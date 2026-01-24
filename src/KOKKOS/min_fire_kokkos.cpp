// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "min_fire_kokkos.h"
#include "atom_kokkos.h"
#include "atom_masks.h"
#include "force.h"
#include "update.h"
#include "timer.h"
#include "error.h"
#include "universe.h"
#include <cmath>

using namespace LAMMPS_NS;

static constexpr double EPS_ENERGY = 1.0e-8;

MinFireKokkos::MinFireKokkos(LAMMPS *lmp) : MinKokkos(lmp) {
  atomKK = (AtomKokkos *) atom;
  kokkosable = 1;
}

void MinFireKokkos::init() {
  MinKokkos::init();

  if (tmax < tmin) error->all(FLERR, "tmax has to be larger than tmin");
  if (dtgrow < 1.0) error->all(FLERR, "dtgrow has to be larger than 1.0");
  if (dtshrink > 1.0) error->all(FLERR, "dtshrink has to be smaller than 1.0");

  alpha = alpha0;
  last_negative = ntimestep_start = update->ntimestep;
  vdotf_negatif = 0;
}

void MinFireKokkos::setup_style() {
  atomKK->sync(Device, V_MASK);
  auto l_v = atomKK->k_v.view_device();
  int nlocal = atom->nlocal;

  Kokkos::parallel_for("min_fire/zero_v", nlocal, LAMMPS_LAMBDA(const int i) {
    l_v(i,0) = l_v(i,1) = l_v(i,2) = 0.0;
  });
  atomKK->modified(Device, V_MASK);
  flagv0 = 1;
}

void MinFireKokkos::reset_vectors() {
  nvec = 3 * atom->nlocal;
  if (nvec) {
    auto d_x = atomKK->k_x.view_device();
    auto d_f = atomKK->k_f.view_device();
    xvec = DAT::t_kkfloat_1d(d_x.data(), nvec);
    fvec = DAT::t_kkacc_1d(d_f.data(), nvec);
  }
}

int MinFireKokkos::iterate(int maxiter) {
  if (integrator == EULERIMPLICIT) {
    return abcflag ? run_iterate<EULERIMPLICIT, true>(maxiter) : run_iterate<EULERIMPLICIT, false>(maxiter);
  } else if (integrator == VERLET) {
    return abcflag ? run_iterate<VERLET, true>(maxiter) : run_iterate<VERLET, false>(maxiter);
  } else if (integrator == LEAPFROG) {
    return abcflag ? run_iterate<LEAPFROG, true>(maxiter) : run_iterate<LEAPFROG, false>(maxiter);
  } else if (integrator == EULEREXPLICIT) {
    return abcflag ? run_iterate<EULEREXPLICIT, true>(maxiter) : run_iterate<EULEREXPLICIT, false>(maxiter);
  }
  return MAXITER;
}

template <int INTEGRATOR, bool ABCFLAG>
int MinFireKokkos::run_iterate(int maxiter) {
  double vdotf_local, vdotfall, vdotv_local, vdotvall, fdotf_local, fdotfall;
  double scale1, scale2;
  double dtv;

  alpha_final = 0.0;
  auto l_dt = update->dt;

  if (INTEGRATOR == LEAPFROG) {
    energy_force(0);
    neval++;
    atomKK->sync(Device, V_MASK | F_MASK);
    auto l_v = atomKK->k_v.view_device();
    auto l_f = atomKK->k_f.view_device();
    auto l_rmass = atomKK->k_rmass.view_device();
    auto l_mass = atomKK->k_mass.view_device();
    auto l_type = atomKK->k_type.view_device();
    auto dtf = -0.5 * l_dt * force->ftm2v;

    Kokkos::parallel_for("min_fire/leapfrog_init", atom->nlocal, LAMMPS_LAMBDA(const int i) {
      double dtfm = dtf / (l_rmass.data() ? l_rmass(i) : l_mass(l_type(i)));
      l_v(i,0) = dtfm * l_f(i,0);
      l_v(i,1) = dtfm * l_f(i,1);
      l_v(i,2) = dtfm * l_f(i,2);
    });
    atomKK->modified(Device, V_MASK);
  }

  for (int iter = 0; iter < maxiter; iter++) {
    if (timer->check_timeout(niter)) return TIMEOUT;

    bigint ntimestep = ++update->ntimestep;
    niter++;

    atomKK->sync(Device, X_MASK | V_MASK | F_MASK);
    auto l_v = atomKK->k_v.view_device();
    auto l_f = atomKK->k_f.view_device();
    auto l_x = atomKK->k_x.view_device();
    auto l_rmass = atomKK->k_rmass.view_device();
    auto l_mass = atomKK->k_mass.view_device();
    auto l_type = atomKK->k_type.view_device();
    auto dtmin = tmin * l_dt;
    auto dtmax = tmax * l_dt;
  
    int nlocal = atom->nlocal;

    vdotf_local = 0.0;
    Kokkos::parallel_reduce("min_fire/vdotf", nlocal, LAMMPS_LAMBDA(const int i, double &vdf) {
      vdf += l_v(i,0)*l_f(i,0) + l_v(i,1)*l_f(i,1) + l_v(i,2)*l_f(i,2);
    }, vdotf_local);
    MPI_Allreduce(&vdotf_local, &vdotfall, 1, MPI_DOUBLE, MPI_SUM, world);

    if (update->multireplica == 1) {
      vdotf_local = vdotfall;
      MPI_Allreduce(&vdotf_local, &vdotfall, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
    }

    if (vdotfall > 0.0) {
      vdotv_local = 0.0;
      fdotf_local = 0.0;
      Kokkos::parallel_reduce("min_fire/norms", nlocal, LAMMPS_LAMBDA(const int i, double &vv) {
        vv += l_v(i,0)*l_v(i,0) + l_v(i,1)*l_v(i,1) + l_v(i,2)*l_v(i,2);
      }, vdotv_local);
      Kokkos::parallel_reduce("min_fire/fnorms", nlocal, LAMMPS_LAMBDA(const int i, double &ff) {
        ff += l_f(i,0)*l_f(i,0) + l_f(i,1)*l_f(i,1) + l_f(i,2)*l_f(i,2);
      }, fdotf_local);

      MPI_Allreduce(&vdotv_local, &vdotvall, 1, MPI_DOUBLE, MPI_SUM, world);
      MPI_Allreduce(&fdotf_local, &fdotfall, 1, MPI_DOUBLE, MPI_SUM, world);

      if (update->multireplica == 1) {
        double dot_tmp = vdotvall;
        MPI_Allreduce(&dot_tmp, &vdotvall, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
        dot_tmp = fdotfall;
        MPI_Allreduce(&dot_tmp, &fdotfall, 1, MPI_DOUBLE, MPI_SUM, universe->uworld);
      }

      if (ABCFLAG) {
        if (alpha < 1e-10) alpha = 1e-10;
        double abc = (1.0 - pow(1.0 - alpha, (double)(ntimestep - last_negative)));
        scale1 = (1.0 - alpha) / abc;
        scale2 = (fdotfall <= 1e-20) ? 0.0 : (alpha * sqrt(vdotvall / fdotfall)) / abc;
      } else {
        scale1 = 1.0 - alpha;
        scale2 = (fdotfall <= 1e-20) ? 0.0 : alpha * sqrt(vdotvall / fdotfall);
      }

      if (ntimestep - last_negative > delaystep) {
        l_dt = fmin(l_dt * dtgrow, dtmax);
        update->dt = l_dt;
        alpha *= alphashrink;
      }
      vdotf_negatif = 0;
    } else {
      last_negative = ntimestep;
      if (ntimestep - ntimestep_start >= delaystep || !delaystep_start_flag) {
        alpha = alpha0;
        if (l_dt * dtshrink >= dtmin) {
          l_dt *= dtshrink;
          update->dt = l_dt;
        }
      }

      vdotf_negatif++;
      if (max_vdotf_negatif > 0 && vdotf_negatif > max_vdotf_negatif) return MAXVDOTF;

      Kokkos::parallel_for("min_fire/inertia_reset", nlocal, LAMMPS_LAMBDA(const int i) {
        if (halfstepback_flag) {
          l_x(i,0) -= 0.5 * l_dt * l_v(i,0);
          l_x(i,1) -= 0.5 * l_dt * l_v(i,1);
          l_x(i,2) -= 0.5 * l_dt * l_v(i,2);
        }
        l_v(i,0) = l_v(i,1) = l_v(i,2) = 0.0;
      });
      atomKK->modified(Device, X_MASK | V_MASK);
      flagv0 = 1;
    }

    if (!ABCFLAG && flagv0) {
      energy_force(0);
      neval++;
      double dtf_init = l_dt * force->ftm2v;
      Kokkos::parallel_for("min_fire/v_init", nlocal, LAMMPS_LAMBDA(const int i) {
        double dtfm = dtf_init / (l_rmass.data() ? l_rmass(i) : l_mass(l_type(i)));
        l_v(i,0) = dtfm * l_f(i,0);
        l_v(i,1) = dtfm * l_f(i,1);
        l_v(i,2) = dtfm * l_f(i,2);
      });
      atomKK->modified(Device, V_MASK);
    }

    double dtvone = l_dt;
    if (!ABCFLAG) {
      Kokkos::parallel_reduce("min_fire/dtv_limit", nlocal, LAMMPS_LAMBDA(const int i, double &dtmin_local) {
        double vmax = fmax(fabs(l_v(i,0)), fmax(fabs(l_v(i,1)), fabs(l_v(i,2))));
        if (dtmin_local * vmax > dmax) dtmin_local = dmax / vmax;
      }, Kokkos::Min<double>(dtvone));
    }
    MPI_Allreduce(&dtvone, &dtv, 1, MPI_DOUBLE, MPI_MIN, world);
    if (update->multireplica == 1) {
      double dtv_tmp = dtv;
      MPI_Allreduce(&dtv_tmp, &dtv, 1, MPI_DOUBLE, MPI_MIN, universe->uworld);
    }

    if (flagv0) {
      Kokkos::parallel_for("min_fire/final_v_zero", nlocal, LAMMPS_LAMBDA(const int i) {
        l_v(i,0) = l_v(i,1) = l_v(i,2) = 0.0;
      });
    }

    double dtf_final = dtv * force->ftm2v;
    double dtf_half = 0.5 * dtf_final;

    Kokkos::parallel_for("min_fire/integrate", nlocal, LAMMPS_LAMBDA(const int i) {
      double mass_val = (l_rmass.data() ? l_rmass(i) : l_mass(l_type(i)));
      double dtfm = dtf_final / mass_val;
      double dtfm_half = dtf_half / mass_val;

      if (INTEGRATOR == EULERIMPLICIT || INTEGRATOR == LEAPFROG) {
        l_v(i,0) += dtfm * l_f(i,0);
        l_v(i,1) += dtfm * l_f(i,1);
        l_v(i,2) += dtfm * l_f(i,2);
        if (vdotfall > 0.0) {
          l_v(i,0) = scale1 * l_v(i,0) + scale2 * l_f(i,0);
          l_v(i,1) = scale1 * l_v(i,1) + scale2 * l_f(i,1);
          l_v(i,2) = scale1 * l_v(i,2) + scale2 * l_f(i,2);
        }
        l_x(i,0) += dtv * l_v(i,0);
        l_x(i,1) += dtv * l_v(i,1);
        l_x(i,2) += dtv * l_v(i,2);
      } else if (INTEGRATOR == VERLET) {
        l_v(i,0) += dtfm_half * l_f(i,0);
        l_v(i,1) += dtfm_half * l_f(i,1);
        l_v(i,2) += dtfm_half * l_f(i,2);
        if (vdotfall > 0.0) {
          l_v(i,0) = scale1 * l_v(i,0) + scale2 * l_f(i,0);
          l_v(i,1) = scale1 * l_v(i,1) + scale2 * l_f(i,1);
          l_v(i,2) = scale1 * l_v(i,2) + scale2 * l_f(i,2);
        }
        l_x(i,0) += dtv * l_v(i,0);
        l_x(i,1) += dtv * l_v(i,1);
        l_x(i,2) += dtv * l_v(i,2);
      } else if (INTEGRATOR == EULEREXPLICIT) {
        if (vdotfall > 0.0) {
          l_v(i,0) = scale1 * l_v(i,0) + scale2 * l_f(i,0);
          l_v(i,1) = scale1 * l_v(i,1) + scale2 * l_f(i,1);
          l_v(i,2) = scale1 * l_v(i,2) + scale2 * l_f(i,2);
        }
        l_x(i,0) += dtv * l_v(i,0);
        l_x(i,1) += dtv * l_v(i,1);
        l_x(i,2) += dtv * l_v(i,2);
        l_v(i,0) += dtfm * l_f(i,0);
        l_v(i,1) += dtfm * l_f(i,1);
        l_v(i,2) += dtfm * l_f(i,2);
      }
    });

    atomKK->modified(Device, X_MASK | V_MASK);
    eprevious = ecurrent;
    ecurrent = energy_force(0);
    neval++;

    if (INTEGRATOR == VERLET) {
      atomKK->sync(Device, V_MASK | F_MASK);
      Kokkos::parallel_for("min_fire/verlet_v_final", nlocal, LAMMPS_LAMBDA(const int i) {
        double dtfm_half = dtf_half / (l_rmass.data() ? l_rmass(i) : l_mass(l_type(i)));
        l_v(i,0) += dtfm_half * l_f(i,0);
        l_v(i,1) += dtfm_half * l_f(i,1);
        l_v(i,2) += dtfm_half * l_f(i,2);
      });
      atomKK->modified(Device, V_MASK);
    }

    flagv0 = 0;

    // Criteria checks
    if (update->etol > 0.0 && ntimestep - last_negative > delaystep) {
      if (fabs(ecurrent - eprevious) < update->etol * 0.5 * (fabs(ecurrent) + fabs(eprevious) + EPS_ENERGY)) {
        if (update->multireplica == 0) return ETOL;
        int flag = 0, flagall;
        MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, universe->uworld);
        if (flagall == 0) return ETOL;
      }
    }

    if (update->ftol > 0.0) {
      double fdotf = (normstyle == MAX) ? fnorm_max() : (normstyle == INF ? fnorm_inf() : fnorm_sqr());
      if (fdotf < update->ftol * update->ftol) {
        if (update->multireplica == 0) return FTOL;
        int flag = 0, flagall;
        MPI_Allreduce(&flag, &flagall, 1, MPI_INT, MPI_SUM, universe->uworld);
        if (flagall == 0) return FTOL;
      }
    }

    if (output->next == ntimestep) {
      timer->stamp();
      output->write(ntimestep);
      timer->stamp(Timer::OUTPUT);
    }
  }
  return MAXITER;
}
