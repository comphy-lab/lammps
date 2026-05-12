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

#include "fix_nh_sphere_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

template <class DeviceType>
FixNHSphereKokkos<DeviceType>::FixNHSphereKokkos(LAMMPS *lmp, int narg, char **arg) :
    FixNHKokkos<DeviceType>(lmp, narg, arg)
{
  if (!atom->omega_flag)
    error->all(FLERR, "Fix {} requires atom attribute omega", style);
  if (!atom->radius_flag)
    error->all(FLERR, "Fix {} requires atom attribute radius", style);

  // inertia = moment of inertia prefactor for sphere or disc

  inertia = 0.4;

  int iarg = 3;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "disc") == 0) {
      inertia = 0.5;
      if (domain->dimension != 2)
        error->all(FLERR, "Fix {} disc option requires 2d simulation", style);
    }
    iarg++;
  }
}

/* ---------------------------------------------------------------------- */

template <class DeviceType>
void FixNHSphereKokkos<DeviceType>::init()
{
  // check that all group particles are finite-size

  double *radius = atom->radius;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit)
      if (radius[i] == 0.0) error->one(FLERR, "Fix {} requires extended particles", style);

  FixNHKokkos<DeviceType>::init();
}

/* ----------------------------------------------------------------------
   perform half-step update of translational and rotational velocities
-----------------------------------------------------------------------*/

template <class DeviceType>
void FixNHSphereKokkos<DeviceType>::nve_v()
{
  // translational velocity update (Kokkos kernel in FixNHKokkos)
  // also sets this->rmass and this->mask views

  FixNHKokkos<DeviceType>::nve_v();

  // rotational velocity (omega) update

  atomKK->sync(execution_space, OMEGA_MASK | TORQUE_MASK | RADIUS_MASK);

  omega_kk  = atomKK->k_omega.view<DeviceType>();
  torque_kk = atomKK->k_torque.view<DeviceType>();
  radius_kk = atomKK->k_radius.view<DeviceType>();

  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNHSphere_nve_v_omega>(0, nlocal), *this);
  copymode = 0;

  atomKK->modified(execution_space, OMEGA_MASK);
}

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixNHSphereKokkos<DeviceType>::operator()(TagFixNHSphere_nve_v_omega, const int &i) const
{
  if (mask(i) & groupbit) {
    const KK_FLOAT dtfrotate  = dtf / inertia;
    const KK_FLOAT dtirotate  = dtfrotate / (radius_kk(i) * radius_kk(i) * rmass(i));
    omega_kk(i, 0) += dtirotate * torque_kk(i, 0);
    omega_kk(i, 1) += dtirotate * torque_kk(i, 1);
    omega_kk(i, 2) += dtirotate * torque_kk(i, 2);
  }
}

/* ----------------------------------------------------------------------
   perform full-step update of positions, and dipole (if requested)
-----------------------------------------------------------------------*/

template <class DeviceType>
void FixNHSphereKokkos<DeviceType>::nve_x()
{
  // position update (Kokkos kernel in FixNHKokkos)
  // also sets this->x, this->v, this->mask views

  FixNHKokkos<DeviceType>::nve_x();

  // dipole orientation update (simple cross-product integrator only)

  if (dipole_flag) {
    if (dlm_flag)
      error->all(FLERR, "Fix {}: DLM dipole integrator not supported in /kk mode", style);

    atomKK->sync(execution_space, MU_MASK | OMEGA_MASK | MASK_MASK);
    mu_kk    = atomKK->k_mu.view<DeviceType>();
    omega_kk = atomKK->k_omega.view<DeviceType>();

    int nlocal = atomKK->nlocal;
    if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

    copymode = 1;
    Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNHSphere_nve_x_dipole>(0, nlocal),
                         *this);
    copymode = 0;

    atomKK->modified(execution_space, MU_MASK);
  }
}

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixNHSphereKokkos<DeviceType>::operator()(TagFixNHSphere_nve_x_dipole, const int &i) const
{
  if (mask(i) & groupbit && mu_kk(i, 3) > 0.0) {
    const KK_FLOAT g0  = mu_kk(i, 0) + dtv * (omega_kk(i, 1) * mu_kk(i, 2) - omega_kk(i, 2) * mu_kk(i, 1));
    const KK_FLOAT g1  = mu_kk(i, 1) + dtv * (omega_kk(i, 2) * mu_kk(i, 0) - omega_kk(i, 0) * mu_kk(i, 2));
    const KK_FLOAT g2  = mu_kk(i, 2) + dtv * (omega_kk(i, 0) * mu_kk(i, 1) - omega_kk(i, 1) * mu_kk(i, 0));
    const KK_FLOAT msq = g0 * g0 + g1 * g1 + g2 * g2;
    const KK_FLOAT scale = mu_kk(i, 3) / Kokkos::sqrt(msq);
    mu_kk(i, 0) = g0 * scale;
    mu_kk(i, 1) = g1 * scale;
    mu_kk(i, 2) = g2 * scale;
  }
}

/* ----------------------------------------------------------------------
   perform half-step thermostat scaling of translational and rotational velocities
-----------------------------------------------------------------------*/

template <class DeviceType>
void FixNHSphereKokkos<DeviceType>::nh_v_temp()
{
  // translational velocity scaling (handles bias removal/restore)

  FixNHKokkos<DeviceType>::nh_v_temp();

  // rotational velocity scaling by same factor_eta

  atomKK->sync(execution_space, OMEGA_MASK | MASK_MASK);
  omega_kk = atomKK->k_omega.view<DeviceType>();

  int nlocal = atomKK->nlocal;
  if (igroup == atomKK->firstgroup) nlocal = atomKK->nfirst;

  copymode = 1;
  Kokkos::parallel_for(Kokkos::RangePolicy<DeviceType, TagFixNHSphere_nh_v_temp_omega>(0, nlocal),
                       *this);
  copymode = 0;

  atomKK->modified(execution_space, OMEGA_MASK);
}

template <class DeviceType>
// NOLINTNEXTLINE
KOKKOS_INLINE_FUNCTION
void FixNHSphereKokkos<DeviceType>::operator()(TagFixNHSphere_nh_v_temp_omega, const int &i) const
{
  if (mask(i) & groupbit) {
    omega_kk(i, 0) *= factor_eta;
    omega_kk(i, 1) *= factor_eta;
    omega_kk(i, 2) *= factor_eta;
  }
}

namespace LAMMPS_NS {
template class FixNHSphereKokkos<LMPDeviceType>;
#ifdef LMP_KOKKOS_GPU
template class FixNHSphereKokkos<LMPHostType>;
#endif
}
