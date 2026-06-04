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

/* ----------------------------------------------------------------------
   Contributing authors: Axel Kohlmeyer (Temple)
------------------------------------------------------------------------- */

#include "ewald_gpu.h"

#include "atom.h"
#include "error.h"
#include "gpu_extra.h"
#include "math_const.h"
#include "neighbor.h"

using namespace LAMMPS_NS;
using namespace MathConst;

// external functions from gpu library

#define EWALD_GPU_API(api)  ewald_gpu_ ## api ## _d

void EWALD_GPU_API(init)(const int nlocal, const int nall, FILE *screen,
                         int &success);
void EWALD_GPU_API(setup)(const int kmax, const int kcount, int *kxvecs,
                          int *kyvecs, int *kzvecs, double *ug, double **eg,
                          double **vg, double *unitk, int &success);
int EWALD_GPU_API(structure)(const int ago, const int nlocal, const int nall,
                             double **host_x, int *host_type, double *host_q,
                             double *host_sfacrl, double *host_sfacim,
                             bool &success);
void EWALD_GPU_API(compute)(double *host_sfacrl_all, double *host_sfacim_all,
                            const double qscale, const int slabflag,
                            const int eflag_atom, const int vflag_atom,
                            double *host_eatom, double **host_vatom,
                            bool &success);
void EWALD_GPU_API(clear)(const double cpu_time);
double EWALD_GPU_API(bytes)();

/* ---------------------------------------------------------------------- */

EwaldGPU::EwaldGPU(LAMMPS *lmp) : Ewald(lmp)
{
  GPU_EXTRA::gpu_ready(lmp->modify, lmp->error);
}

/* ----------------------------------------------------------------------
   free all memory
------------------------------------------------------------------------- */

EwaldGPU::~EwaldGPU()
{
  EWALD_GPU_API(clear)(0.0);
}

/* ----------------------------------------------------------------------
   called once before run
------------------------------------------------------------------------- */

void EwaldGPU::init()
{
  // initialize the GPU device and atom storage first, so that the device is
  // ready when Ewald::init() -> setup() uploads the k-space coefficients

  int success = 0;
  EWALD_GPU_API(init)(atom->nlocal, atom->nlocal+atom->nghost, screen, success);
  GPU_EXTRA::check_flag(success, error, world);

  Ewald::init();
}

/* ----------------------------------------------------------------------
   called whenever the box changes: refresh the k-space coefficients on
   the device after the base class recomputes them
------------------------------------------------------------------------- */

void EwaldGPU::setup()
{
  Ewald::setup();

  int success = 0;
  EWALD_GPU_API(setup)(kmax, kcount, kxvecs, kyvecs, kzvecs, ug, eg, vg, unitk,
                       success);
  GPU_EXTRA::check_flag(success, error, world);
}

/* ----------------------------------------------------------------------
   compute the EwaldGPU long-range force, energy, virial
------------------------------------------------------------------------- */

void EwaldGPU::compute(int eflag, int vflag)
{
  ev_init(eflag, vflag);

  if (atom->natoms != natoms_original) {
    qsum_qsq();
    natoms_original = atom->natoms;
  }
  if (qsqsum == 0.0) return;

  // triclinic boxes use the host path; Ewald::compute() still gets its
  // structure factors from the device via the overridden eik_dot_r()

  if (triclinic) {
    Ewald::compute(eflag, vflag);
    return;
  }

  const int nlocal = atom->nlocal;
  const int nall = atom->nlocal + atom->nghost;
  double *q = atom->q;

  // structure factors on the device

  bool success = true;
  EWALD_GPU_API(structure)(neighbor->ago, nlocal, nall, atom->x, atom->type,
                           atom->q, sfacrl, sfacim, success);
  if (!success)
    error->one(FLERR, "Insufficient memory on accelerator for ewald/gpu");

  MPI_Allreduce(sfacrl, sfacrl_all, kcount, MPI_DOUBLE, MPI_SUM, world);
  MPI_Allreduce(sfacim, sfacim_all, kcount, MPI_DOUBLE, MPI_SUM, world);

  const double qscale = qqrd2e * scale;

  // per-atom field/force on the device (force queued for fix gpu to merge
  // with the pair force; raw per-atom energy/virial copied into eatom/vatom)

  EWALD_GPU_API(compute)(sfacrl_all, sfacim_all, qscale, slabflag, eflag_atom,
                         vflag_atom, eatom, vatom, success);
  if (!success)
    error->one(FLERR, "Insufficient memory on accelerator for ewald/gpu");

  // global energy across k-vectors plus self-energy and volume corrections

  if (eflag_global) {
    for (int k = 0; k < kcount; k++)
      energy += ug[k] * (sfacrl_all[k]*sfacrl_all[k] +
                         sfacim_all[k]*sfacim_all[k]);
    energy -= g_ewald*qsqsum/MY_PIS +
      MY_PI2*qsum*qsum / (g_ewald*g_ewald*volume);
    energy *= qscale;
  }

  // global virial

  if (vflag_global) {
    double uk;
    for (int k = 0; k < kcount; k++) {
      uk = ug[k] * (sfacrl_all[k]*sfacrl_all[k] + sfacim_all[k]*sfacim_all[k]);
      for (int j = 0; j < 6; j++) virial[j] += uk*vg[k][j];
    }
    for (int j = 0; j < 6; j++) virial[j] *= qscale;
  }

  // per-atom energy/virial: the device returned the raw k-vector sums; apply
  // the q_i factor, self-energy correction, and qscale here to match the CPU

  if (evflag_atom) {
    if (eflag_atom) {
      for (int i = 0; i < nlocal; i++) {
        eatom[i] = q[i]*eatom[i] - (g_ewald*q[i]*q[i]/MY_PIS +
          MY_PI2*q[i]*qsum / (g_ewald*g_ewald*volume));
        eatom[i] *= qscale;
      }
    }
    if (vflag_atom)
      for (int i = 0; i < nlocal; i++)
        for (int j = 0; j < 6; j++) vatom[i][j] *= q[i]*qscale;
  }

  // 2d slab correction (host; adds to atom->f and the energy)

  if (slabflag == 1) slabcorr();
}

/* ----------------------------------------------------------------------
   structure factors on the device, used by the host field/force loop in
   Ewald::compute() for the per-atom / triclinic fallback path
------------------------------------------------------------------------- */

void EwaldGPU::eik_dot_r()
{
  // host cs/sn are still needed by the host field/force loop in the fallback
  Ewald::eik_dot_r();

  // recompute the (local) structure factors on the device and overwrite the
  // host values, so the GPU result drives the energy and forces

  bool success = true;
  EWALD_GPU_API(structure)(neighbor->ago, atom->nlocal,
                           atom->nlocal+atom->nghost, atom->x, atom->type,
                           atom->q, sfacrl, sfacim, success);
  if (!success)
    error->one(FLERR, "Insufficient memory on accelerator for ewald/gpu");
}

/* ---------------------------------------------------------------------- */

double EwaldGPU::memory_usage()
{
  double bytes = Ewald::memory_usage();
  bytes += EWALD_GPU_API(bytes)();
  return bytes;
}
