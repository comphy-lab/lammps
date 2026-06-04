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
#include "neighbor.h"

using namespace LAMMPS_NS;

// external functions from gpu library

#define EWALD_GPU_API(api)  ewald_gpu_ ## api ## _d

void EWALD_GPU_API(init)(const int nlocal, const int nall, FILE *screen,
                         int &success);
void EWALD_GPU_API(setup)(const int kmax, const int kcount, int *kxvecs,
                          int *kyvecs, int *kzvecs, double *unitk,
                          int &success);
int EWALD_GPU_API(structure)(const int ago, const int nlocal, const int nall,
                             double **host_x, int *host_type, double *host_q,
                             double *host_sfacrl, double *host_sfacim,
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
  EWALD_GPU_API(setup)(kmax, kcount, kxvecs, kyvecs, kzvecs, unitk, success);
  GPU_EXTRA::check_flag(success, error, world);
}

/* ----------------------------------------------------------------------
   compute the EwaldGPU long-range force, energy, virial
------------------------------------------------------------------------- */

void EwaldGPU::compute(int eflag, int vflag)
{
  // Ewald::compute() dispatches to the virtual eik_dot_r() below for the
  // structure factors (orthogonal boxes) or eik_dot_r_triclinic() on the
  // host (triclinic fallback); the remaining field/force/energy work still
  // runs on the host in this stage.

  Ewald::compute(eflag, vflag);
}

/* ----------------------------------------------------------------------
   structure factors: computed on the device, then used by the host
   field/force loop in Ewald::compute()
------------------------------------------------------------------------- */

void EwaldGPU::eik_dot_r()
{
  // host cs/sn are still needed by the host field/force loop in this stage
  Ewald::eik_dot_r();

  // recompute the (local) structure factors on the device and overwrite
  // the host values, so the GPU path drives the energy and forces

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
