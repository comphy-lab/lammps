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

using namespace LAMMPS_NS;

// external functions from gpu library

#define EWALD_GPU_API(api)  ewald_gpu_ ## api ## _d

void EWALD_GPU_API(init)(const int nlocal, const int nall, FILE *screen,
                         int &success);
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
  Ewald::init();

  // GPU device init (plumbing only in this stage; the reciprocal-space
  // computation still runs on the host through Ewald::compute())

  int success = 0;
  EWALD_GPU_API(init)(atom->nlocal, atom->nlocal+atom->nghost, screen, success);
  GPU_EXTRA::check_flag(success, error, world);
}

/* ----------------------------------------------------------------------
   compute the EwaldGPU long-range force, energy, virial
------------------------------------------------------------------------- */

void EwaldGPU::compute(int eflag, int vflag)
{
  Ewald::compute(eflag, vflag);
}

/* ---------------------------------------------------------------------- */

double EwaldGPU::memory_usage()
{
  double bytes = Ewald::memory_usage();
  bytes += EWALD_GPU_API(bytes)();
  return bytes;
}
