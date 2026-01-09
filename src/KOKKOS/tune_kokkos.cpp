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
   Contributing author: Trung Nguyen (U Chicago)
------------------------------------------------------------------------- */

#include "tune_kokkos.h"

#include "comm.h"
#include "lammps.h"
#include "kokkos.h"
#include "error.h"
#include "timer.h"
#include "update.h"

using namespace LAMMPS_NS;

//#define TUNE_DEBUG

/* ---------------------------------------------------------------------- */

TuneKokkos::TuneKokkos(LAMMPS *lmp, int nevery) : Pointers(lmp),
  performance(nullptr), interval(nevery), scanning_completed(0)
{
  nparams = 0;
  allocated = 0;
  firststep = 1;
}

/* ---------------------------------------------------------------------- */

TuneKokkos::~TuneKokkos()
{
  if (performance) delete[] performance;
}

/* ---------------------------------------------------------------------- */

void TuneKokkos::allocate(int num_params)
{
  team_sizes.clear();
  vector_sizes.clear();

  #if defined(KOKKOS_ENABLE_HIP)
  int max_vectorsize = 64;
  #else
  int max_vectorsize = 32;
  #endif

  for (int ts = 64; ts <= 512; ts += 32)
    team_sizes.push_back(ts);

  if (num_params == 2) {
    for (int vs = 1; vs <= max_vectorsize; vs *= 2)
      vector_sizes.push_back(vs);
  }

  // compute the total number of parameter combinations
  // cols = team sizes (pair/team/size = 64, 96, ..., 512)
  // rows = vector sizes (threads/per/atom = 1, 2, 4, 8, 16, max_vectorsize)

  int num_team_sizes = team_sizes.size();
  int num_vector_sizes;
  if (num_params == 2) num_vector_sizes = vector_sizes.size();
  else num_vector_sizes = 1;

  nparams = num_team_sizes * num_vector_sizes;

  // allocate performance array

  performance = new double[nparams];
  param_idx = 0;

  scanning_completed = 0;
  allocated = 1;
}

/* ----------------------------------------------------------------------
   figure out CPU time per timestep since last time checked
------------------------------------------------------------------------- */

double TuneKokkos::get_timing_info()
{
  double dvalue;
  double new_cpu;
  bigint new_step = update->ntimestep;

  if (firststep == 1) {
    new_cpu = 0.0;
    dvalue = 0.0;
    firststep = 0;
  } else {
    new_cpu = timer->elapsed(Timer::TOTAL);
    double cpu_diff = new_cpu - last_spcpu;
    bigint step_diff = new_step - last_step;
    if (step_diff > 0.0) dvalue = cpu_diff/step_diff;
    else dvalue = 0.0;
  }

  last_step = new_step;
  last_spcpu = new_cpu;

  return dvalue;
}

/* ---------------------------------------------------------------------- */

void TuneKokkos::tuning_kernel_params(class Pair *pair)
{
  if (!allocated) allocate(2);

  // ensure that relevant kokkos parameters are allowed to be specified

  if (!lmp->kokkos->neigh_thread_set) lmp->kokkos->neigh_thread_set = 1;
  if (!lmp->kokkos->threads_per_atom_set) lmp->kokkos->threads_per_atom_set = 1;
  if (!lmp->kokkos->pair_team_size_set) lmp->kokkos->pair_team_size_set = 1;

  if (!scanning_completed) {

    if (update->ntimestep % interval == 0) {

      double tps = get_timing_info();
      if (tps == 0.0) return;

      performance[param_idx] = 1.0 / tps;
      int num_team_sizes = team_sizes.size();
      int pair_team_size = team_sizes[param_idx % num_team_sizes];
      int threads_per_atom = vector_sizes[param_idx / num_team_sizes];

      #ifdef TUNE_DEBUG
      if (comm->me == 0)
        printf("param_idx = %d: pair team size = %d: tpa = %d: perf = %f TPS\n",
          param_idx, pair_team_size, threads_per_atom, performance[param_idx]);
      #endif

      lmp->kokkos->pair_team_size = pair_team_size;
      lmp->kokkos->threads_per_atom = threads_per_atom;

      param_idx++;

      // suppose that interval is sufficiently long to get a stable TPS
      // so that we only need a single pass over all the parameter combinations

      if (param_idx >= nparams) scanning_completed = 1;
    }

  } else {

    // find the optimal parameter set

    double max_perf = performance[0];
    int best_idx = 0;
    for (int i = 1; i < nparams; i++) {
      if (performance[i] > max_perf) {
        max_perf = performance[i];
        best_idx = i;
      }
    }

    int num_team_sizes = team_sizes.size();
    int pair_team_size_opt = team_sizes[best_idx % num_team_sizes];
    int threads_per_atom_opt = vector_sizes[best_idx / num_team_sizes];

    // set the optimal kernel parameters

    lmp->kokkos->pair_team_size = pair_team_size_opt;
    lmp->kokkos->threads_per_atom = threads_per_atom_opt;

    #ifdef TUNE_DEBUG
    if (update->ntimestep % interval == 0) {
      double tps = get_timing_info();
      double perf = 1.0 / tps;
      if (comm->me == 0)
        printf("optimal pair team size = %d: tpa = %d: tps = %f\n",
          lmp->kokkos->pair_team_size, lmp->kokkos->threads_per_atom, perf);

    }
    #endif
  }
}




