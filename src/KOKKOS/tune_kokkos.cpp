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

TuneKokkos::TuneKokkos(LAMMPS *lmp, int nevery, int _nparams) : Pointers(lmp),
  performance(nullptr), interval(nevery), scanning_completed(0)
{
  ncombinations = 0;
  allocated = 0;
  firststep = 1;
  opt_perf = 0.0;
  relative_tolerance = 0.20; // 20% performance degradation allowed

  num_params = _nparams;
  allocate(num_params);
}

/* ---------------------------------------------------------------------- */

TuneKokkos::~TuneKokkos()
{
  if (performance) delete[] performance;
}

/* ----------------------------------------------------------------------
   allocate arrays to store the parameter combinations and performance data
   num_params = 1: only team size (bond compute and nbin kernels)
                2: team size and threads per atom (typical pair compute kernels)
------------------------------------------------------------------------- */

void TuneKokkos::allocate(int num_params)
{
  if (num_params < 1 || num_params > 2)
    error->all(FLERR,"Illegal number of parameters for Kokkos kernel tuning");

  team_sizes.clear();
  vector_sizes.clear();

  // maximum team size is often limited by the GPU architecture, 512 is a safe choice

  int max_team_size = 512;
  for (int ts = 64; ts <= max_team_size; ts += 32)
    team_sizes.push_back(ts);

  #if defined(KOKKOS_ENABLE_HIP)
  int max_vectorsize = 64;
  #else
  int max_vectorsize = 32;
  #endif

  if (num_params == 2) {
    for (int vs = 1; vs <= max_vectorsize; vs *= 2)
      vector_sizes.push_back(vs);
  } else {
    if (lmp->kokkos->threads_per_atom_set)
      vector_sizes.push_back(lmp->kokkos->threads_per_atom);
    else
      vector_sizes.push_back(1);
  }

  // compute the total number of parameter combinations
  // cols = team sizes (pair/team/size = 64, 96, ..., 512)
  // rows = vector sizes (threads/per/atom = 1, 2, 4, 8, 16, max_vectorsize)

  int num_team_sizes = team_sizes.size();
  int num_vector_sizes = vector_sizes.size();

  ncombinations = num_team_sizes * num_vector_sizes;

  // allocate performance array

  if (performance) delete[] performance;
  performance = new double[ncombinations];
  combination_idx = 0;

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
    double cpu_diff = new_cpu - last_cpu;
    bigint step_diff = new_step - last_step;
    if (step_diff > 0.0) dvalue = cpu_diff / step_diff;
    else dvalue = 0.0;
  }

  last_step = new_step;
  last_cpu = new_cpu;

  return dvalue;
}

/* ----------------------------------------------------------------------
   tuning the pair compute kernel parameters
   this function is called by the /kk pair style at every timestep
   if auto-tuning is enabled
------------------------------------------------------------------------- */

void TuneKokkos::tuning_kernel_params(class Pair *pair)
{
  // ensure that relevant kokkos parameters are allowed to be specified

  if (!lmp->kokkos->neigh_thread_set) {
    lmp->kokkos->neigh_thread_set = 1;
    lmp->kokkos->neigh_thread = 1;
  }
  if (!lmp->kokkos->threads_per_atom_set) lmp->kokkos->threads_per_atom_set = 1;
  if (!lmp->kokkos->pair_team_size_set) lmp->kokkos->pair_team_size_set = 1;

  if (!scanning_completed) {

    // retrieve the current parameter set from combination_idx

    int num_team_sizes = team_sizes.size();
    int pair_team_size = team_sizes[combination_idx % num_team_sizes];
    int threads_per_atom = vector_sizes[combination_idx / num_team_sizes];

    // set the KOKKOS kernel parameters

    lmp->kokkos->pair_team_size = pair_team_size;
    lmp->kokkos->threads_per_atom = threads_per_atom;

    // wait for interval timesteps to collect timing info

    if (update->ntimestep % interval == 0) {

      double tps = get_timing_info();

      // skip for the first call when no timing info is available

      if (tps == 0.0) return;

      // store the performance of the current parameter set

      performance[combination_idx] = 1.0 / tps;

      #ifdef TUNE_DEBUG
      if (comm->me == 0) {
        std::string mesg = fmt::format("combination_idx {}: pair team size = {} threads per atom = {} ",
                          combination_idx, lmp->kokkos->pair_team_size, lmp->kokkos->threads_per_atom);
        mesg += fmt::format("perf = {:.1f} TPS\n", performance[combination_idx]);
        utils::logmesg(lmp,mesg);
      }
      #endif

      // move to the next parameter set

      combination_idx++;

      // suppose that interval is sufficiently long to get a stable TPS
      // so that we only need a single pass over all the parameter combinations

      if (combination_idx >= ncombinations) scanning_completed = 1;
    }
  }

  // if scanning just completed, find the parameter set with the optimal performance

  if (scanning_completed && combination_idx == ncombinations) {

    int opt_idx = get_optimal_combination_idx();
    int num_team_sizes = team_sizes.size();
    int pair_team_size_opt = team_sizes[opt_idx % num_team_sizes];
    int threads_per_atom_opt = vector_sizes[opt_idx / num_team_sizes];

    // set the optimal kernel parameters

    lmp->kokkos->pair_team_size = pair_team_size_opt;
    lmp->kokkos->threads_per_atom = threads_per_atom_opt;

    // reset combination_idx to zero to be ready for another scan if needed
    // and to avoid repeating this block

    combination_idx = 0;
  }

  // check if the performance is within acceptable range of the optimal performance
  // if not, re-trigger the scanning process

  regular_performance_check();
}

/* ----------------------------------------------------------------------
   tuning the bond compute kernel parameters
   this function is called by the /kk bond style at every timestep
   if auto-tuning is enabled
------------------------------------------------------------------------- */

void TuneKokkos::tuning_kernel_params(class Bond *bond)
{
  // ensure that relevant kokkos parameters are allowed to be specified

  if (!lmp->kokkos->bond_block_size_set) lmp->kokkos->bond_block_size_set = 1;

  if (!scanning_completed) {

    // retrieve the current parameter set from combination_idx

    int bond_block_size = team_sizes[combination_idx];

    // set the KOKKOS kernel parameters

    lmp->kokkos->bond_block_size = bond_block_size;

    // wait for interval timesteps to collect timing info

    if (update->ntimestep % interval == 0) {

      double tps = get_timing_info();

      // skip for the first call when no timing info is available

      if (tps == 0.0) return;

      // store the performance of the current parameter set

      performance[combination_idx] = 1.0 / tps;

      #ifdef TUNE_DEBUG
      if (comm->me == 0) {
        std::string mesg = fmt::format("combination_idx {}: bond block size = {} ",
                          combination_idx, lmp->kokkos->bond_block_size);
        mesg += fmt::format("perf = {:.1f} TPS\n", performance[combination_idx]);
        utils::logmesg(lmp,mesg);
      }
      #endif

      // move to the next parameter set

      combination_idx++;

      // suppose that interval is sufficiently long to get a stable TPS
      // so that we only need a single pass over all the parameter combinations

      if (combination_idx >= ncombinations) scanning_completed = 1;
    }
  }

  // if scanning just completed, find the parameter set with the optimal performance

  if (scanning_completed && combination_idx == ncombinations) {

    int opt_idx = get_optimal_combination_idx();
    int bond_block_size = team_sizes[opt_idx];

    // set the optimal kernel parameters

    lmp->kokkos->bond_block_size = bond_block_size;

    // reset combination_idx to zero to avoid repeating this block

    combination_idx = 0;
  }

  // check if the performance is within acceptable range of the optimal performance
  // if not, re-trigger the scanning process

  regular_performance_check();
}

/* ----------------------------------------------------------------------
   regularly monitor the current performance after scanning is completed
   if performance degraded beyond acceptable threshold,
     re-trigger the scanning process
   otherwise, continue using the optimal parameter set,
     including with a new run
------------------------------------------------------------------------- */

void TuneKokkos::regular_performance_check()
{
  if (!scanning_completed || update->ntimestep % interval != 0) return;
  if (update->ntimestep == update->beginstep) return;

  double tps = get_timing_info();
  if (tps <= 0.0) return;
  double perf = 1.0 / tps;

  #ifdef TUNE_DEBUG
  if (comm->me == 0) {
    std::string mesg = fmt::format("Using the optimal params: ");
    mesg += fmt::format("pair team size = {} threads per atom = {} ",
                      lmp->kokkos->pair_team_size,
                      lmp->kokkos->threads_per_atom);
    mesg += fmt::format(" current perf = {:.1f} TPS\n", perf);
    utils::logmesg(lmp,mesg);
  }
  #endif

  // compute the relative performance difference wrt the optimal performance

  double diff;
  if (opt_perf > 0.0)
    diff = (opt_perf - perf) / opt_perf;
  else
    diff = 0.0;

  if (diff > relative_tolerance) {
    scanning_completed = 0;
    combination_idx = 0;
    firststep = 1;

    #ifdef TUNE_DEBUG
    if (comm->me == 0) {
      std::string mesg = fmt::format("Performance degraded by {:.2f}%, re-triggering scan\n",
                        diff * 100.0);
      utils::logmesg(lmp,mesg);
    }
    #endif
  }
}

/* ----------------------------------------------------------------------
   find the optimal performance
   return the index of the optimal parameter set
------------------------------------------------------------------------- */

int TuneKokkos::get_optimal_combination_idx()
{
  if (ncombinations == 0 || performance == nullptr)
    error->all(FLERR,"No performance data available for Kokkos kernel tuning");

  opt_perf = performance[0];
  int opt_idx = 0;
  for (int i = 1; i < ncombinations; i++) {
    if (performance[i] > opt_perf) {
      opt_perf = performance[i];
      opt_idx = i;
    }
  }
  return opt_idx;
}





