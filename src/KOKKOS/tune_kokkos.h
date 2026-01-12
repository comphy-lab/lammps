// clang-format off
/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_TUNE_KOKKOS_H
#define LMP_TUNE_KOKKOS_H

#include "pointers.h"
#include <mpi.h>
#include <vector>

namespace LAMMPS_NS {

class TuneKokkos : protected Pointers {
 public:
  TuneKokkos(class LAMMPS *, int kernel_type, int nevery, int nparams=2);
  ~TuneKokkos() override;
  void allocate(int);
  void tuning_kernel_params();

  enum { PAIR, BOND, NBIN, FIX, COMPUTE };

  int kernel_type;  // type of kernel being tuned: PAIR, BOND, NBIN, FIX, COMPUTE
  int interval;     // # of timesteps to run the simulation with a given parameter set
                    //   should be sufficiently large (~100) for estimating the overall performance in TPS
                    //   should be small enough to avoid significant overhead during auto-tuning

  std::vector<int> team_sizes;   // parameter values for the team size (typically, thread block size)
  std::vector<int> vector_sizes; // parameter values for the vector size (the 2nd dimension of thread block)

  int num_params;            // number of parameters to tune: 1 (team size only) or 2 (team size and threads per atom)
  double* performance;       // array to store the performance data for each parameter set
  int ncombinations;         // total number of parameter combinations
  int combination_idx;       // current combination index during scanning
  int scanning_completed;    // 0 if still scanning, 1 if scanning completed
  int allocated;             // 1 if the performance array is allocated and param values set up

  double opt_perf;           // stored the optimal performance
  double relative_tolerance; // acceptable threshold for performance degradation wrt opt_perf

  bigint last_step;          // last timestep when timing info was collected
  double last_cpu;           // last CPU time when timing info was collected
  int firststep;             // 1 if first timestep for timing info collection

  double get_timing_info();                  // get the elapsed time from the last call
  void get_current_params(int, auto&, auto&);  // get the team size and vector size for a given combination index
  void set_param_values(int);                // set the KOKKOS kernel parameters based on the combination index
  int get_optimal_combination_idx();         // find the combination index with the best performance
  void regular_performance_check();          // monitor the performance during normal simulation and re-trigger scanning if needed
};

}

#endif

