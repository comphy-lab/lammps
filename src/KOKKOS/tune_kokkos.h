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

  TuneKokkos(class LAMMPS *, int);
  ~TuneKokkos() override;
  void allocate();
  void tuning_kernel_params();

  int interval;  // # of timesteps for measuring performance in TPS of a given parameter set

  std::vector<int> team_sizes;
  std::vector<int> vector_sizes;
  double* performance;
  int nparams;
  int param_idx;
  int scanning_completed;
  int allocated;

  bigint last_step;
  double last_spcpu;
  int firststep;
  double get_timing_info();
};

}

#endif

