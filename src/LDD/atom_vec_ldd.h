/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ------------------------------------------------------
   This file is part of the USER-LDD package for LAMMPS.
   Contributed by Michael R. DeLyser, mrd5285@psu.edu
   The Pennsylvania State University
   ------------------------------------------------------ */

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(ldd,AtomVecLdd);
// clang-format on
#else

#ifndef LMP_ATOM_VEC_LDD_H
#define LMP_ATOM_VEC_LDD_H

#include "atom_vec.h"

namespace LAMMPS_NS {

class AtomVecLdd : public AtomVec {
 public:
  AtomVecLdd(class LAMMPS *);

  void grow_pointers();
  void process_args(int, char**);
 private:
  double **ldd_local_density;
  double **ldd_energy;
  double **ldd_grad_density;
  double **ldd_grad_energy;
  double *ldd_total_energy;
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:


*/
