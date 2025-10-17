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
/* ------------------------------------------------------
   This file is part of the LDD package for LAMMPS.
   Contributed by Michael R. DeLyser, mrd5285@psu.edu
   and Maria C. Lesniewski, mjl6766@psu.edu
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

  void grow_pointers() override;
  void process_args(int, char**) override;
 private:
  double **ldd_local_density; // n_atom x n_atomtype matrix for LD info
  double **ldd_energy; // n_atom x n_atomtype matrix for u_{b|t_I}(rho_I)
  double **ldd_grad_density; // n_atom x 3 * n_atomtype matrix for LD grad info
  double **ldd_grad_energy; // n_atom x n_atomtype matrix for u_{\nabla b|t_I}(rho_I)
  double *ldd_total_energy; // n_atom vector, stores the sum (across all surrounding types) of u_x(rho_I) terms
};

}    // namespace LAMMPS_NS

#endif
#endif

/* ERROR/WARNING messages:


*/
