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

#ifdef PAIR_CLASS
// clang-format off
PairStyle(oxdna3/coaxstk,PairOxdna3Coaxstk);
// clang-format on
#else

#ifndef LMP_PAIR_OXDNA3_COAXSTK_H
#define LMP_PAIR_OXDNA3_COAXSTK_H

#include "pair_oxdna2_coaxstk.h"

namespace LAMMPS_NS {

class PairOxdna3Coaxstk : public PairOxdna2Coaxstk {
 public:
  PairOxdna3Coaxstk(class LAMMPS *lmp);
  void compute_stacking_site(double *, double *, double *, double *) const override;
};

}    // namespace LAMMPS_NS

#endif
#endif
