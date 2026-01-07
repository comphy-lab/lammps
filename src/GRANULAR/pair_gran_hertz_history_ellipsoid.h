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
PairStyle(gran/hertz/history/ellipsoid,PairGranHertzHistoryEllipsoid);
// clang-format on
#else

#ifndef LMP_PAIR_GRAN_HERTZ_HISTORY_ELLIPSOID_H
#define LMP_PAIR_GRAN_HERTZ_HISTORY_ELLIPSOID_H

#include "pair_gran_hooke_history_ellipsoid.h"

namespace LAMMPS_NS {

class PairGranHertzHistoryEllipsoid : public PairGranHookeHistoryEllipsoid {
 public:
  PairGranHertzHistoryEllipsoid(class LAMMPS *);
  void compute(int, int) override;
  void settings(int, char **) override;
  double single(int, int, int, int, double, double, double, double &) override;

 protected:
  int curvature_model; 
  enum { CURV_MEAN, CURV_GAUSSIAN };
};

}    // namespace LAMMPS_NS

#endif
#endif
