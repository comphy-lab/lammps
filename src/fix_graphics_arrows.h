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

#ifdef FIX_CLASS
// clang-format off
FixStyle(graphics/arrows,FixGraphicsArrows);
// clang-format on
#else

#ifndef LMP_FIX_GRAPHICS_ARROWS_H
#define LMP_FIX_GRAPHICS_ARROWS_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGraphicsArrows : public Fix {
 public:
  FixGraphicsArrows(class LAMMPS *, int, char **);
  ~FixGraphicsArrows() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;

  int image(int *&, double **&) override;

 protected:
  int mode;
  int varflag;
  double scale;
  double radius;
  double trans;
  double val[3];
  char *xstr;
  char *ystr;
  char *zstr;
  int xvar, yvar, zvar;
  char *id_chunk;
  char *id_pos;
  char *id_vec;

  int numobjs;
  int *imgobjs;
  double **imgparms;
};
}    // namespace LAMMPS_NS
#endif
#endif
