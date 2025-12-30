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
FixStyle(graphics/surface,FixGraphicsSurface);
// clang-format on
#else

#ifndef LMP_FIX_GRAPHICS_SURFACE_H
#define LMP_FIX_GRAPHICS_SURFACE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixGraphicsSurface : public Fix {
 public:
  FixGraphicsSurface(class LAMMPS *, int, char **);
  ~FixGraphicsSurface() override;
  int setmask() override;
  void init() override;
  void end_of_step() override;

  int image(int *&, double **&) override;

 private:
  int varflag;
  int atype;
  int quality;
  int binary;
  double iso;
  double rad;
  char *rstr;
  int rvar;
  int pad;
  std::string filename;

  double dx, dy, dz;
  int nx, ny, nz;

  int numobjs;
  int *imgobjs;
  double **imgparms;
};
}    // namespace LAMMPS_NS
#endif
#endif
