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
FixStyle(graphics/labels,FixGraphicsLabels);
// clang-format on
#else

#ifndef LMP_FIX_GRAPHICS_LABELS_H
#define LMP_FIX_GRAPHICS_LABELS_H

#include "fix.h"

namespace LAMMPS_NS {
class Compute;
class ComputeChunkAtom;

class FixGraphicsLabels : public Fix {
 public:
  FixGraphicsLabels(class LAMMPS *, int, char **);
  ~FixGraphicsLabels() override;
  int setmask() override;
  void init() override;
  void setup(int) override;
  void end_of_step() override;

  int image(int *&, double **&) override;

 protected:
  struct PixmapInfo {
    double pos[3];
    int width;
    int height;
    unsigned char *pixmap;
    double transcolor[3];
    double scale;
    int xvar, yvar, zvar, svar;
    char *xstr, *ystr, *zstr, *sstr;
  };
  std::vector<PixmapInfo> pixmaps;

  struct TextInfo {
    std::string text;
    double pos[3];
    int width;
    int height;
    unsigned char *pixmap;
    double fontcolor[3];
    double backcolor[3];
    double framecolor[3];
    double transcolor[3];
    double size;
    int xvar, yvar, zvar, svar;
    char *xstr, *ystr, *zstr, *sstr;
  };
  std::vector<TextInfo> texts;

  int varflag;
  int numobjs;
  int *imgobjs;
  double **imgparms;
};
}    // namespace LAMMPS_NS
#endif
#endif
