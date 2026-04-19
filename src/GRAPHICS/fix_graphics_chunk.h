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
FixStyle(graphics/chunk,FixGraphicsChunk);
// clang-format on
#else

#ifndef LMP_FIX_GRAPHICS_CHUNK_H
#define LMP_FIX_GRAPHICS_CHUNK_H

#include "fix.h"

namespace LAMMPS_NS {
class ComputeChunkAtom;

class FixGraphicsChunk : public Fix {
 public:
  FixGraphicsChunk(class LAMMPS *, int, char **);
  ~FixGraphicsChunk() override;

  int setmask() override;
  void init() override;
  void end_of_step() override;

  int image(int *&, double **&) override;

 private:
  double radius;             // global atom radius for hull inflation
  double alpha;              // triangulation algorithm parameter determining curvature
  int quality;               // icosahedron refinement level (0=12pts, 1=42pts, 2=162pts)
  bool has_global_radius;    // override auto-determined radius with global value
  bool smooth;               // smooth vs flat shading
  char *id_chunk;            // compute chunk/atom ID
  ComputeChunkAtom *cchunk;

  int numobjs;
  int *imgobjs;
  double **imgparms;
};
}    // namespace LAMMPS_NS
#endif
#endif
