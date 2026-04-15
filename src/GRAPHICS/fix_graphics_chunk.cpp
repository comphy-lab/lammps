/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include "fix_graphics_chunk.h"

#include "atom.h"
#include "comm.h"
#include "compute_chunk_atom.h"
#include "domain.h"
#include "error.h"
#include "graphics.h"
#include "image_objects.h"
#include "memory.h"
#include "modify.h"
#include "update.h"

#include <cstring>
#include <vector>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGraphicsChunk::FixGraphicsChunk(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), id_chunk(nullptr), cchunk(nullptr), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 5) utils::missing_cmd_args(FLERR, "fix graphics/chunk", error);

  // parse mandatory args

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Fix graphics/chunk nevery value {} must be > 0", nevery);
  global_freq = nevery;
  dynamic_group_allow = 1;

  id_chunk = utils::strdup(arg[4]);
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(id_chunk));
  if (!cchunk)
    error->all(FLERR, 4, "Chunk/atom compute {} does not exist or is incorrect style for fix {}",
               id_chunk, style);

  // defaults
  numobjs = 0;
  radius = 0.0;
  smooth = true;

  // parse optional args

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "radius") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/chunk radius", error);
      radius = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (radius < 0.0)
        error->all(FLERR, iarg + 1, "Fix graphics/chunk radius value must be >= 0");
      iarg += 2;
    } else if (strcmp(arg[iarg], "shading") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/chunk shading", error);
      if (strcmp(arg[iarg + 1], "smooth") == 0) {
        smooth = true;
      } else if (strcmp(arg[iarg + 1], "flat") == 0) {
        smooth = false;
      } else {
        error->all(FLERR, iarg + 1, "Unknown fix graphics/chunk shading setting {}", arg[iarg + 1]);
      }
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix graphics/chunk keyword {}", arg[iarg]);
    }
  }

  if (domain->dimension == 2)
    error->all(FLERR, "Fix graphics/chunk is currently not compatible with 2d systems");
}

/* ---------------------------------------------------------------------- */

FixGraphicsChunk::~FixGraphicsChunk()
{
  delete[] id_chunk;
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsChunk::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsChunk::init()
{
  cchunk = dynamic_cast<ComputeChunkAtom *>(modify->get_compute_by_id(id_chunk));
  if (!cchunk)
    error->all(FLERR, Error::NOLASTLINE,
               "Chunk/atom compute {} does not exist or is incorrect style for fix {}", id_chunk,
               style);

  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsChunk::end_of_step()
{
  using ImageObjects::vec3;

  memory->destroy(imgobjs);
  memory->destroy(imgparms);
  imgobjs = nullptr;
  imgparms = nullptr;
  numobjs = 0;

  // invoke chunk/atom compute and get chunk assignments

  modify->clearstep_compute();

  int nchunk = cchunk->setup_chunks();
  cchunk->compute_ichunk();

  const int nlocal = atom->nlocal;
  const int *const mask = atom->mask;
  const int *const type = atom->type;
  const double *const *const x = atom->x;
  const int *const ichunk = cchunk->ichunk;

  // determine per-atom radius: use per-atom property if available, else use fixed radius

  double atom_radius = radius;
  bool has_peratom_radius = (atom->radius != nullptr);

  // gather points per chunk along with their atom types
  // chunks are 1-based, ichunk[i] == 0 means atom is not in any chunk

  struct AtomInfo {
    vec3 pos;
    int atype;
    double aradius;
  };
  std::vector<std::vector<AtomInfo>> chunk_atoms(nchunk);

  for (int i = 0; i < nlocal; ++i) {
    if (!(mask[i] & groupbit)) continue;
    int ic = ichunk[i];
    if (ic < 1 || ic > nchunk) continue;
    double r = has_peratom_radius ? atom->radius[i] : atom_radius;
    chunk_atoms[ic - 1].push_back({{x[i][0], x[i][1], x[i][2]}, type[i], r});
  }

  // build convex hulls for each chunk and collect all TRINORM objects

  struct TriData {
    int objtype;
    double type0, type1, type2;
    vec3 v0, v1, v2;
    vec3 n0, n1, n2;
  };
  std::vector<TriData> all_tris;

  ImageObjects::ConvexHullObj hull;

  for (int c = 0; c < nchunk; ++c) {
    const auto &atoms = chunk_atoms[c];
    if (atoms.empty()) continue;

    // collect positions and determine effective radius
    std::vector<vec3> pts;
    pts.reserve(atoms.size());
    double max_radius = 0.0;
    for (const auto &ai : atoms) {
      pts.push_back(ai.pos);
      if (ai.aradius > max_radius) max_radius = ai.aradius;
    }

    // build convex hull with radius inflation
    hull.build(pts, max_radius, smooth);

    const auto &tris = hull.get_triangles();
    const auto &norms = hull.get_normals();
    const auto &cidx = hull.get_color_indices();

    // map color indices to atom types
    for (size_t t = 0; t < tris.size(); ++t) {
      TriData td;
      td.objtype = Graphics::TRINORM;

      // get atom type for each vertex color based on closest atom index
      int ci0 = (cidx[t][0] >= 0 && cidx[t][0] < (int) atoms.size()) ? cidx[t][0] : 0;
      int ci1 = (cidx[t][1] >= 0 && cidx[t][1] < (int) atoms.size()) ? cidx[t][1] : 0;
      int ci2 = (cidx[t][2] >= 0 && cidx[t][2] < (int) atoms.size()) ? cidx[t][2] : 0;

      td.type0 = atoms[ci0].atype;
      td.type1 = atoms[ci1].atype;
      td.type2 = atoms[ci2].atype;
      td.v0 = tris[t][0];
      td.v1 = tris[t][1];
      td.v2 = tris[t][2];
      td.n0 = norms[t][0];
      td.n1 = norms[t][1];
      td.n2 = norms[t][2];
      all_tris.push_back(td);
    }
  }

  // allocate and fill imgobjs and imgparms arrays

  numobjs = static_cast<int>(all_tris.size());
  if (numobjs > 0) {
    memory->create(imgobjs, numobjs, "fix_graphics_chunk:imgobjs");
    memory->create(imgparms, numobjs, 21, "fix_graphics_chunk:imgparms");

    for (int n = 0; n < numobjs; ++n) {
      const auto &td = all_tris[n];
      imgobjs[n] = td.objtype;
      imgparms[n][0] = td.type0;
      imgparms[n][1] = td.type1;
      imgparms[n][2] = td.type2;
      imgparms[n][3] = td.v0[0];
      imgparms[n][4] = td.v0[1];
      imgparms[n][5] = td.v0[2];
      imgparms[n][6] = td.v1[0];
      imgparms[n][7] = td.v1[1];
      imgparms[n][8] = td.v1[2];
      imgparms[n][9] = td.v2[0];
      imgparms[n][10] = td.v2[1];
      imgparms[n][11] = td.v2[2];
      imgparms[n][12] = td.n0[0];
      imgparms[n][13] = td.n0[1];
      imgparms[n][14] = td.n0[2];
      imgparms[n][15] = td.n1[0];
      imgparms[n][16] = td.n1[1];
      imgparms[n][17] = td.n1[2];
      imgparms[n][18] = td.n2[0];
      imgparms[n][19] = td.n2[1];
      imgparms[n][20] = td.n2[2];
    }
  }

  modify->addstep_compute((update->ntimestep / nevery) * nevery + nevery);
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsChunk::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
