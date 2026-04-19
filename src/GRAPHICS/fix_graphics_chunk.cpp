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
#include "force.h"
#include "graphics.h"
#include "image_objects.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <map>
#include <utility>

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {

using ImageObjects::triangle;
using ImageObjects::vec3;

// Define vertices of an icosahedron to approximate a sphere of radius 1 around the origin.
// A and B are the normalized coordinates derived from the golden ratio:
// phi = (1 + sqrt(5)) / 2; A = 1 / sqrt(1 + phi^2); B = phi / sqrt(1 + phi^2)
constexpr double A = 0.5257311121191336;
constexpr double B = 0.8506508083520399;

// Base icosahedron: 12 vertices, 20 faces
const std::vector<vec3> ico_vertices = {
    {-A, B, 0.0},  {A, B, 0.0},  {-A, -B, 0.0}, {A, -B, 0.0}, {0.0, -A, B},  {0.0, A, B},
    {0.0, -A, -B}, {0.0, A, -B}, {B, 0.0, -A},  {B, 0.0, A},  {-B, 0.0, -A}, {-B, 0.0, A},
};

struct IcoFace {
  int v[3];
};
const std::vector<IcoFace> ico_faces = {
    {0, 5, 11},  {0, 1, 5},  {0, 7, 1},  {0, 10, 7}, {0, 11, 10}, {1, 9, 5},  {5, 4, 11},
    {11, 2, 10}, {10, 6, 7}, {7, 8, 1},  {3, 4, 9},  {3, 2, 4},   {3, 6, 2},  {3, 8, 6},
    {3, 9, 8},   {4, 5, 9},  {2, 11, 4}, {6, 10, 2}, {8, 7, 6},   {9, 1, 8},
};

// Generate refined icosahedron points at given level (0=12 pts, 1=42, 2=162).
// Each refinement splits every triangle into 4 sub-triangles and projects
// the new midpoint vertices onto the unit sphere.

std::vector<vec3> generate_sphere_points(int level)
{
  // Start with icosahedron vertices
  std::vector<vec3> verts = ico_vertices;
  std::vector<IcoFace> faces(ico_faces.begin(), ico_faces.end());

  for (int lev = 0; lev < level; ++lev) {
    // midpoint cache: edge -> new vertex index
    std::map<std::pair<int, int>, int> edge_map;
    std::vector<IcoFace> new_faces;
    new_faces.reserve(faces.size() * 4);

    for (const auto &f : faces) {
      int mid[3];
      for (int e = 0; e < 3; ++e) {
        int a = f.v[e], b = f.v[(e + 1) % 3];
        auto key = std::make_pair(std::min(a, b), std::max(a, b));
        auto it = edge_map.find(key);
        if (it != edge_map.end()) {
          mid[e] = it->second;
        } else {
          // Compute midpoint and project onto unit sphere
          vec3 mp = {(verts[a][0] + verts[b][0]) * 0.5, (verts[a][1] + verts[b][1]) * 0.5,
                     (verts[a][2] + verts[b][2]) * 0.5};
          double len =
              std::sqrt(mp[0] * mp[0] + mp[1] * mp[1] + mp[2] * mp[2]);
          if (len > 0.0) {
            mp[0] /= len;
            mp[1] /= len;
            mp[2] /= len;
          }
          mid[e] = static_cast<int>(verts.size());
          edge_map[key] = mid[e];
          verts.push_back(mp);
        }
      }
      new_faces.push_back({f.v[0], mid[0], mid[2]});
      new_faces.push_back({mid[0], f.v[1], mid[1]});
      new_faces.push_back({mid[2], mid[0], mid[1]});
      new_faces.push_back({mid[2], mid[1], f.v[2]});
    }
    faces = std::move(new_faces);
  }
  return verts;
}

}    // namespace

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
  alpha = 0.0;
  quality = 0;
  has_global_radius = false;
  smooth = true;

  // parse optional args

  int iarg = 5;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "radius") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/chunk radius", error);
      radius = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (radius < 0.0) error->all(FLERR, iarg + 1, "Fix graphics/chunk radius value must be >= 0");
      has_global_radius = true;
      iarg += 2;
    } else if (strcmp(arg[iarg], "alpha") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/chunk alpha", error);
      alpha = utils::numeric(FLERR, arg[iarg + 1], false, lmp);
      if (alpha < 0.0) error->all(FLERR, iarg + 1, "Fix graphics/chunk alpha value must be >= 0");
      iarg += 2;
    } else if (strcmp(arg[iarg], "quality") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/chunk quality", error);
      quality = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (quality < 0 || quality > 2)
        error->all(FLERR, iarg + 1,
                   "Fix graphics/chunk quality value must be 0, 1, or 2");
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
  const imageint *const image = atom->image;
  const double *const *const x = atom->x;
  const int *const ichunk = cchunk->ichunk;

  // determine per-atom radius: use per-atom property if available, otherwise get sigma from potential
  // else use fixed radius from input

  int dim = 0;
  double **sigma = nullptr;
  if (force->pair) {
    sigma = (double **) force->pair->extract("sigma", dim);
    if (dim != 2) sigma = nullptr;
  }

  double atom_radius = radius;
  bool has_peratom_radius = (atom->radius != nullptr);
  bool has_pertype_radius = (sigma != nullptr);

  // gather points per chunk along with their atom types
  // chunks are 1-based, ichunk[i] == 0 means atom is not in any chunk

  struct AtomInfo {
    vec3 pos;
    double rad;
    double atype;
  };
  std::vector<std::vector<AtomInfo>> chunk_atoms(nchunk);
  vec3 unwrapped;

  for (int i = 0; i < nlocal; ++i) {
    if (!(mask[i] & groupbit)) continue;
    int ic = ichunk[i];
    if (ic < 1 || ic > nchunk) continue;

    // use global atom radius if set, or else try per-atom or per-atomtype value
    if (!has_global_radius) {
      if (has_peratom_radius)
        atom_radius = atom->radius[i];
      else if (has_pertype_radius)
        atom_radius = 0.5 * sigma[type[i]][type[i]];
    }
    domain->unmap(x[i], image[i], unwrapped.data());
    chunk_atoms[ic - 1].push_back({unwrapped, atom_radius, double(type[i])});
  }

  // build hulls for each chunk and collect all TRINORM objects

  struct ObjData {
    double type0, type1, type2;
    vec3 v0, v1, v2;
    vec3 n0, n1, n2;
  };

  // Generate sphere sample points based on quality level:
  // quality 0: 12 points (icosahedron)
  // quality 1: 42 points (one refinement)
  // quality 2: 162 points (two refinements)
  const std::vector<vec3> sphere_pts = generate_sphere_points(quality);
  const int nsp = static_cast<int>(sphere_pts.size());

  std::vector<ObjData> all_objs;
  ImageObjects::ConvexHullObj hull;
  for (int c = 0; c < nchunk; ++c) {
    const auto &iatoms = chunk_atoms[c];
    const auto natoms = iatoms.size();
    if (iatoms.empty()) continue;

    vec3 center{0.0, 0.0, 0.0};

    for (const auto &ai : iatoms) {
      center[0] += ai.pos[0];
      center[1] += ai.pos[1];
      center[2] += ai.pos[2];
    }
    center[0] /= double(natoms);
    center[1] /= double(natoms);
    center[2] /= double(natoms);

    vec3 wrapped{center};
    domain->remap(wrapped.data());
    vec3 offset{center[0] - wrapped[0], center[1] - wrapped[1], center[2] - wrapped[2]};

    // replace atom positions with sphere sample points scaled to radius
    std::vector<vec3> pts;
    pts.reserve(nsp * natoms);
    for (const auto &ai : iatoms) {
      for (const auto &sp : sphere_pts) {
        pts.push_back({ai.rad * sp[0] + ai.pos[0] - offset[0],
                       ai.rad * sp[1] + ai.pos[1] - offset[1],
                       ai.rad * sp[2] + ai.pos[2] - offset[2]});
      }
    }

    // build hull
    hull.build(pts, smooth, alpha);

    const auto &tris = hull.get_triangles();
    const auto &norms = hull.get_normals();
    const auto &cidx = hull.get_color_indices();

    for (size_t t = 0; t < tris.size(); ++t) {

      // get index into original atoms array for access to atom type
      int ci0 = cidx[t][0] / nsp;
      int ci1 = cidx[t][1] / nsp;
      int ci2 = cidx[t][2] / nsp;
      if ((ci0 < 0) || (ci0 >= (int) iatoms.size())) ci0 = 0;
      if ((ci1 < 0) || (ci1 >= (int) iatoms.size())) ci1 = 0;
      if ((ci2 < 0) || (ci2 >= (int) iatoms.size())) ci2 = 0;

      all_objs.push_back({iatoms[ci0].atype, iatoms[ci1].atype, iatoms[ci2].atype, tris[t][0],
                          tris[t][1], tris[t][2], norms[t][0], norms[t][1], norms[t][2]});
    }
  }

  // allocate and fill imgobjs and imgparms arrays

  numobjs = static_cast<int>(all_objs.size());
  if (numobjs > 0) {
    memory->create(imgobjs, numobjs, "fix_graphics_chunk:imgobjs");
    memory->create(imgparms, numobjs, 21, "fix_graphics_chunk:imgparms");

    for (int n = 0; n < numobjs; ++n) {
      const auto &od = all_objs[n];
      imgobjs[n] = Graphics::TRINORM;
      imgparms[n][0] = od.type0;
      imgparms[n][1] = od.type1;
      imgparms[n][2] = od.type2;
      imgparms[n][3] = od.v0[0];
      imgparms[n][4] = od.v0[1];
      imgparms[n][5] = od.v0[2];
      imgparms[n][6] = od.v1[0];
      imgparms[n][7] = od.v1[1];
      imgparms[n][8] = od.v1[2];
      imgparms[n][9] = od.v2[0];
      imgparms[n][10] = od.v2[1];
      imgparms[n][11] = od.v2[2];
      imgparms[n][12] = od.n0[0];
      imgparms[n][13] = od.n0[1];
      imgparms[n][14] = od.n0[2];
      imgparms[n][15] = od.n1[0];
      imgparms[n][16] = od.n1[1];
      imgparms[n][17] = od.n1[2];
      imgparms[n][18] = od.n2[0];
      imgparms[n][19] = od.n2[1];
      imgparms[n][20] = od.n2[2];
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
