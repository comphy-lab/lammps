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

#include "fix_graphics_replica.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "dump_image.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "universe.h"
#include "update.h"

#include <algorithm>
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixGraphicsReplica::FixGraphicsReplica(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 4) error->universe_all(FLERR, "Too few arguments for fix graphics/replica");

  // parse mandatory arg

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->universe_all(FLERR, "Illegal fix graphics/replica nevery value");
  global_freq = nevery;
  dynamic_group_allow = 1;

  dflag = false;
  aflag = false;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "display") == 0) {
      if (iarg + 4 > narg)
        error->universe_all(FLERR, "Too few arguments for fix graphics/replica display");
      dflag = true;
      iarg += 4;
    } else if (strcmp(arg[iarg], "average") == 0) {
      if (iarg + 4 > narg)
        error->universe_all(FLERR, "Too few arguments for fix graphics/replica average");
      aflag = true;
      iarg += 4;
    } else {
      error->universe_all(FLERR, std::string("Unknown fix graphics/replica keyword: ") + arg[iarg]);
    }
  }

  // require atom map to sort atoms by ID

  if (atom->map_style == Atom::MAP_NONE)
    error->universe_all(FLERR, "Fix graphics/replica requires an atom map, see atom_modify");

  numobjs = 0;
}

/* ---------------------------------------------------------------------- */

FixGraphicsReplica::~FixGraphicsReplica()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsReplica::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsReplica::init()
{
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsReplica::end_of_step()
{
  const int me = comm->me;
  const int nprocs = comm->nprocs;

  memory->destroy(imgobjs);
  memory->destroy(imgparms);

  // count atoms in group and across replica

  bigint nper = 0;
  bigint nall = 0;
  if (me == 0) nper = group->count(igroup);
  MPI_Allreduce(&nper, &nall, 1, MPI_LMP_BIGINT, MPI_SUM, universe->uworld);

  // ensure the group has the same number of atoms on each replica

  bigint nminmax = 0;
  MPI_Allreduce(&nper, &nminmax, 1, MPI_LMP_BIGINT, MPI_MIN, universe->uworld);
  if (nminmax != nper)
    error->universe_all(FLERR, "Fix group must have the same number of atoms for each replica");
  MPI_Allreduce(&nper, &nminmax, 1, MPI_LMP_BIGINT, MPI_MAX, universe->uworld);
  if (nminmax != nper)
    error->universe_all(FLERR, "Fix group must have the same number of atoms for each replica");

  // determine number of spheres to draw and check for overflow

  bigint numtotal = 0;
  if (dflag) numtotal += nall;
  if (aflag) numtotal += nper;
  if (numtotal >= MAXSMALLINT) error->universe_all(FLERR, "Too many graphics objects");
  numobjs = (int) numtotal;

  // create sorted map of atom-IDs

  const auto *const *const x = atom->x;
  const auto *const tag = atom->tag;
  const auto *const mask = atom->mask;
  const auto *const type = atom->type;
  const auto *const image = atom->image;
  const auto nlocal = atom->nlocal;

  std::vector<int> recvcounts(nprocs, 0);
  std::vector<int> displs(nprocs, 0);
  int ngroup = 0;
  for (int i = 0; i < nlocal; ++i)
    if (mask[i] & groupbit) ++ngroup;
  MPI_Allgather(&ngroup, 1, MPI_INT, recvcounts.data(), 1, MPI_INT, world);
  for (int i = 1; i < nprocs; ++i) displs[i] = displs[i - 1] + recvcounts[i - 1];

  std::vector<tagint> tagme;
  std::vector<tagint> taglist(nper, 0);
  for (int i = 0; i < nlocal; ++i)
    if (mask[i] & groupbit) tagme.emplace_back(tag[i]);
  MPI_Allgatherv(tagme.data(), tagme.size(), MPI_LMP_TAGINT, taglist.data(), recvcounts.data(),
                 displs.data(), MPI_LMP_TAGINT, world);
  std::sort(taglist.begin(), taglist.end());

  // collect unwrapped positions and atom type data
  // for selected atoms sorted by ID on the root of each replica

  std::vector<double> coords(3 * nper, 0.0);
  std::vector<int> types(nper, 0);
  int n = 0;
  for (const auto id : taglist) {
    int i = atom->map(id);
    if ((i >= 0) && (i < nlocal)) {
      double tmp[3] = {x[i][0], x[i][1], x[i][2]};
      domain->unmap(tmp, image[i]);
      coords[3 * n] = tmp[0];
      coords[3 * n + 1] = tmp[1];
      coords[3 * n + 2] = tmp[2];
      types[n] = type[i];
    }
    ++n;
  }

  MPI_Reduce(MPI_IN_PLACE, types.data(), 3 * nper, MPI_DOUBLE, MPI_SUM, 0, world);
  MPI_Reduce(MPI_IN_PLACE, coords.data(), 3 * nper, MPI_DOUBLE, MPI_SUM, 0, world);

  // now we are ready to create the graphics items
  // only universe root creates the objects

  if (universe->me == 0) {
    memory->create(imgobjs, numobjs, "fix_graphics:imgobjs");
    memory->create(imgparms, numobjs, 5, "fix_graphics:imgparms");
    int n = 0;
    // first process our own data;
    std::vector<double> buf(coords);
    if (dflag) {
      for (int i = 0; i < nper; ++i) {
        imgobjs[n] = DumpImage::SPHERE;
        imgparms[n][0] = type[i];
        domain->remap(buf.data() + 3 * i);
        imgparms[n][1] = buf[3 * i];
        imgparms[n][2] = buf[3 * i + 1];
        imgparms[n][3] = buf[3 * i + 2];
        imgparms[n][4] = dradius;
        imgparms[n][5] = dtrans;
        ++n;
      }
    }

    // now get data from other replicas

    for (int j = 1; j < universe->nworlds; ++j) {
      MPI_Recv(buf.data(), 3 * nper, MPI_DOUBLE, MPI_ANY_SOURCE, 0, universe->uworld,
               MPI_STATUS_IGNORE);
      if (dflag) {
        for (int i = 0; i < nper; ++i) {
          imgobjs[n] = DumpImage::SPHERE;
          imgparms[n][0] = type[i];
          domain->remap(buf.data() + 3 * i);
          imgparms[n][1] = buf[3 * i];
          imgparms[n][2] = buf[3 * i + 1];
          imgparms[n][3] = buf[3 * i + 2];
          imgparms[n][4] = dradius;
          imgparms[n][5] = dtrans;
          ++n;
        }
      }
    }
    fprintf(stderr, "output %d objects of %d\n", n, numobjs);
  } else {
    if (me == 0) MPI_Send(coords.data(), 3 * nper, MPI_DOUBLE, 0, 0, universe->uworld);
  }
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image on universe root only
------------------------------------------------------------------------- */

int FixGraphicsReplica::image(int *&objs, double **&parms)
{
  if (universe->me == 0) {
    objs = imgobjs;
    parms = imgparms;
    return numobjs;
  }
  return 0;
}
