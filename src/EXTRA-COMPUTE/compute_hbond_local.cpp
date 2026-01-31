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

#include "compute_hbond_local.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "graphics.h"
#include "group.h"
#include "math_const.h"
#include "math_extra.h"
#include "memory.h"
#include "molecule.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

#include <algorithm>
#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;
using MathConst::DEG2RAD;
using MathConst::MY_PI2;
using MathConst::RAD2DEG;

static constexpr int DELTA = 10000;
static constexpr double EPSILON = 1.0e-10;

enum { DONOR = 0, ACCEPTOR, HYDROGEN, DIST, ANGLE, HDIST, ENGPOT, FORCE, MAXVAL };

/* ---------------------------------------------------------------------- */

ComputeHBondLocal::ComputeHBondLocal(LAMMPS *lmp, int narg, char **arg) :
    Compute(lmp, narg, arg), list(nullptr), alocal(nullptr), imgobjs(nullptr), imgparms(nullptr)
{
  if (atom->molecular == Atom::ATOMIC)
    error->all(FLERR, "Cannot (yet) use compute hbond/local with non-molecular system");
  if (atom->map_style == Atom::MAP_NONE)
    error->all(FLERR, "Compute hbond/local requires an atom map, see atom_modify");

  if (narg < 8) utils::missing_cmd_args(FLERR, "compute hbond/local", error);

  local_flag = scalar_flag = image_flag = 1;
  extscalar = 1;

  ncount = singleflag = hdistflag = numobjs = 0;
  hydrogenmask = donormask = acceptormask = 0;

  distcutoff = utils::numeric(FLERR, arg[3], false, lmp);
  if (distcutoff <= 0.0) error->all(FLERR, 3, "Compute hbond/local distance cutoff must be > 0.0");
  distcutoffsq = distcutoff * distcutoff;
  anglecutoff = DEG2RAD * utils::numeric(FLERR, arg[4], false, lmp);
  if (anglecutoff <= 0.0) error->all(FLERR, 4, "Compute hbond/local angle cutoff must be > 0.0");
  if (anglecutoff > MY_PI2)
    error->all(FLERR, 4, "Compute hbond/local angle cutoff must not be > 90.0 degrees");

  donormask = group->get_bitmask_by_id(FLERR, arg[5], "compute hbond/local donor");
  acceptormask = group->get_bitmask_by_id(FLERR, arg[6], "compute hbond/local acceptor");
  hydrogenmask = group->get_bitmask_by_id(FLERR, arg[7], "compute hbond/local hydrogen");

  // fist three elements of a local vector row are always set
  vflag.resize(MAXVAL);
  vflag[DONOR] = DONOR;
  vflag[ACCEPTOR] = ACCEPTOR;
  vflag[HYDROGEN] = HYDROGEN;

  int nvalues = 3;    // always store 3 atom IDs
  for (int iarg = 8; iarg < narg; iarg++) {
    if (strcmp(arg[iarg], "dist") == 0) {
      vflag[nvalues++] = DIST;
    } else if (strcmp(arg[iarg], "angle") == 0) {
      vflag[nvalues++] = ANGLE;
    } else if (strcmp(arg[iarg], "hdist") == 0) {
      hdistflag = 1;
      vflag[nvalues++] = HDIST;
    } else if (strcmp(arg[iarg], "engpot") == 0) {
      hdistflag = 1;
      singleflag = 1;
      vflag[nvalues++] = ENGPOT;
    } else if (strcmp(arg[iarg], "force") == 0) {
      hdistflag = 1;
      singleflag = 1;
      vflag[nvalues++] = FORCE;
    } else {
      error->all(FLERR, iarg, "Unknown compute hbond/local property {}", arg[iarg]);
    }
  }
  // reset to actual size
  vflag.resize(nvalues);

  if (singleflag && (!force->pair || !force->pair->single_enable))
    error->all(FLERR, "Computation of hydrogen bond energy or force not supported by pair style");

  // initialize output settings

  size_local_cols = nvalues;
  nmax = -1;
  alocal = nullptr;
}

/* ---------------------------------------------------------------------- */

ComputeHBondLocal::~ComputeHBondLocal()
{
  memory->destroy(alocal);
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

void ComputeHBondLocal::init()
{
  // need an occasional full neighbor list

  if (neighbor->style == Neighbor::MULTI)
    error->all(FLERR, Error::NOLASTLINE,
               "Compute hbond/local requires neighbor style 'bin' or 'nsq'");
  auto *req = neighbor->add_request(this, NeighConst::REQ_FULL | NeighConst::REQ_OCCASIONAL);
  req->set_cutoff(distcutoff);

  // do initial memory allocation assuming all donors have two hydrogen bonds

  ncount = 0;
  const auto *const mask = atom->mask;
  const auto nlocal = atom->nlocal;
  for (int i = 0; i < nlocal; ++i) {
    if ((mask[i] & groupbit) && (mask[i] & donormask)) ncount += 2;
  }

  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
}

/* ---------------------------------------------------------------------- */

void ComputeHBondLocal::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

double ComputeHBondLocal::compute_scalar()
{
  invoked_scalar = update->ntimestep;
  if (invoked_local != update->ntimestep) compute_local();
  double local_scalar = ncount;
  MPI_Allreduce(&local_scalar, &scalar, 1, MPI_DOUBLE, MPI_SUM, world);
  return scalar;
}

/* ---------------------------------------------------------------------- */

void ComputeHBondLocal::compute_local()
{
  invoked_local = update->ntimestep;
  neighbor->build_one(list);

  // count local entries and compute bond info

  ncount = compute_hbonds(0);
  if (ncount > nmax) reallocate(ncount);
  size_local_rows = ncount;
  ncount = compute_hbonds(1);
}

/* ---------------------------------------------------------------------- */

int ComputeHBondLocal::compute_image(int *&objs, double **&parms)
{
  if (invoked_image != update->ntimestep) {
    invoked_image = update->ntimestep;
    if (invoked_local != update->ntimestep) compute_local();

    memory->destroy(imgobjs);
    memory->destroy(imgparms);
    memory->create(imgobjs, ncount, "hbond/local:imgobjs");
    memory->create(imgparms, ncount, 10, "hbond/local:imgparms");

    const auto *const *const x = atom->x;
    const auto *const type = atom->type;
    double mid[3], vec[3];
    numobjs = 0;
    for (int i = 0; i < ncount; ++i) {
      int idonor = atom->map(alocal[i][DONOR]);    // always a local atom
      int iacceptor = domain->closest_image(idonor, atom->map(alocal[i][ACCEPTOR]));
      int ihydrogen = domain->closest_image(idonor, atom->map(alocal[i][HYDROGEN]));
      // this should not happen, but better safe than sorry and create a crash
      if ((idonor < 0) || (iacceptor < 0) || (ihydrogen < 0)) continue;

      imgobjs[numobjs] = Graphics::ARROW;
      imgparms[numobjs][0] = type[idonor];
      imgparms[numobjs][8] = 0.0;
      imgparms[numobjs][9] = 0.2;
      MathExtra::add3(x[iacceptor], x[ihydrogen], vec);
      MathExtra::scale3(0.5, vec, mid);
      MathExtra::sub3(x[iacceptor], x[ihydrogen], vec);
      imgparms[numobjs][7] = MathExtra::len3(vec);
      MathExtra::norm3(vec);
      imgparms[numobjs][1] = mid[0];
      imgparms[numobjs][2] = mid[1];
      imgparms[numobjs][3] = mid[2];
      imgparms[numobjs][4] = vec[0];
      imgparms[numobjs][5] = vec[1];
      imgparms[numobjs][6] = vec[2];
      ++numobjs;
    }
  }
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
/* ----------------------------------------------------------------------
   determine hydrogen bonds and compute hbond info on this proc
------------------------------------------------------------------------- */

int ComputeHBondLocal::compute_hbonds(int flag)
{
  int nhb = 0;

  auto inum = list->inum;
  auto *ilist = list->ilist;
  auto *numneigh = list->numneigh;
  auto **firstneigh = list->firstneigh;

  const auto *const *const x = atom->x;
  const auto *const tag = atom->tag;
  const auto *const type = atom->type;
  const auto *const mask = atom->mask;
  const auto *const num_bond = atom->num_bond;
  const auto *const *const bond_atom = atom->bond_atom;
  const auto nlocal = atom->nlocal;
  const auto molecular = atom->molecular;
  const auto *const molindex = atom->molindex;
  const auto *const molatom = atom->molatom;
  const auto *const *const onemols = atom->avec->onemols;

  // loop over donors, then acceptors, then hydrogens attached to donor

  for (int ii = 0; ii < inum; ++ii) {
    int i = ilist[ii];
    if ((mask[i] & groupbit) && (mask[i] & donormask)) {
      int numbonds = 0;
      int imol = -1;
      int iatom = -1;
      // TODO: find bonded hydrogens for atomic systems from neighbor list
      if (molecular == Atom::MOLECULAR) {
        numbonds = num_bond[i];
      } else {
        if (molindex[i] >= 0) {
          imol = molindex[i];
          iatom = molatom[i];
          numbonds = onemols[imol]->num_bond[iatom];
        }
      }
      if (numbonds == 0) continue;

      // loop over hydrogens in group.

      for (int kk = 0; kk < numbonds; ++kk) {
        int k = -1;
        if (molecular == Atom::MOLECULAR) {
          k = atom->map(bond_atom[i][kk]);
        } else {
          auto tagprev = tag[i] - iatom - 1;
          k = atom->map(onemols[imol]->bond_atom[iatom][kk] + tagprev);
        }
        k = domain->closest_image(i, k);
        if ((k < 0) || !((mask[k] & groupbit) && (mask[k] & hydrogenmask))) continue;

        const auto xtmp = x[i][0];
        const auto ytmp = x[i][1];
        const auto ztmp = x[i][2];
        const auto *jlist = firstneigh[i];
        const auto jnum = numneigh[i];
        for (int jj = 0; jj < jnum; ++jj) {
          int j = NEIGHMASK & jlist[jj];
          if ((mask[j] & groupbit) && (mask[j] & acceptormask)) {
            double dx1 = x[j][0] - xtmp;
            double dy1 = x[j][1] - ytmp;
            double dz1 = x[j][2] - ztmp;
            double distsq = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
            if (distsq <= distcutoffsq) {
              double dx2 = x[k][0] - xtmp;
              double dy2 = x[k][1] - ytmp;
              double dz2 = x[k][2] - ztmp;

              double r1 = sqrt(distsq);
              double r2 = sqrt(dx2 * dx2 + dy2 * dy2 + dz2 * dz2);
              if ((r1 < EPSILON) || (r2 < EPSILON)) continue;

              double c = std::clamp((dx1 * dx2 + dy1 * dy2 + dz1 * dz2) / r1 * r2, -1.0, 1.0);
              double theta = acos(c);
              if (theta < anglecutoff) {
                if (flag) {
                  double fpair = 0.0;
                  double epot = 0.0;
                  double hdist = 0.0;
                  if (hdistflag) {
                    double dx = x[k][0] - x[j][0];
                    double dy = x[k][1] - x[j][1];
                    double dz = x[k][2] - x[j][2];
                    hdist = sqrt(dx * dx + dy * dy + dz * dz);
                  }
                  if (singleflag)
                    epot =
                        force->pair->single(j, k, type[j], type[k], hdist * hdist, 1.0, 1.0, fpair);

                  int m = 0;
                  for (const auto &val : vflag) {
                    switch (val) {
                      case DONOR:
                        alocal[nhb][m] = tag[i];
                        break;
                      case ACCEPTOR:
                        alocal[nhb][m] = tag[j];
                        break;
                      case HYDROGEN:
                        alocal[nhb][m] = tag[k];
                        break;
                      case DIST:
                        alocal[nhb][m] = r1;
                        break;
                      case ANGLE:
                        alocal[nhb][m] = theta * RAD2DEG;
                        break;
                      case HDIST:
                        alocal[nhb][m] = hdist;
                        break;
                      case ENGPOT:
                        alocal[nhb][m] = epot;
                        break;
                      case FORCE:
                        alocal[nhb][m] = -hdist * fpair;
                        break;
                      default:
                        alocal[nhb][m] = 0.0;
                        break;
                    }
                    ++m;
                  }
                }
                ++nhb;
              }
            }
          }
        }
      }
    }
  }
  return nhb;
}

/* ---------------------------------------------------------------------- */

void ComputeHBondLocal::reallocate(int n)
{
  // grow array_local

  while (nmax < n) nmax += DELTA;

  memory->destroy(alocal);
  memory->create(alocal, nmax, vflag.size(), "hbond/local:array_local");
  array_local = alocal;
}

/* ----------------------------------------------------------------------
   memory usage of local data
------------------------------------------------------------------------- */

double ComputeHBondLocal::memory_usage()
{
  double bytes = (double) nmax * vflag.size() * sizeof(double);
  return bytes;
}
