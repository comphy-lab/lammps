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

#include "fix_oxdna_lrf.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

FixOxdnaLRF::FixOxdnaLRF(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), nx(nullptr), ny(nullptr), nz(nullptr)
{
  comm_forward = 9;

  int nmax = atom->nmax;
  memory->create(nx, nmax, 3, "FixOxdnaLRF:nx");
  memory->create(ny, nmax, 3, "FixOxdnaLRF:ny");
  memory->create(nz, nmax, 3, "FixOxdnaLRF:nz");
  atom->add_callback(0);
}

/* ---------------------------------------------------------------------- */

FixOxdnaLRF::~FixOxdnaLRF()
{
  atom->delete_callback(id, 0);
  memory->destroy(nx);
  memory->destroy(ny);
  memory->destroy(nz);
}

/* ---------------------------------------------------------------------- */

int FixOxdnaLRF::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::init()
{
  // assume force->pair, since fix is initialised through pair oxdna/excv
  auto req = neighbor->add_request(this, NeighConst::REQ_DEFAULT);
  req->set_cutoff(force->pair->cutforce);

  double cutghost;    // as computed by Neighbor and Comm
  cutghost = MAX(force->pair->cutforce + neighbor->skin, comm->cutghostuser);

  if (force->pair->cutforce > cutghost) comm->cutghostuser = force->pair->cutforce + neighbor->skin;
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::setup_pre_force(int vflag)
{
  pre_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::pre_force(int /*vflag*/)
{
  compute_lrf();
}

/* ---------------------------------------------------------------------- */

double FixOxdnaLRF::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 0.0;
  bytes += nmax * 9 * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::grow_arrays(int nmax)
{
  memory->grow(nx, nmax, 3, "FixOxdnaLRF:nx");
  memory->grow(ny, nmax, 3, "FixOxdnaLRF:ny");
  memory->grow(nz, nmax, 3, "FixOxdnaLRF:nz");
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::copy_arrays(int i, int j, int delflag)
{
  memcpy(nx[j], nx[i], sizeof(double) * 3);
  memcpy(ny[j], ny[i], sizeof(double) * 3);
  memcpy(nz[j], nz[i], sizeof(double) * 3);
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::set_arrays(int i)
{
  memset(nx[i], 0, sizeof(double) * 3);
  memset(ny[i], 0, sizeof(double) * 3);
  memset(nz[i], 0, sizeof(double) * 3);
}

/* ---------------------------------------------------------------------- */

int FixOxdnaLRF::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = nx[j][0];
    buf[m++] = nx[j][1];
    buf[m++] = nx[j][2];
    buf[m++] = ny[j][0];
    buf[m++] = ny[j][1];
    buf[m++] = ny[j][2];
    buf[m++] = nz[j][0];
    buf[m++] = nz[j][1];
    buf[m++] = nz[j][2];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    nx[i][0] = buf[m++];
    nx[i][1] = buf[m++];
    nx[i][2] = buf[m++];
    ny[i][0] = buf[m++];
    ny[i][1] = buf[m++];
    ny[i][2] = buf[m++];
    nz[i][0] = buf[m++];
    nz[i][1] = buf[m++];
    nz[i][2] = buf[m++];
  }
}

/* ---------------------------------------------------------------------- */

void FixOxdnaLRF::compute_lrf()
{
  using namespace MathExtra;    // q_to_exyz

  int nlocal = atom->nlocal;
  auto *avec = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::Bonus *bonus = avec->bonus;
  int *ellipsoid = atom->ellipsoid;

  int inum = list->inum;
  int *ilist = list->ilist;

  // loop over all local atoms, calculation of local reference frame
  for (int i = 0; i < nlocal; ++i) {
    if (atom->mask[i] & groupbit) {

      int n = ilist[i];
      // quaternion and Cartesian unit vectors in lab frame
      double *qn, nx_temp[3], ny_temp[3], nz_temp[3];

      if (ellipsoid[n] < 0) continue;    // skip non-ellipsoid atoms
      qn = bonus[ellipsoid[n]].quat;

      q_to_exyz(qn, nx_temp, ny_temp, nz_temp);

      nx[n][0] = nx_temp[0];
      nx[n][1] = nx_temp[1];
      nx[n][2] = nx_temp[2];
      ny[n][0] = ny_temp[0];
      ny[n][1] = ny_temp[1];
      ny[n][2] = ny_temp[2];
      nz[n][0] = nz_temp[0];
      nz[n][1] = nz_temp[1];
      nz[n][2] = nz_temp[2];
    }
  }

  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

void *FixOxdnaLRF::extract(const char *str, int &dim)
{
  dim = 2;

  if (strcmp(str, "nx") == 0) return (void *) nx;
  if (strcmp(str, "ny") == 0) return (void *) ny;
  if (strcmp(str, "nz") == 0) return (void *) nz;
  return nullptr;
}
