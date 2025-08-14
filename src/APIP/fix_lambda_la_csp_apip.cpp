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
/* ----------------------------------------------------------------------
   Contributing author: David Immel (d.immel@fz-juelich.de, FZJ, Germany)
------------------------------------------------------------------------- */

#include "fix_lambda_la_csp_apip.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "fix_store_atom.h"
#include "force.h"
#include "group.h"
#include "memory.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "modify.h"
#include "pair.h"
#include "update.h"

#include <algorithm>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixLambdaLACSPAPIP::FixLambdaLACSPAPIP(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), ngh_pairs(nullptr), list(nullptr), distsq(nullptr), nearest(nullptr), fixstore(nullptr), f_lambda(nullptr)
{
  comm_reverse = 2;
  comm_forward = 2;

  ngh_pairs_size = 0;
  maxneigh = 0;

  peratom_flag = 1;
  peratom_freq = 1;
  size_peratom_cols = 3;

  // set defaults
  cut_lo = cut_hi = -1;
  threshold_lo = threshold_hi = -1;
  lambda_non_group = 1; // fast
  const_ngh_flag = false;

  tags_stored = false;
  counter_changed_csp_nghs = 0;
  scalar_flag = 1;
  extscalar = 1;

  if (narg < 10) error->all(FLERR, "fix lambda/la/csp/apip requires six arguments");

  threshold_lo = utils::numeric(FLERR, arg[3], false, lmp);
  threshold_hi = utils::numeric(FLERR, arg[4], false, lmp);
  threshold_width = threshold_hi - threshold_lo;

  cut_lo = utils::numeric(FLERR, arg[5], false, lmp);
  cut_hi = utils::numeric(FLERR, arg[6], false, lmp);
  cut_width = cut_hi - cut_lo;
  cut_hi_sq = cut_hi * cut_hi;

  if (strcmp(arg[7], "fcc") == 0)
    nnn = 12;
  else if (strcmp(arg[7], "bcc") == 0)
    nnn = 14;
  else
    nnn = utils::inumeric(FLERR, arg[7], false, lmp);
  csp_cutsq = pow(utils::numeric(FLERR, arg[8], false, lmp), 2);
  cutsq_combined = csp_cutsq > cut_hi_sq ? csp_cutsq : cut_hi_sq;

  if (strcmp(arg[9], "dyn_csp_ngh") == 0)
    const_ngh_flag = false;
  else if (strcmp(arg[9], "const_csp_ngh") == 0)
    const_ngh_flag = true;
  else
    error->all(FLERR, "expected dyn_csp_ngh or const_csp_ngh instead of {}", arg[9]);

  if (narg == 11) lambda_non_group = utils::numeric(FLERR, arg[10], false, lmp);

  // verify arguments
  if (cut_lo > cut_hi || cut_lo < 0) error->all(FLERR, "0 <= cut_lo <= cut_hi required");
  if (threshold_lo > threshold_hi || threshold_lo < 0) error->all(FLERR, "0 <= threshold_lo <= threshold_hi required");
  if (lambda_non_group < 0 || lambda_non_group > 1) error->all(FLERR, "0 <= lambda_non_group <= 1 required");



  if (!atom->apip_lambda_flag) { error->all(FLERR, "fix lambda/la/csp/apip requires atomic style with lambda."); }
  if (!atom->apip_la_inp_flag || !atom->apip_la_avg_flag || !atom->apip_la_norm_flag) { error->all(FLERR, "fix lambda/la/csp/apip requires atomic style with csp, csp_avg and csp_norm."); }

  size_f_lambda = atom->nlocal;
  memory->create(f_lambda, size_f_lambda, size_peratom_cols, "pair:lambda:la:csp:apip:f:lambda");
  array_atom = f_lambda;
  for (int i = 0; i < atom->nlocal; i++) {
    for (int j = 0; j < size_peratom_cols; j ++) f_lambda[i][j] = 0;
  }

}

/* ---------------------------------------------------------------------- */

FixLambdaLACSPAPIP::~FixLambdaLACSPAPIP()
{
  memory->destroy(ngh_pairs);
  memory->destroy(f_lambda);
  if (fixstore && modify->nfix) modify->delete_fix(fixstore->id);
  fixstore = nullptr;
}

/**
  * allocate per-particle storage for the last ngh_pairs
  * fix could already be allocated if fix lambda/la/csp/apip is re-specified
  */

void FixLambdaLACSPAPIP::post_constructor()
{
  std::string cmd;
  cmd = id;
  cmd += "LAST_CSP_NGH_PAIR";

  // delete existing fix store if existing
  fixstore = dynamic_cast<FixStoreAtom *>(modify->get_fix_by_id(cmd));
  // check nfix in case all fixes have already been deleted
  if (fixstore && modify->nfix) modify->delete_fix(fixstore->id);
  fixstore = nullptr;

  char str_values[40];
  sprintf(str_values, "%d", nnn);

  // arguments of peratom:
  // first: 1 -> store in restart file
  // second: number of doubles to store per atom
  cmd += " all STORE/ATOM ";
  cmd += str_values;      // n1
  cmd += " 0 0 1";                // n2 gflag rflag
  fixstore = dynamic_cast<FixStoreAtom *>(modify->add_fix(cmd));

  // carry weights with atoms during normal atom migration
  fixstore->disable = 0;
}

/* ---------------------------------------------------------------------- */

int FixLambdaLACSPAPIP::setmask()
{
  int mask = 0;
  mask |= PRE_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixLambdaLACSPAPIP::init()
{
  // full neighbour list for thermostating
  auto req = neighbor->add_request(this, NeighConst::REQ_FULL);
  req->set_cutoff(sqrt(cutsq_combined));

  if (atom->tag_enable == 0)
    error->all(FLERR, "fix lambda/la/csp/apip requires atom IDs");

  // only one fix lambda/la/csp/apip
  // TOOD update names
  int count = 0;
  for (int i = 0; i < modify->nfix; i++) {
    if (strcmp(modify->fix[i]->style, "lambda/la/csp/apip") == 0) count++;
  }
  if (count > 1) error->all(FLERR, "More than one fix lambda/la/csp/apip.");

  Pair *pair_tmp;
  pair_tmp = force->pair_match("la/csp/apip", 0);
  if (!pair_tmp) error->warning(FLERR, "fix lambda/la/csp/apip requires a `pair la/csp/apip` to apply forces");

  if (force->pair->cutforce < cut_hi)
    error->all(FLERR, "cutoff of potential ({}) smaller than cutoff of weighting function ({})", force->pair->cutforce, cut_hi);
  if (force->pair->cutforce*force->pair->cutforce < csp_cutsq)
    error->all(FLERR, "cutoff of CSP ({}) smaller than cutoff of weighting function ({})", force->pair->cutforce, sqrt(csp_cutsq));

  if (strcmp(atom->atom_style, "apip/la")) error->all(FLERR, "fix lambda/la/csp/apip requires atom style apip/la");
}

/**
  *  Init neighbor list.
  */

void FixLambdaLACSPAPIP::init_list(int /*id*/, NeighList *ptr)
{
  list = ptr;
}

/**
  *  Compute lambda, csp, csp_avg, csp_norm for all local atoms.
  */

void FixLambdaLACSPAPIP::setup_pre_force(int vflag)
{
  if (const_ngh_flag && atom->map_style == Atom::MAP_NONE) {
    atom->map_init();
    atom->map_set();
  }

  // TODO allow other mechanisms to define pairs

  if (!const_ngh_flag || (const_ngh_flag && !tags_stored)) pre_force_dyn_pairs();
  else pre_force_const_pairs();

  comm->forward_comm(this);
}

void FixLambdaLACSPAPIP::pre_force(int /*vflag*/) {
  if (const_ngh_flag) pre_force_const_pairs();
  else pre_force_dyn_pairs();
}

/**
  *  Compute lambda, csp, csp_avg, csp_norm for all local atoms with dynamic neighbour pairs.
  */

void FixLambdaLACSPAPIP::pre_force_dyn_pairs()
{
  int i, j, k, ii, jj, kk, n_nearest, n_pairs, inum, jnum;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, value, weight;
  double *lambda, *csp, *csp_avg, *csp_norm;

  tagint *tag = atom->tag;

  lambda = atom->apip_lambda;
  csp = atom->apip_la_inp;
  csp_avg = atom->apip_la_avg;
  csp_norm = atom->apip_la_norm;

  int nlocal = atom->nlocal;
  int nmax = atom->nmax;

  // grow ngh_pairs array if necessary

  if (atom->nlocal > ngh_pairs_size) {
    memory->destroy(ngh_pairs);
    ngh_pairs_size = atom->nlocal;
    memory->create(ngh_pairs, ngh_pairs_size, nnn, "fix_lambda_la_csp_apip:ngh_pairs");
  }

  double ** stored_tags = fixstore->astore;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // npairs = number of unique pairs

  int nhalf = nnn / 2;
  int npairs = nnn * (nnn - 1) / 2;
  auto pairs_value = new double[npairs];
  auto pairs_j = new int[npairs];
  auto pairs_k = new int[npairs];
  auto pairs_index = new int[npairs];
  auto pairs_used_now = new CSPpairAPIP[nhalf];
  auto pairs_used_prev = new CSPpairAPIP[nhalf];

  // compute centro-symmetry parameter for each atom in group

  double **x = atom->x;

  // zero values
  for (i = 0; i < nlocal; i++) csp[i] = 0;
  for (i = 0; i < nmax; i++) csp_norm[i] = csp_avg[i] = 0;

  counter_changed_csp_nghs = 0;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];


    // 1. compute csp[i]


    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    // ensure distsq and nearest arrays are long enough

    if (jnum > maxneigh) {
      memory->destroy(distsq);
      memory->destroy(nearest);
      maxneigh = jnum;
      memory->create(distsq, maxneigh, "centro/atom:distsq");
      memory->create(nearest, maxneigh, "centro/atom:nearest");
    }

    // loop over list of all neighbors within force cutoff
    // distsq[] = distance sq to each
    // nearest[] = atom indices of neighbors

    n_nearest = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq < cutsq_combined) {
        distsq[n_nearest] = rsq;
        nearest[n_nearest++] = j;
      }
    }

    if (n_nearest < nnn) {
      error->all(FLERR, "Atom has only {} neighbours, but {} are required to calculate the CSP.", n_nearest, nnn);
    }

    // store nnn nearest neighs in 1st nnn locations of distsq and nearest

    select2(nnn, n_nearest, distsq, nearest);

    // R = Ri + Rj for each of npairs i,j pairs among nnn neighbors
    // pairs = squared length of each R

    n_pairs = 0;
    for (j = 0; j < nnn; j++) {
      jj = nearest[j];
      for (k = j + 1; k < nnn; k++) {
        kk = nearest[k];
        delx = x[jj][0] + x[kk][0] - 2.0 * xtmp;
        dely = x[jj][1] + x[kk][1] - 2.0 * ytmp;
        delz = x[jj][2] + x[kk][2] - 2.0 * ztmp;
        pairs_value[n_pairs] = delx * delx + dely * dely + delz * delz;
        pairs_j[n_pairs] = jj;
        pairs_k[n_pairs] = kk;
        pairs_index[n_pairs++] = n_pairs;
      }
    }

    // store nhalf smallest pair distances in 1st nhalf locations of pairs

    select2(nhalf, npairs, pairs_value, pairs_index);

    // centrosymmetry = sum of nhalf smallest squared values
    // compute csp with detected neighbour pairs
    value = 0.0;
    for (j = 0; j < nhalf; j++) {
      value += pairs_value[j];
      // store used neighbor pairs
      ngh_pairs[i][j] = pairs_j[pairs_index[j]];
      ngh_pairs[i][j + nhalf] = pairs_k[pairs_index[j]];

      // get new tags
      pairs_used_now[j] =  CSPpairAPIP(tag[ngh_pairs[i][j]], tag[ngh_pairs[i][j + nhalf]]);

    }

    csp[i] = value;

    // get previous tags
    if (tags_stored) {
      for (j = 0; j < nhalf; j++) {
        pairs_used_prev[j] = CSPpairAPIP(stored_tags[i][j], stored_tags[i][j + nhalf]);
      }

      // sort tag arrays
      std::sort(pairs_used_now, pairs_used_now + nhalf);
      std::sort(pairs_used_prev, pairs_used_prev + nhalf);

      // count number of different pairs
      for (j = 0; j < nhalf; j++) {
        if (pairs_used_now[j] != pairs_used_prev[j]) {
          counter_changed_csp_nghs++;
          break;
        }
      }
    }

    // store current tags
    for (j = 0; j < nnn; j++) {
      stored_tags[i][j] = tag[ngh_pairs[i][j]];
    }

    // 2. compute csp avg
    // csp_avg_i += csp_i * w(r_ij)
    // csp_norm_i += w(r_ij)

    // own particle
    weight = weighting_function_poly(0);
    csp_avg[i] += value * weight;
    csp_norm[i] += weight;

    // neighbouring particles
    for (j = 0; j < n_nearest; j++) {
      if (distsq[j] <= cut_hi_sq) {
        weight = weighting_function_poly(sqrt(distsq[j]));
        jj = nearest[j];
        csp_avg[jj] += value * weight;
        csp_norm[jj] += weight;
      }
    }

  }
  tags_stored = true;


  // reverse communication of csp_norm and csp_avgs of ghost atoms
  comm->reverse_comm(this);

  for (i = 0; i < nlocal; i++) {
    csp_avg[i] /= csp_norm[i];
    lambda[i] = switching_function_poly(csp_avg[i]);

    // TOOD consider lambda_non_group here
  }

  delete[] pairs_value;
  delete[] pairs_index;
}

/**
  *  Compute lambda, csp, csp_avg, csp_norm for all local atoms with constant neighbour pairs.
  */

void FixLambdaLACSPAPIP::pre_force_const_pairs()
{
  int i, j, ii, jj, kk, inum, jnum;
  int *ilist, *jlist, *numneigh, **firstneigh;
  double xtmp, ytmp, ztmp, delx, dely, delz, rsq, value, weight, delx_j, dely_j, delz_j, delx_k, dely_k, delz_k;
  double *lambda, *csp, *csp_avg, *csp_norm, **x, **stored_tags;

  int nlocal = atom->nlocal;
  int nmax = atom->nmax;
  int nhalf = nnn / 2;

  tagint *tag = atom->tag;

  lambda = atom->apip_lambda;
  csp = atom->apip_la_inp;
  csp_avg = atom->apip_la_avg;
  csp_norm = atom->apip_la_norm;
  x = atom->x;
  stored_tags = fixstore->astore;

  // grow ngh_pairs array if necessary

  if (atom->nlocal > ngh_pairs_size) {
    memory->destroy(ngh_pairs);
    ngh_pairs_size = atom->nlocal;
    memory->create(ngh_pairs, ngh_pairs_size, nnn, "fix_lambda_la_csp_apip:ngh_pairs");
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;

  // zero values
  for (i = 0; i < nlocal; i++) csp[i] = 0;
  for (i = 0; i < nmax; i++) csp_norm[i] = csp_avg[i] = 0;

  // 1. compute csp[i]

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    value = 0;
    for (j = 0; j < nhalf; j++) {

      jj = atom->map((tagint) stored_tags[i][j]);
      while (jj >= 0) {
        // calculate distance to jj
        delx_j = x[jj][0] - xtmp;
        dely_j = x[jj][1] - ytmp;
        delz_j = x[jj][2] - ztmp;
        // correct periodic image?
        if (fabs(delx_j) < domain->prd_half[0] && fabs(dely_j) < domain->prd_half[1] && fabs(delz_j) < domain->prd_half[2]) {
          // correct periodic image
          break;
        } else {
          // next periodic image
          jj = atom->sametag[jj];
        }
      }
      if (jj == -1) error->one(FLERR, "atom ID {} csp neighbour with ID {} not present", atom->tag[i], stored_tags[i][j]);
      ngh_pairs[i][j] = jj;

      kk = atom->map((tagint) stored_tags[i][j + nhalf]);
      while (kk >= 0) {
        // calculate distance to kk
        delx_k = x[kk][0] - xtmp;
        dely_k = x[kk][1] - ytmp;
        delz_k = x[kk][2] - ztmp;
        // correct periodic image?
        if (fabs(delx_k) < domain->prd_half[0] && fabs(dely_k) < domain->prd_half[1] && fabs(delz_k) < domain->prd_half[2]) {
          // correct periodic image
          break;
        } else {
          // next periodic image
          kk = atom->sametag[kk];
        }
      }
      if (kk == -1) error->one(FLERR, "atom ID {} csp neighbour with ID {} not present", atom->tag[i], stored_tags[i][j + nhalf]);
      ngh_pairs[i][j + nhalf] = kk;

      delx = delx_j + delx_k;
      dely = dely_j + dely_k;
      delz = delz_j + delz_k;

      value += delx * delx + dely * dely + delz * delz;
    }
    csp[i] = value;

  }

  // 2. compute csp avg
  // csp_avg_i += csp_i * w(r_ij)
  // csp_norm_i += w(r_ij)

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];

    // own particle
    weight = weighting_function_poly(0);
    csp_avg[i] += csp[i] * weight;
    csp_norm[i] += weight;

    // neighbouring particles
    jlist = firstneigh[i];
    jnum = numneigh[i];
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      // this is a full neighbour list
      // avoid double counting
      if (i > j) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      if (rsq <= cut_hi_sq) {
        weight = weighting_function_poly(sqrt(rsq));

        csp_avg[j] += csp[i] * weight;
        csp_norm[j] += weight;

        if (j >= nlocal) continue;

        csp_norm[i] += weight;
        csp_avg[i] += csp[j] * weight;
      }
    }
  }

  // reverse communication of csp_norm and csp_avgs of ghost atoms
  comm->reverse_comm(this);

  for (i = 0; i < nlocal; i++) {
    csp_avg[i] /= csp_norm[i];
    lambda[i] = switching_function_poly(csp_avg[i]);

    // TOOD consider lambda_non_group here
  }
}

/* ----------------------------------------------------------------------
   select routine from Numerical Recipes (slightly modified)
   find k smallest values in array of length n
   and sort auxiliary array at same time
------------------------------------------------------------------------- */

void FixLambdaLACSPAPIP::select2(int k, int n, double *arr, int *iarr)
{
  int i, ir, j, l, mid, ia;
  double a;

  arr--;
  iarr--;
  l = 1;
  ir = n;
  while (true) {
    if (ir <= l + 1) {
      if (ir == l + 1 && arr[ir] < arr[l]) {
        std::swap(arr[l], arr[ir]);
        std::swap(iarr[l], iarr[ir]);
      }
      return;
    } else {
      mid = (l + ir) >> 1;
      std::swap(arr[mid], arr[l + 1]);
      std::swap(iarr[mid], iarr[l + 1]);
      if (arr[l] > arr[ir]) {
        std::swap(arr[l], arr[ir]);
        std::swap(iarr[l], iarr[ir]);
      }
      if (arr[l + 1] > arr[ir]) {
        std::swap(arr[l + 1], arr[ir]);
        std::swap(iarr[l + 1], iarr[ir]);
      }
      if (arr[l] > arr[l + 1]) {
        std::swap(arr[l], arr[l + 1]);
        std::swap(iarr[l], iarr[l + 1]);
      }
      i = l + 1;
      j = ir;
      a = arr[l + 1];
      ia = iarr[l + 1];
      while (true) {
        do i++;
        while (arr[i] < a);
        do j--;
        while (arr[j] > a);
        if (j < i) break;
        std::swap(arr[i], arr[j]);
        std::swap(iarr[i], iarr[j]);
      }
      arr[l + 1] = arr[j];
      arr[j] = a;
      iarr[l + 1] = iarr[j];
      iarr[j] = ia;
      if (j >= k) ir = j - 1;
      if (j <= k) l = i;
    }
  }
}

// helper function
// similar to cutoff_func_poly in ace_radial.cpp
// compare Phys Rev Mat 6, 013804 (2022) APPENDIX C: RADIAL AND CUTOFF FUNCTIONS 2. Cutoff function
// the first two derivatives of the switching function lambda vanishes at the boundaries of the switching region
double FixLambdaLACSPAPIP::switching_function_poly(double input)
{
  // calculate lambda
  if (input <= threshold_lo) {
    return 1;
  } else if (input >= threshold_hi) {
    return 0;
  } else {
    double tmp = 1 - 2 * (1 + (input - threshold_hi) / threshold_width);
    return 0.5 + 7.5 / 2. * (tmp / 4. - pow(tmp, 3) / 6. + pow(tmp, 5) / 20.);
  }
}

double FixLambdaLACSPAPIP::der_switching_function_poly(double input)
{
  // calculate lambda
  if (input <= threshold_lo || input >= threshold_hi) {
    return 0;
  } else {
    double tmp = 1 - 2 * (1 + (input - threshold_hi) / threshold_width);
    return -1.875 / threshold_width * (1 - 2 * tmp * tmp + pow(tmp, 4));
  }
}

double FixLambdaLACSPAPIP::weighting_function_poly(double r)
{
  // calculate lambda
  if (r <= cut_lo) {
    return 1;
  } else if (r >= cut_hi) {
    return 0;
  } else {
    double tmp = 1 - 2 * (1 + (r - cut_hi) / cut_width);
    return 0.5 + 7.5 / 2. * (tmp / 4. - pow(tmp, 3) / 6. + pow(tmp, 5) / 20.);
  }
}

double FixLambdaLACSPAPIP::der_weighting_function_poly(double r)
{
  // calculate lambda
  if (r <= cut_lo || r >= cut_hi) {
    return 0;
  } else {
    double tmp = 1 - 2 * (1 + (r - cut_hi) / cut_width);
    return -1.875 / cut_width * (1 - 2 * tmp * tmp + pow(tmp, 4));
  }
}

/**
  * extract lambda(time averaged lambda_input) and lambda_input_history_len
  */

void *FixLambdaLACSPAPIP::extract(const char *str, int &dim)
{
  dim = 2;
  if (strcmp(str, "fix_lambda_la_csp_apip:ngh_pairs") == 0) { return (void *) ngh_pairs; }
  dim = 0;
  if (strcmp(str, "fix_lambda_la_csp_apip:cut_lo") == 0) { return &cut_lo; }
  if (strcmp(str, "fix_lambda_la_csp_apip:cut_hi") == 0) { return &cut_hi; }
  if (strcmp(str, "fix_lambda_la_csp_apip:cut_hi_sq") == 0) { return &cut_hi_sq; }
  if (strcmp(str, "fix_lambda_la_csp_apip:cut_width") == 0) { return &cut_width; }
  if (strcmp(str, "fix_lambda_la_csp_apip:lambda_non_group") == 0) { return &lambda_non_group; }
  if (strcmp(str, "fix_lambda_la_csp_apip:cutsq_combined") == 0) { return &cutsq_combined; }
  if (strcmp(str, "fix_lambda_la_csp_apip:nnn") == 0) { return &nnn; }
  return nullptr;
}

/* ---------------------------------------------------------------------- */
int FixLambdaLACSPAPIP::pack_reverse_comm(int n, int first, double *buf)
{
  int i,m,last;
  double *avg = atom->apip_la_avg;
  double *norm = atom->apip_la_norm;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    buf[m++] = avg[i];
    buf[m++] = norm[i];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void FixLambdaLACSPAPIP::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i,j,m;
  double *avg = atom->apip_la_avg;
  double *norm = atom->apip_la_norm;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    avg[j] += buf[m++];
    norm[j] += buf[m++];
  }
}

/**
  * Send csp to neighbours.
  */

int FixLambdaLACSPAPIP::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i, j, m;
  double *la_inp = atom->apip_la_inp;
  double *lambda = atom->apip_lambda;
  m = 0;

  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = la_inp[j];
    buf[m++] = lambda[j];
  }

  return m;
}

/**
  * Recv csp from neighbours.
  */

void FixLambdaLACSPAPIP::unpack_forward_comm(int n, int first, double *buf)
{
  int i, m, last;
  double *la_inp = atom->apip_la_inp;
  double *lambda = atom->apip_lambda;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    la_inp[i] = buf[m++];
    lambda[i] = buf[m++];
  }
}


/**
  * Initial f_lambda.
  * This function is called from the force calculation routine.
  */

void FixLambdaLACSPAPIP::store_f_lambda_before()
{
  int nlocal = atom->nlocal;
  double **f = atom->f;

  if (nlocal > size_f_lambda) {
    memory->destroy(f_lambda);
    size_f_lambda = nlocal;
    memory->create(f_lambda, size_f_lambda, size_peratom_cols, "pair:lambda:la:csp:apip:f:lambda");
    array_atom = f_lambda;
  }

  for (int i = 0; i < nlocal; i++) {
    f_lambda[i][0] = - f[i][0];
    f_lambda[i][1] = - f[i][1];
    f_lambda[i][2] = - f[i][2];
  }
}

/**
  * Calculate f_lambda.
  * This function is called from the force calculation routine.
  */

void FixLambdaLACSPAPIP::store_f_lambda_after()
{
  int nlocal = atom->nlocal;
  double **f = atom->f;

  for (int i = 0; i < nlocal; i++) {
    f_lambda[i][0] += f[i][0];
    f_lambda[i][1] += f[i][1];
    f_lambda[i][2] += f[i][2];
  }
}



CSPpairAPIP::CSPpairAPIP(int i0, int i1) {
  if (i0 < i1 ) {
    tag_smaller = i0;
    tag_larger = i1;
  } else {
    tag_smaller = i1;
    tag_larger = i0;
  }
}

double FixLambdaLACSPAPIP::compute_scalar()
{
  int counter;
  MPI_Allreduce(&counter_changed_csp_nghs, &counter, 1, MPI_INT, MPI_SUM, world);
  return counter;
}
