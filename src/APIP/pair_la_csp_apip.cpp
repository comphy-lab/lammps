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

#include "pair_la_csp_apip.h"

#include "atom.h"
#include "atom_vec_apip.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "force.h"
#include "math_const.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "neighbor.h"
#include "pair_hybrid.h"
#include "update.h"

#include <exception>

#include "ace-evaluator/ace_c_basis.h"
#include "ace-evaluator/ace_evaluator.h"
#include "ace-evaluator/ace_recursive.h"
#include "ace-evaluator/ace_version.h"
#include "ace/ace_b_basis.h"

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */
PairLACSPAPIP::PairLACSPAPIP(LAMMPS *lmp) : Pair(lmp), fix_la_csp(nullptr),
  prefactor1(nullptr), prefactor2(nullptr)
{
  prefactor1_size = 0;
  prefactor2_size = 0;
  cutsq_combined = 0;
  fix_la_csp_w_cutsq = 0;
  nnn_half = 0;

  single_enable = 0;
  restartinfo = 0;
  one_coeff = 1;
  manybody_flag = 1;

  scale = nullptr;

  comm_forward = 1;

  // start of adaptive-precision modifications by DI
  n_computations_accumulated = 0;
  time_wall_accumulated = 0;
  time_per_atom = -1;
  // end of adaptive-precision modifications by DI
}

/* ----------------------------------------------------------------------
   check if allocated, since class can be destructed when incomplete
------------------------------------------------------------------------- */

PairLACSPAPIP::~PairLACSPAPIP()
{
  if (copymode) return;

  memory->destroy(prefactor1);
  memory->destroy(prefactor2);

  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(scale);
  }
}

/* ---------------------------------------------------------------------- */

void PairLACSPAPIP::compute(int eflag, int vflag)
{
  // start timers
  double time_wall_start = platform::walltime();
  int n_computations = 0;

  int i, j, ii, jj, inum, jnum, k, i_pair, pair_ngh_0, pair_ngh_1;
  int *ilist, *jlist, *numneigh, **firstneigh, **ngh_pairs;
  double **x, **f, *lambda, *csp, *csp_avg, *csp_norm, *e_fast, *e_precise;
  double xtmp, ytmp, ztmp, lambdatmp, fpair, delx, dely, delz, r, rsq, cspavgtmp, prefactortmp, delx0, dely0, delz0, delx1, dely1, delz1;

  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;

  ev_init(eflag, vflag);

  x = atom->x;
  f = atom->f;
  lambda = atom->apip_lambda;
  csp = atom->apip_la_inp;
  csp_norm = atom->apip_la_norm;
  csp_avg = atom->apip_la_avg;
  e_fast = atom->apip_e_fast;
  e_precise = atom->apip_e_precise;

  ngh_pairs = (int **) fix_la_csp->extract("fix_lambda_la_csp_apip:ngh_pairs", i);

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;


  if (atom->nmax > prefactor1_size) {
    memory->destroy(prefactor1);
    prefactor1_size = atom->nmax;
    memory->create(prefactor1,prefactor1_size,"pair:la:csp:apip:prefactor1");
  }
  if (nlocal > prefactor2_size) {
    memory->destroy(prefactor2);
    prefactor2_size = nlocal;
    memory->create(prefactor2,prefactor2_size,"pair:la:csp:apip:prefactor2");
  }


  // e_fast and e_precise are known only for own atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    prefactor1[i] = fix_la_csp->der_switching_function_poly(csp_avg[i]);
    // e_fast and e_precise are computed only for lambda in (0,1)
    // lambda in (0,1) implies that the derivative of the switching function is non-zero
    if (prefactor1[i] != 0) prefactor1[i] *= (e_fast[i] - e_precise[i]) / csp_norm[i];
    prefactor2[i] = prefactor1[i] * fix_la_csp->weighting_function_poly(0);
  }

  // communicate prefactor1 to neighbouring processors
  // and get prefactor1 of ghosts
  comm->forward_comm(this);


  // store all forces before force computation
  fix_la_csp->store_f_lambda_before();


  // compute derivative of the radial weight function first
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    lambdatmp = lambda[i];

    prefactortmp = prefactor1[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    cspavgtmp = csp_avg[i];

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      j &= NEIGHMASK;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      if (rsq >= fix_la_csp_w_cutsq) continue;

      r = sqrt(rsq);
      prefactor2[i] += prefactor1[j] * fix_la_csp->weighting_function_poly(r);
      fpair = prefactortmp * (csp[j] - cspavgtmp) * fix_la_csp->der_weighting_function_poly(r) / r;

      if (fpair == 0) continue;

      f[i][0] -= delx * fpair;
      f[i][1] -= dely * fpair;
      f[i][2] -= delz * fpair;
      f[j][0] += delx * fpair;
      f[j][1] += dely * fpair;
      f[j][2] += delz * fpair;

      if (vflag) ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
    }
  }


  // compute derivative of the CSP
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    fpair = -2 * prefactor2[i];

    if (fpair == 0) continue;

    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];


    for (i_pair = 0; i_pair < nnn_half; i_pair++) {
      pair_ngh_0 = ngh_pairs[i][i_pair];
      pair_ngh_1 = ngh_pairs[i][i_pair + nnn_half];

      delx0 = xtmp - x[pair_ngh_0][0];
      delx1 = xtmp - x[pair_ngh_1][0];
      dely0 = ytmp - x[pair_ngh_0][1];
      dely1 = ytmp - x[pair_ngh_1][1];
      delz0 = ztmp - x[pair_ngh_0][2];
      delz1 = ztmp - x[pair_ngh_1][2];

      delx = delx0 + delx1;
      dely = dely0 + dely1;
      delz = delz0 + delz1;

      f[i][0] += 2 * delx * fpair;
      f[i][1] += 2 * dely * fpair;
      f[i][2] += 2 * delz * fpair;

      f[pair_ngh_0][0] -= delx * fpair;
      f[pair_ngh_0][1] -= dely * fpair;
      f[pair_ngh_0][2] -= delz * fpair;

      f[pair_ngh_1][0] -= delx * fpair;
      f[pair_ngh_1][1] -= dely * fpair;
      f[pair_ngh_1][2] -= delz * fpair;

      if (vflag_either) {
        ev_tally_xyz(i, pair_ngh_0, nlocal, newton_pair, 0.0, 0.0, delx*fpair, dely*fpair, delz*fpair, delx0, dely0, delz0);
        ev_tally_xyz(i, pair_ngh_1, nlocal, newton_pair, 0.0, 0.0, delx*fpair, dely*fpair, delz*fpair, delx1, dely1, delz1);
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();

  // calculate force due to gradient lambda by comparison with before stored force
  fix_la_csp->store_f_lambda_after();

  // stop timers
  time_wall_accumulated += platform::walltime() - time_wall_start;
  n_computations_accumulated += n_computations;
}

/* ---------------------------------------------------------------------- */

void PairLACSPAPIP::allocate()
{
  allocated = 1;
  int n = atom->ntypes + 1;

  memory->create(setflag, n, n, "pair:setflag");
  for (int i = 1; i < n; i++)
    for (int j = i; j < n; j++) setflag[i][j] = 0;
  memory->create(cutsq, n, n, "pair:cutsq");
  memory->create(scale, n, n, "pair:scale");
  map = new int[n];
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairLACSPAPIP::settings(int narg, char **arg)
{
  if (narg != 0) utils::missing_cmd_args(FLERR, "pair_style la/csp/apip", error);
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairLACSPAPIP::coeff(int narg, char **arg)
{

  if (!allocated) allocate();

  if (cutsq_combined == 0) {
    // find fix lambda/la/csp/apip
    int count = 0;
    for (int i = 0; i < modify->nfix; i++) {
      if (strcmp(modify->fix[i]->style, "lambda/la/csp/apip") == 0) {
        fix_la_csp = (FixLambdaLACSPAPIP *) modify->fix[i];
        count++;
      }
    }
    if (count != 1) error->all(FLERR, "Exact one fix lambda/la/csp/apip required");
    int dim = 0;
    cutsq_combined = *((double *) fix_la_csp->extract("fix_la_csp_apip:cutsq_combined", dim));
  }

  const int n = atom->ntypes;
  // initialize scale factor
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      scale[i][j] = 1.0;
      setflag[i][j] = 1;
      cutsq[i][j] = cutsq_combined;
    }
  }
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairLACSPAPIP::init_style()
{
  if (strcmp(atom->atom_style, "apip/la")) error->all(FLERR, "atom style apip/la required");

  if (fix_la_csp == nullptr) {
    // find fix lambda/la/csp/apip
    int count = 0;
    for (int i = 0; i < modify->nfix; i++) {
      if (strcmp(modify->fix[i]->style, "lambda/la/csp/apip") == 0) {
        fix_la_csp = (FixLambdaLACSPAPIP *) modify->fix[i];
        count++;
      }
    }
    if (count != 1) error->all(FLERR, "Exact one fix lambda/la/csp/apip required");
    int dim = 0;
    cutsq_combined = *((double *) fix_la_csp->extract("fix_la_csp_apip:cutsq_combined", dim));
  }

  if (strcmp(force->pair_style, "hybrid/overlay")) error->all(FLERR, "pair la/csp/apip needs to be used with pair hybrid/overlay");
  auto hybrid = dynamic_cast<PairHybrid *>(force->pair);
  if (hybrid->nstyles < 3) error->all(FLERR, "pair la/csp/apip relies on the computation of e_fast and e_precise from other hybrid substyles");
  if (strcmp(hybrid->keywords[hybrid->nstyles -1], "la/csp/apip")) error->all(FLERR, "pair la/csp/apip needs to be the last hybrid substyle so that e_precise and e_fast are computed before la/csp/apip is used");

  int dim2 = 0;
  fix_la_csp_w_cutsq = *((double *) fix_la_csp->extract("fix_la_csp_apip:cut_hi_sq", dim2));
  nnn_half = (*((int *) fix_la_csp->extract("fix_la_csp_apip:nnn", dim2))) / 2;

  // request a full neighbor list
  auto req = neighbor->add_request(this, NeighConst::REQ_FULL);
  req->set_cutoff(sqrt(cutsq_combined));
}

/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairLACSPAPIP::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    setflag[i][j] = 1;
    cutsq[i][j] = cutsq_combined;
    scale[i][j] = 1.0;
  }
  scale[j][i] = scale[i][j];
  return sqrt(cutsq_combined);
}

/* ---------------------------------------------------------------------- */

int PairLACSPAPIP::pack_forward_comm(int n, int *list, double *buf, int /*pbc_flag*/, int * /*pbc*/)
{
  int i,j,m;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    buf[m++] = prefactor1[j];
  }
  return m;
}

/* ---------------------------------------------------------------------- */

void PairLACSPAPIP::unpack_forward_comm(int n, int first, double *buf)
{
  int i,m,last;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) prefactor1[i] = buf[m++];
}


/**
  * set return values for timers and number of computed particles
  */

// written by DI. This function is required for the adaptive-precision.
void PairLACSPAPIP::calculate_time_per_atom()
{
  if (n_computations_accumulated > 0)
    time_per_atom = time_wall_accumulated / n_computations_accumulated;
  else
    time_per_atom = -1;

  // reset
  time_wall_accumulated = 0;
  n_computations_accumulated = 0;
}

/* ----------------------------------------------------------------------
    extract method for extracting value of scale variable
 ---------------------------------------------------------------------- */
void *PairLACSPAPIP::extract(const char *str, int &dim)
{
  dim = 0;
  if (strcmp(str, "la/csp/apip:time_per_atom") == 0) {
    calculate_time_per_atom();
    return (void *) &time_per_atom;
  }

  dim = 2;
  if (strcmp(str, "scale") == 0) return (void *) scale;
  return nullptr;
}
