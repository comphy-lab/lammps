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
/* ------------------------------------------------------
    This file is part of the BOCS package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    and Maria C. Lesniewski, mjl6766@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */

#include "pair_ldd.h"

#include "atom.h"
#include "citeme.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "integrate.h"
#include "memory.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"
#include "utils.h"

#include "ldd_indicator.h"
#include "ldd_potential.h"
#include "the_ldd_indicator_types.h"
#include "the_ldd_potential_types.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

#define DOT_PRODUCT(a, b) (a[0] * b[0] + a[1] * b[1] + a[2] * b[2])
#define DOT_PROD_GRAD(a, ta, b, tb) \
  (a[3 * ta] * b[3 * tb] + a[3 * ta + 1] * b[3 * tb + 1] + a[3 * ta + 2] * b[3 * tb + 2])
#define GRADTYPE(a) (3 * a)

namespace {
constexpr int MAXLINE = 1024;

constexpr char KEY_LDD_IND[] = "indicator";
constexpr char KEY_LDD_POTL[] = "potential";
constexpr char KEY_LDD_GRAD[] = "gradient";
constexpr char KEY_LDD_SELF[] = "self";
constexpr char KEY_LDD_IGNORE[] = "ignore";

constexpr char cite_pair_ldd1_c[] = "pair ldd command: doi:10.1063/1.5128665\n\n"
                                    "@Article{DeLyser1,\n"
                                    " author = {Michael R. DeLyser and W. G. Noid},\n"
                                    " title = {Analysis of local density potentials},\n"
                                    " journal = {The journal of chemical physics},\n"
                                    " year =    2019,\n"
                                    " volume =  151,\n"
                                    " pages =   {22:224106}\n"
                                    "}\n\n";

constexpr char cite_pair_ldd2_c[] =
    "pair ldd command gradient keyword: doi:10.1063/5.0075291\n\n"
    "@Article{DeLyser2,\n"
    " author = {Michael R. DeLyser and W. G. Noid},\n"
    " title = {Coarse-grained models for local density gradients},\n"
    " journal = {The Journal of Chemical Physics},\n"
    " year =    2021,\n"
    " volume =  156,\n"
    " pages =   {034106}\n"
    "}\n\n";
}    // namespace

/* ---------------------------------------------------------------------- */

PairLdd::PairLdd(LAMMPS *lmp) : Pair(lmp)
{
  if (lmp->citeme) lmp->citeme->add(cite_pair_ldd1_c);
  if (lmp->citeme) lmp->citeme->add(cite_pair_ldd2_c);

  writedata = 1;

  restartinfo = 0;

  // MCL 07.30.25, while we have a single routine to compute the pair contribution
  // given info about the LD, we have no way to set up the LD info without compute at present.
  // So I disabled this.
  single_enable = 0;
  one_coeff = 1;
  // We pass the local densities & 3 components of the gradients for each type
  comm_forward = 4 * atom->ntypes;
  comm_reverse = 4 * atom->ntypes;

  // Initialize these to NULL
  Inds = nullptr;
  Potls = nullptr;
  GradPotls = nullptr;

  // Somewhere we need to check to make sure the user is using the correct
  // atom_style. I did it here because it's the first thing called when
  // the user tries to use this pair type.
  // This can be moved if there's a better place for it elsewhere.
  if (!atom->ldd_big_flag) error->all(FLERR, "atomstyle ldd must be used with pair style ldd");

  if (atom->ldd_ntypes != atom->ntypes) error->all(FLERR, "ldd_ntypes doesn't match ntypes");

  // This is the same as is done in force.cpp
  char *str = (char *) "none";
  int n = strlen(str) + 1;
  indicator_style = new char[n];
  strcpy(indicator_style, str);

  potential_style = new char[n];
  strcpy(potential_style, str);
  map = new int[atom->ntypes + 1];

  LDD_factory();
}

// analogous to void _noopt Force::create_factories()

void PairLdd::LDD_factory()
{

  indicator_map = new IndicatorCreatorMap();

#define LDD_INDICATOR_CLASS
#define LddIndicatorStyle(key, Class) (*indicator_map)[#key] = &indicator_creator<Class>;

#include "the_ldd_indicator_types.h"

#undef LddIndicatorStyle
#undef LDD_INDICATOR_CLASS

  potential_map = new PotentialCreatorMap();

#define LDD_POTENTIAL_CLASS
#define LddPotentialStyle(key, Class) (*potential_map)[#key] = &potential_creator<Class>;

#include "the_ldd_potential_types.h"

#undef LddPotentialStyle
#undef LDD_POTENTIAL_CLASS
}

/* ---------------------------------------------------------------------- */

PairLdd::~PairLdd()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(ignore_pair);
    memory->destroy(ignore_me);
    memory->destroy(bGradient);
    memory->destroy(self_interaction);
    memory->destroy(Inds);
    memory->destroy(Potls);
    memory->destroy(GradPotls);
    allocated = 0;
  }
}

/* ---------------------------------------------------------------------- */

void PairLdd::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  memory->create(setflag, n + 1, n + 1, "pair:setflag");
  for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) setflag[i][j] = 0;
  }

  memory->create(cutsq, n + 1, n + 1, "LDD:cutsq");
  memory->create(cut, n + 1, n + 1, "LDD:cut");
  memory->create(ignore_me, n + 1, "LDD:ignore_me");
  memory->create(ignore_pair, n + 1, n + 1, "LDD:ignore_pair");
  memory->create(bGradient, n + 1, n + 1, "LDD:bGradient");
  memory->create(self_interaction, n + 1, n + 1, "LDD:self_interaction");
  memory->create(Inds, n + 1, n + 1, 1, "LDD:indicators");
  memory->create(Potls, n + 1, n + 1, 1, "LDD:potentials");
  memory->create(GradPotls, n + 1, n + 1, 1, "LDD:gradpotls");

  /* MCL 10.04.24
  * This seems like the wrong way to use this function
  * [Interaction A B and B A point to the same address when we do this]
  * Instead I called it over in atom_vec_ldd when the atom_style is called
  * This seems to fix the address issues
  if (n > 1)
  {
    atom->add_peratom_change_columns("ldd_local_density",n+1);
    atom->add_peratom_change_columns("ldd_energy",n+1);
    atom->add_peratom_change_columns("ldd_grad_density",3*(n+1));
    atom->add_peratom_change_columns("ldd_grad_energy",n+1);
  }
  */

  // MCL ignore solution fix - we want to ignore a type if it has no specified ldds
  // 10.04.24                - so I set the default to ignore and change it only if spec/
  for (int i = 0; i <= n; i++) ignore_me[i] = true;
}

/* ---------------------------------------------------------------------- */

void PairLdd::compute(int eflag, int vflag)
{
  // Standard stuff from other force file
  int i, j, ii, jj, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz, evdwl, fpair;
  double rsq;
  double fpair_LD;

  evdwl = 0.0;
  ev_init(eflag,
          vflag);    // MCL 09.24.25, this lets per atom energies be set up so we can talk to them
  // OLD MRD way
  //if (eflag || vflag) ev_setup(eflag,vflag);
  //else evflag = vflag_fdotr = 0;

  const double *const *const x = atom->x;
  double *const *const f = atom->f;
  const int *const type = atom->type;
  const int nlocal = atom->nlocal;
  double fxtmp, fytmp, fztmp;
  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;
  int newton_pair = force->newton_pair;

  // ldd stuff
  double r_pair, norm_fact;
  double **const local_dens = atom->ldd_local_density;
  double *const LD_ttl_nrg = atom->ldd_total_energy;

  // gradient stuff
  // gf1 = grad force 1. used for the part of the force in the ij dir
  // gf2 = grad force 2. used for the part of the force in the ld gradient dir
  double gf1, gf2[3], eij[3];
  double **const grad_dens = atom->ldd_grad_density;

  LDD_calculate_LDs();
  LDD_calculate_energies();

  // loop over this processor's atoms
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    const int itype = type[i];
    // if this type of atom doesn't have any ldd potentials, skip this
    if (!ignore_me[itype]) {
      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      fxtmp = fytmp = fztmp = 0.0;

      const int *const jlist = firstneigh[i];
      const int jnum = numneigh[i];
      // loop over this atom's neighbors
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        jtype = type[j];
        // make sure there's a ld potential for this pair type
        if ((!ignore_pair[itype][jtype]) || (!ignore_pair[jtype][itype])) {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;

          if (rsq < cutsq[itype][jtype]) {
            r_pair = sqrt(rsq);
            // The pair force decomposition for ldd potentials, divided by r_pair
            fpair_LD = 0.0;
            if (!ignore_pair[itype][jtype]) {
              fpair_LD +=
                  Inds[itype][jtype]->wp(r_pair) * Potls[itype][jtype]->f(local_dens[i][jtype]);
            }
            if (!ignore_pair[jtype][itype]) {
              fpair_LD +=
                  Inds[jtype][itype]->wp(r_pair) * Potls[jtype][itype]->f(local_dens[j][itype]);
            }
            fpair_LD /= r_pair;

            f[i][0] += delx * fpair_LD;
            f[i][1] += dely * fpair_LD;
            f[i][2] += delz * fpair_LD;

            if (newton_pair || j < nlocal) {
              f[j][0] -= delx * fpair_LD;
              f[j][1] -= dely * fpair_LD;
              f[j][2] -= delz * fpair_LD;
            }
            // set fpair as the ld pair force.
            // needed for virial/pressure calculation later.
            // will be overwritten if there's a gradient force.
            fpair = fpair_LD;
            // check if there's a gradient force for this pair of atom types
            if (bGradient[itype][jtype] || bGradient[jtype][itype]) {
              // calculate unit vector from atom j to atom i
              eij[0] = delx / r_pair;
              eij[1] = dely / r_pair;
              eij[2] = delz / r_pair;
              gf1 = gf2[0] = gf2[1] = gf2[2] = 0.0;

              if (bGradient[itype][jtype]) {
                gf1 += GradPotls[itype][jtype]->f(local_dens[i][jtype]) *
                        DOT_PROD_GRAD(grad_dens[i], jtype, grad_dens[i], jtype) *
                        //DOT_PRODUCT(grad_dens[i][jtype],grad_dens[i][jtype]) *
                        Inds[itype][jtype]->wp(r_pair) -
                    2.0 * GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
                        (Inds[itype][jtype]->wp2(r_pair) -
                         Inds[itype][jtype]->wp(r_pair) / r_pair) *
                        //DOT_PRODUCT(eij, grad_dens[i][jtype]);
                        DOT_PROD_GRAD(eij, 0, grad_dens[i], jtype);

                gf2[0] -= 2.0 * GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
                    Inds[itype][jtype]->wp(r_pair) / r_pair * grad_dens[i][GRADTYPE(jtype)];
                //                       grad_dens[i][jtype][0];i
                gf2[1] -= 2.0 * GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
                    Inds[itype][jtype]->wp(r_pair) / r_pair * grad_dens[i][GRADTYPE(jtype) + 1];
                //                       grad_dens[i][jtype][1];

                gf2[2] -= 2.0 * GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
                    Inds[itype][jtype]->wp(r_pair) / r_pair * grad_dens[i][GRADTYPE(jtype) + 2];
                //                       grad_dens[i][jtype][2];
              }
              if (bGradient[jtype][itype]) {
                gf1 += GradPotls[jtype][itype]->f(local_dens[j][itype]) *
                        DOT_PROD_GRAD(grad_dens[j], itype, grad_dens[j], itype) *
                        //DOT_PRODUCT(grad_dens[j][itype],grad_dens[j][itype]) *
                        Inds[jtype][itype]->wp(r_pair) +
                    2.0 * GradPotls[jtype][itype]->u(local_dens[j][itype]) *
                        (Inds[jtype][itype]->wp2(r_pair) -
                         Inds[jtype][itype]->wp(r_pair) / r_pair) *
                        DOT_PROD_GRAD(eij, 0, grad_dens[j], itype);

                //DOT_PRODUCT(eij, grad_dens[j][itype]);
                gf2[0] += 2.0 * GradPotls[jtype][itype]->u(local_dens[j][itype]) *
                    Inds[jtype][itype]->wp(r_pair) / r_pair * grad_dens[j][GRADTYPE(itype)];
                //grad_dens[j][itype][0];
                gf2[1] += 2.0 * GradPotls[jtype][itype]->u(local_dens[j][itype]) *
                    Inds[jtype][itype]->wp(r_pair) / r_pair * grad_dens[j][GRADTYPE(itype) + 1];
                //grad_dens[j][itype][1];
                gf2[2] += 2.0 * GradPotls[jtype][itype]->u(local_dens[j][itype]) *
                    Inds[jtype][itype]->wp(r_pair) / r_pair * grad_dens[j][GRADTYPE(itype) + 2];
                //grad_dens[j][itype][2];
              }
              // total gradient force broken down into components
              fxtmp = gf1 * eij[0] + gf2[0];
              fytmp = gf1 * eij[1] + gf2[1];
              fztmp = gf1 * eij[2] + gf2[2];

              f[i][0] += fxtmp;
              f[i][1] += fytmp;
              f[i][2] += fztmp;
              if (newton_pair || j < nlocal) {
                f[j][0] -= fxtmp;
                f[j][1] -= fytmp;
                f[j][2] -= fztmp;
              }
              // add original ld force components
              fxtmp += fpair_LD * delx;
              fytmp += fpair_LD * dely;
              fztmp += fpair_LD * delz;
              // The gradient force contains components not directed along r_ij.
              // Accordingly, we have to call ev_tally_xyz to properly update the virial tensor.
              if (evflag)
                ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fxtmp, fytmp, fztmp, delx, dely,
                             delz);
            }
            // If there's no gradient term, call ev_tally like normal.
            else {
              if (evflag) {
                ev_tally(i, j, nlocal, newton_pair, 0.0, 0.0, fpair, delx, dely, delz);
              }
            }
          }
        }
      }
      // Increment van der Waals' energy with ld potl energy for this atom
      evdwl = LD_ttl_nrg[i];
      if (eflag_atom)    //(eatom != NULL) // We need to talk to the per atom energy like we would've in ev_tally so that compute peratom pe works
      {
        eatom[i] += LD_ttl_nrg[i];
      }
      if (evflag) eng_vdwl += evdwl;
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ---------------------------------------------------------------------- */

void PairLdd::settings(int narg, char **arg)
{
  if (narg != 1) error->all(FLERR, "Need one argument for pair_style pair_ldd");

  cut_LDD_global = utils::numeric(FLERR, arg[0], false, lmp);
  neighbor->cutneighmax = cut_LDD_global;
  if (allocated) {
    int i, j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i + 1; j <= atom->ntypes; j++)
        if (setflag[i][j]) cut[i][j] = cut_LDD_global;    // make sure cut_global is declared
  }
}

/* ---------------------------------------------------------------------- */

// Again, these are analogous to what is done in force.h
LddIndicator *PairLdd::new_indicator(const std::string &wtype)
{
  if (indicator_map->find(wtype) != indicator_map->end()) {
    IndicatorCreator indicator_creator = (*indicator_map)[wtype];
    return indicator_creator(lmp);
  }
  error->all(FLERR, utils::check_packages_for_style("LddIndicator", wtype, lmp));
  return nullptr;
}

template <typename T> LddIndicator *PairLdd::indicator_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

LddPotential *PairLdd::new_potential(const std::string &ptype)
{
  if (potential_map->find(ptype) != potential_map->end()) {
    PotentialCreator potential_creator = (*potential_map)[ptype];
    return potential_creator(lmp);
  }
  error->all(FLERR, utils::check_packages_for_style("LddPotential", ptype, lmp));
  return nullptr;
}

template <typename T> LddPotential *PairLdd::potential_creator(LAMMPS *lmp)
{
  return new T(lmp);
}

/* ---------------------------------------------------------------------- */

void PairLdd::ErrorDoubleKeyword(const char *keyword)
{
  //  char *errmsg = (char *) calloc(100,sizeof(char));
  //  sprintf(errmsg,"Found keyword %s twice!\n"
  //                 "Error: double keyword\n",keyword);
  std::string errmsg = fmt::format("Found keyword {} twice!\nError: double keyword\n", keyword);
  error->all(FLERR, errmsg);
}

void PairLdd::ErrorNumKeywordArgs(const char *keyword, const char *arglist)
{
  std::string errmsg = fmt::format(
      "Expected {} to be followed by: {}\n Error: Invalid or missing arguments to keyword {}\n",
      keyword, arglist, keyword);
  error->all(FLERR, errmsg);
}

/* ---------------------------------------------------------------------- */

void PairLdd::coeff_ldd(int narg, char **arg)
{
  /* Line should look like:
 * pair_coeff      i   j (ldd)
 *      indicator wtype    r0    rc
 *      self self_term
 *      potential  potl_type   *potl_coeffs*
 *      gradient grad_type *grad_coeffs*
 *
 * ldd only included if using hybrid/overlay, which is probable
 * indicator, self, and potential are all required.
 * gradient is optional.
 * potl_coeffs and grad_coeffs depend on the type of potential used
 */
  int me;    // Use lammps comm->me instead?
  MPI_Comm_rank(world, &me);
  if (narg < 2) error->all(FLERR, "You must list coefficients for the ldd pair interaction");
  //if (!allocated) allocate();
  int n = atom->ntypes;
  int i, j, ilo, ihi, jlo, jhi;
  int iarg = 0;
  int ignore = 0;

  utils::bounds(FLERR, arg[0], 1, atom->ntypes, ilo, ihi, error);
  utils::bounds(FLERR, arg[1], 1, atom->ntypes, jlo, jhi, error);
  iarg += 2;
  int the_ind = -1, the_potl = -1, the_grad = -1, n_potl_coeffs = 0;
  double *potl_coeffs;

  bool bSelf = false, bIgnore = false;

  // turn these true after we find the keyword.
  // used to check for double occurrences.
  bool bkInd = false, bkPotl = false, bkSelf = false, bkGrad = false;

  // process keywords
  // all keywords are #defined at the top of this file
  while (iarg < narg) {
    if (strcmp(arg[iarg], KEY_LDD_IND) == 0) {
      /* check for double keyword */
      if (bkInd) ErrorDoubleKeyword(KEY_LDD_IND);
      bkInd = true;
      /* make sure the proper number of arguments for this keyword are present */
      if (iarg + 3 >= narg) ErrorNumKeywordArgs(KEY_LDD_IND, "wtype r0 rc");

      for (i = ilo; i <= ihi; ++i) {
        for (j = jlo; j <= jhi; ++j) {
          // make a new indicator function of the proper type
          Inds[i][j] = new_indicator(arg[iarg + 1]);
          // set coefficients for indicator function
          Inds[i][j]->init_coeffs(utils::numeric(FLERR, arg[iarg + 2], false, lmp),
                                  utils::numeric(FLERR, arg[iarg + 3], false, lmp),
                                  domain->dimension);
        }
      }
      // increment word counter
      iarg += 4;
    } else if (strcmp(arg[iarg], KEY_LDD_SELF) == 0) {
      /* check for double keyword */
      if (bkSelf) ErrorDoubleKeyword(KEY_LDD_SELF);
      bkSelf = true;
      /* make sure the proper number of arguments for this keyword are present */
      if (iarg + 1 >= narg) ErrorNumKeywordArgs(KEY_LDD_SELF, "yes/no");

      if (strcmp(arg[iarg + 1], "yes") == 0)
        bSelf = true;
      else if (strcmp(arg[iarg + 1], "no") == 0)
        bSelf = false;
      else {
        std::string errmsg =
            fmt::format("Expected to find either \"yes\"/\"no\" to follow keyword {} instead, we "
                        "found {}\n Error: Invalid argument to keyword\n",
                        KEY_LDD_SELF, arg[iarg + 1]);
        error->all(FLERR, errmsg);
      }
      iarg += 2;
    } else if (strcmp(arg[iarg], KEY_LDD_POTL) == 0) {
      /* check for double keyword */
      if (bkPotl) ErrorDoubleKeyword(KEY_LDD_POTL);
      bkPotl = true;
      /* make sure the proper number of arguments for this keyword are present */
      if (iarg + 1 >= narg) ErrorNumKeywordArgs(KEY_LDD_POTL, "type *args*");

      for (i = ilo; i <= ihi; ++i) {
        for (j = jlo; j <= jhi; ++j) {
          Potls[i][j] = new_potential(arg[iarg + 1]);
          // pass whole line so we can extract potential type arguments
          Potls[i][j]->setup_potl(iarg, narg, arg);
        }
      }
      iarg += (Potls[ilo][jlo]->n_coeffs + 2);
    } else if (strcmp(arg[iarg], KEY_LDD_GRAD) == 0) {
      /* check for double keyword */
      if (bkGrad) ErrorDoubleKeyword(KEY_LDD_GRAD);
      bkGrad = true;
      /* make sure the proper number of arguments for this keyword are present */
      if (iarg + 1 >= narg) ErrorNumKeywordArgs(KEY_LDD_GRAD, "type *args*");

      for (i = ilo; i <= ihi; ++i) {
        for (j = jlo; j <= jhi; ++j) {
          GradPotls[i][j] = new_potential(arg[iarg + 1]);
          GradPotls[i][j]->setup_potl(iarg, narg, arg);
        }
      }
      iarg += (GradPotls[ilo][jlo]->n_coeffs + 2);
    } else if (strcmp(arg[iarg], KEY_LDD_IGNORE) == 0) {
      bIgnore = true;
      iarg += 1;
    } else {
      std::string errmsg = fmt::format("Recognized keywords for pair_coeff ldd: {} {} {} {} {}\n"
                                       "However, we found unrecognized keyword: {}\n"
                                       "Invalid ldd keyword\n",
                                       KEY_LDD_IND, KEY_LDD_SELF, KEY_LDD_POTL, KEY_LDD_GRAD,
                                       KEY_LDD_IGNORE, arg[iarg]);
      error->all(FLERR, errmsg);
    }
  }
  for (i = ilo; i <= ihi; ++i) {
    for (j = jlo; j <= jhi; ++j) {
      bGradient[i][j] = bkGrad;
      ignore_pair[i][j] = bIgnore;
    }
  }

  if (!bIgnore) {
    // Error checks
    if (!bkInd) {
      std::string errmsg = std::string("We never found required keyword ") +
          std::string(KEY_LDD_IND) + std::string("\n");
      error->all(FLERR, errmsg);
    }
    if (Inds[ilo][jlo]->r0 >= Inds[ilo][jlo]->rc) {
      std::string errmsg = fmt::format("r0 must be less than rC. However, you specified r0 = {},"
                                       " rC = {}\n",
                                       Inds[ilo][jlo]->r0, Inds[ilo][jlo]->rc);
      error->all(FLERR, errmsg);
    }

    if (!bkSelf) {
      std::string errmsg = std::string("We never found required keyword ") +
          std::string(KEY_LDD_SELF) + std::string("\n");
      error->all(FLERR, errmsg);
    }
    if (!bkPotl) {
      std::string errmsg = std::string("We never found required keyword ") +
          std::string(KEY_LDD_POTL) + std::string("\n");
      error->all(FLERR, errmsg);
    }
    for (i = ilo; i <= ihi; ++i) {
      for (j = jlo; j <= jhi; ++j) {
        self_interaction[i][j] = false;
        if (bSelf) {
          if (i == j) {
            self_interaction[i][j] = true;
          } else {
            std::string warnmsg =
                fmt::format("WARNING: you said to include the self interaction "
                            "for itype: {} jtype: {}\nHOWEVER, you can only include the "
                            "self interaction for i == j\nAccordingly, we are "
                            "turning this off\n",
                            i, j);
            error->warning(FLERR, warnmsg);
          }
        }
      }
    }
  }

  double cut_one = cut_LDD_global;

  for (int i = ilo; i <= ihi; ++i) {
    bool ignore = true;
    for (int j = jlo; j <= jhi; ++j) {
      cut[i][j] = Inds[i][j]->rc;
      setflag[i][j] = 1;
      if (!ignore_pair[i][j]) {
        ignore = false;
        ignore_me[j] = (ignore && ignore_me[j]);
      }
    }
    ignore_me[i] = (ignore && ignore_me[i]);
  }

  MPI_Barrier(MPI_COMM_WORLD);
}

/* ---------------------------------------------------------------------- */

double PairLdd::init_one(int i, int j)    // perform initializaion for one i,j type pair
{
  return cut_LDD_global;
}

/* ---------------------------------------------------------------------- */

void PairLdd::init_style()    //initialization specific to this pair style
{
  if (force->newton_pair == 0) error->all(FLERR, "Pair style ldd requires newton pair on");
  neighbor->request(this, instance_me);
}

/* ---------------------------------------------------------------------- */

void PairLdd::write_restart_settings(FILE *fp)    // writes global settings to restart files
{
  fwrite(&cut_LDD_global, sizeof(double), 1, fp);
}

/* ---------------------------------------------------------------------- */

void PairLdd::read_restart_settings(FILE *fp)    // reads global settings from restart files
{
  int me = comm->me;
  if (me == 0) fread(&cut_LDD_global, sizeof(double), 1, fp);
  MPI_Bcast(&cut_LDD_global, 1, MPI_DOUBLE, 0, world);
}

/* ---------------------------------------------------------------------- */

double PairLdd::single(int i, int j, int itype, int jtype, double rsq, double, double factor_lj,
                       double &fforce)
{
  fforce = 0.0;
  if ((ignore_pair[itype][jtype]) && (ignore_pair[jtype][itype])) return 0.0;
  double r_pair = sqrt(rsq);
  const double *const *const x = atom->x;
  double **const local_dens = atom->ldd_local_density;
  double **const grad_dens = atom->ldd_grad_density;
  double energy = 0.0;
  double gf1 = 0.0, gf2[3], eij[3];
  gf2[0] = gf2[1] = gf2[2] = 0.0;
  double delx, dely, delz;
  delx = x[i][0] - x[j][0];
  dely = x[i][1] - x[j][1];
  delz = x[i][2] - x[j][2];
  eij[0] = delx / r_pair;
  eij[1] = dely / r_pair;
  eij[2] = delz / r_pair;

  if (!ignore_pair[itype][jtype]) {
    if ((Inds[itype][jtype]->r0 <= r_pair) && (Inds[itype][jtype]->rc >= r_pair)) {
      fforce += Inds[itype][jtype]->wp(r_pair) * Potls[itype][jtype]->f(local_dens[i][jtype]);
      energy += Potls[itype][jtype]->u(local_dens[i][jtype]);
    }
    if (bGradient[itype][jtype]) {
      gf1 += GradPotls[itype][jtype]->f(local_dens[i][jtype]) *
              DOT_PROD_GRAD(grad_dens[i], jtype, grad_dens[i], jtype) *
              //DOT_PRODUCT(grad_dens[i][jtype],grad_dens[i][jtype]) *
              Inds[itype][jtype]->wp(r_pair) -
          2.0 * GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
              (Inds[itype][jtype]->wp2(r_pair) - Inds[itype][jtype]->wp(r_pair) / r_pair) *
              DOT_PROD_GRAD(eij, 0, grad_dens[i], jtype);
      //DOT_PRODUCT(eij,grad_dens[i][jtype]);
      gf2[0] -= 2.0 * GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
          Inds[itype][jtype]->wp(r_pair) / r_pair * grad_dens[i][GRADTYPE(jtype)];
      //grad_dens[i][jtype][0];
      gf2[1] -= 2.0 * GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
          Inds[itype][jtype]->wp(r_pair) / r_pair * grad_dens[i][GRADTYPE(jtype) + 1];
      //grad_dens[i][jtype][1];
      gf2[2] -= 2.0 * GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
          Inds[itype][jtype]->wp(r_pair) / r_pair * grad_dens[i][GRADTYPE(jtype) + 2];
      //grad_dens[i][jtype][2];
      //
      energy += GradPotls[itype][jtype]->u(local_dens[i][jtype]) *
          DOT_PROD_GRAD(grad_dens[i], jtype, grad_dens[i], jtype);
      //DOT_PRODUCT(grad_dens[i][jtype],grad_dens[i][jtype]);
    }
  }

  if (!ignore_pair[jtype][itype]) {
    if ((Inds[jtype][itype]->r0 <= r_pair) && (Inds[jtype][itype]->rc >= r_pair)) {
      fforce += Inds[jtype][itype]->wp(r_pair) * Potls[jtype][itype]->f(local_dens[j][itype]);
      energy += Potls[jtype][itype]->u(local_dens[j][itype]);
    }
    if (bGradient[jtype][itype]) {
      gf1 += GradPotls[jtype][itype]->f(local_dens[j][itype]) *
              DOT_PROD_GRAD(grad_dens[j], itype, grad_dens[j], itype) *
              //DOT_PRODUCT(grad_dens[j][itype],grad_dens[j][itype]) *
              Inds[jtype][itype]->wp(r_pair) +
          2.0 * GradPotls[jtype][itype]->u(local_dens[j][itype]) *
              (Inds[jtype][itype]->wp2(r_pair) - Inds[jtype][itype]->wp(r_pair) / r_pair) *
              DOT_PROD_GRAD(eij, 0, grad_dens[j], itype);
      //DOT_PRODUCT(eij, grad_dens[j][itype]);
      gf2[0] += 2.0 * GradPotls[jtype][itype]->u(local_dens[j][itype]) *
          Inds[jtype][itype]->wp(r_pair) / r_pair * grad_dens[j][GRADTYPE(itype)];
      //grad_dens[j][itype][0];
      gf2[1] += 2.0 * GradPotls[jtype][itype]->u(local_dens[j][itype]) *
          Inds[jtype][itype]->wp(r_pair) / r_pair * grad_dens[j][GRADTYPE(itype) + 1];
      //grad_dens[j][itype][1];
      gf2[2] += 2.0 * GradPotls[jtype][itype]->u(local_dens[j][itype]) *
          Inds[jtype][itype]->wp(r_pair) / r_pair * grad_dens[j][GRADTYPE(itype) + 2];
      //grad_dens[j][itype][2];
      energy += GradPotls[jtype][itype]->u(local_dens[j][itype]) *
          DOT_PROD_GRAD(grad_dens[j], itype, grad_dens[j], itype);
      //DOT_PRODUCT(grad_dens[j][itype],grad_dens[j][itype]);
    }
  }
  double fx = (fforce + gf1) * delx / r_pair + gf2[0];
  double fy = (fforce + gf1) * dely / r_pair + gf2[1];
  double fz = (fforce + gf1) * delz / r_pair + gf2[2];
  fforce = sqrt(fx * fx + fy * fy + fz * fz);
  return energy;
}

/* ---------------------------------------------------------------------- */

void PairLdd::LDD_calculate_LDs()    //
{
  int i, j, m, ii, jj, jtype;
  double xtmp, ytmp, ztmp, delx, dely, delz;
  double rsq;

  const double *const *const x = atom->x;
  const int *const type = atom->type;
  const int nlocal = atom->nlocal;
  const int ntypes = atom->ntypes;
  double r_pair, wprime;    //MINE
  const int inum = list->inum;
  const int *const ilist = list->ilist;
  const int *const numneigh = list->numneigh;
  const int *const *const firstneigh = list->firstneigh;
  int newton_pair = force->newton_pair;

  double **const local_dens = atom->ldd_local_density;
  double **const grad_dens = atom->ldd_grad_density;

  if (newton_pair) {
    m = nlocal + atom->nghost;
    for (i = 0; i < m; i++) {
      for (int tidx = 0; tidx <= ntypes; tidx++) {
        local_dens[i][tidx] = 0.0;
        grad_dens[i][GRADTYPE(tidx)] = 0.0;
        grad_dens[i][GRADTYPE(tidx) + 1] = 0.0;
        grad_dens[i][GRADTYPE(tidx) + 2] = 0.0;
      }
    }
  } else {
    for (i = 0; i < nlocal; i++) {
      for (int tidx = 0; tidx <= ntypes; tidx++) {
        local_dens[i][tidx] = 0.0;
        grad_dens[i][GRADTYPE(tidx)] = 0.0;
        grad_dens[i][GRADTYPE(tidx) + 1] = 0.0;
        grad_dens[i][GRADTYPE(tidx) + 2] = 0.0;
      }
    }
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    const int itype = type[i];
    if (!ignore_me[itype]) {

      xtmp = x[i][0];
      ytmp = x[i][1];
      ztmp = x[i][2];
      const int *const jlist = firstneigh[i];
      const int jnum = numneigh[i];
      // self interaction
      if (self_interaction[itype][itype]) local_dens[i][itype] += Inds[itype][itype]->invnorm;
      for (jj = 0; jj < jnum; jj++) {
        j = jlist[jj];
        j &= NEIGHMASK;
        jtype = type[j];
        // MCL We need the displacements if relevent to either interaction
        if ((!ignore_pair[itype][jtype]) || (!ignore_pair[jtype][itype])) {
          delx = xtmp - x[j][0];
          dely = ytmp - x[j][1];
          delz = ztmp - x[j][2];
          rsq = delx * delx + dely * dely + delz * delz;

          if (rsq < cutsq[itype][jtype] && (!ignore_pair[itype][jtype])) {
            // Compute local density
            r_pair = sqrt(rsq);
            local_dens[i][jtype] += Inds[itype][jtype]->w(r_pair);
            wprime = Inds[itype][jtype]->wp(r_pair);
            grad_dens[i][GRADTYPE(jtype)] += wprime * delx / r_pair;
            grad_dens[i][GRADTYPE(jtype) + 1] += wprime * dely / r_pair;
            grad_dens[i][GRADTYPE(jtype) + 2] += wprime * delz / r_pair;
          }
          if (newton_pair || j < nlocal) {
            if (!ignore_pair[jtype][itype]) {
              if (rsq < cutsq[jtype][itype]) {
                r_pair = sqrt(rsq);
                local_dens[j][itype] += Inds[jtype][itype]->w(r_pair);
                wprime = Inds[jtype][itype]->wp(r_pair);
                grad_dens[j][GRADTYPE(itype)] -= wprime * delx / r_pair;
                grad_dens[j][GRADTYPE(itype) + 1] -= wprime * dely / r_pair;
                grad_dens[j][GRADTYPE(itype) + 2] -= wprime * delz / r_pair;
              }
            }
          }
        }
      }
    }
  }

  if (newton_pair) comm->reverse_comm(this);    //MCL 4.10.10 on LAMMPS DEV
  comm->forward_comm(this);
}

/* ---------------------------------------------------------------------- */

void PairLdd::LDD_calculate_energies()
{
  int i, ii;

  const int *const type = atom->type;
  const int ntypes = atom->ntypes;
  const int inum = list->inum;
  const int *const ilist = list->ilist;

  double **const local_dens = atom->ldd_local_density;
  double **const grad_dens = atom->ldd_grad_density;
  double **const ld_nrg = atom->ldd_energy;
  double **const ld_grad_nrg = atom->ldd_grad_energy;
  double *const ld_ttl_nrg = atom->ldd_total_energy;
  double Ai2;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    const int itype = type[i];
    ld_ttl_nrg[i] = 0;
    if (!ignore_me[itype]) {
      for (int tidx = 1; tidx <= ntypes; tidx++) {
        ld_nrg[i][tidx] = 0;
        ld_grad_nrg[i][tidx] = 0;
        if (!ignore_pair[itype][tidx]) {
          ld_nrg[i][tidx] = Potls[itype][tidx]->u(local_dens[i][tidx]);
          if (bGradient[itype][tidx]) {
            //Ai2 = DOT_PRODUCT(grad_dens[i][tidx],grad_dens[i][tidx]);
            Ai2 = DOT_PROD_GRAD(grad_dens[i], tidx, grad_dens[i], tidx);
            ld_grad_nrg[i][tidx] = Ai2 * GradPotls[itype][tidx]->u(local_dens[i][tidx]);
          } else {
            ld_grad_nrg[i][tidx] = 0;
          }

          ld_ttl_nrg[i] += ld_nrg[i][tidx] + ld_grad_nrg[i][tidx];
        }
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairLdd::pack_forward_comm(int n, int *list, double *buf, int pbc_flag, int *pbc)
{
  int i, j, k, m;
  double **ld = atom->ldd_local_density;
  double **ldg = atom->ldd_grad_density;
  const int ntypes = atom->ntypes;

  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 1; k <= ntypes; k++) buf[m++] = ld[j][k];
    for (k = 1; k <= ntypes; k++) {
      buf[m++] = ldg[j][GRADTYPE(k)];
      buf[m++] = ldg[j][GRADTYPE(k) + 1];
      buf[m++] = ldg[j][GRADTYPE(k) + 2];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairLdd::unpack_forward_comm(int n, int first, double *buf)
{
  int i, k, m, last;
  double **ld = atom->ldd_local_density;
  double **ldg = atom->ldd_grad_density;
  const int ntypes = atom->ntypes;

  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (k = 1; k <= ntypes; k++) ld[i][k] = buf[m++];
    for (k = 1; k <= ntypes; ++k) {
      ldg[i][GRADTYPE(k)] = buf[m++];
      ldg[i][GRADTYPE(k) + 1] = buf[m++];
      ldg[i][GRADTYPE(k) + 2] = buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

int PairLdd::pack_reverse_comm(int n, int first, double *buf)
{
  int i, k, m, last;
  double **ld = atom->ldd_local_density;
  double **ldg = atom->ldd_grad_density;
  const int ntypes = atom->ntypes;
  m = 0;
  last = first + n;
  for (i = first; i < last; i++) {
    for (k = 1; k <= ntypes; k++) buf[m++] = ld[i][k];
    for (k = 1; k <= ntypes; k++) {
      buf[m++] = ldg[i][GRADTYPE(k)];
      buf[m++] = ldg[i][GRADTYPE(k) + 1];
      buf[m++] = ldg[i][GRADTYPE(k) + 2];
    }
  }

  return m;
}

/* ---------------------------------------------------------------------- */

void PairLdd::unpack_reverse_comm(int n, int *list, double *buf)
{
  int i, j, k, m;
  double **ld = atom->ldd_local_density;
  double **ldg = atom->ldd_grad_density;
  const int ntypes = atom->ntypes;
  m = 0;
  for (i = 0; i < n; i++) {
    j = list[i];
    for (k = 1; k <= ntypes; k++) ld[j][k] += buf[m++];
    for (k = 1; k <= ntypes; k++) {
      ldg[j][GRADTYPE(k)] += buf[m++];
      ldg[j][GRADTYPE(k) + 1] += buf[m++];
      ldg[j][GRADTYPE(k) + 2] += buf[m++];
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairLdd::coeff(int narg, char **arg)
{
  if (!allocated) allocate();
  // then the command probably looked like pair_coeff * * ele1 ele2 etc.
  if (narg > 3) {
    map_element2type(narg - 3, arg + 3);
    if (nelements != atom->ntypes) {
      error->one(FLERR,
                 "Need unique element map for all atom types in system. nelements: {} ntypes: {}",
                 nelements, atom->ntypes);
    }
  }
  // parse pair ij != ji syntax in a different file to get around messing with the main-line of lammps
  read_file(arg[2], nelements);
}

/* ---------------------------------------------------------------------- */

void PairLdd::read_file(char *filename, int nelements)
{
  FILE *lddinp_fp = nullptr;
  if (comm->me == 0) {
    lddinp_fp = utils::open_potential(filename, lmp, 0);
    if (!lddinp_fp) {
      error->one(FLERR, "Cannot open ldd input file {}: {}", filename, utils::getsyserror());
    }
  }

  char line_buf[MAXLINE];
  std::vector<std::string> ldd_arg_string;
  int arg_string_size = 0;
  char **ldd_arg_chars;
  int num_words = 0;
  int bdone = 0;
  // file reader/arg parser loop
  while (!bdone) {
    utils::read_lines_from_file(lddinp_fp, 1, MAXLINE, line_buf, comm->me,
                                world);    // should broadcast to all

    if (comm->me == 0) {
      if (feof(lddinp_fp)) bdone = 1;
    }    // But only 0 will know if done
    MPI_Bcast(&bdone, 1, MPI_INT, 0, world);
    if (bdone) continue;

    num_words = utils::trim_and_count_words(line_buf, " ");
    if (num_words - 1 <= 0) {
      MPI_Barrier(MPI_COMM_WORLD);
      continue;
    } else {
      ldd_arg_string = utils::split_words(line_buf);
    }    // break line into args
    if (!utils::strsame(ldd_arg_string[0].c_str(), "pair_coeff")) {
      error->all(
          FLERR,
          "ERROR:ldd input file {} only accepts lines leading with \"pair_coeff\" commands not: {}\
                                adjust this line: {} in {}",
          filename, ldd_arg_string[0].c_str(), line_buf, filename);
      num_words = 0;
    }
    // then the user is probably passing e.g. pair_coeff A B instead of 1 2
    // coeff_ldd likes 1 2 so we'll change it out before we pass into there.
    if (nelements > 0) {
      for (int k = 0; k < nelements; k++) {
        if (utils::strsame(ldd_arg_string[1], elements[k])) {
          ldd_arg_string[1] = std::to_string(k + 1);
        }
        if (utils::strsame(ldd_arg_string[2], elements[k])) {
          ldd_arg_string[2] = std::to_string(k + 1);
        }
      }
    }

    ldd_arg_chars = new char *[ldd_arg_string.size() - 1];
    for (int i = 1; i < ldd_arg_string.size(); i++) {
      ldd_arg_chars[i - 1] = const_cast<char *>(ldd_arg_string[i].c_str());
      std::strcpy(ldd_arg_chars[i - 1], ldd_arg_string[i].c_str());
    }
    MPI_Barrier(MPI_COMM_WORLD);
    coeff_ldd(ldd_arg_string.size() - 1, ldd_arg_chars);

    delete[] ldd_arg_chars;

  }    // end parser loop
  MPI_Barrier(MPI_COMM_WORLD);
}
