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
   Contributing author: Jacopo Bilotto (EPFL), Jibril B. Coulibaly
------------------------------------------------------------------------- */

#include "pair_gran_hertz_history_ellipsoid.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "error.h"
#include "fix.h"
#include "fix_dummy.h"
#include "fix_neigh_history.h"
#include "force.h"
#include "math_extra.h"    // probably needed for some computations
#include "math_extra_superellipsoids.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "update.h"

#include <cmath>
#include <cstring>

using namespace LAMMPS_NS;

static constexpr int NUMSTEP_INITIAL_GUESS = 5;

/* ---------------------------------------------------------------------- */

PairGranHertzHistoryEllipsoid::PairGranHertzHistoryEllipsoid(LAMMPS *lmp) :
    PairGranHookeHistoryEllipsoid(lmp)
{
}

/* ---------------------------------------------------------------------- */

void PairGranHertzHistoryEllipsoid::compute(int eflag, int vflag)
{
  int i, j, ii, jj, inum, jnum;
  double xtmp, ytmp, ztmp, delx, dely, delz, fx, fy, fz;
  double radi, radj, radsum, rsq, r, rinv, rsqinv, factor_lj;
  double vr1, vr2, vr3, vnnr, vn1, vn2, vn3, vt1, vt2, vt3;
  double wr1, wr2, wr3;
  double vtr1, vtr2, vtr3, vrel;
  double mi, mj, meff, damp, ccel, tor1, tor2, tor3;
  double fn, fs, fs1, fs2, fs3;
  double shrmag, rsht, polyhertz;
  int *ilist, *jlist, *numneigh, **firstneigh;
  int *touch, **firsttouch;
  double *shear, *X0_prev, *history, *allhistory, **firsthistory;

  double shapex, shapey, shapez;    // ellipsoid shape params
  double quat1, quat2, quat3, quat4;
  double block1, block2;

  double X0[4], nij[3], shapei[3], blocki[3], shapej[3], blockj[3], Ri[3][3], Rj[3][3], overlap1,
      overlap2, omegai[3], omegaj[3];
  AtomVecEllipsoid::BlockType flagi, flagj;

  ev_init(eflag, vflag);

  int shearupdate = 1;
  if (update->setupflag) shearupdate = 0;

  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body

  if (fix_rigid && neighbor->ago == 0) {
    int tmp;
    int *body = (int *) fix_rigid->extract("body", tmp);
    auto *mass_body = (double *) fix_rigid->extract("masstotal", tmp);
    if (atom->nmax > nmax) {
      memory->destroy(mass_rigid);
      nmax = atom->nmax;
      memory->create(mass_rigid, nmax, "pair:mass_rigid");
    }
    int nlocal = atom->nlocal;
    for (i = 0; i < nlocal; i++)
      if (body[i] >= 0)
        mass_rigid[i] = mass_body[body[i]];
      else
        mass_rigid[i] = 0.0;
    comm->forward_comm(this);
  }

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  double **angmom = atom->angmom;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;

  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double *special_lj = force->special_lj;
  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::BonusSuper *bonus = avec_ellipsoid->bonus_super;
  int *ellipsoid = atom->ellipsoid;

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  firsttouch = fix_history->firstflag;
  firsthistory = fix_history->firstvalue;

  // loop over neighbors of my atoms

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    touch = firsttouch[i];
    allhistory = firsthistory[i];
    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];
      factor_lj = special_lj[sbmask(j)];
      j &= NEIGHMASK;

      if (factor_lj == 0) continue;

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;
      radj = radius[j];
      radsum = radi + radj;

      X0_prev = &allhistory[3 + size_history * jj];
      int ref_index = (atom->tag[i] < atom->tag[j]) ? i : j;

      // TODO: Below could be a `touch()` function
      bool touching = false;
      if (rsq >= radsum * radsum) {
        touching = false;
      } else {
        MathExtra::copy3(bonus[ellipsoid[i]].shape, shapei);
        MathExtra::copy3(bonus[ellipsoid[j]].shape, shapej);
        MathExtra::copy3(bonus[ellipsoid[i]].block, blocki);
        MathExtra::copy3(bonus[ellipsoid[j]].block, blockj);
        MathExtra::quat_to_mat(bonus[ellipsoid[i]].quat, Ri);
        MathExtra::quat_to_mat(bonus[ellipsoid[j]].quat, Rj);
        bool skip_contact_detection(false);
        if (bounding_box) {
          int cached_axis = (int) (allhistory[7 + size_history * jj]);
          int new_axis = MathExtraSuperellipsoids::check_oriented_bounding_boxes(
              x[i], Ri, shapei, x[j], Rj, shapej, cached_axis);
          if (new_axis != -1) {
            skip_contact_detection = true;
            allhistory[7 + size_history * jj] = (double) new_axis;
          }
        }
        if (skip_contact_detection)
          touching = false;
        else {
          // superellipsoid contact detection between atoms i and j
          flagi = bonus[ellipsoid[i]].type;
          flagj = bonus[ellipsoid[j]].type;
          if (touch[jj] == 1) {
            // Continued contact: use grain true shape and last contact point with respect to grain i
            X0[0] = x[ref_index][0] + X0_prev[0];
            X0[1] = x[ref_index][1] + X0_prev[1];
            X0[2] = x[ref_index][2] + X0_prev[2];
            X0[3] = X0_prev[3];
            int status = MathExtraSuperellipsoids::determine_contact_point(
                x[i], Ri, shapei, blocki, flagi, x[j], Rj, shapej, blockj, flagj, X0, nij,
                contact_formulation);
            if (status == 0)
              touching = true;
            else if (status == 1)
              touching = false;
            else
              error->warning(FLERR,
                             "Ellipsoid contact detection (old contact) failed"
                             "between particle {} and particle {}",
                             atom->tag[i], atom->tag[j]);
          } else {
            // New contact: Build initial guess incrementally by morphing the particles from spheres to actual shape

            // There might be better heuristic for the "volume equivalent spheres" suggested in the paper
            // but this is good enough. We might even be able to use radi and radj which is cheaper
            // MathExtra::scaleadd3(radj / radsum, x[i], radi /radsum, x[j], X0);

            double reqi = std::cbrt(shapei[0] * shapei[1] * shapei[2]);
            double reqj = std::cbrt(shapej[0] * shapej[1] * shapej[2]);
            MathExtra::scaleadd3(reqj / (reqi + reqj), x[i], reqi / (reqi + reqj), x[j], X0);
            X0[3] = reqj / reqi;    // Lagrange multiplier mu^2
            for (int iter_ig = 1; iter_ig <= NUMSTEP_INITIAL_GUESS; iter_ig++) {
              double frac = iter_ig / double(NUMSTEP_INITIAL_GUESS);
              shapei[0] = shapei[1] = shapei[2] = reqi;
              shapej[0] = shapej[1] = shapej[2] = reqj;
              MathExtra::scaleadd3(1.0 - frac, shapei, frac, bonus[ellipsoid[i]].shape, shapei);
              MathExtra::scaleadd3(1.0 - frac, shapej, frac, bonus[ellipsoid[j]].shape, shapej);
              blocki[0] = 2.0 + frac * (bonus[ellipsoid[i]].block[0] - 2.0);
              blocki[1] = 2.0 + frac * (bonus[ellipsoid[i]].block[1] - 2.0);
              blockj[0] = 2.0 + frac * (bonus[ellipsoid[j]].block[0] - 2.0);
              blockj[1] = 2.0 + frac * (bonus[ellipsoid[j]].block[1] - 2.0);

              // force ellipsoid flag for first initial guess iteration.
              // Avoid incorrect values of n1/n2 - 2 in second derivatives.
              int status = MathExtraSuperellipsoids::determine_contact_point(
                  x[i], Ri, shapei, blocki,
                  iter_ig == 1 ? AtomVecEllipsoid::BlockType::ELLIPSOID : flagi, x[j], Rj, shapej,
                  blockj, iter_ig == 1 ? AtomVecEllipsoid::BlockType::ELLIPSOID : flagj, X0, nij,
                  contact_formulation);
              if (status == 0)
                touching = true;
              else if (status == 1)
                touching = false;
              else if (iter_ig == NUMSTEP_INITIAL_GUESS) {
                // keep trying until last iteration to avoid erroring out too early
                error->warning(FLERR,
                               "Ellipsoid contact detection (new contact) failed"
                               "between particle {} and particle {}",
                               atom->tag[i], atom->tag[j]);
              }
            }
          }
        }
      }

      if (!touching) {
        // unset non-touching neighbors

        touch[jj] = 0;
        history = &allhistory[size_history * jj];
        for (int k = 0; k < size_history; k++) history[k] = 0.0;
      } else {
        // Store contact point with respect to grain i for next time step
        // This is crucial for periodic BCs when grains can move by large amount in one time step
        // Keeping the previous contact point relative to global frame would lead to bad initial guess
        X0_prev[0] = X0[0] - x[ref_index][0];
        X0_prev[1] = X0[1] - x[ref_index][1];
        X0_prev[2] = X0[2] - x[ref_index][2];
        X0_prev[3] = X0[3];

        double nji[3] = {-nij[0], -nij[1], -nij[2]};
        // compute overlap depth along normal direction for each grain
        // overlap is positive for both grains
        overlap1 = MathExtraSuperellipsoids::compute_overlap_distance(shapei, blocki, Ri, flagi, X0,
                                                                      nij, x[i]);
        overlap2 = MathExtraSuperellipsoids::compute_overlap_distance(shapej, blockj, Rj, flagj, X0,
                                                                      nji, x[j]);

        double surf_point_i[3], surf_point_j[3], curvature_i, curvature_j;
        MathExtra::scaleadd3(overlap1, nij, X0, surf_point_i);
        MathExtra::scaleadd3(overlap2, nji, X0, surf_point_j);

        if (curvature_model == MathExtraSuperellipsoids::CURV_MEAN) {
          curvature_i = MathExtraSuperellipsoids::mean_curvature_superellipsoid(
              shapei, blocki, flagi, Ri, surf_point_i, x[i]);
          curvature_j = MathExtraSuperellipsoids::mean_curvature_superellipsoid(
              shapej, blockj, flagj, Rj, surf_point_j, x[j]);
        } else {
          curvature_i = MathExtraSuperellipsoids::gaussian_curvature_superellipsoid(
              shapei, blocki, flagi, Ri, surf_point_i, x[i]);
          curvature_j = MathExtraSuperellipsoids::gaussian_curvature_superellipsoid(
              shapej, blockj, flagj, Rj, surf_point_j, x[j]);
        }

        polyhertz = sqrt((overlap1 + overlap2) /
                         (curvature_i + curvature_j));    // hertzian contact radius approximation

        // branch vectors
        double cr1[3], cr2[3];
        MathExtra::sub3(X0, x[i], cr1);
        MathExtra::sub3(X0, x[j], cr2);

        // we need to take the cross product of omega

        double ex_space[3], ey_space[3], ez_space[3];
        MathExtra::q_to_exyz(bonus[ellipsoid[i]].quat, ex_space, ey_space, ez_space);
        MathExtra::angmom_to_omega(angmom[i], ex_space, ey_space, ez_space,
                                   bonus[ellipsoid[i]].inertia, omegai);
        MathExtra::q_to_exyz(bonus[ellipsoid[j]].quat, ex_space, ey_space, ez_space);
        MathExtra::angmom_to_omega(angmom[j], ex_space, ey_space, ez_space,
                                   bonus[ellipsoid[j]].inertia, omegaj);

        double omega_cross_r1[3], omega_cross_r2[3];
        MathExtra::cross3(omegai, cr1, omega_cross_r1);
        MathExtra::cross3(omegaj, cr2, omega_cross_r2);

        // relative translational velocity
        // compute directly the sum of relative translational velocity at contact point
        // since rotational velocity contribution is different for superellipsoids
        double cv1[3], cv2[3];

        cv1[0] = v[i][0] + omega_cross_r1[0];
        cv1[1] = v[i][1] + omega_cross_r1[1];
        cv1[2] = v[i][2] + omega_cross_r1[2];

        cv2[0] = v[j][0] + omega_cross_r2[0];
        cv2[1] = v[j][1] + omega_cross_r2[1];
        cv2[2] = v[j][2] + omega_cross_r2[2];

        // total relavtive velocity at contact point
        vr1 = cv1[0] - cv2[0];
        vr2 = cv1[1] - cv2[1];
        vr3 = cv1[2] - cv2[2];

        // normal component

        vn1 = nij[0] * vr1;    // dot product
        vn2 = nij[1] * vr2;
        vn3 = nij[2] * vr3;

        vnnr = vr1 * nij[0] + vr2 * nij[1] + vr3 * nij[2];    // magnitude

        // tangential component

        vtr1 = vr1 - vnnr * nij[0];
        vtr2 = vr2 - vnnr * nij[1];
        vtr3 = vr3 - vnnr * nij[2];

        vrel = vtr1 * vtr1 + vtr2 * vtr2 + vtr3 * vtr3;
        vrel = sqrt(vrel);

        // meff = effective mass of pair of particles
        // if I or J part of rigid body, use body mass
        // if I or J is frozen, meff is other particle

        mi = rmass[i];
        mj = rmass[j];
        if (fix_rigid) {
          if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
          if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
        }

        meff = mi * mj / (mi + mj);
        if (mask[i] & freeze_group_bit) meff = mj;
        if (mask[j] & freeze_group_bit) meff = mi;

        // normal forces = Hertzian contact + normal velocity damping

        damp = meff * gamman * vnnr;
        ccel = kn * (overlap1 + overlap2) + damp;    // assuming we get the overlap depth
        ccel *= polyhertz;
        if (limit_damping && (ccel < 0.0)) ccel = 0.0;

        // shear history effects

        touch[jj] = 1;
        shear = &allhistory[size_history * jj];

        if (shearupdate) {
          shear[0] += vtr1 * dt;
          shear[1] += vtr2 * dt;
          shear[2] += vtr3 * dt;
        }
        shrmag = sqrt(shear[0] * shear[0] + shear[1] * shear[1] + shear[2] * shear[2]);

        if (shearupdate) {

          // rotate shear displacements

          rsht = shear[0] * nij[0] + shear[1] * nij[1] + shear[2] * nij[2];
          shear[0] -= rsht * nij[0];
          shear[1] -= rsht * nij[1];
          shear[2] -= rsht * nij[2];
        }

        // tangential forces = shear + tangential velocity damping

        fs1 = -polyhertz * (kt * shear[0] + meff * gammat * vtr1);
        fs2 = -polyhertz * (kt * shear[1] + meff * gammat * vtr2);
        fs3 = -polyhertz * (kt * shear[2] + meff * gammat * vtr3);

        // rescale frictional displacements and forces if needed

        fs = sqrt(fs1 * fs1 + fs2 * fs2 + fs3 * fs3);
        fn = xmu * fabs(ccel);

        if (fs > fn) {
          if (shrmag != 0.0) {
            shear[0] =
                (fn / fs) * (shear[0] + meff * gammat * vtr1 / kt) - meff * gammat * vtr1 / kt;
            shear[1] =
                (fn / fs) * (shear[1] + meff * gammat * vtr2 / kt) - meff * gammat * vtr2 / kt;
            shear[2] =
                (fn / fs) * (shear[2] + meff * gammat * vtr3 / kt) - meff * gammat * vtr3 / kt;
            fs1 *= fn / fs;
            fs2 *= fn / fs;
            fs3 *= fn / fs;
          } else
            fs1 = fs2 = fs3 = 0.0;
        }

        // forces & torques

        fx = nji[0] * ccel + fs1;
        fy = nji[1] * ccel + fs2;
        fz = nji[2] * ccel + fs3;
        fx *= factor_lj;    // I think factor lj is just 1 except for special bonds
        fy *= factor_lj;
        fz *= factor_lj;
        f[i][0] += fx;
        f[i][1] += fy;
        f[i][2] += fz;

        // torques are cross prodcuts of branch vector with the entire force at contact point

        tor1 = cr1[1] * fz - cr1[2] * fy;
        tor2 = cr1[2] * fx - cr1[0] * fz;
        tor3 = cr1[0] * fy - cr1[1] * fx;

        torque[i][0] += tor1;
        torque[i][1] += tor2;
        torque[i][2] += tor3;

        if (newton_pair || j < nlocal) {
          f[j][0] -= fx;
          f[j][1] -= fy;
          f[j][2] -= fz;

          tor1 = cr2[1] * fz - cr2[2] * fy;
          tor2 = cr2[2] * fx - cr2[0] * fz;
          tor3 = cr2[0] * fy - cr2[1] * fx;

          torque[j][0] -= tor1;
          torque[j][1] -= tor2;
          torque[j][2] -= tor3;
        }

        if (evflag)
          ev_tally_xyz(i, j, nlocal, newton_pair, 0.0, 0.0, fx, fy, fz, delx, dely,
                       delz);    // Correct even for non-spherical particles
      }
    }
  }

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairGranHertzHistoryEllipsoid::settings(int narg, char **arg)
{
  if (narg < 6) error->all(FLERR, "Illegal pair_style command");

  kn = utils::numeric(FLERR, arg[0], false, lmp);
  if (strcmp(arg[1], "NULL") == 0)
    kt = kn * 2.0 / 7.0;
  else
    kt = utils::numeric(FLERR, arg[1], false, lmp);

  gamman = utils::numeric(FLERR, arg[2], false, lmp);
  if (strcmp(arg[3], "NULL") == 0)
    gammat = 0.5 * gamman;
  else
    gammat = utils::numeric(FLERR, arg[3], false, lmp);

  xmu = utils::numeric(FLERR, arg[4], false, lmp);
  dampflag = utils::inumeric(FLERR, arg[5], false, lmp);
  if (dampflag == 0) gammat = 0.0;

  limit_damping = 0;
  bounding_box = 0;
  curvature_model = MathExtraSuperellipsoids::CURV_MEAN;    // Default to Mean curvature

  for (int iarg = 6; iarg < narg; iarg++) {
    if (strcmp(arg[iarg], "limit_damping") == 0)
      limit_damping = 1;
    else if (strcmp(arg[iarg], "bounding_box") == 0)
      bounding_box = 1;
    else if (strcmp(arg[iarg], "geometric") == 0)
      contact_formulation = MathExtraSuperellipsoids::FORMULATION_GEOMETRIC;
    else if (strcmp(arg[iarg], "curvature_gaussian") == 0)
      curvature_model = MathExtraSuperellipsoids::CURV_GAUSSIAN;
    else
      error->all(FLERR, "Illegal pair_style command");
  }

  size_history = 8;    // reset to default size
  if (bounding_box == 0) size_history--;

  if (kn < 0.0 || kt < 0.0 || gamman < 0.0 || gammat < 0.0 || xmu < 0.0 || xmu > 10000.0 ||
      dampflag < 0 || dampflag > 1)
    error->all(FLERR, "Illegal pair_style command");

  // convert Kn and Kt from pressure units to force/distance^2

  kn /= force->nktv2p;
  kt /= force->nktv2p;
}

/* ---------------------------------------------------------------------- */

double PairGranHertzHistoryEllipsoid::single(int i, int j, int /*itype*/, int /*jtype*/, double rsq,
                                             double /*factor_coul*/, double /*factor_lj*/,
                                             double &fforce)
{
  double radi, radj, radsum;
  double vr1, vr2, vr3, vnnr, vn1, vn2, vn3, vt1, vt2, vt3;
  double mi, mj, meff, damp, ccel;
  double vtr1, vtr2, vtr3, vrel, shrmag, polyhertz;
  double fs1, fs2, fs3, fs, fn;

  double *radius = atom->radius;
  radi = radius[i];
  radj = radius[j];
  radsum = radi + radj;

  double **x = atom->x;

  // history effects
  // neighprev = index of found neigh on previous call
  // search entire jnum list of neighbors of I for neighbor J
  // start from neighprev, since will typically be next neighbor
  // reset neighprev to 0 as necessary
  int jnum = list->numneigh[i];
  int *jlist = list->firstneigh[i];
  int *touch = fix_history->firstflag[i];
  double *allhistory = fix_history->firstvalue[i];
  for (int jj = 0; jj < jnum; jj++) {
    neighprev++;
    if (neighprev >= jnum) neighprev = 0;
    if (jlist[neighprev] == j) break;
  }

  if (rsq >= radsum * radsum) {
    fforce = 0.0;
    for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
    return 0.0;
  }
  auto avec_ellipsoid = dynamic_cast<AtomVecEllipsoid *>(atom->style_match("ellipsoid"));
  AtomVecEllipsoid::BonusSuper *bonus = avec_ellipsoid->bonus_super;
  int *ellipsoid = atom->ellipsoid;
  double shapei[3], blocki[3], shapej[3], blockj[3], Ri[3][3], Rj[3][3];
  MathExtra::copy3(bonus[ellipsoid[i]].shape, shapei);
  MathExtra::copy3(bonus[ellipsoid[j]].shape, shapej);
  MathExtra::copy3(bonus[ellipsoid[i]].block, blocki);
  MathExtra::copy3(bonus[ellipsoid[j]].block, blockj);
  MathExtra::quat_to_mat(bonus[ellipsoid[i]].quat, Ri);
  MathExtra::quat_to_mat(bonus[ellipsoid[j]].quat, Rj);
  bool skip_contact_detection(false);
  if (bounding_box) {
    int cached_axis = (int)
        (allhistory[7 + size_history * neighprev]);    // Copy: no update of history in single
    int new_axis = MathExtraSuperellipsoids::check_oriented_bounding_boxes(
        x[i], Ri, shapei, x[j], Rj, shapej, cached_axis);
    if (new_axis !=-1) skip_contact_detection = true;
    if (skip_contact_detection) {
      fforce = 0.0;
      for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
      return 0.0;
    }
  }
  // superellipsoid contact detection between atoms i and j
  double X0[4], nij[3];
  AtomVecEllipsoid::BlockType flagi, flagj;
  flagi = bonus[ellipsoid[i]].type;
  flagj = bonus[ellipsoid[j]].type;
  double *X0_prev = &allhistory[3 + size_history * neighprev];
  if (touch[neighprev] == 1) {
    int ref_index = (atom->tag[i] < atom->tag[j]) ? i : j;
    // Continued contact: use grain true shape and last contact point
    X0[0] = X0_prev[0] + x[ref_index][0];
    X0[1] = X0_prev[1] + x[ref_index][1];
    X0[2] = X0_prev[2] + x[ref_index][2];
    X0[3] = X0_prev[3];
    int status = MathExtraSuperellipsoids::determine_contact_point(x[i], Ri, shapei, blocki, flagi,
                                                                   x[j], Rj, shapej, blockj, flagj,
                                                                   X0, nij, contact_formulation);
    if (status == 1) {
      fforce = 0.0;
      for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
      return 0.0;
    }
    if (status != 0)
      error->all(FLERR,
                 "Ellipsoid contact detection (old contact) failed"
                 "between particle {} and particle {}",
                 atom->tag[i], atom->tag[j]);
  } else {
    double reqi = std::cbrt(shapei[0] * shapei[1] * shapei[2]);
    double reqj = std::cbrt(shapej[0] * shapej[1] * shapej[2]);
    MathExtra::scaleadd3(reqj / (reqi + reqj), x[i], reqi / (reqi + reqj), x[j], X0);
    X0[3] = reqj / reqi;    // Lagrange multiplier mu^2
    for (int iter_ig = 1; iter_ig <= NUMSTEP_INITIAL_GUESS; iter_ig++) {
      double frac = iter_ig / double(NUMSTEP_INITIAL_GUESS);
      shapei[0] = shapei[1] = shapei[2] = reqi;
      shapej[0] = shapej[1] = shapej[2] = reqj;
      MathExtra::scaleadd3(1.0 - frac, shapei, frac, bonus[ellipsoid[i]].shape, shapei);
      MathExtra::scaleadd3(1.0 - frac, shapej, frac, bonus[ellipsoid[j]].shape, shapej);
      blocki[0] = 2.0 + frac * (bonus[ellipsoid[i]].block[0] - 2.0);
      blocki[1] = 2.0 + frac * (bonus[ellipsoid[i]].block[1] - 2.0);
      blockj[0] = 2.0 + frac * (bonus[ellipsoid[j]].block[0] - 2.0);
      blockj[1] = 2.0 + frac * (bonus[ellipsoid[j]].block[1] - 2.0);

      // force ellipsoid flag for first initial guess iteration.
      // Avoid incorrect values of n1/n2 - 2 in second derivatives.
      int status = MathExtraSuperellipsoids::determine_contact_point(
          x[i], Ri, shapei, blocki, iter_ig == 1 ? AtomVecEllipsoid::BlockType::ELLIPSOID : flagi,
          x[j], Rj, shapej, blockj, iter_ig == 1 ? AtomVecEllipsoid::BlockType::ELLIPSOID : flagj,
          X0, nij, contact_formulation);
      if (status == 1) {
        fforce = 0.0;
        for (int m = 0; m < single_extra; m++) svector[m] = 0.0;
        return 0.0;
      }
      if (status != 0)
        error->all(FLERR,
                   "Ellipsoid contact detection (new contact) failed"
                   "between particle {} and particle {}",
                   atom->tag[i], atom->tag[j]);
    }
  }
  double overlap1, overlap2, omegai[3], omegaj[3];
  double nji[3] = {-nij[0], -nij[1], -nij[2]};
  overlap1 =
      MathExtraSuperellipsoids::compute_overlap_distance(shapei, blocki, Ri, flagi, X0, nij, x[i]);
  overlap2 =
      MathExtraSuperellipsoids::compute_overlap_distance(shapej, blockj, Rj, flagj, X0, nji, x[j]);

  double surf_point_i[3], surf_point_j[3], curvature_i, curvature_j;
  MathExtra::scaleadd3(overlap1, nij, X0, surf_point_i);
  MathExtra::scaleadd3(overlap2, nji, X0, surf_point_j);

  if (curvature_model == MathExtraSuperellipsoids::CURV_MEAN) {
    curvature_i = MathExtraSuperellipsoids::mean_curvature_superellipsoid(shapei, blocki, flagi, Ri,
                                                                          surf_point_i, x[i]);
    curvature_j = MathExtraSuperellipsoids::mean_curvature_superellipsoid(shapej, blockj, flagj, Rj,
                                                                          surf_point_j, x[j]);
  } else {
    curvature_i = MathExtraSuperellipsoids::gaussian_curvature_superellipsoid(
        shapei, blocki, flagi, Ri, surf_point_i, x[i]);
    curvature_j = MathExtraSuperellipsoids::gaussian_curvature_superellipsoid(
        shapej, blockj, flagj, Rj, surf_point_j, x[j]);
  }

  polyhertz = sqrt((overlap1 + overlap2) /
                   (curvature_i + curvature_j));    // hertzian contact radius approximation

  double cr1[3], cr2[3];
  MathExtra::sub3(X0, x[i], cr1);
  MathExtra::sub3(X0, x[j], cr2);

  double ex_space[3], ey_space[3], ez_space[3];
  double **angmom = atom->angmom;
  MathExtra::q_to_exyz(bonus[ellipsoid[i]].quat, ex_space, ey_space, ez_space);
  MathExtra::angmom_to_omega(angmom[i], ex_space, ey_space, ez_space, bonus[ellipsoid[i]].inertia,
                             omegai);
  MathExtra::q_to_exyz(bonus[ellipsoid[j]].quat, ex_space, ey_space, ez_space);
  MathExtra::angmom_to_omega(angmom[j], ex_space, ey_space, ez_space, bonus[ellipsoid[j]].inertia,
                             omegaj);

  double omega_cross_r1[3], omega_cross_r2[3];
  MathExtra::cross3(omegai, cr1, omega_cross_r1);
  MathExtra::cross3(omegaj, cr2, omega_cross_r2);

  // relative translational velocity
  // compute directly the sum of relative translational velocity at contact point
  // since rotational velocity contribution is different for superellipsoids

  double **v = atom->v;
  double cv1[3], cv2[3];

  cv1[0] = v[i][0] + omega_cross_r1[0];
  cv1[1] = v[i][1] + omega_cross_r1[1];
  cv1[2] = v[i][2] + omega_cross_r1[2];

  cv2[0] = v[j][0] + omega_cross_r2[0];
  cv2[1] = v[j][1] + omega_cross_r2[1];
  cv2[2] = v[j][2] + omega_cross_r2[2];

  // total relavtive velocity at contact point

  vr1 = cv1[0] - cv2[0];
  vr2 = cv1[1] - cv2[1];
  vr3 = cv1[2] - cv2[2];

  // normal component

  vn1 = nij[0] * vr1;    // dot product
  vn2 = nij[1] * vr2;
  vn3 = nij[2] * vr3;

  vnnr = vr1 * nij[0] + vr2 * nij[1] + vr3 * nij[2];    // magnitude

  // tangential component

  vtr1 = vr1 - vnnr * nij[0];
  vtr2 = vr2 - vnnr * nij[1];
  vtr3 = vr3 - vnnr * nij[2];

  vrel = vtr1 * vtr1 + vtr2 * vtr2 + vtr3 * vtr3;
  vrel = sqrt(vrel);

  // meff = effective mass of pair of particles
  // if I or J part of rigid body, use body mass
  // if I or J is frozen, meff is other particle
  double *rmass = atom->rmass;
  int *mask = atom->mask;

  mi = rmass[i];
  mj = rmass[j];
  if (fix_rigid) {
    if (mass_rigid[i] > 0.0) mi = mass_rigid[i];
    if (mass_rigid[j] > 0.0) mj = mass_rigid[j];
  }

  meff = mi * mj / (mi + mj);
  if (mask[i] & freeze_group_bit) meff = mj;
  if (mask[j] & freeze_group_bit) meff = mi;

  // normal forces = Hookian contact + normal velocity damping

  damp = meff * gamman * vnnr;
  ccel = kn * (overlap1 + overlap2) + damp;    // assuming we get the overlap depth
  ccel *= polyhertz;
  if (limit_damping && (ccel < 0.0)) ccel = 0.0;

  double *shear = &allhistory[size_history * neighprev];
  shrmag = sqrt(shear[0] * shear[0] + shear[1] * shear[1] + shear[2] * shear[2]);

  // tangential forces = shear + tangential velocity damping

  fs1 = -polyhertz * (kt * shear[0] + meff * gammat * vtr1);
  fs2 = -polyhertz * (kt * shear[1] + meff * gammat * vtr2);
  fs3 = -polyhertz * (kt * shear[2] + meff * gammat * vtr3);

  // rescale frictional displacements and forces if needed

  fs = sqrt(fs1 * fs1 + fs2 * fs2 + fs3 * fs3);
  fn = xmu * fabs(ccel);

  if (fs > fn) {
    if (shrmag != 0.0) {
      fs1 *= fn / fs;
      fs2 *= fn / fs;
      fs3 *= fn / fs;
      fs *= fn / fs;
    } else
      fs1 = fs2 = fs3 = 0.0;
  }

  // set force (normalized by r) and return no energy

  fforce = ccel / sqrt(rsq);

  // set single_extra quantities

  svector[0] = fs1;
  svector[1] = fs2;
  svector[2] = fs3;
  svector[3] = fs;
  svector[4] = vn1;
  svector[5] = vn2;
  svector[6] = vn3;
  svector[7] = vtr1;
  svector[8] = vtr2;
  svector[9] = vtr3;

  return 0.0;
}
