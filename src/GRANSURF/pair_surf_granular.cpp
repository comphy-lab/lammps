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

#include "pair_surf_granular.h"

#include "atom.h"
#include "atom_vec_line.h"
#include "atom_vec_tri.h"
#include "comm.h"
#include "granular_model.h"
#include "gran_sub_mod.h"
#include "error.h"
#include "fix.h"
#include "fix_dummy.h"
#include "fix_neigh_history.h"
#include "fix_surface_local.h"
#include "fix_surface.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "neigh_list.h"
#include "neighbor.h"
#include "surf_extra.h"
#include "update.h"

#include <algorithm>
#include <map>
#include <unordered_set>
#include <vector>

using namespace LAMMPS_NS;
using namespace Granular_NS;
using namespace MathExtra;
using namespace SurfExtra;

enum {NONE, LINE, TRI};
enum{FLAT,FLATEDGE,CONCAVE,CONVEX};
enum{SAME_SIDE,OPPOSITE_SIDE};
enum{INTERNAL = 0, NONFLAT = 1, EXTERNAL = 2};

static constexpr int DELTACONTACTS = 4;
static constexpr double EPSILON = 1e-14;

static inline int FLIPSIDE(int nside) {
  if (nside == OPPOSITE_SIDE) return SAME_SIDE;
  else return OPPOSITE_SIDE;
}

/* ---------------------------------------------------------------------- */

PairSurfGranular::PairSurfGranular(LAMMPS *lmp) : PairGranular(lmp)
{
  single_enable = 0;

  emax = 0;
  endpts = nullptr;
  cmax = 0;
  corners = nullptr;

  contact_surfs = nullptr;
}

/* ---------------------------------------------------------------------- */

PairSurfGranular::~PairSurfGranular()
{
  memory->destroy(endpts);
  memory->destroy(corners);

  memory->sfree(contact_surfs);
}

/* ---------------------------------------------------------------------- */

void PairSurfGranular::compute(int eflag, int vflag)
{
  int i, j, k, a, n, m, nconnect, ii, jj, inum, jnum, itype, jtype;
  int isphere, itri, jflag, kflag, n_contact_surfs, exposed_flag, zero_overlap;
  double xtmp, ytmp, ztmp, radi, delx, dely, delz;
  double rsq, radsum, max_overlap, tmp_max, dot, distance_from_surf;
  double factor_lj, mi, mj, meff;
  double norm[3], dr[3], contact[3], ds[3], xc[3], vc[3], omegac[3], residual[3];
  double *endpt, *corner, *forces, *torquesi, *torquesj, dq;
  double omega0[3] = {0.0, 0.0, 0.0};

  int it, jjtmp, nsidej;
  std::vector<int> *composite_surfs = new std::vector<int>();
  std::unordered_set<int> *processed_contacts = new std::unordered_set<int>();
  std::unordered_set<int> *convex_contacts = new std::unordered_set<int>();
  std::unordered_set<int> *concave_contacts = new std::unordered_set<int>();
  std::map<int, int> *contacts_map = new std::map<int, int>();

  int *ilist, *jlist, *numneigh, **firstneigh;
  int *touch, **firsttouch;
  double *history, *allhistory, **firsthistory;

  bool touchflag = false;
  const bool history_update = update->setupflag == 0;

  class GranularModel* model;

  for (int i = 0; i < nmodels; i++)
    models_list[i]->history_update = history_update;

  ev_init(eflag,vflag);

  // if just reneighbored:
  // update rigid body info for owned & ghost atoms if using FixRigid masses
  // body[i] = which body atom I is in, -1 if none
  // mass_body = mass of each rigid body
  // forward comm mass_rigid so have it for ghost lines
  // also grab current line connectivity info from FixSurfaceLocal

  if (neighbor->ago == 0) {
    if (fix_rigid) {
      int tmp;
      int *body = (int *) fix_rigid->extract("body",tmp);
      double *mass_body = (double *) fix_rigid->extract("masstotal",tmp);
      if (atom->nmax > nmax) {
        memory->destroy(mass_rigid);
        nmax = atom->nmax;
        memory->create(mass_rigid,nmax,"surf/granular:mass_rigid");
      }
      int nlocal = atom->nlocal;
      for (i = 0; i < nlocal; i++)
        if (body[i] >= 0)
          mass_rigid[i] = mass_body[body[i]];
        else
          mass_rigid[i] = 0.0;
      comm->forward_comm(this);
    }

    connect2d = fsl->connect2d;
    connect3d = fsl->connect3d;
  }

  // pre-calculate current end pts of owned+ghost lines
  // only once per reneighbor if surfs not moving

  if (surfmoveflag || neighbor->ago == 0) {
    if (style == LINE) calculate_endpts();
    if (style == TRI) calculate_corners();
  }

  // loop over neighbors of my atoms
  // I is always sphere, J is always line/tri

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *type = atom->type;
  double **omega = atom->omega;
  double **torque = atom->torque;
  double *radius = atom->radius;
  double *rmass = atom->rmass;
  tagint *tag = atom->tag;
  int *line = atom->line;
  int *tri = atom->tri;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  int newton_pair = force->newton_pair;
  double *special_lj = force->special_lj;
  double *heatflow, *temperature;

  if (heat_flag) {
    heatflow = atom->heatflow;
    temperature = atom->temperature;
  }

  inum = list->inum;
  ilist = list->ilist;
  numneigh = list->numneigh;
  firstneigh = list->firstneigh;
  if (use_history) {
    firsttouch = fix_history->firstflag;
    firsthistory = fix_history->firstvalue;
  }

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];
    itype = type[i];
    xtmp = x[i][0];
    ytmp = x[i][1];
    ztmp = x[i][2];
    radi = radius[i];

    if (use_history) {
      touch = firsttouch[i];
      allhistory = firsthistory[i];
    }

    jlist = firstneigh[i];
    jnum = numneigh[i];

    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      // Do we need special bonds?
      // factor_lj = special_lj[sbmask(j)];
      // j &= NEIGHMASK;
      // if (factor_lj == 0) continue;

      // sanity check that neighbor list is built correctly

      if ((style == LINE) && (line[i] >= 0 || line[j] < 0))
        error->one(FLERR,"Pair surf/granular interaction is invalid");

      if ((style == TRI) && (tri[i] >= 0 || tri[j] < 0))
        error->one(FLERR,"Pair surf/granular interaction is invalid");

      delx = xtmp - x[j][0];
      dely = xtmp - x[j][1];
      delz = xtmp - x[j][2];
      rsq = delx * delx + dely * dely + delz * delz;

      // skip contact check if particle/surf are too far apart

      radsum = radi + radius[j];
      if (rsq > radsum * radsum) {
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history * jj];
          for (k = 0; k < size_history; k++) history[k] = 0.0;
        }
        continue;
      }

      // check for overlap of sphere and line segment/triangle
      // for line:
      //   jflag = 0 for no overlap, 1 for interior line pt, -1/-2 for end pts
      // for tri:
      //   jflag = 0 for no overlap, 1 for interior line pt,
      //     -1/-2/-3 for 3 edges, -4/-5/-6 for 3 corner pts
      // if no overlap, just continue
      // for overlap, also return:
      //   contact = nearest point on line/tri to sphere center
      //   dr = vector from contact pt to sphere center
      //   rsq = squared length of dr

      if (style == LINE) {
        endpt = endpts[line[j]];
        jflag = SurfExtra::
          overlap_sphere_line(x[i], radi, &endpt[0], &endpt[3],contact, dr, rsq);

      } else if (style == TRI) {
        corner = corners[tri[j]];
        jflag = SurfExtra::
          overlap_sphere_tri(x[i], radi, &corner[0], &corner[3], &corner[6], &corner[9],
                             contact, dr, rsq);
      }

      // unset non-touching neighbors

      if (!jflag) {
        if (use_history) {
          touch[jj] = 0;
          history = &allhistory[size_history * jj];
          for (k = 0; k < size_history; k++) history[k] = 0.0;
        }
        continue;
      }

      // append surf to list of contacts

      if (n_contact_surfs + 1 > nmax_contact_surfs) {
        nmax_contact_surfs += DELTACONTACTS;
        memory->grow(contact_surfs, nmax_contact_surfs * sizeof(FixSurface::ContactSurf),
                                      "surface/global:contact_surfs");
      }

      // Store which side is in contact relative to normal vector
      exposed_flag = 0;
      if (style == LINE) {
        MathExtra::copy3(&endpts[line[j]][6], norm);
        dot = MathExtra::dot3(norm, dr);
        if (jflag == -1) exposed_flag = connect2d[j].exposed_pt[0];
        if (jflag == -2) exposed_flag = connect2d[j].exposed_pt[1];
      } else {
        MathExtra::copy3(&corners[tri[j]][9], norm);
        dot = MathExtra::dot3(norm, dr);
        if (jflag == -1) exposed_flag = connect3d[j].exposed_edge[0];
        if (jflag == -2) exposed_flag = connect3d[j].exposed_edge[1];
        if (jflag == -3) exposed_flag = connect3d[j].exposed_edge[2];
        if (jflag == -4) exposed_flag = connect3d[j].exposed_pt[0];
        if (jflag == -5) exposed_flag = connect3d[j].exposed_pt[1];
        if (jflag == -6) exposed_flag = connect3d[j].exposed_pt[2];
      }

      if (dot >= 0) nsidej = SAME_SIDE;
      else nsidej = OPPOSITE_SIDE;

      contact_surfs[n_contact_surfs].index = j;
      contact_surfs[n_contact_surfs].neigh_index = jj;
      contact_surfs[n_contact_surfs].type = type[j];
      contact_surfs[n_contact_surfs].flag = jflag;
      contact_surfs[n_contact_surfs].exposed = exposed_flag;
      contact_surfs[n_contact_surfs].nside = nsidej;
      contact_surfs[n_contact_surfs].overlap = radi - sqrt(rsq);
      contact_surfs[n_contact_surfs].norm_def = -1;
      contact_surfs[n_contact_surfs].smooth_ext = 0;
      MathExtra::zero3(contact_surfs[n_contact_surfs].force_norm);
      MathExtra::copy3(norm, contact_surfs[n_contact_surfs].surf_norm);
      MathExtra::copy3(dr, contact_surfs[n_contact_surfs].dr);
      MathExtra::copy3(contact, contact_surfs[n_contact_surfs].contact);
      MathExtra::zero3(contact_surfs[n_contact_surfs].cor_int);
      MathExtra::zero3(contact_surfs[n_contact_surfs].cor_ext);

      // Ensure interior contacts always win in a tie, needed for convex flat structures
      //   that calculate distance to the corner to smooth turning
      if (jflag == 1)
        contact_surfs[n_contact_surfs].overlap += EPSILON;

      n_contact_surfs += 1;
    }

    // Reduce set of contacts

    /*
    For contact in reduced contacts:

      factor_lj = special_lj[sbmask(j)]; // presumably not necessary
      j &= NEIGHMASK;

      if (factor_lj == 0) continue;

      // if any history is needed

      if (use_history) touch[jj] = 1;

      // calculate new data
      // ds = vector from line/tri center to contact pt
      // vs = velocity of contact pt on line/tri, translation + rotation
      // omega for tri was set from angmom by calculate_corners()



      jtype = type[j];
      model = models_list[types_indices[itype][jtype]];

      // reset model and copy initial geometric data

      model->xi = x[i];
      model->xj = x[j];
      model->radi = radius[i];
      model->radj = radius[j];
      if (use_history) model->touch = touch[jj];

      model->rsq = rsq;

      ds[0] = contact[0] - x[j][0];
      ds[1] = contact[1] - x[j][1];
      ds[2] = contact[2] - x[j][2];

      vs[0] = v[j][0] + (omega[j][1] * ds[2] - omega[j][2] * ds[1]);
      vs[1] = v[j][1] + (omega[j][2] * ds[0] - omega[j][0] * ds[2]);
      vs[2] = v[j][2] + (omega[j][0] * ds[1] - omega[j][1] * ds[0]);

      model->dx[0] = dr[0];
      model->dx[1] = dr[1];
      model->dx[2] = dr[2];

      // NOTE: add logic to persist history history if contact has changed

      // NOTE: add logic to check for coupled contacts and weight them

      factor_couple = 1.0;

      // meff = effective mass of sphere and line/tri
      // if I or J is part of rigid body,use body mass
      // if line/tri is not part of rigidbody assume infinite mass

      meff = rmass[i];
      if (fix_rigid) {
        if (mass_rigid[i] > 0.0) meff = mass_rigid[i];
        if (mass_rigid[j] > 0.0) {
          mj = mass_rigid[j];
          meff = meff * mj / (meff + mj);
        }
      }

      // Copy additional information and prepare force calculations
      model->meff = meff;
      model->vi = v[i];
      model->vj = vs;
      model->omegai = omega[i];
      model->omegaj = omega0;

      if (use_history) {
        history = &allhistory[size_history * jj];
        model->history = history;
      }

      if (heat_flag) {
        model->Ti = temperature[i];
        model->Tj = temperature[j];
      }

      model->calculate_forces();

      // need to add support coupled contacts
      // is this just multiplying forces (+torques?) by factor_couple?

      forces = model->forces;
      torquesi = model->torquesi;
      torquesj = model->torquesj;

      // apply forces & torques
      scale3(factor_lj, forces);
      add3(f[i], forces, f[i]);

      scale3(factor_lj, torquesi);
      add3(torque[i], torquesi, torque[i]);

      if (force->newton_pair || j < nlocal) {
        sub3(f[j], forces, f[j]);
        scale3(factor_lj, torquesj);
        add3(torque[j], torquesj, torque[j]);
      }

      if (heat_flag) {
        dq = model->dq;
        heatflow[i] += dq;
        if (force->newton_pair || j < nlocal) heatflow[j] -= dq;
      }

      if (evflag) {
        ev_tally_xyz(i,j,nlocal,force->newton_pair,
          0.0,0.0,forces[0],forces[1],forces[2],model->dx[0],model->dx[1],model->dx[2]);
      }
    */
  }

  // NOTE: should there be virial contributions from boundary tris?

  if (vflag_fdotr) virial_fdotr_compute();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSurfGranular::init_style()
{
  // error and warning checks

  avecline = (AtomVecLine *) atom->style_match("line");
  avectri = (AtomVecTri *) atom->style_match("tri");

  if (avecline) style = LINE;
  if (avectri) style = TRI;

  if (style == NONE)
    error->all(FLERR,"Pair surf/granular requires atom style line or tri");

  if (!atom->radius_flag || !atom->rmass_flag || !atom->omega_flag)
    error->all(FLERR, "Pair surf/granular requires atom attributes radius, rmass, omega");
  if (!force->newton_pair)
    error->all(FLERR,"Pair style surf/granular requires newton pair on");
  if (comm->ghost_velocity == 0)
    error->all(FLERR,"Pair surf/granular requires ghost atoms store velocity");

  if (heat_flag) {
    if (!atom->temperature_flag)
      error->all(FLERR,"Heat conduction in pair surf/granularular requires atom style with temperature property");
    if (!atom->heatflow_flag)
      error->all(FLERR,"Heat conduction in pair surf/granularular requires atom style with heatflow property");
  }

  // allocate history and initialize models
  class GranularModel* model;
  int size_max[NSUBMODELS] = {0};
  for (int n = 0; n < nmodels; n++) {
    model = models_list[n];
    model->contact_type = SURFACE;

    if (model->beyond_contact) {
      beyond_contact = 1;
      use_history = 1; // Need to track if in contact
    }
    if (model->size_history != 0) use_history = 1;

    for (int i = 0; i < NSUBMODELS; i++)
      if (model->sub_models[i]->size_history > size_max[i])
        size_max[i] = model->sub_models[i]->size_history;

    if (model->nondefault_history_transfer) nondefault_history_transfer = 1;
  }

  size_history = 0;
  if (use_history) {
    for (int i = 0; i < NSUBMODELS; i++) size_history += size_max[i];

    // Ensure size history is at least 1 to avoid errors in fix neigh/history
    // This could occur if normal model is beyond_contact but no other entries are required
    // E.g. JKR + linear_nohistory
    size_history = MAX(size_history, 1);
  }

  for (int n = 0; n < nmodels; n++) {
    model = models_list[n];
    int next_index = 0;
    for (int i = 0; i < NSUBMODELS; i++) {
      model->sub_models[i]->history_index = next_index;
      next_index += size_max[i];
    }

    if (model->beyond_contact)
      error->all(FLERR, "Beyond contact models not currenty supported");
  }

  // need a granular neighbor list

  if (use_history)
    neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_ONESIDED |
                          NeighConst::REQ_HISTORY);
  else
    neighbor->add_request(this, NeighConst::REQ_SIZE | NeighConst::REQ_ONESIDED);

  // if history is stored and first init, create Fix to store history
  // it replaces FixDummy, created in the constructor
  // this is so its order in the fix list is preserved

  if (use_history && fix_history == nullptr) {
    fix_history = dynamic_cast<FixNeighHistory *>(modify->replace_fix(id_dummy, fmt::format("{} all NEIGH_HISTORY {} onesided", id_history, size_history),1));
    fix_history->pair = this;
  } else if (use_history) {
    fix_history = dynamic_cast<FixNeighHistory *>(modify->get_fix_by_id(id_history));
    if (!fix_history) error->all(FLERR,"Could not find pair fix neigh history ID");
  }

  // set ptr to FixSurfaceLocal for surf connectivity info

  fsl = nullptr;
  for (int m = 0; m < modify->nfix; m++) {
    if (strcmp(modify->fix[m]->style,"surface/local") == 0) {
      if (fsl)
        error->all(FLERR,"Pair surf/granular requires single fix surface/local");
      fsl = (FixSurfaceLocal *) modify->fix[m];
    }
  }
  if (!fsl) error->all(FLERR,"Pair surf/granular requires a fix surface/local");

  // surfmoveflag = 1 if surfs may move at every step
  // yes if fix move exists and its group includes lines
  // NOTE: are there other conditions, like fix deform or fix npt?

  surfmoveflag = 0;
  for (int m = 0; m < modify->nfix; m++) {
    if (strcmp(modify->fix[m]->style,"move") == 0) {
      int groupbit = modify->fix[m]->groupbit;
      int *mask = atom->mask;
      int nlocal = atom->nlocal;
      int flag = 0;
      if (style == LINE) {
        int *line = atom->line;
        for (int i = 0; i < nlocal; i++) {
          if (line[i] < 0) continue;
          if (mask[i] & groupbit) flag = 1;
        }
      } else {
        int *tri = atom->tri;
        for (int i = 0; i < nlocal; i++) {
          if (tri[i] < 0) continue;
          if (mask[i] & groupbit) flag = 1;
        }
      }
      int any;
      MPI_Allreduce(&flag,&any,1,MPI_INT,MPI_SUM,world);
      if (any) surfmoveflag = 1;
    }
  }

  // check for FixFreeze and set freeze_group_bit

  auto fixlist = modify->get_fix_by_style("^freeze");
  if (fixlist.size() == 0)
    freeze_group_bit = 0;
  else if (fixlist.size() > 1)
    error->all(FLERR, "Only one fix freeze command at a time allowed");
  else
    freeze_group_bit = fixlist.front()->groupbit;

  // check for FixRigid so can extract rigid body masses

  fix_rigid = nullptr;
  for (const auto &ifix : modify->get_fix_list()) {
    if (ifix->rigid_flag) {
      if (fix_rigid)
        error->all(FLERR, "Only one fix rigid command at a time allowed");
      else
        fix_rigid = ifix;
    }
  }

  // check for FixPour and FixDeposit so can extract particle radii

  auto pours = modify->get_fix_by_style("^pour");
  auto deps = modify->get_fix_by_style("^deposit");

  // set maxrad_dynamic and maxrad_frozen for each type
  // include future FixPour and FixDeposit particles as dynamic
  // lines/tris cannot be frozen

  int itype;
  for (int i = 1; i <= atom->ntypes; i++) {
    onerad_dynamic[i] = onerad_frozen[i] = 0.0;
    for (auto &ipour : pours) {
      itype = i;
      double maxrad = *((double *) ipour->extract("radius", itype));
      if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
    }
    for (auto &idep : deps) {
      itype = i;
      double maxrad = *((double *) idep->extract("radius", itype));
      if (maxrad > 0.0) onerad_dynamic[i] = maxrad;
    }
  }

  double *radius = atom->radius;
  int *line = atom->line;
  int *tri = atom->tri;
  int *mask = atom->mask;
  int *type = atom->type;
  int nlocal = atom->nlocal;

  for (int i = 0; i < nlocal; i++) {
    if ((style == LINE) && line[i] >= 0)
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius[i]);
    else if ((style == TRI) && tri[i] >= 0)
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius[i]);
    else {
      if (mask[i] & freeze_group_bit)
        onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]],radius[i]);
      else
        onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]],radius[i]);
    }
  }

  MPI_Allreduce(&onerad_dynamic[1],&maxrad_dynamic[1],atom->ntypes,
                MPI_DOUBLE,MPI_MAX,world);
  MPI_Allreduce(&onerad_frozen[1],&maxrad_frozen[1],atom->ntypes,
                MPI_DOUBLE,MPI_MAX,world);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based arrays
------------------------------------------------------------------------- */

double PairSurfGranular::memory_usage()
{
  double bytes = nmax * sizeof(double);
  if (style == LINE) bytes = emax * 6 * sizeof(double);        // endpts array for line particles
  if (style == TRI) bytes = emax * 12 * sizeof(double);        // corners array for tri particles
  return bytes;
}

/* ----------------------------------------------------------------------
   compute current end points of owned and ghost lines
   nothing computed for particles that are not lines
------------------------------------------------------------------------- */

void PairSurfGranular::calculate_endpts()
{
  int i,m;
  double length,theta,dx,dy;
  double p12[3], norm[3];
  double zunit[3] = {0.0,0.0,1.0};
  double *endpt;

  // realloc endpts array if necssary

  if (fsl->nmax_connect > emax) {
    memory->destroy(endpts);
    emax = fsl->nmax_connect;
    memory->create(endpts,emax,9,"surf/granular:endpts");
  }

  AtomVecLine::Bonus *bonus = avecline->bonus;
  double **x = atom->x;
  int *line = atom->line;
  int n = atom->nlocal + atom->nghost;

  for (i = 0; i < n; i++) {
    if (line[i] < 0) continue;
    m = line[i];
    length = bonus[m].length;
    theta = bonus[m].theta;
    dx = 0.5*length*cos(theta);
    dy = 0.5*length*sin(theta);
    endpt = endpts[m];
    endpt[0] = x[i][0] - dx;
    endpt[1] = x[i][1] - dy;
    endpt[2] = 0.0;
    endpt[3] = x[i][0] + dx;
    endpt[4] = x[i][1] + dy;
    endpt[5] = 0.0;

    MathExtra::sub3(&endpt[3],&endpt[0],p12);
    MathExtra::cross3(zunit,p12,norm);
    MathExtra::normalize3(norm, &endpt[6]);
  }
}

/* ----------------------------------------------------------------------
   compute current corner points and current norm of N triangles
   also compute omega from angmom
   N = nlocal or nlocal+nghost atoms
   nothing computed for particles that are not tris
------------------------------------------------------------------------- */

void PairSurfGranular::calculate_corners()
{
  int i,m;
  double ex[3],ey[3],ez[3],p[3][3];
  double *corner;

  // realloc corners array if necssary

  if (fsl->nmax_connect > cmax) {
    memory->destroy(corners);
    cmax = fsl->nmax_connect;
    memory->create(corners,cmax,12,"surf/granular:corners");
  }

  AtomVecTri::Bonus *bonus = avectri->bonus;
  double **x = atom->x;
  double **omega = atom->omega;
  double **angmom = atom->angmom;
  int *tri = atom->tri;
  int n = atom->nlocal + atom->nghost;

  for (int i = 0; i < n; i++) {
    if (tri[i] < 0) continue;
    m = tri[i];
    corner = corners[m];
    MathExtra::quat_to_mat(bonus[m].quat,p);
    MathExtra::matvec(p,bonus[m].c1,&corner[0]);
    MathExtra::add3(x[i],&corner[0],&corner[0]);
    MathExtra::matvec(p,bonus[m].c2,&corner[3]);
    MathExtra::add3(x[i],&corner[3],&corner[3]);
    MathExtra::matvec(p,bonus[m].c3,&corner[6]);
    MathExtra::add3(x[i],&corner[6],&corner[6]);
    corners2norm(corner,&corner[9]);

    // omega from angmom of tri particles

    if (angmom[i][0] == 0.0 && angmom[i][1] == 0.0 && angmom[i][2] == 0.0) {
      omega[i][0] = omega[i][1] = omega[i][2] = 0.0;
      continue;
    }
    MathExtra::q_to_exyz(bonus[m].quat,ex,ey,ez);
    MathExtra::angmom_to_omega(angmom[i],ex,ey,ez,
                               bonus[m].inertia,omega[i]);
  }
}

/* ----------------------------------------------------------------------
   compute norm of a triangle based on its 3 corner pts
------------------------------------------------------------------------- */

void PairSurfGranular::corners2norm(double *corners, double *norm)
{
  double p12[3],p13[3];

  MathExtra::sub3(&corners[3],&corners[0],p12);
  MathExtra::sub3(&corners[6],&corners[0],p13);
  MathExtra::cross3(p12,p13,norm);
  MathExtra::norm3(norm);
}
