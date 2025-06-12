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
#include "domain.h"
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
static constexpr double BIG = 1.0e20;

static inline int FLIPSIDE(int nside) {
  if (nside == OPPOSITE_SIDE) return SAME_SIDE;
  else return OPPOSITE_SIDE;
}

static inline int EQUAL3(double *pt1, double *pt2) {
  int same = 1;

  if (fabs(pt1[0] - pt2[0]) > EPSILON) same = 0;
  else if (fabs(pt1[1] - pt2[1]) > EPSILON) same = 0;
  else if (fabs(pt1[2] - pt2[2]) > EPSILON) same = 0;

  return same;
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
  int dimension = domain->dimension;

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

    if (n_contact_surfs == 0)
      continue;

    // Sort contacts by overlap and create a map
    std::sort(contact_surfs, contact_surfs + n_contact_surfs,
      [](FixSurface::ContactSurf a, FixSurface::ContactSurf b) {return a.overlap > b.overlap;});
    contacts_map->clear();
    for (n = 0; n < n_contact_surfs; n++)
      (*contacts_map)[contact_surfs[n].index] = n;

    // Initial walk to assign consistent sides of surfaces
    //   Not guaranteed to work for v. complex geometries (e.g. Mobius)

    processed_contacts->clear();
    for (n = 0; n < n_contact_surfs; n++) {
      j = contact_surfs[n].index;
      if (processed_contacts->find(j) != processed_contacts->end())  continue;
      composite_surfs->clear();
      if (dimension == 2) {
        prewalk_connections2d(n, contact_surfs[n].nside,
                              processed_contacts, contacts_map);
      } else {
        prewalk_connections3d(n, contact_surfs[n].nside, composite_surfs,
                              processed_contacts, contacts_map);

        // Closest distance from sphere to surface
        distance_from_surf = BIG;
        if (contact_surfs[0].exposed == INTERNAL) {
          // If the closest point is internal, then the contact is internal
          distance_from_surf = 0.0;
        } else {
          // Otherwise, project towards surface and find the distance
          for (it = 0; it < composite_surfs->size(); it++) {
            m = (*composite_surfs)[it];

            dot = MathExtra::dot3(contact_surfs[m].dr, contact_surfs[m].surf_norm);
            MathExtra::scaleadd3(-dot, contact_surfs[m].surf_norm, contact_surfs[m].dr, residual);
            distance_from_surf = MIN(distance_from_surf, MathExtra::len3(residual));
          }
        }

        // Smooth contributions over 1/4 of radius (arbitrary #)
        for (it = 0; it < composite_surfs->size(); it++) {
          m = (*composite_surfs)[it];
          contact_surfs[m].smooth_ext = MIN(1.0, distance_from_surf / (0.25 * radi));
        }
      }
    }

    processed_contacts->clear();
    convex_contacts->clear();
    concave_contacts->clear();

    for (n = 0; n < n_contact_surfs; n++) {

      j = contact_surfs[n].index;
      if (processed_contacts->find(j) != processed_contacts->end()) continue;


      max_overlap = contact_surfs[n].overlap; // Save here, walking can change
      composite_surfs->clear();
      if (dimension == 2) {
        walk_connections2d(n, composite_surfs, processed_contacts,
                           convex_contacts, contacts_map);
      } else {
        walk_connections3d(n, composite_surfs, processed_contacts,
                           convex_contacts, concave_contacts, contacts_map);
      }

      zero_overlap = rescale_overlaps(max_overlap, composite_surfs);
      if (zero_overlap)
        continue;

      if (concave_contacts->size() != 0)
        process_concave_tris(composite_surfs, concave_contacts);
      if (convex_contacts->size() != 0)
        process_convex_surfs(composite_surfs, convex_contacts);

      max_overlap = 0.0;
      for (it = 0; it < composite_surfs->size(); it++) {
        m = (*composite_surfs)[it];
        max_overlap = MAX(max_overlap, contact_surfs[m].overlap);
      }

      if (max_overlap < EPSILON)
        continue;

      // Calculate geometry of contact
      if (composite_surfs->size() > 1) {

        // Calculate overlap-weighted average normal vector
        MathExtra::zero3(dr);
        for (it = 0; it < composite_surfs->size(); it++) {
          m = (*composite_surfs)[it];
          MathExtra::scaleadd3(contact_surfs[m].overlap, contact_surfs[m].force_norm, dr, dr);
        }

        MathExtra::norm3(dr);
        MathExtra::scale3(radi - max_overlap, dr);
      } else {
        MathExtra::scale3(radi - max_overlap, contact_surfs[n].force_norm, dr);
      }
      //INSERT MORE
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

/* ----------------------------------------------------------------------
   recursively walk through contacting connections and determine side of contact
------------------------------------------------------------------------- */

void PairSurfGranular::prewalk_connections2d(int n, int nsidej, std::unordered_set<int> *processed_contacts, std::map<int, int> *contacts_map)
{
  int j = contact_surfs[n].index;

  processed_contacts->insert(j);
  if (contact_surfs[n].nside == OPPOSITE_SIDE)
    MathExtra::negate3(contact_surfs[n].surf_norm);

  int k, m, nsidek, nconnect;

  for (nconnect = 0; nconnect < (connect2d[j].np1 + connect2d[j].np2); nconnect++) {
    if (nconnect < connect2d[j].np1) {
      k = connect2d[j].neigh_p1[nconnect];
      nsidek = connect2d[j].nside_p1[nconnect];
    } else {
      k = connect2d[j].neigh_p2[nconnect - connect2d[j].np1];
      nsidek = connect2d[j].nside_p2[nconnect - connect2d[j].np1];
    }

    // Skip if not in contact
    if (contacts_map->find(k) == contacts_map->end())
      continue;

    if (processed_contacts->find(k) == processed_contacts->end()) {
      // which side is associated with the initial closest surf
      m = (*contacts_map)[k];
      if (nsidej == OPPOSITE_SIDE)
        nsidek = FLIPSIDE(nsidek);
      contact_surfs[m].nside = nsidek;
      prewalk_connections2d(m, nsidek, processed_contacts, contacts_map);
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairSurfGranular::prewalk_connections3d(int n, int nsidej, std::vector<int> *composite_surfs, std::unordered_set<int> *processed_contacts, std::map<int, int> *contacts_map)
{
  int j = contact_surfs[n].index;

  processed_contacts->insert(j);
  composite_surfs->push_back(n);

  if (contact_surfs[n].nside == OPPOSITE_SIDE)
    MathExtra::negate3(contact_surfs[n].surf_norm);

  int k, m, nsidek, nconnect, nc, ntotal;

  // Loop through edge-connected surfs
  ntotal = connect3d[j].ne1 + connect3d[j].ne2 + connect3d[j].ne3;
  for (nconnect = 0; nconnect < ntotal; nconnect++) {
    if (nconnect < connect3d[j].ne1) {
      nc = nconnect;
      k = connect3d[j].neigh_e1[nc];
      nsidek = connect3d[j].nside_e1[nc];
    } else if (nconnect < connect3d[j].ne1 + connect3d[j].ne2) {
      nc = nconnect - connect3d[j].ne1;
      k = connect3d[j].neigh_e2[nc];
      nsidek = connect3d[j].nside_e2[nc];
    } else {
      nc = nconnect - connect3d[j].ne1 - connect3d[j].ne2;
      k = connect3d[j].neigh_e3[nc];
      nsidek = connect3d[j].nside_e3[nc];
    }

    // Skip if not in contact
    if (contacts_map->find(k) == contacts_map->end())
      continue;

    if (processed_contacts->find(k) == processed_contacts->end()) {
      // which side is associated with the initial closest surf
      m = (*contacts_map)[k];
      if (nsidej == OPPOSITE_SIDE)
        nsidek = FLIPSIDE(nsidek);
      contact_surfs[m].nside = nsidek;
      prewalk_connections3d(m, nsidek, composite_surfs, processed_contacts, contacts_map);
    }
  }

  // Loop through corner-connected surfs
  ntotal = connect3d[j].nc1 + connect3d[j].nc2 + connect3d[j].nc3;
  for (nconnect = 0; nconnect < ntotal; nconnect++) {
    if (nconnect < connect3d[j].nc1) {
      nc = nconnect;
      k = connect3d[j].neigh_c1[nc];
      nsidek = connect3d[j].nside_c1[nc];
    } else if (nconnect < connect3d[j].nc1 + connect3d[j].nc2) {
      nc = nconnect - connect3d[j].nc1;
      k = connect3d[j].neigh_c2[nc];
      nsidek = connect3d[j].nside_c2[nc];
    } else {
      nc = nconnect - connect3d[j].nc1 - connect3d[j].nc2;
      k = connect3d[j].neigh_c3[nc];
      nsidek = connect3d[j].nside_c3[nc];
    }

    // Skip if not in contact
    if (contacts_map->find(k) == contacts_map->end())
      continue;

    if (processed_contacts->find(k) == processed_contacts->end()) {
      // which side is associated with the initial closest surf
      m = (*contacts_map)[k];
      if (nsidej == OPPOSITE_SIDE)
        nsidek = FLIPSIDE(nsidek);
      contact_surfs[m].nside = nsidek;
      prewalk_connections3d(m, nsidek, composite_surfs, processed_contacts, contacts_map);
    }
  }
}

/* ----------------------------------------------------------------------
   recursively walk through flat connections and process any contacts
------------------------------------------------------------------------- */

void PairSurfGranular::walk_connections2d(int n, std::vector<int> *composite_surfs, std::unordered_set<int> *processed_contacts, std::unordered_set<int> *convex_contacts, std::map<int, int> *contacts_map)
{
  int j = contact_surfs[n].index;

  processed_contacts->insert(j);
  composite_surfs->push_back(n);

  int k, m, aflag, which, nconnect, convex_flag, contact_at_joint;
  int jflag = contact_surfs[n].flag;

  for (nconnect = 0; nconnect < (connect2d[j].np1 + connect2d[j].np2); nconnect++) {
    contact_at_joint = 0; // If j's contact is at j-k joint
    if (nconnect < connect2d[j].np1) {
      k = connect2d[j].neigh_p1[nconnect];
      aflag = connect2d[j].aflag_p1[nconnect];
      which = connect2d[j].pwhich_p1[nconnect];
      if (jflag == -1)
        contact_at_joint = 1;
    } else {
      k = connect2d[j].neigh_p2[nconnect - connect2d[j].np1];
      aflag = connect2d[j].aflag_p2[nconnect - connect2d[j].np1];
      which = connect2d[j].pwhich_p2[nconnect - connect2d[j].np1];
      if (jflag == -2)
        contact_at_joint = 2;
    }

    // Skip if not in contact
    if (contacts_map->find(k) == contacts_map->end())
      continue;

    m = (*contacts_map)[k];

    // Different type flat surfs act independently
    if (aflag == FLAT && contact_surfs[n].type == contact_surfs[m].type) {
      // flat, same-type: walk
      if (processed_contacts->find(k) == processed_contacts->end())
        walk_connections2d(m, composite_surfs, processed_contacts, convex_contacts, contacts_map);
    } else {
      convex_flag = 0;
      if ((contact_surfs[n].nside == SAME_SIDE && aflag == CONVEX) ||
          (contact_surfs[n].nside == OPPOSITE_SIDE && aflag == CONCAVE))
        convex_flag = 1;

      if (convex_flag) {
        // convex: process later
        convex_contacts->insert(k);
      } else if (contact_at_joint) {
        // contacting at concave joint: use normal, overriding default use of dr when exposed
        //   unlike 3D, no need to propagate thru corners so process here
        MathExtra::copy3(contact_surfs[n].surf_norm, contact_surfs[n].force_norm);
        contact_surfs[n].norm_def = 1;
      }
    }
  }

  // If norm hasn't been defined
  if (contact_surfs[n].norm_def == -1) {
    if (contact_surfs[n].exposed) {
      // use dr for exposed contacts
      MathExtra::normalize3(contact_surfs[n].dr, contact_surfs[n].force_norm);
      contact_surfs[n].norm_def = 2;
    } else {
      // otherwise default to surface normal
      MathExtra::copy3(contact_surfs[n].surf_norm, contact_surfs[n].force_norm);
      contact_surfs[n].norm_def = 0;
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairSurfGranular::walk_connections3d(int n, std::vector<int> *composite_surfs, std::unordered_set<int> *processed_contacts, std::unordered_set<int> *convex_contacts, std::unordered_set<int> *concave_contacts, std::map<int, int> *contacts_map)
{
  int j = contact_surfs[n].index;

  processed_contacts->insert(j);
  composite_surfs->push_back(n);

  int k, m, aflag, which, nconnect, nc, ntotal, convex_flag, contact_at_joint;
  double r, dot;

  int jflag = contact_surfs[n].flag;

  int needs_correction = 0;
  if (contact_surfs[n].exposed && jflag < -3)
    needs_correction = 1;

  // Whether the composite surface - made up of many flat surfs - is contacted at an exposed point/edge
  int external_to_composite = contact_surfs[n].smooth_ext > EPSILON;

  // Loop through edge-connected surfs
  ntotal = connect3d[j].ne1 + connect3d[j].ne2 + connect3d[j].ne3;
  for (nconnect = 0; nconnect < ntotal; nconnect++) {
    contact_at_joint = 0; // If j's contact is at j-k joint
    if (nconnect < connect3d[j].ne1) {
      // e1 = p1+p2
      nc = nconnect;
      k = connect3d[j].neigh_e1[nc];
      aflag = connect3d[j].aflag_e1[nc];
      which = connect3d[j].ewhich_e1[nc];
      if (jflag == -1 || jflag == -3 || jflag == -4) // do i need corners?
        contact_at_joint = 1;
    } else if (nconnect < connect3d[j].ne1 + connect3d[j].ne2) {
      // e2 = p2+p3
      nc = nconnect - connect3d[j].ne1;
      k = connect3d[j].neigh_e2[nc];
      aflag = connect3d[j].aflag_e2[nc];
      which = connect3d[j].ewhich_e2[nc];
      if (jflag == -2 || jflag == -4 || jflag == -5)
        contact_at_joint = 1;
    } else {
      // e3 = p1+p3
      nc = nconnect - connect3d[j].ne1 - connect3d[j].ne2;
      k = connect3d[j].neigh_e3[nc];
      aflag = connect3d[j].aflag_e3[nc];
      which = connect3d[j].ewhich_e3[nc];
      if (jflag == -3 || jflag == -1 || jflag == -5)
        contact_at_joint = 1;
    }

    // Skip if not in contact
    if (contacts_map->find(k) == contacts_map->end())
      continue;

    m = (*contacts_map)[k];

    // Different type flat surfs act independently
    if (aflag == FLAT && contact_surfs[n].type == contact_surfs[m].type) {
      // flat, same-type: walk
      if (processed_contacts->find(k) == processed_contacts->end())
        walk_connections3d(m, composite_surfs, processed_contacts, convex_contacts, concave_contacts, contacts_map);

      // Adjust dr if j it's an exposed corner
      if (needs_correction && contact_at_joint)
        adjust_exposed_corner_int(j, k, n, m);

    } else {
      convex_flag = 0;
      if ((contact_surfs[n].nside == SAME_SIDE && aflag == CONVEX) ||
          (contact_surfs[n].nside == OPPOSITE_SIDE && aflag == CONCAVE))
        convex_flag = 1;

      if (convex_flag) {
        // convex: process later
        convex_contacts->insert(k);
      } else if (contact_at_joint) {
        // contacting at concave joint: process later
        concave_contacts->insert(k);
      }
    }

    // Adjust dr if j it's an exposed corner and k's exposed for external contacts
    if (needs_correction && external_to_composite && contact_surfs[m].exposed)
      adjust_exposed_corner_ext(j, k, n, m);
  }

  // Loop through corner-connected surfs to find any other flat connections
  ntotal = connect3d[j].nc1 + connect3d[j].nc2 + connect3d[j].nc3;
  for (nconnect = 0; nconnect < ntotal; nconnect++) {
    contact_at_joint = 0; // If j's contact is at j-k joint
    if (nconnect < connect3d[j].nc1) {
      nc = nconnect;
      k = connect3d[j].neigh_c1[nc];
      aflag = connect3d[j].aflag_c1[nc];
      if (jflag == -4)
        contact_at_joint = 1;
    } else if (nconnect < connect3d[j].nc1 + connect3d[j].nc2) {
      nc = nconnect - connect3d[j].nc1;
      k = connect3d[j].neigh_c2[nc];
      aflag = connect3d[j].aflag_c2[nc];
      if (jflag == -5)
        contact_at_joint = 1;
    } else {
      nc = nconnect - connect3d[j].nc1 - connect3d[j].nc2;
      k = connect3d[j].neigh_c3[nc];
      aflag = connect3d[j].aflag_c3[nc];
      if (jflag == -6)
        contact_at_joint = 1;
    }

    // Skip if not in contact
    if (contacts_map->find(k) == contacts_map->end())
      continue;

    m = (*contacts_map)[k];

    if (aflag == FLAT && contact_surfs[n].type == contact_surfs[m].type) {
      if (processed_contacts->find(k) == processed_contacts->end())
        walk_connections3d(m, composite_surfs, processed_contacts, convex_contacts, concave_contacts, contacts_map);

      if (needs_correction && contact_at_joint)
        adjust_exposed_corner_int(j, k, n, m);
    }
    if (needs_correction && external_to_composite && contact_surfs[m].exposed)
      adjust_exposed_corner_ext(j, k, n, m);
  }

  if (needs_correction) {
    // If correction never applied, default to dr
    if (MathExtra::lensq3(contact_surfs[n].cor_int) <= EPSILON)
      MathExtra::copy3(contact_surfs[n].dr, contact_surfs[n].cor_int);
    if (MathExtra::lensq3(contact_surfs[n].cor_ext) <= EPSILON)
      MathExtra::copy3(contact_surfs[n].dr, contact_surfs[n].cor_ext);

    // Smooth between two corrections based on externality of contact
    MathExtra::scaleadd3(contact_surfs[n].smooth_ext, contact_surfs[n].cor_ext, (1.0 - contact_surfs[n].smooth_ext), contact_surfs[n].cor_int, contact_surfs[n].dr);
  }

  // If norm hasn't been defined
  if (contact_surfs[n].norm_def == -1) {
    if (contact_surfs[n].exposed) {
      // use dr for exposed contacts
      MathExtra::normalize3(contact_surfs[n].dr, contact_surfs[n].force_norm);
      contact_surfs[n].norm_def = 2 + needs_correction;
      if (contact_surfs[n].exposed && jflag < -3)
        contact_surfs[n].norm_def = 4 + needs_correction;
    } else {
      // otherwise default to surface normal
      MathExtra::copy3(contact_surfs[n].surf_norm, contact_surfs[n].force_norm);
      contact_surfs[n].norm_def = 0;
    }
  }

  // If contact is on an exposed surface, smooth relative flat contribution
  //   Will rescale overlaps to restore absolute value in compute()
  if (external_to_composite && !contact_surfs[n].exposed)
    contact_surfs[n].overlap *= (1.0 - contact_surfs[n].smooth_ext);
}

/* ----------------------------------------------------------------------
   With an internal contact, correct force normal vector for exposed corners (j),
     relative to their connected neighbor with the highest overlap (k)
     Essential force normal so it doesn't point beyond the two surf norms
------------------------------------------------------------------------- */

void PairSurfGranular::adjust_exposed_corner_int(int j, int k, int n, int m)
{
  // Skip if smoothed to nothing
  if (contact_surfs[n].smooth_ext == 1)
    return;

  // Already adjusted by closer surf
  if (contact_surfs[n].int_overlap >= contact_surfs[m].overlap)
    return;

  double jnorm[3], knorm[3], drnorm[3];
  MathExtra::copy3(contact_surfs[n].surf_norm, jnorm);
  MathExtra::copy3(contact_surfs[m].surf_norm, knorm);
  MathExtra::normalize3(contact_surfs[n].dr, drnorm);

  double dotjr, dotkr;
  dotjr = MathExtra::dot3(jnorm, drnorm);
  dotkr = MathExtra::dot3(knorm, drnorm);

  // Todo: confirm this is smooth around a corner shared by 3+ tris w/
  //   distinct surface normals. However, error no more than flat threshold

  if (contact_surfs[m].flag > -4) {
    // If contacts not all at the corner (i.e. k touches an edge/internal)
    //   -> one surface will dominate
    if (dotjr >= dotkr) {
      // If closer to jnorm, set to jnorm (must be a concave corner)
      MathExtra::copy3(jnorm, contact_surfs[n].cor_int);
    } else {
      // else, set to knorm (must be a convex corner)
      MathExtra::copy3(knorm, contact_surfs[n].cor_int);
    }
  } else {
    // All surfs must touch at the corner (convex corner), just use dr
    MathExtra::copy3(drnorm, contact_surfs[n].cor_int);
  }
  contact_surfs[n].int_overlap = contact_surfs[m].overlap;
}

/* ----------------------------------------------------------------------
   For exposed corner contacts (j) with an external contact, check if connected
     exposed tri (k) contains vector
     if so, adjust dr as necessary to remove any component that lies inside it
------------------------------------------------------------------------- */

void PairSurfGranular::adjust_exposed_corner_ext(int j, int k, int n, int m)
{
  // Skip if smoothed to nothing
  if (contact_surfs[n].smooth_ext == 0)
    return;

  // Already adjusted by closer surf
  if (contact_surfs[n].ext_overlap >= contact_surfs[m].overlap)
    return;

  // Get j's exposed edge vector
  //   Note: if there's two, this will arbitrarily pick one
  //         rare, only if j and k "kiss" at the corner

  int pt = -1;
  int ptj = -1;
  int ptj3;
  if (contact_surfs[n].flag == -4) {
    pt = 0;
    if (connect3d[j].exposed_edge[0]) {
      ptj = 3;
      ptj3 = 6;
    } else if (connect3d[j].exposed_edge[2]) {
      ptj = 6;
      ptj3 = 3;
    }
  } else if (contact_surfs[n].flag == -5) {
    pt = 3;
    if (connect3d[j].exposed_edge[0]) {
      ptj = 0;
      ptj3 = 6;
    } else if (connect3d[j].exposed_edge[1]) {
      ptj = 6;
      ptj3 = 0;
    }
  } else if (contact_surfs[n].flag == -6) {
    pt = 6;
    if (connect3d[j].exposed_edge[1]) {
      ptj = 3;
      ptj3 = 0;
    } else if (connect3d[j].exposed_edge[2]) {
      ptj = 0;
      ptj3 = 3;
    }
  }

  if (ptj == -1) {
    // If a tri that pokes a corner onto perimeter, just remove contribution
    MathExtra::copy3(contact_surfs[n].surf_norm, contact_surfs[n].cor_ext);
    contact_surfs[n].overlap = 0.0;
    return;
  }

  int *tri = atom->tri;
  double jline[3];
  MathExtra::sub3(corners[tri[j]][ptj], corners[tri[j]][pt], jline);
  MathExtra::norm3(jline);

  // If k is an edge, check if its edge-in-contact connects to j's corner
  //   if k is a corner, then find if one of the edges is exposed
  //   if neither, skip

  int ptk = -1;
  int ptk3;
  if (contact_surfs[m].flag == -1) {
    if (EQUAL3(corners[tri[j]][pt], corners[tri[k]][0])) ptk = 3;
    if (EQUAL3(corners[tri[j]][pt], corners[tri[k]][3])) ptk = 0;
    ptk3 = 6;
  } else if (contact_surfs[m].flag == -2) {
    if (EQUAL3(corners[tri[j]][pt], corners[tri[k]][3])) ptk = 6;
    if (EQUAL3(corners[tri[j]][pt], corners[tri[k]][6])) ptk = 3;
    ptk3 = 0;
  } else if (contact_surfs[m].flag == -3) {
    if (EQUAL3(corners[tri[j]][pt], corners[tri[k]][0])) ptk = 6;
    if (EQUAL3(corners[tri[j]][pt], corners[tri[k]][6])) ptk = 0;
    ptk3 = 3;
  } else if (contact_surfs[m].flag == -4) {
    if (connect3d[k].exposed_edge[0]) {
      ptk = 3;
      ptk3 = ;
    } else if (connect3d[k].exposed_edge[2]) {
      ptk = 6;
      ptk3 = ; hmm
    }
  } else if (contact_surfs[m].flag == -5) {
    if (connect3d[k].exposed_edge[0]) ptk = 0;
    if (connect3d[k].exposed_edge[1]) ptk = 6;
  } else if (contact_surfs[m].flag == -6) {
    if (connect3d[k].exposed_edge[1]) ptk = 3;
    if (connect3d[k].exposed_edge[2]) ptk = 0;
  }

  if (ptk == -1)
    return;

  // Calculate k's edge vector

  double kline[3];
  MathExtra::sub3(corners[tri[k]][ptk], corners[tri[j]][pt], kline);
  MathExtra::norm3(kline);

  // Determine if two vectors are concave or convex
  //   first, get normal vector of each edge (perpendicular to surf norm)

  double jnorm[3], knorm[3];
  MathExtra::cross3(corners[tri[j]][9], jline, jnorm);
  MathExtra::norm3(jnorm);
  MathExtra::cross3(corners[tri[k]][9], kline, knorm);
  MathExtra::norm3(knorm);

  // Correct sign to point away from 3rd corner
  double dot, line3[3];
  MathExtra::sub3(points[ptj3].x, points[pt].x, line3);
  dot = MathExtra::dot3(line3, jnorm);
  if (dot > 0) MathExtra::negate3(jnorm);

  if (pt != tris[k].p1 && ptk != tris[k].p1) pt3 = tris[k].p1;
  if (pt != tris[k].p2 && ptk != tris[k].p2) pt3 = tris[k].p2;
  if (pt != tris[k].p3 && ptk != tris[k].p3) pt3 = tris[k].p3;
  MathExtra::sub3(points[pt3].x, points[pt].x, line3);
  dot = MathExtra::dot3(line3, knorm);
  if (dot > 0) MathExtra::negate3(knorm);

  // Convert to 2D problem made up of one of the line's normal vectors and
  //   the vector between the two contacts, rjk
  //   Note: could replace rjk with jline/kline, but find it's generally less smooth

  double rjk[3];
  MathExtra::sub3(contact_surfs[n].contact, contact_surfs[m].contact, rjk);

  // Don't adjust if same contact point, will default cor_ext to dr which
  //   is correct if the sphere is on top of a corner
  if (MathExtra::lensq3(rjk) < EPSILON)
    return;

  // Get normal vector of jnorm and rjk but if aligned, swap to knorm
  double n_plane[3];
  MathExtra::cross3(rjk, jnorm, n_plane);
  if (MathExtra::lensq3(n_plane) < EPSILON)
    MathExtra::cross3(rjk, knorm, n_plane);

  if (MathExtra::lensq3(n_plane) < EPSILON)
    error->one(FLERR, "Bad geometry"); // Can this happen?

  // Project dr into 2D plane
  double drnorm[3], dr_proj[3];
  MathExtra::normalize3(contact_surfs[n].dr, drnorm);
  MathExtra::norm3(n_plane);
  dot = MathExtra::dot3(drnorm, n_plane);
  MathExtra::scaleadd3(-dot, n_plane, drnorm, dr_proj);

  // Calculate nj x nk and dot with the 2d plane vector to identify convex/concave
  double jcrossk[3];
  MathExtra::cross3(jnorm, knorm, jcrossk);
  dot = MathExtra::dot3(n_plane, jcrossk);

  double dr_remove[3];
  if (dot < 0.0 || fabs(dot) < EPSILON) {
    // appears concave (or jnorm = knorm), only preserve component along jnorm
    dot = MathExtra::dot3(dr_proj, jnorm);
    MathExtra::scaleadd3(-dot, jnorm, dr_proj, dr_remove);
  } else {
    // appears convex, skip if both corners
    if (contact_surfs[m].flag < -3)
      return;

    // Project k norm into plane
    double dotk = MathExtra::dot3(n_plane, knorm);
    MathExtra::scaleadd3(-dotk, n_plane, knorm, knorm);
    MathExtra::norm3(knorm);

    double dotjk = MathExtra::dot3(jnorm, knorm);
    double dotjdr = MathExtra::dot3(jnorm, dr_proj);
    double dotkdr = MathExtra::dot3(knorm, dr_proj);

    MathExtra::zero3(dr_remove);
    // Only correct if dr_proj to lies outside of jnorm/knorm sector
    if (dotjk > dotjdr || dotjk > dotkdr) {
      // Correct to whichever is closer (could be outside of both)
      if (dotkdr > dotjdr) {
        dot = MathExtra::dot3(dr_proj, knorm);
        MathExtra::scaleadd3(-dot, knorm, dr_proj, dr_remove);
      } else if (dotjk < dotkdr) {
        dot = MathExtra::dot3(dr_proj, jnorm);
        MathExtra::scaleadd3(-dot, jnorm, dr_proj, dr_remove);
      }
    }
  }

  // get surface normal of j surface in 2d plane
  double jnorm2[3];
  double dot2 = MathExtra::dot3(n_plane, contact_surfs[n].surf_norm);
  MathExtra::scaleadd3(-dot2, n_plane, contact_surfs[n].surf_norm, jnorm2);

  // add back in any component that lies along j's normal
  //   to preserve continuity with edge calculation
  dot2 = MathExtra::dot3(dr_remove, jnorm2);
  MathExtra::scaleadd3(-fabs(dot2), jnorm2, dr_remove, dr_remove);

  MathExtra::sub3(drnorm, dr_remove, contact_surfs[n].cor_ext);

  contact_surfs[n].ext_overlap = contact_surfs[m].overlap;
  return;
}

/* ----------------------------------------------------------------------
   Process overlaps of any tris in a convex geometry
------------------------------------------------------------------------- */

void PairSurfGranular::process_convex_surfs(std::vector<int> *composite_surfs, std::unordered_set<int> *convex_contacts)
{
  int n, j, k, zero_overlap;
  double tmp, dist_convex, min_dist, new_overlap;
  double tmp1[3], tmp2[3];

  // Whether surf & composite contacts are external
  int external = contact_surfs[n].exposed == EXTERNAL && contact_surfs[n].smooth_ext > EPSILON;

  // First, check whether this flat structure is hiding another
  //   If so, more heavily weight surfs abutting convex edge to smooth transition
  //   Preserve max overlap

  // Find minimum distance to convex surf
  double max_overlap = -1;
  if (composite_surfs->size() > 1) {

    min_dist = BIG;
    for (auto it = 0; it < composite_surfs->size(); it++) {
      n = (*composite_surfs)[it];
      j = contact_surfs[n].index;

      max_overlap = MAX(max_overlap, contact_surfs[n].overlap);

      contact_surfs[n].dist = BIG;
      for (auto it2 = convex_contacts->begin(); it2 != convex_contacts->end(); it2++) {
        k = *it2;
        if (dimension == 2) {
          overlap_sphere_line(contact_surfs[n].contact, 0.0,
                              points[lines[k].p1].x, points[lines[k].p2].x,
                              tmp1, tmp2, dist_convex);
        } else {
          overlap_sphere_tri(contact_surfs[n].contact, 0.0,
                             points[tris[k].p1].x, points[tris[k].p2].x,
                             points[tris[k].p3].x, tris[k].norm,
                             tmp1, tmp2, dist_convex);
        }
        contact_surfs[n].dist = MIN(contact_surfs[n].dist, dist_convex);
      }

      contact_surfs[n].dist = sqrt(contact_surfs[n].dist);
      min_dist = MIN(min_dist, contact_surfs[n].dist);
    }

    // De-weight surfs that do not directly touch convex turn
    if (min_dist != BIG) {
      for (auto it = 0; it < composite_surfs->size(); it++) {
        n = (*composite_surfs)[it];

        new_overlap = contact_surfs[n].overlap;
        if (contact_surfs[n].dist != 0)
          new_overlap *= min_dist / contact_surfs[n].dist;

        // If external, preserve part of original overlap
        if (external)
          new_overlap = new_overlap * (1.0 - contact_surfs[n].smooth_ext) + contact_surfs[n].overlap * contact_surfs[n].smooth_ext;
        contact_surfs[n].overlap = new_overlap;
      }

      // Restore absolute max overlap
      zero_overlap = rescale_overlaps(max_overlap, composite_surfs);
    }
  }

  // Secondly, check if this flat structure is being hidden by another
  int hidden = 0;
  for (auto it = 0; it < composite_surfs->size(); it++) {
    n = (*composite_surfs)[it];
    j = contact_surfs[n].index;
    if (convex_contacts->find(j) != convex_contacts->end())
      hidden = 1;
  }

  // If so, remove overlap
  if (hidden) {
    for (auto it = 0; it < composite_surfs->size(); it++) {
      n = (*composite_surfs)[it];

      double weight = 0;
      if (contact_surfs[n].exposed == EXTERNAL) {
        // If external, adjust weight based on proportion of dr that lies outside of the plane of jnorm and j's exposed edge

        int j = contact_surfs[n].index;
        int jflag = contact_surfs[n].flag;

        // If this is a corner-corner connection/contact, will arbitrarily pick external side
        int pt1, pt2;
        pt1 = pt2 = -1;
        if (jflag == -1) {
          pt1 = tris[j].p1;
          pt2 = tris[j].p2;
        } else if (jflag == -2) {
          pt1 = tris[j].p2;
          pt2 = tris[j].p3;
        } else if (jflag == -3) {
          pt1 = tris[j].p1;
          pt2 = tris[j].p3;
        } else if (jflag == -4) {
          pt1 = tris[j].p1;
          if (connect3d[j].exposed_edge[0] == EXTERNAL) pt2 = tris[j].p2;
          if (connect3d[j].exposed_edge[2] == EXTERNAL) pt2 = tris[j].p3;
        } else if (jflag == -5) {
          pt1 = tris[j].p2;
          if (connect3d[j].exposed_edge[0] == EXTERNAL) pt2 = tris[j].p1;
          if (connect3d[j].exposed_edge[1] == EXTERNAL) pt2 = tris[j].p3;
        } else if (jflag == -6) {
          pt1 = tris[j].p3;
          if (connect3d[j].exposed_edge[1] == EXTERNAL) pt2 = tris[j].p2;
          if (connect3d[j].exposed_edge[2] == EXTERNAL) pt2 = tris[j].p1;
        }

        if (pt1 == -1 || pt2 == -1) {
          // If a tri that pokes a corner onto perimeter, just remove contribution
          contact_surfs[n].overlap = 0.0;
          continue;
        }

        double jline[3];
        MathExtra::sub3(points[pt1].x, points[pt2].x, jline);
        MathExtra::norm3(jline);

        double n_plane[3];
        MathExtra::cross3(jline, contact_surfs[n].surf_norm, n_plane);
        MathExtra::norm3(n_plane);

        // Correct sign to point away from 3rd corner
        int pt3;
        double dot, line3[3];
        if (pt1 != tris[j].p1 && pt2 != tris[j].p1) pt3 = tris[j].p1;
        if (pt1 != tris[j].p2 && pt2 != tris[j].p2) pt3 = tris[j].p2;
        if (pt1 != tris[j].p3 && pt2 != tris[j].p3) pt3 = tris[j].p3;
        MathExtra::sub3(points[pt3].x, points[pt1].x, line3);
        dot = MathExtra::dot3(line3, n_plane);
        if (dot > 0) MathExtra::negate3(n_plane);

        double drnorm[3];
        MathExtra::normalize3(contact_surfs[n].dr, drnorm);
        dot = MathExtra::dot3(n_plane, drnorm);
        weight = MAX(0.0, dot);
      }
      contact_surfs[n].overlap *= weight;
      contact_surfs[n].norm_def = 6;
    }
  }
}

/* ----------------------------------------------------------------------
   If a concave tri is found, check all flat tris to find connections
    if internal, use surface normal
    if external, weight with contact direction
------------------------------------------------------------------------- */

void PairSurfGranular::process_concave_tris(std::vector<int> *composite_surfs, std::unordered_set<int> *concave_contacts)
{
  int j, k, n, m, nconnect, nc, ntotal, concave;

  for (auto it = 0; it < composite_surfs->size(); it++) {
    n = (*composite_surfs)[it];
    j = contact_surfs[n].index;
    concave = 0;

    ntotal = connect3d[j].ne1 + connect3d[j].ne2 + connect3d[j].ne3;
    for (nconnect = 0; nconnect < ntotal; nconnect++) {
      if (nconnect < connect3d[j].ne1) {
        nc = nconnect;
        k = connect3d[j].neigh_e1[nc];
      } else if (nconnect < connect3d[j].ne1 + connect3d[j].ne2) {
        nc = nconnect - connect3d[j].ne1;
        k = connect3d[j].neigh_e2[nc];
      } else {
        nc = nconnect - connect3d[j].ne1 - connect3d[j].ne2;
        k = connect3d[j].neigh_e3[nc];
      }

      if (concave_contacts->find(k) != concave_contacts->end()) {
        concave = 1;
        break;
      }
    }

    ntotal = connect3d[j].nc1 + connect3d[j].nc2 + connect3d[j].nc3;
    for (nconnect = 0; nconnect < ntotal; nconnect++) {
      if (nconnect < connect3d[j].nc1) {
        nc = nconnect;
        k = connect3d[j].neigh_c1[nc];
      } else if (nconnect < connect3d[j].nc1 + connect3d[j].nc2) {
        nc = nconnect - connect3d[j].nc1;
        k = connect3d[j].neigh_c2[nc];
      } else {
        nc = nconnect - connect3d[j].nc1 - connect3d[j].nc2;
        k = connect3d[j].neigh_c3[nc];
      }

      if (concave_contacts->find(k) != concave_contacts->end()) {
        concave = 1;
        break;
      }
    }

    if (concave) {
      // Note: dr contains all int/ext adjustments
      MathExtra::norm3(contact_surfs[n].dr);
      MathExtra::scaleadd3(contact_surfs[n].smooth_ext, contact_surfs[n].dr, 1.0 - contact_surfs[n].smooth_ext, contact_surfs[n].surf_norm, contact_surfs[n].force_norm);
      MathExtra::norm3(contact_surfs[n].force_norm);
      contact_surfs[n].norm_def = 1;
    }
  }
}

/* ----------------------------------------------------------------------
   Rescales overlaps in composite_surfs to adjust maximum
------------------------------------------------------------------------- */

int PairSurfGranular::rescale_overlaps(double new_max_overlap, std::vector<int> *composite_surfs)
{
  int n;
  double current_max_overlap = -1;
  for (auto it = 0; it < composite_surfs->size(); it++) {
    n = (*composite_surfs)[it];
    current_max_overlap = MAX(current_max_overlap, contact_surfs[n].overlap);
  }
  if (current_max_overlap < EPSILON) return 1;

  double scale = new_max_overlap / current_max_overlap;
  for (auto it = 0; it < composite_surfs->size(); it++) {
    n = (*composite_surfs)[it];
    contact_surfs[n].overlap *= scale;
  }

  return 0;
}
