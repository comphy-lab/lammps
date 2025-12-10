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

enum{NONE, LINE, TRI};
enum{FLAT,CONCAVE,CONVEX};
enum{SAME_SIDE,OPPOSITE_SIDE};
enum{INTERIOR = 0,EXTERNAL,UNCONNECTED};

static constexpr int DELTACONTACTS = 4;
static constexpr double EPSILON = 1e-12;
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
  nmax_contact_surfs = -1;
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
  int i, j, k, a, n, m, iconnect, jconnect, nconnect, ii, jj;
  int inum, jnum, itype, jtype;
  int isphere, itri, jflag, kflag, external_flag, priority, hidden;
  double xtmp, ytmp, ztmp, radi, delx, dely, delz;
  double rsq, rsq_com, rmag, radsum, max_overlap, dot;
  double factor_lj, mi, mj, meff;
  double norm[3], dr[3], contact[3], ds[3], xc[3], vc[3], omegac[3];
  double *endpt, *corner, *forces, *torquesi, *torquesj, dq;
  double omega0[3] = {0.0, 0.0, 0.0};

  int it, jjtmp, nsidej;
  std::vector<int> *composite_surfs = new std::vector<int>();
  std::unordered_set<int> *processed_contacts = new std::unordered_set<int>();

  int *ilist, *jlist, *numneigh, **firstneigh;
  int *touch, **firsttouch;
  double *history, *allhistory, **firsthistory;

  bool touchflag = false;
  const bool history_update = update->setupflag == 0;

  class GranularModel* model;
  for (n = 0; n < nmodels; n++) {
    model = models_list[n];
    model->history_update = history_update;
    model->radj = 0.0;
  }

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
    atom2connect = fsl->atom2connect;
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
    itype = type[i];
    iconnect = atom2connect[i];

    if (use_history) {
      touch = firsttouch[i];
      allhistory = firsthistory[i];
    }

    jlist = firstneigh[i];
    jnum = numneigh[i];

    n_contact_surfs = 0;
    for (jj = 0; jj < jnum; jj++) {
      j = jlist[jj];

      // Do we need special bonds?
      // factor_lj = special_lj[sbmask(j)];
      // j &= NEIGHMASK;
      // if (factor_lj == 0) continue;

      // sanity check that neighbor list is built correctly

      if ((style == LINE) && (line[i] >= 0 || line[j] < 0))
        error->one(FLERR, "Pair surf/granular interaction is invalid");

      if ((style == TRI) && (tri[i] >= 0 || tri[j] < 0))
        error->one(FLERR, "Pair surf/granular interaction is invalid");

      delx = xtmp - x[j][0];
      dely = ytmp - x[j][1];
      delz = ztmp - x[j][2];
      rsq_com = delx * delx + dely * dely + delz * delz;

      // skip contact check if particle/surf are too far apart

      radsum = radi + radius[j];

      if (rsq_com > radsum * radsum) {
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
          overlap_sphere_line(x[i], radi, &endpt[0], &endpt[3], contact, dr, rsq);

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
                                      "surf/granular:contact_surfs");

      }

      // Find out if contact is on an external edge/corner
      jconnect = atom2connect[j];
      external_flag = INTERIOR;
      if (style == LINE) {
        MathExtra::copy3(&endpts[line[j]][6], norm);
        dot = MathExtra::dot3(norm, dr);
        if (jflag == -1) external_flag = connect2d[jconnect].external_pt[0];
        if (jflag == -2) external_flag = connect2d[jconnect].external_pt[1];
      } else {
        MathExtra::copy3(&corners[tri[j]][9], norm);
        dot = MathExtra::dot3(norm, dr);
        if (jflag == -1) external_flag = connect3d[jconnect].external_edge[0];
        if (jflag == -2) external_flag = connect3d[jconnect].external_edge[1];
        if (jflag == -3) external_flag = connect3d[jconnect].external_edge[2];
        if (jflag == -4) external_flag = connect3d[jconnect].external_pt[0];
        if (jflag == -5) external_flag = connect3d[jconnect].external_pt[1];
        if (jflag == -6) external_flag = connect3d[jconnect].external_pt[2];
      }

      // Store which side is in contact relative to normal vector
      if (dot >= 0) nsidej = SAME_SIDE;
      else nsidej = OPPOSITE_SIDE;

      // Currently does not handle unconnected edges/corners in 3D
      if (dimension == 3) {
        if (external_flag == UNCONNECTED) {
          error->warning(FLERR, "Contact detected with an unconnected tri edge/corner");
          external_flag = EXTERNAL;
        }
      }

      rmag = sqrt(rsq);
      MathExtra::scale3(1.0 / rmag, dr, dr);

      priority = 2;
      if (jflag != 1) priority = 1;
      if (dimension == 3 && jflag < -3) priority = 0;

      contact_surfs[n_contact_surfs].index = j;
      contact_surfs[n_contact_surfs].neigh_index = jj;
      contact_surfs[n_contact_surfs].type = type[j];
      contact_surfs[n_contact_surfs].flag = jflag;
      contact_surfs[n_contact_surfs].external = external_flag;
      contact_surfs[n_contact_surfs].nside = nsidej;
      contact_surfs[n_contact_surfs].overlap = radi - rmag;
      contact_surfs[n_contact_surfs].rmag = rmag;
      contact_surfs[n_contact_surfs].rank_ext = MAXSMALLINT;
      contact_surfs[n_contact_surfs].copy_index_ext = -1;
      contact_surfs[n_contact_surfs].flat_ext = 0;
      contact_surfs[n_contact_surfs].weight_contribution = 1.0;
      contact_surfs[n_contact_surfs].hidden = 0;
      contact_surfs[n_contact_surfs].convex_preceding_contact = -1;
      contact_surfs[n_contact_surfs].rsq_com = rsq_com;
      contact_surfs[n_contact_surfs].priority = priority;

      MathExtra::zero3(contact_surfs[n_contact_surfs].dr_force);
      MathExtra::copy3(norm, contact_surfs[n_contact_surfs].surf_norm);
      MathExtra::copy3(contact, contact_surfs[n_contact_surfs].contact);
      MathExtra::copy3(dr, contact_surfs[n_contact_surfs].dr);
      MathExtra::copy3(dr, contact_surfs[n_contact_surfs].dr_ext);

      n_contact_surfs += 1;
    }

    if (n_contact_surfs == 0)
      continue;

    // Sort contacts by overlap and create a map
    std::sort(contact_surfs, contact_surfs + n_contact_surfs, [](FixSurface::ContactSurf a, FixSurface::ContactSurf b) {
        if (a.overlap > (b.overlap + EPSILON)) return 1; // 1st compare overlaps within epsilon
        if (b.overlap > (a.overlap + EPSILON)) return 0;
        if (a.priority > b.priority) return 1; // 2nd, prioritize interior > edge > corner
        if (b.priority > a.priority) return 0;
        double dota = fabs(MathExtra::dot3(a.surf_norm, a.dr)); // sign may not yet be set
        double dotb = fabs(MathExtra::dot3(b.surf_norm, b.dr));
        if (dota > (dotb + EPSILON)) return 1; // 3rd, prioritize which one aligns best
        if (dotb > (dota + EPSILON)) return 0;
        if (a.rsq_com < (b.rsq_com - EPSILON)) return 1; // 4th, prioritize closer CoM
        if (b.rsq_com < (a.rsq_com - EPSILON)) return 0;
        if (a.index < b.index) return 1;
        else return 0;
      });

    contacts_map.clear();
    for (n = 0; n < n_contact_surfs; n++)
      contacts_map[contact_surfs[n].index] = n;

    // Initial walk to assign consistent sides of surfaces
    //   Not guaranteed to work for v. complex geometries (e.g. Mobius)

    processed_contacts->clear();
    if (dimension == 2) prewalk_connections2d();
    else prewalk_connections3d();

    // Given corrected surface norms, resort contacts
    std::sort(contact_surfs, contact_surfs + n_contact_surfs, [](FixSurface::ContactSurf a, FixSurface::ContactSurf b) {
        if (a.overlap > (b.overlap + EPSILON)) return 1; // 1st compare overlaps within epsilon
        if (b.overlap > (a.overlap + EPSILON)) return 0;
        if (a.priority > b.priority) return 1; // 2nd, prioritize interior > edge > corner
        if (b.priority > a.priority) return 0;
        double dota = MathExtra::dot3(a.surf_norm, a.dr);
        double dotb = MathExtra::dot3(b.surf_norm, b.dr);
        if (dota > (dotb + EPSILON)) return 1; // 3rd, prioritize which one aligns best
        if (dotb > (dota + EPSILON)) return 0;
        if (a.rsq_com < (b.rsq_com - EPSILON)) return 1; // 4th, prioritize closer CoM
        if (b.rsq_com < (a.rsq_com - EPSILON)) return 0;
        if (a.index < b.index) return 1;
        else return 0;
      });

    for (n = 0; n < n_contact_surfs; n++)
      contacts_map[contact_surfs[n].index] = n;

    processed_contacts->clear();
    for (n = 0; n < n_contact_surfs; n++) {

      j = contact_surfs[n].index;
      if (processed_contacts->find(j) != processed_contacts->end()) continue;

      composite_surfs->clear();
      if (dimension == 2) {
        walk_connections2d(n, composite_surfs, processed_contacts);
        calculate_2d_forces(composite_surfs);
      } else {
        walk_connections3d(n, composite_surfs, processed_contacts);
        calculate_3d_forces(composite_surfs);
      }

      max_overlap = -BIG;
      hidden = 0;
      for (auto it = 0; it < composite_surfs->size(); it++) {
        m = (*composite_surfs)[it];
        max_overlap = MAX(max_overlap, contact_surfs[m].overlap);
        if (contact_surfs[m].hidden) {
          hidden = 1;
          break;
        }
      }

      if (hidden || (max_overlap < EPSILON))
        continue;

      // Calculate geometry of contact
      if (composite_surfs->size() > 1) {

        // Calculate overlap-weighted average normal vector
        MathExtra::zero3(dr);
        for (it = 0; it < composite_surfs->size(); it++) {
          m = (*composite_surfs)[it];
          MathExtra::scaleadd3(contact_surfs[m].overlap * contact_surfs[m].weight_contribution, contact_surfs[m].dr_force, dr, dr);
        }

        MathExtra::norm3(dr);
        MathExtra::scale3(radi - max_overlap, dr);

      } else {
        MathExtra::scale3(radi - max_overlap, contact_surfs[n].dr_force, dr);
      }

      for (a = 0; a < 3; a++)
        xc[a] = x[i][a] - dr[a];

      MathExtra::zero3(vc);
      MathExtra::zero3(omegac);

      jtype = type[j];
      model = models_list[types_indices[itype][jtype]];
      model->xi = x[i];
      model->radi = radi;
      model->vi = v[i];
      model->omegai = omega[i];
      if (heat_flag) {
        model->Ti = temperature[i];
        model->Tj = temperature[j];
      }

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
      model->meff = meff;

      // Correct velocity at contact point, extending from closest surf j
      ds[0] = xc[0] - x[j][0];
      ds[1] = xc[1] - x[j][1];
      ds[2] = xc[2] - x[j][2];
      vc[0] = v[j][0] + (omega[j][1] * ds[2] - omega[j][2] * ds[1]);
      vc[1] = v[j][1] + (omega[j][2] * ds[0] - omega[j][0] * ds[2]);
      vc[2] = v[j][2] + (omega[j][0] * ds[1] - omega[j][1] * ds[0]);

      model->xj = xc;
      model->vj = vc;
      model->omegaj = omegac; // Ask Dan

      if (use_history) {
        jj = contact_surfs[n].neigh_index;
        model->touch = touch[jj];
      }

      // guaranteed in contact, but need to calculate intermediate variables
      touchflag = model->check_contact();

      if (use_history) {
        // Check if another flat contact has a stored history
        if (touch[jj] != 1) {
          for (it = 0; it < composite_surfs->size(); it++) {
            m = (*composite_surfs)[it];
            jjtmp = contact_surfs[m].neigh_index;
            if (touch[jjtmp] == 1);
              jj = jjtmp;
          }
        }

        touch[jj] = 1;
        history = &allhistory[size_history * jj];
        model->history = history;
      }

      model->calculate_forces();

      // Sychronize history across flat contacts
      //   can be arbitrary if not all connected flat surfaces are mutually flat
      //   e.g. a hair pin turn where surfs on either end of the 'U' are not flat
      if (use_history) {
        for (it = 0; it < composite_surfs->size(); it++) {
          m = (*composite_surfs)[it];
          jjtmp = contact_surfs[m].neigh_index;
          if (jj != jjtmp) {
            touch[jjtmp] = 1;
            for (k = 0; k < size_history; k++)
              allhistory[size_history * jjtmp + k] = history[k];
          }
        }
      }

      forces = model->forces;
      torquesi = model->torquesi;

      add3(f[i], forces, f[i]);
      add3(torque[i], torquesi, torque[i]);

      if (force->newton_pair || j < nlocal) {
        sub3(f[j], forces, f[j]);
        torquesj = model->torquesj;
        add3(torque[j], torquesj, torque[j]);
      }

      if (heat_flag) {
        heatflow[i] += model->dq;
        if (force->newton_pair || j < nlocal) heatflow[j] -= dq;
      }

      if (evflag) {
        ev_tally_xyz(i,j,nlocal,force->newton_pair,
          0.0,0.0,forces[0],forces[1],forces[2],model->dx[0],model->dx[1],model->dx[2]);
      }
    }
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
    error->all(FLERR, "Pair style surf/granular requires newton pair on");
  if (comm->ghost_velocity == 0)
    error->all(FLERR, "Pair surf/granular requires ghost atoms store velocity");

  if (heat_flag) {
    if (!atom->temperature_flag)
      error->all(FLERR, "Heat conduction in pair surf/granularular requires atom style with temperature property");
    if (!atom->heatflow_flag)
      error->all(FLERR, "Heat conduction in pair surf/granularular requires atom style with heatflow property");
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
    if (!fix_history) error->all(FLERR, "Could not find pair fix neigh history ID");
  }

  // set ptr to FixSurfaceLocal for surf connectivity info

  fsl = nullptr;
  for (int m = 0; m < modify->nfix; m++) {
    if (strcmp(modify->fix[m]->style, "surface/local") == 0) {
      if (fsl)
        error->all(FLERR, "Pair surf/granular requires single fix surface/local");
      fsl = (FixSurfaceLocal *) modify->fix[m];
    }
  }
  if (!fsl) error->all(FLERR, "Pair surf/granular requires a fix surface/local");

  // surfmoveflag = 1 if surfs may move at every step
  // yes if fix move exists and its group includes lines
  // NOTE: are there other conditions, like fix deform or fix npt?

  surfmoveflag = 0;
  for (int m = 0; m < modify->nfix; m++) {
    if (strcmp(modify->fix[m]->style, "move") == 0) {
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
      MPI_Allreduce(&flag, &any, 1, MPI_INT, MPI_SUM, world);
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
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);
    else if ((style == TRI) && tri[i] >= 0)
      onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);
    else {
      if (mask[i] & freeze_group_bit)
        onerad_frozen[type[i]] = MAX(onerad_frozen[type[i]], radius[i]);
      else
        onerad_dynamic[type[i]] = MAX(onerad_dynamic[type[i]], radius[i]);
    }
  }

  MPI_Allreduce(&onerad_dynamic[1], &maxrad_dynamic[1], atom->ntypes,
                MPI_DOUBLE, MPI_MAX, world);
  MPI_Allreduce(&onerad_frozen[1], &maxrad_frozen[1], atom->ntypes,
                MPI_DOUBLE, MPI_MAX, world);
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
    memory->create(endpts, emax, 9, "surf/granular:endpts");
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
    dx = 0.5 * length * cos(theta);
    dy = 0.5 * length * sin(theta);
    endpt = endpts[m];
    endpt[0] = x[i][0] - dx;
    endpt[1] = x[i][1] - dy;
    endpt[2] = 0.0;
    endpt[3] = x[i][0] + dx;
    endpt[4] = x[i][1] + dy;
    endpt[5] = 0.0;

    MathExtra::sub3(&endpt[3], &endpt[0], p12);
    MathExtra::cross3(zunit, p12, norm);
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
  double ex[3], ey[3], ez[3], p[3][3];
  double *corner;

  // realloc corners array if necssary

  if (fsl->nmax_connect > cmax) {
    memory->destroy(corners);
    cmax = fsl->nmax_connect;
    memory->create(corners, cmax, 12, "surf/granular:corners");
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
    MathExtra::quat_to_mat(bonus[m].quat, p);
    MathExtra::matvec(p, bonus[m].c1, &corner[0]);
    MathExtra::add3(x[i], &corner[0], &corner[0]);
    MathExtra::matvec(p, bonus[m].c2, &corner[3]);
    MathExtra::add3(x[i], &corner[3], &corner[3]);
    MathExtra::matvec(p, bonus[m].c3, &corner[6]);
    MathExtra::add3(x[i], &corner[6], &corner[6]);
    corners2norm(corner, &corner[9]);

    // omega from angmom of tri particles

    if (angmom[i][0] == 0.0 && angmom[i][1] == 0.0 && angmom[i][2] == 0.0) {
      omega[i][0] = omega[i][1] = omega[i][2] = 0.0;
      continue;
    }
    MathExtra::q_to_exyz(bonus[m].quat, ex, ey, ez);
    MathExtra::angmom_to_omega(angmom[i], ex, ey, ez,
                               bonus[m].inertia, omega[i]);
  }
}

/* ----------------------------------------------------------------------
   compute norm of a triangle based on its 3 corner pts
------------------------------------------------------------------------- */

void PairSurfGranular::corners2norm(double *corners, double *norm)
{
  double p12[3],p13[3];

  MathExtra::sub3(&corners[3], &corners[0], p12);
  MathExtra::sub3(&corners[6], &corners[0], p13);
  MathExtra::cross3(p12, p13, norm);
  MathExtra::norm3(norm);
}

/* ----------------------------------------------------------------------
   recursively walk through contacting connections and determine side of contact
------------------------------------------------------------------------- */

void PairSurfGranular::prewalk_connections2d()
{
  std::map<int, int> to_walk;
  std::unordered_set<int> walked;

  int j = contact_surfs[0].index;
  to_walk[j] = contact_surfs[0].nside;

  tagint ktag;
  int k, n, m, jconnect, nsidej, nsidek, nconnect, nc;
  std::tuple<int, int> element;
  while (!to_walk.empty()) {
    auto it = to_walk.begin();
    element = *it;
    j = it->first;
    jconnect = atom2connect[j];
    nsidej = it->second;
    to_walk.erase(it);
    walked.insert(j);

    n = contacts_map[j];

    if (nsidej == OPPOSITE_SIDE)
      MathExtra::negate3(contact_surfs[n].surf_norm);

    for (nconnect = 0; nconnect < (connect2d[jconnect].np1 + connect2d[jconnect].np2); nconnect++) {
      if (nconnect < connect2d[jconnect].np1) {
        ktag = connect2d[jconnect].neigh_p1[nconnect];
        nsidek = connect2d[jconnect].nside_p1[nconnect];
      } else {
        ktag = connect2d[jconnect].neigh_p2[nconnect - connect2d[jconnect].np1];
        nsidek = connect2d[jconnect].nside_p2[nconnect - connect2d[jconnect].np1];
      }
      k = atom->map(ktag);
      if (k == -1)
        error->one(FLERR, "Bad atom index");

      // Skip if not in contact
      if (contacts_map.find(k) == contacts_map.end())
        continue;

      if (walked.find(k) == walked.end() && to_walk.find(k) == to_walk.end()) {
        // which side is associated with the initial closest surf
        m = contacts_map[k];
        if (nsidej == OPPOSITE_SIDE)
          nsidek = FLIPSIDE(nsidek);

        contact_surfs[m].nside = nsidek;
        to_walk[k] = nsidek;
      }
    }

    // Check if there is another disconnected surf
    if (to_walk.empty()) {
      for (nc = 0; nc < n_contact_surfs; nc++) {
        j = contact_surfs[nc].index;
        if (walked.find(j) == walked.end())
          to_walk[j] = contact_surfs[nc].nside;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairSurfGranular::prewalk_connections3d()
{
  std::map<int, int> to_walk;
  std::unordered_set<int> walked;

  int j = contact_surfs[0].index;
  to_walk[j] = contact_surfs[0].nside;

  tagint ktag;
  int k, n, m, jconnect, nsidej, nsidek, nconnect, nc, ntotal;
  std::tuple<int, int> element;
  while (!to_walk.empty()) {
    auto it = to_walk.begin();
    element = *it;
    j = it->first;
    jconnect = atom2connect[j];
    nsidej = it->second;
    to_walk.erase(it);
    walked.insert(j);

    n = contacts_map[j];

    if (nsidej == OPPOSITE_SIDE)
      MathExtra::negate3(contact_surfs[n].surf_norm);

    // Loop through edge-connected surfs
    ntotal = connect3d[jconnect].ne1 + connect3d[jconnect].ne2 + connect3d[jconnect].ne3;
    for (nconnect = 0; nconnect < ntotal; nconnect++) {
      if (nconnect < connect3d[jconnect].ne1) {
        nc = nconnect;
        ktag = connect3d[jconnect].neigh_e1[nc];
        nsidek = connect3d[jconnect].nside_e1[nc];
      } else if (nconnect < connect3d[jconnect].ne1 + connect3d[jconnect].ne2) {
        nc = nconnect - connect3d[jconnect].ne1;
        ktag = connect3d[jconnect].neigh_e2[nc];
        nsidek = connect3d[jconnect].nside_e2[nc];
      } else {
        nc = nconnect - connect3d[jconnect].ne1 - connect3d[jconnect].ne2;
        ktag = connect3d[jconnect].neigh_e3[nc];
        nsidek = connect3d[jconnect].nside_e3[nc];
      }
      k = atom->map(ktag);
      if (k == -1)
        error->one(FLERR, "Bad atom index");

      // Skip if not in contact
      if (contacts_map.find(k) == contacts_map.end())
        continue;

      if (walked.find(k) == walked.end() && to_walk.find(k) == to_walk.end()) {
        // which side is associated with the initial closest surf
        m = contacts_map[k];
        if (nsidej == OPPOSITE_SIDE)
          nsidek = FLIPSIDE(nsidek);
        contact_surfs[m].nside = nsidek;
        to_walk[k] = nsidek;
      }
    }

    // Loop through corner-connected surfs
    ntotal = connect3d[jconnect].nc1 + connect3d[jconnect].nc2 + connect3d[jconnect].nc3;
    for (nconnect = 0; nconnect < ntotal; nconnect++) {
      if (nconnect < connect3d[jconnect].nc1) {
        nc = nconnect;
        ktag = connect3d[jconnect].neigh_c1[nc];
        nsidek = connect3d[jconnect].nside_c1[nc];
      } else if (nconnect < connect3d[jconnect].nc1 + connect3d[jconnect].nc2) {
        nc = nconnect - connect3d[jconnect].nc1;
        ktag = connect3d[jconnect].neigh_c2[nc];
        nsidek = connect3d[jconnect].nside_c2[nc];
      } else {
        nc = nconnect - connect3d[jconnect].nc1 - connect3d[jconnect].nc2;
        ktag = connect3d[jconnect].neigh_c3[nc];
        nsidek = connect3d[jconnect].nside_c3[nc];
      }

      k = atom->map(ktag);
      if (k == -1)
        error->one(FLERR, "Bad atom index");

      // Skip if not in contact
      if (contacts_map.find(k) == contacts_map.end())
        continue;

      if (walked.find(k) == walked.end() && to_walk.find(k) == to_walk.end()) {
        // which side is associated with the initial closest surf
        m = contacts_map[k];
        if (nsidej == OPPOSITE_SIDE)
          nsidek = FLIPSIDE(nsidek);
        contact_surfs[m].nside = nsidek;
        to_walk[k] = nsidek;
      }
    }

    // Check if there is another disconnected surf
    if (to_walk.empty()) {
      for (nc = 0; nc < n_contact_surfs; nc++) {
        j = contact_surfs[nc].index;
        if (walked.find(j) == walked.end())
          to_walk[j] = contact_surfs[nc].nside;
      }
    }
  }
}

/* ----------------------------------------------------------------------
   recursively walk through flat connections and process any contacts
------------------------------------------------------------------------- */

void PairSurfGranular::walk_connections2d(int n, std::vector<int> *composite_surfs, std::unordered_set<int> *processed_contacts)
{
  int j = contact_surfs[n].index;
  int jconnect = atom2connect[j];

  processed_contacts->insert(j);
  composite_surfs->push_back(n);

  tagint ktag;
  int k, m, aflag, which, nconnect, nc, convex_flag, contact_at_joint;
  int jflag = contact_surfs[n].flag;

  // Whether surf contact is at an external pt
  int needs_correction = 0;
  if (contact_surfs[n].external && jflag < 0) // (I think jflag < 0 is a given)
    needs_correction = 1;

  for (nconnect = 0; nconnect < (connect2d[jconnect].np1 + connect2d[jconnect].np2); nconnect++) {
    contact_at_joint = 0; // If j's contact is at j-k joint
    if (nconnect < connect2d[jconnect].np1) {
      nc = nconnect;
      ktag = connect2d[jconnect].neigh_p1[nc];
      aflag = connect2d[jconnect].aflag_p1[nc];
      which = connect2d[jconnect].pwhich_p1[nc];
      if (jflag == -1)
        contact_at_joint = 1;
    } else {
      nc = nconnect - connect2d[jconnect].np1;
      ktag = connect2d[jconnect].neigh_p2[nc];
      aflag = connect2d[jconnect].aflag_p2[nc];
      which = connect2d[jconnect].pwhich_p2[nc];
      if (jflag == -2)
        contact_at_joint = 2;
    }
    k = atom->map(ktag);
    if (k == -1)
      error->one(FLERR, "Bad atom index");

    // Skip if not in contact
    if (contacts_map.find(k) == contacts_map.end())
      continue;

    m = contacts_map[k];

    if (aflag == FLAT) {
      // flat, same-type: walk
      if (contact_surfs[n].type == contact_surfs[m].type &&
          processed_contacts->find(k) == processed_contacts->end())
        walk_connections2d(m, composite_surfs, processed_contacts);

      if (needs_correction && contact_at_joint)
        adjust_external_pt_flat_2d(j, k, n, m);
    } else {
      convex_flag = 0;
      if ((contact_surfs[n].nside == SAME_SIDE && aflag == CONVEX) ||
          (contact_surfs[n].nside == OPPOSITE_SIDE && aflag == CONCAVE))
        convex_flag = 1;

      if (convex_flag) {
        contact_surfs[m].convex_preceding_contact = j;
      } else if (contact_at_joint) {
        // Concave, so just use surface normal
        MathExtra::copy3(contact_surfs[n].surf_norm, contact_surfs[n].dr_ext);
        contact_surfs[n].rank_ext = -2; // correction won't be overwritten
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairSurfGranular::walk_connections3d(int n, std::vector<int> *composite_surfs, std::unordered_set<int> *processed_contacts)
{
  int j = contact_surfs[n].index;
  int jconnect = atom2connect[j];

  processed_contacts->insert(j);
  composite_surfs->push_back(n);

  tagint ktag;
  int k, m, aflag, which, nconnect, nc, ntotal, convex_flag, contact_at_joint;
  int jflag = contact_surfs[n].flag;

  // Whether surf contact is at an external pt/edge
  int needs_pt_correction = 0;
  int needs_edge_correction = 0;
  if (contact_surfs[n].external) {
    if (jflag < -3) {
      needs_pt_correction = 1;
    } else if (jflag < 0) {
      needs_edge_correction = 1;
    }
  }

  // Loop through edge-connected surfs
  ntotal = connect3d[jconnect].ne1 + connect3d[jconnect].ne2 + connect3d[jconnect].ne3;
  for (nconnect = 0; nconnect < ntotal; nconnect++) {
    contact_at_joint = 0; // If j's contact is at j-k joint
    if (nconnect < connect3d[jconnect].ne1) {
      // e1 = p1+p2
      nc = nconnect;
      ktag = connect3d[jconnect].neigh_e1[nc];
      aflag = connect3d[jconnect].aflag_e1[nc];
      which = connect3d[jconnect].ewhich_e1[nc];
      if (jflag == -1 || jflag == -4 || jflag == -5)
        contact_at_joint = 1;
    } else if (nconnect < connect3d[jconnect].ne1 + connect3d[jconnect].ne2) {
      // e2 = p2+p3
      nc = nconnect - connect3d[jconnect].ne1;
      ktag = connect3d[jconnect].neigh_e2[nc];
      aflag = connect3d[jconnect].aflag_e2[nc];
      which = connect3d[jconnect].ewhich_e2[nc];
      if (jflag == -2 || jflag == -5 || jflag == -6)
        contact_at_joint = 1;
    } else {
      // e3 = p1+p3
      nc = nconnect - connect3d[jconnect].ne1 - connect3d[jconnect].ne2;
      ktag = connect3d[jconnect].neigh_e3[nc];
      aflag = connect3d[jconnect].aflag_e3[nc];
      which = connect3d[jconnect].ewhich_e3[nc];
      if (jflag == -3 || jflag == -4 || jflag == -6)
        contact_at_joint = 1;
    }

    k = atom->map(ktag);
    if (k == -1)
      error->one(FLERR, "Bad atom index");

    // Skip if not in contact
    if (contacts_map.find(k) == contacts_map.end())
      continue;

    m = contacts_map[k];

    if (aflag == FLAT) {
      // flat, same-type: walk
      if (contact_surfs[n].type == contact_surfs[m].type &&
          processed_contacts->find(k) == processed_contacts->end())
        walk_connections3d(m, composite_surfs, processed_contacts);

      if (needs_edge_correction && contact_at_joint)
        adjust_external_edge_flat_3d(j, k, n, m);

      if (needs_pt_correction && contact_at_joint)
        adjust_external_pt_flat_3d(j, k, n, m);
    } else {
      convex_flag = 0;
      if ((contact_surfs[n].nside == SAME_SIDE && aflag == CONVEX) ||
          (contact_surfs[n].nside == OPPOSITE_SIDE && aflag == CONCAVE))
        convex_flag = 1;

      if (convex_flag) {
        contact_surfs[m].convex_preceding_contact = j;
      } else if (needs_edge_correction && contact_at_joint) {
        // Concave edge, so just use surface normal
        MathExtra::copy3(contact_surfs[n].surf_norm, contact_surfs[n].dr_ext);
        contact_surfs[n].rank_ext = -2; // correction won't be overwritten
      }

      if (needs_pt_correction && contact_at_joint)
        adjust_external_pt_nonflat_3d(j, k, n, m);
    }
  }

  // Loop through corner-connected surfs to find any other flat connections
  ntotal = connect3d[jconnect].nc1 + connect3d[jconnect].nc2 + connect3d[jconnect].nc3;
  for (nconnect = 0; nconnect < ntotal; nconnect++) {
    contact_at_joint = 0;
    if (nconnect < connect3d[jconnect].nc1) {
      nc = nconnect;
      ktag = connect3d[jconnect].neigh_c1[nc];
      aflag = connect3d[jconnect].aflag_c1[nc];
      if (jflag == -4)
        contact_at_joint = 1;
    } else if (nconnect < connect3d[jconnect].nc1 + connect3d[jconnect].nc2) {
      nc = nconnect - connect3d[jconnect].nc1;
      ktag = connect3d[jconnect].neigh_c2[nc];
      aflag = connect3d[jconnect].aflag_c2[nc];
      if (jflag == -5)
        contact_at_joint = 1;
    } else {
      nc = nconnect - connect3d[jconnect].nc1 - connect3d[jconnect].nc2;
      ktag = connect3d[jconnect].neigh_c3[nc];
      aflag = connect3d[jconnect].aflag_c3[nc];
      if (jflag == -6)
        contact_at_joint = 1;
    }

    k = atom->map(ktag);
    if (k == -1)
      error->one(FLERR, "Bad atom index");

    // Skip if not in contact
    if (contacts_map.find(k) == contacts_map.end())
      continue;

    m = contacts_map[k];

    if (aflag == FLAT) {
      if (contact_surfs[n].type == contact_surfs[m].type &&
          processed_contacts->find(k) == processed_contacts->end())
        walk_connections3d(m, composite_surfs, processed_contacts);

      // Edge connections can only correct flat tris
      //   too hard to correct non-flat kissing tris
      if (needs_pt_correction && contact_at_joint)
        adjust_external_pt_flat_3d(j, k, n, m);
    }
  }
}

/* ----------------------------------------------------------------------
   Calculate forces
------------------------------------------------------------------------- */

void PairSurfGranular::calculate_2d_forces(std::vector<int> *composite_surfs)
{
  int n, m, j, k;
  double dot, residual[3];

  // Calculate properties of composite surface

  n = (*composite_surfs)[0];
  double max_overlap = contact_surfs[n].overlap;
  double max_overlap_ext = -BIG;
  double contact_at_max[3];
  MathExtra::copy3(contact_surfs[n].contact, contact_at_max);

  // Check if composite is hidden (convex) and calc min distance to any external points

  int hidden = 0;
  for (auto it = 0; it < composite_surfs->size(); it++) {
    n = (*composite_surfs)[it];
    j = contact_surfs[n].index;

    if (contact_surfs[n].convex_preceding_contact != -1) {
      // Hide composite surface, unless it's self-hidden (e.g. hairpin turns)
      hidden = 1;
      for (auto it2 = 0; it2 < composite_surfs->size(); it2++) {
        m = (*composite_surfs)[it2];
        k = contact_surfs[m].index;
        if (contact_surfs[n].convex_preceding_contact == k) {
          hidden = 0;
          break;
        }
      }

      if (hidden == 1)
        break;
    }

    if (contact_surfs[n].external)
      max_overlap_ext = MAX(max_overlap_ext, contact_surfs[n].overlap);
  }

  if (hidden) {
    for (auto it = 0; it < composite_surfs->size(); it++) {
      n = (*composite_surfs)[it];
      contact_surfs[n].hidden = 1;
    }
    return;
  }

  // Smooth int/ext based on arbitrary ratio of overlaps
  double w_int = 1.0;
  if (max_overlap_ext != -BIG) {
    w_int = max_overlap_ext / max_overlap;
    w_int = 1.0 - w_int * w_int;
  }
  double w_ext = 1.0 - w_int;

  for (auto it = 0; it < composite_surfs->size(); it++) {
    n = (*composite_surfs)[it];
    if (contact_surfs[n].external) {
      // Interpolate between surface normal and dr_ext using smoothing
      MathExtra::scaleadd3(w_ext, contact_surfs[n].dr_ext, w_int, contact_surfs[n].surf_norm, contact_surfs[n].dr_force);
      MathExtra::norm3(contact_surfs[n].dr_force);
    } else {
      MathExtra::copy3(contact_surfs[n].surf_norm, contact_surfs[n].dr_force);
      if (w_ext > EPSILON)
        contact_surfs[n].weight_contribution *= w_int;
    }
  }
}

/* ----------------------------------------------------------------------
   Calculate forces
------------------------------------------------------------------------- */

void PairSurfGranular::calculate_3d_forces(std::vector<int> *composite_surfs)
{
  int n, m, j, k;
  double dot, residual[3];

  // Calculate properties of composite surface

  n = (*composite_surfs)[0];
  double max_overlap = contact_surfs[n].overlap;
  double max_overlap_ext = -BIG;
  double contact_at_max[3];
  MathExtra::copy3(contact_surfs[n].contact, contact_at_max);

  // Check if composite is hidden (convex) and calc min distance to any external features

  int hidden = 0;
  for (auto it = 0; it < composite_surfs->size(); it++) {
    n = (*composite_surfs)[it];
    j = contact_surfs[n].index;

    if (contact_surfs[n].convex_preceding_contact != -1) {
      // Hide composite surface, unless it's self-hidden (e.g. hairpin turns)
      hidden = 1;
      for (auto it2 = 0; it2 < composite_surfs->size(); it2++) {
        m = (*composite_surfs)[it2];
        k = contact_surfs[m].index;
        if (contact_surfs[n].convex_preceding_contact == k) {
          hidden = 0;
          break;
        }
      }

      if (hidden == 1)
        break;
    }

    if (contact_surfs[n].external)
      max_overlap_ext = MAX(max_overlap_ext, contact_surfs[n].overlap);
  }

  if (hidden) {
    for (auto it = 0; it < composite_surfs->size(); it++) {
      n = (*composite_surfs)[it];
      contact_surfs[n].hidden = 1;
    }
    return;
  }

  // Smooth int/ext based on arbitrary ratio of overlaps
  double w_int = 1.0;
  if (max_overlap_ext != -BIG) {
    w_int = max_overlap_ext / max_overlap;
    w_int = 1.0 - w_int * w_int;
  }
  double w_ext = 1.0 - w_int;

  for (auto it = 0; it < composite_surfs->size(); it++) {
    n = (*composite_surfs)[it];

    if (contact_surfs[n].external) {
      // If copying from another surf
      if (contact_surfs[n].copy_index_ext > -1) {
        m = contact_surfs[n].copy_index_ext;
        MathExtra::copy3(contact_surfs[m].dr_ext, contact_surfs[n].dr_ext);
      }

      // Interpolate between surface normal and dr_ext using smoothing
      MathExtra::scaleadd3(w_ext, contact_surfs[n].dr_ext, w_int, contact_surfs[n].surf_norm, contact_surfs[n].dr_force);
      MathExtra::norm3(contact_surfs[n].dr_force);
      contact_surfs[n].weight_contribution *= w_ext;
    } else {
      MathExtra::copy3(contact_surfs[n].surf_norm, contact_surfs[n].dr_force);
      contact_surfs[n].weight_contribution *= w_int;
    }
  }
}

/* ----------------------------------------------------------------------
   For external point contacts (j), adjust dr used for force calculation
------------------------------------------------------------------------- */

void PairSurfGranular::adjust_external_pt_flat_2d(int j, int k, int n, int m)
{
  // Only higher ranked contacts can correct
  if (n < m) return;

  // Already adjusted by higher ranked surf
  if (contact_surfs[n].rank_ext < m)
    return;

  // Otherwise, see if dr has component pointing into other (k) line
  //  If so, remove this component

  int pt, ptk;
  if (contact_surfs[n].flag == -1) {
    pt = 0;
  } else {
    pt = 3;
  }
  int *line = atom->line;
  double *pt_x = &endpts[line[j]][pt];
  double *ptk_x = &endpts[line[k]][0];

  double kline[3];
  MathExtra::sub3(ptk_x, pt_x, kline);

  if (MathExtra::lensq3(kline) < EPSILON) {
    ptk_x = &endpts[line[k]][3];
    MathExtra::sub3(ptk_x, pt_x, kline);
  }

  MathExtra::norm3(kline);

  // remove any component of dr that lies along kline

  double dr[3];
  MathExtra::copy3(contact_surfs[n].dr, dr);
  double dot = MathExtra::dot3(dr, kline);

  if (dot > 0) {
    MathExtra::scaleadd3(-dot, kline, dr, contact_surfs[n].dr_ext);
    MathExtra::norm3(contact_surfs[n].dr_ext);
    contact_surfs[n].rank_ext = m;
  }
}

/* ----------------------------------------------------------------------
   For external edge contacts (j), adjust dr used for force calculation
------------------------------------------------------------------------- */

void PairSurfGranular::adjust_external_edge_flat_3d(int j, int k, int n, int m)
{
  // Only higher priority contacts can correct
  if (n < m) return;

  // Already adjusted by higher ranked surf
  if (contact_surfs[n].rank_ext < m)
    return;

  // If superceding k contact is internal + flat, just use jnorm
  if (!contact_surfs[m].external) {
    MathExtra::copy3(contact_surfs[n].surf_norm, contact_surfs[n].dr_ext);
    contact_surfs[n].rank_ext = m;
    contact_surfs[n].copy_index_ext = -1;
  }
}

/* ----------------------------------------------------------------------
   For external corner contacts (j), adjust dr used for force calculation
------------------------------------------------------------------------- */

void PairSurfGranular::adjust_external_pt_flat_3d(int j, int k, int n, int m)
{
  // Only higher ranked contacts can correct
  if (n < m) return;

  // Already copying from higher ranked flat surf
  if (contact_surfs[n].flat_ext && contact_surfs[n].rank_ext < m)
    return;

  // If superceding k contact is internal + flat, just use jnorm
  if (!contact_surfs[m].external) {
    MathExtra::copy3(contact_surfs[n].surf_norm, contact_surfs[n].dr_ext);
    contact_surfs[n].rank_ext = m;
    contact_surfs[n].copy_index_ext = -1;
    contact_surfs[n].flat_ext = 1;
  }

    // If superceding k contact is a flat external edge, copy it's result
  if (contact_surfs[m].external && contact_surfs[m].flag > -4) {
    contact_surfs[n].copy_index_ext = m;
    contact_surfs[n].rank_ext = m;
    contact_surfs[n].flat_ext = 1;
  }
}

/* ----------------------------------------------------------------------
   For external corner contacts (j), adjust dr used for force calculation
------------------------------------------------------------------------- */

void PairSurfGranular::adjust_external_pt_nonflat_3d(int j, int k, int n, int m)
{
  // Only higher ranked contacts can correct
  if (n < m) return;

  // Already adjusted by higher ranked surf
  if (contact_surfs[n].rank_ext < m)
    return;

  // Skip if corrected by flat surf
  if (contact_surfs[n].flat_ext)
    return;

  // If superceding k contact is not flat, project into j-k norm plane
  //   remove component in that plane not along jnorm
  //   lower priority than FLAT correction

  double jnorm[3], knorm[3];
  MathExtra::copy3(contact_surfs[n].surf_norm, jnorm);
  MathExtra::copy3(contact_surfs[m].surf_norm, knorm);

  double jxk[3];
  MathExtra::cross3(jnorm, knorm, jxk);
  MathExtra::norm3(jxk);

  double dr[3], dr_proj[3], dr_keep[3];
  MathExtra::copy3(contact_surfs[n].dr, dr);
  double dot = MathExtra::dot3(dr, jxk);
  MathExtra::scale3(dot, jxk, dr_keep);
  MathExtra::sub3(dr, dr_keep, dr_proj);

  double magsq_proj = MathExtra::lensq3(dr_proj);

  dot = MathExtra::dot3(jnorm, dr_proj);
  MathExtra::scale3(dot, jnorm, dr_proj);
  double scale = magsq_proj / MathExtra::lensq3(dr_proj);

  double dr_tmp[3];
  MathExtra::scaleadd3(scale, dr_proj, dr_keep, dr_tmp);

  // When inline with edge, smooth it towards dr (like an actual point)
  //   estimate by dotting with edge vector
  //   (ideally would only calc once in calculate_forces)

  int pt = -1;
  int ptj1 = -1;
  int ptj2 = -1;
  if (contact_surfs[n].flag == -4) {
    pt = 0;
    ptj1 = 3;
    ptj2 = 6;
  } else if (contact_surfs[n].flag == -5) {
    pt = 3;
    ptj1 = 0;
    ptj2 = 6;
  } else if (contact_surfs[n].flag == -6) {
    pt = 6;
    ptj1 = 3;
    ptj2 = 0;
  }
  if (pt == -1) error->one(FLERR, "Bad geometry?");

  int *tri = atom->tri;
  double *pt_x = &corners[tri[j]][pt];
  double *ptj1_x = &corners[tri[j]][ptj1];
  double *ptj2_x = &corners[tri[j]][ptj2];

  double jline1[3], jline2[3];
  MathExtra::sub3(pt_x, ptj1_x, jline1);
  MathExtra::sub3(pt_x, ptj2_x, jline2);
  MathExtra::norm3(jline1);
  MathExtra::norm3(jline2);

  double dot1 = MathExtra::dot3(jline1, dr);
  double dot2 = MathExtra::dot3(jline2, dr);
  double w = MIN(fabs(dot1), fabs(dot2));

  MathExtra::scaleadd3(w, dr, 1.0 - w, dr_tmp, contact_surfs[n].dr_ext);
  MathExtra::norm3(contact_surfs[n].dr_ext);

  contact_surfs[n].rank_ext = m;
}
