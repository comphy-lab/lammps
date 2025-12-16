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

#include "fix_graphics.h"

#include "comm.h"
#include "domain.h"
#include "dump_image.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

enum { SPHERE, ARROW, PROGBAR };
enum { X, Y, Z };

/* ---------------------------------------------------------------------- */

FixGraphics::FixGraphics(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 4) utils::missing_cmd_args(FLERR, "fix graphics", error);

  // parse mandatory arg

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix graphics nevery value {}", nevery);
  global_freq = nevery;

  numobjs = 0;
  varflag = 0;

  int iarg = 4;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "sphere") == 0) {
      if (iarg + 6 > narg) utils::missing_cmd_args(FLERR, "fix graphics sphere", error);
      SphereItem sphere{SPHERE, 1, {0.0, 0.0, 0.0}, 0.0, 0, 0, 0, 0, -1, -1, -1, -1};
      sphere.type = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      if (strstr(arg[iarg + 2], "v_") == arg[iarg + 2]) {
        varflag = 1;
        sphere.xstr = utils::strdup(arg[iarg + 2] + 2);
      } else
        sphere.pos[0] = utils::numeric(FLERR, arg[iarg + 2], false, lmp);
      if (strstr(arg[iarg + 3], "v_") == arg[iarg + 3]) {
        varflag = 1;
        sphere.ystr = utils::strdup(arg[iarg + 3] + 2);
      } else
        sphere.pos[1] = utils::numeric(FLERR, arg[iarg + 3], false, lmp);
      if (strstr(arg[iarg + 4], "v_") == arg[iarg + 4]) {
        varflag = 1;
        sphere.zstr = utils::strdup(arg[iarg + 4] + 2);
      } else
        sphere.pos[2] = utils::numeric(FLERR, arg[iarg + 4], false, lmp);
      if (strstr(arg[iarg + 5], "v_") == arg[iarg + 5]) {
        varflag = 1;
        sphere.dstr = utils::strdup(arg[iarg + 5] + 2);
      } else
        sphere.diameter = 2.0 * utils::numeric(FLERR, arg[iarg + 5], false, lmp);
      GraphicsItem g;
      g.sphere = sphere;
      items.emplace_back(g);
      ++numobjs;
      iarg += 6;
    } else {
      error->all(FLERR, iarg, "Unknown fix graphics keyword {}", arg[iarg]);
    }
  }
  memory->create(imgobjs, numobjs, "fix_graphics:imgobjs");
  memory->create(imgparms, numobjs, 10, "fix_graphics:imgparms");
}

/* ---------------------------------------------------------------------- */

FixGraphics::~FixGraphics()
{
  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphics::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphics::init()
{
  int n = 0;
  for (auto &gi : items) {
    if (gi.style == SPHERE) {
      imgobjs[n] = DumpImage::SPHERE;
      imgparms[n][0] = gi.sphere.type;
      if (gi.sphere.xstr) {
        int ivar = input->variable->find(gi.sphere.xstr);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.sphere.xstr);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.sphere.xstr);
        gi.sphere.xvar = ivar;
      }
      if (gi.sphere.ystr) {
        int ivar = input->variable->find(gi.sphere.ystr);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.sphere.ystr);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.sphere.ystr);
        gi.sphere.yvar = ivar;
      }
      if (gi.sphere.zstr) {
        int ivar = input->variable->find(gi.sphere.zstr);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.sphere.zstr);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.sphere.zstr);
        gi.sphere.zvar = ivar;
      }
      if (gi.sphere.dstr) {
        int ivar = input->variable->find(gi.sphere.dstr);
        if (ivar < 0)
          error->all(FLERR, Error::NOLASTLINE, "Variable name {} for fix graphics does not exist",
                     gi.sphere.dstr);
        if (input->variable->equalstyle(ivar) == 0)
          error->all(FLERR, Error::NOLASTLINE,
                     "Fix graphics variable {} is not equal-style variable", gi.sphere.dstr);
        gi.sphere.dvar = ivar;
      }
      ++n;
    } else if (gi.style == ARROW) {
    } else if (gi.style == PROGBAR) {
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixGraphics::end_of_step()
{
  // evaluate variable if necessary, wrap with clear/add

  if (varflag) modify->clearstep_compute();

  int n = 0;
  for (auto &gi : items) {
    if (gi.style == SPHERE) {
      if (gi.sphere.xstr) gi.sphere.pos[0] = input->variable->compute_equal(gi.sphere.xvar);
      if (gi.sphere.ystr) gi.sphere.pos[1] = input->variable->compute_equal(gi.sphere.yvar);
      if (gi.sphere.zstr) gi.sphere.pos[2] = input->variable->compute_equal(gi.sphere.zvar);
      if (gi.sphere.dstr) gi.sphere.diameter = 2.0 * input->variable->compute_equal(gi.sphere.dvar);
      imgparms[n][1] = gi.sphere.pos[0];
      imgparms[n][2] = gi.sphere.pos[1];
      imgparms[n][3] = gi.sphere.pos[2];
      imgparms[n][4] = gi.sphere.diameter;
      ++n;
    } else if (gi.style == ARROW) {
    } else if (gi.style == PROGBAR) {
    }
  }

  if (varflag) modify->addstep_compute((update->ntimestep / nevery) * nevery + nevery);
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphics::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
