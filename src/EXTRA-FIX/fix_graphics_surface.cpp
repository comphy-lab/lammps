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

#include "fix_graphics_surface.h"

#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "dump_image.h"
#include "error.h"
#include "input.h"
#include "math_extra.h"
#include "memory.h"
#include "modify.h"
#include "update.h"
#include "variable.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

namespace {
enum { NONE, SURFLOW, SURFMED, SURFHIGH, SURFMAX };
}    // namespace

/* ---------------------------------------------------------------------- */

FixGraphicsSurface::FixGraphicsSurface(LAMMPS *lmp, int narg, char **arg) :
    Fix(lmp, narg, arg), rstr(nullptr), imgobjs(nullptr), imgparms(nullptr)
{
  if (narg < 6) utils::missing_cmd_args(FLERR, "fix graphics/surface", error);

  global_freq = nevery;
  dynamic_group_allow = 1;

  // defaults
  varflag = 0;
  numobjs = 0;
  quality = SURFMED;
  binary = 0;
  pad = 0;
  rvar = -1;

  // parse mandatory args

  nevery = utils::inumeric(FLERR, arg[3], false, lmp);
  if (nevery <= 0) error->all(FLERR, 3, "Illegal fix graphics/surface nevery value {}", nevery);
  atype = utils::inumeric(FLERR, arg[4], false, lmp);
  if (atype < 0) error->all(FLERR, 4, "Illegal fix graphics/surface type value {}", atype);
  iso = utils::numeric(FLERR, arg[5], false, lmp);
  if (iso <= 0.0) error->all(FLERR, 5, "Illegal fix graphics/surface isovalue {}", iso);
  if (strstr(arg[6], "v_") == arg[6]) {
    varflag = 1;
    rstr = utils::strdup(arg[6] + 2);
  } else {
    rad = utils::numeric(FLERR, arg[6], false, lmp);
    if (rad <= 0) error->all(FLERR, 6, "Illegal fix graphics/surface radius value {}", rad);
  }

  // parse optional args

  int iarg = 7;
  while (iarg < narg) {
    if (strcmp(arg[iarg], "quality") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/surface quality", error);
      iarg += 2;
    } else if (strcmp(arg[iarg], "filename") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/surface filename", error);
      filename = arg[iarg + 1];
      iarg += 2;
    } else if (strcmp(arg[iarg], "binary") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/surface binary", error);
      binary = utils::logical(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else if (strcmp(arg[iarg], "pad") == 0) {
      if (iarg + 2 > narg) utils::missing_cmd_args(FLERR, "fix graphics/surface pad", error);
      pad = utils::inumeric(FLERR, arg[iarg + 1], false, lmp);
      iarg += 2;
    } else {
      error->all(FLERR, iarg, "Unknown fix graphics/surface keyword {}", arg[iarg]);
    }
  }

  memory->create(imgobjs, numobjs, "fix_graphics:imgobjs");
  memory->create(imgparms, numobjs, 8, "fix_graphics:imgparms");
}

/* ---------------------------------------------------------------------- */

FixGraphicsSurface::~FixGraphicsSurface()
{
  delete[] rstr;

  memory->destroy(imgobjs);
  memory->destroy(imgparms);
}

/* ---------------------------------------------------------------------- */

int FixGraphicsSurface::setmask()
{
  return END_OF_STEP;
}

/* ---------------------------------------------------------------------- */

void FixGraphicsSurface::init()
{
  int n = 0;
  if (rstr) {
    int ivar = input->variable->find(rstr);
    if (ivar < 0)
      error->all(FLERR, Error::NOLASTLINE,
                 "Variable name {} for fix graphics/surface does not exist", rstr);
    if ((input->variable->equalstyle(ivar) == 0) && (input->variable->atomstyle(ivar) == 0))
      error->all(FLERR, Error::NOLASTLINE,
                 "fix graphics/surface variable {} is not equal- or atom-style variable", rstr);
    rvar = ivar;
  }
  end_of_step();
}

/* ---------------------------------------------------------------------- */

void FixGraphicsSurface::end_of_step()
{
  // evaluate variable if necessary, wrap with clear/add

  if (varflag) modify->clearstep_compute();

  int n = 0;

  if (varflag) modify->addstep_compute((update->ntimestep / nevery) * nevery + nevery);
}

/* ----------------------------------------------------------------------
   provide graphics information to dump image
------------------------------------------------------------------------- */

int FixGraphicsSurface::image(int *&objs, double **&parms)
{
  objs = imgobjs;
  parms = imgparms;
  return numobjs;
}
