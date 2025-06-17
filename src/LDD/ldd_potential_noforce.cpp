/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ------------------------------------------------------
    This file is part of the USER-LDD package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */
#include "ldd_potential_noforce.h"

#include <cmath>
#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "atom.h"
#include "atom_masks.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "kspace.h"
#include "math_const.h"
#include "memory.h"
#include "neighbor.h"
#include "suffix.h"
#include "update.h"
#include "utils.h"


using namespace LAMMPS_NS;

LddPotentialNoForce::LddPotentialNoForce(class LAMMPS * lmp) : LddPotential(lmp)
{
  n_coeffs = 0;
  ptype_len = 7;
}

LddPotentialNoForce::~LddPotentialNoForce()
{
  if (allocated == 1)
  {
    memory->destroy(ptype);
  }
  allocated = 0;
}

void LddPotentialNoForce::allocate()
{
  memory->create(ptype,ptype_len,"ldd_potential:ptype");
  allocated = 1;
}

void LddPotentialNoForce::setup_potl(int ipt, int narg, char **arg)
{
  if (!allocated) allocate();

  sprintf(ptype,"noforce");
}

/* These last two should never be called in the case of noforce */
double LddPotentialNoForce::u(double rho)
{
  return 0.0;
}

double LddPotentialNoForce::f(double rho)
{
  return 0.0;
}

