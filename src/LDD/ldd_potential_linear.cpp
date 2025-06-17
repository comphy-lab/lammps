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

#include "ldd_potential_linear.h"

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

LddPotentialLinear::LddPotentialLinear(class LAMMPS * lmp) : LddPotential(lmp)
{
  n_coeffs = 2;
  ptype_len = 6;
}

LddPotentialLinear::~LddPotentialLinear()
{
  if (allocated == 1)
  {
    memory->destroy(coeffs);
    memory->destroy(ptype);
  }
  allocated = 0;
}

void LddPotentialLinear::allocate()
{
  memory->create(coeffs,n_coeffs,"ldd_potential:coeffs");
  memory->create(ptype,ptype_len,"ldd_potential:ptype");
  allocated = 1;
}

void LddPotentialLinear::setup_potl(int ipt, int narg, char **arg)
{
  if (!allocated) allocate();

  coeffs[0] = utils::numeric(FLERR,arg[ipt+2],false,lmp);
  coeffs[1] = utils::numeric(FLERR,arg[ipt+3],false,lmp);
  sprintf(ptype,"linear");
}

double LddPotentialLinear::u(double rho)
{
  return (coeffs[0] * rho + coeffs[1]);
}

double LddPotentialLinear::f(double rho)
{
  return -coeffs[0];
}

