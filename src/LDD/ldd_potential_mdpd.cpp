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
/*
 * 
 * This LD potential type is the one employed in many MDPD papers, such as
 * Warren. Phys. Rev. E. 2003. 68, 066702.
 * Ghoufi, Malfreyt. Phys. Rev. E. 2011. 83, 051601.
 *
 */

#include "ldd_potential_mdpd.h"

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
using namespace MathConst;

LddPotentialMdpd::LddPotentialMdpd(class LAMMPS * lmp) : LddPotential(lmp)
{
  n_coeffs = 1;
  ptype_len = 4;
}

LddPotentialMdpd::~LddPotentialMdpd()
{
  if (allocated == 1)
  {
    memory->destroy(coeffs);
    memory->destroy(ptype);
  }
  allocated = 0;
}

void LddPotentialMdpd::allocate()
{
  memory->create(coeffs,n_coeffs,"ldd_potential:coeffs");
  memory->create(ptype,ptype_len,"ldd_potential:ptype");
  allocated = 1;
}

void LddPotentialMdpd::setup_potl(int ipt, int narg, char **arg)
{
  if (!allocated) allocate();

  coeffs[0] = utils::numeric(FLERR,arg[ipt+2],false,lmp); // B
  sprintf(ptype,"mdpd");  
  double rc = -1.0;
  int i;
  // The potential normalization depends on the indicator function
  for (i = 0; i < narg; ++i)
  {
    if (strcmp(arg[i],"indicator") == 0)
    {
      if (i+3 < narg)
      {
        rc = utils::numeric(FLERR,arg[i+3],false,lmp);
      }
    }
  }

  if (rc == -1.0) // Unable to find rc
  {
    error->all(FLERR,"Unable to find rc in LDD pair coeff line.");
  } 
  coeffs[1] = MY_PI * coeffs[0] * pow(rc,4) / 30.0;
}

double LddPotentialMdpd::u(double rho)
{
  return (coeffs[1] * rho * rho);
}

double LddPotentialMdpd::f(double rho)
{
  return (-2.0 * coeffs[1] * rho);
}

