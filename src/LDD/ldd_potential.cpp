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
#include "ldd_potential.h"

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
#include "neigh_list.h"
#include "suffix.h"
#include "update.h"
#include "utils.h"


using namespace LAMMPS_NS;

LddPotential::LddPotential(class LAMMPS * lmp) : Pointers(lmp)
{
  allocated = n_coeffs = ptype_len = 0;
  coeffs = NULL;
  ptype = NULL;
  potl_table.n_pts = 0;
  potl_table.dr = 0.0;
  potl_table.r = NULL;
  potl_table.u = NULL;
  potl_table.u2 = NULL;
  potl_table.f = NULL;
  potl_table.f2 = NULL;
}

LddPotential::~LddPotential()
{
  allocated = 0;
}


void LddPotential::read_table_file(char *fnm, bool bspline)
{
  potl_table.n_pts = 0;
  FILE *fp = fopen(fnm,"r");
  if (!fp)
  {
    char *errmsg = (char *) calloc(100,sizeof(char));
    sprintf(errmsg,"ERROR: unable to open file %s for reading",fnm);
    error->all(FLERR,errmsg);
  }
  char line[1000];
  while (fgets(line,1000,fp) != NULL)
  {
    ++(potl_table.n_pts);
  }
  rewind(fp);

  memory->create(potl_table.r,potl_table.n_pts,"LDpotl:r");
  memory->create(potl_table.u,potl_table.n_pts,"LDpotl:u");
  memory->create(potl_table.f,potl_table.n_pts,"LDpotl:f");
  if (bspline)
  {
    memory->create(potl_table.u2,potl_table.n_pts,"LDpotl:u2");
    memory->create(potl_table.f2,potl_table.n_pts,"LDpotl:f2");
  }

  for (int i = 0; i < potl_table.n_pts; ++i)
  {
    fgets(line,1000,fp);
    float rt, ut, ft;
    int test_sscanf = sscanf(line," %f %f %f ", &rt, &ut, &ft);
    if (test_sscanf != 3)
    {
      char *errmsg = (char *) calloc(100, sizeof(char));
      sprintf(errmsg,"ERROR: unable to read rho u f from line %d in file %s",
                       i,fnm);
      error->all(FLERR,errmsg);
    }
    potl_table.r[i] = (double) rt;
    potl_table.u[i] = (double) ut;
    potl_table.f[i] = (double) ft;
  } 
  fclose(fp);
  potl_table.dr = potl_table.r[1] - potl_table.r[0];
  for (int i = 0; i < potl_table.n_pts - 1; ++i)
  {
    if (fabs(potl_table.r[i+1] - potl_table.r[i] - potl_table.dr) > potl_table.dr / 100.0)
    {
      char *errmsg = (char *) calloc(200, sizeof(char));
      sprintf(errmsg,"ERROR: in file %s\nr[1] - r[0] = %g - %g = %g\nr[%d] - r[%d] = %g - %g = %g",
            fnm,potl_table.r[1],potl_table.r[0],potl_table.dr,
            i+1,i,potl_table.r[i+1],potl_table.r[i],potl_table.r[i+1] - potl_table.r[i]);
      error->all(FLERR,errmsg);
    }
  }
}

int LddPotential::get_table_index(double r)
{
  if (r > potl_table.r[potl_table.n_pts-1])
  {
    char *errmsg = (char *) calloc(1000,sizeof(char));
    sprintf(errmsg,"ERROR: local density = %g > %g = table hi",r,potl_table.r[potl_table.n_pts-1]);
    error->one(FLERR,errmsg);
  }
  if (r < potl_table.r[0])
  {
    char *errmsg = (char *) calloc(1000,sizeof(char));
    sprintf(errmsg,"ERROR: local density = %g < %g = table lo",r,potl_table.r[0]);
    error->one(FLERR,errmsg);
  }
  return floor( (r-potl_table.r[0]) / (potl_table.dr) );
}

double LddPotential::calc_A_table(double r, int idx)
{
  return ( (potl_table.r[idx+1] - r) / (potl_table.dr) );
}
