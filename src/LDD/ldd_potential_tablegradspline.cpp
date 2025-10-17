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
/* ------------------------------------------------------
    This file is part of the LDD package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    and Maria C. Lesniewski, mjl6766@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */
#include "ldd_potential_tablegradspline.h"

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

LddPotentialTableGradSpline::LddPotentialTableGradSpline(class LAMMPS * lmp) : LddPotential(lmp)
{
  n_coeffs = 1;
  ptype_len = 18;
}

LddPotentialTableGradSpline::~LddPotentialTableGradSpline()
{
  if (allocated == 1)
  {
    memory->destroy(coeffs);
    memory->destroy(ptype);
    memory->destroy(table_fnm);
    memory->destroy(potl_table.r);
    memory->destroy(potl_table.u);
    memory->destroy(potl_table.f);
    memory->destroy(potl_table.u2);
    memory->destroy(potl_table.f2);
  }
  allocated = 0;
}

void LddPotentialTableGradSpline::allocate()
{
  memory->create(coeffs,n_coeffs,"ldd_potential:coeffs");
  memory->create(ptype,ptype_len,"ldd_potential:ptype");
  memory->create(table_fnm,100,"ldd_potential:table_fnm");
  allocated = 1;
}

/*******************************
Calculates second derivatives at knot points. stores them in y2[]. pass in
knot points x[] and y[]. pass in the number of knot points n. pass in first
derivative you want for the end points in yp1 and ypn. if yp1 and/or ypn
are > 0.99e30, then second derivative is set to zero for that boundary.
Copied from "Numerical Recipes in C" second edition.
*******************************/
void LddPotentialTableGradSpline::spline(double x[], double y[], int n, double yp1, double ypn, double y2[])
{
  int i, j;
  double p, qn, sig, un, u[n];
  if (yp1 > 0.99e30) y2[0] = u[0] = 0.0;
  else { y2[0] = -0.5; u[0] = (3.0/(x[1]-x[0])) * ((y[1]-y[0])/(x[1]-x[0]) - yp1); }

  for (i=1; i < n-1; i++)
  {
    sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p = sig*y2[i-1]+2.0;
    y2[i] = (sig-1.0)/p;
    u[i] = (y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i] = (6.0*u[i]/(x[i+1]-x[i-1]) - sig*u[i-1]) / p;
  }

  if (ypn > 0.99e30) qn = un = 0.0;
  else { qn = 0.5; un = (3.0/(x[n-1]-x[n-2])) * (ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2])); }

  y2[n-1] = (un - qn * u[n-2])/(qn * y2[n-2] + 1.0);

  for (j=n-2; j>1;j--) { y2[j] = y2[j] * y2[j+1] + u[j]; }
}

void LddPotentialTableGradSpline::setup_potl(int ipt, int narg, char **arg)
{
  if (!allocated) allocate();

  if (narg <= ipt + 1)
  {
    error->all(FLERR,"ERROR: unable to read filename following table_lin");
    exit(EXIT_FAILURE);
  }
  strcpy(table_fnm,arg[ipt+2]);
  read_table_file(arg[ipt+2],true);


  sprintf(ptype,"table/spline");


  spline (&(potl_table.r[0]), &(potl_table.u[0]), potl_table.n_pts, 0.999e30,
                                               0.999e30, &(potl_table.u2[0]));
  spline (&(potl_table.r[0]), &(potl_table.f[0]), potl_table.n_pts, 0.999e30,
                                               0.999e30, &(potl_table.f2[0]));
}

/* Evaluate cubic spline. Modified from "Numerical Recipes in C" Second Edition */
double LddPotentialTableGradSpline::splint(double x0, double x1, double y0,
                          double y1, double y20, double y21,
                          double dr, double x, double a, double b)
{
  if (fabs(x - x0) <= 1e-8) { return y0; }
  if (fabs(x1 - x) <= 1e-8) { return y1; }
  return (a * y0 + b * y1 + ((a*a*a-a) * y20 + (b*b*b-b) * y21) * dr * dr / 6.0);
  //MCL 04.24.25 this is eqn 3.3.3 -> |a is A | b is B| a^3-a/6(dr**2) is C | 1/6(B^3-B)dr is D | y20 is y''_{j} | y21 is y''_{j+1}
}
//MCL adding an interpolator for the derivative
double LddPotentialTableGradSpline::dsplint(double x0, double x1, double y0,
                          double y1, double y20, double y21,
                          double dr, double x, double a, double b)
{
  //MCL -> This will be equation 3.3.5
  double dy_dx = 1/dr*(y1 - y0) - (3*a*a - 1)/6*dr*y20 + (3*b*b-1)/6*dr*y21;
  return dy_dx;

}

double LddPotentialTableGradSpline::u(double rho)
{
  // Handle this case separately
  if (rho == potl_table.r[potl_table.n_pts-1])
  {
    return potl_table.u[potl_table.n_pts-1];
  }
  int idx = get_table_index(rho);
  double A = calc_A_table(rho, idx);

  return splint(potl_table.r[idx],potl_table.r[idx+1],
                potl_table.u[idx],potl_table.u[idx+1],
                potl_table.u2[idx],potl_table.u2[idx+1],
                potl_table.dr,rho,A,1.0-A);
}

double LddPotentialTableGradSpline::f(double rho)
{
  if (rho == potl_table.r[potl_table.n_pts-1])
  {
    return potl_table.f[potl_table.n_pts-1];
  }
  int idx = get_table_index(rho);
  double A = calc_A_table(rho, idx);

  return -dsplint(potl_table.r[idx],potl_table.r[idx+1],
                potl_table.u[idx],potl_table.u[idx+1],
                potl_table.u2[idx],potl_table.u2[idx+1],
                potl_table.dr,rho,A,1.0-A);
}
