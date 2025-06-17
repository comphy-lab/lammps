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
#include "ldd_indicator_smooth.h"

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

LddIndicatorSmooth::LddIndicatorSmooth(class LAMMPS * lmp) : LddIndicator(lmp)
{
  n_coeffs = 7;
}

LddIndicatorSmooth::~LddIndicatorSmooth()
{
  if (allocated == 1)
  {
    memory->destroy(coeffs);
  }
  allocated = 0;
}

void LddIndicatorSmooth::allocate()
{
  memory->create(coeffs,n_coeffs,"ldd_indicator:coeffs");
  allocated = 1;
}

void LddIndicatorSmooth::init_coeffs(double a, double b, int dim)
{
  if (!allocated) allocate();
  r0 = a;
  rc = b;
  coeffs[6] = (pow(r0,5) - pow(rc,5)) / 120.0 -
              (pow(r0,3) - pow(rc,3)) * r0 * rc / 24.0 +
              (r0 - rc) * pow(r0,2) * pow(rc,2) / 12.0;
  coeffs[0] = (- pow(rc,5) / 120.0 + r0 * pow(rc,4) / 24.0
               - pow(r0,2) * pow(rc,3) / 12.0) / coeffs[6];
  coeffs[1] = (pow(r0,2) * pow(rc,2) / 4.0) / coeffs[6];
  coeffs[2] = (-rc * r0 * (r0 + rc) / 4.0) / coeffs[6];
  coeffs[3] = ((pow(r0,2) + 4.0 * r0 * rc + pow(rc,2)) / 12.0) / coeffs[6];
  coeffs[4] = (-(r0 + rc) / 8.0) / coeffs[6];
  coeffs[5] = 0.05 / coeffs[6];
  switch (dim)
  {
  case 2:
    norm = 2.0 * MY_PI * (pow(r0,2) / 2.0  +
            coeffs[5] / 7.0 * (pow(rc,7)-pow(r0,7)) +
            coeffs[4] / 6.0 * (pow(rc,6)-pow(r0,6)) +
            coeffs[3] / 5.0 * (pow(rc,5)-pow(r0,5)) +
            coeffs[2] / 4.0 * (pow(rc,4)-pow(r0,4)) +
            coeffs[1] / 3.0 * (pow(rc,3)-pow(r0,3)) +
            coeffs[0] / 2.0 * (pow(rc,2)-pow(r0,2)));
    break;
  case 3:
    norm = 4.0 * MY_PI * (pow(r0,3) / 3.0  +
            coeffs[5] / 8.0 * (pow(rc,8)-pow(r0,8)) +
            coeffs[4] / 7.0 * (pow(rc,7)-pow(r0,7)) +
            coeffs[3] / 6.0 * (pow(rc,6)-pow(r0,6)) +
            coeffs[2] / 5.0 * (pow(rc,5)-pow(r0,5)) +
            coeffs[1] / 4.0 * (pow(rc,4)-pow(r0,4)) +
            coeffs[0] / 3.0 * (pow(rc,3)-pow(r0,3)));
    break;
  }
  invnorm = 1.0 / norm;               
}

double LddIndicatorSmooth::w(double r)
{
  if (r < r0) { return invnorm; }
  if (r > rc) { return 0.0; }
  return ((coeffs[0] + 
           coeffs[1] * r + 
           coeffs[2] * pow(r,2) + 
           coeffs[3] * pow(r,3) + 
           coeffs[4] * pow(r,4) + 
           coeffs[5] * pow(r,5)) /* coeffs[6]*/) * invnorm;
}

double LddIndicatorSmooth::wp(double r)
{
  if (r < r0) { return 0.0; }
  if (r > rc) { return 0.0; }
  return ((      coeffs[1] + 
           2.0 * coeffs[2] * r +
           3.0 * coeffs[3] * pow(r,2) + 
           4.0 * coeffs[4] * pow(r,3) +
           5.0 * coeffs[5] * pow(r,4)) /* coeffs[6]*/) * invnorm;
}

double LddIndicatorSmooth::wp2(double r)
{
  if (r < r0) { return 0.0; }
  if (r > rc) { return 0.0; }
  return (( 2.0 * coeffs[2] + 
            6.0 * coeffs[3] * r +
           12.0 * coeffs[4] * pow(r,2) + 
           20.0 * coeffs[5] * pow(r,3)) /* coeffs[6]*/ ) * invnorm;
}
