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
#include "ldd_indicator_lucy.h"

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

LddIndicatorLucy::LddIndicatorLucy(class LAMMPS * lmp) : LddIndicator(lmp)
{
  n_coeffs = 5;
}

LddIndicatorLucy::~LddIndicatorLucy()
{
  if (allocated == 1)
  {
    memory->destroy(coeffs);
  }
  allocated = 0;
}

void LddIndicatorLucy::allocate()
{
  memory->create(coeffs,n_coeffs,"ldd_indicator:coeffs");
  allocated = 1;
}

void LddIndicatorLucy::init_coeffs(double a, double b, int dim)
{
  if (!allocated) allocate();
  r0 = a;
  rc = b;
  coeffs[0] = 1.0;
  coeffs[1] = 0.0;
  coeffs[2] = -6.0 / pow(rc,2);
  coeffs[3] = 8.0 / pow(rc,3);
  coeffs[4] = -3.0 / pow(rc,4);
  switch (dim)
  {
    case 2:
      norm = 0.2 * MY_PI * pow(rc,2);
      break;
    case 3:
      norm = 16.0 * MY_PI / 105.0 * pow(rc,3); 
      break;
  }
  invnorm = 1.0 / norm;
}

double LddIndicatorLucy::w(double r)
{
  if (r >= rc) { return 0.0; }
  double w = (coeffs[0] + 
              coeffs[2] * pow(r,2) + 
              coeffs[3] * pow(r,3) + 
              coeffs[4] * pow(r,4)) * invnorm;
/* MRD 5.10.2021
 * I was running into an issue where rc = 3.6,
 * r ~= 3.59999 (i.e. very close to 3.6),
 * and w(r) was getting returned as -9.94151e-18,
 * -1.9883e-17, or -2.98245e-17.
 * when this happened for particles with only one neighbor,
 * the overall LD was negative, causing problems with tabulated potentials
 * that started at 0.
 */
  return (w > 0.0 ? w : 0.0);
}

double LddIndicatorLucy::wp(double r)
{
  if (r >= rc) { return 0.0; }
  return (  2.0 * coeffs[2] * r 
          + 3.0 * coeffs[3] * pow(r,2)
          + 4.0 * coeffs[4] * pow(r,3)) * invnorm;
}

double LddIndicatorLucy::wp2(double r)
{
  if (r > rc) { return 0.0; }
  return (  2.0 * coeffs[2] + 
            6.0 * coeffs[3] * r +
           12.0 * coeffs[4] * pow(r,2)) * invnorm;
}
