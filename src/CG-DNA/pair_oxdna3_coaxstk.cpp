// clang-format off
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
/* ----------------------------------------------------------------------
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "pair_oxdna3_coaxstk.h"
#include "nucleotide_oxdna.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairOxdna3Coaxstk::PairOxdna3Coaxstk(LAMMPS *lmp) : PairOxdna2Coaxstk(lmp)
{

  // sequence-specific coaxial stacking strength
  // A:0 C:1 G:2 T:3, 3'- [i][j] -5'

  eta_cxst[0][0] = 1.00000;
  eta_cxst[1][0] = 1.00000;
  eta_cxst[2][0] = 1.00000;
  eta_cxst[3][0] = 1.00000;

  eta_cxst[0][1] = 1.00000;
  eta_cxst[1][1] = 1.00000;
  eta_cxst[2][1] = 1.00000;
  eta_cxst[3][1] = 1.00000;

  eta_cxst[0][2] = 1.00000;
  eta_cxst[1][2] = 1.00000;
  eta_cxst[2][2] = 1.00000;
  eta_cxst[3][2] = 1.00000;

  eta_cxst[0][3] = 1.00000;
  eta_cxst[1][3] = 1.00000;
  eta_cxst[2][3] = 1.00000;
  eta_cxst[3][3] = 1.00000;

  single_enable = 0;
  writedata = 0;
  trim_flag = 0;
}

/* -----------------------------------------------------------------------
    compute vector COM-stacking interaction site in oxDNA3
-------------------------------------------------------------------------- */
inline void PairOxdna3Coaxstk::compute_stacking_site(double e1[3], double /*e2*/[3],
    double /*e3*/[3], double rstk[3]) const
{
  NucleotideOxdna3 oxdna3;
  oxdna3.stacking_site(e1, NULL, NULL, rstk);
}
