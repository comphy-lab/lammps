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
#ifndef LDD_INDICATOR_H
#define LDD_INDICATOR_H

#include "pointers.h"

namespace LAMMPS_NS {

class LddIndicator : protected Pointers {
 public:
  double r0, rc;  // The min/max length scales (resp) where the indicator is defined [user]
  double rs, rb;  // extra length scales that may be relevent to the indicator [derived] e.g. the radius of small and big sphere inferred from r0 and rc in ldd_indicator_sphere
  double *coeffs; // Coefficients involved in the indicator function w(r) 
  double norm, invnorm; // Normalization [w] and its inverse
  int n_coeffs;  // length of coeffs
  int allocated; // 0 or 1, whether the indicator has been allocated

  LddIndicator(class LAMMPS * );
  virtual ~LddIndicator();
    
  virtual void init_coeffs(double, double, int) {} // Specific coeffs involved in w(r), overriden by type of indicator
  virtual double w(double ) { return 0; }  // w(r) given the coeffs (overriden)
  virtual double wp(double ) { return 0; } // w'(r) 
  virtual double wp2(double ) { return 0; } // w''(r) 

};

}  

#endif
