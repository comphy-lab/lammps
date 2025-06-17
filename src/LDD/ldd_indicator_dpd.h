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
#ifdef LDD_INDICATOR_CLASS
// clang-format off
LddIndicatorStyle(dpd,LddIndicatorDpd);
// clang-format on
#else

#ifndef LDD_INDICATOR_DPD_H
#define LDD_INDICATOR_DPD_H

#include "ldd_indicator.h"

namespace LAMMPS_NS {

class LddIndicatorDpd : public LddIndicator {
  public:
  
    LddIndicatorDpd(class LAMMPS *);
    ~LddIndicatorDpd();
    virtual void init_coeffs(double, double, int);
    virtual double w(double );
    virtual double wp(double );
    virtual double wp2(double );

  protected:
    virtual void allocate();
};

} 

#endif
#endif
