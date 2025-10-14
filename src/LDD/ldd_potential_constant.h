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
#ifdef LDD_POTENTIAL_CLASS
// clang-format off
LddPotentialStyle(constant,LddPotentialConstant);
// clang-format on
#else

#ifndef LDD_POTENTIAL_CONSTANT
#define LDD_POTENTIAL_CONSTANT

#include "ldd_potential.h"

namespace LAMMPS_NS { 

class LddPotentialConstant : public LddPotential {
  public:
 
    LddPotentialConstant(class LAMMPS *);
    ~LddPotentialConstant(); 

    void setup_potl(int, int, char **) override;
    double u(double ) override; 
    double f(double ) override; 

  protected:
    virtual void allocate();

};
}

#endif 
#endif
