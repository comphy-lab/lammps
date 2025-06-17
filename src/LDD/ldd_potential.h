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
#ifndef LDD_POTENTIAL_H
#define LDD_POTENTIAL_H

#include "pointers.h"

namespace LAMMPS_NS {

class LddPotential : protected Pointers {
 public:
  int allocated;
  int n_coeffs;
  double *coeffs;
  int ptype_len;
  char *ptype;
  char *table_fnm;
  
  struct t_table {
    int n_pts;
    double dr;
    double *r;
    double *u;
    double *u2;
    double *f;
    double *f2;
  };

  t_table potl_table;

  LddPotential(class LAMMPS *);
  virtual ~LddPotential();
  
  virtual void setup_potl(int, int, char **) {}
  virtual double u(double ) { return 0; }
  virtual double f(double ) { return 0; }

  void read_table_file(char *, bool );
  int get_table_index(double );
  double calc_A_table(double , int );

};

}  

#endif
