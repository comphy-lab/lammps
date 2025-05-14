/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef CONSTANTS_OXDNA_H
#define CONSTANTS_OXDNA_H

#include "pointers.h"

namespace LAMMPS_NS {

class ConstantsOxdna : protected Pointers {
 public:
  ConstantsOxdna(class LAMMPS *lmp);
  virtual ~ConstantsOxdna(){};

  // oxDNA 1 getters
  static double get_d_cback() { return d_cback; }
  static double get_d_cstack() { return d_cstack; }
  static double get_d_cbase() { return d_cbase; }

  // oxDNA 2 getters
  static double get_d_cback_x() { return d_cback_x; }
  static double get_d_cback_y() { return d_cback_y; }
  static double get_lambda_dh_one_prefactor() { return lambda_dh_one_prefactor; }
  static double get_qeff_dh_pf_one_prefactor() { return qeff_dh_pf_one_prefactor; }

  // oxRNA 2 getters
  static double get_d_cback_z() { return d_cback_z; }
  static double get_d_cstack_x_3p() { return d_cstack_x_3p; }
  static double get_d_cstack_y_3p() { return d_cstack_y_3p; }
  static double get_d_cstack_x_5p() { return d_cstack_x_5p; }
  static double get_d_cstack_y_5p() { return d_cstack_y_5p; }

 private:
  std::string units;
  bool real_flag;
  void set_real_units();

  // oxDNA 1 parameters
  static double d_cback, d_cstack, d_cbase;

  // oxDNA 2 parameters
  static double d_cback_x, d_cback_y;
  static double lambda_dh_one_prefactor, qeff_dh_pf_one_prefactor;

  // oxRNA 2 parameters
  static double d_cback_z;
  static double d_cstack_x_3p, d_cstack_y_3p;
  static double d_cstack_x_5p, d_cstack_y_5p;
};

}    // namespace LAMMPS_NS

#endif
