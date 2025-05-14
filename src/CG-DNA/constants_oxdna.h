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
  static double get_d_cbk() { return d_cbk; }
  static double get_d_cstk() { return d_cstk; }
  static double get_d_cbs() { return d_cbs; }

  // oxDNA 2 getters
  static double get_d_cbk_x() { return d_cbk_x; }
  static double get_d_cbk_y() { return d_cbk_y; }
  static double get_lambda_dh_one_prefactor() { return lambda_dh_one_prefactor; }
  static double get_qeff_dh_pf_one_prefactor() { return qeff_dh_pf_one_prefactor; }

  // oxRNA 2 getters
  static double get_d_cbk_z() { return d_cbk_z; }
  static double get_d_cstk_x_3p() { return d_cstk_x_3p; }
  static double get_d_cstk_y_3p() { return d_cstk_y_3p; }
  static double get_d_cstk_x_5p() { return d_cstk_x_5p; }
  static double get_d_cstk_y_5p() { return d_cstk_y_5p; }

 private:
  std::string units;
  bool real_flag;
  void set_real_units();

  // oxDNA 1 parameters
  static double d_cbk, d_cstk, d_cbs;

  // oxDNA 2 parameters
  static double d_cbk_x, d_cbk_y;
  static double lambda_dh_one_prefactor, qeff_dh_pf_one_prefactor;

  // oxRNA 2 parameters
  static double d_cbk_z;
  static double d_cstk_x_3p, d_cstk_y_3p;
  static double d_cstk_x_5p, d_cstk_y_5p;
};

}    // namespace LAMMPS_NS

#endif
