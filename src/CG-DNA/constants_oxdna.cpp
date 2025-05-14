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
/* ----------------------------------------------------------------------
   Contributing authors: Oliver Henrich (University of Strathclyde, Glasgow)
                         Kierran Falloon (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "constants_oxdna.h"

#include "update.h"

namespace LAMMPS_NS {

ConstantsOxdna::ConstantsOxdna(class LAMMPS *lmp) : Pointers(lmp)
{
  // set oxDNA units
  units = update->unit_style;
  real_flag = utils::strmatch(units.c_str(), "^real");
  if (real_flag) set_real_units();
}

// default to lj units
// oxDNA 1 parameters
double ConstantsOxdna::d_cback = -0.4;
double ConstantsOxdna::d_cstack = +0.34;
double ConstantsOxdna::d_cbase = +0.4;
// oxDNA 2 parameters
double ConstantsOxdna::d_cback_x = -0.34;
double ConstantsOxdna::d_cback_y = +0.3408;
double ConstantsOxdna::lambda_dh_one_prefactor = +0.3616455075438555;      // = C1
double ConstantsOxdna::qeff_dh_pf_one_prefactor = +0.08173808693529228;    // = C2
// oxRNA 2 parameters
double ConstantsOxdna::d_cback_z = +0.2;
double ConstantsOxdna::d_cstack_x_3p = +0.4;
double ConstantsOxdna::d_cstack_y_3p = +0.1;
double ConstantsOxdna::d_cstack_x_5p = +0.124906078525;
double ConstantsOxdna::d_cstack_y_5p = -0.00866274917473;

void ConstantsOxdna::set_real_units()
{
  // oxDNA 1 parameters in real units
  d_cback = -3.4072;
  d_cstack = +2.89612;
  d_cbase = +3.4072;
  // oxDNA 2 parameters in real units
  d_cback_x = -2.89612;
  d_cback_y = +2.9029344;
  lambda_dh_one_prefactor = +0.05624154892;        // = C1 * 8.518 * sqrt(k_B/4.142e-20)
  qeff_dh_pf_one_prefactor = +4.15079634587587;    // = C2 * 5.961689060210325 * 8.518
  // oxRNA 2 parameters in real units
  // d_cback_x = -3.4072 = d_cback for RNA
  d_cback_z = +1.7036;
  d_cstack_x_3p = +3.4072;
  d_cstack_y_3p = +0.8518;
  d_cstack_x_5p = +1.063949977;
  d_cstack_y_5p = -0.07378929747;
};

}    // namespace LAMMPS_NS
