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

#ifndef NUCLEOTIDE_OXDNA_H
#define NUCLEOTIDE_OXDNA_H

# include <cstdio>
#include "pointers.h"
#include "constants_oxdna.h"

namespace LAMMPS_NS {

template <class Derived>
class NucleotideOxdna {
 public:
  NucleotideOxdna(){};
  virtual ~NucleotideOxdna(){};
  void backbone_site() { static_cast<Derived*>(this); }
  void stacking_site() { static_cast<Derived*>(this); }
  template <int N>
  void base_site()     { static_cast<Derived*>(this); }
};

class NucleotideOxdna1 : public NucleotideOxdna<NucleotideOxdna1> {
 public:
  NucleotideOxdna1(){};
  virtual ~NucleotideOxdna1(){};
  void backbone_site(double *, double *, double *, double *);
  void stacking_site();
  template <int N>
  void base_site() {};
};

class NucleotideOxdna2 : public NucleotideOxdna<NucleotideOxdna2> {
 public:
  NucleotideOxdna2(){};
  virtual ~NucleotideOxdna2(){};
  void backbone_site(double *, double *, double *, double *);
  void stacking_site();
  template <int N>
  void base_site() {};
};

class NucleotideOxdna3 : public NucleotideOxdna<NucleotideOxdna3> {
 public:
  NucleotideOxdna3(){};
  virtual ~NucleotideOxdna3(){};
  void backbone_site();
  void stacking_site();
  template <int N>
  void base_site() {};
};

class NucleotideOxrna2 : public NucleotideOxdna<NucleotideOxrna2> {
 public:
  NucleotideOxrna2(){};
  virtual ~NucleotideOxrna2(){};
  void backbone_site(double *, double *, double *, double *);
  void stacking_site();
  template <int N>
  void base_site() {};
};

/* ----------------------------------------------------------------------
   oxDNA1 nucleotide
------------------------------------------------------------------------- */
inline void NucleotideOxdna1::backbone_site(double e1[3], double /*e2*/[3],
  double /*e3*/[3], double rbk[3]) {

  double dx_cbk_oxdna1 = ConstantsOxdna::get_dx_cbk_oxdna1();

  rbk[0] = dx_cbk_oxdna1 * e1[0];
  rbk[1] = dx_cbk_oxdna1 * e1[1];
  rbk[2] = dx_cbk_oxdna1 * e1[2];
}
inline void NucleotideOxdna1::stacking_site() {
  printf("oxDNA1 stacking site\n");
}
template <>
inline void NucleotideOxdna1::base_site<0>() { printf("oxDNA1 base site type 0\n"); }
template <>
inline void NucleotideOxdna1::base_site<1>() { printf("oxDNA1 base site type 1\n"); }
template <>
inline void NucleotideOxdna1::base_site<2>() { printf("oxDNA1 base site type 2\n"); }
template <>
inline void NucleotideOxdna1::base_site<3>() { printf("oxDNA1 base site type 3\n"); }

/* ----------------------------------------------------------------------
   oxDNA2 nucleotide
------------------------------------------------------------------------- */

inline void NucleotideOxdna2::backbone_site(double e1[3], double e2[3],
  double /*e3*/[3], double rbk[3]) {

  double dx_cbk_oxdna2 = ConstantsOxdna::get_dx_cbk_oxdna2();
  double dy_cbk_oxdna2 = ConstantsOxdna::get_dy_cbk_oxdna2();

  rbk[0] = dx_cbk_oxdna2 * e1[0] + dy_cbk_oxdna2 * e2[0];
  rbk[1] = dx_cbk_oxdna2 * e1[1] + dy_cbk_oxdna2 * e2[1];
  rbk[2] = dx_cbk_oxdna2 * e1[2] + dy_cbk_oxdna2 * e2[2];
}
inline void NucleotideOxdna2::stacking_site() {
  printf("oxDNA2 stacking site\n");
}
template <>
inline void NucleotideOxdna2::base_site<0>() { printf("oxDNA2 base site type 0\n"); }
template <>
inline void NucleotideOxdna2::base_site<1>() { printf("oxDNA2 base site type 1\n"); }
template <>
inline void NucleotideOxdna2::base_site<2>() { printf("oxDNA2 base site type 2\n"); }
template <>
inline void NucleotideOxdna2::base_site<3>() { printf("oxDNA2 base site type 3\n"); }

/* ----------------------------------------------------------------------
   oxDNA3 nucleotide
------------------------------------------------------------------------- */

inline void NucleotideOxdna3::backbone_site() {
  printf("oxDNA3 backbone site\n");
}
inline void NucleotideOxdna3::stacking_site() {
  printf("oxDNA3 stacking site\n");
}
template <>
inline void NucleotideOxdna3::base_site<0>() { printf("oxDNA3 base site type 0\n"); }
template <>
inline void NucleotideOxdna3::base_site<1>() { printf("oxDNA3 base site type 1\n"); }
template <>
inline void NucleotideOxdna3::base_site<2>() { printf("oxDNA3 base site type 2\n"); }
template <>
inline void NucleotideOxdna3::base_site<3>() { printf("oxDNA3 base site type 3\n"); }

/* ----------------------------------------------------------------------
   oxRNA2 nucleotide
------------------------------------------------------------------------- */

inline void NucleotideOxrna2::backbone_site(double e1[3], double /*e2*/[3],
  double e3[3], double rbk[3]) {

  double dx_cbk_oxrna2 = ConstantsOxdna::get_dx_cbk_oxrna2();
  double dz_cbk_oxrna2 = ConstantsOxdna::get_dz_cbk_oxrna2();

  rbk[0] = dx_cbk_oxrna2 * e1[0] + dz_cbk_oxrna2 * e3[0];
  rbk[1] = dx_cbk_oxrna2 * e1[1] + dz_cbk_oxrna2 * e3[1];
  rbk[2] = dx_cbk_oxrna2 * e1[2] + dz_cbk_oxrna2 * e3[2];
}
inline void NucleotideOxrna2::stacking_site() {
  printf("oxRNA2 stacking site\n");
}
template <>
inline void NucleotideOxrna2::base_site<0>() { printf("oxRNA2 base site type 0\n"); }
template <>
inline void NucleotideOxrna2::base_site<1>() { printf("oxRNA2 base site type 1\n"); }
template <>
inline void NucleotideOxrna2::base_site<2>() { printf("oxRNA2 base site type 2\n"); }
template <>
inline void NucleotideOxrna2::base_site<3>() { printf("oxRNA2 base site type 3\n"); }

}
#endif
