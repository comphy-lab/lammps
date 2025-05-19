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
  void backbone_site();
  void stacking_site();
  template <int N>
  void base_site() {};
};

class NucleotideOxdna2 : public NucleotideOxdna<NucleotideOxdna2> {
 public:
  NucleotideOxdna2(){};
  virtual ~NucleotideOxdna2(){};
  void backbone_site();
  void stacking_site();
  template <int N>
  void base_site() {};
};

inline void NucleotideOxdna1::backbone_site() {
  printf("oxDNA1 backbone site\n");
}
inline void NucleotideOxdna1::stacking_site() {
  printf("oxDNA1 stacking site\n");
}


inline void NucleotideOxdna2::backbone_site() {
  printf("oxDNA2 backbone site\n");
}

inline void NucleotideOxdna2::stacking_site() {
  printf("oxDNA2 stacking site\n");
}

template <>
inline void NucleotideOxdna1::base_site<0>() { printf("oxDNA1 base site type 0\n"); }

template <>
inline void NucleotideOxdna1::base_site<1>() { printf("oxDNA1 base site type 1\n"); }


template <>
inline void NucleotideOxdna2::base_site<0>() { printf("oxDNA2 base site type 0\n"); }

template <>
inline void NucleotideOxdna2::base_site<1>() { printf("oxDNA2 base site type 1\n"); }


}
#endif
