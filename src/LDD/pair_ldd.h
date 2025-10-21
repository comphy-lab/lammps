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
/* ------------------------------------------------------
    This file is part of the LDD package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    and Maria C. Lesniewski, mjl6766@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(ldd,PairLdd);
// clang-format on
#else

#ifndef LMP_PAIR_LDD_H
#define LMP_PAIR_LDD_H

#include <map>

#include "pair.h"


namespace LAMMPS_NS {

//Forward Declarations
class LddIndicator;
class LddPotential;

class PairLdd: public Pair {
 public:
  PairLdd(class LAMMPS *);
  virtual ~PairLdd();
  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  void write_restart_settings(FILE *) override;
  void read_restart_settings(FILE *) override;
  double single(int, int, int, int, double, double, double, double &) override;

  /* I do this the same way it's done in force.h */

  char *indicator_style; // key for LDD indicator map
  class LddIndicator *indicator; // w(r) and all associated info
  typedef LddIndicator *(*IndicatorCreator)(LAMMPS *);
  typedef std::map<std::string,IndicatorCreator> IndicatorCreatorMap;
  IndicatorCreatorMap *indicator_map;
  class LddIndicator * new_indicator(std::string wtype);

  char *potential_style; // key for LDD potential U_x map
  class LddPotential *potential; // U_x and all associated info
  typedef LddPotential *(*PotentialCreator)(LAMMPS *);
  typedef std::map<std::string,PotentialCreator> PotentialCreatorMap;
  PotentialCreatorMap *potential_map;
  class LddPotential * new_potential (std::string ptype);

  // Functions to calculate the local densities, the gradients,
  // and the associated energies
  void LDD_calculate_LDs();
  void LDD_calculate_energies();

  int pack_forward_comm(int, int *, double *, int, int *);
  void unpack_forward_comm(int, int, double *);
  int pack_reverse_comm(int, int, double *);
  void unpack_reverse_comm(int, int *, double *);

 protected:
  double cut_LDD_global; // longest interaction cutoff obtained from pair_style ldd
  double **cut; // cutoffs passed in pair_coeff ldd commands. dim n_types x n_types
  bool **self_interaction; // Self interaction settings from pair_coeff. dim n_types x n_types
  bool **ignore_pair; // Ignored Potential list. dim n_types x n_types
  bool *ignore_me;  // Totally ignored central type list. dim_ntypes
  bool **bGradient; // List of pair types a surrounded by b that also have a gradient interaction. dim n_types x n_types

  LddIndicator ***Inds; //The address of an n_type x n_type structure, holding all a|b indicator info
  LddPotential ***Potls; //The address of an n_type x n_type structure, holding all a|b U_rho info
  LddPotential ***GradPotls; //The address of an n_type x n_type structure, holding all a|b U_{\nabla} info

  void allocate();
  void ErrorDoubleKeyword(const char *);
  void ErrorNumKeywordArgs(const char *, const char *);

  private:
  // Again, the same as done in force.h
    void LDD_factory();
    template <typename T> static LddIndicator *indicator_creator(LAMMPS *);
    template <typename T> static LddPotential *potential_creator(LAMMPS *);

};

}

#endif
#endif
