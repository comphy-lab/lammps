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
  virtual void compute(int, int);
  void settings(int, char **);
  void coeff(int, char **);
  void init_style();
  double init_one(int, int);
  void write_restart_settings(FILE *);
  void read_restart_settings(FILE *);
  double single(int, int, int, int, double, double, double, double &);

  /* I do this the same way it's done in force.h */

  char *indicator_style;
  class LddIndicator *indicator;
  typedef LddIndicator *(*IndicatorCreator)(LAMMPS *);
  typedef std::map<std::string,IndicatorCreator> IndicatorCreatorMap;
  IndicatorCreatorMap *indicator_map;
  class LddIndicator * new_indicator(std::string wtype);

  char *potential_style;
  class LddPotential *potential;
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
  double cut_LDD_global;
  double **cut;
  bool **self_interaction; 
  bool **ignore_pair;
  bool *ignore_me; 
  bool **bGradient; 

  LddIndicator ***Inds; 
  LddPotential ***Potls; 
  LddPotential ***GradPotls; 
 
  virtual void allocate();
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

/* ERROR/WARNING messages: 
E: Illegal ... command

Self-Explanatory. Check the input script syntax and compare to the documentationf or the command. You can use -echo screen as a command line option when running LAMMPS to see the offending line

E: Incorrect args for pair coefficients

Self-Explanatory. Check the input script of data file.

E: ??

??
*/
