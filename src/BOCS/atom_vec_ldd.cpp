// clang-format off
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

#include "atom_vec_ldd.h"

#include "atom.h"
#include "error.h"
#include "fix.h"
#include "fix_adapt.h"
#include "math_const.h"
#include "modify.h"

#include <cstring>

using namespace LAMMPS_NS;
using namespace MathConst;

/* ---------------------------------------------------------------------- */
namespace {
// With this structure I have now imposed a max of <100k dif ldd atom types that can be used in compute_property_atom
// I hope this is practically beyond the limits of what one might want to simulate.
  enum Ldd_compute_types : int {
          LDD_LOCAL_DENSITY = 100000, // Planning for back 5 digits to tell me the ldd type
          LDD_ENERGY = 200000, // while the 10^6 place tells me the kind of data
          LDD_GRAD_DENSITYX = 300000,
          LDD_GRAD_DENSITYY = 400000,
          LDD_GRAD_DENSITYZ = 500000,
          LDD_GRAD_ENERGY=  600000,
          LDD_TOTAL_ENERGY = 700000}; // Here I'll only take 700k since its not a per_atom_vec

} // End anon NS - Added to hide constants for compute_property_atom / property_atom recognition
AtomVecLdd::AtomVecLdd(LAMMPS *lmp) : AtomVec(lmp)
{
  mass_type = PER_TYPE;
  molecular = Atom::ATOMIC;

  atom->ldd_big_flag = 1;

  // strings with peratom variables to include in each AtomVec method
  // strings cannot contain fields in corresponding AtomVec default strings
  // order of fields in a string does not matter
  // except: fields_data_atom & fields_data_vel must match data file

  fields_grow = {"ldd_local_density", "ldd_energy", "ldd_grad_density", "ldd_grad_energy", "ldd_total_energy"};
  fields_copy = {"ldd_local_density", "ldd_energy", "ldd_grad_density", "ldd_grad_energy", "ldd_total_energy"};
  fields_comm = {};
  fields_comm_vel = {};
  fields_reverse = {};
  fields_border = {};
  fields_border_vel = {};
  fields_exchange = {};
  fields_restart = {};
  fields_create = {"ldd_local_density", "ldd_energy", "ldd_grad_density", "ldd_grad_energy", "ldd_total_energy"};

  fields_data_atom = {"id", "type", "x"};
  fields_data_vel = {"id", "v"};

//  setup_fields();
}

void AtomVecLdd::process_args(int narg, char **arg)
{
  if (narg != 1)
    error->all(FLERR,"Illegal atom_style ldd command");

  atom->ldd_ntypes = utils::inumeric(FLERR,arg[0],true,lmp);

  int n = atom->ldd_ntypes;
  atom->add_peratom_change_columns("ldd_local_density",n+1);
  atom->add_peratom_change_columns("ldd_energy",n+1);
  atom->add_peratom_change_columns("ldd_grad_density",3*(n+1));
  atom->add_peratom_change_columns("ldd_grad_energy",n+1);


  // delay setting up of fields until now
  setup_fields();
}


/* ----------------------------------------------------------------------
   set local copies of all grow ptrs used by this class, except defaults
   needed in replicate when 2 atom classes exist and it calls pack_restart()
------------------------------------------------------------------------- */

void AtomVecLdd::grow_pointers()
{
  ldd_local_density = atom->ldd_local_density;
  ldd_energy = atom->ldd_energy;
  ldd_grad_density = atom->ldd_grad_density;
  ldd_grad_energy = atom->ldd_grad_energy;
  ldd_total_energy = atom->ldd_total_energy;
}


/* Need to make properties available for computes */

/* property_atom gives us index for packing atom style later
 *  The index returned is a big number (see enum at top) so that we can use
 *  certain digits to extract column info in our per_atom vectors*/
int AtomVecLdd::property_atom(const std::string &name)
{

 int typecol; // desired per atom vector column
 int compute_type; // stores the determined Ldd_compute_type
 std::string key; // Holds the root keyword without column wanted

 // first I determine the ldd member wanted
 if (utils::strmatch(name, "ldd_total_energy"))
         return LDD_TOTAL_ENERGY; // scalar per atom, so no added work
 else if (utils::strmatch(name, "ldd_local_density")){
         key = "ldd_local_density"; // for vectors we find colidx wanted
         compute_type = Ldd_compute_types::LDD_LOCAL_DENSITY;
 }
 else if (utils::strmatch(name, "ldd_energy")){
         key = "ldd_energy";
         compute_type = Ldd_compute_types::LDD_ENERGY;
 }
 else if (utils::strmatch(name, "ldd_grad_densityx")){
         key = "ldd_grad_densityx";
         compute_type = Ldd_compute_types::LDD_GRAD_DENSITYX;
 }
 else if (utils::strmatch(name, "ldd_grad_densityy")){
         key = "ldd_grad_densityy";
         compute_type = Ldd_compute_types::LDD_GRAD_DENSITYY;
 }
 else if (utils::strmatch(name, "ldd_grad_densityz")){
         key = "ldd_grad_densityz";
         compute_type = Ldd_compute_types::LDD_GRAD_DENSITYZ;
 }
 else if (utils::strmatch(name, "ldd_gradnrg")){
         key = "ldd_gradnrg";
         compute_type = Ldd_compute_types::LDD_GRAD_ENERGY;
 }
 else {return -1;}

 // Then if its a vector we get the column wanted and pack it into the index
 if (utils::strmatch(name, key)){
         std::string lddid = utils::trim(name.substr(key.size()));
         try {
                 typecol = utils::inumeric(FLERR, lddid, true, lmp);
         }
         catch (...) {
                 error->one(FLERR, "ldd type id {} not valid in property/atom", lddid);
         }
         if (typecol > atom->ldd_ntypes || typecol < 1) {
                 error->one(FLERR, "ldd type id {} not valid in property/atom", typecol);
         }
   return compute_type + typecol;
 }
 return -1;
}

/* pack_property_atom infers the ldd type to pack/(if applicable)the column to pack based on the big
 * index returned by pack_per_atom */
void AtomVecLdd::pack_property_atom(int index, double *buf, int nvalues, int groupbit)
{
  int *mask =  atom->mask;
  int nlocal = atom->nlocal;
  int n = 0;
  double **ldd_vec_variable; // all of them are n_atom x n_type or n_atom x 3 n_type
  int ldd_type_col;
  int ldd_gradxyz_col;

  // Check variable types by greatest to least constant vals
  if (index == LDD_TOTAL_ENERGY) {
          double* total_energy = ldd_total_energy;
          for (int i = 0; i < nlocal; i++){
                  if (mask[i] & groupbit) {
                          buf[n] = total_energy[i];}
                  else {
                          buf[n] = 0.0;}
                  n += nvalues;
          }
          return; // done unpacking in this case
  }
  else if (index > LDD_GRAD_ENERGY) {
          ldd_vec_variable = ldd_grad_energy;
          ldd_type_col = index - LDD_GRAD_ENERGY;
  }
  else if (index > LDD_GRAD_DENSITYZ){
          ldd_vec_variable = ldd_grad_density;
          ldd_type_col = 3*(index - LDD_GRAD_DENSITYZ) + 2;
  }
  else if (index > LDD_GRAD_DENSITYY){
          ldd_vec_variable = ldd_grad_density;
          ldd_type_col = 3*(index - LDD_GRAD_DENSITYY) + 1;
  }
  else if (index > LDD_GRAD_DENSITYX) {
          ldd_vec_variable = ldd_grad_density;
          ldd_type_col = 3*(index - LDD_GRAD_DENSITYX);
  }
  else if (index > LDD_ENERGY){
          ldd_vec_variable = ldd_energy;
          ldd_type_col = index - LDD_ENERGY;
  }
  else if (index > LDD_LOCAL_DENSITY){
          ldd_vec_variable = ldd_local_density;
          ldd_type_col = index - LDD_LOCAL_DENSITY;
  }

  // if it was one of the vectors unpack the column info req.
  for (int i = 0; i < nlocal; i++) {
          if (mask[i] & groupbit){
                  buf[n] = ldd_vec_variable[i][ldd_type_col];}
          else{
                  buf[n] = 0.0;}
            n+=nvalues;
          }
}
