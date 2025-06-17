// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributead under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */
/* ------------------------------------------------------
   This file is part of the USER-LDD package for LAMMPS.
   Contributed by Michael R. DeLyser, mrd5285@psu.edu
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

