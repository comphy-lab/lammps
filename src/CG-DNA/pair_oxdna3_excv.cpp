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
/* ----------------------------------------------------------------------
   Contributing author: Oliver Henrich (University of Strathclyde, Glasgow)
------------------------------------------------------------------------- */

#include "pair_oxdna3_excv.h"

#include "atom.h"
#include "atom_vec_ellipsoid.h"
#include "comm.h"
#include "constants_oxdna.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "mf_oxdna.h"
#include "neigh_list.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>
#include <cassert>

using namespace LAMMPS_NS;
using namespace MFOxdna;

/* ---------------------------------------------------------------
    compute vector COM-hydrogen bonding interaction site in oxDNA3
    A=1, C=2, G=3, T=0
------------------------------------------------------------------ */
template <>
void PairOxdna3Excv::compute_base_site<0>(double e1[3], double /*e2*/[3],
    double /*e3*/[3], double rb[3]) const
{
  double d_cbase = ConstantsOxdna::get_d_cbase();

  rb[0] = d_cbase*e1[0];
  rb[1] = d_cbase*e1[1];
  rb[2] = d_cbase*e1[2];

}

template <>
void PairOxdna3Excv::compute_base_site<1>(double e1[3], double /*e2*/[3],
    double /*e3*/[3], double rb[3]) const
{
  double d_cbase = ConstantsOxdna::get_d_cbase();

  rb[0] = d_cbase*e1[0];
  rb[1] = d_cbase*e1[1];
  rb[2] = d_cbase*e1[2];

}
template <>
void PairOxdna3Excv::compute_base_site<2>(double e1[3], double /*e2*/[3],
    double /*e3*/[3], double rb[3]) const
{
  double d_cbase = ConstantsOxdna::get_d_cbase();

  rb[0] = d_cbase*e1[0];
  rb[1] = d_cbase*e1[1];
  rb[2] = d_cbase*e1[2];

}
template <>
void PairOxdna3Excv::compute_base_site<3>(double e1[3], double /*e2*/[3],
    double /*e3*/[3], double rb[3]) const
{
  double d_cbase = ConstantsOxdna::get_d_cbase();

  rb[0] = d_cbase*e1[0];
  rb[1] = d_cbase*e1[1];
  rb[2] = d_cbase*e1[2];

}

/* ----------------------------------------------------------------------
   set coeffs 
------------------------------------------------------------------------- */

void PairOxdna3Excv::coeff(int narg, char **arg)
{
  int count;

  if (narg != 3 && narg != 11) error->all(FLERR,"Incorrect args for pair coefficients in oxdna3/excv" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi,nlo,nhi;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  assert((ilo == jlo) & (ihi == jhi));
  nlo = ilo;
  nhi = ihi;

  double epsilon_ss_one, sigma_ss_one;
  double cut_ss_ast_one, cut_ss_c_one, b_ss_one;

  double epsilon_sb_one, sigma_sb_one;
  double cut_sb_ast_one, cut_sb_c_one, b_sb_one;

  double epsilon_bb_one, sigma_bb_one;
  double cut_bb_ast_one, cut_bb_c_one, b_bb_one;

  if (narg == 11) {
    // Excluded volume interaction
    // LJ parameters
    epsilon_ss_one = utils::numeric(FLERR,arg[2],false,lmp);
    sigma_ss_one = utils::numeric(FLERR,arg[3],false,lmp);
    cut_ss_ast_one = utils::numeric(FLERR,arg[4],false,lmp);

    // LJ parameters
    epsilon_sb_one = utils::numeric(FLERR,arg[5],false,lmp);
    sigma_sb_one = utils::numeric(FLERR,arg[6],false,lmp);
    cut_sb_ast_one = utils::numeric(FLERR,arg[7],false,lmp);

    // LJ parameters
    epsilon_bb_one = utils::numeric(FLERR,arg[8],false,lmp);
    sigma_bb_one = utils::numeric(FLERR,arg[9],false,lmp);
    cut_bb_ast_one = utils::numeric(FLERR,arg[10],false,lmp);
  }
  else {
    if (comm->me == 0) {
      PotentialFileReader reader(lmp, arg[2], "oxdna3 potential", " (excv)");
      reader.set_bufsize(65336);
      char * line;
      std::string iloc, jloc, potential_name;

      while ((line = reader.next_line())) {
        try {
          ValueTokenizer values(line);
          iloc = values.next_string();
          jloc = values.next_string();
          potential_name = values.next_string();
          if (iloc == arg[0] && jloc == arg[1] && potential_name == "excv") {
            // Excluded volume interaction
            // LJ backbone-backbone parameters
            epsilon_ss_one = values.next_double();
            sigma_ss_one = values.next_double();
            cut_ss_ast_one = values.next_double();

            // LJ backbone-base parameters
            epsilon_sb_one = values.next_double();
            sigma_sb_one = values.next_double();
            cut_sb_ast_one = values.next_double();

            // LJ base-base parameters
            epsilon_bb_one = values.next_double();
            sigma_bb_one = values.next_double();
            cut_bb_ast_one = values.next_double();

            break;
          } else continue;
        } catch (std::exception &e) {
          error->one(FLERR, "Problem parsing oxDNA potential file: {}", e.what());
        }
      }
      if ((iloc != arg[0]) || (jloc != arg[1]) || (potential_name != "excv"))
        error->one(FLERR, "No corresponding excv potential found in file {} for pair type {} {}",
                   arg[2], arg[0], arg[1]);
    }

    MPI_Bcast(&epsilon_ss_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&sigma_ss_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_ss_ast_one, 1, MPI_DOUBLE, 0, world);

    MPI_Bcast(&epsilon_sb_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&sigma_sb_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_sb_ast_one, 1, MPI_DOUBLE, 0, world);

    MPI_Bcast(&epsilon_bb_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&sigma_bb_one, 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&cut_bb_ast_one, 1, MPI_DOUBLE, 0, world);
  }

  // backbone-backbone
  count = 0;

  // smoothing - determined through continuity and differentiability
  b_ss_one = 4.0/sigma_ss_one
      *(6.0*pow(sigma_ss_one/cut_ss_ast_one,7)-12.0*pow(sigma_ss_one/cut_ss_ast_one,13))
      *4.0/sigma_ss_one*(6.0*pow(sigma_ss_one/cut_ss_ast_one,7)-12.0*pow(sigma_ss_one/cut_ss_ast_one,13))
      /4.0/(4.0*(pow(sigma_ss_one/cut_ss_ast_one,12)-pow(sigma_ss_one/cut_ss_ast_one,6)));

  cut_ss_c_one = cut_ss_ast_one
      - 2.0*4.0*(pow(sigma_ss_one/cut_ss_ast_one,12)-pow(sigma_ss_one/cut_ss_ast_one,6))
      /(4.0/sigma_ss_one*(6.0*pow(sigma_ss_one/cut_ss_ast_one,7)-12.0*pow(sigma_ss_one/cut_ss_ast_one,13)));

  // backbone-backbone parameters depending on base step
  for (int i = nlo; i <= nhi; i++) {
    for (int j = nlo; j <= nhi; j++) {
      epsilon_ss[i][j] = epsilon_ss_one;
      sigma_ss[i][j] = sigma_ss_one;
      cut_ss_ast[i][j] = cut_ss_ast_one;
      b_ss[i][j] = b_ss_one;
      cut_ss_c[i][j] = cut_ss_c_one;
      lj1_ss[i][j] = 4.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],12.0);
      lj2_ss[i][j] = 4.0 * epsilon_ss[i][j] * pow(sigma_ss[i][j],6.0);
      cutsq_ss_ast[i][j] = cut_ss_ast[i][j]*cut_ss_ast[i][j];
      cutsq_ss_c[i][j]  = cut_ss_c[i][j]*cut_ss_c[i][j];
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv" + utils::errorurl(21));

  // backbone-base
  count = 0;

  // smoothing - determined through continuity and differentiability
  b_sb_one = 4.0/sigma_sb_one
      *(6.0*pow(sigma_sb_one/cut_sb_ast_one,7)-12.0*pow(sigma_sb_one/cut_sb_ast_one,13))
      *4.0/sigma_sb_one*(6.0*pow(sigma_sb_one/cut_sb_ast_one,7)-12.0*pow(sigma_sb_one/cut_sb_ast_one,13))
      /4.0/(4.0*(pow(sigma_sb_one/cut_sb_ast_one,12)-pow(sigma_sb_one/cut_sb_ast_one,6)));

  cut_sb_c_one = cut_sb_ast_one
      - 2.0*4.0*(pow(sigma_sb_one/cut_sb_ast_one,12)-pow(sigma_sb_one/cut_sb_ast_one,6))
      /(4.0/sigma_sb_one*(6.0*pow(sigma_sb_one/cut_sb_ast_one,7)-12.0*pow(sigma_sb_one/cut_sb_ast_one,13)));

  // backbone-base parameters depending on base step
  for (int i = nlo; i <= nhi; i++) {
    for (int j = nlo; j <= nhi; j++) {
      epsilon_sb[i][j] = epsilon_sb_one;
      sigma_sb[i][j] = sigma_sb_one;
      cut_sb_ast[i][j] = cut_sb_ast_one;
      b_sb[i][j] = b_sb_one;
      cut_sb_c[i][j] = cut_sb_c_one;
      lj1_sb[i][j] = 4.0 * epsilon_sb[i][j] * pow(sigma_sb[i][j],12.0);
      lj2_sb[i][j] = 4.0 * epsilon_sb[i][j] * pow(sigma_sb[i][j],6.0);
      cutsq_sb_ast[i][j] = cut_sb_ast[i][j]*cut_sb_ast[i][j];
      cutsq_sb_c[i][j]  = cut_sb_c[i][j]*cut_sb_c[i][j];
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv" + utils::errorurl(21));

  // base-base
  count = 0;

  // smoothing - determined through continuity and differentiability
  b_bb_one = 4.0/sigma_bb_one
      *(6.0*pow(sigma_bb_one/cut_bb_ast_one,7)-12.0*pow(sigma_bb_one/cut_bb_ast_one,13))
      *4.0/sigma_bb_one*(6.0*pow(sigma_bb_one/cut_bb_ast_one,7)-12.0*pow(sigma_bb_one/cut_bb_ast_one,13))
      /4.0/(4.0*(pow(sigma_bb_one/cut_bb_ast_one,12)-pow(sigma_bb_one/cut_bb_ast_one,6)));

  cut_bb_c_one = cut_bb_ast_one
      - 2.0*4.0*(pow(sigma_bb_one/cut_bb_ast_one,12)-pow(sigma_bb_one/cut_bb_ast_one,6))
      /(4.0/sigma_bb_one*(6.0*pow(sigma_bb_one/cut_bb_ast_one,7)-12.0*pow(sigma_bb_one/cut_bb_ast_one,13)));

  // base-base parameters depending on base step
  for (int i = nlo; i <= nhi; i++) {
    for (int j = nlo; j <= nhi; j++) {
      epsilon_bb[i][j] = epsilon_bb_one;
      sigma_bb[i][j] = sigma_bb_one;
      cut_bb_ast[i][j] = cut_bb_ast_one;
      b_bb[i][j] = b_bb_one;
      cut_bb_c[i][j] = cut_bb_c_one;
      lj1_bb[i][j] = 4.0 * epsilon_bb[i][j] * pow(sigma_bb[i][j],12.0);
      lj2_bb[i][j] = 4.0 * epsilon_bb[i][j] * pow(sigma_bb[i][j],6.0);
      cutsq_bb_ast[i][j] = cut_bb_ast[i][j]*cut_bb_ast[i][j];
      cutsq_bb_c[i][j]  = cut_bb_c[i][j]*cut_bb_c[i][j];
      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv");

  // base-base parameters depending on tetramer
  count = 0;

  for (int i = 0; i <= nhi; i++) { // type 0 for terminal j
    for (int j = nlo; j <= nhi; j++) {
      for (int k = nlo; k <= nhi; k++) {
        for (int l = 0; l <= nhi; l++) { // type 0 for terminal k
          sigma4_bb[i][j][k][l] = sigma_bb_one;
          cut4_bb_ast[i][j][k][l] = cut_bb_ast_one;
          b4_bb[i][j][k][l] = b_bb_one;
          cut4_bb_c[i][j][k][l] = cut_bb_c_one;
          cut4sq_bb_ast[i][j][k][l] = cut4_bb_ast[i][j][k][l]*cut4_bb_ast[i][j][k][l];
          cut4sq_bb_c[i][j][k][l]  = cut4_bb_c[i][j][k][l]*cut4_bb_c[i][j][k][l];
          lj14_bb[i][j][k][l] = 4.0 * epsilon_bb[j][k] * pow(sigma4_bb[i][j][k][l],12.0);
          lj24_bb[i][j][k][l] = 4.0 * epsilon_bb[j][k] * pow(sigma4_bb[i][j][k][l],6.0);
          count++;
       }
      }
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/excv");

}
