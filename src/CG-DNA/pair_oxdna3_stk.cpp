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

#include "pair_oxdna3_stk.h"

#include "atom.h"
#include "comm.h"
#include "constants_oxdna.h"
#include "error.h"
#include "force.h"
#include "math_extra.h"
#include "memory.h"
#include "mf_oxdna.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "potential_file_reader.h"

#include <cmath>
#include <cstring>
#include <cassert>

using namespace LAMMPS_NS;
using namespace MFOxdna;

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairOxdna3Stk::coeff(int narg, char **arg)
{
  int count;

  if (narg != 7) error->all(FLERR,"Incorrect args for pair coefficients in oxdna3/stk, use potential file" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo,ihi,jlo,jhi,nlo,nhi,jmod4,kmod4;
  utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
  utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

  assert((ilo == jlo) & (ihi == jhi));
  nlo = ilo;
  nhi = ihi;

  // stacking interaction
  count = 0;

  double T, epsilon_st_one, xi_st_one, kappa_st_one, a_st_one, b_st_lo_one, b_st_hi_one;
  double cut_st_0_one, cut_st_c_one, cut_st_lo_one, cut_st_hi_one;
  double cut_st_lc_one, cut_st_hc_one, tmp, shift_st_one;

  double a_st4_one, theta_st4_0_one, dtheta_st4_ast_one;
  double b_st4_one, dtheta_st4_c_one;

  double a_st5_one, theta_st5_0_one, dtheta_st5_ast_one;
  double b_st5_one, dtheta_st5_c_one;

  double a_st6_one, theta_st6_0_one, dtheta_st6_ast_one;
  double b_st6_one, dtheta_st6_c_one;

  double a_st1_one, cosphi_st1_ast_one, b_st1_one, cosphi_st1_c_one;
  double a_st2_one, cosphi_st2_ast_one, b_st2_one, cosphi_st2_c_one;

  if (strcmp(arg[2], "seqav") != 0 && strcmp(arg[2], "seqdep") != 0) {
    error->all(FLERR,"Incorrect setting, select seqav or seqdep in oxdna/stk");
  }
  if (strcmp(arg[2],"seqav")  == 0) seqdepflag = 0;
  if (strcmp(arg[2],"seqdep") == 0) seqdepflag = 1;

  T = utils::numeric(FLERR,arg[3],false,lmp);
  xi_st_one = utils::numeric(FLERR,arg[4],false,lmp);
  kappa_st_one = utils::numeric(FLERR,arg[5],false,lmp);

  if (comm->me == 0) {
    PotentialFileReader reader(lmp, arg[6], "oxdna potential", " (stk)");
    reader.set_bufsize(65336);
    char * line;
    std::string iloc, jloc, potential_name;

    while ((line = reader.next_line())) {
      try {
        ValueTokenizer values(line);
        iloc = values.next_string();
        jloc = values.next_string();
        potential_name = values.next_string();
        if (iloc == arg[0] && jloc == arg[1] && potential_name == "stk") {

          a_st_one = values.next_double();
          cut_st_0_one = values.next_double();
          cut_st_c_one = values.next_double();
          cut_st_lo_one = values.next_double();
          cut_st_hi_one = values.next_double();

          a_st4_one = values.next_double();
          theta_st4_0_one = values.next_double();
          dtheta_st4_ast_one = values.next_double();
          a_st5_one = values.next_double();
          theta_st5_0_one = values.next_double();
          dtheta_st5_ast_one = values.next_double();
          a_st6_one = values.next_double();
          theta_st6_0_one = values.next_double();
          dtheta_st6_ast_one = values.next_double();
          a_st1_one = values.next_double();
          cosphi_st1_ast_one = values.next_double();
          a_st2_one = values.next_double();
          cosphi_st2_ast_one = values.next_double();

          break;
        } else continue;
      } catch (std::exception &e) {
        error->one(FLERR, "Problem parsing oxDNA potential file: {}", e.what());
      }
    }
    if ((iloc != arg[0]) || (jloc != arg[1]) || (potential_name != "stk"))
      error->one(FLERR, "No corresponding stk potential found in file {} for pair type {} {}",
                 arg[4], arg[0], arg[1]);

    epsilon_st_one = stacking_strength(xi_st_one, kappa_st_one, T);

  }

  MPI_Bcast(&epsilon_st_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&a_st_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_st_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_st_c_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_st_lo_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cut_st_hi_one, 1, MPI_DOUBLE, 0, world);

  MPI_Bcast(&a_st4_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_st4_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_st4_ast_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&a_st5_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_st5_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_st5_ast_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&a_st6_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&theta_st6_0_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&dtheta_st6_ast_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&a_st1_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cosphi_st1_ast_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&a_st2_one, 1, MPI_DOUBLE, 0, world);
  MPI_Bcast(&cosphi_st2_ast_one, 1, MPI_DOUBLE, 0, world);


  b_st_lo_one = 2*a_st_one*exp(-a_st_one*(cut_st_lo_one-cut_st_0_one))*
        2*a_st_one*exp(-a_st_one*(cut_st_lo_one-cut_st_0_one))*
        (1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))/
        (4*((1-exp(-a_st_one*(cut_st_lo_one -cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))-
        (1-exp(-a_st_one*(cut_st_c_one -cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_c_one-cut_st_0_one)))));

  cut_st_lc_one = cut_st_lo_one - a_st_one*exp(-a_st_one*(cut_st_lo_one-cut_st_0_one))*
        (1-exp(-a_st_one*(cut_st_lo_one-cut_st_0_one)))/b_st_lo_one;

  b_st_hi_one = 2*a_st_one*exp(-a_st_one*(cut_st_hi_one-cut_st_0_one))*
        2*a_st_one*exp(-a_st_one*(cut_st_hi_one-cut_st_0_one))*
        (1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))/
        (4*((1-exp(-a_st_one*(cut_st_hi_one -cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))-
        (1-exp(-a_st_one*(cut_st_c_one -cut_st_0_one)))*
        (1-exp(-a_st_one*(cut_st_c_one-cut_st_0_one)))));

  cut_st_hc_one = cut_st_hi_one - a_st_one*exp(-a_st_one*(cut_st_hi_one-cut_st_0_one))*
        (1-exp(-a_st_one*(cut_st_hi_one-cut_st_0_one)))/b_st_hi_one;

  tmp = 1 - exp(-(cut_st_c_one-cut_st_0_one) * a_st_one);
  shift_st_one = epsilon_st_one * tmp * tmp;

  b_st4_one = a_st4_one*a_st4_one*dtheta_st4_ast_one*dtheta_st4_ast_one/(1-a_st4_one*dtheta_st4_ast_one*dtheta_st4_ast_one);
  dtheta_st4_c_one = 1/(a_st4_one*dtheta_st4_ast_one);

  b_st5_one = a_st5_one*a_st5_one*dtheta_st5_ast_one*dtheta_st5_ast_one/(1-a_st5_one*dtheta_st5_ast_one*dtheta_st5_ast_one);
  dtheta_st5_c_one = 1/(a_st5_one*dtheta_st5_ast_one);

  b_st6_one = a_st6_one*a_st6_one*dtheta_st6_ast_one*dtheta_st6_ast_one/(1-a_st6_one*dtheta_st6_ast_one*dtheta_st6_ast_one);
  dtheta_st6_c_one = 1/(a_st6_one*dtheta_st6_ast_one);

  b_st1_one = a_st1_one*a_st1_one*cosphi_st1_ast_one*cosphi_st1_ast_one/(1-a_st1_one*cosphi_st1_ast_one*cosphi_st1_ast_one);
  cosphi_st1_c_one = 1/(a_st1_one*cosphi_st1_ast_one);

  b_st2_one = a_st2_one*a_st2_one*cosphi_st2_ast_one*cosphi_st2_ast_one/(1-a_st2_one*cosphi_st2_ast_one*cosphi_st2_ast_one);
  cosphi_st2_c_one = 1/(a_st2_one*cosphi_st2_ast_one);

  for (int i = 0; i <= nhi; i++) { // type 0 for terminal j
    for (int j = nlo; j <= nhi; j++) {
      jmod4 = j%4;
      if (jmod4 == 0) jmod4 = 4;

      for (int k = nlo; k <= nhi; k++) {
        kmod4 = k%4;
        if (kmod4 == 0) kmod4 = 4;

        epsilon_st[j][k] = epsilon_st_one;
        if (seqdepflag) epsilon_st[j][k] *= eta_st[jmod4-1][kmod4-1];
        a_st[j][k] = a_st_one;

      // parameters depending on tetramer
        for (int l = 0; l <= nhi; l++) { // type 0 for terminal k
          cut_st_0[i][j][k][l] = cut_st_0_one;
          cut_st_c[i][j][k][l] = cut_st_c_one;
          cut_st_lo[i][j][k][l] = cut_st_lo_one;
          cut_st_hi[i][j][k][l] = cut_st_hi_one;
          cut_st_lc[i][j][k][l] = cut_st_lc_one;
          cut_st_hc[i][j][k][l] = cut_st_hc_one;
          cutsq_st_hc[i][j][k][l] = cut_st_hc[i][j][k][l]*cut_st_hc[i][j][k][l];
          b_st_lo[i][j][k][l] = b_st_lo_one;
          b_st_hi[i][j][k][l] = b_st_hi_one;
          shift_st[i][j][k][l] = shift_st_one;
          if (seqdepflag) {
            shift_st[i][j][k][l] *= eta_st[jmod4-1][kmod4-1];
          }

          a_st4[i][j][k][l] = a_st4_one;
          dtheta_st4_ast[i][j][k][l] = dtheta_st4_ast_one;
          b_st4[i][j][k][l] = b_st4_one;
          dtheta_st4_c[i][j][k][l] = dtheta_st4_c_one;
        }
      }

      theta_st4_0[i][j] = theta_st4_0_one;

      a_st5[i][j] = a_st5_one;
      theta_st5_0[i][j] = theta_st5_0_one;
      dtheta_st5_ast[i][j] = dtheta_st5_ast_one;
      b_st5[i][j] = b_st5_one;
      dtheta_st5_c[i][j] = dtheta_st5_c_one;

      a_st6[i][j] = a_st6_one;
      theta_st6_0[i][j] = theta_st6_0_one;
      dtheta_st6_ast[i][j] = dtheta_st6_ast_one;
      b_st6[i][j] = b_st6_one;
      dtheta_st6_c[i][j] = dtheta_st6_c_one;

      a_st1[i][j] = a_st1_one;
      cosphi_st1_ast[i][j] = cosphi_st1_ast_one;
      b_st1[i][j] = b_st1_one;
      cosphi_st1_c[i][j] = cosphi_st1_c_one;

      a_st2[i][j] = a_st2_one;
      cosphi_st2_ast[i][j] = cosphi_st2_ast_one;
      b_st2[i][j] = b_st2_one;
      cosphi_st2_c[i][j] = cosphi_st2_c_one;

      setflag[i][j] = 1;
      count++;
    }
  }

  if (count == 0) error->all(FLERR,"Incorrect args for pair coefficients in oxdna/stk" + utils::errorurl(21));

}

