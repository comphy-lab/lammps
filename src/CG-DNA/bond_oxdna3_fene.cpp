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

#include "bond_oxdna3_fene.h"

#include "atom.h"
#include "comm.h"
#include "error.h"
#include "potential_file_reader.h"

using namespace LAMMPS_NS;

/* ----------------------------------------------------------------------
   set coeffs for one type
------------------------------------------------------------------------- */

void BondOxdna3Fene::coeff(int narg, char **arg)
{
  if (narg != 2 && narg != 4) error->all(FLERR, "Incorrect args for bond coefficients in oxdna/fene" + utils::errorurl(21));
  if (!allocated) allocate();

  int ilo, ihi;
  utils::bounds(FLERR, arg[0], 1, atom->nbondtypes, ilo, ihi, error);

  if (narg == 4) {
    k[ilo] = utils::numeric(FLERR, arg[1], false, lmp);
    Delta[ilo][0][0][0][0] = utils::numeric(FLERR, arg[2], false, lmp);
    r0[ilo][0][0][0][0] = utils::numeric(FLERR, arg[3], false, lmp);
  } else {
    if (comm->me == 0) { // read values from potential file
      PotentialFileReader reader(lmp, arg[1], "oxdna potential", " (fene)");
      char * line;
      std::string iloc, potential_name;

      while ((line = reader.next_line())) {
        try {
          ValueTokenizer values(line);
          iloc = values.next_string();
          potential_name = values.next_string();
          if (iloc == arg[0] && potential_name == "fene") {
            k[ilo] = values.next_double();
            Delta[ilo][0][0][0][0] = values.next_double();
            r0[ilo][0][0][0][0] = values.next_double();

            break;
          } else continue;
        } catch (std::exception &e) {
          error->one(FLERR, "Problem parsing oxDNA potential file: {}", e.what());
        }
      }
      if ((iloc != arg[0]) || (potential_name != "fene"))
        error->one(FLERR, "No corresponding fene potential found in file {} for bond type {}",
                   arg[1], arg[0]);
    }

    MPI_Bcast(&k[ilo], 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&Delta[ilo][0][0][0][0], 1, MPI_DOUBLE, 0, world);
    MPI_Bcast(&r0[ilo][0][0][0][0], 1, MPI_DOUBLE, 0, world);
  }

  int count = 0;
  int n = atom->ntypes;

  for (int i = ilo; i <= ihi; i++) {
    k[i] = k[ilo];
    for (int n1 = 0; n1 <= n; n1++) {
      for (int n2 = 0; n2 <= n; n2++) {
        for (int n3 = 0; n3 <= n; n3++) {
          for (int n4 = 0; n4 <= n; n4++) {
            Delta[i][n1][n2][n3][n4] = Delta[ilo][0][0][0][0];
            r0[i][n1][n2][n3][n4] = r0[ilo][0][0][0][0];
          }
        }
      }
    }
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR, "Incorrect args for bond coefficients in oxdna/fene" + utils::errorurl(21));
}
