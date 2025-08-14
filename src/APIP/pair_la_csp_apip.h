/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef PAIR_CLASS
// clang-format off
PairStyle(la/csp/apip,PairLACSPAPIP);
// clang-format on
#else

#ifndef LMP_PAIR_LA_CSP_APIP_H
#define LMP_PAIR_LA_CSP_APIP_H

#include "fix_lambda_la_csp_apip.h"
#include "pair.h"

namespace LAMMPS_NS {

class PairLACSPAPIP : public Pair {
 public:
  PairLACSPAPIP(class LAMMPS *);
  ~PairLACSPAPIP() override;

  void compute(int, int) override;
  void settings(int, char **) override;
  void coeff(int, char **) override;
  void init_style() override;
  double init_one(int, int) override;
  int pack_forward_comm(int, int *, double *, int, int *) override;
  void unpack_forward_comm(int, int, double *) override;

  void *extract(const char *, int &) override;

 protected:

  FixLambdaLACSPAPIP * fix_la_csp;

  double cutsq_combined;
  double fix_la_csp_w_cutsq;
  int nnn_half;

  double * prefactor1; // per atom array
  int prefactor1_size; // own + ghosts
  double * prefactor2; // per atom array
  int prefactor2_size; // own


  virtual void allocate();
  double get_cutsq_combined();

  double **scale;

  void calculate_time_per_atom();

  int n_computations_accumulated;    // number of accumulated computations
  double time_wall_accumulated;      // accumulated compute time
  double time_per_atom;              // average time of one computation
};
}    // namespace LAMMPS_NS

#endif
#endif
