/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   BAOAB Langevin integrator (Leimkuhler & Matthews, 2013)

   Implements the BAOAB splitting of the Langevin equation:
     B: half-step velocity kick from conservative forces
     A: half-step position drift
     O: full-step exact Ornstein-Uhlenbeck thermostat
     A: half-step position drift
     B: half-step velocity kick from conservative forces

   This fix performs COMPLETE time integration — do NOT pair with fix nve.

   Usage:
     fix ID group-ID baoab Tstart Tstop damp seed [keyword value ...]

   Optional keywords:
     zero yes/no  — zero total random momentum (default: no)

   Reference:
     B. Leimkuhler and C. Matthews, "Rational Construction of Stochastic
     Numerical Methods for Molecular Sampling", Appl. Math. Res. Express,
     2013(1), 34-56 (2013). https://doi.org/10.1093/amrx/abs010
------------------------------------------------------------------------- */

#ifdef FIX_CLASS
// clang-format off
FixStyle(baoab, FixBAOAB)
// clang-format on
#else

#ifndef LMP_FIX_BAOAB_H
#define LMP_FIX_BAOAB_H

#include "fix.h"

namespace LAMMPS_NS {

class FixBAOAB : public Fix {
 public:
  FixBAOAB(class LAMMPS *, int, char **);
  ~FixBAOAB() override;

  int setmask() override;
  void init() override;
  void setup(int) override;
  void initial_integrate(int) override;
  void final_integrate() override;
  void reset_dt() override;
  void write_restart(FILE *) override;
  void restart(char *) override;
  double compute_scalar() override;
  int modify_param(int, char **) override;

 private:
  double t_start, t_stop;      // temperature ramp endpoints
  double t_target;             // current target temperature
  double gamma;                // friction coefficient = 1/damp
  double damp;                 // damping time (user-specified)
  int seed;                    // RNG seed

  // Precomputed timestep-dependent quantities
  double dtf;                  // 0.5*dt*ftm2v  (force->velocity kick)
  double dtby2;                // 0.5*dt        (half-step drift)
  double c1;                   // exp(-gamma*dt) for O step

  // Energy accounting
  int tallyflag;               // whether to tally thermostat energy
  double energy;               // cumulative thermostat energy
  double energy_onestep;       // thermostat energy this step

  // Options
  int zeroflag;                // zero total random momentum?

  class RanMars *random;

  void compute_target();       // update t_target from ramp
  void recompute_dt_coeffs();  // recompute dtf, dtby2, c1
};

}    // namespace LAMMPS_NS
#endif
#endif
