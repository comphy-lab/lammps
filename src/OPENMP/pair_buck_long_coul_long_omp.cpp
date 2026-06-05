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

#include "pair_buck_long_coul_long_omp.h"

#include "atom.h"
#include "comm.h"
#include "ewald_const.h"
#include "force.h"
#include "math_extra.h"
#include "neigh_list.h"
#include "suffix.h"

#include <cmath>
#include <cstring>

#include "omp_compat.h"
using namespace LAMMPS_NS;
using namespace MathExtra;
using namespace EwaldConst;

// bits in ewald_order / ewald_off flag which interaction order is treated
// with the long-range (k-space) solver: order 1 = Coulomb (1/r),
// order 6 = dispersion (1/r^6).

enum { EWALD_COUL = 1 << 1, EWALD_DISP = 1 << 6 };

/* ---------------------------------------------------------------------- */

PairBuckLongCoulLongOMP::PairBuckLongCoulLongOMP(LAMMPS *lmp) :
  PairBuckLongCoulLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 1;
  cut_respa = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairBuckLongCoulLongOMP::compute(int eflag, int vflag)
{
  ev_init(eflag,vflag);

  const int order1 = ewald_order & EWALD_COUL;
  const int order6 = ewald_order & EWALD_DISP;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (order6) {
      if (order1) {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,0,1,1>(ifrom, ito, thr);
                else eval<1,1,0,0,0,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,0,1,1>(ifrom, ito, thr);
                else eval<1,0,0,0,0,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,0,1,1>(ifrom, ito, thr);
              else eval<0,0,0,0,0,1,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,0,1,1>(ifrom, ito, thr);
                else eval<1,1,0,1,0,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,0,1,1>(ifrom, ito, thr);
                else eval<1,0,0,1,0,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,0,1,1>(ifrom, ito, thr);
              else eval<0,0,0,1,0,1,1>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,1,1,1>(ifrom, ito, thr);
                else eval<1,1,0,0,1,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,1,1,1>(ifrom, ito, thr);
                else eval<1,0,0,0,1,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,1,1,1>(ifrom, ito, thr);
              else eval<0,0,0,0,1,1,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,1,1,1>(ifrom, ito, thr);
                else eval<1,1,0,1,1,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,1,1,1>(ifrom, ito, thr);
                else eval<1,0,0,1,1,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,1,1,1>(ifrom, ito, thr);
              else eval<0,0,0,1,1,1,1>(ifrom, ito, thr);
            }
          }
        }
      } else {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,0,0,1>(ifrom, ito, thr);
                else eval<1,1,0,0,0,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,0,0,1>(ifrom, ito, thr);
                else eval<1,0,0,0,0,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,0,0,1>(ifrom, ito, thr);
              else eval<0,0,0,0,0,0,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,0,0,1>(ifrom, ito, thr);
                else eval<1,1,0,1,0,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,0,0,1>(ifrom, ito, thr);
                else eval<1,0,0,1,0,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,0,0,1>(ifrom, ito, thr);
              else eval<0,0,0,1,0,0,1>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,1,0,1>(ifrom, ito, thr);
                else eval<1,1,0,0,1,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,1,0,1>(ifrom, ito, thr);
                else eval<1,0,0,0,1,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,1,0,1>(ifrom, ito, thr);
              else eval<0,0,0,0,1,0,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,1,0,1>(ifrom, ito, thr);
                else eval<1,1,0,1,1,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,1,0,1>(ifrom, ito, thr);
                else eval<1,0,0,1,1,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,1,0,1>(ifrom, ito, thr);
              else eval<0,0,0,1,1,0,1>(ifrom, ito, thr);
            }
          }
        }
      }
    } else {
      if (order1) {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,0,1,0>(ifrom, ito, thr);
                else eval<1,1,0,0,0,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,0,1,0>(ifrom, ito, thr);
                else eval<1,0,0,0,0,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,0,1,0>(ifrom, ito, thr);
              else eval<0,0,0,0,0,1,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,0,1,0>(ifrom, ito, thr);
                else eval<1,1,0,1,0,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,0,1,0>(ifrom, ito, thr);
                else eval<1,0,0,1,0,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,0,1,0>(ifrom, ito, thr);
              else eval<0,0,0,1,0,1,0>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,1,1,0>(ifrom, ito, thr);
                else eval<1,1,0,0,1,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,1,1,0>(ifrom, ito, thr);
                else eval<1,0,0,0,1,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,1,1,0>(ifrom, ito, thr);
              else eval<0,0,0,0,1,1,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,1,1,0>(ifrom, ito, thr);
                else eval<1,1,0,1,1,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,1,1,0>(ifrom, ito, thr);
                else eval<1,0,0,1,1,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,1,1,0>(ifrom, ito, thr);
              else eval<0,0,0,1,1,1,0>(ifrom, ito, thr);
            }
          }
        }
      } else {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,0,0,0>(ifrom, ito, thr);
                else eval<1,1,0,0,0,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,0,0,0>(ifrom, ito, thr);
                else eval<1,0,0,0,0,0,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,0,0,0,0>(ifrom, ito, thr);
              else eval<0,0,0,0,0,0,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,0,0,0>(ifrom, ito, thr);
                else eval<1,1,0,1,0,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,0,0,0>(ifrom, ito, thr);
                else eval<1,0,0,1,0,0,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,0,0,0>(ifrom, ito, thr);
              else eval<0,0,0,1,0,0,0>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,0,1,0,0>(ifrom, ito, thr);
                else eval<1,1,0,0,1,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,0,1,0,0>(ifrom, ito, thr);
                else eval<1,0,0,0,1,0,0>(ifrom, ito, thr);
              }
            } else {
            if (force->newton_pair) eval<0,0,1,0,1,0,0>(ifrom, ito, thr);
              else eval<0,0,0,0,1,0,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval<1,1,1,1,1,0,0>(ifrom, ito, thr);
                else eval<1,1,0,1,1,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval<1,0,1,1,1,0,0>(ifrom, ito, thr);
                else eval<1,0,0,1,1,0,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval<0,0,1,1,1,0,0>(ifrom, ito, thr);
              else eval<0,0,0,1,1,0,0>(ifrom, ito, thr);
            }
          }
        }
      }
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */

void PairBuckLongCoulLongOMP::compute_inner()
{

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum_inner;
#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(0, 0, nall, nullptr, nullptr, nullptr, thr);
    eval_inner(ifrom, ito, thr);
    thr->timer(Timer::PAIR);

  }  // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

void PairBuckLongCoulLongOMP::compute_middle()
{

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum_middle;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(0, 0, nall, nullptr, nullptr, nullptr, thr);
    eval_middle(ifrom, ito, thr);
    thr->timer(Timer::PAIR);

  }  // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

void PairBuckLongCoulLongOMP::compute_outer(int eflag, int vflag)
{
  ev_init(eflag,vflag);
  const int order1 = ewald_order & EWALD_COUL;
  const int order6 = ewald_order & EWALD_DISP;

  const int nall = atom->nlocal + atom->nghost;
  const int nthreads = comm->nthreads;
  const int inum = list->inum;

#if defined(_OPENMP)
#pragma omp parallel LMP_DEFAULT_NONE LMP_SHARED(eflag,vflag)
#endif
  {
    int ifrom, ito, tid;

    loop_setup_thr(ifrom, ito, tid, inum, nthreads);
    ThrData *thr = fix->get_thr(tid);
    thr->timer(Timer::START);
    ev_setup_thr(eflag, vflag, nall, eatom, vatom, nullptr, thr);

    if (order6) {
      if (order1) {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,0,1,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,0,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,0,1,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,0,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,0,1,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,0,1,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,0,1,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,0,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,0,1,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,0,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,0,1,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,0,1,1>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,1,1,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,1,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,1,1,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,1,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,1,1,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,1,1,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,1,1,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,1,1,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,1,1,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,1,1,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,1,1,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,1,1,1>(ifrom, ito, thr);
            }
          }
        }
      } else {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,0,0,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,0,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,0,0,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,0,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,0,0,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,0,0,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,0,0,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,0,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,0,0,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,0,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,0,0,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,0,0,1>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,1,0,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,1,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,1,0,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,1,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,1,0,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,1,0,1>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,1,0,1>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,1,0,1>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,1,0,1>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,1,0,1>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,1,0,1>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,1,0,1>(ifrom, ito, thr);
            }
          }
        }
      }
    } else {
      if (order1) {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,0,1,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,0,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,0,1,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,0,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,0,1,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,0,1,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,0,1,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,0,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,0,1,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,0,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,0,1,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,0,1,0>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,1,1,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,1,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,1,1,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,1,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,1,1,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,1,1,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,1,1,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,1,1,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,1,1,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,1,1,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,1,1,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,1,1,0>(ifrom, ito, thr);
            }
          }
        }
      } else {
        if (!ndisptablebits) {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,0,0,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,0,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,0,0,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,0,0,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,0,0,0,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,0,0,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,0,0,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,0,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,0,0,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,0,0,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,0,0,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,0,0,0>(ifrom, ito, thr);
            }
          }
        } else {
          if (!ncoultablebits) {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,0,1,0,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,0,1,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,0,1,0,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,0,1,0,0>(ifrom, ito, thr);
              }
            } else {
            if (force->newton_pair) eval_outer<0,0,1,0,1,0,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,0,1,0,0>(ifrom, ito, thr);
            }
          } else {
            if (evflag) {
              if (eflag) {
                if (force->newton_pair) eval_outer<1,1,1,1,1,0,0>(ifrom, ito, thr);
                else eval_outer<1,1,0,1,1,0,0>(ifrom, ito, thr);
              } else {
                if (force->newton_pair) eval_outer<1,0,1,1,1,0,0>(ifrom, ito, thr);
                else eval_outer<1,0,0,1,1,0,0>(ifrom, ito, thr);
              }
            } else {
              if (force->newton_pair) eval_outer<0,0,1,1,1,0,0>(ifrom, ito, thr);
              else eval_outer<0,0,0,1,1,0,0>(ifrom, ito, thr);
            }
          }
        }
      }
    }

    thr->timer(Timer::PAIR);
    reduce_thr(this, eflag, vflag, thr);
  } // end of omp parallel region
}

/* ---------------------------------------------------------------------- */

template < const int EVFLAG, const int EFLAG,
           const int NEWTON_PAIR, const int CTABLE, const int DISPTABLE, const int ORDER1, const int ORDER6 >
void PairBuckLongCoulLongOMP::eval(int iifrom, int iito, ThrData * const thr)
{

  double evdwl,ecoul,fpair;
  evdwl = ecoul = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;

  int *ilist = list->ilist;

  // loop over neighbors of my atoms

  int i, ii, j, jj;
  int typei, typej, ni;
  double qi = 0.0, qri = 0.0, *cutsqi, *cut_bucksqi,
         *buck1i, *buck2i, *buckai, *buckci, *rhoinvi, *offseti;
  double r, rsq, r2inv, force_coul, force_buck;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;
  double xi[3], d[3];

  for (ii = iifrom; ii < iito; ++ii) {                        // loop over my atoms
    i = ilist[ii];
    if (ORDER1) qri = (qi = q[i])*qqrd2e;                // initialize constants
    typei = type[i];
    offseti = offset[typei];
    buck1i = buck1[typei]; buck2i = buck2[typei];
    buckai = buck_a[typei]; buckci = buck_c[typei]; rhoinvi = rhoinv[typei];
    cutsqi = cutsq[typei]; cut_bucksqi = cut_bucksq[typei];
    xi[0] = x[i][0]; xi[1] = x[i][1]; xi[2] = x[i][2];
    int *jlist = list->firstneigh[i];
    int jnum = list->numneigh[i];

    for (jj = 0; jj < jnum; jj++) {                           // loop over neighbors
      j = jlist[jj];
      ni = sbmask(j);
      j &= NEIGHMASK;

      d[0] = xi[0] - x[j][0];                                 // pair vector
      d[1] = xi[1] - x[j][1];
      d[2] = xi[2] - x[j][2];

      rsq = dot3(d, d);
      typej = type[j];
      if (rsq >= cutsqi[typej]) continue;
      r2inv = 1.0/rsq;
      r = sqrt(rsq);

      if (ORDER1 && (rsq < cut_coulsq)) {                // coulombic
        if (!CTABLE || rsq <= tabinnersq) {        // series real space
          double x = g_ewald*r;
          double s = qri*q[j], t = 1.0/(1.0+EWALD_P*x);
          if (ni == 0) {
            s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s;
            if (EFLAG) ecoul = t;
          } else {                                        // special case
            double f = s*(1.0-special_coul[ni])/r;
            s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-f;
            if (EFLAG) ecoul = t-f;
          }
        } else {                                               // table real space
          union_int_float_t t;
          t.f = rsq;
          const int k = (t.i & ncoulmask) >> ncoulshiftbits;
          double f = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
          if (ni == 0) {
            force_coul = qiqj*(ftable[k]+f*dftable[k]);
            if (EFLAG) ecoul = qiqj*(etable[k]+f*detable[k]);
          } else {                                        // special case
            t.f = (1.0-special_coul[ni])*(ctable[k]+f*dctable[k]);
            force_coul = qiqj*(ftable[k]+f*dftable[k]-(double)t.f);
            if (EFLAG) ecoul = qiqj*(etable[k]+f*detable[k]-(double)t.f);
          }
        }
      } else force_coul = ecoul = 0.0;

      if (rsq < cut_bucksqi[typej]) {                        // buckingham
        double rn = r2inv*r2inv*r2inv,
                        expr = exp(-r*rhoinvi[typej]);
        if (ORDER6) {                                        // long-range
          if (!DISPTABLE || rsq <= tabinnerdispsq) {
            double x2 = g2*rsq, a2 = 1.0/x2;
            x2 = a2*exp(-x2)*buckci[typej];
            if (ni == 0) {
              force_buck =
                r*expr*buck1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
              if (EFLAG) evdwl = expr*buckai[typej]-g6*((a2+1.0)*a2+0.5)*x2;
            } else {                                        // special case
              double f = special_lj[ni], t = rn*(1.0-f);
              force_buck = f*r*expr*buck1i[typej]-
                g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*buck2i[typej];
              if (EFLAG) evdwl = f*expr*buckai[typej] -
                           g6*((a2+1.0)*a2+0.5)*x2+t*buckci[typej];
            }
          } else {                                              //table real space
            union_int_float_t disp_t;
            disp_t.f = rsq;
            const int disp_k = (disp_t.i & ndispmask)>>ndispshiftbits;
            double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            if (ni == 0) {
              force_buck = r*expr*buck1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*buckci[typej];
              if (EFLAG) evdwl = expr*buckai[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*buckci[typej];
            } else {                                             //speial case
              double f = special_lj[ni], t = rn*(1.0-f);
              force_buck = f*r*expr*buck1i[typej] -(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*buckci[typej] +t*buck2i[typej];
              if (EFLAG) evdwl = f*expr*buckai[typej] -(edisptable[disp_k]+f_disp*dedisptable[disp_k])*buckci[typej]+t*buckci[typej];
            }
          }
        } else {                                                // cut
          if (ni == 0) {
            force_buck = r*expr*buck1i[typej]-rn*buck2i[typej];
            if (EFLAG) evdwl = expr*buckai[typej] -
                         rn*buckci[typej]-offseti[typej];
          } else {                                        // special case
            double f = special_lj[ni];
            force_buck = f*(r*expr*buck1i[typej]-rn*buck2i[typej]);
            if (EFLAG)
              evdwl = f*(expr*buckai[typej]-rn*buckci[typej]-offseti[typej]);
          }
        }
      } else force_buck = evdwl = 0.0;

      fpair = (force_coul+force_buck)*r2inv;

      double fdx = d[0]*fpair, fdy = d[1]*fpair, fdz = d[2]*fpair;
      f[i][0] += fdx;
      f[i][1] += fdy;
      f[i][2] += fdz;
      if (NEWTON_PAIR || j < nlocal) {
        f[j][0] -= fdx;
        f[j][1] -= fdy;
        f[j][2] -= fdz;
      }

      if (EVFLAG) ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                               evdwl,ecoul,fpair,d[0],d[1],d[2],thr);
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairBuckLongCoulLongOMP::eval_inner(int iifrom, int iito, ThrData * const thr)
{
  double r, rsq, r2inv, force_coul = 0.0, force_buck, fpair;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;

  int *ilist = list->ilist_inner;

  const int newton_pair = force->newton_pair;

  const double cut_out_on = cut_respa[0];
  const double cut_out_off = cut_respa[1];

  const double cut_out_diff = cut_out_off - cut_out_on;
  const double cut_out_on_sq = cut_out_on*cut_out_on;
  const double cut_out_off_sq = cut_out_off*cut_out_off;

  int typei, typej, ni;
  const int order1 = (ewald_order | ~ewald_off) & EWALD_COUL;
  int i, j, ii, jj;
  double qri, *cut_bucksqi, *buck1i, *buck2i, *rhoinvi;
  double xi[3], d[3];

  for (ii = iifrom; ii < iito; ++ii) {                        // loop over my atoms
    i = ilist[ii];
    if (order1) qri = qqrd2e*q[i];
    xi[0] = x[i][0]; xi[1] = x[i][1]; xi[2] = x[i][2];
    typei = type[i];
    cut_bucksqi = cut_bucksq[typei];
    buck1i = buck1[typei]; buck2i = buck2[typei]; rhoinvi = rhoinv[typei];
    int *jlist = list->firstneigh_inner[i];
    int jnum = list->numneigh_inner[i];

    for (jj = 0; jj < jnum; jj++) {                           // loop over neighbors
      j = jlist[jj];
      ni = sbmask(j);
      j &= NEIGHMASK;

      d[0] = xi[0] - x[j][0];                                 // pair vector
      d[1] = xi[1] - x[j][1];
      d[2] = xi[2] - x[j][2];

      rsq = dot3(d, d);
      if (rsq >= cut_out_off_sq) continue;
      r2inv = 1.0/rsq;
      r = sqrt(rsq);

      if (order1 && (rsq < cut_coulsq))                        // coulombic
        force_coul = ni == 0 ?
          qri*q[j]/r : qri*q[j]/r*special_coul[ni];

      typej = type[j];
      if (rsq < cut_bucksqi[typej]) {                          // buckingham
        double rn = r2inv*r2inv*r2inv,
                        expr = exp(-r*rhoinvi[typej]);
        force_buck = ni == 0 ?
          (r*expr*buck1i[typej]-rn*buck2i[typej]) :
          (r*expr*buck1i[typej]-rn*buck2i[typej])*special_lj[ni];
      } else force_buck = 0.0;

      fpair = (force_coul + force_buck) * r2inv;

      if (rsq > cut_out_on_sq) {                        // switching
        double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
        fpair  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
      }

      double fdx = d[0]*fpair, fdy = d[1]*fpair, fdz = d[2]*fpair;
      f[i][0] += fdx;                                          // force update
      f[i][1] += fdy;
      f[i][2] += fdz;
      if (newton_pair || j < nlocal) {
        f[j][0] -= fdx;
        f[j][1] -= fdy;
        f[j][2] -= fdz;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void PairBuckLongCoulLongOMP::eval_middle(int iifrom, int iito, ThrData * const thr)
{
  double r, rsq, r2inv, force_coul = 0.0, force_buck, fpair;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;

  int *ilist = list->ilist_middle;

  const int newton_pair = force->newton_pair;

  const double cut_in_off = cut_respa[0];
  const double cut_in_on = cut_respa[1];
  const double cut_out_on = cut_respa[2];
  const double cut_out_off = cut_respa[3];

  const double cut_in_diff = cut_in_on - cut_in_off;
  const double cut_out_diff = cut_out_off - cut_out_on;
  const double cut_in_off_sq = cut_in_off*cut_in_off;
  const double cut_in_on_sq = cut_in_on*cut_in_on;
  const double cut_out_on_sq = cut_out_on*cut_out_on;
  const double cut_out_off_sq = cut_out_off*cut_out_off;

  int typei, typej, ni;
  const int order1 = (ewald_order | ~ewald_off) & EWALD_COUL;
  int i, j, ii, jj;
  double qri, *cut_bucksqi, *buck1i, *buck2i, *rhoinvi;
  double xi[3], d[3];

  for (ii = iifrom; ii < iito; ++ii) {                        // loop over my atoms
    i = ilist[ii];
    if (order1) qri = qqrd2e*q[i];
    xi[0] = x[i][0]; xi[1] = x[i][1]; xi[2] = x[i][2];
    typei = type[i];
    cut_bucksqi = cut_bucksq[typei];
    buck1i = buck1[typei]; buck2i = buck2[typei]; rhoinvi = rhoinv[typei];
    int *jlist = list->firstneigh_middle[i];
    int jnum = list->numneigh_middle[i];

    for (jj = 0; jj < jnum; jj++) {                           // loop over neighbors
      j = jlist[jj];
      ni = sbmask(j);
      j &= NEIGHMASK;

      d[0] = xi[0] - x[j][0];                                 // pair vector
      d[1] = xi[1] - x[j][1];
      d[2] = xi[2] - x[j][2];

      rsq = dot3(d, d);
      if (rsq >= cut_out_off_sq) continue;
      if (rsq <= cut_in_off_sq) continue;
      r2inv = 1.0/rsq;
      r = sqrt(rsq);

      if (order1 && (rsq < cut_coulsq))                        // coulombic
        force_coul = ni == 0 ?
          qri*q[j]/r : qri*q[j]/r*special_coul[ni];

      typej = type[j];
      if (rsq < cut_bucksqi[typej]) {                          // buckingham
        double rn = r2inv*r2inv*r2inv,
                        expr = exp(-r*rhoinvi[typej]);
        force_buck = ni == 0 ?
          (r*expr*buck1i[typej]-rn*buck2i[typej]) :
          (r*expr*buck1i[typej]-rn*buck2i[typej])*special_lj[ni];
      } else force_buck = 0.0;

      fpair = (force_coul + force_buck) * r2inv;

      if (rsq < cut_in_on_sq) {                                // switching
        double rsw = (sqrt(rsq) - cut_in_off)/cut_in_diff;
        fpair  *= rsw*rsw*(3.0 - 2.0*rsw);
      }
      if (rsq > cut_out_on_sq) {
        double rsw = (sqrt(rsq) - cut_out_on)/cut_out_diff;
        fpair  *= 1.0 + rsw*rsw*(2.0*rsw-3.0);
      }

      double fdx = d[0]*fpair, fdy = d[1]*fpair, fdz = d[2]*fpair;
      f[i][0] += fdx;                                          // force update
      f[i][1] += fdy;
      f[i][2] += fdz;
      if (newton_pair || j < nlocal) {
        f[j][0] -= fdx;
        f[j][1] -= fdy;
        f[j][2] -= fdz;
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

template < const int EVFLAG, const int EFLAG,
           const int NEWTON_PAIR, const int CTABLE, const int DISPTABLE, const int ORDER1, const int ORDER6 >
void PairBuckLongCoulLongOMP::eval_outer(int iiform, int iito, ThrData * const thr)
{
  double evdwl,ecoul,fpair,fvirial;
  evdwl = ecoul = 0.0;

  const double * const * const x = atom->x;
  double * const * const f = thr->get_f();
  const double * const q = atom->q;
  const int * const type = atom->type;
  const int nlocal = atom->nlocal;
  const double * const special_coul = force->special_coul;
  const double * const special_lj = force->special_lj;
  const double qqrd2e = force->qqrd2e;

  int *ilist = list->ilist;

  int i, j, ii, jj;
  int typei, typej, ni, respa_flag;
  double qi = 0.0, qri = 0.0;
  double *cutsqi, *cut_bucksqi, *buck1i, *buck2i, *buckai, *buckci, *rhoinvi, *offseti;
  double r, rsq, r2inv, force_coul, force_buck;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;
  double respa_buck = 0.0, respa_coul = 0.0, frespa = 0.0;
  double xi[3], d[3];

  const double cut_in_off = cut_respa[2];
  const double cut_in_on = cut_respa[3];

  const double cut_in_diff = cut_in_on - cut_in_off;
  const double cut_in_off_sq = cut_in_off*cut_in_off;
  const double cut_in_on_sq = cut_in_on*cut_in_on;

  for (ii = iiform; ii < iito; ++ii) {                        // loop over my atoms
    i = ilist[ii];
    if (ORDER1) qri = (qi = q[i])*qqrd2e;                // initialize constants
    typei = type[i];
    offseti = offset[typei];
    buck1i = buck1[typei]; buck2i = buck2[typei];
    buckai = buck_a[typei]; buckci = buck_c[typei]; rhoinvi = rhoinv[typei];
    cutsqi = cutsq[typei]; cut_bucksqi = cut_bucksq[typei];
    xi[0] = x[i][0]; xi[1] = x[i][1]; xi[2] = x[i][2];
    int *jlist = list->firstneigh[i];
    int jnum = list->numneigh[i];

    for (jj = 0; jj < jnum; jj++) {                           // loop over neighbors
      j = jlist[jj];
      ni = sbmask(j);
      j &= NEIGHMASK;

      d[0] = xi[0] - x[j][0];                                 // pair vector
      d[1] = xi[1] - x[j][1];
      d[2] = xi[2] - x[j][2];

      rsq = dot3(d, d);
      typej = type[j];
      if (rsq >= cutsqi[typej]) continue;
      r2inv = 1.0/rsq;
      r = sqrt(rsq);

      frespa = 1.0;      //check whether and how to compute respa corrections
      respa_coul = 0.0;
      respa_buck = 0.0;
      respa_flag = rsq < cut_in_on_sq ? 1 : 0;
      if (respa_flag && (rsq > cut_in_off_sq)) {
        double rsw = (r-cut_in_off)/cut_in_diff;
        frespa = 1-rsw*rsw*(3.0-2.0*rsw);
      }

      if (ORDER1 && (rsq < cut_coulsq)) {                // coulombic
        if (!CTABLE || rsq <= tabinnersq) {        // series real space
          double s = qri*q[j];
          if (respa_flag)                                // correct for respa
            respa_coul = ni == 0 ? frespa*s/r : frespa*s/r*special_coul[ni];
          double x = g_ewald*r, t = 1.0/(1.0+EWALD_P*x);
          if (ni == 0) {
            s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-respa_coul;
            if (EFLAG) ecoul = t;
          } else {                                        // correct for special
            double ri = s*(1.0-special_coul[ni])/r; s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-ri-respa_coul;
            if (EFLAG) ecoul = t-ri;
          }
        } else {                                             // table real space
          if (respa_flag) {
            double s = qri*q[j];
            respa_coul = ni == 0 ? frespa*s/r : frespa*s/r*special_coul[ni];
          }
          union_int_float_t t;
          t.f = rsq;
          const int k = (t.i & ncoulmask) >> ncoulshiftbits;
          double f = (rsq-rtable[k])*drtable[k], qiqj = qi*q[j];
          if (ni == 0) {
            force_coul = qiqj*(ftable[k]+f*dftable[k]);
            if (EFLAG) ecoul = qiqj*(etable[k]+f*detable[k]);
          } else {                                        // correct for special
            t.f = (1.0-special_coul[ni])*(ctable[k]+f*dctable[k]);
            force_coul = qiqj*(ftable[k]+f*dftable[k]-(double)t.f);
            if (EFLAG) {
              t.f = (1.0-special_coul[ni])*(ptable[k]+f*dptable[k]);
              ecoul = qiqj*(etable[k]+f*detable[k]-(double)t.f);
            }
          }
        }
      } else force_coul = respa_coul = ecoul = 0.0;

      if (rsq < cut_bucksqi[typej]) {                        // buckingham
        double rn = r2inv*r2inv*r2inv,
                        expr = exp(-r*rhoinvi[typej]);
        if (respa_flag) respa_buck = ni == 0 ?                 // correct for respa
            frespa*(r*expr*buck1i[typej]-rn*buck2i[typej]) :
            frespa*(r*expr*buck1i[typej]-rn*buck2i[typej])*special_lj[ni];
        if (ORDER6) {                                        // long-range form
          if (!DISPTABLE || rsq <= tabinnerdispsq) {
            double x2 = g2*rsq, a2 = 1.0/x2;
            x2 = a2*exp(-x2)*buckci[typej];
            if (ni == 0) {
              force_buck =
                r*expr*buck1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq-respa_buck;
              if (EFLAG) evdwl = expr*buckai[typej]-g6*((a2+1.0)*a2+0.5)*x2;
            } else {                                        // correct for special
              double f = special_lj[ni], t = rn*(1.0-f);
              force_buck = f*r*expr*buck1i[typej]-
                g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*buck2i[typej]-respa_buck;
              if (EFLAG) evdwl = f*expr*buckai[typej] -
                           g6*((a2+1.0)*a2+0.5)*x2+t*buckci[typej];
            }
          } else {          // table real space
            union_int_float_t disp_t;
            disp_t.f = rsq;
            const int disp_k = (disp_t.i & ndispmask)>>ndispshiftbits;
            double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            double rn = r2inv*r2inv*r2inv;
            if (ni == 0) {
              force_buck = r*expr*buck1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*buckci[typej]-respa_buck;
              if (EFLAG) evdwl =  expr*buckai[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*buckci[typej];
            } else {                             //special case
              double f = special_lj[ni], t = rn*(1.0-f);
              force_buck = f*r*expr*buck1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*buckci[typej]+t*buck2i[typej]-respa_buck;
              if (EFLAG) evdwl = f*expr*buckai[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*buckci[typej]+t*buckci[typej];
            }
          }
        } else {                                                // cut form
          if (ni == 0) {
            force_buck = r*expr*buck1i[typej]-rn*buck2i[typej]-respa_buck;
            if (EFLAG)
              evdwl = expr*buckai[typej]-rn*buckci[typej]-offseti[typej];
          } else {                                        // correct for special
            double f = special_lj[ni];
            force_buck = f*(r*expr*buck1i[typej]-rn*buck2i[typej])-respa_buck;
            if (EFLAG)
              evdwl = f*(expr*buckai[typej]-rn*buckci[typej]-offseti[typej]);
          }
        }
      } else force_buck = respa_buck = evdwl = 0.0;

      fpair = (force_coul+force_buck)*r2inv;

      double fdx = d[0]*fpair, fdy = d[1]*fpair, fdz = d[2]*fpair;
      f[i][0] += fdx;
      f[i][1] += fdy;
      f[i][2] += fdz;
      if (NEWTON_PAIR || j < nlocal) {
        f[j][0] -= fdx;
        f[j][1] -= fdy;
        f[j][2] -= fdz;
      }

      if (EVFLAG) {
        fvirial = (force_coul + force_buck + respa_coul + respa_buck)*r2inv;
        ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                     evdwl,ecoul,fvirial,d[0],d[1],d[2],thr);
      }
    }
  }
}

