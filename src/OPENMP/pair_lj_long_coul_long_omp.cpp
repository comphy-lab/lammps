// clang-format off
/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   This software is distributed under the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include "pair_lj_long_coul_long_omp.h"

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
// order 6 = LJ dispersion (1/r^6).

enum { EWALD_COUL = 1 << 1, EWALD_DISP = 1 << 6 };

/* ---------------------------------------------------------------------- */

PairLJLongCoulLongOMP::PairLJLongCoulLongOMP(LAMMPS *lmp) :
  PairLJLongCoulLong(lmp), ThrOMP(lmp, THR_PAIR)
{
  suffix_flag |= Suffix::OMP;
  respa_enable = 1;
  cut_respa = nullptr;
}

/* ---------------------------------------------------------------------- */

void PairLJLongCoulLongOMP::compute(int eflag, int vflag)
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

void PairLJLongCoulLongOMP::compute_inner()
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

void PairLJLongCoulLongOMP::compute_middle()
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

void PairLJLongCoulLongOMP::compute_outer(int eflag, int vflag)
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
           const int NEWTON_PAIR, const int CTABLE, const int LJTABLE, const int ORDER1, const int ORDER6 >
void PairLJLongCoulLongOMP::eval(int iifrom, int iito, ThrData * const thr)
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
  double qi, qri, *cutsqi, *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double rsq, r2inv, force_coul, force_lj;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;
  double xi[3], d[3];

  for (ii = iifrom; ii < iito; ++ii) {                        // loop over my atoms
    i = ilist[ii];
    if (ORDER1) qri = (qi = q[i])*qqrd2e;                // initialize constants
    typei = type[i];
    offseti = offset[typei];
    lj1i = lj1[typei]; lj2i = lj2[typei]; lj3i = lj3[typei]; lj4i = lj4[typei];
    cutsqi = cutsq[typei]; cut_ljsqi = cut_ljsq[typei];
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

      if (ORDER1 && (rsq < cut_coulsq)) {                // coulombic
        if (!CTABLE || rsq <= tabinnersq) {        // series real space
          double r = sqrt(rsq), x = g_ewald*r;
          double s = qri*q[j], t = 1.0/(1.0+EWALD_P*x);
          if (ni == 0) {
            s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s;
            if (EFLAG) ecoul = t;
          } else {                                        // special case
            r = s*(1.0-special_coul[ni])/r; s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-r;
            if (EFLAG) ecoul = t-r;
          }
        } else {                                             // table real space
          union_int_float_t t;
          t.f = rsq;
          const int k = (t.i & ncoulmask)>>ncoulshiftbits;
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

      if (rsq < cut_ljsqi[typej]) {                        // lj
        if (ORDER6) {                                        // long-range lj
          if (!LJTABLE || rsq <= tabinnerdispsq) {            //series real space
            double rn = r2inv*r2inv*r2inv;
            double x2 = g2*rsq, a2 = 1.0/x2;
            x2 = a2*exp(-x2)*lj4i[typej];
            if (ni == 0) {
              force_lj =
                (rn*=rn)*lj1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq;
              if (EFLAG)
                evdwl = rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2;
            } else {                                        // special case
              double f = special_lj[ni], t = rn*(1.0-f);
              force_lj = f*(rn *= rn)*lj1i[typej]-
                g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2i[typej];
              if (EFLAG)
                evdwl = f*rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2+t*lj4i[typej];
            }
          } else {                                        // table real space
            union_int_float_t disp_t;
            disp_t.f = rsq;
            const int disp_k = (disp_t.i & ndispmask)>>ndispshiftbits;
            double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            double rn = r2inv*r2inv*r2inv;
            if (ni == 0) {
              force_lj = (rn*=rn)*lj1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[typej];
              if (EFLAG) evdwl = rn*lj3i[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[typej];
            } else {                                      // special case
              double f = special_lj[ni], t = rn*(1.0-f);
              force_lj = f*(rn *= rn)*lj1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[typej]+t*lj2i[typej];
              if (EFLAG) evdwl = f*rn*lj3i[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[typej]+t*lj4i[typej];
            }
          }
        } else {                                                // cut lj
          double rn = r2inv*r2inv*r2inv;
          if (ni == 0) {
            force_lj = rn*(rn*lj1i[typej]-lj2i[typej]);
            if (EFLAG) evdwl = rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej];
          } else {                                        // special case
            double f = special_lj[ni];
            force_lj = f*rn*(rn*lj1i[typej]-lj2i[typej]);
            if (EFLAG)
              evdwl = f * (rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej]);
          }
        }
      } else force_lj = evdwl = 0.0;

      fpair = (force_coul+force_lj)*r2inv;

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

void PairLJLongCoulLongOMP::eval_inner(int iifrom, int iito, ThrData * const thr)
{
  double rsq, r2inv, force_coul = 0.0, force_lj, fpair;

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
  double qri, *cut_ljsqi, *lj1i, *lj2i;
  double xi[3], d[3];

  for (ii = iifrom; ii < iito; ++ii) {                        // loop over my atoms
    i = ilist[ii];
    xi[0] = x[i][0]; xi[1] = x[i][1]; xi[2] = x[i][2];
    typei = type[i];
    cut_ljsqi = cut_ljsq[typei];
    lj1i = lj1[typei]; lj2i = lj2[typei];
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

      if (order1 && (rsq < cut_coulsq)) {                       // coulombic
        qri = qqrd2e*q[i];
        force_coul = ni == 0 ?
          qri*q[j]*sqrt(r2inv) : qri*q[j]*sqrt(r2inv)*special_coul[ni];
      }

      typej = type[j];
      if (rsq < cut_ljsqi[typej]) {                          // lennard-jones
        double rn = r2inv*r2inv*r2inv;
        force_lj = ni == 0 ?
          rn*(rn*lj1i[typej]-lj2i[typej]) :
          rn*(rn*lj1i[typej]-lj2i[typej])*special_lj[ni];
      } else force_lj = 0.0;

      fpair = (force_coul + force_lj) * r2inv;

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

void PairLJLongCoulLongOMP::eval_middle(int iifrom, int iito, ThrData * const thr)
{
  double rsq, r2inv, force_coul = 0.0, force_lj, fpair;

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
  double qri, *cut_ljsqi, *lj1i, *lj2i;
  double xi[3], d[3];

  for (ii = iifrom; ii < iito; ++ii) {                        // loop over my atoms
    i = ilist[ii];
    if (order1) qri = qqrd2e*q[i];
    xi[0] = x[i][0]; xi[1] = x[i][1]; xi[2] = x[i][2];
    typei = type[i];
    cut_ljsqi = cut_ljsq[typei];
    lj1i = lj1[typei]; lj2i = lj2[typei];
    int *jlist = list->firstneigh_middle[i];
    int jnum = list->numneigh_middle[i];

    for (jj = 0; jj < jnum; jj++) {
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

      if (order1 && (rsq < cut_coulsq))                        // coulombic
        force_coul = ni == 0 ?
          qri*q[j]*sqrt(r2inv) : qri*q[j]*sqrt(r2inv)*special_coul[ni];

      typej = type[j];
      if (rsq < cut_ljsqi[typej]) {                          // lennard-jones
        double rn = r2inv*r2inv*r2inv;
        force_lj = ni == 0 ?
          rn*(rn*lj1i[typej]-lj2i[typej]) :
          rn*(rn*lj1i[typej]-lj2i[typej])*special_lj[ni];
      } else force_lj = 0.0;

      fpair = (force_coul + force_lj) * r2inv;

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
           const int NEWTON_PAIR, const int CTABLE, const int LJTABLE, const int ORDER1, const int ORDER6 >
void PairLJLongCoulLongOMP::eval_outer(int iiform, int iito, ThrData * const thr)
{
  double evdwl,ecoul,fvirial,fpair;
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
  double *cutsqi, *cut_ljsqi, *lj1i, *lj2i, *lj3i, *lj4i, *offseti;
  double rsq, r2inv, force_coul, force_lj;
  double g2 = g_ewald_6*g_ewald_6, g6 = g2*g2*g2, g8 = g6*g2;
  double respa_lj = 0.0, respa_coul = 0.0, frespa = 0.0;
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
    lj1i = lj1[typei]; lj2i = lj2[typei]; lj3i = lj3[typei]; lj4i = lj4[typei];
    cutsqi = cutsq[typei]; cut_ljsqi = cut_ljsq[typei];
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

      frespa = 1.0;                                       // check whether and how to compute respa corrections
      respa_coul = 0;
      respa_lj = 0;
      respa_flag = rsq < cut_in_on_sq ? 1 : 0;
      if (respa_flag && (rsq > cut_in_off_sq)) {
        double rsw = (sqrt(rsq)-cut_in_off)/cut_in_diff;
        frespa = 1-rsw*rsw*(3.0-2.0*rsw);
      }

      if (ORDER1 && (rsq < cut_coulsq)) {                // coulombic
        if (!CTABLE || rsq <= tabinnersq) {        // series real space
          double r = sqrt(rsq), s = qri*q[j];
          if (respa_flag)                                // correct for respa
            respa_coul = ni == 0 ? frespa*s/r : frespa*s/r*special_coul[ni];
          double x = g_ewald*r, t = 1.0/(1.0+EWALD_P*x);
          if (ni == 0) {
            s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-respa_coul;
            if (EFLAG) ecoul = t;
          } else {                                        // correct for special
            r = s*(1.0-special_coul[ni])/r; s *= g_ewald*exp(-x*x);
            force_coul = (t *= ((((t*A5+A4)*t+A3)*t+A2)*t+A1)*s/x)+EWALD_F*s-r-respa_coul;
            if (EFLAG) ecoul = t-r;
          }
        } else {                                              // table real space
          if (respa_flag) {
            double r = sqrt(rsq), s = qri*q[j];
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

      if (rsq < cut_ljsqi[typej]) {                        // lennard-jones
        double rn = r2inv*r2inv*r2inv;
        if (respa_flag) respa_lj = ni == 0 ?                 // correct for respa
            frespa*rn*(rn*lj1i[typej]-lj2i[typej]) :
            frespa*rn*(rn*lj1i[typej]-lj2i[typej])*special_lj[ni];
        if (ORDER6) {                                        // long-range form
          if (!LJTABLE || rsq <= tabinnerdispsq) {
            double x2 = g2*rsq, a2 = 1.0/x2;
            x2 = a2*exp(-x2)*lj4i[typej];
            if (ni == 0) {
              force_lj =
                (rn*=rn)*lj1i[typej]-g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq-respa_lj;
              if (EFLAG) evdwl = rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2;
            } else {                                        // correct for special
              double f = special_lj[ni], t = rn*(1.0-f);
              force_lj = f*(rn *= rn)*lj1i[typej]-
                g8*(((6.0*a2+6.0)*a2+3.0)*a2+1.0)*x2*rsq+t*lj2i[typej]-respa_lj;
              if (EFLAG)
                evdwl = f*rn*lj3i[typej]-g6*((a2+1.0)*a2+0.5)*x2+t*lj4i[typej];
            }
          } else {                                                // table real space
            union_int_float_t disp_t;
            disp_t.f = rsq;
            const int disp_k = (disp_t.i & ndispmask)>>ndispshiftbits;
            double f_disp = (rsq-rdisptable[disp_k])*drdisptable[disp_k];
            double rn = r2inv*r2inv*r2inv;
            if (ni == 0) {
              force_lj = (rn*=rn)*lj1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[typej]-respa_lj;
              if (EFLAG) evdwl = rn*lj3i[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[typej];
            } else {                                      // special case
              double f = special_lj[ni], t = rn*(1.0-f);
              force_lj = f*(rn *= rn)*lj1i[typej]-(fdisptable[disp_k]+f_disp*dfdisptable[disp_k])*lj4i[typej]+t*lj2i[typej]-respa_lj;
              if (EFLAG) evdwl = f*rn*lj3i[typej]-(edisptable[disp_k]+f_disp*dedisptable[disp_k])*lj4i[typej]+t*lj4i[typej];
            }
          }
        } else {                                                // cut form
          if (ni == 0) {
            force_lj = rn*(rn*lj1i[typej]-lj2i[typej])-respa_lj;
            if (EFLAG) evdwl = rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej];
          } else {                                        // correct for special
            double f = special_lj[ni];
            force_lj = f*rn*(rn*lj1i[typej]-lj2i[typej])-respa_lj;
            if (EFLAG)
              evdwl = f*(rn*(rn*lj3i[typej]-lj4i[typej])-offseti[typej]);
          }
        }
      } else force_lj = respa_lj = evdwl = 0.0;

      fpair = (force_coul+force_lj)*r2inv;

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
        fvirial = (force_coul + force_lj + respa_coul + respa_lj)*r2inv;
        ev_tally_thr(this,i,j,nlocal,NEWTON_PAIR,
                     evdwl,ecoul,fvirial,d[0],d[1],d[2],thr);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

double PairLJLongCoulLongOMP::memory_usage()
{
  double bytes = memory_usage_thr();
  bytes += PairLJLongCoulLong::memory_usage();

  return bytes;
}
