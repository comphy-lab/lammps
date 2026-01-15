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

#include "atom_vec_hybrid_kokkos.h"

#include "atom_kokkos.h"
#include "atom_masks.h"
#include "domain.h"
#include "error.h"
#include "kokkos.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

AtomVecHybridKokkos::AtomVecHybridKokkos(LAMMPS *lmp) : AtomVec(lmp),
AtomVecKokkos(lmp), AtomVecHybrid(lmp)
{
}

/* ----------------------------------------------------------------------
   process field strings to initialize data structs for all other methods
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::init()
{
  AtomVecHybrid::init();

  set_atom_masks();
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::grow(int n)
{
  for (int k = 0; k < nstyles; k++) styles[k]->grow(n);
  nmax = atomKK->k_x.view_host().extent(0);

  tag = atom->tag;
  type = atom->type;
  mask = atom->mask;
  image = atom->image;
  x = atom->x;
  v = atom->v;
  f = atom->f;
}

// TODO: move dynamic_cast into init

/* ----------------------------------------------------------------------
   sort atom arrays on device
------------------------------------------------------------------------- */

void AtomVecHybridKokkos::sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter)
{
  for (int k = 0; k < nstyles; k++)
    (dynamic_cast<AtomVecKokkos*>(styles[k]))->sort_kokkos(Sorter);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::pack_comm_bonus_kokkos(const int &n, const DAT::tdual_int_1d &list,
                                                 const DAT::tdual_double_2d_lr &buf)
{
  for (int k = 0; k < nstyles; k++)
    (dynamic_cast<AtomVecKokkos*>(styles[k]))->
      pack_comm_bonus_kokkos(n,list,buf);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::unpack_comm_bonus_kokkos(const int &n, const int &nfirst,
                                                   const DAT::tdual_double_2d_lr &buf)
{
  for (int k = 0; k < nstyles; k++)
    (dynamic_cast<AtomVecKokkos*>(styles[k]))->
      unpack_comm_bonus_kokkos(n,nfirst,buf);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::pack_border_bonus_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                                   DAT::tdual_double_2d_lr &buf,
                                                   ExecutionSpace space)
{
  for (int k = 0; k < nstyles; k++)
    (dynamic_cast<AtomVecKokkos*>(styles[k]))->
      pack_border_bonus_kokkos(n,k_sendlist,buf,space);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::unpack_border_bonus_kokkos(const int &n, const int &nfirst,
                                                     const DAT::tdual_double_2d_lr &buf,
                                                     ExecutionSpace space)
{
  for (int k = 0; k < nstyles; k++)
    (dynamic_cast<AtomVecKokkos*>(styles[k]))->
      unpack_border_bonus_kokkos(n,nfirst,buf,space);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::pack_exchange_bonus_kokkos(const int &nsend, DAT::tdual_double_2d_lr &buf,
                                                     DAT::tdual_int_1d k_sendlist,
                                                     DAT::tdual_int_1d k_copylist,
                                                     DAT::tdual_int_1d k_copylist_bonus,
                                                     ExecutionSpace space)
{
  for (int k = 0; k < nstyles; k++)
    (dynamic_cast<AtomVecKokkos*>(styles[k]))->
      pack_exchange_bonus_kokkos(nsend,buf,k_sendlist,k_copylist,
                                 k_copylist_bonus,space);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::unpack_exchange_bonus_kokkos(DAT::tdual_double_2d_lr &k_buf,
                                                       int nrecv,
                                                       ExecutionSpace space,
                                                       DAT::tdual_int_1d &k_indices)
{
  for (int k = 0; k < nstyles; k++)
    (dynamic_cast<AtomVecKokkos*>(styles[k]))->
      unpack_exchange_bonus_kokkos(k_buf,nrecv,space,k_indices);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::set_size_exchange()
{
  AtomVecKokkos::set_size_exchange();

 size_exchange_bonus = 0;
 for (int k = 0; k < nstyles; k++)
    size_exchange_bonus += (dynamic_cast<AtomVecKokkos*>(styles[k]))->
      size_exchange_bonus;

  size_exchange += size_exchange_bonus;

 for (int k = 0; k < nstyles; k++)
    (dynamic_cast<AtomVecKokkos*>(styles[k]))->size_exchange = size_exchange;
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::sync(ExecutionSpace space, uint64_t h_mask)
{
  for (int k = 0; k < nstyles; k++) (dynamic_cast<AtomVecKokkos*>(styles[k]))->sync(space,h_mask);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::sync_pinned(ExecutionSpace space, uint64_t h_mask, int async_flag)
{
  for (int k = 0; k < nstyles; k++) (dynamic_cast<AtomVecKokkos*>(styles[k]))->sync_pinned(space,h_mask,async_flag);
}

/* ---------------------------------------------------------------------- */

void AtomVecHybridKokkos::modified(ExecutionSpace space, uint64_t h_mask)
{
  for (int k = 0; k < nstyles; k++) (dynamic_cast<AtomVecKokkos*>(styles[k]))->modified(space,h_mask);
}
