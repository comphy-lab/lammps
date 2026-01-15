/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef ATOM_CLASS
// clang-format off
AtomStyle(hybrid/kk,AtomVecHybridKokkos);
AtomStyle(hybrid/kk/device,AtomVecHybridKokkos);
AtomStyle(hybrid/kk/host,AtomVecHybridKokkos);
// clang-format on
#else

// clang-format off
#ifndef LMP_ATOM_VEC_HYBRID_KOKKOS_H
#define LMP_ATOM_VEC_HYBRID_KOKKOS_H

#include "atom_vec_kokkos.h"
#include "atom_vec_hybrid.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

class AtomVecHybridKokkos : public AtomVecKokkos, public AtomVecHybrid {
 public:
  AtomVecHybridKokkos(class LAMMPS *);
  void init() override;

  void grow(int) override;
  void sort_kokkos(Kokkos::BinSort<KeyViewType, BinOp> &Sorter) override;

  void pack_comm_bonus_kokkos(const int &n, const DAT::tdual_int_1d &list,
                              const DAT::tdual_double_2d_lr &buf) override;

  void unpack_comm_bonus_kokkos(const int &n, const int &nfirst,
                                const DAT::tdual_double_2d_lr &buf) override;

  void pack_border_bonus_kokkos(int n, DAT::tdual_int_1d k_sendlist,
                                DAT::tdual_double_2d_lr &buf,
                                ExecutionSpace space) override;

  void unpack_border_bonus_kokkos(const int &n, const int &nfirst,
                                  const DAT::tdual_double_2d_lr &buf,
                                  ExecutionSpace space) override;

  void pack_exchange_bonus_kokkos(const int &nsend, DAT::tdual_double_2d_lr &buf,
                                  DAT::tdual_int_1d k_sendlist,
                                  DAT::tdual_int_1d k_copylist,
                                  DAT::tdual_int_1d k_copylist_bonus,
                                  ExecutionSpace space) override;

  void unpack_exchange_bonus_kokkos(DAT::tdual_double_2d_lr &k_buf,
                                    int nrecv,
                                    ExecutionSpace space,
                                    DAT::tdual_int_1d &k_indices) override;

  void set_size_exchange() override;
  void sync(ExecutionSpace space, uint64_t mask) override;
  void modified(ExecutionSpace space, uint64_t mask) override;
  void sync_pinned(ExecutionSpace space, uint64_t mask, int async_flag = 0) override;
};

} // namespace LAMMPS_NS

#endif
#endif
