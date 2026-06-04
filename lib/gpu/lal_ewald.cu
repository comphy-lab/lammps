// **************************************************************************
//                                 ewald.cu
//                             -------------------
//                            W. Michael Brown (ORNL)
//                            Axel Kohlmeyer (Temple)
//
//  Device code for Ewald (reciprocal-space) acceleration
//
// __________________________________________________________________________
//    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
// __________________________________________________________________________
//
//    begin                :
//    email                : developers@lammps.org
// ***************************************************************************

#if defined(NV_KERNEL) || defined(USE_HIP)
#include "lal_aux_fun1.h"
#ifndef _DOUBLE_DOUBLE
_texture( pos_tex,float4);
_texture( q_tex,float);
#else
_texture_2d( pos_tex,int4);
_texture( q_tex,int2);
#endif
#else
#define pos_tex x_
#define q_tex q_
#endif

// Index into the cs/sn arrays.  The arrays store the cosine and sine of
// m*unitk[ic]*r for m in [-kmax,kmax] (full range, negatives included), all
// three Cartesian directions ic, and all local atoms i.  The atom index is
// the fastest-varying dimension so reads in the reduction kernel are coalesced.
#define CS_INDEX(m,ic,i,nlocal,kmax) ((((m)+(kmax))*3+(ic))*(nlocal)+(i))

// ---------------------------------------------------------------------------
// Compute cos(m*unitk[ic]*r_i) and sin(...) for every atom, direction, and
// |m| up to kmax using the standard cosine/sine recurrence.  One thread per
// atom.  Negative m are stored using cos(-m)=cos(m), sin(-m)=-sin(m).
// ---------------------------------------------------------------------------
__kernel void k_ewald_cssn(const __global numtyp4 *restrict x_,
                           __global numtyp *restrict cs,
                           __global numtyp *restrict sn,
                           const numtyp unitkx, const numtyp unitky,
                           const numtyp unitkz, const int kmax,
                           const int nlocal) {
  int i=GLOBAL_ID_X;
  if (i<nlocal) {
    numtyp4 p;
    fetch4(p,i,pos_tex);
    numtyp xx[3], unitk[3];
    xx[0]=p.x; xx[1]=p.y; xx[2]=p.z;
    unitk[0]=unitkx; unitk[1]=unitky; unitk[2]=unitkz;

    for (int ic=0; ic<3; ic++) {
      const numtyp arg=unitk[ic]*xx[ic];
      const numtyp c1=cos(arg);
      const numtyp s1=sin(arg);

      // m = 0
      cs[CS_INDEX(0,ic,i,nlocal,kmax)]=(numtyp)1.0;
      sn[CS_INDEX(0,ic,i,nlocal,kmax)]=(numtyp)0.0;

      // m = 1 and m = -1
      cs[CS_INDEX(1,ic,i,nlocal,kmax)]=c1;
      sn[CS_INDEX(1,ic,i,nlocal,kmax)]=s1;
      cs[CS_INDEX(-1,ic,i,nlocal,kmax)]=c1;
      sn[CS_INDEX(-1,ic,i,nlocal,kmax)]=-s1;

      numtyp csm1=c1, snm1=s1;
      for (int m=2; m<=kmax; m++) {
        const numtyp csm=csm1*c1-snm1*s1;
        const numtyp snm=snm1*c1+csm1*s1;
        cs[CS_INDEX(m,ic,i,nlocal,kmax)]=csm;
        sn[CS_INDEX(m,ic,i,nlocal,kmax)]=snm;
        cs[CS_INDEX(-m,ic,i,nlocal,kmax)]=csm;
        sn[CS_INDEX(-m,ic,i,nlocal,kmax)]=-snm;
        csm1=csm;
        snm1=snm;
      }
    }
  }
}

// ---------------------------------------------------------------------------
// Compute the local (per-MPI-rank) structure factors
//   sfacrl[k] = sum_i q_i cos(k.r_i),  sfacim[k] = sum_i q_i sin(k.r_i)
// One thread-block per k-vector reduces over all local atoms, so the result
// is deterministic (fixed-order tree reduction, no floating-point atomics).
// ---------------------------------------------------------------------------
__kernel void k_ewald_structure(const __global numtyp *restrict q_,
                                const __global numtyp *restrict cs,
                                const __global numtyp *restrict sn,
                                const __global int *restrict kxvecs,
                                const __global int *restrict kyvecs,
                                const __global int *restrict kzvecs,
                                __global acctyp *restrict sfacrl,
                                __global acctyp *restrict sfacim,
                                const int kmax, const int nlocal,
                                const int kcount) {
  __local acctyp red[2][BLOCK_PAIR];

  const int tid=THREAD_ID_X;
  const int k=BLOCK_ID_X;

  acctyp sr=(acctyp)0.0;
  acctyp si=(acctyp)0.0;

  if (k<kcount) {
    const int kx=kxvecs[k];
    const int ky=kyvecs[k];
    const int kz=kzvecs[k];
    const int basex=CS_INDEX(kx,0,0,nlocal,kmax);
    const int basey=CS_INDEX(ky,1,0,nlocal,kmax);
    const int basez=CS_INDEX(kz,2,0,nlocal,kmax);

    for (int i=tid; i<nlocal; i+=BLOCK_SIZE_X) {
      const numtyp csx=cs[basex+i];
      const numtyp snx=sn[basex+i];
      const numtyp csy=cs[basey+i];
      const numtyp sny=sn[basey+i];
      const numtyp csz=cs[basez+i];
      const numtyp snz=sn[basez+i];
      const numtyp cypz=csy*csz-sny*snz;
      const numtyp sypz=sny*csz+csy*snz;
      const numtyp exprl=csx*cypz-snx*sypz;
      const numtyp expim=snx*cypz+csx*sypz;
      numtyp qi;
      fetch(qi,i,q_tex);
      sr+=(acctyp)(qi*exprl);
      si+=(acctyp)(qi*expim);
    }
  }

  block_reduce_add2(simd_size(),red,tid,sr,si);

  if (tid==0 && k<kcount) {
    sfacrl[k]=sr;
    sfacim[k]=si;
  }
}
