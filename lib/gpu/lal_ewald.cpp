/***************************************************************************
                                lal_ewald.cpp
                             -------------------
                            W. Michael Brown (ORNL)
                            Axel Kohlmeyer (Temple)

  Class for Ewald (reciprocal-space) acceleration

 __________________________________________________________________________
    This file is part of the LAMMPS Accelerator Library (LAMMPS_AL)
 __________________________________________________________________________

    begin                :
    email                : developers@lammps.org
 ***************************************************************************/

#if defined(USE_OPENCL)
#include "ewald_cl.h"
#elif defined(USE_CUDART)
const char *ewald=0;
#else
#include "ewald_cubin.h"
#endif

#include "lal_ewald.h"
#include <cassert>
#include <cmath>

namespace LAMMPS_AL {
#define EwaldGPUT EwaldGPU<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> global_device;

template <class numtyp, class acctyp>
EwaldGPUT::EwaldGPU() : _allocated(false), _compiled(false), _max_bytes(0) {
  device=&global_device;
  ans=new Answer<numtyp,acctyp>();
  ewald_program=nullptr;
  ucl_device=nullptr;
  _kmax=0;
  _kcount=0;
  _cs_kmax=-1;
  _cs_nlocal=0;
}

template <class numtyp, class acctyp>
EwaldGPUT::~EwaldGPU() {
  clear(0.0);
  k_cssn.clear();
  k_structure.clear();
  if (ewald_program) delete ewald_program;
  delete ans;
}

template <class numtyp, class acctyp>
int EwaldGPUT::init(const int nlocal, const int nall, FILE *_screen,
                    int &flag) {
  _max_bytes=10;
  screen=_screen;

  flag=device->init(*ans,nlocal,nall);
  if (flag!=0)
    return flag;
  if (device->ptx_arch()>0.0 && device->ptx_arch()<1.1) {
    flag=-4;
    return flag;
  }

  if (ucl_device!=device->gpu) _compiled=false;
  ucl_device=device->gpu;
  atom=&device->atom;

  _block_size=device->pair_block_size();

  compile_kernels(*ucl_device);

  pos_tex.bind_float(atom->x,4);
  q_tex.bind_float(atom->q,1);

  _allocated=true;
  _max_bytes=0;
  _max_an_bytes=ans->gpu_bytes();
  _cs_kmax=-1;
  _cs_nlocal=0;
  _kmax=0;
  _kcount=0;
  return flag;
}

// ---------------------------------------------------------------------------
// Upload the (constant per box) k-vectors and grid parameters
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EwaldGPUT::setup(const int kmax, const int kcount, int *kxvecs,
                      int *kyvecs, int *kzvecs, double *unitk, bool &success) {
  _kmax=kmax;
  _kcount=kcount;
  _unitk[0]=(numtyp)unitk[0];
  _unitk[1]=(numtyp)unitk[1];
  _unitk[2]=(numtyp)unitk[2];

  // (re)allocate and upload the integer k-vectors (length kcount)

  d_kxvecs.clear();
  d_kyvecs.clear();
  d_kzvecs.clear();
  success = success && (d_kxvecs.alloc(kcount,*ucl_device,UCL_READ_ONLY)==
                        UCL_SUCCESS);
  success = success && (d_kyvecs.alloc(kcount,*ucl_device,UCL_READ_ONLY)==
                        UCL_SUCCESS);
  success = success && (d_kzvecs.alloc(kcount,*ucl_device,UCL_READ_ONLY)==
                        UCL_SUCCESS);

  // (re)allocate the structure-factor buffers (host+device, length kcount)

  d_sfacrl.clear();
  d_sfacim.clear();
  success = success && (d_sfacrl.alloc(kcount,*ucl_device,UCL_WRITE_ONLY,
                                       UCL_READ_WRITE)==UCL_SUCCESS);
  success = success && (d_sfacim.alloc(kcount,*ucl_device,UCL_WRITE_ONLY,
                                       UCL_READ_WRITE)==UCL_SUCCESS);
  if (!success)
    return;

  UCL_H_Vec<int> view;
  view.view(kxvecs,kcount,*ucl_device);
  ucl_copy(d_kxvecs,view,false);
  view.view(kyvecs,kcount,*ucl_device);
  ucl_copy(d_kyvecs,view,false);
  view.view(kzvecs,kcount,*ucl_device);
  ucl_copy(d_kzvecs,view,false);

  // force reallocation of the per-atom cs/sn arrays (kmax may have changed)

  _cs_kmax=-1;
  _cs_nlocal=0;
}

// ---------------------------------------------------------------------------
// (Re)allocate the per-atom cos/sin arrays for the current kmax and nlocal
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
void EwaldGPUT::resize_cssn(const int nlocal, bool &success) {
  if (nlocal<=_cs_nlocal && _kmax==_cs_kmax)
    return;

  int newn=static_cast<int>(static_cast<double>(nlocal)*1.10)+1;
  int sz=(2*_kmax+1)*3*newn;

  d_cs.clear();
  d_sn.clear();
  success = success && (d_cs.alloc(sz,*ucl_device)==UCL_SUCCESS);
  success = success && (d_sn.alloc(sz,*ucl_device)==UCL_SUCCESS);
  if (!success)
    return;

  _cs_nlocal=newn;
  _cs_kmax=_kmax;
}

// ---------------------------------------------------------------------------
// Compute the local (per-rank) structure factors on the device
// ---------------------------------------------------------------------------
template <class numtyp, class acctyp>
int EwaldGPUT::structure(const int ago, const int nlocal, const int nall,
                         double **host_x, int *host_type, double *host_q,
                         double *host_sfacrl, double *host_sfacim,
                         bool &success) {
  if (nlocal==0)
    return 0;

  _nlocal=nlocal;
  ans->inum(nlocal);

  if (ago==0) {
    resize_atom(nlocal,nall,success);
    if (!success)
      return 0;
  }
  resize_cssn(nlocal,success);
  if (!success)
    return 0;

  atom->cast_x_data(host_x,host_type);
  atom->cast_q_data(host_q);
  atom->add_x_data(host_x,host_type);
  atom->add_q_data();

  const int BX=_block_size;

  // cos/sin recurrence: one thread per atom
  const int GX=static_cast<int>(ceil(static_cast<double>(nlocal)/BX));
  k_cssn.set_size(GX,BX);
  k_cssn.run(&atom->x, &d_cs, &d_sn, &_unitk[0], &_unitk[1], &_unitk[2],
             &_kmax, &_nlocal);

  // structure factor reduction: one block per k-vector
  k_structure.set_size(_kcount,BX);
  k_structure.run(&atom->q, &d_cs, &d_sn, &d_kxvecs, &d_kyvecs, &d_kzvecs,
                  &d_sfacrl, &d_sfacim, &_kmax, &_nlocal, &_kcount);

  d_sfacrl.update_host(_kcount,false);
  d_sfacim.update_host(_kcount,false);
  for (int k=0; k<_kcount; k++) {
    host_sfacrl[k]=d_sfacrl[k];
    host_sfacim[k]=d_sfacim[k];
  }

  return 0;
}

template <class numtyp, class acctyp>
void EwaldGPUT::clear(const double /*cpu_time*/) {
  if (!_allocated)
    return;
  _allocated=false;

  d_kxvecs.clear();
  d_kyvecs.clear();
  d_kzvecs.clear();
  d_cs.clear();
  d_sn.clear();
  d_sfacrl.clear();
  d_sfacim.clear();

  ans->clear();
}

template <class numtyp, class acctyp>
double EwaldGPUT::host_memory_usage() const {
  return device->atom.host_memory_usage()+
         sizeof(EwaldGPU<numtyp,acctyp>);
}

template <class numtyp, class acctyp>
void EwaldGPUT::compile_kernels(UCL_Device &dev) {
  if (_compiled)
    return;

  std::string flags=device->compile_string();

  if (ewald_program) delete ewald_program;
  ewald_program=new UCL_Program(dev);
  ewald_program->load_string(ewald,flags.c_str(),nullptr,screen);

  k_cssn.set_function(*ewald_program,"k_ewald_cssn");
  k_structure.set_function(*ewald_program,"k_ewald_structure");
  pos_tex.get_texture(*ewald_program,"pos_tex");
  q_tex.get_texture(*ewald_program,"q_tex");

  _compiled=true;
}

template class EwaldGPU<PRECISION,ACC_PRECISION>;
}
