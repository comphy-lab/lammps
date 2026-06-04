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

#include "lal_ewald.h"
#include <cassert>

namespace LAMMPS_AL {
#define EwaldGPUT EwaldGPU<numtyp, acctyp>

extern Device<PRECISION,ACC_PRECISION> global_device;

template <class numtyp, class acctyp>
EwaldGPUT::EwaldGPU() : _allocated(false), _max_bytes(0) {
  device=&global_device;
  ans=new Answer<numtyp,acctyp>();
}

template <class numtyp, class acctyp>
EwaldGPUT::~EwaldGPU() {
  clear(0.0);
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

  ucl_device=device->gpu;
  atom=&device->atom;

  _allocated=true;
  _max_bytes=0;
  _max_an_bytes=ans->gpu_bytes();
  return flag;
}

template <class numtyp, class acctyp>
void EwaldGPUT::clear(const double /*cpu_time*/) {
  if (!_allocated)
    return;
  _allocated=false;

  ans->clear();
}

template <class numtyp, class acctyp>
double EwaldGPUT::host_memory_usage() const {
  return device->atom.host_memory_usage()+
         sizeof(EwaldGPU<numtyp,acctyp>);
}

template class EwaldGPU<PRECISION,ACC_PRECISION>;
}
