/***************************************************************************
                                 lal_ewald.h
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

#ifndef LAL_EWALD_H
#define LAL_EWALD_H

#include "mpi.h"
#include "lal_device.h"

namespace LAMMPS_AL {

template <class numtyp, class acctyp> class Device;

template <class numtyp, class acctyp>
class EwaldGPU {
 public:
  EwaldGPU();
  virtual ~EwaldGPU();

  /// Clear any previous data and set up for a new LAMMPS run
  /** \param success set to the device init error code
    *   -  0 if successful
    *   - -1 if fix gpu not found
    *   - -2 if GPU could not be found
    *   - -3 if there is an out of memory error
    *   - -4 if the GPU library was not compiled for GPU
    *   - -5 if double precision is not supported on card **/
  int init(const int nlocal, const int nall, FILE *screen, int &success);

  /// Clear all host and device data
  /** \note This is called at the beginning of the init() routine **/
  void clear(const double cpu_time);

  /// Total host memory used by library
  double host_memory_usage() const;

  // -------------------------- DEVICE DATA -------------------------

  /// Device Properties and Atom and Neighbor storage
  Device<numtyp,acctyp> *device;

  /// Geryon device
  UCL_Device *ucl_device;

  /// LAMMPS pointer for screen output
  FILE *screen;

  // --------------------------- ATOM DATA --------------------------

  /// Atom Data
  Atom<numtyp,acctyp> *atom;

  // ------------------------ FORCE/ENERGY DATA -----------------------

  Answer<numtyp,acctyp> *ans;

 protected:
  bool _allocated;
  double _max_bytes, _max_an_bytes;
};

}

#endif
