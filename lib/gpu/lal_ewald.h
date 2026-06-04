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

#if defined(USE_OPENCL)
#include "geryon/ocl_texture.h"
#elif defined(USE_CUDART)
#include "geryon/nvc_texture.h"
#elif defined(USE_HIP)
#include "geryon/hip_texture.h"
#else
#include "geryon/nvd_texture.h"
#endif

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

  /// Upload the (constant per box) k-vectors and store grid parameters.
  /** Called from the kspace setup() so it is refreshed whenever the box
    * changes.  (Re)allocates the k-vector and structure-factor buffers. **/
  void setup(const int kmax, const int kcount, int *kxvecs, int *kyvecs,
             int *kzvecs, double **eg, double *unitk, bool &success);

  /// Compute the local (per-rank) structure factors on the device.
  /** Uploads x/q (gated on ago), runs the cs/sn and structure-factor
    * kernels, and copies the kcount structure factors back to the host. **/
  int structure(const int ago, const int nlocal, const int nall,
                double **host_x, int *host_type, double *host_q,
                double *host_sfacrl, double *host_sfacim, bool &success);

  /// K-space field/force: upload the global structure factors and run the
  /// field kernel, reusing the resident cs/sn.  The k-space force is queued
  /// for merging with the pair force by fix gpu.
  void compute_forces(double *host_sfacrl_all, double *host_sfacim_all,
                      const double qscale, const int slabflag, bool &success);

  /// Check if there is enough storage for atom arrays and realloc if not
  inline void resize_atom(const int inum, const int nall, bool &success) {
    if (atom->resize(nall, success)) {
      pos_tex.bind_float(atom->x,4);
      q_tex.bind_float(atom->q,1);
    }
    ans->resize(inum,success);
  }

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

  // ----------------------- RECIPROCAL-SPACE DATA --------------------

  /// k-vector integer indices (length kcount), constant per box
  UCL_D_Vec<int> d_kxvecs, d_kyvecs, d_kzvecs;

  /// per-k field coefficients eg[k][0..2] packed into x,y,z (constant per box)
  UCL_D_Vec<numtyp4> d_eg;

  /// cos/sin of m*unitk[ic]*r for m in [-kmax,kmax], 3 dirs, all local atoms
  UCL_D_Vec<numtyp> d_cs, d_sn;

  /// structure factors: local after k_structure, global (sfac_all) for k_field
  UCL_Vector<acctyp,acctyp> d_sfacrl, d_sfacim;

  // ------------------------- DEVICE KERNELS -------------------------

  UCL_Program *ewald_program;
  UCL_Kernel k_cssn, k_structure, k_field;

  // --------------------------- TEXTURES -----------------------------

  UCL_Texture pos_tex;
  UCL_Texture q_tex;

 protected:
  bool _allocated, _compiled;
  int _block_size;
  double _max_bytes, _max_an_bytes;

  // grid parameters
  int _kmax, _kcount, _nlocal;
  int _cs_kmax, _cs_nlocal;
  numtyp _unitk[3];

  // per-step field/force parameters (kept as members for stable kernel args)
  numtyp _qscale;
  int _slabflag;

  void compile_kernels(UCL_Device &dev);
  void resize_cssn(const int nlocal, bool &success);
};

}

#endif
