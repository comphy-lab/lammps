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
/* ------------------------------------------------------
    This file is part of the LDD package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    and Maria C. Lesniewski, mjl6766@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */

#ifdef DUMP_CLASS
// clang-format off
DumpStyle(ldd,DumpLdd);
// clang-format on
#else

#ifndef LMP_DUMP_LDD_H
#define LMP_DUMP_LDD_H

#include "dump.h"

namespace LAMMPS_NS {

class DumpLdd : public Dump {
 public:
  DumpLdd(LAMMPS *, int, char**);

 protected:
  int scale_flag;            // 1 if atom coords are scaled, 0 if no
  int image_flag;            // 1 if append box count to atom coords, 0 if no

  char *columns;             // column labels

  void init_style() override;
  int modify_param(int, char **) override;
  void write_header(bigint) override;
  void pack(tagint *) override;
  int convert_string(int, double *) override;
  void write_data(int, double *) override;

  typedef void (DumpLdd::*FnPtrHeader)(bigint);
  FnPtrHeader header_choice;           // ptr to write header functions
  void header_binary(bigint);
  void header_binary_triclinic(bigint);
  void header_item(bigint);
  void header_item_triclinic(bigint);

  typedef void (DumpLdd::*FnPtrPack)(tagint *);
  FnPtrPack pack_choice;               // ptr to pack functions
  void pack_scale_image(tagint *);
  void pack_scale_noimage(tagint *);
  void pack_noscale_image(tagint *);
  void pack_noscale_noimage(tagint *);
  void pack_scale_image_triclinic(tagint *);
  void pack_scale_noimage_triclinic(tagint *);

  typedef int (DumpLdd::*FnPtrConvert)(int, double *);
  FnPtrConvert convert_choice;          // ptr to convert data functions
  int convert_image(int, double *);
  int convert_noimage(int, double *);

  typedef void (DumpLdd::*FnPtrWrite)(int, double *);
  FnPtrWrite write_choice;              // ptr to write data functions
  void write_binary(int, double *);
  void write_string(int, double *);
  void write_lines_image(int, double *);
  void write_lines_noimage(int, double *);
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
