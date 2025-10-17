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
/* ------------------------------------------------------
    This file is part of the LDD package for LAMMPS.
    Contributed by Michael R. DeLyser, mrd5285@psu.edu
    and Maria C. Lesniewski, mjl6766@psu.edu
    The Pennsylvania State University
   ------------------------------------------------------ */

#include "dump_ldd.h"

#include <cstring>
#include <cstdio>
#include <cstdlib>

#include "atom.h"
#include "domain.h"
#include "error.h"
#include "group.h"
#include "memory.h"
#include "update.h"

using namespace LAMMPS_NS;

#define ONELINE 256
#define DELTA 1048576

/* ---------------------------------------------------------------------- */

DumpLdd::DumpLdd(LAMMPS *lmp, int narg, char **arg) : Dump(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal dump ldd command");

  scale_flag = 0;
  image_flag = 0;
  buffer_allow = 1;
  buffer_flag = 1;
  format_default = NULL;
}

/* ---------------------------------------------------------------------- */

void DumpLdd::init_style()
{
int chkidx = 0;
  // id type x y z vx vy vz fx fy fz ldds ldnrgs grads gradnrgs ldttlnrg
  int ntypes = atom->ntypes;
  size_one = 12 + 2 * ntypes + 4 * ntypes;
 // 12 = id type x y z vx vy vz fx fy fz ldttlnrg
  // 2 * ntypes = ldds ldnrgs
  // 4 * ntypes = gradx grady gradz gradnrg
  // default format depends on image flags

   if (atom->molecule_flag == 1) { size_one++; }
  delete [] format;

  // 3 chars per value: the % sign, a 'd' or a 'g', and a space
  // we have 12 + 6 * ntypes values = 36 + 18 * ntypes characters.
  // 200 chars works for up to 9 atom types (36 + 18 * 9 = 162 + 36 = 198)
  char *str, str0[200];
  int str_length;
  if (ntypes <= 8) { str = &(str0[0]); str_length = 200; }
  else { memory->create(str,40+18*ntypes,"DumpLDD:str"); str_length = 40+18*ntypes; }
  char * tmpstr;
  memory->create(tmpstr, str_length, "DumpLDD:tmpstr"); //MCL 05.22.25 - Gotta create seperate addresses for src/dest for better defined bhvr
  for (int i = 0; i < str_length; i++)
  {
          str[i] = '\0';
          tmpstr[i] = '\0';
  }

// start with spots for id, (mol,) type, x, y, z, vx, vy, vz, fx, fy, fz
  if (atom->molecule_flag == 1) { sprintf(str,"%s","%d %d %d %g %g %g %g %g %g %g %g %g"); }
  else { sprintf(str,"%s","%d %d %d %g %g %g %g %g %g %g %g %g");}
// for each type, add spots for ldd and ldnrg
  for (int i = 1; i <= ntypes; i++)
  {
          sprintf(tmpstr,"%s %s",str,"%g %g");
          strcpy(str,tmpstr);
  }
// for each type, add spots for gradx, grady, gradz, and gradnrg
  for (int i = 1; i <= ntypes; i++)
  {
          sprintf(tmpstr,"%s %s",str,"%g %g %g %g");
          strcpy(str,tmpstr);
  }
//lastly, add spot for total energy
  sprintf(tmpstr,"%s %s",str,"%g\n");
  strcpy(str,tmpstr);

  format = new char[strlen(str) + 1];
  strcpy(format,str);
  // setup boundary string

  domain->boundary_string(boundstr);

  // setup column string

  // columns = (char *) "id type x y z fx fy fz ";
  // 32 characters for "id type x y z vx vy vz fx fy fz "
  // 9 more for "lddttlnrg"
  // for each type, we need 45 for:
  // "lddens# ldnrg# gradx# grady# gradz# gradnrg# "
  // and that's if each # is 1 digit each.
  // let's make it 57 chars, which is enough for each # to be 3 digits.

  int collength = 41 + 57 * ntypes;

  if (atom->molecule_flag == 1) { collength += 4; }

  char *cols;
  char *tmpcols;
  memory->create(cols,collength,"DumpLDD:cols");
  memory->create(tmpcols, collength, "DumpLDD:tmpcols");
  char ncol[60];
  for (int i=0; i<60; ++i) { ncol[i] = '\0'; }
  for (int i=0; i<collength; ++i) { cols[i] = '\0'; tmpcols[i] = '\0';}
  if (atom->molecule_flag == 1)
  {
          sprintf(cols,"%s","id mol type x y z vx vy vz fx fy fz");
          sprintf(tmpcols,"%s", "id mol type x y z vx vy vz fx fy fz");
  }
  else
  {
          sprintf(cols,"%s","id type x y z vx vy vz fx fy fz");
          sprintf(tmpcols,"%s","id type x y z vx vy vz fx fy fz");
  }
  for (int i=1; i<=ntypes; i++)
  {
    sprintf(ncol,"lddens%d ldnrg%d",i,i);
    sprintf(tmpcols,"%s %s",cols,ncol);
    strcpy(cols,tmpcols);
  }
  for (int i=1; i<=ntypes; ++i)
  {
    sprintf(ncol,"gradx%d grady%d gradz%d gradnrg%d",i,i,i,i);
    sprintf(tmpcols,"%s %s",cols,ncol);
    strcpy(cols,tmpcols);
  }
  sprintf(tmpcols,"%s %s",cols,"lddttlnrg");
  strcpy(cols,tmpcols);
  columns = (char *) cols;

  // setup function ptrs

  header_choice = &DumpLdd::header_item;
  pack_choice = &DumpLdd::pack_noscale_noimage;
  convert_choice = &DumpLdd::convert_noimage;
  write_choice = &DumpLdd::write_string;
  // open single file, one time only

  if (multifile == 0) openfile();

}

/* ---------------------------------------------------------------------- */

int DumpLdd::modify_param(int narg, char **arg)
{
  if (strcmp(arg[0],"scale") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"yes") == 0) scale_flag = 1;
    else if (strcmp(arg[1],"no") == 0) scale_flag = 0;
    else error->all(FLERR,"Illegal dump_modify command");
    return 2;
  } else if (strcmp(arg[0],"image") == 0) {
    if (narg < 2) error->all(FLERR,"Illegal dump_modify command");
    if (strcmp(arg[1],"yes") == 0) image_flag = 1;
    else if (strcmp(arg[1],"no") == 0) image_flag = 0;
    else error->all(FLERR,"Illegal dump_modify command");
    return 2;
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

void DumpLdd::write_header(bigint ndump)
{
  if (multiproc) (this->*header_choice)(ndump);
  else if (me == 0) (this->*header_choice)(ndump);
}

/* ---------------------------------------------------------------------- */

void DumpLdd::pack(tagint *ids)
{
  (this->*pack_choice)(ids);
}

/* ---------------------------------------------------------------------- */

int DumpLdd::convert_string(int n, double *mybuf)
{
  return (this->*convert_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpLdd::write_data(int n, double *mybuf)
{
  (this->*write_choice)(n,mybuf);
}

/* ---------------------------------------------------------------------- */

void DumpLdd::header_binary(bigint ndump)
{
  fwrite(&update->ntimestep,sizeof(bigint),1,fp);
  fwrite(&ndump,sizeof(bigint),1,fp);
  fwrite(&domain->triclinic,sizeof(int),1,fp);
  fwrite(&domain->boundary[0][0],6*sizeof(int),1,fp);
  fwrite(&boxxlo,sizeof(double),1,fp);
  fwrite(&boxxhi,sizeof(double),1,fp);
  fwrite(&boxylo,sizeof(double),1,fp);
  fwrite(&boxyhi,sizeof(double),1,fp);
  fwrite(&boxzlo,sizeof(double),1,fp);
  fwrite(&boxzhi,sizeof(double),1,fp);
  fwrite(&size_one,sizeof(int),1,fp);
  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpLdd::header_binary_triclinic(bigint ndump)
{
  fwrite(&update->ntimestep,sizeof(bigint),1,fp);
  fwrite(&ndump,sizeof(bigint),1,fp);
  fwrite(&domain->triclinic,sizeof(int),1,fp);
  fwrite(&domain->boundary[0][0],6*sizeof(int),1,fp);
  fwrite(&boxxlo,sizeof(double),1,fp);
  fwrite(&boxxhi,sizeof(double),1,fp);
  fwrite(&boxylo,sizeof(double),1,fp);
  fwrite(&boxyhi,sizeof(double),1,fp);
  fwrite(&boxzlo,sizeof(double),1,fp);
  fwrite(&boxzhi,sizeof(double),1,fp);
  fwrite(&boxxy,sizeof(double),1,fp);
  fwrite(&boxxz,sizeof(double),1,fp);
  fwrite(&boxyz,sizeof(double),1,fp);
  fwrite(&size_one,sizeof(int),1,fp);
  if (multiproc) fwrite(&nclusterprocs,sizeof(int),1,fp);
  else fwrite(&nprocs,sizeof(int),1,fp);
}

/* ---------------------------------------------------------------------- */

void DumpLdd::header_item(bigint ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS %s\n",boundstr);
  fprintf(fp,"%g %g\n",boxxlo,boxxhi);
  fprintf(fp,"%g %g\n",boxylo,boxyhi);
  fprintf(fp,"%g %g\n",boxzlo,boxzhi);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpLdd::header_item_triclinic(bigint ndump)
{
  fprintf(fp,"ITEM: TIMESTEP\n");
  fprintf(fp,BIGINT_FORMAT "\n",update->ntimestep);
  fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  fprintf(fp,BIGINT_FORMAT "\n",ndump);
  fprintf(fp,"ITEM: BOX BOUNDS xy xz yz %s\n",boundstr);
  fprintf(fp,"%g %g %g\n",boxxlo,boxxhi,boxxy);
  fprintf(fp,"%g %g %g\n",boxylo,boxyhi,boxxz);
  fprintf(fp,"%g %g %g\n",boxzlo,boxzhi,boxyz);
  fprintf(fp,"ITEM: ATOMS %s\n",columns);
}

/* ---------------------------------------------------------------------- */

void DumpLdd::pack_scale_image(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  imageint *image = atom->image;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double invxprd = 1.0/domain->xprd;
  double invyprd = 1.0/domain->yprd;
  double invzprd = 1.0/domain->zprd;

  m = n = 00;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = (x[i][0] - boxxlo) * invxprd;
      buf[m++] = (x[i][1] - boxylo) * invyprd;
      buf[m++] = (x[i][2] - boxzlo) * invzprd;
      buf[m++] = (image[i] & IMGMASK) - IMGMAX;
      buf[m++] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      buf[m++] = (image[i] >> IMG2BITS) - IMGMAX;
      if (ids) ids[n++] = tag[i];
    }
}

/* ---------------------------------------------------------------------- */

void DumpLdd::pack_scale_noimage(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double invxprd = 1.0/domain->xprd;
  double invyprd = 1.0/domain->yprd;
  double invzprd = 1.0/domain->zprd;

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = (x[i][0] - boxxlo) * invxprd;
      buf[m++] = (x[i][1] - boxylo) * invyprd;
      buf[m++] = (x[i][2] - boxzlo) * invzprd;
      if (ids) ids[n++] = tag[i];
    }
}

/* ---------------------------------------------------------------------- */

void DumpLdd::pack_scale_image_triclinic(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  imageint *image = atom->image;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double lamda[3];

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      domain->x2lamda(x[i],lamda);
      buf[m++] = lamda[0];
      buf[m++] = lamda[1];
      buf[m++] = lamda[2];
      buf[m++] = (image[i] & IMGMASK) - IMGMAX;
      buf[m++] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      buf[m++] = (image[i] >> IMG2BITS) - IMGMAX;
      if (ids) ids[n++] = tag[i];
    }
}

/* ---------------------------------------------------------------------- */

void DumpLdd::pack_scale_noimage_triclinic(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  double lamda[3];

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      domain->x2lamda(x[i],lamda);
      buf[m++] = lamda[0];
      buf[m++] = lamda[1];
      buf[m++] = lamda[2];
      if (ids) ids[n++] = tag[i];
    }
}

/* ---------------------------------------------------------------------- */

void DumpLdd::pack_noscale_image(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  int *type = atom->type;
  imageint *image = atom->image;
  int *mask = atom->mask;
  double **x = atom->x;
  int nlocal = atom->nlocal;

  m = n = 0;
  for (int i = 0; i < nlocal; i++)
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      buf[m++] = (image[i] & IMGMASK) - IMGMAX;
      buf[m++] = (image[i] >> IMGBITS & IMGMASK) - IMGMAX;
      buf[m++] = (image[i] >> IMG2BITS) - IMGMAX;
      if (ids) ids[n++] = tag[i];
    }
}

/* ---------------------------------------------------------------------- */

void DumpLdd::pack_noscale_noimage(tagint *ids)
{
  int m,n;

  tagint *tag = atom->tag;
  int *mol = atom->molecule;
  int *type = atom->type;
  int *mask = atom->mask;
  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int nlocal = atom->nlocal;
  int ntypes = atom->ntypes;
  double **LDs = atom->ldd_local_density;
  double **LDEs = atom->ldd_energy;
  double **LDGrads = atom->ldd_grad_density;
  double **LDGradEs = atom->ldd_grad_energy;
  double *LDTE = atom->ldd_total_energy;

  m = n = 0;

  for (int i = 0; i < nlocal; i++)
  {
    if (mask[i] & groupbit) {
      buf[m++] = tag[i];
      if (atom->molecule_flag == 1) { buf[m++] = mol[i]; }
      buf[m++] = type[i];
      buf[m++] = x[i][0];
      buf[m++] = x[i][1];
      buf[m++] = x[i][2];
      buf[m++] = v[i][0];
      buf[m++] = v[i][1];
      buf[m++] = v[i][2];
      buf[m++] = f[i][0];
      buf[m++] = f[i][1];
      buf[m++] = f[i][2];
      for (int j = 1; j <=ntypes; j++)
      {
        buf[m++] = LDs[i][j];
        buf[m++] = LDEs[i][j];
      }
      for (int j = 1; j <=ntypes; ++j)
      {
        buf[m++] = LDGrads[i][3*j];
        buf[m++] = LDGrads[i][3*j+1];
        buf[m++] = LDGrads[i][3*j+2];
        buf[m++] = LDGradEs[i][j];
      }
      buf[m++] = LDTE[i];
      if (ids) ids[n++] = tag[i];
    }
  }
}

/* ----------------------------------------------------------------------
   convert mybuf of doubles to one big formatted string in sbuf
   return -1 if strlen exceeds an int, since used as arg in MPI calls in Dump
------------------------------------------------------------------------- */

int DumpLdd::convert_image(int n, double *mybuf)
{
  int offset = 0;
  int m = 0;

  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }

    offset += sprintf(&sbuf[offset],format,
                      static_cast<tagint> (mybuf[m]),
                      static_cast<int> (mybuf[m+1]),
                      mybuf[m+2],mybuf[m+3],mybuf[m+4],
                      static_cast<int> (mybuf[m+5]),
                      static_cast<int> (mybuf[m+6]),
                      static_cast<int> (mybuf[m+7]));
    m += size_one;
  }

  return offset;
}

/* ---------------------------------------------------------------------- */

int DumpLdd::convert_noimage(int n, double *mybuf)
{
  int offset = 0;
  int m = 0;
  int ntypes = atom->ntypes;

  // strentries = (char *) "id (mol) type x y z fx fy fz ";
  // 16 characters for id/type/(mol) if we allow 3 digs and spaces
  // 30 characters for xyz positions up to 999.999999 with spaces between (9+1)
  // 36x2=72 characters for v/f entries to allow negatives and spaces in sets of three
  //
  // 39xntypes for gradxyz entries to allow negatives spaces and scientific notations
  // 24xntypes for lddensn/lddensenergy
  // 13xntypes for gradnrgn
  // 12 more for "lddttlnrg"
  // So thats a generous 118(base) + 76*ntypes+ 12(totalld) +1 (terminating)
  int entry_len = 131 + 76*ntypes;

//  char * str = (char *) calloc(800,sizeof(char)); //MCL NAIVE GUESS
   char * str;
   char * tmpstr;
   memory->create(str,entry_len,"DumpLDD:str");
   memory->create(tmpstr,entry_len,"DumpLDD:tmpstr");

  for (int i = 0; i < n; i++) {
    if (offset + ONELINE > maxsbuf) {
      if ((bigint) maxsbuf + DELTA > MAXSMALLINT) return -1;
      maxsbuf += DELTA;
      memory->grow(sbuf,maxsbuf,"dump:sbuf");
    }
    if (atom->molecule_flag == 1)
    {
            sprintf(str,
                    "%d %d %d %g %g %g %g %g %g %g %g %g",
                    static_cast<tagint> (mybuf[m]),
                    static_cast<int> (mybuf[m+1]),
                    static_cast<int> (mybuf[m+2]),
                    mybuf[m+3],
                    mybuf[m+4],
                    mybuf[m+5],
                    mybuf[m+6],
                    mybuf[m+7],
                    mybuf[m+8],
                    mybuf[m+9],
                    mybuf[m+10],
                    mybuf[m+11]);
            m += 11;
            strcpy(tmpstr, str);
    }
    else
    {
            sprintf(str,
                    "%d %d %g %g %g %g %g %g %g %g %g",
                    static_cast<tagint> (mybuf[m]),
                    static_cast<int> (mybuf[m+1]),
                    mybuf[m+2],
                    mybuf[m+3],
                    mybuf[m+4],
                    mybuf[m+5],
                    mybuf[m+6],
                    mybuf[m+7],
                    mybuf[m+8],
                    mybuf[m+9],
                    mybuf[m+10]);
            m += 10;
            strcpy(tmpstr,str);
    }
    for (int j = 1; j <= 2*ntypes; j+=2) // LDs and LDnrgs
    {
      sprintf(tmpstr,"%s %g %g",str,mybuf[m+j],mybuf[m+j+1]);
      strcpy(str,tmpstr);
    }
    m += 2*ntypes;
    for (int j = 1; j <= 4*ntypes; j++) // LDGradx LDGrady LDGradz and LDGradNRG
    {
      sprintf(tmpstr,"%s %g",str,mybuf[m+j]);
      strcpy(str,tmpstr);
    }
    m += 4*ntypes;
    sprintf(tmpstr,"%s %g",str,mybuf[m+1]);
    strcpy(str,tmpstr);
    m += 2;

      offset += sprintf(&sbuf[offset],"%s\n",str);
  }
  return offset;
}

/* ---------------------------------------------------------------------- */

void DumpLdd::write_binary(int n, double *mybuf)
{
  n *= size_one;
  fwrite(&n,sizeof(int),1,fp);
  fwrite(mybuf,sizeof(double),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpLdd::write_string(int n, double *mybuf)
{
  fwrite(mybuf,sizeof(char),n,fp);
}

/* ---------------------------------------------------------------------- */

void DumpLdd::write_lines_image(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
            static_cast<tagint> (mybuf[m]), static_cast<int> (mybuf[m+1]),
            mybuf[m+2],mybuf[m+3],mybuf[m+4], static_cast<int> (mybuf[m+5]),
            static_cast<int> (mybuf[m+6]), static_cast<int> (mybuf[m+7]));
    m += size_one;
  }
}

/* ---------------------------------------------------------------------- */

void DumpLdd::write_lines_noimage(int n, double *mybuf)
{
  int m = 0;
  for (int i = 0; i < n; i++) {
    fprintf(fp,format,
            static_cast<tagint> (mybuf[m]), static_cast<int> (mybuf[m+1]),
            mybuf[m+2],mybuf[m+3],mybuf[m+4]);
    m += size_one;
  }
}
