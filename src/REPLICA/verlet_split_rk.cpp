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

/* ----------------------------------------------------------------------
   Contributing authors: Yuxing Peng and Chris Knight (U Chicago)
------------------------------------------------------------------------- */

#include "verlet_split_rk.h"
#include "universe.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "atom.h"
#include "atom_vec.h"
#include "force.h"
#include "pair.h"
#include "bond.h"
#include "angle.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"
#include "output.h"
#include "update.h"
#include "fix.h"
#include "modify.h"
#include "timer.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

VerletSplitRK::VerletSplitRK(LAMMPS *lmp, int narg, char **arg) :
  Verlet(lmp, narg, arg)
{
  ;
}

/* ---------------------------------------------------------------------- */

VerletSplitRK::~VerletSplitRK()
{
  ;
}

/* ----------------------------------------------------------------------
   initialization before run
------------------------------------------------------------------------- */

void VerletSplitRK::init()
{
  if (comm->style != Comm::BRICK)
    error->universe_all(FLERR,"Verlet/split/rk can only currently be used with comm_style brick");
  if (!force->kspace && comm->me == 0)
    error->warning(FLERR,"A KSpace style must be defined with verlet/split/rk");

  // error for as-yet unsupported verlet/split KSpace options

  int errflag = 0;
  if (!force->kspace->rk_flag) errflag = 1; //Currently requires the use of pppm/rk
  if (force->kspace->tip4pflag) errflag = 1;
  if (force->kspace->dipoleflag) errflag = 1;
  if (force->kspace->spinflag) errflag = 1;

  if (errflag)
    error->all(FLERR,"Verlet/split/rk cannot (yet) be used with kpace style {}", force->kspace_style);


  // invoke parent Verlet init

  Verlet::init();

}

/* ----------------------------------------------------------------------
   setup before run
   kproc partition only sets up KSpace calculation
------------------------------------------------------------------------- */

void VerletSplitRK::setup(int flag)
#if 0
{
  if (comm->me == 0 && screen)
    fprintf(screen,"Setting up Verlet/split/rk run ...\n");

  rproc = (universe->iworld == 0) ? 1 : 0;
  if (!rproc) force->kspace->setup();
  else Verlet::setup(flag);
}
#else
{
  if (comm->me == 0 && screen) {
    fputs("Setting up Verlet/split/rk run ...\n",screen);
    if (flag) {
      utils::print(screen,"  Unit style    : {}\n"
                        "  Current step  : {}\n"
                        "  Time step     : {}\n",
                 update->unit_style,update->ntimestep,update->dt);
      timer->print_timeout(screen);
    }
  }

  if (lmp->kokkos)
    error->all(FLERR,"KOKKOS package requires run_style verlet/kk");

  rproc = (universe->iworld == 0) ? 1 : 0;

  if(rproc){
    update->setupflag = 1;
  
    // setup domain, communication and neighboring
    // acquire ghosts
    // build neighbor lists
  
    atom->setup();
    modify->setup_pre_exchange();
    if (triclinic) domain->x2lamda(atom->nlocal);
    domain->pbc();
    domain->reset_box();
    comm->setup();
    if (neighbor->style) neighbor->setup_bins();
    comm->exchange();
    if (atom->sortfreq > 0) atom->sort();
    comm->borders();
    if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
    domain->image_check();
    domain->box_too_small_check();
    modify->setup_pre_neighbor();
    neighbor->build(1);
    modify->setup_post_neighbor();
    neighbor->ncalls = 0;
  
    // compute all forces
  
    force->setup();
    ev_set(update->ntimestep);
    force_clear();
    modify->setup_pre_force(vflag);
  }

  if (force->kspace){ 
    force->kspace->setup();
    if (kspace_compute_flag) 
    	force->kspace->r2k_comm(eflag,vflag);
  }
  if(rproc){
    if (pair_compute_flag) force->pair->compute(eflag,vflag);
    else if (force->pair) force->pair->compute_dummy(eflag,vflag);

    if (atom->molecular != Atom::ATOMIC) {
      if (force->bond) force->bond->compute(eflag,vflag);
      if (force->angle) force->angle->compute(eflag,vflag);
      if (force->dihedral) force->dihedral->compute(eflag,vflag);
      if (force->improper) force->improper->compute(eflag,vflag);
    }
  }
  if (force->kspace) {
    if (kspace_compute_flag){
	if(!rproc) force->kspace->compute_grid_potentials(eflag,vflag);

    	force->kspace->k2r_comm(eflag,vflag);
    }
    else force->kspace->compute_dummy(eflag,vflag);
  }


  if (rproc){
    modify->setup_pre_reverse(eflag,vflag);
    if (force->newton) comm->reverse_comm();

    modify->setup(vflag);
    output->setup(flag);
    update->setupflag = 0;
  }
}
#endif

/* ----------------------------------------------------------------------
   setup without output
   flag = 0 = just force calculation
   flag = 1 = reneighbor and force calculation
   kproc partition only sets up KSpace calculation
------------------------------------------------------------------------- */

void VerletSplitRK::setup_minimal(int flag)
#if 0
{
  rproc = (universe->iworld == 0) ? 1 : 0;
  if (!rproc) force->kspace->setup();
  else Verlet::setup_minimal(flag);
}
#else
{
  rproc = (universe->iworld == 0) ? 1 : 0;
  if(rproc){
    update->setupflag = 1;
  
    // setup domain, communication and neighboring
    // acquire ghosts
    // build neighbor lists
  
    if (flag) {
      modify->setup_pre_exchange();
      if (triclinic) domain->x2lamda(atom->nlocal);
      domain->pbc();
      domain->reset_box();
      comm->setup();
      if (neighbor->style) neighbor->setup_bins();
      comm->exchange();
      comm->borders();
      if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
      domain->image_check();
      domain->box_too_small_check();
      modify->setup_pre_neighbor();
      neighbor->build(1);
      modify->setup_post_neighbor();
      neighbor->ncalls = 0;
    }
  
    // compute all forces
  
    ev_set(update->ntimestep);
    force_clear();
    modify->setup_pre_force(vflag);
  }//rproc

  if (force->kspace){ 
    force->kspace->setup();
    if (kspace_compute_flag) 
    	force->kspace->r2k_comm(eflag,vflag);
  }

  if(rproc){
    if (pair_compute_flag) force->pair->compute(eflag,vflag);
    else if (force->pair) force->pair->compute_dummy(eflag,vflag);
  
    if (atom->molecular != Atom::ATOMIC) {
      if (force->bond) force->bond->compute(eflag,vflag);
      if (force->angle) force->angle->compute(eflag,vflag);
      if (force->dihedral) force->dihedral->compute(eflag,vflag);
      if (force->improper) force->improper->compute(eflag,vflag);
    }
  }//rproc

  if (force->kspace) {
    if (kspace_compute_flag){
	if(!rproc) force->kspace->compute_grid_potentials(eflag,vflag);

    	force->kspace->k2r_comm(eflag,vflag);
    }
    else force->kspace->compute_dummy(eflag,vflag);
  }

  if(rproc){
    modify->setup_pre_reverse(eflag,vflag);
    if (force->newton) comm->reverse_comm();
  
    modify->setup(vflag);
    update->setupflag = 0;
  }//rproc
}
#endif

/* ----------------------------------------------------------------------
   run for N steps
   rproc partition does everything but Kspace
   kproc partition does just Kspace
   communicate back and forth every step:
     atom coords from rproc -> kproc
     kspace forces from kproc -> rproc
     also box bounds from rproc -> kproc if necessary
------------------------------------------------------------------------- */

void VerletSplitRK::run(int n)
{
  bigint ntimestep;
  int nflag,sortflag;

  // sync both partitions before start timer

  MPI_Barrier(universe->uworld);
  timer->init();
  timer->barrier_start();

  // setup initial Rspace <-> Kspace comm params
  //rk_setup();

  // check if OpenMP support fix defined

  Fix *fix_omp;
  int ifix = modify->find_fix("package_omp");
  if (ifix < 0) fix_omp = nullptr;
  else fix_omp = modify->fix[ifix];

  // flags for timestepping iterations

  int n_post_integrate = modify->n_post_integrate;
  int n_pre_exchange = modify->n_pre_exchange;
  int n_pre_neighbor = modify->n_pre_neighbor;
  int n_pre_force = modify->n_pre_force;
  int n_pre_reverse = modify->n_pre_reverse;
  int n_post_force = modify->n_post_force_any;
  int n_end_of_step = modify->n_end_of_step;

  if (atom->sortfreq > 0) sortflag = 1;
  else sortflag = 0;

  for (int i = 0; i < n; i++) {

    ntimestep = ++update->ntimestep;
    ev_set(ntimestep);

    // initial time integration

    if (rproc) {
      modify->initial_integrate(vflag);
      if (n_post_integrate) modify->post_integrate();
#if 0
    }

    // regular communication vs neighbor list rebuild

    if (rproc) nflag = neighbor->decide();
    MPI_Bcast(&nflag,1,MPI_INT,1,block);

    if (rproc) {
#else
      nflag = neighbor->decide();
#endif
      if (nflag == 0) {
        timer->stamp();
        comm->forward_comm();
        timer->stamp(Timer::COMM);
      } else {
        if (n_pre_exchange) modify->pre_exchange();
        if (triclinic) domain->x2lamda(atom->nlocal);
        domain->pbc();
        if (domain->box_change) {
          domain->reset_box();
          comm->setup();
          if (neighbor->style) neighbor->setup_bins();
        }
        timer->stamp();
        comm->exchange();
        if (sortflag && ntimestep >= atom->nextsort) atom->sort();
        comm->borders();
        if (triclinic) domain->lamda2x(atom->nlocal+atom->nghost);
        timer->stamp(Timer::COMM);
        if (n_pre_neighbor) modify->pre_neighbor();
        neighbor->build(1);
        timer->stamp(Timer::NEIGH);
      }
    }

    // if reneighboring occurred, re-setup Rspace <-> Kspace comm params
    // comm Rspace atom coords to Kspace procs

#if 0
    if (nflag) rk_setup();
#endif
    force->kspace->r2k_comm(eflag,vflag);

#if 0
    force_clear();
#endif
    // force computations
    if (rproc) {
#if 1
      force_clear();
#endif
      if (n_pre_force) modify->pre_force(vflag);


      timer->stamp();
      if (force->pair) {
        force->pair->compute(eflag,vflag);
        timer->stamp(Timer::PAIR);
      }

      if (atom->molecular != Atom::ATOMIC) {
        if (force->bond) force->bond->compute(eflag,vflag);
        if (force->angle) force->angle->compute(eflag,vflag);
        if (force->dihedral) force->dihedral->compute(eflag,vflag);
        if (force->improper) force->improper->compute(eflag,vflag);
        timer->stamp(Timer::BOND);
      }

#if 0
      if (n_pre_reverse) {
        modify->pre_reverse(eflag,vflag);
        timer->stamp(Timer::MODIFY);
      }
      if (force->newton) {
        comm->reverse_comm();
        timer->stamp(Timer::COMM);
      }
#endif

    } else {
#if 0
      // run FixOMP as sole pre_force fix, if defined
      if (fix_omp) fix_omp->pre_force(vflag);

      if (force->kspace) {
        timer->stamp();
        force->kspace->compute(eflag,vflag);
        timer->stamp(Timer::KSPACE);
      }

      if (n_pre_reverse) {
        modify->pre_reverse(eflag,vflag);
        timer->stamp(Timer::MODIFY);
      }

      // TIP4P PPPM puts forces on ghost atoms, so must reverse_comm()

      if (tip4pflag && force->newton) {
        comm->reverse_comm();
        timer->stamp(Timer::COMM);
      }
#else
      timer->stamp();
      //Compute_grid_potentials(eflag,vflag);
      force->kspace->compute_grid_potentials(eflag,vflag);
      timer->stamp(Timer::KSPACE);
#endif
    }

    // comm and sum Kspace forces back to Rspace procs

    force->kspace->k2r_comm(eflag,vflag);
    //TODO
    if(rproc){
#if 1
      if (n_pre_reverse) {
        modify->pre_reverse(eflag,vflag);
        timer->stamp(Timer::MODIFY);
      }
      if (force->newton) {
        comm->reverse_comm();
        timer->stamp(Timer::COMM);
      }
#endif
    // force modifications, final time integration, diagnostics
    // all output
      timer->stamp();
      if (n_post_force) modify->post_force(vflag);
      modify->final_integrate();
      if (n_end_of_step) modify->end_of_step();
      timer->stamp(Timer::MODIFY);

      if (ntimestep == output->next) {
        timer->stamp();
        output->write(ntimestep);
        timer->stamp(Timer::OUTPUT);
      }
    } //rproc
  } //for i
}


