#include "fix_gemc.h"

#include "atom.h"
#include "atom_vec.h"
#include "comm.h"
#include "compute.h"
#include "domain.h"
#include "error.h"
#include "force.h"
#include "group.h"
#include "input.h"
#include "irregular.h"
#include "memory.h"
#include "modify.h"
#include "molecule.h"
#include "neighbor.h"
#include "pair.h"
#include "random_park.h"
#include "universe.h"
#include "update.h"
#include "variable.h"

// for molecule
#include "angle.h"
#include "bond.h"
#include "dihedral.h"
#include "improper.h"
#include "kspace.h"

#include <cstring>

using namespace LAMMPS_NS;

// this must be lower than MAXENERGYSIGNAL
// by a large amount, so that it is still
// less than total energy when negative
// energy contributions are added to MAXENERGYSIGNAL

static constexpr double MAXENERGYTEST = 1.0e50;

static constexpr double BUFFACTOR = 1.2;
static constexpr int BUFMIN = 1024;

static constexpr int NO_TAG = 0;

// translation boxes

/* ----------------------------------------------------------------------
  Shrink/expand boxes (always requires full energy)
------------------------------------------------------------------------- */
void FixGEMC::attempt_volume_change_full()
{
  nvolume_attempts++;

  // current volume

  double Lx = xhi-xlo;
  double Ly = yhi-ylo;
  double Lz = zhi-zlo;
  double i_vol = Lx*Ly*Lz;
  double i_mass = group->mass(0);


  // sample volume change from world 0 comm 0
  double dvolume;
  if (mycomm == 0) {
    // each needs to sample volume which doesn't cause other to have
    // ... too large of density (only need max_mass)
    // have world 1 send current vol and mass to world 0

    double j_vol, j_mass;

    // replace with below after updating mpi
    //MPI_Sendrecv(&i_mass, 1, MPI_DOUBLE, 0, NO_TAG,
    //             &j_mass, 1, MPI_DOUBLE, 1, NO_TAG, comm_replica, NULL);
    //MPI_Sendrecv(&i_vol, 1, MPI_DOUBLE, 0, NO_TAG,
    //             &j_vol, 1, MPI_DOUBLE, 1, NO_TAG, comm_replica, NULL);

    if (myworld == 1) {
      MPI_Send(&i_mass, 1, MPI_DOUBLE, 0, NO_TAG, comm_replica);
      MPI_Send(&i_vol, 1, MPI_DOUBLE, 0, NO_TAG, comm_replica);
    } else {
      MPI_Recv(&j_mass, 1, MPI_DOUBLE, 1, NO_TAG, comm_replica, NULL);
      MPI_Recv(&j_vol, 1, MPI_DOUBLE, 1, NO_TAG, comm_replica, NULL);
    }

    // just have world 0 sample all the volume changes
    // make sure volume change less than available volume

    if (myworld == 0) {
      double i_rho, j_rho;
      while (1) {
        dvolume = (random_proc->uniform()-0.5)*max_volume;
        i_rho = i_mass/(i_vol+dvolume);
        j_rho = j_mass/(j_vol-dvolume); // world 1 will swap sign
        if (MAX(i_rho,j_rho) < max_rho && i_rho > 0 && j_rho > 0) break;
      }
    }

    MPI_Bcast(&dvolume, 1, MPI_DOUBLE, 0, comm_replica);
    if (myworld == 1) dvolume *= -1.0;
  }

  // broadcast volume change to all worlds

  MPI_Bcast(&dvolume, 1, MPI_DOUBLE, 0, world);

  // attempt to change volume

  double fvolume = (i_vol+dvolume)/i_vol;
  if (fvolume < 0) error->universe_one(FLERR,"Negative volume found in fix gemc");
  double scale_length = pow(fvolume, 1.0/domain->dimension);

  // convert to lamda coords so they get scaled
  // TODO : If this is only natom_local, then there is huge spike in energy
  // .... but not with natom_total. So ghost atoms also need to be shifted?

  domain->x2lamda(natom_total);
  for (auto &ifix : rfix) ifix->deform(0);

  // shrink box toward lower corner
  // lower box coordinates always same
  // find center point of box

  xhi_tmp = xlo + Lx*scale_length;
  yhi_tmp = ylo + Ly*scale_length;
  zhi_tmp = zlo + Lz*scale_length;

  // set temporarily

  domain->boxhi[0] = xhi_tmp;
  domain->boxhi[1] = yhi_tmp;
  domain->boxhi[2] = zhi_tmp;

  // reset box and subbox dimensions

  domain->set_global_box();
  domain->set_local_box(); // reassigns sub domains

  // positions are scaled now

  domain->lamda2x(natom_total);
  for (auto &ifix : rfix) ifix->deform(1);

  // remap call

  domain->remap_all();

  // change in potential due to volume change

  double dU_volume = atom->natoms*force->boltz*box_temp*log(fvolume);

  // current system energy

  double energy_before = energy_stored;

  // (possible) future system energy

  double energy_after = energy_full();

  // get total energy from each partition

  double dU;
  if (mycomm == 0) {
    //printf("energy: %g -> %g\n", energy_before, energy_after);
    double idU = energy_after-energy_before-dU_volume;
    // sum change in full energy across each box
    MPI_Allreduce(&idU, &dU, 1, MPI_DOUBLE, MPI_SUM, comm_replica);
  }

  // bcast potential change to rest of world

  MPI_Bcast(&dU, 1, MPI_DOUBLE, 0, world);

  // evaluate probability

  double prob = MIN(exp(-beta*dU),1.0);
  double rf = random_universe->uniform();

  // volume change rejected -> revert atom positions

  if (prob < rf || energy_after >= MAXENERGYTEST) {

    domain->x2lamda(natom_total);
    for (auto &ifix : rfix) ifix->deform(0);

    domain->boxhi[0] = xhi;
    domain->boxhi[1] = yhi;
    domain->boxhi[2] = zhi;

    // reset box and subbox dimensions

    domain->set_global_box();
    domain->set_local_box(); // reassigns sub domains

    //irregular->migrate_atoms();

    domain->lamda2x(natom_total);
    for (auto &ifix : rfix) ifix->deform(1);

    // remap call (may lose atoms if no remap)

    domain->remap_all(); // maybe?

    // build neighbor
    neighbor->build(1);

  // acccept volume change

  } else {
    nvolume_successes++;

    // store new energy

    energy_stored = energy_after;

    // reacquire upper domain bounds

    xhi = domain->boxhi[0];
    yhi = domain->boxhi[1];
    zhi = domain->boxhi[2];

    // reacquire subdomain bounds

    if (triclinic_flag) {
      sublo = domain->sublo_lamda;
      subhi = domain->subhi_lamda;
    } else {
      sublo = domain->sublo;
      subhi = domain->subhi;
    }
  }

  if (energy_stored > MAXENERGYTEST) {
    printf("[%i] volume - %g -> %g\n",
      myworld, energy_before, energy_after);
    error->universe_one(FLERR,"bad energy");
  }
}

/* ----------------------------------------------------------------------
  Attempt atom exchange
------------------------------------------------------------------------- */

// TODO: do we need the force->kspace and force->pair->tail_flag?
void FixGEMC::attempt_atomic_exchange_full()
{
  nexchange_attempts++;

  // Choose sender and receiver

  int sender;
  if (mycomm == 0) {
    double drand = random_proc->uniform();
    double dmean;
    MPI_Allreduce(&drand, &dmean, 1, MPI_DOUBLE, MPI_SUM, comm_replica);
    dmean *= 0.5;
    if (drand > dmean) sender = 1;
    else sender = 0;
  }
  MPI_Bcast(&sender, 1, MPI_INT, 0, world);

  // Setup buffer for atom exchange

  int maxbuf_tmp = init_exchange() + BUFMIN;
  if (maxbuf_tmp > maxbuf) {
    maxbuf = maxbuf_tmp * BUFFACTOR;
    memory->grow(buf,maxbuf,"fix_gemc:buf_send");
  }

  // atom to delete/insert

  int iatom = -1;
  int tmp_mask;
  double q_tmp;

  // save old coordinates in case exchange rejected

  double old_coord[3];

  // pick atom to send

  if (sender) {

    // pick one atom randomly from all atoms in box
    // only one proc will actually delete atom

    iatom = pick_random_gas_atom();
    if (iatom >= 0) {
      // save old coordinates

      old_coord[0] = atom->x[iatom][0];
      old_coord[1] = atom->x[iatom][1];
      old_coord[2] = atom->x[iatom][2];

      // pack atom (only one atom sent per move)

      atom->avec->pack_exchange(iatom,&buf[0]);

      // temporarily set mask to exclusion for full energy later

      tmp_mask = atom->mask[iatom];

      // temporarily zero out charge for kspace later)

      if (q_flag) {
        q_tmp = atom->q[iatom];
        atom->q[iatom] = 0.0;
      }
      atom->mask[iatom] = exclusion_group_bit;
    }

    // send buffer from to comm 0 with
    // each exchange move will only have two procs send/recv per box
    // don't need mpi_barrier
    // exclude case where comm 0 already has the information

    if (iatom >= 0 && mycomm != 0) {
      MPI_Send(&buf[0], maxbuf, MPI_DOUBLE, 0, 0, world);
    } else if (iatom < 0 && mycomm == 0) {
      MPI_Recv(&buf[0], maxbuf, MPI_DOUBLE, MPI_ANY_SOURCE,
               0, world, MPI_STATUS_IGNORE);
    }
  }

  // send over atom thru comm 0's

  if (mycomm == 0) {
    // send buffer from sender to receiver
    // there's only two procs in comm_replica, so other is always 1-myrank

    if (sender) MPI_Send(&buf[0], maxbuf, MPI_DOUBLE,
                         1-myrank_replica, 0, comm_replica);
    else MPI_Recv(&buf[0], maxbuf, MPI_DOUBLE,
                  MPI_ANY_SOURCE, 0, comm_replica, MPI_STATUS_IGNORE);
  }

  // for now bcast buf to all procs

  if (!sender) MPI_Bcast(&buf[0], maxbuf, MPI_DOUBLE, 0, world);

  // pick random proc to place atom in

  int proc_flag = 0;
  if (!sender) {

    // sample random point in box

    double lamda[3], coord[3];
    if (mycomm == 0) {
      if (triclinic_flag) {
        lamda[0] = random_proc->uniform();
        lamda[1] = random_proc->uniform();
        lamda[2] = random_proc->uniform();

        // wasteful, but necessary

        if (lamda[0] == 1.0) lamda[0] = 0.0;
        if (lamda[1] == 1.0) lamda[1] = 0.0;
        if (lamda[2] == 1.0) lamda[2] = 0.0;

        domain->lamda2x(lamda,coord);
      } else {
        coord[0] = xlo + random_proc->uniform() * (xhi-xlo);
        coord[1] = ylo + random_proc->uniform() * (yhi-ylo);
        coord[2] = zlo + random_proc->uniform() * (zhi-zlo);
      }
    } // END mycomm

    // find proc that contains coordinate

    MPI_Bcast(&coord, 3, MPI_DOUBLE, 0, world);
    if (triclinic_flag) {
      if (lamda[0] >= sublo[0] && lamda[0] < subhi[0] &&
          lamda[1] >= sublo[1] && lamda[1] < subhi[1] &&
          lamda[2] >= sublo[2] && lamda[2] < subhi[2]) proc_flag = 1;
    } else {
      domain->remap(coord);
      if (!domain->inside(coord))
        error->universe_one(FLERR,"Fix gemc put atom outside box");
      if (coord[0] >= sublo[0] && coord[0] < subhi[0] &&
          coord[1] >= sublo[1] && coord[1] < subhi[1] &&
          coord[2] >= sublo[2] && coord[2] < subhi[2]) proc_flag = 1;
    } // END if triclinic

    // unpack atom here (only one atom should be received per move)
    // this will also create an atom and add to list (only to nlocal)

    if (proc_flag) {
      // confirmed that charge is stored in avec

      atom->avec->unpack_exchange(&buf[0]);

      int m = atom->nlocal - 1;

      // overwrite coordinates with new ones

      atom->x[m][0] = coord[0];
      atom->x[m][1] = coord[1];
      atom->x[m][2] = coord[2];
    }

    atom->natoms++;
    if (atom->tag_enable) {
      atom->tag_extend();
      if (atom->map_style != Atom::MAP_NONE) atom->map_init();
    }
  } // END if sender

  update_gas_atoms_list();

  // evalute probability for exchange

  if (triclinic_flag) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic_flag) domain->lamda2x(atom->nlocal+atom->nghost);

  if (force->kspace) force->kspace->qsum_qsq();
  if (force->pair->tail_flag) force->pair->reinit();
  double energy_before = energy_stored;
  double energy_after = energy_full();

  int success;
  if (mycomm == 0) {
    double all_dU;
    double idU = energy_after-energy_before;
    MPI_Allreduce(&idU,&all_dU,1,MPI_DOUBLE,MPI_SUM,comm_replica);

    double volume = (xhi-xlo)*(yhi-ylo)*(zhi-zlo);
    double NV;

    if (sender) NV = volume/(natom_total-1);
    else NV = natom_total/volume;
    double allNV;
    MPI_Allreduce(&NV,&allNV,1,MPI_DOUBLE,MPI_PROD,comm_replica);

    all_dU += (box_temp*force->boltz*log(allNV));
    double prob = MIN(exp(-beta*all_dU),1.0);

    if (prob > random_proc->uniform() && energy_after < MAXENERGYTEST)
      success = 1;
    else
      success = 0;

    MPI_Bcast(&success, 1, MPI_INT, 0, comm_replica);
  }
  MPI_Bcast(&success, 1, MPI_INT, 0, world);

  if (success && energy_after > MAXENERGYTEST) {
    printf("bad exchange\n");
    printf("%i -- %i\n", myworld, atom->natoms);
    printf("%g %g\n", energy_before, energy_after);
    error->universe_one(FLERR,"bad energy");
  }

  // handle deletion/insertions or revert

  if (sender) {
    // delete iatom

    if (success) {
      nexchange_successes++;
      if (iatom >= 0) {
        // overwrite iatom with last atom details

        atom->avec->copy(atom->nlocal-1,iatom,1);
        atom->nlocal--;
      }
      atom->natoms--;
      if (atom->map_style != Atom::MAP_NONE) atom->map_init();
      energy_stored = energy_after;

    // packing does not delete atom, just need to revert mask
    } else {
      // revert mask and charge if deletion rejected

      if (iatom >= 0) {
        atom->mask[iatom] = tmp_mask;
        if (q_flag) atom->q[iatom] = q_tmp;
      }
      if (force->kspace) force->kspace->qsum_qsq();
      if (force->pair->tail_flag) force->pair->reinit();
    }
  } else {
    // accept newly inserted iatom there

    if (success) {
      nexchange_successes++;
      energy_stored = energy_after;

    // remove newly inserted iatom (it was added to end)
    } else {
      atom->natoms--;
      if (proc_flag) atom->nlocal--;
      if (force->kspace) force->kspace->qsum_qsq();
      if (force->pair->tail_flag) force->pair->reinit();
    }
  }

  // update counts

  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

void FixGEMC::attempt_atomic_translation_full()
{
  ntranslation_attempts++;

  if (natom_total == 0) return;

  double energy_before = energy_stored;

  int iatom = pick_random_gas_atom();

  double **x = atom->x;
  double xtmp[3];

  xtmp[0] = xtmp[1] = xtmp[2] = 0.0;

  tagint tmptag = 0;

  double rx,ry,rz,rsq;
  rx = ry = rz = 0.0;
  if (iatom >= 0) {
    while(1) {
      rsq = 1.1;
      rx = ry = rz = 0.0;
      while (rsq > 1.0) {
        rx = 2*random_proc->uniform() - 1.0;
        ry = 2*random_proc->uniform() - 1.0;
        rz = 2*random_proc->uniform() - 1.0;
        rsq = rx*rx + ry*ry + rz*rz;
      }

      double coord[3];
      coord[0] = x[iatom][0] + displace*rx;
      coord[1] = x[iatom][1] + displace*ry;
      coord[2] = x[iatom][2] + displace*rz;
      if (coord[0]>xlo && coord[0]<xhi &&
          coord[1]>ylo && coord[1]<yhi &&
          coord[2]>zlo && coord[2]<zhi) break;
    }

    x[iatom][0] += displace*rx;
    x[iatom][1] += displace*ry;
    x[iatom][2] += displace*rz;

    tmptag = atom->tag[iatom];
  }

  if (triclinic_flag) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic_flag) domain->lamda2x(atom->nlocal+atom->nghost);
  double energy_after = energy_full();

  double prob = MIN(1.0,exp(beta*(energy_before - energy_after)));
  int success = 0;
  if (energy_after < MAXENERGYTEST && random_world->uniform() < prob)
    success = 1;

  if (success) {
    energy_stored = energy_after;
    ntranslation_successes++;
  } else {

    tagint tmptag_all;
    MPI_Allreduce(&tmptag,&tmptag_all,1,MPI_LMP_TAGINT,MPI_MAX,world);

    double rx_all, ry_all, rz_all;
    MPI_Allreduce(&rx,&rx_all,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&ry,&ry_all,1,MPI_DOUBLE,MPI_SUM,world);
    MPI_Allreduce(&rz,&rz_all,1,MPI_DOUBLE,MPI_SUM,world);

    if (iatom >= 0) {
      x[iatom][0] -= displace*rx_all;
      x[iatom][1] -= displace*ry_all;
      x[iatom][2] -= displace*rz_all;
    }

    energy_stored = energy_before;
  }

  if (triclinic_flag) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic_flag) domain->lamda2x(atom->nlocal+atom->nghost);

  update_gas_atoms_list();
}

/* ----------------------------------------------------------------------
------------------------------------------------------------------------- */

int FixGEMC::pick_random_gas_atom()
{
  int i = -1;
  int iwhichglobal = static_cast<int> (natom_total*random_world->uniform());
  if ((iwhichglobal >= natom_lower) &&
      (iwhichglobal < natom_lower + natom_local)) {
    i = iwhichglobal - natom_lower;
  }

  return i;
}
