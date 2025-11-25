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
   Contributing author: Andrew Hong, Aidan Thompson (SNL)
------------------------------------------------------------------------- */

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
using namespace FixConst;

// large energy value used to signal overlap

static constexpr double MAXENERGYSIGNAL = 1.0e100;

// this must be lower than MAXENERGYSIGNAL
// by a large amount, so that it is still
// less than total energy when negative
// energy contributions are added to MAXENERGYSIGNAL

static constexpr double MAXENERGYTEST = 1.0e50;

static constexpr double COMMBUFFACTOR = 1.2;
static constexpr int COMMBUFMIN = 1024;

/* ---------------------------------------------------------------------- */

FixGEMC::FixGEMC(LAMMPS *lmp, int narg, char **arg) : Fix(lmp, narg, arg)
{
  if (narg != 11) utils::missing_cmd_args(FLERR, "fix gemc", error);

  // must have only two boxes

  if (universe->nworlds != 2) error->universe_all(FLERR, "Must use exactly two partitions with fix gemc");

  molecule_flag = atom->molecule_flag;
  if (molecule_flag) error->universe_all(FLERR, "Using fix gemc with molecules not (yet) supported");

  // various fix flags (partial copy from gcmc)

  time_integrate = 0; // do not time integrate (use only MC moves)
  global_freq = 1;
  time_depend = 1;
  restart_global = 1;

  // box size changes with volume MC moves

  box_change |= BOX_CHANGE_X;
  box_change |= BOX_CHANGE_Y;
  box_change |= BOX_CHANGE_Z;

  // set up reneighboring

  force_reneighbor = 1;
  next_reneighbor = update->ntimestep + 1;
  nrotate = 0;

  // required user args

  nevery = utils::inumeric(FLERR,     arg[3], false, lmp);
  ntranslate = utils::inumeric(FLERR, arg[4], false, lmp);
  nexchange = utils::inumeric(FLERR,  arg[5], false, lmp);
  nvolume = utils::inumeric(FLERR,    arg[6], false, lmp);
  box_temp = utils::numeric(FLERR,    arg[7], false, lmp);
  displace = utils::numeric(FLERR,    arg[8], false, lmp);
  max_volume = utils::numeric(FLERR,  arg[9], false, lmp);
  seed = utils::inumeric(FLERR,       arg[10], false, lmp);

  // set up comm_replica = communicator between proc 0s across boxes

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);
  myworld = universe->iworld;

  MPI_Comm_split(universe->uworld, me, 0, &comm_replica);

  // use same RNG for each replica for volume MC moves
  // unique to proc

  random_proc = new RanPark(lmp,seed+3.0*(universe->iworld+1)+7.0*(me+1));

  // sync between procs

  random_world = new RanPark(lmp,seed+13.0*(universe->iworld+1));

  // sync between universes

  random_universe = new RanPark(lmp,seed);

  // detect if any rigid fixes exist so rigid bodies move when box is remapped

  rfix.clear();
  for (auto &ifix : modify->get_fix_list())
    if (ifix->rigid_flag) rfix.push_back(ifix);

  gemc_nmax = 0;
  local_gas_list = nullptr;
  commbuf = nullptr;
  maxcommbuf = 0;

}

/* ---------------------------------------------------------------------- */

FixGEMC::~FixGEMC()
{
  memory->destroy(local_gas_list);
  MPI_Comm_free(&comm_replica);
  memory->destroy(commbuf);
}

/* ---------------------------------------------------------------------- */

int FixGEMC::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixGEMC::init()
{
  progress = 0;

  // determine probability of each step during pre_exchange

  // set probabilities for MC moves

  if (!molecule_flag) nrotate = 0;

  // total moves is a double to avoid type casting later

  nmoves = nvolume + nexchange + ntranslate + nrotate;

  double d_nmoves = static_cast<double>(nmoves);
  double p_exchange  = nexchange/d_nmoves;
  double p_volume    = nvolume/d_nmoves;
  double p_translate = ntranslate/d_nmoves;
  double p_rotate    = nrotate/d_nmoves;

  // normalize probabilities

  double p_total = p_exchange + p_volume + p_translate + p_rotate;

  // compute cumulative probabilities

  pc_exchange = p_exchange/p_total;
  pc_volume = (p_volume+p_exchange)/p_total;
  if (molecule_flag) {
    pc_translate = (p_volume+p_exchange+p_translate)/p_total;
    pc_rotate = 1.0;
  } else {
    pc_translate = 1.0;
    pc_rotate = 0.0;
  }

  // for full energy

  c_pe = modify->compute[modify->find_compute("thermo_pe")];

  // check if atoms charged

  q_flag = atom->q_flag;

  // pre compute scaled temperature

  beta = 1.0/(force->boltz*box_temp);

  // get domain dim

  triclinic_flag = domain->triclinic;

  // get domain dims

  xlo = domain->boxlo[0];
  xhi = domain->boxhi[0];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];

  // get subdomain

  if (triclinic_flag) {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  } else {
    sublo = domain->sublo;
    subhi = domain->subhi;
  }

  // create unique group name for atoms to be excluded for particle exchange
  // keeps temporarily deleted particles from being added in potential energy calc

  // id from fix

  auto group_id = std::string("FixGEMC:gemc_exclusion_group:") + id;
  group->assign(group_id + " subtract all all");
  exclusion_group = group->find(group_id);
  if (exclusion_group == -1)
    error->universe_all(FLERR,"Could not find fix gemc exclusion group ID");
  exclusion_group_bit = group->bitmask[exclusion_group];

  // neighbor list exclusion setup
  // turn off interactions between group all and the exclusion group

  neighbor->modify_params(fmt::format("exclude group {} all",group_id));

  // allocate communication buffer for exchange moves

  grow_commbuf();

  groupbitall = 1 | groupbit;

  ntranslation_attempts = ntranslation_successes = 0.0;
  nvolume_attempts = nvolume_successes = 0.0;
  nexchange_attempts = nexchange_successes = 0.0;
}

/* ----------------------------------------------------------------------
   attempt Monte Carlo translations, rotations, insertions, and deletions
   done before exchange, borders, reneighbor
   so that ghost atoms and neighbor lists will be correct

   gcmc + extra volume step and modified insertion/deletion
------------------------------------------------------------------------- */

void FixGEMC::pre_exchange()
{
  // just return if should not be called on this timestep

  if (next_reneighbor != update->ntimestep) return;

  // get domain dims

  xlo = domain->boxlo[0];
  xhi = domain->boxhi[0];
  ylo = domain->boxlo[1];
  yhi = domain->boxhi[1];
  zlo = domain->boxlo[2];
  zhi = domain->boxhi[2];

  // get subdomain

  if (triclinic_flag) {
    sublo = domain->sublo_lamda;
    subhi = domain->subhi_lamda;
  } else {
    sublo = domain->sublo;
    subhi = domain->subhi;
  }

  // three steps in GEMC:
  // 1) translate particles within each box
  // 2) exchange particles between boxes
  // 3) change box volume + scale particle positions

  // update next time to call

  next_reneighbor = update->ntimestep + nevery;

  double imove;
  update_gas_atoms_list();

  energy_stored = energy_full();
  int prev_step = 0;

   for (int i = 0; i < nmoves; i++) {
    imove = random_universe->uniform();

    if (imove < pc_exchange) attempt_atomic_exchange_full();
    else if (imove < pc_volume) attempt_volume_change_full();
    else attempt_atomic_translation_full();

  }

  // print progress info to universe screen/logfile

  if (universe->me == 0) {
    double delta = update->ntimestep - update->beginstep;
    if ((delta != 0.0) && (update->beginstep != update->endstep))
      delta /= update->endstep - update->beginstep;
    int status = static_cast<int>(delta * 100.0);
    if (status > progress) {
      progress = status;
      auto msg = fmt::format(
          " GEMC run progress: {:>3d}% \n  Trans: {:g}/{:g}\n"
          "  Vol: {:g}/{:g}\n  Ex: {:g}/{:g}\n",
          progress,
          ntranslation_successes, ntranslation_attempts,
          nvolume_successes, nvolume_attempts,
          nexchange_successes, nexchange_attempts);
      if (universe->uscreen) utils::print(universe->uscreen, msg);
      if (universe->ulogfile) utils::print(universe->ulogfile, msg);
      ntranslation_attempts = ntranslation_successes = 0.0;
      nvolume_attempts = nvolume_successes = 0.0;
      nexchange_attempts = nexchange_successes = 0.0;
    }
  }
}

/* ----------------------------------------------------------------------
   update per-proc atom count (same for molecules)
   assume all atoms are candidates for MC moves
------------------------------------------------------------------------- */

void FixGEMC::update_gas_atoms_list()
{
  int nlocal = atom->nlocal;
  int *mask = atom->mask;

  natom_local = 0;
  for (int i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {
      natom_local++;
    }
  }

  // ngas is total atoms in whole system

  MPI_Allreduce(&natom_local,&natom_total,1,MPI_INT,MPI_SUM,world);
  MPI_Scan(&natom_local,&natom_lower,1,MPI_INT,MPI_SUM,world);
  natom_lower -= natom_local;
}

/* ----------------------------------------------------------------------
   compute system potential energy
------------------------------------------------------------------------- */

double FixGEMC::energy_full()
{
  if (triclinic_flag) domain->x2lamda(atom->nlocal);
  domain->pbc();
  comm->exchange();
  atom->nghost = 0;
  comm->borders();
  if (triclinic_flag) domain->lamda2x(atom->nlocal+atom->nghost);
  if (modify->n_pre_neighbor) modify->pre_neighbor();
  neighbor->build(1);
  int eflag = 1;
  int vflag = 0;

  // clear forces so they don't accumulate over multiple
  // calls within fix gcmc timestep, e.g. for fix shake

  size_t nbytes = sizeof(double) * (atom->nlocal + atom->nghost);
  if (nbytes) memset(&atom->f[0][0],0,3*nbytes);

  if (modify->n_pre_force) modify->pre_force(vflag);

  if (force->pair) force->pair->compute(eflag,vflag);

  if (atom->molecular != Atom::ATOMIC) {
    if (force->bond) force->bond->compute(eflag,vflag);
    if (force->angle) force->angle->compute(eflag,vflag);
    if (force->dihedral) force->dihedral->compute(eflag,vflag);
    if (force->improper) force->improper->compute(eflag,vflag);
  }

  if (force->kspace) force->kspace->compute(eflag,vflag);

  if (modify->n_post_force_any) modify->post_force(vflag);

  // NOTE: all fixes with energy_global_flag set and which
  //   operate at pre_force() or post_force()
  //   and which user has enabled via fix_modify energy yes,
  //   will contribute to total MC energy via pe->compute_scalar()

  update->eflag_global = update->ntimestep;
  double total_energy = c_pe->compute_scalar();

  return total_energy;
}

/* ----------------------------------------------------------------------
   set bufextra based on AtomVec and fixes
   does not include base data to exchange
   similar to Comm::init_exchange()
   *** this is formally an error, since base data is also packed ***
------------------------------------------------------------------------- */

void FixGEMC::grow_commbuf()
{
  int ntmp = 0;
  for (auto &ifix : modify->get_fix_list())
    ntmp = MAX(ntmp, ifix->maxexchange);
  ntmp += atom->avec->maxexchange;
  ntmp += COMMBUFMIN;
  ntmp *= COMMBUFFACTOR;

  if (maxcommbuf < ntmp) {
    maxcommbuf = ntmp;
    memory->grow(commbuf,maxcommbuf,"fix_gemc:commbuf");
  }
}

/* ----------------------------------------------------------------------
   pack entire state of Fix into one write
------------------------------------------------------------------------- */

void FixGEMC::write_restart(FILE *fp)
{
  int n = 0;
  double list[13];
  list[n++] = random_proc->state();
  list[n++] = random_world->state();
  list[n++] = random_universe->state();
  list[n++] = ubuf(next_reneighbor).d;
  list[n++] = ntranslation_attempts;
  list[n++] = ntranslation_successes;
  list[n++] = nrotation_attempts;
  list[n++] = nrotation_successes;
  list[n++] = nexchange_attempts;
  list[n++] = nexchange_successes;
  list[n++] = nvolume_attempts;
  list[n++] = nvolume_successes;
  list[n++] = ubuf(update->ntimestep).d;

  if (comm->me == 0) {
    int size = n * sizeof(double);
    fwrite(&size,sizeof(int),1,fp);
    fwrite(list,sizeof(double),n,fp);
  }
}

/* ----------------------------------------------------------------------
   use state info from restart file to restart the Fix
------------------------------------------------------------------------- */

void FixGEMC::restart(char *buf)
{
  int n = 0;
  auto *list = (double *) buf;

  seed = static_cast<int> (list[n++]);
  random_proc->reset(seed);

  seed = static_cast<int> (list[n++]);
  random_world->reset(seed);

  seed = static_cast<int> (list[n++]);
  random_universe->reset(seed);

  next_reneighbor = (bigint) ubuf(list[n++]).i;

  ntranslation_attempts  = list[n++];
  ntranslation_successes = list[n++];
  nrotation_attempts     = list[n++];
  nrotation_successes    = list[n++];
  nexchange_attempts     = list[n++];
  nexchange_successes    = list[n++];
  nvolume_attempts    = list[n++];
  nvolume_successes   = list[n++];

  bigint ntimestep_restart = (bigint) ubuf(list[n++]).i;
  if (ntimestep_restart != update->ntimestep)
    error->all(FLERR,"Must not reset timestep when restarting fix gemc");
}
