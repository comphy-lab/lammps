.. index:: pair_style runner

pair_style runner command
===========================

Syntax
------

.. code-block:: LAMMPS

   pair_style runner keyword value ...

* zero or more keyword/value pairs may be appended
* keyword = *dir* or *cflength* or *cfenergy* or *f_comm* or *q_comm* or *committee_size* or *total_charge* or *use_prev_q* or *check_extrap* or *max_extrap* or *show_ew* or *sum_ew_freq* or *reset_ew_freq*
* value depends on the preceding keyword:

.. parsed-literal::

    *dir* value = Path to RuNNer configuration files (default: *./*)

    *cflength* value = length
          Length unit conversion factor (default: *1.0*)
    *cfenergy* value = energy
          Unit conversion factor (default: *1.0*)

    *committee_size* value = number
          number = Number of committee members (default: *1*)
    *f_comm* value = *yes* or *no* 
          Write individual committee forces to a custom per-atom array (default: *no*)
    *q_comm* value = *yes* or *no* 
          Write individual committee charges to a custom per-atom array (default: *no*)

    *total_charge* value = charge
          charge = Total system charge for electrostatics/QEq (default: *0.0*)
    *use_prev_q* value = *yes* or *no* 
          Use predicted charges from previous time step as initial guess for QEq (default: *no*)

    *check_extrap* value = *yes* or *no*
          Enable monitoring of feature extrapolation (default: *no*)
    *max_extrap* value = threshold
          threshold = Stop simulation if EW count exceeds this value (default: *100*)
    *show_ew* value = *yes* or *no*
          Write EWs to the log (default: *no*)
    *sum_ew_freq* value = summary
          summary = Write EW summary every this many time steps (default: *0*)
    *reset_ew_freq* value = frequency
          frequency = Reset extrapolation counters every this many steps (default: *0*)

Examples
""""""""

.. code-block:: LAMMPS

   pair_style runner dir "./model_files"
   pair_coeff * * 1 6 8

**4G potential:**
.. code-block:: lammps
   # Using 4G-HDNNP with charge equilibration and committee force output
   fix 1 all property/atom ghost f_comm 12 # 4 committee members * 3 components
   pair_style runner dir ./nnp_model/ committee_size 4 f_comm yes total_charge 0.0
   pair_coeff * * 1 6 8


**Active Learning Setup (8-member Committee with output):**

.. code-block:: lammps

   # 1. Define storage for individual committee forces (8 members * 3 components = 24)
   fix 1 all property/atom f_comm 24 ghost yes

   # 2. Setup potential with units conversion (e.g., metal to atomic)
   pair_style runner dir "pot" cflength 1.8897 cfenergy 0.0367 committee_size 8 f_comm yes
   pair_coeff * * 1 8 14 26

   # 3. Define compute to access individual member energies
   compute e_comm all pair runner

   # 4. Output results
   thermo_style custom step temp epair c_e_comm[1] c_e_comm[2] c_e_comm[3]
   dump 1 all custom 1000 run.lammpstrj id element x y z fx fy fz f_comm[1] f_comm[2] f_comm[3]

Description
-----------

The *runner* pair style provides an interface to the **RuNNer** (Ruhr University Neural Network Energy Representation) library. It implements High-Dimensional Neural Network Potentials (HDNNPs) that represent the total energy of a system as a sum of environment-dependent atomic contributions.

RuNNer supports several "generations" of HDNNPs:

* **Second-generation (2G):** Short-range many-body potentials where the total energy is the sum of atomic energies predicted from local chemical environments.
* **Third-generation (3G):** Extends 2G by adding long-range electrostatic interactions based on environment-dependent atomic charges.
* **Fourth-generation (4G):** Includes global charge equilibration (QEq) based on environment-dependent electronegativities (and optionally hardness). These charges are fed back into the energy model, providing a global descriptor for the atomic energies.

Additionally, all generations can be augmented with:
* **Hirshfeld-based dispersion:** Long-range vdW interactions.
* **Repulsive potentials:** Short-range pairwise repulsive terms.

---

**Committee Approach:**

The pair style supports running multiple neural network models simultaneously. The forces, energies, and virials propagated in the simulation are the average of all committee members. This is useful for active learning and uncertainty estimation.

**Energy Output (via Compute):**
The individual potential energies of each committee member are accessed using the :doc:`compute pair <compute_pair>` command.

* Example: `compute ec all pair runner`
* Access: The energies are stored in a global vector of length `committee_size`. They are accessed as `c_ec[1]`, `c_ec[2]`, up to `c_ec[N]`.
* In `thermo_style`, these represent the potential energy of the system according to each specific model in the committee.

**Force Output (via Fix Property/Atom):**
To export forces of individual committee members (e.g., to compute force variance per atom), you must define a custom per-atom array **before** the `pair_style` command.

* **Name:** Must be `f_comm`.
* **Size:** Must be exactly $3 \times committee\_size$. (e.g., for a size 8 committee, the size is 24).
* **Ghost:** Must set `ghost yes` to ensure forces are correctly communicated.
* **Indexing:** The forces are stored sequentially. 
    * `f_comm[1], f_comm[2], f_comm[3]` = x, y, z forces for Member 1.
    * `f_comm[4], f_comm[5], f_comm[6]` = x, y, z forces for Member 2.
* **Example Dump:** `dump 1 all custom 100 out.dump id f_comm[1] f_comm[2] f_comm[3] f_comm[4] f_comm[5] f_comm[6]`

**Charge Output (4G only):**
If `q_comm` is set to `yes`, a per-atom array `fix q_comm all property/atom q_comm N ghost yes` (where N is `committee_size`) can be used to extract individual member charges.

---

**Extrapolation Monitoring:**

Since machine learning potentials are most reliable within their training data range, the *runner* style can monitor if the local atomic environments (features) extrapolate beyond the training set.
* `check_extrap`: Toggles monitoring.
* `show_ew`: Logs specific atoms and timesteps where extrapolation occurs.
* `max_extrap`: Provides a safety shutoff if the potential becomes unreliable.

---

Unit Conversion
---------------

The RuNNer library is unit-agnostic. Use `cflength` and `cfenergy` to scale LAMMPS coordinates and energies to the units the potential was trained in. For example, if your LAMMPS simulation uses `metal` units (Angstroms, eV) and your RuNNer model was trained in Bohr and Hartrees:
* `cflength`: 1.889726 (Angstrom to Bohr)
* `cfenergy`: 0.036749 (eV to Hartree)

---

Parallelization and Performance
-------------------------------

*   **MPI:** 2G models scale linearly. 3G/4G models require global structure collection on a single process for electrostatic/QEq calculations, which creates an MPI bottleneck.
*   **OpenMP:** RuNNer is heavily optimized for OpenMP. For 3G and 4G potentials, it is highly recommended to use a small number of MPI tasks and a large number of OpenMP threads per task.

Restrictions
------------

* This pair style is part of the ML-RUNNER package. It is only enabled if LAMMPS was built with that package.
* **Newton Pair:** This style requires `newton on`.
* **Single Instance:** Currently, only one instance of `pair_style runner` can be initialized per simulation.
* **Periodicity:** 3G and 4G HDNNPs require either full 3D periodicity or no periodicity. Partial periodicity (e.g., `boundary p p f`) is not supported for long-range electrostatics.
* **Charge Neutrality:** For periodic systems using 3G/4G models, the system must be charge neutral (`total_charge 0.0`).
* **mixing: ** This pair style does not support mixing.
* **pair_modify:** This style does not support `pair_modify` shift, table, or tail options.
* **Restart:** This pair style does not require its own restart files; the potential is re-initialized from the input files in the specified directory.

Related commands
----------------

:doc:`pair_coeff <pair_coeff>`, :doc:`fix property/atom <fix_property_atom>`

References
----------

**(Behler)** J. Behler and M. Parrinello, Phys. Rev. Lett. 98, 146401 (2007).
**(Behler)** J. Behler, J. Chem. Phys. 134, 074106 (2011).
**(Ko)** T. W. Ko, J. A. Finkler, S. Goedecker, and J. Behler, Nat. Commun. 12, 398 (2021). (4G-HDNNPs)
