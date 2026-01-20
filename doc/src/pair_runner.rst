
`pair_style runner` command
===========================

Syntax
------

.. code-block:: lammps

   pair_style runner keyword value ...

* one or more keywords must be appended
* keyword = *dir* or *cflength* or *cfenergy* or *f_comm* or *q_comm* or *committee_size* or *total_charge* or *use_prev_q* or *check_extrap* or *max_extrap* or *show_ew* or *sum_ew_freq* or *reset_ew_freq*

.. code-block:: lammps

    dir path = path to the directory containing RuNNer potential files (default: ./)
    cflength value = length unit conversion factor (default: 1.0)
    cfenergy value = energy unit conversion factor (default: 1.0)
    f_comm yes/no = optionally write individual committee forces to a custom per-atom array (default: no)
    q_comm yes/no = optionally write individual committee charges to a custom per-atom array (default: no)
    committee_size value = number of models in the committee (default: 1)
    total_charge value = total system charge for electrostatics/QEq (default: 0.0)
    use_prev_q yes/no = use charges from previous step as starting guess for QEq (default: no)
    check_extrap yes/no = enable monitoring of feature extrapolation (default: no)
    max_extrap value = stop simulation if extrapolation count exceeds this value (default: 100)
    show_ew yes/no = print extrapolation warning messages to the log (default: no)
    sum_ew_freq value = frequency to print summary of extrapolations (default: 0)
    reset_ew_freq value = frequency to reset extrapolation counters (default: 0)

Examples
--------

.. code-block:: lammps

   pair_style runner dir ./potential_files/ committee_size 4 check_extrap yes
   pair_coeff * * 1 8 14

   # Using 4G-HDNNP with charge equilibration and committee force output
   fix 1 all property/atom ghost f_comm 12 # 4 committee members * 3 components
   pair_style runner dir ./nnp_model/ committee_size 4 f_comm yes total_charge 0.0
   pair_coeff * * 1 6 8

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

If `f_comm` or `q_comm` are set to `yes`, the individual forces or charges for every committee member can be stored. This requires the user to define a custom per-atom array using :doc:`fix property/atom <fix_property_atom>` with the name `f_comm` (size ``3 * committee_size``) or `q_comm` (size ``committee_size``).

**Extrapolation Monitoring:**

Since machine learning potentials are most reliable within their training data range, the *runner* style can monitor if the local atomic environments (features) extrapolate beyond the training set.
* `check_extrap`: Toggles monitoring.
* `show_ew`: Logs specific atoms and timesteps where extrapolation occurs.
* `max_extrap`: Provides a safety shutoff if the potential becomes unreliable.

**Unit Conversion:**

The RuNNer library is unit-agnostic. Use `cflength` and `cfenergy` to scale LAMMPS coordinates and energies to the units the potential was trained in. For example, if your LAMMPS simulation uses `metal` units (Angstroms, eV) and your RuNNer model was trained in Bohr and Hartrees:
* `cflength`: 1.889726 (Angstrom to Bohr)
* `cfenergy`: 0.036749 (eV to Hartree)

---

Mixing, shift, table, tail correction, restart, r_cut
-----------------------------------------------------

This pair style does not support mixing. It does not support `pair_modify` shift, table, or tail options. It does not write its own restart files; the potential is re-initialized from the input files in the specified directory.

Restrictions
------------

* This pair style is part of the ML-RUNNER package. It is only enabled if LAMMPS was built with that package.
* **Newton Pair:** This style requires `newton on`.
* **Single Instance:** Currently, only one instance of `pair_style runner` can be initialized per simulation.
* **Periodicity:** 3G and 4G HDNNPs require either full 3D periodicity or no periodicity. Partial periodicity (e.g., `boundary p p f`) is not supported for long-range electrostatics.
* **Charge Neutrality:** For periodic systems using 3G/4G models, the system must be charge neutral (`total_charge 0.0`).

Related commands
----------------

:doc:`pair_coeff <pair_coeff>`, :doc:`fix property/atom <fix_property_atom>`

References
----------

**(Behler)** J. Behler and M. Parrinello, Phys. Rev. Lett. 98, 146401 (2007).
**(Behler)** J. Behler, J. Chem. Phys. 134, 074106 (2011).
**(Ko)** T. W. Ko, J. A. Finkler, S. Goedecker, and J. Behler, Nat. Commun. 12, 398 (2021). (4G-HDNNPs)
