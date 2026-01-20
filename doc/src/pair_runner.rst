.. index:: pair_style runner

pair_style runner command
=========================

Syntax
------

.. code-block:: LAMMPS

   pair_style runner keyword value ...

* zero or more keyword/value pairs may be appended
* keyword = *dir* or *cflength* or *cfenergy* or *f_comm* or *q_comm* or *committee_size* or *total_charge* or *use_prev_q* or *check_extrap* or *max_extrap* or *show_ew* or *sum_ew_freq* or *reset_ew_freq*
* value depends on the preceding keyword:

.. list-table::
   :widths: 20 20 60
   :header-rows: 1

   * - Keyword
     - Value
     - Description
   * - *dir*
     - path
     - Path to RuNNer configuration files
   * - *cflength*
     - factor
     - Length unit conversion factor
   * - *cfenergy*
     - factor
     - Energy unit conversion factor
   * - *committee_size*
     - N
     - Number of committee members
   * - *f_comm*
     - *yes* or *no*
     - Write individual committee forces to per-atom array
   * - *q_comm*
     - *yes* or *no*
     - Write individual committee charges to per-atom array
   * - *total_charge*
     - Q
     - Total system charge for 4G/QEq
   * - *use_prev_q*
     - *yes* or *no*
     - Use charges from previous step as QEq guess
   * - *check_extrap*
     - *yes* or *no*
     - Enable extrapolation warning (EW) monitoring
   * - *max_extrap*
     - N
     - Stop simulation if EW count exceeds N
   * - *show_ew*
     - *yes* or *no*
     - Write EWs to the log file
   * - *sum_ew_freq*
     - N
     - Write EW summary every N steps
   * - *reset_ew_freq*
     - N
     - Reset EW counters every N steps

Examples
--------

.. code-block:: LAMMPS

   pair_style runner dir "./potential_files"
   pair_coeff * * 1 8

   fix 1 all property/atom d2_f_comm 24 ghost yes
   pair_style runner dir "./potential_files" cflength 1.889726124626 &
      cfenergy 0.036749322175655 committee_size 8 f_comm yes
   pair_coeff * * 1 3 8 25

   fix 1 all property/atom d2_q_comm 4 ghost yes
   pair_style runner dir "./potential_files" committee_size 4 q_comm yes total_charge 0.0
   pair_coeff * * 1 8 78

Description
-----------

This pair style provides an interface to the `RuNNer 2 <https://gitlab.com/runner-suite/runner2>`_ (Ruhr University Neural Network Energy Representation) library. It implements High-Dimensional Neural Network Potentials (HDNNPs) as introduced in :ref:`(Behler and Parrinello 2007) <Behler_Parrinello_2007>`. HDNNPs are machine learning potentials that represent the total energy of a system as a sum of environment-dependent atomic contributions.



The pair style supports several "generations" of HDNNPs as proposed in :ref:`(Behler 2021) <Behler_2021>`:

* **Second-generation (2G):** Short-range many-body potentials where the total energy is the sum of atomic energies predicted from local chemical environments :ref:`(Behler and Parrinello 2007) <Behler_Parrinello_2007>`.
* **Third-generation (3G):** Extends 2Gs by adding explicit long-range electrostatic interactions based on environment-dependent atomic charges :ref:`(Artrith, Morawietz and Behler 2011) <Artrith_Morawietz_Behler_2011>`.
* **Fourth-generation (4G):** Includes global charge equilibration (QEq) based on environment-dependent electronegativities (and optionally hardness). These charges are fed back into the energy model, providing a global descriptor for the atomic energies :ref:`(Ko et al 2021) <Ko_Finkler_Goedecker_Behler_2021>`.

Additionally, all generations can be augmented with:

* **Hirshfeld-based dispersion:** Long-range dispersion interactions based on the Tkatchenko-Scheffler dispersion model :ref:`(Tkatchenko and Scheffler 2009) <Tkatchenko_Scheffler_2009>`.
* **Repulsive potentials:** Ziegler-Biersack-Littmark-based short-range pairwise repulsive potential.

Only a single :doc:`pair_coeff <pair_coeff>` command with two asterisk wildcards is used with this pair style. Its additional arguments define the mapping of LAMMPS atom types to RuNNer atomic numbers.

.. code-block:: LAMMPS

   pair_coeff * * 1 8

The example above maps LAMMPS atom types 1 and 2 to atomic numbers 1 ("H") and 8 ("O") in RuNNer.

----

General
^^^^^^^

Use the *dir* keyword to specify the directory containing the RuNNer configuration files. The directory must contain ``input.nn`` with the HDNNP architecture and feature map info, ``scaling_?.data`` with feature map scaling data, and ``weights_?.???.data`` with parameters for each element.

The RuNNer library is unit-agnostic. Use *cflength* and *cfenergy* to scale LAMMPS coordinates and energies to the units in which the potential was trained. If the HDNNP was trained in Bohr and Hartree and the LAMMPS simulation uses *metal* units (Angstroms, eV), then *cflength* and *cfenergy* must be the multiplicative factors required to convert LAMMPS units to the respective quantities in native HDNNP units:

.. code-block:: LAMMPS

   cflength 1.8897261328
   cfenergy 0.0367493254

Since machine learning potentials are most reliable within their training data range, the *runner* pair style can monitor whether the features representing the local atomic environments extrapolate beyond their training range. Set *check_extrap* to *yes* to enable monitoring. The keyword *show_ew* enables the writing of extrapolation warnings (EWs) to the log. With *sum_ew_freq*, you can specify whether a summary should be written at specific intervals instead of writing each EW to the log as it occurs. The *max_extrap* threshold allows termination of a simulation when the provided number of EWs is exceeded. Setting *max_extrap* to a negative number disables the termination threshold. With *reset_ew_freq*, the EW counters can be reset at specific intervals.

Committees
^^^^^^^^^^

The pair style supports **Committees**, where multiple HDNNPs sharing atomic descriptors are evaluated simultaneously. The forces, energies, and virials used to propagate the simulation are the average of all committee members. This is useful for Query-by-Committee-based Active Learning approaches and uncertainty estimation for production simulations.

In the case of a *committee_size* greater than 1, *dir* must point to a directory that contains ``input.nn`` and ``scaling_?.data``, as well as *committee_size* many subdirectories named 1 to *N*, which contain the ``weights_?.???.data`` of the respective members.

.. code-block:: text

   dir/
   ├── 1/
   │   ├── weights_short.001.data
   │   └── weights_short.008.data
   ├── 2/
   │   ├── weights_short.001.data
   │   └── weights_short.008.data
   ├── input.nn
   └── scaling.data

The individual potential energies of each committee member can be accessed using the :doc:`compute pair <compute_pair>` command:

.. code-block:: LAMMPS

   compute e_comm all pair runner

Here, the energies are stored in a global vector *e_comm* of length *committee_size*. They can be accessed as follows:

.. code-block:: LAMMPS

   c_e_comm[1] # total energy member 1
   c_e_comm[2] # total energy member 2
   c_e_comm[N] # total energy member N

The individual forces of the committee members (e.g., to compute force variance per atom) can be accessed by defining a custom per-atom array using the :doc:`fix property/atom <fix_property_atom>` command **before** the *pair_style* command. The array **must** be a floating-point array (type ``d2``) with the name *f_comm* and 3 times *committee_size* columns. It is necessary for ghost atom info to be communicated.

.. code-block:: LAMMPS

   # Define custom per-atom force array for committee_size 8 (8 * 3 = 24 columns)
   fix 1 all property/atom d2_f_comm 24 ghost yes

The committee forces are stored sequentially and can be accessed as:

.. code-block:: LAMMPS

   d2_f_comm[1] # fx member 1
   d2_f_comm[2] # fy member 1
   d2_f_comm[3] # fz member 1
   d2_f_comm[4] # fx member 2

For 3G and 4G potentials, *q_comm* can be set to *yes* to extract individual member charges. A custom per-atom array *q_comm* must be specified **before** the *pair_style* command. The array **must** be a floating-point array (type ``d2``) with the name *q_comm* and *committee_size* columns. Ghost atom communication must be enabled.

.. code-block:: LAMMPS

   # Define custom per-atom charge array for committee_size 2
   fix 2 all property/atom d2_q_comm 2 ghost yes 

The committee charges can be accessed as:

.. code-block:: LAMMPS

   d2_q_comm[1] # q member 1
   d2_q_comm[2] # q member 2


3G / 4G only
^^^^^^^^^^^^

For 3G and 4G HDNNPs, the total charge of the system can be specified with the *total_charge* keyword. For periodic systems, the system must be charge neutral.

For 4G HDNNPs only, when setting *use_prev_q* to *yes*, the predicted charges from the previous time step are used as an initial guess for QEq in the current time step.

.. note::

   3G and 4G HDNNPs require either full 3D periodicity or no periodicity. Partial periodicity (e.g., ``boundary p p f``) is not supported for long-range electrostatics.

.. note::

   Long-range electrostatics in RuNNer currently require global structure 
   collection. To avoid MPI bottlenecks, it is highly recommended to use 
   **OpenMP threading** (few MPI tasks, many OpenMP threads).

Mixing, shift, table, tail correction, restart, rRESPA info
-----------------------------------------------------------

This style does not support mixing. The :doc:`pair_coeff <pair_coeff>` command should only be invoked with asterisk wildcards (as shown above).

This style does not support the :doc:`pair_modify <pair_modify>` shift, table, and tail options.

This style does not write information to :doc:`binary restart files <restart>`. Thus, you must re-specify the *pair_style* and *pair_coeff* commands in any input script that reads a restart file.

This style can only be used via the *pair* keyword of the :doc:`run_style respa <run_style>` command. It does not support the *inner*, *middle*, or *outer* keywords.

Restrictions
------------

This pair style is part of the ML-RUNNER package. It is only enabled if LAMMPS was built with that package. See the :doc:`Build package <Build_package>` page for more info.

Currently, only one instance of ``pair_style runner`` can be initialized per simulation.

Related commands
----------------

:doc:`pair_coeff <pair_coeff>`, :doc:`fix property/atom <fix_property_atom>`, :doc:`compute pair <compute_pair>`

Default
-------

The default options are *dir* = "./", *cflength* = 1.0, *cfenergy* = 1.0, *committee_size* = 1, *f_comm* = no, *q_comm* = no, *total_charge* = 0.0, *use_prev_q* = no, *check_extrap* = no, *max_extrap* = 100, *show_ew* = no, *sum_ew_freq* = 0, and *reset_ew_freq* = 0.

---

References
----------

.. _Behler_Parrinello_2007:

**(Behler and Parrinello 2007)** Behler, J.; Parrinello, M. Phys. Rev. Lett. 2007, 98 (14), 146401.

.. _Tkatchenko_Scheffler_2009:

**(Tkatchenko and Scheffler 2009)** Tkatchenko, A.; Scheffler, M., Phys. Rev. Lett. 2009, 102, 073005.

.. _Behler_2011:

**(Behler 2011)** Behler, J., J. Chem. Phys. 2011, 134, 074106.

.. _Artrith_Morawietz_Behler_2011:

**(Artrith, Morawietz and Behler 2011)** Artrith, N.; Morawietz, T.; Behler, J., Phys. Rev. B 2011 83, 153101.

.. _Ko_Finkler_Goedecker_Behler_2021:

**(Ko et al 2021)** Ko, T. W.; Finkler, J. A.; Goedecker, S.; Behler, J, Nat. Commun. 2021 12, 398.

.. _Behler_2021:

**(Behler 2021)** Behler, J., Chem. Rev. 2021, 121, 16, 10037–10072.
