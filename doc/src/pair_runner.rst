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
          Length unit conversion factor
    *cfenergy* value = energy
          Unit conversion factor

    *committee_size* value = number
          number = Number of committee members
    *f_comm* value = *yes* or *no* 
          Write individual committee forces to a custom per-atom array
    *q_comm* value = *yes* or *no* 
          Write individual committee charges to a custom per-atom array

    *total_charge* value = charge
          charge = Total system charge for electrostatics/QEq
    *use_prev_q* value = *yes* or *no* 
          Use predicted charges from previous time step as initial guess for QEq

    *check_extrap* value = *yes* or *no*
          Enable monitoring of feature extrapolation
    *max_extrap* value = threshold
          threshold = Stop simulation if EW count exceeds this value
    *show_ew* value = *yes* or *no*
          Write EWs to the log
    *sum_ew_freq* value = summary
          summary = Write EW summary every this many time steps
    *reset_ew_freq* value = frequency
          frequency = Reset EW counters every this many steps

Examples
""""""""

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
"""""""""""

This pair style provides an interface to the `**RuNNer 2** <https://gitlab.com/runner-suite/runner2>` (Ruhr University Neural Network Energy Representation) library. It implements High-Dimensional Neural Network Potentials (HDNNPs) as introduced in :ref:`(Behler and Parrinello 2007) <Behler_Parrinello_2007>`. HDNNPs are machine learning potentials, which represent the total energy of a system as a sum of environment-dependent atomic contributions.

The pair style supports several "generations" of HDNNPs as proposed in :ref:`(Behler 2021) <Behler_2021>`.

* **Second-generation (2G):** Short-range many-body potentials where the total energy is the sum of atomic energies predicted from local chemical environments :ref:`(Behler and Parrinello 2007) <Behler_Parrinello_2007>`.
* **Third-generation (3G):** Extends 2Gs by adding explicit long-range electrostatic interactions based on environment-dependent atomic charges :ref:`(Artrith, Morawietz and Behler 2011) <Artrith_Morawietz_Behler_2011>`.
* **Fourth-generation (4G):** Includes global charge equilibration (QEq) based on environment-dependent electronegativities (and optionally hardness). These charges are fed back into the energy model, providing a global descriptor for the atomic energies :ref:`(Ko et al 2021) <Ko_Finkler_Goedecker_Behler_2021>`.

Additionally, all generations can be augmented with:

* **Hirshfeld-based dispersion:** Long-range dispersion interactions based on the Tkatchenko-Scheffler dispersion model :ref:`(Tkatchenko and Scheffler 2009) <Tkatchenko_Scheffler_2009>`.
* **Repulsive potentials:** Ziegler-Biersack-Littmark-based short-range pairwise repulsive potential.

Only a single *pair_coeff* command with two asterisk wild-cards is used with this
pair style. Its additional arguments define the mapping of LAMMPS atom types to
RuNNer atomic number.

.. code-block:: LAMMPS

   pair_coeff * * 1 8

maps atom types 1 and 2 to the atomic number 1 ("H") and 8 ("O") in RuNNer.

----

General
"""""""""""

Use the *dir* keyword to specify the directory containing the RuNNer configuration
files. The directory must contain ``input.nn`` with HDNNP architecture and feature map info, ``scaling_?.data`` with feature map scaling data and ``weights_?.???.data`` with parameters for each element.

The RuNNer library is unit-agnostic. Use *cflength* and *cfenergy* to scale LAMMPS coordinates and energies to the units the potential was trained in. If the HDNNP was trained in Bohr and Hartree and the LAMMPS simulation uses *metal* units (Angstroms, eV), then `cflength` and `cfenergy` must be the multiplicative factors required to convert LAMMPS units to the respective quantities in native HDNNP units.

.. code-block:: LAMMPS

   cflength 1.8897261328
   cfenergy 0.0367493254

Since machine learning potentials are most reliable within their training data range, the *runner* pair style can monitor if the features, representing the local atomic environments, extrapolate beyond their training range. Set *check_extrap* to *yes* to toggle on monitoring. The keyword *show_ew* enables writing of EWs to the log. With *sum_ew_freq* it can be specified whether a summary should be written in specific intervals instead of writing the EW to the log when it occurs.
The *max_extrap* threshold allows termination of a simulation when the provided number of EWs is exceeded. Setting *max_extrap* to a negative number disables the termination threshold. With *reset_ew_freq* the EW counters can be reset in specific intervals.

Committees
"""""""""""

The pair style supports **Committees** where multiple HDNNPs, sharing atomic descriptors, are evaluated simultaneously. The forces, energies, and virials used to propagate the simulation are the average of all committee members.
This is useful for Query-by-Committee-based Active Learning approaches and uncertainty estimation of production simulations.
In the case of a *committee_size* greater than 1, *dir* must point to a directory, which contains ``input.nn`` and ``scaling_?.data`` and *committe_size* many directories named 1 to *committee_size*, which contain the ``weights_?.???.data`` of the respective member.

.. code-block:: LAMMPS

   # dir tree for 2G with committee_size 2
   - 1
     - weights_short.001.data
     - weights_short.008.data
   - 2
     - weights_short.001.data
     - weights_short.008.data
   - input.nn
   - scaling.data

The individual potential energies of each committee member can be accessed using the :doc:`compute pair <compute_pair>` command.

.. code-block:: LAMMPS

   compute e_comm all pair runner

Here, the energies are stored in the global *e_comm* vector of length *committee_size*. They can be accessed as

.. code-block:: LAMMPS

   c_ec[1] # total energy member 1
   c_ec[2] # total energy member 1
   c_ec[N] # total energy member N

The individual forces of the committee members (e.g., to compute force variance per atom) can be accessed by defining a custom per-atom array using the :doc:`fix property/atom <fix_property_atom>` command **before** the *pair_style* command.
The array **must** be a floating-point array *d2_name* with name *f_comm* and three times *committee_size* columns. It is necessary for ghost atom info to be communicated.

.. code-block:: LAMMPS

   fix 1 all property/atom d2_f_comm 24 ghost yes # define custom per-atom force array for committee_size 8

The committee forces are stored sequentially and can be accessed as 

.. code-block:: LAMMPS

   d2_f_comm[1] # fx member 1
   d2_f_comm[2] # fy member 1
   d2_f_comm[3] # fz member 1
   d2_f_comm[4] # fx member 2

For 3G and 4G potentials *q_comm* can be set to *yes* to extract individual member charges. 
A custom per-atom array *q_comm* needs to specified **before** the *pair_style* command. The array **must** be a floating-point array *d2_name* with name *q_comm* and *committee_size* columns. It is necessary for ghost atom info to be communicated.

.. code-block:: LAMMPS

   fix 2 all property/atom d2_q_comm 2 ghost yes # define custom per-atom charge array for committee_size 2

The committee charges can be accessed as

.. code-block :: LAMMPS

   d2_q_comm[1] # q member 1
   d2_q_comm[2] # q member 2


3G / 4G only
"""""""""""

For 3G and 4G HDNNPs, the total charge of the system can be specified with the *total_charge* keyword.  For periodic systems, the system must be charge neutral. 

For 4G HDNNPs only, when toggling *use_prev_q* to *yes* the predicted charges from the previous time step are used as initial guess for QEq of the current time step. 

.. note::
   3G and 4G HDNNPs require either full 3D periodicity or no periodicity. Partial periodicity e.g., ``boundary p p f`` is not supported for long-range electrostatics.

.. note::
   3G and 4G HDNNPs require global structure collection on a single process for the long-range electrostatics / global QEq calculations, which creates an MPI bottleneck and subotpimal scaling with the number of MPI tasks. 3G and 4G HDNNPs are heavily optimized for OpenMP in RuNNer. Therefore, it is highly recommended to use only a small number of MPI tasks for 3G and 4G HDNNPs and a large number of OpenMP threads per task.

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This style does not support mixing. The :doc:`pair_coeff <pair_coeff>` command
should only be invoked with asterisk wild cards (see above).

This style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This style does not write information to :doc:`binary restart files <restart>`.
Thus, you need to re-specify the pair_style and pair_coeff commands in an input
script that reads a restart file.

This style can only be used via the *pair* keyword of the :doc:`run_style respa
<run_style>` command.  It does not support the *inner*, *middle*, *outer*
keywords.

Restrictions
""""""""""""

This pair style is part of the ML-RUNNER package. It is only enabled if LAMMPS was built with that package. See the :doc:`Build package <Build_package>` doc page for more info.

Currently, only one instance of `pair_style runner` can be initialized per simulation.

Related commands
""""""""""""

:doc:`pair_coeff <pair_coeff>`, :doc:`fix property/atom <fix_property_atom>`, :doc:`compute pair <compute_pair>`

Default
""""""""""""

The default options are *dir* = "./", *cflength* = 1.0, *cfenergy* = 1.0, *committee_size* = 1, *f_comm* = no, *q_comm* = no, *total_charge* = 0.0, *use_prev_q* = no, *check_extrap* = no, *max_extrap* = 100, *show_ew* = no, *sum_ew_freq* = 0, *reset_ew_freq* = 0.

----

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
