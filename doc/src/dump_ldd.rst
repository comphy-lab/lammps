.. index:: dump ldd

dump ldd command
================

Syntax
""""""

.. code-block:: LAMMPS

   dump ID group-ID ldd filename

* ID = user-assigned name for the dump
* group-ID = ID of the group of atoms to be dumped
* ldd = style of dump command (other styles such as *atom* or *cfg* or *dcd* or *xtc* or *xyz* or *local* or *custom* are discussed on the :doc:`dump <dump>` doc page)
* N = dump every this many timesteps
* file = name of file to write dump info to
* args = none`

Examples
""""""""

.. code-block:: LAMMPS

   dump 4 all ldd 500 dump.txt

Description
"""""""""""

An optional dump style for when :doc:`atom_style ldd <atom_style>` is used in the :ref:`LDD <PKG-LDD>` package.
Dump a snapshot of atom quantities to a file every :math:`N`
timesteps in a text format readable by the `Bottom Up Open Source Coarse Graining Package (BOCS)
<https://github.com/noid-group/BOCS>`_  for constructing CG models and translating trajectories to other file formats. 

The output of this dump style is similar to a :doc:`dump_style custom <dump>` but 
it includes per atom local density and gradient of the local density/ ldd energy information for each type of particle for each particle in the system. 
Each frame of the ldd dumped text file contains 4 ITEM sections: the timestep, the number of atoms, the box dimensions, and the ATOMS sections. 
The first three sections are 1-3 lines containing the time, total number of atoms in the system and box bounds respectively.

Each line in the atoms section contains per atom data fields in an order which is labeled by a headereach frame.

below is a table of ldd header field labels, definitions in the order they appear


+------------+-----------------------------------------------------------------+
| **label**  |  **Definition**                                                 |
+============+=================================================================+
| id         | The atom index, :math:`I`                                       |
+------------+-----------------------------------------------------------------+
| mol        | The molecule index (only appears is atom_style ldd is listed    |
|            | an atom_style hybrid list with a molecular atom_style)          |
+------------+-----------------------------------------------------------------+
| type       | The atom type index, :math:`t_I \in` {1 ... n_types}            |
+------------+-----------------------------------------------------------------+
| x y z      | The x y and z components of the atom's position                 |
+------------+-----------------------------------------------------------------+
| vx vy vz   | The x y and z components of the atom's velocity                 |
+------------+-----------------------------------------------------------------+
| fx fy fz   | The x y and z components of the atom's net force                |
+------------+-----------------------------------------------------------------+
| lddensn    | n = {1 ... ntypes}, the local density of particle type n around | 
|            | the atom,  :math:`\rho_{n|I}`                                   |
+------------+-----------------------------------------------------------------+
| ldnrgn     | n = {1 ... ntypes}, the energy contribution of the local        |
|            |  density of type n around this particle to the total energy     |
|            | :math:`u_{n|t_I}(\rho_{n|I})`                                   |
+------------+-----------------------------------------------------------------+
| gradxn     | n = {1 ... ntypes), the x y and z components of the gradient of |
| gradyn     | the local density of type n surrounding the particle.           |
| gradzn     | :math:`\nabla_{I} \rho_{n|I}`                                   |
+------------+-----------------------------------------------------------------+
| gradnrgn   | n = {1 ... ntypes}, the energy contribution of this particle    |
|            |  surrounded by type n to the square gradient potential.         |
|            | :math:`u_{\nabla; n|t_I}(\rho_{n|I})|\nabla_{I}\rho_{n|I}|^2`   |
+------------+-----------------------------------------------------------------+
| lddttlnrg  | The total energy contribution of this particle to LD and SG     |
|            | potentials.                                                     |
|            | :math:`sum_{\text{n}_{\text{types}}} \text{gradnrgn + ldnrgn}`  |
+------------+-----------------------------------------------------------------+


If local densities and square gradient potentials are not defined with :doc:`pair_style ldd <pair_ldd>` 
then the relevant LD and SG fields will be zeroed out in the output.


Only information for atoms in the specified group is dumped.

.. warning::

   Because periodic boundary conditions are enforced only
   on timesteps when neighbor lists are rebuilt, the coordinates of an
   atom written to a dump file may be slightly outside the simulation
   box.

.. warning::

   Unless the :doc:`dump_modify sort <dump_modify>` option is invoked,
   the lines of atom information written to dump files will be in an
   indeterminate order for each snapshot.  This is even true when
   running on a single processor, if the :doc:`atom_modify sort
   <atom_modify>` option is on, which it is by default.  In this case
   atoms are re-ordered periodically during a simulation, due to spatial
   sorting.  It is also true when running in parallel, because data for
   a single snapshot is collected from multiple processors, each of
   which owns a subset of the atoms.

For the *ldd* style, sorting is off by default. See the
:doc:`dump_modify <dump_modify>` page for details.

----------

Dumps are performed on timesteps that are a multiple of N (including
timestep 0) and on the last timestep of a minimization if the
minimization converges.  Note that this means a dump will not be
performed on the initial timestep after the dump command is invoked,
if the current timestep is not a multiple of N.  This behavior can be
changed via the :doc:`dump_modify first <dump_modify>` command, which
can also be useful if the dump command is invoked after a minimization
ended on an arbitrary timestep.  N can be changed between runs by
using the :doc:`dump_modify every <dump_modify>` command.
The :doc:`dump_modify every <dump_modify>` command
also allows a variable to be used to determine the sequence of
timesteps on which dump files are written.  In this mode a dump on the
first timestep of a run will also not be written unless the
:doc:`dump_modify first <dump_modify>` command is used.

Dump filenames can contain two wildcard characters.  If a "\*"
character appears in the filename, then one file per snapshot is
written and the "\*" character is replaced with the timestep value.
For example, tmp.dump\*.txt becomes tmp.dump0.txt, tmp.dump10000.txt,
tmp.dump20000.txt, etc.  Note that the :doc:`dump_modify pad <dump_modify>`
command can be used to ensure all timestep numbers are the same length
(e.g. 00010), which can make it easier to read a series of dump files
in order with some post-processing tools.

If a "%" character appears in the filename, then each of P processors
writes a portion of the dump file, and the "%" character is replaced
with the processor ID from 0 to P-1 preceded by an underscore character.
For example, tmp.dump%.txt becomes tmp.dump_0.txt, tmp.dump_1.txt, ...
tmp.dump_P-1.txt, etc.  This creates smaller files and can be a fast
mode of output on parallel machines that support parallel I/O for output.

By default, P = the number of processors meaning one file per
processor, but P can be set to a smaller value via the *nfile* or
*fileper* keywords of the :doc:`dump_modify <dump_modify>` command.
These options can be the most efficient way of writing out dump files
when running on large numbers of processors.

Note that using the "\*" and "%" characters together can produce a
large number of small dump files!

----------

Restrictions
""""""""""""


The *ldd* dump style is part of the LDD package. It is only
enabled if LAMMPS was built with that package. See the :doc:`Build package <Build_package>` page for more info.

The *ldd* dump style is only supported when the atom_style ldd is used. See :doc:`Howto_ldd <Howto_ldd>` for more details. 

The *ldd* dump style supports neither buffering or custom format
strings.

Related commands
""""""""""""""""

:doc:`dump <dump>`, :doc:`LDD <Howto_ldd>`, :doc:`pair_style ldd <pair_ldd>`, 
:doc:`dump_modify <dump_modify>`, :doc:`undump <undump>`

Default
"""""""

By default, files are written in ASCII format.
