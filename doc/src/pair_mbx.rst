.. index:: pair_style mbx

pair_style mbx command
======================

Syntax
""""""

.. code-block:: LAMMPS

    pair_style mbx cutoff

* cutoff = real-space cutoff for MBX. 9.0 is usually a safe value.


Examples
""""""""

.. code-block:: LAMMPS
    pair_style mbx 9.0
    pair_coeff      * * 0.0 0.
    compute         mbx all pair mbx


    pair_style      hybrid/overlay mbx 9.0 lj/cut 9.0 coul/exclude 9.0
    pair_coeff      * * mbx  0.0 0.0
    pair_coeff      1*11 1*11 coul/exclude  # 
    compute         mbx all pair mbx


Description
"""""""""""




For hybrid simulations involving MB-nrg and non-MB-nrg molecules in the
same simulation, one can use hybrid/overlay

Do note that all electrostatics must be computed within MBX, so the
coul/exclude pair_style must be applied on the non-MB-nrg molecules.
See  `mof <>` for a complete hybrid example. 


Restrictions
""""""""""""

This pair_style is part of the MBX package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info. This pair_style also relies on the
presence of `fix mbx <fix_mbx>` command.

Due to the usage of Partridge and Schwenke charges for MB-pol,
all electrostatic interactions are calculated internally in MBX.
Therefore one should avoid calculating coulombic interactions in
LAMMPS such as using `coul/cut` or `coul/long`.

