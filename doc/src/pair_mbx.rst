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

    # For a system involving water (atom types O=1, H=2)
    pair_style mbx 9.0
    pair_coeff      * * 0.0 0.0
    compute         mbx all pair mbx
    fix             mbx_fix all mbx 1 h2o 1 2 3 1 2 2 json mbx.json


    # For a system involving ch4 (atom types C=1, H=2) and
    # water (atom types O=3, H=4)
    pair_style mbx 9.0
    pair_coeff      * * 0.0 0.0
    compute         mbx all pair mbx
    fix             mbx_fix all mbx 2 ch4 1 2 5 1 2 2 2 2 h2o 3 4 3 3 4 4 json mbx.json

    # For a system involving water (atom types 0=12, H=13) in a hybrid simulation
    pair_style      hybrid/overlay mbx 9.0 lj/cut 9.0 coul/exclude 9.0
    pair_coeff      * * mbx  0.0 0.0
    pair_coeff      1*11 1*11 coul/exclude
    compute         mbx all pair mbx
    fix             mbx_fix all mbx 2 dp1 1 11 1 1 h2o 12 13 3 12 13 13 json mbx.json

See ``examples/PACKAGES/mbx`` for additional examples of how to use MBX in LAMMPS.


Description
"""""""""""

This pair_style instructs LAMMPS to call the `MBX library <_mbxwebsite>`_
in order to simulate MB-nrg models such as MB-pol. This pair_style must be
used in conjunction with the :doc:`fix mbx <fix_mbx>` command.


The MBX library code development is available at
`https://github.com/paesanilab/MBX <https://github.com/paesanilab/MBX>`_.
A detailed discussion of the code can be found in the manuscript :ref:`(Riera) <Riera>`.

If you have questions not answered by this documentation, please reach
out to us at `https://groups.google.com/g/mbx-users <https://groups.google.com/g/mbx-users>`_


See ``examples/PACKAGES/mbx`` for complete examples of how to use
this fix command.

For hybrid simulations involving MB-nrg and non-MB-nrg molecules in the
same simulation, one can use :doc:`pair_style hybrid/overlay <pair_hybrid_overlay>`
to combine

For hybrid simulations involving MB-nrg and non-MB-nrg molecules in the
same simulation, one can use :doc:`pair_style hybrid/overlay <pair_hybrid_overlay>`
to combine the MB-nrg molecules with other pair styles, such as
:doc:`lj/cut <pair_lj_cut>`. Do note that all electrostatics must be computed within MBX, so the
coul/exclude pair_style must be applied on the non-MB-nrg molecules.
See  ``examples/PACKAGES/mbx`` for a complete hybrid example. 


Restrictions
""""""""""""

This pair_style is part of the MBX package.  It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info. This pair_style also relies on the
presence of :doc:`fix mbx <fix_mbx>` command.

Due to the usage of Partridge and Schwenke charges for MB-pol,
all electrostatic interactions are calculated internally in MBX.
Therefore one should never calculate coulombic interactions in
LAMMPS such as using `coul/cut` or `coul/long`.



Related commands
""""""""""""""""

:doc:`fix mbx <fix_mbx>`
:doc:`pair_style hybrid/overlay <pair_hybrid_overlay>`
:doc:`coul/exclude <pair_coul_exclude>`


-----------

.. _Riera:

**(Riera)** M. Riera, C. Knight, E. Bull-Vulpe, X. Zhu, H. Agnew, D. Smith, A. Simmonett, F. Paesani, J. Chem. Phys. 159, 054802 (2023)

.. _mbxwebsite: https://mbxsimulations.com
