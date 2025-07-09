Local Density Dependent Potentials
==================================

**Overview:**

The LDD package offers local density dependent potentials that are either 
a function of the local density (LD) around a central particle
or a function of the square gradient (SG) of the local density around a central particle.
Such potentials are useful for coarse graining and can be constructed from trajectory data using
`BOCS version 5 and higher <https://github.com/noid-group/BOCS>`_.

**Local Density Definitions**

The local density, :math:`\rho_{\beta|I}` of site types :math:`\beta` around a central particle of type :math:`t_I = \alpha` is defined:

.. math::

   \rho_{\beta|I} = \sum_{J\in S_{\beta|I}} \bar{w}_{\beta|\alpha}(r_{IJ})
 
where :math:`S_{\beta|I}` is the set of particles of type :math:`\beta` not excluded in pairwise interactions involving site I. 
:math:`\bar{w}_{\beta|\alpha}(r_{IJ})` is a normalized indicator function that determines the weight of each particle of type :math:`\beta` that contributes to the :math:`beta` LD around particular particle I of type :math:`t_I = \alpha`.
It is defined by dividing a continuous/differentiably nonincreasing function :math:`w(r)` by its spatial integral.

.. math::

   \bar{w}_{\beta|\alpha}(r_{IJ}) = \frac{w(r)}{\int_{0}^{r_{c}} {w(r)}}


where in the above :math:`w(r=0) = 1` and :math:`w(r=r_c) = 0`.
The local density of particles of type :math:`\beta` surrounding a given particle I is then the sum of particles of type :math:`\beta` that are within (inclusive of) the cutoff :math:`r_c`.
If particle I is of type :math:`t_I = \beta`, it's contribution to the local density can be included or excluded from the sum without changing the net force on a pair of particles. 
Practically, the only difference is a horizontal shift in the LD potential :math:`u_rho`

For molecular systems, particles included in the local density around a given particle I will be based on the exclusion list for the intramolecular nb interactions. 
For example, consider a 5 beaded chain simulated with standard 1-4 exclusions: 

:math:`A_{1}-B_{2}-C_{3}-D_{4}-E_{5}` 

The local density of particles of type E surrounding site 1, :math:`\rho_{E|1}` will include 
the contribution from the end chain particle :math:`\bar{w}_{E|A}(r_{15})`. 
Conversely, the local density of particles of type B surrounding site 1, :math:`\rho_{B|1}` will only include intermolecular contributions to the local density.
The local density of particles of Type A surrounding site 1, :math:`\rho_{A|1}` will include the *self* intramolecular contribution :math:`\bar{w}_{A|A}(0)` based on whether *self yes* or *self no* is listed by the user in the ldd pair_coeff command. 

The choice of LD definition used for a particular set of types :math:`\beta|\alpha` (not nec eq to :math:`\alpha|\beta` ) interaction is set by the user in via the ldd pair_coeff command.

**Local Density/ Square Gradient Potentials**

The LDD package implements forces from two kinds of local density dependent potentials. 

LD potentials defined by:

.. math::
   
   U_{\rho}(\mathbf{R}) = \sum_{I} \sum_{\beta} u_{\beta|t_I}(\rho_{\beta|I})

where the first sum is over all sites, I, and the second sum is over all site types, :math:`\beta`.
The form of :math:`u_{\beta|\alpha}` is set by the user in the ldd pair_coeff command using the *potential* keyword. 

SG potentials are defined by:

.. math::

   U_{SG}(\mathbf{R}) = \sum_{I} \sum_{\beta} u_{\nabla; \beta|t_I} (\rho_{\beta|I}) \left | \nabla_{I} \rho_{\beta|I} \right|^2

where the first sum is over all sites, I, the second sum is over all site types, :math:`\beta` and the gradient :math:`\nabla_{I} = \frac{d}{d\mathbf{R_I}}` is with respect to the position of particle I. 
The form of :math:`u_{\nabla; \beta|t_I}` is specified by the user in the ldd pair_coeff command using the (optional) *gradient* keyword. 


For any given site I, once the local densities :math:`\rho_{\beta|I}` and the gradients of the local density :math:`\nabla_{I} \rho_{\beta|I}` have been determined, the corresponding forces are pairwise decomposable with pair forces analagous to those stated in `Delyser 2021 <https://pubs.aip.org/aip/jcp/article/156/3/034106/2839866/Coarse-grained-models-for-local-density-gradients>`_. 

Because of this, the LDD package must be used with its own custom :doc:`atom_style <atom_style>`, *ldd*, where a list of all local densities and gradients of local densities are added as a per_atom data field to the most basic atomic atom_style. 
Also because of this, interactions may be specified in a pair-typewise manner using the regular pair_coeff commands, where in this instance e.g. pair_coeff 1 2 is treated distinctly from pair_coeff 2 1. 


**Putting it all together, input Overview**

Examples
"""""""""

A simple atomic input example using only tabulated LDD potentials

.. code-block:: LAMMPS
 
   ## Init system/units
   units real
   dimension 3
   atom_style ldd 2  # The total number of particle types

   newton on
   timestep 1

   read_data my_mix.data # System initialization, must pre-ceed pair style init
   velocity        all create 300.0 22345 dist gaussian
    
   ## pair_style ldd, must be used with atomstyle ldd

   pair_style ldd 6.5 # Longest cutoff of all LD interactions

   pair_coeff 1 1 indicator dpd 0.0 6.5 self yes potential table/lin LD_table.1.1.dat # type 1 surrounded by type 1

   pair_coeff 1 2 ignore # type 1 surrounded by type 2 
   # n.b. All 2^ntypes interactions must be specified even if its to ignore them see keyword ignore for details in pair_style doc
   pair_coeff 2 1 indicator dpd 0.0 5.5 self no potential table/lin LD_table.2.1..dat # type 2 surrounded by type 1

   pair_coeff 2 2 indicator dpd 0.0 5.5 self yes potential table/line LD_table.2.2.dat # type 2 surrounded by type 2

   ## Run / Output
   run_style verlet
   neighbor 15.0 bin
   thermo 500

   fix 1 all nvt temp 300.0 300.0 100.0 
   dump 1 all ldd 500 dump.txt #A lammps trajectory file with LD info

   run 10000

Note in the above example, only pair types 12, 11 and 22 will interact according to LD potentials. 


Molecular input example that layers tabulated pair, LD and SG interactions via pair_style hybrid

.. code-block:: LAMMPS

   ## Init system /units

   units real
   dimension 3
   atom_style hybrid full ldd 2 # Molecular topology accessible through hybrid atom_style

   newton on
   timestep 1

   bond_style harmonic 

   region my_box block 0 40 0 40 0 40
   create_box 2 my_box bond/types 1 # System init, must pre-ceed pair style ldd

   mass 1 59.044
   mass 2 59.044

   pair_style hybrid/overlay ldd 6.5 table spline 3000

   molecule DI DIMER.mol 
   create_atoms 0 random 250 885401 my_box mol DI 454353

   include ff_bonded_params.lammps
   velocity all create 300.0 22345 dist gaussian

   ## Pair styles

   pair_coeff * * table lammps_nb_ALL.table nb_All 15.0 # Pair interaction cutoff is independent of LDD package cutoffs

   pair_coeff 1 1 ldd indicator lucy 0.0 6.5 self no potential table/spline LD.1.1.dat gradient table/gradspline SG.1.1.dat
   
   pair_coeff 1 2 ldd indicator lucy 0.0 6.5 self no potential table/spline LD.1.2.dat
   pair_coeff 2 1 ldd indicator lucy 0.0 6.5 self no potential table/spline LD.2.1.dat

   pair_coeff 2 2 ldd indicator lucy 0.0 6.5 self no potential table/spline LD.2.2.dat gradient table/gradspline SG.2.2.dat

   ## Run / Output

   run_style verlet
   neighbor 15.0 bin
   thermo 500

   fix 1 all nvt temp 300.0 300.0 100.0 
   dump 1 all ldd 500 dump.txt #A lammps trajectory file with LD info

   run 10000
 
Note in the above example LD cross interactions are both specified, but not necessarily symetric.  
Also cross interactions opt not to use the optional gradient keyword. 
Conversely homo-interactions (e.g. 11 22) layer SG potentials on top of the already defined LD and pair interactions. 

See :doc:`pair_style ldd <pair_ldd>` for all pair_coeff *args* options that exist, all restrictions, and more examples.


**atom_style ldd read_data input file format**

Atom style ldd is a basic atomic atom_style with per-atom fields added for local densities, gradients of local densities, LD energy contributions and SG energy contributions. 
These can be reported using :doc:`dump ldd <dump_ldd>`, but are not required for starting simulations. 
read_data input therefore can follow usual atomic read_data input formats, and when hybridized with other atom styles, the .data follw is the same as atomic hybridized with those styles.

The following is an example for when the atom_style ldd is used by itself

.. code-block:: LAMMPS

   Topology of CG water translated from gromacs input

   5000 atoms
   0 bonds
   0 angles
   0 dihedrals
   0 impropers

   1 atom types
   0 bond types
   0 angle types
   0 dihedral types
   0 improper types

   0.000000 53.222800 xlo xhi
   0.000000 53.222800 ylo yhi
   0.000000 53.222800 zlo zhi

   Masses

   1 18.0154 # SOL

   Atoms

   1 1 46.050000 19.400000 9.760000 #atidx typeidx x y z
   2 1 7.860000 36.680000 8.090000
   3 1 42.340000 39.600000 8.310000
   4 1 22.950000 37.630000 5.880000

   .
   .
   .
   .
   4998 1 15.120000 8.330000 19.900000
   4999 1 5.300000 42.540000 30.500000
   5000 1 44.600000 18.310000 51.990000


The following is an example for the Atoms section of the read_data file when the atom_style ldd is hybridized with full. Bonds/Angles/Dihedral syntax are standard as listed in :doc:`read_data <read_data>`. 

.. code-block:: LAMMPS

   LAMMPS data file via write_data, version 12 Jun 2025, timestep = 985, units = real

   500 atoms
   4 atom types
   375 bonds
   1 bond types
   250 angles
   1 angle types
   125 dihedrals
   1 dihedral types

   0 40 xlo xhi
   0 40 ylo yhi
   0 40 zlo zhi

   Masses

   1 59.0448
   2 59.0448
   3 59.0448
   4 59.0448

   Atoms # hybrid

   19 3 6.337397372531347 1.5408419900080126 7.155488269009223 5 0 0 0 0 #id type x y z molid q nx ny nz
   20 4 4.504049144010551 1.3471156110429994 9.244879846603403 5 0 0 0 0
   129 1 6.919214369251321 4.484508069873466 10.334931011393717 33 0 0 0 0
   130 2 4.485538449090976 5.13320081653011 11.674419677953514 33 0 0 0 0
   131 3 4.2382325348268 7.693584581555052 10.175944460105603 33 0 0 0 0
   132 4 2.4530920970231946 8.798676084329882 12.317895781012147 33 0 0 0 0

   .
   .
   .



**dump ldd output**

Full LD/SG simulation statistics for each pair of types where an LDD type interaction can be calculated/dumped in the ldd dump command, which is essentially a custom lammps dump trajectory output with local density information. 

See :doc:`dump ldd <dump_ldd>` for details.

This trajectory type is natively compatible with the `Bottom-up Open-source Coarse-graining Software <https://github.com/noid-group/BOCS>`_ which can be used to construct LD/SG potentials from atomistic data, as well as to convert lammps trajectories to .trr files for analysis with gromacs(2019.6) tools.


**When to use LD/SG potentials:**

This package has been primarily developed in the context of constructing/simulating 
bottom up coarse grained models. 
We have found that LD/SG potentials are particularly useful for simulating CG models in the NPT 
and in different interfacial environments. 
SG potentials are slower to simulate than LD potentials alone, so where applicable for coarse graining we recommend trying to construct/simulate LD potentials first and then [if necessary] adding SG potentials to refine further. 
A number of works have been published using this package in its dev stage for different CG studies. 

**References**

.. _DeLyser1:

**Michael R. DeLyser, W. G. Noid. "Extending Pressure-Matching to Inhomogenous Systems via Local-Density Potentials." The Journal of Chemical Physics, 147, no. 13: 134111 (2017)**


.. _DeLyser2:

**Michael R. DeLyser, W. G. Noid. "Analysis of local density potentials." The Journal of Chemical Physics, 151, no. 22:224106 (2019)**

.. _DeLyser3:
   
**Michael R. DeLyser, W. G. Noid "Bottom-up coarse-grained models for external fields and interfaces" The Journal of Chemical Physics 153, 224103 (2020)**

.. _DeLyser4:

**Michael R. DeLyser, W. G. Noid "Coarse-grained models for local density gradients" The Journal of Chemical Physics 156, 034106 (2021)**


.. _Szukalo:

**Ryan K/ Szukalo, W.G. Noid "A temperature-dependent length-scale for transferable local density potentials". The Journal of Chemical Physics, 159, 074104 (2023)**

.. _Lesniewski1:

**Maria C. Lesniewski, R. C. Remsing, W. G. Noid "Coarse Graining the Hydrophobic Effect." The Journal of Chemical Physics, X, X (2025)**

.. _Dutta:

**Sayan Dutta, Muhammad Nawaz Qaisrani, Maria C. Lesniewski, W. G. Noid, Denis Andrienko, Arash Nikoubashman "Many body effect and optimized mapping scheme for systematic coarse-graining". X, X, X, (2025)**

.. _Lesniewski2:

**Maria C. Lesniewski, Michael R. DeLyser, Nicholas J. H. Dunn, Joseph R. Rudzinski, Ryan J. Szukalo, W. G. Noid "BOCS version 5: Bottom Up Coarse Graining with Local Densities" X, X, X (2025)**


