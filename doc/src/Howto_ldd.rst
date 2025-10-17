Local Density Dependent Potentials via PKG-LDD
===============================================
**Overview:**

The LDD package implements potentials that are a function of the local density, :math:`\rho_{\beta|I}`, of type :math:`\beta` particles around a central particle, :math:`I`.
There are two types of local density dependent potentials implemented.
There are Local Density (LD) potentials :math:`U_{LD}` , which are direct functions of :math:`\rho_{\beta|I}`.
There are also Square Gradient (SG) potentials :math:`U_{SG}`.
SG potentials use a direct function of the local density, :math:`U_\nabla(\rho_{\beta|I})`, as a coefficient of the square gradient, :math:`|\frac{\partial \rho_{\beta|I}}{\partial R_I}|^2` of a particular local density.
Such potentials are useful for coarse grained models and can be parameterized from trajectory data using
`BOCS version 5 and higher. <https://github.com/noid-group/BOCS>`_

**Local Density Definitions**

Each site (atom) in a molecular dynamics simulation :math:`I` is assigned a type :math:`t_I`.
These site types may be labeled e.g. :math:`\alpha`, :math:`\beta`, ..., etc.
The local density, :math:`\rho_{\beta|I}`, of :math:`\beta` particles around a central particle, :math:`I`, of type :math:`t_I = \alpha` is defined:

.. math::

   \rho_{\beta|I} = \sum_{J\in S_{\beta|I}} \bar{w}_{\beta|\alpha}(r_{IJ})

where :math:`S_{\beta|I}` is the set of particles of type :math:`\beta` that are not excluded from pairwise interactions involving site :math:`I`.
:math:`\bar{w}_{\beta|\alpha}(r_{IJ})` is a normalized indicator function, it determines the pair weighted contribution from each particle :math:`J`, of type :math:`t_J = \beta` to the :math:`\beta` local density around a central particle :math:`I` of type :math:`t_I = \alpha`.

The normalized indicator function is defined by dividing a continuous/differentiable non-negative non-increasing indicator function :math:`w(r)` by its spatial integral :math:`[w]`.
Therefore, :math:`\bar{w}(r)` is defined for each :math:`\beta|\alpha` pair as:

.. math::

   \bar{w}(r) = \frac{w(r)}{[w]} \theta(r - r_c) \\
   [w] = \int_{0}^{r_{c}} {4\pi r^2 w(r)} dr


where :math:`\theta` is the Heaviside function, :math:`r_c` is the local density length-scale, and :math:`\bar{w}(r)` is defined such that :math:`\bar{w}(r) = 0`, if  :math:`r > r_c`.
The local density of :math:`\beta` particles that surround a given particle :math:`I` is defined by the sum of contributions from :math:`\beta`  particles that are within (<=) the cutoff :math:`r_c`.
The choice of indicator function is specified by the user in the *indicator* argument of the :doc:`ldd pair_coeff <pair_ldd>` command. Different options are listed under different pages indexed on the :doc:`pair ldd <pair_ldd>` page.
Note that the weighting function and its length-scale can be defined independently for each pair of site types.

If particle :math:`I` is of type :math:`t_I = \beta`, its contribution to the local density :math:`\rho_{t_I|I}` can be included or excluded from the sum without changing the net force on a pair of particles.
Practically, the only difference is a horizontal shift in the argument of the LD potential function :math:`U_{\beta|\beta} (\rho)`.
Whether or not to include this contribution is set by the user in the :doc:`ldd pair_coeff <pair_ldd>` command via the *self yes* or *self no* arguments.


For molecular systems, particles included in the local density around a given particle :math:`I` are identified by the neighbor list for the non-bonded interactions.
That is, a particle excluded by the neighbor list for a central particle will also be excluded from the sum defining its local density.

For example, consider a 5 bead chain that is simulated with default 1-2, 1-3, and 1-4 non-bonded exclusions:

:math:`A_{1}-B_{2}-C_{3}-D_{4}-E_{5}`

The local density of particles of type E surrounding site 1, :math:`\rho_{E|1}` will include
the contribution from the end chain particle :math:`\bar{w}_{E|A}(r_{15})`.
Conversely, the local density of particles of type B surrounding site 1, :math:`\rho_{B|1}` will only include intermolecular contributions to the local density.
Whether or not the local density of A particles surrounding site 1, :math:`\rho_{A|1}`, will include the *self* intramolecular contribution :math:`\bar{w}_{A|A}(0)` is determined by whether *self yes* or *self no* is listed as an argument by the user in the :doc:`ldd pair_coeff <pair_ldd>` command.

The choice of local density definition used for a particular interaction between types, :math:`\beta|\alpha` is set by the user via the :doc:`ldd pair_coeff <pair_ldd>` command with the different available options for indicator functions.
Note, that the weighting function, :math:`\bar{w}_{\beta|\alpha}`, and local density potential, :math:`U_{\beta|\alpha}` are not symmetric with respect to :math:`\alpha` and :math:`\beta`.
For instance, if :math:`\alpha = v` corresponds to solvents and :math:`\beta = u` corresponds to solutes, then the package allows the user to specify distinct weighting functions, :math:`\bar{w}_{v|u} \neq \bar{w}_{u|v}` and distinct LD potentials :math:`U_{v|u} \neq U_{u|v}`, for treating solvent around solute and solute around solvent respectively.

-----------------------------------------------------------------------------------------------------------------------------------------

**Local Density/ Square Gradient Potentials**

The LDD package implements forces for two kinds of local density dependent potentials.

LD (Local Density) potentials are defined:

.. math::

   U_{LD}(\mathbf{R}) = \sum_{I} \sum_{\beta} U_{\beta|t_I}(\rho_{\beta|I}),

where the first sum is over all sites, :math:`I` , and the second sum is over all site types, :math:`\beta`.
The form of :math:`U_{\beta|\alpha}` is specified by the user in the :doc:`ldd pair_coeff <pair_ldd>` command using the *potential* keyword.
For instance, :math:`U_{\beta|\alpha}(\rho_{\beta|I})` can be defined as :doc:`zero <ldd_potential_noforce>` or as a :doc:`constant <ldd_potential_constant>`, :doc:`linear <ldd_potential_linear>`, :doc:`quadratic <ldd_potential_quadratic>` function of the local density.
Alternatively, :math:`U_{\beta|\alpha}(\rho_{\beta|I})` can be specified by a table of values using the :doc:`tabulated <ldd_potential_table>` keywords.
Additionally the *potential* arg :doc:`mdpd <ldd_potential_mdpd>` is a quadratic function in local density that can be combined with the dpd pair_style to simulate pairwise mdpd interactions in the LD framework.

SG (Square Gradient) potentials are defined:

.. math::

   U_{SG}(\mathbf{R}) = \sum_{I} \sum_{\beta} U_{\nabla; \beta|t_I} (\rho_{\beta|I}) \left | \nabla_{I} \rho_{\beta|I} \right|^2

where the first sum is over all sites, :math:`I`, the second sum is over all site types, :math:`\beta` and the gradient :math:`\nabla_{I} = \frac{\partial}{\partial\mathbf{R_I}}` is evaluated with respect to the position of particle :math:`I`.
The form of the coefficient, :math:`U_{\nabla; \beta|t_I}`, is specified by the user in the :doc:`ldd pair_coeff <pair_ldd>` command using the (optional) *gradient* keyword.
Many of the available functional forms for the LD potential are also available for the SG coefficients, see :doc:`pair_ldd <pair_ldd>` for details.

-----------------------------------------------------------------------------------------------------------------------------------------

For any given site :math:`I`, once the local densities :math:`\rho_{\beta|I}` and the gradients of the local density :math:`\nabla_{I} \rho_{\beta|I}` have been determined, the corresponding forces are pairwise decomposable with pair forces analogous to those stated in `Delyser 2021 <https://pubs.aip.org/aip/jcp/article/156/3/034106/2839866/Coarse-grained-models-for-local-density-gradients>`_ and `Delyser 2019 <https://pubs.aip.org/aip/jcp/article/151/22/224106/198468/Analysis-of-local-density-potentials>`_.

Because of this, the LDD package must be used with its own custom :doc:`atom_style ldd <atom_style>`, which is basically the *atomic* atom_style, with extra fields for all local densities and local density gradients.

This implementation allows interactions to be specified in a pairtypewise manner using the typical pair_coeff commands for each of :math:`2^{\text{n}_{\text{types}}}` interactions.
In particular, pair_coeff 1 2 specifies the interaction for atoms of type 1 surrounded by atoms of type 2, and is distinct from pair_coeff 2 1, which corresponds to the interaction for atoms of type 2 surrounded by atoms of type 1.
Note that the pair_coeff input (e.g. pair_coeff 1 2) uses a "surrounded by" convention. This describes the local densities for particles type 1 "surrounded by" type 2.
Conversely, often in published works and equations for local densities we tend to use a "given" convention to describe the local density.
In this convention the equivalent interaction is denoted as 2|1, or type 2 given a central particle of type 1.

-----------------------------------------------------------------------------------------------------------------------------------------

Input File Examples
""""""""""""""""""""

Example 1) A simple atomic input example using only tabulated LD potentials

.. code-block:: LAMMPS

   ## Init system/units
   units real
   dimension 3
   atom_style ldd 2  # The total number of particle types

   newton on
   timestep 1

   read_data my_mix.data # System initialization, must preceed pair style init or ldd will exit with a warning
   velocity all create 300.0 22345 dist gaussian

   ## pair_style ldd, must be used with atomstyle ldd

   pair_style ldd 6.5 # Longest cutoff of all LD interactions

   # pair_coeff x y indicator keyword r0 rc self arg potential keyword args
   pair_coeff 1 1 indicator dpd 0.0 6.5 self yes potential table/lin LD_table.1.1.dat # type 1 surrounded by type 1, 1|1

   pair_coeff 1 2 ignore # type 1 surrounded by type 2 , 2|1
   # n.b. As the above line shows, all 2^ntypes interactions must be specified even if its to note that no LD interaction should be used for this type, see keyword ignore for details in pair_style doc
   pair_coeff 2 1 indicator dpd 0.0 5.5 self no potential table/lin LD_table.2.1..dat # type 2 surrounded by type 1, 1|2

   pair_coeff 2 2 indicator dpd 0.0 5.5 self yes potential table/line LD_table.2.2.dat # type 2 surrounded by type 2, 2|2

   ## Run / Output
   run_style verlet
   neighbor 15.0 bin
   thermo 500

   fix 1 all nvt temp 300.0 300.0 100.0
   dump 1 all ldd 500 dump.txt #A lammps trajectory file with LD info

   run 10000

Note in the above example, only pair types 1|2, 1|1 and 2|2 will interact according to LD potentials.

Example 2) Molecular input example that layers tabulated pair, LD and SG interactions via :doc:`pair_style hybrid <pair_hybrid>`

.. code-block:: LAMMPS

   ## Initialize system /units

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

   pair_coeff * * table lammps_nb_ALL.table nb_All 15.0 # Pair interaction cutoff should be >= pair_coeff ldd cutoff

   #pair_coeff x y ldd indicator keyword r0 rc self arg potential keyword arg gradient keyword arg
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

Note in the above example LD cross interactions are both specified, but not necessarily symmetric.
Also in this example, cross interactions opt not to use the optional *gradient* keyword.
Conversely homo-interactions (e.g. 1|1 2|2) layer SG potentials on top of the already defined LD and pair interactions.

See ``examples/PACKAGES/ldd`` for full input files for a 4-site chain example (topology like A-B-C-D).

Note that in all input examples the ldd package is used in the following order: 1) Define the ldd atom style 2) Define the number of atom types in the simulation 3) Call the ldd pair style 4) Call the coeffs for  all types of particles 5) define output 6) run.
If the simulation box is not defined with the same number of atom types that atom_style ldd is defined with, the LDD package will throw an error and close the program.


See :doc:`pair_style ldd <pair_ldd>` for all pair_coeff *args* options that exist, all restrictions, and more examples.

-----------------------------------------------------------------------------------------------------------------------------------------

read_data file format examples
"""""""""""""""""""""""""""""""

Atom style ldd is a basic atomic atom_style with per-atom fields added for local densities, gradients of local densities, LD energy contributions and SG energy contributions.
These can be reported using :doc:`dump ldd <dump_ldd>`, but this information is calculated just based on configurational data, and thus is not used for starting simulations.
read_data input therefore can follow usual atomic read_data input formats, and when hybridized with other atom styles, the .data file is the same as atomic hybridized with those styles.

Example 1) for when the :doc:`atom_style ldd <atom_style>` is used by itself (read_data is atomic)

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


Example 2) for the .data file up to the "Atoms" section of the read_data file when the atom_style ldd is hybridized with full. Bonds/Angles/Dihedral syntax are standard as listed in :doc:`read_data <read_data>`.

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


-----------------------------------------------------------------------------------------------------------------------------------------


**dump ldd output**

For each central particle :math:`I` and each particle type :math:`\beta`, there is a local density of :math:`\beta` particles that surround :math:`I` :math:`\rho_{\beta|I}`, and a corresponding gradient of that local density, :math:`\frac{\partial \rho_{\beta|I}}{\partial \boldsymbol{R}_I}`.
If for example :math:`t_I = \alpha`, and an LDD interaction has been defined for the :math:`\beta|\alpha` local densities, then the :doc:`dump ldd <dump_ldd>` command will report :math:`\rho_{\beta|I}` and :math:`\frac{\partial \rho_{\beta|I}}{\partial \boldsymbol{R}_I}` in addition to the simulation x, v, and f information.
:doc:`dump ldd <dump_ldd>` is essentially a custom lammps dump trajectory output with local density information.

See :doc:`dump ldd <dump_ldd>` for details.

This trajectory type is natively compatible with the `Bottom-up Open-source Coarse-graining Software <https://github.com/noid-group/BOCS>`_ which can be used to parameterize LD/SG potentials from atomistic data, as well as to convert these lammps trajectories to .trr files for analysis with `gromacs <https://www.gromacs.org/>`_ tools.

-----------------------------------------------------------------------------------------------------------------------------------------

**When to use LD/SG potentials:**

This package has been primarily developed in the context of parameterizing/simulating
bottom up coarse grained models.
We have found that LD/SG potentials are particularly useful for supplementing pair potentials when simulating CG models in the NPT ensemble
and in different interfacial environments.
In particular, such potentials dramatically improve CG descriptions of the pressure-density equation of state, interfacial profiles, liquid-vapor coexistence, and transferability between bulk and interfacial systems.
SG potentials are slower to simulate than LD potentials alone, so where applicable for coarse graining we recommend trying to construct/simulate LD potentials first and then [if necessary] adding SG potentials to refine further.
CG force fields with LD and SG potentials can be parameterized via the multi-scale coarse graining method implemented in `BOCSv5 and higher <https://github.com/noid-group/BOCS>`_.


A number of works have been published using this package (with `BOCS <https://github.com/noid-group/BOCS>`_) in its development stage for different CG studies and are listed below:

.. _DeLyser1:

**(DeLyser)** DeLyser, Noid. "Extending Pressure-Matching to Inhomogeneous Systems via Local-Density Potentials." The Journal of Chemical Physics, 147, no. 13: 134111 (2017)

.. _DeLyser2:

**(DeLyser)** DeLyser, Noid. "Analysis of local density potentials." The Journal of Chemical Physics, 151, no. 22:224106 (2019)

.. _DeLyser3:

**(DeLyser)** DeLyser, Noid "Bottom-up coarse-grained models for external fields and interfaces" The Journal of Chemical Physics 153, 224103 (2020)

.. _DeLyser4:

**(DeLyser)** DeLyser, Noid "Coarse-grained models for local density gradients" The Journal of Chemical Physics 156, 034106 (2021)

.. _Szukalo:

**(Szukalo)** Szukalo, Noid "A temperature-dependent length-scale for transferable local density potentials". The Journal of Chemical Physics, 159, 074104 (2023)

.. _Lesniewski1:

**(Lesniewski)** Lesniewski, Remsing, Noid "Coarse Graining the Hydrophobic Effect." The Journal of Chemical Physics, In Prep 2025

.. _Dutta:

**(Dutta)** Dutta, Lesniewski, Qaisrani, Noid, Andrienko, Nikoubashman. "Accurate coarse-graining of small organic molecules in melts and thin films using density-dependent potentials". In Prep 2025

.. _Lesniewski2:

**(Lesniewski)** Lesniewski, DeLyser, W. G. Noid "Progress toward a better BOCS: systematic coarse-graining with local density potentials" In Prep 2025
