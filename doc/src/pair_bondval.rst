.. _pair_bondval:

pair_style bondval, bondval/vec (Bond-Valence Potential)
=========================================================

Syntax
------

The bond-valence potential is typically used with ``hybrid/overlay`` in order
to combine Coulomb, Lennard-Jones repulsion, bond-valence, and bond-valence
vector contributions.

.. code-block:: LAMMPS

   pair_style hybrid/overlay lj/cut/coul/long cutoff1 cutoff2 bondval cutoff1 cutoff2 bondval/vec cutoff1 cutoff2

Pair coefficients are specified per atom-type pair:

.. code-block:: LAMMPS

   pair_coeff i j lj/cut/coul/long cutoff1 Bij
   pair_coeff i j bondval r0,ij Cij Si V0,ij rcut
   pair_coeff i j bondval/vec r0,ij Cij Di W0,i rcut

where:

* ``i,j`` = atom type indices
* ``cutoff1`` = global cutoff for LJ(distance units)
* ``Bij`` = short range cutoff parameter (distance units)
* ``r0,ij`` = first of Brown's empirical parameters for bond valence (distance units)
* ``Cij`` = second of Brown's empirical parameters for bond valence
* ``Si`` = penalty for deviating from ideal bond valence (energy units)
* ``V0,ij`` = ideal bond valence for given atom type
* ``rcut`` = global cutoff for coulombic term of LJ (distance units)
* ``Di`` = penalty for deviating from ideal bond valence vector (energy units)
* ``W0,i`` = ideal bond valence vector for given atom type


Examples
--------

Example for a three-element perovskite BaTiO3 system:

.. code-block:: LAMMPS

   pair_style hybrid/overlay lj/cut/coul/long 8.0 8.0 bondval 2.0 8.0 bondval/vec 2.0 8.0

   read_data btolammps.data

   # lj/cut/coul/long parameters: cutoff1 Bij
   pair_coeff 1 1 lj/cut/coul/long 2.0 2.44805
   pair_coeff 1 2 lj/cut/coul/long 2.0 2.32592
   pair_coeff 1 3 lj/cut/coul/long 2.0 1.98792
   pair_coeff 2 2 lj/cut/coul/long 2.0 2.73825
   pair_coeff 2 3 lj/cut/coul/long 2.0 1.37741
   pair_coeff 3 3 lj/cut/coul/long 2.0 1.99269

   # bondval parameters: r0,ij Cij Si V0,ij rcut
   pair_coeff 1 1 bondval 0.0000000 5.0000000 0.5973900 2.0000000 8.0000000
   pair_coeff 1 2 bondval 0.0000000 5.0000000 0.0000000 0.0000000 8.0000000
   pair_coeff 1 3 bondval 2.2900000 8.9400000 0.0000000 0.0000000 8.0000000
   pair_coeff 2 2 bondval 0.0000000 5.0000000 0.1653300 4.0000000 8.0000000
   pair_coeff 2 3 bondval 1.7980000 5.2000000 0.0000000 0.0000000 8.0000000
   pair_coeff 3 3 bondval 0.0000000 5.0000000 0.9306300 2.0000000 8.0000000

   # bondval/vec parameters: r0,ij Cij Di W0,i rcut
   pair_coeff 1 1 bondval/vec 0.0000000 5.0000000 0.0842900 0.1156100 8.0000000
   pair_coeff 1 2 bondval/vec 0.0000000 5.0000000 0.0000000 0.0000000 8.0000000
   pair_coeff 1 3 bondval/vec 2.1430000 8.9400000 0.0000000 0.0000000 8.0000000
   pair_coeff 2 2 bondval/vec 0.0000000 5.0000000 0.8248400 0.3943700 8.0000000
   pair_coeff 2 3 bondval/vec 1.7980000 5.2000000 0.0000000 0.0000000 8.0000000
   pair_coeff 3 3 bondval/vec 0.0000000 5.0000000 0.2800600 0.3165100 8.0000000

Description
-----------

The Bond-Valence potential is an empirical potential based on the conservation
principles of the bond-valence (bondval) and bond-valence vector (bondval/vec)
that is fitted to DFT calculations for a given bulk semiconductor.

The bond-valence for a given pair of atoms (:math:`V_{ij}`) is a measure of the bonding
strength and can be calculated from the lengths of its bonds (:math:`r_{ij}`) by:

.. math::
   V_{ij} = \left( \frac{r_{0,ij}}{r_{ij}} \right)^{c_{ij}}

where :math:`r_{0,ij}` and :math:`c_{ij}` are Brown’s empirical parameters for bond-valence. The
bond-valence vector is defined as a vector lying along the
bond :math:`\vec{V}_{ij} = V_{ij} \hat{R}_{ij}`,
where :math:`\hat{R}_{ij}` is the unit vector
pointing from atom i to atom j. The total energy (:math:`E`) consists of the Coulombic
energy (:math:`E_c`), the short-range repulsive energy (:math:`E_r`),
the bond-valence energy (:math:`E_{BV}`), the bond-valence vector energy (:math:`E_{BVV}`),
and the angle potential (:math:`E_a`):

.. math::
   E = E_c + E_r + E_{BV} + E_{BVV} + E_a

.. math::
   E_c = \sum_{i<j} \frac{q_i q_j}{r_{ij}}

.. math::
   E_r = \sum_{i<j} \left(\frac{B_{ij}}{r_{ij}}\right)^{12}

.. math::
   E_{BV} = \sum_i S_i \left(V_i - V_{0,i}\right)^2

.. math::
   E_{BVV} =\sum_i D_i \left(W_i^2-W_{0,i}^2\right)^2

.. math::
   E_a = k \sum_i^{N_{\mathrm{oxygen}}} \left(\theta_i - 180^{\circ}\right)^2

Where :math:`V_i = \sum_{j \neq i} V_{ij}` is the bond-valence sum(BVS),
:math:`W_i = \sum_{j \neq i} V_{ij}` is the bond-valence vector sum(BVVS),
:math:`q_i` is the ionic charge, :math:`B_{ij}` is the short-range repulsion
parameter, :math:`S_i` and :math:`D_i` are scaling parameters with the unit of
energy, k is the spring constant, and :math:`\theta_i` is the O-O-O angle along
the common axis of two adjacent oxygen octahedra.

The values (:math:`r_{ij},V_i,W_i,\theta_i`) are calculated at each timestep from
the positions of the atoms in the simulation; all other terms need to be defined
by the user, shown in the Syntax part.

The values come from training the potential against DFT energies for a given system
type. There exist multiple examples of trained and tested parameters for bulk perovskite
oxide materials published and publicly available by the group of Andrew M. Rappe, but
the potential can theoretically be trained and applied to any bulk crystalline semiconductor
system. For reported values of the trained parameters, the LAMMPS metal units are used.
The atom charges and angle potential are also fitted parameters that are required to be
defined by the user. Charges can be set via read_data or set. If no angle potential is
defined, it is assumed to be zero. Otherwise, it can be set with the following command:

.. code-block:: LAMMPS

   Angle Coeffs
   1 15.0 180   # type of angle potential, k(eV/rad^2), degree


Restrictions
------------

The ``lj_cut_coul_long`` is a part of the KSPACE package. Since we only fitted the
repulsion part of the LJ potential to DFT, users need to manually set lj2 and lj4 to be
zero in ``lj_cut_coul_long``. For the code to work with the right physics, ``bondval`` ,
``bondval/vec``, ``lj_cut_coul_long`` must be all used together with ``hybrid/overlay``.
They cannot be used as standalone styles.

Related commands
----------------

* :doc:`pair_style hybrid/overlay <pair_hybrid>`
* :doc:`angle_style <angle_style>`
* :doc:`kspace_style <kspace_style>`

Default
-------

None


----------

**(Grinberg)** I. Grinberg, V. R. Cooper, and A. M. Rappe, *Nature*, **419**, 909 (2002).

**(Shin)** Y.-H. Shin, J.-Y. Son, B.-J. Lee, I. Grinberg, and A. M. Rappe, *J. Phys.: Condens. Matter*, **20**, 015224 (2008).

**(Shin)** Y.-H. Shin, V. R. Cooper, I. Grinberg, and A. M. Rappe, *Phys. Rev. B*, **71**, 054104 (2005).

**(Liu)** S. Liu, I. Grinberg, and A. M. Rappe, *J. Phys.: Condens. Matter*, **25**, 102202 (2013).

**(Liu)** S. Liu, I. Grinberg, H. Takenaka, and A. M. Rappe, *Phys. Rev. B*, **88**, 104102 (2013).

**(Brown)** I. D. Brown, *Chem. Rev.*, **109**, 6858 (2009).

**(Brown)** I. Brown and R. Shannon, *Acta Crystallogr. A*, **29**, 266 (1973).

**(Brown)** I. Brown and K. K. Wu, *Acta Crystallogr. B*, **32**, 1957 (1976).
