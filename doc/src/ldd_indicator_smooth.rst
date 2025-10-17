.. index:: pair_coeff ldd indicator smooth

ldd indicator smooth command
==============================

Syntax
""""""

.. code-block:: LAMMPS

   atom_style ldd ntypes
   pair_style ldd rc
   pair_coeff i j keyword values ...

Examples
""""""""

.. code-block:: LAMMPS

   pair_style ldd 9.6
   pair_coeff 1 1 indicator smooth 6.32 9.6 self yes potential table_spline ldtable.1.1.dat

Description
"""""""""""

The *smooth* indicator style employs a smoothed Heaviside function that is one for :math:`r \leq r_0`, zero for :math:`r \geq r_c`, and smoothly interpolates between zero and one in the interval :math:`r_0 \leq r \leq r_c`.
This function is similar to the indicator function introduced by Sanyal and Shell, but has a continuous second derivative.

.. math::
   w(r) =
   \begin{cases}
   1 & r \leq r_{0} \\
   \sum_{i=0}^{5} c_{i}r^{i} &r_{0} \leq r \leq r_{c} \\
   0 & r \geq r_{c}
   \end{cases}

.. math::
   [w] &= 4 \pi ( r_{0}^{3}/3 + \sum_{i=0}^{5} c_{i} / (i+3) * (r_{c}^{i+3}-r_{0}^{i+3})) \\
   d &= (r_{0}^{5} - r_{c}^{5})/120 - (r_{0}^{3} - r_{c}^{3})r_{0}r_{c}/24 + (r_{0}-r_{c})r_{0}^{2}r_{c}^{2}/12 \\
   c_{0} &= (-r_{c}^{5}/120 + r_{0} r_{c}^{4}/24 - r_{0}^{2}r_{c}^{3}/12)/d \\
   c_{1} &= (r_{0}^{2} r_{c}^{2}) / (4d) \\
   c_{2} &= -(r_{0}r_{c}(r_{0}+r_{c}))/(4d) \\
   c_{3} &= (r_{0}^{2} + 4 r_{0}r_{c} + r_{c}^{2})/(12d) \\
   c_{4} &= -(r_{0} + r_{c})/(8d) \\
   c_{5} &= 1/(20d)

Following *smooth*, you must specify values for :math:`r_{0}` and :math:`r_{c}`.


Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`


**(Sanyal & Shell)** Tanmoy Sanyal, M. Scott Shell "Coarse-grained models using local-density potentials optimized with the relative entropy: Application to implicit solvation.", J Chem Phys, 145, 034109 (2016).

**(DeLyser)** Michael R. DeLyser, W. G. Noid "Extending pressure-matching to inhomogeneous systems via local-density potentials" J Chem Phys, 147, 134111 (2017)

