.. index:: ldd indicator shell

ldd indicator shell command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   pair_style ldd rc
   pair_coeff i j keyword values ...

Examples
""""""""

.. code-block:: LAMMPS

   pair_style ldd 9.6
   pair_coeff 1 1 indicator shell 6.32 9.6 self yes potential table_spline ldtable.1.1.dat

Description
"""""""""""

The *shell* indicator style employs a smoothed Heaviside indicator function originally parameterized in :ref:`(Sanyal)<Sanyal>`:

.. math::
   w(r) =
   \begin{cases}
   1 & r \leq r_{0} \\
   c_{0} + c_{2} r^{2} + c_{4} r^{4} + c_{6}r^{6} & r_{0} \leq r \leq r_{c} \\
   0 & r \geq r_{c}
   \end{cases}

.. math::
   [w] = 4 \pi ( r_{0}^{3}/3 &+ c_{0} / 3 * (r_{c}^{3}-r_{0}^{3}) \\
                           &+ c_{2} / 5 * (r_{c}^{5}-r_{0}^{5}) \\
                           &+ c_{4} / 7 * (r_{c}^{7}-r_{0}^{7}) \\
                           &+ c_{6} / 9 * (r_{c}^{9}-r_{0}^{9}) )

.. math:: 
   d &= (1 - r_{0}^{2}/r_{c}^{2})^{3} \\
   c_{0} &= (1 - 3 r_{0}^{2}/r_{c}^{2})/d \\
   c_{2} &= 6 r_{0}^{2} / r_{c}^{4} / d \\
   c_{4} &= -3(1+r_{0}^{2}/r_{c}^{2})/r_{c}^{4} / d \\
   c_{6} &= 2/r_{c}^{6}/d

Following *shell*, you must specify values for :math:`r_{0}` and :math:`r_{c}`.

Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`

----------

.. _Sanyal:

**(Sanyal)** T. Sanyal, M.S. Shell. "Coarse-grained models using local-density potentials optimized with the relative entropy: Application to implicit solvation.", J. Chem. Phys., 145, 034109 (2016).

