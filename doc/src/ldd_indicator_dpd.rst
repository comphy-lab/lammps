.. index:: pair_coeff ldd indicator dpd

ldd indicator dpd command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   atom_style ldd ntypes
   pair_style ldd rc
   pair_coeff i j keyword values ...

Examples
""""""""

.. code-block:: LAMMPS

   pair_style ldd 1.0
   pair_coeff 1 1 indicator dpd 0.0 1.0 self yes potential table_spline ldtable.1.1.dat

Description
"""""""""""

The *dpd* indicator style employs a simple quadratic indicator function:

.. math::
   w(r) &= (1-r/r_{c})^{2} \\
   [w] &= 2 \pi r_{c}^{3} / 15 \\
   \bar{w}(r) &= 15 (1-r/r_{c})^{2} \theta(r-r_{c}) /(2 \pi r_{c}^{3})


Where :math:`\theta(x)` is the Heaviside function, with :math:`\theta(x)=1` if :math:`x \geq 0` and :math:`\theta(x)=0` otherwise.
Following *dpd*, you must provide values for :math:`r_{0}` and :math:`r_{c}`. Note, you must set :math:`r_{0}=0.0` for the *dpd* indicator type.

Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`

**(Groot & Warren)** R. D. Groot, Patrick B. Warren "Dissipative particle dynamics: Bridging the gap between atomistic and mesoscopic simulation.", J Chem Phys, 107, 11 (1997).

