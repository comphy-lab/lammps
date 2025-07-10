.. index:: pair_coeff ldd/indicator/lucy

ldd indicator lucy command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   atom_style ldd ntypes
   pair_style ldd rc
   pair_coeff i j keyword values ...

Examples
""""""""

.. code-block:: LAMMPS

   pair_style ldd 7.2
   pair_coeff 1 1 indicator lucy 0.0 7.2 self yes potential table_spline ldtable.1.1.dat

Description
"""""""""""

The *lucy* indicator style employs the Lucy function popularized from smoothed particle hydrodynamics as :math:`w(r)`:

.. math::
   w(r) &= (1-r/r_{c})^{3}(1+3r/r_{c}) \\
   [w] &= 16 \pi r_{c}^{3} / 105 \\
   \bar{w}(r) &= 105 (1-r/r_{c})^{3}(1+3r/r_{c})/(16 \pi r_{c}^{3})
 
Following the *lucy* argument after the indicator keyword, you must specify the values for :math:`r_{0}` and :math:`r_{c}`. For the *lucy* indicator type, you must set :math:`r_{0}=0.0`. 

Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`


