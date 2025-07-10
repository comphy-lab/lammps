.. index:: pair_coeff ldd indicator sphere

ldd indicator sphere command
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
   pair_coeff 1 1 indicator sphere 6.32 9.6 self yes potential table_spline ldtable.1.1.dat

Description
"""""""""""

The *sphere* indicator style employs a smoothed heaviside indicator function. If you consider a large sphere of radius :math:`X` and a small sphere of radius :math:`x`, with :math:`r_{0} = X - x` and :math:`r_{c} = X + x`, then :math:`w(r)` is the fraction of the small sphere contained in the large sphere when they are separated by a distance :math:`r`.

.. math::
   w(r) =
   \begin{cases}
   1 & r \leq r_{0} \\
   c_{3}r^{3} + c_{1} r + c_{0} + c_{-1}/r  &r_{0} \leq r \leq r_{c} \\
   0 & r \geq r_{c}
   \end{cases}

.. math::
   [w] &= 4 \pi ( r_{0}^{3} + \sum_{i}^{-1,0,1,3} c_{i} / (i+3) * (r_{c}^{i+3}-r_{0}^{i+3})) \\
   c_{-1} &= -3 (X^{4}+x^{4})/(16x^{3}) + 3X^{2}/(8x) \\
   c_{0} &= (X^{3}+x^{3})/(2x^{3}) \\
   c_{1} &= -3(X^{2}+x^{2})/(8x^{3}) \\
   c_{3} &= (16x^{3})^{-1} 

Following *sphere*, you must provide values for :math:`r_{0}` and :math:`r_{c}`

Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`



