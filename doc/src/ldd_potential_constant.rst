.. index:: pair_coeff ldd potential constant

ldd potential/gradient constant command
=====================---------==========

Syntax
""""""

.. code-block:: LAMMPS

   atom_style ldd ntypes
   pair_style ldd rc
   pair_coeff i j keyword value ...

Examples
""""""""

.. code-block:: LAMMPS

   pair_style ldd 7.2
   pair_coeff 1 1 indicator lucy 0.0 7.2 self yes potential table/spline ldtable.1.1.dat gradient constant 0.01

.. parsed-literal::

    keyword = *gradient* 
      *gradient* value = constant args
        *constant* args = c 

where :math:`c` is the value of the potential at all values of the local density

Description
"""""""""""
This option can follow the *gradient* or *potential* keywords in the ldd pair_coeff command.
Following the *potential* keyword, it controls the functional form for :math:`U_{\rho}` and :math:`F_{\rho}`
Following the *gradient* keyword, it controls the function form for :math:`U_{\nabla}` and :math:`F_{\nabla}`.
For generality we note each case with a dummy :math:`X` below.


Style *constant* applies a constant value as the potential. 
Note that using the constant type for the LD potential won't actually change the behavior of the system. 
However, using the constant type for the gradient potential will change the behavior of the system.

.. math::
   U_{X}(\rho) &= c \\
   F_{X}(\rho) &= 0


Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`
   

