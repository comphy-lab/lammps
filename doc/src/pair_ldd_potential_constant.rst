.. index:: pair_coeff ldd potential constant

ldd potential constant command
==============================

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

Style *constant* applies a constant value as the potential. Using the constant type for the LD potential won't actually change the behavior of the system. However, using the constant type for the gradient potential will change the behavior of the system.

.. math::
   u_{\rho}(\rho) &= c \\
   f_{\rho}(\rho) &= 0


Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`
   

