.. index:: pair_coeff ldd potential quadratic

ldd potential quadratic command
===============================

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
   pair_coeff 1 1 indicator lucy 0.0 7.2 self yes potential quadratic 2.0 -0.5 1.0

.. parsed-literal::

    keyword = *potential* 
      *potental* value = quadratic args
        *quadratic* args = a b c

where :math:`a`, :math:`b`, and  :math:`c` are the three coefficients in a quadratic polynomial.

Description
"""""""""""

Style *quadratic* applies a quadratic potential: 

.. math::
   u_{\rho}(\rho) &= a\rho^{2} + b\rho + c \\
   f_{\rho}(\rho) &= -2a\rho - b

Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`
   

