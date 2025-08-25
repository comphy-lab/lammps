.. index:: pair_coeff ldd potential linear

ldd potential/gradient linear command
======================================

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
   pair_coeff 1 1 indicator lucy 0.0 7.2 self yes potential linear 2.0 -0.5

.. parsed-literal::

    keyword = *potential* 
      *potental* value = linear args
        *linear* args = m b

where :math:`m` and :math:`b` are the slope and intercept of the LD potential.

Description
"""""""""""
This option can follow the *gradient* or *potential* keywords in the ldd pair_coeff command.
Following the *potential* keyword, it controls the functional form for :math:`U_{\rho}` and :math:`F_{\rho}`
Following the *gradient* keyword, it controls the function form for :math:`U_{\nabla}` and :math:`F_{\nabla}`.
For generality we note each case with a dummy :math:`X` below.

Style *linear* applies the form: 

.. math::
   U_{X}(\rho) &= m\rho + b \\
   F_{X}(\rho) &= -m


Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`
   

