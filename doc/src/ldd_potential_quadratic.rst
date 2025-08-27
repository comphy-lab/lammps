.. index:: pair_coeff ldd potential quadratic

ldd potential/gradient quadratic command
=========================================

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
This option can follow the *gradient* or *potential* keywords in the ldd pair_coeff command.
Following the *potential* keyword, it controls the functional form for :math:`U_{\rho}` and :math:`F_{\rho}`
Following the *gradient* keyword, it controls the function form for :math:`U_{\nabla}` and :math:`F_{\nabla}`.
For generality we note each case with a dummy :math:`X` below.

Style *quadratic* applies a quadratic potential: 

.. math::
   U_{X}(\rho) &= a\rho^{2} + b\rho + c \\
   F_{X}(\rho) &= -2a\rho - b

Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`
   
------------

.. _DeLyser:

**(DeLyser)** M.R. DeLyser, W.G. Noid, "Coarse-grained models for local density gradients." J. Chem. Phys., 156, 034106 (2021).

.. _Warren:

**(Warren)** P. B. Warren, "Vapor-liquid coexistence in many-body dissipative particle dynamics." Phys. Rev. E., 68, 066702 (2003). 
