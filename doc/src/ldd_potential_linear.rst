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

where :math:`m` and :math:`b` are the slope and intercept of the relevant LDD function.

Description
"""""""""""
This option can follow the *gradient* or *potential* keywords in the ldd pair_coeff command.
Following the *potential* keyword, it controls the functional form for :math:`U_{\rho}` and :math:`F_{\rho}  = -\frac{dU_{\rho}}{d\rho}`
Following the *gradient* keyword, it controls the function form for :math:`U_{\nabla}` and :math:`F_{\nabla}  = -\frac{dU_{\nabla}}{d\rho}`.
For generality we note each case with a dummy :math:`X` below.

Style *linear* applies the form:

.. math::
   U_{X}(\rho) &= m\rho + b \\
   F_{X}(\rho) &= -m , \text{ where } X=\rho \text{ or } X=\nabla. \\


Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`

------------

.. .. _DeLyser:

.. see :ref:`DeLyser` for more

**(DeLyser)** M.R. DeLyser, W.G. Noid, "Coarse-grained models for local density gradients." J. Chem. Phys., 156, 034106 (2021).

