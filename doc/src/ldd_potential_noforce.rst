.. index:: pair_coeff ldd potential noforce

ldd potential noforce command
=============================

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
   pair_coeff 1 1 indicator lucy 0.0 7.2 self yes potential noforce


.. parsed-literal::

    keyword = *potential* 
      *potential* value = noforce args
        *noforce* args = none

Examples
""""""""

.. code-block:: LAMMPS

  pair_coeff 1 1 indicator dpd 0.0 1.0 self no potential noforce

Description
"""""""""""
This option can follow the *potential* keyword in the ldd pair_coeff command.
Following the *potential* keyword, it controls the functional form for :math:`U_{\rho}` and :math:`F_{\rho}`.

Operationally, *noforce* is meant to be used to calculate local densities and gradients of the local density without applying a potential to your system. 
This could alternatively be calculated by post-processing the simulation with tools like `BOCS <https://github.com/noid-group/BOCS>`_.

In the code we implement this by setting the functional form of :math:`U_{\rho}` and :math:`F_{\rho}` to zero.

.. math::
   U_{\rho}(\rho) &= 0 \\
   F_{\rho}(\rho) &= 0

Note that since gradients of the local density are calculated when the local density is calculated, it is not necessary to use this option after the *gradient* keyword is used.

Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`
   

