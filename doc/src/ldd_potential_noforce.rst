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

Style *noforce* does not actually apply a local density potential to your system. It is be used to calculate the local densities during a simulation (which could alternatively be calculated by post-processing the simulation).

.. math::
   u_{\rho}(\rho) &= 0 \\
   f_{\rho}(\rho) &= 0


Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`
   

