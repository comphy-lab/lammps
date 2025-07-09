.. index:: pair_coeff ldd potential linear

ldd potential linear command
============================

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

Style *linear* applies a linear potential: 

.. math::
   u_{\rho}(\rho) &= m\rho + b \\
   f_{\rho}(\rho) &= -m


Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`
   

