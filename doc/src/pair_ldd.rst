.. index:: pair_style ldd

pair_style ldd command
======================

Syntax
""""""

.. code-block:: LAMMPS

   atom_style ldd ntypes
   pair_style ldd rc
   pair_coeff i j keyword value ...

* atom_style arg = ntypes (number of types of particles)
* pair_style arg = rc (cutoff of longest indicator function)
* one or more keyword/value pairs must be appended

.. parsed-literal::

    keyword = *indicator* or *self* or *potential* or *gradient* or *ignore*
      *indicator* values = type r0 rc
       type = dpd, lucy, shell, sphere, or smooth
       r0 = start of indicator function's decay to 0. must be 0.0 for dpd and lucy.
       rc = range of indicator function.
      *self* values = val
       val = yes or no based on if you want to include the self term in the LD calculation.
      *potential* values = type args
      *gradient* values = type args
       potential and gradient args are specific to the type, see its doc page
      *ignore* values = none

Examples
""""""""

.. code-block:: LAMMPS
   
  atom_style ldd 1
  pair_style ldd 7.2
  pair_coeff 1 1 indicator lucy 0.0 7.2 self yes potential table_spline ldtable.1.1.dat

  atom_style ldd 2
  pair_style ldd 1.0
  pair_coeff 1 1 indicator dpd 0.0 1.0 self no potential mdpd 25.0
  pair_coeff 1 2 indicator dpd 0.0 1.0 self no potential table_spline ldtable.1.2.dat gradient table_spline gradtable.1.2.dat
  pair_coeff 2 2 ignore

Description
"""""""""""

Style *ldd* implements the local density potential as first described by Pagonabarraga and Frenkel :ref:`(Pagonabarraga)<Pagonabarraga>`. Given an indicator function :math:`w(r)` such that :math:`w(0)=1`, :math:`w(r_{c})=0`, and :math:`w(r)` continuously and differentiably transitions from 1 to 0 across :math:`0 \leq r \leq r_{c}`, we define :math:`[w]` as the spatial integral of :math:`w(r)` and the normalized indicator function as :math:`\bar{w}(r) = w(r)/[w]`. We then define the local density around a given particle as :math:`\rho_{i} = \sum_{j} \bar{w}(r_{ij})`, where the sum over :math:`j` can either include or exclude :math:`i`, depending on whether the argument following the self keyword is yes or no. Practically, the only difference is a horizontal shift in the potential by :math:`1/[w]`. The local density potential and corresponding pair forces are given by:

.. math::
   U_{\rho}(\mathbf{r}) &= \sum_{i} u_{\rho}(\rho_{i}) \\
   F_{ij;\rho}(\mathbf{r}) &= (f_{\rho}(\rho_{i}) + f_{\rho}(\rho_{j}))\bar{w}'(r_{ij}) \mathbf{e}_{ij}
   
where :math:`f_{\rho}(x) = -du_{\rho}(x)/dx`.

The optional gradient keyword implements the gradient expansion of the local density potential as described in :ref:`(DeLyser)<DeLyser>`. These simulations add an extra term to the overall potential:

.. math::
   U_{\nabla}(\mathbf{r}) = \sum_{i} u_{\nabla}(\rho_{i}) (\nabla_{i} \rho_{i})^{2}



The *ignore* keyword is used in simulations with mutliple particle types where only some of the type pairs have local density potentials acting between them, as in the example above.

See individual files for the different forms of the available indicator functions and potentials.

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This pair style does not support automatic mixing. 

This pair style does not support the :doc:`pair_modify <pair_modify>` shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`.
Therefore, you must re-specify the pair_style and pair_coeff commands in an input script that reads a restart file.


Restrictions
""""""""""""

This pair style must be used with the atom_style ldd. This atom style requires an argument of ntypes, which is the number of particle types you will use in the simulation.

To save the properties associated with the local density, use dump style ldd.

The *indicator*, *self*, and *potential* keywords are mandatory, unless the *ignore* keyword is provided. The *gradient* keyword is optional.

Any value that can follow the *potential* keyword can also follow the *gradient* keyword, with the same arguments to that value. However, not all of the available styles should follow *potential*. For example, a constant ld potential will not change the simulation behavior. Simiarly, not all of the available styles should follow *gradient*. For example, the potential noforce style is used to calculate the local densities of particles in a simulation without actually applying a force. Since the gradient keyword is optional, you should just omit it instead of specifying gradient noforce.

Related commands
""""""""""""""""

:doc:`pair_ldd_indcator_dpd <pair_ldd_indicator_dpd>`, :doc:`pair_ldd_indicator_lucy <pair_ldd_indicator_lucy>`, 
:doc:`pair_ldd_indicator_shell <pair_ldd_indicator_shell>`, :doc:`pair_ldd_indicator_smooth <pair_ldd_indicator_smooth>`, 
:doc:`pair_ldd_indicator_sphere <pair_ldd_indicator_sphere>`, 
:doc:`pair_ldd_potential_noforce <pair_ldd_potential_noforce>`, :doc:`pair_ldd_potential_constant <pair_ldd_potential_constant>`, 
:doc:`pair_ldd_potential_linear <pair_ldd_potential_linear>`, :doc:`pair_ldd_potential_quadratic <pair_ldd_potential_quadratic>`, 
:doc:`pair_ldd_potential_mdpd <pair_ldd_potential_mdpd>`, :doc:`pair_ldd_potential_table <pair_ldd_potential_table>`


----------

.. _Pagonabarraga:

**(Pagonabarraga)** I. Pagonabarraga, D. Frenkel. "Dissipative particle dynamics for interacting systems.", J. Chem. Phys., 115, 5015 (2001).

.. _DeLyser:

**(DeLyser)** M.R. DeLyser, W.G. Noid. "Corase-grained models for local density gradients." J. Chem. Phys., ???, ??? (2022).
