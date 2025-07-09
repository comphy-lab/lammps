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
       type = quadratic, linear, constant, mdpd, table/spline, table/lin, table/gradlin, table/gradspline, noforce
       potential/gradient args are specific to the type, see each ldd_potential doc page
      *gradient* values = type args
       type = quadratic, linear, constant, table/spline, table/lin, table/gradlin, table/gradspline
       potential/gradient args are specific to the type, see each ldd_potential doc page
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


The (req.) *potential* keyword defines the form for :math:`U_{\rho}`. See each ldd_potential doc page for details.
The (opt.) *gradient* keyword defines the form for :math:`U{\nabla}`. See each ldd_potential doc page for details.
The (req.) *indicator* keyword defines the form for :math:`w(r)`. See each ldd_indicator doc page for details.

The *ignore* keyword is used in simulations with mutliple particle types where only some of the type pairs have local density potentials acting between them, as in the example above.

The *self* argument indicates whether the particle i=j term is included in the local densities and gradients calculated. 
Note that for heterogenous LD types [e.g. type A surrounded by B] the self term is automatically excluded. 
This is automatically enforced by the ldd package 
and a warning will be issued to the user to note that a requested self term in this case has been turned off.


.. toctree::
   :maxdepth: 1
   :caption: indicator options
   
   ldd_indicator_lucy
   ldd_indicator_dpd
   ldd_indicator_sphere
   ldd_indicator_smooth
   ldd_indicator_shell



.. toctree::
   :maxdepth: 1
   :caption: potential and gradient keywords

   ldd_potential_noforce
   ldd_potential_constant
   ldd_potential_linear
   ldd_potential_quadratic
   ldd_potential_table
   ldd_potential_mdpd

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This pair style does not support automatic mixing. 

This pair style does not support the :doc:`pair_modify <pair_modify>` shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart files <restart>`.
Therefore, you must re-specify the pair_style and pair_coeff commands in an input script that reads a restart file.


Restrictions
""""""""""""

This pair style must be used with the atom_style ldd or atom_style hybrid ldd etc. 
This atom style requires an argument of ntypes, which is the number of particle types you will use in the simulation.

To save the properties associated with the local density, use dump style ldd.

The *indicator*, *self*, and *potential* keywords are mandatory, unless the *ignore* keyword is provided. The *gradient* keyword is optional.

For all :math:`2^{n_{\text{types}}}` possible local density interactions that are initialized when atom_style ldd is used, the user must specify what kind of/if a local density interaction should be 
defined for each pair of types. 
The *ignore* keyword is used for not setting a local density interaction/not calculating local densities or gradients for that kind of interaction. 
Once the *ignore* keyword is specified, the pair will be ignored regardless of future pair_coeff commands, so use with care.
Conversely the *noforce* keyword will set up a potential with a constant 0 force, which will calculate local densities and gradients for that kind of interaction but not change statistics of the simulation. 

Note, not all of the available potential styles should follow *gradient*. 
For example, the potential noforce style is used to calculate the local densities/square gradients of particles in a simulation without actually applying a force. 
Since the gradient keyword is optional, you should just omit it instead of specifying gradient noforce.

Related commands
""""""""""""""""

:doc:`ldd_indicator_dpd <ldd_indicator_dpd>`, :doc:`ldd_indicator_lucy <ldd_indicator_lucy>`, 
:doc:`ldd_indicator_shell <ldd_indicator_shell>`, :doc:`ldd_indicator_smooth <ldd_indicator_smooth>`, 
:doc:`ldd_indicator_sphere <ldd_indicator_sphere>`, 
:doc:`ldd_potential_noforce <ldd_potential_noforce>`, :doc:`ldd_potential_constant <ldd_potential_constant>`, 
:doc:`ldd_potential_linear <ldd_potential_linear>`, :doc:`ldd_potential_quadratic <ldd_potential_quadratic>`, 
:doc:`ldd_potential_mdpd <ldd_potential_mdpd>`, :doc:`ldd_potential_table <ldd_potential_table>`


----------

.. _Pagonabarraga:

**(Pagonabarraga)** I. Pagonabarraga, D. Frenkel. "Dissipative particle dynamics for interacting systems.", J. Chem. Phys., 115, 5015 (2001).

.. _DeLyser:

**(DeLyser)** M.R. DeLyser, W.G. Noid. "Corase-grained models for local density gradients." J. Chem. Phys., 156, 034106 (2021).
