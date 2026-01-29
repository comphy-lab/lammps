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
       val = yes or no specifies whether to include the self term in defining the local density.
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

  ######### Example 1 : 1 type
  #atom_style ldd ntypes
  atom_style ldd 1

  #pair_style ldd maxcut
  pair_style ldd 7.2

  #pair_coeff x surrounded by y indicator type r0 rC self args potential type args
  pair_coeff 1 1 indicator lucy 0.0 7.2 self yes potential table_spline ldtable.1.1.dat # 1|1

  #########  Example 2 : 2 types

  #atom_style ldd ntypes
  atom_style ldd 2
  #pair_style ldd maxcut
  pair_style ldd 7.2

  #pair_coeff x surrounded by y indicator type r0 rC self args potential type args
  pair_coeff 1 1 indicator dpd 0.0 1.0 self no potential mdpd 25.0 # 1|1
  pair_coeff 1 2 indicator dpd 0.0 1.0 self no potential table_spline ldtable.1.2.dat gradient table_spline gradtable.1.2.dat # 2|1
  pair_coeff 2 1 indicator lucy 0.0 7.2 self no potential linear 1 1 # 1|2
  pair_coeff 2 2 ignore # 2|2

Description
"""""""""""

Style *ldd* implements the local density potential as first described by
Pagonabarraga and Frenkel :ref:`(Pagonabarraga)<Pagonabarraga>` and
additionally the square gradient of local densities first introduced by
:ref:`(DeLyser)<DeLyser>`.  The pair_style *ldd* is compatible with a
variety of molecular and atomic topologies, (See :doc:`Howto_ldd
<Howto_ldd>`) for details) and offers a variety of options for how to
define the local density.

Here for notational simplicity we outline the theory for style *ldd*
potentials with just 1 type of particle. :doc:`Howto_ldd <Howto_ldd>`
explains the more general :math:`n_{\text{type}}` case.

Consider an indicator function of pair distance, :math:`w(r)`, where
:math:`w(0)=1`, :math:`w(r_{c})=0`, and :math:`w(r)` is
continuous/differentiable for :math:`0 \leq r \leq r_{c}`.  We define
:math:`[w]` as the spatial integral of :math:`w(r)` and the normalized
indicator function as :math:`\bar{w}(r) = w(r)/[w]`.  The local density
around a given particle is then defined as :math:`\rho_{I} = \sum_{J}
\bar{w}(r_{IJ})`, where the sum over :math:`J` can either include or
exclude :math:`I`, depending on whether the argument following the
*self* keyword is yes or no.  Practically, the only difference is a
horizontal shift in the argument of the potential by :math:`1/[w]`. The
local density potential and corresponding pair forces are given by:

.. math::
   U_{LD}(\mathbf{R}) &= \sum_{I} U_{\rho}(\rho_{I}) \\
   F_{IJ;\rho}(\mathbf{R}) &= (F_{\rho}(\rho_{I}) + F_{\rho}(\rho_{J}))\bar{w}'(r_{IJ}) \mathbf{e}_{IJ}

where :math:`F_{\rho}(x) = -dU_{\rho}(x)/dx`.

The optional *gradient* keyword implements a potential that is a
function of the square gradient of the local density as described in
:ref:`(DeLyser)<DeLyser>`.  This keyword adds an extra square gradient
(SG) term to the overall potential:

.. math::
   U_{SG}(\mathbf{R}) = \sum_{I} U_{\nabla}(\rho_{I}) \| \nabla_{I} \rho_{I} \|^{2}

The (req.) *potential* keyword defines the form for
:math:`U_{\rho}`. See each :ref:`ldd_potential <ldd_potential>` doc page
for details.  The (opt.) *gradient* keyword defines the form for
:math:`U_{\nabla}`. See each :ref:`ldd_potential <ldd_potential>` doc
page for details.  The BOCS package offers a library of functional forms
:math:`g(\rho)` that may be implemented to define either
:math:`U_{\rho}` or :math:`U_{\nabla}` for LD and SG potentials.  These
functional forms are specified by the *args* that follow each respective
*keyword*.

The (req.) *indicator* keyword defines the form for :math:`w(r)`. See
each :ref:`ldd_indicator <ldd_indicator>` doc page for details.

The *ignore* keyword is used in simulations with multiple particle types
where only some of the type pairs have local density potentials acting
between them, as in the example above.

The *self* argument indicates whether the particle :math:`I=J` term is
included in the local densities and gradients calculated.  Note that if
site :math:`I` is of type :math:`\alpha`, then it cannot contribute to
the local density of :math:`\beta` sites around :math:`I` when
:math:`\beta \neq \alpha`.  Thus in this case the self term is
automatically excluded and a warning will be generated to note that a
requested self term in this case has been turned off.

.. _ldd_indicator:
.. toctree::
   :maxdepth: 1
   :caption: indicator options

   lucy <ldd_indicator_lucy>
   dpd <ldd_indicator_dpd>
   sphere <ldd_indicator_sphere>
   smooth <ldd_indicator_smooth>
   shell <ldd_indicator_shell>


.. _ldd_potential:
.. toctree::
   :maxdepth: 1
   :caption: potential/gradient keywords for defining functions of the local density

   noforce <ldd_potential_noforce>
   constant <ldd_potential_constant>
   linear <ldd_potential_linear>
   quadratic <ldd_potential_quadratic>
   tabulated <ldd_potential_table>
   mdpd <ldd_potential_mdpd>

Mixing, shift, table, tail correction, restart, rRESPA info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This pair style does not support automatic mixing.

This pair style does not support the :doc:`pair_modify <pair_modify>`
shift, table, and tail options.

This pair style does not write its information to :doc:`binary restart
files <restart>`.  Therefore, you must re-specify the pair_style and
pair_coeff commands in an input script that reads a restart file.


Restrictions
""""""""""""

This pair_style is part of the BOCS package. It is only enabled if
LAMMPS was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This pair style must be used with the :doc:`atom_style ldd <atom_style>`
or :doc:`atom_style hybrid <atom_style>` with ldd listed as an arg.
This atom style requires an argument of ntypes, which is the number of
particle types used in the simulation.

To save the properties associated with the local density, use :doc:`dump
style ldd <dump_ldd>`.

The *indicator*, *self*, and *potential* keywords are mandatory, unless
the *ignore* keyword is provided. The *gradient* keyword is optional.

For all :math:`2^{n_{\text{types}}}` of local densities that are
possible, the user must specify what kind of/if a local density
interaction should be defined.  The *ignore* keyword is used for turning
local density interactions off and not calculating local densities or
gradients for that kind of interaction.  Once the *ignore* keyword is
specified, the pair will be ignored regardless of future pair_coeff
commands, so use with care.  Conversely the *noforce* keyword will set
up a potential with a constant 0 force.  This option turns local density
interactions off, but allows local densities and gradients to be
calculated during the simulation without changing the statistics of the
simulation.  See :doc:`noforce <ldd_potential_noforce>` for details.

Note, that not all of the available potential styles should follow
*gradient*.  For example, the potential *noforce* arg is used to
calculate the local densities/square gradients of particles in a
simulation without actually applying a force.  Since the gradient
keyword is optional, and gradients are calculated at the same time as
local densities, you should just omit the gradient keyword instead of
specifying gradient noforce if you want LD/SG stats without changing the
forces in the simulation.

Related commands
""""""""""""""""

:doc:`atom_style ldd <atom_style>`, :doc:`dump ldd <dump_ldd>`, :doc:`Howto_ldd <Howto_ldd>`

----------

.. _Pagonabarraga:

**(Pagonabarraga)** I. Pagonabarraga, D. Frenkel. "Dissipative particle dynamics for interacting systems.", J. Chem. Phys., 115, 5015 (2001).

.. _DeLyser:

**(DeLyser)** M.R. DeLyser, W.G. Noid. "Coarse-grained models for local density gradients." J. Chem. Phys., 156, 034106 (2021).
