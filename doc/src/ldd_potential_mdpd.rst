.. index:: pair_coeff ldd potential mdpd

ldd potential mdpd command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   atom_style ldd ntypes
   pair_style ldd rc
   pair_coeff i j keyword value ...

Examples
""""""""

.. code-block:: LAMMPS

   pair_style hybrid/overlay dpd 1.0 1.0 16802 ldd 0.75
   pair_coeff 1 1 dpd -40.0 4.5 1.0
   pair_coeff 1 1 ldd indicator dpd 0.0 1.0 self no potential mdpd 25.0

.. parsed-literal::

    keyword = *potential*
      *potental* value = mdpd args
        *mdpd* args = B

Description
"""""""""""
This option can follow the *gradient* or *potential* keywords in the ldd pair_coeff command.
Following the *potential* keyword, it controls the functional form for :math:`U_{\rho}` and :math:`F_{\rho}`
Following the *gradient* keyword, it controls the function form for :math:`U_{\nabla}` and :math:`F_{\nabla}`.
For generality we note each case with a dummy :math:`X` below.

.. math::
   U_{X}(\rho) = B \pi r_{c}^{4} \rho^{2} / 30 \\
   F_{X}(\rho) = - B \pi r_{c}^{4} \rho / 15 \text{ where } X=\rho \text{ or } X=\nabla.\\

Note that when used after the *potential* keyword, style *mdpd* applies a quadratic local density potential used in many MDPD studies, e.g. :ref:`(Warren)<Warren>` or :ref:`(Ghoufi)<Ghoufi>`.
It is traditionally used to supplement the standard dpd pair force.
As described by :ref:`(Ghoufi)<Ghoufi>`, the desired conservative force is given by:

.. math::
   F^{C}_{ij} & = A \omega_{C}(r_{ij}) \mathbf{e}_{IJ} + B(\rho_{i} + \rho_{j}) \omega_{d}(r_{ij}) \mathbf{e}_{ij} \\
   \omega_{d}(r) & = 1 - r/r_{d}

In practice, to obtain the desired force with the LDD package, we use the :doc:`dpd indicator style <ldd_indicator_dpd>` combined with the mdpd *potential* style above to obtain:

.. math::
   w(r) &= (1-r/r_{c})^{2} \\
   [w] &= 2 \pi r_{c}^{3} / 15 \\
   \bar{w}(r) &= 15 (1-r/r_{c})^{2} / (2 \pi r_{c}^{3}) \\
   U(\rho) &= B \pi r_{c}^{4} \rho^{2} / 30 \\
   F(\rho) &= - B \pi r_{c}^{4} \rho / 15 \\
   \mathbf{F}_{ij} &= (f(\rho_{i}) + F(\rho_{j}))\bar{w}'(r_{ij}) \mathbf{e}_{ij} \\
   &= B(\rho_{i} + \rho_{j})(1-r/r_{c})\mathbf{e}_{ij}


Restrictions
""""""""""""

The mdpd potential style is only equivalent to the traditional MDPD force above when hybridized with the :doc:`pair_style dpd <pair_dpd>` and used with the :doc:`dpd indicator style <ldd_indicator_dpd>`.
The basic pair style required is part of the DPD-BASIC package, while the ldd potential/gradient style is part of the pair style ldd implemented in the LDD package.
They are both only enabled if LAMMPS was built with each package.
See the :doc:`Build package
<Build_package>` page for more info.



Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`

----------

.. _Warren:

**(Warren)** P. Warren. "Vapor-liquid coexistence in many-body dissipative particle dynamics.", Phys. Rev. E., 68, 066702 (2003).

.. _Ghoufi:

**(Ghoufi)** A. Ghoufi, P. Malfreyt. "Mesoscale modeling of the water liquid-vapor interface: A surface tension calculation.", Phys. Rev. E., 83, 051601 (2011).
