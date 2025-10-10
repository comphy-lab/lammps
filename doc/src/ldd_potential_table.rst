.. index:: pair_coeff ldd potential table/lin
.. index:: pair_coeff ldd potential table/spline

ldd potential/gradient table/X command
========================================

Syntax
""""""

.. code-block:: LAMMPS

   atom_style ldd ntypes
   pair_style ldd rc
   pair_coeff i j keyword value ...

Examples
""""""""

.. code-block:: LAMMPS

   pair_coeff 1 1 indicator lucy 0.0 7.2 self yes potential table/lin ldtable.1.1.dat
   pair_coeff 1 1 indicator lucy 0.0 7.2 self yes potential table/spline ldtable.1.1.dat
   pair_coeff 2 1 indicator lucy 0.0 7.2 self no potential table/lin ldtable.2.1.dat gradient table/gradline sgtable.2.1.dat
   pair_coeff 1 2 indicator lucy 0.0 7.2 self no potential table/spline ldtable.1.2.dat gradient table/gradspline sgtable.1.2.dat

.. parsed-literal::

    keyword = *potential* 
      *potental* value = table/XXX args
        *table/XXX* args = filename


The only argument following *table/lin* or *table/spline* or *table/gradlin* or *table/gradspline* is the filename for the table.
Example tables are used for *potential table/x* and *gradient table/gradx* in the ``examples/PACKAGES/ldd`` system.

Description
"""""""""""
This family of arguments can follow the *gradient* or *potential* keywords in the ldd pair_coeff command. 
Following the *potential* keyword, these options control the functional form for :math:`U_{\rho}` and :math:`F_{\rho}`
Following the *gradient* keyword, these options control the function form for :math:`U_{\nabla}` and :math:`F_{\nabla}`.

Args *table/lin* and *table/spline* read in a file containing three columns: :math:`\rho`, :math:`U`, :math:`F`. 
E.g. For a LD potential these three columns should be :math:`\rho` :math:`U_{\rho}(\rho)` :math:`F_{\rho}(\rho)` where the columns are separated by whitespace and :math:`F(\rho) = -dU(\rho)/d\rho`. 
Note, the grid domain for the table spacing must be uniform. 
If the simulation ever encounters a particle with a local density outside the domain of values provided in the table, the simulation will exit with an error. 
Accordingly, we advise providing values starting at a much lower/higher :math:`\rho` value than you expect to sample.

Args *table/gradlin* and *table/gradspline* read in a file containing two columns.
E.g. for a SG potential, these two columns should be :math:`\rho` :math:`U_{\nabla}(\rho)` where the columns are separated by whitespace. 
If a table formatted for *table/lin* or *table/spline* is passed under this command, only the first two columns will affect the simulation. 
As above, the table spacing must be uniform. 
If the simulation ever encounters a particle with a local density outside the domain of values provided in the table, the simulation will exit with an error. 
Accordingly, we advise providing values starting at a much lower/higher domain than you expect to sample. 

Arg *table/lin*  interpolates between entries in the the :math:`U` and :math:`F` columns linearly. 
Similarly arg *table/spline* constructs a cubic spline for interpolation between grid points in each column.

Arg *table/gradlin* interpolates between the :math:`U` column values linearly, and constructs a :math:`F(\rho) = -dU/d\rho` with the corresponding delta interpolation. 
Similarly, Arg *table/gradlin* constructs a cubic spline for interpolating between entries in the :math:`U` column provided, and then constructs a quadratic interpolation for the corresponding :math:`F(\rho) = -dU/d\rho`. 

Note that for SG potentials, the expression for the pair additive force contains both the coefficient :math:`U_{\nabla}(\rho)` and its derivative :math:`F_{\nabla}(\rho)`.
Thus using arg *table/lin* and *table/spline* after the *gradient* keyword will not calculate the forces consistently with (DeLyser 2021). 
Thus, while e.g. *table/lin* and *table/gradlin* can be used somewhat interchangeably to generate forces following the *potential* keyword, the *gradient* keyword should only be used with the *table/gradlin* and *table/gradspline* keywords. 

Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`


