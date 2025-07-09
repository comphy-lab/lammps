.. index:: pair_coeff ldd potential table/lin
.. index:: pair_coeff ldd potential table/spline

ldd potential table command
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

   pair_coeff 1 1 indicator lucy 0.0 7.2 self yes potential table/lin ldtable.1.1.dat
   pair_coeff 1 1 indicator lucy 0.0 7.2 self yes potential table/spline ldtable.1.1.dat
   pair_coeff 2 1 indicator lucy 0.0 7.2 self no potential table/lin ldtable.2.1.dat gradient table/gradline sgtable.2.1.dat
   pair_coeff 1 2 indicator lucy 0.0 7.2 self no potential table/spline ldtable.1.2.dat gradient table/gradspline sgtable.1.2.dat

.. parsed-literal::

    keyword = *potential* 
      *potental* value = table/XXX args
        *table/XXX* args = filename

The only argument following *table/lin* or *table/spline* or *table/gradlin* or *table/gradspline* is the filename for the table.

Description
"""""""""""

Styles *table/lin* and *table/spline* read in a file containing three columns: :math:`\rho` :math:`u_{\rho}(\rho)` :math:`f_{\rho}(\rho)` where the columns are separated by whitespace and :math:`f(\rho) = -du(\rho)/d\rho`. Note, the table spacing must be uniform. If the simulation ever encounters a particle with a local density outside the domain of values provided in the table, the simulation will crash and produce an error. Accordingly, we advise providing values starting at a much lower/higher :math:`\rho` value than you expect to sample.

Styles *table/gradlin* and *table/gradspline* read in a file containing two columns :math:`\rho` :math:`u_{\nabla}(\rho)` where the columns are seperated by whitespace. 
As above, the table spacing must be uniform. If the simulation ever encounters a particle with a local density outside the domain of values provided in the table, the simulation will crash and produce an error. 
Accordingly, we advise providing values starting at a much lower/higher domain than you expect to sample. 

Style *table/lin* interpolates the potential and force between grid points linearly. Style *table/spline* constructs a cubic spline for interpolation between grid points.

Style *table/lin* interpolates between the :math:`u_{\nabla}` values linearly, and constructs :math:`f_{\nabla}(\rho) = -du_{\nabla}/d\rho` with the analogous delta extrapolation. 
Style *table/gradlin* constructs a cubic spline for interpolation between `u_{\nabla}` and the corresponding quadratically interpolated `f_{\nabla}(\rho) = -du_{\nabla}/d\rho`. 

Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`


