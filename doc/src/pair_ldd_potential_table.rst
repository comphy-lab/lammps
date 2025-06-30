.. index:: pair_coeff ldd potential table/lin
.. index:: pair_coeff ldd potential table/spline

ldd potential table command
===========================

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

.. parsed-literal::

    keyword = *potential* 
      *potental* value = table/XXX args
        *table/XXX* args = filename

The only argument following *table/lin* or *table/spline* is the filename for the table.

Description
"""""""""""

Styles *table/lin* and *table/spline* read in a file containing three columns: :math:`\rho` :math:`u(\rho)` :math:`f(\rho)` where the columns are separated by whitespace and :math:`f(\rho) = -du(\rho)/d\rho`. Note, the table spacing must be uniform. If the simulation ever encounters a particle with a local density outside the domain of values provided in the table, the simulation will crash and produce an error. Accordingly, we advise providing values starting at :math:`\rho = 0` and ending at a value much higher than you expect to sample.

Style *table/lin* interpolates the potential and force between grid points linearly. Style *table/spline* constructs a cubic spline for interpolation between grid points.
   

Related commands
""""""""""""""""

:doc:`pair_ldd <pair_ldd>`


