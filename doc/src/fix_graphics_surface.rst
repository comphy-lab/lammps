.. index:: fix graphics/surface

fix graphics/surface command
============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID graphics/surface Nevery isovalue radius keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* graphics/surface = style name of this fix command
* Nevery = update graphics information every this many time steps
* isovalue = isovalue for the particle property surface selection
* radius = radius describing the spread of the atoms to the density grid (distance units)
* one or more keyword/args pairs may be appended
* keyword = *quality* or *property* or *filename* or *binary* or *pad*

  .. parsed-literal::

     *quality* keyword = isosurface grid resolurion setting
        keyword = one of *min*, *low*, *med*, *high*, or *max*
     *property* value = per-atom property used to create the isosurface grid
        value = *none*, *mass*,  c_ID, c_ID[i], f_ID, f_ID[i], v_name
           *none* = 1.0 for all atoms
           *mass* = mass of the atoms
           c_ID = per-atom vector calculated by a compute with ID
           c_ID[I] = Ith column of per-atom array calculated by a compute with ID
           f_ID = per-atom vector calculated by a fix with ID
           f_ID[I] = Ith column of per-atom array calculated by a fix with ID
           v_name = per-atom vector calculated by an atom-style variable with name
     *filename* name = name pattern for output of a sequence of STL format mesh files (must contain a \* character to be replaced by the timestep number)
     *binary* logical = select whether to output a binary STL file (default is text mode)
     *pad* number = pad the timestep in the output file name with zeroes to have this many digits (default is 0)

Examples
""""""""

.. code-block:: LAMMPS

   fix sf1 water graphics/surface 200 0.1 2.5 quality high property mass
   fix stl water graphics/surface 200 0.01 1.5 filename water-surface-*.stl pad 5

Description
"""""""""""

.. versionadded:: TBD

This fix allows to add an isosurface graphics object representing the
triangulated surface at a given isovalue on a grid to images rendered
with :doc:`dump image <dump_image>` using the *fix* keyword and
optionally to output the computed mesh as a series of STL format files
for external processing.

The *group-ID* sets the group ID of the atoms selected to be represented
by the surface.  This may be a dynamic group.

The *Nevery* keyword determines how often the surface graphics data is
updated.  This should be the same value as the corresponding *N*
parameter of the :doc:`dump <dump>` image command.  LAMMPS will stop
with an error message if the settings for this fix and the dump command
are not compatible.

The surface objects will be colored by the atom type that is closest to
the isosurface grid cell when the *type* coloring scheme is used in the
:doc:`dump image fix <dump_image>` command is used.  The color may also
be that of the atom type's element with the *element* scheme, or just a
globally set constant color for the whole surface with the *const*
scheme.  Similarly, the transparency settings are inherited from the
corresponding atom types as well.  The globally set color and
transparency can be changed with the *fcolor* and *ftrans* keywords of
the :doc:`dump modify fcolor <dump_image>` command.

The *isovalue* argument sets the isovalue used to compute the surface.
The optimum value depends on the property that is being used and may
require some experimentation in combination with varying the *radius*
setting.  For the default setting of using the number density, a good
starting point should be a small positive value like 0.1.

The *radius* argument sets the width of the gaussian distribution
function used to distribute the per-particle data across the grid.  Its
value controls the smoothness of the isosurface and may need some
experimentation in combination with the choice of isovalue to achieve
the best output.

The *quality* keyword can have any of these words as argument: "min",
"low", "med", "high", or "max", and selects the grid resolution used
for the isosurface.  The actual grid dimensions depend on the geometry
of the simulation cell.

-----------

Dump image info
"""""""""""""""

.. versionadded:: TBD

Fix graphics/surface is designed to be used with the *fix* keyword of
:doc:`dump image <dump_image>`.  The fix will construct an isosurface
based on the positions and radii of the atoms in the fix group and pass
the graphics geometry information about it to *dump image* so that it is
included in the rendered image.

The *fflag1* setting of *dump image fix* determines whether the surface
will be rendered as a set of connected triangles (1) or as a mesh of
cylinders (2).

If using a mesh of cylinders, the *fflag2* setting determines the
diameter of the cylinders.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

None of the :doc:`fix_modify <fix_modify>` options apply to this fix.

Restrictions
""""""""""""

This fix is part of the EXTRA-FIX package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

Related commands
""""""""""""""""


Defaults
""""""""

quality = medium, property = none, binary = no, pad = 0, filename = none
