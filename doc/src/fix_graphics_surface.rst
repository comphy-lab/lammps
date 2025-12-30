.. index:: fix graphics/surface

fix graphics/surface command
============================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID graphics/surface Nevery type isovalue radius keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* graphics/surface = style name of this fix command
* Nevery = update graphics information every this many time steps
* type = an atom type value to select the color of the surface
* isovalue = isovalue for the surface selection (unitless)
* radius = radius of the atoms for computing the isosurface density (may be a variable)
* one or more keyword/args pairs may be appended
* keyword = *quality* or *filename* or *binary* or *pad*

  .. parsed-literal::

     *quality* keyword = surface algorithm quality setting
        keyword = one of *max*, *high*, *medium*, or *low*
     *filename* name = name pattern for output of a sequence of STL format mesh files (must contain a \* character to be replaced by the timestep number)
     *binary* logical = select whether to output a binary STL file (default is text mode)
     *pad* number = pad the timestep in the output file name with zeroes to have this many digits (default is 0)

Examples
""""""""

.. code-block:: LAMMPS

   fix sf1 water graphics/surface 200 5 isovalue 0.1 radius 1.5 quality high
   fix stl water graphics/surface 200 5 filename water-surface-*.stl pad 10

Description
"""""""""""

.. versionadded:: TBD

This fix allows to add an isosurface graphics object to images rendered with
:doc:`dump image <dump_image>` using the *fix* keyword and optionally to output
the computed mesh as a series of STL format files for external processing.

The *group-ID* sets the group ID of the atoms selected to be represented
by the surface.  This may be a dynamic group.

The *Nevery* keyword determines how often the surface graphics data is
updated.  This should be the same value as the corresponding *N*
parameter of the :doc:`dump <dump>` image command.  LAMMPS will stop
with an error message if the settings for this fix and the dump command
are not compatible.

The *type* quantity determines the color of the object.  Its represents
an *atom* type and the object will be colored the same as the
corresponding atom type when the *type* coloring scheme is used in the
:doc:`dump image fix <dump_image>` command is used.  The color may also
be that of the atom type's element or just a globally set constant color
for *all* objects of this fix instance, which can be changed using a
:doc:`dump modify fcolor <dump_image>` command.

The *isovalue* keyword argument sets the isovalue used to compute the surface.

*R* sets the radius for the atoms.  It may be a variable reference (*v_name*)
for a variable (called "name") that can be either an equal-style or an
atom-style or compatible variable.

The *quality* keyword can have any of these words as argument: "max", "high",
"medium", or "low" and selects the smoothness and resolution of the surface
graphics object.

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

In case of using triangles, the *fflag2* setting determines the
transparency of the triangles and must use a value between 0.0
(invisible) and 1.0 (fully opaque).  If using a mesh of cylinders, the
*fflag2* setting determines the diameter of the cylinders.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

None of the :doc:`fix_modify <fix_modify>` options apply to this fix.

Restrictions
""""""""""""

none

Related commands
""""""""""""""""


Defaults
""""""""

quality = medium, binary = no, pad = 0
