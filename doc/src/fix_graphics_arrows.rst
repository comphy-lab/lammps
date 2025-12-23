.. index:: fix graphics/arrows

fix graphics/arrows command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID graphics/arrows Nevery mode keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* graphics/arrows = style name of this fix command
* Nevery = update graphics information every this many time steps
* mode = one of the following modes *dipole* or *force* or *velocity* or *variable* or *chunk*

  .. parsed-literal::

     *dipole* args = scale radius transparency
       scale = scale factor for the dipole moment to determine the arrow length
       radius = radius for arrows
       transparency = transparency setting for the arrows, a value from 0 (invisible) to 1 (fully opaque)
     *force* args = scale radius transparency
       scale = scale factor for the force vector to determine the arrow length
       radius = radius for arrows
       transparency = transparency setting for the arrows, a value from 0 (invisible) to 1 (fully opaque)
     *velocity* args = scale radius transparency
       scale = scale factor for the velocity vector to determine the arrow length
       radius = radius for arrows
       transparency = transparency setting for the arrow, a value from 0 (invisible) to 1 (fully opaque)
     *variable* args = xval yval zval radius transparency
       xval = xvalue for arrow vector
       yval = xvalue for arrow vector
       zval = xvalue for arrow vector
       radius = radius for arrows
       transparency = transparency setting for the arrow, a value from 0 (invisible) to 1 (fully opaque)
     *chunk* args = chunk-ID posval vecval scale radius transparency
       chunk-ID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command
       posval = ID of a per-chunk compute or fix that computes the position vector for the chunk
       vecval = ID of a per-chunk compute or fix that computes the arrow vector for the chunk
       yval = xvalue for arrow vector
       zval = xvalue for arrow vector
       radius = radius for arrows
       transparency = transparency setting for the arrow, a value from 0 (invisible) to 1 (fully opaque)

Examples
""""""""

.. code-block:: LAMMPS

   fix sf1 water graphics/arrows 200 dipole 100.0 0.1 1.0

Description
"""""""""""

.. versionadded:: TBD

This fix allows to add arrows to images rendered with :doc:`dump image
<dump_image>` using the *fix* keyword to represent vector properties
with arrows atoms from all atoms in the fix group.

The *group-ID* sets the group ID of the atoms selected to have the selected
property represented.  This may be a dynamic group.

The *Nevery* keyword determines how often the arrows graphics data is
updated.  This should be the same value as the corresponding *N*
parameter of the :doc:`dump <dump>` image command.  LAMMPS will stop
with an error message if the settings for this fix and the dump command
are not compatible.

There are two keywords available that determine what is shown: *display*
and *average*.  With *display* all atoms from all arrows and are in the
fix group will be displayed.  With *average* only the average position
of the atoms with the same atom-ID across all arrows will be shown.

The *type* quantity determines the color of the objects.  Its represents
an *atom* type and the atoms will be colored the same as the
corresponding atom type when the *type* coloring scheme is used in the
:doc:`dump image fix <dump_image>` command.  If the value of *type* is
0 then the atom type of the individual atoms is used.

The *radius* quantity determines the radius of the atoms.  The value of
*radius* is 0 then largest distance of an atom to the average position
from all arrowss is used.

The *transparency* quantity determines the transparency of the objects.
Its value must be between 0 (invisible) and 1 (fully opaque).

-----------

Dump image info
"""""""""""""""

.. versionadded:: TBD

Fix graphics/arrows is designed to be used with the *fix* keyword of
:doc:`dump image <dump_image>`.  The fix will add spheres based on the
atoms in the fix group across all arrows to *dump image* so that they
are included in the rendered image.

The *fflag1* setting of *dump image fix* are currently ignored.

and *fflag2* setting of *dump image fix* is used as an adjustment
to the radius of the rendered sphere.  Since the radius is already
determined by this fix, it is recommended to set this flag to 0.0.

Restart, fix_modify, output, run start/stop, minimize info
==========================================================

No information about this fix is written to :doc:`binary restart files
<restart>`.

None of the :doc:`fix_modify <fix_modify>` options apply to this fix.

Restrictions
""""""""""""

None

Related commands
""""""""""""""""

:doc:`fix_graphics <fix_graphics>`

Default
"""""""

none
