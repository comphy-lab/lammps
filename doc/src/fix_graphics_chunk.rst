.. index:: fix graphics/chunk

fix graphics/chunk command
==========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID graphics/chunk Nevery chunkID keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* graphics/chunk = style name of this fix command
* Nevery = update graphics information every this many time steps
* chunkID = ID of :doc:`compute chunk/atom <compute_chunk_atom>` command
* zero or more keyword/args pairs may be appended
* keyword = *radius* or *shading*

  .. parsed-literal::

     *radius* value = radius for hull inflation (distance units, default 0.0)
     *shading* value = *smooth* or *flat*
        *smooth* = compute per-vertex normals for smooth shading (default)
        *flat* = use face normals for flat shading

Examples
""""""""

.. code-block:: LAMMPS

   compute cc1 all chunk/atom molecule
   fix hull all graphics/chunk 100 cc1
   fix hull all graphics/chunk 100 cc1 radius 1.0 shading smooth
   fix hull all graphics/chunk 100 cc1 radius 0.5 shading flat

Description
"""""""""""

.. versionadded:: TBD

This fix generates convex hull graphics objects from chunks of atoms
defined by the :doc:`compute chunk/atom <compute_chunk_atom>` command.
One convex hull is created for each chunk that contains atoms.  The
resulting triangle mesh is passed to :doc:`dump image <dump_image>` for
rendering via the *fix* keyword.

The *group-ID* selects the atoms included in the hull computation.  Only
atoms that belong to the specified group **and** are assigned to a chunk
are considered.

The *Nevery* keyword determines how often the hull geometry is
recomputed.  It should match the dump frequency of the corresponding
:doc:`dump image <dump_image>` command.

Special cases are handled automatically:

* Chunks with a single atom are rendered as a sphere.
* Chunks with exactly two atoms are rendered as a capped cylinder.
* Chunks with three or more atoms produce a full 3-D convex hull.

For the special cases of one or two atoms, a minimum visible radius of
0.1 (in distance units) is used if the effective radius would otherwise
be zero (i.e., when neither the *radius* keyword is set nor per-atom
radii exist).  This ensures that point particles are still visible.

The convex hull triangles are colored per vertex using the atom type of
the closest atom (from the point set that built the hull) when the
*type* or *element* coloring scheme is selected in :doc:`dump image
<dump_image>`.  With the *const* coloring scheme a uniform color is used
instead.  The color can be set with the *fcolor* keyword of the
:doc:`dump modify <dump_image>` command.

The optional *radius* keyword inflates the hull outward from the
centroid of each chunk by the specified distance.  If per-atom radii
exist (e.g. from :doc:`atom_style sphere <atom_style>`), the maximum
per-atom radius across all atoms in each chunk is used as the inflation
distance for that chunk.  The *radius* keyword value serves as the
default when no per-atom radii are available.  The default value of 0
corresponds to point particles with no inflation.

The optional *shading* keyword selects the normal computation used for
rendering.  *smooth* (the default) computes averaged per-vertex normals
so that adjacent triangles blend smoothly.  *flat* uses the face normal
for all three corners of each triangle, giving the hull a faceted
appearance.

----------

Dump image info
"""""""""""""""

Fix graphics/chunk is designed to be used with the *fix* keyword of
:doc:`dump image <dump_image>`.  The fix constructs one convex hull per
chunk and passes the geometry as triangles with per-vertex normals
(TRINORM objects) to the image renderer.

The *fflag1* setting of *dump image fix* determines whether the hull is
rendered as connected rounded triangles (1) or as a wireframe mesh of
cylinders (2).

If using a wireframe mesh, the *fflag2* setting determines the diameter
of the cylinders, otherwise it is ignored.

Restart, fix_modify, output, run start/stop, minimize info
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

No information about this fix is written to :doc:`binary restart files
<restart>`.

None of the :doc:`fix_modify <fix_modify>` options apply to this fix.

Restrictions
""""""""""""

This fix is part of the GRAPHICS package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

This fix is not compatible with 2d simulations.

Related commands
""""""""""""""""

:doc:`compute chunk/atom <compute_chunk_atom>`,
:doc:`fix graphics/arrows <fix_graphics_arrows>`,
:doc:`fix graphics/isosurface <fix_graphics_isosurface>`,
:doc:`fix graphics/labels <fix_graphics_labels>`,
:doc:`fix graphics/lines <fix_graphics_lines>`,
:doc:`fix graphics/objects <fix_graphics_objects>`,
:doc:`fix graphics/periodic <fix_graphics_periodic>`,
:doc:`dump image <dump_image>`

Defaults
""""""""

radius = 0.0, shading = smooth
