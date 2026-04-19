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
* keyword = *radius* or *shading* or *alpha*

  .. parsed-literal::

     *radius* value = override per-atom or per-type radius if > 0.0 (distance units, default 0.0)
     *alpha* value  = override multiplier for alpha shape extraction (distance units)
     *shading* value = *smooth* or *flat*
        *smooth* = compute per-vertex normals for smooth shading (default)
        *flat* = use face normals for flat shading

Examples
""""""""

.. code-block:: LAMMPS

   compute cc1 all chunk/atom molecule
   fix hull all graphics/chunk 100 cc1
   fix hull all graphics/chunk 100 cc1 radius 1.0 shading smooth
   fix hull all graphics/chunk 100 cc1 shading flat alpha 10.0

Description
"""""""""""

.. versionadded:: TBD

This fix generates graphics objects from chunks of atoms defined by the
:doc:`compute chunk/atom <compute_chunk_atom>` command.  Each chunk is
represented by a point cloud created by replacing each atom position
with the corners of an octahedron scaled to the radius of the atom.  A
triangulated surface is created from that point cloud using a 3-D
`Delaunay triangulation
<https://en.wikipedia.org/wiki/Delaunay_triangulation>`_ combined with
`alpha shape <https://en.wikipedia.org/wiki/Alpha_shape>`_ extraction.
This allows the resulting surface to follow concave features of the
molecular geometry rather than always producing a convex hull.  The
resulting list of graphics objects is passed to :doc:`dump image
<dump_image>` for rendering via the *fix* keyword.

The positions used for the generation of the graphics are based on
unwrapped coordinates which are then mapped back into the simulation
cell based on the position of the first atom in the cluster.  When a
cluster will straddle a periodic boundary it should be drawn only on one
side of the boundary.

The *group-ID* selects the atoms included in the hull computation.  Only
atoms that belong to the specified group **and** are assigned to a chunk
are considered.

The *Nevery* keyword determines how often the list of the graphics
objects is recomputed.  It should match the dump frequency of the
corresponding :doc:`dump image <dump_image>` command.

The color of the graphics objects depends on the coloring scheme
selected in :doc:`dump image <dump_image>` command.  With the *type* or
*element* coloring scheme the color is based on atom type as described
below, with the *const* coloring scheme a uniform color is used instead.
This color can be set with the *fcolor* keyword of the :doc:`dump modify
<dump_image>` command.  When using atom type based colors the vertices
of the surface are colored using the atom type of the closest atom
and the color between vertices is interpolated.

The optional *radius* keyword allows to override the radius value used
to determine the size of the represented graphics by scaling the
octahedron positions that represents each atom for computing the
surface.  If available, the per-atom radius (e.g. for simulations using
:doc:`atom style sphere <atom_style>`) is used, otherwise half of the
value of the Lennard-Jones *sigma* parameter for the atom type is used.

The optional *alpha* keyword allows to adjust the alpha shape extraction
algorithm which determines how closely the generated triangulation
follows the shape of chunks of atoms.  It should be at least about 3x
the average distance of closest neighbors.  For larger values, the
generated shape will become mostly a conventional convex hull. A value
of 0.0 (the default) triggers an estimation of a suitable value.

The optional *shading* keyword selects how triangle normals are
determined for rendering surfaces.  The *smooth* setting (the default)
computes averaged per-vertex normals so that adjacent triangles appear
curved and blend smoothly (except for sharp edges).  The *flat* uses
the face normal for all three corners of each triangle, giving the
surface a faceted appearance.

----------

Dump image info
"""""""""""""""

Fix graphics/chunk is designed to be used with the *fix* keyword of
:doc:`dump image <dump_image>`.  The fix constructs a list of graphics
objects based on the size and geometry of the chunks in the fix group
and passes the information to the image renderer.

The *fflag1* setting of *dump image fix* determines whether the surface is
rendered as connected rounded triangles (1) or as a wireframe mesh of
cylinders (2).  If using a wireframe mesh, the *fflag2* setting
determines the diameter of the cylinders.

----------

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

When running in parallel, the chunks and corresponding graphics objects
are currently computed separately for each subdomain, so that the graphics
will be different for all chunks that are distributed across sub-domains
depending on the number of processors used.

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

radius = 0.0, alpha = 0.0, shading = smooth
