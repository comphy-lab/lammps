Granular surfaces
=================

.. versionadded:: TBD

As explained on the :doc:`Howto granular <Howto_granular>` doc page,
granular systems are composed of spherical particles with a diameter,
as opposed to point particles.  This means they have an angular
velocity and torque can be imparted to them to cause them to rotate.

The :doc:`Howto granular <Howto_granular>` doc page lists various
atom, pair, fix, and compute styles useful for creating granular
models for systems of interacting particles.

This page explains how you can also define granular surfaces which are
a collection of triangles (3d systems) or line segments (2d systems),
which act as boundaries interacting with the particles.  Different
kinds of particle/surface interactions can be specified with similar
options as the :doc:`granular pair style <pair_granular>` command.

----------

Global versus local surfaces
""""""""""""""""""""""""""""

A key point to understand is that LAMMPS supports two forms of
granular surfaces.  You cannot use both in the same simulation.

The first is *global* which means that each processor stores a copy of
all the triangles/lines.  This is suitable when a modest number of
triangles/lines is needed.  They can be large triangles/lines, any of
which span a significant fraction of the simulation box size in one or
more dimensions.

The second is *local* which means that the collection of
triangles/lines is distributed across processors in the same manner
that particles are distributed.  Each processor is assigned to a
sub-domain of the simulation box and owns whichever triangles/lines
have their center point in the processor's sub-domain.  Similar to
particles, each processor may also own ghost copies of triangles/lines
whose finite size overlaps with the processor's sub-domain.  The total
number of triangles/lines in the system can now be very large.  For
effective distribution and minimal communication, all the
triangles/lines should be small, no more than a few particle diameters
in size.  If even one larger triangle or line is defined then the
neighbor list cutoff and communication cutoff will be set
correspondingly larger, which can slow down the simulation.  Note that
a large triangle are line can instead be defined as multiple smaller
trianges or lines without changing the topology of the collective
surface.

One of these two commands must be specified to use *global* or *local*
surfaces in your granular simulation:

* :doc:`fix surface/global <fix_surface_global>`
* :doc:`fix surface/local <fix_surface_local>`

The :doc:`fix surface/global <fix_surface_global>` command reads in
the global surfaces in one of two ways.  The first option is from a
molecule file(s) previously defined by the :doc:`molecule <molecule>`
command.  The file should define triangles or lines with header
keywords and a Triangles or Lines section.  The second option is from
a text or binary STL (stereolithography) file which defines a set of
triangles.  It can only be used with 3d simulations.

The :doc:`fix surface/local <fix_surface_local>` command defines local
surface in one of three ways.  The first two options are the same
molecule and STL files explained in the previous paragraph.  In this
case, the list of triangles/lines is distributed across processors
based on the center point of each triangle/line.  The third option is
to include them in a LAMMPS data file which has been previously read
in via the :doc:`read_data <read_data>` command.  If the file has a
Triangles or Lines section, then triangles/lines will be read in and
distributed along with any particles the data file includes, assuming
an appropriate :doc:`atom_style <atom_style>` has been specified, as
explained below.

----------

Surface attributes
""""""""""""""""""

For both global and local surfaces, each triangle/line is assigned a
*type* and a *molecule ID*.  This is done when surfaces are read-in
from a molecule, STL, or data file.  Since STL files do not define
types or molecule IDs, the :doc:`fix surface/global
<fix_surface_global>` and :doc:`fix surface/local <fix_surface_local>`
commands specify the type and molecule ID that will be assigned to the
read-in triangles.  The :doc:`fix surface/global <fix_surface_global>`
command also allows use of the :doc:`fix_modify type/region <fix_modify>`
command to assign types based on a geometric region.  Since local
surfaces are effectively particles, the :doc:`set <set>` command can
be used to alter the *type* or *molecule ID* of any triangle or line.

For both global and local surfaces, types are used to define the style
of granular interactions for individual triangles/lines.  Different
styles can be used within a single object consisting of connected
triangles/lines.  See the Surface Connectivity section below.

Molecule IDs are used to determine which triangles/lines are connected.
They are therefore intended to be assigned uniquely to each
inter-connected set of triangles/lines, as if each object were a
"molecule".

For local surfaces, the molecule ID can be used to define groups which
enables assignment of different motions to different surface objects.
See the Surface Motion section below.  Various other LAMMPS commands
operate on groups or molecules and can thus be used to gather
statistics about or output information about individual surface
objects.

----------

Atom styles for granular surfaces
"""""""""""""""""""""""""""""""""

For all three ways of defining *local* surfaces, the triangles/lines
are stored internally in LAMMPS as triangle-style or line-style
particles.  This means you must use a hybrid atom style which includes
one of these two sub-styles (for 3d or 2d):

* :doc:`atom_style tri <atom_style>` for 3d simulations
* :doc:`atom_style line <atom_style>` for 2d simulations

The atom_style hybrid command must also define a :doc:`atom_style
sphere <atom_style>` sub-style for the granular particles which
interact with the surfaces.

Note that for molecule or STL file input, the :doc:`fix surface/local
<fix_surface_local>` command reads the file(s) and uses the values for
each surface to create a single new triangle or line particle.  For
data file input, the triangle/line particles are created when the data
file is read.

For granular simulations with *global* surfaces, a hybrid atom style
which defines triangle-style or line-style particles should NOT be
used.  Typically only an :doc:`atom_style sphere <atom_style>` command
is needed to define the properties of particles in the simulation.
The triangles/lines are stored by the :doc:`fix surface/global
<fix_surface_global>` command and not as triangle-style or line-style
particles.

----------

Rules for surface topology
""""""""""""""""""""""""""

For both *global* and *local* surfaces, granular particles interact
with both sides of each triangle or line segment.  This means a surface
such as a mixer blade can be infinitely thin.

Triangles and line segments can be "connected" to form a contiguous
surface if they share common edges or corner point (triangles) or end
points (line segments).  Each triangle edge or corner point can be
shared by multiple adjacent triangles.  A triangle edge is shared by
two triangles if both the end points of the edge are corner points of
both triangles.  A triangle corner point is shared by two triangles if
it is a corner point of both.  Likewise a line segment end point is
shared by two line segments if it is an end point of both segments.

NOTE: say something about epsilon criterion for "shared" ?

There is no requirement that a triangle edge or triangle corner point
or line segment end point be connected to another triangle or line
segment.  If an edge or point has no connection, it is a free
(unconnected) edge or point.  Particles interact with the free edge
or corner point in a manner consistent with forces generated by
particles overlapping with the interior of a triangle or line segment.

NOTE: need an explanation of PBC

No check is made to see if two triangles or line segments intersect
each other; this is allowed if it makes sense for the geometry of the
collection of surfaces. However, intersections can cause issues if the
missing connectivity leads to inaccurate forces.

As an example of a valid intersection, consider a 2d simulation which
mixes a container of granular particles.  *Global* line segments are
used to define both the box-shaped container and the mixer in the
center.  The 4 mixer blades are in the shape of a large X and are made
to rotate using the :doc:`fix_modify <fix_modify>` command (see
below).

.. figure:: img/gransurf_mixer.jpg
            :figwidth: 33%
            :align: right
            :target: _images/gransurf_mixer.jpg

The 2 blades could be defined by 2 line segments which cross each
other at their centers (left).  Or the 2 blades could be defined by 4
line segments, all of which have a common endpoint at the center of the
mixer (middle).  Or the 2 blades could be defined by 4 non-touching line
segments, all of which have a distinct endpoint near the center of the
mixer, but displaced from it by a distance less than the radius of a
granular particle (right).  In any of these cases, when a particle gets
very close to the center of the mixer it will interact with both nearby
line segments as expected.

As an example of an invalid intersection, consider a 2d simulation of a
T shaped object defined by 2 line segments (next figure, left).  The
vertical line segment (blue) ends at the midpoint of the horizontal line
segment (red).  Without proper connectivity, there is no way to censor a
force from the geometrically hidden vertical segment as a particle (gray)
moves horizontally across the top of the T. In contrast, if the T shape is
defined by 3 line segments that all share a common endpoint at the center
of the top of the T (right), then the connectivity would censor the force
from the vertical segment (blue).

.. figure:: img/gransurf_T_geom.jpg
            :figwidth: 33%
            :align: right
            :target: _images/gransurf_T_geom.jpg

See the next section on connectivity for how two triangles or line
segments are treated if they share a common edge (triangle) or point
(triangle or line).

----------

Surface connectivity
""""""""""""""""""""

If multiple triangles/lines are used to define a contiguous surface
which is flat or gently curved or has sharp edges or corners, LAMMPS
will detect when two or more line segments (2d) in the same molecule
share the same endpoint.  Or when two or more triangles (3d) in the
same molecule share the same edge or same corner point.

This connectivity is stored internally and is used when appropriate to
calculate accurate forces on particles which simultaneously overlap
with 2 or more connected triangles or line segments.

Consider the simulation model of the previous section for a 2d mixer
now defined by *local* line segments.  The flat surface of each mixer
blade (and container box faces) is defined by multiple small line
segments.  It is important that these line segments be "connected" so
that when a particle contacts two adjacent line segments at the same
time, the resulting force on the particle is the same as it would be
if it were contacting the middle of a single long line segment.

Here is how to ensure that LAMMPS detects the appropriate connections.

For either *global* or *local* surfaces, if the triangles/lines are
defined in a molecule or STL file, then 3 corner points (triangle) or
2 end points (line) will be listed for each triangle/line in the file.
LAMMPS will only make a connection between 2 triangles or lines if a
shared point is EXACTLY the same in both.  This is a single point in
both for a corner point or end point connection.  It is two points in
both triangles for an edge connection.

For *local* surfaces, if the triangles/lines are defined in a data
file, then 3 corner points (triangle) or 2 end points (line) will be
listed for each triangle/line in the file.  However in this case,
LAMMPS will allow for an INEXACT match of a shared point to make a
connection between 2 triangles or lines.  Again, this is a single
point in both for a corner point or end point connection.  It is two
points in both triangles for an edge connection.

An INEXACT match means that the two points can be EPSILON apart.
EPSILON is defined as a tiny fraction (1.0e-4) of the size of
the smallest triangle or line in the system.

The reason INEXACT matches are allowed is that data files can be
created in a variety of manners, including by LAMMPS itself as a
simulation runs via the :doc:`write_data <write_data>` command.
Internally, triangle-style and line-style particles do not store their
corner points directly.  Instead, the center point of the
triangle/line is stored, along with an orientation of the
triangle/line and a displacement vector from the center point for each
corner point.  This means that when new corner points values are
written to a data file for two different triangles/line, they may
differ by epsilon due to round-offs in finite-precision arithmetic.

Note that due to how connectivity is defined, two triangles/lines will
not be connected if their corner points are separated by even small
distances (greater than EPSILON).  Likewise they will not be connected
if the corner point of one triangle/line is very close to (or even on)
the surface of another triangle or middle of another line segment.  In
general these kinds of granular surfaces could be problematic and
should be avoided, but LAMMPS does not check for these conditions.

In addition, note that connectivity is only defined between two
triangles/lines of the same molecule ID. This way surfaces of two
molecules can move independently, as described in the following
section.

Note that if a triangle or line segment has a free edge or free
corner/end point (not connected to any other triangle/line), granular
particles will still interact with the triangle/line if the nearest
contact point to the spherical particle center is on the free edge or
is the free corner/end point.

----------

Surface motion
""""""""""""""

By default, surface triangles/lines are motionless during a
simulation, whether they are *global* or *local*.  Triangles/lines
impart forces and torques to granular particles, but the inverse
forces/torques on the triangles/lines do not cause them to move.

However, triangles/lines can be made to move in a prescribed manner.
E.g. the rotation of 2d mixer blades in the example described above.
These two commands can be used for that purpose:

* :doc:`fix_modify move <fix_modify>` for *global* surfaces
* :doc:`fix move <fix_move>` for *local* surfaces

For *global* surfaces, the :doc:`fix_modify move <fix_modify>` command
can move a specified subset of the triangles/lines in various ways
(translation, rotation, etc).  Which triangles move is specified based
on the *molecule ID* of each triangle.  Molecule IDs are specified when
surfaces are defined by the :doc:`fix surface/global <fix_surface_global>`
command. They can also be defined by the :doc:`fix_modify mol/region
<fix_modify>` command.

For *local* surfaces, the :doc:`fix move <fix_move>` command can move
a specified subset of the triangles/lines in various ways
(translation, rotation, etc).  Which triangles move is specified based
on the group-ID argument to the :doc:`fix move <fix_move>` command.
Groups of *local* surfaces can be defined by the :doc:`group <group>`
command.

.. note::

   For an object defined by two or more connected triangles/lines, it is
   an error to assign a motion and not include all the connected
   triangles/lines, since this would break the connections.  LAMMPS checks
   this for *global* surfaces but only checks that the fix group of any
   instances of :doc:`fix move <fix_move>` include all or none of a
   set of connnected *local* triangles/lines.

----------

Calculation of forces
"""""""""""""""""""""

After generating the surface connectivity, LAMMPS classifies each connection
as being flat or non-flat based on the angular difference between normal vectors.
For each non-flat connection, the two sides are then classified as being concave
or convex.

Each point or edge of a line or triangle are then classified as being internal,
external, or unconnected based on the connectivity. For lines, an end point
is internal if it only has flat connections, external if it has at least one
non-flat (concave or convex) connection, and unconnected if it has no connections.
The same is true for edges on a triangle. Corners on triangles inherit their
classification from the two edges that meet at it. If either edge is unconnected,
the corner is unconnected. Otherwise, the corner is external if either edge is
external, and internal if both edges are internal.

If a particle is in contact with a single line or triangle, the force is applied
in the direction from the contact point to the particle center, :math:`\hat{n}_c`.
The magnitude of the force then depends on the details of the pair style, but may
be a function of the overlap distance :math:`\delta_f`, the relative (angular)
velocity, and the type of the sphere/surface.

If a particle is in contact with two or more connected lines or triangles,
:math:`\hat{n}_f` depends on details of the connection and the classification of
any points, edges, or corners involved in the contact. At a high level, this
process combines all connected flat surfaces in contact with the particle into a
single effective compsite surface that applies a single force on the particle.
Forces are also adjusted near convex or concave turns or unconnected points/edges
to prevent aphysical forces. For instance, surfaces on the other side of a convex
turn may be hidden such that do not apply a force. Details follow below.

To do this, LAMMPS first identifies which surface is closest to the particle. This
primary contact is effectively the surface with the largest overlap with the
particle, however, other criteria are used to break ties to ensure it is consistent
and does not depend on machine precision. LAMMPS then identifies which side of the
primary contact is in contact with the particle and recursively walks all
connections to determine which side of all other surfaces are in contact. If
there are unchecked surfaces, LAMMPS finds the next unchecked surface with the
largest overlap and repeats the process.

Now that all contacting surfaces know which side should be touching the particle,
LAMMPS performs a second iteration through contacting surfaces to calculate forces.
Starting from the primary contact, LAMMPS checks all connections. Flat connections
are added to a growing set of composite surfaces and added to a list of surfaces
which need to be walked. Concave connections are skipped while convex connections
are flagged to be hidden if the convex surface has a smaller overlap. LAMMPS then
iterates through the list of surfaces to be walked, effectively searching all
1st flat neighbors of the primary contact, before repeating the process with their
flat connections (2nd neighbors) and so forth until there are no more flat connections
in contact with the particle. LAMMPS then calculates a force from the composite
surface consisting of all flat connections (described in the following paragraph)
before finding the next surface with the largest overlap and repeats the process with
that surface. Continuing with this until all contacting surfaces have been processed
and have had forces calculated. While checking for flat connections to a surface that
was flagged as hidden due to a convex turn, LAMMPS passes along this flag to any such
flat connections that have not already been processed. Thus a set of flat surfaces
behind a convex turn can all be hidden as if they were a single surface. When
calculating forces, surfaces are hidden by simply zeroing their overlap
:math:`\delta_f`.

Given a set of mutually flat surfaces that make up a single composite surface,
LAMMPS calculates a single force using the maximum overlap across constituent
surfaces :math:`\delta_f = \delta_\mathrm{max}` and an effective direction
:math:`\hat{n}_f` of the contact force. This unit vector is a weighted average
across all surfaces :math:`i` in the set :math:`\hat{n}_f = A \sum W_i \hat{n}_{f,i}`
where :math:`\hat{n}_{f,i}` is the calculated direction of force for that surface
and :math:`A` is a normalization factor.

Before calculating contributions and weights for individual surfaces, LAMMPS first
calculates a weight for externally and internally contacted surfaces defined as
:math:`W_\mathrm{ext} \equiv (\delta_\mathrm{max,ext} / \delta_\mathrm{max})^2` and
:math:`W_\mathrm{int} \equiv 1 - W_\mathrm{ext}`, respectively, where
:math:`\delta_\mathrm{max,ext}` is the maximum overlap with an external boundary
of a surface. As a reminder, external boundaries are points in 2D or edges in 3D
between two non-flat surfaces. Corners in 3D inherit an external status from
either of the two edges that meet at it. This weighting is used amplify contributions
from externally-contacted surfaces as the particle moves around a concave or
convex turn to smoothly interpolate :math:`\hat{n}_f`.

In 2D, a surface's weight :math:`W_i` is simply the particle's overlap with
that line :math:`delta_i` multiplied by either :math:`W_\mathrm{ext}` or
or :math:`W_\mathrm{int}`. By default, :math:`\hat{n}_{f,i}` is simply the
direction from that contact point to the particle :math:`\hat{n}_{r,i}`.
If the contact point is inside of the line, :math:`\hat{n}_{r,i}` is
equivalent to the surface normal :math:`\hat{n}_{s,i}`. However, if the
contact point is at a concave corner, then
:math:`\hat{n}_{f,i} = \hat{n}_{s,i}`. If the contact point is at a convex
corner and :math:`\hat{n}_{r,i}` has a component pointing into the the
adjacent line :math:`j`, then set :math:`\hat{n}_{f,i} = \hat{n}_{s,j}`.

.. figure:: img/gransurf_nonflat_turn.jpg
            :figwidth: 33%
            :align: right
            :target: _images/gransurf_nonflat_turn.jpg

Contributions in 3D (no free edges)

Contributions in 3D (with free edges)

----------

Recommendations for geometries
""""""""""""""""""""""""""""""

When designing geometries for granular surfaces, there are several things
to keep in mind to avoid unintended force contributions and to improve
efficiency.

While convex corners are identified and used to censor forces from
physically hidden sections of a wall, if a particle is not in contact
with the entirety of a convex turn, then forces cannot be properly
censored. For example, consider a 2d system with a U shaped wall
defined by 3 line segments (see figure). If the width of the U is wider
than the typical particle-wall overlap (right), no issues are anticipated.
However, if the width of the U is very thin relative to the typical
particle-wall overlap (middle, right), then a particle could potentially
be in contact with both vertical legs of the U. If the particle is also
in contact with the base of the U (middle, green), then it can identify
the convex turn and censor forces from the rightmost vertical leg (blue).
However, if the particle is towards the top of the U (right) and not in
contact with the base (green), then there is no way for the particle to
identify the convex turn and censor forces from the far vertical wall
(blue). Therefore, walls with very thin features separated by a gap less
than the expected overlap distance between a particle and a surface can
lead to unintended additional forces.

.. figure:: img/gransurf_U_geom.jpg
            :figwidth: 33%
            :align: right
            :target: _images/gransurf_U_geom.jpg

As mentioned in the above section, forces resulting from contact with
unconnected endpoints of lines do point in the expected direction and
experience no discontinuities as the particle moves around the
endpoint. However, in 3D, contacts with unconnected edges only produce
reasonably directed forces oriented away from the edge. However, the
exact direction of a force can wobble as the contact moves across a
series of disconnected edges and convex turns may not be appropriately
censored as LAMMPS does not currently construct the proper quasi-2d
connectivity of 2D features. Therefore, it is recommended to avoid
complex geometries along unconnected boundaries such as rapid
oscillations in- or out-of-plane such as pleats or sawteeth, relative
to the length of an edge.

To build a neighborlist between particles and lines/triangles, LAMMPS
assigns a radius to each line/triangle that corresponds to the radius
of the circle/sphere that encloses the object. Therefore, one must be
aware that systems with large disparities between triangle/line and
partile radii may slow down the neighbor list build for and could
benefit from :doc:`the multi neighbor style <neighbor>` for *local*
surfaces. Furthermore, triangles with very high aspect ratios should
generally be avoided as they can lead to large neighbor lists
containing many particles which are not actually close to being
contact with the triangle (false positives).

----------

Example scripts
"""""""""""""""

The ``examples/gransurf`` directory has example input scripts which use
both *global* and *local* surfaces.  Both 2d and 3d models are included.

Each script produces a series of snapshot images using the :doc:`dump
image <dump_image>` command.  The snapshots visualize both the particles
and granular surfaces.  The snapshots can be animated to view a movie of
the simulation.
