.. index:: fix graphics/labels

fix graphics/labels command
===========================

Syntax
""""""

.. code-block:: LAMMPS

   fix ID group-ID graphics/labels Nevery mode keyword args ...

* ID, group-ID are documented in :doc:`fix <fix>` command
* graphics/labels = style name of this fix command
* Nevery = update graphics information every this many time steps
* zero or more keyword/args groups may be appended
* keyword = *image* or *text*

  .. parsed-literal::

     *image* filename x y z keyword args = display image in visualization
        filename = name of the image file
        x, y, z  = position where the center of the image is located in the visualization
        any of x, y, or z can be a variable (see below)
        one or more keyword/arg pairs may be appended
        keyword = *scale* or *transcolor*
          *scale* value = the image is scaled by this value (default 1.0), can be a variable (see below)
          *transcolor* arg = select color for transparency: *auto* or *none* or *r/g/b*
             *auto* = uses the color in the lower left corner of the image for transparency
             *none* = disables transparency
             *r/g/b* = provide three integers in the range 0 to 255 to select transparancy color in RGB color space

     *text* labeltext x y z keyword args = display text in visualization
        labeltext = text for the label, must be quoted if it contains whitespace
        x, y, z  = position where the center of the text is located in the visualization
        any of x, y, or z can be a variable (see below)

        keyword = *fontcolor* or *framecolor* or *backcolor* or *transcolor* or *size*
          *fontcolor* arg = select color for text: *white* (default) or *black* or *r/g/b*
             *white* = uses white
             *black* = uses black
             *r/g/b* = provide three integers in the range 0 to 255
          *framecolor* arg = select color for frame around text: *silver* (default) or *darkgray* or *white* or *black* or *r/g/b*
             *silver* = uses a very light gray
             *darkgray* = uses a very dark gray
             *white* = uses white
             *black* = uses black
             *r/g/b* = provide three integers in the range 0 to 255
          *backcolor* arg = select color for background of the text: *silver* (default) or *darkgray* or *white* or *black* *r/g/b*
             *silver* = uses a very light gray
             *darkgray* = uses a very dark gray
             *white* = uses white
             *black* = uses black
             *r/g/b* = provide three integers in the range 0 to 255
          *transcolor* arg = select color for transparency: *silver* (default) or *darkgray* or *white* or *black* or *none* or *r/g/b*
             *silver* = uses a very light gray
             *darkgray* = uses a very dark gray
             *white* = uses white
             *black* = uses black
             *none* = disables transparency
             *r/g/b* = provide three integers in the range 0 to 255
          *size* value = set the size of the characters (default 24), can be a variable (see below)

Examples
""""""""

.. code-block:: LAMMPS

   fix pix all graphics/labels 100 image teapot.png 5.0 -1.0 -2.0 transcolor auto
   fix pot all graphics/labels 100 image teapot.ppm 1.0 v_ypos v_zpos scale v_prog transcolor 19/92/192
   fix lbl all graphics/labels 1000 text "LAMMPS graphics demo" 5.0 -1.0 -2.0 backcolor darkgray framecolor black
   fix info all graphics/labels 1000 text "Step: $(step)  Angle: ${rot}" 5.0 -1.0 -2.0 size 32

Description
"""""""""""

.. versionadded:: TBD

This fix allows to add either images or text as "labels" to :doc:`dump
image <dump_image>` created images by using the *fix* keyword.  This can
be useful to augment images with additional graphics or text directly
and without having to post-process the images.  Since the positioning
uses the coordinate system of the simulation and because the graphics
objects use the depth buffer of the image rasterizer, atoms and other
graphics in the "scene" can be located before or behind any text or
image label.

The *group-id* is ignored by this fix.

The *Nevery* keyword determines how often the graphics data is updated.
This should be the same value as the corresponding *N* parameter of the
:doc:`dump <dump>` image command.  LAMMPS will stop with an error
message if the settings for this fix and the dump command are not
compatible.

The *image* keyword reads an image file and adds it to the visualization
centered around the provided position and optionally scaled by the
provided scale factor.  The filename suffix determines whether LAMMPS
will try to read a file in JPEG, PNG, or PPM format.  If the suffix is
".jpg" or ".jpeg", then LAMMPS attempts to read the image in `JPEG
format <jpeg_format_>`_, if the suffix is ".png", then lammps attempts
to read the image in `PNG format <png_format_>`_.  Otherwise LAMMPS will
try to read the image in `ppm (aka netpbm) format <ppm_format_>`_.  Not
all variants of those file formats are compatible with image reader code
in LAMMPS.  If LAMMPS encounters an incompatible or unrecognizable file
format or a corrupted file, it will stop with an error.

When using the *image* keyword, the name of the image file and its position
in the "scene" are required arguments.  Optional keyword / value pairs
may be added:

  The *scale* value determines if the image is scaled before it is added
  to the :doc:`dump image <dump_image>` output.  LAMMPS currently
  employs a bilinear scaling algorithm.

  The *transcolor* value selects a color for transparency.  All pixels
  in the image with that color will be skipped when the image is
  rendered.  The color is specified as an R/G/B triple with values
  ranging from 0 to 255 for each channel.  There are also two special
  arguments: *auto* will pick the color of the pixel in the lower left
  corner as transparency color and *none* will disable all transparency
  processing (this is the default).


The *text* keyword will process a provided text into a pixmap and adds
it to the visualization centered around the provided position in a
similar fashion as with the *image* keyword.  The requirements for the
text argument are the same as in the :doc:`fix print <fix_print>`
command: it must be a single argument, so text with whitespace must be
quoted; and the text may contain equal style or immediate variables
using the ``${name}`` or ``$(expression)`` format.  The variables are
evaluated and expanded at every *Nevery* time step.

When using the *text* keyword, the text and its position in the "scene"
are required arguments.  Optional keyword / value pairs may be added:

  The *size* value determines the size of the letters in the text in
  pixels (approximately) and values between 4 and 512 are accepted.
  The default value is 24.

  There are four color settings: *fontcolor* or *framecolor* or
  *backcolor* or *transcolor*\ .  The color can be specified for all of
  those either as an R/G/B triple with values ranging from 0 to 255 for
  each channel (e.g. yellow would be "255/255/0").  There are also a few
  shortcuts for common choices: *silver*, *darkgray*, *white*, *black*.
  The default *backcolor* value is *silver*.

  - *fontcolor* selects the color for the text, default is *white*
  - *backcolor* selects the color for the background, default is
    *silver*
  - *framecolor* selects the color for the frame around the background,
    default is *silver*.
  - *transcolor* value selects a color for transparency, default is
    *silver*.  If this color is the same as any of the other color
    settings, those pixels are not drawn.  Thus with the default
    settings, the text will be rendered in white without background or
    frame.  The *none* setting for *transcolor* disables transparency
    processing.

  When rendering text with transparent background it is recommended to
  select a similar color but slightly darker or brighter color as background.
  This will reduce unwanted color effects at the edges due to anti-aliasing.

There may be multiple *image* or *text* keywords with their arguments
in a single fix *graphics/labels* command.

The arguments for the positions of an *image* or *text* and the *scale*
factor of an *image* or the *size* of a *text* can be specified as
equal-style :doc:`variables <variable>`, namely *x*, *y*, *z*, *scale*,
or *size*.  If any of these values is a variable, it should be specified
as `v_name`, where `name` is the variable name.  In this case, the
variable will be evaluated each *nevery* timestep, and its value used to
position and resize the image or text.  Please see the documentation of
the :doc:`fix graphics/objects <fix_graphics_objects>` command for a
more detailed discussion on using variables with graphics objects.

.. _jpeg_format: https://jpeg.org/jpeg/
.. _png_format: https://en.wikipedia.org/wiki/portable_network_graphics
.. _ppm_format: https://en.wikipedia.org/wiki/netpbm

-----------

Dump image info
"""""""""""""""

The fix graphics/labels command is designed to be used with the *fix*
keyword of :doc:`dump image <dump_image>`.  The fix adds images or text
to the visualization.

The color style setting for the fix in the :doc:`dump image
<dump_image>` has no effect on either image or text labels. The
transparency is by default fully opaque and can be changed with
*dump\_modify ftrans*\ .

The *fflag1* and *fflag2* settings of *dump image fix* are ignored.

Restart, fix_modify, output, run start/stop, minimize info
==========================================================

No information about this fix is written to :doc:`binary restart files
<restart>`.

None of the :doc:`fix_modify <fix_modify>` options apply to this fix.

Restrictions
""""""""""""

This fix is part of the GRAPHICS package.  It is only enabled if LAMMPS
was built with that package.  See the :doc:`Build package
<Build_package>` page for more info.

To read JPEG or PNG format images, support for the corresponding
graphics libraries must have been compiled and linked into LAMMPS.
Please see the :ref:`instructions for building LAMMPS with the GRAPHICS
package <graphics>` for more information on how to do that.

Related commands
""""""""""""""""

:doc:`fix print <fix_print>`,
:doc:`fix graphics/arrows <fix_graphics_arrows>`,
:doc:`fix graphics/isosurface <fix_graphics_isosurface>`,
:doc:`fix graphics/objects <fix_graphics_objects>`,
:doc:`fix graphics/periodic <fix_graphics_periodic>`

Default
"""""""

transcolor = "none" for *image* and "silver" for *text*, scale = 1.0, fontcolor = white, backcolor = silver, framecolor = silver, size = 24
