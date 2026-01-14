/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://www.lammps.org/, Sandia National Laboratories
   LAMMPS development team: developers@lammps.org

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_SCALABLE_FONT_H
#define LMP_SCALABLE_FONT_H

#include <exception>
#include <string>

#define SSFN_DATA_MAX 65536

/* rendering modes */
#define SSFN_MODE_NONE 0    /* just select the font to get the glyph (for kerning) */
#define SSFN_MODE_OUTLINE 1 /* return the glyph's outlines */
#define SSFN_MODE_BITMAP 2  /* render into bitmap */

/* font family group */
#define SSFN_FAMILY_SERIF 0
#define SSFN_FAMILY_SANS 1
#define SSFN_FAMILY_DECOR 2
#define SSFN_FAMILY_MONOSPACE 3
#define SSFN_FAMILY_HAND 4

/* font style flags */
#define SSFN_STYLE_REGULAR 0
#define SSFN_STYLE_BOLD 1
#define SSFN_STYLE_ITALIC 2

namespace SSFN {

using ssfn_font_t = struct _ssfn_font_t {
  uint8_t magic[4];  /* SSFN magic bytes */
  uint32_t size;     /* total size in bytes */
  uint8_t family;    /* font family group */
  uint8_t style;     /* font style, zero or OR'd SSFN_STYLE_BOLD and SSFN_STYLE_ITALIC */
  uint8_t quality;   /* quality, defines grid size, 0 - 8 */
  uint8_t features;  /* feature flags, OR'd SSFN_FEAT_* */
  uint8_t revision;  /* format revision, must be zero */
  uint8_t reserved0; /* must be zero */
  uint16_t reserved1;
  uint16_t baseline;  /* horizontal baseline in grid pixels */
  uint16_t underline; /* position of under line in grid pixels */
  uint16_t bbox_left; /* overall bounding box for all glyphs in grid pixels */
  uint16_t bbox_top;
  uint16_t bbox_right;
  uint16_t bbox_bottom;
  uint32_t fragments_offs;  /* offset of fragments table relative to magic */
  uint32_t characters_offs; /* offset of characters tables relative to magic */
  uint32_t kerning_offs;    /* kerning table offset relative to magic */
};

/* returned bitmap struct */
using ssfn_glyph_t = struct _ssfn_glyph_t {
  uint8_t mode;                /* returned glyph's data format */
  uint8_t baseline;            /* baseline of glyph, scaled to size */
  uint8_t w;                   /* width */
  uint8_t h;                   /* height */
  uint8_t adv_x;               /* advance x */
  uint8_t adv_y;               /* advance y */
  uint16_t pitch;              /* data buffer bytes per line */
  uint32_t *cmap;              /* pointer to color map */
  uint8_t data[SSFN_DATA_MAX]; /* data buffer */
};

/* renderer context */

using ssfn_t = struct _ssfn_t {
  const ssfn_font_t **fnt[5];             /* font registry */
  const ssfn_font_t *s;                   /* explicitly selected font */
  const ssfn_font_t *f;                   /* font selected by best match */
  ssfn_glyph_t *ret;                      /* glyph to return */
  uint16_t *p;                            /* outline points */
  uint16_t *r[256];                       /* raster for scanlines */
  uint16_t *h;                            /* auto hinting grid */
  int len[5];                             /* number of fonts in registry */
  int mp;                                 /* memory allocated for points */
  int np;                                 /* how many points actually are there */
  int nr[256];                            /* number of coordinates in each raster line */
  int family;                             /* required family */
  int style;                              /* required style */
  int size;                               /* required size */
  int mode;                               /* required mode */
  int g;                                  /* shift value for grid size */
  int m, ix, u, uix, uax, lx, ly, mx, my; /* helper variables */
};

/* normal font renderer API */
// add SSFN font to context
extern void ssfn_load(ssfn_t *ctx, const ssfn_font_t *font);
// select font to use
extern void ssfn_select(ssfn_t *ctx, int family, int style, int size);
// return allocated glyph bitmap
extern ssfn_glyph_t *ssfn_render(ssfn_t *ctx, uint32_t unicode);
// get kerning values
extern int ssfn_kern(ssfn_t *ctx, uint32_t unicode, uint32_t nextunicode, int *x, int *y);
// get bounding box of a rendered string
extern void ssfn_bbox(ssfn_t *ctx, char *str, int usekern, int *w, int *h);
// free context
extern void ssfn_free(ssfn_t *ctx);

/** Font renderer exception class */
class SSFNException : public std::exception {
 public:
  SSFNException() = delete;
  /** Thrown during font processing
   *
   * \param   file    source file where exception was thrown
   * \param   line    line in source file where exception was thrown
   * \param   flag    select error message */
  explicit SSFNException(const std::string &file, int line, int flag);

  /** Retrieve message describing the thrown exception
   *
   * This function provides the message that can be retrieved when the corresponding
   * exception is caught.
   *
   * \return  String with error message */
  const char *what() const noexcept override { return message.c_str(); }

 private:
  std::string message;
};

extern const ssfn_font_t *const ssfn_sans_font;
}    // namespace SSFN

#endif
