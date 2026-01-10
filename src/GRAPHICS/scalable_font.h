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

/* error codes */
#define SSFN_OK 0            /* success */
#define SSFN_ERR_ALLOC 1     /* allocation error */
#define SSFN_ERR_NOFACE 2    /* no font face selected */
#define SSFN_ERR_INVINP 3    /* invalid input */
#define SSFN_ERR_BADFILE 4   /* bad SSFN file format */
#define SSFN_ERR_BADSTYLE 5  /* bad style */
#define SSFN_ERR_BADSIZE 6   /* bad size */
#define SSFN_ERR_BADMODE 7   /* bad mode */
#define SSFN_ERR_NOGLYPH 8   /* glyph (or kerning info) not found */

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
  uint32_t fragments_offs; /* offset of fragments table relative to magic */
  uint32_t characters_offs; /* offset of characters tables relative to magic */
  uint32_t kerning_offs;  /* kerning table offset relative to magic */
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
  int err;                                /* returned error code */
  int family;                             /* required family */
  int style;                              /* required style */
  int size;                               /* required size */
  int mode;                               /* required mode */
  int g;                                  /* shift value for grid size */
  int m, ix, u, uix, uax, lx, ly, mx, my; /* helper variables */
};

/* normal renderer */
extern int ssfn_load(ssfn_t *ctx, const ssfn_font_t *font);           /* add an SSFN to context */
extern int ssfn_select(ssfn_t *ctx, int family, int style, int size); /* select font to use */
extern uint32_t ssfn_utf8(char **str);             /* decode UTF-8 sequence */
extern ssfn_glyph_t *ssfn_render(ssfn_t *ctx, uint32_t unicode); /* return allocated glyph bitmap */
extern int ssfn_kern(ssfn_t *ctx, uint32_t unicode, uint32_t nextunicode, int *x,
                     int *y); /* get kerning values */
extern int ssfn_bbox(ssfn_t *ctx, char *str, int usekern, int *w,
                     int *h);                   /* get bounding box of a rendered string */
extern int ssfn_mem(ssfn_t *ctx);               /* return how much memory is used */
extern void ssfn_free(ssfn_t *ctx);             /* free context */
#define ssfn_lasterr(ctx) ((ssfn_t *) ctx)->err /* return last error code */
#define ssfn_error(err) \
  (err >= 0 && err <= 9 ? ssfn_errstr[err] : "Unknown error") /* return string for error code */
extern const char *ssfn_errstr[];
extern const ssfn_font_t *const ssfn_sans_font;
}    // namespace SSFN

#endif
