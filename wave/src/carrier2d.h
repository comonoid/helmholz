#ifndef HZ_CARRIER2D_H
#define HZ_CARRIER2D_H

#include <complex.h>

/* 2D carrier (Trefftz) basis:
 *   B(x,y) = phi(x/W - nx) * phi(y/W - ny) * exp(i(kx(x - nx W) + ky(y - ny W)))
 *
 * WHY A TENSOR PRODUCT OF PLANE WAVES AND NOT SPHERICAL HARMONICS. The whole
 * assembly rests on separability: e^{i(kx x + ky y)} FACTORISES, so every 2D
 * integral is a product of two 1D integrals that are already closed-form and
 * already validated (hz_phi_prod_integral_osc, 980/980). Spherical harmonics
 * span the same physics and destroy exactly that.
 *
 * PHASE REFERENCE IS LOCAL on both axes - see the note in carrier.h. With a
 * real carrier that is a constant phase; with a complex one it is the
 * difference between a usable basis and an unusable one.
 *
 * BILINEAR, NOT SESQUILINEAR, like the 1D form: the test function is NOT
 * conjugated, so the exponents simply ADD and each axis reduces to one
 * oscillatory integral at omega = k_i + k_j. */
typedef struct {
  double W;          /* element width, same on both axes */
  int nx, ny;        /* node indices; support is W*(n-2) .. W*(n+2) per axis */
  double complex kx; /* signed carrier wavenumber, x component */
  double complex ky; /* signed carrier wavenumber, y component */
} hz_carrier2d;

/* Piecewise-constant medium on an axis-aligned rectangle: ABSOLUTE k^2.
 * Rectangles are what makes the medium term separable INSIDE a cell, which is
 * what keeps the assembly a product of 1D integrals. A quadtree hands over
 * exactly this. */
typedef struct {
  double x0, x1, y0, y1;
  double complex k2;
} hz_medrect2d;

double complex hz_carrier2d_val(hz_carrier2d b, double x, double y);

/* Galerkin entry  Int B_i (d2/dx2 + d2/dy2 + k^2(x,y)) B_j dA over the domain
 * [0,domx] x [0,domy].
 *   = Ixx*Iy0 + Ix0*Iyy + sum_cells (k2_cell - kx_j^2 - ky_j^2) * Ix0c * Iy0c
 * The medium rectangles must cover the domain; anything not covered is treated
 * as the carrier's own k^2, i.e. contributes nothing (the Trefftz cancellation).
 * ncalls (may be NULL) counts 1D closed-form integral evaluations. */
double complex hz_carrier2d_entry(hz_carrier2d bi, hz_carrier2d bj, double domx, double domy,
                                  const hz_medrect2d *rects, int nrect, long *ncalls);

/* Int B_i * f dA for a separable source f(x,y) = phi((x-xs)/ws) * phi((y-ys)/ws). */
double complex hz_carrier2d_rhs(hz_carrier2d bi, double domx, double domy, double xs, double ys,
                                double ws, long *ncalls);

#endif
