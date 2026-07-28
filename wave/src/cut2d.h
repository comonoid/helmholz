#ifndef HZ_CUT2D_H
#define HZ_CUT2D_H

#include "phi.h"
#include <complex.h>

/* EXACT integration of (piecewise polynomial) x (plane wave) over a rectangle
 * cut by an OBLIQUE strip. This is the machinery a cut costs in 2D.
 *
 * WHY IT CANNOT BE QUADRATURE. The integrand oscillates as e^{i(k_i+k_j).x}
 * across a support of 4W; at W = 1e5 lambda that is 1e6 periods, and W >> lambda
 * is precisely the regime the carrier basis exists for. Sampling is therefore
 * not an option in the assembly (M14 audit, R3) -- only as a witness at small W.
 *
 * WHY IT IS STILL CLOSED FORM. Separability is lost when the cut is oblique,
 * but ANALYTICITY is not: with F = (W(x) B(y) e^{i om_y y}, 0) and
 * W' + i om_x W = A(x), the divergence theorem turns the area integral into
 * edge integrals of (polynomial) x (exponential), each of which is closed form.
 * Horizontal edges carry n_x dl = dy = 0 and drop out, so a rectangle cut by a
 * strip leaves at most four contributing edges.
 * Cost is O(edges x breakpoints), independent of omega -- the same property that
 * makes hz_phi_prod_integral_osc affordable.
 *
 * The 3D case is the same theorem applied twice (polygon -> polyhedron); not
 * implemented here. */

/* Product of up to two phi factors on one axis; nf == 0 means the constant 1
 * (needed for plane-wave x plane-wave integrals: field norms and cross terms). */
typedef struct {
  hz_phi_factor f[2];
  int nf;
} hz_axis2;

#define HZ_CUT_INF 1e300

/* Strip  slo <= x*ca + y*sa <= shi  in the direction of the interface normal.
 * A half-plane is a strip with one end at +-HZ_CUT_INF. Two parallel cuts (the
 * true medium boundary and a DISPLACED basis cut) intersect to a strip, which
 * is why this is the primitive rather than a half-plane. */
typedef struct {
  double ca, sa;
  double slo, shi;
} hz_strip2;

/* Integral over {[xlo,xhi] x [ylo,yhi]} ^ {strip} of
 *     Fx(x) * Fy(y) * exp(i(omx*x + omy*y)) dA.
 * Exact to rounding. Returns 0 for an empty or degenerate region. */
double complex hz_cut2d_integral(double xlo, double xhi, double ylo, double yhi, hz_axis2 fx,
                                 hz_axis2 fy, double complex omx, double complex omy, hz_strip2 st);

/* Half-plane {x*ca + y*sa <= c} when keep_le, else {>= c}. */
typedef struct {
  double ca, sa, c;
  int keep_le;
} hz_half2;

#define HZ_CUT_MAXHP 96

/* The same integral over the rectangle intersected with ANY convex region given
 * as a list of half-planes — a faceted material boundary is exactly this. The
 * region OUTSIDE a convex body is not convex, so it is taken as the difference
 * of two calls (whole rectangle minus inside); that is linear and exact.
 * nhp <= HZ_CUT_MAXHP. */
double complex hz_cut2d_poly(double xlo, double xhi, double ylo, double yhi, hz_axis2 fx,
                             hz_axis2 fy, double complex omx, double complex omy,
                             const hz_half2 *hp, int nhp);

/* Half-planes of a regular m-gon approximating the circle of radius r about
 * (cx,cy). fit = 0 inscribes it (every facet falls inside the circle, so the
 * displacement is one-signed); fit = 1 pushes the facets out by the mean sagitta
 * so the displacement has ZERO MEAN — M14 [3] measured the error as
 * 0.125 (k d)^2 sqrt(d/W) for a UNIFORM offset, and a mean-zero one is the
 * cheaper half of that trade. Returns m. */
int hz_cut2d_disc(double cx, double cy, double r, int m, int fit, hz_half2 *hp);

/* The VERTICES of the very same m-gon, counter-clockwise, vx/vy of length m.
 * Facet i (half-plane i above) runs from vertex i to vertex i+1.
 * WHY THIS LIVES HERE. A Nitsche term integrates the jump ALONG the cut locus,
 * so the polygon it walks and the polygon the basis is cut by must be the same
 * curve to the last bit; a second copy of the "d = r cos(pi/m)" logic elsewhere
 * would let the two drift and the jump would then be measured where the
 * function is continuous. Returns m. */
int hz_cut2d_disc_verts(double cx, double cy, double r, int m, int fit, double *vx, double *vy);

/* mom[n] = integral over [-1,1] of u^n exp(i a u) du, n = 0..N, N <= 30.
 * Exported for src/nitsche2d.c, which integrates the same class of integrand
 * (piecewise polynomial x plane wave) along a LINE instead of over an area and
 * must not grow its own copy of the series/recurrence seam. */
void hz_cut2d_osc_moments(double complex a, int N, double complex *mom);

#endif
