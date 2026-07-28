#ifndef HZ_PHI_H
#define HZ_PHI_H

#include <complex.h>

/* 1D potential phi: C^1 piecewise-quadratic bump, integer knots, support [-2,2].
 *   phi(t) = 2 - t^2        on [-1,1]
 *   phi(t) = (|t| - 2)^2    on [1,2] and [-2,-1]
 *   phi(t) = 0              outside
 * Properties (tested in tests/test_phi.c): sum of integer translates == 4
 * (partition of unity), alternating sum == 0 (Nyquist null mode), refinable:
 * phi(x) = 1/4*[phi(2x+2) + 2phi(2x+1) + 2phi(2x) + 2phi(2x-1) + phi(2x-2)].
 * 3D potentials are tensor products of scaled translates phi(x/h - n). */

double hz_phi(double t);
double hz_phi_d1(double t);
double hz_phi_d2(double t);

/* One factor of a 1D product integrand: phi^(deriv)(x/h - n).
 * h > 0 is the level spacing (dyadic power of two in real use, any positive
 * value works), n the grid offset, deriv 0, 1 or 2 (derivative d^k/dx^k, which
 * brings a 1/h^k factor). The plain Helmholtz operator needs only 0 and 2;
 * deriv 1 appears in the carrier (modulated) basis, where
 * (d^2/dx^2 + k^2)[phi*exp(ikx)] leaves a phi' term (see M9, carrier basis). */
typedef struct {
  double h;
  double n;
  int deriv;
} hz_phi_factor;

/* Coefficients of ONE factor on the piece containing x = m, in the local
 * coordinate s = x - m:  f(x) = q[0] + q[1]*s + q[2]*s^2. m must not sit exactly
 * on a knot (callers pass piece midpoints). Local coordinates keep the
 * expansion well-conditioned for narrow pieces far from the origin.
 * Exposed for src/cut2d.c, which needs the polynomial itself rather than an
 * integral of it; duplicating the piece logic there would let the two copies
 * drift apart. */
void hz_phi_local_poly(hz_phi_factor f, double m, double q[3]);

/* Exact integral over [a,b] of f1(x)*f2(x). Piecewise-polynomial closed form:
 * the integrand is a polynomial of degree <= 4 between knots of the two
 * factors; no quadrature involved. Returns exactly 0.0 when supports do not
 * intersect [a,b]. */
double hz_phi_prod_integral(double a, double b, hz_phi_factor f1, hz_phi_factor f2);

/* Exact integral over [a,b] of f1(x)*f2(x)*exp(i*omega*x) -- same piecewise
 * decomposition, but each piece is integrated against the oscillation in closed
 * form instead of being sampled. This is what makes coarse (>> wavelength)
 * elements affordable: assembly cost no longer scales with omega, so no
 * wavelength-sized grid sneaks back in through quadrature (M9 audit, A1).
 * Returns exactly 0.0 when supports do not intersect [a,b]; omega == 0.0
 * reproduces hz_phi_prod_integral. */
double complex hz_phi_prod_integral_osc(double a, double b, hz_phi_factor f1, hz_phi_factor f2,
                                        double complex omega);

/* Exact integral over [a,b] of f(x)*exp(i*omega*x) — ONE factor, not a pair.
 * Needed wherever a potential meets a plane wave rather than another potential:
 * the projection of a known field onto the carrier basis, and the closed-form
 * field of an extended source (a partial integral of the source against
 * exp(i k y)). Same piecewise machinery, integrand degree <= 2. */
double complex hz_phi_integral_osc(double a, double b, hz_phi_factor f, double complex omega);

#endif
