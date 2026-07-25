#ifndef HZ_CARRIER_H
#define HZ_CARRIER_H

#include <complex.h>

/* Carrier (Trefftz) basis element:  phi(x/W - n) * exp(i * kx * x),
 * optionally CUT by a material boundary.
 *
 * WHY THE MODULATION. (d2/dx2 + k^2)[phi e^{ikx}] = (phi'' + 2ik phi'
 * + (k^2 - kx^2) phi) e^{ikx}: when the carrier matches the local medium the
 * k^2 term cancels ANALYTICALLY, which is what makes an element thousands of
 * wavelengths wide a legal trial function. A plain smooth bump cannot do that.
 *
 * WHY A SIGNED kx RATHER THAN A cos/sin PAIR (M9 question 12). The pair
 * (phi cos, phi sin) spans the same space, but it MIXES +k and -k: you cannot
 * drop one propagation direction without constraining a pair of coefficients.
 * With one signed kx per function each direction is a separate basis function
 * and can be selected or dropped individually — which is exactly what M9c has
 * to measure. In 3D +k and -k are "light towards the camera" and "away from
 * it": physically distinct, so keeping them in one function is wrong.
 * kx also carries the MEDIUM: it is the local wavenumber, not a global one
 * (measured: a mismatched carrier costs 100% error at W >= 10 lambda).
 *
 * WHY THE CUT. A smooth element straddling a material boundary is spoiled over
 * its WHOLE support (measured: ~20% floor at any width), because one smooth
 * envelope cannot hold different waves on the two sides. Cutting gives each
 * side its own coefficient and its own carrier, and restores exact
 * representation (measured: machine zero up to 1e5 lambda). Cutting, NOT
 * shrinking: forbidding elements to cross boundaries would bound their size by
 * the distance to the nearest surface and destroy LOD. */
/* PHASE REFERENCE IS LOCAL: the function is phi(x/W - n) * exp(i kx (x - n W)),
 * NOT exp(i kx x). With a real carrier the difference is a constant phase and
 * nothing else. With a COMPLEX carrier it is the difference between a usable
 * basis and an unusable one: a global reference scales an element sitting at
 * x by exp(-Im(kx) x), so an absorbing layer a thousand wavelengths from the
 * origin arrives with amplitude e^-175 and the system cannot be solved at all.
 * Measured — the absorbing-layer scheme returned a zero field at full rank
 * until this was fixed. */
typedef struct {
  double W;          /* element width */
  int n;             /* grid offset; support is W*(n-2) .. W*(n+2) */
  double complex kx; /* SIGNED carrier wavenumber. COMPLEX: an absorbing medium
                      * has a complex k, and a real carrier cannot match it at
                      * all — measured, 95% residual at alpha = 0.05 over a
                      * 10-lambda element. Carries direction, medium and loss. */
  int side;          /* cut: -1 = whole support, 0 = left of xcut, 1 = right of xcut */
} hz_carrier;

/* Piecewise-constant medium: ABSOLUTE k^2 on [a,b). Absolute rather than a
 * correction to some global background, because there is no global background
 * any more — each element measures the medium against ITS OWN carrier. */
typedef struct {
  double a, b;
  double complex k2;
} hz_medseg;

double complex hz_carrier_val(hz_carrier b, double x, double xcut);

/* Intersection of the two supports, clipped to [0,dom] and to both cut sides.
 * Empty => lo >= hi (functions on opposite sides of a cut never overlap). */
void hz_carrier_overlap(hz_carrier bi, hz_carrier bj, double dom, double xcut, double *lo,
                        double *hi);

/* Galerkin entry Int B_i * (d2/dx2 + k^2(x)) B_j dx.
 *
 * NOTE — BILINEAR, NOT SESQUILINEAR: the test function is NOT conjugated. That
 * keeps this an exact change of coordinates from the (phi cos, phi sin) form
 * (the two are related by T^T A T, which conjugation would break), and keeps
 * the system complex-symmetric. The L2 inner product used for PROJECTION is a
 * different object and does conjugate — do not confuse them.
 * Consequence of the bilinear form: exponentials simply ADD, so every term is
 * one oscillatory integral at omega = kx_i + kx_j, and the four cos/sin cases
 * collapse into a single rule.
 * ncalls (may be NULL) counts closed-form integral evaluations. */
double complex hz_carrier_entry(hz_carrier bi, hz_carrier bj, double dom, double xcut,
                                const hz_medseg *segs, int nseg, long *ncalls);

/* Int B_i * f dx for a source f(x) = phi((x - xsrc)/ws). */
double complex hz_carrier_rhs(hz_carrier bi, double dom, double xcut, double xsrc, double ws,
                              long *ncalls);

#endif
