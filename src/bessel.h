#ifndef HZ_BESSEL_H
#define HZ_BESSEL_H

#include <complex.h>

/* Cylinder functions for the 2D exterior: the DtN map on a circle and the
 * exact Mie reference both live here (SLICE_PLAN step 2a).
 *
 * WHY THIS IS ITS OWN FILE WITH ITS OWN TESTS. The DtN operator is diagonal in
 * angular harmonics with symbol k*H'_m(kR)/H_m(kR), and the whole termination
 * rests on that one number being right for every m up to ~kR + margin. libm's
 * jn/yn are accurate in the range this project needs but not everywhere, so the
 * Wronskian and the recurrence are checked as gates rather than assumed
 * (SLICE_PLAN doubt 4: "written as one line, as if trivial. It is not"). */

/* J_m(x), Y_m(x) and their derivatives, x > 0. Derivatives from the standard
 * relation f'_m = (f_{m-1} - f_{m+1})/2. */
void hz_bessel_jy(int m, double x, double *j, double *y, double *dj, double *dy);

/* H^(1)_m'(x) / H^(1)_m(x) — the DtN symbol, without the leading k.
 * Finite for every m >= 0, x > 0 within the validated range (see
 * tests/test_bessel.c: |m| <= 3x, x up to 1e3). Beyond it Y_m overflows and the
 * function returns the large-order limit -m/x, which is the correct evanescent
 * behaviour but only to leading order — the test pins where that starts. */
double complex hz_hankel_ratio(int m, double x);

/* Wronskian residual |J_m(x)Y'_m(x) - J'_m(x)Y_m(x) - 2/(pi x)| / (2/(pi x)).
 * Exposed because it is the cheapest independent statement about the accuracy
 * of the pair at a given (m, x), and the slice bench prints it for the harmonic
 * range it actually uses. */
double hz_bessel_wronskian_err(int m, double x);

#endif
