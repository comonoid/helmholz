#ifndef HZ_DTN2D_H
#define HZ_DTN2D_H

#include "carrier2d.h"
#include <complex.h>

/* Exact radiation condition on a circle, for the 2D carrier basis
 * (SLICE_PLAN step 2d).
 *
 * The exterior of a circle of radius R is spanned by outgoing cylinder waves,
 * so the map from u to du/dr is DIAGONAL in angular harmonics:
 *     d/dr u_m(R) = k0 * H'_m(k0 R) / H_m(k0 R) * u_m(R).
 * TERMINATION_REPORT recorded "in 2D the DtN is nonlocal" and stopped there;
 * nonlocal in SPACE it is, but diagonal in harmonics, which is what matters.
 *
 * WHAT THE PLAN DID NOT SAY AND THIS HEADER DOES. Applying the symbol is one
 * line; COUPLING it to a tensor-product basis is not:
 *  - the rows are EXTRA rows, not a boundary term. The Galerkin form here is
 *    strong (no integration by parts), so no boundary term arises to hang the
 *    condition on — exactly as the 1D bench used two extra rows;
 *  - each row is DENSE in the basis functions that reach the circle;
 *  - the harmonic of one basis function is not band-limited (phi is only C^1),
 *    but their SUM is analytic in theta, so the sampling error cancels in any
 *    combination that represents a smooth field. ntheta must still clear
 *    2*(k0 R + mmax);
 *  - harmonics above mmax are left UNCONSTRAINED and can collect rubbish; the
 *    caller must check the computed tail has decayed. */

typedef struct {
  double R;   /* radius of the truncation circle */
  double k0;  /* exterior wavenumber */
  int mmax;   /* harmonics |m| <= mmax */
  int ntheta; /* sample points on the circle */
} hz_dtn;

/* The m-th angular harmonic of one basis function on the circle, and of its
 * radial derivative. Exposed separately because the solver needs the harmonics
 * themselves — the exterior field, and with it the far-field pattern, is
 * reconstructed from them — and because a test that only ever sees the finished
 * row cannot tell a wrong symbol from a wrong derivative. */
void hz_dtn_harm(const hz_dtn *d, int m, hz_carrier2d b, double complex *uh, double complex *duh);

/* Row entry for harmonic m and basis function b:
 *     (d/dr - Z_m) applied to the m-th harmonic of B.
 * The whole row is the vector over basis functions; a coefficient vector that
 * describes a purely outgoing field annihilates every row. */
double complex hz_dtn_row(const hz_dtn *d, int m, hz_carrier2d b);

/* The symbol itself, k0 * H'_m(k0 R)/H_m(k0 R). */
double complex hz_dtn_symbol(const hz_dtn *d, int m);

#endif
