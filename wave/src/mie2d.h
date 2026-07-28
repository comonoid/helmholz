#ifndef HZ_MIE2D_H
#define HZ_MIE2D_H

#include <complex.h>

/* Exact scattering by a RADIALLY SYMMETRIC scatterer in 2D — the reference for
 * gate D2 of the vertical slice (SLICE_PLAN step 2b).
 *
 * ONE CODE PATH FOR BOTH SCENES. Separation of variables is exact for any
 * radial k(r): the exterior is always
 *     u = u_inc + sum_m i^m c_m H^(1)_m(k0 r) e^{i m theta},
 * and everything the scatterer does enters through ONE number per harmonic —
 * the interior logarithmic derivative L_m = R'_m(a)/R_m(a):
 *     c_m = [L J_m(k0 a) - k0 J'_m(k0 a)] / [k0 H'_m(k0 a) - L H_m(k0 a)].
 * For a sharp cylinder L is closed form (Bessel of the interior wavenumber);
 * for a graded profile it comes from integrating a 1D ODE to machine accuracy.
 * A 1D ODE is NOT a grid in the sense audit A4 of M9 forbids: its cost does not
 * grow with the domain and its accuracy is checkable by step halving.
 *
 * NO GRID ANYWHERE, so the reference stays valid at any ka — which is the whole
 * reason the cylinder was chosen as the gate in the first place. */

typedef enum {
  HZ_MIE_SHARP = 0,    /* k^2 = k0^2 * n^2 for r < a, jump at r = a */
  HZ_MIE_GRADED = 1,   /* k^2 = k0^2 * (1 + d*(1-(r/a)^2)^3), C^2 at r = a */
  HZ_MIE_SHARP_ODE = 2 /* the sharp profile through the ODE path: exists only so
                        * that the closed form and the integrator can be made to
                        * disagree, which is the only way to trust either */
} hz_mie_profile;

typedef struct {
  hz_mie_profile prof;
  double a;     /* radius */
  double k0;    /* background wavenumber */
  double param; /* SHARP: refractive index n. GRADED: peak d in k^2/k0^2 - 1 */
  int mmax;     /* harmonics kept, |m| <= mmax */
  int nr;       /* radial steps for the ODE (GRADED only) */
} hz_mie;

/* k^2(r) of the profile, exposed because the solver's medium must be the SAME
 * function the reference used — a reference and a scene that disagree about the
 * scatterer is artefact 12 waiting to happen. */
double hz_mie_k2(const hz_mie *s, double r);

/* Scattering coefficients c_m, m = 0..mmax (c_{-m} = c_m for axial incidence).
 * Returns 0 on success. */
int hz_mie_coeffs(const hz_mie *s, double complex *c);

/* Far-field pattern f(theta) = sum_m c_m e^{i m theta}, with u_s ~
 * sqrt(2/(pi k r)) e^{i(kr - pi/4)} f(theta). */
double complex hz_mie_far(const hz_mie *s, const double complex *c, double theta);

/* Scattered field at a point OUTSIDE the scatterer (r >= a). */
double complex hz_mie_scat(const hz_mie *s, const double complex *c, double x, double y);

/* First Born approximation of the same far-field pattern — an INDEPENDENT
 * reference, valid only at weak contrast. The two must not agree to a fixed
 * tolerance (Born's error is O(contrast) by construction); what must hold is
 * that the difference falls LINEARLY with the contrast. That rate is the real
 * test, and it is stronger than any magnitude threshold. */
double complex hz_mie_born(const hz_mie *s, double theta);

#endif
