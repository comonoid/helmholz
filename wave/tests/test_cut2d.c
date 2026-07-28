/* Unit tests for src/cut2d.c — the exact oblique-cut integral.
 *
 * WHY THIS FILE IS LONGER THAN THE CODE IT TESTS. M14 asks whether a cut
 * survives the move from 1D to 2D, and its headline branch ("exact cut =>
 * machine zero") is an ALGEBRAIC IDENTITY that holds iff these integrals are
 * exact — so a green M14 would be a statement about this file, not about the
 * basis (M14 audit, R2). The integrator therefore gets its own falsifiers:
 *   V1  full strip must reproduce the separable product from src/phi.c, which
 *       is independently tested (980/980) — this checks polygon, antiderivative
 *       continuity, edge splitting and moments at once;
 *   V2  an axis-aligned cut must reproduce the separable product with a clipped
 *       range — the degenerate case of the oblique code path;
 *   V3  additivity across the cut: the two sides must sum to the whole;
 *   V4  the genuinely OBLIQUE case against an independent witness (exact in x,
 *       Gauss in y) — the only check that does not share code with the object
 *       under test on the geometry;
 *   V5  transpose invariance, which exercises BOTH primitive-axis branches on
 *       the same mathematical quantity;
 *   V6  strict zero when the strip misses the rectangle.
 *
 * THE WITNESS IS THE WEAKER PARTY, MEASURED NOT ASSUMED. Where the strip
 * boundary crosses a corner of the rectangle the witness integrand has a kink,
 * and its Gauss-in-y error falls off as roughly NSUB^-2.5: on the plane-wave
 * band it read 2.1e-7 / 3.7e-8 / 1.2e-9 / 1.1e-10 at NSUB = 200 / 800 / 3200 /
 * 12800, converging onto the closed form. So V4 bounds the witness, not the
 * code, and NSUB is set where that bound is comfortably below the threshold. */
#include "cut2d.h"
#include "phi.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>

static int pass_count = 0, fail_count = 0;

static void check(const char *name, int ok, double v) {
  if (ok) {
    pass_count++;
    printf("  ok   %-52s %.3e\n", name, v);
    return;
  }
  fail_count++;
  printf("  FAIL %-52s %.3e\n", name, v);
}

/* Error normalised by the UNOSCILLATED size of the integral, not by the result.
 * At omega*W >> 1 the result is small purely by oscillatory cancellation (1.1e-7
 * against a natural scale of 15.7 in the om=40 row below), so a relative-to-
 * result metric reports the dynamic range of the integrand rather than an error
 * — the same trap as the annihilation residuals in test_carrier_op. Measured on
 * the first run: the om=40 row was 8.2e-12 "relative" and 5.8e-20 to scale. */
static double relsc(double complex a, double complex b, double scale) {
  return cabs(a - b) / scale;
}

static double axis_scale(double lo, double hi, hz_axis2 f) {
  if (f.nf == 2) return fabs(hz_phi_prod_integral(lo, hi, f.f[0], f.f[1]));
  if (f.nf == 1) return fabs(creal(hz_phi_integral_osc(lo, hi, f.f[0], 0.0)));
  return hi - lo;
}

/* ---- independent witness: exact in x, Gauss-Legendre in y --------------- */
enum { NG = 8, NSUB = 3200 };
static const double GX[NG] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                              -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                              0.7966664774136267,  0.9602898564975363};
static const double GW[NG] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                              0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                              0.2223810344533745, 0.1012285362903763};

static double complex xint(double a, double b, hz_axis2 fx, double complex om) {
  if (!(a < b)) return 0.0;
  if (fx.nf == 2) return hz_phi_prod_integral_osc(a, b, fx.f[0], fx.f[1], om);
  if (fx.nf == 1) return hz_phi_integral_osc(a, b, fx.f[0], om);
  if (cabs(om) < 1e-12) return b - a;
  double complex i = CMPLX(0.0, 1.0);
  return (cexp(i * om * b) - cexp(i * om * a)) / (i * om);
}

static double complex witness(double xlo, double xhi, double ylo, double yhi, hz_axis2 fx,
                              hz_axis2 fy, double complex omx, double complex omy, hz_strip2 st) {
  double complex tot = 0.0;
  for (int s = 0; s < NSUB; s++) {
    double ya = ylo + (yhi - ylo) * (double)s / (double)NSUB;
    double yb = ylo + (yhi - ylo) * (double)(s + 1) / (double)NSUB;
    double hw = 0.5 * (yb - ya), mid = 0.5 * (ya + yb);
    for (int g = 0; g < NG; g++) {
      double y = mid + hw * GX[g];
      double fyv = 1.0;
      for (int i = 0; i < fy.nf; i++)
        fyv *= hz_phi(y / fy.f[i].h - fy.f[i].n);
      double xa = xlo, xb = xhi;
      if (fabs(st.ca) > 1e-14) {
        double lo = (st.slo - y * st.sa) / st.ca, hi = (st.shi - y * st.sa) / st.ca;
        if (st.ca < 0.0) {
          double t = lo;
          lo = hi;
          hi = t;
        }
        if (lo > xa) xa = lo;
        if (hi < xb) xb = hi;
      } else {
        double sv = y * st.sa;
        if (sv < st.slo || sv > st.shi) continue;
      }
      tot += GW[g] * hw * fyv * cexp(CMPLX(0.0, 1.0) * omy * y) * xint(xa, xb, fx, omx);
    }
  }
  return tot;
}

int main(void) {
  printf("test_cut2d: exact integration over a rectangle cut by an oblique strip\n");
  double W = 1.0;
  hz_phi_factor a0 = {W, 0.0, 0}, a1 = {W, 1.0, 0}, b0 = {W, 0.0, 0}, b1 = {W, -1.0, 0};
  hz_axis2 fx = {{a0, a1}, 2}, fy = {{b0, b1}, 2};
  hz_axis2 f1 = {{a0, a0}, 1}, f0 = {{a0, a0}, 0};
  double xlo = -2.0 * W, xhi = 3.0 * W, ylo = -3.0 * W, yhi = 2.0 * W;
  hz_strip2 all = {1.0, 0.0, -HZ_CUT_INF, HZ_CUT_INF};
  double sc = axis_scale(xlo, xhi, fx) * axis_scale(ylo, yhi, fy);
  double sc0 = axis_scale(xlo, xhi, f0) * axis_scale(ylo, yhi, f0);

  /* --- V1: no cut at all must reproduce the separable closed form --------- */
  {
    static const double OM[5] = {0.0, 1e-6, 0.7, 5.0, 40.0};
    for (int i = 0; i < 5; i++) {
      double complex omx = OM[i], omy = 0.3 * OM[i] - 0.11;
      double complex got = hz_cut2d_integral(xlo, xhi, ylo, yhi, fx, fy, omx, omy, all);
      double complex want = hz_phi_prod_integral_osc(xlo, xhi, a0, a1, omx) *
                            hz_phi_prod_integral_osc(ylo, yhi, b0, b1, omy);
      char nm[64];
      snprintf(nm, sizeof nm, "V1 full strip == separable, om=%.0e", OM[i]);
      check(nm, relsc(got, want, sc) < 1e-14, relsc(got, want, sc));
    }
    double complex g1 = hz_cut2d_integral(xlo, xhi, ylo, yhi, f1, f0, 2.0, 3.0, all);
    double complex w1 = hz_phi_integral_osc(xlo, xhi, a0, 2.0) *
                        (cexp(CMPLX(0.0, 3.0) * yhi) - cexp(CMPLX(0.0, 3.0) * ylo)) /
                        CMPLX(0.0, 3.0);
    double s1 = axis_scale(xlo, xhi, f1) * axis_scale(ylo, yhi, f0);
    check("V1 nf=1 x nf=0 == closed form", relsc(g1, w1, s1) < 1e-14, relsc(g1, w1, s1));
  }

  /* --- V2: an axis-aligned cut is the separable form with a clipped range -- */
  {
    double xc = 0.37 * W;
    hz_strip2 left = {1.0, 0.0, -HZ_CUT_INF, xc};
    double complex omx = 3.3, omy = -1.7;
    double complex got = hz_cut2d_integral(xlo, xhi, ylo, yhi, fx, fy, omx, omy, left);
    double complex want = hz_phi_prod_integral_osc(xlo, xc, a0, a1, omx) *
                          hz_phi_prod_integral_osc(ylo, yhi, b0, b1, omy);
    check("V2 axis-aligned cut (x) == clipped separable", relsc(got, want, sc) < 1e-14,
          relsc(got, want, sc));
    hz_strip2 down = {0.0, 1.0, -HZ_CUT_INF, xc};
    got = hz_cut2d_integral(xlo, xhi, ylo, yhi, fx, fy, omx, omy, down);
    want = hz_phi_prod_integral_osc(xlo, xhi, a0, a1, omx) *
           hz_phi_prod_integral_osc(ylo, xc, b0, b1, omy);
    check("V2 axis-aligned cut (y) == clipped separable", relsc(got, want, sc) < 1e-14,
          relsc(got, want, sc));
  }

  /* --- V3: the two sides of an oblique cut must add up to the whole ------- */
  {
    double ca = cos(0.6), sa = sin(0.6), c = 0.41 * W;
    hz_strip2 lo = {ca, sa, -HZ_CUT_INF, c}, hi = {ca, sa, c, HZ_CUT_INF};
    static const double OM[3] = {0.0, 2.5, 30.0};
    for (int i = 0; i < 3; i++) {
      double complex omx = OM[i] + 0.3, omy = -0.6 * OM[i];
      double complex whole = hz_cut2d_integral(xlo, xhi, ylo, yhi, fx, fy, omx, omy, all);
      double complex s1 = hz_cut2d_integral(xlo, xhi, ylo, yhi, fx, fy, omx, omy, lo);
      double complex s2 = hz_cut2d_integral(xlo, xhi, ylo, yhi, fx, fy, omx, omy, hi);
      char nm[64];
      snprintf(nm, sizeof nm, "V3 oblique halves sum to whole, om=%.1f", OM[i]);
      check(nm, relsc(s1 + s2, whole, sc) < 1e-14, relsc(s1 + s2, whole, sc));
    }
  }

  /* --- V4: the oblique case against the independent witness --------------- */
  {
    static const double AL[3] = {0.3, 0.785398163397448, 1.2};
    for (int i = 0; i < 3; i++) {
      double ca = cos(AL[i]), sa = sin(AL[i]);
      hz_strip2 st = {ca, sa, -HZ_CUT_INF, 0.23 * W};
      double complex omx = 2.1, omy = -3.4;
      double complex got = hz_cut2d_integral(xlo, xhi, ylo, yhi, fx, fy, omx, omy, st);
      double complex ref = witness(xlo, xhi, ylo, yhi, fx, fy, omx, omy, st);
      char nm[64];
      snprintf(nm, sizeof nm, "V4 oblique vs witness, alpha=%.2f", AL[i]);
      check(nm, relsc(got, ref, sc) < 1e-10, relsc(got, ref, sc));
    }
    /* a genuine STRIP (two parallel cuts), which is what a displaced cut needs */
    double ca = cos(0.9), sa = sin(0.9);
    hz_strip2 band = {ca, sa, -0.3 * W, 0.5 * W};
    double complex got = hz_cut2d_integral(xlo, xhi, ylo, yhi, fx, fy, 1.3, 0.9, band);
    double complex ref = witness(xlo, xhi, ylo, yhi, fx, fy, 1.3, 0.9, band);
    check("V4 oblique BAND vs witness", relsc(got, ref, sc) < 1e-10, relsc(got, ref, sc));
    /* plane wave only (nf=0 on both axes): the pure geometric case */
    got = hz_cut2d_integral(xlo, xhi, ylo, yhi, f0, f0, 1.3, 0.9, band);
    ref = witness(xlo, xhi, ylo, yhi, f0, f0, 1.3, 0.9, band);
    check("V4 plane wave over a band vs witness", relsc(got, ref, sc0) < 1e-9,
          relsc(got, ref, sc0));
  }

  /* --- V5: transpose invariance exercises both primitive-axis branches ---- */
  {
    double ca = cos(0.7), sa = sin(0.7);
    hz_strip2 st = {ca, sa, -HZ_CUT_INF, 0.19 * W};
    hz_strip2 sw = {sa, ca, -HZ_CUT_INF, 0.19 * W};
    double complex g1 = hz_cut2d_integral(xlo, xhi, ylo, yhi, fx, fy, 0.4, 9.0, st);
    double complex g2 = hz_cut2d_integral(ylo, yhi, xlo, xhi, fy, fx, 9.0, 0.4, sw);
    check("V5 transpose invariance", relsc(g1, g2, sc) < 1e-14, relsc(g1, g2, sc));
  }

  /* --- V6: strictly zero when the strip misses the rectangle -------------- */
  {
    hz_strip2 far = {1.0, 0.0, 100.0 * W, 200.0 * W};
    double complex z = hz_cut2d_integral(xlo, xhi, ylo, yhi, fx, fy, 1.0, 2.0, far);
    check("V6 strictly zero when the strip misses", !(cabs(z) > 0.0), cabs(z));
    hz_strip2 empty = {1.0, 0.0, 1.0, 1.0};
    z = hz_cut2d_integral(xlo, xhi, ylo, yhi, fx, fy, 1.0, 2.0, empty);
    check("V6 strictly zero for an empty strip", !(cabs(z) > 0.0), cabs(z));
  }

  /* --- V7: a FACETED DISC against an exact analytic integral --------------
   * Int over |x| < a of e^{i q.x} dA = 2 pi a J1(qa)/q — a closed form that
   * knows nothing about polygons, so it measures the geometry alone. This is
   * the design rule of PLAN врезка A4 in its purest form: how many facets buy
   * how much accuracy, with no basis and no solver in the way. */
  {
    double ad = 1.0, q = 4.0;
    double exact = 2.0 * M_PI * ad * j1(q * ad) / q;
    /* normalised by the disc AREA, not by the integral: at qa = 4 the integral is
     * 0.104 against an area of 3.14, so a relative-to-result number would
     * overstate the geometry error thirtyfold — the same trap as relsc above */
    double area = M_PI * ad * ad;
    double prev0 = 0.0, prev1 = 0.0, rate0 = 0.0, rate1 = 0.0;
    printf("       facets   inscribed        mid-fit     (relative to the exact disc)\n");
    for (int mm = 8; mm <= 64; mm *= 2) {
      hz_half2 hp[HZ_CUT_MAXHP];
      double e[2];
      for (int fit = 0; fit < 2; fit++) {
        int nh = hz_cut2d_disc(0.0, 0.0, ad, mm, fit, hp);
        double complex got =
            hz_cut2d_poly(-1.2 * ad, 1.2 * ad, -1.2 * ad, 1.2 * ad, f0, f0, q, 0.0, hp, nh);
        e[fit] = fabs(creal(got) - exact) / area;
      }
      printf("       %6d   %10.3e   %10.3e\n", mm, e[0], e[1]);
      if (mm > 8) {
        rate0 = prev0 / e[0];
        rate1 = prev1 / e[1];
      }
      prev0 = e[0];
      prev1 = e[1];
    }
    check("V7 inscribed disc converges as 1/m^2", rate0 > 3.5 && rate0 < 4.6, rate0);
    check("V7 mid-fit disc converges at least as fast", rate1 > 3.5, rate1);
    /* the outside of a convex body is the whole minus the inside */
    hz_half2 hp[HZ_CUT_MAXHP];
    int nh = hz_cut2d_disc(0.1, -0.2, 0.7, 16, 1, hp);
    double complex whole = hz_cut2d_poly(xlo, xhi, ylo, yhi, fx, fy, 1.7, -0.9, hp, 0);
    double complex ins = hz_cut2d_poly(xlo, xhi, ylo, yhi, fx, fy, 1.7, -0.9, hp, nh);
    double complex ref = hz_cut2d_integral(xlo, xhi, ylo, yhi, fx, fy, 1.7, -0.9, all);
    check("V7 zero half-planes == the whole rectangle", relsc(whole, ref, sc) < 1e-14,
          relsc(whole, ref, sc));
    check("V7 inside is a strict part of the whole", cabs(ins) < cabs(whole),
          cabs(ins) / cabs(whole));
  }

  printf("test_cut2d: %d/%d passed\n", pass_count, pass_count + fail_count);
  return fail_count == 0 ? 0 : 1;
}
