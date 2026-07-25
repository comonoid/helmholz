/* Unit tests for src/phi.c — the M1 acceptance list from PLAN.md. */
#include "phi.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>

static int g_fail = 0;
static int g_total = 0;

static void check_close(double got, double want, double tol, const char *what) {
  g_total++;
  if (!(fabs(got - want) <= tol)) {
    g_fail++;
    printf("FAIL: %s: got %.17g want %.17g (tol %g)\n", what, got, want, tol);
  }
}

/* --- independent comparator: Gauss-Legendre 5 per knot-free piece ---------
 * Exact for polynomials of degree <= 9, so exact for our degree-4 integrand.
 * Independent of the production path: evaluates phi pointwise instead of
 * expanding piece polynomials. */

static double feval(hz_phi_factor f, double x) {
  double t = x / f.h - f.n;
  if (f.deriv == 2) return hz_phi_d2(t) / (f.h * f.h);
  if (f.deriv == 1) return hz_phi_d1(t) / f.h;
  return hz_phi(t);
}

static const double GL5_X[5] = {0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640,
                                0.9061798459386640};
static const double GL5_W[5] = {0.5688888888888889, 0.4786286704993665, 0.4786286704993665,
                                0.2369268850561891, 0.2369268850561891};

/* Breakpoints of [a,b] at the knots of both factors — the test's own copy, so
 * a bug in the production splitter cannot hide by being shared. */
static int split_knots(double a, double b, hz_phi_factor f1, hz_phi_factor f2, double *pts) {
  const hz_phi_factor fs[2] = {f1, f2};
  int np = 0;
  pts[np++] = a;
  pts[np++] = b;
  static const double knot_t[4] = {-2.0, -1.0, 1.0, 2.0};
  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 4; k++) {
      double x = fs[i].h * (fs[i].n + knot_t[k]);
      if (x > a && x < b) pts[np++] = x;
    }
  }
  for (int i = 1; i < np; i++) {
    double x = pts[i];
    int j = i - 1;
    while (j >= 0 && pts[j] > x) {
      pts[j + 1] = pts[j];
      j--;
    }
    pts[j + 1] = x;
  }
  return np;
}

static double gl5_prod(double a, double b, hz_phi_factor f1, hz_phi_factor f2) {
  double pts[10];
  int np = split_knots(a, b, f1, f2, pts);
  double total = 0.0;
  for (int p = 0; p + 1 < np; p++) {
    double half = 0.5 * (pts[p + 1] - pts[p]);
    double mid = 0.5 * (pts[p] + pts[p + 1]);
    if (half <= 0.0) continue;
    double s = 0.0;
    for (int g = 0; g < 5; g++) {
      double x = mid + half * GL5_X[g];
      s += GL5_W[g] * feval(f1, x) * feval(f2, x);
    }
    total += s * half;
  }
  return total;
}

/* phi evaluated in long double, written out rather than calling hz_phi — the
 * oscillatory comparator must not share a single line with the production
 * path. */
static long double lfeval(hz_phi_factor f, long double x) {
  long double t = x / (long double)f.h - (long double)f.n;
  long double a = fabsl(t);
  if (a >= 2.0L) return 0.0L;
  if (f.deriv == 2) return (a <= 1.0L ? -2.0L : 2.0L) / ((long double)f.h * (long double)f.h);
  if (f.deriv == 1) {
    long double d = (a <= 1.0L) ? -2.0L * t : 2.0L * (t < 0.0L ? -1.0L : 1.0L) * (a - 2.0L);
    return d / (long double)f.h;
  }
  if (a <= 1.0L) return 2.0L - t * t;
  long double d = a - 2.0L;
  return d * d;
}

static const long double GLL_X[5] = {0.0L, -0.5384693101056830910363144L,
                                     0.5384693101056830910363144L, -0.9061798459386639927976269L,
                                     0.9061798459386639927976269L};
static const long double GLL_W[5] = {0.5688888888888888888888889L, 0.4786286704993664680412915L,
                                     0.4786286704993664680412915L, 0.2369268850561890875142640L,
                                     0.2369268850561890875142640L};

/* Oscillatory comparator: the same GL5, with the oscillation resolved by brute
 * subdivision until omega*subwidth <= 0.1 rad. Deliberately dumb and slow —
 * this is exactly the quadrature the production path must NOT need (M9 audit,
 * A1), kept only as an independent witness.
 * It runs in long double for a measured reason: in double it sums ~1e5 terms
 * whose cos/sin arguments reach ~4e3, and the per-term rounding alone (not the
 * summation — Kahan does not help) caps what it can witness at ~1.4e-16
 * absolute. That is coarser than the analytic path's own error, so a double
 * witness cannot judge it. In long double the witness lands at ~3e-19. */
static double complex gl5_prod_osc(double a, double b, hz_phi_factor f1, hz_phi_factor f2,
                                   double omega) {
  double pts[10];
  int np = split_knots(a, b, f1, f2, pts);
  long double complex total = 0.0L;
  for (int p = 0; p + 1 < np; p++) {
    long double lo = (long double)pts[p], hi = (long double)pts[p + 1];
    if (hi <= lo) continue;
    long ns = (long)(10.0 * (pts[p + 1] - pts[p]) * fabs(omega)) + 1;
    long double sub = (hi - lo) / (long double)ns;
    for (long k = 0; k < ns; k++) {
      long double half = 0.5L * sub;
      long double mid = lo + ((long double)k + 0.5L) * sub;
      long double complex s = 0.0L;
      for (int g = 0; g < 5; g++) {
        long double x = mid + half * GLL_X[g];
        long double ph = (long double)omega * x;
        s += GLL_W[g] * lfeval(f1, x) * lfeval(f2, x) * (cosl(ph) + I * sinl(ph));
      }
      total += s * half;
    }
  }
  return CMPLX((double)creall(total), (double)cimagl(total));
}

/* --- tests ---------------------------------------------------------------- */

static void test_partition_of_unity(void) {
  for (int i = 0; i <= 200; i++) {
    double x = -1.0 + 0.01 * (double)i;
    double sum = 0.0, alt = 0.0;
    for (int n = -4; n <= 4; n++) {
      double v = hz_phi(x - (double)n);
      sum += v;
      alt += (n % 2 == 0) ? v : -v;
    }
    check_close(sum, 4.0, 1e-14, "partition of unity");
    check_close(alt, 0.0, 1e-14, "alternating null mode");
  }
}

static void test_c1_continuity(void) {
  /* One-sided samples at distance eps from a knot differ by ~2*|phi'|*eps
   * (<= 4*eps: |phi'| <= 2) for phi, and by ~2*|phi''|*eps (= 4*eps) for phi';
   * 5*eps covers both plus rounding. */
  static const double knots[5] = {-2.0, -1.0, 0.0, 1.0, 2.0};
  double eps = 1e-7;
  for (int k = 0; k < 5; k++) {
    check_close(hz_phi(knots[k] - eps), hz_phi(knots[k] + eps), 5.0 * eps, "phi continuity");
    check_close(hz_phi_d1(knots[k] - eps), hz_phi_d1(knots[k] + eps), 5.0 * eps, "phi' continuity");
  }
}

static void test_two_scale(void) {
  /* phi(x) = 1/4 * [phi(2x+2) + 2phi(2x+1) + 2phi(2x) + 2phi(2x-1) + phi(2x-2)] */
  static const double w[5] = {1.0, 2.0, 2.0, 2.0, 1.0};
  for (int i = 0; i <= 400; i++) {
    double x = -2.5 + 0.0125 * (double)i;
    double fine = 0.0;
    for (int n = -2; n <= 2; n++)
      fine += w[n + 2] * hz_phi(2.0 * x - (double)n);
    check_close(0.25 * fine, hz_phi(x), 1e-14, "two-scale relation");
  }
}

static void test_closed_forms(void) {
  hz_phi_factor p = {1.0, 0.0, 0};
  hz_phi_factor pdd = {1.0, 0.0, 2};
  /* <phi,phi> = 92/15, <phi'',phi> = -Int(phi'^2) = -16/3 */
  check_close(hz_phi_prod_integral(-2.0, 2.0, p, p), 92.0 / 15.0, 1e-13, "<phi,phi>");
  check_close(hz_phi_prod_integral(-2.0, 2.0, pdd, p), -16.0 / 3.0, 1e-13, "<phi'',phi>");
  /* clipping to a subinterval: [0,2] gives exactly half of each by symmetry */
  check_close(hz_phi_prod_integral(0.0, 2.0, p, p), 46.0 / 15.0, 1e-13, "<phi,phi> half");
}

static void test_disjoint_supports(void) {
  hz_phi_factor f1 = {1.0, 0.0, 0};
  hz_phi_factor f2 = {1.0, 4.0, 0}; /* supports [-2,2] and [2,6]: touch at a point */
  g_total++;
  double v = hz_phi_prod_integral(-10.0, 10.0, f1, f2);
  if (v < 0.0 || v > 0.0) {
    g_fail++;
    printf("FAIL: disjoint supports must give exact 0, got %.17g\n", v);
  }
}

static void test_vs_quadrature(void) {
  /* Cross-level pairs, both derivs, asymmetric windows — the real assembly
   * shapes. Fixed deterministic list. */
  static const struct {
    hz_phi_factor f1, f2;
    double a, b;
  } cases[] = {
      {{1.0, 0.0, 0}, {1.0, 1.0, 0}, -2.0, 3.0},       /* same level, neighbors */
      {{1.0, 0.0, 0}, {1.0, 2.0, 0}, -2.0, 4.0},       /* same level, distance 2 */
      {{1.0, 0.0, 2}, {1.0, 1.0, 0}, -2.0, 3.0},       /* operator x value */
      {{1.0, 0.0, 2}, {1.0, 1.0, 2}, -2.0, 3.0},       /* operator x operator */
      {{2.0, 0.0, 0}, {0.5, 3.0, 0}, -4.0, 4.0},       /* level diff 2 */
      {{2.0, -1.0, 2}, {0.5, 2.0, 0}, -6.0, 3.0},      /* level diff 2, deriv */
      {{8.0, 0.0, 0}, {0.25, 10.0, 2}, 0.0, 16.0},     /* level diff 5 */
      {{16.0, 1.0, 0}, {0.125, 100.0, 0}, 10.0, 14.0}, /* level diff 7, far window */
      {{1.0, 0.0, 0}, {1.0, 0.5, 0}, -2.0, 2.5},       /* non-integer offset */
      {{4.0, 0.25, 2}, {0.5, 1.0, 2}, -8.0, 8.0},      /* both derivs, cross level */
  };
  int ncases = (int)(sizeof(cases) / sizeof(cases[0]));
  for (int i = 0; i < ncases; i++) {
    double got = hz_phi_prod_integral(cases[i].a, cases[i].b, cases[i].f1, cases[i].f2);
    double want = gl5_prod(cases[i].a, cases[i].b, cases[i].f1, cases[i].f2);
    double scale = fabs(want) > 1.0 ? fabs(want) : 1.0;
    check_close(got, want, 1e-12 * scale, "analytic vs GL5 quadrature");
  }
}

/* --- M9a: deriv 1 and the oscillatory weight -------------------------------- */

static void check_cclose(double complex got, double complex want, double tol, const char *what) {
  g_total++;
  if (!(cabs(got - want) <= tol)) {
    g_fail++;
    printf("FAIL: %s: got (%.17g,%.17g) want (%.17g,%.17g) (tol %g)\n", what, creal(got),
           cimag(got), creal(want), cimag(want), tol);
  }
}

/* Deterministic case list reused by the deriv-1 and oscillatory tests: the
 * shapes the carrier assembly actually produces (phi', cross-level pairs). */
static const struct {
  hz_phi_factor f1, f2;
  double a, b;
} OSC_CASES[] = {
    {{1.0, 0.0, 0}, {1.0, 1.0, 0}, -2.0, 3.0},   /* value x value, neighbors */
    {{1.0, 0.0, 1}, {1.0, 1.0, 0}, -2.0, 3.0},   /* phi' x value — the carrier term */
    {{1.0, 0.0, 1}, {1.0, 1.0, 1}, -2.0, 3.0},   /* phi' x phi' */
    {{1.0, 0.0, 2}, {1.0, 1.0, 1}, -2.0, 3.0},   /* phi'' x phi' */
    {{2.0, 0.0, 1}, {0.5, 3.0, 0}, -4.0, 4.0},   /* level diff 2, deriv 1 */
    {{2.0, -1.0, 2}, {0.5, 2.0, 1}, -6.0, 3.0},  /* level diff 2, mixed derivs */
    {{8.0, 0.0, 1}, {0.25, 10.0, 2}, 0.0, 16.0}, /* level diff 5 */
    {{1.0, 0.0, 0}, {1.0, 0.5, 1}, -2.0, 2.5},   /* non-integer offset */
    {{4.0, 0.25, 1}, {0.5, 1.0, 1}, -8.0, 8.0},  /* both deriv 1, cross level */
};
static const int NOSC = (int)(sizeof(OSC_CASES) / sizeof(OSC_CASES[0]));

static void test_deriv1_vs_quadrature(void) {
  for (int i = 0; i < NOSC; i++) {
    double got =
        hz_phi_prod_integral(OSC_CASES[i].a, OSC_CASES[i].b, OSC_CASES[i].f1, OSC_CASES[i].f2);
    double want = gl5_prod(OSC_CASES[i].a, OSC_CASES[i].b, OSC_CASES[i].f1, OSC_CASES[i].f2);
    double scale = fabs(want) > 1.0 ? fabs(want) : 1.0;
    check_close(got, want, 1e-12 * scale, "deriv1 analytic vs GL5");
  }
}

static void test_osc_vs_quadrature(void) {
  /* omega spread across both branches of osc_moments: the series regime
   * (|omega*w| < 4), the seam, and the closed-form regime. */
  static const double omegas[7] = {0.3, 1.7, 3.9, 4.1, 8.1, 25.0, 130.0};
  for (int i = 0; i < NOSC; i++) {
    for (int j = 0; j < 7; j++) {
      double complex got = hz_phi_prod_integral_osc(OSC_CASES[i].a, OSC_CASES[i].b, OSC_CASES[i].f1,
                                                    OSC_CASES[i].f2, omegas[j]);
      double complex want =
          gl5_prod_osc(OSC_CASES[i].a, OSC_CASES[i].b, OSC_CASES[i].f1, OSC_CASES[i].f2, omegas[j]);
      double scale = cabs(want) > 1.0 ? cabs(want) : 1.0;
      check_close(cabs(got - want), 0.0, 1e-11 * scale, "osc analytic vs GL5");
    }
  }
}

static void test_osc_zero_omega(void) {
  /* omega = 0 must reproduce the non-oscillatory path (not bitwise: the two
   * group the same polynomial moments differently). */
  for (int i = 0; i < NOSC; i++) {
    double complex got = hz_phi_prod_integral_osc(OSC_CASES[i].a, OSC_CASES[i].b, OSC_CASES[i].f1,
                                                  OSC_CASES[i].f2, 0.0);
    double want =
        hz_phi_prod_integral(OSC_CASES[i].a, OSC_CASES[i].b, OSC_CASES[i].f1, OSC_CASES[i].f2);
    double scale = fabs(want) > 1.0 ? fabs(want) : 1.0;
    check_cclose(got, want, 1e-15 * scale, "osc at omega=0 == plain");
  }
}

static void test_osc_seam(void) {
  /* Crossing |omega*w| = 4 must not produce a step. Pieces of the h=1 potential
   * have half-widths 1.0 ([-1,1]) and 0.5 ([-2,-1], [1,2]), so the seam is hit
   * at omega = 4 and omega = 8; step over both. */
  hz_phi_factor f1 = {1.0, 0.0, 0}, f2 = {1.0, 1.0, 1};
  static const double om[6] = {3.99, 4.0, 4.01, 7.99, 8.0, 8.01};
  for (int j = 0; j < 6; j++) {
    double complex got = hz_phi_prod_integral_osc(-2.0, 3.0, f1, f2, om[j]);
    double complex want = gl5_prod_osc(-2.0, 3.0, f1, f2, om[j]);
    check_close(cabs(got - want), 0.0, 1e-11, "osc across the series/closed-form seam");
  }
}

static void test_osc_large_omega(void) {
  /* omega*h up to 1e3 — the coarse-element regime M9 is built for. The integral
   * itself decays hard (smooth envelope against fast oscillation: down to
   * ~3.5e-9 from an O(1) integrand, six orders of cancellation between pieces),
   * so the meaningful assertion is ABSOLUTE at the integrand scale. Measured
   * against the long double witness: the analytic path stays within 3e-19, so
   * the 5e-17 floor below carries a ~150x margin for libm differences. */
  static const double omegas[3] = {200.0, 500.0, 1000.0};
  hz_phi_factor cases[3][2] = {
      {{1.0, 0.0, 0}, {1.0, 1.0, 0}},
      {{1.0, 0.0, 1}, {1.0, 1.0, 0}},
      {{2.0, 0.0, 2}, {1.0, 1.0, 1}},
  };
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      double complex got = hz_phi_prod_integral_osc(-4.0, 4.0, cases[i][0], cases[i][1], omegas[j]);
      double complex want = gl5_prod_osc(-4.0, 4.0, cases[i][0], cases[i][1], omegas[j]);
      check_close(cabs(got - want), 0.0, 1e-12 * cabs(want) + 5e-17, "osc at omega*h ~ 1e3");
    }
  }
}

static void test_osc_symmetry(void) {
  /* A configuration symmetric about x=0 has an even integrand, so the sin part
   * cancels exactly and the result must be real. Independent of the comparator:
   * it tests the moment parity that osc_moments does not hardcode. */
  hz_phi_factor f = {1.0, 0.0, 0}, fdd = {1.0, 0.0, 2};
  static const double om[4] = {0.7, 3.0, 9.0, 300.0};
  for (int j = 0; j < 4; j++) {
    double complex a = hz_phi_prod_integral_osc(-2.0, 2.0, f, f, om[j]);
    double complex b = hz_phi_prod_integral_osc(-2.0, 2.0, fdd, f, om[j]);
    check_close(cimag(a), 0.0, 1e-14, "symmetric <phi,phi> osc is real");
    check_close(cimag(b), 0.0, 1e-14, "symmetric <phi'',phi> osc is real");
  }
}

/* --- M9b/proj: single-factor oscillatory integral ---------------------------- */

/* Long double witness for Int f(x) e^{i omega x} dx, brute subdivision. */
static double complex gll_single_osc(double a, double b, hz_phi_factor f, double omega) {
  double pts[10];
  int np = split_knots(a, b, f, f, pts);
  long double complex total = 0.0L;
  for (int p = 0; p + 1 < np; p++) {
    long double lo = (long double)pts[p], hi = (long double)pts[p + 1];
    if (hi <= lo) continue;
    long ns = (long)(10.0 * (pts[p + 1] - pts[p]) * fabs(omega)) + 1;
    long double sub = (hi - lo) / (long double)ns;
    for (long k = 0; k < ns; k++) {
      long double half = 0.5L * sub, mid = lo + ((long double)k + 0.5L) * sub;
      long double complex s = 0.0L;
      for (int g = 0; g < 5; g++) {
        long double x = mid + half * GLL_X[g];
        long double ph = (long double)omega * x;
        s += GLL_W[g] * lfeval(f, x) * (cosl(ph) + I * sinl(ph));
      }
      total += s * half;
    }
  }
  return CMPLX((double)creall(total), (double)cimagl(total));
}

static void test_single_osc(void) {
  static const struct {
    hz_phi_factor f;
    double a, b;
  } cases[] = {
      {{1.0, 0.0, 0}, -3.0, 3.0},   {{1.0, 0.0, 1}, -3.0, 3.0}, {{1.0, 0.0, 2}, -3.0, 3.0},
      {{4.0, 1.0, 0}, -10.0, 30.0}, {{0.25, 7.0, 1}, 0.0, 4.0}, {{16.0, -1.0, 0}, -40.0, 40.0},
      {{2.0, 0.5, 2}, -6.0, 6.0},
  };
  int nc = (int)(sizeof(cases) / sizeof(cases[0]));
  static const double om[6] = {0.4, 2.0, 3.9, 4.1, 30.0, 400.0};
  for (int i = 0; i < nc; i++)
    for (int j = 0; j < 6; j++) {
      double complex got = hz_phi_integral_osc(cases[i].a, cases[i].b, cases[i].f, om[j]);
      double complex want = gll_single_osc(cases[i].a, cases[i].b, cases[i].f, om[j]);
      double scale = cabs(want) > 1.0 ? cabs(want) : 1.0;
      check_close(cabs(got - want), 0.0, 1e-11 * scale + 5e-17, "single-factor osc vs witness");
    }

  /* omega = 0 closed forms: Int phi = 4 over its support (partition of unity
   * has the same constant), Int phi' = 0 by oddness, Int phi'' = 0 because
   * phi' vanishes at both ends. Scaling by h multiplies by h^(1-deriv). */
  hz_phi_factor p = {1.0, 0.0, 0}, p1 = {1.0, 0.0, 1}, p2 = {1.0, 0.0, 2};
  check_close(creal(hz_phi_integral_osc(-2.0, 2.0, p, 0.0)), 4.0, 1e-14, "Int phi = 4");
  check_close(creal(hz_phi_integral_osc(-2.0, 2.0, p1, 0.0)), 0.0, 1e-14, "Int phi' = 0");
  check_close(creal(hz_phi_integral_osc(-2.0, 2.0, p2, 0.0)), 0.0, 1e-14, "Int phi'' = 0");
  hz_phi_factor ph8 = {8.0, 0.0, 0};
  check_close(creal(hz_phi_integral_osc(-16.0, 16.0, ph8, 0.0)), 32.0, 1e-13, "Int phi(x/8) = 32");

  /* Scaling identity: Int phi(x/h - n) e^{i w x} dx == h e^{i w h n} F(w h),
   * where F is the same integral at h=1, n=0. Ties different levels together
   * without going through the witness. */
  double h = 4.0, n = 3.0, w = 0.7;
  hz_phi_factor fs = {h, n, 0}, f0 = {1.0, 0.0, 0};
  double complex lhs = hz_phi_integral_osc(-100.0, 100.0, fs, w);
  double complex rhs =
      h * CMPLX(cos(w * h * n), sin(w * h * n)) * hz_phi_integral_osc(-2.0, 2.0, f0, w * h);
  check_close(cabs(lhs - rhs), 0.0, 1e-13 * cabs(rhs), "single-factor scaling identity");

  /* outside the support: exact zero */
  hz_phi_factor far = {1.0, 20.0, 0};
  g_total++;
  double complex v = hz_phi_integral_osc(-5.0, 5.0, far, 1.0);
  if (creal(v) < 0.0 || creal(v) > 0.0 || cimag(v) < 0.0 || cimag(v) > 0.0) {
    g_fail++;
    printf("FAIL: single-factor osc outside support must be exact 0\n");
  }
}

static void test_osc_disjoint(void) {
  hz_phi_factor f1 = {1.0, 0.0, 0};
  hz_phi_factor f2 = {1.0, 4.0, 0}; /* supports [-2,2] and [2,6]: touch at a point */
  g_total++;
  double complex v = hz_phi_prod_integral_osc(-10.0, 10.0, f1, f2, 3.0);
  if (creal(v) < 0.0 || creal(v) > 0.0 || cimag(v) < 0.0 || cimag(v) > 0.0) {
    g_fail++;
    printf("FAIL: osc disjoint supports must give exact 0, got (%.17g,%.17g)\n", creal(v),
           cimag(v));
  }
}

int main(void) {
  test_partition_of_unity();
  test_c1_continuity();
  test_two_scale();
  test_closed_forms();
  test_disjoint_supports();
  test_vs_quadrature();
  test_deriv1_vs_quadrature();
  test_osc_vs_quadrature();
  test_osc_zero_omega();
  test_osc_seam();
  test_osc_large_omega();
  test_osc_symmetry();
  test_single_osc();
  test_osc_disjoint();
  printf("test_phi: %d/%d passed\n", g_total - g_fail, g_total);
  return g_fail > 0 ? 1 : 0;
}
