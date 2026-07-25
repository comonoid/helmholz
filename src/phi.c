#include "phi.h"
#include <math.h>

double hz_phi(double t) {
  double a = fabs(t);
  if (a >= 2.0) return 0.0;
  if (a <= 1.0) return 2.0 - t * t;
  double d = a - 2.0;
  return d * d;
}

double hz_phi_d1(double t) {
  double a = fabs(t);
  if (a >= 2.0) return 0.0;
  if (a <= 1.0) return -2.0 * t;
  double s = t < 0.0 ? -1.0 : 1.0;
  return 2.0 * s * (a - 2.0);
}

double hz_phi_d2(double t) {
  double a = fabs(t);
  if (a >= 2.0) return 0.0;
  return a <= 1.0 ? -2.0 : 2.0;
}

/* Coefficients of phi on the piece containing t: phi(t) = c0 + c1*t + c2*t^2.
 * t must not sit exactly on a piece boundary (callers pass piece midpoints). */
static void piece_coef(double t, double c[3]) {
  if (t <= -2.0 || t >= 2.0) {
    c[0] = 0.0;
    c[1] = 0.0;
    c[2] = 0.0;
  } else if (t < -1.0) { /* (t+2)^2 */
    c[0] = 4.0;
    c[1] = 4.0;
    c[2] = 1.0;
  } else if (t <= 1.0) { /* 2 - t^2 */
    c[0] = 2.0;
    c[1] = 0.0;
    c[2] = -1.0;
  } else { /* (t-2)^2 */
    c[0] = 4.0;
    c[1] = -4.0;
    c[2] = 1.0;
  }
}

/* Polynomial of one factor in the local coordinate s = x - m (m = piece
 * midpoint): q[0] + q[1]*s + q[2]*s^2. Local coordinates keep the expansion
 * well-conditioned for narrow pieces far from the origin. */
static void local_poly(hz_phi_factor f, double m, double q[3]) {
  double u = 1.0 / f.h;
  double tm = m * u - f.n;
  double c[3];
  piece_coef(tm, c);
  if (f.deriv == 2) {
    q[0] = 2.0 * c[2] * u * u;
    q[1] = 0.0;
    q[2] = 0.0;
  } else if (f.deriv == 1) {
    /* d/ds of the deriv-0 polynomial below, times the chain factor u */
    q[0] = (c[1] + 2.0 * c[2] * tm) * u;
    q[1] = 2.0 * c[2] * u * u;
    q[2] = 0.0;
  } else {
    q[0] = c[0] + (c[1] + c[2] * tm) * tm;
    q[1] = (c[1] + 2.0 * c[2] * tm) * u;
    q[2] = c[2] * u * u;
  }
}

static void dsort(double *v, int n) {
  for (int i = 1; i < n; i++) {
    double x = v[i];
    int j = i - 1;
    while (j >= 0 && v[j] > x) {
      v[j + 1] = v[j];
      j--;
    }
    v[j + 1] = x;
  }
}

/* Clip [a,b] to both supports and split at every interior knot of both factors
 * (piece boundaries of phi at t in {-2,-1,1,2}; t=0 is not a knot, [-1,1] is
 * one parabola). Writes the sorted breakpoints to pts and returns their count,
 * 0 when the supports do not meet [a,b]. pts must hold 10 doubles: 2 ends plus
 * at most 4 interior knots per factor. */
static int clip_split(double a, double b, const hz_phi_factor fs[2], double pts[10]) {
  double lo = a, hi = b;
  for (int i = 0; i < 2; i++) {
    double s0 = fs[i].h * (fs[i].n - 2.0);
    double s1 = fs[i].h * (fs[i].n + 2.0);
    if (s0 > lo) lo = s0;
    if (s1 < hi) hi = s1;
  }
  if (lo >= hi) return 0;

  int np = 0;
  pts[np++] = lo;
  pts[np++] = hi;
  static const double knot_t[4] = {-2.0, -1.0, 1.0, 2.0};
  for (int i = 0; i < 2; i++) {
    for (int k = 0; k < 4; k++) {
      double x = fs[i].h * (fs[i].n + knot_t[k]);
      if (x > lo && x < hi) pts[np++] = x;
    }
  }
  dsort(pts, np);
  return np;
}

double hz_phi_prod_integral(double a, double b, hz_phi_factor f1, hz_phi_factor f2) {
  const hz_phi_factor fs[2] = {f1, f2};
  double pts[10];
  int np = clip_split(a, b, fs, pts);
  if (np == 0) return 0.0;

  /* On each piece the integrand is a degree<=4 polynomial; integrate in the
   * midpoint-centered coordinate, where odd powers vanish. */
  double total = 0.0;
  for (int p = 0; p + 1 < np; p++) {
    double w = 0.5 * (pts[p + 1] - pts[p]);
    if (w <= 0.0) continue;
    double m = 0.5 * (pts[p] + pts[p + 1]);
    double q1[3], q2[3];
    local_poly(fs[0], m, q1);
    local_poly(fs[1], m, q2);
    double r0 = q1[0] * q2[0];
    double r2 = q1[0] * q2[2] + q1[1] * q2[1] + q1[2] * q2[0];
    double r4 = q1[2] * q2[2];
    double w2 = w * w;
    total += 2.0 * w * (r0 + w2 * (r2 / 3.0 + w2 * r4 / 5.0));
  }
  return total;
}

/* Series/closed-form seam for the oscillatory moments below. The Taylor series
 * terms peak at m ~ |a| with size |a|^m/m!, so the cancellation it costs is
 * bounded by that peak divided by the O(1) result: under one decimal digit at
 * |a| = 4 (peak 4^4/4! = 10.7), growing fast beyond. The closed forms are the
 * mirror image -- they divide by powers of a, so they are clean for large |a|
 * and lose digits for small. At |a| = 4 both branches are accurate, and
 * test_phi checks they agree across the seam. */
static const double OSC_SERIES_MAX = 4.0;

/* mom[n] = integral over [-1,1] of u^n * exp(i*a*u) du, for n = 0..4.
 * Purely real for even n, purely imaginary for odd n (parity of the integrand);
 * the code does not exploit that, so a violation shows up as a test failure
 * rather than being silently forced. */
static void osc_moments(double complex a, double complex mom[5]) {
  if (cabs(a) <= OSC_SERIES_MAX) {
    /* Expand exp(i*a*u) and integrate term by term: the u^(n+m) moment over
     * [-1,1] is 2/(n+m+1) for n+m even and vanishes otherwise. */
    for (int n = 0; n < 5; n++)
      mom[n] = 0.0;
    double complex coef = 1.0; /* (i*a)^m / m! */
    for (int m = 0; m < 48; m++) {
      if (m > 0) coef *= CMPLX(0.0, 1.0) * a / (double)m;
      for (int n = m & 1; n < 5; n += 2)
        mom[n] += 2.0 * coef / (double)(n + m + 1);
      /* past the peak the terms fall factorially; 1e-19 is below the double
       * epsilon of the O(1) result, so the tail cannot move it */
      if ((double)m > cabs(a) && cabs(coef) < 1e-19) break;
    }
    return;
  }
  /* Repeated integration by parts, exact. */
  double complex s = csin(a), c = ccos(a);
  double complex a1 = 1.0 / a, a2 = a1 * a1, a3 = a2 * a1, a4 = a3 * a1, a5 = a4 * a1;
  mom[0] = 2.0 * s * a1;
  mom[1] = CMPLX(0.0, 2.0) * (s * a2 - c * a1);
  mom[2] = 2.0 * (s * a1 + 2.0 * c * a2 - 2.0 * s * a3);
  mom[3] = CMPLX(0.0, 2.0) * (-c * a1 + 3.0 * s * a2 + 6.0 * c * a3 - 6.0 * s * a4);
  mom[4] = 2.0 * (s * a1 + 4.0 * c * a2 - 12.0 * s * a3 - 24.0 * c * a4 + 24.0 * s * a5);
}

double complex hz_phi_integral_osc(double a, double b, hz_phi_factor f, double complex omega) {
  /* clip_split takes a pair; passing the same factor twice gives exactly this
   * factor's support and knots, which is what a single-factor integral needs */
  const hz_phi_factor fs[2] = {f, f};
  double pts[10];
  int np = clip_split(a, b, fs, pts);
  if (np == 0) return 0.0;

  double complex total = 0.0;
  for (int p = 0; p + 1 < np; p++) {
    double w = 0.5 * (pts[p + 1] - pts[p]);
    if (w <= 0.0) continue;
    double m = 0.5 * (pts[p] + pts[p + 1]);
    double q[3];
    local_poly(f, m, q);
    double complex mom[5];
    osc_moments(omega * w, mom);
    double complex acc = 0.0;
    double wp = w;
    for (int n = 0; n < 3; n++) {
      acc += q[n] * mom[n] * wp;
      wp *= w;
    }
    total += cexp(CMPLX(0.0, 1.0) * omega * m) * acc;
  }
  return total;
}

double complex hz_phi_prod_integral_osc(double a, double b, hz_phi_factor f1, hz_phi_factor f2,
                                        double complex omega) {
  const hz_phi_factor fs[2] = {f1, f2};
  double pts[10];
  int np = clip_split(a, b, fs, pts);
  if (np == 0) return 0.0;

  double complex total = 0.0;
  for (int p = 0; p + 1 < np; p++) {
    double w = 0.5 * (pts[p + 1] - pts[p]);
    if (w <= 0.0) continue;
    double m = 0.5 * (pts[p] + pts[p + 1]);
    double q1[3], q2[3];
    local_poly(fs[0], m, q1);
    local_poly(fs[1], m, q2);
    /* product of the two local quadratics, degree <= 4 in s = x - m */
    double r[5];
    r[0] = q1[0] * q2[0];
    r[1] = q1[0] * q2[1] + q1[1] * q2[0];
    r[2] = q1[0] * q2[2] + q1[1] * q2[1] + q1[2] * q2[0];
    r[3] = q1[1] * q2[2] + q1[2] * q2[1];
    r[4] = q1[2] * q2[2];

    /* integral of s^n * exp(i*omega*s) over [-w,w] is w^(n+1)*mom[n](omega*w);
     * the piece midpoint contributes the exp(i*omega*m) carrier phase */
    double complex mom[5];
    osc_moments(omega * w, mom);
    double complex acc = 0.0;
    double wp = w;
    for (int n = 0; n < 5; n++) {
      acc += r[n] * mom[n] * wp;
      wp *= w;
    }
    total += cexp(CMPLX(0.0, 1.0) * omega * m) * acc;
  }
  return total;
}
