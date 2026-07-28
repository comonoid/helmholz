#include "mie2d.h"
#include "bessel.h"
#include <math.h>

/* The INTERIOR formula, continued analytically to r = a and beyond. The ODE
 * below needs exactly this: its last RK4 stage samples r = a, and if that
 * returned the exterior value the final step would be wrong by O(h) times the
 * contrast. Measured before the split: the sharp profile disagreed with its own
 * closed form by 4.768e-04, and the number did not move across three different
 * integrators — constancy under the parameter you are changing has meant
 * systematics every time in this project (artefact 15). */
static double k2_interior(const hz_mie *s, double r) {
  double k02 = s->k0 * s->k0;
  if (s->prof != HZ_MIE_GRADED) return k02 * s->param * s->param;
  double t = 1.0 - (r / s->a) * (r / s->a);
  return k02 * (1.0 + s->param * t * t * t);
}

double hz_mie_k2(const hz_mie *s, double r) {
  if (r >= s->a) return s->k0 * s->k0;
  return k2_interior(s, r);
}

/* Interior logarithmic derivative L(a) = R'(a)/R(a) of the regular solution of
 *     R'' + R'/r + (k^2(r) - m^2/r^2) R = 0,
 * returned as a FRACTION Ln/Ld so that a zero of R (which happens inside the
 * scatterer for every m below ka) is not a division by zero.
 *
 * SUBSTITUTION R = r^m u, WITHOUT WHICH THIS DOES NOT WORK. In R the m^2/r^2
 * term forces a step much smaller than r/m, and the integration starts at r = h
 * where that is violated outright — the first run showed it as a step-INDEPENDENT
 * relative error of 5e-4 (8000 and 32000 steps disagreed by as much as either
 * disagreed with the closed form). Substituting kills the singular term exactly:
 *     u'' + ((2m+1)/r) u' + k^2 u = 0,
 * whose regular solution has u(0) = 1, u'(0) = 0 and u' ~ r, so the surviving
 * (2m+1)u'/r is finite at the origin and nothing is stiff. As a bonus u stays
 * O(1), so no renormalisation is needed either. */
static void interior_logderiv(const hz_mie *s, int m, double *Ln, double *Ld) {
  if (s->prof == HZ_MIE_SHARP) {
    double k1 = s->k0 * s->param;
    double j, y, dj, dy;
    hz_bessel_jy(m, k1 * s->a, &j, &y, &dj, &dy);
    *Ln = k1 * dj;
    *Ld = j;
    return;
  }
  int n = s->nr;
  double a = s->a, h = a / (double)n;
  double dm = (double)m, kap2 = k2_interior(s, 0.0);
  double r = h;
  double u = 1.0 - kap2 * h * h / (4.0 * (dm + 1.0));
  double v = -kap2 * h / (2.0 * (dm + 1.0));
  for (int i = 1; i < n; i++) {
    double k1u, k1v, k2u, k2v, k3u, k3v, k4u, k4v, rr, uu, vv;
    rr = r;
    uu = u;
    vv = v;
    k1u = vv;
    k1v = -(2.0 * dm + 1.0) / rr * vv - k2_interior(s, rr) * uu;
    rr = r + 0.5 * h;
    uu = u + 0.5 * h * k1u;
    vv = v + 0.5 * h * k1v;
    k2u = vv;
    k2v = -(2.0 * dm + 1.0) / rr * vv - k2_interior(s, rr) * uu;
    uu = u + 0.5 * h * k2u;
    vv = v + 0.5 * h * k2v;
    k3u = vv;
    k3v = -(2.0 * dm + 1.0) / rr * vv - k2_interior(s, rr) * uu;
    rr = r + h;
    uu = u + h * k3u;
    vv = v + h * k3v;
    k4u = vv;
    k4v = -(2.0 * dm + 1.0) / rr * vv - k2_interior(s, rr) * uu;
    u += h / 6.0 * (k1u + 2.0 * k2u + 2.0 * k3u + k4u);
    v += h / 6.0 * (k1v + 2.0 * k2v + 2.0 * k3v + k4v);
    r += h;
  }
  *Ln = dm / a * u + v; /* L = m/a + u'/u */
  *Ld = u;
}

int hz_mie_coeffs(const hz_mie *s, double complex *c) {
  double ka = s->k0 * s->a;
  for (int m = 0; m <= s->mmax; m++) {
    double Ln, Ld;
    interior_logderiv(s, m, &Ln, &Ld);
    double j, y, dj, dy;
    hz_bessel_jy(m, ka, &j, &y, &dj, &dy);
    double complex H = CMPLX(j, y), dH = CMPLX(dj, dy);
    /* multiplied through by Ld so a zero of the interior solution is harmless */
    double complex num = Ln * j - s->k0 * dj * Ld;
    double complex den = s->k0 * dH * Ld - Ln * H;
    c[m] = num / den;
  }
  return 0;
}

double complex hz_mie_far(const hz_mie *s, const double complex *c, double theta) {
  double complex f = c[0];
  for (int m = 1; m <= s->mmax; m++)
    f += 2.0 * c[m] * cos((double)m * theta);
  return f;
}

double complex hz_mie_scat(const hz_mie *s, const double complex *c, double x, double y) {
  double r = sqrt(x * x + y * y), th = atan2(y, x);
  if (r < s->a) return 0.0; /* exterior only; the interior needs R_m(r) */
  double complex u = 0.0, ipow = 1.0;
  for (int m = 0; m <= s->mmax; m++) {
    double jj, yy, dj, dy;
    hz_bessel_jy(m, s->k0 * r, &jj, &yy, &dj, &dy);
    double complex H = CMPLX(jj, yy);
    double w = (m == 0) ? 1.0 : 2.0 * cos((double)m * th);
    u += ipow * c[m] * H * w;
    ipow *= CMPLX(0.0, 1.0);
  }
  return u;
}

double complex hz_mie_born(const hz_mie *s, double theta) {
  /* f_Born = (i/4) * FT of the contrast at q = 2 k sin(theta/2), which for a
   * radial profile is a Hankel transform: 2 pi Int V(r) J0(q r) r dr. */
  double q = 2.0 * s->k0 * fabs(sin(0.5 * theta));
  int n = 4000;
  double h = s->a / (double)n, acc = 0.0;
  for (int i = 0; i < n; i++) {
    double r = h * ((double)i + 0.5);
    double V = hz_mie_k2(s, r) - s->k0 * s->k0;
    acc += V * j0(q * r) * r * h;
  }
  return CMPLX(0.0, 0.25) * 2.0 * M_PI * acc;
}
