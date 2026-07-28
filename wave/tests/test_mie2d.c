/* Gate D2-3 of the vertical slice: the exact scattering reference.
 *
 * ARTEFACT 12 COST THIS PROJECT A SESSION — a reference taken for the wrong
 * source, whose error stayed politely constant across every domain size. The
 * defence here is that no single statement is trusted:
 *   [1] the closed form and the ODE integrator are made to compute the SAME
 *       quantity by different routes;
 *   [2] the ODE is checked by step halving (must be 4th order);
 *   [3] UNITARITY, |1 + 2 c_m| = 1, which needs no second reference at all: for
 *       a lossless radial scatterer each partial wave is a pure phase shift, so
 *       this one identity validates the Bessel functions, the matching and the
 *       ODE at once;
 *   [4] the first Born approximation as an INDEPENDENT reference. Note the
 *       criterion: not agreement to a fixed tolerance (Born is wrong by
 *       O(contrast) by construction — the threshold in the audit was simply
 *       wrong), but agreement whose error falls LINEARLY with the contrast. A
 *       rate is a stronger statement than a magnitude. */
#include "bessel.h"
#include "mie2d.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>

enum { MM = 200 };
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

int main(void) {
  printf("test_mie2d: exact 2D scattering by a radial scatterer\n");
  double k0 = 2.0 * M_PI, a = 1.0; /* lambda = 1, ka = 6.28 */
  double complex c[MM + 1], c2[MM + 1];

  /* --- [1] closed form vs the ODE on the SAME sharp profile ---------------- */
  {
    hz_mie sc = {HZ_MIE_SHARP, a, k0, 1.5, 24, 8000};
    hz_mie so = {HZ_MIE_SHARP_ODE, a, k0, 1.5, 24, 8000};
    hz_mie_coeffs(&sc, c);
    hz_mie_coeffs(&so, c2);
    double worst = 0.0;
    for (int m = 0; m <= sc.mmax; m++) {
      double d = cabs(c[m] - c2[m]);
      if (d > worst) worst = d;
    }
    check("[1] closed form == ODE on the sharp profile", worst < 1e-9, worst);
  }

  /* --- [2] the ODE is 4th order: halving the step cuts the error by 16 ----- */
  {
    hz_mie s1 = {HZ_MIE_GRADED, a, k0, 0.3225, 24, 500};
    hz_mie s2 = {HZ_MIE_GRADED, a, k0, 0.3225, 24, 1000};
    hz_mie s4 = {HZ_MIE_GRADED, a, k0, 0.3225, 24, 2000};
    double complex ca[MM + 1], cb[MM + 1], cc[MM + 1];
    hz_mie_coeffs(&s1, ca);
    hz_mie_coeffs(&s2, cb);
    hz_mie_coeffs(&s4, cc);
    double e1 = 0.0, e2 = 0.0;
    for (int m = 0; m <= s1.mmax; m++) {
      double d1 = cabs(ca[m] - cc[m]), d2 = cabs(cb[m] - cc[m]);
      if (d1 > e1) e1 = d1;
      if (d2 > e2) e2 = d2;
    }
    double rate = (e2 > 0.0) ? e1 / e2 : 1e9;
    check("[2] ODE convergence rate (want ~16, 4th order)", rate > 8.0, rate);
  }

  /* --- [3] unitarity: |1 + 2 c_m| = 1 for a lossless scatterer ------------- */
  {
    hz_mie sh = {HZ_MIE_SHARP, a, k0, 1.5, 40, 8000};
    hz_mie gr = {HZ_MIE_GRADED, a, k0, 0.3225, 40, 8000};
    hz_mie_coeffs(&sh, c);
    hz_mie_coeffs(&gr, c2);
    double w1 = 0.0, w2 = 0.0;
    for (int m = 0; m <= sh.mmax; m++) {
      double d1 = fabs(cabs(1.0 + 2.0 * c[m]) - 1.0);
      double d2 = fabs(cabs(1.0 + 2.0 * c2[m]) - 1.0);
      if (d1 > w1) w1 = d1;
      if (d2 > w2) w2 = d2;
    }
    check("[3] unitarity |1+2c_m| = 1, sharp", w1 < 1e-10, w1);
    check("[3] unitarity |1+2c_m| = 1, graded", w2 < 1e-9, w2);
  }

  /* --- [4] the harmonic tail must have decayed at m = mmax ---------------- */
  {
    hz_mie sh = {HZ_MIE_SHARP, a, k0, 1.5, 40, 8000};
    hz_mie_coeffs(&sh, c);
    double mx = 0.0;
    for (int m = 0; m <= sh.mmax; m++)
      if (cabs(c[m]) > mx) mx = cabs(c[m]);
    double tail = cabs(c[sh.mmax]) / mx;
    check("[4] tail |c_M|/max|c| at M = ka + 34", tail < 1e-6, tail);
  }

  /* --- [5] zero contrast scatters nothing -------------------------------- */
  {
    hz_mie sh = {HZ_MIE_SHARP, a, k0, 1.0, 24, 8000};
    hz_mie gr = {HZ_MIE_GRADED, a, k0, 0.0, 24, 8000};
    hz_mie_coeffs(&sh, c);
    hz_mie_coeffs(&gr, c2);
    double w = 0.0;
    for (int m = 0; m <= 24; m++) {
      if (cabs(c[m]) > w) w = cabs(c[m]);
      if (cabs(c2[m]) > w) w = cabs(c2[m]);
    }
    check("[5] zero contrast => c_m = 0", w < 1e-12, w);
  }

  /* --- [6] Born: the DIFFERENCE must fall linearly with the contrast ------ */
  {
    double prev = 0.0, rate = 0.0;
    printf("       contrast   |f_exact|      |f_exact - f_Born| / |f_exact|\n");
    for (int i = 0; i < 4; i++) {
      double d = 0.08 / pow(2.0, (double)i);
      hz_mie gr = {HZ_MIE_GRADED, a, k0, d, 40, 8000};
      hz_mie_coeffs(&gr, c);
      double num = 0.0, den = 0.0;
      for (int t = 0; t < 32; t++) {
        double th = 2.0 * M_PI * (double)t / 32.0;
        double complex fe = hz_mie_far(&gr, c, th), fb = hz_mie_born(&gr, th);
        num += cabs(fe - fb) * cabs(fe - fb);
        den += cabs(fe) * cabs(fe);
      }
      double rel = sqrt(num / den);
      printf("       %8.4f   %10.4e   %12.4e\n", d, sqrt(den / 32.0), rel);
      if (i > 0) rate = prev / rel;
      prev = rel;
    }
    check("[6] Born error halves with the contrast (rate ~2)", rate > 1.7 && rate < 2.4, rate);
  }

  /* --- [7] the far field agrees with the near field pushed out ------------- */
  {
    hz_mie sh = {HZ_MIE_SHARP, a, k0, 1.5, 40, 8000};
    hz_mie_coeffs(&sh, c);
    double r = 4000.0 * a, th = 0.7;
    double complex us = hz_mie_scat(&sh, c, r * cos(th), r * sin(th));
    double complex asym = sqrt(2.0 / (M_PI * k0 * r)) * cexp(CMPLX(0.0, k0 * r - M_PI / 4.0)) *
                          hz_mie_far(&sh, c, th);
    double e = cabs(us - asym) / cabs(us);
    check("[7] near field -> far-field pattern at kr = 2.5e4", e < 1e-3, e);
  }

  printf("test_mie2d: %d/%d passed\n", pass_count, pass_count + fail_count);
  return fail_count == 0 ? 0 : 1;
}
