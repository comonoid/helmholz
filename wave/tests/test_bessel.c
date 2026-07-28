/* Gate D2-0 of the vertical slice: the DtN symbol H'_m/H_m.
 *
 * SLICE_PLAN doubt 4 said the harmonic order was written into the plan as one
 * line, as if trivial. The audit downgraded the danger (no overflow below
 * m ~ 3.6x) but kept the gate, because the whole termination is this one number
 * and a wrong Y_m would look like "the basis is bad near the boundary".
 * Three independent statements, none of which shares a code path with the
 * others: the Wronskian (an identity between J, Y and their derivatives), the
 * three-term recurrence (an identity in m alone), and the evanescent limit
 * (an asymptotic statement about large m). */
#include "bessel.h"
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

int main(void) {
  printf("test_bessel: the DtN symbol on a circle\n");
  static const double X[5] = {1.0, 10.0, 31.4, 200.0, 1000.0};

  /* --- [1] Wronskian over the range the slice actually uses --------------- */
  for (int ix = 0; ix < 5; ix++) {
    double x = X[ix], worst = 0.0;
    int mworst = 0, mmax = (int)(3.0 * x) + 10;
    for (int m = 0; m <= mmax; m++) {
      double e = hz_bessel_wronskian_err(m, x);
      if (!isfinite(e)) break; /* Y_m overflowed: [3] pins where */
      if (e > worst) {
        worst = e;
        mworst = m;
      }
    }
    char nm[64];
    snprintf(nm, sizeof nm, "[1] Wronskian x=%.1f, worst m=%d", x, mworst);
    check(nm, worst < 1e-10, worst);
  }

  /* --- [2] three-term recurrence, an identity in m alone ------------------ */
  for (int ix = 0; ix < 5; ix++) {
    double x = X[ix], worst = 0.0;
    for (int m = 1; m <= (int)x + 20; m++) {
      double jm = jn(m, x), jp = jn(m + 1, x), jmm = jn(m - 1, x);
      double ym = yn(m, x), yp = yn(m + 1, x), ymm = yn(m - 1, x);
      double sj = fabs(jmm + jp - 2.0 * (double)m / x * jm);
      double sy = fabs(ymm + yp - 2.0 * (double)m / x * ym);
      double scale = fabs(jm) + fabs(jp) + fabs(jmm);
      double scy = fabs(ym) + fabs(yp) + fabs(ymm);
      if (scale > 0.0 && sj / scale > worst) worst = sj / scale;
      if (scy > 0.0 && sy / scy > worst) worst = sy / scy;
    }
    char nm[64];
    snprintf(nm, sizeof nm, "[2] recurrence x=%.1f", x);
    check(nm, worst < 1e-10, worst);
  }

  /* --- [3] where does Y_m overflow, and is the fallback the right limit? --- */
  {
    double x = 200.0;
    int mov = 0;
    for (int m = 0; m < 4000; m++)
      if (!isfinite(yn(m, x))) {
        mov = m;
        break;
      }
    printf("  note x=200: Y_m first overflows at m=%d (m/x = %.2f)\n", mov, (double)mov / x);
    check("[3] overflow only well above m = 3x", (double)mov > 3.0 * x, (double)mov / x);
    /* below the overflow the evanescent limit must already be approached */
    double complex r = hz_hankel_ratio(500, x);
    double relerr = fabs(creal(r) + 500.0 / x) / (500.0 / x);
    check("[3] H'/H -> -m/x for m >> x (m=500, x=200)", relerr < 0.35, relerr);
  }

  /* --- [4] the symbol is a RADIATION condition: Re part of the flux ------- */
  /* For a purely outgoing wave the radial flux must be OUT, i.e. the symbol has
   * a positive imaginary part for propagating orders (m < x). This is the sign
   * convention the DtN rows depend on; getting it backwards is the single most
   * likely error in the termination (slice NC-B). */
  {
    double x = 31.4;
    int bad = 0;
    double worst = 1.0;
    for (int m = 0; m < (int)x; m++) {
      double complex r = hz_hankel_ratio(m, x);
      if (!(cimag(r) > 0.0)) bad++;
      if (cimag(r) < worst) worst = cimag(r);
    }
    check("[4] Im(H'/H) > 0 for every propagating order", bad == 0, worst);
    /* and for m >> x the symbol is real and negative (evanescent, no flux) */
    double complex re = hz_hankel_ratio(100, x);
    check("[4] evanescent order carries no flux", fabs(cimag(re)) < 1e-6 * fabs(creal(re)),
          fabs(cimag(re)) / fabs(creal(re)));
  }

  /* --- [5] large-argument asymptotics: H'/H -> i for x >> m --------------- */
  {
    double complex r = hz_hankel_ratio(0, 1e5);
    double e = cabs(r - CMPLX(0.0, 1.0));
    check("[5] H'/H -> i as x -> infinity (m=0, x=1e5)", e < 1e-4, e);
  }

  printf("test_bessel: %d/%d passed\n", pass_count, pass_count + fail_count);
  return fail_count == 0 ? 0 : 1;
}
