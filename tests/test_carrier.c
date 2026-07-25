/* Acceptance for the directional-carrier refactor (PLAN queue step 2).
 *
 * (a) SPAN UNCHANGED. The directional form {phi e^{+ikx}, phi e^{-ikx}} and the
 *     old quadrature form {phi cos, phi sin} must span the SAME space, so the
 *     best L2 approximation of any field is identical. This is the criterion
 *     that a "pure change of coordinates" actually is one.
 * (b) DIRECTIONS ARE SEPARABLE. This is the whole reason for the refactor: a
 *     pure e^{+ikx} projected on the directional basis must put ~nothing on the
 *     -k functions, and dropping them must not change the result. With cos/sin
 *     that statement cannot even be expressed.
 * (c) OPERATOR CANCELLATION. With the carrier matched to the medium, the k^2
 *     term must cancel analytically: the entry must equal the pure
 *     phi''/phi' part with no medium contribution.
 * (d) CUT. Functions on opposite sides of a cut have disjoint supports. */
#include "carrier.h"
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int g_fail = 0;
static int g_total = 0;

static void check(int ok, const char *what, double got) {
  g_total++;
  if (!ok) {
    g_fail++;
    printf("FAIL: %s (got %.6e)\n", what, got);
  }
}

enum { NB = 64, NSAMP = 4000 };
static const double LAM = 16.0;
static const double RCOND = 1e-12;

/* Least-squares fit of a sampled field by a set of sampled basis functions;
 * returns the relative residual and writes the coefficients. Deliberately
 * sampled and dumb: this test must not share machinery with what it checks. */
static double fit(const double complex *B, int nb, const double complex *u, int ns,
                  double complex *coef) {
  double complex *A = calloc((size_t)ns * (size_t)nb, sizeof(double complex));
  double complex *b = calloc((size_t)(ns > nb ? ns : nb), sizeof(double complex));
  double *sv = calloc((size_t)nb, sizeof(double));
  if (!A || !b || !sv) {
    free(A);
    free(b);
    free(sv);
    return -1.0;
  }
  for (int p = 0; p < ns; p++) {
    for (int i = 0; i < nb; i++)
      A[(size_t)p * (size_t)nb + (size_t)i] = B[(size_t)i * (size_t)ns + (size_t)p];
    b[p] = u[p];
  }
  lapack_int rank = 0;
  lapack_int info = LAPACKE_zgelsd(LAPACK_ROW_MAJOR, ns, nb, 1, A, nb, b, 1, sv, RCOND, &rank);
  double num = 0.0, den = 0.0;
  if (info == 0) {
    for (int p = 0; p < ns; p++) {
      double complex v = 0.0;
      for (int i = 0; i < nb; i++)
        v += b[i] * B[(size_t)i * (size_t)ns + (size_t)p];
      num += creal((u[p] - v) * conj(u[p] - v));
      den += creal(u[p] * conj(u[p]));
    }
    for (int i = 0; i < nb; i++)
      coef[i] = b[i];
  }
  free(A);
  free(b);
  free(sv);
  return info == 0 ? sqrt(num / den) : -1.0;
}

int main(void) {
  double k = 2.0 * M_PI / LAM;
  double W = 4.0 * LAM; /* deliberately coarse: 4 wavelengths per element */
  double x0 = -2.0 * W, x1 = 2.0 * W;

  static double complex Bdir[NB][NSAMP], Bquad[NB][NSAMP], u[NSAMP];
  static double complex cd[NB], cq[NB];
  double xs[NSAMP];
  for (int p = 0; p < NSAMP; p++)
    xs[p] = x0 + (x1 - x0) * ((double)p + 0.5) / (double)NSAMP;

  /* --- (a) span unchanged --------------------------------------------- */
  int nd = 0, nq = 0;
  for (int n = -4; n <= 4; n++) {
    for (int s = 0; s < 2; s++) {
      hz_carrier c = {W, n, s ? -k : k, -1};
      for (int p = 0; p < NSAMP; p++)
        Bdir[nd][p] = hz_carrier_val(c, xs[p], 0.0);
      nd++;
    }
    for (int q = 0; q < 2; q++) {
      for (int p = 0; p < NSAMP; p++) {
        double t = xs[p] / W - (double)n;
        double ph = fabs(t) < 2.0 ? hz_phi(t) : 0.0;
        Bquad[nq][p] = q ? ph * sin(k * xs[p]) : ph * cos(k * xs[p]);
      }
      nq++;
    }
  }
  /* a field that is NOT exactly in the span, so the two fits have something to
   * disagree about: a wave at 1.3k, deliberately mismatched */
  for (int p = 0; p < NSAMP; p++)
    u[p] = CMPLX(cos(1.3 * k * xs[p]), sin(1.3 * k * xs[p]));
  double ed = fit(&Bdir[0][0], nd, u, NSAMP, cd);
  double eq = fit(&Bquad[0][0], nq, u, NSAMP, cq);
  check(nd == nq, "same dimension", (double)(nd - nq));
  check(ed >= 0.0 && eq >= 0.0 && fabs(ed - eq) <= 1e-12 * (eq > 1.0 ? eq : 1.0),
        "span unchanged: same best-approximation error", fabs(ed - eq));

  /* --- (b) directions separable --------------------------------------- */
  for (int p = 0; p < NSAMP; p++)
    u[p] = CMPLX(cos(k * xs[p]), sin(k * xs[p])); /* pure +k */
  double efull = fit(&Bdir[0][0], nd, u, NSAMP, cd);
  double amp_plus = 0.0, amp_minus = 0.0;
  for (int i = 0; i < nd; i++) {
    double a = cabs(cd[i]);
    if (i % 2 == 0)
      amp_plus += a;
    else
      amp_minus += a;
  }
  check(amp_minus <= 1e-8 * amp_plus, "pure +k puts nothing on -k functions",
        amp_minus / (amp_plus > 0.0 ? amp_plus : 1.0));

  /* dropping the -k half must not degrade the fit at all */
  static double complex Bplus[NB][NSAMP];
  int np = 0;
  for (int i = 0; i < nd; i += 2) {
    for (int p = 0; p < NSAMP; p++)
      Bplus[np][p] = Bdir[i][p];
    np++;
  }
  double ehalf = fit(&Bplus[0][0], np, u, NSAMP, cd);
  check(fabs(ehalf - efull) <= 1e-10, "dropping -k costs nothing for a +k field",
        fabs(ehalf - efull));

  /* --- (c) operator cancellation -------------------------------------- */
  hz_medseg seg = {-1e9, 1e9, CMPLX(k * k, 0.0)}; /* medium matches the carrier */
  hz_carrier bi = {W, 0, k, -1}, bj = {W, 1, k, -1};
  double complex with_med = hz_carrier_entry(bi, bj, 1e9, 0.0, &seg, 1, NULL);
  double complex no_med = hz_carrier_entry(bi, bj, 1e9, 0.0, NULL, 0, NULL);
  check(cabs(with_med - no_med) <= 1e-12 * (cabs(no_med) > 1.0 ? cabs(no_med) : 1.0),
        "matched carrier: k^2 term cancels, medium contributes nothing", cabs(with_med - no_med));
  /* and a MISMATCHED medium must contribute */
  hz_medseg seg2 = {-1e9, 1e9, CMPLX(2.25 * k * k, 0.0)};
  double complex mismatched = hz_carrier_entry(bi, bj, 1e9, 0.0, &seg2, 1, NULL);
  check(cabs(mismatched - no_med) > 1e-6, "mismatched medium does contribute",
        cabs(mismatched - no_med));

  /* --- (d) cut ---------------------------------------------------------- */
  hz_carrier left = {W, 0, k, 0}, right = {W, 0, k, 1};
  double lo = 0.0, hi = 0.0;
  hz_carrier_overlap(left, right, 1e9, 0.0, &lo, &hi);
  check(lo >= hi, "opposite sides of a cut do not overlap", hi - lo);
  double complex e_cut = hz_carrier_entry(left, right, 1e9, 0.0, &seg, 1, NULL);
  check(!(cabs(e_cut) > 0.0), "entry across a cut is exactly zero", cabs(e_cut));

  printf("test_carrier: %d/%d passed\n", g_total - g_fail, g_total);
  return g_fail > 0 ? 1 : 0;
}
