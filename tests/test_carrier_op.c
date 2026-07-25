/* Operator consistency for the carrier basis — exact, and with NO BOUNDARY.
 *
 * Every solver bench in this project has been spoiled by its domain
 * termination. This test avoids the problem entirely instead of fighting it:
 * in a homogeneous medium a plane wave is an EXACT solution of the homogeneous
 * equation, and by the partition of unity (sum of phi translates == 4) it lies
 * in the span of the carrier basis with CONSTANT coefficients c = 1/4. So the
 * discrete operator must ANNIHILATE it:
 *     sum_j A_ij c_j = 0   for every i whose support is interior.
 * That single identity checks the assembly, the analytic k^2 cancellation, the
 * partition of unity and the medium term at once, exactly, with no boundary,
 * no solve and no conditioning.
 *
 * It also has a negative control (a medium that does NOT match the carrier must
 * give a non-zero residual) — without one, a bug that returns zero everywhere
 * would pass. */
#include "carrier.h"
#include "phi.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>

static int g_fail = 0;
static int g_total = 0;

static void check(int ok, const char *what, double got) {
  g_total++;
  if (!ok) {
    g_fail++;
    printf("FAIL: %s (got %.6e)\n", what, got);
  }
}

enum { NB = 900 };
static const double LAM = 16.0;
static const double PU = 4.0; /* sum of phi translates */

/* Largest relative residual over elements whose support lies well inside the
 * domain. Relative to the row's own magnitude, so it is scale-free. */
/* The basis must hold BOTH directions, and the reason is the bilinear form.
 * <B_i, L B_j> carries e^{i(kx_i + kx_j)x}: with both functions at +k that is
 * e^{2ikx}, whose integral is exponentially small once W >> lambda. So a +k
 * test function barely sees a +k wave, and a basis built from +k alone makes
 * the WHOLE matrix negligible — the residual then comes out near zero for any
 * medium at all, and even the negative control passes. The test function that
 * actually pairs with a +k wave is the -k one (omega = 0). An earlier version
 * of this test had exactly that hole.
 * NELEM is in ELEMENTS, not wavelengths: sizing the domain in lambda made the
 * element count blow up at small W and silently truncated the neighbours of the
 * "interior" elements. */
enum { NELEM = 60 };

static double annihilation(double W, double complex kx, const hz_medseg *segs, int nseg) {
  static hz_carrier b[NB];
  double dom = (double)NELEM * W;
  int dim = 0;
  for (int n = -3; n <= NELEM + 3 && dim + 2 <= NB; n++) {
    b[dim++] = (hz_carrier){W, n, kx, -1};
    b[dim++] = (hz_carrier){W, n, -kx, -1};
  }
  /* the plane wave e^{i kx x}: constant coefficient on the +kx family only */
  double worst = 0.0;
  for (int i = 0; i < dim; i++) {
    double lo = W * ((double)b[i].n - 2.0), hi = W * ((double)b[i].n + 2.0);
    if (lo < 4.0 * W || hi > dom - 4.0 * W) continue;
    double complex acc = 0.0;
    double scale = 0.0;
    for (int j = 0; j < dim; j++) {
      double complex a = hz_carrier_entry(b[i], b[j], dom, 0.0, segs, nseg, NULL);
      /* the plane wave in a LOCALLY referenced basis needs coefficients
       * e^{i kx n W} / 4, not constants: each function already carries
       * e^{i kx (x - nW)}, so the node phase has to be put back */
      double complex cj = 0.0;
      if ((creal(b[j].kx) > 0.0) == (creal(kx) > 0.0))
        cj = cexp(CMPLX(0.0, 1.0) * b[j].kx * b[j].W * (double)b[j].n) / PU;
      acc += a * cj;
      scale += cabs(a) * cabs(cj);
    }
    if (scale <= 0.0) continue;
    double rel = cabs(acc) / scale;
    if (rel > worst) worst = rel;
  }
  return worst;
}

int main(void) {
  double k = 2.0 * M_PI / LAM;
  /* (1) homogeneous, lossless, carrier matched — must annihilate exactly */
  hz_medseg homo = {-1e9, 1e9, CMPLX(k * k, 0.0)};
  for (int e = 0; e < 4; e++) {
    double W = LAM * pow(10.0, (double)e - 1.0); /* 0.1, 1, 10, 100 lambda */
    /* tolerance grows with W: at 100 lambda the entries span a huge dynamic
     * range (the omega=2k parts are ~1e-15 of the omega=0 ones), so the
     * relative residual accumulates rounding. Measured: 1e-15 at 0.1-1 lambda,
     * 3e-10 at 10, 5e-7 at 100. */
    double tol = 5e-11 * (W / LAM) * (W / LAM);
    if (tol < 1e-12) tol = 1e-12;
    double r = annihilation(W, k, &homo, 1);
    check(r < tol, "plane wave annihilated (+k, homogeneous)", r);
    double rm = annihilation(W, -k, &homo, 1);
    check(rm < tol, "plane wave annihilated (-k, homogeneous)", rm);
    printf("  W = %6.1f lambda: residual %.2e (tol %.1e)\n", W / LAM, r, tol);
  }

  /* (2) with absorption: the carrier must be the COMPLEX local wavenumber for
   * the cancellation to survive. k2 = k^2 (1 + i a), so kx = k sqrt(1+ia). */
  double alpha = 0.05;
  double complex k2c = CMPLX(k * k, k * k * alpha);
  hz_medseg lossy = {-1e9, 1e9, k2c};
  double complex kx_full = csqrt(k2c);
  double r_lossy_real = annihilation(LAM * 10.0, creal(kx_full), &lossy, 1);
  /* THE POINT OF MAKING kx COMPLEX: with the full complex wavenumber the
   * cancellation is restored and the operator annihilates again. */
  double r_lossy_cplx = annihilation(LAM * 10.0, kx_full, &lossy, 1);
  /* a real carrier cannot match a complex medium exactly — this is expected to
   * be small but NOT machine zero, and it is recorded rather than asserted away */
  /* FINDING, not a failure: hz_carrier.kx is a double, so an ABSORBING medium
   * cannot be matched at all — its wavenumber is complex. With alpha = 0.05 the
   * envelope decays by ~2x across a 10-lambda element, which constant
   * coefficients cannot represent, and the residual is ~95%%. Two consequences:
   * (a) the API needs a COMPLEX kx before absorption can be handled;
   * (b) with the complex carrier there is NO representational size law from
   *     absorption at all: a decaying wave e^{i k~ x} is exactly in the span
   *     with constant coefficients, just like a lossless one. The only limit is
   *     FLOATING-POINT DYNAMIC RANGE across the support, exp(-Im(k)*4W): the
   *     sweep below shows the residual creeping up to 4e-5 at Im(k)*W = 10
   *     (span e^-40, the edge of double) and then becoming meaningless because
   *     the field itself has underflowed. An earlier comment here claimed
   *     "Im(k)*W <~ 1" as a representational law — that was wrong.
   * Still direct support for PLAN question 5: keeping the interior real avoids
   * the dynamic-range issue entirely and lets the interior operator be real
   * symmetric (MINRES). */
  printf("FINDING: real carrier vs lossy medium (alpha=%.2f): residual %.3e\n", alpha,
         r_lossy_real);
  printf("  => hz_carrier.kx must become complex before absorption is supported;\n");
  printf("     the limit is floating-point dynamic range exp(-Im(k)*4W), not representability.\n");
  check(r_lossy_real > 0.1, "lossy medium with a REAL carrier does NOT annihilate (expected)",
        r_lossy_real);
  printf("  with the COMPLEX carrier: residual %.3e\n", r_lossy_cplx);
  check(r_lossy_cplx < 5e-9, "lossy medium: COMPLEX carrier restores annihilation", r_lossy_cplx);
  /* and the size law Im(k)*W <~ 1: widen the element until the decay across it
   * stops being representable */
  printf("\ndynamic range vs Im(k)*W (complex carrier, alpha=%.2f):\n", alpha);
  for (int e = 0; e < 5; e++) {
    double Wl = LAM * pow(4.0, (double)e);
    double rr2 = annihilation(Wl, kx_full, &lossy, 1);
    printf("  W = %8.1f lambda   Im(k)*W = %7.3f   residual %.2e\n", Wl / LAM, cimag(kx_full) * Wl,
           rr2);
  }

  /* (3) NEGATIVE CONTROL: a medium that does not match the carrier must leave a
   * substantial residual, or the test above proves nothing */
  hz_medseg wrong = {-1e9, 1e9, CMPLX(2.25 * k * k, 0.0)};
  double rw = annihilation(LAM * 10.0, k, &wrong, 1);
  check(rw > 0.1, "negative control: mismatched medium leaves a residual", rw);

  /* (4) PARTITION OF UNITY UNDER A CUT — the reason a cut basis is not
   * automatically conforming. Truncated translates cannot sum to 4 near the
   * cut, so the exact solution is NOT reproduced there by constant
   * coefficients, however well the operator behaves. Quantified, not asserted:
   * this is the size of the problem PLAN question 10 has to solve. */
  double W = 10.0 * LAM;
  printf("\npartition of unity near a cut at x=0 (sum of truncated translates, exact = %.0f):\n",
         PU);
  for (int i = 0; i <= 6; i++) {
    double x = -1.5 * W + 0.5 * W * (double)i;
    double sum_left = 0.0, sum_all = 0.0;
    for (int n = -20; n <= 20; n++) {
      double t = x / W - (double)n;
      if (fabs(t) >= 2.0) continue;
      double v = hz_phi(t);
      sum_all += v;
      if ((double)n * W < 0.0) sum_left += v; /* only elements owning the left side */
    }
    printf("  x/W = %+5.2f   all %.4f   left-side only %.4f   deficit %.4f\n", x / W, sum_all,
           sum_left, sum_all - sum_left);
  }
  printf("READ: the deficit is what a cut removes from the left-hand span; it is\n");
  printf("zero far from the cut and large at it. Constant coefficients therefore\n");
  printf("stop reproducing a plane wave there — the non-conformity is structural,\n");
  printf("not a solver artefact.\n");

  printf("\ntest_carrier_op: %d/%d passed\n", g_total - g_fail, g_total);
  return g_fail > 0 ? 1 : 0;
}
