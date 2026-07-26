/* GATE D1 of the vertical slice (SLICE_PLAN.md step 1): the 2D carrier assembly.
 *
 * THE ACCEPTANCE TRICK, LIFTED FROM 1D (tests/test_carrier_op.c, 11/11). In a
 * homogeneous medium a plane wave is an EXACT solution, and by the partition of
 * unity it lies in the span of the carrier basis with KNOWN coefficients.
 * Therefore the discrete operator must ANNIHILATE it. That checks the assembly,
 * the k^2 cancellation, the partition of unity and the medium term at once -
 * exactly, with no boundary and no solve. It is the cheapest decisive test this
 * project has, and it is the reason step 1 has a hard gate.
 *
 * THE COEFFICIENTS. Sum_n phi(t - n) = 4, so on a tensor grid Sum_{nx,ny} = 16.
 * With the phase referenced to the ELEMENT (mandatory once k is complex, see
 * carrier.h), reproducing e^{i(kx x + ky y)} needs
 *      c_{nx,ny} = e^{i(kx nx W + ky ny W)} / 16
 * for the basis functions whose carrier equals (kx,ky), and 0 for every other
 * direction.
 *
 * THE HOLE THE 1D VERSION HAD, AND WHY IT MATTERS HERE TOO: the form is
 * BILINEAR, so <B_i, L B_j> with both carriers at +k carries e^{2ikx} and is
 * exponentially small when W >> lambda. The test function that actually SEES a
 * +k wave is the one at -k. A basis holding one direction only makes the whole
 * matrix negligible and the negative control passes falsely. So the direction
 * fan here always contains the opposite direction as well. */
#include "carrier2d.h"
#include "phi.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum {
  ND = 8,    /* directions per element; EVEN, so the fan is closed under negation */
  NN = 11,   /* nodes per axis */
  PMARG = 4, /* interior margin: see below */
  MAXB = ND * NN * NN
};

/* WHY THE MARGIN IS 4 AND NOT 2 - the first run of this gate got it wrong and
 * failed with a residual of 1.1e-2 that was almost INDEPENDENT of W, which in
 * this project has meant "bench artefact" every single time (artefact 12 was a
 * constant 2.8481 at every domain size). supp(phi) is 4W wide, so reproducing
 * the plane wave over the support of row i needs every translate from nx-4 to
 * nx+4 to be present. With margin 2 the partition of unity is truncated inside
 * the very rows being tested, and the rim of the patch gets measured instead of
 * the assembly. With NN = 7 and margin 4 there is no interior row at all -
 * hence NN = 11. */

static int pass_count = 0, fail_count = 0;

static void check(const char *name, int ok, double v) {
  if (ok) {
    pass_count++;
    printf("  ok   %-44s %.3e\n", name, v);
    return;
  }
  fail_count++;
  printf("  FAIL %-44s %.3e\n", name, v);
}

/* Build a fan of ND directions on the circle plus the full node grid. The fan
 * is closed under negation (ND even), which is what keeps the bilinear form
 * from being blind - see the header note. */
static int build(hz_carrier2d *b, double W, double k, int n0, int n1, double *thetas) {
  int d = 0;
  for (int nx = n0; nx <= n1; nx++)
    for (int ny = n0; ny <= n1; ny++)
      for (int t = 0; t < ND; t++) {
        b[d].W = W;
        b[d].nx = nx;
        b[d].ny = ny;
        b[d].kx = CMPLX(k * cos(thetas[t]), 0.0);
        b[d].ky = CMPLX(k * sin(thetas[t]), 0.0);
        d++;
      }
  return d;
}

/* Residual ||A c|| / (||A|| ||c||) for the plane-wave coefficient vector, over
 * the rows of elements whose support lies strictly INSIDE the patch (only there
 * is the partition of unity complete - at the rim it is truncated, and a
 * nonzero residual there is the boundary of the patch, not an assembly error). */
static double annihilation(const hz_carrier2d *b, int dim, double W, int n0, int n1,
                           const double complex *c, const hz_medrect2d *rects, int nrect,
                           double domx, double domy) {
  double num = 0.0, den = 0.0;
  for (int i = 0; i < dim; i++) {
    if (b[i].nx < n0 + PMARG || b[i].nx > n1 - PMARG) continue;
    if (b[i].ny < n0 + PMARG || b[i].ny > n1 - PMARG) continue;
    double complex s = 0.0, scale = 0.0;
    for (int j = 0; j < dim; j++) {
      double complex a = hz_carrier2d_entry(b[i], b[j], domx, domy, rects, nrect, NULL);
      s += a * c[j];
      scale += cabs(a) * cabs(c[j]);
    }
    num += creal(s * conj(s));
    den += creal(scale * conj(scale));
  }
  (void)W;
  return (den > 0.0) ? sqrt(num / den) : -1.0;
}

int main(void) {
  printf("test_carrier2d: gate D1 of the vertical slice\n");

  static double thetas[ND];
  for (int t = 0; t < ND; t++)
    thetas[t] = 2.0 * M_PI * (double)t / (double)ND;

  static hz_carrier2d b[MAXB];
  static double complex c[MAXB];

  /* --- [1] separability sanity: the value of a basis function factorises --- */
  {
    hz_carrier2d e = {3.0, 2, 1, CMPLX(0.3, 0.0), CMPLX(-0.2, 0.0)};
    double x = 4.1, y = 2.7;
    double complex v = hz_carrier2d_val(e, x, y);
    double complex want = hz_phi(x / e.W - 2.0) * hz_phi(y / e.W - 1.0) *
                          cexp(CMPLX(0.0, 1.0) * (e.kx * (x - 3.0 * 2.0) + e.ky * (y - 3.0 * 1.0)));
    check("value factorises with local phase", cabs(v - want) < 1e-15, cabs(v - want));
    check("strictly zero outside the support",
          !(cabs(hz_carrier2d_val(e, 100.0, 2.7)) > 0.0) &&
              !(cabs(hz_carrier2d_val(e, 4.1, -100.0)) > 0.0),
          0.0);
  }

  /* --- [2] THE GATE: annihilation of a plane wave, several element widths ---
   * W is given in wavelengths, so this also answers "how wide may an element
   * be" - the 1D answer was "arbitrarily, the wave is in the span exactly". */
  static const double LAMS[4] = {0.5, 2.0, 8.0, 32.0}; /* W / lambda */
  double lam = 16.0, k = 2.0 * M_PI / lam;
  for (int w = 0; w < 4; w++) {
    double W = LAMS[w] * lam;
    int n0 = 0, n1 = NN - 1;
    double domx = W * ((double)NN + 4.0), domy = domx;
    int dim = build(b, W, k, n0, n1, thetas);
    hz_medrect2d vac = {0.0, domx, 0.0, domy, CMPLX(k * k, 0.0)};

    int dir = 1; /* not axis-aligned, so both axes are exercised */
    for (int j = 0; j < dim; j++) {
      int t = j % ND;
      if (t != dir) {
        c[j] = 0.0;
        continue;
      }
      double complex ph =
          CMPLX(0.0, 1.0) * (b[j].kx * W * (double)b[j].nx + b[j].ky * W * (double)b[j].ny);
      c[j] = cexp(ph) / 16.0;
    }
    double r = annihilation(b, dim, W, n0, n1, c, &vac, 1, domx, domy);
    char nm[64];
    snprintf(nm, sizeof nm, "annihilation, W = %.1f lambda", LAMS[w]);
    check(nm, r >= 0.0 && r < 1e-9, r);
  }

  /* --- [3] NEGATIVE CONTROL: a medium that does not match the carrier must
   * leave a visible residual. Without this, an assembly that returns garbage of
   * small magnitude would pass [2] trivially. --- */
  {
    double W = 2.0 * lam;
    int n0 = 0, n1 = NN - 1;
    double domx = W * ((double)NN + 4.0), domy = domx;
    int dim = build(b, W, k, n0, n1, thetas);
    double kk = 1.5 * k;
    hz_medrect2d bad = {0.0, domx, 0.0, domy, CMPLX(kk * kk, 0.0)};
    int dir = 1;
    for (int j = 0; j < dim; j++) {
      int t = j % ND;
      if (t != dir) {
        c[j] = 0.0;
        continue;
      }
      double complex ph =
          CMPLX(0.0, 1.0) * (b[j].kx * W * (double)b[j].nx + b[j].ky * W * (double)b[j].ny);
      c[j] = cexp(ph) / 16.0;
    }
    double r = annihilation(b, dim, W, n0, n1, c, &bad, 1, domx, domy);
    check("negative control: mismatched medium", r > 1e-3, r);
  }

  /* --- [4] the source integral factorises and is strictly zero when the
   * supports miss each other --- */
  {
    hz_carrier2d e = {4.0, 3, 3, CMPLX(0.1, 0.0), CMPLX(0.1, 0.0)};
    double complex q = hz_carrier2d_rhs(e, 100.0, 100.0, 12.0, 12.0, 1.0, NULL);
    check("source overlap is nonzero when it should be", cabs(q) > 1e-12, cabs(q));
    double complex z = hz_carrier2d_rhs(e, 100.0, 100.0, 90.0, 12.0, 1.0, NULL);
    check("source overlap strictly 0 when disjoint", !(cabs(z) > 0.0), cabs(z));
  }

  printf("test_carrier2d: %d/%d passed\n", pass_count, pass_count + fail_count);
  return fail_count == 0 ? 0 : 1;
}
