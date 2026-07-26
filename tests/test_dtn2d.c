/* GATE D2-1 of the vertical slice: the DtN rows, checked WITHOUT a solver.
 *
 * The audited plan asked for "project the exact scattered field onto the basis
 * and check the rows are satisfied to within the projection error", which needs
 * a Gram matrix and a least-squares solve at dim ~2000. There is a cheaper and
 * STRICTER check, and it is the D1 trick again: a PLANE WAVE lies in the span
 * exactly, with coefficients known in closed form (partition of unity, /16), and
 * its angular harmonics are Bessel functions. So the assembled rows applied to
 * that known vector must equal a number one can write down.
 *
 * AND THE NUMBER IS NOT ZERO, WHICH IS THE POINT. A plane wave is not outgoing:
 * with 2 J_m = H1_m + H2_m and Z_m = H1'_m/H1_m,
 *     J'_m - Z_m J_m = (H2'_m - Z_m H2_m) / 2,
 * so the row measures exactly the INCOMING half of the plane wave. A test whose
 * expected value is zero can pass by everything being zero — this one cannot.
 * It checks the boundary sampling, the DFT normalisation, the radial derivative,
 * the symbol and the local phase convention at once, against a closed form.
 *
 * NC-B (the slice's own negative control): replace the outgoing symbol by the
 * incoming one. The rows must then annihilate the INCOMING half and return the
 * outgoing one instead — a different, also computable number. Getting the sign
 * of the radiation condition backwards is the single most likely error in a
 * termination, and this pins it. */
#include "bessel.h"
#include "carrier2d.h"
#include "dtn2d.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>

enum { NN = 9, MAXB = (2 * NN + 1) * (2 * NN + 1) };
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

/* H^(2)_m and its derivative, from the same J, Y as everything else */
static void h2(int m, double x, double complex *h, double complex *dh) {
  double j, y, dj, dy;
  hz_bessel_jy(m, x, &j, &y, &dj, &dy);
  *h = CMPLX(j, -y);
  *dh = CMPLX(dj, -dy);
}

int main(void) {
  printf("test_dtn2d: DtN rows on a circle, against a closed form\n");
  double lam = 1.0, k0 = 2.0 * M_PI / lam, W = 1.0 * lam, R = 5.0 * lam;
  hz_dtn d = {R, k0, 40, 512};

  /* the plane wave travelling along +x, and its exact expansion coefficients */
  double kx = k0, ky = 0.0;
  static hz_carrier2d b[MAXB];
  static double complex c[MAXB];
  int dim = 0;
  for (int nx = -NN; nx <= NN; nx++)
    for (int ny = -NN; ny <= NN; ny++) {
      b[dim] = (hz_carrier2d){W, nx, ny, CMPLX(kx, 0.0), CMPLX(ky, 0.0)};
      c[dim] = cexp(CMPLX(0.0, 1.0) * (kx * W * (double)nx + ky * W * (double)ny)) / 16.0;
      dim++;
    }

  /* --- [1] the rows reproduce the incoming half of the plane wave --------- */
  {
    double worst = 0.0;
    int mworst = 0;
    for (int m = 0; m <= d.mmax; m++) {
      double complex got = 0.0;
      for (int j = 0; j < dim; j++)
        got += c[j] * hz_dtn_row(&d, m, b[j]);
      double complex H2, dH2;
      h2(m, k0 * R, &H2, &dH2);
      double complex ipow = cexp(CMPLX(0.0, 1.0) * M_PI * 0.5 * (double)m);
      double complex want = ipow * k0 * (dH2 - hz_dtn_symbol(&d, m) / k0 * H2) * 0.5;
      double sc = k0 * (cabs(dH2) + cabs(hz_dtn_symbol(&d, m) / k0 * H2)) * 0.5;
      double e = cabs(got - want) / sc;
      if (e > worst) {
        worst = e;
        mworst = m;
      }
    }
    char nm[64];
    snprintf(nm, sizeof nm, "[1] rows == incoming half, worst m=%d", mworst);
    check(nm, worst < 1e-9, worst);
  }

  /* --- [1a,1b] the two ASSEMBLED harmonics, each against its closed form ---
   * Without these, test [1] alone cannot tell a wrong derivative from a wrong
   * symbol: the row is one number and the two errors could cancel in it. */
  {
    double wu = 0.0, wd = 0.0;
    for (int m = 0; m <= d.mmax; m++) {
      double complex uh = 0.0, duh = 0.0;
      for (int j = 0; j < dim; j++) {
        double complex a, b2;
        hz_dtn_harm(&d, m, b[j], &a, &b2);
        uh += c[j] * a;
        duh += c[j] * b2;
      }
      double jj, y, dj, dy;
      hz_bessel_jy(m, k0 * R, &jj, &y, &dj, &dy);
      double complex ipow = cexp(CMPLX(0.0, 1.0) * M_PI * 0.5 * (double)m);
      double eu = cabs(uh - ipow * jj), ed = cabs(duh - ipow * k0 * dj) / k0;
      if (eu > wu) wu = eu;
      if (ed > wd) wd = ed;
    }
    check("[1a] assembled harmonic == i^m J_m(k0 R)", wu < 1e-12, wu);
    check("[1b] assembled radial derivative == i^m k0 J'_m", wd < 1e-12, wd);
  }

  /* --- [2] NC-B: the incoming symbol must keep the OTHER half ------------- */
  {
    double worst = 0.0;
    for (int m = 0; m <= 12; m++) {
      double complex H2, dH2;
      double j, y, dj, dy;
      h2(m, k0 * R, &H2, &dH2);
      hz_bessel_jy(m, k0 * R, &j, &y, &dj, &dy);
      double complex H1 = CMPLX(j, y), dH1 = CMPLX(dj, dy);
      double complex Z2 = k0 * dH2 / H2;
      /* the row with the WRONG (incoming) symbol, built from the SAME assembled
       * harmonics — not from the analytic ones, or this would be an algebraic
       * consequence of [1] rather than a control */
      double complex uh = 0.0, duh = 0.0;
      for (int jj = 0; jj < dim; jj++) {
        double complex a, b2;
        hz_dtn_harm(&d, m, b[jj], &a, &b2);
        uh += c[jj] * a;
        duh += c[jj] * b2;
      }
      double complex row2 = duh - Z2 * uh;
      double complex ipow = cexp(CMPLX(0.0, 1.0) * M_PI * 0.5 * (double)m);
      double complex want2 = ipow * (k0 * dH1 - Z2 * H1) * 0.5;
      double sc = cabs(k0 * dH1) + cabs(Z2 * H1);
      double e = cabs(row2 - want2) / sc;
      if (e > worst) worst = e;
    }
    check("[2] NC-B: incoming symbol keeps the outgoing half", worst < 1e-9, worst);
  }

  /* --- [3] a purely OUTGOING field is annihilated ------------------------- */
  /* Not representable in the basis, so it is fed directly: the row functional
   * applied to u = H1_m e^{i m theta} must vanish. This is the statement the
   * whole termination rests on, and it is independent of the basis. */
  {
    double worst = 0.0;
    for (int m = 0; m <= 20; m++) {
      double j, y, dj, dy;
      hz_bessel_jy(m, k0 * R, &j, &y, &dj, &dy);
      double complex H1 = CMPLX(j, y), dH1 = CMPLX(dj, dy);
      double complex res = k0 * dH1 - hz_dtn_symbol(&d, m) * H1;
      double e = cabs(res) / (cabs(k0 * dH1) + cabs(hz_dtn_symbol(&d, m) * H1));
      if (e > worst) worst = e;
    }
    check("[3] outgoing harmonic annihilated by its own symbol", worst < 1e-14, worst);
  }

  /* --- [4] sampling: doubling ntheta must not move the answer ------------- */
  {
    hz_dtn d2 = {R, k0, 40, 1024};
    double worst = 0.0;
    for (int m = 0; m <= d.mmax; m += 8) {
      double complex g1 = 0.0, g2 = 0.0;
      for (int j = 0; j < dim; j++) {
        g1 += c[j] * hz_dtn_row(&d, m, b[j]);
        g2 += c[j] * hz_dtn_row(&d2, m, b[j]);
      }
      double e = cabs(g1 - g2) / (cabs(g1) + cabs(g2));
      if (e > worst) worst = e;
    }
    check("[4] ntheta 512 vs 1024 agree", worst < 1e-9, worst);
  }

  /* --- [5] a basis function far from the circle contributes strictly 0 ---- */
  {
    hz_carrier2d far = {W, 0, 0, CMPLX(k0, 0.0), CMPLX(0.0, 0.0)};
    double complex r = hz_dtn_row(&d, 3, far);
    check("[5] element at the centre never reaches the circle", !(cabs(r) > 0.0), cabs(r));
  }

  printf("test_dtn2d: %d/%d passed\n", pass_count, pass_count + fail_count);
  return fail_count == 0 ? 0 : 1;
}
