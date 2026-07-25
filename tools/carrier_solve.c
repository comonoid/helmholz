/* THE OPERATOR SOLVE — everything measured separately, assembled and run.
 *
 * Pieces, each validated on its own before being trusted here:
 *   carrier basis, complex kx, phase referenced to the element  (test_carrier_op)
 *   dyadic shells, redundancy tracks overlap not depth          (SHELLS_REPORT)
 *   cascade removes inter-level redundancy, full rank per level (CASCADE_REPORT)
 *   DtN termination, 0.7-5% with elements up to 200 lambda      (TERMINATION_REPORT)
 *   the assembly annihilates a plane wave to 1e-10              (test_carrier_op)
 *
 * WHAT THIS MEASURES AND WHY IT IS THE LAST ONE. Every result so far has been
 * about REPRESENTATION — best approximation, a lower bound on what any Galerkin
 * solve can reach. Representability does not imply the solve gets there:
 * Helmholtz is indefinite and quasi-optimality is not free. Falsifier 5 of M9c
 * is the only one still open, and it asks for the NUMBER OF SWEEPS to a given
 * accuracy OF THE FIELD (not of the residual — M2 measured that the raw
 * residual stagnates on a boundary mode that does not hurt the field).
 * Without that number "faster than path tracing" is unverifiable.
 *
 * Reference: the exact outgoing 1D Green's function, valid at any distance. */
#include "carrier.h"
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum { MAXB = 1000, NLEVMAX = 24, NSWEEP = 8 };
static const double LAM = 16.0;
static const double W0 = 2.0;
static const double EPS = 0.125;
static const double RCOND = 1e-12;

static double complex rhs_src(hz_carrier b, double xsrc) {
  hz_phi_factor fi = {b.W, (double)b.n, 0}, fs = {1.0, xsrc, 0};
  double lo = b.W * ((double)b.n - 2.0), hi = b.W * ((double)b.n + 2.0);
  if (lo < xsrc - 2.0) lo = xsrc - 2.0;
  if (hi > xsrc + 2.0) hi = xsrc + 2.0;
  if (lo >= hi) return 0.0;
  return cexp(CMPLX(0.0, -1.0) * b.kx * b.W * (double)b.n) *
         hz_phi_prod_integral_osc(lo, hi, fi, fs, b.kx);
}

/* value and derivative of a basis function at x (local phase reference) */
static void bval_der(hz_carrier b, double x, double complex *v, double complex *d) {
  double t = x / b.W - (double)b.n;
  *v = 0.0;
  *d = 0.0;
  if (fabs(t) >= 2.0) return;
  double complex e = cexp(CMPLX(0.0, 1.0) * b.kx * (x - b.W * (double)b.n));
  *v = hz_phi(t) * e;
  *d = (hz_phi_d1(t) / b.W + CMPLX(0.0, 1.0) * b.kx * hz_phi(t)) * e;
}

int main(void) {
  double k = 2.0 * M_PI / LAM;
  double complex ii = CMPLX(0.0, 1.0);
  /* Green amplitude for the source ACTUALLY used, f = phi(x - xs), not for a
   * delta: u = phihat(k) * e^{ik|x-xs|} / (2ik). Dropping phihat(k) makes the
   * reference wrong by a constant factor — the first run of this bench reported
   * exactly 2.8481 for every domain size, which is |phihat(k) - 1|, not an
   * error of the solve. */
  double complex phihat = 0.0;
  {
    static const double GX[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                                 -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                                 0.7966664774136267,  0.9602898564975363};
    static const double GW[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                                 0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                                 0.2223810344533745, 0.1012285362903763};
    for (int p = -2; p < 2; p++)
      for (int g = 0; g < 8; g++) {
        double t = (double)p + 0.5 + 0.5 * GX[g];
        phihat += 0.5 * GW[g] * hz_phi(t) * cexp(CMPLX(0.0, -1.0) * k * t);
      }
  }
  double complex amp = phihat / CMPLX(0.0, 2.0 * k);

  printf("Operator solve: carrier basis + dyadic shells + cascade + DtN\n");
  printf("Falsifier 5 of M9c: sweeps to a given accuracy OF THE FIELD.\n\n");
  printf("  %10s %6s %6s %10s %10s %10s %10s\n", "domain/lam", "lev", "dim", "direct", "sweep 1",
         "sweep 2", "sweep 4");

  static const double DOM[4] = {200.0, 2000.0, 20000.0, 200000.0};
  for (int id = 0; id < 4; id++) {
    double dom = DOM[id] * LAM;
    double xsrc = 0.5 * dom;

    static hz_carrier b[MAXB];
    static int lev0[NLEVMAX + 1];
    int dim = 0, nlev = 0;
    for (int j = 0; j < NLEVMAX; j++) {
      double W = W0 * pow(2.0, (double)j);
      if (W > 0.02 * dom) break;
      lev0[nlev] = dim;
      double alo = (j == 0) ? 0.0 : W / EPS, ahi = 2.0 * W / EPS;
      if (W0 * pow(2.0, (double)(j + 1)) > 0.02 * dom) ahi = 2.0 * dom;
      alo -= 2.0 * W;
      ahi += 2.0 * W;
      if (alo < 0.0) alo = 0.0;
      int n0 = (int)((xsrc - ahi) / W) - 2, n1 = (int)((xsrc + ahi) / W) + 2;
      int added = 0;
      for (int n = n0; n <= n1 && dim + 2 <= MAXB; n++) {
        double xc = (double)n * W, d = fabs(xc - xsrc);
        if (xc < -2.0 * W || xc > dom + 2.0 * W || d < alo || d >= ahi) continue;
        b[dim++] = (hz_carrier){W, n, CMPLX(k, 0.0), -1};
        b[dim++] = (hz_carrier){W, n, CMPLX(-k, 0.0), -1};
        added = 1;
      }
      if (added) nlev++;
    }
    lev0[nlev] = dim;

    hz_medseg vac = {0.0, dom, CMPLX(k * k, 0.0)};
    int nrow = dim + 2;
    double complex *A = calloc((size_t)nrow * (size_t)dim, sizeof(double complex));
    double complex *rr = calloc((size_t)nrow, sizeof(double complex));
    double *sv = calloc((size_t)dim, sizeof(double));
    static double complex coef[MAXB];
    if (!A || !rr || !sv) return 1;

    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++)
        A[(size_t)i * (size_t)dim + (size_t)j] =
            hz_carrier_entry(b[i], b[j], dom, 0.0, &vac, 1, NULL);
      rr[i] = rhs_src(b[i], xsrc);
    }
    /* DtN rows, scaled to the operator rows so they neither vanish nor dominate */
    double scl = 0.0;
    for (int j = 0; j < dim; j++)
      scl += cabs(A[(size_t)(dim / 2) * (size_t)dim + (size_t)j]);
    if (scl <= 0.0) scl = 1.0;
    for (int e = 0; e < 2; e++) {
      double xb = (e == 0) ? 0.0 : dom;
      double complex sgn = (e == 0) ? -ii * k : ii * k;
      for (int j = 0; j < dim; j++) {
        double complex v, d;
        bval_der(b[j], xb, &v, &d);
        A[(size_t)(dim + e) * (size_t)dim + (size_t)j] = scl * (d - sgn * v);
      }
      rr[dim + e] = 0.0;
    }

    /* field error against the analytic Green's function, away from the source
     * and away from the ends */
    double ma = xsrc + 0.1 * dom, mb = xsrc + 0.3 * dom;
    double err[5];
    for (int q = 0; q < 5; q++)
      err[q] = -1.0;

    /* --- direct solve, the accuracy yardstick --- */
    {
      double complex *Aw = calloc((size_t)nrow * (size_t)dim, sizeof(double complex));
      double complex *rw = calloc((size_t)nrow, sizeof(double complex));
      if (!Aw || !rw) return 1;
      for (size_t t = 0; t < (size_t)nrow * (size_t)dim; t++)
        Aw[t] = A[t];
      for (int t = 0; t < nrow; t++)
        rw[t] = rr[t];
      lapack_int rank = 0;
      if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, nrow, dim, 1, Aw, dim, rw, 1, sv, RCOND, &rank) == 0) {
        double num = 0.0, den = 0.0;
        for (int p = 0; p < 3000; p++) {
          double x = ma + (mb - ma) * ((double)p + 0.5) / 3000.0;
          double complex v = 0.0;
          for (int j = 0; j < dim; j++)
            v += rw[j] * hz_carrier_val(b[j], x, 0.0);
          double complex u = amp * cexp(ii * k * (x - xsrc));
          num += creal((v - u) * conj(v - u));
          den += creal(u * conj(u));
        }
        err[0] = sqrt(num / den);
      }
      free(Aw);
      free(rw);
    }

    /* --- cascade: block Gauss-Seidel over levels, coarse to fine --- */
    for (int i = 0; i < dim; i++)
      coef[i] = 0.0;
    for (int sw = 1; sw <= 4; sw++) {
      for (int L = nlev - 1; L >= 0; L--) { /* coarse first: level nlev-1 is widest */
        int a = lev0[L], nl = lev0[L + 1] - a;
        if (nl <= 0) continue;
        int nr2 = nl + 2;
        double complex *Al = calloc((size_t)nr2 * (size_t)nl, sizeof(double complex));
        double complex *rl = calloc((size_t)nr2, sizeof(double complex));
        double *svl = calloc((size_t)nl, sizeof(double));
        if (!Al || !rl || !svl) return 1;
        for (int i = 0; i < nl; i++) {
          for (int j = 0; j < nl; j++)
            Al[(size_t)i * (size_t)nl + (size_t)j] =
                A[(size_t)(a + i) * (size_t)dim + (size_t)(a + j)];
          double complex t = rr[a + i];
          for (int m = 0; m < dim; m++) {
            if (m >= a && m < a + nl) continue;
            if (!(cabs(coef[m]) > 0.0)) continue;
            t -= coef[m] * A[(size_t)(a + i) * (size_t)dim + (size_t)m];
          }
          rl[i] = t;
        }
        for (int e = 0; e < 2; e++) {
          double complex t = rr[dim + e];
          for (int j = 0; j < nl; j++)
            Al[(size_t)(nl + e) * (size_t)nl + (size_t)j] =
                A[(size_t)(dim + e) * (size_t)dim + (size_t)(a + j)];
          for (int m = 0; m < dim; m++) {
            if (m >= a && m < a + nl) continue;
            if (!(cabs(coef[m]) > 0.0)) continue;
            t -= coef[m] * A[(size_t)(dim + e) * (size_t)dim + (size_t)m];
          }
          rl[nl + e] = t;
        }
        lapack_int rk = 0;
        if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, nr2, nl, 1, Al, nl, rl, 1, svl, RCOND, &rk) == 0)
          for (int i = 0; i < nl; i++)
            coef[a + i] = rl[i];
        free(Al);
        free(rl);
        free(svl);
      }
      if (sw == 1 || sw == 2 || sw == 4) {
        double num = 0.0, den = 0.0;
        for (int p = 0; p < 3000; p++) {
          double x = ma + (mb - ma) * ((double)p + 0.5) / 3000.0;
          double complex v = 0.0;
          for (int j = 0; j < dim; j++)
            v += coef[j] * hz_carrier_val(b[j], x, 0.0);
          double complex u = amp * cexp(ii * k * (x - xsrc));
          num += creal((v - u) * conj(v - u));
          den += creal(u * conj(u));
        }
        err[sw == 1 ? 1 : (sw == 2 ? 2 : 3)] = sqrt(num / den);
      }
    }

    printf("  %10.0f %6d %6d %10.4f %10.4f %10.4f %10.4f\n", DOM[id], nlev, dim, err[0], err[1],
           err[2], err[3]);
    free(A);
    free(rr);
    free(sv);
  }
  printf("\nREAD: 'direct' is what the basis can do; the sweep columns are what the\n");
  printf("cascade reaches and how fast. Falsifier 5 asks whether a small, bounded\n");
  printf("number of sweeps gets there — and whether that number stays bounded as\n");
  printf("the domain grows by three orders of magnitude.\n");
  return 0;
}
