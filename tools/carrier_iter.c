/* FALSIFIER 5 of M9c, measured directly: ITERATIONS TO A GIVEN FIELD ACCURACY,
 * and whether that count stays bounded as the domain grows.
 *
 * Order matters here and the audit of this bench reversed it. The obvious move
 * was "measure the spectrum, then pick a method", but for this operator that is
 * naive: the matrix is complex-symmetric, hence NON-NORMAL, and for non-normal
 * operators Krylov convergence is not governed by eigenvalues at all. Worse,
 * the basis is degenerate BY CONSTRUCTION (Nyquist null mode plus inter-level
 * redundancy), so its exact zeros have nothing to do with the physics and a
 * naive "deflate the small modes" would deflate the basis, not the operator.
 * So: measure the thing the falsifier actually asks for — the iteration count —
 * and reach for spectra only to explain a failure.
 *
 * METHOD: LSQR (Paige-Saunders) on the rectangular system, which is what the
 * plan has specified from the start (min-norm, never form A^H A). It needs only
 * A*v and A^H*u, so it says nothing about the assembly and everything about the
 * operator.
 * ACCURACY IS OF THE FIELD, not of the residual: M2 measured that the raw
 * residual stagnates on a boundary mode that does not hurt the field, so a
 * residual-based count would be meaningless here. */
#include "carrier.h"
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum { MAXB = 1000, NLEVMAX = 24, MAXIT = 400, NPROBE = 1500 };
static const double LAM = 16.0;
static const double W0 = 2.0;
static const double EPS = 0.125;
static const double RCOND = 1e-12;
static const double TARGET = 0.01;

static double complex rhs_src(hz_carrier b, double xsrc) {
  hz_phi_factor fi = {b.W, (double)b.n, 0}, fs = {1.0, xsrc, 0};
  double lo = b.W * ((double)b.n - 2.0), hi = b.W * ((double)b.n + 2.0);
  if (lo < xsrc - 2.0) lo = xsrc - 2.0;
  if (hi > xsrc + 2.0) hi = xsrc + 2.0;
  if (lo >= hi) return 0.0;
  return cexp(CMPLX(0.0, -1.0) * b.kx * b.W * (double)b.n) *
         hz_phi_prod_integral_osc(lo, hi, fi, fs, b.kx);
}

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

  printf("Falsifier 5: LSQR iterations to %.0f%% FIELD accuracy vs domain size.\n", TARGET * 100.0);
  printf("The question is not whether it converges but whether the count stays\n");
  printf("bounded while the domain grows by three orders of magnitude.\n\n");
  printf("  %10s %5s %6s %10s %8s %10s\n", "domain/lam", "lev", "dim", "direct err", "iters",
         "err@iters");

  static const double DOM[4] = {200.0, 2000.0, 20000.0, 200000.0};
  for (int id = 0; id < 4; id++) {
    double dom = DOM[id] * LAM, xsrc = 0.5 * dom;
    static hz_carrier b[MAXB];
    int dim = 0, nlev = 0;
    for (int j = 0; j < NLEVMAX; j++) {
      double W = W0 * pow(2.0, (double)j);
      if (W > 0.02 * dom) break;
      double alo = (j == 0) ? 0.0 : W / EPS, ahi = 2.0 * W / EPS;
      if (W0 * pow(2.0, (double)(j + 1)) > 0.02 * dom) ahi = 2.0 * dom;
      alo -= 2.0 * W;
      ahi += 2.0 * W;
      if (alo < 0.0) alo = 0.0;
      int n0 = (int)((xsrc - ahi) / W) - 2, n1 = (int)((xsrc + ahi) / W) + 2, added = 0;
      for (int n = n0; n <= n1 && dim + 2 <= MAXB; n++) {
        double xc = (double)n * W, d = fabs(xc - xsrc);
        if (xc < -2.0 * W || xc > dom + 2.0 * W || d < alo || d >= ahi) continue;
        b[dim++] = (hz_carrier){W, n, CMPLX(k, 0.0), -1};
        b[dim++] = (hz_carrier){W, n, CMPLX(-k, 0.0), -1};
        added = 1;
      }
      if (added) nlev++;
    }

    hz_medseg vac = {0.0, dom, CMPLX(k * k, 0.0)};
    int nrow = dim + 2;
    double complex *A = calloc((size_t)nrow * (size_t)dim, sizeof(double complex));
    double complex *rhs = calloc((size_t)nrow, sizeof(double complex));
    if (!A || !rhs) return 1;
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++)
        A[(size_t)i * (size_t)dim + (size_t)j] =
            hz_carrier_entry(b[i], b[j], dom, 0.0, &vac, 1, NULL);
      rhs[i] = rhs_src(b[i], xsrc);
    }
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
      rhs[dim + e] = 0.0;
    }

    /* probe points and the exact field there */
    static double px[NPROBE];
    static double complex pu[NPROBE];
    double ma = xsrc + 0.1 * dom, mb = xsrc + 0.3 * dom, den = 0.0;
    for (int p = 0; p < NPROBE; p++) {
      px[p] = ma + (mb - ma) * ((double)p + 0.5) / (double)NPROBE;
      pu[p] = amp * cexp(ii * k * (px[p] - xsrc));
      den += creal(pu[p] * conj(pu[p]));
    }

    /* direct solve: the accuracy floor LSQR is aiming at */
    double direct = -1.0;
    {
      double complex *Aw = calloc((size_t)nrow * (size_t)dim, sizeof(double complex));
      double complex *rw = calloc((size_t)nrow, sizeof(double complex));
      double *sv = calloc((size_t)dim, sizeof(double));
      if (!Aw || !rw || !sv) return 1;
      for (size_t t = 0; t < (size_t)nrow * (size_t)dim; t++)
        Aw[t] = A[t];
      for (int t = 0; t < nrow; t++)
        rw[t] = rhs[t];
      lapack_int rank = 0;
      if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, nrow, dim, 1, Aw, dim, rw, 1, sv, RCOND, &rank) == 0) {
        double num = 0.0;
        for (int p = 0; p < NPROBE; p++) {
          double complex v = 0.0;
          for (int j = 0; j < dim; j++)
            v += rw[j] * hz_carrier_val(b[j], px[p], 0.0);
          num += creal((v - pu[p]) * conj(v - pu[p]));
        }
        direct = sqrt(num / den);
      }
      free(Aw);
      free(rw);
      free(sv);
    }

    /* --- LSQR (Paige-Saunders), complex, no normal equations formed --- */
    double complex *u = calloc((size_t)nrow, sizeof(double complex));
    double complex *v = calloc((size_t)dim, sizeof(double complex));
    double complex *w = calloc((size_t)dim, sizeof(double complex));
    double complex *x = calloc((size_t)dim, sizeof(double complex));
    double complex *tmp = calloc((size_t)(nrow > dim ? nrow : dim), sizeof(double complex));
    if (!u || !v || !w || !x || !tmp) return 1;

    double beta = 0.0;
    for (int i = 0; i < nrow; i++) {
      u[i] = rhs[i];
      beta += creal(u[i] * conj(u[i]));
    }
    beta = sqrt(beta);
    for (int i = 0; i < nrow; i++)
      u[i] /= beta;
    double alpha = 0.0;
    for (int j = 0; j < dim; j++) {
      double complex s = 0.0;
      for (int i = 0; i < nrow; i++)
        s += conj(A[(size_t)i * (size_t)dim + (size_t)j]) * u[i];
      v[j] = s;
      alpha += creal(s * conj(s));
    }
    alpha = sqrt(alpha);
    for (int j = 0; j < dim; j++) {
      v[j] /= alpha;
      w[j] = v[j];
      x[j] = 0.0;
    }
    double phibar = beta, rhobar = alpha;

    int iters = -1;
    double errit = -1.0;
    for (int it = 1; it <= MAXIT; it++) {
      for (int i = 0; i < nrow; i++) {
        double complex s = 0.0;
        for (int j = 0; j < dim; j++)
          s += A[(size_t)i * (size_t)dim + (size_t)j] * v[j];
        tmp[i] = s - alpha * u[i];
      }
      double bnew = 0.0;
      for (int i = 0; i < nrow; i++)
        bnew += creal(tmp[i] * conj(tmp[i]));
      bnew = sqrt(bnew);
      if (bnew > 0.0)
        for (int i = 0; i < nrow; i++)
          u[i] = tmp[i] / bnew;
      for (int j = 0; j < dim; j++) {
        double complex s = 0.0;
        for (int i = 0; i < nrow; i++)
          s += conj(A[(size_t)i * (size_t)dim + (size_t)j]) * u[i];
        tmp[j] = s - bnew * v[j];
      }
      double anew = 0.0;
      for (int j = 0; j < dim; j++)
        anew += creal(tmp[j] * conj(tmp[j]));
      anew = sqrt(anew);
      if (anew > 0.0)
        for (int j = 0; j < dim; j++)
          v[j] = tmp[j] / anew;

      double rho = sqrt(rhobar * rhobar + bnew * bnew);
      double c = rhobar / rho, s2 = bnew / rho;
      double theta = s2 * anew;
      rhobar = -c * anew;
      double phi = c * phibar;
      phibar = s2 * phibar;
      for (int j = 0; j < dim; j++) {
        x[j] += (phi / rho) * w[j];
        w[j] = v[j] - (theta / rho) * w[j];
      }
      alpha = anew;

      if (it % 5 == 0 || it == 1) {
        double num = 0.0;
        for (int p = 0; p < NPROBE; p++) {
          double complex vv = 0.0;
          for (int j = 0; j < dim; j++)
            vv += x[j] * hz_carrier_val(b[j], px[p], 0.0);
          num += creal((vv - pu[p]) * conj(vv - pu[p]));
        }
        errit = sqrt(num / den);
        if (errit <= TARGET) {
          iters = it;
          break;
        }
      }
    }
    printf("  %10.0f %5d %6d %10.2e %8d %10.2e\n", DOM[id], nlev, dim, direct, iters, errit);
    free(A);
    free(rhs);
    free(u);
    free(v);
    free(w);
    free(x);
    free(tmp);
  }
  printf("\nREAD: iters = -1 means %d LSQR iterations did not reach the target.\n", MAXIT);
  printf("A count that grows with the domain kills the cold-solve real-time claim\n");
  printf("and routes the argument to the incremental path (PLAN M11) instead.\n");
  return 0;
}
