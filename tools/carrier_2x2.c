/* M9c, the 2x2 matrix with its negative control — the last structural
 * falsifier of the milestone.
 *
 *              uniform fine (lambda/8)      LOD (coarse with distance)
 *   plain             A                              B   <-- MUST FAIL
 *   carrier           C                              D
 *
 * Cell B is the point. The architecture's whole claim is that LOD is only
 * affordable BECAUSE the basis carries the oscillation; if plain potentials
 * survived LOD too, the carrier would be an optimisation rather than a
 * condition, and the plan's central assertion would be wrong. A negative
 * control that is PREDICTED to fail is what makes the other three cells mean
 * anything — without it a bench that passes everything proves nothing.
 *
 * Everything is closed-form: a plain potential is just a carrier with kx = 0,
 * so all four cells go through the same code path and differ only in the
 * wavenumber attached to each function. The field is the exact outgoing 1D
 * Green's function of a point source — constant amplitude, valid at any
 * distance, no grid (the reference rule this project has followed throughout).
 */
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum { MAXB = 2400 };
static const double LAM = 16.0;
static const double W0 = 2.0;    /* lambda/8 */
static const double EPS = 0.125; /* LOD rule, the value SHELLS_REPORT settled on */
static const double RCOND = 1e-12;

typedef struct {
  double W;
  int n;
  double kx; /* 0 = plain potential, +-k = carrier */
} bfn;

static double complex gram(bfn a, bfn b) {
  hz_phi_factor fa = {a.W, (double)a.n, 0}, fb = {b.W, (double)b.n, 0};
  double lo = a.W * ((double)a.n - 2.0), hi = a.W * ((double)a.n + 2.0);
  double lo2 = b.W * ((double)b.n - 2.0), hi2 = b.W * ((double)b.n + 2.0);
  if (lo2 > lo) lo = lo2;
  if (hi2 < hi) hi = hi2;
  if (lo >= hi) return 0.0;
  return hz_phi_prod_integral_osc(lo, hi, fa, fb, b.kx - a.kx);
}

/* <u, B_i> for u = A e^{i k (x - xs)} on x > xs and A e^{-i k (x - xs)} below */
static double complex rhs_of(bfn a, double k, double xs, double complex amp) {
  hz_phi_factor fa = {a.W, (double)a.n, 0};
  double lo = a.W * ((double)a.n - 2.0), hi = a.W * ((double)a.n + 2.0);
  double complex t = 0.0;
  if (hi > xs) {
    double l = lo > xs ? lo : xs;
    t += amp * cexp(CMPLX(0.0, -1.0) * k * xs) * hz_phi_integral_osc(l, hi, fa, k - a.kx);
  }
  if (lo < xs) {
    double h = hi < xs ? hi : xs;
    t += amp * cexp(CMPLX(0.0, 1.0) * k * xs) * hz_phi_integral_osc(lo, h, fa, -k - a.kx);
  }
  return t;
}

int main(void) {
  double k = 2.0 * M_PI / LAM;
  double complex amp = 1.0 / CMPLX(0.0, 2.0 * k); /* 1D Green's amplitude */
  double dom = 2000.0 * LAM;
  double xs = 0.5 * dom;

  printf("M9c 2x2 with the negative control: does LOD survive WITHOUT a carrier?\n");
  printf("domain %.0f lambda, source in the middle, exact 1D Green reference\n\n", dom / LAM);
  printf("  %-26s %8s %10s\n", "basis", "dim", "rel err");

  static bfn b[MAXB];
  for (int cell = 0; cell < 4; cell++) {
    int carrier = cell % 2; /* 0 = plain, 1 = carrier */
    int lod = cell / 2;     /* 0 = uniform fine, 1 = LOD */
    int dim = 0;

    if (!lod) {
      /* uniform lambda/8 over a window around the source — the whole domain at
       * this spacing would be 32000 elements, far past any sane cap, which is
       * itself the v1 cost statement */
      double half = 60.0 * LAM;
      int n0 = (int)((xs - half) / W0), n1 = (int)((xs + half) / W0);
      for (int n = n0; n <= n1 && dim + 2 <= MAXB; n++) {
        if (carrier) {
          b[dim++] = (bfn){W0, n, k};
          b[dim++] = (bfn){W0, n, -k};
        } else {
          b[dim++] = (bfn){W0, n, 0.0};
        }
      }
    } else {
      /* dyadic shells around the source, W_j = W0 * 2^j on [W_j/eps, 2W_j/eps) */
      for (int j = 0; j < 20; j++) {
        double W = W0 * pow(2.0, (double)j);
        if (W > 0.05 * dom) break;
        double alo = (j == 0) ? 0.0 : W / EPS;
        double ahi = 2.0 * W / EPS;
        alo -= 2.0 * W;
        ahi += 2.0 * W;
        if (alo < 0.0) alo = 0.0;
        int n0 = (int)((xs - ahi) / W) - 2, n1 = (int)((xs + ahi) / W) + 2;
        for (int n = n0; n <= n1 && dim + 2 <= MAXB; n++) {
          double xc = (double)n * W, d = fabs(xc - xs);
          if (xc < 0.0 || xc > dom || d < alo || d >= ahi) continue;
          if (carrier) {
            b[dim++] = (bfn){W, n, k};
            b[dim++] = (bfn){W, n, -k};
          } else {
            b[dim++] = (bfn){W, n, 0.0};
          }
        }
      }
    }

    /* measure over a window well away from the source, where the field is a
     * clean outgoing wave and the only question is whether the basis holds it */
    double ma = xs + 20.0 * LAM, mb = xs + 40.0 * LAM;
    double complex *G = calloc((size_t)dim * (size_t)dim, sizeof(double complex));
    double complex *r = calloc((size_t)dim, sizeof(double complex));
    double *sv = calloc((size_t)dim, sizeof(double));
    if (!G || !r || !sv) return 1;
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++)
        G[(size_t)i * (size_t)dim + (size_t)j] = gram(b[i], b[j]);
      r[i] = rhs_of(b[i], k, xs, amp);
    }
    static double complex rr[MAXB];
    for (int i = 0; i < dim; i++)
      rr[i] = r[i];
    lapack_int rank = 0;
    lapack_int info =
        LAPACKE_zgelsd(LAPACK_ROW_MAJOR, dim, dim, 1, G, dim, rr, 1, sv, RCOND, &rank);

    double err = -1.0;
    if (info == 0) {
      /* error over the measurement window, sampled: the basis functions there
       * are a small subset, so this is cheap and independent of the closed
       * forms used to build the system */
      double num = 0.0, den = 0.0;
      int np = 4000;
      for (int p = 0; p < np; p++) {
        double x = ma + (mb - ma) * ((double)p + 0.5) / (double)np;
        double complex v = 0.0;
        for (int i = 0; i < dim; i++) {
          double t = x / b[i].W - (double)b[i].n;
          if (fabs(t) >= 2.0) continue;
          v += rr[i] * hz_phi(t) * cexp(CMPLX(0.0, 1.0) * b[i].kx * x);
        }
        double complex u = amp * cexp(CMPLX(0.0, 1.0) * k * (x - xs));
        num += creal((v - u) * conj(v - u));
        den += creal(u * conj(u));
      }
      err = sqrt(num / den);
    }
    printf("  %-26s %8d %10.4f%s\n",
           carrier ? (lod ? "D carrier + LOD" : "C carrier + uniform fine")
                   : (lod ? "B plain + LOD  <-- control" : "A plain + uniform fine"),
           dim, err, (!carrier && lod) ? "   (MUST FAIL)" : "");
    free(G);
    free(r);
    free(sv);
  }
  printf("\nREAD: B must fail and D must not. If B succeeded, the carrier would be\n");
  printf("an optimisation rather than a condition, and the plan's central claim\n");
  printf("would be wrong. D at a small dim is the claim itself.\n");
  return 0;
}
