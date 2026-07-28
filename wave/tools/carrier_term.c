/* PLAN question 1: HOW TO TERMINATE THE DOMAIN WHEN ELEMENTS ARE >> LAMBDA.
 *
 * This question was deferred four times and spoiled four measurements in a row;
 * it is now on the critical path. Three schemes were listed: an absorbing layer
 * (PML-like), a DtN / exact radiation condition, and Green's-function elements.
 *
 * HONEST SCOPE, STATED UP FRONT: in 1D the third scheme DEGENERATES INTO THE
 * SECOND. Beyond the last scatterer the exterior solution is exactly one
 * outgoing plane wave, so "represent the exterior analytically" and "impose
 * u' = i k u at the boundary" are the same statement. They part company only in
 * 2D/3D, where the exterior is not a single plane wave and DtN becomes
 * non-local. So this bench measures TWO schemes and says so, rather than
 * pretending to three.
 *
 * WHAT IS ACTUALLY NEW HERE. The absorbing layer used to be unusable with large
 * elements: its decay length is ~lambda whatever its thickness, so an element
 * of hundreds of wavelengths saw a step, not a ramp, and the boundary simply
 * reflected. With the carrier now COMPLEX (see tests/test_carrier_op.c) the
 * layer's own decaying wave is exactly in the span, so the layer should be
 * representable at ANY element size. That is the measurement.
 *
 * Reference is the analytic outgoing Green's function of the uniform medium —
 * exact at any distance, no grid, valid at any lambda (M9 audit rule A4). */
#include "carrier.h"
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum { MAXB = 1200, MAXSEG = 64, NSAMP = 2000 };
static const double LAM = 16.0;
static const double RCOND = 1e-12;
/* The absorbing scheme gets its own, much smaller threshold on purpose. The
 * radiation condition enters an absorbing formulation ONLY through a
 * near-null direction of the operator (a plane wave is annihilated exactly in
 * the lossless interior and merely ALMOST annihilated once the layer damps it),
 * so a min-norm solver that truncates small singular values throws the
 * radiation condition away and returns the silent solution. That is the same
 * "the minimiser prefers silence" mechanism M2 found for LSQ formulations. If
 * the diagnosis is right, lowering the threshold must revive the field. */
static const double RCOND_ABS = 1e-16;
/* Absorbing layer: NSUB geometrically graded sub-layers. Graded because a step
 * in k^2 reflects; geometric because that is what keeps each sub-layer's
 * contrast modest while the total absorption grows fast. */
enum { NSUB = 6 };
/* The layer is scaled to a fixed OPTICAL DEPTH, not to a fixed sigma. A fixed
 * sigma is what a lambda-sized layer wants, but here the layer is sized in
 * ELEMENTS, so at 200 lambda per element it is 1200 lambda thick and a sigma of
 * order 1 decays by e^-3700 — pure underflow, and the solve simply fails. What
 * absorption needs is total depth Int Im(k) dx of a few, whatever the
 * thickness. With the quartic profile Int sigma dx = sigma_max*layer/5 and
 * Im(k) ~ k sigma/2, so sigma_max = 10*DEPTH/(k*layer). */
static const double DEPTH = 4.0;

typedef struct {
  double a, b;
  double complex k2;
  double complex kx; /* carrier matched to THIS region */
} region;

/* Int phi(x/W - n) * phi(y source) e^{i kx x}: right-hand side for a source
 * f(x) = phi(x - xsrc) of unit width. */
static double complex rhs_src(hz_carrier b, double xsrc) {
  hz_phi_factor fi = {b.W, (double)b.n, 0}, fs = {1.0, xsrc, 0};
  double lo = b.W * ((double)b.n - 2.0), hi = b.W * ((double)b.n + 2.0);
  if (lo < xsrc - 2.0) lo = xsrc - 2.0;
  if (hi > xsrc + 2.0) hi = xsrc + 2.0;
  if (lo >= hi) return 0.0;
  return hz_phi_prod_integral_osc(lo, hi, fi, fs, b.kx);
}

int main(void) {
  double k = 2.0 * M_PI / LAM;
  double complex ii = CMPLX(0.0, 1.0);

  printf("Domain termination with elements >> lambda (PLAN question 1)\n");
  printf("NB: in 1D the Green's-function scheme degenerates into DtN — the\n");
  printf("exterior is one outgoing plane wave. Two schemes measured, not three.\n\n");
  printf("  %8s %8s %6s %14s %14s\n", "W/lam", "layer/lam", "dim", "err ABSORBING", "err DtN");

  static const double WL[5] = {0.5, 2.0, 10.0, 50.0, 200.0};
  for (int iw = 0; iw < 5; iw++) {
    double W = WL[iw] * LAM;
    /* domain sized in elements so the comparison is like-for-like */
    double core = 40.0 * W; /* physical region */
    double layer = 6.0 * W; /* absorbing layer, one sub-layer per element */
    double dom = core + 2.0 * layer;
    double xsrc = 0.35 * core + layer;

    /* regions: graded layer | core | graded layer, each with its own carrier */
    static region rg[MAXSEG];
    int nr = 0;
    for (int s = 0; s < NSUB; s++) {
      double f0 = (double)s / (double)NSUB, f1 = (double)(s + 1) / (double)NSUB;
      double sig_max = 10.0 * DEPTH / (k * layer);
      double sg = sig_max * pow(1.0 - 0.5 * (f0 + f1), 4.0); /* quartic, peak at the wall */
      double complex k2 = CMPLX(k * k, k * k * sg);
      rg[nr++] = (region){f0 * layer, f1 * layer, k2, csqrt(k2)};
    }
    rg[nr++] = (region){layer, dom - layer, CMPLX(k * k, 0.0), CMPLX(k, 0.0)};
    for (int s = 0; s < NSUB; s++) {
      double f0 = (double)s / (double)NSUB, f1 = (double)(s + 1) / (double)NSUB;
      double sig_max = 10.0 * DEPTH / (k * layer);
      double sg = sig_max * pow(0.5 * (f0 + f1), 4.0);
      double complex k2 = CMPLX(k * k, k * k * sg);
      rg[nr++] = (region){dom - layer + f0 * layer, dom - layer + f1 * layer, k2, csqrt(k2)};
    }

    static hz_medseg segs[MAXSEG];
    for (int s = 0; s < nr; s++)
      segs[s] = (hz_medseg){rg[s].a, rg[s].b, rg[s].k2};

    /* basis: two directions per node; the carrier of a node is that of the
     * region its CENTRE falls in — with the layer graded per element, an
     * element never spans two regions by construction */
    static hz_carrier b[MAXB];
    int dim = 0;
    int n1 = (int)(dom / W) + 3;
    for (int n = -3; n <= n1 && dim + 2 <= MAXB; n++) {
      double xc = (double)n * W;
      double complex kc = CMPLX(k, 0.0);
      for (int s = 0; s < nr; s++)
        if (xc >= rg[s].a && xc < rg[s].b) kc = rg[s].kx;
      b[dim++] = (hz_carrier){W, n, kc, -1};
      b[dim++] = (hz_carrier){W, n, -kc, -1};
    }

    /* two systems: rows = Galerkin equations; DtN adds two constraint rows */
    for (int scheme = 0; scheme < 2; scheme++) {
      int nrow = dim + (scheme == 1 ? 2 : 0);
      double complex *A = calloc((size_t)nrow * (size_t)dim, sizeof(double complex));
      double complex *r = calloc((size_t)(nrow > dim ? nrow : dim), sizeof(double complex));
      double *sv = calloc((size_t)dim, sizeof(double));
      if (!A || !r || !sv) return 1;

      /* scheme 0: absorbing layer everywhere, no boundary condition.
       * scheme 1: NO absorption at all (lossless), DtN rows instead. */
      int nseg_use = scheme == 0 ? nr : 1;
      hz_medseg vac = {0.0, dom, CMPLX(k * k, 0.0)};
      const hz_medseg *sg_use = scheme == 0 ? segs : &vac;

      for (int i = 0; i < dim; i++) {
        hz_carrier bi = b[i];
        if (scheme == 1) bi.kx = creal(bi.kx) > 0.0 ? CMPLX(k, 0.0) : CMPLX(-k, 0.0);
        for (int j = 0; j < dim; j++) {
          hz_carrier bj = b[j];
          if (scheme == 1) bj.kx = creal(bj.kx) > 0.0 ? CMPLX(k, 0.0) : CMPLX(-k, 0.0);
          A[(size_t)i * (size_t)dim + (size_t)j] =
              hz_carrier_entry(bi, bj, dom, 0.0, sg_use, nseg_use, NULL);
        }
        r[i] = rhs_src(bi, xsrc);
      }
      if (scheme == 1) {
        /* exact 1D radiation condition: u' = +ik u at the right end,
         * u' = -ik u at the left. Added as two least-squares rows, scaled to
         * the magnitude of the operator rows so they are neither ignored nor
         * dominant. */
        double scale = 0.0;
        for (int j = 0; j < dim; j++)
          scale += cabs(A[(size_t)(dim / 2) * (size_t)dim + (size_t)j]);
        if (scale <= 0.0) scale = 1.0;
        for (int e = 0; e < 2; e++) {
          double xb = (e == 0) ? 0.0 : dom;
          double complex sgn = (e == 0) ? -ii * k : ii * k;
          for (int j = 0; j < dim; j++) {
            hz_carrier bj = b[j];
            bj.kx = creal(bj.kx) > 0.0 ? CMPLX(k, 0.0) : CMPLX(-k, 0.0);
            double t = xb / bj.W - (double)bj.n;
            double complex val = 0.0, der = 0.0;
            if (fabs(t) < 2.0) {
              double complex e_i = cexp(ii * bj.kx * xb);
              val = hz_phi(t) * e_i;
              der = (hz_phi_d1(t) / bj.W + ii * bj.kx * hz_phi(t)) * e_i;
            }
            A[(size_t)(dim + e) * (size_t)dim + (size_t)j] = scale * (der - sgn * val);
          }
          r[dim + e] = 0.0;
        }
      }

      lapack_int rank = 0;
      double rc = (scheme == 0) ? RCOND_ABS : RCOND;
      lapack_int info = LAPACKE_zgelsd(LAPACK_ROW_MAJOR, nrow, dim, 1, A, dim, r, 1, sv, rc, &rank);

      /* compare with the analytic outgoing Green's function in the core */
      double num = 0.0, den = 0.0;
      if (info == 0) {
        double complex ph = 0.0;
        {
          /* phihat(k): Int phi(t) e^{-ikt} dt, 8-point GL per unit piece */
          static const double GX[8] = {
              -0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
              0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363};
          static const double GW[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                                       0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                                       0.2223810344533745, 0.1012285362903763};
          for (int p = -2; p < 2; p++)
            for (int g = 0; g < 8; g++) {
              double t = (double)p + 0.5 + 0.5 * GX[g];
              ph += 0.5 * GW[g] * hz_phi(t) * cexp(-ii * k * t);
            }
        }
        double x0 = xsrc + 10.0 * LAM, x1 = dom - layer - 10.0 * LAM;
        if (x1 > x0)
          for (int p = 0; p < NSAMP; p++) {
            double x = x0 + (x1 - x0) * ((double)p + 0.5) / (double)NSAMP;
            double complex v = 0.0;
            for (int j = 0; j < dim; j++) {
              hz_carrier bj = b[j];
              if (scheme == 1) bj.kx = creal(bj.kx) > 0.0 ? CMPLX(k, 0.0) : CMPLX(-k, 0.0);
              v += r[j] * hz_carrier_val(bj, x, 0.0);
            }
            double complex u = ph * cexp(ii * k * (x - xsrc)) / (2.0 * ii * k);
            num += creal((v - u) * conj(v - u));
            den += creal(u * conj(u));
          }
      }
      double err = den > 0.0 ? sqrt(num / den) : -1.0;
      if (scheme == 0)
        printf("  %8.1f %8.1f %6d %11.3e (rank %3d)", WL[iw], layer / LAM, dim, err, (int)rank);
      else
        printf(" %11.3e (rank %3d)\n", err, (int)rank);
      free(A);
      free(r);
      free(sv);
    }
  }
  printf("\nREAD: the absorbing layer used to be unusable past a few lambda per\n");
  printf("element because its decay length is ~lambda however thick it is. With a\n");
  printf("COMPLEX carrier the layer's own decaying wave is exactly in the span, so\n");
  printf("the layer should now work at any element size. DtN is exact in 1D and\n");
  printf("serves as the yardstick.\n");
  return 0;
}
