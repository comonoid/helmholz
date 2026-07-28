/* PLAN question 4, first half: DOES THE CASCADE REMOVE INTER-LEVEL REDUNDANCY?
 *
 * The plan has asserted since day one that "the cascade takes out the
 * inter-level redundancy, so no preconditioner is needed". That has never been
 * tested with the carrier basis — and the one time all levels were thrown into
 * a single solve, the system went to cond 1e10 with 27 levels and the accuracy
 * collapsed. Every solver bench since has been measuring THE COST OF NOT HAVING
 * A CASCADE rather than any property of the method.
 *
 * Tested here in the PROJECTION setting: no operator, no domain boundary, no
 * absorbing layer — none of the six bench artefacts this project has produced
 * can reach in. If the cascade fails to remove the redundancy even in plain
 * approximation, it certainly will not in a solve. Necessary condition, cheap.
 *
 * MONOLITHIC: fit all levels at once, min-norm.
 * CASCADIC:   fit the coarsest level; subtract; fit the next level to what is
 *             left; and so on. Each level then faces a SINGLE-LEVEL system,
 *             which is only mildly degenerate (the Nyquist null mode of phi
 *             translates) instead of massively so.
 * SWEEPS:     E3 measured in 1D that one coarse-to-fine pass does not deliver a
 *             local refinement to the sensor — the coarse level was solved
 *             before the fine one existed and still carries the old radiation.
 *             So a second pass is checked here too, as standard rather than
 *             optional.
 *
 * The target field is multi-scale on purpose: a carrier plus beats at 1e-1 and
 * 1e-2 of k, the finer beats confined near the origin. Its envelope therefore
 * has structure at ~10 lambda near the centre and ~100 lambda further out,
 * which is exactly what a shell ladder is supposed to capture level by level.
 * Being a piecewise sum of plane waves, every integral stays closed-form. */
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum { MAXB = 900, MAXSEG = 16, NLEV = 6 };
static const double LAM = 16.0;
static const double W0 = 2.0;    /* finest element = lambda/8 */
static const double EPS = 0.125; /* shell rule: overlap 1.0 (SHELLS_REPORT) */
static const double RCOND = 1e-12;

typedef struct {
  double x0, x1;
  double complex amp;
  double kap;
} pw;

typedef struct {
  double W;
  int n;
  double kx;
  int lev;
} bfn;

/* <B_i, B_j> = Int phi_i phi_j e^{i(kx_j - kx_i)x} dx (conjugated: L2) */
static double complex gram(bfn a, bfn b) {
  hz_phi_factor fa = {a.W, (double)a.n, 0}, fb = {b.W, (double)b.n, 0};
  double lo = a.W * ((double)a.n - 2.0), hi = a.W * ((double)a.n + 2.0);
  double lo2 = b.W * ((double)b.n - 2.0), hi2 = b.W * ((double)b.n + 2.0);
  if (lo2 > lo) lo = lo2;
  if (hi2 < hi) hi = hi2;
  if (lo >= hi) return 0.0;
  return hz_phi_prod_integral_osc(lo, hi, fa, fb, b.kx - a.kx);
}

/* <u, B_i> for a piecewise plane-wave field */
static double complex rhs_of(const pw *seg, int ns, bfn a) {
  hz_phi_factor fa = {a.W, (double)a.n, 0};
  double lo = a.W * ((double)a.n - 2.0), hi = a.W * ((double)a.n + 2.0);
  double complex t = 0.0;
  for (int s = 0; s < ns; s++) {
    double l = lo > seg[s].x0 ? lo : seg[s].x0;
    double h = hi < seg[s].x1 ? hi : seg[s].x1;
    if (l >= h) continue;
    t += seg[s].amp * hz_phi_integral_osc(l, h, fa, seg[s].kap - a.kx);
  }
  return t;
}

static double field_norm2(const pw *seg, int ns, double a, double b) {
  double t = 0.0;
  for (int p = 0; p < ns; p++)
    for (int q = 0; q < ns; q++) {
      double l = a > seg[p].x0 ? a : seg[p].x0, h = b < seg[p].x1 ? b : seg[p].x1;
      double l2 = seg[q].x0, h2 = seg[q].x1;
      if (l2 > l) l = l2;
      if (h2 < h) h = h2;
      if (l >= h) continue;
      double w = seg[p].kap - seg[q].kap;
      double complex ii;
      if (fabs(w) < 1e-14)
        ii = h - l;
      else
        ii = (CMPLX(cos(w * h), sin(w * h)) - CMPLX(cos(w * l), sin(w * l))) / CMPLX(0.0, w);
      t += creal(seg[p].amp * conj(seg[q].amp) * ii);
    }
  return t;
}

/* min-norm solve of G c = r, returns rank */
static int solve(double complex *G, double complex *r, int n, double *sv) {
  lapack_int rank = 0;
  if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, n, n, 1, G, n, r, 1, sv, RCOND, &rank) != 0) return -1;
  return (int)rank;
}

int main(void) {
  double k = 2.0 * M_PI / LAM;
  double dom = 4000.0 * LAM;
  double xs = 0.5 * dom; /* the structure sits mid-domain */

  /* MULTI-SCALE TARGET, AND IT MUST BE A PHYSICAL FIELD.
   * The first version of this bench used components at 1.1k and 1.01k. That was
   * wrong: this is a TREFFTZ basis, its functions satisfy the homogeneous
   * equation locally, and a plane wave with |kappa| != k does not solve it at
   * all. In a homogeneous medium every physical field is a superposition of
   * plane waves of FIXED MAGNITUDE differing only in DIRECTION — in 1D that is
   * exactly +k and -k, both of which the basis already holds. The 22% error the
   * first run reported was therefore the price of an unphysical target, not a
   * property of the cascade. (It also kills the idea of "letting the carrier
   * magnitude be fitted": in a piecewise-constant medium the magnitude is fixed
   * by the medium. Only for gradient media, via the WKB phase, does it vary.)
   *
   * So: a right-going wave everywhere, plus a left-going one whose amplitude is
   * piecewise constant on DYADIC annuli around the centre. Each shell of the
   * ladder then meets a feature at its own scale, which is what the ladder is
   * for, and every piece is still a sum of +-k plane waves, so the integrals
   * stay closed-form. */
  pw seg[MAXSEG];
  int ns = 0;
  seg[ns++] = (pw){0.0, dom, 1.0, k};
  double sig = 25.0 * LAM;
  for (int m = 0; m < 5 && ns + 2 < MAXSEG; m++) {
    double r_in = (m == 0) ? 0.0 : sig * pow(2.0, (double)(m - 1));
    double r_out = sig * pow(2.0, (double)m);
    double amp = (m % 2 == 0) ? 0.5 : -0.3;
    seg[ns++] = (pw){xs - r_out, xs - r_in, amp, -k};
    seg[ns++] = (pw){xs + r_in, xs + r_out, amp, -k};
  }

  static bfn b[MAXB];
  static int lev0[NLEV + 1];
  int dim = 0;
  for (int j = 0; j < NLEV; j++) {
    lev0[j] = dim;
    double W = W0 * pow(2.0, (double)(NLEV - 1 - j)) * 16.0; /* coarse first */
    double alo = (j == NLEV - 1) ? 0.0 : W / EPS;
    double ahi = (j == 0) ? 2.0 * dom : 2.0 * W / EPS;
    alo -= 2.0 * W;
    ahi += 2.0 * W;
    if (alo < 0.0) alo = 0.0;
    int n0 = (int)((xs - ahi) / W) - 2, n1 = (int)((xs + ahi) / W) + 2;
    for (int n = n0; n <= n1 && dim + 2 <= MAXB; n++) {
      double xc = (double)n * W;
      if (xc < 0.0 || xc > dom) continue;
      double d = fabs(xc - xs);
      if (d < alo || d >= ahi) continue;
      b[dim++] = (bfn){W, n, k, j};
      b[dim++] = (bfn){W, n, -k, j};
    }
  }
  lev0[NLEV] = dim;

  double nu = field_norm2(seg, ns, 0.0, dom);
  printf("Cascade vs monolithic, projection setting (no operator, no boundary)\n");
  printf("domain %.0f lambda, %d levels, dim %d, |u|^2 = %.4e\n\n", dom / LAM, NLEV, dim, nu);

  /* ---------- monolithic ---------- */
  {
    double complex *G = calloc((size_t)dim * (size_t)dim, sizeof(double complex));
    double complex *r = calloc((size_t)dim, sizeof(double complex));
    double *sv = calloc((size_t)dim, sizeof(double));
    if (!G || !r || !sv) return 1;
    for (int i = 0; i < dim; i++) {
      for (int j = 0; j < dim; j++)
        G[(size_t)i * (size_t)dim + (size_t)j] = gram(b[i], b[j]);
      r[i] = rhs_of(seg, ns, b[i]);
    }
    static double complex rr[MAXB];
    for (int i = 0; i < dim; i++)
      rr[i] = r[i];
    int rank = solve(G, rr, dim, sv);
    /* error: |u|^2 - 2Re<u,v> + <v,v>, all coefficients now in rr */
    double uv = 0.0, vv = 0.0;
    for (int i = 0; i < dim; i++) {
      uv += creal(conj(rr[i]) * r[i]);
      for (int j = 0; j < dim; j++)
        vv += creal(rr[i] * conj(rr[j]) * gram(b[i], b[j]));
    }
    double e2 = nu - 2.0 * uv + vv;
    if (e2 < 0.0) e2 = 0.0;
    printf("[monolithic] dim %d  rank %d (%.0f%%)  rel err %.4e\n", dim, rank,
           100.0 * (double)rank / (double)dim, sqrt(e2 / nu));
    free(G);
    free(r);
    free(sv);
  }

  /* ---------- cascadic, with a second sweep ---------- */
  static double complex coef[MAXB];
  for (int i = 0; i < dim; i++)
    coef[i] = 0.0;
  for (int sweep = 0; sweep < 2; sweep++) {
    printf("\n[cascadic] sweep %d\n", sweep + 1);
    for (int j = 0; j < NLEV; j++) {
      int a = lev0[j], nlev = lev0[j + 1] - a;
      if (nlev <= 0) continue;
      double complex *G = calloc((size_t)nlev * (size_t)nlev, sizeof(double complex));
      double complex *r = calloc((size_t)nlev, sizeof(double complex));
      double *sv = calloc((size_t)nlev, sizeof(double));
      if (!G || !r || !sv) return 1;
      for (int i = 0; i < nlev; i++) {
        for (int l = 0; l < nlev; l++)
          G[(size_t)i * (size_t)nlev + (size_t)l] = gram(b[a + i], b[a + l]);
        /* residual right-hand side: the field minus what all OTHER levels
         * currently hold. On sweep 2 that includes the finer levels, which is
         * exactly what E3 said a single pass cannot do. */
        double complex t = rhs_of(seg, ns, b[a + i]);
        for (int m = 0; m < dim; m++) {
          if (m >= a && m < a + nlev) continue;
          if (!(cabs(coef[m]) > 0.0)) continue;
          t -= coef[m] * gram(b[a + i], b[m]);
        }
        r[i] = t;
      }
      int rank = solve(G, r, nlev, sv);
      for (int i = 0; i < nlev; i++)
        coef[a + i] = r[i];
      printf("  level %d: dim %4d  rank %4d (%3.0f%%)\n", j, nlev, rank,
             100.0 * (double)rank / (double)nlev);
      free(G);
      free(r);
      free(sv);
    }
    double uv = 0.0, vv = 0.0;
    for (int i = 0; i < dim; i++) {
      uv += creal(conj(coef[i]) * rhs_of(seg, ns, b[i]));
      for (int j = 0; j < dim; j++)
        vv += creal(coef[i] * conj(coef[j]) * gram(b[i], b[j]));
    }
    double e2 = nu - 2.0 * uv + vv;
    if (e2 < 0.0) e2 = 0.0;
    printf("  rel err after sweep %d: %.4e\n", sweep + 1, sqrt(e2 / nu));
  }
  printf("\nREAD: per-level rank near 100%% means the cascade really does remove\n");
  printf("the inter-level redundancy; the monolithic rank is the price of not\n");
  printf("doing it. Sweep 2 vs 1 is the E3 effect.\n");
  return 0;
}
