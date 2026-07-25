/* M9b — carrier (modulated) potentials on ANALYTIC assembly.
 *
 * Basis: two real quadrature families per node,
 *     B_n^c = phi(x/W - n) cos(kx),   B_n^s = phi(x/W - n) sin(kx)
 * Why coarse cells become legal: (d2/dx2 + k^2)[phi cos] = phi''/W^2 cos
 * - 2k phi'/W sin — the k^2 term cancels, which a plain smooth bump can never
 * produce.
 *
 * WHY THIS FILE WAS REWRITTEN (M9 audit, defect A1): the previous version
 * assembled by Gauss-Legendre quadrature over unit cells, i.e. it sampled the
 * oscillation at ~lambda/16 — a wavelength-sized grid sneaking back in through
 * the assembly. Its "4.2x fewer unknowns" therefore said nothing about COST.
 * Here every entry is closed-form: the carrier products expand by
 *     cos*cos = (1 + cos2kx)/2,  sin*sin = (1 - cos2kx)/2,  cos*sin = sin2kx/2
 * so each entry is a combination of  Int phi^(a) phi^(b) dx  and
 * Int phi^(a) phi^(b) e^{i 2k x} dx — both closed-form (M9a). The NUMBER OF
 * CALLS DOES NOT DEPEND ON k. The quadrature path stays as an independent
 * witness (cross-checked, same role GL5 plays in test_phi).
 *
 * Measures (PLAN.md M9b): cost quadruple (unknowns / structural nnz /
 * integral calls / wall time), error-vs-distance and error-vs-element-size
 * slices, seam cost between element sizes, conditioning vs element size. */
#include "helm1d.h"
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

enum { NCELL = 1024, PER_CELL = 8, GLN = 8, RAMP = 32, MAXDIM = 8192 };
static const double XSRC = 128.0;
static const double LAM = 16.0;
/* windows for the error profile: a growing profile means accumulated phase
 * (dispersion), a flat one means amplitude */
static const double DIST[4] = {5.0, 15.0, 30.0, 42.0};
/* rcond for zgelsd: the system is exactly rank-deficient by construction
 * (Nyquist null mode), 1e-10 is the value M2 established as working */
static const double RCOND = 1e-10;
/* Two medium cells join one segment when their k2 agree to rounding: the
 * medium is piecewise constant by construction, so anything below this is
 * float noise from a shared expression, not physics. */
static const double MEDSEG_REL_TOL = 1e-15;

typedef struct {
  double W;
  int n;
  int quad; /* 0 = cos, 1 = sin */
} basis1d;

typedef struct {
  double a, b;
  double complex dk2;
} medseg;

typedef struct {
  double wfine;    /* fine width, 0 = none */
  double wcoarse;  /* coarse width left of xseam */
  double wcoarse2; /* coarse width right of xseam, 0 = same as wcoarse */
  double xseam;
  int mixed; /* 1 = fine only in the fine zones */
} cfg;

typedef struct {
  int dim;
  long nnz;     /* structurally nonzero entries (overlapping supports) */
  long ncalls;  /* closed-form integral evaluations */
  double t_asm; /* seconds */
  double t_solve;
  double cond; /* sv[0]/sv[rank-1] over the numerically nonzero part */
  int rank;
  double err[4];
} runres;

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static int in_support(basis1d b, double x) {
  return fabs(x / b.W - (double)b.n) < 2.0;
}

static double bval(basis1d b, double x, double k) {
  if (!in_support(b, x)) return 0.0;
  double p = hz_phi(x / b.W - (double)b.n);
  return b.quad ? p * sin(k * x) : p * cos(k * x);
}

/* (d2/dx2 + k_bg^2) B — the background term cancels analytically */
static double bop(basis1d b, double x, double k) {
  double t = x / b.W - (double)b.n;
  double d1 = hz_phi_d1(t) / b.W, d2 = hz_phi_d2(t) / (b.W * b.W);
  double c = cos(k * x), s = sin(k * x);
  return b.quad ? (d2 * s + 2.0 * k * d1 * c) : (d2 * c - 2.0 * k * d1 * s);
}

static int in_fine_zone(double x) {
  if (x < (double)RAMP + 4.0) return 1;
  if (x > (double)(NCELL - RAMP) - 4.0) return 1;
  if (x > XSRC - LAM && x < XSRC + LAM) return 1;
  return 0;
}

static int build(basis1d *bs, cfg c) {
  int d = 0;
  double wc2 = c.wcoarse2 > 0.0 ? c.wcoarse2 : c.wcoarse;
  if (c.wcoarse > 0.0) {
    /* left of the seam at wcoarse, right of it at wcoarse2 */
    int n1 = (int)((double)NCELL / c.wcoarse) + 2;
    for (int n = -2; n <= n1; n++) {
      if ((double)n * c.wcoarse >= c.xseam) continue;
      for (int q = 0; q < 2; q++)
        bs[d++] = (basis1d){c.wcoarse, n, q};
    }
    int n2 = (int)((double)NCELL / wc2) + 2;
    for (int n = -2; n <= n2; n++) {
      if ((double)n * wc2 < c.xseam) continue;
      for (int q = 0; q < 2; q++)
        bs[d++] = (basis1d){wc2, n, q};
    }
  }
  if (c.wfine > 0.0) {
    int n1 = (int)((double)NCELL / c.wfine) + 2;
    for (int n = -2; n <= n1; n++) {
      if (c.mixed && !in_fine_zone((double)n * c.wfine)) continue;
      for (int q = 0; q < 2; q++)
        bs[d++] = (basis1d){c.wfine, n, q};
    }
  }
  return d;
}

/* Maximal intervals of constant dk2 = k2(x) - k_bg^2. The interior of this
 * scene is uniform (global absorption), so this collapses ~1024 cells to the
 * ramp cells plus one big segment — the medium term costs O(ramp), not O(N). */
static int build_segments(const hz_med1d *med, double kbg2, medseg *segs) {
  int ns = 0;
  double complex cur = med->k2[0] - kbg2;
  double start = 0.0;
  for (int c = 1; c < NCELL; c++) {
    double complex v = med->k2[c] - kbg2;
    if (cabs(v - cur) > MEDSEG_REL_TOL * kbg2) {
      segs[ns++] = (medseg){start, (double)c, cur};
      start = (double)c;
      cur = v;
    }
  }
  segs[ns++] = (medseg){start, (double)NCELL, cur};
  return ns;
}

/* Int over [x0,x1] of Phi_a^(da) * Phi_b^(db) * (trig_a * trig_b), closed form.
 * The trig product collapses to a DC part and a 2k part, nothing else. */
static double trigpair(double x0, double x1, hz_phi_factor fa, hz_phi_factor fb, int qa, int qb,
                       double k, long *ncalls) {
  if (x0 >= x1) return 0.0;
  *ncalls += 2;
  double i0 = hz_phi_prod_integral(x0, x1, fa, fb);
  double complex iw = hz_phi_prod_integral_osc(x0, x1, fa, fb, 2.0 * k);
  if (qa == qb) return qa == 0 ? 0.5 * (i0 + creal(iw)) : 0.5 * (i0 - creal(iw));
  return 0.5 * cimag(iw);
}

/* support of phi(x/W - n) is [W(n-2), W(n+2)], clipped to the domain */
static void overlap(basis1d bi, basis1d bj, double *lo, double *hi) {
  double l1 = bi.W * ((double)bi.n - 2.0), h1 = bi.W * ((double)bi.n + 2.0);
  double l2 = bj.W * ((double)bj.n - 2.0), h2 = bj.W * ((double)bj.n + 2.0);
  double l = l1 > l2 ? l1 : l2, h = h1 < h2 ? h1 : h2;
  if (l < 0.0) l = 0.0;
  if (h > (double)NCELL) h = (double)NCELL;
  *lo = l;
  *hi = h;
}

static double complex entry_analytic(basis1d bi, basis1d bj, double k, const medseg *segs, int nseg,
                                     long *ncalls) {
  double lo, hi;
  overlap(bi, bj, &lo, &hi);
  if (lo >= hi) return 0.0;
  hz_phi_factor fi0 = {bi.W, (double)bi.n, 0};
  hz_phi_factor fj0 = {bj.W, (double)bj.n, 0};
  hz_phi_factor fj1 = {bj.W, (double)bj.n, 1};
  hz_phi_factor fj2 = {bj.W, (double)bj.n, 2};

  /* (d2/dx2 + k_bg^2) B_j = Phi_j'' * t_j  +/-  2k Phi_j' * t_j^perp
   * (deriv factors already carry their 1/W^d, see hz_phi_factor) */
  double v = trigpair(lo, hi, fi0, fj2, bi.quad, bj.quad, k, ncalls);
  double sgn = bj.quad ? 1.0 : -1.0;
  v += sgn * 2.0 * k * trigpair(lo, hi, fi0, fj1, bi.quad, 1 - bj.quad, k, ncalls);

  double complex total = v;
  for (int s = 0; s < nseg; s++) {
    double a = lo > segs[s].a ? lo : segs[s].a;
    double b = hi < segs[s].b ? hi : segs[s].b;
    if (a >= b) continue;
    total += segs[s].dk2 * trigpair(a, b, fi0, fj0, bi.quad, bj.quad, k, ncalls);
  }
  return total;
}

static void assemble_analytic(const basis1d *bs, int dim, double k, const medseg *segs, int nseg,
                              double complex *A, double complex *b, long *nnz, long *ncalls) {
  hz_phi_factor fsrc = {1.0, XSRC, 0};
  for (int i = 0; i < dim; i++) {
    hz_phi_factor fi0 = {bs[i].W, (double)bs[i].n, 0};
    for (int j = 0; j < dim; j++) {
      double lo, hi;
      overlap(bs[i], bs[j], &lo, &hi);
      if (lo >= hi) continue;
      (*nnz)++;
      A[(size_t)i * (size_t)dim + (size_t)j] = entry_analytic(bs[i], bs[j], k, segs, nseg, ncalls);
    }
    double l1 = bs[i].W * ((double)bs[i].n - 2.0), h1 = bs[i].W * ((double)bs[i].n + 2.0);
    double lo = l1 > XSRC - 2.0 ? l1 : XSRC - 2.0, hi = h1 < XSRC + 2.0 ? h1 : XSRC + 2.0;
    if (lo < 0.0) lo = 0.0;
    if (hi > (double)NCELL) hi = (double)NCELL;
    if (lo < hi) {
      (*ncalls)++;
      double complex r = hz_phi_prod_integral_osc(lo, hi, fi0, fsrc, k);
      b[i] = bs[i].quad ? cimag(r) : creal(r);
    }
  }
}

/* Independent witness: the quadrature assembly this file used before M9b.
 * Kept to cross-check the closed forms, NOT used for the cost numbers. */
static void assemble_quad(const basis1d *bs, int dim, double k, const hz_med1d *med, double kbg2,
                          double complex *A, double complex *b) {
  static const double GX[GLN] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                                 -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                                 0.7966664774136267,  0.9602898564975363};
  static const double GW[GLN] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                                 0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                                 0.2223810344533745, 0.1012285362903763};
  for (int cell = 0; cell < NCELL; cell++) {
    double complex dk2 = med->k2[cell] - kbg2;
    for (int g = 0; g < GLN; g++) {
      double x = (double)cell + 0.5 + 0.5 * GX[g];
      double w = 0.5 * GW[g];
      for (int i = 0; i < dim; i++) {
        if (!in_support(bs[i], x)) continue;
        double vi = bval(bs[i], x, k);
        for (int j = 0; j < dim; j++) {
          if (!in_support(bs[j], x)) continue;
          A[(size_t)i * (size_t)dim + (size_t)j] +=
              w * vi * (bop(bs[j], x, k) + dk2 * bval(bs[j], x, k));
        }
        b[i] += w * vi * hz_phi(x - XSRC);
      }
    }
  }
}

static void field_errors(const basis1d *bs, int dim, const double complex *coef, double k,
                         const double complex *ufd, double *err) {
  for (int w = 0; w < 4; w++) {
    double c0 = XSRC + DIST[w] * LAM, c1 = c0 + 2.0 * LAM;
    double num = 0.0, den = 0.0;
    for (int p = 0; p <= NCELL * PER_CELL; p++) {
      double x = (double)p / (double)PER_CELL;
      if (x < c0 || x > c1) continue;
      double complex s = 0.0;
      for (int i = 0; i < dim; i++)
        s += coef[i] * bval(bs[i], x, k);
      double complex dd = s - ufd[p];
      num += creal(dd * conj(dd));
      den += creal(ufd[p] * conj(ufd[p]));
    }
    err[w] = sqrt(num / den);
  }
}

/* use_quad != 0 selects the witness path (cost numbers then meaningless) */
static int run(const char *tag, cfg c, const hz_med1d *med, const double complex *ufd, double k,
               const medseg *segs, int nseg, int use_quad, runres *out, int quiet) {
  basis1d *bs = calloc(MAXDIM, sizeof(basis1d));
  if (bs == NULL) return -1;
  int dim = build(bs, c);
  if (dim <= 0 || dim > MAXDIM) {
    free(bs);
    return -1;
  }
  double complex *A = calloc((size_t)dim * (size_t)dim, sizeof(double complex));
  double complex *b = calloc((size_t)dim, sizeof(double complex));
  double *sv = calloc((size_t)dim, sizeof(double));
  if (!A || !b || !sv) {
    free(A);
    free(b);
    free(sv);
    free(bs);
    return -1;
  }

  runres r = {0};
  r.dim = dim;
  double t0 = now_s();
  if (use_quad)
    assemble_quad(bs, dim, k, med, k * k, A, b);
  else
    assemble_analytic(bs, dim, k, segs, nseg, A, b, &r.nnz, &r.ncalls);
  r.t_asm = now_s() - t0;

  t0 = now_s();
  lapack_int rank = 0;
  lapack_int info = LAPACKE_zgelsd(LAPACK_ROW_MAJOR, dim, dim, 1, A, dim, b, 1, sv, RCOND, &rank);
  r.t_solve = now_s() - t0;
  r.rank = (int)rank;
  r.cond = (rank > 0 && sv[rank - 1] > 0.0) ? sv[0] / sv[rank - 1] : -1.0;
  for (int w = 0; w < 4; w++)
    r.err[w] = -1.0;
  if (info == 0) field_errors(bs, dim, b, k, ufd, r.err);

  if (!quiet) {
    if (use_quad)
      printf("  %-24s dim=%-5d [quadrature witness]        err@5=%.3f @15=%.3f @30=%.3f @42=%.3f\n",
             tag, dim, r.err[0], r.err[1], r.err[2], r.err[3]);
    else
      printf("  %-24s dim=%-5d nnz/row=%-6.1f calls=%-9ld asm=%6.2fs solve=%6.2fs  "
             "rank=%d/%d cond=%.2e  err@5=%.3f @15=%.3f @30=%.3f @42=%.3f\n",
             tag, dim, (double)r.nnz / (double)dim, r.ncalls, r.t_asm, r.t_solve, r.rank, dim,
             r.cond, r.err[0], r.err[1], r.err[2], r.err[3]);
  }
  if (out) *out = r;
  free(A);
  free(b);
  free(sv);
  free(bs);
  return 0;
}

int main(void) {
  double k = 2.0 * M_PI / LAM;
  hz_med1d med;
  hz_med1d_init(&med, NCELL, k, 0.02, RAMP);
  double complex *ufd = calloc((size_t)NCELL * PER_CELL + 1, sizeof(double complex));
  medseg *segs = calloc(NCELL + 1, sizeof(medseg));
  if (ufd == NULL || segs == NULL || hz_fd_reference(&med, PER_CELL, XSRC, ufd) != 0) {
    free(ufd);
    free(segs);
    hz_med1d_free(&med);
    return 1;
  }
  int nseg = build_segments(&med, k * k, segs);

  printf("M9b: carrier basis on ANALYTIC assembly\n");
  printf("domain %d cells = %.0f lambda, source at %.0f, windows at %.0f/%.0f/%.0f/%.0f lambda\n",
         NCELL, (double)NCELL / LAM, XSRC, DIST[0], DIST[1], DIST[2], DIST[3]);
  printf("medium segments of constant dk2: %d (vs %d cells) — ramps + one interior\n\n", nseg,
         NCELL);

  /* --- 1. agreement with the quadrature witness ------------------------- */
  printf("[1] analytic vs quadrature witness (must agree to 1e-3)\n");
  cfg small = {2.0, 32.0, 0.0, 0.0, 1};
  runres ra, rq;
  if (run("analytic", small, &med, ufd, k, segs, nseg, 0, &ra, 0) != 0) return 1;
  if (run("quadrature", small, &med, ufd, k, segs, nseg, 1, &rq, 0) != 0) return 1;
  double dmax = 0.0;
  for (int w = 0; w < 4; w++) {
    double d = fabs(ra.err[w] - rq.err[w]);
    if (d > dmax) dmax = d;
  }
  printf("  max |err_analytic - err_quad| = %.2e  -> %s\n\n", dmax,
         dmax < 1e-3 ? "AGREE" : "MISMATCH");

  /* --- 2. cost quadruple + error vs element size ------------------------ */
  printf("[2] cost quadruple and error vs element size (fine = lambda/8 baseline)\n");
  cfg fine = {2.0, 0.0, 0.0, 0.0, 0};
  runres rfine;
  if (run("FINE uniform lambda/8", fine, &med, ufd, k, segs, nseg, 0, &rfine, 0) != 0) return 1;
  runres rmix[4];
  for (int i = 0; i < 4; i++) {
    double wc = LAM * (double)(1 << i);
    char t[64];
    snprintf(t, sizeof(t), "MIXED coarse %.0f lambda", wc / LAM);
    cfg c = {2.0, wc, 0.0, 0.0, 1};
    if (run(t, c, &med, ufd, k, segs, nseg, 0, &rmix[i], 0) != 0) return 1;
  }
  printf("\n  cost ratio vs FINE (lower is better):\n");
  for (int i = 0; i < 4; i++)
    printf("    coarse %.0f lambda: dim x%.2f  calls x%.3f  asm x%.3f  solve x%.4f\n",
           (double)(1 << i), (double)rmix[i].dim / (double)rfine.dim,
           (double)rmix[i].ncalls / (double)rfine.ncalls, rmix[i].t_asm / rfine.t_asm,
           rmix[i].t_solve / rfine.t_solve);

  /* --- 3. seam cost between element sizes (input to PLAN question 3) ---- */
  printf("\n[3] seam cost: element size changes mid free-propagation at x=%.0f\n",
         (double)NCELL / 2.0);
  printf("    near zone stays 1 lambda, far zone is coarsened; the reference\n");
  printf("    seamless runs are in [2] above (1 lambda: %.3f, 2 lambda: %.3f @42)\n",
         rmix[0].err[3], rmix[1].err[3]);
  runres rseam[4];
  for (int i = 0; i < 4; i++) {
    double wfar = LAM * (double)(1 << (i + 1));
    char t[64];
    snprintf(t, sizeof(t), "SEAM 1 -> %.0f lambda", wfar / LAM);
    cfg c = {2.0, LAM, wfar, (double)NCELL / 2.0, 1};
    if (run(t, c, &med, ufd, k, segs, nseg, 0, &rseam[i], 0) != 0) return 1;
  }
  printf("\n  far-zone coarsening, error @30/@42 (both downstream of the seam):\n");
  printf("    seamless 1 lambda everywhere : %.3f / %.3f  (dim %d)\n", rmix[0].err[2],
         rmix[0].err[3], rmix[0].dim);
  for (int i = 0; i < 4; i++)
    printf("    far %2.0f lambda                : %.3f / %.3f  (dim %d)\n", (double)(1 << (i + 1)),
           rseam[i].err[2], rseam[i].err[3], rseam[i].dim);
  printf("  reading: if far-zone coarsening is ~free while UNIFORM coarsening is not,\n");
  printf("  the element-size limit comes from STRUCTURE (the source), not from transport.\n");

  free(ufd);
  free(segs);
  hz_med1d_free(&med);
  return 0;
}
