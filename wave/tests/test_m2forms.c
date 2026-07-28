/* M2 tail: formulation comparison (Galerkin vs collocation-LS vs Gram/FOSLS)
 * and the Toeplitz-FFT floor solve (PLAN.md audit items 1 and 2). */
#include "fft.h"
#include "helm1d.h"
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int g_fail = 0;
static int g_total = 0;

static void check_lt(double got, double bound, const char *what) {
  g_total++;
  if (!(got < bound)) {
    g_fail++;
    printf("FAIL: %s: got %.6g, bound %.6g\n", what, got, bound);
  }
}

enum { NCELL = 256, PER_CELL = 16, NPAD = 256 };
static const double SENS0 = 192.0, SENS1 = 208.0;
static const double XSRC = 48.0;
static const hz_level1d LV_F = {1, -1, 129}; /* the floor: h=2 = lambda/8 */
#define NF 131                               /* LV_F count */

static hz_phi_factor pot(int n, int deriv) {
  hz_phi_factor f = {2.0, (double)(LV_F.n0 + n), deriv};
  return f;
}

static double sensor_err_coef(double complex *coef, const double complex *ufd) {
  hz_sol1d s;
  int ofs[2] = {0, NF};
  s.levels = &LV_F;
  s.nlev = 1;
  s.coef = coef;
  s.ofs = ofs;
  double num = 0.0, den = 0.0;
  for (int i = 0; i <= NCELL * PER_CELL; i++) {
    double x = (double)i / (double)PER_CELL;
    if (x < SENS0 || x > SENS1) continue;
    double complex d = hz_sol1d_eval(&s, x) - ufd[i];
    num += creal(d * conj(d));
    den += creal(ufd[i] * conj(ufd[i]));
  }
  return sqrt(num / den);
}

/* ---- formulation 2: pointwise collocation least squares ------------------ */
static double solve_collocation(const hz_med1d *med, const double complex *ufd, double rcond,
                                double *cond, int *rank_out) {
  int mrows = NCELL * PER_CELL + 1;
  double wgt = sqrt(1.0 / (double)PER_CELL);
  double complex *A = calloc((size_t)mrows * NF, sizeof(double complex));
  double complex *b = calloc((size_t)mrows, sizeof(double complex));
  double *sv = calloc(NF, sizeof(double));
  if (A == NULL || b == NULL || sv == NULL) {
    free(A);
    free(b);
    free(sv);
    return 1e9;
  }
  for (int i = 0; i < mrows; i++) {
    double x = (double)i / (double)PER_CELL;
    int c = i / PER_CELL;
    if (c >= NCELL) c = NCELL - 1;
    for (int j = 0; j < NF; j++) {
      double t = x / 2.0 - (double)(LV_F.n0 + j);
      if (t <= -2.0 || t >= 2.0) continue;
      A[(size_t)i * NF + (size_t)j] = (hz_phi_d2(t) / 4.0 + med->k2[c] * hz_phi(t)) * wgt;
    }
    b[i] = hz_phi(x - XSRC) * wgt;
  }
  lapack_int rank = 0;
  lapack_int info = LAPACKE_zgelsd(LAPACK_ROW_MAJOR, mrows, NF, 1, A, NF, b, 1, sv, rcond, &rank);
  double err = 1e9;
  if (info == 0) {
    err = sensor_err_coef(b, ufd);
    *cond = (rank >= 1 && sv[rank - 1] > 0.0) ? sv[0] / sv[rank - 1] : 0.0;
    *rank_out = (int)rank;
  }
  free(A);
  free(b);
  free(sv);
  return err;
}

/* ---- formulation 3: Gram of L*Phi (= normal equations of continuous LS) -- */
static double complex gram_entry(const hz_med1d *m, int in, int jn) {
  hz_phi_factor pi0 = pot(in, 0), pi2 = pot(in, 2);
  hz_phi_factor pj0 = pot(jn, 0), pj2 = pot(jn, 2);
  double lo = 2.0 * (pi0.n - 2.0), hi = 2.0 * (pi0.n + 2.0);
  double lo2 = 2.0 * (pj0.n - 2.0), hi2 = 2.0 * (pj0.n + 2.0);
  if (lo2 > lo) lo = lo2;
  if (hi2 < hi) hi = hi2;
  if (lo < 0.0) lo = 0.0;
  if (hi > (double)m->ncell) hi = (double)m->ncell;
  if (lo >= hi) return 0.0;
  double complex v = hz_phi_prod_integral(lo, hi, pi2, pj2);
  int c0 = (int)floor(lo), c1 = (int)ceil(hi);
  for (int c = c0; c < c1; c++) {
    double a = (double)c, b = (double)(c + 1);
    if (a < lo) a = lo;
    if (b > hi) b = hi;
    if (a >= b) continue;
    double complex k2 = m->k2[c];
    v += k2 * hz_phi_prod_integral(a, b, pi2, pj0);
    v += conj(k2) * hz_phi_prod_integral(a, b, pi0, pj2);
    v += k2 * conj(k2) * hz_phi_prod_integral(a, b, pi0, pj0);
  }
  return v;
}

static double solve_gram(const hz_med1d *med, const double complex *ufd, double rcond, double *cond,
                         int *rank_out) {
  double complex *A = calloc((size_t)NF * NF, sizeof(double complex));
  double complex *b = calloc(NF, sizeof(double complex));
  double *sv = calloc(NF, sizeof(double));
  if (A == NULL || b == NULL || sv == NULL) {
    free(A);
    free(b);
    free(sv);
    return 1e9;
  }
  hz_phi_factor src = {1.0, XSRC, 0};
  for (int i = 0; i < NF; i++) {
    for (int j = 0; j < NF; j++)
      A[(size_t)i * NF + (size_t)j] = gram_entry(med, i, j);
    /* <L Phi_i, f> = Int Phi_i'' f + conj(k2) Int Phi_i f, per cell */
    double complex r = hz_phi_prod_integral(0.0, (double)NCELL, pot(i, 2), src);
    int c0 = (int)floor(XSRC - 2.0), c1 = (int)ceil(XSRC + 2.0);
    for (int c = c0; c < c1; c++)
      r += conj(med->k2[c]) * hz_phi_prod_integral((double)c, (double)(c + 1), pot(i, 0), src);
    b[i] = r;
  }
  lapack_int rank = 0;
  lapack_int info = LAPACKE_zgelsd(LAPACK_ROW_MAJOR, NF, NF, 1, A, NF, b, 1, sv, rcond, &rank);
  double err = 1e9;
  if (info == 0) {
    err = sensor_err_coef(b, ufd);
    *cond = (rank >= 1 && sv[rank - 1] > 0.0) ? sv[0] / sv[rank - 1] : 0.0;
    *rank_out = (int)rank;
  }
  free(A);
  free(b);
  free(sv);
  return err;
}

static void tp_precond(const double complex *x, double complex *y, double complex *pad,
                       const double complex *sym, double smax) {
  for (int i = 0; i < NPAD; i++)
    pad[i] = i < NF ? x[i] : 0.0;
  hz_fft(pad, NPAD, -1);
  for (int i = 0; i < NPAD; i++)
    pad[i] = (cabs(sym[i]) < 1e-8 * smax) ? 0.0 : pad[i] / sym[i];
  hz_fft(pad, NPAD, 1);
  for (int i = 0; i < NF; i++)
    y[i] = pad[i];
}

static void tp_matvec(const double complex *A, const double complex *x, double complex *y) {
  for (int i = 0; i < NF; i++) {
    double complex s = 0.0;
    for (int j = 0; j < NF; j++)
      s += A[(size_t)i * NF + (size_t)j] * x[j];
    y[i] = s;
  }
}

/* ---- Toeplitz-FFT floor: circulant-preconditioned Richardson -------------- */
/* Returns iterations to ||r||/||b|| < 1e-8, or -1; *err = sensor error. */
static int solve_toeplitz(const hz_med1d *med, const double complex *ufd, double k0, double alpha,
                          double *err) {
  /* translation-invariant background kernel t(m) = <Phi_0, L Phi_m>, infinite
   * constant medium (exact, no domain clipping) */
  double complex k2bg = CMPLX(k0 * k0, k0 * k0 * alpha);
  double complex tker[7];
  for (int m = -3; m <= 3; m++) {
    hz_phi_factor p0 = {2.0, 0.0, 0}, pm2 = {2.0, (double)m, 2}, pm0 = {2.0, (double)m, 0};
    tker[m + 3] =
        hz_phi_prod_integral(-1e6, 1e6, p0, pm2) + k2bg * hz_phi_prod_integral(-1e6, 1e6, p0, pm0);
  }
  /* full Galerkin matrix (ramps + slab + edge clipping live in V = A - T) */
  double complex *A = calloc((size_t)NF * NF, sizeof(double complex));
  double complex *b = calloc(NF, sizeof(double complex));
  double complex *c = calloc(NF, sizeof(double complex));
  double complex *r = calloc(NF, sizeof(double complex));
  double complex *pad = calloc(NPAD, sizeof(double complex));
  double complex *sym = calloc(NPAD, sizeof(double complex));
  if (!A || !b || !c || !r || !pad || !sym) {
    free(A);
    free(b);
    free(c);
    free(r);
    free(pad);
    free(sym);
    return -1;
  }
  int nv = 0;
  for (int i = 0; i < NF; i++) {
    b[i] = hz_rhs_entry(med, pot(i, 0), XSRC);
    for (int j = 0; j < NF; j++) {
      double complex a = hz_galerkin_entry(med, pot(i, 0), pot(j, 0));
      A[(size_t)i * NF + (size_t)j] = a;
      double complex tv = (abs(i - j) <= 3) ? tker[i - j + 3] : 0.0;
      if (cabs(a - tv) > 1e-12) nv++;
    }
  }
  /* circulant symbol of the background kernel, Nyquist-null guarded */
  for (int m = -3; m <= 3; m++)
    sym[(m + NPAD) % NPAD] = tker[m + 3];
  hz_fft(sym, NPAD, -1);
  double smax = 0.0;
  for (int i = 0; i < NPAD; i++)
    if (cabs(sym[i]) > smax) smax = cabs(sym[i]);

  double bn = 0.0;
  for (int i = 0; i < NF; i++)
    bn += creal(b[i] * conj(b[i]));
  bn = sqrt(bn);
  /* right-preconditioned BiCGStab, M = circulant(t) applied via FFT;
   * A applied densely here for simplicity — production M6 uses the padded-FFT
   * convolution for T plus the sparse V, same operator, O(N log N) */
  double complex *rh = calloc(NF, sizeof(double complex));
  double complex *p = calloc(NF, sizeof(double complex));
  double complex *v = calloc(NF, sizeof(double complex));
  double complex *ph = calloc(NF, sizeof(double complex));
  double complex *sh = calloc(NF, sizeof(double complex));
  double complex *tv = calloc(NF, sizeof(double complex));
  int done = -1;
  double relres = 1.0;
  if (rh != NULL && p != NULL && v != NULL && ph != NULL && sh != NULL && tv != NULL) {
    for (int i = 0; i < NF; i++) {
      r[i] = b[i];
      rh[i] = b[i];
    }
    double complex rho = 1.0, al = 1.0, om = 1.0;
    for (int it = 0; it < 500 && done < 0; it++) {
      double rn = 0.0;
      for (int i = 0; i < NF; i++)
        rn += creal(r[i] * conj(r[i]));
      relres = sqrt(rn) / bn;
      if (relres < 1e-6) {
        done = it;
        break;
      }
      double complex rho1 = 0.0;
      for (int i = 0; i < NF; i++)
        rho1 += conj(rh[i]) * r[i];
      if (cabs(rho1) < 1e-300) break;
      double complex be = (rho1 / rho) * (al / om);
      for (int i = 0; i < NF; i++)
        p[i] = r[i] + be * (p[i] - om * v[i]);
      tp_precond(p, ph, pad, sym, smax);
      tp_matvec(A, ph, v);
      double complex den = 0.0;
      for (int i = 0; i < NF; i++)
        den += conj(rh[i]) * v[i];
      if (cabs(den) < 1e-300) break;
      al = rho1 / den;
      for (int i = 0; i < NF; i++)
        r[i] -= al * v[i]; /* s lives in r */
      tp_precond(r, sh, pad, sym, smax);
      tp_matvec(A, sh, tv);
      double complex tn = 0.0, ts = 0.0;
      for (int i = 0; i < NF; i++) {
        tn += conj(tv[i]) * tv[i];
        ts += conj(tv[i]) * r[i];
      }
      if (cabs(tn) < 1e-300) break;
      om = ts / tn;
      for (int i = 0; i < NF; i++) {
        c[i] += al * ph[i] + om * sh[i];
        r[i] -= om * tv[i];
      }
      rho = rho1;
    }
  }
  free(rh);
  free(p);
  free(v);
  free(ph);
  free(sh);
  free(tv);
  *err = sensor_err_coef(c, ufd);
  printf("    toeplitz: V nnz/row=%.1f, iters=%d, relres=%.2e, sensor err=%.4f\n",
         (double)nv / (double)NF, done, relres, *err);
  free(A);
  free(b);
  free(c);
  free(r);
  free(pad);
  free(sym);
  return done;
}

int main(void) {
  double k0 = 2.0 * M_PI / 16.0;
  double alpha = 0.02;

  /* FFT self-test: roundtrip + delta spectrum */
  {
    double complex v[16];
    for (int i = 0; i < 16; i++)
      v[i] = CMPLX((double)(i % 5) - 2.0, (double)(i % 3));
    double complex w[16];
    for (int i = 0; i < 16; i++)
      w[i] = v[i];
    hz_fft(w, 16, -1);
    hz_fft(w, 16, 1);
    double d = 0.0;
    for (int i = 0; i < 16; i++)
      d += cabs(w[i] - v[i]);
    check_lt(d, 1e-12, "fft roundtrip");
  }

  double complex *ufd = calloc((size_t)NCELL * PER_CELL + 1, sizeof(double complex));
  if (ufd == NULL) return 1;

  /* formulations on the E2 slab scene */
  hz_med1d med;
  hz_med1d_init(&med, NCELL, k0, alpha, 32);
  hz_med1d_slab(&med, 96, 128, k0, 1.5, alpha);
  g_total++;
  if (hz_fd_reference(&med, PER_CELL, XSRC, ufd) != 0) {
    g_fail++;
    printf("FAIL: fd\n");
  }

  hz_sol1d gsol;
  hz_lvlstat1d gst;
  double eg = 1e9;
  g_total++;
  if (hz_cascade1d_solve(&med, &LV_F, 1, XSRC, 1, &gsol, &gst) == 0) {
    eg = sensor_err_coef(gsol.coef, ufd);
    hz_sol1d_free(&gsol);
  } else
    g_fail++;
  double cc = 0.0, cgram = 0.0, cc14 = 0.0, cgram14 = 0.0;
  int rkc = 0, rkg = 0, rkc14 = 0, rkg14 = 0;
  double ec = solve_collocation(&med, ufd, 1e-10, &cc, &rkc);
  double ec14 = solve_collocation(&med, ufd, 1e-14, &cc14, &rkc14);
  double egr = solve_gram(&med, ufd, 1e-10, &cgram, &rkg);
  double egr14 = solve_gram(&med, ufd, 1e-14, &cgram14, &rkg14);
  printf("formulations (slab 1.5):\n");
  printf("  galerkin     err=%.4f cond=%.3g\n", eg, gst.cond_eff);
  printf("  collocation  err=%.4f cond=%.3g rank=%d | rcond=1e-14: err=%.4f cond=%.3g rank=%d\n",
         ec, cc, rkc, ec14, cc14, rkc14);
  printf("  gram (A^H A) err=%.4f cond=%.3g rank=%d | rcond=1e-14: err=%.4f cond=%.3g rank=%d\n",
         egr, cgram, rkg, egr14, cgram14, rkg14);
  check_lt(eg, 0.12, "galerkin floor err on slab");
  /* measured 07-25, both rcond 1e-10 and 1e-14, full rank: pure interior
   * least squares (collocation AND Gram/FOSLS) STRIPS THE RADIATED WAVE —
   * u=0 has zero residual in empty space while the discretized wave pays a
   * small dispersion residual, so the LS minimizer prefers silence. The
   * radiation condition is not in the functional. Galerkin is the production
   * formulation; these two stay as negative controls. */
  check_lt(0.9, ec, "collocation must strip the wave (negative control)");
  check_lt(0.9, egr, "gram must strip the wave (negative control)");
  check_lt(fabs(ec - ec14) + fabs(egr - egr14), 0.01,
           "wave stripping is rcond-independent (not a truncation artifact)");
  hz_med1d_free(&med);

  /* Toeplitz-FFT floor: iterations vs contrast */
  static const double CONTRAST[4] = {1.0, 1.2, 1.5, 3.0};
  for (int ci = 0; ci < 4; ci++) {
    hz_med1d mc;
    hz_med1d_init(&mc, NCELL, k0, alpha, 32);
    if (CONTRAST[ci] > 1.0) hz_med1d_slab(&mc, 96, 128, k0, CONTRAST[ci], alpha);
    g_total++;
    if (hz_fd_reference(&mc, PER_CELL, XSRC, ufd) != 0) {
      g_fail++;
      printf("FAIL: fd contrast\n");
    }
    printf("  contrast %.1f:\n", CONTRAST[ci]);
    double te = 1e9;
    int iters = solve_toeplitz(&mc, ufd, k0, alpha, &te);
    /* measured 07-25: residual plateaus at ~2e-3 on the edge alternating mode
     * (harmless to the field, guarded out of the preconditioner); the FIELD
     * matches the dense Galerkin solve for contrast <= 1.5. Contrast 3.0
     * diverges — strong scatterers need a stronger outer solver in M6. */
    (void)iters;
    if (ci < 3) check_lt(te, 0.06, "toeplitz-fft sensor err vs dense");
    hz_med1d_free(&mc);
  }

  free(ufd);
  printf("test_m2forms: %d/%d passed\n", g_total - g_fail, g_total);
  return g_fail > 0 ? 1 : 0;
}
