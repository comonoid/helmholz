#include "helm1d.h"
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdlib.h>

/* Singular values below RCOND*s_max are treated as null modes (the alternating
 * mode of phi translates is an exact null of the infinite grid; boundary
 * truncation leaves it tiny but nonzero). */
static const double RCOND = 1e-10;

void hz_med1d_init(hz_med1d *m, int ncell, double k0, double alpha, int ramp_cells) {
  m->ncell = ncell;
  m->k2 = calloc((size_t)ncell, sizeof(double complex));
  if (m->k2 == NULL) {
    m->ncell = 0;
    return;
  }
  double b = k0 * k0;
  for (int c = 0; c < ncell; c++)
    m->k2[c] = CMPLX(b, b * alpha);
  for (int i = 0; i < ramp_cells && i < ncell; i++) {
    double t = (double)(ramp_cells - i) / (double)ramp_cells;
    double a =
        alpha + (2.0 - alpha) * t * t * t *
                    t; /* smooth quartic, deep peak: kill the wall return without ramp reflection */
    m->k2[i] = CMPLX(b, b * a);
    m->k2[ncell - 1 - i] = CMPLX(b, b * a);
  }
}

void hz_med1d_slab(hz_med1d *m, int c0, int c1, double k0, double kfac, double alpha) {
  double kk = k0 * kfac;
  for (int c = c0; c < c1 && c < m->ncell; c++)
    if (c >= 0) m->k2[c] = CMPLX(kk * kk, kk * kk * alpha);
}

void hz_med1d_free(hz_med1d *m) {
  free(m->k2);
  m->k2 = NULL;
  m->ncell = 0;
}

/* ---- FD reference: tridiagonal complex Thomas solve ---------------------- */

int hz_fd_reference(const hz_med1d *m, int per_cell, double xsrc, double complex *u) {
  if (m->ncell < 1 || per_cell < 1) return 1;
  int npts = m->ncell * per_cell + 1;
  double dx = 1.0 / (double)per_cell;
  double complex *dia = calloc((size_t)npts, sizeof(double complex));
  double complex *rhs = calloc((size_t)npts, sizeof(double complex));
  if (dia == NULL || rhs == NULL) {
    free(dia);
    free(rhs);
    return 1;
  }
  double idx2 = 1.0 / (dx * dx);
  for (int i = 0; i < npts; i++) {
    double x = (double)i * dx;
    int c = i / per_cell;
    if (c >= m->ncell) c = m->ncell - 1;
    dia[i] = -2.0 * idx2 + m->k2[c];
    rhs[i] = hz_phi(x - xsrc);
  }
  /* Dirichlet ends: u[0] = u[npts-1] = 0 */
  u[0] = 0.0;
  u[npts - 1] = 0.0;
  int n = npts - 2; /* interior unknowns u[1..npts-2] */
  /* Thomas on constant off-diagonal idx2 */
  double complex *cp = calloc((size_t)n, sizeof(double complex));
  double complex *dp = calloc((size_t)n, sizeof(double complex));
  if (cp == NULL || dp == NULL) {
    free(dia);
    free(rhs);
    free(cp);
    free(dp);
    return 1;
  }
  cp[0] = idx2 / dia[1];
  dp[0] = rhs[1] / dia[1];
  for (int i = 1; i < n; i++) {
    double complex den = dia[i + 1] - idx2 * cp[i - 1];
    cp[i] = idx2 / den;
    dp[i] = (rhs[i + 1] - idx2 * dp[i - 1]) / den;
  }
  u[n] = dp[n - 1];
  for (int i = n - 2; i >= 0; i--)
    u[i + 1] = dp[i] - cp[i] * u[i + 2];
  free(dia);
  free(rhs);
  free(cp);
  free(dp);
  return 0;
}

/* ---- Galerkin assembly over the potential basis -------------------------- */

static hz_phi_factor pot_of(const hz_level1d *L, int idx, int deriv) {
  hz_phi_factor f;
  f.h = (double)(1 << L->lvl);
  f.n = (double)(L->n0 + idx);
  f.deriv = deriv;
  return f;
}

/* <Phi_i, L Phi_j> over [0, ncell]: exact phi'' part + per-cell k2-weighted
 * mass part. Both factors deriv=0 in `a`/`b`. */
double complex hz_galerkin_entry(const hz_med1d *m, hz_phi_factor pi, hz_phi_factor pj) {
  double lo = pi.h * (pi.n - 2.0);
  double hi = pi.h * (pi.n + 2.0);
  double lo2 = pj.h * (pj.n - 2.0);
  double hi2 = pj.h * (pj.n + 2.0);
  if (lo2 > lo) lo = lo2;
  if (hi2 < hi) hi = hi2;
  if (lo < 0.0) lo = 0.0;
  if (hi > (double)m->ncell) hi = (double)m->ncell;
  if (lo >= hi) return 0.0;
  hz_phi_factor pj_dd = pj;
  pj_dd.deriv = 2;
  double complex v = hz_phi_prod_integral(lo, hi, pi, pj_dd);
  int c0 = (int)floor(lo);
  int c1 = (int)ceil(hi);
  for (int c = c0; c < c1; c++) {
    double a = (double)c, b = (double)(c + 1);
    if (a < lo) a = lo;
    if (b > hi) b = hi;
    if (a >= b) continue;
    v += m->k2[c] * hz_phi_prod_integral(a, b, pi, pj);
  }
  return v;
}

/* <Phi_i, f> with f = phi(x - xsrc) */
double complex hz_rhs_entry(const hz_med1d *m, hz_phi_factor pi, double xsrc) {
  hz_phi_factor src = {1.0, xsrc, 0};
  double lo = 0.0, hi = (double)m->ncell;
  return hz_phi_prod_integral(lo, hi, pi, src);
}

/* <Phi_i, f - L u_total> against the full current solution */
static double complex resid_entry(const hz_med1d *m, hz_phi_factor pi, double xsrc,
                                  const hz_sol1d *sol) {
  double complex r = hz_rhs_entry(m, pi, xsrc);
  for (int p = 0; p < sol->nlev; p++) {
    int np = sol->levels[p].n1 - sol->levels[p].n0 + 1;
    for (int j = 0; j < np; j++) {
      double complex c = sol->coef[sol->ofs[p] + j];
      if (cabs(c) > 0.0) r -= c * hz_galerkin_entry(m, pi, pot_of(&sol->levels[p], j, 0));
    }
  }
  return r;
}

int hz_cascade1d_solve(const hz_med1d *m, const hz_level1d *levels, int nlev, double xsrc,
                       int nsweeps, hz_sol1d *sol, hz_lvlstat1d *stats) {
  sol->levels = levels;
  sol->nlev = nlev;
  sol->ofs = calloc((size_t)nlev + 1, sizeof(int));
  if (sol->ofs == NULL) return 1;
  int total = 0;
  for (int l = 0; l < nlev; l++) {
    sol->ofs[l] = total;
    total += levels[l].n1 - levels[l].n0 + 1;
  }
  sol->ofs[nlev] = total;
  if (total < 1) {
    free(sol->ofs);
    sol->ofs = NULL;
    return 1;
  }
  sol->coef = calloc((size_t)total, sizeof(double complex));
  if (sol->coef == NULL) {
    free(sol->ofs);
    sol->ofs = NULL;
    return 1;
  }

  for (int sweep = 0; sweep < nsweeps; sweep++)
    for (int l = 0; l < nlev; l++) {
      int n = levels[l].n1 - levels[l].n0 + 1;
      double complex *A = calloc((size_t)n * (size_t)n, sizeof(double complex));
      double complex *b = calloc((size_t)n, sizeof(double complex));
      double *sv = calloc((size_t)n, sizeof(double));
      if (A == NULL || b == NULL || sv == NULL) {
        free(A);
        free(b);
        free(sv);
        return 1;
      }
      for (int i = 0; i < n; i++) {
        hz_phi_factor pi = pot_of(&levels[l], i, 0);
        for (int j = 0; j < n; j++)
          A[(size_t)i * (size_t)n + (size_t)j] = hz_galerkin_entry(m, pi, pot_of(&levels[l], j, 0));
        /* correction rhs: residual of the FULL current solution (block
         * Gauss-Seidel over levels; extra sweeps feed local refinements back
         * into the transport levels — one-way cascade cannot, see E3) */
        b[i] = resid_entry(m, pi, xsrc, sol);
      }
      double rnorm0 = 0.0;
      for (int i = 0; i < n; i++)
        rnorm0 += creal(b[i] * conj(b[i]));

      lapack_int rank = 0;
      lapack_int info = LAPACKE_zgelsd(LAPACK_ROW_MAJOR, n, n, 1, A, n, b, 1, sv, RCOND, &rank);
      if (info != 0) {
        free(A);
        free(b);
        free(sv);
        return 2;
      }
      for (int i = 0; i < n; i++)
        sol->coef[sol->ofs[l] + i] += b[i];

      if (stats != NULL) {
        stats[l].dim = n;
        stats[l].nnull = n - (int)rank;
        double smin = sv[n - 1] > 0.0 ? sv[n - 1] : RCOND * sv[0];
        stats[l].cond = sv[0] > 0.0 ? sv[0] / smin : 0.0;
        int rk = (int)rank;
        stats[l].cond_eff = (rk >= 1 && sv[rk - 1] > 0.0) ? sv[0] / sv[rk - 1] : 0.0;
        /* residual drop on this level's test space: recompute <Phi_i, f - L u> */
        double rnorm1 = 0.0;
        for (int i = 0; i < n; i++) {
          double complex r = resid_entry(m, pot_of(&levels[l], i, 0), xsrc, sol);
          rnorm1 += creal(r * conj(r));
        }
        stats[l].res_drop = rnorm0 > 0.0 ? sqrt(rnorm1 / rnorm0) : 0.0;
      }
      free(A);
      free(b);
      free(sv);
    }
  return 0;
}

double complex hz_sol1d_eval(const hz_sol1d *s, double x) {
  double complex u = 0.0;
  for (int l = 0; l < s->nlev; l++) {
    double h = (double)(1 << s->levels[l].lvl);
    int nlo = (int)floor(x / h) - 1;
    for (int n = nlo; n <= nlo + 3; n++) {
      if (n < s->levels[l].n0 || n > s->levels[l].n1) continue;
      double complex c = s->coef[s->ofs[l] + (n - s->levels[l].n0)];
      u += c * hz_phi(x / h - (double)n);
    }
  }
  return u;
}

void hz_sol1d_free(hz_sol1d *s) {
  free(s->coef);
  free(s->ofs);
  s->coef = NULL;
  s->ofs = NULL;
  s->nlev = 0;
}
