#include "solver3d.h"
#include "fft.h"
#include "phi.h"
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

enum { DENSE_CAP = 1500, BICG_MAX = 400 };
static double g_bicg_tol = 1e-6;
void hz_solver3d_set_tol(double t) {
  g_bicg_tol = t;
}
static const double RCOND = 1e-10;

static double now_ms(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec * 1e3 + (double)ts.tv_nsec * 1e-6;
}

/* ---- per-level solver state ---------------------------------------------- */

typedef struct {
  int lvl, m[3], dim;
  int P[3], pn; /* fft pad dims and volume */
  double complex *kspec;
  double complex *pad; /* scratch */
  int64_t *rowptr;     /* V csr (dim+1) or NULL */
  int32_t *col;
  double complex *val;
  double complex *dense; /* dense A (dim x dim) or NULL */
  double kscale;         /* |T(0)| for tolerances */
} lstate;

static int next_pow2(int v) {
  int p = 1;
  while (p < v)
    p *= 2;
  return p;
}

static int lidx(const int m[3], int x, int y, int z) {
  return (z * m[1] + y) * m[0] + x;
}

/* 3D FFT over padded buffer (in place), sign as in hz_fft */
/* deterministic under OpenMP: every line transform is independent, no shared
 * accumulation anywhere (the `line` parameter kept for call compatibility). */
static void fft3(double complex *a, const int P[3], int sign, double complex *line) {
  (void)line;
  size_t P0 = (size_t)P[0], P1 = (size_t)P[1];
  int lmax = P[0] > P[1] ? (P[0] > P[2] ? P[0] : P[2]) : (P[1] > P[2] ? P[1] : P[2]);
#pragma omp parallel
  {
    double complex *tl = malloc((size_t)lmax * sizeof(double complex));
    if (tl != NULL) {
#pragma omp for
      for (int zy = 0; zy < P[2] * P[1]; zy++)
        hz_fft(a + (size_t)zy * P0, P[0], sign);
#pragma omp for
      for (int zx = 0; zx < P[2] * P[0]; zx++) {
        int z = zx / P[0], x = zx % P[0];
        for (int y = 0; y < P[1]; y++)
          tl[y] = a[((size_t)z * P1 + (size_t)y) * P0 + (size_t)x];
        hz_fft(tl, P[1], sign);
        for (int y = 0; y < P[1]; y++)
          a[((size_t)z * P1 + (size_t)y) * P0 + (size_t)x] = tl[y];
      }
#pragma omp for
      for (int yx = 0; yx < P[1] * P[0]; yx++) {
        int y = yx / P[0], x = yx % P[0];
        for (int z = 0; z < P[2]; z++)
          tl[z] = a[((size_t)z * P1 + (size_t)y) * P0 + (size_t)x];
        hz_fft(tl, P[2], sign);
        for (int z = 0; z < P[2]; z++)
          a[((size_t)z * P1 + (size_t)y) * P0 + (size_t)x] = tl[z];
      }
    }
    free(tl);
  }
}

/* infinite-background Toeplitz kernel t(mx,my,mz), support |m|<=3 */
static void kernel_build(lstate *L, double complex k2bg) {
  int h = 1 << L->lvl;
  double D[7], M[7];
  for (int m = -3; m <= 3; m++) {
    hz_phi_factor f0 = {(double)h, 0.0, 0};
    hz_phi_factor fm0 = {(double)h, (double)m, 0};
    hz_phi_factor fm2 = {(double)h, (double)m, 2};
    double lo = -4.0 * (double)h, hi = 4.0 * (double)h;
    D[m + 3] = hz_phi_prod_integral(lo, hi, f0, fm2);
    M[m + 3] = hz_phi_prod_integral(lo, hi, f0, fm0);
  }
  memset(L->kspec, 0, (size_t)L->pn * sizeof(double complex));
  for (int mz = -3; mz <= 3; mz++)
    for (int my = -3; my <= 3; my++)
      for (int mx = -3; mx <= 3; mx++) {
        double complex t = D[mx + 3] * M[my + 3] * M[mz + 3] + M[mx + 3] * D[my + 3] * M[mz + 3] +
                           M[mx + 3] * M[my + 3] * D[mz + 3] +
                           k2bg * M[mx + 3] * M[my + 3] * M[mz + 3];
        int ix = (mx + L->P[0]) % L->P[0];
        int iy = (my + L->P[1]) % L->P[1];
        int iz = (mz + L->P[2]) % L->P[2];
        L->kspec[((size_t)iz * (size_t)L->P[1] + (size_t)iy) * (size_t)L->P[0] + (size_t)ix] = t;
      }
  L->kscale = cabs(L->kspec[0]);
  double complex *line =
      malloc((size_t)next_pow2(L->P[0] > L->P[1] ? (L->P[0] > L->P[2] ? L->P[0] : L->P[2])
                                                 : (L->P[1] > L->P[2] ? L->P[1] : L->P[2])) *
             sizeof(double complex));
  if (line != NULL) {
    fft3(L->kspec, L->P, -1, line);
    free(line);
  }
}

/* y = T x via padded FFT (exact linear convolution: pad >= m+3 per axis) */
static void matvec_T(lstate *L, const double complex *x, double complex *y, double complex *line) {
  memset(L->pad, 0, (size_t)L->pn * sizeof(double complex));
  for (int z = 0; z < L->m[2]; z++)
    for (int yy = 0; yy < L->m[1]; yy++)
      memcpy(L->pad + ((size_t)z * (size_t)L->P[1] + (size_t)yy) * (size_t)L->P[0],
             x + (size_t)lidx(L->m, 0, yy, z), (size_t)L->m[0] * sizeof(double complex));
  fft3(L->pad, L->P, -1, line);
  for (int i = 0; i < L->pn; i++)
    L->pad[i] *= L->kspec[i];
  fft3(L->pad, L->P, 1, line);
  for (int z = 0; z < L->m[2]; z++)
    for (int yy = 0; yy < L->m[1]; yy++)
      memcpy(y + (size_t)lidx(L->m, 0, yy, z),
             L->pad + ((size_t)z * (size_t)L->P[1] + (size_t)yy) * (size_t)L->P[0],
             (size_t)L->m[0] * sizeof(double complex));
}

static void matvec_A(lstate *L, const double complex *x, double complex *y, double complex *line) {
  matvec_T(L, x, y, line);
  if (L->rowptr != NULL) {
/* rows are independent; per-row sums keep a fixed order — deterministic */
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L->dim; i++) {
      double complex s = 0.0;
      for (int64_t k = L->rowptr[i]; k < L->rowptr[i + 1]; k++)
        s += L->val[k] * x[L->col[k]];
      y[i] += s;
    }
  }
}

/* z = P^-1 r: circulant division with a guard on near-null symbol bins */
static void precond(lstate *L, const double complex *r, double complex *z, double complex *line) {
  memset(L->pad, 0, (size_t)L->pn * sizeof(double complex));
  for (int zz = 0; zz < L->m[2]; zz++)
    for (int yy = 0; yy < L->m[1]; yy++)
      memcpy(L->pad + ((size_t)zz * (size_t)L->P[1] + (size_t)yy) * (size_t)L->P[0],
             r + (size_t)lidx(L->m, 0, yy, zz), (size_t)L->m[0] * sizeof(double complex));
  fft3(L->pad, L->P, -1, line);
  for (int i = 0; i < L->pn; i++)
    L->pad[i] = (cabs(L->kspec[i]) < 1e-8 * L->kscale) ? 0.0 : L->pad[i] / L->kspec[i];
  fft3(L->pad, L->P, 1, line);
  for (int zz = 0; zz < L->m[2]; zz++)
    for (int yy = 0; yy < L->m[1]; yy++)
      memcpy(z + (size_t)lidx(L->m, 0, yy, zz),
             L->pad + ((size_t)zz * (size_t)L->P[1] + (size_t)yy) * (size_t)L->P[0],
             (size_t)L->m[0] * sizeof(double complex));
}

/* ---- V: deviation of the exact Galerkin matrix from the Toeplitz kernel --- */

typedef struct {
  char *mark;
  const int *m;
  int lvl, dom0, dom1, dom2;
} markctx;

static int mark_leaf_cb(void *vctx, const int blo[3], const int bhi[3], double complex k2) {
  (void)k2;
  markctx *mc = vctx;
  int h = 1 << mc->lvl;
  int lo[3], hi[3];
  for (int a = 0; a < 3; a++) {
    lo[a] = blo[a] / h - 2;
    hi[a] = bhi[a] / h + 2;
    lo[a] += 1; /* index offset: n0 = -1 */
    hi[a] += 1;
    if (lo[a] < 0) lo[a] = 0;
    if (hi[a] > mc->m[a] - 1) hi[a] = mc->m[a] - 1;
  }
  for (int z = lo[2]; z <= hi[2]; z++)
    for (int y = lo[1]; y <= hi[1]; y++)
      for (int x = lo[0]; x <= hi[0]; x++)
        mc->mark[lidx(mc->m, x, y, z)] = 1;
  return 0;
}

typedef struct {
  markctx mc;
  const hz_octree *tree;
  double complex k2bg;
} devmark;

static int dev_mark_cb(void *vctx, const int blo[3], const int bhi[3], double complex k2) {
  devmark *dm = vctx;
  if (cabs(k2 - dm->k2bg) > 1e-14 * (cabs(dm->k2bg) + 1.0)) mark_leaf_cb(&dm->mc, blo, bhi, k2);
  return 0;
}

/* kernel value for offset d (or 0 outside the band) */
static double complex kernel_at(const double complex *tk, const int d[3]) {
  for (int a = 0; a < 3; a++)
    if (d[a] < -3 || d[a] > 3) return 0.0;
  return tk[((d[2] + 3) * 7 + (d[1] + 3)) * 7 + (d[0] + 3)];
}

static int build_V(lstate *L, const hz_scene3 *sc) {
  int h = 1 << L->lvl;
  /* raw band kernel (not the spectrum) for subtraction */
  double complex tk[343];
  {
    double D[7], M[7];
    for (int m = -3; m <= 3; m++) {
      hz_phi_factor f0 = {(double)h, 0.0, 0};
      hz_phi_factor fm0 = {(double)h, (double)m, 0};
      hz_phi_factor fm2 = {(double)h, (double)m, 2};
      double lo = -4.0 * (double)h, hi = 4.0 * (double)h;
      D[m + 3] = hz_phi_prod_integral(lo, hi, f0, fm2);
      M[m + 3] = hz_phi_prod_integral(lo, hi, f0, fm0);
    }
    for (int mz = -3; mz <= 3; mz++)
      for (int my = -3; my <= 3; my++)
        for (int mx = -3; mx <= 3; mx++)
          tk[((mz + 3) * 7 + (my + 3)) * 7 + (mx + 3)] =
              D[mx + 3] * M[my + 3] * M[mz + 3] + M[mx + 3] * D[my + 3] * M[mz + 3] +
              M[mx + 3] * M[my + 3] * D[mz + 3] + sc->k2bg * M[mx + 3] * M[my + 3] * M[mz + 3];
  }

  char *mark = calloc((size_t)L->dim, 1);
  if (mark == NULL) return 1;
  /* domain-edge clipping rows */
  for (int z = 0; z < L->m[2]; z++)
    for (int y = 0; y < L->m[1]; y++)
      for (int x = 0; x < L->m[0]; x++) {
        int n[3] = {x - 1, y - 1, z - 1};
        int clip = 0;
        for (int a = 0; a < 3; a++)
          if (h * (n[a] - 2) < 0 || h * (n[a] + 2) > sc->dom[a]) clip = 1;
        if (clip) mark[lidx(L->m, x, y, z)] = 1;
      }
  /* medium deviation rows */
  devmark dm = {{mark, L->m, L->lvl, sc->dom[0], sc->dom[1], sc->dom[2]}, sc->tree, sc->k2bg};
  int zlo[3] = {0, 0, 0};
  int zhi[3] = {sc->dom[0], sc->dom[1], sc->dom[2]};
  hz_oct_visit(sc->tree, zlo, zhi, dev_mark_cb, &dm);

  /* count and fill CSR: two passes over independent rows, parallel and
   * deterministic (pass 1 writes into disjoint precomputed row ranges) */
  double tol = 1e-13 * L->kscale;
  L->rowptr = calloc((size_t)L->dim + 1, sizeof(int64_t));
  if (L->rowptr == NULL) {
    free(mark);
    return 1;
  }
#pragma omp parallel for schedule(dynamic, 8)
  for (int i = 0; i < L->dim; i++) {
    if (!mark[i]) continue;
    int x = i % L->m[0], y = (i / L->m[0]) % L->m[1], z = i / (L->m[0] * L->m[1]);
    hz_pot3 a = {L->lvl, {x - 1, y - 1, z - 1}};
    int64_t cnt = 0;
    for (int dz = -3; dz <= 3; dz++)
      for (int dy = -3; dy <= 3; dy++)
        for (int dx = -3; dx <= 3; dx++) {
          int jx = x + dx, jy = y + dy, jz = z + dz;
          if (jx < 0 || jy < 0 || jz < 0 || jx >= L->m[0] || jy >= L->m[1] || jz >= L->m[2])
            continue;
          hz_pot3 b = {L->lvl, {jx - 1, jy - 1, jz - 1}};
          int d[3] = {dx, dy, dz};
          if (cabs(hz_entry3d(sc->tree, sc->dom, a, b) - kernel_at(tk, d)) > tol) cnt++;
        }
    L->rowptr[i + 1] = cnt;
  }
  for (int i = 0; i < L->dim; i++)
    L->rowptr[i + 1] += L->rowptr[i];
  int64_t nnz = L->rowptr[L->dim];
  L->col = malloc((size_t)nnz * sizeof(int32_t));
  L->val = malloc((size_t)nnz * sizeof(double complex));
  if (L->col == NULL || L->val == NULL) {
    free(mark);
    return 1;
  }
#pragma omp parallel for schedule(dynamic, 8)
  for (int i = 0; i < L->dim; i++) {
    if (!mark[i]) continue;
    int x = i % L->m[0], y = (i / L->m[0]) % L->m[1], z = i / (L->m[0] * L->m[1]);
    hz_pot3 a = {L->lvl, {x - 1, y - 1, z - 1}};
    int64_t w = L->rowptr[i];
    for (int dz = -3; dz <= 3; dz++)
      for (int dy = -3; dy <= 3; dy++)
        for (int dx = -3; dx <= 3; dx++) {
          int jx = x + dx, jy = y + dy, jz = z + dz;
          if (jx < 0 || jy < 0 || jz < 0 || jx >= L->m[0] || jy >= L->m[1] || jz >= L->m[2])
            continue;
          hz_pot3 b = {L->lvl, {jx - 1, jy - 1, jz - 1}};
          int d[3] = {dx, dy, dz};
          double complex v = hz_entry3d(sc->tree, sc->dom, a, b) - kernel_at(tk, d);
          if (cabs(v) <= tol) continue;
          L->col[w] = (int32_t)lidx(L->m, jx, jy, jz);
          L->val[w] = v;
          w++;
        }
  }
  free(mark);
  return 0;
}

/* ---- BiCGStab -------------------------------------------------------------- */

static int bicgstab(lstate *L, const double complex *b, double complex *x, double *relres_out) {
  int n = L->dim;
  double complex *r = calloc((size_t)n, sizeof(double complex));
  double complex *rh = calloc((size_t)n, sizeof(double complex));
  double complex *p = calloc((size_t)n, sizeof(double complex));
  double complex *v = calloc((size_t)n, sizeof(double complex));
  double complex *ph = calloc((size_t)n, sizeof(double complex));
  double complex *sh = calloc((size_t)n, sizeof(double complex));
  double complex *tv = calloc((size_t)n, sizeof(double complex));
  double complex *line =
      calloc((size_t)next_pow2(L->P[0] | L->P[1] | L->P[2]), sizeof(double complex));
  int done = -1;
  double relres = 1.0;
  if (r && rh && p && v && ph && sh && tv && line) {
    double bn = 0.0;
    for (int i = 0; i < n; i++) {
      r[i] = b[i];
      rh[i] = b[i];
      x[i] = 0.0;
      bn += creal(b[i] * conj(b[i]));
    }
    bn = sqrt(bn);
    if (bn > 0.0) {
      double complex rho = 1.0, al = 1.0, om = 1.0;
      double best = 1e300;
      int best_it = 0;
      for (int it = 0; it < BICG_MAX && done < 0; it++) {
        double rn = 0.0;
        for (int i = 0; i < n; i++)
          rn += creal(r[i] * conj(r[i]));
        relres = sqrt(rn) / bn;
        if (relres < g_bicg_tol) {
          done = it;
          break;
        }
        if (relres < 0.7 * best) {
          best = relres;
          best_it = it;
        } else if (it - best_it > 40) {
          /* stagnation plateau (edge near-null modes) — the field is done
           * even though the raw residual is not (measured in 1D, M2) */
          done = it;
          break;
        }
        double complex rho1 = 0.0;
        for (int i = 0; i < n; i++)
          rho1 += conj(rh[i]) * r[i];
        if (cabs(rho1) < 1e-290) break;
        double complex be = (rho1 / rho) * (al / om);
        for (int i = 0; i < n; i++)
          p[i] = r[i] + be * (p[i] - om * v[i]);
        precond(L, p, ph, line);
        matvec_A(L, ph, v, line);
        double complex den = 0.0;
        for (int i = 0; i < n; i++)
          den += conj(rh[i]) * v[i];
        if (cabs(den) < 1e-290) break;
        al = rho1 / den;
        for (int i = 0; i < n; i++)
          r[i] -= al * v[i];
        precond(L, r, sh, line);
        matvec_A(L, sh, tv, line);
        double complex tn = 0.0, ts = 0.0;
        for (int i = 0; i < n; i++) {
          tn += conj(tv[i]) * tv[i];
          ts += conj(tv[i]) * r[i];
        }
        if (cabs(tn) < 1e-290) break;
        om = ts / tn;
        for (int i = 0; i < n; i++) {
          x[i] += al * ph[i] + om * sh[i];
          r[i] -= om * tv[i];
        }
        rho = rho1;
      }
    } else
      done = 0;
  }
  free(r);
  free(rh);
  free(p);
  free(v);
  free(ph);
  free(sh);
  free(tv);
  free(line);
  *relres_out = relres;
  return done;
}

/* ---- two-scale transfer ---------------------------------------------------- */

static const double W5[5] = {0.25, 0.5, 0.5, 0.5, 0.25}; /* [1,2,2,2,1]/4 */

/* add prolongation of coarse (mc, offsets -1) into fine (mf): fine_abs = 2*coarse_abs + s */
static void prolong_add(const int mc[3], const double complex *cc, const int mf[3],
                        double complex *cf) {
  for (int z = 0; z < mc[2]; z++)
    for (int y = 0; y < mc[1]; y++)
      for (int x = 0; x < mc[0]; x++) {
        double complex v = cc[lidx(mc, x, y, z)];
        if (!(cabs(v) > 0.0)) continue;
        int ax = x - 1, ay = y - 1, az = z - 1; /* absolute coarse offsets */
        for (int sz = -2; sz <= 2; sz++)
          for (int sy = -2; sy <= 2; sy++)
            for (int sx = -2; sx <= 2; sx++) {
              int fx = 2 * ax + sx + 1, fy = 2 * ay + sy + 1, fz = 2 * az + sz + 1;
              if (fx < 0 || fy < 0 || fz < 0 || fx >= mf[0] || fy >= mf[1] || fz >= mf[2]) continue;
              cf[lidx(mf, fx, fy, fz)] += v * W5[sx + 2] * W5[sy + 2] * W5[sz + 2];
            }
      }
}

/* coarse = restriction of fine (same weights, transposed) */
static void restrict_to(const int mf[3], const double complex *rf, const int mc[3],
                        double complex *rc) {
  for (int z = 0; z < mc[2]; z++)
    for (int y = 0; y < mc[1]; y++)
      for (int x = 0; x < mc[0]; x++) {
        int ax = x - 1, ay = y - 1, az = z - 1;
        double complex s = 0.0;
        for (int sz = -2; sz <= 2; sz++)
          for (int sy = -2; sy <= 2; sy++)
            for (int sx = -2; sx <= 2; sx++) {
              int fx = 2 * ax + sx + 1, fy = 2 * ay + sy + 1, fz = 2 * az + sz + 1;
              if (fx < 0 || fy < 0 || fz < 0 || fx >= mf[0] || fy >= mf[1] || fz >= mf[2]) continue;
              s += rf[lidx(mf, fx, fy, fz)] * W5[sx + 2] * W5[sy + 2] * W5[sz + 2];
            }
        rc[lidx(mc, x, y, z)] = s;
      }
}

/* ---- driver ----------------------------------------------------------------- */

int hz_solve3d(const hz_scene3 *sc, int top_lvl, int floor_lvl, int nsweeps, int verbose,
               hz_sol3 *out) {
  int nlev = top_lvl - floor_lvl + 1;
  if (nlev < 1) return 1;
  for (int a = 0; a < 3; a++)
    if (sc->dom[a] % (1 << top_lvl) != 0) return 1;

  out->nlev = nlev;
  out->dom[0] = sc->dom[0];
  out->dom[1] = sc->dom[1];
  out->dom[2] = sc->dom[2];
  out->lv = calloc((size_t)nlev, sizeof(hz_slevel));
  lstate *st = calloc((size_t)nlev, sizeof(lstate));
  if (out->lv == NULL || st == NULL) {
    free(out->lv);
    free(st);
    return 1;
  }

  for (int l = 0; l < nlev; l++) {
    int lvl = top_lvl - l;
    int h = 1 << lvl;
    out->lv[l].lvl = lvl;
    st[l].lvl = lvl;
    st[l].dim = 1;
    for (int a = 0; a < 3; a++) {
      out->lv[l].m[a] = sc->dom[a] / h + 3;
      st[l].m[a] = out->lv[l].m[a];
      st[l].dim *= st[l].m[a];
    }
    out->lv[l].coef = calloc((size_t)st[l].dim, sizeof(double complex));
    if (out->lv[l].coef == NULL) return 1;
  }
  int FL = nlev - 1; /* floor index */

  double t0 = now_ms();
  /* per-level operators */
  for (int l = 0; l < nlev; l++) {
    lstate *L = &st[l];
    if (L->dim <= DENSE_CAP && l != FL) {
      L->dense = malloc((size_t)L->dim * (size_t)L->dim * sizeof(double complex));
      if (L->dense == NULL) return 1;
#pragma omp parallel for schedule(dynamic, 4)
      for (int i = 0; i < L->dim; i++) {
        int x = i % L->m[0], y = (i / L->m[0]) % L->m[1], z = i / (L->m[0] * L->m[1]);
        hz_pot3 a = {L->lvl, {x - 1, y - 1, z - 1}};
        for (int j = 0; j < L->dim; j++) {
          int xx = j % L->m[0], yy = (j / L->m[0]) % L->m[1], zz = j / (L->m[0] * L->m[1]);
          hz_pot3 b = {L->lvl, {xx - 1, yy - 1, zz - 1}};
          L->dense[(size_t)i * (size_t)L->dim + (size_t)j] = hz_entry3d(sc->tree, sc->dom, a, b);
        }
      }
    } else {
      for (int a = 0; a < 3; a++)
        L->P[a] = next_pow2(L->m[a] + 6);
      L->pn = L->P[0] * L->P[1] * L->P[2];
      L->kspec = malloc((size_t)L->pn * sizeof(double complex));
      L->pad = malloc((size_t)L->pn * sizeof(double complex));
      if (L->kspec == NULL || L->pad == NULL) return 1;
      kernel_build(L, sc->k2bg);
      if (build_V(L, sc)) return 1;
    }
    if (verbose)
      printf("  level h=%-3d dim=%-8d %s nnzV=%lld\n", 1 << L->lvl, L->dim,
             L->dense ? "dense" : "fft+V", L->rowptr ? (long long)L->rowptr[L->dim] : 0LL);
  }
  double t1 = now_ms();

  /* floor rhs */
  lstate *F = &st[FL];
  double complex *bf = calloc((size_t)F->dim, sizeof(double complex));
  double complex *img = calloc((size_t)F->dim, sizeof(double complex));
  double complex *rf = calloc((size_t)F->dim, sizeof(double complex));
  double complex *line =
      calloc((size_t)next_pow2(F->P[0] | F->P[1] | F->P[2]), sizeof(double complex));
  double complex **tmp = calloc((size_t)nlev, sizeof(double complex *));
  if (bf == NULL || img == NULL || rf == NULL || line == NULL || tmp == NULL) return 1;
  for (int l = 0; l < nlev; l++) {
    tmp[l] = calloc((size_t)st[l].dim, sizeof(double complex));
    if (tmp[l] == NULL) return 1;
  }
#pragma omp parallel for schedule(static)
  for (int i = 0; i < F->dim; i++) {
    int x = i % F->m[0], y = (i / F->m[0]) % F->m[1], z = i / (F->m[0] * F->m[1]);
    hz_pot3 a = {F->lvl, {x - 1, y - 1, z - 1}};
    double complex s = 0.0;
    for (int q = 0; q < sc->nsrc; q++)
      s += hz_rhs3d(sc->dom, a, sc->src[q].p, sc->src[q].amp);
    bf[i] = s;
  }

  double relres = 1.0;
  for (int sweep = 0; sweep < nsweeps; sweep++) {
    for (int l = 0; l < nlev; l++) {
      /* floor image of the whole current solution */
      memcpy(img, out->lv[FL].coef, (size_t)F->dim * sizeof(double complex));
      for (int c = 0; c < FL; c++) {
        memcpy(tmp[c], out->lv[c].coef, (size_t)st[c].dim * sizeof(double complex));
        for (int step = c; step < FL; step++) {
          memset(tmp[step + 1], 0, (size_t)st[step + 1].dim * sizeof(double complex));
          prolong_add(st[step].m, tmp[step], st[step + 1].m, tmp[step + 1]);
        }
        for (int i = 0; i < F->dim; i++)
          img[i] += tmp[FL][i];
      }
      /* floor residual */
      matvec_A(F, img, rf, line);
      double rn = 0.0, bn = 0.0;
      for (int i = 0; i < F->dim; i++) {
        rf[i] = bf[i] - rf[i];
        rn += creal(rf[i] * conj(rf[i]));
        bn += creal(bf[i] * conj(bf[i]));
      }
      relres = bn > 0.0 ? sqrt(rn / bn) : 0.0;

      /* restrict to level l */
      memcpy(tmp[FL], rf, (size_t)F->dim * sizeof(double complex));
      for (int step = FL; step > l; step--)
        restrict_to(st[step].m, tmp[step], st[step - 1].m, tmp[step - 1]);

      /* solve the level for a correction */
      lstate *L = &st[l];
      if (L->dense != NULL) {
        double complex *A = malloc((size_t)L->dim * (size_t)L->dim * sizeof(double complex));
        double *sv = malloc((size_t)L->dim * sizeof(double));
        if (A == NULL || sv == NULL) {
          free(A);
          free(sv);
          return 1;
        }
        memcpy(A, L->dense, (size_t)L->dim * (size_t)L->dim * sizeof(double complex));
        lapack_int rank = 0;
        lapack_int info = LAPACKE_zgelsd(LAPACK_ROW_MAJOR, L->dim, L->dim, 1, A, L->dim, tmp[l], 1,
                                         sv, RCOND, &rank);
        free(A);
        free(sv);
        if (info != 0) return 2;
        for (int i = 0; i < L->dim; i++)
          out->lv[l].coef[i] += tmp[l][i];
      } else {
        double complex *dx = calloc((size_t)L->dim, sizeof(double complex));
        if (dx == NULL) return 1;
        double rr = 0.0;
        int it = bicgstab(L, tmp[l], dx, &rr);
        if (verbose)
          printf("  sweep %d level h=%-3d bicg iters=%d relres=%.2e floor_relres=%.2e\n", sweep,
                 1 << L->lvl, it, rr, relres);
        for (int i = 0; i < L->dim; i++)
          out->lv[l].coef[i] += dx[i];
        free(dx);
      }
    }
  }

  /* final floor residual for the report */
  memcpy(img, out->lv[FL].coef, (size_t)F->dim * sizeof(double complex));
  for (int c = 0; c < FL; c++) {
    memcpy(tmp[c], out->lv[c].coef, (size_t)st[c].dim * sizeof(double complex));
    for (int step = c; step < FL; step++) {
      memset(tmp[step + 1], 0, (size_t)st[step + 1].dim * sizeof(double complex));
      prolong_add(st[step].m, tmp[step], st[step + 1].m, tmp[step + 1]);
    }
    for (int i = 0; i < F->dim; i++)
      img[i] += tmp[FL][i];
  }
  matvec_A(F, img, rf, line);
  double rn = 0.0, bn = 0.0;
  for (int i = 0; i < F->dim; i++) {
    rn += creal((bf[i] - rf[i]) * conj(bf[i] - rf[i]));
    bn += creal(bf[i] * conj(bf[i]));
  }
  out->final_relres = bn > 0.0 ? sqrt(rn / bn) : 0.0;
  if (verbose)
    printf("  setup %.0f ms, solve %.0f ms, final floor relres %.2e\n", t1 - t0, now_ms() - t1,
           out->final_relres);

  for (int l = 0; l < nlev; l++) {
    free(tmp[l]);
    free(st[l].kspec);
    free(st[l].pad);
    free(st[l].rowptr);
    free(st[l].col);
    free(st[l].val);
    free(st[l].dense);
  }
  free(tmp);
  free(bf);
  free(img);
  free(rf);
  free(line);
  free(st);
  return 0;
}

double complex hz_sol3_eval(const hz_sol3 *s, const double p[3]) {
  double complex u = 0.0;
  for (int l = 0; l < s->nlev; l++) {
    const hz_slevel *L = &s->lv[l];
    double h = (double)(1 << L->lvl);
    int lo[3];
    for (int a = 0; a < 3; a++)
      lo[a] = (int)floor(p[a] / h) - 1;
    for (int dz = 0; dz < 4; dz++)
      for (int dy = 0; dy < 4; dy++)
        for (int dx = 0; dx < 4; dx++) {
          int n[3] = {lo[0] + dx, lo[1] + dy, lo[2] + dz};
          int ok = 1;
          for (int a = 0; a < 3; a++)
            if (n[a] < -1 || n[a] > lo[a] + 3 || n[a] + 1 >= L->m[a] || n[a] + 1 < 0) ok = 0;
          if (!ok) continue;
          double v = 1.0;
          for (int a = 0; a < 3; a++)
            v *= hz_phi(p[a] / h - (double)n[a]);
          if (fabs(v) > 0.0)
            u += L->coef[((n[2] + 1) * L->m[1] + (n[1] + 1)) * L->m[0] + (n[0] + 1)] * v;
        }
  }
  return u;
}

void hz_sol3_free(hz_sol3 *s) {
  for (int l = 0; l < s->nlev; l++)
    free(s->lv[l].coef);
  free(s->lv);
  s->lv = NULL;
  s->nlev = 0;
}

int hz_dbg_matvec(const hz_scene3 *sc, int lvl, const double complex *x, double complex *y) {
  lstate L;
  memset(&L, 0, sizeof(L));
  L.lvl = lvl;
  L.dim = 1;
  for (int a = 0; a < 3; a++) {
    L.m[a] = sc->dom[a] / (1 << lvl) + 3;
    L.P[a] = next_pow2(L.m[a] + 6);
    L.dim *= L.m[a];
  }
  L.pn = L.P[0] * L.P[1] * L.P[2];
  L.kspec = malloc((size_t)L.pn * sizeof(double complex));
  L.pad = malloc((size_t)L.pn * sizeof(double complex));
  double complex *line =
      malloc((size_t)next_pow2(L.P[0] | L.P[1] | L.P[2]) * sizeof(double complex));
  if (L.kspec == NULL || L.pad == NULL || line == NULL) {
    free(L.kspec);
    free(L.pad);
    free(line);
    return 1;
  }
  kernel_build(&L, sc->k2bg);
  if (build_V(&L, sc)) return 1;
  matvec_A(&L, x, y, line);
  free(L.kspec);
  free(L.pad);
  free(L.rowptr);
  free(L.col);
  free(L.val);
  free(line);
  return 0;
}
