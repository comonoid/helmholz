/* rte2d.c — T2: рассеяние, DG1, ограничитель. См. rte2d.h. */
#include "rte2d.h"
#include "quad.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Итерация по рассеянию сходится как ρⁿ по альбедо (аудит С1). TR2_TOL —
 * относительное изменение φ; 1e-15 это несколько ulp, ниже итерация упирается
 * в округление и не останавливается. TR2_MAXIT: при альбедо 1 и τ~2 утечка
 * даёт спектральный радиус ~0.7 ⇒ ~100 итераций; при 0.99 счёт идёт на сотни
 * (С1), поэтому 20000 — предохранитель, а не рабочая величина. */
#define TR2_TOL 1e-15
#define TR2_MAXIT 20000

/* ============================================================ направления */

int tr2_dirs_quadrant_gauss(tr2_dirs *d, int nper) {
  d->n = 0;
  d->ox = d->oy = d->w = NULL;
  if (nper < 1) return -1;

  double *gx = calloc((size_t)nper, sizeof *gx);
  double *gw = calloc((size_t)nper, sizeof *gw);
  if (!gx || !gw) {
    free(gx);
    free(gw);
    return -1;
  }
  tr_gauss_legendre(nper, gx, gw);

  int n = 4 * nper;
  double *ox = calloc((size_t)n, sizeof *ox);
  double *oy = calloc((size_t)n, sizeof *oy);
  double *w = calloc((size_t)n, sizeof *w);
  if (!ox || !oy || !w) {
    free(gx);
    free(gw);
    free(ox);
    free(oy);
    free(w);
    return -1;
  }
  /* квадрант q занимает азимут [q·π/2, (q+1)·π/2]; полуширина π/4, вес wg·π/4
   * ⇒ Σ по квадранту = π/2, по всем четырём = 2π */
  const double halfw = M_PI / 4.0;
  for (int q = 0; q < 4; q++) {
    double mid = ((double)q + 0.5) * (M_PI / 2.0);
    for (int i = 0; i < nper; i++) {
      double th = mid + halfw * gx[i];
      int k = q * nper + i;
      ox[k] = cos(th);
      oy[k] = sin(th);
      w[k] = halfw * gw[i];
    }
  }
  free(gx);
  free(gw);
  d->n = n;
  d->ox = ox;
  d->oy = oy;
  d->w = w;
  return 0;
}

void tr2_dirs_free(tr2_dirs *d) {
  free(d->ox);
  free(d->oy);
  free(d->w);
  d->ox = d->oy = d->w = NULL;
  d->n = 0;
}

void tr2_phase_isotropic(tr2_phase *ph, const tr2_dirs *d) {
  double s = 0.0;
  for (int m = 0; m < d->n; m++)
    s += d->w[m];
  ph->p0 = 1.0 / s; /* НЕ 1/(2π): так рассеяние сохраняет энергию точно */
}

/* ================================================================= сетка */

int tr2_mesh_init(tr2_mesh *m, int nx, int ny, double ux, double uy, int ccap) {
  memset(m, 0, sizeof *m);
  if (nx < 1 || ny < 1 || ccap < 1) return -1;
  tr2_cell *cell = calloc((size_t)ccap, sizeof *cell);
  if (!cell) return -1;
  m->nx = nx;
  m->ny = ny;
  m->ux = ux;
  m->uy = uy;
  m->ccap = ccap;
  m->cell = cell;
  return 0;
}

int tr2_mesh_add(tr2_mesh *m, int ix0, int iy0, int ix1, int iy1) {
  if (m->ncell >= m->ccap) return -1;
  if (ix0 < 0 || iy0 < 0 || ix1 > m->nx || iy1 > m->ny || ix1 <= ix0 || iy1 <= iy0) return -1;
  m->cell[m->ncell].ix0 = ix0;
  m->cell[m->ncell].iy0 = iy0;
  m->cell[m->ncell].ix1 = ix1;
  m->cell[m->ncell].iy1 = iy1;
  m->ncell++;
  return 0;
}

void tr2_mesh_free(tr2_mesh *m) {
  free(m->cell);
  free(m->face);
  free(m->coff);
  free(m->cfac);
  memset(m, 0, sizeof *m);
}

/* полуграни: ячейка + плоскость, к которой она примыкает.
 * side 0 — ячейка с − стороны плоскости (плоскость = её верхняя грань),
 * side 1 — ячейка с + стороны (плоскость = её нижняя грань). */
typedef struct {
  int axis, ipos, it0, it1, cell, side;
} halfface;

static int hf_cmp(const void *pa, const void *pb) {
  const halfface *a = pa, *b = pb;
  if (a->axis != b->axis) return a->axis < b->axis ? -1 : 1;
  if (a->ipos != b->ipos) return a->ipos < b->ipos ? -1 : 1;
  if (a->side != b->side) return a->side < b->side ? -1 : 1;
  if (a->it0 != b->it0) return a->it0 < b->it0 ? -1 : 1;
  return 0;
}

/* Кто из полуграней группы [g0,g1) с нужной стороной покрывает [p,q). Отрезки
 * одной стороны в группе не пересекаются, поэтому такой не более одного. */
static int covering(const halfface *hf, int g0, int g1, int side, int p, int q) {
  for (int i = g0; i < g1; i++)
    if (hf[i].side == side && hf[i].it0 <= p && hf[i].it1 >= q) return hf[i].cell;
  return -1;
}

/* Один проход построения граней: считает (out == NULL) либо заполняет. */
static int build_faces(const halfface *hf, int nhf, tr2_face *out) {
  int nf = 0;
  int g0 = 0;
  while (g0 < nhf) {
    int g1 = g0;
    while (g1 < nhf && hf[g1].axis == hf[g0].axis && hf[g1].ipos == hf[g0].ipos)
      g1++;
    /* элементарные отрезки = все концы группы, отсортированные и уникальные */
    int nb = 0;
    int bmax = 2 * (g1 - g0);
    int *b = calloc((size_t)bmax, sizeof *b);
    if (!b) return -1;
    for (int i = g0; i < g1; i++) {
      b[nb++] = hf[i].it0;
      b[nb++] = hf[i].it1;
    }
    /* вставками: групп мало и они короткие */
    for (int i = 1; i < nb; i++) {
      int v = b[i], j = i - 1;
      while (j >= 0 && b[j] > v) {
        b[j + 1] = b[j];
        j--;
      }
      b[j + 1] = v;
    }
    for (int i = 0; i < nb - 1; i++) {
      if (b[i] == b[i + 1]) continue;
      int lo = covering(hf, g0, g1, 0, b[i], b[i + 1]);
      int hi = covering(hf, g0, g1, 1, b[i], b[i + 1]);
      if (lo < 0 && hi < 0) continue;
      if (out) {
        out[nf].axis = hf[g0].axis;
        out[nf].ipos = hf[g0].ipos;
        out[nf].it0 = b[i];
        out[nf].it1 = b[i + 1];
        out[nf].lo = lo;
        out[nf].hi = hi;
      }
      nf++;
    }
    free(b);
    g0 = g1;
  }
  return nf;
}

int tr2_mesh_finish(tr2_mesh *m) {
  if (m->ncell < 1) return -1;

  /* --- замощение ровно один раз, целочисленно --- */
  int *cov = calloc((size_t)m->nx * (size_t)m->ny, sizeof *cov);
  if (!cov) return -1;
  for (int c = 0; c < m->ncell; c++)
    for (int ix = m->cell[c].ix0; ix < m->cell[c].ix1; ix++)
      for (int iy = m->cell[c].iy0; iy < m->cell[c].iy1; iy++)
        cov[ix * m->ny + iy]++;
  int bad = 0;
  for (int i = 0; i < m->nx * m->ny; i++)
    if (cov[i] != 1) bad = 1;
  free(cov);
  if (bad) return -1;

  /* --- полуграни --- */
  int nhf = 4 * m->ncell;
  halfface *hf = calloc((size_t)nhf, sizeof *hf);
  if (!hf) return -1;
  for (int c = 0; c < m->ncell; c++) {
    const tr2_cell *k = &m->cell[c];
    halfface *h = &hf[4 * c];
    h[0] = (halfface){0, k->ix1, k->iy0, k->iy1, c, 0}; /* верхняя грань по x */
    h[1] = (halfface){0, k->ix0, k->iy0, k->iy1, c, 1}; /* нижняя грань по x  */
    h[2] = (halfface){1, k->iy1, k->ix0, k->ix1, c, 0};
    h[3] = (halfface){1, k->iy0, k->ix0, k->ix1, c, 1};
  }
  qsort(hf, (size_t)nhf, sizeof *hf, hf_cmp);

  int nf = build_faces(hf, nhf, NULL);
  if (nf < 0) {
    free(hf);
    return -1;
  }
  tr2_face *face = calloc((size_t)nf, sizeof *face);
  int *coff = calloc((size_t)m->ncell + 1, sizeof *coff);
  if (!face || !coff) {
    free(hf);
    free(face);
    free(coff);
    return -1;
  }
  if (build_faces(hf, nhf, face) != nf) {
    free(hf);
    free(face);
    free(coff);
    return -1;
  }
  free(hf);

  /* --- CSR граней по ячейкам --- */
  for (int f = 0; f < nf; f++) {
    if (face[f].lo >= 0) coff[face[f].lo + 1]++;
    if (face[f].hi >= 0) coff[face[f].hi + 1]++;
  }
  for (int c = 0; c < m->ncell; c++)
    coff[c + 1] += coff[c];
  int *cfac = calloc((size_t)coff[m->ncell], sizeof *cfac);
  int *fill = calloc((size_t)m->ncell, sizeof *fill);
  if (!cfac || !fill) {
    free(face);
    free(coff);
    free(cfac);
    free(fill);
    return -1;
  }
  for (int f = 0; f < nf; f++) {
    if (face[f].lo >= 0) cfac[coff[face[f].lo] + fill[face[f].lo]++] = f;
    if (face[f].hi >= 0) cfac[coff[face[f].hi] + fill[face[f].hi]++] = f;
  }
  free(fill);

  m->nface = nf;
  m->face = face;
  m->coff = coff;
  m->cfac = cfac;
  return 0;
}

int tr2_mesh_nsplit(const tr2_mesh *m) {
  int ns = 0;
  for (int c = 0; c < m->ncell; c++) {
    const tr2_cell *k = &m->cell[c];
    int len[2] = {k->ix1 - k->ix0, k->iy1 - k->iy0};
    for (int j = m->coff[c]; j < m->coff[c + 1]; j++) {
      const tr2_face *f = &m->face[m->cfac[j]];
      if (f->it1 - f->it0 < len[1 - f->axis]) ns++;
    }
  }
  return ns;
}

int tr2_mesh_uniform(tr2_mesh *m, int nx, int ny, double ux, double uy) {
  if (tr2_mesh_init(m, nx, ny, ux, uy, nx * ny) != 0) return -1;
  for (int ix = 0; ix < nx; ix++)
    for (int iy = 0; iy < ny; iy++)
      if (tr2_mesh_add(m, ix, iy, ix + 1, iy + 1) != 0) return -1;
  return tr2_mesh_finish(m);
}

/* Желаемый размер блока: чем дальше от центра, тем крупнее. Целочисленно. */
static int want_size(int nunit, int ix, int iy, int s) {
  int cx = nunit / 2, cy = nunit / 2;
  int d = abs(ix + s / 2 - cx) + abs(iy + s / 2 - cy);
  int t = 1;
  while (t * 4 <= d && t < nunit)
    t *= 2;
  return t;
}

static int grade_rec(tr2_mesh *m, int nunit, int ix, int iy, int s) {
  if (s > 1 && s > want_size(nunit, ix, iy, s)) {
    int h = s / 2;
    if (grade_rec(m, nunit, ix, iy, h) != 0) return -1;
    if (grade_rec(m, nunit, ix + h, iy, h) != 0) return -1;
    if (grade_rec(m, nunit, ix, iy + h, h) != 0) return -1;
    if (grade_rec(m, nunit, ix + h, iy + h, h) != 0) return -1;
    return 0;
  }
  return tr2_mesh_add(m, ix, iy, ix + s, iy + s);
}

int tr2_mesh_graded(tr2_mesh *m, int nunit, double ux, double uy) {
  if (tr2_mesh_init(m, nunit, nunit, ux, uy, nunit * nunit) != 0) return -1;
  if (grade_rec(m, nunit, 0, 0, nunit) != 0) return -1;
  return tr2_mesh_finish(m);
}

int tr2_mesh_aniso(tr2_mesh *m, int nunit, double ux, double uy, int w) {
  if (w < 1 || nunit % (2 * w) != 0) return -1;
  if (tr2_mesh_init(m, nunit, nunit, ux, uy, nunit * nunit) != 0) return -1;
  int half = nunit / 2;
  for (int ix = 0; ix < half; ix += w) /* слева: вытянуты по x */
    for (int iy = 0; iy < nunit; iy++)
      if (tr2_mesh_add(m, ix, iy, ix + w, iy + 1) != 0) return -1;
  for (int ix = half; ix < nunit; ix++) /* справа: вытянуты по y */
    for (int iy = 0; iy < nunit; iy += w)
      if (tr2_mesh_add(m, ix, iy, ix + 1, iy + w) != 0) return -1;
  return tr2_mesh_finish(m);
}

/* ============================================================== ядро DG1 */

/* Локальные координаты ячейки ξ,η ∈ [−1/2,1/2], L = c0 + c1·ξ + c2·η.
 * Базис ортогонален: ∫1 = A, ∫ξ = ∫η = ∫ξη = 0, ∫ξ² = ∫η² = A/12.
 * Слабая форма (поток, с противопоточным следом L̂ на гранях):
 *   −∫L(ω·∇v) + ∮(ω·n)L̂·v + ∫σLv = ∫qv,   v ∈ {1, ξ, η}
 * При v = 1 это ТОЧНОЕ уравнение баланса ячейки, поэтому сумма по ячейкам
 * телескопируется в глобальный баланс — отсюда консервативность. */

static int solve3(double a[3][3], const double b[3], double c[3]) {
  double mm[3][4];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++)
      mm[i][j] = a[i][j];
    mm[i][3] = b[i];
  }
  for (int k = 0; k < 3; k++) {
    int p = k;
    for (int i = k + 1; i < 3; i++)
      if (fabs(mm[i][k]) > fabs(mm[p][k])) p = i;
    if (fabs(mm[p][k]) <= 0.0) return -1;
    if (p != k)
      for (int j = k; j < 4; j++) {
        double t = mm[k][j];
        mm[k][j] = mm[p][j];
        mm[p][j] = t;
      }
    for (int i = k + 1; i < 3; i++) {
      double r = mm[i][k] / mm[k][k];
      for (int j = k; j < 4; j++)
        mm[i][j] -= r * mm[k][j];
    }
  }
  for (int i = 2; i >= 0; i--) {
    double s = mm[i][3];
    for (int j = i + 1; j < 3; j++)
      s -= mm[i][j] * c[j];
    c[i] = s / mm[i][i];
  }
  return 0;
}

/* Всё, что нужно про одну ячейку в одном направлении. */
typedef struct {
  double h[2], ctr[2], area;
  int i0[2], i1[2];
} cellgeo;

static void geom(const tr2_mesh *m, int c, cellgeo *g) {
  const tr2_cell *k = &m->cell[c];
  g->i0[0] = k->ix0;
  g->i1[0] = k->ix1;
  g->i0[1] = k->iy0;
  g->i1[1] = k->iy1;
  /* развёрнуто, а не циклом по осям: gcc -fanalyzer не видит, что цикл из двух
   * шагов заполняет ОБА элемента, и сообщает про h[1] как про неинициализированный */
  g->h[0] = (double)(k->ix1 - k->ix0) * m->ux;
  g->h[1] = (double)(k->iy1 - k->iy0) * m->uy;
  g->ctr[0] = 0.5 * (double)(k->ix0 + k->ix1) * m->ux;
  g->ctr[1] = 0.5 * (double)(k->iy0 + k->iy1) * m->uy;
  g->area = g->h[0] * g->h[1];
}

int tr2_solve(const tr2_problem *pr, tr2_stats *st) {
  const tr2_mesh *m = pr->mesh;
  const tr2_dirs *d = pr->dirs;
  const int nc = m->ncell, nd = d->n;
  const double p0 = pr->phase->p0;
  const double tol = pr->tol > 0.0 ? pr->tol : TR2_TOL;
  const int maxit = pr->one_pass ? 1 : (pr->maxit > 0 ? pr->maxit : TR2_MAXIT);
  double u[2] = {m->ux, m->uy};

  memset(st, 0, sizeof *st);
  double *phi = calloc((size_t)nc * 3, sizeof *phi);
  double *phin = calloc((size_t)nc * 3, sizeof *phin);
  double *gphi = calloc((size_t)nc * 3, sizeof *gphi);
  double *lc = calloc((size_t)nc * 3, sizeof *lc);
  int *indeg = calloc((size_t)nc, sizeof *indeg);
  int *order = calloc((size_t)nc, sizeof *order);
  if (!phi || !phin || !gphi || !lc || !indeg || !order) {
    free(phi);
    free(phin);
    free(gphi);
    free(lc);
    free(indeg);
    free(order);
    return -1;
  }
  if (pr->phi0) memcpy(phi, pr->phi0, (size_t)nc * 3 * sizeof *phi);

  /* влёт через границу считается один раз: от решения не зависит */
  double pin = 0.0;
  for (int f = 0; f < m->nface; f++) {
    const tr2_face *fc = &m->face[f];
    if (fc->lo >= 0 && fc->hi >= 0) continue;
    int a = fc->axis, t = 1 - a;
    int side = 2 * a + (fc->hi < 0 ? 1 : 0);
    double len = (double)(fc->it1 - fc->it0) * u[t];
    double tmid = 0.5 * (double)(fc->it0 + fc->it1) * u[t];
    double lbar = pr->lb[side] + pr->lbt[side] * tmid; /* среднее по грани */
    /* внешняя нормаль области на этой грани: +e_a если ячейка снизу (hi<0) */
    double nsign = fc->hi < 0 ? 1.0 : -1.0;
    for (int mi = 0; mi < nd; mi++) {
      double on = nsign * (a == 0 ? d->ox[mi] : d->oy[mi]);
      if (on < 0.0) pin += d->w[mi] * (-on) * lbar * len;
    }
  }
  st->pin = pin;

  int it = 0;
  double resid = 0.0;
  for (it = 0; it < maxit; it++) {
    memset(phin, 0, (size_t)nc * 3 * sizeof *phin);
    memset(gphi, 0, (size_t)nc * 3 * sizeof *gphi);
    st->pout = 0.0;
    st->pemit = 0.0;
    st->nlimited = 0;
    st->nneg = 0;
    st->lmin = INFINITY;
    st->lmax = -INFINITY;
    st->cmin = INFINITY;

    for (int mi = 0; mi < nd; mi++) {
      double o[2] = {d->ox[mi], d->oy[mi]};
      double gm = pr->gdir ? pr->gdir[mi] : 1.0;

      /* --- порядок обхода: топологическая сортировка по потоку (Кан) --- */
      memset(indeg, 0, (size_t)nc * sizeof *indeg);
      for (int f = 0; f < m->nface; f++) {
        const tr2_face *fc = &m->face[f];
        if (fc->lo < 0 || fc->hi < 0) continue;
        if (o[fc->axis] > 0.0)
          indeg[fc->hi]++;
        else
          indeg[fc->lo]++;
      }
      int head = 0, tail = 0;
      for (int c = 0; c < nc; c++)
        if (indeg[c] == 0) order[tail++] = c;
      while (head < tail) {
        int c = order[head++];
        for (int j = m->coff[c]; j < m->coff[c + 1]; j++) {
          const tr2_face *fc = &m->face[m->cfac[j]];
          if (fc->lo < 0 || fc->hi < 0) continue;
          int down = o[fc->axis] > 0.0 ? fc->hi : fc->lo;
          if (down == c) continue; /* эта грань втекает в c, а не вытекает */
          if (--indeg[down] == 0) order[tail++] = down;
        }
      }
      if (tail != nc) { /* цикл обхода: развёртка неприменима, молча считать нельзя */
        free(phi);
        free(phin);
        free(gphi);
        free(lc);
        free(indeg);
        free(order);
        return -2;
      }

      /* --- сама развёртка --- */
      for (int oi = 0; oi < nc; oi++) {
        int c = order[oi];
        cellgeo g;
        geom(m, c, &g);
        double sig = gm * pr->sig_t[c];
        double mtx[3][3] = {{0}}, rhs[3] = {0}, cf[3] = {0};

        /* объёмные члены */
        mtx[0][0] += sig * g.area;
        mtx[1][1] += sig * g.area / 12.0;
        mtx[2][2] += sig * g.area / 12.0;
        mtx[1][0] += -(o[0] / g.h[0]) * g.area;
        mtx[2][0] += -(o[1] / g.h[1]) * g.area;
        double ss = pr->sig_s[c] * p0;
        double q[3] = {ss * phi[3 * c], ss * phi[3 * c + 1], ss * phi[3 * c + 2]};
        if (pr->emit) {
          const double *e = &pr->emit[((size_t)c * (size_t)nd + (size_t)mi) * 3];
          q[0] += e[0];
          q[1] += e[1];
          q[2] += e[2];
        }
        rhs[0] += g.area * q[0];
        rhs[1] += g.area / 12.0 * q[1];
        rhs[2] += g.area / 12.0 * q[2];
        st->pemit += d->w[mi] * g.area *
                     (pr->emit ? pr->emit[((size_t)c * (size_t)nd + (size_t)mi) * 3] : 0.0);

        for (int j = m->coff[c]; j < m->coff[c + 1]; j++) {
          int fi = m->cfac[j];
          const tr2_face *fc = &m->face[fi];
          int a = fc->axis, t = 1 - a;
          int own_lo = (fc->lo == c);
          double wn = own_lo ? o[a] : -o[a];
          double xif = own_lo ? 0.5 : -0.5;
          double s0 = (double)(fc->it0 - g.i0[t]) / (double)(g.i1[t] - g.i0[t]) - 0.5;
          double s1 = (double)(fc->it1 - g.i0[t]) / (double)(g.i1[t] - g.i0[t]) - 0.5;
          double E0 = s1 - s0;
          double E1 = 0.5 * (s1 * s1 - s0 * s0);
          double E2 = (s1 * s1 * s1 - s0 * s0 * s0) / 3.0;
          double sc = wn * g.h[t];

          if (wn > 0.0) { /* выток: след свой, идёт в матрицу */
            mtx[0][0] += sc * E0;
            mtx[0][1 + a] += sc * xif * E0;
            mtx[0][1 + t] += sc * E1;
            mtx[1 + a][0] += sc * xif * E0;
            mtx[1 + a][1 + a] += sc * xif * xif * E0;
            mtx[1 + a][1 + t] += sc * xif * E1;
            mtx[1 + t][0] += sc * E1;
            mtx[1 + t][1 + a] += sc * xif * E1;
            mtx[1 + t][1 + t] += sc * E2;
          } else { /* вток: след известен, L̂ = P + Q·τ в СВОИХ координатах */
            double P, Q;
            int nb = own_lo ? fc->hi : fc->lo;
            if (nb < 0) {
              int side = 2 * a + (own_lo ? 1 : 0);
              /* L̂ = lb + lbt·t_физ, а t_физ = ctr[t] + τ·h[t] ⇒ P берётся при
               * τ = 0, то есть в ЦЕНТРЕ ячейки по тангенциали, а не в центре
               * грани: интегралы E0,E1,E2 живут в τ, отсчитанном от центра. */
              P = pr->lb[side] + pr->lbt[side] * g.ctr[t];
              Q = pr->lbt[side] * g.h[t];
            } else {
              cellgeo gn;
              geom(m, nb, &gn);
              double beta = g.h[t] / gn.h[t];
              double alpha = (g.ctr[t] - gn.ctr[t]) / gn.h[t];
              const double *dc = &lc[3 * nb];
              P = dc[0] + dc[1 + a] * (-xif) + dc[1 + t] * alpha;
              Q = dc[1 + t] * beta;
            }
            double J0 = P * E0 + Q * E1, J1 = P * E1 + Q * E2;
            rhs[0] -= sc * J0;
            rhs[1 + a] -= sc * xif * J0;
            rhs[1 + t] -= sc * J1;
          }
        }

        if (solve3(mtx, rhs, cf) != 0) {
          free(phi);
          free(phin);
          free(gphi);
          free(lc);
          free(indeg);
          free(order);
          return -3;
        }

        /* --- ограничитель: ТОЛЬКО ПОЛОЖИТЕЛЬНОСТЬ ---
         * Обычный ограничитель монотонности берёт границы со ВСЕХ соседей, но
         * в развёртке доступны только ВЕРХНИЕ ПО ПОТОКУ, а односторонняя
         * граница обнуляет наклон у любого гладкого поля (у растущего вдоль
         * потока поля все верхние соседи ниже, и запас вверх равен нулю).
         * Нижние соседи по потоку — это поле ПРЕДЫДУЩЕЙ итерации для ЭТОГО
         * направления, то есть полные ND ординат в памяти, что запрещено Р1.
         * Значит в развёртке остаётся положительность — и это ровно то, чего
         * требуют К6/T5: около границы тени радианс не должен уходить в минус.
         * Константу она сохраняет тривиально (наклон нулевой, срезать нечего). */
        double dev0 = 0.5 * (fabs(cf[1]) + fabs(cf[2])); /* минимум по углам = c0 − dev0 */
        if (pr->lim != TR2_LIM_OFF && dev0 > 0.0) {
          double alpha = 1.0;
          if (pr->lim == TR2_LIM_CONSERVATIVE) {
            /* НАКЛОН И СРЕДНЕЕ ИЩУТСЯ СОВМЕСТНО, А НЕ ПО ОЧЕРЕДИ. Срезав
             * наклон и пересчитав среднее из баланса, положительность
             * приходится проверять заново, и поочерёдные проходы к ней лишь
             * ПРИБЛИЖАЮТСЯ: измерено, что после трёх проходов оставалось
             * −3.4e-6 от L_b — на десять порядков выше округления. Но среднее
             * есть АФФИННАЯ функция коэффициента срезки α:
             *     c0(α) = (rhs0 − α·K)/M00,   K = M01·c1 + M02·c2,
             * а требование «минимум по углам ≥ 0» это c0(α) ≥ α·dev0, откуда
             *     α ≤ rhs0 / (K + dev0·M00).
             * Замкнутая форма, один шаг — и баланс, и положительность точны
             * ОБА. Порога при этом не появилось: α берётся из уравнения. */
            double K = mtx[0][1] * cf[1] + mtx[0][2] * cf[2];
            double den = K + dev0 * mtx[0][0];
            if (den > 0.0) {
              double amax = rhs[0] / den;
              if (amax < alpha) alpha = amax;
            }
            if (alpha < 0.0) alpha = 0.0;
            if (alpha < 1.0) {
              cf[1] *= alpha;
              cf[2] *= alpha;
              cf[0] = (rhs[0] - mtx[0][1] * cf[1] - mtx[0][2] * cf[2]) / mtx[0][0];
              st->nlimited++;
            }
          } else { /* SLOPE: срезка без пересчёта среднего ⇒ баланс ломается */
            if (cf[0] <= 0.0)
              alpha = 0.0;
            else if (cf[0] < dev0)
              alpha = cf[0] / dev0;
            if (alpha < 1.0) {
              cf[1] *= alpha;
              cf[2] *= alpha;
              st->nlimited++;
            }
          }
        }

        lc[3 * c] = cf[0];
        lc[3 * c + 1] = cf[1];
        lc[3 * c + 2] = cf[2];
        for (int k = 0; k < 3; k++) {
          phin[3 * c + k] += d->w[mi] * cf[k];
          gphi[3 * c + k] += d->w[mi] * gm * cf[k];
        }
        if (cf[0] < st->lmin) st->lmin = cf[0];
        if (cf[0] > st->lmax) st->lmax = cf[0];
        double corner = cf[0] - 0.5 * (fabs(cf[1]) + fabs(cf[2]));
        if (corner < st->cmin) st->cmin = corner;
        if (corner < 0.0) st->nneg++;

        /* выток через границу области */
        for (int j = m->coff[c]; j < m->coff[c + 1]; j++) {
          const tr2_face *fc = &m->face[m->cfac[j]];
          if (fc->lo >= 0 && fc->hi >= 0) continue;
          int a = fc->axis, t = 1 - a;
          int own_lo = (fc->lo == c);
          double wn = own_lo ? o[a] : -o[a];
          if (wn <= 0.0) continue;
          double xif = own_lo ? 0.5 : -0.5;
          double s0 = (double)(fc->it0 - g.i0[t]) / (double)(g.i1[t] - g.i0[t]) - 0.5;
          double s1 = (double)(fc->it1 - g.i0[t]) / (double)(g.i1[t] - g.i0[t]) - 0.5;
          double E0 = s1 - s0, E1 = 0.5 * (s1 * s1 - s0 * s0);
          st->pout += d->w[mi] * wn * g.h[t] * ((cf[0] + cf[1 + a] * xif) * E0 + cf[1 + t] * E1);
        }
      }
    }

    /* поглощение: σ_t от НОВОГО прохода, σ_s·φ от ТОГО ЖЕ φ, что дал источник */
    double pabs = 0.0;
    for (int c = 0; c < nc; c++) {
      cellgeo g;
      geom(m, c, &g);
      pabs += g.area * (pr->sig_t[c] * gphi[3 * c] - pr->sig_s[c] * phi[3 * c]);
    }
    st->pabs = pabs;
    st->balance = st->pin + st->pemit - st->pout - st->pabs;

    double dn = 0.0, nn = 0.0;
    for (int i = 0; i < 3 * nc; i++) {
      double dd = phin[i] - phi[i];
      dn += dd * dd;
      nn += phin[i] * phin[i];
    }
    resid = nn > 0.0 ? sqrt(dn / nn) : 0.0;
    memcpy(phi, phin, (size_t)nc * 3 * sizeof *phi);
    if (resid < tol) {
      it++;
      break;
    }
  }
  st->iters = it;
  st->resid = resid;
  st->phi = phi;
  free(phin);
  free(gphi);
  free(lc);
  free(indeg);
  free(order);
  return 0;
}

void tr2_stats_free(tr2_stats *st) {
  free(st->phi);
  st->phi = NULL;
}
