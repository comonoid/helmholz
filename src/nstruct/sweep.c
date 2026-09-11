/* sweep.c — свип на новых массивах: §835 (скалярный луч), §836-§837
 * (колонки/линии), §838 (материализованный обход + трафик). См. sweep.h.
 */
#include "sweep.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define HZ_SW_W (4.0 * 3.14159265358979323846 / 26.0) /* вес направления, изотропия 4π/N */
#define HZ_SW_SLABBITS 21 /* слэб и позиция в ключе; оси разумных сцен ≤ 2^21 */
#define HZ_SW_ND 26

/* знаки 26 направлений: 6 осей, 12 рёбер (√2), 8 углов (√3); ЦЕЛЫЕ —
 * сравнение с нулём без float-equal, нормировка отдельной таблицей */
static const int sw_sv[HZ_SW_ND][3] = {
    {1, 0, 0},  {-1, 0, 0},  {0, 1, 0},   {0, -1, 0},  {0, 0, 1},   {0, 0, -1}, {1, 1, 0},
    {1, -1, 0}, {-1, 1, 0},  {-1, -1, 0}, {1, 0, 1},   {1, 0, -1},  {-1, 0, 1}, {-1, 0, -1},
    {0, 1, 1},  {0, 1, -1},  {0, -1, 1},  {0, -1, -1}, {1, 1, 1},   {1, 1, -1}, {1, -1, 1},
    {-1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}};
static const double sw_norm[HZ_SW_ND] = {1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.7320508075688772,
                                         1.7320508075688772,
                                         1.7320508075688772,
                                         1.7320508075688772,
                                         1.7320508075688772,
                                         1.7320508075688772,
                                         1.7320508075688772,
                                         1.7320508075688772};

typedef struct {
  int64_t key;
  int32_t li;
} hz_sw_pair;

static int sw_cmp_pair(const void *a, const void *b) {
  int64_t x = ((const hz_sw_pair *)a)->key, y = ((const hz_sw_pair *)b)->key;
  return (x > y) - (x < y);
}

/* линия и слэб клетки id на направление d. ЛИНИЯ — поперечные разности
 * индексов со знаками (А1481); слэб — знаковая проекция, в линии строг.
 * Знак направления сидит в sv. Ключ умещается в int64 при осях ≤ 2^21. */
static void sw_line_slab(const hz_pyr *py, int d, int64_t id, int64_t *line, int64_t *slab) {
  const int *sv = sw_sv[d];
  int64_t ii[3], len[3];
  int na[3] = {0, 0, 0}, nnz = 0, ax;
  len[0] = py->nx;
  len[1] = py->ny;
  len[2] = py->nz;
  ii[0] = id % py->nx;
  ii[1] = (id / py->nx) % py->ny;
  ii[2] = id / ((int64_t)py->nx * py->ny);
  for (ax = 0; ax < 3; ax++)
    if (sv[ax] != 0) na[nnz++] = ax;
  if (nnz == 1) {
    int c = (na[0] + 1) % 3, dd = (na[0] + 2) % 3;
    *line = ii[c] + len[c] * ii[dd];
    *slab = (int64_t)sv[na[0]] * ii[na[0]];
  } else if (nnz == 2) {
    int a = na[0], b = na[1], c = 3 - a - b;
    /* a ≠ b по построению (оси уникальны); анализатор этого не выводит */
    if (c < 0) c = 0;
    if (c > 2) c = 2;
    int64_t o1 = (int64_t)sv[b] * ii[b] - (int64_t)sv[a] * ii[a] + len[a];
    *line = o1 + (len[a] + len[b]) * ii[c];
    *slab = (int64_t)sv[a] * ii[a] + (int64_t)sv[b] * ii[b];
  } else {
    int a = na[0], b = na[1], c = na[2];
    int64_t o1 = (int64_t)sv[b] * ii[b] - (int64_t)sv[a] * ii[a] + len[a];
    int64_t o2 = (int64_t)sv[c] * ii[c] - (int64_t)sv[b] * ii[b] + len[b];
    *line = o1 + (len[a] + len[b]) * o2;
    *slab = (int64_t)sv[a] * ii[a] + (int64_t)sv[b] * ii[b] + (int64_t)sv[c] * ii[c];
  }
}

/* ПОХОДНЫЙ МАССИВ направления (§838): листья в порядке обхода, куски —
 * в порядке обхода; cn и wn посчитаны ОДИН раз (А1488: память nd×nt —
 * шаг-измерение, постоянное решение — Morton-раскладка базы). */
typedef struct {
  int32_t npcs;
  int64_t first; /* офсет в wcn/wn/wpid */
  int64_t line;  /* линия этого листа; смена — тёмный вход луча */
} hz_sw_wleaf;

typedef struct {
  hz_sw_wleaf *leaf;
  double *wcn; /* area·|ω·n| */
  double *wn;  /* |ω·n| */
  double *kd;  /* альбедо куска — из походного слота: минус случайный поток */
  int32_t *wpid;
  int32_t n; /* листьев в походе — только С КУСКАМИ (§838: пустой лист в
              * модели слоя не делает ничего, его визит — чистые 24 Б) */
} hz_sw_walk;

static int sw_build_walk(const hz_pyr *py, int d, const int32_t *order, const double *area,
                         const double *nrm, const double *om, const double *kd, hz_sw_walk *w) {
  int32_t k, wk = 0;
  int64_t pos = 0;
  w->leaf = (hz_sw_wleaf *)malloc((size_t)py->nleaf * sizeof *w->leaf);
  w->wcn = (double *)malloc((size_t)py->nt * sizeof *w->wcn);
  w->wn = (double *)malloc((size_t)py->nt * sizeof *w->wn);
  w->kd = (double *)malloc((size_t)py->nt * sizeof *w->kd);
  w->wpid = (int32_t *)malloc((size_t)py->nt * sizeof *w->wpid);
  if (!w->leaf || !w->wcn || !w->wn || !w->kd || !w->wpid) return 2;
  for (k = 0; k < py->nleaf; k++) {
    const hz_pyr_leaf *lf = &py->leaf[order[k]];
    int32_t u;
    if (lf->npcs == 0) continue; /* пустой лист — не участник переноса */
    {
      int64_t slab_dummy;
      sw_line_slab(py, d, py->leaf_id[order[k]], &w->leaf[wk].line, &slab_dummy);
    }
    w->leaf[wk].npcs = lf->npcs;
    w->leaf[wk].first = pos;
    for (u = 0; u < lf->npcs; u++) {
      int32_t p = py->csr[lf->pcs_first + u];
      double d0 = om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                  om[2] * nrm[3 * (int64_t)p + 2];
      double an = fabs(d0);
      w->wcn[pos] = area[p] * an;
      w->wn[pos] = an;
      w->kd[pos] = kd[p];
      w->wpid[pos] = p;
      pos++;
    }
    wk++;
  }
  w->n = wk;
  return 0;
}

static void sw_free_walk(hz_sw_walk *w) {
  free(w->leaf);
  free(w->wcn);
  free(w->wn);
  free(w->kd);
  free(w->wpid);
  w->leaf = NULL;
  w->wcn = NULL;
  w->wn = NULL;
  w->kd = NULL;
  w->wpid = NULL;
}

int hz_sw_run(hz_pyr *py, int32_t nt, const double *area, const double *nrm, const double *kd,
              const hz_sw_opts *o, hz_sw_stat *st, double *e_hist) {
  double *Ed = NULL, *Eprev = NULL;
  int32_t **order = NULL;
  hz_sw_walk *walks = NULL;
  int32_t *vindex = NULL, *counts = NULL;
  uint64_t *bits = NULL;
  hz_sw_pair *pairs = NULL;
  int nd, d, it, rc = 0;
  double csec;

  memset(st, 0, sizeof *st);
  if (!py || !area || !nrm || !kd || !o) return 1;
  if (o->iters <= 0) return 1;
  if (!py->pcs || nt != py->nt) return 1;
  nd = (o->ndirs == 26) ? HZ_SW_ND : 6;
  csec = py->cell * py->cell; /* A_⊥ луча в клетке листа */

  Ed = (double *)calloc((size_t)nt, sizeof *Ed);
  Eprev = (double *)calloc((size_t)nt, sizeof *Eprev);
  vindex = (int32_t *)malloc((size_t)py->nleaf * sizeof *vindex);
  bits = (uint64_t *)malloc((size_t)((py->nleaf + 63) >> 6) * sizeof *bits);
  pairs = (hz_sw_pair *)malloc((size_t)py->nleaf * sizeof *pairs);
  order = (int32_t **)calloc((size_t)nd, sizeof *order);
  counts = (int32_t *)calloc((size_t)nd, sizeof *counts);
  walks = (hz_sw_walk *)calloc((size_t)nd, sizeof *walks);
  if (!Ed || !Eprev || !vindex || !bits || !pairs || !order || !counts || !walks) {
    rc = 2;
    goto done;
  }

  /* порядок листов на направление. scalar §835: марш (DFS пирамиды), 6 осей.
   * col §836-§838: ключ (линия, слэб) — сортировка. */
  for (d = 0; d < nd; d++) {
    int32_t li;
    order[d] = (int32_t *)malloc((size_t)py->nleaf * sizeof *order[d]);
    if (!order[d]) {
      rc = 2;
      goto done;
    }
    if (o->mode == 1) {
      for (li = 0; li < py->nleaf; li++) {
        int64_t line, slab;
        sw_line_slab(py, d, py->leaf_id[li], &line, &slab);
        pairs[li].key = line * ((int64_t)1 << HZ_SW_SLABBITS) + slab;
        pairs[li].li = li;
      }
      qsort(pairs, (size_t)py->nleaf, sizeof *pairs, sw_cmp_pair);
      for (li = 0; li < py->nleaf; li++)
        order[d][li] = pairs[li].li;
    } else {
      hz_pyr_march_stat ms;
      double om[3] = {(double)sw_sv[d][0] / sw_norm[d], (double)sw_sv[d][1] / sw_norm[d],
                      (double)sw_sv[d][2] / sw_norm[d]};
      hz_pyr_march(py, om, 0, &ms, bits, vindex);
      for (li = 0; li < py->nleaf; li++)
        order[d][vindex[li] - 1] = li;
    }
    counts[d] = (o->mode == 1) ? 0 : py->nleaf; /* col: заполнится после сборки похода */
  }

  /* §838: материализация походов (режим col) — cn/wn один раз на направление */
  if (o->mode == 1) {
    for (d = 0; d < nd; d++) {
      double om[3] = {(double)sw_sv[d][0] / sw_norm[d], (double)sw_sv[d][1] / sw_norm[d],
                      (double)sw_sv[d][2] / sw_norm[d]};
      rc = sw_build_walk(py, d, order[d], area, nrm, om, kd, &walks[d]);
      if (rc != 0) goto done;
      counts[d] = walks[d].n;
    }
  }

  for (it = 0; it < o->iters; it++) {
    double emitted = 0, absorbed = 0, lost = 0, recycled = 0, e_sum = 0, area_sum = 0;
    int64_t traffic = 0;
    int32_t p;
    memset(Ed, 0, (size_t)nt * sizeof *Ed);
    for (p = 0; p < nt; p++)
      Eprev[p] = py->pcs[p].e;
    for (d = 0; d < nd; d++) {
      double om[3] = {(double)sw_sv[d][0] / sw_norm[d], (double)sw_sv[d][1] / sw_norm[d],
                      (double)sw_sv[d][2] / sw_norm[d]};
      double L = 0.0; /* луч входит в домён тёмным */
      int64_t prevline = -1;
      int32_t k;
      if (o->mode == 1) {
        const hz_sw_walk *w = &walks[d];
        for (k = 0; k < counts[d]; k++) {
          const hz_sw_wleaf *lf = &w->leaf[k];
          double cntot = 0.0, tau, T, Lsurf = 0.0;
          int64_t e = lf->first;
          int32_t u;
          if (k + 8 < counts[d]) __builtin_prefetch(&w->leaf[k + 8]);
          if (w->leaf[k].line != prevline) { /* новая линия — тёмный вход */
            L = 0.0;
            prevline = w->leaf[k].line;
          }
          traffic += 24; /* походный лист: npcs 4 + first 8 + line 8 + паддинг */
          if (o->noprop) L = 0.0;
          for (u = 0; u < lf->npcs; u++, e++) {
            traffic += 28; /* cn 8 + wn 8 + p 4 + E r/w 8 */
            cntot += w->wcn[e];
          }
          tau = cntot / csec;
          T = o->tau0 ? 1.0 : exp(-tau);
          if (tau > 0.0) {
            e = lf->first;
            for (u = 0; u < lf->npcs; u++, e++) {
              p = w->wpid[e];
              double rho = o->rho < 0 ? w->kd[e] : o->rho;
              Lsurf += (o->le + rho * Eprev[p] / (2.0 * M_PI)) * w->wcn[e];
              Ed[p] += HZ_SW_W * L * (1.0 - T) * w->wn[e];
            }
            Lsurf /= (double)cntot;
            absorbed += HZ_SW_W * L * (1.0 - T) * csec;
            emitted += HZ_SW_W * o->le * (1.0 - T) * csec;            /* внешний вход Le */
            recycled += HZ_SW_W * (Lsurf - o->le) * (1.0 - T) * csec; /* ρ·E_prev */
            L = L * T + Lsurf * (1.0 - T);
          }
          if (o->noprop) L = 0.0;
        }
        lost += L * csec;
      } else {
        /* scalar §835 — эталонный контраст, старый путь по csr/нормалям */
        for (k = 0; k < counts[d]; k++) {
          const hz_pyr_leaf *lf = &py->leaf[order[d][k]];
          double cntot = 0.0, tau, T, Lsurf = 0.0;
          int32_t u;
          if (o->noprop) L = 0.0;
          for (u = 0; u < lf->npcs; u++) {
            p = py->csr[lf->pcs_first + u];
            double d0 = om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                        om[2] * nrm[3 * (int64_t)p + 2];
            cntot += area[p] * fabs(d0);
          }
          tau = cntot / csec;
          T = o->tau0 ? 1.0 : exp(-tau);
          if (tau > 0.0) {
            for (u = 0; u < lf->npcs; u++) {
              p = py->csr[lf->pcs_first + u];
              double rho = o->rho < 0 ? kd[p] : o->rho;
              double cnm = fabs(om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                                om[2] * nrm[3 * (int64_t)p + 2]);
              Lsurf += (o->le + rho * Eprev[p] / (2.0 * M_PI)) * area[p] * cnm;
              Ed[p] += HZ_SW_W * L * (1.0 - T) * cnm;
            }
            Lsurf /= cntot;
            absorbed += HZ_SW_W * L * (1.0 - T) * csec;
            emitted += HZ_SW_W * o->le * (1.0 - T) * csec;
            recycled += HZ_SW_W * (Lsurf - o->le) * (1.0 - T) * csec;
            L = L * T + Lsurf * (1.0 - T);
          }
          if (o->noprop) L = 0.0;
        }
        lost += L * csec;
      }
    }
    for (p = 0; p < nt; p++) {
      py->pcs[p].e = (float)Ed[p];
      e_sum += Ed[p] * area[p];
      area_sum += area[p];
    }
    st->e_avg = e_sum / area_sum;
    st->emitted = emitted;
    st->recycled = recycled;
    st->absorbed = absorbed;
    st->lost = lost;
    st->traffic = traffic;
    if (e_hist) e_hist[it] = st->e_avg;
  }

done:
  free(Ed);
  free(Eprev);
  free(vindex);
  free(bits);
  free(pairs);
  if (walks) {
    for (d = 0; d < nd; d++)
      sw_free_walk(&walks[d]);
    free(walks);
  }
  if (order) {
    int dd;
    for (dd = 0; dd < nd; dd++)
      free((void *)order[dd]);
    free(order);
  }
  free(counts);
  return rc;
}
