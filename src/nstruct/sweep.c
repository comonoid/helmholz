/* sweep.c — свип на новых массивах: §835 (скалярный луч), §836-§838
 * (колонки/линии, материализованный обход), §841 (гигиена цикла).
 * См. sweep.h.
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

/* линия и слэб клетки id на направление d (А1481). Знак направления в sv. */
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

/* ПОХОДНЫЙ ЛИСТ направления (§841): cntot и T предвычислены — итерация
 * однопроходна. Только кусковые листья. */
typedef struct {
  int32_t npcs;
  int64_t first; /* офсет в wcn/wn/kd/wpid */
  int64_t line;  /* линия; смена — тёмный вход луча */
  double cntot;  /* Σ area·|ω·n| по кускам листа */
  double T;      /* exp(−τ) */
} hz_sw_wleaf;

typedef struct {
  hz_sw_wleaf *leaf;
  double *wcn, *wn, *kd;
  int32_t *wpid;
  int32_t n; /* походных листьев (с кусками) */
} hz_sw_walk;

static void sw_free_walk(hz_sw_walk *w) {
  free(w->leaf);
  free(w->wcn);
  free(w->wn);
  free(w->kd);
  free(w->wpid);
  memset(w, 0, sizeof *w);
}

/* §841: поход строится ПРЯМО — фильтр кусковых листьев, сортировка ТОЛЬКО
 * их (полная сортировка 14.9M пустых не нужна вовсе), cntot/T сразу. */
static int sw_build_walk(const hz_pyr *py, int d, const double *area, const double *nrm,
                         const double *kd, const double *om, int tau0, hz_sw_walk *w) {
  int32_t li, k;
  int64_t pos = 0;
  hz_sw_pair *pairs = (hz_sw_pair *)malloc((size_t)py->nleaf * sizeof *pairs);
  w->leaf = (hz_sw_wleaf *)malloc((size_t)py->nleaf * sizeof *w->leaf);
  w->wcn = (double *)malloc((size_t)py->nt * sizeof *w->wcn);
  w->wn = (double *)malloc((size_t)py->nt * sizeof *w->wn);
  w->kd = (double *)malloc((size_t)py->nt * sizeof *w->kd);
  w->wpid = (int32_t *)malloc((size_t)py->nt * sizeof *w->wpid);
  if (!pairs || !w->leaf || !w->wcn || !w->wn || !w->kd || !w->wpid) {
    free(pairs);
    sw_free_walk(w);
    return 2;
  }
  {
    int32_t np = 0;
    for (li = 0; li < py->nleaf; li++) {
      if (py->leaf[li].npcs == 0) continue; /* пустой лист — не участник переноса */
      sw_line_slab(py, d, py->leaf_id[li], &pairs[np].key, &(int64_t){0});
      pairs[np].li = li;
      np++;
    }
    w->n = np;
  }
  qsort(pairs, (size_t)w->n, sizeof *pairs, sw_cmp_pair);
  for (k = 0; k < w->n; k++) {
    const hz_pyr_leaf *lf = &py->leaf[pairs[k].li];
    int32_t u;
    double cntot = 0.0;
    {
      int64_t slab_dummy;
      sw_line_slab(py, d, py->leaf_id[pairs[k].li], &w->leaf[k].line, &slab_dummy);
    }
    w->leaf[k].npcs = lf->npcs;
    w->leaf[k].first = pos;
    for (u = 0; u < lf->npcs; u++) {
      int32_t p = py->csr[lf->pcs_first + u];
      double d0 = om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                  om[2] * nrm[3 * (int64_t)p + 2];
      double an = fabs(d0);
      w->wcn[pos] = area[p] * an;
      w->wn[pos] = an;
      w->kd[pos] = kd[p];
      w->wpid[pos] = p;
      cntot += w->wcn[pos];
      pos++;
    }
    w->leaf[k].cntot = cntot;
    w->leaf[k].T = (tau0 || cntot <= 0.0) ? 1.0 : exp(-cntot / (py->cell * py->cell));
  }
  free(pairs);
  return 0;
}

int hz_sw_run(hz_pyr *py, int32_t nt, const double *area, const double *nrm, const double *kd,
              const hz_sw_opts *o, hz_sw_stat *st, double *e_hist) {
  double *Ed = NULL, *Eprev = NULL;
  int32_t *order = NULL, *vindex = NULL;
  uint64_t *bits = NULL;
  hz_sw_walk *walks = NULL;
  int nd, d, it, rc = 0;
  double csec;

  memset(st, 0, sizeof *st);
  if (!py || !area || !nrm || !kd || !o) return 1;
  if (o->iters <= 0) return 1;
  if (!py->pcs || nt != py->nt) return 1;
  nd = (o->ndirs == 26) ? HZ_SW_ND : 6;
  csec = py->cell * py->cell;

  Ed = (double *)calloc((size_t)nt, sizeof *Ed);
  Eprev = (double *)calloc((size_t)nt, sizeof *Eprev);
  order = (int32_t *)calloc(
      (size_t)py->nleaf,
      sizeof *order); /* calloc: анализатор видит инициализацию (FP-класс diam 07-24) */
  vindex = (int32_t *)malloc((size_t)py->nleaf * sizeof *vindex);
  bits = (uint64_t *)malloc((size_t)((py->nleaf + 63) >> 6) * sizeof *bits);
  walks = (hz_sw_walk *)calloc((size_t)nd, sizeof *walks);
  if (!Ed || !Eprev || !order || !vindex || !bits || !walks) {
    rc = 2;
    goto done;
  }

  if (o->mode == 1) {
    for (d = 0; d < nd; d++) {
      double om[3] = {(double)sw_sv[d][0] / sw_norm[d], (double)sw_sv[d][1] / sw_norm[d],
                      (double)sw_sv[d][2] / sw_norm[d]};
      rc = sw_build_walk(py, d, area, nrm, kd, om, o->tau0, &walks[d]);
      if (rc != 0) goto done;
    }
  } else {
    /* scalar §835: марш по всем листьям, 6 осей */
    {
      int32_t li2;
      for (li2 = 0; li2 < py->nleaf; li2++)
        order[li2] = li2; /* анализатор: полная инициализация */
    }
    hz_pyr_march_stat ms;
    double om[3] = {1.0, 0.0, 0.0};
    int32_t li;
    hz_pyr_march(py, om, 0, &ms, bits, vindex);
    for (li = 0; li < py->nleaf; li++)
      order[vindex[li] - 1] = li;
  }

  for (it = 0; it < o->iters; it++) {
    double emitted = 0, absorbed = 0, lost = 0, recycled = 0, e_sum = 0, area_sum = 0;
    int32_t p;
    memset(Ed, 0, (size_t)nt * sizeof *Ed);
    for (p = 0; p < nt; p++)
      Eprev[p] = py->pcs[p].e;
    for (d = 0; d < nd; d++) {
      double om[3] = {(double)sw_sv[d][0] / sw_norm[d], (double)sw_sv[d][1] / sw_norm[d],
                      (double)sw_sv[d][2] / sw_norm[d]};
      double L = 0.0;
      int64_t prevline = -1;
      if (o->mode == 1) {
        const hz_sw_walk *w = &walks[d];
        int32_t k;
        for (k = 0; k < w->n; k++) {
          const hz_sw_wleaf *lf = &w->leaf[k];
          double Lsurf = 0.0, tau0T = lf->T;
          int64_t e = lf->first, eend = lf->first + lf->npcs;
          int64_t pf = lf->first + 8 < nt ? lf->first + 8 : nt - 1;
          if (lf->line != prevline) { /* новая линия — тёмный вход */
            L = 0.0;
            prevline = lf->line;
          }
          if (o->noprop) L = 0.0;
          if (lf->cntot <= 0.0) continue; /* куски встык лучу — слоя нет */
          __builtin_prefetch(&w->leaf[k + 8 < w->n ? k + 8 : w->n - 1]);
          __builtin_prefetch(&Eprev[w->wpid[pf]]);
          for (e = lf->first; e < eend; e++) {
            p = w->wpid[e];
            double rho = o->rho < 0 ? w->kd[e] : o->rho;
            Lsurf += (o->le + rho * Eprev[p] / (2.0 * M_PI)) * w->wcn[e];
            Ed[p] += HZ_SW_W * L * (1.0 - tau0T) * w->wn[e];
          }
          Lsurf /= lf->cntot;
          absorbed += HZ_SW_W * L * (1.0 - tau0T) * csec;
          emitted += HZ_SW_W * o->le * (1.0 - tau0T) * csec;
          recycled += HZ_SW_W * (Lsurf - o->le) * (1.0 - tau0T) * csec;
          L = L * tau0T + Lsurf * (1.0 - tau0T);
          if (o->noprop) L = 0.0;
        }
        lost += L * csec;
      } else {
        int32_t k;
        for (k = 0; k < py->nleaf; k++) {
          const hz_pyr_leaf *lf = &py->leaf[order[k]];
          double cntot = 0.0, T, Lsurf = 0.0;
          int32_t u;
          if (o->noprop) L = 0.0;
          for (u = 0; u < lf->npcs; u++) {
            p = py->csr[lf->pcs_first + u];
            double d0 = om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                        om[2] * nrm[3 * (int64_t)p + 2];
            cntot += area[p] * fabs(d0);
          }
          T = o->tau0 ? 1.0 : exp(-cntot / csec);
          if (cntot > 0.0) {
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
    st->traffic = (int64_t)nd * ((int64_t)nt * 36 + (int64_t)walks[0].n * 40);
    if (e_hist) e_hist[it] = st->e_avg;
  }

done:
  free(Ed);
  free(Eprev);
  free(order);
  free(vindex);
  free(bits);
  if (walks) {
    for (d = 0; d < nd; d++)
      sw_free_walk(&walks[d]);
    free(walks);
  }
  return rc;
}
