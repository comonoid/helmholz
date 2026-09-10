/* sweep.c — свип на новых массивах: §835 (скалярный луч) и §836 (колонки).
 * См. sweep.h.
 *
 * ЭНЕРГЕТИКА ИТЕРАЦИИ (тождество схемы): на направлении луч входит тёмным,
 * на каждом слое излучённый поток w·Lsurf·(1−T)·A_⊥ равен поглощённому
 * w·L·(1−T)·A_⊥ плюс прирост луча; остаток за последним листом — потери.
 *
 * §836, режим col: занятые листья направления группируются в КОЛОНКИ
 * (поперечная пара координат), внутри колонки — обход по оси со знаком;
 * L — скаляр колонки. Поперечные колонки независимы: кусок получает свет
 * только своих колонок-направлений. Связь направлений — через E_prev
 * (источниковые итерации, Жакоби). Марш §834 в этом режиме не участвует —
 * порядок строится сортировкой; марш остаётся машинерией диагоналей и
 * агрегатов.
 */
#include "sweep.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#define HZ_SW_W 2.09439510239319549 /* 2π/3 — вес осевого направления, изотропия 4π/6 */
#define HZ_SW_POSBITS 21            /* позиций на оси < 2^21 — потолок сетки листа 2^23 */

typedef struct {
  int64_t key;
  int32_t li;
} hz_sw_pair;

static int sw_cmp_pair(const void *a, const void *b) {
  int64_t x = ((const hz_sw_pair *)a)->key, y = ((const hz_sw_pair *)b)->key;
  return (x > y) - (x < y);
}

int hz_sw_run(hz_pyr *py, int32_t nt, const double *area, const double *nrm, const double *kd,
              const hz_sw_opts *o, hz_sw_stat *st, double *e_hist) {
  static const double dirs[6][3] = {{1, 0, 0},  {-1, 0, 0}, {0, 1, 0},
                                    {0, -1, 0}, {0, 0, 1},  {0, 0, -1}};
  double *Ed = NULL, *Eprev = NULL, *cn = NULL;
  int32_t **order = NULL;
  hz_sw_pair *pairs = NULL;
  int32_t *vindex = NULL, *counts = NULL;
  uint64_t *bits = NULL;
  int d, it, rc = 0;
  double csec;

  memset(st, 0, sizeof *st);
  if (!py || !area || !nrm || !kd || !o || o->iters <= 0) return 1;
  if (!py->pcs || nt != py->nt) return 1;
  csec = py->cell * py->cell; /* A_⊥ осевого луча в клетке листа */

  Ed = (double *)calloc((size_t)nt, sizeof *Ed);
  Eprev = (double *)calloc((size_t)nt, sizeof *Eprev);
  cn = (double *)malloc((size_t)nt * sizeof *cn); /* area·|ω·n| на направление, scratch */
  vindex = (int32_t *)malloc((size_t)py->nleaf * sizeof *vindex);
  bits = (uint64_t *)malloc((size_t)((py->nleaf + 63) >> 6) * sizeof *bits);
  pairs = (hz_sw_pair *)malloc((size_t)py->nleaf * sizeof *pairs);
  order = (int32_t **)calloc(6, sizeof *order);
  counts = (int32_t *)calloc(6, sizeof *counts);
  if (!Ed || !Eprev || !cn || !vindex || !bits || !pairs || !order || !counts) {
    rc = 2;
    goto done;
  }

  /* порядок листов на направление. scalar §835: марш (DFS пирамиды).
   * col §836: ключ (колонка, позиция по оси со знаком) — сортировка. */
  for (d = 0; d < 6; d++) {
    int32_t li;
    order[d] = (int32_t *)malloc((size_t)py->nleaf * sizeof *order[d]);
    if (!order[d]) {
      rc = 2;
      goto done;
    }
    if (o->mode == 1) {
      int axis = d / 2, sign = (d % 2 == 0) ? 1 : -1;
      int64_t s = (axis == 0) ? 1 : (axis == 1 ? (int64_t)py->nx : (int64_t)py->nx * py->ny);
      int64_t len = (axis == 0) ? py->nx : (axis == 1 ? py->ny : py->nz);
      for (li = 0; li < py->nleaf; li++) {
        int64_t id = py->leaf_id[li];
        int64_t col = id / s, pos = id % s;
        pairs[li].key = col * ((int64_t)1 << HZ_SW_POSBITS) + (sign > 0 ? pos : (len - 1 - pos));
        pairs[li].li = li;
      }
      qsort(pairs, (size_t)py->nleaf, sizeof *pairs, sw_cmp_pair);
      for (li = 0; li < py->nleaf; li++)
        order[d][li] = pairs[li].li;
    } else {
      hz_pyr_march_stat ms;
      hz_pyr_march(py, dirs[d], 0, &ms, bits, vindex);
      for (li = 0; li < py->nleaf; li++)
        order[d][vindex[li] - 1] = li;
    }
    counts[d] = py->nleaf;
  }

  for (it = 0; it < o->iters; it++) {
    double emitted = 0, absorbed = 0, lost = 0, recycled = 0, e_sum = 0, area_sum = 0;
    int32_t p;
    memset(Ed, 0, (size_t)nt * sizeof *Ed);
    for (p = 0; p < nt; p++)
      Eprev[p] = py->pcs[p].e;
    for (d = 0; d < 6; d++) {
      const double *om = dirs[d];
      double L = 0.0; /* луч входит в домён тёмным */
      int64_t prevcol = -1;
      int axis = d / 2;
      int64_t s = (axis == 0) ? 1 : (axis == 1 ? (int64_t)py->nx : (int64_t)py->nx * py->ny);
      int32_t k;
      for (k = 0; k < counts[d]; k++) {
        const hz_pyr_leaf *lf;
        double cntot = 0.0, tau, T, Lsurf = 0.0, ab_flux;
        int32_t u;
        if (o->mode == 1) { /* новая колонка — новый тёмный вход */
          int64_t col = py->leaf_id[order[d][k]] / s;
          if (col != prevcol) {
            L = 0.0;
            prevcol = col;
          }
        }
        lf = &py->leaf[order[d][k]];
        if (o->noprop) L = 0.0;
        for (u = 0; u < lf->npcs; u++) {
          p = py->csr[lf->pcs_first + u];
          cn[p] = area[p] * fabs(om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                                 om[2] * nrm[3 * (int64_t)p + 2]);
          cntot += cn[p];
        }
        tau = cntot / csec;
        T = o->tau0 ? 1.0 : exp(-tau);
        if (tau > 0.0) {
          for (u = 0; u < lf->npcs; u++) {
            p = py->csr[lf->pcs_first + u];
            double rho = o->rho < 0 ? kd[p] : o->rho;
            double cnm = fabs(om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                              om[2] * nrm[3 * (int64_t)p + 2]);
            Lsurf += (o->le + rho * Eprev[p] / (2.0 * M_PI)) * cn[p];
            /* косинусная проекция радианса на кусок; Σ ΔE·area = absorbed */
            Ed[p] += HZ_SW_W * L * (1.0 - T) * cnm;
          }
          Lsurf /= cntot;
          ab_flux = HZ_SW_W * L * (1.0 - T) * csec;
          absorbed += ab_flux;
          emitted += HZ_SW_W * o->le * (1.0 - T) * csec;            /* внешний вход Le */
          recycled += HZ_SW_W * (Lsurf - o->le) * (1.0 - T) * csec; /* ρ·E_prev */
          L = L * T + Lsurf * (1.0 - T);
        }
        if (o->noprop) L = 0.0;
      }
      lost += L * csec; /* остаток луча за границей домена */
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
    if (e_hist) e_hist[it] = st->e_avg;
  }

done:
  free(Ed);
  free(Eprev);
  free(cn);
  free(vindex);
  free(bits);
  free(pairs);
  if (order) {
    for (d = 0; d < 6; d++)
      free((void *)order[d]);
    free(order);
  }
  free(counts);
  return rc;
}
