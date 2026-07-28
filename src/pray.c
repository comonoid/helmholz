/* Луч по полигонам. Разбор и оговорки — в `pray.h`. */

#include "pray.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Мировая коробка полигона по его краю. */
static void poly_box(const hz_polyset *ps, int32_t k, double lo[3], double hi[3]) {
  const hz_poly *P = &ps->p[k];
  for (int a = 0; a < 3; a++) {
    lo[a] = 1e300;
    hi[a] = -1e300;
  }
  for (int32_t l = P->l0; l < P->l0 + P->nloop; l++)
    for (int32_t b = ps->loop[l]; b < ps->loop[l + 1]; b++) {
      double x[3];
      hz_poly_world(P, ps->bv[(size_t)b * 2], ps->bv[(size_t)b * 2 + 1], x);
      for (int a = 0; a < 3; a++) {
        if (x[a] < lo[a]) lo[a] = x[a];
        if (x[a] > hi[a]) hi[a] = x[a];
      }
    }
}

void hz_pray_free(hz_pray *g) {
  free(g->start);
  free(g->idx);
  memset(g, 0, sizeof *g);
}

int hz_pray_build(hz_pray *g, const hz_polyset *ps, double target) {
  memset(g, 0, sizeof *g);
  g->ps = ps;
  if (ps->np <= 0) return 1;
  for (int a = 0; a < 3; a++) {
    g->lo[a] = 1e300;
    g->hi[a] = -1e300;
  }
  double *bb = malloc((size_t)ps->np * 6 * sizeof *bb);
  if (bb == NULL) return 2;
  for (int32_t k = 0; k < ps->np; k++) {
    double lo[3], hi[3];
    poly_box(ps, k, lo, hi);
    for (int a = 0; a < 3; a++) {
      bb[(size_t)k * 6 + (size_t)a] = lo[a];
      bb[(size_t)k * 6 + 3 + (size_t)a] = hi[a];
      if (lo[a] < g->lo[a]) g->lo[a] = lo[a];
      if (hi[a] > g->hi[a]) g->hi[a] = hi[a];
    }
  }
  /* Число ячеек ~ np/target, распределённое по осям пропорционально габариту.
   * Порога здесь нет: и `target`, и габарит — данные, а не константы. */
  double ext[3], vol = 1.0;
  for (int a = 0; a < 3; a++) {
    ext[a] = g->hi[a] - g->lo[a];
    if (!(ext[a] > 0.0)) ext[a] = 1e-9;
    vol *= ext[a];
  }
  double want = (double)ps->np / (target > 0.0 ? target : 2.0);
  if (want < 1.0) want = 1.0;
  double s = cbrt(vol / want);
  int64_t tot = 1;
  for (int a = 0; a < 3; a++) {
    double n = floor(ext[a] / s) + 1.0;
    if (n > 512.0) n = 512.0;
    g->nc[a] = (int32_t)n;
    g->cs[a] = ext[a] / (double)g->nc[a];
    tot *= g->nc[a];
  }
  g->ncell = tot;
  g->start = calloc((size_t)tot + 1, sizeof *g->start);
  if (g->start == NULL) {
    free(bb);
    return 2;
  }
  /* два прохода: счёт, потом заполнение — растущих массивов нет */
  for (int pass = 0; pass < 2; pass++) {
    if (pass == 1) {
      /* Счётчики лежат в `start[c+1]`, поэтому префиксная сумма идёт ВПЕРЁД по
       * тому же массиву: `start[c+1] += start[c]`. Первая редакция считала
       * сумму со сдвигом на ячейку, курсоры заполнения писали в чужой пробег,
       * и `idx` наполнялся мусором — ASan поймал это как чтение по мёртвому
       * адресу в `poly_hit`, то есть за три вызова от места ошибки. */
      for (int64_t c = 0; c < tot; c++)
        g->start[c + 1] += g->start[c];
      int64_t acc = g->start[tot];
      g->nent = acc;
      g->idx = malloc((size_t)(acc > 0 ? acc : 1) * sizeof *g->idx);
      if (g->idx == NULL) {
        free(bb);
        return 2;
      }
    }
    for (int32_t k = 0; k < ps->np; k++) {
      int32_t c0[3], c1[3];
      for (int a = 0; a < 3; a++) {
        double t0 = (bb[(size_t)k * 6 + (size_t)a] - g->lo[a]) / g->cs[a];
        double t1 = (bb[(size_t)k * 6 + 3 + (size_t)a] - g->lo[a]) / g->cs[a];
        int32_t i0 = (int32_t)floor(t0), i1 = (int32_t)floor(t1);
        if (i0 < 0) i0 = 0;
        if (i1 > g->nc[a] - 1) i1 = g->nc[a] - 1;
        if (i1 < i0) i1 = i0;
        c0[a] = i0;
        c1[a] = i1;
      }
      for (int32_t z = c0[2]; z <= c1[2]; z++)
        for (int32_t y = c0[1]; y <= c1[1]; y++)
          for (int32_t x = c0[0]; x <= c1[0]; x++) {
            int64_t c = ((int64_t)z * g->nc[1] + y) * g->nc[0] + x;
            if (pass == 0)
              g->start[c + 1]++;
            else
              g->idx[g->start[c]++] = k;
          }
    }
  }
  /* курсоры сдвинулись — вернуть начала */
  for (int64_t c = tot; c > 0; c--)
    g->start[c] = g->start[c - 1];
  g->start[0] = 0;
  free(bb);
  return 0;
}

/* Пересечение луча с полигоном: плоскость, затем край. */
static int poly_hit(const hz_polyset *ps, int32_t k, const double o[3], const double d[3],
                    double *tout) {
  const hz_poly *P = &ps->p[k];
  if (P->nloop == 0) return 0;
  double dn = d[0] * P->n[0] + d[1] * P->n[1] + d[2] * P->n[2];
  if (!(fabs(dn) > 0.0)) return 0;
  double t = (P->off - (o[0] * P->n[0] + o[1] * P->n[1] + o[2] * P->n[2])) / dn;
  double q[3];
  for (int a = 0; a < 3; a++)
    q[a] = o[a] + t * d[a] - P->org[a];
  double u = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
  double v = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
  if (u < P->uvlo[0] || u > P->uvhi[0] || v < P->uvlo[1] || v > P->uvhi[1]) return 0;
  if (!hz_poly_inside(ps, P, u, v)) return 0;
  *tout = t;
  return 1;
}

/* Обход сеткой по алгоритму Amanatides–Woo; выход из коробки — конец. */
typedef struct {
  int32_t c[3], step[3];
  double tmax[3], tdelta[3], t;
  int alive;
} ray_walk;

static void walk_init(const hz_pray *g, const double o[3], const double d[3], double tmin,
                      ray_walk *w) {
  memset(w, 0, sizeof *w);
  /* вход в коробку сетки */
  double t0 = tmin, t1 = 1e300;
  for (int a = 0; a < 3; a++) {
    if (fabs(d[a]) > 0.0) {
      double ta = (g->lo[a] - o[a]) / d[a], tb = (g->hi[a] - o[a]) / d[a];
      if (ta > tb) {
        double s = ta;
        ta = tb;
        tb = s;
      }
      if (ta > t0) t0 = ta;
      if (tb < t1) t1 = tb;
    } else if (o[a] < g->lo[a] || o[a] > g->hi[a]) {
      return;
    }
  }
  if (t0 > t1) return;
  w->t = t0;
  for (int a = 0; a < 3; a++) {
    double p = o[a] + t0 * d[a];
    int32_t i = (int32_t)floor((p - g->lo[a]) / g->cs[a]);
    if (i < 0) i = 0;
    if (i > g->nc[a] - 1) i = g->nc[a] - 1;
    w->c[a] = i;
    if (d[a] > 0.0) {
      w->step[a] = 1;
      w->tmax[a] = t0 + (g->lo[a] + (double)(i + 1) * g->cs[a] - p) / d[a];
      w->tdelta[a] = g->cs[a] / d[a];
    } else if (d[a] < 0.0) {
      w->step[a] = -1;
      w->tmax[a] = t0 + (g->lo[a] + (double)i * g->cs[a] - p) / d[a];
      w->tdelta[a] = -g->cs[a] / d[a];
    } else {
      w->step[a] = 0;
      w->tmax[a] = 1e300;
      w->tdelta[a] = 1e300;
    }
  }
  w->alive = 1;
}

static int walk_next(const hz_pray *g, ray_walk *w) {
  int a = 0;
  if (w->tmax[1] < w->tmax[a]) a = 1;
  if (w->tmax[2] < w->tmax[a]) a = 2;
  w->t = w->tmax[a];
  w->c[a] += w->step[a];
  w->tmax[a] += w->tdelta[a];
  if (w->step[a] == 0 || w->c[a] < 0 || w->c[a] >= g->nc[a]) {
    w->alive = 0;
    return 0;
  }
  return 1;
}

int32_t hz_pray_hit(const hz_pray *g, const double o[3], const double d[3], double tmin,
                    double *tout) {
  ray_walk w;
  walk_init(g, o, d, tmin, &w);
  int32_t best = -1;
  double bt = 1e300;
  while (w.alive) {
    int64_t c = ((int64_t)w.c[2] * g->nc[1] + w.c[1]) * g->nc[0] + w.c[0];
    for (int32_t i = g->start[c]; i < g->start[c + 1]; i++) {
      double t;
      if (!poly_hit(g->ps, g->idx[i], o, d, &t)) continue;
      if (t > tmin && t < bt) {
        bt = t;
        best = g->idx[i];
      }
    }
    /* Найденное в этой ячейке может лежать ДАЛЬШЕ её границы — тогда впереди
     * может быть ближе. Останавливаемся только когда точка попадания заведомо
     * внутри пройденного. */
    double tnext = w.tmax[0];
    if (w.tmax[1] < tnext) tnext = w.tmax[1];
    if (w.tmax[2] < tnext) tnext = w.tmax[2];
    if (best >= 0 && bt <= tnext) break;
    if (!walk_next(g, &w)) break;
  }
  if (best >= 0) *tout = bt;
  return best;
}

int hz_pray_occluded(const hz_pray *g, const double o[3], const double d[3], double tmin,
                     double tmax) {
  ray_walk w;
  walk_init(g, o, d, tmin, &w);
  while (w.alive) {
    int64_t c = ((int64_t)w.c[2] * g->nc[1] + w.c[1]) * g->nc[0] + w.c[0];
    for (int32_t i = g->start[c]; i < g->start[c + 1]; i++) {
      double t;
      if (!poly_hit(g->ps, g->idx[i], o, d, &t)) continue;
      if (t > tmin && t < tmax) return 1;
    }
    if (w.t > tmax) break;
    if (!walk_next(g, &w)) break;
  }
  return 0;
}
