/* pgrid.c — реализация; довод, зачем этот файл отдельно, — в `pgrid.h`.
 *
 * ПЕРЕНОС ИЗ `tools/psil.c` СДЕЛАН БЕЗ ЕДИНОЙ ПРАВКИ АРИФМЕТИКИ, и это
 * проверяется числом: `./build/psil city noocc n=32` обязан дать ту же медиану
 * `20518` (§172, пункт 0). Переименованы только имена; порядок операций,
 * константы и ветвления оставлены дословно.
 */
#include "pgrid.h"
#include <math.h>
#include <stdlib.h>

static void pgrid_tri_box(const hz_objmesh *m, int32_t t, double blo[3], double bhi[3]) {
  for (int a = 0; a < 3; a++) {
    blo[a] = 1e300;
    bhi[a] = -1e300;
  }
  for (int i = 0; i < 3; i++) {
    const double *p = m->v + 3 * (size_t)m->f[(size_t)t * 3 + (size_t)i];
    for (int a = 0; a < 3; a++) {
      if (p[a] < blo[a]) blo[a] = p[a];
      if (p[a] > bhi[a]) bhi[a] = p[a];
    }
  }
}

int hz_pgrid_build(hz_pgrid *g, const hz_objmesh *m) {
  g->start = NULL;
  g->idx = NULL;
  g->ncell = 0;
  g->nref = 0;
  double asum = 0.0;
  for (int32_t t = 0; t < m->nt; t++)
    asum += hz_obj_tri_area(m, t);
  double cs = HZ_PGRID_CELL_MUL * sqrt(2.0 * asum / (double)m->nt);
  for (int a = 0; a < 3; a++) {
    /* Коробка расширяется на ячейку: точка ровно на грани сцены обязана иметь
     * ячейку, а не выпадать из сетки. */
    g->lo[a] = m->lo[a] - cs;
    g->hi[a] = m->hi[a] + cs;
  }
  for (;;) {
    int64_t n = 1;
    for (int a = 0; a < 3; a++) {
      double w = (g->hi[a] - g->lo[a]) / cs;
      int32_t k = (int32_t)floor(w) + 1;
      g->nc[a] = k < 1 ? 1 : k;
      n *= g->nc[a];
    }
    if (n <= HZ_PGRID_CELL_MAX) {
      g->ncell = n;
      break;
    }
    cs *= 2.0;
  }
  g->cs = cs;
  g->start = calloc((size_t)g->ncell + 2, sizeof *g->start);
  if (g->start == NULL) return 2;
  for (int pass = 0; pass < 2; pass++) {
    if (pass == 1) {
      for (int64_t i = 0; i < g->ncell; i++)
        g->start[i + 1] += g->start[i];
      g->nref = g->start[g->ncell];
      g->idx = malloc((size_t)(g->nref > 0 ? g->nref : 1) * sizeof *g->idx);
      if (g->idx == NULL) return 2;
      for (int64_t i = g->ncell; i > 0; i--)
        g->start[i] = g->start[i - 1];
      g->start[0] = 0;
    }
    for (int32_t t = 0; t < m->nt; t++) {
      double blo[3], bhi[3];
      pgrid_tri_box(m, t, blo, bhi);
      /* Инициализация нулём не «на всякий случай»: без неё gcc-analyzer не
       * связывает цикл по трём осям с последующим использованием и даёт
       * ложное «use of uninitialized value». */
      int32_t c0[3] = {0, 0, 0}, c1[3] = {0, 0, 0};
      for (int a = 0; a < 3; a++) {
        int32_t i0 = (int32_t)floor((blo[a] - g->lo[a]) / cs);
        int32_t i1 = (int32_t)floor((bhi[a] - g->lo[a]) / cs);
        c0[a] = i0 < 0 ? 0 : (i0 >= g->nc[a] ? g->nc[a] - 1 : i0);
        c1[a] = i1 < 0 ? 0 : (i1 >= g->nc[a] ? g->nc[a] - 1 : i1);
      }
      for (int32_t z = c0[2]; z <= c1[2]; z++)
        for (int32_t y = c0[1]; y <= c1[1]; y++)
          for (int32_t x = c0[0]; x <= c1[0]; x++) {
            int64_t c = (int64_t)x + (int64_t)g->nc[0] * ((int64_t)y + (int64_t)g->nc[1] * z);
            if (pass == 0)
              g->start[c + 1]++;
            else
              g->idx[g->start[c + 1]++] = t;
          }
    }
  }
  return 0;
}

void hz_pgrid_free(hz_pgrid *g) {
  free(g->start);
  free(g->idx);
  g->start = NULL;
  g->idx = NULL;
}

int hz_pgrid_tri_hit(const hz_objmesh *m, int32_t t, const double o[3], const double d[3],
                     double *tout) {
  const double *p0 = m->v + 3 * (size_t)m->f[(size_t)t * 3 + 0];
  const double *p1 = m->v + 3 * (size_t)m->f[(size_t)t * 3 + 1];
  const double *p2 = m->v + 3 * (size_t)m->f[(size_t)t * 3 + 2];
  double e1[3], e2[3], pv[3], tv[3], qv[3];
  for (int a = 0; a < 3; a++) {
    e1[a] = p1[a] - p0[a];
    e2[a] = p2[a] - p0[a];
  }
  pv[0] = d[1] * e2[2] - d[2] * e2[1];
  pv[1] = d[2] * e2[0] - d[0] * e2[2];
  pv[2] = d[0] * e2[1] - d[1] * e2[0];
  double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
  if (fabs(det) < HZ_PGRID_DET_TINY) return 0; /* луч лежит в плоскости */
  double inv = 1.0 / det;
  for (int a = 0; a < 3; a++)
    tv[a] = o[a] - p0[a];
  double u = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
  if (u < 0.0 || u > 1.0) return 0;
  qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
  qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
  qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
  double v = (d[0] * qv[0] + d[1] * qv[1] + d[2] * qv[2]) * inv;
  if (v < 0.0 || u + v > 1.0) return 0;
  *tout = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
  return 1;
}

int hz_pgrid_trace(const hz_pgrid *g, const hz_objmesh *m, const double o[3], const double dir[3],
                   double tmin, double tmax, const int32_t *skip, int nskip, int anyhit,
                   double *thit) {
  return hz_pgrid_trace_tri(g, m, o, dir, tmin, tmax, skip, nskip, anyhit, thit, NULL);
}

int hz_pgrid_trace_tri(const hz_pgrid *g, const hz_objmesh *m, const double o[3],
                       const double dir[3], double tmin, double tmax, const int32_t *skip,
                       int nskip, int anyhit, double *thit, int32_t *tri) {
  return hz_pgrid_trace_cnt(g, m, o, dir, tmin, tmax, skip, nskip, anyhit, thit, tri, NULL, NULL,
                            NULL, -1);
}

int hz_pgrid_trace_cnt(const hz_pgrid *g, const hz_objmesh *m, const double o[3],
                       const double dir[3], double tmin, double tmax, const int32_t *skip,
                       int nskip, int anyhit, double *thit, int32_t *tri, int64_t *ncell,
                       int64_t *ntest, const int32_t *tag, int32_t tagskip) {
  double t0 = tmin, t1 = tmax;
  for (int a = 0; a < 3; a++) {
    if (fabs(dir[a]) < HZ_PGRID_DTINY) {
      if (o[a] < g->lo[a] || o[a] > g->hi[a]) return 0;
    } else {
      double ta = (g->lo[a] - o[a]) / dir[a], tb = (g->hi[a] - o[a]) / dir[a];
      if (ta > tb) {
        double s = ta;
        ta = tb;
        tb = s;
      }
      if (ta > t0) t0 = ta;
      if (tb < t1) t1 = tb;
    }
  }
  if (t0 > t1) return 0;
  int32_t ix[3], st[3];
  double tnext[3], tdel[3];
  for (int a = 0; a < 3; a++) {
    double p = o[a] + t0 * dir[a];
    int32_t k = (int32_t)floor((p - g->lo[a]) / g->cs);
    ix[a] = k < 0 ? 0 : (k >= g->nc[a] ? g->nc[a] - 1 : k);
    if (dir[a] > HZ_PGRID_DTINY) {
      st[a] = 1;
      tdel[a] = g->cs / dir[a];
      tnext[a] = (g->lo[a] + (double)(ix[a] + 1) * g->cs - o[a]) / dir[a];
    } else if (dir[a] < -HZ_PGRID_DTINY) {
      st[a] = -1;
      tdel[a] = -g->cs / dir[a];
      tnext[a] = (g->lo[a] + (double)ix[a] * g->cs - o[a]) / dir[a];
    } else {
      st[a] = 0;
      tdel[a] = 1e300;
      tnext[a] = 1e300;
    }
  }
  double best = tmax;
  int found = 0;
  for (;;) {
    int64_t c = (int64_t)ix[0] + (int64_t)g->nc[0] * ((int64_t)ix[1] + (int64_t)g->nc[1] * ix[2]);
    if (ncell != NULL) (*ncell)++;
    for (int64_t i = g->start[c]; i < g->start[c + 1]; i++) {
      int32_t t = g->idx[i];
      int sk = 0;
      for (int q = 0; q < nskip; q++)
        if (skip[q] == t) sk = 1;
      if (sk) continue;
      /* Треугольник с тегом `tagskip` заслоном не считается — маска самозаслона
       * (§237). Проверка стоит ДО пересечения: она дешевле его. */
      if (tag != NULL && tag[t] == tagskip) continue;
      if (ntest != NULL) (*ntest)++;
      double th;
      if (!hz_pgrid_tri_hit(m, t, o, dir, &th)) continue;
      if (th <= tmin || th >= best) continue;
      if (anyhit) {
        *thit = th;
        if (tri != NULL) *tri = t;
        return 1;
      }
      best = th;
      found = 1;
      if (tri != NULL) *tri = t;
    }
    int a =
        (tnext[0] < tnext[1]) ? ((tnext[0] < tnext[2]) ? 0 : 2) : ((tnext[1] < tnext[2]) ? 1 : 2);
    if (tnext[a] > t1 || (found && tnext[a] > best)) break;
    ix[a] += st[a];
    if (ix[a] < 0 || ix[a] >= g->nc[a]) break;
    tnext[a] += tdel[a];
  }
  if (found) *thit = best;
  return found;
}
