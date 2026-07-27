/* PLAN_TRANSPORT.md, проход 3 со ЗЕРКАЛАМИ (К5 и К3). Разбор — в gather3.h. */

#include "transport/gather3.h"
#include <math.h>
#include <string.h>

/* Значение DG1-поля, хранимого в базисе ячейки `c`, в мировой точке `x`.
 * ЧИТАТЬ НУЛЕВОЙ КОЭФФИЦИЕНТ КАК РАДИАНС НЕЛЬЗЯ (К39): базис центрирован на
 * ЯЧЕЙКЕ, а элемент лежит на её краю, поэтому у постоянного поля `c`
 * представитель минимальной нормы равен (0.8c, 0.4c, 0, 0). */
static double dg1_at(const tr3_mesh *m, int32_t c, const double coef[4], const double x[3]) {
  double s = (double)m->csize[c], v = coef[0];
  for (int a = 0; a < 3; a++) {
    double xu = (x[a] - m->fr.o[a]) / m->fr.u[a];
    v += coef[a + 1] * ((xu - ((double)m->clo[c][a] + 0.5 * s)) / s);
  }
  return v;
}

/* Поверхностный элемент, в который попал луч: БЛИЖАЙШИЙ ПО НОРМАЛИ.
 * Искать по номеру фасета нельзя — в примитивном пути (К15) номер не
 * проставляется вовсе, потому что попадание считает примитив, и фасет к ответу
 * отношения не имеет. */
static int32_t pick_selem(const tr3_cut *cut, int32_t mc, const double n[3]) {
  double best = -2.0;
  int32_t bi = -1;
  for (int32_t k = cut->sestart[mc]; k < cut->sestart[mc + 1]; k++) {
    int32_t e = cut->selist[k];
    double dp = cut->se[e].n[0] * n[0] + cut->se[e].n[1] * n[1] + cut->se[e].n[2] * n[2];
    if (dp > best) {
      best = dp;
      bi = e;
    }
  }
  return bi;
}

/* Граничная грань, накрывающая точку выхода луча из куба. */
static int32_t wall_face(const tr3_gather *g, const double o[3], const double d[3], double hp[3]) {
  if (g->wallidx == NULL) return -1;
  double tex = 1e300;
  int wall = -1;
  for (int a = 0; a < 3; a++) {
    if (!(fabs(d[a]) > 0.0)) continue;
    double lim = d[a] > 0.0 ? (double)g->nwall : 0.0;
    double tt = (lim - o[a]) / d[a];
    if (tt > 0.0 && tt < tex) {
      tex = tt;
      wall = 2 * a + (d[a] > 0.0 ? 1 : 0);
    }
  }
  if (wall < 0) return -1;
  for (int a = 0; a < 3; a++)
    hp[a] = o[a] + tex * d[a];
  int axis = wall / 2, u = (axis + 1) % 3, v = (axis + 2) % 3;
  if (u > v) {
    int t = u;
    u = v;
    v = t;
  }
  int32_t iu = (int32_t)floor(hp[u]), iv = (int32_t)floor(hp[v]);
  if (iu < 0) iu = 0;
  if (iu >= g->nwall) iu = g->nwall - 1;
  if (iv < 0) iv = 0;
  if (iv >= g->nwall) iv = g->nwall - 1;
  return g->wallidx[((int32_t)wall * g->nwall + iu) * g->nwall + iv];
}

double tr3_gather_ray(const tr3_gather *g, const double o[3], const double d[3], int *nbounce) {
  double p[3], dir[3], thr = 1.0, val = 0.0;
  memcpy(p, o, sizeof p);
  memcpy(dir, d, sizeof dir);
  int nb = 0;
  for (;;) {
    tr3_hit h;
    tr3_march(g->sc, p, dir, -1.0, &h);
    if (!h.hit) { /* ушёл в стенку куба */
      double hp[3];
      int32_t f = wall_face(g, p, dir, hp);
      if (f >= 0) val += thr * dg1_at(g->m, g->m->f[f].ca, g->bout + (size_t)f * 4, hp);
      break;
    }
    int32_t mc = g->m->cellof[h.cell];
    if (mc < 0) break;
    int32_t e = pick_selem(g->cut, mc, h.n);
    if (e < 0) break;
    double spec = 0.0;
    if (g->facet_spec != NULL && g->cut->se[e].facet < g->nfacet)
      spec = g->facet_spec[g->cut->se[e].facet];
    if (spec > 0.0 && nb < g->maxbounce) {
      /* ЗАКОН ОТРАЖЕНИЯ по нормали ПРИМИТИВА: d' = d − 2(d·n)n.
       * Точка старта берётся РОВНО на поверхности, без смещения по нормали:
       * самопопадание исключено строгостью неравенств в марше (ray3.h). */
      double dn = dir[0] * h.n[0] + dir[1] * h.n[1] + dir[2] * h.n[2];
      for (int a = 0; a < 3; a++) {
        p[a] = h.p[a];
        dir[a] -= 2.0 * dn * h.n[a];
      }
      thr *= spec;
      nb++;
      continue;
    }
    val += thr * dg1_at(g->m, mc, g->sout + (size_t)e * 4, h.p);
    break;
  }
  if (nbounce != NULL) *nbounce = nb;
  return val > 0.0 ? val : 0.0;
}
