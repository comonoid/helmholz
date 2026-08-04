/* ptrace.c — реализация; разбор, зачем этот файл, — в `ptrace.h`. */
#include "ptrace.h"
#include "pgrid.h"
#include <math.h>
#include <stddef.h>

/* Состояние одного луча. Собрано в структуру, чтобы рекурсия несла один
 * указатель, а не двенадцать аргументов. */
typedef struct {
  const hz_ptree *t;
  const hz_objmesh *m;
  const double *o, *d;
  const int32_t *skip;
  double tmin, best;
  int64_t *nnode, *ntest;
  int32_t tri;
  int nskip, anyhit, found, mask;
} pt_ray;

/* Пересечение луча с коробкой. Взято дословно из `hz_pgrid_trace_tri` (отсечение
 * луча габаритом сетки), включая обработку нулевой компоненты направления через
 * `HZ_PGRID_DTINY`: у двух трассировщиков договор о том, что считается лучом
 * ВДОЛЬ грани, обязан быть один. */
static int pt_box(const double lo[3], const double hi[3], const double o[3], const double d[3],
                  double t0in, double t1in, double *tent) {
  double t0 = t0in, t1 = t1in;
  for (int a = 0; a < 3; a++) {
    if (fabs(d[a]) < HZ_PGRID_DTINY) {
      if (o[a] < lo[a] || o[a] > hi[a]) return 0;
    } else {
      double ta = (lo[a] - o[a]) / d[a], tb = (hi[a] - o[a]) / d[a];
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
  *tent = t0;
  return 1;
}

/* Спуск. Возврат 1 — ответ готов (только при `anyhit`), 0 — продолжать. */
static int pt_walk(pt_ray *r, int32_t ni, double tent) {
  if (tent >= r->best) return 0;
  const hz_ptnode *nd = &r->t->nd[ni];
  if (r->nnode != NULL) (*r->nnode)++;
  for (int32_t i = nd->t0; i < nd->t0 + nd->ntri; i++) {
    int32_t tr = r->t->ref[i];
    int sk = 0;
    for (int q = 0; q < r->nskip; q++)
      if (r->skip[q] == tr) sk = 1;
    if (sk) continue;
    if (r->ntest != NULL) (*r->ntest)++;
    double th;
    if (!hz_pgrid_tri_hit(r->m, tr, r->o, r->d, &th)) continue;
    if (th <= r->tmin || th >= r->best) continue;
    if (r->anyhit) {
      r->best = th;
      r->tri = tr;
      r->found = 1;
      return 1;
    }
    r->best = th;
    r->tri = tr;
    r->found = 1;
  }
  if (nd->child < 0) return 0;
  int32_t c0 = nd->child;
  for (int j = 0; j < 8; j++) {
    int32_t ch = c0 + (j ^ r->mask);
    const hz_ptnode *C = &r->t->nd[ch];
    /* ПУСТОТА ОТСЕИВАЕТСЯ ДО СПУСКА (§209.2). Пустой лист — это и есть «пустое
     * пространство одним узлом»: дробление останавливается при `n <= leafmax`,
     * поэтому лестницы из пустых узлов не бывает. */
    if (C->ntri == 0 && C->child < 0) continue;
    double te;
    if (!pt_box(C->lo, C->hi, r->o, r->d, r->tmin, r->best, &te)) continue;
    if (pt_walk(r, ch, te)) return 1;
  }
  return 0;
}

int hz_ptrace_cnt(const hz_ptree *t, const hz_objmesh *m, const double o[3], const double dir[3],
                  double tmin, double tmax, const int32_t *skip, int nskip, int anyhit,
                  double *thit, int32_t *tri, int64_t *nnode, int64_t *ntest) {
  if (t->nnd <= 0) return 0;
  double tent;
  if (!pt_box(t->nd[0].lo, t->nd[0].hi, o, dir, tmin, tmax, &tent)) return 0;
  pt_ray r;
  r.t = t;
  r.m = m;
  r.o = o;
  r.d = dir;
  r.skip = skip;
  r.tmin = tmin;
  r.best = tmax;
  r.nnode = nnode;
  r.ntest = ntest;
  r.tri = -1;
  r.nskip = nskip;
  r.anyhit = anyhit;
  r.found = 0;
  r.mask = 0;
  for (int c = 0; c < 3; c++)
    if (dir[c] < 0.0) r.mask |= (1 << c);
  pt_walk(&r, 0, tent);
  if (r.found) {
    *thit = r.best;
    if (tri != NULL) *tri = r.tri;
  }
  return r.found;
}
