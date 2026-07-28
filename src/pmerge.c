/* Огрубление слиянием. Разбор и оговорки — в `pmerge.h`. */

#include "pmerge.h"
#include "cut/qef.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Проб на полигон при проверке НЕПЕРЕКРЫТИЯ. Число выведено, а не подобрано:
 * перекрытие, занимающее меньше 1/16 площади, слиянию не мешает (ошибка
 * площади ниже допуска на неё), а 16 проб по каждой оси такое перекрытие
 * обнаруживают заведомо. */
#define PM_OVSAMP 8

typedef struct {
  int32_t a, b;
  double err;
} pm_pair;

static int cmp_pair(const void *x, const void *y) {
  const pm_pair *p = x, *q = y;
  return (p->err < q->err) ? -1 : ((p->err > q->err) ? 1 : 0);
}

static int32_t uf_find(int32_t *p, int32_t x) {
  while (p[x] != x) {
    p[x] = p[p[x]];
    x = p[x];
  }
  return x;
}

/* Мировая коробка полигона по краю. */
static void pbox(const hz_polyset *ps, int32_t k, double lo[3], double hi[3]) {
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

/* Общая плоскость СПИСКА полигонов и МАКСИМАЛЬНОЕ отклонение их вершин от неё —
 * то есть `dmax`, а не среднеквадратичное (Г40/Г44).
 *
 * ПРОВЕРЯЕТСЯ ОБЪЕДИНЕНИЕ ГРУПП, А НЕ ПАРА, и это исправление настоящей ошибки.
 * Первая редакция проверяла пару и сливала транзитивно: `A~B` и `B~C` проходили
 * по отдельности, а у `A∪B∪C` отклонение не ограничивалось ничем. Замер поймал
 * это сразу — `dmax` дорастал до `44 δ` при допуске `1 δ`, — и поймал его
 * ИМЕННО тот пересчёт `dmax`, которого требует Ш7 («ложный ноль опаснее
 * UNKNOWN»). Метрика картинки к этому слепа: при уехавшей на два метра
 * плоскости она даже УЛУЧШАЛАСЬ, потому что сравнивает поле с точным светом в
 * той же уехавшей точке (узор К13/К40/К94). */
static double group_plane(const hz_polyset *ps, const int32_t *mem, int32_t nm, double n[3],
                          double *off) {
  double s = 0.0, ns[3] = {0, 0, 0}, org[3] = {0, 0, 0};
  for (int32_t i = 0; i < nm; i++) {
    const hz_poly *P = &ps->p[mem[i]];
    double w = P->area;
    s += w;
    for (int c = 0; c < 3; c++) {
      ns[c] += w * P->n[c];
      org[c] += w * P->org[c];
    }
  }
  if (!(s > 0.0)) return 1e300;
  double nn = sqrt(ns[0] * ns[0] + ns[1] * ns[1] + ns[2] * ns[2]);
  if (!(nn > 0.0)) return 1e300;
  for (int c = 0; c < 3; c++) {
    n[c] = ns[c] / nn;
    org[c] /= s;
  }
  *off = n[0] * org[0] + n[1] * org[1] + n[2] * org[2];

  double dmax = 0.0;
  for (int32_t i = 0; i < nm; i++) {
    const hz_poly *P = &ps->p[mem[i]];
    for (int32_t l = P->l0; l < P->l0 + P->nloop; l++)
      for (int32_t q = ps->loop[l]; q < ps->loop[l + 1]; q++) {
        double x[3];
        hz_poly_world(P, ps->bv[(size_t)q * 2], ps->bv[(size_t)q * 2 + 1], x);
        double d = fabs(x[0] * n[0] + x[1] * n[1] + x[2] * n[2] - *off);
        if (d > dmax) dmax = d;
      }
  }
  return dmax;
}

/* Ворота по QEF: эрмитовы образцы обоих краёв. Невязка есть
 * СРЕДНЕКВАДРАТИЧНОЕ, поэтому она годится ТОЛЬКО как дешёвый отсев — связывает
 * решение точный `dmax` выше (Г40/Г44). */
static double qef_gate(const hz_polyset *ps, int32_t a, int32_t b) {
  hz_qef q;
  hz_qef_zero(&q);
  double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
  for (int side = 0; side < 2; side++) {
    int32_t k = side ? b : a;
    const hz_poly *P = &ps->p[k];
    for (int32_t l = P->l0; l < P->l0 + P->nloop; l++)
      for (int32_t i = ps->loop[l]; i < ps->loop[l + 1]; i++) {
        double x[3];
        hz_poly_world(P, ps->bv[(size_t)i * 2], ps->bv[(size_t)i * 2 + 1], x);
        hz_qef_add_sample(&q, x, P->n);
        for (int c = 0; c < 3; c++) {
          if (x[c] < lo[c]) lo[c] = x[c];
          if (x[c] > hi[c]) hi[c] = x[c];
        }
      }
  }
  if (q.n < 3) return 1e300;
  double v[3], resid = 0.0;
  if (hz_qef_solve(&q, lo, hi, v, &resid) == HZ_QEF_EMPTY) return 1e300;
  return sqrt(fabs(resid) / (double)q.n);
}

/* ПЕРЕКРЫВАЮТСЯ ЛИ ПРОЕКЦИИ. Добавлено по замеру Ш3: критерий «все точки в
 * пределах δ от общей плоскости» разрешает слить две ОДИНАКОВО СМОТРЯЩИЕ
 * поверхности ближе δ, и на грубом δ это съедает до 10.6% площади сцены, чего
 * баланс энергии не ловит вовсе. */
static int overlaps(const hz_polyset *ps, int32_t a, int32_t b) {
  const hz_poly *A = &ps->p[a], *B = &ps->p[b];
  int inside = 0, tried = 0;
  for (int i = 0; i < PM_OVSAMP; i++)
    for (int j = 0; j < PM_OVSAMP; j++) {
      double u = B->uvlo[0] + ((double)i + 0.5) / PM_OVSAMP * (B->uvhi[0] - B->uvlo[0]);
      double v = B->uvlo[1] + ((double)j + 0.5) / PM_OVSAMP * (B->uvhi[1] - B->uvlo[1]);
      if (!hz_poly_inside(ps, B, u, v)) continue;
      tried++;
      double x[3], q[3];
      hz_poly_world(B, u, v, x);
      for (int c = 0; c < 3; c++)
        q[c] = x[c] - A->org[c];
      double ua = q[0] * A->eu[0] + q[1] * A->eu[1] + q[2] * A->eu[2];
      double va = q[0] * A->ev[0] + q[1] * A->ev[1] + q[2] * A->ev[2];
      if (ua < A->uvlo[0] || ua > A->uvhi[0] || va < A->uvlo[1] || va > A->uvhi[1]) continue;
      if (hz_poly_inside(ps, A, ua, va)) inside++;
    }
  /* Соседние по ребру полигоны дают единичные попадания на самой кромке;
   * перекрытием считается заметная ДОЛЯ, а не факт попадания. */
  return (tried > 0) && (inside * 16 > tried);
}

int hz_merge(hz_pseglist *so, const hz_objmesh *m, const hz_pseglist *si, const hz_polyset *ps,
             const double *E, const double *rho, const hz_mergecfg *cfg, hz_mergestat *st) {
  memset(so, 0, sizeof *so);
  memset(st, 0, sizeof *st);
  st->nseg_in = si->nseg;
  /* ТОЛЬКО ПОЛИГОНЫ ИЗ УЧАСТКОВ. В наборе за ними идут ИСТОЧНИКИ, добавленные
   * `hz_poly_add_quad`: у них нет треугольников в сетке, сливать их не с чем и
   * незачем, а группа без треугольников даёт пустой участок на выходе. */
  const int32_t np = (ps->np < si->nseg) ? ps->np : si->nseg;

  double *bb = malloc((size_t)np * 6 * sizeof *bb);
  int32_t *par = malloc((size_t)np * sizeof *par);
  if (bb == NULL || par == NULL) {
    free(bb);
    free(par);
    return 2;
  }
  double emean = 0.0, atot = 0.0;
  for (int32_t k = 0; k < np; k++) {
    pbox(ps, k, bb + (size_t)k * 6, bb + (size_t)k * 6 + 3);
    par[k] = k;
    atot += ps->p[k].area;
    if (E != NULL) emean += fabs(E[(size_t)k * 3]) * ps->p[k].area;
  }
  emean = (atot > 0.0) ? emean / atot : 1.0;
  if (!(emean > 0.0)) emean = 1.0;

  /* --- кандидаты: коробки ближе δ --- */
  int64_t pcap = 4096, pn = 0;
  pm_pair *pl = malloc((size_t)pcap * sizeof *pl);
  if (pl == NULL) {
    free(bb);
    free(par);
    return 2;
  }
  uint64_t rnd = 0x9E3779B97F4A7C15ULL;
  for (int32_t a = 0; a < np; a++)
    for (int32_t b = a + 1; b < np; b++) {
      const double *A = bb + (size_t)a * 6, *B = bb + (size_t)b * 6;
      int near = 1;
      for (int c = 0; c < 3 && near; c++) {
        if (A[c] - B[3 + c] > cfg->delta) near = 0;
        if (B[c] - A[3 + c] > cfg->delta) near = 0;
      }
      if (!near) continue;
      st->npair++;

      double err = 0.0;
      /* ГЕОМЕТРИЯ ЗДЕСЬ НЕ ПРОВЕРЯЕТСЯ — она проверяется при СЛИЯНИИ, по
       * объединению групп (см. `group_plane`). Здесь только парные и потому
       * законно парные члены: поле, материал, перекрытие. */
      if (cfg->random) {
        /* НЕГАТИВНЫЙ КОНТРОЛЬ: метрика СЛУЧАЙНАЯ. Ошибка обязана стать O(1);
         * если не стала — конъюнкция ничего не решает и мерили не её. */
        rnd = rnd * 6364136223846793005ULL + 1442695040888963407ULL;
        err = (double)(rnd >> 11) / 9007199254740992.0;
      } else {
        if (cfg->use_geom) {
          /* Дешёвые ворота по QEF (среднеквадратичное — только отсев). */
          double qg = qef_gate(ps, a, b);
          if (qg > cfg->delta) {
            st->nrej_geom++;
            continue;
          }
          err += qg / cfg->delta;
        }
        if (cfg->use_rad && E != NULL) {
          double d = fabs(E[(size_t)a * 3] - E[(size_t)b * 3]) / emean;
          if (!(d < cfg->ltol)) {
            st->nrej_rad++;
            continue;
          }
          err += d / cfg->ltol;
        }
        if (cfg->use_mtl && rho != NULL) {
          double d = fabs(rho[a] - rho[b]);
          if (!(d < cfg->rtol)) {
            st->nrej_mtl++;
            continue;
          }
          err += d / (cfg->rtol > 0.0 ? cfg->rtol : 1.0);
        }
        if (cfg->use_overlap && overlaps(ps, a, b)) {
          st->nrej_overlap++;
          continue;
        }
      }
      if (pn >= pcap) {
        int64_t nc = pcap * 2;
        pm_pair *q = realloc(pl, (size_t)nc * sizeof *q);
        if (q == NULL) {
          free(pl);
          free(bb);
          free(par);
          return 2;
        }
        pl = q;
        pcap = nc;
      }
      pl[pn].a = a;
      pl[pn].b = b;
      pl[pn].err = err;
      pn++;
    }

  /* --- слияние по возрастанию ошибки (Т2: приоритет — качество) ---
   * ГЕОМЕТРИЯ ПРОВЕРЯЕТСЯ ЗДЕСЬ, ПО ОБЪЕДИНЕНИЮ ГРУПП. Списки членов ведутся
   * односвязно: `head[корень] -> nxt`. Слияние принимается, только если у
   * ОБЪЕДИНЁННОЙ группы `dmax < δ`; иначе пара отбрасывается, а группы живут
   * дальше и могут слиться с другими. */
  qsort(pl, (size_t)pn, sizeof *pl, cmp_pair);
  int32_t *head = malloc((size_t)np * sizeof *head);
  int32_t *nxt = malloc((size_t)np * sizeof *nxt);
  int32_t *mem = malloc((size_t)np * sizeof *mem);
  if (head == NULL || nxt == NULL || mem == NULL) {
    free(head);
    free(nxt);
    free(mem);
    free(pl);
    free(bb);
    free(par);
    return 2;
  }
  for (int32_t k = 0; k < np; k++) {
    head[k] = k;
    nxt[k] = -1;
  }
  int32_t nleft = np;
  for (int64_t i = 0; i < pn; i++) {
    if (cfg->target > 0 && nleft <= cfg->target) break;
    int32_t ra = uf_find(par, pl[i].a), rb = uf_find(par, pl[i].b);
    if (ra == rb) continue;
    int32_t nm = 0;
    for (int32_t x = head[ra]; x >= 0 && nm < np; x = nxt[x])
      mem[nm++] = x;
    for (int32_t x = head[rb]; x >= 0 && nm < np; x = nxt[x])
      mem[nm++] = x;
    if (cfg->use_geom && !cfg->random) {
      double n[3], off;
      if (!(group_plane(ps, mem, nm, n, &off) < cfg->delta)) {
        st->nrej_geom++;
        continue;
      }
    }
    /* сцепить списки */
    int32_t tail = head[ra];
    while (nxt[tail] >= 0)
      tail = nxt[tail];
    nxt[tail] = head[rb];
    par[rb] = ra;
    head[ra] = head[ra];
    nleft--;
    st->nmerged++;
  }
  free(head);
  free(nxt);
  free(mem);
  free(pl);
  free(bb);

  /* --- новая разметка: треугольник наследует корень своего полигона --- */
  int32_t *rank = malloc((size_t)np * sizeof *rank);
  so->label = malloc((size_t)m->nt * sizeof *so->label);
  so->seg = malloc((size_t)np * sizeof *so->seg);
  if (rank == NULL || so->label == NULL || so->seg == NULL) {
    free(rank);
    free(par);
    hz_seg_free(so);
    return 2;
  }
  for (int32_t k = 0; k < np; k++)
    rank[k] = -1;
  int32_t nn = 0;
  for (int32_t k = 0; k < np; k++) {
    int32_t r = uf_find(par, k);
    if (rank[r] < 0) {
      rank[r] = nn;
      so->seg[nn] = si->seg[si->label[ps->tri[ps->p[r].t0]]];
      so->seg[nn].area = 0.0;
      so->seg[nn].ntri = 0;
      so->seg[nn].dmax = 0.0;
      nn++;
    }
  }
  /* ПЛОСКОСТЬ ГРУППЫ ПЕРЕСЧИТЫВАЕТСЯ, а не берётся у первого: слияние меняет
   * и нормаль, и смещение, а `dmax` после него ДРУГОЙ. Ложный ноль опаснее
   * UNKNOWN — тот останавливает потребителя, этот пропускает. */
  double *acc = calloc((size_t)nn * 8, sizeof *acc);
  if (acc == NULL) {
    free(rank);
    free(par);
    hz_seg_free(so);
    return 2;
  }
  for (int32_t k = 0; k < np; k++) {
    int32_t g = rank[uf_find(par, k)];
    const hz_poly *P = &ps->p[k];
    double w = P->area;
    for (int c = 0; c < 3; c++) {
      acc[(size_t)g * 8 + (size_t)c] += w * P->n[c];
      acc[(size_t)g * 8 + 3 + (size_t)c] += w * P->org[c];
    }
    acc[(size_t)g * 8 + 6] += w;
  }
  for (int32_t g = 0; g < nn; g++) {
    double w = acc[(size_t)g * 8 + 6];
    if (!(w > 0.0)) continue;
    double n[3], nl = 0.0;
    for (int c = 0; c < 3; c++) {
      n[c] = acc[(size_t)g * 8 + (size_t)c] / w;
      nl += n[c] * n[c];
    }
    nl = sqrt(nl);
    if (!(nl > 0.0)) continue;
    double org[3];
    for (int c = 0; c < 3; c++) {
      n[c] /= nl;
      org[c] = acc[(size_t)g * 8 + 3 + (size_t)c] / w;
      so->seg[g].n[c] = n[c];
    }
    so->seg[g].off = n[0] * org[0] + n[1] * org[1] + n[2] * org[2];
  }
  free(acc);

  /* Треугольник знает свой ИСХОДНЫЙ участок; полигон — тот же индекс, потому
   * что `hz_poly_build` нумерует полигоны участками один к одному. */
  for (int32_t t = 0; t < m->nt; t++) {
    int32_t k = si->label[t];
    int32_t g = rank[uf_find(par, k)];
    so->label[t] = g;
    so->seg[g].ntri++;
    so->seg[g].area += hz_obj_tri_area(m, t);
    double p[3][3];
    hz_obj_tri(m, t, p);
    for (int i = 0; i < 3; i++) {
      double d = fabs(p[i][0] * so->seg[g].n[0] + p[i][1] * so->seg[g].n[1] +
                      p[i][2] * so->seg[g].n[2] - so->seg[g].off);
      if (d > so->seg[g].dmax) so->seg[g].dmax = d;
    }
  }
  for (int32_t g = 0; g < nn; g++)
    if (so->seg[g].dmax > st->dmax_worst) st->dmax_worst = so->seg[g].dmax;
  so->nseg = nn;
  so->delta = si->delta;
  st->nseg_out = nn;
  free(rank);
  free(par);
  return 0;
}
