/* pvert.c — неизвестные на вершинах границы. Разбор — в `pvert.h`. */
#include "pvert.h"
#include "phcube.h"
#include "ptree.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

int hz_vset_build(hz_vset *V, const hz_polyset *ps) {
  memset(V, 0, sizeof *V);
  if (ps->nbv <= 0) return 1;
  /* Единый номер: сварной, если есть, иначе собственный. `bw` — номер в таблице
   * сварки полигонизатора, и его диапазон СВОЙ, не `nbv` (А252). */
  int32_t wmax = 0;
  for (int32_t b = 0; b < ps->nbv; b++)
    if (ps->bw != NULL && ps->bw[b] + 1 > wmax) wmax = ps->bw[b] + 1;
  if (wmax < 0) wmax = 0;
  int32_t nid = wmax + ps->nbv;
  int32_t *map = malloc((size_t)nid * sizeof *map);
  V->id = malloc((size_t)ps->nbv * sizeof *V->id);
  if (map == NULL || V->id == NULL) {
    free(map);
    free(V->id);
    memset(V, 0, sizeof *V);
    return 2;
  }
  for (int32_t q = 0; q < nid; q++)
    map[q] = -1;
  int32_t nv = 0;
  for (int32_t b = 0; b < ps->nbv; b++) {
    int32_t raw = (ps->bw != NULL && ps->bw[b] >= 0) ? ps->bw[b] : wmax + b;
    if (map[raw] < 0) map[raw] = nv++;
    V->id[b] = map[raw];
  }
  free(map);
  V->nv = nv;
  V->x = calloc(3 * (size_t)nv, sizeof *V->x);
  V->n = calloc(3 * (size_t)nv, sizeof *V->n);
  V->area = calloc((size_t)nv, sizeof *V->area);
  V->nown = calloc((size_t)nv, sizeof *V->nown);
  if (V->x == NULL || V->n == NULL || V->area == NULL || V->nown == NULL) {
    hz_vset_free(V);
    return 2;
  }
  for (int32_t k = 0; k < ps->np; k++) {
    const hz_poly *p = &ps->p[k];
    for (int32_t b = ps->loop[p->l0]; b < ps->loop[p->l0 + p->nloop]; b++) {
      int32_t v = V->id[b];
      double u = ps->bv[2 * b], w = ps->bv[2 * b + 1];
      /* Положение восстанавливается из рамы владельца; у всех владельцев одна
       * точка, потому и сварилось. Пишется каждым — значение одно. */
      for (int c = 0; c < 3; c++)
        V->x[3 * v + c] = p->org[c] + u * p->eu[c] + w * p->ev[c];
      /* Нормаль — средневзвешенная по площади: вершина принадлежит нескольким
       * граням, и своей нормали у неё нет. */
      for (int c = 0; c < 3; c++)
        V->n[3 * v + c] += p->n[c] * p->area;
      V->area[v] += p->area;
      V->nown[v]++;
    }
  }
  for (int32_t v = 0; v < nv; v++) {
    double l = 0.0;
    for (int c = 0; c < 3; c++)
      l += V->n[3 * v + c] * V->n[3 * v + c];
    l = sqrt(l);
    if (l > 0.0)
      for (int c = 0; c < 3; c++)
        V->n[3 * v + c] /= l;
    if (V->nown[v] == 1)
      V->n_one++;
    else if (V->nown[v] > 1)
      V->n_many++;
  }
  return 0;
}

void hz_vset_free(hz_vset *V) {
  free(V->id);
  free(V->x);
  free(V->n);
  free(V->area);
  free(V->nown);
  memset(V, 0, sizeof *V);
}

/* Предок полигона по тому же критерию, что в элементной сборке: пиксель голосует
 * за тот уровень, на котором излучатель уже не разрешается полукубом. */
static int32_t vt_pick(const hz_lod *L, int32_t poly, double r, double eps) {
  if (poly < 0 || poly >= L->np) return -1;
  int32_t nd = L->lab[poly];
  if (!(r > 0.0)) return nd;
  for (;;) {
    int32_t p = L->nd[nd].parent;
    if (p < 0 || p >= L->nnd) break;
    if (sqrt(L->nd[p].area_surf) / r >= eps) break;
    nd = p;
  }
  return nd;
}

int hz_vert_build(hz_linkset *S, const hz_vset *V, const hz_scene *sc, const hz_linkcfg *cfg,
                  int R) {
  memset(S, 0, sizeof *S);
  const hz_lod *L = sc->L;
  const hz_polyset *ps = sc->ps;
  if (ps == NULL || sc->m == NULL) return 1;
  int32_t *t2p = malloc((size_t)sc->m->nt * sizeof *t2p);
  if (t2p == NULL) return 2;
  for (int32_t k = 0; k < ps->np; k++)
    for (int32_t t = ps->p[k].t0; t < ps->p[k].t0 + ps->p[k].ntri; t++)
      t2p[ps->tri[t]] = k;
  /* ВЛАДЕЛЬЦЫ ВЕРШИНЫ — их полигоны в полукуб не рисуются: это правило §88
   * («элемент не затеняет сам себя») на вершине. Список хранится сжатым. */
  int32_t *ostart = calloc((size_t)V->nv + 1, sizeof *ostart);
  if (ostart == NULL) {
    free(t2p);
    return 2;
  }
  for (int32_t k = 0; k < ps->np; k++)
    for (int32_t b = ps->loop[ps->p[k].l0]; b < ps->loop[ps->p[k].l0 + ps->p[k].nloop]; b++)
      ostart[V->id[b] + 1]++;
  for (int32_t v = 0; v < V->nv; v++)
    ostart[v + 1] += ostart[v];
  int32_t *own = malloc((size_t)(ostart[V->nv] > 0 ? ostart[V->nv] : 1) * sizeof *own);
  int32_t *pos = malloc((size_t)V->nv * sizeof *pos);
  if (own == NULL || pos == NULL) {
    free(t2p);
    free(ostart);
    free(own);
    free(pos);
    return 2;
  }
  memcpy(pos, ostart, (size_t)V->nv * sizeof *pos);
  for (int32_t k = 0; k < ps->np; k++)
    for (int32_t b = ps->loop[ps->p[k].l0]; b < ps->loop[ps->p[k].l0 + ps->p[k].nloop]; b++)
      own[pos[V->id[b]]++] = k;
  free(pos);

  hz_ptree T;
  if (hz_ptree_build(&T, sc->m, cfg->ptleaf, 0) != 0) {
    free(t2p);
    free(ostart);
    free(own);
    return 2;
  }
  double eps = cfg->eps, epx = 2.0 / (double)R;
  if (eps < epx) eps = epx;
  int64_t npix = 0, nmiss = 0, nsetup = 0;
  int rc = 0;
#pragma omp parallel reduction(+ : npix, nmiss, nsetup)
  {
    hz_hcube h;
    double *acc = calloc((size_t)L->nnd, sizeof *acc);
    int32_t *touch = malloc((size_t)L->nnd * sizeof *touch);
    int32_t *stamp = calloc((size_t)sc->m->nt, sizeof *stamp);
    signed char *skip = calloc((size_t)ps->np, 1);
    hz_link *loc = NULL;
    int64_t nloc = 0, cloc = 0;
    hz_hcube_stat hst;
    memset(&hst, 0, sizeof hst);
    int ok = (acc != NULL && touch != NULL && stamp != NULL && skip != NULL &&
              hz_hcube_init(&h, R) == 0);
#pragma omp for schedule(dynamic, 64)
    for (int32_t v = 0; v < V->nv; v++) {
      if (!ok) continue;
      for (int32_t q = ostart[v]; q < ostart[v + 1]; q++)
        skip[own[q]] = 1;
      hz_hcube_draw_tree(&h, sc->m, t2p, &T, V->x + 3 * v, V->n + 3 * v, -1, stamp, v + 1,
                         cfg->nozb, &hst);
      int32_t nt2 = 0;
      for (int i = 0; i < h.npix; i++) {
        double w = h.dff[i];
        if (!(w > 0.0)) continue;
        npix++;
        int32_t pid = h.id[i];
        if (pid < 0 || skip[pid]) {
          nmiss++;
          continue;
        }
        int32_t e = vt_pick(L, pid, h.depth[i], eps);
        if (e < 0) continue;
        if (!(acc[e] > 0.0)) touch[nt2++] = e;
        acc[e] += w;
      }
      for (int32_t q = 0; q < nt2; q++) {
        int32_t e = touch[q];
        if (nloc >= cloc) {
          int64_t nc = (cloc > 0) ? cloc * 2 : 8192;
          hz_link *nl = realloc(loc, (size_t)nc * sizeof *nl);
          if (nl == NULL) {
            ok = 0;
            break;
          }
          loc = nl;
          cloc = nc;
        }
        loc[nloc].i = v;
        loc[nloc].j = e;
        loc[nloc].f = acc[e];
        nloc++;
        acc[e] = 0.0;
      }
      for (int32_t q = ostart[v]; q < ostart[v + 1]; q++)
        skip[own[q]] = 0;
    }
    nsetup += hst.nsetup;
#pragma omp critical
    {
      if (!ok) rc = 2;
      for (int64_t q = 0; q < nloc && rc == 0; q++) {
        if (S->n >= S->cap) {
          int64_t nc = (S->cap > 0) ? S->cap * 2 : 1 << 20;
          hz_link *nl = realloc(S->l, (size_t)nc * sizeof *nl);
          if (nl == NULL) {
            rc = 2;
            break;
          }
          S->l = nl;
          S->cap = nc;
        }
        S->l[S->n++] = loc[q];
      }
    }
    free(acc);
    free(touch);
    free(stamp);
    free(skip);
    free(loc);
    hz_hcube_free(&h);
  }
  S->nray = npix;
  S->nzero = nmiss;
  S->nrefine = nsetup;
  hz_ptree_free(&T);
  free(t2p);
  free(ostart);
  free(own);
  return rc;
}

int hz_vert_solve(const hz_linkset *S, const hz_vset *V, const hz_polyset *ps, const hz_lod *L,
                  const double *Le, const double *rho, double *Bv, double *Bp, double *Bn,
                  double tol, int maxit, double *reshist) {
  double *acc = malloc((size_t)V->nv * sizeof *acc);
  double *prev = malloc((size_t)V->nv * sizeof *prev);
  double *warea = calloc((size_t)L->nnd, sizeof *warea);
  if (acc == NULL || prev == NULL || warea == NULL) {
    free(acc);
    free(prev);
    free(warea);
    return -1;
  }
  for (int32_t v = 0; v < V->nv; v++)
    Bv[v] = Le[v];
  for (int32_t k = 0; k < ps->np; k++)
    warea[L->lab[k]] += ps->p[k].area;
  int it = 0;
  for (; it < maxit; it++) {
    memcpy(prev, Bv, (size_t)V->nv * sizeof *prev);
    /* ВЕРШИНЫ → ПОЛИГОНЫ: радианс полигона есть среднее его вершин. */
    for (int32_t k = 0; k < ps->np; k++) {
      const hz_poly *p = &ps->p[k];
      int32_t b0 = ps->loop[p->l0], b1 = ps->loop[p->l0 + p->nloop];
      double s = 0.0;
      int32_t n = 0;
      for (int32_t b = b0; b < b1; b++) {
        s += Bv[V->id[b]];
        n++;
      }
      Bp[k] = (n > 0) ? s / n : 0.0;
    }
    /* ПОЛИГОНЫ → УЗЛЫ: средневзвешенное по площади. */
    for (int32_t q = 0; q < L->nnd; q++)
      Bn[q] = 0.0;
    for (int32_t k = 0; k < ps->np; k++)
      Bn[L->lab[k]] += Bp[k] * ps->p[k].area;
    for (int32_t q = 0; q < L->nnd; q++)
      if (warea[q] > 0.0) Bn[q] /= warea[q];
    /* Грубые узлы: снизу вверх, средневзвешенное по площади детей. */
    for (int lev = 1; lev < L->nlev; lev++) {
      for (int32_t q = 0; q < L->nnd; q++)
        if (L->nd[q].level == lev) Bn[q] = 0.0;
      for (int32_t q = 0; q < L->nnd; q++) {
        if (L->nd[q].level != lev - 1) continue;
        int32_t p = L->nd[q].parent;
        if (p >= 0 && p < L->nnd) Bn[p] += Bn[q] * L->nd[q].area_surf;
      }
      for (int32_t q = 0; q < L->nnd; q++)
        if (L->nd[q].level == lev && L->nd[q].area_surf > 0.0) Bn[q] /= L->nd[q].area_surf;
    }
    for (int32_t v = 0; v < V->nv; v++)
      acc[v] = 0.0;
    for (int64_t l = 0; l < S->n; l++)
      acc[S->l[l].i] += S->l[l].f * Bn[S->l[l].j];
    double dmx = 0.0, bmx = 0.0;
    for (int32_t v = 0; v < V->nv; v++) {
      Bv[v] = Le[v] + rho[v] * acc[v];
      double d = fabs(Bv[v] - prev[v]), b = fabs(Bv[v]);
      if (d > dmx) dmx = d;
      if (b > bmx) bmx = b;
    }
    double res = (bmx > 0.0) ? dmx / bmx : 0.0;
    if (reshist != NULL) reshist[it] = res;
    if (res < tol) {
      it++;
      break;
    }
  }
  free(acc);
  free(prev);
  free(warea);
  return it;
}
