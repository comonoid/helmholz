/* pocc.c — занятость узлов дерева. Разбор — в `pocc.h`. */
#include "pocc.h"
#include "pclip.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
  const hz_objmesh *m;
  const hz_ptree *T;
  const double *tlo, *thi, *tpl;
  unsigned char *occ, *bocc;
  int fail;
} occctx;

static int box_hits(const double *bl, const double *bh, const double *clo, const double *chi) {
  for (int c = 0; c < 3; c++)
    if (!(bl[c] < chi[c]) || !(bh[c] >= clo[c])) return 0;
  return 1;
}

static void occ_walk(occctx *X, int32_t nid, const int32_t *list, int32_t n) {
  if (X->fail) return;
  const hz_ptnode *N = &X->T->nd[nid];
  if (n > 0 && X->bocc != NULL) X->bocc[nid] = 1;
  if (N->child < 0) {
    const double *lo = N->lo, *hi = N->hi;
    for (int32_t i = 0; i < n; i++) {
      int32_t t = list[i];
      const double *tp = X->tpl + 4 * (size_t)t;
      double ctr[3], hlf[3], sd = 0.0, rr = 0.0;
      for (int c = 0; c < 3; c++) {
        ctr[c] = 0.5 * (lo[c] + hi[c]);
        hlf[c] = 0.5 * (hi[c] - lo[c]);
        sd += tp[c] * ctr[c];
        rr += fabs(tp[c]) * hlf[c];
      }
      sd -= tp[3];
      /* Дешёвый отказ по плоскости: коробка целиком по одну сторону — куска нет
       * ТОЧНО, и это не приближение. */
      if (fabs(sd) > rr) continue;
      const double *A = X->m->v + 3 * (size_t)X->m->f[3 * (size_t)t + 0];
      const double *B = X->m->v + 3 * (size_t)X->m->f[3 * (size_t)t + 1];
      const double *C = X->m->v + 3 * (size_t)X->m->f[3 * (size_t)t + 2];
      hz_pclip_poly P;
      hz_pclip_tri(A, B, C, lo, hi, &P);
      if (P.nv >= 3 && hz_pclip_area(&P) > 0.0) {
        X->occ[nid] = 1;
        return;
      }
    }
    return;
  }
  int32_t c0 = N->child;
  int32_t *sub = malloc((size_t)(n > 0 ? n : 1) * sizeof *sub);
  if (sub == NULL) {
    X->fail = 1;
    return;
  }
  for (int k = 0; k < 8; k++) {
    const double *clo = X->T->nd[c0 + k].lo, *chi = X->T->nd[c0 + k].hi;
    int32_t ns = 0;
    for (int32_t i = 0; i < n; i++) {
      const double *bl = X->tlo + 3 * (size_t)list[i], *bh = X->thi + 3 * (size_t)list[i];
      if (box_hits(bl, bh, clo, chi)) sub[ns++] = list[i];
    }
    occ_walk(X, c0 + k, sub, ns);
  }
  free(sub);
}

int hz_pocc_build(const hz_objmesh *m, const hz_ptree *T, const double *tlo, const double *thi,
                  const double *tpl, unsigned char *occ, unsigned char *bocc) {
  int32_t *list = malloc((size_t)m->nt * sizeof *list);
  if (list == NULL) return 1;
  for (int32_t t = 0; t < m->nt; t++)
    list[t] = t;
  occctx X;
  X.m = m;
  X.T = T;
  X.tlo = tlo;
  X.thi = thi;
  X.tpl = tpl;
  X.occ = occ;
  X.bocc = bocc;
  X.fail = 0;
  occ_walk(&X, 0, list, m->nt);
  free(list);
  if (X.fail) return 1;
  /* Поднять вверх: если кусок есть у ребёнка, он есть и у родителя. Обратное
   * верно потому, что дети покрывают родителя целиком. Узлы дописываются после
   * родителя, поэтому обратный порядок индексов и есть порядок «снизу вверх». */
  for (int32_t i = T->nnd - 1; i >= 0; i--) {
    int32_t c0 = T->nd[i].child;
    if (c0 < 0) continue;
    for (int k = 0; k < 8; k++) {
      occ[i] = (unsigned char)(occ[i] | occ[c0 + k]);
      if (bocc != NULL) bocc[i] = (unsigned char)(bocc[i] | bocc[c0 + k]);
    }
  }
  return 0;
}
