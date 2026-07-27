/* PLAN_TRANSPORT.md, шаг B. Сетка граней над октодеревом: грань есть
 * ПЕРЕСЕЧЕНИЕ граней соседей, висячего узла как случая нет. */

#include "transport/mesh3.h"
#include <stdlib.h>
#include <string.h>

void tr3_mesh_free(tr3_mesh *m) {
  free(m->node);
  free(m->cellof);
  free(m->clo);
  free(m->csize);
  free(m->f);
  free(m->fstart);
  free(m->flist);
  memset(m, 0, sizeof *m);
}

/* --- перечисление листьев --------------------------------------------------- */

static void leaves_rec(tr3_mesh *m, int32_t ni, const int32_t lo[3], int32_t size) {
  if (m->tree->nodes[ni].child0 >= 0) {
    int32_t half = size / 2;
    for (int i = 0; i < 8; i++) {
      int32_t clo[3] = {lo[0] + ((i & 1) ? half : 0), lo[1] + ((i & 2) ? half : 0),
                        lo[2] + ((i & 4) ? half : 0)};
      leaves_rec(m, m->tree->nodes[ni].child0 + i, clo, half);
    }
    return;
  }
  int32_t c = m->ncell++;
  m->node[c] = ni;
  m->cellof[ni] = c;
  for (int a = 0; a < 3; a++)
    m->clo[c][a] = lo[a];
  m->csize[c] = size;
}

/* Листья, пересекающие коробку [blo, bhi). Собираются в out. */
static void collect_rec(const hz_octree *t, int32_t ni, const int32_t nlo[3], int32_t size,
                        const int32_t blo[3], const int32_t bhi[3], int32_t *out, int32_t *n,
                        int32_t max) {
  for (int a = 0; a < 3; a++)
    if (nlo[a] + size <= blo[a] || bhi[a] <= nlo[a]) return;
  if (t->nodes[ni].child0 < 0) {
    if (*n < max) out[(*n)++] = ni;
    return;
  }
  int32_t half = size / 2;
  for (int i = 0; i < 8; i++) {
    int32_t clo[3] = {nlo[0] + ((i & 1) ? half : 0), nlo[1] + ((i & 2) ? half : 0),
                      nlo[2] + ((i & 4) ? half : 0)};
    collect_rec(t, t->nodes[ni].child0 + i, clo, half, blo, bhi, out, n, max);
  }
}

/* --- построение ------------------------------------------------------------- */

static int push_face(tr3_mesh *m, int32_t ca, int32_t cb, int axis, int32_t pos,
                     const int32_t lo[2], const int32_t hi[2]) {
  if (m->nf >= m->fcap) {
    int32_t nc = m->fcap * 2;
    tr3_face *nf = realloc(m->f, (size_t)nc * sizeof(tr3_face));
    if (nf == NULL) return 1;
    m->f = nf;
    m->fcap = nc;
  }
  tr3_face *f = &m->f[m->nf++];
  f->ca = ca;
  f->cb = cb;
  f->axis = axis;
  f->pos = pos;
  for (int k = 0; k < 2; k++) {
    f->lo[k] = lo[k];
    f->hi[k] = hi[k];
  }
  return 0;
}

int tr3_mesh_build(tr3_mesh *m, const hz_octree *t, const hz_frame *fr) {
  memset(m, 0, sizeof *m);
  m->tree = t;
  m->fr = *fr;
  int32_t n = (int32_t)1 << t->log2size;

  m->node = calloc((size_t)t->n, sizeof(int32_t));
  m->cellof = calloc((size_t)t->n, sizeof(int32_t));
  m->clo = calloc((size_t)t->n, sizeof(int32_t[3]));
  m->csize = calloc((size_t)t->n, sizeof(int32_t));
  if (m->node == NULL || m->cellof == NULL || m->clo == NULL || m->csize == NULL) {
    tr3_mesh_free(m);
    return 1;
  }
  for (int32_t i = 0; i < t->n; i++)
    m->cellof[i] = -1;
  int32_t zero[3] = {0, 0, 0};
  leaves_rec(m, 0, zero, n);

  m->fcap = 64;
  m->f = calloc((size_t)m->fcap, sizeof(tr3_face));
  if (m->f == NULL) {
    tr3_mesh_free(m);
    return 1;
  }
  int32_t *buf = calloc((size_t)t->n, sizeof(int32_t));
  if (buf == NULL) {
    tr3_mesh_free(m);
    return 1;
  }

  for (int32_t c = 0; c < m->ncell; c++) {
    int32_t s = m->csize[c];
    for (int a = 0; a < 3; a++) {
      int u = (a + 1) % 3, v = (a + 2) % 3;
      if (u > v) {
        int tmp = u;
        u = v;
        v = tmp;
      }
      int32_t rlo[2] = {m->clo[c][u], m->clo[c][v]};
      int32_t rhi[2] = {rlo[0] + s, rlo[1] + s};

      /* МИНУС-сторона: только внешняя граница. Внутренние грани заводит сосед,
       * стоящий с МЕНЬШЕЙ стороны, — так каждая внутренняя грань появляется
       * РОВНО ОДИН раз, и это свойство построения, а не проверка потом. */
      if (m->clo[c][a] == 0) {
        if (push_face(m, c, ~(int32_t)(2 * a), a, 0, rlo, rhi)) {
          free(buf);
          tr3_mesh_free(m);
          return 1;
        }
      }

      int32_t pos = m->clo[c][a] + s;
      if (pos == n) { /* внешняя граница с ПЛЮС стороны */
        if (push_face(m, c, ~(int32_t)(2 * a + 1), a, pos, rlo, rhi)) {
          free(buf);
          tr3_mesh_free(m);
          return 1;
        }
        continue;
      }
      /* соседи: листья в тонком слое толщиной 1 за плоскостью */
      int32_t blo[3], bhi[3];
      for (int k = 0; k < 3; k++) {
        blo[k] = m->clo[c][k];
        bhi[k] = m->clo[c][k] + s;
      }
      blo[a] = pos;
      bhi[a] = pos + 1;
      int32_t nb = 0;
      collect_rec(t, 0, zero, n, blo, bhi, buf, &nb, t->n);
      for (int32_t i = 0; i < nb; i++) {
        int32_t cb = m->cellof[buf[i]];
        if (cb < 0) continue;
        /* ПЕРЕСЕЧЕНИЕ прямоугольников: если сосед крупнее, выйдет грань этой
         * ячейки; если мельче — его. Обе стороны получат ОДИН И ТОТ ЖЕ. */
        int32_t flo[2], fhi[2];
        int ok = 1;
        int ax[2] = {u, v};
        for (int k = 0; k < 2; k++) {
          int32_t alo = m->clo[c][ax[k]], ahi = alo + s;
          int32_t nlo2 = m->clo[cb][ax[k]], nhi2 = nlo2 + m->csize[cb];
          flo[k] = alo > nlo2 ? alo : nlo2;
          fhi[k] = ahi < nhi2 ? ahi : nhi2;
          if (fhi[k] <= flo[k]) ok = 0;
        }
        if (!ok) continue;
        if (push_face(m, c, cb, a, pos, flo, fhi)) {
          free(buf);
          tr3_mesh_free(m);
          return 1;
        }
      }
    }
  }
  free(buf);

  /* CSR: грани, инцидентные ячейке */
  m->fstart = calloc((size_t)m->ncell + 1, sizeof(int32_t));
  if (m->fstart == NULL) {
    tr3_mesh_free(m);
    return 1;
  }
  for (int32_t i = 0; i < m->nf; i++) {
    m->fstart[m->f[i].ca + 1]++;
    if (m->f[i].cb >= 0) m->fstart[m->f[i].cb + 1]++;
  }
  for (int32_t c = 0; c < m->ncell; c++)
    m->fstart[c + 1] += m->fstart[c];
  /* пустая сетка возможна (дерево без листьев не бывает, но грани могут не
   * появиться при вырожденных параметрах), а calloc(0) — неопределённость */
  m->flist = calloc((size_t)(m->fstart[m->ncell] > 0 ? m->fstart[m->ncell] : 1), sizeof(int32_t));
  if (m->flist == NULL) {
    tr3_mesh_free(m);
    return 1;
  }
  int32_t *fill = calloc((size_t)m->ncell, sizeof(int32_t));
  if (fill == NULL) {
    tr3_mesh_free(m);
    return 1;
  }
  for (int32_t i = 0; i < m->nf; i++) {
    int32_t a = m->f[i].ca;
    m->flist[m->fstart[a] + fill[a]++] = i;
    if (m->f[i].cb >= 0) {
      int32_t b = m->f[i].cb;
      m->flist[m->fstart[b] + fill[b]++] = i;
    }
  }
  free(fill);
  return 0;
}

void tr3_face_corners(const tr3_mesh *m, int32_t fi, double v[4][3]) {
  const tr3_face *f = &m->f[fi];
  int a = f->axis, u = (a + 1) % 3, w = (a + 2) % 3;
  if (u > w) {
    int t = u;
    u = w;
    w = t;
  }
  const int32_t cu[4] = {f->lo[0], f->hi[0], f->hi[0], f->lo[0]};
  const int32_t cw[4] = {f->lo[1], f->lo[1], f->hi[1], f->hi[1]};
  for (int k = 0; k < 4; k++) {
    double p[3];
    p[a] = (double)f->pos;
    p[u] = (double)cu[k];
    p[w] = (double)cw[k];
    for (int c = 0; c < 3; c++)
      v[k][c] = m->fr.o[c] + m->fr.u[c] * p[c];
  }
}

double tr3_face_area(const tr3_mesh *m, int32_t fi) {
  const tr3_face *f = &m->f[fi];
  int a = f->axis, u = (a + 1) % 3, w = (a + 2) % 3;
  if (u > w) {
    int t = u;
    u = w;
    w = t;
  }
  return (double)(f->hi[0] - f->lo[0]) * m->fr.u[u] * (double)(f->hi[1] - f->lo[1]) * m->fr.u[w];
}
