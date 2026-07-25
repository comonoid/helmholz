#include "octree.h"
#include <complex.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/* upper bound for hz_oct_load: 2^26 nodes = 1.5 GiB of arena — beyond any
 * sane scene, and it bounds the attacker-controlled allocation (CWE-789) */
#define HZ_OCT_MAX_NODES (1 << 26)

static int32_t node_alloc8(hz_octree *t, double complex k2) {
  if (t->n + 8 > t->cap) {
    int32_t ncap = t->cap * 2;
    if (ncap < t->n + 8) ncap = t->n + 8;
    hz_onode *nn = realloc(t->nodes, (size_t)ncap * sizeof(hz_onode));
    if (nn == NULL) return -1;
    t->nodes = nn;
    t->cap = ncap;
  }
  int32_t base = t->n;
  for (int i = 0; i < 8; i++) {
    t->nodes[base + i].child0 = -1;
    t->nodes[base + i].k2 = k2;
  }
  t->n += 8;
  return base;
}

int hz_oct_init(hz_octree *t, int log2size, double complex background) {
  t->log2size = log2size;
  t->cap = 64;
  t->n = 0;
  t->nodes = calloc((size_t)t->cap, sizeof(hz_onode));
  if (t->nodes == NULL) return 1;
  t->nodes[0].child0 = -1;
  t->nodes[0].k2 = background;
  t->n = 1;
  return 0;
}

void hz_oct_free(hz_octree *t) {
  free(t->nodes);
  t->nodes = NULL;
  t->n = t->cap = 0;
}

/* --- set_box ------------------------------------------------------------- */

static int box_disjoint(const int alo[3], const int ahi[3], const int blo[3], const int bhi[3]) {
  for (int a = 0; a < 3; a++)
    if (ahi[a] <= blo[a] || bhi[a] <= alo[a]) return 1;
  return 0;
}

static int box_contains(const int outlo[3], const int outhi[3], const int inlo[3],
                        const int inhi[3]) {
  for (int a = 0; a < 3; a++)
    if (inlo[a] < outlo[a] || inhi[a] > outhi[a]) return 0;
  return 1;
}

static int set_box_rec(hz_octree *t, int32_t ni, const int nlo[3], int size, const int lo[3],
                       const int hi[3], double complex k2) {
  int nhi[3] = {nlo[0] + size, nlo[1] + size, nlo[2] + size};
  if (box_disjoint(nlo, nhi, lo, hi)) return 0;
  if (box_contains(lo, hi, nlo, nhi)) {
    t->nodes[ni].child0 = -1; /* pruned subtree stays in the arena (compaction TBD) */
    t->nodes[ni].k2 = k2;
    return 0;
  }
  if (t->nodes[ni].child0 < 0) {
    int32_t c = node_alloc8(t, t->nodes[ni].k2);
    if (c < 0) return 1;
    t->nodes[ni].child0 = c;
  }
  int half = size / 2;
  for (int i = 0; i < 8; i++) {
    int clo[3] = {nlo[0] + (i & 1 ? half : 0), nlo[1] + (i & 2 ? half : 0),
                  nlo[2] + (i & 4 ? half : 0)};
    if (set_box_rec(t, t->nodes[ni].child0 + i, clo, half, lo, hi, k2)) return 1;
  }
  return 0;
}

int hz_oct_set_box(hz_octree *t, const int lo[3], const int hi[3], double complex k2) {
  int size = 1 << t->log2size;
  int clo[3], chi[3];
  for (int a = 0; a < 3; a++) {
    clo[a] = lo[a] < 0 ? 0 : lo[a];
    chi[a] = hi[a] > size ? size : hi[a];
    if (clo[a] >= chi[a]) return 0;
  }
  int zlo[3] = {0, 0, 0};
  return set_box_rec(t, 0, zlo, size, clo, chi, k2);
}

/* --- set_ball ------------------------------------------------------------ */

static int ball_rec(hz_octree *t, int32_t ni, const int nlo[3], int size, const double c[3],
                    double r, double complex k2) {
  /* farthest / nearest corner distances to the ball center */
  double dmin2 = 0.0, dmax2 = 0.0;
  for (int a = 0; a < 3; a++) {
    double lo = (double)nlo[a], hi = (double)(nlo[a] + size);
    double dlo = c[a] - lo, dhi = hi - c[a];
    double far = dlo > dhi ? dlo : dhi;
    dmax2 += far * far;
    double near = 0.0;
    if (c[a] < lo)
      near = lo - c[a];
    else if (c[a] > hi)
      near = c[a] - hi;
    dmin2 += near * near;
  }
  if (dmin2 >= r * r) return 0; /* disjoint */
  if (dmax2 <= r * r) {         /* fully inside */
    t->nodes[ni].child0 = -1;
    t->nodes[ni].k2 = k2;
    return 0;
  }
  if (size == 1) { /* partial unit cell: center-inside test */
    double d2 = 0.0;
    for (int a = 0; a < 3; a++) {
      double d = (double)nlo[a] + 0.5 - c[a];
      d2 += d * d;
    }
    if (d2 < r * r) {
      t->nodes[ni].child0 = -1;
      t->nodes[ni].k2 = k2;
    }
    return 0;
  }
  if (t->nodes[ni].child0 < 0) {
    int32_t ch = node_alloc8(t, t->nodes[ni].k2);
    if (ch < 0) return 1;
    t->nodes[ni].child0 = ch;
  }
  int half = size / 2;
  for (int i = 0; i < 8; i++) {
    int clo[3] = {nlo[0] + (i & 1 ? half : 0), nlo[1] + (i & 2 ? half : 0),
                  nlo[2] + (i & 4 ? half : 0)};
    if (ball_rec(t, t->nodes[ni].child0 + i, clo, half, c, r, k2)) return 1;
  }
  return 0;
}

int hz_oct_set_ball(hz_octree *t, const double c[3], double r, double complex k2) {
  int zlo[3] = {0, 0, 0};
  return ball_rec(t, 0, zlo, 1 << t->log2size, c, r, k2);
}

/* --- queries ------------------------------------------------------------- */

double complex hz_oct_at(const hz_octree *t, int x, int y, int z) {
  int size = 1 << t->log2size;
  if (x < 0 || y < 0 || z < 0 || x >= size || y >= size || z >= size) return 0.0;
  int32_t ni = 0;
  int lo[3] = {0, 0, 0};
  while (t->nodes[ni].child0 >= 0) {
    size /= 2;
    int idx = 0;
    if (x >= lo[0] + size) {
      idx |= 1;
      lo[0] += size;
    }
    if (y >= lo[1] + size) {
      idx |= 2;
      lo[1] += size;
    }
    if (z >= lo[2] + size) {
      idx |= 4;
      lo[2] += size;
    }
    ni = t->nodes[ni].child0 + idx;
  }
  return t->nodes[ni].k2;
}

static int visit_rec(const hz_octree *t, int32_t ni, const int nlo[3], int size, const int lo[3],
                     const int hi[3], hz_oct_cb cb, void *ctx) {
  int nhi[3] = {nlo[0] + size, nlo[1] + size, nlo[2] + size};
  if (box_disjoint(nlo, nhi, lo, hi)) return 0;
  if (t->nodes[ni].child0 < 0) {
    int blo[3], bhi[3];
    for (int a = 0; a < 3; a++) {
      blo[a] = nlo[a] > lo[a] ? nlo[a] : lo[a];
      bhi[a] = nhi[a] < hi[a] ? nhi[a] : hi[a];
    }
    return cb(ctx, blo, bhi, t->nodes[ni].k2);
  }
  int half = size / 2;
  for (int i = 0; i < 8; i++) {
    int clo[3] = {nlo[0] + (i & 1 ? half : 0), nlo[1] + (i & 2 ? half : 0),
                  nlo[2] + (i & 4 ? half : 0)};
    int rc = visit_rec(t, t->nodes[ni].child0 + i, clo, half, lo, hi, cb, ctx);
    if (rc) return rc;
  }
  return 0;
}

int hz_oct_visit(const hz_octree *t, const int lo[3], const int hi[3], hz_oct_cb cb, void *ctx) {
  int zlo[3] = {0, 0, 0};
  return visit_rec(t, 0, zlo, 1 << t->log2size, lo, hi, cb, ctx);
}

/* --- serialization -------------------------------------------------------- */

int hz_oct_save(const hz_octree *t, FILE *f) {
  int32_t hdr[2] = {(int32_t)t->log2size, t->n};
  if (fwrite(hdr, sizeof(hdr), 1, f) != 1) return 1;
  if (fwrite(t->nodes, sizeof(hz_onode), (size_t)t->n, f) != (size_t)t->n) return 1;
  return 0;
}

int hz_oct_load(hz_octree *t, FILE *f) {
  int32_t hdr[2];
  if (fread(hdr, sizeof(hdr), 1, f) != 1) return 1;
  if (hdr[0] < 0 || hdr[0] > 20 || hdr[1] < 1 || hdr[1] > HZ_OCT_MAX_NODES) return 1;
  t->log2size = hdr[0];
  t->n = t->cap = hdr[1];
  t->nodes = calloc((size_t)t->cap, sizeof(hz_onode));
  if (t->nodes == NULL) return 1;
  if (fread(t->nodes, sizeof(hz_onode), (size_t)t->n, f) != (size_t)t->n) {
    hz_oct_free(t);
    return 1;
  }
  /* child indices must stay inside the arena */
  for (int32_t i = 0; i < t->n; i++)
    if (t->nodes[i].child0 >= 0 && t->nodes[i].child0 + 8 > t->n) {
      hz_oct_free(t);
      return 1;
    }
  return 0;
}
