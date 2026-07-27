#include "octree.h"
#include <complex.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

int32_t hz_oct_leaf(const hz_octree *t, int x, int y, int z) {
  int size = 1 << t->log2size;
  if (x < 0 || y < 0 || z < 0 || x >= size || y >= size || z >= size) return -1;
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
  return ni;
}

double complex hz_oct_at(const hz_octree *t, int x, int y, int z) {
  int32_t ni = hz_oct_leaf(t, x, y, z);
  /* CMPLX and not the bare 0.0: CBMC models double complex as a struct and its
   * symex asserts type equality on assignment, so an implicit real->complex
   * conversion aborts the checker outright (tests/cbmc_octree.c) */
  return ni < 0 ? CMPLX(0.0, 0.0) : t->nodes[ni].k2;
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

/* --- validation (pure, integer; PLAN_CUT.md Р-2) -------------------------- */

/* depth scratch marker; real depths are <= HZ_OCT_MAX_LOG2SIZE = 20 */
#define HZ_OCT_UNSEEN 0xFFu

int hz_oct_validate(const hz_octree *t) {
  if (t->log2size < 0 || t->log2size > HZ_OCT_MAX_LOG2SIZE) return HZ_OCT_E_LOG2SIZE;
  /* (n-1) % 8: node 0 is the root and node_alloc8 hands out blocks of 8 after
   * it, so any other count is a file that this code could not have written. */
  if (t->n < 1 || t->n > HZ_OCT_MAX_NODES || (t->n - 1) % 8 != 0) return HZ_OCT_E_NCOUNT;

  /* inline malloc, no allocation wrapper: gcc -fanalyzer loses the
   * capacity<->count link through wrappers and reports phantom overflows
   * (CLAUDE.md, diam audit 07-24) */
  unsigned char *depth = malloc((size_t)t->n);
  if (depth == NULL) return HZ_OCT_E_IO;
  memset(depth, (int)HZ_OCT_UNSEEN, (size_t)t->n);
  depth[0] = 0; /* the root is reachable by definition */

  int rc = HZ_OCT_OK;
  for (int32_t i = 0; i < t->n && rc == HZ_OCT_OK; i++) {
    if (depth[i] == HZ_OCT_UNSEEN) {
      rc = HZ_OCT_E_UNREACHABLE;
      break;
    }
    int32_t c = t->nodes[i].child0;
    if (c == -1) continue; /* leaf */
    if (c < 0) {
      rc = HZ_OCT_E_CHILD_NEG;
      break;
    }
    if (c <= i) {
      rc = HZ_OCT_E_CYCLE;
      break;
    }
    if (c > t->n - 8) { /* i.e. c + 8 > n, written so it cannot overflow */
      rc = HZ_OCT_E_RANGE;
      break;
    }
    if ((c - 1) % 8 != 0) {
      rc = HZ_OCT_E_ALIGN;
      break;
    }
    if ((int)depth[i] + 1 > t->log2size) { /* a node of size 1 has no children */
      rc = HZ_OCT_E_DEPTH;
      break;
    }
    for (int j = 0; j < 8; j++) {
      if (depth[c + j] != HZ_OCT_UNSEEN) {
        rc = HZ_OCT_E_TWOPARENT;
        break;
      }
      depth[c + j] = (unsigned char)(depth[i] + 1);
    }
  }
  free(depth);
  return rc;
}

/* --- serialization -------------------------------------------------------- */

static const char hz_oct_magic[8] = {'H', 'Z', 'O', 'C', 'T', 'R', 'E', '\0'};
#define HZ_OCT_ENDIAN_WITNESS 0x04030201u
/* bit pattern of the IEEE754 double 1.0; the writer stores an actual 1.0, so a
 * build with another floating-point representation fails HERE instead of
 * reading the node values as garbage */
#define HZ_OCT_CANARY_BITS UINT64_C(0x3FF0000000000000)

/* elements per bulk read/write: 1024 int32 + 1024 double = 12 KiB of stack,
 * inside a 32 KiB L1, and it keeps peak memory at the arena — gathering all
 * three arrays at once would cost 1.3 GB on top of a 1.5 GB arena at n = 2^26 */
enum { HZ_OCT_IOCHUNK = 1024 };

/* One array per call, one chunked loop per call. THE SPLIT IS NOT COSMETIC:
 * with all three loops in hz_oct_load, gcc -fanalyzer widens after the first
 * chunk iteration, loses the calloc-capacity <-> loop-bound link and reports a
 * phantom heap overflow (the class CLAUDE.md records from the diam audit). One
 * loop per function keeps each state space small enough to stay exact. */
static int read_child0(FILE *f, hz_onode *nodes, int32_t n) {
  int32_t buf[HZ_OCT_IOCHUNK];
  for (int32_t i = 0; i < n;) {
    int32_t m = n - i;
    if (m > HZ_OCT_IOCHUNK) m = HZ_OCT_IOCHUNK;
    if (fread(buf, sizeof buf[0], (size_t)m, f) != (size_t)m) return HZ_OCT_E_IO;
    for (int32_t j = 0; j < m; j++)
      nodes[i + j].child0 = buf[j];
    i += m;
  }
  return HZ_OCT_OK;
}

/* part 0 = real, 1 = imaginary */
static int read_part(FILE *f, hz_onode *nodes, int32_t n, int part) {
  double buf[HZ_OCT_IOCHUNK];
  for (int32_t i = 0; i < n;) {
    int32_t m = n - i;
    if (m > HZ_OCT_IOCHUNK) m = HZ_OCT_IOCHUNK;
    if (fread(buf, sizeof buf[0], (size_t)m, f) != (size_t)m) return HZ_OCT_E_IO;
    for (int32_t j = 0; j < m; j++)
      nodes[i + j].k2 =
          part == 0 ? CMPLX(buf[j], cimag(nodes[i + j].k2)) : CMPLX(creal(nodes[i + j].k2), buf[j]);
    i += m;
  }
  return HZ_OCT_OK;
}

static int write_child0(FILE *f, const hz_onode *nodes, int32_t n) {
  int32_t buf[HZ_OCT_IOCHUNK];
  for (int32_t i = 0; i < n;) {
    int32_t m = n - i;
    if (m > HZ_OCT_IOCHUNK) m = HZ_OCT_IOCHUNK;
    for (int32_t j = 0; j < m; j++)
      buf[j] = nodes[i + j].child0;
    if (fwrite(buf, sizeof buf[0], (size_t)m, f) != (size_t)m) return HZ_OCT_E_IO;
    i += m;
  }
  return HZ_OCT_OK;
}

static int write_part(FILE *f, const hz_onode *nodes, int32_t n, int part) {
  double buf[HZ_OCT_IOCHUNK];
  for (int32_t i = 0; i < n;) {
    int32_t m = n - i;
    if (m > HZ_OCT_IOCHUNK) m = HZ_OCT_IOCHUNK;
    for (int32_t j = 0; j < m; j++)
      buf[j] = part == 0 ? creal(nodes[i + j].k2) : cimag(nodes[i + j].k2);
    if (fwrite(buf, sizeof buf[0], (size_t)m, f) != (size_t)m) return HZ_OCT_E_IO;
    i += m;
  }
  return HZ_OCT_OK;
}

int hz_oct_hdr_decode(const unsigned char *h, int *log2size, int32_t *n) {
  if (memcmp(h + HZ_OCT_OFF_MAGIC, hz_oct_magic, sizeof hz_oct_magic) != 0) return HZ_OCT_E_MAGIC;
  uint32_t u;
  memcpy(&u, h + HZ_OCT_OFF_VERSION, sizeof u);
  if (u != HZ_OCT_VERSION) return HZ_OCT_E_VERSION;
  memcpy(&u, h + HZ_OCT_OFF_FLAGS, sizeof u);
  if (u != 0u) return HZ_OCT_E_FLAGS; /* fail closed on any unknown bit */
  memcpy(&u, h + HZ_OCT_OFF_ENDIAN, sizeof u);
  if (u != HZ_OCT_ENDIAN_WITNESS) return HZ_OCT_E_ENDIAN;
  uint64_t cb;
  memcpy(&cb, h + HZ_OCT_OFF_CANARY, sizeof cb);
  if (cb != HZ_OCT_CANARY_BITS) return HZ_OCT_E_CANARY;
  int32_t rsv;
  memcpy(&rsv, h + HZ_OCT_OFF_RESERVED, sizeof rsv);
  if (rsv != 0) return HZ_OCT_E_RESERVED;
  int32_t l2, nn;
  memcpy(&l2, h + HZ_OCT_OFF_LOG2SIZE, sizeof l2);
  memcpy(&nn, h + HZ_OCT_OFF_N, sizeof nn);
  if (l2 < 0 || l2 > HZ_OCT_MAX_LOG2SIZE) return HZ_OCT_E_LOG2SIZE;
  if (nn < 1 || nn > HZ_OCT_MAX_NODES || (nn - 1) % 8 != 0) return HZ_OCT_E_NCOUNT;
  *log2size = (int)l2;
  *n = nn;
  return HZ_OCT_OK;
}

int hz_oct_save(const hz_octree *t, FILE *f) {
  unsigned char h[HZ_OCT_HDRSIZE];
  memset(h, 0, sizeof h);
  memcpy(h + HZ_OCT_OFF_MAGIC, hz_oct_magic, sizeof hz_oct_magic);
  uint32_t u = HZ_OCT_VERSION;
  memcpy(h + HZ_OCT_OFF_VERSION, &u, sizeof u);
  u = 0u; /* flags */
  memcpy(h + HZ_OCT_OFF_FLAGS, &u, sizeof u);
  u = HZ_OCT_ENDIAN_WITNESS;
  memcpy(h + HZ_OCT_OFF_ENDIAN, &u, sizeof u);
  int32_t v = (int32_t)t->log2size;
  memcpy(h + HZ_OCT_OFF_LOG2SIZE, &v, sizeof v);
  v = t->n;
  memcpy(h + HZ_OCT_OFF_N, &v, sizeof v);
  v = 0; /* reserved */
  memcpy(h + HZ_OCT_OFF_RESERVED, &v, sizeof v);
  double canary = 1.0;
  memcpy(h + HZ_OCT_OFF_CANARY, &canary, sizeof canary);
  if (fwrite(h, sizeof h, 1, f) != 1) return HZ_OCT_E_IO;

  int rc = write_child0(f, t->nodes, t->n);
  if (rc == HZ_OCT_OK) rc = write_part(f, t->nodes, t->n, 0);
  if (rc == HZ_OCT_OK) rc = write_part(f, t->nodes, t->n, 1);
  return rc;
}

int hz_oct_load(hz_octree *t, FILE *f) {
  t->nodes = NULL;
  t->n = t->cap = 0;
  t->log2size = 0;

  unsigned char h[HZ_OCT_HDRSIZE];
  if (fread(h, sizeof h, 1, f) != 1) return HZ_OCT_E_IO;
  int log2size = 0;
  int32_t n = 0;
  int rc = hz_oct_hdr_decode(h, &log2size, &n);
  if (rc != HZ_OCT_OK) return rc;

  /* n is inside HZ_OCT_MAX_NODES before this point, so the product cannot
   * overflow size_t on any platform where the arena could be allocated at all */
  hz_onode *nodes = calloc((size_t)n, sizeof(hz_onode));
  if (nodes == NULL) return HZ_OCT_E_IO;

  /* real parts first, then imaginary: the file keeps the two apart, so the node
   * is filled in two passes rather than the file carrying an interleaved copy
   * of whatever the compiler chose for double complex */
  if (read_child0(f, nodes, n) != HZ_OCT_OK || read_part(f, nodes, n, 0) != HZ_OCT_OK ||
      read_part(f, nodes, n, 1) != HZ_OCT_OK) {
    free(nodes);
    return HZ_OCT_E_IO;
  }

  t->nodes = nodes;
  t->n = t->cap = n;
  t->log2size = log2size;
  rc = hz_oct_validate(t);
  if (rc != HZ_OCT_OK) {
    hz_oct_free(t);
    return rc;
  }
  return HZ_OCT_OK;
}
