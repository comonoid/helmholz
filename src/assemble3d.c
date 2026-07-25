#include "assemble3d.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ---- memoized 1D integrals ----------------------------------------------
 * key: canonicalized (lvlA,dA, lvlB,nB',dB, interval') after translating the
 * larger-h factor to offset 0. All quantities are integers in cell units, so
 * the key is exact. Open-addressing hash, single-threaded. */

#define MEMO_BITS 17
#define MEMO_SIZE (1u << MEMO_BITS)

typedef struct {
  int64_t key[3]; /* full key, no lossy folding */
  double val;
  int used;
} memo_ent;

/* thread-local: hz_entry3d runs under OpenMP; every thread owns a table
 * (values are identical by construction, duplication only costs memory) */
static _Thread_local memo_ent *g_memo;
static _Thread_local long g_hits, g_misses;

void hz_asm3d_memo_stats(long *hits, long *misses) {
  *hits = g_hits;
  *misses = g_misses;
}

static double int1d(int la, int na, int da, int lb, int nb, int db, int c0, int c1) {
  /* canonical translation: shift so the larger-h factor sits at offset 0.
   * The shift s is a multiple of BOTH spacings (the larger h is a power-of-two
   * multiple of the smaller), so the translated offsets stay integers and the
   * key is exact. */
  int ha = 1 << la, hb = 1 << lb;
  int s = (ha >= hb) ? ha * na : hb * nb;
  int fa_n = na - s / ha;
  int fb_n = nb - s / hb;
  int k0 = c0 - s, k1 = c1 - s;

  if (g_memo == NULL) g_memo = calloc(MEMO_SIZE, sizeof(memo_ent));
  if (g_memo != NULL) {
    /* deriv needs 2 bits (values 0 and 2): (lvl<<2)|deriv, never (lvl<<1) —
     * the 1-bit packing aliased "level l, phi''" with "level l+1, phi" */
    int64_t ka = (((int64_t)la << 2) | da) << 32 | (int64_t)(uint32_t)fa_n;
    int64_t kb = (((int64_t)lb << 2) | db) << 32 | (int64_t)(uint32_t)fb_n;
    int64_t kc = ((int64_t)k0 << 32) | (int64_t)(uint32_t)k1;
    uint64_t h = ((uint64_t)ka * 0x9E3779B97F4A7C15ULL) ^ ((uint64_t)kb * 0xC2B2AE3D27D4EB4FULL) ^
                 ((uint64_t)kc * 0x165667B19E3779F9ULL);
    uint32_t idx = (uint32_t)(h >> (64 - MEMO_BITS));
    for (uint32_t probe = 0; probe < 64; probe++) {
      memo_ent *e = &g_memo[(idx + probe) & (MEMO_SIZE - 1)];
      if (!e->used) {
        hz_phi_factor fa = {(double)ha, (double)fa_n, da};
        hz_phi_factor fb = {(double)hb, (double)fb_n, db};
        double v = hz_phi_prod_integral((double)k0, (double)k1, fa, fb);
        e->key[0] = ka;
        e->key[1] = kb;
        e->key[2] = kc;
        e->val = v;
        e->used = 1;
        g_misses++;
        return v;
      }
      if (e->key[0] == ka && e->key[1] == kb && e->key[2] == kc) {
        g_hits++;
        return e->val;
      }
    }
  }
  /* table degraded or full neighborhood: compute directly */
  hz_phi_factor fa = {(double)ha, (double)fa_n, da};
  hz_phi_factor fb = {(double)hb, (double)fb_n, db};
  g_misses++;
  return hz_phi_prod_integral((double)k0, (double)k1, fa, fb);
}

/* per-axis overlap of two potentials clipped to [0, dom] (in cells) */
static int axis_overlap(int la, int na, int lb, int nb, int dom, int *c0, int *c1) {
  int ha = 1 << la, hb = 1 << lb;
  int lo = ha * (na - 2), hi = ha * (na + 2);
  int lo2 = hb * (nb - 2), hi2 = hb * (nb + 2);
  if (lo2 > lo) lo = lo2;
  if (hi2 < hi) hi = hi2;
  if (lo < 0) lo = 0;
  if (hi > dom) hi = dom;
  *c0 = lo;
  *c1 = hi;
  return lo < hi;
}

typedef struct {
  hz_pot3 a, b;
  double complex acc;
} k2ctx;

static int k2_cb(void *vctx, const int blo[3], const int bhi[3], double complex k2) {
  k2ctx *c = vctx;
  double m = 1.0;
  for (int ax = 0; ax < 3; ax++) {
    m *= int1d(c->a.lvl, c->a.n[ax], 0, c->b.lvl, c->b.n[ax], 0, blo[ax], bhi[ax]);
    if (!(fabs(m) > 0.0)) return 0; /* early out: no overlap on this axis */
  }
  c->acc += k2 * m;
  return 0;
}

double complex hz_entry3d(const hz_octree *t, const int dom[3], hz_pot3 a, hz_pot3 b) {
  int c0[3], c1[3];
  for (int ax = 0; ax < 3; ax++)
    if (!axis_overlap(a.lvl, a.n[ax], b.lvl, b.n[ax], dom[ax], &c0[ax], &c1[ax])) return 0.0;

  /* Laplacian: sum over axes of D_ax * prod of masses on other axes */
  double D[3], M[3];
  for (int ax = 0; ax < 3; ax++) {
    D[ax] = int1d(a.lvl, a.n[ax], 0, b.lvl, b.n[ax], 2, c0[ax], c1[ax]);
    M[ax] = int1d(a.lvl, a.n[ax], 0, b.lvl, b.n[ax], 0, c0[ax], c1[ax]);
  }
  double lap = D[0] * M[1] * M[2] + M[0] * D[1] * M[2] + M[0] * M[1] * D[2];

  /* k2 part: exact sum over octree leaves intersecting the overlap box */
  k2ctx ctx = {a, b, 0.0};
  hz_oct_visit(t, c0, c1, k2_cb, &ctx);
  return lap + ctx.acc;
}

double complex hz_rhs3d(const int dom[3], hz_pot3 a, hz_pot3 src, double complex amp) {
  double m = 1.0;
  for (int ax = 0; ax < 3; ax++) {
    int c0, c1;
    if (!axis_overlap(a.lvl, a.n[ax], src.lvl, src.n[ax], dom[ax], &c0, &c1)) return 0.0;
    m *= int1d(a.lvl, a.n[ax], 0, src.lvl, src.n[ax], 0, c0, c1);
  }
  return amp * m;
}
