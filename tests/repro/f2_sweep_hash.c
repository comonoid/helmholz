/* F2 probe: sweep.c:3521 "The right operand of '^' is a garbage value".
 * Verbatim (reduced) sw_line_slab + sw_travel_rank + sw_build_walk body from
 * src/nstruct/sweep.c with stub types.  Point: the filter loop at line ~393
 * asserts leaf.npcs != 0, but the opaque call sw_line_slab() inside the write
 * loop drops field precision, so the analyzer may re-read lf->npcs as 0 and
 * skip the wpid[] fill loop -> "garbage value" on the order_hash read. */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  int32_t npcs, pcs_first;
} hz_pyr_leaf;
typedef struct {
  int64_t nx, ny, nz;
  int32_t nleaf, nt, ncentleaf;
  hz_pyr_leaf *leaf;
  int64_t *leaf_id;
  int32_t *csr;
} hz_pyr;
typedef struct {
  int64_t key, key2;
  int32_t li;
} hz_sw_pair;
typedef struct {
  int32_t npcs;
  int64_t first, line, cid;
  double cntot;
} hz_sw_wleaf;
typedef struct {
  hz_sw_wleaf *leaf;
  double *wcn, *wn, *kd;
  int32_t *wpid;
  int32_t n;
} hz_sw_walk;

static void sw_free_walk(hz_sw_walk *w) {
  free(w->leaf);
  free(w->wcn);
  free(w->wn);
  free(w->kd);
  free(w->wpid);
  memset(w, 0, sizeof *w);
}
static int sw_cmp_pair(const void *a, const void *b) {
  const hz_sw_pair *x = (const hz_sw_pair *)a, *y = (const hz_sw_pair *)b;
  if (x->key != y->key) return (x->key > y->key) - (x->key < y->key);
  return (x->key2 > y->key2) - (x->key2 < y->key2);
}

/* ---- verbatim sw_line_slab (sweep.c:158-193) ---- */
static void sw_line_slab(const hz_pyr *py, const int sg[3], int64_t id, int64_t *line,
                         int64_t *slab) {
  const int *sv = sg;
  int64_t ii[3], len[3];
  int na[3] = {0, 0, 0}, nnz = 0, ax;
  len[0] = py->nx;
  len[1] = py->ny;
  len[2] = py->nz;
  ii[0] = id % py->nx;
  ii[1] = (id / py->nx) % py->ny;
  ii[2] = id / ((int64_t)py->nx * py->ny);
  for (ax = 0; ax < 3; ax++)
    if (sv[ax] != 0) na[nnz++] = ax;
  if (nnz == 1) {
    int c = (na[0] + 1) % 3, dd = (na[0] + 2) % 3;
    *line = ii[c] + len[c] * ii[dd];
    *slab = (int64_t)sv[na[0]] * ii[na[0]];
  } else if (nnz == 2) {
    int a = na[0], b = na[1], c = 3 - a - b;
    if (c < 0) c = 0;
    if (c > 2) c = 2;
    int64_t o1 = (int64_t)sv[b] * ii[b] - (int64_t)sv[a] * ii[a] + len[a];
    *line = o1 + (len[a] + len[b]) * ii[c];
    *slab = (int64_t)sv[a] * ii[a] + (int64_t)sv[b] * ii[b];
  } else {
    int a = na[0], b = na[1], c = na[2];
    int64_t o1 = (int64_t)sv[b] * ii[b] - (int64_t)sv[a] * ii[a] + len[a];
    int64_t o2 = (int64_t)sv[c] * ii[c] - (int64_t)sv[b] * ii[b] + len[b];
    *line = o1 + (len[a] + len[b]) * o2;
    *slab = (int64_t)sv[a] * ii[a] + (int64_t)sv[b] * ii[b] + (int64_t)sv[c] * ii[c];
  }
}

/* ---- verbatim sw_travel_rank (sweep.c:225-237) ---- */
static int64_t sw_travel_rank(const hz_pyr *py, const int sg[3], int64_t id) {
  int64_t ii[3], len[3], u = 0;
  int ax;
  len[0] = py->nx;
  len[1] = py->ny;
  len[2] = py->nz;
  ii[0] = id % py->nx;
  ii[1] = (id / py->nx) % py->ny;
  ii[2] = id / ((int64_t)py->nx * py->ny);
  for (ax = 0; ax < 3; ax++)
    u += (sg[ax] >= 0) ? ii[ax] : len[ax] - 1 - ii[ax];
  return u;
}

/* ---- verbatim sw_build_walk body (sweep.c:376-427) ---- */
static int sw_build_walk(const hz_pyr *py, const int sg[3], const double *om, const double *area,
                         const double *nrm, const double *kd, hz_sw_walk *w) {
  int32_t li, k;
  int64_t pos = 0;
  hz_sw_pair *pairs = (hz_sw_pair *)malloc((size_t)py->nleaf * sizeof *pairs);
  w->leaf = (hz_sw_wleaf *)malloc((size_t)py->nleaf * sizeof *w->leaf);
  w->wcn = (double *)malloc((size_t)py->nt * sizeof *w->wcn);
  w->wn = (double *)malloc((size_t)py->nt * sizeof *w->wn);
  w->kd = (double *)malloc((size_t)py->nt * sizeof *w->kd);
  w->wpid = (int32_t *)malloc((size_t)py->nt * sizeof *w->wpid);
  if (!pairs || !w->leaf || !w->wcn || !w->wn || !w->kd || !w->wpid) {
    free(pairs);
    sw_free_walk(w);
    return 2;
  }
  {
    int32_t np = 0;
    for (li = 0; li < py->nleaf; li++) {
      if (py->leaf[li].npcs == 0) continue; /* пустой лист — не участник переноса */
      sw_line_slab(py, sg, py->leaf_id[li], &pairs[np].key, &(int64_t){0});
      pairs[np].key2 = sw_travel_rank(py, sg, py->leaf_id[li]);
      pairs[np].li = li;
      np++;
    }
    w->n = np;
  }
  qsort(pairs, (size_t)w->n, sizeof *pairs, sw_cmp_pair);
  for (k = 0; k < w->n; k++) {
    const hz_pyr_leaf *lf = &py->leaf[pairs[k].li];
    int32_t u;
    double cntot = 0.0;
    {
      int64_t slab_dummy;
      sw_line_slab(py, sg, py->leaf_id[pairs[k].li], &w->leaf[k].line, &slab_dummy);
    }
    w->leaf[k].npcs = lf->npcs;
    w->leaf[k].first = pos;
    for (u = 0; u < lf->npcs; u++) {
      int32_t p = py->csr[lf->pcs_first + u];
      double d0 = om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                  om[2] * nrm[3 * (int64_t)p + 2];
      double an = fabs(d0);
      w->wcn[pos] = area[p] * an;
      w->wn[pos] = an;
      w->kd[pos] = kd[p];
      w->wpid[pos] = p;
      cntot += w->wcn[pos];
      pos++;
    }
    w->leaf[k].cntot = cntot;
  }
  free(pairs);
  return 0;
}

/* caller of the SAME shape as hz_sw_run(): external, not called inside the
 * TU -> clang-analyzer analyses it as a TOP-LEVEL root with symbolic args,
 * exactly as it does hz_sw_run() inside sweep.c. */
int probe_run(const hz_pyr *py, const double *area, const double *nrm, const double *kd, int nd) {
  hz_sw_walk *walks;
  uint64_t order_hash = 1469598103934665603ULL;
  double om[3] = {1.0, 0.0, 0.0};
  int sg[3] = {1, 1, 1};
  int d;
  walks = (hz_sw_walk *)malloc((size_t)nd * sizeof *walks);
  if (!walks) return 2;
  for (d = 0; d < nd; d++) {
    if (sw_build_walk(py, sg, om, area, nrm, kd, &walks[d]) != 0) return 3;
    for (int32_t k = 0; k < walks[d].n; k++)
      order_hash = order_hash * 1099511628211ULL ^
                   (uint64_t)(uint32_t)walks[d].wpid[walks[d].leaf[k].first];
  }
  return (int)(order_hash & 1u);
}

int main(void) { return 0; }
