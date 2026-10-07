/* F3 probe (replaces build/probe_triage/repro_1194_fp.c, which no longer
 * reproduces): tools/pgather.c "Assigned value is uninitialized" at
 * `g->tnext[t] = g->dhead[s];`. Shape of pg924_grid_build(): ncell from
 * ternary-clamped spans, dhead=malloc(ncell), init loop `s < ncell`,
 * then dhead[s] is read. Top-level root -> symbolic fields. */
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
  const double *vx;
  uint32_t ntris;
} scene_t;
typedef struct {
  int64_t ncell, gnx, gny, gnz;
  double gorg[3], ginv;
  double *rmax, gmaxr;
  int32_t *dhead, *tnext;
} grid_t;

int probe_grid(const scene_t *S) {
  grid_t g;
  double lo[3] = {0, 0, 0}, hi[3] = {1, 1, 1};
  int64_t span[3];
  int a;
  if (S->ntris == 0) return 1;
  for (a = 0; a < 3; a++) g.gorg[a] = lo[a];
  g.ginv = 1.0;
  for (a = 0; a < 3; a++) span[a] = (int64_t)((hi[a] - g.gorg[a]) * g.ginv) + 2;
  g.gnx = span[0] > 1 ? span[0] : 1;
  g.gny = span[1] > 1 ? span[1] : 1;
  g.gnz = span[2] > 1 ? span[2] : 1;
  g.ncell = g.gnx * g.gny * g.gnz;
  g.dhead = (int32_t *)malloc((size_t)g.ncell * sizeof *g.dhead);
  g.tnext = (int32_t *)malloc((size_t)S->ntris * sizeof *g.tnext);
  if (!g.dhead || !g.tnext) {
    free(g.dhead);
    free(g.tnext);
    return 2;
  }
  if (g.ncell <= 0) return 3; /* VARIANT: re-assert capacity */
  for (int64_t s = 0; s < g.ncell; s++) g.dhead[s] = -1;
  g.rmax = (double *)malloc((size_t)S->ntris * sizeof *g.rmax);
  if (g.rmax == NULL) {
    free(g.dhead);
    free(g.tnext);
    return 2;
  }
  for (uint32_t t = 0; t < S->ntris; t++) {
    double cx = fabs(S->vx[t]);
    g.rmax[t] = cx;
    int64_t ix = (int64_t)((cx - g.gorg[0]) * g.ginv), iy = 0, iz = 0;
    if (ix < 0) ix = 0;
    if (iy < 0) iy = 0;
    if (iz < 0) iz = 0;
    if (ix >= g.gnx) ix = g.gnx - 1;
    if (iy >= g.gny) iy = g.gny - 1;
    if (iz >= g.gnz) iz = g.gnz - 1;
    int64_t s = ix + g.gnx * (iy + g.gny * iz);
    g.tnext[t] = g.dhead[s];
    g.dhead[s] = (int32_t)t;
  }
  g.gmaxr = 0.0;
  for (uint32_t t = 0; t < S->ntris; t++)
    if (g.rmax[t] > g.gmaxr) g.gmaxr = g.rmax[t];
  return 0;
}
int main(void) { return 0; }
