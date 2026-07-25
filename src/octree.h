#ifndef HZ_OCTREE_H
#define HZ_OCTREE_H

/* Stage-1 octree: cubic voxels, piecewise-constant complex k^2 per leaf.
 * Cube [0, 2^log2size)^3 in unit cells. Nodes live in one arena, children
 * referenced by index (8 contiguous) — no per-node malloc (CLAUDE.md: alloc
 * wrappers vs gcc-analyzer, plus locality). Node: child0 < 0 = leaf. */

#include <complex.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
  int32_t child0; /* index of first of 8 children, or -1 for a leaf */
  double complex k2;
} hz_onode;

typedef struct {
  hz_onode *nodes;
  int32_t n, cap;
  int log2size;
} hz_octree;

int hz_oct_init(hz_octree *t, int log2size, double complex background);
void hz_oct_free(hz_octree *t);

/* set k2 on the cell box [lo, hi) (clipped to the cube); returns 0 on ok */
int hz_oct_set_box(hz_octree *t, const int lo[3], const int hi[3], double complex k2);
/* set k2 on cells whose centers lie inside the ball */
int hz_oct_set_ball(hz_octree *t, const double c[3], double r, double complex k2);

double complex hz_oct_at(const hz_octree *t, int x, int y, int z);

/* visit leaves intersecting [lo, hi): cb gets the CLIPPED box [blo, bhi).
 * Nonzero from cb aborts and is returned. */
typedef int (*hz_oct_cb)(void *ctx, const int blo[3], const int bhi[3], double complex k2);
int hz_oct_visit(const hz_octree *t, const int lo[3], const int hi[3], hz_oct_cb cb, void *ctx);

int hz_oct_save(const hz_octree *t, FILE *f);
int hz_oct_load(hz_octree *t, FILE *f); /* into an uninitialized tree */

#endif
