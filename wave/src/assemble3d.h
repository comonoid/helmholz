#ifndef HZ_ASSEMBLE3D_H
#define HZ_ASSEMBLE3D_H

/* M5: matrix-free Galerkin entries over tensor-product potentials.
 * Phi_{l,n}(x,y,z) = phi(x/h - nx) phi(y/h - ny) phi(z/h - nz), h = 2^l.
 * <Phi_a, L Phi_b> = Laplacian part (separable, background-free) +
 * sum over octree leaves of k2 * product of per-axis integrals.
 * All 1D integrals are exact (src/phi.c) and memoized with translation
 * canonicalization — the hot path does hash lookups, not integration. */

#include "octree.h"
#include "phi.h"
#include <complex.h>

typedef struct {
  int lvl;
  int n[3];
} hz_pot3;

/* <Phi_a, (Lap + k2(x)) Phi_b> over [0,dom]^3 (dom in cells, per axis) */
double complex hz_entry3d(const hz_octree *t, const int dom[3], hz_pot3 a, hz_pot3 b);

/* <Phi_a, f> for a blob source f = amp * Phi_src (tensor potential shape) */
double complex hz_rhs3d(const int dom[3], hz_pot3 a, hz_pot3 src, double complex amp);

/* memo statistics (for the M8 report) */
void hz_asm3d_memo_stats(long *hits, long *misses);

#endif
