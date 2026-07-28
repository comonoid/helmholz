#ifndef HZ_SOLVER3D_H
#define HZ_SOLVER3D_H

/* M6: cascade solver. Levels coarse->fine on regular dyadic grids covering
 * the whole domain (index range [-1, dom/h+1] per axis). Small levels are
 * solved densely (zgelsd, min-norm); large levels structurally: Toeplitz
 * background applied via padded 3D FFT (exact linear convolution) + sparse
 * deviation matrix V (medium leaves != background, domain-edge clipping),
 * BiCGStab with a guarded circulant preconditioner. Cross-level coupling is
 * exact through the two-scale relation: coarse levels prolong onto the floor
 * grid, residuals restrict back (block Gauss-Seidel sweeps, >= 2 per E3). */

#include "assemble3d.h"
#include "octree.h"
#include <complex.h>

typedef struct {
  hz_pot3 p;
  double complex amp;
} hz_src3;

typedef struct {
  const hz_octree *tree;
  int dom[3]; /* domain in cells, each divisible by the top-level h */
  double complex k2bg;
  const hz_src3 *src;
  int nsrc;
} hz_scene3;

typedef struct {
  int lvl;
  int m[3]; /* index counts per axis, offsets start at -1 */
  double complex *coef;
} hz_slevel;

typedef struct {
  hz_slevel *lv;
  int nlev; /* coarse -> fine; lv[nlev-1] is the floor */
  int dom[3];
  double final_relres; /* floor-level residual after the last sweep */
} hz_sol3;

/* returns 0 on success; verbose!=0 prints per-level dims/iterations/timings */
int hz_solve3d(const hz_scene3 *sc, int top_lvl, int floor_lvl, int nsweeps, int verbose,
               hz_sol3 *out);

/* BiCGStab relative-residual target (default 1e-6; pictures are fine at 1e-3
 * — the field error is dominated by discretization far above that) */
void hz_solver3d_set_tol(double t);

/* test hook: y = (T_fft + V) x on one level — must equal the exact dense
 * Galerkin matrix times x (tests/test_solver3d.c) */
int hz_dbg_matvec(const hz_scene3 *sc, int lvl, const double complex *x, double complex *y);
double complex hz_sol3_eval(const hz_sol3 *s, const double p[3]);
void hz_sol3_free(hz_sol3 *s);

#endif
