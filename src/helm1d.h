#ifndef HZ_HELM1D_H
#define HZ_HELM1D_H

/* M2: 1D Helmholtz test bench. u'' + k2(x)*u = f on [0, ncell], k2 constant
 * per unit cell (the 1D analog of stage-1 cubic voxels), complex with a small
 * absorbing imaginary part everywhere and absorbing ramps near both ends
 * (radiation condition surrogate). Source f = one potential blob phi(x-xsrc).
 * Solved two ways: independent fine-grid FD reference (Thomas), and the
 * cascade over potential levels (Galerkin + min-norm SVD per level). */

#include "phi.h"
#include <complex.h>

typedef struct {
  int ncell;
  double complex *k2; /* per unit cell [c, c+1) */
} hz_med1d;

/* background k0^2*(1 + i*alpha); ramps: Im part grows quadratically to
 * k0^2*1.0 over ramp_cells at each end */
void hz_med1d_init(hz_med1d *m, int ncell, double k0, double alpha, int ramp_cells);
void hz_med1d_slab(hz_med1d *m, int c0, int c1, double k0, double kfac, double alpha);
void hz_med1d_free(hz_med1d *m);

/* FD reference on grid x_i = i/per_cell, i = 0 .. ncell*per_cell (Dirichlet
 * ends; the ramps absorb before the walls matter). u has npts entries. */
int hz_fd_reference(const hz_med1d *m, int per_cell, double xsrc, double complex *u);

/* One cascade level: potentials phi(x/h - n), h = 2^lvl, n = n0..n1. */
typedef struct {
  int lvl;
  int n0, n1;
} hz_level1d;

typedef struct {
  const hz_level1d *levels;
  int nlev;
  double complex *coef; /* concatenated per level */
  int *ofs;             /* nlev+1 prefix offsets into coef */
} hz_sol1d;

/* Per-level diagnostics (PLAN M2: measure, don't guess). */
typedef struct {
  int dim;         /* system size */
  double cond;     /* s_max/s_min from SVD (includes null modes) */
  double cond_eff; /* s_max/s_rank — conditioning of the solved subspace */
  int nnull;       /* singular values below RCOND*s_max (null modes) */
  double res_drop; /* ||r_after|| / ||r_before|| on this level's test space */
} hz_lvlstat1d;

/* Cascade solve, coarse -> fine, Galerkin test space = level being solved,
 * dense min-norm via LAPACK zgelsd; nsweeps > 1 repeats the level loop as
 * block Gauss-Seidel — required whenever local refinements must feed back
 * into the transport levels (see E3 in tests). sol->coef/ofs are allocated
 * (caller frees with hz_sol1d_free). stats may be NULL (filled by the last
 * sweep). Returns 0 on success. */
int hz_cascade1d_solve(const hz_med1d *m, const hz_level1d *levels, int nlev, double xsrc,
                       int nsweeps, hz_sol1d *sol, hz_lvlstat1d *stats);

double complex hz_sol1d_eval(const hz_sol1d *s, double x);
void hz_sol1d_free(hz_sol1d *s);

/* Assembly entries, exported for the M2 formulation experiments:
 * <Phi_i, L Phi_j> and <Phi_i, f> over [0, ncell], f = phi(x - xsrc).
 * pi/pj must have deriv = 0 (the operator is applied inside). */
double complex hz_galerkin_entry(const hz_med1d *m, hz_phi_factor pi, hz_phi_factor pj);
double complex hz_rhs_entry(const hz_med1d *m, hz_phi_factor pi, double xsrc);

#endif
