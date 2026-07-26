/* TREFFTZ-DG: the reformulation the slab2d bench forced.
 *
 * WHAT slab2d PROVED. With a CONTINUOUS (partition-of-unity) carrier basis the
 * volume Galerkin residual is VACUOUS at k W >> 1: the carrier cancels k^2
 * analytically, so a plane wave — the whole solution — sits in the kernel of the
 * volume operator. Measured: the exact solution satisfies the assembled system
 * to 3e-13 while LSQR drives the residual to 9e-8 and walks AWAY from it (field
 * error 0.84 -> 1.13 as the residual falls). Independent of the box size
 * (NE = 8/16/24/32 gave 0.88/1.03/1.05/1.05) and reproduced in the simplest
 * possible setting — one medium, one plane wave, one true direction per element.
 *
 * WHY IT IS STRUCTURAL AND NOT A BUG. A C^1 basis has continuity BUILT IN, so
 * continuity yields no equations; and the volume residual yields none either,
 * because every basis combination nearly solves the PDE. That leaves only the
 * outer boundary: NE^2 * ND unknowns against O(NE * ND) conditions. A
 * DISCONTINUOUS basis is what supplies interior equations — one set per face,
 * ~2 NE^2 faces — and that is why every Trefftz method in the literature is a DG
 * method.
 *
 * THIS BENCH. Cells do not overlap; on each cell the basis is PURE plane waves,
 * so every basis function solves the Helmholtz equation exactly and there is no
 * volume integral at all. All equations live on the skeleton:
 *     interior face:  [u] = 0   and   [dn u]/(i k) = 0
 *     boundary face:  dn u - i k u = h                    (impedance, data exact)
 * The two interior conditions carry the SAME units, so the relative weight is
 * derived rather than tuned — there is no penalty parameter anywhere.
 *
 * Every integral is elementary: a segment integral of exp(i(k_j - k_m).x). No
 * phi, no polygons, no quadrature.
 *
 * TEST FUNCTIONS ARE CONJUGATED, unlike the rest of the project's bilinear
 * assembly. Here the tests on a face are the traces of the same plane waves that
 * make up the jump; testing bilinearly would give a symmetric (not Hermitian)
 * Gram which can be singular, and "all rows zero" would then not mean "the jump
 * is zero" — the condition would be vacuous in a new way. */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { MAXCELL = 4096, MAXDIR = 64, MAXDIM = 60000 };
static const double LAM = 1.0;

typedef struct {
  double W, k0, theta;
  int nd, ne;
  int oracle; /* 1 = the exact direction is in the fan */
} cfg;

/* Integral over the segment P0->P1 of exp(i g.x) ds, closed form. */
static double complex seg_exp(double gx, double gy, double x0, double y0, double x1, double y1) {
  double dx = x1 - x0, dy = y1 - y0;
  double len = sqrt(dx * dx + dy * dy);
  if (!(len > 0.0)) return 0.0;
  double gam = gx * dx + gy * dy; /* phase across the whole segment */
  double complex e0 = cexp(CMPLX(0.0, 1.0) * (gx * x0 + gy * y0));
  if (fabs(gam) < 1e-8) {
    /* series, because (e^{ig}-1)/(ig) loses everything at small g */
    double complex s = 1.0, t = 1.0;
    for (int m = 1; m < 12; m++) {
      t *= CMPLX(0.0, 1.0) * gam / (double)(m + 1);
      s += t;
    }
    return len * e0 * s;
  }
  return len * e0 * (cexp(CMPLX(0.0, 1.0) * gam) - 1.0) / (CMPLX(0.0, 1.0) * gam);
}

int main(int argc, char **argv) {
  static const char *const KEYS[] = {"W", "nd", "ne", "th", "it", "oracle", NULL};
  for (int i = 1; i < argc; i++) {
    const char *eq = strchr(argv[i], '=');
    int ok = 0;
    if (eq)
      for (int k = 0; KEYS[k]; k++)
        if ((size_t)(eq - argv[i]) == strlen(KEYS[k]) &&
            !strncmp(argv[i], KEYS[k], strlen(KEYS[k])))
          ok = 1;
    if (!ok) {
      printf("tdg2d: unknown argument '%s'; keys are", argv[i]);
      for (int k = 0; KEYS[k]; k++)
        printf(" %s", KEYS[k]);
      printf("\n");
      return 1;
    }
  }
  cfg c = {1.0 * LAM, 2.0 * M_PI / LAM, 0.3, 8, 8, 0};
  int itmax = 4000;
  for (int i = 1; i < argc; i++) {
    double v = atof(strchr(argv[i], '=') + 1);
    if (!strncmp(argv[i], "W=", 2)) c.W = v * LAM;
    if (!strncmp(argv[i], "nd=", 3)) c.nd = (int)v;
    if (!strncmp(argv[i], "ne=", 3)) c.ne = (int)v;
    if (!strncmp(argv[i], "th=", 3)) c.theta = v;
    if (!strncmp(argv[i], "it=", 3)) itmax = (int)v;
    if (!strncmp(argv[i], "oracle=", 7)) c.oracle = (int)v;
  }
  if (c.nd > MAXDIR || c.ne * c.ne > MAXCELL) return 1;
  int ncell = c.ne * c.ne, dim = ncell * c.nd;
  if (dim > MAXDIM) {
    printf("tdg2d: dim %d too large\n", dim);
    return 1;
  }
  double L = 0.5 * (double)c.ne * c.W;
  printf("tdg2d: W=%.4g lam  kW=%.4g  ND=%d  NE=%d  box=%.4g lam  dim=%d\n", c.W / LAM, c.k0 * c.W,
         c.nd, c.ne, 2.0 * L / LAM, dim);

  /* directions: a uniform fan, offset so the exact direction is NOT in it unless
   * asked for — otherwise the exact solution is reproduced trivially */
  static double dkx[MAXDIR], dky[MAXDIR];
  for (int d = 0; d < c.nd; d++) {
    double th = 2.0 * M_PI * ((double)d + 0.5) / (double)c.nd;
    if (c.oracle && d == 0) th = c.theta;
    dkx[d] = c.k0 * cos(th);
    dky[d] = c.k0 * sin(th);
  }
  double ekx = c.k0 * cos(c.theta), eky = c.k0 * sin(c.theta);

  /* cell centres; the basis of cell q is exp(i k_d . (x - x_q)) */
  static double cx[MAXCELL], cy[MAXCELL];
  for (int a = 0; a < c.ne; a++)
    for (int bb = 0; bb < c.ne; bb++) {
      int q = a * c.ne + bb;
      cx[q] = -L + ((double)a + 0.5) * c.W;
      cy[q] = -L + ((double)bb + 0.5) * c.W;
    }

  size_t cap = 32u << 20;
  int *ja = malloc(cap * sizeof(int));
  double complex *va = malloc(cap * sizeof(double complex));
  size_t *rp = malloc(200000 * sizeof(size_t));
  double complex *rhs = calloc(200000, sizeof(double complex));
  if (!ja || !va || !rp || !rhs) return 1;
  size_t nnz = 0;
  int nrow = 0;
  double complex i1 = CMPLX(0.0, 1.0);

  /* ---- faces ------------------------------------------------------------
   * For each face: the two adjacent cells (or one, at the boundary), the
   * segment, and the normal pointing from cell A to cell B. */
  for (int dir = 0; dir < 2; dir++)
    for (int a = 0; a < c.ne; a++)
      for (int bb = 0; bb <= c.ne; bb++) {
        int qa = -1, qb = -1;
        double x0, y0, x1, y1, nx, ny;
        if (dir == 0) { /* vertical face at x = -L + bb*W */
          double xf = -L + (double)bb * c.W;
          x0 = xf;
          y0 = -L + (double)a * c.W;
          x1 = xf;
          y1 = y0 + c.W;
          nx = 1.0;
          ny = 0.0;
          if (bb > 0) qa = (bb - 1) * c.ne + a;
          if (bb < c.ne) qb = bb * c.ne + a;
        } else { /* horizontal face at y = -L + bb*W */
          double yf = -L + (double)bb * c.W;
          x0 = -L + (double)a * c.W;
          y0 = yf;
          x1 = x0 + c.W;
          y1 = yf;
          nx = 0.0;
          ny = 1.0;
          if (bb > 0) qa = a * c.ne + (bb - 1);
          if (bb < c.ne) qb = a * c.ne + bb;
        }
        /* tests on this face: the traces of every basis function of the
         * adjacent cells, conjugated */
        int qs[2] = {qa, qb};
        int interior = (qa >= 0 && qb >= 0);
        for (int ti = 0; ti < 2; ti++) {
          if (qs[ti] < 0) continue;
          for (int tm = 0; tm < c.nd; tm++) {
            int qt = qs[ti];
            /* conj(T) = exp(-i k_tm .(x - x_qt)) */
            double tgx = -dkx[tm], tgy = -dky[tm];
            double complex tph = cexp(-i1 * (tgx * cx[qt] + tgy * cy[qt]));
            int nc = interior ? 2 : 1;
            for (int cond = 0; cond < nc; cond++) {
              rp[nrow] = nnz;
              double complex acc = 0.0;
              for (int ti2 = 0; ti2 < 2; ti2++) {
                int q = qs[ti2];
                if (q < 0) continue;
                double sgn = interior ? (ti2 == 0 ? -1.0 : 1.0) : 1.0; /* [.] = B - A */
                for (int d = 0; d < c.nd; d++) {
                  double gx = dkx[d] + tgx, gy = dky[d] + tgy;
                  double complex ph = cexp(-i1 * (dkx[d] * cx[q] + dky[d] * cy[q]));
                  double complex Iv = seg_exp(gx, gy, x0, y0, x1, y1) * ph * tph;
                  double complex w;
                  if (!interior) {
                    /* impedance: dn u - i k u */
                    w = (i1 * (dkx[d] * nx + dky[d] * ny) - i1 * c.k0) * Iv;
                  } else if (cond == 0) {
                    w = sgn * Iv; /* [u] */
                  } else {
                    w = sgn * i1 * (dkx[d] * nx + dky[d] * ny) * Iv / (i1 * c.k0);
                  }
                  if (!(cabs(w) > 0.0) || nnz >= cap) continue;
                  ja[nnz] = q * c.nd + d;
                  va[nnz++] = w;
                }
              }
              (void)acc;
              if (!interior) {
                /* drive: the same impedance combination of the exact wave */
                double gx = ekx + tgx, gy = eky + tgy;
                double complex Iv = seg_exp(gx, gy, x0, y0, x1, y1) * tph;
                rhs[nrow] = (i1 * (ekx * nx + eky * ny) - i1 * c.k0) * Iv;
              } else {
                rhs[nrow] = 0.0;
              }
              nrow++;
            }
          }
        }
      }
  rp[nrow] = nnz;
  printf("  rows=%d  unknowns=%d  (%.1fx overdetermined)  nnz=%zu\n", nrow, dim,
         (double)nrow / (double)dim, nnz);

  /* row equilibration */
  for (int r = 0; r < nrow; r++) {
    double s = 0.0;
    for (size_t p = rp[r]; p < rp[r + 1]; p++)
      s += cabs(va[p]) * cabs(va[p]);
    s = sqrt(s);
    if (!(s > 0.0)) continue;
    for (size_t p = rp[r]; p < rp[r + 1]; p++)
      va[p] /= s;
    rhs[r] /= s;
  }

  /* the exact solution IS in the span when oracle=1: coefficient 1 on direction
   * 0 of every cell, with the cell's phase reference */
  double complex *ce = calloc((size_t)dim, sizeof(double complex));
  if (!ce) return 1;
  if (c.oracle)
    for (int q = 0; q < ncell; q++)
      ce[q * c.nd + 0] = cexp(i1 * (ekx * cx[q] + eky * cy[q]));
  if (c.oracle) {
    double e2 = 0.0, nb = 0.0;
    for (int r = 0; r < nrow; r++) {
      double complex s = 0.0;
      for (size_t p = rp[r]; p < rp[r + 1]; p++)
        s += va[p] * ce[ja[p]];
      e2 += cabs(s - rhs[r]) * cabs(s - rhs[r]);
      nb += cabs(rhs[r]) * cabs(rhs[r]);
    }
    printf("  [diag] exact solution in assembled rows: |Ac-b| = %.3e of |b| = %.3e\n", sqrt(e2),
           sqrt(nb));
  }

  /* ---- LSQR ------------------------------------------------------------- */
  double complex *xs = calloc((size_t)dim, sizeof(double complex));
  double complex *uu = calloc((size_t)nrow, sizeof(double complex));
  double complex *vv = calloc((size_t)dim, sizeof(double complex));
  double complex *ww = calloc((size_t)dim, sizeof(double complex));
  double complex *tv = calloc((size_t)dim, sizeof(double complex));
  if (!xs || !uu || !vv || !ww || !tv) return 1;
  for (int r = 0; r < nrow; r++)
    uu[r] = rhs[r];
  double beta = 0.0;
  for (int r = 0; r < nrow; r++)
    beta += cabs(uu[r]) * cabs(uu[r]);
  beta = sqrt(beta);
  for (int r = 0; r < nrow; r++)
    uu[r] /= beta;
  for (int j = 0; j < dim; j++)
    vv[j] = 0.0;
  for (int r = 0; r < nrow; r++)
    for (size_t p = rp[r]; p < rp[r + 1]; p++)
      vv[ja[p]] += conj(va[p]) * uu[r];
  double alpha = 0.0;
  for (int j = 0; j < dim; j++)
    alpha += cabs(vv[j]) * cabs(vv[j]);
  alpha = sqrt(alpha);
  for (int j = 0; j < dim; j++) {
    vv[j] /= alpha;
    ww[j] = vv[j];
  }
  double phibar = beta, rhobar = alpha, b0 = beta;
  int n2 = 0, n3 = 0, n4 = 0, n6 = 0, n8 = 0;
  for (int it = 0; it < itmax; it++) {
    for (int r = 0; r < nrow; r++) {
      double complex s = 0.0;
      for (size_t p = rp[r]; p < rp[r + 1]; p++)
        s += va[p] * vv[ja[p]];
      uu[r] = s - alpha * uu[r];
    }
    beta = 0.0;
    for (int r = 0; r < nrow; r++)
      beta += cabs(uu[r]) * cabs(uu[r]);
    beta = sqrt(beta);
    if (!(beta > 0.0)) break;
    for (int r = 0; r < nrow; r++)
      uu[r] /= beta;
    for (int j = 0; j < dim; j++)
      tv[j] = 0.0;
    for (int r = 0; r < nrow; r++)
      for (size_t p = rp[r]; p < rp[r + 1]; p++)
        tv[ja[p]] += conj(va[p]) * uu[r];
    alpha = 0.0;
    for (int j = 0; j < dim; j++) {
      vv[j] = tv[j] - beta * vv[j];
      alpha += cabs(vv[j]) * cabs(vv[j]);
    }
    alpha = sqrt(alpha);
    if (!(alpha > 0.0)) break;
    for (int j = 0; j < dim; j++)
      vv[j] /= alpha;
    double rho = sqrt(rhobar * rhobar + beta * beta);
    double cs = rhobar / rho, sn = beta / rho;
    double theta = sn * alpha;
    rhobar = -cs * alpha;
    double phi = cs * phibar;
    phibar = sn * phibar;
    for (int j = 0; j < dim; j++) {
      xs[j] += (phi / rho) * ww[j];
      ww[j] = vv[j] - (theta / rho) * ww[j];
    }
    double rr = phibar / b0;
    if (!n2 && rr <= 1e-2) n2 = it + 1;
    if (!n3 && rr <= 1e-3) n3 = it + 1;
    if (!n4 && rr <= 1e-4) n4 = it + 1;
    if (!n6 && rr <= 1e-6) n6 = it + 1;
    if (!n8 && rr <= 1e-8) n8 = it + 1;
  }
  printf("  LSQR |r|/|b| = %.3e   iters to 1e-2/1e-3/1e-4/1e-6/1e-8: %d/%d/%d/%d/%d\n", phibar / b0,
         n2, n3, n4, n6, n8);

  /* ---- error against the exact plane wave, cell centres and quarter points */
  double num = 0.0, den = 0.0;
  for (int q = 0; q < ncell; q++)
    for (int s = 0; s < 9; s++) {
      double ox = ((double)(s % 3) - 1.0) * 0.3 * c.W, oy = ((double)(s / 3) - 1.0) * 0.3 * c.W;
      double x = cx[q] + ox, y = cy[q] + oy;
      double complex got = 0.0;
      for (int d = 0; d < c.nd; d++)
        got += xs[q * c.nd + d] * cexp(i1 * (dkx[d] * (x - cx[q]) + dky[d] * (y - cy[q])));
      double complex ex = cexp(i1 * (ekx * x + eky * y));
      num += cabs(got - ex) * cabs(got - ex);
      den += cabs(ex) * cabs(ex);
    }
  printf("  FIELD ERROR = %.4e\n", sqrt(num / den));

  free(ja);
  free(va);
  free(rp);
  free(rhs);
  free(ce);
  free(xs);
  free(uu);
  free(vv);
  free(ww);
  free(tv);
  return 0;
}
