/* TREFFTZ-DG at the element size the project actually targets, with a material
 * interface. Derivation and numbers: PLAN.md, "ПОСТАНОВКА СМЕНЕНА НА ТРЕФТЦ-DG".
 *
 * WHAT slab2d PROVED AND WHY THIS EXISTS. With a CONTINUOUS (partition-of-unity)
 * carrier basis the volume Galerkin residual is VACUOUS at k W >> 1: the carrier
 * cancels k^2 analytically, so a plane wave — the whole solution — sits in the
 * kernel of the volume operator, while a C^1 basis has continuity built in and so
 * yields no equations either. That leaves NE^2 * ND unknowns against O(NE * ND)
 * boundary conditions. Measured: the exact solution satisfied the assembled
 * system to 3e-13 while LSQR drove the residual to 9e-8 and walked AWAY from it.
 *
 * HERE cells do not overlap and the basis is PURE plane waves, so every basis
 * function solves the Helmholtz equation exactly on a region of ANY shape and
 * there is no volume integral anywhere. Everything lives on the skeleton:
 *     interior face:  [u] = 0   and   [dn u]/(i k) = 0
 *     box wall:       dn u - i k u = h     (impedance; the data come from the
 *                                           exact solution, so termination is
 *                                           NOT what this bench tests)
 * The two interior conditions carry the same units, so there is no penalty
 * parameter anywhere in the formulation.
 *
 * THE MATERIAL INTERFACE IS NOT A SPECIAL CASE. A cell the interface crosses is
 * split into two sub-cells carrying their own wavenumbers, and the interface
 * becomes an ORDINARY FACE between them with the very same two conditions —
 * which are exactly the physical transmission conditions. Nothing is "cut": the
 * apparatus a continuous basis needed (exact phi x phi x plane-wave integration
 * over polygons, a Nitsche term on the cut locus) is not required at all, because
 * there are no volume integrals to restrict. A cut costs only the clipping of
 * SEGMENTS, since every integral is a face integral.
 *
 * TEST FUNCTIONS ARE CONJUGATED, unlike the project's bilinear volume assembly:
 * the tests on a face are traces of the same plane waves that make up the jump,
 * and a bilinear test would give a symmetric rather than Hermitian Gram, so "all
 * rows zero" would stop meaning "the jump is zero". */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { MAXCELL = 4096, MAXSUB = 8192, MAXDIR = 64, MAXDIM = 80000, MAXROW = 400000 };
static const double LAM = 1.0;

typedef struct {
  double W, k0, theta, n, alpha, dth;
  int nd, ne, oracle, lev;
  double lodk;
  double nhx, nhy, thx, thy;
  double complex amp[3];
  double px[3], py[3];
} cfg;

/* --- exact Fresnel solution: incident + reflected outside, transmitted inside */
static void fresnel(cfg *c) {
  c->nhx = cos(c->alpha);
  c->nhy = sin(c->alpha);
  c->thx = -sin(c->alpha);
  c->thy = cos(c->alpha);
  double k1 = c->k0, k2 = c->k0 * c->n;
  double kt = k1 * sin(c->theta);
  double c1 = sqrt(k1 * k1 - kt * kt), c2 = sqrt(k2 * k2 - kt * kt);
  c->amp[0] = 1.0;
  c->amp[1] = (c1 - c2) / (c1 + c2);
  c->amp[2] = 2.0 * c1 / (c1 + c2);
  c->px[0] = kt * c->thx - c1 * c->nhx;
  c->py[0] = kt * c->thy - c1 * c->nhy;
  c->px[1] = kt * c->thx + c1 * c->nhx;
  c->py[1] = kt * c->thy + c1 * c->nhy;
  c->px[2] = kt * c->thx - c2 * c->nhx;
  c->py[2] = kt * c->thy - c2 * c->nhy;
}

static double sdist(const cfg *c, double x, double y) {
  return x * c->nhx + y * c->nhy;
}

static double complex uexact(const cfg *c, double x, double y) {
  double complex i1 = CMPLX(0.0, 1.0);
  if (sdist(c, x, y) > 0.0)
    return c->amp[0] * cexp(i1 * (c->px[0] * x + c->py[0] * y)) +
           c->amp[1] * cexp(i1 * (c->px[1] * x + c->py[1] * y));
  return c->amp[2] * cexp(i1 * (c->px[2] * x + c->py[2] * y));
}

/* Integral over the segment P0->P1 of exp(i g.x) ds, closed form. */
static double complex seg_exp(double gx, double gy, double x0, double y0, double x1, double y1) {
  double dx = x1 - x0, dy = y1 - y0;
  double len = sqrt(dx * dx + dy * dy);
  if (!(len > 0.0)) return 0.0;
  double gam = gx * dx + gy * dy;
  double complex e0 = cexp(CMPLX(0.0, 1.0) * (gx * x0 + gy * y0));
  if (fabs(gam) < 1e-8) { /* (e^{ig}-1)/(ig) loses everything at small g */
    double complex s = 1.0, t = 1.0;
    for (int m = 1; m < 12; m++) {
      t *= CMPLX(0.0, 1.0) * gam / (double)(m + 1);
      s += t;
    }
    return len * e0 * s;
  }
  return len * e0 * (cexp(CMPLX(0.0, 1.0) * gam) - 1.0) / (CMPLX(0.0, 1.0) * gam);
}

/* Clip a segment to one side of the interface: side 0 = medium (s<0), 1 = the
 * rest. This is the ONLY geometry a cut costs here, because every integral in the
 * formulation is a face integral. */
static int clip_side(const cfg *c, int side, double *ax, double *ay, double *bx, double *by) {
  double s0 = sdist(c, *ax, *ay), s1 = sdist(c, *bx, *by);
  int in0 = (side == 0) ? (s0 < 0.0) : (s0 > 0.0);
  int in1 = (side == 0) ? (s1 < 0.0) : (s1 > 0.0);
  if (!in0 && !in1) return 0;
  if (in0 && in1) return 1;
  double t = s0 / (s0 - s1);
  double mx = *ax + t * (*bx - *ax), my = *ay + t * (*by - *ay);
  if (in0) {
    *bx = mx;
    *by = my;
  } else {
    *ax = mx;
    *ay = my;
  }
  return 1;
}

typedef struct {
  int cell, side; /* side 0 = medium, 1 = the rest */
  double cx, cy;  /* phase reference: the cell centre */
  int base, nd;   /* first unknown, and how many directions */
  double kx[MAXDIR], ky[MAXDIR];
} sub;

static sub S[MAXSUB];
static int subof[MAXCELL][2];
static double cellx[MAXCELL], celly[MAXCELL], cellh[MAXCELL];
static int celld[MAXCELL];

int main(int argc, char **argv) {
  static const char *const KEYS[] = {"W",     "nd",  "ne",     "th",  "it",   "n",
                                     "alpha", "dth", "oracle", "lev", "lodk", NULL};
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
  cfg c;
  memset(&c, 0, sizeof c);
  c.W = 1.0 * LAM;
  c.k0 = 2.0 * M_PI / LAM;
  c.theta = 0.3;
  c.n = 1.0;
  c.alpha = 0.4;
  c.nd = 8;
  c.ne = 8;
  int itmax = 4000;
  for (int i = 1; i < argc; i++) {
    double v = atof(strchr(argv[i], '=') + 1);
    if (!strncmp(argv[i], "W=", 2)) c.W = v * LAM;
    if (!strncmp(argv[i], "nd=", 3)) c.nd = (int)v;
    if (!strncmp(argv[i], "ne=", 3)) c.ne = (int)v;
    if (!strncmp(argv[i], "th=", 3)) c.theta = v;
    if (!strncmp(argv[i], "n=", 2)) c.n = v;
    if (!strncmp(argv[i], "alpha=", 6)) c.alpha = v;
    if (!strncmp(argv[i], "dth=", 4)) c.dth = v;
    if (!strncmp(argv[i], "it=", 3)) itmax = (int)v;
    if (!strncmp(argv[i], "oracle=", 7)) c.oracle = (int)v;
    if (!strncmp(argv[i], "lev=", 4)) c.lev = (int)v;
    if (!strncmp(argv[i], "lodk=", 5)) c.lodk = v;
  }
  if (c.nd > MAXDIR || c.ne * c.ne > MAXCELL) return 1;
  fresnel(&c);
  int ncell = c.ne * c.ne;
  double L = 0.5 * (double)c.ne * c.W;
  int contrast = fabs(c.n - 1.0) > 1e-12;
  printf("tdg2d: W=%.4g lam  kW=%.4g  ND=%d  NE=%d  box=%.4g lam  n=%.3f alpha=%.2f th=%.2f\n",
         c.W / LAM, c.k0 * c.W, c.nd, c.ne, 2.0 * L / LAM, c.n, c.alpha, c.theta);

  /* --- CELLS: a graded ladder, not a grid --------------------------------
   * The architecture sizes elements by LOD, so the mesh has to carry cells of
   * different sizes side by side. In DG that costs nothing structurally: a face
   * between a coarse and a fine cell is just the segment they share, and a
   * coarse cell simply has several faces along the side a fine neighbour
   * touches. "Hanging nodes" are not a concept here.
   * Refinement rule: split while the cell is larger than lodk times its distance
   * to the interface — regulator 2 of the plan (a FLOOR on size driven by
   * geometry), capped by lev. lev = 0 gives back the uniform grid. */
  int ncellv = 0;
  {
    int head = 0;
    for (int a = 0; a < c.ne; a++)
      for (int bb = 0; bb < c.ne; bb++) {
        cellx[ncellv] = -L + ((double)a + 0.5) * c.W;
        celly[ncellv] = -L + ((double)bb + 0.5) * c.W;
        cellh[ncellv] = c.W;
        celld[ncellv] = 0;
        ncellv++;
      }
    while (head < ncellv) {
      int q = head++;
      if (celld[q] >= c.lev) continue;
      double d = fabs(sdist(&c, cellx[q], celly[q]));
      if (!(cellh[q] > c.lodk * d)) continue;
      if (ncellv + 3 >= MAXCELL) break;
      double h2 = 0.5 * cellh[q], x0 = cellx[q], y0 = celly[q];
      int dep = celld[q] + 1;
      cellh[q] = h2;
      celld[q] = dep;
      cellx[q] = x0 - 0.5 * h2;
      celly[q] = y0 - 0.5 * h2;
      const double ox[3] = {0.5, -0.5, 0.5}, oy[3] = {-0.5, 0.5, 0.5};
      for (int t = 0; t < 3; t++) {
        cellx[ncellv] = x0 + ox[t] * h2;
        celly[ncellv] = y0 + oy[t] * h2;
        cellh[ncellv] = h2;
        celld[ncellv] = dep;
        ncellv++;
      }
      head = 0; /* the split cell and its siblings must be reconsidered */
    }
  }
  ncell = ncellv;
  {
    double hmin = 1e300, hmax = 0.0;
    for (int q = 0; q < ncell; q++) {
      if (cellh[q] < hmin) hmin = cellh[q];
      if (cellh[q] > hmax) hmax = cellh[q];
    }
    printf("  cells=%d  sizes %.4g..%.4g lam  (ratio %.0f)\n", ncell, hmin / LAM, hmax / LAM,
           hmax / hmin);
  }

  /* --- sub-cells ---------------------------------------------------------- */
  int nsub = 0, dim = 0;
  {
    for (int q = 0; q < ncell; q++) {
      subof[q][0] = subof[q][1] = -1;
      double smin = 1e300, smax = -1e300;
      for (int i = 0; i < 4; i++) {
        double x = cellx[q] + ((i & 1) ? 0.5 : -0.5) * cellh[q];
        double y = celly[q] + ((i & 2) ? 0.5 : -0.5) * cellh[q];
        double s = sdist(&c, x, y);
        if (s < smin) smin = s;
        if (s > smax) smax = s;
      }
      int cut = contrast && smin < 0.0 && smax > 0.0;
      int own = (contrast && sdist(&c, cellx[q], celly[q]) < 0.0) ? 0 : 1;
      for (int sd = 0; sd < 2; sd++) {
        if (!cut && sd != own) continue;
        if (nsub >= MAXSUB) return 1;
        sub *s = &S[nsub];
        s->cell = q;
        s->side = sd;
        s->cx = cellx[q];
        s->cy = celly[q];
        double km = (sd == 0 && contrast) ? c.k0 * c.n : c.k0;
        if (c.oracle) {
          /* the field's OWN directions: transmitted inside, incident + reflected
           * outside; with no contrast the two outside waves coincide */
          if (!contrast) {
            /* the incident wave's OWN direction, which fresnel() defines against
             * the interface normal — not the raw angle theta. dth rotates it
             * deliberately, which is how the direction tolerance is measured. */
            double a0 = atan2(c.py[0], c.px[0]) + c.dth;
            s->nd = 1;
            s->kx[0] = c.k0 * cos(a0);
            s->ky[0] = c.k0 * sin(a0);
          } else if (sd == 0) {
            s->nd = 1;
            s->kx[0] = c.px[2];
            s->ky[0] = c.py[2];
          } else {
            s->nd = 2;
            for (int d = 0; d < 2; d++) {
              s->kx[d] = c.px[d];
              s->ky[d] = c.py[d];
            }
          }
        } else {
          s->nd = c.nd;
          for (int d = 0; d < c.nd; d++) {
            double th = 2.0 * M_PI * ((double)d + 0.5) / (double)c.nd;
            s->kx[d] = km * cos(th);
            s->ky[d] = km * sin(th);
          }
        }
        s->base = dim;
        dim += s->nd;
        subof[q][sd] = nsub++;
      }
    }
  }
  if (dim > MAXDIM) {
    printf("tdg2d: dim %d too large\n", dim);
    return 1;
  }
  printf("  cells=%d  sub-cells=%d  unknowns=%d\n", ncell, nsub, dim);

  size_t cap = 48u << 20;
  int *ja = malloc(cap * sizeof(int));
  double complex *va = malloc(cap * sizeof(double complex));
  size_t *rp = malloc((MAXROW + 1) * sizeof(size_t));
  double complex *rhs = calloc(MAXROW, sizeof(double complex));
  if (!ja || !va || !rp || !rhs) return 1;
  size_t nnz = 0;
  int nrow = 0;
  double complex i1 = CMPLX(0.0, 1.0);

  /* One face. sb < 0 marks a box wall, where the impedance condition and its
   * analytic drive replace the two jump conditions; wmask selects which exact
   * waves live on that side. */
#define ADDFACE(SA, SB, X0, Y0, X1, Y1, NX, NY, WMASK)                                             \
  do {                                                                                             \
    int scs[2] = {(SA), (SB)};                                                                     \
    int isint = ((SB) >= 0);                                                                       \
    for (int ti = 0; ti < 2; ti++) {                                                               \
      if (scs[ti] < 0) continue;                                                                   \
      const sub *ST = &S[scs[ti]];                                                                 \
      for (int tm = 0; tm < ST->nd; tm++) {                                                        \
        double tgx = -ST->kx[tm], tgy = -ST->ky[tm];                                               \
        double complex tph = cexp(-i1 * (tgx * ST->cx + tgy * ST->cy));                            \
        for (int cond = 0; cond < (isint ? 2 : 1); cond++) {                                       \
          if (nrow >= MAXROW) return 1;                                                            \
          rp[nrow] = nnz;                                                                          \
          for (int t2 = 0; t2 < 2; t2++) {                                                         \
            if (scs[t2] < 0) continue;                                                             \
            const sub *SU = &S[scs[t2]];                                                           \
            double sgn = isint ? (t2 == 0 ? -1.0 : 1.0) : 1.0;                                     \
            for (int d = 0; d < SU->nd; d++) {                                                     \
              double gx = SU->kx[d] + tgx, gy = SU->ky[d] + tgy;                                   \
              double complex ph = cexp(-i1 * (SU->kx[d] * SU->cx + SU->ky[d] * SU->cy));           \
              double complex Iv = seg_exp(gx, gy, (X0), (Y0), (X1), (Y1)) * ph * tph;              \
              double dn = SU->kx[d] * (NX) + SU->ky[d] * (NY);                                     \
              double complex w;                                                                    \
              if (!isint)                                                                          \
                w = (i1 * dn - i1 * c.k0) * Iv;                                                    \
              else if (cond == 0)                                                                  \
                w = sgn * Iv;                                                                      \
              else                                                                                 \
                w = sgn * i1 * dn * Iv / (i1 * c.k0);                                              \
              if (!(cabs(w) > 0.0) || nnz >= cap) continue;                                        \
              ja[nnz] = SU->base + d;                                                              \
              va[nnz++] = w;                                                                       \
            }                                                                                      \
          }                                                                                        \
          double complex dr = 0.0;                                                                 \
          if (!isint)                                                                              \
            for (int q2 = 0; q2 < 3; q2++) {                                                       \
              if (!((WMASK) & (1 << q2))) continue;                                                \
              double gx = c.px[q2] + tgx, gy = c.py[q2] + tgy;                                     \
              double complex Iv = seg_exp(gx, gy, (X0), (Y0), (X1), (Y1)) * tph;                   \
              dr += c.amp[q2] * (i1 * (c.px[q2] * (NX) + c.py[q2] * (NY)) - i1 * c.k0) * Iv;       \
            }                                                                                      \
          rhs[nrow] = dr;                                                                          \
          nrow++;                                                                                  \
        }                                                                                          \
      }                                                                                            \
    }                                                                                              \
  } while (0)

  /* INTER-CELL FACES BY GEOMETRIC ADJACENCY, not by grid indices. A face is the
   * OVERLAP of two cell sides that lie on the same line, so a coarse cell facing
   * several fine ones simply yields several faces along that side, each of the
   * fine cell's length. Emitted once per pair: only the +x and +y sides of the
   * cell behind the normal are scanned. */
  int nfint = 0, nfwall = 0;
  for (int qa = 0; qa < ncell; qa++)
    for (int ax2 = 0; ax2 < 2; ax2++) {
      double ha = 0.5 * cellh[qa];
      double coord = (ax2 == 0) ? cellx[qa] + ha : celly[qa] + ha; /* the +x / +y side */
      double lo = (ax2 == 0) ? celly[qa] - ha : cellx[qa] - ha, hi = lo + cellh[qa];
      double nx = (ax2 == 0) ? 1.0 : 0.0, ny = (ax2 == 0) ? 0.0 : 1.0;
      /* the box wall: no neighbour there by construction */
      if (fabs(coord - L) < 1e-9 * c.W) {
        for (int sd = 0; sd < 2; sd++) {
          if (!contrast && sd == 0) continue;
          double px0 = (ax2 == 0) ? coord : lo, py0 = (ax2 == 0) ? lo : coord;
          double px1 = (ax2 == 0) ? coord : hi, py1 = (ax2 == 0) ? hi : coord;
          if (contrast && !clip_side(&c, sd, &px0, &py0, &px1, &py1)) continue;
          int sa = subof[qa][sd];
          if (sa < 0) continue;
          int mask = !contrast ? 1 : (sd == 0 ? 4 : 3);
          nfwall++;
          ADDFACE(sa, -1, px0, py0, px1, py1, nx, ny, mask);
        }
        continue;
      }
      for (int qb = 0; qb < ncell; qb++) {
        if (qb == qa) continue;
        double hb = 0.5 * cellh[qb];
        double cb = (ax2 == 0) ? cellx[qb] - hb : celly[qb] - hb; /* B's -x / -y side */
        if (fabs(cb - coord) > 1e-9 * c.W) continue;
        double blo = (ax2 == 0) ? celly[qb] - hb : cellx[qb] - hb, bhi = blo + cellh[qb];
        double olo = lo > blo ? lo : blo, ohi = hi < bhi ? hi : bhi;
        if (!(ohi - olo > 1e-9 * c.W)) continue;
        for (int sd = 0; sd < 2; sd++) {
          if (!contrast && sd == 0) continue;
          double px0 = (ax2 == 0) ? coord : olo, py0 = (ax2 == 0) ? olo : coord;
          double px1 = (ax2 == 0) ? coord : ohi, py1 = (ax2 == 0) ? ohi : coord;
          if (contrast && !clip_side(&c, sd, &px0, &py0, &px1, &py1)) continue;
          int sa = subof[qa][sd], sb = subof[qb][sd];
          int mask = !contrast ? 1 : (sd == 0 ? 4 : 3);
          if (sa >= 0 && sb >= 0) {
            nfint++;
            ADDFACE(sa, sb, px0, py0, px1, py1, nx, ny, mask);
          }
        }
      }
    }
  /* the -x and -y sides that lie on the box wall */
  for (int qa = 0; qa < ncell; qa++)
    for (int ax2 = 0; ax2 < 2; ax2++) {
      double ha = 0.5 * cellh[qa];
      double coord = (ax2 == 0) ? cellx[qa] - ha : celly[qa] - ha;
      if (fabs(coord + L) > 1e-9 * c.W) continue;
      double lo = (ax2 == 0) ? celly[qa] - ha : cellx[qa] - ha, hi = lo + cellh[qa];
      double nx = (ax2 == 0) ? -1.0 : 0.0, ny = (ax2 == 0) ? 0.0 : -1.0;
      for (int sd = 0; sd < 2; sd++) {
        if (!contrast && sd == 0) continue;
        double px0 = (ax2 == 0) ? coord : lo, py0 = (ax2 == 0) ? lo : coord;
        double px1 = (ax2 == 0) ? coord : hi, py1 = (ax2 == 0) ? hi : coord;
        if (contrast && !clip_side(&c, sd, &px0, &py0, &px1, &py1)) continue;
        int sa = subof[qa][sd];
        if (sa < 0) continue;
        int mask = !contrast ? 1 : (sd == 0 ? 4 : 3);
        nfwall++;
        ADDFACE(sa, -1, px0, py0, px1, py1, nx, ny, mask);
      }
    }

  /* THE MATERIAL INTERFACE, as an ordinary face between the two sub-cells of a
   * cut cell. Same two conditions as every interior face — they ARE the physical
   * transmission conditions. */
  int nif = 0;
  if (contrast)
    for (int q = 0; q < ncell; q++) {
      if (subof[q][0] < 0 || subof[q][1] < 0) continue;
      double sc = sdist(&c, cellx[q], celly[q]);
      double ox = cellx[q] - sc * c.nhx, oy = celly[q] - sc * c.nhy;
      double t0 = -1e300, t1 = 1e300;
      const double lo[2] = {cellx[q] - 0.5 * c.W, celly[q] - 0.5 * c.W};
      const double hi[2] = {cellx[q] + 0.5 * c.W, celly[q] + 0.5 * c.W};
      const double dd[2] = {c.thx, c.thy}, oo[2] = {ox, oy};
      int empty = 0;
      for (int ax2 = 0; ax2 < 2; ax2++) {
        if (fabs(dd[ax2]) < 1e-300) {
          if (oo[ax2] < lo[ax2] || oo[ax2] > hi[ax2]) empty = 1;
          continue;
        }
        double ta = (lo[ax2] - oo[ax2]) / dd[ax2], tb = (hi[ax2] - oo[ax2]) / dd[ax2];
        if (ta > tb) {
          double sw = ta;
          ta = tb;
          tb = sw;
        }
        if (ta > t0) t0 = ta;
        if (tb < t1) t1 = tb;
      }
      if (empty || !(t0 < t1)) continue;
      double ax = ox + c.thx * t0, ay = oy + c.thy * t0;
      double bx = ox + c.thx * t1, by = oy + c.thy * t1;
      nif++;
      ADDFACE(subof[q][0], subof[q][1], ax, ay, bx, by, c.nhx, c.nhy, 0);
    }
  rp[nrow] = nnz;
  printf("  rows=%d  (%.1fx overdetermined)  nnz=%zu  interface faces=%d\n", nrow,
         (double)nrow / (double)dim, nnz, nif);

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

  /* With the field's own directions and no deliberate angular error the exact
   * solution is IN THE SPAN exactly, so the assembled rows must reproduce it. */
  if (c.oracle && !(fabs(c.dth) > 0.0)) {
    double complex *ce = calloc((size_t)dim, sizeof(double complex));
    if (ce) {
      for (int s = 0; s < nsub; s++) {
        const sub *SU = &S[s];
        for (int d = 0; d < SU->nd; d++) {
          int q2 = (contrast && SU->side == 0) ? 2 : (contrast ? d : 0);
          ce[SU->base + d] = c.amp[q2] * cexp(i1 * (SU->kx[d] * SU->cx + SU->ky[d] * SU->cy));
        }
      }
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
      free(ce);
    }
  }

  /* --- LSQR --------------------------------------------------------------- */
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
  if (!(beta > 0.0)) {
    printf("  ABORT: zero drive\n");
    return 1;
  }
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
    double th2 = sn * alpha;
    rhobar = -cs * alpha;
    double phi = cs * phibar;
    phibar = sn * phibar;
    for (int j = 0; j < dim; j++) {
      xs[j] += (phi / rho) * ww[j];
      ww[j] = vv[j] - (th2 / rho) * ww[j];
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

  /* --- error against Fresnel, sampled inside each sub-cell's own region ---- */
  double num = 0.0, den = 0.0;
  int nsamp = 0;
  for (int s = 0; s < nsub; s++) {
    const sub *SU = &S[s];
    for (int p = 0; p < 25; p++) {
      double ox = ((double)(p % 5) - 2.0) * 0.2 * c.W, oy = ((double)(p / 5) - 2.0) * 0.2 * c.W;
      double x = SU->cx + ox, y = SU->cy + oy;
      if (contrast && ((SU->side == 0) != (sdist(&c, x, y) < 0.0))) continue;
      double complex got = 0.0;
      for (int d = 0; d < SU->nd; d++)
        got += xs[SU->base + d] * cexp(i1 * (SU->kx[d] * (x - SU->cx) + SU->ky[d] * (y - SU->cy)));
      double complex ex = uexact(&c, x, y);
      num += cabs(got - ex) * cabs(got - ex);
      den += cabs(ex) * cabs(ex);
      nsamp++;
    }
  }
  printf("  FIELD ERROR = %.4e   (%d samples)\n", sqrt(num / den), nsamp);

  free(ja);
  free(va);
  free(rp);
  free(rhs);
  free(xs);
  free(uu);
  free(vv);
  free(ww);
  free(tv);
  return 0;
}
