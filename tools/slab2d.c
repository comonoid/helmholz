/* THE BENCH AT THE ELEMENT SIZE THE PROJECT ACTUALLY TARGETS.
 *
 * WHY A NEW BENCH. Everything measured so far ran at W = lambda/2 on an object
 * one wavelength across, because that is what the Mie reference and the DtN
 * circle could afford. But the target regime is W >~ 1000 lambda on objects of
 * 10^4..10^6 lambda (a table, a wall, a chair), and the two are not the same
 * problem: the pathologies fought at W = lambda/2 (a kernel carrying field,
 * 10^4 LSQR iterations, cond ~ 1e4) all come from NEARLY PARALLEL directions,
 * and two directions are distinguishable to the operator only when
 * k W dtheta >~ 1. At W = lambda/2 with 8 directions that product is 0.8 — the
 * fan is invisible to the operator. At W = 1000 lambda it is 1600. So the
 * measurement had to move, not the method.
 *
 * WHAT IS TESTED HERE: the assembled SYSTEM (operator + cut + Nitsche + solver)
 * at k W up to 1e5, against a closed-form reference.
 * WHAT IS NOT: domain termination. The data on the box walls come from the exact
 * solution, so this says nothing about a radiation condition — and that is
 * deliberate, because the exact DtN on a circle has rank ~2kR and is hopeless at
 * room scale (26 000 harmonics at R = 2000 lambda, 1e8 in a room). Termination
 * needs a different mechanism and its own bench.
 *
 * GEOMETRY: one straight material interface through a box of NE x NE elements —
 * a table top, a wall. The exact solution is Fresnel: incident + reflected
 * outside, transmitted inside, all closed form, valid at any kW.
 *
 * BOUNDARY: IMPEDANCE, not Dirichlet. A Dirichlet box at kL ~ 1e4 sits in a
 * dense spectrum of interior resonances and would be measuring the resonance
 * rather than the method; the impedance problem has none. It also needs no
 * penalty parameter at all: with dn u = i k u + h the boundary term drops
 * straight out of Green's identity.
 *
 *   row_i:  Strong_ij - Int_dbox B_i (dn B_j - i k B_j)  =  -Int_dbox B_i h
 *
 * and h = dn u_exact - i k u_exact is known analytically. */
#include "carrier2d.h"
#include "cut2d.h"
#include "nitsche2d.h"
#include "phi.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { MAXDIM = 40000 };
static const double LAM = 1.0;
/* L2 inverse inequality for a quadratic on a piece of length W (M15 audit) */
static const double CINV = 7.745966692414834;

typedef struct {
  int nx, ny;
  int side; /* -1 uncut, 0 inside (medium n), 1 outside */
  double complex kx, ky;
  double km2, kc2;
} elem;

typedef struct {
  double W, k0, n, alpha, theta;
  int nd, ne;
  double pen, beta;
  int oracle;                  /* 1 = put the TRUE Fresnel directions in the fan */
  int nocut;                   /* 1 = suppress the cut entirely (isolation of the other blocks) */
  double L;                    /* half-side of the box */
  double nhx, nhy;             /* interface normal, pointing OUT of the medium */
  double thx, thy;             /* interface tangent */
  double complex amp[3];       /* 1, r, t */
  double complex px[3], py[3]; /* wave vectors of incident, reflected, transmitted */
} cfg;

/* --- exact Fresnel solution ---------------------------------------------- */
static void fresnel(cfg *c) {
  c->nhx = cos(c->alpha);
  c->nhy = sin(c->alpha);
  c->thx = -sin(c->alpha);
  c->thy = cos(c->alpha);
  double k1 = c->k0, k2 = c->k0 * c->n;
  double kt = k1 * sin(c->theta);
  double c1 = sqrt(k1 * k1 - kt * kt), c2 = sqrt(k2 * k2 - kt * kt);
  double complex r = (c1 - c2) / (c1 + c2), t = 2.0 * c1 / (c1 + c2);
  c->amp[0] = 1.0;
  c->amp[1] = r;
  c->amp[2] = t;
  /* incident travels toward the medium (-n), reflected away (+n), transmitted
   * into it (-n) with the interior normal component */
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

/* --- basis ---------------------------------------------------------------- */
static int straddles(const cfg *c, int nx, int ny) {
  if (c->nocut) return 0; /* isolation: no cut at all, leaving only volume + impedance */
  double d = fabs(sdist(c, c->W * (double)nx, c->W * (double)ny));
  return d < 2.0 * c->W * sqrt(2.0);
}

static int build(const cfg *c, elem *b) {
  int nmax = (int)(c->L / c->W) + 3, dim = 0;
  for (int nx = -nmax; nx <= nmax; nx++)
    for (int ny = -nmax; ny <= nmax; ny++) {
      double px = c->W * (double)nx, py = c->W * (double)ny;
      /* the support must meet the box */
      if (fabs(px) > c->L + 2.0 * c->W || fabs(py) > c->L + 2.0 * c->W) continue;
      int str = straddles(c, nx, ny);
      for (int pass = 0; pass < 2; pass++) {
        int side, inside;
        if (str) {
          side = pass;
          inside = (pass == 0);
        } else {
          if (pass == 1) break;
          side = -1;
          inside = sdist(c, px, py) < 0.0;
        }
        double km = inside ? c->k0 * c->n : c->k0;
        /* ORACLE: the true directions instead of a uniform fan. At kW >> 1 a
         * carrier direction must be right to dtheta ~ 1/(kW) — 1.6e-4 rad at
         * W = 1000 lam — and no uniform fan of a few directions can be. This
         * separates "the method fails at large W" from "a uniform fan fails at
         * large W", which are entirely different statements. */
        int nds = c->oracle ? (inside ? 1 : 2) : c->nd;
        for (int d = 0; d < nds; d++) {
          /* HALF-STEP OFFSET ON PURPOSE: if a true Fresnel direction were in the
           * fan the exact solution would be reproduced exactly and the bench
           * would measure nothing (M14 artefact 16 — breaking both sides of a
           * comparison at once). The fan must be a genuine approximation. */
          double th = 2.0 * M_PI * ((double)d + 0.5) / (double)c->nd;
          double dkx = km * cos(th), dky = km * sin(th);
          if (c->oracle) {
            int q = inside ? 2 : d;
            dkx = creal(c->px[q]);
            dky = creal(c->py[q]);
          }
          if (dim >= MAXDIM) return dim;
          b[dim].nx = nx;
          b[dim].ny = ny;
          b[dim].side = side;
          b[dim].kx = dkx;
          b[dim].ky = dky;
          b[dim].km2 = km * km;
          b[dim].kc2 = km * km;
          dim++;
        }
      }
    }
  return dim;
}

/* --- assembly ------------------------------------------------------------- */
static double pair_k2(elem ei, elem ej) {
  return ei.side >= 0 ? ei.km2 : ej.km2;
}

/* Does a half-plane actually cut this rectangle? Deep inside the box neither the
 * walls nor (for an uncut pair) the interface do, and then the region is the
 * plain support rectangle and the integral factorises. Keeping every pair on the
 * polygon path costs the whole assembly for nothing. */
/* THREE cases, and collapsing them into two was the bug. Signed distance of the
 * rectangle's corners to the plane, keeping {s <= 0}:
 *   lo >= 0  the rectangle is wholly OUTSIDE — the region is EMPTY;
 *   hi >  0  the plane genuinely cuts — it must be carried;
 *   else     the rectangle is wholly inside — the plane is a no-op and dropping
 *            it is what buys the separable fast path.
 * The first and second used to be tested by one predicate, so a support lying
 * beyond a wall and merely TOUCHING it (lo = 0, hi > 0 — no strict crossing) was
 * read as "no cut" and integrated over its whole rectangle, which lies outside
 * the domain entirely. Elements have nodes up to 2W beyond the box, so this is
 * not a corner case, it is a whole ring of them. */
enum { PLANE_EMPTY = -1, PLANE_SKIP = 0, PLANE_CUTS = 1 };
static int classify(hz_half2 h, double xlo, double xhi, double ylo, double yhi) {
  double lo = 1e300, hi = -1e300;
  const double xs[2] = {xlo, xhi}, ys[2] = {ylo, yhi};
  for (int a = 0; a < 2; a++)
    for (int d = 0; d < 2; d++) {
      double s = xs[a] * h.ca + ys[d] * h.sa - h.c;
      if (s < lo) lo = s;
      if (s > hi) hi = s;
    }
  if (lo >= 0.0) return PLANE_EMPTY;
  return hi > 0.0 ? PLANE_CUTS : PLANE_SKIP;
}

/* Returns the number of half-planes that matter, filling hp; -1 if the region is
 * empty. */
static int region_planes(const cfg *c, const hz_half2 *box, int cut, double xlo, double xhi,
                         double ylo, double yhi, hz_half2 *hp, int *nbox) {
  int n = 0;
  for (int i = 0; i < 4; i++) {
    int cl = classify(box[i], xlo, xhi, ylo, yhi);
    if (cl == PLANE_EMPTY) return -1;
    if (cl == PLANE_CUTS) hp[n++] = box[i];
  }
  *nbox = n;
  if (cut) {
    hp[n].ca = c->nhx;
    hp[n].sa = c->nhy;
    hp[n].c = 0.0;
    hp[n].keep_le = 1; /* the medium is s <= 0 */
    n++;
  }
  return n;
}

static double complex region_int(const hz_half2 *hp, int nbox, int nall, int sidei, int sidej,
                                 double xlo, double xhi, double ylo, double yhi, hz_axis2 fx,
                                 hz_axis2 fy, double complex omx, double complex omy) {
  int in = (sidei == 0) || (sidej == 0);
  int out = (sidei == 1) || (sidej == 1);
  if (in && out) return 0.0;
  if (in) return hz_cut2d_poly(xlo, xhi, ylo, yhi, fx, fy, omx, omy, hp, nall);
  if (out)
    return hz_cut2d_poly(xlo, xhi, ylo, yhi, fx, fy, omx, omy, hp, nbox) -
           hz_cut2d_poly(xlo, xhi, ylo, yhi, fx, fy, omx, omy, hp, nall);
  return hz_cut2d_poly(xlo, xhi, ylo, yhi, fx, fy, omx, omy, hp, nbox);
}

/* Separable fast path: nothing cuts the support rectangle, so the region is the
 * rectangle itself and the whole entry is six 1D closed forms. */
static double complex entry_sep(double W, elem ei, elem ej, double xlo, double xhi, double ylo,
                                double yhi) {
  double complex omx = ei.kx + ej.kx, omy = ei.ky + ej.ky, i1 = CMPLX(0.0, 1.0);
  hz_phi_factor fi = {W, (double)ei.nx, 0}, gi = {W, (double)ei.ny, 0};
  hz_phi_factor f0 = {W, (double)ej.nx, 0}, f1 = {W, (double)ej.nx, 1}, f2 = {W, (double)ej.nx, 2};
  hz_phi_factor g0 = {W, (double)ej.ny, 0}, g1 = {W, (double)ej.ny, 1}, g2 = {W, (double)ej.ny, 2};
  double complex X0 = hz_phi_prod_integral_osc(xlo, xhi, fi, f0, omx);
  double complex X1 = hz_phi_prod_integral_osc(xlo, xhi, fi, f1, omx);
  double complex X2 = hz_phi_prod_integral_osc(xlo, xhi, fi, f2, omx);
  double complex Y0 = hz_phi_prod_integral_osc(ylo, yhi, gi, g0, omy);
  double complex Y1 = hz_phi_prod_integral_osc(ylo, yhi, gi, g1, omy);
  double complex Y2 = hz_phi_prod_integral_osc(ylo, yhi, gi, g2, omy);
  double complex pref = cexp(-i1 * (ei.kx * W * (double)ei.nx + ej.kx * W * (double)ej.nx +
                                    ei.ky * W * (double)ei.ny + ej.ky * W * (double)ej.ny));
  double dk2 = pair_k2(ei, ej) - ej.kc2;
  return pref *
         ((X2 + 2.0 * i1 * ej.kx * X1) * Y0 + X0 * (Y2 + 2.0 * i1 * ej.ky * Y1) + dk2 * X0 * Y0);
}

static double trace_alpha(int side) {
  return side < 0 ? 0.0 : (side == 1 ? 1.0 : -1.0);
}
static double trace_beta(int side) {
  return side < 0 ? 1.0 : 0.5;
}

/* Clip a segment to one side of the interface. side 0 = the medium (s < 0),
 * side 1 = outside. Returns 0 when nothing of the segment survives. */
static int clip_side(const cfg *c, int side, double *ax, double *ay, double *bx, double *by) {
  double s0 = sdist(c, *ax, *ay), s1 = sdist(c, *bx, *by);
  int k0in = (side == 0) ? (s0 < 0.0) : (s0 > 0.0);
  int k1in = (side == 0) ? (s1 < 0.0) : (s1 > 0.0);
  if (!k0in && !k1in) return 0;
  if (k0in && k1in) return 1;
  double t = s0 / (s0 - s1);
  double mx = *ax + t * (*bx - *ax), my = *ay + t * (*by - *ay);
  if (k0in) {
    *bx = mx;
    *by = my;
  } else {
    *ax = mx;
    *ay = my;
  }
  return 1;
}

/* the interface chord clipped to the box, as a segment */
static int iface_seg(const cfg *c, double *x0, double *y0, double *x1, double *y1) {
  double t0 = -1e300, t1 = 1e300;
  /* point on the line: origin; direction: the tangent */
  double dx = c->thx, dy = c->thy;
  const double lim[2] = {-c->L, c->L};
  for (int ax = 0; ax < 2; ax++) {
    double d = ax ? dy : dx, p = 0.0;
    for (int s = 0; s < 2; s++) {
      if (fabs(d) < 1e-300) {
        if (p < lim[0] || p > lim[1]) return 0;
        continue;
      }
      double tt = (lim[s] - p) / d;
      if (d > 0.0) {
        if (s == 0 && tt > t0) t0 = tt;
        if (s == 1 && tt < t1) t1 = tt;
      } else {
        if (s == 0 && tt < t1) t1 = tt;
        if (s == 1 && tt > t0) t0 = tt;
      }
    }
  }
  if (!(t0 < t1)) return 0;
  *x0 = dx * t0;
  *y0 = dy * t0;
  *x1 = dx * t1;
  *y1 = dy * t1;
  return 1;
}

int main(int argc, char **argv) {
  static const char *const KEYS[] = {"W",   "nd", "ne",     "n",     "alpha", "th",
                                     "pen", "it", "oracle", "nocut", NULL};
  for (int i = 1; i < argc; i++) {
    const char *eq = strchr(argv[i], '=');
    int ok = 0;
    if (eq)
      for (int k = 0; KEYS[k]; k++)
        if ((size_t)(eq - argv[i]) == strlen(KEYS[k]) &&
            !strncmp(argv[i], KEYS[k], strlen(KEYS[k])))
          ok = 1;
    if (!ok) {
      printf("slab2d: unknown argument '%s'; keys are", argv[i]);
      for (int k = 0; KEYS[k]; k++)
        printf(" %s", KEYS[k]);
      printf("\n");
      return 1;
    }
  }
  cfg c;
  memset(&c, 0, sizeof c);
  c.k0 = 2.0 * M_PI / LAM;
  c.W = 1.0 * LAM;
  c.nd = 8;
  c.ne = 8;
  c.n = 1.5;
  c.alpha = 0.4;
  c.theta = 0.3;
  c.pen = 1.0;
  int itmax = 4000;
  int oracle = 0;
  for (int i = 1; i < argc; i++) {
    const char *eq = strchr(argv[i], '=');
    double v = atof(eq + 1);
    if (!strncmp(argv[i], "W=", 2)) c.W = v * LAM;
    if (!strncmp(argv[i], "nd=", 3)) c.nd = (int)v;
    if (!strncmp(argv[i], "ne=", 3)) c.ne = (int)v;
    if (!strncmp(argv[i], "n=", 2)) c.n = v;
    if (!strncmp(argv[i], "alpha=", 6)) c.alpha = v;
    if (!strncmp(argv[i], "th=", 3)) c.theta = v;
    if (!strncmp(argv[i], "pen=", 4)) c.pen = v;
    if (!strncmp(argv[i], "it=", 3)) itmax = (int)v;
    if (!strncmp(argv[i], "oracle=", 7)) oracle = (int)v;
    if (!strncmp(argv[i], "nocut=", 6)) c.nocut = (int)v;
  }
  c.oracle = oracle;
  c.L = 0.5 * (double)c.ne * c.W;
  fresnel(&c);
  double kmax = c.k0 * (c.n > 1.0 ? c.n : 1.0);
  c.beta = c.pen * (kmax + CINV / c.W);
  printf("slab2d: W=%.4g lam  kW=%.4g  ND=%d  NE=%d  box=%.4g lam  n=%.2f  alpha=%.2f th=%.2f\n",
         c.W / LAM, c.k0 * c.W, c.nd, c.ne, 2.0 * c.L / LAM, c.n, c.alpha, c.theta);

  static elem b[MAXDIM];
  int dim = build(&c, b);
  printf("  dim=%d  (elements x directions; INDEPENDENT of W by construction)\n", dim);
  if (dim >= MAXDIM) {
    printf("  ABORT: basis truncated\n");
    return 1;
  }

  hz_half2 box[4];
  const double bn[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
  for (int i = 0; i < 4; i++) {
    box[i].ca = bn[i][0];
    box[i].sa = bn[i][1];
    box[i].c = c.L;
    box[i].keep_le = 1;
  }
  /* wall segments, counter-clockwise, with outward normals */
  double wx0[4] = {-c.L, c.L, c.L, -c.L}, wy0[4] = {-c.L, -c.L, c.L, c.L};
  double wx1[4] = {c.L, c.L, -c.L, -c.L}, wy1[4] = {-c.L, c.L, c.L, -c.L};
  double ix0, iy0, ix1, iy1;
  int have_iface = iface_seg(&c, &ix0, &iy0, &ix1, &iy1);

  size_t cap = 16u << 20;
  int *ja = malloc(cap * sizeof(int));
  double complex *va = malloc(cap * sizeof(double complex));
  size_t *rp = calloc((size_t)dim + 1, sizeof(size_t));
  double complex *rhs = calloc((size_t)dim, sizeof(double complex));
  if (!ja || !va || !rp || !rhs) return 1;
  size_t nnz = 0;
  double complex i1 = CMPLX(0.0, 1.0);

  /* The exact solution's coefficients, built BEFORE assembly so each block can be
   * tested against it as it is formed. With oracle directions the partition of
   * unity holds across the whole box (elements out to L+2W are kept), so this
   * vector represents the exact field with no approximation at all. */
  double complex *ce = NULL, *dvol = NULL, *dnit = NULL, *dbnd = NULL;
  if (oracle) {
    ce = calloc((size_t)dim, sizeof(double complex));
    dvol = calloc((size_t)dim, sizeof(double complex));
    dnit = calloc((size_t)dim, sizeof(double complex));
    dbnd = calloc((size_t)dim, sizeof(double complex));
    if (!ce || !dvol || !dnit || !dbnd) exit(1);
    for (int j = 0; j < dim; j++) {
      int q = (b[j].side == 0 || (b[j].side < 0 && sdist(&c, c.W * b[j].nx, c.W * b[j].ny) < 0.0))
                  ? 2
                  : -1;
      if (q < 0) /* outside: two waves, match this column to its own direction */
        for (int t = 0; t < 2; t++)
          if (cabs(b[j].kx - creal(c.px[t])) + cabs(b[j].ky - creal(c.py[t])) < 1e-9 * c.k0) q = t;
      if (q < 0) continue;
      double xn = c.W * (double)b[j].nx, yn = c.W * (double)b[j].ny;
      ce[j] = c.amp[q] * cexp(i1 * (c.px[q] * xn + c.py[q] * yn)) / 16.0;
    }
  }

  for (int i = 0; i < dim; i++) {
    rp[i] = nnz;
    for (int j = 0; j < dim; j++) {
      if (abs(b[i].nx - b[j].nx) > 3 || abs(b[i].ny - b[j].ny) > 3) continue;
      double W = c.W;
      double xlo = W * ((double)(b[i].nx > b[j].nx ? b[i].nx : b[j].nx) - 2.0);
      double xhi = W * ((double)(b[i].nx < b[j].nx ? b[i].nx : b[j].nx) + 2.0);
      double ylo = W * ((double)(b[i].ny > b[j].ny ? b[i].ny : b[j].ny) - 2.0);
      double yhi = W * ((double)(b[i].ny < b[j].ny ? b[i].ny : b[j].ny) + 2.0);
      if (!(xlo < xhi) || !(ylo < yhi)) continue;
      double complex omx = b[i].kx + b[j].kx, omy = b[i].ky + b[j].ky;
      double complex pref =
          cexp(-i1 * (b[i].kx * W * (double)b[i].nx + b[j].kx * W * (double)b[j].nx +
                      b[i].ky * W * (double)b[i].ny + b[j].ky * W * (double)b[j].ny));
      hz_phi_factor fi0 = {W, (double)b[i].nx, 0}, fj0 = {W, (double)b[j].nx, 0};
      hz_phi_factor gi0 = {W, (double)b[i].ny, 0}, gj0 = {W, (double)b[j].ny, 0};
      hz_phi_factor fj1 = {W, (double)b[j].nx, 1}, fj2 = {W, (double)b[j].nx, 2};
      hz_phi_factor gj1 = {W, (double)b[j].ny, 1}, gj2 = {W, (double)b[j].ny, 2};
      hz_axis2 X0 = {{fi0, fj0}, 2}, X1 = {{fi0, fj1}, 2}, X2 = {{fi0, fj2}, 2};
      hz_axis2 Y0 = {{gi0, gj0}, 2}, Y1 = {{gi0, gj1}, 2}, Y2 = {{gi0, gj2}, 2};
      int cutpair = (b[i].side >= 0 || b[j].side >= 0);
      hz_half2 hp[5];
      int nbox = 0;
      int nall = region_planes(&c, box, cutpair, xlo, xhi, ylo, yhi, hp, &nbox);
      if (nall < 0) continue; /* support lies wholly outside the box */
      double dk2 = pair_k2(b[i], b[j]) - b[j].kc2;
      double complex vvol;
      if (nall == 0) {
        vvol = entry_sep(W, b[i], b[j], xlo, xhi, ylo, yhi);
      } else {
        double complex t = 0.0;
        t += region_int(hp, nbox, nall, b[i].side, b[j].side, xlo, xhi, ylo, yhi, X2, Y0, omx, omy);
        t += 2.0 * i1 * b[j].kx *
             region_int(hp, nbox, nall, b[i].side, b[j].side, xlo, xhi, ylo, yhi, X1, Y0, omx, omy);
        t += region_int(hp, nbox, nall, b[i].side, b[j].side, xlo, xhi, ylo, yhi, X0, Y2, omx, omy);
        t += 2.0 * i1 * b[j].ky *
             region_int(hp, nbox, nall, b[i].side, b[j].side, xlo, xhi, ylo, yhi, X0, Y1, omx, omy);
        if (fabs(dk2) > 0.0)
          t += dk2 * region_int(hp, nbox, nall, b[i].side, b[j].side, xlo, xhi, ylo, yhi, X0, Y0,
                                omx, omy);
        vvol = pref * t;
      }

      hz_carrier2d ci = {W, b[i].nx, b[i].ny, b[i].kx, b[i].ky};
      hz_carrier2d cj = {W, b[j].nx, b[j].ny, b[j].kx, b[j].ky};
      /* Nitsche on the material interface (M15) */
      double complex vnit = 0.0;
      if (have_iface && cutpair) {
        hz_nit2d s = hz_nitsche2d_seg(ci, cj, ix0, iy0, ix1, iy1, c.nhx, c.nhy);
        double ai = trace_alpha(b[i].side), aj = trace_alpha(b[j].side), bi = trace_beta(b[i].side);
        vnit = (bi * aj) * (s.t1j - s.t1i) - c.beta * (ai * aj) * s.t0;
      }
      /* Impedance on the box walls: -Int B_i (dn B_j - i k B_j).
       * A CUT function is zero on the far side of the interface, and the line
       * integrator knows nothing about sides, so the wall has to be clipped to
       * the side the pair shares — otherwise the two halves of a cut element
       * each contribute their whole trace and the wall term is counted twice. */
      double complex vbnd = 0.0;
      int sidepair = (b[i].side >= 0) ? b[i].side : b[j].side;
      if (!(b[i].side >= 0 && b[j].side >= 0 && b[i].side != b[j].side)) {
        for (int w = 0; w < 4; w++) {
          double ax = wx0[w], ay = wy0[w], bx = wx1[w], by = wy1[w];
          if (sidepair >= 0 && !clip_side(&c, sidepair, &ax, &ay, &bx, &by)) continue;
          double dx = wx1[w] - wx0[w], dy = wy1[w] - wy0[w];
          double len = sqrt(dx * dx + dy * dy);
          hz_nit2d s = hz_nitsche2d_seg(ci, cj, ax, ay, bx, by, dy / len, -dx / len);
          vbnd -= s.t1j - i1 * c.k0 * s.t0;
        }
      }
      /* EACH BLOCK KEPT APART FOR THE DIAGNOSTIC. The exact solution must be
       * annihilated by the volume block on its own (it satisfies the PDE
       * pointwise), by the Nitsche block on its own (all its jumps vanish), and
       * the boundary block on its own must reproduce the drive. One number per
       * block instead of one number for the sum. */
      if (dvol) {
        dvol[i] += vvol * ce[j];
        dnit[i] += vnit * ce[j];
        dbnd[i] += vbnd * ce[j];
      }
      double complex v = vvol + vnit + vbnd;
      if (!(cabs(v) > 0.0) || nnz >= cap) continue;
      ja[nnz] = j;
      va[nnz++] = v;
    }
    /* drive: -Int_dbox B_i h,  h = dn u_exact - i k u_exact */
    double complex acc = 0.0;
    hz_carrier2d ci = {c.W, b[i].nx, b[i].ny, b[i].kx, b[i].ky};
    for (int w = 0; w < 4; w++) {
      double dx = wx1[w] - wx0[w], dy = wy1[w] - wy0[w];
      double len = sqrt(dx * dx + dy * dy);
      double nx = dy / len, ny = -dx / len;
      /* the wall lies wholly outside the medium only if the interface misses it;
       * in general each wall may straddle, so both branches contribute on their
       * own side. The exact field is a plane wave on each side, so the integral
       * splits at the interface crossing. */
      for (int q = 0; q < 3; q++) {
        /* TWO clips, and the second is the one that was missing. The datum h is
         * a different plane wave on each side, so the wall splits at the
         * interface — that is the first clip. But the TEST function may itself
         * be cut, and then it is identically zero on the far side, so the wall
         * splits again by ITS side — and without that a cut test function was
         * being driven by the wave living on the other side of the interface. */
        double ax = wx0[w], ay = wy0[w], bx = wx1[w], by = wy1[w];
        if (!clip_side(&c, (q != 2) ? 1 : 0, &ax, &ay, &bx, &by)) continue;
        if (b[i].side >= 0 && !clip_side(&c, b[i].side, &ax, &ay, &bx, &by)) continue;
        double complex i0v, i1v;
        hz_nitsche2d_seg_pw(ci, c.px[q], c.py[q], ax, ay, bx, by, nx, ny, &i0v, &i1v);
        (void)i1v;
        double complex pn = c.px[q] * nx + c.py[q] * ny;
        acc += c.amp[q] * (i1 * pn - i1 * c.k0) * i0v;
      }
    }
    rhs[i] = -acc;
  }
  rp[dim] = nnz;
  printf("  nnz=%zu (%.0f per row)\n", nnz, (double)nnz / (double)dim);

  /* DIAGNOSTIC BEFORE ANY SOLVE: with oracle directions the exact solution is in
   * the span EXACTLY (partition of unity holds across the whole box, since
   * elements out to L+2W are kept), so the assembled system must reproduce it.
   * Whichever block does not is the defective one — one run instead of a hunt. */
  if (oracle) {
    /* Is the COEFFICIENT VECTOR right? A separate question from whether the
     * matrix is, and conflating the two is how a hunt starts. */
    double e = 0.0, dn = 0.0;
    for (int a = 0; a < 21; a++)
      for (int d2 = 0; d2 < 21; d2++) {
        /* the WHOLE box, not the central half: the volume rows integrate over all
         * of it, so a partition of unity that fails only near the wall would be
         * invisible to a central-half check and would show up as a volume block
         * that refuses to annihilate */
        double x = -0.999 * c.L + 1.998 * c.L * (double)a / 20.0;
        double y = -0.999 * c.L + 1.998 * c.L * (double)d2 / 20.0;
        int in2 = sdist(&c, x, y) < 0.0;
        double complex g = 0.0;
        for (int j = 0; j < dim; j++) {
          if (b[j].side == 0 && !in2) continue;
          if (b[j].side == 1 && in2) continue;
          hz_carrier2d cb = {c.W, b[j].nx, b[j].ny, b[j].kx, b[j].ky};
          g += ce[j] * hz_carrier2d_val(cb, x, y);
        }
        double complex ex2 = uexact(&c, x, y);
        e += cabs(g - ex2) * cabs(g - ex2);
        dn += cabs(ex2) * cabs(ex2);
      }
    printf("  [diag] exact COEFFICIENT VECTOR reproduces the field to %.3e\n", sqrt(e / dn));
    double v2 = 0.0, n2 = 0.0, b2 = 0.0, nb = 0.0, tot = 0.0;
    for (int r = 0; r < dim; r++) {
      v2 += cabs(dvol[r]) * cabs(dvol[r]);
      n2 += cabs(dnit[r]) * cabs(dnit[r]);
      b2 += cabs(dbnd[r] - rhs[r]) * cabs(dbnd[r] - rhs[r]);
      nb += cabs(rhs[r]) * cabs(rhs[r]);
      double complex t = dvol[r] + dnit[r] + dbnd[r] - rhs[r];
      tot += cabs(t) * cabs(t);
    }
    printf("  [diag] BLOCK BY BLOCK on the exact solution, against |b| = %.3e:\n", sqrt(nb));
    printf("         volume   |Av| = %.3e   (must annihilate: it satisfies the PDE)\n", sqrt(v2));
    printf("         Nitsche  |An| = %.3e   (must annihilate: every jump is zero)\n", sqrt(n2));
    printf("         boundary |Ab-b| = %.3e (must reproduce the drive)\n", sqrt(b2));
    printf("         total    |Ac-b| = %.3e\n", sqrt(tot));
  }
  free(ce);
  free(dvol);
  free(dnit);
  free(dbnd);

  /* row equilibration, once */
  for (int r = 0; r < dim; r++) {
    double s = 0.0;
    for (size_t p = rp[r]; p < rp[r + 1]; p++)
      s += cabs(va[p]) * cabs(va[p]);
    s = sqrt(s);
    if (!(s > 0.0)) continue;
    for (size_t p = rp[r]; p < rp[r + 1]; p++)
      va[p] /= s;
    rhs[r] /= s;
  }

  /* --- LSQR ------------------------------------------------------------- */
  double complex *xs = calloc((size_t)dim, sizeof(double complex));
  double complex *uu = calloc((size_t)dim, sizeof(double complex));
  double complex *vv = calloc((size_t)dim, sizeof(double complex));
  double complex *ww = calloc((size_t)dim, sizeof(double complex));
  double complex *tv = calloc((size_t)dim, sizeof(double complex));
  if (!xs || !uu || !vv || !ww || !tv) return 1;
  for (int r = 0; r < dim; r++)
    uu[r] = rhs[r];
  double beta = 0.0;
  for (int r = 0; r < dim; r++)
    beta += cabs(uu[r]) * cabs(uu[r]);
  beta = sqrt(beta);
  if (!(beta > 0.0)) {
    printf("  ABORT: zero right-hand side\n");
    return 1;
  }
  for (int r = 0; r < dim; r++)
    uu[r] /= beta;
  for (int j = 0; j < dim; j++)
    vv[j] = 0.0;
  for (int r = 0; r < dim; r++)
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
  int n2 = 0, n3 = 0, n4 = 0, n6 = 0;
  for (int it = 0; it < itmax; it++) {
    for (int r = 0; r < dim; r++) {
      double complex s = 0.0;
      for (size_t p = rp[r]; p < rp[r + 1]; p++)
        s += va[p] * vv[ja[p]];
      uu[r] = s - alpha * uu[r];
    }
    beta = 0.0;
    for (int r = 0; r < dim; r++)
      beta += cabs(uu[r]) * cabs(uu[r]);
    beta = sqrt(beta);
    if (!(beta > 0.0)) break;
    for (int r = 0; r < dim; r++)
      uu[r] /= beta;
    for (int j = 0; j < dim; j++)
      tv[j] = 0.0;
    for (int r = 0; r < dim; r++)
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
  }
  printf("  LSQR: |r|/|b| = %.3e   iters to 1e-2/1e-3/1e-4/1e-6: %d/%d/%d/%d\n", phibar / b0, n2,
         n3, n4, n6);

  /* --- error against Fresnel, on the CENTRAL HALF of the box -------------- */
  double num = 0.0, den = 0.0;
  int NS = 61;
  for (int a = 0; a < NS; a++)
    for (int d = 0; d < NS; d++) {
      double x = -0.5 * c.L + c.L * (double)a / (double)(NS - 1);
      double y = -0.5 * c.L + c.L * (double)d / (double)(NS - 1);
      int in = sdist(&c, x, y) < 0.0;
      double complex got = 0.0;
      for (int j = 0; j < dim; j++) {
        if (b[j].side == 0 && !in) continue;
        if (b[j].side == 1 && in) continue;
        hz_carrier2d cb = {c.W, b[j].nx, b[j].ny, b[j].kx, b[j].ky};
        got += xs[j] * hz_carrier2d_val(cb, x, y);
      }
      double complex ex = uexact(&c, x, y);
      num += cabs(got - ex) * cabs(got - ex);
      den += cabs(ex) * cabs(ex);
    }
  printf("  FIELD ERROR (central half) = %.4e\n", sqrt(num / den));

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
