/* M9b, scale study: how far can an element be coarsened, in units of lambda?
 *
 * carrier1d.c answered "the element-size limit comes from structure, not from
 * transport" — but only up to 16 lambda, because its reference is a finite
 * difference solve, which needs ~128 points per wavelength and therefore
 * cannot follow the domain past a few hundred lambda. The architecture needs
 * L/lambda of 1e3..1e7 (a 10 m room at 550 nm: LOD elements are sub-millimetre
 * at a few metres, i.e. ~1e3 lambda, and much larger further out).
 *
 * Two things had to change to get there, both of them from the M9 audit's own
 * rule (A4: analytic references are valid at any lambda):
 *   1. reference = the exact 1D Green's function of the uniform medium,
 *      u(x) = phihat(k~) * e^{i k~ (x - xsrc)} / (2 i k~) for x past the source
 *      support. No grid, exact at any distance.
 *   2. the medium is a LIST OF SEGMENTS, not an array of cells, so neither the
 *      scene nor the assembly grows with the domain.
 * The absorbing ramps stay in lambda units and keep fine elements — the
 * unsolved domain-termination question (PLAN open question 1) is deliberately
 * NOT what this measures; its cost is reported separately as the fine-element
 * count, which stays constant while the domain grows by seven decades. */
#include "carrier.h"
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum { MAXDIM = 4096, RAMPSEG = 32, NSAMP = 400, LEVFIX = 8 };
static const double LAM = 16.0;  /* cells per wavelength */
static const double WFINE = 2.0; /* fine element = lambda/8, at structure */
/* Carrier magnitude of the background medium. In the directional form every
 * basis function carries a SIGNED wavenumber, so a node contributes +KCAR and
 * -KCAR instead of a cos/sin pair — same span, but the two directions are now
 * separate functions and can be selected individually (PLAN question 12). */
static double KCAR = 0.0;

static const double RAMPEL = 2.0; /* absorbing layer thickness, in ELEMENTS */
/* Cap on the coarsest element, as a fraction of the domain. Not cosmetic: the
 * absorbing layer is RAMPEL elements thick, so an unconstrained ladder (whose
 * top element is a sizeable fraction of the domain) makes the absorber a
 * quarter of the world and BURIES THE SOURCE INSIDE IT. That is what produced
 * err = 1.0 (field annihilated at birth) in the first ladder runs — a bench
 * artefact, not a property of the basis. With 1/64 the absorber is ~3%. */
static const double WMAX_FRAC = 1.0 / 64.0;

static const double RCOND = 1e-10;   /* M2 value; system is exactly degenerate */
static const double RAMP_PEAK = 2.0; /* Im k2 peak in units of k0^2 (M2 spec) */

/* Gauss-Legendre 8, exact to degree 15 — used ONCE for the source spectrum
 * phihat(k~), whose integrand is quadratic times a barely-oscillating
 * exponential (|k~| ~ 0.4 over a support of 4). Not in any hot path, and
 * independent of the domain scale. */
static const double GX[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                             -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                             0.7966664774136267,  0.9602898564975363};
static const double GW[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                             0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                             0.2223810344533745, 0.1012285362903763};

/* Int over [-2,2] of phi(t) e^{-i xi t} dt */
static double complex phihat(double complex xi) {
  double complex acc = 0.0;
  for (int p = -2; p < 2; p++) {
    double m = (double)p + 0.5, h = 0.5;
    for (int g = 0; g < 8; g++) {
      double t = m + h * GX[g];
      acc += h * GW[g] * hz_phi(t) * cexp(CMPLX(0.0, -1.0) * xi * t);
    }
  }
  return acc;
}

/* Exact outgoing solution of u'' + k~^2 u = phi((x - xsrc)/ws), valid past the
 * source support. ph must be ws * phihat(kt * ws) — the source spectrum scales
 * with the source, which matters for the similarity sweep in [3]. */
static double complex uref(double x, double complex kt, double xsrc, double complex ph) {
  return ph * cexp(CMPLX(0.0, 1.0) * kt * (x - xsrc)) / (CMPLX(0.0, 2.0) * kt);
}

/* Absorbing ramps at both ends, quartic to RAMP_PEAK * k0^2, plus one interior
 * segment carrying the small global absorption.
 * RAMP LENGTH IS PASSED IN, IN CELLS, AND SIZED IN ELEMENTS BY THE CALLER —
 * not in wavelengths as v1 specified. A layer the basis cannot resolve cannot
 * absorb: with lambda-sized ramps and elements thousands of wavelengths wide
 * the boundary simply reflects, which is what the first run of this study was
 * really measuring. (Input to PLAN open question 1.) */
static int build_segments(double dom, double k0sq, double alpha, double rl, hz_medseg *segs) {
  int ns = 0;
  for (int i = 0; i < RAMPSEG; i++) {
    double f0 = (double)i / (double)RAMPSEG, f1 = (double)(i + 1) / (double)RAMPSEG;
    double t = 1.0 - 0.5 * (f0 + f1); /* 1 at the wall, 0 at the inner edge */
    double a = alpha + (RAMP_PEAK - alpha) * t * t * t * t;
    segs[ns++] = (hz_medseg){f0 * rl, f1 * rl, CMPLX(k0sq, k0sq * a)};
  }
  segs[ns++] = (hz_medseg){rl, dom - rl, CMPLX(k0sq, k0sq * alpha)};
  for (int i = 0; i < RAMPSEG; i++) {
    double f0 = (double)i / (double)RAMPSEG, f1 = (double)(i + 1) / (double)RAMPSEG;
    double t = 0.5 * (f0 + f1);
    double a = alpha + (RAMP_PEAK - alpha) * t * t * t * t;
    segs[ns++] = (hz_medseg){dom - rl + f0 * rl, dom - rl + f1 * rl, CMPLX(k0sq, k0sq * a)};
  }
  return ns;
}

/* LOD ladder: dyadic levels W_j = WFINE * 2^j, level j active where the local
 * element size rule W = eps * dist puts it, i.e. dist in [W_j/(2 eps), W_j/eps).
 * The finest level also covers everything nearer than that, the coarsest
 * everything further. This is the L = eps*R rule of the architecture with the
 * source standing in for the camera (structure is what refinement tracks).
 *
 * The first version of this file jumped straight from lambda/8 to the coarse
 * size with NOTHING in between, and measured 25-56% error that was flat across
 * six decades of scale — i.e. it was measuring the missing ladder, not the
 * basis. Levels are what the cascade is made of; skipping them is not a
 * cheaper cascade, it is no cascade. */
static int build_lod(hz_carrier *bs, double dom, double xsrc, double eps, double wfine, int cap,
                     double *wmax, int *nlev) {
  int d = 0;
  *wmax = wfine;
  *nlev = 0;
  for (int j = 0; j < 64; j++) {
    double W = wfine * pow(2.0, (double)j);
    if (W > WMAX_FRAC * dom) break;
    int top = (wfine * pow(2.0, (double)(j + 1)) > WMAX_FRAC * dom);
    double dhi = top ? dom * 2.0 : W / eps;
    double dlo = (j == 0) ? 0.0 : W / (2.0 * eps);
    /* MARGIN: phi has support 4W, so the partition of unity (sum of translates
     * == 4) only holds where FOUR translates overlap. Selecting exactly the
     * nodes whose centre falls in the band leaves the band edges with one or
     * two translates and the basis dips there — that alone cost 30% error in
     * the first run of this ladder. Widen by 2 nodes on each side so every
     * point of the band is covered by a full set of translates. */
    dlo -= 2.0 * W;
    dhi += 2.0 * W;
    if (dlo < 0.0) dlo = 0.0;
    int n0 = (int)((xsrc - dhi) / W) - 2, n1 = (int)((xsrc + dhi) / W) + 2;
    int added = 0;
    for (int n = n0; n <= n1; n++) {
      double xc = (double)n * W;
      if (xc < -2.0 * W || xc > dom + 2.0 * W) continue;
      double dist = fabs(xc - xsrc);
      if (dist < dlo || dist >= dhi) continue;
      if (d + 2 > cap) return d;
      for (int q = 0; q < 2; q++)
        bs[d++] = (hz_carrier){W, n, q ? -KCAR : KCAR, -1};
      added = 1;
    }
    if (added) {
      *wmax = W;
      (*nlev)++;
    }
  }
  return d;
}

/* Relative L2 against the analytic reference in a 4-lambda window at x0. */
static double window_err(const hz_carrier *bs, int dim, const double complex *c, double x0,
                         double complex kt, double xsrc, double complex ph) {
  double num = 0.0, den = 0.0;
  double w = 4.0 * LAM;
  for (int p = 0; p < NSAMP; p++) {
    double x = x0 + w * ((double)p + 0.5) / (double)NSAMP;
    double complex s = 0.0;
    for (int i = 0; i < dim; i++)
      s += c[i] * hz_carrier_val(bs[i], x, 0.0);
    double complex d = s - uref(x, kt, xsrc, ph);
    num += creal(d * conj(d));
    den += creal(uref(x, kt, xsrc, ph) * conj(uref(x, kt, xsrc, ph)));
  }
  return sqrt(num / den);
}

typedef struct {
  int dim, nlev, rank;
  double wmax, wfine, cond, nnzrow, smin_rel;
  long calls;
  double err[3];
} caseres;

static int run_case(double dom, double eps, double rcond, double wfine, hz_carrier *bs,
                    hz_medseg *segs, double k0, caseres *out) {
  double xsrc = 0.15 * dom;
  /* the source is half the finest element, so a similarity sweep can scale
   * every length in the problem together (phase [3]) */
  double ws = 0.5 * wfine;
  /* keep the total absorption across the domain comparable at every scale */
  double alpha = 2.0 / (k0 * dom);
  caseres r = {0};
  r.wmax = wfine;
  r.wfine = wfine;
  r.dim = build_lod(bs, dom, xsrc, eps, wfine, MAXDIM, &r.wmax, &r.nlev);
  /* the absorbing layer must be several ELEMENTS thick to be representable */
  double rl = RAMPEL * r.wmax;
  int nseg = build_segments(dom, k0 * k0, alpha, rl, segs);
  int dim = r.dim;

  double complex *A = calloc((size_t)dim * (size_t)dim, sizeof(double complex));
  double complex *b = calloc((size_t)dim, sizeof(double complex));
  double *sv = calloc((size_t)dim, sizeof(double));
  if (!A || !b || !sv) {
    free(A);
    free(b);
    free(sv);
    return -1;
  }
  long nnz = 0;
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      double lo, hi;
      hz_carrier_overlap(bs[i], bs[j], dom, 0.0, &lo, &hi);
      if (lo >= hi) continue;
      nnz++;
      A[(size_t)i * (size_t)dim + (size_t)j] =
          hz_carrier_entry(bs[i], bs[j], dom, 0.0, segs, nseg, &r.calls);
    }
    b[i] = hz_carrier_rhs(bs[i], dom, 0.0, xsrc, ws, &r.calls);
  }
  r.nnzrow = (double)nnz / (double)dim;

  lapack_int rank = 0;
  lapack_int info = LAPACKE_zgelsd(LAPACK_ROW_MAJOR, dim, dim, 1, A, dim, b, 1, sv, rcond, &rank);
  r.rank = (int)rank;
  r.cond = (rank > 0 && sv[rank - 1] > 0.0) ? sv[0] / sv[rank - 1] : -1.0;
  /* smallest singular value of the WHOLE matrix relative to the largest —
   * unlike cond it does not depend on where the truncation fell */
  r.smin_rel = (sv[0] > 0.0) ? sv[dim - 1] / sv[0] : -1.0;

  double complex kt = k0 * csqrt(CMPLX(1.0, alpha));
  double complex ph = ws * phihat(kt * ws);
  for (int w = 0; w < 3; w++)
    r.err[w] = -1.0;
  if (info == 0) {
    static const double FRAC[3] = {0.25, 0.5, 0.75};
    for (int w = 0; w < 3; w++) {
      double x0 = xsrc + FRAC[w] * (dom - rl - xsrc);
      r.err[w] = window_err(bs, dim, b, x0, kt, xsrc, ph);
    }
  }
  free(A);
  free(b);
  free(sv);
  *out = r;
  return 0;
}

int main(void) {
  double k0 = 2.0 * M_PI / LAM;
  KCAR = k0;
  printf("M9b scale study: element size in wavelengths, analytic Green reference\n");
  printf("LOD ladder: element size = eps * distance from the structure; ramp = %.0f elements\n",
         RAMPEL);
  printf("absorption scaled so the field decays ~e^-1 across the domain\n\n");
  printf("  domain(lam)  levels  Wmax/lam     dim  nnz/row     calls  rank    cond   "
         "err@25%%  @50%%  @75%%\n");

  hz_carrier *bs = calloc(MAXDIM, sizeof(hz_carrier));
  hz_medseg *segs = calloc(2 * RAMPSEG + 2, sizeof(hz_medseg));
  if (bs == NULL || segs == NULL) {
    free(bs);
    free(segs);
    return 1;
  }

  /* EPS sweep: the level bands are [W/(2 eps), W/eps), i.e. W/(2 eps) wide,
   * while phi's support is 4W. Bands narrower than the support cannot be
   * separated — the partition-of-unity margin then swallows the lower edge and
   * every level ends up active everywhere. That needs eps <= 1/8, and it is a
   * REAL CONSTRAINT ON THE LOD RULE, not a tuning knob: the potential's support
   * width bounds how fast element size may grow with distance. */
  static const double EPSV[4] = {0.25, 0.125, 0.0625, 0.03125};
  static const double DLAM[4] = {640.0, 6.4e4, 6.4e6, 6.4e8};
  for (int ie = 0; ie < 4; ie++) {
    printf("\n  --- eps = 1/%.0f (band width %.1f W vs support 4 W) ---\n", 1.0 / EPSV[ie],
           0.5 / EPSV[ie]);
    for (int e = 0; e < 4; e++) {
      caseres r;
      if (run_case(DLAM[e] * LAM, EPSV[ie], RCOND, WFINE, bs, segs, k0, &r) != 0) break;
      printf("  %11.0f  %6d  %8.0f  %6d  %7.1f  %8ld  %4d  %.1e   %6.4f %6.4f %6.4f\n", DLAM[e],
             r.nlev, r.wmax / LAM, r.dim, r.nnzrow, r.calls, r.rank, r.cond, r.err[0], r.err[1],
             r.err[2]);
    }
  }

  /* --- is the accuracy loss TRUNCATION or genuine degeneracy? --------------
   * The two have opposite cures. If some rcond recovers the field, the modes
   * were there and the threshold was wrong — a solver fix. If no rcond does,
   * the ladder is genuinely degenerate and the basis needs orthogonalisation
   * between levels. smin/smax is printed because, unlike cond, it does not
   * depend on where the truncation happened to fall. */
  printf("\n[2] rcond sweep on the configurations that lost accuracy above\n");
  static const double RC[6] = {1e-4, 1e-6, 1e-8, 1e-10, 1e-12, 1e-14};
  static const double CASE_D[3] = {6.4e4, 6.4e6, 6.4e8};
  static const double CASE_E[3] = {0.0625, 0.0625, 0.03125};
  for (int c = 0; c < 3; c++) {
    printf("\n  domain %.0e lambda, eps = 1/%.0f\n", CASE_D[c], 1.0 / CASE_E[c]);
    printf("      rcond    rank/dim   smin/smax      cond   err@25%%  @50%%  @75%%\n");
    for (int i = 0; i < 6; i++) {
      caseres r;
      if (run_case(CASE_D[c] * LAM, CASE_E[c], RC[i], WFINE, bs, segs, k0, &r) != 0) break;
      printf("     %.0e   %5d/%-5d  %.2e  %.2e   %6.4f %6.4f %6.4f\n", RC[i], r.rank, r.dim,
             r.smin_rel, r.cond, r.err[0], r.err[1], r.err[2]);
    }
  }

  /* --- WIDTH or DEPTH? -----------------------------------------------------
   * In [1] element width and ladder depth grow together, so they cannot be
   * told apart there. Here EVERY length scales by the same factor (domain,
   * finest element, coarsest element, source width), so the ladder keeps a
   * FIXED number of levels and the only thing that changes is L/lambda.
   * Flat error across the sweep => absolute width is innocent and depth is the
   * culprit. Rising error => width itself costs. */
  printf("\n[3] similarity sweep: %d levels held fixed, every length scaled together\n", LEVFIX);
  printf(
      "  domain(lam)   Wfine/lam   Wmax/lam  levels    dim      cond   err@25%%  @50%%  @75%%\n");
  for (int e = 0; e < 6; e++) {
    double dom = 640.0 * pow(10.0, (double)e) * LAM;
    /* choose the finest element so that exactly LEVFIX levels fit under the cap */
    double wfine = WMAX_FRAC * dom / pow(2.0, (double)(LEVFIX - 1));
    caseres r;
    if (run_case(dom, 0.0625, 1e-8, wfine, bs, segs, k0, &r) != 0) break;
    printf("  %11.0f  %10.1f  %9.0f  %6d %6d  %.2e   %6.4f %6.4f %6.4f\n", dom / LAM, r.wfine / LAM,
           r.wmax / LAM, r.nlev, r.dim, r.cond, r.err[0], r.err[1], r.err[2]);
  }

  free(bs);
  free(segs);
  return 0;
}
