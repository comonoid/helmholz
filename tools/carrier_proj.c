/* M9b/proj — CAN the carrier basis express the field at all? Representation
 * only: no solver, no domain truncation, no cascade.
 *
 * Why this exists. Four benches in a row failed at their boundary or their
 * ladder, and each time the failure was reported as a property of the method
 * before the bench was trustworthy. Best-approximation error has none of those
 * moving parts: project a KNOWN exact field onto the basis and measure the
 * residue. It is a LOWER BOUND on what any Galerkin solve can reach — a failure
 * here is decisive, a success here proves nothing about the solve (Helmholtz is
 * indefinite; quasi-optimality is not free). That asymmetry is the point.
 *
 * Design follows the audit of this experiment:
 *  A1  ONE LEVEL at a time. Projecting onto the multilevel ladder would need
 *      its Gram matrix — the very degenerate object under investigation — and
 *      would smuggle the same confound back in.
 *  A2  the measurement window STRADDLES the structure. A plane wave lies in the
 *      span exactly at any W, so measuring in free space is a vacuous test that
 *      returns green regardless.
 *  A3  the carrier is LOCAL TO THE MEDIUM. Field F2 below shows what a global
 *      carrier costs at a material interface.
 *  B3  everything closed-form, so W/lambda can go to 1e6 without sampling.
 *  B4  a sampled witness at small W, since the analytic path shares src/phi.c
 *      with the production assembly and a shared bug would agree with itself. */
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>

enum { MAXDIM = 2048, MAXSEG = 4, NSAMP = 20000 };
static const double LAM = 16.0; /* cells per wavelength */
/* Gram of a single level is still exactly degenerate: the alternating sum of
 * phi translates vanishes identically (Nyquist null mode, M1). Min-norm is
 * therefore mandatory even here; any minimiser gives the same projection. */
static const double RCOND = 1e-12;

/* A field made of plane-wave pieces: amp * exp(i*kap*x) on [x0,x1). Everything
 * measured here is of this form, which is what makes the closed forms possible:
 * an incident/reflected pair, a transmitted wave past an interface, or the two
 * halves of an outgoing point-source field. */
typedef struct {
  double x0, x1;
  double complex amp;
  double kap;
} pwseg;

typedef struct {
  const char *name;
  pwseg seg[MAXSEG];
  int nseg;
  double kcar;  /* carrier for elements left of xsplit */
  double kcar2; /* carrier for elements right of xsplit; 0 = same as kcar */
  double xsplit;
} field;

/* Int over [a,b] of exp(i*kap*x) * phi(x/W-n) * {cos,sin}(kcar*x) dx */
static double complex proj_rhs_seg(double a, double b, hz_phi_factor f, double kap, double kcar,
                                   int quad) {
  double complex sp = hz_phi_integral_osc(a, b, f, kap + kcar);
  double complex sm = hz_phi_integral_osc(a, b, f, kap - kcar);
  if (quad == 0) return 0.5 * (sp + sm);
  return (sp - sm) / CMPLX(0.0, 2.0);
}

/* <B_i, B_j> over [a,b] with a carrier PER ELEMENT (ki may differ from kj —
 * that is the point of the local-carrier variant). The trig product collapses
 * to a difference term and a sum term:
 *   cos*cos = [cos(D) + cos(S)]/2      sin*sin = [cos(D) - cos(S)]/2
 *   cos_i*sin_j = [sin(S) - sin(D)]/2  sin_i*cos_j = [sin(S) + sin(D)]/2
 * with D = (ki-kj)x, S = (ki+kj)x. Equal carriers reduce to the DC + 2k form. */
static double gram_entry(double a, double b, hz_phi_factor fi, hz_phi_factor fj, int qi, int qj,
                         double ki, double kj) {
  double complex id = hz_phi_prod_integral_osc(a, b, fi, fj, ki - kj);
  double complex is = hz_phi_prod_integral_osc(a, b, fi, fj, ki + kj);
  if (qi == 0 && qj == 0) return 0.5 * (creal(id) + creal(is));
  if (qi == 1 && qj == 1) return 0.5 * (creal(id) - creal(is));
  if (qi == 0) return 0.5 * (cimag(is) - cimag(id));
  return 0.5 * (cimag(is) + cimag(id));
}

static double complex field_at(const field *fl, double x) {
  for (int s = 0; s < fl->nseg; s++)
    if (x >= fl->seg[s].x0 && x < fl->seg[s].x1)
      return fl->seg[s].amp * CMPLX(cos(fl->seg[s].kap * x), sin(fl->seg[s].kap * x));
  return 0.0;
}

/* Int |u|^2 over [a,b]. Pieces are disjoint, and each has |amp*exp(i*kap*x)|^2
 * constant for real kap, so this is exact. */
static double field_norm2(const field *fl, double a, double b) {
  double t = 0.0;
  for (int s = 0; s < fl->nseg; s++) {
    double lo = a > fl->seg[s].x0 ? a : fl->seg[s].x0;
    double hi = b < fl->seg[s].x1 ? b : fl->seg[s].x1;
    if (lo < hi)
      t += (creal(fl->seg[s].amp) * creal(fl->seg[s].amp) +
            cimag(fl->seg[s].amp) * cimag(fl->seg[s].amp)) *
           (hi - lo);
  }
  return t;
}

/* <u, B_i> over [a,b], summed over the field's plane-wave pieces */
static double complex field_rhs(const field *fl, double a, double b, hz_phi_factor f, int quad,
                                double kcar) {
  double complex t = 0.0;
  for (int s = 0; s < fl->nseg; s++) {
    double lo = a > fl->seg[s].x0 ? a : fl->seg[s].x0;
    double hi = b < fl->seg[s].x1 ? b : fl->seg[s].x1;
    if (lo >= hi) continue;
    t += fl->seg[s].amp * proj_rhs_seg(lo, hi, f, fl->seg[s].kap, kcar, quad);
  }
  return t;
}

/* Half-line a basis function is restricted to when the element is CUT by the
 * material boundary: -1 = uncut (whole support), 0 = left of xsplit, 1 = right.
 * Cutting is how "the basis must be discontinuous AT the boundary" is realised
 * without shrinking anything (PLAN, section on cuts). */
static void side_range(int s, double xsplit, double *lo, double *hi) {
  double big = 1e30;
  *lo = (s == 1) ? xsplit : -big;
  *hi = (s == 0) ? xsplit : big;
}

/* Relative L2 error of the best approximation of fl by one level of width W.
 * Fit over [fa,fb], measure over the inner half [ma,mb] so that the elements
 * hanging off the ends of the fit window cannot fake the result.
 * cut != 0 splits every element at fl->xsplit into two independent functions.
 * jump_out (may be NULL) receives |v(x+) - v(x-)| / |u| at the boundary — the
 * projection is in L2, which does not see continuity, so a cut basis will
 * happily produce a jump. That number is the direct measure of the
 * non-conformity flagged as PLAN open question 10; it is NOT fixed by cutting. */
static double project_cut(const field *fl, double W, double fa, double fb, double ma, double mb,
                          int cut, int *rank_out, int *dim_out, double *jump_out) {
  int n0 = (int)(fa / W) - 3, n1 = (int)(fb / W) + 3;
  int dim = 0;
  static hz_phi_factor fac[MAXDIM];
  static int quad[MAXDIM];
  static double kc[MAXDIM];
  static int sd[MAXDIM];
  double k2 = fl->kcar2 > 0.0 ? fl->kcar2 : fl->kcar;
  int nside = cut ? 2 : 1;
  for (int n = n0; n <= n1; n++)
    for (int s = 0; s < nside; s++)
      for (int q = 0; q < 2; q++) {
        if (dim >= MAXDIM) break;
        fac[dim] = (hz_phi_factor){W, (double)n, 0};
        quad[dim] = q;
        sd[dim] = cut ? s : -1;
        /* the carrier follows the medium this piece sits in (audit A3) */
        kc[dim] = cut ? (s == 0 ? fl->kcar : k2) : (((double)n * W < fl->xsplit) ? fl->kcar : k2);
        dim++;
      }
  double complex *G = calloc((size_t)dim * (size_t)dim, sizeof(double complex));
  double complex *r = calloc((size_t)dim, sizeof(double complex));
  double *sv = calloc((size_t)dim, sizeof(double));
  if (!G || !r || !sv) {
    free(G);
    free(r);
    free(sv);
    return -1.0;
  }
  for (int i = 0; i < dim; i++) {
    double li, hri;
    side_range(sd[i], fl->xsplit, &li, &hri);
    for (int j = 0; j < dim; j++) {
      double lj, hrj;
      side_range(sd[j], fl->xsplit, &lj, &hrj);
      double a = fa > li ? fa : li, b = fb < hri ? fb : hri;
      a = a > lj ? a : lj;
      b = b < hrj ? b : hrj;
      if (a >= b) continue; /* opposite sides of the cut: disjoint supports */
      G[(size_t)i * (size_t)dim + (size_t)j] =
          gram_entry(a, b, fac[i], fac[j], quad[i], quad[j], kc[i], kc[j]);
    }
    double a = fa > li ? fa : li, b = fb < hri ? fb : hri;
    if (a < b) r[i] = field_rhs(fl, a, b, fac[i], quad[i], kc[i]);
  }
  lapack_int rank = 0;
  lapack_int info = LAPACKE_zgelsd(LAPACK_ROW_MAJOR, dim, dim, 1, G, dim, r, 1, sv, RCOND, &rank);
  double err = -1.0;
  if (info == 0) {
    /* ||u-v||^2 = ||u||^2 - 2 Re<u,v> + <v,v>, ALL over the measure window.
     * v = sum c_i B_i with B_i real, so <u,v> = sum conj(c_i) * Int_M u B_i and
     * <v,v> = sum_ij c_i conj(c_j) Int_M B_i B_j. The solve left c in r[]. */
    double nu = field_norm2(fl, ma, mb);
    double complex uv = 0.0;
    for (int i = 0; i < dim; i++) {
      double li, hri;
      side_range(sd[i], fl->xsplit, &li, &hri);
      double a = ma > li ? ma : li, b = mb < hri ? mb : hri;
      if (a < b) uv += conj(r[i]) * field_rhs(fl, a, b, fac[i], quad[i], kc[i]);
    }
    double vv = 0.0;
    for (int i = 0; i < dim; i++) {
      double li, hri;
      side_range(sd[i], fl->xsplit, &li, &hri);
      for (int j = 0; j < dim; j++) {
        double lj, hrj;
        side_range(sd[j], fl->xsplit, &lj, &hrj);
        double a = ma > li ? ma : li, b = mb < hri ? mb : hri;
        a = a > lj ? a : lj;
        b = b < hrj ? b : hrj;
        if (a >= b) continue;
        vv += creal(r[i] * conj(r[j])) *
              gram_entry(a, b, fac[i], fac[j], quad[i], quad[j], kc[i], kc[j]);
      }
    }
    double e2 = nu - 2.0 * creal(uv) + vv;
    if (e2 < 0.0) e2 = 0.0;
    err = sqrt(e2 / nu);

    if (jump_out != NULL) {
      /* Jump of the approximation across the boundary, normalised by |u|.
       * Evaluated AT xsplit exactly, one-sided sums: each cut function is
       * defined up to (or from) the boundary, so no epsilon offset is needed.
       * An earlier version used x = xsplit +/- 1e-6*W, which made the offset
       * grow with the element and reported the SLOPE discontinuity times that
       * offset instead of the value jump — the column then rose linearly with W
       * and looked like catastrophic non-conformity at large elements. */
      double complex vm = 0.0, vp = 0.0;
      double x = fl->xsplit;
      for (int i = 0; i < dim; i++) {
        double t = x / fac[i].h - fac[i].n;
        if (fabs(t) >= 2.0) continue;
        double ph = hz_phi(t);
        double complex c = r[i] * (quad[i] ? ph * sin(kc[i] * x) : ph * cos(kc[i] * x));
        if (sd[i] == 1)
          vp += c;
        else
          vm += c;
      }
      *jump_out = cabs(vp - vm) / cabs(field_at(fl, x));
    }
  }
  *rank_out = (int)rank;
  *dim_out = dim;
  free(G);
  free(r);
  free(sv);
  return err;
}

/* Independent witness: same projection, but the error integral is sampled
 * rather than assembled from closed forms (B4). Only affordable at small W. */
static double project_sampled(const field *fl, double W, double fa, double fb, double ma,
                              double mb) {
  int n0 = (int)(fa / W) - 3, n1 = (int)(fb / W) + 3;
  int dim = 0;
  static hz_phi_factor fac[MAXDIM];
  static int quad[MAXDIM];
  static double kc[MAXDIM];
  double k2 = fl->kcar2 > 0.0 ? fl->kcar2 : fl->kcar;
  for (int n = n0; n <= n1; n++)
    for (int q = 0; q < 2; q++) {
      if (dim >= MAXDIM) break;
      fac[dim] = (hz_phi_factor){W, (double)n, 0};
      quad[dim] = q;
      /* the carrier follows the medium the element sits in (audit A3) */
      kc[dim] = ((double)n * W < fl->xsplit) ? fl->kcar : k2;
      dim++;
    }
  double complex *G = calloc((size_t)dim * (size_t)dim, sizeof(double complex));
  double complex *r = calloc((size_t)dim, sizeof(double complex));
  double *sv = calloc((size_t)dim, sizeof(double));
  if (!G || !r || !sv) {
    free(G);
    free(r);
    free(sv);
    return -1.0;
  }
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++)
      G[(size_t)i * (size_t)dim + (size_t)j] =
          gram_entry(fa, fb, fac[i], fac[j], quad[i], quad[j], kc[i], kc[j]);
    r[i] = field_rhs(fl, fa, fb, fac[i], quad[i], kc[i]);
  }
  lapack_int rank = 0;
  if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, dim, dim, 1, G, dim, r, 1, sv, RCOND, &rank) != 0) {
    free(G);
    free(r);
    free(sv);
    return -1.0;
  }
  double num = 0.0, den = 0.0;
  for (int p = 0; p < NSAMP; p++) {
    double x = ma + (mb - ma) * ((double)p + 0.5) / (double)NSAMP;
    double complex v = 0.0;
    for (int i = 0; i < dim; i++) {
      double t = x / fac[i].h - fac[i].n;
      if (fabs(t) >= 2.0) continue;
      double ph = hz_phi(t);
      v += r[i] * (quad[i] ? ph * sin(kc[i] * x) : ph * cos(kc[i] * x));
    }
    double complex u = field_at(fl, x);
    num += creal((u - v) * conj(u - v));
    den += creal(u * conj(u));
  }
  free(G);
  free(r);
  free(sv);
  return sqrt(num / den);
}

int main(void) {
  double k = 2.0 * M_PI / LAM;
  double kp = k * 1.5; /* the medium on the far side of an interface: n = 1.5 */
  double X0 = 0.0;     /* the structure sits at the origin of every window */

  printf("M9b/proj: best-approximation error of ONE level, window straddling the structure\n");
  printf("lambda = %.0f cells, interface contrast n = 1.5\n\n", LAM);

  /* F1 control: a plane wave at the carrier. Must be exact at ANY W — this is
   *     the claim the whole architecture rests on, so it is checked, not assumed.
   * F2 mismatched carrier: a plane wave at 1.5k represented by a carrier at k.
   *     This is what a GLOBAL carrier meets on the far side of an interface.
   * F3 point-source kink: exp(i k |x|), a genuine singularity — the outgoing
   *     field of a point source, worst case for a smooth envelope.
   * F4 same kink but with the two halves at the two media wavenumbers, carrier
   *     local to the left medium — an interface with a HALF-correct carrier. */
  double big = 1e9;
  field F[5] = {
      {"F1 plane wave, carrier matched", {{-big, big, 1.0, 0.0}}, 1, 0.0, 0.0, 0.0},
      {"F2 plane wave, carrier MISMATCHED", {{-big, big, 1.0, 0.0}}, 1, 0.0, 0.0, 0.0},
      {"F3 point-source kink exp(ik|x|)",
       {{-big, X0, 1.0, 0.0}, {X0, big, 1.0, 0.0}},
       2,
       0.0,
       0.0,
       0.0},
      {"F4 interface, ONE global carrier",
       {{-big, X0, 1.0, 0.0}, {X0, big, 1.0, 0.0}},
       2,
       0.0,
       0.0,
       0.0},
      {"F5 interface, carrier LOCAL to medium",
       {{-big, X0, 1.0, 0.0}, {X0, big, 1.0, 0.0}},
       2,
       0.0,
       0.0,
       0.0},
  };
  F[0].seg[0].kap = k;
  F[0].kcar = k;
  F[1].seg[0].kap = kp;
  F[1].kcar = k;
  F[2].seg[0].kap = -k;
  F[2].seg[1].kap = k;
  F[2].kcar = k;
  F[3].seg[0].kap = k;
  F[3].seg[1].kap = kp;
  F[3].kcar = k;
  /* F5: same field, but every element carries the wavenumber of the medium it
   * sits in — the architectural fix the audit demanded (A3) */
  F[4].seg[0].kap = k;
  F[4].seg[1].kap = kp;
  F[4].kcar = k;
  F[4].kcar2 = kp;
  F[4].xsplit = X0;

  printf("  %-38s %10s %8s %8s %10s\n", "field", "W/lambda", "dim", "rank", "rel L2 err");
  static const double WL[7] = {0.125, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5};
  for (int f = 0; f < 5; f++) {
    for (int i = 0; i < 7; i++) {
      double W = WL[i] * LAM;
      /* fit over 8 elements centred on the structure, measure over the inner 4;
       * the window always contains the structure (audit A2) */
      double fa = X0 - 4.0 * W, fb = X0 + 4.0 * W;
      double ma = X0 - 2.0 * W, mb = X0 + 2.0 * W;
      int rank = 0, dim = 0;
      double e = project_cut(&F[f], W, fa, fb, ma, mb, 0, &rank, &dim, NULL);
      printf("  %-38s %10.3f %8d %8d %10.3e\n", i == 0 ? F[f].name : "", WL[i], dim, rank, e);
    }
    printf("\n");
  }

  /* --- PLAN question 13: the phase-accuracy law -------------------------
   * Predicted: the carrier phase must stay within ~1 rad across an element, so
   * (k'-k)*W ~ 1, i.e. W_max ~ lambda/(2*pi*delta) with delta = dk/k. If that
   * is the law, W_max*delta/lambda is a CONSTANT down the table. This is the
   * number that says how accurately the medium has to be known — an input
   * parameter for gradient media and for geometric tolerances. */
  printf("[phase law] W at which a mismatched carrier crosses 10%% error\n");
  printf("  %10s  %12s  %14s\n", "dk/k", "W_max/lambda", "W_max*dk/lambda");
  static const double DK[6] = {0.5, 0.1, 0.01, 1e-3, 1e-4, 1e-5};
  for (int d = 0; d < 6; d++) {
    field G = {"", {{-1e9, 1e9, 1.0, 0.0}}, 1, k, 0.0, 0.0};
    G.seg[0].kap = k * (1.0 + DK[d]);
    double wlo = 1e-4, whi = 1e7; /* in lambda */
    for (int it = 0; it < 40; it++) {
      double wm = sqrt(wlo * whi), W = wm * LAM;
      int rk = 0, dm = 0;
      double e = project_cut(&G, W, -4.0 * W, 4.0 * W, -2.0 * W, 2.0 * W, 0, &rk, &dm, NULL);
      if (e < 0.1)
        wlo = wm;
      else
        whi = wm;
    }
    printf("  %10.0e  %12.4g  %14.4g\n", DK[d], wlo, wlo * DK[d]);
  }
  printf("  constant last column => law confirmed; drifting => it is not this law\n\n");

  /* --- PLAN question 10: does CUTTING repair the straddling element? -----
   * F5 saturates at ~20% however wide the element, because one smooth envelope
   * cannot hold "wave left" on one side and "wave right" on the other; the
   * misfit smears over the whole support. Cutting gives each side its own
   * coefficient AND its own carrier. The jump column is the price: L2
   * projection does not see continuity, so the cut basis is non-conforming —
   * that is a real open problem (question 10), not an artefact. */
  printf("[cut] straddling element, smooth vs cut at the boundary\n");
  printf("  %10s %10s %10s %12s %10s\n", "W/lambda", "err smooth", "err cut", "dim smooth/cut",
         "jump/|u|");
  for (int i = 0; i < 7; i++) {
    double W = WL[i] * LAM;
    double fa = X0 - 4.0 * W, fb = X0 + 4.0 * W;
    double ma = X0 - 2.0 * W, mb = X0 + 2.0 * W;
    int r1 = 0, d1 = 0, r2 = 0, d2 = 0;
    double jump = -1.0;
    double es = project_cut(&F[4], W, fa, fb, ma, mb, 0, &r1, &d1, NULL);
    double ec = project_cut(&F[4], W, fa, fb, ma, mb, 1, &r2, &d2, &jump);
    printf("  %10.3f %10.3e %10.3e %6d/%-5d %10.3e\n", WL[i], es, ec, d1, d2, jump);
  }
  printf("\n");

  /* --- PLAN question 10, SECOND half: CONFORMITY -------------------------
   * The [cut] phase above showed that cutting restores exact representability,
   * but it could NOT judge conformity: each side there was exactly
   * representable, so v == u and the continuity was inherited from the field
   * rather than produced by the basis — there was nothing for conformity to
   * fail on. Here the right-hand side is only APPROXIMATELY representable (the
   * field reverses direction partway past the cut, which one carrier cannot
   * hold), so the two independent one-sided projections have no reason to agree
   * at the boundary. The jump they leave is the real measure of the structural
   * non-conformity that a Nitsche / constraint / Heaviside treatment would have
   * to remove. */
  printf("[conformity] cut, field only APPROXIMATELY representable on one side\n");
  printf("  %10s %12s %12s\n", "W/lambda", "err cut", "jump/|u|");
  field FC = {
      "reversal past the cut",
      {{-big, X0, 1.0, 0.0}, {X0, X0 + 40.0 * LAM, 1.0, 0.0}, {X0 + 40.0 * LAM, big, 1.0, 0.0}},
      3,
      0.0,
      0.0,
      0.0};
  FC.seg[0].kap = k;
  FC.seg[1].kap = kp;
  FC.seg[2].kap = -kp;
  FC.kcar = k;
  FC.kcar2 = kp;
  FC.xsplit = X0;
  for (int i = 0; i < 5; i++) {
    double W = WL[i] * LAM;
    double fa = X0 - 4.0 * W, fb = X0 + 4.0 * W;
    double ma = X0 - 2.0 * W, mb = X0 + 2.0 * W;
    int rr = 0, dd = 0;
    double jj = -1.0;
    double ee = project_cut(&FC, W, fa, fb, ma, mb, 1, &rr, &dd, &jj);
    printf("  %10.3f %12.4e %12.4e\n", WL[i], ee, jj);
  }
  printf("  READ: a non-zero jump here is structural non-conformity — the price\n");
  printf("  a Nitsche / constraint / Heaviside treatment would have to pay.\n\n");

  printf("[witness] sampled vs closed-form at small W (audit B4)\n");
  for (int f = 0; f < 5; f++) {
    double W = 1.0 * LAM;
    double fa = X0 - 4.0 * W, fb = X0 + 4.0 * W;
    double ma = X0 - 2.0 * W, mb = X0 + 2.0 * W;
    int rank = 0, dim = 0;
    double ea = project_cut(&F[f], W, fa, fb, ma, mb, 0, &rank, &dim, NULL);
    double es = project_sampled(&F[f], W, fa, fb, ma, mb);
    printf("  %-38s analytic %.4e  sampled %.4e  diff %.1e\n", F[f].name, ea, es, fabs(ea - es));
  }
  return 0;
}
