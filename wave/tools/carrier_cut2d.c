/* M14 — does the CUT survive the move from 1D to 2D on an OBLIQUE boundary?
 * Audit, falsifiers and predicted failures: result/M14_AUDIT.md (written first).
 *
 * The whole of "LOD on the basis is free" rests on the cut: without one, an
 * element straddling a material boundary sits on a 20.6% floor that does not
 * fall with width (M9b F5). Cutting removed that floor to machine zero — but in
 * 1D, where the boundary is a POINT and the cut is a split of an interval. Here
 * the boundary is a LINE at an arbitrary angle, the cut breaks separability, and
 * nobody has measured whether the 1D result transfers.
 *
 * Projection only: no operator, no domain boundary, no solver of the physics —
 * the arrangement of tools/carrier_proj.c (M9b §5-§8) lifted into 2D. Failure
 * here is decisive; success promises nothing about a Galerkin solve.
 *
 * The exact field is a ROTATED 1D layered problem, so the reference is closed
 * form at any angle and any element width:
 *   side 1:  e^{i K_i.x} + R e^{i K_r.x}      side 2:  T e^{i K_t.x}
 * with the tangential wavenumber conserved and R, T from continuity of u and
 * du/dn. NC3 below drives alpha and theta to zero, where the problem becomes
 * y-independent and the answer must reproduce the published 1D column. */
#include "cut2d.h"
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum { MAXDIM = 1400, NNODE = 7, MAXPIECE = 5 };
static const double LAM = 16.0;
/* Single-level Gram, still exactly degenerate (the alternating sum of phi
 * translates vanishes identically, M1), so min-norm is mandatory. The value is
 * the one carrier_proj.c used; F5 of the audit sweeps it. */
static const double RCOND_DEF = 1e-12;

typedef struct {
  double x, y;
} vec2;

typedef struct {
  int nx, ny;
  int side; /* -1 = uncut (whole support), 0 = s < scut, 1 = s > scut */
  vec2 K;   /* carrier wavevector */
} elem;

typedef struct {
  double complex amp;
  vec2 K;
  int side; /* 0 = s < 0, 1 = s > 0 — the TRUE medium boundary */
} piece;

typedef struct {
  double ca, sa;      /* interface normal */
  double scut;        /* where the BASIS is cut (0 = on the boundary) */
  piece pc[MAXPIECE]; /* the exact field */
  int npc;
  vec2 Ki, Kr, Kt; /* incident, reflected, transmitted (the FIELD) */
  vec2 KtB;        /* what the BASIS uses on side + (NC1 breaks only this) */
  double k1, k2;
} scene;

/* --- geometry helpers ---------------------------------------------------- */
static void side_strip(int side, double ca, double sa, double scut, hz_strip2 *st) {
  st->ca = ca;
  st->sa = sa;
  st->slo = (side == 1) ? scut : -HZ_CUT_INF;
  st->shi = (side == 0) ? scut : HZ_CUT_INF;
}

/* intersection of two strips with the SAME normal */
static hz_strip2 strip_and(hz_strip2 a, hz_strip2 b) {
  hz_strip2 r = a;
  if (b.slo > r.slo) r.slo = b.slo;
  if (b.shi < r.shi) r.shi = b.shi;
  return r;
}

static void support(int n, double W, double lo, double hi, double *a, double *b) {
  *a = W * ((double)n - 2.0);
  *b = W * ((double)n + 2.0);
  if (*a < lo) *a = lo;
  if (*b > hi) *b = hi;
}

/* --- integrals ----------------------------------------------------------- */
/* <B_a, B_b> over the window square [lo,hi]^2, sesquilinear (the projection
 * conjugates; the Galerkin form does NOT — see PLAN, API note). */
static double complex gram(const elem *a, const elem *b, double W, double lo, double hi,
                           const scene *sc) {
  double xa, xb, ya, yb, ta, tb;
  support(a->nx, W, lo, hi, &xa, &xb);
  support(b->nx, W, xa, xb, &ta, &tb);
  support(a->ny, W, lo, hi, &ya, &yb);
  double ua, ub;
  support(b->ny, W, ya, yb, &ua, &ub);
  if (!(ta < tb) || !(ua < ub)) return 0.0;
  hz_strip2 sa2, sb2;
  side_strip(a->side, sc->ca, sc->sa, sc->scut, &sa2);
  side_strip(b->side, sc->ca, sc->sa, sc->scut, &sb2);
  hz_strip2 st = strip_and(sa2, sb2);
  if (!(st.slo < st.shi)) return 0.0;
  hz_phi_factor fax = {W, (double)a->nx, 0}, fbx = {W, (double)b->nx, 0};
  hz_phi_factor fay = {W, (double)a->ny, 0}, fby = {W, (double)b->ny, 0};
  hz_axis2 fx = {{fax, fbx}, 2}, fy = {{fay, fby}, 2};
  double cax = W * (double)a->nx, cay = W * (double)a->ny;
  double cbx = W * (double)b->nx, cby = W * (double)b->ny;
  double complex ph =
      cexp(CMPLX(0.0, 1.0) * (a->K.x * cax + a->K.y * cay - b->K.x * cbx - b->K.y * cby));
  return ph * hz_cut2d_integral(ta, tb, ua, ub, fx, fy, b->K.x - a->K.x, b->K.y - a->K.y, st);
}

/* <B_a, u> over the window square */
static double complex rhs(const elem *a, double W, double lo, double hi, const scene *sc) {
  double xa, xb, ya, yb;
  support(a->nx, W, lo, hi, &xa, &xb);
  support(a->ny, W, lo, hi, &ya, &yb);
  if (!(xa < xb) || !(ya < yb)) return 0.0;
  hz_strip2 sa2;
  side_strip(a->side, sc->ca, sc->sa, sc->scut, &sa2);
  hz_phi_factor fax = {W, (double)a->nx, 0}, fay = {W, (double)a->ny, 0};
  hz_axis2 fx = {{fax, fax}, 1}, fy = {{fay, fay}, 1};
  double cax = W * (double)a->nx, cay = W * (double)a->ny;
  double complex ph = cexp(CMPLX(0.0, 1.0) * (a->K.x * cax + a->K.y * cay));
  double complex tot = 0.0;
  for (int p = 0; p < sc->npc; p++) {
    hz_strip2 fs;
    side_strip(sc->pc[p].side, sc->ca, sc->sa, 0.0, &fs);
    hz_strip2 st = strip_and(sa2, fs);
    if (!(st.slo < st.shi)) continue;
    tot += sc->pc[p].amp * hz_cut2d_integral(xa, xb, ya, yb, fx, fy, sc->pc[p].K.x - a->K.x,
                                             sc->pc[p].K.y - a->K.y, st);
  }
  return ph * tot;
}

/* ||u||^2 over the window square */
static double field_norm2(double lo, double hi, const scene *sc) {
  hz_phi_factor dummy = {1.0, 0.0, 0};
  hz_axis2 f0 = {{dummy, dummy}, 0};
  double complex tot = 0.0;
  for (int p = 0; p < sc->npc; p++)
    for (int q = 0; q < sc->npc; q++) {
      if (sc->pc[p].side != sc->pc[q].side) continue;
      hz_strip2 st;
      side_strip(sc->pc[p].side, sc->ca, sc->sa, 0.0, &st);
      tot += conj(sc->pc[p].amp) * sc->pc[q].amp *
             hz_cut2d_integral(lo, hi, lo, hi, f0, f0, sc->pc[q].K.x - sc->pc[p].K.x,
                               sc->pc[q].K.y - sc->pc[p].K.y, st);
    }
  return creal(tot);
}

/* --- basis construction --------------------------------------------------- */
/* mode: 0 = uncut, carrier from the element centre's medium
 *       1 = uncut, ALL directions on straddling elements ("poor man's cut")
 *       2 = cut, each side its own directions
 *       3 = cut, ONE global carrier magnitude (k1) on both sides */
static int build(elem *b, double W, int mode, const scene *sc) {
  int dim = 0;
  for (int nx = -NNODE; nx <= NNODE; nx++)
    for (int ny = -NNODE; ny <= NNODE; ny++) {
      double cx = W * (double)nx, cy = W * (double)ny;
      double s = cx * sc->ca + cy * sc->sa;
      double reach = 2.0 * W * (fabs(sc->ca) + fabs(sc->sa));
      int straddles = fabs(s - sc->scut) < reach;
      vec2 Kt2 = sc->KtB;
      if (mode == 3) { /* same direction, magnitude of medium 1 (M9b F4) */
        double nrm = sqrt(sc->KtB.x * sc->KtB.x + sc->KtB.y * sc->KtB.y);
        Kt2.x = sc->KtB.x * sc->k1 / nrm;
        Kt2.y = sc->KtB.y * sc->k1 / nrm;
      }
      if (mode <= 1) {
        int minus = s < 0.0;
        if (mode == 1 && straddles) {
          b[dim++] = (elem){nx, ny, -1, sc->Ki};
          b[dim++] = (elem){nx, ny, -1, sc->Kr};
          b[dim++] = (elem){nx, ny, -1, Kt2};
        } else if (minus) {
          b[dim++] = (elem){nx, ny, -1, sc->Ki};
          b[dim++] = (elem){nx, ny, -1, sc->Kr};
        } else {
          b[dim++] = (elem){nx, ny, -1, Kt2};
        }
      } else {
        if (straddles || s < sc->scut) {
          b[dim++] = (elem){nx, ny, straddles ? 0 : -1, sc->Ki};
          b[dim++] = (elem){nx, ny, straddles ? 0 : -1, sc->Kr};
        }
        if (straddles || s > sc->scut) b[dim++] = (elem){nx, ny, straddles ? 1 : -1, Kt2};
      }
      if (dim > MAXDIM - 8) return dim;
    }
  return dim;
}

/* --- conformity: the jump the CUT leaves on the interface ------------------
 * The cut makes the basis discontinuous on purpose, and an L2 projection does
 * not penalise a jump at all — so the two one-sided expansions have no reason
 * to agree. In the exact-field runs above they must agree anyway, because each
 * side is representable EXACTLY and continuity is inherited from the field
 * rather than produced by the basis; that is the hole M9b §7 found in itself
 * and closed in §8 with a field only APPROXIMATELY representable. Here the same
 * device: a ghost plane wave whose direction the basis does not carry, added on
 * BOTH sides so the field itself stays continuous (adding it on one side only
 * would make u jump, and the measurement would be meaningless).
 * Sampled on the cut line inside the measurement window; uncut functions are
 * continuous and cancel in the difference, so only cut ones are summed. */
static double cut_jump(const elem *b, int dim, const double complex *c, double W, const scene *sc) {
  double tx = -sc->sa, ty = sc->ca;
  double worst = 0.0;
  for (int p = 0; p <= 64; p++) {
    double t = -1.5 * W + 3.0 * W * (double)p / 64.0;
    double x = t * tx + sc->scut * sc->ca, y = t * ty + sc->scut * sc->sa;
    double complex vm = 0.0, vp = 0.0, u = 0.0;
    for (int i = 0; i < dim; i++) {
      if (b[i].side < 0) continue;
      double ax = x / W - (double)b[i].nx, ay = y / W - (double)b[i].ny;
      if (fabs(ax) >= 2.0 || fabs(ay) >= 2.0) continue;
      double complex v = c[i] * hz_phi(ax) * hz_phi(ay) *
                         cexp(CMPLX(0.0, 1.0) * (b[i].K.x * (x - W * (double)b[i].nx) +
                                                 b[i].K.y * (y - W * (double)b[i].ny)));
      if (b[i].side == 1)
        vp += v;
      else
        vm += v;
    }
    for (int q = 0; q < sc->npc; q++)
      if (sc->pc[q].side == 0)
        u += sc->pc[q].amp * cexp(CMPLX(0.0, 1.0) * (sc->pc[q].K.x * x + sc->pc[q].K.y * y));
    double j = cabs(vp - vm) / cabs(u);
    if (j > worst) worst = j;
  }
  return worst;
}

/* --- one projection ------------------------------------------------------- */
static double project(double W, int mode, const scene *sc, double rcond, int *dim_out,
                      int *rank_out, int on_fit, double *jump_out) {
  static elem b[MAXDIM];
  int dim = build(b, W, mode, sc);
  double fa = -4.0 * W, fb = 4.0 * W, ma = -2.0 * W, mb = 2.0 * W;
  if (on_fit) { /* NC4: measure where the fit was made */
    ma = fa;
    mb = fb;
  }
  double complex *G = calloc((size_t)dim * (size_t)dim, sizeof(double complex));
  double complex *r = calloc((size_t)dim, sizeof(double complex));
  double complex *GM = calloc((size_t)dim * (size_t)dim, sizeof(double complex));
  double complex *rM = calloc((size_t)dim, sizeof(double complex));
  double *sv = calloc((size_t)dim, sizeof(double));
  if (!G || !r || !GM || !rM || !sv) {
    free(G);
    free(r);
    free(GM);
    free(rM);
    free(sv);
    return -1.0;
  }
  for (int i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++) {
      G[(size_t)i * (size_t)dim + (size_t)j] = gram(&b[i], &b[j], W, fa, fb, sc);
      GM[(size_t)i * (size_t)dim + (size_t)j] = gram(&b[i], &b[j], W, ma, mb, sc);
    }
    r[i] = rhs(&b[i], W, fa, fb, sc);
    rM[i] = rhs(&b[i], W, ma, mb, sc);
  }
  lapack_int rank = 0;
  lapack_int info = LAPACKE_zgelsd(LAPACK_ROW_MAJOR, dim, dim, 1, G, dim, r, 1, sv, rcond, &rank);
  double err = -1.0;
  if (info == 0) {
    double nu = field_norm2(ma, mb, sc);
    double complex uv = 0.0;
    for (int i = 0; i < dim; i++)
      uv += r[i] * conj(rM[i]);
    double complex vv = 0.0;
    for (int i = 0; i < dim; i++)
      for (int j = 0; j < dim; j++)
        vv += conj(r[i]) * r[j] * GM[(size_t)i * (size_t)dim + (size_t)j];
    double e2 = nu - 2.0 * creal(uv) + creal(vv);
    if (e2 < 0.0) e2 = 0.0;
    err = sqrt(e2 / nu);
  }
  *dim_out = dim;
  if (jump_out != NULL) *jump_out = cut_jump(b, dim, r, W, sc);
  *rank_out = (int)rank;
  free(G);
  free(r);
  free(GM);
  free(rM);
  free(sv);
  return err;
}

/* --- scene ---------------------------------------------------------------- */
static scene make_scene(double alpha, double theta, double contrast, double scut, int snell_break,
                        double ghost) {
  scene sc;
  sc.ca = cos(alpha);
  sc.sa = sin(alpha);
  sc.scut = scut;
  sc.k1 = 2.0 * M_PI / LAM;
  sc.k2 = sc.k1 * contrast;
  double tx = -sin(alpha), ty = cos(alpha);
  double kt = sc.k1 * sin(theta), k1n = sc.k1 * cos(theta);
  double arg = sc.k2 * sc.k2 - kt * kt;
  double k2n = arg > 0.0 ? sqrt(arg) : 0.0;
  double complex R = (k1n - k2n) / (k1n + k2n), T = 1.0 + R;
  sc.Ki = (vec2){sc.k1 * sin(theta) * tx + k1n * sc.ca, sc.k1 * sin(theta) * ty + k1n * sc.sa};
  sc.Kr = (vec2){sc.k1 * sin(theta) * tx - k1n * sc.ca, sc.k1 * sin(theta) * ty - k1n * sc.sa};
  sc.Kt = (vec2){kt * tx + k2n * sc.ca, kt * ty + k2n * sc.sa};
  sc.KtB = sc.Kt;
  if (snell_break) {
    /* NC1 rotates the direction the BASIS carries on side +, and ONLY that.
     * The first version rotated the field's transmitted wave as well, which
     * left the two consistent: the control returned 1.16e-7, exactly the same
     * as the correct cut, and would have been read as "the cut works". Caught
     * because a negative control has no business landing on the metric floor. */
    double th2 = atan2(kt, k2n) + 10.0 * M_PI / 180.0;
    double ktb = sc.k2 * sin(th2), k2nb = sc.k2 * cos(th2);
    sc.KtB = (vec2){ktb * tx + k2nb * sc.ca, ktb * ty + k2nb * sc.sa};
  }
  sc.pc[0] = (piece){1.0, sc.Ki, 0};
  sc.pc[1] = (piece){R, sc.Kr, 0};
  sc.pc[2] = (piece){T, sc.Kt, 1};
  sc.npc = 3;
  if (ghost > 0.0) {
    /* a plane wave the basis does not carry, present on BOTH sides so that u
     * itself stays continuous across the interface (see cut_jump) */
    double gth = theta + 40.0 * M_PI / 180.0;
    vec2 Kg = {sc.k1 * (sin(gth) * tx + cos(gth) * sc.ca),
               sc.k1 * (sin(gth) * ty + cos(gth) * sc.sa)};
    sc.pc[3] = (piece){ghost, Kg, 0};
    sc.pc[4] = (piece){ghost, Kg, 1};
    sc.npc = 5;
  }
  return sc;
}

static const char *MODE_NAME[4] = {"1a uncut, carrier by centre", "1b uncut, ALL directions",
                                   "2  CUT, local carriers", "4  CUT, one global carrier"};

int main(void) {
  static const double WL[7] = {0.125, 1.0, 10.0, 100.0, 1e3, 1e4, 1e5};
  printf("M14: the cut on an OBLIQUE material boundary, projection only\n");
  printf("lambda = %.0f cells, contrast n = 1.5, incidence 30 deg\n\n", LAM);

  /* --- NC3 (instrument): alpha = 0, theta = 0 must reproduce the 1D column -- */
  printf("[NC3] alpha=0, theta=0: y-independent, must reproduce M9b (0.0034/0.040/0.206)\n");
  printf("  %10s %14s %14s %10s\n", "W/lambda", "1a uncut", "2 cut", "dim/rank");
  for (int i = 0; i < 5; i++) {
    scene sc = make_scene(0.0, 0.0, 1.5, 0.0, 0, 0.0);
    int d1 = 0, r1 = 0, d2 = 0, r2 = 0;
    double e1 = project(WL[i] * LAM, 0, &sc, RCOND_DEF, &d1, &r1, 0, NULL);
    double e2 = project(WL[i] * LAM, 2, &sc, RCOND_DEF, &d2, &r2, 0, NULL);
    printf("  %10.3f %14.4e %14.4e %6d/%-5d\n", WL[i], e1, e2, d2, r2);
  }

  /* --- the four branches at 30 degrees ------------------------------------ */
  double alpha = 30.0 * M_PI / 180.0, theta = 30.0 * M_PI / 180.0;
  printf("\n[branches] interface at 30 deg, incidence 30 deg\n");
  printf("  %-28s %10s %12s %8s %8s\n", "branch", "W/lambda", "rel L2 err", "dim", "rank");
  for (int m = 0; m < 4; m++) {
    for (int i = 0; i < 7; i++) {
      scene sc = make_scene(alpha, theta, 1.5, 0.0, 0, 0.0);
      int d = 0, rk = 0;
      double e = project(WL[i] * LAM, m == 3 ? 3 : m, &sc, RCOND_DEF, &d, &rk, 0, NULL);
      printf("  %-28s %10.3f %12.4e %8d %8d\n", i == 0 ? MODE_NAME[m] : "", WL[i], e, d, rk);
    }
    printf("\n");
  }

  /* --- NC1: broken Snell in the cut basis --------------------------------- */
  {
    scene sc = make_scene(alpha, theta, 1.5, 0.0, 1, 0.0);
    int d = 0, rk = 0;
    double e = project(10.0 * LAM, 2, &sc, RCOND_DEF, &d, &rk, 0, NULL);
    printf("[NC1] cut with the transmitted direction rotated 10 deg, W=10lam: %.4e\n", e);
  }
  /* --- NC2: zero contrast: every branch must be exact ---------------------- */
  {
    for (int m = 0; m < 3; m++) {
      scene sc = make_scene(alpha, theta, 1.0, 0.0, 0, 0.0);
      int d = 0, rk = 0;
      double e = project(10.0 * LAM, m, &sc, RCOND_DEF, &d, &rk, 0, NULL);
      printf("[NC2] zero contrast, branch %d, W=10lam: %.4e\n", m, e);
    }
  }
  /* --- NC4: does the measurement window matter? --------------------------- */
  {
    printf("[NC4] same runs measured on the FIT window +-4W instead of +-2W\n");
    for (int m = 0; m < 3; m += 2) {
      scene sc = make_scene(alpha, theta, 1.5, 0.0, 0, 0.0);
      int d = 0, rk = 0;
      double ein = project(10.0 * LAM, m, &sc, RCOND_DEF, &d, &rk, 0, NULL);
      double eout = project(10.0 * LAM, m, &sc, RCOND_DEF, &d, &rk, 1, NULL);
      printf("       %-28s inner %.4e  fit-window %.4e\n", MODE_NAME[m], ein, eout);
    }
  }

  /* --- F3: the law of a DISPLACED cut ------------------------------------- */
  printf("\n[F3] cut displaced by delta from the true boundary (worst case: uniform)\n");
  printf("  %10s %12s %12s %12s %12s\n", "W/lambda", "k*d=0.01", "0.03", "0.1", "0.3");
  static const double KD[4] = {0.01, 0.03, 0.1, 0.3};
  for (int i = 1; i < 5; i++) {
    printf("  %10.3f", WL[i]);
    for (int j = 0; j < 4; j++) {
      double k1 = 2.0 * M_PI / LAM;
      scene sc = make_scene(alpha, theta, 1.5, KD[j] / k1, 0, 0.0);
      int d = 0, rk = 0;
      double e = project(WL[i] * LAM, 2, &sc, RCOND_DEF, &d, &rk, 0, NULL);
      printf(" %12.4e", e);
    }
    printf("\n");
  }

  /* --- F5: is the answer a truncation artefact? --------------------------- */
  printf("\n[F5] RCOND sweep at W=10lam (audit F5: must not move by more than 2x)\n");
  static const double RC[4] = {1e-14, 1e-12, 1e-10, 1e-8};
  for (int m = 0; m < 3; m += 2) {
    printf("  %-28s", MODE_NAME[m]);
    for (int i = 0; i < 4; i++) {
      scene sc = make_scene(alpha, theta, 1.5, 0.0, 0, 0.0);
      int d = 0, rk = 0;
      double e = project(10.0 * LAM, m, &sc, RC[i], &d, &rk, 0, NULL);
      printf(" %11.4e", e);
    }
    printf("\n");
  }

  /* --- conformity: the jump left by the cut (M9b §8 lifted into 2D) -------
   * PREDICTION, RECORDED BEFORE THE RUN: with the exact field both sides are
   * representable exactly, so the jump stays at the metric floor and says
   * nothing (that is the hole §7 had). With a ghost wave the basis cannot
   * carry, the approximation error jumps to ~0.1-0.3 while the JUMP stays two
   * to three orders below it, as in 1D (0.03-0.3% against 12-17%). */
  printf("\n[conformity] jump across the cut, max over the interface in the window\n");
  printf("  %10s %10s %14s %14s %10s\n", "W/lambda", "ghost", "err", "jump/|u|", "jump/err");
  static const double GH[3] = {0.0, 0.1, 0.3};
  for (int g = 0; g < 3; g++)
    for (int i = 1; i < 4; i++) {
      scene sc = make_scene(alpha, theta, 1.5, 0.0, 0, GH[g]);
      int d = 0, rk = 0;
      double jmp = -1.0;
      double e = project(WL[i] * LAM, 2, &sc, RCOND_DEF, &d, &rk, 0, &jmp);
      printf("  %10.3f %10.2f %14.4e %14.4e %10.2e\n", WL[i], GH[g], e, jmp,
             e > 0.0 ? jmp / e : -1.0);
    }

  /* --- angle sweep at W = 10 lambda --------------------------------------- */
  printf("\n[angles] W = 10 lambda\n");
  printf("  %10s %14s %14s\n", "alpha,deg", "1a uncut", "2 cut");
  static const double AL[4] = {0.0, 15.0, 30.0, 45.0};
  for (int i = 0; i < 4; i++) {
    scene sc = make_scene(AL[i] * M_PI / 180.0, theta, 1.5, 0.0, 0, 0.0);
    int d = 0, rk = 0;
    double e1 = project(10.0 * LAM, 0, &sc, RCOND_DEF, &d, &rk, 0, NULL);
    double e2 = project(10.0 * LAM, 2, &sc, RCOND_DEF, &d, &rk, 0, NULL);
    printf("  %10.1f %14.4e %14.4e\n", AL[i], e1, e2);
  }
  return 0;
}
