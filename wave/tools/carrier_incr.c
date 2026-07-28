/* M11: CORRECTION ITERATIONS WHEN THE CAMERA MOVES.
 *
 * WHY THIS IS THE CRITICAL PATH. SOLVE_REPORT measured that the COLD solve
 * does not give real time: LSQR converges, but 165 -> 320 iterations while the
 * domain grows by 1000x (~0.4*dim). At the 3D estimate of 1e8 unknowns that is
 * ~5 s per frame. The real-time claim therefore rests on INCREMENTALITY alone,
 * and this bench is the only place where incrementality becomes a number.
 *
 * The full pre-run audit is result/M11_AUDIT.md. The six things it rejected,
 * because each of them would have produced a green result that means nothing:
 *
 * R1. The existing operator benches (carrier_solve.c, carrier_iter.c) anchor
 *     the shells at the SOURCE, not at a camera. "Moving the camera" there
 *     moves the refinement around the source: a different phenomenon under the
 *     right name. Here the camera is its own anchor and the source never moves.
 * R2. A vacuum scene self-confirms. A plane wave lies in the span of the
 *     carrier basis EXACTLY at any element width (M9b F1, zero to 1e5 lambda),
 *     so re-projecting the old solution would be near-exact for a reason that
 *     has nothing to do with incrementality, and the bench would print "0
 *     iterations" as a property of the method. The scene must carry structure
 *     the basis holds only approximately. Vacuum is still run - as the FLOOR of
 *     the instrument, not as the result.
 * R3. Convergence judged on the residual. M2's lesson: the raw residual
 *     stagnates on a boundary mode that does not hurt the field. FIELD only.
 * R4. Measuring the field only at the camera, where the basis is finest BY
 *     CONSTRUCTION. Both are measured: at the camera (the renderer's output)
 *     and over the whole domain (honesty).
 * R5. Letting the DtN row scaling follow the basis. It is derived from the
 *     assembled matrix, so it would differ between the old and the new camera
 *     position and the difference would be charged to the camera. Fixed once,
 *     at frame 0, and reused.
 * R6. Measuring a single step. One small shift is cheap almost by definition.
 *     The claim is 30 frames a second indefinitely, and partial corrections
 *     ACCUMULATE, so falsifier 4 (drift over 64 frames at a fixed budget) was
 *     added before the run.
 *
 * DEVIATION FROM THE PLAN, STATED OUT LOUD: this is 1D, where PLAN M11 asks for
 * 2D. All four falsifiers are expressible in 1D (the camera is the point R is
 * measured from in L = eps*R; "no camera in 1D" was about ANGULAR structure).
 * A 2D operator does not exist yet and its termination is an open construction
 * (TERMINATION_REPORT: DtN is nonlocal in 2D). What 1D cannot see: whether
 * camera motion re-selects carrier DIRECTIONS. Necessary, not sufficient. */
#include "carrier.h"
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum {
  MAXB = 2400,     /* basis capacity; measured worst case ~1050 */
  NLEVMAX = 26,    /* level ladder cap */
  MARGIN = 2,      /* shell margin in elements (SHELLS_REPORT rule) */
  SRCLEV = -1,     /* source ladder depth; -1 = full dyadic ladder, like the
                    * camera one. A 3-level patch was tried first and FAILED:
                    * beyond it a single coarse element straddles the kink of
                    * the emitter and is spoiled over its whole support (M9b F5
                    * with the jump supplied by the source instead of the
                    * medium) - global error 15-25% in VACUUM, where the field
                    * is otherwise exact. The source is geometry, and regulator
                    * 2 of the plan says geometry sets a camera-INDEPENDENT
                    * floor: a full ladder, not a patch. */
  MAXSLAB = 6,     /* thin slabs making up the depth spread */
  MAXSEG = 16,     /* medium segments = 2*MAXSLAB+2, rounded up */
  NPROBE = 400,    /* probe points at the camera (checked every iteration) */
  NPROBEG = 1600,  /* probe points over the whole domain */
  MAXIT = 400,     /* LSQR cap; cold needs 320 at the largest domain */
  NFRAME = 64,     /* frames in the drift test */
  DRIFT_DIRECT = 8 /* recompute the direct-solve floor every N frames */
};

static const double LAM = 16.0;       /* cells per wavelength */
static const double W0 = 2.0;         /* finest element = lambda/8 */
static const double EPS = 0.125;      /* LOD constant; M9b: eps <= 1/8 or the shell
                                       * band is narrower than the support */
static const double WMAX_FRAC = 0.02; /* coarsest element vs domain (M9b: an
                                       * unbounded ladder swallows the source) */
static const double RCOND = 1e-12;
static const double TARGET = 0.01; /* the pre-registered 1% of the FIELD */
/* ARTEFACT 13, CAUGHT BY THIS BENCH ITSELF. At TARGET = 1e-2 the measurement is
 * DEGENERATE: the cold solve settles at 0.0089-0.0099 because it stops the
 * instant it crosses the threshold, and the transferred warm start arrives at
 * 0.011-0.023 - just above the same threshold. "Iterations to 1%" then counts
 * the shaving of a factor 1.2, not convergence, and the count is threshold-
 * crossing noise: the teleport control returned ZERO iterations at the two
 * large domains purely because its start happened to land at 9.9e-3, which
 * would have read as "a teleport is free" - a negative control passing by
 * accident. A target is a measurement only if it sits well BELOW what the warm
 * start already reaches and well ABOVE the direct-solve floor. 2e-3 is 5.5x
 * below the warm start and 5.3x above the floor. BOTH are reported: the
 * degenerate one to show the trap, the working one to answer the question. */
static const double TARGET2 = 2e-3;
static const double DIRECT_MAX = 6e-4; /* bench validity: the direct solve must
                                        * clear TARGET2 by 3x, else no
                                        * iteration count means anything */
/* Slab contrast dk2/k2 and thickness. Weak on purpose: cuts are NOT exercised
 * here, and an uncut element straddling a strong boundary has a known ~20%
 * floor (M9b F5). Strong contrast is a different measurement.
 * HOW THIS NUMBER WAS SET, and it was set by the pre-registered remedy, not by
 * taste: dk2/k2 = 0.25 gave a direct-solve error of 3.0e-2, ten times over the
 * validity gate of 3e-3, entirely from uncut elements straddling the slab
 * faces (the same run in VACUUM gives 7e-8). The audit fixed the remedy in
 * advance - lower the contrast until the direct solve clears the gate, and say
 * so - and weak scattering is linear in the contrast, so 0.25 -> 0.02.
 * The price is stated, not hidden: with a weak scene the structure the basis
 * cannot represent exactly is only a few percent of the field, so struct_strength()
 * is printed next to every result. If it ever falls far below the 1% target we
 * are back in the vacuum trap R2 and the numbers mean nothing. */
static const double SLAB_DK2 = 0.02;
static const double SLAB_TH = 4.0; /* lambda/4 */
static const double SLAB_AT[MAXSLAB] = {0.30, 0.36, 0.42, 0.48, 0.80, 0.88};

static const double DOMS[4] = {200.0, 2000.0, 20000.0, 200000.0};
static const double CAM_AT = 0.60; /* camera start, fraction of the domain */
static const double SRC_AT = 0.25; /* source, fraction of the domain */

/* ------------------------------------------------------------------ medium */

/* Piecewise-constant k^2 over [0,dom]: vacuum with MAXSLAB thin slabs. Returns
 * the number of segments; segments are contiguous and cover the whole domain,
 * which is what hz_carrier_entry expects when it clips the medium integral. */
static int build_medium(double dom, double k, int with_slabs, hz_medseg *segs) {
  double complex k2v = CMPLX(k * k, 0.0);
  if (!with_slabs) {
    segs[0] = (hz_medseg){0.0, dom, k2v};
    return 1;
  }
  double complex k2s = CMPLX(k * k * (1.0 + SLAB_DK2), 0.0);
  int ns = 0;
  double x = 0.0;
  for (int s = 0; s < MAXSLAB; s++) {
    double a = SLAB_AT[s] * dom, b = a + SLAB_TH;
    segs[ns++] = (hz_medseg){x, a, k2v};
    segs[ns++] = (hz_medseg){a, b, k2s};
    x = b;
  }
  segs[ns++] = (hz_medseg){x, dom, k2v};
  return ns;
}

/* --------------------------------------------------------------- reference */

/* Exact layered Green's function. u = u_L(x_<) u_R(x_>) / W, with u_L outgoing
 * at the left end, u_R outgoing at the right, built by transfer matrices; W is
 * the (constant) Wronskian. Analytic at any lambda, no grid anywhere - A4 of
 * the M9 audit forbids a numerical reference in this regime.
 *
 * THE TWELFTH ARTEFACT OF THE LAST SESSION LIVED HERE: the reference had been
 * taken for a DELTA source while the source is phi(x - xs), and the missing
 * phihat(k) showed up as a constant 2.8481 error at every domain size. Hence
 * the three self-tests in ref_selftest(). */
typedef struct {
  int nreg;
  double xb[MAXSEG];                     /* xb[s] = right end of region s */
  double complex kap[MAXSEG];            /* wavenumber in region s */
  double complex AL[MAXSEG], BL[MAXSEG]; /* u_L = AL e^{i kap x} + BL e^{-i kap x} */
  double complex AR[MAXSEG], BR[MAXSEG];
  double complex wro;
  double complex ovl_L, ovl_R; /* Int u_L f dx, Int u_R f dx */
  double xs;
  double wro_spread; /* self-test: relative variation of W across regions */
} layref;

static double complex reg_val(double complex a, double complex b, double complex kap, double x) {
  double complex i = CMPLX(0.0, 1.0);
  return a * cexp(i * kap * x) + b * cexp(-i * kap * x);
}

/* Int_{-2}^{2} phi(t) e^{-i om t} dt, by Gauss-Legendre on the four unit cells.
 * phi is even so the result is real; ng is the rule order (8 and 16 are
 * compared in the self-test, they must agree to 1e-14). */
static double phihat(double om, int ng) {
  static const double G8X[8] = {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290,
                                -0.1834346424956498, 0.1834346424956498,  0.5255324099163290,
                                0.7966664774136267,  0.9602898564975363};
  static const double G8W[8] = {0.1012285362903763, 0.2223810344533745, 0.3137066458778873,
                                0.3626837833783620, 0.3626837833783620, 0.3137066458778873,
                                0.2223810344533745, 0.1012285362903763};
  static const double G16X[16] = {
      -0.9894009349916499, -0.9445750230732326, -0.8656312023878318, -0.7554044083550030,
      -0.6178762444026438, -0.4580167776572274, -0.2816035507792589, -0.0950125098376374,
      0.0950125098376374,  0.2816035507792589,  0.4580167776572274,  0.6178762444026438,
      0.7554044083550030,  0.8656312023878318,  0.9445750230732326,  0.9894009349916499};
  static const double G16W[16] = {
      0.0271524594117541, 0.0622535239386479, 0.0951585116824928, 0.1246289712555339,
      0.1495959888165767, 0.1691565193950025, 0.1826034150449236, 0.1894506104550685,
      0.1894506104550685, 0.1826034150449236, 0.1691565193950025, 0.1495959888165767,
      0.1246289712555339, 0.0951585116824928, 0.0622535239386479, 0.0271524594117541};
  const double *gx = (ng == 8) ? G8X : G16X;
  const double *gw = (ng == 8) ? G8W : G16W;
  double sum = 0.0;
  for (int p = -2; p < 2; p++)
    for (int g = 0; g < ng; g++) {
      double t = (double)p + 0.5 + 0.5 * gx[g];
      sum += 0.5 * gw[g] * hz_phi(t) * cos(om * t);
    }
  return sum;
}

static void ref_build(layref *r, const hz_medseg *segs, int nseg, double xs) {
  double complex i = CMPLX(0.0, 1.0);
  r->nreg = nseg;
  r->xs = xs;
  for (int s = 0; s < nseg; s++) {
    r->xb[s] = segs[s].b;
    r->kap[s] = csqrt(segs[s].k2);
  }
  /* u_L: e^{-i kap x} in region 0, propagated rightwards */
  r->AL[0] = 0.0;
  r->BL[0] = 1.0;
  for (int s = 0; s + 1 < nseg; s++) {
    double xi = r->xb[s];
    double complex k1 = r->kap[s], k2 = r->kap[s + 1];
    double complex p = r->AL[s] * cexp(i * k1 * xi), q = r->BL[s] * cexp(-i * k1 * xi);
    double complex u = p + q, du = k1 * (p - q);
    r->AL[s + 1] = cexp(-i * k2 * xi) * (u + du / k2) * 0.5;
    r->BL[s + 1] = cexp(i * k2 * xi) * (u - du / k2) * 0.5;
  }
  /* u_R: e^{+i kap x} in the last region, propagated leftwards */
  r->AR[nseg - 1] = 1.0;
  r->BR[nseg - 1] = 0.0;
  for (int s = nseg - 2; s >= 0; s--) {
    double xi = r->xb[s];
    double complex k1 = r->kap[s], k2 = r->kap[s + 1];
    double complex p = r->AR[s + 1] * cexp(i * k2 * xi), q = r->BR[s + 1] * cexp(-i * k2 * xi);
    double complex u = p + q, du = k2 * (p - q);
    r->AR[s] = cexp(-i * k1 * xi) * (u + du / k1) * 0.5;
    r->BR[s] = cexp(i * k1 * xi) * (u - du / k1) * 0.5;
  }
  /* W = u_L u_R' - u_L' u_R = 2 i kap (BL AR - AL BR), constant in x within a
   * region and constant ACROSS regions - which is self-test 2. */
  double lo = 1e300, hi = 0.0;
  for (int s = 0; s < nseg; s++) {
    double complex w = CMPLX(0.0, 2.0) * r->kap[s] * (r->BL[s] * r->AR[s] - r->AL[s] * r->BR[s]);
    if (s == 0) r->wro = w;
    double m = cabs(w);
    if (m < lo) lo = m;
    if (m > hi) hi = m;
  }
  r->wro_spread = (hi > 0.0) ? (hi - lo) / hi : 0.0;

  /* Int u f dx with f = phi(x - xs): the source support sits inside one
   * homogeneous region, so u = A e^{i kap x} + B e^{-i kap x} there and the
   * overlap is A e^{i kap xs} phihat(-kap) + B e^{-i kap xs} phihat(kap). */
  int s0 = 0;
  for (int s = 0; s < nseg; s++)
    if (xs >= segs[s].a && xs < segs[s].b) s0 = s;
  double kk = creal(r->kap[s0]);
  double ph = phihat(kk, 16); /* even in om, so phihat(+kap) = phihat(-kap) */
  r->ovl_L = ph * (r->AL[s0] * cexp(i * r->kap[s0] * xs) + r->BL[s0] * cexp(-i * r->kap[s0] * xs));
  r->ovl_R = ph * (r->AR[s0] * cexp(i * r->kap[s0] * xs) + r->BR[s0] * cexp(-i * r->kap[s0] * xs));
}

static double complex ref_u(const layref *r, double x) {
  int s = 0;
  for (int t = 0; t < r->nreg; t++)
    if (x >= r->xb[t]) s = t + 1;
  if (s >= r->nreg) s = r->nreg - 1;
  if (x >= r->xs) return reg_val(r->AR[s], r->BR[s], r->kap[s], x) * r->ovl_L / r->wro;
  return reg_val(r->AL[s], r->BL[s], r->kap[s], x) * r->ovl_R / r->wro;
}

/* Self-tests of the reference itself. Without these the whole bench measures
 * the reference, as it did once already. */
static int ref_selftest(double k) {
  int ok = 1;
  double e8 = phihat(k, 8), e16 = phihat(k, 16);
  double dq = fabs(e8 - e16) / fabs(e16);
  printf("  [ref 3/3] phihat GL8 vs GL16   rel %8.1e   %s\n", dq, dq < 1e-14 ? "PASS" : "FAIL");
  if (!(dq < 1e-14)) ok = 0;

  /* Test 1 is run at FOUR domain sizes on purpose. A wrong reference (the
   * source taken as a delta) showed up last session as an error that was
   * CONSTANT across domain sizes; pure phase rounding, in contrast, must grow
   * like eps*k*dom, because the transfer matrices carry e^{i k x} at absolute x
   * while the closed form carries e^{i k |x - xs|}, and the large phase is
   * cancelled between two independently rounded numbers. So the test is on the
   * error DIVIDED by eps*k*dom: a bug makes that ratio explode as the domain
   * grows, rounding keeps it O(1). */
  double complex amp = e16 / CMPLX(0.0, 2.0 * k);
  double worst_norm = 0.0;
  hz_medseg segs[MAXSEG];
  layref r;
  for (int id = 0; id < 4; id++) {
    double dom = DOMS[id] * LAM, xs = SRC_AT * dom;
    int nseg = build_medium(dom, k, 0, segs);
    ref_build(&r, segs, nseg, xs);
    double worst = 0.0;
    for (int p = 0; p < 500; p++) {
      double x = 0.05 * dom + 0.9 * dom * ((double)p + 0.5) / 500.0;
      if (fabs(x - xs) < 4.0) continue;
      double complex u = amp * cexp(CMPLX(0.0, 1.0) * k * fabs(x - xs));
      double d = cabs(ref_u(&r, x) - u) / cabs(u);
      if (d > worst) worst = d;
    }
    double nrm = worst / (2.2e-16 * k * dom);
    if (nrm > worst_norm) worst_norm = nrm;
    printf("  [ref 1/3] no slabs vs Green    dom %7.0f lam  rel %8.1e  / (eps k dom) %6.2f\n",
           DOMS[id], worst, nrm);
  }
  printf("            verdict: %s (rounding-limited iff the ratio stays O(1))\n",
         worst_norm < 10.0 ? "PASS" : "FAIL");
  if (!(worst_norm < 10.0)) ok = 0;

  double dom = 200.0 * LAM, xs = SRC_AT * dom;
  int nseg = build_medium(dom, k, 1, segs);
  ref_build(&r, segs, nseg, xs);
  printf("  [ref 2/3] Wronskian spread     rel %8.1e   %s (%d regions)\n", r.wro_spread,
         r.wro_spread < 1e-12 ? "PASS" : "FAIL", nseg);
  if (!(r.wro_spread < 1e-12)) ok = 0;
  return ok;
}

/* -------------------------------------------------------------------- basis */

/* An element is identified by (level, node, direction) and NOT by its width as
 * a float: the node lattice is fixed in space (xc = n*W), so two bases taken at
 * different camera positions share elements EXACTLY, and a surviving element
 * keeps its coefficient verbatim. That exactness is the whole mechanism of the
 * increment - if elements drifted in size with R, nothing would transfer. */
typedef struct {
  int j;
  int n;
  int dir;
} bkey;

static int key_find(const bkey *set, int n, bkey q) {
  for (int i = 0; i < n; i++)
    if (set[i].j == q.j && set[i].n == q.n && set[i].dir == q.dir) return i;
  return -1;
}

/* Shells anchored at xa: level j covers R in [W_j/eps, 2 W_j/eps), widened by
 * MARGIN elements so the partition of unity holds across the shell rather than
 * dipping at its edges. jmax < 0 means "the full ladder up to WMAX_FRAC*dom,
 * with the finest owning everything nearer and the coarsest everything beyond";
 * jmax >= 0 means a bounded patch (used for the source). */
static int add_shells(bkey *set, int cap, int nhave, double xa, double dom, int jmax) {
  int d = nhave;
  for (int j = 0; j < NLEVMAX; j++) {
    double W = W0 * pow(2.0, (double)j);
    if (jmax < 0) {
      if (W > WMAX_FRAC * dom) break;
    } else if (j > jmax) {
      break;
    }
    double alo = (j == 0) ? 0.0 : W / EPS, ahi = 2.0 * W / EPS;
    if (jmax < 0 && W0 * pow(2.0, (double)(j + 1)) > WMAX_FRAC * dom) ahi = 2.0 * dom;
    alo -= (double)MARGIN * W;
    ahi += (double)MARGIN * W;
    if (alo < 0.0) alo = 0.0;
    int n0 = (int)((xa - ahi) / W) - 2, n1 = (int)((xa + ahi) / W) + 2;
    for (int n = n0; n <= n1; n++) {
      double xc = (double)n * W, dist = fabs(xc - xa);
      if (xc < -2.0 * W || xc > dom + 2.0 * W) continue;
      if (dist < alo || dist >= ahi) continue;
      for (int s = 0; s < 2; s++) {
        bkey q = {j, n, s == 0 ? 1 : -1};
        if (key_find(set, d, q) >= 0)
          continue; /* dedupe: duplicate columns
                     * would add degeneracy that has
                     * nothing to do with the case */
        if (d >= cap) return d;
        set[d++] = q;
      }
    }
  }
  return d;
}

/* camera shells (LOD ceiling) + a camera-INDEPENDENT patch at the source
 * (geometry floor, regulator 2 of the plan). Without the patch a distant camera
 * puts the emitter inside an element thousands of lambda wide and the bench
 * measures the failure to represent a source, not the increment. */
static int build_basis(bkey *set, int cap, double xcam, double xs, double dom) {
  int d = add_shells(set, cap, 0, xcam, dom, -1);
  d = add_shells(set, cap, d, xs, dom, SRCLEV);
  return d;
}

static void keys_to_carriers(const bkey *set, int n, double k, hz_carrier *b) {
  for (int i = 0; i < n; i++) {
    double W = W0 * pow(2.0, (double)set[i].j);
    b[i] = (hz_carrier){W, set[i].n, CMPLX((double)set[i].dir * k, 0.0), -1};
  }
}

/* ----------------------------------------------------------------- assembly */

static void assemble(const hz_carrier *b, int dim, double dom, const hz_medseg *segs, int nseg,
                     double xs, double scl, double k, double complex *A, double complex *rhs) {
  int i;
#pragma omp parallel for schedule(dynamic, 8) private(i)
  for (i = 0; i < dim; i++) {
    for (int j = 0; j < dim; j++)
      A[(size_t)i * (size_t)dim + (size_t)j] =
          hz_carrier_entry(b[i], b[j], dom, 0.0, segs, nseg, NULL);
    rhs[i] = hz_carrier_rhs(b[i], dom, 0.0, xs, 1.0, NULL);
  }
  /* DtN rows: u' = -i k u at 0, u' = +i k u at dom. Two rows, and they are the
   * whole radiation condition (TERMINATION_REPORT: DtN works at any element
   * size, absorption does not). scl is FIXED at frame 0 - see R5. */
  double complex ii = CMPLX(0.0, 1.0);
  for (int e = 0; e < 2; e++) {
    double xb = (e == 0) ? 0.0 : dom;
    double complex sgn = (e == 0) ? -ii * k : ii * k;
    for (int j = 0; j < dim; j++) {
      double t = xb / b[j].W - (double)b[j].n;
      double complex v = 0.0, d = 0.0;
      if (fabs(t) < 2.0) {
        double complex ex = cexp(ii * b[j].kx * (xb - b[j].W * (double)b[j].n));
        v = hz_phi(t) * ex;
        d = (hz_phi_d1(t) / b[j].W + ii * b[j].kx * hz_phi(t)) * ex;
      }
      A[(size_t)(dim + e) * (size_t)dim + (size_t)j] = scl * (d - sgn * v);
    }
    rhs[dim + e] = 0.0;
  }
}

/* the scaling recipe of carrier_solve.c, kept verbatim for continuity, but
 * evaluated ONCE and then frozen */
static double dtn_scale(const double complex *A, int dim) {
  double s = 0.0;
  for (int j = 0; j < dim; j++)
    s += cabs(A[(size_t)(dim / 2) * (size_t)dim + (size_t)j]);
  return (s > 0.0) ? s : 1.0;
}

/* --------------------------------------------------------------- field/LSQR */

typedef struct {
  int np;
  double x[NPROBEG];
  double complex u[NPROBEG];
  double den;
} probes;

static void probes_fill(probes *p, const layref *r, double lo, double hi, int np, double xs) {
  p->np = 0;
  p->den = 0.0;
  for (int q = 0; q < np; q++) {
    double x = lo + (hi - lo) * ((double)q + 0.5) / (double)np;
    if (fabs(x - xs) < 4.0)
      continue; /* inside the source support the field is
                 * not a pure outgoing wave */
    p->x[p->np] = x;
    p->u[p->np] = ref_u(r, x);
    p->den += creal(p->u[p->np] * conj(p->u[p->np]));
    p->np++;
  }
}

/* P[q][j] = B_j(x_q): the field at the probes is one dense mat-vec, so it can
 * be checked EVERY iteration (carrier_iter checked every 5th, too coarse for
 * warm counts of a few) */
static void probe_matrix(const probes *p, const hz_carrier *b, int dim, double complex *P) {
  int q;
#pragma omp parallel for schedule(static) private(q)
  for (q = 0; q < p->np; q++)
    for (int j = 0; j < dim; j++)
      P[(size_t)q * (size_t)dim + (size_t)j] = hz_carrier_val(b[j], p->x[q], 0.0);
}

static double field_err(const probes *p, const double complex *P, int dim,
                        const double complex *x) {
  double num = 0.0;
  for (int q = 0; q < p->np; q++) {
    double complex v = 0.0;
    for (int j = 0; j < dim; j++)
      v += P[(size_t)q * (size_t)dim + (size_t)j] * x[j];
    num += creal((v - p->u[q]) * conj(v - p->u[q]));
  }
  return sqrt(num / p->den);
}

typedef struct {
  double complex *u, *v, *w, *d, *tmp, *xw;
} lsqr_ws;

/* one owner for the six work vectors: cppcheck was right that the scattered
 * "return 0" paths leaked them */
static int ws_alloc(lsqr_ws *ws, size_t nrow) {
  size_t nb = (nrow > (size_t)MAXB ? nrow : (size_t)MAXB) + (size_t)MAXB;
  ws->u = calloc(nb, sizeof(double complex));
  ws->v = calloc((size_t)MAXB, sizeof(double complex));
  ws->w = calloc((size_t)MAXB, sizeof(double complex));
  ws->d = calloc((size_t)MAXB, sizeof(double complex));
  ws->tmp = calloc(nb, sizeof(double complex));
  ws->xw = calloc((size_t)MAXB, sizeof(double complex));
  return ws->u != NULL && ws->v != NULL && ws->w != NULL && ws->d != NULL && ws->tmp != NULL &&
         ws->xw != NULL;
}

static void ws_free(lsqr_ws *ws) {
  free(ws->u);
  free(ws->v);
  free(ws->w);
  free(ws->d);
  free(ws->tmp);
  free(ws->xw);
  ws->u = NULL;
  ws->v = NULL;
  ws->w = NULL;
  ws->d = NULL;
  ws->tmp = NULL;
  ws->xw = NULL;
}

/* LSQR (Paige-Saunders) on A dx = rhs - A x0, then x = x0 + dx. Algebraically a
 * warm start; LSQR itself takes no initial guess, so the correction system is
 * the honest way to give it one.
 * target > 0: stop at the first iteration whose FIELD error is at or below it,
 * return the count. target <= 0: run exactly maxit iterations (drift test). */
static int lsqr_solve(const double complex *A, int nrow, int dim, const double complex *rhs,
                      const double complex *x0, const probes *p, const double complex *P,
                      double target, int maxit, lsqr_ws *ws, double *err_out, double complex *x_out,
                      double *err0_out) {
  double complex *u = ws->u, *v = ws->v, *w = ws->w, *d = ws->d, *tmp = ws->tmp, *xw = ws->xw;
  for (int i = 0; i < nrow; i++) {
    double complex s = rhs[i];
    if (x0 != NULL)
      for (int j = 0; j < dim; j++)
        s -= A[(size_t)i * (size_t)dim + (size_t)j] * x0[j];
    u[i] = s;
  }
  double beta = 0.0;
  for (int i = 0; i < nrow; i++)
    beta += creal(u[i] * conj(u[i]));
  beta = sqrt(beta);
  for (int j = 0; j < dim; j++) {
    d[j] = 0.0;
    xw[j] = (x0 != NULL) ? x0[j] : 0.0;
  }
  /* The field error of the STARTING point. For a warm start this is the number
   * that explains everything: LSQR needs roughly rate * log(e0/target)
   * iterations, so a warm start buys the LOG factor and nothing else - the rate
   * is a property of the operator and the basis, not of the starting point. */
  double e0 = field_err(p, P, dim, xw);
  if (err0_out != NULL) *err0_out = e0;
  if (err_out != NULL) *err_out = e0;
  if (x_out != NULL)
    for (int j = 0; j < dim; j++)
      x_out[j] = xw[j];
  if (!(beta > 0.0)) return 0;
  if (target > 0.0 && e0 <= target) return 0;
  for (int i = 0; i < nrow; i++)
    u[i] /= beta;
  double alpha = 0.0;
  for (int j = 0; j < dim; j++) {
    double complex s = 0.0;
    for (int i = 0; i < nrow; i++)
      s += conj(A[(size_t)i * (size_t)dim + (size_t)j]) * u[i];
    v[j] = s;
    alpha += creal(s * conj(s));
  }
  alpha = sqrt(alpha);
  if (!(alpha > 0.0)) return 0;
  for (int j = 0; j < dim; j++) {
    v[j] /= alpha;
    w[j] = v[j];
  }
  double phibar = beta, rhobar = alpha;

  for (int it = 1; it <= maxit; it++) {
    int i;
#pragma omp parallel for schedule(static) private(i)
    for (i = 0; i < nrow; i++) {
      double complex s = 0.0;
      for (int j = 0; j < dim; j++)
        s += A[(size_t)i * (size_t)dim + (size_t)j] * v[j];
      tmp[i] = s - alpha * u[i];
    }
    double bnew = 0.0;
    for (i = 0; i < nrow; i++)
      bnew += creal(tmp[i] * conj(tmp[i]));
    bnew = sqrt(bnew);
    if (bnew > 0.0)
      for (i = 0; i < nrow; i++)
        u[i] = tmp[i] / bnew;
    int j;
#pragma omp parallel for schedule(static) private(j)
    for (j = 0; j < dim; j++) {
      double complex s = 0.0;
      for (int t = 0; t < nrow; t++)
        s += conj(A[(size_t)t * (size_t)dim + (size_t)j]) * u[t];
      tmp[j] = s - bnew * v[j];
    }
    double anew = 0.0;
    for (j = 0; j < dim; j++)
      anew += creal(tmp[j] * conj(tmp[j]));
    anew = sqrt(anew);
    if (anew > 0.0)
      for (j = 0; j < dim; j++)
        v[j] = tmp[j] / anew;

    double rho = sqrt(rhobar * rhobar + bnew * bnew);
    double c = rhobar / rho, s2 = bnew / rho;
    double theta = s2 * anew;
    rhobar = -c * anew;
    double phi = c * phibar;
    phibar = s2 * phibar;
    for (j = 0; j < dim; j++) {
      d[j] += (phi / rho) * w[j];
      w[j] = v[j] - (theta / rho) * w[j];
    }
    alpha = anew;

    for (j = 0; j < dim; j++)
      xw[j] = ((x0 != NULL) ? x0[j] : 0.0) + d[j];
    double e = field_err(p, P, dim, xw);
    if (err_out != NULL) *err_out = e;
    if (x_out != NULL)
      for (j = 0; j < dim; j++)
        x_out[j] = xw[j];
    if (target > 0.0 && e <= target) return it;
  }
  return -1;
}

/* --------------------------------------------------------------- experiment */

/* Everything needed to pose the problem at one camera position. */
typedef struct {
  bkey key[MAXB];
  hz_carrier b[MAXB];
  int dim;
  double complex *A;
  double complex *rhs;
} sysview;

static int sys_build(sysview *s, double xcam, double xs, double dom, const hz_medseg *segs,
                     int nseg, double k, double scl) {
  s->dim = build_basis(s->key, MAXB, xcam, xs, dom);
  if (s->dim <= 0) return 0;
  keys_to_carriers(s->key, s->dim, k, s->b);
  int nrow = s->dim + 2;
  s->A = calloc((size_t)nrow * (size_t)s->dim, sizeof(double complex));
  s->rhs = calloc((size_t)nrow, sizeof(double complex));
  if (s->A == NULL || s->rhs == NULL) return 0;
  /* scl < 0 asks for the unscaled matrix so that dtn_scale() can be read off it
   * once; the caller then applies the frozen value with sys_scale_dtn(). */
  assemble(s->b, s->dim, dom, segs, nseg, xs, scl < 0.0 ? 1.0 : scl, k, s->A, s->rhs);
  return 1;
}

static void sys_scale_dtn(sysview *s, double scl) {
  for (int e = 0; e < 2; e++)
    for (int j = 0; j < s->dim; j++)
      s->A[(size_t)(s->dim + e) * (size_t)s->dim + (size_t)j] *= scl;
}

static void sys_free(sysview *s) {
  free(s->A);
  free(s->rhs);
  s->A = NULL;
  s->rhs = NULL;
}

/* Copy the coefficients of the elements that survived; new elements start at
 * zero; the content of the elements that LEFT is honestly lost. There is no
 * exact fine->coarse restriction (refinability only runs the other way), so the
 * plan's phrase "departing elements are removed by exact restriction" does not
 * apply and is not used. */
static void transfer(const sysview *old, const double complex *xold, const sysview *nw,
                     double complex *xnew, int *nkept) {
  int kept = 0;
  for (int i = 0; i < nw->dim; i++) {
    int q = key_find(old->key, old->dim, nw->key[i]);
    if (q >= 0) {
      xnew[i] = xold[q];
      kept++;
    } else {
      xnew[i] = 0.0;
    }
  }
  *nkept = kept;
}

/* Falsifier 1: is the starting residual carried by the rows of elements that
 * were touched by the camera move? "Touched" = the element is new, or its
 * support meets the support of an element that appeared or disappeared. Both
 * the residual share AND the concentration ratio are reported: a share of 0.75
 * means nothing if 90% of the rows are touched. */
static void residual_locality(const sysview *old, const sysview *nw, const double complex *x0,
                              double *share, double *touched_frac) {
  static double chlo[MAXB * 2], chhi[MAXB * 2];
  int nch = 0;
  for (int i = 0; i < nw->dim && nch < MAXB * 2; i++)
    if (key_find(old->key, old->dim, nw->key[i]) < 0) {
      chlo[nch] = nw->b[i].W * ((double)nw->b[i].n - 2.0);
      chhi[nch] = nw->b[i].W * ((double)nw->b[i].n + 2.0);
      nch++;
    }
  for (int i = 0; i < old->dim && nch < MAXB * 2; i++)
    if (key_find(nw->key, nw->dim, old->key[i]) < 0) {
      chlo[nch] = old->b[i].W * ((double)old->b[i].n - 2.0);
      chhi[nch] = old->b[i].W * ((double)old->b[i].n + 2.0);
      nch++;
    }

  double tot = 0.0, hit = 0.0;
  int nhit = 0;
  for (int i = 0; i < nw->dim; i++) {
    double complex s = nw->rhs[i];
    for (int j = 0; j < nw->dim; j++)
      s -= nw->A[(size_t)i * (size_t)nw->dim + (size_t)j] * x0[j];
    double m = creal(s * conj(s));
    tot += m;
    double lo = nw->b[i].W * ((double)nw->b[i].n - 2.0);
    double hi = nw->b[i].W * ((double)nw->b[i].n + 2.0);
    int touched = 0;
    for (int c = 0; c < nch; c++)
      if (lo < chhi[c] && hi > chlo[c]) {
        touched = 1;
        break;
      }
    if (touched) {
      hit += m;
      nhit++;
    }
  }
  *share = (tot > 0.0) ? sqrt(hit / tot) : 0.0;
  *touched_frac = (nw->dim > 0) ? (double)nhit / (double)nw->dim : 0.0;
}

/* How much of the field is NOT a plain outgoing wave: ||u_scene - u_vacuum|| /
 * ||u_vacuum|| over the probe window. This is the number that says whether the
 * bench escaped the vacuum trap (audit R2). In 1D the only sources of structure
 * are the medium and the source kink - any field in a homogeneous region is
 * A e^{ikx} + B e^{-ikx} with CONSTANT coefficients, i.e. exactly in the span
 * at any element width. So this quantity is the whole margin the measurement
 * has, and it belongs next to every result rather than in a footnote. */
static double struct_strength(const probes *p, const layref *vac) {
  double num = 0.0, den = 0.0;
  for (int q = 0; q < p->np; q++) {
    double complex uv = ref_u(vac, p->x[q]);
    num += creal((p->u[q] - uv) * conj(p->u[q] - uv));
    den += creal(uv * conj(uv));
  }
  return (den > 0.0) ? sqrt(num / den) : 0.0;
}

static double direct_err(const sysview *s, const probes *p, const double complex *P) {
  int nrow = s->dim + 2;
  double complex *Aw = calloc((size_t)nrow * (size_t)s->dim, sizeof(double complex));
  double complex *rw = calloc((size_t)nrow, sizeof(double complex));
  double *sv = calloc((size_t)s->dim, sizeof(double));
  if (Aw == NULL || rw == NULL || sv == NULL) {
    free(Aw);
    free(rw);
    free(sv);
    return -1.0;
  }
  for (size_t t = 0; t < (size_t)nrow * (size_t)s->dim; t++)
    Aw[t] = s->A[t];
  for (int t = 0; t < nrow; t++)
    rw[t] = s->rhs[t];
  lapack_int rank = 0;
  double e = -1.0;
  if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, nrow, s->dim, 1, Aw, s->dim, rw, 1, sv, RCOND, &rank) == 0)
    e = field_err(p, P, s->dim, rw);
  free(Aw);
  free(rw);
  free(sv);
  return e;
}

/* deterministic scramble for NC2: same norm as the transferred start, no
 * information about the solution */
static void scramble(double complex *x, int dim, const double complex *ref) {
  double nrm = 0.0;
  for (int j = 0; j < dim; j++)
    nrm += creal(ref[j] * conj(ref[j]));
  nrm = sqrt(nrm);
  unsigned long st = 88172645463325252UL;
  double acc = 0.0;
  for (int j = 0; j < dim; j++) {
    st ^= st << 13;
    st ^= st >> 7;
    st ^= st << 17;
    double a = (double)(st >> 11) / 9007199254740992.0 - 0.5;
    st ^= st << 13;
    st ^= st >> 7;
    st ^= st << 17;
    double b = (double)(st >> 11) / 9007199254740992.0 - 0.5;
    x[j] = CMPLX(a, b);
    acc += a * a + b * b;
  }
  acc = sqrt(acc);
  if (acc > 0.0)
    for (int j = 0; j < dim; j++)
      x[j] *= nrm / acc;
}

/* shifts in units of R_near = W0/EPS, the inner radius of the first shell */
static const double SHIFTS[6] = {0.0, 0.0625, 0.125, 0.25, 0.5, 2.0};

static int experiment_a(double k, int with_slabs) {
  double rnear = W0 / EPS;
  printf("\n=== EXPERIMENT A: one camera step, %s ===\n",
         with_slabs ? "SCENE WITH SLABS (the result)" : "VACUUM (instrument floor, see R2)");
  printf("  %8s %6s %9s %6s %6s | %7s %6s %6s %6s %7s %9s %5s %5s\n", "dom/lam", "dim", "direct",
         "cold1%", "cold.2%", "d/Rnear", "chg", "resid", "touch", "conc", "err0", "w1%", "w.2%");

  for (int id = 0; id < 4; id++) {
    double dom = DOMS[id] * LAM, xs = SRC_AT * dom, xcam = CAM_AT * dom;
    hz_medseg segs[MAXSEG], vsegs[MAXSEG];
    int nseg = build_medium(dom, k, with_slabs, segs);
    layref ref, vref;
    ref_build(&ref, segs, nseg, xs);
    int vnseg = build_medium(dom, k, 0, vsegs); /* separate statement: passing
                                                 * vsegs and filling it in the same argument list is
                                                 * unsequenced */
    ref_build(&vref, vsegs, vnseg, xs);

    static sysview s0, s1;
    if (!sys_build(&s0, xcam, xs, dom, segs, nseg, k, -1.0)) return 0;
    double scl = dtn_scale(s0.A, s0.dim);
    sys_scale_dtn(&s0, scl);

    static probes pc;
    probes_fill(&pc, &ref, xcam - 4.0 * LAM, xcam + 4.0 * LAM, NPROBE, xs);
    double complex *P0 = calloc((size_t)pc.np * (size_t)s0.dim, sizeof(double complex));
    if (P0 == NULL) return 0;
    probe_matrix(&pc, s0.b, s0.dim, P0);

    /* R4: the camera window is where the basis is finest BY CONSTRUCTION, so a
     * global window is measured next to it and both are reported. */
    static probes pg;
    probes_fill(&pg, &ref, 0.05 * dom, 0.95 * dom, NPROBEG, xs);
    double complex *PG = calloc((size_t)pg.np * (size_t)s0.dim, sizeof(double complex));
    if (PG == NULL) return 0;
    probe_matrix(&pg, s0.b, s0.dim, PG);

    double de = direct_err(&s0, &pc, P0);
    double deg = direct_err(&s0, &pg, PG);
    int nrow0 = s0.dim + 2;
    lsqr_ws ws;
    static double complex xcold[MAXB], xwarm[MAXB];
    if (!ws_alloc(&ws, (size_t)nrow0)) {
      ws_free(&ws);
      return 0;
    }
    double ecold = 0.0;
    int cold = lsqr_solve(s0.A, nrow0, s0.dim, s0.rhs, NULL, &pc, P0, TARGET, MAXIT, &ws, &ecold,
                          xcold, NULL);
    int cold2 = lsqr_solve(s0.A, nrow0, s0.dim, s0.rhs, NULL, &pc, P0, TARGET2, MAXIT, &ws, &ecold,
                           xcold, NULL);
    double ecoldg = field_err(&pg, PG, s0.dim, xcold);
    printf("  %8.0f %6d %9.2e %6d %6d |  R_near %.0f  direct(global) %8.2e  cold(global) %8.2e  "
           "structure %5.2f%%\n",
           DOMS[id], s0.dim, de, cold, cold2, rnear, deg, ecoldg,
           100.0 * struct_strength(&pc, &vref));
    if (!(de >= 0.0) || de > DIRECT_MAX)
      printf("      !! BENCH INVALID: direct error %.2e exceeds %.1e, the 1%% target is not\n"
             "         reachable in this basis and no iteration count below means anything\n",
             de, DIRECT_MAX);
    fflush(stdout);

    for (int is = 0; is < 6; is++) {
      double delta = SHIFTS[is] * rnear;
      if (!sys_build(&s1, xcam + delta, xs, dom, segs, nseg, k, scl)) return 0;
      int kept = 0;
      transfer(&s0, xcold, &s1, xwarm, &kept);
      int chg = (s0.dim - kept) + (s1.dim - kept);
      double share = 0.0, tf = 0.0;
      residual_locality(&s0, &s1, xwarm, &share, &tf);
      static probes pc1;
      probes_fill(&pc1, &ref, xcam + delta - 4.0 * LAM, xcam + delta + 4.0 * LAM, NPROBE, xs);
      double complex *P1 = calloc((size_t)pc1.np * (size_t)s1.dim, sizeof(double complex));
      if (P1 == NULL) return 0;
      probe_matrix(&pc1, s1.b, s1.dim, P1);
      double ew = 0.0;
      double e0 = 0.0;
      int warm = lsqr_solve(s1.A, s1.dim + 2, s1.dim, s1.rhs, xwarm, &pc1, P1, TARGET, MAXIT, &ws,
                            &ew, NULL, &e0);
      int warm2 = lsqr_solve(s1.A, s1.dim + 2, s1.dim, s1.rhs, xwarm, &pc1, P1, TARGET2, MAXIT, &ws,
                             &ew, NULL, NULL);
      printf("  %8s %6s %9s %6s %6s | %7.4f %6.3f %6.2f %6.2f %7.2f %9.2e %5d %5d\n", "", "", "",
             "", "", SHIFTS[is], (double)chg / (double)s1.dim, share, tf,
             tf > 0.0 ? share / tf : 0.0, e0, warm, warm2);
      fflush(stdout);
      free(P1);
      sys_free(&s1);
    }

    /* --- negative controls, failure PREDICTED (see M11_AUDIT) --- */
    double dnc = 0.5 * rnear;
    /* NC1: teleport. predicted warm ~ cold */
    {
      double delta = 0.25 * dom;
      if (!sys_build(&s1, xcam + delta, xs, dom, segs, nseg, k, scl)) return 0;
      int kept = 0;
      transfer(&s0, xcold, &s1, xwarm, &kept);
      static probes pt;
      probes_fill(&pt, &ref, xcam + delta - 4.0 * LAM, xcam + delta + 4.0 * LAM, NPROBE, xs);
      double complex *P1 = calloc((size_t)pt.np * (size_t)s1.dim, sizeof(double complex));
      if (P1 == NULL) return 0;
      probe_matrix(&pt, s1.b, s1.dim, P1);
      double ew = 0.0;
      double e0 = 0.0;
      int warm = lsqr_solve(s1.A, s1.dim + 2, s1.dim, s1.rhs, xwarm, &pt, P1, TARGET2, MAXIT, &ws,
                            &ew, NULL, &e0);
      int cw = lsqr_solve(s1.A, s1.dim + 2, s1.dim, s1.rhs, NULL, &pt, P1, TARGET2, MAXIT, &ws, &ew,
                          NULL, NULL);
      printf("      NC1 teleport d=0.25*dom   kept %5.3f  warm %4d  cold %4d  err0 %8.2e  "
             "(predict warm~cold)\n",
             (double)kept / (double)s1.dim, warm, cw, e0);
      fflush(stdout);
      free(P1);
      sys_free(&s1);
    }
    /* NC2: scrambled warm start at the same small shift. predicted ~cold */
    {
      if (!sys_build(&s1, xcam + dnc, xs, dom, segs, nseg, k, scl)) return 0;
      int kept = 0;
      transfer(&s0, xcold, &s1, xwarm, &kept);
      static double complex xscr[MAXB];
      scramble(xscr, s1.dim, xwarm);
      static probes ps;
      probes_fill(&ps, &ref, xcam + dnc - 4.0 * LAM, xcam + dnc + 4.0 * LAM, NPROBE, xs);
      double complex *P1 = calloc((size_t)ps.np * (size_t)s1.dim, sizeof(double complex));
      if (P1 == NULL) return 0;
      probe_matrix(&ps, s1.b, s1.dim, P1);
      double ew = 0.0;
      double e0 = 0.0;
      int scr = lsqr_solve(s1.A, s1.dim + 2, s1.dim, s1.rhs, xscr, &ps, P1, TARGET2, MAXIT, &ws,
                           &ew, NULL, &e0);
      printf("      NC2 scrambled start       d/Rnear 0.5   iters %4d  err0 %8.2e  "
             "(predict ~cold2 %d)\n",
             scr, e0, cold2);
      fflush(stdout);
      free(P1);
      sys_free(&s1);
    }
    /* NC3: move the SOURCE by the same delta, camera fixed. The field really
     * changes; predicted much more expensive than moving the camera. */
    {
      double xs2 = xs + dnc;
      layref ref2;
      ref_build(&ref2, segs, nseg, xs2);
      if (!sys_build(&s1, xcam, xs2, dom, segs, nseg, k, scl)) return 0;
      int kept = 0;
      transfer(&s0, xcold, &s1, xwarm, &kept);
      static probes p3;
      probes_fill(&p3, &ref2, xcam - 4.0 * LAM, xcam + 4.0 * LAM, NPROBE, xs2);
      double complex *P1 = calloc((size_t)p3.np * (size_t)s1.dim, sizeof(double complex));
      if (P1 == NULL) return 0;
      probe_matrix(&p3, s1.b, s1.dim, P1);
      double ew = 0.0;
      double e0 = 0.0;
      int mv = lsqr_solve(s1.A, s1.dim + 2, s1.dim, s1.rhs, xwarm, &p3, P1, TARGET2, MAXIT, &ws,
                          &ew, NULL, &e0);
      printf("      NC3 source moves instead  d/Rnear 0.5   iters %4d  kept %5.3f  err0 %8.2e  "
             "(predict >> camera)\n",
             mv, (double)kept / (double)s1.dim, e0);
      fflush(stdout);
      free(P1);
      sys_free(&s1);
    }

    free(P0);
    free(PG);
    ws_free(&ws);
    sys_free(&s0);
  }
  return 1;
}

/* Falsifier 4 (added by the audit, R6): the camera walks NFRAME frames, each
 * frame gets a FIXED small budget of iterations, and we watch whether the field
 * error stays put or drifts. A scheme that looks perfect for one step and
 * accumulates over a minute is not a real-time scheme. */
static int experiment_b(double k, double domlam, int budget) {
  double dom = domlam * LAM, xs = SRC_AT * dom, xcam = CAM_AT * dom;
  double step = 0.25 * (W0 / EPS);
  hz_medseg segs[MAXSEG];
  int nseg = build_medium(dom, k, 1, segs);
  layref ref;
  ref_build(&ref, segs, nseg, xs);

  printf(
      "\n=== EXPERIMENT B: drift, %.0f lambda, %d frames, %d iters/frame, step %.2f R_near ===\n",
      domlam, NFRAME, budget, 0.25);

  static sysview s0, s1;
  if (!sys_build(&s0, xcam, xs, dom, segs, nseg, k, -1.0)) return 0;
  double scl = dtn_scale(s0.A, s0.dim);
  sys_scale_dtn(&s0, scl);

  static probes pc;
  probes_fill(&pc, &ref, xcam - 4.0 * LAM, xcam + 4.0 * LAM, NPROBE, xs);
  double complex *P = calloc((size_t)pc.np * (size_t)MAXB, sizeof(double complex));
  if (P == NULL) return 0;
  probe_matrix(&pc, s0.b, s0.dim, P);

  lsqr_ws ws;
  static double complex xcur[MAXB], xnext[MAXB];
  if (!ws_alloc(&ws, (size_t)MAXB)) {
    ws_free(&ws);
    return 0;
  }

  double e = 0.0;
  int cold = lsqr_solve(s0.A, s0.dim + 2, s0.dim, s0.rhs, NULL, &pc, P, TARGET, MAXIT, &ws, &e,
                        xcur, NULL);
  double floor0 = direct_err(&s0, &pc, P);
  printf("  frame 0: cold %d iters, err %.3e, direct floor %.3e\n", cold, e, floor0);
  printf("  %6s %8s %9s %9s %7s\n", "frame", "x/dom", "err", "floor", "chg");
  fflush(stdout);

  for (int f = 1; f <= NFRAME; f++) {
    double xc = xcam + step * (double)f;
    if (!sys_build(&s1, xc, xs, dom, segs, nseg, k, scl)) return 0;
    int kept = 0;
    transfer(&s0, xcur, &s1, xnext, &kept);
    int chg = (s0.dim - kept) + (s1.dim - kept);
    probes_fill(&pc, &ref, xc - 4.0 * LAM, xc + 4.0 * LAM, NPROBE, xs);
    probe_matrix(&pc, s1.b, s1.dim, P);
    lsqr_solve(s1.A, s1.dim + 2, s1.dim, s1.rhs, xnext, &pc, P, -1.0, budget, &ws, &e, xcur, NULL);
    double fl = -1.0;
    if (f % DRIFT_DIRECT == 0 || f == NFRAME) fl = direct_err(&s1, &pc, P);
    if (f % 4 == 0 || f <= 2) {
      printf("  %6d %8.4f %9.3e", f, xc / dom, e);
      if (fl >= 0.0)
        printf(" %9.3e", fl);
      else
        printf(" %9s", "-");
      printf(" %7.3f\n", (double)chg / (double)s1.dim);
      fflush(stdout);
    }
    sys_free(&s0);
    s0 = s1; /* struct copy takes the arrays and the heap pointers */
    s1.A = NULL;
    s1.rhs = NULL;
  }
  free(P);
  ws_free(&ws);
  sys_free(&s0);
  return 1;
}

int main(int argc, char **argv) {
  int mode = 0; /* 0 = everything, 1 = vacuum only, 2 = slabs only */
  if (argc > 1 && argv[1][0] == 'v') mode = 1;
  if (argc > 1 && argv[1][0] == 's') mode = 2;
  double k = 2.0 * M_PI / LAM;
  printf("M11: correction iterations when the camera moves.\n");
  printf("Falsifiers and negative controls fixed in advance: result/M11_AUDIT.md\n\n");
  printf("Reference self-tests (the 12th artefact of the last session was a wrong\n");
  printf("reference that held a CONSTANT error at every domain size):\n");
  if (!ref_selftest(k)) {
    printf("\nREFERENCE SELF-TEST FAILED - nothing below would mean anything.\n");
    return 1;
  }

  if (mode != 1 && !experiment_a(k, 1)) return 1;
  if (mode != 2 && !experiment_a(k, 0)) return 1;
  if (mode == 0) {
    if (!experiment_b(k, 2000.0, 2)) return 1;
    if (!experiment_b(k, 2000.0, 5)) return 1;
    if (!experiment_b(k, 2000.0, 20)) return 1;
  }

  printf("\nREAD: warm = LSQR iterations on the CORRECTION after the camera step,\n");
  printf("cold = the same target from zero. -1 means the cap was hit. The claim\n");
  printf("under test is a CLASS difference (units vs hundreds), not a factor of 2.\n");
  return 0;
}
