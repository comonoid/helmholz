/* THE VERTICAL SLICE, step 3: the whole 2D pipeline in one place —
 * carrier basis + oblique cuts on a faceted cylinder + exact DtN on a circle,
 * solved and measured against the exact Mie series. Plan: SLICE_PLAN.md.
 *
 * TOTAL FIELD, NOT SCATTERED, AND THE MEDIUM TERM THEN VANISHES. With the
 * scattered field the interior would need BOTH wavenumbers (u_s = u - u_inc
 * carries k1 from u and k0 from u_inc) and a volume source integral. Solving
 * for the TOTAL field instead:
 *   - every element carries the wavenumber of the medium it sits in, and every
 *     element touching the boundary is CUT, so on each entry's region the
 *     carrier matches the medium and the term (k^2 - k_j^2) is IDENTICALLY
 *     ZERO. The medium enters only through the geometry of the cut and the
 *     choice of carrier — which is what "the carrier is local to the medium"
 *     means operationally;
 *   - the operator rows are homogeneous, so there is no volume integral of a
 *     source anywhere;
 *   - the incident wave enters ONLY through the DtN rows, whose right-hand side
 *     is the incoming half of the plane wave — the very quantity gate D2-1
 *     verified to 7.8e-15.
 * So the whole physics of the scene sits in three places, each already gated.
 *
 * Rows are equilibrated to unit norm ONCE (M11 R5: a scale taken from the
 * assembled matrix changes with the basis, and then two configurations solve
 * slightly different least-squares problems and the difference is blamed on the
 * physics). This is a REFORMULATION, not preconditioning (M12 R4). */
#include "bessel.h"
#include "carrier2d.h"
#include "cut2d.h"
#include "dtn2d.h"
#include "mie2d.h"
#include "nitsche2d.h"
#include "phi.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { MAXDIM = 30000, MAXM = 120 };
static const double LAM = 1.0;

/* L2 inverse inequality for a quadratic on a piece of length W: the extremal
 * polynomial on [-1,1] is u^2 - 1/3 and gives ||p'||/||p|| = sqrt(15), i.e.
 * 2 sqrt(15)/W after scaling. This is what bounds the Nitsche penalty FROM
 * BELOW; it is computed, not swept (M15 audit §2). */
static const double CINV = 7.745966692414834;

typedef struct {
  int nx, ny;
  int side; /* -1 uncut, 0 inside the scatterer, 1 outside */
  double complex kx, ky;
  double kc2; /* |k|^2 of the CARRIER: equal to km2 while the carrier is local
               * to the medium, exactly 0 with the carrier off. Stored rather than
               * recomputed as kx^2+ky^2 so that the Trefftz cancellation is exact
               * in floating point and not merely small. */
  double km2; /* k^2 of the MEDIUM on this element's region — not of its carrier.
               * The two are equal while the carrier is local to the medium (the
               * Trefftz cancellation), and they part company the moment the
               * carrier is switched off, which is exactly the baseline run. */
} elem;

/* Nitsche variants. FULL is the method; the other three exist because the audit
 * demands negative controls whose failure is predicted, not because anything is
 * meant to be chosen among them. */
enum { NIT_FULL = 0, NIT_PENALTY_ONLY = 1, NIT_FLIPPED = 2, NIT_OFF = 3 };

typedef struct {
  double a, k0, n;  /* scatterer radius, background k, index */
  double R;         /* DtN circle */
  double W;         /* element width */
  int nd;           /* directions per element */
  int nfacet;       /* facets of the cylinder */
  int mmax, ntheta; /* DtN truncation and sampling */
  double dtnw;      /* weight of the DtN block, see main() */
  double pen;       /* gamma: penalty = pen * (k_max + CINV/W), pen >= 1 by the bound */
  double beta;      /* the penalty itself, formed once in main() */
  int nitmode;      /* NIT_* above */
  int forcecut;     /* cut the basis even at zero contrast (falsifier F-N1) */
  int carrier;      /* 0 = cost baseline: plain partition-of-unity basis, no plane wave */
  int nv;           /* vertices of the cut polygon = facets */
  double vx[HZ_CUT_MAXHP], vy[HZ_CUT_MAXHP];
} cfg;

/* --- geometry ------------------------------------------------------------ */
static int straddles(const cfg *c, int nx, int ny) {
  /* No contrast, nothing to cut. In vacuum the cut is geometry without physics,
   * and it doubles the degrees of freedom on a whole ring of elements for no
   * reason — which is exactly the kind of free redundancy that fed the non-zero
   * field kernel (PLAN A5). Keeping it there would confound the direction fan
   * with the cut.
   * forcecut turns it back on WITHOUT the contrast, and that is the sharpest
   * control this bench has: the exact solution is then the incident plane wave,
   * it lies in the span exactly and with EQUAL coefficients on the two sides, so
   * every jump is zero and the Nitsche term must be invisible. Measured before
   * Nitsche: 0.197..0.405 where the uncut basis gives 0.0073..0.0166. */
  if (!c->forcecut && fabs(c->n - 1.0) < 1e-12) return 0;
  double px = c->W * (double)nx, py = c->W * (double)ny;
  double d = fabs(sqrt(px * px + py * py) - c->a);
  return d < 2.0 * c->W * sqrt(2.0); /* support half-diagonal */
}

static int build(const cfg *c, elem *b) {
  int nmax = (int)((c->R + 2.0 * c->W) / c->W) + 2, dim = 0;
  for (int nx = -nmax; nx <= nmax; nx++)
    for (int ny = -nmax; ny <= nmax; ny++) {
      double px = c->W * (double)nx, py = c->W * (double)ny;
      double r = sqrt(px * px + py * py);
      if (r > c->R + 4.0 * c->W) continue; /* columns: rows need n+-4 neighbours */
      int str = straddles(c, nx, ny);
      for (int pass = 0; pass < 2; pass++) {
        int side, inside;
        if (str) {
          side = pass; /* 0 = inside, 1 = outside */
          inside = (pass == 0);
        } else {
          if (pass == 1) break;
          side = -1;
          inside = (r < c->a);
        }
        double km = inside ? c->k0 * c->n : c->k0;
        for (int d = 0; d < c->nd; d++) {
          /* carrier = 0 is the COST BASELINE: the same bench, the same boundary,
           * the same reference and the same error metric, with the one thing
           * under test removed. What is left is an ordinary quadratic
           * partition-of-unity discretisation, so the comparison isolates the
           * basis and nothing else. */
          /* NO half-step offset: the drive is a plane wave along +x, and offsetting
           * the fan puts it MAXIMALLY far from every candidate (22.5 deg at
           * ND=8, i.e. 1.2 rad of phase across an element at W=lambda/2). The
           * incident direction must BE in the set. */
          double th = 2.0 * M_PI * (double)d / (double)c->nd;
          if (dim >= MAXDIM) return dim;
          b[dim].nx = nx;
          b[dim].ny = ny;
          b[dim].side = side;
          b[dim].kx = c->carrier ? km * cos(th) : 0.0;
          b[dim].ky = c->carrier ? km * sin(th) : 0.0;
          b[dim].km2 = km * km;
          b[dim].kc2 = c->carrier ? km * km : 0.0;
          dim++;
        }
      }
    }
  return dim;
}

/* One term of the Galerkin entry over the region of the pair. side handling:
 * inside is the convex polygon, outside is "whole minus inside". */
static double complex region_int(const hz_half2 *hp, int nhp, int sidei, int sidej, double xlo,
                                 double xhi, double ylo, double yhi, hz_axis2 fx, hz_axis2 fy,
                                 double complex omx, double complex omy) {
  int in = (sidei == 0) || (sidej == 0);
  int out = (sidei == 1) || (sidej == 1);
  if (in && out) return 0.0; /* opposite sides of the cut: disjoint */
  if (in) return hz_cut2d_poly(xlo, xhi, ylo, yhi, fx, fy, omx, omy, hp, nhp);
  if (out)
    return hz_cut2d_poly(xlo, xhi, ylo, yhi, fx, fy, omx, omy, hp, 0) -
           hz_cut2d_poly(xlo, xhi, ylo, yhi, fx, fy, omx, omy, hp, nhp);
  return hz_cut2d_poly(xlo, xhi, ylo, yhi, fx, fy, omx, omy, hp, 0);
}

/* k^2 of the MEDIUM on the region a pair shares.
 * If either function is cut, its side IS the region. If both are uncut then
 * neither support reaches the interface (straddles() is conservative), so an
 * overlapping pair lies wholly on one side and the two agree. */
static double pair_k2(elem ei, elem ej) {
  return ei.side >= 0 ? ei.km2 : ej.km2;
}

/* Fast path for a pair that neither cut touches: the region is the plain
 * support rectangle, so the integral factorises and costs six 1D closed forms
 * instead of four polygon traversals. Verified against the polygon path to
 * 2e-13, and against hz_carrier2d_entry (gate D1) to the same. */
static double complex entry_sep(double W, elem ei, elem ej) {
  double xlo = W * ((double)(ei.nx > ej.nx ? ei.nx : ej.nx) - 2.0);
  double xhi = W * ((double)(ei.nx < ej.nx ? ei.nx : ej.nx) + 2.0);
  double ylo = W * ((double)(ei.ny > ej.ny ? ei.ny : ej.ny) - 2.0);
  double yhi = W * ((double)(ei.ny < ej.ny ? ei.ny : ej.ny) + 2.0);
  if (!(xlo < xhi) || !(ylo < yhi)) return 0.0;
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
  /* THE MEDIUM TERM, WHICH IS IDENTICALLY ZERO WHENEVER THE CARRIER IS LOCAL TO
   * THE MEDIUM. entry() computes Int B_i (Lap + |k_j|^2) B_j, because the -k^2
   * of the carrier cancels against the +k^2 of the operator; what the equation
   * actually needs is the medium's k^2, so the difference has to be carried
   * explicitly. With the carrier on it is exactly 0 (|k_j|^2 == km2 by
   * construction) and nothing changes; with the carrier off it is the whole
   * Helmholtz term. */
  double dk2 = pair_k2(ei, ej) - ej.kc2;
  return pref *
         ((X2 + 2.0 * i1 * ej.kx * X1) * Y0 + X0 * (Y2 + 2.0 * i1 * ej.ky * Y1) + dk2 * X0 * Y0);
}

/* Trace coefficients of one basis function on the cut locus:
 *   [v] = alpha * v|_G ,   {v} = beta * v|_G .
 * A cut function is identically zero on the far side, so its jump is its own
 * trace up to the sign of the side and its average is half of it; an uncut
 * function is continuous, so it has no jump and its average is the trace. */
static double trace_alpha(int side) {
  if (side < 0) return 0.0;        /* uncut: continuous, no jump */
  return (side == 1) ? 1.0 : -1.0; /* [.] = outside - inside */
}

static double trace_beta(int side) {
  return side < 0 ? 1.0 : 0.5;
}

/* The interface part of the entry. Derivation in src/nitsche2d.h; with
 *   Strong_ij = -a_broken - Int_G ([B_i]{dnB_j} + {B_i}[dnB_j])
 * and a_broken symmetric, the assembled matrix becomes MINUS the symmetric
 * Nitsche form once these terms are added:
 *   A_ij = Strong_ij + beta_i alpha_j (t1j - t1i) - beta_pen alpha_i alpha_j t0
 * The overall sign is free — operator rows have an identically zero right-hand
 * side — but the signs WITHIN the row are not.
 * THE SIGN OF THE CONSISTENCY PAIR IS THE ONE THING THE SYMMETRY CHECK CANNOT
 * SEE, because flipping a symmetric pair leaves the matrix symmetric. It is
 * fixed instead by consistency: for the exact solution [u] = [dn u] = 0, yet
 * a_broken(u,v) = -Int_G [v]{dn u} is NOT zero (the TEST function still jumps),
 * so the interface term must CANCEL it rather than double it. Written with the
 * other sign the assembled operator failed the known plane wave by |Ac| = 0.40
 * (F-N4) while still passing the symmetry check at 8e-16 — which is exactly why
 * both falsifiers were registered. */
static double complex nitsche_term(const cfg *c, elem ei, elem ej) {
  if (c->nitmode == NIT_OFF || (ei.side < 0 && ej.side < 0)) return 0.0;
  double ai = trace_alpha(ei.side), aj = trace_alpha(ej.side), bi = trace_beta(ei.side);
  if (c->nitmode == NIT_FLIPPED) { /* NC-N1: break ONE side of the comparison */
    ai = -ai;
    aj = -aj;
  }
  hz_carrier2d ci = {c->W, ei.nx, ei.ny, ei.kx, ei.ky};
  hz_carrier2d cj = {c->W, ej.nx, ej.ny, ej.kx, ej.ky};
  hz_nit2d t = hz_nitsche2d_poly(ci, cj, c->vx, c->vy, c->nv);
  double complex v = -c->beta * (ai * aj) * t.t0;
  if (c->nitmode != NIT_PENALTY_ONLY) /* NC-N2 drops exactly the consistency pair */
    v += (bi * aj) * (t.t1j - t.t1i);
  return v;
}

static double complex entry(const cfg *c, const hz_half2 *hp, int nhp, elem ei, elem ej) {
  double W = c->W;
  if (ei.side < 0 && ej.side < 0) return entry_sep(W, ei, ej);
  double xlo = W * ((double)(ei.nx > ej.nx ? ei.nx : ej.nx) - 2.0);
  double xhi = W * ((double)(ei.nx < ej.nx ? ei.nx : ej.nx) + 2.0);
  double ylo = W * ((double)(ei.ny > ej.ny ? ei.ny : ej.ny) - 2.0);
  double yhi = W * ((double)(ei.ny < ej.ny ? ei.ny : ej.ny) + 2.0);
  if (!(xlo < xhi) || !(ylo < yhi)) return 0.0;
  double complex omx = ei.kx + ej.kx, omy = ei.ky + ej.ky; /* bilinear: add */
  double complex i1 = CMPLX(0.0, 1.0);
  double complex pref = cexp(-i1 * (ei.kx * W * (double)ei.nx + ej.kx * W * (double)ej.nx +
                                    ei.ky * W * (double)ei.ny + ej.ky * W * (double)ej.ny));
  hz_phi_factor fi0 = {W, (double)ei.nx, 0}, fj0 = {W, (double)ej.nx, 0};
  hz_phi_factor gi0 = {W, (double)ei.ny, 0}, gj0 = {W, (double)ej.ny, 0};
  hz_phi_factor fj1 = {W, (double)ej.nx, 1}, fj2 = {W, (double)ej.nx, 2};
  hz_phi_factor gj1 = {W, (double)ej.ny, 1}, gj2 = {W, (double)ej.ny, 2};
  hz_axis2 X0 = {{fi0, fj0}, 2}, X1 = {{fi0, fj1}, 2}, X2 = {{fi0, fj2}, 2};
  hz_axis2 Y0 = {{gi0, gj0}, 2}, Y1 = {{gi0, gj1}, 2}, Y2 = {{gi0, gj2}, 2};
  double complex t = 0.0;
  t += region_int(hp, nhp, ei.side, ej.side, xlo, xhi, ylo, yhi, X2, Y0, omx, omy);
  t += 2.0 * i1 * ej.kx *
       region_int(hp, nhp, ei.side, ej.side, xlo, xhi, ylo, yhi, X1, Y0, omx, omy);
  t += region_int(hp, nhp, ei.side, ej.side, xlo, xhi, ylo, yhi, X0, Y2, omx, omy);
  t += 2.0 * i1 * ej.ky *
       region_int(hp, nhp, ei.side, ej.side, xlo, xhi, ylo, yhi, X0, Y1, omx, omy);
  /* Medium term: identically zero while the carrier is local to the medium (see
   * entry_sep), the whole Helmholtz term once the carrier is switched off. */
  double dk2 = pair_k2(ei, ej) - ej.kc2;
  if (fabs(dk2) > 0.0)
    t += dk2 * region_int(hp, nhp, ei.side, ej.side, xlo, xhi, ylo, yhi, X0, Y0, omx, omy);
  return pref * t + nitsche_term(c, ei, ej);
}

/* Point against the FACETED boundary — the same curve the basis is cut by, not
 * the circle it approximates. Counter-clockwise ring, so "inside" is "left of
 * every edge". */
static int poly_inside(const cfg *c, double x, double y) {
  for (int e = 0; e < c->nv; e++) {
    int f = (e + 1) % c->nv;
    double dx = c->vx[f] - c->vx[e], dy = c->vy[f] - c->vy[e];
    if (dx * (y - c->vy[e]) - dy * (x - c->vx[e]) < 0.0) return 0;
  }
  return 1;
}

static double complex field_at(const cfg *c, const elem *b, int dim, const double complex *coef,
                               double x, double y) {
  int in = poly_inside(c, x, y);
  double complex got = 0.0;
  for (int j = 0; j < dim; j++) {
    if (b[j].side == 0 && !in) continue; /* a cut function lives on one side only */
    if (b[j].side == 1 && in) continue;
    hz_carrier2d cb = {c->W, b[j].nx, b[j].ny, b[j].kx, b[j].ky};
    got += coef[j] * hz_carrier2d_val(cb, x, y);
  }
  return got;
}

/* THE OTHER HALF OF GATE D2: the SCATTERING DIAGRAM.
 * The field error is normalised by the TOTAL field, which is everywhere of order
 * one because of the incident wave; the diagram is normalised by the SCATTERED
 * field alone, so it removes that dilution — the same correction M14 had to make
 * with its +-2W window rule. Nothing new is approximated here: the solution is
 * sampled on the measurement circle, the incident wave is subtracted, and the
 * result is projected onto outgoing harmonics by a DFT and one Hankel division,
 *     u_s = sum_m A_m H_|m|(k0 r) e^{i m th},   A_m = i^|m| c_|m|
 * so c_m comes straight out. The exact side is hz_mie_far, already gated (8/8).
 *
 * m AND -m ARE KEPT APART AND NEVER SYMMETRISED. The scene is symmetric about
 * the x axis (axial incidence; the facet ring has a vertex at theta = 0), so the
 * EXACT c_m obey c_{-m} = c_{+m} while the numerical ones need not: the basis
 * grid is square and its direction fan is its own. The gap between them is
 * therefore a control that does not know the reference exists — it measures the
 * bench. Symmetrising would hide exactly that half of the error. */
static double diagram_err(const cfg *c, const elem *b, int dim, const double complex *coef,
                          const hz_mie *mie, const double complex *mc, double rm, double *asym) {
  enum { NT = 512 };
  static double complex us[NT];
  for (int p = 0; p < NT; p++) {
    double th = 2.0 * M_PI * (double)p / (double)NT;
    double x = rm * cos(th), y = rm * sin(th);
    us[p] = field_at(c, b, dim, coef, x, y) - cexp(CMPLX(0.0, 1.0) * c->k0 * x);
  }
  double complex cp[MAXM + 1], cm[MAXM + 1];
  double cmax = 0.0, dmax = 0.0;
  for (int m = 0; m <= c->mmax; m++) {
    double complex up = 0.0, um = 0.0;
    for (int p = 0; p < NT; p++) {
      double th = 2.0 * M_PI * (double)p / (double)NT;
      up += us[p] * cexp(CMPLX(0.0, -1.0) * (double)m * th);
      um += us[p] * cexp(CMPLX(0.0, 1.0) * (double)m * th);
    }
    up /= (double)NT;
    um /= (double)NT;
    double jj, yy, dj, dy;
    hz_bessel_jy(m, c->k0 * rm, &jj, &yy, &dj, &dy);
    double complex den = CMPLX(jj, yy) * cexp(CMPLX(0.0, 1.0) * M_PI * 0.5 * (double)m);
    cp[m] = up / den;
    cm[m] = um / den;
    if (cabs(cp[m]) > cmax) cmax = cabs(cp[m]);
    if (cabs(cp[m] - cm[m]) > dmax) dmax = cabs(cp[m] - cm[m]);
  }
  *asym = cmax > 0.0 ? dmax / cmax : 0.0;
  double num = 0.0, den = 0.0;
  for (int p = 0; p < NT; p++) {
    double th = 2.0 * M_PI * (double)p / (double)NT;
    double complex fn = cp[0];
    for (int m = 1; m <= c->mmax; m++)
      fn += cp[m] * cexp(CMPLX(0.0, 1.0) * (double)m * th) +
            cm[m] * cexp(CMPLX(0.0, -1.0) * (double)m * th);
    double complex fe = hz_mie_far(mie, mc, th);
    num += cabs(fn - fe) * cabs(fn - fe);
    den += cabs(fe) * cabs(fe);
  }
  return den > 0.0 ? sqrt(num / den) : 0.0;
}

/* Field error on one circle — the quantity the gate is stated in. */
static double field_err(const cfg *c, const elem *b, int dim, const double complex *coef,
                        const hz_mie *mie, const double complex *mc, double rr) {
  double num = 0.0, den = 0.0;
  for (int p = 0; p < 256; p++) {
    double th = 2.0 * M_PI * (double)p / 256.0;
    double x = rr * cos(th), y = rr * sin(th);
    double complex got = field_at(c, b, dim, coef, x, y);
    double complex ex = cexp(CMPLX(0.0, 1.0) * c->k0 * x) + hz_mie_scat(mie, mc, x, y);
    num += cabs(got - ex) * cabs(got - ex);
    den += cabs(ex) * cabs(ex);
  }
  return sqrt(num / den);
}

/* KEY=VALUE arguments, and an unknown key is an ABORT rather than a default.
 * Thirteen positional arguments were ready ground for the seventeenth artefact
 * of this project: one misplaced position in a sweep is indistinguishable from
 * physics in the output. */
static const char *const KEYS[] = {"W",   "nd",  "facets", "R",        "mmax",    "dtnw",    "it",
                                   "n",   "pen", "nit",    "forcecut", "rrow",    "ngam",    "gamw",
                                   "sym", "a5",  "rm",     "plateau",  "carrier", "ftarget", NULL};

static double argval(int argc, char **argv, const char *key, double def) {
  size_t kl = strlen(key);
  for (int i = 1; i < argc; i++)
    if (!strncmp(argv[i], key, kl) && argv[i][kl] == '=') return atof(argv[i] + kl + 1);
  return def;
}

static int args_ok(int argc, char **argv) {
  for (int i = 1; i < argc; i++) {
    const char *eq = strchr(argv[i], '=');
    int ok = 0;
    if (eq)
      for (int k = 0; KEYS[k]; k++)
        if ((size_t)(eq - argv[i]) == strlen(KEYS[k]) &&
            !strncmp(argv[i], KEYS[k], strlen(KEYS[k])))
          ok = 1;
    if (!ok) {
      printf("slice2d: unknown argument '%s'; keys are", argv[i]);
      for (int k = 0; KEYS[k]; k++)
        printf(" %s", KEYS[k]);
      printf("\n");
      return 0;
    }
  }
  return 1;
}

int main(int argc, char **argv) {
  cfg c = {0};
  c.a = 1.0 * LAM;
  c.k0 = 2.0 * M_PI / LAM;
  c.ntheta = 512;
  if (!args_ok(argc, argv)) return 1;
  c.W = argval(argc, argv, "W", 1.0) * LAM;
  c.nd = (int)argval(argc, argv, "nd", 8);
  c.nfacet = (int)argval(argc, argv, "facets", 32);
  c.R = argval(argc, argv, "R", 2.5) * LAM;
  c.mmax = (int)argval(argc, argv, "mmax", 25);
  c.dtnw = argval(argc, argv, "dtnw", 1.0);
  c.n = argval(argc, argv, "n", 1.5);
  c.pen = argval(argc, argv, "pen", 1.0);
  c.nitmode = (int)argval(argc, argv, "nit", NIT_FULL);
  c.forcecut = (int)argval(argc, argv, "forcecut", 0);
  c.carrier = (int)argval(argc, argv, "carrier", 1);
  if (!c.carrier && c.nd != 1) {
    /* every direction collapses to the same column without a carrier */
    printf("  carrier=0: ND forced from %d to 1 (directions are identical without a carrier)\n",
           c.nd);
    c.nd = 1;
  }
  /* THE PENALTY IS FORMED FROM THE BOUND, NOT FROM A SWEEP (M15 audit §2):
   * beta >= k_max + CINV/W, and pen is the dimensionless factor on top of it
   * whose only legitimate value is "anywhere on the plateau". */
  double kmax = c.k0 * (c.n > 1.0 ? c.n : 1.0);
  c.beta = c.pen * (kmax + CINV / c.W);
  c.nv = hz_cut2d_disc_verts(0.0, 0.0, c.a, c.nfacet, 1, c.vx, c.vy);
  printf("slice2d: ka=%.2f n=%.2f R=%.2flam W=%.3flam ND=%d facets=%d mmax=%d\n", c.k0 * c.a, c.n,
         c.R / LAM, c.W / LAM, c.nd, c.nfacet, c.mmax);
  printf("  nitsche: mode=%d pen=%.4g beta=%.4f (kmax=%.3f + %.3f/W) forcecut=%d dtnw=%.4g\n",
         c.nitmode, c.pen, c.beta, kmax, CINV, c.forcecut, c.dtnw);

  static elem b[MAXDIM];
  int dim = build(&c, b);
  int ngam0 = (int)argval(argc, argv, "ngam", 0);
  int nrow = dim + 2 * c.mmax + 1 + 2 * ngam0;
  printf("  dim=%d  operator rows=%d  DtN rows=%d  total rows=%d\n", dim, dim, 2 * c.mmax + 1,
         nrow);
  if (dim >= MAXDIM) {
    printf("  ABORT: basis truncated\n");
    return 1;
  }

  hz_half2 hp[HZ_CUT_MAXHP];
  int nhp = hz_cut2d_disc(0.0, 0.0, c.a, c.nfacet, 1, hp);
  hz_dtn dt = {c.R, c.k0, c.mmax, c.ntheta};

  /* SPARSE (CSR). The dense matrix is not an option any more: the geometry that
   * makes the system determined needs R >~ 20 W and W >~ lambda/2, hence
   * dim >~ 12000, and dense would be terabytes. Operator rows touch only
   * |dn| <= 3 per axis (the support is 4W), DtN rows only elements that reach
   * the circle. */
  size_t cap = 8u << 20;
  int *ja = malloc(cap * sizeof(int));
  double complex *va = malloc(cap * sizeof(double complex));
  size_t *rp = calloc((size_t)nrow + 1, sizeof(size_t));
  double complex *rhs = calloc((size_t)nrow, sizeof(double complex));
  double *rowsc = calloc((size_t)nrow, sizeof(double));
  if (!ja || !va || !rp || !rhs || !rowsc) return 1;
  size_t nnz = 0;

  /* OPERATOR ROWS ONLY WHERE THE NEIGHBOURHOOD IS COMPLETE.
   * The support of phi is 4W, so a test function needs trial elements out to
   * n +- 4 before a plane wave can be reproduced across its support (artefact
   * 15). At the outer edge of the column set those neighbours do not exist, and
   * such a row can only be satisfied by a field that VANISHES there — a
   * homogeneous Dirichlet condition smuggled in by the truncation, which then
   * propagates inward and kills the solution. Measured before this guard:
   * norm(c) = 2.0 while the field it produced was 1e-11, i.e. the solver
   * returned a null mode and the boundary drive did nothing at all.
   * Rows are therefore written only for elements whose whole support sits
   * inside the column set; the outer ring keeps its COLUMNS (it must, or the
   * interior loses its neighbours) and simply carries no equation. */
  int nop = 0;
  /* ROWS FOR EVERY COLUMN. The earlier cutoff at R-2W was wrong: <B_i, Lu> = 0
   * is an INTEGRAL IDENTITY satisfied by the true solution for ANY test
   * function, complete neighbourhood or not, so dropping those rows removed
   * legitimate equations and left the system underdetermined (PLAN A5).
   * rrow= restores the cutoff for comparison. */
  double rrow = argval(argc, argv, "rrow", 0.0) * c.W;
  if (!(rrow > 0.0)) rrow = c.R;
  for (int i = 0; i < dim; i++) {
    rp[i] = nnz;
    double px = c.W * (double)b[i].nx, py = c.W * (double)b[i].ny;
    if (sqrt(px * px + py * py) > rrow) continue;
    nop++;
    for (int j = 0; j < dim; j++) {
      if (abs(b[i].nx - b[j].nx) > 3 || abs(b[i].ny - b[j].ny) > 3) continue;
      double complex v = entry(&c, hp, nhp, b[i], b[j]);
      if (!(cabs(v) > 0.0) || nnz >= cap) continue;
      ja[nnz] = j;
      va[nnz++] = v;
    }
  }
  printf("  operator rows with a complete neighbourhood: %d of %d, nnz=%zu\n", nop, dim, nnz);
  double reach = 2.0 * c.W * sqrt(2.0);
  for (int mi = 0; mi <= 2 * c.mmax; mi++) {
    int m = mi - c.mmax;
    size_t row = (size_t)(dim + mi);
    rp[row] = nnz;
    for (int j = 0; j < dim; j++) {
      double px = c.W * (double)b[j].nx, py = c.W * (double)b[j].ny;
      if (fabs(sqrt(px * px + py * py) - c.R) > reach) continue;
      hz_carrier2d cb = {c.W, b[j].nx, b[j].ny, b[j].kx, b[j].ky};
      double complex v = hz_dtn_row(&dt, m, cb);
      if (!(cabs(v) > 0.0) || nnz >= cap) continue;
      ja[nnz] = j;
      va[nnz++] = v;
    }
    /* RHS: the incident wave is not outgoing, so (d/dr - Z) applied to it does
     * not vanish; that residual IS the drive. Closed form, gate D2-1. */
    int am = m < 0 ? -m : m;
    double j0v, y0v, dj0, dy0;
    hz_bessel_jy(am, c.k0 * c.R, &j0v, &y0v, &dj0, &dy0);
    double complex ipow = cexp(CMPLX(0.0, 1.0) * M_PI * 0.5 * (double)am);
    rhs[row] = ipow * (c.k0 * dj0 - hz_dtn_symbol(&dt, m) * j0v);
  }
  /* --- CONTINUITY ROWS ON THE CUT (what the strong form cannot see) -------
   * A cut basis function has a JUMP, and the strong form integrates each side
   * separately, so the distributional terms on the cut locus (a delta and a
   * delta-prime) are silently dropped: the discrete operator does not see the
   * jump AT ALL. Measured — with the cut on, the vacuum error was 0.197..0.405;
   * with it off, 0.0073..0.0166. The jump was free, so every combination
   * a*B_in + b*B_out with a != b sat in the kernel WITH a non-zero field.
   * The physical solution of a piecewise-constant Helmholtz problem has BOTH u
   * and du/dn continuous across the material interface — the cut exists so each
   * side can carry its own wavenumber, not so the field can jump. So the
   * missing information is added as explicit constraints at sample points on
   * the interface (M9b §8: "Nitsche / constraint / Heaviside"); this is the
   * constraint form, the cheapest of the three.
   * SUPERSEDED BY THE NITSCHE TERM and kept only as NC-N3: the sweep over the
   * number of points was a PIT (0.697/0.624/0.503/0.181/0.095/0.131/0.129 at
   * 0/8/16/24/32/48/96), i.e. a knob, and a number found at the bottom of a knob
   * sweep is not a result. It stays here so that "pit versus plateau" can be
   * measured by ONE binary on ONE configuration instead of against numbers from
   * a previous session. Off unless ngam= is given. */
  int ngam = ngam0;
  double gamw = argval(argc, argv, "gamw", 1.0);
  for (int g = 0; g < 2 * ngam; g++) {
    int ip = g / 2, deriv = g % 2;
    double th = 2.0 * M_PI * (double)ip / (double)ngam;
    double gx = c.a * cos(th), gy = c.a * sin(th);
    size_t row = (size_t)(dim + 2 * c.mmax + 1 + g);
    rp[row] = nnz;
    for (int j = 0; j < dim; j++) {
      if (b[j].side < 0) continue; /* uncut functions are continuous already */
      double tx = gx / c.W - (double)b[j].nx, ty = gy / c.W - (double)b[j].ny;
      if (fabs(tx) >= 2.0 || fabs(ty) >= 2.0) continue;
      double complex i1 = CMPLX(0.0, 1.0);
      double complex e = cexp(
          i1 * (b[j].kx * (gx - c.W * (double)b[j].nx) + b[j].ky * (gy - c.W * (double)b[j].ny)));
      double px = hz_phi(tx), py = hz_phi(ty);
      double complex v;
      if (deriv == 0) {
        v = px * py * e;
      } else {
        double complex bx = (hz_phi_d1(tx) / c.W + i1 * b[j].kx * px) * py * e;
        double complex by = px * (hz_phi_d1(ty) / c.W + i1 * b[j].ky * py) * e;
        v = (gx * bx + gy * by) / c.a;
      }
      v *= (b[j].side == 1) ? 1.0 : -1.0; /* [.] = outside - inside */
      if (!(cabs(v) > 0.0) || nnz >= cap) continue;
      ja[nnz] = j;
      va[nnz++] = v;
    }
    rhs[row] = 0.0;
  }
  rp[nrow] = nnz;
  if (nnz >= cap) {
    printf("  ABORT: nnz cap\n");
    return 1;
  }

  /* --- F-N3: SYMMETRY OF THE OPERATOR BLOCK, before equilibration ---------
   * a_broken is symmetric and so is the Nitsche form built on it, but the STRONG
   * form is not: integrating by parts over each side of a cut leaves an
   * uncancelled surface term. So the assembled block must be ASYMMETRIC without
   * the interface terms and symmetric with them, and that is an independent
   * witness of the derivation and of every sign in it — free, no solve.
   * (Equilibration destroys symmetry, hence "before".) */
  if (argval(argc, argv, "sym", 0.0) > 0.0) {
    double amax = 0.0, dmax = 0.0;
    for (int i = 0; i < dim; i++)
      for (size_t p = rp[i]; p < rp[i + 1]; p++) {
        int j = ja[p];
        if (rp[j + 1] <= rp[j]) continue; /* row j not assembled: no pair */
        size_t lo = rp[j], hi = rp[j + 1];
        double complex aji = 0.0;
        while (lo < hi) { /* columns within a row are ascending */
          size_t mid = lo + (hi - lo) / 2;
          if (ja[mid] < i)
            lo = mid + 1;
          else if (ja[mid] > i)
            hi = mid;
          else {
            aji = va[mid];
            break;
          }
        }
        if (cabs(va[p]) > amax) amax = cabs(va[p]);
        if (cabs(va[p] - aji) > dmax) dmax = cabs(va[p] - aji);
      }
    printf("  [F-N3] operator block asymmetry max|Aij-Aji|/max|Aij| = %.3e  (max|Aij|=%.3e)\n",
           amax > 0.0 ? dmax / amax : 0.0, amax);
  }

  /* Row equilibration, once — and then the DtN block is WEIGHTED UP.
   * WHY, AND IT IS NOT A TUNING KNOB. The radiation condition is a CONSTRAINT,
   * not one more residual to trade against: with 1744 operator rows and 41
   * boundary rows at equal weight, the minimiser buys operator residual by
   * throwing the drive away and returns u = 0 — measured, error exactly 1.0000.
   * This is the v1 failure mode verbatim (PLAN: "pure interior least squares
   * CUTS OUT RADIATION; u=0 gives zero residual in vacuum, the minimiser
   * prefers silence"), reappearing because the boundary rows are outnumbered.
   * The weight must therefore be swept and the answer must PLATEAU in it; a
   * result that moves with the weight is a tuned number, not a solution. */
  double rnmin = 1e300, rnmax = 0.0;
  for (size_t r = 0; r < (size_t)nrow; r++) {
    double s = 0.0;
    for (size_t p = rp[r]; p < rp[r + 1]; p++)
      s += cabs(va[p]) * cabs(va[p]);
    s = sqrt(s);
    if (!(s > 0.0)) continue;
    if (r < (size_t)dim) {
      if (s < rnmin) rnmin = s;
      if (s > rnmax) rnmax = s;
    }
    double w = 1.0;
    if (r >= (size_t)(dim + 2 * c.mmax + 1))
      w = gamw;
    else if (r >= (size_t)dim)
      w = c.dtnw;
    for (size_t p = rp[r]; p < rp[r + 1]; p++)
      va[p] *= w / s;
    rhs[r] *= w / s;
    rowsc[r] = w / s;
  }
  /* Sliver watch (M15 audit §6): an element whose region on one side is a thin
   * shard has a nearly zero row, and equilibration then promotes its noise to
   * the weight of a real equation. Alarm below 1e-10. */
  printf("  operator row norms before equilibration: min=%.3e max=%.3e ratio=%.2e\n", rnmin, rnmax,
         rnmax > 0.0 ? rnmin / rnmax : 0.0);

  /* --- DIAGNOSTIC: feed the KNOWN vector into the assembled rows ----------
   * In vacuum the exact solution is the incident plane wave, and its expansion
   * coefficients are closed form (partition of unity, /16) on the d = 0
   * direction. Gate D1 says the operator block must annihilate it; gate D2-1
   * says the DtN block must reproduce b exactly. Whichever block does not is
   * the one with the defect — one run instead of a knob hunt. */
  double complex *ctrue = calloc((size_t)dim, sizeof(double complex));
  if (!ctrue) return 1;
  /* Both sides of a cut element get the SAME coefficient: B_in + B_out is the
   * uncut function, so the plane wave is reproduced exactly and every jump is
   * identically zero. That is what makes this vector a test of the Nitsche term
   * rather than of the basis. */
  for (int j = 0; j < dim; j++) {
    if (cabs(b[j].ky) > 1e-9 || creal(b[j].kx) < 0.0) continue; /* d = 0 only */
    ctrue[j] = cexp(CMPLX(0.0, 1.0) * (creal(b[j].kx) * c.W * (double)b[j].nx)) / 16.0;
  }
  {
    double eop = 0.0, edt = 0.0, nbd = 0.0;
    for (size_t r = 0; r < (size_t)nrow; r++) {
      double complex s = 0.0;
      for (size_t p = rp[r]; p < rp[r + 1]; p++)
        s += va[p] * ctrue[ja[p]];
      double d2 = cabs(s - rhs[r]) * cabs(s - rhs[r]);
      if (r < (size_t)dim)
        eop += d2;
      else {
        edt += d2;
        nbd += cabs(rhs[r]) * cabs(rhs[r]);
      }
    }
    /* Only meaningful in vacuum: with a scatterer present the plane wave is not
     * the solution, so |Ac| there is a number about nothing. Printing it anyway
     * would be an invitation to read it as a defect. */
    if (fabs(c.n - 1.0) < 1e-12)
      printf("  [F-N4] plane-wave vector: operator block |Ac|=%.3e   DtN block "
             "|Ac-b|=%.3e of |b|=%.3e\n",
             sqrt(eop), sqrt(edt), sqrt(nbd));
  }

  /* The reference is built BEFORE the solve, because the plateau curve needs it:
   * §8 of SLICE_PLAN asks for the discretisation floor as a CURVE of error
   * against iterations, not as one number, and a curve assembled from separate
   * runs cannot distinguish "converged" from "stopped". */
  hz_mie mie = {HZ_MIE_SHARP, c.a, c.k0, c.n, c.mmax, 8000};
  static double complex mc[MAXM + 1];
  hz_mie_coeffs(&mie, mc);
  if (cabs(mc[0]) > 0.0)
    printf("  Mie tail |c_M|/max = %.2e\n", cabs(mc[c.mmax]) / cabs(mc[0]));
  else
    printf("  Mie coefficients identically zero (no contrast): exact = incident wave\n");
  double rm = argval(argc, argv, "rm", 2.0) * LAM; /* gate window, NEVER swept */
  int plateau = (int)argval(argc, argv, "plateau", 0.0);

  /* --- LSQR (Paige-Saunders), min-norm least squares, matrix by CSR ------- */
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
  for (size_t r = 0; r < (size_t)nrow; r++)
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
  enum { NTOL = 3 };
  static const double restol[NTOL] = {1e-2, 1e-3, 1e-4};
  int nres[NTOL] = {0, 0, 0}, nfmin = 0, nftgt = 0;
  double ftarget = argval(argc, argv, "ftarget", 1e-2);
  double fmin = 1e300;
  int itmax = (int)argval(argc, argv, "it", 4000);
  int itdone = 0, stop = 0; /* 1 = beta breakdown, 2 = alpha breakdown */
  for (int it = 0; it < itmax; it++) {
    itdone = it + 1;
    for (size_t r = 0; r < (size_t)nrow; r++) {
      double complex s = 0.0;
      for (size_t p = rp[r]; p < rp[r + 1]; p++)
        s += va[p] * vv[ja[p]];
      uu[r] = s - alpha * uu[r];
    }
    beta = 0.0;
    for (int r = 0; r < nrow; r++)
      beta += cabs(uu[r]) * cabs(uu[r]);
    beta = sqrt(beta);
    if (!(beta > 0.0)) {
      stop = 1;
      break;
    }
    for (int r = 0; r < nrow; r++)
      uu[r] /= beta;
    for (int j = 0; j < dim; j++)
      tv[j] = 0.0;
    for (size_t r = 0; r < (size_t)nrow; r++)
      for (size_t p = rp[r]; p < rp[r + 1]; p++)
        tv[ja[p]] += conj(va[p]) * uu[r];
    alpha = 0.0;
    for (int j = 0; j < dim; j++) {
      vv[j] = tv[j] - beta * vv[j];
      alpha += cabs(vv[j]) * cabs(vv[j]);
    }
    alpha = sqrt(alpha);
    if (!(alpha > 0.0)) {
      stop = 2;
      break;
    }
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
    /* M17: iterations needed, by the two criteria that can disagree. The
     * residual one is pure solver; the field one is what the plan actually
     * demands ("convergence judged on the FIELD, not the raw residual" — in v1
     * the residual stagnated on a boundary mode that did not harm the field). */
    for (int t = 0; t < NTOL; t++)
      if (nres[t] == 0 && phibar / b0 <= restol[t]) nres[t] = it + 1;
    if (plateau > 0 && (it + 1) % plateau == 0) {
      double as = 0.0;
      double fe = field_err(&c, b, dim, xs, &mie, mc, rm);
      if (nftgt == 0 && fe <= ftarget) nftgt = it + 1;
      if (fe < fmin) {
        fmin = fe;
        nfmin = it + 1;
      }
      printf("    it %6d  |r|/|b| = %.4e  field(r=%.2f) = %.4e  diagram = %.4e\n", it + 1,
             phibar / b0, rm, fe, diagram_err(&c, b, dim, xs, &mie, mc, rm, &as));
    } else if (plateau == 0 && (it + 1) % 500 == 0) {
      printf("    it %5d  |r|/|b| = %.4e\n", it + 1, phibar / b0);
    }
  }
  {
    double nc = 0.0;
    for (int j = 0; j < dim; j++)
      nc += cabs(xs[j]) * cabs(xs[j]);
    printf("  LSQR done: |c| = %.4e  |r|/|b| = %.4e  iters = %d%s\n", sqrt(nc), phibar / b0, itdone,
           stop == 1   ? " (beta breakdown: Krylov space exhausted)"
           : stop == 2 ? " (alpha breakdown)"
                       : " (cap reached)");
  }

  /* COST, in the only unit that compares two different bases honestly: work is
   * iterations x nnz, because a row of the carrier basis is far denser than a
   * row of the plain one and bare iteration counts would flatter it. */
  printf("  [COST] dim=%d nnz=%zu  iters to |r|/|b|<1e-2/1e-3/1e-4: %d/%d/%d   field min %.4e at "
         "it %d   iters to field<%.1e: %d   work(it x nnz) = %.3e\n",
         dim, nnz, nres[0], nres[1], nres[2], fmin, nfmin, ftarget, nftgt,
         (double)nftgt * (double)nnz);

  /* --- error against the exact Mie series on a circle -------------------- */
  {
    double asym = 0.0;
    double dg = diagram_err(&c, b, dim, xs, &mie, mc, rm, &asym);
    printf("  [GATE D2] field(r=%.2flam) = %.4e   DIAGRAM = %.4e   (threshold 5e-2)\n", rm / LAM,
           field_err(&c, b, dim, xs, &mie, mc, rm), dg);
    /* Reference-free: the exact solution has c_{-m} = c_{+m}; the basis grid has
     * no reason to. If this is not well below the diagram error, part of what is
     * being called method error is the bench's own asymmetry. */
    printf("  [control, no reference] harmonic asymmetry |c_+m - c_-m|/max|c| = %.3e (%.0f%% of "
           "the diagram error)\n",
           asym, dg > 0.0 ? 100.0 * asym / dg : 0.0);
  }

  /* NC-D IS PRINTED ON EVERY ROW, NOT ON REQUEST. It is the error of the
   * TRIVIAL answer "there is no scatterer", i.e. the size of the scattered field
   * itself, and it is the denominator of every claim made here: the plan's
   * unconditional rule is that a result must beat its negative control at least
   * threefold. A field error of 1% means nothing until it is set against it. */
  printf("  %8s %14s %14s %14s %14s\n", "r/a", "rel L2 err", "NC-D no-scat", "|u| exact",
         "|u| got");
  for (int ir = 0; ir < 4; ir++) {
    double rr = c.a * (1.5 + 0.5 * (double)ir);
    if (rr > c.R - 2.0 * c.W) break;
    double num = 0.0, den = 0.0, ngot = 0.0, ncd = 0.0;
    for (int p = 0; p < 256; p++) {
      double th = 2.0 * M_PI * (double)p / 256.0;
      double x = rr * cos(th), y = rr * sin(th);
      double complex got = field_at(&c, b, dim, xs, x, y);
      double complex inc = cexp(CMPLX(0.0, 1.0) * c.k0 * x);
      double complex ex = inc + hz_mie_scat(&mie, mc, x, y);
      ngot += cabs(got) * cabs(got);
      num += cabs(got - ex) * cabs(got - ex);
      ncd += cabs(inc - ex) * cabs(inc - ex);
      den += cabs(ex) * cabs(ex);
    }
    printf("  %8.2f %14.4e %14.4e %14.4e %14.4e\n", rr / c.a, sqrt(num / den), sqrt(ncd / den),
           sqrt(den / 256.0), sqrt(ngot / 256.0));
  }

  /* --- F-N5: does the KERNEL lie in the ZERO FIELD? (PLAN врезка A5) -------
   * In vacuum both the LSQR answer and the plane wave satisfy the system, so
   * their difference d is a direction of the (numerical) kernel. The criterion
   * the architecture rests on — "an overcomplete basis is harmless because every
   * kernel direction carries no field" — is then a measurement, not an opinion:
   * compute the field of d. Measured before Nitsche: |Ad|/(|A||d|) = 8.4e-7 with
   * a field of 0.41/0.35/0.31/0.28 of the incident amplitude, i.e. the criterion
   * was VIOLATED. */
  if (argval(argc, argv, "a5", 0.0) > 0.0) {
    double complex *d = calloc((size_t)dim, sizeof(double complex));
    if (d) {
      double nd2 = 0.0, af2 = 0.0, ad2 = 0.0;
      for (int j = 0; j < dim; j++) {
        d[j] = xs[j] - ctrue[j];
        nd2 += cabs(d[j]) * cabs(d[j]);
      }
      for (size_t r = 0; r < (size_t)nrow; r++) {
        double complex s = 0.0;
        for (size_t p = rp[r]; p < rp[r + 1]; p++) {
          s += va[p] * d[ja[p]];
          af2 += cabs(va[p]) * cabs(va[p]);
        }
        ad2 += cabs(s) * cabs(s);
      }
      printf("  [F-N5] kernel test: |Ad|/(|A||d|) = %.3e   |d| = %.3e\n",
             sqrt(ad2) / (sqrt(af2) * sqrt(nd2)), sqrt(nd2));
      printf("  %8s %14s\n", "r/a", "field of d");
      for (int ir = 0; ir < 4; ir++) {
        double rr = c.a * (1.5 + 0.5 * (double)ir);
        if (rr > c.R - 2.0 * c.W) break;
        double s = 0.0;
        for (int p = 0; p < 256; p++) {
          double th = 2.0 * M_PI * (double)p / 256.0;
          double complex v = field_at(&c, b, dim, d, rr * cos(th), rr * sin(th));
          s += cabs(v) * cabs(v);
        }
        printf("  %8.2f %14.4e\n", rr / c.a, sqrt(s / 256.0));
      }
      free(d);
    }
  }

  free(ctrue);
  free(rhs);
  free(xs);
  free(uu);
  free(vv);
  free(ww);
  free(tv);
  free(ja);
  free(va);
  free(rp);
  free(rowsc);
  return 0;
}
