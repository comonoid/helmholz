/* Unit tests for src/nitsche2d.c — the interface (Nitsche) line integrals.
 *
 * WHY THESE COME BEFORE ANY PHYSICS. The whole point of the Nitsche term is that
 * it removes a sampling knob; if its integrals are wrong, the run that follows
 * will still produce a plausible number and it will be read as a statement about
 * the method. So the three integrals get their own falsifiers first
 * (M15 audit §3, F-N6/F-N7):
 *   N1  the polygon walked here is the SAME polygon the basis is cut by —
 *       otherwise the jump is measured where the function is continuous;
 *   N2  against an independent witness (composite Gauss on sampled phi and
 *       phi', sharing no code with the piecewise/moment machinery);
 *   N3  additivity of a split facet — the only check that survives at kW >> 1,
 *       where quadrature cannot testify at all (M9a);
 *   N4  strict zero when the facet misses the supports;
 *   N5  GREEN'S IDENTITY against src/cut2d.c: the asymmetry the strong form
 *       leaves on a cut pair must equal exactly the interface term. This is the
 *       derivation itself under test, and it is the one check that fixes the
 *       SIGNS and the direction of the normal;
 *   N6  parametrisation: reversing a facet must change nothing, flipping its
 *       normal must flip the two normal-derivative terms and leave the first.
 *
 * ERRORS ARE NORMALISED BY THE NATURAL SCALE, never by the result: at k W >> 1
 * these integrals are small by oscillatory cancellation, and dividing by the
 * result would report the dynamic range of the integrand instead of an error
 * (the trap that cost M14 a reading). */
#include "cut2d.h"
#include "nitsche2d.h"
#include "phi.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>

static int pass_count = 0, fail_count = 0;

static void check(const char *name, int ok, double v) {
  if (ok) {
    pass_count++;
    printf("  ok   %-54s %.3e\n", name, v);
    return;
  }
  fail_count++;
  printf("  FAIL %-54s %.3e\n", name, v);
}

/* --- independent witness: sample the basis and integrate with composite Gauss */
static double complex bval(hz_carrier2d b, double x, double y) {
  double tx = x / b.W - (double)b.nx, ty = y / b.W - (double)b.ny;
  double complex ph = b.kx * (x - b.W * (double)b.nx) + b.ky * (y - b.W * (double)b.ny);
  return hz_phi(tx) * hz_phi(ty) * cexp(CMPLX(0.0, 1.0) * ph);
}

static void bgrad(hz_carrier2d b, double x, double y, double complex *gx, double complex *gy) {
  double tx = x / b.W - (double)b.nx, ty = y / b.W - (double)b.ny;
  double complex ph = b.kx * (x - b.W * (double)b.nx) + b.ky * (y - b.W * (double)b.ny);
  double complex e = cexp(CMPLX(0.0, 1.0) * ph);
  double px = hz_phi(tx), py = hz_phi(ty);
  double complex i1 = CMPLX(0.0, 1.0);
  *gx = (hz_phi_d1(tx) / b.W + i1 * b.kx * px) * py * e;
  *gy = px * (hz_phi_d1(ty) / b.W + i1 * b.ky * py) * e;
}

/* 4-point Gauss-Legendre, NPAN panels. phi is only C^1, so a panel containing a
 * knot converges at O(h^3) and not at the Gauss rate; NPAN is set so that even
 * those panels sit far below the thresholds used here. */
enum { NPAN = 20000 };
static void witness(hz_carrier2d bi, hz_carrier2d bj, double x0, double y0, double x1, double y1,
                    double nx, double ny, hz_nit2d *out) {
  static const double gx[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
                               0.8611363115940526};
  static const double gw[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
                               0.3478548451374538};
  double dx = x1 - x0, dy = y1 - y0;
  double L = sqrt(dx * dx + dy * dy);
  out->t0 = out->t1i = out->t1j = 0.0;
  for (int p = 0; p < NPAN; p++) {
    double ta = (double)p / (double)NPAN, tb = (double)(p + 1) / (double)NPAN;
    double c = 0.5 * (ta + tb), h = 0.5 * (tb - ta);
    for (int q = 0; q < 4; q++) {
      double t = c + h * gx[q], w = gw[q] * h * L;
      double x = x0 + dx * t, y = y0 + dy * t;
      double complex vi = bval(bi, x, y), vj = bval(bj, x, y);
      double complex gix, giy, gjx, gjy;
      bgrad(bi, x, y, &gix, &giy);
      bgrad(bj, x, y, &gjx, &gjy);
      out->t0 += w * vi * vj;
      out->t1j += w * vi * (nx * gjx + ny * gjy);
      out->t1i += w * vj * (nx * gix + ny * giy);
    }
  }
}

/* --- the strong-form Galerkin entry over one side of the cut, via src/cut2d.c.
 * Deliberately written from the definition here rather than shared with the
 * bench: N5 wants an INDEPENDENT witness for the volume half of Green's
 * identity, and the cut2d path it uses is separately gated (23/23). */
static double complex strong(hz_carrier2d bi, hz_carrier2d bj, int inside, const hz_half2 *hp,
                             int nhp) {
  double W = bi.W;
  double xlo = W * ((double)(bi.nx > bj.nx ? bi.nx : bj.nx) - 2.0);
  double xhi = W * ((double)(bi.nx < bj.nx ? bi.nx : bj.nx) + 2.0);
  double ylo = W * ((double)(bi.ny > bj.ny ? bi.ny : bj.ny) - 2.0);
  double yhi = W * ((double)(bi.ny < bj.ny ? bi.ny : bj.ny) + 2.0);
  if (!(xlo < xhi) || !(ylo < yhi)) return 0.0;
  double complex omx = bi.kx + bj.kx, omy = bi.ky + bj.ky, i1 = CMPLX(0.0, 1.0);
  hz_phi_factor fi0 = {W, (double)bi.nx, 0}, fj0 = {W, (double)bj.nx, 0};
  hz_phi_factor gi0 = {W, (double)bi.ny, 0}, gj0 = {W, (double)bj.ny, 0};
  hz_phi_factor fj1 = {W, (double)bj.nx, 1}, fj2 = {W, (double)bj.nx, 2};
  hz_phi_factor gj1 = {W, (double)bj.ny, 1}, gj2 = {W, (double)bj.ny, 2};
  hz_axis2 X0 = {{fi0, fj0}, 2}, X1 = {{fi0, fj1}, 2}, X2 = {{fi0, fj2}, 2};
  hz_axis2 Y0 = {{gi0, gj0}, 2}, Y1 = {{gi0, gj1}, 2}, Y2 = {{gi0, gj2}, 2};
  hz_axis2 ax[4] = {X2, X1, X0, X0}, ay[4] = {Y0, Y0, Y2, Y1};
  double complex cf[4] = {1.0, 2.0 * i1 * bj.kx, 1.0, 2.0 * i1 * bj.ky};
  double complex t = 0.0;
  for (int k = 0; k < 4; k++) {
    double complex v;
    if (inside) {
      v = hz_cut2d_poly(xlo, xhi, ylo, yhi, ax[k], ay[k], omx, omy, hp, nhp);
    } else {
      v = hz_cut2d_poly(xlo, xhi, ylo, yhi, ax[k], ay[k], omx, omy, hp, 0) -
          hz_cut2d_poly(xlo, xhi, ylo, yhi, ax[k], ay[k], omx, omy, hp, nhp);
    }
    t += cf[k] * v;
  }
  double complex pref = cexp(-i1 * (bi.kx * W * (double)bi.nx + bj.kx * W * (double)bj.nx +
                                    bi.ky * W * (double)bi.ny + bj.ky * W * (double)bj.ny));
  return pref * t;
}

int main(void) {
  printf("test_nitsche2d: interface integrals for the cut\n");

  /* ---------------- N1: the interface polygon IS the cut polygon ---------- */
  {
    double worst_on = 0.0, worst_in = -1e300, worst_out = -1e300;
    for (int fit = 0; fit < 2; fit++)
      for (int m = 3; m <= 48; m += 5) {
        hz_half2 hp[HZ_CUT_MAXHP];
        double vx[HZ_CUT_MAXHP], vy[HZ_CUT_MAXHP];
        int nh = hz_cut2d_disc(0.3, -0.4, 1.7, m, fit, hp);
        int nv = hz_cut2d_disc_verts(0.3, -0.4, 1.7, m, fit, vx, vy);
        if (nv != nh) {
          check("N1 vertex count matches facet count", 0, (double)(nv - nh));
          continue;
        }
        for (int j = 0; j < nv; j++) {
          for (int i = 0; i < nh; i++) {
            double s = vx[j] * hp[i].ca + vy[j] * hp[i].sa - hp[i].c;
            int adj = (i == j) || (i == (j + nh - 1) % nh);
            if (adj) {
              if (fabs(s) > worst_on) worst_on = fabs(s);
            } else if (s > worst_out) {
              worst_out = s;
            }
          }
          /* edge midpoints must be strictly inside every other facet */
          int f = (j + 1) % nv;
          double mx = 0.5 * (vx[j] + vx[f]), my = 0.5 * (vy[j] + vy[f]);
          double s = mx * hp[j].ca + my * hp[j].sa - hp[j].c;
          if (s > worst_in) worst_in = s;
          /* the outward normal (dy,-dx) must point away from the centre */
          double dx = vx[f] - vx[j], dy = vy[f] - vy[j];
          double L = sqrt(dx * dx + dy * dy);
          double dot = (mx - 0.3) * (dy / L) + (my + 0.4) * (-dx / L);
          if (dot <= 0.0) check("N1 edge normal points outward", 0, dot);
        }
      }
    check("N1 vertices lie ON their two facets", worst_on < 1e-14, worst_on);
    check("N1 vertices lie inside all other facets", worst_out < -1e-9, worst_out);
    check("N1 edge midpoint lies ON its own facet", fabs(worst_in) < 1e-15, fabs(worst_in));
    check("N1 outward normals verified on every edge", 1, 0.0);
  }

  /* ---------------- N2: against the independent witness ------------------ */
  {
    double W = 1.0;
    /* an oblique facet that genuinely crosses both supports and several knots */
    double x0 = -1.35, y0 = -0.9, x1 = 1.15, y1 = 1.6;
    double dx = x1 - x0, dy = y1 - y0, L = sqrt(dx * dx + dy * dy);
    double nx = dy / L, ny = -dx / L;
    const double kk[3] = {0.0, 1.3, 4.0};
    for (int c = 0; c < 3; c++) {
      double k = kk[c];
      hz_carrier2d bi = {W, 0, 0, k * cos(0.4), k * sin(0.4)};
      hz_carrier2d bj = {W, 1, -1, k * cos(2.1), k * sin(2.1)};
      hz_nit2d got = hz_nitsche2d_seg(bi, bj, x0, y0, x1, y1, nx, ny);
      hz_nit2d ref;
      witness(bi, bj, x0, y0, x1, y1, nx, ny, &ref);
      double sc0 = 16.0 * L, sc1 = 16.0 * L * (k + 1.0 / W);
      char nm[64];
      snprintf(nm, sizeof nm, "N2 t0  vs quadrature, kW = %.1f", k * W);
      check(nm, cabs(got.t0 - ref.t0) / sc0 < 1e-12, cabs(got.t0 - ref.t0) / sc0);
      snprintf(nm, sizeof nm, "N2 t1j vs quadrature, kW = %.1f", k * W);
      check(nm, cabs(got.t1j - ref.t1j) / sc1 < 1e-12, cabs(got.t1j - ref.t1j) / sc1);
      snprintf(nm, sizeof nm, "N2 t1i vs quadrature, kW = %.1f", k * W);
      check(nm, cabs(got.t1i - ref.t1i) / sc1 < 1e-12, cabs(got.t1i - ref.t1i) / sc1);
    }
  }

  /* ---------------- N3: additivity of a split facet ---------------------- */
  {
    double W = 1.0;
    double x0 = -1.35, y0 = -0.9, x1 = 1.15, y1 = 1.6;
    double dx = x1 - x0, dy = y1 - y0, L = sqrt(dx * dx + dy * dy);
    double nx = dy / L, ny = -dx / L;
    const double kk[4] = {0.0, 2.0, 20.0, 200.0};
    for (int c = 0; c < 4; c++) {
      double k = kk[c];
      hz_carrier2d bi = {W, 0, 0, k * cos(0.4), k * sin(0.4)};
      hz_carrier2d bj = {W, 1, -1, k * cos(2.1), k * sin(2.1)};
      double tsp = 0.37; /* not a knot crossing, on purpose */
      double mx = x0 + dx * tsp, my = y0 + dy * tsp;
      hz_nit2d w = hz_nitsche2d_seg(bi, bj, x0, y0, x1, y1, nx, ny);
      hz_nit2d a = hz_nitsche2d_seg(bi, bj, x0, y0, mx, my, nx, ny);
      hz_nit2d b = hz_nitsche2d_seg(bi, bj, mx, my, x1, y1, nx, ny);
      double sc0 = 16.0 * L, sc1 = 16.0 * L * (k + 1.0 / W);
      double e0 = cabs(w.t0 - a.t0 - b.t0) / sc0;
      double e1 = (cabs(w.t1j - a.t1j - b.t1j) + cabs(w.t1i - a.t1i - b.t1i)) / sc1;
      char nm[64];
      snprintf(nm, sizeof nm, "N3 additivity of t0,  kW = %.0f", k * W);
      check(nm, e0 < 1e-14, e0);
      snprintf(nm, sizeof nm, "N3 additivity of t1,  kW = %.0f", k * W);
      check(nm, e1 < 1e-14, e1);
    }
  }

  /* ---------------- N4: strict zero off the supports --------------------- */
  {
    double W = 1.0;
    hz_carrier2d bi = {W, 0, 0, 1.0, 0.5};
    hz_carrier2d bj = {W, 1, 0, -0.7, 1.1};
    hz_nit2d far = hz_nitsche2d_seg(bi, bj, 9.0, -3.0, 9.0, 3.0, 1.0, 0.0);
    double sfar = cabs(far.t0) + cabs(far.t1i) + cabs(far.t1j);
    check("N4 facet outside both supports gives exact 0", !(sfar > 0.0), sfar);
    hz_carrier2d bk = {W, 9, 0, 0.3, 0.2}; /* supports do not overlap at all */
    hz_nit2d dis =
        hz_nitsche2d_seg(bi, bk, -3.0, -3.0, 3.0, 3.0, 0.7071067811865476, -0.7071067811865476);
    double sdis = cabs(dis.t0) + cabs(dis.t1i) + cabs(dis.t1j);
    check("N4 disjoint supports give exact 0", !(sdis > 0.0), sdis);
  }

  /* ---------------- N5: Green's identity against src/cut2d.c ------------- */
  {
    /* The strong form loses exactly one surface term on a cut pair:
     *   Strong_ij = -a_broken(B_i,B_j) - Int_G ([B_i]{dn B_j} + {B_i}[dn B_j])
     * with a_broken symmetric. With traces [v] = alpha v|G, {v} = beta v|G
     * (alpha,beta) = (sigma,1/2) for a cut function, that says
     *   Strong_ij - Strong_ji = K (t1i - t1j),  K = alpha_i beta_j + beta_i alpha_j
     * i.e. K = -1 for a pair on the inside and +1 for a pair on the outside.
     * Nothing here is fitted: both sides are computed and compared. */
    double W = 1.0;
    hz_half2 hp[HZ_CUT_MAXHP];
    double vx[HZ_CUT_MAXHP], vy[HZ_CUT_MAXHP];
    int nh = hz_cut2d_disc(0.15, -0.2, 1.4, 12, 1, hp);
    hz_cut2d_disc_verts(0.15, -0.2, 1.4, 12, 1, vx, vy);
    const double kk[3] = {0.0, 1.7, 12.0};
    const int dn[3][2] = {{0, 0}, {1, 0}, {1, -1}};
    for (int c = 0; c < 3; c++)
      for (int d = 0; d < 3; d++)
        for (int side = 0; side < 2; side++) {
          double k = kk[c];
          /* SAME |k| on both functions: they share a region, so they share the
           * medium, and only then is the k^2 term of the strong form symmetric */
          hz_carrier2d bi = {W, 0, 0, k * cos(0.4), k * sin(0.4)};
          hz_carrier2d bj = {W, dn[d][0], dn[d][1], k * cos(2.3), k * sin(2.3)};
          double complex sij = strong(bi, bj, side == 0, hp, nh);
          double complex sji = strong(bj, bi, side == 0, hp, nh);
          hz_nit2d t = hz_nitsche2d_poly(bi, bj, vx, vy, nh);
          double K = (side == 0) ? -1.0 : 1.0;
          double sc = 16.0 * 4.0 * W * (k + 1.0 / W);
          double e = cabs((sij - sji) - K * (t.t1i - t.t1j)) / sc;
          char nm[80];
          snprintf(nm, sizeof nm, "N5 Green, %s, dn=(%d,%d), kW=%.0f", side ? "outside" : "inside",
                   dn[d][0], dn[d][1], k * W);
          check(nm, e < 1e-12, e);
        }
  }

  /* ---------------- N6: parametrisation and normal ----------------------- */
  {
    double W = 1.0, k = 3.0;
    hz_carrier2d bi = {W, 0, 0, k * cos(0.4), k * sin(0.4)};
    hz_carrier2d bj = {W, 1, -1, k * cos(2.1), k * sin(2.1)};
    double x0 = -1.35, y0 = -0.9, x1 = 1.15, y1 = 1.6;
    double dx = x1 - x0, dy = y1 - y0, L = sqrt(dx * dx + dy * dy);
    double nx = dy / L, ny = -dx / L;
    hz_nit2d a = hz_nitsche2d_seg(bi, bj, x0, y0, x1, y1, nx, ny);
    hz_nit2d r = hz_nitsche2d_seg(bi, bj, x1, y1, x0, y0, nx, ny);
    hz_nit2d f = hz_nitsche2d_seg(bi, bj, x0, y0, x1, y1, -nx, -ny);
    double sc0 = 16.0 * L, sc1 = 16.0 * L * (k + 1.0 / W);
    double er = (cabs(a.t0 - r.t0) / sc0 + cabs(a.t1j - r.t1j) / sc1 + cabs(a.t1i - r.t1i) / sc1);
    double ef = (cabs(a.t0 - f.t0) / sc0 + cabs(a.t1j + f.t1j) / sc1 + cabs(a.t1i + f.t1i) / sc1);
    check("N6 reversing the facet changes nothing", er < 1e-15, er);
    check("N6 flipping the normal flips t1 only", ef < 1e-15, ef);
    /* i<->j swap must exchange the two normal-derivative terms */
    hz_nit2d s = hz_nitsche2d_seg(bj, bi, x0, y0, x1, y1, nx, ny);
    double es = (cabs(a.t0 - s.t0) / sc0 + cabs(a.t1j - s.t1i) / sc1 + cabs(a.t1i - s.t1j) / sc1);
    check("N6 swapping i and j exchanges t1i and t1j", es < 1e-15, es);
  }

  printf("test_nitsche2d: %d/%d passed\n", pass_count, pass_count + fail_count);
  return fail_count == 0 ? 0 : 1;
}
