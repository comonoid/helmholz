#include "nitsche2d.h"
#include "cut2d.h"
#include "phi.h"
#include <complex.h>
#include <math.h>

/* Four phi factors of degree 2 each => degree 8; the derivative factors are one
 * lower and simply carry a zero top coefficient. */
enum { DMAX = 8, NPT = 20 };

static void pmul(const double complex *a, int da, const double complex *b, int db,
                 double complex *o) {
  for (int i = 0; i <= da + db; i++)
    o[i] = 0.0;
  for (int i = 0; i <= da; i++)
    for (int j = 0; j <= db; j++)
      o[i + j] += a[i] * b[j];
}

/* One phi factor as a polynomial in u on [-1,1], where the point on the segment
 * is coord = m + scale*u. hz_phi_local_poly gives it in s = coord - m, so the
 * k-th coefficient picks up scale^k. Composing into the SUB-PIECE parameter
 * rather than into a global coordinate is what keeps this conditioned (M14: the
 * same expansion taken about an origin a few units away read 1.1e-13 where the
 * local one reads 1e-15). */
static void fac_u(hz_phi_factor f, double m, double scale, double complex *p) {
  double q[3];
  hz_phi_local_poly(f, m, q);
  p[0] = q[0];
  p[1] = q[1] * scale;
  p[2] = q[2] * scale * scale;
}

/* Intersect the parameter window [*lo,*hi] with {a < p0 + dir*s < b}. */
static int clip_ax(double p0, double dir, double a, double b, double *lo, double *hi) {
  if (!(fabs(dir) > 0.0)) return p0 > a && p0 < b;
  double s0 = (a - p0) / dir, s1 = (b - p0) / dir;
  if (s0 > s1) {
    double t = s0;
    s0 = s1;
    s1 = t;
  }
  if (s0 > *lo) *lo = s0;
  if (s1 < *hi) *hi = s1;
  return 1;
}

/* Parameter values at which the factor's piece boundaries are crossed. phi has
 * knots at t = -2,-1,1,2 (t = 0 is not one: [-1,1] is a single parabola). */
static void add_knots(double W, int n, double p0, double dir, double lo, double hi, double *pt,
                      int *np) {
  static const double knot_t[4] = {-2.0, -1.0, 1.0, 2.0};
  if (!(fabs(dir) > 0.0)) return;
  for (int k = 0; k < 4; k++) {
    double s = (W * ((double)n + knot_t[k]) - p0) / dir;
    if (s > lo && s < hi && *np < NPT) pt[(*np)++] = s;
  }
}

hz_nit2d hz_nitsche2d_seg(hz_carrier2d bi, hz_carrier2d bj, double x0, double y0, double x1,
                          double y1, double nx, double ny) {
  hz_nit2d out = {0.0, 0.0, 0.0};
  double W = bi.W;
  double dx = x1 - x0, dy = y1 - y0;
  double L = sqrt(dx * dx + dy * dy);
  if (!(L > 0.0) || !(W > 0.0)) return out;
  double tx = dx / L, ty = dy / L;

  double lo = 0.0, hi = L;
  if (!clip_ax(x0, tx, W * ((double)bi.nx - 2.0), W * ((double)bi.nx + 2.0), &lo, &hi)) return out;
  if (!clip_ax(x0, tx, W * ((double)bj.nx - 2.0), W * ((double)bj.nx + 2.0), &lo, &hi)) return out;
  if (!clip_ax(y0, ty, W * ((double)bi.ny - 2.0), W * ((double)bi.ny + 2.0), &lo, &hi)) return out;
  if (!clip_ax(y0, ty, W * ((double)bj.ny - 2.0), W * ((double)bj.ny + 2.0), &lo, &hi)) return out;
  if (!(lo < hi)) return out;

  double pt[NPT];
  int np = 0;
  pt[np++] = lo;
  pt[np++] = hi;
  add_knots(W, bi.nx, x0, tx, lo, hi, pt, &np);
  add_knots(W, bj.nx, x0, tx, lo, hi, pt, &np);
  add_knots(W, bi.ny, y0, ty, lo, hi, pt, &np);
  add_knots(W, bj.ny, y0, ty, lo, hi, pt, &np);
  for (int i = 1; i < np; i++) { /* insertion sort */
    double v = pt[i];
    int j = i - 1;
    while (j >= 0 && pt[j] > v) {
      pt[j + 1] = pt[j];
      j--;
    }
    pt[j + 1] = v;
  }

  double complex omx = bi.kx + bj.kx, omy = bi.ky + bj.ky;
  for (int p = 0; p + 1 < np; p++) {
    double hw = 0.5 * (pt[p + 1] - pt[p]);
    if (!(hw > 0.0)) continue;
    double sm = 0.5 * (pt[p] + pt[p + 1]);
    double xm = x0 + tx * sm, ym = y0 + ty * sm;
    /* A segment parallel to an axis has scale = 0 on that axis, and then xm (or
     * ym) may land exactly on a knot. Harmless: phi is C^1 there, so the
     * deriv-0 and deriv-1 coefficients agree from both sides, and the only
     * coefficient that differs multiplies scale^2 == 0. */
    double complex Xi[3], Xj[3], Yi[3], Yj[3], Xid[3], Xjd[3], Yid[3], Yjd[3];
    hz_phi_factor f;
    f.h = W;
    f.n = (double)bi.nx;
    f.deriv = 0;
    fac_u(f, xm, tx * hw, Xi);
    f.deriv = 1;
    fac_u(f, xm, tx * hw, Xid);
    f.n = (double)bj.nx;
    f.deriv = 0;
    fac_u(f, xm, tx * hw, Xj);
    f.deriv = 1;
    fac_u(f, xm, tx * hw, Xjd);
    f.n = (double)bi.ny;
    f.deriv = 0;
    fac_u(f, ym, ty * hw, Yi);
    f.deriv = 1;
    fac_u(f, ym, ty * hw, Yid);
    f.n = (double)bj.ny;
    f.deriv = 0;
    fac_u(f, ym, ty * hw, Yj);
    f.deriv = 1;
    fac_u(f, ym, ty * hw, Yjd);

    double complex A[5], B[5], Axj[5], Axi[5], Byj[5], Byi[5];
    double complex P00[DMAX + 1], Pxj[DMAX + 1], Pyj[DMAX + 1], Pxi[DMAX + 1], Pyi[DMAX + 1];
    pmul(Xi, 2, Xj, 2, A);
    pmul(Yi, 2, Yj, 2, B);
    pmul(Xi, 2, Xjd, 2, Axj);
    pmul(Xid, 2, Xj, 2, Axi);
    pmul(Yi, 2, Yjd, 2, Byj);
    pmul(Yid, 2, Yj, 2, Byi);
    pmul(A, 4, B, 4, P00);
    pmul(Axj, 4, B, 4, Pxj);
    pmul(Axi, 4, B, 4, Pxi);
    pmul(A, 4, Byj, 4, Pyj);
    pmul(A, 4, Byi, 4, Pyi);

    double complex mom[DMAX + 1];
    hz_cut2d_osc_moments((omx * tx + omy * ty) * hw, DMAX, mom);
    double complex i00 = 0.0, ixj = 0.0, iyj = 0.0, ixi = 0.0, iyi = 0.0;
    for (int k = 0; k <= DMAX; k++) {
      i00 += P00[k] * mom[k];
      ixj += Pxj[k] * mom[k];
      iyj += Pyj[k] * mom[k];
      ixi += Pxi[k] * mom[k];
      iyi += Pyi[k] * mom[k];
    }
    /* Phase kept in the LOCAL reference of each factor (|x - n W| <= 2W), never
     * as a global omega*x that would be a huge angle far from the origin. */
    double complex ph = bi.kx * (xm - W * (double)bi.nx) + bj.kx * (xm - W * (double)bj.nx) +
                        bi.ky * (ym - W * (double)bi.ny) + bj.ky * (ym - W * (double)bj.ny);
    double complex E = hw * cexp(CMPLX(0.0, 1.0) * ph);
    double complex i1 = CMPLX(0.0, 1.0);
    out.t0 += E * i00;
    out.t1j += E * (nx * (ixj + i1 * bj.kx * i00) + ny * (iyj + i1 * bj.ky * i00));
    out.t1i += E * (nx * (ixi + i1 * bi.kx * i00) + ny * (iyi + i1 * bi.ky * i00));
  }
  return out;
}

hz_nit2d hz_nitsche2d_poly(hz_carrier2d bi, hz_carrier2d bj, const double *vx, const double *vy,
                           int nv) {
  hz_nit2d out = {0.0, 0.0, 0.0};
  double W = bi.W;
  double xlo = W * ((double)(bi.nx > bj.nx ? bi.nx : bj.nx) - 2.0);
  double xhi = W * ((double)(bi.nx < bj.nx ? bi.nx : bj.nx) + 2.0);
  double ylo = W * ((double)(bi.ny > bj.ny ? bi.ny : bj.ny) - 2.0);
  double yhi = W * ((double)(bi.ny < bj.ny ? bi.ny : bj.ny) + 2.0);
  if (!(xlo < xhi) || !(ylo < yhi)) return out;
  for (int e = 0; e < nv; e++) {
    int f = (e + 1) % nv;
    double ax = vx[e], ay = vy[e], bx = vx[f], by = vy[f];
    /* box reject before any arithmetic on the facet: most facets of the ring
     * are nowhere near the overlap of two supports */
    if ((ax < xlo && bx < xlo) || (ax > xhi && bx > xhi)) continue;
    if ((ay < ylo && by < ylo) || (ay > yhi && by > yhi)) continue;
    double dx = bx - ax, dy = by - ay;
    double L = sqrt(dx * dx + dy * dy);
    if (!(L > 0.0)) continue;
    hz_nit2d s = hz_nitsche2d_seg(bi, bj, ax, ay, bx, by, dy / L, -dx / L);
    out.t0 += s.t0;
    out.t1j += s.t1j;
    out.t1i += s.t1i;
  }
  return out;
}
