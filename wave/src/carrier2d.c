#include "carrier2d.h"
#include "phi.h"
#include <math.h>
#include <stddef.h>

/* One axis of the tensor product. Returns Int phi^(di)(x/Wi - ni)
 * phi^(dj)(x/Wj - nj) e^{i om x} dx over [lo,hi], with the LOCAL phase
 * prefactor left to the caller (it is common to all four terms). */
static double complex axis_int(double lo, double hi, double Wi, int ni, int di, double Wj, int nj,
                               int dj, double complex om, long *ncalls) {
  if (lo >= hi) return 0.0;
  hz_phi_factor fi = {Wi, (double)ni, di};
  hz_phi_factor fj = {Wj, (double)nj, dj};
  if (ncalls != NULL) *ncalls += 1;
  return hz_phi_prod_integral_osc(lo, hi, fi, fj, om);
}

/* Overlap of the two supports on one axis, clipped to [0,dom]. */
static void axis_overlap(double Wi, int ni, double Wj, int nj, double dom, double *lo, double *hi) {
  double l1 = Wi * ((double)ni - 2.0), h1 = Wi * ((double)ni + 2.0);
  double l2 = Wj * ((double)nj - 2.0), h2 = Wj * ((double)nj + 2.0);
  double l = l1 > l2 ? l1 : l2;
  double h = h1 < h2 ? h1 : h2;
  if (l < 0.0) l = 0.0;
  if (h > dom) h = dom;
  *lo = l;
  *hi = h;
}

double complex hz_carrier2d_val(hz_carrier2d b, double x, double y) {
  double tx = x / b.W - (double)b.nx, ty = y / b.W - (double)b.ny;
  if (fabs(tx) >= 2.0 || fabs(ty) >= 2.0) return 0.0;
  double complex i = CMPLX(0.0, 1.0);
  return hz_phi(tx) * hz_phi(ty) *
         cexp(i * (b.kx * (x - b.W * (double)b.nx) + b.ky * (y - b.W * (double)b.ny)));
}

double complex hz_carrier2d_entry(hz_carrier2d bi, hz_carrier2d bj, double domx, double domy,
                                  const hz_medrect2d *rects, int nrect, long *ncalls) {
  double xlo, xhi, ylo, yhi;
  axis_overlap(bi.W, bi.nx, bj.W, bj.nx, domx, &xlo, &xhi);
  axis_overlap(bi.W, bi.ny, bj.W, bj.ny, domy, &ylo, &yhi);
  if (xlo >= xhi || ylo >= yhi) return 0.0;

  double complex omx = bi.kx + bj.kx; /* bilinear form: exponents add */
  double complex omy = bi.ky + bj.ky;
  double complex i = CMPLX(0.0, 1.0);
  /* local phase reference of both factors on both axes, common prefactor */
  double complex pref = cexp(-i * (bi.kx * bi.W * (double)bi.nx + bj.kx * bj.W * (double)bj.nx +
                                   bi.ky * bi.W * (double)bi.ny + bj.ky * bj.W * (double)bj.ny));

  /* (d2/dx2 + d2/dy2 + k^2) B_j, with the k^2 term measured against THIS
   * ELEMENT'S OWN carrier - that is what "carrier local to the medium" means
   * operationally, and why there is no global background wavenumber here. */
  double complex Ix0 = axis_int(xlo, xhi, bi.W, bi.nx, 0, bj.W, bj.nx, 0, omx, ncalls);
  double complex Iy0 = axis_int(ylo, yhi, bi.W, bi.ny, 0, bj.W, bj.ny, 0, omy, ncalls);
  double complex Ix1 = axis_int(xlo, xhi, bi.W, bi.nx, 0, bj.W, bj.nx, 1, omx, ncalls);
  double complex Iy1 = axis_int(ylo, yhi, bi.W, bi.ny, 0, bj.W, bj.ny, 1, omy, ncalls);
  double complex Ix2 = axis_int(xlo, xhi, bi.W, bi.nx, 0, bj.W, bj.nx, 2, omx, ncalls);
  double complex Iy2 = axis_int(ylo, yhi, bi.W, bi.ny, 0, bj.W, bj.ny, 2, omy, ncalls);

  /* d2/dx2 [phi(x) e^{i kx x}] = (phi'' + 2 i kx phi' - kx^2 phi) e^{i kx x};
   * the -kx^2 and -ky^2 pieces are folded into the medium term below so that
   * they cancel there exactly when the carrier matches the medium. */
  double complex total = (Ix2 + CMPLX(0.0, 2.0) * bj.kx * Ix1) * Iy0;
  total += Ix0 * (Iy2 + CMPLX(0.0, 2.0) * bj.ky * Iy1);

  double complex kj2 = bj.kx * bj.kx + bj.ky * bj.ky;
  for (int r = 0; r < nrect; r++) {
    double a = xlo > rects[r].x0 ? xlo : rects[r].x0;
    double b = xhi < rects[r].x1 ? xhi : rects[r].x1;
    double c = ylo > rects[r].y0 ? ylo : rects[r].y0;
    double d = yhi < rects[r].y1 ? yhi : rects[r].y1;
    if (a >= b || c >= d) continue;
    double complex dk2 = rects[r].k2 - kj2;
    if (!(cabs(dk2) > 0.0)) continue; /* carrier matches the medium: cancels */
    double complex Jx = axis_int(a, b, bi.W, bi.nx, 0, bj.W, bj.nx, 0, omx, ncalls);
    double complex Jy = axis_int(c, d, bi.W, bi.ny, 0, bj.W, bj.ny, 0, omy, ncalls);
    total += dk2 * Jx * Jy;
  }
  return pref * total;
}

double complex hz_carrier2d_rhs(hz_carrier2d bi, double domx, double domy, double xs, double ys,
                                double ws, long *ncalls) {
  /* the source is a tensor product too, so the overlap is taken per axis
   * against the source support rather than against another element */
  double xlo = bi.W * ((double)bi.nx - 2.0);
  double xhi = bi.W * ((double)bi.nx + 2.0);
  double ylo, yhi;
  if (xlo < xs - 2.0 * ws) xlo = xs - 2.0 * ws;
  if (xhi > xs + 2.0 * ws) xhi = xs + 2.0 * ws;
  if (xlo < 0.0) xlo = 0.0;
  if (xhi > domx) xhi = domx;
  ylo = bi.W * ((double)bi.ny - 2.0);
  yhi = bi.W * ((double)bi.ny + 2.0);
  if (ylo < ys - 2.0 * ws) ylo = ys - 2.0 * ws;
  if (yhi > ys + 2.0 * ws) yhi = ys + 2.0 * ws;
  if (ylo < 0.0) ylo = 0.0;
  if (yhi > domy) yhi = domy;
  if (xlo >= xhi || ylo >= yhi) return 0.0;

  double complex i = CMPLX(0.0, 1.0);
  double complex pref = cexp(-i * (bi.kx * bi.W * (double)bi.nx + bi.ky * bi.W * (double)bi.ny));
  hz_phi_factor fx = {bi.W, (double)bi.nx, 0}, fsx = {ws, xs / ws, 0};
  hz_phi_factor fy = {bi.W, (double)bi.ny, 0}, fsy = {ws, ys / ws, 0};
  if (ncalls != NULL) *ncalls += 2;
  return pref * hz_phi_prod_integral_osc(xlo, xhi, fx, fsx, bi.kx) *
         hz_phi_prod_integral_osc(ylo, yhi, fy, fsy, bi.ky);
}
