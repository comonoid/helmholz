#include "carrier.h"
#include "phi.h"
#include <math.h>
#include <stddef.h>

/* Half-line this function lives on, given its cut side. */
static void side_range(int side, double xcut, double *lo, double *hi) {
  double big = 1e300;
  *lo = (side == 1) ? xcut : -big;
  *hi = (side == 0) ? xcut : big;
}

double complex hz_carrier_val(hz_carrier b, double x, double xcut) {
  double lo, hi;
  side_range(b.side, xcut, &lo, &hi);
  if (x < lo || x > hi) return 0.0;
  double t = x / b.W - (double)b.n;
  if (fabs(t) >= 2.0) return 0.0;
  double p = hz_phi(t);
  /* phase referenced to the ELEMENT, not to the origin — see the note in
   * carrier.h on why this is mandatory once kx is complex */
  return p * cexp(CMPLX(0.0, 1.0) * b.kx * (x - b.W * (double)b.n));
}

void hz_carrier_overlap(hz_carrier bi, hz_carrier bj, double dom, double xcut, double *lo,
                        double *hi) {
  double l1 = bi.W * ((double)bi.n - 2.0), h1 = bi.W * ((double)bi.n + 2.0);
  double l2 = bj.W * ((double)bj.n - 2.0), h2 = bj.W * ((double)bj.n + 2.0);
  double si_lo, si_hi, sj_lo, sj_hi;
  side_range(bi.side, xcut, &si_lo, &si_hi);
  side_range(bj.side, xcut, &sj_lo, &sj_hi);

  double l = l1 > l2 ? l1 : l2;
  if (si_lo > l) l = si_lo;
  if (sj_lo > l) l = sj_lo;
  if (l < 0.0) l = 0.0;

  double h = h1 < h2 ? h1 : h2;
  if (si_hi < h) h = si_hi;
  if (sj_hi < h) h = sj_hi;
  if (h > dom) h = dom;

  *lo = l;
  *hi = h;
}

double complex hz_carrier_entry(hz_carrier bi, hz_carrier bj, double dom, double xcut,
                                const hz_medseg *segs, int nseg, long *ncalls) {
  double lo, hi;
  hz_carrier_overlap(bi, bj, dom, xcut, &lo, &hi);
  if (lo >= hi) return 0.0;
  double complex om = bi.kx + bj.kx; /* bilinear form: the exponents just add */
  /* Local phase reference: each factor carries e^{i kx (x - n W)}, so the pair
   * contributes a constant prefactor on top of the e^{i om x} integral. For a
   * real carrier this is a pure phase and harmless; for a COMPLEX one it is a
   * rescaling by e^{Im(kx) n W} that keeps every element O(1) instead of
   * e^{-Im(kx) x}. Without it an absorbing layer at x ~ 1000 lambda scales its
   * own basis functions by e^{-175} and the system is unusable. */
  double complex pref =
      cexp(CMPLX(0.0, -1.0) * (bi.kx * bi.W * (double)bi.n + bj.kx * bj.W * (double)bj.n));
  hz_phi_factor fi0 = {bi.W, (double)bi.n, 0};
  hz_phi_factor fj0 = {bj.W, (double)bj.n, 0};
  hz_phi_factor fj1 = {bj.W, (double)bj.n, 1};
  hz_phi_factor fj2 = {bj.W, (double)bj.n, 2};

  /* (d2/dx2 + k^2) B_j = [phi'' + 2i kx phi' + (k^2(x) - kx^2) phi] e^{i kx x};
   * the deriv factors already carry their 1/W^d. The medium term is measured
   * against THIS ELEMENT'S OWN kx^2 — that is what "carrier local to the
   * medium" means operationally, and it is why there is no global background
   * wavenumber anywhere in this file. */
  if (ncalls != NULL) *ncalls += 2;
  double complex total = hz_phi_prod_integral_osc(lo, hi, fi0, fj2, om);
  total += CMPLX(0.0, 2.0) * bj.kx * hz_phi_prod_integral_osc(lo, hi, fi0, fj1, om);

  double complex kx2 = bj.kx * bj.kx;
  for (int s = 0; s < nseg; s++) {
    double a = lo > segs[s].a ? lo : segs[s].a;
    double b = hi < segs[s].b ? hi : segs[s].b;
    if (a >= b) continue;
    double complex dk2 = segs[s].k2 - kx2;
    if (!(cabs(dk2) > 0.0)) continue; /* carrier matches the medium: term cancels */
    if (ncalls != NULL) *ncalls += 1;
    total += dk2 * hz_phi_prod_integral_osc(a, b, fi0, fj0, om);
  }
  return pref * total;
}

double complex hz_carrier_rhs(hz_carrier bi, double dom, double xcut, double xsrc, double ws,
                              long *ncalls) {
  hz_phi_factor fi0 = {bi.W, (double)bi.n, 0};
  hz_phi_factor fsrc = {ws, xsrc / ws, 0};
  double l1 = bi.W * ((double)bi.n - 2.0), h1 = bi.W * ((double)bi.n + 2.0);
  double si_lo, si_hi;
  side_range(bi.side, xcut, &si_lo, &si_hi);
  double lo = l1 > xsrc - 2.0 * ws ? l1 : xsrc - 2.0 * ws;
  double hi = h1 < xsrc + 2.0 * ws ? h1 : xsrc + 2.0 * ws;
  if (si_lo > lo) lo = si_lo;
  if (si_hi < hi) hi = si_hi;
  if (lo < 0.0) lo = 0.0;
  if (hi > dom) hi = dom;
  if (lo >= hi) return 0.0;
  if (ncalls != NULL) *ncalls += 1;
  return cexp(CMPLX(0.0, -1.0) * bi.kx * bi.W * (double)bi.n) *
         hz_phi_prod_integral_osc(lo, hi, fi0, fsrc, bi.kx);
}
