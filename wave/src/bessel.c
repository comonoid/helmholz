#include "bessel.h"
#include <math.h>

void hz_bessel_jy(int m, double x, double *j, double *y, double *dj, double *dy) {
  double jm = jn(m, x), ym = yn(m, x);
  double jp = jn(m + 1, x), yp = yn(m + 1, x);
  /* f'_m = f_{m-1} - (m/x) f_m, which avoids naming f_{-1} for m = 0 (there
   * f'_0 = -f_1 comes out of the same expression). */
  double jmm = (m == 0) ? -jp : jn(m - 1, x);
  double ymm = (m == 0) ? -yp : yn(m - 1, x);
  *j = jm;
  *y = ym;
  *dj = (m == 0) ? jmm : 0.5 * (jmm - jp);
  *dy = (m == 0) ? ymm : 0.5 * (ymm - yp);
}

double complex hz_hankel_ratio(int m, double x) {
  double j, y, dj, dy;
  hz_bessel_jy(m, x, &j, &y, &dj, &dy);
  double complex h = CMPLX(j, y), dh = CMPLX(dj, dy);
  /* Y_m overflows for m >> x (around m ~ 3.6x at x = 200). There H^(1)_m is
   * pure evanescent and the ratio tends to -m/x; the test says where this
   * branch starts, so it is never entered silently inside the validated range. */
  if (!isfinite(y) || !isfinite(dy) || !(cabs(h) > 0.0)) return -(double)m / x;
  return dh / h;
}

double hz_bessel_wronskian_err(int m, double x) {
  double j, y, dj, dy;
  hz_bessel_jy(m, x, &j, &y, &dj, &dy);
  double want = 2.0 / (M_PI * x);
  double got = j * dy - dj * y;
  if (!isfinite(got)) return INFINITY;
  return fabs(got - want) / fabs(want);
}
