#include "dtn2d.h"
#include "bessel.h"
#include "phi.h"
#include <math.h>

double complex hz_dtn_symbol(const hz_dtn *d, int m) {
  return d->k0 * hz_hankel_ratio(m < 0 ? -m : m, d->k0 * d->R);
}

/* Value and radial derivative of one basis function at a point on the circle.
 * B = phi(x/W - nx) phi(y/W - ny) exp(i K.(x - c)), so
 *   dB/dx = [phi'(.)/W + i kx phi(.)] phi(.) exp(...),
 * and d/dr = (x dx + y dy)/r. The phase is referenced to the ELEMENT, which is
 * not cosmetic: with a complex carrier a global reference is the difference
 * between a usable basis and an unusable one (TERMINATION_REPORT). */
static void bval(hz_carrier2d b, double x, double y, double r, double complex *v,
                 double complex *dv) {
  double tx = x / b.W - (double)b.nx, ty = y / b.W - (double)b.ny;
  if (fabs(tx) >= 2.0 || fabs(ty) >= 2.0) {
    *v = 0.0;
    *dv = 0.0;
    return;
  }
  double px = hz_phi(tx), py = hz_phi(ty);
  double dpx = hz_phi_d1(tx) / b.W, dpy = hz_phi_d1(ty) / b.W;
  double complex i = CMPLX(0.0, 1.0);
  double complex e = cexp(i * (b.kx * (x - b.W * (double)b.nx) + b.ky * (y - b.W * (double)b.ny)));
  *v = px * py * e;
  double complex bx = (dpx + i * b.kx * px) * py * e;
  double complex by = px * (dpy + i * b.ky * py) * e;
  *dv = (x * bx + y * by) / r;
}

void hz_dtn_harm(const hz_dtn *d, int m, hz_carrier2d b, double complex *uho,
                 double complex *duho) {
  double complex uh = 0.0, duh = 0.0;
  double complex i = CMPLX(0.0, 1.0);
  for (int p = 0; p < d->ntheta; p++) {
    double th = 2.0 * M_PI * (double)p / (double)d->ntheta;
    double x = d->R * cos(th), y = d->R * sin(th);
    double complex v, dv;
    bval(b, x, y, d->R, &v, &dv);
    if (!(cabs(v) > 0.0) && !(cabs(dv) > 0.0)) continue;
    double complex w = cexp(-i * (double)m * th);
    uh += v * w;
    duh += dv * w;
  }
  *uho = uh / (double)d->ntheta;
  *duho = duh / (double)d->ntheta;
}

double complex hz_dtn_row(const hz_dtn *d, int m, hz_carrier2d b) {
  double complex uh, duh;
  hz_dtn_harm(d, m, b, &uh, &duh);
  return duh - hz_dtn_symbol(d, m) * uh;
}
