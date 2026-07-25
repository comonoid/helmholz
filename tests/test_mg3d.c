/* diagnostic: matvec vs dense; restriction vs direct rhs */
#include "solver3d.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static const double K0 = 2.0 * M_PI / 16.0;
static const double ALPHA = 0.05;

int main(void) {
  int DOM = 16;
  hz_octree t;
  hz_oct_init(&t, 4, CMPLX(K0 * K0, K0 * K0 * ALPHA));
  /* a small block so V has medium content too */
  int b0[3] = {6, 6, 6}, b1[3] = {10, 9, 8};
  hz_oct_set_box(&t, b0, b1, CMPLX(2.0 * K0 * K0, 0.1));
  hz_src3 src = {{1, {4, 4, 4}}, CMPLX(1.0, 0.0)};
  hz_scene3 sc = {&t, {DOM, DOM, DOM}, CMPLX(K0 * K0, K0 * K0 * ALPHA), &src, 1};

  int m = DOM / 2 + 3, dim = m * m * m; /* lvl 1 */
  int dom[3] = {DOM, DOM, DOM};
  double complex *x = calloc((size_t)dim, sizeof(double complex));
  double complex *y = calloc((size_t)dim, sizeof(double complex));
  double complex *yd = calloc((size_t)dim, sizeof(double complex));
  if (x == NULL || y == NULL || yd == NULL) {
    printf("test_mg3d: 0/2 passed (alloc)\n");
    return 1;
  }
  for (int i = 0; i < dim; i++)
    x[i] = CMPLX((double)((i * 37) % 11) - 5.0, (double)((i * 53) % 7) - 3.0);
  hz_dbg_matvec(&sc, 1, x, y);
  /* dense via exact entries */
  double num = 0.0, den = 0.0;
  for (int z = 0; z < m; z++)
    for (int yy = 0; yy < m; yy++)
      for (int xx = 0; xx < m; xx++) {
        int i = (z * m + yy) * m + xx;
        hz_pot3 a = {1, {xx - 1, yy - 1, z - 1}};
        double complex s = 0.0;
        for (int zz = 0; zz < m; zz++)
          for (int y2 = 0; y2 < m; y2++)
            for (int x2 = 0; x2 < m; x2++) {
              hz_pot3 b = {1, {x2 - 1, y2 - 1, zz - 1}};
              double complex e = hz_entry3d(&t, dom, a, b);
              if (cabs(e) > 0.0) s += e * x[(zz * m + y2) * m + x2];
            }
        yd[i] = s;
        num += creal((y[i] - s) * conj(y[i] - s));
        den += creal(s * conj(s));
      }
  double rel1 = sqrt(num / den);
  printf("matvec (T_fft + V) vs dense: rel = %.3e\n", rel1);

  /* restriction identity: <Phi_coarse, f> == restrict(<Phi_floor, f>) */
  int mc = DOM / 4 + 3;
  static const double W5[5] = {0.25, 0.5, 0.5, 0.5, 0.25};
  num = den = 0.0;
  for (int z = 0; z < mc; z++)
    for (int yy = 0; yy < mc; yy++)
      for (int xx = 0; xx < mc; xx++) {
        hz_pot3 c = {2, {xx - 1, yy - 1, z - 1}};
        double complex direct = hz_rhs3d(dom, c, src.p, src.amp);
        double complex viaR = 0.0;
        for (int sz = -2; sz <= 2; sz++)
          for (int sy = -2; sy <= 2; sy++)
            for (int sx = -2; sx <= 2; sx++) {
              int fx = 2 * (xx - 1) + sx, fy = 2 * (yy - 1) + sy, fz = 2 * (z - 1) + sz;
              if (fx < -1 || fy < -1 || fz < -1 || fx > DOM / 2 + 1 || fy > DOM / 2 + 1 ||
                  fz > DOM / 2 + 1)
                continue;
              hz_pot3 fpot = {1, {fx, fy, fz}};
              viaR += W5[sx + 2] * W5[sy + 2] * W5[sz + 2] * hz_rhs3d(dom, fpot, src.p, src.amp);
            }
        num += creal((direct - viaR) * conj(direct - viaR));
        den += creal(direct * conj(direct)) + 1e-30;
      }
  double rel2 = sqrt(num / den);
  printf("restriction identity: rel = %.3e\n", rel2);
  free(x);
  free(y);
  free(yd);
  hz_oct_free(&t);
  int npass = (rel1 < 1e-12 ? 1 : 0) + (rel2 < 1e-12 ? 1 : 0);
  printf("test_mg3d: %d/2 passed\n", npass);
  return npass == 2 ? 0 : 1;
}
