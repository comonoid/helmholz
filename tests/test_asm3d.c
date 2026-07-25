/* M5 acceptance: hz_entry3d vs independent 3D Gauss-Legendre quadrature. */
#include "assemble3d.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int g_fail = 0;
static int g_total = 0;

static void check_close(double complex got, double complex want, double tol, const char *what) {
  g_total++;
  if (!(cabs(got - want) <= tol)) {
    g_fail++;
    printf("FAIL: %s: got %.12g%+.12gi want %.12g%+.12gi\n", what, creal(got), cimag(got),
           creal(want), cimag(want));
  }
}

static const double GLX[5] = {0.0, -0.5384693101056831, 0.5384693101056831, -0.9061798459386640,
                              0.9061798459386640};
static const double GLW[5] = {0.5688888888888889, 0.4786286704993665, 0.4786286704993665,
                              0.2369268850561891, 0.2369268850561891};

static double potv(hz_pot3 p, int ax, double x, int deriv) {
  double h = (double)(1 << p.lvl);
  double t = x / h - (double)p.n[ax];
  if (deriv == 2) return hz_phi_d2(t) / (h * h);
  return hz_phi(t);
}

/* quadrature over unit cells: exact per knot-free cell (knots on integers) */
static double complex brute_entry(const hz_octree *t, const int dom[3], hz_pot3 a, hz_pot3 b) {
  int lo[3], hi[3];
  for (int ax = 0; ax < 3; ax++) {
    int ha = 1 << a.lvl, hb = 1 << b.lvl;
    lo[ax] = ha * (a.n[ax] - 2) > hb * (b.n[ax] - 2) ? ha * (a.n[ax] - 2) : hb * (b.n[ax] - 2);
    hi[ax] = ha * (a.n[ax] + 2) < hb * (b.n[ax] + 2) ? ha * (a.n[ax] + 2) : hb * (b.n[ax] + 2);
    if (lo[ax] < 0) lo[ax] = 0;
    if (hi[ax] > dom[ax]) hi[ax] = dom[ax];
    if (lo[ax] >= hi[ax]) return 0.0;
  }
  double complex total = 0.0;
  for (int cx = lo[0]; cx < hi[0]; cx++)
    for (int cy = lo[1]; cy < hi[1]; cy++)
      for (int cz = lo[2]; cz < hi[2]; cz++) {
        double complex k2 = hz_oct_at(t, cx, cy, cz);
        double complex cell = 0.0;
        for (int i = 0; i < 5; i++)
          for (int j = 0; j < 5; j++)
            for (int k = 0; k < 5; k++) {
              double x = (double)cx + 0.5 + 0.5 * GLX[i];
              double y = (double)cy + 0.5 + 0.5 * GLX[j];
              double z = (double)cz + 0.5 + 0.5 * GLX[k];
              double A = potv(a, 0, x, 0) * potv(a, 1, y, 0) * potv(a, 2, z, 0);
              double lapB = potv(b, 0, x, 2) * potv(b, 1, y, 0) * potv(b, 2, z, 0) +
                            potv(b, 0, x, 0) * potv(b, 1, y, 2) * potv(b, 2, z, 0) +
                            potv(b, 0, x, 0) * potv(b, 1, y, 0) * potv(b, 2, z, 2);
              double B0 = potv(b, 0, x, 0) * potv(b, 1, y, 0) * potv(b, 2, z, 0);
              cell += GLW[i] * GLW[j] * GLW[k] * A * (lapB + k2 * B0);
            }
        total += cell * 0.125; /* (0.5)^3 jacobian */
      }
  return total;
}

int main(void) {
  hz_octree t;
  if (hz_oct_init(&t, 4, CMPLX(0.154, 0.003)) != 0) return 1; /* 16^3 */
  int w0[3] = {0, 6, 0}, w1[3] = {16, 8, 16};
  hz_oct_set_box(&t, w0, w1, CMPLX(1.39, 0.07));
  double c[3] = {10.0, 11.0, 6.0};
  hz_oct_set_ball(&t, c, 3.2, CMPLX(0.62, 0.01));
  int dom[3] = {16, 16, 16};

  /* sym=0 for boundary-clipped pairs: the strong form <Phi_a, L Phi_b> keeps
   * a boundary term after clipping, so exact symmetry holds only in the
   * interior (the clipped shell lies inside the absorbing region anyway) */
  static const struct {
    hz_pot3 a, b;
    int sym;
    const char *name;
  } cases[] = {
      {{0, {5, 6, 5}}, {0, {5, 6, 5}}, 1, "h1 self, on wall"},
      {{0, {5, 6, 5}}, {0, {6, 7, 4}}, 1, "h1 diag neighbors"},
      {{1, {3, 3, 3}}, {0, {7, 6, 6}}, 1, "h2 x h1 cross-level"},
      {{2, {2, 2, 1}}, {0, {9, 10, 6}}, 0, "h4 x h1 across ball"},
      {{1, {0, 3, 3}}, {1, {1, 3, 3}}, 0, "h2 pair at boundary clip"},
      {{2, {1, 1, 1}}, {2, {2, 2, 2}}, 0, "h4 pair spanning wall"},
      {{0, {14, 14, 14}}, {0, {15, 15, 15}}, 0, "h1 corner clip"},
  };
  int nc = (int)(sizeof(cases) / sizeof(cases[0]));
  for (int i = 0; i < nc; i++) {
    double complex got = hz_entry3d(&t, dom, cases[i].a, cases[i].b);
    double complex want = brute_entry(&t, dom, cases[i].a, cases[i].b);
    double scale = cabs(want) > 1.0 ? cabs(want) : 1.0;
    check_close(got, want, 1e-11 * scale, cases[i].name);
    /* the reversed order must also match ITS quadrature */
    double complex got2 = hz_entry3d(&t, dom, cases[i].b, cases[i].a);
    double complex want2 = brute_entry(&t, dom, cases[i].b, cases[i].a);
    check_close(got2, want2, 1e-11 * scale, "reversed vs quadrature");
    if (cases[i].sym) check_close(got2, got, 1e-11 * scale, "interior symmetry");
  }

  /* rhs: blob source = potential shape, quadrature-free closed form check */
  hz_pot3 aa = {1, {3, 3, 3}}, ss = {1, {3, 3, 3}};
  double complex r = hz_rhs3d(dom, aa, ss, CMPLX(2.0, 0.0));
  /* <phi,phi>^3 per axis at h=2: (2 * 92/15)^3, amp 2 */
  double m1 = 2.0 * 92.0 / 15.0;
  check_close(r, 2.0 * m1 * m1 * m1, 1e-10 * m1 * m1 * m1, "rhs blob self");

  long hits = 0, misses = 0;
  hz_asm3d_memo_stats(&hits, &misses);
  printf("memo: %ld hits / %ld misses\n", hits, misses);
  g_total++;
  if (misses > hits) {
    g_fail++;
    printf("FAIL: memo must actually help\n");
  }

  hz_oct_free(&t);
  printf("test_asm3d: %d/%d passed\n", g_total - g_fail, g_total);
  return g_fail > 0 ? 1 : 0;
}
