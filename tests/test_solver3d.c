/* M6 acceptance: empty-scene transport vs the free-space Green's function.
 * u_ref(x) = -Int G(|x-y|) f(y) dy, G = e^{ikr}/(4 pi r), k = k0 sqrt(1+i a).
 * The solver sees the same medium plus a quartic absorbing shell; probes stay
 * in the clean interior. */
#include "solver3d.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int g_fail = 0;
static int g_total = 0;

static void check_lt(double got, double bound, const char *what) {
  g_total++;
  if (!(got < bound)) {
    g_fail++;
    printf("FAIL: %s: got %.6g bound %.6g\n", what, got, bound);
  }
}

enum { DOM = 48 };
static const double K0 = 2.0 * M_PI / 8.0; /* lambda = 8 cells */
static const double ALPHA = 0.12;

/* quartic absorbing shell, thickness SH cells, peak extra Im = PEAK*k0^2 */
static void add_shell(hz_octree *t, int dom, double k0) {
  const int SH = 8;
  const double PEAK = 1.2;
  double b = k0 * k0;
  for (int s = 0; s < SH; s++) {
    double tt = (double)(SH - s) / (double)SH;
    double a = ALPHA + (PEAK - ALPHA) * tt * tt * tt * tt;
    double complex k2 = CMPLX(b, b * a);
    int lo0[3] = {s, s, s}, hi0[3] = {dom - s, dom - s, dom - s};
    int lo1[3] = {s + 1, s + 1, s + 1}, hi1[3] = {dom - s - 1, dom - s - 1, dom - s - 1};
    /* one-cell thick shell layer: set the box, then the inner box will be
     * overwritten by the next (less absorbing) layer */
    (void)lo1;
    (void)hi1;
    hz_oct_set_box(t, lo0, hi0, k2);
  }
  int lo[3] = {SH, SH, SH}, hi[3] = {dom - SH, dom - SH, dom - SH};
  hz_oct_set_box(t, lo, hi, CMPLX(b, b * ALPHA));
}

/* reference: -G * f over the blob support, per-piece Gauss-Legendre 4 */
static double complex green_ref(const double p[3], const double s[3]) {
  static const double gx[4] = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563,
                               0.8611363115940526};
  static const double gw[4] = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461,
                               0.3478548451374538};
  double complex k = K0 * csqrt(CMPLX(1.0, ALPHA));
  double complex acc = 0.0;
  /* blob = phi(x/2 - n): support s +- 4, pieces of width 2 */
  for (int px = -2; px < 2; px++)
    for (int py = -2; py < 2; py++)
      for (int pz = -2; pz < 2; pz++)
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
            for (int l = 0; l < 4; l++) {
              double y0 = s[0] + 2.0 * px + 1.0 + gx[i];
              double y1 = s[1] + 2.0 * py + 1.0 + gx[j];
              double y2 = s[2] + 2.0 * pz + 1.0 + gx[l];
              double f =
                  hz_phi((y0 - s[0]) / 2.0) * hz_phi((y1 - s[1]) / 2.0) * hz_phi((y2 - s[2]) / 2.0);
              double dx = p[0] - y0, dy = p[1] - y1, dz = p[2] - y2;
              double r = sqrt(dx * dx + dy * dy + dz * dz);
              double complex G = cexp(CMPLX(0.0, 1.0) * k * r) / (4.0 * M_PI * r);
              acc += gw[i] * gw[j] * gw[l] * f * G;
            }
  return -acc; /* (Lap+k^2)G = -delta  =>  u = -G*f */
}

int main(void) {
  hz_octree t;
  if (hz_oct_init(&t, 6, CMPLX(K0 * K0, K0 * K0 * ALPHA)) != 0) return 1;
  add_shell(&t, DOM, K0);

  hz_src3 src = {{1, {12, 12, 12}}, CMPLX(1.0, 0.0)}; /* blob centered at 24 */
  hz_scene3 sc = {&t, {DOM, DOM, DOM}, CMPLX(K0 * K0, K0 * K0 * ALPHA), &src, 1};

  hz_sol3 sol;
  g_total++;
  if (hz_solve3d(&sc, 3, 0, 2, 1, &sol) != 0) {
    g_fail++;
    printf("FAIL: solve\n");
    printf("test_solver3d: %d/%d passed\n", g_total - g_fail, g_total);
    return 1;
  }

  double sctr[3] = {24.0, 24.0, 24.0};
  static const double probes[5][3] = {{34.0, 24.0, 24.0},
                                      {24.0, 33.0, 24.0},
                                      {24.0, 24.0, 34.0},
                                      {30.0, 30.0, 24.0},
                                      {30.0, 30.0, 30.0}};
  double num = 0.0, den = 0.0;
  for (int p = 0; p < 5; p++) {
    double complex got = hz_sol3_eval(&sol, probes[p]);
    double complex want = green_ref(probes[p], sctr);
    printf("probe %d: got %+.5e%+.5ei  want %+.5e%+.5ei  |rel| %.3f\n", p, creal(got), cimag(got),
           creal(want), cimag(want), cabs(got - want) / cabs(want));
    num += creal((got - want) * conj(got - want));
    den += creal(want * conj(want));
  }
  double rel = sqrt(num / den);
  printf("green check: rel L2 over probes = %.4f, floor relres = %.2e\n", rel, sol.final_relres);
  check_lt(rel, 0.2, "3D transport vs Green function");

  hz_sol3_free(&sol);
  hz_oct_free(&t);
  printf("test_solver3d: %d/%d passed\n", g_total - g_fail, g_total);
  return g_fail > 0 ? 1 : 0;
}
