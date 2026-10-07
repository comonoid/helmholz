/* CBMC formal check for finding 1: sw_gauss_legendre (verbatim copy from
 * src/nstruct/sweep.c:78-100) fills EVERY slot for even n.  Two independent
 * buffers, same n: an unwritten slot stays nondet and differs -> assert fails. */
#include <assert.h>
#include <math.h>
#define HZ_SW_GL_EPS 1e-15
#define HZ_SW_GL_ITMAX 4 /* "узел Гаусса сходится за ~5" (sweep.c:17) */
static void sw_gauss_legendre(int n, double *x, double *w) {
  for (int i = 0; i < n / 2; i++) {
    double z = cos(M_PI * (i + 0.75) / (n + 0.5));
    double pp = 0.0;
    for (int it = 0; it < HZ_SW_GL_ITMAX; it++) {
      double p0 = 1.0, p1 = z;
      for (int k = 2; k <= n; k++) {
        double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / (double)k;
        p0 = p1;
        p1 = p2;
      }
      pp = (double)n * (z * p1 - p0) / (z * z - 1.0);
      double dz = p1 / pp;
      z -= dz;
      if (fabs(dz) <= HZ_SW_GL_EPS) break;
    }
    x[i] = -z;
    w[i] = 2.0 / ((1.0 - z * z) * pp * pp);
    x[n - 1 - i] = z;
    w[n - 1 - i] = w[i];
  }
}
int main(void) {
  int n;
  double x1[4], x2[4], w1[4], w2[4];
  __CPROVER_assume(n >= 2 && n <= 4 && (n % 2) == 0);
  sw_gauss_legendre(n, x1, w1);
  sw_gauss_legendre(n, x2, w2);
  for (int i = 0; i < n; i++) {
    assert(x1[i] == x2[i]);
    assert(w1[i] == w2[i]);
  }
  return 0;
}
