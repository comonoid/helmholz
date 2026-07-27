/* quad.c — узлы Гаусса–Лежандра ньютоном по P_n. См. quad.h. */
#include "quad.h"
#include <math.h>

/* Ньютон по P_n сходится квадратично от чебышёвского начального приближения:
 * из 1e-1 до 1e-15 это ~5 шагов. TR_GL_TOL = 1e-15 — примерно 4 ulp на [−1,1],
 * практический пол двойной точности для корня; ниже итерация зациклится на
 * округлении. TR_GL_ITER — предохранитель, а не рабочая величина. */
#define TR_GL_TOL 1e-15
#define TR_GL_ITER 100

void tr_gauss_legendre(int n, double *x, double *w) {
  int half = (n + 1) / 2; /* узлы симметричны, считаем половину */
  for (int i = 0; i < half; i++) {
    /* начальное приближение: нули P_n близки к косинусам (Chebyshev-like) */
    double z = cos(M_PI * ((double)i + 0.75) / ((double)n + 0.5));
    double pp = 0.0;
    for (int it = 0; it < TR_GL_ITER; it++) {
      double p1 = 1.0, p2 = 0.0;
      for (int j = 0; j < n; j++) { /* рекуррентность Лежандра */
        double p3 = p2;
        p2 = p1;
        p1 = ((2.0 * (double)j + 1.0) * z * p2 - (double)j * p3) / ((double)j + 1.0);
      }
      pp = (double)n * (z * p1 - p2) / (z * z - 1.0); /* P_n'(z) */
      double dz = p1 / pp;
      z -= dz;
      if (fabs(dz) < TR_GL_TOL) break;
    }
    x[i] = -z; /* x[i] — отрицательный конец, отражение — положительный */
    x[n - 1 - i] = z;
    w[i] = 2.0 / ((1.0 - z * z) * pp * pp);
    w[n - 1 - i] = w[i];
  }
}
