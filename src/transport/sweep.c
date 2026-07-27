/* sweep.c — T1: ядро развёртки, чистое поглощение. См. sweep.h. */
#include "sweep.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

/* Ньютон по P_n сходится квадратично от чебышёвского начального приближения:
 * из 1e-1 до 1e-15 это ~5 шагов. TR_GL_TOL = 1e-15 — примерно 4 ulp на [−1,1],
 * практический пол двойной точности для корня; ниже итерация зациклится на
 * округлении. TR_GL_ITER — предохранитель, а не рабочая величина. */
#define TR_GL_TOL 1e-15
#define TR_GL_ITER 100

/* Узлы и веса Гаусса–Лежандра на [−1,1], n узлов, по возрастанию. */
static void gauss_legendre(int n, double *x, double *w) {
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
    /* x[i] — отрицательный конец, отражение — положительный */
    x[i] = -z;
    x[n - 1 - i] = z;
    w[i] = 2.0 / ((1.0 - z * z) * pp * pp);
    w[n - 1 - i] = w[i];
  }
}

int tr_ordinates_double_gauss(tr_ordinates *o, int nhalf) {
  o->n = 0;
  o->mu = NULL;
  o->w = NULL;
  if (nhalf < 1) return -1;

  double *gx = calloc((size_t)nhalf, sizeof *gx);
  double *gw = calloc((size_t)nhalf, sizeof *gw);
  if (!gx || !gw) {
    free(gx);
    free(gw);
    return -1;
  }
  gauss_legendre(nhalf, gx, gw);

  int n = 2 * nhalf;
  double *mu = calloc((size_t)n, sizeof *mu);
  double *w = calloc((size_t)n, sizeof *w);
  if (!mu || !w) {
    free(gx);
    free(gw);
    free(mu);
    free(w);
    return -1;
  }
  /* [−1,1] → [0,1]: μ = (1+g)/2, вес g/2 ⇒ Σ по полусфере = 1, по сфере = 2 */
  for (int i = 0; i < nhalf; i++) {
    double m = 0.5 * (1.0 + gx[i]);
    double ww = 0.5 * gw[i];
    mu[nhalf + i] = m;
    w[nhalf + i] = ww;
    mu[nhalf - 1 - i] = -m;
    w[nhalf - 1 - i] = ww;
  }
  free(gx);
  free(gw);

  o->n = n;
  o->mu = mu;
  o->w = w;
  return 0;
}

void tr_ordinates_free(tr_ordinates *o) {
  free(o->mu);
  free(o->w);
  o->mu = NULL;
  o->w = NULL;
  o->n = 0;
}

int tr_slab_alloc(tr_slab *s, int ncell) {
  s->ncell = 0;
  s->z = NULL;
  s->sig_lo = NULL;
  s->sig_hi = NULL;
  if (ncell < 1) return -1;

  double *z = calloc((size_t)ncell + 1, sizeof *z);
  double *lo = calloc((size_t)ncell, sizeof *lo);
  double *hi = calloc((size_t)ncell, sizeof *hi);
  if (!z || !lo || !hi) {
    free(z);
    free(lo);
    free(hi);
    return -1;
  }
  s->ncell = ncell;
  s->z = z;
  s->sig_lo = lo;
  s->sig_hi = hi;
  return 0;
}

void tr_slab_free(tr_slab *s) {
  free(s->z);
  free(s->sig_lo);
  free(s->sig_hi);
  s->z = NULL;
  s->sig_lo = NULL;
  s->sig_hi = NULL;
  s->ncell = 0;
}

void tr_slab_grade(tr_slab *s, double zmax, double stretch) {
  int n = s->ncell;
  assert(stretch > 0.0);
  /* Δz_i ∝ q^i с q = stretch^{−1/(n−1)}: отношение крайних ячеек равно stretch
   * при ЛЮБОМ n, поэтому измельчение уменьшает каждую ячейку ∝ 1/n (см. sweep.h,
   * почему шаг-на-ячейку здесь не годится). Сумма прогрессии нормируется на
   * zmax, так что z[n] = zmax до округления. */
  double q = (n > 1) ? pow(stretch, -1.0 / (double)(n - 1)) : 1.0;
  double step = 1.0, total = 0.0;
  for (int i = 0; i < n; i++) {
    total += step;
    step *= q;
  }
  step = 1.0;
  double acc = 0.0;
  s->z[0] = 0.0;
  for (int i = 0; i < n; i++) {
    acc += step;
    step *= q;
    s->z[i + 1] = zmax * acc / total;
  }
  s->z[n] = zmax; /* убрать дрейф накопления на последней грани */
}

void tr_slab_set_linear(tr_slab *s, double a, double b) {
  for (int i = 0; i < s->ncell; i++) {
    s->sig_lo[i] = a + b * s->z[i];
    s->sig_hi[i] = a + b * s->z[i + 1];
  }
}

double tr_cell_absorb(double sig_in, double sig_out, double len, double l_in) {
  /* Трапеция для ЛИНЕЙНОЙ σ_t точна ⇒ это замкнутая форма, не квадратура. */
  double tau = 0.5 * (sig_in + sig_out) * len;
  return l_in * exp(-tau);
}

double tr_sweep_absorb(const tr_slab *s, double mu, double g, double l_in, double *l_face) {
  double amu = fabs(mu);
  assert(amu > 0.0); /* μ = 0 в двойном Гауссе не бывает; хорда не определена */
  double l = l_in;
  if (mu > 0.0) {
    if (l_face) l_face[0] = l;
    for (int i = 0; i < s->ncell; i++) {
      double len = (s->z[i + 1] - s->z[i]) / amu;
      l = tr_cell_absorb(g * s->sig_lo[i], g * s->sig_hi[i], len, l);
      if (l_face) l_face[i + 1] = l;
    }
  } else {
    if (l_face) l_face[s->ncell] = l;
    for (int i = s->ncell - 1; i >= 0; i--) {
      double len = (s->z[i + 1] - s->z[i]) / amu;
      /* вток сверху: входная грань ячейки — верхняя */
      l = tr_cell_absorb(g * s->sig_hi[i], g * s->sig_lo[i], len, l);
      if (l_face) l_face[i] = l;
    }
  }
  return l;
}

double tr_moment0(const tr_ordinates *o, const double *l) {
  double sum = 0.0;
  for (int m = 0; m < o->n; m++)
    sum += o->w[m] * l[m];
  return 2.0 * M_PI * sum;
}

double tr_moment1(const tr_ordinates *o, const double *l) {
  double sum = 0.0;
  for (int m = 0; m < o->n; m++)
    sum += o->w[m] * o->mu[m] * l[m];
  return 2.0 * M_PI * sum;
}
