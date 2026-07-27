/* PLAN_TRANSPORT.md, шаг A. Ординаты на сфере: двойной Гаусс по μ (К10) на
 * произведение с гауссом по квадрантам азимута. */

#include "transport/dirs3.h"
#include "transport/quad.h"
#include <math.h>
#include <stdlib.h>

void tr3_dirs_free(tr3_dirs *d) {
  free(d->ox);
  free(d->oy);
  free(d->oz);
  free(d->w);
  d->ox = d->oy = d->oz = d->w = NULL;
  d->n = 0;
}

int tr3_dirs_product(tr3_dirs *d, int nmu, int nphi) {
  if (nmu < 1 || nphi < 1) return 1;
  int nm = 2 * nmu, np = 4 * nphi, n = nm * np;
  d->n = n;
  d->ox = calloc((size_t)n, sizeof(double));
  d->oy = calloc((size_t)n, sizeof(double));
  d->oz = calloc((size_t)n, sizeof(double));
  d->w = calloc((size_t)n, sizeof(double));
  double *gx = calloc((size_t)(nmu > nphi ? nmu : nphi), sizeof(double));
  double *gw = calloc((size_t)(nmu > nphi ? nmu : nphi), sizeof(double));
  if (d->ox == NULL || d->oy == NULL || d->oz == NULL || d->w == NULL || gx == NULL || gw == NULL) {
    tr3_dirs_free(d);
    free(gx);
    free(gw);
    return 1;
  }

  /* μ: ДВОЙНОЙ Гаусс — полусферы набираются РАЗДЕЛЬНО (К10) */
  double *mu = calloc((size_t)nm, sizeof(double));
  double *wmu = calloc((size_t)nm, sizeof(double));
  /* φ: Гаусс по КВАДРАНТАМ — узлы строго внутри, осей в наборе нет */
  double *ph = calloc((size_t)np, sizeof(double));
  double *wph = calloc((size_t)np, sizeof(double));
  if (mu == NULL || wmu == NULL || ph == NULL || wph == NULL) {
    tr3_dirs_free(d);
    free(gx);
    free(gw);
    free(mu);
    free(wmu);
    free(ph);
    free(wph);
    return 1;
  }

  tr_gauss_legendre(nmu, gx, gw);
  for (int h = 0; h < 2; h++) {
    double a = h ? 0.0 : -1.0, b = h ? 1.0 : 0.0;
    for (int i = 0; i < nmu; i++) {
      mu[h * nmu + i] = 0.5 * (b - a) * gx[i] + 0.5 * (a + b);
      wmu[h * nmu + i] = 0.5 * (b - a) * gw[i];
    }
  }
  tr_gauss_legendre(nphi, gx, gw);
  for (int q = 0; q < 4; q++) {
    double a = 0.5 * M_PI * (double)q, b = a + 0.5 * M_PI;
    for (int i = 0; i < nphi; i++) {
      ph[q * nphi + i] = 0.5 * (b - a) * gx[i] + 0.5 * (a + b);
      wph[q * nphi + i] = 0.5 * (b - a) * gw[i];
    }
  }

  int k = 0;
  for (int i = 0; i < nm; i++) {
    double st = sqrt(1.0 - mu[i] * mu[i]);
    for (int j = 0; j < np; j++) {
      d->ox[k] = st * cos(ph[j]);
      d->oy[k] = st * sin(ph[j]);
      d->oz[k] = mu[i];
      d->w[k] = wmu[i] * wph[j];
      k++;
    }
  }
  free(gx);
  free(gw);
  free(mu);
  free(wmu);
  free(ph);
  free(wph);
  return 0;
}
