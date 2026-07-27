/* PLAN_TRANSPORT.md T3. Слой в оптической глубине, изотропное рассеяние,
 * источник первого столкновения, H-функция как эталон. */

#include "transport/scat1d.h"
#include "transport/quad.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* --- H-функция: эталон, вычисляемый независимо от развёртки ---------------- */

double tr_hfunc_moment0(double c) {
  return 2.0 / c * (1.0 - sqrt(1.0 - c));
}

int tr_hfunc(double c, int nquad, int maxit, double tol, const double *mu, int nmu, double *h) {
  if (nquad < 2 || nmu < 1) return 1;
  double *x = calloc((size_t)nquad, sizeof(double));
  double *w = calloc((size_t)nquad, sizeof(double));
  double *hq = calloc((size_t)nquad, sizeof(double));
  double *hn = calloc((size_t)nquad, sizeof(double));
  if (x == NULL || w == NULL || hq == NULL || hn == NULL) {
    free(x);
    free(w);
    free(hq);
    free(hn);
    return 1;
  }
  /* узлы Гаусса на [0,1]: из [−1,1] переносом */
  tr_gauss_legendre(nquad, x, w);
  for (int i = 0; i < nquad; i++) {
    x[i] = 0.5 * (x[i] + 1.0);
    w[i] *= 0.5;
    hq[i] = 1.0;
  }
  for (int it = 0; it < maxit; it++) {
    double dmax = 0.0;
    for (int i = 0; i < nquad; i++) {
      double s = 0.0;
      for (int j = 0; j < nquad; j++)
        s += w[j] * hq[j] / (x[i] + x[j]);
      /* H = 1/(1 − (c/2)μ∫H/(μ+μ')) — форма, у которой неподвижная точка
       * устойчива и при c → 1, в отличие от прямой подстановки */
      double v = 1.0 / (1.0 - 0.5 * c * x[i] * s);
      double d = fabs(v - hq[i]);
      if (d > dmax) dmax = d;
      hn[i] = v;
    }
    for (int i = 0; i < nquad; i++)
      hq[i] = hn[i];
    if (dmax < tol) break;
  }
  for (int k = 0; k < nmu; k++) {
    double s = 0.0;
    for (int j = 0; j < nquad; j++)
      s += w[j] * hq[j] / (mu[k] + x[j]);
    h[k] = 1.0 / (1.0 - 0.5 * c * mu[k] * s);
  }
  /* моментное тождество считается на СВОЕЙ квадратуре и возвращается вызывающему
   * через tr_hfunc_moment0 — здесь только сама H */
  free(x);
  free(w);
  free(hq);
  free(hn);
  return 0;
}

/* --- развёртка по слою ----------------------------------------------------- */

/* Максимум ячеек — страховка от вырожденных параметров, а не физический предел. */
#define TR_SCAT1D_MAXCELL 200000

int tr_scat1d_solve(const tr_scat1d *p, const tr_ordinates *o, int maxit, double tol,
                    tr_scat1d_out *out) {
  int nm = o->n;
  if (nm < 2 || !(p->mu0 > 0.0) || !(p->dtau0 > 0.0) || !(p->growth >= 1.0)) return 1;
  if (!(p->taumax > p->dtau0)) return 1;

  /* ГРАДУИРОВАННАЯ СЕТКА: шаг растёт от поверхности вглубь. Число ячеек — не
   * параметр, а следствие: свип по taumax при этом НЕ трогает разрешение у
   * поверхности, и К22 меряет одну вещь, а не две. */
  double *dt = calloc(TR_SCAT1D_MAXCELL, sizeof(double));
  if (dt == NULL) return 1;
  int nc = 0;
  double acc = 0.0, step = p->dtau0;
  while (acc < p->taumax && nc < TR_SCAT1D_MAXCELL) {
    double s = step;
    if (acc + s > p->taumax) s = p->taumax - acc;
    dt[nc++] = s;
    acc += s;
    step *= p->growth;
  }
  if (nc < 2) {
    free(dt);
    return 1;
  }

  double *S = calloc((size_t)nc, sizeof(double));  /* источник по ячейкам */
  double *J = calloc((size_t)nc, sizeof(double));  /* средняя интенсивность */
  double *Q1 = calloc((size_t)nc, sizeof(double)); /* источник первого столкновения */
  double *Jn = calloc((size_t)nc, sizeof(double));
  if (S == NULL || J == NULL || Q1 == NULL || Jn == NULL) {
    free(dt);
    free(S);
    free(J);
    free(Q1);
    free(Jn);
    return 1;
  }

  /* К23: незатенённая часть — ТОЧНО, и в правую часть идёт только однократно
   * рассеянная компонента. Среднее по ячейке берётся ТОЧНЫМ интегралом
   * экспоненты, а не значением в центре: иначе при dtau ~ mu0 появится ошибка,
   * которую потом припишут угловой сетке. */
  double tacc = 0.0;
  for (int i = 0; i < nc; i++) {
    double ta = tacc, tb = ta + dt[i];
    tacc = tb;
    double mean = p->mu0 / dt[i] * (exp(-ta / p->mu0) - exp(-tb / p->mu0));
    Q1[i] = p->beam_as_bc ? 0.0 : 0.25 * p->c * p->f * mean;
  }

  double *iface = calloc((size_t)(nc + 1), sizeof(double));
  if (iface == NULL) {
    free(dt);
    free(S);
    free(J);
    free(Q1);
    free(Jn);
    return 1;
  }
  memset(out->iout, 0, (size_t)nm * sizeof(double));

  int it = 0;
  double resid = 0.0;
  for (it = 0; it < maxit; it++) {
    for (int i = 0; i < nc; i++)
      S[i] = p->pscale * p->c * J[i] + Q1[i];
    for (int i = 0; i < nc; i++)
      Jn[i] = 0.0;

    for (int m = 0; m < nm; m++) {
      double mu = o->mu[m], amu = fabs(mu);
      double l;
      if (mu < 0.0) { /* вниз: с верхней границы */
        l = 0.0;
        /* НЕГАТИВНЫЙ КОНТРОЛЬ К23: пучок влётом по ближайшей ординате */
        if (p->beam_as_bc) {
          int best = -1;
          double bd = 1e300;
          for (int k = 0; k < nm; k++) {
            if (o->mu[k] >= 0.0) continue;
            double d = fabs(-o->mu[k] - p->mu0);
            if (d < bd) {
              bd = d;
              best = k;
            }
          }
          /* тот же поток, что у пучка: 2π·w·|μ|·I = π·f·μ0 */
          if (m == best) l = p->f * p->mu0 / (2.0 * o->w[m] * amu);
        }
        iface[0] = l;
        for (int i = 0; i < nc; i++) {
          double ex = exp(-dt[i] / amu);
          l = l * ex + S[i] * (1.0 - ex);
          iface[i + 1] = l;
        }
      } else { /* вверх: снизу, полубесконечность = нет влёта */
        l = 0.0;
        iface[nc] = l;
        for (int i = nc - 1; i >= 0; i--) {
          double ex = exp(-dt[i] / amu);
          l = l * ex + S[i] * (1.0 - ex);
          iface[i] = l;
        }
        out->iout[m] = l; /* радианс на ВЕРХНЕЙ границе */
      }
      /* среднее по ячейке для J: точное среднее характеристики
       * L(s) = L_вх e^{−s} + S(1−e^{−s}) по ячейке */
      for (int i = 0; i < nc; i++) {
        double tcell = dt[i] / amu, ex = exp(-tcell);
        double lin = mu < 0.0 ? iface[i] : iface[i + 1];
        double avg = S[i] + (lin - S[i]) * (1.0 - ex) / tcell;
        Jn[i] += 0.5 * o->w[m] * avg;
      }
    }

    resid = 0.0;
    for (int i = 0; i < nc; i++) {
      double d = fabs(Jn[i] - J[i]);
      if (d > resid) resid = d;
      J[i] = Jn[i];
    }
    if (resid < tol) break;
  }

  /* плоское альбедо: отражённый поток на падающий. Падающий на ГОРИЗОНТАЛЬНУЮ
   * площадку есть F·μ0 в той же нормировке, что и I. */
  double fup = 0.0;
  for (int m = 0; m < nm; m++)
    if (o->mu[m] > 0.0) fup += o->w[m] * o->mu[m] * out->iout[m];
  out->albedo = 2.0 * fup / (p->f * p->mu0);
  out->iters = it;
  out->resid = resid;

  free(dt);
  free(S);
  free(J);
  free(Q1);
  free(Jn);
  free(iface);
  return 0;
}
