/* PLAN_CUT.md Р-5б. Квадратичная форма эрмитовых образцов: накопление, перенос,
 * сложение, минимум с усечением вырожденных направлений. */

#include "cut/qef.h"
#include <float.h>
#include <math.h>
#include <string.h>

void hz_qef_zero(hz_qef *q) {
  /* memset, а не по полям, И ЭТО НЕ ЛЕНЬ: формы копируются присваиванием
   * структуры и сравниваются ПОБИТОВО (Г49, починка против перестройки). Байты
   * выравнивания при заполнении по полям остались бы неопределёнными — ровно
   * дыра (b) из Р-2, где сырой дамп узла писал в файл набивку. Нули дают +0.0 и
   * целый ноль, то есть то же, что и по полям. */
  memset(q, 0, sizeof *q);
}

void hz_qef_add_sample(hz_qef *q, const double p[3], const double nrm[3]) {
  double d = nrm[0] * p[0] + nrm[1] * p[1] + nrm[2] * p[2];
  q->a[0] += nrm[0] * nrm[0];
  q->a[1] += nrm[0] * nrm[1];
  q->a[2] += nrm[0] * nrm[2];
  q->a[3] += nrm[1] * nrm[1];
  q->a[4] += nrm[1] * nrm[2];
  q->a[5] += nrm[2] * nrm[2];
  for (int k = 0; k < 3; k++) {
    q->b[k] += nrm[k] * d;
    q->m[k] += p[k];
  }
  q->c += d * d;
  q->n++;
}

/* y = A x, A из шести чисел */
static void amul(const double a[6], const double x[3], double y[3]) {
  y[0] = a[0] * x[0] + a[1] * x[1] + a[2] * x[2];
  y[1] = a[1] * x[0] + a[3] * x[1] + a[4] * x[2];
  y[2] = a[2] * x[0] + a[4] * x[1] + a[5] * x[2];
}

void hz_qef_add_shifted(hz_qef *dst, const hz_qef *src, const double t[3]) {
  double at[3];
  amul(src->a, t, at);
  double tb = t[0] * src->b[0] + t[1] * src->b[1] + t[2] * src->b[2];
  double tat = t[0] * at[0] + t[1] * at[1] + t[2] * at[2];
  for (int i = 0; i < 6; i++)
    dst->a[i] += src->a[i];
  for (int k = 0; k < 3; k++) {
    dst->b[k] += src->b[k] + at[k];
    dst->m[k] += src->m[k] + t[k] * (double)src->n;
  }
  dst->c += src->c + 2.0 * tb + tat;
  dst->n += src->n;
}

double hz_qef_eval(const hz_qef *q, const double x[3]) {
  double ax[3];
  amul(q->a, x, ax);
  double xax = x[0] * ax[0] + x[1] * ax[1] + x[2] * ax[2];
  double bx = q->b[0] * x[0] + q->b[1] * x[1] + q->b[2] * x[2];
  return xax - 2.0 * bx + q->c;
}

void hz_qef_jacobi3(const double a[6], double val[3], double vec[3][3]) {
  double m[3][3] = {{a[0], a[1], a[2]}, {a[1], a[3], a[4]}, {a[2], a[4], a[5]}};
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      vec[i][j] = i == j ? 1.0 : 0.0;

  double nrm2 = 0.0;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      nrm2 += m[i][j] * m[i][j];
  /* выход по норме внедиагонали относительно нормы матрицы: дальше вращения
   * ничего не меняют в двойной точности, а HZ_QEF_SWEEPS — только страховка */
  double stop = DBL_EPSILON * DBL_EPSILON * nrm2;

  static const int P[3] = {0, 0, 1}, Q[3] = {1, 2, 2};
  for (int sweep = 0; sweep < HZ_QEF_SWEEPS; sweep++) {
    double off = m[0][1] * m[0][1] + m[0][2] * m[0][2] + m[1][2] * m[1][2];
    if (!(off > stop)) break;
    for (int k = 0; k < 3; k++) {
      int p = P[k], q = Q[k];
      double apq = m[p][q];
      if (!(fabs(apq) > 0.0)) continue;
      double theta = (m[q][q] - m[p][p]) / (2.0 * apq);
      double t = (theta >= 0.0 ? 1.0 : -1.0) / (fabs(theta) + sqrt(theta * theta + 1.0));
      double c = 1.0 / sqrt(t * t + 1.0), s = t * c;
      double app = m[p][p], aqq = m[q][q];
      m[p][p] = app - t * apq;
      m[q][q] = aqq + t * apq;
      m[p][q] = m[q][p] = 0.0;
      for (int r = 0; r < 3; r++) {
        if (r == p || r == q) continue;
        double arp = m[r][p], arq = m[r][q];
        m[r][p] = m[p][r] = c * arp - s * arq;
        m[r][q] = m[q][r] = c * arq + s * arp;
      }
      for (int r = 0; r < 3; r++) {
        double vrp = vec[r][p], vrq = vec[r][q];
        vec[r][p] = c * vrp - s * vrq;
        vec[r][q] = s * vrp + c * vrq;
      }
    }
  }
  for (int i = 0; i < 3; i++)
    val[i] = m[i][i];
}

int hz_qef_solve_tau(const hz_qef *q, double tau_rel, double x[3], double *resid) {
  if (q->n <= 0) {
    for (int k = 0; k < 3; k++)
      x[k] = 0.0;
    *resid = 0.0;
    return HZ_QEF_EMPTY;
  }
  double xb[3];
  for (int k = 0; k < 3; k++)
    xb[k] = q->m[k] / (double)q->n;

  double axb[3], r[3];
  amul(q->a, xb, axb);
  for (int k = 0; k < 3; k++)
    r[k] = q->b[k] - axb[k];

  double val[3], vec[3][3];
  hz_qef_jacobi3(q->a, val, vec);
  double lmax = 0.0;
  for (int i = 0; i < 3; i++)
    if (val[i] > lmax) lmax = val[i];
  double tau = tau_rel * lmax;

  for (int k = 0; k < 3; k++)
    x[k] = xb[k];
  for (int i = 0; i < 3; i++) {
    if (!(val[i] > tau)) continue; /* направление, о котором данные молчат */
    double d = (vec[0][i] * r[0] + vec[1][i] * r[1] + vec[2][i] * r[2]) / val[i];
    for (int k = 0; k < 3; k++)
      x[k] += d * vec[k][i];
  }
  double e = hz_qef_eval(q, x);
  /* Невязка неотрицательна по построению; отрицательное значение здесь —
   * сокращение близких чисел у идеально подогнанной формы, а не ошибка. */
  *resid = e > 0.0 ? e : 0.0;
  return HZ_QEF_OK;
}

int hz_qef_solve(const hz_qef *q, const double lo[3], const double hi[3], double x[3],
                 double *resid) {
  int rc = hz_qef_solve_tau(q, HZ_QEF_TAU_REL, x, resid);
  if (rc != HZ_QEF_OK) return rc;
  int clamped = 0;
  for (int k = 0; k < 3; k++) {
    if (x[k] < lo[k]) {
      x[k] = lo[k];
      clamped = 1;
    } else if (x[k] > hi[k]) {
      x[k] = hi[k];
      clamped = 1;
    }
  }
  if (clamped) {
    double e = hz_qef_eval(q, x);
    *resid = e > 0.0 ? e : 0.0;
    return HZ_QEF_CLAMPED;
  }
  return HZ_QEF_OK;
}
