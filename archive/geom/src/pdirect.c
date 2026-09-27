/* Прямой свет от компактных источников. Разбор и оговорки — в `pdirect.h`. */

#include "pdirect.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PD_MAXV 16 /* угол источника + прибавки от отсечения горизонтом */

double hz_direct_unocc(const hz_polyset *ps, int32_t src, double Le, const double x[3],
                       const double n[3]) {
  const hz_poly *S = &ps->p[src];
  if (S->nloop == 0 || !(Le > 0.0)) return 0.0;
  /* ИСТОЧНИК ОДНОСТОРОННИЙ: излучает только со стороны `+n`. Приёмник с
   * изнанки не получает ничего, и это проверяется ЗДЕСЬ, а не знаком суммы —
   * знак суммы задаётся обходом петли, а не физикой. */
  double side = 0.0;
  for (int a = 0; a < 3; a++)
    side += S->n[a] * (x[a] - S->org[a]);
  if (!(side > 0.0)) return 0.0;
  int32_t b = ps->loop[S->l0], e = ps->loop[S->l0 + 1];
  int32_t nv = e - b;
  if (nv < 3 || nv > PD_MAXV - 2) return 0.0;

  /* вершины источника в мире, относительно x */
  double p[PD_MAXV][3], q[PD_MAXV][3];
  for (int32_t i = 0; i < nv; i++) {
    double w[3];
    hz_poly_world(S, ps->bv[(size_t)(b + i) * 2], ps->bv[(size_t)(b + i) * 2 + 1], w);
    for (int a = 0; a < 3; a++)
      p[i][a] = w[a] - x[a];
  }
  /* ОТСЕЧЕНИЕ ПО ГОРИЗОНТУ ПРИЁМНИКА: держим {n·p > 0}. */
  int32_t m = 0;
  for (int32_t i = 0; i < nv; i++) {
    const double *A = p[i], *B = p[(i + 1) % nv];
    double da = A[0] * n[0] + A[1] * n[1] + A[2] * n[2];
    double db = B[0] * n[0] + B[1] * n[1] + B[2] * n[2];
    if (da > 0.0) {
      if (m >= PD_MAXV) return 0.0;
      for (int a = 0; a < 3; a++)
        q[m][a] = A[a];
      m++;
    }
    if ((da > 0.0) != (db > 0.0)) {
      double s = da / (da - db);
      if (m >= PD_MAXV) return 0.0;
      for (int a = 0; a < 3; a++)
        q[m][a] = A[a] + s * (B[a] - A[a]);
      m++;
    }
  }
  if (m < 3) return 0.0;

  double sum = 0.0;
  for (int32_t i = 0; i < m; i++) {
    const double *A = q[i], *B = q[(i + 1) % m];
    double la = sqrt(A[0] * A[0] + A[1] * A[1] + A[2] * A[2]);
    double lb = sqrt(B[0] * B[0] + B[1] * B[1] + B[2] * B[2]);
    if (!(la > 0.0) || !(lb > 0.0)) continue;
    double c = (A[0] * B[0] + A[1] * B[1] + A[2] * B[2]) / (la * lb);
    if (c > 1.0) c = 1.0;
    if (c < -1.0) c = -1.0;
    double th = acos(c);
    double cr[3] = {A[1] * B[2] - A[2] * B[1], A[2] * B[0] - A[0] * B[2],
                    A[0] * B[1] - A[1] * B[0]};
    double cn = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
    if (!(cn > 0.0)) continue;
    sum += th * (cr[0] * n[0] + cr[1] * n[1] + cr[2] * n[2]) / cn;
  }
  /* МОДУЛЬ, А НЕ ЗАЖИМ НУЛЁМ. После отсечения по горизонту каждая оставшаяся
   * точка источника лежит выше горизонта приёмника, то есть вносит
   * ПОЛОЖИТЕЛЬНЫЙ вклад; знак суммы Ламберта при этом целиком задан НАПРАВЛЕНИЕМ
   * ОБХОДА петли источника, а оно — соглашение сборки полигона, не физика.
   * Первая редакция зажимала нулём и получала тождественный ноль по всей сцене:
   * петля светильника обходится так, что сумма выходит отрицательной. */
  double E = 0.5 * Le * sum;
  return (E > 0.0) ? E : -E;
}

/* Точка на источнике по параметрам (a, b) ∈ [0,1]² — по его uv-габариту.
 * Для прямоугольного светильника это точная равномерная выборка. */
static void src_point(const hz_polyset *ps, int32_t src, double a, double b, double x[3]) {
  const hz_poly *S = &ps->p[src];
  double u = S->uvlo[0] + a * (S->uvhi[0] - S->uvlo[0]);
  double v = S->uvlo[1] + b * (S->uvhi[1] - S->uvlo[1]);
  hz_poly_world(S, u, v, x);
}

int hz_direct_add(hz_ptrans *t, const hz_pray *g, const int32_t *src, int nsrc, const double *Le,
                  int nvis, double hsamp, int cap, hz_pstats *st) {
  const hz_polyset *ps = t->ps;
  if (nvis < 1) nvis = 1;
  if (cap < 4) cap = (cap <= 0) ? 64 : 4;
  int64_t nray = 0, nsample = 0;

#pragma omp parallel for schedule(dynamic, 8) reduction(+ : nray, nsample)
  for (int32_t k = 0; k < ps->np; k++) {
    const hz_poly *P = &ps->p[k];
    if (P->nloop == 0 || !(P->mom[0] > 0.0)) continue;
    int is_src = 0;
    for (int s = 0; s < nsrc; s++)
      if (src[s] == k) is_src = 1;
    if (is_src) continue; /* источник сам себя не освещает */

    /* Регулярная выборка по uv-габариту с отбором по краю. Вес каждой принятой
     * точки — не площадь ячейки, а `mom[0]/принято`: так суммарный вес РАВЕН
     * точной площади полигона, и ошибка отбора не течёт в энергию. */
    double du = P->uvhi[0] - P->uvlo[0], dv = P->uvhi[1] - P->uvlo[1];
    int nu = (int)(du / hsamp) + 1, nvv = (int)(dv / hsamp) + 1;
    if (nu < 4) nu = 4;
    if (nvv < 4) nvv = 4;
    /* ПОТОЛОК ВЫБОРКИ. Подгоняются ТРИ коэффициента гладкой функции, и `64×64`
     * точек для этого с огромным запасом. Прежний потолок `256` давал крупному
     * полигону 65 536 точек, то есть 18.9 млн теневых лучей на ОДИН полигон при
     * восьми источниках и 36 пробах, — час на конфигурацию и ни одной цифры
     * сверх того, что даёт `64`. */
    if (nu > cap) nu = cap;
    if (nvv > cap) nvv = cap;
    double s0 = 0.0, s1 = 0.0, s2 = 0.0;
    int32_t nin = 0;
    for (int iu = 0; iu < nu; iu++)
      for (int iv = 0; iv < nvv; iv++) {
        double u = P->uvlo[0] + ((double)iu + 0.5) / nu * du;
        double v = P->uvlo[1] + ((double)iv + 0.5) / nvv * dv;
        if (!hz_poly_inside(ps, P, u, v)) continue;
        nin++;
        double x[3];
        hz_poly_world(P, u, v, x);
        /* точку чуть поднять над плоскостью — иначе теневой луч ловит сам себя */
        double eps = 1e-7 * (1.0 + fabs(x[0]) + fabs(x[1]) + fabs(x[2]));
        double xo[3];
        for (int a = 0; a < 3; a++)
          xo[a] = x[a] + eps * P->n[a];
        double Etot = 0.0;
        for (int s = 0; s < nsrc; s++) {
          double Eu = hz_direct_unocc(ps, src[s], Le[src[s]], x, P->n);
          if (!(Eu > 0.0)) continue;
          int vis = 0;
          for (int a = 0; a < nvis; a++)
            for (int bq = 0; bq < nvis; bq++) {
              double sp[3];
              src_point(ps, src[s], ((double)a + 0.5) / nvis, ((double)bq + 0.5) / nvis, sp);
              double d[3];
              double len = 0.0;
              for (int c = 0; c < 3; c++) {
                d[c] = sp[c] - xo[c];
                len += d[c] * d[c];
              }
              len = sqrt(len);
              if (!(len > 0.0)) continue;
              for (int c = 0; c < 3; c++)
                d[c] /= len;
              nray++;
              if (!hz_pray_occluded(g, xo, d, 0.0, len * (1.0 - 1e-6))) vis++;
            }
          Etot += Eu * (double)vis / (double)(nvis * nvis);
        }
        s0 += Etot;
        s1 += Etot * u;
        s2 += Etot * v;
      }
    if (nin == 0) continue;
    nsample += nin;
    double sc = P->mom[0] / (double)nin;
    t->acc[(size_t)k * 3 + 0] += sc * s0;
    t->acc[(size_t)k * 3 + 1] += sc * s1;
    t->acc[(size_t)k * 3 + 2] += sc * s2;
  }
  if (st != NULL) {
    st->npair += nsample;
    st->nfrag += nray;
  }
  return 0;
}
