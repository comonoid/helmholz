/* Подгонка вершин края в 3D. Разбор и оговорки — в `pvfit.h`. */

#include "pvfit.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Решение системы `n×n` (n ≤ 4) методом Гаусса с частичным выбором. 0 — успех. */
static int vf_solve(double A[4][4], double b[4], double x[4], int n) {
  /* Граница ЯВНАЯ: без неё статический анализ не видит, что `n ≤ 4`, и
   * справедливо ругается на выход за массив (гейт cppcheck). */
  if (n < 1 || n > 4) return 1;
  for (int c = 0; c < n; c++) {
    int piv = c;
    for (int r = c + 1; r < n; r++)
      if (fabs(A[r][c]) > fabs(A[piv][c])) piv = r;
    if (!(fabs(A[piv][c]) > 0.0)) return 1;
    if (piv != c) {
      for (int k = 0; k < n; k++) {
        double t = A[c][k];
        A[c][k] = A[piv][k];
        A[piv][k] = t;
      }
      double t = b[c];
      b[c] = b[piv];
      b[piv] = t;
    }
    for (int r = 0; r < n; r++) {
      if (r == c) continue;
      double f = A[r][c] / A[c][c];
      for (int k = 0; k < n; k++)
        A[r][k] -= f * A[c][k];
      b[r] -= f * b[c];
    }
  }
  for (int c = 0; c < n; c++)
    x[c] = b[c] / A[c][c];
  return 0;
}

/* Худшая невязка точки по набору плоскостей. */
static double vf_worst(const double *pl, int m, const double v[3]) {
  double w = 0.0;
  for (int i = 0; i < m; i++) {
    double d = fabs(pl[(size_t)i * 4] * v[0] + pl[(size_t)i * 4 + 1] * v[1] +
                    pl[(size_t)i * 4 + 2] * v[2] - pl[(size_t)i * 4 + 3]);
    if (d > w) w = d;
  }
  return w;
}

/* МИНИМАКС ПО ЧЕТВЁРКАМ ЗНАКОВЫХ ОГРАНИЧЕНИЙ — ТОЧНО, ПЕРЕБОРОМ.
 *
 * Задача: минимизировать `t` при `|n_i·v − d_i| ≤ t`. Это линейная программа с
 * четырьмя неизвестными `(v, t)`, и её оптимум достигается в вершине, то есть на
 * четырёх активных ограничениях. Плоскостей у вершины единицы (замер А122: у
 * города `47.6%` вершин имеют ровно три группы, `15.3%` — четыре, `8.2%` — пять),
 * поэтому перебор всех четвёрок ЗНАКОВЫХ ограничений дёшев и даёт ТОЧНЫЙ ответ:
 * при `m = 6` это `C(12,4) = 495` систем `4×4`.
 * Так честнее, чем итеративный решатель: у обменного алгоритма в `pmerge.c`
 * сходимость не доказана и замерена как `74%` недосходов (А89), а здесь ответ
 * точный по построению. */
static int vf_minimax(const double *pl, int m, double v[3], double *pt) {
  double best = 1e300, bv[3] = {0, 0, 0};
  int ok = 0;
  for (int i = 0; i < m; i++)
    for (int j = i + 1; j < m; j++)
      for (int k = j + 1; k < m; k++)
        for (int l = k + 1; l < m + 1; l++) {
          /* Четвёртое ограничение может отсутствовать (m == 3): тогда система
           * трёх плоскостей с `t = 0`, то есть точное пересечение. */
          int id[4] = {i, j, k, (l < m) ? l : -1};
          int nn = (l < m) ? 4 : 3;
          for (int sg = 0; sg < (1 << nn); sg++) {
            double A[4][4], b[4], x[4];
            memset(A, 0, sizeof A);
            for (int r = 0; r < nn; r++) {
              const double *P = pl + (size_t)id[r] * 4;
              double s = ((sg >> r) & 1) ? 1.0 : -1.0;
              A[r][0] = P[0];
              A[r][1] = P[1];
              A[r][2] = P[2];
              A[r][3] = -s;
              b[r] = P[3];
            }
            if (nn == 3) {
              /* `t` не участвует:три плоскости пересекаются в точке. */
              double A3[4][4], b3[4], x3[4];
              memcpy(A3, A, sizeof A3);
              memcpy(b3, b, sizeof b3);
              if (vf_solve(A3, b3, x3, 3) != 0) continue;
              double vv[3] = {x3[0], x3[1], x3[2]};
              double t = vf_worst(pl, m, vv);
              if (t < best) {
                best = t;
                bv[0] = vv[0];
                bv[1] = vv[1];
                bv[2] = vv[2];
                ok = 1;
              }
              break; /* знаки при `t = 0` не важны */
            }
            if (vf_solve(A, b, x, 4) != 0) continue;
            if (!(x[3] >= 0.0)) continue;
            double vv[3] = {x[0], x[1], x[2]};
            double t = vf_worst(pl, m, vv);
            if (t < best) {
              best = t;
              bv[0] = vv[0];
              bv[1] = vv[1];
              bv[2] = vv[2];
              ok = 1;
            }
          }
        }
  if (!ok) return 1;
  for (int c = 0; c < 3; c++)
    v[c] = bv[c];
  *pt = best;
  return 0;
}

/* НАИМЕНЬШИЕ КВАДРАТЫ БЕЗ ОГРАНИЧЕНИЙ — второй этап §41, и он здесь на своём
 * месте: минимакс уже задал гарантию `t*`, а L2 выбирает внутри неё гладко.
 * Решатель безусловный (как `src/cut/qef.c`), поэтому результат ПРОВЕРЯЕТСЯ на
 * попадание в `t*`; не попал — берётся минимаксная точка, и это считается
 * (А123). */
static int vf_lsq(const double *pl, int m, double v[3]) {
  double A[4][4], b[4], x[4];
  memset(A, 0, sizeof A);
  memset(b, 0, sizeof b);
  for (int i = 0; i < m; i++) {
    const double *P = pl + (size_t)i * 4;
    for (int r = 0; r < 3; r++) {
      for (int c = 0; c < 3; c++)
        A[r][c] += P[r] * P[c];
      b[r] += P[r] * P[3];
    }
  }
  if (vf_solve(A, b, x, 3) != 0) return 1;
  for (int c = 0; c < 3; c++)
    v[c] = x[c];
  return 0;
}

int hz_poly_vfit(hz_polyset *ps, double dmove, hz_vfitstat *st) {
  memset(st, 0, sizeof *st);
  const int32_t nw = (int32_t)ps->nweld;
  if (nw <= 0 || ps->nbv <= 0) return 0;
  /* Владелец краевой вершины: полигон, которому принадлежит её петля. */
  int32_t *own = malloc((size_t)ps->nbv * sizeof *own);
  int32_t *cnt = calloc((size_t)nw + 1, sizeof *cnt);
  if (own == NULL || cnt == NULL) {
    free(own);
    free(cnt);
    return 2;
  }
  for (int32_t k = 0; k < ps->np; k++)
    for (int32_t l = ps->p[k].l0; l < ps->p[k].l0 + ps->p[k].nloop; l++)
      for (int32_t b = ps->loop[l]; b < ps->loop[l + 1]; b++)
        own[b] = k;
  for (int32_t b = 0; b < ps->nbv; b++) {
    int32_t w = ps->bw[b];
    if (w >= 0 && w < nw) cnt[w + 1]++;
  }
  for (int32_t i = 0; i < nw; i++)
    cnt[i + 1] += cnt[i];
  int32_t *lst = malloc((size_t)(cnt[nw] > 0 ? cnt[nw] : 1) * sizeof *lst);
  int32_t *cur = malloc((size_t)nw * sizeof *cur);
  if (lst == NULL || cur == NULL) {
    free(own);
    free(cnt);
    free(lst);
    free(cur);
    return 2;
  }
  memcpy(cur, cnt, (size_t)nw * sizeof *cur);
  for (int32_t b = 0; b < ps->nbv; b++) {
    int32_t w = ps->bw[b];
    if (w >= 0 && w < nw) lst[cur[w]++] = b;
  }
  free(cur);

  enum { VF_MAXPL = 16 };
  double pl[VF_MAXPL * 4];
  int32_t pk[VF_MAXPL];
  for (int32_t w = 0; w < nw; w++) {
    int32_t b0 = cnt[w], b1 = cnt[w + 1];
    if (b1 <= b0) continue;
    st->nvert++;
    /* Различные плоскости — по различным полигонам-владельцам. */
    int m = 0;
    for (int32_t i = b0; i < b1 && m < VF_MAXPL; i++) {
      int32_t k = own[lst[i]];
      int seen = 0;
      for (int j = 0; j < m; j++)
        if (pk[j] == k) seen = 1;
      if (seen) continue;
      pk[m] = k;
      const hz_poly *P = &ps->p[k];
      pl[(size_t)m * 4 + 0] = P->n[0];
      pl[(size_t)m * 4 + 1] = P->n[1];
      pl[(size_t)m * 4 + 2] = P->n[2];
      pl[(size_t)m * 4 + 3] = P->off;
      m++;
    }
    if (m > 8) st->nvert_many++;
    /* А121: ДВИГАЕМ ТОЛЬКО ВЕРШИНЫ С ТРЕМЯ И БОЛЕЕ ПЛОСКОСТЯМИ. При одной
     * свободны два направления в плоскости, при двух — одно вдоль их линии
     * пересечения; ограничение «не дальше `t*` от прежнего места» оставило бы
     * вершине круг или отрезок, и край поехал бы вдоль поверхности. */
    if (m < 3) {
      st->nvert_fixed++;
      continue;
    }
    /* Исходное положение вершины — из первого владельца, ДО перезаписи. */
    double v0[3];
    {
      int32_t b = lst[b0];
      const hz_poly *P = &ps->p[own[b]];
      hz_poly_world(P, ps->bv[(size_t)b * 2], ps->bv[(size_t)b * 2 + 1], v0);
    }
    double v[3], t = 0.0;
    if (vf_minimax(pl, m, v, &t) != 0) {
      st->nvert_fail++;
      continue;
    }
    /* ДВИГАТЬ ТОЛЬКО ЕСЛИ СТАНОВИТСЯ ЛУЧШЕ В ТОЙ ЖЕ НОРМЕ (А125). Порога тут не
     * нужно: сравниваются худшая невязка в НОВОЙ точке и в СТАРОЙ. Если старая
     * уже не хуже, двигать нечего — а именно так и бывает у почти параллельных
     * плоскостей, где пересечение уезжает далеко при крошечной невязке. */
    {
      double t0 = vf_worst(pl, m, v0);
      if (!(t < t0)) {
        st->nvert_noimp++;
        continue;
      }
    }
    double v2[3];
    if (vf_lsq(pl, m, v2) == 0) {
      /* L2 принимается ТОЛЬКО внутри гарантии минимакса (§41). */
      double t2 = vf_worst(pl, m, v2);
      if (t2 <= t) {
        for (int c = 0; c < 3; c++)
          v[c] = v2[c];
      } else {
        st->nvert_l2out++;
      }
    }
    /* ПРЕДЕЛ СМЕЩЕНИЯ (А124): `t*` ограничивает расстояние до плоскостей, а не
     * перемещение. У почти параллельных плоскостей пересечение уезжает далеко
     * при крошечной невязке — и край рвётся. Отказ, а не подтягивание: подтянуть
     * значило бы выдать за минимакс точку, которая им не является. */
    {
      double dd = 0.0;
      for (int c = 0; c < 3; c++)
        dd += (v[c] - v0[c]) * (v[c] - v0[c]);
      dd = sqrt(dd);
      if (dd > dmove) {
        st->nvert_far++;
        continue;
      }
      if (dd > st->dmax_move) st->dmax_move = dd;
    }
    st->nvert_moved++;
    if (t > st->tmax) st->tmax = t;
    /* Гистограмма `t*` — по долям допуска не считаем (допуск сюда не приходит),
     * а по абсолютной величине с логарифмическими корзинами: от `1e−6` м. */
    {
      int bk = 0;
      if (t > 0.0) {
        double lg = log10(t) + 6.0;
        bk = (int)lg;
        if (bk < 0) bk = 0;
        if (bk > 11) bk = 11;
      }
      st->t_hist[bk]++;
    }
    /* Новое положение кладётся в КАЖДЫЙ полигон своей проекцией: край хранится
     * в двумерных `(u, v)`, и полигон обязан остаться плоским (§42, п. 3).
     * Остаточная щель в этой вершине не больше `2 t*`. */
    for (int32_t i = b0; i < b1; i++) {
      int32_t b = lst[i];
      const hz_poly *P = &ps->p[own[b]];
      double q[3];
      for (int c = 0; c < 3; c++)
        q[c] = v[c] - P->org[c];
      ps->bv[(size_t)b * 2] = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
      ps->bv[(size_t)b * 2 + 1] = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
    }
  }
  /* Габарит края пересчитывается: по нему идёт отсев в луче (`pray.c`) и в
   * проверке перекрытия, и оставить его старым значило бы терять попадания. */
  for (int32_t k = 0; k < ps->np; k++) {
    hz_poly *P = &ps->p[k];
    if (P->nloop <= 0) continue;
    double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
    for (int32_t l = P->l0; l < P->l0 + P->nloop; l++)
      for (int32_t b = ps->loop[l]; b < ps->loop[l + 1]; b++)
        for (int c = 0; c < 2; c++) {
          double x = ps->bv[(size_t)b * 2 + (size_t)c];
          if (x < lo[c]) lo[c] = x;
          if (x > hi[c]) hi[c] = x;
        }
    if (lo[0] <= hi[0]) {
      for (int c = 0; c < 2; c++) {
        P->uvlo[c] = lo[c];
        P->uvhi[c] = hi[c];
      }
    }
  }
  free(own);
  free(cnt);
  free(lst);
  return 0;
}
