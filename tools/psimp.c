/* psimp — УПРОЩЕНИЕ КРАЯ ПРИБЛИЖЕНИЕМ (третий ход между слиянием и подгонкой).
 *
 * PLAN_ELEMENTS.md, Ш7. Слияние оставляет край КАК ЕСТЬ — объединение краёв
 * входных полигонов, со всеми их вершинами (у зала 57 на полигон). Подгонка
 * (edge-collapse) край упростила бы, но Ш7 показал, что по ЧИСЛУ ПОЛИГОНОВ она
 * ничего не даёт: дно задаёт топология поверхности. Остаётся третий ход —
 * ПРИБЛИЗИТЬ КРАЙ, не трогая ни плоскость, ни разбиение.
 *
 * ЗАКОННОСТЬ ТА ЖЕ, ЧТО У ПЛОСКОСТИ. Вершина края, сдвинутая меньше чем на `δ`,
 * укладывается в тот же допуск, которым уже пользуются `dmax` и сегментация.
 * Отдельного порога не заводится.
 *
 * ЗАЧЕМ ЭТО НУЖНО — ЦЕНА ЛУЧА, А НЕ ПАМЯТЬ. `hz_poly_inside` перебирает ВСЕ
 * вершины края на КАЖДЫЙ тест луча с полигоном, то есть край стоит во
 * внутреннем цикле трассировщика. Ш6 намерил кадр `48.2` с против `7.1` при
 * сопоставимом числе полигонов, и разница шла от формы, а не от счёта.
 *
 * Метод — Дуглас–Пекер по замкнутой петле. Он НЕ двигает вершины (в отличие от
 * edge-collapse), а только ВЫБРАСЫВАЕТ те, что лежат ближе `δ` к хорде, поэтому
 * ошибка ограничена `δ` по построению и проверять её отдельно не нужно —
 * достаточно проверить, что площадь не уехала.
 */

#include "poly_seg.h"
#include "polygon.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Расстояние от точки до отрезка в (u,v). */
static double seg_dist(const double *p, const double *a, const double *b) {
  double du = b[0] - a[0], dv = b[1] - a[1];
  double L2 = du * du + dv * dv;
  double t = 0.0;
  if (L2 > 0.0) {
    t = ((p[0] - a[0]) * du + (p[1] - a[1]) * dv) / L2;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
  }
  double qu = a[0] + t * du - p[0], qv = a[1] + t * dv - p[1];
  return sqrt(qu * qu + qv * qv);
}

/* Дуглас–Пекер на открытой цепи [i0, i1] массива `v`; помечает `keep`. */
static void dp(const double *v, int32_t i0, int32_t i1, double tol, unsigned char *keep) {
  if (i1 <= i0 + 1) return;
  double dmax = -1.0;
  int32_t im = i0;
  for (int32_t i = i0 + 1; i < i1; i++) {
    double d = seg_dist(v + (size_t)i * 2, v + (size_t)i0 * 2, v + (size_t)i1 * 2);
    if (d > dmax) {
      dmax = d;
      im = i;
    }
  }
  if (dmax <= tol) return;
  keep[im] = 1;
  dp(v, i0, im, tol, keep);
  dp(v, im, i1, tol, keep);
}

int main(int argc, char **argv) {
  double delta = (argc > 1) ? strtod(argv[1], NULL) : 0.045;
  hz_objmesh m;
  if (hz_obj_load(&m, HZ_CFG_HALL_OBJ, HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет %s\n", HZ_CFG_HALL_OBJ);
    return 1;
  }
  hz_pseglist sg;
  if (hz_seg_planar(&sg, &m, delta) != 0) return 1;
  hz_polyset ps;
  if (hz_poly_build(&ps, &m, &sg) != 0) return 1;

  printf("== упрощение края приближением: зал, δ = %g м, полигонов %d\n", delta, ps.np);
  printf("   %-10s %10s %10s %9s %12s %12s %11s\n", "допуск, м", "вершин", "на полигон", "осталось",
         "Σ шнуровка", "было", "отн. ошибка");

  /* Исходная площадь по шнуровке — эталон: она же служит проверкой, что
   * упрощение не съело поверхность. */
  double lace0 = 0.0;
  int32_t nv0 = ps.nbv;
  for (int32_t k = 0; k < ps.np; k++) {
    const hz_poly *P = &ps.p[k];
    for (int32_t l = P->l0; l < P->l0 + P->nloop; l++) {
      int32_t b = ps.loop[l], e = ps.loop[l + 1], n = e - b;
      for (int32_t q = 0; q < n; q++) {
        const double *A = ps.bv + (size_t)(b + q) * 2;
        const double *B = ps.bv + (size_t)(b + (q + 1) % n) * 2;
        lace0 += 0.5 * (A[0] * B[1] - B[0] * A[1]);
      }
    }
  }

  const double tols[6] = {1e-9, 0.001, 0.005, 0.02, 0.045, 0.15};
  unsigned char *keep = malloc((size_t)ps.nbv);
  if (keep == NULL) return 1;
  for (int it = 0; it < 6; it++) {
    double tol = tols[it];
    int32_t nkeep = 0;
    double lace = 0.0;
    for (int32_t k = 0; k < ps.np; k++) {
      const hz_poly *P = &ps.p[k];
      for (int32_t l = P->l0; l < P->l0 + P->nloop; l++) {
        int32_t b = ps.loop[l], e = ps.loop[l + 1], n = e - b;
        if (n < 3) continue;
        memset(keep + b, 0, (size_t)n);
        /* ЗАМКНУТАЯ петля: два опорных узла — самая дальняя пара по первой
         * вершине, чтобы результат не зависел от того, где петля начинается. */
        int32_t far = 0;
        double dm = -1.0;
        for (int32_t q = 1; q < n; q++) {
          double du = ps.bv[(size_t)(b + q) * 2] - ps.bv[(size_t)b * 2];
          double dv = ps.bv[(size_t)(b + q) * 2 + 1] - ps.bv[(size_t)b * 2 + 1];
          double d = du * du + dv * dv;
          if (d > dm) {
            dm = d;
            far = q;
          }
        }
        keep[b] = 1;
        keep[b + far] = 1;
        /* КОНЕЦ ВТОРОЙ ЦЕПИ НАДО ПОМЕТИТЬ ЯВНО: Дуглас–Пекер метит только
         * ВНУТРЕННИЕ узлы, опорные обязан задать вызывающий. Без этого
         * последняя вершина петли терялась вместе со всем, что лежало в
         * допуске от её хорды. */
        keep[b + n - 1] = 1;
        dp(ps.bv, b, b + far, tol, keep + 0);
        dp(ps.bv, b + far, b + n - 1, tol, keep + 0);
        /* последний отрезок замыкает петлю на b */
        int32_t prev = -1, first = -1;
        for (int32_t q = 0; q < n; q++)
          if (keep[b + q]) {
            nkeep++;
            if (first < 0) first = b + q;
            if (prev >= 0) {
              const double *A = ps.bv + (size_t)prev * 2;
              const double *B = ps.bv + (size_t)(b + q) * 2;
              lace += 0.5 * (A[0] * B[1] - B[0] * A[1]);
            }
            prev = b + q;
          }
        if (first >= 0 && prev >= 0) {
          const double *A = ps.bv + (size_t)prev * 2;
          const double *B = ps.bv + (size_t)first * 2;
          lace += 0.5 * (A[0] * B[1] - B[0] * A[1]);
        }
      }
    }
    printf("   %-10g %10d %10.1f %8.1f%% %12.3f %12.3f %11.3e\n", tol, nkeep, (double)nkeep / ps.np,
           100.0 * nkeep / (double)nv0, lace, lace0, fabs(lace - lace0) / fabs(lace0));
  }
  free(keep);
  hz_poly_free(&ps);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
