/* pff.c — реализация; довод, зачем отдельный файл, и вся физика — в `pff.h`.
 *
 * ПЕРЕНОС ИЗ `src/plink.c` СДЕЛАН БЕЗ ЕДИНОЙ ПРАВКИ АРИФМЕТИКИ, и это
 * проверяется числом: самопроверка обязана давать ту же худшую невязку
 * `1.110e-16` (`./build/prad`, первая строка вывода).
 */
#include "pff.h"
#include <math.h>

#define PFF_PI 3.14159265358979323846

double hz_pff_loop(const double *R, int n, const double nv[3]) {
  double s = 0.0;
  for (int k = 0; k < n; k++) {
    const double *a = R + 3 * k, *b = R + 3 * ((k + 1) % n);
    double c[3] = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    double cl = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    if (!(cl > 0.0)) continue; /* совпавшие или противоположные лучи угла не задают */
    /* УГОЛ БЕРЁТСЯ ЧЕРЕЗ `atan2`, А НЕ ЧЕРЕЗ `acos`, И ЭТО НЕ КОСМЕТИКА.
     *
     * `acos(a·b/|a||b|)` теряет половину разрядов при МАЛОМ угле: производная
     * `acos` у единицы бесконечна, и погрешность аргумента `ε` даёт погрешность
     * угла `√(2ε)`. Работаем мы как раз там — солнце занимает `0.533°`, у
     * далёкого источника соседние лучи почти совпадают.
     *
     * `atan2(|a×b|, a·b)` = `atan2(|a||b|sinθ, |a||b|cosθ)` устойчива при ЛЮБОМ
     * угле: длины сокращаются сами, поэтому не нужны ни `|a|`, ни `|b|`, ни
     * деление, ни зажим аргумента в `[-1, 1]` — три источника ошибки и одна
     * ветвь исчезают вместе с ними. `|a×b|` уже посчитан как `cl`.
     *
     * Замерено самопроверкой (`hz_pff_selftest`): худшая невязка была
     * `1.110e-16` на `acos`, стала `1.110e-16` на `atan2` — на ЭТОМ случае
     * разницы нет, потому что квадрат `a = h` даёт углы порядка радиана. Ход
     * сделан ради случая МАЛОГО угла, который самопроверкой не покрыт, и это
     * записано здесь, чтобы «числа не изменились» не прочли как «правка
     * бесполезна». */
    double beta = atan2(cl, a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
    s += beta * (nv[0] * c[0] + nv[1] * c[1] + nv[2] * c[2]) / cl;
  }
  /* Минус — довод о знаке в `pff.h`, и он не косметический. */
  return -s / (2.0 * PFF_PI);
}

int hz_pff_clip_half(const double *R, int n, const double nv[3], double *out) {
  int m = 0;
  for (int k = 0; k < n; k++) {
    const double *a = R + 3 * k, *b = R + 3 * ((k + 1) % n);
    double da = nv[0] * a[0] + nv[1] * a[1] + nv[2] * a[2];
    double db = nv[0] * b[0] + nv[1] * b[1] + nv[2] * b[2];
    if (da > 0.0) {
      for (int c = 0; c < 3; c++)
        out[3 * m + c] = a[c];
      m++;
    }
    if ((da > 0.0) != (db > 0.0)) {
      double t = da / (da - db);
      for (int c = 0; c < 3; c++)
        out[3 * m + c] = a[c] + t * (b[c] - a[c]);
      m++;
    }
  }
  return m;
}

double hz_pff_point_poly(const double x[3], const double nrec[3], const double *P, int n,
                         double *w) {
  if (n < 3) return 0.0;
  /* Вершины сдвигаются в раму приёмника ПРЯМО В БУФЕР ОТСЕЧЕНИЯ нельзя: у
   * Сазерленда–Ходжмена вход и выход обязаны быть разными массивами. Поэтому
   * сдвиг делается на месте в `w` только после того, как отсечение прочло
   * исходные вершины, — здесь проще: сдвинутые кладём в `w`, отсечённые в
   * хвост того же `w` за границей `3·n`. */
  double *rel = w;
  double *clp = w + 3 * n;
  for (int k = 0; k < n; k++)
    for (int c = 0; c < 3; c++)
      rel[3 * k + c] = P[3 * k + c] - x[c];
  int m = hz_pff_clip_half(rel, n, nrec, clp);
  if (m < 3) return 0.0;
  return hz_pff_loop(clp, m, nrec);
}

double hz_pff_selftest(double *f_big, double *f_unit, double *ref_unit) {
  double x[3] = {0.0, 0.0, 0.0}, nrec[3] = {0.0, 0.0, 1.0};
  /* Излучатель смотрит ВНИЗ, на приёмника: `n_j = (0,0,−1)`, и рама берётся с
   * `eu × ev = n_j`, как у настоящего полигона. */
  double eu[3] = {1.0, 0.0, 0.0}, ev[3] = {0.0, -1.0, 0.0};
  double P[12], w[27], worst = 0.0;
  const double h = 1.0;
  for (int cs = 0; cs < 2; cs++) {
    double a = (cs == 0) ? 1e6 : 1.0;
    double org[3] = {0.0, 0.0, h};
    double uv[8] = {-a, -a, a, -a, a, a, -a, a}; /* против часовой в раме (eu,ev) */
    for (int k = 0; k < 4; k++)
      for (int c = 0; c < 3; c++)
        P[3 * k + c] = org[c] + uv[2 * k] * eu[c] + uv[2 * k + 1] * ev[c];
    double f = hz_pff_point_poly(x, nrec, P, 4, w);
    double X = a / h, t = X / sqrt(1.0 + X * X);
    double ref = (4.0 / PFF_PI) * t * atan(t);
    if (cs == 0) {
      *f_big = f;
    } else {
      *f_unit = f;
      *ref_unit = ref;
    }
    double e = fabs(f - ref);
    if (e > worst) worst = e;
  }
  return worst;
}
