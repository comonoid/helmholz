/* scenechk — ПРИГОДНОСТЬ СЦЕНЫ, ОДНОЙ КОМАНДОЙ (PLAN_ELEMENTS.md, §105).
 *
 * ЗАЧЕМ. Два дня ушло на борьбу с дефектами, которые оказались свойствами
 * ВХОДА, а не схемы: разброс площади треугольника в 250 млн раз (от 0.08 мм² до
 * 19.7 м²) убил подряд ограничение размера элемента, пространственное дерево и
 * нижнюю границу элемента; 12.63 %% рёбер с одним владельцем сделали
 * непрерывное поле невозможным вдоль трети границ. Ни одно из этих чисел не
 * стоило больше двадцати строк, и ни одно не было снято вовремя.
 *
 * ЧТО МЕРИТСЯ И ПОЧЕМУ ИМЕННО ЭТО:
 *   - РЁБРА С ОДНИМ ВЛАДЕЛЬЦЕМ: открытый край. Вдоль него непрерывность
 *     невозможна — сшивать не с чем;
 *   - РЁБРА БОЛЕЕ ЧЕМ С ДВУМЯ: неманифолдный стык, «сосед» неоднозначен;
 *   - РАЗБРОС ПЛОЩАДИ ТРЕУГОЛЬНИКА: нижняя граница размера элемента задаётся
 *     входом, и ограничение размера бессильно против одного крупного
 *     треугольника (§92.1.4);
 *   - ГАБАРИТ и ЧИСЛО: цена сборки растёт как квадрат числа элементов.
 *
 * Запуск: `scenechk путь.obj масштаб`. */
#include <math.h>
#include "scene_obj.h"
#include <stdio.h>
#include <stdlib.h>
static int cmpk(const void *a, const void *b) {
  int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}
int main(int argc, char **argv) {
  if (argc < 3) return 1;
  hz_objmesh m;
  if (hz_obj_load(&m, argv[1], atof(argv[2])) != 0) {
    printf("нет сцены %s\n", argv[1]);
    return 1;
  }
  int64_t ne = (int64_t)m.nt * 3;
  int64_t *e = malloc((size_t)ne * sizeof *e);
  if (e == NULL) return 1;
  for (int32_t t = 0; t < m.nt; t++)
    for (int i = 0; i < 3; i++) {
      int64_t a = m.f[(size_t)t * 3 + (size_t)i], b = m.f[(size_t)t * 3 + (size_t)((i + 1) % 3)];
      int64_t lo = a < b ? a : b, hi = a < b ? b : a;
      e[(size_t)t * 3 + (size_t)i] = lo * 16000000LL + hi;
    }
  qsort(e, (size_t)ne, sizeof *e, cmpk);
  int64_t c1 = 0, c2 = 0, cm = 0, tot = 0;
  for (int64_t i = 0; i < ne;) {
    int64_t j = i;
    while (j < ne && e[j] == e[i])
      j++;
    int64_t c = j - i;
    tot++;
    if (c == 1)
      c1++;
    else if (c == 2)
      c2++;
    else
      cm++;
    i = j;
  }
  double amin = 1e300, amax = 0.0, asum = 0.0;
  for (int32_t t = 0; t < m.nt; t++) {
    const double *p0 = m.v + 3 * (size_t)m.f[(size_t)t * 3 + 0];
    const double *p1 = m.v + 3 * (size_t)m.f[(size_t)t * 3 + 1];
    const double *p2 = m.v + 3 * (size_t)m.f[(size_t)t * 3 + 2];
    double u[3], v[3], c[3];
    for (int a = 0; a < 3; a++) {
      u[a] = p1[a] - p0[a];
      v[a] = p2[a] - p0[a];
    }
    c[0] = u[1] * v[2] - u[2] * v[1];
    c[1] = u[2] * v[0] - u[0] * v[2];
    c[2] = u[0] * v[1] - u[1] * v[0];
    double ar = 0.5 * sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    if (ar < amin) amin = ar;
    if (ar > amax) amax = ar;
    asum += ar;
  }
  printf("%s: треугольников %d, вершин %d, нормалей %d, материалов %d\n", argv[1], m.nt, m.nv,
         m.nvn, m.nmtl);
  printf("  рёбер %lld; с ОДНИМ владельцем %lld (%.2f %%), с двумя %lld, больше двух %lld\n",
         (long long)tot, (long long)c1, 100.0 * (double)c1 / (double)tot, (long long)c2,
         (long long)cm);
  printf("  площадь треугольника: min %.3e, max %.3e, средняя %.3e м²; вся %.1f м²\n", amin, amax,
         asum / (double)m.nt, asum);
  printf("  габарит: %.2f x %.2f x %.2f м\n", m.hi[0] - m.lo[0], m.hi[1] - m.lo[1],
         m.hi[2] - m.lo[2]);
  free(e);
  hz_obj_free(&m);
  return 0;
}
