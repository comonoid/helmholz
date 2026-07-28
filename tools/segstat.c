/* segstat — Ш1: ЧИСЛО ПОЛИГОНОВ. Свип по δ на скачанных сценах.
 *
 * PLAN_ELEMENTS.md, Ш1. Печатает статистику входа и, для каждого δ, число
 * плоских участков, найденных сегментатором пространственной смежности.
 *
 * Запуск: build/segstat FILE.obj SCALE [delta ...]
 */

#include "poly_seg.h"
#include "polygon.h"
#include "scene_obj.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: %s FILE.obj SCALE [delta ...]\n", argv[0]);
    return 2;
  }
  double scale = strtod(argv[2], NULL);
  double t0 = now_s();
  hz_objmesh m;
  int rc = hz_obj_load(&m, argv[1], scale);
  if (rc != 0) {
    fprintf(stderr, "hz_obj_load: rc=%d\n", rc);
    return 1;
  }
  double t1 = now_s();

  printf("== %s (масштаб ×%g), разбор %.2f с\n", argv[1], scale, t1 - t0);
  printf("   v = %d, vn = %d, треугольников = %d (вырожденных отброшено %lld,\n"
         "   граней >3 вершин %lld), материалов = %d\n",
         m.nv, m.nvn, m.nt, (long long)m.ndegen, (long long)m.nquad, m.nmtl);
  printf("   габарит = %.3f × %.3f × %.3f м, начало (%.3f, %.3f, %.3f)\n", m.hi[0] - m.lo[0],
         m.hi[1] - m.lo[1], m.hi[2] - m.lo[2], m.lo[0], m.lo[1], m.lo[2]);

  /* Распределение РАЗМЕРА треугольника: от него зависит шаг сетки поиска, и
   * брать среднее нельзя — у него хвост. */
  double *ext = malloc((size_t)m.nt * sizeof *ext);
  if (ext == NULL) return 1;
  double area = 0.0;
  for (int32_t t = 0; t < m.nt; t++) {
    double p[3][3];
    hz_obj_tri(&m, t, p);
    double e = 0.0;
    for (int a = 0; a < 3; a++) {
      double lo = p[0][a], hi = p[0][a];
      for (int i = 1; i < 3; i++) {
        if (p[i][a] < lo) lo = p[i][a];
        if (p[i][a] > hi) hi = p[i][a];
      }
      if (hi - lo > e) e = hi - lo;
    }
    ext[t] = e;
    area += hz_obj_tri_area(&m, t);
  }
  qsort(ext, (size_t)m.nt, sizeof *ext, cmp_d);
  printf("   площадь = %.1f м²; габарит треугольника, м: p10 %.4g  p50 %.4g  p90 %.4g\n"
         "   p99 %.4g  p99.9 %.4g  max %.4g\n",
         area, ext[m.nt / 10], ext[m.nt / 2], ext[m.nt * 9 / 10], ext[(int32_t)(m.nt * 0.99)],
         ext[(int32_t)(m.nt * 0.999)], ext[m.nt - 1]);
  free(ext);

  /* --- свип по δ ------------------------------------------------------------ */
  if (argc > 3) {
    printf("\n      δ, м    участков     ср.тр/уч   max dmax/δ  раунд  отпущ.  "
           "межур.пар   канд/тр    время, с\n");
    double d0 = 0.0, n0 = 0.0, d1 = 0.0, n1 = 0.0;
    for (int i = 3; i < argc; i++) {
      double delta = strtod(argv[i], NULL);
      hz_pseglist s;
      double ta = now_s();
      int src = hz_seg_planar(&s, &m, delta);
      double tb = now_s();
      if (src != 0) {
        fprintf(stderr, "hz_seg_planar(%g): rc=%d\n", delta, src);
        continue;
      }
      /* ИНВАРИАНТЫ РАЗБИЕНИЯ. Без них число участков — просто число: разметка,
       * потерявшая треугольники, дала бы участков МЕНЬШЕ и выглядела бы лучше. */
      int64_t sumtri = 0;
      double sumarea = 0.0;
      int bad = 0;
      for (int32_t k = 0; k < s.nseg; k++) {
        sumtri += s.seg[k].ntri;
        sumarea += s.seg[k].area;
      }
      for (int32_t t = 0; t < m.nt; t++)
        if (s.label[t] < 0 || s.label[t] >= s.nseg) bad++;
      if (bad != 0 || sumtri != m.nt || fabs(sumarea - area) > 1e-6 * area) {
        printf("  !! РАЗБИЕНИЕ НЕВЕРНО при δ=%g: вне диапазона %d, Σтр %lld против %d, "
               "Σплощадь %.6f против %.6f\n",
               delta, bad, (long long)sumtri, m.nt, sumarea, area);
        hz_seg_free(&s);
        continue;
      }
      double worst = 0.0;
      for (int32_t k = 0; k < s.nseg; k++)
        if (s.seg[k].dmax / delta > worst) worst = s.seg[k].dmax / delta;
      /* Сборка полигонов гоняется здесь же: на масштабной сцене надо знать не
       * только СКОЛЬКО участков, но и во что обходится край. */
      hz_polyset ps;
      double tc = now_s();
      int prc = hz_poly_build(&ps, &m, &s);
      double td = now_s();
      printf("  %9.4g   %9d   %10.1f   %10.4f  %5lld  %6lld  %9lld  %8.1f  %10.2f\n", delta, s.nseg,
             (double)m.nt / (double)s.nseg, worst, (long long)s.nrounds, (long long)s.nreleased,
             (long long)s.ncross, (double)s.ncand / (double)m.nt, tb - ta);
      if (prc == 0) {
        int32_t noloop = 0;
        for (int32_t k = 0; k < ps.np; k++)
          if (ps.p[k].nloop == 0) noloop++;
        printf("             край: %.2f с, петель %d, вершин края %d (%.1f на полигон), "
               "без петли %d, не замкнулось %lld, развилок %lld;\n"
               "             полурёбер без пары %.3f%%, полигонов на шве материалов %lld\n",
               td - tc, ps.nloopall, ps.nbv, (double)ps.nbv / (double)ps.np, noloop,
               (long long)ps.nopen, (long long)ps.nfork,
               100.0 * (double)ps.nhe_open / (double)ps.nhe, (long long)ps.nmixed);
        hz_poly_free(&ps);
      } else {
        printf("             край: сборка полигонов не прошла (rc=%d)\n", prc);
      }
      if (n0 <= 0.0) {
        d0 = delta;
        n0 = s.nseg;
      }
      d1 = delta;
      n1 = s.nseg;
      hz_seg_free(&s);
    }
    if (n0 > 0.0 && n1 > 0.0 && d1 > 0.0 && d0 > 0.0 && fabs(log(d1 / d0)) > 0.0)
      printf("  показатель по крайним точкам: N ∝ δ^%.2f\n", log(n1 / n0) / log(d1 / d0));
  }

  hz_obj_free(&m);
  return 0;
}
