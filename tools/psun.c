/* psun.c — ШАГ О73 (Ф3): ФРОНТ ОТ СОЛНЦА С СОСТОЯНИЕМ НА ГРАНЯХ.
 * План §266, аудит §267, поправки 08-05 (огрублять, а не отсекать).
 *
 * ЗАПУСК:
 *   psun ФАЙЛ.obj МАСШТАБ [leaf=N] [lev=N] [grade=0|1] [px=9.1]
 *        [sun=dx,dy,dz] [shadow=0]
 *   `shadow=0` — НЕГАТИВНЫЙ КОНТРОЛЬ: ни заслонения, ни огрубления по ровности.
 *   Обязан воспроизвести числа О72 ДО ЕДИНИЦЫ; расхождение означает, что
 *   изменилось что-то ещё.
 */
#include "pclip.h"
#include "pfront.h"
#include "pocc.h"
#include "ptree.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: psun ФАЙЛ.obj МАСШТАБ [leaf=N] [lev=N] [grade=0|1] [px=9.1] "
                    "[sun=dx,dy,dz] [shadow=0]\n");
    return 1;
  }
  int leafmax = 0, maxlev = 0, grade = 1, shadow = 1;
  double px = 9.1, sdir[3] = {0.3, -0.9, 0.3};
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "leaf=", 5) == 0) leafmax = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "lev=", 4) == 0) maxlev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "grade=", 6) == 0) grade = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "px=", 3) == 0) px = strtod(argv[i] + 3, NULL);
    if (strncmp(argv[i], "shadow=", 7) == 0) shadow = (int)strtol(argv[i] + 7, NULL, 10);
    if (strncmp(argv[i], "sun=", 4) == 0) {
      char *e = NULL;
      sdir[0] = strtod(argv[i] + 4, &e);
      if (e != NULL && *e == ',') sdir[1] = strtod(e + 1, &e);
      if (e != NULL && *e == ',') sdir[2] = strtod(e + 1, NULL);
    }
  }
  hz_objmesh m;
  double t0 = now_s();
  if (hz_obj_load(&m, argv[1], strtod(argv[2], NULL)) != 0) {
    fprintf(stderr, "нет сцены %s\n", argv[1]);
    return 1;
  }
  printf("== СЦЕНА %s: треугольников %d, габарит %.2f x %.2f x %.2f м, чтение %.1f с\n", argv[1],
         m.nt, m.hi[0] - m.lo[0], m.hi[1] - m.lo[1], m.hi[2] - m.lo[2], now_s() - t0);
  hz_ptree T;
  t0 = now_s();
  if (hz_ptree_build_cut(&T, &m, leafmax, maxlev, grade) != 0) {
    fprintf(stderr, "отказ дерева\n");
    hz_obj_free(&m);
    return 2;
  }
  printf("== ДЕРЕВО за %.2f с: узлов %d, листьев %lld, глубина %d\n", now_s() - t0, T.nnd,
         (long long)T.nleaf, T.depth);

  double *tlo = malloc(3 * (size_t)m.nt * sizeof *tlo);
  double *thi = malloc(3 * (size_t)m.nt * sizeof *thi);
  double *tpl = malloc(4 * (size_t)m.nt * sizeof *tpl);
  int32_t *list = malloc((size_t)m.nt * sizeof *list);
  unsigned char *occ = calloc((size_t)T.nnd, 1);
  if (tlo == NULL || thi == NULL || tpl == NULL || list == NULL || occ == NULL) {
    free(tlo);
    free(thi);
    free(tpl);
    free(list);
    free(occ);
    hz_ptree_free(&T);
    hz_obj_free(&m);
    return 2;
  }
  for (int32_t t = 0; t < m.nt; t++) {
    const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
    const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
    const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
    double u[3], w[3], nr[3];
    for (int c = 0; c < 3; c++) {
      u[c] = B[c] - A[c];
      w[c] = C[c] - A[c];
    }
    nr[0] = u[1] * w[2] - u[2] * w[1];
    nr[1] = u[2] * w[0] - u[0] * w[2];
    nr[2] = u[0] * w[1] - u[1] * w[0];
    double L = sqrt(nr[0] * nr[0] + nr[1] * nr[1] + nr[2] * nr[2]);
    if (!(L > 0.0)) L = 1.0;
    double *pl = tpl + 4 * (size_t)t;
    for (int c = 0; c < 3; c++)
      pl[c] = nr[c] / L;
    pl[3] = pl[0] * A[0] + pl[1] * A[1] + pl[2] * A[2];
    double *bl = tlo + 3 * (size_t)t, *bh = thi + 3 * (size_t)t;
    for (int c = 0; c < 3; c++) {
      bl[c] = 1e300;
      bh[c] = -1e300;
    }
    for (int i = 0; i < 3; i++) {
      const double *p = m.v + 3 * (size_t)m.f[3 * (size_t)t + (size_t)i];
      for (int c = 0; c < 3; c++) {
        if (p[c] < bl[c]) bl[c] = p[c];
        if (p[c] > bh[c]) bh[c] = p[c];
      }
    }
    list[t] = t;
  }
  t0 = now_s();
  if (hz_pocc_build(&m, &T, tlo, thi, tpl, occ, NULL) != 0) {
    fprintf(stderr, "отказ занятости\n");
    free(tlo);
    free(thi);
    free(tpl);
    free(list);
    free(occ);
    hz_ptree_free(&T);
    hz_obj_free(&m);
    return 2;
  }
  int64_t nocc = 0;
  for (int32_t i = 0; i < T.nnd; i++)
    if (occ[i]) nocc++;
  printf("== ЗАНЯТОСТЬ за %.2f с: узлов с геометрией %lld (%.2f %%)\n", now_s() - t0,
         (long long)nocc, 100.0 * (double)nocc / (double)T.nnd);

  hz_pfront_ctx X;
  memset(&X, 0, sizeof X);
  X.T = &T;
  X.m = &m;
  X.tlo = tlo;
  X.thi = thi;
  X.occ = occ;
  X.noshadow = !shadow;
  double L = sqrt(sdir[0] * sdir[0] + sdir[1] * sdir[1] + sdir[2] * sdir[2]);
  if (!(L > 0.0)) L = 1.0;
  for (int c = 0; c < 3; c++)
    X.dir[c] = sdir[c] / L;
  X.pxeps2 = (px * HZ_CFG_EPS) * (px * HZ_CFG_EPS);
  X.t_entry = 1e300;
  for (int k = 0; k < 8; k++) {
    double t = 0.0;
    for (int c = 0; c < 3; c++)
      t += ((k & (1 << c)) ? T.nd[0].hi[c] : T.nd[0].lo[c]) * X.dir[c];
    if (t < X.t_entry) X.t_entry = t;
  }
  printf("== ФРОНТ: направление %.3f,%.3f,%.3f, пол %g px ОТ ИСТОЧНИКА, потолок огрубления %d, "
         "заслонение %s\n",
         X.dir[0], X.dir[1], X.dir[2], px, HZ_PFRONT_COARSEN_MAX, shadow ? "ВКЛ" : "ВЫКЛ");
  printf("== ОГОВОРКА ОПЕРАТОРА: %s\n", hz_pfront_note());

  hz_pfront_face in[3], out[3];
  for (int a = 0; a < 3; a++) {
    in[a].c0 = 1.0; /* на входе в сцену диск источника открыт целиком */
    in[a].cu = 0.0;
    in[a].cv = 0.0;
  }
  t0 = now_s();
  hz_pfront_walk(&X, 0, T.nd[0].lo, T.nd[0].hi, in, out, list, m.nt, 0);
  double secs = now_s() - t0;
  if (X.fail) {
    fprintf(stderr, "отказ обхода\n");
    free(tlo);
    free(thi);
    free(tpl);
    free(list);
    free(occ);
    hz_ptree_free(&T);
    hz_obj_free(&m);
    return 2;
  }
  printf("   ЯЧЕЕК %lld (с геометрией %lld, пустых %lld)\n", (long long)X.ncell,
         (long long)X.ncell_geo, (long long)X.ncell_void);
  printf("   спуск остановлен: ПУСТОТОЙ %lld (%.2f %%), ЛИСТОМ %lld (%.2f %%), ПОЛОМ %lld "
         "(%.2f %%), ТЕМНЫМ %lld (%.2f %%)\n",
         (long long)X.stop_void, 100.0 * (double)X.stop_void / (double)X.ncell,
         (long long)X.stop_leaf, 100.0 * (double)X.stop_leaf / (double)X.ncell,
         (long long)X.stop_floor, 100.0 * (double)X.stop_floor / (double)X.ncell,
         (long long)X.stop_flat, 100.0 * (double)X.stop_flat / (double)X.ncell);
  printf("   ячеек ОСВЕЩЁННЫХ %lld, В ТЕНИ %lld; потолок огрубления сработал %lld раз; предел "
         "кусков %lld раз\n",
         (long long)X.nlit, (long long)X.nshadow, (long long)X.ncap, (long long)X.npclip);
  printf("   ВРЕМЯ %.2f с, на ячейку %.1f нс\n", secs,
         X.ncell > 0 ? 1e9 * secs / (double)X.ncell : 0.0);
  free(tlo);
  free(thi);
  free(tpl);
  free(list);
  free(occ);
  hz_ptree_free(&T);
  hz_obj_free(&m);
  return 0;
}
