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
#include "pmark.h"
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

#define HZ_REF_RAYS 16 /* лучей на точку выборки: диск солнца мал, шестнадцать дают шаг ~2 % */

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
  int leafmax = 0, maxlev = 0, grade = 1, shadow = 1, cam = 0, nref = 0, cliplev = -1, camfull = 0,
      usemark = 0, relax = HZ_PMARK_OUT_LEVELS;
  double camo[3] = {0.0, 0.0, 0.0}, camf[3] = {0.0, 0.0, 1.0};
  double px = 9.1, sdir[3] = {0.3, -0.9, 0.3};
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "leaf=", 5) == 0) leafmax = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "lev=", 4) == 0) maxlev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "grade=", 6) == 0) grade = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "px=", 3) == 0) px = strtod(argv[i] + 3, NULL);
    if (strncmp(argv[i], "shadow=", 7) == 0) shadow = (int)strtol(argv[i] + 7, NULL, 10);
    if (strncmp(argv[i], "cam=", 4) == 0) cam = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "ref=", 4) == 0) nref = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "clip=", 5) == 0) cliplev = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "camfull=", 8) == 0) camfull = (int)strtol(argv[i] + 8, NULL, 10);
    if (strncmp(argv[i], "mark=", 5) == 0) usemark = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "relax=", 6) == 0) relax = (int)strtol(argv[i] + 6, NULL, 10);
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
  X.cliplev = cliplev;
  X.frustfull = camfull;
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
  /* ПИРАМИДА КАМЕРЫ (ключ `cam=1`). Глаз и цель — из `scene_cfg.h`, где камеры
   * НАЙДЕНЫ ЗАМЕРОМ, а не назначены (§187); поле зрения и разрешение оттуда же.
   * Плоскости строятся нормалями ВНУТРЬ. */
  if (cam || usemark) {
    double e[3], at[3];
    if (strstr(argv[1], "conference") != NULL) {
      double a1[3] = HZ_CFG_HALL_EYE, a2[3] = HZ_CFG_HALL_AT;
      memcpy(e, a1, sizeof e);
      memcpy(at, a2, sizeof at);
    } else if (strstr(argv[1], "rungholt") != NULL) {
      double a1[3] = HZ_CFG_CITY_EYE, a2[3] = HZ_CFG_CITY_AT;
      memcpy(e, a1, sizeof e);
      memcpy(at, a2, sizeof at);
    } else {
      double a1[3] = HZ_CFG_MIGUEL_EYE, a2[3] = HZ_CFG_MIGUEL_AT;
      memcpy(e, a1, sizeof e);
      memcpy(at, a2, sizeof at);
    }
    double f[3], up[3] = {0.0, 1.0, 0.0}, r[3], u2[3];
    double ln = 0.0;
    for (int c = 0; c < 3; c++) {
      f[c] = at[c] - e[c];
      ln += f[c] * f[c];
    }
    ln = sqrt(ln);
    if (!(ln > 0.0)) ln = 1.0;
    for (int c = 0; c < 3; c++)
      f[c] /= ln;
    r[0] = f[1] * up[2] - f[2] * up[1];
    r[1] = f[2] * up[0] - f[0] * up[2];
    r[2] = f[0] * up[1] - f[1] * up[0];
    ln = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    if (!(ln > 0.0)) ln = 1.0;
    for (int c = 0; c < 3; c++)
      r[c] /= ln;
    u2[0] = r[1] * f[2] - r[2] * f[1];
    u2[1] = r[2] * f[0] - r[0] * f[2];
    u2[2] = r[0] * f[1] - r[1] * f[0];
    double th = 0.5 * HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0;
    double cs = cos(th), sn = sin(th);
    /* Четыре боковые: нормаль внутрь есть cs·(ось) ± sn·(вперёд) с точностью до
     * знака; ближняя — плоскость через глаз вперёд. */
    double nn[6][3];
    for (int c = 0; c < 3; c++) {
      nn[0][c] = cs * r[c] + sn * f[c];
      nn[1][c] = -cs * r[c] + sn * f[c];
      nn[2][c] = cs * u2[c] + sn * f[c];
      nn[3][c] = -cs * u2[c] + sn * f[c];
      nn[4][c] = f[c];
      nn[5][c] = -f[c];
    }
    for (int k = 0; k < 6; k++) {
      double d = 0.0;
      for (int c = 0; c < 3; c++) {
        X.fr[k][c] = nn[k][c];
        d += nn[k][c] * e[c];
      }
      X.fr[k][3] = -d;
    }
    /* Дальняя плоскость: сцена целиком, то есть отодвинута за диагональ. */
    double diag = 0.0;
    for (int c = 0; c < 3; c++) {
      double s = T.nd[0].hi[c] - T.nd[0].lo[c];
      diag += s * s;
    }
    X.fr[5][3] += sqrt(diag);
    memcpy(camo, e, sizeof camo);
    memcpy(camf, f, sizeof camf);
    /* Пирамида нужна и проходу ПОМЕТОК, и прямой проверке в обходе. Считается
     * она в обоих случаях, а ВКЛЮЧАЕТСЯ в обходе только по `cam=1`: иначе
     * пирамида и пометка сделали бы одну работу дважды. */
    X.usefrustum = cam;
    printf("== ПИРАМИДА КАМЕРЫ: глаз %.2f,%.2f,%.2f, цель %.2f,%.2f,%.2f, поле %.0f°;\n"
           "   ячейка ЦЕЛИКОМ снаружи берётся ОДНОЙ (огрубление), а не выбрасывается\n",
           e[0], e[1], e[2], at[0], at[1], at[2], HZ_CFG_FOV_DEG);
  }
  printf("== ФРОНТ: направление %.3f,%.3f,%.3f, пол %g px ОТ ИСТОЧНИКА, потолок огрубления %d, "
         "заслонение %s\n",
         X.dir[0], X.dir[1], X.dir[2], px, HZ_PFRONT_COARSEN_MAX, shadow ? "ВКЛ" : "ВЫКЛ");
  printf("== ОГОВОРКА ОПЕРАТОРА: %s\n", hz_pfront_note());

  /* ПРОХОД ПОМЕТОК ОТ КАМЕРЫ (ключ `mark=1`). Дёшев: идёт только по крупным
   * узлам. Кладёт каждому предельный уровень, который фронт читает третьим
   * условием остановки. */
  unsigned char *dep = calloc((size_t)T.nnd, 1);
  if (dep == NULL) return 2;
  hz_pocc_depth(&T, dep);
  X.depth = dep;
  unsigned char *mk = NULL;
  if (usemark) {
    /* ПРОХОД ОТ КАМЕРЫ, ищущий НЕВИДИМОЕ НАПРЯМУЮ. Это ТОТ ЖЕ фронт, только
     * пущенный от камеры и остановленный на крупных уровнях: где он приходит
     * тёмным, там от камеры не видно, и там ставится пометка. Направление берётся
     * взглядом камеры — приближение параллельным пучком, точное в середине кадра
     * и грубеющее к краям; конус при этом учтён отдельно, пирамидой. */
    mk = calloc((size_t)T.nnd, 1);
    if (mk == NULL) return 2;
    for (int32_t i = 0; i < T.nnd; i++)
      mk[i] = 255; /* 255 — пометки нет */
    hz_pfront_ctx CX = X;
    CX.mark = NULL;
    CX.markout = mk;
    CX.markoutlev = HZ_PMARK_LEVEL;
    CX.markrelax = relax;
    CX.depth = dep;
    CX.cliplev = HZ_PMARK_LEVEL;
    CX.usefrustum = 1;
    CX.nmarkset = 0;
    for (int c = 0; c < 3; c++)
      CX.dir[c] = camf[c];
    CX.t_entry = 1e300;
    for (int k = 0; k < 8; k++) {
      double tt = 0.0;
      for (int c = 0; c < 3; c++)
        tt += ((k & (1 << c)) ? T.nd[0].hi[c] : T.nd[0].lo[c]) * CX.dir[c];
      if (tt < CX.t_entry) CX.t_entry = tt;
    }
    hz_pfront_face ci[3], co[3];
    for (int a = 0; a < 3; a++) {
      ci[a].c0 = 1.0;
      ci[a].cu = 0.0;
      ci[a].cv = 0.0;
    }
    double tmk = now_s();
    hz_pfront_walk(&CX, 0, T.nd[0].lo, T.nd[0].hi, ci, co, list, m.nt, 0);
    printf("== ПРОХОД ОТ КАМЕРЫ за %.2f с: ячеек %lld, ПОМЕТОК ПОСТАВЛЕНО %lld (невидимое, "
           "огрубление на %d уровня)\n",
           now_s() - tmk, (long long)CX.ncell, (long long)CX.nmarkset, relax);
    X.mark = mk;
    X.markcur = -1;
    X.markrelax = relax;
  }
  hz_pfront_face in[3], out[3];
  for (int a = 0; a < 3; a++) {
    in[a].c0 = 1.0; /* на входе в сцену диск источника открыт целиком */
    in[a].cu = 0.0;
    in[a].cv = 0.0;
  }
  double *refpt = NULL, *refvis = NULL;
  if (nref > 0) {
    refpt = malloc(3 * (size_t)nref * sizeof *refpt);
    refvis = malloc((size_t)nref * sizeof *refvis);
    X.refpt = refpt;
    X.refvis = refvis;
    X.refcap = nref;
    X.refstride = 997; /* простое число: выборка не попадает в такт обхода */
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
  printf("   вне пирамиды с ОТПУЩЕННЫМ полом %lld ячеек\n", (long long)X.noutside);
  printf("   остановок ПО ПОМЕТКЕ %lld\n", (long long)X.nmarkstop);
  printf("   ЯЧЕЕК С ГЕОМЕТРИЕЙ ОДНОЙ ПЛОСКОСТИ %lld из %lld (%.2f %%) — их дробить незачем ни при "
         "каком поле\n",
         (long long)X.nflat1, (long long)X.nflatgeo,
         X.nflatgeo > 0 ? 100.0 * (double)X.nflat1 / (double)X.nflatgeo : 0.0);
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
