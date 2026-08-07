/* psun.c — ШАГ О73 (Ф3): ФРОНТ ОТ СОЛНЦА С СОСТОЯНИЕМ НА ГРАНЯХ.
 * План §266, аудит §267, поправки 08-05 (огрублять, а не отсекать).
 *
 * ЗАПУСК:
 *   psun ФАЙЛ.obj МАСШТАБ [leaf=N] [lev=N] [grade=0|1] [px=9.1]
 *        [sun=dx,dy,dz] [shadow=0] [cam=1] [mark=1|2] [buf=N] [loose=1|2]
 *        [ref=N] [eta=0]
 *   `shadow=0` — НЕГАТИВНЫЙ КОНТРОЛЬ: ни заслонения, ни огрубления по ровности.
 *   Обязан воспроизвести числа О72 ДО ЕДИНИЦЫ; расхождение означает, что
 *   изменилось что-то ещё.
 *   `mark=1` — пометки невидимости буфером дальности (§278, односторонние);
 *   `mark=2` — ПРЕЖНИЙ проход фронтом, приёмку НЕ проходящий (оставлен ради
 *              воспроизводимости §271/§272 и как предмет сравнения);
 *   `buf=N`  — сторона буфера дальности, по умолчанию сторона кадра;
 *   `loose=` — негативные контроли буфера, разбор в `pcull.h`;
 *   `eta=0`  — без эталона лучами по диску: он стоит `8 × 8` лучей на выборочную
 *              ячейку перебором ВСЕХ треугольников и на Сан-Мигеле идёт
 *              десятки минут, тогда как приёмка пометок — один луч на ячейку.
 */
#include "image.h"
#include "pclip.h"
#include "pcull.h"
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
      usemark = 0, relax = HZ_PMARK_OUT_LEVELS, refinemax = 0, nk = 1;
  double camo[3] = {0.0, 0.0, 0.0}, camf[3] = {0.0, 0.0, 1.0}, camat[3] = {0.0, 0.0, 1.0};
  int loose = 0, bufside = HZ_CFG_W, eta = 1, mlev = HZ_PMARK_LEVEL, img = 0;
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
    if (strncmp(argv[i], "loose=", 6) == 0) loose = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "buf=", 4) == 0) bufside = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "eta=", 4) == 0) eta = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "mlev=", 5) == 0) mlev = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "img=", 4) == 0) img = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "cover=", 6) == 0)
      hz_pfront_cover_sum = (int)strtol(argv[i] + 6, NULL, 10);
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
  /* ВЫБОРКИ ПО ДИСКУ ИСТОЧНИКА (Ф4, §275). `K = 1` — в точности О73 и он же
   * негативный контроль. Выборка — спираль Ферма, детерминированная; у эталона
   * она ДРУГАЯ и вчетверо гуще (А533), иначе сверка сойдётся по построению. */
  X.nk = (nk < 1) ? 1 : (nk > HZ_PFRONT_MAXK ? HZ_PFRONT_MAXK : nk);
  {
    double half = 0.5 * 9.3e-3, t1[3] = {0.0, 0.0, 1.0}, q1[3], q2[3];
    if (fabs(X.dir[2]) > 0.9) {
      t1[0] = 1.0;
      t1[2] = 0.0;
    }
    q1[0] = X.dir[1] * t1[2] - X.dir[2] * t1[1];
    q1[1] = X.dir[2] * t1[0] - X.dir[0] * t1[2];
    q1[2] = X.dir[0] * t1[1] - X.dir[1] * t1[0];
    double lq = sqrt(q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2]);
    if (!(lq > 0.0)) lq = 1.0;
    for (int c = 0; c < 3; c++)
      q1[c] /= lq;
    q2[0] = X.dir[1] * q1[2] - X.dir[2] * q1[1];
    q2[1] = X.dir[2] * q1[0] - X.dir[0] * q1[2];
    q2[2] = X.dir[0] * q1[1] - X.dir[1] * q1[0];
    for (int s = 0; s < X.nk; s++) {
      double rr = (X.nk == 1) ? 0.0 : half * sqrt(((double)s + 0.5) / (double)X.nk);
      double ph = 2.39996322972865332 * (double)s;
      double nn = 0.0;
      for (int c = 0; c < 3; c++) {
        X.sdir[s][c] = X.dir[c] + rr * (cos(ph) * q1[c] + sin(ph) * q2[c]);
        nn += X.sdir[s][c] * X.sdir[s][c];
      }
      nn = sqrt(nn);
      if (!(nn > 0.0)) nn = 1.0;
      for (int c = 0; c < 3; c++)
        X.sdir[s][c] /= nn;
    }
  }
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
  /* Камера нужна и картинке (`img`), а не только пометкам: без этого условия
   * прогон с `mark=0` снимал кадр из точки (0,0,0) — поймано сравнением
   * картинок, где различие вышло во весь сигнал. */
  if (cam || usemark || img) {
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
    memcpy(camat, at, sizeof camat);
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
    mk = calloc((size_t)T.nnd, 1);
    if (mk == NULL) return 2;
    for (int32_t i = 0; i < T.nnd; i++)
      mk[i] = 255; /* 255 — пометки нет */
  }
  if (usemark == 1) {
    /* ПОМЕТКИ БУФЕРОМ ДАЛЬНОСТИ ОТ ГЛАЗА (§278, `pcull.h`). Здесь камера —
     * ТОЧКА, а обе границы (заслонитель по дальней точке накрытого ЦЕЛИКОМ
     * пикселя, ячейка по ближней своей точке против максимума по всей проекции)
     * повёрнуты в сторону «пометок меньше». Приёмка обязана дать НОЛЬ. */
    double diag = 0.0;
    for (int c = 0; c < 3; c++) {
      double s = T.nd[0].hi[c] - T.nd[0].lo[c];
      diag += s * s;
    }
    diag = sqrt(diag);
    double upv[3] = HZ_CFG_UP;
    hz_pcull C;
    double tmk = now_s();
    if (hz_pcull_init(&C, camo, camat, upv, HZ_CFG_FOV_DEG, bufside, diag) != 0) {
      fprintf(stderr, "отказ буфера дальности\n");
      return 2;
    }
    C.loose = loose;
    /* Пирамида кадра отдаётся `pcull` ВСЕГДА, а не по ключу `cam`: она отвечает
     * на вопрос «виден ли этот кусок пространства», а не служит послаблением
     * пола. Вклады «вне кадра» и «заслонено» печатаются порознь. */
    memcpy(C.fr, X.fr, sizeof C.fr);
    C.usefr = 1;
    hz_pcull_draw(&C, &m);
    hz_pcull_pyramid(&C);
    int64_t nset = 0;
    hz_pcull_marks(&C, &T, mlev, mk);
    /* ДОЛЯ СЦЕНЫ, ЗАСЛОНЁННАЯ ОТ КАМЕРЫ, — величина про ПРЕДСТАВЛЕНИЕ, а не про
     * фронт (замечание пользователя 08-07). §279 мерил пометки против пола
     * СОЛНЕЧНОГО прохода и получил полтора процента; но это пересечение двух
     * разных множеств, а не размер заслонённого. Здесь считается прямо: сколько
     * ЗАНЯТЫХ ЛИСТЬЕВ — то есть самой мелкой части представления — лежит под
     * заслонёнными узлами. Лист взят мерой сознательно: он и есть то, что
     * предлагается огрублять.
     *
     * ПИРАМИДА КАДРА В ЭТО ЧИСЛО НЕ ВХОДИТ: `pcull` судит невидимость ТОЛЬКО по
     * заслонению одних тел другими, а всё, что за краем кадра, он не помечает
     * вовсе. Значит доля НЕ ЗАВЫШЕНА за счёт того, что просто не попало в кадр. */
    {
      int64_t leaf_all = 0, leaf_hid = 0, occ_all = 0, occ_hid = 0, leaf_occl = 0, leaf_frust = 0;
      int32_t *st = malloc((size_t)T.nnd * sizeof *st);
      unsigned char *hid = malloc((size_t)T.nnd);
      if (st != NULL && hid != NULL) {
        int32_t sp = 0;
        st[sp] = 0;
        hid[0] = (unsigned char)(mk[0] != 255 ? mk[0] : 0);
        sp = 1;
        while (sp > 0) {
          int32_t nid = st[--sp];
          int h = hid[nid];
          if (occ[nid]) {
            occ_all++;
            if (h) occ_hid++;
            if (T.nd[nid].child < 0) {
              leaf_all++;
              if (h) leaf_hid++;
              if (h == 1) leaf_occl++;
              if (h == 2) leaf_frust++;
            }
          }
          if (T.nd[nid].child >= 0)
            for (int k = 0; k < 8; k++) {
              int32_t c = T.nd[nid].child + k;
              /* Причина наследуется вниз: заслонённый кусок пространства не
               * становится видимым от того, что его поделили. `1` — заслонено
               * телами, `2` — вне кадра; при споре побеждает та, что выше. */
              hid[c] = (unsigned char)(h != 0 ? h : (mk[c] != 255 ? mk[c] : 0));
              st[sp++] = c;
            }
        }
        printf(
            "== НЕ ВИДНО ИЗ КАМЕРЫ — ДОЛЯ ПРЕДСТАВЛЕНИЯ (занятый лист есть самое мелкое,\n"
            "   что в дереве есть, и ровно то, что предлагается огрублять):\n"
            "   ВСЕГО %lld из %lld (%.1f %%); из них ЗАСЛОНЕНО ТЕЛАМИ %lld (%.1f %%), вне "
            "кадра %lld (%.1f %%)\n"
            "   занятых узлов %lld из %lld (%.1f %%)\n",
            (long long)leaf_hid, (long long)leaf_all,
            leaf_all > 0 ? 100.0 * (double)leaf_hid / (double)leaf_all : 0.0, (long long)leaf_occl,
            leaf_all > 0 ? 100.0 * (double)leaf_occl / (double)leaf_all : 0.0,
            (long long)leaf_frust,
            leaf_all > 0 ? 100.0 * (double)leaf_frust / (double)leaf_all : 0.0, (long long)occ_hid,
            (long long)occ_all, occ_all > 0 ? 100.0 * (double)occ_hid / (double)occ_all : 0.0);
      }
      free(st);
      free(hid);
    }
    for (int32_t i = 0; i < T.nnd; i++)
      if (mk[i] != 255) nset++;
    printf("== БУФЕР ДАЛЬНОСТИ ОТ ГЛАЗА за %.2f с (%d²%s): пикселей занято %lld, треугольников "
           "мимо кадра %lld, у ближней плоскости %lld\n",
           now_s() - tmk, bufside,
           loose == 1   ? ", НЕГАТИВНЫЙ КОНТРОЛЬ: по СЕРЕДИНЕ пикселя"
           : loose == 2 ? ", НЕГАТИВНЫЙ КОНТРОЛЬ: БЛИЖНЯЯ точка пикселя"
                        : "",
           (long long)C.npix, (long long)C.ntri_off, (long long)C.ntri_near);
    printf("== ПОМЕТКИ: узлов проверено %lld, НЕ ВИДНО %lld (из них ВНЕ КАДРА %lld, ЗАСЛОНЕНО "
           "%lld), пометок %lld; отказов: за краем буфера %lld, у ближней плоскости %lld; "
           "огрубление на %d уровня\n",
           (long long)C.ntest, (long long)(C.nfrust + C.nhidden), (long long)C.nfrust,
           (long long)C.nhidden, (long long)nset, (long long)C.noff, (long long)C.nnear, relax);
    hz_pcull_free(&C);
    X.mark = mk;
    X.markcur = -1;
    X.markrelax = relax;
  } else if (usemark) {
    /* ПРЕЖНИЙ ПРОХОД ФРОНТОМ (`mark=2`) — оставлен ради воспроизводимости чисел
     * §271/§272 и как то, ЧТО ИМЕННО не проходит приёмку. Это ТОТ ЖЕ фронт,
     * пущенный от камеры: где он приходит тёмным, там ставится пометка.
     * Направление берётся взглядом камеры — приближение параллельным пучком. */
    hz_pfront_ctx CX = X;
    CX.nk = 1; /* пометки ставит ОДНА выборка: невидимость двоична */
    for (int c = 0; c < 3; c++)
      CX.sdir[0][c] = camf[c];
    CX.mark = NULL;
    CX.markout = mk;
    CX.markoutlev = mlev;
    CX.markrelax = relax;
    CX.depth = dep;
    CX.cliplev = mlev;
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
    hz_pfront_face ci[3 * HZ_PFRONT_MAXK], co[3 * HZ_PFRONT_MAXK];
    for (int a = 0; a < 3 * CX.nk; a++) {
      ci[a].c0 = 1.0;
      ci[a].cu = 0.0;
      ci[a].cv = 0.0;
    }
    double tmk = now_s();
    hz_pfront_walk(&CX, 0, T.nd[0].lo, T.nd[0].hi, ci, co, list, m.nt, 0, 0);
    printf("== ПРОХОД ОТ КАМЕРЫ за %.2f с: ячеек %lld, ПОМЕТОК ПОСТАВЛЕНО %lld (невидимое, "
           "огрубление на %d уровня)\n",
           now_s() - tmk, (long long)CX.ncell, (long long)CX.nmarkset, relax);
    X.mark = mk;
    X.markcur = -1;
    X.markrelax = relax;
  }
  /* СРЕЗ КАДРА — ЧТО ИМЕННО ОГРУБЛЕНИЕ НЕВИДИМОГО ЭКОНОМИТ (ключ `cut=1`).
   *
   * ЗАЧЕМ ОТДЕЛЬНО ОТ ФРОНТА. §279 мерил пометки против фронта и получил
   * `1.016×`, а §283 показал, что невидимо `56…88 %` представления. Оба числа
   * верны, и противоречия нет: **фронт останавливается КРУПНЕЕ той зернистости,
   * на которой невидимость вообще доказуема.** Пометки живут на уровнях `8…12`,
   * а солнечный фронт при своём поле встаёт на `5…7` — и внутрь помеченного
   * узла просто не спускается («остановок ПО ПОМЕТКЕ `358` из `195 364`»).
   * Значит мерить надо потребителя, работающего НА ЗЕРНИСТОСТИ ЛИСТА, а таков
   * срез кадра: набор ячеек, который кадру нужен от представления.
   *
   * ЧТО СЧИТАЕТСЯ. Спуск с полом камеры `4r² ≤ (px·ε)²·d²`, где `d` — расстояние
   * ОТ ГЛАЗА (а не вдоль солнца, как у фронта). Пустой узел берётся целиком
   * (§241.4). Три варианта: без пометок, с огрублением невидимого на `relax`
   * ступеней и с полным его схлопыванием — последнее есть ПОТОЛОК, и держать
   * его рядом обязательно (§279: потолок мерить до доводки).
   *
   * ОТСЕЧЕНИЯ ЗДЕСЬ НЕТ НИ В ОДНОМ ВАРИАНТЕ: невидимая ячейка БЕРЁТСЯ ЦЕЛИКОМ,
   * а не выбрасывается — свет приходит из-за кадра и от невидимых поверхностей
   * (указание пользователя 08-05, §259 граница 2). */
  if (usemark && mk != NULL) {
    static const char *NM[3] = {"без пометок", "невидимое грубее на relax", "невидимое ЦЕЛИКОМ"};
    double eps2 = (px * HZ_CFG_EPS) * (px * HZ_CFG_EPS);
    int64_t cnt[3] = {0, 0, 0};
    double tcut = now_s();
    for (int mode = 0; mode < 3; mode++) {
      int32_t *st = malloc((size_t)T.nnd * sizeof *st);
      unsigned char *hd = malloc((size_t)T.nnd);
      if (st == NULL || hd == NULL) {
        free(st);
        free(hd);
        break;
      }
      int32_t sp = 0;
      st[sp] = 0;
      hd[0] = (unsigned char)(mk[0] != 255);
      sp = 1;
      while (sp > 0) {
        int32_t nid = st[--sp];
        int h = hd[nid];
        const hz_ptnode *N = &T.nd[nid];
        int stop = 0;
        if (!occ[nid] || N->child < 0) {
          stop = 1;
        } else if (h && mode == 2) {
          stop = 1;
        } else {
          double r2 = 0.0, d2 = 0.0;
          for (int c = 0; c < 3; c++) {
            double hh = 0.5 * (N->hi[c] - N->lo[c]);
            double dc = 0.5 * (N->lo[c] + N->hi[c]) - camo[c];
            r2 += hh * hh;
            d2 += dc * dc;
          }
          double e2 = eps2;
          if (h && mode == 1)
            for (int k = 0; k < relax && k < 16; k++)
              e2 *= 4.0;
          if (d2 > r2 && 4.0 * r2 <= e2 * d2) stop = 1;
        }
        if (stop) {
          cnt[mode]++;
          continue;
        }
        for (int k = 0; k < 8; k++) {
          int32_t c = N->child + k;
          hd[c] = (unsigned char)(h || mk[c] != 255);
          st[sp++] = c;
        }
      }
      free(st);
      free(hd);
    }
    printf("== СРЕЗ КАДРА (пол камеры %g px ОТ ГЛАЗА, ячеек в представлении, "
           "нужном кадру), за %.2f с:\n",
           px, now_s() - tcut);
    for (int mode = 0; mode < 3; mode++)
      printf("   %-28s %10lld\n", NM[mode], (long long)cnt[mode]);
    if (cnt[1] > 0 && cnt[2] > 0)
      printf("   выигрыш: огрубление на %d ступени — %.2fx, ПОТОЛОК — %.2fx\n", relax,
             (double)cnt[0] / (double)cnt[1], (double)cnt[0] / (double)cnt[2]);
  }

  hz_pfront_face in[3 * HZ_PFRONT_MAXK], out[3 * HZ_PFRONT_MAXK];
  for (int a = 0; a < 3 * X.nk; a++) {
    in[a].c0 = 1.0; /* на входе в сцену диск источника открыт целиком */
    in[a].cu = 0.0;
    in[a].cv = 0.0;
  }
  double *refpt = NULL, *refvis = NULL, *refh = NULL;
  if (nref > 0) {
    refpt = malloc(3 * (size_t)nref * sizeof *refpt);
    refvis = malloc((size_t)nref * sizeof *refvis);
    refh = malloc((size_t)nref * sizeof *refh);
    X.refh = refh;
    X.refpt = refpt;
    X.refvis = refvis;
    X.refcap = nref;
    X.refstride = 997; /* простое число: выборка не попадает в такт обхода */
  }
  t0 = now_s();
  X.refinemax = refinemax;
  /* Место, куда фронт положит свой ответ, — только если он кому-то нужен. */
  if (img) {
    X.cellf = malloc((size_t)T.nnd * sizeof *X.cellf);
    if (X.cellf != NULL)
      for (int32_t i = 0; i < T.nnd; i++)
        X.cellf[i] = -1.0f;
  }
  hz_pfront_walk(&X, 0, T.nd[0].lo, T.nd[0].hi, in, out, list, m.nt, 0, 0);
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

  /* ПРОХОД 3 — СБОР ПО ПИКСЕЛЮ, ТО ЕСТЬ КАРТИНКА (ключ `img=1`).
   *
   * ЗАЧЕМ ОН ЗДЕСЬ И ПОЧЕМУ ТОЛЬКО ТЕПЕРЬ. У линии фронта не было ни одной
   * картинки: фронт считал долю открытого диска и ВЫБРАСЫВАЛ её, поэтому всякий
   * выигрыш оставался числом в счётчике. Замечание пользователя §277 («эти
   * оптимизации ничего не позволили нового») упирается ровно в это. Теперь фронт
   * оставляет ответ в `cellf`, и его можно посмотреть глазами.
   *
   * УСТРОЙСТВО, ПРОСТЕЙШЕЕ ИЗ ЗАКОННЫХ (§7: «камера на первое время как угодно»):
   * обычный z-буфер с номером треугольника даёт видимую поверхность в пикселе,
   * точка попадания — из дальности, ячейка — спуском по дереву до той, где фронт
   * остановился, яркость — `ρ·f·|n·ω|/π`.
   *
   * ЛУЧЕЙ ЗДЕСЬ НЕТ И БЫТЬ НЕ МОЖЕТ: дерево строится укладкой резкой и ссылок на
   * треугольники не хранит вовсе (А472), так что `ptrace` по нему не пойдёт.
   * Проекция отвечает на тот же вопрос разом для всех пикселей. */
  if (img) {
    double timg = now_s();
    size_t np = (size_t)bufside * (size_t)bufside;
    float *zb = malloc(np * sizeof *zb);
    int32_t *ib = malloc(np * sizeof *ib);
    double *pix = malloc(np * sizeof *pix);
    /* Доля открытого диска ОТДЕЛЬНО от яркости: эталон спрашивает «видно ли
     * солнце», а не «насколько ярко», и мешать в одну величину косинус с
     * альбедо значило бы сверять три вещи разом. */
    double *fpix = malloc(np * sizeof *fpix);
    if (zb != NULL && ib != NULL && pix != NULL && fpix != NULL) {
      double diag = 0.0;
      for (int c = 0; c < 3; c++) {
        double s = T.nd[0].hi[c] - T.nd[0].lo[c];
        diag += s * s;
      }
      double upv[3] = HZ_CFG_UP;
      hz_pcull C2;
      if (hz_pcull_init(&C2, camo, camat, upv, HZ_CFG_FOV_DEG, bufside, sqrt(diag)) == 0) {
        hz_pcull_shot(&C2, &m, zb, ib);
        int64_t nsky = 0, nlit2 = 0, nsh2 = 0;
        double sum = 0.0;
        for (size_t p = 0; p < np; p++) {
          pix[p] = 0.0;
          fpix[p] = 0.0;
          if (ib[p] < 0) {
            nsky++;
            continue;
          }
          double d[3];
          hz_pcull_ray(&C2, (int)(p % (size_t)bufside), (int)(p / (size_t)bufside), d);
          double P[3];
          for (int c = 0; c < 3; c++)
            P[c] = camo[c] + (double)zb[p] * d[c];
          /* Ячейка, в которой фронт остановился: спуск, пока ответа нет. */
          int32_t nid = 0;
          while (X.cellf[nid] < 0.0f && T.nd[nid].child >= 0) {
            double mid[3];
            int k = 0;
            for (int c = 0; c < 3; c++) {
              mid[c] = 0.5 * (T.nd[nid].lo[c] + T.nd[nid].hi[c]);
              if (P[c] >= mid[c]) k |= (1 << c);
            }
            nid = T.nd[nid].child + k;
          }
          double f = (X.cellf[nid] >= 0.0f) ? (double)X.cellf[nid] : 0.0;
          /* Косинус на поверхности — из нормали видимого треугольника. */
          int32_t t = ib[p];
          const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
          const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
          const double *Cc = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
          double u1[3], u2[3], nn[3];
          for (int c = 0; c < 3; c++) {
            u1[c] = B[c] - A[c];
            u2[c] = Cc[c] - A[c];
          }
          nn[0] = u1[1] * u2[2] - u1[2] * u2[1];
          nn[1] = u1[2] * u2[0] - u1[0] * u2[2];
          nn[2] = u1[0] * u2[1] - u1[1] * u2[0];
          double ln2 = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
          if (!(ln2 > 0.0)) ln2 = 1.0;
          double cs = 0.0;
          for (int c = 0; c < 3; c++)
            cs -= (nn[c] / ln2) * X.dir[c]; /* `dir` — направление ЛУЧЕЙ источника */
          if (cs < 0.0) cs = -cs;           /* нормаль .obj может смотреть внутрь */
          double kd = m.mtl[m.fm[t]].kd;
          pix[p] = kd * f * cs;
          fpix[p] = f;
          sum += pix[p];
          if (f > 0.5)
            nlit2++;
          else
            nsh2++;
        }
        /* ЭТАЛОН ПО ПИКСЕЛЯМ — НАСКОЛЬКО КАРТИНКА НЕВЕРНА. До сих пор у линии
         * фронта сверялись ЯЧЕЙКИ (§272, выборка по обходу), а картинки не было
         * вовсе; теперь спрашивается прямо: в этом пикселе солнце видно или нет?
         * Эталон — теневой луч из точки попадания перебором ВСЕХ треугольников,
         * то есть заведомо проще проверяемого (правило эталона).
         *
         * ЛУЧ ОДИН, А НЕ ШЕСТНАДЦАТЬ, И ЭТО НЕ ЭКОНОМИЯ. При `nk = 1` фронт
         * несёт видимость ЦЕНТРА диска, то есть точечное солнце; сверять его с
         * усреднением по диску значило бы мерить полутень, которой в этой
         * постановке нет (§274). Диск вернётся вместе с `nk > 1`. */
        if (nref > 0) {
          double tref2 = now_s();
          int32_t step = (int32_t)(np / (size_t)(nref > 0 ? nref : 1));
          if (step < 1) step = 1;
          int64_t nch = 0, nlitbad = 0, nshbad = 0;
          double sumd = 0.0, maxd = 0.0;
          for (size_t p = 0; p < np; p += (size_t)step) {
            if (ib[p] < 0) continue;
            double d[3];
            hz_pcull_ray(&C2, (int)(p % (size_t)bufside), (int)(p / (size_t)bufside), d);
            double P[3];
            for (int c = 0; c < 3; c++)
              P[c] = camo[c] + (double)zb[p] * d[c];
            /* Отступ от поверхности вдоль луча НА СОЛНЦЕ, а не по нормали:
             * нормаль у .obj может смотреть внутрь, и отступ по ней уводил бы
             * точку под поверхность ровно в половине случаев. */
            double ofs = 1e-5 * (fabs(P[0]) + fabs(P[1]) + fabs(P[2]) + 1.0);
            double O[3], S[3];
            for (int c = 0; c < 3; c++) {
              S[c] = -X.dir[c];
              O[c] = P[c] + ofs * S[c];
            }
            int blk = 0;
            for (int32_t t = 0; t < m.nt && !blk; t++) {
              const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
              const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
              const double *Cc = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
              double q1[3], q2[3], pv[3], tv[3], qv[3];
              for (int c = 0; c < 3; c++) {
                q1[c] = B[c] - A[c];
                q2[c] = Cc[c] - A[c];
              }
              pv[0] = S[1] * q2[2] - S[2] * q2[1];
              pv[1] = S[2] * q2[0] - S[0] * q2[2];
              pv[2] = S[0] * q2[1] - S[1] * q2[0];
              double det = q1[0] * pv[0] + q1[1] * pv[1] + q1[2] * pv[2];
              if (det > -1e-12 && det < 1e-12) continue;
              double inv = 1.0 / det;
              for (int c = 0; c < 3; c++)
                tv[c] = O[c] - A[c];
              double uu = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
              if (uu < 0.0 || uu > 1.0) continue;
              qv[0] = tv[1] * q1[2] - tv[2] * q1[1];
              qv[1] = tv[2] * q1[0] - tv[0] * q1[2];
              qv[2] = tv[0] * q1[1] - tv[1] * q1[0];
              double vv = (S[0] * qv[0] + S[1] * qv[1] + S[2] * qv[2]) * inv;
              if (vv < 0.0 || uu + vv > 1.0) continue;
              double tt = (q2[0] * qv[0] + q2[1] * qv[1] + q2[2] * qv[2]) * inv;
              if (tt > 0.0) blk = 1;
            }
            double ftrue = blk ? 0.0 : 1.0;
            double four = fpix[p];
            double dd = fabs(four - ftrue);
            nch++;
            sumd += dd;
            if (dd > maxd) maxd = dd;
            if (ftrue < 0.5 && four > 0.5) nlitbad++; /* светло у нас, тень у эталона */
            if (ftrue > 0.5 && four < 0.5) nshbad++;  /* тень у нас, светло у эталона */
          }
          printf("   ЭТАЛОН ПО ПИКСЕЛЯМ (теневой луч перебором, %lld пикселей за %.1f с):\n"
                 "     |Δдоли| среднее %.4f, наибольшее %.4f; РАЗОШЛИСЬ ЗНАКОМ: у нас светло "
                 "а в тени %lld (%.2f %%), у нас тень а на свету %lld (%.2f %%)\n",
                 (long long)nch, now_s() - tref2, nch > 0 ? sumd / (double)nch : 0.0, maxd,
                 (long long)nlitbad, nch > 0 ? 100.0 * (double)nlitbad / (double)nch : 0.0,
                 (long long)nshbad, nch > 0 ? 100.0 * (double)nshbad / (double)nch : 0.0);
        }
        char path[256];
        const char *base = strrchr(argv[1], '/');
        snprintf(path, sizeof path, "img/psun_%s_mark%d.ppm", base ? base + 1 : argv[1], usemark);
        hz_ppm_write(path, pix, bufside, bufside);
        /* PFM рядом с PPM — чтобы сравнивать две картинки ЧИСЛАМИ (`pfmdiff`), а
         * не глазами: приёмка огрубления невидимого есть «видимое не изменилось»,
         * и глаз на такой вопрос не отвечает. */
        snprintf(path, sizeof path, "img/psun_%s_mark%d.pfm", base ? base + 1 : argv[1], usemark);
        hz_pfm_write(path, pix, pix, pix, bufside, bufside);
        printf("== КАРТИНКА (проход 3, сбор по пикселю) за %.2f с: %s\n"
               "   пикселей неба %lld, освещённых %lld, в тени %lld; средняя яркость %.6f\n",
               now_s() - timg, path, (long long)nsky, (long long)nlit2, (long long)nsh2,
               sum / (double)np);
        hz_pcull_free(&C2);
      }
    }
    free(zb);
    free(ib);
    free(pix);
    free(fpix);
  }

  /* ПРИЁМКА ЭТОЙ ЛИНИИ ОПТИМИЗАЦИЙ — ОДНОСТОРОННЯЯ (указание пользователя 08-06):
   * «нужно лишь, чтобы мы не огрубили ВИДИМЫЕ части; вполне допустимо, что часть
   * невидимых останется неогрублённой». Значит проверять надо НЕ расхождение с
   * эталоном по яркости — то физика тени, а у нас тень в смысле НЕВИДИМОСТИ, —
   * а одно: нет ли среди помеченных ячейки, до которой камера ДОСТАЁТ.
   * Проверка независима от фронта: прямой луч от глаза к центру ячейки перебором
   * треугольников. Нарушений обязано быть НОЛЬ. */
  if (usemark && mk != NULL && nref > 0) {
    int32_t nchk = 0, nbadv = 0, nbadbox = 0;
    double tchk = now_s();
    /* РАСПРЕДЕЛЕНИЕ НАРУШЕНИЙ ПО УГЛУ ОТ ОСИ ВЗГЛЯДА — РАЗЛИЧИТЕЛЬ ПРИЧИН.
     * Параллельное приближение камеры ошибается тем сильнее, чем дальше от оси:
     * у края кадра расхождение равно половине поля. Если нарушения сидят у края
     * — виновато оно; если разбросаны по всему полю (и есть у самой оси) —
     * причина другая, и параллельностью её не объяснить. Четыре корзины по
     * долям полуполя. */
    int32_t bin_chk[4] = {0, 0, 0, 0}, bin_bad[4] = {0, 0, 0, 0};
    double halffov = 0.5 * HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0;
    /* Перебирать надо ПОМЕЧЕННЫЕ, а не сетку по всем узлам: пометок 47 тыс. из
     * 34 млн, и шаг по всем узлам не попал ни в одну (проверено 0). Считаем их
     * сперва, потом берём каждую `stride`-ю ИЗ НИХ. */
    int32_t nmk = 0;
    for (int32_t i = 0; i < T.nnd; i++)
      if (mk[i] != 255) nmk++;
    int32_t stride = (nmk > nref) ? nmk / nref : 1;
    int32_t seen = 0;
    for (int32_t i = 0; i < T.nnd; i++) {
      if (mk[i] == 255) continue;
      if ((seen++ % stride) != 0) continue;
      double P0[3], d2[3], len = 0.0;
      for (int c = 0; c < 3; c++) {
        P0[c] = 0.5 * (T.nd[i].lo[c] + T.nd[i].hi[c]);
        d2[c] = P0[c] - camo[c];
        len += d2[c] * d2[c];
      }
      len = sqrt(len);
      if (!(len > 0.0)) continue;
      for (int c = 0; c < 3; c++)
        d2[c] /= len;
      /* ВНУТРИ ли пирамиды — И ЗДЕСЬ ВАЖНО, ЧТО ИМЕННО ПРОВЕРЯЕТСЯ. Луч пускается
       * в ЦЕНТР ячейки, значит и на видимость проверять надо ЭТУ ТОЧКУ, а не
       * коробку: точка вне кадра лучом достижима, но не видима, и нарушением не
       * является. Прежде проверялась КОРОБКА (перекрывает ли она пирамиду), и
       * пара «коробка × луч в центр» была рассогласована.
       *
       * ЧЕСТНОСТЬ ТРЕБУЕТ ПЕЧАТАТЬ ОБА ЧИСЛА, потому что правило пометки в этом
       * же шаге стало отвечать «не видно» и за краем кадра: если бы приёмка
       * смягчалась вместе с правилом, она перестала бы что-либо проверять.
       * Поэтому считаются нарушения ПО ТОЧКЕ (строгая пара) и ПО КОРОБКЕ
       * (прежний, более широкий охват), и печатаются рядом. */
      int inbox = 1, inpt = 1;
      for (int kf = 0; kf < 6; kf++) {
        const double *PL = X.fr[kf];
        double sf = PL[3], sp = PL[3];
        for (int c = 0; c < 3; c++) {
          sf += (PL[c] > 0.0 ? T.nd[i].hi[c] : T.nd[i].lo[c]) * PL[c];
          sp += P0[c] * PL[c];
        }
        if (sf < 0.0) inbox = 0;
        if (sp < 0.0) inpt = 0;
      }
      if (!inbox) continue;
      nchk++;
      int blocked = 0;
      for (int32_t t = 0; t < m.nt && !blocked; t++) {
        const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
        const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
        const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
        double e1[3], e2[3], pv[3], tv[3], qv[3];
        for (int c = 0; c < 3; c++) {
          e1[c] = B[c] - A[c];
          e2[c] = C[c] - A[c];
        }
        pv[0] = d2[1] * e2[2] - d2[2] * e2[1];
        pv[1] = d2[2] * e2[0] - d2[0] * e2[2];
        pv[2] = d2[0] * e2[1] - d2[1] * e2[0];
        double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
        if (det > -1e-12 && det < 1e-12) continue;
        double inv = 1.0 / det;
        for (int c = 0; c < 3; c++)
          tv[c] = camo[c] - A[c];
        double uu = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
        if (uu < 0.0 || uu > 1.0) continue;
        qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
        qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
        qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
        double vv = (d2[0] * qv[0] + d2[1] * qv[1] + d2[2] * qv[2]) * inv;
        if (vv < 0.0 || uu + vv > 1.0) continue;
        double tt = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
        if (tt > 1e-6 && tt < len - 1e-6) blocked = 1;
      }
      double ca = 0.0;
      for (int c = 0; c < 3; c++)
        ca += d2[c] * camf[c];
      if (ca > 1.0) ca = 1.0;
      if (ca < -1.0) ca = -1.0;
      int bin = (int)(4.0 * acos(ca) / halffov);
      if (bin < 0) bin = 0;
      if (bin > 3) bin = 3;
      bin_chk[bin]++;
      if (!blocked) {
        nbadbox++;
        if (inpt) {
          nbadv++;
          bin_bad[bin]++;
        }
      }
    }
    printf("   ПРИЁМКА: помеченных проверено %d за %.1f с; ВИДИМЫХ СРЕДИ НИХ %d (обязано 0); "
           "по прежнему, более широкому критерию коробки — %d\n",
           nchk, now_s() - tchk, nbadv, nbadbox);
    printf("   по углу от оси взгляда (доля полуполя): ");
    for (int b = 0; b < 4; b++)
      printf("%d/4 %d из %d | ", b + 1, bin_bad[b], bin_chk[b]);
    printf("\n");
  }

  /* ЭТАЛОН ЛУЧАМИ (А519). Считает ТО ЖЕ — долю открытого диска источника, — но
   * другим способом: прямым перебором треугольников. Перебор выбран сознательно:
   * эталон обязан быть ПРОЩЕ проверяемого, иначе сверяются две одинаковые
   * ошибки. Отсюда и ограничение по сцене — на зале это 124 тыс. треугольников
   * на луч, и дальше растёт линейно.
   * Диск источника выбирается спиралью Ферма: ДЕТЕРМИНИРОВАННО, без случайности,
   * чтобы повтор давал те же числа. */
  /* ЭТАЛОН ЛУЧАМИ И ПРИЁМКА ПОМЕТОК — РАЗНЫЕ ПРОВЕРКИ, И ЦЕНА У НИХ РАЗНАЯ.
   * Эталон стреляет `8 × 8` лучей на выборочную ячейку перебором ВСЕХ
   * треугольников: на зале это две секунды, на Сан-Мигеле — `600 × 64 × 6.7e6`,
   * то есть десятки минут (замерено 08-06 прогоном, который пришлось прервать).
   * Приёмка невидимости стоит один луч на ячейку и ни от чего этого не зависит.
   * Поэтому эталон включается отдельным ключом `eta=1`. */
  if (eta && nref > 0 && X.refn > 0) {
    double half = 0.5 * 9.3e-3; /* угловой радиус солнца, §241.4 */
    double t1[3] = {0.0, 0.0, 1.0}, r1[3], r2[3];
    if (fabs(X.dir[2]) > 0.9) {
      t1[0] = 1.0;
      t1[2] = 0.0;
    }
    r1[0] = X.dir[1] * t1[2] - X.dir[2] * t1[1];
    r1[1] = X.dir[2] * t1[0] - X.dir[0] * t1[2];
    r1[2] = X.dir[0] * t1[1] - X.dir[1] * t1[0];
    double ln2 = sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
    if (!(ln2 > 0.0)) ln2 = 1.0;
    for (int c = 0; c < 3; c++)
      r1[c] /= ln2;
    r2[0] = X.dir[1] * r1[2] - X.dir[2] * r1[1];
    r2[1] = X.dir[2] * r1[0] - X.dir[0] * r1[2];
    r2[2] = X.dir[0] * r1[1] - X.dir[1] * r1[0];
    /* РАЗДЕЛЕНИЕ ПО СОСТОЯНИЮ (разбор 08-06). Две причины расхождения известны
     * заранее и лечатся по-разному: перекрытие по ОДНОМУ куску (А517) недооценит
     * тень ВЕЗДЕ, а линейность состояния (Р1) промахнётся только на КРАЮ тени.
     * Поэтому ошибка считается порознь в трёх корзинах: глубокая тень, открытый
     * свет и край. Если промахи только на краю — виновата линейность и это цена
     * представления; если и в тени — виновато перекрытие. */
    double sum3[3] = {0.0, 0.0, 0.0}, wor3[3] = {0.0, 0.0, 0.0};
    int64_t cnt3[3] = {0, 0, 0}, bad3[3] = {0, 0, 0};
    double sum = 0.0, worst = 0.0, tref = now_s();
    int64_t nbad = 0;
    for (int32_t s = 0; s < X.refn; s++) {
      /* СВЕРЯТЬ НАДО СРЕДНЕЕ СО СРЕДНИМ (разбор 08-06). Наше состояние есть
       * ЛИНЕЙНАЯ функция по грани, и `c0` — её среднее; эталон же стрелял из
       * ЦЕНТРА, то есть сравнивал точку с проекцией. На краю тени это давало
       * расхождение 0.248 при 100 %% ячеек вне приёмки — артефакт сверки, а не
       * ошибка фронта. Теперь эталон берёт восемь точек по ячейке (углы
       * полуразмера) и усредняет: величина становится той же, что у нас. */
      const double *Pc = refpt + 3 * (size_t)s;
      double hh = (refh != NULL) ? 0.5 * refh[s] : 0.0;
      int open = 0, ntot = 0;
      for (int pt = 0; pt < 8; pt++) {
        double P0[3];
        for (int c = 0; c < 3; c++)
          P0[c] = Pc[c] + ((pt & (1 << c)) ? hh : -hh);
        for (int k = 0; k < HZ_REF_RAYS / 2; k++) {
          ntot++;
          double rr = half * sqrt(((double)k + 0.5) / (double)(HZ_REF_RAYS / 2));
          double ph = 2.39996322972865332 * (double)k;
          double d2[3];
          for (int c = 0; c < 3; c++)
            d2[c] = -X.dir[c] + rr * (cos(ph) * r1[c] + sin(ph) * r2[c]);
          int hit = 0;
          for (int32_t t = 0; t < m.nt && !hit; t++) {
            const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
            const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
            const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
            double e1[3], e2[3], pv[3], tv[3], qv[3];
            for (int c = 0; c < 3; c++) {
              e1[c] = B[c] - A[c];
              e2[c] = C[c] - A[c];
            }
            pv[0] = d2[1] * e2[2] - d2[2] * e2[1];
            pv[1] = d2[2] * e2[0] - d2[0] * e2[2];
            pv[2] = d2[0] * e2[1] - d2[1] * e2[0];
            double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
            if (det > -1e-12 && det < 1e-12) continue;
            double inv = 1.0 / det;
            for (int c = 0; c < 3; c++)
              tv[c] = P0[c] - A[c];
            double uu = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
            if (uu < 0.0 || uu > 1.0) continue;
            qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
            qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
            qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
            double vv = (d2[0] * qv[0] + d2[1] * qv[1] + d2[2] * qv[2]) * inv;
            if (vv < 0.0 || uu + vv > 1.0) continue;
            double tt = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
            if (tt > 1e-6) hit = 1;
          }
          if (!hit) open++;
        }
      }
      double vref = (ntot > 0) ? (double)open / (double)ntot : 0.0;
      double d = refvis[s] - vref;
      if (d < 0.0) d = -d;
      sum += d;
      if (d > worst) worst = d;
      if (d > 0.02) nbad++;
      int b = (refvis[s] < 0.05) ? 0 : ((refvis[s] > 0.95) ? 1 : 2);
      cnt3[b]++;
      sum3[b] += d;
      if (d > wor3[b]) wor3[b] = d;
      if (d > 0.02) bad3[b]++;
    }
    printf("   ЭТАЛОН ЛУЧАМИ (%d лучей по диску, перебор треугольников): выборка %d ячеек за "
           "%.1f с\n",
           HZ_REF_RAYS, X.refn, now_s() - tref);
    printf("   расхождение среднее %.4f, наибольшее %.4f; выше приёмки 2 %% — %lld ячеек "
           "(%.2f %%)\n",
           sum / (double)X.refn, worst, (long long)nbad, 100.0 * (double)nbad / (double)X.refn);
    {
      static const char *nm[3] = {"глубокая тень", "открытый свет", "КРАЙ тени"};
      for (int b = 0; b < 3; b++)
        if (cnt3[b] > 0)
          printf("     %-14s ячеек %4lld, среднее %.4f, наибольшее %.4f, выше 2 %%%% — %lld "
                 "(%.1f %%%%)\n",
                 nm[b], (long long)cnt3[b], sum3[b] / (double)cnt3[b], wor3[b], (long long)bad3[b],
                 100.0 * (double)bad3[b] / (double)cnt3[b]);
    }
  }
  free(refpt);
  free(refh);
  free(refvis);
  /* Пометки и высоты поддеревьев жили до конца обхода; санитайзер 08-06
   * показал, что их не освобождали ВООБЩЕ — при выходе безвредно, но течь на
   * гейте видна, и починить дешевле, чем оговаривать. */
  free(dep);
  free(mk);
  free(X.cellf);
  free(tlo);
  free(thi);
  free(tpl);
  free(list);
  free(occ);
  hz_ptree_free(&T);
  hz_obj_free(&m);
  return 0;
}
