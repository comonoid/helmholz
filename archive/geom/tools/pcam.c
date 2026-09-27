/* pcam — КАМЕРА СЦЕНЫ НАХОДИТСЯ ЗАМЕРОМ, А НЕ НАЗНАЧАЕТСЯ.
 *
 * ЗАЧЕМ. Обе камеры проекта были назначены на глаз и обе оказались негодными,
 * причём выяснилось это не сразу: глаз зала из §2 стоял ВНУТРИ шкафа (всё
 * видимое в 4…160 см), глаз города — ПОД ЗЕМЛЁЙ (доля неба ровно 0.0000). Ни
 * то, ни другое не было заподозрено — оба поймались числами. Третья сцена
 * (San Miguel, 08-04) заводится сразу с замером, чтобы не повторять этот путь
 * в третий раз.
 *
 * КРИТЕРИЙ РАНГА НАЗВАН ДО ПРОГОНА: `p90/p10` по дальностям видимого из точки.
 * То есть побеждает точка, из которой в ОДНОМ кадре видно и близкое, и
 * далёкое. Взят он не за красоту, а потому что это ровно то свойство, ради
 * которого San Miguel и заводился: разброс размеров детали (пункт 6 реестра
 * §170). Ранг по «глубине каньона», которым нашлась камера города (§110),
 * здесь неприменим — Сан-Мигель есть двор, а не улица.
 *
 * ОГРАНИЧЕНИЯ ОТБОРА, тоже до прогона:
 *   до ближайшей геометрии ≥ `PCAM_NEAR_MIN` — иначе глаз внутри предмета,
 *     ровно ошибка зала;
 *   доля неба в `[PCAM_SKY_LO, PCAM_SKY_HI]` — снизу отсекает подвал и
 *     закрытый объём, сверху — чистое поле, где смотреть не на что.
 *
 * КОНТРОЛЬ ЭКВИВАЛЕНТНОСТИ (ключ `check`): тот же счёт в ИЗВЕСТНОЙ точке
 * города §110 обязан воспроизвести уже опубликованные числа — доля неба
 * `0.1501`, горизонт `p10 5.45` м, `p50 8.10` м. Не воспроизвёл — инструмент
 * мерит не то, и найденной камере верить нельзя. Это НЕ негативный контроль, а
 * проверка эквивалентности, и называется так прямо (правило А10).
 *
 * Запуск: `pcam ФАЙЛ.obj МАСШТАБ [ndir=N] [grid=N] [eye=X,Y,Z] [check]`.
 */
#include "pgrid.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Высота глаза над полом. Человеческий рост — величина сцены, а не порог
 * схемы; берётся тот же `1.6` м, что стоит у камеры зала в §2. */
#define PCAM_EYE_H 1.6
/* Ближе этого до геометрии точка считается стоящей ВНУТРИ предмета. Полметра —
 * это половина шага человека, и меньше брать нельзя: глаз зала §2 отстоял от
 * геометрии на `0.04` м и был признан негодным. */
#define PCAM_NEAR_MIN 0.5
/* Доля неба. Снизу — чтобы точка не оказалась в закрытом объёме (у глаза
 * города §2 было ровно `0.0000`); сверху — чтобы не выйти в чистое поле, где
 * разброса дальностей нет по построению. */
/* ЧИСТОЕ ПОЛЕ ВБЛИЗИ. Ранг ОТНОШЕНИЕМ (`p90/p10`, потом `p90/p50`) оказался
 * дефектным дважды по одной причине: любое отношение максимизируется МАЛЫМ
 * знаменателем, то есть награждает точку, прижатую к предмету. Замерено: первый
 * победитель имел `p10 = 0.60` м, второй `p50 = 1.94` м. Лечение не в выборе
 * перцентиля, а в отказе от отношения: ранг есть `p90` (насколько далеко видно)
 * ПРИ УСЛОВИИ, что ближнее поле чисто. Два метра — расстояние вытянутой руки с
 * запасом; меньше брать нельзя, иначе условие перестаёт что-либо отсекать. */
#define PCAM_NEAR_CLEAR 2.0
#define PCAM_SKY_LO 0.02
#define PCAM_SKY_HI 0.50
/* Пол сцены определяется МЕДИАНОЙ высот вершин, как это уже сделано в
 * `pfront`: медиана устойчива к одиночному выбросу, среднее — нет. */
#define PCAM_FLOOR_STEP 97

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

/* Направления — узлы Фибоначчиевой сферы. Детерминированно (то есть прогон
 * воспроизводится побитово) и без сгущения у полюсов, которое даёт наивная
 * сетка по углам. */
static void dir_fib(int i, int n, double d[3]) {
  double ga = M_PI * (3.0 - sqrt(5.0));
  double y = 1.0 - 2.0 * ((double)i + 0.5) / (double)n;
  double r = sqrt(y * y < 1.0 ? 1.0 - y * y : 0.0);
  double th = ga * (double)i;
  d[0] = r * cos(th);
  d[1] = y;
  d[2] = r * sin(th);
}

/* Замер в одной точке. ТРИ РАЗНЫЕ ВЕЛИЧИНЫ, и путать их нельзя — на этом
 * инструмент уже провалил свой же контроль эквивалентности:
 *
 *   ДОЛЯ НЕБА — КОСИНУСНО ВЗВЕШЕННАЯ доля ВЕРХНЕЙ полусферы, как её считает
 *     полукуб в §110 (`Σ dff` по пустым пикселям). Равномерная доля по ПОЛНОЙ
 *     сфере — другая величина: на городской камере она дала `0.0430` против
 *     опубликованных `0.1501`, потому что половина направлений смотрит в землю,
 *     а косинус их и так почти гасит. Здесь косинусное распределение берётся
 *     выборкой (`r = √u`), и тогда простая доля ушедших И ЕСТЬ `Σ dff`.
 *   ГОРИЗОНТ `p10/p50` — только СТРОГО ГОРИЗОНТАЛЬНЫЕ направления, как нижняя
 *     строка боковых граней полукуба. Никакого порога «что считать горизонтом»
 *     не вводится: азимуты берутся при `y = 0` точно.
 *   РАНГ `p90/p10` — по ПОЛНОЙ сфере, как и было объявлено до прогона. Он
 *     остаётся тем, чем назван, и правка определений его не трогает. */
static void probe(const hz_pgrid *g, const hz_objmesh *m, const double o[3], int ndir, double far,
                  double *dist, double *nearest, double *p10, double *p50, double *p90, double *sky,
                  double *h10, double *h50) {
  /* --- полная сфера: ранг и ближайшее --- */
  int nh = 0;
  for (int i = 0; i < ndir; i++) {
    double d[3] = {0.0, 0.0, 0.0};
    dir_fib(i, ndir, d);
    double t = 0.0;
    if (hz_pgrid_trace(g, m, o, d, 0.0, far, NULL, 0, 0, &t)) dist[nh++] = t;
  }
  if (nh == 0) {
    *nearest = far;
    *p10 = *p50 = *p90 = far;
  } else {
    qsort(dist, (size_t)nh, sizeof *dist, cmp_d);
    *nearest = dist[0];
    *p10 = dist[(nh * 10) / 100];
    *p50 = dist[nh / 2];
    *p90 = dist[(nh * 90) / 100 < nh ? (nh * 90) / 100 : nh - 1];
  }

  /* --- верхняя полусфера, косинусное распределение: доля неба --- */
  {
    int nsky = 0;
    for (int i = 0; i < ndir; i++) {
      double u1 = ((double)i + 0.5) / (double)ndir;
      /* Тот же радикальный обратный по основанию 2, что в `pfront`:
       * детерминированно и без полос. */
      uint32_t b = (uint32_t)i;
      b = (b << 16) | (b >> 16);
      b = ((b & 0x55555555u) << 1) | ((b & 0xAAAAAAAAu) >> 1);
      b = ((b & 0x33333333u) << 2) | ((b & 0xCCCCCCCCu) >> 2);
      b = ((b & 0x0F0F0F0Fu) << 4) | ((b & 0xF0F0F0F0u) >> 4);
      b = ((b & 0x00FF00FFu) << 8) | ((b & 0xFF00FF00u) >> 8);
      double u2 = (double)b * 2.3283064365386963e-10;
      double r = sqrt(u1), ph = 2.0 * M_PI * u2;
      double dy = sqrt(1.0 - u1 > 0.0 ? 1.0 - u1 : 0.0);
      double d[3] = {r * cos(ph), dy, r * sin(ph)};
      double t = 0.0;
      if (!hz_pgrid_trace(g, m, o, d, 0.0, far, NULL, 0, 0, &t)) nsky++;
    }
    *sky = (double)nsky / (double)ndir;
  }

  /* --- строго горизонтальные направления: горизонт --- */
  {
    int nz = 0;
    for (int i = 0; i < ndir; i++) {
      double a = 2.0 * M_PI * ((double)i + 0.5) / (double)ndir;
      double d[3] = {cos(a), 0.0, sin(a)};
      double t = 0.0;
      dist[nz++] = hz_pgrid_trace(g, m, o, d, 0.0, far, NULL, 0, 0, &t) ? t : far;
    }
    qsort(dist, (size_t)nz, sizeof *dist, cmp_d);
    *h10 = dist[(nz * 10) / 100];
    *h50 = dist[nz / 2];
  }
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "pcam ФАЙЛ.obj МАСШТАБ [ndir=N] [grid=N] [eye=X,Y,Z] [check]\n");
    return 1;
  }
  int ndir = 512, ngrid = 24, docheck = 0, haseye = 0;
  double eye0[3] = {0.0, 0.0, 0.0};
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "ndir=", 5) == 0) ndir = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "grid=", 5) == 0) ngrid = (int)strtol(argv[i] + 5, NULL, 10);
    if (strcmp(argv[i], "check") == 0) docheck = 1;
    if (strncmp(argv[i], "eye=", 4) == 0) {
      char *e = NULL;
      eye0[0] = strtod(argv[i] + 4, &e);
      if (e != NULL && *e == ',') eye0[1] = strtod(e + 1, &e);
      if (e != NULL && *e == ',') eye0[2] = strtod(e + 1, NULL);
      haseye = 1;
    }
  }
  if (ndir < 32) ndir = 32;
  if (ngrid < 2) ngrid = 2;

  hz_objmesh m;
  double t0 = now_s();
  if (hz_obj_load(&m, argv[1], atof(argv[2])) != 0) {
    fprintf(stderr, "нет сцены %s\n", argv[1]);
    return 1;
  }
  printf("== СЦЕНА %s: треугольников %d, габарит %.2f x %.2f x %.2f м, чтение %.1f с\n", argv[1],
         m.nt, m.hi[0] - m.lo[0], m.hi[1] - m.lo[1], m.hi[2] - m.lo[2], now_s() - t0);

  hz_pgrid g;
  t0 = now_s();
  if (hz_pgrid_build(&g, &m) != 0) {
    fprintf(stderr, "отказ сетки\n");
    hz_pgrid_free(&g);
    return 2;
  }
  printf("== СЕТКА за %.1f с\n", now_s() - t0);

  double diag = 0.0;
  for (int c = 0; c < 3; c++) {
    double d2 = m.hi[c] - m.lo[c];
    diag += d2 * d2;
  }
  diag = sqrt(diag);
  double far = 4.0 * diag;

  /* ПОЛ БЕРЁТСЯ ЛУЧОМ ВНИЗ В КАЖДОЙ ТОЧКЕ, А НЕ МЕДИАНОЙ ВЫСОТ ВЕРШИН.
   * Медиана вершин — НЕ земля: её тянут вверх крыши и стены. Замерено на
   * Rungholt: медиана даёт `13.00` м, тогда как §110 нашёл землю на `4.00` м,
   * то есть ошибка в девять метров, и глаз оказался бы на уровне крыш. Луч
   * вниз из-под потолка габарита даёт МЕСТНЫЙ пол, что и нужно: во дворе с
   * галереями он разный на первом и втором этаже. */
  printf("== ПОЛ: лучом вниз в каждой точке (медиана вершин негодна: на Rungholt она даёт 13.00 м "
         "против замеренных §110 4.00 м)\n");

  double *dist = malloc((size_t)ndir * sizeof *dist);
  if (dist == NULL) {
    hz_pgrid_free(&g);
    return 2;
  }

  /* --- КОНТРОЛЬ ЭКВИВАЛЕНТНОСТИ ЛИБО ОДНА ТОЧКА ------------------------- */
  if (haseye || docheck) {
    double o[3] = {eye0[0], eye0[1], eye0[2]};
    if (docheck && !haseye) {
      double e2[3] = HZ_CFG_CITY_EYE;
      for (int c = 0; c < 3; c++)
        o[c] = e2[c];
    }
    double nr = 0, p10 = 0, p50 = 0, p90 = 0, sky = 0, h10 = 0, h50 = 0;
    probe(&g, &m, o, ndir, far, dist, &nr, &p10, &p50, &p90, &sky, &h10, &h50);
    printf("== ТОЧКА (%.2f, %.2f, %.2f): доля неба %.4f; дальности min %.2f, p10 %.2f, p50 %.2f, "
           "p90 %.2f м; ранг p90/p10 = %.2f; ГОРИЗОНТ p10 %.2f, p50 %.2f м\n",
           o[0], o[1], o[2], sky, nr, p10, p50, p90, (p10 > 0.0) ? p90 / p10 : -1.0, h10, h50);
    if (docheck)
      printf("   КОНТРОЛЬ ЭКВИВАЛЕНТНОСТИ (§110: небо 0.1501, p10 5.45, p50 8.10): %s\n",
             (fabs(sky - 0.1501) < 0.02 && fabs(h10 - 5.45) < 1.0 && fabs(h50 - 8.10) < 1.0)
                 ? "СОШЛОСЬ"
                 : "НЕ СОШЛОСЬ — инструмент мерит не то");
    free(dist);
    hz_pgrid_free(&g);
    return 0;
  }

  /* --- СКАН ПО СЕТКЕ ----------------------------------------------------- */
  printf("== СКАН: %d x %d точек, %d направлений; отбор: ближайшее ≥ %.2f м, небо в %.2f…%.2f\n",
         ngrid, ngrid, ndir, PCAM_NEAR_MIN, PCAM_SKY_LO, PCAM_SKY_HI);
  double best = -1.0, beye[3] = {0.0, 0.0, 0.0};
  double bnr = 0, b10 = 0, b50 = 0, b90 = 0, bsky = 0;
  int64_t nok = 0, nnear = 0, nsky_out = 0;
  double t_scan = now_s();
  for (int ix = 0; ix < ngrid; ix++) {
    for (int iz = 0; iz < ngrid; iz++) {
      double o[3];
      o[0] = m.lo[0] + (m.hi[0] - m.lo[0]) * ((double)ix + 0.5) / (double)ngrid;
      o[2] = m.lo[2] + (m.hi[2] - m.lo[2]) * ((double)iz + 0.5) / (double)ngrid;
      /* Местный пол: луч вниз из-под потолка габарита. Не нашёл — под этой
       * точкой пола нет вовсе, и стоять там нельзя. */
      {
        double top[3] = {o[0], m.hi[1] - 1e-3, o[2]};
        double dn[3] = {0.0, -1.0, 0.0};
        double t = 0.0;
        if (!hz_pgrid_trace(&g, &m, top, dn, 0.0, far, NULL, 0, 0, &t)) {
          nnear++;
          continue;
        }
        o[1] = top[1] - t + PCAM_EYE_H;
      }
      double nr = 0, p10 = 0, p50 = 0, p90 = 0, sky = 0, h10 = 0, h50 = 0;
      probe(&g, &m, o, ndir, far, dist, &nr, &p10, &p50, &p90, &sky, &h10, &h50);
      if (nr < PCAM_NEAR_MIN) {
        nnear++;
        continue;
      }
      if (sky < PCAM_SKY_LO || sky > PCAM_SKY_HI) {
        nsky_out++;
        continue;
      }
      nok++;
      /* РАНГ ИСПРАВЛЕН: `p90/p50`, а не `p90/p10`. Объявленный до прогона
       * `p90/p10` оказался ДЕФЕКТНЫМ, и это видно числом, а не на вкус: он
       * растёт, когда `p10` мал, то есть награждает точку, стоящую ВПЛОТНУЮ к
       * предмету. У победителя первого прогона было `p10 = 0.60` м — в десятой
       * части направлений что-то в шестидесяти сантиметрах, то есть глаз в
       * углу, а вовсе не «видно близкое и далёкое». Медиана одиночным близким
       * предметом не сдвигается. Замена сделана ПО ДЕФЕКТУ МЕТРИКИ, а не по
       * тому, какая точка победила. */
      if (p10 < PCAM_NEAR_CLEAR) {
        nnear++;
        continue;
      }
      double rank = p90;
      if (rank > best) {
        best = rank;
        for (int c = 0; c < 3; c++)
          beye[c] = o[c];
        bnr = nr;
        b10 = p10;
        b50 = p50;
        b90 = p90;
        bsky = sky;
      }
    }
  }
  printf("== ОТБОР: годных точек %lld из %d; отброшено внутри геометрии %lld, по небу %lld; "
         "скан %.1f с\n",
         (long long)nok, ngrid * ngrid, (long long)nnear, (long long)nsky_out, now_s() - t_scan);
  if (best < 0.0) {
    printf("== КАМЕРЫ НЕ НАЙДЕНО: ни одна точка сетки не прошла отбор. Сетка груба либо сцена "
           "закрыта — разбирать вход, а не ослаблять отбор.\n");
    free(dist);
    hz_pgrid_free(&g);
    return 3;
  }
  printf("== КАМЕРА (максимум p90 = %.2f м при чистом ближнем поле p10 ≥ %.1f м):\n", best,
         PCAM_NEAR_CLEAR);
  printf("   ГЛАЗ  {%.2f, %.2f, %.2f}\n", beye[0], beye[1], beye[2]);
  printf("   доля неба %.4f; дальности min %.2f, p10 %.2f, p50 %.2f, p90 %.2f м\n", bsky, bnr, b10,
         b50, b90);

  /* ЦЕЛЬ — В СТОРОНУ НАИБОЛЬШЕЙ ГЛУБИНЫ, а не в центр габарита: центр может
   * лежать внутри геометрии, и тогда камера смотрит в стену. Берётся то из
   * `ndir` направлений, что дало наибольшую дальность при неотрицательной
   * высоте (смотреть в пол незачем). */
  double bat[3] = {0.0, 0.0, 0.0};
  {
    double bestt = -1.0, bd[3] = {1.0, 0.0, 0.0};
    for (int i = 0; i < ndir; i++) {
      double d[3] = {0.0, 0.0, 0.0};
      dir_fib(i, ndir, d);
      if (d[1] < -0.1) continue;
      double t = 0.0;
      if (!hz_pgrid_trace(&g, &m, beye, d, 0.0, far, NULL, 0, 0, &t)) continue;
      if (t > bestt) {
        bestt = t;
        for (int c = 0; c < 3; c++)
          bd[c] = d[c];
      }
    }
    for (int c = 0; c < 3; c++)
      bat[c] = beye[c] + bd[c] * bestt * 0.5;
    printf("   ЦЕЛЬ  {%.2f, %.2f, %.2f}  (направление наибольшей глубины, %.2f м)\n", bat[0],
           bat[1], bat[2], bestt);
  }

  /* СНИМОК ИЗ НАЙДЕННОЙ ТОЧКИ. Камера, принятая по числам без взгляда, — это
   * ровно та ошибка, ради которой инструмент писался: у зала числа 'min 0.04'
   * тоже были числами. Снимок диагностический: яркость по `|n·взгляд|`,
   * притемнение по дальности; освещения тут нет и быть не должно. */
  {
    int W = 512, H = 512;
    unsigned char *rgb = malloc((size_t)W * (size_t)H * 3);
    if (rgb != NULL) {
      double fwd[3], right[3], upv[3], up[3] = HZ_CFG_UP;
      for (int c = 0; c < 3; c++)
        fwd[c] = bat[c] - beye[c];
      double fl = sqrt(fwd[0] * fwd[0] + fwd[1] * fwd[1] + fwd[2] * fwd[2]);
      for (int c = 0; c < 3; c++)
        fwd[c] /= fl;
      right[0] = fwd[1] * up[2] - fwd[2] * up[1];
      right[1] = fwd[2] * up[0] - fwd[0] * up[2];
      right[2] = fwd[0] * up[1] - fwd[1] * up[0];
      double rl = sqrt(right[0] * right[0] + right[1] * right[1] + right[2] * right[2]);
      for (int c = 0; c < 3; c++)
        right[c] /= rl;
      upv[0] = right[1] * fwd[2] - right[2] * fwd[1];
      upv[1] = right[2] * fwd[0] - right[0] * fwd[2];
      upv[2] = right[0] * fwd[1] - right[1] * fwd[0];
      double th2 = tan(0.5 * HZ_CFG_FOV_DEG * M_PI / 180.0);
      for (int j = 0; j < H; j++)
        for (int i = 0; i < W; i++) {
          double sx = (2.0 * ((double)i + 0.5) / W - 1.0) * th2;
          double sy = (1.0 - 2.0 * ((double)j + 0.5) / H) * th2;
          double d[3] = {0.0, 0.0, 0.0};
          for (int c = 0; c < 3; c++)
            d[c] = fwd[c] + sx * right[c] + sy * upv[c];
          double dl = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
          for (int c = 0; c < 3; c++)
            d[c] /= dl;
          double t = 0.0;
          int32_t tri = -1;
          int64_t q = (int64_t)j * W + i;
          double v = 0.55;
          if (hz_pgrid_trace_tri(&g, &m, beye, d, 0.0, far, NULL, 0, 0, &t, &tri) && tri >= 0) {
            const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)tri];
            const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)tri + 1];
            const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)tri + 2];
            double e1[3], e2[3], n[3];
            for (int c = 0; c < 3; c++) {
              e1[c] = B[c] - A[c];
              e2[c] = C[c] - A[c];
            }
            n[0] = e1[1] * e2[2] - e1[2] * e2[1];
            n[1] = e1[2] * e2[0] - e1[0] * e2[2];
            n[2] = e1[0] * e2[1] - e1[1] * e2[0];
            double nl = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
            double cs = (nl > 0.0) ? fabs((n[0] * d[0] + n[1] * d[1] + n[2] * d[2]) / nl) : 0.0;
            v = (0.15 + 0.85 * cs) * exp(-t / (4.0 * b90));
          }
          unsigned char u =
              (unsigned char)(255.0 * pow(v < 0.0 ? 0.0 : (v > 1.0 ? 1.0 : v), 1.0 / 2.2) + 0.5);
          rgb[3 * q] = u;
          rgb[3 * q + 1] = u;
          rgb[3 * q + 2] = u;
        }
      FILE *f = fopen("img/pcam_scene.ppm", "wb");
      if (f != NULL) {
        fprintf(f, "P6\n%d %d\n255\n", W, H);
        fwrite(rgb, 1, (size_t)W * (size_t)H * 3, f);
        fclose(f);
        printf("   СНИМОК img/pcam_scene.ppm (диагностический: наклон грани и дальность, света "
               "нет)\n");
      }
      free(rgb);
    }
  }

  free(dist);
  hz_pgrid_free(&g);
  return 0;
}
