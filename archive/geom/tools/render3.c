/* render3.c — ПОЛНОЦЕННЫЙ РЕНДЕР линии переноса: сфера в освещённой комнате,
 * глобальное освещение развёрткой по ординатам.
 *
 * PLAN_TRANSPORT.md, все три прохода вместе:
 *   развёртка   — полное решение уравнения переноса с многократными отражениями;
 *   поверхности — граничное условие, исходящий радианс хранится DG1 по положению;
 *   сбор        — луч на пиксель до первого попадания, чтение хранимого В ТОЧКЕ.
 *
 * КАРТА ПОЗИЦИОННЫХ КЛЮЧЕЙ (argv, порядок жёсткий):
 *   1 W  2 H  3 nmu  4 out.ppm  5 spec  6 maxbounce  7 nframe  8 facet_only
 *   9 log2n  10 limiter  11 trace  12 maxit  13 warmdelta  14 krylov  15 stol
 *   16 albs  17 nch  18 chsame  19 chdark  20 chdbl  21 lodeps(LOD ε)  22 cutfine
 *   23 raster  24 fsub  25 xwarm(1=П5 перенос, 2=орбита+стабильность, 3=орбита без
 *   гистерезис-состояния)  26 bprof(затухание по отскокам)  27 visstop(видимый
 *   порог reltol=HZ_SWEEP_VISIBLE_RELTOL)  28 twophase(грубый хвост; отклонён
 *   при видимом пороге, прибор ниши альбедо→1).
 *   Примеры: кадр по умолчанию — 640 480 2 img/room.ppm; LOD+видимый порог —
 *   640 480 2 img/v.ppm 0 8 1 0 5 1 0 4000 0 0 1e-9 1 1 0 -1 -1 0.15 0 0 2 0 0 1.
 *
 * СИЛУЭТ СФЕРЫ ТОЧЕН: попадание считается по АНАЛИТИЧЕСКОМУ ПРИМИТИВУ (К15),
 * фасеты работают только пространственным индексом и режут ячейку. Фасетной
 * огранки на картинке поэтому нет вовсе — она была бы, если бы луч спрашивал
 * фасеты.
 *
 * ДВОЙНОГО УЧЁТА НЕТ ПО ПОСТРОЕНИЮ (К9/К17): источник здесь — ИЗЛУЧАЮЩАЯ
 * ПОВЕРХНОСТЬ, развёртка считает всё, а сбор только ЧИТАЕТ. Прямой свет
 * отдельным проходом не добавляется, потому что он уже внутри.
 */

#include "cut/surf.h"
#include "image.h"
#include "octree.h"
#include "transport/cam3.h"
#include "transport/cut3.h"
#include "transport/dirs3.h"
#include "transport/gather3.h"
#include "transport/krylov3.h"
#include "transport/mesh3.h"
#include "transport/raster3.h"
#include "transport/ray3.h"
#include "transport/xfer3.h"
#include "transport/sweep3.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int orbit_fail(int line) {
  fprintf(stderr, "orbit: fail at line %d\n", line);
  return 1;
}

static double now(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* РАЗМЕР КОМНАТЫ В МИРЕ — ФИКСИРОВАН, а число ячеек параметр (К68).
 * Свип по размеру элемента обязан менять ТОЛЬКО сетку: если менять `log2n` не
 * трогая кадр, вырастет сама комната, и сравнивались бы РАЗНЫЕ сцены. Поэтому
 * единица дерева есть `ROOM / nc`, и мир остаётся тем же при любом `log2n`.
 * Уровень фасетизации сфер тоже держится постоянным — иначе за один свип
 * менялись бы две вещи сразу, и это ровно ошибка К53. */
#define ROOM 16.0

#define ROOM 16.0

/* УЧЁТ ХОЛОДНОГО СТАРТА ПО ЭТАПАМ, И ЭТО ТРЕБОВАНИЕ, А НЕ УДОБСТВО.
 * Прежде «холодным стартом» звалась ОДНА развёртка, а вся расстановка сцены —
 * дерево, фасетизация, боковая таблица, сетка граней, матрицы масс разреза —
 * не мерилась нигде и потому молча выпадала из отчёта. Правило CLAUDE.md
 * («докладывать холодный старт и кадр РАЗДЕЛЬНО») этим нарушалось в свою
 * пользу: цифра выходила меньше настоящей.
 *
 * Столбец «на кадр» отвечает на вопрос, ради которого учёт и заведён: ЧТО ИЗ
 * ЭТОГО ПОВТОРИТСЯ, ЕСЛИ ДВИНУТЬ КАМЕРУ. Сегодня — только сбор, потому что ни
 * одна структура выше камеры не читает (`tr3_problem` её не содержит вовсе).
 * ОГОВОРКА, БЕЗ КОТОРОЙ ЭТО ХВАСТОВСТВО: так выходит ещё и потому, что LOD НЕ
 * СДЕЛАН — сетка равномерная. С камерно-привязанным LOD движение камеры меняет
 * набор элементов, и часть расстановки вернётся в кадр. */
typedef struct {
  const char *name;
  double t;
  int per_frame; /* 1 — повторяется при движении камеры */
} stage;

/* Запись боковой таблицы: ячейка и её веер фасетов. Вынесена на уровень файла
 * ради компаратора `qsort` — вставками сортировка была `O(n^2)` и на мелкой
 * сетке съедала десятки минут ПЕРЕД развёрткой. */
typedef struct {
  int32_t cell, f[HZ_P3_MAXH], nf;
} rec_t;

/* захват тактов для bprof (К102-класс инструмента; λ-функций в C нет) */
static double *g_bp_sn;
static int g_bp_cap, g_bp_n;
static void bp_capture(void *ctx, int pass, const double *phi) {
  (void)ctx;
  (void)pass;
  if (g_bp_n < 12 && g_bp_sn != NULL) {
    memcpy(g_bp_sn + (size_t)g_bp_n * (size_t)g_bp_cap, phi, (size_t)g_bp_cap * sizeof(double));
    g_bp_n++; /* растёт ТОЛЬКО при записи: nsnaps == число валидных тактов */
  }
}

static int rec_cmp(const void *a, const void *b) {
  const rec_t *x = a, *y = b;
  return x->cell < y->cell ? -1 : (x->cell > y->cell ? 1 : 0);
}

/* ЭТАП A, ШАГ 1: ДЕРЕВО ПО ПРАВИЛУ `L = εR` СО СХЛОПЫВАНИЕМ ПУСТОТЫ.
 *
 * Разделять этап A надвое стало можно после того, как К52 закрыта замером:
 * при схлопнутой пустоте число ячеек ВДОЛЬ ЛУЧА почти не растёт (`18 → 21 → 22`
 * при учетверении `ε`), тогда как «`23 → 7500`» принадлежало ОБЪЁМНОМУ закону.
 * Значит LOD можно ставить при СТАРОМ, маршевом сборе — и тогда сбор остаётся
 * эталоном, а меняется ровно одна вещь.
 *
 * Дробится узел, только если в нём ЕСТЬ поверхность (стенка куба либо шар) и он
 * крупнее `εR` в БЛИЖНЕЙ своей точке. Ближняя, а не центр: правило есть потолок
 * на размер элемента, и нарушать его хоть где-то внутри узла нельзя. */
/* Возврат: 0 — поверхности нет, 1 — стенка куба, 2 — ФАСЕТЫ тела.
 *
 * Тело спрашивается ТЕМ ЖЕ `hz_facets_for_box`, которым потом строится боковая
 * таблица (К95). Спрашивать аналитический ПРИМИТИВ нельзя: фасеты торчат наружу
 * него на `dmax`, и тогда критерий дробления и таблица говорят о РАЗНЫХ телах —
 * коробка чуть снаружи сферы не дробится, но фасеты получает. */
static int box_has_surface(const int lo[3], int size, int world, const hz_facettab *ft,
                           const int32_t *f0, const int32_t *nfac, int nb) {
  int wall = 0;
  for (int a = 0; a < 3; a++)
    if (lo[a] == 0 || lo[a] + size == world) wall = 1;
  int32_t blo[3] = {lo[0], lo[1], lo[2]};
  int32_t bhi[3] = {lo[0] + size, lo[1] + size, lo[2] + size};
  int32_t sel[HZ_P3_MAXH];
  for (int b = 0; b < nb; b++)
    if (hz_facets_for_box(ft, f0[b], nfac[b], blo, bhi, sel, HZ_P3_MAXH) != HZ_BOX_OUTSIDE)
      return 2;
  return wall;
}

static void lod_build(hz_octree *t, const int lo[3], int size, int world, double u,
                      const double eye[3], double eps, const hz_facettab *ft, const int32_t *f0,
                      const int32_t *nfac, int nb, int cutfine) {
  int hi[3] = {lo[0] + size, lo[1] + size, lo[2] + size};
  int kind = box_has_surface(lo, size, world, ft, f0, nfac, nb);
  if (kind == 0) {
    hz_oct_set_box(t, lo, hi, 1.0);
    return;
  }
  /* РАЗДЕЛИТЕЛЬ К94: при `cutfine` РАЗРЕЗАННЫЕ ячейки дробятся до предела. */
  if (kind == 2 && cutfine && size > 1) {
    int h3 = size / 2;
    for (int k = 0; k < 8; k++) {
      int c[3] = {lo[0] + ((k & 1) ? h3 : 0), lo[1] + ((k & 2) ? h3 : 0),
                  lo[2] + ((k & 4) ? h3 : 0)};
      lod_build(t, c, h3, world, u, eye, eps, ft, f0, nfac, nb, cutfine);
    }
    return;
  }
  double near2 = 0.0;
  for (int a = 0; a < 3; a++) {
    double l = (double)lo[a] * u, h = (double)hi[a] * u, dd = 0.0;
    if (eye[a] < l)
      dd = l - eye[a];
    else if (eye[a] > h)
      dd = eye[a] - h;
    near2 += dd * dd;
  }
  if (size <= 1) {
    hz_oct_set_box(t, lo, hi, 1.0);
    return;
  }
  if (size * u <= eps * sqrt(near2)) {
    hz_oct_set_box(t, lo, hi, 1.0);
    return;
  }
  int h2 = size / 2;
  for (int k = 0; k < 8; k++) {
    int c[3] = {lo[0] + ((k & 1) ? h2 : 0), lo[1] + ((k & 2) ? h2 : 0), lo[2] + ((k & 4) ? h2 : 0)};
    lod_build(t, c, h2, world, u, eye, eps, ft, f0, nfac, nb, cutfine);
  }
}

/* ЭТАП A, ШАГ 3: МИР = ВСЁ, ЧТО ЗАВИСИТ ОТ СЕТКИ. Дерево, боковая таблица,
 * сетка граней, маска полных ячеек и разрез строятся ОДНОЙ функцией, потому что
 * при якоре на камере их придётся строить БОЛЬШЕ ОДНОГО РАЗА: сдвиг камеры
 * меняет `L = εR`, значит меняет дерево, а с ним всё над ним. Выделение —
 * чистая перестановка кода: на той же геометрии тот же вызов обязан дать те же
 * биты (проверено прогоном до/после рефакторинга, см. доклад шага 3). */
typedef struct {
  hz_octree t;
  hz_cutmap cmap;
  tr3_mesh mesh;
  uint8_t *solid;
  tr3_cut cut;
  int nboth, nover; /* ячейки с несколькими телами и из них не влезшие */
} world;

static void world_free(world *w) {
  tr3_cut_free(&w->cut);
  tr3_mesh_free(&w->mesh);
  hz_cutmap_free(&w->cmap);
  free(w->solid);
  hz_oct_free(&w->t);
}

static int world_build(world *w, int log2n, int nc, const hz_frame *fr, const hz_surftab *stab,
                       const hz_facettab *ftab, const int32_t *f0, const int32_t *nfac, int nb,
                       double lodeps, const double eye[3], int cutfine) {
  memset(w, 0, sizeof *w);
  if (hz_oct_init(&w->t, log2n, 0.0)) return 1;
  if (hz_cutmap_init(&w->cmap)) return 1;

  /* РАВНОМЕРНОЕ дробление либо LOD по `L = εR` (этап A, шаг 1). Критерий
   * «есть ли здесь поверхность» спрашивает ТОТ ЖЕ `hz_facets_for_box`, что и
   * боковая таблица (К95) — иначе тела расходятся и условие 1:1 ломается. */
  if (lodeps > 0.0) {
    int z0[3] = {0, 0, 0};
    lod_build(&w->t, z0, nc, nc, ROOM / nc, eye, lodeps, ftab, f0, nfac, nb, cutfine);
  } else {
    for (int x = 0; x < nc; x++)
      for (int y = 0; y < nc; y++)
        for (int z = 0; z < nc; z++) {
          int lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
          hz_oct_set_box(&w->t, lo, hi, 1.0);
        }
  }

  /* боковая таблица: ключи ОБЯЗАНЫ идти по возрастанию (Г45) */
  rec_t *recs = calloc((size_t)nc * (size_t)nc * (size_t)nc, sizeof(rec_t));
  if (recs == NULL) return 1;
  int nrec = 0;
  /* ОБХОД ИДЁТ ПО ЛИСТЬЯМ, А НЕ ПО ЕДИНИЧНЫМ ЯЧЕЙКАМ — К93. Отбор фасетов
   * ведётся по КОРОБКЕ ЛИСТА, а не по единичной ячейке. */
  {
    int32_t stack[64][4]; /* lo[3] и size; глубина дерева заведомо меньше */
    int sp = 0;
    stack[sp][0] = stack[sp][1] = stack[sp][2] = 0;
    stack[sp][3] = nc;
    sp = 1;
    while (sp > 0) {
      sp--;
      int32_t blo[3] = {stack[sp][0], stack[sp][1], stack[sp][2]};
      int32_t bsz = stack[sp][3];
      int32_t rlo[3], rsz = 1;
      int32_t ni = hz_oct_leaf_box(&w->t, blo[0], blo[1], blo[2], rlo, &rsz);
      if (ni < 0) continue;
      if (rsz < bsz) { /* лист мельче — спускаемся */
        int32_t h2 = bsz / 2;
        for (int k = 0; k < 8; k++) {
          stack[sp][0] = blo[0] + ((k & 1) ? h2 : 0);
          stack[sp][1] = blo[1] + ((k & 2) ? h2 : 0);
          stack[sp][2] = blo[2] + ((k & 4) ? h2 : 0);
          stack[sp][3] = h2;
          sp++;
        }
        continue;
      }
      int32_t lo[3] = {rlo[0], rlo[1], rlo[2]};
      int32_t hi[3] = {rlo[0] + rsz, rlo[1] + rsz, rlo[2] + rsz};
      int32_t sel[HZ_P3_MAXH];
      /* ХВОСТ «ДВА ТЕЛА»: в веер ячейки идут фасеты ВСЕХ тел, чьи куски её
       * задевают; не влезшее помечается, а не молчит. */
      int ns = 0, nbody = 0, nover_b = 0, first_b = -1;
      for (int b = 0; b < nb; b++) {
        int32_t s2[HZ_P3_MAXH];
        int k2 = hz_facets_for_box(ftab, f0[b], nfac[b], lo, hi, s2, HZ_P3_MAXH);
        if (k2 <= 0) continue;
        nbody++;
        if (first_b < 0) first_b = b;
        if (ns + k2 > HZ_P3_MAXH) {
          nover_b = 1;
          continue;
        }
        for (int j = 0; j < k2; j++)
          sel[ns++] = s2[j];
      }
      if (nbody > 1) {
        w->nboth++;
        if (nover_b) w->nover++;
      }
      if (ns <= 0) continue;
      recs[nrec].cell = ni;
      recs[nrec].nf = ns;
      for (int j = 0; j < ns; j++)
        recs[nrec].f[j] = sel[j];
      nrec++;
    }
  }
  /* СОРТИРОВКА `qsort`, А НЕ ВСТАВКАМИ (К68-оснастка): ключи различны, порядок
   * совпадает с прежним. */
  qsort(recs, (size_t)nrec, sizeof(rec_t), rec_cmp);
  for (int i = 0; i < nrec; i++)
    if (hz_cutmap_add(&w->cmap, recs[i].cell, recs[i].f, recs[i].nf) != 0) {
      /* К93: МОЛЧА УМИРАТЬ НЕЛЬЗЯ. */
      fprintf(stderr, "боковая таблица отвергла запись %d (ячейка %d): ключи не возрастают\n", i,
              recs[i].cell);
      free(recs);
      return 1;
    }
  free(recs);

  if (tr3_mesh_build(&w->mesh, &w->t, fr)) return 1;
  /* Г38: маска полных ячеек строится ОТДЕЛЬНО — в боковой таблице их нет.
   * Обход тоже по ЛИСТЬЯМ (К93). */
  w->solid = calloc((size_t)w->mesh.ncell, 1);
  if (w->solid == NULL) return 1;
  for (int32_t ci = 0; ci < w->mesh.ncell; ci++) {
    int32_t lo[3] = {w->mesh.clo[ci][0], w->mesh.clo[ci][1], w->mesh.clo[ci][2]};
    int32_t hi[3] = {lo[0] + w->mesh.csize[ci], lo[1] + w->mesh.csize[ci],
                     lo[2] + w->mesh.csize[ci]};
    int32_t sel[HZ_P3_MAXH];
    for (int b = 0; b < nb; b++)
      if (hz_facets_for_box(ftab, f0[b], nfac[b], lo, hi, sel, HZ_P3_MAXH) == 0) w->solid[ci] = 1;
  }
  /* К50: ТЕЛО = ПРИМИТИВ (см. доклад); сбор и развёртка смотрят на ОДНО тело. */
  if (tr3_cut_build3(&w->cut, &w->mesh, ftab, &w->cmap, w->solid, NULL, NULL, stab)) return 1;
  return 0;
}

#define NSTAGE 12
static stage g_st[NSTAGE];
static int g_ns = 0;
static double g_mark = 0.0;

static void stage_mark(void);
static void stage_add(const char *name, int per_frame);

static void stage_mark(void) {
  g_mark = now();
}

static void stage_add(const char *name, int per_frame) {
  double t = now();
  if (g_ns < NSTAGE) {
    g_st[g_ns].name = name;
    g_st[g_ns].t = t - g_mark;
    g_st[g_ns].per_frame = per_frame;
    g_ns++;
  }
  g_mark = t;
}

int main(int argc, char **argv) {
  int W = argc > 1 ? atoi(argv[1]) : 640;
  int H = argc > 2 ? atoi(argv[2]) : 480;
  int nmu = argc > 3 ? atoi(argv[3]) : 2;
  const char *out = argc > 4 ? argv[4] : "img/room_sphere.ppm";
  /* ЗЕРКАЛЬНАЯ ДОЛЯ ШАРОВ. 0 — диффузные, как было. >0 — зеркала, и тогда
   * диффузная доля обнуляется: складывать их без модели Френеля значило бы
   * учесть энергию дважды. */
  double spec = argc > 5 ? atof(argv[5]) : 0.0;
  int maxbounce = argc > 6 ? atoi(argv[6]) : 8;
  /* сколько КАДРОВ отрисовать движущейся камерой при неподвижной сцене */
  int nframe = argc > 7 ? atoi(argv[7]) : 1;
  /* ОГРАНЁННОЕ ТЕЛО: примитив выключен, тело есть многогранник (см. ray3.h) */
  int facet_only = argc > 8 ? atoi(argv[8]) : 0;
  /* К68: РАЗМЕР ЭЛЕМЕНТА — ПАРАМЕТР, А НЕ КОНСТАНТА. Мир при этом не меняется
   * (см. ROOM выше), поэтому свип по `log2n` меряет ровно сетку. */
  int log2n = argc > 9 ? atoi(argv[9]) : 4;
  /* К65: ограничитель положительности параметром — надо измерить, активен ли он
   * на СХОДИМОСТИ рендерной сцены и меняет ли ответ. От этого зависит, нужен ли
   * Крылову внешний нелинейный цикл (этап B плана перехода). */
  int limiter = argc > 10 ? atoi(argv[10]) : 1;
  /* К38: печатать историю невязки каждые `trace` итераций; 0 — молчать */
  int trace = argc > 11 ? atoi(argv[11]) : 0;
  /* предел итераций развёртки: диагностический прогон не обязан ждать сходимости */
  int maxit = argc > 12 ? atoi(argv[12]) : 4000;
  /* К76: возмущение света для опыта «тёплый старт против холодного»; 0 — не ставить */
  double warmdelta = argc > 13 ? atof(argv[13]) : 0.0;
  /* ЭТАП B: сравнить Крылов с рядом Неймана на этой сцене; 0 — не ставить */
  int krylov = argc > 14 ? atoi(argv[14]) : 0;
  /* ЭТАП E: число спектральных каналов, 1 или 3 */
  int nch = argc > 17 ? atoi(argv[17]) : 1;
  /* фальсификатор №1: три ОДИНАКОВЫХ канала */
  int chsame = argc > 18 ? atoi(argv[18]) : 0;
  /* фальсификатор №2: канал без источников; -1 — нет такого */
  int chdark = argc > 19 ? atoi(argv[19]) : -1;
  /* фальсификатор №3: канал, излучение которого УДВОЕНО; -1 — нет такого.
   * Задача линейна по источнику, поэтому ответ обязан удвоиться ТОЧНО. */
  int chdbl = argc > 20 ? atoi(argv[20]) : -1;
  /* ДОПУСК РЕШАТЕЛЯ — ПАРАМЕТР, А НЕ КОНСТАНТА. `1e-9` при `φ ≈ 50` есть `2e-11`
   * относительных, тогда как ошибка дискретизации на `16³` измерена в 3%
   * (К68). То есть система решается на девять порядков точнее, чем имеет смысл,
   * и цена этого — итерации. Число назначено, а не выведено: тот же класс, что
   * `ε` в К68. */
  double stol = argc > 15 ? atof(argv[15]) : 1e-9;
  /* МНОЖИТЕЛЬ АЛЬБЕДО СТЕНОК: настоящий довод этапа B не в множителе на этой
   * сцене, а в НЕЧУВСТВИТЕЛЬНОСТИ к альбедо. Ряд Неймана идёт как `ρⁿ` и при
   * `ρ → 1` встаёт (С1, К37); Крылов на печи держал шесть проходов при любой
   * толщине. Проверяется это только свипом по альбедо. */
  double albs = argc > 16 ? atof(argv[16]) : 1.0;
  /* ЭТАП A, ШАГ 2: сверить ПРОЕКЦИОННЫЙ сбор с маршем НА ТОЙ ЖЕ СЕТКЕ (К53) */
  int raster = argc > 23 ? atoi(argv[23]) : 0;
  /* К59: уровень фасетизации сфер параметром. Ошибка силуэта у проекции обязана
   * ПАДАТЬ с ним — иначе растеризуется не то, что думаем. */
  int fsub = argc > 24 ? atoi(argv[24]) : 2;
  /* ЭТАП A, ШАГ 1: угловой размер пикселя для правила `L = εR`; 0 — равномерная
   * сетка, как прежде. Разделять A надвое стало можно после К52. */
  double lodeps = argc > 21 ? atof(argv[21]) : 0.0;
  /* К94: дробить РАЗРЕЗАННЫЕ ячейки до предела — разделитель причины */
  int cutfine = argc > 22 ? atoi(argv[22]) : 0;
  /* ЭТАП A, ШАГ 3 (П5): тёплый старт через ПЕРЕНОС на сетке, переякоренной к
   * сдвинутой камере. 0 — блок выключен, прежний мир побитово прежний. */
  int xwarm = argc > 25 ? atoi(argv[25]) : 0;
  /* ШАГ 3-е: прибор затухания по отскокам — серия холодных решений с растущим
   * maxit (stol недостижим), разности φ_k − φ_{k−1} = вклад k-го такта. */
  int bprof = argc > 26 ? atoi(argv[26]) : 0;
  /* ШАГ 3-е (2): видимый порог остановки — reltol вместо абсолютного допуска. */
  int visstop = argc > 27 ? atoi(argv[27]) : 0;
  /* ШАГ 3-е (3): двухфазный решатель — такты 1..K* на точной сетке, хвост на
   * грубой (соседний этаж лестницы), поле переносится tr3_xfer. */
  int twophase = argc > 28 ? atoi(argv[28]) : 0;
  if (log2n < 1 || log2n > 8) {
    fprintf(stderr, "log2n вне [1,8]\n");
    return 1;
  }
  const int nc = 1 << log2n;

  stage_mark();
  hz_frame fr = {{0, 0, 0}, {ROOM / nc, ROOM / nc, ROOM / nc}};

  /* ФАСЕТИЗАЦИЯ ИДЁТ ПЕРЕД ДЕРЕВОМ, И ЭТО НЕ ПЕРЕСТАНОВКА РАДИ УДОБСТВА — К95.
   * Критерий дробления обязан спрашивать ТО ЖЕ ТЕЛО, что и боковая таблица.
   * Прежде он проверял аналитическую СФЕРУ, а таблица — ФАСЕТЫ, которые торчат
   * наружу примитива на `dmax`: коробка чуть снаружи сферы не дробилась, но
   * фасеты получала, становилась КРУПНОЙ разрезанной ячейкой рядом с мелкими и
   * ломала условие 1:1. Это дословно К21 («один индекс обслуживает два разных
   * тела»), только на другом конце. */
  hz_surftab stab;
  hz_facettab ftab;
  if (hz_surftab_init(&stab) || hz_facettab_init(&ftab)) return 1;
  /* ДВА ТЕЛА: тени одного на другом и переотражения между ними — то, ради чего
   * развёртка и нужна. Материал ячейки, куда попадают фасеты ОБОИХ тел, есть
   * ОБЪЕДИНЕНИЕ: флюид считает cut3 (cu->nunion), марш видит оба тела (К27).
   * Класть оба тела в веер раньше было нельзя — модель реза умела только
   * пересечение; теперь можно, и тела допустимо даже касатьcя. */
  const int NB = 2;
  double sc[2][3] = {{5.6, 9.2, 4.2}, {10.9, 7.4, 3.1}};
  double sr[2] = {3.0, 1.9};
  int32_t f0[2], nfac[2];
  for (int b = 0; b < NB; b++) {
    hz_surf sp = {HZ_SURF_SPHERE, {sc[b][0], sc[b][1], sc[b][2], sr[b], 0, 0, 0}, 1, 0};
    int32_t si = hz_surftab_add(&stab, &sp);
    f0[b] = 0;
    nfac[b] = hz_surf_facet_sphere(&ftab, &fr, sc[b], sr[b], fsub, HZ_FIT_MEAN_SAGITTA, si, &f0[b]);
    if (nfac[b] <= 0) return 1;
  }

  stage_add("фасетизация примитивов", 0);

  /* МИР №1: сетка, боковая таблица, меш, маска, разрез — одной функцией
   * (этап A, шаг 3: то же самое делается ВТОРОЙ раз под сдвинутую камеру). */
  world w1;
  {
    double leye[3] = {8.0, 0.6, 7.2};
    if (world_build(&w1, log2n, nc, &fr, &stab, &ftab, f0, nfac, NB, lodeps, leye, cutfine))
      return 1;
  }
  tr3_mesh *mesh = &w1.mesh;
  tr3_cut *cut = &w1.cut;
  hz_cutmap *cmap = &w1.cmap;
  int nboth = w1.nboth, nover = (int)w1.nover; /* solid читает только world_build */
  stage_add("дерево, боковая таблица, меш, маска, РАЗРЕЗ", 0);
  tr3_dirs dirs;
  if (tr3_dirs_product(&dirs, nmu, nmu)) return 1;
  stage_add("набор ординат", 0);

  double *sig_t = calloc((size_t)mesh->ncell, sizeof(double));
  double *sig_s = calloc((size_t)mesh->ncell, sizeof(double));
  if (sig_t == NULL || sig_s == NULL) return 1;
  for (int32_t c = 0; c < mesh->ncell; c++) {
    sig_t[c] = 0.015; /* лёгкая дымка: среда участвует, но не глушит */
    sig_s[c] = 0.012; /* альбедо 0.8, а не 1: при 1 ряд по рассеянию сходится как
                       * ρⁿ и упирается ровно в то, ради чего С1 и заводил DSA */
    ;
  }
  /* КОМНАТА ПО КАНАЛАМ (этап E). Прежде считался ОДИН спектральный канал, и
   * «красноватая стена» была красноватой только по яркости — а на картинке
   * стояла ЛОЖНАЯ РАСКРАСКА скалярной величины палитрой.
   *
   * Каналы суть НЕЗАВИСИМЫЕ задачи переноса: развёртка гоняется по разу на
   * канал со своими свойствами. Смешаться они не могут ПО ПОСТРОЕНИЮ — общего
   * состояния у них нет вовсе, — и это довод в пользу такой раскладки, а не
   * следствие лени.
   *
   * ПРИ `nch = 1` ВСЁ ОБЯЗАНО ОСТАТЬСЯ ПРЕЖНИМ ПОБИТОВО: нулевой канал несёт
   * ровно те числа, что стояли здесь до правки. */
  double wrc[TR3_MAXCH][6] = {{0.72, 0.35, 0.72, 0.72, 0.65, 0.05},
                              {0.72, 0.20, 0.55, 0.72, 0.65, 0.05},
                              {0.72, 0.18, 0.42, 0.72, 0.65, 0.05}};
  double wec[TR3_MAXCH][6] = {{0, 0, 0, 0, 0, 6.0}, {0, 0, 0, 0, 0, 5.4}, {0, 0, 0, 0, 0, 4.5}};
  /* ТРИ ОДИНАКОВЫХ КАНАЛА — режим фальсификатора №1: буферы обязаны совпасть с
   * одноканальным прогоном ПОБИТОВО, иначе многоканальность что-то трогает. */
  if (chsame)
    for (int c = 1; c < TR3_MAXCH; c++)
      for (int i = 0; i < 6; i++) {
        wrc[c][i] = wrc[0][i];
        wec[c][i] = wec[0][i];
      }
  /* КАНАЛ БЕЗ ИСТОЧНИКОВ — режим фальсификатора №2: обязан выйти РОВНО нулём. */
  if (chdark >= 0 && chdark < TR3_MAXCH)
    for (int i = 0; i < 6; i++)
      wec[chdark][i] = 0.0;
  if (chdbl >= 0 && chdbl < TR3_MAXCH)
    for (int i = 0; i < 6; i++)
      wec[chdbl][i] *= 2.0;
  for (int c = 0; c < TR3_MAXCH; c++)
    for (int i = 0; i < 6; i++) {
      wrc[c][i] *= albs;
      if (wrc[c][i] > 0.999) wrc[c][i] = 0.999;
    }
  double *wr = wrc[0];
  double *we = wec[0];
  double *frho = calloc((size_t)ftab.n, sizeof(double));
  double *femit = calloc((size_t)ftab.n, sizeof(double));
  if (frho == NULL || femit == NULL) return 1;
  double *fspec = calloc((size_t)ftab.n, sizeof(double));
  if (fspec == NULL) return 1;
  for (int32_t i = 0; i < ftab.n; i++) {
    /* К5: для РАЗВЁРТКИ зеркало есть ЧЁРНОЕ тело — отражённое направление в
     * наборе ординат отсутствует, развёртка на такой границе обрывается, и
     * зеркальная энергия из объёмного решения выпадает. Комната с зеркалами
     * выходит темнее физической ровно на эту долю; это граница метода, а не
     * недосмотр, и она записана в gather3.h. */
    frho[i] = spec > 0.0 ? 0.0 : 0.78;
    fspec[i] = spec;
  }

  tr3_problem prob = {.m = mesh,
                      .d = &dirs,
                      .cut = cut,
                      .facet_rho = frho,
                      .facet_emit = femit,
                      .nfacet = ftab.n,
                      .sig_t = sig_t,
                      .sig_s = sig_s,
                      .wall_rho = wr,
                      .wall_emit = we,
                      .limiter = limiter,
                      .trace = trace};
  double *phi = calloc((size_t)mesh->ncell * 4, sizeof(double));
  if (phi == NULL) return 1;
  tr3_stats st;
  memset(&st, 0, sizeof st);
  printf("ячеек %d, граней %d, поверхностных элементов %d, направлений %d, 1:1 нарушений %d, "
         "ячеек с НЕСКОЛЬКИМИ телами %d (объединение: %d, nbad: %d, веер не влез: %d)\n",
         mesh->ncell, mesh->nf, cut->nse, dirs.n, cut->nbad, nboth, cut->nunion, cut->nbad, nover);
  /* ХРАНИМОЕ ПО КАНАЛАМ, СПЛОШНЫМ МАССИВОМ [nch][4·nf] — та же раскладка, что
   * ждёт сбор (gather3.h). Развёртка о каналах не знает вовсе: она вызывается
   * по разу со своими свойствами, и это ровно то, что делает смешение
   * невозможным по построению. */
  double *boutc = calloc((size_t)nch * (size_t)mesh->nf * 4, sizeof(double));
  double *soutc = calloc((size_t)nch * (size_t)(cut->nse > 0 ? cut->nse : 1) * 4, sizeof(double));
  if (boutc == NULL || soutc == NULL) return 1;
  stage_mark();
  int rc = 0;
  double t_sweep = 0.0;
  memset(&st, 0, sizeof st);
  for (int ch = 0; ch < nch; ch++) {
    tr3_problem pc = prob;
    pc.wall_rho = wrc[ch];
    pc.wall_emit = wec[ch];
    if (visstop) pc.reltol = HZ_SWEEP_VISIBLE_RELTOL; /* ШАГ 3-е (2), см. sweep3.h */
    tr3_stats sc2;
    double tc0 = now();
    int rcc = tr3_sweep_solve(&pc, maxit, stol, phi, &sc2);
    t_sweep += now() - tc0;
    if (rcc != 0) rc = rcc;
    memcpy(boutc + (size_t)ch * (size_t)mesh->nf * 4, sc2.bout,
           (size_t)mesh->nf * 4 * sizeof(double));
    if (cut->nse > 0)
      memcpy(soutc + (size_t)ch * (size_t)cut->nse * 4, sc2.sout,
             (size_t)cut->nse * 4 * sizeof(double));
    if (ch + 1 == nch) {
      st = sc2; /* последний канал: баланс и статистика докладываются по нему */
    } else {
      free(sc2.bout);
      free(sc2.sout);
    }
    if (nch > 1)
      printf("  канал %d: итераций %d, невязка %.2e, баланс на втекшее %.3e\n", ch, sc2.iters,
             sc2.resid, sc2.pin > 0.0 ? fabs(sc2.balance) / sc2.pin : 0.0);
  }
  stage_add("РАЗВЁРТКА (итерация по рассеянию)", 0);
  printf("развёртка (ХОЛОДНЫЙ СТАРТ, каналов %d): код %d, итераций %d, невязка %.2e, срезок %d, "
         "откатов %d, "
         "%.2f с (%.2f мс на итерацию)\n",
         nch, rc, st.iters, st.resid, st.nclip, st.nfallback, t_sweep,
         1e3 * t_sweep / (double)(st.iters + 1) / (double)nch);
  printf("баланс: втекло %.4f, вытекло %.4f, поглощено средой %.4f, ушло в поверхности %.4f, "
         "отдано ими %.4f\n",
         st.pin, st.pout, st.pabs, st.psin, st.psout);
  printf("        невязка %.3e, она же на втекшее %.3e\n", st.balance,
         st.pin > 0.0 ? fabs(st.balance) / st.pin : 0.0);
  /* РАЗДЕЛИТЕЛЬ «ВРЁТ РЕШЕНИЕ ИЛИ ВРЁТ СБОР» (К94). Печатается максимум
   * ХРАНИМОГО радианса в углах граней и на поверхностных элементах. Если он
   * разумен, а картинка показывает больше — виноват СБОР, а не развёртка, и
   * искать надо в чтении DG1, а не в схеме. */
  {
    double bmax = 0.0, smax = 0.0;
    for (int32_t f = 0; f < mesh->nf; f++) {
      if (mesh->f[f].cb >= 0) continue;
      int32_t cc = mesh->f[f].ca;
      double s = (double)mesh->csize[cc];
      (void)s;
      for (int k = 0; k < 8; k++) {
        double v = st.bout[f * 4];
        for (int a = 0; a < 3; a++)
          v += ((k >> a) & 1 ? 0.5 : -0.5) * st.bout[f * 4 + 1 + a];
        if (v > bmax) bmax = v;
      }
    }
    for (int32_t e = 0; e < cut->nse; e++)
      for (int k = 0; k < 8; k++) {
        double v = st.sout[e * 4];
        for (int a = 0; a < 3; a++)
          v += ((k >> a) & 1 ? 0.5 : -0.5) * st.sout[e * 4 + 1 + a];
        if (v > smax) smax = v;
      }
    printf("К94: max ХРАНИМОГО в углах — стенки %.4f, поверхности %.4f "
           "(излучение стенки %.2f)\n",
           bmax, smax, we[5]);
  }
  /* К65: САМА ВЕЛИЧИНА, РАДИ КОТОРОЙ ОГРАНИЧИТЕЛЬ СТОИТ. Картинка её не видит:
   * она читает хранимое на ПОВЕРХНОСТЯХ, а ограничитель следит за полем в
   * ОБЪЁМЕ. Печатается минимум φ по углам ячеек — то самое, что уходит в минус
   * около границы тени, — и минимум хранимого. Без этих двух чисел вопрос
   * «нужен ли Крылову внешний нелинейный цикл» решается на глаз. */
  {
    double phimin = 1e300;
    int nneg = 0;
    for (int32_t c = 0; c < mesh->ncell; c++) {
      for (int k = 0; k < 8; k++) {
        double v = phi[4 * c];
        for (int a = 0; a < 3; a++)
          v += ((k >> a) & 1 ? 0.5 : -0.5) * phi[4 * c + 1 + a];
        if (v < phimin) phimin = v;
        if (v < 0.0) {
          nneg++;
          break;
        }
      }
    }
    printf("К65: min φ по углам ячеек %.6e, ячеек с отрицательным углом %d из %d\n", phimin, nneg,
           mesh->ncell);
  }
  /* --- ШАГ 3-е (bprof): ЗАТУХАНИЕ ВКЛАДОВ ПО ОТСКОКАМ ---
   * d_k = φ_k − φ_{k−1} есть вклад k-го такта (первый такт смешивает прямой
   * свет и первое рассеяние). Предсказания до прогона: П1 d_k убывает с
   * отношением ≈ ρ=0.63; П2 отношение (наклоны/константа) внутри такта убывает
   * с k — такт глаже предыдущего; П3 (НК) первый такт наибольший. */
  if (bprof) {
    enum { NB2 = 12 };
    double *sn = calloc((size_t)(NB2 + 1) * (size_t)mesh->ncell * 4, sizeof(double));
    if (!sn) return 1;
    g_bp_sn = sn;
    g_bp_cap = mesh->ncell * 4;
    g_bp_n = 0;
    tr3_problem pbp = prob;
    pbp.limiter = 0;  /* чистый ряд, без нелинейных веток */
    pbp.maxp_off = 1; /* иначе такты — не элементы ряда */
    pbp.pass_fn = bp_capture;
    tr3_stats sbp;
    double *phib = calloc((size_t)mesh->ncell * 4, sizeof(double));
    if (!phib) return 1;
    if (tr3_sweep_solve(&pbp, NB2 + 1, 1e-13, phib, &sbp) != 0) return 1;
    int nsnaps = g_bp_n;
    free(sbp.bout);
    free(sbp.sout);
    free(sbp.eirr);
    free(phib);
    double *zeros = calloc((size_t)mesh->ncell * 4, sizeof(double));
    if (!zeros) return 1;
    double prev_amp = 0.0;
    for (int n = 1; n <= nsnaps; n++) {
      const double *cur = sn + (size_t)(n - 1) * (size_t)mesh->ncell * 4;
      const double *prev = n == 1 ? zeros : sn + (size_t)(n - 2) * (size_t)mesh->ncell * 4;
      double amp = 0.0, sl = 0.0, cn = 0.0;
      for (int32_t c = 0; c < mesh->ncell; c++) {
        double dv[4];
        for (int i = 0; i < 4; i++)
          dv[i] = cur[4 * c + i] - prev[4 * c + i];
        if (fabs(dv[0]) > amp) amp = fabs(dv[0]);
        cn += fabs(dv[0]);
        for (int j = 1; j < 4; j++)
          sl += fabs(dv[j]);
      }
      double slope_mean = sl / (double)mesh->ncell / 3.0;
      double const_mean = cn / (double)mesh->ncell;
      double ratio = const_mean > 0.0 ? slope_mean / const_mean : 0.0;
      printf("  [bprof %2d] d_k max %.4e, сред.конст %.4e, накл/конст %.3f", n, amp, const_mean,
             ratio);
      if (n > 1 && prev_amp > 0.0) printf(", к предыдущему: %.3f", amp / prev_amp);
      printf("\n");
      prev_amp = amp;
    }
    free(zeros);
    free(sn);
  }
  if (rc != 0) return 1;

  /* --- ЭТАП B: КРЫЛОВ ПРОТИВ РЯДА НЕЙМАНА, НА ТОЙ ЖЕ СЦЕНЕ ---
   *
   * Единица сравнения — ПРОХОДЫ РАЗВЁРТКИ, а не итерации метода: BiCGStab
   * делает два прохода на итерацию, и считать итерации значило бы сравнивать
   * разное.
   *
   * СРАВНЕНИЕ ИДЁТ С РЯДОМ НЕЙМАНА БЕЗ ОГРАНИЧИТЕЛЯ, и это не поблажка, а
   * требование: Крылову нужен ЛИНЕЙНЫЙ оператор, поэтому ограничитель у него
   * выключен, а К65 измерила, что с ним и без него неподвижные точки РАЗНЫЕ
   * (0.2% в среднем, 6.8% в худшем пикселе). Сравнивать с ограниченным рядом
   * значило бы приписать методу чужую разницу. */
  if (krylov) {
    double *phi_n = calloc((size_t)mesh->ncell * 4, sizeof(double));
    double *phi_k = calloc((size_t)mesh->ncell * 4, sizeof(double));
    double *bk = calloc((size_t)mesh->nf * 4, sizeof(double));
    double *sk = calloc((size_t)(cut->nse > 0 ? cut->nse : 1) * 4, sizeof(double));
    if (phi_n == NULL || phi_k == NULL || bk == NULL || sk == NULL) return 1;

    tr3_problem pn = prob;
    pn.limiter = 0;
    /* К102: сравнение обязано быть ЛИНЕЙНЫМ С ОБОИХ СТОРОН. maxp_off не был
     * выставлен — развёртка применяла §735-проекции принципа максимума
     * (нелинейные клампы), её ответ не удовлетворял линейному оператору
     * стенда (невязка 14–114), и объёмные φ расходились с Крыловым на сотни.
     * Крылов внутри всегда линеен (krylov3.c: limiter=0, maxp_off=1). */
    pn.maxp_off = 1;
    tr3_stats stn;
    double t0n = now();
    int rcn = tr3_sweep_solve(&pn, maxit, stol, phi_n, &stn);
    double tn = now() - t0n;

    tr3_kstats stk;
    double t0k = now();
    int rck = tr3_krylov_solve(&prob, 200, stol, phi_k, bk, sk, &stk);
    double tk = now() - t0k;

    double dmax = 0.0, pmax = 0.0, dbm = 0.0, bmax = 0.0;
    for (int32_t i = 0; i < mesh->ncell * 4; i++) {
      double d = fabs(phi_n[i] - phi_k[i]);
      if (d > dmax) dmax = d;
      if (fabs(phi_n[i]) > pmax) pmax = fabs(phi_n[i]);
    }
    for (int32_t i = 0; i < mesh->nf * 4; i++) {
      double d = fabs(stn.bout[i] - bk[i]);
      if (d > dbm) dbm = d;
      if (fabs(stn.bout[i]) > bmax) bmax = fabs(stn.bout[i]);
    }
    printf("ЭТАП B, КРЫЛОВ ПРОТИВ РЯДА НЕЙМАНА (оба БЕЗ ограничителя):\n");
    printf("     НЕЙМАН:  код %d, ПРОХОДОВ %d, невязка %.2e, застой %d, %.2f с\n", rcn, stn.iters,
           stn.resid, stn.stalled, tn);
    printf("     КРЫЛОВ:  код %d, ПРОХОДОВ %ld (итераций %d), невязка %.2e, %s, %.2f с  (×%.2f)\n",
           rck, stk.npass, stk.iters, stk.resid, stk.why, tk, tk > 0.0 ? tn / tk : 0.0);
    /* ПРОВЕРКА ПРИБОРА ПЕРЕД СПОРОМ О ЧИСЛАХ: сошедшийся ответ ряда Неймана
     * обязан давать под тем же оператором НУЛЕВУЮ невязку. Если не даёт —
     * неверен оператор или правая часть, и сравнивать проходы бессмысленно. */
    double relN = -1.0, relK = -1.0;
    tr3_krylov_residual(&prob, phi_n, stn.bout, stn.sout, &relN);
    tr3_krylov_residual(&prob, phi_k, bk, sk, &relK);
    printf("     ПРИБОР: невязка ответа НЕЙМАНА под оператором Крылова %.3e; ответа Крылова %.3e\n",
           relN, relK);
    printf("     СОВПАДЕНИЕ ОТВЕТОВ: max|Δφ| %.3e при |φ|max %.3e; "
           "max|Δbout| %.3e при |bout|max %.3e\n",
           dmax, pmax, dbm, bmax);
    free(stn.bout);
    free(stn.sout);
    free(stn.eirr);
    /* stk — tr3_kstats: bout/sout живут в bk/sk вызывающего и freed выше */
    free(phi_n);
    free(phi_k);
    free(bk);
    free(sk);
  }

  /* --- К76: ТЁПЛЫЙ СТАРТ ПРОТИВ ХОЛОДНОГО ПРИ ДВИЖУЩЕМСЯ СВЕТЕ ---
   *
   * Свет меняется на `warmdelta`, и та же задача решается ДВАЖДЫ: от нуля и от
   * поля предыдущего кадра (включая `bout`/`sout` — без них тёплый старт на
   * сцене с отражением тёплым не является).
   *
   * ПРЕДСКАЗАНИЕ ЗАПИСАНО ДО ПРОГОНА и выведено, а не угадано: число итераций
   * есть `ln(tol/e₀)/ln(ρ)`, тёплый старт уменьшает ТОЛЬКО `e₀`. При `ρ = 0.63`,
   * `tol = 1e-9` и `φ ≈ 16` возмущение `δ` даёт `e₀ ≈ δ·φ`, то есть экономию
   * `ln(37/(δ·16))/ln(1/0.63)` итераций: около 7 при `δ = 0.1`, 12 при 0.01 и
   * 17 при 0.001. Ждать `×20` неоткуда.
   *
   * НЕГАТИВНЫЙ КОНТРОЛЬ ВСТРОЕН: оба решения обязаны СОВПАСТЬ. Если тёплое
   * отличается от холодного больше допуска, тёплый старт меняет ОТВЕТ, а не
   * только путь к нему, и тогда он незаконен. */
  if (warmdelta > 0.0) {
    double we2[6];
    for (int i = 0; i < 6; i++)
      we2[i] = we[i] * (1.0 + warmdelta);
    double *phi_c = calloc((size_t)mesh->ncell * 4, sizeof(double));
    double *phi_w = calloc((size_t)mesh->ncell * 4, sizeof(double));
    if (phi_c == NULL || phi_w == NULL) return 1;
    memcpy(phi_w, phi, (size_t)mesh->ncell * 4 * sizeof(double));

    tr3_problem p2 = prob;
    p2.wall_emit = we2;
    tr3_stats sc_st, sw_st;
    double t1 = now();
    int rc1 = tr3_sweep_solve(&p2, maxit, 1e-9, phi_c, &sc_st);
    double tc = now() - t1;

    p2.warm_start = 1;
    p2.bout_in = st.bout;
    p2.sout_in = st.sout;
    t1 = now();
    int rc2 = tr3_sweep_solve(&p2, maxit, 1e-9, phi_w, &sw_st);
    double tw = now() - t1;

    double dmax = 0.0, pmax = 0.0;
    for (int32_t i = 0; i < mesh->ncell * 4; i++) {
      double d = fabs(phi_c[i] - phi_w[i]);
      if (d > dmax) dmax = d;
      if (fabs(phi_c[i]) > pmax) pmax = fabs(phi_c[i]);
    }
    printf("К76 ТЁПЛЫЙ СТАРТ, возмущение света %.4g:\n", warmdelta);
    printf("     ХОЛОДНЫЙ: код %d, итераций %d, невязка %.2e, застой %d, %.2f с\n", rc1,
           sc_st.iters, sc_st.resid, sc_st.stalled, tc);
    printf("     ТЁПЛЫЙ:   код %d, итераций %d, невязка %.2e, застой %d, %.2f с  (×%.2f)\n", rc2,
           sw_st.iters, sw_st.resid, sw_st.stalled, tw, tc > 0.0 ? tc / tw : 0.0);
    printf("     СОВПАДЕНИЕ ОТВЕТОВ (негативный контроль): max|Δφ| %.3e при |φ|max %.3e\n", dmax,
           pmax);
    /* ШАГ 3/5: ЦЕНА ОТСУТСТВИЯ ПЕРЕНОСА bout/sout. Тот же тёплый старт, но
     * БЕЗ поверхностного хранимого (bout_in = NULL): разница итераций с
     * строкой «ТЁПЛЫЙ» выше есть ВЕРХНЯЯ граница того, что мог бы дать
     * перенос bout/sout между сетками — на той же сетке поверхностное поле
     * доступно ТОЧНО, между сетками оно было бы лишь приближённым. */
    {
      double *phi_x = calloc((size_t)mesh->ncell * 4, sizeof(double));
      if (phi_x == NULL) return 1;
      memcpy(phi_x, phi, (size_t)mesh->ncell * 4 * sizeof(double));
      tr3_problem p3 = p2;
      p3.warm_start = 1;
      p3.bout_in = NULL;
      p3.sout_in = NULL;
      tr3_stats sx;
      double t3 = now();
      int rc3 = tr3_sweep_solve(&p3, maxit, 1e-9, phi_x, &sx);
      printf("     ТЁПЛЫЙ БЕЗ bout/sout: код %d, итераций %d, %.2f с — разница с ТЁПЛЫМ "
             "есть цена поверхностного хранимого\n",
             rc3, sx.iters, now() - t3);
      double dx = 0.0;
      for (int32_t i = 0; i < mesh->ncell * 4; i++) {
        double d = fabs(phi_x[i] - phi_w[i]);
        if (d > dx) dx = d;
      }
      printf("     и против ТЁПЛОГО с bout/sout: max|Δφ| %.3e\n", dx);
      free(sx.bout);
      free(sx.sout);
      free(sx.eirr);
      free(phi_x);
    }
    free(sc_st.bout);
    free(sc_st.sout);
    free(sw_st.bout);
    free(sw_st.sout);
    free(phi_c);
    free(phi_w);
  }

  /* --- ШАГ 3-е (3): ДВУХФАЗНЫЙ РЕШАТЕЛЬ ---
   * Фаза 1: такты 1..K* на ТОЧНОМ мире (прямой свет и первый отскок несут
   * угловую структуру — лепестки, тени; bprof: отношения тактов садятся к
   * третьему). Фаза 2: хвост ряда на ГРУБОМ мире (соседний этаж лестницы,
   * CMULT — минимальный шаг, не подбор), поле переносится tr3_xfer
   * (ограничение точно в моменты). Сравнение — с односеточным видимым
   * решением на точной сетке (обратная пролонгация Крылову не нужна:
   * собираем поле на точной сетке). Бюджет ошибки: хвост фазы 2 амплитудой
   * <= ρ^{K*}·|φ| считается на грубой сетке. */
  if (twophase) {
    const int KSTAR = 2;
    const double CMULT = 2.0;
    /* база: односеточное ВИДИМОЕ решение на точной сетке */
    double *phi_vis = calloc((size_t)mesh->ncell * 4, sizeof(double));
    double *phi_p1 = calloc((size_t)mesh->ncell * 4, sizeof(double));
    if (!phi_vis || !phi_p1) return 1;
    /* Все фазы — PRODUCTION-настройки (limiter + maxp): это Нейман, линейность
     * ему не нужна (К102 требовала её только для Крылова), а без maxp грубая
     * сетка даёт DG1-выбросы (замерено: |φ|=713 при норме 85 — §734-класс). */
    tr3_problem pv = prob;
    pv.reltol = HZ_SWEEP_VISIBLE_RELTOL;
    tr3_stats stv;
    double tv0 = now();
    if (tr3_sweep_solve(&pv, maxit, stol, phi_vis, &stv) != 0) return 1;
    double tvis = now() - tv0;
    /* фаза 1 */
    tr3_problem p1 = pv;
    p1.reltol = 0.0;
    p1.warm_start = 0;
    tr3_stats st1;
    double tf0 = now();
    if (tr3_sweep_solve(&p1, KSTAR, 0.0, phi_p1, &st1) != 0) return 1;
    double tfine = now() - tf0;
    /* грубый мир и перенос */
    world wc;
    {
      double leye[3] = {8.0, 0.6, 7.2};
      if (world_build(&wc, log2n, nc, &fr, &stab, &ftab, f0, nfac, NB,
                      lodeps > 0.0 ? lodeps * CMULT : 0.0, leye, cutfine))
        return orbit_fail(__LINE__);
    }
    int32_t ncs = wc.mesh.ncell;
    double *phi_c = calloc((size_t)ncs * 4, sizeof(double));
    double *phi_back = calloc((size_t)mesh->ncell * 4, sizeof(double));
    double *sigct = calloc((size_t)ncs, sizeof(double));
    double *sigcs = calloc((size_t)ncs, sizeof(double));
    if (!phi_c || !phi_back || !sigct || !sigcs) return orbit_fail(__LINE__);
    for (int32_t c = 0; c < ncs; c++) {
      sigct[c] = 0.015;
      sigcs[c] = 0.012;
    }
    double tx0 = now();
    if (tr3_xfer(&wc.mesh, phi_c, mesh, phi_p1)) return orbit_fail(__LINE__);
    double txfer = now() - tx0;
    /* фаза 2: хвост на грубой сетке до видимого порога */
    tr3_problem pc2 = {.m = &wc.mesh,
                       .d = &dirs,
                       .cut = &wc.cut,
                       .facet_rho = frho,
                       .facet_emit = femit,
                       .nfacet = ftab.n,
                       .sig_t = sigct,
                       .sig_s = sigcs,
                       .wall_rho = wr,
                       .wall_emit = we,
                       .limiter = limiter,
                       .reltol = HZ_SWEEP_VISIBLE_RELTOL,
                       .warm_start = 1,
                       .trace = 0};
    tr3_stats st2;
    double tc0 = now();
    if (tr3_sweep_solve(&pc2, maxit, stol, phi_c, &st2) != 0) return orbit_fail(__LINE__);
    double tcoarse = now() - tc0;
    /* обратная пролонгация на точную сетку для сверки */
    if (tr3_xfer(mesh, phi_back, &wc.mesh, phi_c)) return orbit_fail(__LINE__);
    double dvis = 0.0, dvis_sum = 0.0, dfull = 0.0, pmax = 0.0, srel = 0.0;
    for (int32_t i = 0; i < mesh->ncell * 4; i++) {
      double a = fabs(phi_back[i] - phi_vis[i]);
      if (a > dvis) dvis = a;
      dvis_sum += a;
      double b = fabs(phi_back[i] - phi[i]);
      if (b > dfull) dfull = b;
      if (fabs(phi_vis[i]) > pmax) pmax = fabs(phi_vis[i]);
      srel += fabs(phi_vis[i]);
    }
    double dvis_mean = dvis_sum / (double)(mesh->ncell * 4);
    printf("ШАГ 3-е (3) ДВУХФАЗНЫЙ РЕШАТЕЛЬ (ε точной сетки ×%.1f для хвоста):\n", CMULT);
    printf("     одно-сеточно видимый: проходов %d, %.2f с\n", stv.iters, tvis);
    printf("     двухфазно: фаза1 %d тактов %.2f с + перенос %.4f с + фаза2 %d тактов %.2f с; "
           "итого решатель %.2f с (×%.2f)\n",
           st1.iters, tfine, txfer, st2.iters, tcoarse, tfine + txfer + tcoarse,
           tvis > 0.0 ? tvis / (tfine + txfer + tcoarse) : 0.0);
    printf("     φ: двухфазно-пролонгировано против одно-сеточного: max|Δ| %.3e, сред |Δ| %.3e "
           "при |φ|max %.3e (средн. |φ| %.3e)\n",
           dvis, dvis_mean, pmax, srel / (double)(mesh->ncell * 4));
    printf("     против ПОЛНОГО решения: max|Δ| %.3e\n", dfull);
    free(sigcs);
    free(sigct);
    free(phi_back);
    free(phi_c);
    world_free(&wc);
    free(phi_p1);
    free(phi_vis);
  }

  /* --- сбор по пикселю --- */ /* --- сбор по пикселю --- */
  tr3_scene scn = {.tree = &w1.t,
                   .fr = fr,
                   .sigma = NULL,
                   .st = &stab,
                   .ft = &ftab,
                   .cm = cmap,
                   .facet_only = facet_only};
  tr3_camera cam;
  double eye[3] = {8.0, 0.6, 7.2}, at[3] = {8.2, 9.5, 3.6}, up[3] = {0, 0, 1};
  if (tr3_camera_look(&cam, eye, at, up, 1.3, W, H)) return 1;
  /* ТРИ БУФЕРА, А НЕ ОДИН: `buf` есть нулевой канал, и одноканальный путь
   * пользуется им как прежде. */
  double *buf = calloc((size_t)W * (size_t)H * (size_t)TR3_MAXCH, sizeof(double));
  if (buf == NULL) return 1;

  /* ИНДЕКС ГРАНИЧНЫХ ГРАНЕЙ ПО СЕТКЕ. Первая редакция искала грань ЛИНЕЙНЫМ
   * перебором всех 13056 граней на КАЖДЫЙ пиксель — это 2.7e10 сравнений на
   * кадр, и сбор стоил 3930 нс на луч при разумных 100-200. Грани лежат на
   * ЦЕЛОЧИСЛЕННОЙ сетке, поэтому индекс прямой: (стенка, u, v) -> грань. */
  int32_t *wallidx = calloc((size_t)6 * (size_t)nc * (size_t)nc, sizeof(int32_t));
  if (wallidx == NULL) return 1;
  for (int32_t i = 0; i < 6 * (int32_t)nc * (int32_t)nc; i++)
    wallidx[i] = -1;
  for (int32_t f = 0; f < mesh->nf; f++) {
    const tr3_face *ff = &mesh->f[f];
    if (ff->cb >= 0) continue;
    int wl = (int)(~ff->cb);
    for (int32_t a = ff->lo[0]; a < ff->hi[0]; a++)
      for (int32_t b = ff->lo[1]; b < ff->hi[1]; b++)
        wallidx[((int32_t)wl * nc + a) * nc + b] = f;
  }

  stage_add("индекс граничных граней", 0);
  double t1 = now();
  /* СБОР ВЫНЕСЕН В МОДУЛЬ (gather3.c): зеркало добавляет в него ЦИКЛ, а цикл
   * надо фальсифицировать, чего внутри main было негде делать. При spec = 0
   * результат обязан совпасть с прежним ПОБИТОВО — это перестановка кода, а не
   * изменение расчёта, и это проверено. */
  tr3_gather gg = {.sc = &scn,
                   .m = mesh,
                   .cut = cut,
                   .bout = boutc,
                   .sout = soutc,
                   .nch = nch,
                   .facet_spec = spec > 0.0 ? fspec : NULL,
                   .nfacet = ftab.n,
                   .wallidx = wallidx,
                   .nwall = nc,
                   .maxbounce = maxbounce,
                   /* К45: среда вдоль луча — ослабление и вклад рассеяния */
                   .sig_t = sig_t,
                   .sig_s = sig_s,
                   .phi = phi};
  long nbtot = 0;
  int nbmax = 0;
  /* ДВИЖУЩАЯСЯ КАМЕРА ПРИ НЕПОДВИЖНОЙ СЦЕНЕ — ЗАМЕР, А НЕ ВЫКЛАДКА.
   *
   * Вопрос «сколько стоит КАДР» нельзя честно ответить одним прогоном: в нём
   * цена кадра неотличима от разовой расстановки. Поэтому камера двигается по
   * дуге, а сцена остаётся на месте, и печатается время КАЖДОГО кадра. Что при
   * этом НЕ пересчитывается, видно прямо по коду: между кадрами меняется только
   * `cam`, а `gg` (сетка, разрез, хранимое поле, индекс граней) не трогается
   * ВООБЩЕ. Развёртка позади цикла и в него не входит.
   *
   * Первый кадр отделён от остальных: он греет кэш, и мешать его с
   * установившимися значило бы завысить цену движения. */
  double t_first = 0.0, t_rest = 0.0;
  for (int fri = 0; fri < nframe; fri++) {
    if (fri > 0) {
      /* дуга вокруг центра комнаты: меняется ТОЛЬКО камера */
      double a = 0.35 * (double)fri / (double)(nframe > 1 ? nframe - 1 : 1);
      double ex = 8.0 + 6.0 * sin(a), ey = 0.6 - 6.0 * (1.0 - cos(a));
      double e2[3] = {ex, ey, 7.2};
      if (tr3_camera_look(&cam, e2, at, up, 1.3, W, H)) return 1;
    }
    double tf = now();
    nbtot = 0;
    nbmax = 0;
    for (int py = 0; py < H; py++)
      for (int px = 0; px < W; px++) {
        double o[3], d[3];
        tr3_camera_ray(&cam, px, py, o, d);
        int nb = 0;
        double px3[TR3_MAXCH] = {0, 0, 0};
        tr3_gather_ray(&gg, o, d, &nb, px3);
        for (int ch = 0; ch < nch; ch++)
          buf[(size_t)ch * (size_t)W * (size_t)H + (size_t)py * (size_t)W + (size_t)px] = px3[ch];
        nbtot += nb;
        if (nb > nbmax) nbmax = nb;
      }
    double dt = now() - tf;
    if (nframe > 1)
      printf("  кадр %2d: %.3f с (%.0f нс на луч, %.2f кадра в секунду)\n", fri, dt,
             1e9 * dt / ((double)W * (double)H), 1.0 / dt);
    if (fri == 0)
      t_first = dt;
    else
      t_rest += dt;
  }
  if (nframe > 1)
    printf("КАДР ПРИ НЕПОДВИЖНОЙ СЦЕНЕ: первый %.3f с, установившийся %.3f с в среднем "
           "по %d кадрам (%.2f кадра в секунду). Между кадрами меняется ТОЛЬКО камера.\n",
           t_first, t_rest / (double)(nframe - 1), nframe - 1,
           (double)(nframe - 1) / (t_rest > 0.0 ? t_rest : 1.0));

  /* В УЧЁТ ИДЁТ УСТАНОВИВШИЙСЯ КАДР, А НЕ СУММА ПО ВСЕМ. Иначе строка «сбор» и
   * доля «на каждый кадр» растут вместе с числом кадров, а это бессмыслица:
   * кадр стоит столько, сколько стоит ОДИН кадр. */
  double t_gather = nframe > 1 ? t_rest / (double)(nframe - 1) : now() - t1;
  (void)t1;
  g_mark = now() - t_gather; /* в этап идёт УСТАНОВИВШИЙСЯ кадр, см. выше */
  stage_add("СБОР ПО ПИКСЕЛЮ (установившийся кадр)", 1);
  /* РАЗДЕЛЬНЫЙ ДОКЛАД ХОЛОДНОГО СТАРТА И КАДРА — требование CLAUDE.md. Поле от
   * камеры не зависит, поэтому поворот камеры стоит ТОЛЬКО сбора. */
  printf("СБОР ПО ПИКСЕЛЮ (кадр): %.3f с на %dx%d = %.2f млн лучей, %.0f нс на луч\n", t_gather, W,
         H, 1e-6 * (double)W * (double)H, 1e9 * t_gather / ((double)W * (double)H));
  {
    double tot = 0.0, per = 0.0;
    for (int i = 0; i < g_ns; i++) {
      tot += g_st[i].t;
      if (g_st[i].per_frame) per += g_st[i].t;
    }
    printf("\n--- ХОЛОДНЫЙ СТАРТ ПО ЭТАПАМ (всё, от расстановки сцены) ---\n");
    for (int i = 0; i < g_ns; i++)
      printf("  %-52s %7.3f с  %5.1f%%  %s\n", g_st[i].name, g_st[i].t, 100.0 * g_st[i].t / tot,
             g_st[i].per_frame ? "НА КАЖДЫЙ КАДР" : "один раз на сцену");
    printf("  %-52s %7.3f с\n", "ИТОГО холодный старт", tot);
    printf("  %-52s %7.3f с   (отношение %.1f)\n", "из них ПОВТОРИТСЯ при движении камеры", per,
           per > 0.0 ? tot / per : 0.0);
  }
  /* ДВА ФАЙЛА, И ЭТО НЕ ИЗБЫТОЧНОСТЬ. Канонический — PFM: РАДИАНС КАК ЕСТЬ,
   * 32-битный float, без нормировки, гаммы и раскраски. По нему можно сверять
   * числа; по PPM нельзя ничего (К19: та же ошибка после тон-маппинга доходит
   * до 255 уровней из 255). PPM пишется рядом и только для глаз.
   * Имя PFM получается заменой расширения у заданного пути. */
  /* --- ЭТАП A, ШАГ 2: ПРОЕКЦИЯ ПРОТИВ МАРША НА ОДНОЙ И ТОЙ ЖЕ СЕТКЕ (К53) ---
   *
   * Сетка, хранимое поле и камера — одни и те же; меняется РОВНО СБОР. Это и
   * есть изолирующая проверка, которую формула «A только вместе» запрещала, а
   * К52 разрешила, показав, что марш под LOD не дорожает. */
  if (raster) {
    double *rb = calloc((size_t)W * (size_t)H * (size_t)TR3_MAXCH, sizeof(double));
    if (rb == NULL) return 1;
    tr3_raster rr = {.m = mesh, .cut = cut, .bout = boutc, .sout = soutc, .nch = nch};
    tr3_rstats rs;
    double tr0 = now();
    if (tr3_raster_render(&rr, &cam, rb, &rs)) return 1;
    double t_ras = now() - tr0;
    printf("ЭТАП A ШАГ 2, ПРОЕКЦИЯ: элементов %ld, пропущено %ld, пикселей %ld, ПУСТЫХ %ld, "
           "%.3f с (%.0f нс на пиксель)\n",
           rs.nelem, rs.nskip, rs.npix, rs.nempty, t_ras,
           1e9 * t_ras / (double)((size_t)W * (size_t)H));
    double dmax2 = 0.0, dsum = 0.0, vmax = 0.0;
    long nin = 0;
    for (size_t i = 0; i < (size_t)W * (size_t)H; i++) {
      double a = buf[i], b2 = rb[i];
      if (fabs(a) > vmax) vmax = fabs(a);
      double e = fabs(a - b2);
      dsum += e;
      nin++;
      if (e > dmax2) dmax2 = e;
    }
    printf("     ПРОТИВ МАРША: |Δ| сред %.3e, макс %.3e при |L|max %.3e\n",
           nin > 0 ? dsum / (double)nin : 0.0, dmax2, vmax);
    {
      char rp[512];
      size_t ln2 = strlen(out);
      if (ln2 + 6 < sizeof rp) {
        memcpy(rp, out, ln2 + 1);
        char *dt = strrchr(rp, (int)0x2E);
        size_t ex = dt != NULL ? (size_t)(dt - rp) : ln2;
        memcpy(rp + ex, "_r.pfm", 7);
        const double *pg2 = nch > 1 ? rb + (size_t)W * (size_t)H : rb;
        const double *pb2 = nch > 2 ? rb + (size_t)2 * (size_t)W * (size_t)H : rb;
        hz_pfm_write(rp, rb, pg2, pb2, W, H);
      }
    }
    free(rb);
  }

  /* --- ЭТАП A, ШАГ 3 (П5): ТЁПЛЫЙ СТАРТ ЧЕРЕЗ ПЕРЕНОС МЕЖДУ LOD-СЕТКАМИ ---
   *
   * Камера сдвигается по той же дуге, что и в кадре, сетка ПЕРЕСТРАИВАЕТСЯ под
   * её новую позицию (`L = εR` якорится в eye), и задача решается ТРИЖДЫ:
   * холодным стартом на старой сетке (уже решено выше), холодным и ТЁПЛЫМ на
   * новой. Тёплый начальное поле получает ПЕРЕНОСОМ `tr3_xfer` со старой
   * сетки, а не нулём. Вопросы, на которые блок отвечает числами:
   *   (а) сколько итераций экономит перенос против холодного старта;
   *   (б) не меняет ли тёплый старт ОТВЕТ (негативный контроль, как в К76);
   *   (в) как далек перенос от ответа, то есть размер начальной ошибки e0.
   * `bout_in` не передаётся сознательно: грани двух сеток не совпадают, их
   * перенос — отдельная работа, и честный замер ценности φ-переноса идёт без него. */
  if (xwarm == 1) {
    world w2;
    double a2 = 0.35, eye2[3] = {8.0 + 6.0 * sin(a2), 0.6 - 6.0 * (1.0 - cos(a2)), 7.2};
    double tb0 = now();
    if (world_build(&w2, log2n, nc, &fr, &stab, &ftab, f0, nfac, NB, lodeps, eye2, cutfine))
      return 1;
    double tbuild = now() - tb0;
    int32_t n2 = w2.mesh.ncell;
    double *phi2c = calloc((size_t)n2 * 4, sizeof(double));
    double *phi2w = calloc((size_t)n2 * 4, sizeof(double));
    double *sig2t = calloc((size_t)n2, sizeof(double));
    double *sig2s = calloc((size_t)n2, sizeof(double));
    if (!phi2c || !phi2w || !sig2t || !sig2s) return 1;
    for (int32_t c = 0; c < n2; c++) {
      sig2t[c] = 0.015;
      sig2s[c] = 0.012;
    }
    double tx0 = now();
    if (tr3_xfer(&w2.mesh, phi2w, mesh, phi)) return 1;
    double txfer = now() - tx0; /* А3: перенос измеряется, а не печатается нулём */
    /* К21-ПРИБОР: доля ячеек новой сетки, у которых в старом дереве есть лист
     * С ТОЙ ЖЕ коробкой. Ровно они получают поле ПОБИТОВОЙ копией (П1 переноса),
     * по ним же считается потолок переиспользования узлов при инкрементальности:
     * пережившая ячейка — кандидат на переживший узел дерева, разрез и элемент. */
    long nsurv = 0;
    for (int32_t c = 0; c < n2; c++) {
      int32_t rlo[3], rsz = 0;
      int32_t ni = hz_oct_leaf_box(&w1.t, w2.mesh.clo[c][0], w2.mesh.clo[c][1], w2.mesh.clo[c][2],
                                   rlo, &rsz);
      if (ni >= 0 && rsz == w2.mesh.csize[c] && rlo[0] == w2.mesh.clo[c][0] &&
          rlo[1] == w2.mesh.clo[c][1] && rlo[2] == w2.mesh.clo[c][2])
        nsurv++;
    }
    printf("     К21: пережило ячеек %ld из %d (%.1f%%) — потолок переиспользования узлов\n", nsurv,
           n2, 100.0 * (double)nsurv / (double)n2);
    /* (в): размер начальной ошибки против сошедшегося ответа, ПО РАЗРЯДАМ:
     * ячейки У ПОВЕРХНОСТИ (есть элементы разреза) против ОБЪЁМНЫХ. Если e0
     * держится поверхностью — виновата коробка против флюида в моментной
     * модели переноса; если и объём — границы теней (граница метода). */
    double e0 = 0.0, pmax2 = 0.0;
    double e0_surf = 0.0, e0_vol = 0.0, sum_surf = 0.0, sum_vol = 0.0;
    long n_surf = 0, n_vol = 0;
    tr3_problem pc2 = {.m = &w2.mesh,
                       .d = &dirs,
                       .cut = &w2.cut,
                       .facet_rho = frho,
                       .facet_emit = femit,
                       .nfacet = ftab.n,
                       .sig_t = sig2t,
                       .sig_s = sig2s,
                       .wall_rho = wr,
                       .wall_emit = we,
                       .limiter = limiter,
                       .trace = 0};
    tr3_stats stc2;
    double tk0 = now();
    int rcc2 = tr3_sweep_solve(&pc2, maxit, stol, phi2c, &stc2);
    double tcool = now() - tk0;
    for (int32_t c = 0; c < n2; c++) {
      double dmaxc = 0.0;
      for (int i = 0; i < 4; i++) {
        double d = fabs(phi2w[4 * c + i] - phi2c[4 * c + i]);
        if (d > dmaxc) dmaxc = d;
      }
      if (w2.cut.sestart[c + 1] > w2.cut.sestart[c]) {
        n_surf++;
        sum_surf += dmaxc;
        if (dmaxc > e0_surf) e0_surf = dmaxc;
      } else {
        n_vol++;
        sum_vol += dmaxc;
        if (dmaxc > e0_vol) e0_vol = dmaxc;
      }
    }
    tr3_stats stw2;
    pc2.warm_start = 1; /* φ уже в phi2w — перенос и есть начальное поле */
    tk0 = now();
    int rcw2 = tr3_sweep_solve(&pc2, maxit, stol, phi2w, &stw2);
    double twarm = now() - tk0;
    double dw = 0.0;
    for (int32_t i = 0; i < n2 * 4; i++) {
      double d = fabs(phi2w[i] - phi2c[i]);
      if (d > dw) dw = d;
    }
    printf("ЭТАП A ШАГ 3, ТЁПЛЫЙ СТАРТ ЧЕРЕЗ ПЕРЕНОС (сетка под сдвинутой камерой):\n");
    printf("     сеток: старая %d ячеек, новая %d ячеек; построение новой %.2f с; перенос %.3f с\n",
           mesh->ncell, n2, tbuild, txfer);
    printf("     ХОЛОДНЫЙ: код %d, итераций %d, %.2f с\n", rcc2, stc2.iters, tcool);
    printf("     ТЁПЛЫЙ:   код %d, итераций %d, %.2f с  (×%.2f)\n", rcw2, stw2.iters, twarm,
           tcool > 0.0 ? tcool / twarm : 0.0);
    printf("     начальная ошибка переноса e0 = %.3e при |φ| = %.3e;\n", e0, pmax2);
    printf("     разряды e0: у поверхности max %.3e, сред %.3e (яч. %ld); в объёме max %.3e, "
           "сред %.3e (яч. %ld)\n",
           e0_surf, n_surf ? sum_surf / (double)n_surf : 0.0, n_surf, e0_vol,
           n_vol ? sum_vol / (double)n_vol : 0.0, n_vol);
    printf("     СОВПАДЕНИЕ ОТВЕТОВ (НК, как К76): max|Δφ| %.3e\n", dw);
    free(stc2.bout);
    free(stc2.sout);
    free(stc2.eirr);
    free(stw2.bout);
    free(stw2.sout);
    free(stw2.eirr);
    free(sig2s);
    free(sig2t);
    free(phi2w);
    free(phi2c);
    world_free(&w2);
  }

  /* --- ЭТАП A, ШАГ 3 (П5-КАДР, xwarm = 2): ЦЕПОЧКА КАДРОВ И СТАБИЛЬНОСТЬ ---
   *
   * Кадровый фальсификатор движения: K позиций вдоль дуги камеры; на КАЖДОЙ
   * строится своё дерево (`L = εR` от текущего глаза), хранимое переносится
   * с предыдущей сетки и тёплым стартом доводится до сходимости. В КОНЦЕ
   * траектории задача решается ещё раз ХОЛОДНЫМ стартом, и оба кадра
   * собираются из одной точки. Картинки ОБЯЗАНЫ совпасть в допуске
   * дискретизации вне кромок (К20) — расхождение сверх этого значило бы, что
   * цепочка переносов уводит решение, и движение камеры незаконно. */
  if (xwarm == 2 || xwarm == 3) {
    const int K = 4;
    /* миры — из кучи, обмен УКАЗАТЕЛЯМИ (копия структуры делает mesh.tree
     * висячим адресом &wnext.t — найдено этим прогоном) */
    world *wprev = malloc(sizeof(world)), *wnext = malloc(sizeof(world));
    if (!wprev || !wnext) return orbit_fail(__LINE__);
    {
      double leye0[3] = {8.0, 0.6, 7.2};
      if (world_build(wprev, log2n, nc, &fr, &stab, &ftab, f0, nfac, NB, lodeps, leye0, cutfine))
        return orbit_fail(__LINE__);
    }
    double *last_bout = NULL, *last_sout = NULL, *last_eirr = NULL;
    double *phicur = calloc((size_t)wprev->mesh.ncell * 4, sizeof(double));
    double *sigcur_t = calloc((size_t)wprev->mesh.ncell, sizeof(double));
    double *sigcur_s = calloc((size_t)wprev->mesh.ncell, sizeof(double));
    if (!phicur || !sigcur_t || !sigcur_s) return orbit_fail(__LINE__);
    for (int32_t c = 0; c < wprev->mesh.ncell; c++) {
      sigcur_t[c] = 0.015;
      sigcur_s[c] = 0.012;
    }
    tr3_problem pc0 = {.m = &wprev->mesh,
                       .d = &dirs,
                       .cut = &wprev->cut,
                       .facet_rho = frho,
                       .facet_emit = femit,
                       .nfacet = ftab.n,
                       .sig_t = sigcur_t,
                       .sig_s = sigcur_s,
                       .wall_rho = wr,
                       .wall_emit = we,
                       .limiter = limiter,
                       .trace = 0};
    tr3_stats st0;
    if (tr3_sweep_solve(&pc0, maxit, stol, phicur, &st0) != 0) return orbit_fail(__LINE__);
    /* Г1: наборы коробок (lo, size) всех миров цепочки */
    enum { MAXB = 40000, MAXSTEP = 4 };
    static int32_t boxes[MAXSTEP + 1][MAXB][4];
    static int nbox[MAXSTEP + 1];
    {
      world *ws[4] = {NULL, NULL, NULL, NULL};
      (void)ws;
    }
    nbox[0] = 0;
    for (int32_t c = 0; c < wprev->mesh.ncell && nbox[0] < MAXB; c++) {
      boxes[0][nbox[0]][0] = wprev->mesh.clo[c][0];
      boxes[0][nbox[0]][1] = wprev->mesh.clo[c][1];
      boxes[0][nbox[0]][2] = wprev->mesh.clo[c][2];
      boxes[0][nbox[0]][3] = wprev->mesh.csize[c];
      nbox[0]++;
    }
    free(st0.bout);
    free(st0.sout);
    free(st0.eirr);
    for (int k = 1; k <= K; k++) {
      double a2 = 0.35 * (double)k / (double)K;
      double eye2[3] = {8.0 + 6.0 * sin(a2), 0.6 - 6.0 * (1.0 - cos(a2)), 7.2};
      if (world_build(wnext, log2n, nc, &fr, &stab, &ftab, f0, nfac, NB, lodeps, eye2, cutfine))
        return orbit_fail(__LINE__);
      double *phinext = calloc((size_t)wnext->mesh.ncell * 4, sizeof(double));
      if (!phinext) return orbit_fail(__LINE__);
      if (tr3_xfer(&wnext->mesh, phinext, &wprev->mesh, phicur)) return orbit_fail(__LINE__);
      double *sgt = calloc((size_t)wnext->mesh.ncell, sizeof(double));
      double *sgs = calloc((size_t)wnext->mesh.ncell, sizeof(double));
      if (!sgt || !sgs) return orbit_fail(__LINE__);
      for (int32_t c = 0; c < wnext->mesh.ncell; c++) {
        sgt[c] = 0.015;
        sgs[c] = 0.012;
      }
      tr3_problem pcn = {.m = &wnext->mesh,
                         .d = &dirs,
                         .cut = &wnext->cut,
                         .facet_rho = frho,
                         .facet_emit = femit,
                         .nfacet = ftab.n,
                         .sig_t = sgt,
                         .sig_s = sgs,
                         .wall_rho = wr,
                         .wall_emit = we,
                         .limiter = limiter,
                         .warm_start = 1,
                         .trace = 0};
      tr3_stats stn2;
      double tk = now();
      if (tr3_sweep_solve(&pcn, maxit, stol, phinext, &stn2) != 0) return orbit_fail(__LINE__);
      long nsurv = 0;
      for (int32_t c = 0; c < wnext->mesh.ncell; c++) {
        int32_t rlo[3], rsz = 0;
        int32_t ni = hz_oct_leaf_box(&wprev->t, wnext->mesh.clo[c][0], wnext->mesh.clo[c][1],
                                     wnext->mesh.clo[c][2], rlo, &rsz);
        if (ni >= 0 && rsz == wnext->mesh.csize[c] && rlo[0] == wnext->mesh.clo[c][0] &&
            rlo[1] == wnext->mesh.clo[c][1] && rlo[2] == wnext->mesh.clo[c][2])
          nsurv++;
      }
      printf("  орбита %d/%d: ячеек %d (пережило %ld, %.1f%%), тёплых итераций %d, %.2f с\n", k, K,
             wnext->mesh.ncell, nsurv, 100.0 * (double)nsurv / (double)wnext->mesh.ncell,
             stn2.iters, now() - tk);
      free(last_bout);
      free(last_sout);
      free(last_eirr);
      nbox[k] = 0;
      for (int32_t c = 0; c < wnext->mesh.ncell && nbox[k] < MAXB; c++) {
        boxes[k][nbox[k]][0] = wnext->mesh.clo[c][0];
        boxes[k][nbox[k]][1] = wnext->mesh.clo[c][1];
        boxes[k][nbox[k]][2] = wnext->mesh.clo[c][2];
        boxes[k][nbox[k]][3] = wnext->mesh.csize[c];
        nbox[k]++;
      }
      last_bout = stn2.bout;
      last_sout = stn2.sout;
      last_eirr = stn2.eirr; /* сошедшееся хранимое последнего звена — кадру */
      /* кадр на этой позиции собирается из ТЁПЛОГО поля цепочки */
      free(sigcur_t);
      free(sigcur_s);
      free(phicur);
      phicur = phinext;
      sigcur_t = sgt;
      sigcur_s = sgs;
      world_free(wprev);
      {
        world *tmp = wprev;
        wprev = wnext;
        wnext = tmp; /* буфер прежнего мира переиспользуется, а не копируется */
      }
    }
    /* ХОЛОДНЫЙ контроль в КОНЕЧНОЙ точке: та же сетка, нулевое начало */
    double *phicold = calloc((size_t)wprev->mesh.ncell * 4, sizeof(double));
    if (!phicold) return orbit_fail(__LINE__);
    tr3_problem pcc = {.m = &wprev->mesh,
                       .d = &dirs,
                       .cut = &wprev->cut,
                       .facet_rho = frho,
                       .facet_emit = femit,
                       .nfacet = ftab.n,
                       .sig_t = sigcur_t,
                       .sig_s = sigcur_s,
                       .wall_rho = wr,
                       .wall_emit = we,
                       .limiter = limiter,
                       .trace = 0};
    tr3_stats stc3;
    double tcc = now();
    if (tr3_sweep_solve(&pcc, maxit, stol, phicold, &stc3) != 0) return orbit_fail(__LINE__);
    double dcold = 0.0, pm3 = 0.0;
    for (int32_t i = 0; i < wprev->mesh.ncell * 4; i++) {
      double d = fabs(phicur[i] - phicold[i]);
      if (d > dcold) dcold = d;
      if (fabs(phicold[i]) > pm3) pm3 = fabs(phicold[i]);
    }
    printf("ШАГ 3, ЦЕПОЧКА КАДРОВ: цепочка из %d переякорений, холодный в конце %d итераций, "
           "%.2f с\n",
           K, stc3.iters, now() - tcc);
    { /* Г1: осцилляция — коробка ЕСТЬ в мире k, НЕТ в k+1, СНОВА ЕСТЬ в k+2 */
      long osc = 0, gone = 0;
      for (int k = 0; k + 2 <= K; k++) {
        for (int i = 0; i < nbox[k]; i++) {
          int in_k1 = 0, in_k2 = 0;
          for (int q = 0; q < nbox[k + 1] && !in_k1; q++)
            if (boxes[k + 1][q][0] == boxes[k][i][0] && boxes[k + 1][q][1] == boxes[k][i][1] &&
                boxes[k + 1][q][2] == boxes[k][i][2] && boxes[k + 1][q][3] == boxes[k][i][3])
              in_k1 = 1;
          if (in_k1) continue;
          gone++;
          for (int q = 0; q < nbox[k + 2] && !in_k2; q++)
            if (boxes[k + 2][q][0] == boxes[k][i][0] && boxes[k + 2][q][1] == boxes[k][i][1] &&
                boxes[k + 2][q][2] == boxes[k][i][2] && boxes[k + 2][q][3] == boxes[k][i][3])
              in_k2 = 1;
          if (in_k2) osc++;
        }
      }
      printf("     Г1: коробок, исчезнувших на шаге, %ld; из них ВЕРНУЛОСЬ через шаг (осцилляция) "
             "%ld (%.1f%%)\n",
             gone, osc, gone ? 100.0 * (double)osc / (double)gone : 0.0);
    }
    printf("     φ: max|Δ| цепочка против холодного %.3e при |φ|max %.3e\n", dcold, pm3);
    /* КАДРЫ из одной точки: тёплый (цепочка) против холодного */
    {
      double aK = 0.35, eyeK[3] = {8.0 + 6.0 * sin(aK), 0.6 - 6.0 * (1.0 - cos(aK)), 7.2};
      tr3_scene scn2 = {.tree = &wprev->t,
                        .fr = fr,
                        .sigma = NULL,
                        .st = &stab,
                        .ft = &ftab,
                        .cm = &wprev->cmap,
                        .facet_only = facet_only};
      tr3_camera cam2;
      if (tr3_camera_look(&cam2, eyeK, at, up, 1.3, W, H)) return 1;
      int32_t *wi2 = calloc((size_t)6 * (size_t)nc * (size_t)nc, sizeof(int32_t));
      double *bufA = calloc((size_t)W * (size_t)H * TR3_MAXCH, sizeof(double));
      double *bufB = calloc((size_t)W * (size_t)H * TR3_MAXCH, sizeof(double));
      if (!wi2 || !bufA || !bufB) return 1;
      for (int32_t i = 0; i < 6 * (int32_t)nc * (int32_t)nc; i++)
        wi2[i] = -1;
      for (int32_t f = 0; f < wprev->mesh.nf; f++) {
        const tr3_face *ff = &wprev->mesh.f[f];
        if (ff->cb >= 0) continue;
        int wl = (int)(~ff->cb);
        for (int32_t a2 = ff->lo[0]; a2 < ff->hi[0]; a2++)
          for (int32_t b2 = ff->lo[1]; b2 < ff->hi[1]; b2++)
            wi2[((int32_t)wl * nc + a2) * nc + b2] = f;
      }
      /* хранимое поверхностей — СОШЕДШЕЕСЯ в обоих случаях: у тёплой цепочки
       * это финальный такт последнего звена (last_*), у холодного — stc3.
       * Первая редакция гоняла один проход развёртки для «тёплого» хранимого —
       * кадры расходились на 53% именно поэтому: дефект был в конструкции
       * теста, а не в переносе. */
      double *bouts = last_bout, *souts = last_sout;
      double *bouts2 = stc3.bout, *souts2 = stc3.sout;
      tr3_gather gA = {.sc = &scn2,
                       .m = &wprev->mesh,
                       .cut = &wprev->cut,
                       .bout = bouts,
                       .sout = souts,
                       .nch = nch,
                       .facet_spec = spec > 0.0 ? fspec : NULL,
                       .nfacet = ftab.n,
                       .wallidx = wi2,
                       .nwall = nc,
                       .maxbounce = maxbounce,
                       .sig_t = sigcur_t,
                       .sig_s = sigcur_s,
                       .phi = phicur};
      tr3_gather gB = {.sc = &scn2,
                       .m = &wprev->mesh,
                       .cut = &wprev->cut,
                       .bout = bouts2,
                       .sout = souts2,
                       .nch = nch,
                       .facet_spec = spec > 0.0 ? fspec : NULL,
                       .nfacet = ftab.n,
                       .wallidx = wi2,
                       .nwall = nc,
                       .maxbounce = maxbounce,
                       .sig_t = sigcur_t,
                       .sig_s = sigcur_s,
                       .phi = phicold};
      double dsum = 0.0, dmx = 0.0, ssum = 0.0;
      for (int py = 0; py < H; py++)
        for (int px = 0; px < W; px++) {
          double o[3], d[3], va = 0, vb = 0;
          tr3_camera_ray(&cam2, px, py, o, d);
          va = tr3_gather_ray(&gA, o, d, NULL, NULL);
          vb = tr3_gather_ray(&gB, o, d, NULL, NULL);
          dsum += fabs(va - vb);
          if (fabs(va - vb) > dmx) dmx = fabs(va - vb);
          ssum += vb;
        }
      printf("     КАДР из конечной точки: |Δ| сред %.3e (%.3f%% среднего), макс %.3e\n",
             dsum / (double)(W * H), 100.0 * dsum / (ssum > 0.0 ? ssum : 1.0), dmx);
      /* приёмка: среднее расхождение обязано сидеть в допуске дискретизации
       * (0.5% среднего радианса — тот же порядок, что сверка проекции с
       * маршем вне кромок, 0.20%); провал значил бы, что цепочка переносов
       * уводит решение. */
      if (!(dsum / (double)(W * H) <= 0.005 * (ssum / (double)(W * H)))) {
        fprintf(stderr, "П5-КАДР: картинки цепочки и холодного старта РАСШОДИТСЯ сверх допуска\n");
        return orbit_fail(__LINE__);
      }
      free(stc3.eirr);
      free(bufB);
      free(bufA);
      free(wi2);
    }
    free(phicold);
    free(last_bout);
    free(last_sout);
    free(last_eirr);
    free(stc3.bout);
    free(stc3.sout);
    free(sigcur_s);
    free(sigcur_t);
    free(phicur);
    world_free(wprev);
    free(wprev);
    free(wnext);
  }

  char pfm[512];
  size_t ln = strlen(out);
  if (ln + 5 < sizeof pfm) {
    memcpy(pfm, out, ln + 1);
    char *dot = strrchr(pfm, '.');
    size_t ext = dot != NULL ? (size_t)(dot - pfm) : ln;
    memcpy(pfm + ext, ".pfm", 5);
    /* R = G = B: спектральный канал ОДИН, и файл об этом не врёт */
    /* НАСТОЯЩИЙ RGB, если каналов три; при одном — `R = G = B`, и файл об этом
     * не врёт (формат взят трёхканальным сознательно, см. доклад 07-27). */
    const double *pg = nch > 1 ? buf + (size_t)W * (size_t)H : buf;
    const double *pb = nch > 2 ? buf + (size_t)2 * (size_t)W * (size_t)H : buf;
    if (hz_pfm_write(pfm, buf, pg, pb, W, H) != 0) {
      fprintf(stderr, "не записалось: %s\n", pfm);
      return 1;
    }
  } else {
    pfm[0] = '\0';
  }
  if (hz_ppm_write(out, buf, W, H) != 0) {
    fprintf(stderr, "не записалось: %s\n", out);
    return 1;
  }
  if (spec > 0.0)
    printf("ЗЕРКАЛА: доля %.2f, предел отскоков %d, отскоков всего %ld (до %d на луч)\n", spec,
           maxbounce, nbtot, nbmax);
  double lo = 1e300, hi2 = -1e300;
  for (size_t i = 0; i < (size_t)W * (size_t)H; i++) {
    if (buf[i] < lo) lo = buf[i];
    if (buf[i] > hi2) hi2 = buf[i];
  }
  printf("картинки: %s — РАДИАНС как есть, float RGB (R=G=B, канал один);\n"
         "          %s — тон-маппинг для глаз, сверять по нему НЕЛЬЗЯ (К19)\n"
         "          %dx%d, радианс от %.6f до %.6f\n",
         pfm[0] ? pfm : "(не записан)", out, W, H, lo, hi2);
  free(wallidx);
  free(w1.solid);
  free(buf);
  free(phi);
  free(st.bout);
  free(st.sout);
  free(st.eirr); /* §744: экспорт облучённости — живёт в stats, freeing на вызывающем */
  free(boutc);
  free(soutc);
  free(sig_t);
  free(sig_s);
  free(frho);
  free(fspec);
  free(femit);
  tr3_dirs_free(&dirs);
  tr3_cut_free(cut);
  tr3_mesh_free(mesh);
  hz_cutmap_free(cmap);
  hz_facettab_free(&ftab);
  hz_surftab_free(&stab);
  hz_oct_free(&w1.t);
  /* §829: пассивная печать приборов карты/марша — НК-б плана §829 (зубы
   * прибора: здесь марш по построению в каждом пикселе, счётчики обязаны быть
   * ненулевыми). Физику не читает и не меняет. */
  hz_cutmap_stats_print("render3");
  hz_ray3_stats_print("render3");
  return 0;
}
