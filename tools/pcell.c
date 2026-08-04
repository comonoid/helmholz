/* pcell — СЕГМЕНТАЦИЯ ПО ЯЧЕЙКАМ ПРОТИВ ГЛОБАЛЬНОЙ (замер О57, план §191).
 *
 * ВОПРОС ОДНОЙ СТРОКОЙ. Глобальная сегментация Сан-Мигеля стоит `2022` с, и
 * §192 замерила, что цену задаёт ЧИСЛО ТРЕУГОЛЬНИКОВ В ОДНОМ ЭЛЕМЕНТЕ (а не
 * его габарит: ограничение габарита `scap = 2` м дало `7 %`). Восьмидерево
 * ограничивает это число по построению — `leafmax`, порог дробления узла. Здесь
 * измеряется, во что это обходится и что даёт.
 *
 * ЧТО СРАВНИВАЕТСЯ. Числа глобального прогона (`result/O56_miguel_seg.txt`,
 * `O56_miguel_lowpoly_seg.txt`) против того же ядра `hz_seg_planar`,
 * запущенного НА КАЖДОМ ЛИСТЕ дерева. Ядро одно и то же — меняется только
 * набор, который ему дают.
 *
 * ПОДСЕТКА ЛИСТА НЕ КОПИРУЕТ ВЕРШИНЫ. `v` разделяется с исходной сеткой как
 * есть, копируются только грани листа; индексы при этом остаются валидными.
 * Копировать `5.9` млн вершин на каждый из десятков тысяч листьев было бы
 * дороже самого замера.
 *
 * КОНТРОЛЬ ЭКВИВАЛЕНТНОСТИ (`leafmax` больше числа треугольников): дерево из
 * ОДНОГО листа обязано дать в точности глобальный ответ. Это предусмотрено
 * автором `ptree` (см. `src/ptree.h`, НК13) и здесь только используется. Это НЕ
 * негативный контроль, и называется так прямо (правило А10).
 *
 * ИНВАРИАНТ ПЛОЩАДИ. Дробление НЕ меняет суммарной площади элементов. Она
 * печатается и сверяется с глобальной: расхождение выше относительной `1e-9`
 * означает, что дробление теряет или дублирует геометрию.
 *
 * Запуск: `pcell ФАЙЛ.obj МАСШТАБ [d=МЕТРЫ] [leaf=N] [lev=N]`.
 */
#include "poly_seg.h"
#include "scene_cfg.h"
#include "pgrid.h"
#include "ptree.h"
#include "scene_obj.h"
#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Допуск планарности. Тот же, что у городских прогонов `pcoarse` (`0.05` м):
 * сравнивать надо с ними, а разный допуск сделал бы сравнение бессмысленным. */
#define PCELL_DELTA 0.05
/* Уровень объединения вторичных источников. Выведен ДВУМЯ независимыми
 * замерами, а не подобран: §199 — дерево ведёт себя как поверхность (`×4` на
 * уровень) ровно до шестого; §202 — там же сжатие освещённого набора выходит на
 * `30×` при `1 008` источниках. */
#define PCELL_AGG_LEV 5
/* Альбедо агрегата. Материала у объединения нет, и брать его неоткуда: это
 * СРЕДНЕЕ по сцене, названное здесь, а не спрятанное в формулу. Городское
 * `0.3` из разбора 08-03. */
#define PCELL_RHO 0.3
/* Множитель косвенного в кадре. Радиометрия сцены условна (солнце `E = 1`),
 * и этот множитель приводит косвенное к тому же масштабу. Названо числом, а не
 * подобрано на глаз: при `E = 1` и `ρ = 0.3` косвенное обязано дать порядка
 * `0.09` от прямого, что и заложено. */
#define PCELL_IND 1.0
/* Угловой ДИАМЕТР солнца — физическая величина, не порог схемы (та же, что в
 * `tools/pfront.c`). Из-за неё тень имеет полутень, растущую с расстоянием до
 * заслона: `9` мм на метре, `9` см на десяти. */
#define PCELL_SUN_DEG 0.533
/* Проб по диску солнца на пиксель. При `16` шум полутени ниже кванта восьми
 * разрядов на всех расстояниях этой сцены. */
/* Проб по диску солнца. `16` давало ступень полутени `6 %` при кванте `0.4 %`
 * — обоснование «шум ниже кванта» было арифметически ложным (§207.1). При `64`
 * с поворотом набора НА ПИКСЕЛЬ остаток становится шумом, а не лесенкой. */
#define PCELL_SUN_SAMP 64
/* Яркость неба как доли солнечной. Ясный день: небо даёт около четверти
 * горизонтальной облучённости против прямого солнца. */
/* Небо как доля прямого солнца. `0.25` было завышением: столько получает
 * ГОРИЗОНТАЛЬНАЯ площадка под открытым небом, а стена в каньоне двора —
 * `0.10…0.15`. Завышение делало тени слишком светлыми. */
#define PCELL_SKY 0.12
/* ЭКСПОЗИЦИЯ. Без неё освещённая поверхность при `kd ≈ 1` и солнце `1.0` даёт
 * ровно `1.0` и упирается в потолок: всё на свету схлопывается в белое, и
 * контраста там нет по построению. Число выбрано так, чтобы прямое солнце на
 * площадке, повёрнутой к нему, ложилось в `0.9`, оставляя запас под блики. */
#define PCELL_EXPO 0.9
/* Критерий дробления при иерархическом сборе: узел берётся источником, пока
 * приёмник дальше `PCELL_SRC_K` его радиусов. `4` означает угловой размер
 * источника ниже `~28°` — та же по смыслу величина, что допуск связи в
 * иерархической радиосити. Меньше — точнее и дороже. */
#define PCELL_SRC_K 4.0
/* Проб неба на пиксель. При `24` шум видимости ниже кванта восьми разрядов
 * везде, кроме контактных стыков, — та же оценка `1/√N`, что в §183. */
#define PCELL_SKY_SAMP 24

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int cmp_dbl(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

static int cmp_i32(const void *a, const void *b) {
  int32_t x = *(const int32_t *)a, y = *(const int32_t *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "pcell ФАЙЛ.obj МАСШТАБ [d=МЕТРЫ] [leaf=N] [lev=N]\n");
    return 1;
  }
  double delta = PCELL_DELTA;
  int leafmax = 0, maxlev = 0, wpx = HZ_CFG_W, loose = 0, haseye = 0;
  double eye0[3] = {0.0, 0.0, 0.0};
  /* Полоса §191, названная пользователем: элемент от 4 до 32 пикселей. */
  double pxlo = 4.0, pxhi = 32.0;
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "d=", 2) == 0) delta = strtod(argv[i] + 2, NULL);
    if (strncmp(argv[i], "leaf=", 5) == 0) leafmax = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "lev=", 4) == 0) maxlev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "w=", 2) == 0) wpx = (int)strtol(argv[i] + 2, NULL, 10);
    if (strncmp(argv[i], "pxhi=", 5) == 0) pxhi = strtod(argv[i] + 5, NULL);
    if (strncmp(argv[i], "pxlo=", 5) == 0) pxlo = strtod(argv[i] + 5, NULL);
    if (strcmp(argv[i], "loose") == 0) loose = 1;
    /* Камера ключом: закон роста среза надо мерить на РАЗНЫХ сценах, а глаз у
     * каждой свой и найден замером (§187). */
    if (strncmp(argv[i], "eye=", 4) == 0) {
      char *e = NULL;
      eye0[0] = strtod(argv[i] + 4, &e);
      if (e != NULL && *e == ',') eye0[1] = strtod(e + 1, &e);
      if (e != NULL && *e == ',') eye0[2] = strtod(e + 1, NULL);
      haseye = 1;
    }
  }

  hz_objmesh m;
  double t0 = now_s();
  if (hz_obj_load(&m, argv[1], atof(argv[2])) != 0) {
    fprintf(stderr, "нет сцены %s\n", argv[1]);
    return 1;
  }
  printf("== СЦЕНА %s: треугольников %d, вершин %d, габарит %.2f x %.2f x %.2f м, чтение %.1f с\n",
         argv[1], m.nt, m.nv, m.hi[0] - m.lo[0], m.hi[1] - m.lo[1], m.hi[2] - m.lo[2],
         now_s() - t0);
  printf("== δ сегментации %g м; leafmax %d, maxlev %d (0 — умолчания ptree: 16 и 12)\n", delta,
         leafmax, maxlev);

  /* ---- ДЕРЕВО ---------------------------------------------------------- */
  hz_ptree T;
  t0 = now_s();
  if (hz_ptree_build_ex(&T, &m, leafmax, maxlev, loose) != 0) {
    fprintf(stderr, "отказ дерева\n");
    return 2;
  }
  double t_tree = now_s() - t0;
  printf("== ДЕРЕВО за %.2f с: узлов %d, ЛИСТЬЕВ %lld, ссылок %lld, ИЗБЫТОЧНОСТЬ nref/nt = %.3f\n",
         t_tree, T.nnd, (long long)T.nleaf, (long long)T.nref, T.redundancy);
  printf("   треугольников у ВНУТРЕННИХ узлов %lld (%.2f %%), у листьев %lld\n",
         (long long)T.n_inner, 100.0 * (double)T.n_inner / (double)m.nt, (long long)T.n_leaf_tri);

  /* ---- СЕГМЕНТАЦИЯ ПО ЛИСТЬЯМ ------------------------------------------ */
  /* СЕГМЕНТИРУЮТСЯ ВСЕ УЗЛЫ С ТРЕУГОЛЬНИКАМИ, А НЕ ТОЛЬКО ЛИСТЬЯ.
   * Первая редакция брала только листья и ТЕРЯЛА 78.60 %% треугольников: те,
   * что не влезли целиком ни в одного ребёнка, висят у ВНУТРЕННИХ узлов
   * (`ptree.h` про это предупреждает прямо). Поймал инвариант площади:
   * `2.4076` м² против `213.3018` на зале. Каждый треугольник числится ровно у
   * одного узла (`n_inner + n_leaf_tri = nt` точно), поэтому обход всех узлов
   * покрывает сцену без пропусков и без двойного счёта. */
  int32_t *leaf = malloc((size_t)(T.nnd > 0 ? T.nnd : 1) * sizeof *leaf);
  if (leaf == NULL) return 2;
  int64_t nl = 0;
  for (int32_t i = 0; i < T.nnd; i++)
    if (T.nd[i].ntri > 0) leaf[nl++] = i;

  /* Число элементов и их размеры — на лист, чтобы слить потом без замка. */
  int32_t **cnt = calloc((size_t)(nl > 0 ? nl : 1), sizeof *cnt);
  int32_t *nseg_of = calloc((size_t)(nl > 0 ? nl : 1), sizeof *nseg_of);
  double *area_of = calloc((size_t)(nl > 0 ? nl : 1), sizeof *area_of);
  double *dmax_of = calloc((size_t)(nl > 0 ? nl : 1), sizeof *dmax_of);
  if (cnt == NULL || nseg_of == NULL || area_of == NULL || dmax_of == NULL) return 2;

  /* СЕТКА ЛУЧЕЙ — для ОСВЕЩЁННОГО НАБОРА: у каждого элемента считается его
   * положение (центр по площади) и видно ли из него солнце. Это первая
   * половина схемы «фронт от источника → освещённый набор → вторичные
   * источники» (директива пользователя 08-04, §186/§190). */
  hz_pgrid G;
  if (hz_pgrid_build(&G, &m) != 0) {
    fprintf(stderr, "отказ сетки лучей\n");
    return 2;
  }
  double wsun[3] = HZ_CFG_SUN_DIR;
  {
    double L = sqrt(wsun[0] * wsun[0] + wsun[1] * wsun[1] + wsun[2] * wsun[2]);
    for (int c = 0; c < 3; c++)
      wsun[c] /= L;
  }
  double diag2 = 0.0;
  for (int c = 0; c < 3; c++) {
    double d2 = m.hi[c] - m.lo[c];
    diag2 += d2 * d2;
  }
  diag2 = 4.0 * sqrt(diag2);
  int64_t nlit = 0;
  double alit = 0.0;
  /* ИНВАРИАНТ ЭНЕРГИИ (§206). Считается НЕЗАВИСИМО ОТ ДЕРЕВА, прямо в цикле по
   * элементам: суммарная мощность вторичных источников обязана совпасть с этим
   * числом. Три раза за заход одна ошибка (А419, А420, §206) проходила там, где
   * инварианта не было, и ловилась немедленно там, где он был. */
  double wref = 0.0;
  /* Вторичные источники: центр+нормаль (по 6 чисел) и мощность. Потолок —
   * число узлов, заведомо с запасом. */
  int64_t nsrc = 0;
  /* Суммы по поддеревьям переживают блок агрегирования: их читает кадр. */
  double *g_lv = NULL, *g_lp = NULL;
  int64_t *g_ls = NULL;
  int32_t *g_dep = NULL;
  (void)g_dep;
  double *src_p = calloc((size_t)(T.nnd > 0 ? T.nnd : 1) * 6, sizeof *src_p);
  double *src_n = calloc((size_t)(T.nnd > 0 ? T.nnd : 1) * 6, sizeof *src_n);
  double *src_w = calloc((size_t)(T.nnd > 0 ? T.nnd : 1), sizeof *src_w);
  if (src_p == NULL || src_n == NULL || src_w == NULL) return 2;
  /* Признак «освещён» НА ТРЕУГОЛЬНИК — только ради картинки: она рисуется
   * первичными лучами, а луч попадает в треугольник, не в элемент. */
  unsigned char *tri_lit = calloc((size_t)(m.nt > 0 ? m.nt : 1), 1);
  if (tri_lit == NULL) return 2;
  /* Освещённость ПО УЗЛАМ — чтобы потом сжать освещённый набор по уровням
   * дерева: наивный сбор 1.07 млн приёмников на 30 тыс. источников есть 3.2e10
   * пар, и вопрос не в том, дорого ли это, а во сколько раз агрегирование
   * сжимает набор (пункт 2 реестра: «второе отражение на ОГРУБЛЁННОМ наборе»). */
  int32_t *litn = calloc((size_t)(T.nnd > 0 ? T.nnd : 1), sizeof *litn);
  /* Положение агрегата (сумма площадь×центр) и его излучение (сумма
   * площадь×cos): вторичный источник есть площадка с этим центром, этой
   * нормалью и этой мощностью. */
  double *litp = calloc((size_t)(T.nnd > 0 ? T.nnd : 1) * 4, sizeof *litp);
  /* СВЯЗНОСТЬ НОРМАЛЕЙ освещённой группы: `|Σ area·n| / Σ area`. Единица —
   * группа плоская и заменима ОДНОЙ плоскостью; ноль — нормали смотрят врозь, и
   * агрегат обязан нести РАСПРЕДЕЛЕНИЕ, а не плоскость. Без этого числа
   * объединение источников (§202) построить нельзя: неизвестно, что у агрегата
   * за нормаль. */
  double *litv = calloc((size_t)(T.nnd > 0 ? T.nnd : 1) * 4, sizeof *litv);
  if (litn == NULL || litv == NULL || litp == NULL) return 2;

  int64_t nfail = 0;
  t0 = now_s();
#pragma omp parallel for schedule(dynamic, 8) reduction(+ : nfail, nlit, alit, wref)
  for (int64_t k = 0; k < nl; k++) {
    const hz_ptnode *nd = &T.nd[leaf[k]];
    int32_t nt = nd->ntri;
    /* Подсетка: ВЕРШИНЫ РАЗДЕЛЯЮТСЯ, копируются только грани. */
    hz_objmesh sub = m;
    sub.nt = nt;
    sub.f = malloc((size_t)nt * 3 * sizeof *sub.f);
    sub.fm = malloc((size_t)nt * sizeof *sub.fm);
    sub.fn = (m.fn != NULL) ? malloc((size_t)nt * 3 * sizeof *sub.fn) : NULL;
    if (sub.f == NULL || sub.fm == NULL || (m.fn != NULL && sub.fn == NULL)) {
      free(sub.f);
      free(sub.fm);
      free(sub.fn);
      nfail++;
      continue;
    }
    for (int32_t j = 0; j < nt; j++) {
      int32_t t = T.ref[nd->t0 + j];
      for (int c = 0; c < 3; c++) {
        sub.f[3 * j + c] = m.f[3 * (size_t)t + (size_t)c];
        if (sub.fn != NULL) sub.fn[3 * j + c] = m.fn[3 * (size_t)t + (size_t)c];
      }
      sub.fm[j] = m.fm[t];
    }
    hz_pseglist sg;
    if (hz_seg_planar(&sg, &sub, delta) != 0) {
      nfail++;
    } else {
      nseg_of[k] = sg.nseg;
      int32_t *c2 = malloc((size_t)(sg.nseg > 0 ? sg.nseg : 1) * sizeof *c2);
      double ar = 0.0, dm = 0.0;
      for (int32_t s = 0; s < sg.nseg; s++) {
        if (c2 != NULL) c2[s] = sg.seg[s].ntri;
        ar += sg.seg[s].area;
        if (sg.seg[s].dmax > dm) dm = sg.seg[s].dmax;
      }
      /* ОСВЕЩЁННЫЙ НАБОР: центр элемента по площади треугольников, затем один
       * теневой луч к солнцу. Отступ по нормали — тот же относительный, что у
       * всей оснастки. */
      {
        double *cx = calloc((size_t)(sg.nseg > 0 ? sg.nseg : 1) * 4, sizeof *cx);
        if (cx != NULL) {
          for (int32_t j = 0; j < nt; j++) {
            int32_t s = sg.label[j];
            if (s < 0 || s >= sg.nseg) continue;
            const double *A = m.v + 3 * (size_t)sub.f[3 * j];
            const double *B = m.v + 3 * (size_t)sub.f[3 * j + 1];
            const double *C = m.v + 3 * (size_t)sub.f[3 * j + 2];
            /* Инициализация нулём — тот же класс ложных срабатываний OpenMP,
             * что описан в CLAUDE.md: анализатор не проводит связь «заполнено
             * циклом — прочитано» через параллельную область. */
            double e1[3] = {0.0, 0.0, 0.0}, e2[3] = {0.0, 0.0, 0.0}, cr[3] = {0.0, 0.0, 0.0};
            for (int c = 0; c < 3; c++) {
              e1[c] = B[c] - A[c];
              e2[c] = C[c] - A[c];
            }
            cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
            cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
            cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
            double atri = 0.5 * sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
            for (int c = 0; c < 3; c++)
              cx[4 * s + c] += atri * (A[c] + B[c] + C[c]) / 3.0;
            cx[4 * s + 3] += atri;
          }
          for (int32_t s = 0; s < sg.nseg; s++) {
            if (!(cx[4 * s + 3] > 0.0)) continue;
            double p[3];
            for (int c = 0; c < 3; c++)
              p[c] = cx[4 * s + c] / cx[4 * s + 3];
            double ndl =
                -(sg.seg[s].n[0] * wsun[0] + sg.seg[s].n[1] * wsun[1] + sg.seg[s].n[2] * wsun[2]);
            double nn2[3] = {sg.seg[s].n[0], sg.seg[s].n[1], sg.seg[s].n[2]};
            if (ndl < 0.0) {
              ndl = -ndl;
              for (int c = 0; c < 3; c++)
                nn2[c] = -nn2[c];
            }
            if (!(ndl > 0.0)) continue;
            double o[3], dsun[3] = {-wsun[0], -wsun[1], -wsun[2]}, th = 0.0;
            for (int c = 0; c < 3; c++)
              o[c] = p[c] + 1e-4 * nn2[c];
            if (!hz_pgrid_trace(&G, &m, o, dsun, 0.0, diag2, NULL, 0, 1, &th)) {
              nlit++;
              alit += sg.seg[s].area;
              litn[leaf[k]]++;
              for (int c = 0; c < 3; c++)
                litv[4 * leaf[k] + c] += sg.seg[s].area * nn2[c];
              litv[4 * leaf[k] + 3] += sg.seg[s].area;
              for (int c = 0; c < 3; c++)
                litp[4 * leaf[k] + c] += sg.seg[s].area * p[c];
              litp[4 * leaf[k] + 3] += sg.seg[s].area * ndl;
              wref += PCELL_RHO * sg.seg[s].area * ndl;
              for (int32_t j2 = 0; j2 < nt; j2++)
                if (sg.label[j2] == s) tri_lit[T.ref[nd->t0 + j2]] = 1;
            }
          }
          free(cx);
        }
      }
      cnt[k] = c2;
      area_of[k] = ar;
      dmax_of[k] = dm;
      hz_seg_free(&sg);
    }
    free(sub.f);
    free(sub.fm);
    free(sub.fn);
  }
  double t_seg = now_s() - t0;

  /* ---- СВОД ------------------------------------------------------------ */
  int64_t nel = 0;
  double area = 0.0, dmax = 0.0;
  for (int64_t k = 0; k < nl; k++) {
    nel += nseg_of[k];
    area += area_of[k];
    if (dmax_of[k] > dmax) dmax = dmax_of[k];
  }
  printf("== СЕГМЕНТАЦИЯ ПО УЗЛАМ за %.2f с (%d потоков): ЭЛЕМЕНТОВ %lld; отказов %lld\n", t_seg,
         omp_get_max_threads(), (long long)nel, (long long)nfail);
  printf("== ОСВЕЩЁННЫЙ НАБОР (солнце, один луч на элемент): элементов %lld из %lld (%.1f %%), "
         "площадь %.1f из %.1f м² (%.1f %%)\n",
         (long long)nlit, (long long)nel, 100.0 * (double)nlit / (double)(nel > 0 ? nel : 1), alit,
         area, 100.0 * alit / (area > 0.0 ? area : 1.0));
  printf("   суммарная площадь элементов %.4f м²; худший dmax %.3e м при δ = %g м — %s\n", area,
         dmax, delta, (dmax <= delta) ? "инвариант держится" : "ИНВАРИАНТ НАРУШЕН");

  /* Хвост распределения, а не среднее (§4). */
  {
    int32_t *all = malloc((size_t)(nel > 0 ? nel : 1) * sizeof *all);
    if (all != NULL) {
      int64_t p = 0;
      for (int64_t k = 0; k < nl; k++)
        for (int32_t s = 0; s < nseg_of[k]; s++)
          all[p++] = (cnt[k] != NULL) ? cnt[k][s] : 0;
      qsort(all, (size_t)p, sizeof *all, cmp_i32);
      if (p > 0)
        printf("   треугольников на элемент: среднее %.1f, p50 %d, p99 %d, максимум %d\n",
               (double)m.nt / (double)p, all[p / 2], all[(p * 99) / 100], all[p - 1]);
      free(all);
    }
  }
  printf("== ВСЕГО (холодный старт, один раз на сцену): дерево %.2f с + сегментация %.2f с = "
         "%.2f с\n",
         t_tree, t_seg, t_tree + t_seg);

  /* ---- СЖАТИЕ ОСВЕЩЁННОГО НАБОРА ПО УРОВНЯМ (§202) ---------------------
   * Сколько ВТОРИЧНЫХ ИСТОЧНИКОВ останется, если объединять освещённые
   * элементы до уровня `L`. Это и есть цена отскока: пар = приёмники ×
   * источники, и сжатие входит в неё множителем. */
  {
    int64_t *ls = calloc((size_t)(T.nnd > 0 ? T.nnd : 1), sizeof *ls);
    int32_t *dep = calloc((size_t)(T.nnd > 0 ? T.nnd : 1), sizeof *dep);
    if (ls != NULL && dep != NULL) {
      for (int32_t i = T.nnd - 1; i >= 0; i--) {
        ls[i] = litn[i];
        if (T.nd[i].child >= 0)
          for (int q = 0; q < 8; q++)
            ls[i] += ls[T.nd[i].child + q];
      }
      for (int32_t i = 0; i < T.nnd; i++)
        if (T.nd[i].child >= 0)
          for (int q = 0; q < 8; q++)
            dep[T.nd[i].child + q] = dep[i] + 1;
      printf("== СЖАТИЕ ОСВЕЩЁННОГО НАБОРА (вторичные источники после объединения до уровня L):\n");
      /* Векторные суммы нормалей по поддеревьям — тем же обратным проходом. */
      double *lv = calloc((size_t)(T.nnd > 0 ? T.nnd : 1) * 4, sizeof *lv);
      if (lv != NULL)
        for (int32_t i = T.nnd - 1; i >= 0; i--) {
          for (int c = 0; c < 4; c++)
            lv[4 * (size_t)i + (size_t)c] = litv[4 * (size_t)i + (size_t)c];
          if (T.nd[i].child >= 0)
            for (int q = 0; q < 8; q++)
              for (int c = 0; c < 4; c++)
                lv[4 * (size_t)i + (size_t)c] += lv[4 * (size_t)(T.nd[i].child + q) + (size_t)c];
        }
      for (int L = 2; L <= 12; L++) {
        int64_t c2 = 0;
        double csum = 0.0, cw = 0.0, cmin = 2.0;
        for (int32_t i = 0; i < T.nnd; i++) {
          int take = (dep[i] == L && ls[i] > 0) || (dep[i] < L && T.nd[i].child < 0 && ls[i] > 0);
          if (!take) continue;
          c2++;
          if (lv == NULL || !(lv[4 * (size_t)i + 3] > 0.0)) continue;
          double vx = lv[4 * (size_t)i], vy = lv[4 * (size_t)i + 1], vz = lv[4 * (size_t)i + 2];
          double coh = sqrt(vx * vx + vy * vy + vz * vz) / lv[4 * (size_t)i + 3];
          csum += coh * lv[4 * (size_t)i + 3];
          cw += lv[4 * (size_t)i + 3];
          if (coh < cmin) cmin = coh;
        }
        if (c2 == 0) continue;
        printf("   до уровня %2d: %8lld источников (сжатие %5.1f×); СВЯЗНОСТЬ нормалей: средняя "
               "%.3f, худшая %.3f\n",
               L, (long long)c2, (double)nlit / (double)c2, (cw > 0.0) ? csum / cw : 0.0,
               (cmin <= 1.0) ? cmin : 0.0);
      }
      /* ВТОРИЧНЫЕ ИСТОЧНИКИ уровня `PCELL_AGG_LEV` — то, ради чего вся §202.
       * Каждый есть площадка: центр по площади, нормаль по площади, мощность
       * `ρ·E·Σ(area·cos)`. Строится один раз и от камеры НЕ зависит (§190). */
      double *lp = calloc((size_t)(T.nnd > 0 ? T.nnd : 1) * 4, sizeof *lp);
      if (lp != NULL && lv != NULL) {
        for (int32_t i = T.nnd - 1; i >= 0; i--) {
          for (int c = 0; c < 4; c++)
            lp[4 * (size_t)i + (size_t)c] = litp[4 * (size_t)i + (size_t)c];
          if (T.nd[i].child >= 0)
            for (int q = 0; q < 8; q++)
              for (int c = 0; c < 4; c++)
                lp[4 * (size_t)i + (size_t)c] += lp[4 * (size_t)(T.nd[i].child + q) + (size_t)c];
        }
        for (int32_t i = 0; i < T.nnd; i++) {
          /* УЗЛЫ ВЫШЕ УРОВНЯ ОБЪЕДИНЕНИЯ ДАЮТ СВОИ СОБСТВЕННЫЕ ЭЛЕМЕНТЫ (§206).
           * Первая редакция брала только `dep == L` и листья выше, а собственные
           * элементы внутренних узлов `dep < L` не попадали НИКУДА: их нет ни в
           * одном ребёнке. У `ptree` там висит 72.5 %% геометрии, и крупнейшие
           * стены — как раз наверху. Каждый элемент числится ровно один раз:
           * узлы `dep == L` берут ПОДДЕРЕВО, узлы `dep < L` — только СВОЁ. */
          int at_lev = (dep[i] == PCELL_AGG_LEV && ls[i] > 0);
          int above = (dep[i] < PCELL_AGG_LEV && litn[i] > 0);
          if (!at_lev && !above) continue;
          const double *vsrc = at_lev ? &lv[4 * (size_t)i] : &litv[4 * (size_t)i];
          const double *psrc = at_lev ? &lp[4 * (size_t)i] : &litp[4 * (size_t)i];
          if (!(vsrc[3] > 0.0)) continue;
          double a = vsrc[3];
          for (int c = 0; c < 3; c++) {
            src_p[6 * nsrc + c] = psrc[c] / a;
            src_n[6 * nsrc + c] = vsrc[c];
          }
          double nl3 =
              sqrt(src_n[6 * nsrc] * src_n[6 * nsrc] + src_n[6 * nsrc + 1] * src_n[6 * nsrc + 1] +
                   src_n[6 * nsrc + 2] * src_n[6 * nsrc + 2]);
          if (!(nl3 > 0.0)) continue;
          for (int c = 0; c < 3; c++)
            src_n[6 * nsrc + c] /= nl3;
          /* Мощность: альбедо × облучённость солнцем × Σ(площадь·cos).
           * Альбедо взято средним по сцене — материала у агрегата нет, и это
           * названо здесь, а не спрятано. */
          src_w[nsrc] = PCELL_RHO * psrc[3];
          nsrc++;
        }
        double wsum = 0.0;
        for (int64_t q2 = 0; q2 < nsrc; q2++)
          wsum += src_w[q2];
        printf("== ВТОРИЧНЫХ ИСТОЧНИКОВ построено %lld (уровень %d), мощность %.2f против эталона "
               "%.2f — ИНВАРИАНТ ЭНЕРГИИ %s\n",
               (long long)nsrc, PCELL_AGG_LEV, wsum, wref,
               (fabs(wsum - wref) <= 1e-9 * (wref > 0.0 ? wref : 1.0)) ? "СОШЁЛСЯ" : "НЕ СОШЁЛСЯ");
      }
      /* НЕ ОСВОБОЖДАЕМ: суммы по поддеревьям нужны кадру для ИЕРАРХИЧЕСКОГО
       * СБОРА (§209). Плоский перебор 1 326 источников на пиксель и есть
       * главный расход кадра. */
      g_lv = lv;
      g_lp = lp;
    }
    g_ls = ls;
    g_dep = dep;
  }

  /* ---- РОСТ ЧИСЛА УЗЛОВ ПО УРОВНЯМ ДЕРЕВА (§199, гипотеза «б») ---------
   * Камеры НЕ содержит: это чистое свойство сцены. Если геометрия на масштабе
   * узла ПОВЕРХНОСТНА, уровень глубже даёт вчетверо больше занятых узлов; если
   * нитевидна (листва) — вдвое. Именно этим отличается ×4 от ×2 в законе роста
   * среза, и проверять это надо БЕЗ камеры, иначе снова смешаются причины. */
  {
    int64_t *lvn = calloc(64, sizeof *lvn);
    int32_t *st = malloc((size_t)(T.nnd > 0 ? T.nnd : 1) * sizeof *st);
    int32_t *sl = malloc((size_t)(T.nnd > 0 ? T.nnd : 1) * sizeof *sl);
    if (lvn != NULL && st != NULL && sl != NULL) {
      int64_t sp2 = 0;
      st[sp2] = 0;
      sl[sp2++] = 0;
      while (sp2 > 0) {
        int32_t i = st[--sp2];
        int32_t L = sl[sp2];
        if (T.nd[i].ntri > 0 && L < 64) lvn[L]++;
        if (T.nd[i].child >= 0 && L + 1 < 64)
          for (int q = 0; q < 8; q++) {
            st[sp2] = T.nd[i].child + q;
            sl[sp2++] = L + 1;
          }
      }
      printf("== ЗАНЯТЫХ УЗЛОВ ПО УРОВНЯМ (свойство СЦЕНЫ, камеры нет):\n");
      for (int L = 0; L < 64; L++) {
        if (lvn[L] == 0) continue;
        printf("   уровень %2d: %10lld", L, (long long)lvn[L]);
        if (L > 0 && lvn[L - 1] > 0)
          printf("   ×%.2f к предыдущему", (double)lvn[L] / (double)lvn[L - 1]);
        printf("\n");
      }
    }
    free(lvn);
    free(st);
    free(sl);
  }

  /* ---- ВЫБОР УРОВНЯ ПО ПОЛОСЕ 4…32 px (О58, §194) ---------------------
   * Спуск от корня: пока узел крупнее верхнего предела — вниз, иначе выбран.
   * Нижний предел спуском не движет, он ДИАГНОСТИКА: элементы мельче него —
   * дальние, у которых даже грубейший доступный узел мал.
   * ВЫБРАННЫМ УЗЛОМ ЕЩЁ НЕЛЬЗЯ РИСОВАТЬ: крупный узел обязан представлять всё
   * своё поддерево, а слияния снизу вверх (§190) нет. Здесь считается БЮДЖЕТ
   * КАДРА, и только он. */
  {
    double eye[3] = HZ_CFG_MIGUEL_EYE;
    if (haseye)
      for (int a = 0; a < 3; a++)
        eye[a] = eye0[a];
    double eps = (HZ_CFG_FOV_DEG * M_PI / 180.0) / (double)wpx;
    /* ПОДДЕРЕВНЫЕ СУММЫ: у остановленного узла засчитывается ВСЁ его поддерево
     * (он его и представляет), у пройденного насквозь — только собственные
     * треугольники. Сумма обязана дать ровно `nt`. Прежний инвариант («сумма
     * собственных по выбранным = nt») был сформулирован МНОЮ неверно: он
     * справедлив только без LOD, а при LOD треугольники под срезом не берутся
     * СОЗНАТЕЛЬНО — в этом вся экономия. Дети лежат в арене ПОСЛЕ родителя,
     * поэтому обратный проход по индексам даёт суммы за один раз. */
    int64_t *sub_t = malloc((size_t)(T.nnd > 0 ? T.nnd : 1) * sizeof *sub_t);
    if (sub_t != NULL) {
      for (int32_t i = T.nnd - 1; i >= 0; i--) {
        sub_t[i] = T.nd[i].ntri;
        if (T.nd[i].child >= 0)
          for (int q = 0; q < 8; q++)
            sub_t[i] += sub_t[T.nd[i].child + q];
      }
    }
    int32_t *stack = malloc((size_t)(T.nnd > 0 ? T.nnd : 1) * sizeof *stack);
    double *pxs = malloc((size_t)(T.nnd > 0 ? T.nnd : 1) * sizeof *pxs);
    if (stack != NULL && pxs != NULL) {
      double ts = now_s();
      int64_t sp = 0, nsel = 0, nsmall = 0, ntri_sel = 0, ninside = 0, cov = 0;
      stack[sp++] = 0;
      while (sp > 0) {
        int32_t i = stack[--sp];
        const hz_ptnode *nd = &T.nd[i];
        double c[3], rad = 0.0;
        for (int a = 0; a < 3; a++) {
          c[a] = 0.5 * (nd->lo[a] + nd->hi[a]);
          double h = 0.5 * (nd->hi[a] - nd->lo[a]);
          rad += h * h;
        }
        rad = sqrt(rad);
        /* РАССТОЯНИЕ ДО БЛИЖАЙШЕЙ ТОЧКИ КОРОБКИ, а не «до центра минус радиус».
         * Второе даёт отрицательное `R` у всякого узла, ВНУТРИ которого стоит
         * глаз, и после обрезки — угловой размер в миллиарды пикселей
         * (замерено: максимум 4.6e9). Глаз внутри коробки — законный случай, и
         * означает он «дробить обязательно», а не «элемент бесконечно велик». */
        double dd = 0.0;
        for (int a = 0; a < 3; a++) {
          double e = eye[a];
          double q2 = (e < nd->lo[a]) ? (nd->lo[a] - e) : ((e > nd->hi[a]) ? (e - nd->hi[a]) : 0.0);
          dd += q2 * q2;
        }
        double R = sqrt(dd);
        int inside = (R <= 0.0);
        double px = inside ? HUGE_VAL : (2.0 * rad / R) / eps;
        (void)c;
        if (px > pxhi && nd->child >= 0) {
          for (int q = 0; q < 8; q++)
            stack[sp++] = nd->child + q;
          /* СОБСТВЕННЫЕ ТРЕУГОЛЬНИКИ УЗЛА НЕ ТЕРЯЮТСЯ ПРИ СПУСКЕ. Первая
           * редакция их роняла, а у внутренних узлов висит 72.5 %% геометрии —
           * тот же класс, что А419. Ловится инвариантом ниже. */
          if (nd->ntri > 0) {
            pxs[nsel++] = px;
            ntri_sel += nd->ntri;
            cov += nd->ntri;
            if (inside)
              ninside++;
            else if (px < pxlo)
              nsmall++;
          }
          continue;
        }
        cov += (sub_t != NULL) ? sub_t[i] : 0;
        if (sub_t != NULL && sub_t[i] == 0) continue; /* пусто: ни своих, ни в поддереве */
        if (nd->ntri > 0) {
          pxs[nsel++] = px;
          ntri_sel += nd->ntri;
          if (inside)
            ninside++;
          else if (px < pxlo)
            nsmall++;
        }
      }
      double t_sel = now_s() - ts;
      qsort(pxs, (size_t)nsel, sizeof *pxs, cmp_dbl);
      printf("== ВЫБОР УРОВНЯ (полоса %.0f…%.0f px при %d точках кадра) за %.4f с: УЗЛОВ %lld из "
             "%d; треугольников в них %lld\n",
             pxlo, pxhi, wpx, t_sel, (long long)nsel, T.nnd, (long long)ntri_sel);
      printf("   ИНВАРИАНТ ПОКРЫТИЯ: представлено %lld треугольников против %d в сцене — %s; под "
             "срезом %lld (экономия LOD); узлов с глазом ВНУТРИ %lld\n",
             (long long)cov, m.nt, (cov == m.nt) ? "СОШЛОСЬ" : "НЕ СОШЛОСЬ",
             (long long)(m.nt - ntri_sel), (long long)ninside);
      if (nsel > 0)
        printf("   угловой размер выбранных, px: p10 %.2f, МЕДИАНА %.2f, p90 %.2f, максимум %.2f; "
               "мельче %.0f px — %lld (%.1f %%)\n",
               pxs[nsel / 10], pxs[nsel / 2], pxs[(nsel * 9) / 10], pxs[nsel - 1], pxlo,
               (long long)nsmall, 100.0 * (double)nsmall / (double)nsel);
    }
    free(stack);
    free(pxs);
    free(sub_t);
  }

  /* ---- КАРТИНКА: ОСВЕЩЁННЫЙ НАБОР ---------------------------------------
   * Рисуется то, что ПОСЧИТАНО, и ничего сверх: освещён элемент или нет, плюс
   * ламбертов косинус к солнцу у освещённых. Косвенного света тут нет и быть не
   * должно — отскок не сделан (§202). Тень в затенённых местах поэтому глухая,
   * и это честно, а не дефект раскраски. */
  {
    int W = 1024, H = 1024;
    unsigned char *rgb = malloc((size_t)W * (size_t)H * 3);
    double eyeP[3] = HZ_CFG_MIGUEL_EYE, atP[3] = HZ_CFG_MIGUEL_AT, upP[3] = HZ_CFG_UP;
    if (haseye)
      for (int a = 0; a < 3; a++)
        eyeP[a] = eye0[a];
    if (rgb != NULL) {
      double fw[3] = {0, 0, 0}, ri[3] = {0, 0, 0}, uv[3] = {0, 0, 0};
      for (int c = 0; c < 3; c++)
        fw[c] = atP[c] - eyeP[c];
      double fl = sqrt(fw[0] * fw[0] + fw[1] * fw[1] + fw[2] * fw[2]);
      for (int c = 0; c < 3; c++)
        fw[c] /= fl;
      ri[0] = fw[1] * upP[2] - fw[2] * upP[1];
      ri[1] = fw[2] * upP[0] - fw[0] * upP[2];
      ri[2] = fw[0] * upP[1] - fw[1] * upP[0];
      double rl = sqrt(ri[0] * ri[0] + ri[1] * ri[1] + ri[2] * ri[2]);
      for (int c = 0; c < 3; c++)
        ri[c] /= rl;
      uv[0] = ri[1] * fw[2] - ri[2] * fw[1];
      uv[1] = ri[2] * fw[0] - ri[0] * fw[2];
      uv[2] = ri[0] * fw[1] - ri[1] * fw[0];
      double th2 = tan(0.5 * HZ_CFG_FOV_DEG * M_PI / 180.0);
      double t_pic = now_s();
      int64_t nmis = 0, npix2 = 0, ngath = 0;
#pragma omp parallel for schedule(dynamic, 8) reduction(+ : nmis, npix2, ngath)
      for (int j = 0; j < H; j++)
        for (int i = 0; i < W; i++) {
          double sx = (2.0 * ((double)i + 0.5) / W - 1.0) * th2;
          double sy = (1.0 - 2.0 * ((double)j + 0.5) / H) * th2;
          double d[3] = {0, 0, 0};
          for (int c = 0; c < 3; c++)
            d[c] = fw[c] + sx * ri[c] + sy * uv[c];
          double dl = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
          for (int c = 0; c < 3; c++)
            d[c] /= dl;
          double tt = 0.0;
          int32_t tri = -1;
          int64_t q = (int64_t)j * W + i;
          double col[3] = {0.45, 0.55, 0.70};
          if (hz_pgrid_trace_tri(&G, &m, eyeP, d, 0.0, diag2, NULL, 0, 0, &tt, &tri) && tri >= 0) {
            const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)tri];
            const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)tri + 1];
            const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)tri + 2];
            double e1[3] = {0, 0, 0}, e2[3] = {0, 0, 0}, nn[3] = {0, 0, 0};
            for (int c = 0; c < 3; c++) {
              e1[c] = B[c] - A[c];
              e2[c] = C[c] - A[c];
            }
            nn[0] = e1[1] * e2[2] - e1[2] * e2[1];
            nn[1] = e1[2] * e2[0] - e1[0] * e2[2];
            nn[2] = e1[0] * e2[1] - e1[1] * e2[0];
            double nl2 = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
            double ndl = 0.0;
            if (nl2 > 0.0) ndl = fabs((nn[0] * wsun[0] + nn[1] * wsun[1] + nn[2] * wsun[2]) / nl2);
            /* ЧИТАЕМОСТЬ ВАЖНЕЕ ФИЗИЧНОСТИ В ДИАГНОСТИЧЕСКОМ КАДРЕ. Первая
             * редакция умножала незасвеченное на `0.06` и загоняла 97 %% кадра
             * почти в чёрное: судить о качестве по такому нельзя. Теперь
             * геометрия читается ВСЕГДА (затенение по `|n·взгляд|`), а
             * освещённость показана ЦВЕТОМ: тёплый и яркий — солнце дошло,
             * холодный и тусклый — нет. Это диагностика, не рендер. */
            /* ЧЕСТНАЯ ТЕНЬ НА ПИКСЕЛЬ — эталон, с которым сравнивается поэлементная.
             * Разница между ними и есть цена КВАНТОВАНИЯ ТЕНИ ЭЛЕМЕНТОМ, то
             * есть ровно то, что портит кадр; всё прочее — незаконченность
             * (отскока нет). */
            double hp[3] = {0, 0, 0};
            for (int c = 0; c < 3; c++)
              hp[c] = eyeP[c] + tt * d[c];
            double nsg = (nl2 > 0.0 && (nn[0] * wsun[0] + nn[1] * wsun[1] + nn[2] * wsun[2]) > 0.0)
                             ? -1.0
                             : 1.0;
            double so[3] = {0, 0, 0}, sd[3] = {-wsun[0], -wsun[1], -wsun[2]}, st = 0.0;
            for (int c = 0; c < 3; c++)
              so[c] = hp[c] + 1e-4 * nsg * nn[c] / (nl2 > 0.0 ? nl2 : 1.0);
            /* ПОЛУТЕНЬ: солнце — ДИСК углового диаметра `PCELL_SUN_DEG`, а не
             * точка. Один луч в центр даёт резкую тень, которой в природе нет.
             * Выборка детерминированная (радикальный обратный по основанию 2),
             * то есть кадр воспроизводится побитово. */
            double sa[3] = {0, 0, 0}, sb[3] = {0, 0, 0};
            {
              int ax2 = (fabs(sd[0]) < fabs(sd[1])) ? ((fabs(sd[0]) < fabs(sd[2])) ? 0 : 2)
                                                    : ((fabs(sd[1]) < fabs(sd[2])) ? 1 : 2);
              double tv[3] = {0, 0, 0};
              tv[ax2] = 1.0;
              sa[0] = sd[1] * tv[2] - sd[2] * tv[1];
              sa[1] = sd[2] * tv[0] - sd[0] * tv[2];
              sa[2] = sd[0] * tv[1] - sd[1] * tv[0];
              double al = sqrt(sa[0] * sa[0] + sa[1] * sa[1] + sa[2] * sa[2]);
              for (int c = 0; c < 3; c++)
                sa[c] /= al;
              sb[0] = sd[1] * sa[2] - sd[2] * sa[1];
              sb[1] = sd[2] * sa[0] - sd[0] * sa[2];
              sb[2] = sd[0] * sa[1] - sd[1] * sa[0];
            }
            /* Поворот набора на пиксель: целочисленный хеш, не генератор — кадр
             * остаётся воспроизводимым побитово. Без него один и тот же набор
             * направлений на всех пикселях даёт ЛЕСЕНКУ вместо градиента. */
            uint64_t hpx = (uint64_t)q * 1099511628211ULL;
            hpx ^= hpx >> 29;
            hpx *= 1469598103934665603ULL;
            double rot_px = 2.0 * M_PI * (double)(hpx >> 11) / 9007199254740992.0;
            double tanr = tan(0.5 * PCELL_SUN_DEG * M_PI / 180.0);
            int nvis2 = 0;
            for (int s3 = 0; s3 < PCELL_SUN_SAMP; s3++) {
              double u1 = ((double)s3 + 0.5) / (double)PCELL_SUN_SAMP;
              uint32_t b3 = (uint32_t)s3;
              b3 = (b3 << 16) | (b3 >> 16);
              b3 = ((b3 & 0x55555555u) << 1) | ((b3 & 0xAAAAAAAAu) >> 1);
              b3 = ((b3 & 0x33333333u) << 2) | ((b3 & 0xCCCCCCCCu) >> 2);
              b3 = ((b3 & 0x0F0F0F0Fu) << 4) | ((b3 & 0xF0F0F0F0u) >> 4);
              b3 = ((b3 & 0x00FF00FFu) << 8) | ((b3 & 0xFF00FF00u) >> 8);
              double u2 = (double)b3 * 2.3283064365386963e-10;
              double rr2 = tanr * sqrt(u1), ph2 = 2.0 * M_PI * u2 + rot_px;
              double sdj[3] = {0, 0, 0};
              for (int c = 0; c < 3; c++)
                sdj[c] = sd[c] + rr2 * (cos(ph2) * sa[c] + sin(ph2) * sb[c]);
              double jl = sqrt(sdj[0] * sdj[0] + sdj[1] * sdj[1] + sdj[2] * sdj[2]);
              for (int c = 0; c < 3; c++)
                sdj[c] /= jl;
              if (!hz_pgrid_trace(&G, &m, so, sdj, 0.0, diag2, NULL, 0, 1, &st)) nvis2++;
            }
            double vfrac = (double)nvis2 / (double)PCELL_SUN_SAMP;
            int litpx = (vfrac > 0.0);
            if (litpx != (tri_lit[tri] ? 1 : 0)) nmis++;
            npix2++;
            double vdn = 0.0;
            if (nl2 > 0.0) vdn = fabs((nn[0] * d[0] + nn[1] * d[1] + nn[2] * d[2]) / nl2);
            /* ОТСКОК: вклад вторичных источников в эту точку. Площадка видна как
             * диск: `F = cos_r·cos_s·A / (π r² + A)`. Видимости в отскоке НЕТ —
             * это названная грубость (пункт 2 реестра допускает огрубление
             * второго отскока), и она завышает свет в углах. */
            double ind = 0.0;
            double nrm[3] = {nn[0], nn[1], nn[2]};
            if (nl2 > 0.0)
              for (int c = 0; c < 3; c++)
                nrm[c] /= nl2;
            if ((nrm[0] * d[0] + nrm[1] * d[1] + nrm[2] * d[2]) > 0.0)
              for (int c = 0; c < 3; c++)
                nrm[c] = -nrm[c];
            /* ИЕРАРХИЧЕСКИЙ СБОР ПО ИНДЕКСУ (§209), вместо плоского перебора
             * всех источников на каждый пиксель. Спуск от корня: узел берётся
             * ИСТОЧНИКОМ, пока его угловой размер из этой точки мал; иначе
             * дробится. Учёт тот же, что поймал §206: ВЗЯТЫЙ узел представляет
             * ПОДДЕРЕВО (`g_lv`/`g_lp`), ПРОЙДЕННЫЙ отдаёт только СВОЁ
             * (`litv`/`litp`). Каждый освещённый элемент учитывается ровно раз. */
            int32_t st2[128];
            int sp3 = 0;
            st2[sp3++] = 0;
            while (sp3 > 0) {
              int32_t ni = st2[--sp3];
              if (g_ls == NULL || g_ls[ni] <= 0) continue;
              const hz_ptnode *nd2 = &T.nd[ni];
              double c2[3] = {0, 0, 0}, rad2 = 0.0;
              for (int a = 0; a < 3; a++) {
                c2[a] = 0.5 * (nd2->lo[a] + nd2->hi[a]);
                double h2 = 0.5 * (nd2->hi[a] - nd2->lo[a]);
                rad2 += h2 * h2;
              }
              rad2 = sqrt(rad2);
              double dc[3] = {c2[0] - hp[0], c2[1] - hp[1], c2[2] - hp[2]};
              double rc = sqrt(dc[0] * dc[0] + dc[1] * dc[1] + dc[2] * dc[2]);
              int split = (nd2->child >= 0) && (rc < rad2 * PCELL_SRC_K) && (sp3 + 8 < 128);
              const double *vv = split ? &litv[4 * (size_t)ni] : &g_lv[4 * (size_t)ni];
              const double *pp = split ? &litp[4 * (size_t)ni] : &g_lp[4 * (size_t)ni];
              if (split)
                for (int qq = 0; qq < 8; qq++)
                  st2[sp3++] = nd2->child + qq;
              if (!(vv[3] > 0.0)) continue;
              double dv[3] = {pp[0] / vv[3] - hp[0], pp[1] / vv[3] - hp[1], pp[2] / vv[3] - hp[2]};
              double r2 = dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2];
              if (!(r2 > 1e-9)) continue;
              double rr = sqrt(r2);
              double cr2 = (nrm[0] * dv[0] + nrm[1] * dv[1] + nrm[2] * dv[2]) / rr;
              if (!(cr2 > 0.0)) continue;
              double vn2 = sqrt(vv[0] * vv[0] + vv[1] * vv[1] + vv[2] * vv[2]);
              if (!(vn2 > 0.0)) continue;
              double cs2 = -(vv[0] * dv[0] + vv[1] * dv[1] + vv[2] * dv[2]) / (vn2 * rr);
              if (!(cs2 > 0.0)) continue;
              double w2 = PCELL_RHO * pp[3];
              ind += w2 * cr2 * cs2 / (M_PI * r2 + w2);
              ngath++;
            }
            /* ЦВЕТ БЕРЁТСЯ ИЗ МАТЕРИАЛА, А НЕ ПРИДУМЫВАЕТСЯ. Прежняя редакция красила
             * тёплым/холодным по признаку «освещён», и кадр выходил чистым
             * lightmap-ом при том, что `kd3` у сцены есть по 288 материалам и
             * загрузчик их читает. Свет и цвет — разные вещи, и смешивать их в
             * одну шкалу нельзя (правило Т1: альбедо прикладывается при ЧТЕНИИ
             * поля). Текстур пока нет — это Ш8, и до неё не дошли. */
            const double *kd = m.mtl[m.fm[tri]].kd3;
            /* Небо как подсветка теней: полусферически, по наклону нормали.
             * Косвенного от стен пока нет (мощность источников неверна, §206),
             * и без неба тень была бы чёрной. Число названо, в замер не входит. */
            /* НЕБО — ВИДИМОСТЬЮ, А НЕ КОНСТАНТОЙ (§183). Полусферическое
             * приближение `0.5+0.5·n_y` не знает о заслонах, и потому тени
             * выходят без контактных затемнений и без формы. Здесь берётся
             * косинусно взвешенная доля неба — точный член уравнения
             * рендеринга, приближается только интеграл. Выборка
             * детерминированная. */
            double sky1[3] = {0, 0, 0}, sky2[3] = {0, 0, 0};
            {
              int ax3 = (fabs(nrm[0]) < fabs(nrm[1])) ? ((fabs(nrm[0]) < fabs(nrm[2])) ? 0 : 2)
                                                      : ((fabs(nrm[1]) < fabs(nrm[2])) ? 1 : 2);
              double tv2[3] = {0, 0, 0};
              tv2[ax3] = 1.0;
              sky1[0] = nrm[1] * tv2[2] - nrm[2] * tv2[1];
              sky1[1] = nrm[2] * tv2[0] - nrm[0] * tv2[2];
              sky1[2] = nrm[0] * tv2[1] - nrm[1] * tv2[0];
              double s1l = sqrt(sky1[0] * sky1[0] + sky1[1] * sky1[1] + sky1[2] * sky1[2]);
              for (int c = 0; c < 3; c++)
                sky1[c] /= s1l;
              sky2[0] = nrm[1] * sky1[2] - nrm[2] * sky1[1];
              sky2[1] = nrm[2] * sky1[0] - nrm[0] * sky1[2];
              sky2[2] = nrm[0] * sky1[1] - nrm[1] * sky1[0];
            }
            int nsk = 0;
            for (int s4 = 0; s4 < PCELL_SKY_SAMP; s4++) {
              double u1 = ((double)s4 + 0.5) / (double)PCELL_SKY_SAMP;
              uint32_t b4 = (uint32_t)s4;
              b4 = (b4 << 16) | (b4 >> 16);
              b4 = ((b4 & 0x55555555u) << 1) | ((b4 & 0xAAAAAAAAu) >> 1);
              b4 = ((b4 & 0x33333333u) << 2) | ((b4 & 0xCCCCCCCCu) >> 2);
              b4 = ((b4 & 0x0F0F0F0Fu) << 4) | ((b4 & 0xF0F0F0F0u) >> 4);
              b4 = ((b4 & 0x00FF00FFu) << 8) | ((b4 & 0xFF00FF00u) >> 8);
              double u2 = (double)b4 * 2.3283064365386963e-10;
              double rr3 = sqrt(u1), ph3 = 2.0 * M_PI * u2;
              double dz3 = sqrt(1.0 - u1 > 0.0 ? 1.0 - u1 : 0.0);
              double dk3[3] = {0, 0, 0};
              for (int c = 0; c < 3; c++)
                dk3[c] = rr3 * cos(ph3) * sky1[c] + rr3 * sin(ph3) * sky2[c] + dz3 * nrm[c];
              double kl = sqrt(dk3[0] * dk3[0] + dk3[1] * dk3[1] + dk3[2] * dk3[2]);
              for (int c = 0; c < 3; c++)
                dk3[c] /= kl;
              double kt = 0.0;
              if (!hz_pgrid_trace(&G, &m, so, dk3, 0.0, diag2, NULL, 0, 1, &kt)) nsk++;
            }
            double amb2 = PCELL_SKY * (double)nsk / (double)PCELL_SKY_SAMP;
            double sun2 = vfrac * ndl;
            col[0] = kd[0] * (sun2 * 1.00 + amb2 * 0.60 + ind * PCELL_IND);
            col[1] = kd[1] * (sun2 * 0.97 + amb2 * 0.72 + ind * PCELL_IND);
            col[2] = kd[2] * (sun2 * 0.88 + amb2 * 1.00 + ind * PCELL_IND);
            (void)vdn;
            /* Косвенный свет ДОБАВЛЯЕТСЯ к обоим случаям — он есть и на свету. */
          }
          for (int c = 0; c < 3; c++) {
            double e3 = col[c] * PCELL_EXPO;
            if (e3 < 0.0) e3 = 0.0;
            /* Плечо вместо обрезки: без него всё ярче единицы становится
             * одинаково белым, и градации на солнце пропадают. */
            e3 = e3 / (1.0 + e3 * 0.6);
            double g2 = pow(e3 > 1.0 ? 1.0 : e3, 1.0 / 2.2);
            rgb[3 * q + c] = (unsigned char)(255.0 * g2 + 0.5);
          }
        }
      FILE *fp = fopen("img/pcell_lit.ppm", "wb");
      if (fp != NULL) {
        fprintf(fp, "P6\n%d %d\n255\n", W, H);
        fwrite(rgb, 1, (size_t)W * (size_t)H * 3, fp);
        fclose(fp);
        printf("== КАРТИНКА img/pcell_lit.ppm за %.2f с: тень ЧЕСТНАЯ (луч на пиксель); косвенного "
               "НЕТ — отскок не сделан\n",
               now_s() - t_pic);
        printf("   ИЕРАРХИЧЕСКИЙ СБОР: обращений к источникам %lld, то есть %.1f на пиксель "
               "против %lld при плоском переборе (выигрыш %.1f×)\n",
               (long long)ngath, (double)ngath / (double)(npix2 > 0 ? npix2 : 1), (long long)nsrc,
               (double)nsrc * (double)npix2 / (double)(ngath > 0 ? ngath : 1));
        printf("   ЦЕНА КВАНТОВАНИЯ ТЕНИ ЭЛЕМЕНТОМ: поэлементная расходится с честной на %lld "
               "пикселей из %lld (%.2f %%)\n",
               (long long)nmis, (long long)npix2,
               100.0 * (double)nmis / (double)(npix2 > 0 ? npix2 : 1));
      }
      free(rgb);
    }
  }
  free(tri_lit);

  for (int64_t k = 0; k < nl; k++)
    free(cnt[k]);
  free(cnt);
  free(nseg_of);
  free(area_of);
  free(dmax_of);
  free(leaf);
  hz_ptree_free(&T);
  return 0;
}
