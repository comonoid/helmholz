/* render3.c — ПОЛНОЦЕННЫЙ РЕНДЕР линии переноса: сфера в освещённой комнате,
 * глобальное освещение развёрткой по ординатам.
 *
 * PLAN_TRANSPORT.md, все три прохода вместе:
 *   развёртка   — полное решение уравнения переноса с многократными отражениями;
 *   поверхности — граничное условие, исходящий радианс хранится DG1 по положению;
 *   сбор        — луч на пиксель до первого попадания, чтение хранимого В ТОЧКЕ.
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
#include "transport/ray3.h"
#include "transport/sweep3.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

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

static int rec_cmp(const void *a, const void *b) {
  const rec_t *x = a, *y = b;
  return x->cell < y->cell ? -1 : (x->cell > y->cell ? 1 : 0);
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
  if (log2n < 1 || log2n > 8) {
    fprintf(stderr, "log2n вне [1,8]\n");
    return 1;
  }
  const int nc = 1 << log2n;

  stage_mark();
  hz_frame fr = {{0, 0, 0}, {ROOM / nc, ROOM / nc, ROOM / nc}};
  hz_octree t;
  if (hz_oct_init(&t, log2n, 0.0)) return 1;
  /* РАВНОМЕРНОЕ дробление: условие 1:1 у разрезанных ячеек (cut3) требует, чтобы
   * грань сетки совпадала с гранью коробки. Градуированную сетку у поверхности
   * пришлось бы ещё и обрезать прямоугольником — это отдельная работа. */
  for (int x = 0; x < nc; x++)
    for (int y = 0; y < nc; y++)
      for (int z = 0; z < nc; z++) {
        int lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        hz_oct_set_box(&t, lo, hi, 1.0);
      }

  stage_add("октодерево", 0);

  /* сфера как АНАЛИТИЧЕСКИЙ примитив плюс её фасеты для разреза */
  hz_surftab stab;
  hz_facettab ftab;
  hz_cutmap cmap;
  if (hz_surftab_init(&stab) || hz_facettab_init(&ftab) || hz_cutmap_init(&cmap)) return 1;
  /* ДВА ТЕЛА: тени одного на другом и переотражения между ними — то, ради чего
   * развёртка и нужна. Материал в модели разреза есть ПЕРЕСЕЧЕНИЕ полуплоскостей,
   * то есть ОДНО выпуклое тело на ячейку; два тела в одной ячейке дали бы
   * пересечение вместо объединения. Тела разнесены, и это ПРОВЕРЯЕТСЯ. */
  const int NB = 2;
  double sc[2][3] = {{5.6, 9.2, 4.2}, {10.9, 7.4, 3.1}};
  double sr[2] = {3.0, 1.9};
  int32_t f0[2], nfac[2];
  for (int b = 0; b < NB; b++) {
    hz_surf sp = {HZ_SURF_SPHERE, {sc[b][0], sc[b][1], sc[b][2], sr[b], 0, 0, 0}, 1, 0};
    int32_t si = hz_surftab_add(&stab, &sp);
    f0[b] = 0;
    nfac[b] = hz_surf_facet_sphere(&ftab, &fr, sc[b], sr[b], 2, HZ_FIT_MEAN_SAGITTA, si, &f0[b]);
    if (nfac[b] <= 0) return 1;
  }

  stage_add("фасетизация примитивов", 0);

  /* боковая таблица: ключи ОБЯЗАНЫ идти по возрастанию (Г45) */
  rec_t *recs = calloc((size_t)nc * (size_t)nc * (size_t)nc, sizeof(rec_t));
  if (recs == NULL) return 1;
  int nrec = 0, nboth = 0;
  for (int x = 0; x < nc; x++)
    for (int y = 0; y < nc; y++)
      for (int z = 0; z < nc; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        int used = -1, ns = 0;
        for (int b = 0; b < NB; b++) {
          int32_t s2[HZ_P3_MAXH];
          int k2 = hz_facets_for_box(&ftab, f0[b], nfac[b], lo, hi, s2, HZ_P3_MAXH);
          if (k2 <= 0) continue;
          if (used >= 0) {
            nboth++;
            continue;
          }
          used = b;
          ns = k2;
          for (int j = 0; j < k2; j++)
            sel[j] = s2[j];
        }
        if (ns <= 0) continue;
        int32_t ni = hz_oct_leaf(&t, x, y, z);
        recs[nrec].cell = ni;
        recs[nrec].nf = ns;
        for (int j = 0; j < ns; j++)
          recs[nrec].f[j] = sel[j];
        nrec++;
      }
  /* СОРТИРОВКА `qsort`, А НЕ ВСТАВКАМИ, И ЭТО НЕ УКРАШЕНИЕ.
   * Вставками было `O(n^2)` при записи в 392 байта. На сетке 16^3 записей около
   * тысячи, и это не мешало; на 64^3 их десятки тысяч, и построитель сцены
   * вставал на десятки минут ПЕРЕД развёрткой — то есть свип по размеру ячейки
   * (К68) упирался в оснастку, а не в схему. Ключи различны (одна запись на
   * ячейку), поэтому порядок совпадает с прежним, а не просто похож. */
  qsort(recs, (size_t)nrec, sizeof(rec_t), rec_cmp);
  for (int i = 0; i < nrec; i++)
    if (hz_cutmap_add(&cmap, recs[i].cell, recs[i].f, recs[i].nf) != 0) return 1;
  free(recs);

  stage_add("боковая таблица (отбор фасетов по ячейкам)", 0);

  tr3_mesh mesh;
  if (tr3_mesh_build(&mesh, &t, &fr)) return 1;
  stage_add("сетка ГРАНЕЙ над деревом", 0);
  /* Г38: маска полных ячеек строится ОТДЕЛЬНО — в боковой таблице их нет */
  uint8_t *solid = calloc((size_t)mesh.ncell, 1);
  if (solid == NULL) return 1;
  for (int x = 0; x < nc; x++)
    for (int y = 0; y < nc; y++)
      for (int z = 0; z < nc; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        for (int b = 0; b < NB; b++)
          if (hz_facets_for_box(&ftab, f0[b], nfac[b], lo, hi, sel, HZ_P3_MAXH) == 0)
            solid[mesh.cellof[hz_oct_leaf(&t, x, y, z)]] = 1;
      }
  stage_add("маска полных ячеек (Г38)", 0);
  tr3_cut cut;
  if (tr3_cut_build(&cut, &mesh, &ftab, &cmap, solid)) return 1;
  stage_add("РАЗРЕЗ: флюид, матрицы масс, поверхностные элементы", 0);
  tr3_dirs dirs;
  if (tr3_dirs_product(&dirs, nmu, nmu)) return 1;
  stage_add("набор ординат", 0);

  double *sig_t = calloc((size_t)mesh.ncell, sizeof(double));
  double *sig_s = calloc((size_t)mesh.ncell, sizeof(double));
  if (sig_t == NULL || sig_s == NULL) return 1;
  for (int32_t c = 0; c < mesh.ncell; c++) {
    sig_t[c] = 0.015; /* лёгкая дымка: среда участвует, но не глушит */
    sig_s[c] = 0.012; /* альбедо 0.8, а не 1: при 1 ряд по рассеянию сходится как
                       * ρⁿ и упирается ровно в то, ради чего С1 и заводил DSA */
    ;
  }
  /* комната: потолок светит, пол и стены серые, одна стена красноватая по
   * яркости (цвета нет — считается один спектральный канал) */
  double wr[6] = {0.72, 0.35, 0.72, 0.72, 0.65, 0.05};
  for (int i = 0; i < 6; i++) {
    wr[i] *= albs;
    if (wr[i] > 0.999) wr[i] = 0.999;
  }
  double we[6] = {0, 0, 0, 0, 0, 6.0};
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

  tr3_problem prob = {.m = &mesh,
                      .d = &dirs,
                      .cut = &cut,
                      .facet_rho = frho,
                      .facet_emit = femit,
                      .nfacet = ftab.n,
                      .sig_t = sig_t,
                      .sig_s = sig_s,
                      .wall_rho = wr,
                      .wall_emit = we,
                      .limiter = limiter,
                      .trace = trace};
  double *phi = calloc((size_t)mesh.ncell * 4, sizeof(double));
  if (phi == NULL) return 1;
  tr3_stats st;
  printf("ячеек %d, граней %d, поверхностных элементов %d, направлений %d, 1:1 нарушений %d, "
         "ячеек с ДВУМЯ телами %d (обязано быть 0)\n",
         mesh.ncell, mesh.nf, cut.nse, dirs.n, cut.nbad, nboth);
  stage_mark();
  double t0 = now();
  int rc = tr3_sweep_solve(&prob, maxit, stol, phi, &st);
  double t_sweep = now() - t0;
  stage_add("РАЗВЁРТКА (итерация по рассеянию)", 0);
  printf("развёртка (ХОЛОДНЫЙ СТАРТ): код %d, итераций %d, невязка %.2e, срезок %d, %.2f с "
         "(%.2f мс на итерацию)\n",
         rc, st.iters, st.resid, st.nclip, t_sweep, 1e3 * t_sweep / (double)(st.iters + 1));
  printf("баланс: втекло %.4f, вытекло %.4f, поглощено средой %.4f, ушло в поверхности %.4f, "
         "отдано ими %.4f\n",
         st.pin, st.pout, st.pabs, st.psin, st.psout);
  printf("        невязка %.3e, она же на втекшее %.3e\n", st.balance,
         st.pin > 0.0 ? fabs(st.balance) / st.pin : 0.0);
  /* К65: САМА ВЕЛИЧИНА, РАДИ КОТОРОЙ ОГРАНИЧИТЕЛЬ СТОИТ. Картинка её не видит:
   * она читает хранимое на ПОВЕРХНОСТЯХ, а ограничитель следит за полем в
   * ОБЪЁМЕ. Печатается минимум φ по углам ячеек — то самое, что уходит в минус
   * около границы тени, — и минимум хранимого. Без этих двух чисел вопрос
   * «нужен ли Крылову внешний нелинейный цикл» решается на глаз. */
  {
    double phimin = 1e300;
    int nneg = 0;
    for (int32_t c = 0; c < mesh.ncell; c++) {
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
           mesh.ncell);
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
    double *phi_n = calloc((size_t)mesh.ncell * 4, sizeof(double));
    double *phi_k = calloc((size_t)mesh.ncell * 4, sizeof(double));
    double *bk = calloc((size_t)mesh.nf * 4, sizeof(double));
    double *sk = calloc((size_t)(cut.nse > 0 ? cut.nse : 1) * 4, sizeof(double));
    if (phi_n == NULL || phi_k == NULL || bk == NULL || sk == NULL) return 1;

    tr3_problem pn = prob;
    pn.limiter = 0;
    tr3_stats stn;
    double t0n = now();
    int rcn = tr3_sweep_solve(&pn, maxit, stol, phi_n, &stn);
    double tn = now() - t0n;

    tr3_kstats stk;
    double t0k = now();
    int rck = tr3_krylov_solve(&prob, 200, stol, phi_k, bk, sk, &stk);
    double tk = now() - t0k;

    double dmax = 0.0, pmax = 0.0, dbm = 0.0, bmax = 0.0;
    for (int32_t i = 0; i < mesh.ncell * 4; i++) {
      double d = fabs(phi_n[i] - phi_k[i]);
      if (d > dmax) dmax = d;
      if (fabs(phi_n[i]) > pmax) pmax = fabs(phi_n[i]);
    }
    for (int32_t i = 0; i < mesh.nf * 4; i++) {
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
    double *phi_c = calloc((size_t)mesh.ncell * 4, sizeof(double));
    double *phi_w = calloc((size_t)mesh.ncell * 4, sizeof(double));
    if (phi_c == NULL || phi_w == NULL) return 1;
    memcpy(phi_w, phi, (size_t)mesh.ncell * 4 * sizeof(double));

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
    for (int32_t i = 0; i < mesh.ncell * 4; i++) {
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
    free(sc_st.bout);
    free(sc_st.sout);
    free(sw_st.bout);
    free(sw_st.sout);
    free(phi_c);
    free(phi_w);
  }

  /* --- сбор по пикселю --- */
  tr3_scene scn = {.tree = &t,
                   .fr = fr,
                   .sigma = NULL,
                   .st = &stab,
                   .ft = &ftab,
                   .cm = &cmap,
                   .facet_only = facet_only};
  tr3_camera cam;
  double eye[3] = {8.0, 0.6, 7.2}, at[3] = {8.2, 9.5, 3.6}, up[3] = {0, 0, 1};
  if (tr3_camera_look(&cam, eye, at, up, 1.3, W, H)) return 1;
  double *buf = calloc((size_t)W * (size_t)H, sizeof(double));
  if (buf == NULL) return 1;

  /* ИНДЕКС ГРАНИЧНЫХ ГРАНЕЙ ПО СЕТКЕ. Первая редакция искала грань ЛИНЕЙНЫМ
   * перебором всех 13056 граней на КАЖДЫЙ пиксель — это 2.7e10 сравнений на
   * кадр, и сбор стоил 3930 нс на луч при разумных 100-200. Грани лежат на
   * ЦЕЛОЧИСЛЕННОЙ сетке, поэтому индекс прямой: (стенка, u, v) -> грань. */
  int32_t *wallidx = calloc((size_t)6 * (size_t)nc * (size_t)nc, sizeof(int32_t));
  if (wallidx == NULL) return 1;
  for (int32_t i = 0; i < 6 * (int32_t)nc * (int32_t)nc; i++)
    wallidx[i] = -1;
  for (int32_t f = 0; f < mesh.nf; f++) {
    const tr3_face *ff = &mesh.f[f];
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
                   .m = &mesh,
                   .cut = &cut,
                   .bout = st.bout,
                   .sout = st.sout,
                   .facet_spec = spec > 0.0 ? fspec : NULL,
                   .nfacet = ftab.n,
                   .wallidx = wallidx,
                   .nwall = nc,
                   .maxbounce = maxbounce};
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
        buf[(size_t)py * (size_t)W + (size_t)px] = tr3_gather_ray(&gg, o, d, &nb);
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
  char pfm[512];
  size_t ln = strlen(out);
  if (ln + 5 < sizeof pfm) {
    memcpy(pfm, out, ln + 1);
    char *dot = strrchr(pfm, '.');
    size_t ext = dot != NULL ? (size_t)(dot - pfm) : ln;
    memcpy(pfm + ext, ".pfm", 5);
    /* R = G = B: спектральный канал ОДИН, и файл об этом не врёт */
    if (hz_pfm_write(pfm, buf, buf, buf, W, H) != 0) {
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
  free(solid);
  free(buf);
  free(phi);
  free(st.bout);
  free(st.sout);
  free(sig_t);
  free(sig_s);
  free(frho);
  free(fspec);
  free(femit);
  tr3_dirs_free(&dirs);
  tr3_cut_free(&cut);
  tr3_mesh_free(&mesh);
  hz_cutmap_free(&cmap);
  hz_facettab_free(&ftab);
  hz_surftab_free(&stab);
  hz_oct_free(&t);
  return 0;
}
