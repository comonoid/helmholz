/* psimp — УПРОЩЕНИЕ КРАЯ: СОГЛАСОВАННОЕ ПРОТИВ НЕЗАВИСИМОГО.
 *
 * PLAN_ELEMENTS.md, §10, открытый пункт 1. Ш7 (четвёртая поправка) намерил, что
 * край жмётся приближением `2.4×` даром, и там же назвал предел: каждый полигон
 * жмётся НЕЗАВИСИМО, общее ребро двух соседей перестаёт совпадать, и после 5 мм
 * появляются щели. Здесь мерится, что даёт согласование (`src/pedge.c`).
 *
 * ТРИ МЕРЫ, И КАЖДАЯ ОТВЕЧАЕТ ЗА СВОЁ.
 *   1. ВЕРШИН НА ПОЛИГОН — цена во внутреннем цикле трассировщика, ради которой
 *      всё и делается (`hz_poly_inside` перебирает весь край на КАЖДЫЙ тест).
 *   2. ВЕРШИНЫ, ВЫБРОШЕННЫЕ ЧАСТЬЮ ВЛАДЕЛЬЦЕВ — прямая мера щели, не требующая
 *      ни лучей, ни картинки, и это ОПРЕДЕЛЕНИЕ рассогласования: вершина,
 *      которую один сосед выбросил, а другой оставил, и есть щель. Считается
 *      по числу вхождений сварной вершины до и после: `0 < после < до`.
 *      (Несопряжённые рёбра — та же болезнь в другой мере, докладываются
 *      долей; АБСОЛЮТНОЕ их число для сравнения не годится, потому что
 *      упрощение уменьшает и общее число рёбер.)
 *   3. ЛУЧИ — физическое следствие. Щель видна как ПРОСКОК: луч, попадавший в
 *      поверхность, пролетает мимо и садится на дальнюю. Мерится ХВОСТ по |Δt|,
 *      а не среднее (§4), и отдельно — доля лучей, потерявших попадание вовсе.
 *
 * ЭТАЛОН ЛУЧЕЙ — НЕУПРОЩЁННЫЙ НАБОР, а не треугольники: сравнивается ровно то,
 * что упрощение меняет, и ошибка планарного представления (уже измеренная в Ш2)
 * в сравнение не входит.
 *
 * НЕГАТИВНЫЙ КОНТРОЛЬ ЗДЕСЬ — НЕЗАВИСИМЫЙ РЕЖИМ, и он обязан провалиться: тот
 * же метод, тот же допуск, та же мера расстояния, разная только единица
 * решения. Не провалился — значит согласование ничего не решает и мерили не то.
 */

#include "pedge.h"
#include "poly_seg.h"
#include "polygon.h"
#include "pray.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include "transport/cam3.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Лучей по стороне кадра. 512² = 262 144 луча по сетке `hz_pray` — десятки
 * секунд, и хвоста хватает на доли процента (Ш2 снимал 128²). */
#define NRAY 512
/* ДВА ПОРОГА, И ОНИ МЕРЯЮТ РАЗНОЕ.
 * `JUMP_M` — луч сдвинулся дальше сантиметра. Отделить этим законное движение
 * края от щели НЕЛЬЗЯ: на скользящем падении сдвиг края усиливается как
 * `1/|n·d|`, и Ш2 намерил усиление `42×`, так что миллиметровый допуск даёт
 * сантиметры законно. Величина докладывается как есть.
 * `THRU_M` — луч сел на ДРУГУЮ поверхность (дальше десяти сантиметров) либо
 * потерял попадание вовсе. Вот это уже проход СКВОЗЬ, и он законным быть не
 * может ни при каком допуске из свипа. */
#define JUMP_M 0.01
#define THRU_M 0.10

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

/* Площадь по ШНУРОВКЕ края (теорема Грина) — она и меряет, что упрощение
 * съело: моменты полигона считаны по ТРЕУГОЛЬНИКАМ и упрощением не меняются. */
static double lace_total(const hz_polyset *ps) {
  double s = 0.0;
  for (int32_t k = 0; k < ps->np; k++) {
    const hz_poly *P = &ps->p[k];
    for (int32_t l = P->l0; l < P->l0 + P->nloop; l++) {
      int32_t b = ps->loop[l], e = ps->loop[l + 1], n = e - b;
      for (int32_t q = 0; q < n; q++) {
        const double *A = ps->bv + (size_t)(b + q) * 2;
        const double *B = ps->bv + (size_t)(b + (q + 1) % n) * 2;
        s += 0.5 * (A[0] * B[1] - B[0] * A[1]);
      }
    }
  }
  return s;
}

int main(int argc, char **argv) {
  double delta = (argc > 1) ? strtod(argv[1], NULL) : 0.045;
  /* Город — второй аргумент. Нужен не ради ещё одной таблицы, а ради проверки
   * ЦЕНЫ: рабочие массивы `pedge` заведены по числу СВАРНЫХ вершин, и у
   * Rungholt их миллионы. */
  int city = (argc > 2 && strcmp(argv[2], "city") == 0);
  const double eyeh[3] = HZ_CFG_HALL_EYE, ath[3] = HZ_CFG_HALL_AT;
  const double eyec[3] = HZ_CFG_CITY_EYE, atc[3] = HZ_CFG_CITY_AT;
  hz_objmesh m;
  if (hz_obj_load(&m, city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ,
                  city ? HZ_CFG_CITY_SCALE : HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  hz_pseglist sg;
  if (hz_seg_planar(&sg, &m, delta) != 0) return 1;

  /* --- эталон: неупрощённый набор + лучи по нему --- */
  hz_polyset ref;
  if (hz_poly_build(&ref, &m, &sg) != 0) return 1;
  double lace0 = lace_total(&ref);
  int64_t ne0 = 0, ns0 = 0;
  hz_edge_unpaired(&ref, &ne0, &ns0);

  /* Вхождения сварной вершины в край ДО упрощения. Вершина с двумя и более
   * вхождениями делится соседями, и только у неё рассогласование возможно. */
  const int32_t nw = (int32_t)ref.nweld;
  int32_t *occ0 = calloc((size_t)nw, sizeof *occ0);
  int32_t *occ = calloc((size_t)nw, sizeof *occ);
  if (occ0 == NULL || occ == NULL) {
    free(occ0);
    free(occ);
    return 1;
  }
  for (int32_t i = 0; i < ref.nbv; i++)
    if (ref.bw[i] >= 0) occ0[ref.bw[i]]++;
  int64_t nshared = 0;
  for (int32_t v = 0; v < nw; v++)
    if (occ0[v] > 1) nshared++;

  const double up[3] = HZ_CFG_UP;
  tr3_camera cam;
  if (tr3_camera_look(&cam, city ? eyec : eyeh, city ? atc : ath, up, HZ_CFG_FOV_DEG * M_PI / 180.0,
                      NRAY, NRAY) != 0) {
    free(occ0);
    free(occ);
    return 1;
  }

  double *t0 = malloc((size_t)NRAY * NRAY * sizeof *t0);
  int32_t *w0 = malloc((size_t)NRAY * NRAY * sizeof *w0);
  double *dd = malloc((size_t)NRAY * NRAY * sizeof *dd);
  if (t0 == NULL || w0 == NULL || dd == NULL) {
    free(t0);
    free(w0);
    free(dd);
    free(occ0);
    free(occ);
    return 1;
  }
  {
    hz_pray g;
    if (hz_pray_build(&g, &ref, 4.0) != 0) {
      free(t0);
      free(w0);
      free(dd);
      free(occ0);
      free(occ);
      return 1;
    }
    for (int py = 0; py < NRAY; py++)
      for (int px = 0; px < NRAY; px++) {
        double o[3], d[3], t;
        tr3_camera_ray(&cam, px, py, o, d);
        int32_t who = hz_pray_hit(&g, o, d, 1e-6, &t);
        t0[py * NRAY + px] = (who >= 0) ? t : -1.0;
        w0[py * NRAY + px] = who;
      }
    hz_pray_free(&g);
  }
  int64_t nhit0 = 0;
  for (int i = 0; i < NRAY * NRAY; i++)
    if (w0[i] >= 0) nhit0++;

  printf("== упрощение края: %s, δ = %g м, полигонов %d, вершин края %d\n", city ? "ГОРОД" : "зал",
         delta, ref.np, ref.nbv);
  printf("   эталон: рёбер края %lld, одиночных %lld (%.3f%%), Σ шнуровка %.4f, "
         "лучей с попаданием %lld из %d\n",
         (long long)ne0, (long long)ns0, 100.0 * (double)ns0 / (double)ne0, lace0, (long long)nhit0,
         NRAY * NRAY);
  printf("   вершин края, делимых соседями: %lld из %lld сварных\n\n", (long long)nshared,
         (long long)ref.nweld);
  printf("   %-6s %-12s %8s %7s %11s %8s %9s %9s %9s %8s %9s\n", "допуск", "режим", "вершин",
         "на пол", "рассогл.", "одиноч.", "Δплощ", "p99|Δt|", "p99.9|Δt|", "проскок", "сквозь");

  const double tols[7] = {1e-9, 0.001, 0.005, 0.01, 0.02, 0.045, 0.15};
  for (int it = 0; it < 7; it++) {
    static const hz_edgemode md[3] = {HZ_EDGE_SHARED, HZ_EDGE_SHARED_INNER, HZ_EDGE_INDEP};
    static const char *const mn[3] = {"согласован.", "только шов", "независимо"};
    for (int mode = 0; mode < 3; mode++) {
      hz_polyset ps;
      if (hz_poly_build(&ps, &m, &sg) != 0) return 1;
      hz_edgestat es;
      double ta = now_s();
      if (hz_edge_simplify(&ps, tols[it], md[mode], &es) != 0) return 1;
      double tb = now_s();
      int64_t ne = 0, ns = 0;
      hz_edge_unpaired(&ps, &ne, &ns);
      double lace = lace_total(&ps);
      memset(occ, 0, (size_t)nw * sizeof *occ);
      for (int32_t i = 0; i < ps.nbv; i++)
        if (ps.bw[i] >= 0) occ[ps.bw[i]]++;
      int64_t ndis = 0;
      for (int32_t v = 0; v < nw; v++)
        if (occ[v] > 0 && occ[v] < occ0[v]) ndis++;

      hz_pray g;
      if (hz_pray_build(&g, &ps, 4.0) != 0) return 1;
      int64_t njump = 0, nlost = 0, nboth = 0, nthru = 0;
      for (int py = 0; py < NRAY; py++)
        for (int px = 0; px < NRAY; px++) {
          int i = py * NRAY + px;
          if (w0[i] < 0) continue;
          double o[3], d[3], t;
          tr3_camera_ray(&cam, px, py, o, d);
          int32_t who = hz_pray_hit(&g, o, d, 1e-6, &t);
          if (who < 0) {
            nlost++;
            nthru++;
            continue;
          }
          double e = fabs(t - t0[i]);
          dd[nboth++] = e;
          if (e > JUMP_M) njump++;
          if (e > THRU_M) nthru++;
        }
      hz_pray_free(&g);
      qsort(dd, (size_t)nboth, sizeof *dd, cmp_d);
      double p99 = (nboth > 0) ? dd[(int64_t)((double)nboth * 0.99)] : 0.0;
      double p999 = (nboth > 0) ? dd[(int64_t)((double)nboth * 0.999)] : 0.0;

      printf("   %-6g %-12s %8d %7.1f %11lld %7.2f%% %9.2e %9.2e %9.2e %7.3f%% %8.4f%%\n", tols[it],
             mn[mode], ps.nbv, (double)ps.nbv / ps.np, (long long)ndis,
             100.0 * (double)ns / (double)ne, fabs(lace - lace0) / fabs(lace0), p99, p999,
             100.0 * (double)njump / (double)nhit0, 100.0 * (double)nthru / (double)nhit0);
      if (mode == 0 && it == 5)
        printf("        (цепей %lld, из них замкнутых %lld; заперто %lld сварных вершин из %lld; "
               "петель ниже трёх вершин %lld; %.3f с)\n",
               (long long)es.nchain, (long long)es.ncycle, (long long)es.nlock,
               (long long)es.nlock + (long long)es.nfree, (long long)es.nloop_short, tb - ta);
      hz_poly_free(&ps);
    }
  }
  free(t0);
  free(w0);
  free(dd);
  free(occ0);
  free(occ);
  hz_poly_free(&ref);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
