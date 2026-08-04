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
  int leafmax = 0, maxlev = 0, wpx = HZ_CFG_W;
  /* Полоса §191, названная пользователем: элемент от 4 до 32 пикселей. */
  double pxlo = 4.0, pxhi = 32.0;
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "d=", 2) == 0) delta = strtod(argv[i] + 2, NULL);
    if (strncmp(argv[i], "leaf=", 5) == 0) leafmax = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "lev=", 4) == 0) maxlev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "w=", 2) == 0) wpx = (int)strtol(argv[i] + 2, NULL, 10);
    if (strncmp(argv[i], "pxhi=", 5) == 0) pxhi = strtod(argv[i] + 5, NULL);
    if (strncmp(argv[i], "pxlo=", 5) == 0) pxlo = strtod(argv[i] + 5, NULL);
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
  if (hz_ptree_build(&T, &m, leafmax, maxlev) != 0) {
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

  int64_t nfail = 0;
  t0 = now_s();
#pragma omp parallel for schedule(dynamic, 8) reduction(+ : nfail)
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

  /* ---- ВЫБОР УРОВНЯ ПО ПОЛОСЕ 4…32 px (О58, §194) ---------------------
   * Спуск от корня: пока узел крупнее верхнего предела — вниз, иначе выбран.
   * Нижний предел спуском не движет, он ДИАГНОСТИКА: элементы мельче него —
   * дальние, у которых даже грубейший доступный узел мал.
   * ВЫБРАННЫМ УЗЛОМ ЕЩЁ НЕЛЬЗЯ РИСОВАТЬ: крупный узел обязан представлять всё
   * своё поддерево, а слияния снизу вверх (§190) нет. Здесь считается БЮДЖЕТ
   * КАДРА, и только он. */
  {
    double eye[3] = HZ_CFG_MIGUEL_EYE;
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
