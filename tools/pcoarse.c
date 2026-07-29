/* pcoarse — ЦЕНА ОГРУБЛЕНИЯ И ГОРОД (PLAN_ELEMENTS.md, §10, пункты 2 и 6).
 *
 * ДВА ВОПРОСА, И ОБА ПРО МАСШТАБ, А НЕ ПРО КАЧЕСТВО.
 *   (2) `hz_merge` перебирал пары, то есть стоил `O(np²)`. На зале (993
 *       полигона) это секунды, на городе (301 тыс.) — `4.5e10` пар, то есть
 *       город был недостижим не по существу, а по устройству цикла.
 *   (6) Весь Ш7 снят на зале. Огрубление на ГОРОДЕ не пробовано вовсе, а
 *       именно город — заявленная цель (CLAUDE.md: «сцена — город, не комната»).
 *
 * КАРТИНКИ ЗДЕСЬ НЕТ, И ЭТО НЕ УПУЩЕНИЕ. Метрика картинки для огрубления уже
 * признана НЕДЕЙСТВИТЕЛЬНОЙ (Ш7): при разрушенной геометрии она даже
 * улучшается, потому что эталон считает точный свет в той же уехавшей точке.
 * Связывающая величина — `dmax`, и она считается здесь. Качество картинки
 * меряет `prender ... c` на зале, где эталон осмыслен.
 */

#include "pedge.h"
#include "pmerge.h"
#include "poly_seg.h"
#include "polygon.h"
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

static void run(const hz_objmesh *m, const hz_pseglist *si, const hz_polyset *ps, double eps,
                int32_t target, int random, int brute, double delta, int ov) {
  hz_mergecfg mc;
  memset(&mc, 0, sizeof mc);
  mc.delta = eps;
  mc.target = target;
  mc.use_geom = 1;
  mc.use_overlap = ov && !random;
  mc.random = random;
  mc.brute = brute;
  hz_pseglist so;
  hz_mergestat st;
  double t0 = now_s();
  if (hz_merge(&so, m, si, ps, &mc, &st) != 0) {
    printf("   ОТКАЗ hz_merge\n");
    return;
  }
  double t1 = now_s();
  printf("   %5.2f %8d %7s %8d %9.3e %6.2f %10lld %9lld %9lld %8lld %10lld %8.2f\n", eps, target,
         random ? "случ" : (brute ? "перебор" : "сетка"), st.nseg_out, st.dmax_worst,
         st.dmax_worst / delta, (long long)st.ncand, (long long)st.npair, (long long)st.nmerged,
         (long long)st.ngridcell, (long long)st.ngrident, t1 - t0);
  /* ГДЕ ИМЕННО ВРЕМЯ. «Огрубление идёт семь минут» — жалоба, а не адрес. */
  printf("        по фазам, с: сетка %.2f, ВОРОТА %.2f, слияние %.2f, разметка %.2f "
         "(на пару-кандидата %.2f мкс)\n",
         st.t_grid, st.t_gate, st.t_merge, st.t_label,
         (st.ncand > 0) ? 1e6 * st.t_gate / (double)st.ncand : 0.0);
  /* ИНВАРИАНТ, А НЕ УКРАШЕНИЕ: критерий обязан держать `dmax` под допуском
   * ОГРУБЛЕНИЯ. Нарушение печатается с УЛИКАМИ — сколько групп вышло за допуск и
   * из скольких ЧЛЕНОВ состоит худшая, — потому что «сколько-то больше» не
   * говорит, откуда: одиночная группа означала бы ошибку сегментации, крупная —
   * ошибку проверки при слиянии. */
  if (st.dmax_worst > eps) {
    int32_t *nmem = calloc((size_t)so.nseg, sizeof *nmem);
    int32_t *seen = malloc((size_t)si->nseg * sizeof *seen);
    if (nmem != NULL && seen != NULL) {
      for (int32_t k = 0; k < si->nseg; k++)
        seen[k] = -1;
      for (int32_t t = 0; t < m->nt; t++) {
        int32_t k = si->label[t];
        if (seen[k] < 0) {
          seen[k] = so.label[t];
          nmem[so.label[t]]++;
        }
      }
      int32_t worst = 0, nover = 0;
      for (int32_t g = 0; g < so.nseg; g++) {
        if (so.seg[g].dmax > eps) nover++;
        if (so.seg[g].dmax > so.seg[worst].dmax) worst = g;
      }
      printf("        ИНВАРИАНТ НАРУШЕН: групп с dmax > %g — %d из %d; худшая: членов %d, "
             "треугольников %d, площадь %.2f м², dmax %.3f м\n",
             eps, nover, so.nseg, nmem[worst], so.seg[worst].ntri, so.seg[worst].area,
             so.seg[worst].dmax);
      /* УЛИКИ ПО ХУДШЕЙ ГРУППЕ. Вопрос ровно один: проверка при слиянии
       * пропустила эту группу — или проверила ДРУГУЮ плоскость? Поэтому плоскость
       * пересчитывается здесь ЗАНОВО, той же формулой, и печатаются обе. */
      double s = 0.0, ns[3] = {0, 0, 0}, org[3] = {0, 0, 0};
      printf("        члены худшей группы (участок: площадь, нормаль, off):\n");
      for (int32_t k = 0; k < si->nseg; k++) {
        if (seen[k] != worst) continue;
        const hz_poly *P = &ps->p[k];
        s += P->area;
        for (int c = 0; c < 3; c++) {
          ns[c] += P->area * P->n[c];
          org[c] += P->area * P->org[c];
        }
        printf("          %7d: A %8.4f  n (%+.4f %+.4f %+.4f)  off %+9.4f  ntri %d\n", k, P->area,
               P->n[0], P->n[1], P->n[2], P->off, P->ntri);
      }
      double nn = sqrt(ns[0] * ns[0] + ns[1] * ns[1] + ns[2] * ns[2]);
      printf("        |Σ A·n| = %.6f при Σ A = %.6f  (отношение %.6f — единица значит "
             "согласованные нормали, ноль — взаимно погасившиеся)\n",
             nn, s, (s > 0.0) ? nn / s : 0.0);
      if (nn > 0.0 && s > 0.0) {
        double n2[3], off2 = 0.0;
        for (int c = 0; c < 3; c++) {
          n2[c] = ns[c] / nn;
          org[c] /= s;
          off2 += n2[c] * org[c];
        }
        double dm2 = 0.0;
        for (int32_t t = 0; t < m->nt; t++) {
          if (so.label[t] != worst) continue;
          double p[3][3];
          hz_obj_tri(m, t, p);
          for (int j = 0; j < 3; j++) {
            double dv = fabs(p[j][0] * n2[0] + p[j][1] * n2[1] + p[j][2] * n2[2] - off2);
            if (dv > dm2) dm2 = dv;
          }
        }
        printf("        плоскость группы: доклад n (%+.4f %+.4f %+.4f) off %+9.4f; "
               "пересчёт n (%+.4f %+.4f %+.4f) off %+9.4f; dmax по пересчёту %.4f м\n",
               so.seg[worst].n[0], so.seg[worst].n[1], so.seg[worst].n[2], so.seg[worst].off, n2[0],
               n2[1], n2[2], off2, dm2);
      }
    }
    free(nmem);
    free(seen);
  }
  fflush(stdout);
  hz_seg_free(&so);
}

int main(int argc, char **argv) {
  int city = (argc > 1 && strcmp(argv[1], "city") == 0);
  double delta = (argc > 2) ? strtod(argv[2], NULL) : (city ? 0.05 : 0.045);
  hz_objmesh m;
  if (hz_obj_load(&m, city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ,
                  city ? HZ_CFG_CITY_SCALE : HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  double t0 = now_s();
  hz_pseglist sg;
  if (hz_seg_planar(&sg, &m, delta) != 0) return 1;
  double t1 = now_s();
  hz_polyset ps;
  if (hz_poly_build(&ps, &m, &sg) != 0) return 1;
  double t2 = now_s();
  printf("== огрубление: %s, δ сегментации %g м, треугольников %d, полигонов %d\n",
         city ? "ГОРОД" : "зал", delta, m.nt, ps.np);
  /* ПОЛИГОНЫ БЕЗ КРАЯ — не любопытство, а то, на чём села первая редакция
   * критерия: `hz_poly_build` отбрасывает петли короче трёх вершин, и такой
   * полигон не даёт НИ ОДНОЙ пробы при проверке по краю, оставаясь в группе
   * своими треугольниками. Поэтому их считают и печатают. */
  int32_t noloop = 0;
  for (int32_t k = 0; k < ps.np; k++)
    if (ps.p[k].nloop == 0) noloop++;
  printf("   сегментация %.2f с, сборка полигонов %.2f с, вершин края %d, "
         "полигонов БЕЗ КРАЯ %d (%.2f%%)\n\n",
         t1 - t0, t2 - t1, ps.nbv, noloop, 100.0 * noloop / (double)ps.np);
  printf("   %5s %8s %7s %8s %9s %6s %10s %9s %9s %8s %10s %8s\n", "δ огр", "цель", "отбор",
         "полиг", "dmax", "/δ", "кандидат", "пар", "слито", "ячеек", "вхожд", "с");

  /* КРАЙ, СЖАТЫЙ СОГЛАСОВАННО (О1), — и это не украшение таблицы. Проверка
   * НЕПЕРЕКРЫТИЯ ставит 64 пробы «точка внутри края» на каждую пару, а
   * `hz_poly_inside` обходит ВЕСЬ край; значит цена огрубления линейна по числу
   * вершин края, и упрощение края обязано её уронить во столько же раз, во
   * сколько уронило край. Допуск вчетверо мельче δ сегментации — та же точка,
   * что в О1. */
  if (argc > 3 && strcmp(argv[3], "simp") == 0) {
    hz_edgestat es;
    double ta = now_s();
    if (hz_edge_simplify(&ps, 0.25 * delta, HZ_EDGE_SHARED, &es) != 0) return 1;
    printf("   край сжат согласованно за %.3f с: %lld -> %lld вершин (%.1f на полигон)\n\n",
           now_s() - ta, (long long)es.nbv_in, (long long)es.nbv_out, (double)es.nbv_out / ps.np);
  }

  if (city) {
    /* Допуски огрубления заданы масштабом города, а не залом: у Rungholt
     * сегментация при δ ≤ 0.2 м не аппроксимирует НИЧЕГО (Ш1: `max dmax = 0`
     * точно), поэтому огрублять имеет смысл начиная с метров. */
    run(&m, &sg, &ps, 0.5, 150000, 0, 0, delta, 1);
    run(&m, &sg, &ps, 2.0, 20000, 0, 0, delta, 1);
    run(&m, &sg, &ps, 8.0, 20000, 0, 0, delta, 1);
    /* ТА ЖЕ ТОЧКА БЕЗ ПРОВЕРКИ ПЕРЕКРЫТИЯ — чтобы разделить цену ворот надвое. */
    run(&m, &sg, &ps, 0.5, 150000, 0, 0, delta, 0);
  } else {
    run(&m, &sg, &ps, 0.3, 250, 0, 0, delta, 0);
    const int32_t tg[4] = {700, 450, 250, 150};
    for (int j = 0; j < 4; j++)
      run(&m, &sg, &ps, 0.3, tg[j], 0, 0, delta, 1);
    /* ЭТАЛОН СЕТКИ — ПЕРЕБОР. На зале он ещё выполним, и только поэтому
     * равенство проверяемо вообще. */
    for (int j = 0; j < 4; j++)
      run(&m, &sg, &ps, 0.3, tg[j], 0, 1, delta, 1);
  }
  hz_poly_free(&ps);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
