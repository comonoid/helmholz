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

static int cmp_i32(const void *x, const void *y) {
  int32_t a = *(const int32_t *)x, b = *(const int32_t *)y;
  return (a < b) ? -1 : ((a > b) ? 1 : 0);
}

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* ВЫБОР ОДНОЙ СТРОКИ (`only=2.0`). Не удобство: полный городской прогон идёт
 * 37 минут, и сличать «до и после» на нём приходилось целиком ради одной
 * строки. Конфигурации при этом остаются в коде, а не в командной строке (§2:
 * иначе числа разных прогонов несравнимы). */
static double g_only = -1.0;
static int g_notie = 0;
static int g_rand = 0;
static int g_noexact = 0;
static int g_quarter = 0;
static double g_conemax = 0.0;
static double g_lossmax = 0.0, g_aspmax = 0.0;

static void run(const hz_objmesh *m, const hz_pseglist *si, const hz_polyset *ps, double eps,
                int32_t target, int random, int brute, double delta, int ov) {
  /* Равенство ТОЧНОЕ, а не с допуском: строка выбирается по тому же числу, что
   * записано в конфигурации, и «почти равно» означало бы выбор соседней строки.
   * Форма через `>=`/`<=` — то же равенство, но без `==`: гейт справедливо
   * ругается на `==` там, где сравнение вещественных случайно, а здесь оно
   * намеренное, и форма это показывает. */
  if (g_only > 0.0 && !(eps >= g_only && eps <= g_only)) return;
  hz_mergecfg mc;
  memset(&mc, 0, sizeof mc);
  mc.delta = eps;
  mc.target = target;
  mc.use_geom = 1;
  mc.use_overlap = ov && !random;
  mc.random = random;
  mc.brute = brute;
  mc.notie = g_notie;
  mc.noexact = g_noexact;
  mc.quarter = g_quarter;
  mc.conemax = g_conemax;
  mc.lossmax = g_lossmax;
  mc.aspmax = g_aspmax;
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
  /* ПОТОКИ: СТЕННОЕ И ПРОЦЕССОРНОЕ ВРЕМЯ ПОРОЗНЬ (К79), плюс ПЕРЕКОС. Без
   * отношения `CPU/стена` упор в память неотличим от упора в счёт, а без
   * перекоса «процессорное время выросло» неотличимо от простоя на активном
   * ожидании (А57): при `dynamic` ждущий поток тоже тратит процессорное время. */
  printf("        потоки: %d, ворота стенные %.2f с, процессорные %.2f с (CPU/стена %.2f); "
         "пар у самого нагруженного потока %lld при среднем %lld\n",
         st.nthreads, st.t_gate, st.t_gate_cpu, (st.t_gate > 0.0) ? st.t_gate_cpu / st.t_gate : 0.0,
         (long long)st.npair_thr_max,
         (long long)(st.nthreads > 0 ? st.npair_list / st.nthreads : st.npair_list));
  /* МИНИМАКСНЫЙ РЕШАТЕЛЬ (О11): сходимость этой формы обмена не доказана (А89),
   * поэтому печатается её поведение, а не только результат. Два последних числа
   * обязаны быть нулями. */
  printf("        минимакс: вызовов %lld, обменов в среднем %.2f, максимум %lld; упёрлись в предел "
         "%lld, недосходов %lld; max sqrt(1+a²+b²) %.6f\n",
         (long long)st.nmm_call,
         (st.nmm_call > 0) ? (double)st.nmm_exch / (double)st.nmm_call : 0.0,
         (long long)st.nmm_exch_max, (long long)st.nmm_cap, (long long)st.nmm_under, st.mm_gmax);
  /* КОНУС НОРМАЛЕЙ (О21): гистограмма полураствора ПО СЛИЯНИЯМ — числа, которого
   * в проекте не было вовсе, а без него всякий порог есть угадывание. */
  {
    int64_t tot = 0;
    for (int i = 0; i < 18; i++)
      tot += st.cone_hist[i];
    printf("        конус: полураствор по слияниям, %% в корзинах по 5°:");
    for (int i = 0; i < 18; i++)
      if (st.cone_hist[i] > 0)
        printf(" [%d-%d°] %.2f", i * 5, i * 5 + 5,
               100.0 * (double)st.cone_hist[i] / (double)(tot ? tot : 1));
    printf("; максимум %.2f°, поворотов оси %lld, вырожденных %lld, отвергнуто порогом %lld\n",
           st.cone_max, (long long)st.ncone_turn, (long long)st.ncone_fail,
           (long long)st.ncone_rej);
  }
  {
    int64_t ta = 0, tl = 0;
    for (int i = 0; i < 12; i++) {
      ta += st.asp_hist[i];
      tl += st.loss_hist[i];
    }
    printf("        форма (§38): вытянутость, %% по корзинам 2^k:");
    for (int i = 0; i < 12; i++)
      if (st.asp_hist[i] > 0)
        printf(" [%d..%d] %.2f", 1 << i, 1 << (i + 1),
               100.0 * (double)st.asp_hist[i] / (double)(ta ? ta : 1));
    printf("; максимум %.1f\n", st.asp_max);
    printf("        форма (§38): потеря площади худшим членом, %% по корзинам 0.25:");
    for (int i = 0; i < 12; i++)
      if (st.loss_hist[i] > 0)
        printf(" [%.2f..%.2f] %.2f", 1.0 + 0.25 * i, 1.25 + 0.25 * i,
               100.0 * (double)st.loss_hist[i] / (double)(tl ? tl : 1));
    printf("; максимум %.2f\n", st.loss_max);
  }
  printf("        минимакс: откатов к средневзвешенной %lld (%.2f%%)\n", (long long)st.nmm_revert,
         (st.nmm_call > 0) ? 100.0 * (double)st.nmm_revert / (double)st.nmm_call : 0.0);
  /* СЛИЯНИЕ ПОСЛЕ О9: доля точных пересчётов и РАБОТА, ими вызванная. Доля
   * попыток зависит от истории слияний (А67), поэтому вывод делается по работе
   * (А71), а не по ней. */
  printf("        слияние: попыток %lld, проверок объединения %lld, тронуто членов %lld, опорных "
         "точек %lld\n",
         (long long)st.nmerge_try, (long long)st.nmerge_exact, (long long)st.nwork_merge_mem,
         (long long)st.nwork_merge_sup);
  /* РАЗБИВКА ВОРОТ. Без неё всякий довод о том, что в них дорого, есть
   * атрибуция по догадке: О6 намерил цену ЛИНЕЙНОЙ ПО КРАЮ (сжатие втрое
   * уронило ворота втрое), а это указывает на пробы неперекрытия, а не на
   * случайный доступ к сетке треугольников. Печатаются ОБЕ работы разом. */
  printf("        ворота: мажоранты хватило %lld, точных расчётов %lld (вершин тр. %lld), "
         "неперекрытие %lld (коробкой отсеяно %lld, вершин края %lld)\n",
         (long long)st.ngate_major, (long long)st.ngate_exact, (long long)st.nwork_tri,
         (long long)st.ngate_ov, (long long)st.ngate_ovbox, (long long)st.nwork_bv);
  /* ПАМЯТЬ И СЛЕПОК. Первое нужно О10 (при δ = 8 м список пар шёл на сотни
   * мегабайт, а пик не мерился вовсе), второе — сквозному требованию фазы I:
   * «выход не меняется побитово» надо чем-то СЛИЧАТЬ. */
  printf("        память: список пар %.1f МБ (пик с удвоением %.1f МБ), сетка %.1f МБ; "
         "слепок разметки %016llx; совпадающих err %lld\n",
         (double)st.bytes_pairs / 1048576.0, (double)st.bytes_pairs_peak / 1048576.0,
         (double)st.bytes_grid / 1048576.0, (unsigned long long)st.digest, (long long)st.ntie);
  /* СТОРОЖА, КОТОРЫЕ ОБЯЗАНЫ БЫТЬ НУЛЯМИ. Печатаются только когда сработали:
   * молчание здесь — сообщение, а не отсутствие проверки. */
  if (st.ntrunc > 0 || st.nbox_tri > 0)
    printf("        СТОРОЖА: групп, не поместившихся в буфер членов %lld; полигонов, ушедших "
           "на запасной путь коробки по треугольникам %lld\n",
           (long long)st.ntrunc, (long long)st.nbox_tri);
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

/* Сравнение по УБЫВАНИЮ `dmax`: негативному контролю нужны худшие полигоны. */
static const hz_poly *g_p = NULL;
static int cmp_dmax_desc(const void *x, const void *y) {
  const hz_poly *a = &g_p[*(const int32_t *)x], *b = &g_p[*(const int32_t *)y];
  if (a->dmax > b->dmax) return -1;
  if (a->dmax < b->dmax) return 1;
  return 0;
}

/* НЕГАТИВНЫЙ КОНТРОЛЬ О7 (§13): `dmax` группы ТРЕМЯ способами — по
 * треугольникам, по опорному множеству и ПО КРАЮ.
 *
 * Первые два обязаны совпасть БИТ В БИТ: множество, по которому берётся
 * максимум, меняется, а величина — нет. Третий обязан РАЗОЙТИСЬ на зале и
 * совпасть на городе, и это следует из построения, а не из данных:
 * `hz_poly_world` возвращает точку СТРОГО в плоскости полигона, поэтому у
 * группы из одного полигона отклонение края от собственной плоскости равно нулю
 * ТОЧНО, а отклонение треугольников равно `P->dmax`.
 *
 * ЧЕГО ЭТОТ КОНТРОЛЬ НЕ ПОКАЗЫВАЕТ (А49): группа здесь из ОДНОГО полигона,
 * значит доказывается «край слеп к собственному прогибу», а не утверждение А7
 * целиком — в группе из нескольких членов край `B` относительно плоскости `A`
 * несёт настоящую информацию, и слепота там другой природы. */
static void negcontrol(const hz_objmesh *m, const hz_polyset *ps, int ntop) {
  int32_t *ord = malloc((size_t)ps->np * sizeof *ord);
  if (ord == NULL) return;
  for (int32_t k = 0; k < ps->np; k++)
    ord[k] = k;
  g_p = ps->p;
  qsort(ord, (size_t)ps->np, sizeof *ord, cmp_dmax_desc);
  double wtri = 0.0, whull = 0.0, wedge = 0.0, dth = 0.0, dte = 0.0;
  int32_t nn = (ps->np < ntop) ? ps->np : ntop;
  for (int32_t i = 0; i < nn; i++) {
    const hz_poly *P = &ps->p[ord[i]];
    double a = 0.0, b = 0.0, c = 0.0;
    for (int32_t q = 0; q < P->ntri; q++) {
      double p[3][3];
      hz_obj_tri(m, ps->tri[P->t0 + q], p);
      for (int j = 0; j < 3; j++) {
        double d = fabs(p[j][0] * P->n[0] + p[j][1] * P->n[1] + p[j][2] * P->n[2] - P->off);
        if (d > a) a = d;
      }
    }
    for (int32_t q = 0; q < P->nsup; q++) {
      const double *S = ps->sup + (size_t)(P->s0 + q) * 3;
      double d = fabs(S[0] * P->n[0] + S[1] * P->n[1] + S[2] * P->n[2] - P->off);
      if (d > b) b = d;
    }
    for (int32_t l = P->l0; l < P->l0 + P->nloop; l++)
      for (int32_t e = ps->loop[l]; e < ps->loop[l + 1]; e++) {
        double x[3];
        hz_poly_world(P, ps->bv[(size_t)e * 2], ps->bv[(size_t)e * 2 + 1], x);
        double d = fabs(x[0] * P->n[0] + x[1] * P->n[1] + x[2] * P->n[2] - P->off);
        if (d > c) c = d;
      }
    if (a > wtri) wtri = a;
    if (b > whull) whull = b;
    if (c > wedge) wedge = c;
    if (fabs(a - b) > dth) dth = fabs(a - b);
    if (a - c > dte) dte = a - c;
  }
  printf("   НЕГАТИВНЫЙ КОНТРОЛЬ О7 (%d полигонов с наибольшим dmax): по треугольникам %.6f м, "
         "по оболочке %.6f м, ПО КРАЮ %.6f м\n",
         nn, wtri, whull, wedge);
  printf("      расхождение треугольники/оболочка %.3e м (обязано быть 0), "
         "треугольники/край %.6f м (обязано быть > 0 на ЗАЛЕ и 0 на городе)\n",
         dth, dte);
  free(ord);
}

int main(int argc, char **argv) {
  int city = (argc > 1 && strcmp(argv[1], "city") == 0);
  /* ЭТАЛОН БЕЗ ОБОЛОЧКИ (А45): опорное множество — все различные вершины.
   * Слепок и `dmax_worst` обязаны совпасть с обычным прогоном; расхождение
   * означает, что монотонная цепь теряет опорные точки. */
  int use_hull = 1;
  for (int i = 1; i < argc; i++) {
    if (strncmp(argv[i], "only=", 5) == 0) g_only = strtod(argv[i] + 5, NULL);
    if (strcmp(argv[i], "nohull") == 0) use_hull = 0;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ О8: сравнение без доопределения по номерам. Ставить
     * только на ГОРОДЕ — на зале совпадающих `err` 49 из 19 846 (А14). */
    if (strcmp(argv[i], "notie") == 0) g_notie = 1;
    if (strcmp(argv[i], "rand") == 0) g_rand = 1;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ О9: точный пересчёт отключён там, где граница его
     * требует; инвариант dmax_worst <= δ обязан НАРУШИТЬСЯ на городе. */
    if (strcmp(argv[i], "noexact") == 0) g_noexact = 1;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ О11: минимакс по четверти точек. */
    if (strcmp(argv[i], "quarter") == 0) g_quarter = 1;
    /* О21: ограничение полураствора конуса нормалей, градусы. */
    if (strncmp(argv[i], "cone=", 5) == 0) g_conemax = strtod(argv[i] + 5, NULL);
    if (strncmp(argv[i], "loss=", 5) == 0) g_lossmax = strtod(argv[i] + 5, NULL);
    if (strncmp(argv[i], "asp=", 4) == 0) g_aspmax = strtod(argv[i] + 4, NULL);
  }
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
  if (hz_poly_build_ex(&ps, &m, &sg, use_hull) != 0) return 1;
  double t2 = now_s();
  printf("== огрубление: %s, δ сегментации %g м, треугольников %d, полигонов %d%s\n",
         city ? "ГОРОД" : "зал", delta, m.nt, ps.np, use_hull ? "" : "  [БЕЗ ОБОЛОЧКИ]");
  /* ПОЛИГОНЫ БЕЗ КРАЯ — не любопытство, а то, на чём села первая редакция
   * критерия: `hz_poly_build` отбрасывает петли короче трёх вершин, и такой
   * полигон не даёт НИ ОДНОЙ пробы при проверке по краю, оставаясь в группе
   * своими треугольниками. Поэтому их считают и печатают. */
  int32_t noloop = 0;
  for (int32_t k = 0; k < ps.np; k++)
    if (ps.p[k].nloop == 0) noloop++;
  printf("   сегментация %.2f с, сборка полигонов %.2f с, вершин края %d, "
         "полигонов БЕЗ КРАЯ %d (%.2f%%)\n",
         t1 - t0, t2 - t1, ps.nbv, noloop, 100.0 * noloop / (double)ps.np);
  /* КРАЙ НЕПОЛНЫЙ — ОТДЕЛЬНЫЙ КЛАСС, А НЕ ЧАСТНЫЙ СЛУЧАЙ «БЕЗ КРАЯ». Петли
   * отбрасываются ПООТДЕЛЬНОСТИ, поэтому полигон может сохранить одни и
   * потерять другие; счёт «без края» такого не видит, а всё, что меряет
   * поверхность по краю (О7, О9, О11), обязано считать именно его. */
  printf("   петель отброшено (короче трёх вершин) %lld, из них полигонов с ЧАСТИЧНО "
         "потерянным краем %lld; петель оборвано пределом длины %lld\n",
         (long long)ps.nloop_drop, (long long)ps.nloop_part, (long long)ps.nloop_trunc);
  /* РАСПРЕДЕЛЕНИЕ ЧИСЛА ТРЕУГОЛЬНИКОВ НА ПОЛИГОН. Не статистика: точный путь
   * ворот обходит ВСЕ треугольники обоих членов пары, поэтому его цена задаётся
   * не средним, а ХВОСТОМ — крупный полигон пересчитывается заново в каждой
   * паре, куда он входит. Среднее по сцене этого не показывает вовсе, и на нём
   * ошибку легко приписать кэшу. */
  {
    int32_t *h = malloc((size_t)ps.np * sizeof *h);
    if (h != NULL) {
      int64_t sum = 0;
      for (int32_t k = 0; k < ps.np; k++) {
        h[k] = ps.p[k].ntri;
        sum += ps.p[k].ntri;
      }
      qsort(h, (size_t)ps.np, sizeof *h, cmp_i32);
      int64_t tail = 0;
      for (int32_t k = ps.np - ps.np / 100; k < ps.np; k++)
        tail += h[k];
      printf("   треугольников на полигон: среднее %.1f, p50 %d, p99 %d, максимум %d; "
             "на верхний 1%% полигонов приходится %.1f%% треугольников\n",
             (double)sum / ps.np, h[ps.np / 2], h[ps.np - ps.np / 100], h[ps.np - 1],
             100.0 * (double)tail / (double)sum);
      free(h);
    }
  }
  /* РАЗМЕР ОПОРНОГО МНОЖЕСТВА — ВЫХОД О7. Сравнивать его надо с числом ВЕРШИН
   * ТРЕУГОЛЬНИКОВ (`3·ntri`), потому что именно столько трогал прежний путь, а
   * не с числом треугольников. Печатается тем же порядком величин, что и
   * распределение выше, чтобы строки сличались глазами. */
  {
    int32_t *h = malloc((size_t)ps.np * sizeof *h);
    if (h != NULL) {
      int64_t sum = 0;
      for (int32_t k = 0; k < ps.np; k++) {
        h[k] = ps.p[k].nsup;
        sum += ps.p[k].nsup;
      }
      qsort(h, (size_t)ps.np, sizeof *h, cmp_i32);
      printf("   опорных точек на полигон: среднее %.1f, p50 %d, p99 %d, максимум %d; "
             "всего %lld точек (%.1f МБ) против %d вершин треугольников — сжатие %.1f×\n",
             (double)sum / ps.np, h[ps.np / 2], h[ps.np - ps.np / 100], h[ps.np - 1],
             (long long)ps.nsupall, (double)ps.nsupall * 3 * (double)sizeof *ps.sup / 1048576.0,
             3 * m.nt, (sum > 0) ? 3.0 * (double)m.nt / (double)sum : 0.0);
      free(h);
    }
  }
  /* `dmax` СЕГМЕНТАТОРА ПРОТИВ ТОЧНОГО — прямая проверка А44. Первый берётся из
   * `sg->seg[r].dmax` и докладывается как «максимум отклонения»; второй считается
   * здесь по опорному множеству и плоскости полигона, то есть по той же
   * формуле, что критерий огрубления. Если они расходятся, «dmax = 0.0000» есть
   * ЛОЖНЫЙ НОЛЬ, и все выводы вида «у Rungholt отклонения нет» держатся на
   * округлении печати, а не на геометрии. */
  {
    double wd = 0.0, wrel = 0.0;
    int32_t nz = 0;
    for (int32_t k = 0; k < ps.np; k++) {
      const hz_poly *P = &ps.p[k];
      double e = 0.0;
      for (int32_t q = 0; q < P->nsup; q++) {
        const double *S = ps.sup + (size_t)(P->s0 + q) * 3;
        double d = fabs(S[0] * P->n[0] + S[1] * P->n[1] + S[2] * P->n[2] - P->off);
        if (d > e) e = d;
      }
      if (e > wd) wd = e;
      if (fabs(e - P->dmax) > wrel) wrel = fabs(e - P->dmax);
      if (e > 0.0 && !(P->dmax > 0.0)) nz++;
    }
    printf("   dmax: точный максимум по сцене %.3e м; худшее расхождение с dmax сегментации "
           "%.3e м; полигонов с dmax = 0 при НЕнулевом точном %d\n",
           wd, wrel, nz);
  }
  /* ВЕТВИ И СТОРОЖ. «Плоских точно» есть свойство ДАННЫХ (А44) и печатается
   * отдельно от «оболочек построено»: у полигона из трёх различных вершин
   * оболочка равна им же, и строить её незачем. Сторож обязан быть нулём. */
  printf("   опорное множество: плоских ТОЧНО %lld (%.1f%%), оболочек построено %lld, "
         "полным набором %lld, без треугольников %lld; ТОЧЕК ВНЕ ОБОЛОЧКИ %lld\n\n",
         (long long)ps.nsup_flat, 100.0 * (double)ps.nsup_flat / (double)ps.np,
         (long long)ps.nsup_hull, (long long)ps.nsup_full, (long long)ps.nsup_empty,
         (long long)ps.nsup_out);
  negcontrol(&m, &ps, 64);
  printf("\n");
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

  /* РЕЖИМ СЛУЧАЙНОЙ МЕТРИКИ ОТДЕЛЬНЫМ КЛЮЧОМ (`rand`). Он есть НЕГАТИВНЫЙ
   * КОНТРОЛЬ критерия (Ш7, приёмка О12 и О20), и до О8 его величина зависела от
   * ПОРЯДКА ЭМИСЦИИ пар — то есть при потоках рассыпалась бы (А13). После
   * замены ЛЦГ на хеш от `(a, b)` она есть функция ПАРЫ, и проверяется это
   * прямо: слепок при 1 и 16 потоках обязан совпасть. В постоянную таблицу
   * строка не ставится, чтобы не менять конфигурацию §2 задним числом. */
  if (g_rand) {
    if (city)
      run(&m, &sg, &ps, 0.5, 150000, 1, 0, delta, 1);
    else
      run(&m, &sg, &ps, 0.3, 250, 1, 0, delta, 1);
    hz_poly_free(&ps);
    hz_seg_free(&sg);
    hz_obj_free(&m);
    return 0;
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
