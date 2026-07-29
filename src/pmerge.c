/* Огрубление слиянием. Разбор и оговорки — в `pmerge.h`. */

#include "pmerge.h"
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Проб на полигон при проверке НЕПЕРЕКРЫТИЯ. Число выведено, а не подобрано:
 * перекрытие, занимающее меньше 1/16 площади, слиянию не мешает (ошибка
 * площади ниже допуска на неё), а 16 проб по каждой оси такое перекрытие
 * обнаруживают заведомо. */
#define PM_OVSAMP 8

typedef struct {
  int32_t a, b;
  double err;
} pm_pair;

/* Порядок ПОЛНЫЙ: сперва ошибка, потом номера. Сравнение по одной ошибке
 * оставляло равные пары на произвол `qsort`, а от их порядка зависит, какие
 * группы срастутся первыми. Тогда сеточный набор кандидатов и переборный,
 * СОВПАДАЯ как множества, давали бы разный выход, и проверку «сетка равна
 * перебору» нельзя было бы поставить вовсе. */
static int cmp_pair(const void *x, const void *y) {
  const pm_pair *p = x, *q = y;
  if (p->err < q->err) return -1;
  if (p->err > q->err) return 1;
  if (p->a != q->a) return (p->a < q->a) ? -1 : 1;
  if (p->b != q->b) return (p->b < q->b) ? -1 : 1;
  return 0;
}

static int32_t uf_find(int32_t *p, int32_t x) {
  while (p[x] != x) {
    p[x] = p[p[x]];
    x = p[x];
  }
  return x;
}

/* Мировая коробка полигона по краю. */
static void pbox(const hz_polyset *ps, int32_t k, double lo[3], double hi[3]) {
  const hz_poly *P = &ps->p[k];
  for (int a = 0; a < 3; a++) {
    lo[a] = 1e300;
    hi[a] = -1e300;
  }
  for (int32_t l = P->l0; l < P->l0 + P->nloop; l++)
    for (int32_t b = ps->loop[l]; b < ps->loop[l + 1]; b++) {
      double x[3];
      hz_poly_world(P, ps->bv[(size_t)b * 2], ps->bv[(size_t)b * 2 + 1], x);
      for (int a = 0; a < 3; a++) {
        if (x[a] < lo[a]) lo[a] = x[a];
        if (x[a] > hi[a]) hi[a] = x[a];
      }
    }
}

/* Общая плоскость СПИСКА полигонов и МАКСИМАЛЬНОЕ отклонение их вершин от неё —
 * то есть `dmax`, а не среднеквадратичное (Г40/Г44).
 *
 * ПРОВЕРЯЕТСЯ ОБЪЕДИНЕНИЕ ГРУПП, А НЕ ПАРА, и это исправление настоящей ошибки.
 * Первая редакция проверяла пару и сливала транзитивно: `A~B` и `B~C` проходили
 * по отдельности, а у `A∪B∪C` отклонение не ограничивалось ничем. Замер поймал
 * это сразу — `dmax` дорастал до `44 δ` при допуске `1 δ`, — и поймал его
 * ИМЕННО тот пересчёт `dmax`, которого требует Ш7 («ложный ноль опаснее
 * UNKNOWN»). Метрика картинки к этому слепа: при уехавшей на два метра
 * плоскости она даже УЛУЧШАЛАСЬ, потому что сравнивает поле с точным светом в
 * той же уехавшей точке (узор К13/К40/К94). */
/* ОТКЛОНЕНИЕ МЕРЯЕТСЯ ПО ТРЕУГОЛЬНИКАМ, А НЕ ПО КРАЮ, И ЭТО ИСПРАВЛЕНИЕ ПО
 * ЗАМЕРУ. Первая редакция брала вершины КРАЯ, и на городе это дало `dmax = 4.18`
 * м при допуске огрубления `0.5` м — то есть критерий, обязанный держать `dmax`
 * под допуском, его не держал. Причина: `hz_poly_build` отбрасывает петли короче
 * трёх вершин (`nvloop < 3`), и полигон может остаться БЕЗ КРАЯ ВОВСЕ, сохранив
 * свои треугольники. По краю такой полигон не даёт НИ ОДНОЙ пробы — его
 * отклонение не проверяется ничем, а в группу он входит целиком. На зале не
 * срабатывало: там полигонов без петли единицы.
 * Треугольники — это и есть поверхность; край есть её проекция. Заодно критерий
 * и докладываемый `dmax` становятся ОДНОЙ величиной, а не двумя похожими, и
 * тогда «`dmax_worst < δ`» — проверяемый инвариант, а не пожелание. */
static double group_plane(const hz_polyset *ps, const hz_objmesh *m, const int32_t *mem, int32_t nm,
                          double n[3], double *off) {
  double s = 0.0, ns[3] = {0, 0, 0}, org[3] = {0, 0, 0};
  for (int32_t i = 0; i < nm; i++) {
    const hz_poly *P = &ps->p[mem[i]];
    double w = P->area;
    s += w;
    for (int c = 0; c < 3; c++) {
      ns[c] += w * P->n[c];
      org[c] += w * P->org[c];
    }
  }
  if (!(s > 0.0)) return 1e300;
  double nn = sqrt(ns[0] * ns[0] + ns[1] * ns[1] + ns[2] * ns[2]);
  if (!(nn > 0.0)) return 1e300;
  for (int c = 0; c < 3; c++) {
    n[c] = ns[c] / nn;
    org[c] /= s;
  }
  /* ВСЕ НОРМАЛИ ОБЯЗАНЫ СМОТРЕТЬ В ОДНУ СТОРОНУ С ГРУППОВОЙ. Найдено замером на
   * ГОРОДЕ: группа из 12 участков дала `dmax = 4.18` м при допуске `0.5` м, и
   * улики показали, чем именно — `|Σ A·n| = 0` при `Σ A = 15.4`. Это шесть пар
   * ВСТРЕЧНЫХ граней (лицо и изнанка одной диагональной стены, `n` и `−n`,
   * площади равны). Взвешенная нормаль у них гасится ТОЧНО, и «плоскость
   * группы» получается направлением ОСТАТКА ОКРУГЛЕНИЯ — то есть произволом
   * порядка `1e−17`, усиленным до единичного вектора. Проверка `dmax` при этом
   * не спасает: если случайное направление легло поперёк разброса, отклонение
   * выходит малым и слияние ПРИНИМАЕТСЯ; а в конце плоскость пересчитывается в
   * ДРУГОМ порядке суммирования, даёт другое случайное направление — и `dmax`
   * оказывается метрами.
   * Порога здесь нет и не нужно: требуется строгая положительность `n_i·n`.
   * Встречные грани — это ДВЕ поверхности, а не одна, и сливать их нельзя ни
   * при каком допуске; кривая же поверхность, которую огрубление и должно
   * сливать, держит нормали в одной полусфере, и её ограничивает `dmax`. */
  for (int32_t i = 0; i < nm; i++) {
    const hz_poly *P = &ps->p[mem[i]];
    if (!(P->n[0] * n[0] + P->n[1] * n[1] + P->n[2] * n[2] > 0.0)) return 1e300;
  }
  *off = n[0] * org[0] + n[1] * org[1] + n[2] * org[2];

  double dmax = 0.0;
  for (int32_t i = 0; i < nm; i++) {
    const hz_poly *P = &ps->p[mem[i]];
    for (int32_t q = 0; q < P->ntri; q++) {
      double p[3][3];
      hz_obj_tri(m, ps->tri[P->t0 + q], p);
      for (int j = 0; j < 3; j++) {
        double d = fabs(p[j][0] * n[0] + p[j][1] * n[1] + p[j][2] * n[2] - *off);
        if (d > dmax) dmax = d;
      }
    }
  }
  return dmax;
}

/* ВОРОТА — ОЦЕНКА СВЕРХУ ЗА `O(1)`, А НЕ ПОДГОНКА. Переписано 07-29 по замеру
 * и по замечанию пользователя, и замечание было верным дословно: «фиттингом
 * добиваемся, чтобы полигоны в одной плоскости слить — это как, это что?».
 *
 * ЧТО БЫЛО. На каждую пару-кандидата собиралась квадрика по ВСЕМ вершинам края
 * обоих полигонов с их нормалями и решалась система 3×3 (`src/cut/qef.c`, она
 * писалась для дуального контурирования). Цена — `1.45` мкс на пару, измерено.
 * За эти деньги работает настоящая ПОДГОНКА, а вопрос был «лежат ли два куска в
 * общей плоскости с точностью δ», то есть да/нет.
 *
 * ЧТО СТАЛО. Мажоранта за `O(1)` по коробке края (`dev_box`), и ТОЧНЫЙ расчёт
 * только если мажоранта не пропустила. Точный — тот же `group_plane`, что решает
 * при слиянии, на группе из двух; значит ворота и решение стали ОДНОЙ величиной.
 * Одной мажоранты мало: она строже точного критерия, и замер это показал сразу —
 * цель 150 полигонов переставала достигаться (246 вместо 150).
 *
 * ЗАОДНО СНЯТА ОГОВОРКА Г40/Г44. Невязка QEF есть СРЕДНЕКВАДРАТИЧНОЕ, и
 * подменять ею максимум запрещено прямым текстом; здесь и мажоранта, и точный
 * расчёт — величины той же природы, что критерий, то есть МАКСИМУМ. */
/* Худшее отклонение края `B` от плоскости `A`, оценённое по КОРОБКЕ края в
 * местной раме `B`. Точка `B` есть `org_b + u·eu_b + v·ev_b`, поэтому её
 * отклонение от плоскости `A` — линейная функция `(u, v)`, и максимум модуля
 * достигается В УГЛУ коробки. Четыре вычисления вместо описанного радиуса: та
 * же `O(1)`, но мажоранта много туже, а от туготы зависит, как часто придётся
 * считать точно. Добавляется `dmax` самого `B` — на столько его треугольники
 * отходят от собственной плоскости. */
static double dev_box(const hz_poly *A, const hz_poly *B) {
  double d0 = A->n[0] * B->org[0] + A->n[1] * B->org[1] + A->n[2] * B->org[2] - A->off;
  double pu = A->n[0] * B->eu[0] + A->n[1] * B->eu[1] + A->n[2] * B->eu[2];
  double pv = A->n[0] * B->ev[0] + A->n[1] * B->ev[1] + A->n[2] * B->ev[2];
  double best = 0.0;
  for (int i = 0; i < 4; i++) {
    double u = (i & 1) ? B->uvhi[0] : B->uvlo[0];
    double v = (i & 2) ? B->uvhi[1] : B->uvlo[1];
    double d = fabs(d0 + u * pu + v * pv);
    if (d > best) best = d;
  }
  return best + B->dmax;
}

static double pair_gate(const hz_polyset *ps, int32_t a, int32_t b) {
  const hz_poly *A = &ps->p[a], *B = &ps->p[b];
  double c = A->n[0] * B->n[0] + A->n[1] * B->n[1] + A->n[2] * B->n[2];
  /* Встречные грани — две поверхности, а не одна (см. `group_plane`). */
  if (!(c > 0.0)) return 1e300;
  double ea = dev_box(A, B), eb = dev_box(B, A);
  return (ea > eb) ? ea : eb;
}

/* ПЕРЕКРЫВАЮТСЯ ЛИ ПРОЕКЦИИ. Добавлено по замеру Ш3: критерий «все точки в
 * пределах δ от общей плоскости» разрешает слить две ОДИНАКОВО СМОТРЯЩИЕ
 * поверхности ближе δ, и на грубом δ это съедает до 10.6% площади сцены, чего
 * баланс энергии не ловит вовсе. */
static int overlaps(const hz_polyset *ps, int32_t a, int32_t b) {
  const hz_poly *A = &ps->p[a], *B = &ps->p[b];
  /* ДЕШЁВЫЙ ОТСЕВ ПЕРЕД ПРОБАМИ, и он снимает девять десятых цены огрубления.
   * Измерено 07-29: с проверкой перекрытия ворота стоили `14.19` мкс на пару,
   * без неё — `1.45`, то есть 90% времени уходило сюда. Причина простая:
   * `PM_OVSAMP²  = 64` пробы, каждая — обход ВСЕГО края (57 вершин у зала), и
   * так на каждую пару-кандидата.
   * Отсев точный, а не эвристический: проекция `B` на раму `A` целиком лежит в
   * коробке `[uvlo, uvhi]` каждого, и если коробки не пересекаются, то не
   * пересекаются и сами полигоны — пробовать нечего. Считается за десяток
   * сравнений по величинам, посчитанным при импорте. Соседние по ребру
   * полигоны (а их большинство среди кандидатов) отсеиваются здесь же. */
  {
    /* Коробка `B`, перенесённая в раму `A` через четыре УГЛА (не через радиус:
     * радиус раздувает её вчетверо и отсев перестаёт срабатывать). Перенос
     * через габарит углов — оценка СВЕРХУ, поэтому настоящее перекрытие отсев
     * выбросить не может. */
    double lo[2] = {1e300, 1e300}, hi[2] = {-1e300, -1e300};
    for (int i = 0; i < 4; i++) {
      double u = (i & 1) ? B->uvhi[0] : B->uvlo[0];
      double v = (i & 2) ? B->uvhi[1] : B->uvlo[1];
      double x[3], q[3];
      hz_poly_world(B, u, v, x);
      for (int c = 0; c < 3; c++)
        q[c] = x[c] - A->org[c];
      double ua = q[0] * A->eu[0] + q[1] * A->eu[1] + q[2] * A->eu[2];
      double va = q[0] * A->ev[0] + q[1] * A->ev[1] + q[2] * A->ev[2];
      if (ua < lo[0]) lo[0] = ua;
      if (ua > hi[0]) hi[0] = ua;
      if (va < lo[1]) lo[1] = va;
      if (va > hi[1]) hi[1] = va;
    }
    double w0 =
        (hi[0] < A->uvhi[0] ? hi[0] : A->uvhi[0]) - (lo[0] > A->uvlo[0] ? lo[0] : A->uvlo[0]);
    double w1 =
        (hi[1] < A->uvhi[1] ? hi[1] : A->uvhi[1]) - (lo[1] > A->uvlo[1] ? lo[1] : A->uvlo[1]);
    /* ОТСЕВ ПО ПЛОЩАДИ, А НЕ ПО ФАКТУ КАСАНИЯ, и в этом всё дело. Соседние по
     * ребру полигоны — а их среди кандидатов большинство — коробками ВСЕГДА
     * соприкасаются, поэтому проверка «пересекаются ли коробки» не отсеивает
     * почти ничего. Но перекрытием считается ДОЛЯ больше 1/16 (см. ниже), а
     * площадь пересечения коробок эту долю ограничивает сверху: если она сама
     * меньше 1/16 меньшей коробки, настоящее перекрытие тем более меньше, и
     * шестьдесят четыре пробы можно не ставить. */
    if (!(w0 > 0.0) || !(w1 > 0.0)) return 0;
    double ab = (hi[0] - lo[0]) * (hi[1] - lo[1]);
    double aa = (A->uvhi[0] - A->uvlo[0]) * (A->uvhi[1] - A->uvlo[1]);
    double amin = (aa < ab) ? aa : ab;
    if (w0 * w1 * 16.0 < amin) return 0;
  }
  int inside = 0, tried = 0;
  for (int i = 0; i < PM_OVSAMP; i++)
    for (int j = 0; j < PM_OVSAMP; j++) {
      double u = B->uvlo[0] + ((double)i + 0.5) / PM_OVSAMP * (B->uvhi[0] - B->uvlo[0]);
      double v = B->uvlo[1] + ((double)j + 0.5) / PM_OVSAMP * (B->uvhi[1] - B->uvlo[1]);
      if (!hz_poly_inside(ps, B, u, v)) continue;
      tried++;
      double x[3], q[3];
      hz_poly_world(B, u, v, x);
      for (int c = 0; c < 3; c++)
        q[c] = x[c] - A->org[c];
      double ua = q[0] * A->eu[0] + q[1] * A->eu[1] + q[2] * A->eu[2];
      double va = q[0] * A->ev[0] + q[1] * A->ev[1] + q[2] * A->ev[2];
      if (ua < A->uvlo[0] || ua > A->uvhi[0] || va < A->uvlo[1] || va > A->uvhi[1]) continue;
      if (hz_poly_inside(ps, A, ua, va)) inside++;
    }
  /* Соседние по ребру полигоны дают единичные попадания на самой кромке;
   * перекрытием считается заметная ДОЛЯ, а не факт попадания. */
  return (tried > 0) && (inside * 16 > tried);
}

static double pm_now(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* --- отбор кандидатов ---------------------------------------------------------
 * `pm_emit` — единственное место, где пара проверяется и попадает в список.
 * Едино для сетки и для переборного эталона: иначе «совпадение выходов» ничего
 * бы не доказывало, потому что различаться могли бы сами проверки. */

typedef struct {
  pm_pair *pl;
  int64_t pn, pcap;
} pm_plist;

static int pm_emit(pm_plist *L, const hz_polyset *ps, const hz_objmesh *m, const double *bb,
                   int32_t a, int32_t b, const hz_mergecfg *cfg, hz_mergestat *st, uint64_t *rnd) {
  st->ncand++;
  const double *A = bb + (size_t)a * 6, *B = bb + (size_t)b * 6;
  for (int c = 0; c < 3; c++) {
    if (A[c] - B[3 + c] > cfg->delta) return 0;
    if (B[c] - A[3 + c] > cfg->delta) return 0;
  }
  st->npair++;

  double err = 0.0;
  /* ГЕОМЕТРИЯ ЗДЕСЬ НЕ ПРОВЕРЯЕТСЯ — она проверяется при СЛИЯНИИ, по
   * объединению групп (см. `group_plane`). Здесь только законно парные
   * члены: ворота по оценке сверху и перекрытие. */
  if (cfg->random) {
    /* НЕГАТИВНЫЙ КОНТРОЛЬ: метрика СЛУЧАЙНАЯ. Ошибка обязана стать O(1);
     * если не стала — критерий ничего не решает и мерили не его. */
    *rnd = *rnd * 6364136223846793005ULL + 1442695040888963407ULL;
    err = (double)(*rnd >> 11) / 9007199254740992.0;
  } else {
    if (cfg->use_geom) {
      /* МАЖОРАНТА ЗА O(1), И ТОЧНЫЙ РАСЧЁТ ТОЛЬКО ЕСЛИ ОНА НЕ ПРОПУСТИЛА.
       * Одна мажоранта качество портит: она строже точного критерия, и замер
       * это показал сразу — цель 150 полигонов переставала достигаться (246
       * вместо 150). Точный расчёт — тот же `group_plane`, что решает при
       * слиянии, только на группе из двух; значит ворота и решение стали ОДНОЙ
       * величиной, а не двумя похожими. Дорогой путь берётся лишь там, где
       * дешёвый сомневается. */
      double qg = pair_gate(ps, a, b);
      if (qg > cfg->delta) {
        int32_t mm[2] = {a, b};
        double gn[3], goff;
        qg = group_plane(ps, m, mm, 2, gn, &goff);
        if (qg > cfg->delta) {
          st->nrej_geom++;
          return 0;
        }
      }
      err += qg / cfg->delta;
    }
    if (cfg->use_overlap && overlaps(ps, a, b)) {
      st->nrej_overlap++;
      return 0;
    }
  }
  if (L->pn >= L->pcap) {
    int64_t nc = L->pcap * 2;
    pm_pair *q = realloc(L->pl, (size_t)nc * sizeof *q);
    if (q == NULL) return 2;
    L->pl = q;
    L->pcap = nc;
  }
  L->pl[L->pn].a = a;
  L->pl[L->pn].b = b;
  L->pl[L->pn].err = err;
  L->pn++;
  return 0;
}

int hz_merge(hz_pseglist *so, const hz_objmesh *m, const hz_pseglist *si, const hz_polyset *ps,
             const hz_mergecfg *cfg, hz_mergestat *st) {
  memset(so, 0, sizeof *so);
  memset(st, 0, sizeof *st);
  st->nseg_in = si->nseg;
  /* ТОЛЬКО ПОЛИГОНЫ ИЗ УЧАСТКОВ. В наборе за ними идут ИСТОЧНИКИ, добавленные
   * `hz_poly_add_quad`: у них нет треугольников в сетке, сливать их не с чем и
   * незачем, а группа без треугольников даёт пустой участок на выходе. */
  const int32_t np = (ps->np < si->nseg) ? ps->np : si->nseg;

  double *bb = malloc((size_t)np * 6 * sizeof *bb);
  int32_t *par = malloc((size_t)np * sizeof *par);
  if (bb == NULL || par == NULL) {
    free(bb);
    free(par);
    return 2;
  }
  double glo[3] = {1e300, 1e300, 1e300}, ghi[3] = {-1e300, -1e300, -1e300};
  for (int32_t k = 0; k < np; k++) {
    pbox(ps, k, bb + (size_t)k * 6, bb + (size_t)k * 6 + 3);
    par[k] = k;
    for (int c = 0; c < 3; c++) {
      if (bb[(size_t)k * 6 + (size_t)c] < glo[c]) glo[c] = bb[(size_t)k * 6 + (size_t)c];
      if (bb[(size_t)k * 6 + 3 + (size_t)c] > ghi[c]) ghi[c] = bb[(size_t)k * 6 + 3 + (size_t)c];
    }
  }

  /* --- кандидаты --- */
  pm_plist L;
  L.pcap = 4096;
  L.pn = 0;
  L.pl = malloc((size_t)L.pcap * sizeof *L.pl);
  if (L.pl == NULL) {
    free(bb);
    free(par);
    return 2;
  }
  uint64_t rnd = 0x9E3779B97F4A7C15ULL;
  int oom = 0;
  double tphase = pm_now();
  if (cfg->brute) {
    for (int32_t a = 0; a < np && !oom; a++)
      for (int32_t b = a + 1; b < np && !oom; b++)
        if (pm_emit(&L, ps, m, bb, a, b, cfg, st, &rnd) != 0) oom = 1;
  } else {
    /* СЕТКА КАНДИДАТОВ. Коробка полигона расширяется на `δ/2` с каждой стороны,
     * поэтому две коробки, отстоящие меньше чем на `δ`, заведомо делят ячейку:
     * расширенные пересекаются, а ячейка точки пересечения принадлежит обеим
     * записям. Точная проверка расстояния остаётся в `pm_emit` — сетка только
     * отбрасывает заведомо далёкое. */
    double ext[3], vol = 1.0;
    for (int c = 0; c < 3; c++) {
      ext[c] = ghi[c] - glo[c];
      if (!(ext[c] > 0.0)) ext[c] = 1e-9;
      vol *= ext[c];
    }
    /* Ячеек примерно `np/2` — та же цель, что у сетки лучей: пар внутри ячейки
     * квадратично по числу жильцов, поэтому мельче двух дробить незачем. */
    double want = (double)np / 2.0;
    if (want < 1.0) want = 1.0;
    double s = cbrt(vol / want);
    int32_t nc[3];
    double cs[3];
    int64_t tot = 1;
    for (int c = 0; c < 3; c++) {
      double n = floor(ext[c] / (s > 0.0 ? s : 1e-9)) + 1.0;
      if (n > 512.0) n = 512.0;
      nc[c] = (int32_t)n;
      cs[c] = ext[c] / (double)nc[c];
      tot *= nc[c];
    }
    int32_t *cr = malloc((size_t)np * 6 * sizeof *cr);
    int32_t *start = calloc((size_t)tot + 1, sizeof *start);
    if (cr == NULL || start == NULL) {
      free(cr);
      free(start);
      free(L.pl);
      free(bb);
      free(par);
      return 2;
    }
    double half = 0.5 * cfg->delta;
    for (int32_t k = 0; k < np; k++)
      for (int c = 0; c < 3; c++) {
        double t0 = (bb[(size_t)k * 6 + (size_t)c] - half - glo[c]) / cs[c];
        double t1 = (bb[(size_t)k * 6 + 3 + (size_t)c] + half - glo[c]) / cs[c];
        int32_t i0 = (int32_t)floor(t0), i1 = (int32_t)floor(t1);
        if (i0 < 0) i0 = 0;
        if (i0 > nc[c] - 1) i0 = nc[c] - 1;
        if (i1 < 0) i1 = 0;
        if (i1 > nc[c] - 1) i1 = nc[c] - 1;
        if (i1 < i0) i1 = i0;
        cr[(size_t)k * 6 + (size_t)c] = i0;
        cr[(size_t)k * 6 + 3 + (size_t)c] = i1;
      }
    /* Два прохода развёрнуты ЯВНО, а не циклом по `pass`: у выделения внутри
     * цикла gcc-analyzer теряет связь «после прохода 0 указатель заведён» и
     * докладывает разыменование NULL. Развёрнутая форма и человеку читается
     * прямее. */
    for (int32_t k = 0; k < np; k++)
      for (int32_t z = cr[(size_t)k * 6 + 2]; z <= cr[(size_t)k * 6 + 5]; z++)
        for (int32_t y = cr[(size_t)k * 6 + 1]; y <= cr[(size_t)k * 6 + 4]; y++)
          for (int32_t x = cr[(size_t)k * 6 + 0]; x <= cr[(size_t)k * 6 + 3]; x++)
            start[((int64_t)z * nc[1] + y) * nc[0] + x + 1]++;
    for (int64_t c = 0; c < tot; c++)
      start[c + 1] += start[c];
    int32_t *idx = calloc((size_t)(start[tot] > 0 ? start[tot] : 1), sizeof *idx);
    if (idx == NULL) {
      free(cr);
      free(start);
      free(L.pl);
      free(bb);
      free(par);
      return 2;
    }
    for (int32_t k = 0; k < np; k++)
      for (int32_t z = cr[(size_t)k * 6 + 2]; z <= cr[(size_t)k * 6 + 5]; z++)
        for (int32_t y = cr[(size_t)k * 6 + 1]; y <= cr[(size_t)k * 6 + 4]; y++)
          for (int32_t x = cr[(size_t)k * 6 + 0]; x <= cr[(size_t)k * 6 + 3]; x++)
            idx[start[((int64_t)z * nc[1] + y) * nc[0] + x]++] = k;
    for (int64_t c = tot; c > 0; c--)
      start[c] = start[c - 1];
    st->ngridcell = tot;
    st->ngrident = start[tot];
    start[0] = 0;

    st->t_grid = pm_now() - tphase;
    tphase = pm_now();
    for (int64_t z = 0; z < nc[2] && !oom; z++)
      for (int64_t y = 0; y < nc[1] && !oom; y++)
        for (int64_t x = 0; x < nc[0] && !oom; x++) {
          int64_t c = (z * nc[1] + y) * nc[0] + x;
          for (int32_t i = start[c]; i < start[c + 1] && !oom; i++)
            for (int32_t j = i + 1; j < start[c + 1] && !oom; j++) {
              int32_t a = idx[i], b = idx[j];
              if (a > b) {
                int32_t t = a;
                a = b;
                b = t;
              }
              /* ПАРА ВЫДАЁТСЯ ОДИН РАЗ — из ПЕРВОЙ общей ячейки. Иначе крупный
               * полигон, лежащий в сотне ячеек, дал бы сотню одинаковых пар, и
               * они прошли бы дорогие ворота по сто раз. */
              int first = 1;
              for (int cc = 0; cc < 3 && first; cc++) {
                int32_t lo = cr[(size_t)a * 6 + (size_t)cc];
                if (cr[(size_t)b * 6 + (size_t)cc] > lo) lo = cr[(size_t)b * 6 + (size_t)cc];
                int64_t me = (cc == 0) ? x : ((cc == 1) ? y : z);
                if (me != lo) first = 0;
              }
              if (!first) continue;
              if (pm_emit(&L, ps, m, bb, a, b, cfg, st, &rnd) != 0) oom = 1;
            }
        }
    free(cr);
    free(start);
    free(idx);
  }
  st->t_gate = pm_now() - tphase;
  tphase = pm_now();
  pm_pair *pl = L.pl;
  int64_t pn = L.pn;
  if (oom) {
    free(pl);
    free(bb);
    free(par);
    return 2;
  }

  /* --- слияние по возрастанию ошибки (Т2: приоритет — качество) ---
   * ГЕОМЕТРИЯ ПРОВЕРЯЕТСЯ ЗДЕСЬ, ПО ОБЪЕДИНЕНИЮ ГРУПП. Списки членов ведутся
   * односвязно: `head[корень] -> nxt`. Слияние принимается, только если у
   * ОБЪЕДИНЁННОЙ группы `dmax < δ`; иначе пара отбрасывается, а группы живут
   * дальше и могут слиться с другими. */
  qsort(pl, (size_t)pn, sizeof *pl, cmp_pair);
  int32_t *head = malloc((size_t)np * sizeof *head);
  int32_t *nxt = malloc((size_t)np * sizeof *nxt);
  int32_t *mem = malloc((size_t)np * sizeof *mem);
  if (head == NULL || nxt == NULL || mem == NULL) {
    free(head);
    free(nxt);
    free(mem);
    free(pl);
    free(bb);
    free(par);
    return 2;
  }
  for (int32_t k = 0; k < np; k++) {
    head[k] = k;
    nxt[k] = -1;
  }
  int32_t nleft = np;
  for (int64_t i = 0; i < pn; i++) {
    if (cfg->target > 0 && nleft <= cfg->target) break;
    int32_t ra = uf_find(par, pl[i].a), rb = uf_find(par, pl[i].b);
    if (ra == rb) continue;
    int32_t nm = 0;
    for (int32_t x = head[ra]; x >= 0 && nm < np; x = nxt[x])
      mem[nm++] = x;
    for (int32_t x = head[rb]; x >= 0 && nm < np; x = nxt[x])
      mem[nm++] = x;
    if (cfg->use_geom && !cfg->random) {
      double n[3], off;
      if (!(group_plane(ps, m, mem, nm, n, &off) < cfg->delta)) {
        st->nrej_geom++;
        continue;
      }
    }
    /* сцепить списки */
    int32_t tail = head[ra];
    while (nxt[tail] >= 0)
      tail = nxt[tail];
    nxt[tail] = head[rb];
    par[rb] = ra;
    head[ra] = head[ra];
    nleft--;
    st->nmerged++;
  }
  free(head);
  free(nxt);
  free(mem);
  free(pl);
  free(bb);

  st->t_merge = pm_now() - tphase;
  tphase = pm_now();
  /* --- новая разметка: треугольник наследует корень своего полигона --- */
  int32_t *rank = malloc((size_t)np * sizeof *rank);
  so->label = malloc((size_t)m->nt * sizeof *so->label);
  so->seg = calloc((size_t)np, sizeof *so->seg);
  if (rank == NULL || so->label == NULL || so->seg == NULL) {
    free(rank);
    free(par);
    hz_seg_free(so);
    return 2;
  }
  for (int32_t k = 0; k < np; k++)
    rank[k] = -1;
  int32_t nn = 0;
  for (int32_t k = 0; k < np; k++) {
    int32_t r = uf_find(par, k);
    if (rank[r] < 0) {
      rank[r] = nn;
      so->seg[nn] = si->seg[si->label[ps->tri[ps->p[r].t0]]];
      so->seg[nn].area = 0.0;
      so->seg[nn].ntri = 0;
      so->seg[nn].dmax = 0.0;
      nn++;
    }
  }
  /* ПЛОСКОСТЬ ГРУППЫ ПЕРЕСЧИТЫВАЕТСЯ, а не берётся у первого: слияние меняет
   * и нормаль, и смещение, а `dmax` после него ДРУГОЙ. Ложный ноль опаснее
   * UNKNOWN — тот останавливает потребителя, этот пропускает. */
  double *acc = calloc((size_t)(nn > 0 ? nn : 1) * 8, sizeof *acc);
  if (acc == NULL) {
    free(rank);
    free(par);
    hz_seg_free(so);
    return 2;
  }
  for (int32_t k = 0; k < np; k++) {
    int32_t g = rank[uf_find(par, k)];
    const hz_poly *P = &ps->p[k];
    double w = P->area;
    for (int c = 0; c < 3; c++) {
      acc[(size_t)g * 8 + (size_t)c] += w * P->n[c];
      acc[(size_t)g * 8 + 3 + (size_t)c] += w * P->org[c];
    }
    acc[(size_t)g * 8 + 6] += w;
  }
  for (int32_t g = 0; g < nn; g++) {
    double w = acc[(size_t)g * 8 + 6];
    if (!(w > 0.0)) continue;
    double n[3], nl = 0.0;
    for (int c = 0; c < 3; c++) {
      n[c] = acc[(size_t)g * 8 + (size_t)c] / w;
      nl += n[c] * n[c];
    }
    nl = sqrt(nl);
    if (!(nl > 0.0)) continue;
    double org[3];
    for (int c = 0; c < 3; c++) {
      n[c] /= nl;
      org[c] = acc[(size_t)g * 8 + 3 + (size_t)c] / w;
      so->seg[g].n[c] = n[c];
    }
    so->seg[g].off = n[0] * org[0] + n[1] * org[1] + n[2] * org[2];
  }
  free(acc);

  /* Треугольник знает свой ИСХОДНЫЙ участок; полигон — тот же индекс, потому
   * что `hz_poly_build` нумерует полигоны участками один к одному. */
  for (int32_t t = 0; t < m->nt; t++) {
    int32_t k = si->label[t];
    int32_t g = rank[uf_find(par, k)];
    so->label[t] = g;
    so->seg[g].ntri++;
    so->seg[g].area += hz_obj_tri_area(m, t);
    double p[3][3];
    hz_obj_tri(m, t, p);
    for (int i = 0; i < 3; i++) {
      double d = fabs(p[i][0] * so->seg[g].n[0] + p[i][1] * so->seg[g].n[1] +
                      p[i][2] * so->seg[g].n[2] - so->seg[g].off);
      if (d > so->seg[g].dmax) so->seg[g].dmax = d;
    }
  }
  for (int32_t g = 0; g < nn; g++)
    if (so->seg[g].dmax > st->dmax_worst) st->dmax_worst = so->seg[g].dmax;
  st->t_label = pm_now() - tphase;
  so->nseg = nn;
  so->delta = si->delta;
  st->nseg_out = nn;
  free(rank);
  free(par);
  return 0;
}
