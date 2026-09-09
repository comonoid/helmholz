/* PLAN_TRANSPORT.md T5а. Марш по лучу: накопление τ и попадание в границу. */

#include "transport/ray3.h"
#include <math.h>
#include <stdatomic.h>
#include <stdio.h>
#include <string.h>

/* §829: прибор ДОЛИ МАРША В КОНВЕЙЕРЕ — счётчик вызовов и Σ шагов по ячейкам.
 * Пассивные атомики с relaxed (вызовы возможны из omp-регионов, А1036); физика
 * их не читает. Долю ВРЕМЕНИ даёт perf по символу tr3_march (план §829). */
static _Atomic long long g_stat_march = 0;
static _Atomic long long g_stat_march_steps = 0;

void hz_ray3_stats_get(long long *calls, long long *steps) {
  *calls = atomic_load_explicit(&g_stat_march, memory_order_relaxed);
  *steps = atomic_load_explicit(&g_stat_march_steps, memory_order_relaxed);
}

void hz_ray3_stats_print(const char *tag) {
  printf("   §829 МАРШ (%s): вызовов %lld, шагов %lld\n", tag,
         atomic_load_explicit(&g_stat_march, memory_order_relaxed),
         atomic_load_explicit(&g_stat_march_steps, memory_order_relaxed));
}

static int tr3_march_body(const tr3_scene *sc, const double o[3], const double d[3], double tmax,
                          tr3_hit *h);

/* Пересечение луча со СЛОЕМ [lo, hi] по оси a. Возвращает 0, если пусто. */
static int slab(double p0, double dd, double lo, double hi, double *t0, double *t1) {
  if (!(fabs(dd) > 0.0)) return p0 >= lo && p0 <= hi;
  double a = (lo - p0) / dd, b = (hi - p0) / dd;
  if (a > b) {
    double s = a;
    a = b;
    b = s;
  }
  if (a > *t0) *t0 = a;
  if (b < *t1) *t1 = b;
  return *t0 <= *t1;
}

/* Уточнение попадания по ПРИМИТИВУ (К15): фасеты режут ЯЧЕЙКУ, а отвечать на луч
 * обязан примитив — иначе картинка получит огранку на геометрии, заданной точно.
 *
 * КОРЕНЬ ВЫБИРАЕТСЯ БЛИЖАЙШИЙ К ФАСЕТНОМУ ОТВЕТУ, а не «попавший в отрезок
 * ячейки», и это правка по результату прогона. Требование «корень внутри
 * [t, t_вых]» отвергало 9 попаданий из 169 при k = 1: фасет отстоит от примитива
 * на dmax, и при dmax ~ 0.29 против ячейки в 1 единицу корень законно уезжает в
 * соседнюю ячейку. Окно шириной dmax тоже не годится — на скользящем луче
 * смещение вдоль луча есть dmax/|n·d| и не ограничено. Ближайший корень
 * однозначен, потому что фасетное тело отличается от примитива не более чем на
 * dmax, а корней всего два и они разнесены на хорду. */
static int sphere_hit(const hz_surf *s, const double o[3], const double d[3], double t0, double t1,
                      double *t, double n[3]) {
  if (s->kind != HZ_SURF_SPHERE) return 0;
  double m[3];
  for (int k = 0; k < 3; k++)
    m[k] = o[k] - s->p[k];
  double r = s->p[3];
  double b = m[0] * d[0] + m[1] * d[1] + m[2] * d[2];
  double c = m[0] * m[0] + m[1] * m[1] + m[2] * m[2] - r * r;
  double disc = b * b - c;
  if (disc < 0.0) return 0;
  double sq = sqrt(disc), root[2] = {-b - sq, -b + sq};
  for (int i = 0; i < 2; i++) {
    if (root[i] < t0 || root[i] > t1) continue;
    double nn[3], nd = 0.0;
    for (int k = 0; k < 3; k++) {
      nn[k] = (o[k] + root[i] * d[k] - s->p[k]) / r;
      nd += nn[k] * d[k];
    }
    /* ПОПАДАНИЕ ЕСТЬ ВХОД, а не касание: требуется n·d < 0 СТРОГО. Отсюда же
     * теневой луч, выходящий РОВНО с поверхности, сам себя не затеняет — на
     * своём корне у него n·d > 0. ε не нужен, ровно как у полуплоскостей. */
    if (!(nd < 0.0)) continue;
    *t = root[i];
    for (int k = 0; k < 3; k++)
      n[k] = nn[k];
    return 1;
  }
  return 0;
}

int tr3_march(const tr3_scene *sc, const double o[3], const double d[3], double tmax, tr3_hit *h) {
  /* §829: счёт — одной точкой, в обёртке; сумма шагов читается из готового hit,
   * поэтому место возврата не важно. */
  atomic_fetch_add_explicit(&g_stat_march, 1, memory_order_relaxed);
  int rc = tr3_march_body(sc, o, d, tmax, h);
  atomic_fetch_add_explicit(&g_stat_march_steps, h->nsteps, memory_order_relaxed);
  return rc;
}

static int tr3_march_body(const tr3_scene *sc, const double o[3], const double d[3], double tmax,
                          tr3_hit *h) {
  memset(h, 0, sizeof *h);
  h->facet = -1;
  h->surf = -1;
  h->cell = -1;

  /* К28: длина `sigma` СВЕРЯЕТСЯ, а не подразумевается. Проверка стоит одно
   * сравнение на весь марш и ловит ровно тот случай, ради которого заведена, —
   * дерево доросло после выделения массива. Fail closed: марш не выполняется. */
  if (sc->sigma != NULL && sc->nsigma < sc->tree->n) return TR3_MARCH_ESIGMA;

  int32_t n = (int32_t)1 << sc->tree->log2size;
  /* Луч в ЕДИНИЦАХ дерева. Параметр t остаётся МИРОВОЙ длиной, поэтому длины
   * отрезков сразу годятся для τ и делить их обратно не надо. */
  double p0[3], du[3];
  for (int a = 0; a < 3; a++) {
    p0[a] = (o[a] - sc->fr.o[a]) / sc->fr.u[a];
    du[a] = d[a] / sc->fr.u[a];
  }

  double t0 = 0.0, t1 = tmax >= 0.0 ? tmax : 1e300;
  for (int a = 0; a < 3; a++)
    if (!slab(p0[a], du[a], 0.0, (double)n, &t0, &t1)) return 0; /* куб не задет */
  if (t0 < 0.0) t0 = 0.0;
  if (t1 <= t0) return 0;

  double t = t0;
  hz_hspace hs[HZ_P3_MAXH];
  int32_t hid[HZ_P3_MAXH];

  /* ТЕКУЩАЯ ЯЧЕЙКА ВЕДЁТСЯ ЦЕЛЫМИ ЧИСЛАМИ, а не пересчитывается из точки на
   * каждом шаге. Пересчёт был бы ровно тем местом, где заводится ε: на выходной
   * грани при d < 0 координата точки РОВНО целая, floor даёт ту же ячейку, и
   * обход зациклился бы. Здесь шаг делается по целым, и такого случая нет. */
  int c[3];
  for (int a = 0; a < 3; a++) {
    double x = p0[a] + t * du[a];
    int v = (int)floor(x);
    if (v < 0) v = 0;
    if (v >= (int)n) v = (int)n - 1;
    c[a] = v;
  }

  while (t < t1) {
    int blo[3], bsize = 1;
    int32_t ni = hz_oct_leaf_box(sc->tree, c[0], c[1], c[2], blo, &bsize);
    if (ni < 0) break;
    h->nsteps++;

    /* выход из коробки: параметр и ОСЬ */
    double texit = t1;
    int axis = -1;
    for (int a = 0; a < 3; a++) {
      if (!(fabs(du[a]) > 0.0)) continue;
      double face = du[a] > 0.0 ? (double)(blo[a] + bsize) : (double)blo[a];
      double te = (face - p0[a]) / du[a];
      if (te < texit) {
        texit = te;
        axis = a;
      }
    }
    if (texit > t1) texit = t1;

    /* --- граница материала внутри этой ячейки --- */
    const hz_cutrec *rec = sc->cm != NULL ? hz_cutmap_find(sc->cm, ni) : NULL;
    int nh = rec != NULL ? hz_cutmap_hspaces(sc->ft, sc->cm, rec, hs, hid, NULL, HZ_P3_MAXH) : -1;
    /* К15 ЦЕЛИКОМ: если веер ячейки ссылается на ПРИМИТИВ, то и решение
     * «попал / не попал», и точка, и нормаль берутся у примитива, а фасеты
     * работают только пространственным признаком «эта поверхность рядом».
     * Первая редакция уточняла примитивом лишь ТОЧКУ, оставив решение фасетам, —
     * и силуэт с кромкой тени уезжали на O(dmax), то есть на несколько пикселей
     * при k = 1. Это измерено, а не предположено. */
    /* ВСЕ РАЗЛИЧНЫЕ ПРИМИТИВЫ ВЕЕРА, А НЕ ПЕРВЫЙ (К27).
     *
     * Прежде бралcя первый же фасет с `surf ≥ 0`, и на этом перебор кончался.
     * Если в веере ячейки встречались ДВА разных тела — две сферы, сфера и
     * цилиндр, — второе не проверялось вовсе, и луч проходил сквозь него. Тест
     * T5а этот случай не просто не покрывал, а ОБХОДИЛ проверкой `both == 0`.
     *
     * Берётся БЛИЖАЙШЕЕ попадание, а не первое найденное: на стыке двух тел
     * порядок фасетов в веере произволен, и «первое» дало бы то дальнее тело,
     * то ближнее в зависимости от того, как лёг отбор. */
    int32_t plist[HZ_P3_MAXH];
    int np = 0;
    if (!sc->facet_only) /* facet_only: тело есть МНОГОГРАННИК, примитива нет */
      for (int j = 0; j < nh; j++) {
        int32_t fi = hid[j] >= 0 ? hid[j] : ~hid[j];
        int32_t sf = sc->ft->f[fi].surf;
        if (sf < 0) continue;
        int seen = 0;
        for (int q = 0; q < np && !seen; q++)
          if (plist[q] == sf) seen = 1;
        if (!seen) plist[np++] = sf;
      }
    if (np > 0 && sc->st != NULL) {
      double bt = 0.0, bn[3] = {0, 0, 0};
      int32_t bs = -1;
      for (int q = 0; q < np; q++) {
        if (plist[q] >= sc->st->n) continue;
        double tp, npv[3];
        if (!sphere_hit(&sc->st->s[plist[q]], o, d, t, texit, &tp, npv)) continue;
        if (bs < 0 || tp < bt) {
          bt = tp;
          bs = plist[q];
          for (int a = 0; a < 3; a++)
            bn[a] = npv[a];
        }
      }
      if (bs >= 0) {
        h->hit = 1;
        h->refined = 1;
        h->t = bt;
        h->cell = ni;
        h->surf = bs;
        for (int a = 0; a < 3; a++) {
          h->n[a] = bn[a];
          h->p[a] = o[a] + bt * d[a];
        }
        h->tau += (sc->sigma != NULL ? sc->sigma[ni] : 0.0) * (bt - t);
        return 0;
      }
      nh = 0; /* ни один примитив эту ячейку не задел — поверхности здесь нет */
    }
    if (nh > 0) {
      double lo_t = t, hi_t = texit;
      int32_t enter_facet = -1;
      int empty = 0;
      for (int j = 0; j < nh && !empty; j++) {
        double s0 = -hs[j].off, sd = 0.0;
        for (int a = 0; a < 3; a++) {
          s0 += hs[j].n[a] * (p0[a] + t * du[a]);
          sd += hs[j].n[a] * du[a];
        }
        if (!(fabs(sd) > 0.0)) {
          if (s0 > 0.0) empty = 1; /* весь отрезок снаружи этой полуплоскости */
          continue;
        }
        double te = t - s0 / sd;
        if (sd < 0.0) { /* полупространство ВХОДИТСЯ */
          if (te > lo_t) {
            lo_t = te;
            enter_facet = hid[j] >= 0 ? hid[j] : ~hid[j];
          }
        } else { /* ПОКИДАЕТСЯ: сюда же попадает старт ровно на поверхности */
          if (te < hi_t) hi_t = te;
        }
      }
      if (!empty && lo_t <= hi_t && lo_t > t && enter_facet >= 0) {
        h->hit = 1;
        h->t = lo_t;
        h->cell = ni;
        h->facet = enter_facet;
        const hz_facet *f = &sc->ft->f[enter_facet];
        h->surf = f->surf;
        /* нормаль фасета — В ЕДИНИЦАХ; в мир она переносится делением на u и
         * пере-нормировкой (аффинное сжатие углов не сохраняет, PLAN_CUT Г21) */
        double nw[3], nm = 0.0;
        for (int a = 0; a < 3; a++) {
          nw[a] = f->n[a] / sc->fr.u[a];
          nm += nw[a] * nw[a];
        }
        nm = sqrt(nm);
        for (int a = 0; a < 3; a++)
          h->n[a] = nw[a] / nm;
        /* Уточнять примитивом здесь нечего: ячейки со ссылкой на примитив
         * разобраны ВЫШЕ и сюда не доходят. Сюда попадают только плоскости,
         * заданные прямо, и фасеты dual contouring — у них примитива нет по
         * определению (surf < 0), и фасет ЕСТЬ поверхность, а не приближение. */
        for (int a = 0; a < 3; a++)
          h->p[a] = o[a] + h->t * d[a];
        h->tau += (sc->sigma != NULL ? sc->sigma[ni] : 0.0) * (h->t - t);
        return 0;
      }
    }

    /* --- среда --- Та же замкнутая форма, что в T1: трапеция для линейной σ
     * точна, а при постоянной по ячейке σ она вырождается в σ·длину. Это
     * частный случай уже проверенного, а не новая формула. */
    if (sc->sigma != NULL) h->tau += sc->sigma[ni] * (texit - t);

    if (axis < 0) break; /* луч параллелен всем осям — невозможно при |d| = 1 */
    /* ЦЕЛОЧИСЛЕННЫЙ ШАГ ЧЕРЕЗ ГРАНЬ, без ε. По оси выхода координата шагает
     * через грань; по остальным берётся из точки выхода с ЗАЖИМОМ в текущую
     * коробку — тогда луч, прошедший ровно через ребро, входит в соседа через
     * шаг, а не теряет ячейку и не создаёт щели. */
    for (int a = 0; a < 3; a++) {
      if (a == axis) continue;
      double x = p0[a] + texit * du[a];
      int v = (int)floor(x);
      if (v < blo[a]) v = blo[a];
      if (v > blo[a] + bsize - 1) v = blo[a] + bsize - 1;
      c[a] = v;
    }
    c[axis] = du[axis] > 0.0 ? blo[axis] + bsize : blo[axis] - 1;
    if (c[axis] < 0 || c[axis] >= (int)n) break;
    t = texit;
    if (t >= t1) break;
  }
  return 0;
}
