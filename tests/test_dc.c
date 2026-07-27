/* Фальсификаторы эрмитова дерева и QEF-иерархии (PLAN_CUT.md, Р-5б; Г13, Г14,
 * Г41, Г47, Г49). ПРЕДСКАЗАНИЯ ПЕЧАТАЮТСЯ ДО РЕЗУЛЬТАТОВ. */

#include "cut/dc.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fail = 0, g_total = 0;

static void check(int ok, const char *what) {
  g_total++;
  if (!ok) {
    g_fail++;
    printf("FAIL: %s\n", what);
  }
}

/* --------------------------------------------------------- источники данных */

typedef struct {
  double c[3], r;
} sph;

static int sph_sign(void *ctx, const int32_t p[3]) {
  const sph *s = ctx;
  double d = 0.0;
  for (int k = 0; k < 3; k++) {
    double q = (double)p[k] - s->c[k];
    d += q * q;
  }
  return d <= s->r * s->r;
}

static int sph_cross(void *ctx, const int32_t p[3], int axis, double *t, double nrm[3]) {
  const sph *s = ctx;
  double d[3];
  for (int k = 0; k < 3; k++)
    d[k] = (double)p[k] - s->c[k];
  double b = d[axis], c = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - s->r * s->r;
  double disc = b * b - c;
  if (disc < 0.0) return 1;
  double sq = sqrt(disc), s0 = -b - sq, s1 = -b + sq, ss = -1.0;
  if (s0 >= 0.0 && s0 <= 1.0)
    ss = s0;
  else if (s1 >= 0.0 && s1 <= 1.0)
    ss = s1;
  if (ss < 0.0) return 1;
  *t = ss;
  d[axis] += ss;
  for (int k = 0; k < 3; k++)
    nrm[k] = d[k] / s->r;
  return 0;
}

/* ТО ЖЕ пересечение, посчитанное ОТ ВЕРХНЕГО конца ребра: сегмент тот же,
 * параметризация обратная. Алгебраически одно и то же, в битах — нет. Так
 * выглядела бы независимая подгонка на ячейку (Г3, Г39). */
static void sph_cross_rev(const sph *s, const int32_t p[3], int axis, double out[3]) {
  double q[3] = {(double)p[0], (double)p[1], (double)p[2]};
  q[axis] += 1.0;
  double d[3];
  for (int k = 0; k < 3; k++)
    d[k] = q[k] - s->c[k];
  double b = -d[axis], c = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - s->r * s->r;
  double disc = b * b - c, sq = sqrt(disc);
  double s0 = -b - sq, s1 = -b + sq, ss = (s0 >= 0.0 && s0 <= 1.0) ? s0 : s1;
  for (int k = 0; k < 3; k++)
    out[k] = q[k];
  out[axis] -= ss;
}

/* Выпуклое тело — пересечение полупространств {n_i·x <= d_i}. Пересечение ребра
 * считается ТОЧНО (вдоль ребра всё линейно), нормаль берётся у активной в этой
 * точке плоскости: получается НАСТОЯЩЕЕ острое ребро, а не сглаженное. */
typedef struct {
  int np;
  double n[4][3], d[4];
} conv;

static double conv_f(const conv *v, const double x[3]) {
  double f = -1e300;
  for (int i = 0; i < v->np; i++) {
    double g = v->n[i][0] * x[0] + v->n[i][1] * x[1] + v->n[i][2] * x[2] - v->d[i];
    if (g > f) f = g;
  }
  return f;
}

static int conv_sign(void *ctx, const int32_t p[3]) {
  const conv *v = ctx;
  double x[3] = {(double)p[0], (double)p[1], (double)p[2]};
  return conv_f(v, x) <= 0.0;
}

static int conv_cross(void *ctx, const int32_t p[3], int axis, double *t, double nrm[3]) {
  const conv *v = ctx;
  double x0[3] = {(double)p[0], (double)p[1], (double)p[2]};
  double best = -1.0;
  int bi = -1;
  for (int i = 0; i < v->np; i++) {
    double g0 = v->n[i][0] * x0[0] + v->n[i][1] * x0[1] + v->n[i][2] * x0[2] - v->d[i];
    double g1 = g0 + v->n[i][axis];
    if (!(fabs(g1 - g0) > 0.0)) continue;
    double s = -g0 / (g1 - g0);
    if (s < 0.0 || s > 1.0) continue;
    double y[3] = {x0[0], x0[1], x0[2]};
    y[axis] += s;
    if (fabs(conv_f(v, y)) > 1e-12) continue; /* корень другой плоскости, не границы */
    if (s > best) {
      best = s;
      bi = i;
    }
  }
  if (bi < 0) return 1;
  *t = best;
  for (int k = 0; k < 3; k++)
    nrm[k] = v->n[bi][k];
  return 0;
}

/* НЕГАТИВНЫЙ КОНТРОЛЬ: ЗАНЯТОСТЬ. Известно только, что ребро пересечено; точка
 * берётся посередине, нормаль — с сетки. Ровно то представление, из которого
 * получается лесенка (CLAUDE.md, Г40). */
static int occ_cross(void *ctx, const int32_t p[3], int axis, double *t, double nrm[3]) {
  const conv *v = ctx;
  double x0[3] = {(double)p[0], (double)p[1], (double)p[2]};
  *t = 0.5;
  for (int k = 0; k < 3; k++)
    nrm[k] = 0.0;
  nrm[axis] = conv_f(v, x0) <= 0.0 ? 1.0 : -1.0;
  return 0;
}

/* ------------------------------------------------------- обход узлов с геометрией */

typedef void (*nodecb)(void *ctx, const hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size,
                       int lv);

static void for_nodes(const hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size, int lv,
                      nodecb cb, void *ctx) {
  cb(ctx, t, ni, lo, size, lv);
  int32_t c0 = t->nd[ni].child0;
  if (c0 < 0) return;
  int32_t half = size / 2;
  for (int i = 0; i < 8; i++) {
    int32_t clo[3];
    for (int a = 0; a < 3; a++)
      clo[a] = lo[a] + (((i >> a) & 1) ? half : 0);
    for_nodes(t, c0 + i, clo, half, lv + 1, cb, ctx);
  }
}

static void walk_tree(const hz_dctree *t, nodecb cb, void *ctx) {
  int32_t zero[3] = {0, 0, 0};
  for_nodes(t, 0, zero, (int32_t)1 << t->log2size, 0, cb, ctx);
}

/* --------------------------------------------------------------- 1. сфера */

typedef struct {
  const sph *s;
  int finest;
  double dworst;    /* max |расстояние вершины до сферы| на самом мелком уровне */
  double eworst[8]; /* max невязка по уровню */
  int cnt[8], nleaf, nvert;
  double dconv; /* то же, но при общепринятом пороге усечения */
  int nclamp_conv;
} sphstat;

static void sph_cb(void *ctx, const hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size,
                   int lv) {
  sphstat *st = ctx;
  const hz_dcnode *nd = &t->nd[ni];
  if (!(nd->flags & HZ_DC_HASVERT)) return;
  st->cnt[lv]++;
  st->nvert++;
  if (nd->err > st->eworst[lv]) st->eworst[lv] = nd->err;
  double d = 0.0;
  for (int k = 0; k < 3; k++) {
    double q = nd->vx[k] - st->s->c[k];
    d += q * q;
  }
  double e = fabs(sqrt(d) - st->s->r);
  if (size == 1) {
    st->nleaf++;
    if (e > st->dworst) st->dworst = e;
    /* ТА ЖЕ форма при ОБЩЕПРИНЯТОМ пороге 0.01 по собственным — чтобы сравнить
     * выведенный порог не с воздухом, а с тем, что делают все */
    double x[3], r;
    hz_qef_solve_tau(&nd->q, 0.01, x, &r);
    int cl = 0;
    for (int k = 0; k < 3; k++) {
      if (x[k] < 0.0) {
        x[k] = 0.0;
        cl = 1;
      } else if (x[k] > (double)size) {
        x[k] = (double)size;
        cl = 1;
      }
    }
    st->nclamp_conv += cl;
    double dd = 0.0;
    for (int k = 0; k < 3; k++) {
      double q = (double)lo[k] + x[k] - st->s->c[k];
      dd += q * q;
    }
    double ee = fabs(sqrt(dd) - st->s->r);
    if (ee > st->dconv) st->dconv = ee;
  }
}

static void sphere_at(int L, sphstat *st, sph *s, int *nclamped, int *nnodes) {
  hz_signgrid g = {NULL, 0, 0};
  hz_htab ht;
  hz_htab_init(&ht);
  check(hz_dc_sample(&g, &ht, L, sph_sign, sph_cross, s) == HZ_DC_OK, "опрос источника");
  hz_dctree t;
  hz_dc_init(&t, L);
  check(hz_dc_build(&t, &g, &ht) == HZ_DC_OK, "сборка дерева");
  memset(st, 0, sizeof *st);
  st->s = s;
  st->finest = L;
  walk_tree(&t, sph_cb, st);
  *nclamped = t.nclamped;
  *nnodes = t.n;
  hz_dc_free(&t);
  hz_htab_free(&ht);
  hz_signgrid_free(&g);
}

static void t_sphere(void) {
  double dprev = 0.0;
  int mono_err = 1;
  for (int L = 4; L <= 5; L++) {
    sph s = {{(double)(1 << L) / 2 + 0.3, (double)(1 << L) / 2 - 0.3, (double)(1 << L) / 2 + 0.1},
             (double)(1 << L) * 0.36};
    sphstat st;
    int ncl = 0, nn = 0;
    sphere_at(L, &st, &s, &ncl, &nn);
    printf("  [сфера L=%d, r=%.1f] узлов %d, вершин %d, листьев %d; max |d − r| на листе %.4f; "
           "зажатых %d (%.1f%%)\n",
           L, s.r, nn, st.nvert, st.nleaf, st.dworst, ncl, 100.0 * (double)ncl / (double)st.nleaf);
    printf("    невязка по уровням: ");
    for (int i = 0; i <= L; i++)
      printf("%d:%.2e ", i, st.eworst[i]);
    printf("\n");
    printf("    ТОТ ЖЕ набор при общепринятом пороге 0.01: max |d − r| %.4f, зажатых %d\n",
           st.dconv, st.nclamp_conv);
    for (int i = 1; i <= L; i++)
      if (st.eworst[i] > st.eworst[i - 1]) mono_err = 0;
    printf("    сагитта ячейки h²/(8r) = %.4f при h = 1\n", 1.0 / (8.0 * s.r));
    if (L == 5) {
      /* ЗАКОН ЗДЕСЬ h²/r, А НЕ h². Ячейка всегда единичная — мельче минимальной
       * не бывает, — поэтому «удвоить глубину» значит удвоить радиус в единицах,
       * то есть уполовинить h/r. Предсказание «×4» было записано по невнимании:
       * при фиксированном h и удвоенном r сагитта падает ровно ВДВОЕ. */
      double ratio = dprev / st.dworst;
      printf("    сходимость: %.4f -> %.4f (отношение %.2f; закон h²/(8r) даёт 2)\n", dprev,
             st.dworst, ratio);
      check(ratio > 1.6 && ratio < 2.6, "ошибка вершины падает как h²/r — вдвое на удвоение r");
      check(st.dworst < 4.0 / (8.0 * s.r),
            "и держится в пределах нескольких сагитт (max берётся по ХУДШЕЙ ячейке)");
    }
    dprev = st.dworst;
  }
  check(mono_err, "невязка падает с уровнем — значит по ней МОЖНО выбирать уровень");
}

/* ------------------------------------------- 2. острое ребро и занятость */

typedef struct {
  double ex, ey;
  double worst;
  int n;
} edgestat;

static void edge_cb(void *ctx, const hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size,
                    int lv) {
  edgestat *es = ctx;
  const hz_dcnode *nd = &t->nd[ni];
  (void)lv;
  if (size != 1 || nd->child0 >= 0 || !(nd->flags & HZ_DC_HASVERT)) return;
  /* только ячейки, СКВОЗЬ которые проходит линия ребра */
  if (es->ex < (double)lo[0] || es->ex > (double)lo[0] + 1.0) return;
  if (es->ey < (double)lo[1] || es->ey > (double)lo[1] + 1.0) return;
  double d = hypot(nd->vx[0] - es->ex, nd->vx[1] - es->ey);
  if (d > es->worst) es->worst = d;
  es->n++;
}

static void t_sharp(void) {
  const int L = 4;
  /* Клин с ребром вдоль z: два полупространства под 90° плюс крышки по z.
   *
   * РЕБРО ПОСТАВЛЕНО ВНУТРЬ ЯЧЕЙКИ НАМЕРЕННО, и первая попытка этого не сделала.
   * Там вершина клина приходилась на ячейку, ВСЕ ВОСЕМЬ углов которой лежали
   * снаружи: материал заходил в неё тонким носом, ни одного узла сетки не задев.
   * Сетка знаков такого не видит вовсе — и это не дефект DC, а свойство
   * ОБЪЁМНОГО ВХОДА (Г14: подпиксельное тело в знаковое поле не влезает).
   * Мерить на такой конфигурации остроту ребра нельзя: мерить было бы нечего. */
  const double ax = 8.63, ay = 8.39, s2 = 0.70710678118654752;
  conv v = {4,
            {{s2, s2, 0.0}, {-s2, s2, 0.0}, {0.0, 0.0, 1.0}, {0.0, 0.0, -1.0}},
            {s2 * (ax + ay), s2 * (ay - ax), 14.0, -2.0}};
  double ex = ax, ey = ay;
  for (int mode = 0; mode < 2; mode++) {
    hz_signgrid g = {NULL, 0, 0};
    hz_htab ht;
    hz_htab_init(&ht);
    check(hz_dc_sample(&g, &ht, L, conv_sign, mode ? occ_cross : conv_cross, &v) == HZ_DC_OK,
          "опрос источника (клин)");
    hz_dctree t;
    hz_dc_init(&t, L);
    hz_dc_build(&t, &g, &ht);
    edgestat es = {ex, ey, 0.0, 0};
    walk_tree(&t, edge_cb, &es);
    printf("  [%s] ячеек НА ребре %d, max расстояние вершины от линии %.3e\n",
           mode ? "НК занятость" : "эрмитовы данные", es.n, es.worst);
    if (mode == 0) {
      check(es.n >= 4, "у ребра есть что мерить");
      check(es.worst < 1e-9, "эрмитовы нормали: вершины ЛЕЖАТ на ребре");
    } else {
      check(es.worst > 0.05, "негативный контроль: занятость ребро СКРУГЛЯЕТ");
    }
    hz_dc_free(&t);
    hz_htab_free(&ht);
    hz_signgrid_free(&g);
  }
}

/* ------------------------------------- 3. общее ребро: одна запись на четверых */

static void t_shared_edge(void) {
  sph s = {{16.3, 15.7, 16.1}, 11.5};
  hz_signgrid g = {NULL, 0, 0};
  hz_htab ht;
  hz_htab_init(&ht);
  check(hz_dc_sample(&g, &ht, 5, sph_sign, sph_cross, &s) == HZ_DC_OK, "опрос");
  check(hz_htab_sort(&ht) == HZ_DC_OK, "в таблице НЕТ повторов: ребро занесено ОДИН раз");

  /* НЕГАТИВНЫЙ КОНТРОЛЬ (предсказан провал): та же точка, посчитанная НЕЗАВИСИМО
   * от другого конца ребра. Алгебраически одно и то же; в битах — нет, и вот
   * ровно поэтому ребро есть ОБЩИЙ ОБЪЕКТ, а не формула, повторённая в каждой
   * из четырёх ячеек (Г3 на уровне данных, Г39 на уровне ядра). */
  int diff = 0;
  double maxd = 0.0;
  for (int32_t i = 0; i < ht.n; i++) {
    const hz_hedge *x = &ht.e[i];
    double a[3], b[3];
    hz_hedge_point(x, a);
    sph_cross_rev(&s, x->p, x->axis, b);
    if (memcmp(a, b, sizeof a) != 0) {
      diff++;
      double d = fabs(a[x->axis] - b[x->axis]);
      if (d > maxd) maxd = d;
    }
  }
  printf("  [НК независимый счёт] разошлись %d из %d записей, max на %.2e\n", diff, ht.n, maxd);
  check(diff > 0, "негативный контроль: независимый счёт того же пересечения ломает биты");
  check(maxd < 1e-12, "и расхождение именно в последних битах, а не в геометрии");

  hz_htab_free(&ht);
  hz_signgrid_free(&g);
}

/* ------------------------------- 4. Г13/Г49: починка против перестройки */

static int nodes_equal(const hz_dctree *a, const hz_dctree *b, int32_t *first_diff) {
  if (a->n != b->n) return 0;
  for (int32_t i = 0; i < a->n; i++)
    if (memcmp(&a->nd[i], &b->nd[i], sizeof(hz_dcnode)) != 0) {
      if (first_diff != NULL) *first_diff = i;
      return 0;
    }
  return 1;
}

static void t_repair(void) {
  const int L = 4;
  sph s = {{8.3, 7.7, 8.1}, 5.5};
  hz_signgrid g = {NULL, 0, 0};
  hz_htab ht;
  hz_htab_init(&ht);
  check(hz_dc_sample(&g, &ht, L, sph_sign, sph_cross, &s) == HZ_DC_OK, "опрос");

  hz_dctree t0;
  hz_dc_init(&t0, L);
  hz_dc_build(&t0, &g, &ht);

  /* «скол»: правится ОДНА эрмитова запись — точка съехала, нормаль повернулась.
   * Знаки не трогаются, поэтому форма дерева та же (иначе это перестройка). */
  int32_t k = ht.n / 3;
  hz_hedge *e = &ht.e[k];
  e->t = e->t * 0.5 + 0.25;
  double nn[3] = {e->nrm[0] + 0.3, e->nrm[1] - 0.2, e->nrm[2] + 0.1};
  double m = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
  for (int a = 0; a < 3; a++)
    e->nrm[a] = nn[a] / m;

  hz_dctree tr;
  hz_dc_init(&tr, L);
  hz_dc_build(&tr, &g, &ht); /* ПОЛНАЯ ПЕРЕСТРОЙКА с правленой таблицей */

  /* ПОЧИНКА: ребро принадлежит ЧЕТЫРЁМ ячейкам, чинить надо все четыре. */
  int u = (e->axis + 1) % 3, v = (e->axis + 2) % 3;
  int nrep = 0;
  for (int du = 0; du <= 1; du++)
    for (int dv = 0; dv <= 1; dv++) {
      int32_t cell[3] = {e->p[0], e->p[1], e->p[2]};
      cell[u] -= du;
      cell[v] -= dv;
      if (cell[u] < 0 || cell[v] < 0) continue;
      int rc = hz_dc_repair(&t0, &ht, cell);
      check(rc == HZ_DC_OK || rc == HZ_DC_ETOPO, "починка вернула определённый код");
      if (rc == HZ_DC_OK) nrep++;
    }

  int32_t fd = -1;
  int same = nodes_equal(&t0, &tr, &fd);
  printf("  [Г13] узлов %d, починено ячеек %d; починка против перестройки побитово=%d "
         "(первое расхождение %d)\n",
         tr.n, nrep, same, fd);
  check(same, "Г49: ПОЧИНКА даёт ПОБИТОВО то же, что полная перестройка");

  /* НЕГАТИВНЫЙ КОНТРОЛЬ (предсказан провал): дерево, собранное по ИСХОДНОЙ
   * таблице, с правкой только в листе и без пересчёта предков. */
  hz_dctree t1;
  hz_dc_init(&t1, L);
  hz_htab ht0;
  hz_htab_init(&ht0);
  hz_signgrid g0 = {NULL, 0, 0};
  hz_dc_sample(&g0, &ht0, L, sph_sign, sph_cross, &s);
  hz_dc_build(&t1, &g0, &ht0);
  {
    int32_t cell[3] = {e->p[0], e->p[1], e->p[2]};
    int32_t ni = 0, size = (int32_t)1 << L, lo[3] = {0, 0, 0};
    while (size > 1 && t1.nd[ni].child0 >= 0) {
      int32_t half = size / 2;
      int bit = 0;
      for (int a = 0; a < 3; a++)
        if (cell[a] >= lo[a] + half) {
          bit |= 1 << a;
          lo[a] += half;
        }
      ni = t1.nd[ni].child0 + bit;
      size = half;
    }
    if (size == 1) hz_dc_repair(&t1, &ht0, cell); /* правим лист по СТАРЫМ данным */
  }
  int32_t fd2 = -1;
  int same2 = nodes_equal(&t1, &tr, &fd2);
  printf("  [НК старые данные] побитово=%d (первое расхождение %d)\n", same2, fd2);
  check(!same2, "негативный контроль: без правки данных совпадения быть НЕ МОЖЕТ");

  hz_dc_free(&t0);
  hz_dc_free(&tr);
  hz_dc_free(&t1);
  hz_htab_free(&ht);
  hz_htab_free(&ht0);
  hz_signgrid_free(&g);
  hz_signgrid_free(&g0);
}

/* ------------------------------------------------ 5. Г47: неманифолдность */

static int two_balls_sign(void *ctx, const int32_t p[3]) {
  const double *r = ctx;
  double n = (double)((int32_t)1 << 4);
  double d0 = 0.0, d1 = 0.0;
  for (int k = 0; k < 3; k++) {
    d0 += (double)p[k] * (double)p[k];
    double q = (double)p[k] - n;
    d1 += q * q;
  }
  return d0 <= *r * *r || d1 <= *r * *r;
}

static int two_balls_cross(void *ctx, const int32_t p[3], int axis, double *t, double nrm[3]) {
  const double *r = ctx;
  double n = (double)((int32_t)1 << 4);
  for (int b = 0; b < 2; b++) {
    sph s;
    for (int k = 0; k < 3; k++)
      s.c[k] = b ? n : 0.0;
    s.r = *r;
    if (sph_cross(&s, p, axis, t, nrm) == 0) return 0;
  }
  return 1;
}

static void t_manifold(void) {
  check(hz_dc_manifold(0x00u) == 1, "маска 0x00: поверхности нет");
  check(hz_dc_manifold(0xFFu) == 1, "маска 0xFF: поверхности нет");
  check(hz_dc_manifold(0x01u) == 1, "маска 0x01: один угол внутри — манифолдна");
  check(hz_dc_manifold(0x0Fu) == 1, "маска 0x0F: половина куба — манифолдна");
  check(hz_dc_manifold(0x81u) == 0, "маска 0x81: ДВА диагональных угла — два листа (Г47)");
  check(hz_dc_manifold(0x18u) == 0, "маска 0x18: то же на другой диагонали");

  /* Два шара в противоположных углах: у КРУПНОГО узла маска диагональная, и
   * одна вершина его не представляет. Отказ обязан быть виден снаружи. */
  double r = 4.5;
  hz_signgrid g = {NULL, 0, 0};
  hz_htab ht;
  hz_htab_init(&ht);
  check(hz_dc_sample(&g, &ht, 4, two_balls_sign, two_balls_cross, &r) == HZ_DC_OK, "опрос");
  hz_dctree t;
  hz_dc_init(&t, 4);
  int rc = hz_dc_build(&t, &g, &ht);
  printf("  [Г47] два шара по диагонали: код сборки %d, неманифолдных узлов %d, "
         "у корня вершина=%d\n",
         rc, t.nmulti, (t.nd[0].flags & HZ_DC_HASVERT) ? 1 : 0);
  check(rc == HZ_DC_EMULTI, "сборка СООБЩАЕТ о неманифолдной конфигурации, а не молчит");
  check(t.nmulti > 0, "и считает такие узлы");
  check(!(t.nd[0].flags & HZ_DC_HASVERT),
        "у корня вершины НЕТ — fail closed, а не одна на два листа");
  hz_dc_free(&t);
  hz_htab_free(&ht);
  hz_signgrid_free(&g);
}

/* ------------------------------------------- 6. Г14: тонкое тело — адрес Г12 */

typedef struct {
  int nleaf;
} leafcount;

static void leaf_cb(void *ctx, const hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size,
                    int lv) {
  leafcount *lc = ctx;
  (void)lo;
  (void)lv;
  if (size == 1 && (t->nd[ni].flags & HZ_DC_HASVERT)) lc->nleaf++;
}

static void t_thin(void) {
  const int L = 5;
  /* клин с очень острым носом: полураствор ~0.05 рад */
  conv v = {4,
            {{-0.0499, 0.9988, 0.0}, {-0.0499, -0.9988, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}},
            {-0.0499 * 4.0 + 0.9988 * 16.0, -0.0499 * 4.0 - 0.9988 * 16.0, 28.0, 28.0}};
  hz_signgrid g = {NULL, 0, 0};
  hz_htab ht;
  hz_htab_init(&ht);
  check(hz_dc_sample(&g, &ht, L, conv_sign, conv_cross, &v) == HZ_DC_OK, "опрос (тонкий клин)");
  hz_dctree t;
  hz_dc_init(&t, L);
  int rc = hz_dc_build(&t, &g, &ht);
  leafcount lc = {0};
  walk_tree(&t, leaf_cb, &lc);
  printf("  [Г14] тонкий клин: код %d, листьев %d, зажатых %d (%.1f%%), неманифолдных %d\n", rc,
         lc.nleaf, t.nclamped, 100.0 * (double)t.nclamped / (double)lc.nleaf, t.nmulti);
  check(t.nclamped > 0, "Г12/Г14: на ТОНКОМ теле вершина уезжает из узла и её зажимает");
  hz_dc_free(&t);
  hz_htab_free(&ht);
  hz_signgrid_free(&g);
}

/* -------------------------------------------------------------------------- */

int main(void) {
  printf("=== ПРЕДСКАЗАНИЯ (до единого результата) ===\n");
  printf("  1 сфера: ошибка вершины падает как h² при удвоении глубины, а НЕВЯЗКА\n");
  printf("    падает с уровнем — значит по ней можно выбирать уровень (Г11)\n");
  printf("  2 острое ребро: вершины ЛЕЖАТ на линии пересечения (<=1e-9)\n");
  printf("  3 эрмитово ребро в таблице ОДНО (повторов нет)\n");
  printf("  4 Г13/Г49: ПОЧИНКА даёт ПОБИТОВО то же, что полная перестройка\n");
  printf("  5 Г47: диагональная маска -> отказ и счёт, а не вершина на два листа\n");
  printf("  6 Г14: на тонком теле зажим срабатывает — это и есть адрес оговорок Г12\n");
  printf("=== НЕГАТИВНЫЕ КОНТРОЛИ, каждый с ПРЕДСКАЗАННЫМ провалом ===\n");
  printf("  занятость вместо нормалей -> ребро ОБЯЗАНО скруглиться (>0.05)\n");
  printf("  независимый счёт того же пересечения -> биты ОБЯЗАНЫ разойтись\n");
  printf("  починка по СТАРЫМ данным -> совпадения с перестройкой БЫТЬ НЕ МОЖЕТ\n");
  printf("=== РЕЗУЛЬТАТЫ ===\n");

  t_sphere();
  t_sharp();
  t_shared_edge();
  t_repair();
  t_manifold();
  t_thin();

  printf("%s: %d/%d\n", g_fail ? "FAILURES" : "ok", g_total - g_fail, g_total);
  return g_fail != 0;
}
