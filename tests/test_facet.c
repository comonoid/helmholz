/* Фальсификаторы фасетизации примитива (PLAN_CUT.md, Р-5а; Г35, Г37, Г38).
 * ПРЕДСКАЗАНИЯ ПЕЧАТАЮТСЯ ДО РЕЗУЛЬТАТОВ.
 *
 * Метрика взята НЕ на полной сфере (Г35): подгон, который масштабирует все
 * смещения на общий множитель, делает «объём многогранника = объёму сферы»
 * тождеством при любом, сколь угодно грубом k. Мерится объём ЯЧЕЙКИ, отсечённой
 * сферой, против аналитического шарового сегмента — эту величину общий
 * множитель не спасает. */
#include "cut/surf.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

static int g_fail = 0, g_total = 0;

static void check(int ok, const char *what) {
  g_total++;
  if (!ok) {
    g_fail++;
    printf("FAIL: %s\n", what);
  }
}

static const hz_frame FR = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};

/* Шаровой сегмент {|x| <= r} ^ {x0 >= a}, a < r: V = pi (r-a)^2 (2r+a) / 3 */
static double cap_volume(double r, double a) {
  return M_PI * (r - a) * (r - a) * (2.0 * r + a) / 3.0;
}

/* Объём тела в коробке по ОТОБРАННЫМ фасетам. sel<0 — коробка вне тела. */
static double body_in_box(const hz_facettab *ft, int32_t f0, int32_t nf, const int32_t lo[3],
                          const int32_t hi[3], int *status) {
  int32_t sel[HZ_P3_MAXH];
  int ns = hz_facets_for_box(ft, f0, nf, lo, hi, sel, HZ_P3_MAXH);
  *status = ns;
  if (ns < 0) return 0.0;
  /* ноль-инициализация: при ns == 0 (ячейка целиком внутри тела, Г38) цикл ниже
   * не выполняется, и cppcheck справедливо не может доказать, что ядро туда не
   * заглянет */
  hz_hspace h[HZ_P3_MAXH] = {{{0, 0, 0}, 0}};
  int32_t hid[HZ_P3_MAXH] = {0};
  for (int j = 0; j < ns; j++) {
    memcpy(h[j].n, ft->f[sel[j]].n, sizeof h[j].n);
    h[j].off = ft->f[sel[j]].off;
    hid[j] = sel[j];
  }
  hz_poly3 p;
  if (hz_poly3_cut(&p, lo, hi, h, hid, ns) != HZ_P3_OK) return 0.0;
  return hz_poly3_volume(&p, &FR);
}

/* ------------------------------- 1. сходимость и Г35 (метрика на сегменте) */

static void t_convergence(void) {
  const double r = 3.0, a = 2.0;
  int32_t lo[3] = {2, -4, -4}, hi[3] = {4, 4, 4}; /* содержит весь сегмент x>=2 */
  double c[3] = {0, 0, 0};
  double exact = cap_volume(r, a);

  printf("  [Р-5а] сегмент r=%.1f a=%.1f, точный объём %.6f\n", r, a, exact);
  double prev_ins = 0.0, prev_fit = 0.0;
  int mono_ins = 1, mono_fit = 1, fit_better = 1, ins_one_signed = 1, kmax = -1;
  for (int k = 0; k <= HZ_ICO_MAXK; k++) {
    int toomany = 0;
    double e[2];
    int nfac[2];
    for (int fit = 0; fit <= 1; fit++) {
      hz_facettab ft;
      check(hz_facettab_init(&ft) == 0, k == 0 && fit == 0 ? "init" : "");
      int32_t f0 = 0;
      int32_t n = hz_surf_facet_sphere(&ft, &FR, c, r, k,
                                       fit ? HZ_FIT_MEAN_SAGITTA : HZ_FIT_INSCRIBED, -1, &f0);
      check(n == 20 * (1 << (2 * k)), k == 0 && fit == 0 ? "число граней = 20*4^k" : "");
      int st = 0;
      double v = body_in_box(&ft, f0, n, lo, hi, &st);
      nfac[fit] = st;
      e[fit] = v - exact;
      hz_facettab_free(&ft);
      if (st == HZ_BOX_TOOMANY) toomany = 1;
    }
    if (toomany) {
      printf("    k=%d (%4d граней): режущих БОЛЬШЕ HZ_P3_MAXH=%d — предел ядра,\n"
             "      а не ошибка подгона; одноуровневая фасетизация упирается сюда\n",
             k, 20 * (1 << (2 * k)), HZ_P3_MAXH);
      break;
    }
    kmax = k;
    if (e[0] >= 0.0) ins_one_signed = 0; /* вписанный обязан НЕДОдавать */
    if (k > 0 && fabs(e[0]) > prev_ins) mono_ins = 0;
    if (k > 0 && fabs(e[1]) > prev_fit) mono_fit = 0;
    if (fabs(e[1]) >= fabs(e[0])) fit_better = 0;
    prev_ins = fabs(e[0]);
    prev_fit = fabs(e[1]);
    printf("    k=%d (%4d граней, режут ячейку %3d/%3d): вписанный %+.6f, средняя сагитта %+.6f\n",
           k, 20 * (1 << (2 * k)), nfac[0], nfac[1], e[0], e[1]);
  }
  check(kmax >= 1, "свип дошёл хотя бы до k=1");
  check(ins_one_signed, "вписанный даёт ОДНОЗНАКОВУЮ ошибку (тело всюду меньше)");
  check(mono_ins, "ошибка вписанного падает с k");
  check(mono_fit, "ошибка средней сагитты падает с k");
  check(fit_better, "средняя сагитта ТОЧНЕЕ вписанного на каждом k");
}

/* --------------------------------- 2. Г37: отбор ТОЖДЕСТВЕН, а не приближён */

/* тот же отбор, но по ЦЕНТРУ коробки вместо опорной вершины — так делать нельзя */
static int facets_by_center(const hz_facettab *ft, int32_t f0, int32_t nf, const int32_t lo[3],
                            const int32_t hi[3], int32_t *out, int max) {
  int ns = 0;
  for (int32_t j = 0; j < nf; j++) {
    const hz_facet *f = &ft->f[f0 + j];
    double v = 0.0;
    for (int a = 0; a < 3; a++)
      v += f->n[a] * 0.5 * ((double)lo[a] + (double)hi[a]);
    if (v <= f->off) continue; /* «центр внутри — плоскость лишняя»: НЕВЕРНО */
    if (ns >= max) return HZ_BOX_TOOMANY;
    out[ns++] = f0 + j;
  }
  return ns;
}

static double vol_from_set(const hz_facettab *ft, const int32_t *sel, int ns, const int32_t lo[3],
                           const int32_t hi[3]) {
  /* ноль-инициализация: при ns == 0 (ячейка целиком внутри тела, Г38) цикл ниже
   * не выполняется, и cppcheck справедливо не может доказать, что ядро туда не
   * заглянет */
  hz_hspace h[HZ_P3_MAXH] = {{{0, 0, 0}, 0}};
  int32_t hid[HZ_P3_MAXH] = {0};
  for (int j = 0; j < ns; j++) {
    memcpy(h[j].n, ft->f[sel[j]].n, sizeof h[j].n);
    h[j].off = ft->f[sel[j]].off;
    hid[j] = sel[j];
  }
  hz_poly3 p;
  if (hz_poly3_cut(&p, lo, hi, h, hid, ns) != HZ_P3_OK) return 0.0;
  return hz_poly3_volume(&p, &FR);
}

static void t_selection_exact(void) {
  hz_facettab ft;
  check(hz_facettab_init(&ft) == 0, "init");
  double c[3] = {0, 0, 0};
  int32_t f0 = 0;
  int32_t n = hz_surf_facet_sphere(&ft, &FR, c, 3.0, 1, HZ_FIT_MEAN_SAGITTA, -1, &f0);
  check(n == 80, "80 фасетов");
  int32_t lo[3] = {1, -1, -1}, hi[3] = {3, 1, 1};

  /* полный набор — но он больше HZ_P3_MAXH, поэтому берём k=1 (80 <= 96) */
  int32_t full[HZ_P3_MAXH] = {
      0}; /* см. ноль-инициализацию выше: n доказуемо > 0 только нам, не cppcheck */
  check(n <= HZ_P3_MAXH, "полный набор влезает в ядро (иначе сравнивать не с чем)");
  for (int32_t j = 0; j < n; j++)
    full[j] = f0 + j;
  double v_full = vol_from_set(&ft, full, (int)n, lo, hi);

  int32_t sel[HZ_P3_MAXH];
  int ns = hz_facets_for_box(&ft, f0, n, lo, hi, sel, HZ_P3_MAXH);
  check(ns > 0 && ns < n, "отбор что-то отбросил и что-то оставил");
  double v_sel = vol_from_set(&ft, sel, ns, lo, hi);
  check(memcmp(&v_sel, &v_full, sizeof v_sel) == 0,
        "Г37: отбор ТОЖДЕСТВЕН — объём совпадает ПОБИТОВО с полным набором");

  /* НЕГАТИВНЫЙ КОНТРОЛЬ (предсказан провал): отбор по центру коробки */
  int32_t bad[HZ_P3_MAXH];
  int nb = facets_by_center(&ft, f0, n, lo, hi, bad, HZ_P3_MAXH);
  double v_bad = vol_from_set(&ft, bad, nb, lo, hi);
  check(memcmp(&v_bad, &v_full, sizeof v_bad) != 0,
        "негативный контроль: отбор по центру ОБЯЗАН испортить объём");
  printf("  [Г37] полный %d -> отбор %d: V=%.9f побитово; отбор по центру %d: V=%.9f\n", (int)n, ns,
         v_sel, nb, v_bad);
  hz_facettab_free(&ft);
}

/* ------------------------- 3. Г38: полная и пустая ячейка различимы */

static void t_inside_outside(void) {
  hz_facettab ft;
  check(hz_facettab_init(&ft) == 0, "init");
  double c[3] = {0, 0, 0};
  int32_t f0 = 0;
  int32_t n = hz_surf_facet_sphere(&ft, &FR, c, 8.0, 1, HZ_FIT_MEAN_SAGITTA, -1, &f0);
  check(n == 80, "80 фасетов, r=8");
  int32_t sel[HZ_P3_MAXH];

  int32_t ilo[3] = {-1, -1, -1}, ihi[3] = {1, 1, 1}; /* глубоко внутри шара */
  int in = hz_facets_for_box(&ft, f0, n, ilo, ihi, sel, HZ_P3_MAXH);
  check(in == 0, "Г38: ячейка внутри тела -> НОЛЬ отобранных");

  int32_t olo[3] = {20, 20, 20}, ohi[3] = {22, 22, 22}; /* далеко снаружи */
  int out = hz_facets_for_box(&ft, f0, n, olo, ohi, sel, HZ_P3_MAXH);
  check(out == HZ_BOX_OUTSIDE, "Г38: ячейка снаружи -> HZ_BOX_OUTSIDE, а НЕ ноль");
  check(in != out, "два случая РАЗЛИЧИМЫ — это разница между стеной и воздухом");

  int32_t clo[3] = {5, -1, -1}, chi[3] = {9, 1, 1}; /* режется */
  int cut = hz_facets_for_box(&ft, f0, n, clo, chi, sel, HZ_P3_MAXH);
  check(cut > 0, "режущаяся ячейка -> положительное число");
  printf("  [Г38] внутри=%d, снаружи=%d, режется=%d\n", in, out, cut);

  /* переполнение отбора не молчит */
  int over = hz_facets_for_box(&ft, f0, n, clo, chi, sel, 1);
  check(over == HZ_BOX_TOOMANY, "переполнение отбора возвращает HZ_BOX_TOOMANY, а не обрезает");
  hz_facettab_free(&ft);
}

/* ------------------------------------ 4. dmax заполнен и открывает ворота Г25 */

static void t_dmax(void) {
  hz_facettab ft;
  check(hz_facettab_init(&ft) == 0, "init");
  double c[3] = {0, 0, 0}, r = 3.0;
  int32_t f0 = 0;
  int32_t n0 = hz_surf_facet_sphere(&ft, &FR, c, r, 1, HZ_FIT_INSCRIBED, 0, &f0);
  double dins = 0.0;
  for (int32_t j = 0; j < n0; j++)
    if (ft.f[f0 + j].dmax > dins) dins = ft.f[f0 + j].dmax;
  int32_t f1 = 0;
  int32_t n1 = hz_surf_facet_sphere(&ft, &FR, c, r, 1, HZ_FIT_MEAN_SAGITTA, 0, &f1);
  double dfit = 0.0;
  for (int32_t j = 0; j < n1; j++)
    if (ft.f[f1 + j].dmax > dfit) dfit = ft.f[f1 + j].dmax;
  check(n0 == 80 && n1 == 80, "оба набора построены");
  check(dins > 0.0 && dfit > 0.0, "dmax заполнен, а не оставлен неизвестным");
  check(dfit < dins, "средняя сагитта уменьшает dmax");
  check(fabs(dfit / dins - 0.75) < 0.05, "и уменьшает его примерно на четверть (3/8 против 1/2)");

  double err = 0.0;
  check(hz_facet_error(&ft.f[f1], 1e6, 1.0, &err) == 0,
        "Г25: у фасетизированной сферы ворота ОТКРЫТЫ (dmax известен)");
  int32_t unk = hz_facettab_add_plane(&ft, &FR, (double[]){0, 0, 1}, 1.0, 0, HZ_FACET_DMAX_UNKNOWN);
  check(hz_facet_error(&ft.f[unk], 1e6, 1.0, &err) != 0,
        "и по-прежнему ЗАКРЫТЫ у фасета без вычисленного смещения");
  printf("  [Р-5а] dmax: вписанный %.6f, средняя сагитта %.6f (отношение %.3f)\n", dins, dfit,
         dfit / dins);
  hz_facettab_free(&ft);
}

/* --------------------------- 5. водонепроницаемость на фасетах сферы */

static void t_watertight(void) {
  hz_facettab ft;
  check(hz_facettab_init(&ft) == 0, "init");
  double c[3] = {0.37, 0.11, -0.23}, r = 3.3;
  int32_t f0 = 0;
  int32_t n = hz_surf_facet_sphere(&ft, &FR, c, r, 1, HZ_FIT_MEAN_SAGITTA, -1, &f0);
  check(n == 80, "80 фасетов");
  int32_t alo[3] = {-1, -3, -3}, ahi[3] = {1, 3, 3};
  int32_t blo[3] = {1, -3, -3}, bhi[3] = {3, 3, 3};
  int sa = 0, sb = 0;
  double va = body_in_box(&ft, f0, n, alo, ahi, &sa);
  double vb = body_in_box(&ft, f0, n, blo, bhi, &sb);
  check(sa > 0 && sb > 0, "обе ячейки режутся");

  /* точки на общей стене x=1 обязаны совпасть побитово */
  int32_t sel[HZ_P3_MAXH];
  /* ноль-инициализация: при ns == 0 (ячейка целиком внутри тела, Г38) цикл ниже
   * не выполняется, и cppcheck справедливо не может доказать, что ядро туда не
   * заглянет */
  hz_hspace h[HZ_P3_MAXH] = {{{0, 0, 0}, 0}};
  int32_t hid[HZ_P3_MAXH] = {0};
  hz_poly3 pa, pb;
  int ns = hz_facets_for_box(&ft, f0, n, alo, ahi, sel, HZ_P3_MAXH);
  for (int j = 0; j < ns; j++) {
    memcpy(h[j].n, ft.f[sel[j]].n, sizeof h[j].n);
    h[j].off = ft.f[sel[j]].off;
    hid[j] = sel[j];
  }
  check(hz_poly3_cut(&pa, alo, ahi, h, hid, ns) == HZ_P3_OK, "многогранник A");
  ns = hz_facets_for_box(&ft, f0, n, blo, bhi, sel, HZ_P3_MAXH);
  for (int j = 0; j < ns; j++) {
    memcpy(h[j].n, ft.f[sel[j]].n, sizeof h[j].n);
    h[j].off = ft.f[sel[j]].off;
    hid[j] = sel[j];
  }
  check(hz_poly3_cut(&pb, blo, bhi, h, hid, ns) == HZ_P3_OK, "многогранник B");

  int pairs = 0, exact = 1;
  for (int32_t i = 0; i < pa.nv; i++) {
    if (fabs(pa.v[i][0] - 1.0) > 0.0) continue;
    int best = -1;
    double bd = 1e300;
    for (int32_t j = 0; j < pb.nv; j++) {
      if (fabs(pb.v[j][0] - 1.0) > 0.0) continue;
      double d = fabs(pa.v[i][1] - pb.v[j][1]) + fabs(pa.v[i][2] - pb.v[j][2]);
      if (d < bd) {
        bd = d;
        best = j;
      }
    }
    if (best < 0) continue;
    pairs++;
    if (memcmp(pa.v[i], pb.v[best], 3 * sizeof(double)) != 0) exact = 0;
  }
  check(pairs >= 3, "на общей стене есть что сравнивать");
  check(exact, "фасеты сферы: точки на общей стене совпадают ПОБИТОВО");
  printf("  [Г3+Р-5а] V(A)=%.6f V(B)=%.6f, пар на стене %d, побитово=%d\n", va, vb, pairs, exact);
  hz_facettab_free(&ft);
}

/* -------------------------------------------------------------------------- */

int main(void) {
  printf("=== ПРЕДСКАЗАНИЯ (до единого результата) ===\n");
  printf("  1 метрика — объём ЯЧЕЙКИ против аналитического сегмента, НЕ полной\n");
  printf("    сферы (Г35: общий множитель сделал бы её тождеством). Ошибка\n");
  printf("    вписанного ОДНОЗНАКОВА (недодаёт) и падает с k; средняя сагитта\n");
  printf("    точнее вписанного на каждом k\n");
  printf("  2 Г37: отбор фасетов ТОЖДЕСТВЕН — объём совпадает с полным набором\n");
  printf("    ПОБИТОВО, а не в пределах допуска\n");
  printf("  3 Г38: ячейка внутри -> 0 отобранных; снаружи -> HZ_BOX_OUTSIDE;\n");
  printf("    это РАЗНЫЕ ответы, потому что это стена против воздуха\n");
  printf("  4 dmax заполнен; средняя сагитта уменьшает его в ~0.75 (3/8 против 1/2)\n");
  printf("  5 фасеты сферы водонепроницаемы: общая стена совпадает ПОБИТОВО\n");
  printf("=== НЕГАТИВНЫЕ КОНТРОЛИ, каждый с ПРЕДСКАЗАННЫМ провалом ===\n");
  printf("  Г37 отбор по ЦЕНТРУ коробки вместо опорной вершины\n");
  printf("      -> объём ОБЯЗАН разойтись с полным набором\n");
  printf("  Г25 фасет без вычисленного dmax обязан по-прежнему ЗАКРЫВАТЬ ворота\n");
  printf("=== РЕЗУЛЬТАТЫ ===\n");

  t_convergence();
  t_selection_exact();
  t_inside_outside();
  t_dmax();
  t_watertight();

  printf("%s: %d/%d\n", g_fail ? "FAILURES" : "ok", g_total - g_fail, g_total);
  return g_fail != 0;
}
