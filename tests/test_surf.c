/* Фальсификаторы таблиц общего слоя (PLAN_CUT.md, Р-3; Г3, Г4, Г20, Г25).
 * ПРЕДСКАЗАНИЯ ПЕЧАТАЮТСЯ ДО РЕЗУЛЬТАТОВ.
 *
 * Главное, что здесь проверяется, — не «складываются ли записи в массив».
 * Проверяются три свойства, каждое из которых при неверной реализации
 * молчаливо: водонепроницаемость держится на ССЫЛКЕ НА ОДИН ОБЪЕКТ (Г3), ключ
 * не должен знать про дерево (Г20), а невычисленное смещение обязано ЗАКРЫВАТЬ
 * дорогу, а не подставлять ноль (Г25). */
#include "cut/surf.h"
#include "octree.h"
#include <complex.h>
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

/* ------------------------------------------ 1. сторона кодируется ~, не - */

static void t_side_encoding(void) {
  hz_facettab ft;
  check(hz_facettab_init(&ft) == 0, "facettab init");
  double n[3] = {0.0, 0.0, 1.0};
  int32_t f0 = hz_facettab_add_plane(&ft, &FR, n, 2.0, -1, 0.0);
  check(f0 == 0, "ФАСЕТ НОЛЬ — именно тот случай, на котором ломается минус");

  hz_cutmap m;
  check(hz_cutmap_init(&m) == 0, "cutmap init");
  int32_t le = f0, ge = ~f0; /* ~0 == -1, а -0 == 0 — вот в чём разница */
  check(ge == -1, "~0 == -1: у нулевого фасета вторая сторона представима");
  check(hz_cutmap_add(&m, 10, &le, 1) == 0, "запись со стороной <=");
  check(hz_cutmap_add(&m, 11, &ge, 1) == 0, "запись со стороной >=");

  hz_hspace h_le, h_ge;
  int32_t id_le, id_ge;
  check(hz_cutmap_hspaces(&ft, &m, hz_cutmap_find(&m, 10), &h_le, &id_le, NULL, 1) == 1,
        "разворот стороны <=");
  check(hz_cutmap_hspaces(&ft, &m, hz_cutmap_find(&m, 11), &h_ge, &id_ge, NULL, 1) == 1,
        "разворот стороны >=");
  int opposite = 1;
  for (int a = 0; a < 3; a++)
    if (memcmp(&h_le.n[a], &(double){-h_ge.n[a]}, sizeof(double)) != 0) opposite = 0;
  check(opposite, "две стороны нулевого фасета — точное отрицание друг друга");
  check(id_le == 0 && id_ge == -1, "hid несёт сам fref, то есть фасет И сторону");
  printf("  [Р-3] фасет 0: <= даёт n_z=%+.1f, >= даёт n_z=%+.1f\n", h_le.n[2], h_ge.n[2]);

  hz_cutmap_free(&m);
  hz_facettab_free(&ft);
}

/* --------------------- 2. Г3: водонепроницаемость от ОБЩЕЙ ссылки */

/* Объём материала в ячейке [lo,hi), отсечённой одним фасетом со стороной. */
static int cell_poly(hz_poly3 *p, const hz_facettab *ft, const hz_cutmap *m, int32_t cell,
                     const int32_t lo[3], const int32_t hi[3]) {
  const hz_cutrec *r = hz_cutmap_find(m, cell);
  if (r == NULL) return -1;
  hz_hspace h[HZ_P3_MAXH];
  int32_t hid[HZ_P3_MAXH];
  int nh = hz_cutmap_hspaces(ft, m, r, h, hid, NULL, HZ_P3_MAXH);
  if (nh < 0) return -1;
  return hz_poly3_cut(p, lo, hi, h, hid, nh);
}

/* сколько вершин p лежат ровно на x = xc и совпадают ли ПОБИТОВО с вершинами q */
static int wall_exact(const hz_poly3 *p, const hz_poly3 *q, double xc, int *pairs) {
  int exact = 1;
  *pairs = 0;
  for (int32_t i = 0; i < p->nv; i++) {
    if (fabs(p->v[i][0] - xc) > 0.0) continue;
    int best = -1;
    double bd = 1e300;
    for (int32_t j = 0; j < q->nv; j++) {
      if (fabs(q->v[j][0] - xc) > 0.0) continue;
      double d = fabs(p->v[i][1] - q->v[j][1]) + fabs(p->v[i][2] - q->v[j][2]);
      if (d < bd) {
        bd = d;
        best = j;
      }
    }
    if (best < 0) continue;
    (*pairs)++;
    if (memcmp(p->v[i], q->v[best], 3 * sizeof(double)) != 0) exact = 0;
  }
  return exact;
}

static void t_shared_reference(void) {
  int32_t alo[3] = {0, 0, 0}, ahi[3] = {4, 4, 4};
  int32_t blo[3] = {4, 0, 0}, bhi[3] = {8, 4, 4};
  double n[3] = {0.3137, 1.0, 0.4271};
  double off = 5.137;

  /* ОБЩАЯ ссылка: один фасет, две ячейки */
  hz_facettab ft;
  hz_cutmap m;
  check(hz_facettab_init(&ft) == 0 && hz_cutmap_init(&m) == 0, "init (общая ссылка)");
  int32_t f = hz_facettab_add_plane(&ft, &FR, n, off, -1, 0.0);
  check(hz_cutmap_add(&m, 100, &f, 1) == 0 && hz_cutmap_add(&m, 200, &f, 1) == 0, "две ячейки");
  hz_poly3 pa, pb;
  check(cell_poly(&pa, &ft, &m, 100, alo, ahi) == HZ_P3_OK, "ячейка A");
  check(cell_poly(&pb, &ft, &m, 200, blo, bhi) == HZ_P3_OK, "ячейка B");
  int pairs = 0;
  int exact = wall_exact(&pa, &pb, 4.0, &pairs);
  check(pairs >= 2, "на общей стене есть что сравнивать");
  check(exact, "Г3: при ОБЩЕЙ ссылке точки на стене совпадают ПОБИТОВО");

  /* НЕГАТИВНЫЙ КОНТРОЛЬ (предсказан провал): у каждой ячейки СВОЙ фасет, и его
   * смещение отличается на один ulp — ровно то, что даёт независимая подгонка */
  hz_facettab ft2;
  hz_cutmap m2;
  check(hz_facettab_init(&ft2) == 0 && hz_cutmap_init(&m2) == 0, "init (свой фасет)");
  int32_t g0 = hz_facettab_add_plane(&ft2, &FR, n, off, -1, 0.0);
  int32_t g1 = hz_facettab_add_plane(&ft2, &FR, n, nextafter(off, 1e300), -1, 0.0);
  check(g0 == 0 && g1 == 1, "два независимых фасета заведены");
  check(hz_cutmap_add(&m2, 100, &g0, 1) == 0 && hz_cutmap_add(&m2, 200, &g1, 1) == 0, "две ячейки");
  hz_poly3 qa, qb;
  check(cell_poly(&qa, &ft2, &m2, 100, alo, ahi) == HZ_P3_OK, "ячейка A (свой фасет)");
  check(cell_poly(&qb, &ft2, &m2, 200, blo, bhi) == HZ_P3_OK, "ячейка B (свой фасет)");
  int pairs2 = 0;
  int exact2 = wall_exact(&qa, &qb, 4.0, &pairs2);
  check(pairs2 >= 2, "и там есть что сравнивать");
  check(!exact2, "негативный контроль: независимые фасеты ОБЯЗАНЫ разойтись в битах");
  printf("  [Г3] общая ссылка: побитово=%d; независимые фасеты (+1 ulp): побитово=%d\n", exact,
         exact2);

  hz_cutmap_free(&m);
  hz_facettab_free(&ft);
  hz_cutmap_free(&m2);
  hz_facettab_free(&ft2);
}

/* ------------------------- 3. Г20: ключ не знает про дерево */

static void t_host_agnostic(void) {
  hz_facettab ft;
  check(hz_facettab_init(&ft) == 0, "init");
  double n[3] = {0.0, 0.0, 1.0};
  int32_t f = hz_facettab_add_plane(&ft, &FR, n, 2.0, -1, 0.0);
  int32_t lo[3] = {0, 0, 0}, hi[3] = {4, 4, 4};

  /* хозяин 1: октодерево, ключ = индекс УЗЛА (hz_oct_leaf) */
  hz_octree t;
  check(hz_oct_init(&t, 3, CMPLX(1.0, 0.0)) == 0, "octree init");
  int blo[3] = {0, 0, 0}, bhi[3] = {4, 4, 4};
  check(hz_oct_set_box(&t, blo, bhi, CMPLX(2.0, 0.0)) == 0, "octree set_box");
  int32_t node = hz_oct_leaf(&t, 1, 1, 1);
  check(node >= 0, "лист найден");
  hz_cutmap m1;
  check(hz_cutmap_init(&m1) == 0, "cutmap (октодерево)");
  check(hz_cutmap_add(&m1, node, &f, 1) == 0, "запись по индексу УЗЛА");
  hz_poly3 p1;
  check(cell_poly(&p1, &ft, &m1, node, lo, hi) == HZ_P3_OK, "рез по ключу-узлу");
  double v1 = hz_poly3_volume(&p1, &FR);

  /* хозяин 2: сетка переноса, ключ = ПОРЯДКОВЫЙ номер ячейки. Ни строки про
   * октодерево — вот содержание Г20, и оно проверяемо, а не декларативно. */
  hz_cutmap m2;
  check(hz_cutmap_init(&m2) == 0, "cutmap (сетка)");
  int32_t mesh_cell = 7;
  check(hz_cutmap_add(&m2, mesh_cell, &f, 1) == 0, "запись по индексу ЯЧЕЙКИ СЕТКИ");
  hz_poly3 p2;
  check(cell_poly(&p2, &ft, &m2, mesh_cell, lo, hi) == HZ_P3_OK, "рез по ключу-ячейке");
  double v2 = hz_poly3_volume(&p2, &FR);

  check(memcmp(&v1, &v2, sizeof v1) == 0, "Г20: оба хозяина дают ПОБИТОВО один результат");
  printf("  [Г20] ключ-узел %d и ключ-ячейка %d: V=%.6f обе\n", node, mesh_cell, v1);

  hz_cutmap_free(&m1);
  hz_cutmap_free(&m2);
  hz_facettab_free(&ft);
  hz_oct_free(&t);
}

/* --------------------------- 4. курсор против двоичного поиска */

static void t_cursor(void) {
  hz_cutmap m;
  check(hz_cutmap_init(&m) == 0, "init");
  int32_t f = 0;
  enum { N = 500 };
  for (int32_t k = 0; k < N; k++)
    check(hz_cutmap_add(&m, k * 3, &f, 1) == 0, k == 0 ? "заполнение" : "");
  g_total -= (N - 1); /* считаем заполнение одной проверкой, а не пятьюстами */

  check(hz_cutmap_add(&m, 0, &f, 1) == 2, "порядок ключей — ИНВАРИАНТ, нарушение отвергается");

  hz_cutcur c;
  hz_cutcur_begin(&c, &m);
  int agree = 1;
  for (int32_t cell = 0; cell < N * 3; cell++) {
    const hz_cutrec *a = hz_cutcur_next(&c, cell);
    const hz_cutrec *b = hz_cutmap_find(&m, cell);
    if (a != b) agree = 0;
  }
  check(agree, "курсор и двоичный поиск дают ОДНО И ТО ЖЕ на всём диапазоне");
  check(c.i <= m.nr, "курсор прошёл массив не более одного раза (слияние, а не поиск)");
  printf("  [Р-3] %d записей, курсор сдвинулся на %d — O(1) на ячейку\n", m.nr, c.i);
  hz_cutmap_free(&m);
}

/* ------------------------------ 5. Г25: невычисленное смещение закрывает */

static void t_dmax_fail_closed(void) {
  hz_facettab ft;
  check(hz_facettab_init(&ft) == 0, "init");
  double n[3] = {0.0, 0.0, 1.0};
  int32_t exact_f = hz_facettab_add_plane(&ft, &FR, n, 2.0, -1, 0.0);
  int32_t derived = hz_facettab_add_plane(&ft, &FR, n, 2.0, 3, HZ_FACET_DMAX_UNKNOWN);
  int32_t fitted = hz_facettab_add_plane(&ft, &FR, n, 2.0, 3, 1e-4);

  double err = -12345.0, keep = err;
  check(hz_facet_error(&ft.f[derived], 1e7, 1.0, &err) != 0,
        "Г25: невычисленное смещение — ОТКАЗ, а не ноль");
  check(memcmp(&err, &keep, sizeof err) == 0, "Г25: при отказе *err не трогается вовсе");

  /* НЕГАТИВНЫЙ КОНТРОЛЬ: отказ обязан быть ИЗБИРАТЕЛЬНЫМ, иначе проверка
   * проходится реализацией `return 1;` */
  check(hz_facet_error(&ft.f[exact_f], 1e7, 1.0, &err) == 0,
        "точный фасет (dmax=0) — пропускается");
  check(err <= 0.0, "у точного фасета ошибка нулевая");
  check(hz_facet_error(&ft.f[fitted], 1e7, 1.0, &err) == 0, "подогнанный фасет — пропускается");
  double expect = 0.125 * (1e7 * 1e-4) * (1e7 * 1e-4) * sqrt(1e-4 / 1.0);
  check(fabs(err - expect) <= 1e-12 * expect, "закон 0.125(kd)^2 sqrt(d/W) воспроизведён");
  printf("  [Г25] dmax неизвестен -> отказ; dmax=0 -> err=0; dmax=1e-4 -> err=%.6f\n", err);
  hz_facettab_free(&ft);
}

/* --------------------- 6. сквозной прогон: материал + дополнение */

static void t_end_to_end(void) {
  hz_surftab st;
  hz_facettab ft;
  hz_cutmap m;
  check(hz_surftab_init(&st) == 0 && hz_facettab_init(&ft) == 0 && hz_cutmap_init(&m) == 0, "init");

  /* поверхность несёт ДВЕ НЕПРОЗРАЧНЫЕ метки; общий слой в них не смотрит */
  double pn[3] = {1.0, 1.0, 1.0};
  hz_surf s = {HZ_SURF_PLANE, {1.0, 1.0, 1.0, 2.0, 0, 0, 0}, 77, 88};
  int32_t si = hz_surftab_add(&st, &s);
  check(si == 0, "поверхность заведена");
  int32_t f = hz_facettab_add_plane(&ft, &FR, pn, 2.0, si, 0.0);
  check(hz_cutmap_add(&m, 42, &f, 1) == 0, "ячейка 42 несёт границу");

  int32_t lo[3] = {0, 0, 0}, hi[3] = {4, 4, 4};
  const hz_cutrec *r = hz_cutmap_find(&m, 42);
  hz_hspace h[HZ_P3_MAXH];
  int32_t hid[HZ_P3_MAXH], hidf[HZ_P3_MAXH];
  int nh = hz_cutmap_hspaces(&ft, &m, r, h, hid, hidf, HZ_P3_MAXH);
  check(nh == 1, "развёрнуто одно полупространство");

  hz_poly3 inside;
  check(hz_poly3_cut(&inside, lo, hi, h, hid, nh) == HZ_P3_OK, "материал построен");
  double vin = hz_poly3_volume(&inside, &FR);
  hz_poly3 outside[HZ_P3_MAXH];
  int nout = 0;
  check(hz_poly3_complement(outside, HZ_P3_MAXH, &nout, lo, hi, h, hid, hidf, nh) == HZ_P3_OK,
        "внешность построена");
  double vout = 0.0;
  for (int j = 0; j < nout; j++)
    vout += hz_poly3_volume(&outside[j], &FR);
  check(fabs(vin + vout - 64.0) <= 1e-14 * 64.0, "материал + внешность = ячейка");

  /* метка стороны берётся из ПОВЕРХНОСТИ, а не из ядра */
  int ff = -1;
  for (int32_t k = 0; k < inside.nf; k++)
    if (inside.fsrc[k] == f) ff = (int)k;
  check(ff >= 0, "грань материала найдена по fsrc = fref");
  const hz_facet *fa = &ft.f[f];
  check(fa->surf == si, "фасет помнит свой примитив");
  check(st.s[fa->surf].tag_le == 77 && st.s[fa->surf].tag_ge == 88,
        "метки сторон доехали нетронутыми — общий слой в них не смотрел");
  printf("  [сквозной] V_материала=%.6f, кусков внешности %d, ΣV=%.6f, метки %d/%d\n", vin, nout,
         vin + vout, st.s[fa->surf].tag_le, st.s[fa->surf].tag_ge);

  hz_cutmap_free(&m);
  hz_facettab_free(&ft);
  hz_surftab_free(&st);
}

/* ------------------------------ 7. перевод в единицы делается и делается раз */

static void t_frame_once(void) {
  hz_facettab ft;
  check(hz_facettab_init(&ft) == 0, "init");
  double n[3] = {0.0, 0.0, 1.0};
  hz_frame iso = {{0, 0, 0}, {1, 1, 1}};
  hz_frame ani = {{0, 0, 10.0}, {1, 1, 8}};
  int32_t a = hz_facettab_add_plane(&ft, &iso, n, 2.0, -1, 0.0);
  int32_t b = hz_facettab_add_plane(&ft, &ani, n, 2.0, -1, 0.0);
  check(memcmp(&ft.f[a].off, &ft.f[b].off, sizeof(double)) != 0,
        "кадр ДЕЙСТВУЕТ: тот же мировой off даёт разное в единицах");
  double eoff = 2.0 - 10.0, en = 8.0;
  check(memcmp(&ft.f[b].off, &eoff, sizeof(double)) == 0, "off' = off - n·o");
  check(memcmp(&ft.f[b].n[2], &en, sizeof(double)) == 0, "n'_z = n_z·u_z");
  printf("  [Р-4] кадр (o_z=10,u_z=8): мировая z<=2 стала %.1f·z <= %.1f\n", ft.f[b].n[2],
         ft.f[b].off);
  hz_facettab_free(&ft);
}

/* -------------------------------------------------------------------------- */

int main(void) {
  printf("=== ПРЕДСКАЗАНИЯ (до единого результата) ===\n");
  printf("  1 сторона кодируется ~fref: у ФАСЕТА НОЛЬ обе стороны представимы\n");
  printf("    (с минусом -0 == 0 и сторона терялась бы)\n");
  printf("  2 Г3: две ячейки по ОБЩЕЙ ссылке дают на стене ПОБИТОВО одни точки\n");
  printf("  3 Г20: тот же результат по ключу-узлу дерева и по ключу-ячейке сетки,\n");
  printf("    и в surf.c нет ни строки про октодерево\n");
  printf("  4 курсор и двоичный поиск совпадают, курсор проходит массив однажды\n");
  printf("  5 Г25: dmax не вычислен -> ОТКАЗ и *err не трогается\n");
  printf("  6 сквозной: V_материала + ΣV_внешности = V_ячейки, метки сторон целы\n");
  printf("=== НЕГАТИВНЫЕ КОНТРОЛИ, каждый с ПРЕДСКАЗАННЫМ провалом ===\n");
  printf("  Г3  у каждой ячейки СВОЙ фасет, off отличается на 1 ulp\n");
  printf("      -> точки на стене ОБЯЗАНЫ разойтись в битах\n");
  printf("  Г25 фасеты с dmax=0 и dmax=1e-4 ОБЯЗАНЫ проходить\n");
  printf("      -> иначе проверку удовлетворяет `return 1;`\n");
  printf("=== РЕЗУЛЬТАТЫ ===\n");

  t_side_encoding();
  t_shared_reference();
  t_host_agnostic();
  t_cursor();
  t_dmax_fail_closed();
  t_end_to_end();
  t_frame_once();

  printf("%s: %d/%d\n", g_fail ? "FAILURES" : "ok", g_total - g_fail, g_total);
  return g_fail != 0;
}
