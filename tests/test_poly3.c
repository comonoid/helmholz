/* Фальсификаторы общего ядра отсечения (PLAN_CUT.md, Г8 в редакции после
 * Г21…Г24, Г30). ПРЕДСКАЗАНИЯ ПЕЧАТАЮТСЯ ДО РЕЗУЛЬТАТОВ.
 *
 * Главное, что здесь проверяется, — не «сходится ли сумма». Сумма объёмов
 * сходится тождественно, при любой плоскости, включая непрочитанную (Г22),
 * поэтому у каждой измеримой величины есть свой негативный контроль с
 * ПРЕДСКАЗАННЫМ провалом, и он гоняется рядом. */
#include "cut/poly3.h"
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

/* Пол округления. Не 1e-15: объём собирается теоремой о дивергенции суммой
 * ~10^2 знаковых тетраэдров, и несколько десятков ulp здесь штатны. Порог
 * назван и обоснован, а не подобран (CLAUDE.md: магических порогов нет). */
static const double TOL = 1e-14;

static int near_rel(double a, double b, double scale) {
  return fabs(a - b) <= TOL * scale;
}

static const hz_frame FR_UNIT = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};

/* ---------------------------------------------------------------- helpers */

/* Объём внутреннего куска для плоскости n·x <= off над коробкой [lo,hi). */
static double vol_in(const int32_t lo[3], const int32_t hi[3], const double n[3], double off) {
  hz_poly3 p;
  hz_hspace h;
  memcpy(h.n, n, sizeof h.n);
  h.off = off;
  int32_t id = 0;
  int rc = hz_poly3_cut(&p, lo, hi, &h, &id, 1);
  if (rc != HZ_P3_OK) return 0.0;
  return hz_poly3_volume(&p, &FR_UNIT);
}

/* Индекс грани, лежащей на полуплоскости hid (не на коробке). */
static int face_of(const hz_poly3 *p, int32_t hid) {
  for (int f = 0; f < p->nf; f++)
    if (p->fsrc[f] == hid) return f;
  return -1;
}

/* ------------------------------------------------- 1. Σ объёмов и свип */

static void t_volume_sum(void) {
  int32_t lo[3] = {0, 0, 0}, hi[3] = {4, 4, 4};
  double n[3] = {1.0, 1.0, 1.0};
  hz_hspace h;
  memcpy(h.n, n, sizeof h.n);
  h.off = 2.0;
  int32_t id = 0, idf = 1;

  hz_poly3 in;
  check(hz_poly3_cut(&in, lo, hi, &h, &id, 1) == HZ_P3_OK, "cut: угловой тетраэдр построен");
  double vin = hz_poly3_volume(&in, &FR_UNIT);
  check(near_rel(vin, 8.0 / 6.0, 64.0), "V_внутр = a^3/6 замкнутой формы");

  hz_poly3 out[HZ_P3_MAXH];
  int nout = 0;
  check(hz_poly3_complement(out, HZ_P3_MAXH, &nout, lo, hi, &h, &id, &idf, 1) == HZ_P3_OK,
        "complement построен");
  double vout = 0.0;
  for (int j = 0; j < nout; j++)
    vout += hz_poly3_volume(&out[j], &FR_UNIT);
  check(near_rel(vin + vout, 64.0, 64.0), "V_внутр + ΣV_внешн = V_коробки");

  /* Г22: сумма сходится тождественно, поэтому одна она НЕ доказывает, что
   * плоскость вообще прочитана. Свип — вот что доказывает. */
  double prev = -1.0;
  int monotone = 1, moved = 0;
  for (int k = 0; k <= 24; k++) {
    double off = 12.0 * (double)k / 24.0; /* 0 .. 12 = 3*4, вся диагональ */
    double v = vol_in(lo, hi, n, off);
    if (v < prev - TOL * 64.0) monotone = 0;
    if (k > 0 && fabs(v - prev) > 1e-9) moved = 1;
    prev = v;
  }
  check(monotone, "свип off: V_внутр монотонна");
  check(moved, "свип off: V_внутр ВООБЩЕ меняется");
  check(near_rel(vol_in(lo, hi, n, 0.0), 0.0, 64.0), "свип: при off=0 объёма нет");
  check(near_rel(vol_in(lo, hi, n, 12.0), 64.0, 64.0), "свип: при off=12 вся коробка");

  /* НЕГАТИВНЫЙ КОНТРОЛЬ (предсказан провал): читать не off, а константу ⇒
   * зависимость обязана исчезнуть, А СУММА ВСЁ РАВНО СОЙДЁТСЯ. */
  double frozen_spread = 0.0, v0 = vol_in(lo, hi, n, 6.0);
  for (int k = 0; k <= 8; k++) {
    double v = vol_in(lo, hi, n, 6.0); /* off игнорируется — «замороженное чтение» */
    frozen_spread += fabs(v - v0);
  }
  check(frozen_spread <= 0.0, "негативный контроль: при замороженном off разброс НУЛЕВОЙ");
  printf("  [Г22] свип живой: V(0)=%.3f V(6)=%.3f V(12)=%.3f | замороженный разброс %.1e\n",
         vol_in(lo, hi, n, 0.0), v0, vol_in(lo, hi, n, 12.0), frozen_spread);
}

/* ------------------------------- 2a. побитовость на общей грани соседей */

/* Тот же рез в ячейке [lo,hi), но посчитанный в ЛОКАЛЬНЫХ координатах:
 * коробка сдвинута в ноль, плоскость сдвинута вместе с ней, вершины возвращены
 * прибавлением lo. Математически то же самое, в битах — нет. */
static int cut_local(hz_poly3 *p, const int32_t lo[3], const int32_t hi[3], const double n[3],
                     double off) {
  int32_t zlo[3] = {0, 0, 0}, zhi[3];
  double shift = 0.0;
  for (int a = 0; a < 3; a++) {
    zhi[a] = hi[a] - lo[a];
    shift += n[a] * (double)lo[a];
  }
  hz_hspace h;
  memcpy(h.n, n, sizeof h.n);
  h.off = off - shift;
  int32_t id = 0;
  int rc = hz_poly3_cut(p, zlo, zhi, &h, &id, 1);
  if (rc != HZ_P3_OK) return rc;
  for (int32_t i = 0; i < p->nv; i++)
    for (int a = 0; a < 3; a++)
      p->v[i][a] += (double)lo[a];
  return HZ_P3_OK;
}

/* Сколько вершин p лежат ровно на x = xc, и совпадают ли они ПОБИТОВО с
 * вершинами q на той же плоскости. Возвращает число сопоставленных пар. */
static int match_on_wall(const hz_poly3 *p, const hz_poly3 *q, double xc, int *exact) {
  int pairs = 0;
  *exact = 1;
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
    pairs++;
    if (memcmp(p->v[i], q->v[best], 3 * sizeof(double)) != 0) *exact = 0;
  }
  return pairs;
}

static void t_bitwise(void) {
  /* нормаль намеренно «некруглая»: на осевой плоскости совпало бы что угодно */
  double n[3] = {0.3137, 1.0, 0.4271};
  double off = 5.137;
  int32_t alo[3] = {0, 0, 0}, ahi[3] = {4, 4, 4};
  int32_t blo[3] = {4, 0, 0}, bhi[3] = {8, 4, 4};

  hz_poly3 a, b, bloc;
  hz_hspace h;
  memcpy(h.n, n, sizeof h.n);
  h.off = off;
  int32_t id = 0;
  check(hz_poly3_cut(&a, alo, ahi, &h, &id, 1) == HZ_P3_OK, "сосед A построен");
  check(hz_poly3_cut(&b, blo, bhi, &h, &id, 1) == HZ_P3_OK, "сосед B построен");

  int exact = 0;
  int pairs = match_on_wall(&a, &b, 4.0, &exact);
  check(pairs >= 2, "на общей грани есть что сравнивать");
  check(exact, "2a: точки на общем ребре соседей совпадают ПОБИТОВО");

  /* НЕГАТИВНЫЙ КОНТРОЛЬ (предсказан провал): локальные координаты */
  check(cut_local(&bloc, blo, bhi, n, off) == HZ_P3_OK, "локальный вариант построен");
  int exact_loc = 0;
  int pairs_loc = match_on_wall(&a, &bloc, 4.0, &exact_loc);
  check(pairs_loc >= 2, "локальный вариант даёт те же точки геометрически");
  check(!exact_loc, "негативный контроль: в локальных координатах биты ОБЯЗАНЫ разойтись");
  printf("  [Г23] пар на стене %d, побитово=%d; локальные координаты: побитово=%d\n", pairs, exact,
         exact_loc);
}

/* ---------------------------------------- 4. вырождения (Г8 №4, Г30) */

static void t_degenerate(void) {
  int32_t lo[3] = {0, 0, 0}, hi[3] = {4, 4, 4};
  int32_t id = 0;
  hz_poly3 p;
  struct {
    double n[3], off;
    const char *what;
    int expect_ok;
    double vol;
  } c[] = {
      {{1, 1, 1}, 0.0, "плоскость точно через УГОЛ", 0, 0.0},
      {{1, 1, 0}, 0.0, "плоскость точно через РЕБРО", 0, 0.0},
      {{0, 0, 1}, 0.0, "плоскость точно по ГРАНИ (снаружи)", 0, 0.0},
      {{0, 0, 1}, 4.0, "плоскость точно по ГРАНИ (внутри)", 1, 64.0},
      {{0, 0, 1}, 100.0, "плоскость целиком вне ячейки (всё внутри)", 1, 64.0},
      {{0, 0, 1}, -100.0, "плоскость целиком вне ячейки (всё снаружи)", 0, 0.0},
  };
  for (size_t k = 0; k < sizeof c / sizeof c[0]; k++) {
    hz_hspace h;
    memcpy(h.n, c[k].n, sizeof h.n);
    h.off = c[k].off;
    int rc = hz_poly3_cut(&p, lo, hi, &h, &id, 1);
    if (c[k].expect_ok) {
      check(rc == HZ_P3_OK, c[k].what);
      if (rc == HZ_P3_OK) check(near_rel(hz_poly3_volume(&p, &FR_UNIT), c[k].vol, 64.0), c[k].what);
    } else {
      /* куска нулевого объёма НЕ порождается: либо EMPTY, либо нулевой объём */
      check(rc == HZ_P3_EMPTY, c[k].what);
    }
  }

  /* скачок в ОБЕ стороны от вырождения, а не залипание */
  double n[3] = {1, 1, 1};
  double vm = vol_in(lo, hi, n, -0.001), v0 = vol_in(lo, hi, n, 0.0), vp = vol_in(lo, hi, n, 0.001);
  check(vm <= 0.0 && v0 <= 0.0 && vp > 0.0, "у вырождения объём трогается с нуля в ОДНУ сторону");
  printf("  [Г8-4] через угол: V(-e)=%.3g V(0)=%.3g V(+e)=%.3g\n", vm, v0, vp);
}

/* ------------------------------ 5. площадь фасета + 6. перестановка */

static void t_area_and_permute(void) {
  int32_t lo[3] = {0, 0, 0}, hi[3] = {4, 4, 4};
  double n[3] = {1, 1, 1};
  hz_hspace h;
  memcpy(h.n, n, sizeof h.n);
  h.off = 2.0;
  int32_t id = 7;
  hz_poly3 p;
  check(hz_poly3_cut(&p, lo, hi, &h, &id, 1) == HZ_P3_OK, "фасет: рез построен");
  int f = face_of(&p, id);
  check(f >= 0, "фасет найден по fsrc");
  double area = f >= 0 ? hz_poly3_area(&p, f, &FR_UNIT) : 0.0;
  double exact = sqrt(3.0) / 2.0 * 4.0; /* (√3/2)a² при a=2 — сечение куба */
  check(near_rel(area, exact, exact), "5: площадь фасета = замкнутой форме");

  /* НЕГАТИВНЫЙ КОНТРОЛЬ (Г24, предсказан провал): сбить нормаль на 1% */
  double n2[3] = {1.01, 1.0, 1.0};
  hz_hspace h2;
  memcpy(h2.n, n2, sizeof h2.n);
  h2.off = 2.0;
  hz_poly3 p2;
  check(hz_poly3_cut(&p2, lo, hi, &h2, &id, 1) == HZ_P3_OK, "фасет: возмущённый рез построен");
  int f2 = face_of(&p2, id);
  double area2 = f2 >= 0 ? hz_poly3_area(&p2, f2, &FR_UNIT) : 0.0;
  double rel = fabs(area2 - area) / area;
  check(rel > 1e-3, "негативный контроль: сбитая на 1% нормаль двигает площадь на O(1%)");
  printf("  [Г24] площадь %.6f (замкнутая %.6f); при нормали +1%%: %.6f, сдвиг %.3f%%\n", area,
         exact, area2, 100.0 * rel);

  /* 6 (Г30): перестановка полуплоскостей — тот же объём ПОБИТОВО */
  hz_hspace hh[3];
  int32_t ids[3] = {10, 11, 12};
  double na[3][3] = {{1, 0, 0}, {0, 1, 0}, {1, 1, 1}};
  double offs[3] = {3.0, 3.0, 5.0};
  for (int j = 0; j < 3; j++) {
    memcpy(hh[j].n, na[j], sizeof hh[j].n);
    hh[j].off = offs[j];
  }
  hz_poly3 q;
  check(hz_poly3_cut(&q, lo, hi, hh, ids, 3) == HZ_P3_OK, "три полуплоскости");
  double vq = hz_poly3_volume(&q, &FR_UNIT);

  int perm[6][3] = {{0, 1, 2}, {0, 2, 1}, {1, 0, 2}, {1, 2, 0}, {2, 0, 1}, {2, 1, 0}};
  int same = 1;
  for (int k = 0; k < 6; k++) {
    hz_hspace pp[3];
    int32_t pid[3];
    for (int j = 0; j < 3; j++) {
      pp[j] = hh[perm[k][j]];
      pid[j] = ids[perm[k][j]];
    }
    hz_poly3 r;
    if (hz_poly3_cut(&r, lo, hi, pp, pid, 3) != HZ_P3_OK) {
      same = 0;
    } else {
      /* ПОБИТОВО, поэтому memcmp, а не ==: «совпало с точностью до» тут было бы
       * слабее заявленного, и -Wfloat-equal ругался бы по делу */
      double vr = hz_poly3_volume(&r, &FR_UNIT);
      if (memcmp(&vr, &vq, sizeof vr) != 0) same = 0;
    }
  }
  check(same, "6: перестановка полуплоскостей даёт тот же объём ПОБИТОВО");
  printf("  [Г30] V=%.9f при всех 6 перестановках\n", vq);
}

/* --------------------------------------- Г21: метрика на анизотропии */

static void t_anisotropy(void) {
  int32_t lo[3] = {0, 0, 0}, hi[3] = {4, 4, 4};
  double n[3] = {1, 1, 1};
  hz_hspace h;
  memcpy(h.n, n, sizeof h.n);
  h.off = 2.0;
  int32_t id = 5;
  hz_poly3 p;
  check(hz_poly3_cut(&p, lo, hi, &h, &id, 1) == HZ_P3_OK, "анизотропия: рез построен");
  int f = face_of(&p, id);
  check(f >= 0, "анизотропия: фасет найден");
  if (f < 0) return;

  hz_frame iso = {{0, 0, 0}, {1, 1, 1}};
  hz_frame ani = {{0, 0, 0}, {1, 1, 8}};
  double a_iso = hz_poly3_area(&p, f, &iso);
  double a_ani = hz_poly3_area(&p, f, &ani);
  double v_iso = hz_poly3_volume(&p, &iso);
  double v_ani = hz_poly3_volume(&p, &ani);

  /* объём переносится умножением ТОЧНО, площадь — нет */
  check(near_rel(v_ani, v_iso * 8.0, v_iso * 8.0), "объём в мир переносится множителем u0u1u2");
  check(fabs(a_ani - a_iso) / a_iso > 0.1,
        "негативный контроль Г21: на ячейке 1x1x8 площадь ОБЯЗАНА отличаться от единичной");

  double nrm[3];
  hz_poly3_normal(&p, f, &ani, nrm);
  /* у плоскости x+y+z=2 в кадре с u=(1,1,8) нормаль в мире ∝ (8,8,1)/... */
  double ex[3] = {8.0, 8.0, 1.0};
  double m = sqrt(ex[0] * ex[0] + ex[1] * ex[1] + ex[2] * ex[2]);
  int nok = 1;
  for (int a = 0; a < 3; a++)
    if (!near_rel(nrm[a], ex[a] / m, 1.0)) nok = 0;
  check(nok, "нормаль в анизотропном кадре считается в МИРЕ, а не в единицах");
  printf("  [Г21] площадь: единичный кадр %.6f, кадр 1x1x8 %.6f (отношение %.3f)\n", a_iso, a_ani,
         a_ani / a_iso);
}

/* ------------------------------------------------------------------------ */

int main(void) {
  printf("=== ПРЕДСКАЗАНИЯ (до единого результата) ===\n");
  printf("  1  V_внутр + ΣV_внешн = V_коробки, пол округления %.0e отн.\n", TOL);
  printf("     ТОЖДЕСТВО: одно оно ничего не проверяет, поэтому рядом свип off —\n");
  printf("     V_внутр обязана монотонно пройти 0 -> V_коробки (Г22)\n");
  printf("  2a точки на общем ребре соседей ОДНОГО уровня совпадают ПОБИТОВО;\n");
  printf("     2b (разные уровни) здесь НЕ проверяется — это шов пункта 5 (Г23)\n");
  printf("  4  вырождения (угол/ребро/грань/вне): куска нулевого объёма НЕ\n");
  printf("     порождается, число кусков считается только для ВНУТРЕННЕГО (Г30)\n");
  printf("  5  площадь фасета = (sqrt(3)/2)a^2 замкнутой формы\n");
  printf("  6  перестановка полуплоскостей: объём совпадает ПОБИТОВО (Г30)\n");
  printf("=== НЕГАТИВНЫЕ КОНТРОЛИ, каждый с ПРЕДСКАЗАННЫМ провалом ===\n");
  printf("  Г22 замороженное чтение off  -> разброс V_внутр обязан стать НУЛЕВЫМ,\n");
  printf("      а сумма объёмов при этом всё равно сойдётся\n");
  printf("  Г23 счёт в локальных координатах -> биты ОБЯЗАНЫ разойтись\n");
  printf("  Г24 нормаль сбита на 1%%       -> площадь обязана сдвинуться на O(1%%)\n");
  printf("  Г21 кадр 1x1x8 против 1x1x1   -> площадь обязана отличиться, объём нет\n");
  printf("=== РЕЗУЛЬТАТЫ ===\n");

  t_volume_sum();
  t_bitwise();
  t_degenerate();
  t_area_and_permute();
  t_anisotropy();

  printf("%s: %d/%d\n", g_fail ? "FAILURES" : "ok", g_total - g_fail, g_total);
  return g_fail != 0;
}
