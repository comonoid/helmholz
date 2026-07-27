/* Фальсификаторы интегратора переноса над разрезанной ячейкой
 * (PLAN_TRANSPORT.md T4, PLAN_CUT.md пункт 6). ПРЕДСКАЗАНИЯ ДО РЕЗУЛЬТАТОВ. */

#include "cut/surf.h"
#include "transport/tet3.h"
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

/* --------------------------------------- 1. объём: тетраэдры против ядра */

static void t_volume(void) {
  const hz_frame FR = {{0, 0, 0}, {1, 1, 1}};
  const hz_frame FRA = {{0, 0, 0}, {1, 1, 8}}; /* анизотропный: Г21 */
  int32_t lo[3] = {0, 0, 0}, hi[3] = {4, 4, 4};
  double worst = 0.0;
  int ncut = 0;
  /* набор косых плоскостей, режущих коробку по-разному */
  for (int k = 0; k < 6; k++) {
    double nn[3] = {0.3 + 0.1 * k, 0.7 - 0.05 * k, 0.5 + 0.07 * k};
    double m = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
    for (int a = 0; a < 3; a++)
      nn[a] /= m;
    hz_hspace h[2];
    int32_t hid[2] = {0, 1};
    memcpy(h[0].n, nn, sizeof nn);
    h[0].off = nn[0] * 2.1 + nn[1] * 1.9 + nn[2] * 2.3;
    double n2[3] = {-0.4, 0.6 + 0.03 * k, -0.5};
    m = sqrt(n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]);
    for (int a = 0; a < 3; a++)
      n2[a] /= m;
    memcpy(h[1].n, n2, sizeof n2);
    h[1].off = n2[0] * 1.7 + n2[1] * 2.2 + n2[2] * 1.6;

    hz_poly3 p;
    if (hz_poly3_cut(&p, lo, hi, h, hid, 2) != HZ_P3_OK) continue;
    ncut++;
    tr3_tet tet[TR3_MAXTET];
    int nt = tr3_tet_from_poly(&p, tet, TR3_MAXTET);
    check(nt > 0, "разбиение непусто");
    double vsum = 0.0, vmin = 1e300;
    for (int t = 0; t < nt; t++) {
      double v = tr3_tet_volume(&tet[t]);
      vsum += v;
      if (v < vmin) vmin = v;
    }
    double vk = hz_poly3_volume(&p, &FR);
    double e = fabs(vsum - vk) / fabs(vk);
    if (e > worst) worst = e;
    check(vmin > 0.0, "все тетраэдры ОДНОГО знака — разбиение ориентировано верно");
    /* та же ячейка в анизотропном кадре: объём переносится множителем ТОЧНО */
    double m00[4][4], cc[3] = {2, 2, 2};
    tr3_mass_matrix(&p, &FRA, cc, 4.0, m00);
    double vka = hz_poly3_volume(&p, &FRA);
    check(fabs(m00[0][0] - vka) < 1e-12 * fabs(vka),
          "M00 матрицы масс ЕСТЬ объём, и на анизотропном кадре тоже");
  }
  printf("  [объём] %d конфигураций: max отн. расхождение тетраэдров с ядром %.2e\n", ncut, worst);
  check(ncut >= 4, "конфигураций достаточно");
  check(worst < 1e-13, "сумма объёмов тетраэдров совпадает с ядром");
}

/* ------------------------------- 2. моменты: против замкнутой формы на КОРОБКЕ */

static void t_moments(void) {
  const hz_frame FR = {{0, 0, 0}, {1, 1, 1}};
  int32_t lo[3] = {0, 0, 0}, hi[3] = {2, 2, 2};
  hz_poly3 p;
  check(hz_poly3_box(&p, lo, hi) == HZ_P3_OK, "коробка");
  double c[3] = {1, 1, 1}, h = 2.0;
  double m[4][4];
  check(tr3_mass_matrix(&p, &FR, c, h, m) == 0, "матрица масс посчитана");

  /* Замкнутая форма для куба со стороной h и центром c при ξ=(x−c)/h:
   *   ∫1 = h³,  ∫ξ = 0,  ∫ξ² = h³/12,  ∫ξη = 0 */
  double v = h * h * h;
  double worst = 0.0;
  double exact[4][4] = {{v, 0, 0, 0}, {0, v / 12, 0, 0}, {0, 0, v / 12, 0}, {0, 0, 0, v / 12}};
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      double e = fabs(m[i][j] - exact[i][j]) / v;
      if (e > worst) worst = e;
    }
  printf("  [моменты] коробка: max |M − замкнутая форма| / V = %.2e\n", worst);
  check(worst < 1e-14, "матрица масс на коробке совпадает с замкнутой формой");

  /* линейный интеграл: ∫(a0 + a·ξ) = a0·V, наклоны обязаны сократиться */
  double a[4] = {3.0, 5.0, -7.0, 11.0};
  double il = tr3_int_linear(&p, &FR, c, h, a);
  printf("  [линейный] ∫(3 + 5ξ − 7η + 11ζ) = %.15f, замкнутая форма %.15f\n", il, 3.0 * v);
  check(fabs(il - 3.0 * v) < 1e-13 * v, "наклоны на симметричной ячейке сокращаются ТОЧНО");
}

/* ------------------- 3. разрезанная ячейка против аналитического сегмента */

static void t_cut_vs_analytic(void) {
  const hz_frame FR = {{0, 0, 0}, {1, 1, 1}};
  /* ячейка, отсечённая ОДНОЙ плоскостью: объём и центроид — замкнутая форма */
  int32_t lo[3] = {0, 0, 0}, hi[3] = {2, 2, 2};
  hz_hspace h1;
  double nn[3] = {0, 0, 1};
  memcpy(h1.n, nn, sizeof nn);
  h1.off = 1.3; /* материал z ≤ 1.3 */
  int32_t hid = 0;
  hz_poly3 p;
  check(hz_poly3_cut(&p, lo, hi, &h1, &hid, 1) == HZ_P3_OK, "отсечение");
  double c[3] = {1, 1, 1}, hh = 2.0;
  double m[4][4];
  tr3_mass_matrix(&p, &FR, c, hh, m);
  double vex = 2.0 * 2.0 * 1.3; /* плита 2×2×1.3 */
  /* центроид по z: 0.65, значит ∫ζ = V·(0.65−1)/2 */
  double zex = vex * (0.65 - 1.0) / 2.0;
  printf("  [срез] V = %.15f против %.15f;  ∫ζ = %.15f против %.15f\n", m[0][0], vex, m[0][3], zex);
  check(fabs(m[0][0] - vex) < 1e-13 * vex, "объём среза точен");
  check(fabs(m[0][3] - zex) < 1e-13 * vex, "первый момент среза точен — центроид на месте");
}

/* --------------------------- 4. негативный контроль: чужая вершина веера */

static void t_control(void) {
  const hz_frame FR = {{0, 0, 0}, {1, 1, 1}};
  int32_t lo[3] = {0, 0, 0}, hi[3] = {4, 4, 4};
  hz_hspace h1;
  double nn[3] = {0.5773502691896258, 0.5773502691896258, 0.5773502691896258};
  memcpy(h1.n, nn, sizeof nn);
  h1.off = nn[0] * 5.0;
  int32_t hid = 0;
  hz_poly3 p;
  check(hz_poly3_cut(&p, lo, hi, &h1, &hid, 1) == HZ_P3_OK, "отсечение");
  double vk = hz_poly3_volume(&p, &FR);
  tr3_tet tet[TR3_MAXTET];
  int nt = tr3_tet_from_poly(&p, tet, TR3_MAXTET);
  double vs = 0.0;
  for (int t = 0; t < nt; t++)
    vs += tr3_tet_volume(&tet[t]);
  printf("  [контроль] правильный веер: %d тетраэдров, V = %.15f против ядра %.15f\n", nt, vs, vk);
  check(fabs(vs - vk) < 1e-13 * vk, "правильный веер даёт объём ядра");

  /* НЕГАТИВНЫЙ КОНТРОЛЬ (предсказан провал): веер от ТОЧКИ ВНЕ тела. Объём со
   * знаком обязан устоять (теорема о дивергенции его не замечает), а вот
   * ТЕТРАЭДРЫ ОБЯЗАНЫ СТАТЬ РАЗНОЗНАКОВЫМИ — то есть разбиение перестанет быть
   * разбиением. Признак ловится именно знаком, а не суммой. */
  double apex[3] = {-5.0, -5.0, -5.0};
  double vs2 = 0.0;
  int neg = 0, pos = 0;
  for (int32_t f = 0; f < p.nf; f++) {
    int32_t b0 = p.floff[f], k = p.floff[f + 1] - b0;
    for (int32_t e = 1; e + 1 < k; e++) {
      tr3_tet t;
      for (int cc = 0; cc < 3; cc++) {
        t.v[0][cc] = apex[cc];
        t.v[1][cc] = p.v[p.fl[b0]][cc];
        t.v[2][cc] = p.v[p.fl[b0 + e]][cc];
        t.v[3][cc] = p.v[p.fl[b0 + e + 1]][cc];
      }
      double v = tr3_tet_volume(&t);
      vs2 += v;
      if (v < 0.0)
        neg++;
      else
        pos++;
    }
  }
  printf("    веер от точки ВНЕ тела: V = %.15f (сумма устояла), но знаков +%d/−%d\n", vs2, pos,
         neg);
  check(fabs(vs2 - vk) < 1e-12 * vk, "сумма со знаком устояла — она этого не ловит");
  check(neg > 0 && pos > 0,
        "негативный контроль: у чужой вершины тетраэдры РАЗНОЗНАКОВЫ, и ловится это знаком");
}

/* -------------------------------------------------------------------------- */

int main(void) {
  printf("=== ПРЕДСКАЗАНИЯ (до единого результата) ===\n");
  printf("  1 сумма объёмов тетраэдров = объём ядра, ≤1e-13, и все ОДНОГО знака\n");
  printf("  2 матрица масс на коробке = замкнутая форма (V, V/12 по диагонали)\n");
  printf("  3 срез плоскостью: объём и первый момент точны — центроид на месте\n");
  printf("  4 M00 ЕСТЬ объём, в том числе на анизотропном кадре (Г21)\n");
  printf("=== НЕГАТИВНЫЙ КОНТРОЛЬ, с ПРЕДСКАЗАННЫМ провалом ===\n");
  printf("  веер от точки ВНЕ тела: СУММА устоит (дивергенция её не различает),\n");
  printf("  а знаки тетраэдров ОБЯЗАНЫ разойтись — вот чем это ловится\n");
  printf("=== РЕЗУЛЬТАТЫ ===\n");

  t_volume();
  t_moments();
  t_cut_vs_analytic();
  t_control();

  printf("%s: %d/%d\n", g_fail ? "FAILURES" : "ok", g_total - g_fail, g_total);
  return g_fail != 0;
}
