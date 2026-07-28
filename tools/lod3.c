/* lod3.c — ПРОВЕРКА ЗАКОНА LOD СЧЁТОМ ЯЧЕЕК, БЕЗ ВСЯКОГО РЕШЕНИЯ.
 *
 * PLAN_TRANSPORT.md. Центральное утверждение архитектуры звучит так: размер
 * элемента задаётся угловым размером пикселя, `L = εR`, и тогда на КАЖДУЮ
 * ОКТАВУ расстояния приходится ОДНО И ТО ЖЕ число элементов, `~4π/ε²`. Отсюда
 * полное число растёт как `log₂(R_max/R_min)`, а не как объём, и цена
 * становится ВЫХОДНО-ОГРАНИЧЕННОЙ — определяется числом пикселей, а не
 * содержимым сцены.
 *
 * Это утверждение до сих пор НЕ ПРОВЕРЯЛОСЬ ничем: во всех прогонах сетка была
 * равномерной. Здесь оно проверяется самым дешёвым способом, каким возможно, —
 * дерево строится, ячейки считаются, решение не запускается вовсе. Если закон
 * неверен, это видно на счётчике, и незачем платить за развёртку, чтобы узнать.
 *
 * ЧТО ИМЕННО ФАЛЬСИФИЦИРУЕТСЯ:
 *   1. число элементов ∝ 1/ε² — вдвое мельче пиксель даёт вчетверо больше;
 *   2. число элементов ∝ log₂(протяжённости), а НЕ ∝ объёму: удвоение сцены по
 *      каждой оси обязано дать ПРИБАВКУ, а не восьмикратный рост;
 *   3. выигрыш против равномерной сетки того же мелкого шага.
 *
 * Признак артефакта, за которым тут надо следить особо (первая подпись): если
 * число ячеек НЕ МЕНЯЕТСЯ при смене ε, значит правило не применяется вовсе.
 */

#include "octree.h"
#include "transport/ray3.h"
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ЕСТЬ ЛИ В УЗЛЕ ПОВЕРХНОСТЬ. Модель сцены здесь нарочно грубая — шесть стенок
 * куба и шар посередине, — потому что проверяется ЗАКОН СЧЁТА, а не геометрия.
 * Нужна она затем, что заполнять ячейками ВЕСЬ ОБЪЁМ и дробить только у
 * ПОВЕРХНОСТЕЙ — это два разных закона: `1/ε³` на октаву против `1/ε²`. Вся
 * ставка архитектуры на втором, и держится он на СХЛОПЫВАНИИ ПУСТОТЫ. */
static int has_surface(const int lo[3], int size, int world, double sc[3], double sr) {
  for (int a = 0; a < 3; a++) /* стенка куба проходит по узлу */
    if (lo[a] == 0 || lo[a] + size == world) return 1;
  double dmin2 = 0.0, dmax2 = 0.0;
  for (int a = 0; a < 3; a++) {
    double l = (double)lo[a], h = (double)(lo[a] + size);
    double dl = sc[a] - l, dh = h - sc[a];
    double far = dl > dh ? dl : dh;
    dmax2 += far * far;
    double near = 0.0;
    if (sc[a] < l)
      near = l - sc[a];
    else if (sc[a] > h)
      near = sc[a] - h;
    dmin2 += near * near;
  }
  return dmin2 <= sr * sr && dmax2 >= sr * sr; /* сфера пересекает коробку */
}

/* То же дробление, но пустой узел НЕ ДРОБИТСЯ ВОВСЕ: поверхности в нём нет,
 * поле сквозь него просто течёт, и мельчить незачем. */
static void lod_surf_rec(hz_octree *t, const int lo[3], int size, const double eye[3], double eps,
                         int minsize, int world, double sc[3], double sr, long *ncell) {
  if (!has_surface(lo, size, world, sc, sr)) {
    int hi[3] = {lo[0] + size, lo[1] + size, lo[2] + size};
    hz_oct_set_box(t, lo, hi, 1.0);
    (*ncell)++;
    return;
  }
  double near2 = 0.0;
  for (int a = 0; a < 3; a++) {
    double l = (double)lo[a], h = (double)(lo[a] + size), d = 0.0;
    if (eye[a] < l)
      d = l - eye[a];
    else if (eye[a] > h)
      d = eye[a] - h;
    near2 += d * d;
  }
  if (size <= minsize || (double)size <= eps * sqrt(near2)) {
    int hi[3] = {lo[0] + size, lo[1] + size, lo[2] + size};
    hz_oct_set_box(t, lo, hi, 1.0);
    (*ncell)++;
    return;
  }
  int h2 = size / 2;
  for (int k = 0; k < 8; k++) {
    int c[3] = {lo[0] + ((k & 1) ? h2 : 0), lo[1] + ((k & 2) ? h2 : 0), lo[2] + ((k & 4) ? h2 : 0)};
    lod_surf_rec(t, c, h2, eye, eps, minsize, world, sc, sr, ncell);
  }
}

/* Рекурсивный построитель: узел дробится, пока его размер БОЛЬШЕ, чем `εR` в
 * ближайшей его точке. Берётся именно БЛИЖНЯЯ точка, а не центр: правило есть
 * потолок на размер элемента, и нарушать его хоть где-то внутри узла нельзя. */
static void lod_rec(hz_octree *t, const int lo[3], int size, const double eye[3], double eps,
                    int minsize, long *ncell) {
  double near2 = 0.0;
  for (int a = 0; a < 3; a++) {
    double l = (double)lo[a], h = (double)(lo[a] + size), d = 0.0;
    if (eye[a] < l)
      d = l - eye[a];
    else if (eye[a] > h)
      d = eye[a] - h;
    near2 += d * d;
  }
  double want = eps * sqrt(near2);
  if (size <= minsize || (double)size <= want) {
    int hi[3] = {lo[0] + size, lo[1] + size, lo[2] + size};
    hz_oct_set_box(t, lo, hi, 1.0);
    (*ncell)++;
    return;
  }
  int h = size / 2;
  for (int k = 0; k < 8; k++) {
    int c[3] = {lo[0] + ((k & 1) ? h : 0), lo[1] + ((k & 2) ? h : 0), lo[2] + ((k & 4) ? h : 0)};
    lod_rec(t, c, h, eye, eps, minsize, ncell);
  }
}

/* Число ячеек при заданных ε и глубине дерева. minsize — предел дробления,
 * он ЕСТЬ ЧАСТЬ ЗАМЕРА: без него дерево у самой камеры уходит в бесконечность,
 * и это не свойство закона, а свойство того, что камера стоит внутри сцены. */
static long count_lod(int log2size, const double eye[3], double eps, int minsize) {
  hz_octree t;
  if (hz_oct_init(&t, log2size, 0.0)) return -1;
  long n = 0;
  int lo[3] = {0, 0, 0};
  lod_rec(&t, lo, 1 << log2size, eye, eps, minsize, &n);
  hz_oct_free(&t);
  return n;
}

/* ---------------- К52: СКОЛЬКО ЯЧЕЕК ВДОЛЬ ЛУЧА, И ЭТО РАЗНЫЕ ЗАКОНЫ ------
 *
 * План обосновывал формулу «A только вместе» тем, что LOD без проекции
 * превращает марш из 23 ячеек в 7500, и число это выведено как
 * `ln(R_max/R_min)/ε`. Но такой счёт предполагает, что ячейками заполнен ВЕСЬ
 * ОБЪЁМ, а К49 измерила обратное: закон ПОВЕРХНОСТНЫЙ, и держится он на
 * СХЛОПЫВАНИИ ПУСТОТЫ. Луч большую часть пути идёт по пустым КРУПНЫМ узлам.
 *
 * Здесь марш идёт по НАСТОЯЩЕМУ маршу (`tr3_march` без геометрии), а не по
 * переписанному обходу: переписанный проверял бы себя. Считаются два дерева при
 * ОДНОМ `ε` — со схлопнутой пустотой и заполняющее объём, — и вся разница между
 * ними есть цена схлопывания. */
static void ray_steps(const hz_octree *t, int world, const double eye[3], int nray, double *med,
                      double *mx, double *mean) {
  hz_frame fr = {{0, 0, 0}, {1, 1, 1}};
  tr3_scene sc;
  memset(&sc, 0, sizeof sc);
  sc.tree = t;
  sc.fr = fr;
  int *hist = calloc((size_t)nray, sizeof(int));
  if (hist == NULL) return;
  long sum = 0;
  int worst = 0, nn = 0;
  /* ДЕТЕРМИНИРОВАННЫЙ веер направлений: спираль Фибоначчи по сфере. Случайности
   * здесь не нужно, а воспроизводимость нужна. */
  const double ga = 2.39996322972865332;
  for (int i = 0; i < nray; i++) {
    double z = 1.0 - 2.0 * ((double)i + 0.5) / (double)nray;
    double r = sqrt(1.0 - z * z), a = ga * (double)i;
    double d[3] = {r * cos(a), r * sin(a), z};
    tr3_hit h;
    if (tr3_march(&sc, eye, d, -1.0, &h) != 0) continue;
    hist[nn++] = h.nsteps;
    sum += h.nsteps;
    if (h.nsteps > worst) worst = h.nsteps;
  }
  (void)world;
  for (int i = 1; i < nn; i++) { /* сортировка вставками: nray мал */
    int v = hist[i], j = i - 1;
    while (j >= 0 && hist[j] > v) {
      hist[j + 1] = hist[j];
      j--;
    }
    hist[j + 1] = v;
  }
  *med = nn > 0 ? (double)hist[nn / 2] : 0.0;
  *mx = (double)worst;
  *mean = nn > 0 ? (double)sum / (double)nn : 0.0;
  free(hist);
}

int main(void) {
  printf("=== ПРЕДСКАЗАНИЯ (до единого результата) ===\n");
  printf("  П1 число элементов ∝ 1/ε²: вдвое мельче пиксель -> ВЧЕТВЕРО больше\n");
  printf("  П2 удвоение сцены по каждой оси даёт ПРИБАВКУ (логарифм), а не ×8\n");
  printf("  П3 выигрыш против равномерной сетки того же шага — многие порядки\n");
  printf("  ПРИЗНАК АРТЕФАКТА: если число НЕ меняется при смене ε — правило не\n");
  printf("  применяется вовсе\n");
  printf("  П5 ячеек ВДОЛЬ ЛУЧА при схлопнутой пустоте — ДЕСЯТКИ, и растут\n");
  printf("     логарифмически с размером сцены\n");
  printf("  П6 у дерева, дробящего ВЕСЬ ОБЪЁМ, — ТЫСЯЧИ, порядка ln(Rmax/Rmin)/ε;\n");
  printf("     значит 7500 из плана принадлежит ОБЪЁМНОМУ закону (К52)\n");
  printf("=== РЕЗУЛЬТАТЫ ===\n");

  /* П1: ЗАВИСИМОСТЬ ОТ ε ПРИ ФИКСИРОВАННОЙ СЦЕНЕ */
  const int L = 9; /* куб 512 единиц */
  double eye[3] = {256.0, 8.0, 256.0};
  printf("\n  [П1] куб %d, камера внутри, предел дробления 1:\n", 1 << L);
  long prev = 0;
  for (int k = 0; k < 5; k++) {
    double eps = 0.4 / (double)(1 << k);
    long n = count_lod(L, eye, eps, 1);
    printf("     ε = %.4f : %9ld ячеек", eps, n);
    if (prev > 0) printf("   отношение к предыдущему %.2f (ожидается 4)", (double)n / (double)prev);
    printf("\n");
    prev = n;
  }

  /* П2: ЗАВИСИМОСТЬ ОТ ПРОТЯЖЁННОСТИ СЦЕНЫ ПРИ ФИКСИРОВАННОМ ε.
   * Камера стоит в одном и том же месте у пола, сцена растёт вокруг. Если закон
   * верен, прибавка на каждое удвоение ПОСТОЯННА — это и есть логарифм. */
  printf("\n  [П2] ε = 0.05 фиксирован, сцена растёт вдвое по каждой оси:\n");
  long p2 = 0;
  for (int lg = 5; lg <= 11; lg++) {
    double e2[3] = {(double)(1 << lg) * 0.5, 8.0, (double)(1 << lg) * 0.5};
    long n = count_lod(lg, e2, 0.05, 1);
    printf("     куб %5d : %9ld ячеек", 1 << lg, n);
    if (p2 > 0)
      printf("   прибавка %+8ld   (при росте по ОБЪЁМУ было бы ×8 = %ld)", n - p2, p2 * 8);
    printf("\n");
    p2 = n;
  }

  /* П3: ПРОТИВ РАВНОМЕРНОЙ СЕТКИ ТОГО ЖЕ МЕЛКОГО ШАГА */
  printf("\n  [П3] выигрыш против равномерной сетки с тем же мелким шагом:\n");
  for (int lg = 7; lg <= 10; lg++) {
    double e2[3] = {(double)(1 << lg) * 0.5, 8.0, (double)(1 << lg) * 0.5};
    long n = count_lod(lg, e2, 0.05, 1);
    double unif = pow((double)(1 << lg), 3.0);
    printf("     куб %5d : LOD %9ld против равномерной %.3e   выигрыш ×%.3e\n", 1 << lg, n, unif,
           unif / (double)n);
  }
  /* П4: НАСТОЯЩИЙ ЗАКОН — ДРОБИТЬ ТОЛЬКО У ПОВЕРХНОСТЕЙ, ПУСТОТУ СХЛОПЫВАТЬ.
   * Здесь ожидается уже `1/ε²`, то есть ×4 на уполовинивание, а не ×8: элементы
   * живут на ДВУМЕРНОМ множестве, а не в объёме. Разница между этим прогоном и
   * [П1] и есть цена того, что пустота не схлопнута. */
  printf("\n  [П4] дробление ТОЛЬКО у поверхностей (стенки куба + шар), пустота схлопнута:\n");
  {
    const int LG = 9, WORLD = 1 << LG;
    double sc[3] = {(double)WORLD * 0.5, (double)WORLD * 0.5, (double)WORLD * 0.5};
    double sr = (double)WORLD * 0.15;
    double e3[3] = {(double)WORLD * 0.5, 8.0, (double)WORLD * 0.5};
    long pv = 0;
    for (int k = 0; k < 5; k++) {
      double eps = 0.4 / (double)(1 << k);
      hz_octree t;
      if (hz_oct_init(&t, LG, 0.0)) break;
      long n = 0;
      int lo[3] = {0, 0, 0};
      lod_surf_rec(&t, lo, WORLD, e3, eps, 1, WORLD, sc, sr, &n);
      hz_oct_free(&t);
      printf("     ε = %.4f : %9ld ячеек", eps, n);
      if (pv > 0)
        printf("   отношение %.2f (ожидается 4 у поверхностного закона, 8 у объёмного)",
               (double)n / (double)pv);
      printf("\n");
      pv = n;
    }
  }

  /* [П5/П6] К52: СКОЛЬКО ЯЧЕЕК ВДОЛЬ ЛУЧА У ДВУХ ЗАКОНОВ ПРИ ОДНОМ `ε`. */
  printf("\n  [П5/П6] ячеек ВДОЛЬ ЛУЧА, 4096 лучей из камеры (К52):\n");
  {
    const int LG = 9, WORLD = 1 << LG;
    double sc[3] = {(double)WORLD * 0.5, (double)WORLD * 0.5, (double)WORLD * 0.5};
    double sr = (double)WORLD * 0.15;
    double e3[3] = {(double)WORLD * 0.5, 8.0, (double)WORLD * 0.5};
    for (int k = 0; k < 3; k++) {
      double eps = 0.1 / (double)(1 << k);
      double md, mx, mn;
      hz_octree t1;
      if (hz_oct_init(&t1, LG, 0.0)) break;
      long n1 = 0;
      int lo[3] = {0, 0, 0};
      lod_surf_rec(&t1, lo, WORLD, e3, eps, 1, WORLD, sc, sr, &n1);
      ray_steps(&t1, WORLD, e3, 4096, &md, &mx, &mn);
      printf("     ε = %.4f  ПУСТОТА СХЛОПНУТА: ячеек %8ld, вдоль луча медиана %.0f, "
             "среднее %.1f, макс %.0f\n",
             eps, n1, md, mn, mx);
      hz_oct_free(&t1);

      hz_octree t2;
      if (hz_oct_init(&t2, LG, 0.0)) break;
      long n2 = 0;
      lod_rec(&t2, lo, WORLD, e3, eps, 1, &n2);
      ray_steps(&t2, WORLD, e3, 4096, &md, &mx, &mn);
      printf("     ε = %.4f  ВЕСЬ ОБЪЁМ:       ячеек %8ld, вдоль луча медиана %.0f, "
             "среднее %.1f, макс %.0f\n",
             eps, n2, md, mn, mx);
      hz_oct_free(&t2);
    }
  }
  return 0;
}
