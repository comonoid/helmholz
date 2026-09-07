/* pfield — ИСТОЧНИК ЭРМИТОВА ПОЛЯ ИЗ МЕША (PLAN_ELEMENTS.md §363, §393; Р1, Р7).
 *
 * ЗАЧЕМ. Решение §363: представление сцены — эрмитово поле на восьмидереве, а не
 * индекс над треугольниками.
 *
 * ЗНАКА «ВНУТРИ/СНАРУЖИ» ЗДЕСЬ НЕТ ВОВСЕ, И ЭТО ЗАМЕР, А НЕ ВКУС (Р7, §391).
 * До 08-10 знак строился заливкой пустоты от границы, и она отказывала так:
 * пока ячейка крупна, заливка в зал не проникает, зал считается «телом», стены
 * рисуются, а предметы ВНУТРИ зала контраста не имеют и пропадают; как ячейка
 * мельчает, заливка протекает сквозь щели меша, зал становится «наружей», и
 * односторонняя стена, у которой «наружа» с обеих сторон, пропадает тоже.
 * Уровня, на котором сцена цела, не существовало ни одного. Причина глубже
 * щелей: у односторонней поверхности нулевой толщины ВНУТРЕННОСТИ НЕТ.
 *
 * ЧТО ВМЕСТО. Поверхность опознаётся РЕБРОМ С ПЕРЕСЕЧЕНИЕМ. Сцена строит и
 * ОТДАЁТ (Р1) две вещи: занятость (пирамида битовых карт, по ней идёт спуск) и
 * таблицу эрмитовых рёбер (точка на ребре плюс нормаль грани). Обратных вызовов
 * «знак угла» и «пересечение ребра» больше нет; остался один — «есть ли
 * геометрия в коробке», и он про коробку, а не про поле.
 *
 * ЧЕГО ЭТОТ СТЕНД НЕ ПОКАЗЫВАЕТ. Ошибка поверхности ОДНОСТОРОННЯЯ (А598): ловит
 * потерю геометрии, не ловит появление лишней, и популяция у неё — ВЫХОД.
 * Приёмка потери считается отдельно, и её популяция ВХОД: занятые ячейки,
 * поклеточно сверенные с независимым эталоном `tools/poccref.c` (Ш0, §387).
 * Ни отскока, ни тени, ни цвета материала здесь нет: картинка проверяет
 * ГЕОМЕТРИЮ поля и только её.
 */
#include "cut/dc.h"
#include "cut/dcslice.h"
#include "image.h"
#include "occmap.h"
#include "pclip.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include "transport/cam3.h"
#include "transport/cut3.h"
#include "transport/dirs3.h"
#include "transport/mesh3.h"
#include "transport/sweep3.h"
#include <math.h>
#include <omp.h>
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

static int cmp_i32(const void *x, const void *y) {
  int32_t a = *(const int32_t *)x, b = *(const int32_t *)y;
  return a < b ? -1 : (a > b ? 1 : 0);
}

static int cmp_f(const void *x, const void *y) {
  float a = *(const float *)x, b = *(const float *)y;
  return (a < b) ? -1 : ((a > b) ? 1 : 0);
}

static int cmp_dev699(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return x < y ? -1 : (x > y ? 1 : 0);
}

static int cmp_d(const void *x, const void *y) {
  double a = *(const double *)x, b = *(const double *)y;
  return (a < b) ? -1 : ((a > b) ? 1 : 0);
}

/* --- общая рама сцены ------------------------------------------------------ */

typedef struct {
  double org[3]; /* нижний угол кубической рамы */
  double h;      /* сторона мелкой ячейки */
  int32_t n;     /* ячеек по стороне = 1 << lev */
  int lev;
} frame;

/* --- 1. ЗАНЯТОСТЬ: ПИРАМИДА БИТОВЫХ КАРТ (Р7, §393) ------------------------ */

/* ЗАЧЕМ ПИРАМИДА, А НЕ ДЕРЕВО. Спуску нужен один ответ — «есть ли геометрия в
 * коробке `[lo, lo+size)`», — и коробки эти ВЫРОВНЕНЫ по своему размеру. Для
 * выровненной коробки ответ есть ОДИН БИТ на соответствующем уровне пирамиды,
 * то есть чтение по адресу, а не поиск по дереву. Второе дерево (§381, 68.9 МБ)
 * заводилось ровно затем, чтобы отвечать на этот вопрос поиском; здесь его нет
 * и завестись оно не может — массив запросов по координате не принимает.
 * Уровень L при L=9 стоит 16.8 МБ, вся пирамида 8/7 от этого.
 *
 * ОПРЕДЕЛЕНИЕ ЗАНЯТОСТИ ТО ЖЕ, ЧТО У ЭТАЛОНА Ш0 (§385, `tools/poccref.c`):
 * кусок треугольника, отсечённый коробкой ячейки, имеет положительную площадь.
 * Формулы коробки и диапазона ячеек скопированы дословно — разойдись они,
 * сверка с эталоном мерила бы разницу формул. */
typedef struct {
  unsigned char *b[HZ_DC_MAX_LOG2SIZE + 1];
  int lev;
} opyr;

static void opyr_free(opyr *P) {
  for (int l = 0; l <= P->lev; l++)
    free(P->b[l]);
}

static void tri_verts(const hz_objmesh *m, int32_t t, const double **A, const double **B,
                      const double **C) {
  *A = m->v + 3 * (size_t)m->f[3 * (size_t)t + 0];
  *B = m->v + 3 * (size_t)m->f[3 * (size_t)t + 1];
  *C = m->v + 3 * (size_t)m->f[3 * (size_t)t + 2];
}

static void tri_cells(const hz_objmesh *m, const frame *fr, int32_t t, int64_t i0[3],
                      int64_t i1[3]) {
  for (int c = 0; c < 3; c++) {
    double a = 1e300, b = -1e300;
    for (int v = 0; v < 3; v++) {
      double x = m->v[3 * (size_t)m->f[3 * (size_t)t + (size_t)v] + (size_t)c];
      if (x < a) a = x;
      if (x > b) b = x;
    }
    /* РАСШИРЕНИЕ НА ОДНУ ЯЧЕЙКУ В КАЖДУЮ СТОРОНУ (§443). Прежний диапазон
     * назывался ТОЧНЫМ отсевом, и на зале это подтверждалось; на ОСЕВОЙ
     * геометрии он терял 52 % занятых ячеек: габарит треугольника вырожден по
     * оси, диапазон сжимается в один слой, а плоскость касается двух. Точность
     * обеспечивает точное отсечение, которое всё равно вызывается; габарит
     * обязан лишь НЕ ТЕРЯТЬ, и теперь он консервативен. */
    i0[c] = (int64_t)floor((a - fr->org[c]) / fr->h) - 1;
    i1[c] = (int64_t)floor((b - fr->org[c]) / fr->h) + 1;
    if (i0[c] < 0) i0[c] = 0;
    if (i1[c] >= fr->n) i1[c] = fr->n - 1;
  }
}

static void cell_box(const frame *fr, int64_t x, int64_t y, int64_t z, double cl[3], double ch[3]) {
  cl[0] = fr->org[0] + (double)x * fr->h;
  cl[1] = fr->org[1] + (double)y * fr->h;
  cl[2] = fr->org[2] + (double)z * fr->h;
  for (int c = 0; c < 3; c++)
    ch[c] = cl[c] + fr->h;
}

static int64_t opyr_build(opyr *P, const hz_objmesh *m, const frame *fr) {
  P->lev = fr->lev;
  for (int l = 0; l <= P->lev; l++) {
    int32_t n = (int32_t)1 << l;
    P->b[l] = calloc(hz_occ_bytes((size_t)n * (size_t)n * (size_t)n), 1);
    if (P->b[l] == NULL) exit(1);
  }
  int64_t nocc = 0;
  unsigned char *bt = P->b[P->lev];
  for (int32_t t = 0; t < m->nt; t++) {
    int64_t i0[3] = {0, 0, 0}, i1[3] = {0, 0, 0};
    tri_cells(m, fr, t, i0, i1);
    const double *A, *B, *C;
    tri_verts(m, t, &A, &B, &C);
    for (int64_t z = i0[2]; z <= i1[2]; z++)
      for (int64_t y = i0[1]; y <= i1[1]; y++)
        for (int64_t x = i0[0]; x <= i1[0]; x++) {
          double cl[3], ch[3];
          cell_box(fr, x, y, z, cl, ch);
          hz_pclip_poly Q;
          if (hz_pclip_tri(A, B, C, cl, ch, &Q) < 3) continue;
          if (!(hz_pclip_area(&Q) > 0.0)) continue;
          size_t ci = hz_occ_index(fr->n, x, y, z);
          if (!hz_occ_get(bt, ci)) {
            hz_occ_set(bt, ci);
            nocc++;
          }
        }
  }
  /* Подъём: бит родителя = ИЛИ по восьми детям. Тождество, а не приближение. */
  for (int l = P->lev; l > 0; l--) {
    int32_t n = (int32_t)1 << l, n2 = n / 2;
    for (int32_t z = 0; z < n; z++)
      for (int32_t y = 0; y < n; y++)
        for (int32_t x = 0; x < n; x++)
          if (hz_occ_get(P->b[l], hz_occ_index(n, x, y, z)))
            hz_occ_set(P->b[l - 1], hz_occ_index(n2, x / 2, y / 2, z / 2));
  }
  return nocc;
}

/* ЕДИНСТВЕННЫЙ ВОПРОС, КОТОРЫЙ ОБЩИЙ СЛОЙ ЗАДАЁТ СЦЕНЕ. Ответ ТОЧНЫЙ, а не
 * консервативный, — это сильнее контракта и потому законно. */
static int u_occ(void *ctx, const int32_t lo[3], int32_t size) {
  const opyr *P = (const opyr *)ctx;
  int sh = 0;
  for (int32_t s = size; s > 1; s >>= 1)
    sh++;
  int lvl = P->lev - sh;
  int32_t n = (int32_t)1 << lvl;
  return hz_occ_get(P->b[lvl], hz_occ_index(n, lo[0] >> sh, lo[1] >> sh, lo[2] >> sh));
}

/* --- 2. ЭРМИТОВЫ РЁБРА: СЦЕНА ИХ СТРОИТ И ОТДАЁТ (Р1) ---------------------- */

/* Обратных вызовов `hz_dc_cross` больше нет. Сцена проходит по треугольникам и
 * заносит в таблицу пересечение с каждым сеточным ребром, которое треугольник
 * задевает; сведение по ближайшему — `hz_htab_reduce_min`. Отсюда же исчезает
 * четырёхкратный опрос ребра (§379: 1.1 млн вызовов на 277 тыс. рёбер): ребро
 * рождается от ТРЕУГОЛЬНИКА, а не от каждой из четырёх своих ячеек.
 *
 * НОРМАЛЬ — ГРАНИ (А656, гранёный материал). Вершинные нормали модели в поле не
 * попадают и здесь: гладкий режим Р2 требует бита материала, а материалов в
 * структуре ещё нет. Говорится прямо, чтобы пропуск не читался сделанным. */
/* СЕТКА, СОВПАВШАЯ С ГЕОМЕТРИЕЙ, СЧИТАЕТ ПОВЕРХНОСТЬ ДВАЖДЫ (§451). Осевая стена,
 * легшая в плоскость сеточных вершин, пересекается ОБОИМИ смежными рёбрами: снизу
 * с `t = 1−δ`, сверху с `t = +δ'`, и оба СТРОГО ВНУТРИ отрезка. Полуоткрытый
 * отрезок `t ∈ [0, 1)` этого НЕ ЛЕЧИТ — проверено замером: записей осталось те же
 * `8 649`, потому что `t` не равен единице ни разу, он равен `1 − 1e-16`.
 * Лечится ЕДИНСТВЕННЫМ местом — положением рамы (`hz_gridalign` ниже).
 * Ключ `gridalign` возвращает прежнюю раму: это НЕГАТИВНЫЙ КОНТРОЛЬ, и отношение
 * площади обязано вернуться к `1.5000` ТОЧНО. */
/* НЕГАТИВНЫЙ КОНТРОЛЬ §581: один поток принудительно — все времена обязаны
 * вернуться к однопоточным. Объявлен рано, потому что читается прагмами. */
static int g_omp1 = 0;

static int hz_gridalign = 0;
/* НЕГАТИВНЫЙ КОНТРОЛЬ §463 (`flipall`): обратить порядок вершин У ВСЕХ выданных
 * треугольников. Доля вывернутых и по числу, и ПО ПЛОЩАДИ обязана стать почти
 * единицей, а гистограмма косинуса — зеркально перевернуться. Не перевернётся —
 * измеритель меряет не то, что называет. */
static int hz_flipall = 0;
/* НЕГАТИВНЫЙ КОНТРОЛЬ Ш10 (`longdiag`): резать квад по ДЛИННОЙ диагонали. Доля
 * вывернутых по площади обязана ВЫРАСТИ; не шелохнётся — выбор диагонали ни на
 * что не влияет. */
static int hz_longdiag = 0;

/* Дерево по ЗАНЯТОСТИ для стыка с переносом (Ш13): пустая коробка остаётся
 * крупным листом, занятая дробится до предела. Тот же предикат, что у
 * `hz_dc_shape_occ`, поэтому сетка переноса и дерево DC согласованы по
 * построению, а не по совпадению. */
static void oct_by_occ(hz_octree *t, const opyr *P, int lev, int32_t x, int32_t y, int32_t z,
                       int32_t size) {
  int lv = lev;
  int32_t s = size;
  while (s > 1) {
    lv--;
    s /= 2;
  }
  int occ = hz_occ_get(P->b[lv], hz_occ_index((int32_t)1 << lv, x >> (lev - lv), y >> (lev - lv),
                                              z >> (lev - lv))) != 0;
  if (!occ || size == 1) {
    int lo[3] = {(int)x, (int)y, (int)z},
        hi[3] = {(int)(x + size), (int)(y + size), (int)(z + size)};
    hz_oct_set_box(t, lo, hi, 1.0);
    return;
  }
  int32_t h = size / 2;
  for (int k = 0; k < 8; k++)
    oct_by_occ(t, P, lev, x + ((k & 1) ? h : 0), y + ((k & 2) ? h : 0), z + ((k & 4) ? h : 0), h);
}
/* П2 §467: пометка «треугольник родом из СЛОЖЕННОГО квада». Объявлена здесь,
 * потому что ставит её обход, а читает укладка треугольника. */
static int g_curfold = 0;

static int seg_tri(const double *A, const double *B, const double *C, const double P0[3], int axis,
                   double h, double *tt, double nrm[3]) {
  double e1[3], e2[3], pv[3], tv[3], qv[3], dir[3] = {0.0, 0.0, 0.0};
  dir[axis] = h;
  for (int c = 0; c < 3; c++) {
    e1[c] = B[c] - A[c];
    e2[c] = C[c] - A[c];
  }
  pv[0] = dir[1] * e2[2] - dir[2] * e2[1];
  pv[1] = dir[2] * e2[0] - dir[0] * e2[2];
  pv[2] = dir[0] * e2[1] - dir[1] * e2[0];
  double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
  if (det > -1e-15 && det < 1e-15) return 0;
  double inv = 1.0 / det;
  for (int c = 0; c < 3; c++)
    tv[c] = P0[c] - A[c];
  double uu = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
  if (uu < 0.0 || uu > 1.0) return 0;
  qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
  qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
  qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
  double vv = (dir[0] * qv[0] + dir[1] * qv[1] + dir[2] * qv[2]) * inv;
  if (vv < 0.0 || uu + vv > 1.0) return 0;
  double t = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
  if (t < 0.0 || t > 1.0) return 0;
  double cr[3];
  cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
  cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
  cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
  double l2 = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
  if (!(l2 > 0.0)) return 0;
  *tt = t;
  for (int c = 0; c < 3; c++)
    nrm[c] = cr[c] / l2;
  return 1;
}

static int64_t edges_build(hz_htab *ht, const hz_objmesh *m, const frame *fr, int nonrm) {
  int64_t nhit = 0;
  for (int32_t t = 0; t < m->nt; t++) {
    int64_t i0[3] = {0, 0, 0}, i1[3] = {0, 0, 0};
    tri_cells(m, fr, t, i0, i1);
    const double *A, *B, *C;
    tri_verts(m, t, &A, &B, &C);
    /* Решётка на единицу шире коробки в обе стороны: ребро `[p, p+1]` вдоль оси
     * может задеть треугольник и при `p = i0-1`. Лишние узлы отсеет сам тест. */
    for (int64_t z = i0[2] - 1; z <= i1[2] + 1; z++) {
      if (z < 0 || z > fr->n) continue;
      for (int64_t y = i0[1] - 1; y <= i1[1] + 1; y++) {
        if (y < 0 || y > fr->n) continue;
        for (int64_t x = i0[0] - 1; x <= i1[0] + 1; x++) {
          if (x < 0 || x > fr->n) continue;
          double P0[3];
          cell_box(fr, x, y, z, P0, (double[3]){0, 0, 0});
          for (int a = 0; a < 3; a++) {
            int32_t p[3] = {(int32_t)x, (int32_t)y, (int32_t)z};
            if (p[a] >= fr->n) continue;
            double tv, nr[3];
            if (!seg_tri(A, B, C, P0, a, fr->h, &tv, nr)) continue;
            if (nonrm) { /* НЕГАТИВНЫЙ КОНТРОЛЬ §393: нормаль осевая, как заплата */
              nr[0] = nr[1] = nr[2] = 0.0;
              nr[a] = 1.0;
            }
            if (hz_htab_add(ht, a, p, tv, nr, 0) != 0) exit(1);
            nhit++;
          }
        }
      }
    }
  }
  return nhit;
}

/* --- 2а. ИНДЕКС ЗАНЯТЫХ ЯЧЕЕК: РАНГ ПО БИТОВОЙ КАРТЕ ----------------------- */

/* Списки треугольников нужны ДВУМ потребителям: ошибке поверхности и разбору
 * потерянных. Держать их на все `n³` ячеек нельзя (при L=9 это 134 млн int32),
 * поэтому индекс идёт по РАНГУ в битовой карте занятости: занятых 817 тыс., и
 * массивы считаются от них. Ранг — сумма по восьмибайтовым блокам плюс popcount
 * хвоста; это тот же приём, что prefix-sum, только без массива на ячейку. */
typedef struct {
  const unsigned char *b;
  int32_t *pre; /* pre[i] — единиц в байтах [0, 8i) */
  int64_t nw;
  int32_t n;
} rankmap;

static int64_t rm_rank(const rankmap *R, size_t bit) {
  size_t byte = bit >> 3;
  int64_t r = R->pre[byte >> 3];
  for (size_t k = (byte & ~(size_t)7); k < byte; k++)
    r += __builtin_popcount(R->b[k]);
  r += __builtin_popcount((unsigned)(R->b[byte] & (unsigned char)((1u << (bit & 7u)) - 1u)));
  return r;
}

static void rm_build(rankmap *R, const unsigned char *b, int32_t n) {
  R->b = b;
  R->n = n;
  size_t nb = hz_occ_bytes((size_t)n * (size_t)n * (size_t)n);
  R->nw = (int64_t)((nb + 7) / 8) + 1;
  R->pre = malloc((size_t)R->nw * sizeof *R->pre);
  if (R->pre == NULL) exit(1);
  int32_t acc = 0;
  for (int64_t w = 0; w < R->nw; w++) {
    R->pre[w] = acc;
    for (size_t k = (size_t)w * 8; k < (size_t)(w + 1) * 8 && k < nb; k++)
      acc += __builtin_popcount(b[k]);
  }
}

typedef struct {
  rankmap rm;
  int32_t *cnt; /* nocc+1, префиксные суммы */
  int32_t *lst;
  const frame *fr;
} celltris;

static int32_t ct_list(void *ctx, const int32_t cell[3], const int32_t **out) {
  celltris *T = (celltris *)ctx;
  int32_t n = T->fr->n;
  if (cell[0] < 0 || cell[1] < 0 || cell[2] < 0 || cell[0] >= n || cell[1] >= n || cell[2] >= n)
    return 0;
  size_t ci = hz_occ_index(n, cell[0], cell[1], cell[2]);
  if (!hz_occ_get(T->rm.b, ci)) return 0;
  int64_t r = rm_rank(&T->rm, ci);
  *out = T->lst + T->cnt[r];
  return T->cnt[r + 1] - T->cnt[r];
}

static void ct_build(celltris *T, const hz_objmesh *m, const frame *fr, const unsigned char *occ,
                     int64_t nocc) {
  T->fr = fr;
  rm_build(&T->rm, occ, fr->n);
  T->cnt = calloc((size_t)nocc + 2, sizeof *T->cnt);
  if (T->cnt == NULL) exit(1);
  for (int pass = 0; pass < 2; pass++) {
    for (int32_t t = 0; t < m->nt; t++) {
      int64_t i0[3] = {0, 0, 0}, i1[3] = {0, 0, 0};
      tri_cells(m, fr, t, i0, i1);
      const double *A, *B, *C;
      tri_verts(m, t, &A, &B, &C);
      for (int64_t z = i0[2]; z <= i1[2]; z++)
        for (int64_t y = i0[1]; y <= i1[1]; y++)
          for (int64_t x = i0[0]; x <= i1[0]; x++) {
            size_t ci = hz_occ_index(fr->n, x, y, z);
            if (!hz_occ_get(occ, ci)) continue;
            double cl[3], ch[3];
            cell_box(fr, x, y, z, cl, ch);
            hz_pclip_poly Q;
            if (hz_pclip_tri(A, B, C, cl, ch, &Q) < 3) continue;
            if (!(hz_pclip_area(&Q) > 0.0)) continue;
            int64_t r = rm_rank(&T->rm, ci);
            if (pass == 0)
              T->cnt[r + 1]++;
            else
              T->lst[T->cnt[r]++] = t;
          }
    }
    if (pass == 0) {
      for (int64_t i = 0; i < nocc; i++)
        T->cnt[i + 1] += T->cnt[i];
      T->lst = malloc((size_t)(T->cnt[nocc] > 0 ? T->cnt[nocc] : 1) * sizeof *T->lst);
      if (T->lst == NULL) exit(1);
    }
  }
  /* Второй проход сдвинул cnt на единицу вперёд — вернуть на место. */
  for (int64_t i = nocc; i > 0; i--)
    T->cnt[i] = T->cnt[i - 1];
  T->cnt[0] = 0;
}

static void ct_free(celltris *T) {
  free(T->rm.pre);
  free(T->cnt);
  free(T->lst);
}
/* --- 3. ОШИБКА ПОВЕРХНОСТИ -------------------------------------------------- */

/* Список треугольников ячейки берётся из `celltris` — ранга по битовой карте
 * занятости. Источник один, второго нет, и сверять их не нужно: прежние
 * `dense_tris`/`sparse_tris` существовали затем, чтобы сличать плотный путь с
 * ленивым, а ленивого больше нет. */
typedef int32_t (*trilist_fn)(void *ctx, const int32_t cell[3], const int32_t **out);

typedef struct {
  int32_t lo[3], size;
  double vx[3];
} vrec;

static void collect_verts(const hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size,
                          vrec *out, int64_t cap, int64_t *n, int leaves) {
  /* `leaves` РАЗДЕЛЯЕТ ДВЕ РАЗНЫЕ ВЕЛИЧИНЫ, которые до сих пор считались одной и
   * той же (найдено замером среза, §399): вершина есть У КАЖДОГО узла дерева,
   * включая внутренние, а поверхность при полной глубине строится по ЛИСТЬЯМ.
   * «Вершин выдано» без этого различения завышает счёт на все внутренние узлы —
   * при `L = 6` это `13 554` против `10 624`, то есть `27 %`. */
  if ((t->nd[ni].flags & HZ_DC_HASVERT) && (!leaves || t->nd[ni].child0 < 0)) {
    /* ЁМКОСТЬ ПЕРЕДАЁТСЯ ЯВНО, а не подразумевается из первого прохода. Проход
     * первый считает, второй пишет, и связь между ними держалась ТОЛЬКО тем, что
     * дерево между вызовами не менялось. Здесь она стала проверяемой: анализатор
     * эту связь через рекурсию теряет (класс из CLAUDE.md), а расхождение двух
     * проходов было бы настоящим дефектом, и вызывающий его ловит сверкой
     * `k == nv`. */
    if (out != NULL && *n < cap) {
      for (int a = 0; a < 3; a++) {
        out[*n].lo[a] = lo[a];
        out[*n].vx[a] = (double)lo[a] + (double)t->nd[ni].vq[a] * (double)size / 65535.0;
      }
      out[*n].size = size;
    }
    (*n)++;
  }
  if (t->nd[ni].child0 < 0) return;
  int32_t half = size / 2;
  for (int i = 0; i < 8; i++) {
    int32_t clo[3];
    for (int a = 0; a < 3; a++)
      clo[a] = lo[a] + (((i >> a) & 1) ? half : 0);
    collect_verts(t, t->nd[ni].child0 + i, clo, half, out, cap, n, leaves);
  }
}

/* §654: СКОЛЬКО ЛИСТОВ ПОВЕРХНОСТИ ВИДНО В МАСКЕ РЁБЕР ЯЧЕЙКИ.
 *
 * ЗАЧЕМ. Вопрос «стоит ли менять ядро на Cubical Marching Squares» (§653)
 * упирается в одно: много ли у нас ячеек, несущих БОЛЬШЕ ОДНОГО куска
 * поверхности. Dual contouring кладёт в ячейку ОДНУ вершину и такую ячейку
 * представить не может; CMS выделяет петли и может. Счётчик `nmulti` на этот
 * вопрос не отвечает: он про маску УГЛОВ, а в беззнаковом пути её нет (§654).
 *
 * КАК. Ровно первым шагом CMS: на каждой из шести граней пересечённые рёбра
 * соединяются попарно (двумерный marching squares), после чего число связных
 * компонент среди пересечённых рёбер и есть число петель. Соединение делается
 * объединением множеств по двенадцати рёбрам — целочисленно, без геометрии.
 *
 * ЧЕГО ЭТА ВЕЛИЧИНА НЕ ЗНАЕТ (А1062): два ПАРАЛЛЕЛЬНЫХ листа, пересекающих одни
 * и те же рёбра, дают ту же маску, что один лист. Значит это НИЖНЯЯ оценка
 * числа листов, и называть её надо так.
 *
 * НЕОДНОЗНАЧНАЯ ГРАНЬ (А1064): при ЧЕТЫРЁХ пересечениях на грани двумерный
 * marching squares выбирает соединение ПО ЗНАКАМ УГЛОВ, а знака у нас нет
 * вовсе. Такие грани СЧИТАЮТСЯ отдельно, а соединяются произвольно, но
 * ОДИНАКОВО у обоих соседей — то есть щели это не даёт, а правильность не
 * обещана. */
static int ecross_loops(uint16_t m, int *nambig) {
  int par[12];
  for (int i = 0; i < 12; i++)
    par[i] = i;
  for (int f = 0; f < 3; f++)
    for (int s = 0; s < 2; s++) {
      int e[4], ne = 0;
      for (int i = 0; i < 12; i++) {
        if (!((m >> i) & 1u)) continue;
        int a = i / 4, k = i % 4;
        if (a == f) continue; /* ребро идёт вдоль оси грани — на грани не лежит */
        int u = (a + 1) % 3;
        int off = (f == u) ? (k & 1) : ((k >> 1) & 1);
        if (off != s) continue;
        if (ne < 4) e[ne++] = i;
      }
      if (ne == 4 && nambig != NULL) (*nambig)++;
      for (int j = 0; j + 1 < ne; j += 2) {
        int x = e[j], y = e[j + 1];
        while (par[x] != x)
          x = par[x];
        while (par[y] != y)
          y = par[y];
        if (x != y) par[x] = y;
      }
    }
  int nloop = 0;
  for (int i = 0; i < 12; i++) {
    if (!((m >> i) & 1u)) continue;
    int x = i;
    while (par[x] != x)
      x = par[x];
    if (x == i) nloop++;
  }
  return nloop;
}

typedef struct {
  int64_t npop;     /* листья размера 1 с непустой маской рёбер — ЗНАМЕНАТЕЛЬ */
  int64_t nloop[4]; /* петель 0 / 1 / 2 / >=3 */
  int64_t nambigc;  /* ячеек, где есть неоднозначная грань */
  int64_t nambigf;  /* самих неоднозначных граней */
  int64_t nopp;     /* ячеек со ВСТРЕЧНЫМИ эрмитовыми нормалями (n_i·n_j < 0) */
} sheetstat;

static void sheet_rec(const hz_dctree *t, const hz_htab *ht, int32_t ni, const int32_t lo[3],
                      int32_t size, sheetstat *S) {
  if (t->nd[ni].child0 < 0) {
    /* ПОПУЛЯЦИЯ ОГРАНИЧЕНА ЯВНО (А1063): только лист размера 1. У крупного узла
     * `ecross` есть ИЛИ по детям (А709), и петли по нему считать бессмысленно. */
    if (size != 1) return;
    uint16_t m = t->nd[ni].ecross;
    if (m == 0) return;
    S->npop++;
    int amb = 0;
    int nl = ecross_loops(m, &amb);
    if (amb > 0) {
      S->nambigc++;
      S->nambigf += amb;
    }
    S->nloop[nl <= 0 ? 0 : (nl >= 3 ? 3 : nl)]++;
    /* ВСТРЕЧНЫЕ НОРМАЛИ — признак ИЗ ДРУГОЙ ПРИРОДЫ, чем петли: он видит два
     * листа даже там, где маска одна. Но он ШИРЕ: складка острее 90 градусов
     * даёт встречные нормали без всяких двух листов (П5 это и говорит). */
    double nr[12][3];
    int nn = 0;
    for (int i = 0; i < 12; i++) {
      int a = i / 4, k = i % 4;
      int u = (a + 1) % 3, v = (a + 2) % 3;
      int32_t p[3] = {lo[0], lo[1], lo[2]};
      p[u] += k & 1;
      p[v] += (k >> 1) & 1;
      const hz_hedge *e = hz_htab_find(ht, a, p);
      if (e == NULL || e->in_lo == HZ_HEDGE_ERASED) continue;
      for (int c = 0; c < 3; c++)
        nr[nn][c] = e->nrm[c];
      nn++;
    }
    for (int i = 0; i < nn; i++)
      for (int j = i + 1; j < nn; j++)
        if (nr[i][0] * nr[j][0] + nr[i][1] * nr[j][1] + nr[i][2] * nr[j][2] < 0.0) {
          S->nopp++;
          i = nn;
          break;
        }
    return;
  }
  int32_t half = size / 2;
  for (int i = 0; i < 8; i++) {
    int32_t clo[3];
    for (int a = 0; a < 3; a++)
      clo[a] = lo[a] + (((i >> a) & 1) ? half : 0);
    sheet_rec(t, ht, t->nd[ni].child0 + i, clo, half, S);
  }
}

static void surf_err(const hz_dctree *T, const frame *fr, const hz_objmesh *m, trilist_fn tl,
                     void *tctx, const char *tag) {
  int64_t nv = 0;
  collect_verts(T, 0, (int32_t[3]){0, 0, 0}, fr->n, NULL, 0, &nv, 0);
  if (nv == 0) return;
  double *er = malloc((size_t)nv * sizeof *er);
  vrec *vr = malloc((size_t)nv * sizeof *vr);
  if (er == NULL || vr == NULL) {
    free(er);
    free(vr);
    return;
  }
  int64_t k = 0;
  collect_verts(T, 0, (int32_t[3]){0, 0, 0}, fr->n, vr, nv, &k, 0);
  if (k != nv)
    printf("   %s ДВА ПРОХОДА РАЗОШЛИСЬ: %lld против %lld\n", tag, (long long)k, (long long)nv);
  int64_t ne = 0;
  for (int64_t i = 0; i < k; i++) {
    double V[3];
    for (int c = 0; c < 3; c++)
      V[c] = fr->org[c] + vr[i].vx[c] * fr->h;
    int32_t ic[3];
    int okc = 1;
    for (int c = 0; c < 3; c++) {
      double f = floor((V[c] - fr->org[c]) / fr->h);
      if (!(f >= 0.0) || !(f < (double)fr->n)) okc = 0;
      ic[c] = okc ? (int32_t)f : 0;
    }
    if (!okc) continue;
    double bd = 1e300;
    for (int dz = -1; dz <= 1; dz++)
      for (int dy = -1; dy <= 1; dy++)
        for (int dx = -1; dx <= 1; dx++) {
          int32_t c3[3] = {ic[0] + dx, ic[1] + dy, ic[2] + dz};
          if (c3[0] < 0 || c3[1] < 0 || c3[2] < 0 || c3[0] >= fr->n || c3[1] >= fr->n ||
              c3[2] >= fr->n)
            continue;
          const int32_t *ls = NULL;
          int32_t nl = tl(tctx, c3, &ls);
          for (int32_t q = 0; q < nl; q++) {
            const double *A, *B, *C;
            tri_verts(m, ls[q], &A, &B, &C);
            double e1[3], e2[3], cr[3];
            for (int c = 0; c < 3; c++) {
              e1[c] = B[c] - A[c];
              e2[c] = C[c] - A[c];
            }
            cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
            cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
            cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
            double l2 = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
            if (!(l2 > 0.0)) continue;
            double d =
                fabs((V[0] - A[0]) * cr[0] + (V[1] - A[1]) * cr[1] + (V[2] - A[2]) * cr[2]) / l2;
            if (d < bd) bd = d;
          }
        }
    if (bd < 1e299) er[ne++] = bd;
  }
  if (ne > 0) {
    qsort(er, (size_t)ne, sizeof *er, cmp_d);
    printf("   %s ОШИБКА ПОВЕРХНОСТИ (ОДНОСТОРОННЯЯ, А598): p50 %.6f p90 %.6f max %.6f м, по "
           "%lld вершинам\n",
           tag, er[ne / 2], er[(ne * 9) / 10], er[ne - 1], (long long)ne);
  }
  free(er);
  free(vr);
}

/* --- 3a. ПОКРЫТИЕ: ВТОРАЯ СТОРОНА ОШИБКИ (§380) ---------------------------- */

/* ЗАЧЕМ, И ПОЧЕМУ ЭТО НАДО БЫЛО ПИСАТЬ ПЕРВЫМ. У `surf_err` ПОПУЛЯЦИЯ ЕСТЬ
 * ВЫХОД: он берёт ВЫДАННЫЕ вершины и меряет расстояние до ближайшего
 * треугольника. Величина, чья популяция — выход, ОТСУТСТВИЯ НЕ ОБНАРУЖИВАЕТ
 * НИКОГДА: не выдали половину сцены — не с чего и мерить, медиана останется
 * нулём. Это третья подпись артефакта §4 («ошибка нормирована на сокращённый
 * результат»), и она была в правилах всё это время.
 *
 * ЗДЕСЬ ПОПУЛЯЦИЯ — ВХОД: перебираются ЗАНЯТЫЕ ячейки, то есть те, в которых
 * ЕСТЬ треугольник, и спрашивается, выдало ли поле вершину в этой ячейке или
 * хотя бы в одной из 26 соседних. Доля ячеек, где не выдало, И ЕСТЬ ДОЛЯ СЦЕНЫ,
 * ПОТЕРЯННОЙ ПОЛЕМ. Она не может «сесть на пол» от того, что мы построили
 * меньше: чем меньше построено, тем она ХУЖЕ. */
/* Вершина именно В ЭТОЙ ячейке, то есть у ЛИСТА. Засчитывать вершину предка
 * нельзя: у корня она есть всегда, и проверка стала бы тождеством — что и
 * случилось при первом прогоне (0.0 % на всех уровнях, §380). */
static int dc_vert_at(const hz_dctree *t, const int32_t cell[3]) {
  int32_t ni = 0, size = (int32_t)1 << t->log2size, lo[3] = {0, 0, 0};
  for (;;) {
    if (t->nd[ni].child0 < 0) return (t->nd[ni].flags & HZ_DC_HASVERT) ? 1 : 0;
    int32_t half = size / 2;
    int bit = 0;
    for (int a = 0; a < 3; a++)
      if (cell[a] >= lo[a] + half) {
        bit |= 1 << a;
        lo[a] += half;
      }
    ni = t->nd[ni].child0 + bit;
    size = half;
  }
}

/* РАЗБОР ПОТЕРЯННЫХ (§393, А705). ДВЕ величины, а не одна, и вторая строже:
 *   ОКНО 27  вершина в ячейке ИЛИ в любой из 26 соседних — так считал §380, и
 *            только ради сравнимости с ним она и остаётся;
 *   ОКНО 1   вершина В САМОЙ ячейке. После беззнакового перехода вершины
 *            появляются почти в каждой ячейке с геометрией, и широкое окно
 *            село бы на пол не потому, что поле полно, а потому, что окно
 *            широко. Честная мера — эта.
 * КЛАСС потерянной теперь один и он геометрический: пересечено ли хоть одно из
 * 12 рёбер самой ячейки. Класс «все восемь углов одного знака» исчез вместе со
 * знаком (§391 им и пользовался). */
typedef struct {
  int64_t nocc, lost27, lost1, lost_noedge;
} covstat;

static void occ_cover(const hz_dctree *T, const frame *fr, const unsigned char *occ,
                      const hz_htab *ht, covstat *st) {
  memset(st, 0, sizeof *st);
  for (int64_t z = 0; z < fr->n; z++)
    for (int64_t y = 0; y < fr->n; y++)
      for (int64_t x = 0; x < fr->n; x++) {
        if (!hz_occ_get(occ, hz_occ_index(fr->n, x, y, z))) continue;
        st->nocc++;
        int32_t c0[3] = {(int32_t)x, (int32_t)y, (int32_t)z};
        int here = dc_vert_at(T, c0);
        if (!here) st->lost1++;
        int near = here;
        for (int dz = -1; dz <= 1 && !near; dz++)
          for (int dy = -1; dy <= 1 && !near; dy++)
            for (int dx = -1; dx <= 1 && !near; dx++) {
              int32_t c[3] = {c0[0] + dx, c0[1] + dy, c0[2] + dz};
              if (c[0] < 0 || c[1] < 0 || c[2] < 0 || c[0] >= fr->n || c[1] >= fr->n ||
                  c[2] >= fr->n)
                continue;
              if (dc_vert_at(T, c)) near = 1;
            }
        if (!near) st->lost27++;
        if (here) continue;
        /* Почему вершины нет: ни одно ребро ячейки не пересечено. Это и есть
         * предсказанный остаток (А596) — кусок поверхности целиком внутри
         * ячейки, не задевший её рёбер. */
        int any = 0;
        for (int i = 0; i < 12 && !any; i++) {
          int a = i / 4, k = i % 4;
          int u = (a + 1) % 3, v = (a + 2) % 3;
          int32_t p[3] = {c0[0], c0[1], c0[2]};
          p[u] += k & 1;
          p[v] += (k >> 1) & 1;
          if (hz_htab_find(ht, a, p) != NULL) any = 1;
        }
        if (!any) st->lost_noedge++;
      }
}

/* --- 3b. КАРТИНКА ПОЛЯ (первая; §377) -------------------------------------- */

/* ЗАЧЕМ. До сих пор поле проверялось только числами, а глазами смотрели на
 * картинки ЛУЧЕВОГО пути (`prad`, полигональная сегментация). Пока поля не
 * видно, «одно представление» проверяется вслепую, и §368 п. 4 («сегментация
 * снимается, когда поле даёт картинку не хуже prad») не с чем сравнивать.
 *
 * ЧТО ЭТО НЕ ЕСТЬ. Это НЕ рендер переноса и не фронт: ни одного отскока, ни
 * тени, ни цвета материала. Полигоны выдаёт `hz_dc_walk`, они растеризуются
 * z-буфером, яркость — |n·l| от налобного источника. Единственное, что картинка
 * проверяет, — ГЕОМЕТРИЮ поля, и говорить о ней надо только это. */
typedef struct {
  const tr3_camera *cam;
  const frame *fr;
  double *z, *v;
  int64_t ntri, nclip;
} rastctx;

static void rast_tri(rastctx *R, const double a[3], const double b[3], const double c[3]) {
  const tr3_camera *cm = R->cam;
  double p[3][3] = {{a[0], a[1], a[2]}, {b[0], b[1], b[2]}, {c[0], c[1], c[2]}};
  double sx[3], sy[3], sz[3];
  for (int k = 0; k < 3; k++) {
    double d[3];
    for (int q = 0; q < 3; q++)
      d[q] = p[k][q] - cm->eye[q];
    double zz = d[0] * cm->fwd[0] + d[1] * cm->fwd[1] + d[2] * cm->fwd[2];
    if (!(zz > 1e-6)) { /* за камерой или в ней: отсекать надо, а не сдвигать */
      R->nclip++;
      return;
    }
    double rr = d[0] * cm->right[0] + d[1] * cm->right[1] + d[2] * cm->right[2];
    double uu = d[0] * cm->up[0] + d[1] * cm->up[1] + d[2] * cm->up[2];
    sz[k] = zz;
    sx[k] = ((rr / (zz * cm->tanx)) + 1.0) * 0.5 * (double)cm->w;
    sy[k] = (1.0 - uu / (zz * cm->tany)) * 0.5 * (double)cm->h;
  }
  /* нормаль СЧИТАЕТСЯ по мировым вершинам, а не берётся у ячейки: картинка
   * обязана показывать ту поверхность, которую выдал обход */
  double e1[3], e2[3], nn[3];
  for (int q = 0; q < 3; q++) {
    e1[q] = b[q] - a[q];
    e2[q] = c[q] - a[q];
  }
  nn[0] = e1[1] * e2[2] - e1[2] * e2[1];
  nn[1] = e1[2] * e2[0] - e1[0] * e2[2];
  nn[2] = e1[0] * e2[1] - e1[1] * e2[0];
  double m = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
  if (!(m > 0.0)) return;
  double lx = a[0] - cm->eye[0], ly = a[1] - cm->eye[1], lz = a[2] - cm->eye[2];
  double lm = sqrt(lx * lx + ly * ly + lz * lz);
  if (!(lm > 0.0)) return;
  double cosv = fabs((nn[0] * lx + nn[1] * ly + nn[2] * lz) / (m * lm));
  double shade = 0.15 + 0.85 * cosv;

  double x0 = sx[0], x1 = sx[0], y0 = sy[0], y1 = sy[0];
  for (int k = 1; k < 3; k++) {
    if (sx[k] < x0) x0 = sx[k];
    if (sx[k] > x1) x1 = sx[k];
    if (sy[k] < y0) y0 = sy[k];
    if (sy[k] > y1) y1 = sy[k];
  }
  int ix0 = (int)floor(x0), ix1 = (int)ceil(x1), iy0 = (int)floor(y0), iy1 = (int)ceil(y1);
  if (ix0 < 0) ix0 = 0;
  if (iy0 < 0) iy0 = 0;
  if (ix1 > cm->w - 1) ix1 = cm->w - 1;
  if (iy1 > cm->h - 1) iy1 = cm->h - 1;
  double d21x = sx[1] - sx[0], d21y = sy[1] - sy[0], d31x = sx[2] - sx[0], d31y = sy[2] - sy[0];
  double det = d21x * d31y - d21y * d31x;
  if (!(fabs(det) > 0.0)) return;
  R->ntri++;
  for (int py = iy0; py <= iy1; py++)
    for (int px = ix0; px <= ix1; px++) {
      double qx = (double)px + 0.5 - sx[0], qy = (double)py + 0.5 - sy[0];
      double u = (qx * d31y - qy * d31x) / det, v = (qy * d21x - qx * d21y) / det;
      if (u < 0.0 || v < 0.0 || u + v > 1.0) continue;
      double zz = sz[0] + u * (sz[1] - sz[0]) + v * (sz[2] - sz[0]);
      size_t k = (size_t)py * (size_t)cm->w + (size_t)px;
      if (zz >= R->z[k]) continue;
      R->z[k] = zz;
      R->v[k] = shade;
    }
}

static int rast_poly(void *ctx, const hz_dcref *ref, const double (*v)[3], int nv) {
  rastctx *R = (rastctx *)ctx;
  (void)ref;
  double w[4][3];
  for (int i = 0; i < nv; i++)
    for (int c = 0; c < 3; c++)
      w[i][c] = R->fr->org[c] + v[i][c] * R->fr->h;
  for (int i = 1; i + 1 < nv; i++)
    rast_tri(R, w[0], w[i], w[i + 1]);
  return 0;
}

/* --- 3в. ЭТАЛОН РЕГРЕССА: ВЫДАННЫЕ МНОГОУГОЛЬНИКИ (А666) ------------------- */

/* Ш2 меняет РАСКЛАДКУ, а не геометрию, поэтому его приёмка — совпадение выхода.
 * Сравнивать с тем, что было ДО Ш1б, нельзя: Ш1б законно добавляет потерянное.
 * Значит эталон снимается ЗДЕСЬ. Пишется текстом, а не двоичным: файл читают
 * глазами при разборе расхождения, и точность `%.17g` возвращает double
 * побитово. */
typedef struct {
  FILE *f;
  int64_t n;
} polydumpctx;

static int poly_dump(void *ctx, const hz_dcref *ref, const double (*v)[3], int nv) {
  polydumpctx *D = (polydumpctx *)ctx;
  (void)ref;
  fprintf(D->f, "%d", nv);
  for (int i = 0; i < nv; i++)
    for (int c = 0; c < 3; c++)
      fprintf(D->f, " %.17g", v[i][c]);
  fputc('\n', D->f);
  D->n++;
  return 0;
}

/* --- 3г. КРИТЕРИЙ АДАПТИВНОГО СРЕЗА (Ш3, §401) ----------------------------- */

/* «НЕВЯЗКА × ПРОЕКЦИЯ», ПОРОГ В ПИКСЕЛЯХ. Порог в метрах был бы магическим:
 * он зависел бы от масштаба сцены. В пикселях — нет, и потому число в нём
 * названо тем, что оно есть, — допуском наблюдателя.
 * ЗАВИСИМОСТЬ ОТ КАМЕРЫ ЗДЕСЬ ЗАКОННА И ЕДИНСТВЕННА ВО ВСЁМ УСТРОЙСТВЕ (А730):
 * LOD привязан к наблюдателю по построению, иначе цена становится input-bounded.
 * Сама  при этом чисто геометрическая и от камеры не зависит. */
typedef struct {
  double eye[3]; /* камера в координатах ДЕРЕВА, то есть в ячейках */
  double pxrad;  /* радиан на пиксель */
  double thr;    /* порог по НЕВЯЗКЕ, пикселей — ПОЛ по геометрии */
  /* ПОТОЛОК по угловому размеру ЯЧЕЙКИ, пикселей (§567). Задаётся отдельным
   * числом, а не равенством `thr`: это разные регуляторы, и связывать их одной
   * ручкой значит потерять возможность мерить их порознь. Ключ `ceil=`;
   * `0` выключает потолок и возвращает прежнее поведение (негативный контроль). */
  double thrcell;
  /* Пол по экранному размеру ячейки, пикселей; 0 — выключен. */
  double thrpoly;
  int64_t hist[HZ_DC_MAX_LOG2SIZE + 2];
  /* §580: отсечение по пирамиде видимости ВО ВРЕМЯ СПУСКА. Ставится ТОЛЬКО на
   * обход растеризатора: для переноса оно незаконно (свет приходит извне
   * кадра), для кадра — законно по определению пикселя. */
  int cull;
  const tr3_camera *cam;
  double org[3], h;
} lodctx;

static int lod_stop(void *ctx, const hz_dctree *t, const hz_dcref *r) {
  lodctx *L = (lodctx *)ctx;
  if (t->nd[r->ni].child0 < 0) return 1;
  /* Расстояние до БЛИЖАЙШЕЙ точки коробки (А733): у крупного узла оно разное в
   * разных его точках, а консервативно то, что даёт более мелкое дробление.
   * Камера внутри узла даёт ноль, и тогда спуск безусловен — деления на ноль не
   * возникает по построению, а не по проверке. */
  double d2 = 0.0;
  for (int a = 0; a < 3; a++) {
    double lo = (double)r->lo[a], hi = lo + (double)r->size, e = L->eye[a];
    double dd = e < lo ? lo - e : (e > hi ? e - hi : 0.0);
    d2 += dd * dd;
  }
  /* ОТСЕЧЕНИЕ ПО ПИРАМИДЕ ВИДИМОСТИ ВО ВРЕМЯ СПУСКА (§580). Узел, целиком
   * лежащий вне пирамиды, объявляется ЛИСТОМ — и всё его поддерево не
   * обходится вовсе. «Пропустить поддерево» в контракте `hz_dc_walk` не
   * предусмотрено, но «считать листом» даёт то же: лист выдаст пару полигонов,
   * их отбросит `lit_cull`, а спуск на тысячи узлов не пойдёт.
   *
   * ЭТО ЗАКОННО ЗДЕСЬ И НЕЗАКОННО В ПЕРЕНОСЕ, И РАЗНИЦУ НАДО НАЗВАТЬ. Проект
   * запрещает фрустумное отсечение для СВЕТА: он приходит извне кадра и от
   * невидимых поверхностей (`CLAUDE.md`, «FRUSTUM/VISIBILITY CULLING IS STILL
   * NOT AVAILABLE»). Здесь же делается КАДР, то есть пиксели: то, чего камера
   * не видит, в пиксель не попадает по определению. Поэтому флаг ставится
   * ТОЛЬКО на обход растеризатора, а срез для переноса строится тем же
   * `lod_stop` с флагом выключенным.
   *
   * ЗАМЕРЕНО ДО ПРАВКИ: обход выдавал `252 955` многоугольников, из которых
   * `181 626` (`71.8 %`) отбрасывались УЖЕ ПОСТРОЕННЫМИ. */
  if (L->cull && L->cam != NULL) {
    const tr3_camera *cm = L->cam;
    int outn = 0, outl = 0, outr = 0, outb = 0, outt = 0;
    for (int k = 0; k < 8; k++) {
      double p[3], d[3];
      for (int a = 0; a < 3; a++)
        p[a] = ((double)r->lo[a] + ((k >> a) & 1 ? (double)r->size : 0.0)) * L->h + L->org[a];
      for (int a = 0; a < 3; a++)
        d[a] = p[a] - cm->eye[a];
      double zz = d[0] * cm->fwd[0] + d[1] * cm->fwd[1] + d[2] * cm->fwd[2];
      double rr = d[0] * cm->right[0] + d[1] * cm->right[1] + d[2] * cm->right[2];
      double uu = d[0] * cm->up[0] + d[1] * cm->up[1] + d[2] * cm->up[2];
      if (!(zz > 1e-6)) outn++;
      if (rr < -zz * cm->tanx) outl++;
      if (rr > zz * cm->tanx) outr++;
      if (uu < -zz * cm->tany) outb++;
      if (uu > zz * cm->tany) outt++;
    }
    if (outn == 8 || outl == 8 || outr == 8 || outb == 8 || outt == 8) return 1;
  }
  if (!(d2 > 0.0)) return 0;
  /* ПОТОЛОК ПО УГЛОВОМУ РАЗМЕРУ ЯЧЕЙКИ (§567, замечание пользователя 08-12).
   * До этой правки в критерии стояла ТОЛЬКО невязка, то есть угловой размер
   * ОШИБКИ. У шумной геометрии (листва, орнамент) невязка держится порядка
   * полуячейки, и условие вырождалось в «дробить, пока ячейка не станет
   * `2·thr` пикселя» — то есть до предела разрешения, ровно там, где деталь
   * показать нечем. Замерено: у Сан-Мигеля выходило `0.95` ячейки НА ПИКСЕЛЬ
   * против `0.067` у комнаты, при вдвое большей дальности и трети кадра.
   * `CLAUDE.md` называет ТРИ регулятора, и первый из них — «LOD `L = εR` —
   * CEILING on element size». В коде был только второй (пол по геометрии).
   * ПРОВЕРЯЕТСЯ ПЕРВЫМ, и это не порядок ради порядка: потолок есть ГРАНИЦА
   * ПРИМЕНИМОСТИ невязки, а не ещё одно условие рядом. Невязку, не видимую в
   * кадре, мерить бессмысленно — её там нельзя ни показать, ни отличить. */
  double pxcell = ((double)r->size / sqrt(d2)) / L->pxrad;
  if (pxcell <= L->thrcell) return 1;
  /* ПОЛ ПО ЭКРАННОМУ РАЗМЕРУ (замечание пользователя 08-13). Элемент несёт не
   * только геометрию, но и СВЕТОВУЮ ПРОБУ: облучённость хранится в ячейке, а
   * растеризатор её интерполирует. Значит на полигоне в тысячу пикселей тень
   * не изобразится НИКАК — там три значения на всю площадь.
   * Прежний критерий (невязка) у плоской стены равен нулю и останавливал
   * дробление немедленно. ЗАМЕРЕНО: полигоны крупнее 64 пикселей кроют 45 %
   * экрана у Bistro и 43 % у Сан-Мигеля.
   * Поэтому: пока ячейка КРУПНЕЕ порога в пикселях, дробим независимо от того,
   * насколько она плоская. Это не поправка к невязке, а ДРУГОЕ ограничение —
   * по частоте СВЕТА, а не геометрии. */
  if (L->thrpoly > 0.0 && pxcell > L->thrpoly) return 0;
  double px = (hz_dc_rms(t, r->ni) / sqrt(d2)) / L->pxrad;
  return px <= L->thr;
}

/* СКОЛЬКО РАБОТЫ ПРИХОДИТ В СТАДИЮ, А НЕ ТОЛЬКО ЧТО ИЗ НЕЁ ВЫХОДИТ (А791, А807;
 * §457). Прежде чем судить о ПРАВИЛЕ выбора `edir` у крупного узла, надо знать,
 * СКОЛЬКО РАЗ оно вообще неоднозначно: «правило ни при чём» и «данные его ни
 * разу не проверили» — разные утверждения, и без этого счёта они неотличимы.
 *
 * Для каждого КРУПНОГО листа среза перебираются его 12 рёбер, и на каждом —
 * все накрытые сеточные рёбра из таблицы. Считается:
 *   - крупных рёбер с ХОТЯ БЫ ОДНИМ пересечением;
 *   - из них с БОЛЕЕ ЧЕМ ОДНИМ;
 *   - из них с РАЗНЫМИ знаками проекции нормали на ось (только здесь правило и
 *     значит что-нибудь);
 *   - на скольких ПРЕЖНЕЕ правило («первое пересечение от нижнего конца», в
 *     точности то, во что разворачивается «бит первого ребёнка с пересечением»)
 *     даёт бит, ОТЛИЧНЫЙ от знака суммы проекций.
 * Последнее число и есть верхняя оценка того, что вообще могла бы изменить
 * правка Ш7. Ноль в нём означает, что менять нечего. */
typedef struct {
  const hz_htab *ht;
  const hz_dctree *t;
  hz_dc_stop stop;
  void *sctx;
  int64_t nbig, nedge1, nedgem, ndis, nchange, nzero;
} ambctx;

static void amb_rec(ambctx *C, int32_t ni, const int32_t lo[3], int32_t size) {
  hz_dcref r = {ni, {lo[0], lo[1], lo[2]}, size};
  int leaf = C->t->nd[ni].child0 < 0 || (C->stop != NULL && C->stop(C->sctx, C->t, &r));
  if (!leaf) {
    int32_t half = size / 2;
    for (int k = 0; k < 8; k++) {
      int32_t clo[3];
      for (int a = 0; a < 3; a++)
        clo[a] = lo[a] + (((k >> a) & 1) ? half : 0);
      amb_rec(C, C->t->nd[ni].child0 + k, clo, half);
    }
    return;
  }
  if (size == 1) return;
  C->nbig++;
  for (int i = 0; i < 12; i++) {
    /* Нумерация рёбер — та же, что у `unit_edge` в dc.c и у сборки среза. */
    int a = i / 4, kk = i % 4, u = (a + 1) % 3, v = (a + 2) % 3;
    int32_t p0[3] = {lo[0], lo[1], lo[2]};
    p0[u] += (kk & 1) * size;
    p0[v] += ((kk >> 1) & 1) * size;
    int64_t npos = 0, nneg = 0;
    double sum = 0.0;
    int first = 0, got = 0;
    for (int32_t s = 0; s < size; s++) {
      int32_t p[3] = {p0[0], p0[1], p0[2]};
      p[a] += s;
      const hz_hedge *e = hz_htab_find(C->ht, a, p);
      if (e == NULL || e->in_lo == HZ_HEDGE_ERASED) continue;
      sum += e->nrm[a];
      if (e->nrm[a] > 0.0)
        npos++;
      else
        nneg++;
      if (!got) {
        got = 1;
        first = e->nrm[a] > 0.0;
      }
    }
    if (!got) continue;
    C->nedge1++;
    if (npos + nneg > 1) C->nedgem++;
    if (npos > 0 && nneg > 0) C->ndis++;
    /* ТОЧНЫЙ ноль, а не «почти»: порога здесь быть не может, ничья — это
     * ровное сокращение. Записано так, а не `sum == 0.0`, только чтобы не
     * будить -Wfloat-equal: сравнение остаётся точным. */
    if (!(sum < 0.0) && !(sum > 0.0)) C->nzero++;
    if ((sum > 0.0) != (first != 0)) C->nchange++;
  }
}

/* ПЛОЩАДЬ ВЫДАННОЙ ПОВЕРХНОСТИ (§446). Мера площадки `h²` на ЯЧЕЙКУ не может
 * быть верной для листа нулевой толщины: ячеек по обе стороны листа вдвое
 * больше, чем листов (§445). Площадь по ВЫДАННЫМ МНОГОУГОЛЬНИКАМ двойного счёта
 * не должна иметь — обход выдаёт многоугольник на ПЕРЕСЕЧЁННОЕ РЕБРО, а ребро у
 * листа одно, с какой бы стороны ни стояли ячейки. Здесь это проверяется числом.
 *
 * ДВА ВЕЕРА, А НЕ ОДИН (А796): площадь неплоского четырёхугольника от веера
 * ЗАВИСИТ — складка идёт по разной диагонали. Считаются оба, и печатается
 * разность: она и есть цена неоднозначности, а не погрешность. */
/* ГИСТОГРАММА ПО ПЛОСКОСТЯМ (§450, П2). `3 × 2 + 3 × 1 = 9` есть вывод из
 * арифметики, а не измерение: надо увидеть САМИ листы — сколько их, где они и на
 * каком расстоянии друг от друга. Ячейка делится на `HZ_PLBIN` долей, чтобы
 * различить два слоя, отстоящих на одну ячейку, и заодно увидеть смещение листа
 * ВНУТРИ ячейки. Корзины заводит вызывающий; NULL — гистограммы нет. */
#define HZ_PLBIN 8

typedef struct {
  double fan0, fan1; /* веер от v0 и от v1, в ЯЧЕЙКАХ² */
  double ax[3];      /* площадь по главной оси нормали (без знака), веер от v0 */
  int64_t px[3];     /* многоугольников по той же оси — различитель А798 */
  int64_t npoly, ntri, ndeg;
  double flatmax;  /* максимум неплоскостности в долях ячейки */
  double *pl_area; /* [3][nbin]: площадь по плоскостям, м² считает вызывающий */
  int64_t *pl_cnt; /* [3][nbin]: многоугольников по плоскостям */
  int32_t nbin;
  /* РАЗЛИЧИТЕЛЬ, БЕЗ КОТОРОГО ГИСТОГРАММА ПО ПЛОСКОСТЯМ НЕ РАЗДЕЛЯЕТ ДВЕ
   * ГИПОТЕЗЫ (П2). Вершина DC ставится НА ПОВЕРХНОСТЬ, поэтому два слоя ячеек по
   * обе стороны листа дадут вершины в ОДНОЙ плоскости — геометрически дубль
   * неотличим от «одно ребро учтено дважды». Отличает их только ЯЧЕЙКА: здесь
   * многоугольники раскладываются по координате своей ячейки вдоль главной оси. */
  double *cl_area; /* [3][ncell]: площадь по СЛОЯМ ЯЧЕЕК */
  int64_t *cl_cnt;
  int32_t ncell;
} areacnt;

static double tri_area2(const double a[3], const double b[3], const double c[3], double n[3]) {
  double u[3], v[3];
  for (int k = 0; k < 3; k++) {
    u[k] = b[k] - a[k];
    v[k] = c[k] - a[k];
  }
  n[0] = u[1] * v[2] - u[2] * v[1];
  n[1] = u[2] * v[0] - u[0] * v[2];
  n[2] = u[0] * v[1] - u[1] * v[0];
  return 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
}

/* Веер от вершины `b`: треугольники (b, b+i, b+i+1) по кругу. Для тройки оба
 * веера совпадают тождественно, для четвёрки — нет. */
static double fan_area(const double (*v)[3], int nv, int b, double nsum[3]) {
  double s = 0.0;
  if (nsum != NULL) nsum[0] = nsum[1] = nsum[2] = 0.0;
  for (int i = 1; i + 1 < nv; i++) {
    double n[3];
    double a = tri_area2(v[b], v[(b + i) % nv], v[(b + i + 1) % nv], n);
    s += a;
    if (nsum != NULL)
      for (int k = 0; k < 3; k++)
        nsum[k] += n[k];
  }
  return s;
}

static int area_emit(void *ctx, const hz_dcref *ref, const double (*v)[3], int nv) {
  areacnt *A = (areacnt *)ctx;
  double nsum[3];
  double s0 = fan_area(v, nv, 0, nsum);
  A->fan0 += s0;
  A->fan1 += fan_area(v, nv, 1, NULL);
  A->npoly++;
  A->ntri += nv - 2;
  /* ВЫРОЖДЕННЫЕ СЧИТАЮТСЯ ОТДЕЛЬНО (А797): у них направление нормали есть шум
   * округления, и в разбивку по осям их пускать нельзя. Порог — ТОЧНЫЙ ноль,
   * а не подобранная малость. */
  double ln = sqrt(nsum[0] * nsum[0] + nsum[1] * nsum[1] + nsum[2] * nsum[2]);
  if (!(ln > 0.0)) {
    A->ndeg++;
    return 0;
  }
  int ax = 0;
  for (int k = 1; k < 3; k++)
    if (fabs(nsum[k]) > fabs(nsum[ax])) ax = k;
  A->ax[ax] += s0;
  A->px[ax]++;
  if (A->pl_area != NULL) {
    /* Плоскость листа — координата ЦЕНТРОИДА многоугольника вдоль его главной
     * оси. Центроид, а не вершина: у наклонного многоугольника вершины разъедутся
     * по корзинам, а центроид — нет. */
    double c = 0.0;
    for (int i = 0; i < nv; i++)
      c += v[i][ax];
    c /= (double)nv;
    long long b = llround(c * (double)HZ_PLBIN);
    if (b >= 0 && b < (long long)A->nbin) {
      A->pl_area[(size_t)ax * (size_t)A->nbin + (size_t)b] += s0;
      A->pl_cnt[(size_t)ax * (size_t)A->nbin + (size_t)b]++;
    }
  }
  if (A->cl_area != NULL) {
    /* Ячейка берётся ПЕРВАЯ из четырёх (`ref[0]`): у многоугольника вокруг
     * ребра все четыре лежат в ОДНОМ слое вдоль оси ребра, а ось ребра совпадает
     * с главной осью нормали листа. Слой поэтому определён однозначно. */
    int32_t c0 = ref[0].lo[ax];
    if (c0 >= 0 && c0 < A->ncell) {
      A->cl_area[(size_t)ax * (size_t)A->ncell + (size_t)c0] += s0;
      A->cl_cnt[(size_t)ax * (size_t)A->ncell + (size_t)c0]++;
    }
  }
  /* Неплоскостность четвёрки: расстояние `v3` до плоскости первых трёх, в
   * долях ячейки (координаты обхода — в ячейках рамы). */
  if (nv == 4) {
    double n[3];
    double a012 = tri_area2(v[0], v[1], v[2], n);
    double l = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (a012 > 0.0 && l > 0.0) {
      double d = 0.0;
      for (int k = 0; k < 3; k++)
        d += n[k] * (v[3][k] - v[0][k]);
      d = fabs(d) / l;
      if (d > A->flatmax) A->flatmax = d;
    }
  }
  return 0;
}

/* ВЫДАННАЯ СЕТКА, СОБРАННАЯ ОБХОДОМ (А735). Треугольники складываются в один
 * массив и привязываются к КАЖДОЙ из четырёх своих ячеек односвязным списком:
 * поиск ближайшего к точке идёт тогда по кандидатам ЭТОЙ ячейки, а не по всей
 * сетке. Список, а не массив на ячейку, потому что число многоугольников на
 * ячейку заранее не известно и хвост у него длинный.
 *
 * ЗДЕСЬ ЖЕ СЧИТАЕТСЯ ОБРАЩЁННЫЙ ОБХОД (А729) — та проверка, без которой выбор
 * `edir` у крупного узла остаётся ничем не подкреплённым: нормаль
 * многоугольника (по Ньюэллу) сравнивается со средней нормалью его ячеек, и
 * несовпадение знака СЧИТАЕТСЯ. Адаптивный срез делает перепад уровней
 * повсеместным, поэтому если правило негодно, счётчик обязан это показать. */
/* ИНДЕКС ПО ГРУБОЙ СЕТКЕ, А НЕ ПО УЗЛУ — ПОПРАВКА К СЕБЕ (А738). Первая
 * редакция привязывала треугольник к его четырём ЯЧЕЙКАМ и искала ближайший
 * только среди них. Для крупной ячейки это опять мажоранта: точка у её края
 * ближе всего к треугольнику, привязанному к СОСЕДУ, а он в список не попадал.
 * Здесь треугольник кладётся во все клетки грубой сетки, которые задевает его
 * габарит, и поиск идёт по клетке точки плюс кольцо соседей — тогда величина
 * есть настоящее расстояние до выданной поверхности, а не оценка сверху. */
#define HZ_EGRID 6 /* сторона сетки индекса = 1 << HZ_EGRID */

typedef struct {
  double *v; /* 9 на треугольник: три вершины в координатах ДЕРЕВА */
  /* КЛАСС ТРЕУГОЛЬНИКА ПО ТОМУ, УЧАСТВОВАЛО ЛИ ПРАВИЛО (Ш7, А823). Обход берёт
   * бит ориентации у САМОЙ МЕЛКОЙ из четырёх ячеек ребра, поэтому свёртка
   * `edir` по детям участвует тогда и только тогда, когда минимальная из
   * четырёх — КРУПНАЯ. 0 = класс F (min size == 1, бит прямо из таблицы,
   * правкой не затронут), 1 = класс G (min size > 1). Класс F служит
   * ВНУТРЕННИМ эталоном в том же прогоне: доля обращённых в нём от правки
   * зависеть не может. */
  uint8_t *cls;
  int32_t ntri, cap;
  int32_t *head; /* на клетку индекса: первое вхождение, -1 — нет */
  int32_t *nxt;  /* следующее вхождение */
  int32_t *tri;  /* какому треугольнику принадлежит вхождение */
  int32_t nn, ncap;
  int32_t gn; /* сторона сетки индекса */
  int lev;
} emesh;

static int emesh_link(emesh *E, int32_t gi, int32_t ti) {
  if (E->nn >= E->ncap) {
    int32_t nc = E->ncap > 0 ? E->ncap * 2 : 4096;
    int32_t *a1 = realloc(E->nxt, (size_t)nc * sizeof *a1);
    if (a1 == NULL) return -1;
    E->nxt = a1;
    int32_t *a2 = realloc(E->tri, (size_t)nc * sizeof *a2);
    if (a2 == NULL) return -1;
    E->tri = a2;
    E->ncap = nc;
  }
  E->nxt[E->nn] = E->head[gi];
  E->tri[E->nn] = ti;
  E->head[gi] = E->nn;
  E->nn++;
  return 0;
}

static int emesh_push(emesh *E, const double a[3], const double b[3], const double c[3],
                      const hz_dcref *ref, int nv) {
  if (E->ntri >= E->cap) {
    int32_t nc = E->cap * 2;
    double *nv2 = realloc(E->v, (size_t)nc * 9 * sizeof *nv2);
    if (nv2 == NULL) return -1;
    E->v = nv2;
    uint8_t *nc2 = realloc(E->cls, (size_t)nc * sizeof *nc2);
    if (nc2 == NULL) return -1;
    E->cls = nc2;
    E->cap = nc;
  }
  double *d = E->v + 9 * (size_t)E->ntri;
  for (int k = 0; k < 3; k++) {
    d[k] = a[k];
    d[3 + k] = hz_flipall ? c[k] : b[k]; /* НК §463: обмотка обращена у ВСЕХ */
    d[6 + k] = hz_flipall ? b[k] : c[k];
  }
  int32_t ti = E->ntri++;
  /* ТРИ КЛАССА, А НЕ ДВА (правка по первому прогону, §457). Минимальная из
   * ячеек решает, ОТКУДА ВЗЯТ БИТ (А823); максимальная решает, ОГРУБЛЕНА ЛИ
   * ГЕОМЕТРИЯ полигона. Это разные вопросы, и слив их в один класс делает
   * величину неразличающей:
   *   0 = ВСЕ ЧЕТЫРЕ МЕЛКИЕ — ни правило, ни огрубление не участвуют (эталон);
   *   1 = min > 1 — бит из свёртки по детям, то есть предмет Ш7;
   *   2 = СМЕШАННЫЙ (min == 1, max > 1) — бит из таблицы, но геометрия
   *       огрублена: этот класс и отделяет вклад ПРАВИЛА от вклада ОГРУБЛЕНИЯ. */
  int32_t msz = ref[0].size, xsz = ref[0].size;
  for (int k = 1; k < nv; k++) {
    if (ref[k].size < msz) msz = ref[k].size;
    if (ref[k].size > xsz) xsz = ref[k].size;
  }
  /* Биты 0-1 — класс уровня (§457), бит 2 — «родом из сложенного квада» (П2). */
  E->cls[ti] = (uint8_t)((msz > 1 ? 1 : (xsz > 1 ? 2 : 0)) | (g_curfold ? 4 : 0));
  /* Габарит треугольника в клетках индекса. Координаты — в ячейках сетки уровня
   * `lev`, клетка индекса шире в `1 << (lev - HZ_EGRID)` раз. */
  int sh = E->lev - HZ_EGRID;
  if (sh < 0) sh = 0;
  int32_t g0[3] = {0, 0, 0}, g1[3] = {0, 0, 0};
  for (int k = 0; k < 3; k++) {
    double lo = d[k], hi = d[k];
    for (int q = 1; q < 3; q++) {
      if (d[3 * q + k] < lo) lo = d[3 * q + k];
      if (d[3 * q + k] > hi) hi = d[3 * q + k];
    }
    int32_t i0 = (int32_t)floor(lo) >> sh, i1 = (int32_t)floor(hi) >> sh;
    if (i0 < 0) i0 = 0;
    if (i1 >= E->gn) i1 = E->gn - 1;
    g0[k] = i0;
    g1[k] = i1;
  }
  for (int32_t z = g0[2]; z <= g1[2]; z++)
    for (int32_t y = g0[1]; y <= g1[1]; y++)
      for (int32_t x = g0[0]; x <= g1[0]; x++)
        if (emesh_link(E, (int32_t)hz_occ_index(E->gn, x, y, z), ti) != 0) return -1;
  return 0;
}

/* Нормаль треугольника по трём точкам — векторное произведение, НЕ нормированное:
 * дальше нужен только знак и длина как площадь. */
static void tri_nrm(const double a[3], const double b[3], const double c[3], double n[3]) {
  double e1[3], e2[3];
  for (int k = 0; k < 3; k++) {
    e1[k] = b[k] - a[k];
    e2[k] = c[k] - a[k];
  }
  n[0] = e1[1] * e2[2] - e1[2] * e2[1];
  n[1] = e1[2] * e2[0] - e1[0] * e2[2];
  n[2] = e1[0] * e2[1] - e1[1] * e2[0];
}

/* Ш10 (§467): СЛОЖЕННЫЕ КВАДЫ И ВЫБОР ДИАГОНАЛИ. Считаются глобально, потому
 * что обход зовёт `emesh_emit` через указатель и второго канала для чисел нет. */
static int64_t g_quad = 0, g_fold = 0;

static int emesh_emit(void *ctx, const hz_dcref *ref, const double (*v)[3], int nv) {
  emesh *E = (emesh *)ctx;
  int i0 = 0;
  g_curfold = 0;
  if (nv == 4) {
    /* СЛОЖЕННОСТЬ — свойство самого квада, эталон для неё не нужен: две
     * половины веера сравниваются между собой. Знак берётся у веера от `v[0]`,
     * то есть у того разбиения, которое было до этого шага. */
    double n1[3], n2[3];
    tri_nrm(v[0], v[1], v[2], n1);
    tri_nrm(v[0], v[2], v[3], n2);
    g_quad++;
    if (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2] < 0.0) {
      g_fold++;
      g_curfold = 1;
    }
    /* ДИАГОНАЛЬ ПО КОРОТКОЙ. Веер от `v[0]` режет по `v0–v2`, веер от `v[1]` —
     * по `v1–v3`; берётся та, что короче. У неплоского квада это не косметика:
     * длинная диагональ проходит дальше от поверхности, и её половины сильнее
     * расходятся по нормали (та же причина, по которой А796 требовал считать
     * ОБА веера для площади). */
    double d02 = 0.0, d13 = 0.0;
    for (int k = 0; k < 3; k++) {
      double x = v[0][k] - v[2][k], y = v[1][k] - v[3][k];
      d02 += x * x;
      d13 += y * y;
    }
    int shortis1 = d13 < d02;
    i0 = (hz_longdiag ? !shortis1 : shortis1) ? 1 : 0;
  }
  double a[3] = {v[i0][0], v[i0][1], v[i0][2]};
  for (int i = 1; i + 1 < nv; i++) {
    const double *b = v[(i0 + i) % nv], *c = v[(i0 + i + 1) % nv];
    if (emesh_push(E, a, b, c, ref, nv) != 0) return 1;
  }
  return 0;
}

/* Ближайшая точка треугольника: квадрат расстояния от точки до треугольника.
 * Разбор по областям Вороного (Эриксон, «Real-Time Collision Detection»);
 * величина ТОЧНАЯ, а не приближение по вершинам, — в этом весь смысл правки
 * А735. */
static double pt_tri_d2(const double p[3], const double a[3], const double b[3],
                        const double c[3]) {
  double ab[3], ac[3], ap[3];
  for (int k = 0; k < 3; k++) {
    ab[k] = b[k] - a[k];
    ac[k] = c[k] - a[k];
    ap[k] = p[k] - a[k];
  }
  double d1 = ab[0] * ap[0] + ab[1] * ap[1] + ab[2] * ap[2];
  double d2 = ac[0] * ap[0] + ac[1] * ap[1] + ac[2] * ap[2];
  double q[3];
  if (d1 <= 0.0 && d2 <= 0.0) {
    for (int k = 0; k < 3; k++)
      q[k] = a[k];
  } else {
    double bp[3];
    for (int k = 0; k < 3; k++)
      bp[k] = p[k] - b[k];
    double d3 = ab[0] * bp[0] + ab[1] * bp[1] + ab[2] * bp[2];
    double d4 = ac[0] * bp[0] + ac[1] * bp[1] + ac[2] * bp[2];
    if (d3 >= 0.0 && d4 <= d3) {
      for (int k = 0; k < 3; k++)
        q[k] = b[k];
    } else {
      double cp[3];
      for (int k = 0; k < 3; k++)
        cp[k] = p[k] - c[k];
      double d5 = ab[0] * cp[0] + ab[1] * cp[1] + ab[2] * cp[2];
      double d6 = ac[0] * cp[0] + ac[1] * cp[1] + ac[2] * cp[2];
      double vc = d1 * d4 - d3 * d2, vb = d5 * d2 - d1 * d6, va = d3 * d6 - d5 * d4;
      if (d6 >= 0.0 && d5 <= d6) {
        for (int k = 0; k < 3; k++)
          q[k] = c[k];
      } else if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        double t = d1 / (d1 - d3);
        for (int k = 0; k < 3; k++)
          q[k] = a[k] + t * ab[k];
      } else if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        double t = d2 / (d2 - d6);
        for (int k = 0; k < 3; k++)
          q[k] = a[k] + t * ac[k];
      } else if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0) {
        double t = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        for (int k = 0; k < 3; k++)
          q[k] = b[k] + t * (c[k] - b[k]);
      } else {
        double den = 1.0 / (va + vb + vc);
        double u = vb * den, w = vc * den;
        for (int k = 0; k < 3; k++)
          q[k] = a[k] + ab[k] * u + ac[k] * w;
      }
    }
  }
  double s = 0.0;
  for (int k = 0; k < 3; k++) {
    double d = p[k] - q[k];
    s += d * d;
  }
  return s;
}

/* Ближайшее расстояние от точки до ВЫДАННОЙ поверхности. Ищется по клетке точки
 * и кольцу соседей: радиус растёт, пока не найдено хоть что-то, но не дальше
 * половины сетки — иначе на пустой сцене цикл стал бы полным перебором молча. */
static double emesh_nearest(const emesh *E, const frame *fr, const double p[3]) {
  int sh = E->lev - HZ_EGRID;
  if (sh < 0) sh = 0;
  int32_t g[3];
  for (int k = 0; k < 3; k++) {
    double f = floor((p[k] - fr->org[k]) / fr->h);
    int32_t i = (int32_t)f >> sh;
    if (i < 0) i = 0;
    if (i >= E->gn) i = E->gn - 1;
    g[k] = i;
  }
  double best = 1e300;
  for (int32_t r = 0; r <= E->gn / 2; r++) {
    for (int32_t dz = -r; dz <= r; dz++)
      for (int32_t dy = -r; dy <= r; dy++)
        for (int32_t dx = -r; dx <= r; dx++) {
          /* только КОЖУРА кольца: внутренность просмотрена на прошлом радиусе */
          if (r > 0 && dx > -r && dx < r && dy > -r && dy < r && dz > -r && dz < r) continue;
          int32_t x = g[0] + dx, y = g[1] + dy, z = g[2] + dz;
          if (x < 0 || y < 0 || z < 0 || x >= E->gn || y >= E->gn || z >= E->gn) continue;
          for (int32_t e = E->head[hz_occ_index(E->gn, x, y, z)]; e >= 0; e = E->nxt[e]) {
            const double *tv = E->v + 9 * (size_t)E->tri[e];
            double w[3][3];
            for (int q = 0; q < 3; q++)
              for (int c = 0; c < 3; c++)
                w[q][c] = fr->org[c] + tv[3 * q + c] * fr->h;
            double d = pt_tri_d2(p, w[0], w[1], w[2]);
            if (d < best) best = d;
          }
        }
    /* Найденное на радиусе r могло оказаться дальше, чем есть на r+1, поэтому
     * после первой находки просматривается ещё одно кольцо — и только потом
     * выход. Без этого величина была бы «первое попавшееся», а не ближайшее. */
    if (best < 1e299 && r > 0) break;
  }
  return best;
}
/* Ключ ячейки среза для сличения ДВУХ срезов множествами (не построчно, А720). */
static uint64_t cellkey(const hz_dccell *c) {
  return ((uint64_t)c->lvl << 48) | ((uint64_t)c->lo[0] << 32) | ((uint64_t)c->lo[1] << 16) |
         (uint64_t)c->lo[2];
}

/* Пара «ключ ячейки — её номер в срезе» для сортировки индекса цвета (§590).
 * Сравнение вторым полем делает порядок ОДНОЗНАЧНЫМ при равных ключах, то есть
 * заменяет устойчивость, которой у `qsort` нет. */
typedef struct {
  uint64_t k;
  int32_t o;
} cellkv;

static int cmp_cellkv(const void *x, const void *y) {
  const cellkv *a = (const cellkv *)x, *b = (const cellkv *)y;
  if (a->k != b->k) return a->k < b->k ? -1 : 1;
  return a->o < b->o ? -1 : (a->o > b->o ? 1 : 0);
}

static int cmp_u64(const void *x, const void *y) {
  uint64_t a = *(const uint64_t *)x, b = *(const uint64_t *)y;
  return a < b ? -1 : (a > b ? 1 : 0);
}

/* ОБРАЩЁННЫЙ ОБХОД (А729) — проверка выбора `edir` у крупного узла. Нормаль
 * выданного треугольника сравнивается с нормалью ГРАНИ ближайшего исходного
 * треугольника в той же ячейке.
 *
 * ОГОВОРКА, БЕЗ КОТОРОЙ ЧИСЛО НЕ ЧИТАЕТСЯ: обход исходной модели сам может быть
 * несогласован (суп треугольников), и часть несовпадений придёт оттуда, а не от
 * нас. Поэтому величина берётся ДИФФЕРЕНЦИАЛЬНО: та же доля считается при
 * ПОЛНОЙ ГЛУБИНЕ, где перепада уровней нет вовсе, и сравниваются ДВЕ доли.
 * Разность и есть вклад перепада; общий уровень — свойство модели. */
/* Разбор вывернутости (§463). Складывается рядом с самим счётом, чтобы
 * популяция была та же самая до последнего треугольника. `NULL` = не считать. */
typedef struct {
  int64_t hist[7], histf[7], nzn, n, cap, nfoldtri, nfoldfl;
  double asum, afl;
  double *ar, *as;
  uint8_t *fl;
} flipdiag;

static void flipdiag_init(flipdiag *D, int32_t ntri) {
  memset(D, 0, sizeof *D);
  D->cap = ntri > 0 ? ntri : 1;
  D->ar = malloc((size_t)D->cap * sizeof *D->ar);
  D->as = malloc((size_t)D->cap * sizeof *D->as);
  D->fl = malloc((size_t)D->cap * sizeof *D->fl);
  if (D->ar == NULL || D->as == NULL || D->fl == NULL) exit(1);
}

/* Медиана по подмножеству (только вывернутые, либо все) — медиана, а не среднее:
 * у площадей длинный хвост, и среднее по нему сказало бы о хвосте, а не о том,
 * каков типичный треугольник. */
static double flipdiag_med(const flipdiag *D, const double *v, int onlyflip) {
  int64_t k = 0;
  double *tmp = malloc((size_t)(D->n > 0 ? D->n : 1) * sizeof *tmp);
  if (tmp == NULL) exit(1);
  for (int64_t i = 0; i < D->n; i++)
    if (!onlyflip || D->fl[i]) tmp[k++] = v[i];
  double r = 0.0;
  if (k > 0) {
    qsort(tmp, (size_t)k, sizeof *tmp, cmp_d);
    r = tmp[k / 2];
  }
  free(tmp);
  return r;
}

static void flipdiag_report(flipdiag *D, const char *what) {
  static const char *nm[7] = {"[-1,-.9)", "[-.9,-.5)", "[-.5,-.1)", "[-.1,.1)",
                              "[.1,.5)",  "[.5,.9)",   "[.9,1]"};
  int64_t nfl = 0;
  for (int b = 0; b < 7; b++)
    nfl += D->histf[b];
  printf("      §463 РАЗБОР ВЫВЕРНУТОСТИ (%s): доля ПО ПЛОЩАДИ %.3f %% (по числу считалась выше); "
         "медиана площади вывернутых %.4f h² против всех %.4f h²; медиана АСПЕКТА вывернутых "
         "%.4f против всех %.4f; нормаль ровно нулевая у %lld\n",
         what, 100.0 * D->afl / (D->asum > 0.0 ? D->asum : 1.0), flipdiag_med(D, D->ar, 1),
         flipdiag_med(D, D->ar, 0), flipdiag_med(D, D->as, 1), flipdiag_med(D, D->as, 0),
         (long long)D->nzn);
  printf("         П2: треугольников из СЛОЖЕННЫХ квадов %lld, из них вывернутых %lld (%.1f %%); а "
         "среди ВСЕХ вывернутых доля родом из сложенных %.1f %%\n",
         (long long)D->nfoldtri, (long long)D->nfoldfl,
         100.0 * (double)D->nfoldfl / (double)(D->nfoldtri ? D->nfoldtri : 1),
         100.0 * (double)D->nfoldfl / (double)(nfl ? nfl : 1));
  printf("         КОСИНУС вывернутых по корзинам:");
  for (int b = 0; b < 7; b++)
    if (D->histf[b] > 0)
      printf(" %s %lld (%.1f %%)", nm[b], (long long)D->histf[b],
             100.0 * (double)D->histf[b] / (double)(nfl ? nfl : 1));
  printf("\n");
  free(D->ar);
  free(D->as);
  free(D->fl);
  memset(D, 0, sizeof *D);
}

static int64_t emesh_flips(const emesh *E, celltris *CT, const frame *fr, const hz_objmesh *m,
                           int64_t *ncmp, int64_t nfc[3], int64_t ncc[3], flipdiag *D) {
  int64_t nf = 0, nc = 0;
  for (int k = 0; k < 3; k++)
    nfc[k] = ncc[k] = 0;
  for (int32_t i = 0; i < E->ntri; i++) {
    const double *tv = E->v + 9 * (size_t)i;
    double w[3][3], ctr[3] = {0, 0, 0};
    for (int q = 0; q < 3; q++)
      for (int c = 0; c < 3; c++) {
        w[q][c] = fr->org[c] + tv[3 * q + c] * fr->h;
        ctr[c] += w[q][c] / 3.0;
      }
    double e1[3], e2[3], nn[3];
    for (int c = 0; c < 3; c++) {
      e1[c] = w[1][c] - w[0][c];
      e2[c] = w[2][c] - w[0][c];
    }
    nn[0] = e1[1] * e2[2] - e1[2] * e2[1];
    nn[1] = e1[2] * e2[0] - e1[0] * e2[2];
    nn[2] = e1[0] * e2[1] - e1[1] * e2[0];
    int32_t cell[3];
    int ok = 1;
    for (int c = 0; c < 3; c++) {
      double f = floor((ctr[c] - fr->org[c]) / fr->h);
      if (!(f >= 0.0) || !(f < (double)fr->n)) ok = 0;
      cell[c] = ok ? (int32_t)f : 0;
    }
    if (!ok) continue;
    const int32_t *ls = NULL;
    int32_t nls = ct_list(CT, cell, &ls);
    if (nls == 0) continue;
    /* ЭТАЛОННЫЙ ТРЕУГОЛЬНИК — БЛИЖАЙШИЙ, А НЕ ПЕРВЫЙ (§457). Первый в списке —
     * произвольный: в приграничной ячейке сферы лежат и передняя, и задняя
     * грань, и их нормали смотрят ВРОЗЬ, так что «первый» давал бы вывернутость
     * там, где её нет. Это ровно класс «подмена величины»: мажоранта вместо
     * точного. Ближайший ищется по расстоянию от ЦЕНТРА выданного треугольника
     * до исходного — той же точной формулой, что и ошибка поверхности. */
    int32_t bi = ls[0];
    if (nls > 1) {
      double bd = 1e300;
      for (int32_t q = 0; q < nls; q++) {
        const double *Aq, *Bq, *Cq;
        tri_verts(m, ls[q], &Aq, &Bq, &Cq);
        double dq = pt_tri_d2(ctr, Aq, Bq, Cq);
        if (dq < bd) {
          bd = dq;
          bi = ls[q];
        }
      }
    }
    const double *A2, *B2, *C2;
    tri_verts(m, bi, &A2, &B2, &C2);
    double f1[3], f2[3], fn[3];
    for (int c = 0; c < 3; c++) {
      f1[c] = B2[c] - A2[c];
      f2[c] = C2[c] - A2[c];
    }
    fn[0] = f1[1] * f2[2] - f1[2] * f2[1];
    fn[1] = f1[2] * f2[0] - f1[0] * f2[2];
    fn[2] = f1[0] * f2[1] - f1[1] * f2[0];
    double d = nn[0] * fn[0] + nn[1] * fn[1] + nn[2] * fn[2];
    int cl = (int)(E->cls[i] & 3u);
    if (cl > 2) cl = 0;
    if (D != NULL && (E->cls[i] & 4u)) {
      D->nfoldtri++;
      if (d < 0.0) D->nfoldfl++;
    }
    nc++;
    ncc[cl]++;
    if (d < 0.0) {
      nf++;
      nfc[cl]++;
    }
    if (D == NULL) continue;
    /* §463: РАЗДЕЛИТЬ НАСТОЯЩИЙ ПЕРЕВОРОТ И ВЫРОЖДЕННЫЙ ЧЕТЫРЁХУГОЛЬНИК.
     * Косинус НОРМИРУЕТСЯ (А848): у знака длина не нужна, у гистограммы —
     * нужна. Точное вырождение (нулевая нормаль) идёт в СВОЙ счётчик, а не в
     * среднюю корзину, иначе исход (Б) подтверждался бы тем, что в него же и
     * записано. */
    double ln = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
    double lf = sqrt(fn[0] * fn[0] + fn[1] * fn[1] + fn[2] * fn[2]);
    double ar = 0.5 * ln / (fr->h * fr->h); /* площадь в долях h² */
    D->asum += ar;
    if (d < 0.0) D->afl += ar;
    if (!(ln > 0.0) || !(lf > 0.0)) {
      D->nzn++;
      continue;
    }
    double cs = d / (ln * lf);
    if (cs > 1.0) cs = 1.0;
    if (cs < -1.0) cs = -1.0;
    static const double edge[6] = {-0.9, -0.5, -0.1, 0.1, 0.5, 0.9};
    int b = 0;
    while (b < 6 && cs >= edge[b])
      b++;
    D->hist[b]++;
    if (d < 0.0) D->histf[b]++;
    /* АСПЕКТ (А850) — признак ИЗ ДРУГОЙ ПРИРОДЫ, чем площадь: у иглы он около
     * нуля при любом масштабе. `4A/(√3 L²)`, где L — длиннейшая сторона; у
     * равностороннего равен 1. */
    double e3[3];
    for (int c = 0; c < 3; c++)
      e3[c] = w[2][c] - w[1][c];
    double l2m = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
    double t2 = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
    double t3 = e3[0] * e3[0] + e3[1] * e3[1] + e3[2] * e3[2];
    if (t2 > l2m) l2m = t2;
    if (t3 > l2m) l2m = t3;
    double asp = l2m > 0.0 ? (2.0 * ln) / (1.7320508075688772 * l2m) : 0.0;
    if (D->n < D->cap) {
      D->ar[D->n] = ar;
      D->as[D->n] = asp;
      D->fl[D->n] = d < 0.0 ? 1 : 0;
      D->n++;
    }
  }
  if (ncmp != NULL) *ncmp = nc;
  return nf;
}

/* --- 3д. РАЗРУШЕНИЕ (Ш4, Р8) ----------------------------------------------- */

/* УДАЛЕНИЕ ПОВЕРХНОСТИ СФЕРОЙ. По Р8 это одна из двух определённых операций:
 * стираются эрмитовы образцы рёбер, ЦЕЛИКОМ лежащих внутри сферы. «Внутри»
 * принадлежит ПРИМИТИВУ — у него знак свой и аналитический, поле в этом не
 * участвует, и потому отсутствие знака у поля здесь ничему не мешает.
 *
 * КРАЙ ДЫРЫ ВЫЙДЕТ СТУПЕНЧАТЫМ, и это сознательно (А743): ребро либо целиком
 * внутри, либо нет. Доля рёбер, накрытых ЧАСТИЧНО, считается и печатается —
 * если она велика, гладкий край становится обязательным.
 *
 * Стирание — пометкой в поле in_lo, а не удалением из массива: таблица
 * отсортирована по ключу, и вырезание записи сдвинуло бы всё за ней, то есть
 * стоило бы O(таблицы) вместо O(затронутых). Поиск `hz_htab_find` при этом
 * обязан отдавать помеченные как «нет ребра» — это делает вызывающий. */
typedef struct {
  int64_t nedge, npart, ncell;
  int32_t lo[3], hi[3]; /* коробка затронутых ЯЧЕЕК */
} hitstat;

static int hit_sphere(hz_htab *ht, const frame *fr, const double c[3], double rad, hitstat *S) {
  memset(S, 0, sizeof *S);
  for (int a = 0; a < 3; a++) {
    S->lo[a] = fr->n;
    S->hi[a] = -1;
  }
  for (int32_t i = 0; i < ht->n; i++) {
    hz_hedge *e = &ht->e[i];
    if (e->in_lo == HZ_HEDGE_ERASED) continue; /* уже стёрто */
    /* Оба конца ребра в мировых координатах. */
    double p0[3], p1[3];
    for (int k = 0; k < 3; k++) {
      p0[k] = fr->org[k] + (double)e->p[k] * fr->h;
      p1[k] = p0[k];
    }
    p1[e->axis] += fr->h;
    double d0 = 0.0, d1 = 0.0;
    for (int k = 0; k < 3; k++) {
      double a0 = p0[k] - c[k], a1 = p1[k] - c[k];
      d0 += a0 * a0;
      d1 += a1 * a1;
    }
    int in0 = d0 <= rad * rad, in1 = d1 <= rad * rad;
    if (!in0 && !in1) continue;
    if (in0 != in1) {
      S->npart++;
      continue; /* частично накрытое ребро НЕ трогаем — А743 */
    }
    e->in_lo = HZ_HEDGE_ERASED;
    S->nedge++;
    for (int a = 0; a < 3; a++) {
      int32_t l = e->p[a] - 1, h = e->p[a];
      if (l < 0) l = 0;
      if (h >= fr->n) h = fr->n - 1;
      if (l < S->lo[a]) S->lo[a] = l;
      if (h > S->hi[a]) S->hi[a] = h;
    }
  }
  return 0;
}

/* УДАР ПО ЗАНЯТОСТИ (§524, пункт 2; находка А918). Стирание рёбер эрмитовой
 * таблицы до свипа НЕ ДОХОДИТ: дерево свипа строится из ПИРАМИДЫ (`stree_rec`
 * читает `P->b`, заслон листа — тоже из неё). Без этой функции доля затронутых
 * листьев вышла бы `0.000 %` не потому, что грязь локальна, а потому, что её
 * нет вовсе, — то есть ложным нулём, четвёртым за проект.
 *
 * ПРАВИЛО ТО ЖЕ, ЧТО У `hit_sphere` (А743): гасится ячейка, коробка которой
 * ЦЕЛИКОМ внутри сферы; частично накрытая не трогается, и край дыры выходит
 * ступенчатым сознательно. Число частично накрытых считается и печатается.
 *
 * ПОДЪЁМ ПО ПИРАМИДЕ — ТОЛЬКО НАД ЗАТРОНУТОЙ КОРОБКОЙ: родитель есть ИЛИ восьми
 * детей (тождество `opyr_build`), и измениться могут лишь предки погашенных.
 * Полный пересчёт стоил бы O(объёма) на каждый удар и мерил бы не то. */
static int64_t hit_occ(opyr *P, const frame *fr, const double c[3], double rad, int64_t *npart) {
  int32_t lo[3], hi[3];
  for (int a = 0; a < 3; a++) {
    double l = floor((c[a] - rad - fr->org[a]) / fr->h) - 1.0;
    double h = floor((c[a] + rad - fr->org[a]) / fr->h) + 1.0;
    if (l < 0.0) l = 0.0;
    if (h > (double)(fr->n - 1)) h = (double)(fr->n - 1);
    lo[a] = (int32_t)l;
    hi[a] = (int32_t)h;
  }
  int64_t nclr = 0;
  *npart = 0;
  unsigned char *bt = P->b[P->lev];
  for (int32_t z = lo[2]; z <= hi[2]; z++)
    for (int32_t y = lo[1]; y <= hi[1]; y++)
      for (int32_t x = lo[0]; x <= hi[0]; x++) {
        size_t ci = hz_occ_index(fr->n, x, y, z);
        if (!hz_occ_get(bt, ci)) continue;
        /* Дальний и ближний углы коробки ячейки от центра сферы: целиком внутри
         * тогда и только тогда, когда ДАЛЬНИЙ угол ближе радиуса. */
        double cl[3], ch[3];
        cell_box(fr, x, y, z, cl, ch);
        double dfar = 0.0, dnear = 0.0;
        for (int a = 0; a < 3; a++) {
          double d0 = fabs(cl[a] - c[a]), d1 = fabs(ch[a] - c[a]);
          double dm = d0 > d1 ? d0 : d1;
          double dn = c[a] < cl[a] ? cl[a] - c[a] : (c[a] > ch[a] ? c[a] - ch[a] : 0.0);
          dfar += dm * dm;
          dnear += dn * dn;
        }
        if (dfar <= rad * rad) {
          bt[ci >> 3] = (unsigned char)(bt[ci >> 3] & ~(1u << (ci & 7)));
          nclr++;
        } else if (dnear <= rad * rad)
          (*npart)++;
      }
  /* Подъём: на каждом уровне пересчитать родителей затронутой коробки. */
  int32_t plo[3] = {lo[0], lo[1], lo[2]}, phi[3] = {hi[0], hi[1], hi[2]};
  for (int l = P->lev; l > 0; l--) {
    int32_t n = (int32_t)1 << l, n2 = n / 2;
    for (int a = 0; a < 3; a++) {
      plo[a] >>= 1;
      phi[a] >>= 1;
    }
    for (int32_t z = plo[2]; z <= phi[2]; z++)
      for (int32_t y = plo[1]; y <= phi[1]; y++)
        for (int32_t x = plo[0]; x <= phi[0]; x++) {
          int any = 0;
          for (int k = 0; k < 8 && !any; k++)
            any = hz_occ_get(P->b[l], hz_occ_index(n, 2 * x + (k & 1), 2 * y + ((k >> 1) & 1),
                                                   2 * z + ((k >> 2) & 1))) != 0;
          size_t pi = hz_occ_index(n2, x, y, z);
          if (any)
            hz_occ_set(P->b[l - 1], pi);
          else
            P->b[l - 1][pi >> 3] = (unsigned char)(P->b[l - 1][pi >> 3] & ~(1u << (pi & 7)));
        }
  }
  return nclr;
}

/* --- 3е. ПРЯМОЙ СВЕТ ПО СРЕЗУ (Ш5, часть первая) --------------------------- */

/* ЧЕСТНОЕ ОТСТУПЛЕНИЕ ОТ ПЛАНА, ЗАПИСАННОЕ, А НЕ СДЕЛАННОЕ МОЛЧА. §383 требует
 * СВИПА ПОТОКОМ в октантном порядке. Здесь сделан не свип, а ЗАТЕНЕНИЕ ЛУЧОМ по
 * пирамиде занятости: от ячейки к образцам площадки. Причина — не удобство:
 * свип переносит РАДИАНТНОСТЬ по направлениям и требует соседа в срезе, то есть
 * индекса, которого ещё нет; луч же пользуется УЖЕ ПОСТРОЕННОЙ пирамидой (Р7) и
 * даёт прямой свет с тенями немедленно. Цена названа: это `O(ячейки × образцы ×
 * шаги)`, то есть по построению дороже свипа, и в приёмку П5.1 (нс на ячейку НА
 * СВИП) оно НЕ ЗАСЧИТЫВАЕТСЯ. Отскока здесь нет вовсе — средство 5 по-прежнему
 * не начато, и одна эта часть его не закрывает.
 *
 * ИСТОЧНИК ПРОТЯЖЁННЫЙ И ПОСТАВЛЕН ОСМЫСЛЕННО: площадка под потолком зала.
 * Налобный источник для теней негоден по построению — он их не отбрасывает.
 * ЦВЕТ ТРЁХКАНАЛЬНЫЙ И НЕ ФИКТИВНЫЙ (А757): источник ОКРАШЕН, иначе три
 * одинаковых числа выдавались бы за RGB. */
/* АЛЬБЕДО ПО КАНАЛАМ ИЗ ТАБЛИЦЫ МАТЕРИАЛОВ. `kd3` читается загрузчиком из
 * `.mtl` и до сих пор в поле не доходил вовсе (А761: «цвет формально
 * трёхканальный, фактически однотонный»). Здесь он доходит. */
/* ФАЛЬСИФИКАТОР §529 (`albone`): все альбедо РОВНО единица. Тогда лишний
 * множитель `alb(j)` есть умножение на `1`, и ключ `emitalb2` обязан не изменить
 * НИЧЕГО. Изменит — значит лишний множитель не альбедо, и разбор §529 неверен.
 * Назван фальсификатором, а не негативным контролем: он предсказан УСПЕШНЫМ, то
 * есть это проверка эквивалентности (А10), и путать их нельзя. */
static int g_albone = 0;

static double alb(const hz_objmesh *m, uint8_t mi, int k) {
  if (g_albone) return 1.0;
  if (m->mtl == NULL || mi >= m->nmtl) return 0.5;
  double a = m->mtl[mi].kd3[k];
  return a >= 0.0 && a <= 1.0 ? a : 0.5;
}

typedef struct {
  double c[3]; /* центр площадки */

  double u[3], v[3]; /* полуоси */
  double rgb[3];     /* сила по каналам */
} arealight;

/* ПОСТАНОВКА ПЛОЩАДКИ — ОДНА ФУНКЦИЯ НА ВСЕХ ПОТРЕБИТЕЛЕЙ (§524, пункт 1). Была
 * встроена в ветвь `lit`; ветви замера грязи (4г2) нужен ТОТ ЖЕ источник, а
 * вторая редакция того же правила означала бы, что замер меряет разницу формул.
 * Тело перенесено дословно; выход ветви `lit` обязан остаться побитово тем же.
 * `eyeg` — глаз в единицах ячеек поля (точка, заведомо лежащая В ПОЛОСТИ). */
static void hall_light(arealight *AL, const opyr *P, const frame *fr, const double lo[3],
                       const double hi[3], const double eyeg[3], int verbose) {
  for (int k = 0; k < 3; k++)
    AL->c[k] = 0.5 * (lo[k] + hi[k]);
  /* ПОТОЛОК ПОЛОСТИ, А НЕ ВЕРХ ГАБАРИТА (§472). Прежнее «на 10 см ниже hi[1]»
   * верно только для сцены-оболочки без толщины: у замкнутой сцены с толстыми
   * стенами там материал, и лампа оказывается замурованной (замерено:
   * освещённых ячеек 0, кадр чёрный целиком).
   * СЧИТАТЬ НАДО ОТ ТОЧКИ, ЗАВЕДОМО ЛЕЖАЩЕЙ В ПОЛОСТИ, А НЕ СВЕРХУ: спуск
   * сверху упирается не в потолок полости, а в ПУСТОТУ ВНУТРИ ПЛИТЫ — занятость
   * метит ячейки, ЗАДЕТЫЕ ТРЕУГОЛЬНИКАМИ, а внутренность сплошного тела
   * треугольников не содержит. Замерено: правило сверху дало потолок на 2.385 м
   * при настоящем 2.2 м. Точка внутри полости у нас есть по построению — камера. */
  int32_t ye = (int32_t)eyeg[1], xe = (int32_t)eyeg[0], ze = (int32_t)eyeg[2];
  if (ye < 0) ye = 0;
  if (ye >= fr->n) ye = fr->n - 1;
  if (xe < 0) xe = 0;
  if (xe >= fr->n) xe = fr->n - 1;
  if (ze < 0) ze = 0;
  if (ze >= fr->n) ze = fr->n - 1;
  int32_t ytop = -1;
  for (int32_t y = ye; y < fr->n; y++)
    if (hz_occ_get(P->b[fr->lev], hz_occ_index(fr->n, xe, y, ze)) != 0) {
      ytop = y - 1;
      break;
    }
  /* ОТСТУП 0.10 М ПОД ПОТОЛКОМ. Прежний комментарий здесь объяснял решётку пятен
   * на потолке «артефактом выборки» и объявлял её вылеченной подъёмом панели.
   * ОБА УТВЕРЖДЕНИЯ БЫЛИ НЕВЕРНЫ (А958): решётка была точным изображением
   * `N×N` ТОЧЕЧНЫХ источников, потому что в `g` не было ни площади, ни косинуса
   * излучателя; панель при этом так и осталась висеть в 10 см. Теперь площадка
   * есть площадка и односторонняя, потолок над ней не освещён по построению, а
   * сам отступ безвреден и оставлен как есть. */
  AL->c[1] = ytop >= 0 ? fr->org[1] + ((double)ytop + 0.5) * fr->h - 0.10 : hi[1] - 0.10;
  double half = 0.25 * (hi[0] - lo[0]);
  AL->u[0] = half;
  AL->u[1] = 0.0;
  AL->u[2] = 0.0;
  AL->v[0] = 0.0;
  AL->v[1] = 0.0;
  AL->v[2] = 0.25 * (hi[2] - lo[2]);
  /* ОКРАШЕН СОЗНАТЕЛЬНО (А757): три одинаковых числа выдавать за RGB нельзя. */
  AL->rgb[0] = 1.00;
  AL->rgb[1] = 0.92;
  AL->rgb[2] = 0.78;
  if (verbose)
    printf("   ИСТОЧНИК: потолок полости на y = %.3f м, площадка на %.3f м, ПЛОЩАДЬ %.4f м² "
           "(%.3f × %.3f)\n",
           ytop >= 0 ? fr->org[1] + ((double)ytop + 0.5) * fr->h : hi[1], AL->c[1],
           4.0 * half * AL->v[2], 2.0 * half, 2.0 * AL->v[2]);
}

/* Затенён ли путь от точки к точке. Марш по ЗАНЯТОСТИ мелкого уровня с шагом в
 * пол-ячейки: занятая ячейка на пути — заслон. Концы исключаются, иначе сама
 * поверхность закрывала бы себя. */
static int shadowed(const opyr *P, const frame *fr, const double a[3], const double b[3],
                    double stepfrac) {
  double d[3], len = 0.0;
  for (int k = 0; k < 3; k++) {
    d[k] = b[k] - a[k];
    len += d[k] * d[k];
  }
  len = sqrt(len);
  if (!(len > 0.0)) return 0;
  double step = fr->h * stepfrac;
  int ns = (int)(len / step);
  if (ns > 16384) ns = 16384;
  /* ИСКЛЮЧЕНИЕ НАЧАЛА — ПО РАССТОЯНИЮ, А НЕ ПО ЧИСЛУ ОБРАЗЦОВ (найдено А776).
   * Прежде пропускались первые ДВА образца: при шаге в полячейки это отступ в
   * ЦЕЛУЮ ячейку, а при шаге в четверть — уже в половину. Тогда мельчение шага
   * начинало ловить СОБСТВЕННУЮ ячейку поверхности как заслон, и треть
   * освещённых ячеек пропадала — проверка сходимости мерила самозатенение, а не
   * заслоны. Отступ `1.5` ячейки назван числом: он выводит луч за пределы своей
   * ячейки и её соседа по грани, и от шага марша больше не зависит. */
  double skip = 1.5 * fr->h;
  for (int i = 0; i < ns - 1; i++) {
    double t = (double)i / (double)ns;
    if (t * len < skip) continue;
    int32_t c[3];
    int ok = 1;
    for (int k = 0; k < 3; k++) {
      double w = a[k] + d[k] * t;
      double f = floor((w - fr->org[k]) / fr->h);
      if (!(f >= 0.0) || !(f < (double)fr->n)) ok = 0;
      c[k] = ok ? (int32_t)f : 0;
    }
    if (!ok) continue;
    if (hz_occ_get(P->b[fr->lev], hz_occ_index(fr->n, c[0], c[1], c[2]))) return 1;
  }
  return 0;
}

/* ЗАТЕНЕНИЕ С ПОДЪЁМОМ В ГРУБЫЕ ВЕТКИ (А784). Плоский марш идёт шагом в
 * полячейки и платит по РАССТОЯНИЮ; здесь на каждом шаге ищется САМЫЙ КРУПНЫЙ
 * пустой узел, накрывающий текущую точку, и он пересекается ЦЕЛИКОМ. Большой
 * пустой куб стоит один шаг, а не свою сторону в ячейках.
 *
 * ПРЕДИКАТ ТОТ ЖЕ, что у плоского марша: заслон — занятая ячейка мелкого уровня.
 * Значит множество затенённых обязано СОВПАСТЬ, и это приёмка, а не пожелание:
 * пропуск пустоты ТОЧЕН по построению (пустой узел не содержит геометрии вовсе),
 * и любое расхождение означает ошибку в подъёме, а не приближение. */
static int shadowed_h(const opyr *P, const frame *fr, const double a[3], const double b[3],
                      int64_t *nstep) {
  double d[3], len = 0.0;
  for (int k = 0; k < 3; k++) {
    d[k] = b[k] - a[k];
    len += d[k] * d[k];
  }
  len = sqrt(len);
  if (!(len > 0.0)) return 0;
  for (int k = 0; k < 3; k++)
    d[k] /= len;
  double skip = 1.5 * fr->h; /* тот же отступ начала, что и у плоского (А776) */
  /* ОТСТУП И В КОНЦЕ — НАЙДЕНО §617. Плоский марш исключает последние образцы
   * (`i < ns − 1`), а здесь этого не было, и марш упирался в СОБСТВЕННУЮ ячейку
   * цели. Пока единственным потребителем был прямой свет, пробел не проявлялся:
   * там цель — солнце, точка далеко ВНЕ геометрии, и конец луча ни на что не
   * попадает. У связи «поверхность -> поверхность» цель ЛЕЖИТ НА ПОВЕРХНОСТИ, и
   * каждая связь объявлялась заслонённой. Замерено: `Σ b_1` падало в `4.8` раза.
   * Довод тот же, что у А776, и симметричный ему: луч обязан выйти из своей
   * ячейки и её соседа по грани НА ОБОИХ концах. */
  double tend = len - skip;
  double t = skip;
  while (t < tend) {
    int32_t c[3];
    int ok = 1;
    for (int k = 0; k < 3; k++) {
      double w = a[k] + d[k] * t;
      double f = floor((w - fr->org[k]) / fr->h);
      if (!(f >= 0.0) || !(f < (double)fr->n)) ok = 0;
      c[k] = ok ? (int32_t)f : 0;
    }
    if (!ok) return 0;
    if (nstep != NULL) (*nstep)++;
    if (hz_occ_get(P->b[fr->lev], hz_occ_index(fr->n, c[0], c[1], c[2]))) return 1;
    /* ПОДЪЁМ: самый крупный ПУСТОЙ узел, накрывающий точку. Пока предок пуст —
     * поднимаемся; шаг равен стороне найденного узла. */
    int l = fr->lev;
    while (l > 0) {
      int32_t nl = (int32_t)1 << (l - 1);
      if (hz_occ_get(P->b[l - 1],
                     hz_occ_index(nl, c[0] >> (fr->lev - l + 1), c[1] >> (fr->lev - l + 1),
                                  c[2] >> (fr->lev - l + 1))))
        break;
      l--;
    }
    /* ШАГ — ДО ВЫХОДА ИЗ УЗЛА, А НЕ НА ЕГО СТОРОНУ. Первая редакция прибавляла
     * сторону узла от ТЕКУЩЕЙ точки и потому перескакивала за его дальнюю грань:
     * `958` ячеек разошлись с плоским маршем — заслон сразу за узлом
     * проглатывался. Здесь считается расстояние до ближайшей из трёх дальних
     * граней (обычная плитовая проба), и оно всегда меньше стороны. */
    int sh2 = fr->lev - l;
    double side = fr->h * (double)((int32_t)1 << sh2);
    double texit = 1e300;
    for (int k = 0; k < 3; k++) {
      if (!(fabs(d[k]) > 1e-300)) continue;
      int32_t nodelo = (c[k] >> sh2) << sh2;
      double lo0 = fr->org[k] + (double)nodelo * fr->h;
      double bnd = d[k] > 0.0 ? lo0 + side : lo0;
      double tk = (bnd - a[k]) / d[k];
      if (tk > t && tk < texit) texit = tk;
    }
    /* Отступ в четверть ячейки выводит точку ЗА грань: без него следующая проба
     * попадала бы ровно на границу и топталась на месте. */
    t = (texit < 1e299 ? texit : t + side) + fr->h * 0.25;
  }
  return 0;
}

/* Прямая облучённость ячейки среза по трём каналам.
 *
 * ЧИСЛО ПРОБ ПЛОЩАДКИ ВЫВЕДЕНО ИЗ РАЗРЯДНОСТИ ВЫВОДА, А НЕ ПОДОБРАНО (правка
 * 08-11 по замечанию пользователя «рваные тени»). Прежние ЧЕТЫРЕ пробы стояли
 * по углам площадки, то есть полутень имела ровно ПЯТЬ ступеней — и они видны
 * на картинке как ступеньки, а не как градиент. Сетка `N×N` даёт `N²+1`
 * ступеней; чтобы ступенька была ниже кванта восьмибитного вывода (1/255),
 * нужно `N² >= 255`, то есть `N >= 16`. Это дорого (цена прямого света линейна
 * по пробам), поэтому берётся `N = 8`: `65` ступеней, квант `1.5 %` — на глаз
 * уже градиент, а цена растёт вчетверо, а не в шестнадцать. Число названо
 * вместе с ценой и с тем, чего оно НЕ ДАЁТ: полной гладкости полутени.
 *
 * ПРОБЫ СТОЯТ В ЦЕНТРАХ ЯЧЕЕК сетки, а не по углам: угловая сетка смещена
 * наружу и переоценивает полутень у самого края площадки. */
#define HZ_LIGHT_NS 8
#define HZ_LIGHT_SAMPLES (HZ_LIGHT_NS * HZ_LIGHT_NS)

/* Смещение `s`-й пробы в долях полуоси, в `[-1, 1]`. */
static double lsamp(int i) {
  return (2.0 * ((double)i + 0.5) / (double)HZ_LIGHT_NS) - 1.0;
}

/* --- Ф10' (§541): ПЛОЩАДКА КАК ПЛОЩАДКА, А НЕ КАК N ТОЧЕК ------------------ */

/* ЧТО ЗДЕСЬ БЫЛО НЕ ТАК ДЕСЯТЬ ДНЕЙ (А958). Стояло `g = cos θ_r / r² / N` — без
 * площади, без косинуса излучателя и без дискового члена. Это не квадратура
 * интеграла по площадке, а сумма `N` ТОЧЕЧНЫХ источников силы `L_e/N`, то есть
 * ровно то, что §261 запрещает («точечного не заводить ВОВСЕ, он сингулярность
 * модели цены») и что §212 уже ловил 08-04 в другом месте. Решётка `8×8`
 * светлых пятен на потолке во всех картинках проекта есть их точное
 * изображение.
 *
 * ЧТО СТОИТ ТЕПЕРЬ. Проба несёт подплощадь `a = A/N`:
 *
 *     E_s = L_e · cos θ_s · cos θ_r · π·a / (π r² + a)
 *
 * Дальнее поле даёт `L_e cos θ_s cos θ_r · a/r²` — обычная квадратура. Ближнее
 * ограничено, и ограничивает его именно `+a`: это форм-фактор на ДИСК площади
 * `a`, а не на точку.
 *
 * ОДНА ФУНКЦИЯ НА ОБА ПУТИ (А963). `front_direct` и `front_sweep` считают одно и
 * то же; будь множитель написан в каждой отдельно, они разошлись бы, как
 * разошлись соглашения об альбедо в §531. Здесь им разойтись нечем. */
static int g_ptlight = 0; /* НК §541: вернуть прежнюю точечную формулу */

/* СОЛНЦЕ (§570): направленный источник вместо площадки. Наружной сцене площадка
 * под «потолком полости» не годится по построению — у двора это правило дало
 * светящийся потолок `465.7` м² размером с половину сцены, и `613.7` мс прямого
 * света были ценой именно его, а не алгоритма.
 * У направленного источника углового размера нет, поэтому проба ОДНА, а не
 * `64`: облучённость есть `L_e·cos θ_r`, без `1/r²` и без площади. Тень —
 * тот же марш, но к точке, отнесённой на габарит сцены. */
static int g_sun = 0;
static double g_sundir[3] = {-0.3, -1.0, -0.2}; /* §2: город, ω = normalize(−0.3,−1,−0.2) */

static double alight_g(const arealight *L, const double w[3], double r2, double cosr) {
  if (g_ptlight) return cosr / r2 / (double)HZ_LIGHT_SAMPLES;
  double r = sqrt(r2);
  /* Косинус ИЗЛУЧАТЕЛЯ: нормаль площадки есть `u × v`, направление на приёмник
   * противоположно `w` (тот идёт от приёмника к источнику). */
  double ns[3] = {L->u[1] * L->v[2] - L->u[2] * L->v[1], L->u[2] * L->v[0] - L->u[0] * L->v[2],
                  L->u[0] * L->v[1] - L->u[1] * L->v[0]};
  double nl = sqrt(ns[0] * ns[0] + ns[1] * ns[1] + ns[2] * ns[2]);
  if (!(nl > 0.0)) return 0.0;
  /* ПЛОЩАДКА ОДНОСТОРОННЯЯ, и это не выбор, а то, что §472 уже объявил физикой:
   * «у панели заподлицо потолок получает cos ~ 0 и не светится вовсе». Намерение
   * было записано, а косинуса излучателя в коде не было вовсе — потому потолок и
   * светился. Двусторонняя площадка светила бы вверх, В МАТЕРИАЛ потолка. */
  double coss = -(w[0] * ns[0] + w[1] * ns[1] + w[2] * ns[2]) / (r * nl);
  if (!(coss > 0.0)) return 0.0;
  /* `A = 4|u||v|`, потому что пробы бегут по `(−1, 1)` в обеих полуосях. */
  double a = 4.0 * nl / (double)HZ_LIGHT_SAMPLES;
  return coss * cosr * 3.14159265358979323846 * a / (3.14159265358979323846 * r2 + a);
}

static void front_direct(const hz_dcslice *S, const frame *fr, const opyr *P, const arealight *L,
                         float *irr, double stepfrac, int hier, int64_t *nstep,
                         const hz_objmesh *A) {
  double su[HZ_LIGHT_SAMPLES], sv[HZ_LIGHT_SAMPLES];
  /* СОЛНЦЕ: одна проба вместо `64`, и она не на площадке, а «за горизонтом» —
   * точка, отнесённая вдоль `−ω` на габарит сцены. Заслон считается тем же
   * маршем, что и у площадки, поэтому сравнение идёт схема в схему. */
  int nsmp = g_sun ? 1 : HZ_LIGHT_SAMPLES;
  double sunfar = 0.0;
  for (int a = 0; a < 3; a++)
    sunfar += (double)fr->n * fr->h * (double)fr->n * fr->h;
  sunfar = sqrt(sunfar);
  for (int a = 0; a < HZ_LIGHT_NS; a++)
    for (int b = 0; b < HZ_LIGHT_NS; b++) {
      su[a * HZ_LIGHT_NS + b] = lsamp(a);
      sv[a * HZ_LIGHT_NS + b] = lsamp(b);
    }
  /* Р1 (§581): ДЕЛЕНИЕ ПО ПРИЁМНИКАМ. Записи независимы (`irr[i]`), чтения
   * общие и неизменные, порядок сложения ВНУТРИ приёмника не меняется — значит
   * ответ побитово тот же, и это проверяется приёмкой, а не предполагается. */
#pragma omp parallel for schedule(dynamic, 256) if (!g_omp1)
  for (int32_t i = 0; i < S->n; i++) {
    double p[3], n[3];
    hz_slice_vertex(S, i, p);
    for (int k = 0; k < 3; k++)
      p[k] = fr->org[k] + p[k] * fr->h;
    hz_slice_normal(S, i, n);
    double acc[3] = {0, 0, 0};
    for (int s = 0; s < nsmp; s++) {
      double q[3], w[3], r2 = 0.0;
      for (int k = 0; k < 3; k++) {
        q[k] = g_sun ? p[k] - g_sundir[k] * sunfar : L->c[k] + L->u[k] * su[s] + L->v[k] * sv[s];
        w[k] = q[k] - p[k];
        r2 += w[k] * w[k];
      }
      if (!(r2 > 0.0)) continue;
      double r = sqrt(r2);
      double cosr = (w[0] * n[0] + w[1] * n[1] + w[2] * n[2]) / r;
      /* Поверхность односторонняя, и знака у неё нет (Р7): освещённой считается
       * та сторона, к которой нормаль обращена, а модуль брать нельзя — иначе
       * стена светилась бы с обратной стороны. */
      if (!(cosr > 0.0)) continue;
      /* `stepfrac < 0` — БЕЗ ЗАТЕНЕНИЯ ВОВСЕ: это знаменатель доли
       * открытости, а не режим рендера. */
      if (stepfrac > 0.0) {
        int sh = hier ? shadowed_h(P, fr, p, q, nstep) : shadowed(P, fr, p, q, stepfrac);
        if (sh) continue;
      }
      /* У направленного источника ни площади, ни `1/r²`: облучённость есть
       * `L_e·cos θ_r`, и это не упрощение, а определение. */
      double g = g_sun ? cosr : alight_g(L, w, r2, cosr);
      for (int k = 0; k < 3; k++)
        acc[k] += L->rgb[k] * g * alb(A, S->c[i].mat, k);
    }
    for (int k = 0; k < 3; k++)
      irr[3 * (size_t)i + (size_t)k] = (float)acc[k];
  }
}

/* §520: допуск на разброс нормалей — отдельно от углового. */
static double g_hspread = 0.25;

/* НЕГАТИВНЫЙ КОНТРОЛЬ §520: все нормали в ОДНУ корзину — усреднение через
 * складку, как было до §519. Энергия обязана уехать вдвое. */
static int g_onebin = 0;

/* --- Ф6' (§516): ДЕРЕВО ИЗЛУЧАТЕЛЕЙ ---------------------------------------- */

/* ЗАЧЕМ. Гатер «каждая с каждой» стоит `N²` пар (`3.7e8` на комнате) и занимает
 * 94 % кадра. Но приёмнику не нужен каждый излучатель по отдельности: дальняя
 * стена светит как ОДНА площадка. Огрубление идёт по ВЗАИМНОМУ расстоянию, а не
 * по расстоянию до глаза — §510 проверил обратное и был наказан замером
 * (энергия упала вдвое при падении числа излучателей на треть).
 *
 * УЗЕЛ НЕСЁТ ПОТОК, А НЕ ЯРКОСТЬ: `flux = Σ радианс·площадь` по детям. Тогда
 * сложение точное, а яркость восстанавливается делением на площадь там, где
 * нужна. Складывать яркости было бы неверно — они не аддитивны. */
/* КОРЗИНЫ ПО НАПРАВЛЕНИЮ (§519). Узел агрегируется не как одно целое, а по
 * СЕМЕЙСТВАМ поверхностей: шесть корзин по знаковой главной оси нормали.
 * Тогда столешница (+y) и ножка (±x) — разные семейства, у каждого разброс
 * ничтожен, и оба агрегируются точно. Стена с толщиной: внутренняя сторона в
 * одной корзине, наружная в противоположной, и они не смешиваются.
 * ШЕСТЬ, А НЕ БОЛЬШЕ: внутри корзины нормали лежат в конусе 90°, то есть
 * разброс ограничен `1 − cos 45° = 0.29` сверху ПО ПОСТРОЕНИЮ, а у плоскости он
 * ноль. Больше корзин — точнее конус, но дороже узел; число проверяется
 * замером, а не назначается навсегда. */
typedef struct {
  float c[3], n[3], flux[3];
  float area, nsum, rad;
  signed char q; /* номер корзины: хранятся только НЕПУСТЫЕ */
} ebin;

/* §529/§531 (ПРОВЕРКА ПРИБОРА ДО Ф8', правило А824). АЛЬБЕДО ИЗЛУЧАТЕЛЯ ВХОДИЛО
 * ДВАЖДЫ, И ЭТО ИСПРАВЛЕНО ЗДЕСЬ. `irr` уже содержит `alb` ПРИЁМНИКА (кладут
 * `front_direct` и `front_sweep`), поэтому отскок, умножавший `irr[j]` ещё раз
 * на `alb(j)`, нёс `rho_j^2` вместо `rho_j`. Радиосити требует `B_j = rho_j E_j`
 * ровно один раз.
 * ЗАМЕРЕНО (§531): отношение `E_ind/E_dir` на комнате `0.4755 -> 0.6601` у
 * прямого гатера и `0.4643 -> 0.6430` у иерархического — два НЕЗАВИСИМЫХ пути
 * сдвинулись в `1.388` и `1.385` раза, при обратном средневзвешенном альбедо
 * излучателей `1/0.7342 = 1.362`.
 * КЛЮЧ `emitalb2` ВОЗВРАЩАЕТ ПРЕЖНЕЕ (НЕВЕРНОЕ) ПОВЕДЕНИЕ — он нужен затем,
 * чтобы старые доклады оставались ВОСПРОИЗВОДИМЫМИ, а не просто объявленными
 * недействительными. Числа под ним цитировать нельзя. */
static int g_emitalb2 = 0;

/* УЗЕЛ — ЗАГОЛОВОК, А НЕ КОНТЕЙНЕР (§521, замечание пользователя «СКОЛЬКО?!»).
 * Первая редакция держала ШЕСТЬ корзин у КАЖДОГО узла и всё в `double`: 616 Б на
 * узел при том, что непустых корзин обычно одна-две. Это нарушение правила,
 * записанного мною же сутки назад (§498: считать байты на элемент ДО того, как
 * структура написана). Теперь корзины лежат общим массивом, у узла — диапазон;
 * полезная нагрузка во `float` (координаты и потоки ограничены по диапазону, а
 * относительная точность 1e-7 на порядок ниже допуска спуска). */
typedef struct {
  int32_t b0;    /* первая корзина в общем массиве */
  int32_t ch[8]; /* дети; -1 там, где пусто. НЕ подряд (§512). */
  signed char nb, nch;
} enode;

typedef struct {
  enode *e;
  /* ГРУППЫ ЛЕЖАТ ОБЩИМ МАССИВОМ, У УЗЛА — ДИАПАЗОН (§521, замечание
   * пользователя). Первая редакция держала ШЕСТЬ корзин у КАЖДОГО узла и всё в
   * `double` — 616 Б на узел при том, что непустых обычно одна-две. Это
   * нарушение правила, записанного сутки назад (§498). Теперь узел есть
   * ЗАГОЛОВОК: диапазон групп плюс дети. */
  ebin *b;
  int32_t nb, bcap;
  int32_t n, cap;
} etree;

static void etree_free(etree *T) {
  free(T->e);
  free(T->b);
  T->b = NULL;
  T->nb = T->bcap = 0;
  T->e = NULL;
  T->n = T->cap = 0;
}

static int32_t etree_alloc(etree *T, int32_t k) {
  if (T->n + k > T->cap) {
    int32_t nc = T->cap > 0 ? T->cap * 2 : 4096;
    if (nc < T->n + k) nc = T->n + k;
    enode *nn = realloc(T->e, (size_t)nc * sizeof *nn);
    if (nn == NULL) exit(1);
    T->e = nn;
    T->cap = nc;
  }
  int32_t b = T->n;
  T->n += k;
  return b;
}

/* Дерево строится ПРЯМО ПО СРЕЗУ, потому что срез лежит в МОРТОНОВОМ порядке
 * (dcslice.h): у любого узла его ячейки образуют НЕПРЕРЫВНЫЙ отрезок массива, и
 * разбиение на восьмерых детей есть разрезание отрезка по биту координаты. Ни
 * хеша, ни второго индекса не нужно — и это прямая выгода от решения хранить
 * срез массивом, а не деревом (§382). */
static int32_t ebin_alloc(etree *T, int k) {
  if (T->nb + k > T->bcap) {
    int32_t nc = T->bcap > 0 ? T->bcap * 2 : 8192;
    if (nc < T->nb + k) nc = T->nb + k;
    ebin *nn = realloc(T->b, (size_t)nc * sizeof *nn);
    if (nn == NULL) exit(1);
    T->b = nn;
    T->bcap = nc;
  }
  int32_t r = T->nb;
  T->nb += k;
  return r;
}

static int32_t etree_build(etree *T, const hz_dcslice *S, const frame *fr, const float *irr,
                           const hz_objmesh *m, int32_t a, int32_t b, int lvl, int lev) {
  int32_t me = etree_alloc(T, 1);
  {
    /* Шесть корзин набираются во ВРЕМЕННЫХ, а хранятся только непустые. */
    double c[6][3] = {{0}}, n[6][3] = {{0}}, fl[6][3] = {{0}}, ar[6] = {0};
    double lo[6][3], hi[6][3];
    for (int q = 0; q < 6; q++)
      for (int k = 0; k < 3; k++) {
        lo[q][k] = 1e300;
        hi[q][k] = -1e300;
      }
    for (int32_t i = a; i < b; i++) {
      double p[3], nn[3];
      hz_slice_vertex(S, i, p);
      for (int k = 0; k < 3; k++)
        p[k] = fr->org[k] + p[k] * fr->h;
      hz_slice_normal(S, i, nn);
      int ax = 0;
      for (int k = 1; k < 3; k++)
        if (fabs(nn[k]) > fabs(nn[ax])) ax = k;
      int q = g_onebin ? 0 : 2 * ax + (nn[ax] > 0.0 ? 1 : 0);
      double side = fr->h * (double)((int32_t)1 << (lev - (int)S->c[i].lvl));
      double aa = side * side;
      ar[q] += aa;
      for (int k = 0; k < 3; k++) {
        c[q][k] += p[k] * aa;
        n[q][k] += nn[k] * aa;
        if (p[k] < lo[q][k]) lo[q][k] = p[k];
        if (p[k] > hi[q][k]) hi[q][k] = p[k];
        fl[q][k] += (double)irr[3 * (size_t)i + (size_t)k] * aa *
                    (g_emitalb2 ? alb(m, S->c[i].mat, k) : 1.0);
      }
    }
    int nbq = 0;
    for (int q = 0; q < 6; q++)
      if (ar[q] > 0.0) nbq++;
    int32_t b0 = ebin_alloc(T, nbq);
    T->e[me].b0 = b0;
    T->e[me].nb = (signed char)nbq;
    T->e[me].nch = 0;
    for (int k = 0; k < 8; k++)
      T->e[me].ch[k] = -1;
    int j = 0;
    for (int q = 0; q < 6; q++) {
      if (!(ar[q] > 0.0)) continue;
      ebin *bb = &T->b[b0 + j++];
      memset(bb, 0, sizeof *bb);
      bb->q = (signed char)q;
      bb->area = (float)ar[q];
      double nl = 0.0;
      for (int k = 0; k < 3; k++) {
        bb->c[k] = (float)(c[q][k] / ar[q]);
        bb->flux[k] = (float)fl[q][k];
        double nk = n[q][k] / ar[q];
        nl += nk * nk;
      }
      nl = sqrt(nl);
      bb->nsum = (float)(nl * ar[q]);
      for (int k = 0; k < 3; k++)
        bb->n[k] = (float)(nl > 0.0 ? (n[q][k] / ar[q]) / nl : 0.0);
      double r2 = 0.0;
      for (int k = 0; k < 3; k++) {
        double d = 0.5 * (hi[q][k] - lo[q][k]);
        r2 += d * d;
      }
      bb->rad = (float)sqrt(r2);
    }
  }
  if (b - a <= 1 || lvl >= lev) return me;
  int sh = lev - lvl - 1;
  int32_t bnd[9];
  bnd[0] = a;
  int32_t cur = a;
  for (int k = 1; k <= 8; k++) {
    while (cur < b) {
      const hz_dccell *cc = &S->c[cur];
      int bit = (((cc->lo[0] >> sh) & 1) | (((cc->lo[1] >> sh) & 1) << 1) |
                 (((cc->lo[2] >> sh) & 1) << 2));
      if (bit >= k) break;
      cur++;
    }
    bnd[k] = cur;
  }
  int nc2 = 0;
  int32_t ch2[8];
  for (int k = 0; k < 8; k++) {
    if (bnd[k + 1] <= bnd[k]) continue;
    ch2[nc2++] = etree_build(T, S, fr, irr, m, bnd[k], bnd[k + 1], lvl + 1, lev);
  }
  for (int k = 0; k < nc2; k++)
    T->e[me].ch[k] = ch2[k];
  T->e[me].nch = (signed char)nc2;
  return me;
}

/* Спуск ведётся ПО ОДНОЙ КОРЗИНЕ; корзины узла лежат подряд, и пустых среди них
 * нет — поэтому обход читает ровно то, что нужно, и ни байтом больше. */
/* §600: НК1 — прежнее ТОЧЕЧНОЕ ядро; НК2 — поправка в 1000 раз (обязана уехать
 * в ДАЛЁКОМ поле, чем и проверяет, что П2 не слепа). */
static int g_gpoint = 0, g_gwide = 0;

/* ---- ХРАНИМЫЕ СВЯЗИ СБОРА (§613) -----------------------------------------
 *
 * ЗАЧЕМ. Видимость есть свойство ГЕОМЕТРИИ: ни яркость излучателей, ни BRDF, ни
 * камера в неё не входят. Замерено (§612): заслоны поднимают сбор с `11.1` до
 * `59.3` с. Платить это на КАЖДОМ отскоке — значит сделать многократный отскок
 * невозможным ровно тогда, когда он стал осмысленным.
 *
 * ЧТО ХРАНИТСЯ. На приёмник — список `(номер корзины, g)`, и только тот, что
 * ПРОШЁЛ заслон. Дальше `E_i = Σ g·flux[корзина]`, обхода дерева нет вовсе.
 *
 * ПОЧЕМУ БУФЕРЫ ПОПОТОЧНЫЕ. Приёмники делятся по потокам, и каждый поток пишет
 * свой буфер подряд. Атомарных операций нет, а ПОРЯДОК СЛОЖЕНИЯ ВНУТРИ приёмника
 * не меняется — потому первый отскок обязан совпасть с необходимым путём
 * ПОБИТОВО, и это проверяется, а не предполагается.
 *
 * НОМЕРА КОРЗИН УСТОЙЧИВЫ МЕЖДУ ОТСКОКАМИ, потому что `etree_build`
 * детерминирована при том же срезе: дерево перестраивается каждый отскок
 * (`1.05` с), но раскладка корзин повторяется. Переписывать построение ради
 * обновления одних потоков — отдельная работа, и `1.05` с против `59` не та
 * статья, ради которой стоит рисковать. */
typedef struct {
  int32_t bin;
  float g;
} glink;

typedef struct {
  glink **buf;  /* [поток] растущий массив */
  int64_t *n;   /* сколько занято */
  int64_t *cap; /* сколько выделено */
  int64_t *off; /* [приёмник] смещение в буфере своего потока */
  int32_t *cnt; /* [приёмник] сколько связей */
  int8_t *th;   /* [приёмник] чей это буфер */
  int nth;
  int32_t nrecv;
} linkcache;

static void lc_free(linkcache *L) {
  if (L->buf != NULL)
    for (int t = 0; t < L->nth; t++)
      free(L->buf[t]);
  free(L->buf);
  free(L->n);
  free(L->cap);
  free(L->off);
  free(L->cnt);
  free(L->th);
  memset(L, 0, sizeof *L);
}

static int lc_init(linkcache *L, int nth, int32_t nrecv) {
  memset(L, 0, sizeof *L);
  L->nth = nth;
  L->nrecv = nrecv;
  L->buf = calloc((size_t)nth, sizeof *L->buf);
  L->n = calloc((size_t)nth, sizeof *L->n);
  L->cap = calloc((size_t)nth, sizeof *L->cap);
  L->off = calloc((size_t)nrecv, sizeof *L->off);
  L->cnt = calloc((size_t)nrecv, sizeof *L->cnt);
  L->th = calloc((size_t)nrecv, sizeof *L->th);
  if (L->buf == NULL || L->n == NULL || L->cap == NULL || L->off == NULL || L->cnt == NULL ||
      L->th == NULL) {
    lc_free(L);
    return 0;
  }
  return 1;
}
/* §611: НК1 — без нормировки суммы формфакторов. */
static int g_gnonorm = 0;
/* §613, НК2: не хранить связи — обходить дерево на каждом отскоке. Без этого
 * контроля «стало быстрее» неотличимо от «стало меньше работы». */
static int g_nolinkcache = 0;
/* §616: прежний ПЛОСКИЙ марш заслона — эталон и негативный контроль. */
static int g_gflatvis = 0;
/* §621: `hangle` — прежний УГЛОВОЙ критерий дробления (негативный контроль и
 * путь отката); `hcontrib` — доля от масштаба сцены для отбора по вкладу. */
static int g_hangle = 0;
static double g_hcontrib = 1e-3;
/* §618, НК: не читать `Ke`. Без солнца `Σ b_0` обязана стать РОВНО НОЛЬ. */
static int g_nokemit = 0;
/* §618: `b0` ТОЛЬКО из `Ke` — ни солнца, ни площадки. */
static int g_keonly = 0;
/* §620: шагов иерархического марша заслона — чтобы делить время замером, а не
 * догадкой. Сводится в фиксированном порядке из попоточных. */
static int64_t g_march_steps = 0;

static void hgather_rec(const etree *T, int32_t ni, int q, const double pi[3], const double ni_[3],
                        double eps, double tau, double rrecv, const opyr *P, const frame *fr,
                        int vis, double out[3], int64_t *nlink, double *sthru, double *sall,
                        int64_t *nnear, double *ffsum, linkcache *lc, int lcth, int64_t *nmarch) {
  const enode *e = &T->e[ni];
  const ebin *bb = NULL;
  for (int k = 0; k < e->nb; k++)
    if (T->b[e->b0 + k].q == q) {
      bb = &T->b[e->b0 + k];
      break;
    }
  if (bb == NULL) return;
  /* нулевая инициализация — гейт (§794-класс): gcc-analyzer на пути из
   * efc_run теряет запись w[k] в цикле ниже; путь ложный, семантика та же */
  double w[3] = {0, 0, 0}, r2 = 0.0;
  for (int k = 0; k < 3; k++) {
    w[k] = (double)bb->c[k] - pi[k];
    r2 += w[k] * w[k];
  }
  if (!(r2 > 0.0)) return;
  /* Корень — только на принятой связи: критерий спуска через КВАДРАТЫ. */
  double d4 = 4.0 * (double)bb->rad * (double)bb->rad;
  double spread = 1.0 - (double)bb->nsum / (double)bb->area;
  /* §621: ДРОБЛЕНИЕ ПО ВКЛАДУ, А НЕ ПО УГЛОВОМУ РАЗМЕРУ.
   *
   * ЗАМЕРЕНО (§620): по угловому критерию величина НЕ СХОДИТСЯ — приращения при
   * измельчении `eps` идут `+7`, `+25`, `+14`, `+13 %` и не убывают, при том что
   * цена строго линейна по числу связей. Признак отбирал связи не по тому, что
   * важно: `d4 > eps·r2` смотрит на УГЛОВОЙ РАЗМЕР излучателя и НЕ смотрит на
   * его ЯРКОСТЬ, поэтому тусклая близкая площадка дробилась как яркая.
   *
   * ОЦЕНКА ВКЛАДА — КОНСЕРВАТИВНАЯ, СВЕРХУ: косинусы не больше единицы, значит
   * вклад не больше `flux_max/(π r² + A)`. Так критерий никогда не НЕДОдробит, а
   * косинусы для решения не считаются — оно остаётся дешевле самого сложения.
   * ПОРОГ — ДОЛЯ ОТ МАСШТАБА СЦЕНЫ (`tau`, см. `gather_run`), а не абсолютное
   * число: абсолютных порогов проект не допускает. */
  int split;
  if (g_hangle) {
    split = d4 > eps * r2 && 2.0 * (double)bb->rad > rrecv;
  } else {
    double fmax = (double)bb->flux[0];
    if ((double)bb->flux[1] > fmax) fmax = (double)bb->flux[1];
    if ((double)bb->flux[2] > fmax) fmax = (double)bb->flux[2];
    double cmax = fmax / (3.14159265358979323846 * r2 + (double)bb->area);
    split = cmax > tau && 2.0 * (double)bb->rad > rrecv;
  }
  if (e->nch > 0 && (spread > g_hspread || split)) {
    for (int k = 0; k < e->nch; k++)
      hgather_rec(T, e->ch[k], q, pi, ni_, eps, tau, rrecv, P, fr, vis, out, nlink, sthru, sall,
                  nnear, ffsum, lc, lcth, nmarch);
    return;
  }
  double di = w[0] * ni_[0] + w[1] * ni_[1] + w[2] * ni_[2];
  double dj = -(w[0] * (double)bb->n[0] + w[1] * (double)bb->n[1] + w[2] * (double)bb->n[2]);
  if (!(di > 0.0) || !(dj > 0.0)) return; /* отсев БЕЗ корня */
  double rr = sqrt(r2);
  (*nlink)++;
  /* ОГРАНИЧЕННОЕ ЯДРО ВМЕСТО ТОЧЕЧНОГО (§600). Точечное `cos·cos/(π r²)`
   * расходится как `1/r²`, и условие дробления выше при КРУПНОМ приёмнике
   * ложно — связь принимается вплотную. Замерено (§599): на Bistro отношение
   * косвенного к прямому вышло `1 712 335` вместо порядка единицы.
   *
   * ФИЗИКА, КОТОРУЮ ТОЧЕЧНОЕ ЯДРО НЕ ЗНАЕТ: облучённость от полусферы яркости
   * `B` равна `B`, и больше не бывает НИКОГДА. Прибавка `A_j` в знаменателе
   * даёт `E = B·cos·cos·A/(π r² + A) ≤ B` — граница по ПОСТРОЕНИЮ, а не по
   * порогу. На далёких связях `A ≪ π r²`, поправка исчезает, и далёкое поле
   * обязано остаться прежним (это проверяется П2 §600 на комнате).
   *
   * `bb->area` — СУММА площадей ячеек корзины, а не площадь пятна, которое они
   * занимают (А1034). Граница `E ≤ B` верна при любом соотношении, но для
   * рассеянной корзины `π·rad²` ограничивало бы сильнее. Это НЕ ПРОВЕРЕННАЯ
   * альтернатива и записанный замерный долг, а не молчаливый выбор. */
  double soft = g_gpoint ? 0.0 : (double)bb->area * (g_gwide ? 1000.0 : 1.0);
  double g = (di / rr) * (dj / rr) / (3.14159265358979323846 * r2 + soft);
  /* А1036: доля связей, где поправка ЗНАЧИМА. Счётчик ПОПОТОЧНЫЙ, как и
   * остальные: гонка в горячем цикле недопустима, а редукция OpenMP отдала бы
   * порядок планировщику. */
  if (nnear != NULL && soft > 0.01 * 3.14159265358979323846 * r2) (*nnear)++;
  int blocked = 0;
  double cw[3] = {(double)bb->c[0], (double)bb->c[1], (double)bb->c[2]};
  /* §616: ИЕРАРХИЧЕСКИЙ МАРШ, А НЕ ПЛОСКИЙ. Плоский идёт шагом в полячейки и
   * платит ПО РАССТОЯНИЮ: на Bistro (сцена `512` ячеек) длинная связь стоила
   * порядка тысячи проб. `shadowed_h` пропускает крупные пустые узлы ЦЕЛИКОМ и
   * заявлен тем же предикатом (§426, А784) — прямой свет считает им давно.
   * Замечание пользователя 08-13: «все заслоны уже есть в архитектуре». */
  if (vis) blocked = g_gflatvis ? shadowed(P, fr, pi, cw, 0.5) : shadowed_h(P, fr, pi, cw, nmarch);
  /* §609: СУММА ФОРМФАКТОРОВ. `F_ij = cos_i cos_j A/(π r² + A)` есть в точности
   * `g · A_j`. Физика: `Σ_j F_ij ≤ 1` у полностью замкнутой точки и СТРОГО МЕНЬШЕ
   * у открытой. Всякий приёмник выше единицы — доказательство завышения.
   * СЧИТАЕТСЯ ПОСЛЕ ЗАСЛОНА И ТОЛЬКО ПО ПРОШЕДШИМ СВЯЗЯМ (§612): иначе замер
   * слеп к тому, чинят ли заслоны сумму, — а именно это и надо проверить. */
  if (ffsum != NULL && !(blocked && vis)) *ffsum += g * (double)bb->area;
  /* §613: связь ЗАПОМИНАЕТСЯ, если прошла заслон. Индекс корзины устойчив
   * между отскоками (см. шапку `linkcache`). */
  if (lc != NULL && !(blocked && vis)) {
    if (lc->n[lcth] >= lc->cap[lcth]) {
      int64_t nc = lc->cap[lcth] > 0 ? lc->cap[lcth] * 2 : 1 << 16;
      glink *nb2 = realloc(lc->buf[lcth], (size_t)nc * sizeof *nb2);
      if (nb2 == NULL) exit(1);
      lc->buf[lcth] = nb2;
      lc->cap[lcth] = nc;
    }
    lc->buf[lcth][lc->n[lcth]].bin = (int32_t)(bb - T->b);
    lc->buf[lcth][lc->n[lcth]].g = (float)g;
    lc->n[lcth]++;
  }
  for (int k = 0; k < 3; k++) {
    double v = (double)bb->flux[k] * g;
    if (sall != NULL) {
      *sall += v;
      if (blocked && sthru != NULL) *sthru += v;
    }
    if (!(blocked && vis)) out[k] += v;
  }
}

/* --- 3ж. §812: БЛИКИ — АНАЛИТИЧЕСКИЙ ФОНГ-ЧЛЕН В РАСТРЕ ---------------------
 *
 * Раскол §625, измеренный планом §812: свип на канон-ординатах ND 32 несёт
 * лепестки шире ~0.9 рад (§538), Bistro-глянец живёт в 0.1…1 рад — ОБРАЗ в
 * свипе недостижим (ND 64…256 — такт ×2…×8, линейность А1358), аналитике
 * число ординат не нужно. Альбедо свипа идёт ОТ kd (§739-скоринг), Ks-доля
 * сегодня ПРОПАДАЕТ — член её добавляет, двойного счёта нет; глянцевое
 * косвенное (Ks²) — названная потеря.
 *
 * Машина: etree излучателей §796 (efc_build, светящиеся элементы) + СВОЙ
 * обход с ядром Фонга (А1367). Ядро (А1366): L_o = f_r·E_bin,
 * f_r = ks·(ns+2)/(2π)·cosψ^ns, E_bin = Φ·cosθ_i·cosθ_lamp/(πr² + A) — та же
 * форма, что у диффузного члена (при ns = 0, ks = ρ ядро ТОЧНО ламбертово,
 * поэтому ns ≤ 0 пропускается — дубль диффуза, А1365). Лепесток режет обход:
 * вне cosψ < cext вклад нуль, конус ~sqrt(ln(1/ε)/ns) узок. Белая точка —
 * ДИФФУЗНАЯ (А1364): блик клипуется насыщением, а не тонет в перцентиле. */
static int g_xgloss = 0;
static etree *g_glet = NULL;
static const opyr *g_glP = NULL;
static const frame *g_glfr = NULL;
static double g_gltau = 0.0, g_glcext = 0.0, g_glbeta = 0.0;
static int g_glks0 = 0;     /* НК-а: ks := 0 — кадр обязан стать посимвольным */
static int g_glnorm812 = 0; /* НК-б: нормаль перевёрнута у каждого 100-го треугольника */
static float *g_glspec = NULL;
static int64_t g_glnpx1 = 0, g_glnpx10 = 0, g_glclip = 0, g_glmarch = 0, g_glprune = 0,
               g_glleaf = 0;
static double g_glsum = 0.0, g_glt812 = 0.0;

/* Потоковые счётчики обхода: в OMP-регионе у каждого потока свои, сведение —
 * критической секцией после цикла (гонка в горячем цикле недопустима, А1036). */
typedef struct {
  int64_t prune, leaf, march;
  double sum, t;
} glstat812;

static void gloss812_rec(const etree *T, int32_t ni, int q, const double pi[3],
                         const double nsurf[3], const double wr[3], double ks, double nsp, int vis,
                         double out[3], glstat812 *st) {
  const enode *e = &T->e[ni];
  const ebin *bb = NULL;
  for (int k = 0; k < e->nb; k++)
    if (T->b[e->b0 + k].q == q) {
      bb = &T->b[e->b0 + k];
      break;
    }
  if (bb == NULL) return;
  double w[3] = {0, 0, 0}, r2 = 0.0;
  for (int k = 0; k < 3; k++) {
    w[k] = (double)bb->c[k] - pi[k];
    r2 += w[k] * w[k];
  }
  if (!(r2 > 0.0)) return;
  double r = sqrt(r2);
  double ci8 = (w[0] * nsurf[0] + w[1] * nsurf[1] + w[2] * nsurf[2]) / r;
  double cj8 = -(w[0] * (double)bb->n[0] + w[1] * (double)bb->n[1] + w[2] * (double)bb->n[2]) / r;
  if (!(ci8 > 0.0) || !(cj8 > 0.0)) {
    st->prune++;
    return;
  }
  double cps = (w[0] * wr[0] + w[1] * wr[1] + w[2] * wr[2]) / r;
  /* ЛЕПЕСТОК: (cosψ)^ns < 1e-3 вне конуса; радиус корзины — консервативно */
  if (cps + (double)bb->rad / r < g_glcext) {
    st->prune++;
    return;
  }
  /* консервативная оценка вклада СВЕРХУ (урок §621: дробление по ВКЛАДУ, а не
   * по угловому размеру: косинусы ≤ 1, значит вклад ≤ ядро_макс·Φ/(πr²)) */
  double fmax = (double)bb->flux[0];
  if ((double)bb->flux[1] > fmax) fmax = (double)bb->flux[1];
  if ((double)bb->flux[2] > fmax) fmax = (double)bb->flux[2];
  double cmax =
      ks * (nsp + 2.0) / 6.28318530717958647692 * fmax * ci8 * cj8 / (3.14159265358979323846 * r2);
  if (!(cmax > g_gltau)) {
    st->prune++;
    return;
  }
  /* дробление: ядро меняется внутри корзины — пока угловой радиус больше
   * полуширины лепестка на полувысоте β = sqrt(2·ln2/ns)/2 */
  if (e->nch > 0 && (double)bb->rad / r > 0.5 * g_glbeta) {
    for (int k = 0; k < e->nch; k++)
      gloss812_rec(T, e->ch[k], q, pi, nsurf, wr, ks, nsp, vis, out, st);
    return;
  }
  st->leaf++;
  int blocked = 0;
  if (vis) {
    double cw[3] = {(double)bb->c[0], (double)bb->c[1], (double)bb->c[2]};
    blocked = shadowed_h(g_glP, g_glfr, pi, cw, &st->march);
  }
  if (blocked) return;
  double g8 = ci8 * cj8 / (3.14159265358979323846 * r2 + (double)bb->area);
  double fr8 = ks * (nsp + 2.0) / 6.28318530717958647692;
  if (cps < 0.0) cps = 0.0;
  double k8 = fr8 * pow(cps, nsp) * g8;
  for (int c = 0; c < 3; c++)
    out[c] += k8 * (double)bb->flux[c];
}

/* Спек-член одного пикселя: ks/ns — материала победившего треугольника
 * (смотрит ВЫЗЫВАЮЩИЙ, §667: материал живёт у элемента); dr — направление
 * ГЛАЗ→пиксель, исходящее к камере V = −dr; отражение R = 2(N·V)N − V. */
static void gloss812_pixel(const etree *T, const double wp[3], const double nn[3],
                           const double dr[3], double ks, double nsp, int vis, double out[3],
                           glstat812 *st) {
  out[0] = out[1] = out[2] = 0.0;
  if (!(ks > 0.0) || !(nsp > 0.0)) return; /* А1365: ns ≤ 0 — дубль диффуза */
  double vv[3];
  for (int c = 0; c < 3; c++)
    vv[c] = -dr[c];
  double dn = vv[0] * nn[0] + vv[1] * nn[1] + vv[2] * nn[2];
  if (!(dn > 0.0)) return; /* пиксель видит тыл — лоба нет */
  double wr[3];
  for (int c = 0; c < 3; c++)
    wr[c] = 2.0 * dn * nn[c] - vv[c];
  for (int q = 0; q < 6; q++)
    gloss812_rec(T, 0, q, wp, nn, wr, ks, nsp, vis, out, st);
}

/* САМОТЕСТ А1366: одна корзина ровно над точкой, лампа перпендикулярно,
 * отражение совпадает с направлением на лампу (cps = 1): полученное против
 * аналитики f_r·Φ/(πr² + A). Возвращает 0 при совпадении до 1e-12. */
static int gloss812_selftest(void) {
  etree T;
  memset(&T, 0, sizeof T);
  if (etree_alloc(&T, 1) < 0) return 1;
  if (ebin_alloc(&T, 1) < 0) {
    etree_free(&T);
    return 1;
  }
  T.e[0].b0 = 0;
  T.e[0].nb = 1;
  T.e[0].nch = 0;
  for (int k = 0; k < 8; k++)
    T.e[0].ch[k] = -1;
  ebin *bb = &T.b[0];
  memset(bb, 0, sizeof *bb);
  bb->q = 5; /* нормаль (0,0,−1): ось 2, минус */
  bb->c[0] = 0.0f;
  bb->c[1] = 0.0f;
  bb->c[2] = 1.0f;
  bb->n[2] = -1.0f;
  bb->area = 0.01f;
  bb->flux[0] = bb->flux[1] = bb->flux[2] = 2.0f;
  bb->rad = 0.05f;
  /* камера НАД точкой смотрит вниз: dr = глаз→пиксель = (0,0,−1); исходящее к
   * камере V = (0,0,1); отражение R = V — лампа прямо в лоб, cps = 1 */
  double pi[3] = {0, 0, 0}, nn[3] = {0, 0, 1}, dr[3] = {0, 0, -1};
  double out[3] = {0, 0, 0};
  glstat812 st;
  memset(&st, 0, sizeof st);
  double tau_sv = g_gltau, beta_sv = g_glbeta, cext_sv = g_glcext;
  g_gltau = 0.0;
  g_glbeta = 1e9;  /* не дробить (детей и так нет) */
  g_glcext = -1.0; /* лепесток не отсекает */
  double ks = 0.3, nsp = 8.0;
  gloss812_pixel(&T, pi, nn, dr, ks, nsp, 0, out, &st);
  g_gltau = tau_sv;
  g_glbeta = beta_sv;
  g_glcext = cext_sv;
  double r = 1.0;
  /* ожидание — из ТЕХ ЖЕ float-констант корзины (урок А1344: предсказывать от
   * носителя; расхождение 2e-10 было float-округлением 0.01f против 0.01) */
  double A8 = (double)bb->area, F8 = (double)bb->flux[0];
  double expect =
      ks * (nsp + 2.0) / 6.28318530717958647692 * F8 / (3.14159265358979323846 * r * r + A8);
  int bad = fabs(out[0] - expect) > 1e-12 * (fabs(expect) > 0.0 ? fabs(expect) : 1.0);
  printf("   §812 САМОТЕСТ ЯДРА: одна лампа над точкой, cps = 1: получено %.12e против "
         "аналитики %.12e — %s\n",
         out[0], expect, bad ? "РАСХОЖДЕНИЕ" : "ТОЖДЕСТВО");
  etree_free(&T);
  return bad;
}

/* --- 3з. СВИП ПОТОКОМ (Ш5, §383) ------------------------------------------- */

/* ЧТО ЗДЕСЬ ДЕЛАЕТСЯ И ЧЕМ ЭТО ОТЛИЧАЕТСЯ ОТ ЛУЧА. Луч платит `O(шаги)` за
 * КАЖДУЮ пару «ячейка — образец источника»: марш повторяется для каждой ячейки
 * заново, хотя соседние ячейки идут почти одним и тем же путём. Свип платит
 * `O(1)` на ячейку: видимость ячейки выводится из видимости её соседа СО
 * СТОРОНЫ ИСТОЧНИКА, уже посчитанной, потому что обход идёт в порядке удаления
 * от источника. Работа, которую луч делает заново каждый раз, здесь делится
 * между всеми ячейками — в этом и весь выигрыш, а не в мелкой оптимизации.
 *
 * ПОРЯДОК ОБХОДА И ЕСТЬ ОКТАНТНЫЙ: знак `(c − s)` по каждой оси задаёт октант,
 * и внутри октанта координата обходится ОТ источника, поэтому сосед со стороны
 * источника заведомо посчитан. Октантов восемь, и они покрывают сетку целиком.
 *
 * ЧЕГО ЭТОТ СВИП НЕ ДЕЛАЕТ — СКАЗАНО ЗДЕСЬ, ЧТОБЫ НЕ ЧИТАЛОСЬ СДЕЛАННЫМ. Он
 * несёт ВИДИМОСТЬ (затенение), а не радиантность: перенос энергии с отскоком —
 * следующая часть. Тень получается ПОЛУТЕНЬЮ только за счёт нескольких образцов
 * площадки, как и у луча.
 *
 * ПРИБЛИЖЕНИЕ, КОТОРОЕ НАДО НАЗВАТЬ: сосед берётся по ТРЁМ осям, и видимость
 * наследуется как максимум по ним. Это распространение вдоль ступенчатого пути,
 * а не вдоль прямой, поэтому у длинных косых теней край поедет. Насколько —
 * МЕРИТСЯ поячеечно против луча (А763), а не оценивается на глаз. */

/* Сетка свипа грубее сетки поля: тень не обязана иметь разрешение поверхности, а
 * цена свипа кубична по стороне. Уровень назван числом и проверяется замером
 * расхождения с эталоном-лучом. */
#ifndef HZ_SWEEP_DROP
#define HZ_SWEEP_DROP 2
#endif

/* РАЗРЕШЕНИЕ ПРОФИЛЯ НА ГРАНИ (Ф1', §502). `1` — прежнее поведение: одно число
 * на ячейку. Больше — профиль `P×P` подпроб, несущий ФОРМУ тени. Число НЕ
 * выводится из допуска (это следующий шаг), а перебирается замером: §425
 * требует сперва кривую «цена — точность» при фиксированных `P`. */
#ifndef HZ_SWEEP_P
#define HZ_SWEEP_P 1
#endif

/* НЕГАТИВНЫЙ КОНТРОЛЬ Ф1' (§502): профиль БЕЗ СДВИГА — все подпробы берутся из
 * одной клетки соседа. Расхождение обязано вернуться к уровню P = 1. */
static int g_noshift = 0;
/* Ф2' (§505): перенос вдоль ПРЯМОЙ вместо трёх осевых соседей. */
static int g_raysweep = 0;
/* §507: радиус засева окрестности источника маршем, в ячейках сетки свипа. */
static int32_t g_seed = 0;
/* §508: сколько проходов свипа подряд без обнуления поля. */
static int g_passes = 1;
/* Ф3' (§510): порог LOD для среза ИЗЛУЧАТЕЛЕЙ; 0 — излучатели те же, что приёмники. */
static double g_emitthr = 0.0;
/* Ф4' (§511): свип по дереву вместо плоской сетки. */
static int g_treesweep = 0;
/* Ф5. (§514): порог дробления дерева свипа и камера для него. */
static double g_sweepthr = 0.0, g_sweeppx = 1.0, g_sweepeye[3] = {0, 0, 0};
/* Ф6. (§516): угловой порог иерархического отскока; 0 — прежний гатер N². */
static double g_hgather = 0.0;

/* ---- НОСИТЕЛЬ КОСВЕННОГО СВЕТА — УЗЛЫ ДЕРЕВА (§597) ----------------------
 *
 * ЗАЧЕМ. Облучённость лежала в ячейках СРЕЗА, а срез камерозависим и строится
 * заново каждый кадр. Значит косвенный свет был привязан к камере ПО ПОСТРОЕНИЮ
 * и выбрасывался при каждом её движении — отсюда `5254` мс «вне кадра» (§593).
 * Узел дерева камеры не содержит, поэтому ядро считается ОДИН раз и живёт.
 *
 * ПОЧЕМУ МАССИВ ЗДЕСЬ, А НЕ В `hz_dctree`. `src/cut/` — ОБЩИЙ СЛОЙ ГЕОМЕТРИИ.
 * Светимость — не геометрия, и класть её в дерево разреза значило бы сращивать
 * слои, которые проект держит порознь сознательно. Индекс — номер узла. */
static float (*g_indnode)[3];
static int32_t g_indnode_n;
/* §762: флаг «узел заполнен раскладкой §597 или несёт заполненного потомка» —
 * различает легальный тёмный ноль от «данных нет» (класс А917). */
static uint8_t *g_indfill = NULL;

/* §744: агрегаты свипа для поячеечного сличения с ядром §597. Заполняются в
 * блоке развёртки (E_ind = полный − прямой, площадно-взвешенно по ячейке),
 * читаются после ind_core_build. Центры излучающих ячеек — для корзин по
 * расстоянию. */
/* §750: живая машинерия ядра §597 для двухсеточного цикла (keep-флаг велит
 * ind_core_build не освобождать срез, связи и ffv, а отдать сюда). */
static int g_dsa_keep = 0, g_dsa_ready = 0;
static hz_dcslice g_dsa_SF;
static linkcache g_dsa_LC;
static double *g_dsa_ffv = NULL;

/* §758: свиповое поле на УЗЛАХ дерева для кадра — num/den (E·площадь и
 * площадь), подъём суммами поддеревьев: кадровый срез с LOD читает любой
 * уровень как площадно-взвешенное среднее поддерева. */
static double *g_swEn = NULL, *g_swEd = NULL;
/* §798: второй узловой канал кадра — ПРЯМОЙ свет из ядрового сбора (гибрид
 * К-а §797); num/den той же механики, что g_swEn/g_swEd. */
static double *g_fcdn = NULL, *g_fcdd = NULL;
/* §798: режим прямого канала кадра для печатей: 0 — без fc, 1 — xfcelem
 * (прямой из элементов, §796), 2 — гибрид (прямой из ядрового сбора). */
static int g_fcmode = 0;
/* §798: подъём и спуск обобщены МАССИВАМИ — каналов стало два, а механика
 * (суммы поддеревьев, спуск к глубочайшему заполненному предку) одна. */
static void swE_lift(const hz_dctree *t, int32_t ni, double *En, double *Ed) {
  if (t->nd[ni].child0 < 0) return;
  for (int k = 0; k < 8; k++) {
    int32_t c = t->nd[ni].child0 + k;
    swE_lift(t, c, En, Ed);
    En[ni] += En[c];
    Ed[ni] += Ed[c];
  }
}

/* §760: спуск с запоминанием ГЛУБОЧАЙШЕГО предка, чьё поддерево несёт
 * свиповое поле; uplev — на сколько уровней выше запрошенного взято (0 —
 * дыры не было). Лечит и оборванные спуски. */
static int32_t node_swE_best(const hz_dctree *t, int lev, const hz_dccell *c, int *uplev,
                             const double *Ed) {
  int lvl = (int)c->lvl;
  if (lvl < 0 || lvl > lev) return -1;
  int32_t ni = 0, best = Ed[0] > 0.0 ? 0 : -1;
  int bestd = 0;
  for (int d = 0; d < lvl; d++) {
    if (t->nd[ni].child0 < 0) break;
    int bit = 0;
    for (int a = 0; a < 3; a++)
      if (((int32_t)c->lo[a] >> (lev - 1 - d)) & 1) bit |= 1 << a;
    ni = t->nd[ni].child0 + bit;
    if (Ed[ni] > 0.0) {
      best = ni;
      bestd = d + 1;
    }
  }
  *uplev = lvl - bestd;
  return best;
}

/* §772: ПРЕДИКАТ ОГРУБЛЕНИЯ ПРИЁМНИКОВ — лестница L = εR с fail-closed к
 * тонкому. Узел размера S листов берётся ячейкой, если он дальше D0·S от
 * глаза, его плоскости влезают в бюджет и ОБЪЁМНАЯ СВЕРКА сходится: флюид
 * грубого реза равен сумме листовых (невыпуклый материал в коробке даёт
 * макроскопическое расхождение и узел остаётся тонким). */
typedef struct {
  const hz_octree *t;
  const hz_cutmap *cm;
  const hz_facettab *ft;
  const hz_frame *fr;
  const unsigned char *occ;
  const uint8_t *smask;
  int occn, innerfluid;
  double eye[3], d0;
  int64_t nacc, nacc_empty, nrej_budget;
} c772;

/* §782: листовая сплошность для агрегации грубой ячейки — ТА ЖЕ логика, что
 * маска solid у заливки (А1253: рассинхрон = двойной учёт объёма). Класс
 * однороден по листу — хватает одной пробы (§772-прецедент). */
static int c772_leaf_solid(void *vc, const int32_t lo[3], int32_t size) {
  c772 *cx = vc;
  (void)size;
  size_t k = hz_occ_index(cx->occn, lo[0], lo[1], lo[2]);
  unsigned cls = cx->smask != NULL ? cx->smask[k] : 1u;
  return cx->smask != NULL && (cls == 0u || (cls == 2u && !cx->innerfluid)) &&
         !hz_occ_get(cx->occ, k);
}

/* §782: предикат огрубления похудел до ДАЛЬНОСТИ и БЮДЖЕТА КУСКОВ: объёмная
 * сверка §772 (пробный union-рез) не нужна — объём грубой ячейки собирается
 * ТОЧНОЙ агрегацией по листьям в cut3, предел §773 обойдён кусками Р-8.
 * Fail closed остаётся на бюджете: не влезли куски — спуск. */
static int c772_stop(void *vc, int32_t ni, const int32_t lo[3], int32_t size) {
  c772 *cx = vc;
  double dmin2 = 0.0;
  for (int a = 0; a < 3; a++) {
    double blo = cx->fr->o[a] + cx->fr->u[a] * (double)lo[a];
    double bhi = cx->fr->o[a] + cx->fr->u[a] * (double)(lo[a] + size);
    double d = cx->eye[a] < blo ? blo - cx->eye[a] : (cx->eye[a] > bhi ? cx->eye[a] - bhi : 0.0);
    dmin2 += d * d;
  }
  if (sqrt(dmin2) < cx->d0 * (double)size) return 0; /* близко — спуск, не отказ */
  int32_t refs[HZ_P3_MAXH];
  int nr = tr3_cut_subtree_refs(cx->t, cx->cm, ni, refs, HZ_P3_MAXH);
  if (nr < 0) {
    cx->nrej_budget++;
    return 0;
  }
  if (nr == 0) cx->nacc_empty++;
  cx->nacc++;
  return 1;
}

/* §794: провайдер треугольников сцены для hz_dc_facets_tris — перевод вершин
 * из мира в ЕДИНИЦЫ кадра; вырожденные отсеет сам мост (нулевая нормаль). */
typedef struct {
  const hz_objmesh *om;
  const frame *fr;
} otri794;

static int otri_get794(void *vc, int32_t i, double tv[3][3]) {
  otri794 *c = vc;
  const double *A, *B, *C;
  tri_verts(c->om, i, &A, &B, &C);
  for (int k = 0; k < 3; k++) {
    tv[0][k] = (A[k] - c->fr->org[k]) / c->fr->h;
    tv[1][k] = (B[k] - c->fr->org[k]) / c->fr->h;
    tv[2][k] = (C[k] - c->fr->org[k]) / c->fr->h;
  }
  return 1;
}

/* ---- §800: КЛАСТЕРИЗАЦИЯ КУСКОВ (А1285) -----------------------------------
 *
 * Жадное слияние смежных треугольников сцены в ВЫПУКЛЫЕ копланарные
 * куски-полигоны (≤ HZ_FACET_TVMAX вершин). Критерий слияния ОДИН и связывает
 * нормаль с офсетом (урок §5.1): все вершины кандидата в пределах δ от
 * плоскости СЕМЕНИ; плюс сонаправленность нормалей (двусторонние листы не
 * сливаются с изнанкой) и выпуклость в проекции на доминантную ось — знаком,
 * БЕЗ допуска (А1317). dmax куска — ФАКТИЧЕСКИЙ максимум |отклонения| (А1315).
 * Смежность — по ИНДЕКСАМ вершин OBJ; дубли вершин и T-стыки рвут её и режут
 * сжатие (А1316) — доля рёбер без пары печатается. */

/* δ по умолчанию, МЕТРЫ: на три порядка меньше листа lev=7 (0.9 м) и на
 * порядок больше float-шума вершин (~1e-5 м на сцене 115 м). Ключ xclustol=. */
#define HZ_CLUS_TOL 1e-4

typedef struct {
  uint64_t key; /* (min(u,v) << 32) | max(u,v) */
  int32_t tri;
  int32_t e;
} cedge800;

static int cmp_cedge800(const void *a, const void *b) {
  const cedge800 *x = a, *y = b;
  if (x->key != y->key) return x->key < y->key ? -1 : 1;
  if (x->tri != y->tri) return x->tri < y->tri ? -1 : 1;
  return (x->e > y->e) - (x->e < y->e);
}

typedef struct {
  int32_t *pv;  /* [np][HZ_FACET_TVMAX] индексы вершин сцены */
  int8_t *pnv;  /* [np] */
  float *pdmax; /* [np] метры */
  int32_t np;
} clus800;

/* выпуклость кандидата в проекции на доминантную ось нормали: все повороты
 * одного знака со знаком n[d] либо ровно ноль (коллинеарные вершины законны) */
static int clus800_convex(const hz_objmesh *om, const int32_t *vi, int nv, const double n[3]) {
  int d = 0;
  for (int a = 1; a < 3; a++)
    if (fabs(n[a]) > fabs(n[d])) d = a;
  int a1 = (d + 1) % 3, a2 = (d + 2) % 3;
  double sgn = n[d] > 0.0 ? 1.0 : -1.0;
  for (int j = 0; j < nv; j++) {
    const double *p0 = &om->v[3 * (size_t)vi[j]];
    const double *p1 = &om->v[3 * (size_t)vi[(j + 1) % nv]];
    const double *p2 = &om->v[3 * (size_t)vi[(j + 2) % nv]];
    double cr = (p1[a1] - p0[a1]) * (p2[a2] - p1[a2]) - (p1[a2] - p0[a2]) * (p2[a1] - p1[a1]);
    if (sgn * cr < 0.0) return 0;
  }
  return 1;
}

/* поворот в вершине j: 0.0 — битовая коллинеарность (кандидат на удаление) */
static double clus800_turn(const hz_objmesh *om, const int32_t *vi, int nv, const double n[3],
                           int j) {
  int d = 0;
  for (int a = 1; a < 3; a++)
    if (fabs(n[a]) > fabs(n[d])) d = a;
  int a1 = (d + 1) % 3, a2 = (d + 2) % 3;
  const double *p0 = &om->v[3 * (size_t)vi[(j + nv - 1) % nv]];
  const double *p1 = &om->v[3 * (size_t)vi[j]];
  const double *p2 = &om->v[3 * (size_t)vi[(j + 1) % nv]];
  return (p1[a1] - p0[a1]) * (p2[a2] - p1[a2]) - (p1[a2] - p0[a2]) * (p2[a1] - p1[a1]);
}

static void clus800_build(const hz_objmesh *om, double tol, clus800 *cl) {
  double t0 = now_s();
  int32_t nt = om->nt;
  cedge800 *ed = malloc(3 * (size_t)nt * sizeof *ed);
  int32_t *adj = malloc(3 * (size_t)nt * sizeof *adj);
  double *tn = malloc(3 * (size_t)nt * sizeof *tn);
  uint8_t *used = calloc((size_t)nt, 1);
  cl->pv = malloc((size_t)nt * HZ_FACET_TVMAX * sizeof *cl->pv);
  cl->pnv = malloc((size_t)nt * sizeof *cl->pnv);
  cl->pdmax = malloc((size_t)nt * sizeof *cl->pdmax);
  cl->np = 0;
  if (ed == NULL || adj == NULL || tn == NULL || used == NULL || cl->pv == NULL ||
      cl->pnv == NULL || cl->pdmax == NULL)
    exit(1);
  /* нормали (единичные) и рёбра; вырожденные помечаются использованными */
  for (int32_t t = 0; t < nt; t++) {
    const double *A = &om->v[3 * (size_t)om->f[3 * (size_t)t]];
    const double *B = &om->v[3 * (size_t)om->f[3 * (size_t)t + 1]];
    const double *C = &om->v[3 * (size_t)om->f[3 * (size_t)t + 2]];
    double e1[3], e2[3], nn[3];
    for (int k = 0; k < 3; k++) {
      e1[k] = B[k] - A[k];
      e2[k] = C[k] - A[k];
    }
    nn[0] = e1[1] * e2[2] - e1[2] * e2[1];
    nn[1] = e1[2] * e2[0] - e1[0] * e2[2];
    nn[2] = e1[0] * e2[1] - e1[1] * e2[0];
    double ml = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
    if (!(ml > 0.0)) used[t] = 1;
    for (int k = 0; k < 3; k++)
      tn[3 * (size_t)t + (size_t)k] = ml > 0.0 ? nn[k] / ml : 0.0;
    for (int e = 0; e < 3; e++) {
      int32_t u = om->f[3 * (size_t)t + (size_t)e];
      int32_t v = om->f[3 * (size_t)t + (size_t)((e + 1) % 3)];
      uint64_t lo9 = (uint64_t)(uint32_t)(u < v ? u : v), hi9 = (uint64_t)(uint32_t)(u < v ? v : u);
      ed[3 * (size_t)t + (size_t)e].key = (lo9 << 32) | hi9;
      ed[3 * (size_t)t + (size_t)e].tri = t;
      ed[3 * (size_t)t + (size_t)e].e = e;
      adj[3 * (size_t)t + (size_t)e] = -1;
    }
  }
  qsort(ed, 3 * (size_t)nt, sizeof *ed, cmp_cedge800);
  int64_t nopen = 0, nnonm = 0;
  {
    int64_t i = 0, n3 = 3 * (int64_t)nt;
    while (i < n3) {
      int64_t j = i;
      while (j < n3 && ed[j].key == ed[i].key)
        j++;
      if (j - i == 2) {
        adj[3 * (size_t)ed[i].tri + (size_t)ed[i].e] = ed[i + 1].tri;
        adj[3 * (size_t)ed[i + 1].tri + (size_t)ed[i + 1].e] = ed[i].tri;
      } else if (j - i == 1)
        nopen++;
      else
        nnonm += j - i;
      i = j;
    }
  }
  free(ed);
  /* жадный рост: полигон = список индексов вершин + владелец каждого
   * граничного ребра (tri*4 + e; −1 — ребро закрыто удалением коллинеарной) */
  int64_t nvhist[HZ_FACET_TVMAX + 1];
  memset(nvhist, 0, sizeof nvhist);
  double dmx_all = 0.0;
  int64_t nmerged = 0;
  for (int32_t t = 0; t < nt; t++) {
    if (used[t]) continue;
    used[t] = 1;
    int32_t vi[HZ_FACET_TVMAX];
    int64_t own[HZ_FACET_TVMAX];
    int nv = 3;
    for (int e = 0; e < 3; e++) {
      vi[e] = om->f[3 * (size_t)t + (size_t)e];
      own[e] = 4 * (int64_t)t + e;
    }
    const double *n0 = &tn[3 * (size_t)t];
    double off0 = n0[0] * om->v[3 * (size_t)vi[0]] + n0[1] * om->v[3 * (size_t)vi[0] + 1] +
                  n0[2] * om->v[3 * (size_t)vi[0] + 2];
    double dmx = 0.0;
    int progress = 1;
    while (progress && nv < HZ_FACET_TVMAX) {
      progress = 0;
      for (int j = 0; j < nv; j++) {
        if (own[j] < 0) continue;
        int32_t ot = (int32_t)(own[j] >> 2);
        int oe = (int)(own[j] & 3);
        int32_t nb = adj[3 * (size_t)ot + (size_t)oe];
        if (nb < 0 || used[nb]) continue;
        const double *nnb = &tn[3 * (size_t)nb];
        if (!(nnb[0] * n0[0] + nnb[1] * n0[1] + nnb[2] * n0[2] > 0.0)) continue;
        /* ребро полигона (u, v); у согласованно обмотанного соседа оно (v, u) */
        int32_t u = vi[j], v = vi[(j + 1) % nv];
        int enb = -1;
        for (int e = 0; e < 3; e++)
          if (om->f[3 * (size_t)nb + (size_t)e] == v &&
              om->f[3 * (size_t)nb + (size_t)((e + 1) % 3)] == u)
            enb = e;
        if (enb < 0) continue; /* несогласованная обмотка — не сливать */
        int32_t w = om->f[3 * (size_t)nb + (size_t)((enb + 2) % 3)];
        int dup = 0;
        for (int q = 0; q < nv; q++)
          if (vi[q] == w) dup = 1;
        if (dup) continue; /* не-простой полигон */
        double dw = fabs(n0[0] * om->v[3 * (size_t)w] + n0[1] * om->v[3 * (size_t)w + 1] +
                         n0[2] * om->v[3 * (size_t)w + 2] - off0);
        if (!(dw <= tol)) continue;
        int32_t cand[HZ_FACET_TVMAX + 1];
        for (int q = 0; q <= j; q++)
          cand[q] = vi[q];
        cand[j + 1] = w;
        for (int q = j + 1; q < nv; q++)
          cand[q + 1] = vi[q];
        if (!clus800_convex(om, cand, nv + 1, n0)) continue;
        /* принять: w после j; владельцы новых рёбер — два других ребра nb */
        for (int q = nv; q > j + 1; q--) {
          vi[q] = vi[q - 1];
          own[q] = own[q - 1];
        }
        vi[j + 1] = w;
        own[j] = 4 * (int64_t)nb + ((enb + 1) % 3);     /* (u, w) */
        own[j + 1] = 4 * (int64_t)nb + ((enb + 2) % 3); /* (w, v) */
        nv++;
        used[nb] = 1;
        nmerged++;
        if (dw > dmx) dmx = dw;
        /* битово-коллинеарные вершины выбрасываются (страйпы не упираются в
         * потолок вершин); ребро слитого отрезка закрывается (own = −1) —
         * смежность через него в сетке не существует (А1317) */
        for (int q = 0; q < nv && nv > 3; q++) {
          /* битовый ноль поворота — без ==: гейт запрещает равенство плавучих */
          if (fabs(clus800_turn(om, vi, nv, n0, q)) > 0.0) continue;
          int p = (q + nv - 1) % nv;
          for (int r = q; r < nv - 1; r++) {
            vi[r] = vi[r + 1];
            own[r] = own[r + 1];
          }
          nv--;
          own[p >= nv ? nv - 1 : p] = -1;
          q--;
        }
        progress = 1;
        break;
      }
    }
    int32_t *dst = &cl->pv[(size_t)cl->np * HZ_FACET_TVMAX];
    for (int q = 0; q < nv; q++)
      dst[q] = vi[q];
    cl->pnv[cl->np] = (int8_t)nv;
    cl->pdmax[cl->np] = (float)dmx;
    cl->np++;
    nvhist[nv]++;
    if (dmx > dmx_all) dmx_all = dmx;
  }
  free(adj);
  free(tn);
  free(used);
  printf("   §800 КЛАСТЕРЫ: %d треугольников -> %d кусков (×%.2f, поглощено %lld); nv: ", om->nt,
         cl->np, (double)om->nt / (double)(cl->np > 0 ? cl->np : 1), (long long)nmerged);
  for (int q = 3; q <= HZ_FACET_TVMAX; q++)
    printf("%d:%lld ", q, (long long)nvhist[q]);
  printf("; dmax макс %.3e м (допуск %.1e); рёбер без пары %lld, немногообразных %lld; %.2f с\n",
         dmx_all, tol, (long long)nopen, (long long)nnonm, now_s() - t0);
}

static void clus800_free(clus800 *cl) {
  free(cl->pv);
  free(cl->pnv);
  free(cl->pdmax);
  memset(cl, 0, sizeof *cl);
}

/* провайдер полигонов для hz_dc_facets_polys: мир -> единицы кадра */
typedef struct {
  const hz_objmesh *om;
  const frame *fr;
  const clus800 *cl;
} opoly800;

static int opoly_get800(void *vc, int32_t i, double pv[][3], double *dmax) {
  opoly800 *c = vc;
  int nv = c->cl->pnv[i];
  const int32_t *vi = &c->cl->pv[(size_t)i * HZ_FACET_TVMAX];
  for (int q = 0; q < nv; q++)
    for (int k = 0; k < 3; k++)
      pv[q][k] = (c->om->v[3 * (size_t)vi[q] + (size_t)k] - c->fr->org[k]) / c->fr->h;
  *dmax = (double)c->cl->pdmax[i];
  return nv;
}

/* §778: клип треугольника к коробке — Сазерленд–Ходжман по шести полуплоскостям.
 * ДИАГНОСТИКА (ключ xleak): эталон покрытия, независимый от рабочего пути cut3;
 * плавучка здесь законна — рабочая геометрия этим не пользуется. */
static int leak_clip_ax(const double (*in)[3], int nin, double (*out)[3], int ax, double val,
                        int keep_less) {
  int nout = 0;
  for (int i = 0; i < nin; i++) {
    const double *a = in[i], *b = in[(i + 1) % nin];
    double da = keep_less ? val - a[ax] : a[ax] - val;
    double db = keep_less ? val - b[ax] : b[ax] - val;
    if (da >= 0.0) memcpy(out[nout++], a, 3 * sizeof(double));
    if ((da >= 0.0) != (db >= 0.0)) {
      double t = da / (da - db);
      for (int k = 0; k < 3; k++)
        out[nout][k] = a[k] + t * (b[k] - a[k]);
      nout++;
    }
  }
  return nout;
}

/* Площадь куска треугольника в коробке [lo, hi]. Ёмкость 16: выпуклый клип
 * добавляет не больше одной вершины на полуплоскость, 3 + 6 = 9. */
static double leak_tri_box_area(const double tri[3][3], const double lo[3], const double hi[3]) {
  double A[16][3], B[16][3];
  memcpy(A, tri, 9 * sizeof(double));
  int n = 3;
  for (int ax = 0; ax < 3; ax++) {
    n = leak_clip_ax((const double (*)[3])A, n, B, ax, lo[ax], 0);
    if (n < 3) return 0.0;
    n = leak_clip_ax((const double (*)[3])B, n, A, ax, hi[ax], 1);
    if (n < 3) return 0.0;
  }
  double sx = 0.0, sy = 0.0, sz = 0.0;
  for (int i = 1; i + 1 < n; i++) {
    double u[3], w[3];
    for (int k = 0; k < 3; k++) {
      u[k] = A[i][k] - A[0][k];
      w[k] = A[i + 1][k] - A[0][k];
    }
    sx += u[1] * w[2] - u[2] * w[1];
    sy += u[2] * w[0] - u[0] * w[2];
    sz += u[0] * w[1] - u[1] * w[0];
  }
  return 0.5 * sqrt(sx * sx + sy * sy + sz * sz);
}

/* §762: подъём флага заполнения предкам и спуск с last-good по флагу. */
static int fill_lift(const hz_dctree *t, int32_t ni) {
  int any = g_indfill[ni];
  if (t->nd[ni].child0 >= 0)
    for (int k = 0; k < 8; k++)
      any |= fill_lift(t, t->nd[ni].child0 + k);
  g_indfill[ni] = (uint8_t)(any ? 1 : 0);
  return any;
}

static int32_t node_fill_best(const hz_dctree *t, int lev, const hz_dccell *c, int *uplev) {
  int lvl = (int)c->lvl;
  if (lvl < 0 || lvl > lev) return -1;
  int32_t ni = 0, best = g_indfill != NULL && g_indfill[0] ? 0 : -1;
  int bestd = 0;
  for (int d = 0; d < lvl; d++) {
    if (t->nd[ni].child0 < 0) break;
    int bit = 0;
    for (int a = 0; a < 3; a++)
      if (((int32_t)c->lo[a] >> (lev - 1 - d)) & 1) bit |= 1 << a;
    ni = t->nd[ni].child0 + bit;
    if (g_indfill != NULL && g_indfill[ni]) {
      best = ni;
      bestd = d + 1;
    }
  }
  *uplev = lvl - bestd;
  return best;
}

/* §774: потолок q̂ и минимум тактов для его оценки. Значения выведены в §774:
 * выше замеренного диапазона сжатия 0.61…0.89 замыкание умножает шум последнего
 * приращения как 1/(1−q); для оценки нужны ДВА приращения, и такт 1 — весь
 * прямой свет от нуля, вне геометрического режима. Файловая область видимости:
 * лестницей пользуется и стенд §808 (база кадра 2 — та же лестница). */
static const double HZ_TAILQ_MAX = 0.95;
static const int HZ_TAILN_MIN = 3;

/* ---- §808: СТЕНД ИНКРЕМЕНТНОСТИ СВЕТА (xwarm=, xscramble) ------------------
 *
 * Двухкадровый стенд (план §808, аудит А1346…А1349). Кадр 1 — холодная лестница
 * xbounce при камере cam; стенд сдвигает глаз на (DX,DY,DZ) и строит кадр 2
 * тёплым путём: пересборка грейда+реза+Ke (ОБЩАЯ часть кадра 2 — А1346: она
 * нужна и тёплому пути, и холодной базе, и в выигрыш инкрементности не
 * входит), ремап состояния, ОДНА-ДВЕ корректирующие прокидки К76 + замыкание
 * хвоста; рядом — же-процессная холодная база: та же лестница на сетке кадра
 * 2. Метрики времени печатаются раздельно.
 *
 * РЕМАП ТРЁХ НОСИТЕЛЕЙ РАЗДЕЛЁН (А1348): φ — вложенно через ОБЩЕЕ октодерево
 * (обе сетки суть наборы узлов одного ot, а два узла октодерева либо вложены,
 * либо disjoint — соответствие есть лес по построению); sout — рестрикцией по
 * ПЛОЩАДИ с контролем Σ sout·area (носитель — радианс, сохранять надо
 * энергию, не значение); bout НЕ ремапится — грань новой сетки не имеет
 * вложенного прообраза, он восстанавливается первой прокидкой из φ (в каноне
 * bout и вовсе не читается: wall_rho NULL, solid_rho 0). Поэтому варианта
 * «прокидок ноль» не существует — минимум одна, что план и предписывал.
 *
 * ДЕТЕКТОР СМЕНЫ ПОЛОСТИ (А1326): заливка идёт волной от ячейки КАМЕРЫ (Ш12),
 * и телепорт в другую полость меняет её целиком. Проверка дёшевa и точна:
 * класс ячейки НОВОГО глаза в маске кадра 1. «Полость» — та же компонента
 * связности (заливка отмечает ровно одну); «тело»/«наружное» — полость иная
 * или камера в материале: тёплый путь запрещён, кадр 2 считается холодным. */
typedef struct {
  /* сцена и дерево — камера-НЕзависимые, общие для обоих кадров */
  const hz_octree *ot;
  const hz_frame *ofr;
  const hz_cutmap *cm;
  const hz_facettab *ft;
  const hz_objmesh *om;
  void *ctctx; /* celltris* — ct_list принимает void* */
  const opyr *occ;
  const unsigned char *smask; /* заливка кадра 1 (Ш13) */
  int occn, innerfluid, lev;
  double xcoarse;
  int xfernosolid, xmatrho, xmatfar;
  double xrho, xrhoscale, xthin, xsemin, xtailq;
  /* кадр 1: сетка, рез, закрытое состояние */
  const tr3_mesh *m1;
  const tr3_cut *c1;
  const double *phi1;
  double *sout1; /* НЕ const: НК-скрэмбл портит хранилище (А1349) */
  double qhat1;  /* q̂ лестницы кадра 1; ≤0 — замыкания не было */
  double emitpow1;
  /* Р3 §810/А1357: альбедо фасетов кадра 1 для ПЕРЕИСПОЛЬЗОВАНИЯ (свойство
   * фасета, не сетки; уже умножено на xrhoscale в кадре 1 — повторное
   * умножение запрещено). NULL — прежний скоринг по элементам сетки 2. */
  const double *frho1;
  int xrhocal; /* калибровка: построить ОБЕ версии, сверить, напечатать */
  /* задача-шаблон: ординаты, лимитер, reltol — всё камера-независимое */
  const tr3_problem *prob;
  double eye1[3], shift[3];
  int nwarm, xbounce, scramble;
} warm808;

/* Узел ot на глубине depth, содержащий коробку с углом clo (листовая решётка
 * 2^log2size). −1 — дерево на этом пути мельче. Та же нумерация детей, что в
 * node_of_cell: бит оси a — разряд log2(size) координаты a. */
static int32_t warm808_node_at(const hz_octree *t, int depth, const int32_t clo[3]) {
  int32_t ni = 0, size = (int32_t)1 << t->log2size;
  for (int d = 0; d < depth; d++) {
    if (t->nodes[ni].child0 < 0) return -1;
    size >>= 1;
    int bit = 0;
    for (int a = 0; a < 3; a++)
      if (clo[a] & size) bit |= 1 << a;
    ni = t->nodes[ni].child0 + bit;
  }
  return ni;
}

/* ЯЧЕЙКИ СЕТКИ 1 под узлом ni: спуск останавливается на них — глубже сетки 1
 * нет. Переполнение считается (fail closed: усредняется собранное, печатается). */
static void warm808_collect(const hz_octree *t, const tr3_mesh *m1, int32_t ni, const int32_t lo[3],
                            int32_t size, int32_t *out, int32_t *n, int32_t cap, int32_t *nover) {
  if (m1->cellof[ni] >= 0) {
    if (*n < cap)
      out[(*n)++] = m1->cellof[ni];
    else
      (*nover)++;
    return;
  }
  if (t->nodes[ni].child0 < 0) return;
  int32_t half = size / 2;
  for (int k = 0; k < 8; k++) {
    int32_t clo[3] = {lo[0] + ((k & 1) ? half : 0), lo[1] + ((k & 2) ? half : 0),
                      lo[2] + ((k & 4) ? half : 0)};
    warm808_collect(t, m1, t->nodes[ni].child0 + k, clo, half, out, n, cap, nover);
  }
}

/* Ячейка СЕТКИ 1, содержащая мировую точку p: спуск от корня, останов на
 * первой ячейке сетки 1. −1 — вне куба либо спуск ушёл в лист, не встретив
 * ячейки сетки 1 (у ремапа sout это «источника нет», а не ошибка). */
static int32_t warm808_cell_at(const hz_octree *t, const tr3_mesh *m1, const hz_frame *ofr,
                               const double p[3]) {
  int32_t pi[3];
  for (int a = 0; a < 3; a++) {
    double f = floor((p[a] - ofr->o[a]) / ofr->u[a]);
    if (!(f >= 0.0) || f > (double)((int32_t)1 << t->log2size) - 1.0) return -1;
    pi[a] = (int32_t)f;
  }
  int32_t ni = 0, size = (int32_t)1 << t->log2size;
  while (m1->cellof[ni] < 0 && t->nodes[ni].child0 >= 0) {
    size >>= 1;
    int bit = 0;
    for (int a = 0; a < 3; a++)
      if (pi[a] & size) bit |= 1 << a;
    ni = t->nodes[ni].child0 + bit;
  }
  return m1->cellof[ni] >= 0 ? m1->cellof[ni] : -1;
}

/* Русская множественная для строк прибора («1 прокидка», «2 прокидки»). */
static const char *warm808_plu(int n, const char *one, const char *few, const char *many) {
  int n10 = n % 10, n100 = n % 100;
  if (n10 == 1 && n100 != 11) return one;
  if (n10 >= 2 && n10 <= 4 && (n100 < 12 || n100 > 14)) return few;
  return many;
}

/* Пара (значение, элемент) для отбора top-1 % по |sout| (А1349) — детерминизм
 * без ГСЧ: сортировка по значению, при равных по индексу. */
typedef struct {
  double v;
  int32_t e;
} warm808_pv;

static int warm808_pv_cmp(const void *a, const void *b) {
  const warm808_pv *x = a, *y = b;
  if (x->v > y->v) return -1; /* по убыванию яркости; равные значения — по индексу */
  if (x->v < y->v) return 1;
  return (x->e > y->e) - (x->e < y->e);
}

static void warm808_run(const warm808 *cx) {
  const double t808 = now_s();
  double eye2[3];
  int zero808 = 1;
  for (int a = 0; a < 3; a++) {
    eye2[a] = cx->eye1[a] + cx->shift[a];
    if (fabs(cx->shift[a]) > 0.0) zero808 = 0; /* xwarm=0,0,0 даёт ровно нули */
  }
  printf(
      "   §808 СДВИГ КАМЕРЫ: (%+.2f, %+.2f, %+.2f) м — глаз (%.2f %.2f %.2f) -> (%.2f %.2f %.2f); "
      "%s\n",
      cx->shift[0], cx->shift[1], cx->shift[2], cx->eye1[0], cx->eye1[1], cx->eye1[2], eye2[0],
      eye2[1], eye2[2],
      cx->nwarm == 1 ? "одна корректирующая прокидка"
                     : "предсказанная ветка: две прокидки (А1348)");
  /* НК-1 (А1349): порча ЗНАКА хранения кадра 1 — до ремапа, чтобы тёплый путь
   * унаследовал её, а холодная база (с нуля) — нет. Переворачиваются ВСЕ
   * ЧЕТЫРЕ коэффициента выбранных элементов: портится элемент, не число. */
  if (cx->scramble && cx->c1->nse > 0) {
    warm808_pv *pv = malloc((size_t)cx->c1->nse * sizeof *pv);
    if (pv == NULL) exit(1);
    int32_t npos = 0;
    double pabs = 0.0, paff = 0.0;
    for (int32_t e = 0; e < cx->c1->nse; e++) {
      if (!(cx->c1->se[e].area > 0.0)) continue;
      double v = fabs(cx->sout1[4 * (size_t)e]);
      pv[npos].v = v;
      pv[npos].e = e;
      npos++;
      pabs += v * cx->c1->se[e].area;
    }
    qsort(pv, (size_t)npos, sizeof *pv, warm808_pv_cmp);
    int32_t ktop = npos / 100;
    if (ktop < 1) ktop = 1;
    if (ktop > npos) ktop = npos;
    for (int32_t j = 0; j < ktop; j++) {
      int32_t e = pv[j].e;
      paff += pv[j].v * cx->c1->se[e].area;
      for (int q = 0; q < 4; q++)
        cx->sout1[4 * (size_t)e + (size_t)q] *= -1.0;
    }
    printf("   §808 НК-СКРЭМБЛ: перевернут знак у %d из %d элементов хранения (top-1 %% по "
           "|sout|), доля |мощности| %.2f %%\n",
           (int)ktop, (int)npos, 100.0 * paff / (pabs > 0.0 ? pabs : 1.0));
    free(pv);
  }
  /* ДЕТЕКТОР ПОЛОСТИ (А1326). Класс ячейки нового глаза в маске кадра 1. */
  {
    int32_t e0[3];
    for (int a = 0; a < 3; a++) {
      double f = floor((eye2[a] - cx->ofr->o[a]) / cx->ofr->u[a]);
      if (f < 0.0) f = 0.0;
      if (f > (double)(cx->occn - 1)) f = (double)(cx->occn - 1);
      e0[a] = (int32_t)f;
    }
    unsigned cls = cx->smask[hz_occ_index(cx->occn, e0[0], e0[1], e0[2])];
    if (cls != 1u) {
      printf("   §808 ЗАЛИВКА СМЕНИЛАСЬ: класс ячейки нового глаза в маске кадра 1 — %s; "
             "тёплый путь ЗАПРЕЩЁН, кадр 2 — холодный пересчёт (стендом не повторяется: "
             "холодная лестница уже измерена как база)\n",
             cls == 0u ? "ТЕЛО" : "НАРУЖНОЕ");
      printf("   §808 ИНКРЕМЕНТНОСТЬ: холод —, тёплый —; заливка СМЕНИЛАСЬ — тёплый путь "
             "запрещён\n");
      return;
    }
    printf("   §808 ЗАЛИВКА: ячейка нового глаза — ПОЛОСТЬ в маске кадра 1 (та же компонента "
           "связности)\n");
  }
  /* ---- ОБЩАЯ ЧАСТЬ КАДРА 2: грейд + рез + раздача Ke (А1346) ---- */
  double tgc = now_s();
  c772 cb2;
  memset(&cb2, 0, sizeof cb2);
  cb2.t = cx->ot;
  cb2.cm = cx->cm;
  cb2.ft = cx->ft;
  cb2.fr = cx->ofr;
  cb2.occ = cx->occ->b[cx->lev];
  cb2.smask = cx->smask;
  cb2.occn = cx->occn;
  cb2.innerfluid = cx->innerfluid;
  for (int a = 0; a < 3; a++)
    cb2.eye[a] = eye2[a];
  cb2.d0 = cx->xcoarse;
  tr3_mesh mesh2;
  if (tr3_mesh_build_lod(&mesh2, cx->ot, cx->ofr, c772_stop, &cb2) != 0) exit(1);
  uint8_t *solid2 = calloc((size_t)mesh2.ncell, 1);
  if (solid2 == NULL) exit(1);
  for (int32_t ci = 0; ci < mesh2.ncell; ci++) {
    /* та же проба листьев, что у кадра 1 (Ш16/§782): сплошная = все листы сплошные */
    int32_t s9 = mesh2.csize[ci];
    int allsolid = 1;
    for (int32_t iz = 0; iz < s9 && allsolid; iz++)
      for (int32_t iy = 0; iy < s9 && allsolid; iy++)
        for (int32_t ix = 0; ix < s9 && allsolid; ix++) {
          size_t k = hz_occ_index(cx->occn, mesh2.clo[ci][0] + ix, mesh2.clo[ci][1] + iy,
                                  mesh2.clo[ci][2] + iz);
          if (!((cx->smask[k] == 0u || (cx->smask[k] == 2u && !cx->innerfluid)) &&
                !hz_occ_get(cx->occ->b[cx->lev], k)))
            allsolid = 0;
        }
    if (allsolid) solid2[ci] = 1u;
  }
  tr3_cut cut2;
  int crc2 = tr3_cut_build2(&cut2, &mesh2, cx->ft, cx->cm, cx->xfernosolid ? NULL : solid2,
                            cx->xcoarse > 0.0 ? c772_leaf_solid : NULL, &cb2);
  if (crc2 != 0) {
    fprintf(stderr, "§808: рез сетки кадра 2 не собрался (код %d)\n", crc2);
    exit(1);
  }
  /* §717/§711 — те же постфильтры, что у кадра 1, в том же порядке */
  if (cx->xsemin > 0.0)
    for (int32_t k = 0; k < cut2.nse; k++) {
      double h2 = cx->ofr->u[0] * cx->ofr->u[1];
      if (cut2.se[k].area > 0.0 && cut2.se[k].area < cx->xsemin * h2) {
        cut2.se[k].area = 0.0;
        memset(cut2.se[k].m, 0, sizeof cut2.se[k].m);
      }
    }
  if (cx->xthin > 0.0)
    for (int32_t ci = 0; ci < mesh2.ncell; ci++) {
      double vfl = cut2.mvol[ci][0][0];
      if (!(vfl > 0.0) || cut2.solid[ci]) continue;
      double s3 = (double)mesh2.csize[ci];
      double V = s3 * s3 * s3 * cx->ofr->u[0] * cx->ofr->u[1] * cx->ofr->u[2];
      if (V > 0.0 && vfl < cx->xthin * V) {
        cut2.solid[ci] = 1;
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
            cut2.mvol[ci][i][j] = 0.0;
      }
    }
  /* АЛЬБЕДО ФАСЕТОВ сетки 2. Р3 §810/А1357: по умолчанию ПЕРЕИСПОЛЬЗУЕТСЯ
   * frho кадра 1 — альбедо есть свойство ФАСЕТА (общая таблица ft), а
   * скоринг — лишь атрибуция по данной сетке (разброс А1135); frho1 уже несёт
   * множитель xrhoscale, повторное умножение запрещено. Калибровка xrhocal
   * строит скоринг-версию и печатает расхождение (порог 5 % — информация,
   * не переключатель: переключать поведение в прогоне значило бы мерить
   * другую физику, чем напечатана). */
  const int32_t nse2 = cut2.nse;
  double *frho2 = calloc((size_t)(cx->ft->n > 0 ? cx->ft->n : 1), sizeof *frho2);
  double *eemit2 = calloc((size_t)(nse2 > 0 ? nse2 : 1), sizeof *eemit2);
  double *sigt2 = calloc((size_t)mesh2.ncell, sizeof *sigt2);
  double *sigs2 = calloc((size_t)mesh2.ncell, sizeof *sigs2);
  if (frho2 == NULL || eemit2 == NULL || sigt2 == NULL || sigs2 == NULL) exit(1);
  int frho_reused = 0;
  if (cx->frho1 != NULL) {
    memcpy(frho2, cx->frho1, (size_t)cx->ft->n * sizeof *frho2);
    frho_reused = 1;
  } else {
    for (int32_t i = 0; i < cx->ft->n; i++)
      frho2[i] = cx->xrho;
  }
  if (cx->xmatrho) {
    double *frnum = calloc((size_t)cx->ft->n, sizeof *frnum);
    double *frden = calloc((size_t)cx->ft->n, sizeof *frden);
    if (frnum == NULL || frden == NULL) exit(1);
    enum { LS808 = 4096 };
    int32_t ls808[LS808];
    int64_t nlsclip808 = 0;
    for (int32_t k = 0; k < nse2; k++) {
      if (cut2.se[k].nv <= 0) continue;
      int32_t ci = cut2.se[k].cell;
      if (ci < 0 || ci >= mesh2.ncell) continue;
      const int32_t *ls = NULL;
      int32_t nls = 0;
      if (mesh2.csize[ci] > 1 && cut2.se[k].facet >= 0 && cut2.se[k].facet < cx->ft->n &&
          cx->ft->f[cut2.se[k].facet].bounded) {
        const hz_facet *fp8 = &cx->ft->f[cut2.se[k].facet];
        int32_t blo8[3], bhi8[3];
        int ntv8 = fp8->tnv >= 3 && fp8->tnv <= HZ_FACET_TVMAX ? fp8->tnv : 3;
        for (int a = 0; a < 3; a++) {
          double mn = fp8->tv[0][a], mx = fp8->tv[0][a];
          for (int q2 = 1; q2 < ntv8; q2++) {
            if (fp8->tv[q2][a] < mn) mn = fp8->tv[q2][a];
            if (fp8->tv[q2][a] > mx) mx = fp8->tv[q2][a];
          }
          blo8[a] = (int32_t)floor(mn);
          bhi8[a] = (int32_t)floor(mx);
          if (blo8[a] < 0) blo8[a] = 0;
          if (bhi8[a] > cx->occn - 1) bhi8[a] = cx->occn - 1;
        }
        int32_t nu = 0;
        for (int32_t iz = blo8[2]; iz <= bhi8[2]; iz++)
          for (int32_t iy = blo8[1]; iy <= bhi8[1]; iy++)
            for (int32_t ix = blo8[0]; ix <= bhi8[0]; ix++) {
              int32_t cl8[3] = {ix, iy, iz};
              const int32_t *l2 = NULL;
              int32_t n2 = ct_list(cx->ctctx, cl8, &l2);
              for (int32_t q2 = 0; q2 < n2; q2++) {
                if (nu >= LS808) {
                  nlsclip808++;
                  break;
                }
                ls808[nu++] = l2[q2];
              }
            }
        ls = ls808;
        nls = nu;
      } else {
        int32_t cellc[3] = {mesh2.clo[ci][0], mesh2.clo[ci][1], mesh2.clo[ci][2]};
        nls = ct_list(cx->ctctx, cellc, &ls);
      }
      if (nls == 0) continue;
      double ec[3] = {0, 0, 0};
      for (int q2 = 0; q2 < cut2.se[k].nv; q2++)
        for (int a = 0; a < 3; a++)
          ec[a] += cut2.se[k].v[q2][a] / (double)cut2.se[k].nv;
      const double *en = cut2.se[k].n;
      int32_t tbest = ls[0];
      double sbest = cx->xmatfar ? -1.0 : 1e300, cbest = 1e300;
      for (int32_t q2 = 0; q2 < nls; q2++) {
        const double *A3, *B3, *C3;
        tri_verts(cx->om, ls[q2], &A3, &B3, &C3);
        double e1[3] = {B3[0] - A3[0], B3[1] - A3[1], B3[2] - A3[2]};
        double e2[3] = {C3[0] - A3[0], C3[1] - A3[1], C3[2] - A3[2]};
        double nt[3];
        nt[0] = e1[1] * e2[2] - e1[2] * e2[1];
        nt[1] = e1[2] * e2[0] - e1[0] * e2[2];
        nt[2] = e1[0] * e2[1] - e1[1] * e2[0];
        double nl = sqrt(nt[0] * nt[0] + nt[1] * nt[1] + nt[2] * nt[2]);
        if (!(nl > 0.0)) continue;
        double dist = 0.0, dot = 0.0, cd = 0.0;
        for (int a = 0; a < 3; a++) {
          dist += nt[a] / nl * (ec[a] - A3[a]);
          dot += nt[a] / nl * en[a];
          double tc = (A3[a] + B3[a] + C3[a]) / 3.0 - ec[a];
          cd += tc * tc;
        }
        double sc = fabs(dist) + cx->ofr->u[0] * (1.0 - fabs(dot));
        if (cx->xmatfar) {
          if (sc > sbest) {
            sbest = sc;
            tbest = ls[q2];
          }
        } else if (sc < sbest - 1e-6 * cx->ofr->u[0] ||
                   (sc < sbest + 1e-6 * cx->ofr->u[0] && cd < cbest)) {
          sbest = sc;
          cbest = cd;
          tbest = ls[q2];
        }
      }
      int32_t mi2 = cx->om->fm != NULL ? cx->om->fm[tbest] : 0;
      if (mi2 < 0 || mi2 >= cx->om->nmtl) mi2 = 0;
      double kd = cx->om->mtl[mi2].kd;
      if (!(kd >= 0.0 && kd <= 1.0)) kd = 0.5;
      int32_t fi2 = cut2.se[k].facet;
      double ar = cut2.se[k].area;
      frnum[fi2] += kd * ar;
      frden[fi2] += ar;
    }
    for (int32_t i = 0; i < cx->ft->n; i++)
      frho2[i] = frden[i] > 0.0 ? frnum[i] / frden[i] : 0.5;
    free(frnum);
    free(frden);
    if (nlsclip808 > 0)
      printf("   §808 СКОРИНГ АЛЬБЕДО: усечено списков клеток %lld (бюджет %d треугольников, "
             "класс А1217)\n",
             (long long)nlsclip808, (int)LS808);
  }
  if (!frho_reused)
    for (int32_t i = 0; i < cx->ft->n; i++)
      frho2[i] *= cx->xrhoscale;
  /* Калибровка Р3 (xrhocal): построить скоринг-версию тем же кодом и сверить
   * с переиспользованной — расхождение печатается ДОРОГОЙ ОДНОКРАТНОЙ
   * проверкой, в рабочий путь не входит. */
  if (cx->xrhocal && frho_reused) {
    double *frs = calloc((size_t)(cx->ft->n > 0 ? cx->ft->n : 1), sizeof *frs);
    double *frnum = calloc((size_t)cx->ft->n, sizeof *frnum);
    double *frden = calloc((size_t)cx->ft->n, sizeof *frden);
    if (frs == NULL || frnum == NULL || frden == NULL) exit(1);
    for (int32_t i = 0; i < cx->ft->n; i++)
      frs[i] = cx->xrho;
    if (cx->xmatrho) {
      enum { LSC8 = 4096 };
      int32_t lsc8[LSC8];
      int64_t nclipc8 = 0;
      for (int32_t k = 0; k < nse2; k++) {
        if (cut2.se[k].nv <= 0) continue;
        int32_t ci = cut2.se[k].cell;
        if (ci < 0 || ci >= mesh2.ncell) continue;
        const int32_t *ls = NULL;
        int32_t nls = 0;
        if (mesh2.csize[ci] > 1 && cut2.se[k].facet >= 0 && cut2.se[k].facet < cx->ft->n &&
            cx->ft->f[cut2.se[k].facet].bounded) {
          const hz_facet *fp8 = &cx->ft->f[cut2.se[k].facet];
          int32_t blo8[3], bhi8[3];
          int ntv8 = fp8->tnv >= 3 && fp8->tnv <= HZ_FACET_TVMAX ? fp8->tnv : 3;
          for (int a = 0; a < 3; a++) {
            double mn = fp8->tv[0][a], mx = fp8->tv[0][a];
            for (int q2 = 1; q2 < ntv8; q2++) {
              if (fp8->tv[q2][a] < mn) mn = fp8->tv[q2][a];
              if (fp8->tv[q2][a] > mx) mx = fp8->tv[q2][a];
            }
            blo8[a] = (int32_t)floor(mn);
            bhi8[a] = (int32_t)floor(mx);
            if (blo8[a] < 0) blo8[a] = 0;
            if (bhi8[a] > cx->occn - 1) bhi8[a] = cx->occn - 1;
          }
          int32_t nu = 0;
          for (int32_t iz = blo8[2]; iz <= bhi8[2]; iz++)
            for (int32_t iy = blo8[1]; iy <= bhi8[1]; iy++)
              for (int32_t ix = blo8[0]; ix <= bhi8[0]; ix++) {
                int32_t cl8[3] = {ix, iy, iz};
                const int32_t *l2 = NULL;
                int32_t n2 = ct_list(cx->ctctx, cl8, &l2);
                for (int32_t q2 = 0; q2 < n2; q2++) {
                  if (nu >= LSC8) {
                    nclipc8++;
                    break;
                  }
                  lsc8[nu++] = l2[q2];
                }
              }
          ls = lsc8;
          nls = nu;
        } else {
          int32_t cellc[3] = {mesh2.clo[ci][0], mesh2.clo[ci][1], mesh2.clo[ci][2]};
          nls = ct_list(cx->ctctx, cellc, &ls);
        }
        if (nls == 0) continue;
        double ec[3] = {0, 0, 0};
        for (int q2 = 0; q2 < cut2.se[k].nv; q2++)
          for (int a = 0; a < 3; a++)
            ec[a] += cut2.se[k].v[q2][a] / (double)cut2.se[k].nv;
        const double *en = cut2.se[k].n;
        int32_t tbest = ls[0];
        double sbest = cx->xmatfar ? -1.0 : 1e300, cbest = 1e300;
        for (int32_t q2 = 0; q2 < nls; q2++) {
          const double *A3, *B3, *C3;
          tri_verts(cx->om, ls[q2], &A3, &B3, &C3);
          double e1[3] = {B3[0] - A3[0], B3[1] - A3[1], B3[2] - A3[2]};
          double e2[3] = {C3[0] - A3[0], C3[1] - A3[1], C3[2] - A3[2]};
          double nt[3];
          nt[0] = e1[1] * e2[2] - e1[2] * e2[1];
          nt[1] = e1[2] * e2[0] - e1[0] * e2[2];
          nt[2] = e1[0] * e2[1] - e1[1] * e2[0];
          double nl = sqrt(nt[0] * nt[0] + nt[1] * nt[1] + nt[2] * nt[2]);
          if (!(nl > 0.0)) continue;
          double dist = 0.0, dot = 0.0, cd = 0.0;
          for (int a = 0; a < 3; a++) {
            dist += nt[a] / nl * (ec[a] - A3[a]);
            dot += nt[a] / nl * en[a];
            double tc = (A3[a] + B3[a] + C3[a]) / 3.0 - ec[a];
            cd += tc * tc;
          }
          double sc = fabs(dist) + cx->ofr->u[0] * (1.0 - fabs(dot));
          if (cx->xmatfar) {
            if (sc > sbest) {
              sbest = sc;
              tbest = ls[q2];
            }
          } else if (sc < sbest - 1e-6 * cx->ofr->u[0] ||
                     (sc < sbest + 1e-6 * cx->ofr->u[0] && cd < cbest)) {
            sbest = sc;
            cbest = cd;
            tbest = ls[q2];
          }
        }
        int32_t mi2 = cx->om->fm != NULL ? cx->om->fm[tbest] : 0;
        if (mi2 < 0 || mi2 >= cx->om->nmtl) mi2 = 0;
        double kd = cx->om->mtl[mi2].kd;
        if (!(kd >= 0.0 && kd <= 1.0)) kd = 0.5;
        frnum[cut2.se[k].facet] += kd * cut2.se[k].area;
        frden[cut2.se[k].facet] += cut2.se[k].area;
      }
      for (int32_t i = 0; i < cx->ft->n; i++)
        frs[i] = frden[i] > 0.0 ? cx->xrhoscale * frnum[i] / frden[i] : cx->xrhoscale * 0.5;
      (void)nclipc8;
    } else
      for (int32_t i = 0; i < cx->ft->n; i++)
        frs[i] = cx->xrho * cx->xrhoscale;
    int64_t ndiff8 = 0;
    double dmax8 = 0.0;
    for (int32_t i = 0; i < cx->ft->n; i++) {
      double d8 = fabs(frs[i] - frho2[i]);
      if (d8 > 0.01) ndiff8++;
      if (d8 > dmax8) dmax8 = d8;
    }
    printf("   §808/Р3 КАЛИБРОВКА АЛЬБЕДО: фасетов с |frho1−frho2| > 0.01 — %lld из %d "
           "(%.2f %%), макс |Δ| %.4f%s\n",
           (long long)ndiff8, cx->ft->n,
           100.0 * (double)ndiff8 / (double)(cx->ft->n > 0 ? cx->ft->n : 1), dmax8,
           100.0 * (double)ndiff8 / (double)(cx->ft->n > 0 ? cx->ft->n : 1) > 5.0
               ? " — ПРЕВЫШЕН ПОРОГ 5 %%, НАХОДКА"
               : "");
    free(frs);
    free(frnum);
    free(frden);
  }
  /* §786-раздача Ke на сетке 2 — та же машина: энергия ячейки точным клипом
   * треугольников к листовым клеткам, раздача элементам по площадям. */
  double emitpow2b = 0.0;
  {
    int32_t *cellat6 =
        malloc((size_t)cx->occn * (size_t)cx->occn * (size_t)cx->occn * sizeof *cellat6);
    double *Ecell6 = calloc((size_t)mesh2.ncell, sizeof *Ecell6);
    double *Sarea6 = calloc((size_t)mesh2.ncell, sizeof *Sarea6);
    if (cellat6 == NULL || Ecell6 == NULL || Sarea6 == NULL) exit(1);
    for (size_t g = 0; g < (size_t)cx->occn * (size_t)cx->occn * (size_t)cx->occn; g++)
      cellat6[g] = -1;
    for (int32_t ci = 0; ci < mesh2.ncell; ci++)
      for (int32_t iz = 0; iz < mesh2.csize[ci]; iz++)
        for (int32_t iy = 0; iy < mesh2.csize[ci]; iy++)
          for (int32_t ix = 0; ix < mesh2.csize[ci]; ix++)
            cellat6[hz_occ_index(cx->occn, mesh2.clo[ci][0] + ix, mesh2.clo[ci][1] + iy,
                                 mesh2.clo[ci][2] + iz)] = ci;
    for (int32_t k = 0; k < nse2; k++)
      if (cut2.se[k].area > 0.0) Sarea6[cut2.se[k].cell] += cut2.se[k].area;
    for (int32_t t9 = 0; t9 < cx->om->nt; t9++) {
      int32_t mt = cx->om->fm != NULL ? cx->om->fm[t9] : 0;
      if (mt < 0 || mt >= cx->om->nmtl) mt = 0;
      const double *ke3 = cx->om->mtl[mt].ke3;
      double ke = (ke3[0] + ke3[1] + ke3[2]) / 3.0;
      if (!(ke > 0.0)) continue;
      /* нулевая инициализация — гейт: анализатор теряет заполнение через
       * вложенный цикл по q2/a (ложный класс, что у cellc/e1 выше) */
      double tri[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
      for (int q2 = 0; q2 < 3; q2++)
        for (int a = 0; a < 3; a++)
          tri[q2][a] = cx->om->v[3 * (size_t)cx->om->f[3 * (size_t)t9 + (size_t)q2] + (size_t)a];
      int32_t blo9[3] = {0, 0, 0}, bhi9[3] = {0, 0, 0};
      for (int a = 0; a < 3; a++) {
        double mn = tri[0][a], mx = tri[0][a];
        for (int q2 = 1; q2 < 3; q2++) {
          if (tri[q2][a] < mn) mn = tri[q2][a];
          if (tri[q2][a] > mx) mx = tri[q2][a];
        }
        blo9[a] = (int32_t)floor((mn - cx->ofr->o[a]) / cx->ofr->u[a]);
        bhi9[a] = (int32_t)floor((mx - cx->ofr->o[a]) / cx->ofr->u[a]);
        if (blo9[a] < 0) blo9[a] = 0;
        if (bhi9[a] > cx->occn - 1) bhi9[a] = cx->occn - 1;
      }
      for (int32_t iz = blo9[2]; iz <= bhi9[2]; iz++)
        for (int32_t iy = blo9[1]; iy <= bhi9[1]; iy++)
          for (int32_t ix = blo9[0]; ix <= bhi9[0]; ix++) {
            double blo[3] = {cx->ofr->o[0] + cx->ofr->u[0] * (double)ix,
                             cx->ofr->o[1] + cx->ofr->u[1] * (double)iy,
                             cx->ofr->o[2] + cx->ofr->u[2] * (double)iz};
            double bhi[3] = {blo[0] + cx->ofr->u[0], blo[1] + cx->ofr->u[1],
                             blo[2] + cx->ofr->u[2]};
            double A = leak_tri_box_area((const double (*)[3])tri, blo, bhi);
            if (!(A > 0.0)) continue;
            int32_t ci = cellat6[hz_occ_index(cx->occn, ix, iy, iz)];
            if (ci >= 0 && Sarea6[ci] > 0.0) Ecell6[ci] += ke * A;
          }
    }
    for (int32_t k = 0; k < nse2; k++) {
      int32_t ci = cut2.se[k].cell;
      if (!(cut2.se[k].area > 0.0) || !(Ecell6[ci] > 0.0)) continue;
      eemit2[k] = Ecell6[ci] / (3.14159265358979323846 * Sarea6[ci]);
      emitpow2b += Ecell6[ci] * cut2.se[k].area / Sarea6[ci];
    }
    free(cellat6);
    free(Ecell6);
    free(Sarea6);
  }
  double t_gc = now_s() - tgc;
  printf("   §808 ГРЕЙД2+РЕЗ2 (общая часть кадра 2, обеим веткам): ячеек %d (кадр 1: %d), "
         "элементов %d (было %d), роздано Ke %.4f (кадр 1: %.4f, Δ %+.2f %%); %.2f с\n",
         mesh2.ncell, cx->m1->ncell, nse2, cx->c1->nse, emitpow2b, cx->emitpow1,
         100.0 * (emitpow2b - cx->emitpow1) / (cx->emitpow1 > 0.0 ? cx->emitpow1 : 1.0), t_gc);
  /* ---- РЕМАП (А1348) ---- */
  const double trm0 = now_s();
  const size_t nph2 = (size_t)mesh2.ncell * 4, nso2 = (size_t)(nse2 > 0 ? nse2 : 1) * 4,
               nbo2 = (size_t)mesh2.nf * 4;
  double *phi2 = calloc(nph2, sizeof *phi2);
  double *sout2 = calloc(nso2, sizeof *sout2);
  int32_t *buf8 = malloc((size_t)4096 * sizeof *buf8);
  if (phi2 == NULL || sout2 == NULL || buf8 == NULL) exit(1);
  /* масса φ кадра 1 — до ремапа (контроль УБИВАЕТ-в) */
  double mass1 = 0.0, massabs1 = 0.0;
  for (int32_t c = 0; c < cx->m1->ncell; c++) {
    double w = cx->c1->mvol[c][0][0];
    mass1 += cx->phi1[4 * (size_t)c] * w;
    massabs1 += fabs(cx->phi1[4 * (size_t)c]) * w;
  }
  int64_t nid8 = 0, npr8 = 0, nre8 = 0, nov8 = 0;
  /* гистограмма глубины: d = уровень_новой − уровень_источника (0 — тождество,
   * + — продление вниз, − — рестрикция вверх); ±6 и дальше — крайние корзины.
   * ГЛУБИНА УЗЛА = lev − log2(csize): считаем через log2, а не спуском —
   * иначе перепутываются размер и глубина (поймано прогоном §808: «подъём»
   * проверял один корень). */
  int64_t dhist[13] = {0};
  for (int32_t ci = 0; ci < mesh2.ncell; ci++) {
    int32_t csize2 = mesh2.csize[ci];
    int k2 = 0;
    for (int32_t s = csize2; s > 1; s >>= 1)
      k2++;
    int depth2 = cx->lev - k2; /* глубина узла ni2 от корня */
    int32_t ni2 = mesh2.node[ci];
    int32_t found = -1, fdepth = -1;
    for (int d = depth2; d >= 0; d--) {
      int32_t ni = (d == depth2) ? ni2 : warm808_node_at(cx->ot, d, mesh2.clo[ci]);
      if (ni >= 0 && cx->m1->cellof[ni] >= 0) {
        found = cx->m1->cellof[ni];
        fdepth = d;
        break;
      }
    }
    int drec = 0;
    if (found >= 0 && fdepth == depth2) {
      /* тождество: тот же узел ot — ячейка есть у обеих сеток */
      memcpy(&phi2[4 * (size_t)ci], &cx->phi1[4 * (size_t)found], 4 * sizeof(double));
      nid8++;
    } else if (found >= 0) {
      /* ПРОДЛЕНИЕ: старая ячейка крупнее; поле линейно в её базисе
       * φ = m + Σ s_k·(x_k−c_k)/h_k, среднее в центре новой, наклоны ×h2/h1 */
      int32_t cj = found;
      const double *s1 = &cx->phi1[4 * (size_t)cj];
      double h1 = (double)cx->m1->csize[cj], h2d = (double)csize2;
      double m = s1[0];
      for (int k = 1; k < 4; k++) {
        double dc = ((double)mesh2.clo[ci][k - 1] + 0.5 * h2d) -
                    ((double)cx->m1->clo[cj][k - 1] + 0.5 * h1);
        m += s1[k] * dc / h1;
      }
      phi2[4 * (size_t)ci] = m;
      for (int k = 1; k < 4; k++)
        phi2[4 * (size_t)ci + (size_t)k] = s1[k] * (h2d / h1);
      npr8++;
      drec = depth2 - fdepth;
    } else {
      /* РЕСТРИКЦИЯ: сетка 1 мельче внутри; флюидо-взвешенное среднее её ячеек,
       * наклоны пересчитываются в базис крупной (×h_старой/h_новой) */
      int32_t ncol = 0, novc = 0;
      warm808_collect(cx->ot, cx->m1, ni2, mesh2.clo[ci], csize2, buf8, &ncol, 4096, &novc);
      nov8 += novc;
      double wsum = 0.0, msum = 0.0, ssum[3] = {0, 0, 0};
      double wdom = -1.0;
      int ddom = 0;
      for (int32_t j = 0; j < ncol; j++) {
        int32_t cj = buf8[j];
        double w = cx->c1->mvol[cj][0][0];
        if (!(w > 0.0)) continue;
        double ratio = (double)cx->m1->csize[cj] / (double)csize2;
        wsum += w;
        msum += w * cx->phi1[4 * (size_t)cj];
        for (int k = 1; k < 4; k++)
          ssum[k - 1] += w * cx->phi1[4 * (size_t)cj + (size_t)k] * ratio;
        int dj = 0;
        for (int32_t s = cx->m1->csize[cj]; s > 1; s >>= 1)
          dj++;
        if (w > wdom) {
          wdom = w;
          ddom = cx->lev - dj; /* глубина источника-доминанты */
        }
      }
      if (wsum > 0.0) {
        phi2[4 * (size_t)ci] = msum / wsum;
        for (int k = 1; k < 4; k++)
          phi2[4 * (size_t)ci + (size_t)k] = ssum[k - 1] / wsum;
      }
      nre8++;
      drec = depth2 - ddom;
    }
    if (drec < -6) drec = -6;
    if (drec > 6) drec = 6;
    dhist[drec + 6]++;
  }
  double mass2 = 0.0, massabs2 = 0.0;
  for (int32_t ci = 0; ci < mesh2.ncell; ci++) {
    double w = cut2.mvol[ci][0][0];
    mass2 += phi2[4 * (size_t)ci] * w;
    massabs2 += fabs(phi2[4 * (size_t)ci]) * w;
  }
  /* sout — рестрикция по площади: веерная триангуляция многоугольника,
   * центроиды треугольников с весами площадей ищут старый элемент ТОГО ЖЕ
   * фасета (номер фасета — индекс в общей таблице ft, сетки его не меняют).
   * Контроль — Σ sout·area до/после (энергия, УБИКАЕТ-в). */
  double spow1 = 0.0, spowabs1 = 0.0;
  for (int32_t e = 0; e < cx->c1->nse; e++) {
    if (!(cx->c1->se[e].area > 0.0)) continue;
    spow1 += cx->sout1[4 * (size_t)e] * cx->c1->se[e].area;
    spowabs1 += fabs(cx->sout1[4 * (size_t)e]) * cx->c1->se[e].area;
  }
  int64_t ndark8 = 0;
  double adark8 = 0.0, atot2 = 0.0;
  for (int32_t e2 = 0; e2 < nse2; e2++) {
    double a2 = cut2.se[e2].area;
    if (!(a2 > 0.0)) continue;
    atot2 += a2;
    int32_t nv = cut2.se[e2].nv;
    if (nv < 3) {
      ndark8++;
      adark8 += a2;
      continue;
    }
    double acc0 = 0.0, accs[3] = {0, 0, 0}, wsum = 0.0;
    for (int q = 1; q + 1 < nv; q++) {
      const double *A = cut2.se[e2].v[0], *B = cut2.se[e2].v[q], *C = cut2.se[e2].v[q + 1];
      double u[3], w[3], nx[3];
      for (int k = 0; k < 3; k++) {
        u[k] = B[k] - A[k];
        w[k] = C[k] - A[k];
      }
      nx[0] = u[1] * w[2] - u[2] * w[1];
      nx[1] = u[2] * w[0] - u[0] * w[2];
      nx[2] = u[0] * w[1] - u[1] * w[0];
      double ar = 0.5 * sqrt(nx[0] * nx[0] + nx[1] * nx[1] + nx[2] * nx[2]);
      if (!(ar > 0.0)) continue;
      double p[3];
      for (int k = 0; k < 3; k++)
        p[k] = (A[k] + B[k] + C[k]) / 3.0;
      int32_t c1i = warm808_cell_at(cx->ot, cx->m1, cx->ofr, p);
      if (c1i < 0) continue;
      int32_t e1 = -1;
      for (int32_t kk = cx->c1->sestart[c1i]; kk < cx->c1->sestart[c1i + 1]; kk++) {
        int32_t ej = cx->c1->selist[kk];
        if (cx->c1->se[ej].facet == cut2.se[e2].facet) {
          e1 = ej;
          break;
        }
      }
      if (e1 < 0) continue;
      double ratio =
          (double)cx->m1->csize[cx->c1->se[e1].cell] / (double)mesh2.csize[cut2.se[e2].cell];
      acc0 += ar * cx->sout1[4 * (size_t)e1];
      for (int k = 1; k < 4; k++)
        accs[k - 1] += ar * cx->sout1[4 * (size_t)e1 + (size_t)k] * ratio;
      wsum += ar;
    }
    if (wsum > 0.0) {
      sout2[4 * (size_t)e2] = acc0 / wsum;
      for (int k = 1; k < 4; k++)
        sout2[4 * (size_t)e2 + (size_t)k] = accs[k - 1] / wsum;
    } else {
      ndark8++;
      adark8 += a2;
    }
  }
  double spow2 = 0.0, spowabs2 = 0.0;
  for (int32_t e2 = 0; e2 < nse2; e2++) {
    if (!(cut2.se[e2].area > 0.0)) continue;
    spow2 += sout2[4 * (size_t)e2] * cut2.se[e2].area;
    spowabs2 += fabs(sout2[4 * (size_t)e2]) * cut2.se[e2].area;
  }
  double t_rm = now_s() - trm0;
  double poured = 100.0 * (double)(npr8 + nre8) / (double)(mesh2.ncell > 0 ? mesh2.ncell : 1);
  {
    int64_t tot8 = nid8 + npr8 + nre8, acc8 = 0;
    int p998 = 0, pset8 = 0;
    for (int d = 0; d <= 6; d++) {
      int64_t c8 = d == 0 ? dhist[6] : dhist[6 - d] + dhist[6 + d];
      acc8 += c8;
      if (!pset8 && acc8 >= (tot8 * 99) / 100) {
        p998 = d;
        pset8 = 1;
      }
    }
    if (!pset8) p998 = 6;
    printf(
        "   §808 РЕМАП: φ — тождество %lld, продление %lld, рестрикция %lld (перелито %.1f %%, "
        "глубина p99 %d ур., гистограмма d [−6..+6]: %lld %lld %lld %lld %lld %lld | %lld | "
        "%lld %lld %lld %lld %lld %lld, переполнений сбора %lld); Σ|φ|·V %.6g -> %.6g (Δ %+.3f %%; "
        "знаковая %.6g -> %.6g); sout — без источника %lld элементов (%.2f %% площади), "
        "Σ sout·area %.6g -> %.6g (Δ %+.3f %%), Σ|sout|·area %.6g -> %.6g (Δ %+.3f %%); %.2f с\n",
        (long long)nid8, (long long)npr8, (long long)nre8, poured, p998, (long long)dhist[0],
        (long long)dhist[1], (long long)dhist[2], (long long)dhist[3], (long long)dhist[4],
        (long long)dhist[5], (long long)dhist[6], (long long)dhist[7], (long long)dhist[8],
        (long long)dhist[9], (long long)dhist[10], (long long)dhist[11], (long long)dhist[12],
        (long long)nov8, massabs1, massabs2,
        100.0 * (massabs2 - massabs1) / (massabs1 > 0.0 ? massabs1 : 1.0), mass1, mass2,
        (long long)ndark8, 100.0 * adark8 / (atot2 > 0.0 ? atot2 : 1.0), spow1, spow2,
        100.0 * (spow2 - spow1) / (fabs(spow1) > 0.0 ? fabs(spow1) : 1.0), spowabs1, spowabs2,
        100.0 * (spowabs2 - spowabs1) / (spowabs1 > 0.0 ? spowabs1 : 1.0), t_rm);
  }
  /* ---- ТЁПЛЫЙ ПУТЬ: прокидки К76 поверх ремапа + замыкание ---- */
  tr3_problem pw = *cx->prob;
  pw.m = &mesh2;
  pw.cut = &cut2;
  pw.facet_rho = frho2;
  pw.elem_emit = eemit2;
  pw.sig_t = sigt2;
  pw.sig_s = sigs2;
  const int nw = cx->nwarm;
  double *phip = malloc(nph2 * sizeof *phip);
  double *sop = malloc(nso2 * sizeof *sop);
  double *dphi = malloc(nph2 * sizeof *dphi);
  double *dso = malloc(nso2 * sizeof *dso);
  double *dbo = malloc(nbo2 * sizeof *dbo);
  if (phip == NULL || sop == NULL || dphi == NULL || dso == NULL || dbo == NULL) exit(1);
  memcpy(phip, phi2, nph2 * sizeof *phip);
  memcpy(sop, sout2, nso2 * sizeof *sop);
  for (size_t i = 0; i < nbo2; i++)
    dbo[i] = 0.0; /* байт-нулевая база: bout первой прокидкой восстанавливается */
  double *ubw = NULL, *usw = NULL, *eirr_w = NULL, *eirr_wp = NULL;
  double *bop = NULL;
  if (nw >= 2) {
    bop = calloc(nbo2, sizeof *bop);
    if (bop == NULL) exit(1);
  }
  double qnumw = 0.0, qdenw = 0.0, denprevw = 0.0;
  double psinw = 0.0, psinwp = 0.0;
  for (int t = 1; t <= nw; t++) {
    tr3_problem pt = pw;
    pt.warm_start = 1;
    pt.bout_in = t > 1 ? ubw : NULL; /* bout НЕ ремапится (А1348) */
    pt.sout_in = t > 1 ? usw : sout2;
    tr3_stats stt;
    memset(&stt, 0, sizeof stt);
    double tt808 = now_s();
    if (tr3_sweep_solve(&pt, 1, 0.0, phi2, &stt) != 0) exit(1);
    double num_t = 0.0, dencur = 0.0, dmx8 = 0.0;
    for (size_t i = 0; i < nph2; i++) {
      double dc = phi2[i] - phip[i];
      if (t > 1) num_t += dc * dphi[i];
      if (fabs(dc) > dmx8) dmx8 = fabs(dc);
      dencur += dc * dc;
      dphi[i] = dc;
      phip[i] = phi2[i];
    }
    for (size_t i = 0; i < nso2; i++) {
      double dc = stt.sout[i] - sop[i];
      if (t > 1) num_t += dc * dso[i];
      dencur += dc * dc;
      dso[i] = dc;
      sop[i] = stt.sout[i];
    }
    if (nw >= 2 && t > 1)
      for (size_t i = 0; i < nbo2; i++) {
        double dc = stt.bout[i] - bop[i];
        num_t += dc * dbo[i];
        dencur += dc * dc;
        dbo[i] = dc;
        bop[i] = stt.bout[i];
      }
    else if (nw >= 2)
      for (size_t i = 0; i < nbo2; i++) {
        dbo[i] = stt.bout[i]; /* первая прокидка bout строит, не приращивает */
        bop[i] = stt.bout[i];
      }
    if (t > 1) {
      qnumw = num_t;
      qdenw = denprevw;
    }
    denprevw = dencur;
    free(ubw);
    free(usw);
    free(eirr_wp);
    eirr_wp = eirr_w;
    ubw = stt.bout;
    usw = stt.sout;
    eirr_w = stt.eirr;
    psinwp = psinw;
    psinw = stt.psin;
    printf("   §808 тёплый такт %d: psin %.6g, |Δφ|∞ %.3e, %.2f с\n", t, stt.psin, dmx8,
           now_s() - tt808);
  }
  (void)psinwp;
  /* Замыкание: прирост последней прокидки × q̂/(1−q̂). q̂ при nw ≥ 2 измеряется
   * по двум приростам; при nw = 1 наследуется у лестницы кадра 1 — оператор
   * (сцена, альбедо) не менялся, скаляр сжатия глобальный. bout/eirr/psin при
   * nw = 1 не замыкаются: их прирост НЕ в режиме сжатия (bout строится с нуля,
   * у eirr нет базиса), замыкание завышало бы их в q̂/(1−q̂) раз. */
  double qw = -1.0, multw = 0.0;
  const char *qsrcw = "нет (чистое усечение)";
  if (nw >= 2 && qdenw > 0.0) qw = qnumw / qdenw;
  if (qw > 0.0 && qw <= HZ_TAILQ_MAX) {
    multw = qw / (1.0 - qw);
    qsrcw = "измерен (тёплый)";
  } else if (cx->qhat1 > 0.0 && cx->qhat1 <= HZ_TAILQ_MAX) {
    multw = cx->qhat1 / (1.0 - cx->qhat1);
    qsrcw = "унаследован у кадра 1";
  }
  if (multw > 0.0) {
    for (size_t i = 0; i < nph2; i++)
      phi2[i] += multw * dphi[i];
    for (size_t i = 0; i < nso2; i++)
      usw[i] += multw * dso[i];
    if (nw >= 2) {
      for (size_t i = 0; i < nbo2; i++)
        ubw[i] += multw * dbo[i];
      if (eirr_w != NULL && eirr_wp != NULL)
        for (int32_t e = 0; e < nse2; e++)
          eirr_w[e] += multw * (eirr_w[e] - eirr_wp[e]);
      psinw += multw * (psinw - psinwp);
    }
  }
  double t_warm = now_s() - trm0; /* ремап + прокидки + замыкание (А1346-б) */
  printf("   §808 ТЁПЛЫЙ: %d %s, замыкание q̂ %.4f (%s), psin %.6g, %.2f с (из них ремап %.2f с)\n",
         nw, warm808_plu(nw, "прокидка", "прокидки", "прокидок"),
         qw > 0.0 ? qw : (cx->qhat1 > 0.0 ? cx->qhat1 : 0.0), qsrcw, psinw, t_warm, t_rm);
  /* ---- ХОЛОДНАЯ БАЗА: та же лестница на сетке кадра 2, с нуля ---- */
  double tb0 = now_s();
  double *phiB = calloc(nph2, sizeof *phiB);
  double *pphiB = calloc(nph2, sizeof *pphiB);
  double *pboB = calloc(nbo2, sizeof *pboB);
  double *psoB = calloc(nso2, sizeof *psoB);
  double *dpphiB = calloc(nph2, sizeof *dpphiB);
  double *dpboB = calloc(nbo2, sizeof *dpboB);
  double *dpsoB = calloc(nso2, sizeof *dpsoB);
  double *peirrB = calloc((size_t)(nse2 > 0 ? nse2 : 1), sizeof *peirrB);
  if (phiB == NULL || pphiB == NULL || pboB == NULL || psoB == NULL || dpphiB == NULL ||
      dpboB == NULL || dpsoB == NULL || peirrB == NULL)
    exit(1);
  double *ubB = NULL, *usB = NULL;
  double qnumB = 0.0, qdenB = 0.0, denprevB = 0.0;
  double scprevB[6] = {0, 0, 0, 0, 0, 0};
  tr3_stats stB;
  memset(&stB, 0, sizeof stB);
  for (int t = 1; t <= cx->xbounce; t++) {
    tr3_problem pl = pw;
    if (t > 1) {
      pl.warm_start = 1;
      pl.bout_in = ubB;
      pl.sout_in = usB;
    }
    tr3_stats stt;
    memset(&stt, 0, sizeof stt);
    if (tr3_sweep_solve(&pl, 1, 0.0, phiB, &stt) != 0) exit(1);
    double num_t = 0.0, dencur = 0.0;
    for (size_t i = 0; i < nph2; i++) {
      double dc = phiB[i] - pphiB[i];
      num_t += dc * dpphiB[i];
      dencur += dc * dc;
      dpphiB[i] = dc;
      pphiB[i] = phiB[i];
    }
    for (size_t i = 0; i < nbo2; i++) {
      double dc = stt.bout[i] - pboB[i];
      num_t += dc * dpboB[i];
      dencur += dc * dc;
      dpboB[i] = dc;
      pboB[i] = stt.bout[i];
    }
    for (size_t i = 0; i < nso2; i++) {
      double dc = stt.sout[i] - psoB[i];
      num_t += dc * dpsoB[i];
      dencur += dc * dc;
      dpsoB[i] = dc;
      psoB[i] = stt.sout[i];
    }
    if (t > 1) {
      qnumB = num_t;
      qdenB = denprevB;
    }
    denprevB = dencur;
    free(ubB);
    free(usB);
    ubB = stt.bout;
    usB = stt.sout;
    if (t < cx->xbounce) {
      if (stt.eirr != NULL)
        for (int32_t e = 0; e < nse2; e++)
          peirrB[e] = stt.eirr[e];
      scprevB[0] = stt.pin;
      scprevB[1] = stt.pout;
      scprevB[2] = stt.pabs;
      scprevB[3] = stt.psin;
      scprevB[4] = stt.psout;
      scprevB[5] = stt.psolid;
      free(stt.eirr);
    } else
      stB = stt;
  }
  double qB = -1.0;
  if (cx->xtailq > 0.0)
    qB = cx->xtailq;
  else if (cx->xtailq < 0.0 && cx->xbounce >= HZ_TAILN_MIN && qdenB > 0.0)
    qB = qnumB / qdenB;
  if (qB > 0.0 && qB <= HZ_TAILQ_MAX) {
    double mult = qB / (1.0 - qB);
    for (size_t i = 0; i < nph2; i++)
      phiB[i] += mult * dpphiB[i];
    for (size_t i = 0; i < nbo2; i++)
      ubB[i] += mult * dpboB[i];
    for (size_t i = 0; i < nso2; i++)
      usB[i] += mult * dpsoB[i];
    if (stB.eirr != NULL)
      for (int32_t e = 0; e < nse2; e++)
        stB.eirr[e] += mult * (stB.eirr[e] - peirrB[e]);
    stB.pin += mult * (stB.pin - scprevB[0]);
    stB.pout += mult * (stB.pout - scprevB[1]);
    stB.pabs += mult * (stB.pabs - scprevB[2]);
    stB.psin += mult * (stB.psin - scprevB[3]);
    stB.psout += mult * (stB.psout - scprevB[4]);
    stB.psolid += mult * (stB.psolid - scprevB[5]);
  }
  double t_base = now_s() - tb0; /* прокидки + замыкание (А1346-в) */
  printf("   §808 БАЗА: %d %s с нуля, замыкание q̂ %.4f (%s), psin %.6g, %.2f с\n", cx->xbounce,
         warm808_plu(cx->xbounce, "прокидка", "прокидки", "прокидок"), qB > 0.0 ? qB : 0.0,
         cx->xtailq > 0.0 ? "ФОРСИРОВАН — НК" : (qB > 0.0 ? "измерен" : "нет — усечение"), stB.psin,
         t_base);
  /* ---- МЕТРИКИ: зона от НОВОГО глаза, поэлементный хвост ---- */
  double dps = stB.psin > 0.0 ? 100.0 * (psinw - stB.psin) / stB.psin : 0.0;
  double zw8 = 0.0, zb8 = 0.0;
  int64_t nzone8 = 0, nuse8 = 0;
  double *rel8 = malloc((size_t)(nse2 > 0 ? nse2 : 1) * sizeof *rel8);
  if (rel8 == NULL) exit(1);
  for (int32_t e = 0; e < nse2; e++) {
    double a9 = cut2.se[e].area;
    if (!(a9 > 0.0) || eirr_w == NULL || stB.eirr == NULL) continue;
    double c9[3] = {0, 0, 0};
    for (int q = 0; q < cut2.se[e].nv; q++)
      for (int a = 0; a < 3; a++)
        c9[a] += cut2.se[e].v[q][a] / (double)(cut2.se[e].nv > 0 ? cut2.se[e].nv : 1);
    double d2 = 0.0, dot = 0.0;
    for (int a = 0; a < 3; a++) {
      double dd = c9[a] - eye2[a];
      d2 += dd * dd;
      dot += cut2.se[e].n[a] * dd;
    }
    int inzone = sqrt(d2) < 15.0 && !(dot > 0.0);
    double ww8 = fabs(eirr_w[e]) * a9, wb8 = fabs(stB.eirr[e]) * a9;
    if (inzone) {
      zw8 += ww8;
      zb8 += wb8;
      nzone8++;
    }
    if (fabs(stB.eirr[e]) > 0.0) rel8[nuse8++] = fabs(eirr_w[e] - stB.eirr[e]) / fabs(stB.eirr[e]);
  }
  double dz = 100.0 * (zw8 - zb8) / (zb8 > 0.0 ? zb8 : 1.0);
  printf("   §808 ИНКРЕМЕНТНОСТЬ: холод %.2f с (%d %s), тёплый %.2f с (%d %s + ремап %.2f с, "
         "перелито %.1f %%), Δpsin %+.3f %%, Δ зоны %+.3f %%; заливка %s\n",
         t_base, cx->xbounce, warm808_plu(cx->xbounce, "прокидка", "прокидки", "прокидок"), t_warm,
         nw, warm808_plu(nw, "прокидка", "прокидки", "прокидок"), t_rm, poured, dps, dz,
         "полость та же");
  if (nuse8 > 0) {
    qsort(rel8, (size_t)nuse8, sizeof *rel8, cmp_dev699);
    printf("   §808 ПОЭЛЕМЕНТНО |Δeirr|/eirr (покрытие %lld из %d, база > 0): медиана %.4g, "
           "p90 %.4g, p99 %.4g, макс %.4g\n",
           (long long)nuse8, nse2, rel8[nuse8 / 2], rel8[(nuse8 * 9) / 10],
           rel8[(nuse8 * 99) / 100], rel8[nuse8 - 1]);
  }
  printf("   §808 ЗОНА: Σ|E·area| (<15 м, лицевые, %lld элементов) %.6g против %.6g\n",
         (long long)nzone8, zw8, zb8);
  if (zero808) {
    /* НК-3: сетки совпадают — база обязана быть ПОБИТОВО лестницей кадра 1
     * (одно и то же вычисление), ремап — тождеством; тёплый путь отличается
     * ровно на одну лишнюю эффективную итерацию (аудит исполнения скажет,
     * посимвольно ли это равенство — план §808 НК-3 утверждал «да»). */
    int bitsame = mesh2.ncell == cx->m1->ncell &&
                  memcmp(phiB, cx->phi1, (size_t)cx->m1->ncell * 4 * sizeof(double)) == 0;
    printf("   §808 НК-3 (сдвиг 0): ремап перелито %.1f %%, БАЗА == КАДР 1 побитово (φ): %s\n",
           poured, bitsame ? "ДА" : "нет");
  }
  printf("   §808 СТЕНД всего: %.2f с (общая часть %.2f, тёплый %.2f, база %.2f)\n", now_s() - t808,
         t_gc, t_warm, t_base);
  free(rel8);
  free(phi2);
  free(sout2);
  free(buf8);
  free(phip);
  free(sop);
  free(dphi);
  free(dso);
  free(dbo);
  free(bop);
  free(ubw);
  free(usw);
  free(eirr_w);
  free(eirr_wp);
  free(phiB);
  free(pphiB);
  free(pboB);
  free(psoB);
  free(dpphiB);
  free(dpboB);
  free(dpsoB);
  free(peirrB);
  free(ubB);
  free(usB);
  free(stB.eirr);
  free(eemit2);
  free(frho2);
  free(sigt2);
  free(sigs2);
  free(solid2);
  tr3_cut_free(&cut2);
  tr3_mesh_free(&mesh2);
}

/* ---- §810: РАСПИСАНИЕ ОТСКОКОВ — лестница на СМЕШАННЫХ ординатах -----------
 *
 * Директива §768-приписки: «отскок n глаже отскока n−1 ⇒ грубить с номером;
 * агрессивность от ШИРИНЫ ЛЕПЕСТКА; первый отскок не грубить (тени)». Ось
 * этого шага — УГЛОВАЯ: такт 1 идёт набором d1 (рабочие ординаты), такты ≥2 —
 * набором d2 (реже: после одного диффузного отскока угловой спектр — свёртка
 * с ламбертовым лепестком). Состояние (φ/bout/sout) угловой оси НЕ ИМЕЕТ
 * (φ — скалярный поток), наборы меняются между тактами без преобразований.
 * q̂ по двум последним приращениям: при N ≥ 3 оба — из ОДНОГО оператора S_{d2}
 * (А1356); переходная Δ_1 в оценку не входит. Замыкание — то же, что §774.
 * ПРОСТРАНСТВЕННАЯ ось (грейд на отскок) закрыта экономикой плана §810:
 * пересборка грейда [1.0, 1.3] с против экономии такта ≤ 0.5 с. */
static void sched810_ladder(const tr3_problem *pb, const tr3_dirs *d1, const tr3_dirs *d2, int nmu1,
                            int nmu2, int nb, double xtailq, double *phi, tr3_stats *stout,
                            double *qout, double *tout) {
  const int32_t nse4 = pb->cut != NULL ? pb->cut->nse : 0;
  const size_t nph = (size_t)pb->m->ncell * 4, nbo = (size_t)pb->m->nf * 4,
               nso = (size_t)(nse4 > 0 ? nse4 : 1) * 4;
  double *pphi = calloc(nph, sizeof *pphi);
  double *pbo = calloc(nbo, sizeof *pbo);
  double *pso = calloc(nso, sizeof *pso);
  double *dpphi = calloc(nph, sizeof *dpphi);
  double *dpbo = calloc(nbo, sizeof *dpbo);
  double *dpso = calloc(nso, sizeof *dpso);
  double *peirr = calloc((size_t)(nse4 > 0 ? nse4 : 1), sizeof *peirr);
  if (pphi == NULL || pbo == NULL || pso == NULL || dpphi == NULL || dpbo == NULL || dpso == NULL ||
      peirr == NULL)
    exit(1);
  double *ub = NULL, *us = NULL;
  double qnum = 0.0, qden = 0.0, denprev = 0.0;
  double scprev[6] = {0, 0, 0, 0, 0, 0};
  memset(stout, 0, sizeof *stout);
  double t0 = now_s();
  for (int t = 1; t <= nb; t++) {
    tr3_problem pl = *pb;
    pl.d = t == 1 ? d1 : d2; /* §810: первый отскок не грубится */
    if (t > 1) {
      pl.warm_start = 1;
      pl.bout_in = ub;
      pl.sout_in = us;
    }
    tr3_stats stt;
    memset(&stt, 0, sizeof stt);
    double tt = now_s();
    if (tr3_sweep_solve(&pl, 1, 0.0, phi, &stt) != 0) exit(1);
    double num_t = 0.0, dencur = 0.0, dmx = 0.0;
    for (size_t i = 0; i < nph; i++) {
      double dc = phi[i] - pphi[i];
      num_t += dc * dpphi[i];
      if (fabs(dc) > dmx) dmx = fabs(dc);
      dencur += dc * dc;
      dpphi[i] = dc;
      pphi[i] = phi[i];
    }
    for (size_t i = 0; i < nbo; i++) {
      double dc = stt.bout[i] - pbo[i];
      num_t += dc * dpbo[i];
      dencur += dc * dc;
      dpbo[i] = dc;
      pbo[i] = stt.bout[i];
    }
    for (size_t i = 0; i < nso; i++) {
      double dc = stt.sout[i] - pso[i];
      num_t += dc * dpso[i];
      dencur += dc * dc;
      dpso[i] = dc;
      pso[i] = stt.sout[i];
    }
    if (t > 1) {
      qnum = num_t;
      qden = denprev;
    }
    denprev = dencur;
    printf("   §810 такт %2d (nmu %d, ND %d): |Δφ|∞ %.3e, psin %.6g, %.2f с\n", t,
           t == 1 ? nmu1 : nmu2, t == 1 ? d1->n : d2->n, dmx, stt.psin, now_s() - tt);
    free(ub);
    free(us);
    ub = stt.bout;
    us = stt.sout;
    if (t < nb) {
      if (stt.eirr != NULL)
        for (int32_t e = 0; e < nse4; e++)
          peirr[e] = stt.eirr[e];
      scprev[0] = stt.pin;
      scprev[1] = stt.pout;
      scprev[2] = stt.pabs;
      scprev[3] = stt.psin;
      scprev[4] = stt.psout;
      scprev[5] = stt.psolid;
      free(stt.eirr);
    } else
      *stout = stt;
  }
  *tout = now_s() - t0;
  if (ub == NULL || us == NULL) exit(1); /* nb ≥ 1: решатель обязан оставить состояние */
  double q = -1.0;
  if (xtailq > 0.0)
    q = xtailq;
  else if (xtailq < 0.0 && nb >= HZ_TAILN_MIN && qden > 0.0)
    q = qnum / qden;
  if (q > 0.0 && q <= HZ_TAILQ_MAX) {
    double mult = q / (1.0 - q);
    for (size_t i = 0; i < nph; i++)
      phi[i] += mult * dpphi[i];
    for (size_t i = 0; i < nbo; i++)
      ub[i] += mult * dpbo[i];
    for (size_t i = 0; i < nso; i++)
      us[i] += mult * dpso[i];
    if (stout->eirr != NULL)
      for (int32_t e = 0; e < nse4; e++)
        stout->eirr[e] += mult * (stout->eirr[e] - peirr[e]);
    stout->pin += mult * (stout->pin - scprev[0]);
    stout->pout += mult * (stout->pout - scprev[1]);
    stout->pabs += mult * (stout->pabs - scprev[2]);
    stout->psin += mult * (stout->psin - scprev[3]);
    stout->psout += mult * (stout->psout - scprev[4]);
    stout->psolid += mult * (stout->psolid - scprev[5]);
  }
  *qout = q;
  printf("   §810 ЛЕСТНИЦА: %d тактов за %.2f с, q̂ %.4f (%s)\n", nb, *tout, q > 0.0 ? q : 0.0,
         xtailq > 0.0 ? "ФОРСИРОВАН — НК" : (q > 0.0 ? "измерен" : "нет — усечение"));
  free(pphi);
  free(pbo);
  free(pso);
  free(dpphi);
  free(dpbo);
  free(dpso);
  free(peirr);
}

static int32_t g_xcmp_n = 0, g_xcmp_nemit = 0;
static double *g_xcmp_E = NULL, *g_xcmp_Edir = NULL, *g_xcmp_aw = NULL;
static int32_t (*g_xcmp_lo)[3] = NULL;
static int32_t *g_xcmp_sz = NULL;
static double (*g_xcmp_ec)[3] = NULL;
/* НК1 §597: прежний путь — сбор в срезе, каждый кадр. */
static int g_indslice = 0;
/* НК2 §597: подъём по иерархии выключен. Крупные узлы получают ноль, и крупные
 * ячейки кадра обязаны почернеть — это контроль ЛОЖНОГО НУЛЯ (А917). */
static int g_indnolift = 0;
/* НК3 §598/А1028: подъём НЕВЗВЕШЕННЫЙ. Если разницы с взвешенным нет, значит
 * площади у детей одинаковы, и это надо ЗНАТЬ, а не предполагать. */
static int g_indflat = 0;
/* §607: число отскоков в ядре (ряд Неймана). Потолок назван константой, чтобы
 * массив сумм не заводился по вводу пользователя. */
#define HZ_BOUNCE_MAX 8
static int g_hbounce = 1;
/* Ф8' (§532): отскок направленным свипом. `nmu,nphi` — набор §2 (1,1 / 2,2 /
 * 2,4 / 4,4 даёт ND = 8, 32, 64, 128); `dpass` — число отскоков (довод №3 §523:
 * матрицы нет, есть ещё один проход); `dirsall` — НЕГАТИВНЫЙ КОНТРОЛЬ, излучать
 * во все стороны вместо наружной полусферы. */

/* СБОР ПО ИЕРАРХИИ ИЗЛУЧАТЕЛЕЙ — ОДНА РЕАЛИЗАЦИЯ НА ОБА ПУТИ (§597).
 *
 * Вынесено из кадрового блока БЕЗ изменения смысла: те же `hgather_rec`, то же
 * деление по приёмникам, тот же фиксированный порядок сведения счётчиков.
 * Копии не заводится СОЗНАТЕЛЬНО: копия — это два места, где чинить один и тот
 * же промах, а прежний путь остаётся эталоном (`0.7460` на комнате).
 *
 * `ind` — выход, `3` float на ячейку среза; `nlink`, `sthru`, `sall` могут быть
 * NULL. Альбедо приёмника применяется здесь же, как и было. */
static void gather_run(const etree *ET, const hz_dcslice *S, const frame *fr, const opyr *P,
                       const hz_objmesh *m, int lev, double eps, int blockvis, int alb0, float *ind,
                       int64_t *nlink, double *sthru, double *sall, int64_t *nnear, double *ffout,
                       linkcache *lc) {
  /* §621: МАСШТАБ СЦЕНЫ ДЛЯ ПОРОГА ВКЛАДА — средневзвешенная по площади
   * радиосить излучателей, взятая с КОРНЯ дерева. Порог задаётся её долей, а не
   * абсолютным числом: абсолютных порогов проект не допускает, и здесь величина
   * выводится из самой сцены. */
  double bsum = 0.0, asum = 0.0;
  for (int k = 0; k < ET->e[0].nb; k++) {
    const ebin *rb = &ET->b[ET->e[0].b0 + k];
    asum += (double)rb->area;
    for (int c = 0; c < 3; c++)
      bsum += (double)rb->flux[c] / 3.0;
  }
  double tau = g_hcontrib * (asum > 0.0 ? bsum / asum : 0.0);
  int nth = g_omp1 ? 1 : omp_get_max_threads();
  int64_t *plink = calloc((size_t)nth, sizeof *plink);
  double *pthru = calloc((size_t)nth, sizeof *pthru);
  double *pall = calloc((size_t)nth, sizeof *pall);
  int64_t *pnear = calloc((size_t)nth, sizeof *pnear);
  int64_t *pmar = calloc((size_t)nth, sizeof *pmar);
  if (plink == NULL || pthru == NULL || pall == NULL || pnear == NULL || pmar == NULL) exit(1);
#pragma omp parallel for schedule(dynamic, 64) if (!g_omp1)
  for (int32_t i = 0; i < S->n; i++) {
    int th = g_omp1 ? 0 : omp_get_thread_num();
    double pi[3], nn2[3], acc2[3] = {0, 0, 0}, ffacc = 0.0;
    int64_t lc0 = lc != NULL ? lc->n[th] : 0;
    hz_slice_vertex(S, i, pi);
    for (int k = 0; k < 3; k++)
      pi[k] = fr->org[k] + pi[k] * fr->h;
    hz_slice_normal(S, i, nn2);
    double rrecv = fr->h * (double)((int32_t)1 << (lev - (int)S->c[i].lvl));
    for (int q2 = 0; q2 < 6; q2++)
      hgather_rec(ET, 0, q2, pi, nn2, eps, tau, rrecv, P, fr, blockvis, acc2, &plink[th],
                  &pthru[th], &pall[th], &pnear[th], &ffacc, lc, th, &pmar[th]);
    /* §611: НОРМИРОВКА. `Σ_j F_ij` физически не больше единицы; замерено, что у
     * `68.9 %` приёмников она больше (§610, медиана `1.8367`). Деление на
     * `max(1, Σ F)` делает оператор СЖАТИЕМ по построению.
     * ЧЕСТНО: это заставляет энергию сохраняться, но перекрытие шести корзин НЕ
     * чинит — ошибка перераспределяется по направлениям. Диагностика `Σ F`
     * печатается ПОСЛЕ нормировки и обязана показывать прежние `1.84`. */
    double fnorm = (!g_gnonorm && ffacc > 1.0) ? 1.0 / ffacc : 1.0;
    for (int k = 0; k < 3; k++)
      ind[3 * (size_t)i + (size_t)k] =
          (float)(acc2[k] * fnorm * (alb0 ? 0.0 : alb(m, S->c[i].mat, k)));
    if (ffout != NULL) ffout[i] = ffacc;
    if (lc != NULL) {
      lc->off[i] = lc0;
      lc->cnt[i] = (int32_t)(lc->n[th] - lc0);
      lc->th[i] = (int8_t)th;
    }
  }
  /* Счётчики сводятся в ФИКСИРОВАННОМ порядке: редукция OpenMP отдала бы его
   * планировщику, и число поехало бы от запуска к запуску. */
  for (int t4 = 0; t4 < nth; t4++) {
    if (nlink != NULL) *nlink += plink[t4];
    if (nnear != NULL) *nnear += pnear[t4];
    g_march_steps += pmar[t4];
    if (sthru != NULL) *sthru += pthru[t4];
    if (sall != NULL) *sall += pall[t4];
  }
  free(plink);
  free(pthru);
  free(pall);
  free(pnear);
  free(pmar);
}

/* ПРИМЕНИТЬ ХРАНИМЫЕ СВЯЗИ — БЕЗ ОБХОДА ДЕРЕВА (§613).
 *
 * `E_i = Σ g·flux[корзина]`, затем нормировка §611 и альбедо приёмника — всё
 * ровно как в `gather_run`, иначе первый и второй отскок считались бы разными
 * формулами. Сумма формфакторов берётся ГОТОВОЙ (`ffv`): она геометрическая и
 * между отскоками не меняется. */
static void gather_apply(const linkcache *lc, const etree *T, const hz_dcslice *S,
                         const hz_objmesh *m, const double *ffv, float *ind) {
#pragma omp parallel for schedule(dynamic, 256) if (!g_omp1)
  for (int32_t i = 0; i < lc->nrecv; i++) {
    const glink *lk = lc->buf[lc->th[i]] + lc->off[i];
    int32_t nl = lc->cnt[i];
    double acc[3] = {0, 0, 0};
    for (int32_t j = 0; j < nl; j++) {
      const ebin *bb = &T->b[lk[j].bin];
      for (int k = 0; k < 3; k++)
        acc[k] += (double)bb->flux[k] * (double)lk[j].g;
    }
    double fnorm = (!g_gnonorm && ffv[i] > 1.0) ? 1.0 / ffv[i] : 1.0;
    for (int k = 0; k < 3; k++)
      ind[3 * (size_t)i + (size_t)k] = (float)(acc[k] * fnorm * alb(m, S->c[i].mat, k));
  }
}
/* ---- §796: FIRST-COLLISION SOURCE ----------------------------------------
 *
 * Прямой свет Ke-ламп до ЭЛЕМЕНТОВ РАЗВЁРТКИ считается НЕПРЕРЫВНЫМ по
 * направлениям сбором — той же машинерией, что ядро §597 (etree + hgather_rec
 * + заслон shadowed_h), а в развёртку инъецируется ПЕРВОЕ ОТРАЖЕНИЕ
 * ρ·E_fc/hsum; свип разносит только рассеянный хвост. Лечение А1286: «чешуя»
 * кадра L9 — ray effect дискретных ординат (компактная лампа рисует звёзды по
 * ND направлениям); у сбора направления непрерывны, и класс уходит на ЛЮБЫХ
 * кусках, включая дешёвые DC. Включается КЛЮЧОМ `xfc` (исследовательский:
 * приёмка §796-П3 провалена — чешуя при fc осталась, УБИВАЕТ сработал);
 * без ключа — прежняя инъекция эмиссией ПОСИМВОЛЬНО (НК §796).
 *
 * Марши заслона fc-сбора НЕ попадают в `g_march_steps` (А1290): та строка
 * §620 принадлежит ядру §597, и смешение изменило бы её в комбинированных
 * прогонах. Счётчик у fc свой и печатается своей строкой. */

/* Излучатель fc-дерева: светящийся элемент свипа. `flux` = π·L·area — ровно
 * слагаемое `emitpow2`, поэтому Σ flux по корню обязана сойтись с печатью
 * «Σ π·L·площадь»入 развёртки. `key` — мортон-код клетки (те же биты и тот же
 * порядок осей, что у детей etree_build), `idx` — разрыв ничьих: сортировка
 * становится тотальной, и порядок сложения в корзинах детерминирован. */
typedef struct {
  double c[3], n[3], flux, area;
  uint64_t key;
  int32_t idx;
} efcsrc;

static int efc_cmp(const void *a, const void *b) {
  const efcsrc *x = a, *y = b;
  if (x->key != y->key) return x->key < y->key ? -1 : 1;
  return (x->idx > y->idx) - (x->idx < y->idx);
}

/* Клон etree_build над массивом излучателей-элементов: те же корзины по
 * квадранту нормали, тот же спуск по битам мортон-кода. Копия сознательная:
 * etree_build читает hz_dcslice (ячейки среза), а здесь элементы разреза —
 * общей структуры входа у них нет, а менять §597-путь ради §796 значило бы
 * трогать эталон гладкости, которым §796 принимается. */
static int32_t efc_build(etree *T, const efcsrc *s, int32_t a, int32_t b, int lvl, int lev) {
  int32_t me = etree_alloc(T, 1);
  {
    double c[6][3] = {{0}}, n[6][3] = {{0}}, fl[6] = {0}, ar[6] = {0};
    double lo[6][3], hi[6][3];
    for (int q = 0; q < 6; q++)
      for (int k = 0; k < 3; k++) {
        lo[q][k] = 1e300;
        hi[q][k] = -1e300;
      }
    for (int32_t i = a; i < b; i++) {
      int ax = 0;
      for (int k = 1; k < 3; k++)
        if (fabs(s[i].n[k]) > fabs(s[i].n[ax])) ax = k;
      int q = 2 * ax + (s[i].n[ax] > 0.0 ? 1 : 0);
      double aa = s[i].area;
      ar[q] += aa;
      fl[q] += s[i].flux;
      for (int k = 0; k < 3; k++) {
        c[q][k] += s[i].c[k] * aa;
        n[q][k] += s[i].n[k] * aa;
        if (s[i].c[k] < lo[q][k]) lo[q][k] = s[i].c[k];
        if (s[i].c[k] > hi[q][k]) hi[q][k] = s[i].c[k];
      }
    }
    int nbq = 0;
    for (int q = 0; q < 6; q++)
      if (ar[q] > 0.0) nbq++;
    int32_t b0 = ebin_alloc(T, nbq);
    T->e[me].b0 = b0;
    T->e[me].nb = (signed char)nbq;
    T->e[me].nch = 0;
    for (int k = 0; k < 8; k++)
      T->e[me].ch[k] = -1;
    int j = 0;
    for (int q = 0; q < 6; q++) {
      if (!(ar[q] > 0.0)) continue;
      ebin *bb = &T->b[b0 + j++];
      memset(bb, 0, sizeof *bb);
      bb->q = (signed char)q;
      bb->area = (float)ar[q];
      double nl = 0.0;
      for (int k = 0; k < 3; k++) {
        bb->c[k] = (float)(c[q][k] / ar[q]);
        /* свип монохромный: три канала корзины несут одно и то же число */
        bb->flux[k] = (float)fl[q];
        double nk = n[q][k] / ar[q];
        nl += nk * nk;
      }
      nl = sqrt(nl);
      bb->nsum = (float)(nl * ar[q]);
      for (int k = 0; k < 3; k++)
        bb->n[k] = (float)(nl > 0.0 ? (n[q][k] / ar[q]) / nl : 0.0);
      double r2 = 0.0;
      for (int k = 0; k < 3; k++) {
        double d = 0.5 * (hi[q][k] - lo[q][k]);
        r2 += d * d;
      }
      bb->rad = (float)sqrt(r2);
    }
  }
  if (b - a <= 1 || lvl >= lev) return me;
  int sh = 3 * (lev - lvl - 1);
  int32_t bnd[9];
  bnd[0] = a;
  int32_t cur = a;
  for (int k = 1; k <= 8; k++) {
    while (cur < b && (int)((s[cur].key >> sh) & 7u) < k)
      cur++;
    bnd[k] = cur;
  }
  int nc2 = 0;
  int32_t ch2[8];
  for (int k = 0; k < 8; k++) {
    if (bnd[k + 1] <= bnd[k]) continue;
    ch2[nc2++] = efc_build(T, s, bnd[k], bnd[k + 1], lvl + 1, lev);
  }
  for (int k = 0; k < nc2; k++)
    T->e[me].ch[k] = ch2[k];
  T->e[me].nch = (signed char)nc2;
  return me;
}

/* Центроид элемента разреза В МИРЕ: среднее вершин многоугольника; при
 * nv = 0 (не поместился, считан в nsebig) — центр клетки. Возврат 0 — откат
 * к центру клетки, вызывающий его считает. */
static int efc_center(const tr3_cut *cu, const tr3_mesh *ms, const frame *fr, int32_t e,
                      double pi[3]) {
  const tr3_selem *se = &cu->se[e];
  if (se->nv > 0) {
    for (int a = 0; a < 3; a++) {
      pi[a] = 0.0;
      for (int q = 0; q < se->nv; q++)
        pi[a] += se->v[q][a];
      pi[a] /= (double)se->nv;
    }
    return 1;
  }
  for (int a = 0; a < 3; a++)
    pi[a] = fr->org[a] + ((double)ms->clo[se->cell][a] + 0.5 * (double)ms->csize[se->cell]) * fr->h;
  return 0;
}

/* Угловой допуск спуска fc-сбора. Читается только под `hangle` (НК-режим
 * прежнего углового критерия); рабочий спуск идёт по вкладу §621. Значение —
 * рабочее у ядра §597 (`hgather=0.5`, канон §793). */
#define HZ_FC_EPS 0.5

/* Порог вклада fc-сбора, доля средней радиосити ламп (§621). Умолчание —
 * КАНОН ядра §597 (`hcontrib=1e-1`, §793/§650-команды): при общем умолчании
 * 1e-3 сбор на конвейере D0=10 стоил 9.35 с (замерено §796, 168 М связей);
 * явный `hcontrib=` его перебивает. */
#define HZ_FC_CONTRIB 1e-1

/* Сбор E_fc по всем элементам. Возвращает счётчики фиксированным порядком
 * сведения (§770-дисциплина): связи, поток весь/заслонённый, шаги марша. */
static void efc_run(const etree *T, const tr3_cut *cu, const tr3_mesh *ms, const frame *fr,
                    const opyr *P, double hcontrib, double *efc, int64_t *nlink, double *sall,
                    double *sthru, int64_t *nmarch) {
  /* порог вклада §621: доля средневзвешенной радиосити ЛАМП с корня */
  double bsum = 0.0, asum = 0.0;
  for (int k = 0; k < T->e[0].nb; k++) {
    const ebin *rb = &T->b[T->e[0].b0 + k];
    asum += (double)rb->area;
    for (int c = 0; c < 3; c++)
      bsum += (double)rb->flux[c] / 3.0;
  }
  double tau = hcontrib * (asum > 0.0 ? bsum / asum : 0.0);
  int nth = g_omp1 ? 1 : omp_get_max_threads();
  int64_t *plink = calloc((size_t)nth, sizeof *plink);
  double *pthru = calloc((size_t)nth, sizeof *pthru);
  double *pall = calloc((size_t)nth, sizeof *pall);
  int64_t *pmar = calloc((size_t)nth, sizeof *pmar);
  if (plink == NULL || pthru == NULL || pall == NULL || pmar == NULL) exit(1);
#pragma omp parallel for schedule(dynamic, 64) if (!g_omp1)
  for (int32_t e = 0; e < cu->nse; e++) {
    int th = g_omp1 ? 0 : omp_get_thread_num();
    const tr3_selem *se = &cu->se[e];
    efc[e] = 0.0;
    if (!(se->area > 0.0)) continue;
    /* нулевая инициализация — гейт: gcc-analyzer теряет запись pi[a] в цикле
     * efc_center через границу OMP-замыкания (класс §794-FP), путь ложный */
    double pi[3] = {0, 0, 0};
    efc_center(cu, ms, fr, e, pi);
    /* А1289: протяжённость ПРИЁМНИКА — сам элемент (√area), не клетка:
     * на грубой сетке §782 клетка в метры заморозила бы дробление ламп. */
    double rrecv = sqrt(se->area);
    double acc[3] = {0, 0, 0}, ffacc = 0.0;
    for (int q = 0; q < 6; q++)
      hgather_rec(T, 0, q, pi, se->n, HZ_FC_EPS, tau, rrecv, P, fr, 1, acc, &plink[th], &pthru[th],
                  &pall[th], NULL, &ffacc, NULL, 0, &pmar[th]);
    /* нормировка §611 — как у сбора: оператор сжатие по построению */
    double fnorm = (!g_gnonorm && ffacc > 1.0) ? 1.0 / ffacc : 1.0;
    efc[e] = acc[0] * fnorm;
  }
  for (int t = 0; t < nth; t++) {
    *nlink += plink[t];
    *sall += pall[t];
    *sthru += pthru[t];
    *nmarch += pmar[t];
  }
  free(plink);
  free(pthru);
  free(pall);
  free(pmar);
}

/* §796: пары «ключ квантованной плоскости — элемент» прибора чешуи и их
 * порядок (по ключу, ничьи по элементу — сортировка тотальна, группы
 * детерминированы). */
typedef struct {
  uint64_t key;
  int32_t e;
} cvpair796;

static int cmp_cvpair796(const void *a, const void *b) {
  const cvpair796 *x = a, *y = b;
  if (x->key != y->key) return x->key < y->key ? -1 : 1;
  return (x->e > y->e) - (x->e < y->e);
}

static int cmp_i32_796(const void *a, const void *b) {
  int32_t x = *(const int32_t *)a, y = *(const int32_t *)b;
  return (x > y) - (x < y);
}

/* УЗЕЛ ПО ЯЧЕЙКЕ СРЕЗА (§597, Р3). Ячейка несёт `lo` в сетке САМОГО МЕЛКОГО
 * уровня и свой уровень `lvl`; номера узла у неё нет. Спуск целочисленный: на
 * шаге `d` бит ребёнка берётся из разряда `lev−1−d` координаты — та же нумерация
 * детей, что в `cell_proc` (бит `a` оси `a`).
 *
 * ОСТАНОВКА РОВНО НА ГЛУБИНЕ `lvl`, И ЭТО ТРУДНОЕ МЕСТО (А1029): ошибка на
 * единицу уводит в ребёнка, косвенный свет уезжает в соседнюю ячейку, а на
 * картинке это выглядит лёгким шумом — то есть НЕ бросается в глаза. Поэтому
 * возвращается `−1` при любом несоответствии, а вызывающий ОБЯЗАН считать
 * отказы и печатать их число. */
static int32_t node_of_cell(const hz_dctree *t, int lev, const hz_dccell *c) {
  int32_t ni = 0;
  int lvl = (int)c->lvl;
  if (lvl < 0 || lvl > lev) return -1;
  for (int d = 0; d < lvl; d++) {
    if (t->nd[ni].child0 < 0) return -1; /* дерево мельче, чем просит срез */
    int bit = 0;
    for (int a = 0; a < 3; a++)
      if (((int32_t)c->lo[a] >> (lev - 1 - d)) & 1) bit |= 1 << a;
    ni = t->nd[ni].child0 + bit;
  }
  return ni;
}

/* КАКОЙ ТРЕУГОЛЬНИК ПРИНАДЛЕЖИТ ЯЧЕЙКЕ СРЕЗА (§597; правило из Ш5б, §430).
 *
 * ОДНА РЕАЛИЗАЦИЯ НА ОБА ПУТИ — кадровый и ядро. Прежде правило было записано
 * прямо в кадровом блоке, и ядро §597 без него получило `mat = 0` у ВСЕХ ячеек:
 * `Σ E_ind` вышло `1.18` против `Σ E_dir = 2.7e5`, то есть косвенный свет был
 * НУЛЕВОЙ. Копию заводить нельзя — это два места, где чинить один промах.
 *
 * ИСКАТЬ НАДО ПО ВЕРШИНЕ, А НЕ ПО УГЛУ ЯЧЕЙКИ (замерено 08-11): угол лежит
 * где угодно — внутри тела, в пустоте, на соседнем предмете, — и материал
 * оттуда берётся чужой (назначено `31 %` ячеек, шар вышел пятнистым). Вершина
 * же лежит НА поверхности по построению DC.
 *
 * Возвращает индекс треугольника или `−1`. */
static int32_t slice_cell_tri(const hz_dcslice *S, int32_t i, celltris *CT, const frame *fr,
                              double vw[3]) {
  hz_slice_vertex(S, i, vw);
  int32_t cell[3];
  for (int a = 0; a < 3; a++) {
    double f = floor(vw[a]);
    if (f < 0.0) f = 0.0;
    if (f > (double)(fr->n - 1)) f = (double)(fr->n - 1);
    cell[a] = (int32_t)f;
  }
  const int32_t *ls = NULL;
  if (ct_list(CT, cell, &ls) == 0) return -1;
  return ls[0];
}

/* Материал ячейки по тому же правилу. Возвращает число назначенных. */
static int64_t slice_assign_mat(hz_dcslice *S, celltris *CT, const frame *fr, const hz_objmesh *m) {
  int64_t n = 0;
  for (int32_t i = 0; i < S->n; i++) {
    double vw[3];
    int32_t tri = slice_cell_tri(S, i, CT, fr, vw);
    if (tri < 0) continue;
    int32_t mi = m->fm != NULL ? m->fm[tri] : 0;
    if (mi < 0 || mi >= m->nmtl) mi = 0;
    S->c[i].mat = (uint8_t)(mi < 255 ? mi : 255);
    n++;
  }
  return n;
}

/* ПОДЪЁМ ПО ИЕРАРХИИ (§597, Р4; правка А1028).
 *
 * ЗАЧЕМ. Ядро считается на срезе ПОЛНОЙ глубины, а камера может взять КРУПНЫЙ
 * узел, которого в полном срезе нет вовсе. Без подъёма он получил бы ноль — а
 * «нет пометки» есть НЕИЗВЕСТНО, а не ноль (А917, четвёртый случай за проект).
 *
 * ВЗВЕШИВАЕТСЯ ЧИСЛОМ ЛИСТЬЕВ С ВЕРШИНОЙ, А НЕ ПРОСТЫМ СРЕДНИМ ПО ВОСЬМИ, И
 * РАЗНИЦА НЕ КОСМЕТИЧЕСКАЯ (А1028). Облучённость есть ПЛОТНОСТЬ, а не
 * аддитивная величина; у детей разное количество поверхности внутри, и простое
 * среднее по восьми дало бы крупному узлу не ту величину, которую он несёт.
 * Число листьев с вершиной — заместитель площади: точной площади под рукой нет,
 * и подменять её единицей было бы той же ошибкой, только молча.
 * `indflat` (НК3) возвращает невзвешенное среднее.
 *
 * Возвращает число листьев с вершиной в поддереве. */
static int32_t ind_lift(const hz_dctree *t, int32_t ni) {
  if (t->nd[ni].child0 < 0) {
    /* Лист: значение уже разложено (либо ноль, если вершины нет). */
    return hz_dc_hasvert(t, ni) ? 1 : 0;
  }
  double acc[3] = {0, 0, 0};
  double wsum = 0.0;
  int32_t nleaf = 0;
  for (int k = 0; k < 8; k++) {
    int32_t ci = t->nd[ni].child0 + k;
    int32_t nc = ind_lift(t, ci);
    nleaf += nc;
    if (nc == 0) continue;
    double w = g_indflat ? 1.0 : (double)nc;
    for (int a = 0; a < 3; a++)
      acc[a] += (double)g_indnode[ci][a] * w;
    wsum += w;
  }
  if (!g_indnolift && wsum > 0.0)
    for (int a = 0; a < 3; a++)
      g_indnode[ni][a] = (float)(acc[a] / wsum);
  return nleaf;
}

/* ЯДРО КОСВЕННОГО СВЕТА — СЧИТАЕТСЯ ОДИН РАЗ, ЖИВЁТ МЕЖДУ КАДРАМИ (§597).
 *
 * Приёмники — срез ПОЛНОЙ ГЛУБИНЫ (`stop = NULL`): он камеронезависим ПО
 * ОПРЕДЕЛЕНИЮ, а не по удачно выбранному порогу. Отвергнутая альтернатива
 * (считать на кадровом срезе первой камеры) записана в А1030: она оставила бы
 * скрытую привязку к точке запуска — «свет хороший там, откуда я вышел».
 *
 * ИСТОЧНИК ОБЯЗАН БЫТЬ КАМЕРОНЕЗАВИСИМЫМ, иначе всё это бессмысленно.
 * Проверено чтением (А1033): под `sun` в расчёт входят только точка приёмника,
 * фиксированное направление, диагональ сцены и постоянный цвет. А вот
 * `hall_light` ищет потолок полости СПУСКОМ ОТ КАМЕРЫ — с площадным источником
 * ядро несовместимо, и это записанный долг, а не забытое. */
/* §744: статистика одной пары полей — счётчики по прецеденту §559 (порог от
 * данных, выброшенное считается, А999), медиана отношения, доля вне ×2 от
 * медианы, корзины по расстоянию до ближайшего излучателя (кратные листу L —
 * не подбор), топ-3 расхождений в обе стороны. */
static void xcmp_stats(const char *tag, const double *A, const double *B, const double (*cen)[3],
                       int64_t n, double Lh) {
  double *tmp = malloc((size_t)(n > 0 ? n : 1) * sizeof *tmp);
  double *rt = malloc((size_t)(n > 0 ? n : 1) * sizeof *rt);
  double *dst = malloc((size_t)(n > 0 ? n : 1) * sizeof *dst);
  int64_t *idx = malloc((size_t)(n > 0 ? n : 1) * sizeof *idx);
  if (tmp == NULL || rt == NULL || dst == NULL || idx == NULL) exit(1);
  int64_t nb = 0;
  for (int64_t i = 0; i < n; i++)
    if (B[i] > 0.0) tmp[nb++] = B[i];
  if (nb == 0) {
    printf("   §744 %s: знаменатель пуст, сличать нечего\n", tag);
    free(tmp);
    free(rt);
    free(dst);
    free(idx);
    return;
  }
  qsort(tmp, (size_t)nb, sizeof *tmp, cmp_d);
  double thr = 1e-3 * tmp[nb / 2];
  int64_t nboth0 = 0, ngz = 0, nsz = 0, nuse = 0;
  for (int64_t i = 0; i < n; i++) {
    double a = A[i], b = B[i];
    if (!(a > 0.0) && !(b > thr)) {
      nboth0++;
      continue;
    }
    if (!(b > thr)) {
      ngz++;
      continue;
    }
    if (!(a > 0.0)) {
      nsz++;
      continue;
    }
    rt[nuse] = a / b;
    /* расстояние до ближайшего излучателя — по центрам излучающих ячеек */
    double dmin = 1e300;
    for (int32_t j = 0; j < g_xcmp_nemit; j++) {
      double d2 = 0.0;
      for (int k = 0; k < 3; k++) {
        double dd = cen[i][k] - g_xcmp_ec[j][k];
        d2 += dd * dd;
      }
      if (d2 < dmin) dmin = d2;
    }
    dst[nuse] = sqrt(dmin);
    idx[nuse] = i;
    nuse++;
  }
  if (nuse == 0) {
    printf("   §744 %s: пар нет (оба нуля %lld, ядро~0 %lld, свип~0 %lld)\n", tag,
           (long long)nboth0, (long long)ngz, (long long)nsz);
    free(tmp);
    free(rt);
    free(dst);
    free(idx);
    return;
  }
  memcpy(tmp, rt, (size_t)nuse * sizeof *tmp);
  qsort(tmp, (size_t)nuse, sizeof *tmp, cmp_d);
  double med = tmp[nuse / 2], p10 = tmp[(nuse * 10) / 100], p90 = tmp[(nuse * 90) / 100];
  int64_t nout2 = 0;
  for (int64_t i = 0; i < nuse; i++)
    if (rt[i] < 0.5 * med || rt[i] > 2.0 * med) nout2++;
  printf("   §744 %s: пар %lld (оба нуля %lld, ядро~0 %lld, свип~0 %lld); медиана r %.4g, "
         "p10 %.4g, p90 %.4g; ВНЕ x2 от медианы: %lld (%.1f %%)\n",
         tag, (long long)nuse, (long long)nboth0, (long long)ngz, (long long)nsz, med, p10, p90,
         (long long)nout2, 100.0 * (double)nout2 / (double)nuse);
  static const double BINL[4] = {2.0, 8.0, 32.0, 1e30};
  double blo = 0.0;
  for (int bq = 0; bq < 4; bq++) {
    double bhi = BINL[bq] < 1e29 ? BINL[bq] * Lh : 1e30;
    int64_t cnt = 0;
    for (int64_t i = 0; i < nuse; i++)
      if (dst[i] >= blo && dst[i] < bhi) tmp[cnt++] = rt[i];
    if (cnt > 0) {
      qsort(tmp, (size_t)cnt, sizeof *tmp, cmp_d);
      char hib[32];
      if (BINL[bq] < 1e29)
        snprintf(hib, sizeof hib, "%.1f", bhi);
      else
        snprintf(hib, sizeof hib, "inf");
      printf("      §744 %s корзина [%.1f, %s м): пар %lld, медиана r/med %.3f\n", tag, blo, hib,
             (long long)cnt, tmp[cnt / 2] / med);
    }
    blo = bhi;
  }
  for (int side = 0; side < 2; side++) {
    for (int t3 = 0; t3 < 3; t3++) {
      double best = -1.0;
      int64_t bi = -1;
      /* взятое помечается знаком (rt < 0) и исключается условием rt > 0 */
      for (int64_t i = 0; i < nuse; i++) {
        if (!(rt[i] > 0.0)) continue;
        double v = side == 0 ? rt[i] / med : med / rt[i];
        if (v > best) {
          best = v;
          bi = i;
        }
      }
      if (bi < 0) break;
      printf("      §744 %s топ-%d %s: r/med %.3g в (%.1f, %.1f, %.1f) м, до излучателя %.1f м\n",
             tag, t3 + 1, side == 0 ? "ВВЕРХ" : "ВНИЗ", side == 0 ? rt[bi] / med : med / rt[bi],
             cen[idx[bi]][0], cen[idx[bi]][1], cen[idx[bi]][2], dst[bi]);
      rt[bi] = -rt[bi]; /* пометка взятого: отрицательное исключается условием rt > 0 */
    }
    for (int64_t i = 0; i < nuse; i++)
      if (rt[i] < 0.0) rt[i] = -rt[i];
  }
  free(tmp);
  free(rt);
  free(dst);
  free(idx);
}

/* §744: мост — ячейки свипа к узлам ядра §597 спуском node_of_cell; три
 * сличения: главное, тождество (НК-а) и заведомо разные поля (НК-б). */
static void xcmp_report(const hz_dctree *T, int lev, const frame *fr) {
  if (g_xcmp_n <= 0 || g_indnode == NULL) return;
  int32_t nc = g_xcmp_n;
  double *Esw = malloc((size_t)nc * sizeof *Esw);
  double *Edir = malloc((size_t)nc * sizeof *Edir);
  double *Eg = malloc((size_t)nc * sizeof *Eg);
  double (*cen)[3] = malloc((size_t)nc * sizeof *cen);
  int8_t *up = malloc((size_t)nc * sizeof *up); /* §762: глубина подъёма пары */
  if (Esw == NULL || Edir == NULL || Eg == NULL || cen == NULL || up == NULL) exit(1);
  int64_t nv = 0, ndesc = 0, nel = 0, nlift = 0;
  for (int32_t ci = 0; ci < nc; ci++) {
    if (!(g_xcmp_aw[ci] > 0.0)) continue;
    nel++;
    int32_t szc = g_xcmp_sz[ci];
    int lvl = lev;
    while (szc > 1) {
      szc >>= 1;
      lvl--;
    }
    hz_dccell q;
    memset(&q, 0, sizeof q);
    q.lvl = (uint8_t)lvl;
    for (int a = 0; a < 3; a++)
      q.lo[a] = (uint16_t)g_xcmp_lo[ci][a];
    /* §762: спуск с подъёмом к ЗАПОЛНЕННОМУ узлу (мост А1154); прежнее
     * население (uplev = 0) отделяется — сдвиг от починки не смешивается со
     * сдвигом населения (А1195/А6). */
    int uplev = 0;
    int32_t ni = g_indfill != NULL ? node_fill_best(T, lev, &q, &uplev) : node_of_cell(T, lev, &q);
    if (ni < 0 || ni >= g_indnode_n) {
      ndesc++;
      continue;
    }
    double eg = 0.0;
    for (int a = 0; a < 3; a++)
      eg += (double)g_indnode[ni][a];
    Esw[nv] = g_xcmp_E[ci];
    Edir[nv] = g_xcmp_Edir[ci];
    Eg[nv] = eg;
    up[nv] = (int8_t)(uplev > 3 ? 3 : uplev);
    if (uplev > 0) nlift++;
    for (int a = 0; a < 3; a++)
      cen[nv][a] = fr->org[a] + ((double)g_xcmp_lo[ci][a] + 0.5 * (double)g_xcmp_sz[ci]) * fr->h;
    nv++;
  }
  printf("   §744 МОСТ: ячеек с элементами %lld, отказов спуска %lld (%.2f %%), пар для сличения "
         "%lld (из них ПОДНЯТЫХ к заполненному %lld); излучающих ячеек %d\n",
         (long long)nel, (long long)ndesc, 100.0 * (double)ndesc / (double)(nel > 0 ? nel : 1),
         (long long)nv, (long long)nlift, g_xcmp_nemit);
  xcmp_stats("СВИП/ЯДРО (все пары)", Esw, Eg, (const double (*)[3])cen, nv, fr->h);
  /* §762: прежнее население — только точный спуск */
  {
    double *E2 = malloc((size_t)(nv > 0 ? nv : 1) * sizeof *E2);
    double *G2 = malloc((size_t)(nv > 0 ? nv : 1) * sizeof *G2);
    double (*c2)[3] = malloc((size_t)(nv > 0 ? nv : 1) * sizeof *c2);
    if (E2 == NULL || G2 == NULL || c2 == NULL) exit(1);
    int64_t n2 = 0;
    for (int64_t i = 0; i < nv; i++)
      if (up[i] == 0) {
        E2[n2] = Esw[i];
        G2[n2] = Eg[i];
        for (int a = 0; a < 3; a++)
          c2[n2][a] = cen[i][a];
        n2++;
      }
    xcmp_stats("СВИП/ЯДРО (точный спуск)", E2, G2, (const double (*)[3])c2, n2, fr->h);
    free(E2);
    free(G2);
    free(c2);
  }
  xcmp_stats("НК-а ТОЖДЕСТВО", Esw, Esw, (const double (*)[3])cen, nv, fr->h);
  xcmp_stats("НК-б ПРЯМОЙ/ЯДРО", Edir, Eg, (const double (*)[3])cen, nv, fr->h);
  free(Esw);
  free(Edir);
  free(Eg);
  free(cen);
  free(up);
}

static void ind_core_build(const hz_dctree *T, const hz_htab *ht, const frame *fr, const opyr *P,
                           const arealight *AL, const hz_objmesh *m, celltris *CT, int lev,
                           double eps, int blockvis) {
  double t0 = now_s();
  g_indnode_n = T->n;
  g_indnode = calloc((size_t)T->n, sizeof *g_indnode);
  free(g_indfill);
  g_indfill = calloc((size_t)T->n, sizeof *g_indfill); /* §762 */
  if (g_indnode == NULL) exit(1);

  hz_dcslice SF;
  if (hz_slice_init(&SF, lev) != HZ_DC_OK) exit(1);
  if (hz_slice_build(&SF, T, ht, NULL, NULL) != HZ_DC_OK) exit(1);
  /* МАТЕРИАЛ ЯЧЕЙКАМ ЯДРА — ТЕМ ЖЕ ПРАВИЛОМ, ЧТО В КАДРЕ. Без него альбедо
   * приёмника берётся у материала `0`, и косвенный свет выходит НУЛЕВЫМ
   * (замерено: `Σ E_ind = 1.18` против `Σ E_dir = 2.7e5`). */
  int64_t nmatc = slice_assign_mat(&SF, CT, fr, m);
  double t_slice = now_s() - t0;

  double t1 = now_s();
  float *irrF = malloc(3 * (size_t)SF.n * sizeof *irrF);
  float *indF = calloc(3 * (size_t)SF.n, sizeof *indF);
  if (irrF == NULL || indF == NULL) exit(1);
  /* §618, `keonly`: НИКАКОГО внешнего источника. `b0` целиком из `Ke`, то есть
   * сцена светит сама — ровно §1.7. Площадка `hall_light` при этом не зовётся
   * вовсе, и её камерозависимость (§597) перестаёт мешать. */
  if (g_keonly)
    memset(irrF, 0, 3 * (size_t)SF.n * sizeof *irrF);
  else
    front_direct(&SF, fr, P, AL, irrF, 0.5, 1, NULL, m);
  /* §618: ИЗЛУЧАЮЩИЕ ПОВЕРХНОСТИ. §1.7: «источник есть полигон с заданным
   * `L_e`, отдельной машинерии нет». Прибавляется к `b0`, а не заменяет его:
   * сцена может светиться И быть освещённой. Солнце и площадка становятся
   * необязательными. */
  int64_t nemit = 0, nkemat = 0;
  if (!g_nokemit) {
    for (int32_t q = 0; q < m->nmtl; q++)
      if (m->mtl[q].ke3[0] > 0.0 || m->mtl[q].ke3[1] > 0.0 || m->mtl[q].ke3[2] > 0.0) nkemat++;
    for (int32_t i = 0; i < SF.n; i++) {
      const double *ke = m->mtl[SF.c[i].mat < m->nmtl ? SF.c[i].mat : 0].ke3;
      if (!(ke[0] > 0.0 || ke[1] > 0.0 || ke[2] > 0.0)) continue;
      for (int k = 0; k < 3; k++)
        irrF[3 * (size_t)i + (size_t)k] += (float)ke[k];
      nemit++;
    }
  }
  double t_dir = now_s() - t1;
  printf("   §618 ИЗЛУЧАЮЩИЕ ПОВЕРХНОСТИ: материалов с Ke > 0 — %lld из %d; ячеек среза со "
         "свечением %lld из %d (%.2f %%)\n",
         (long long)nkemat, m->nmtl, (long long)nemit, SF.n,
         100.0 * (double)nemit / (double)(SF.n ? SF.n : 1));

  /* РЯД НЕЙМАНА: `b0` — прямой свет, `b_{k+1} = K·b_k`, ответ `Σ b_k` (§607).
   *
   * ДЕРЕВО ИЗЛУЧАТЕЛЕЙ СТРОИТСЯ ИЗ `b_k`, А НЕ ИЗ НАКОПЛЕННОЙ СУММЫ, и это не
   * мелочь: из суммы члены ряда сложились бы ДВАЖДЫ, и «второй отскок» дал бы
   * завышение, неотличимое на глаз от настоящего переотражения.
   *
   * ЗАЧЕМ МЕРИТЬ ОТНОШЕНИЯ `Σ b_{k+1} / Σ b_k`. Оператор переноса есть СЖАТИЕ:
   * для замкнутой сцены с альбедо `ρ < 1` обязано быть `Σ b_{k+1} ≤ ρ_max·Σ b_k`.
   * Это ВНУТРЕННИЙ эталон — он не требует ни другой реализации, ни другой сцены,
   * а §599/§604 записали, что внешнего эталона отскока на Bistro нет. */
  double t_tree = 0.0, t_gath = 0.0;
  int64_t nlink = 0, nnear = 0;
  /* §609: сумма формфакторов копится ТОЛЬКО на первом отскоке — она свойство
   * ГЕОМЕТРИИ и от яркости излучателей не зависит вовсе. */
  double *ffv = calloc((size_t)SF.n, sizeof *ffv);
  if (ffv == NULL) exit(1);
  linkcache LC2;
  if (!lc_init(&LC2, g_omp1 ? 1 : omp_get_max_threads(), SF.n)) exit(1);
  float *bcur = malloc(3 * (size_t)SF.n * sizeof *bcur);
  float *bnext = calloc(3 * (size_t)SF.n, sizeof *bnext);
  if (bcur == NULL || bnext == NULL) exit(1);
  memcpy(bcur, irrF, 3 * (size_t)SF.n * sizeof *bcur);
  double tbk[HZ_BOUNCE_MAX];
  double sb[HZ_BOUNCE_MAX + 1];
  sb[0] = 0.0;
  for (int32_t i = 0; i < 3 * SF.n; i++)
    sb[0] += (double)irrF[i];
  int nb = g_hbounce < 1 ? 1 : (g_hbounce > HZ_BOUNCE_MAX ? HZ_BOUNCE_MAX : g_hbounce);
  for (int k = 0; k < nb; k++) {
    double ta1 = now_s();
    etree ETk;
    memset(&ETk, 0, sizeof ETk);
    etree_build(&ETk, &SF, fr, bcur, m, 0, SF.n, 0, lev);
    t_tree += now_s() - ta1;
    double ta2 = now_s();
    int64_t nl = 0;
    if (k == 0 || g_nolinkcache) {
      gather_run(&ETk, &SF, fr, P, m, lev, eps, blockvis, 0, bnext, &nl, NULL, NULL, &nnear,
                 k == 0 ? ffv : NULL, (k == 0 && !g_nolinkcache) ? &LC2 : NULL);
    } else {
      /* §613: ВТОРОЙ И ДАЛЬШЕ — БЕЗ ОБХОДА. Видимость уже оплачена. */
      gather_apply(&LC2, &ETk, &SF, m, ffv, bnext);
    }
    t_gath += now_s() - ta2;
    if (k == 0) nlink = nl;
    etree_free(&ETk);
    sb[k + 1] = 0.0;
    for (int32_t i = 0; i < 3 * SF.n; i++) {
      sb[k + 1] += (double)bnext[i];
      indF[i] += bnext[i];
    }
    tbk[k] = now_s() - ta1;
    float *sw = bcur;
    bcur = bnext;
    bnext = sw;
  }
  free(bcur);
  free(bnext);
  {
    int64_t nlc = 0;
    for (int t = 0; t < LC2.nth; t++)
      nlc += LC2.n[t];
    if (!g_nolinkcache)
      printf("   §613 ХРАНИМЫЕ СВЯЗИ: прошло заслон %lld из %lld геометрических (%.1f %%), "
             "память %.2f ГБ\n",
             (long long)nlc, (long long)nlink, 100.0 * (double)nlc / (double)(nlink ? nlink : 1),
             (double)nlc * sizeof(glink) / 1073741824.0);
  }
  if (g_dsa_keep)
    g_dsa_LC = LC2; /* §750: связи остаются жить для грубого звена */
  else
    lc_free(&LC2);
  /* НАИБОЛЬШЕЕ АЛЬБЕДО СЦЕНЫ — граница сжатия. Без него «больше единицы» не с
   * чем сравнивать: граница есть `ρ_max`, а не `1` (Р3 §607). */
  double rhomax = 0.0;
  for (int32_t q = 0; q < m->nmtl; q++)
    for (int k = 0; k < 3; k++) {
      double a = alb(m, (uint8_t)(q < 255 ? q : 255), k);
      if (a > rhomax) rhomax = a;
    }
  /* §609: ПОЯЧЕЕЧНЫЙ ФАЛЬСИФИКАТОР. `Σ_j F_ij ≤ 1` у замкнутой точки и СТРОГО
   * МЕНЬШЕ у открытой; Bistro — улица с небом. Всякий приёмник выше единицы есть
   * доказательство завышения, и это не сводная величина, а поячеечная. */
  {
    double *fs = malloc((size_t)SF.n * sizeof *fs);
    if (fs == NULL) exit(1);
    memcpy(fs, ffv, (size_t)SF.n * sizeof *fs);
    qsort(fs, (size_t)SF.n, sizeof *fs, cmp_d);
    int64_t nover = 0;
    for (int32_t i = 0; i < SF.n; i++)
      if (ffv[i] > 1.0) nover++;
    printf("   §609 СУММА ФОРМФАКТОРОВ Σ_j F_ij (физика: < 1 у открытой точки):\n"
           "      медиана %.4f, p90 %.4f, МАКСИМУМ %.4f; ВЫШЕ ЕДИНИЦЫ %lld из %d (%.1f %%)\n",
           fs[SF.n / 2], fs[(int32_t)((double)SF.n * 0.9)], fs[SF.n - 1], (long long)nover, SF.n,
           100.0 * (double)nover / (double)(SF.n ? SF.n : 1));
    free(fs);
  }
  printf("   §620 ЦЕНА ЗАСЛОНА: шагов иерархического марша %lld на %lld связей (%.1f на связь)\n",
         (long long)g_march_steps, (long long)nlink,
         (double)g_march_steps / (double)(nlink ? nlink : 1));
  printf("   §607 РЯД НЕЙМАНА: отскоков %d, ρ_max = %.4f\n", nb, rhomax);
  for (int k = 0; k < nb; k++)
    printf("      Σ b_%d = %.6e -> Σ b_%d = %.6e; ОТНОШЕНИЕ %.4f %s; отскок %.2f с\n", k, sb[k],
           k + 1, sb[k + 1], sb[k + 1] / (sb[k] > 0.0 ? sb[k] : 1.0),
           sb[k + 1] / (sb[k] > 0.0 ? sb[k] : 1.0) <= rhomax ? "(сжатие)" : "<- ВЫШЕ ГРАНИЦЫ",
           tbk[k]);

  /* РАСКЛАДКА ПО УЗЛАМ. Отказы СЧИТАЮТСЯ и печатаются: спуск на глубину `lvl` —
   * трудное место (А1029), и молчаливый промах выглядел бы лёгким шумом. */
  double t4 = now_s();
  int64_t nbad = 0, nput = 0;
  double sd = 0.0, si = 0.0;
  for (int32_t i = 0; i < SF.n; i++) {
    int32_t ni = node_of_cell(T, lev, &SF.c[i]);
    if (ni < 0 || !hz_dc_hasvert(T, ni)) {
      nbad++;
      continue;
    }
    for (int a = 0; a < 3; a++)
      g_indnode[ni][a] = indF[3 * (size_t)i + (size_t)a];
    g_indfill[ni] = 1; /* §762 */
    nput++;
    for (int a = 0; a < 3; a++) {
      sd += (double)irrF[3 * (size_t)i + (size_t)a];
      si += (double)indF[3 * (size_t)i + (size_t)a];
    }
  }
  ind_lift(T, 0);
  fill_lift(T, 0); /* §762: флаг заполнения — предкам */
  double t_put = now_s() - t4;

  /* А1032: ДВОЙНОЙ УЧЁТ ловится ОТНОШЕНИЕМ, а не «стало ярче». На комнате
   * эталон отношения — `0.7460` при `hgather = 0.5`. */
  printf("   §597 ЯДРО КОСВЕННОГО: приёмников %d (полная глубина), связей %lld, узлов %d, "
         "память %.1f МБ, материал назначен %lld\n"
         "      срез %.2f с + прямой %.2f с + дерево %.2f с + СБОР %.2f с + раскладка %.2f с "
         "= %.2f с; РАСКЛАДКА: узлов заполнено %lld, ОТКАЗОВ %lld\n"
         "      Σ E_dir = %.6e, Σ E_ind = %.6e, отношение = %.4f (двойной учёт дал бы вдвое)\n",
         SF.n, (long long)nlink, T->n, (double)T->n * sizeof *g_indnode / 1048576.0,
         (long long)nmatc, t_slice, t_dir, t_tree, t_gath, t_put, now_s() - t0, (long long)nput,
         (long long)nbad, sd, si, si / (sd > 0.0 ? sd : 1.0));

  free(irrF);
  free(indF);
  if (g_dsa_keep) {
    /* §750: срез и ffv переезжают в статики; освободит их некому — живут до
     * конца процесса, как g_indnode */
    g_dsa_SF = SF;
    g_dsa_ffv = ffv;
    g_dsa_ready = 1;
  } else {
    hz_slice_free(&SF);
    free(ffv); /* утечка была всегда; вскрыта анализатором при §750-ветвлении */
  }
}
static int g_dsweep = 0, g_dnmu = 2, g_dnphi = 2, g_dpass = 1, g_dirsall = 0, g_dblkopen = 0;
/* Ф9. (§536): замер границы узости доли; `dlobeflat` — негативный контроль. */
static int g_lobetest = 0, g_dlobeflat = 0;
/* Ф11. (§545): `ffull` — вернуть булев заслон (сведение); `fzero` — НК. */
static int g_ffull = 0, g_fzero = 0, g_fnotrans = 0;
/* Ф12. (§549): сторона коробки стенда; 0 — стенд не гоняется. */
static int g_diffbench = 0, g_onenb = 0, g_schar = 0, g_scharax = 0;
/* Ф14. (§557): поячеечное сличение свипа с гатером; `cmpself` — подсунуть один
 * и тот же массив дважды (проверка на ложный ноль, А1002). */
static int g_cmpcell = 0, g_cmpself = 0, g_fnoclamp = 0;
/* РЕЖИМ ТОЛЬКО РЕНДЕРА (08-12, выбор пользователя): пропустить диагностические
 * проходы, СОХРАНИВ все тайминги рабочих стадий. Прибор не должен мешать
 * измерять то, ради чего он заведён. */
static int g_render = 0;
/* §567: потолок LOD по угловому размеру ячейки, пикселей. По умолчанию равен
 * порогу невязки — иначе правка вводила бы ДВА новых числа сразу. */
static double g_lodceil = 1.0;
/* §584: пол по экранному размеру полигона, пикселей. 0 — прежнее поведение. */
static double g_lodpoly = 0.0;
/* §572: размер таблицы гаммы (`HZ_GAMN` объявлен у растеризатора; здесь стоит то
 * же число, потому что макрос определён ниже по файлу). Негативный контроль
 * `gam16` ломает её нарочно — картинка обязана пойти полосами. */
static int g_gamn = 4096;
/* НЕГАТИВНЫЙ КОНТРОЛЬ §575: вернуть постоянное альбедо вместо выборки. Картинка
 * обязана вернуться к СЕРОЙ побитово. */
static int g_texflat = 0, g_texnomip = 0;
/* НЕГАТИВНЫЙ КОНТРОЛЬ §580: отсечение по пирамиде выключено — обход и число
 * отброшенных обязаны вернуться к прежним. */
static int g_nofrustum = 0;

/* ХОДЬБА ПО СЦЕНЕ (§589). Ключ `walk`: вместо одного кадра — цикл, камера
 * правится с клавиатуры и мыши. Постройка сцены (`38` с) от этого не меняется:
 * она разовая, и цикл начинается ПОСЛЕ неё.
 *
 * ПОЧЕМУ ЭТО ПРАВКА `pfield.c`, А НЕ НОВЫЙ ИНСТРУМЕНТ. Растеризатор, прямой
 * свет, срез и таблица треугольников ячейки — статические функции ЭТОГО файла.
 * Новый инструмент значил бы копию тысячи строк, а копия — это два места, где
 * чинить один и тот же промах. */
static int g_walk = 0;

/* §804: PFM-ВЫХОД КАДРА — умолчание ВКЛ (image.h: PFM — канонический выход,
 * радианс как есть до тон-маппа), НК-ключ xnopfm возвращает мир без PFM.
 * Пишется тем же ветвлением, что и PPM (в ходьбе не пишется вовсе). */
static int g_nopfm = 0;

/* ТЕКСТУРЫ ЖИВУТ МЕЖДУ КАДРАМИ. Их загрузка — работа РАЗОВАЯ (`105` файлов,
 * `24.5` МБ), и в цикле ходьбы она платилась бы каждый кадр. Поэтому массивы
 * вынесены из кадрового контекста в файловые: `litctx` получает УКАЗАТЕЛИ на
 * них, а владение остаётся здесь. Пиксельные поля (`defuv`, `defmat`) НЕ
 * вынесены сознательно: они размером с кадр и к нему же относятся. */
static unsigned char **g_texrgb;
static int *g_texw, *g_texh;
static unsigned char ***g_texmip;
static int **g_mipw, **g_miph;
static int *g_nmip;
static int g_texloaded = 0;

/* Секундомеры сводки: заполняются рабочими стадиями по ходу. */
static double g_t0 = 0.0, g_t_frame = 0.0, g_t_fslice = 0.0, g_t_fdir = 0.0, g_t_fras = 0.0,
              g_t_bounce = 0.0;
/* КАМЕРА ДЛЯ ЛЮБОЙ СЦЕНЫ (08-12). Прежде все восемь мест читали `HZ_CFG_HALL_EYE`
 * напрямую, и никакая сцена кроме комнаты не рендерилась вовсе. `cam=` задаёт
 * шесть чисел явно; `camauto` ставит взгляд снаружи габарита — для предметов
 * (шары, дом), а для города и интерьера камеру надо задавать руками. */
static double g_eye[3] = HZ_CFG_HALL_EYE, g_at[3] = HZ_CFG_HALL_AT;
static int g_camauto = 0, g_caminside = 0;
/* А1001: прореживание излучателей гатера — ручкой, чтобы мерить ЕГО собственный
 * разброс тем же прибором. 0 — прежний автоматический выбор. */
static int g_gstride = 0;
/* Перебивка узости и зеркальной доли для СВИПА без правки сцены; -1 у `dks` —
 * «не перебивать», брать из материала. */
static double g_dns = 0.0, g_dks = -1.0;
/* Ф8'-0 (§524, А921): доля затронутых листьев печатается как ФУНКЦИЯ допуска, а
 * не при одном пороге — иначе порог был бы магическим. Пять уровней, первый
 * (0.0) есть точное неравенство, то есть отсутствие порога вовсе. */
#define HZ_HS_NTOL 5
/* Ф11. (§545): сколько значений `f` держать для распределения. Степень двойки —
 * чтобы выборка бралась маской, а не делением; это не порог, а размер буфера. */
#define HZ_FSTAT_CAP 65536
/* НЕГАТИВНЫЙ КОНТРОЛЬ §520: все нормали в ОДНУ корзину — то есть усреднение
 * через складку, как было до §519. Энергия обязана уехать вдвое. */

/* --- Ф4' (§511): СВИП ПО ДЕРЕВУ ------------------------------------------- */

/* ЗАЧЕМ. Плоский свип ходит по всем ячейкам сетки; адаптивное покрытие того же
 * уровня на комнате вдесятеро меньше (`200 593` против `2 097 152`), потому что
 * пустой блок не дробится и фронт обязан проносить его ЗА ОДИН ШАГ.
 *
 * УСТРОЙСТВО. Узлы перечисляются спуском по пирамиде занятости: пустой блок —
 * один узел, занятый дробится до уровня свипа. Индекс «ячейка → узел» строится
 * ОДИН РАЗ на кадр и делает поиск соседа O(1); он стоит O(объёма), но не на
 * образец, а на кадр, и потому дешёв.
 *
 * ПОРЯДОК ОБХОДА — рекурсивный спуск с детьми в порядке октанта: ближний к
 * источнику раньше. Условие свипа («сосед со стороны источника посчитан») тогда
 * выполняется по построению, и сортировать ничего не надо. */
typedef struct {
  int32_t lo[3], size;
  int32_t child0; /* индекс первого из 8 детей, -1 у листа */
  /* ССЫЛКИ НА СОСЕДЕЙ ПО ШЕСТИ ГРАНЯМ, предвычисленные при постройке (§513).
   * Иначе обход на каждого соседа считает координату и лезет в индекс размером
   * с сетку (8 МБ при 128³) — случайный доступ мимо кэша, и он же оказался
   * главной статьёй 43 нс на лист. Порядок: -x, +x, -y, +y, -z, +z. */
  int32_t nb[6];
} snode;

typedef struct {
  snode *nd;
  int32_t n, cap;
  int32_t *idx; /* ячейка сетки свипа -> ЛИСТ */
  int32_t nleaf;
  /* Ф5. (§514): дробление по ДОПУСКУ, а не только по занятости. Порог — тот же
   * угловой, что у среза, чтобы подробность тени и подробность поверхности были
   * согласованы. Невидимое огрубляется тем же правилом: далёкое и мелкое в
   * кадре получает крупный узел. */
  double eye[3], pxrad, thr;
  unsigned char *occl;            /* лист занят (заслон) — считается при постройке */
  int32_t *nbstart, *nblist, nnb; /* CSR: соседи листа по шести граням (§514) */
  float *nbw;                     /* доля общей площади грани */
  float *open;                    /* открытость на УЗЕЛ */
  int32_t gn;                     /* сторона сетки свипа */
} stree;

static void stree_free(stree *T) {
  free(T->nd);
  free(T->idx);
  free(T->open);
  free(T->occl);
  free(T->nbstart);
  free(T->nblist);
  free(T->nbw);
  T->occl = NULL;
  T->nbstart = NULL;
  T->nblist = NULL;
  T->nbw = NULL;
  T->nd = NULL;
  T->idx = NULL;
  T->open = NULL;
  T->n = T->cap = 0;
}

/* Восемь детей выделяются ПОДРЯД, и только потом каждый достраивается. Первая
 * редакция создавала их рекурсивно, и они ложились вразнобой, а обход считал
 * `child0 + k` — то есть ходил не туда. Замер поймал это ложным нулём: свип
 * выдал 0 освещённых ячеек за 3.8 мс. */
static int32_t stree_alloc8(stree *T) {
  if (T->n + 8 > T->cap) {
    int32_t nc = T->cap > 0 ? T->cap * 2 : 4096;
    if (nc < T->n + 8) nc = T->n + 8;
    snode *nn = realloc(T->nd, (size_t)nc * sizeof *nn);
    if (nn == NULL) exit(1);
    T->nd = nn;
    T->cap = nc;
  }
  int32_t b = T->n;
  T->n += 8;
  return b;
}

static void stree_rec(stree *T, const opyr *P, int lev, int drop, int32_t me) {
  int32_t x = T->nd[me].lo[0], y = T->nd[me].lo[1], z = T->nd[me].lo[2], size = T->nd[me].size;
  int lv = lev - drop, sh = 0;
  int32_t s2 = size;
  while (s2 > 1) {
    lv--;
    sh++;
    s2 >>= 1;
  }
  int occ = hz_occ_get(P->b[lv], hz_occ_index((int32_t)1 << lv, x >> sh, y >> sh, z >> sh)) != 0;
  T->nd[me].child0 = -1;
  if (!occ || size == 1) return;
  if (T->thr > 0.0) {
    /* Расстояние до БЛИЖАЙШЕЙ точки коробки — та же осторожная мера, что у
     * `lod_stop` (А733): у крупного узла оно разное в разных его точках. */
    double d2 = 0.0;
    for (int a = 0; a < 3; a++) {
      double lo2 = (double)T->nd[me].lo[a], hi2 = lo2 + (double)size, e2 = T->eye[a];
      double dd2 = e2 < lo2 ? lo2 - e2 : (e2 > hi2 ? e2 - hi2 : 0.0);
      d2 += dd2 * dd2;
    }
    if (d2 > 0.0) {
      double px = ((double)size / sqrt(d2)) / T->pxrad;
      if (px <= T->thr) return; /* мельче допуска в кадре — лист */
    }
  }
  int32_t h = size >> 1;
  int32_t c0 = stree_alloc8(T);
  T->nd[me].child0 = c0;
  for (int k = 0; k < 8; k++) {
    T->nd[c0 + k].lo[0] = x + ((k & 1) ? h : 0);
    T->nd[c0 + k].lo[1] = y + ((k & 2) ? h : 0);
    T->nd[c0 + k].lo[2] = z + ((k & 4) ? h : 0);
    T->nd[c0 + k].size = h;
    T->nd[c0 + k].child0 = -1;
  }
  for (int k = 0; k < 8; k++)
    stree_rec(T, P, lev, drop, c0 + k);
}

static void stree_build(stree *T, const opyr *P, int lev, int drop, int32_t gn, const double eye[3],
                        double pxrad, double thr) {
  memset(T, 0, sizeof *T);
  T->gn = gn;
  for (int k = 0; k < 3; k++)
    T->eye[k] = eye[k];
  T->pxrad = pxrad;
  T->thr = thr;
  T->nd = malloc(sizeof *T->nd);
  if (T->nd == NULL) exit(1);
  T->cap = 1;
  T->n = 1;
  T->nd[0].lo[0] = T->nd[0].lo[1] = T->nd[0].lo[2] = 0;
  T->nd[0].size = gn;
  T->nd[0].child0 = -1;
  stree_rec(T, P, lev, drop, 0);
  size_t nc = (size_t)gn * (size_t)gn * (size_t)gn;
  T->idx = malloc(nc * sizeof *T->idx);
  T->open = malloc((size_t)T->n * sizeof *T->open);
  T->occl = calloc((size_t)T->n, 1);
  if (T->idx == NULL || T->open == NULL || T->occl == NULL) exit(1);
  /* Ссылки на соседей: лист, накрывающий ЦЕНТР соответствующей грани. Для
   * крупного узла сосед через грань может быть не один; берётся тот же, что
   * брала прежняя редакция по координате, — значит замер сравнивает СКОРОСТЬ,
   * а не схему. */
  /* Индексируются ТОЛЬКО листья: внутренние узлы нужны обходу, но накрывают те
   * же ячейки, и запись их поверх листьев испортила бы индекс. */
  T->nleaf = 0;
  for (int32_t i = 0; i < T->n; i++) {
    const snode *s = &T->nd[i];
    if (s->child0 >= 0) continue;
    T->nleaf++;
    /* Занятость листа решается ЗДЕСЬ, один раз: у листа размера 1 — по пирамиде
     * мелкого уровня, у крупного она ноль по построению (он пуст, иначе бы
     * дробился). */
    T->occl[i] = (s->size == 1 &&
                  hz_occ_get(P->b[lev - drop], hz_occ_index(gn, s->lo[0], s->lo[1], s->lo[2])))
                     ? 1u
                     : 0u;
    for (int32_t z = s->lo[2]; z < s->lo[2] + s->size; z++)
      for (int32_t y = s->lo[1]; y < s->lo[1] + s->size; y++)
        for (int32_t x = s->lo[0]; x < s->lo[0] + s->size; x++)
          T->idx[hz_occ_index(gn, x, y, z)] = i;
  }
  for (int32_t i = 0; i < T->n; i++) {
    snode *s = &T->nd[i];
    if (s->child0 >= 0) continue;
    double c[3] = {(double)s->lo[0] + 0.5 * (double)s->size,
                   (double)s->lo[1] + 0.5 * (double)s->size,
                   (double)s->lo[2] + 0.5 * (double)s->size};
    for (int f = 0; f < 6; f++) {
      int ax = f / 2, sg = (f & 1) ? 1 : -1;
      double p[3] = {c[0], c[1], c[2]};
      p[ax] += (double)sg * (0.5 * (double)s->size + 0.5);
      int32_t q[3];
      int ok = 1;
      for (int a = 0; a < 3; a++) {
        double f2 = floor(p[a]);
        if (!(f2 >= 0.0) || !(f2 < (double)gn)) ok = 0;
        q[a] = ok ? (int32_t)f2 : 0;
      }
      s->nb[f] = ok ? T->idx[hz_occ_index(gn, q[0], q[1], q[2])] : -1;
    }
  }
}

typedef struct {
  unsigned char *vis; /* на ячейку: булева видимость (осевое и направленное правила) */
  float *open;        /* Ш5а2: ДОЛЯ ОТКРЫТОСТИ, переносимая с весами граней */
  /* Ф1': профиль открытости, `P²` подпроб на ячейку, `uint8` (§425/§498: свип
   * упирается в память, а квант 1/255 на порядок ниже порога приёмки 0.05). */
  unsigned char *prof;
  unsigned char *seeded; /* §507: ячейка засеяна маршем, свипом не пересчитывается */
  int32_t n;             /* сторона грубой сетки */
  int drop;              /* lev − log2(n) */
  int axis;              /* НЕГАТИВНЫЙ КОНТРОЛЬ: прежнее осевое наследование */
  int frac;              /* Ш5а2: дробная открытость вместо булевой */
  int round01;           /* НЕГАТИВНЫЙ КОНТРОЛЬ Ш5а2: округлять F до 0/1 */
  stree *tree;           /* Ф4': свип по дереву; NULL — плоская сетка */
} sweepgrid;

/* Обход дерева в порядке октанта: дети, ближние к источнику, раньше. Условие
 * свипа выполняется по построению — сортировать нечего. Правило переноса ТО ЖЕ,
 * что на плоской сетке (взвешенное среднее по трём верхним соседям), чтобы
 * сравнение шло схема в схему. */
/* СПИСКИ СОСЕДЕЙ ПО ПЛОЩАДИ ГРАНИ (Ф5', §514). У градуированного дерева сосед
 * через грань не один: крупный лист граничит с несколькими мелкими, и читать
 * только того, кто накрывает центр грани, значит терять три четверти потока.
 * Здесь для каждого листа и каждой из шести граней строится СПИСОК соседей с
 * долями общей площади. Строится раз на кадр, стоит O(суммарной площади
 * листьев), то есть по построению меньше объёма. */
static void stree_links(stree *T, int32_t gn) {
  free(T->nbstart);
  free(T->nblist);
  free(T->nbw);
  T->nbstart = malloc(((size_t)T->n * 6 + 1) * sizeof *T->nbstart);
  if (T->nbstart == NULL) exit(1);
  int32_t cap = T->n * 6, cnt = 0;
  T->nblist = malloc((size_t)cap * sizeof *T->nblist);
  T->nbw = malloc((size_t)cap * sizeof *T->nbw);
  if (T->nblist == NULL || T->nbw == NULL) exit(1);
  for (int32_t i = 0; i < T->n; i++)
    for (int f = 0; f < 6; f++) {
      T->nbstart[(size_t)i * 6 + (size_t)f] = cnt;
      const snode *s = &T->nd[i];
      if (s->child0 >= 0) continue;
      int ax = f / 2, sg = (f & 1) ? 1 : -1;
      int u = (ax + 1) % 3, v = (ax + 2) % 3;
      int32_t pos = sg > 0 ? s->lo[ax] + s->size : s->lo[ax] - 1;
      if (pos < 0 || pos >= gn) continue;
      /* Перебор ячеек грани в МЕЛКОЙ сетке; одинаковые соседи схлопываются.
       * Площадь считается в мелких ячейках, доля — от площади грани листа. */
      double tot = (double)s->size * (double)s->size;
      /* ДЛИННАЯ СЕРИЯ ОДНОГО СОСЕДА СХЛОПЫВАЕТСЯ ПО ПОИСКУ, НО НЕ ПО ПЛОЩАДИ
       * (А942). Прежняя редакция писала `if (nj == seen) continue;` и тем
       * ВЫБРАСЫВАЛА площадь всех повторов: крупный сосед получал вес одной
       * мелкой клетки вместо своей доли, и суммы весов выходили много меньше
       * единицы. Найдено первым же потребителем этих списков (Ф8', §534);
       * до него их не читал никто, поэтому ошибка и жила. Здесь повтор
       * по-прежнему не ищется заново — но вес ему добавляется. */
      int32_t seen = -1, seent = -1;
      for (int32_t a = 0; a < s->size; a++)
        for (int32_t b = 0; b < s->size; b++) {
          int32_t q[3];
          q[ax] = pos;
          q[u] = s->lo[u] + a;
          q[v] = s->lo[v] + b;
          int32_t nj = T->idx[hz_occ_index(gn, q[0], q[1], q[2])];
          if (nj == seen && seent >= 0) {
            T->nbw[seent] += 1.0f / (float)tot;
            continue;
          }
          int found = 0;
          for (int32_t t = T->nbstart[(size_t)i * 6 + (size_t)f]; t < cnt; t++)
            if (T->nblist[t] == nj) {
              T->nbw[t] += 1.0f / (float)tot;
              found = 1;
              seent = t;
              break;
            }
          if (!found) {
            if (cnt >= cap) {
              int32_t nc = cap * 2;
              int32_t *l2 = realloc(T->nblist, (size_t)nc * sizeof *l2);
              float *w2 = realloc(T->nbw, (size_t)nc * sizeof *w2);
              if (l2 == NULL || w2 == NULL) exit(1);
              T->nblist = l2;
              T->nbw = w2;
              cap = nc;
            }
            T->nblist[cnt] = nj;
            T->nbw[cnt] = 1.0f / (float)tot;
            seent = cnt;
            cnt++;
          }
          seen = nj;
        }
    }
  T->nbstart[(size_t)T->n * 6] = cnt;
  T->nnb = cnt;
}

/* --- Ф8' (§532): РАДИАНС ПО НАПРАВЛЕНИЯМ НА ТОМ ЖЕ ДЕРЕВЕ ------------------ */

/* ЧТО ЗДЕСЬ НОВОГО И ЧТО СТАРОГО. Правило переноса ТО ЖЕ, что у открытости:
 * взвешенное среднее по трём ВХОДНЫМ граням, вес оси `|ω_a|`, внутри грани —
 * доли площади из готового CSR (§514). Новая только НЕИЗВЕСТНАЯ: вместо доли
 * видимости несётся радианс, и на поверхности вместо булева заслона стоит
 * граничное условие.
 *
 * ND-СОСТОЯНИЯ В ОБЪЁМЕ НЕТ (§530, возражение пользователя «узел вырастет
 * ойойойййй»). Свип идёт ПО ОДНОМУ НАПРАВЛЕНИЮ ЗА РАЗ, поэтому узел несёт три
 * float, а не `ND × 3`: при `ND = 128` это 12 Б против 1 536 Б, то есть разница
 * между 1.8 МБ и 65 МБ на дерево. Собранная облучённость копится ПО ХОДУ
 * прохода — радианс аддитивен, это довод №3 самого §523.
 *
 * ПОВЕРХНОСТЬ — ЭТО ЛЮБОЙ ЗАНЯТЫЙ ЛИСТ, А НЕ ТОЛЬКО ТОТ, КУДА ПОПАЛА ЯЧЕЙКА
 * СРЕЗА (А935). Срез огрублён по камере, и дальняя стена представлена немногими
 * крупными ячейками; считай поверхностью только их — и свет пошёл бы сквозь
 * стену там, где срез редок. Здесь два класса врозь и оба считаются:
 *     лист с `Bs` — светит наружу своей радиосити;
 *     лист занят, но `Bs` нет — ЧЁРНАЯ стена: гасит, но не светит.
 * Второе физически честнее дыры и печатается счётчиком (А891). */
typedef struct {
  float *Ld; /* 3 на УЗЕЛ: радианс текущего направления */
  float *Bs; /* 3 на УЗЕЛ: исходящая радиосити поверхности */
  float *Bn; /* 3 на УЗЕЛ: нормаль поверхности */
  /* Ф9. (§536): УЗКАЯ ЧАСТЬ ИЗЛУЧЕНИЯ. Диффузная часть `Bs` изотропна и хранится
   * одним числом; узкая зависит от направления и потому хранится ТРЕМЯ вещами:
   * амплитудой `Bsp` (= rho_s * E_dir), направлением ПРИХОДА света `Bwi` и
   * показателем `Bns`. Это ОДНА доля на лист, а не распределение по ND, и
   * оговорка записана: первый отскок несёт узость честно, второй и дальше —
   * только диффузно (§538). */
  float *Bsp, *Bwi, *Bns;
  unsigned char *srf; /* 1 — есть Bs; 2 — занят без Bs (чёрная стена) */
  /* Лист уже посчитан в этом направлении. Ровно та же оговорка, что у
   * открытости (`open < 0` — «ещё не посчитан, не наш порядок»): у
   * градуированного дерева сосед через грань может лежать ниже по потоку, и
   * читать его нельзя. Такой сосед ИСКЛЮЧАЕТСЯ ИЗ ВЕСА, а не берётся нулём —
   * иначе схема теряла бы энергию на каждом перепаде уровня. */
  unsigned char *vis;
  int32_t *cstart, *clist;       /* CSR: лист -> ячейки среза */
  double *Eind;                  /* 3 на ячейку среза, копится по ходу */
  const double *snx, *sny, *snz; /* нормали ячеек среза */
  int32_t nsl;
  /* ВИЛКА ВМЕСТО ОДНОГО ЧИСЛА (§534). Занятый лист без ячейки среза — это
   * поверхность, которую камера не разрешила: радиосити ей взять НЕОТКУДА, и
   * обе крайности неверны. Гасить (по умолчанию) — НИЖНЯЯ граница: такой лист
   * поглощает и не светит. Пропускать (`dblkopen`) — ВЕРХНЯЯ: свет идёт сквозь
   * стену. Истина между, и замер обязан дать обе стороны, а не одну. */
  int blkopen;
  /* Ф11. (§545): доля перекрытия. `Ap` — площадь площадки в листе, `cw` —
   * сторона ячейки сетки свипа в метрах; `fone`/`fzero` — сведение к прежнему
   * правилу и негативный контроль; `fstat` — выборка `f` для распределения. */
  float *Ap;
  double cw;
  int fone, fzero, fnotrans;
  /* НЕГАТИВНЫЙ КОНТРОЛЬ §549: брать ОДНОГО входного соседа по главной оси —
   * поперечного перемешивания нет вовсе. Разброс обязан остаться нулём, а
   * коридор — пропустить почти всё. */
  int onenb;
  /* Ф13. (§553): короткие характеристики; `scharax` — НК (шаг назад вдоль оси,
   * а не вдоль ω); счётчики опор — сколько раз сработала перенормировка. */
  int schar, scharax;
  /* Ф15. (§561): бюджет излучения. `fnoclamp` — снять обрезку ТОЛЬКО в
   * излучающем члене (А1010); `emitact`/`emitcut` — излучено и обрезано;
   * `frawv` — выборка сырых `f` среди упёршихся в единицу. */
  int fnoclamp;
  double *emitact, *emitcut, *frawv;
  int64_t nfraw;
  int64_t nschar, nrenorm;
  double *fstat;
  int64_t nfstat, nfone;
} dfield;

/* --- Ф12. (§549): СТЕНД НА ДИФФУЗИЮ ПРАВИЛА ПЕРЕНОСА ---------------------- */

/* ЗАЧЕМ СТЕНД, А НЕ СЦЕНА. На комнате в одном числе смешаны четыре механизма
 * (угловая дискретизация, толщина заслона, поглощение, диффузия), и А977 отверг
 * три из них замером, оставив четвёртый непроверенным. Здесь ответ ИЗВЕСТЕН:
 * в пустоте радианс постоянен вдоль луча, пучок не расплывается, поток
 * сохраняется. Всё, что схема сделает сверх этого, и есть её диффузия.
 *
 * ПРАВИЛО БЕРЁТСЯ НАСТОЯЩЕЕ. Дерево строится синтетически (пирамида занятости —
 * все единицы, отчего листья выходят размера 1 и равномерными), а зовутся те же
 * `stree_links` и `dsweep_rec`. Копия правила разошлась бы с оригиналом, и цена
 * этому в проекте уже заплачена (§531, §541).
 *
 * КОНСЕРВАТИВНОСТЬ ДОКАЗЫВАЕТСЯ БАЛАНСОМ, А НЕ ПОТОКОМ ЧЕРЕЗ СЛОЙ (А979): у
 * косого `ω` часть пучка уходит в боковую грань раньше, чем доходит до слоя, и
 * падение потока через слой смешало бы «схема теряет» с «пучок вышел». Поэтому
 * считается вытекшее ЧЕРЕЗ ВСЕ ГРАНИ и сверяется с втекшим. */
/* --- Ф9. (§536): ШИРИНА ДОЛИ РАССЕЯНИЯ ИЗ МАТЕРИАЛА ------------------------ */

/* ЧТО ЗАДАЁТСЯ ОДНИМ ЧИСЛОМ. `f_r = ρ_d/π + ρ_s·(s+2)/(2π)·cos^s α`, где `α` —
 * угол между исходящим направлением и зеркальным отражением входящего
 * относительно нормали ЯЧЕЙКИ. При `s = 0` доля обращается в `1/π`, то есть в
 * диффузную: диапазон «диффузное … зеркало» непрерывен по построению, и
 * отдельного переключателя между режимами нет.
 *
 * КРИВЫЕ ЗЕРКАЛА ЛОЖАТСЯ САМИ: нормаль берётся у ячейки среза (DC), то есть у
 * поверхности, а не у ординаты, — кривизна несётся сеткой.
 *
 * ЗЕРКАЛЬНЫЙ ВЕКТОР НЕ СТРОИТСЯ. `cos α = ω_e·ω_d − 2(ω_d·n)(ω_e·n)` —
 * тождество, а не приближение: подставить `mirror(ω_d) = ω_d − 2(ω_d·n)n` в
 * `ω_e·mirror(ω_d)` и раскрыть. Два скалярных произведения вместо вектора.
 *
 * ГДЕ ГРАНИЦА, СКАЗАНО ДО КОДА (К3/К5/К56): отражённого направления в наборе
 * ординат нет, интерполировать между ними нельзя, поэтому доля у́же шага
 * `2π/ND` ложится на сетку — ЭНЕРГИЯ верна, ОБРАЗ разрешается лишь до шага.
 * Обе половины границы меряются порознь (А947): `A` — энергетическая, угловая
 * ошибка пика — образная. */
static double lobe(double cosa, double ns) {
  if (g_dlobeflat) return 1.0 / (2.0 * 3.14159265358979323846);
  if (!(cosa > 0.0)) return 0.0;
  if (!(ns > 0.0)) return 1.0 / 3.14159265358979323846;
  return (ns + 2.0) / (2.0 * 3.14159265358979323846) * pow(cosa, ns);
}

/* Один лист: посчитать входящий радианс, собрать его в приёмники, выставить
 * исходящий. Порядок именно такой — сбор читает ВХОДЯЩЕЕ, иначе ячейка собирала
 * бы собственное излучение. */
static void dsweep_leaf(const stree *T, dfield *D, int32_t ni, const double om[3], double wd,
                        int dirsall) {
  /* УЖЕ ПОСЧИТАН — не трогать. В рабочем пути этого не случается никогда
   * (`vis` обнуляется перед каждым направлением, а обход посещает лист ровно
   * раз), и потому поведение не меняется; проверяется побитовостью (§549 П5).
   * Нужно это СТЕНДУ (§551): им он зажигает входную грань и защищает её от
   * пересчёта, пользуясь НАСТОЯЩИМ правилом, а не копией. */
  if (D->vis[ni]) return;
  double lin[3] = {0.0, 0.0, 0.0}, wsum = 0.0;
  int amj = 0;
  for (int a = 1; a < 3; a++)
    if (fabs(om[a]) > fabs(om[amj])) amj = a;
  /* --- КОРОТКИЕ ХАРАКТЕРИСТИКИ (Ф13', §553) --------------------------------
   * Трёхгранное среднее раздаёт значение всем трём нисходящим соседям, отчего
   * носитель растёт на ячейку за шаг — конус, а не диффузия (§551). Здесь
   * вместо среднего берётся значение В ТОЧКЕ, ОТКУДА ЛУЧ ПРИШЁЛ: из центра
   * листа шагнуть назад вдоль `ω` до входной грани и интерполировать там
   * билинейно по четырём соседям. Направление тогда помнится геометрически, а
   * не «в среднем».
   * ЗАСЛОНЁННЫЕ И НЕПОСЧИТАННЫЕ ОПОРЫ исключаются из веса, а сумма
   * нормируется на принятые — то же правило, без которого §534 потерял
   * тридцатикратно (А943). Сколько раз перенормировка сработала, СЧИТАЕТСЯ
   * (А990): если часто, у результата есть названная оговорка. */
  if (D->schar) {
    int u = (amj + 1) % 3, v = (amj + 2) % 3;
    double s = (double)T->nd[ni].size;
    double c[3];
    for (int a = 0; a < 3; a++)
      c[a] = (double)T->nd[ni].lo[a] + 0.5 * s;
    /* Шаг назад до входной грани. `scharax` — НЕГАТИВНЫЙ КОНТРОЛЬ: шагать вдоль
     * главной ОСИ, то есть брать точку не на луче. Снос центроида обязан
     * вырасти линейно, и если он не вырастет — стенд к направлению слеп. */
    /* ШАГ НАЗАД — ДО ПЛОСКОСТИ ЦЕНТРОВ СОСЕДЕЙ, А НЕ ДО ГРАНИ. Первая редакция
     * брала `0.5·s/|ω|` (до грани), а опорные значения при этом лежат на
     * ПОЛКЛЕТКИ дальше — в центрах соседних ячеек. Рассогласование в половину
     * шага давало ЛИНЕЙНЫЙ СНОС пучка: замерено `3.896` ячейки к слою 32 при
     * счётном `0.5·|ω_⊥|/|ω_amj|·32 = 3.89`. Поймано ровно тем, что стенд
     * печатает снос ОТДЕЛЬНО от разброса (А992) — одной величиной это читалось
     * бы как «схема стала резче». */
    double t = s / fabs(om[amj]);
    double p[3];
    for (int a = 0; a < 3; a++)
      p[a] = D->scharax ? c[a] : c[a] - om[a] * t;
    p[amj] = om[amj] > 0.0 ? (double)T->nd[ni].lo[amj] - 0.5 : (double)T->nd[ni].lo[amj] + s + 0.5;
    double fu = p[u] - 0.5, fv = p[v] - 0.5;
    double bu = floor(fu), bv = floor(fv);
    double gu = fu - bu, gv = fv - bv;
    double acc2[3] = {0.0, 0.0, 0.0}, aw2 = 0.0;
    int dropped = 0;
    for (int du = 0; du < 2; du++)
      for (int dv = 0; dv < 2; dv++) {
        double ww = (du ? gu : 1.0 - gu) * (dv ? gv : 1.0 - gv);
        if (!(ww > 0.0)) continue;
        int32_t q[3];
        q[amj] = (int32_t)floor(p[amj]);
        q[u] = (int32_t)bu + du;
        q[v] = (int32_t)bv + dv;
        int ok = 1;
        for (int a = 0; a < 3; a++)
          if (q[a] < 0 || q[a] >= T->gn) ok = 0;
        if (!ok) {
          dropped = 1;
          continue;
        }
        int32_t nj = T->idx[hz_occ_index(T->gn, q[0], q[1], q[2])];
        if (nj < 0 || !D->vis[nj]) {
          dropped = 1;
          continue;
        }
        aw2 += ww;
        for (int k = 0; k < 3; k++)
          acc2[k] += ww * (double)D->Ld[3 * (size_t)nj + (size_t)k];
      }
    if (aw2 > 0.0)
      for (int k = 0; k < 3; k++)
        lin[k] = acc2[k] / aw2;
    if (dropped) D->nrenorm++;
    D->nschar++;
    wsum = 1.0; /* реконструкция закончена: дальше — только граничное условие */
  }
  for (int a = 0; a < 3 && !D->schar; a++) {
    if (D->onenb && a != amj) continue;
    double w = fabs(om[a]);
    if (!(w > 0.0)) continue;
    /* Входная грань по оси `a`: свет идёт в сторону `sign(om[a])`, значит
     * входит через грань с противоположной стороны. */
    int f = 2 * a + (om[a] > 0.0 ? 0 : 1);
    int32_t b0 = T->nbstart[(size_t)ni * 6 + (size_t)f];
    int32_t b1 = T->nbstart[(size_t)ni * 6 + (size_t)f + 1];
    double acc[3] = {0.0, 0.0, 0.0}, aw = 0.0;
    for (int32_t t = b0; t < b1; t++) {
      int32_t nj = T->nblist[t];
      if (nj < 0) continue;
      if (!D->vis[nj]) continue; /* ещё не посчитан — не наш порядок */
      double fw = (double)T->nbw[t];
      aw += fw;
      for (int k = 0; k < 3; k++)
        acc[k] += fw * (double)D->Ld[3 * (size_t)nj + (size_t)k];
    }
    if (!(aw > 0.0)) continue; /* грань наружу сетки: тьма, вклада нет */
    /* ДЕЛИТЬ НА `aw` ОБЯЗАТЕЛЬНО. Доли площади `nbw` суммируются в единицу лишь
     * когда ВСЕ соседи грани уже посчитаны; на перепаде уровня часть их лежит
     * ниже по потоку и исключается. Без деления недостача уходила бы прямо в
     * потерю энергии, и она уходила: первая редакция дала отношение `0.0059`
     * вместо ожидавшихся десятых долей (§534). */
    for (int k = 0; k < 3; k++)
      lin[k] += w * acc[k] / aw;
    wsum += w;
  }
  if (wsum > 0.0)
    for (int k = 0; k < 3; k++)
      lin[k] /= wsum;
  /* СБОР ПО ХОДУ: приёмник берёт входящее с косинусом, вес квадратуры — `wd`. */
  for (int32_t t = D->cstart[ni]; t < D->cstart[ni + 1]; t++) {
    int32_t i = D->clist[t];
    double cs = -(om[0] * D->snx[i] + om[1] * D->sny[i] + om[2] * D->snz[i]);
    if (!(cs > 0.0)) continue;
    for (int k = 0; k < 3; k++)
      D->Eind[3 * (size_t)i + (size_t)k] += wd * cs * lin[k];
  }
  /* ИСХОДЯЩЕЕ. Поверхность с радиосити светит в НАРУЖНУЮ полусферу и гасит
   * приходящее с изнанки; чёрная стена гасит всё. `dirsall` — негативный
   * контроль: излучать во все стороны, игнорируя нормаль. */
  if (D->srf[ni] == 1) {
    double nn[3] = {(double)D->Bn[3 * (size_t)ni], (double)D->Bn[3 * (size_t)ni + 1],
                    (double)D->Bn[3 * (size_t)ni + 2]};
    double dn = om[0] * nn[0] + om[1] * nn[1] + om[2] * nn[2];
    int out = dirsall || dn > 0.0;
    /* ДОЛЯ ПЕРЕКРЫТИЯ (Ф11', §545). Прежде лист гасил направление ЦЕЛИКОМ, хотя
     * поверхность занимает в нём не весь объём: заслон выходил толщиной в лист
     * (`0.066` м) вместо толщины поверхности, и луч, проходящий в двух
     * сантиметрах от стены, гиб. Здесь считается отношение ПАРАЛЛЕЛЬНЫХ
     * ПРОЕКЦИЙ на плоскость ⊥ω:
     *     площадка   `A_p·|ω·n|`
     *     сам лист   `h²·(|ω_x| + |ω_y| + |ω_z|)`
     * то есть доля лучей пучка, встречающих площадку. При `f = 1` правило
     * тождественно прежнему, и это проверяется приёмкой.
     * ОГОВОРКИ НАЗВАНЫ И ИЗМЕРЕНЫ ОТДЕЛЬНО: `A_p` есть сумма площадей ГРАНЕЙ
     * ячеек среза, а не наклонённой поверхности (А969, занижает `f`); проекции
     * нескольких площадок в одном листе могут перекрываться (А970, завышает);
     * нормаль усреднена и на складке не значит направления (А971). */
    double f = 1.0;
    if (!D->fone) {
      double sab = fabs(om[0]) + fabs(om[1]) + fabs(om[2]);
      double h = (double)T->nd[ni].size * D->cw;
      double cell = h * h * sab;
      f = cell > 0.0 ? (double)D->Ap[ni] * fabs(dn) / cell : 1.0;
      /* БЮДЖЕТ ИЗЛУЧЕНИЯ (Ф15', §561). Обрезка `min(1, f)` выбрасывает энергию
       * там, где ячейка среза шире листа и вся её площадь свалена в один узел.
       * Считается ДВЕ величины порознь (А1009): сколько обрезано и сколько
       * площади лежит сверх собственного сечения листа — вторая и есть перекос
       * раскладки, первая лишь его следствие. */
      double fraw = f;
      if (f > 1.0) {
        if (D->emitcut != NULL)
          *D->emitcut += wd * (fraw - 1.0) * (double)D->Bs[3 * (size_t)ni] * cell;
        if (D->frawv != NULL && D->nfraw < (int64_t)HZ_FSTAT_CAP) D->frawv[D->nfraw++] = fraw;
        /* А1010: обрезка снимается ТОЛЬКО в излучающем члене; иначе проходящий
         * член `(1−f)` стал бы отрицательным, и опыт мерил бы бессмыслицу. */
        if (!D->fnoclamp) f = 1.0;
      }
      if (D->emitact != NULL)
        *D->emitact += wd * (f > 1.0 ? f : f) * (double)D->Bs[3 * (size_t)ni] * cell;
      if (D->fzero) f = 0.0;
      if (D->fstat != NULL) {
        D->fstat[D->nfstat & (HZ_FSTAT_CAP - 1)] = f;
        D->nfstat++;
        if (f >= 1.0) D->nfone++;
      }
    }
    /* УЗКАЯ ЧАСТЬ. Зеркальный вектор не строится: `cos α = ω_e·ω_in −
     * 2(ω_in·n)(ω_e·n)` — тождество. Кривизна входит через `n`, взятую у ЯЧЕЙКИ
     * СРЕЗА, поэтому произвольная кривая зеркальная поверхность работает без
     * отдельной машинерии. */
    double lv = 0.0;
    if (out && D->Bns != NULL && (double)D->Bns[ni] > 0.0) {
      double wi[3] = {(double)D->Bwi[3 * (size_t)ni], (double)D->Bwi[3 * (size_t)ni + 1],
                      (double)D->Bwi[3 * (size_t)ni + 2]};
      double wn = wi[0] * nn[0] + wi[1] * nn[1] + wi[2] * nn[2];
      double ca = om[0] * wi[0] + om[1] * wi[1] + om[2] * wi[2] - 2.0 * wn * dn;
      lv = lobe(ca, (double)D->Bns[ni]);
    }
    /* ВЫПУКЛАЯ КОМБИНАЦИЯ, а не сумма (А973): перекрытая доля потока заменяется
     * излучением поверхности, неперекрытая проходит насквозь. Энергия не
     * рождается — значение лежит между `B_out` и `L_in`. */
    for (int k = 0; k < 3; k++) {
      double bo = out ? (double)D->Bs[3 * (size_t)ni + (size_t)k] +
                            lv * (double)D->Bsp[3 * (size_t)ni + (size_t)k]
                      : 0.0;
      /* РАЗДЕЛЕНИЕ ДВУХ ПОЛОВИН ПРАВИЛА (§547). `f` меняет СРАЗУ ДВЕ вещи:
       * ослабляет излучение (перекрыта лишь доля сечения) и пропускает остаток
       * насквозь. Чтобы сказать, какая половина что делает, `fnotrans`
       * оставляет первую и выключает вторую. */
      D->Ld[3 * (size_t)ni + (size_t)k] =
          (float)(f * bo + (D->fnotrans ? 0.0 : (f < 1.0 ? 1.0 - f : 0.0) * lin[k]));
    }
  } else if (D->srf[ni] == 2 && !D->blkopen) {
    for (int k = 0; k < 3; k++)
      D->Ld[3 * (size_t)ni + (size_t)k] = 0.0f;
  } else {
    for (int k = 0; k < 3; k++)
      D->Ld[3 * (size_t)ni + (size_t)k] = (float)lin[k];
  }
  D->vis[ni] = 1u;
}

/* Обход в октантном порядке ДЛЯ НАПРАВЛЕНИЯ: по оси `a` свет идёт в сторону
 * `sign(om[a])`, значит первым обходится ребёнок с той стороны, ОТКУДА свет
 * приходит. Условие свипа тогда выполняется по построению — как и у открытости,
 * сортировать нечего. */
static void dsweep_rec(const stree *T, dfield *D, int32_t ni, const double om[3], double wd,
                       int dirsall) {
  const snode *nd = &T->nd[ni];
  if (nd->child0 >= 0) {
    int bx = om[0] > 0.0 ? 0 : 1, by = om[1] > 0.0 ? 0 : 1, bz = om[2] > 0.0 ? 0 : 1;
    for (int i = 0; i < 8; i++) {
      int kx = (i & 1) ? 1 - bx : bx, ky = (i & 2) ? 1 - by : by, kz = (i & 4) ? 1 - bz : bz;
      dsweep_rec(T, D, nd->child0 + (kx | (ky << 1) | (kz << 2)), om, wd, dirsall);
    }
    return;
  }
  dsweep_leaf(T, D, ni, om, wd, dirsall);
}

/* ОТСТУПЛЕНИЕ ОТ §549, ЗАПИСАННОЕ, А НЕ СДЕЛАННОЕ МОЛЧА. План велел доказывать
 * консервативность балансом «втекло = вышло сбоку + дошло» (А979). Здесь стоит
 * проверка СТРОЖЕ и проще: РАВНОМЕРНЫЙ ВТОК. Если вся входная граница горит
 * `L = 1`, точное решение есть `L ≡ 1` во всём объёме, бокового вытока нет по
 * построению, и любое отклонение от единицы — ошибка схемы, без примесей.
 * Возражение А979 при этом снимается само: терять некуда. Карандаш остаётся —
 * он меряет не потери, а РАСПЛЫВАНИЕ, и для него баланс не нужен. */
static void diffbench(int nb, int corridor) {
  int lev = 0;
  while ((1 << lev) < nb)
    lev++;
  opyr P;
  memset(&P, 0, sizeof P);
  P.lev = lev;
  for (int l = 0; l <= lev; l++) {
    int32_t n2 = (int32_t)1 << l;
    size_t byc = hz_occ_bytes((size_t)n2 * (size_t)n2 * (size_t)n2);
    P.b[l] = malloc(byc);
    if (P.b[l] == NULL) exit(1);
    memset(P.b[l], 0xFF, byc); /* всё занято ⇒ дерево дробится до листа-ячейки */
  }
  stree T;
  double eye0[3] = {0.0, 0.0, 0.0};
  stree_build(&T, &P, lev, 0, nb, eye0, 1.0, 0.0);
  stree_links(&T, nb);
  dfield D;
  memset(&D, 0, sizeof D);
  D.Ld = calloc(3 * (size_t)T.n, sizeof *D.Ld);
  D.Bs = calloc(3 * (size_t)T.n, sizeof *D.Bs);
  D.Bn = calloc(3 * (size_t)T.n, sizeof *D.Bn);
  D.Bsp = calloc(3 * (size_t)T.n, sizeof *D.Bsp);
  D.Bwi = calloc(3 * (size_t)T.n, sizeof *D.Bwi);
  D.Bns = calloc((size_t)T.n, sizeof *D.Bns);
  D.Ap = calloc((size_t)T.n, sizeof *D.Ap);
  D.srf = calloc((size_t)T.n, 1);
  D.vis = calloc((size_t)T.n, 1);
  D.cstart = calloc((size_t)T.n + 1, sizeof *D.cstart);
  D.cw = 1.0;
  D.onenb = g_onenb;
  D.schar = g_schar;
  D.scharax = g_scharax;
  if (D.Ld == NULL || D.srf == NULL || D.vis == NULL || D.cstart == NULL || D.Ap == NULL) exit(1);
  /* НАПРАВЛЕНИЯ ВЫБИРАЮТСЯ ЧИСЛОМ (А981): «почти осевое» — с наибольшим
   * `max|ω_a| / Σ|ω_a|`, «диагональное» — с наименьшим. Осевых в наборе нет по
   * построению (`dirs3.h`), и на глаз их не отобрать. */
  tr3_dirs DR;
  if (tr3_dirs_product(&DR, 2, 4) != 0) exit(1);
  int dax = 0, ddg = 0;
  double bax = -1.0, bdg = 2.0;
  for (int d = 0; d < DR.n; d++) {
    double ax = fabs(DR.ox[d]), ay = fabs(DR.oy[d]), az = fabs(DR.oz[d]);
    double mx = ax > ay ? (ax > az ? ax : az) : (ay > az ? ay : az);
    double q = mx / (ax + ay + az);
    if (q > bax) {
      bax = q;
      dax = d;
    }
    if (q < bdg) {
      bdg = q;
      ddg = d;
    }
  }
  printf("== Ф12' СТЕНД НА ДИФФУЗИЮ: коробка %d³, листьев %d\n"
         "   почти осевое ω = (%.4f, %.4f, %.4f), max|ω|/Σ|ω| = %.4f\n"
         "   диагональное ω = (%.4f, %.4f, %.4f), max|ω|/Σ|ω| = %.4f\n",
         nb, T.nleaf, DR.ox[dax], DR.oy[dax], DR.oz[dax], bax, DR.ox[ddg], DR.oy[ddg], DR.oz[ddg],
         bdg);
  int dsel[2] = {dax, ddg};
  const char *dnm[2] = {"почти осевое", "диагональное"};
  for (int t = 0; t < 2; t++) {
    int d = dsel[t];
    double om[3] = {DR.ox[d], DR.oy[d], DR.oz[d]};
    /* --- A1. РАВНОМЕРНЫЙ ВТОК: точный ответ L ≡ 1 --- */
    memset(D.vis, 0, (size_t)T.n);
    memset(D.Ld, 0, 3 * (size_t)T.n * sizeof *D.Ld);
    for (int32_t i = 0; i < T.n; i++) {
      if (T.nd[i].child0 >= 0) continue;
      int on = 0;
      for (int a = 0; a < 3; a++) {
        int32_t c = T.nd[i].lo[a];
        if (om[a] > 0.0 ? (c == 0) : (c == nb - 1)) on = 1;
      }
      if (!on) continue;
      D.vis[i] = 1u;
      for (int k = 0; k < 3; k++)
        D.Ld[3 * (size_t)i + (size_t)k] = 1.0f;
    }
    dsweep_rec(&T, &D, 0, om, 1.0, 0);
    double wmin = 1e300, wmax = -1e300;
    int64_t nin = 0;
    for (int32_t i = 0; i < T.n; i++) {
      if (T.nd[i].child0 >= 0) continue;
      double v = (double)D.Ld[3 * (size_t)i];
      if (v < wmin) wmin = v;
      if (v > wmax) wmax = v;
      nin++;
    }
    printf("   A1 РАВНОМЕРНЫЙ ВТОК (%s): L ∈ [%.9f, %.9f] при точном 1, "
           "макс отклонение %.3e по %lld листьям\n",
           dnm[t], wmin, wmax,
           fabs(wmax - 1.0) > fabs(1.0 - wmin) ? fabs(wmax - 1.0) : fabs(1.0 - wmin),
           (long long)nin);
    /* --- A2. КАРАНДАШ: точный ответ — не расплывается --- */
    memset(D.vis, 0, (size_t)T.n);
    memset(D.Ld, 0, 3 * (size_t)T.n * sizeof *D.Ld);
    int32_t c0[3];
    for (int a = 0; a < 3; a++)
      c0[a] = om[a] > 0.0 ? 0 : nb - 1;
    for (int32_t i = 0; i < T.n; i++) {
      if (T.nd[i].child0 >= 0) continue;
      int on = 0;
      for (int a = 0; a < 3; a++)
        if (T.nd[i].lo[a] == c0[a]) on = 1;
      if (on) D.vis[i] = 1u; /* граница втока: вносит ноль, а не выпадает из веса */
    }
    int32_t pc[3] = {nb / 2, nb / 2, nb / 2};
    for (int a = 0; a < 3; a++)
      pc[a] = om[a] > 0.0 ? 1 : nb - 2;
    /* Карандаш ставится на оси коробки по тем осям, вдоль которых он идёт вглубь. */
    int amaj = 0;
    for (int a = 1; a < 3; a++)
      if (fabs(om[a]) > fabs(om[amaj])) amaj = a;
    for (int a = 0; a < 3; a++)
      if (a != amaj) pc[a] = nb / 2;
    int32_t pl = T.idx[hz_occ_index(nb, pc[0], pc[1], pc[2])];
    D.vis[pl] = 1u;
    for (int k = 0; k < 3; k++)
      D.Ld[3 * (size_t)pl + (size_t)k] = 1.0f;
    dsweep_rec(&T, &D, 0, om, 1.0, 0);
    /* СНОС И РАЗБРОС ПОРОЗНЬ (А992). В §551 они были смешаны: у `onenb` пучок
     * не расплывается вовсе, но целиком уезжает с луча, и одна величина
     * показывала это как «разброс». Снос — смещение ЦЕНТРОИДА от точного луча;
     * разброс — среднеквадратичное ВОКРУГ ЦЕНТРОИДА. */
    printf("   A2 КАРАНДАШ (%s): слой |  max L  | Σ L поток | СНОС | разброс | ненулевых\n",
           dnm[t]);
    for (int32_t kk = 4; kk <= nb - 4; kk *= 2) {
      int u2 = (amaj + 1) % 3, v2 = (amaj + 2) % 3;
      double eu = (double)pc[u2] + (double)kk * om[u2] / fabs(om[amaj]);
      double ev = (double)pc[v2] + (double)kk * om[v2] / fabs(om[amaj]);
      double mx2 = 0.0, s0 = 0.0, cu = 0.0, cv = 0.0, s2 = 0.0;
      int64_t nnz = 0;
      for (int pass2 = 0; pass2 < 2; pass2++) {
        for (int32_t i = 0; i < T.n; i++) {
          if (T.nd[i].child0 >= 0) continue;
          int32_t cm = T.nd[i].lo[amaj];
          int32_t off = om[amaj] > 0.0 ? cm - pc[amaj] : pc[amaj] - cm;
          if (off != kk) continue;
          double vv = (double)D.Ld[3 * (size_t)i];
          if (!(vv > 0.0)) continue;
          double du = (double)T.nd[i].lo[u2], dv = (double)T.nd[i].lo[v2];
          if (pass2 == 0) {
            if (vv > mx2) mx2 = vv;
            s0 += vv;
            cu += vv * du;
            cv += vv * dv;
            nnz++;
          } else {
            double au = du - cu, av = dv - cv;
            s2 += vv * (au * au + av * av);
          }
        }
        if (pass2 == 0 && s0 > 0.0) {
          cu /= s0;
          cv /= s0;
        }
      }
      printf("      %4d  | %.6f | %.6f | %.3f | %.3f | %lld\n", kk, mx2, s0,
             sqrt((cu - eu) * (cu - eu) + (cv - ev) * (cv - ev)), s0 > 0.0 ? sqrt(s2 / s0) : 0.0,
             (long long)nnz);
    }
  }
  if (corridor) {
    /* --- B. КОРИДОР: стены поглощают, поток обязан дойти целиком --- */
    int d = dax;
    double om[3] = {DR.ox[d], DR.oy[d], DR.oz[d]};
    int amaj = 0;
    for (int a = 1; a < 3; a++)
      if (fabs(om[a]) > fabs(om[amaj])) amaj = a;
    printf("   B КОРИДОР (%s, ось %d): ширина | доля дошедшего потока на длине %d\n", dnm[0], amaj,
           nb - 8);
    for (int W = 4; W <= 16; W *= 2) {
      memset(D.srf, 0, (size_t)T.n);
      memset(D.vis, 0, (size_t)T.n);
      memset(D.Ld, 0, 3 * (size_t)T.n * sizeof *D.Ld);
      for (int32_t i = 0; i < T.n; i++) {
        if (T.nd[i].child0 >= 0) continue;
        int wall = 0;
        for (int a = 0; a < 3; a++) {
          if (a == amaj) continue;
          if (labs((long)T.nd[i].lo[a] - (long)(nb / 2)) > W / 2) wall = 1;
        }
        if (wall) D.srf[i] = 2u; /* ЧЁРНАЯ стена: гасит, не светит */
      }
      /* РАЗМЕТКА ВТОКА. Первая редакция метила посчитанным всё, у чего ХОТЬ ОДНА
       * координата не дошла до плоскости источника, — то есть почти всю
       * коробку, и свип не делал ничего (все три ширины дали ровно ноль).
       * Верно так: посчитанными метятся ТОЛЬКО грани втока (по каждой оси — та,
       * с которой приходит свет), а источник — плоскость `lo[amaj] == c1`
       * внутри коридора. */
      int32_t c1 = om[amaj] > 0.0 ? 2 : nb - 3;
      int64_t nsrc = 0;
      for (int32_t i = 0; i < T.n; i++) {
        if (T.nd[i].child0 >= 0) continue;
        int on = 0;
        for (int a = 0; a < 3; a++) {
          int32_t c = T.nd[i].lo[a];
          if (om[a] > 0.0 ? (c == 0) : (c == nb - 1)) on = 1;
        }
        if (on) D.vis[i] = 1u;
      }
      for (int32_t i = 0; i < T.n; i++) {
        if (T.nd[i].child0 >= 0 || D.srf[i] != 0) continue;
        if (T.nd[i].lo[amaj] != c1) continue;
        D.vis[i] = 1u;
        for (int k = 0; k < 3; k++)
          D.Ld[3 * (size_t)i + (size_t)k] = 1.0f;
        nsrc++;
      }
      dsweep_rec(&T, &D, 0, om, 1.0, 0);
      int32_t c2 = om[amaj] > 0.0 ? nb - 6 : 5;
      double sin2 = (double)nsrc, sout = 0.0;
      for (int32_t i = 0; i < T.n; i++) {
        if (T.nd[i].child0 >= 0 || D.srf[i] != 0) continue;
        if (T.nd[i].lo[amaj] != c2) continue;
        sout += (double)D.Ld[3 * (size_t)i];
      }
      printf("      %4d  | %.4f  (втекло по %lld ячейкам)\n", W, sin2 > 0.0 ? sout / sin2 : 0.0,
             (long long)nsrc);
    }
  }
  free(D.Ld);
  free(D.Bs);
  free(D.Bn);
  free(D.Bsp);
  free(D.Bwi);
  free(D.Bns);
  free(D.Ap);
  free(D.srf);
  free(D.vis);
  free(D.cstart);
  tr3_dirs_free(&DR);
  stree_free(&T);
  opyr_free(&P);
}

static void tsweep_rec(stree *T, const opyr *P, int lev, int drop, const double sc[3], int32_t ni) {
  const snode *nd = &T->nd[ni];
  if (nd->child0 >= 0) {
    /* Дети в порядке октанта: ближний к источнику раньше. Условие свипа тогда
     * выполняется по построению, и сортировать нечего. */
    double cx = (double)nd->lo[0] + 0.5 * (double)nd->size;
    double cy = (double)nd->lo[1] + 0.5 * (double)nd->size;
    double cz = (double)nd->lo[2] + 0.5 * (double)nd->size;
    int bx = sc[0] > cx ? 1 : 0, by = sc[1] > cy ? 1 : 0, bz = sc[2] > cz ? 1 : 0;
    for (int i = 0; i < 8; i++) {
      int kx = (i & 1) ? 1 - bx : bx, ky = (i & 2) ? 1 - by : by, kz = (i & 4) ? 1 - bz : bz;
      tsweep_rec(T, P, lev, drop, sc, nd->child0 + (kx | (ky << 1) | (kz << 2)));
    }
    return;
  }
  if (T->open[ni] >= 0.0f) return; /* засеяно маршем либо уже посчитано */
  double c0[3] = {(double)nd->lo[0] + 0.5 * (double)nd->size,
                  (double)nd->lo[1] + 0.5 * (double)nd->size,
                  (double)nd->lo[2] + 0.5 * (double)nd->size};
  double dd[3], sabs = 0.0;
  for (int k = 0; k < 3; k++) {
    dd[k] = sc[k] - c0[k];
    sabs += fabs(dd[k]);
  }
  if (!(sabs > 0.0)) {
    T->open[ni] = 1.0f;
    return;
  }
  double acc = 0.0, wsum = 0.0;
  for (int k = 0; k < 3; k++) {
    double w = fabs(dd[k]) / sabs;
    if (!(w > 0.0)) continue;
    /* СОСЕД ПО ПРЕДВЫЧИСЛЕННОЙ ССЫЛКЕ (§513): ни координаты, ни индекса, ни
     * запроса пирамиды — одно чтение. Заслон опознаётся тем, что у занятого
     * листа открытость нулевая по построению (он ничего не пропускает). */
    int f = 2 * k + (dd[k] > 0.0 ? 1 : 0);
    int32_t nj = nd->nb[f];
    if (nj < 0) continue;
    if (T->occl[nj]) continue;        /* заслон: направление исключается */
    if (T->open[nj] < 0.0f) continue; /* ещё не посчитан — не наш порядок */
    acc += w * (double)T->open[nj];
    wsum += w;
  }
  T->open[ni] = (float)(wsum > 0.0 ? acc / wsum : 0.0);
}

static void tsweep_light(stree *T, const opyr *P, const frame *fr, int lev, int drop,
                         const double q[3]) {
  for (int32_t i = 0; i < T->n; i++)
    T->open[i] = -1.0f;
  double sc[3];
  for (int k = 0; k < 3; k++)
    sc[k] = (q[k] - fr->org[k]) / (fr->h * (double)((int32_t)1 << drop));
  /* Узел источника открыт по определению. */
  int32_t s[3];
  for (int k = 0; k < 3; k++) {
    double f = floor(sc[k]);
    if (f < 0.0) f = 0.0;
    if (f > (double)(T->gn - 1)) f = (double)(T->gn - 1);
    s[k] = (int32_t)f;
  }
  T->open[T->idx[hz_occ_index(T->gn, s[0], s[1], s[2])]] = 1.0f;
  /* §507: засев окрестности источника маршем — тот же, что у плоского свипа. */
  if (g_seed > 0) {
    double cw = fr->h * (double)((int32_t)1 << drop);
    for (int32_t dz2 = -g_seed; dz2 <= g_seed; dz2++)
      for (int32_t dy2 = -g_seed; dy2 <= g_seed; dy2++)
        for (int32_t dx2 = -g_seed; dx2 <= g_seed; dx2++) {
          int32_t cx2 = s[0] + dx2, cy2 = s[1] + dy2, cz2 = s[2] + dz2;
          if (cx2 < 0 || cy2 < 0 || cz2 < 0 || cx2 >= T->gn || cy2 >= T->gn || cz2 >= T->gn)
            continue;
          double pw[3] = {fr->org[0] + ((double)cx2 + 0.5) * cw,
                          fr->org[1] + ((double)cy2 + 0.5) * cw,
                          fr->org[2] + ((double)cz2 + 0.5) * cw};
          T->open[T->idx[hz_occ_index(T->gn, cx2, cy2, cz2)]] =
              shadowed(P, fr, pw, q, 0.5) ? 0.0f : 1.0f;
        }
  }
  tsweep_rec(T, P, lev, drop, sc, 0);
  for (int32_t i = 0; i < T->n; i++)
    if (T->open[i] < 0.0f) T->open[i] = 0.0f;
}

static void sweep_free(sweepgrid *G) {
  free(G->vis);
  free(G->open);
  free(G->prof);
  free(G->seeded);
  G->vis = NULL;
  G->open = NULL;
  G->prof = NULL;
  G->seeded = NULL;
}

/* Один образец источника: заполнить видимость на всей грубой сетке. */
static void sweep_light(sweepgrid *G, const opyr *P, const frame *fr, const double q[3]) {
  int32_t n = G->n;
  size_t nc = (size_t)n * (size_t)n * (size_t)n;
  memset(G->vis, 0, nc);
  for (size_t i = 0; i < nc; i++)
    G->open[i] = 0.0f;
  if (G->seeded != NULL) memset(G->seeded, 0, nc);
  /* Ячейка источника видима по определению — с неё начинается всякий путь. */
  int32_t s[3] = {0, 0, 0};
  for (int k = 0; k < 3; k++) {
    double f = floor((q[k] - fr->org[k]) / (fr->h * (double)((int32_t)1 << G->drop)));
    if (f < 0.0) f = 0.0;
    if (f > (double)(n - 1)) f = (double)(n - 1);
    s[k] = (int32_t)f;
  }
  G->vis[hz_occ_index(n, s[0], s[1], s[2])] = 1u;
  size_t sidx = hz_occ_index(n, s[0], s[1], s[2]);
  G->open[sidx] = 1.0f;
  /* §507: ЗАСЕВ ОКРЕСТНОСТИ ИСТОЧНИКА МАРШЕМ. Образец источника есть ТОЧКА, а
   * стартовая ячейка сетки свипа имеет сторону `2^drop` ячеек поля и «светит»
   * всем своим объёмом во все стороны — включая те, куда из точки ничего не
   * видно. Смещение от этого ПОСТОЯННО и не убывает ни с разрешением профиля,
   * ни со схемой пути, ни с измельчением сетки (§503, §506). Здесь первый слой
   * радиуса `g_seed` считается прямым маршем — ячеек там `(2R+1)³`, то есть
   * десятки, — и свип стартует с верной угловой картины. */
  if (g_seed > 0) {
    double cw = fr->h * (double)((int32_t)1 << G->drop);
    for (int32_t dz2 = -g_seed; dz2 <= g_seed; dz2++)
      for (int32_t dy2 = -g_seed; dy2 <= g_seed; dy2++)
        for (int32_t dx2 = -g_seed; dx2 <= g_seed; dx2++) {
          int32_t cx2 = s[0] + dx2, cy2 = s[1] + dy2, cz2 = s[2] + dz2;
          if (cx2 < 0 || cy2 < 0 || cz2 < 0 || cx2 >= n || cy2 >= n || cz2 >= n) continue;
          size_t k2 = hz_occ_index(n, cx2, cy2, cz2);
          double pw[3] = {fr->org[0] + ((double)cx2 + 0.5) * cw,
                          fr->org[1] + ((double)cy2 + 0.5) * cw,
                          fr->org[2] + ((double)cz2 + 0.5) * cw};
          G->open[k2] = shadowed(P, fr, pw, q, 0.5) ? 0.0f : 1.0f;
          G->vis[k2] = G->open[k2] > 0.0f ? 1u : 0u;
          G->seeded[k2] = 1u;
        }
  }
  /* Источник в координатах ГРУБОЙ сетки — к нему и строится направление. */
  double sc[3];
  for (int k = 0; k < 3; k++)
    sc[k] = (q[k] - fr->org[k]) / (fr->h * (double)((int32_t)1 << G->drop));
  /* ВОСЕМЬ ОКТАНТОВ. Внутри октанта каждая ось идёт ОТ источника, поэтому сосед
   * со стороны источника уже посчитан — это и есть условие свипа. */
  /* §508: НЕСКОЛЬКО ПРОХОДОВ БЕЗ ОБНУЛЕНИЯ. Октантный свип обходит ячейку РОВНО
   * ОДИН РАЗ, а в невыпуклой сцене свет может приходить путём, который в этом
   * порядке ещё не посчитан. Второй проход по готовому полю это показывает: если
   * значения меняются, виноват порядок; если нет — остаток есть диффузия самой
   * схемы, и патчить её нечем. */
  for (int pass = 0; pass < (g_passes > 0 ? g_passes : 1); pass++)
    for (int oct = 0; oct < 8; oct++) {
      int dx = (oct & 1) ? 1 : -1, dy = (oct & 2) ? 1 : -1, dz = (oct & 4) ? 1 : -1;
      int32_t x0 = dx > 0 ? s[0] : s[0], y0 = dy > 0 ? s[1] : s[1], z0 = dz > 0 ? s[2] : s[2];
      for (int32_t z = z0; z >= 0 && z < n; z += dz)
        for (int32_t y = y0; y >= 0 && y < n; y += dy)
          for (int32_t x = x0; x >= 0 && x < n; x += dx) {
            size_t ci = hz_occ_index(n, x, y, z);
            if (G->frac) {
              /* ЯЧЕЙКА ИСТОЧНИКА НЕ ПЕРЕСЧИТЫВАЕТСЯ. Без этой оговорки свип
               * обнулял сам источник первым же шагом: у булевых правил его
               * защищал , а дробный путь идёт мимо него. */
              if (ci == sidx) continue;
              if (G->seeded != NULL && G->seeded[ci]) continue; /* засеяно маршем (§507) */
              /* Ш5а2: ПЕРЕНОС ЧЕРЕЗ ГРАНИ С ВЕСАМИ (§423, А780). Открытость есть
               * взвешенное среднее открытостей входных соседей; веса — доли потока
               * через соответствующие грани. Ни максимума (оптимизм А778), ни
               * одного пути (пессимизм А778) — первый порядок переноса. */
              /* ЗАСЛОНОМ СЛУЖИТ ПРЕДШЕСТВЕННИК, А НЕ САМА ЯЧЕЙКА. Первая редакция
               * обнуляла открытость у занятой ячейки — и тем гасила ровно те
               * ячейки, ради которых всё считается: поверхность И ЕСТЬ занятые
               * ячейки. Булево правило этой ошибки не имело, потому что проверяло
               * занятость СОСЕДА. */
              double c0[3] = {(double)x + 0.5, (double)y + 0.5, (double)z + 0.5};
              double dd[3], sabs = 0.0;
              for (int k = 0; k < 3; k++) {
                dd[k] = sc[k] - c0[k];
                sabs += fabs(dd[k]);
              }
              if (!(sabs > 0.0)) {
                G->open[ci] = 1.0f;
                continue;
              }
              if (g_raysweep) {
                /* Ф2' (§505): ПЕРЕНОС ВДОЛЬ ПРЯМОЙ. Шаг назад делается только по
                 * ГЛАВНОЙ оси, а поперечное смещение точки входа берётся из самого
                 * направления — то есть путь есть прямая, а не лесенка. Четыре
                 * вкладчика лежат в ОДНОЙ плоскости, поперёк одного луча.
                 * ЗАНЯТЫЙ ВКЛАДЧИК ВХОДИТ НУЛЁМ: здесь он значит «эта доля сечения
                 * пучка заслонена», а не «этот путь закрыт», — и потому обнулять
                 * его физически верно (разбор в §505). */
                int km = 0;
                for (int k = 1; k < 3; k++)
                  if (fabs(dd[k]) > fabs(dd[km])) km = k;
                int u2 = (km + 1) % 3, v2 = (km + 2) % 3;
                double inv = fabs(dd[km]) > 0.0 ? 1.0 / fabs(dd[km]) : 0.0;
                double fu = (double)(km == 0 ? x : (km == 1 ? y : z));
                (void)fu;
                double pu =
                    (double)(u2 == 0 ? x : (u2 == 1 ? y : z)) + (g_noshift ? 0.0 : dd[u2] * inv);
                double pv =
                    (double)(v2 == 0 ? x : (v2 == 1 ? y : z)) + (g_noshift ? 0.0 : dd[v2] * inv);
                int32_t bu = (int32_t)floor(pu), bv = (int32_t)floor(pv);
                double tu = pu - (double)bu, tv = pv - (double)bv;
                int32_t base[3] = {x, y, z};
                base[km] += (dd[km] > 0.0 ? 1 : -1);
                double accr = 0.0;
                for (int au = 0; au < 2; au++)
                  for (int av = 0; av < 2; av++) {
                    int32_t q2[3] = {base[0], base[1], base[2]};
                    q2[u2] = bu + au;
                    q2[v2] = bv + av;
                    double w2 = (au ? tu : 1.0 - tu) * (av ? tv : 1.0 - tv);
                    if (!(w2 > 0.0)) continue;
                    if (q2[0] < 0 || q2[1] < 0 || q2[2] < 0 || q2[0] >= n || q2[1] >= n ||
                        q2[2] >= n)
                      continue;
                    size_t qi = hz_occ_index(n, q2[0], q2[1], q2[2]);
                    if (hz_occ_get(P->b[P->lev - G->drop], qi)) continue; /* заслонено: ноль */
                    accr += w2 * (double)G->open[qi];
                  }
                if (G->round01) accr = accr >= 0.5 ? 1.0 : 0.0;
                G->open[ci] = (float)accr;
                continue;
              }
              double acc = 0.0, wsum = 0.0;
              int32_t pp[3] = {x, y, z};
              for (int k = 0; k < 3; k++) {
                double w = fabs(dd[k]) / sabs;
                if (!(w > 0.0)) continue;
                int32_t save = pp[k];
                pp[k] = save + (dd[k] > 0.0 ? 1 : -1);
                if (pp[k] >= 0 && pp[k] < n) {
                  size_t pi = hz_occ_index(n, pp[0], pp[1], pp[2]);
                  /* ЗАКРЫТОЕ НАПРАВЛЕНИЕ ИСКЛЮЧАЕТСЯ ИЗ СРЕДНЕГО, А НЕ ВХОДИТ В
                   * НЕГО НУЛЁМ. Иначе у стены одно направление из трёх всегда
                   * закрыто, среднее падает на каждом шаге, и открытость ТАЕТ
                   * вдоль стены — это ложное затухание, а не тень. Замерено на
                   * первой редакции: средняя открытость 0.0605 против 0.3995 у
                   * эталона. */
                  if (hz_occ_get(P->b[P->lev - G->drop], pi)) continue;
                  acc += w * (double)G->open[pi];
                  wsum += w;
                }
                pp[k] = save;
              }
              double f = wsum > 0.0 ? acc / wsum : 0.0;
#if HZ_SWEEP_P > 1
              /* Ф1' (§502): ПРОФИЛЬ НА ГРАНИ ВМЕСТО ОДНОГО ЧИСЛА. Ячейка несёт
               * `P×P` подпроб открытости в плоскости, поперечной ГЛАВНОЙ оси
               * направления на источник. Подпроба берётся у верхнего по потоку
               * соседа ПО ГЛАВНОЙ ОСИ со сдвигом `P·d_u/|d_k|` — это и есть форма
               * тени, переносимая вдоль луча, а не размазанная средним.
               *
               * ПРИБЛИЖЕНИЕ, НАЗВАННОЕ ЗДЕСЬ: два ПОБОЧНЫХ соседа отдают своё
               * СРЕДНЕЕ, а не профиль. Их профили лежат в других плоскостях, и
               * честный перенос потребовал бы поворота выборки. Главная ось несёт
               * границу тени (по ней идёт основной поток), побочные подмешивают
               * фон — то есть приближение бьёт по фону, а не по границе. Если
               * замер покажет, что этого мало, следующий ход — общий профиль в
               * плоскости, поперечной СРЕДНЕМУ направлению, а не по осям. */
              int kmax = 0;
              for (int k = 1; k < 3; k++)
                if (fabs(dd[k]) > fabs(dd[kmax])) kmax = k;
              int uu = (kmax + 1) % 3, vv = (kmax + 2) % 3;
              unsigned char *pc = G->prof + ci * (size_t)(HZ_SWEEP_P * HZ_SWEEP_P);
              int32_t up[3] = {x, y, z};
              up[kmax] += (dd[kmax] > 0.0 ? 1 : -1);
              int haveup = 0;
              const unsigned char *pu = NULL;
              if (up[0] >= 0 && up[0] < n && up[1] >= 0 && up[1] < n && up[2] >= 0 && up[2] < n) {
                size_t ui = hz_occ_index(n, up[0], up[1], up[2]);
                if (!hz_occ_get(P->b[P->lev - G->drop], ui)) {
                  pu = G->prof + ui * (size_t)(HZ_SWEEP_P * HZ_SWEEP_P);
                  haveup = 1;
                }
              }
              /* Сдвиг в подпробах: пересечение ячейки вдоль `d` смещает луч на
               * `d_u/|d_k|` ячейки поперёк, то есть на `P·d_u/|d_k|` подпроб. */
              double shu = (double)HZ_SWEEP_P * dd[uu] / fabs(dd[kmax]);
              double shv = (double)HZ_SWEEP_P * dd[vv] / fabs(dd[kmax]);
              double bg = f * 255.0; /* фон — тот же взвешенный ответ, что и раньше */
              double psum = 0.0;
              for (int a = 0; a < HZ_SWEEP_P; a++)
                for (int b = 0; b < HZ_SWEEP_P; b++) {
                  double val = bg;
                  if (haveup) {
                    int sa = (int)lround((double)a + (g_noshift ? 0.0 : shu));
                    int sb = (int)lround((double)b + (g_noshift ? 0.0 : shv));
                    if (sa >= 0 && sa < HZ_SWEEP_P && sb >= 0 && sb < HZ_SWEEP_P)
                      val = (double)pu[sa * HZ_SWEEP_P + sb];
                    /* Вышли за грань — луч пришёл из СОСЕДНЕЙ ячейки того же
                     * уровня; её профиля здесь нет, и берётся взвешенный фон.
                     * Это та же диффузия, но только на краю профиля, а не везде. */
                  }
                  pc[a * HZ_SWEEP_P + b] =
                      (unsigned char)(val < 0.0     ? 0
                                      : val > 255.0 ? 255
                                                    : (unsigned char)lround(val));
                  psum += (double)pc[a * HZ_SWEEP_P + b];
                }
              f = psum / (255.0 * (double)(HZ_SWEEP_P * HZ_SWEEP_P));
#endif
              /* НЕГАТИВНЫЙ КОНТРОЛЬ (§423): округление до 0/1 обязано вернуть
               * смещения булевых правил. */
              if (G->round01) f = f >= 0.5 ? 1.0 : 0.0;
              G->open[ci] = (float)f;
              continue;
            }
            if (G->vis[ci]) continue;
            int v = 0;
            if (G->axis) {

              /* НЕГАТИВНЫЙ КОНТРОЛЬ (§419): прежнее правило — максимум по трём
               * ОСЕВЫМ соседям. Путь ступенчатый, свет заворачивает за угол, и
               * расхождение с эталоном обязано вернуться к `18 %`. */
              int32_t px = x - dx, py = y - dy, pz = z - dz;
              if (px >= 0 && px < n) {
                size_t pi = hz_occ_index(n, px, y, z);
                if (G->vis[pi] && !hz_occ_get(P->b[P->lev - G->drop], pi)) v = 1;
              }
              if (!v && py >= 0 && py < n) {
                size_t pi = hz_occ_index(n, x, py, z);
                if (G->vis[pi] && !hz_occ_get(P->b[P->lev - G->drop], pi)) v = 1;
              }
              if (!v && pz >= 0 && pz < n) {
                size_t pi = hz_occ_index(n, x, y, pz);
                if (G->vis[pi] && !hz_occ_get(P->b[P->lev - G->drop], pi)) v = 1;
              }
            } else {
              /* Ш5а1: ПРЕДШЕСТВЕННИК ВДОЛЬ НАСТОЯЩЕГО НАПРАВЛЕНИЯ НА ИСТОЧНИК.
               * Берётся ячейка, содержащая точку `центр − шаг·d`, где `d` —
               * единичное направление на источник. Предшественник ОДИН, а не
               * максимум по трём, — и потому свет не может «свернуть за угол»:
               * путь идёт по прямой, а не ступенькой.
               * Шаг — сторона ячейки: меньший шаг дал бы ту же ячейку, больший
               * перепрыгнул бы заслон. */
              double c0[3] = {(double)x + 0.5, (double)y + 0.5, (double)z + 0.5};
              double d[3], dl = 0.0;
              for (int k = 0; k < 3; k++) {
                d[k] = sc[k] - c0[k];
                dl += d[k] * d[k];
              }
              dl = sqrt(dl);
              if (!(dl > 0.0)) {
                G->vis[ci] = 1u;
                continue;
              }
              int32_t qc[3];
              int ok = 1;
              for (int k = 0; k < 3; k++) {
                double w = c0[k] + d[k] / dl;
                int32_t iw = (int32_t)floor(w);
                if (iw < 0 || iw >= n) ok = 0;
                qc[k] = ok ? iw : 0;
              }
              if (ok) {
                size_t pi = hz_occ_index(n, qc[0], qc[1], qc[2]);
                if (G->vis[pi] && !hz_occ_get(P->b[P->lev - G->drop], pi)) v = 1;
              }
            }
            if (v) G->vis[ci] = 1u;
          }
    }
}

static double sweep_vis(const sweepgrid *G, const frame *fr, const double p[3]) {
  int32_t c[3];
  for (int k = 0; k < 3; k++) {
    double f = floor((p[k] - fr->org[k]) / (fr->h * (double)((int32_t)1 << G->drop)));
    if (!(f >= 0.0) || !(f < (double)G->n)) return 0.0;
    c[k] = (int32_t)f;
  }
  size_t ci = hz_occ_index(G->n, c[0], c[1], c[2]);
  if (G->tree != NULL) return (double)G->tree->open[G->tree->idx[ci]];
  return G->frac ? (double)G->open[ci] : (G->vis[ci] ? 1.0 : 0.0);
}

/* Прямая облучённость СВИПОМ. Отличие от `front_direct` только в том, откуда
 * берётся затенение; геометрия (косинус, `1/r²`, образцы площадки) та же — иначе
 * сверка мерила бы разницу формул, а не разницу механизмов. */
static void front_sweep(const hz_dcslice *S, const frame *fr, const opyr *P, const arealight *L,
                        float *irr, const hz_objmesh *A, double *t_sweep, double *t_gather,
                        int axismode, int fracmode, int round01) {
  double su[HZ_LIGHT_SAMPLES], sv[HZ_LIGHT_SAMPLES];
  for (int a = 0; a < HZ_LIGHT_NS; a++)
    for (int b = 0; b < HZ_LIGHT_NS; b++) {
      su[a * HZ_LIGHT_NS + b] = lsamp(a);
      sv[a * HZ_LIGHT_NS + b] = lsamp(b);
    }
  sweepgrid G;
  memset(&G, 0, sizeof G);
  G.drop = HZ_SWEEP_DROP;
  G.axis = axismode;
  G.frac = fracmode;
  G.round01 = round01;
  G.n = (int32_t)1 << (fr->lev - G.drop);
  size_t gcells = (size_t)G.n * (size_t)G.n * (size_t)G.n;
  G.vis = malloc(gcells);
  G.open = malloc(gcells * sizeof *G.open);
  G.prof = HZ_SWEEP_P > 1 ? malloc(gcells * (size_t)(HZ_SWEEP_P * HZ_SWEEP_P)) : NULL;
  G.seeded = calloc(gcells, 1);
  if (G.vis == NULL || G.open == NULL) exit(1);
  stree TR;
  if (g_treesweep) {
    double tt = now_s();
    double eyeg[3];
    for (int k = 0; k < 3; k++)
      eyeg[k] = g_sweepeye[k] / (double)((int32_t)1 << G.drop);
    stree_build(&TR, P, fr->lev, G.drop, G.n, eyeg, g_sweeppx, g_sweepthr);
    G.tree = &TR;
    stree_links(&TR, G.n);
    printf("   Ф4' ДЕРЕВО СВИПА: узлов %d против %lld ячеек плоской сетки (в %.1f раза меньше), "
           "перечень и индекс за %.1f мс\n",
           TR.nleaf, (long long)gcells, (double)gcells / (double)(TR.nleaf ? TR.nleaf : 1),
           (now_s() - tt) * 1e3);
  }
  for (int32_t i = 0; i < 3 * S->n; i++)
    irr[i] = 0.0f;
  *t_sweep = 0.0;
  *t_gather = 0.0;
  for (int sm = 0; sm < HZ_LIGHT_SAMPLES; sm++) {
    double q[3];
    for (int k = 0; k < 3; k++)
      q[k] = L->c[k] + L->u[k] * su[sm] + L->v[k] * sv[sm];
    double ta = now_s();
    if (g_treesweep)
      tsweep_light(&TR, P, fr, fr->lev, G.drop, q);
    else
      sweep_light(&G, P, fr, q);
    *t_sweep += now_s() - ta;
    ta = now_s();
    for (int32_t i = 0; i < S->n; i++) {
      double p[3], n[3];
      hz_slice_vertex(S, i, p);
      for (int k = 0; k < 3; k++)
        p[k] = fr->org[k] + p[k] * fr->h;
      hz_slice_normal(S, i, n);
      double w[3], r2 = 0.0;
      for (int k = 0; k < 3; k++) {
        w[k] = q[k] - p[k];
        r2 += w[k] * w[k];
      }
      if (!(r2 > 0.0)) continue;
      double r = sqrt(r2);
      double cosr = (w[0] * n[0] + w[1] * n[1] + w[2] * n[2]) / r;
      if (!(cosr > 0.0)) continue;
      double vis = sweep_vis(&G, fr, p);
      if (!(vis > 0.0)) continue;
      double g = vis * alight_g(L, w, r2, cosr);
      for (int k = 0; k < 3; k++)
        irr[3 * (size_t)i + (size_t)k] += (float)(L->rgb[k] * g * alb(A, S->c[i].mat, k));
    }
    *t_gather += now_s() - ta;
  }
  if (g_treesweep) stree_free(&TR);
  sweep_free(&G);
}
/* ИСПОЛНИТЕЛЬ ЗАПРЕТА §261, КОТОРОГО НЕ БЫЛО ДЕСЯТЬ ДНЕЙ (§212 записал этот долг
 * 08-04 и он не был заплачен; А958). Проверяется ЗАКОНОМ, а не допуском, и
 * гоняется ДО всякого счёта — как `hz_pff_selftest` в линии переноса.
 *
 * ТРИ ЧАСТИ, И КАЖДАЯ ЛОВИТ СВОЁ (А959, А962):
 *   ДАЛЬНЕЕ ПОЛЕ ПОД УГЛОМ — сходится к `L·A·cos θ_s·cos θ_r / r²`. Под углом, а
 *       НЕ на оси: на оси `cos θ_s = 1`, и проверка была бы слепа ровно к тому
 *       множителю, ради которого ставится.
 *   БЛИЖНЕЕ ПОЛЕ — `E ≤ π·L` при ЛЮБОМ `r`, включая `r → 0`. Это физический
 *       потолок: облучённость от источника с радиансом `L` не может превысить
 *       `πL` ни при какой геометрии. Порога нет — есть закон, и точечная
 *       формула его нарушает на порядки.
 *   ЛИНЕЙНОСТЬ ПО ПЛОЩАДИ — вдвое большая площадка при том же радиансе и том же
 *       дальнем расстоянии даёт ровно вдвое большую облучённость. Ноль при
 *       `A = 0` этого НЕ ловит: `a²` вместо `a` дал бы тот же ноль.
 *
 * Провал — `exit`, а не предупреждение: запрет, который можно не заметить, и
 * есть та самая «текстовая» форма, из-за которой ошибка прожила десять дней. */
static double alight_sum(const arealight *L, const double p[3], const double nr[3]) {
  double acc = 0.0;
  for (int a = 0; a < HZ_LIGHT_NS; a++)
    for (int b = 0; b < HZ_LIGHT_NS; b++) {
      /* Явный ноль — не перестраховка: gcc-analyzer не доказывает, что цикл ниже
       * заполняет все три компоненты, и без инициализации гейт краснеет. */
      double q[3] = {0.0, 0.0, 0.0}, w[3] = {0.0, 0.0, 0.0}, r2 = 0.0;
      for (int k = 0; k < 3; k++) {
        q[k] = L->c[k] + L->u[k] * lsamp(a) + L->v[k] * lsamp(b);
        w[k] = q[k] - p[k];
        r2 += w[k] * w[k];
      }
      if (!(r2 > 0.0)) continue;
      double r = sqrt(r2);
      double cosr = (w[0] * nr[0] + w[1] * nr[1] + w[2] * nr[2]) / r;
      if (!(cosr > 0.0)) continue;
      acc += alight_g(L, w, r2, cosr);
    }
  return acc;
}

/* САМОПРОВЕРКА БЮДЖЕТА ИЗЛУЧЕНИЯ (§561, требование А1008). Прежде чем читать
 * невязку на сцене, надо убедиться, что `Φ_факт` и `Φ_аналит` — одна и та же
 * величина, а не спор о `π`. Ответ известен точно: площадка `A` с радиосити `B`
 * излучает в полусферу поток `A·B`.
 * `Φ_факт` в коде складывается как `Σ_d w_d · f · Bs · h²Σ|ω_a|`, где
 * `Bs = B/π` — радианс, `f = A|ω·n|/(h²Σ|ω_a|)`. Подставив, получаем
 * `Σ_d w_d · A|ω·n| · B/π` — а это в точности `A·B/π · ∫|cos| dω` по полусфере,
 * то есть `A·B/π · π = A·B`. Проверяется квадратурой, а не рассуждением. */
static void emitbudget_selftest(void) {
  const double PI = 3.14159265358979323846;
  tr3_dirs DR;
  if (tr3_dirs_product(&DR, 4, 4) != 0) exit(1);
  double A = 0.37, B = 2.4, h = 1.0; /* площадка меньше сечения листа: обрезки нет */
  double n[3] = {0.0, 1.0, 0.0};
  double acc = 0.0;
  int64_t ncut = 0;
  for (int d = 0; d < DR.n; d++) {
    double om[3] = {DR.ox[d], DR.oy[d], DR.oz[d]};
    double dn = om[0] * n[0] + om[1] * n[1] + om[2] * n[2];
    if (!(dn > 0.0)) continue; /* наружная полусфера */
    double sab = fabs(om[0]) + fabs(om[1]) + fabs(om[2]);
    double cell = h * h * sab;
    double f = A * fabs(dn) / cell;
    if (f > 1.0) {
      f = 1.0;
      ncut++;
    }
    acc += DR.w[d] * f * (B / PI) * cell;
  }
  printf("== БЮДЖЕТ ИЗЛУЧЕНИЯ, САМОПРОВЕРКА (§561): Φ_факт = %.6f против A·B = %.6f, "
         "расхождение %.3f %% (порог 2 %%); обрезано направлений %lld\n",
         acc, A * B, 100.0 * fabs(acc - A * B) / (A * B), (long long)ncut);
  if (!(fabs(acc - A * B) / (A * B) < 0.02)) {
    printf("   ПРОВАЛ: Φ_факт и Φ_аналит — разные величины, невязку на сцене читать нельзя.\n");
    exit(3);
  }
  tr3_dirs_free(&DR);
}

static void alight_selftest(void) {
  const double PI = 3.14159265358979323846;
  arealight L;
  memset(&L, 0, sizeof L);
  L.u[0] = 1.0;
  L.v[2] = 1.0; /* нормаль u×v = (0,−1,0), площадка смотрит ВНИЗ */
  double A = 4.0 * 1.0 * 1.0;
  /* 1. ДАЛЬНЕЕ ПОЛЕ ПОД УГЛОМ 60°: cos θ_s = 0.5, cos θ_r = 1. */
  double rr = 10.0 * sqrt(A);
  double d[3] = {0.8660254037844386, -0.5, 0.0};
  double p[3] = {L.c[0] + rr * d[0], L.c[1] + rr * d[1], L.c[2] + rr * d[2]};
  double nr[3] = {-d[0], -d[1], -d[2]};
  double got = alight_sum(&L, p, nr);
  double want = A * 0.5 / (rr * rr);
  double dev = 100.0 * fabs(got - want) / want;
  /* 2. БЛИЖНЕЕ ПОЛЕ. ПРИЁМНИК СТАВИТСЯ ПОД ПРОБУ, А НЕ ПО ОСИ (А965), и это не
   * мелочь: по оси ближайшая проба всё равно отстоит вбок на полшага сетки
   * проб, расстояние до неё снизу не убывает, косинус приёмника стремится к
   * нулю — и особенность точечной формулы НЕ ПРОЯВЛЯЕТСЯ. Первая редакция
   * проверки так и промолчала на негативном контроле. Под пробой же `r → 0`
   * по-настоящему: верная формула упирается в `π·a/a = π`, то есть ровно в
   * потолок закона, а точечная уходит в бесконечность как `1/r²`. */
  double worst = 0.0, rworst = 0.0;
  double ps[3];
  for (int k = 0; k < 3; k++)
    ps[k] = L.c[k] + L.u[k] * lsamp(0) + L.v[k] * lsamp(0);
  for (double t = 1.0; t >= 1e-4; t *= 0.5) {
    double rn = t * sqrt(A);
    double p2[3] = {ps[0], ps[1] - rn, ps[2]};
    double nr2[3] = {0.0, 1.0, 0.0};
    double e = alight_sum(&L, p2, nr2) / PI;
    if (e > worst) {
      worst = e;
      rworst = rn;
    }
  }
  /* 3. ЛИНЕЙНОСТЬ ПО ПЛОЩАДИ: вдвое шире — вдвое ярче. */
  arealight L2 = L;
  L2.u[0] = 2.0;
  double got2 = alight_sum(&L2, p, nr);
  double lin = got2 / (got > 0.0 ? got : 1.0);
  printf("== ИСТОЧНИК, САМОПРОВЕРКА ЗАКОНОМ (§541, исполнитель запрета §261)\n"
         "   ДАЛЬНЕЕ ПОЛЕ под 60°: %.6e против A·cos/r² = %.6e, расхождение %.3f %% "
         "(порог 1 %%)\n"
         "   БЛИЖНЕЕ ПОЛЕ: max E/(π·L) = %.4e при r = %.4f м (закон: ≤ 1)\n"
         "   ЛИНЕЙНОСТЬ ПО ПЛОЩАДИ: вдвое шире дало ×%.4f (закон: ×2)\n",
         got, want, dev, worst, rworst, lin);
  if (g_ptlight) {
    printf("   [НК ptlight: прежняя точечная формула — закон обязан быть НАРУШЕН]\n");
    return;
  }
  if (!(dev < 1.0)) {
    printf("   ПРОВАЛ: дальнее поле не сходится к площадке.\n");
    exit(3);
  }
  if (!(worst <= 1.0)) {
    printf("   ПРОВАЛ: облучённость превысила π·L — источник ТОЧЕЧНЫЙ (§261).\n");
    exit(3);
  }
  if (!(fabs(lin - 2.0) < 0.02)) {
    printf("   ПРОВАЛ: площадь входит не первой степенью.\n");
    exit(3);
  }
}

/* --- 3ж. РАСТЕРИЗАЦИЯ СО СВЕТОМ (Ш5) --------------------------------------- */

/* ЦВЕТ ИНТЕРПОЛИРУЕТСЯ ПО МНОГОУГОЛЬНИКУ, а не берётся плоским на треугольник —
 * ровно как записано в А663: угол многоугольника ЕСТЬ дуальная вершина ячейки,
 * значит облучённость в этом углу — облучённость ТОЙ ЯЧЕЙКИ, и никаких
 * вершинных величин заводить не надо.
 *
 * БЕЛАЯ ТОЧКА — не подобранная константа, а перцентиль по ЯЧЕЙКАМ СРЕЗА (то же
 * правило, что в `hz_ppm_write`: одиночный яркий блик не должен утопить кадр).
 * Гамма `1/2.2`. Ложноцветной палитры здесь нет: она годится полю интенсивности,
 * а на геометрии делает картинку нечитаемой. */
struct littri {
  double p[3][3], col[3][3], uv[3][2];
  int mat;
  /* Габарит по строкам, посчитанный ОДИН раз при сборе: без него каждая полоса
   * перепроецировала бы все треугольники заново, и деление на потоки не давало
   * ничего (замерено: 218 -> 234 мс, то есть хуже). */
  int iy0, iy1;
};

typedef struct {
  const hz_dcslice *S;
  const uint64_t *key;
  const int32_t *ord;
  const float *irr;
  const frame *fr;
  const tr3_camera *cam;
  double *z;
  int w, h;
  unsigned char *rgb;
  double white;
  int nocull;           /* НЕГАТИВНЫЙ КОНТРОЛЬ: отсечение выключено */
  int64_t nseen, ncull; /* сколько многоугольников пришло и сколько отброшено */
  uint64_t polysum;     /* порядковая сумма потока многоугольников (А1013) */
  /* Ш8 (§575): текстуры. Таблица по МАТЕРИАЛУ; `uv` и материал кладутся в
   * отложенный буфер вместе с цветом и выбираются ОДИН раз на видимый пиксель.
   * Т1: альбедо применяется при ЧТЕНИИ поля, поэтому в перенос текстура не
   * входит вовсе — и цветного непрямого света от неё не будет (§575). */
  /* Ш8б (08-12): МИП-ПИРАМИДА на материал. Точечная выборка рябила на листве и
   * дальних поверхностях; уровень выбирается по следу пикселя в текселях, а
   * внутри уровня и между уровнями идёт линейная интерполяция (трилинейная).
   * Память растёт на треть (сумма 1+1/4+1/16+… = 4/3) — это цена, названная
   * до кода. */
  unsigned char **texrgb; /* [материал][уровень] */
  int *texw, *texh;       /* размеры УРОВНЯ 0 */
  unsigned char ***texmip;
  int **mipw, **miph;
  int *nmip;
  const float *uvs;
  float *defuv;
  unsigned char *defmat;
  int64_t ntexpx;
  /* Замечание пользователя 08-12: координата берётся С ПОВЕРХНОСТИ, а не с
   * нашего многоугольника. Для этого нужны сама сетка и указатель ячейка ->
   * треугольники; `nsurfuv` считает, скольким пикселям это удалось. */
  const hz_objmesh *mesh;
  celltris *ct;
  int nmtl;
  int64_t nsurfuv;
  double pxrad; /* радиан на пиксель — нужен для выбора уровня пирамиды */
  int nomip;    /* НК: точечная выборка нулевого уровня */
  /* Р3 (§572): буфер отложенного затенения — радианс на пиксель, три канала.
   * `nfrag` — сколько фрагментов прошло z; отношение к числу закрытых пикселей
   * есть ГЛУБИНА ПЕРЕКРЫТИЯ, и от неё прямо зависит выигрыш Р3. */
  float *defcol;
  int64_t nfrag;
  /* Р3 (§581): СБОР ТРЕУГОЛЬНИКОВ, потом отрисовка по полосам. Обход остаётся
   * однопоточным (у него общий выход), а рисование делится: каждый поток берёт
   * свою полосу строк и трогает только свои пиксели. */
  struct littri *tris;
  int64_t ntris, captris;
  /* Гистограмма площади полигона на экране: корзина `k` — площадь `4^k…4^(k+1)`
   * пикселей. Степень четвёрки, потому что дробление узла делит площадь на 4. */
  int64_t *areahist;
  double *areapix;
} litctx;

static int lit_find(const litctx *L, const hz_dcref *r) {
  hz_dccell c;
  int sh = 0;
  for (int32_t s = r->size; s > 1; s >>= 1)
    sh++;
  c.lvl = (uint8_t)(L->S->lev - sh);
  for (int a = 0; a < 3; a++)
    c.lo[a] = (uint16_t)r->lo[a];
  uint64_t k = cellkey(&c);
  int32_t lo = 0, hi = L->S->n - 1;
  while (lo <= hi) {
    int32_t mid = lo + (hi - lo) / 2;
    if (L->key[mid] == k) return L->ord[mid];
    if (L->key[mid] < k)
      lo = mid + 1;
    else
      hi = mid - 1;
  }
  return -1;
}

/* Р3 (§581): `by0..by1` — ПОЛОСА ЭКРАНА, за которую отвечает поток; `by0 > by1`
 * значит без ограничения. Пиксель принадлежит РОВНО ОДНОЙ полосе, поэтому
 * z-буфер идёт без гонок и без атомарных операций, а порядок детерминирован —
 * это и даёт побитовость. */
static void lit_tri(litctx *L, const double p[3][3], const double col[3][3], const double uv[3][2],
                    int mat, int by0, int by1) {
  const tr3_camera *cm = L->cam;
  double sx[3], sy[3], sz[3];
  for (int k = 0; k < 3; k++) {
    double d[3];
    for (int q = 0; q < 3; q++)
      d[q] = p[k][q] - cm->eye[q];
    double zz = d[0] * cm->fwd[0] + d[1] * cm->fwd[1] + d[2] * cm->fwd[2];
    if (!(zz > 1e-6)) return;
    double rr = d[0] * cm->right[0] + d[1] * cm->right[1] + d[2] * cm->right[2];
    double uu = d[0] * cm->up[0] + d[1] * cm->up[1] + d[2] * cm->up[2];
    sz[k] = zz;
    sx[k] = ((rr / (zz * cm->tanx)) + 1.0) * 0.5 * (double)cm->w;
    sy[k] = (1.0 - uu / (zz * cm->tany)) * 0.5 * (double)cm->h;
  }
  double x0 = sx[0], x1 = sx[0], y0 = sy[0], y1 = sy[0];
  for (int k = 1; k < 3; k++) {
    if (sx[k] < x0) x0 = sx[k];
    if (sx[k] > x1) x1 = sx[k];
    if (sy[k] < y0) y0 = sy[k];
    if (sy[k] > y1) y1 = sy[k];
  }
  int ix0 = (int)floor(x0), ix1 = (int)ceil(x1), iy0 = (int)floor(y0), iy1 = (int)ceil(y1);
  if (ix0 < 0) ix0 = 0;
  if (iy0 < 0) iy0 = 0;
  if (ix1 >= L->w) ix1 = L->w - 1;
  if (iy1 >= L->h) iy1 = L->h - 1;
  if (by0 <= by1) {
    if (iy0 < by0) iy0 = by0;
    if (iy1 > by1) iy1 = by1;
    if (iy0 > iy1) return;
  }
  double d21x = sx[1] - sx[0], d21y = sy[1] - sy[0];
  double d31x = sx[2] - sx[0], d31y = sy[2] - sy[0];
  double det = d21x * d31y - d21y * d31x;
  if (!(fabs(det) > 0.0)) return;
  /* Р1 (§572): `1/det` ВЫНЕСЕНО. Прежде на каждый пиксель ограничивающей
   * коробки — включая отвергнутые — приходилось ДВА деления; теперь два
   * умножения. */
  double inv = 1.0 / det;
  double dz1 = sz[1] - sz[0], dz2 = sz[2] - sz[0];
  for (int py = iy0; py <= iy1; py++)
    for (int px = ix0; px <= ix1; px++) {
      double qx = (double)px + 0.5 - sx[0], qy = (double)py + 0.5 - sy[0];
      double u = (qx * d31y - qy * d31x) * inv, v = (qy * d21x - qx * d21y) * inv;
      if (u < 0.0 || v < 0.0 || u + v > 1.0) continue;
      double zz = sz[0] + u * dz1 + v * dz2;
      size_t k = (size_t)py * (size_t)L->w + (size_t)px;
      if (zz >= L->z[k]) continue;
      L->z[k] = zz;
      L->nfrag++;
      /* Р3 (§572): ОТЛОЖЕННОЕ ЗАТЕНЕНИЕ. Здесь только запоминается, ЧЕМ пиксель
       * закрыт; цвет и гамма считаются ОДИН раз на видимый пиксель после
       * обхода. Прежде гамма платилась за каждый прошедший z фрагмент, а
       * большая часть их затиралась следующими треугольниками — то есть работа
       * делалась и выбрасывалась. */
      for (int c = 0; c < 3; c++)
        L->defcol[3 * k + (size_t)c] =
            (float)(col[0][c] + u * (col[1][c] - col[0][c]) + v * (col[2][c] - col[0][c]));
      if (L->defuv != NULL) {
        for (int c = 0; c < 2; c++)
          L->defuv[2 * k + (size_t)c] =
              (float)(uv[0][c] + u * (uv[1][c] - uv[0][c]) + v * (uv[2][c] - uv[0][c]));
        L->defmat[k] = (unsigned char)mat;
      }
    }
}

/* Р2 (§572): ГАММА ТАБЛИЦЕЙ. Выход всё равно байт, поэтому `4096` шагов дают
 * ошибку ниже половины кванта ПО ПОСТРОЕНИЮ, а не по замеру. `pow` при этом
 * зовётся `4096` раз на кадр вместо трёх раз на фрагмент.
 * `HZ_GAMN` — размер таблицы, а не порог: он назван здесь, потому что от него
 * зависит точность, и негативный контроль (`gam16`) её ломает нарочно. */
#define HZ_GAMN 4096
/* §791: предвычисленные данные треугольника КЛЕТКИ для surf-uv резолва.
 * Соседние пиксели строки почти всегда падают в ту же клетку — нормали,
 * барицентрические скаляры и uv-плотность считаются один раз на смену
 * клетки, а не на каждый пиксель × каждый треугольник (замерено §790-доп:
 * 98 % растра сидело именно здесь, texflat 8998 → 197 мс). */
typedef struct {
  double n[3], A[3], e1[3], e2[3];
  double d11, d12, d22, dn3;
  int32_t q0, q1, q2, mt;
  double uvdens;
} uvtri791;

static void lit_resolve(litctx *L, int gamn) {
  /* БЕЛАЯ ТОЧКА — СВОЙСТВО КАДРА, А НЕ СЦЕНЫ (08-12). Прежде она бралась
   * перцентилем по ВСЕМ ячейкам среза, включая невидимые и залитые солнцем
   * снаружи; при камере внутри двора это давило видимый интерьер в чёрное —
   * замерено: кадр выходил сплошь тёмным при `E_ind/E_dir = 1.59`.
   * С отложенным буфером (Р3) правильная величина под рукой: перцентиль по
   * ЗАКРЫТЫМ пикселям. Перцентиль, а не максимум, — по той же причине, что и
   * раньше: одиночный блик не должен утопить кадр. */
  size_t npx = (size_t)L->w * (size_t)L->h;
  {
    float *v = malloc(npx * sizeof *v);
    if (v == NULL) exit(1);
    size_t nv = 0;
    for (size_t k = 0; k < npx; k++) {
      if (L->z[k] >= 1e299) continue;
      float mx = L->defcol[3 * k];
      for (int c = 1; c < 3; c++)
        if (L->defcol[3 * k + (size_t)c] > mx) mx = L->defcol[3 * k + (size_t)c];
      v[nv++] = mx;
    }
    if (nv > 0) {
      qsort(v, nv, sizeof *v, cmp_f);
      double w995 = (double)v[(size_t)((double)nv * 0.995)];
      if (w995 > 0.0) L->white = w995;
    }
    free(v);
  }
  /* Ш8 (§575): ВЫБОРКА ТЕКСТУРЫ — ОДИН РАЗ НА ВИДИМЫЙ ПИКСЕЛЬ. Отложенный
   * буфер (Р3 §572), сделанный ради скорости, здесь окупается второй раз:
   * перекрытые фрагменты текстуру не читают вовсе.
   * Правило Т1: текстура умножается на радианс ПРИ ЧТЕНИИ поля. В перенос она
   * не входит, поэтому цветного непрямого света от неё не будет — цена названа
   * в §575 и платится сознательно.
   * `v` растёт вверх в OBJ и вниз в изображении — отсюда `1 − v`. */
  if (L->texrgb != NULL && L->defuv != NULL) {
    int64_t nsurf791 = 0, ntex791 = 0;
#pragma omp parallel reduction(+ : nsurf791, ntex791)
    {
      /* §791: кэш треугольников ТЕКУЩЕЙ КЛЕТКИ потока — соседние пиксели
       * строки почти всегда в той же клетке; пересборка только на смене. */
      int32_t ccl[3] = {0, 0, 0};
      int cvalid = 0, cnt791 = 0, ccap = 0;
      uvtri791 *ctri = NULL;
      /* §812: потоковые счётчики блик-канала (сведение — критической секцией
       * ниже; гонка в горячем цикле недопустима, А1036) */
      glstat812 st812;
      memset(&st812, 0, sizeof st812);
#pragma omp for schedule(static)
      for (int py2 = 0; py2 < L->h; py2++) {
        for (int px2 = 0; px2 < L->w; px2++) {
          size_t k = (size_t)py2 * (size_t)L->w + (size_t)px2;
          if (L->z[k] >= 1e299) continue;
          int mt = L->defmat[k];
          double uu = (double)L->defuv[2 * k + 0], vv0 = (double)L->defuv[2 * k + 1];
          double uvdens = 0.0;
          /* КООРДИНАТА БЕРЁТСЯ С САМОЙ ПОВЕРХНОСТИ, А НЕ С НАШЕГО МНОГОУГОЛЬНИКА
           * (замечание пользователя 08-12): по глубине восстанавливается мировая
           * точка пикселя, по ней — ЯЧЕЙКА ПОЛЯ и её треугольники, uv — у того,
           * к чьей плоскости точка ближе. §791: та же семантика argmin, но
           * пер-треугольная тригонометрия предвычислена на клетку, барицентрика
           * считается только победителю. */
          if (L->mesh != NULL && L->ct != NULL) {
            const tr3_camera *cm2 = L->cam;
            double ax = ((double)px2 + 0.5) / (double)L->w * 2.0 - 1.0;
            double ay = 1.0 - ((double)py2 + 0.5) / (double)L->h * 2.0;
            double dr[3], wp[3];
            for (int c = 0; c < 3; c++)
              dr[c] = cm2->fwd[c] + cm2->right[c] * ax * cm2->tanx + cm2->up[c] * ay * cm2->tany;
            for (int c = 0; c < 3; c++)
              wp[c] = cm2->eye[c] + dr[c] * L->z[k];
            int32_t cl3[3] = {0, 0, 0}; /* явный ноль: анализатор не видит
                                         * заполнения через тернарник ok3 */
            int ok3 = 1;
            for (int c = 0; c < 3; c++) {
              double f3 = floor((wp[c] - L->fr->org[c]) / L->fr->h);
              if (!(f3 >= 0.0) || !(f3 < (double)L->fr->n)) ok3 = 0;
              cl3[c] = ok3 ? (int32_t)f3 : 0;
            }
            if (ok3 && (!cvalid || cl3[0] != ccl[0] || cl3[1] != ccl[1] || cl3[2] != ccl[2])) {
              const int32_t *ls3 = NULL;
              int32_t nls = ct_list(L->ct, cl3, &ls3);
              cnt791 = 0;
              for (int32_t t3 = 0; t3 < nls; t3++) {
                const double *A3, *B3, *C3;
                tri_verts(L->mesh, ls3[t3], &A3, &B3, &C3);
                /* нулевая инициализация — ложный класс gcc-analyzer: внутри
                 * omp-региона он теряет безусловное присваивание через
                 * указатели tri_verts (CLAUDE.md, класс diam 07-24) */
                double e1[3] = {0, 0, 0}, e2[3] = {0, 0, 0}, nn3[3];
                for (int c = 0; c < 3; c++) {
                  e1[c] = B3[c] - A3[c];
                  e2[c] = C3[c] - A3[c];
                }
                nn3[0] = e1[1] * e2[2] - e1[2] * e2[1];
                nn3[1] = e1[2] * e2[0] - e1[0] * e2[2];
                nn3[2] = e1[0] * e2[1] - e1[1] * e2[0];
                double nl3 = sqrt(nn3[0] * nn3[0] + nn3[1] * nn3[1] + nn3[2] * nn3[2]);
                if (!(nl3 > 0.0)) continue; /* вырожденный: плоскости нет */
                const int32_t *ft3 = L->mesh->ft;
                if (ft3 == NULL) continue;
                int32_t q0 = ft3[3 * (size_t)ls3[t3] + 0], q1 = ft3[3 * (size_t)ls3[t3] + 1],
                        q2 = ft3[3 * (size_t)ls3[t3] + 2];
                if (q0 < 0 || q1 < 0 || q2 < 0) continue;
                double d11 = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
                double d12 = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
                double d22 = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
                double dn3 = d11 * d22 - d12 * d12;
                if (!(fabs(dn3) > 0.0)) continue;
                if (cnt791 >= ccap) {
                  ccap = ccap ? ccap * 2 : 64;
                  uvtri791 *nc9 = realloc(ctri, (size_t)ccap * sizeof *nc9);
                  if (nc9 == NULL) exit(1);
                  ctri = nc9;
                }
                uvtri791 *u9 = &ctri[cnt791++];
                for (int c = 0; c < 3; c++) {
                  u9->n[c] = nn3[c] / nl3;
                  u9->A[c] = A3[c];
                  u9->e1[c] = e1[c];
                  u9->e2[c] = e2[c];
                }
                /* §812 НК-б: перевёрнутая нормаль у каждого 100-го
                 * треугольника. Расстояние argmin не меняется (|dot|), uv и
                 * текстура не меняются — шевелится ТОЛЬКО блик-канал. */
                if (g_glnorm812 && (ls3[t3] % 100) == 0)
                  for (int c = 0; c < 3; c++)
                    u9->n[c] = -u9->n[c];
                u9->d11 = d11;
                u9->d12 = d12;
                u9->d22 = d22;
                u9->dn3 = dn3;
                u9->q0 = q0;
                u9->q1 = q1;
                u9->q2 = q2;
                u9->mt = L->mesh->fm != NULL ? L->mesh->fm[ls3[t3]] : -1;
                /* ПЛОТНОСТЬ ТЕКСЕЛЕЙ НА МЕТР — константа треугольника: отношение
                 * площади в uv к площади в мире; предвычисляется здесь же. */
                double au = L->mesh->vt[2 * (size_t)q1 + 0] - L->mesh->vt[2 * (size_t)q0 + 0];
                double av = L->mesh->vt[2 * (size_t)q1 + 1] - L->mesh->vt[2 * (size_t)q0 + 1];
                double bu2 = L->mesh->vt[2 * (size_t)q2 + 0] - L->mesh->vt[2 * (size_t)q0 + 0];
                double bv2 = L->mesh->vt[2 * (size_t)q2 + 1] - L->mesh->vt[2 * (size_t)q0 + 1];
                double auv = fabs(au * bv2 - av * bu2);
                double aw3 = sqrt(dn3);
                u9->uvdens = (aw3 > 0.0) ? sqrt(auv / aw3) : 0.0;
              }
              ccl[0] = cl3[0];
              ccl[1] = cl3[1];
              ccl[2] = cl3[2];
              cvalid = 1;
            }
            int best9 = -1;
            double bestd3 = 1e300;
            if (ok3 && cvalid)
              for (int j9 = 0; j9 < cnt791; j9++) {
                const uvtri791 *u8 = &ctri[j9];
                double dd3 = fabs((wp[0] - u8->A[0]) * u8->n[0] + (wp[1] - u8->A[1]) * u8->n[1] +
                                  (wp[2] - u8->A[2]) * u8->n[2]);
                if (dd3 < bestd3) {
                  bestd3 = dd3;
                  best9 = j9;
                }
              }
            if (best9 >= 0) { /* барицентрика — ТОЛЬКО победителю */
              const uvtri791 *u9 = &ctri[best9];
              double vp[3];
              for (int c = 0; c < 3; c++)
                vp[c] = wp[c] - u9->A[c];
              double dp1 = vp[0] * u9->e1[0] + vp[1] * u9->e1[1] + vp[2] * u9->e1[2];
              double dp2 = vp[0] * u9->e2[0] + vp[1] * u9->e2[1] + vp[2] * u9->e2[2];
              double bu = (u9->d22 * dp1 - u9->d12 * dp2) / u9->dn3;
              double bv = (u9->d11 * dp2 - u9->d12 * dp1) / u9->dn3;
              uu =
                  L->mesh->vt[2 * (size_t)u9->q0 + 0] +
                  bu * (L->mesh->vt[2 * (size_t)u9->q1 + 0] - L->mesh->vt[2 * (size_t)u9->q0 + 0]) +
                  bv * (L->mesh->vt[2 * (size_t)u9->q2 + 0] - L->mesh->vt[2 * (size_t)u9->q0 + 0]);
              vv0 =
                  L->mesh->vt[2 * (size_t)u9->q0 + 1] +
                  bu * (L->mesh->vt[2 * (size_t)u9->q1 + 1] - L->mesh->vt[2 * (size_t)u9->q0 + 1]) +
                  bv * (L->mesh->vt[2 * (size_t)u9->q2 + 1] - L->mesh->vt[2 * (size_t)u9->q0 + 1]);
              mt = u9->mt;
              if (mt < 0 || mt >= L->nmtl) mt = L->defmat[k];
              uvdens = u9->uvdens;
              nsurf791++;
              /* §812: СПЕК-ЧЛЕН — по авторской нормали и материалу победившего
               * треугольника, ДО текстурного continue (материал без текстуры
               * бликует тоже). Буфер глобальный, индексация как defcol;
               * сведение — в критической секции после цикла. */
              if (g_glet != NULL && g_glspec != NULL) {
                double ks812 = 0.0, ns812 = 0.0;
                if (mt >= 0 && mt < L->nmtl) {
                  const double *ks3 = L->mesh->mtl[mt].ks3;
                  ks812 = (ks3[0] + ks3[1] + ks3[2]) / 3.0;
                  ns812 = L->mesh->mtl[mt].ns;
                }
                if (g_glks0) ks812 = 0.0;
                if (ks812 > 0.0 && ns812 > 0.0) {
                  double t812a = now_s();
                  double out812[3];
                  gloss812_pixel(g_glet, wp, u9->n, dr, ks812, ns812, 1, out812, &st812);
                  for (int c = 0; c < 3; c++)
                    g_glspec[3 * k + (size_t)c] = (float)out812[c];
                  st812.t += now_s() - t812a;
                }
              }
            }
          }
          if (L->texrgb[mt] == NULL) continue;
          double vv = 1.0 - vv0;
          /* УРОВЕНЬ ПИРАМИДЫ по следу пикселя в текселях. Ширина пикселя на
           * глубине z есть z·pxrad; умноженная на плотность текселей и размер
           * текстуры, она даёт след; log2 — уровень. */
          double lodf = 0.0;
          if (!L->nomip && uvdens > 0.0 && L->pxrad > 0.0) {
            double foot = L->z[k] * L->pxrad * uvdens * (double)L->texw[mt];
            if (foot > 1.0) lodf = log2(foot);
          }
          int nl4 = L->nmip[mt] > 0 ? L->nmip[mt] : 1;
          if (lodf > (double)(nl4 - 1)) lodf = (double)(nl4 - 1);
          int l0 = (int)lodf, l1 = l0 + 1 < nl4 ? l0 + 1 : l0;
          double fl = lodf - (double)l0;
          /* ТРИЛИНЕЙНО: билинейно внутри двух уровней и линейно между ними. */
          double acc3[3] = {0.0, 0.0, 0.0};
          for (int s6 = 0; s6 < 2; s6++) {
            int lv = s6 ? l1 : l0;
            double wl = s6 ? fl : 1.0 - fl;
            if (!(wl > 0.0)) continue;
            const unsigned char *tx = L->texmip[mt][lv];
            int tw = L->mipw[mt][lv], th = L->miph[mt][lv];
            double fu = uu - floor(uu), fv = vv - floor(vv);
            double gx = fu * (double)tw - 0.5, gy = fv * (double)th - 0.5;
            int x0 = (int)floor(gx), y0 = (int)floor(gy);
            double tx0 = gx - (double)x0, ty0 = gy - (double)y0;
            for (int dy = 0; dy < 2; dy++)
              for (int dx = 0; dx < 2; dx++) {
                int xx = x0 + dx, yy = y0 + dy;
                /* Повтор по краю: текстуры тайловые, обрезка дала бы шов. */
                xx = ((xx % tw) + tw) % tw;
                yy = ((yy % th) + th) % th;
                double wq = (dx ? tx0 : 1.0 - tx0) * (dy ? ty0 : 1.0 - ty0) * wl;
                const unsigned char *px = tx + 3 * ((size_t)yy * (size_t)tw + (size_t)xx);
                for (int c = 0; c < 3; c++)
                  acc3[c] += wq * (double)px[c];
              }
          }
          for (int c = 0; c < 3; c++)
            L->defcol[3 * k + (size_t)c] *= (float)(acc3[c] / 255.0);
          ntex791++;
        }
      }
      free(ctri);
      /* §812: сведение потоковых счётчиков блик-канала (по одному входу на
       * поток — не в горячем цикле) */
      if (g_glet != NULL) {
#pragma omp critical
        {
          g_glprune += st812.prune;
          g_glleaf += st812.leaf;
          g_glmarch += st812.march;
          g_glt812 += st812.t;
        }
      }
    }
    L->nsurfuv += nsurf791;
    L->ntexpx += ntex791;
  }
  unsigned char *lut = malloc((size_t)gamn);
  if (lut == NULL) exit(1);
  for (int i = 0; i < gamn; i++) {
    double t = ((double)i + 0.5) / (double)gamn;
    lut[i] = (unsigned char)(pow(t, 1.0 / 2.2) * 255.0 + 0.5);
  }
  size_t np = (size_t)L->w * (size_t)L->h;
  for (size_t k = 0; k < np; k++) {
    if (L->z[k] >= 1e299) continue; /* пиксель не закрыт ничем */
    /* §812: спек — ПОСЛЕ текстуры (блик не крашится ею), белая точка
     * ДИФФУЗНАЯ (А1364): блик честно клипуется насыщением. Доли и клип —
     * в прибор (сериальный цикл, глобальные счётчики без гонок). */
    if (g_glspec != NULL) {
      double sp8 =
          (double)g_glspec[3 * k] + (double)g_glspec[3 * k + 1] + (double)g_glspec[3 * k + 2];
      double dd8 =
          (double)L->defcol[3 * k] + (double)L->defcol[3 * k + 1] + (double)L->defcol[3 * k + 2];
      if (sp8 > 0.01 * (dd8 + sp8)) g_glnpx1++;
      if (sp8 > 0.1 * (dd8 + sp8)) g_glnpx10++;
      g_glsum += sp8;
    }
    for (int c = 0; c < 3; c++) {
      /* без ключа — ИСХОДНОЕ выражение без добавок: мир посимвольный */
      double t = (double)L->defcol[3 * k + (size_t)c] / L->white;
      if (g_glspec != NULL) t += (double)g_glspec[3 * k + (size_t)c] / L->white;
      if (t < 0.0) t = 0.0;
      if (t > 1.0) {
        t = 1.0;
        /* клип считается только там, где спек ДОБАВИЛ (иначе это диффузный
         * клип, не блик) */
        if (g_glspec != NULL && g_glspec[3 * k + (size_t)c] > 0.0f) g_glclip++;
      }
      int ix = (int)(t * (double)gamn);
      if (ix >= gamn) ix = gamn - 1;
      L->rgb[3 * k + (size_t)c] = lut[ix];
    }
  }
  free(lut);
}

/* ОТСЕЧЕНИЕ ДО ПРОЕКЦИИ (Ш5в). Срез строится на ПОЛНЫЙ ШАР — так и задумано
 * (§383: свет приходит и из-за спины), но КАМЕРНЫЙ проход обязан брать только
 * видимое. Прежде каждый многоугольник проецировался целиком, и лишь потом
 * выяснялось, что он за камерой или за краем экрана.
 *
 * Проба — по четырём плоскостям пирамиды видимости плюс ближняя: если ВСЕ углы
 * снаружи одной и той же плоскости, многоугольник отбрасывается. Это ТОЧНОЕ
 * отсечение в одну сторону: отбрасывается только заведомо невидимое, поэтому
 * картинка обязана совпасть ПОБИТОВО — она и есть приёмка. */
static int lit_cull(const litctx *L, const double w[4][3], int nv) {
  const tr3_camera *cm = L->cam;
  int out_near = 0, out_l = 0, out_r = 0, out_b = 0, out_t = 0;
  for (int i = 0; i < nv; i++) {
    double d[3];
    for (int c = 0; c < 3; c++)
      d[c] = w[i][c] - cm->eye[c];
    double zz = d[0] * cm->fwd[0] + d[1] * cm->fwd[1] + d[2] * cm->fwd[2];
    double rr = d[0] * cm->right[0] + d[1] * cm->right[1] + d[2] * cm->right[2];
    double uu = d[0] * cm->up[0] + d[1] * cm->up[1] + d[2] * cm->up[2];
    if (!(zz > 1e-6)) out_near++;
    if (rr < -zz * cm->tanx) out_l++;
    if (rr > zz * cm->tanx) out_r++;
    if (uu < -zz * cm->tany) out_b++;
    if (uu > zz * cm->tany) out_t++;
  }
  return out_near == nv || out_l == nv || out_r == nv || out_b == nv || out_t == nv;
}

/* §573: пустой обработчик — им меряется цена САМОГО обхода дерева и критерия
 * LOD, отдельно от растеризации. */
static int lit_none(void *ctx, const hz_dcref *ref, const double (*v)[3], int nv) {
  (void)ctx;
  (void)ref;
  (void)v;
  (void)nv;
  return 0;
}

/* КОНТРОЛЬНАЯ СУММА ПОТОКА МНОГОУГОЛЬНИКОВ (А1013). Картинка к ПОРЯДКУ выдачи
 * слепа: треугольники собираются в массив и рисуются по полосам, а z-буфер
 * порядок не различает, пока глубины не равны. Значит «картинка побитово та же»
 * доказывает МЕНЬШЕ, чем кажется, — а правка Р2 §585 трогает как раз то место,
 * где порядок мог бы поехать. Поэтому здесь считается ПОРЯДКОВАЯ сумма: FNV-1a
 * по битовым образцам координат В ПОРЯДКЕ ВЫДАЧИ. Считается ДО отсечения — она
 * должна отвечать за выход ОБХОДА, а не за работу растеризатора. */
#define HZ_FNV_PRIME UINT64_C(1099511628211)
#define HZ_FNV_BASIS UINT64_C(14695981039346656037)

static void polysum_add(uint64_t *h, const double (*v)[3], int nv) {
  unsigned char b[8];
  *h = (*h ^ (uint64_t)nv) * HZ_FNV_PRIME;
  for (int i = 0; i < nv; i++)
    for (int c = 0; c < 3; c++) {
      memcpy(b, &v[i][c], sizeof b);
      for (int k = 0; k < 8; k++)
        *h = (*h ^ (uint64_t)b[k]) * HZ_FNV_PRIME;
    }
}

static int lit_poly(void *ctx, const hz_dcref *ref, const double (*v)[3], int nv) {
  litctx *L = (litctx *)ctx;
  double w[4][3], col[4][3];

  polysum_add(&L->polysum, v, nv);
  for (int i = 0; i < nv; i++)
    for (int c = 0; c < 3; c++)
      w[i][c] = L->fr->org[c] + v[i][c] * L->fr->h;
  L->nseen++;
  if (!L->nocull && lit_cull(L, w, nv)) {
    L->ncull++;
    return 0;
  }
  double uvv[4][2];
  int matp = 0;
  for (int i = 0; i < nv; i++) {
    int idx = lit_find(L, &ref[i]);
    for (int c = 0; c < 3; c++)
      col[i][c] = idx >= 0 ? (double)L->irr[3 * (size_t)idx + (size_t)c] : 0.0;
    uvv[i][0] = uvv[i][1] = 0.0;
    if (idx >= 0 && L->uvs != NULL) {
      uvv[i][0] = (double)L->uvs[2 * (size_t)idx + 0];
      uvv[i][1] = (double)L->uvs[2 * (size_t)idx + 1];
      matp = L->S->c[idx].mat;
    }
  }
  for (int i = 1; i + 1 < nv; i++) {
    double p3[3][3], c3[3][3], u3[3][2];
    for (int c = 0; c < 3; c++) {
      p3[0][c] = w[0][c];
      p3[1][c] = w[i][c];
      p3[2][c] = w[i + 1][c];
      c3[0][c] = col[0][c];
      c3[1][c] = col[i][c];
      c3[2][c] = col[i + 1][c];
    }
    for (int c = 0; c < 2; c++) {
      u3[0][c] = uvv[0][c];
      u3[1][c] = uvv[i][c];
      u3[2][c] = uvv[i + 1][c];
    }
    if (L->tris != NULL) {
      if (L->ntris >= L->captris) {
        int64_t nc2 = L->captris > 0 ? L->captris * 2 : 65536;
        struct littri *nt = realloc(L->tris, (size_t)nc2 * sizeof *nt);
        if (nt == NULL) exit(1);
        L->tris = nt;
        L->captris = nc2;
      }
      struct littri *dst = &L->tris[L->ntris++];
      memcpy(dst->p, p3, sizeof p3);
      memcpy(dst->col, c3, sizeof c3);
      memcpy(dst->uv, u3, sizeof u3);
      dst->mat = matp;
      {
        const tr3_camera *cm3 = L->cam;
        double y0f = 1e300, y1f = -1e300;
        int okp = 1;
        for (int q3 = 0; q3 < 3; q3++) {
          double d3[3];
          for (int c3i = 0; c3i < 3; c3i++)
            d3[c3i] = p3[q3][c3i] - cm3->eye[c3i];
          double zz3 = d3[0] * cm3->fwd[0] + d3[1] * cm3->fwd[1] + d3[2] * cm3->fwd[2];
          if (!(zz3 > 1e-6)) {
            okp = 0;
            break;
          }
          double uu3 = d3[0] * cm3->up[0] + d3[1] * cm3->up[1] + d3[2] * cm3->up[2];
          double sy3 = (1.0 - uu3 / (zz3 * cm3->tany)) * 0.5 * (double)cm3->h;
          if (sy3 < y0f) y0f = sy3;
          if (sy3 > y1f) y1f = sy3;
        }
        dst->iy0 = okp ? (int)floor(y0f) : 0;
        dst->iy1 = okp ? (int)ceil(y1f) : L->h - 1;
        /* РАЗМЕР ПОЛИГОНА НА ЭКРАНЕ — гистограмма по площади в пикселях.
         * Замечание пользователя 08-13: дефекты видны там, где полигон КРУПНЫЙ,
         * и вопрос «дробить или интерполировать тоньше» решается этим числом, а
         * не на глаз. Считается ЗДЕСЬ, при сборе: в отрисовке по полосам один
         * треугольник попадает в несколько потоков, и счёт был бы и гонкой, и
         * многократным. */
        if (L->areahist != NULL && okp) {
          double x0f = 1e300, x1f = -1e300;
          for (int q3 = 0; q3 < 3; q3++) {
            double d3[3];
            for (int c3i = 0; c3i < 3; c3i++)
              d3[c3i] = p3[q3][c3i] - cm3->eye[c3i];
            double zz3 = d3[0] * cm3->fwd[0] + d3[1] * cm3->fwd[1] + d3[2] * cm3->fwd[2];
            double rr3 = d3[0] * cm3->right[0] + d3[1] * cm3->right[1] + d3[2] * cm3->right[2];
            double sx3 = ((rr3 / (zz3 * cm3->tanx)) + 1.0) * 0.5 * (double)cm3->w;
            if (sx3 < x0f) x0f = sx3;
            if (sx3 > x1f) x1f = sx3;
          }
          double ar4 = 0.5 * (x1f - x0f) * (y1f - y0f);
          if (ar4 > 0.0) {
            int b4 = 0;
            double t4 = ar4;
            while (b4 < 9 && t4 >= 4.0) {
              t4 /= 4.0;
              b4++;
            }
            L->areahist[b4]++;
            L->areapix[b4] += ar4;
          }
        }
      }
    } else
      lit_tri(L, p3, c3, u3, matp, 0, -1);
  }
  return 0;
}
/* --- 4. главная ------------------------------------------------------------ */

/* ======================= ХОДЬБА ПО СЦЕНЕ (§589) =========================
 *
 * ОКНО И ВВОД ЖИВУТ ТОЛЬКО ПОД `-DHZ_SDL`, И ЭТО СОЗНАТЕЛЬНО. Рабочий
 * `build/pfield` остаётся без единой внешней зависимости — он замерный
 * инструмент, и тащить в него оконную библиотеку значит менять условия всех
 * прежних прогонов. Ходилка собирается отдельной целью `build/pwalk` из ТОГО ЖЕ
 * файла: одна реализация растеризатора, а не две.
 *
 * РАСКЛАДКА — DESCENT, шесть степеней свободы, как названо пользователем:
 *     мышь        поворот: вправо/влево — рыскание, вверх/вниз — тангаж
 *     W A S D     СКОЛЬЖЕНИЕ в плоскости вида: вверх, влево, вниз, вправо
 *     пробел      тяга ВПЕРЁД
 *     левый Shift тяга НАЗАД
 *     Q E         КРЕН относительно вида
 *     колесо      скорость (шаг вдвое), Tab — вернуть 6 м/с
 *     Esc         выход, ` (тильда) — отпустить/схватить мышь
 *
 * ПОЧЕМУ БАЗИС ХРАНИТСЯ ЦЕЛИКОМ, А НЕ УГЛАМИ. При крене «верх мира» перестаёт
 * быть верхом камеры, и пара (углы Эйлера + фиксированный up) шесть степеней
 * свободы не выражает вовсе: у неё нет крена по построению. Поэтому хранится
 * тройка ортонормированных векторов, повороты — их вращения, а ортогональность
 * восстанавливается КАЖДЫЙ кадр (иначе накопление ошибки за тысячу кадров
 * уводит базис, и это видно как медленный завал горизонта).
 */
#ifdef HZ_SDL
#include <SDL2/SDL.h>

static SDL_Window *g_win;
static SDL_Renderer *g_ren;
static SDL_Texture *g_tex;
static int g_texres;
static int g_texresh;
/* Базис камеры. Заводится из `eyec/atc/upc` при первом кадре. */
static double g_wpos[3], g_wfwd[3], g_wup[3], g_wright[3];
static int g_winit;
static double g_wspeed = 6.0; /* м/с; сцена Bistro 108 м поперёк — 18 с на проход */
static int g_wgrab = 1;

/* Скорость поворота мыши: радиан на пиксель. Число не магическое — это поле
 * зрения, делённое на сторону окна, то есть «пиксель мыши = пиксель экрана». */
#define HZ_WALK_MOUSE_RAD_PER_PX (HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0 / 512.0)
/* Крен клавишей: радиан в секунду. */
#define HZ_WALK_ROLL_RAD_PER_S 1.5

static void wnorm(double v[3]) {
  double s = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  if (s > 0.0)
    for (int a = 0; a < 3; a++)
      v[a] /= s;
}

static void wcross(const double a[3], const double b[3], double o[3]) {
  o[0] = a[1] * b[2] - a[2] * b[1];
  o[1] = a[2] * b[0] - a[0] * b[2];
  o[2] = a[0] * b[1] - a[1] * b[0];
}

/* Поворот вектора `v` вокруг оси `k` (единичной) на угол `t` — формула Родрига.
 * Одна формула на все три поворота; отдельных матриц для рыскания, тангажа и
 * крена не заводится, потому что разница между ними только в оси. */
static void wrot(double v[3], const double k[3], double t) {
  double c = cos(t), s = sin(t), kv[3], d = 0.0;
  wcross(k, v, kv);
  for (int a = 0; a < 3; a++)
    d += k[a] * v[a];
  for (int a = 0; a < 3; a++)
    v[a] = v[a] * c + kv[a] * s + k[a] * d * (1.0 - c);
}

static int walk_open(int resw, int resh) {
  if (SDL_Init(SDL_INIT_VIDEO) != 0) {
    fprintf(stderr, "SDL_Init: %s\n", SDL_GetError());
    return 0;
  }
  /* ОКНО ТЯНЕТСЯ, И РАЗРЕШЕНИЕ ИДЁТ ЗА НИМ ПО-НАСТОЯЩЕМУ (§591). Без
   * `RESIZABLE` окно было фиксированным, а растянутое средствами оконного
   * менеджера давало ПИКСЕЛИ: текстура оставалась `res × res`, её просто
   * масштабировали. Теперь смена размера меняет `res` СЛЕДУЮЩЕГО кадра, то есть
   * считается настоящий кадр нового разрешения.
   * СТОРОНА БЕРЁТСЯ МЕНЬШАЯ: камера пока квадратная (`tr3_camera_look` получает
   * `res, res`), и растягивать её на неквадратное окно значило бы врать про поле
   * зрения. Неквадратный кадр — отдельная правка камеры, а не подгонка здесь. */
  g_win =
      SDL_CreateWindow("helmholz — ходьба по сцене", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                       resw, resh, SDL_WINDOW_SHOWN | SDL_WINDOW_RESIZABLE);
  if (g_win == NULL) {
    fprintf(stderr, "SDL_CreateWindow: %s\n", SDL_GetError());
    return 0;
  }
  g_ren = SDL_CreateRenderer(g_win, -1, SDL_RENDERER_ACCELERATED);
  if (g_ren == NULL) g_ren = SDL_CreateRenderer(g_win, -1, 0);
  if (g_ren == NULL) {
    fprintf(stderr, "SDL_CreateRenderer: %s\n", SDL_GetError());
    return 0;
  }
  g_tex = SDL_CreateTexture(g_ren, SDL_PIXELFORMAT_RGB24, SDL_TEXTUREACCESS_STREAMING, resw, resh);
  if (g_tex == NULL) {
    fprintf(stderr, "SDL_CreateTexture: %s\n", SDL_GetError());
    return 0;
  }
  g_texres = resw;
  g_texresh = resh;
  SDL_SetRelativeMouseMode(SDL_TRUE);
  return 1;
}

static void walk_close(void) {
  if (g_tex != NULL) SDL_DestroyTexture(g_tex);
  if (g_ren != NULL) SDL_DestroyRenderer(g_ren);
  if (g_win != NULL) SDL_DestroyWindow(g_win);
  SDL_Quit();
}

/* Показать кадр и принять ввод. Возвращает 0, если пора выходить.
 * `dt` — сколько заняло ПРЕДЫДУЩЕЕ построение кадра: движение считается по
 * времени, а не по кадрам, иначе скорость ходьбы зависела бы от того, куда
 * смотришь (у нас кадр от 0.15 до 0.9 с — разница втрое). */
static int walk_present(const unsigned char *rgb, int resw, int resh, int *resw_next,
                        int *resh_next, double eyec[3], double atc[3], double upc[3], double dt) {
  if (!g_winit) {
    g_winit = 1;
    for (int a = 0; a < 3; a++) {
      g_wpos[a] = eyec[a];
      g_wfwd[a] = atc[a] - eyec[a];
      g_wup[a] = upc[a];
    }
    wnorm(g_wfwd);
    wcross(g_wfwd, g_wup, g_wright);
    wnorm(g_wright);
    wcross(g_wright, g_wfwd, g_wup);
    wnorm(g_wup);
  }
  /* Кадр пришёл не того размера, что текстура (окно потянули, пока он считался)
   * — пересоздать под кадр, а не выходить: терять ходьбу из-за движения мышью по
   * рамке нельзя. */
  if (resw != g_texres || resh != g_texresh) {
    SDL_Texture *nt =
        SDL_CreateTexture(g_ren, SDL_PIXELFORMAT_RGB24, SDL_TEXTUREACCESS_STREAMING, resw, resh);
    if (nt == NULL) return 0;
    SDL_DestroyTexture(g_tex);
    g_tex = nt;
    g_texres = resw;
    g_texresh = resh;
  }
  SDL_UpdateTexture(g_tex, NULL, rgb, resw * 3);
  SDL_RenderClear(g_ren);
  SDL_RenderCopy(g_ren, g_tex, NULL, NULL);
  SDL_RenderPresent(g_ren);

  SDL_Event e;
  double dyaw = 0.0, dpitch = 0.0;
  int nmot = 0;
  while (SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT) return 0;
    /* СМЕНА РАЗМЕРА ОКНА = СМЕНА РАЗРЕШЕНИЯ СЛЕДУЮЩЕГО КАДРА. Текстура
     * пересоздаётся здесь же, потому что нынешний кадр в неё уже показан. */
    if (e.type == SDL_WINDOWEVENT && e.window.event == SDL_WINDOWEVENT_SIZE_CHANGED) {
      /* §631: кадр теперь НЕКВАДРАТНЫЙ, и окно берётся как есть — камера это
       * умела всегда (`tany = tanx·h/w`), упиралось в одно число разрешения. */
      int ww = e.window.data1, wh = e.window.data2;
      if (ww < 64) ww = 64;
      if (wh < 64) wh = 64;
      if (ww != g_texres || wh != g_texresh) {
        SDL_Texture *nt =
            SDL_CreateTexture(g_ren, SDL_PIXELFORMAT_RGB24, SDL_TEXTUREACCESS_STREAMING, ww, wh);
        if (nt != NULL) {
          SDL_DestroyTexture(g_tex);
          g_tex = nt;
          g_texres = ww;
          g_texresh = wh;
          *resw_next = ww;
          *resh_next = wh;
          printf("   разрешение -> %d x %d\n", ww, wh);
          fflush(stdout);
        }
      }
    }
    if (e.type == SDL_MOUSEMOTION && g_wgrab) {
      dyaw += (double)e.motion.xrel * HZ_WALK_MOUSE_RAD_PER_PX;
      dpitch += (double)e.motion.yrel * HZ_WALK_MOUSE_RAD_PER_PX;
      nmot++;
    }
    if (e.type == SDL_MOUSEWHEEL) {
      if (e.wheel.y > 0) g_wspeed *= 2.0;
      if (e.wheel.y < 0) g_wspeed *= 0.5;
      printf("   скорость %.2f м/с\n", g_wspeed);
      fflush(stdout);
    }
    if (e.type == SDL_KEYDOWN) {
      if (e.key.keysym.sym == SDLK_ESCAPE) return 0;
      if (e.key.keysym.sym == SDLK_TAB) g_wspeed = 6.0;
      if (e.key.keysym.sym == SDLK_BACKQUOTE) {
        g_wgrab = !g_wgrab;
        SDL_SetRelativeMouseMode(g_wgrab ? SDL_TRUE : SDL_FALSE);
      }
    }
  }

  /* ПОВОРОТЫ. Рыскание — вокруг СОБСТВЕННОГО верха, а не вокруг верха мира:
   * иначе при крене поворот мыши уводил бы взгляд вбок от того, что видно. */
  /* Сравнений с нулём у плавучки здесь нет СОЗНАТЕЛЬНО (гейт `-Wfloat-equal`):
   * «было ли движение» — это ЦЕЛЫЙ признак события, а не свойство числа. */
  if (nmot > 0) {
    wrot(g_wfwd, g_wup, -dyaw);
    wrot(g_wright, g_wup, -dyaw);
    wrot(g_wfwd, g_wright, -dpitch);
    wrot(g_wup, g_wright, -dpitch);
  }

  const Uint8 *ks = SDL_GetKeyboardState(NULL);
  int roll = 0;
  if (ks[SDL_SCANCODE_Q]) roll += 1;
  if (ks[SDL_SCANCODE_E]) roll -= 1;
  if (roll != 0) {
    double t = (double)roll * HZ_WALK_ROLL_RAD_PER_S * dt;
    wrot(g_wup, g_wfwd, t);
    wrot(g_wright, g_wfwd, t);
  }

  /* ОРТОГОНАЛИЗАЦИЯ КАЖДЫЙ КАДР — см. шапку: без неё базис уводит. */
  wnorm(g_wfwd);
  wcross(g_wfwd, g_wup, g_wright);
  wnorm(g_wright);
  wcross(g_wright, g_wfwd, g_wup);
  wnorm(g_wup);

  /* ТЯГА И СКОЛЬЖЕНИЕ. */
  double step = g_wspeed * dt;
  double mv[3] = {0.0, 0.0, 0.0};
  if (ks[SDL_SCANCODE_SPACE])
    for (int a = 0; a < 3; a++)
      mv[a] += g_wfwd[a];
  if (ks[SDL_SCANCODE_LSHIFT])
    for (int a = 0; a < 3; a++)
      mv[a] -= g_wfwd[a];
  if (ks[SDL_SCANCODE_W])
    for (int a = 0; a < 3; a++)
      mv[a] += g_wup[a];
  if (ks[SDL_SCANCODE_S])
    for (int a = 0; a < 3; a++)
      mv[a] -= g_wup[a];
  if (ks[SDL_SCANCODE_D])
    for (int a = 0; a < 3; a++)
      mv[a] += g_wright[a];
  if (ks[SDL_SCANCODE_A])
    for (int a = 0; a < 3; a++)
      mv[a] -= g_wright[a];
  double ml = sqrt(mv[0] * mv[0] + mv[1] * mv[1] + mv[2] * mv[2]);
  if (ml > 0.0)
    for (int a = 0; a < 3; a++)
      g_wpos[a] += mv[a] / ml * step;

  for (int a = 0; a < 3; a++) {
    eyec[a] = g_wpos[a];
    atc[a] = g_wpos[a] + g_wfwd[a];
    upc[a] = g_wup[a];
  }
  return 1;
}
#else
/* БЕЗ SDL ХОДЬБА НЕ МОЛЧИТ, А ОТКАЗЫВАЕТ. Тихо отрисовать один кадр вместо
 * запрошенного цикла — худший вид отказа: выглядит как работа. */
static int walk_open(int resw, int resh) {
  (void)resw;
  (void)resh;
  fprintf(stderr, "pfield: ключ `walk` требует сборки с SDL — собирайте `make build/pwalk`\n");
  return 0;
}
static void walk_close(void) {}
static int walk_present(const unsigned char *rgb, int resw, int resh, int *resw_next,
                        int *resh_next, double eyec[3], double atc[3], double upc[3], double dt) {
  (void)rgb;
  (void)resw;
  (void)resh;
  (void)resw_next;
  (void)resh_next;
  (void)eyec;
  (void)atc;
  (void)upc;
  (void)dt;
  return 0;
}
#endif

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "pfield ФАЙЛ.obj МАСШТАБ [lev=N] [nonrm] [occdump=ПУТЬ] [polydump=ПУТЬ]\n");
    return 2;
  }
  /* §631: РАЗРЕШЕНИЕ ПО УМОЛЧАНИЮ — FULL HD, И КАДР НЕКВАДРАТНЫЙ. Камера это
   * умела всегда (`tr3_camera_look` берёт `w` и `h` порознь и держит пиксель
   * квадратным: `tany = tanx·h/w`); упиралось всё в то, что в этом файле
   * разрешение было ОДНИМ числом. Замечание пользователя 08-13: «у меня 4k, на
   * нём 512x512 смотрится ну такое себе; мы же Full HD делаем». */
  int lev = 6, nonrm = 0, vq1 = 0, nodevq1 = 0, hit = 0, nofix = 0, lit = 0, sweepaxis = 0;
  int resw = 1920, resh = 1080;
  /* Ф8'-0 (§524): замер грязи после удара. `hitrad` — радиус сферы в метрах;
   * `0.20` — то же значение, на котором сняты числа 4г, и негативный контроль
   * меняет именно его. */
  int hitsweep = 0, hitnoocc = 0;
  double hitrad = 0.20;
  /* §613: ЗАСЛОНЫ ПО УМОЛЧАНИЮ. До этого шага `indvis` был ключом ЗАМЕРА, и все
   * числа отскока снимались БЕЗ заслонов — то есть `79.99 %` энергии приходило
   * сквозь стены (§612). Теперь наоборот: заслоны включены, а `noindvis`
   * возвращает прежнее поведение как негативный контроль. */
  int nrmflip = 0, nonsum = 0, indvis = 1, indmeas = 0, dosolid = 0;
  int doxfer = 0, xfernosolid = 0, doxsweep = 0, nmu = 4, xcorner = 0, xinnerfluid = 0;
  /* §667: `xhall` — прежний источник-площадка от камеры, оставлен для
   * воспроизведения старых прогонов; `xrho` — альбедо фасета, число не
   * магическое, оно ключ, и негативный контроль ставит его в ноль. */
  int xhall = 0;
  int xemitfacet = 0;            /* §670 НК: излучение по-старому, ПО ФАСЕТУ */
  int xnolim = 0;                /* §677 НК: выключить ограничитель — оператор станет ЛИНЕЙНЫМ */
  int xunit = 0;                 /* §729: применить оператор к ЕДИНИЧНОМУ состоянию (источники в
                                  * ноль, φ = bout = sout = 1); всё, что вышло > 1 при альбедо
                                  * 0.7, — локальный усилитель с адресом. С зеркалами (wall_spec)
                                  * прокидка состояния незаконна — mspec наружу не выносится
                                  * (А1108); здесь wall_spec не задаётся вовсе. */
  int xmatrho = 0, xrho_set = 0; /* §737: материальные альбедо; явный xrho= несовместим */
  int xmatfar = 0; /* §739 НК: наихудший треугольник по score — мажоранта произвола атрибуции */
  double xrhoscale = 1.0; /* §742: множитель всех альбедо — извлечение ряда отскоков */
  double xrelax = 0.0;    /* §748: демпфирование состояния; 0 — выключено */
  int xdsa = 0;           /* §750: двухсеточный цикл — грубое звено из ядра §597 */
  double xdsagain = 1.0;  /* §750: множитель поправки; 0 и −1 — НК */
  int xdsagal = 0;        /* §752: галёркинское грубое звено из самого свипа */
  int xdsalev = 5;        /* §752: сторона агрегата 2^N листьев */
  int xdsadiff = 0;       /* §754: разностные зонды вокруг рабочего поля */
  double xdsaeps = 1.0;   /* §754: множитель ε (НК: 2) */
  int xdsasplit = 0;      /* §756: элементы и стык — раздельные грубые переменные */
  int xframe = 0;         /* §758: кадр развёрткой — свиповое поле вместо ядра в irr */
  int xcontrib = 0;       /* §768: прибор вклада — перевозмущения ρ→0 по классам */
  double xcoarse = 0.0;   /* §772: метров дальности на лист размера; 0 — выключено */
  int xbounce = 0;        /* §774: лестница N прокидок-отскоков; 0 — выключено */
  /* §796: first-collision source — ПО КЛЮЧУ, а не умолчанием. План §796 менял
   * умолчание (прецедент §780), но приёмка П3 провалена (УБИВАЕТ: чешуя при
   * fc осталась), и переворачивать канон под несработавшее лечение нельзя —
   * ключ остаётся исследовательским до вердикта (класс А1287/xobjpiece). */
  int xfc = 0;
  int xfcelem = 0;              /* §798 НК: прямой канал кадра по-старому — агрегат E_fc+хвост
                                 * по ЭЛЕМЕНТАМ на узлах свипа (§796), ядровый канал не строится */
  int hcontrib_set = 0;         /* §796: задан ли hcontrib= явно (для fc-умолчания) */
  double xtailq = -1.0;         /* §774: q хвоста: <0 — измерить, 0 — усечение, >0 — НК */
  int xbcmp774 = 0;             /* §774: базовый прогон и сравнение в одном процессе */
  int xwarm_set = 0;            /* §808: стенд инкрементности — сдвиг камеры задан */
  double xwarm3[3] = {0, 0, 0}; /* сдвиг глаза кадра 2, МЕТРЫ (А1347) */
  int xwarm_n = 1;              /* корректирующих прокидок: 1 умолчание, 2 — ветка А1348 */
  int xscramble808 = 0;         /* §808 НК-1: знак у top-1 % |sout| хранения кадра 1 */
  int xs810_set = 0;            /* §810: расписание отскоков — ключ задан */
  int xs810_1 = 0, xs810_2 = 0; /* nmu такта 1 и тактов ≥2; 0,0 = ординаты канона */
  int xrhocal810 = 0;           /* §810 Р3: калибровка переиспользования альбедо */
  int xleak = 0;                /* §778: диагностический клип покрытия — адреса дыр */
  int xnopiece = 0;             /* §780 НК: раздача и рез бесконечными плоскостями, как до Р-8 */
  int xobjpiece = 0;            /* §794: куски из ТРЕУГОЛЬНИКОВ СЦЕНЫ (авторские нормали) */
  int xclus = 0;                /* §800: кластеризованные OBJ-куски (полигоны, А1285) */
  double xclustol = -1.0;       /* §800: δ слияния, м; <0 — умолчание HZ_CLUS_TOL */
  int xnomaxp = 0;              /* §735 НК: выключить принцип максимума — вернуть расходимость */
  int xcmp = 0;                 /* §744: поячеечное сличение свипа с ядром §597 */
  int32_t xchain = -1;          /* §731: трасса цепочки к ячейке — только под xunit: пол
                                 * обрыва прогулки есть уровень единичного входа */
  int32_t xcelll[4] = {-1, -1, -1, -1}; /* §733: вскрытие обновления, до 4 ячеек */
  int xdir = -1;                        /* §733: направление вскрытия */
  /* §714: альбедо стыка вынесено в ОТДЕЛЬНЫЙ ключ и по умолчанию ВЫКЛЮЧЕНО:
   * условие §707 усиливает (`bout` доходит до `4.9e+34` за один проход), и
   * держать его рабочим путём нельзя, пока не починено. */
  double xsolidrho = 0.0;
  double xsemin = 0.0; /* §717: доля h², ниже которой элемент убирается из носителя */
  double xthin = 0.0;  /* §711: доля флюида, ниже которой ячейка считается сплошной */
  int xconst = 0;      /* §699: печь на разрезанной геометрии — точное решение известно */
  int xwholemass = 0;  /* §697: подмена матрицы масс целой — различитель, прогон нефизичен */
  /* §674: потолок итераций и допуск были зашиты числами `30` и `1e-4`. Ключи
   * нужны, чтобы отличить «не сошлось» от «не дали сойтись». */
  int xit = 30;
  double xtol = 1e-4;
  /* §804/А1270: относительный критерий останова; 0 (умолчание) — абсолютный,
   * мир без ключей посимвольно прежний. Канон-команда ставит 1e-5: при
   * выросшем ×10 масштабе состояния абсолютный 1e-5 упирается в пол невязки
   * (база §804: 40 тактов при невязке 6.67e-04 и max|состояние| 161.23,
   * rel 4.1e-6 — уже сошлось). */
  double xreltol = 0.0;
  double xrho = 0.7;
  int ss2 = 0, xtrace = 0, qplane = 0;
  int sweepfrac = 1, sweepr01 = 0, nocull = 0, alb0 = 0, area = 0;
  double oven = 0.0, plates = 0.0;
  double lodthr = 1.0;
  const char *occdump = NULL, *polydump = NULL;
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0) lev = (int)strtol(argv[i] + 4, NULL, 10);
    /* НЕГАТИВНЫЙ КОНТРОЛЬ (§393): все нормали рёбер осевые. Вершины обязаны
     * остаться (ребро пересечено — вершина есть), а качество обязано упасть. */
    if (strcmp(argv[i], "nonrm") == 0) nonrm = 1;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ Ш9 (НК-2): подать в измеритель ОБРАЩЁННУЮ нормаль
     * ячейки — угол обязан стать 180° − x. Если не станет, измеритель меряет
     * длину, а не направление. */
    if (strcmp(argv[i], "nrmflip") == 0) nrmflip = 1;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ Ш9 (НК-1): вернуть прежнее поведение — крупной ячейке
     * нормаль не считать вовсе. Доля noct == 0 обязана вернуться к 100 %. */
    if (strcmp(argv[i], "nonsum") == 0) nonsum = 1;
    /* §474: `indmeas` считает долю косвенного, приходящую СКВОЗЬ заслоны;
     * `indvis` вдобавок её отбрасывает. Порознь затем, что первое — замер
     * приближения, а второе — уже другая физика. */
    if (strcmp(argv[i], "indmeas") == 0) indmeas = 1;
    /* Ш12 (§475): заливка внутренностей и сверка со знаковым объёмом. */
    if (strcmp(argv[i], "solid") == 0) dosolid = 1;
    /* Ш13 (§479): стык с переносом. `xfernosolid` — негативный контроль: без
     * маски полных ячеек объём материала обязан рухнуть. */
    if (strncmp(argv[i], "nmu=", 4) == 0) nmu = (int)strtol(argv[i] + 4, NULL, 10);
    /* Ш14 (§482): запустить саму развёртку по ординатам на стыке. */
    /* НК Ш16: чтение радианса по УГЛУ (как до §486) и наружное как ФЛЮИД. */
    /* Ш18 (§494): растр в `res`, на диск вдвое меньше свёрткой 2×2. */
    if (strcmp(argv[i], "ss2") == 0) ss2 = 1;
    if (strcmp(argv[i], "xtrace") == 0) xtrace = 1;
    if (strcmp(argv[i], "xhall") == 0) xhall = 1;
    if (strcmp(argv[i], "xemitfacet") == 0) xemitfacet = 1;
    if (strcmp(argv[i], "xnolim") == 0) xnolim = 1;
    if (strcmp(argv[i], "xunit") == 0) xunit = 1;
    if (strcmp(argv[i], "xnomaxp") == 0) xnomaxp = 1;
    if (strncmp(argv[i], "xchain=", 7) == 0) xchain = (int32_t)strtol(argv[i] + 7, NULL, 10);
    if (strncmp(argv[i], "xcell=", 6) == 0) {
      const char *pl = argv[i] + 6;
      for (int q = 0; q < 4 && *pl != 0; q++) {
        char *end;
        xcelll[q] = (int32_t)strtol(pl, &end, 10);
        pl = (*end == ',') ? end + 1 : end;
      }
    }
    if (strncmp(argv[i], "xdir=", 5) == 0) xdir = (int)strtol(argv[i] + 5, NULL, 10);
    if (strcmp(argv[i], "xwholemass") == 0) xwholemass = 1;
    if (strcmp(argv[i], "xconst") == 0) xconst = 1;
    if (strncmp(argv[i], "xthin=", 6) == 0) xthin = strtod(argv[i] + 6, NULL);
    if (strncmp(argv[i], "xsemin=", 7) == 0) xsemin = strtod(argv[i] + 7, NULL);
    if (strncmp(argv[i], "xsolidrho=", 10) == 0) xsolidrho = strtod(argv[i] + 10, NULL);
    if (strncmp(argv[i], "xit=", 4) == 0) xit = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "xtol=", 5) == 0) xtol = strtod(argv[i] + 5, NULL);
    if (strncmp(argv[i], "xrho=", 5) == 0) {
      xrho = strtod(argv[i] + 5, NULL);
      xrho_set = 1;
    }
    if (strcmp(argv[i], "xmatrho") == 0) xmatrho = 1;
    if (strncmp(argv[i], "xrhoscale=", 10) == 0) xrhoscale = strtod(argv[i] + 10, NULL);
    if (strncmp(argv[i], "xrelax=", 7) == 0) xrelax = strtod(argv[i] + 7, NULL);
    if (strncmp(argv[i], "xreltol=", 8) == 0) xreltol = strtod(argv[i] + 8, NULL);
    if (strcmp(argv[i], "xnopfm") == 0) g_nopfm = 1;
    if (strcmp(argv[i], "xdsa") == 0) {
      xdsa = 1;
      xmatrho = 1;
      doxfer = 1;
      dosolid = 1;
      doxsweep = 1;
    }
    if (strncmp(argv[i], "xdsagain=", 9) == 0) xdsagain = strtod(argv[i] + 9, NULL);
    if (strcmp(argv[i], "xdsagal") == 0) {
      xdsagal = 1;
      xdsa = 1;
      xmatrho = 1;
      doxfer = 1;
      dosolid = 1;
      doxsweep = 1;
    }
    if (strncmp(argv[i], "xdsalev=", 8) == 0) xdsalev = (int)strtol(argv[i] + 8, NULL, 10);
    if (strcmp(argv[i], "xdsadiff") == 0) {
      xdsadiff = 1;
      xdsagal = 1;
      xdsa = 1;
      xmatrho = 1;
      doxfer = 1;
      dosolid = 1;
      doxsweep = 1;
    }
    if (strncmp(argv[i], "xdsaeps=", 8) == 0) xdsaeps = strtod(argv[i] + 8, NULL);
    if (strncmp(argv[i], "xcoarse=", 8) == 0) xcoarse = strtod(argv[i] + 8, NULL);
    if (strncmp(argv[i], "xbounce=", 8) == 0) xbounce = (int)strtol(argv[i] + 8, NULL, 10);
    if (strcmp(argv[i], "xfc") == 0) xfc = 1;
    if (strcmp(argv[i], "xfcelem") == 0) xfcelem = 1;
    if (strncmp(argv[i], "xtailq=", 7) == 0) xtailq = strtod(argv[i] + 7, NULL);
    if (strcmp(argv[i], "xbcmp") == 0) xbcmp774 = 1;
    /* §808: xwarm=DX,DY,DZ[,N] — двухкадровый стенд инкрементности света;
     * необязательный N — число корректирующих прокидок (умолчание 1; 2 —
     * предсказанная ветка А1348). Требует xcoarse>0 и xbounce>0. */
    if (strncmp(argv[i], "xwarm=", 6) == 0) {
      const char *cp = argv[i] + 6;
      char *ep = NULL;
      for (int k = 0; k < 3; k++) {
        xwarm3[k] = strtod(cp, &ep);
        if (ep == cp) break;
        cp = ep;
        if (*cp == ',') cp++;
      }
      if (*cp == ',') cp++;
      if (*cp != '\0') xwarm_n = (int)strtol(cp, NULL, 10);
      xwarm_set = 1;
    }
    if (strcmp(argv[i], "xscramble") == 0) xscramble808 = 1;
    /* §810: xsched=N1,N2 — ординаты (nmu) такта 1 и тактов ≥2 лестницы
     * xbounce; 0 в поле = ордината канона (nmu). Умолчание «канон,канон» —
     * прежний мир; рабочая точка «2,1»; НК — «1,1» (тени грубятся) и «2,3»
     * (хвост дороже — двусторонняя вилка). */
    if (strncmp(argv[i], "xsched=", 7) == 0) {
      xs810_1 = (int)strtol(argv[i] + 7, NULL, 10);
      const char *cp810 = strchr(argv[i] + 7, ',');
      xs810_2 = cp810 != NULL ? (int)strtol(cp810 + 1, NULL, 10) : xs810_1;
      xs810_set = 1;
    }
    if (strcmp(argv[i], "xrhocal") == 0) xrhocal810 = 1;
    /* §812: блики — аналитический Фонг-член в растре (требует кадр со светом
     * и развёртку: дерево ламп строится по eemit). НК-ключи: xgloss0 — ks := 0
     * (кадр посимвольно прежний), xglnorm812 — нормаль перевёрнута у каждого
     * 100-го треугольника (блик-канал обязан поехать). */
    if (strcmp(argv[i], "xgloss") == 0) g_xgloss = 1;
    if (strcmp(argv[i], "xgloss0") == 0) {
      g_xgloss = 1;
      g_glks0 = 1;
    }
    if (strcmp(argv[i], "xglnorm812") == 0) {
      g_xgloss = 1;
      g_glnorm812 = 1;
    }
    if (strcmp(argv[i], "xleak") == 0) xleak = 1;
    if (strcmp(argv[i], "xnopiece") == 0) xnopiece = 1;
    if (strcmp(argv[i], "xobjpiece") == 0) xobjpiece = 1;
    if (strcmp(argv[i], "xclus") == 0) xclus = 1;
    if (strncmp(argv[i], "xclustol=", 9) == 0) xclustol = strtod(argv[i] + 9, NULL);
    if (strcmp(argv[i], "xcontrib") == 0) {
      xcontrib = 1;
      xmatrho = 1;
      doxfer = 1;
      dosolid = 1;
      doxsweep = 1;
    }
    if (strcmp(argv[i], "xframe") == 0) {
      xframe = 1;
      xmatrho = 1;
      doxfer = 1;
      dosolid = 1;
      doxsweep = 1;
    }
    if (strcmp(argv[i], "xdsasplit") == 0) {
      xdsasplit = 1;
      xdsadiff = 1;
      xdsagal = 1;
      xdsa = 1;
      xmatrho = 1;
      doxfer = 1;
      dosolid = 1;
      doxsweep = 1;
    }
    if (strcmp(argv[i], "xcmp") == 0) {
      xcmp = 1;
      xmatrho = 1;
      doxfer = 1;
      dosolid = 1;
      doxsweep = 1;
    }
    if (strcmp(argv[i], "xmatfar") == 0) {
      xmatfar = 1;
      xmatrho = 1;
    }
    if (strcmp(argv[i], "noshift") == 0) g_noshift = 1;
    if (strcmp(argv[i], "raysweep") == 0) g_raysweep = 1;
    if (strncmp(argv[i], "seed=", 5) == 0) g_seed = (int32_t)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "passes=", 7) == 0) g_passes = (int)strtol(argv[i] + 7, NULL, 10);
    if (strncmp(argv[i], "emit=", 5) == 0) g_emitthr = strtod(argv[i] + 5, NULL);
    if (strcmp(argv[i], "treesweep") == 0) g_treesweep = 1;
    if (strncmp(argv[i], "hgather=", 8) == 0) g_hgather = strtod(argv[i] + 8, NULL);
    if (strncmp(argv[i], "hspread=", 8) == 0) g_hspread = strtod(argv[i] + 8, NULL);
    if (strcmp(argv[i], "onebin") == 0) g_onebin = 1;
    if (strncmp(argv[i], "sweepthr=", 9) == 0) {
      g_sweepthr = strtod(argv[i] + 9, NULL);
      g_treesweep = 1;
    }
    /* Р-7а: замер квантования плоскости. */
    if (strcmp(argv[i], "qplane") == 0) {
      qplane = 1;
      doxfer = 1;
      dosolid = 1;
    }
    if (strcmp(argv[i], "xcorner") == 0) xcorner = 1;
    if (strcmp(argv[i], "xinnerfluid") == 0) xinnerfluid = 1;
    if (strcmp(argv[i], "xsweep") == 0) {
      doxfer = 1;
      dosolid = 1;
      doxsweep = 1;
    }
    if (strcmp(argv[i], "xfer") == 0) {
      doxfer = 1;
      dosolid = 1;
    }
    if (strcmp(argv[i], "xfernosolid") == 0) {
      doxfer = 1;
      dosolid = 1;
      xfernosolid = 1;
    }
    if (strcmp(argv[i], "indvis") == 0) indmeas = 1; /* заслоны и так включены (§613) */
    if (strcmp(argv[i], "noindvis") == 0) indvis = 0;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ §397: вершина среза в один бит на ось. */
    if (strcmp(argv[i], "vq1") == 0) vq1 = 1;
    /* ФАЛЬСИФИКАТОР §638: вершина УЗЛА в один бит на ось. Проверяет не точность,
     * а СЛЕПОТУ сверки: срез и дерево читают одну и ту же квантованную вершину,
     * поэтому квант узла в `dmax` входить не обязан. Ключ отличается от `vq1`
     * стороной, которую он огрубляет, и путать их нельзя. */
    if (strcmp(argv[i], "nodevq1") == 0) nodevq1 = 1;
    /* Ш4: удар сферой, два случая врозь (А669). */
    if (strcmp(argv[i], "hit") == 0) hit = 1;
    /* Ф8'-0 (§524): доля листьев свипа, затронутых ударом. `hitr=` — радиус
     * сферы (негативный контроль — `hitr=2.0`, доля обязана уйти в десятки
     * процентов); `hitnoocc` — проверка на ложный ноль, занятость не правится. */
    /* Ф9. (§536): граница узости доли рассеяния — замер без переноса. */
    if (strcmp(argv[i], "lobetest") == 0) {
      g_lobetest = 1;
      lit = 1;
    }
    /* НЕГАТИВНЫЙ КОНТРОЛЬ §536: доля с той же энергией, но БЕЗ зависимости от
     * угла. Отклонение A обязано стать ничтожным при любом s, то есть граница
     * исчезнуть. Не исчезнет — значит меряется не доля, а машинерия. */
    if (strcmp(argv[i], "dlobeflat") == 0) g_dlobeflat = 1;
    /* Ф11. (§545): `ffull` — прежний булев заслон, приёмка сведения;
     * `fzero` — НЕГАТИВНЫЙ КОНТРОЛЬ, поверхности не заслоняют вовсе. */
    if (strcmp(argv[i], "ffull") == 0) g_ffull = 1;
    if (strcmp(argv[i], "fzero") == 0) g_fzero = 1;
    if (strcmp(argv[i], "fnotrans") == 0) g_fnotrans = 1;
    if (strcmp(argv[i], "fnoclamp") == 0) g_fnoclamp = 1;
    if (strcmp(argv[i], "camauto") == 0) g_camauto = 1;
    if (strcmp(argv[i], "caminside") == 0) g_caminside = 1;
    if (strncmp(argv[i], "ceil=", 5) == 0) g_lodceil = strtod(argv[i] + 5, NULL);
    if (strncmp(argv[i], "poly=", 5) == 0) g_lodpoly = strtod(argv[i] + 5, NULL);
    /* §570: СОЛНЦЕ — направленный источник для наружной сцены. */
    if (strcmp(argv[i], "sun") == 0) g_sun = 1;
    if (strcmp(argv[i], "gam16") == 0) g_gamn = 16;
    if (strcmp(argv[i], "texflat") == 0) g_texflat = 1;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ: точечная выборка без пирамиды — рябь обязана
     * вернуться, и расхождение с передискретизованным эталоном вырасти. */
    if (strcmp(argv[i], "texnomip") == 0) g_texnomip = 1;
    if (strcmp(argv[i], "nofrustum") == 0) g_nofrustum = 1;
    if (strcmp(argv[i], "omp1") == 0) g_omp1 = 1;
    /* §585, НК1 и НК2. Оба — ПРИБОР: без первого нечем показать, что выигрыш
     * пришёл от памяти ответа, без второго — что побитовое сличение поломку
     * вообще заметило бы. */
    /* §589: ходьба подразумевает кадр со светом и режим `render` — иначе цикл
     * тащил бы за собой замерную диагностику по нескольку секунд на кадр. */
    if (strcmp(argv[i], "walk") == 0) {
      g_walk = 1;
      g_render = 1;
      lit = 1;
    }
    /* §597: НК1 — прежний путь (сбор в срезе, каждый кадр); НК2 — без подъёма по
     * иерархии; НК3 — подъём невзвешенный. */
    if (strncmp(argv[i], "hbounce=", 8) == 0) g_hbounce = (int)strtol(argv[i] + 8, NULL, 10);
    if (strcmp(argv[i], "keonly") == 0) {
      g_keonly = 1;
      g_render = 1;
      lit = 1;
    }
    if (strcmp(argv[i], "nokemit") == 0) g_nokemit = 1;
    if (strcmp(argv[i], "hangle") == 0) g_hangle = 1;
    if (strncmp(argv[i], "hcontrib=", 9) == 0) {
      g_hcontrib = strtod(argv[i] + 9, NULL);
      hcontrib_set = 1; /* §796: явный порог перебивает fc-умолчание */
    }
    if (strcmp(argv[i], "gflatvis") == 0) g_gflatvis = 1;
    if (strcmp(argv[i], "nolinkcache") == 0) g_nolinkcache = 1;
    if (strcmp(argv[i], "gnonorm") == 0) g_gnonorm = 1;
    if (strcmp(argv[i], "gpoint") == 0) g_gpoint = 1;
    if (strcmp(argv[i], "gwide") == 0) g_gwide = 1;
    if (strcmp(argv[i], "indslice") == 0) g_indslice = 1;
    if (strcmp(argv[i], "indnolift") == 0) g_indnolift = 1;
    if (strcmp(argv[i], "indflat") == 0) g_indflat = 1;
    if (strcmp(argv[i], "nomemo") == 0) hz_dc_walk_memo(HZ_DC_MEMO_OFF);
    if (strcmp(argv[i], "memoscramble") == 0) hz_dc_walk_memo(HZ_DC_MEMO_SCRAMBLE);
    if (strcmp(argv[i], "texflat") == 0) g_texflat = 1;
    if (strncmp(argv[i], "sun=", 4) == 0) {
      const char *sp = argv[i] + 4;
      char *se = NULL;
      double sl = 0.0;
      for (int a = 0; a < 3; a++) {
        g_sundir[a] = strtod(sp, &se);
        sp = (*se == ',') ? se + 1 : se;
        sl += g_sundir[a] * g_sundir[a];
      }
      sl = sqrt(sl);
      if (sl > 0.0)
        for (int a = 0; a < 3; a++)
          g_sundir[a] /= sl;
      g_sun = 1;
    }
    if (strcmp(argv[i], "render") == 0) {
      g_render = 1;
      lit = 1;
    }
    if (strncmp(argv[i], "cam=", 4) == 0) {
      const char *cp = argv[i] + 4;
      char *ep = NULL;
      for (int a = 0; a < 3; a++) {
        g_eye[a] = strtod(cp, &ep);
        cp = (*ep == ',') ? ep + 1 : ep;
      }
      for (int a = 0; a < 3; a++) {
        g_at[a] = strtod(cp, &ep);
        cp = (*ep == ',') ? ep + 1 : ep;
      }
    }
    /* Ф12. (§549): стенд на диффузию правила переноса. */
    if (strncmp(argv[i], "diffbench=", 10) == 0) g_diffbench = (int)strtol(argv[i] + 10, NULL, 10);
    if (strcmp(argv[i], "onenb") == 0) g_onenb = 1;
    /* Ф13. (§553): короткие характеристики; `scharax` — негативный контроль. */
    if (strcmp(argv[i], "schar") == 0) g_schar = 1;
    if (strcmp(argv[i], "cmpcell") == 0) {
      g_cmpcell = 1;
      g_dsweep = 1;
      lit = 1;
    }
    if (strncmp(argv[i], "gstride=", 8) == 0) g_gstride = (int)strtol(argv[i] + 8, NULL, 10);
    if (strcmp(argv[i], "cmpself") == 0) {
      g_cmpself = 1;
      g_cmpcell = 1;
      g_dsweep = 1;
      lit = 1;
    }
    if (strcmp(argv[i], "scharax") == 0) {
      g_schar = 1;
      g_scharax = 1;
    }
    /* НЕГАТИВНЫЙ КОНТРОЛЬ §541: прежняя ТОЧЕЧНАЯ формула. Закон обязан пасть. */
    if (strcmp(argv[i], "ptlight") == 0) g_ptlight = 1;
    if (strncmp(argv[i], "dns=", 4) == 0) g_dns = strtod(argv[i] + 4, NULL);
    if (strncmp(argv[i], "dks=", 4) == 0) g_dks = strtod(argv[i] + 4, NULL);
    /* Ф8. (§532): отскок направленным свипом по дереву. */
    if (strcmp(argv[i], "dsweep") == 0) {
      g_dsweep = 1;
      lit = 1;
    }
    if (strncmp(argv[i], "dnmu=", 5) == 0) {
      g_dnmu = (int)strtol(argv[i] + 5, NULL, 10);
      g_dsweep = 1;
      lit = 1;
    }
    if (strncmp(argv[i], "dnphi=", 6) == 0) {
      g_dnphi = (int)strtol(argv[i] + 6, NULL, 10);
      g_dsweep = 1;
      lit = 1;
    }
    if (strncmp(argv[i], "dpass=", 6) == 0) g_dpass = (int)strtol(argv[i] + 6, NULL, 10);
    /* НЕГАТИВНЫЙ КОНТРОЛЬ §532: излучать во ВСЕ направления, игнорируя нормаль
     * поверхности. Свет пойдёт с изнанки и сквозь тонкую геометрию, и отношение
     * обязано вырасти не менее чем в 1.4 раза. */
    if (strcmp(argv[i], "dblkopen") == 0) {
      g_dblkopen = 1;
      g_dsweep = 1;
      lit = 1;
    }
    if (strcmp(argv[i], "dirsall") == 0) {
      g_dirsall = 1;
      g_dsweep = 1;
      lit = 1;
    }
    /* §529: `albone` — все альбедо единица (фальсификатор разбора). */
    if (strcmp(argv[i], "albone") == 0) g_albone = 1;
    /* §531: `emitalb2` — вернуть ПРЕЖНЕЕ неверное поведение (альбедо излучателя
     * дважды) ради воспроизводимости старых докладов. */
    if (strcmp(argv[i], "emitalb2") == 0) g_emitalb2 = 1;
    if (strcmp(argv[i], "hitsweep") == 0) hitsweep = 1;
    if (strncmp(argv[i], "hitr=", 5) == 0) {
      hitrad = strtod(argv[i] + 5, NULL);
      hitsweep = 1;
    }
    if (strcmp(argv[i], "hitnoocc") == 0) {
      hitnoocc = 1;
      hitsweep = 1;
    }
    if (strcmp(argv[i], "nofix") == 0) {
      hit = 1;
      nofix = 1;
    }
    /* Ш5: кадр со светом. `res=` — разрешение, `thr=` — порог среза в пикселях. */
    if (strcmp(argv[i], "lit") == 0) lit = 1;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ §419: вернуть осевое наследование видимости. */
    /* Ш5а2 включён по умолчанию; ключи возвращают прежние правила. */
    if (strcmp(argv[i], "sweepbool") == 0) sweepfrac = 0;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ Ш5в: отсечение выключено, цена обязана вернуться. */
    if (strcmp(argv[i], "nocull") == 0) {
      lit = 1;
      nocull = 1;
    }
    /* НЕГАТИВНЫЙ КОНТРОЛЬ §430: нулевое альбедо ОТСКОКА — косвенное обязано
     * стать РОВНО нулём, а прямой свет не измениться. */
    if (strcmp(argv[i], "alb0") == 0) {
      lit = 1;
      alb0 = 1;
    }
    /* НЕГАТИВНЫЙ КОНТРОЛЬ §451: рама, совпадающая с габаритом, — площадь
     * коробки обязана вернуться к `1.5000`, записей к `8 649`. */
    if (strcmp(argv[i], "gridalign") == 0) hz_gridalign = 1;
    if (strcmp(argv[i], "flipall") == 0) hz_flipall = 1;
    if (strcmp(argv[i], "longdiag") == 0) hz_longdiag = 1;
    /* ПЛОЩАДЬ ВЫДАННОЙ ПОВЕРХНОСТИ (§446/§450 П3): на любой сцене, без печи. */
    if (strcmp(argv[i], "area") == 0) {
      lit = 1;
      area = 1;
    }
    /* ПЕЧЬ (§430): замкнутая коробка, постоянное альбедо, ответ известен. */
    if (strncmp(argv[i], "oven=", 5) == 0) {
      lit = 1;
      oven = strtod(argv[i] + 5, NULL);
    }
    /* ДВЕ ПЛАСТИНЫ (§450, П6): `E_ind = ρ·E_dir·F`, где `F` — угловой
     * коэффициент двух соосных квадратов, известный в замкнутой форме. Печь
     * ловит СОХРАНЕНИЕ энергии, этот стенд — её РАСПРЕДЕЛЕНИЕ. */
    if (strncmp(argv[i], "plates=", 7) == 0) {
      lit = 1;
      plates = strtod(argv[i] + 7, NULL);
    }
    if (strcmp(argv[i], "sweepr01") == 0) {
      lit = 1;
      sweepr01 = 1;
    }
    if (strcmp(argv[i], "sweepaxis") == 0) {
      sweepfrac = 0;
      lit = 1;
      sweepaxis = 1;
    }
    if (strncmp(argv[i], "res=", 4) == 0) {
      /* `res=WxH` или `res=N` (квадрат). */
      char *e2 = NULL;
      long w2 = strtol(argv[i] + 4, &e2, 10);
      long h2 = (e2 != NULL && (*e2 == 0x78 || *e2 == 0x58)) ? strtol(e2 + 1, NULL, 10) : w2;
      if (w2 > 0 && h2 > 0) {
        resw = (int)w2;
        resh = (int)h2;
      }
    }
    if (strncmp(argv[i], "thr=", 4) == 0) lodthr = strtod(argv[i] + 4, NULL);
    /* Выгрузка занятости для сверки с эталоном Ш0 (`tools/poccref.c`). */
    if (strncmp(argv[i], "occdump=", 8) == 0) occdump = argv[i] + 8;
    /* ЭТАЛОН РЕГРЕССА ДЛЯ Ш2 (А666): выданные многоугольники снимаются здесь. */
    if (strncmp(argv[i], "polydump=", 9) == 0) polydump = argv[i] + 9;
  }
  /* ПОСЛЕ РАЗБОРА КЛЮЧЕЙ, А НЕ ДО (А964). Первая редакция звала самопроверку
   * перед циклом разбора, и негативный контроль `ptlight` до неё не доходил:
   * закон проверял НОВУЮ формулу в обоих прогонах и печатал одинаковые числа.
   * Это третья вакуумная проверка за два дня после А932 и А941 — контроль, до
   * которого не доходит правка, не контроль.
   * И всё же ДО всякого счёта: смысл исполнителя в том, чтобы негодная модель
   * света не доживала до первого замера. */
  g_t0 = now_s();
  alight_selftest();
  emitbudget_selftest();
  if (g_diffbench > 0) {
    diffbench(g_diffbench, 1);
    return 0;
  }
  if (lev < 1 || lev > HZ_DC_MAX_LOG2SIZE) {
    fprintf(stderr, "lev вне разрядного предела (1..%d)\n", HZ_DC_MAX_LOG2SIZE);
    return 2;
  }
  hz_objmesh m;
  if (hz_obj_load(&m, argv[1], strtod(argv[2], NULL)) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  frame fr;
  fr.lev = lev;
  fr.n = (int32_t)1 << lev;
  double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
  for (int32_t k = 0; k < m.nv; k++)
    for (int c = 0; c < 3; c++) {
      double x = m.v[3 * (size_t)k + (size_t)c];
      if (x < lo[c]) lo[c] = x;
      if (x > hi[c]) hi[c] = x;
    }
  /* Кубическая рама: сетка одна на все оси, иначе ключ ребра и координаты DC
   * перестают быть однородными. Габарит расширяется до куба с запасом в ячейку. */
  double side = 0.0;
  for (int c = 0; c < 3; c++)
    if (hi[c] - lo[c] > side) side = hi[c] - lo[c];
  fr.h = side / (double)(fr.n - 2);
  for (int c = 0; c < 3; c++)
    /* ПОЛКЛЕТКИ СМЕЩЕНИЯ — НЕ ПОДГОНКА, А УСЛОВНОСТЬ (§451). Без него
     * `h = сторона/(n−2)` ставит габарит РОВНО на границы ячеек, и всякая осевая
     * стена сцены оказывается в плоскости сеточных вершин. Тогда её пересекают
     * ОБА смежных ребра (`t = 1−δ` снизу и `t = +δ'` сверху, оба строго внутри
     * отрезка), поверхность выдаётся дважды, и площадь выходит в `1.5` раза
     * больше истинной — замерено. Полклетки уводят габарит в СЕРЕДИНУ ячейки,
     * где такого совпадения нет ни у одной из шести граней габарита. Величина
     * `0.5` не настраивается: это центр ячейки, единственная точка, равноудалённая
     * от обеих её границ.
     * ЧЕГО ЭТО НЕ ЛЕЧИТ: внутренняя стена сцены может лечь на плоскость сетки и
     * при сдвинутой раме — случайно, а не систематически. Общее лечение —
     * целочисленный разрез (PLAN_CUT.md, Р-7), и оно сюда не входит. */
    fr.org[c] =
        0.5 * (lo[c] + hi[c]) - 0.5 * (double)fr.n * fr.h + (hz_gridalign ? 0.0 : 0.5 * fr.h);
  printf("== ПОЛЕ (БЕЗЗНАКОВОЕ DC, Р7): %s, треугольников %d, сетка %d^3, ячейка %.4f м%s\n",
         argv[1], m.nt, fr.n, fr.h, nonrm ? "  [НОРМАЛИ ОСЕВЫЕ]" : "");
  /* ГАБАРИТ ПЕЧАТАЕТСЯ ВСЕГДА: без него камеру для новой сцены не назначить, а
   * без камеры сцена не рендерится вовсе (до 08-12 все восемь мест читали
   * координаты комнаты). */
  printf("   ГАБАРИТ: [%.2f %.2f %.2f] … [%.2f %.2f %.2f] м, размах %.2f × %.2f × %.2f\n", lo[0],
         lo[1], lo[2], hi[0], hi[1], hi[2], hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2]);
  if (g_camauto) {
    /* ВЗГЛЯД СНАРУЖИ ГАБАРИТА — годится ПРЕДМЕТУ (шары, дом), но НЕ интерьеру и
     * не городу: там камера обязана стоять внутри, и её задают ключом `cam=`.
     * Правило простое и без подгонки: глаз отодвинут на диагональ габарита от
     * его центра по направлению (−1, +0.6, −1), взгляд — в центр. */
    double c2[3], dg = 0.0;
    for (int a = 0; a < 3; a++) {
      c2[a] = 0.5 * (lo[a] + hi[a]);
      dg += (hi[a] - lo[a]) * (hi[a] - lo[a]);
    }
    dg = sqrt(dg);
    double dv[3] = {-1.0, 0.6, -1.0}, dl2 = sqrt(1.0 + 0.36 + 1.0);
    for (int a = 0; a < 3; a++) {
      g_at[a] = c2[a];
      g_eye[a] = c2[a] + dv[a] / dl2 * dg * 0.75;
    }
  }

  /* ---- 1. занятость ---- */
  double t0 = now_s();
  /* Обнуление ЯВНОЕ: анализатор не видит, что `opyr_build` заполняет все
   * уровни циклом, и это его законное сомнение — массив указателей, частично
   * заполненный, дал бы чтение мусора при первом же изменении порядка. */
  opyr P;
  memset(&P, 0, sizeof P);
  int64_t nocc = opyr_build(&P, &m, &fr);
  double t_occ = now_s() - t0;
  printf(
      "   занятость за %.2f с: занятых ячеек %lld (%.3f %%), пирамида %.1f МБ\n", t_occ,
      (long long)nocc, 100.0 * (double)nocc / (double)((size_t)fr.n * (size_t)fr.n * (size_t)fr.n),
      (double)hz_occ_bytes((size_t)fr.n * (size_t)fr.n * (size_t)fr.n) * (8.0 / 7.0) / 1048576.0);
  /* РЕЖИМ `render` (08-12): диагностические проходы пропускаются целиком. На
   * доме они стоили `> 600` с против `33` с сборки — то есть прибор мешал
   * измерять то, ради чего заведён. Тайминги РАБОЧИХ стадий остаются все. */
  if (!g_render) {
    /* СКОЛЬКО СЕТКИ МОЖНО ПРОПУСТИТЬ: доля ЗАНЯТЫХ ячеек по уровням пирамиды.
     * Если на грубых уровнях занято единицы процентов, то свип, умеющий
     * перепрыгивать пустой узел, платит по ПОВЕРХНОСТИ, а не по объёму. */
    /* СКОЛЬКО УЗЛОВ В АДАПТИВНОМ ПОКРЫТИИ ПРОСТРАНСТВА. Свип по узлам платит по
     * ЭТОМУ числу, а плоский — по числу ячеек сетки. Считается спуском: пустой
     * узел не дробится (сквозь него профиль проносится за один шаг), занятый
     * дробится до предельного уровня. Величина нужна ДО постройки: если узлов
     * немногим меньше, чем ячеек, свип по узлам скорости не даст. */
    for (int lmax = 5; lmax <= 8 && lmax <= lev; lmax++) {
      int64_t nnode = 0;
      int32_t stx[64], sty[64], stz[64], stl[64];
      int sp = 0;
      stx[0] = sty[0] = stz[0] = 0;
      stl[0] = 0;
      sp = 1;
      while (sp > 0) {
        sp--;
        int32_t x = stx[sp], y = sty[sp], z = stz[sp];
        int l = stl[sp];
        nnode++;
        if (l >= lmax) continue;
        if (P.b[l] == NULL) continue;
        int32_t nl = (int32_t)1 << l;
        if (!hz_occ_get(P.b[l], hz_occ_index(nl, x, y, z))) continue;
        for (int k = 0; k < 8 && sp < 60; k++) {
          stx[sp] = 2 * x + (k & 1);
          sty[sp] = 2 * y + ((k >> 1) & 1);
          stz[sp] = 2 * z + ((k >> 2) & 1);
          stl[sp] = l + 1;
          sp++;
        }
      }
      int64_t ncell = (int64_t)1 << (3 * lmax);
      printf("   ПОКРЫТИЕ до уровня L%d: узлов %lld против %lld ячеек плоской сетки "
             "(в %.1f раза меньше)\n",
             lmax, (long long)nnode, (long long)ncell, (double)ncell / (double)(nnode ? nnode : 1));
    }
    printf("   ПИРАМИДА ЗАНЯТОСТИ по уровням (занято/всего):");

    for (int l = 1; l <= lev; l++) {
      /* Проверка на NULL не лишняя: инвариант «opyr_build заполнил все уровни»
       * анализатору не виден, а полагаться на невидимый инвариант в коде,
       * читающем массив указателей, — тот же класс, что чтение мусора. */
      if (P.b[l] == NULL) break;
      int32_t nl = (int32_t)1 << l;
      int64_t tot = (int64_t)nl * nl * nl, occn = 0;
      for (int64_t i = 0; i < tot; i++)
        if (hz_occ_get(P.b[l], (size_t)i)) occn++;
      printf(" L%d %.2f%%", l, 100.0 * (double)occn / (double)tot);
    }
    printf("\n");
  }
  /* ПОКАЗАТЕЛЬ ОГРУБЛЕНИЯ — ГЛАВНЫЙ ДИАГНОЗ ГЕОМЕТРИИ (08-12, замечание
   * пользователя). Во сколько раз падает ЧИСЛО занятых ячеек при подъёме на
   * уровень, говорит, ЧТО в них лежит:
   *     4  — ПОВЕРХНОСТЬ (стена, пол): площадь падает вчетверо на уровень;
   *     8  — ОБЪЁМ (крона дерева, кустарник, любая подъячеечная взвесь):
   *          она заполняет объём, и счёт падает как объём.
   * Восьмёрка есть признак того, что геометрия в этих ячейках НЕ РАЗРЕШАЕТСЯ, а
   * лишь дробится: разрешать её бессмысленно, её надо ГОМОГЕНИЗИРОВАТЬ
   * (`CLAUDE.md`, случай (iii), в коде отсутствует и названо пробелом).
   * ПЕЧАТАЕТСЯ ВСЕГДА, в том числе в режиме `render`: это не диагностика ради
   * диагностики, а величина, по которой решают, где огрублять. */
  {
    printf("   ПОКАЗАТЕЛЬ ОГРУБЛЕНИЯ (во сколько раз падает счёт на уровень; 4 = поверхность, "
           "8 = объём):\n     ");
    int64_t prev = 0;
    for (int l = lev; l >= 1; l--) {
      if (P.b[l] == NULL) break;
      int32_t nl = (int32_t)1 << l;
      int64_t tot = (int64_t)nl * nl * nl, occn = 0;
      for (int64_t i = 0; i < tot; i++)
        if (hz_occ_get(P.b[l], (size_t)i)) occn++;
      if (prev > 0 && occn > 0) printf(" L%d->L%d %.2f", l + 1, l, (double)prev / (double)occn);
      prev = occn;
    }
    printf("\n");
  }
  if (occdump != NULL) {
    int wrc = hz_occ_write(occdump, lev, fr.n, P.b[lev]);
    printf("   ДАМП ЗАНЯТОСТИ -> %s (код %d)\n", occdump, wrc);
  }

  if (g_caminside) {
    /* КАМЕРА ВНУТРИ СЦЕНЫ. Точка и направление ИЩУТСЯ по занятости, а не
     * подбираются координатами: правило обязано работать и на следующей сцене.
     *
     * ЗЕМЛЯ ИЩЕТСЯ В КАЖДОМ СТОЛБЦЕ-КАНДИДАТЕ, А НЕ ОДНА НА СЦЕНУ. Три
     * редакции подряд ошиблись именно датумом высоты, и каждая — по-своему:
     *   от начала РАМЫ  — камера под сценой («до занятой 69.32 м»);
     *   от низа СЦЕНЫ   — на Bistro камера под улицей (y = −3.22 м), потому что
     *                     в габарит входят подвалы и рельеф;
     *   по ЦЕНТРАЛЬНОМУ столбцу — глаз на крыше (y = 20.96 м), потому что в
     *                     центре габарита стоит здание, а не мостовая.
     * Верно так: у КАЖДОГО столбца своя земля (первая занятая ячейка снизу),
     * над нею требуется свободная высота, и уже среди таких выбирается самый
     * открытый. Земля есть то, на что встанешь, — определение, а не эвристика.
     *
     * ЦЕНТРАЛЬНАЯ ПОЛОВИНА плана — чтобы не уйти в чистое поле за сценой: самая
     * открытая точка обычно снаружи, а нужен двор или улица. */
    const double EYEH = 1.6;  /* рост человека — единственное число здесь */
    const double FREEH = 2.5; /* свободная высота над землёй, чтобы не встать под балкой */
    int32_t bx0 = (int32_t)((lo[0] + 0.25 * (hi[0] - lo[0]) - fr.org[0]) / fr.h);
    int32_t bx1 = (int32_t)((lo[0] + 0.75 * (hi[0] - lo[0]) - fr.org[0]) / fr.h);
    int32_t bz0 = (int32_t)((lo[2] + 0.25 * (hi[2] - lo[2]) - fr.org[2]) / fr.h);
    int32_t bz1 = (int32_t)((lo[2] + 0.75 * (hi[2] - lo[2]) - fr.org[2]) / fr.h);
    int32_t y00 = (int32_t)((lo[1] - fr.org[1]) / fr.h);
    int32_t y11 = (int32_t)((hi[1] - fr.org[1]) / fr.h);
    if (y00 < 0) y00 = 0;
    if (y11 >= fr.n) y11 = fr.n - 1;
    int32_t neye = (int32_t)(EYEH / fr.h), nfree = (int32_t)(FREEH / fr.h);
    if (neye < 1) neye = 1;
    if (nfree < 1) nfree = 1;
    /* УРОВЕНЬ ЗЕМЛИ, А НЕ КРЫША. Без этого ограничения «самая открытая точка»
     * встаёт НА КРЫШУ: над нею открыто небо, и по открытости она побеждает
     * всегда — замерено, глаз выходил на `y = 20.96` м при улице около нуля.
     * Берётся нижняя ЧЕТВЕРТЬ распределения высот земли: это и есть уровень
     * улицы или двора. Четверть, а не минимум — минимум пришёлся бы на
     * случайную яму или подвал. */
    int32_t gcap = y11;
    {
      size_t ncol = (size_t)(bx1 - bx0 + 1) * (size_t)(bz1 - bz0 + 1);
      int32_t *gh = malloc((ncol > 0 ? ncol : 1) * sizeof *gh);
      if (gh == NULL) exit(1);
      int64_t ng = 0;
      for (int32_t z = bz0; z <= bz1; z++)
        for (int32_t x = bx0; x <= bx1; x++) {
          if (x < 0 || z < 0 || x >= fr.n || z >= fr.n) continue;
          for (int32_t y = y00; y <= y11; y++)
            if (hz_occ_get(P.b[lev], hz_occ_index(fr.n, x, y, z))) {
              gh[ng++] = y;
              break;
            }
        }
      if (ng > 0) {
        qsort(gh, (size_t)ng, sizeof *gh, cmp_i32);
        gcap = gh[ng / 4] + (int32_t)(1.0 / fr.h);
      }
      free(gh);
    }
    int32_t best[3] = {(bx0 + bx1) / 2, y00 + neye, (bz0 + bz1) / 2};
    double bestd = -1.0;
    int64_t ncand = 0;
    for (int32_t z = bz0; z <= bz1; z++)
      for (int32_t x = bx0; x <= bx1; x++) {
        if (x < 0 || z < 0 || x >= fr.n || z >= fr.n) continue;
        /* Земля этого столбца: первая занятая снизу. */
        int32_t gy = -1;
        for (int32_t y = y00; y <= y11; y++)
          if (hz_occ_get(P.b[lev], hz_occ_index(fr.n, x, y, z))) {
            gy = y;
            break;
          }
        if (gy < 0) continue;    /* столбец пуст насквозь — это не пол, а дыра */
        if (gy > gcap) continue; /* это крыша или балкон, а не уровень улицы */
        int32_t ey = gy + neye;
        if (ey + nfree > y11) continue;
        /* Свободно ли над землёй: иначе камера окажется в толще здания. */
        int freeok = 1;
        for (int32_t y = gy + 1; y <= gy + nfree && freeok; y++)
          if (hz_occ_get(P.b[lev], hz_occ_index(fr.n, x, y, z))) freeok = 0;
        if (!freeok) continue;
        ncand++;
        /* Открытость в горизонтальной плоскости глаза — расширяющимся квадратом
         * с ранним выходом: как только квадрат превысил лучшее, дальше не
         * смотрим. */
        /* ОШИБКА ПЕРВОЙ РЕДАКЦИИ, ПОЙМАННАЯ ЧИСЛОМ: стояло «если радиус не
         * превысил лучшего — бросить», и первый же кандидат ставил `bestd = 1`,
         * после чего ВСЕ следующие обрывались на первом кольце. Замер показал
         * «до ближайшей занятой 0.23 м» при 28 936 кандидатах — то есть
         * побеждал всегда первый. Правильный ранний выход один: остановиться,
         * КОГДА ВСТРЕТИЛИ занятую. Предел кольца задан, чтобы «очень открыто»
         * не стоило квадрата сцены. */
        const int32_t RMAX = 48;
        int32_t rad = 1;
        for (; rad <= RMAX; rad++) {
          int hitocc = 0;
          for (int32_t d = -rad; d <= rad && !hitocc; d++) {
            int32_t qs[4][2] = {
                {x + d, z - rad}, {x + d, z + rad}, {x - rad, z + d}, {x + rad, z + d}};
            for (int q = 0; q < 4 && !hitocc; q++) {
              if (qs[q][0] < 0 || qs[q][1] < 0 || qs[q][0] >= fr.n || qs[q][1] >= fr.n) continue;
              if (hz_occ_get(P.b[lev], hz_occ_index(fr.n, qs[q][0], ey, qs[q][1]))) hitocc = 1;
            }
          }
          if (hitocc) break;
        }
        if ((double)rad > bestd) {
          bestd = (double)rad;
          best[0] = x;
          best[1] = ey;
          best[2] = z;
        }
      }
    g_eye[0] = fr.org[0] + ((double)best[0] + 0.5) * fr.h;
    g_eye[1] = fr.org[1] + ((double)best[1] + 0.5) * fr.h;
    g_eye[2] = fr.org[2] + ((double)best[2] + 0.5) * fr.h;
    /* НАПРАВЛЕНИЕ ВЗГЛЯДА ТОЖЕ ИЩЕТСЯ. Первая редакция смотрела «в дальнюю
     * половину сцены» и упёрлась в стену: кадр вышел одной плоскостью во весь
     * экран. Пробуется 16 азимутов, берётся тот, где свободный ход дальше. */
    double bestlen = -1.0;
    for (int az = 0; az < 16; az++) {
      double ang = 2.0 * 3.14159265358979323846 * (double)az / 16.0;
      double dx = cos(ang), dz = sin(ang), len = 0.0;
      for (double s = fr.h; s < 500.0; s += fr.h) {
        int32_t qx = (int32_t)((g_eye[0] + dx * s - fr.org[0]) / fr.h);
        int32_t qz = (int32_t)((g_eye[2] + dz * s - fr.org[2]) / fr.h);
        if (qx < 0 || qz < 0 || qx >= fr.n || qz >= fr.n) break;
        if (hz_occ_get(P.b[lev], hz_occ_index(fr.n, qx, best[1], qz))) break;
        len = s;
      }
      if (len > bestlen) {
        bestlen = len;
        g_at[0] = g_eye[0] + dx * len;
        g_at[2] = g_eye[2] + dz * len;
      }
    }
    g_at[1] = g_eye[1];
    printf("   КАМЕРА ВНУТРИ: столбцов-кандидатов %lld, глаз на y = %.2f м (земля + %.1f м), "
           "до ближайшей занятой %.2f м, свободный ход взгляда %.2f м\n",
           (long long)ncand, g_eye[1], EYEH, bestd * fr.h, bestlen);
  }
  printf("   КАМЕРА: глаз (%.2f %.2f %.2f) -> (%.2f %.2f %.2f)%s%s\n", g_eye[0], g_eye[1], g_eye[2],
         g_at[0], g_at[1], g_at[2], g_camauto ? "  [camauto]" : "",
         g_caminside ? "  [caminside]" : "");

  /* ---- 2. эрмитовы рёбра: сцена строит и отдаёт ---- */
  hz_htab ht;
  if (hz_htab_init(&ht) != 0) exit(1);
  float *xrad = NULL, **xradlv = NULL; /* Ш15: радианс развёртки на сетке, для картинки */
  unsigned char *solidmask = NULL;     /* Ш12/Ш13: 0 тело, 1 полость, 2 наружное */
  /* ---- 2б. ЗАЛИВКА ВНУТРЕННОСТЕЙ (Ш12, §475) ---- */
  /* ЗАЧЕМ. Развёртке по ординатам нужен вход `solid[ncell]` — «ячейка целиком в
   * материале»; без него свет идёт СКВОЗЬ ТЕЛА (замерено в стыке: объём
   * материала выходил 29.4 вместо 164.6). У беззнакового DC такой величины нет
   * вовсе, и это Р7, а не недоделка: у сцены из односторонних треугольников
   * «внутри» не существует (§391). У ЗАМКНУТОЙ сцены существует, и здесь это
   * проверяется ЧИСЛОМ: знаковый объём меша (теорема о дивергенции) обязан
   * лежать в вилке [ТЕЛО, ТЕЛО + ГРАНИЦА] (А877).
   * Волна идёт от ячейки КАМЕРЫ — точки, заведомо лежащей в полости.
   * ОГРАНИЧЕНИЕ (А878): заливается ТА полость, где камера. На сцене с двумя
   * комнатами вторая будет объявлена сплошным телом, и молча. Названо. */
  if (dosolid) {
    double tsl = now_s();
    size_t nc = (size_t)fr.n * (size_t)fr.n * (size_t)fr.n;
    unsigned char *fl = calloc(nc, 1);
    int32_t *stk = malloc(nc * sizeof *stk > 0 ? (size_t)(1 << 22) * sizeof *stk : 1);
    if (fl == NULL || stk == NULL) exit(1);
    int32_t scap = 1 << 22, ntop = 0;
    double eyec[3] = {g_eye[0], g_eye[1], g_eye[2]};
    int32_t e0[3];
    for (int k = 0; k < 3; k++) {
      double f = floor((eyec[k] - fr.org[k]) / fr.h);
      if (f < 0.0) f = 0.0;
      if (f > (double)(fr.n - 1)) f = (double)(fr.n - 1);
      e0[k] = (int32_t)f;
    }
    int64_t nfluid = 0;
    if (!hz_occ_get(P.b[lev], hz_occ_index(fr.n, e0[0], e0[1], e0[2]))) {
      /* Стек хранит МОРТОНОВ-НЕЗАВИСИМЫЙ линейный индекс; ёмкость растёт. */
      stk[ntop++] = (int32_t)hz_occ_index(fr.n, e0[0], e0[1], e0[2]);
      fl[hz_occ_index(fr.n, e0[0], e0[1], e0[2])] = 1u;
      nfluid = 1;
      while (ntop > 0) {
        size_t ci = (size_t)stk[--ntop];
        /* Раскладка индекса — обратная hz_occ_index: (z·n + y)·n + x. */
        int32_t x = (int32_t)(ci % (size_t)fr.n);
        int32_t y = (int32_t)((ci / (size_t)fr.n) % (size_t)fr.n);
        int32_t z = (int32_t)(ci / ((size_t)fr.n * (size_t)fr.n));
        static const int8_t dxs[6] = {1, -1, 0, 0, 0, 0};
        static const int8_t dys[6] = {0, 0, 1, -1, 0, 0};
        static const int8_t dzs[6] = {0, 0, 0, 0, 1, -1};
        for (int d = 0; d < 6; d++) {
          int32_t nx = x + dxs[d], ny = y + dys[d], nz = z + dzs[d];
          if (nx < 0 || ny < 0 || nz < 0 || nx >= fr.n || ny >= fr.n || nz >= fr.n) continue;
          size_t k2 = hz_occ_index(fr.n, nx, ny, nz);
          if (fl[k2] || hz_occ_get(P.b[lev], k2)) continue;
          fl[k2] = 1u;
          nfluid++;
          if (ntop >= scap) {
            int32_t nsc = scap * 2;
            int32_t *ns = realloc(stk, (size_t)nsc * sizeof *ns);
            if (ns == NULL) exit(1);
            stk = ns;
            scap = nsc;
          }
          stk[ntop++] = (int32_t)k2;
        }
      }
    }
    /* ВТОРАЯ ВОЛНА — ОТ ГРАНИЦЫ КУБА (найдено прогоном, §477). Правило
     * «недостижимо из полости ⇒ тело» неверно: сетка есть КУБ по наибольшему
     * габариту, и всё, что лежит вне сцены (над потолком, за стенами, в запасе
     * по короткой оси), недостижимо из полости, но материалом не является.
     * Замерено до правки: ТЕЛО вышло 491.6 м³ при истинных 36.5.
     * Поэтому классов ЧЕТЫРЕ: полость, оболочка, НАРУЖНОЕ, тело. */
    for (int32_t a = 0; a < fr.n; a++)
      for (int32_t b2 = 0; b2 < fr.n; b2++) {
        static const int8_t f6[6][3] = {{0, 1, 2}, {0, 1, 2}, {1, 0, 2},
                                        {1, 0, 2}, {2, 0, 1}, {2, 1, 0}};
        (void)f6;
        int32_t st[6][3] = {{0, a, b2},        {fr.n - 1, a, b2}, {a, 0, b2},
                            {a, fr.n - 1, b2}, {a, b2, 0},        {a, b2, fr.n - 1}};
        for (int d = 0; d < 6; d++) {
          size_t k2 = hz_occ_index(fr.n, st[d][0], st[d][1], st[d][2]);
          if (fl[k2] || hz_occ_get(P.b[lev], k2)) continue;
          fl[k2] = 2u;
          if (ntop >= scap) {
            int32_t nsc = scap * 2;
            int32_t *ns = realloc(stk, (size_t)nsc * sizeof *ns);
            if (ns == NULL) exit(1);
            stk = ns;
            scap = nsc;
          }
          stk[ntop++] = (int32_t)k2;
        }
      }
    int64_t nout = 0;
    while (ntop > 0) {
      size_t ci = (size_t)stk[--ntop];
      nout++;
      int32_t x = (int32_t)(ci % (size_t)fr.n);
      int32_t y = (int32_t)((ci / (size_t)fr.n) % (size_t)fr.n);
      int32_t z = (int32_t)(ci / ((size_t)fr.n * (size_t)fr.n));
      static const int8_t dxs[6] = {1, -1, 0, 0, 0, 0};
      static const int8_t dys[6] = {0, 0, 1, -1, 0, 0};
      static const int8_t dzs[6] = {0, 0, 0, 0, 1, -1};
      for (int d = 0; d < 6; d++) {
        int32_t nx = x + dxs[d], ny = y + dys[d], nz = z + dzs[d];
        if (nx < 0 || ny < 0 || nz < 0 || nx >= fr.n || ny >= fr.n || nz >= fr.n) continue;
        size_t k2 = hz_occ_index(fr.n, nx, ny, nz);
        if (fl[k2] || hz_occ_get(P.b[lev], k2)) continue;
        fl[k2] = 2u;
        if (ntop >= scap) {
          int32_t nsc = scap * 2;
          int32_t *ns = realloc(stk, (size_t)nsc * sizeof *ns);
          if (ns == NULL) exit(1);
          stk = ns;
          scap = nsc;
        }
        stk[ntop++] = (int32_t)k2;
      }
    }
    int64_t nbnd = 0;
    for (size_t i = 0; i < nc; i++)
      if (hz_occ_get(P.b[lev], i)) nbnd++;
    int64_t nbody = (int64_t)nc - nfluid - nbnd - nout;
    double cv = fr.h * fr.h * fr.h;
    /* ПРОТЕЧКА: флюид на самой границе габарита. Волна пущена изнутри полости,
     * значит любая её ячейка у стенки куба означает дыру в оболочке. */
    int64_t nleak = 0;
    for (int32_t a = 0; a < fr.n; a++)
      for (int32_t b2 = 0; b2 < fr.n; b2++) {
        /* СЧИТАЕТСЯ ТОЛЬКО ВНУТРЕННЯЯ волна (fl == 1): наружная стоит на стенке
         * куба по построению, и путать их значит мерить собственный алгоритм. */
        if (fl[hz_occ_index(fr.n, 0, a, b2)] == 1u) nleak++;
        if (fl[hz_occ_index(fr.n, fr.n - 1, a, b2)] == 1u) nleak++;
        if (fl[hz_occ_index(fr.n, a, 0, b2)] == 1u) nleak++;
        if (fl[hz_occ_index(fr.n, a, fr.n - 1, b2)] == 1u) nleak++;
        if (fl[hz_occ_index(fr.n, a, b2, 0)] == 1u) nleak++;
        if (fl[hz_occ_index(fr.n, a, b2, fr.n - 1)] == 1u) nleak++;
      }
    /* ЭТАЛОН: знаковый объём меша (А868). Считается ЗДЕСЬ, а не однострочником
     * в докладе: у него теперь есть потребитель. Для НЕзамкнутой сцены он
     * бессмыслен (А879), и это сказано рядом с числом. */
    double vsig = 0.0;
    for (int32_t t2 = 0; t2 < m.nt; t2++) {
      const double *A2, *B2, *C2;
      tri_verts(&m, t2, &A2, &B2, &C2);
      vsig += (A2[0] * (B2[1] * C2[2] - B2[2] * C2[1]) - A2[1] * (B2[0] * C2[2] - B2[2] * C2[0]) +
               A2[2] * (B2[0] * C2[1] - B2[1] * C2[0])) /
              6.0;
    }
    printf("   ЗАЛИВКА за %.2f с: ФЛЮИД %lld (%.3f м³), НАРУЖНОЕ %.3f м³, ГРАНИЦА %lld (%.3f м³), "
           "ТЕЛО %lld "
           "(%.3f м³); ВИЛКА [%.3f, %.3f] м³ против знакового объёма меша %.4f м³ — %s; "
           "ПРОТЕЧКА %lld ячеек на стенке куба; маска %.0f МБ\n",
           now_s() - tsl, (long long)nfluid, (double)nfluid * cv, (double)nout * cv,
           (long long)nbnd, (double)nbnd * cv, (long long)nbody, (double)nbody * cv,
           (double)nbody * cv, (double)(nbody + nbnd) * cv, vsig,
           (vsig >= (double)nbody * cv && vsig <= (double)(nbody + nbnd) * cv) ? "В ВИЛКЕ"
                                                                               : "ВНЕ ВИЛКИ",
           (long long)nleak, (double)nc / 1048576.0);
    free(stk);
    solidmask = fl; /* Ш13: маска нужна стыку переноса; освобождается в конце */
  }

  t0 = now_s();
  int64_t nhit = edges_build(&ht, &m, &fr, nonrm);
  int32_t ndrop = 0;
  if (hz_htab_reduce_min(&ht, &ndrop) != HZ_DC_OK) exit(1);
  double t_edges = now_s() - t0;
  printf("   РЁБРА за %.2f с: занесено %lld, УНИКАЛЬНЫХ %d, отброшено дальних %d\n", t_edges,
         (long long)nhit, ht.n, ndrop);

  /* ---- 3. дерево: спуск по занятости, маски рёбер, формы ---- */
  hz_dctree T;
  if (hz_dc_init(&T, lev) != 0) exit(1);
  if (nodevq1) T.vbits = 1; /* ФАЛЬСИФИКАТОР §638: вершина УЗЛА в один бит на ось */
  t0 = now_s();
  int rc = hz_dc_shape_occ(&T, lev, u_occ, &P);
  double t_shape = now_s() - t0;
  if (rc != HZ_DC_OK) {
    fprintf(stderr, "спуск по занятости: код %d\n", rc);
    return 1;
  }
  t0 = now_s();
  hz_dc_masks_occ(&T, &ht);
  double t_masks = now_s() - t0;
  /* НК-1 Ш9: обнулить сложенные суммы — тогда крупная ячейка снова остаётся без
   * нормали, ровно как до правки. Контроль стоит У ПОТРЕБИТЕЛЯ, а не ключом в
   * общем слое: ядро не должно уметь работать неправильно по флагу. */
  if (nonsum && T.nsum != NULL) memset(T.nsum, 0, (size_t)T.nsumcap * sizeof *T.nsum);
  t0 = now_s();
  rc = hz_dc_forms_lazy(&T, &ht);
  double t_dc = now_s() - t0;
  printf("   ДЕРЕВО: спуск %.2f с (узлов %d), маски %.2f с (крупный лист с маской %d — обязан "
         "быть 0), формы %.2f с (код %d; ЗАГНАНО %d)\n",
         t_shape, T.n, t_masks, T.nbigmask, t_dc, rc, T.nclamped);
  /* Г47/§654: отказ выдать вершину при неманифолдной маске УГЛОВ. Печатается С
   * ПОМЕТКОЙ, потому что в БЕЗЗНАКОВОМ пути (`hz_dc_shape_occ`, которым строится
   * всё рабочее) `corner` ставится нулём явно (`dc.c:596`), `hz_dc_manifold(0)`
   * возвращает единицу, и число не может отличаться от нуля ПО ПОСТРОЕНИЮ.
   * Печатать его без этой строки значило бы выдать ложный ноль за измерение. */
  printf("   Г47 ОТКАЗОВ ПО НЕМАНИФОЛДНОЙ МАСКЕ УГЛОВ: %d%s\n", T.nmulti,
         T.unsgn ? "  [БЕЗЗНАКОВЫЙ ПУТЬ: знака углов нет, число есть СТРУКТУРНЫЙ НОЛЬ, "
                   "а не замер — §654]"
                 : "");
  {
    /* САМОПРОВЕРКА ДЕТЕКТОРА ПЕРЕД ЗАМЕРОМ (Р4 §654). Без неё «ноль
     * многопетельных» неотличим от «детектор всегда возвращает единицу».
     * `0x00F0` — плоскость поперёк, режет четыре ребра одной оси: одна петля.
     * `0x999` — срезаны два ПРОТИВОПОЛОЖНЫХ угла, по три ребра у каждого: две. */
    int a1 = 0, a2 = 0;
    int l1 = ecross_loops(0x00F0u, &a1), l2 = ecross_loops(0x999u, &a2);
    printf("   §654 САМОПРОВЕРКА ДЕТЕКТОРА ПЕТЕЛЬ: плоскость -> %d (ждём 1), два угла -> %d "
           "(ждём 2) — %s\n",
           l1, l2, (l1 == 1 && l2 == 2) ? "СОШЛОСЬ" : "ОТКАЗ");
    sheetstat SH;
    memset(&SH, 0, sizeof SH);
    int32_t z0[3] = {0, 0, 0};
    sheet_rec(&T, &ht, 0, z0, fr.n, &SH);
    double den = (double)(SH.npop ? SH.npop : 1);
    printf("   §654 ЛИСТОВ ПОВЕРХНОСТИ В ЯЧЕЙКЕ (НИЖНЯЯ оценка по маске рёбер; популяция — "
           "листья размера 1 с непустой маской, их %lld):\n"
           "      петель 1: %lld (%.3f %%); 2: %lld (%.3f %%); >=3: %lld (%.3f %%); "
           "маска без петель: %lld\n"
           "      неоднозначных граней %lld в %lld ячейках (%.3f %%); ВСТРЕЧНЫЕ НОРМАЛИ у %lld "
           "(%.3f %%)\n",
           (long long)SH.npop, (long long)SH.nloop[1], 100.0 * (double)SH.nloop[1] / den,
           (long long)SH.nloop[2], 100.0 * (double)SH.nloop[2] / den, (long long)SH.nloop[3],
           100.0 * (double)SH.nloop[3] / den, (long long)SH.nloop[0], (long long)SH.nambigf,
           (long long)SH.nambigc, 100.0 * (double)SH.nambigc / den, (long long)SH.nopp,
           100.0 * (double)SH.nopp / den);
  }
  double t_tree = t_shape + t_masks + t_dc;
  int64_t nv = 0;
  int32_t zero[3] = {0, 0, 0};
  int64_t nvl = 0;
  collect_verts(&T, 0, zero, fr.n, NULL, 0, &nv, 0);
  collect_verts(&T, 0, zero, fr.n, NULL, 0, &nvl, 1);
  printf("   ВЕРШИН: У ВСЕХ УЗЛОВ %lld (рёбер/вершину %.3f), У ЛИСТЬЕВ %lld "
         "(рёбер/вершину %.3f)\n",
         (long long)nv, (double)ht.n / (double)(nv ? nv : 1), (long long)nvl,
         (double)ht.n / (double)(nvl ? nvl : 1));
  /* §632: ФОРМЫ СНИМАЮТСЯ ПОСЛЕ ПОСТРОЙКИ. Они весят `112` Б на узел и нужны
   * только сборщику; после решения вершин это мёртвый груз, который иначе
   * таскался бы в каждом обходе. При замере разрушения (`hit`) они остаются:
   * починка без них ОТКАЗЫВАЕТ, и это отказ, а не тишина. */
  double qfmb = (double)T.cap * (double)sizeof(hz_qef) / 1048576.0;
  if (!hit) hz_dc_drop_forms(&T);
  /* СУММЫ НОРМАЛЕЙ ПЕЧАТАЮТСЯ ОТДЕЛЬНОЙ СТРОКОЙ, А НЕ ПРЯЧУТСЯ В «Б/узел»:
   * побочный массив в `sizeof(hz_dcnode)` не входит, и умолчание о нём читалось
   * бы как «правка Ш9 памяти не стоит». */
  printf("   ПАМЯТЬ: дерево DC %.1f МБ (%zu Б/узел), ФОРМЫ %.1f МБ (%s), рёбра %.1f МБ, "
         "пирамида %.1f МБ, СУММЫ НОРМАЛЕЙ %.1f МБ (%.1f Б/узел)\n",
         (double)T.n * (double)sizeof(hz_dcnode) / 1048576.0, sizeof(hz_dcnode), qfmb,
         hit ? "оставлены: нужны починке" : "СНЯТЫ после постройки",
         (double)ht.n * (double)sizeof(hz_hedge) / 1048576.0,
         (double)hz_occ_bytes((size_t)fr.n * (size_t)fr.n * (size_t)fr.n) * (8.0 / 7.0) / 1048576.0,
         (double)T.nsumcap * 24.0 / 1048576.0,
         T.n > 0 ? (double)T.nsumcap * 24.0 / (double)T.n : 0.0);

  /* СПИСОК ТРЕУГОЛЬНИКОВ НА ЯЧЕЙКУ СТРОИТСЯ ЗДЕСЬ, А НЕ ПОСЛЕ СТЫКА (§667).
   * Материал поверхности берётся только отсюда, а стыку с переносом он нужен,
   * чтобы подать `Ke` в источник развёртки. Перенос вверх — правка ПОРЯДКА, не
   * расчёта: ниже по тексту `CT` используется как прежде. */
  celltris CT;
  ct_build(&CT, &m, &fr, P.b[lev], nocc);

  /* ---- 3г. СТЫК С ПЕРЕНОСОМ (Ш13, §479) ---- */
  /* ЧТО ЗДЕСЬ ПРОВЕРЯЕТСЯ. Развёртка по ординатам даёт глобальное освещение С
   * ЗАСЛОНАМИ по построению — ту физику, что маршем стоит 65 с (§474). Её стык
   * с геометрией — `tr3_cut_build`, и принимает он ровно то, что общий слой уже
   * отдаёт: фасеты `hz_dc_facets` плюс маску «ячейка целиком в материале» из
   * заливки Ш12. Здесь развёртка НЕ ЗАПУСКАЕТСЯ: шаг отвечает на один вопрос —
   * принимает ли стык геометрию DC и при каком уровне.
   * ДЕРЕВО СТРОИТСЯ ПО ЗАНЯТОСТИ, А НЕ ПО ОТБОРУ ФАСЕТОВ, как в `render3`: у
   * него фасетов сотни (сфера), у нас их сотни тысяч, и `hz_facets_for_box` на
   * каждую коробку стоил бы O(фасетов). Пирамида занятости отвечает на тот же
   * вопрос за O(1) и по построению согласована с деревом DC. */
  if (doxfer) {
    double tx = now_s();
    hz_facettab ftab;
    hz_cutmap cmap;
    if (hz_facettab_init(&ftab) != 0 || hz_cutmap_init(&cmap) != 0) exit(1);
    /* Р-8/§794: рабочий путь — раздача по КУСКУ (источник кусков: DC-полигоны
     * либо, при xobjpiece, ТРЕУГОЛЬНИКИ СЦЕНЫ — авторские нормали, dmax = 0);
     * xnopiece — прежние бесконечные плоскости (негативный контроль §780). */
    int frc;
    if ((xclus && xobjpiece) || ((xclus || xobjpiece) && xnopiece)) {
      fprintf(stderr, "xclus/xobjpiece/xnopiece взаимоисключающие\n");
      exit(1);
    }
    if (xclus) {
      /* §800: кластеризованные куски-полигоны из треугольников сцены */
      clus800 CL;
      clus800_build(&m, xclustol >= 0.0 ? xclustol : HZ_CLUS_TOL, &CL);
      opoly800 opc = {&m, &fr, &CL};
      frc = hz_dc_facets_polys(&T, opoly_get800, &opc, CL.np, &ftab, &cmap);
      clus800_free(&CL); /* вершины скопированы в ftab при добавлении */
    } else if (xobjpiece) {
      otri794 otc = {&m, &fr};
      frc = hz_dc_facets_tris(&T, otri_get794, &otc, m.nt, &ftab, &cmap);
    } else
      frc = hz_dc_facets2(&T, NULL, NULL, &ftab, &cmap, xnopiece ? 0 : 1);
    printf("   Р-8 РАЗДАЧА (%s): пар %lld, отброшено кусочным отбором %lld, записей опустело "
           "%lld\n",
           xclus ? "КЛАСТЕРЫ §800"
                 : (xobjpiece ? "OBJ-куски §794" : (xnopiece ? "ПЛОСКОСТИ — НК" : "DC-куски")),
           hz_dc_facets_pairs(), hz_dc_facets_dropped(), hz_dc_facets_empty());
    /* Р-7а (§Р-7а): ЗАМЕР КВАНТОВАНИЯ ПЛОСКОСТИ. Вклад в `dmax` считается на
     * ВСЕХ фасетах сцены, радиус — половина диагонали ЕДИНИЧНОЙ ячейки
     * (разрезанная ячейка всегда самого мелкого уровня, условие 1:1). Отдельно
     * считаются ОСЕВЫЕ плоскости: у них нормаль есть орт, и вклад обязан быть
     * РОВНО НОЛЬ — это проверка измерителя, а не измеряемого. */
    if (qplane) {
      hz_frame ofr2 = {{fr.org[0], fr.org[1], fr.org[2]}, {fr.h, fr.h, fr.h}};
      double *qd = malloc((size_t)(ftab.n > 0 ? ftab.n : 1) * sizeof *qd);
      if (qd == NULL) exit(1);
      int64_t nq2 = 0, naxis = 0, naxisbad = 0, nfail = 0;
      double dmaxmax = 0.0;
      for (int32_t i = 0; i < ftab.n; i++) {
        int64_t nqi[3], offqi;
        double add = 0.0;
        if (hz_plane_quant(ftab.f[i].n, ftab.f[i].off, 0.8660254037844386, &ofr2, nqi, &offqi,
                           &add) != 0) {
          nfail++;
          continue;
        }
        int axis = 0;
        for (int k = 0; k < 3; k++)
          if (fabs(fabs(ftab.f[i].n[k]) - 1.0) < 1e-12) axis = 1;
        if (axis) {
          naxis++;
          if (add > 0.0) naxisbad++;
        }
        qd[nq2++] = add;
        if (add > dmaxmax) dmaxmax = add;
      }
      qsort(qd, (size_t)nq2, sizeof *qd, cmp_d);
      printf("      Р-7а КВАНТОВАНИЕ ПЛОСКОСТИ (FRAC %d, NBITS %d): фасетов %lld, вклад в dmax "
             "p50 %.3e p90 %.3e max %.3e м (ячейка %.4f м); ОСЕВЫХ %lld, из них с ненулевым "
             "вкладом %lld; отказов %lld\n",
             HZ_P3_FRAC, HZ_P3_NBITS, (long long)nq2, nq2 ? qd[nq2 / 2] : 0.0,
             nq2 ? qd[(nq2 * 9) / 10] : 0.0, dmaxmax, fr.h, (long long)naxis, (long long)naxisbad,
             (long long)nfail);
      free(qd);
    }
    double t_fac = now_s() - tx;
    printf("   СТЫК: фасетов %d, записей боковой таблицы %d, код %d, за %.2f с\n", ftab.n, cmap.nr,
           frc, t_fac);
    (void)0;
    tx = now_s();
    hz_octree ot;
    if (hz_oct_init(&ot, lev, 0.0) != 0) exit(1);
    oct_by_occ(&ot, &P, lev, 0, 0, 0, fr.n);
    hz_frame ofr = {{fr.org[0], fr.org[1], fr.org[2]}, {fr.h, fr.h, fr.h}};
    tr3_mesh mesh;
    c772 cx772;
    memset(&cx772, 0, sizeof cx772);
    if (xcoarse > 0.0) {
      /* §772: огрубление приёмников — узел размера S листов с дальности D0·S */
      cx772.t = &ot;
      cx772.cm = &cmap;
      cx772.ft = &ftab;
      cx772.fr = &ofr;
      cx772.occ = P.b[lev];
      cx772.smask = solidmask;
      cx772.occn = fr.n;
      cx772.innerfluid = xinnerfluid;
      for (int a = 0; a < 3; a++)
        cx772.eye[a] = g_eye[a];
      cx772.d0 = xcoarse;
      if (tr3_mesh_build_lod(&mesh, &ot, &ofr, c772_stop, &cx772) != 0) exit(1);
      printf("   §772/§782 ОГРУБЛЕНИЕ (D0 = %.1f м/лист): принято %lld узлов (из них пустых "
             "%lld); отказы бюджета %lld; ячеек сетки %d\n",
             xcoarse, (long long)cx772.nacc, (long long)cx772.nacc_empty,
             (long long)cx772.nrej_budget, mesh.ncell);
    } else if (tr3_mesh_build(&mesh, &ot, &ofr) != 0)
      exit(1);
    double t_mesh = now_s() - tx;
    /* МАСКА ПОЛНЫХ ЯЧЕЕК — из заливки. У листа с поверхностью размер 1, у
     * пустого — класс однороден по построению, поэтому хватает ОДНОЙ пробы в
     * его нижнем углу. */
    uint8_t *solid = calloc((size_t)mesh.ncell, 1);
    if (solid == NULL) exit(1);
    int64_t nsolid = 0;
    for (int32_t ci = 0; ci < mesh.ncell; ci++) {
      /* Ш16 (§486): НАРУЖНОЕ тоже сплошное. Снаружи замкнутой комнаты свету
       * взяться неоткуда, а флюид там даёт 458 м³ бесполезной работы и щепки с
       * φ = 1.9e+05. Это утверждение о СЦЕНЕ, а не приближение.
       * §782: у ГРУБОЙ ячейки одной пробы угла МАЛО — узел с поверхностью
       * внутри метился сплошным целиком (замерено: недостача К40 4.6 % на
       * D0=10). Сплошная = ВСЕ листовые клетки сплошные; ранний выход на
       * первом флюиде держит цену (сплошные крупные однородны). */
      int32_t s9 = mesh.csize[ci];
      int allsolid = solidmask != NULL;
      for (int32_t iz = 0; iz < s9 && allsolid; iz++)
        for (int32_t iy = 0; iy < s9 && allsolid; iy++)
          for (int32_t ix = 0; ix < s9 && allsolid; ix++) {
            size_t k = hz_occ_index(fr.n, mesh.clo[ci][0] + ix, mesh.clo[ci][1] + iy,
                                    mesh.clo[ci][2] + iz);
            int cls = solidmask[k];
            if (!((cls == 0u || (cls == 2u && !xinnerfluid)) && !hz_occ_get(P.b[lev], k)))
              allsolid = 0;
          }
      if (allsolid) {
        solid[ci] = 1u;
        nsolid++;
      }
    }
    tx = now_s();
    tr3_cut cut;
    /* §782: листовой каллбэк — только при огрублении (иначе NULL: прежний
     * путь побитово, и НК xcoarse=1e9 остаётся честным) */
    int crc = tr3_cut_build2(&cut, &mesh, &ftab, &cmap, xfernosolid ? NULL : solid,
                             xcoarse > 0.0 ? c772_leaf_solid : NULL, &cx772);
    /* §717: ПЛОЩАДИ ЭЛЕМЕНТОВ И УДАЛЕНИЕ ВЫРОЖДЕННЫХ. Элемент есть сечение
     * коробки плоскостью; у плоскости, чуть задевшей угол, площадь исчезающе
     * мала, а исходящий радианс считается делением на его матрицу масс — тот же
     * механизм, что у ячейки-щепки (§714), только на другом носителе.
     * Удаление идёт ДО подачи `Ke`, чтобы `Σ Ke·площадь` считалась по
     * оставшимся и приёмка «вышло = излучено» осталась верной. */
    if (crc == 0) {
      double h2 = ofr.u[0] * ofr.u[1];
      double *ar = malloc((size_t)(cut.nse > 0 ? cut.nse : 1) * sizeof *ar);
      if (ar != NULL) {
        int64_t na = 0;
        double atot0 = 0.0;
        for (int32_t k = 0; k < cut.nse; k++) {
          if (!(cut.se[k].area > 0.0)) continue;
          ar[na++] = cut.se[k].area / h2;
          atot0 += cut.se[k].area;
        }
        if (na > 0) {
          qsort(ar, (size_t)na, sizeof *ar, cmp_dev699);
          int64_t nb1 = 0, nb2 = 0, nb3 = 0;
          for (int64_t i = 0; i < na; i++) {
            if (ar[i] < 1e-6) nb1++;
            if (ar[i] < 1e-4) nb2++;
            if (ar[i] < 1e-2) nb3++;
          }
          printf("   §717 ПЛОЩАДЬ ЭЛЕМЕНТОВ (в долях h²): элементов %lld; мин %.3e, p1 %.3e, "
                 "медиана %.3e, макс %.3e; ниже 1e-6: %lld (%.2f %%), ниже 1e-4: %lld (%.2f %%), "
                 "ниже 1e-2: %lld (%.2f %%)\n",
                 (long long)na, ar[0], ar[na / 100], ar[na / 2], ar[na - 1], (long long)nb1,
                 100.0 * (double)nb1 / (double)na, (long long)nb2, 100.0 * (double)nb2 / (double)na,
                 (long long)nb3, 100.0 * (double)nb3 / (double)na);
        }
        free(ar);
        /* §778: ПОКРЫТИЕ ПОВЕРХНОСТИ ЭЛЕМЕНТАМИ — дыры покрытия есть канал
         * протечки (§767). Эталон Σ tri.area из OBJ; Σ se.area — ДО фильтра
         * xsemin: фильтр в разложении утечки идёт отдельной осью. */
        {
          double atri778 = 0.0;
          for (int32_t t9 = 0; t9 < m.nt; t9++) {
            double u[3], w[3], cx9[3];
            const double *va = m.v + 3 * (size_t)m.f[3 * (size_t)t9];
            const double *vb = m.v + 3 * (size_t)m.f[3 * (size_t)t9 + 1];
            const double *vc = m.v + 3 * (size_t)m.f[3 * (size_t)t9 + 2];
            for (int k = 0; k < 3; k++) {
              u[k] = vb[k] - va[k];
              w[k] = vc[k] - va[k];
            }
            cx9[0] = u[1] * w[2] - u[2] * w[1];
            cx9[1] = u[2] * w[0] - u[0] * w[2];
            cx9[2] = u[0] * w[1] - u[1] * w[0];
            atri778 += 0.5 * sqrt(cx9[0] * cx9[0] + cx9[1] * cx9[1] + cx9[2] * cx9[2]);
          }
          int64_t nz778 = 0;
          for (int32_t k = 0; k < cut.nse; k++)
            if (!(cut.se[k].area > 0.0)) nz778++;
          printf("   §778 ПОКРЫТИЕ: Σ tri.area %.6f, Σ se.area %.6f, ДЕФИЦИТ %.4f %%; "
                 "элементов с area <= 0: %lld из %d\n",
                 atri778, atot0, 100.0 * (1.0 - atot0 / (atri778 > 0.0 ? atri778 : 1.0)),
                 (long long)nz778, cut.nse);
          if (xleak) {
            /* §778 Р2: эталон по клеткам листовой решётки независимым клипом;
             * невязка (ref − se) с обоими хвостами: дыры И избыток (А1239 —
             * дубли прятали бы дефицит, их выдаёт отрицательный хвост). */
            size_t ng9 = (size_t)fr.n * (size_t)fr.n * (size_t)fr.n;
            double *refA = calloc(ng9, sizeof *refA);
            double *seA = calloc(ng9, sizeof *seA);
            if (refA == NULL || seA == NULL) exit(1);
            for (int32_t t9 = 0; t9 < m.nt; t9++) {
              double tri[3][3];
              for (int q = 0; q < 3; q++)
                for (int k = 0; k < 3; k++)
                  tri[q][k] = m.v[3 * (size_t)m.f[3 * (size_t)t9 + (size_t)q] + (size_t)k];
              int lo9[3], hi9[3];
              for (int a = 0; a < 3; a++) {
                double bl = tri[0][a], bh = tri[0][a];
                for (int q = 1; q < 3; q++) {
                  if (tri[q][a] < bl) bl = tri[q][a];
                  if (tri[q][a] > bh) bh = tri[q][a];
                }
                lo9[a] = (int)floor((bl - ofr.o[a]) / ofr.u[a]);
                hi9[a] = (int)floor((bh - ofr.o[a]) / ofr.u[a]);
                if (lo9[a] < 0) lo9[a] = 0;
                if (hi9[a] > fr.n - 1) hi9[a] = fr.n - 1;
              }
              for (int iz = lo9[2]; iz <= hi9[2]; iz++)
                for (int iy = lo9[1]; iy <= hi9[1]; iy++)
                  for (int ix = lo9[0]; ix <= hi9[0]; ix++) {
                    double blo[3] = {ofr.o[0] + ofr.u[0] * (double)ix,
                                     ofr.o[1] + ofr.u[1] * (double)iy,
                                     ofr.o[2] + ofr.u[2] * (double)iz};
                    double bhi[3] = {blo[0] + ofr.u[0], blo[1] + ofr.u[1], blo[2] + ofr.u[2]};
                    double a9 = leak_tri_box_area((const double (*)[3])tri, blo, bhi);
                    if (a9 > 0.0) refA[hz_occ_index(fr.n, ix, iy, iz)] += a9;
                  }
            }
            int64_t nbig9 = 0;
            for (int32_t k = 0; k < cut.nse; k++) {
              if (!(cut.se[k].area > 0.0)) continue;
              int32_t c9 = cut.se[k].cell;
              if (mesh.csize[c9] != 1) {
                nbig9++; /* элементов в КРУПНЫХ ячейках быть не должно — счёт */
                continue;
              }
              seA[hz_occ_index(fr.n, mesh.clo[c9][0], mesh.clo[c9][1], mesh.clo[c9][2])] +=
                  cut.se[k].area;
            }
            /* допуск невязки клетки: плавучка двух клипов на площади ~h² */
            const double eps9 = 1e-9;
            double reftot = 0.0, refcut = 0.0, secut = 0.0, possum = 0.0, negsum = 0.0;
            int64_t ncut9 = 0, npos = 0, nneg9 = 0;
            enum { NT9 = 10 };
            double tdev[NT9];
            size_t tg[NT9];
            int tn = 0;
            double ndev[3] = {0, 0, 0};
            size_t ngc[3] = {0, 0, 0};
            int nn9 = 0;
            for (size_t g = 0; g < ng9; g++) {
              reftot += refA[g];
              if (!(refA[g] > 0.0) && !(seA[g] > 0.0)) continue;
              ncut9++;
              refcut += refA[g];
              secut += seA[g];
              double dv = refA[g] - seA[g];
              if (dv > eps9) {
                possum += dv;
                npos++;
                int j = tn < NT9 ? tn : NT9 - 1;
                if (tn < NT9 || dv > tdev[NT9 - 1]) {
                  while (j > 0 && tdev[j - 1] < dv) {
                    tdev[j] = tdev[j - 1];
                    tg[j] = tg[j - 1];
                    j--;
                  }
                  tdev[j] = dv;
                  tg[j] = g;
                  if (tn < NT9) tn++;
                }
              } else if (dv < -eps9) {
                negsum -= dv;
                nneg9++;
                int j = nn9 < 3 ? nn9 : 2;
                if (nn9 < 3 || -dv > ndev[2]) {
                  while (j > 0 && ndev[j - 1] < -dv) {
                    ndev[j] = ndev[j - 1];
                    ngc[j] = ngc[j - 1];
                    j--;
                  }
                  ndev[j] = -dv;
                  ngc[j] = g;
                  if (nn9 < 3) nn9++;
                }
              }
            }
            printf("   §778 КЛИП: Σ ref %.6f (самопроверка против Σ tri.area: Δ %.2e отн.); "
                   "клеток с поверхностью %lld; ДЫРЫ Σ(ref−se)+ %.6f (%.3f %% от ref) в %lld "
                   "клетках; ИЗБЫТОК Σ(se−ref)+ %.6f в %lld клетках; se в крупных ячейках %lld\n",
                   reftot, fabs(reftot - atri778) / (atri778 > 0.0 ? atri778 : 1.0),
                   (long long)ncut9, possum, 100.0 * possum / (refcut > 0.0 ? refcut : 1.0),
                   (long long)npos, negsum, (long long)nneg9, (long long)nbig9);
            /* Гистограмма дыр по ЧИСЛУ ФАСЕТОВ записи и по КЛАССУ реза —
             * атрибуция классом, а не десятью топ-клетками (А11-класс). */
            {
              int32_t *cellat = malloc(ng9 * sizeof *cellat);
              if (cellat != NULL) {
                for (size_t g = 0; g < ng9; g++)
                  cellat[g] = -1;
                for (int32_t c9 = 0; c9 < mesh.ncell; c9++)
                  if (mesh.csize[c9] == 1)
                    cellat[hz_occ_index(fr.n, mesh.clo[c9][0], mesh.clo[c9][1], mesh.clo[c9][2])] =
                        c9;
                double hnf[4] = {0, 0, 0, 0}; /* nf: 0-1 / 2-3 / 4-7 / >=8 */
                /* рез: пустой / полный / частичный / нет ячейки / СПЛОШНАЯ —
                 * у сплошной ref есть, а элементов нет НАМЕРЕННО (внутренность
                 * тела); её «дыра» ложная и из суммы механизма исключается */
                double hfr[5] = {0, 0, 0, 0, 0};
                double vb9 = ofr.u[0] * ofr.u[1] * ofr.u[2];
                /* допуск класса реза: 1e-9 объёма коробки — тот же класс, что
                 * у объёмной сверки §772 (макроскопика против плавучки) */
                const double veps = 1e-9 * vb9;
                for (size_t g = 0; g < ng9; g++) {
                  double dv = refA[g] - seA[g];
                  if (!(dv > eps9)) continue;
                  int32_t cm9 = cellat[g];
                  if (cm9 < 0) {
                    hfr[3] += dv;
                    continue;
                  }
                  if (cut.solid[cm9]) {
                    hfr[4] += dv;
                    continue;
                  }
                  const hz_cutrec *r9 = hz_cutmap_find(&cmap, mesh.node[cm9]);
                  int nf9 = r9 != NULL ? (int)r9->nf : 0;
                  hnf[nf9 <= 1 ? 0 : (nf9 <= 3 ? 1 : (nf9 <= 7 ? 2 : 3))] += dv;
                  double vf9 = cut.mvol[cm9][0][0];
                  hfr[vf9 < veps ? 0 : (vf9 > vb9 - veps ? 1 : 2)] += dv;
                }
                printf("   §778 ДЫРЫ ПО ФАСЕТАМ ЗАПИСИ (без сплошных): 0-1: %.4f, 2-3: %.4f, "
                       "4-7: %.4f, >=8: %.4f; ПО КЛАССУ РЕЗА: пустой %.4f, полный %.4f, "
                       "частичный %.4f, без ячейки %.4f, СПЛОШНАЯ (ложная дыра) %.4f\n",
                       hnf[0], hnf[1], hnf[2], hnf[3], hfr[0], hfr[1], hfr[2], hfr[3], hfr[4]);
                free(cellat);
              }
            }
            for (int i = 0; i < tn; i++) {
              size_t g = tg[i];
              int32_t ix = (int32_t)(g % (size_t)fr.n);
              int32_t iy = (int32_t)((g / (size_t)fr.n) % (size_t)fr.n);
              int32_t iz = (int32_t)(g / ((size_t)fr.n * (size_t)fr.n));
              int32_t cm9 = -1;
              for (int32_t c9 = 0; c9 < mesh.ncell; c9++)
                if (mesh.csize[c9] == 1 && mesh.clo[c9][0] == ix && mesh.clo[c9][1] == iy &&
                    mesh.clo[c9][2] == iz) {
                  cm9 = c9;
                  break;
                }
              int64_t nse9 = 0;
              if (cm9 >= 0)
                for (int32_t k = 0; k < cut.nse; k++)
                  if (cut.se[k].cell == cm9 && cut.se[k].area > 0.0) nse9++;
              int nf9 = -1;
              double ffr9 = -1.0;
              if (cm9 >= 0) {
                const hz_cutrec *r9 = hz_cutmap_find(&cmap, mesh.node[cm9]);
                nf9 = r9 != NULL ? (int)r9->nf : 0;
                double vb9 = ofr.u[0] * ofr.u[1] * ofr.u[2];
                ffr9 = vb9 > 0.0 ? cut.mvol[cm9][0][0] / vb9 : -1.0;
              }
              printf("      §778 дыра %2d: клетка (%d,%d,%d) мир (%.3f, %.3f, %.3f); потеря "
                     "%.6f (ref %.6f, se %.6f); элементов %lld, фасетов записи %d, доля флюида "
                     "%.3f; ячейка сетки %s\n",
                     i, ix, iy, iz, ofr.o[0] + ofr.u[0] * ((double)ix + 0.5),
                     ofr.o[1] + ofr.u[1] * ((double)iy + 0.5),
                     ofr.o[2] + ofr.u[2] * ((double)iz + 0.5), tdev[i], refA[g], seA[g],
                     (long long)nse9, nf9, ffr9, cm9 >= 0 ? "есть" : "НЕТ");
            }
            for (int i = 0; i < nn9; i++) {
              size_t g = ngc[i];
              printf("      §778 избыток %d: клетка g=%zu, se−ref %.6f (ref %.6f, se %.6f)\n", i, g,
                     ndev[i], refA[g], seA[g]);
            }
            free(refA);
            free(seA);
          }
        }
        if (xsemin > 0.0) {
          int64_t nrem = 0;
          double arem = 0.0;
          for (int32_t k = 0; k < cut.nse; k++) {
            if (!(cut.se[k].area > 0.0)) continue;
            if (cut.se[k].area >= xsemin * h2) continue;
            arem += cut.se[k].area;
            cut.se[k].area = 0.0;
            memset(cut.se[k].m, 0, sizeof cut.se[k].m);
            nrem++;
          }
          printf("   §717 УДАЛЕНО ВЫРОЖДЕННЫХ ЭЛЕМЕНТОВ (порог %.3e h²): %lld из %d; потеряно "
                 "площади %.4e из %.4e (%.4f %%)\n",
                 xsemin, (long long)nrem, cut.nse, arem, atot0,
                 100.0 * arem / (atot0 > 0.0 ? atot0 : 1.0));
        }
      }
    }
    /* §711: ЩЕПКИ УХОДЯТ ИЗ НОСИТЕЛЯ. Ячейка, где флюида меньше `xthin` доли,
     * объявляется СПЛОШНОЙ: она и так почти целиком материал, а деление на её
     * исчезающий объём даёт радианс на порядки выше соседского (§710). Её грани
     * к флюидным соседям становятся стыком, у которого с §707 есть граничное
     * условие, — энергия не теряется, а отражается.
     * Порог задаётся КЛЮЧОМ и просматривается свипом: назначать его числом в
     * коде значило бы повторить ошибку §690. */
    if (xthin > 0.0 && crc == 0) {
      int64_t nconv = 0;
      double vlost = 0.0, vall = 0.0;
      for (int32_t ci = 0; ci < mesh.ncell; ci++) {
        double vfl = cut.mvol[ci][0][0];
        if (!(vfl > 0.0)) continue;
        vall += vfl;
      }
      for (int32_t ci = 0; ci < mesh.ncell; ci++) {
        double vfl = cut.mvol[ci][0][0];
        if (!(vfl > 0.0) || cut.solid[ci]) continue;
        double s3 = (double)mesh.csize[ci];
        double V = s3 * s3 * s3 * ofr.u[0] * ofr.u[1] * ofr.u[2];
        if (!(V > 0.0) || vfl >= xthin * V) continue;
        cut.solid[ci] = 1;
        vlost += vfl;
        nconv++;
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
            cut.mvol[ci][i][j] = 0.0;
      }
      printf("   §711 ЩЕПКИ ИЗ НОСИТЕЛЯ (порог %.3e): переведено в сплошные %lld ячеек; "
             "потеряно флюидного объёма %.4e из %.4e (%.4f %%)\n",
             xthin, (long long)nconv, vlost, vall, 100.0 * vlost / (vall > 0.0 ? vall : 1.0));
    }
    /* §702: СЧЁТ СТЫКА «ФЛЮИД — СПЛОШНОЕ». Грани флюидных ячеек раскладываются
     * на три класса; сумма обязана дать полное число, и это печатается как
     * проверка разбиения. Класс «стык» делится надвое: есть ли у флюидной ячейки
     * хотя бы один поверхностный элемент. Ячейка без элементов — место, где
     * энергия уходит в материал и нигде не учитывается (`sweep3.c:562`
     * обнуляет радианс сплошной ячейки). */
    if (crc == 0) {
      int64_t nff = 0, nfs = 0, nbd = 0, nfs_bare = 0, ncell_bare = 0, ncell_fl = 0;
      double aff = 0.0, afs = 0.0, abd = 0.0, afs_bare = 0.0, asel = 0.0;
      for (int32_t ci = 0; ci < mesh.ncell; ci++) {
        if (cut.solid[ci]) continue;
        if (!(cut.mvol[ci][0][0] > 0.0)) continue;
        ncell_fl++;
        int has_se = (cut.sestart[ci + 1] > cut.sestart[ci]);
        if (!has_se) ncell_bare++;
        for (int32_t q = mesh.fstart[ci]; q < mesh.fstart[ci + 1]; q++) {
          int32_t fi = mesh.flist[q];
          int32_t nb = (mesh.f[fi].ca == ci) ? mesh.f[fi].cb : mesh.f[fi].ca;
          double ar = cut.farea[fi];
          if (nb < 0) {
            nbd++;
            abd += ar;
          } else if (cut.solid[nb]) {
            nfs++;
            afs += ar;
            if (!has_se) {
              nfs_bare++;
              afs_bare += ar;
            }
          } else {
            nff++;
            aff += ar;
          }
        }
      }
      for (int32_t k = 0; k < cut.nse; k++)
        asel += cut.se[k].area;
      double atot = aff + afs + abd;
      printf("   §702 ГРАНИ ФЛЮИДНЫХ ЯЧЕЕК (%lld ячеек, из них БЕЗ элементов %lld = %.2f %%):\n"
             "      флюид-флюид  %8lld граней, площадь %.4e (%.2f %%)\n"
             "      ФЛЮИД-СПЛОШНОЕ %6lld граней, площадь %.4e (%.2f %%); ИЗ НИХ В ЯЧЕЙКАХ БЕЗ "
             "ЭЛЕМЕНТОВ %lld, площадь %.4e (%.2f %% от всей границы флюида)\n"
             "      граница области %5lld граней, площадь %.4e (%.2f %%)\n"
             "      СУММА КЛАССОВ %lld против полного числа граней флюидных ячеек — обязана "
             "совпасть\n"
             "      площадь поверхностных ЭЛЕМЕНТОВ %.4e; отношение к площади стыка %.4f\n",
             (long long)ncell_fl, (long long)ncell_bare,
             100.0 * (double)ncell_bare / (double)(ncell_fl ? ncell_fl : 1), (long long)nff, aff,
             100.0 * aff / (atot > 0.0 ? atot : 1.0), (long long)nfs, afs,
             100.0 * afs / (atot > 0.0 ? atot : 1.0), (long long)nfs_bare, afs_bare,
             100.0 * afs_bare / (atot > 0.0 ? atot : 1.0), (long long)nbd, abd,
             100.0 * abd / (atot > 0.0 ? atot : 1.0), (long long)(nff + nfs + nbd), asel,
             afs > 0.0 ? asel / afs : -1.0);
    }
    /* §697: СИНТЕТИЧЕСКИЙ КОНТРОЛЬ. Флюидная матрица масс заменяется матрицей
     * ЦЕЛОЙ ячейки — диагональю `diag(V, V/12, V/12, V/12)` в базисе
     * `{1, ξ, η, ζ}`. Прогон НЕФИЗИЧЕН и служит РАЗЛИЧИТЕЛЕМ: он отвечает, сидит
     * ли усиление в объёме разреза. Его числа в доклады о физике не идут. */
    if (xwholemass && crc == 0) {
      for (int32_t ci = 0; ci < mesh.ncell; ci++) {
        double s3 = (double)mesh.csize[ci];
        double V = s3 * s3 * s3 * ofr.u[0] * ofr.u[1] * ofr.u[2];
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
            cut.mvol[ci][i][j] = 0.0;
        cut.mvol[ci][0][0] = V;
        for (int i = 1; i < 4; i++)
          cut.mvol[ci][i][i] = V / 12.0;
      }
      printf("   §697 ПОДМЕНА: флюидная матрица масс заменена ЦЕЛОЙ у всех %d ячеек "
             "(ПРОГОН НЕФИЗИЧЕН, только различитель)\n",
             mesh.ncell);
    }
    double t_cut = now_s() - tx;
    double vfl = 0.0;
    if (crc == 0)
      for (int32_t ci = 0; ci < mesh.ncell; ci++)
        vfl += cut.mvol[ci][0][0];
    double vbox = 0.0;
    for (int32_t ci = 0; ci < mesh.ncell; ci++) {
      double s = (double)mesh.csize[ci] * fr.h;
      vbox += s * s * s;
    }
    printf("   СТЫК: ячеек сетки %d (полных в материале %lld), граней %d; сетка за %.2f с, "
           "разрез за %.2f с, код %d\n",
           mesh.ncell, (long long)nsolid, mesh.nf, t_mesh, t_cut, crc);
    if (crc == 0)
      printf("      nbad %d (нарушивших 1:1), nse %d (поверхностных элементов), nsebig %d "
             "(многоугольник не поместился), РЕЗ ОТКАЗАН (веер > %d) у %d клеток; ОБЪЁМ: куб "
             "%.3f, ФЛЮИД %.3f, материал %.3f м³%s\n",
             cut.nbad, cut.nse, cut.nsebig, HZ_P3_MAXH, cut.nrezover, vbox, vfl, vbox - vfl,
             xfernosolid ? "   [БЕЗ МАСКИ — НК]" : "");
    /* ---- РАЗВЁРТКА ПО ОРДИНАТАМ (Ш14, §482) ---- */
    /* ИСТОЧНИК — САМ ПОТОЛОК, А НЕ ОТДЕЛЬНОЕ ТЕЛО. Лампа в `pfield` есть
     * площадка под потолком, геометрии у неё нет; развёртка же светит
     * ПОВЕРХНОСТЯМИ. Поэтому светящимися объявляются те фасеты, чьи
     * поверхностные элементы лежат в прямоугольнике лампы — физически это
     * люминесцентная панель заподлицо с потолком, и новой геометрии не надо.
     * СТЕНКИ КУБА ЧЁРНЫЕ: они лежат ВНЕ комнаты, свет до них не доходит, а
     * ненулевое отражение на них добавило бы энергию из ниоткуда. */
    if (crc == 0 && doxsweep) {
      double *frho = calloc((size_t)ftab.n, sizeof *frho);
      double *femit = calloc((size_t)ftab.n, sizeof *femit);
      /* §670: излучение ПО ЭЛЕМЕНТУ, в РАДИАНСЕ. Фасетный массив остаётся ради
       * ключа `xemitfacet` (негативный контроль) и ветви `xhall`. */
      double *eemit = calloc((size_t)(cut.nse > 0 ? cut.nse : 1), sizeof *eemit);
      double *sig_t = calloc((size_t)mesh.ncell, sizeof *sig_t);
      double *sig_s = calloc((size_t)mesh.ncell, sizeof *sig_s);
      double *phi = calloc((size_t)mesh.ncell * 4, sizeof *phi);
      if (eemit == NULL) exit(1);
      if (frho == NULL || femit == NULL || sig_t == NULL || sig_s == NULL || phi == NULL) exit(1);
      for (int32_t i = 0; i < ftab.n; i++)
        frho[i] = xrho;
      /* §737: накопители материального альбедо по фасету (только под xmatrho) */
      double *frnum = NULL, *frden = NULL, *frmn = NULL, *frmx = NULL;
      if (xmatrho) {
        if (xrho_set || xhall) {
          fprintf(stderr, "xmatrho несовместим с явным xrho= и с xhall\n");
          exit(1);
        }
        frnum = calloc((size_t)ftab.n, sizeof *frnum);
        frden = calloc((size_t)ftab.n, sizeof *frden);
        frmn = malloc((size_t)ftab.n * sizeof *frmn);
        frmx = malloc((size_t)ftab.n * sizeof *frmx);
        if (frnum == NULL || frden == NULL || frmn == NULL || frmx == NULL) exit(1);
        for (int32_t i = 0; i < ftab.n; i++) {
          frmn[i] = 2.0;
          frmx[i] = -1.0;
        }
      }
      double lc[3], lu = 0.0, lv = 0.0;
      for (int k = 0; k < 3; k++)
        lc[k] = 0.5 * (lo[k] + hi[k]);
      /* ПОТОЛОК ПОЛОСТИ — тем же спуском от КАМЕРЫ вверх, что и у площадки в
       * ветви `lit` (§472): «пусто» в занятости значит «нет поверхности», а не
       * «нет материала», поэтому считать сверху нельзя. */
      double eyew[3] = {g_eye[0], g_eye[1], g_eye[2]};
      int32_t ex = (int32_t)((eyew[0] - fr.org[0]) / fr.h);
      int32_t ez = (int32_t)((eyew[2] - fr.org[2]) / fr.h);
      int32_t ey = (int32_t)((eyew[1] - fr.org[1]) / fr.h);
      if (ex < 0) ex = 0;
      if (ez < 0) ez = 0;
      if (ey < 0) ey = 0;
      if (ex >= fr.n) ex = fr.n - 1;
      if (ez >= fr.n) ez = fr.n - 1;
      if (ey >= fr.n) ey = fr.n - 1;
      double ceilY = hi[1];
      for (int32_t yy = ey; yy < fr.n; yy++)
        if (hz_occ_get(P.b[lev], hz_occ_index(fr.n, ex, yy, ez)) != 0) {
          ceilY = fr.org[1] + ((double)yy + 0.5) * fr.h;
          break;
        }
      lc[1] = ceilY;
      lu = 0.25 * (hi[0] - lo[0]);
      lv = 0.25 * (hi[2] - lo[2]);
      int64_t nlit = 0, nconf = 0, nbigcell = 0, nnomat = 0;
      /* §739: статистика смены kd при атрибуции по близости против ls[0] */
      int64_t nkchg = 0, nkall = 0;
      double dksum = 0.0, dkmax = 0.0;
      double emitpow = 0.0;
      /* §784: раскладка ФАКТИЧЕСКИ ПРИМЕНЁННОЙ эмиссии по материалам — прибор
       * односторонней сверки (А1249) и адреса потерь огрубления (А1257) */
      double *emat784 = calloc((size_t)(m.nmtl > 0 ? m.nmtl : 1), sizeof *emat784);
      if (emat784 == NULL) exit(1);
      if (xhall) {
        /* ПРЕЖНИЙ ИСТОЧНИК (§479): площадка под потолком, найденная спуском ОТ
         * КАМЕРЫ. Оставлен ключом, чтобы прежние прогоны воспроизводились, но
         * умолчанием быть не может: он камерозависим (§597) и на сцене без
         * потолка над камерой даёт РОВНО НОЛЬ светящихся фасетов (§665). */
        for (int32_t k = 0; k < cut.nse; k++) {
          if (cut.se[k].nv <= 0) continue;
          double c[3] = {0, 0, 0};
          for (int q = 0; q < cut.se[k].nv; q++)
            for (int a = 0; a < 3; a++)
              c[a] += cut.se[k].v[q][a] / (double)cut.se[k].nv;
          if (fabs(c[0] - lc[0]) <= lu && fabs(c[2] - lc[2]) <= lv && c[1] >= lc[1] - 2.0 * fr.h &&
              c[1] <= lc[1] + 2.0 * fr.h) {
            femit[cut.se[k].facet] = 1.0;
            nlit++;
          }
        }
      } else {
        /* ИСТОЧНИК ИЗ `Ke` (§667, ровно как §618 сделал для сбора). Путь:
         * элемент -> его ячейка -> треугольники ячейки -> материал -> `ke3`.
         * Свип монохромный, поэтому берётся среднее трёх каналов; единицы —
         * излучённая РАДИОСИТЬ без множителя `π` (`scene_obj.h`).
         *
         * ФАСЕТ ОБЩИЙ У МНОГИХ ЭЛЕМЕНТОВ, А МАТЕРИАЛ ЖИВЁТ У ЭЛЕМЕНТА (§667 Р3):
         * одна плоскость может нести и вывеску, и тёмную стену. Берётся
         * МАКСИМУМ, а число расхождений СЧИТАЕТСЯ и печатается — приближение
         * обязано быть числом, а не умолчанием. */
        /* §782: буфер объединения ct-списков для ГРУБОЙ ячейки — треугольники
         * собираются по листовым клеткам BBOX КУСКА элемента, а не по одному
         * угловому листу (А1217: с угла терялась Ke-эмиссия — замерено −31 %
         * эмиссии на D0=10). Ёмкость: bbox куска ~ листовой масштаб. */
        enum { LS782 = 4096 };
        int32_t ls782[LS782];
        int64_t nlsclip782 = 0;
        for (int32_t k = 0; k < cut.nse; k++) {
          if (cut.se[k].nv <= 0) continue;
          int32_t ci = cut.se[k].cell;
          if (ci < 0 || ci >= mesh.ncell) continue;
          const int32_t *ls = NULL;
          int32_t nls = 0;
          if (mesh.csize[ci] > 1 && cut.se[k].facet >= 0 && cut.se[k].facet < ftab.n &&
              ftab.f[cut.se[k].facet].bounded) {
            nbigcell++;
            const hz_facet *fp8 = &ftab.f[cut.se[k].facet];
            int32_t blo8[3], bhi8[3];
            /* §800: bbox куска — по всем tnv вершинам полигона */
            int ntv8 = fp8->tnv >= 3 && fp8->tnv <= HZ_FACET_TVMAX ? fp8->tnv : 3;
            for (int a = 0; a < 3; a++) {
              double mn = fp8->tv[0][a], mx = fp8->tv[0][a];
              for (int q2 = 1; q2 < ntv8; q2++) {
                if (fp8->tv[q2][a] < mn) mn = fp8->tv[q2][a];
                if (fp8->tv[q2][a] > mx) mx = fp8->tv[q2][a];
              }
              /* кламп — к РЕШЁТКЕ, а не к коробке ячейки: кусок у края ячейки
               * держит треугольники-источники в соседней клетке, и кламп к
               * коробке терял их материал/Ke (часть −8.3 % эмиссии D0=10) */
              blo8[a] = (int32_t)floor(mn);
              bhi8[a] = (int32_t)floor(mx);
              if (blo8[a] < 0) blo8[a] = 0;
              if (bhi8[a] > fr.n - 1) bhi8[a] = fr.n - 1;
            }
            int32_t nu = 0;
            for (int32_t iz = blo8[2]; iz <= bhi8[2]; iz++)
              for (int32_t iy = blo8[1]; iy <= bhi8[1]; iy++)
                for (int32_t ix = blo8[0]; ix <= bhi8[0]; ix++) {
                  int32_t cl8[3] = {ix, iy, iz};
                  const int32_t *l2 = NULL;
                  int32_t n2 = ct_list(&CT, cl8, &l2);
                  for (int32_t q2 = 0; q2 < n2; q2++) {
                    if (nu >= LS782) {
                      nlsclip782++;
                      break;
                    }
                    ls782[nu++] = l2[q2];
                  }
                }
            ls = ls782;
            nls = nu;
          } else {
            int32_t cellc[3];
            for (int a = 0; a < 3; a++)
              cellc[a] = mesh.clo[ci][a];
            nls = ct_list(&CT, cellc, &ls);
          }
          if (nls == 0) {
            nnomat++;
            continue;
          }
          int32_t mi = m.fm != NULL ? m.fm[ls[0]] : 0;
          if (mi < 0 || mi >= m.nmtl) mi = 0;
          /* §737: материальное альбедо, площадно-взвешенно ПО ФАСЕТУ. Путь
           * материала тот же, что у Ke; kd — одноканальное среднее
           * (scene_obj.h: «перенос считает им»), свип монохромный. Разброс kd
           * на общем фасете СЧИТАЕТСЯ (А1135 — класс смазки §670). */
          if (xmatrho && frnum != NULL) {
            /* §739: материал — треугольник, БЛИЖАЙШИЙ К ПЛОСКОСТИ элемента:
             * score = dist(центроид, плоскость треугольника) + h·(1−|n_t·n_e|);
             * при равных (допуск 1e-6·h — плоскостная точность) ближний по
             * центроиду треугольника (А1141: dist слеп к коллинеарным). Ke-путь
             * НЕ трогается — точка отсчёта эмиссии §670. `xmatfar` (НК) берёт
             * НАИХУДШИЙ score — мажоранта произвола атрибуции. */
            double ec[3] = {0, 0, 0};
            for (int q2 = 0; q2 < cut.se[k].nv; q2++)
              for (int a = 0; a < 3; a++)
                ec[a] += cut.se[k].v[q2][a] / (double)cut.se[k].nv;
            const double *en = cut.se[k].n;
            int32_t tbest = ls[0];
            double sbest = xmatfar ? -1.0 : 1e300, cbest = 1e300;
            for (int32_t q2 = 0; q2 < nls; q2++) {
              const double *A3, *B3, *C3;
              tri_verts(&m, ls[q2], &A3, &B3, &C3);
              double e1[3], e2[3], nt[3];
              for (int a = 0; a < 3; a++) {
                e1[a] = B3[a] - A3[a];
                e2[a] = C3[a] - A3[a];
              }
              nt[0] = e1[1] * e2[2] - e1[2] * e2[1];
              nt[1] = e1[2] * e2[0] - e1[0] * e2[2];
              nt[2] = e1[0] * e2[1] - e1[1] * e2[0];
              double nl = sqrt(nt[0] * nt[0] + nt[1] * nt[1] + nt[2] * nt[2]);
              if (!(nl > 0.0)) continue;
              double dist = 0.0, dot = 0.0, cd = 0.0;
              for (int a = 0; a < 3; a++) {
                dist += nt[a] / nl * (ec[a] - A3[a]);
                dot += nt[a] / nl * en[a];
                double tc = (A3[a] + B3[a] + C3[a]) / 3.0 - ec[a];
                cd += tc * tc;
              }
              double sc = fabs(dist) + fr.h * (1.0 - fabs(dot));
              if (xmatfar) {
                if (sc > sbest) {
                  sbest = sc;
                  tbest = ls[q2];
                }
              } else if (sc < sbest - 1e-6 * fr.h || (sc < sbest + 1e-6 * fr.h && cd < cbest)) {
                sbest = sc;
                cbest = cd;
                tbest = ls[q2];
              }
            }
            int32_t mi2 = m.fm != NULL ? m.fm[tbest] : 0;
            if (mi2 < 0 || mi2 >= m.nmtl) mi2 = 0;
            double kd = m.mtl[mi2].kd;
            if (!(kd >= 0.0 && kd <= 1.0)) kd = 0.5; /* как в alb() сбора */
            double kd0 = m.mtl[mi].kd;
            if (!(kd0 >= 0.0 && kd0 <= 1.0)) kd0 = 0.5;
            if (fabs(kd - kd0) > 0.0) {
              nkchg++;
              dksum += fabs(kd - kd0);
              if (fabs(kd - kd0) > dkmax) dkmax = fabs(kd - kd0);
            }
            nkall++;
            int32_t fi2 = cut.se[k].facet;
            double ar = cut.se[k].area;
            frnum[fi2] += kd * ar;
            frden[fi2] += ar;
            if (kd < frmn[fi2]) frmn[fi2] = kd;
            if (kd > frmx[fi2]) frmx[fi2] = kd;
          }
          (void)mi; /* материал ls[0] нужен только статистике смены kd */
        }
        /* §786: ПЛОЩАДНАЯ РАЗДАЧА Ke (А1261). Атрибуция «весь элемент от
         * одного треугольника» давала нули у субвоксельных ламп и переборы
         * ×18…×21 (замерено §784/§785). Здесь энергия ЯЧЕЙКИ считается
         * ТОЧНЫМ клипом треугольников к листовым клеткам (машинерия §778) и
         * раздаётся элементам ячейки по площадям: Σ мощностей = энергии
         * ячейки точно. СЕМАНТИКА ДВУСТОРОННОСТИ — СОХРАНЕНИЕ ЭНЕРГИИ:
         * DC-лист двуслойной шторы несёт сумму обеих сторон (§786).
         * Энергия клеток без элементов НЕ раздаётся и печатается (плафон в
         * сплошном — потеря по-имени, не молчание). §670-единица прежняя:
         * Ke — радиосность, радианс = E/(π·Σ area). */
        {
          double tke0 = now_s();
          int32_t *cellat6 = malloc((size_t)fr.n * (size_t)fr.n * (size_t)fr.n * sizeof *cellat6);
          double *Ecell6 = calloc((size_t)mesh.ncell, sizeof *Ecell6);
          double *Sarea6 = calloc((size_t)mesh.ncell, sizeof *Sarea6);
          if (cellat6 == NULL || Ecell6 == NULL || Sarea6 == NULL) exit(1);
          for (size_t g = 0; g < (size_t)fr.n * (size_t)fr.n * (size_t)fr.n; g++)
            cellat6[g] = -1;
          for (int32_t ci = 0; ci < mesh.ncell; ci++)
            for (int32_t iz = 0; iz < mesh.csize[ci]; iz++)
              for (int32_t iy = 0; iy < mesh.csize[ci]; iy++)
                for (int32_t ix = 0; ix < mesh.csize[ci]; ix++)
                  cellat6[hz_occ_index(fr.n, mesh.clo[ci][0] + ix, mesh.clo[ci][1] + iy,
                                       mesh.clo[ci][2] + iz)] = ci;
          for (int32_t k = 0; k < cut.nse; k++)
            if (cut.se[k].area > 0.0) Sarea6[cut.se[k].cell] += cut.se[k].area;
          double elost_ne = 0.0, elost_nc = 0.0, etri_tot = 0.0;
          for (int32_t t9 = 0; t9 < m.nt; t9++) {
            int32_t mt = m.fm != NULL ? m.fm[t9] : 0;
            if (mt < 0 || mt >= m.nmtl) mt = 0;
            const double *ke3 = m.mtl[mt].ke3;
            double ke = (ke3[0] + ke3[1] + ke3[2]) / 3.0;
            if (!(ke > 0.0)) continue;
            double tri[3][3];
            for (int q2 = 0; q2 < 3; q2++)
              for (int a = 0; a < 3; a++)
                tri[q2][a] = m.v[3 * (size_t)m.f[3 * (size_t)t9 + (size_t)q2] + (size_t)a];
            int32_t blo9[3], bhi9[3];
            for (int a = 0; a < 3; a++) {
              double mn = tri[0][a], mx = tri[0][a];
              for (int q2 = 1; q2 < 3; q2++) {
                if (tri[q2][a] < mn) mn = tri[q2][a];
                if (tri[q2][a] > mx) mx = tri[q2][a];
              }
              blo9[a] = (int32_t)floor((mn - ofr.o[a]) / ofr.u[a]);
              bhi9[a] = (int32_t)floor((mx - ofr.o[a]) / ofr.u[a]);
              if (blo9[a] < 0) blo9[a] = 0;
              if (bhi9[a] > fr.n - 1) bhi9[a] = fr.n - 1;
            }
            for (int32_t iz = blo9[2]; iz <= bhi9[2]; iz++)
              for (int32_t iy = blo9[1]; iy <= bhi9[1]; iy++)
                for (int32_t ix = blo9[0]; ix <= bhi9[0]; ix++) {
                  double blo[3] = {ofr.o[0] + ofr.u[0] * (double)ix,
                                   ofr.o[1] + ofr.u[1] * (double)iy,
                                   ofr.o[2] + ofr.u[2] * (double)iz};
                  double bhi[3] = {blo[0] + ofr.u[0], blo[1] + ofr.u[1], blo[2] + ofr.u[2]};
                  double A = leak_tri_box_area((const double (*)[3])tri, blo, bhi);
                  if (!(A > 0.0)) continue;
                  etri_tot += ke * A;
                  int32_t ci = cellat6[hz_occ_index(fr.n, ix, iy, iz)];
                  if (ci < 0) {
                    elost_nc += ke * A;
                  } else if (!(Sarea6[ci] > 0.0)) {
                    elost_ne += ke * A;
                  } else {
                    Ecell6[ci] += ke * A;
                    emat784[mt] += ke * A;
                  }
                }
          }
          for (int32_t k = 0; k < cut.nse; k++) {
            int32_t ci = cut.se[k].cell;
            if (!(cut.se[k].area > 0.0) || !(Ecell6[ci] > 0.0)) continue;
            eemit[k] = Ecell6[ci] / (3.14159265358979323846 * Sarea6[ci]);
            emitpow += Ecell6[ci] * cut.se[k].area / Sarea6[ci];
            double b6 = 3.14159265358979323846 * eemit[k];
            if (b6 > femit[cut.se[k].facet]) femit[cut.se[k].facet] = b6;
            nlit++;
          }
          printf("   §786 РАЗДАЧА Ke: Σ по треугольникам %.4f = роздано %.4f + потеряно (клетки "
                 "без элементов) %.4f + (вне ячеек) %.4f; %.2f с\n",
                 etri_tot, etri_tot - elost_ne - elost_nc, elost_ne, elost_nc, now_s() - tke0);
          free(cellat6);
          free(Ecell6);
          free(Sarea6);
        }
        /* §784: печать раскладки — доля печатается от СВОЕЙ суммы прибора;
         * сверка с OBJ-аналитикой (awk) идёт по ИМЕНИ (А1260). */
        double etot784 = 0.0;
        for (int32_t mi9 = 0; mi9 < m.nmtl; mi9++)
          etot784 += emat784[mi9];
        printf("   §784 ЭМИССИЯ ПО МАТЕРИАЛАМ (Σ %.4f; тождество против Σ Ke·площадь — "
               "УБИВАЕТ-проверка):\n",
               etot784);
        for (int32_t mi9 = 0; mi9 < m.nmtl; mi9++) {
          double keavg = (m.mtl[mi9].ke3[0] + m.mtl[mi9].ke3[1] + m.mtl[mi9].ke3[2]) / 3.0;
          if (!(keavg > 0.0) && !(emat784[mi9] > 0.0)) continue;
          printf("      §784 %-40s Ke %.3f  Σ Ke·area %.4f\n", m.mtl[mi9].name, keavg,
                 emat784[mi9]);
        }
      }
      free(emat784);
      /* §737: заполнение и печать материального альбедо. Без материала — 0.5,
       * как в alb() сбора; сводка обязана быть напечатана, иначе сверка с
       * эталоном останется без знаменателя (А1063). */
      if (xmatrho && frnum != NULL) {
        int64_t nsup = 0, nnomat2 = 0, nconfr = 0;
        double wnum = 0.0, wden = 0.0, spread = 0.0;
        double *med9 = malloc((size_t)(ftab.n > 0 ? ftab.n : 1) * sizeof *med9);
        if (med9 == NULL) exit(1);
        for (int32_t i = 0; i < ftab.n; i++) {
          if (frden[i] > 0.0) {
            frho[i] = frnum[i] / frden[i];
            med9[nsup++] = frho[i];
            wnum += frnum[i];
            wden += frden[i];
            if (frmx[i] - frmn[i] > 0.0) {
              nconfr++;
              if (frmx[i] - frmn[i] > spread) spread = frmx[i] - frmn[i];
            }
          } else {
            frho[i] = 0.5;
            nnomat2++;
          }
        }
        double wavg = wden > 0.0 ? wnum / wden : -1.0;
        double medv = -1.0;
        if (nsup > 0) {
          qsort(med9, (size_t)nsup, sizeof *med9, cmp_dev699);
          medv = med9[nsup / 2];
        }
        printf("   §737 МАТЕРИАЛЬНЫЕ АЛЬБЕДО: фасетов с носителем %lld из %d, без материала %lld "
               "(-> 0.5); площадь-взвешенное %.4f, медиана %.4f; фасетов с РАЗНЫМИ kd %lld "
               "(макс разброс %.3f)\n",
               (long long)nsup, ftab.n, (long long)nnomat2, wavg, medv, (long long)nconfr, spread);
        printf("   §739 АТРИБУЦИЯ (%s): элементов со сменой kd против ls[0]: %lld из %lld "
               "(%.2f %%), средний |dkd| %.4f, макс %.4f\n",
               xmatfar ? "ДАЛЬНИЙ, НК" : "БЛИЖНИЙ", (long long)nkchg, (long long)nkall,
               100.0 * (double)nkchg / (double)(nkall > 0 ? nkall : 1),
               nkchg > 0 ? dksum / (double)nkchg : 0.0, dkmax);
        free(med9);
      }
      free(frnum);
      free(frden);
      free(frmn);
      free(frmx);
      /* §742: масштаб всех альбедо — извлечение ряда отскоков: доля(s) =
       * s·B1 + s²·B2 + …; при s = 0 psin обязан совпасть с xrho=0 побитово.
       * Умножение безусловное: умолчание 1.0 ничего не меняет, а сравнение
       * плавучих на равенство запрещено гейтом. */
      for (int32_t i = 0; i < ftab.n; i++)
        frho[i] *= xrhoscale;
      /* ВХОД РАЗВЁРТКИ, А НЕ ТОЛЬКО ЕЁ ВЫХОД (А807). Ложный ноль на выходе
       * неотличим от «источник не задан», пока не напечатан сам источник. */
      /* ВТОРАЯ СУММА СЧИТАЕТСЯ ПО САМОМУ МАССИВУ (§670 Р4), а не повторным
       * умножением тех же чисел: иначе она совпала бы тождественно и не могла
       * бы поймать ошибку ЗАПИСИ в `eemit`. */
      double emitpow2 = 0.0;
      for (int32_t k = 0; k < cut.nse; k++)
        emitpow2 += 3.14159265358979323846 * eemit[k] * cut.se[k].area;
      int64_t nlitfac = 0, nlitfluid = 0;
      for (int32_t i = 0; i < ftab.n; i++)
        if (femit[i] > 0.0) nlitfac++;
      for (int32_t k = 0; k < cut.nse; k++)
        if (femit[cut.se[k].facet] > 0.0 && cut.mvol[cut.se[k].cell][0][0] > 0.0) nlitfluid++;
      printf("      ВХОД РАЗВЁРТКИ: светящихся ФАСЕТОВ %lld из %d (%.3f %%), светящихся элементов "
             "во ФЛЮИДНЫХ ячейках %lld из %lld; альбедо фасета %.2f, ординат %d; ИСТОЧНИК %s\n",
             (long long)nlitfac, ftab.n, 100.0 * (double)nlitfac / (double)(ftab.n ? ftab.n : 1),
             (long long)nlitfluid, (long long)nlit, frho[0], 2 * nmu * 4 * nmu,
             xhall ? "площадка от камеры (§479)" : "Ke материалов (§667)");
      if (!xhall)
        printf("      ИСТОЧНИК ИЗ Ke: мощность Σ Ke·площадь %.4e, расхождений материала на общем "
               "фасете %lld, элементов в КРУПНЫХ ячейках %lld, без материала %lld; НОСИТЕЛЬ %s, "
               "Σ π·L·площадь %.4e\n",
               emitpow, (long long)nconf, (long long)nbigcell, (long long)nnomat,
               xemitfacet ? "ФАСЕТ (НК §670)" : "ЭЛЕМЕНТ (§670)", emitpow2);
      /* §699: ПЕЧЬ НА РАЗРЕЗАННОЙ ГЕОМЕТРИИ. Форма задачи взята у `t_furnace` из
       * `tests/test_sweep3.c`, где она проходит на коробках: альбедо РОВНО 1,
       * отражение 1, влёт постоянный. Тогда постоянный радианс есть ТОЧНОЕ
       * решение, и ответ `φ = 4π·binc0` известен без всякой эталонной
       * реализации. Меняется ровно одно — сцена. */
      double bconst = 0.0;
      if (xconst) {
        bconst = 1.0;
        for (int32_t ci = 0; ci < mesh.ncell; ci++) {
          sig_t[ci] = 0.8;
          sig_s[ci] = 0.8;
        }
        for (int32_t i = 0; i < ftab.n; i++)
          frho[i] = 1.0;
        for (int32_t k = 0; k < cut.nse; k++)
          eemit[k] = 0.0;
        printf("   §699 ПЕЧЬ НА РАЗРЕЗЕ: альбедо 1 (σ_t = σ_s = 0.8), отражение 1, влёт %.3f; "
               "точное решение φ = 4π·влёт = %.6f\n",
               bconst, 4.0 * 3.14159265358979323846 * bconst);
      }
      tr3_dirs dirs;
      if (tr3_dirs_product(&dirs, nmu, nmu) != 0) exit(1);
      /* ---- §796: FIRST-COLLISION SOURCE — инъекция первым отражением ----
       * Прямой свет от eemit-ламп собирается непрерывным механизмом (efc_*),
       * elem_emit подменяется на ρ·E_fc/hsum; Ke в свип не инъецируется.
       * Пути со СТАРОЙ семантикой инъекции сохраняют её и под ключом:
       * xhall/xemitfacet (elem_emit там NULL), xconst (печь, eemit = 0),
       * xunit (источники в ноль), xcmp (мост к ядру считает прямой свипом),
       * xcontrib (перевозмущение ρ→0 обязано гасить и первый отскок),
       * xdsa (закрытая линия §757 со старыми единицами моста). */
      double *efc796 = NULL, *einj796 = NULL;
      double injpow796 = 0.0;
      int fcold = !xfc || xhall || xemitfacet || xconst || xunit || xcmp || xcontrib || xdsa;
      if (!fcold) {
        double tfc0 = now_s();
        int32_t nsrc = 0;
        for (int32_t k = 0; k < cut.nse; k++)
          if (eemit[k] > 0.0 && cut.se[k].area > 0.0) nsrc++;
        if (nsrc == 0) {
          printf("   §796 FIRST-COLLISION: светящихся элементов нет — инъекция прежняя\n");
          fcold = 1;
        } else {
          efcsrc *sr6 = malloc((size_t)nsrc * sizeof *sr6);
          efc796 = calloc((size_t)(cut.nse > 0 ? cut.nse : 1), sizeof *efc796);
          einj796 = calloc((size_t)(cut.nse > 0 ? cut.nse : 1), sizeof *einj796);
          if (sr6 == NULL || efc796 == NULL || einj796 == NULL) exit(1);
          int32_t js = 0;
          for (int32_t k = 0; k < cut.nse; k++) {
            if (!(eemit[k] > 0.0 && cut.se[k].area > 0.0)) continue;
            efcsrc *s = &sr6[js];
            efc_center(&cut, &mesh, &fr, k, s->c);
            for (int a = 0; a < 3; a++)
              s->n[a] = cut.se[k].n[a];
            s->area = cut.se[k].area;
            s->flux = 3.14159265358979323846 * eemit[k] * cut.se[k].area;
            s->idx = js;
            s->key = 0;
            for (int sh = lev - 1; sh >= 0; sh--) {
              uint64_t bit = (uint64_t)(((mesh.clo[cut.se[k].cell][0] >> sh) & 1) |
                                        (((mesh.clo[cut.se[k].cell][1] >> sh) & 1) << 1) |
                                        (((mesh.clo[cut.se[k].cell][2] >> sh) & 1) << 2));
              s->key = (s->key << 3) | bit;
            }
            js++;
          }
          qsort(sr6, (size_t)nsrc, sizeof *sr6, efc_cmp);
          etree ET6;
          memset(&ET6, 0, sizeof ET6);
          efc_build(&ET6, sr6, 0, nsrc, 0, lev);
          int64_t nl6 = 0, nm6 = 0;
          double all6 = 0.0, thr6 = 0.0;
          efc_run(&ET6, &cut, &mesh, &fr, &P, hcontrib_set ? g_hcontrib : HZ_FC_CONTRIB, efc796,
                  &nl6, &all6, &thr6, &nm6);
          etree_free(&ET6);
          free(sr6);
          /* перехват и инъекция; сторож самозаслона А1292 — доля E_fc = 0 */
          double icept6 = 0.0;
          int64_t nz6 = 0, nrc6 = 0;
          for (int32_t k = 0; k < cut.nse; k++) {
            if (!(cut.se[k].area > 0.0)) continue;
            nrc6++;
            icept6 += efc796[k] * cut.se[k].area;
            if (!(efc796[k] > 0.0)) nz6++;
            double rho6 =
                (cut.se[k].facet >= 0 && cut.se[k].facet < ftab.n) ? frho[cut.se[k].facet] : 0.0;
            /* К29-нормировка Σ_{ω·n>0} w(ω·n) считается ЗДЕСЬ по набору
             * ординат (поля se.hsum в cut3 НЕТ С §804 — оно не заполнялось
             * нигде, А1294; поймано первым же прогоном: инъекция через него
             * выходила нулевой). Деление
             * на неё, а не на π, делает инъецированную мощность В ЕДИНИЦАХ
             * СХЕМЫ равной ρ·E_fc·area точно. */
            double hs6 = 0.0;
            for (int mm = 0; mm < dirs.n; mm++) {
              double on6 = dirs.ox[mm] * cut.se[k].n[0] + dirs.oy[mm] * cut.se[k].n[1] +
                           dirs.oz[mm] * cut.se[k].n[2];
              if (on6 > 0.0) hs6 += dirs.w[mm] * on6;
            }
            if (hs6 > 0.0) einj796[k] = rho6 * efc796[k] / hs6;
            injpow796 += rho6 * efc796[k] * cut.se[k].area;
          }
          printf("   §796 FIRST-COLLISION: излучателей %d, излучено %.6e, ПЕРЕХВАЧЕНО %.6e "
                 "(%.1f %%), инъекция Σρ·E·area %.6e; связей %lld (заслонённого потока %.1f %%), "
                 "маршей %lld; E_fc = 0 у %lld из %lld элементов (%.1f %%); %.2f с\n",
                 nsrc, emitpow2, icept6, 100.0 * icept6 / (emitpow2 > 0.0 ? emitpow2 : 1.0),
                 injpow796, (long long)nl6, 100.0 * thr6 / (all6 > 0.0 ? all6 : 1.0),
                 (long long)nm6, (long long)nz6, (long long)nrc6,
                 100.0 * (double)nz6 / (double)(nrc6 > 0 ? nrc6 : 1), now_s() - tfc0);
        }
      }
      tr3_problem prob = {.m = &mesh,
                          .d = &dirs,
                          .cut = &cut,
                          .facet_rho = frho,
                          .facet_emit = femit,
                          .elem_emit = (xhall || xemitfacet) ? NULL : (fcold ? eemit : einj796),
                          .nfacet = ftab.n,
                          .solid_rho = xsolidrho,
                          .binc0 = bconst,
                          .sig_t = sig_t,
                          .sig_s = sig_s,
                          .limiter = xnolim ? 0 : 1,
                          .maxp_off = xnomaxp,
                          .relax = xrelax,
                          .chain_cell1 = xchain >= 0 ? xchain + 1 : 0,
                          .trace = xtrace,
                          .reltol = xreltol};
      if (xchain >= 0 && !xunit) {
        fprintf(stderr, "xchain= требует xunit: пол обрыва трассы — уровень единичного входа\n");
        exit(1);
      }
      /* §744: сличение — рабочий режим, не диагностический */
      if (xdsa && (xunit || xconst || xcmp)) {
        fprintf(stderr, "xdsa несовместим с xunit/xconst/xcmp\n");
        exit(1);
      }
      if (xcmp && (xunit || xconst)) {
        fprintf(stderr, "xcmp несовместим с xunit и xconst\n");
        exit(1);
      }
      /* §774: лестница отскоков — рабочий режим поверх обычной задачи; с
       * приборными и двухсеточными режимами не совмещается сознательно. */
      if (xbounce > 0 && (xdsa || xunit || xconst || xcmp || xcontrib)) {
        fprintf(stderr, "xbounce= несовместим с xdsa/xunit/xconst/xcmp/xcontrib\n");
        exit(1);
      }
      /* §733: вскрытие обновления — предсказания калиброваны единичным входом */
      if ((xcelll[0] >= 0 || xdir >= 0) && (!xunit || xcelll[0] < 0 || xdir < 0)) {
        fprintf(stderr, "xcell=/xdir= требуют xunit и друг друга\n");
        exit(1);
      }
      for (int q = 0; q < 4; q++)
        prob.dump_cell1[q] = xcelll[q] >= 0 ? xcelll[q] + 1 : 0;
      prob.dump_dir1 = xdir >= 0 ? xdir + 1 : 0;
      tr3_stats st;
      memset(&st, 0, sizeof st);
      /* §744: ПРЯМОЙ прогон (альбедо 0) до полного — E_dir поэлементно для
       * сличения; полный прогон ниже перезапишет phi (warm_start = 0). */
      double *eind_dir = NULL;
      if (xcmp) {
        double *fz2 = calloc((size_t)ftab.n, sizeof *fz2);
        if (fz2 == NULL) exit(1);
        tr3_problem pd = prob;
        pd.facet_rho = fz2;
        /* §746: «прямой» обязан быть прямым — стык-отскоки принадлежат
         * КОСВЕННОМУ, каким бы ни был ключ xsolidrho */
        pd.solid_rho = 0.0;
        tr3_stats std;
        memset(&std, 0, sizeof std);
        if (tr3_sweep_solve(&pd, xit, xtol, phi, &std) == 0) {
          eind_dir = std.eirr;
          std.eirr = NULL;
        }
        free(std.bout);
        free(std.sout);
        free(std.eirr);
        free(fz2);
      }
      double tsw = now_s();
      int src = 0;
      if (xdsa) {
        /* §750: ДВУХСЕТОЧНЫЙ ЦИКЛ. Такт развёртки (прокидка К76) → приращение
         * sout ограничивается на срез ядра §597 → хвост ряда добирается
         * грубым оператором (etree_build + gather_apply по хранимым связям) →
         * поправка продлевается в sout. Единицы моста: радиосити B = π·радианс
         * (перевод расчётный, γ по умолчанию 1 и НЕ подгоняется — А1170). */
        if (!xdsagal && (!g_keonly || !(g_hgather > 0.0))) {
          fprintf(stderr, "xdsa требует keonly и hgather= (или xdsagal)\n");
          exit(1);
        }
        /* §752: ГАЛЁРКИНСКОЕ ГРУБОЕ ЗВЕНО ИЗ САМОГО СВИПА. Агрегаты — коробки
         * октодерева стороной 2^xdsalev листьев; грубая переменная — площадно-
         * взвешенная исходящая радиансность агрегата (элементы + стык-грани,
         * А1172). A_c = R·S·P зондами: индикатор на агрегат, один такт (ω=1,
         * эмиссия 0), ограничение ответа — столбец. Единицы по построению
         * (радианс туда и обратно, А1170 исчезает). Клипы однородны, но не
         * аддитивны — A_c есть линеаризация вокруг индикаторов (А1176). */
        int32_t nagg = 0;
        int32_t *aggof_e = NULL, *fagg = NULL, *piv9 = NULL;
        double *aden = NULL, *LU9 = NULL, *cagg = NULL, *dagg = NULL, *ubprev = NULL;
        double *farea = NULL;
        /* §754: рабочее поле и базовый ответ разностных зондов */
        double *usB = NULL, *ubB = NULL, *phiB = NULL, *r0agg = NULL;
        double eps754 = 0.0;
        if (xdsagal) {
          int32_t gside = (1 << lev) >> xdsalev;
          if (gside < 1) gside = 1;
          int32_t ngrid = 2 * gside * gside * gside; /* §756: два слота на коробку при xdsasplit */
          int32_t *gid = malloc((size_t)ngrid * sizeof *gid);
          if (gid == NULL) exit(1);
          for (int32_t i = 0; i < ngrid; i++)
            gid[i] = -1;
          aggof_e = malloc((size_t)(cut.nse > 0 ? cut.nse : 1) * sizeof *aggof_e);
          fagg = malloc((size_t)mesh.nf * sizeof *fagg);
          farea = calloc((size_t)mesh.nf, sizeof *farea);
          if (aggof_e == NULL || fagg == NULL || farea == NULL) exit(1);
          /* raw-номер агрегата ячейки */
          int32_t *craw = malloc((size_t)mesh.ncell * sizeof *craw);
          if (craw == NULL) exit(1);
          for (int32_t ci = 0; ci < mesh.ncell; ci++) {
            int32_t gx = mesh.clo[ci][0] >> xdsalev, gy = mesh.clo[ci][1] >> xdsalev,
                    gz = mesh.clo[ci][2] >> xdsalev;
            craw[ci] = (gz * gside + gy) * gside + gx;
          }
          for (int32_t e = 0; e < cut.nse; e++) {
            int32_t ci = cut.se[e].cell;
            aggof_e[e] = -2;
            if (ci < 0 || ci >= mesh.ncell || !(cut.se[e].area > 0.0)) continue;
            int32_t r = xdsasplit ? 2 * craw[ci] : craw[ci]; /* §756 */
            if (gid[r] < 0) gid[r] = nagg++;
            aggof_e[e] = gid[r];
          }
          for (int32_t f = 0; f < mesh.nf; f++) {
            fagg[f] = -2;
            if (mesh.f[f].cb < 0) continue;
            int as9 = cut.solid[mesh.f[f].ca] ? 1 : 0, bs9 = cut.solid[mesh.f[f].cb] ? 1 : 0;
            if (as9 == bs9) continue; /* не стык */
            int32_t cf9 = as9 ? mesh.f[f].cb : mesh.f[f].ca;
            farea[f] = cut.ffm[f][0][0] > 0.0 ? cut.ffm[f][0][0] : cut.ffmb[f][0][0];
            if (!(farea[f] > 0.0)) continue;
            int32_t r = xdsasplit ? 2 * craw[cf9] + 1 : craw[cf9]; /* §756: стык — свой слот */
            if (gid[r] < 0) gid[r] = nagg++;
            fagg[f] = gid[r];
          }
          free(craw);
          free(gid);
          if (nagg < 2 || nagg > 2048) {
            fprintf(stderr, "xdsagal: агрегатов %d — вне [2, 2048]\n", nagg);
            exit(1);
          }
          aden = calloc((size_t)nagg, sizeof *aden);
          double *Ac = calloc((size_t)nagg * (size_t)nagg, sizeof *Ac);
          cagg = calloc((size_t)nagg, sizeof *cagg);
          dagg = calloc((size_t)nagg, sizeof *dagg);
          ubprev = calloc((size_t)mesh.nf, sizeof *ubprev);
          if (aden == NULL || Ac == NULL || cagg == NULL || dagg == NULL || ubprev == NULL) exit(1);
          for (int32_t e = 0; e < cut.nse; e++)
            if (aggof_e[e] >= 0) aden[aggof_e[e]] += cut.se[e].area;
          for (int32_t f = 0; f < mesh.nf; f++)
            if (fagg[f] >= 0) aden[fagg[f]] += farea[f];
          /* зонды */
          double tprobe = now_s();
          double *usP = calloc((size_t)(cut.nse > 0 ? cut.nse : 1) * 4, sizeof *usP);
          double *ubP = calloc((size_t)mesh.nf * 4, sizeof *ubP);
          double *phi9 = calloc((size_t)mesh.ncell * 4, sizeof *phi9);
          if (usP == NULL || ubP == NULL || phi9 == NULL) exit(1);
          /* §754: РАЗНОСТНЫЕ ЗОНДЫ ВОКРУГ РАБОЧЕГО ПОЛЯ (А1180а). x̄ — после
           * NBASE демпфированных тактов; столбец = (R·S(x̄+ε·инд) − R·S(x̄))/ε;
           * эмиссия сокращается разностью точно, клипы — в рабочем режиме.
           * ε = 0.1·max|sout(x̄)| (десятая масштаба поля) × xdsaeps (НК ×2). */

          if (xdsadiff) {
            enum { NBASE754 = 8 }; /* рабочая окрестность; выбор — А1182 */
            phiB = calloc((size_t)mesh.ncell * 4, sizeof *phiB);
            r0agg = calloc((size_t)nagg, sizeof *r0agg);
            if (phiB == NULL || r0agg == NULL) exit(1);
            tr3_problem pdB = prob;
            pdB.relax = xrelax > 0.0 ? xrelax : 0.7;
            pdB.warm_start = 1;
            pdB.trace = 0;
            for (int k = 0; k < NBASE754; k++) {
              pdB.bout_in = ubB;
              pdB.sout_in = usB;
              tr3_stats stB;
              memset(&stB, 0, sizeof stB);
              if (tr3_sweep_solve(&pdB, 1, 0.0, phiB, &stB) != 0) exit(1);
              free(ubB);
              free(usB);
              free(stB.eirr);
              ubB = stB.bout;
              usB = stB.sout;
            }
            for (int32_t e = 0; e < cut.nse; e++)
              if (fabs(usB[4 * (size_t)e]) > eps754) eps754 = fabs(usB[4 * (size_t)e]);
            eps754 *= 0.1 * xdsaeps;
            if (!(eps754 > 0.0)) {
              fprintf(stderr, "xdsadiff: рабочее поле нулевое, ε не определить\n");
              exit(1);
            }
            /* базовый ответ R·S(x̄), ω=1 */
            tr3_problem pd0 = prob;
            pd0.relax = 0.0;
            pd0.warm_start = 1;
            pd0.trace = 0;
            pd0.bout_in = ubB;
            pd0.sout_in = usB;
            memcpy(phi9, phiB, (size_t)mesh.ncell * 4 * sizeof *phi9);
            tr3_stats st0;
            memset(&st0, 0, sizeof st0);
            if (tr3_sweep_solve(&pd0, 1, 0.0, phi9, &st0) != 0) exit(1);
            for (int32_t e = 0; e < cut.nse; e++)
              if (aggof_e[e] >= 0) r0agg[aggof_e[e]] += st0.sout[4 * (size_t)e] * cut.se[e].area;
            for (int32_t f = 0; f < mesh.nf; f++)
              if (fagg[f] >= 0) r0agg[fagg[f]] += st0.bout[4 * (size_t)f] * farea[f];
            free(st0.bout);
            free(st0.sout);
            free(st0.eirr);
            printf("   §754 БАЗА: %d тактов, ε = %.4e (0.1·max|sout| × %.2g)\n", NBASE754, eps754,
                   xdsaeps);
          }
          for (int32_t i = 0; i < nagg; i++) {
            if (xdsadiff) {
              memcpy(usP, usB, (size_t)(cut.nse > 0 ? cut.nse : 1) * 4 * sizeof *usP);
              memcpy(ubP, ubB, (size_t)mesh.nf * 4 * sizeof *ubP);
              memcpy(phi9, phiB, (size_t)mesh.ncell * 4 * sizeof *phi9);
              for (int32_t e = 0; e < cut.nse; e++)
                if (aggof_e[e] == i) usP[4 * (size_t)e] += eps754;
              for (int32_t f = 0; f < mesh.nf; f++)
                if (fagg[f] == i) ubP[4 * (size_t)f] += eps754;
            } else {
              memset(usP, 0, (size_t)(cut.nse > 0 ? cut.nse : 1) * 4 * sizeof *usP);
              memset(ubP, 0, (size_t)mesh.nf * 4 * sizeof *ubP);
              memset(phi9, 0, (size_t)mesh.ncell * 4 * sizeof *phi9);
              for (int32_t e = 0; e < cut.nse; e++)
                if (aggof_e[e] == i) usP[4 * (size_t)e] = 1.0;
              for (int32_t f = 0; f < mesh.nf; f++)
                if (fagg[f] == i) ubP[4 * (size_t)f] = 1.0;
            }
            tr3_problem pd9 = prob;
            if (!xdsadiff) {
              /* индикаторные зонды: эмиссия в ноль (разностные сокращают её
               * разностью и потому идут с ЖИВОЙ эмиссией) */
              pd9.elem_emit = NULL;
              pd9.facet_emit = NULL;
            }
            pd9.relax = 0.0; /* ω=1: зондируется сырой S */
            pd9.warm_start = 1;
            pd9.bout_in = ubP;
            pd9.sout_in = usP;
            pd9.trace = 0;
            tr3_stats st9;
            memset(&st9, 0, sizeof st9);
            if (tr3_sweep_solve(&pd9, 1, 0.0, phi9, &st9) != 0) exit(1);
            for (int32_t e = 0; e < cut.nse; e++)
              if (aggof_e[e] >= 0)
                Ac[(size_t)aggof_e[e] * (size_t)nagg + (size_t)i] +=
                    st9.sout[4 * (size_t)e] * cut.se[e].area;
            for (int32_t f = 0; f < mesh.nf; f++)
              if (fagg[f] >= 0)
                Ac[(size_t)fagg[f] * (size_t)nagg + (size_t)i] +=
                    st9.bout[4 * (size_t)f] * farea[f];
            /* §754: столбец разностью вокруг базового ответа, до нормировки
             * на площадь агрегата (r0agg накоплен теми же весами) */
            if (xdsadiff)
              for (int32_t j = 0; j < nagg; j++)
                Ac[(size_t)j * (size_t)nagg + (size_t)i] =
                    (Ac[(size_t)j * (size_t)nagg + (size_t)i] - r0agg[j]) / eps754;
            free(st9.bout);
            free(st9.sout);
            free(st9.eirr);
          }
          for (int32_t j = 0; j < nagg; j++)
            for (int32_t i = 0; i < nagg; i++)
              Ac[(size_t)j * (size_t)nagg + (size_t)i] /= aden[j] > 0.0 ? aden[j] : 1.0;
          free(usP);
          free(ubP);
          free(phi9);
          /* спектр A_c степенной итерацией — ключевой тест согласованности */
          double rho9 = 0.0;
          {
            double *v9 = malloc((size_t)nagg * sizeof *v9);
            double *w9v = malloc((size_t)nagg * sizeof *w9v);
            if (v9 == NULL || w9v == NULL) exit(1);
            for (int32_t i = 0; i < nagg; i++)
              v9[i] = 1.0;
            for (int itp = 0; itp < 100; itp++) {
              double nn = 0.0;
              for (int32_t j = 0; j < nagg; j++) {
                double s9 = 0.0;
                for (int32_t i = 0; i < nagg; i++)
                  s9 += Ac[(size_t)j * (size_t)nagg + (size_t)i] * v9[i];
                w9v[j] = s9;
                nn += s9 * s9;
              }
              nn = sqrt(nn);
              rho9 = nn / sqrt((double)nagg);
              if (nn > 0.0)
                for (int32_t j = 0; j < nagg; j++)
                  v9[j] = w9v[j] / nn;
              double vn = 0.0;
              for (int32_t j = 0; j < nagg; j++)
                vn += v9[j] * v9[j];
              rho9 = nn; /* после нормировки v: ‖A v‖ при ‖v‖=1 */
              (void)vn;
            }
            free(v9);
            free(w9v);
          }
          /* LU (I − A_c) с частичным выбором */
          LU9 = malloc((size_t)nagg * (size_t)nagg * sizeof *LU9);
          piv9 = malloc((size_t)nagg * sizeof *piv9);
          if (LU9 == NULL || piv9 == NULL) exit(1);
          for (int32_t j = 0; j < nagg; j++)
            for (int32_t i = 0; i < nagg; i++)
              LU9[(size_t)j * (size_t)nagg + (size_t)i] =
                  (i == j ? 1.0 : 0.0) - Ac[(size_t)j * (size_t)nagg + (size_t)i];
          for (int32_t k = 0; k < nagg; k++) {
            int32_t bp = k;
            for (int32_t j = k + 1; j < nagg; j++)
              if (fabs(LU9[(size_t)j * (size_t)nagg + (size_t)k]) >
                  fabs(LU9[(size_t)bp * (size_t)nagg + (size_t)k]))
                bp = j;
            piv9[k] = bp;
            if (bp != k)
              for (int32_t i = 0; i < nagg; i++) {
                double t9s = LU9[(size_t)k * (size_t)nagg + (size_t)i];
                LU9[(size_t)k * (size_t)nagg + (size_t)i] =
                    LU9[(size_t)bp * (size_t)nagg + (size_t)i];
                LU9[(size_t)bp * (size_t)nagg + (size_t)i] = t9s;
              }
            double pv9 = LU9[(size_t)k * (size_t)nagg + (size_t)k];
            if (!(fabs(pv9) > 0.0)) {
              fprintf(stderr, "xdsagal: вырождение LU на шаге %d\n", k);
              exit(1);
            }
            for (int32_t j = k + 1; j < nagg; j++) {
              double m9 = LU9[(size_t)j * (size_t)nagg + (size_t)k] / pv9;
              LU9[(size_t)j * (size_t)nagg + (size_t)k] = m9;
              for (int32_t i = k + 1; i < nagg; i++)
                LU9[(size_t)j * (size_t)nagg + (size_t)i] -=
                    m9 * LU9[(size_t)k * (size_t)nagg + (size_t)i];
            }
          }
          free(Ac);
          printf("   §752 ГАЛЁРКИН: агрегатов %d (коробка %d листьев), зонды %.1f с; "
                 "ρ(A_c) = %.4f\n",
                 nagg, 1 << xdsalev, now_s() - tprobe, rho9);
        }
        int32_t NS = 0;
        if (!xdsagal) {
          arealight AL9;
          memset(&AL9, 0, sizeof AL9); /* keonly: AL ядром не читается */
          g_dsa_keep = 1;
          ind_core_build(&T, &ht, &fr, &P, &AL9, &m, &CT, lev, g_hgather, indvis);
          g_dsa_keep = 0;
          if (!g_dsa_ready) exit(1);
          NS = g_dsa_SF.n;
        }
        /* мосты: узел -> срез, ячейка меша -> срез; отказы — счётчиками (А891) */
        int32_t *n2s = malloc((size_t)T.n * sizeof *n2s);
        int32_t *c2s = malloc((size_t)mesh.ncell * sizeof *c2s);
        if (n2s == NULL || c2s == NULL) exit(1);
        for (int32_t i = 0; i < T.n; i++)
          n2s[i] = -1;
        for (int32_t i = 0; i < NS; i++) {
          int32_t ni = node_of_cell(&T, lev, &g_dsa_SF.c[i]);
          if (ni >= 0 && ni < T.n) n2s[ni] = i;
        }
        for (int32_t ci = 0; ci < mesh.ncell; ci++) {
          int32_t szc = mesh.csize[ci];
          int lvl9 = lev;
          while (szc > 1) {
            szc >>= 1;
            lvl9--;
          }
          hz_dccell q9;
          memset(&q9, 0, sizeof q9);
          q9.lvl = (uint8_t)lvl9;
          for (int a = 0; a < 3; a++)
            q9.lo[a] = (uint16_t)mesh.clo[ci][a];
          int32_t ni = node_of_cell(&T, lev, &q9);
          c2s[ci] = (ni >= 0 && ni < T.n) ? n2s[ni] : -1;
        }
        double acov = 0.0, atot = 0.0;
        for (int32_t e = 0; e < cut.nse; e++) {
          double a = cut.se[e].area;
          if (!(a > 0.0)) continue;
          atot += a;
          int32_t ci = cut.se[e].cell;
          if (ci >= 0 && ci < mesh.ncell && c2s[ci] >= 0) acov += a;
        }
        if (!xdsagal)
          printf("   §750 МОСТ: покрытие элементов срезом %.1f %% площади\n",
                 100.0 * acov / (atot > 0.0 ? atot : 1.0));
        float *bc9 = malloc(3 * (size_t)(NS > 0 ? NS : 1) * sizeof *bc9);
        float *bn9 = malloc(3 * (size_t)(NS > 0 ? NS : 1) * sizeof *bn9);
        double *cslice = calloc((size_t)(NS > 0 ? NS : 1), sizeof *cslice);
        double *dnum = calloc((size_t)(NS > 0 ? NS : 1), sizeof *dnum);
        double *dden = calloc((size_t)(NS > 0 ? NS : 1), sizeof *dden);
        double *usprev = calloc((size_t)(cut.nse > 0 ? cut.nse : 1), sizeof *usprev);
        if (bc9 == NULL || bn9 == NULL || cslice == NULL || dnum == NULL || dden == NULL ||
            usprev == NULL)
          exit(1);
        prob.warm_start = 1;
        if (!(xrelax > 0.0)) prob.relax = 0.7; /* ω из §749 (А1166) */
        double *ub9 = NULL, *us9 = NULL;
        if (xdsagal && xdsadiff && usB != NULL) {
          /* §754: главный цикл ПРОДОЛЖАЕТ с рабочего поля x̄ — базовые такты не
           * выбрасываются; usprev/ubprev от x̄, чтобы первый δ был мал */
          ub9 = ubB;
          us9 = usB;
          ubB = NULL;
          usB = NULL;
          memcpy(phi, phiB, (size_t)mesh.ncell * 4 * sizeof *phi);
          for (int32_t e = 0; e < cut.nse; e++)
            usprev[e] = us9[4 * (size_t)e];
          if (ubprev != NULL)
            for (int32_t f = 0; f < mesh.nf; f++)
              ubprev[f] = ub9[4 * (size_t)f];
        }
        for (int it9 = 0; it9 < xit; it9++) {
          double t9 = now_s();
          prob.bout_in = ub9;
          prob.sout_in = us9;
          memset(&st, 0, sizeof st);
          src = tr3_sweep_solve(&prob, 1, 0.0, phi, &st);
          if (src != 0) break;
          /* §752: галёркинское звено — δ по агрегатам (sout И стык-bout),
           * деление на ω (демпфированное приращение = ω·(S−x)), плотное LU,
           * продление кусочно-постоянно в оба носителя. */
          if (xdsagal) {
            double w9 = prob.relax > 0.0 ? prob.relax : 1.0;
            memset(dagg, 0, (size_t)nagg * sizeof *dagg);
            for (int32_t e = 0; e < cut.nse; e++)
              if (aggof_e[e] >= 0)
                dagg[aggof_e[e]] += (st.sout[4 * (size_t)e] - usprev[e]) * cut.se[e].area;
            for (int32_t f = 0; f < mesh.nf; f++)
              if (fagg[f] >= 0) dagg[fagg[f]] += (st.bout[4 * (size_t)f] - ubprev[f]) * farea[f];
            for (int32_t j = 0; j < nagg; j++)
              dagg[j] /= (aden[j] > 0.0 ? aden[j] : 1.0) * w9;
            /* LU-решение (I − A_c)·c = δ */
            for (int32_t k = 0; k < nagg; k++) {
              if (piv9[k] != k) {
                double t9s = dagg[k];
                dagg[k] = dagg[piv9[k]];
                dagg[piv9[k]] = t9s;
              }
              for (int32_t j = k + 1; j < nagg; j++)
                dagg[j] -= LU9[(size_t)j * (size_t)nagg + (size_t)k] * dagg[k];
            }
            for (int32_t k = nagg - 1; k >= 0; k--) {
              for (int32_t i = k + 1; i < nagg; i++)
                dagg[k] -= LU9[(size_t)k * (size_t)nagg + (size_t)i] * cagg[i];
              cagg[k] = dagg[k] / LU9[(size_t)k * (size_t)nagg + (size_t)k];
            }
            double cmx9 = 0.0;
            for (int32_t e = 0; e < cut.nse; e++)
              if (aggof_e[e] >= 0) {
                double dc = xdsagain * cagg[aggof_e[e]];
                st.sout[4 * (size_t)e] += dc;
                if (fabs(dc) > cmx9) cmx9 = fabs(dc);
              }
            for (int32_t f = 0; f < mesh.nf; f++)
              if (fagg[f] >= 0) st.bout[4 * (size_t)f] += xdsagain * cagg[fagg[f]];
            for (int32_t e = 0; e < cut.nse; e++)
              usprev[e] = st.sout[4 * (size_t)e];
            for (int32_t f = 0; f < mesh.nf; f++)
              ubprev[f] = st.bout[4 * (size_t)f];
            printf("   §752 такт %2d: resid %.3e (rel %.1e, max|состояние| %.3e); max|поправка| "
                   "%.3e; %.1f с\n",
                   it9, st.resid, st.resid_rel, st.statemax, cmx9, now_s() - t9);
            free(ub9);
            free(us9);
            free(st.eirr);
            st.eirr = NULL;
            ub9 = st.bout;
            us9 = st.sout;
            if (st.resid < xtol && it9 > 2) break;
            continue;
          }
          /* приращение -> срез (площадно-взвешенно, радиосити = π·радианс) */
          memset(dnum, 0, (size_t)NS * sizeof *dnum);
          memset(dden, 0, (size_t)NS * sizeof *dden);
          for (int32_t e = 0; e < cut.nse; e++) {
            double a = cut.se[e].area;
            int32_t ci = cut.se[e].cell;
            if (!(a > 0.0) || ci < 0 || ci >= mesh.ncell || c2s[ci] < 0) continue;
            dnum[c2s[ci]] += (st.sout[4 * (size_t)e] - usprev[e]) * a;
            dden[c2s[ci]] += a;
          }
          for (int32_t s = 0; s < NS; s++)
            dnum[s] = dden[s] > 0.0 ? 3.14159265358979323846 * dnum[s] / dden[s] : 0.0;
          /* грубое звено: c <- A(delta + c), тёплый старт с прошлого такта;
           * остановка — сжатие грубой невязки 10x либо потолок 30 (из плана) */
          int mco = 0;
          double cr0 = -1.0, crn = -1.0;
          for (; mco < 30; mco++) {
            for (int32_t s = 0; s < NS; s++) {
              float v = (float)(dnum[s] + cslice[s]);
              bc9[3 * (size_t)s] = v;
              bc9[3 * (size_t)s + 1] = v;
              bc9[3 * (size_t)s + 2] = v;
            }
            etree ET9;
            memset(&ET9, 0, sizeof ET9);
            etree_build(&ET9, &g_dsa_SF, &fr, bc9, &m, 0, NS, 0, lev);
            gather_apply(&g_dsa_LC, &ET9, &g_dsa_SF, &m, g_dsa_ffv, bn9);
            etree_free(&ET9);
            double dr = 0.0;
            for (int32_t s = 0; s < NS; s++) {
              double nv9 = ((double)bn9[3 * (size_t)s] + (double)bn9[3 * (size_t)s + 1] +
                            (double)bn9[3 * (size_t)s + 2]) /
                           3.0;
              double dd = fabs(nv9 - cslice[s]);
              if (dd > dr) dr = dd;
              cslice[s] = nv9;
            }
            if (mco == 0) cr0 = dr;
            crn = dr;
            if (crn <= 0.1 * cr0) {
              mco++;
              break;
            }
          }
          /* продление поправки в sout (нулевой коэффициент) */
          double cmx9 = 0.0;
          for (int32_t e = 0; e < cut.nse; e++) {
            int32_t ci = cut.se[e].cell;
            if (ci < 0 || ci >= mesh.ncell || c2s[ci] < 0) continue;
            double dc = xdsagain * cslice[c2s[ci]] / 3.14159265358979323846;
            st.sout[4 * (size_t)e] += dc;
            if (fabs(dc) > cmx9) cmx9 = fabs(dc);
          }
          for (int32_t e = 0; e < cut.nse; e++)
            usprev[e] = st.sout[4 * (size_t)e];
          printf("   §750 такт %2d: resid %.3e (rel %.1e, max|состояние| %.3e); грубое звено %d "
                 "итер (%.2e -> %.2e); max|поправка| %.3e; %.1f с\n",
                 it9, st.resid, st.resid_rel, st.statemax, mco, cr0, crn, cmx9, now_s() - t9);
          free(ub9);
          free(us9);
          free(st.eirr);
          st.eirr = NULL;
          ub9 = st.bout;
          us9 = st.sout;
          if (st.resid < xtol && it9 > 2) break;
        }
        free(n2s);
        free(aggof_e);
        free(fagg);
        free(farea);
        free(aden);
        free(LU9);
        free(piv9);
        free(cagg);
        free(dagg);
        free(ubprev);
        free(usB);
        free(ubB);
        free(phiB);
        free(r0agg);
        free(c2s);
        free(bc9);
        free(bn9);
        free(cslice);
        free(dnum);
        free(dden);
        free(usprev);
      } else if (xbounce > 0) {
        /* §774: ОГРУБЛЕНИЕ ЧАСТОТЫ ОТСКОКОВ. Отскоки 1…N считаются точно —
         * N прокидок К76 (такт = отскок, А1222: внутри такта марш переносит
         * эмиссию через весь домен, поверхности обновляются после цикла
         * направлений); все поздние отскоки — одним агрегатом, геометрическим
         * замыканием хвоста: x_inf ≈ x_N + Δ_N·q̂/(1−q̂), где q̂ измеряется
         * скалярным произведением двух последних приращений объединённого
         * состояния (φ, bout, sout). Замыкаются состояние, eirr и скаляры К40
         * (все — линейные функционалы поля такта). Клип замыкания не ставится
         * (А1227): отрицательные eirr печатаются счётчиком, а не прячутся. */
        /* Потолок q̂ и минимум тактов — HZ_TAILQ_MAX/HZ_TAILN_MIN у §808-блока
         * (файловая область видимости: та же лестница в стенде кадра 2). */
        const int32_t nse4 = cut.nse > 0 ? cut.nse : 1;
        const size_t nph = (size_t)mesh.ncell * 4, nbo = (size_t)mesh.nf * 4,
                     nso = (size_t)nse4 * 4;
        /* базовый прогон для сравнения — в том же процессе (xbcmp) */
        double *phiB4 = NULL, *eirrB = NULL;
        double psinB = 0.0;
        if (xbcmp774) {
          phiB4 = calloc(nph, sizeof *phiB4);
          if (phiB4 == NULL) exit(1);
          tr3_stats stB4;
          memset(&stB4, 0, sizeof stB4);
          double tb0 = now_s();
          if (tr3_sweep_solve(&prob, xit, xtol, phiB4, &stB4) != 0) exit(1);
          double sphiB = 0.0, mphiB = 0.0;
          for (size_t i = 0; i < nph; i += 4) {
            sphiB += fabs(phiB4[i]);
            if (fabs(phiB4[i]) > mphiB) mphiB = fabs(phiB4[i]);
          }
          /* счётчик отрицательных eirr БАЗЫ: без него отрицательные лестницы
           * неатрибуируемы — замыкание против исходного состояния (А11-класс) */
          int64_t nnegB = 0;
          if (stB4.eirr != NULL)
            for (int32_t e = 0; e < cut.nse; e++)
              if (stB4.eirr[e] < 0.0) nnegB++;
          printf("   §774 БАЗА (xit=%d, xtol=%.1e): итераций %d, невязка %.3e (rel %.1e, "
                 "max|состояние| %.3e), psin %.6g, %.2f с; слепок Σ|φ| %.17g, max|φ| %.17g; "
                 "eirr < 0 у %lld из %d\n",
                 xit, xtol, stB4.iters, stB4.resid, stB4.resid_rel, stB4.statemax, stB4.psin,
                 now_s() - tb0, sphiB, mphiB, (long long)nnegB, cut.nse);
          eirrB = stB4.eirr;
          psinB = stB4.psin;
          free(stB4.bout);
          free(stB4.sout);
        }
        double *pphi = calloc(nph, sizeof *pphi); /* состояние такта t−1 */
        double *pbo = calloc(nbo, sizeof *pbo);
        double *pso = calloc(nso, sizeof *pso);
        double *dpphi = calloc(nph, sizeof *dpphi); /* приращение такта t−1 */
        double *dpbo = calloc(nbo, sizeof *dpbo);
        double *dpso = calloc(nso, sizeof *dpso);
        double *peirr = calloc((size_t)nse4, sizeof *peirr); /* eirr такта t−1 */
        if (pphi == NULL || pbo == NULL || pso == NULL || dpphi == NULL || dpbo == NULL ||
            dpso == NULL || peirr == NULL)
          exit(1);
        double *ubL = NULL, *usL = NULL; /* прокидка К76 */
        double q_num = 0.0, q_den = 0.0, denprev = 0.0;
        double sc_prev[6] = {0, 0, 0, 0, 0, 0}; /* скаляры К40 такта N−1 */
        tr3_stats stfin;
        memset(&stfin, 0, sizeof stfin);
        double tlad = now_s();
        for (int t = 1; t <= xbounce; t++) {
          tr3_problem pl = prob;
          tr3_stats stt;
          memset(&stt, 0, sizeof stt);
          if (t > 1) {
            pl.warm_start = 1;
            pl.bout_in = ubL;
            pl.sout_in = usL;
          }
          double tt0 = now_s();
          if (tr3_sweep_solve(&pl, 1, 0.0, phi, &stt) != 0) exit(1);
          /* приращение такта, его норма и скалярные произведения с прошлым */
          double dmaxphi = 0.0, dmaxso = 0.0, num_t = 0.0, dencur = 0.0;
          for (size_t i = 0; i < nph; i++) {
            double dc = phi[i] - pphi[i];
            if (fabs(dc) > dmaxphi) dmaxphi = fabs(dc);
            num_t += dc * dpphi[i];
            dencur += dc * dc;
            dpphi[i] = dc;
            pphi[i] = phi[i];
          }
          for (size_t i = 0; i < nbo; i++) {
            double dc = stt.bout[i] - pbo[i];
            num_t += dc * dpbo[i];
            dencur += dc * dc;
            dpbo[i] = dc;
            pbo[i] = stt.bout[i];
          }
          for (size_t i = 0; i < nso; i++) {
            double dc = stt.sout[i] - pso[i];
            if (fabs(dc) > dmaxso) dmaxso = fabs(dc);
            num_t += dc * dpso[i];
            dencur += dc * dc;
            dpso[i] = dc;
            pso[i] = stt.sout[i];
          }
          if (t > 1) {
            q_num = num_t;
            q_den = denprev;
          }
          denprev = dencur;
          printf("   §774 такт %2d: |Δφ|∞ %.3e, |Δsout|∞ %.3e, psin %.6g, %.2f с\n", t, dmaxphi,
                 dmaxso, stt.psin, now_s() - tt0);
          free(ubL); /* состояние t−1 больше не нужно: копии лежат в pbo/pso */
          free(usL);
          ubL = stt.bout;
          usL = stt.sout;
          if (t < xbounce) {
            if (stt.eirr != NULL)
              for (int32_t e = 0; e < cut.nse; e++)
                peirr[e] = stt.eirr[e];
            sc_prev[0] = stt.pin;
            sc_prev[1] = stt.pout;
            sc_prev[2] = stt.pabs;
            sc_prev[3] = stt.psin;
            sc_prev[4] = stt.psout;
            sc_prev[5] = stt.psolid;
            free(stt.eirr);
          } else
            stfin = stt; /* bout/sout = ubL/usL, освобождает общий путь */
        }
        st = stfin;
        printf("   §774 ЛЕСТНИЦА: %d тактов за %.2f с\n", xbounce, now_s() - tlad);
        const double t774lad = now_s() - tlad; /* §810: база сравнения — та же
                                                * величина, что напечатана выше */
        /* замыкание хвоста */
        double qhat = -1.0;
        if (xtailq > 0.0)
          qhat = xtailq; /* НК: форсированный множитель */
        else if (xtailq < 0.0 && xbounce >= HZ_TAILN_MIN && q_den > 0.0)
          qhat = q_num / q_den;
        if (qhat > 0.0 && qhat <= HZ_TAILQ_MAX) {
          double mult = qhat / (1.0 - qhat);
          double psin_raw = st.psin;
          for (size_t i = 0; i < nph; i++)
            phi[i] += mult * dpphi[i];
          for (size_t i = 0; i < nbo; i++)
            st.bout[i] += mult * dpbo[i];
          for (size_t i = 0; i < nso; i++)
            st.sout[i] += mult * dpso[i];
          int64_t nneg = 0;
          if (st.eirr != NULL)
            for (int32_t e = 0; e < cut.nse; e++) {
              st.eirr[e] += mult * (st.eirr[e] - peirr[e]);
              if (st.eirr[e] < 0.0) nneg++;
            }
          st.pin += mult * (st.pin - sc_prev[0]);
          st.pout += mult * (st.pout - sc_prev[1]);
          st.pabs += mult * (st.pabs - sc_prev[2]);
          st.psin += mult * (st.psin - sc_prev[3]);
          st.psout += mult * (st.psout - sc_prev[4]);
          st.psolid += mult * (st.psolid - sc_prev[5]);
          printf("   §774 ЗАМЫКАНИЕ: q̂ %.4f (%s), множитель %.3f; psin %.6g -> %.6g; "
                 "eirr < 0 у %lld элементов из %d\n",
                 qhat, xtailq > 0.0 ? "ФОРСИРОВАН — НК" : "измерен", mult, psin_raw, st.psin,
                 (long long)nneg, cut.nse);
        } else if (qhat > HZ_TAILQ_MAX)
          printf("   §774 замыкание ПРОПУЩЕНО: q̂ %.4f вне (0, %.2f] — ЧИСТОЕ УСЕЧЕНИЕ "
                 "(fail closed)\n",
                 qhat, HZ_TAILQ_MAX);
        else
          printf("   §774 замыкание ВЫКЛЮЧЕНО (%s) — ЧИСТОЕ УСЕЧЕНИЕ\n",
                 xtailq < 0.0 ? "тактов меньше минимума либо нулевое приращение"
                              : "xtailq=0 по ключу");
        {
          double sphiL = 0.0, mphiL = 0.0;
          for (size_t i = 0; i < nph; i += 4) {
            sphiL += fabs(phi[i]);
            if (fabs(phi[i]) > mphiL) mphiL = fabs(phi[i]);
          }
          printf("   §774 СЛЕПОК: Σ|φ| %.17g, max|φ| %.17g\n", sphiL, mphiL);
        }
        /* сравнение с базой того же процесса */
        if (xbcmp774) {
          int bitphi = memcmp(phi, phiB4, nph * sizeof(double)) == 0;
          int biteirr = eirrB != NULL && st.eirr != NULL &&
                        memcmp(st.eirr, eirrB, (size_t)cut.nse * sizeof(double)) == 0;
          double dps = psinB > 0.0 ? (st.psin - psinB) / psinB : 0.0;
          /* зона §768: < 15 м от глаза, лицевые; предикат тот же, что у
           * прибора вклада — сравнение читается против А1212 */
          double zb = 0.0, zl = 0.0, tb = 0.0, tl = 0.0;
          int64_t nzone = 0, nuse = 0, nneg2 = 0;
          double *rel = malloc((size_t)nse4 * sizeof *rel);
          if (rel == NULL) exit(1);
          for (int32_t e = 0; e < cut.nse && eirrB != NULL && st.eirr != NULL; e++) {
            double c9[3] = {0, 0, 0};
            for (int q2 = 0; q2 < cut.se[e].nv; q2++)
              for (int a = 0; a < 3; a++)
                c9[a] += cut.se[e].v[q2][a] / (double)(cut.se[e].nv > 0 ? cut.se[e].nv : 1);
            double d2 = 0.0, dot = 0.0;
            for (int a = 0; a < 3; a++) {
              double dd = c9[a] - g_eye[a];
              d2 += dd * dd;
              dot += cut.se[e].n[a] * dd;
            }
            int inzone = sqrt(d2) < 15.0 && !(dot > 0.0);
            double wb = eirrB[e] * cut.se[e].area, wl = st.eirr[e] * cut.se[e].area;
            tb += fabs(wb);
            tl += fabs(wl);
            if (inzone) {
              zb += fabs(wb);
              zl += fabs(wl);
              nzone++;
            }
            if (fabs(eirrB[e]) > 0.0) rel[nuse++] = fabs(st.eirr[e] - eirrB[e]) / fabs(eirrB[e]);
            if (st.eirr[e] < 0.0) nneg2++;
          }
          printf("   §774 СРАВНЕНИЕ С БАЗОЙ: Δpsin %+.3f %%; Σ|E·area| зоны (<15 м, лицевые, "
                 "%lld элементов) %.6g против %.6g (Δ %+.3f %%); всей сцены %.6g против %.6g "
                 "(Δ %+.3f %%); ПОБИТОВО φ: %s, eirr: %s\n",
                 100.0 * dps, (long long)nzone, zl, zb, 100.0 * (zl - zb) / (zb > 0.0 ? zb : 1.0),
                 tl, tb, 100.0 * (tl - tb) / (tb > 0.0 ? tb : 1.0), bitphi ? "ДА" : "нет",
                 biteirr ? "ДА" : "нет");
          if (nuse > 0) {
            qsort(rel, (size_t)nuse, sizeof *rel, cmp_dev699);
            printf("   §774 ПОЭЛЕМЕНТНО |Δeirr|/eirr (покрытие %lld из %d, база > 0): медиана "
                   "%.4g, p90 %.4g, p99 %.4g, макс %.4g; eirr < 0 после замыкания: %lld\n",
                   (long long)nuse, cut.nse, rel[nuse / 2], rel[(nuse * 9) / 10],
                   rel[(nuse * 99) / 100], rel[nuse - 1], (long long)nneg2);
          }
          free(rel);
          free(phiB4);
          free(eirrB);
        }
        /* §810: РАСПИСАНИЕ ОТСКОКОВ. База сравнения — УЖЕ посчитанная лестница
         * кадра («2,2»: ординаты канона); расписание — диагностический прогон
         * на своих наборах по sched810_ladder. Кадр остаётся при ординатах
         * канона до приёмки П1/П2 — переворот умолчания отдельным решением. */
        if (xs810_set) {
          if (xs810_1 < 0 || xs810_1 > 4 || xs810_2 < 0 || xs810_2 > 4) {
            fprintf(stderr, "xsched=: поля %d,%d — вне [0, 4] (0 = ордината канона)\n", xs810_1,
                    xs810_2);
            exit(1);
          }
          int n1_810 = xs810_1 > 0 ? xs810_1 : nmu;
          int n2_810 = xs810_2 > 0 ? xs810_2 : nmu;
          tr3_dirs d1_810, d2_810;
          const tr3_dirs *pd1 = &dirs, *pd2 = &dirs;
          if (n1_810 != nmu) {
            if (tr3_dirs_product(&d1_810, n1_810, n1_810) != 0) exit(1);
            pd1 = &d1_810;
          }
          if (n2_810 != nmu) {
            if (tr3_dirs_product(&d2_810, n2_810, n2_810) != 0) exit(1);
            pd2 = &d2_810;
          }
          const size_t nph810 = (size_t)mesh.ncell * 4;
          double *phiS = calloc(nph810, sizeof *phiS);
          double *rel810 = malloc((size_t)(cut.nse > 0 ? cut.nse : 1) * sizeof *rel810);
          if (phiS == NULL || rel810 == NULL) exit(1);
          tr3_stats stS;
          memset(&stS, 0, sizeof stS);
          double qS = -1.0, tS = 0.0;
          sched810_ladder(&prob, pd1, pd2, n1_810, n2_810, xbounce, xtailq, phiS, &stS, &qS, &tS);
          /* сравнение с лестницей кадра (st — закрытое состояние той же сетки) */
          double dps810 = st.psin > 0.0 ? 100.0 * (stS.psin - st.psin) / st.psin : 0.0;
          double zw810 = 0.0, zb810 = 0.0;
          int64_t nzone810 = 0, nuse810 = 0;
          for (int32_t e = 0; e < cut.nse; e++) {
            double a9 = cut.se[e].area;
            if (!(a9 > 0.0) || stS.eirr == NULL || st.eirr == NULL) continue;
            double c9[3] = {0, 0, 0};
            for (int q = 0; q < cut.se[e].nv; q++)
              for (int a = 0; a < 3; a++)
                c9[a] += cut.se[e].v[q][a] / (double)(cut.se[e].nv > 0 ? cut.se[e].nv : 1);
            double d2 = 0.0, dot = 0.0;
            for (int a = 0; a < 3; a++) {
              double dd = c9[a] - g_eye[a];
              d2 += dd * dd;
              dot += cut.se[e].n[a] * dd;
            }
            int inzone = sqrt(d2) < 15.0 && !(dot > 0.0);
            double wS = fabs(stS.eirr[e]) * a9, wB = fabs(st.eirr[e]) * a9;
            if (inzone) {
              zw810 += wS;
              zb810 += wB;
              nzone810++;
            }
            if (fabs(st.eirr[e]) > 0.0)
              rel810[nuse810++] = fabs(stS.eirr[e] - st.eirr[e]) / fabs(st.eirr[e]);
          }
          double dz810 = 100.0 * (zw810 - zb810) / (zb810 > 0.0 ? zb810 : 1.0);
          double p99810 = 0.0;
          if (nuse810 > 0) {
            qsort(rel810, (size_t)nuse810, sizeof *rel810, cmp_dev699);
            p99810 = rel810[(nuse810 * 99) / 100];
          }
          int bit810 =
              n1_810 == nmu && n2_810 == nmu && memcmp(phiS, phi, nph810 * sizeof(double)) == 0;
          printf("   §810 РАСПИСАНИЕ: такт1 nmu %d (ND %d), хвост nmu %d (ND %d); лестница "
                 "%.2f с против %.2f с (выигрыш ×%.2f); Δpsin %+.3f %%, Δ зоны %+.3f %%; "
                 "q̂ %.4f; eirr p99 %.4g (покрытие %lld); тождество «канон,канон»: %s\n",
                 n1_810, pd1->n, n2_810, pd2->n, tS, t774lad, tS > 0.0 ? t774lad / tS : 0.0, dps810,
                 dz810, qS > 0.0 ? qS : 0.0, p99810, (long long)nuse810, bit810 ? "ДА" : "нет");
          free(stS.bout);
          free(stS.sout);
          free(stS.eirr);
          free(phiS);
          free(rel810);
          if (pd1 == &d1_810) tr3_dirs_free(&d1_810);
          if (pd2 == &d2_810) tr3_dirs_free(&d2_810);
        }
        /* §808: СТЕНД ИНКРЕМЕНТНОСТИ СВЕТА. Стоит ПОСЛЕ лестницы и замыкания
         * кадра 1: его вход — закрытое состояние (phi, st.bout, st.sout) и
         * измеренный q̂. Гварды здесь, а не в разборе ключей: условие — не
         * синтаксис, а семантика конфигурации (fcold считается позже). */
        if (xwarm_set) {
          if (!(xcoarse > 0.0)) {
            fprintf(stderr, "xwarm= требует xcoarse>0: ремап грейда определён только "
                            "при огрублении приёмников (план §808)\n");
            exit(1);
          }
          if (!fcold || xhall || xemitfacet) {
            fprintf(stderr, "xwarm= требует канонную инъекцию Ke по элементам (§786): "
                            "несовместим с xfc/xhall/xemitfacet\n");
            exit(1);
          }
          if (xrelax > 0.0) {
            fprintf(stderr, "xwarm= при xrelax>0: демпфер гасил бы bout первой прокидки "
                            "к нулю (§750)\n");
            exit(1);
          }
          if (xwarm_n < 1 || xwarm_n > 2) {
            fprintf(stderr, "xwarm=: прокидок %d — вне [1, 2] (А1348)\n", xwarm_n);
            exit(1);
          }
          warm808 cx8;
          memset(&cx8, 0, sizeof cx8);
          cx8.ot = &ot;
          cx8.ofr = &ofr;
          cx8.cm = &cmap;
          cx8.ft = &ftab;
          cx8.om = &m;
          cx8.ctctx = &CT;
          cx8.occ = &P;
          cx8.smask = solidmask;
          cx8.occn = fr.n;
          cx8.innerfluid = xinnerfluid;
          cx8.lev = lev;
          cx8.xcoarse = xcoarse;
          cx8.xfernosolid = xfernosolid;
          cx8.xmatrho = xmatrho;
          cx8.xmatfar = xmatfar;
          cx8.xrho = xrho;
          cx8.xrhoscale = xrhoscale;
          cx8.xthin = xthin;
          cx8.xsemin = xsemin;
          cx8.xtailq = xtailq;
          cx8.m1 = &mesh;
          cx8.c1 = &cut;
          cx8.phi1 = phi;
          cx8.sout1 = st.sout;
          cx8.qhat1 = qhat;
          cx8.emitpow1 = emitpow2;
          cx8.frho1 = frho; /* Р3 §810: переиспользование альбедо фасетов */
          cx8.xrhocal = xrhocal810;
          cx8.prob = &prob;
          for (int a = 0; a < 3; a++)
            cx8.eye1[a] = g_eye[a];
          cx8.shift[0] = xwarm3[0];
          cx8.shift[1] = xwarm3[1];
          cx8.shift[2] = xwarm3[2];
          cx8.nwarm = xwarm_n;
          cx8.xbounce = xbounce;
          cx8.scramble = xscramble808;
          warm808_run(&cx8);
        }
        free(pphi);
        free(pbo);
        free(pso);
        free(dpphi);
        free(dpbo);
        free(dpso);
        free(peirr);
      } else if (!xunit) {
        src = tr3_sweep_solve(&prob, xit, xtol, phi, &st);
      } else {
        /* §729: ПРИМЕНИТЬ ОПЕРАТОР К ЕДИНИЧНОМУ СОСТОЯНИЮ. Источники в ноль,
         * φ = bout = sout = 1, и `xit` применений по одному такту с прокидкой
         * поверхностного состояния (К76: `bout_in`/`sout_in`). Всё, что вышло
         * `> 1` при альбедо < 1, — локальный усилитель с адресом. На сценах
         * без рассеяния (σ_s = 0) и с `wall_rho = NULL` состояние — только
         * `sout`; `φ` печатается как показание в нормировке `φ/4π` (единичный
         * радианс со всех сторон даёт скалярный поток 4π, §699). */
        if (xconst) {
          fprintf(stderr, "xunit несовместим с xconst: печь задаёт источники, xunit их обнуляет\n");
          exit(1);
        }
        prob.elem_emit = NULL;
        prob.facet_emit = NULL;
        prob.warm_start = 1;
        for (int32_t ci = 0; ci < mesh.ncell; ci++) {
          phi[4 * (size_t)ci] = 1.0;
          phi[4 * (size_t)ci + 1] = phi[4 * (size_t)ci + 2] = phi[4 * (size_t)ci + 3] = 0.0;
        }
        double *ub = calloc((size_t)mesh.nf * 4, sizeof *ub);
        double *us = calloc((size_t)(cut.nse > 0 ? cut.nse : 1) * 4, sizeof *us);
        double *scr9 = malloc((size_t)(mesh.ncell > cut.nse ? mesh.ncell : cut.nse) * sizeof *scr9);
        if (ub == NULL || us == NULL || scr9 == NULL) exit(1);
        for (int32_t f = 0; f < mesh.nf; f++)
          ub[4 * (size_t)f] = 1.0;
        for (int32_t k = 0; k < cut.nse; k++)
          us[4 * (size_t)k] = 1.0;
        const double fourpi = 4.0 * 3.14159265358979323846;
        for (int ap = 1; ap <= xit; ap++) {
          prob.bout_in = ub;
          prob.sout_in = us;
          memset(&st, 0, sizeof st);
          src = tr3_sweep_solve(&prob, 1, 0.0, phi, &st);
          if (src != 0) break;
          /* элементы: население — ненулевой нулевой коэффициент; нуль здесь
           * означает «не участвует» (за фильтрами xthin/xsemin или без влёта),
           * а не «усиление нулевое» — класс А1097 */
          int64_t nz = 0, ngt1 = 0;
          double emax = 0.0;
          int32_t iemax = -1;
          for (int32_t k = 0; k < cut.nse; k++) {
            double v = fabs(st.sout[4 * (size_t)k]);
            if (!(v > 0.0)) continue;
            scr9[nz++] = v;
            if (v > 1.0) ngt1++;
            if (v > emax) {
              emax = v;
              iemax = k;
            }
          }
          if (nz > 0) {
            qsort(scr9, (size_t)nz, sizeof *scr9, cmp_dev699);
            printf("   §729 ПРИМЕНЕНИЕ %d, ЭЛЕМЕНТЫ (нулевой коэфф.): ненулевых %lld из %d, "
                   "медиана %.4g, p99 %.4g, МАКС %.4g (элемент %d); > 1: %lld (%.2f %%)\n",
                   ap, (long long)nz, cut.nse, scr9[nz / 2], scr9[(nz * 99) / 100], emax, iemax,
                   (long long)ngt1, 100.0 * (double)ngt1 / (double)nz);
          }
          /* ячейки: φ/4π по флюиду; нули включаются — замурованные полости
           * законны (§699). Σφ и max φ в %.17g — слепок проводки П5(б):
           * такт 1 обязан совпасть побитово между xrho=0.7 и xrho=0. */
          int64_t nfl9 = 0, cgt1 = 0;
          double cmax = 0.0, csum = 0.0, phimax = 0.0;
          int32_t icmax = -1;
          for (int32_t ci = 0; ci < mesh.ncell; ci++) {
            if (!(cut.mvol[ci][0][0] > 0.0)) continue;
            double v = phi[4 * (size_t)ci] / fourpi;
            scr9[nfl9++] = fabs(v);
            csum += fabs(phi[4 * (size_t)ci]);
            if (fabs(phi[4 * (size_t)ci]) > phimax) phimax = fabs(phi[4 * (size_t)ci]);
            if (v > 1.0) cgt1++;
            if (fabs(v) > cmax) {
              cmax = fabs(v);
              icmax = ci;
            }
          }
          if (nfl9 > 0) {
            qsort(scr9, (size_t)nfl9, sizeof *scr9, cmp_dev699);
            printf("   §729   ЯЧЕЙКИ φ/4π: флюидных %lld, медиана %.4g, p99 %.4g, МАКС %.4g "
                   "(ячейка %d); > 1: %lld (%.2f %%); слепок: Σ|φ| %.17g, max|φ| %.17g\n",
                   (long long)nfl9, scr9[nfl9 / 2], scr9[(nfl9 * 99) / 100], cmax, icmax,
                   (long long)cgt1, 100.0 * (double)cgt1 / (double)nfl9, csum, phimax);
          }
          double bmax = 0.0;
          int32_t ibmax = -1;
          for (int32_t f = 0; f < mesh.nf; f++)
            if (fabs(st.bout[4 * (size_t)f]) > bmax) {
              bmax = fabs(st.bout[4 * (size_t)f]);
              ibmax = f;
            }
          printf("   §729   ГРАНИ bout (канал инертен при wall_rho = NULL): max %.4g (грань "
                 "%d)\n",
                 bmax, ibmax);
          /* прокидка состояния: владение прежним входом снимается, новым — у
           * последнего `st`, и его освобождает общий путь ниже */
          free(ub);
          free(us);
          ub = st.bout;
          us = st.sout;
        }
        free(scr9);
        if (src != 0) {
          free(ub);
          free(us);
        }
      }
      tsw = now_s() - tsw;
      if (xconst) {
        /* §699: ОТКЛОНЕНИЕ ОТ ТОЧНОГО РЕШЕНИЯ. Считается по ФЛЮИДНЫМ ячейкам:
         * в сплошных поля нет по построению, и включать их значило бы мерить
         * пустоту. Печатается и максимум наклонов — печь проверяет только
         * среднее (К12), а наклоны суть две трети неизвестных. */
        double exact4pi = 4.0 * 3.14159265358979323846 * bconst;
        double worst = 0.0, wslope = 0.0;
        int32_t wc = -1;
        int64_t nfl = 0;
        for (int32_t ci = 0; ci < mesh.ncell; ci++) {
          if (!(cut.mvol[ci][0][0] > 0.0)) continue;
          nfl++;
          double e = fabs(phi[4 * (size_t)ci] - exact4pi) / exact4pi;
          if (e > worst) {
            worst = e;
            wc = ci;
          }
          for (int j = 1; j < 4; j++) {
            double sl = fabs(phi[4 * (size_t)ci + (size_t)j]) / exact4pi;
            if (sl > wslope) wslope = sl;
          }
        }
        double wfrac = -1.0;
        if (wc >= 0) {
          double s3 = (double)mesh.csize[wc];
          double V = s3 * s3 * s3 * ofr.u[0] * ofr.u[1] * ofr.u[2];
          wfrac = V > 0.0 ? cut.mvol[wc][0][0] / V : -1.0;
        }
        /* §699 ПОПРАВКА ПРИБОРА ПО ХОДУ: максимум сам по себе неоднозначен.
         * Отклонение `1.000` значит `φ = 0`, а нулевая ячейка может быть
         * ЗАМУРОВАННОЙ полостью — туда свету взяться неоткуда, и ноль там
         * ЗАКОНЕН. Различить можно только распределением: если тождество держится
         * везде, кроме отрезанных карманов, медиана обязана быть нулевой. */
        {
          double *dv = malloc((size_t)(nfl > 0 ? nfl : 1) * sizeof *dv);
          if (dv != NULL) {
            int64_t nd2 = 0, nzero = 0;
            for (int32_t ci = 0; ci < mesh.ncell; ci++) {
              if (!(cut.mvol[ci][0][0] > 0.0)) continue;
              double v = phi[4 * (size_t)ci];
              if (!(v > 0.0)) nzero++;
              dv[nd2++] = fabs(v - exact4pi) / exact4pi;
            }
            qsort(dv, (size_t)nd2, sizeof *dv, cmp_dev699);
            printf("   §699 РАСПРЕДЕЛЕНИЕ ОТКЛОНЕНИЯ: медиана %.3e, p90 %.3e, p99 %.3e, макс "
                   "%.3e; ячеек РОВНО В НУЛЕ %lld из %lld (%.2f %%)\n",
                   dv[nd2 / 2], dv[(nd2 * 9) / 10], dv[(nd2 * 99) / 100], dv[nd2 - 1],
                   (long long)nzero, (long long)nd2, 100.0 * (double)nzero / (double)nd2);
            free(dv);
          }
        }
        printf("   §699 ПЕЧЬ НА РАЗРЕЗЕ, ИСХОД: флюидных ячеек %lld; МАКС ОТН. ОТКЛОНЕНИЕ %.3e "
               "(ячейка %d, доля флюида %.4g, размер %d); макс наклон %.3e; итераций %d, "
               "невязка %.3e, nclip %d\n",
               (long long)nfl, worst, wc, wfrac, wc >= 0 ? mesh.csize[wc] : -1, wslope, st.iters,
               st.resid, st.nclip);
      }
      /* ГДЕ максимум — в комнате или снаружи. Без этого «φ = 1.9e5» неотличимо
       * от «φ велико в ячейке-щепке вне сцены» (правило А807: печатать вход
       * подозреваемой стадии, а не только её выход). */
      double phimax = 0.0, phisum = 0.0, soutmax = 0.0, phimax_in = 0.0, phimax_out = 0.0;
      double volmin_at_max = 0.0;
      int32_t maxcell = -1;
      for (int32_t ci = 0; ci < mesh.ncell; ci++) {
        double p0 = phi[4 * (size_t)ci];
        size_t kk = hz_occ_index(fr.n, mesh.clo[ci][0], mesh.clo[ci][1], mesh.clo[ci][2]);
        int inner = solidmask != NULL && solidmask[kk] == 1u;
        if (p0 > phimax) {
          phimax = p0;
          volmin_at_max = cut.mvol[ci][0][0];
          maxcell = ci;
        }
        if (inner && p0 > phimax_in) phimax_in = p0;
        if (!inner && p0 > phimax_out) phimax_out = p0;
        phisum += p0;
      }
      if (st.sout != NULL)
        for (int32_t k = 0; k < cut.nse; k++)
          if (st.sout[4 * (size_t)k] > soutmax) soutmax = st.sout[4 * (size_t)k];
      /* §489: ПОЧЕМУ У ЧАСТИ ЯЧЕЕК СРЕЗА НЕТ РАДИАНСА. Два множества по одной
       * сетке: A — ячейки с ВЕРШИНОЙ DC (из них собирается срез), B — ячейки с
       * хотя бы одним ПОВЕРХНОСТНЫМ ЭЛЕМЕНТОМ. Предикаты разные (вершина против
       * «плоскость режет коробку»), и совпадать они не обязаны. Считается
       * прямо, а не оценивается. */
      {
        hz_dcslice SF;
        int64_t na = 0, nb = 0, nab = 0;
        if (hz_slice_init(&SF, lev) != HZ_DC_OK) exit(1);
        if (hz_slice_build(&SF, &T, &ht, NULL, NULL) != HZ_DC_OK) exit(1);
        unsigned char *hasel = calloc((size_t)fr.n * (size_t)fr.n * (size_t)fr.n, 1);
        if (hasel == NULL) exit(1);
        for (int32_t k = 0; k < cut.nse; k++) {
          int32_t ci2 = cut.se[k].cell;
          hasel[hz_occ_index(fr.n, mesh.clo[ci2][0], mesh.clo[ci2][1], mesh.clo[ci2][2])] = 1u;
        }
        for (size_t i2 = 0; i2 < (size_t)fr.n * (size_t)fr.n * (size_t)fr.n; i2++)
          if (hasel[i2]) nb++;
        for (int32_t i2 = 0; i2 < SF.n; i2++) {
          size_t g = hz_occ_index(fr.n, SF.c[i2].lo[0], SF.c[i2].lo[1], SF.c[i2].lo[2]);
          na++;
          if (hasel[g]) nab++;
        }
        /* РЕШАЮЩЕЕ ЧИСЛО (§490): сколько ячеек СЕТКИ ПЕРЕНОСА вообще нашли свою
         * запись в боковой таблице. Ключом там служит индекс узла, и деревьев
         * ДВА — DC и hz_octree, — а нумерация у них своя. Если совпадений
         * заметно меньше числа записей, ключи разные, и фасеты розданы не тем
         * ячейкам. */
        int64_t nrecfound = 0;
        for (int32_t c2 = 0; c2 < mesh.ncell; c2++)
          if (hz_cutmap_find(&cmap, mesh.node[c2]) != NULL) nrecfound++;
        printf("      §490 КЛЮЧ БОКОВОЙ ТАБЛИЦЫ: записей %d, ячеек сетки переноса, нашедших "
               "запись, %lld (%.1f %%)\n",
               cmap.nr, (long long)nrecfound,
               100.0 * (double)nrecfound / (double)(cmap.nr ? cmap.nr : 1));
        printf("      §489 МНОЖЕСТВА: с вершиной DC %lld, с поверхностным элементом %lld, "
               "пересечение %lld (%.1f %% от вершин); вершин без элемента %lld, элементов без "
               "вершины %lld\n",
               (long long)na, (long long)nb, (long long)nab,
               100.0 * (double)nab / (double)(na ? na : 1), (long long)(na - nab),
               (long long)(nb - nab));
        free(hasel);
        hz_slice_free(&SF);
      }
      printf("   РАЗВЁРТКА: направлений %d, светящихся элементов %lld; код %d, итераций %d, "
             "невязка %.2e (rel %.1e, max|состояние| %.3e), за %.2f с\n",
             dirs.n, (long long)nlit, src, st.iters, st.resid, st.resid_rel, st.statemax, tsw);
      printf("      ЭНЕРГИЯ: втекло %.4e, вытекло %.4e, поглощено %.4e, баланс %.2e; в "
             "поверхности %.4e, из них %.4e; max φ %.4e, max исходящий радианс %.4e\n",
             st.pin, st.pout, st.pabs, st.balance, st.psin, st.psout, phimax, soutmax);
      /* ТОЖДЕСТВО К40 ПЕЧАТАЕТСЯ ЯВНО И С ДОЛЕЙ (§674). До этого «баланс» стоял
       * одним числом без знаменателя, и прочесть его было нельзя: `5.61` — это
       * много или мало? Теперь видно, что это `76 %` излучённого.
       *     pin + psout = pout + pabs + psin
       * Недостача считается ЗДЕСЬ, а не берётся из `st.balance`: так видно, из
       * каких именно членов она сложилась. */
      {
        double lack = st.pin + st.psout - st.pout - st.pabs - st.psin - st.psolid;
        double den = st.psout > 0.0 ? st.psout : 1.0;
        printf("      БАЛАНС К40: втекло %.4e + отдано поверхностями %.4e = вытекло %.4e + "
               "поглощено объёмом %.4e + упало на поверхности %.4e + УШЛО В СПЛОШНОЕ %.4e; "
               "НЕДОСТАЧА %.4e = %.1f %% отданного\n",
               st.pin, st.psout, st.pout, st.pabs, st.psin, st.psolid, lack, 100.0 * lack / den);
      }
      /* §739: СВЕТО-ВЗВЕШЕННОЕ ЭФФЕКТИВНОЕ АЛЬБЕДО. psout = отражение +
       * собственное излучение (К40), излучение = Σπ·L·площадь (emitpow2).
       * Прибор проверен на двух известных точках: 0.6969 при xrho=0.7 и
       * 0.9973 при xrho=1 (ручной счёт по логам §737). При xhall излучение
       * задано иначе, и формула не действует. */
      /* §796: под first-collision «эмиссия» формулы — ИНЪЕЦИРОВАННАЯ мощность
       * Σρ·E_fc·area, а не π·L·area ламп: инъецируется первое отражение. */
      double emitref796 = fcold ? emitpow2 : injpow796;
      if (st.psin > 0.0 && !xhall)
        printf("      §739 ρ_eff = (psout − эмиссия)/psin = (%.4e − %.4e)/%.4e = %.4f\n", st.psout,
               emitref796, st.psin, (st.psout - emitref796) / st.psin);
      /* §782: Σ|E·area| БЛИЖНЕЙ ЛИЦЕВОЙ ЗОНЫ — приёмочная величина политики
       * огрубления (тот же предикат, что §768/§774): дальняя зона грубится,
       * и её ошибка не должна доносить сюда (А1212). Печатается всегда. */
      if (st.eirr != NULL) {
        double zs2 = 0.0, ts2 = 0.0;
        int64_t nz2 = 0;
        for (int32_t e = 0; e < cut.nse; e++) {
          double c9[3] = {0, 0, 0};
          for (int q2 = 0; q2 < cut.se[e].nv; q2++)
            for (int a = 0; a < 3; a++)
              c9[a] += cut.se[e].v[q2][a] / (double)(cut.se[e].nv > 0 ? cut.se[e].nv : 1);
          double d2 = 0.0, dot = 0.0;
          for (int a = 0; a < 3; a++) {
            double dd = c9[a] - g_eye[a];
            d2 += dd * dd;
            dot += cut.se[e].n[a] * dd;
          }
          /* §796: зона меряет ПОЛНУЮ облучённость E_fc + рассеянное; при
           * старой инъекции — прежнее eirr (ветвь, не «+0», ради посимвольной
           * воспроизводимости мира без ключа). */
          double e796 = efc796 != NULL ? st.eirr[e] + efc796[e] : st.eirr[e];
          double w = e796 * cut.se[e].area;
          ts2 += fabs(w);
          if (cut.se[e].nv > 0 && sqrt(d2) < 15.0 && !(dot > 0.0)) {
            zs2 += fabs(w);
            nz2++;
          }
        }
        printf("      §782 ЗОНА (<15 м, лицевые): элементов %lld, Σ|E·area| %.6g; всей сцены "
               "%.6g\n",
               (long long)nz2, zs2, ts2);
      }
      /* ---- §796: ПРИБОР ЧЕШУИ ЧИСЛОМ (только в кадровых прогонах) --------
       * Чешуя = вариация полной облучённости между элементами ОДНОЙ плоской
       * поверхности. Группировка по КВАНТОВАННОЙ ПЛОСКОСТИ, а не по фасету:
       * DC-фасеты одноклеточные, и фасетные группы мелки (А1288 — населённость
       * фасетов печатается рядом как замер этого). Ключ группы: нормаль,
       * округлённая к решётке 1/HZ_CV_NQ (4 — разделяет оси и диагонали, но
       * терпит QEF-рябь), плюс корзина офсета шагом HZ_CV_OFFH·h (4h — толще
       * кванта DC-вершины, тоньше межэтажного шага). CV = std/mean по группе
       * из ≥ HZ_CV_MIN элементов (8 — минимум, при котором std читаем). */
      if (xframe && st.eirr != NULL && cut.nse > 0) {
        enum { HZ_CV_NQ = 4, HZ_CV_OFFH = 4, HZ_CV_MIN = 8 };
        cvpair796 *cp = malloc((size_t)cut.nse * sizeof *cp);
        int32_t *fpop = calloc((size_t)(ftab.n > 0 ? ftab.n : 1), sizeof *fpop);
        if (cp == NULL || fpop == NULL) exit(1);
        int32_t ncp = 0;
        for (int32_t e = 0; e < cut.nse; e++) {
          if (!(cut.se[e].area > 0.0)) continue;
          if (cut.se[e].facet >= 0 && cut.se[e].facet < ftab.n) fpop[cut.se[e].facet]++;
          double pc[3];
          efc_center(&cut, &mesh, &fr, e, pc);
          double off = 0.0;
          uint64_t kq = 0;
          for (int a = 0; a < 3; a++) {
            long qn = lround(cut.se[e].n[a] * (double)HZ_CV_NQ);
            kq = (kq << 8) | (uint64_t)(uint8_t)(qn + 16);
            off += cut.se[e].n[a] * pc[a];
          }
          long ob = (long)floor(off / ((double)HZ_CV_OFFH * fr.h));
          cp[ncp].key = (kq << 32) | (uint64_t)(uint32_t)(int32_t)ob;
          cp[ncp].e = e;
          ncp++;
        }
        qsort(cp, (size_t)ncp, sizeof *cp, cmp_cvpair796);
        double *cvv = malloc((size_t)(ncp > 0 ? ncp : 1) * sizeof *cvv);
        double *cvf6 = malloc((size_t)(ncp > 0 ? ncp : 1) * sizeof *cvf6);
        double *cvt6 = malloc((size_t)(ncp > 0 ? ncp : 1) * sizeof *cvt6);
        if (cvv == NULL || cvf6 == NULL || cvt6 == NULL) exit(1);
        int32_t ng = 0, ngbig = 0, nf6 = 0, nt6 = 0;
        int32_t i0 = 0;
        while (i0 < ncp) {
          int32_t i1 = i0;
          while (i1 < ncp && cp[i1].key == cp[i0].key)
            i1++;
          ng++;
          if (i1 - i0 >= HZ_CV_MIN) {
            /* три канала: полный свет, прямой E_fc, рассеянный хвост —
             * компонентная атрибуция ряби (какое поле пятнисто) */
            double s1 = 0.0, s2 = 0.0, f1 = 0.0, f2 = 0.0, t1 = 0.0, t2 = 0.0;
            for (int32_t j = i0; j < i1; j++) {
              int32_t e = cp[j].e;
              double vf = efc796 != NULL ? efc796[e] : 0.0;
              double vt = st.eirr[e];
              double v = vf + vt;
              s1 += v;
              s2 += v * v;
              f1 += vf;
              f2 += vf * vf;
              t1 += vt;
              t2 += vt * vt;
            }
            double nn6 = (double)(i1 - i0);
            double mean = s1 / nn6;
            double var = s2 / nn6 - mean * mean;
            if (mean > 0.0 && var > 0.0) cvv[ngbig++] = sqrt(var) / mean;
            if (efc796 != NULL) {
              double mf = f1 / nn6, vfv = f2 / nn6 - mf * mf;
              double mt = t1 / nn6, vtv = t2 / nn6 - mt * mt;
              if (mf > 0.0 && vfv > 0.0) cvf6[nf6++] = sqrt(vfv) / mf;
              if (mt > 0.0 && vtv > 0.0) cvt6[nt6++] = sqrt(vtv) / mt;
            }
          }
          i0 = i1;
        }
        /* населённость фасетов — замер А1288 (сила фасетной группировки) */
        int32_t nfnz = 0;
        for (int32_t i = 0; i < ftab.n; i++)
          if (fpop[i] > 0) {
            fpop[nfnz++] = fpop[i]; /* уплотнение на месте: fpop дальше не ключ */
          }
        double fp50 = 0.0, fp90 = 0.0;
        int32_t fpmax = 0;
        if (nfnz > 0) {
          qsort(fpop, (size_t)nfnz, sizeof *fpop, cmp_i32_796);
          fp50 = fpop[nfnz / 2];
          fp90 = fpop[(int32_t)((int64_t)nfnz * 9 / 10)];
          fpmax = fpop[nfnz - 1];
        }
        if (ngbig > 0) {
          qsort(cvv, (size_t)ngbig, sizeof *cvv, cmp_dev699);
          printf("   §796 ЧЕШУЯ (%s): групп плоскостей %d (с ≥%d элементами %d); CV p50 %.4f, "
                 "p90 %.4f, макс %.4f; населённость фасетов: p50 %.0f, p90 %.0f, макс %d\n",
                 efc796 != NULL ? "fc" : "без fc", ng, (int)HZ_CV_MIN, ngbig, cvv[ngbig / 2],
                 cvv[(int32_t)((int64_t)ngbig * 9 / 10)], cvv[ngbig - 1], fp50, fp90, fpmax);
          if (nf6 > 0 && nt6 > 0) {
            qsort(cvf6, (size_t)nf6, sizeof *cvf6, cmp_dev699);
            qsort(cvt6, (size_t)nt6, sizeof *cvt6, cmp_dev699);
            printf("   §796 ЧЕШУЯ ПО КОМПОНЕНТАМ: ПРЯМОЙ E_fc CV p50 %.4f, p90 %.4f (групп %d); "
                   "ХВОСТ свипа CV p50 %.4f, p90 %.4f (групп %d)\n",
                   cvf6[nf6 / 2], cvf6[(int32_t)((int64_t)nf6 * 9 / 10)], nf6, cvt6[nt6 / 2],
                   cvt6[(int32_t)((int64_t)nt6 * 9 / 10)], nt6);
          }
        } else
          printf("   §796 ЧЕШУЯ: групп с ≥%d элементами НЕТ (плоскостных групп %d) — прибор "
                 "пуст, приёмка только глазами\n",
                 (int)HZ_CV_MIN, ng);
        free(cvv);
        free(cvf6);
        free(cvt6);
        free(cp);
        free(fpop);
      }
      /* §768: ПРИБОР ВКЛАДА — перевозмущение. Выключаем отражение класса
       * поверхностей (ρ→0 поэлементно; заслон сохраняется — А1210) и меряем
       * изменение поля E в БЛИЖНЕЙ ЛИЦЕВОЙ зоне (прокси видимого: < 15 м от
       * камеры, нормаль к камере; заслонённость не учитывается — А1209,
       * смещение консервативно). Эмиссия не выключается (А1211). */
      if (xcontrib && st.eirr != NULL) {
        double *rbase = malloc((size_t)(cut.nse > 0 ? cut.nse : 1) * sizeof *rbase);
        double *rkill = malloc((size_t)(cut.nse > 0 ? cut.nse : 1) * sizeof *rkill);
        double *phic = malloc((size_t)mesh.ncell * 4 * sizeof *phic);
        double (*ecen)[3] = malloc((size_t)(cut.nse > 0 ? cut.nse : 1) * sizeof *ecen);
        double *edist = malloc((size_t)(cut.nse > 0 ? cut.nse : 1) * sizeof *edist);
        unsigned char *zone = malloc((size_t)(cut.nse > 0 ? cut.nse : 1));
        unsigned char *eback = malloc((size_t)(cut.nse > 0 ? cut.nse : 1));
        if (rbase == NULL || rkill == NULL || phic == NULL || ecen == NULL || edist == NULL ||
            zone == NULL || eback == NULL)
          exit(1);
        for (int32_t e = 0; e < cut.nse; e++) {
          rbase[e] = (frho != NULL && cut.se[e].facet < ftab.n) ? frho[cut.se[e].facet] : 0.0;
          double c9[3] = {0, 0, 0};
          for (int q2 = 0; q2 < cut.se[e].nv; q2++)
            for (int a = 0; a < 3; a++)
              c9[a] += cut.se[e].v[q2][a] / (double)(cut.se[e].nv > 0 ? cut.se[e].nv : 1);
          double d2 = 0.0, dot = 0.0;
          for (int a = 0; a < 3; a++) {
            ecen[e][a] = c9[a];
            double dd = c9[a] - g_eye[a];
            d2 += dd * dd;
            dot += cut.se[e].n[a] * dd; /* > 0: нормаль ОТ камеры (обратная) */
          }
          edist[e] = sqrt(d2);
          eback[e] = dot > 0.0 ? 1 : 0;
          /* зона: ближняя (< 15 м — четверть габарита сцены) и лицевая */
          zone[e] = (edist[e] < 15.0 && !eback[e]) ? 1 : 0;
        }
        double zden = 0.0, tden = 0.0;
        int64_t nzone = 0;
        for (int32_t e = 0; e < cut.nse; e++) {
          double w = st.eirr[e] * cut.se[e].area;
          tden += fabs(w);
          if (zone[e]) {
            zden += fabs(w);
            nzone++;
          }
        }
        printf("   §768 ЗОНА (< 15 м, лицевые): элементов %lld из %d; Σ|E·area| зоны %.4g, всей "
               "сцены %.4g\n",
               (long long)nzone, cut.nse, zden, tden);
        static const double KD[6] = {40.0, 20.0, 10.0, -1.0 /* back */, -2.0 /* inf НК */, 0};
        for (int kv = 0; kv < 5; kv++) {
          int64_t nkill = 0;
          for (int32_t e = 0; e < cut.nse; e++) {
            int kill = 0;
            if (KD[kv] > 0.0)
              kill = edist[e] > KD[kv];
            else if (KD[kv] < -1.5)
              kill = 0; /* НК: никого */
            else
              kill = eback[e]; /* обратные */
            rkill[e] = kill ? 0.0 : rbase[e];
            if (kill) nkill++;
          }
          tr3_problem pk = prob;
          pk.elem_rho = rkill;
          memset(phic, 0, (size_t)mesh.ncell * 4 * sizeof *phic);
          tr3_stats stk;
          memset(&stk, 0, sizeof stk);
          if (tr3_sweep_solve(&pk, xit, xtol, phic, &stk) != 0) exit(1);
          double zsum = 0.0, tsum = 0.0;
          int bit = 1;
          if (stk.eirr != NULL)
            for (int32_t e = 0; e < cut.nse; e++) {
              double dw = (stk.eirr[e] - st.eirr[e]) * cut.se[e].area;
              tsum += fabs(dw);
              if (zone[e]) zsum += fabs(dw);
              if (memcmp(&stk.eirr[e], &st.eirr[e], sizeof(double)) != 0) bit = 0;
            }
          const char *nm = KD[kv] > 0.0 ? "ДАЛЬШЕ" : (KD[kv] < -1.5 ? "НИКОГО (НК)" : "ОБРАТНЫЕ");
          printf("   §768 KILL %s %s%.0f м: выключено %lld элементов; Δ в ближней лицевой зоне "
                 "%.4g (%.2f %% от Σ зоны); Δ всюду %.2f %%; побитово с базой: %s\n",
                 nm, KD[kv] > 0.0 ? "> " : "", KD[kv] > 0.0 ? KD[kv] : 0.0, (long long)nkill, zsum,
                 100.0 * zsum / (zden > 0.0 ? zden : 1.0), 100.0 * tsum / (tden > 0.0 ? tden : 1.0),
                 bit ? "ДА" : "нет");
          free(stk.bout);
          free(stk.sout);
          free(stk.eirr);
        }
        free(rbase);
        free(rkill);
        free(phic);
        free(ecen);
        free(edist);
        free(zone);
        free(eback);
      }
      /* §758: свиповое поле на узлы — для кадра развёрткой. E = eirr
       * (входящая облучённость финального такта, А1187), вес — площадь. */
      if (xframe && st.eirr != NULL) {
        g_swEn = calloc((size_t)T.n, sizeof *g_swEn);
        g_swEd = calloc((size_t)T.n, sizeof *g_swEd);
        if (g_swEn == NULL || g_swEd == NULL) exit(1);
        int64_t nmiss758 = 0;
        for (int32_t e = 0; e < cut.nse; e++) {
          double a = cut.se[e].area;
          int32_t ci = cut.se[e].cell;
          if (!(a > 0.0) || ci < 0 || ci >= mesh.ncell) continue;
          int32_t szc = mesh.csize[ci];
          int lvl8 = lev;
          while (szc > 1) {
            szc >>= 1;
            lvl8--;
          }
          hz_dccell q8;
          memset(&q8, 0, sizeof q8);
          q8.lvl = (uint8_t)lvl8;
          for (int a2 = 0; a2 < 3; a2++)
            q8.lo[a2] = (uint16_t)mesh.clo[ci][a2];
          int32_t ni = node_of_cell(&T, lev, &q8);
          if (ni < 0 || ni >= T.n) {
            nmiss758++;
            continue;
          }
          /* §796/§798: под гибридом узлы свипа несут ТОЛЬКО хвост (прямой
           * идёт вторым каналом из ядрового сбора); под xfcelem (НК §798) —
           * §796-агрегат E_fc + хвост; при старой инъекции — прежнее eirr */
          double ev8 = st.eirr[e];
          if (efc796 != NULL && xfcelem) ev8 += efc796[e];
          g_swEn[ni] += ev8 * a;
          g_swEd[ni] += a;
        }
        g_fcmode = efc796 == NULL ? 0 : (xfcelem ? 1 : 2);
        swE_lift(&T, 0, g_swEn, g_swEd);
        printf("   §758 СВИП->УЗЛЫ: элементов без узла %lld\n", (long long)nmiss758);
      }
      /* ---- §798: ПРЯМОЙ КАНАЛ КАДРА ИЗ ЯДРОВОГО СБОРА (гибрид К-а §797) --
       * Приёмники — ячейки СРЕЗА ПОЛНОЙ ГЛУБИНЫ (одна DC-нормаль на ячейку,
       * двусторонний агрегат А1192 не возникает); b0 — скалярное среднее Ke
       * материала (монохромная семантика свипа §667); сбор hgather_rec с
       * заслонами, БЕЗ альбедо приёмника — узлы получают ОБЛУЧЁННОСТЬ E_dir.
       * Кадр читает сумму двух каналов. Перенос не тронут. */
      if (xframe && st.eirr != NULL && efc796 != NULL && !xfcelem) {
        double tdc0 = now_s();
        hz_dcslice SD;
        if (hz_slice_init(&SD, lev) != HZ_DC_OK) exit(1);
        if (hz_slice_build(&SD, &T, &ht, NULL, NULL) != HZ_DC_OK) exit(1);
        int64_t nmatd = slice_assign_mat(&SD, &CT, &fr, &m);
        float *b0d = calloc(3 * (size_t)(SD.n > 0 ? SD.n : 1), sizeof *b0d);
        double *edir8 = calloc((size_t)(SD.n > 0 ? SD.n : 1), sizeof *edir8);
        if (b0d == NULL || edir8 == NULL) exit(1);
        int64_t nlitd = 0;
        for (int32_t i = 0; i < SD.n; i++) {
          const double *ke = m.mtl[SD.c[i].mat < m.nmtl ? SD.c[i].mat : 0].ke3;
          double kem = (ke[0] + ke[1] + ke[2]) / 3.0;
          if (!(kem > 0.0)) continue;
          for (int k = 0; k < 3; k++)
            b0d[3 * (size_t)i + (size_t)k] = (float)kem;
          nlitd++;
        }
        etree ETD;
        memset(&ETD, 0, sizeof ETD);
        etree_build(&ETD, &SD, &fr, b0d, &m, 0, SD.n, 0, lev);
        /* порог вклада — как у fc-сбора §796 (канон ядра, hcontrib= перебивает) */
        double hcd = hcontrib_set ? g_hcontrib : HZ_FC_CONTRIB;
        double bsumd = 0.0, asumd = 0.0;
        for (int k = 0; k < ETD.e[0].nb; k++) {
          const ebin *rb = &ETD.b[ETD.e[0].b0 + k];
          asumd += (double)rb->area;
          for (int c = 0; c < 3; c++)
            bsumd += (double)rb->flux[c] / 3.0;
        }
        double taud = hcd * (asumd > 0.0 ? bsumd / asumd : 0.0);
        int nthd = g_omp1 ? 1 : omp_get_max_threads();
        int64_t *pld = calloc((size_t)nthd, sizeof *pld);
        double *ptd = calloc((size_t)nthd, sizeof *ptd);
        double *pad = calloc((size_t)nthd, sizeof *pad);
        int64_t *pmd = calloc((size_t)nthd, sizeof *pmd);
        if (pld == NULL || ptd == NULL || pad == NULL || pmd == NULL) exit(1);
#pragma omp parallel for schedule(dynamic, 64) if (!g_omp1)
        for (int32_t i = 0; i < SD.n; i++) {
          int th = g_omp1 ? 0 : omp_get_thread_num();
          double pi[3] = {0, 0, 0}, nn2[3] = {0, 0, 0};
          hz_slice_vertex(&SD, i, pi);
          for (int k = 0; k < 3; k++)
            pi[k] = fr.org[k] + pi[k] * fr.h;
          hz_slice_normal(&SD, i, nn2);
          double rrecv = fr.h * (double)((int32_t)1 << (lev - (int)SD.c[i].lvl));
          double acc[3] = {0, 0, 0}, ffacc = 0.0;
          for (int q2 = 0; q2 < 6; q2++)
            hgather_rec(&ETD, 0, q2, pi, nn2, HZ_FC_EPS, taud, rrecv, &P, &fr, 1, acc, &pld[th],
                        &ptd[th], &pad[th], NULL, &ffacc, NULL, 0, &pmd[th]);
          double fnorm = (!g_gnonorm && ffacc > 1.0) ? 1.0 / ffacc : 1.0;
          edir8[i] = acc[0] * fnorm; /* каналы равны: b0 скалярный */
        }
        int64_t nld = 0, nmd = 0;
        double alld = 0.0, thrd = 0.0;
        for (int t = 0; t < nthd; t++) {
          nld += pld[t];
          alld += pad[t];
          thrd += ptd[t];
          nmd += pmd[t];
        }
        free(pld);
        free(ptd);
        free(pad);
        free(pmd);
        etree_free(&ETD);
        free(b0d);
        /* раскладка на узлы: срез полной глубины — ячейка ↔ узел 1:1; отказы
         * спуска считаются (А1029/А1305), непокрытый узел ≠ тёмный ноль */
        g_fcdn = calloc((size_t)T.n, sizeof *g_fcdn);
        g_fcdd = calloc((size_t)T.n, sizeof *g_fcdd);
        if (g_fcdn == NULL || g_fcdd == NULL) exit(1);
        int64_t nputd = 0, nbadd = 0, nzd = 0;
        double sed8 = 0.0;
        for (int32_t i = 0; i < SD.n; i++) {
          int32_t ni = node_of_cell(&T, lev, &SD.c[i]);
          if (ni < 0 || ni >= T.n) {
            nbadd++;
            continue;
          }
          g_fcdn[ni] += edir8[i];
          g_fcdd[ni] += 1.0;
          nputd++;
          sed8 += edir8[i];
          if (!(edir8[i] > 0.0)) nzd++;
        }
        swE_lift(&T, 0, g_fcdn, g_fcdd);
        printf("   §798 ПРЯМОЙ КАНАЛ (ядровый сбор): ячеек среза %d (материал %lld, светящихся "
               "%lld), связей %lld (заслонённого потока %.1f %%), маршей %lld; узлов заполнено "
               "%lld, отказов %lld; Σ E_dir %.6e, E_dir = 0 у %lld (%.1f %%); %.2f с\n",
               SD.n, (long long)nmatd, (long long)nlitd, (long long)nld,
               100.0 * thrd / (alld > 0.0 ? alld : 1.0), (long long)nmd, (long long)nputd,
               (long long)nbadd, sed8, (long long)nzd,
               100.0 * (double)nzd / (double)(nputd > 0 ? nputd : 1), now_s() - tdc0);
        free(edir8);
        hz_slice_free(&SD);
      }
      /* §744: агрегаты по ячейкам для сличения с ядром §597 — площадно-
       * взвешенная косвенная облучённость (полный − прямой) и центры
       * излучающих ячеек. Читает xcmp_report после ind_core_build. */
      if (xcmp && st.eirr != NULL && eind_dir != NULL) {
        g_xcmp_n = mesh.ncell;
        g_xcmp_E = calloc((size_t)mesh.ncell, sizeof *g_xcmp_E);
        g_xcmp_Edir = calloc((size_t)mesh.ncell, sizeof *g_xcmp_Edir);
        g_xcmp_aw = calloc((size_t)mesh.ncell, sizeof *g_xcmp_aw);
        g_xcmp_lo = malloc((size_t)mesh.ncell * sizeof *g_xcmp_lo);
        g_xcmp_sz = malloc((size_t)mesh.ncell * sizeof *g_xcmp_sz);
        unsigned char *em9 = calloc((size_t)mesh.ncell, 1);
        if (g_xcmp_E == NULL || g_xcmp_Edir == NULL || g_xcmp_aw == NULL || g_xcmp_lo == NULL ||
            g_xcmp_sz == NULL || em9 == NULL)
          exit(1);
        for (int32_t ci = 0; ci < mesh.ncell; ci++) {
          for (int a = 0; a < 3; a++)
            g_xcmp_lo[ci][a] = mesh.clo[ci][a];
          g_xcmp_sz[ci] = mesh.csize[ci];
        }
        for (int32_t e = 0; e < cut.nse; e++) {
          int32_t ci = cut.se[e].cell;
          double a = cut.se[e].area;
          if (ci < 0 || ci >= mesh.ncell || !(a > 0.0)) continue;
          g_xcmp_E[ci] += (st.eirr[e] - eind_dir[e]) * a;
          g_xcmp_Edir[ci] += eind_dir[e] * a;
          g_xcmp_aw[ci] += a;
          if (eemit[e] > 0.0) em9[ci] = 1;
        }
        int32_t nem = 0;
        for (int32_t ci = 0; ci < mesh.ncell; ci++)
          if (em9[ci]) nem++;
        g_xcmp_ec = malloc((size_t)(nem > 0 ? nem : 1) * sizeof *g_xcmp_ec);
        if (g_xcmp_ec == NULL) exit(1);
        g_xcmp_nemit = 0;
        for (int32_t ci = 0; ci < mesh.ncell; ci++) {
          if (g_xcmp_aw[ci] > 0.0) {
            g_xcmp_E[ci] /= g_xcmp_aw[ci];
            g_xcmp_Edir[ci] /= g_xcmp_aw[ci];
          }
          if (em9[ci]) {
            for (int a = 0; a < 3; a++)
              g_xcmp_ec[g_xcmp_nemit][a] =
                  fr.org[a] + ((double)mesh.clo[ci][a] + 0.5 * (double)mesh.csize[ci]) * fr.h;
            g_xcmp_nemit++;
          }
        }
        free(em9);
        printf("   §744 АГРЕГАТЫ СВИПА: ячеек %d, излучающих %d\n", g_xcmp_n, g_xcmp_nemit);
      }
      free(eind_dir);
      /* А1338/§806: эталон «целой ячейки» — ФАКТИЧЕСКИЙ размер ячейки
       * максимума (csize в единицах h), а не захардкоженное 2^(lev-6): в
       * каноне lev=7 у ячейки максимума csize = 1, старая формула печатала
       * (2h)³ — завышение ×8; при lev < 6 сдвиг 1<<(lev-6) — неопределённое
       * поведение (мусор −1.2e+25 на cavity lev=5). maxcell < 0 (все φ ≤ 0)
       * — эталона нет, 0.0. */
      double whole_at_max = maxcell >= 0 ? pow((double)mesh.csize[maxcell] * fr.h, 3.0) : 0.0;
      printf("      ГДЕ МАКСИМУМ: в полости %.4e, вне её %.4e; флюидный объём ячейки с "
             "максимумом %.3e м³ (у целой ячейки %.3e)\n",
             phimax_in, phimax_out, volmin_at_max, whole_at_max);
      /* ГДЕ ИМЕННО эта ячейка — координата, размер, занятость, класс заливки и
       * число поверхностных элементов. Без этих пяти чисел «вне полости»
       * остаётся догадкой: у ПОГРАНИЧНОЙ ячейки половина лежит по одну сторону
       * стены, половина по другую, и классифицировать её по нижнему углу
       * нельзя (найдено этим же прогоном). */
      if (maxcell >= 0) {
        size_t km =
            hz_occ_index(fr.n, mesh.clo[maxcell][0], mesh.clo[maxcell][1], mesh.clo[maxcell][2]);
        printf("      ЯЧЕЙКА МАКСИМУМА: мир (%.3f, %.3f, %.3f), размер %d, занята %d, класс "
               "заливки %d, поверхностных элементов %d\n",
               fr.org[0] + (double)mesh.clo[maxcell][0] * fr.h,
               fr.org[1] + (double)mesh.clo[maxcell][1] * fr.h,
               fr.org[2] + (double)mesh.clo[maxcell][2] * fr.h, mesh.csize[maxcell],
               hz_occ_get(P.b[lev], km) ? 1 : 0, solidmask != NULL ? (int)solidmask[km] : -1,
               cut.sestart[maxcell + 1] - cut.sestart[maxcell]);
      }
      /* Ш15: исходящий радианс поверхностных элементов — на СЕТКУ, чтобы его
       * могла прочесть картинка. Берётся постоянный член DG1; у ячейки с
       * несколькими элементами — наибольший (они суть куски одной поверхности
       * в одной ячейке, и брать среднее значило бы гасить край). */
      if (st.sout != NULL) {
        size_t ng = (size_t)fr.n * (size_t)fr.n * (size_t)fr.n;
        xrad = calloc(3 * ng, sizeof *xrad);
        if (xrad == NULL) exit(1);
        for (int32_t k = 0; k < cut.nse; k++) {
          int32_t ci = cut.se[k].cell;
          size_t g = hz_occ_index(fr.n, mesh.clo[ci][0], mesh.clo[ci][1], mesh.clo[ci][2]);
          float v = (float)st.sout[4 * (size_t)k];
          if (v > xrad[3 * g]) {
            xrad[3 * g] = v;
            xrad[3 * g + 1] = v * 0.92f;
            xrad[3 * g + 2] = v * 0.78f;
          }
        }
        /* ПИРАМИДА (Ш16, §486). Ячейка среза уровнем выше самого мелкого
         * адресуется своим НИЖНИМ УГЛОМ, а элементы живут только на мелком —
         * отсюда 55.4 % чёрных ячеек, читавшихся как тень. Свёртка вверх:
         * СРЕДНЕЕ ПО НЕПУСТЫМ детям. Пустой ребёнок значит «поверхности нет», а
         * не «темно»; включить его нулём — то же ложное затухание, что ловилось
         * в §424. */
        xradlv = calloc((size_t)lev + 1, sizeof *xradlv);
        if (xradlv == NULL) exit(1);
        xradlv[lev] = xrad;
        for (int l = lev - 1; l >= 0; l--) {
          int32_t nl = (int32_t)1 << l;
          xradlv[l] = calloc(3 * (size_t)nl * (size_t)nl * (size_t)nl, sizeof **xradlv);
          if (xradlv[l] == NULL) exit(1);
          const float *up = xradlv[l + 1];
          int32_t nu = nl * 2;
          for (int32_t z = 0; z < nl; z++)
            for (int32_t y = 0; y < nl; y++)
              for (int32_t x = 0; x < nl; x++) {
                double acc[3] = {0, 0, 0};
                int cnt = 0;
                for (int d = 0; d < 8; d++) {
                  size_t gu = hz_occ_index(nu, 2 * x + (d & 1), 2 * y + ((d >> 1) & 1),
                                           2 * z + ((d >> 2) & 1));
                  if (!(up[3 * gu] > 0.0f)) continue;
                  for (size_t c = 0; c < 3; c++)
                    acc[c] += (double)up[3 * gu + c];
                  cnt++;
                }
                if (cnt == 0) continue;
                size_t gl = hz_occ_index(nl, x, y, z);
                for (size_t c = 0; c < 3; c++)
                  xradlv[l][3 * gl + c] = (float)(acc[c] / (double)cnt);
              }
        }
      }
      /* §812: ДЕРЕВО ЛАМП ДЛЯ БЛИКОВ (машина §796: efcsrc/efc_build). Строится
       * ПО развёртке по фактическим eemit и живёт ДО КОНЦА прогона — кадровый
       * путь читает; сведение счётчиков — прибором после кадра. Порог вклада —
       * как у fc-сбора (§621/§796): доля средней радиосити с корня. */
      if (g_xgloss) {
        if (crc != 0 || cut.nse <= 0) {
          fprintf(stderr, "xgloss требует развёртку с элементами (xsweep)\n");
          exit(1);
        }
        int32_t nsrc812 = 0;
        for (int32_t k = 0; k < cut.nse; k++)
          if (eemit[k] > 0.0 && cut.se[k].area > 0.0) nsrc812++;
        if (nsrc812 > 0) {
          efcsrc *sr812 = malloc((size_t)nsrc812 * sizeof *sr812);
          if (sr812 == NULL) exit(1);
          int32_t js812 = 0;
          for (int32_t k = 0; k < cut.nse; k++) {
            if (!(eemit[k] > 0.0 && cut.se[k].area > 0.0)) continue;
            efcsrc *s8 = &sr812[js812];
            efc_center(&cut, &mesh, &fr, k, s8->c);
            for (int a = 0; a < 3; a++)
              s8->n[a] = cut.se[k].n[a];
            s8->area = cut.se[k].area;
            s8->flux = 3.14159265358979323846 * eemit[k] * cut.se[k].area;
            s8->idx = js812;
            s8->key = 0;
            for (int sh = lev - 1; sh >= 0; sh--) {
              uint64_t bit = (uint64_t)(((mesh.clo[cut.se[k].cell][0] >> sh) & 1) |
                                        (((mesh.clo[cut.se[k].cell][1] >> sh) & 1) << 1) |
                                        (((mesh.clo[cut.se[k].cell][2] >> sh) & 1) << 2));
              s8->key = (s8->key << 3) | bit;
            }
            js812++;
          }
          qsort(sr812, (size_t)nsrc812, sizeof *sr812, efc_cmp);
          g_glet = malloc(sizeof *g_glet);
          if (g_glet == NULL) exit(1);
          memset(g_glet, 0, sizeof *g_glet);
          if (efc_build(g_glet, sr812, 0, nsrc812, 0, lev) < 0) exit(1);
          double bsum812 = 0.0, asum812 = 0.0;
          for (int k = 0; k < g_glet->e[0].nb; k++) {
            const ebin *rb = &g_glet->b[g_glet->e[0].b0 + k];
            asum812 += (double)rb->area;
            bsum812 += (double)rb->flux[0];
          }
          g_gltau = (hcontrib_set ? g_hcontrib : HZ_FC_CONTRIB) *
                    (asum812 > 0.0 ? bsum812 / asum812 : 0.0);
          /* полуширина лепестка на полувысоте β = sqrt(2·ln2/ns) и косинус
           * отсечения (cps)^ns < 1e-3 — НОСИТЕЛИ узкие, берутся по медианному
           * ns сцены (приближение счётчика, не порог: точная ширина у каждого
           * материала своя, консервативность даёт оценка cmax) */
          double nsmed812 = 0.0;
          {
            double *nsv = malloc((size_t)(m.nmtl > 0 ? m.nmtl : 1) * sizeof *nsv);
            if (nsv == NULL) exit(1);
            int32_t nnv = 0;
            for (int32_t mi9 = 0; mi9 < m.nmtl; mi9++) {
              const double *ks3 = m.mtl[mi9].ks3;
              double ks9 = (ks3[0] + ks3[1] + ks3[2]) / 3.0;
              if (ks9 > 0.0 && m.mtl[mi9].ns > 0.0) nsv[nnv++] = m.mtl[mi9].ns;
            }
            if (nnv > 0) {
              qsort(nsv, (size_t)nnv, sizeof *nsv, cmp_dev699);
              nsmed812 = nsv[nnv / 2];
            }
            free(nsv);
          }
          g_glbeta = nsmed812 > 0.0 ? sqrt(2.0 * 0.69314718055994530942 / nsmed812) : 1.0;
          g_glcext = nsmed812 > 0.0 ? pow(1e-3, 1.0 / nsmed812) : -1.0;
          g_glP = &P;
          g_glfr = &fr;
          printf("   §812 ДЕРЕВО ЛАМП: излучателей %d, корзин %d; порог вклада %.3e; "
                 "медианный Ns сцены %.1f (β %.3f рад, cext %.3f)\n",
                 (int)nsrc812, (int)g_glet->nb, g_gltau, nsmed812, g_glbeta, g_glcext);
          free(sr812);
        } else
          printf("   §812 ДЕРЕВО ЛАМП: светящихся элементов нет — канал пуст\n");
        if (gloss812_selftest() != 0) exit(1);
      }
      free(st.bout);
      free(st.eirr);
      free(st.sout);
      tr3_dirs_free(&dirs);
      free(frho);
      free(efc796);
      free(einj796);
      free(eemit);
      free(femit);
      free(sig_t);
      free(sig_s);
      free(phi);
    }
    if (crc == 0) tr3_cut_free(&cut);
    free(solid);
    tr3_mesh_free(&mesh);
    hz_oct_free(&ot);
    hz_facettab_free(&ftab);
    hz_cutmap_free(&cmap);
  }

  /* ---- 4. приёмка: потеря поля (популяция — ВХОД) ---- */
  if (!g_render) surf_err(&T, &fr, &m, ct_list, &CT, "");
  if (!g_render) {
    covstat st;
    occ_cover(&T, &fr, P.b[lev], &ht, &st);
    printf("   ПОТЕРЯНО ПОЛЕМ: занятых ячеек %lld; ОКНО 27 — %lld (%.2f %%); ОКНО 1 — %lld "
           "(%.2f %%), из них БЕЗ ЕДИНОГО ПЕРЕСЕЧЁННОГО РЕБРА %lld (%.2f %% от занятых)\n",
           (long long)st.nocc, (long long)st.lost27,
           100.0 * (double)st.lost27 / (double)(st.nocc ? st.nocc : 1), (long long)st.lost1,
           100.0 * (double)st.lost1 / (double)(st.nocc ? st.nocc : 1), (long long)st.lost_noedge,
           100.0 * (double)st.lost_noedge / (double)(st.nocc ? st.nocc : 1));
  }

  /* ---- 5. картинка ---- */
  {
    tr3_camera cam;
    double eye[3] = {g_eye[0], g_eye[1], g_eye[2]}, at[3] = {g_at[0], g_at[1], g_at[2]},
           up[3] = HZ_CFG_UP;
    if (tr3_camera_look(&cam, eye, at, up, HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0, 512,
                        512) == 0) {
      size_t np = (size_t)512u * 512u;
      double *zb = malloc(np * sizeof *zb), *vb = calloc(np, sizeof *vb);
      if (zb == NULL || vb == NULL) exit(1);
      for (size_t i = 0; i < np; i++)
        zb[i] = 1e300;
      rastctx R = {&cam, &fr, zb, vb, 0, 0};
      t0 = now_s();
      int wrc = hz_dc_walk(&T, NULL, NULL, rast_poly, &R);
      double t_img = now_s() - t0;
      int64_t nfill = 0;
      for (size_t i = 0; i < np; i++)
        if (zb[i] < 1e299) nfill++;
      char path[256];
      snprintf(path, sizeof path, "img/pfield_hall_L%d%s.ppm", lev, nonrm ? "_nonrm" : "");
      int prc = hz_ppm_write(path, vb, 512, 512);
      printf("   КАРТИНКА за %.2f с: обход код %d, треугольников %lld (за камерой %lld), "
             "заполнено %lld из %zu пикселей (%.1f %%) -> %s (код %d)\n",
             t_img, wrc, (long long)R.ntri, (long long)R.nclip, (long long)nfill, np,
             100.0 * (double)nfill / (double)np, path, prc);
      free(zb);
      free(vb);
    }
  }

  /* ---- 6. эталон регресса для Ш2 (А666) ---- */
  if (polydump != NULL) {
    FILE *f = fopen(polydump, "wb");
    if (f == NULL)
      printf("   ЭТАЛОН МНОГОУГОЛЬНИКОВ: не открыть %s\n", polydump);
    else {
      polydumpctx D = {f, 0};
      int wrc = hz_dc_walk(&T, NULL, NULL, poly_dump, &D);
      printf("   ЭТАЛОН МНОГОУГОЛЬНИКОВ -> %s: код %d, многоугольников %lld\n", polydump, wrc,
             (long long)D.n);
      fclose(f);
    }
  }

  /* ---- 4б. СРЕЗ (Ш2): построение, замер обхода, сверка вершин ---- */
  if (!g_render) {
    hz_dcslice S;
    if (hz_slice_init(&S, lev) != HZ_DC_OK) exit(1);
    if (vq1) S.vbits = 1; /* НЕГАТИВНЫЙ КОНТРОЛЬ §397: вершина в один бит на ось */
    t0 = now_s();
    int src = hz_slice_build(&S, &T, &ht, NULL, NULL);
    double t_slice = now_s() - t0;
    printf("   СРЕЗ за %.2f с: код %d, ячеек %d, %zu Б/ячейка, %.1f МБ%s\n", t_slice, src, S.n,
           sizeof(hz_dccell), (double)S.n * (double)sizeof(hz_dccell) / 1048576.0,
           vq1 ? "  [ВЕРШИНА В 1 БИТ]" : "");

    /* ЗАМЕР ГОРЯЧЕГО ПУТИ. Обход читает ровно то, что будет читать фронт:
     * вершину и нормаль каждой ячейки. Результат СУММИРУЕТСЯ и печатается —
     * иначе компилятор вправе выбросить цикл целиком, и замер померяет пустоту.
     * Проходов несколько: первый греет кэш, и мешать холодный с установившимся
     * нельзя (правило раздельного доклада). */
    /* ЦЕНА РАЗЛАГАЕТСЯ НА ТРИ, иначе «нс на ячейку» не скажет, ЧТО именно
     * дорого. Вариант 0 — чистый ПОТОК (сложение самих байт, ничего не
     * декодируется): это пол, задаваемый памятью. Вариант 1 добавляет
     * распаковку ВЕРШИНЫ, а с ней разбор мортонова кода (цикл по уровням).
     * Вариант 2 добавляет распаковку НОРМАЛИ (октаэдр, корень). Разность
     * соседних вариантов и есть цена каждой части. */
    double acc[3] = {0.0, 0.0, 0.0}, t_first[3] = {0, 0, 0}, t_best[3] = {1e300, 1e300, 1e300};
    for (int mode = 0; mode < 3; mode++)
      for (int pass = 0; pass < 5; pass++) {
        double ta = now_s();
        double a = 0.0;
        for (int32_t i = 0; i < S.n; i++) {
          if (mode == 0) {
            a += (double)S.c[i].vx[0] + (double)S.c[i].vx[1] + (double)S.c[i].vx[2] +
                 (double)S.c[i].noct + (double)S.c[i].lo[0] + (double)S.c[i].lo[1] +
                 (double)S.c[i].lo[2] + (double)S.c[i].lvl;
            continue;
          }
          double v[3];
          hz_slice_vertex(&S, i, v);
          if (mode == 1) {
            a += v[0] + v[1] + v[2] + (double)S.c[i].noct;
            continue;
          }
          double nn[3];
          hz_slice_normal(&S, i, nn);
          a += v[0] * nn[0] + v[1] * nn[1] + v[2] * nn[2];
        }
        double dt = now_s() - ta;
        if (pass == 0) t_first[mode] = dt;
        if (dt < t_best[mode]) t_best[mode] = dt;
        acc[mode] = a;
      }
    double per = 1e9 / (double)(S.n ? S.n : 1);
    printf("   ОБХОД СРЕЗА (установившийся, нс/ячейка): ПОТОК %.2f; +ВЕРШИНА %.2f; "
           "+НОРМАЛЬ %.2f. ХОЛОДНЫЙ полный %.2f мс (%.2f нс/ячейка)\n",
           t_best[0] * per, t_best[1] * per, t_best[2] * per, t_first[2] * 1e3, t_first[2] * per);
    printf("   контрольные суммы %.6e %.6e %.6e\n", acc[0], acc[1], acc[2]);

    /* ОШИБКА ПЕРЕКОДИРОВКИ ВЕРШИНЫ `16 -> 8` БИТ. Это НЕ «сверка среза с
     * деревом», как здесь было написано до 08-19, и разница не в словах:
     * обе стороны читают ОДНО И ТО ЖЕ квантованное `vq` (`collect_verts` пишет
     * `lo + vq·size/65535`, `hz_slice_vertex` раскодирует то, что `sb_push` взял
     * из того же `hz_dc_vertex`), поэтому величина слепа к кванту УЗЛА
     * СТРУКТУРНО и мерит ровно одно — потерю на перекодировке в срез.
     * ЗАМЕРЕНО, а не выведено (§640, ключ `nodevq1`): при вершине узла в ОДИН
     * бит на ось эта величина печатает `0.000e+00`, идеальное согласие, при
     * ошибке поверхности `p50 0.000000 -> 0.005882` м. Прежнее название и стало
     * причиной того, что П4 §635 был предсказан на неверной посылке.
     * Популяция — ячейки СРЕЗА (их столько же, сколько вершин у листьев дерева):
     * величина ловит потерю точности, а не потерю ячеек. */
    {
      int64_t nvt = 0;
      collect_verts(&T, 0, zero, fr.n, NULL, 0, &nvt, 1);
      double dmax = 0.0;
      vrec *vt = malloc((size_t)(nvt > 0 ? nvt : 1) * sizeof *vt);
      if (vt == NULL) exit(1);
      int64_t kk = 0;
      collect_verts(&T, 0, zero, fr.n, vt, nvt, &kk, 1);
      /* Обход дерева в collect_verts идёт тем же порядком детей 0..7, что и
       * сборка среза, поэтому сопоставление ПОРЯДКОВОЕ и второго индекса не
       * заводит. Если порядки разойдутся, расхождение будет огромным, а не
       * тонким, — то есть проверка не слепа к собственной ошибке. */
      int64_t nn = (int64_t)S.n < kk ? (int64_t)S.n : kk;
      for (int64_t i = 0; i < nn; i++) {
        double v[3];
        hz_slice_vertex(&S, (int32_t)i, v);
        for (int a = 0; a < 3; a++) {
          double d = fabs(v[a] - vt[i].vx[a]) * fr.h;
          if (d > dmax) dmax = d;
        }
      }
      free(vt);
      double quant = fr.h / (double)((1u << S.vbits) - 1u);
      /* ОШИБКА НОРМАЛИ меряется НА КОДЕКЕ, а не на ячейке: берутся нормали
       * таблицы рёбер, кодируются и раскодируются, считается угол. Так и
       * говорится в докладе — это ошибка ОКТАЭДРИЧЕСКОГО КОДА, и она НЕ
       * покрывает ошибки суммирования нормалей по ячейке. */
      double a90 = 0.0, amax = 0.0;
      {
        double *ang = malloc((size_t)(ht.n > 0 ? ht.n : 1) * sizeof *ang);
        if (ang == NULL) exit(1);
        for (int32_t i = 0; i < ht.n; i++) {
          double d[3];
          hz_oct_decode(hz_oct_encode(ht.e[i].nrm), d);
          double c = d[0] * ht.e[i].nrm[0] + d[1] * ht.e[i].nrm[1] + d[2] * ht.e[i].nrm[2];
          if (c > 1.0) c = 1.0;
          if (c < -1.0) c = -1.0;
          ang[i] = acos(c) * 180.0 / 3.14159265358979323846;
        }
        qsort(ang, (size_t)ht.n, sizeof *ang, cmp_d);
        if (ht.n > 0) {
          a90 = ang[(ht.n * 9) / 10];
          amax = ang[ht.n - 1];
        }
        free(ang);
      }
      /* РАЗМЕР ЯЧЕЙКИ ПЕЧАТАЕТСЯ РЯДОМ, И ЭТО НЕ УКРАШЕНИЕ (А1046): отношение
       * `dmax/quant` САМОСОКРАЩАЮЩЕЕСЯ — при `vq1` числитель и знаменатель
       * выросли ОБА в `256` раз, а отношение осталось `0.500` неподвижным.
       * Судить о величине можно только против `h`, которое от кванта среза не
       * зависит. Отношение оставлено, потому что показывает работу округления к
       * ближайшему, и подписано тем, чем является.
       * Метка режима стоит В КОНЦЕ строки сознательно. */
      printf("   ОШИБКА ПЕРЕКОДИРОВКИ ВЕРШИНЫ 16 -> 8 БИТ (в срез; до §642 звалась «СВЕРКА "
             "СРЕЗА С ДЕРЕВОМ»): ячеек среза %d, вершин дерева %lld; МАКС ОШИБКА %.3e м при "
             "ячейке %.3e м, квант среза %.3e м (отношение к кванту %.3f — самосокращающееся); "
             "УГОЛ КОДА НОРМАЛИ p90 %.3f, макс %.3f град%s%s\n",
             S.n, (long long)nvt, dmax, fr.h, quant, dmax / quant, a90, amax,
             vq1 ? "  [ВЕРШИНА СРЕЗА В 1 БИТ]" : "", nodevq1 ? "  [ВЕРШИНА УЗЛА В 1 БИТ]" : "");
    }
    hz_slice_free(&S);
  }

  /* ---- 4в. АДАПТИВНЫЙ СРЕЗ (Ш3): критерий, ошибка ДО ВЫДАННОЙ ПОВЕРХНОСТИ,
   *      ориентация на перепаде уровней, доля среза за кадр ---- */
  if (!g_render) {
    lodctx L;
    memset(&L, 0, sizeof L);
    {
      double eyec[3] = {g_eye[0], g_eye[1], g_eye[2]};
      for (int a = 0; a < 3; a++)
        L.eye[a] = (eyec[a] - fr.org[a]) / fr.h;
      /* А737 ПРОВЕРЯЕТСЯ ПРЯМО: сколько ВНУТРЕННИХ узлов имеют невязку РОВНО НОЛЬ
       * по уровням. Узел уровня 2 накрывает 1/64 сцены и плоским быть не может;
       * ноль у него означает не «идеально подогнано», а сокращение близких больших
       * чисел, о котором предупреждает сам `qef.h`. */
      {
        int64_t nz[HZ_DC_MAX_LOG2SIZE + 2], nt2[HZ_DC_MAX_LOG2SIZE + 2];
        for (int i = 0; i <= HZ_DC_MAX_LOG2SIZE + 1; i++)
          nz[i] = nt2[i] = 0;
        for (int32_t i = 0; i < T.n; i++) {
          if (T.nd[i].child0 < 0) continue;
          if (!hz_dc_hasvert(&T, i)) continue;
          /* уровень узла восстанавливается по числу образцов? нет — по глубине,
           * а глубины у узла нет; поэтому считается обходом ниже */
          (void)0;
        }
        /* обход с глубиной */
        int32_t st[64];
        int32_t sl[64];
        int sp = 0;
        st[sp] = 0;
        sl[sp] = 0;
        sp = 1;
        while (sp > 0) {
          sp--;
          int32_t ni = st[sp];
          int lv = sl[sp];
          if (T.nd[ni].child0 < 0) continue;
          if (hz_dc_hasvert(&T, ni)) {
            nt2[lv]++;
            if (!(hz_dc_rms(&T, ni) > 0.0)) nz[lv]++;
          }
          for (int k = 0; k < 8 && sp < 60; k++) {
            st[sp] = T.nd[ni].child0 + k;
            sl[sp] = lv + 1;
            sp++;
          }
        }
        printf("   А737 НЕВЯЗКА РОВНО НОЛЬ У ВНУТРЕННИХ УЗЛОВ ПО УРОВНЯМ:");
        for (int i = 0; i <= lev; i++)
          if (nt2[i] > 0) printf(" L%d %lld/%lld", i, (long long)nz[i], (long long)nt2[i]);
        printf("\n");
      }
    }
    /* Радиан на пиксель — из поля зрения и разрешения §2, а не константой. */
    L.pxrad = (HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0) / 512.0;

    /* ЭТАЛОН ПОЛНОЙ ГЛУБИНЫ ДЛЯ А729. Без него доля обращённых не читается:
     * часть несовпадений идёт от самой модели (суп треугольников с
     * несогласованным обходом), и отделить её можно только сравнением с
     * режимом, где перепада уровней нет вовсе. */
    {
      emesh E0;
      E0.ntri = 0;
      E0.cap = 1024;
      E0.v = malloc((size_t)E0.cap * 9 * sizeof *E0.v);
      E0.cls = malloc((size_t)E0.cap * sizeof *E0.cls);
      E0.gn = (int32_t)1 << (lev < HZ_EGRID ? lev : HZ_EGRID);
      E0.lev = lev;
      E0.tri = NULL;
      E0.head = malloc((size_t)E0.gn * (size_t)E0.gn * (size_t)E0.gn * sizeof *E0.head);
      E0.nxt = NULL;
      E0.ncap = 0;
      E0.nn = 0;
      if (E0.v == NULL || E0.cls == NULL || E0.head == NULL) exit(1);
      for (int64_t i = 0; i < (int64_t)E0.gn * E0.gn * E0.gn; i++)
        E0.head[i] = -1;
      /* HZ_DC_EMULTI (код 4) НЕ ФАТАЛЕН: обход так помечает ячейки, где вершины
       * нет, и продолжает — «дыра честнее пустого экрана» (dc.c, §377). Считать
       * его отказом значило бы падать на любой сцене, где такие ячейки есть; на
       * зале их 0, на замкнутой коробке — нет, и стенд падал молча. Код
       * печатается, а не глотается. */
      {
        int wrc0 = hz_dc_walk(&T, NULL, NULL, emesh_emit, &E0);
        if (wrc0 != HZ_DC_OK && wrc0 != HZ_DC_EMULTI) exit(1);
        if (wrc0 == HZ_DC_EMULTI) printf("   (обход эталона: код 4, часть ячеек без вершины)\n");
      }
      int64_t nc0 = 0, nfc0[3], ncc0[3];
      flipdiag D0;
      flipdiag_init(&D0, E0.ntri);
      int64_t nf0 = emesh_flips(&E0, &CT, &fr, &m, &nc0, nfc0, ncc0, &D0);
      printf("   ЭТАЛОН ПОЛНОЙ ГЛУБИНЫ: треугольников %lld, ОБРАЩЁННЫХ %lld из %lld (%.3f %%) — "
             "это доля МОДЕЛИ, перепада уровней здесь нет\n",
             (long long)E0.ntri, (long long)nf0, (long long)nc0,
             100.0 * (double)nf0 / (double)(nc0 ? nc0 : 1));
      /* А815: класс F при полной глубине — это ВСЕ полигоны, и печатается он
       * затем, чтобы было с чем сравнивать класс F на срезе. Разойдись они —
       * класс F не представителен, и критерий механизма читать нельзя. */
      printf("      из них класс МЕЛКИЙ %lld из %lld (%.3f %%), класс ПРАВИЛО %lld из %lld, класс "
             "СМЕШАННЫЙ %lld из %lld\n",
             (long long)nfc0[0], (long long)ncc0[0],
             100.0 * (double)nfc0[0] / (double)(ncc0[0] ? ncc0[0] : 1), (long long)nfc0[1],
             (long long)ncc0[1], (long long)nfc0[2], (long long)ncc0[2]);
      flipdiag_report(&D0, "полная глубина");
      printf("      Ш10 СЛОЖЕННЫХ КВАДОВ: %lld из %lld (%.3f %%)%s\n", (long long)g_fold,
             (long long)g_quad, 100.0 * (double)g_fold / (double)(g_quad ? g_quad : 1),
             hz_longdiag ? "   [ДЛИННАЯ ДИАГОНАЛЬ — НК]" : "");
      free(E0.v);
      free(E0.cls);
      free(E0.head);
      free(E0.nxt);
      free(E0.tri);
    }
    static const double thrs[6] = {0.0, 0.25, 0.5, 1.0, 2.0, 1e30};
    hz_dcslice A;
    if (hz_slice_init(&A, lev) != HZ_DC_OK) exit(1);
    for (int k = 0; k < 6; k++) {
      L.thr = thrs[k];
      for (int i = 0; i <= HZ_DC_MAX_LOG2SIZE + 1; i++)
        L.hist[i] = 0;
      t0 = now_s();
      if (hz_slice_build(&A, &T, &ht, lod_stop, &L) != HZ_DC_OK) exit(1);
      double ts = now_s() - t0;
      for (int32_t i = 0; i < A.n; i++)
        L.hist[A.c[i].lvl]++;
      /* А813: ЛОЖНЫЙ НОЛЬ НОРМАЛИ У КРУПНОЙ ЯЧЕЙКИ. `hz_slice_build` набирает
       * сумму нормалей только при size == 1, у крупной ячейки код нормали
       * остаётся нулём — а раскодируется он не в «нормали нет», а в вектор
       * (0, 0, −1). Здесь это СЧИТАЕТСЯ, порознь по мелким и крупным. Чинится
       * не в Ш7 (объём шага): число нужно, чтобы находка не осталась чтением. */
      int64_t nz0f = 0, nz0c = 0, ncoarse = 0;
      for (int32_t i = 0; i < A.n; i++) {
        int fine = A.c[i].lvl == (uint8_t)lev;
        if (!fine) ncoarse++;
        if (A.c[i].noct != 0) continue;
        if (fine)
          nz0f++;
        else
          nz0c++;
      }
      /* УГОЛ НОРМАЛИ ЯЧЕЙКИ ПРОТИВ ИСТИННОЙ (Ш9, §459). Эталон — нормаль
       * БЛИЖАЙШЕГО исходного треугольника к вершине ячейки; «первый в ячейке»
       * запрещён после А824. ОГОВОРКА А835, без неё число не читается: у
       * крупной ячейки сравнивается СРЕДНЯЯ по куску нормаль с МЕСТНОЙ, поэтому
       * метрика меряет «наша ошибка + кривизна куска» и систематически завышает.
       * Мелкие ячейки печатаются рядом как внутренний эталон: их путь не тронут,
       * и их угол обязан не измениться ни у одной.
       * НЕГАТИВНЫЙ КОНТРОЛЬ (флаг `nrmflip`): нормаль ячейки подаётся
       * ОБРАЩЁННОЙ, и угол обязан стать 180° − x. Не станет — метрика меряет
       * длину, а не направление, и всё, что она подтвердила, не подтверждено. */
      {
        double *ang[2];
        int64_t na2[2] = {0, 0}, nomiss = 0;
        for (int cl2 = 0; cl2 < 2; cl2++) {
          ang[cl2] = malloc((size_t)(A.n > 0 ? A.n : 1) * sizeof *ang[cl2]);
          if (ang[cl2] == NULL) exit(1);
        }
        for (int32_t i = 0; i < A.n; i++) {
          double v[3], nn[3];
          hz_slice_vertex(&A, i, v);
          hz_slice_normal(&A, i, nn);
          if (nrmflip)
            for (int c = 0; c < 3; c++)
              nn[c] = -nn[c];
          double w[3];
          for (int c = 0; c < 3; c++)
            w[c] = fr.org[c] + v[c] * fr.h;
          int32_t cell[3];
          int ok = 1;
          for (int c = 0; c < 3; c++) {
            double f = floor((w[c] - fr.org[c]) / fr.h);
            if (!(f >= 0.0) || !(f < (double)fr.n)) ok = 0;
            cell[c] = ok ? (int32_t)f : 0;
          }
          const int32_t *ls = NULL;
          int32_t nls = ok ? ct_list(&CT, cell, &ls) : 0;
          if (nls == 0) {
            nomiss++;
            continue;
          }
          double bd = 1e300, fn[3] = {0.0, 0.0, 0.0};
          for (int32_t q = 0; q < nls; q++) {
            const double *Aq, *Bq, *Cq;
            tri_verts(&m, ls[q], &Aq, &Bq, &Cq);
            double dq = pt_tri_d2(w, Aq, Bq, Cq);
            if (dq >= bd) continue;
            bd = dq;
            double f1[3], f2[3];
            for (int c = 0; c < 3; c++) {
              f1[c] = Bq[c] - Aq[c];
              f2[c] = Cq[c] - Aq[c];
            }
            fn[0] = f1[1] * f2[2] - f1[2] * f2[1];
            fn[1] = f1[2] * f2[0] - f1[0] * f2[2];
            fn[2] = f1[0] * f2[1] - f1[1] * f2[0];
          }
          double fl = sqrt(fn[0] * fn[0] + fn[1] * fn[1] + fn[2] * fn[2]);
          if (!(fl > 0.0)) {
            nomiss++;
            continue;
          }
          double cs = (nn[0] * fn[0] + nn[1] * fn[1] + nn[2] * fn[2]) / fl;
          if (cs > 1.0) cs = 1.0;
          if (cs < -1.0) cs = -1.0;
          int fine = A.c[i].lvl == (uint8_t)lev;
          ang[fine ? 0 : 1][na2[fine ? 0 : 1]++] = acos(cs) * 180.0 / 3.14159265358979323846;
        }
        double q50[2] = {0, 0}, q90[2] = {0, 0}, qmx[2] = {0, 0};
        for (int cl2 = 0; cl2 < 2; cl2++) {
          if (na2[cl2] > 0) {
            qsort(ang[cl2], (size_t)na2[cl2], sizeof *ang[cl2], cmp_d);
            q50[cl2] = ang[cl2][na2[cl2] / 2];
            q90[cl2] = ang[cl2][(na2[cl2] * 9) / 10];
            qmx[cl2] = ang[cl2][na2[cl2] - 1];
          }
          free(ang[cl2]); /* БЕЗУСЛОВНО: пустой класс тоже был выделен (утечка
                           * 192 Б, пойманная санитайзером на зале — А833). */
        }
        printf("      Ш9 УГОЛ НОРМАЛИ ЯЧЕЙКИ (град, эталон — ближайший исходный треугольник%s): "
               "МЕЛКИЕ p50 %.1f p90 %.1f max %.1f (%lld шт); КРУПНЫЕ p50 %.1f p90 %.1f max %.1f "
               "(%lld шт); без эталона %lld\n",
               nrmflip ? ", НОРМАЛЬ ОБРАЩЕНА — НК" : "", q50[0], q90[0], qmx[0], (long long)na2[0],
               q50[1], q90[1], qmx[1], (long long)na2[1], (long long)nomiss);
      }
      int64_t half = A.n / 2, acc2 = 0;
      int lmed = 0, lmin = 99, lmax = -1;
      for (int i = 0; i <= HZ_DC_MAX_LOG2SIZE + 1; i++) {
        if (L.hist[i] == 0) continue;
        if (i < lmin) lmin = i;
        if (i > lmax) lmax = i;
        acc2 += L.hist[i];
        if (acc2 <= half) lmed = i;
      }

      /* ВЫДАННАЯ ПОВЕРХНОСТЬ СОБИРАЕТСЯ ЦЕЛИКОМ (А735). Мерить расстояние до
       * ВЕРШИНЫ накрывающей ячейки — подмена: поверхность есть сетка
       * многоугольников, а не набор вершин, и у крупной ячейки её единственная
       * вершина отстоит от края на размер ячейки. Здесь треугольники обхода
       * складываются в массив и привязываются к КАЖДОЙ из своих ячеек списком,
       * чтобы поиск ближайшего шёл по O(1) кандидатов, а не по всей сетке. */
      emesh EM;
      EM.ntri = 0;
      EM.cap = 1024;
      EM.v = malloc((size_t)EM.cap * 9 * sizeof *EM.v);
      EM.cls = malloc((size_t)EM.cap * sizeof *EM.cls);
      EM.gn = (int32_t)1 << (lev < HZ_EGRID ? lev : HZ_EGRID);
      EM.lev = lev;
      EM.tri = NULL;
      EM.head = malloc((size_t)EM.gn * (size_t)EM.gn * (size_t)EM.gn * sizeof *EM.head);
      EM.nxt = NULL;
      EM.ncap = 0;
      EM.nn = 0;
      if (EM.v == NULL || EM.cls == NULL || EM.head == NULL) exit(1);
      for (int64_t i = 0; i < (int64_t)EM.gn * EM.gn * EM.gn; i++)
        EM.head[i] = -1;
      {
        int wrcm = hz_dc_walk(&T, lod_stop, &L, emesh_emit, &EM);
        if (wrcm != HZ_DC_OK && wrcm != HZ_DC_EMULTI) exit(1);
      }

      /* ОШИБКА: от ИСТИННОЙ поверхности до ВЫДАННОЙ. ПОПУЛЯЦИЯ — ВХОД (занятые
       * ячейки эталона Ш0). Точка на истинной поверхности — центр тяжести куска
       * треугольника в занятой ячейке; расстояние — до ближайшей точки
       * треугольника выданной сетки, привязанного к накрывающему узлу. */
      double *er = malloc((size_t)(nocc > 0 ? nocc : 1) * sizeof *er);
      if (er == NULL) exit(1);
      int64_t ne = 0, nmiss = 0;
      for (int64_t z = 0; z < fr.n; z++)
        for (int64_t y = 0; y < fr.n; y++)
          for (int64_t x = 0; x < fr.n; x++) {
            size_t ci = hz_occ_index(fr.n, x, y, z);
            if (!hz_occ_get(P.b[lev], ci)) continue;
            int32_t cell[3] = {(int32_t)x, (int32_t)y, (int32_t)z};
            const int32_t *ls = NULL;
            if (ct_list(&CT, cell, &ls) == 0) continue;
            double cl[3], ch[3];
            cell_box(&fr, x, y, z, cl, ch);
            const double *Aa, *Bb, *Cc;
            tri_verts(&m, ls[0], &Aa, &Bb, &Cc);
            hz_pclip_poly Q;
            if (hz_pclip_tri(Aa, Bb, Cc, cl, ch, &Q) < 3) continue;
            double tp[3] = {0, 0, 0};
            for (int q = 0; q < Q.nv; q++)
              for (int c = 0; c < 3; c++)
                tp[c] += Q.v[q][c] / (double)Q.nv;
            double best = emesh_nearest(&EM, &fr, tp);
            if (best > 1e299) {
              nmiss++;
              continue;
            }
            er[ne++] = sqrt(best);
          }
      double p50 = 0, p90 = 0, p99 = 0;
      if (ne > 0) {
        qsort(er, (size_t)ne, sizeof *er, cmp_d);
        p50 = er[ne / 2];
        p90 = er[(ne * 9) / 10];
        p99 = er[(ne * 99) / 100];
      }
      free(er);
      {
        ambctx C;
        memset(&C, 0, sizeof C);
        C.ht = &ht;
        C.t = &T;
        C.stop = lod_stop;
        C.sctx = &L;
        int32_t zlo[3] = {0, 0, 0};
        amb_rec(&C, 0, zlo, fr.n);
        printf("      §457 НЕОДНОЗНАЧНОСТЬ ПРАВИЛА: крупных листьев %lld; их рёбер с "
               "пересечением %lld, из них с ДВУМЯ и более %lld, со ВСТРЕЧНЫМИ знаками %lld; "
               "сумма РОВНО ноль %lld; БИТ ИЗМЕНИЛСЯ БЫ на %lld рёбрах\n",
               (long long)C.nbig, (long long)C.nedge1, (long long)C.nedgem, (long long)C.ndis,
               (long long)C.nzero, (long long)C.nchange);
      }
      printf("      А813 НОРМАЛЬ СРЕЗА == 0 (раскодируется в (0,0,-1)): у мелких %lld из %lld, у "
             "КРУПНЫХ %lld из %lld\n",
             (long long)nz0f, (long long)(A.n - ncoarse), (long long)nz0c, (long long)ncoarse);
      printf("   СРЕЗ порог %.2f пикс: ячеек %d, доля %.4f, уровни %d..%d медиана %d, сборка "
             "%.3f с, треугольников %lld; ОШИБКА ДО ВЫДАННОЙ ПОВЕРХНОСТИ p50 %.5f p90 %.5f "
             "p99 %.5f м; без поверхности %lld из %lld\n",
             L.thr, A.n, (double)A.n / (double)(nvl ? nvl : 1), lmin, lmax, lmed, ts,
             (long long)EM.ntri, p50, p90, p99, (long long)nmiss, (long long)(ne + nmiss));
      {
        int64_t nc2 = 0, nfc[3], ncc[3];
        flipdiag D2;
        flipdiag_init(&D2, EM.ntri);
        int64_t nf2 = emesh_flips(&EM, &CT, &fr, &m, &nc2, nfc, ncc, &D2);
        printf("      ОБРАЩЁННЫХ (А729, дифференциально): %lld из %lld (%.3f %%)\n", (long long)nf2,
               (long long)nc2, 100.0 * (double)nf2 / (double)(nc2 ? nc2 : 1));
        /* ТРИ КЛАССА (§457). МЕЛКИЙ — ни правило, ни огрубление не участвуют:
         * эталон в этом же прогоне. ПРАВИЛО — бит из свёртки по детям, то есть
         * предмет Ш7. СМЕШАННЫЙ — бит из таблицы, но геометрия огрублена: он и
         * отделяет вклад ПРАВИЛА от вклада ОГРУБЛЕНИЯ. Без третьего класса
         * первые два неразличимы, и первая редакция замера на этом и сбилась. */
        printf("      МЕЛКИЙ %lld/%lld (%.3f %%); ПРАВИЛО (min>1) %lld/%lld (%.3f %%); СМЕШАННЫЙ "
               "%lld/%lld (%.3f %%); доли по числу %.1f/%.1f/%.1f %%\n",
               (long long)nfc[0], (long long)ncc[0],
               100.0 * (double)nfc[0] / (double)(ncc[0] ? ncc[0] : 1), (long long)nfc[1],
               (long long)ncc[1], 100.0 * (double)nfc[1] / (double)(ncc[1] ? ncc[1] : 1),
               (long long)nfc[2], (long long)ncc[2],
               100.0 * (double)nfc[2] / (double)(ncc[2] ? ncc[2] : 1),
               100.0 * (double)ncc[0] / (double)(nc2 ? nc2 : 1),
               100.0 * (double)ncc[1] / (double)(nc2 ? nc2 : 1),
               100.0 * (double)ncc[2] / (double)(nc2 ? nc2 : 1));
        flipdiag_report(&D2, "срез");
      }
      free(EM.v);
      free(EM.cls);
      free(EM.head);
      free(EM.nxt);
      free(EM.tri);
    }

    /* ДОЛЯ СРЕЗА, МЕНЯЮЩАЯСЯ ЗА КАДР (П3.4, А674). Скорость и частота названы
     * ДО прогона: ходьба 1.4 м/с при 30 кадрах в секунду — шаг 4.667 см.
     * ЭТО ВЕРХНЯЯ ОЦЕНКА ПОЛЬЗЫ инкрементальности, а не её цена (А732):
     * меряется, СКОЛЬКО ячеек изменилось между двумя ПОЛНЫМИ сборками. */
    {
      L.thr = 1.0;
      hz_dcslice B;
      if (hz_slice_init(&B, lev) != HZ_DC_OK) exit(1);
      if (hz_slice_build(&A, &T, &ht, lod_stop, &L) != HZ_DC_OK) exit(1);
      double step_m = 1.4 / 30.0;
      /* шаг вдоль взгляда: камера идёт туда, куда смотрит */
      double eyec[3] = {g_eye[0], g_eye[1], g_eye[2]}, atc[3] = {g_at[0], g_at[1], g_at[2]}, dir[3];
      double dl = 0.0;
      for (int a = 0; a < 3; a++) {
        dir[a] = atc[a] - eyec[a];
        dl += dir[a] * dir[a];
      }
      dl = sqrt(dl);
      for (int a = 0; a < 3; a++)
        L.eye[a] += (dir[a] / dl) * step_m / fr.h;
      if (hz_slice_build(&B, &T, &ht, lod_stop, &L) != HZ_DC_OK) exit(1);
      /* Сличаются МНОЖЕСТВА по ключу (координата, уровень). Оба среза уже в
       * мортоновом порядке, но порядок при разных срезах разный, поэтому
       * ключи сортируются — сверять построчно нельзя (класс А720). */
      int64_t na = A.n, nb = B.n;
      uint64_t *ka = malloc((size_t)(na > 0 ? na : 1) * sizeof *ka);
      uint64_t *kb = malloc((size_t)(nb > 0 ? nb : 1) * sizeof *kb);
      if (ka == NULL || kb == NULL) exit(1);
      for (int64_t i = 0; i < na; i++)
        ka[i] = cellkey(&A.c[i]);
      for (int64_t i = 0; i < nb; i++)
        kb[i] = cellkey(&B.c[i]);
      qsort(ka, (size_t)na, sizeof *ka, cmp_u64);
      qsort(kb, (size_t)nb, sizeof *kb, cmp_u64);
      int64_t i = 0, j = 0, same = 0;
      while (i < na && j < nb) {
        if (ka[i] == kb[j]) {
          same++;
          i++;
          j++;
        } else if (ka[i] < kb[j])
          i++;
        else
          j++;
      }
      printf("   ДОЛЯ СРЕЗА ЗА КАДР (ходьба 1.4 м/с, 30 к/с, шаг %.4f м): было %lld, стало "
             "%lld, общих %lld, ИЗМЕНИЛОСЬ %.2f %% (верхняя оценка пользы, не цена — А732)\n",
             step_m, (long long)na, (long long)nb, (long long)same,
             100.0 * (double)(na + nb - 2 * same) / (double)(na > 0 ? na : 1));
      free(ka);
      free(kb);
      hz_slice_free(&B);
    }
    hz_slice_free(&A);
  }

  /* ---- 4г. РАЗРУШЕНИЕ (Ш4, Р8) ---- */
  if (hit) {
    /* Порог среза ЗАФИКСИРОВАН до прогона (А745): 1 пиксель, тот же, на котором
     * сняты числа §405. Менять его в этом шаге нельзя. */
    lodctx LH;
    memset(&LH, 0, sizeof LH);
    double eyec[3] = {g_eye[0], g_eye[1], g_eye[2]}, atc[3] = {g_at[0], g_at[1], g_at[2]};
    for (int a = 0; a < 3; a++)
      LH.eye[a] = (eyec[a] - fr.org[a]) / fr.h;
    LH.pxrad = (HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0) / 512.0;
    LH.thr = 1.0;

    LH.thrcell = g_lodceil;

    LH.thrpoly = g_lodpoly;

    hz_dcslice S0;
    if (hz_slice_init(&S0, lev) != HZ_DC_OK) exit(1);
    if (hz_slice_build(&S0, &T, &ht, lod_stop, &LH) != HZ_DC_OK) exit(1);

    /* ДВА СЛУЧАЯ ВРОЗЬ (А669), и оба называются: удар В УПОР — в метре перед
     * камерой по взгляду; удар ЗА СПИНОЙ — в метре позади. Одно число здесь
     * было бы подменой величины: доля зависит не от структуры, а от того, где
     * камера относительно удара. */
    double dir[3], dl = 0.0;
    for (int a = 0; a < 3; a++) {
      dir[a] = atc[a] - eyec[a];
      dl += dir[a] * dir[a];
    }
    dl = sqrt(dl);
    for (int a = 0; a < 3; a++)
      dir[a] /= dl;
    static const char *nm[2] = {"В УПОР", "ЗА СПИНОЙ"};
    for (int cse = 0; cse < 2; cse++) {
      /* ЦЕЛИТЬСЯ НАДО В ПОВЕРХНОСТЬ, А НЕ В ВОЗДУХ. Первая редакция ставила
       * сферу в метре по взгляду, попадала в пустоту, стирала НОЛЬ рёбер — и
       * все числа выходили нулями, неотличимыми от «правка не работает».
       * Здесь центр берётся в первой ЗАНЯТОЙ ячейке вдоль луча (шаг h/2 по
       * занятости Ш0): удар В УПОР — вперёд по взгляду, ЗА СПИНОЙ — назад. */
      double c[3], sgn = (cse == 0 ? 1.0 : -1.0);
      int found = 0;
      for (double s = 0.0; s < 20.0 && !found; s += fr.h * 0.5) {
        int32_t cc[3];
        int ok2 = 1;
        for (int a = 0; a < 3; a++) {
          double w = eyec[a] + dir[a] * sgn * s;
          double f = floor((w - fr.org[a]) / fr.h);
          if (!(f >= 0.0) || !(f < (double)fr.n)) ok2 = 0;
          cc[a] = ok2 ? (int32_t)f : 0;
          c[a] = w;
        }
        if (ok2 && hz_occ_get(P.b[lev], hz_occ_index(fr.n, cc[0], cc[1], cc[2]))) found = 1;
      }
      if (!found) {
        printf("   УДАР %s: луч не встретил геометрии — случай не измерен\n", nm[cse]);
        continue;
      }
      /* Копия таблицы: каждый случай бьёт по НЕТРОНУТОМУ полю. */
      hz_htab h2;
      if (hz_htab_init(&h2) != 0) exit(1);
      for (int32_t i = 0; i < ht.n; i++)
        if (hz_htab_add(&h2, ht.e[i].axis, ht.e[i].p, ht.e[i].t, ht.e[i].nrm, ht.e[i].in_lo) != 0)
          exit(1);
      hitstat HS;
      double ta = now_s();
      hit_sphere(&h2, &fr, c, 0.2, &HS);
      double t_edit = now_s() - ta;
      /* ПОЧИНКА: только затронутые ячейки, подъём к корню за O(глубины). */
      ta = now_s();
      int64_t ncell = 0;
      /* НЕГАТИВНЫЙ КОНТРОЛЬ nofix (§407): правка таблицы БЕЗ починки дерева.
       * Срез обязан остаться СТАРЫМ, то есть разойтись с полной пересборкой. */
      for (int32_t z = nofix ? 1 : HS.lo[2]; z <= (nofix ? 0 : HS.hi[2]); z++)
        for (int32_t y = HS.lo[1]; y <= HS.hi[1]; y++)
          for (int32_t x = HS.lo[0]; x <= HS.hi[0]; x++) {
            int32_t cl2[3] = {x, y, z};
            hz_dc_fix_cell(&T, &h2, cl2);
            ncell++;
          }
      double t_fix = now_s() - ta;
      ta = now_s();
      hz_dcslice S1;
      if (hz_slice_init(&S1, lev) != HZ_DC_OK) exit(1);
      if (hz_slice_build(&S1, &T, &h2, lod_stop, &LH) != HZ_DC_OK) exit(1);
      double t_slice2 = now_s() - ta;

      /* Г49 В ВЕРНОЙ ФОРМЕ (А744): сверка не с перестройкой ТОГО ЖЕ поддерева
       * (это тождество), а с ПОЛНОЙ пересборкой всего поля из правленой
       * таблицы — только она ловит недооценку затронутой области. */
      hz_dctree T2;
      if (hz_dc_init(&T2, lev) != HZ_DC_OK) exit(1);
      if (hz_dc_shape_occ(&T2, lev, u_occ, &P) != HZ_DC_OK) exit(1);
      hz_dc_masks_occ(&T2, &h2);
      hz_dc_forms_lazy(&T2, &h2);
      hz_dcslice S2;
      if (hz_slice_init(&S2, lev) != HZ_DC_OK) exit(1);
      if (hz_slice_build(&S2, &T2, &h2, lod_stop, &LH) != HZ_DC_OK) exit(1);
      int64_t nbit = 0;
      if (S1.n == S2.n)
        for (int32_t i = 0; i < S1.n; i++)
          if (memcmp(&S1.c[i], &S2.c[i], sizeof(hz_dccell)) != 0) nbit++;

      /* ДОЛЯ СРЕЗА, изменившаяся от удара: множествами по ключу. */
      int64_t na = S0.n, nb = S1.n, same = 0;
      uint64_t *ka = malloc((size_t)(na > 0 ? na : 1) * sizeof *ka);
      uint64_t *kb = malloc((size_t)(nb > 0 ? nb : 1) * sizeof *kb);
      if (ka == NULL || kb == NULL) exit(1);
      for (int64_t i = 0; i < na; i++)
        ka[i] = cellkey(&S0.c[i]);
      for (int64_t i = 0; i < nb; i++)
        kb[i] = cellkey(&S1.c[i]);
      qsort(ka, (size_t)na, sizeof *ka, cmp_u64);
      qsort(kb, (size_t)nb, sizeof *kb, cmp_u64);
      for (int64_t i = 0, j = 0; i < na && j < nb;) {
        if (ka[i] == kb[j]) {
          same++;
          i++;
          j++;
        } else if (ka[i] < kb[j])
          i++;
        else
          j++;
      }
      printf("   УДАР %s (сфера r=0.20 м): рёбер стёрто %lld, ЧАСТИЧНО накрытых %lld (%.1f %%), "
             "ячеек чинено %lld; ПРАВКА %.3f мс, ПОЧИНКА %.3f мс, СРЕЗ %.3f мс, всего %.3f мс; "
             "срез %lld -> %lld, ИЗМЕНИЛОСЬ %.3f %%; Г49 против ПОЛНОЙ пересборки: ячеек %d "
             "против %d, различий %lld\n",
             nm[cse], (long long)HS.nedge, (long long)HS.npart,
             100.0 * (double)HS.npart / (double)((HS.nedge + HS.npart) ? (HS.nedge + HS.npart) : 1),
             (long long)ncell, t_edit * 1e3, t_fix * 1e3, t_slice2 * 1e3,
             (t_edit + t_fix + t_slice2) * 1e3, (long long)na, (long long)nb,
             100.0 * (double)(na + nb - 2 * same) / (double)(na > 0 ? na : 1), S1.n, S2.n,
             (long long)nbit);
      free(ka);
      free(kb);
      hz_slice_free(&S1);
      hz_slice_free(&S2);
      hz_dc_free(&T2);
      hz_htab_free(&h2);
      /* Дерево испорчено починкой по правленой таблице — вернуть в исходное. */
      hz_dc_masks_occ(&T, &ht);
      hz_dc_forms_lazy(&T, &ht);
    }
    hz_slice_free(&S0);
  }

  /* ---- 4г2. ДОЛЯ ЛИСТЬЕВ СВИПА, ЗАТРОНУТЫХ УДАРОМ (Ф8'-0, §524) ---- */
  /* ВОПРОС, НА КОТОРЫЙ ЭТО ОТВЕЧАЕТ. §523 хочет нести по дереву свипа радианс по
   * направлениям. Возражение (конец §523): направленное поле есть состояние,
   * протянутое ВДОЛЬ ЛУЧЕЙ, и правка геометрии обесценивает не шар вокруг
   * удара, а ПУЧОК за ним. Пока не названа числом доля дерева, которую удар
   * пачкает, идти туда нельзя.
   *
   * ЧТО ИМЕННО МЕРИТСЯ. Открытость считается свипом ДО и ПОСЛЕ удара на двух
   * деревьях, построенных из ДВУХ пирамид; сличение идёт по ЯЧЕЙКАМ мелкой
   * сетки свипа — ровно по тому, что читает `sweep_vis`. Грязь копится
   * ОБЪЕДИНЕНИЕМ по всем образцам площадки: это множество, обесцененное ударом
   * за весь проход прямого света, а не за один образец.
   *
   * ДОЛЯ ПЕЧАТАЕТСЯ КАК ФУНКЦИЯ ДОПУСКА, А НЕ ПРИ ОДНОМ ПОРОГЕ (А921). Один
   * порог был бы магическим: моя первая редакция взяла `1/10` от `9.33 %`, а
   * `9.33 %` есть ДОЛЯ разошедшихся ячеек, то есть частота, а не величина
   * ошибки. Точка К88 называется явно — `1/10 × 1.799e−01 = 0.018` (§515,
   * максимум расхождения с маршем при пороге 32), — и читается по столбцу
   * `1e−2` как ближайшему снизу, то есть в сторону запаса.
   *
   * КАЙМА (А919). Цена инкрементного пересчёта есть доля ПОСЕЩЁННЫХ листьев, а
   * не изменившихся: чтобы узнать, что значение не поехало, лист надо посетить.
   * Кайма — листья, сами не затронутые, но соседние затронутым по `nb[6]`, —
   * считается тут же и даром, и тогда цена зажата с обеих сторон.
   *
   * ЧЕГО ЗДЕСЬ НЕТ (А926): энергии. Открытость безразмерна, `1/r²`, косинуса и
   * альбедо в ней нет, поэтому далёкий тёмный лист весит столько же, сколько
   * ярко освещённый. Значит число — ВЕРХНЯЯ оценка энергетически значимой
   * грязи. Для решения «идти или не идти» это сторона запаса; ЦЕНУ отсюда
   * брать нельзя. */
  if (hitsweep) {
    /* Допуски: 0 — точное неравенство (порога нет вовсе). */
    static const double TOL[HZ_HS_NTOL] = {0.0, 1e-3, 1e-2, 5e-2, 1e-1};
    double eyec[3] = {g_eye[0], g_eye[1], g_eye[2]}, atc[3] = {g_at[0], g_at[1], g_at[2]};
    double eyeg[3];
    for (int a = 0; a < 3; a++)
      eyeg[a] = (eyec[a] - fr.org[a]) / fr.h;
    /* КАМЕРА ДЕРЕВА СВИПА — ТА ЖЕ, ЧТО У ВЕТВИ `lit` (§514), и в тех же
     * единицах: сетка свипа вчетверо грубее поля. Иначе замер мерил бы другое
     * дерево, чем то, которое работает. */
    for (int a = 0; a < 3; a++)
      g_sweepeye[a] = eyeg[a];
    g_sweeppx = (HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0) / (double)resw;
    arealight AL;
    hall_light(&AL, &P, &fr, lo, hi, eyeg, 1);
    double su[HZ_LIGHT_SAMPLES], sv[HZ_LIGHT_SAMPLES];
    for (int a = 0; a < HZ_LIGHT_NS; a++)
      for (int b = 0; b < HZ_LIGHT_NS; b++) {
        su[a * HZ_LIGHT_NS + b] = lsamp(a);
        sv[a * HZ_LIGHT_NS + b] = lsamp(b);
      }
    int32_t gn = (int32_t)1 << (fr.lev - HZ_SWEEP_DROP);
    size_t gcells = (size_t)gn * (size_t)gn * (size_t)gn;
    double eyes[3];
    for (int a = 0; a < 3; a++)
      eyes[a] = eyeg[a] / (double)((int32_t)1 << HZ_SWEEP_DROP);
    printf("   Ф8'-0 ГРЯЗЬ ПОСЛЕ УДАРА: сфера r=%.2f м, сетка свипа %d^3, образцов площадки %d, "
           "порог дерева %.1f%s\n",
           hitrad, gn, HZ_LIGHT_SAMPLES, g_sweepthr, hitnoocc ? ", ЗАНЯТОСТЬ НЕ ПРАВИТСЯ" : "");
    /* Два случая врозь (А669): доля зависит не только от структуры, но и от
     * того, где камера относительно удара, — одно число было бы подменой. */
    double dir[3], dl = 0.0;
    for (int a = 0; a < 3; a++) {
      dir[a] = atc[a] - eyec[a];
      dl += dir[a] * dir[a];
    }
    dl = sqrt(dl);
    for (int a = 0; a < 3; a++)
      dir[a] /= dl;
    static const char *nmh[2] = {"В УПОР", "ЗА СПИНОЙ"};
    for (int cse = 0; cse < 2; cse++) {
      /* ЦЕЛИТЬСЯ В ПОВЕРХНОСТЬ, А НЕ В ВОЗДУХ — то же правило, что в 4г: первая
       * редакция там ставила сферу в метре по взгляду, попадала в пустоту и
       * давала одни нули, неотличимые от «правка не работает». */
      double c[3] = {0.0, 0.0, 0.0}, sgn = (cse == 0 ? 1.0 : -1.0), hdist = 0.0;
      int found = 0;
      for (double s = 0.0; s < 20.0 && !found; s += fr.h * 0.5) {
        hdist = s;
        int32_t cc[3];
        int ok2 = 1;
        for (int a = 0; a < 3; a++) {
          double w = eyec[a] + dir[a] * sgn * s;
          double f = floor((w - fr.org[a]) / fr.h);
          if (!(f >= 0.0) || !(f < (double)fr.n)) ok2 = 0;
          cc[a] = ok2 ? (int32_t)f : 0;
          c[a] = w;
        }
        if (ok2 && hz_occ_get(P.b[lev], hz_occ_index(fr.n, cc[0], cc[1], cc[2]))) found = 1;
      }
      if (!found) {
        printf("      УДАР %s: луч не встретил геометрии — случай не измерен\n", nmh[cse]);
        continue;
      }
      /* КОПИЯ ПИРАМИДЫ: каждый случай бьёт по НЕТРОНУТОМУ полю, и исходная
       * пирамида нужна дальше ветви `lit`. */
      opyr P2;
      P2.lev = P.lev;
      for (int l = 0; l <= P.lev; l++) {
        int32_t n = (int32_t)1 << l;
        size_t nb2 = hz_occ_bytes((size_t)n * (size_t)n * (size_t)n);
        P2.b[l] = malloc(nb2);
        if (P2.b[l] == NULL) exit(1);
        memcpy(P2.b[l], P.b[l], nb2);
      }
      int64_t nclr = 0, npart = 0;
      double ta = now_s();
      /* `hitnoocc` — ПРОВЕРКА НА ЛОЖНЫЙ НОЛЬ, а не негативный контроль (А10):
       * пирамида копируется, но не правится, дерево строится заново из копии, и
       * доля обязана выйти РОВНО нулём по ТОЧНОМУ неравенству. Не выйдет —
       * точный критерий негоден сам по себе, и читать надо столбцы с допуском. */
      if (!hitnoocc) nclr = hit_occ(&P2, &fr, c, hitrad, &npart);
      double t_carve = now_s() - ta;

      stree TA, TB;
      ta = now_s();
      stree_build(&TA, &P, fr.lev, HZ_SWEEP_DROP, gn, eyes, g_sweeppx, g_sweepthr);
      stree_links(&TA, gn);
      stree_build(&TB, &P2, fr.lev, HZ_SWEEP_DROP, gn, eyes, g_sweeppx, g_sweepthr);
      stree_links(&TB, gn);
      double t_build = now_s() - ta;

      /* Грязь по ЯЧЕЙКАМ: бит на допуск. Копится объединением по образцам. */
      unsigned char *dirty = calloc(gcells, 1);
      if (dirty == NULL) exit(1);
      double dmaxv = 0.0;
      ta = now_s();
      for (int sm = 0; sm < HZ_LIGHT_SAMPLES; sm++) {
        double q[3];
        for (int k = 0; k < 3; k++)
          q[k] = AL.c[k] + AL.u[k] * su[sm] + AL.v[k] * sv[sm];
        tsweep_light(&TA, &P, &fr, fr.lev, HZ_SWEEP_DROP, q);
        tsweep_light(&TB, &P2, &fr, fr.lev, HZ_SWEEP_DROP, q);
        for (size_t ci = 0; ci < gcells; ci++) {
          float a2 = TA.open[TA.idx[ci]], b2 = TB.open[TB.idx[ci]];
          /* ПОБИТОВО, а не по допуску: нижняя строка таблицы есть «порога нет
           * вовсе», и сравнение представлений — точно то, что это значит. */
          if (memcmp(&a2, &b2, sizeof a2) == 0) continue;
          double d = fabs((double)a2 - (double)b2);
          if (d > dmaxv) dmaxv = d;
          unsigned m2 = 0;
          for (int t = 0; t < HZ_HS_NTOL; t++)
            if (d > TOL[t]) m2 |= 1u << t;
          dirty[ci] |= (unsigned char)m2;
        }
      }
      double t_sweep = now_s() - ta;

      /* Лист затронут, если накрывает хотя бы одну грязную ячейку. Считается по
       * обоим деревьям: топология у них разная, и доля от этого зависит. */
      unsigned char *lfa = calloc((size_t)TA.n, 1), *lfb = calloc((size_t)TB.n, 1);
      if (lfa == NULL || lfb == NULL) exit(1);
      for (size_t ci = 0; ci < gcells; ci++) {
        unsigned char m2 = dirty[ci];
        if (m2 == 0) continue;
        lfa[TA.idx[ci]] |= m2;
        lfb[TB.idx[ci]] |= m2;
      }
      int64_t na2[HZ_HS_NTOL], nb2c[HZ_HS_NTOL], ncell2[HZ_HS_NTOL], rim[HZ_HS_NTOL];
      for (int t = 0; t < HZ_HS_NTOL; t++)
        na2[t] = nb2c[t] = ncell2[t] = rim[t] = 0;
      for (size_t ci = 0; ci < gcells; ci++)
        for (int t = 0; t < HZ_HS_NTOL; t++)
          if (dirty[ci] & (1u << t)) ncell2[t]++;
      for (int32_t i = 0; i < TA.n; i++)
        for (int t = 0; t < HZ_HS_NTOL; t++)
          if (lfa[i] & (1u << t)) na2[t]++;
      for (int32_t i = 0; i < TB.n; i++)
        for (int t = 0; t < HZ_HS_NTOL; t++)
          if (lfb[i] & (1u << t)) nb2c[t]++;
      /* КАЙМА (А919): лист НЕ затронут, но сосед по одной из шести ссылок —
       * затронут. Это фронт, который пришлось бы посетить, чтобы убедиться, что
       * дальше идти не надо; вместе с затронутыми он и есть цена. */
      for (int32_t i = 0; i < TB.n; i++) {
        if (TB.nd[i].child0 >= 0) continue;
        for (int t = 0; t < HZ_HS_NTOL; t++) {
          if (lfb[i] & (1u << t)) continue;
          for (int f = 0; f < 6; f++) {
            int32_t nj = TB.nd[i].nb[f];
            if (nj >= 0 && (lfb[nj] & (1u << t))) {
              rim[t]++;
              break;
            }
          }
        }
      }
      /* РАССТОЯНИЕ ОТ КАМЕРЫ ДО УДАРА ПЕЧАТАЕТСЯ, А НЕ ПОДРАЗУМЕВАЕТСЯ. Имя «в
       * упор» говорит о НАПРАВЛЕНИИ (вперёд по взгляду), а не о близости: сфера
       * ставится в ПЕРВУЮ ЗАНЯТУЮ ячейку вдоль луча, и она может оказаться в
       * нескольких метрах. Без этого числа доля при пороге дерева читалась бы
       * неверно — дробление задано УГЛОВЫМ размером, то есть расстоянием. */
      printf("      УДАР %s (%.2f м от камеры): ячеек занятости СТЁРТО %lld, частично накрытых "
             "%lld; листьев ДО %d, ПОСЛЕ %d; макс |Дельта открытости| %.4e; правка %.1f мс, "
             "постройка двух деревьев %.1f мс, %d свипов %.0f мс\n",
             nmh[cse], hdist, (long long)nclr, (long long)npart, TA.nleaf, TB.nleaf, dmaxv,
             t_carve * 1e3, t_build * 1e3, 2 * HZ_LIGHT_SAMPLES, t_sweep * 1e3);
      printf("         допуск | лист ДО | лист ПОСЛЕ | кайма ПОСЛЕ | ячейки сетки\n");
      for (int t = 0; t < HZ_HS_NTOL; t++)
        printf("         %-6.0e | %6.3f %% | %8.3f %% | %9.3f %% | %8.4f %%%s\n", TOL[t],
               100.0 * (double)na2[t] / (double)(TA.nleaf ? TA.nleaf : 1),
               100.0 * (double)nb2c[t] / (double)(TB.nleaf ? TB.nleaf : 1),
               100.0 * (double)rim[t] / (double)(TB.nleaf ? TB.nleaf : 1),
               100.0 * (double)ncell2[t] / (double)gcells,
               t == 2 ? "   <- К88: 1/10 x 1.799e-01 = 0.018, читается отсюда" : "");
      free(dirty);
      free(lfa);
      free(lfb);
      stree_free(&TA);
      stree_free(&TB);
      opyr_free(&P2);
    }
  }

  /* ---- 4д. КАДР СО СВЕТОМ (Ш5, часть первая) ---- */
  if (lit) {
    lodctx LL;
    memset(&LL, 0, sizeof LL);
    double eyec[3] = {g_eye[0], g_eye[1], g_eye[2]}, atc[3] = {g_at[0], g_at[1], g_at[2]},
           upc[3] = HZ_CFG_UP;
    /* ЦИКЛ ХОДЬБЫ (§589). Без ключа `walk` тело исполняется РОВНО ОДИН РАЗ и
     * кончается `break` внизу — прежнее поведение слово в слово. С ключом камера
     * правится в конце каждого прохода, и всё, что от неё зависит (срез, свет,
     * растр), считается заново; всё, что не зависит (дерево, рёбра, текстуры,
     * таблица треугольников ячейки), остаётся построенным. */
    /* ЯДРО КОСВЕННОГО СВЕТА СЧИТАЕТСЯ ЗДЕСЬ — ОДИН РАЗ, ДО ЦИКЛА (§597).
     * Камеры в нём нет, поэтому оно переживает и кадр, и ходьбу. Источник
     * обязан быть камеронезависимым: под `sun` это проверено чтением (А1033),
     * с площадным `hall_light` ядро несовместимо и потому не строится. */
    if (g_hgather > 0.0 && !g_indslice) {
      if (!g_sun && !g_keonly) {
        printf("   §597 ЯДРО НЕ СТРОИТСЯ: источник площадной, а `hall_light` ищет потолок "
               "СПУСКОМ ОТ КАМЕРЫ — сценно закреплённое ядро с ним несовместимо. Нужен `sun=` "
               "или `keonly` (§618: светятся сами поверхности).\n");
      } else {
        arealight ALc;
        hall_light(&ALc, &P, &fr, lo, hi, LL.eye, 0);
        ind_core_build(&T, &ht, &fr, &P, &ALc, &m, &CT, lev, g_hgather, indvis);
        xcmp_report(&T, lev, &fr);
      }
    }
    int walk_alive = 1;
    if (g_walk && !walk_open(resw, resh)) return 2;
    double t_prev = 1.0 / 30.0; /* первый шаг движения — как при 30 к/с */
    for (;;) {
      /* СЕКУНДОМЕР НА ВСЮ ИТЕРАЦИЮ, А НЕ НА ТРИ НАЗВАННЫЕ СТАДИИ (§590).
       * Первая редакция ходилки печатала `срез + свет + растр` и называла это
       * кадром. Это ПОДМЕНА ВЕЛИЧИНЫ: в итерации есть работа помимо этих трёх, и
       * пользователь увидел `0.2` к/с там, где строка обещала `5`. Разность
       * теперь печатается ОТДЕЛЬНОЙ статьёй «прочее» — величина, которую никто
       * не мерил, обязана быть видна, а не растворяться. */
      double t_iter0 = now_s();
      for (int a = 0; a < 3; a++)
        LL.eye[a] = (eyec[a] - fr.org[a]) / fr.h;
      LL.pxrad = (HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0) / (double)resw;
      /* Ф5. (§514): камера дерева свипа — та же, что у среза, и в тех же единицах.
       * Разные камеры у тени и у поверхности дали бы несогласованную подробность. */
      for (int a = 0; a < 3; a++)
        g_sweepeye[a] = LL.eye[a];
      g_sweeppx = LL.pxrad;
      LL.thr = lodthr;

      LL.thrcell = g_lodceil;

      LL.thrpoly = g_lodpoly;
      hz_dcslice S;
      if (hz_slice_init(&S, lev) != HZ_DC_OK) exit(1);
      double ta = now_s();
      /* ПЕЧЬ СТРОИТСЯ НА ПОЛНУЮ ГЛУБИНУ, А НЕ ПО СРЕЗУ. Коробка вся ПЛОСКАЯ,
       * невязка QEF на ней ноль, и критерий LOD законно огрубляет её до предела:
       * при камере зала срез вышел в СЕМЬ ячеек, и печь мерила перенос между
       * семью гигантскими площадками. Это был не отказ переноса, а отказ моего
       * замера — величина считалась не на том. */
      if (hz_slice_build(&S, &T, &ht, (oven > 0.0 || plates > 0.0) ? NULL : lod_stop, &LL) !=
          HZ_DC_OK)
        exit(1);
      double t_slice = now_s() - ta;

      /* ИСТОЧНИК: площадка под потолком зала. Габарит сцены известен, потолок —
       * его верх по оси Y; площадка ставится на 10 см ниже и в центре плана.
       * Числа не магические: они выведены из ГАБАРИТА, а не подобраны на глаз.
       * Правило вынесено в `hall_light` (§524): им же пользуется замер грязи. */
      arealight AL;
      hall_light(&AL, &P, &fr, lo, hi, LL.eye, 1);

      /* МАТЕРИАЛ ЯЧЕЙКИ СРЕЗА (Ш5б, §430). Байт `mat` перестаёт быть резервом:
       * берётся материал первого треугольника в ячейке (`celltris` уже построен).
       * Индекс, а не альбедо: в срезе лежит ИНДЕКС, полезная нагрузка — в таблице
       * (Р4). Материалов на сцене десятки, таблица горяча в кэше.
       * ЧЕГО ЭТО НЕ ДЕЛАЕТ: на границе двух материалов ячейка получает ОДИН из
       * них, а не оба — граница пройдёт по ячейкам среза, то есть с точностью
       * LOD. Для отскока это законно (энергия), для резкой границы текстуры —
       * нет; текстур пока и нет. */
      float *uvs = (m.vt != NULL && m.ft != NULL) ? calloc(2 * (size_t)S.n, sizeof *uvs) : NULL;
      double t_mat = now_s();
      {
        int64_t nmat = 0, nuv = 0;
        for (int32_t i = 0; i < S.n; i++) {
          /* ЧИТАТЬ ПО ВЕРШИНЕ, А НЕ ПО УГЛУ (найдено 08-11 картинкой с
           * материалами). У крупной ячейки среза нижний угол лежит ГДЕ УГОДНО —
           * внутри тела, в пустоте, на соседнем предмете, — и материал оттуда
           * либо не находится вовсе, либо берётся чужой. Замерено: назначено
           * `6 084` ячейкам из `19 345` (31 %), а большой шар вышел пятнистым,
           * потому что часть его ячеек получила материал соседнего.
           * ВЕРШИНА ЖЕ ЛЕЖИТ НА ПОВЕРХНОСТИ по построению DC, и её ячейка
           * содержит ровно тот треугольник, который эту вершину и породил.
           * Это тот же класс, что чтение радианса по углу (§486). */
          double vw[3];
          hz_slice_vertex(&S, i, vw);
          int32_t cell[3];
          for (int a = 0; a < 3; a++) {
            double f = floor(vw[a]);
            if (f < 0.0) f = 0.0;
            if (f > (double)(fr.n - 1)) f = (double)(fr.n - 1);
            cell[a] = (int32_t)f;
          }
          const int32_t *ls = NULL;
          if (ct_list(&CT, cell, &ls) == 0) continue;
          int32_t mi = m.fm != NULL ? m.fm[ls[0]] : 0;
          if (mi < 0 || mi >= m.nmtl) mi = 0;
          S.c[i].mat = (uint8_t)(mi < 255 ? mi : 255);
          nmat++;
          /* Ш8 (§575): КООРДИНАТА ТЕКСТУРЫ БЕРЁТСЯ ОТТУДА ЖЕ, ОТКУДА МАТЕРИАЛ —
           * по ТОМУ ЖЕ треугольнику `ls[0]`, барицентрикой в точке вершины DC.
           * Вершина лежит НА поверхности по построению, поэтому проекции не
           * нужно; барицентрика считается в плоскости треугольника, и при выходе
           * за него (вершина ячейки чуть в стороне) координаты не отбрасываются,
           * а ЗАЖИМАЮТСЯ — иначе край поверхности остался бы без текстуры. */
          if (uvs != NULL && m.ft != NULL && m.vt != NULL) {
            const double *A2, *B2, *C2;
            tri_verts(&m, ls[0], &A2, &B2, &C2);
            double e1[3], e2[3], vp[3];
            for (int k = 0; k < 3; k++) {
              e1[k] = B2[k] - A2[k];
              e2[k] = C2[k] - A2[k];
              vp[k] = fr.org[k] + vw[k] * fr.h - A2[k];
            }
            double d11 = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
            double d12 = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
            double d22 = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
            double dp1 = vp[0] * e1[0] + vp[1] * e1[1] + vp[2] * e1[2];
            double dp2 = vp[0] * e2[0] + vp[1] * e2[1] + vp[2] * e2[2];
            double dn = d11 * d22 - d12 * d12;
            if (fabs(dn) > 0.0) {
              double bu = (d22 * dp1 - d12 * dp2) / dn;
              double bv = (d11 * dp2 - d12 * dp1) / dn;
              if (bu < 0.0) bu = 0.0;
              if (bv < 0.0) bv = 0.0;
              if (bu + bv > 1.0) {
                double s2 = bu + bv;
                bu /= s2;
                bv /= s2;
              }
              int32_t q0 = m.ft[3 * (size_t)ls[0] + 0], q1 = m.ft[3 * (size_t)ls[0] + 1],
                      q2 = m.ft[3 * (size_t)ls[0] + 2];
              if (q0 >= 0 && q1 >= 0 && q2 >= 0 && q0 < m.nvt && q1 < m.nvt && q2 < m.nvt) {
                for (int k = 0; k < 2; k++)
                  uvs[2 * (size_t)i + (size_t)k] = (float)(m.vt[2 * (size_t)q0 + (size_t)k] +
                                                           bu * (m.vt[2 * (size_t)q1 + (size_t)k] -
                                                                 m.vt[2 * (size_t)q0 + (size_t)k]) +
                                                           bv * (m.vt[2 * (size_t)q2 + (size_t)k] -
                                                                 m.vt[2 * (size_t)q0 + (size_t)k]));
                nuv++;
              }
            }
          }
        }
        if (!g_walk)
          printf("   МАТЕРИАЛ В СРЕЗЕ: назначен %lld ячейкам из %d, материалов в сцене %d; "
                 "КООРДИНАТА ТЕКСТУРЫ у %lld (%.1f %%)\n",
                 (long long)nmat, S.n, m.nmtl, (long long)nuv,
                 100.0 * (double)nuv / (double)(S.n ? S.n : 1));
      }
      t_mat = now_s() - t_mat;

      float *irr = malloc(3 * (size_t)S.n * sizeof *irr);
      if (irr == NULL) exit(1);
      ta = now_s();
      /* РАБОЧИЙ ПУТЬ — С ПОДЪЁМОМ (§426): тот же предикат, вчетверо дешевле.
       * Плоский марш остаётся АРБИТРОМ и зовётся ниже. */
      int64_t nstep_w = 0;
      /* §618: при `keonly` внешнего источника НЕТ и в кадре — `b0` берётся из `Ke`
       * тех же материалов, что и в ядре. Иначе кадр светился бы площадкой, а
       * косвенный приходил бы от светящихся тел: две разные сцены в одном кадре. */
      if (g_keonly) {
        memset(irr, 0, 3 * (size_t)S.n * sizeof *irr);
        if (!g_nokemit)
          for (int32_t i2 = 0; i2 < S.n; i2++) {
            const double *ke = m.mtl[S.c[i2].mat < m.nmtl ? S.c[i2].mat : 0].ke3;
            for (int k2 = 0; k2 < 3; k2++)
              irr[3 * (size_t)i2 + (size_t)k2] += (float)ke[k2];
          }
      } else
        front_direct(&S, &fr, &P, &AL, irr, 0.5, 1, &nstep_w, &m);
      double t_dir = now_s() - ta;
      /* КОСВЕННЫЙ СВЕТ ИЗ УЗЛОВ (§597, Р5). Сбора в кадре НЕТ: ячейка среза
       * находит свой узел спуском O(глубины) и ЧИТАЕТ готовое значение.
       * Отказы спуска СЧИТАЮТСЯ — молчаливый промах выглядел бы лёгким шумом
       * (А1029), а не отказом. */
      double t_ind = 0.0;
      int64_t nindbad = 0;
      double s_dir = 0.0, s_ind = 0.0;
      if (g_indnode != NULL) {
        double ti0 = now_s();
        for (int32_t i = 0; i < S.n; i++) {
          int32_t ni = node_of_cell(&T, lev, &S.c[i]);
          if (ni < 0 || ni >= g_indnode_n) {
            nindbad++;
            continue;
          }
          for (int a = 0; a < 3; a++) {
            s_dir += (double)irr[3 * (size_t)i + (size_t)a];
            s_ind += (double)g_indnode[ni][a];
            irr[3 * (size_t)i + (size_t)a] += g_indnode[ni][a];
          }
        }
        t_ind = now_s() - ti0;
        /* А1029: отказы спуска ПЕЧАТАЮТСЯ — ноль обязателен, иначе раскладка
         * неверна, и это видно числом, а не глазом.
         * А1032: отношение косвенного к прямому — ловушка на ДВОЙНОЙ УЧЁТ.
         * Печатается один раз (первый кадр), чтобы не засорять цикл ходьбы. */
        static int said = 0;
        if (!said) {
          said = 1;
          printf("   §597 КОСВЕННЫЙ ИЗ УЗЛОВ: %.1f мс на %d ячеек, ОТКАЗОВ СПУСКА %lld; "
                 "Σ E_ind / Σ E_dir = %.4f\n",
                 t_ind * 1e3, S.n, (long long)nindbad, s_ind / (s_dir > 0.0 ? s_dir : 1.0));
        }
      }
      /* §758: КАДР РАЗВЁРТКОЙ. Сличение ДО растра в линейном пространстве
       * кадрового среза (мера кадровая — А1188), затем ПОДМЕНА irr на
       * свиповое E·alb + Ke. Покрытие и отказы — счётчиками (А891). */
      if (xframe && g_swEn != NULL) {
        int64_t ncov8 = 0, nmiss8 = 0, nz8 = 0, nuse8 = 0, nout8 = 0;
        double *rt8 = malloc((size_t)S.n * sizeof *rt8);
        double *tmp8 = malloc((size_t)S.n * sizeof *tmp8);
        float *swv = malloc(3 * (size_t)S.n * sizeof *swv);
        if (rt8 == NULL || tmp8 == NULL || swv == NULL) exit(1);
        int64_t nb8 = 0;
        /* §760: дыры закрываются глубочайшим предком с полем; распределение
         * глубины подъёма — встроенный НК (подъём к корню = заливка) */
        int64_t nup1 = 0, nup2 = 0, nup3 = 0, nuproot = 0, nupall = 0;
        /* §798: читаемое Ei (до альбедо) запоминается для прибора кадровой
         * чешуи; −1 — непокрыто */
        double *eall8 = malloc((size_t)(S.n > 0 ? S.n : 1) * sizeof *eall8);
        if (eall8 == NULL) exit(1);
        for (int32_t i = 0; i < S.n; i++) {
          int uplev8 = 0;
          int32_t ni = node_swE_best(&T, lev, &S.c[i], &uplev8, g_swEd);
          double Ei = (ni >= 0) ? g_swEn[ni] / g_swEd[ni] : -1.0;
          /* §798: гибрид — прямой свет ВТОРЫМ каналом из ядрового сбора */
          if (g_fcdn != NULL) {
            int updir8 = 0;
            int32_t nd8 = node_swE_best(&T, lev, &S.c[i], &updir8, g_fcdd);
            if (nd8 >= 0) Ei = (Ei > 0.0 ? Ei : 0.0) + g_fcdn[nd8] / g_fcdd[nd8];
          }
          eall8[i] = Ei;
          if (Ei < 0.0) {
            nmiss8++;
            for (int k = 0; k < 3; k++)
              swv[3 * (size_t)i + (size_t)k] = 0.0f;
            continue;
          }
          if (uplev8 > 0) {
            nupall++;
            if (uplev8 == 1)
              nup1++;
            else if (uplev8 == 2)
              nup2++;
            else
              nup3++;
            if (uplev8 >= (int)S.c[i].lvl - 2) nuproot++;
          }
          ncov8++;
          const double *ke8 = m.mtl[S.c[i].mat < m.nmtl ? S.c[i].mat : 0].ke3;
          double bsum = 0.0;
          for (int k = 0; k < 3; k++) {
            double v = Ei * alb(&m, S.c[i].mat, k) + (!g_nokemit ? ke8[k] : 0.0);
            swv[3 * (size_t)i + (size_t)k] = (float)v;
            bsum += v;
          }
          if (bsum > 0.0) tmp8[nb8++] = bsum;
        }
        double thr8 = 0.0;
        if (nb8 > 0) {
          qsort(tmp8, (size_t)nb8, sizeof *tmp8, cmp_d);
          thr8 = 1e-3 * tmp8[nb8 / 2]; /* порог от данных, §559 */
        }
        for (int32_t i = 0; i < S.n; i++) {
          double a8 = 0.0, b8 = 0.0;
          for (int k = 0; k < 3; k++) {
            a8 += (double)irr[3 * (size_t)i + (size_t)k];
            b8 += (double)swv[3 * (size_t)i + (size_t)k];
          }
          if (!(b8 > thr8)) {
            if (a8 > 0.0) nz8++;
            continue;
          }
          rt8[nuse8++] = a8 / b8;
        }
        if (nuse8 > 0) {
          memcpy(tmp8, rt8, (size_t)nuse8 * sizeof *tmp8);
          qsort(tmp8, (size_t)nuse8, sizeof *tmp8, cmp_d);
          double med8 = tmp8[nuse8 / 2];
          for (int64_t q8i = 0; q8i < nuse8; q8i++)
            if (rt8[q8i] < 0.5 * med8 || rt8[q8i] > 2.0 * med8) nout8++;
          printf("   §758 КАДР-СЛИЧЕНИЕ (ядро/свип, полный свет, кадровый срез): покрытие %lld "
                 "из %d (отказов %lld, свип~0 при ядре>0: %lld); пар %lld, медиана %.4g, p10 "
                 "%.4g, p90 %.4g; ВНЕ x2: %lld (%.1f %%)\n",
                 (long long)ncov8, S.n, (long long)nmiss8, (long long)nz8, (long long)nuse8, med8,
                 tmp8[(nuse8 * 10) / 100], tmp8[(nuse8 * 90) / 100], (long long)nout8,
                 100.0 * (double)nout8 / (double)nuse8);
        }
        printf("   §760 ДЫРЫ: закрыто подъёмом %lld (1 ур. %lld, 2 ур. %lld, 3+ %lld; с "
               "почти-корня %lld = %.1f %% дыр); осталось непокрытых %lld\n",
               (long long)nupall, (long long)nup1, (long long)nup2, (long long)nup3,
               (long long)nuproot, 100.0 * (double)nuproot / (double)(nupall > 0 ? nupall : 1),
               (long long)nmiss8);
        /* ---- §798: ПРИБОР КАДРОВОЙ ЧЕШУИ — CV читаемого Ei по квантованным
         * плоскостям КАДРОВОГО среза (группировка и константы §796-прибора;
         * значение — Ei до альбедо, канал монохромный). Печатается во всех
         * xframe-режимах, метка называет источник прямого канала. */
        {
          enum { HZ_CV_NQ = 4, HZ_CV_OFFH = 4, HZ_CV_MIN = 8 };
          cvpair796 *cp8 = malloc((size_t)(S.n > 0 ? S.n : 1) * sizeof *cp8);
          if (cp8 == NULL) exit(1);
          int32_t ncp8 = 0;
          for (int32_t i = 0; i < S.n; i++) {
            if (!(eall8[i] >= 0.0)) continue;
            double pw8[3], nn8[3];
            hz_slice_vertex(&S, i, pw8);
            for (int a = 0; a < 3; a++)
              pw8[a] = fr.org[a] + pw8[a] * fr.h;
            hz_slice_normal(&S, i, nn8);
            double off = 0.0;
            uint64_t kq = 0;
            for (int a = 0; a < 3; a++) {
              long qn = lround(nn8[a] * (double)HZ_CV_NQ);
              kq = (kq << 8) | (uint64_t)(uint8_t)(qn + 16);
              off += nn8[a] * pw8[a];
            }
            long ob = (long)floor(off / ((double)HZ_CV_OFFH * fr.h));
            cp8[ncp8].key = (kq << 32) | (uint64_t)(uint32_t)(int32_t)ob;
            cp8[ncp8].e = i;
            ncp8++;
          }
          qsort(cp8, (size_t)ncp8, sizeof *cp8, cmp_cvpair796);
          double *cvv8 = malloc((size_t)(ncp8 > 0 ? ncp8 : 1) * sizeof *cvv8);
          if (cvv8 == NULL) exit(1);
          int32_t ng8 = 0, ngb8 = 0;
          int32_t i0 = 0;
          while (i0 < ncp8) {
            int32_t i1 = i0;
            while (i1 < ncp8 && cp8[i1].key == cp8[i0].key)
              i1++;
            ng8++;
            if (i1 - i0 >= HZ_CV_MIN) {
              double s1 = 0.0, s2 = 0.0;
              for (int32_t j = i0; j < i1; j++) {
                double v = eall8[cp8[j].e];
                s1 += v;
                s2 += v * v;
              }
              double nn6 = (double)(i1 - i0);
              double mean = s1 / nn6;
              double var = s2 / nn6 - mean * mean;
              if (mean > 0.0 && var > 0.0) cvv8[ngb8++] = sqrt(var) / mean;
            }
            i0 = i1;
          }
          if (ngb8 > 0) {
            qsort(cvv8, (size_t)ngb8, sizeof *cvv8, cmp_dev699);
            printf("   §798 ЧЕШУЯ КАДРА (%s): групп %d (с >=%d ячейками %d); CV p50 %.4f, "
                   "p90 %.4f, макс %.4f\n",
                   g_fcmode == 2 ? "гибрид: прямой из ядра"
                                 : (g_fcmode == 1 ? "xfcelem: прямой из элементов" : "без fc"),
                   ng8, (int)HZ_CV_MIN, ngb8, cvv8[ngb8 / 2],
                   cvv8[(int32_t)((int64_t)ngb8 * 9 / 10)], cvv8[ngb8 - 1]);
          } else
            printf("   §798 ЧЕШУЯ КАДРА: групп с >=%d ячейками НЕТ (%d групп)\n", (int)HZ_CV_MIN,
                   ng8);
          free(cvv8);
          free(cp8);
        }
        free(eall8);
        memcpy(irr, swv, 3 * (size_t)S.n * sizeof *irr);
        free(rt8);
        free(tmp8);
        free(swv);
      }
      /* ---- ДИАГНОСТИКА КАДРА: ВСЁ, ЧТО НИЖЕ, К КАРТИНКЕ НЕ ОТНОСИТСЯ (§589) ----
       * Свип, арбитр-марш, эталоны, поячеечные сличения, иерархический отскок —
       * это ЗАМЕРЫ, а не кадр. Одним прогоном они стоят секунды и потому в цикле
       * ходьбы недопустимы: `walk` их выключает целиком.
       * ЧТО ЭТО ЗНАЧИТ ЧЕСТНО: в ходьбе НЕТ КОСВЕННОГО СВЕТА. Виден прямой
       * солнечный плюс тень; отскок (`3.4` с) в кадровый бюджет не входит и здесь
       * не считается вовсе. Это не «упрощение ради скорости», а прямое следствие
       * того, что отскок стоит втрое дороже секунды. */
      float *irr2 = NULL;
      if (!g_walk) {
        /* РЕШЁТКА НА ПОТОЛКЕ — ЧИСЛОМ, А НЕ ГЛАЗОМ (§541; К19 запрещает мерить
         * картинкой). Берутся ячейки среза с нормалью ВНИЗ и лежащие ПОД САМОЙ
         * ПЛОЩАДКОЙ (А961): по всему потолку метрика смешала бы решётку с законным
         * краевым спадом. `max/p50` при точечных источниках велик — каждая проба
         * ставит пятно; у настоящей площадки потолок под ней не освещён вовсе
         * (`cos θ_s = 0` у односторонней панели), и метрика вырождается. */
        {
          double sde = 0.0;
          int64_t nlit = 0;
          for (int32_t i = 0; i < S.n; i++) {
            for (int k = 0; k < 3; k++)
              sde += (double)irr[3 * (size_t)i + (size_t)k];
            if (irr[3 * (size_t)i] > 0.0f) nlit++;
          }
          printf(
              "   §541 ПРЯМОЙ СВЕТ АБСОЛЮТНО: Σ E_dir по срезу %.6e, освещённых ячеек %lld из %d\n",
              sde, (long long)nlit, S.n);
        }
        {
          double *ce = malloc((size_t)(S.n > 0 ? S.n : 1) * sizeof *ce);
          if (ce == NULL) exit(1);
          int64_t nce = 0;
          double au = sqrt(AL.u[0] * AL.u[0] + AL.u[1] * AL.u[1] + AL.u[2] * AL.u[2]);
          double av = sqrt(AL.v[0] * AL.v[0] + AL.v[1] * AL.v[1] + AL.v[2] * AL.v[2]);
          for (int32_t i = 0; i < S.n; i++) {
            double nn4[3], pw4[3];
            hz_slice_normal(&S, i, nn4);
            if (!(nn4[1] < -0.9)) continue; /* нормаль вниз — это потолок */
            hz_slice_vertex(&S, i, pw4);
            for (int k = 0; k < 3; k++)
              pw4[k] = fr.org[k] + pw4[k] * fr.h;
            if (fabs(pw4[0] - AL.c[0]) > au || fabs(pw4[2] - AL.c[2]) > av) continue;
            ce[nce++] = (double)irr[3 * (size_t)i];
          }
          if (nce > 0) {
            qsort(ce, (size_t)nce, sizeof *ce, cmp_d);
            double p50 = ce[nce / 2], mx = ce[nce - 1];
            printf("   §541 РЕШЁТКА НА ПОТОЛКЕ ЧИСЛОМ: ячеек под площадкой %lld, E_dir p50 %.4e, "
                   "max %.4e, max/p50 = %.2f\n",
                   (long long)nce, p50, mx, mx / (p50 > 0.0 ? p50 : 1.0));
          } else
            printf("   §541 РЕШЁТКА НА ПОТОЛКЕ: ячеек под площадкой НЕТ — метрика не снята\n");
          free(ce);
        }
        /* Ш15 (§484): КАРТИНКА ИЗ РАЗВЁРТКИ. Если сработал `xsweep`, облучённость
         * берётся не маршем и не гатером, а ИСХОДЯЩИМ РАДИАНСОМ поверхностных
         * элементов, посчитанным уравнением переноса. Заслоны там по построению —
         * то, что маршем стоит 65 с (§474).
         * СВЯЗЬ ЯЧЕЕК ПРЯМАЯ: и срез, и сетка переноса адресуются одной сеткой
         * уровня `lev`, поэтому переклад — это чтение по координате, а не поиск. */
        /* §798/А1307: ПРИ xframe ЭТА ПЕРЕЗАПИСЬ ПРОПУСКАЕТСЯ. Ш15-путь кладёт
         * в irr ИСХОДЯЩИЙ РАДИАНС sout (константный член, максимум по
         * элементам клетки, псевдоцвет ×1/0.92/0.78) и стоит НИЖЕ подмены
         * §758 — поэтому «кадр развёрткой» с §758 по §797 рисовался из sout,
         * а канал E·alb + Ke до пикселей не доходил вовсе (вскрыто §798:
         * гибрид/НК дали ПОБИТОВО один кадр при разных Ei). xframe читает
         * узловые каналы; без xframe Ш15-картинка остаётся как была. */
        if (xrad != NULL && !(xframe && g_swEn != NULL)) {
          int64_t nfound = 0;
          for (int32_t i = 0; i < S.n; i++) {
            /* Ячейка читает СВОЙ уровень пирамиды, а не угол на самом мелком. НК
             * `xcorner` возвращает прежнее чтение. */
            int lv = xcorner ? lev : (int)S.c[i].lvl;
            int32_t nl = (int32_t)1 << lv;
            int sh = lev - lv;
            size_t k = hz_occ_index(nl, S.c[i].lo[0] >> sh, S.c[i].lo[1] >> sh, S.c[i].lo[2] >> sh);
            const float *src = xradlv != NULL ? xradlv[lv] : xrad;
            for (int c = 0; c < 3; c++)
              irr[3 * (size_t)i + (size_t)c] = src[3 * k + (size_t)c];
            if (src[3 * k] > 0.0f) nfound++;
          }
          printf("   КАРТИНКА ИЗ РАЗВЁРТКИ: ячеек среза с радиансом %lld из %d (%.1f %%)\n",
                 (long long)nfound, S.n, 100.0 * (double)nfound / (double)(S.n ? S.n : 1));
        }
        /* А772/А775: ЭТАЛОН ПРОВЕРЯЕТСЯ САМ. Марш идёт шагом , и тонкий заслон
         * он может проскочить. Пересчёт вдвое мельче: если множество затенённых
         * почти не изменилось, эталон устойчив, и доли расхождения со свипом
         * говорят про свип. Если изменилось — все эти доли наполовину про эталон. */
        {
          float *irrf = malloc(3 * (size_t)S.n * sizeof *irrf);
          if (irrf == NULL) exit(1);
          double tf = now_s();
          front_direct(&S, &fr, &P, &AL, irrf, 0.25, 0, NULL, &m);
          tf = now_s() - tf;
          int64_t nd = 0, na = 0, nb = 0;
          for (int32_t i = 0; i < S.n; i++) {
            int a1 = irr[3 * (size_t)i] > 0.0f, b1 = irrf[3 * (size_t)i] > 0.0f;
            na += a1;
            nb += b1;
            nd += (a1 != b1);
          }
          printf("   ЭТАЛОН ПРИ ПОЛОВИННОМ ШАГЕ (%.1f мс): освещённых h/2 %lld, h/4 %lld, "
                 "РАЗОШЛИСЬ %lld (%.3f %% от среза)\n",
                 tf * 1e3, (long long)na, (long long)nb, (long long)nd,
                 100.0 * (double)nd / (double)(S.n ? S.n : 1));
          free(irrf);
        }

        /* ЭТАЛОН С ПОДЪЁМОМ (А784). Предикат тот же, значит множество затенённых
         * обязано СОВПАСТЬ побитово; падают только шаги и время. */
        {
          float *irrh = malloc(3 * (size_t)S.n * sizeof *irrh);
          if (irrh == NULL) exit(1);
          int64_t nsteph = 0;
          double th = now_s();
          /* АРБИТР — ПЛОСКИЙ марш: рабочий путь стал иерархическим, и сравнивать
           * его с самим собою значило бы печатать ложный ноль. */
          front_direct(&S, &fr, &P, &AL, irrh, 0.5, 0, &nsteph, &m);
          th = now_s() - th;
          int64_t nd2 = 0;
          for (int32_t i = 0; i < S.n; i++)
            if ((irr[3 * (size_t)i] > 0.0f) != (irrh[3 * (size_t)i] > 0.0f)) nd2++;
          printf("   АРБИТР (плоский марш): %.1f мс против %.1f мс рабочего с подъёмом (в %.2f "
                 "раза), шагов %lld (%.1f на луч); РАСХОЖДЕНИЕ %lld ячеек\n",
                 th * 1e3, t_dir * 1e3, t_dir / (th > 0.0 ? th : 1.0), (long long)nsteph,
                 (double)nsteph / (double)((int64_t)S.n * HZ_LIGHT_SAMPLES), (long long)nd2);
          free(irrh);
        }

        /* СВИП — то, что Ш5 обязан измерить
    ; луч выше остаётся ЭТАЛОНОМ (А763), и
         * сверка идёт ПОЯЧЕЕЧНО, а не по картинке. */
        irr2 = malloc(3 * (size_t)S.n * sizeof *irr2);
        if (irr2 == NULL) exit(1);
        double t_sw = 0.0, t_ga = 0.0;
        /* Ф5. (§514): камера дерева свипа ставится ЗДЕСЬ, непосредственно перед
         * прогоном. Прежде она ставилась в ветви `lit`, которая идёт ПОЗЖЕ, и
         * дерево строилось с камерой в нуле — то есть порог не действовал вовсе.
         * Замер это и показал: число листьев не менялось ни при каком пороге. */
        for (int a2 = 0; a2 < 3; a2++)
          g_sweepeye[a2] = LL.eye[a2];
        g_sweeppx = LL.pxrad;
        front_sweep(&S, &fr, &P, &AL, irr2, &m, &t_sw, &t_ga, sweepaxis, sweepfrac, sweepr01);
        /* ДОЛЯ ОТКРЫТОСТИ, А НЕ ОБЛУЧЁННОСТЬ (§423, А781). Приёмка задана на долю,
         * поэтому нужен знаменатель — облучённость БЕЗ всякого затенения. Считается
         * третьим проходом, дешёвым (теста заслона в нём нет вовсе). */
        float *irro = malloc(3 * (size_t)S.n * sizeof *irro);
        if (irro == NULL) exit(1);
        front_direct(&S, &fr, &P, &AL, irro, -1.0, 0, NULL, &m);
        {
          double *dv = malloc((size_t)S.n * sizeof *dv);
          if (dv == NULL) exit(1);
          int64_t nv2 = 0;
          for (int32_t i = 0; i < S.n; i++) {
            double den = (double)irro[3 * (size_t)i];
            if (!(den > 0.0)) continue;
            double a = (double)irr[3 * (size_t)i] / den, b = (double)irr2[3 * (size_t)i] / den;
            dv[nv2++] = fabs(a - b);
          }
          double q90 = 0.0, q99 = 0.0, qmax = 0.0;
          if (nv2 > 0) {
            qsort(dv, (size_t)nv2, sizeof *dv, cmp_d);
            q90 = dv[(nv2 * 9) / 10];
            q99 = dv[(nv2 * 99) / 100];
            qmax = dv[nv2 - 1];
          }
          {
            double ma = 0.0, mb = 0.0;
            int64_t nm = 0;
            for (int32_t i = 0; i < S.n; i++) {
              double den = (double)irro[3 * (size_t)i];
              if (!(den > 0.0)) continue;
              ma += (double)irr[3 * (size_t)i] / den;
              mb += (double)irr2[3 * (size_t)i] / den;
              nm++;
            }
            printf("   СРЕДНЯЯ ОТКРЫТОСТЬ: эталон %.4f, свип %.4f (по %lld ячейкам)\n",
                   ma / (double)(nm ? nm : 1), mb / (double)(nm ? nm : 1), (long long)nm);
          }
          printf("   ПРИЁМКА Ш5а2 (доля открытости против эталона, популяция — ячейки среза с "
                 "ненулевым знаменателем %lld): p90 %.4f, p99 %.4f, макс %.4f\n",
                 (long long)nv2, q90, q99, qmax);
          free(dv);
        }
        free(irro);
        {
          int64_t ndiff = 0, nlit_r = 0, nlit_s = 0, nsum_r = 0;
          double emax = 0.0, sum_r = 0.0;
          for (int32_t i = 0; i < S.n; i++) {
            double a = irr[3 * (size_t)i], b = irr2[3 * (size_t)i];
            if (a > 0.0) nlit_r++;
            if (b > 0.0) nlit_s++;
            if ((a > 0.0) != (b > 0.0)) ndiff++;
            double d = fabs(a - b);
            if (d > emax) emax = d;
            /* §506: РАЗОШЛИСЬ ПО МНОЖЕСТВУ СЧИТАЕТСЯ ПО «> 0», А У ДИФФУЗНОЙ СХЕМЫ
             * ХВОСТЫ НИКОГДА НЕ ДОХОДЯТ ДО НУЛЯ. Значит счётчик может мерить не
             * обтекание заслона, а экспоненциально малый шум. Здесь то же множество
             * считается по долям от СРЕДНЕЙ облучённости — порога нет, есть
             * зависимость, и она сама скажет, шум это или свет. */
            sum_r += a;
            nsum_r++;
          }
          double mean_r = nsum_r > 0 ? sum_r / (double)nsum_r : 0.0;
          int64_t ndf[3] = {0, 0, 0};
          static const double frac3[3] = {0.01, 0.05, 0.20};
          for (int32_t i = 0; i < S.n; i++) {
            double a = irr[3 * (size_t)i], b = irr2[3 * (size_t)i];
            for (int t2 = 0; t2 < 3; t2++) {
              double th = frac3[t2] * mean_r;
              if ((a > th) != (b > th)) ndf[t2]++;
            }
          }
          printf(
              "      §506 РАЗОШЛИСЬ ПО ПОРОГУ ОТ СРЕДНЕЙ: 1 %% — %lld (%.2f %%), 5 %% — %lld (%.2f "
              "%%), 20 %% — %lld (%.2f %%)\n",
              (long long)ndf[0], 100.0 * (double)ndf[0] / (double)(S.n ? S.n : 1),
              (long long)ndf[1], 100.0 * (double)ndf[1] / (double)(S.n ? S.n : 1),
              (long long)ndf[2], 100.0 * (double)ndf[2] / (double)(S.n ? S.n : 1));
          printf(
              "   СВИП: %.1f мс на %d образцов; сетка %d^3 = %lld ячеек, %.2f нс НА ЯЧЕЙКУ СЕТКИ; "
              "%.0f нс на ячейку СРЕЗА (это цена для БЮДЖЕТА, а не цена обработки); сбор %.1f мс; "
              "освещённых ячеек "
              "луч %lld, свип %lld, РАЗОШЛИСЬ %lld (%.2f %%), макс расхождение %.3e\n",
              t_sw * 1e3, HZ_LIGHT_SAMPLES, (int)1 << (lev - HZ_SWEEP_DROP),
              (long long)((int64_t)1 << (3 * (lev - HZ_SWEEP_DROP))),
              t_sw * 1e9 / (double)((int64_t)1 << (3 * (lev - HZ_SWEEP_DROP))) /
                  (double)HZ_LIGHT_SAMPLES,
              t_sw * 1e9 / (double)(S.n ? S.n : 1) / (double)HZ_LIGHT_SAMPLES, t_ga * 1e3,
              (long long)nlit_r, (long long)nlit_s, (long long)ndiff,
              100.0 * (double)ndiff / (double)(S.n ? S.n : 1), emax);
        }

        /* ---- ПЕЧЬ (Ш5б, §430): ЗАМКНУТАЯ ФОРМА ПРОТИВ ИТЕРАЦИИ ---- */
        /* В замкнутой полости с ПОСТОЯННЫМ альбедо и постоянной эмиссией угловые
         * коэффициенты каждой площадки суммируются в единицу, поэтому радиозность
         * удовлетворяет `B = E + ρ·B`, то есть `B = E/(1−ρ)` ТОЧНО. Ответ не зависит
         * ни от формы полости, ни от разбиения — потому это и приёмка: всякое
         * отклонение есть УТЕЧКА (или приток) энергии в моём переносе, а не
         * погрешность геометрии.
         * ЗДЕСЬ ПРОВЕРЯЕТСЯ МОЙ ГАТЕР, А НЕ АРИФМЕТИКА: угловые коэффициенты
         * считаются тем же кодом, что и отскок на сцене. */
        /* ---- ПЛОЩАДЬ ВЫДАННОЙ ПОВЕРХНОСТИ (§446): ПЕРВОЕ ДЕЙСТВИЕ, НАЗНАЧЕННОЕ
         * ЗАРАНЕЕ В §445. КЛЮЧ `area` — ЧТОБЫ ЗАМЕР ШЁЛ НА ЛЮБОЙ СЦЕНЕ (§450, П3а):
         * величина `1.5` снята только на коробке, а правило А793 требует гонять
         * эталон НА КАЖДОЙ НОВОЙ СЦЕНЕ. ---- */
        if (oven > 0.0 || area) {
          /* ИСТИНА БЕРЁТСЯ ИЗ САМОГО МЕША, А НЕ КОНСТАНТОЙ `6.0`. Тогда замер
           * остаётся верным при любом масштабе — и негативный контроль НК-1
           * (масштаб `2.0`) проверяет себя сам, а не сверяется с вписанным числом. */
          double atrue = 0.0;
          for (int32_t t3 = 0; t3 < m.nt; t3++) {
            const double *p0 = m.v + 3 * (size_t)m.f[3 * (size_t)t3];
            const double *p1 = m.v + 3 * (size_t)m.f[3 * (size_t)t3 + 1];
            const double *p2 = m.v + 3 * (size_t)m.f[3 * (size_t)t3 + 2];
            double nn[3];
            atrue += tri_area2(p0, p1, p2, nn);
          }
          /* ТАБЛИЦА РЁБЕР ДО ВСЯКОГО ОБХОДА: сколько пересечений в каждом СЛОЕ и
           * сколько из них сидят на КОНЦАХ ребра (`t` у нуля или у единицы). Обход и
           * вершины стоят ниже по цепочке, и судить по ним о причине — та же ошибка,
           * что А791 (искать вниз, когда надо вверх). */
          {
            int64_t *el = calloc(3 * (size_t)(fr.n + 1), sizeof *el);
            if (el == NULL) exit(1);
            int64_t nt0 = 0, nt1 = 0;
            for (int32_t i = 0; i < ht.n; i++) {
              int a = ht.e[i].axis;
              int32_t p = ht.e[i].p[a];
              if (p >= 0 && p <= fr.n) el[(size_t)a * (size_t)(fr.n + 1) + (size_t)p]++;
              /* Порог — не магический: это шаг ребра в долях, при котором точка
               * неотличима от конца в двойной точности на сетке в 2^lev ячеек. */
              if (ht.e[i].t < 1e-12) nt0++;
              if (ht.e[i].t > 1.0 - 1e-12) nt1++;
            }
            printf("   ТАБЛИЦА РЁБЕР: записей %d; `t` у НУЛЯ %lld, `t` у ЕДИНИЦЫ %lld "
                   "(поверхность на границе ячеек даёт и то и другое)\n",
                   ht.n, (long long)nt0, (long long)nt1);
            for (int a = 0; a < 3; a++) {
              int64_t tot = 0;
              for (int32_t p = 0; p <= fr.n; p++)
                tot += el[(size_t)a * (size_t)(fr.n + 1) + (size_t)p];
              if (tot == 0) continue;
              printf("   РЁБЕРА ПО СЛОЯМ ось %d (всего %lld):", a, (long long)tot);
              for (int32_t p = 0; p <= fr.n; p++) {
                int64_t c = el[(size_t)a * (size_t)(fr.n + 1) + (size_t)p];
                if (c * 100 < tot) continue;
                printf(" | слой %d: %lld", p, (long long)c);
              }
              printf("\n");
            }
            free(el);
          }

          /* ДВА ОГРАНИЧИТЕЛЯ ПОРОЗНЬ (§450, П3.2): ПОЛНАЯ ГЛУБИНА и СРЕЗ. Площадь
           * одной поверхности нельзя сравнивать с ячейками другой, а срез огрубляет
           * — значит числа разные, и печатать их надо врозь, а не одно за оба. */
          int32_t nbin = HZ_PLBIN * fr.n + 2;
          for (int pass = 0; pass < 2; pass++) {
            areacnt A;
            memset(&A, 0, sizeof A);
            A.nbin = nbin;
            A.pl_area = calloc(3 * (size_t)nbin, sizeof *A.pl_area);
            A.pl_cnt = calloc(3 * (size_t)nbin, sizeof *A.pl_cnt);
            A.ncell = fr.n + 1;
            A.cl_area = calloc(3 * (size_t)A.ncell, sizeof *A.cl_area);
            A.cl_cnt = calloc(3 * (size_t)A.ncell, sizeof *A.cl_cnt);
            if (A.pl_area == NULL || A.pl_cnt == NULL || A.cl_area == NULL || A.cl_cnt == NULL)
              exit(1);
            int32_t nskip = 0;
            int wrca = pass == 0 ? hz_dc_walk_stats(&T, NULL, NULL, area_emit, &A, &nskip)
                                 : hz_dc_walk_stats(&T, lod_stop, &LL, area_emit, &A, &nskip);
            double h2 = fr.h * fr.h;
            const char *tag = pass == 0 ? "ПОЛНАЯ ГЛУБИНА" : "СРЕЗ";
            printf("   ПЛОЩАДЬ ВЫДАННОЙ ПОВЕРХНОСТИ [%s] (§446): веер от v0 %.5f м², от v1 %.5f м² "
                   "(разность %.3e); ИСТИННАЯ по мешу %.5f м², отношение %.4f\n",
                   tag, A.fan0 * h2, A.fan1 * h2, fabs(A.fan0 - A.fan1) * h2, atrue,
                   A.fan0 * h2 / (atrue > 0.0 ? atrue : 1.0));
            printf(
                "   ПО ОСЯМ [%s] (главная ось нормали, БЕЗ знака): площадь %.5f / %.5f / %.5f м²; "
                "многоугольников %lld / %lld / %lld\n",
                tag, A.ax[0] * h2, A.ax[1] * h2, A.ax[2] * h2, (long long)A.px[0],
                (long long)A.px[1], (long long)A.px[2]);
            printf("   ОБХОД [%s]: многоугольников %lld, треугольников %lld, вырожденных %lld, "
                   "неплоскостность макс %.3e ячейки, ПРОПУЩЕНО полигонов %lld (код %d)\n",
                   tag, (long long)A.npoly, (long long)A.ntri, (long long)A.ndeg, A.flatmax,
                   (long long)nskip, wrca);
            /* ГИСТОГРАММА ПО ПЛОСКОСТЯМ (П2). Печатаются корзины, несущие не менее
             * сотой доли площади своей оси: иначе список утонет в хвосте из
             * единичных многоугольников на стыках стен. Отсечённая доля печатается,
             * чтобы «показано не всё» не читалось как «больше ничего нет». */
            for (int ax = 0; ax < 3; ax++) {
              if (!(A.ax[ax] > 0.0)) continue;
              printf("   ПЛОСКОСТИ [%s] ось %d (площадь оси %.5f м²):", tag, ax, A.ax[ax] * h2);
              double shown = 0.0;
              int nsh = 0;
              for (int32_t b = 0; b < nbin; b++) {
                double a = A.pl_area[(size_t)ax * (size_t)nbin + (size_t)b];
                if (!(a > 0.01 * A.ax[ax])) continue;
                printf(" | %.4f м: %.5f м² (%lld мн-ков)",
                       fr.org[ax] + (double)b / (double)HZ_PLBIN * fr.h, a * h2,
                       (long long)A.pl_cnt[(size_t)ax * (size_t)nbin + (size_t)b]);
                shown += a;
                nsh++;
              }
              printf(" || показано %d корзин, %.1f %% площади оси\n", nsh,
                     100.0 * shown / A.ax[ax]);
              printf("   СЛОИ ЯЧЕЕК [%s] ось %d:", tag, ax);
              int nsh2 = 0;
              double shown2 = 0.0;
              for (int32_t c = 0; c < A.ncell; c++) {
                double a = A.cl_area[(size_t)ax * (size_t)A.ncell + (size_t)c];
                if (!(a > 0.01 * A.ax[ax])) continue;
                printf(" | слой %d: %.5f м² (%lld мн-ков)", c, a * h2,
                       (long long)A.cl_cnt[(size_t)ax * (size_t)A.ncell + (size_t)c]);
                shown2 += a;
                nsh2++;
              }
              printf(" || показано %d слоёв, %.1f %% площади оси\n", nsh2,
                     100.0 * shown2 / A.ax[ax]);
            }
            if (pass == 0) {
              /* ЭТАЛОННЫЙ ОБХОД ТОЛЬКО НА ПОЛНОЙ ГЛУБИНЕ: на срезе отображение
               * «ребро -> полигон» у него отсутствует по построению (Г42). */
              areacnt R;
              memset(&R, 0, sizeof R);
              int wrcr = hz_dc_walk_ref(&T, &ht, NULL, NULL, area_emit, &R);
              printf("   ЭТАЛОННЫЙ ОБХОД (Г42, независимая реализация, код %d): %.5f м², "
                     "многоугольников %lld; РАСХОЖДЕНИЕ с рабочим %.3e отн.\n",
                     wrcr, R.fan0 * h2, (long long)R.npoly,
                     fabs(R.fan0 - A.fan0) / (A.fan0 > 0.0 ? A.fan0 : 1.0));
            }
            free(A.pl_area);
            free(A.pl_cnt);
            free(A.cl_area);
            free(A.cl_cnt);
          }
        }

        if (oven > 0.0) {
          double rho = oven, Le = 1.0;
          float *B = malloc(3 * (size_t)S.n * sizeof *B);
          float *Bn = malloc(3 * (size_t)S.n * sizeof *Bn);
          if (B == NULL || Bn == NULL) exit(1);
          for (int32_t i = 0; i < 3 * S.n; i++)
            B[i] = (float)Le;
          double exact = Le / (1.0 - rho);
          /* ДВЕНАДЦАТИ ОТСКОКОВ МАЛО ПРИ ВЫСОКОМ АЛЬБЕДО, И ЭТО АРИФМЕТИКА, А НЕ
           * догадка: невязка итерации есть `ρ^n`, то есть при `ρ = 0.7` и `n = 12`
           * она `1.4 %` — сравнима с тем систематическим смещением, которое печь и
           * должна измерять. Двадцать четыре дают `0.02 %` и разделяют их. */
          const int OVEN_ITERS = 24;
          for (int it = 1; it <= OVEN_ITERS; it++) {
            for (int32_t i = 0; i < 3 * S.n; i++)
              Bn[i] = (float)Le;
            for (int32_t j = 0; j < S.n; j++) {
              double pj[3], nj[3];
              hz_slice_vertex(&S, j, pj);
              for (int k = 0; k < 3; k++)
                pj[k] = fr.org[k] + pj[k] * fr.h;
              hz_slice_normal(&S, j, nj);
              double cs = fr.h * (double)((int32_t)1 << (lev - (int)S.c[j].lvl));
              double aj = cs * cs;
              for (int32_t i = 0; i < S.n; i++) {
                if (i == j) continue;
                double pi[3], ni[3], w[3], r2 = 0.0;
                hz_slice_vertex(&S, i, pi);
                for (int k = 0; k < 3; k++)
                  pi[k] = fr.org[k] + pi[k] * fr.h;
                hz_slice_normal(&S, i, ni);
                for (int k = 0; k < 3; k++) {
                  w[k] = pj[k] - pi[k];
                  r2 += w[k] * w[k];
                }
                if (!(r2 > 0.0)) continue;
                double r = sqrt(r2);
                double ci = (w[0] * ni[0] + w[1] * ni[1] + w[2] * ni[2]) / r;
                double cj = -(w[0] * nj[0] + w[1] * nj[1] + w[2] * nj[2]) / r;
                if (!(ci > 0.0) || !(cj > 0.0)) continue;
                double ff = ci * cj * aj / (3.14159265358979323846 * r2);
                for (int k = 0; k < 3; k++)
                  Bn[3 * (size_t)i + (size_t)k] +=
                      (float)(rho * (double)B[3 * (size_t)j + (size_t)k] * ff);
              }
            }
            /* ДИАГНОЗ §435, ПУНКТ (а) и (б): сколько пар прошло оба `cos > 0` и
             * чему равна сумма угловых коэффициентов ОДНОЙ площадки. В замкнутой
             * полости вторая обязана быть `1`; отклонение и есть мера того,
             * насколько гатер теряет энергию. */
            if (it == 1) {
              /* СУММА ПО ВСЕМ ПЛОЩАДКАМ, А НЕ ПО ОДНОЙ (А790). Площадка `0` —
               * угловая, и она видит меньше типичной; одно число с неё мерой
               * потери гатера не является. Здесь считается РАСПРЕДЕЛЕНИЕ. */
              {
                double *fs = malloc((size_t)S.n * sizeof *fs);
                if (fs == NULL) exit(1);
                for (int32_t jj = 0; jj < S.n; jj++) {
                  double pa[3], na[3];
                  hz_slice_vertex(&S, jj, pa);
                  for (int k = 0; k < 3; k++)
                    pa[k] = fr.org[k] + pa[k] * fr.h;
                  hz_slice_normal(&S, jj, na);
                  double acc2 = 0.0;
                  for (int32_t ii = 0; ii < S.n; ii++) {
                    if (ii == jj) continue;
                    double pb[3], nb[3], ww[3], rr2 = 0.0;
                    hz_slice_vertex(&S, ii, pb);
                    for (int k = 0; k < 3; k++)
                      pb[k] = fr.org[k] + pb[k] * fr.h;
                    hz_slice_normal(&S, ii, nb);
                    for (int k = 0; k < 3; k++) {
                      ww[k] = pb[k] - pa[k];
                      rr2 += ww[k] * ww[k];
                    }
                    if (!(rr2 > 0.0)) continue;
                    double rr = sqrt(rr2);
                    double caa = (ww[0] * na[0] + ww[1] * na[1] + ww[2] * na[2]) / rr;
                    double cbb = -(ww[0] * nb[0] + ww[1] * nb[1] + ww[2] * nb[2]) / rr;
                    if (!(caa > 0.0) || !(cbb > 0.0)) continue;
                    double csb = fr.h * (double)((int32_t)1 << (lev - (int)S.c[ii].lvl));
                    acc2 += caa * cbb * csb * csb / (3.14159265358979323846 * rr2);
                  }
                  fs[jj] = acc2;
                }
                qsort(fs, (size_t)S.n, sizeof *fs, cmp_d);
                double mean = 0.0;
                for (int32_t jj = 0; jj < S.n; jj++)
                  mean += fs[jj];
                mean /= (double)(S.n ? S.n : 1);
                printf(
                    "   СУММА УГЛОВЫХ КОЭФФИЦИЕНТОВ ПО ВСЕМ %d ПЛОЩАДКАМ: среднее %.4f, p10 %.4f, "
                    "p50 %.4f, p90 %.4f (обязана быть 1)\n",
                    S.n, mean, fs[S.n / 10], fs[S.n / 2], fs[(S.n * 9) / 10]);
                free(fs);
              }
              int64_t npair = 0;
              double ffsum = 0.0;
              int32_t j0 = 0;
              double pj0[3], nj0[3];
              hz_slice_vertex(&S, j0, pj0);
              for (int k = 0; k < 3; k++)
                pj0[k] = fr.org[k] + pj0[k] * fr.h;
              hz_slice_normal(&S, j0, nj0);
              for (int32_t i = 0; i < S.n; i++) {
                if (i == j0) continue;
                double pi[3], ni[3], w[3], r2 = 0.0;
                hz_slice_vertex(&S, i, pi);
                for (int k = 0; k < 3; k++)
                  pi[k] = fr.org[k] + pi[k] * fr.h;
                hz_slice_normal(&S, i, ni);
                for (int k = 0; k < 3; k++) {
                  w[k] = pi[k] - pj0[k];
                  r2 += w[k] * w[k];
                }
                if (!(r2 > 0.0)) continue;
                double r = sqrt(r2);
                double cj = (w[0] * nj0[0] + w[1] * nj0[1] + w[2] * nj0[2]) / r;
                double ci = -(w[0] * ni[0] + w[1] * ni[1] + w[2] * ni[2]) / r;
                if (!(ci > 0.0) || !(cj > 0.0)) continue;
                double cs2 = fr.h * (double)((int32_t)1 << (lev - (int)S.c[i].lvl));
                npair++;
                ffsum += ci * cj * cs2 * cs2 / (3.14159265358979323846 * r2);
              }
              {
                double atot = 0.0;
                for (int32_t i = 0; i < S.n; i++) {
                  double cs3 = fr.h * (double)((int32_t)1 << (lev - (int)S.c[i].lvl));
                  atot += cs3 * cs3;
                }
                {
                  /* СКОЛЬКО ЯЧЕЕК НА ГРАНЬ. У единичной коробки при lev=5 стена
                   * занимает 30x30 = 900 клеток; если ячеек среза меньше, значит
                   * вершина выдана не в каждой клетке стены — это и есть недосчёт
                   * площади. Группировка по ГЛАВНОЙ оси нормали: шесть граней. */
                  int64_t hg[6] = {0, 0, 0, 0, 0, 0};
                  for (int32_t i = 0; i < S.n; i++) {
                    double nn2[3];
                    hz_slice_normal(&S, i, nn2);
                    int ax = 0;
                    for (int k = 1; k < 3; k++)
                      if (fabs(nn2[k]) > fabs(nn2[ax])) ax = k;
                    hg[2 * ax + (nn2[ax] > 0.0 ? 1 : 0)]++;
                  }
                  printf(
                      "   ЯЧЕЕК НА ГРАНЬ (по главной оси нормали): %lld %lld %lld %lld %lld %lld "
                      "против 900 клеток стены\n",
                      (long long)hg[0], (long long)hg[1], (long long)hg[2], (long long)hg[3],
                      (long long)hg[4], (long long)hg[5]);
                }
                printf("   ПЛОЩАДЬ СРЕЗА: %.5f м² против истинной 6.00000 м² у единичной коробки "
                       "(отношение %.4f)\n",
                       atot, atot / 6.0);
              }
              printf("   ДИАГНОЗ ПЕЧИ: у площадки 0 видимых партнёров %lld из %d; СУММА УГЛОВЫХ "
                     "КОЭФФИЦИЕНТОВ %.5f (в замкнутой полости обязана быть 1)\n",
                     (long long)npair, S.n - 1, ffsum);
            }
            double sum = 0.0;
            for (int32_t i = 0; i < S.n; i++)
              sum += (double)Bn[3 * (size_t)i];

            double mean = sum / (double)(S.n ? S.n : 1);
            printf("   ПЕЧЬ ρ=%.2f, отскок %2d: средняя B = %.5f против замкнутой формы %.5f "
                   "(отклонение %.2f %%)\n",
                   rho, it, mean, exact, 100.0 * (mean - exact) / exact);
            for (int32_t i = 0; i < 3 * S.n; i++)
              B[i] = Bn[i];
          }
          free(B);
          free(Bn);
        }

        /* ---- ДВЕ ПЛАСТИНЫ: `E_ind = ρ·E_dir·F` (Ш5б, §430 П5б.2; §450 П6) ---- */
        /* ПЕЧЬ ЛОВИТ СОХРАНЕНИЕ ЭНЕРГИИ, ЭТОТ СТЕНД — ЕЁ РАСПРЕДЕЛЕНИЕ. В замкнутой
         * полости сумма угловых коэффициентов равна единице при ЛЮБОМ разумном ядре,
         * лишь бы оно было симметрично и нормировано; отдельные коэффициенты она не
         * проверяет. Два соосных квадрата проверяют именно их: `F` известен в
         * замкнутой форме (каталог Хауэлла C-11), и косинусы с `1/r²` входят в него
         * порознь.
         * ЧЕГО ЭТОТ СТЕНД НЕ ПРОВЕРЯЕТ, И ЭТО СКАЗАНО ЗДЕСЬ, А НЕ В ДОКЛАДЕ: ЗАСЛОНЫ.
         * Гатер незаслонённый (приближение (1) §432), значит требование А756 —
         * «печь не ловит тени» — этим стендом ТОЖЕ не закрывается. Обе половины
         * приёмки §430 меряют неэкранированный перенос, и тени остаются
         * неизмеренными вовсе. */
        if (plates > 0.0) {
          double rho = plates;
          /* Габарит пластин и зазор берутся ИЗ СЦЕНЫ, а не вписываются: стенд обязан
           * оставаться верным, если пластины подвинут. */
          double side_a = hi[0] - lo[0], gap = hi[1] - lo[1], mid = 0.5 * (lo[1] + hi[1]);
          double X = side_a / gap;
          double X2 = X * X, s = sqrt(1.0 + X2);
          double F = (2.0 / (3.14159265358979323846 * X2)) *
                     (0.5 * log((1.0 + X2) * (1.0 + X2) / (1.0 + 2.0 * X2)) +
                      2.0 * X * s * atan(X / s) - 2.0 * X * atan(X));
          double acc = 0.0;
          int64_t nup = 0, nlo = 0;
          for (int32_t i = 0; i < S.n; i++) {
            double pi[3], ni[3];
            hz_slice_vertex(&S, i, pi);
            for (int k = 0; k < 3; k++)
              pi[k] = fr.org[k] + pi[k] * fr.h;
            if (pi[1] < mid) {
              nlo++;
              continue;
            }
            nup++;
            hz_slice_normal(&S, i, ni);
            double sum = 0.0;
            for (int32_t j = 0; j < S.n; j++) {
              double pj[3], nj[3], w[3], r2 = 0.0;
              hz_slice_vertex(&S, j, pj);
              for (int k = 0; k < 3; k++)
                pj[k] = fr.org[k] + pj[k] * fr.h;
              if (pj[1] >= mid) continue; /* излучает только НИЖНЯЯ пластина */
              hz_slice_normal(&S, j, nj);
              for (int k = 0; k < 3; k++) {
                w[k] = pj[k] - pi[k];
                r2 += w[k] * w[k];
              }
              if (!(r2 > 0.0)) continue;
              double r = sqrt(r2);
              double ci = (w[0] * ni[0] + w[1] * ni[1] + w[2] * ni[2]) / r;
              double cj = -(w[0] * nj[0] + w[1] * nj[1] + w[2] * nj[2]) / r;
              if (!(ci > 0.0) || !(cj > 0.0)) continue;
              double cs = fr.h * (double)((int32_t)1 << (lev - (int)S.c[j].lvl));
              sum += ci * cj * cs * cs / (3.14159265358979323846 * r2);
            }
            acc += sum;
          }
          double mean = acc / (double)(nup ? nup : 1);
          printf("   ПЛАСТИНЫ: сторона %.4f м, зазор %.4f м, X = %.3f; ЯЧЕЕК верх %lld, низ %lld\n",
                 side_a, gap, X, (long long)nup, (long long)nlo);
          printf("   `E_ind = ρ·E_dir·F` при ρ=%.2f: замерено %.5f, замкнутая форма %.5f "
                 "(F = %.5f), ОТКЛОНЕНИЕ %.2f %%\n",
                 rho, rho * mean, rho * F, F, 100.0 * (mean - F) / F);
        }

        /* ---- ОДИН ОТСКОК (Ш5б, §430) ---- */
        /* ПРИБЛИЖЕНИЯ НАЗЫВАЮТСЯ ЗДЕСЬ, А НЕ В ДОКЛАДЕ ЗАДНИМ ЧИСЛОМ.
         *   (1) ВИДИМОСТИ МЕЖДУ ЯЧЕЙКАМИ НЕТ: перенос идёт по незаслонённому
         *       угловому коэффициенту. Значит свет проходит сквозь стены, и на
         *       сцене с комнатами это ВИДНО. Взято сознательно: с видимостью цена
         *       умножается на марш (~4.8 шага), а замер физики от заслонов не
         *       зависит — `E_ind = ρ·E_dir` проверяется на ПЛОСКОЙ стене, где
         *       заслонов нет вовсе.
         *   (2) ИЗЛУЧАТЕЛИ ПРОРЕЖЕНЫ шагом `stride`: берётся каждый `stride`-й, а
         *       вклад умножается на `stride`. Это несмещённая оценка суммы, но с
         *       разбросом; разброс НЕ ИЗМЕРЕН и в приёмку не входит.
         *   (3) ПЛОЩАДЬ ЯЧЕЙКИ взята как площадь её грани `(h·2^(lev-lvl))²` —
         *       поверхность внутри ячейки наклонена и её площадь больше; это
         *       систематическая недооценка, названная и не исправленная.
         * ЦЕНА ОЖИДАЕТСЯ ПЛОХОЙ И ПРЕДСКАЗАНА ДО ПРОГОНА (§430, П5б.3): это замер,
         * обосновывающий свип, а не попытка уложиться в бюджет. */
        float *ind = calloc(3 * (size_t)S.n, sizeof *ind);
        float *indsw = NULL;
        if (ind == NULL) exit(1);
        /* ---- Ф9' (§536): САМАЯ УЗКАЯ ДОЛЯ, КОТОРУЮ НЕСЁТ ДАННОЕ ND ---- */
        /* ЗАМЕР БЕЗ ПЕРЕНОСА, И ЭТО СОЗНАТЕЛЬНО. Граница «до какой узости ординаты
         * несут отражение» есть свойство НАБОРА ОРДИНАТ и геометрии, а не света:
         * ни источника, ни камеры в ней нет. Значит мерить её надо отдельно от
         * транспорта, иначе три разные ошибки сложатся в одно число (А934).
         *
         * ДВЕ ПОЛОВИНЫ ГРАНИЦЫ МЕРЯЮТСЯ ПОРОЗНЬ (А947), потому что К56 говорит
         * ровно о том, что они РАЗНЫЕ:
         *   A — квадратурное альбедо доли `Σ_e w_e f_r cos θ_e`. Это ЭНЕРГИЯ.
         *       Сравнивается не с `ρ_s` (нормировка Фонга точна лишь при нормальном
         *       падении и спутала бы свою погрешность с квадратурной), а с ТЕМ ЖЕ
         *       интегралом на заведомо избыточном наборе.
         *   УГОЛ ПИКА — на какую ординату легла вершина доли против истинного
         *       зеркального направления. Это ОБРАЗ, и просил пользователь именно
         *       его. `A` к нему слепа: интеграл сходится и тогда, когда пик уехал.
         *
         * НОРМАЛИ БЕРУТСЯ ИЗ СЦЕНЫ, А НЕ ПРИДУМЫВАЮТСЯ: выборка по срезу с шагом,
         * число печатается. */
        if (g_lobetest) {
          static const double SS[] = {0.0, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 20.0, 30.0, 50.0, 200.0};
          static const int NMU[] = {1, 2, 2, 4}, NPHI[] = {1, 2, 4, 4};
          tr3_dirs RF;
          /* ЭТАЛОННЫЙ НАБОР — тот же механизм, вчетверо гуще самого густого из
           * испытуемых по каждой оси. Своей аналитики здесь нет намеренно: она
           * внесла бы вторую формулу, и замер мерил бы разницу формул. */
          if (tr3_dirs_product(&RF, 16, 16) != 0) exit(1);
          int32_t nstep = S.n / 256 > 0 ? S.n / 256 : 1;
          int64_t nsmp = 0;
          for (int32_t i = 0; i < S.n; i += nstep)
            nsmp++;
          printf("   Ф9' ГРАНИЦА УЗОСТИ: нормалей из среза %lld (шаг %d из %d), эталон ND %d\n",
                 (long long)nsmp, nstep, S.n, RF.n);
          printf("        ND  полуугол |  s  | ОТКЛОНЕНИЕ A, %% (p50/p90/max) | УГОЛ ПИКА, град "
                 "(p50/p90/max)\n");
          for (int ci = 0; ci < 4; ci++) {
            tr3_dirs DT;
            if (tr3_dirs_product(&DT, NMU[ci], NPHI[ci]) != 0) exit(1);
            double thnd = acos(1.0 - 2.0 / (double)DT.n) * 180.0 / 3.14159265358979323846;
            for (size_t si = 0; si < sizeof SS / sizeof SS[0]; si++) {
              double ns2 = SS[si];
              double *ea = malloc((size_t)(nsmp * DT.n) * sizeof *ea);
              double *pa = malloc((size_t)(nsmp * DT.n) * sizeof *pa);
              if (ea == NULL || pa == NULL) exit(1);
              int64_t ne = 0;
              for (int32_t i = 0; i < S.n; i += nstep) {
                double nn3[3];
                hz_slice_normal(&S, i, nn3);
                for (int d = 0; d < DT.n; d++) {
                  double od[3] = {DT.ox[d], DT.oy[d], DT.oz[d]};
                  double dn2 = od[0] * nn3[0] + od[1] * nn3[1] + od[2] * nn3[2];
                  if (!(dn2 < 0.0)) continue; /* не падает на эту сторону */
                  /* Испытуемый набор: альбедо доли и ординату пика. */
                  double at = 0.0, pk = -1.0;
                  int be = -1;
                  for (int e = 0; e < DT.n; e++) {
                    double oe[3] = {DT.ox[e], DT.oy[e], DT.oz[e]};
                    double en2 = oe[0] * nn3[0] + oe[1] * nn3[1] + oe[2] * nn3[2];
                    if (!(en2 > 0.0)) continue;
                    double ca = oe[0] * od[0] + oe[1] * od[1] + oe[2] * od[2] - 2.0 * dn2 * en2;
                    double fv = lobe(ca, ns2) * DT.w[e] * en2;
                    at += fv;
                    if (fv > pk) {
                      pk = fv;
                      be = e;
                    }
                  }
                  /* Эталон: тот же интеграл на избыточном наборе. */
                  double ar = 0.0;
                  for (int e = 0; e < RF.n; e++) {
                    double oe[3] = {RF.ox[e], RF.oy[e], RF.oz[e]};
                    double en2 = oe[0] * nn3[0] + oe[1] * nn3[1] + oe[2] * nn3[2];
                    if (!(en2 > 0.0)) continue;
                    double ca = oe[0] * od[0] + oe[1] * od[1] + oe[2] * od[2] - 2.0 * dn2 * en2;
                    ar += lobe(ca, ns2) * RF.w[e] * en2;
                  }
                  if (!(ar > 0.0)) continue;
                  ea[ne] = 100.0 * fabs(at - ar) / ar;
                  /* Угол между ординатой пика и ИСТИННЫМ зеркальным направлением. */
                  double mr[3] = {od[0] - 2.0 * dn2 * nn3[0], od[1] - 2.0 * dn2 * nn3[1],
                                  od[2] - 2.0 * dn2 * nn3[2]};
                  if (be >= 0) {
                    double cm = DT.ox[be] * mr[0] + DT.oy[be] * mr[1] + DT.oz[be] * mr[2];
                    if (cm > 1.0) cm = 1.0;
                    if (cm < -1.0) cm = -1.0;
                    pa[ne] = acos(cm) * 180.0 / 3.14159265358979323846;
                  } else
                    pa[ne] = 180.0;
                  ne++;
                }
              }
              if (ne > 0) {
                qsort(ea, (size_t)ne, sizeof *ea, cmp_d);
                qsort(pa, (size_t)ne, sizeof *pa, cmp_d);
                printf("       %4d  %6.1f°  |%5.0f| %8.2f %8.2f %8.2f      | %8.2f %8.2f %8.2f\n",
                       DT.n, thnd, ns2, ea[ne / 2], ea[(ne * 9) / 10], ea[ne - 1], pa[ne / 2],
                       pa[(ne * 9) / 10], pa[ne - 1]);
              }
              free(ea);
              free(pa);
            }
            tr3_dirs_free(&DT);
          }
          tr3_dirs_free(&RF);
        }

        /* ---- Ф8' (§532): ОТСКОК НАПРАВЛЕННЫМ СВИПОМ ПО ДЕРЕВУ ---- */
        if (g_dsweep) {
          tr3_dirs DR;
          if (tr3_dirs_product(&DR, g_dnmu, g_dnphi) != 0) exit(1);
          stree TD;
          double tb3 = now_s();
          double eyes2[3];
          for (int a = 0; a < 3; a++)
            eyes2[a] = LL.eye[a] / (double)((int32_t)1 << HZ_SWEEP_DROP);
          int32_t gn2 = (int32_t)1 << (fr.lev - HZ_SWEEP_DROP);
          stree_build(&TD, &P, fr.lev, HZ_SWEEP_DROP, gn2, eyes2, g_sweeppx, g_sweepthr);
          stree_links(&TD, gn2);
          dfield D;
          memset(&D, 0, sizeof D);
          D.nsl = S.n;
          D.blkopen = g_dblkopen;
          D.fone = g_ffull;
          D.fzero = g_fzero;
          D.fnotrans = g_fnotrans;
          D.fnoclamp = g_fnoclamp;
          double emact = 0.0, emcut = 0.0;
          D.emitact = &emact;
          D.emitcut = &emcut;
          D.frawv = malloc((size_t)HZ_FSTAT_CAP * sizeof *D.frawv);
          if (D.frawv == NULL) exit(1);
          D.schar = g_schar;
          D.scharax = g_scharax;
          D.onenb = g_onenb;
          D.cw = fr.h * (double)((int32_t)1 << HZ_SWEEP_DROP);
          D.Ap = calloc((size_t)TD.n, sizeof *D.Ap);
          D.fstat = malloc((size_t)HZ_FSTAT_CAP * sizeof *D.fstat);
          if (D.Ap == NULL || D.fstat == NULL) exit(1);
          D.Ld = calloc(3 * (size_t)TD.n, sizeof *D.Ld);
          D.Bs = calloc(3 * (size_t)TD.n, sizeof *D.Bs);
          D.Bn = calloc(3 * (size_t)TD.n, sizeof *D.Bn);
          D.Bsp = calloc(3 * (size_t)TD.n, sizeof *D.Bsp);
          D.Bwi = calloc(3 * (size_t)TD.n, sizeof *D.Bwi);
          D.Bns = calloc((size_t)TD.n, sizeof *D.Bns);
          D.srf = calloc((size_t)TD.n, 1);
          D.vis = calloc((size_t)TD.n, 1);
          D.cstart = calloc((size_t)TD.n + 1, sizeof *D.cstart);
          D.clist = malloc((size_t)(S.n > 0 ? S.n : 1) * sizeof *D.clist);
          D.Eind = calloc(3 * (size_t)S.n, sizeof *D.Eind);
          double *snx = malloc((size_t)(S.n > 0 ? S.n : 1) * sizeof *snx);
          double *sny = malloc((size_t)(S.n > 0 ? S.n : 1) * sizeof *sny);
          double *snz = malloc((size_t)(S.n > 0 ? S.n : 1) * sizeof *snz);
          double *sar = malloc((size_t)(S.n > 0 ? S.n : 1) * sizeof *sar);
          int32_t *slf = malloc((size_t)(S.n > 0 ? S.n : 1) * sizeof *slf);
          double *wchk = calloc((size_t)(S.n > 0 ? S.n : 1), sizeof *wchk);
          if (D.Ld == NULL || D.Bs == NULL || D.Bn == NULL || D.Bsp == NULL || D.Bwi == NULL ||
              D.Bns == NULL || D.srf == NULL || D.vis == NULL || D.cstart == NULL ||
              D.clist == NULL || D.Eind == NULL || snx == NULL || sny == NULL || snz == NULL ||
              sar == NULL || slf == NULL || wchk == NULL)
            exit(1);
          D.snx = snx;
          D.sny = sny;
          D.snz = snz;
          /* ОТОБРАЖЕНИЕ «ЯЧЕЙКА СРЕЗА -> ЛИСТ» тем же правилом, что у `sweep_vis`:
           * иначе перенос и сбор читали бы разные ячейки. */
          double cwid = fr.h * (double)((int32_t)1 << HZ_SWEEP_DROP);
          int64_t nmap = 0;
          for (int32_t i = 0; i < S.n; i++) {
            double p2[3], n2[3];
            hz_slice_vertex(&S, i, p2);
            for (int k = 0; k < 3; k++)
              p2[k] = fr.org[k] + p2[k] * fr.h;
            hz_slice_normal(&S, i, n2);
            snx[i] = n2[0];
            sny[i] = n2[1];
            snz[i] = n2[2];
            double cs4 = fr.h * (double)((int32_t)1 << (lev - (int)S.c[i].lvl));
            sar[i] = cs4 * cs4;
            int32_t cc2[3];
            int ok3 = 1;
            for (int k = 0; k < 3; k++) {
              double f3 = floor((p2[k] - fr.org[k]) / cwid);
              if (!(f3 >= 0.0) || !(f3 < (double)gn2)) ok3 = 0;
              cc2[k] = ok3 ? (int32_t)f3 : 0;
            }
            slf[i] = ok3 ? TD.idx[hz_occ_index(gn2, cc2[0], cc2[1], cc2[2])] : -1;
            if (slf[i] >= 0) {
              nmap++;
              D.cstart[slf[i] + 1]++;
            }
          }
          for (int32_t i = 0; i < TD.n; i++)
            D.cstart[i + 1] += D.cstart[i];
          {
            int32_t *cur = malloc((size_t)TD.n * sizeof *cur);
            if (cur == NULL) exit(1);
            memcpy(cur, D.cstart, (size_t)TD.n * sizeof *cur);
            for (int32_t i = 0; i < S.n; i++)
              if (slf[i] >= 0) D.clist[cur[slf[i]]++] = i;
            free(cur);
          }
          /* РАДИОСИТИ ЛИСТА — СРЕДНЕЕ ПО ПЛОЩАДИ, А НЕ СУММА (А936): радианс есть
           * величина УДЕЛЬНАЯ, и две ячейки среза в одном листе не светят вдвое
           * ярче. Делится на `π`, потому что диффузная поверхность с радиосити `B`
           * имеет радианс `B/π`. */
          {
            double *aw2 = calloc((size_t)TD.n, sizeof *aw2);
            if (aw2 == NULL) exit(1);
            for (int32_t i = 0; i < S.n; i++) {
              int32_t l2 = slf[i];
              if (l2 < 0) continue;
              aw2[l2] += sar[i];
              /* УЗОСТЬ И ЗЕРКАЛЬНАЯ ДОЛЯ ИЗ МАТЕРИАЛА (Ф9', §536). `Ns` и `Ks` уже
               * разбираются `scene_obj.c` и до сих пор не читались никем; ключи
               * `dns=`/`dks=` перебивают их глобально — для СВИПА по узости без
               * правки сцены.
               * `Ns = 0` ЧИТАЕТСЯ КАК «ДОЛИ НЕТ», а не как `s = 0` (А949): в `.mtl`
               * ноль пишут и диффузным, и по умолчанию, и принять молчание формата
               * за значение — та же ошибка, что А917.
               * `E_dir` ВОССТАНАВЛИВАЕТСЯ ДЕЛЕНИЕМ на диффузное альбедо, потому что
               * `irr` хранит уже `rho_d·E`. При `rho_d = 0` восстановить нечего, и
               * узкая часть там просто не заводится — оговорка, а не молчание. */
              double nsi =
                  g_dns > 0.0 ? g_dns
                              : (m.mtl != NULL && S.c[i].mat < m.nmtl ? m.mtl[S.c[i].mat].ns : 0.0);
              double wi2[3] = {0.0, 0.0, 0.0}, wl = 0.0;
              for (int k = 0; k < 3; k++) {
                double pw2[3];
                hz_slice_vertex(&S, i, pw2);
                wi2[k] = (fr.org[k] + pw2[k] * fr.h) - AL.c[k];
                wl += wi2[k] * wi2[k];
              }
              wl = sqrt(wl);
              for (int k = 0; k < 3; k++) {
                D.Bs[3 * (size_t)l2 + (size_t)k] +=
                    (float)(sar[i] * (double)irr[3 * (size_t)i + (size_t)k]);
                D.Bn[3 * (size_t)l2 + (size_t)k] +=
                    (float)(sar[i] * (k == 0 ? snx[i] : (k == 1 ? sny[i] : snz[i])));
                if (nsi > 0.0 && wl > 0.0) {
                  double rd = alb(&m, S.c[i].mat, k);
                  double rs =
                      g_dks >= 0.0
                          ? g_dks
                          : (m.mtl != NULL && S.c[i].mat < m.nmtl ? m.mtl[S.c[i].mat].ks3[k] : 0.0);
                  double ed = rd > 1e-6 ? (double)irr[3 * (size_t)i + (size_t)k] / rd : 0.0;
                  D.Bsp[3 * (size_t)l2 + (size_t)k] += (float)(sar[i] * rs * ed);
                  D.Bwi[3 * (size_t)l2 + (size_t)k] += (float)(sar[i] * wi2[k] / wl);
                }
              }
              if (nsi > 0.0) D.Bns[l2] += (float)(sar[i] * nsi);
            }
            int64_t nsrf = 0, nblk = 0, nglos = 0, nfold = 0, nmulti = 0;
            double *foldv = malloc((size_t)HZ_FSTAT_CAP * sizeof *foldv);
            if (foldv == NULL) exit(1);
            for (int32_t l2 = 0; l2 < TD.n; l2++) {
              if (TD.nd[l2].child0 >= 0) continue;
              if (aw2[l2] > 0.0) {
                double nl = 0.0;
                for (int k = 0; k < 3; k++) {
                  D.Bs[3 * (size_t)l2 + (size_t)k] =
                      (float)((double)D.Bs[3 * (size_t)l2 + (size_t)k] / aw2[l2] /
                              3.14159265358979323846);
                  nl += (double)D.Bn[3 * (size_t)l2 + (size_t)k] *
                        (double)D.Bn[3 * (size_t)l2 + (size_t)k];
                }
                nl = sqrt(nl);
                if (nl > 0.0)
                  for (int k = 0; k < 3; k++)
                    D.Bn[3 * (size_t)l2 + (size_t)k] =
                        (float)((double)D.Bn[3 * (size_t)l2 + (size_t)k] / nl);
                /* Узкая часть: амплитуда — среднее по площади (А936, тот же довод,
                 * что у диффузной); направление прихода нормируется; показатель —
                 * среднее по площади. Делить на `π` здесь НЕ надо: нормировка
                 * `(s+2)/2π` сидит в самой доле. */
                {
                  double wl2 = 0.0;
                  for (int k = 0; k < 3; k++) {
                    D.Bsp[3 * (size_t)l2 + (size_t)k] =
                        (float)((double)D.Bsp[3 * (size_t)l2 + (size_t)k] / aw2[l2]);
                    wl2 += (double)D.Bwi[3 * (size_t)l2 + (size_t)k] *
                           (double)D.Bwi[3 * (size_t)l2 + (size_t)k];
                  }
                  wl2 = sqrt(wl2);
                  if (wl2 > 0.0)
                    for (int k = 0; k < 3; k++)
                      D.Bwi[3 * (size_t)l2 + (size_t)k] =
                          (float)((double)D.Bwi[3 * (size_t)l2 + (size_t)k] / wl2);
                  D.Bns[l2] = (float)((double)D.Bns[l2] / aw2[l2]);
                  if (D.Bns[l2] > 0.0f) nglos++;
                }
                D.Ap[l2] = (float)aw2[l2];
                /* А971: СКЛАДЧАТОСТЬ ЛИСТА `1 − |Σ A_i n_i| / Σ A_i`. У плоской
                 * площадки ноль, у угла — заметно больше нуля, и тогда усреднённая
                 * нормаль не значит направления поверхности (§518/§519). Печатается,
                 * чтобы оговорка была числом, а не словом. */
                if (nfold < (int64_t)HZ_FSTAT_CAP) foldv[nfold++] = 1.0 - nl / aw2[l2];
                D.srf[l2] = 1;
                nsrf++;
              } else {
                /* А935: ЛЮБОЙ занятый лист заслоняет, даже если среза в нём нет —
                 * иначе там, где срез огрублён, свет пошёл бы сквозь стену. Такой
                 * лист есть ЧЁРНАЯ стена: гасит, но не светит. */
                int32_t lof[3] = {TD.nd[l2].lo[0] << HZ_SWEEP_DROP,
                                  TD.nd[l2].lo[1] << HZ_SWEEP_DROP,
                                  TD.nd[l2].lo[2] << HZ_SWEEP_DROP};
                if (u_occ(&P, lof, TD.nd[l2].size << HZ_SWEEP_DROP)) {
                  D.srf[l2] = 2;
                  nblk++;
                }
              }
            }
            /* А970: сколько листьев несут ДВЕ и более ячейки среза — там проекции
             * площадок могут перекрываться, и `f` завышена. */
            for (int32_t l2 = 0; l2 < TD.n; l2++)
              if (D.srf[l2] == 1 && D.cstart[l2 + 1] - D.cstart[l2] >= 2) nmulti++;
            if (nfold > 0) {
              qsort(foldv, (size_t)nfold, sizeof *foldv, cmp_d);
              printf(
                  "      А971 СКЛАДЧАТОСТЬ ЛИСТА (1 − |ΣAn|/ΣA): p50 %.4f, p90 %.4f, max %.4f по "
                  "%lld листьям; А970 листьев с 2+ ячейками среза %lld (%.1f %%)\n",
                  foldv[nfold / 2], foldv[(nfold * 9) / 10], foldv[nfold - 1], (long long)nfold,
                  (long long)nmulti, 100.0 * (double)nmulti / (double)(nsrf ? nsrf : 1));
            }
            free(foldv);
            free(aw2);
            printf("   Ф8' ПОКРЫТИЕ: ячеек среза отображено %lld из %d; листьев-ИЗЛУЧАТЕЛЕЙ %lld, "
                   "листьев-ЗАСЛОНОВ без среза %lld, всего листьев %d; из излучателей ГЛЯНЦЕВЫХ "
                   "%lld\n",
                   (long long)nmap, S.n, (long long)nsrf, (long long)nblk, TD.nleaf,
                   (long long)nglos);
          }
          double t_build3 = now_s() - tb3;
          /* ПРОХОДЫ ПО НАПРАВЛЕНИЯМ. Ld обнуляется на каждое направление: они
           * независимы, и `vis` метит уже посчитанные — как `open < 0` у открытости. */
          double si3 = 0.0, sd3 = 0.0;
          double t_pass = 0.0;
          for (int pass = 0; pass < g_dpass; pass++) {
            for (int32_t i = 0; i < 3 * S.n; i++)
              D.Eind[i] = 0.0;
            for (int32_t i = 0; i < S.n; i++)
              wchk[i] = 0.0;
            double tp = now_s();
            for (int d = 0; d < DR.n; d++) {
              double om[3] = {DR.ox[d], DR.oy[d], DR.oz[d]};
              memset(D.Ld, 0, 3 * (size_t)TD.n * sizeof *D.Ld);
              memset(D.vis, 0, (size_t)TD.n);
              dsweep_rec(&TD, &D, 0, om, DR.w[d], g_dirsall);
              /* А937: ПРОВЕРКА НОРМИРОВКИ ЗАКОНОМ. `Σ_d w_d max(0, −ω·n)` обязана
               * быть `π` — иначе ошибка множителя смешается с физикой и проживёт,
               * как прожила ошибка §531. */
              for (int32_t i = 0; i < S.n; i++) {
                double cs = -(om[0] * snx[i] + om[1] * sny[i] + om[2] * snz[i]);
                if (cs > 0.0) wchk[i] += DR.w[d] * cs;
              }
            }
            t_pass = now_s() - tp;
            si3 = sd3 = 0.0;
            for (int32_t i = 0; i < S.n; i++)
              for (int k = 0; k < 3; k++) {
                double e3 =
                    D.Eind[3 * (size_t)i + (size_t)k] * (alb0 ? 0.0 : alb(&m, S.c[i].mat, k));
                if (pass + 1 == g_dpass) ind[3 * (size_t)i + (size_t)k] = (float)e3;
                si3 += e3;
                sd3 += (double)irr[3 * (size_t)i + (size_t)k];
              }
            double wmn = 1e300, wmx = -1e300, wav = 0.0;
            for (int32_t i = 0; i < S.n; i++) {
              if (wchk[i] < wmn) wmn = wchk[i];
              if (wchk[i] > wmx) wmx = wchk[i];
              wav += wchk[i];
            }
            wav /= (double)(S.n ? S.n : 1);
            printf(
                "   Ф8' СВИП ПО НАПРАВЛЕНИЯМ: ND %d (nmu %d, nphi %d), отскок %d; дерево+раскладка "
                "%.1f мс, проход %.1f мс (%.0f нс на лист-направление); СУММА косвенного / прямого "
                "= %.4f%s\n",
                DR.n, g_dnmu, g_dnphi, pass + 1, t_build3 * 1e3, t_pass * 1e3,
                t_pass * 1e9 / ((double)DR.n * (double)(TD.nleaf ? TD.nleaf : 1)),
                si3 / (sd3 > 0.0 ? sd3 : 1.0),
                g_dirsall ? "  [НК dirsall]"
                          : (g_dblkopen ? "  [ВЕРХНЯЯ граница: заслоны без среза ПРОЗРАЧНЫ]" : ""));
            printf("      А937 НОРМИРОВКА: Σ w·max(0,−ω·n) = %.5f…%.5f, среднее %.5f (обязана быть "
                   "π = %.5f, отклонение среднего %.2f %%)\n",
                   wmn, wmx, wav, 3.14159265358979323846,
                   100.0 * (wav - 3.14159265358979323846) / 3.14159265358979323846);
            /* Ф11' (§545): РАСПРЕДЕЛЕНИЕ ДОЛИ ПЕРЕКРЫТИЯ. Без него «правило
             * изменило ответ» не отличить от «правило почти не сработало». */
            if (D.nfstat > 0) {
              int64_t nf = D.nfstat < (int64_t)HZ_FSTAT_CAP ? D.nfstat : (int64_t)HZ_FSTAT_CAP;
              double *fc = malloc((size_t)nf * sizeof *fc);
              if (fc == NULL) exit(1);
              memcpy(fc, D.fstat, (size_t)nf * sizeof *fc);
              qsort(fc, (size_t)nf, sizeof *fc, cmp_d);
              printf(
                  "      Ф11' ДОЛЯ ПЕРЕКРЫТИЯ f: p10 %.4f, p50 %.4f, p90 %.4f; в единицу упёрлось "
                  "%.2f %% случаев (выборка %lld из %lld)\n",
                  fc[nf / 10], fc[nf / 2], fc[(nf * 9) / 10],
                  100.0 * (double)D.nfone / (double)(D.nfstat ? D.nfstat : 1), (long long)nf,
                  (long long)D.nfstat);
              free(fc);
            }
            /* Ф15' (§561): БЮДЖЕТ ИЗЛУЧЕНИЯ. `Φ_аналит` — сколько поверхность
             * обязана излучить в полусферу (`Σ A_i·irr_i`); `Φ_факт` — сколько
             * излучено; `Φ_обрезано` — съеденное обрезкой `min(1, f)`.
             * ПЕРЕКОС РАСКЛАДКИ (А1009) считается ОТДЕЛЬНО: `Σ max(0, A_p − h²)` —
             * сколько площади лежит в листьях сверх их собственного сечения. Это и
             * есть болезнь; обрезка — лишь её следствие. */
            {
              double phan = 0.0, over = 0.0, atot2 = 0.0;
              for (int32_t i = 0; i < S.n; i++)
                for (int k = 0; k < 3; k++)
                  phan += sar[i] * (double)irr[3 * (size_t)i + (size_t)k];
              for (int32_t l3 = 0; l3 < TD.n; l3++) {
                if (TD.nd[l3].child0 >= 0 || D.srf[l3] != 1) continue;
                double hl = (double)TD.nd[l3].size * D.cw;
                atot2 += (double)D.Ap[l3];
                if ((double)D.Ap[l3] > hl * hl) over += (double)D.Ap[l3] - hl * hl;
              }
              printf("      Ф15' БЮДЖЕТ ИЗЛУЧЕНИЯ: Φ_аналит %.4e, Φ_факт %.4e, Φ_обрезано %.4e; "
                     "невязка %.2f %%%s\n",
                     phan, emact * 3.0, emcut * 3.0,
                     100.0 * (phan - emact * 3.0 - emcut * 3.0) / (phan > 0.0 ? phan : 1.0),
                     g_fnoclamp ? "  [fnoclamp: обрезка СНЯТА]" : "");
              printf("      Ф15' ПЕРЕКОС РАСКЛАДКИ (А1009): площади сверх сечения листа %.4e из "
                     "%.4e (%.1f %%)\n",
                     over, atot2, 100.0 * over / (atot2 > 0.0 ? atot2 : 1.0));
              if (D.nfraw > 0) {
                qsort(D.frawv, (size_t)D.nfraw, sizeof *D.frawv, cmp_d);
                printf("      Ф15' СЫРОЕ f СРЕДИ УПЁРШИХСЯ: p50 %.3f, p90 %.3f, max %.3f по %lld "
                       "случаям\n",
                       D.frawv[D.nfraw / 2], D.frawv[(D.nfraw * 9) / 10], D.frawv[D.nfraw - 1],
                       (long long)D.nfraw);
              }
              emact = emcut = 0.0;
              D.nfraw = 0;
            }
            D.nfstat = 0;
            D.nfone = 0;
            /* МНОГОКРАТНЫЕ ОТРАЖЕНИЯ БЕЗ МАТРИЦЫ (довод №3 §523): следующая
             * радиосити есть собранная облучённость на альбедо. Матрицы нет, есть
             * ещё один проход по тем же направлениям. */
            if (pass + 1 < g_dpass) {
              double *aw3 = calloc((size_t)TD.n, sizeof *aw3);
              if (aw3 == NULL) exit(1);
              memset(D.Bs, 0, 3 * (size_t)TD.n * sizeof *D.Bs);
              for (int32_t i = 0; i < S.n; i++) {
                int32_t l2 = slf[i];
                if (l2 < 0) continue;
                aw3[l2] += sar[i];
                for (int k = 0; k < 3; k++)
                  D.Bs[3 * (size_t)l2 + (size_t)k] +=
                      (float)(sar[i] * D.Eind[3 * (size_t)i + (size_t)k] * alb(&m, S.c[i].mat, k));
              }
              for (int32_t l2 = 0; l2 < TD.n; l2++)
                if (D.srf[l2] == 1 && aw3[l2] > 0.0)
                  for (int k = 0; k < 3; k++)
                    D.Bs[3 * (size_t)l2 + (size_t)k] =
                        (float)((double)D.Bs[3 * (size_t)l2 + (size_t)k] / aw3[l2] /
                                3.14159265358979323846);
              free(aw3);
            }
          }
          /* Ф14. (§557): при сличении свип не вливается в `irr` — иначе гатер, идущий
           * следом, считал бы отскок от уже подсвеченной поверхности. */
          if (g_cmpcell) {
            indsw = malloc(3 * (size_t)S.n * sizeof *indsw);
            if (indsw == NULL) exit(1);
            memcpy(indsw, ind, 3 * (size_t)S.n * sizeof *indsw);
            memset(ind, 0, 3 * (size_t)S.n * sizeof *ind);
          } else
            for (int32_t i = 0; i < 3 * S.n; i++)
              irr[i] += ind[i];
          free(D.Ld);
          free(D.Bs);
          free(D.Bn);
          free(D.Bsp);
          free(D.Bwi);
          free(D.Bns);
          free(D.Ap);
          free(D.fstat);
          free(D.frawv);
          free(D.srf);
          free(D.vis);
          free(D.cstart);
          free(D.clist);
          free(D.Eind);
          free(snx);
          free(sny);
          free(snz);
          free(sar);
          free(slf);
          free(wchk);
          stree_free(&TD);
          tr3_dirs_free(&DR);
        }
        /* §597: ПРЕЖНИЙ ПУТЬ ИДЁТ ТОЛЬКО ПРИ `indslice` (НК1). Иначе косвенный
         * уже пришёл из узлов, и второй сбор дал бы ДВОЙНОЙ УЧЁТ — а он
         * выглядит как «стало ярче», то есть как улучшение (А1032). */
        if (g_hgather > 0.0 && g_indslice) {
          /* Ф6. (§516): ОТСКОК ПО ИЕРАРХИИ ИЗЛУЧАТЕЛЕЙ. */
          etree ET;
          memset(&ET, 0, sizeof ET);
          double tb2 = now_s();
          etree_build(&ET, &S, &fr, irr, &m, 0, S.n, 0, lev);
          double t_build = now_s() - tb2;
          tb2 = now_s();
          int64_t nlink = 0;
          double sthru2 = 0.0, sall2 = 0.0;
          /* Р2 (§581): ДЕЛЕНИЕ ПО ПРИЁМНИКАМ — самое большое число в системе
           * (`16.8` с). Каждый приёмник пишет свой `ind[i]`, дерево излучателей
           * читается всеми и не меняется. Порядок сложения ВНУТРИ приёмника не
           * меняется, значит ответ побитово тот же.
           * СЧЁТЧИКИ — ПО ПОТОКАМ, А СВОДЯТСЯ В ФИКСИРОВАННОМ ПОРЯДКЕ: редукция
           * OpenMP отдала бы порядок планировщику, и число поехало бы от запуска к
           * запуску. Здесь оно воспроизводимо. */
          int nth = g_omp1 ? 1 : omp_get_max_threads();
          int64_t *plink = calloc((size_t)nth, sizeof *plink);
          int64_t *pnv = calloc((size_t)nth, sizeof *pnv);
          double *pthru = calloc((size_t)nth, sizeof *pthru);
          double *pall = calloc((size_t)nth, sizeof *pall);
          if (plink == NULL || pthru == NULL || pall == NULL || pnv == NULL) exit(1);
#pragma omp parallel for schedule(dynamic, 64) if (!g_omp1)
          for (int32_t i = 0; i < S.n; i++) {
            int th = g_omp1 ? 0 : omp_get_thread_num();
            double pi[3], nn2[3], acc2[3] = {0, 0, 0};
            hz_slice_vertex(&S, i, pi);
            for (int k = 0; k < 3; k++)
              pi[k] = fr.org[k] + pi[k] * fr.h;
            hz_slice_normal(&S, i, nn2);
            double rrecv = fr.h * (double)((int32_t)1 << (lev - (int)S.c[i].lvl));
            for (int q2 = 0; q2 < 6; q2++)
              hgather_rec(&ET, 0, q2, pi, nn2, g_hgather, 0.0, rrecv, &P, &fr, indvis, acc2,
                          &plink[th], &pthru[th], &pall[th], &pnv[th], NULL, NULL, 0, NULL);
            for (int k = 0; k < 3; k++)
              ind[3 * (size_t)i + (size_t)k] =
                  (float)(acc2[k] * (alb0 ? 0.0 : alb(&m, S.c[i].mat, k)));
          }
          int64_t pnear2 = 0;
          for (int t4 = 0; t4 < nth; t4++) {
            pnear2 += pnv[t4];
            nlink += plink[t4];
            sthru2 += pthru[t4];
            sall2 += pall[t4];
          }
          free(plink);
          free(pthru);
          free(pall);
          free(pnv);
          double t_g2 = now_s() - tb2;
          double sd2 = 0.0, si2 = 0.0, sda = 0.0, sia = 0.0, sarea = 0.0;
          for (int32_t i = 0; i < S.n; i++) {
            /* ВЗВЕШЕННОЕ ПЛОЩАДЬЮ — ЕДИНСТВЕННОЕ, ЧТО СРАВНИМО МЕЖДУ УРОВНЯМИ
             * (§603). Сумма ПО ЯЧЕЙКАМ инвариантом не является: число ячеек
             * меняется с уровнем вчетверо, и отношение поехало бы просто от
             * смены населения, а не от физики. `Σ E·A / Σ E·A` приближает
             * `∫E dA / ∫E dA`, а это свойство СЦЕНЫ. */
            double sidew = fr.h * (double)((int32_t)1 << (lev - (int)S.c[i].lvl));
            double aa = sidew * sidew;
            sarea += aa;
            for (int k = 0; k < 3; k++) {
              sd2 += (double)irr[3 * (size_t)i + (size_t)k];
              si2 += (double)ind[3 * (size_t)i + (size_t)k];
              sda += (double)irr[3 * (size_t)i + (size_t)k] * aa;
              sia += (double)ind[3 * (size_t)i + (size_t)k] * aa;
            }
          }
          printf(
              "   Ф6. ИЕРАРХИЧЕСКИЙ ОТСКОК: eps %.3f, узлов дерева %d, СВЯЗЕЙ %lld против %lld пар "
              "(в %.0f раз меньше); дерево %.1f мс, сбор %.1f мс; СУММА косвенного / прямого = "
              "%.4f%s\n",
              g_hgather, ET.n, (long long)nlink, (long long)S.n * (long long)S.n,
              (double)((long long)S.n * (long long)S.n) / (double)(nlink ? nlink : 1),
              t_build * 1e3, t_g2 * 1e3, si2 / (sd2 > 0.0 ? sd2 : 1.0),
              indvis ? " [С ЗАСЛОНАМИ]" : "");
          printf("      §603 ВЗВЕШЕННОЕ ПЛОЩАДЬЮ (сравнимо между уровнями): Σ E_ind·A / Σ E_dir·A "
                 "= %.4f; площадь среза %.1f м²\n",
                 sia / (sda > 0.0 ? sda : 1.0), sarea);
          /* А1036: доля связей, где поправка ближней зоны ЗНАЧИМА. Без неё «комната
           * почти не сдвинулась» смешивает «далёкое поле цело» с «ближнее чуть
           * уменьшилось». Счётчик ПОПОТОЧНЫЙ и сводится в фиксированном
           * порядке — как и остальные здесь. Гонку в горячем цикле оставлять
           * нельзя даже в диагностике: число стало бы невоспроизводимым. */
          printf("      §600 БЛИЖНЯЯ ЗОНА: связей с ЗНАЧИМОЙ поправкой (A > 0.01 pi r^2) около "
                 "%lld из %lld (%.2f %%)\n",
                 (long long)pnear2, (long long)nlink,
                 100.0 * (double)pnear2 / (double)(nlink ? nlink : 1));
          g_t_bounce = t_build + t_g2;
          if (indmeas || indvis)
            printf("      §474 СКВОЗЬ ЗАСЛОНЫ: %.2f %%\n",
                   100.0 * sthru2 / (sall2 > 0.0 ? sall2 : 1.0));
          for (int32_t i = 0; i < 3 * S.n; i++)
            irr[i] += ind[i];
          etree_free(&ET);
        } else {
          /* Ф3. (§510): ИЗЛУЧАТЕЛИ С ОГРУБЛЁННОГО СРЕЗА. Приёмнику нужна
           * подробность, излучателю — нет: дальняя стена светит как ОДНА площадка
           * со своей средней яркостью. Второй срез того же дерева с бо́льшим порогом
           * и есть эта огрублённая раздача; прямой свет на нём считается тем же
           * `front_direct` (ячеек мало, цена ничтожна). */
          hz_dcslice SE;
          float *irre = irr;
          const hz_dcslice *SRC2 = &S;
          if (g_emitthr > 0.0) {
            lodctx LE = LL;
            LE.thr = g_emitthr;
            if (hz_slice_init(&SE, lev) != HZ_DC_OK) exit(1);
            if (hz_slice_build(&SE, &T, &ht, lod_stop, &LE) != HZ_DC_OK) exit(1);
            for (int32_t i2 = 0; i2 < SE.n; i2++) {
              double vw2[3];
              hz_slice_vertex(&SE, i2, vw2);
              int32_t cl2[3];
              for (int a2 = 0; a2 < 3; a2++) {
                double f2 = floor(vw2[a2]);
                if (f2 < 0.0) f2 = 0.0;
                if (f2 > (double)(fr.n - 1)) f2 = (double)(fr.n - 1);
                cl2[a2] = (int32_t)f2;
              }
              const int32_t *ls2 = NULL;
              if (ct_list(&CT, cl2, &ls2) == 0) continue;
              int32_t mi2 = m.fm != NULL ? m.fm[ls2[0]] : 0;
              SE.c[i2].mat = (uint8_t)(mi2 < 255 ? mi2 : 255);
            }
            irre = malloc(3 * (size_t)SE.n * sizeof *irre);
            if (irre == NULL) exit(1);
            front_direct(&SE, &fr, &P, &AL, irre, 0.5, 1, NULL, &m);
            SRC2 = &SE;
            printf("   Ф3' ОГРУБЛЁННЫЕ ИЗЛУЧАТЕЛИ: порог %.2f, ячеек %d против %d приёмников\n",
                   g_emitthr, SE.n, S.n);
          }
          int32_t stride = 1;
          while ((int64_t)(SRC2->n / (stride > 0 ? stride : 1)) * (int64_t)S.n > 200000000LL)
            stride *= 2;
          if (g_gstride > 0) stride = g_gstride;
          double tb = now_s();
          int64_t nemit = 0;
          double sthru = 0.0, sall = 0.0;
          double walb = 0.0, wtot = 0.0;
          for (int32_t j = 0; j < SRC2->n; j += stride) {
            double ej[3] = {(double)irre[3 * (size_t)j], (double)irre[3 * (size_t)j + 1],
                            (double)irre[3 * (size_t)j + 2]};
            if (!(ej[0] + ej[1] + ej[2] > 0.0)) continue;
            nemit++;
            /* §529: средневзвешенное альбедо ИЗЛУЧАТЕЛЕЙ, взвешенное их же потоком.
             * Печатается затем, что предсказание П1 сделано именно через него: если
             * второй множитель лишний, отношение обязано вырасти ровно в `1/⟨ρ⟩`. */
            for (int k = 0; k < 3; k++) {
              wtot += ej[k];
              walb += ej[k] * alb(&m, S.c[j].mat, k);
            }
            double pj[3], nj[3];
            hz_slice_vertex(&S, j, pj);
            for (int k = 0; k < 3; k++)
              pj[k] = fr.org[k] + pj[k] * fr.h;
            hz_slice_normal(&S, j, nj);
            double cside = fr.h * (double)((int32_t)1 << (lev - (int)SRC2->c[j].lvl));
            double aj = cside * cside * (double)stride;
            for (int32_t i = 0; i < S.n; i++) {
              if (i == j) continue;
              double pi[3], ni[3], w[3], r2 = 0.0;
              hz_slice_vertex(&S, i, pi);
              for (int k = 0; k < 3; k++)
                pi[k] = fr.org[k] + pi[k] * fr.h;
              hz_slice_normal(&S, i, ni);
              for (int k = 0; k < 3; k++) {
                w[k] = pj[k] - pi[k];
                r2 += w[k] * w[k];
              }
              if (!(r2 > 0.0)) continue;
              double r = sqrt(r2);
              double ci = (w[0] * ni[0] + w[1] * ni[1] + w[2] * ni[2]) / r;
              double cj = -(w[0] * nj[0] + w[1] * nj[1] + w[2] * nj[2]) / r;
              if (!(ci > 0.0) || !(cj > 0.0)) continue;
              /* §474: СКОЛЬКО КОСВЕННОГО ПРИХОДИТ СКВОЗЬ СТЕНЫ. Приближение (1)
               * названо в коде с самого начала, но НЕ ИЗМЕРЕНО ни разу; в закрытой
               * комнате оно перестаёт быть безобидным — наружная сторона стены
               * светит внутрь. Здесь тем же маршем, что у прямого света, считается
               * доля энергии, чей путь пересекает занятую ячейку. Ключ `indvis`
               * её ЗАСЛОНЯЕТ, `indmeas` — только считает. */
              double ff = ci * cj * aj / (3.14159265358979323846 * r2);
              int blocked = 0;
              if (indmeas) blocked = shadowed(&P, &fr, pi, pj, 0.5);
              for (int k = 0; k < 3; k++) {
                double e = ej[k] * (alb0 ? 0.0 : (g_emitalb2 ? alb(&m, S.c[j].mat, k) : 1.0)) * ff *
                           alb(&m, S.c[i].mat, k);
                if (indmeas) {
                  sall += e;
                  if (blocked) sthru += e;
                }
                if (blocked && indvis) continue;
                ind[3 * (size_t)i + (size_t)k] += (float)e;
              }
            }
          }
          tb = now_s() - tb;
          /* ПРИЁМКА: отношение косвенного к прямому обязано быть порядка альбедо. */
          double sd = 0.0, si = 0.0;
          for (int32_t i = 0; i < S.n; i++)
            for (int k = 0; k < 3; k++) {
              sd += (double)irr[3 * (size_t)i + (size_t)k];
              si += (double)ind[3 * (size_t)i + (size_t)k];
            }
          printf("   ОТСКОК: %.1f с (в %.0f раз дороже прямого света), излучателей %lld из %d "
                 "(прореживание %d); СУММА косвенного / прямого = %.4f\n",
                 tb, tb / (t_dir > 0.0 ? t_dir : 1.0), (long long)nemit, S.n, stride,
                 si / (sd > 0.0 ? sd : 1.0));
          printf("      §529 АЛЬБЕДО ИЗЛУЧАТЕЛЕЙ, взвешенное потоком: %.4f%s\n",
                 walb / (wtot > 0.0 ? wtot : 1.0),
                 g_emitalb2 ? "; ВТОРОЙ множитель ВОЗВРАЩЁН (emitalb2, СТАРОЕ НЕВЕРНОЕ)" : "");
          if (indmeas)
            printf("      §474 СКВОЗЬ ЗАСЛОНЫ: %.2f %% энергии отскока идёт путём, пересекающим "
                   "занятую ячейку%s\n",
                   100.0 * sthru / (sall > 0.0 ? sall : 1.0), indvis ? " (и ОТБРОШЕНА)" : "");
          /* ---- Ф14' (§557): ПОЯЧЕЕЧНОЕ СЛИЧЕНИЕ СВИПА С ГАТЕРОМ ---- */
          /* СУММАРНОЕ ЧИСЛО НЕ РАЗЛИЧАЕТ ДВЕ БОЛЕЗНИ: постоянный множитель (тогда
           * это ошибка нормировки, и схема ни при чём) и потерю с расстоянием
           * (тогда виноват перенос). Спутать их — потерять целый шаг на постройку
           * схемы, которая не нужна; в проекте так уже выходило трижды (А933,
           * А955, А994).
           * ОТНОШЕНИЕ С НУЛЯМИ — НЕ ВЕЛИЧИНА (А999): берутся ячейки, где ГАТЕР выше
           * порога от собственной медианы, а выброшенное считается тремя
           * счётчиками, а не прячется.
           * ГАТЕР — НЕ ИСТИНА, А ВТОРАЯ СХЕМА (А1001): его собственный разброс
           * меряется тем же прибором через `hgather=` и `stride`. */
          if (g_cmpcell && indsw != NULL) {
            const float *A1 = indsw, *B1 = g_cmpself ? indsw : ind;
            double *gv = malloc((size_t)S.n * sizeof *gv);
            double *rt = malloc((size_t)S.n * sizeof *rt);
            if (gv == NULL || rt == NULL) exit(1);
            int64_t ng = 0;
            for (int32_t i = 0; i < S.n; i++) {
              double b = 0.0;
              for (int k = 0; k < 3; k++)
                b += (double)B1[3 * (size_t)i + (size_t)k];
              if (b > 0.0) gv[ng++] = b;
            }
            double gmed = 0.0;
            if (ng > 0) {
              qsort(gv, (size_t)ng, sizeof *gv, cmp_d);
              gmed = gv[ng / 2];
            }
            /* Порог назван ОТ ДАННЫХ, а не с потолка: тысячная медианы гатера. */
            double thr2 = 1e-3 * gmed;
            int64_t nboth0 = 0, ngz = 0, nsz = 0, nuse = 0;
            /* Корзины по расстоянию ДО БЛИЖАЙШЕГО ИЗЛУЧАТЕЛЯ ПО ПРЯМОЙ. Величина
             * названа честно (А1000): за стеной она даёт НИЖНЮЮ оценку длины пути,
             * и если зависимость на такой оси найдётся — вывод тем крепче. */
            static const double DB[4] = {0.5, 1.5, 3.0, 1e9};
            double bs[4] = {0, 0, 0, 0}, bn[4] = {0, 0, 0, 0};
            for (int32_t i = 0; i < S.n; i++) {
              double a = 0.0, b = 0.0;
              for (int k = 0; k < 3; k++) {
                a += (double)A1[3 * (size_t)i + (size_t)k];
                b += (double)B1[3 * (size_t)i + (size_t)k];
              }
              if (!(a > 0.0) && !(b > 0.0)) {
                nboth0++;
                continue;
              }
              if (!(b > thr2)) {
                if (a > 0.0) ngz++;
                continue;
              }
              if (!(a > 0.0)) nsz++;
              rt[nuse++] = a / b;
              double pw[3];
              hz_slice_vertex(&S, i, pw);
              for (int k = 0; k < 3; k++)
                pw[k] = fr.org[k] + pw[k] * fr.h;
              double dmin = 1e300;
              for (int32_t j = 0; j < S.n; j += 16) {
                if (!(irr[3 * (size_t)j] > 0.0f)) continue;
                double pj2[3], d2 = 0.0;
                hz_slice_vertex(&S, j, pj2);
                for (int k = 0; k < 3; k++) {
                  double dd = fr.org[k] + pj2[k] * fr.h - pw[k];
                  d2 += dd * dd;
                }
                if (d2 < dmin) dmin = d2;
              }
              dmin = sqrt(dmin);
              for (int bq = 0; bq < 4; bq++)
                if (dmin < DB[bq]) {
                  bs[bq] += a / b;
                  bn[bq] += 1.0;
                  break;
                }
            }
            if (nuse > 0) {
              qsort(rt, (size_t)nuse, sizeof *rt, cmp_d);
              double p10 = rt[nuse / 10], p50 = rt[nuse / 2], p90 = rt[(nuse * 9) / 10];
              printf("   Ф14' СВИП / ГАТЕР ПОЯЧЕЕЧНО%s: p10 %.4f, p50 %.4f, p90 %.4f, "
                     "p90/p10 = %.2f по %lld ячейкам\n",
                     g_cmpself ? " [cmpself: обязано быть 1.0000]" : "", p10, p50, p90,
                     p10 > 0.0 ? p90 / p10 : 0.0, (long long)nuse);
              printf("      ВЫБРОШЕНО: обе нули %lld, гатер ниже порога при ненулевом свипе %lld, "
                     "свип ноль при живом гатере %lld (порог %.3e = 1e-3 медианы)\n",
                     (long long)nboth0, (long long)ngz, (long long)nsz, thr2);
              printf("      ПО РАССТОЯНИЮ ДО БЛИЖАЙШЕГО ИЗЛУЧАТЕЛЯ (по прямой, НИЖНЯЯ оценка "
                     "пути):\n");
              static const char *DN[4] = {"< 0.5 м", "0.5…1.5 м", "1.5…3 м", "> 3 м"};
              for (int bq = 0; bq < 4; bq++)
                printf("         %-10s среднее отношение %.4f по %.0f ячейкам\n", DN[bq],
                       bn[bq] > 0.0 ? bs[bq] / bn[bq] : 0.0, bn[bq]);
            }
            free(gv);
            free(rt);
          }
          for (int32_t i = 0; i < 3 * S.n; i++)
            irr[i] += ind[i];
          if (indsw != NULL) {
            for (int32_t i = 0; i < 3 * S.n; i++)
              irr[i] += indsw[i];
            free(indsw);
            indsw = NULL;
          }
          if (g_emitthr > 0.0) {
            free(irre);
            hz_slice_free(&SE);
          }
        }
        free(ind);
      } /* конец диагностики кадра (§589) */

      /* ЦВЕТ ЯЧЕЙКИ КЛАДЁТСЯ В ИНДЕКС ПО КЛЮЧУ, чтобы растеризатор мог его взять
       * по ячейке многоугольника. Индекс ПЛОСКИЙ (отсортированные ключи +
       * двоичный поиск), а не дерево: у него нет ни спуска, ни владения. */
      double t_sort = now_s();
      uint64_t *key = malloc((size_t)S.n * sizeof *key);
      int32_t *ord = malloc((size_t)S.n * sizeof *ord);
      if (key == NULL || ord == NULL) exit(1);
      for (int32_t i = 0; i < S.n; i++) {
        key[i] = cellkey(&S.c[i]);
        ord[i] = i;
      }
      /* СОРТИРОВКА ЗА `n log n`, А НЕ ЗА `n²` (§590).
       *
       * ЗДЕСЬ СТОЯЛА ВСТАВОЧНАЯ СОРТИРОВКА, и в комментарии рядом было честно
       * написано, почему вход к ней не отсортирован: срез идёт в МОРТОНОВОМ
       * порядке, а ключ лексикографичен по `(lvl, x, y, z)` — это разные
       * порядки, да ещё `lvl` в старших битах. То есть вход был практически
       * случайным, и вставка вырождалась в `n²`.
       * ЗАМЕРЕНО (§590, Bistro, `256²`, `poly=0`, `129 836` ячеек): `1828` мс из
       * `2084` мс кадра — `88 %`. В однокадровом прогоне это тонуло в `38` с
       * постройки и потому не было видно ни разу.
       *
       * ТАЙ-БРЕЙК ПО `ord` — НЕ УКРАШЕНИЕ. Вставочная сортировка УСТОЙЧИВА, а
       * `qsort` нет; при равных ключах порядок решал бы, чей цвет достанется
       * ячейке, и картинка поехала бы. Сравнение вторым ключом по возрастанию
       * `ord` даёт РОВНО тот же результат, что устойчивая сортировка исходно
       * возрастающего `ord`, — то есть побитовость сохраняется по построению, а
       * не по надежде. */
      {
        cellkv *kv = malloc((size_t)S.n * sizeof *kv);
        if (kv == NULL) exit(1);
        for (int32_t i = 0; i < S.n; i++) {
          kv[i].k = key[i];
          kv[i].o = ord[i];
        }
        qsort(kv, (size_t)S.n, sizeof *kv, cmp_cellkv);
        for (int32_t i = 0; i < S.n; i++) {
          key[i] = kv[i].k;
          ord[i] = kv[i].o;
        }
        free(kv);
      }

      /* БЕЛАЯ ТОЧКА: перцентиль 99.5 по ЯЧЕЙКАМ СРЕЗА — то же правило, что в
       * hz_ppm_write, но применённое к населению, у которого оно осмысленно. */
      double white = 1.0;
      {
        float *tmpw = malloc((size_t)S.n * sizeof *tmpw);
        if (tmpw == NULL) exit(1);
        for (int32_t i = 0; i < S.n; i++) {
          float mx = irr[3 * (size_t)i];
          if (irr[3 * (size_t)i + 1] > mx) mx = irr[3 * (size_t)i + 1];
          if (irr[3 * (size_t)i + 2] > mx) mx = irr[3 * (size_t)i + 2];
          tmpw[i] = mx;
        }
        qsort(tmpw, (size_t)S.n, sizeof *tmpw, cmp_f);
        double w995 = (double)tmpw[(size_t)((double)S.n * 0.995)];
        if (w995 > 0.0) white = w995;
        free(tmpw);
      }
      t_sort = now_s() - t_sort;
      litctx LC;
      memset(&LC, 0, sizeof LC);
      LC.polysum = HZ_FNV_BASIS; /* начальное значение FNV-1a */
      LC.S = &S;
      LC.key = key;
      LC.ord = ord;
      LC.irr = irr;
      LC.fr = &fr;
      LC.white = white;
      LC.nocull = nocull;
      LC.uvs = uvs;
      LC.mesh = &m;
      LC.ct = &CT;
      LC.nmtl = m.nmtl;
      LC.pxrad = LL.pxrad;
      LC.nomip = g_texnomip;
      tr3_camera cam;
      if (tr3_camera_look(&cam, eyec, atc, upc, HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0,
                          resw, resh) == 0) {
        size_t np = (size_t)resw * (size_t)resh;
        double *zb = malloc(np * sizeof *zb);
        unsigned char *rgb = calloc(np * 3, 1);
        if (zb == NULL || rgb == NULL) exit(1);
        for (size_t i = 0; i < np; i++)
          zb[i] = 1e300;
        LC.cam = &cam;
        LC.z = zb;
        LC.rgb = rgb;
        LC.w = resw;
        LC.h = resh;
        int64_t ahist[10] = {0};
        double apix[10] = {0};
        LC.areahist = ahist;
        LC.areapix = apix;
        LC.defcol = calloc(np * 3, sizeof *LC.defcol);
        if (LC.defcol == NULL) exit(1);
        /* §812: буфер спек-члена той же кадровой размерности; в цикле ходьбы
         * переиспользуется с обнулением (кадровое поле, не сцено-якорное) */
        if (g_glet != NULL) {
          if (g_glspec == NULL) {
            g_glspec = calloc(np * 3, sizeof *g_glspec);
            if (g_glspec == NULL) exit(1);
          } else
            memset(g_glspec, 0, np * 3 * sizeof *g_glspec);
          g_glnpx1 = g_glnpx10 = g_glclip = g_glprune = g_glleaf = g_glmarch = 0;
          g_glsum = 0.0;
          g_glt812 = 0.0;
        }
        /* Ш8 (§575): ЗАГРУЗКА ТЕКСТУР. Имя из `map_Kd`, каталог — `ppm256` рядом
         * с текстурами сцены (готовит `scripts/tex_prep.sh`). Отсутствующая
         * текстура НЕ ошибка: материал остаётся одноцветным, и число таких
         * ПЕЧАТАЕТСЯ, а не замалчивается (А891). */
        if (uvs != NULL && !g_texflat) {
          /* ПОЛЯ РАЗМЕРОМ С КАДР — КАЖДЫЙ КАДР. Они и есть кадр, а не сцена. */
          LC.defuv = calloc(np * 2, sizeof *LC.defuv);
          LC.defmat = calloc(np, 1);
          if (LC.defuv == NULL || LC.defmat == NULL) exit(1);
        }
        /* А ЗАГРУЗКА ТЕКСТУР — ОДИН РАЗ ЗА ЗАПУСК (§589). В цикле ходьбы иначе
         * платились бы `105` файлов на кадр. */
        if (uvs != NULL && !g_texflat && !g_texloaded) {
          g_texloaded = 1;
          g_texrgb = calloc((size_t)m.nmtl, sizeof *g_texrgb);
          g_texw = calloc((size_t)m.nmtl, sizeof *g_texw);
          g_texh = calloc((size_t)m.nmtl, sizeof *g_texh);
          g_texmip = calloc((size_t)m.nmtl, sizeof *g_texmip);
          g_mipw = calloc((size_t)m.nmtl, sizeof *g_mipw);
          g_miph = calloc((size_t)m.nmtl, sizeof *g_miph);
          g_nmip = calloc((size_t)m.nmtl, sizeof *g_nmip);
          if (g_texmip == NULL || g_mipw == NULL || g_miph == NULL || g_nmip == NULL) exit(1);
          if (g_texrgb == NULL || g_texw == NULL || g_texh == NULL) exit(1);
          char dir[512];
          snprintf(dir, sizeof dir, "%s", argv[1]);
          char *sl = strrchr(dir, '/');
          if (sl != NULL)
            *sl = '\0';
          else
            dir[0] = '\0';
          int64_t nload = 0, nmiss = 0, tbytes = 0;
          double tt0 = now_s();
          for (int32_t mi2 = 0; mi2 < m.nmtl; mi2++) {
            if (m.mtl[mi2].tex[0] == '\0') continue;
            /* ПУТЬ В `map_Kd` ОТБРАСЫВАЕТСЯ, И РАЗДЕЛИТЕЛЬ ТАМ ОБРАТНЫЙ. У
             * Сан-Мигеля стоит `textures\individual_b.png` — экспортёр писал под
             * Windows. Берём только имя файла: каталог у нас свой (`ppm256`), и
             * доверять пути из чужого файла на недоверенном входе нельзя тем
             * более. Замерено: без этого нашлось `0` текстур из `271`. */
            const char *nm2 = m.mtl[mi2].tex;
            for (const char *s2 = nm2; *s2 != '\0'; s2++)
              if (*s2 == '/' || *s2 == '\\') nm2 = s2 + 1;
            char base[160];
            snprintf(base, sizeof base, "%s", nm2);
            char *dot = strrchr(base, '.');
            if (dot != NULL) *dot = '\0';
            /* ДВА МЕСТА ПОИСКА, И ЭТО НЕ ПЕРЕСТРАХОВКА. У Сан-Мигеля текстуры
             * лежат в `textures/`, у Bistro — в нескольких каталогах рядом со
             * сценой (`BuildingTextures`, `Street`, `Natural`…), и `map_Kd` там
             * ссылается через `..\`. Мы кладём переведённые в ОДИН плоский
             * `ppm256` у сцены, поэтому ищем по имени файла в обоих местах. */
            char path2[900];
            snprintf(path2, sizeof path2, "%s/textures/ppm256/%s.ppm", dir, base);
            if (hz_ppm_read(path2, &g_texrgb[mi2], &g_texw[mi2], &g_texh[mi2]) != 0)
              snprintf(path2, sizeof path2, "%s/ppm256/%s.ppm", dir, base);
            if (hz_ppm_read(path2, &g_texrgb[mi2], &g_texw[mi2], &g_texh[mi2]) == 0) {
              nload++;
              tbytes += (int64_t)g_texw[mi2] * g_texh[mi2] * 3;
              /* МИП-ПИРАМИДА строится сразу: коробка 2×2 на уровень, пока сторона
               * не станет единицей. Коробка, а не что-то умнее, — потому что
               * уровень всё равно интерполируется линейно, и лишняя точность
               * ниже кванта байта. */
              int lv2 = 1, w2 = g_texw[mi2], h2 = g_texh[mi2];
              while (w2 > 1 || h2 > 1) {
                w2 = w2 > 1 ? w2 / 2 : 1;
                h2 = h2 > 1 ? h2 / 2 : 1;
                lv2++;
              }
              g_nmip[mi2] = lv2;
              g_texmip[mi2] = calloc((size_t)lv2, sizeof *g_texmip[mi2]);
              g_mipw[mi2] = calloc((size_t)lv2, sizeof *g_mipw[mi2]);
              g_miph[mi2] = calloc((size_t)lv2, sizeof *g_miph[mi2]);
              if (g_texmip[mi2] == NULL || g_mipw[mi2] == NULL || g_miph[mi2] == NULL) exit(1);
              g_texmip[mi2][0] = g_texrgb[mi2];
              g_mipw[mi2][0] = g_texw[mi2];
              g_miph[mi2][0] = g_texh[mi2];
              for (int l2 = 1; l2 < lv2; l2++) {
                int pw = g_mipw[mi2][l2 - 1], ph = g_miph[mi2][l2 - 1];
                int cw = pw > 1 ? pw / 2 : 1, ch = ph > 1 ? ph / 2 : 1;
                unsigned char *dst = malloc((size_t)cw * (size_t)ch * 3);
                if (dst == NULL) exit(1);
                const unsigned char *src = g_texmip[mi2][l2 - 1];
                for (int y2 = 0; y2 < ch; y2++)
                  for (int x2 = 0; x2 < cw; x2++)
                    for (int c2 = 0; c2 < 3; c2++) {
                      int sx0 = (pw > 1) ? 2 * x2 : 0, sy0 = (ph > 1) ? 2 * y2 : 0;
                      int sx1 = (pw > 1) ? sx0 + 1 : sx0, sy1 = (ph > 1) ? sy0 + 1 : sy0;
                      size_t o0 = 3 * ((size_t)sy0 * (size_t)pw + (size_t)sx0) + (size_t)c2;
                      size_t o1 = 3 * ((size_t)sy0 * (size_t)pw + (size_t)sx1) + (size_t)c2;
                      size_t o2 = 3 * ((size_t)sy1 * (size_t)pw + (size_t)sx0) + (size_t)c2;
                      size_t o3 = 3 * ((size_t)sy1 * (size_t)pw + (size_t)sx1) + (size_t)c2;
                      unsigned s5 = (unsigned)src[o0] + src[o1] + src[o2] + src[o3];
                      dst[3 * ((size_t)y2 * (size_t)cw + (size_t)x2) + (size_t)c2] =
                          (unsigned char)(s5 / 4);
                    }
                g_texmip[mi2][l2] = dst;
                g_mipw[mi2][l2] = cw;
                g_miph[mi2][l2] = ch;
                tbytes += (int64_t)cw * ch * 3;
              }
            } else
              nmiss++;
          }
          printf("   Ш8 ТЕКСТУРЫ: загружено %lld, не найдено %lld, память %.1f МБ, за %.2f с\n",
                 (long long)nload, (long long)nmiss, (double)tbytes / 1048576.0, now_s() - tt0);
        }
        if (!g_walk) printf("   §600 БЕЛАЯ ТОЧКА (перцентиль 99.5 по срезу): %.6e\n", white);
        /* Кадровый контекст получает УКАЗАТЕЛИ на разово загруженное. Владения он
         * не берёт: освобождать их в конце кадра значило бы грузить их заново. */
        if (uvs != NULL && !g_texflat) {
          LC.texrgb = g_texrgb;
          LC.texw = g_texw;
          LC.texh = g_texh;
          LC.texmip = g_texmip;
          LC.mipw = g_mipw;
          LC.miph = g_miph;
          LC.nmip = g_nmip;
        }
        ta = now_s();
        /* §573: РАЗДЕЛЕНИЕ ОБХОДА И РАСТЕРИЗАЦИИ. «УБИВАЕТ» §572 сработало
         * (126 мс против порога 100), и условие требует профилировать, а не
         * догадываться. Здесь тот же обход гоняется с ПУСТЫМ обработчиком: его
         * время есть цена обхода дерева и критерия LOD, а разность — цена самой
         * растеризации. */
        /* §580: КОПИЯ КОНТЕКСТА С ОТСЕЧЕНИЕМ. Оригинал `LL` идёт в срез, который
         * кормит ПЕРЕНОС, и там фрустумное отсечение запрещено. */
        lodctx LLc = LL;
        LLc.cull = !g_nofrustum;
        LLc.cam = &cam;
        LLc.h = fr.h;
        for (int a = 0; a < 3; a++)
          LLc.org[a] = fr.org[a];
        /* §590: В ХОДЬБЕ ПРОФИЛЬНОГО ОБХОДА НЕТ. Это ЗАМЕР (§573) — второй
         * полный обход дерева с пустым обработчиком, нужный только чтобы отделить
         * цену обхода от цены рисования. Замерено: `141` мс на кадр при `512²`,
         * то есть `13 %` кадра платились за прибор. В однокадровом прогоне он
         * остаётся: там он и осмыслен. */
        double ta_w = now_s();
        int wrc0 = 0;
        if (!g_walk) wrc0 = hz_dc_walk(&T, lod_stop, &LLc, lit_none, &LC);
        double t_walk = now_s() - ta_w;
#ifdef HZ_DC_COUNT
        {
          extern long long hz_dc_n_cell, hz_dc_n_face, hz_dc_n_edge, hz_dc_n_leafish, hz_dc_n_stop,
              hz_dc_n_proc;
          printf(
              "      СЧЁТ ОБХОДА: cellProc %lld, faceProc %lld, edgeProc %lld, process_edge %lld; "
              "leafish %lld (из них до lod_stop дошло %lld)\n",
              hz_dc_n_cell, hz_dc_n_face, hz_dc_n_edge, hz_dc_n_proc, hz_dc_n_leafish,
              hz_dc_n_stop);
        }
#endif
        ta = now_s();
        /* Р3 (§581): обход СОБИРАЕТ, отрисовка идёт ПО ПОЛОСАМ параллельно. */
        if (!g_omp1) {
          LC.captris = 65536;
          LC.tris = malloc((size_t)LC.captris * sizeof *LC.tris);
          if (LC.tris == NULL) exit(1);
        }
        int wrc = hz_dc_walk(&T, lod_stop, &LLc, lit_poly, &LC);
        if (LC.tris != NULL) {
          int nb2 = omp_get_max_threads();
          int bh = (resh + nb2 - 1) / nb2;
#pragma omp parallel for schedule(static)
          for (int b2 = 0; b2 < nb2; b2++) {
            int y0b = b2 * bh, y1b = y0b + bh - 1;
            if (y1b >= resh) y1b = resh - 1;
            for (int64_t t5 = 0; t5 < LC.ntris; t5++) {
              /* Отсев по предвычисленному габариту строк — два сравнения вместо
               * повторного проецирования трёх вершин. Без него деление на полосы
               * не давало ничего: замерено `218 → 234` мс, то есть ХУЖЕ
               * однопоточного, потому что каждая полоса перепроецировала все
               * треугольники заново. */
              if (LC.tris[t5].iy1 < y0b || LC.tris[t5].iy0 > y1b) continue;
              lit_tri(&LC, LC.tris[t5].p, LC.tris[t5].col, LC.tris[t5].uv, LC.tris[t5].mat, y0b,
                      y1b);
            }
          }
          free(LC.tris);
          LC.tris = NULL;
        }
        lit_resolve(&LC, g_gamn);
        if (g_glet != NULL) {
          int64_t ncovg8 = 0;
          size_t npg8 = (size_t)resw * (size_t)resh;
          for (size_t i8 = 0; i8 < npg8; i8++)
            if (zb[i8] < 1e299) ncovg8++;
          printf("   §812 БЛИКИ: канал %.2f с; пикселей со спеком > 1 %% — %.2f %%, > 10 %% — "
                 "%.2f %% (закрытых %lld); в клипе каналов %lld; корзин-листьев %lld, маршей "
                 "%lld, отсечений %lld; Σ спек %.4e\n",
                 g_glt812, 100.0 * (double)g_glnpx1 / (double)(ncovg8 > 0 ? ncovg8 : 1),
                 100.0 * (double)g_glnpx10 / (double)(ncovg8 > 0 ? ncovg8 : 1), (long long)ncovg8,
                 (long long)g_glclip, (long long)g_glleaf, (long long)g_glmarch,
                 (long long)g_glprune, g_glsum);
        }
        printf("      §573 ПРОФИЛЬ РАСТРА: обход дерева с ПУСТЫМ обработчиком %.1f мс (код %d), "
               "обход+растр %.1f мс — значит сама растеризация %.1f мс\n",
               t_walk * 1e3, wrc0, (now_s() - ta) * 1e3, (now_s() - ta - t_walk) * 1e3);
        /* §585: механизм проверяется СЧЁТОМ, а выгода — временем, и путать их
         * нельзя. `узлов посчитано` — это число РАЗЛИЧНЫХ узлов, у которых ответ
         * критерия вычислен полностью; при `nomemo` считать нечем, и печатается
         * ноль, а не подставленное число вызовов. */
        printf("      §585 ПАМЯТЬ ОТВЕТА: КРИТЕРИЙ СЧИТАН ПОЛНОСТЬЮ %lld раз (различных узлов "
               "тронуто %lld), перевёрнуто ответов %lld; ПОРЯДКОВАЯ СУММА ПОТОКА МНОГОУГОЛЬНИКОВ "
               "%016llx\n",
               hz_dc_walk_memo_stops(), hz_dc_walk_memo_evals(), hz_dc_walk_memo_flips(),
               (unsigned long long)LC.polysum);
        double t_rast = now_s() - ta;
        int64_t ncov = 0;
        for (size_t i2 = 0; i2 < np; i2++)
          if (zb[i2] < 1e299) ncov++;
        printf("      §572 РАСТР: фрагментов прошло z %lld на %lld закрытых пикселей, ГЛУБИНА "
               "ПЕРЕКРЫТИЯ %.2f; таблица гаммы %d входов; КООРДИНАТА С ПОВЕРХНОСТИ у %lld "
               "пикселей (%.1f %%)\n",
               (long long)LC.nfrag, (long long)ncov, (double)LC.nfrag / (double)(ncov ? ncov : 1),
               g_gamn, (long long)LC.nsurfuv,
               100.0 * (double)LC.nsurfuv / (double)(ncov ? ncov : 1));
        {
          double tot4 = 0.0;
          for (int b5 = 0; b5 < 10; b5++)
            tot4 += apix[b5];
          printf("      РАЗМЕР ПОЛИГОНА НА ЭКРАНЕ (площадь в пикселях -> сколько их, и какую долю "
                 "экрана они кроют):\n        ");
          int lo4 = 1;
          for (int b5 = 0; b5 < 10; b5++) {
            if (ahist[b5] > 0)
              printf("%d..%d: %lld шт / %.0f %%   ", lo4, lo4 * 4 - 1, (long long)ahist[b5],
                     100.0 * apix[b5] / (tot4 > 0.0 ? tot4 : 1.0));
            lo4 *= 4;
          }
          printf("\n");
        }
        char path[256];
        /* Ш18 (§494): СГЛАЖИВАНИЕ. Растр идёт в `res`, а на диск пишется вдвое
         * меньше со свёрткой коробкой 2×2 — четыре пробы на пиксель. Это НЕ
         * полноценное сглаживание: края ГЕОМЕТРИИ остаются ступенчатыми на уровне
         * ячейки, сглаживается только край многоугольника. Так и называется. */
        int outw = ss2 ? resw / 2 : resw, outh = ss2 ? resh / 2 : resh;
        unsigned char *outrgb = rgb;
        if (ss2) {
          outrgb = malloc((size_t)outw * (size_t)outh * 3);
          if (outrgb == NULL) exit(1);
          for (int y = 0; y < outh; y++)
            for (int x = 0; x < outw; x++)
              for (int c = 0; c < 3; c++) {
                unsigned s4 = 0;
                for (int dy = 0; dy < 2; dy++)
                  for (int dx = 0; dx < 2; dx++)
                    s4 += rgb[3 * ((size_t)(2 * y + dy) * (size_t)resw + (size_t)(2 * x + dx)) +
                              (size_t)c];
                outrgb[3 * ((size_t)y * (size_t)outw + (size_t)x) + (size_t)c] =
                    (unsigned char)((s4 + 2u) / 4u);
              }
        }
        snprintf(path, sizeof path, "img/pfield_lit_L%d_%dx%d%s.ppm", lev, outw, outh,
                 xframe ? "_sweep" : ""); /* §758 */
        /* В ХОДЬБЕ КАДР НЕ ПИШЕТСЯ НА ДИСК: `786` КБ на кадр — это и лишняя
         * работа, и мусор в `img/`. Снимок делает отдельный запуск без `walk`. */
        int prc = g_walk ? 0 : hz_ppm_write_rgb(path, outrgb, outw, outh);
        /* §804/А1327: PFM — ТОТ ЖЕ КАДР ДО ТОН-МАППА. Источник — defcol
         * растеризатора (радианса как есть); тон-мапп (белая точка + гамма)
         * необратим, а PFM существует затем, чтобы диффузный инструмент читал
         * ЛИНЕЙНЫЕ данные (формат — канонический выход по image.h). При ss2 —
         * то же коробочное 2×2 среднее, что у PPM, но во FLOAT и без
         * округления до байта: квантование до усреднения теряло бы ровно тот
         * диапазон, ради которого PFM заводился. Незакрытые пиксели — нули из
         * calloc. Умолчание ВКЛ, НК-ключ xnopfm; в ходьбе не пишется, как и
         * PPM. */
        if (!g_walk && !g_nopfm) {
          double *pchan[3];
          for (int c = 0; c < 3; c++) {
            pchan[c] = malloc((size_t)outw * (size_t)outh * sizeof *pchan[c]);
            if (pchan[c] == NULL) exit(1);
          }
          for (int y = 0; y < outh; y++)
            for (int x = 0; x < outw; x++)
              for (int c = 0; c < 3; c++) {
                double s5 = 0.0;
                int ns5 = 0;
                for (int dy = 0; dy < 2; dy++)
                  for (int dx = 0; dx < 2; dx++) {
                    /* без ss2 четыре пробы читают ОДИН пиксель — среднее
                     * равно ему самому, отдельная ветка не нужна */
                    int rx = ss2 ? 2 * x + dx : x, ry = ss2 ? 2 * y + dy : y;
                    s5 +=
                        (double)LC.defcol[3 * ((size_t)ry * (size_t)resw + (size_t)rx) + (size_t)c];
                    ns5++;
                  }
                pchan[c][(size_t)y * (size_t)outw + (size_t)x] = s5 / (double)ns5;
              }
          char pathp[256];
          snprintf(pathp, sizeof pathp, "img/pfield_lit_L%d_%dx%d%s.pfm", lev, outw, outh,
                   xframe ? "_sweep" : "");
          int pfrc = hz_pfm_write(pathp, pchan[0], pchan[1], pchan[2], outw, outh);
          printf("      PFM: радианса ДО тон-маппа, %dx%d -> %s (код %d)\n", outw, outh, pathp,
                 pfrc);
          for (int c = 0; c < 3; c++)
            free(pchan[c]);
        }
        if (g_walk) {
          double tf = t_slice + t_dir + t_rast;
          /* ПОЛОЖЕНИЕ ПЕЧАТАЕТСЯ, А НЕ ПОДРАЗУМЕВАЕТСЯ: без него «управление не
           * работает» и «работает, но смотрю в стену» неотличимы, а проверять
           * придётся первым делом. */
          double t_wall = now_s() - t_iter0;
          double t_other = t_wall - tf - t_mat - t_sort - t_walk - t_ind;
          printf(
              "   кадр %6.0f мс (%.2f к/с) = срез %.0f + свет %.0f + растр %.0f + материал %.0f "
              "+ сортировка %.0f + КОСВЕННЫЙ %.0f + профильный обход %.0f + прочее %.0f; ячеек %d, "
              "многоугольников %lld; глаз (%.2f %.2f %.2f) взгляд (%.2f %.2f %.2f)\n",
              t_wall * 1e3, 1.0 / (t_wall > 0.0 ? t_wall : 1.0), t_slice * 1e3, t_dir * 1e3,
              t_rast * 1e3, t_mat * 1e3, t_sort * 1e3, t_ind * 1e3, t_walk * 1e3, t_other * 1e3,
              S.n, (long long)LC.nseen, eyec[0], eyec[1], eyec[2], atc[0] - eyec[0],
              atc[1] - eyec[1], atc[2] - eyec[2]);
          tf = t_wall; /* движение считается по НАСТОЯЩЕМУ времени кадра */
          fflush(stdout);
          walk_alive = walk_present(outrgb, outw, outh, &resw, &resh, eyec, atc, upc, t_prev);
          t_prev = tf;
        }
        if (ss2) free(outrgb);
        if (!g_walk) {
          printf("   ОТСЕЧЕНИЕ: пришло %lld многоугольников, отброшено %lld (%.1f %%)\n",
                 (long long)LC.nseen, (long long)LC.ncull,
                 100.0 * (double)LC.ncull / (double)(LC.nseen ? LC.nseen : 1));
          printf(
              "   КАДР СО СВЕТОМ %dx%d: срез %.1f мс (ячеек %d), ПРЯМОЙ СВЕТ %.1f мс (%.0f нс на "
              "ячейку), растеризация %.1f мс (код %d), ВСЕГО %.1f мс -> %s (код %d)\n",
              outw, outh, t_slice * 1e3, S.n, t_dir * 1e3, t_dir * 1e9 / (double)(S.n ? S.n : 1),
              t_rast * 1e3, wrc, (t_slice + t_dir + t_rast) * 1e3, path, prc);
        }
        g_t_fslice = t_slice;
        g_t_fdir = t_dir;
        g_t_fras = t_rast;
        g_t_frame = t_slice + t_dir + t_rast;
        free(zb);
        free(rgb);
        free(LC.defcol);
      }
      free(key);
      free(ord);
      free(irr2);
      free(irr);
      free(uvs);
      hz_slice_free(&S);
      /* ХВОСТ ЦИКЛА. Всё кадровое освобождено ВЫШЕ — иначе тысяча кадров съела бы
       * память срезом по `460` тыс. ячеек на каждый. `uvs` добавлен сюда же: в
       * однократном прогоне он не освобождался вовсе, и это была утечка, безвредная
       * ровно потому, что программа тут же кончалась (§589, А1029). */
      if (!g_walk || !walk_alive) break;
    }
    if (g_walk) walk_close();
  }
  /* СВОДКА ПО СТАДИЯМ (08-12, требование пользователя «тайминги считать надо»).
   * Печатается ВСЕГДА, а не только в режиме `render`: до неё числа лежали
   * россыпью по два десятка строк, и сложить их глазом было нельзя.
   * ЧТО ЗДЕСЬ ЕСТЬ И ЧЕГО НЕТ. `ПОСТРОЙКА` — разовая работа над сценой
   * (занятость, рёбра, дерево); `КАДР` — то, что платится КАЖДЫЙ раз. Смешивать
   * их — то же, что смешивать холодный старт с установившимся, а это запрещено
   * прямо (`CLAUDE.md`). Отскок стоит отдельной строкой: он платится за кадр
   * только у динамического света. */
  {
    double t_build_all = t_occ + t_edges + t_tree;
    printf("== СВОДКА ВРЕМЕНИ (%s%s)\n"
           "   ПОСТРОЙКА СЦЕНЫ: занятость %.2f с + рёбра %.2f с + дерево %.2f с = %.2f с\n"
           "   ЗА КАДР: %.1f мс (срез %.1f + прямой свет %.1f + растеризация %.1f)\n"
           "   ОТСКОК:  %.1f мс (за кадр только при динамическом свете)\n"
           "   ВСЕГО ОТ ЗАПУСКА ДО КАДРА: %.2f с\n",
           g_render ? "режим render" : "полный, с диагностикой", nonrm ? ", НОРМАЛИ ОСЕВЫЕ" : "",
           t_occ, t_edges, t_tree, t_build_all, g_t_frame * 1e3, g_t_fslice * 1e3, g_t_fdir * 1e3,
           g_t_fras * 1e3, g_t_bounce * 1e3, now_s() - g_t0);
  }
  ct_free(&CT);
  hz_dc_free(&T);
  hz_htab_free(&ht);
  opyr_free(&P);
  hz_obj_free(&m);
  /* §806/А1339-А1340: exit-путь. Маска заливки жила до конца main (комментарий
   * на месте присвоения solidmask = fl обещал это освобождение с самого
   * заведения, но free не было); Ш15-пирамида xradlv — та же история: последний
   * читатель — секция картинки. xradlv[lev] == xrad, отдельных free не надо. */
  free(solidmask);
  if (xradlv != NULL) {
    for (int l = 0; l <= lev; l++)
      free(xradlv[l]);
    free(xradlv);
  }
  return 0;
}
