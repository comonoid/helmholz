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

static int cmp_f(const void *x, const void *y) {
  float a = *(const float *)x, b = *(const float *)y;
  return (a < b) ? -1 : ((a > b) ? 1 : 0);
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
        out[*n].vx[a] = t->nd[ni].vx[a];
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
  double thr;    /* порог, пикселей */
  int64_t hist[HZ_DC_MAX_LOG2SIZE + 2];
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
  if (!(d2 > 0.0)) return 0;
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
static double alb(const hz_objmesh *m, uint8_t mi, int k) {
  if (m->mtl == NULL || mi >= m->nmtl) return 0.5;
  double a = m->mtl[mi].kd3[k];
  return a >= 0.0 && a <= 1.0 ? a : 0.5;
}

typedef struct {
  double c[3]; /* центр площадки */

  double u[3], v[3]; /* полуоси */
  double rgb[3];     /* сила по каналам */
} arealight;

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
  double t = skip;
  while (t < len) {
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

/* Прямая облучённость ячейки среза по трём каналам. Площадка берётся четырьмя
 * образцами — число НАЗВАНО, а не подобрано: это углы, то есть худший случай для
 * полутени, и увеличение его только сгладит край. */
#define HZ_LIGHT_SAMPLES 4

static void front_direct(const hz_dcslice *S, const frame *fr, const opyr *P, const arealight *L,
                         float *irr, double stepfrac, int hier, int64_t *nstep,
                         const hz_objmesh *A) {
  static const double su[HZ_LIGHT_SAMPLES] = {-0.5, 0.5, -0.5, 0.5};
  static const double sv[HZ_LIGHT_SAMPLES] = {-0.5, -0.5, 0.5, 0.5};
  for (int32_t i = 0; i < S->n; i++) {
    double p[3], n[3];
    hz_slice_vertex(S, i, p);
    for (int k = 0; k < 3; k++)
      p[k] = fr->org[k] + p[k] * fr->h;
    hz_slice_normal(S, i, n);
    double acc[3] = {0, 0, 0};
    for (int s = 0; s < HZ_LIGHT_SAMPLES; s++) {
      double q[3], w[3], r2 = 0.0;
      for (int k = 0; k < 3; k++) {
        q[k] = L->c[k] + L->u[k] * su[s] + L->v[k] * sv[s];
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
      double g = cosr / r2 / (double)HZ_LIGHT_SAMPLES;
      for (int k = 0; k < 3; k++)
        acc[k] += L->rgb[k] * g * alb(A, S->c[i].mat, k);
    }
    for (int k = 0; k < 3; k++)
      irr[3 * (size_t)i + (size_t)k] = (float)acc[k];
  }
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
#define HZ_SWEEP_DROP 2

typedef struct {
  unsigned char *vis; /* на ячейку: булева видимость (осевое и направленное правила) */
  float *open;        /* Ш5а2: ДОЛЯ ОТКРЫТОСТИ, переносимая с весами граней */
  int32_t n;          /* сторона грубой сетки */
  int drop;           /* lev − log2(n) */
  int axis;           /* НЕГАТИВНЫЙ КОНТРОЛЬ: прежнее осевое наследование */
  int frac;           /* Ш5а2: дробная открытость вместо булевой */
  int round01;        /* НЕГАТИВНЫЙ КОНТРОЛЬ Ш5а2: округлять F до 0/1 */
} sweepgrid;

static void sweep_free(sweepgrid *G) {
  free(G->vis);
  free(G->open);
  G->vis = NULL;
  G->open = NULL;
}

/* Один образец источника: заполнить видимость на всей грубой сетке. */
static void sweep_light(sweepgrid *G, const opyr *P, const frame *fr, const double q[3]) {
  int32_t n = G->n;
  size_t nc = (size_t)n * (size_t)n * (size_t)n;
  memset(G->vis, 0, nc);
  for (size_t i = 0; i < nc; i++)
    G->open[i] = 0.0f;
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
  /* Источник в координатах ГРУБОЙ сетки — к нему и строится направление. */
  double sc[3];
  for (int k = 0; k < 3; k++)
    sc[k] = (q[k] - fr->org[k]) / (fr->h * (double)((int32_t)1 << G->drop));
  /* ВОСЕМЬ ОКТАНТОВ. Внутри октанта каждая ось идёт ОТ источника, поэтому сосед
   * со стороны источника уже посчитан — это и есть условие свипа. */
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
  return G->frac ? (double)G->open[ci] : (G->vis[ci] ? 1.0 : 0.0);
}

/* Прямая облучённость СВИПОМ. Отличие от `front_direct` только в том, откуда
 * берётся затенение; геометрия (косинус, `1/r²`, образцы площадки) та же — иначе
 * сверка мерила бы разницу формул, а не разницу механизмов. */
static void front_sweep(const hz_dcslice *S, const frame *fr, const opyr *P, const arealight *L,
                        float *irr, const hz_objmesh *A, double *t_sweep, double *t_gather,
                        int axismode, int fracmode, int round01) {
  static const double su[HZ_LIGHT_SAMPLES] = {-0.5, 0.5, -0.5, 0.5};
  static const double sv[HZ_LIGHT_SAMPLES] = {-0.5, -0.5, 0.5, 0.5};
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
  if (G.vis == NULL || G.open == NULL) exit(1);
  for (int32_t i = 0; i < 3 * S->n; i++)
    irr[i] = 0.0f;
  *t_sweep = 0.0;
  *t_gather = 0.0;
  for (int sm = 0; sm < HZ_LIGHT_SAMPLES; sm++) {
    double q[3];
    for (int k = 0; k < 3; k++)
      q[k] = L->c[k] + L->u[k] * su[sm] + L->v[k] * sv[sm];
    double ta = now_s();
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
      double g = vis * cosr / r2 / (double)HZ_LIGHT_SAMPLES;
      for (int k = 0; k < 3; k++)
        irr[3 * (size_t)i + (size_t)k] += (float)(L->rgb[k] * g * alb(A, S->c[i].mat, k));
    }
    *t_gather += now_s() - ta;
  }
  sweep_free(&G);
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

static void lit_tri(litctx *L, const double p[3][3], const double col[3][3]) {
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
  double d21x = sx[1] - sx[0], d21y = sy[1] - sy[0];
  double d31x = sx[2] - sx[0], d31y = sy[2] - sy[0];
  double det = d21x * d31y - d21y * d31x;
  if (!(fabs(det) > 0.0)) return;
  for (int py = iy0; py <= iy1; py++)
    for (int px = ix0; px <= ix1; px++) {
      double qx = (double)px + 0.5 - sx[0], qy = (double)py + 0.5 - sy[0];
      double u = (qx * d31y - qy * d31x) / det, v = (qy * d21x - qx * d21y) / det;
      if (u < 0.0 || v < 0.0 || u + v > 1.0) continue;
      double zz = sz[0] + u * (sz[1] - sz[0]) + v * (sz[2] - sz[0]);
      size_t k = (size_t)py * (size_t)L->w + (size_t)px;
      if (zz >= L->z[k]) continue;
      L->z[k] = zz;
      for (int c = 0; c < 3; c++) {
        double e = col[0][c] + u * (col[1][c] - col[0][c]) + v * (col[2][c] - col[0][c]);
        double t = e / L->white;
        if (t < 0.0) t = 0.0;
        if (t > 1.0) t = 1.0;
        double g = pow(t, 1.0 / 2.2);
        L->rgb[3 * k + (size_t)c] = (unsigned char)(g * 255.0 + 0.5);
      }
    }
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

static int lit_poly(void *ctx, const hz_dcref *ref, const double (*v)[3], int nv) {
  litctx *L = (litctx *)ctx;
  double w[4][3], col[4][3];

  for (int i = 0; i < nv; i++)
    for (int c = 0; c < 3; c++)
      w[i][c] = L->fr->org[c] + v[i][c] * L->fr->h;
  L->nseen++;
  if (!L->nocull && lit_cull(L, w, nv)) {
    L->ncull++;
    return 0;
  }
  for (int i = 0; i < nv; i++) {
    int idx = lit_find(L, &ref[i]);
    for (int c = 0; c < 3; c++)
      col[i][c] = idx >= 0 ? (double)L->irr[3 * (size_t)idx + (size_t)c] : 0.0;
  }
  for (int i = 1; i + 1 < nv; i++) {
    double p3[3][3], c3[3][3];
    for (int c = 0; c < 3; c++) {
      p3[0][c] = w[0][c];
      p3[1][c] = w[i][c];
      p3[2][c] = w[i + 1][c];
      c3[0][c] = col[0][c];
      c3[1][c] = col[i][c];
      c3[2][c] = col[i + 1][c];
    }
    lit_tri(L, p3, c3);
  }
  return 0;
}
/* --- 4. главная ------------------------------------------------------------ */

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "pfield ФАЙЛ.obj МАСШТАБ [lev=N] [nonrm] [occdump=ПУТЬ] [polydump=ПУТЬ]\n");
    return 2;
  }
  int lev = 6, nonrm = 0, vq1 = 0, hit = 0, nofix = 0, lit = 0, res = 512, sweepaxis = 0;
  int nrmflip = 0, nonsum = 0, indvis = 0, indmeas = 0, dosolid = 0;
  int doxfer = 0, xfernosolid = 0, doxsweep = 0, nmu = 4, xcorner = 0, xinnerfluid = 0;
  int ss2 = 0;
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
    if (strcmp(argv[i], "indvis") == 0) {
      indmeas = 1;
      indvis = 1;
    }
    /* НЕГАТИВНЫЙ КОНТРОЛЬ §397: вершина среза в один бит на ось. */
    if (strcmp(argv[i], "vq1") == 0) vq1 = 1;
    /* Ш4: удар сферой, два случая врозь (А669). */
    if (strcmp(argv[i], "hit") == 0) hit = 1;
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
    if (strncmp(argv[i], "res=", 4) == 0) res = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "thr=", 4) == 0) lodthr = strtod(argv[i] + 4, NULL);
    /* Выгрузка занятости для сверки с эталоном Ш0 (`tools/poccref.c`). */
    if (strncmp(argv[i], "occdump=", 8) == 0) occdump = argv[i] + 8;
    /* ЭТАЛОН РЕГРЕССА ДЛЯ Ш2 (А666): выданные многоугольники снимаются здесь. */
    if (strncmp(argv[i], "polydump=", 9) == 0) polydump = argv[i] + 9;
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
  {
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
  if (occdump != NULL) {
    int wrc = hz_occ_write(occdump, lev, fr.n, P.b[lev]);
    printf("   ДАМП ЗАНЯТОСТИ -> %s (код %d)\n", occdump, wrc);
  }

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
    double eyec[3] = HZ_CFG_HALL_EYE;
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
  int64_t nv = 0;
  int32_t zero[3] = {0, 0, 0};
  int64_t nvl = 0;
  collect_verts(&T, 0, zero, fr.n, NULL, 0, &nv, 0);
  collect_verts(&T, 0, zero, fr.n, NULL, 0, &nvl, 1);
  printf("   ВЕРШИН: У ВСЕХ УЗЛОВ %lld (рёбер/вершину %.3f), У ЛИСТЬЕВ %lld "
         "(рёбер/вершину %.3f)\n",
         (long long)nv, (double)ht.n / (double)(nv ? nv : 1), (long long)nvl,
         (double)ht.n / (double)(nvl ? nvl : 1));
  /* СУММЫ НОРМАЛЕЙ ПЕЧАТАЮТСЯ ОТДЕЛЬНОЙ СТРОКОЙ, А НЕ ПРЯЧУТСЯ В «Б/узел»:
   * побочный массив в `sizeof(hz_dcnode)` не входит, и умолчание о нём читалось
   * бы как «правка Ш9 памяти не стоит». */
  printf("   ПАМЯТЬ: дерево DC %.1f МБ (%zu Б/узел), рёбра %.1f МБ, пирамида %.1f МБ, СУММЫ "
         "НОРМАЛЕЙ %.1f МБ (%.1f Б/узел)\n",
         (double)T.n * (double)sizeof(hz_dcnode) / 1048576.0, sizeof(hz_dcnode),
         (double)ht.n * (double)sizeof(hz_hedge) / 1048576.0,
         (double)hz_occ_bytes((size_t)fr.n * (size_t)fr.n * (size_t)fr.n) * (8.0 / 7.0) / 1048576.0,
         (double)T.nsumcap * 24.0 / 1048576.0,
         T.n > 0 ? (double)T.nsumcap * 24.0 / (double)T.n : 0.0);

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
    int frc = hz_dc_facets(&T, NULL, NULL, &ftab, &cmap);
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
    if (tr3_mesh_build(&mesh, &ot, &ofr) != 0) exit(1);
    double t_mesh = now_s() - tx;
    /* МАСКА ПОЛНЫХ ЯЧЕЕК — из заливки. У листа с поверхностью размер 1, у
     * пустого — класс однороден по построению, поэтому хватает ОДНОЙ пробы в
     * его нижнем углу. */
    uint8_t *solid = calloc((size_t)mesh.ncell, 1);
    if (solid == NULL) exit(1);
    int64_t nsolid = 0;
    for (int32_t ci = 0; ci < mesh.ncell; ci++) {
      size_t k = hz_occ_index(fr.n, mesh.clo[ci][0], mesh.clo[ci][1], mesh.clo[ci][2]);
      /* Ш16 (§486): НАРУЖНОЕ тоже сплошное. Снаружи замкнутой комнаты свету
       * взяться неоткуда, а флюид там даёт 458 м³ бесполезной работы и щепки с
       * φ = 1.9e+05. Это утверждение о СЦЕНЕ, а не приближение. */
      int cls = solidmask != NULL ? solidmask[k] : 1u;
      if (solidmask != NULL && (cls == 0u || (cls == 2u && !xinnerfluid)) &&
          !hz_occ_get(P.b[lev], k)) {
        solid[ci] = 1u;
        nsolid++;
      }
    }
    tx = now_s();
    tr3_cut cut;
    int crc = tr3_cut_build(&cut, &mesh, &ftab, &cmap, xfernosolid ? NULL : solid);
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
             "(многоугольник не поместился); ОБЪЁМ: куб %.3f, ФЛЮИД %.3f, материал %.3f м³%s\n",
             cut.nbad, cut.nse, cut.nsebig, vbox, vfl, vbox - vfl,
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
      double *sig_t = calloc((size_t)mesh.ncell, sizeof *sig_t);
      double *sig_s = calloc((size_t)mesh.ncell, sizeof *sig_s);
      double *phi = calloc((size_t)mesh.ncell * 4, sizeof *phi);
      if (frho == NULL || femit == NULL || sig_t == NULL || sig_s == NULL || phi == NULL) exit(1);
      for (int32_t i = 0; i < ftab.n; i++)
        frho[i] = 0.7;
      double lc[3], lu = 0.0, lv = 0.0;
      for (int k = 0; k < 3; k++)
        lc[k] = 0.5 * (lo[k] + hi[k]);
      /* ПОТОЛОК ПОЛОСТИ — тем же спуском от КАМЕРЫ вверх, что и у площадки в
       * ветви `lit` (§472): «пусто» в занятости значит «нет поверхности», а не
       * «нет материала», поэтому считать сверху нельзя. */
      double eyew[3] = HZ_CFG_HALL_EYE;
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
      int64_t nlit = 0;
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
      /* ВХОД РАЗВЁРТКИ, А НЕ ТОЛЬКО ЕЁ ВЫХОД (А807). Ложный ноль на выходе
       * неотличим от «источник не задан», пока не напечатан сам источник. */
      int64_t nlitfac = 0, nlitfluid = 0;
      for (int32_t i = 0; i < ftab.n; i++)
        if (femit[i] > 0.0) nlitfac++;
      for (int32_t k = 0; k < cut.nse; k++)
        if (femit[cut.se[k].facet] > 0.0 && cut.mvol[cut.se[k].cell][0][0] > 0.0) nlitfluid++;
      printf("      ВХОД РАЗВЁРТКИ: светящихся ФАСЕТОВ %lld из %d, светящихся элементов во "
             "ФЛЮИДНЫХ ячейках %lld из %lld; альбедо фасета %.2f, ординат %d\n",
             (long long)nlitfac, ftab.n, (long long)nlitfluid, (long long)nlit, frho[0],
             2 * nmu * 4 * nmu);
      tr3_dirs dirs;
      if (tr3_dirs_product(&dirs, nmu, nmu) != 0) exit(1);
      tr3_problem prob = {.m = &mesh,
                          .d = &dirs,
                          .cut = &cut,
                          .facet_rho = frho,
                          .facet_emit = femit,
                          .nfacet = ftab.n,
                          .sig_t = sig_t,
                          .sig_s = sig_s,
                          .limiter = 1};
      tr3_stats st;
      memset(&st, 0, sizeof st);
      double tsw = now_s();
      int src = tr3_sweep_solve(&prob, 30, 1e-4, phi, &st);
      tsw = now_s() - tsw;
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
             "невязка %.2e, за %.2f с\n",
             dirs.n, (long long)nlit, src, st.iters, st.resid, tsw);
      printf("      ЭНЕРГИЯ: втекло %.4e, вытекло %.4e, поглощено %.4e, баланс %.2e; в "
             "поверхности %.4e, из них %.4e; max φ %.4e, max исходящий радианс %.4e\n",
             st.pin, st.pout, st.pabs, st.balance, st.psin, st.psout, phimax, soutmax);
      printf("      ГДЕ МАКСИМУМ: в полости %.4e, вне её %.4e; флюидный объём ячейки с "
             "максимумом %.3e м³ (у целой ячейки %.3e)\n",
             phimax_in, phimax_out, volmin_at_max, pow((double)(1 << (lev - 6)) * fr.h, 3.0));
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
      free(st.bout);
      free(st.sout);
      tr3_dirs_free(&dirs);
      free(frho);
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
  celltris CT;
  ct_build(&CT, &m, &fr, P.b[lev], nocc);
  surf_err(&T, &fr, &m, ct_list, &CT, "");
  {
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
    double eye[3] = HZ_CFG_HALL_EYE, at[3] = HZ_CFG_HALL_AT, up[3] = HZ_CFG_UP;
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
  {
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

    /* СВЕРКА С ДЕРЕВОМ: срез обязан нести ТЕ ЖЕ вершины с точностью кванта.
     * Популяция — ячейки СРЕЗА (их столько же, сколько вершин у дерева), и это
     * сказано прямо: сверка ловит потерю точности, а не потерю ячеек. */
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
      printf("   СВЕРКА СРЕЗА С ДЕРЕВОМ: ячеек среза %d, вершин дерева %lld; МАКС СМЕЩЕНИЕ "
             "%.3e м при кванте %.3e м (отношение %.3f); УГОЛ КОДА НОРМАЛИ p90 %.3f, макс "
             "%.3f град\n",
             S.n, (long long)nvt, dmax, quant, dmax / quant, a90, amax);
    }
    hz_slice_free(&S);
  }

  /* ---- 4в. АДАПТИВНЫЙ СРЕЗ (Ш3): критерий, ошибка ДО ВЫДАННОЙ ПОВЕРХНОСТИ,
   *      ориентация на перепаде уровней, доля среза за кадр ---- */
  {
    lodctx L;
    memset(&L, 0, sizeof L);
    {
      double eyec[3] = HZ_CFG_HALL_EYE;
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
      double eyec[3] = HZ_CFG_HALL_EYE, atc[3] = HZ_CFG_HALL_AT, dir[3];
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
    double eyec[3] = HZ_CFG_HALL_EYE, atc[3] = HZ_CFG_HALL_AT;
    for (int a = 0; a < 3; a++)
      LH.eye[a] = (eyec[a] - fr.org[a]) / fr.h;
    LH.pxrad = (HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0) / 512.0;
    LH.thr = 1.0;

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

  /* ---- 4д. КАДР СО СВЕТОМ (Ш5, часть первая) ---- */
  if (lit) {
    lodctx LL;
    memset(&LL, 0, sizeof LL);
    double eyec[3] = HZ_CFG_HALL_EYE, atc[3] = HZ_CFG_HALL_AT, upc[3] = HZ_CFG_UP;
    for (int a = 0; a < 3; a++)
      LL.eye[a] = (eyec[a] - fr.org[a]) / fr.h;
    LL.pxrad = (HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0) / (double)res;
    LL.thr = lodthr;
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
     * Числа не магические: они выведены из ГАБАРИТА, а не подобраны на глаз. */
    arealight AL;
    for (int k = 0; k < 3; k++)
      AL.c[k] = 0.5 * (lo[k] + hi[k]);
    /* ПОТОЛОК ПОЛОСТИ, А НЕ ВЕРХ ГАБАРИТА (§472). Прежнее «на 10 см ниже
     * hi[1]» верно только для сцены-оболочки без толщины: у зала верх габарита
     * и есть потолок. У ЗАМКНУТОЙ сцены с толстыми стенами там материал, и
     * лампа оказывается замурованной — замерено: освещённых ячеек 0, кадр
     * чёрный целиком.
     * Здесь потолок ИЩЕТСЯ ПО ЗАНЯТОСТИ: спуск по центральному столбцу от верха
     * габарита до первой ПУСТОЙ ячейки. Порога нет, число не подбирается —
     * читается сетка. Если пустых нет вовсе (сцена сплошная), остаётся прежнее
     * правило, и это честнее, чем молча повесить лампу в никуда. */
    {
      int32_t cx = (int32_t)((AL.c[0] - fr.org[0]) / fr.h);
      int32_t cz = (int32_t)((AL.c[2] - fr.org[2]) / fr.h);
      if (cx < 0) cx = 0;
      if (cx >= fr.n) cx = fr.n - 1;
      if (cz < 0) cz = 0;
      if (cz >= fr.n) cz = fr.n - 1;
      /* СЧИТАТЬ НАДО ОТ ТОЧКИ, ЗАВЕДОМО ЛЕЖАЩЕЙ В ПОЛОСТИ, А НЕ СВЕРХУ.
       * Спуск сверху упирается не в потолок полости, а в ПУСТОТУ ВНУТРИ ПЛИТЫ:
       * занятость метит ячейки, ЗАДЕТЫЕ ТРЕУГОЛЬНИКАМИ, а внутренность
       * сплошного тела треугольников не содержит и потому «пуста». Замерено:
       * правило сверху дало потолок на 2.385 м при настоящем 2.2 м.
       * Точка внутри полости у нас есть по построению — это КАМЕРА. */
      int32_t ye = (int32_t)((LL.eye[1]));
      if (ye < 0) ye = 0;
      if (ye >= fr.n) ye = fr.n - 1;
      int32_t xe = (int32_t)LL.eye[0], ze = (int32_t)LL.eye[2];
      if (xe < 0) xe = 0;
      if (xe >= fr.n) xe = fr.n - 1;
      if (ze < 0) ze = 0;
      if (ze >= fr.n) ze = fr.n - 1;
      cx = xe;
      cz = ze;
      int32_t ytop = -1;
      for (int32_t y = ye; y < fr.n; y++)
        if (hz_occ_get(P.b[lev], hz_occ_index(fr.n, cx, y, cz)) != 0) {
          ytop = y - 1;
          break;
        }
      AL.c[1] = ytop >= 0 ? fr.org[1] + ((double)ytop + 0.5) * fr.h - 0.10 : hi[1] - 0.10;
      printf("   ИСТОЧНИК: потолок полости найден по занятости на y = %.3f м, площадка на "
             "%.3f м\n",
             ytop >= 0 ? fr.org[1] + ((double)ytop + 0.5) * fr.h : hi[1], AL.c[1]);
    }
    double half = 0.25 * (hi[0] - lo[0]);
    AL.u[0] = half;
    AL.u[1] = 0.0;
    AL.u[2] = 0.0;
    AL.v[0] = 0.0;
    AL.v[1] = 0.0;
    AL.v[2] = 0.25 * (hi[2] - lo[2]);
    /* ОКРАШЕН СОЗНАТЕЛЬНО (А757): три одинаковых числа выдавать за RGB нельзя. */
    AL.rgb[0] = 1.00;
    AL.rgb[1] = 0.92;
    AL.rgb[2] = 0.78;

    /* МАТЕРИАЛ ЯЧЕЙКИ СРЕЗА (Ш5б, §430). Байт `mat` перестаёт быть резервом:
     * берётся материал первого треугольника в ячейке (`celltris` уже построен).
     * Индекс, а не альбедо: в срезе лежит ИНДЕКС, полезная нагрузка — в таблице
     * (Р4). Материалов на сцене десятки, таблица горяча в кэше.
     * ЧЕГО ЭТО НЕ ДЕЛАЕТ: на границе двух материалов ячейка получает ОДИН из
     * них, а не оба — граница пройдёт по ячейкам среза, то есть с точностью
     * LOD. Для отскока это законно (энергия), для резкой границы текстуры —
     * нет; текстур пока и нет. */
    {
      int64_t nmat = 0;
      for (int32_t i = 0; i < S.n; i++) {
        int32_t cell[3] = {(int32_t)S.c[i].lo[0], (int32_t)S.c[i].lo[1], (int32_t)S.c[i].lo[2]};
        const int32_t *ls = NULL;
        if (ct_list(&CT, cell, &ls) == 0) continue;
        int32_t mi = m.fm != NULL ? m.fm[ls[0]] : 0;
        if (mi < 0 || mi >= m.nmtl) mi = 0;
        S.c[i].mat = (uint8_t)(mi < 255 ? mi : 255);
        nmat++;
      }
      printf("   МАТЕРИАЛ В СРЕЗЕ: назначен %lld ячейкам из %d, материалов в сцене %d\n",
             (long long)nmat, S.n, m.nmtl);
    }

    float *irr = malloc(3 * (size_t)S.n * sizeof *irr);
    if (irr == NULL) exit(1);
    ta = now_s();
    /* РАБОЧИЙ ПУТЬ — С ПОДЪЁМОМ (§426): тот же предикат, вчетверо дешевле.
     * Плоский марш остаётся АРБИТРОМ и зовётся ниже. */
    int64_t nstep_w = 0;
    front_direct(&S, &fr, &P, &AL, irr, 0.5, 1, &nstep_w, &m);
    double t_dir = now_s() - ta;
    /* Ш15 (§484): КАРТИНКА ИЗ РАЗВЁРТКИ. Если сработал `xsweep`, облучённость
     * берётся не маршем и не гатером, а ИСХОДЯЩИМ РАДИАНСОМ поверхностных
     * элементов, посчитанным уравнением переноса. Заслоны там по построению —
     * то, что маршем стоит 65 с (§474).
     * СВЯЗЬ ЯЧЕЕК ПРЯМАЯ: и срез, и сетка переноса адресуются одной сеткой
     * уровня `lev`, поэтому переклад — это чтение по координате, а не поиск. */
    if (xrad != NULL) {
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
    float *irr2 = malloc(3 * (size_t)S.n * sizeof *irr2);
    if (irr2 == NULL) exit(1);
    double t_sw = 0.0, t_ga = 0.0;
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
      int64_t ndiff = 0, nlit_r = 0, nlit_s = 0;
      double emax = 0.0;
      for (int32_t i = 0; i < S.n; i++) {
        double a = irr[3 * (size_t)i], b = irr2[3 * (size_t)i];
        if (a > 0.0) nlit_r++;
        if (b > 0.0) nlit_s++;
        if ((a > 0.0) != (b > 0.0)) ndiff++;
        double d = fabs(a - b);
        if (d > emax) emax = d;
      }
      printf("   СВИП: %.1f мс на %d образцов; сетка %d^3 = %lld ячеек, %.2f нс НА ЯЧЕЙКУ СЕТКИ; "
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
        if (A.pl_area == NULL || A.pl_cnt == NULL || A.cl_area == NULL || A.cl_cnt == NULL) exit(1);
        int32_t nskip = 0;
        int wrca = pass == 0 ? hz_dc_walk_stats(&T, NULL, NULL, area_emit, &A, &nskip)
                             : hz_dc_walk_stats(&T, lod_stop, &LL, area_emit, &A, &nskip);
        double h2 = fr.h * fr.h;
        const char *tag = pass == 0 ? "ПОЛНАЯ ГЛУБИНА" : "СРЕЗ";
        printf("   ПЛОЩАДЬ ВЫДАННОЙ ПОВЕРХНОСТИ [%s] (§446): веер от v0 %.5f м², от v1 %.5f м² "
               "(разность %.3e); ИСТИННАЯ по мешу %.5f м², отношение %.4f\n",
               tag, A.fan0 * h2, A.fan1 * h2, fabs(A.fan0 - A.fan1) * h2, atrue,
               A.fan0 * h2 / (atrue > 0.0 ? atrue : 1.0));
        printf("   ПО ОСЯМ [%s] (главная ось нормали, БЕЗ знака): площадь %.5f / %.5f / %.5f м²; "
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
          printf(" || показано %d корзин, %.1f %% площади оси\n", nsh, 100.0 * shown / A.ax[ax]);
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
          printf(" || показано %d слоёв, %.1f %% площади оси\n", nsh2, 100.0 * shown2 / A.ax[ax]);
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
            printf("   СУММА УГЛОВЫХ КОЭФФИЦИЕНТОВ ПО ВСЕМ %d ПЛОЩАДКАМ: среднее %.4f, p10 %.4f, "
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
              printf("   ЯЧЕЕК НА ГРАНЬ (по главной оси нормали): %lld %lld %lld %lld %lld %lld "
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
    if (ind == NULL) exit(1);
    {
      int32_t stride = 1;
      while ((int64_t)(S.n / (stride > 0 ? stride : 1)) * (int64_t)S.n > 200000000LL)
        stride *= 2;
      double tb = now_s();
      int64_t nemit = 0;
      double sthru = 0.0, sall = 0.0;
      for (int32_t j = 0; j < S.n; j += stride) {
        double ej[3] = {(double)irr[3 * (size_t)j], (double)irr[3 * (size_t)j + 1],
                        (double)irr[3 * (size_t)j + 2]};
        if (!(ej[0] + ej[1] + ej[2] > 0.0)) continue;
        nemit++;
        double pj[3], nj[3];
        hz_slice_vertex(&S, j, pj);
        for (int k = 0; k < 3; k++)
          pj[k] = fr.org[k] + pj[k] * fr.h;
        hz_slice_normal(&S, j, nj);
        double cside = fr.h * (double)((int32_t)1 << (lev - (int)S.c[j].lvl));
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
            double e = ej[k] * (alb0 ? 0.0 : alb(&m, S.c[j].mat, k)) * ff * alb(&m, S.c[i].mat, k);
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
      if (indmeas)
        printf("      §474 СКВОЗЬ ЗАСЛОНЫ: %.2f %% энергии отскока идёт путём, пересекающим "
               "занятую ячейку%s\n",
               100.0 * sthru / (sall > 0.0 ? sall : 1.0), indvis ? " (и ОТБРОШЕНА)" : "");
      for (int32_t i = 0; i < 3 * S.n; i++)
        irr[i] += ind[i];
    }
    free(ind);

    /* ЦВЕТ ЯЧЕЙКИ КЛАДЁТСЯ В ИНДЕКС ПО КЛЮЧУ, чтобы растеризатор мог его взять
     * по ячейке многоугольника. Индекс ПЛОСКИЙ (отсортированные ключи +
     * двоичный поиск), а не дерево: у него нет ни спуска, ни владения. */
    uint64_t *key = malloc((size_t)S.n * sizeof *key);
    int32_t *ord = malloc((size_t)S.n * sizeof *ord);
    if (key == NULL || ord == NULL) exit(1);
    for (int32_t i = 0; i < S.n; i++) {
      key[i] = cellkey(&S.c[i]);
      ord[i] = i;
    }
    /* Срез уже в мортоновом порядке; ключ (lvl, lo) монотонен по нему не всегда,
     * поэтому сортируется явно. */
    for (int32_t i = 1; i < S.n; i++) {
      uint64_t k = key[i];
      int32_t o = ord[i];
      int32_t j = i - 1;
      while (j >= 0 && key[j] > k) {
        key[j + 1] = key[j];
        ord[j + 1] = ord[j];
        j--;
      }
      key[j + 1] = k;
      ord[j + 1] = o;
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
    litctx LC = {&S, key, ord, irr, &fr, NULL, NULL, 0, 0, NULL, white, nocull, 0, 0};
    tr3_camera cam;
    if (tr3_camera_look(&cam, eyec, atc, upc, HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0, res,
                        res) == 0) {
      size_t np = (size_t)res * (size_t)res;
      double *zb = malloc(np * sizeof *zb);
      unsigned char *rgb = calloc(np * 3, 1);
      if (zb == NULL || rgb == NULL) exit(1);
      for (size_t i = 0; i < np; i++)
        zb[i] = 1e300;
      LC.cam = &cam;
      LC.z = zb;
      LC.rgb = rgb;
      LC.w = res;
      LC.h = res;
      ta = now_s();
      int wrc = hz_dc_walk(&T, lod_stop, &LL, lit_poly, &LC);
      double t_rast = now_s() - ta;
      char path[256];
      /* Ш18 (§494): СГЛАЖИВАНИЕ. Растр идёт в `res`, а на диск пишется вдвое
       * меньше со свёрткой коробкой 2×2 — четыре пробы на пиксель. Это НЕ
       * полноценное сглаживание: края ГЕОМЕТРИИ остаются ступенчатыми на уровне
       * ячейки, сглаживается только край многоугольника. Так и называется. */
      int outres = ss2 ? res / 2 : res;
      unsigned char *outrgb = rgb;
      if (ss2) {
        outrgb = malloc((size_t)outres * (size_t)outres * 3);
        if (outrgb == NULL) exit(1);
        for (int y = 0; y < outres; y++)
          for (int x = 0; x < outres; x++)
            for (int c = 0; c < 3; c++) {
              unsigned s4 = 0;
              for (int dy = 0; dy < 2; dy++)
                for (int dx = 0; dx < 2; dx++)
                  s4 += rgb[3 * ((size_t)(2 * y + dy) * (size_t)res + (size_t)(2 * x + dx)) +
                            (size_t)c];
              outrgb[3 * ((size_t)y * (size_t)outres + (size_t)x) + (size_t)c] =
                  (unsigned char)((s4 + 2u) / 4u);
            }
      }
      snprintf(path, sizeof path, "img/pfield_lit_L%d_%d.ppm", lev, outres);
      int prc = hz_ppm_write_rgb(path, outrgb, outres, outres);
      if (ss2) free(outrgb);
      printf("   ОТСЕЧЕНИЕ: пришло %lld многоугольников, отброшено %lld (%.1f %%)\n",
             (long long)LC.nseen, (long long)LC.ncull,
             100.0 * (double)LC.ncull / (double)(LC.nseen ? LC.nseen : 1));
      printf("   КАДР СО СВЕТОМ %d²: срез %.1f мс (ячеек %d), ПРЯМОЙ СВЕТ %.1f мс (%.0f нс на "
             "ячейку), растеризация %.1f мс (код %d), ВСЕГО %.1f мс -> %s (код %d)\n",
             res, t_slice * 1e3, S.n, t_dir * 1e3, t_dir * 1e9 / (double)(S.n ? S.n : 1),
             t_rast * 1e3, wrc, (t_slice + t_dir + t_rast) * 1e3, path, prc);
      free(zb);
      free(rgb);
    }
    free(key);
    free(ord);
    free(irr2);
    free(irr);
    hz_slice_free(&S);
  }
  ct_free(&CT);
  hz_dc_free(&T);
  hz_htab_free(&ht);
  opyr_free(&P);
  hz_obj_free(&m);
  return 0;
}
