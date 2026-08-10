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
    i0[c] = (int64_t)floor((a - fr->org[c]) / fr->h);
    i1[c] = (int64_t)floor((b - fr->org[c]) / fr->h);
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
    E->cap = nc;
  }
  double *d = E->v + 9 * (size_t)E->ntri;
  for (int k = 0; k < 3; k++) {
    d[k] = a[k];
    d[3 + k] = b[k];
    d[6 + k] = c[k];
  }
  int32_t ti = E->ntri++;
  (void)ref;
  (void)nv;
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

static int emesh_emit(void *ctx, const hz_dcref *ref, const double (*v)[3], int nv) {
  emesh *E = (emesh *)ctx;
  double a[3] = {v[0][0], v[0][1], v[0][2]};
  for (int i = 1; i + 1 < nv; i++)
    if (emesh_push(E, a, v[i], v[i + 1], ref, nv) != 0) return 1;
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
static int64_t emesh_flips(const emesh *E, celltris *CT, const frame *fr, const hz_objmesh *m,
                           int64_t *ncmp) {
  int64_t nf = 0, nc = 0;
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
    if (ct_list(CT, cell, &ls) == 0) continue;
    const double *A2, *B2, *C2;
    tri_verts(m, ls[0], &A2, &B2, &C2);
    double f1[3], f2[3], fn[3];
    for (int c = 0; c < 3; c++) {
      f1[c] = B2[c] - A2[c];
      f2[c] = C2[c] - A2[c];
    }
    fn[0] = f1[1] * f2[2] - f1[2] * f2[1];
    fn[1] = f1[2] * f2[0] - f1[0] * f2[2];
    fn[2] = f1[0] * f2[1] - f1[1] * f2[0];
    double d = nn[0] * fn[0] + nn[1] * fn[1] + nn[2] * fn[2];
    nc++;
    if (d < 0.0) nf++;
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
typedef struct {
  double c[3];       /* центр площадки */
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

/* Прямая облучённость ячейки среза по трём каналам. Площадка берётся четырьмя
 * образцами — число НАЗВАНО, а не подобрано: это углы, то есть худший случай для
 * полутени, и увеличение его только сгладит край. */
#define HZ_LIGHT_SAMPLES 4

static void front_direct(const hz_dcslice *S, const frame *fr, const opyr *P, const arealight *L,
                         float *irr, double stepfrac) {
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
      if (stepfrac > 0.0 && shadowed(P, fr, p, q, stepfrac)) continue;
      double g = cosr / r2 / (double)HZ_LIGHT_SAMPLES;
      for (int k = 0; k < 3; k++)
        acc[k] += L->rgb[k] * g;
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
                        float *irr, double *t_sweep, double *t_gather, int axismode, int fracmode,
                        int round01) {
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
        irr[3 * (size_t)i + (size_t)k] += (float)(L->rgb[k] * g);
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

static int lit_poly(void *ctx, const hz_dcref *ref, const double (*v)[3], int nv) {
  litctx *L = (litctx *)ctx;
  double w[4][3], col[4][3];
  for (int i = 0; i < nv; i++) {
    for (int c = 0; c < 3; c++)
      w[i][c] = L->fr->org[c] + v[i][c] * L->fr->h;
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
  int sweepfrac = 1, sweepr01 = 0;
  double lodthr = 1.0;
  const char *occdump = NULL, *polydump = NULL;
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0) lev = (int)strtol(argv[i] + 4, NULL, 10);
    /* НЕГАТИВНЫЙ КОНТРОЛЬ (§393): все нормали рёбер осевые. Вершины обязаны
     * остаться (ребро пересечено — вершина есть), а качество обязано упасть. */
    if (strcmp(argv[i], "nonrm") == 0) nonrm = 1;
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
    fr.org[c] = 0.5 * (lo[c] + hi[c]) - 0.5 * (double)fr.n * fr.h;
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
  printf("   ПАМЯТЬ: дерево DC %.1f МБ (%zu Б/узел), рёбра %.1f МБ, пирамида %.1f МБ\n",
         (double)T.n * (double)sizeof(hz_dcnode) / 1048576.0, sizeof(hz_dcnode),
         (double)ht.n * (double)sizeof(hz_hedge) / 1048576.0,
         (double)hz_occ_bytes((size_t)fr.n * (size_t)fr.n * (size_t)fr.n) * (8.0 / 7.0) /
             1048576.0);

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
      E0.gn = (int32_t)1 << (lev < HZ_EGRID ? lev : HZ_EGRID);
      E0.lev = lev;
      E0.tri = NULL;
      E0.head = malloc((size_t)E0.gn * (size_t)E0.gn * (size_t)E0.gn * sizeof *E0.head);
      E0.nxt = NULL;
      E0.ncap = 0;
      E0.nn = 0;
      if (E0.v == NULL || E0.head == NULL) exit(1);
      for (int64_t i = 0; i < (int64_t)E0.gn * E0.gn * E0.gn; i++)
        E0.head[i] = -1;
      if (hz_dc_walk(&T, NULL, NULL, emesh_emit, &E0) != HZ_DC_OK) exit(1);
      int64_t nc0 = 0, nf0 = emesh_flips(&E0, &CT, &fr, &m, &nc0);
      printf("   ЭТАЛОН ПОЛНОЙ ГЛУБИНЫ: треугольников %lld, ОБРАЩЁННЫХ %lld из %lld (%.3f %%) — "
             "это доля МОДЕЛИ, перепада уровней здесь нет\n",
             (long long)E0.ntri, (long long)nf0, (long long)nc0,
             100.0 * (double)nf0 / (double)(nc0 ? nc0 : 1));
      free(E0.v);
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
      EM.gn = (int32_t)1 << (lev < HZ_EGRID ? lev : HZ_EGRID);
      EM.lev = lev;
      EM.tri = NULL;
      EM.head = malloc((size_t)EM.gn * (size_t)EM.gn * (size_t)EM.gn * sizeof *EM.head);
      EM.nxt = NULL;
      EM.ncap = 0;
      EM.nn = 0;
      if (EM.v == NULL || EM.head == NULL) exit(1);
      for (int64_t i = 0; i < (int64_t)EM.gn * EM.gn * EM.gn; i++)
        EM.head[i] = -1;
      if (hz_dc_walk(&T, lod_stop, &L, emesh_emit, &EM) != HZ_DC_OK) exit(1);

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
      printf("   СРЕЗ порог %.2f пикс: ячеек %d, доля %.4f, уровни %d..%d медиана %d, сборка "
             "%.3f с, треугольников %lld; ОШИБКА ДО ВЫДАННОЙ ПОВЕРХНОСТИ p50 %.5f p90 %.5f "
             "p99 %.5f м; без поверхности %lld из %lld\n",
             L.thr, A.n, (double)A.n / (double)(nvl ? nvl : 1), lmin, lmax, lmed, ts,
             (long long)EM.ntri, p50, p90, p99, (long long)nmiss, (long long)(ne + nmiss));
      {
        int64_t nc2 = 0;
        int64_t nf2 = emesh_flips(&EM, &CT, &fr, &m, &nc2);
        printf("      ОБРАЩЁННЫХ (А729, дифференциально): %lld из %lld (%.3f %%)\n", (long long)nf2,
               (long long)nc2, 100.0 * (double)nf2 / (double)(nc2 ? nc2 : 1));
      }
      free(EM.v);
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
    if (hz_slice_build(&S, &T, &ht, lod_stop, &LL) != HZ_DC_OK) exit(1);
    double t_slice = now_s() - ta;

    /* ИСТОЧНИК: площадка под потолком зала. Габарит сцены известен, потолок —
     * его верх по оси Y; площадка ставится на 10 см ниже и в центре плана.
     * Числа не магические: они выведены из ГАБАРИТА, а не подобраны на глаз. */
    arealight AL;
    for (int k = 0; k < 3; k++)
      AL.c[k] = 0.5 * (lo[k] + hi[k]);
    AL.c[1] = hi[1] - 0.10;
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

    float *irr = malloc(3 * (size_t)S.n * sizeof *irr);
    if (irr == NULL) exit(1);
    ta = now_s();
    front_direct(&S, &fr, &P, &AL, irr, 0.5);
    double t_dir = now_s() - ta;
    /* А772/А775: ЭТАЛОН ПРОВЕРЯЕТСЯ САМ. Марш идёт шагом , и тонкий заслон
     * он может проскочить. Пересчёт вдвое мельче: если множество затенённых
     * почти не изменилось, эталон устойчив, и доли расхождения со свипом
     * говорят про свип. Если изменилось — все эти доли наполовину про эталон. */
    {
      float *irrf = malloc(3 * (size_t)S.n * sizeof *irrf);
      if (irrf == NULL) exit(1);
      double tf = now_s();
      front_direct(&S, &fr, &P, &AL, irrf, 0.25);
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

    /* СВИП — то, что Ш5 обязан измерить; луч выше остаётся ЭТАЛОНОМ (А763), и
     * сверка идёт ПОЯЧЕЕЧНО, а не по картинке. */
    float *irr2 = malloc(3 * (size_t)S.n * sizeof *irr2);
    if (irr2 == NULL) exit(1);
    double t_sw = 0.0, t_ga = 0.0;
    front_sweep(&S, &fr, &P, &AL, irr2, &t_sw, &t_ga, sweepaxis, sweepfrac, sweepr01);
    /* ДОЛЯ ОТКРЫТОСТИ, А НЕ ОБЛУЧЁННОСТЬ (§423, А781). Приёмка задана на долю,
     * поэтому нужен знаменатель — облучённость БЕЗ всякого затенения. Считается
     * третьим проходом, дешёвым (теста заслона в нём нет вовсе). */
    float *irro = malloc(3 * (size_t)S.n * sizeof *irro);
    if (irro == NULL) exit(1);
    front_direct(&S, &fr, &P, &AL, irro, -1.0);
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
    litctx LC = {&S, key, ord, irr, &fr, NULL, NULL, 0, 0, NULL, white};
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
      snprintf(path, sizeof path, "img/pfield_lit_L%d_%d.ppm", lev, res);
      int prc = hz_ppm_write_rgb(path, rgb, res, res);
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
