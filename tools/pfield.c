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
                          vrec *out, int64_t cap, int64_t *n) {
  if (t->nd[ni].flags & HZ_DC_HASVERT) {
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
    collect_verts(t, t->nd[ni].child0 + i, clo, half, out, cap, n);
  }
}

static void surf_err(const hz_dctree *T, const frame *fr, const hz_objmesh *m, trilist_fn tl,
                     void *tctx, const char *tag) {
  int64_t nv = 0;
  collect_verts(T, 0, (int32_t[3]){0, 0, 0}, fr->n, NULL, 0, &nv);
  if (nv == 0) return;
  double *er = malloc((size_t)nv * sizeof *er);
  vrec *vr = malloc((size_t)nv * sizeof *vr);
  if (er == NULL || vr == NULL) {
    free(er);
    free(vr);
    return;
  }
  int64_t k = 0;
  collect_verts(T, 0, (int32_t[3]){0, 0, 0}, fr->n, vr, nv, &k);
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
    printf("   %s ОШИБКА ПОВЕРХНОСТИ (ОДНОСТОРОННЯЯ, А598): p50 %.4f p90 %.4f max %.4f м, по "
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

/* --- 4. главная ------------------------------------------------------------ */

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "pfield ФАЙЛ.obj МАСШТАБ [lev=N] [nonrm] [occdump=ПУТЬ] [polydump=ПУТЬ]\n");
    return 2;
  }
  int lev = 6, nonrm = 0;
  const char *occdump = NULL, *polydump = NULL;
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0) lev = (int)strtol(argv[i] + 4, NULL, 10);
    /* НЕГАТИВНЫЙ КОНТРОЛЬ (§393): все нормали рёбер осевые. Вершины обязаны
     * остаться (ребро пересечено — вершина есть), а качество обязано упасть. */
    if (strcmp(argv[i], "nonrm") == 0) nonrm = 1;
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
  opyr P;
  int64_t nocc = opyr_build(&P, &m, &fr);
  double t_occ = now_s() - t0;
  printf(
      "   занятость за %.2f с: занятых ячеек %lld (%.3f %%), пирамида %.1f МБ\n", t_occ,
      (long long)nocc, 100.0 * (double)nocc / (double)((size_t)fr.n * (size_t)fr.n * (size_t)fr.n),
      (double)hz_occ_bytes((size_t)fr.n * (size_t)fr.n * (size_t)fr.n) * (8.0 / 7.0) / 1048576.0);
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
  collect_verts(&T, 0, zero, fr.n, NULL, 0, &nv);
  printf("   ВЕРШИН ВЫДАНО %lld; рёбер/вершину %.3f\n", (long long)nv,
         (double)ht.n / (double)(nv ? nv : 1));
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

  ct_free(&CT);
  hz_dc_free(&T);
  hz_htab_free(&ht);
  opyr_free(&P);
  hz_obj_free(&m);
  return 0;
}
