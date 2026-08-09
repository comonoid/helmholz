/* pfield — ИСТОЧНИК ЭРМИТОВА ПОЛЯ ИЗ МЕША (PLAN_ELEMENTS.md §363, §365, §366,
 * §370, §371).
 *
 * ЗАЧЕМ. Решение §363: представление сцены — эрмитово поле на восьмидереве, а не
 * индекс над треугольниками. Готовая машинерия (`src/cut/dc.c`) лежала без дела
 * ровно по одной причине: ей нужен ЗНАК «внутри/снаружи», а наши сцены —
 * незамкнутый суп треугольников. Знак строится заливкой пустоты от границы.
 *
 * ЗДЕСЬ ДВА НЕЗАВИСИМЫХ ИСТОЧНИКА ОДНОГО И ТОГО ЖЕ ПОЛЯ, И ЭТО НЕ ИЗБЫТОЧНОСТЬ
 * (Г42, §370). ПЛОТНЫЙ: три массива по `n³` и `3·(n+1)³`, опрос всех узлов
 * сетки, цена `O(8^L)`, предел `L = 6`. РАЗРЕЖЁННЫЙ: дерево занятости, в котором
 * узел дробится, только если его коробка, РАСШИРЕННАЯ НА ОДНУ ЯЧЕЙКУ, задета
 * треугольником; заливка идёт по граням между листьями этого дерева; знак,
 * пересечение и предикат «есть ли смена знака в коробке» отвечаются спуском.
 * Цена — по числу ПОСЕЩЁННЫХ УЗЛОВ. Плотный путь оставлен ЭТАЛОНОМ: при `L = 6`
 * ключ `both` гоняет оба и сличает вершины ПОБИТОВО, сопоставляя их по КОРОБКЕ
 * узла, а не по индексу (А613 — индекс есть функция формы дерева, а формы
 * разные).
 *
 * ЗАЛИВКА ПЕРЕКРЫВАЕТСЯ ЯЧЕЙКОЙ С ГЕОМЕТРИЕЙ, А НЕ ПЕРЕСЕЧЁННЫМ РЕБРОМ (А596).
 * Поверхность, целиком лежащая внутри ячейки и не задевающая её рёбер, рёбер не
 * пересекает и заливку не остановила бы: она протекла бы сквозь стену и объявила
 * комнату наружей. Это отказ в пользу пустоты, что запрещено (Г25).
 *
 * ЧЕГО ЭТОТ СТЕНД НЕ ПОКАЗЫВАЕТ. Ошибка поверхности ОДНОСТОРОННЯЯ (А598): ловит
 * потерю геометрии, не ловит появление лишней; её медиана сидит на полу метрики
 * (А615), поэтому приёмкой служит побитовое равенство вершин, а не перцентили.
 * Сверка двух путей возможна только там, где есть плотный эталон, то есть до
 * `L = 6` (А614).
 */
#include "cut/dc.h"
#include "image.h"
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

/* --- 1. ПЛОТНЫЙ ИСТОЧНИК (эталон, L <= 6) ---------------------------------- */

typedef struct {
  int32_t n;          /* сторона сетки ячеек */
  unsigned char *sgn; /* (n+1)^3: 1 — угол ВНУТРИ тела */
  double *et;         /* 3*(n+1)^3: доля вдоль ребра, <0 — пересечения нет */
  double *enrm;       /* 9*(n+1)^3 */
} src_ctx;

static int src_sign(void *ctx, const int32_t c[3]) {
  const src_ctx *S = (const src_ctx *)ctx;
  int32_t n1 = S->n + 1;
  if (c[0] < 0 || c[1] < 0 || c[2] < 0 || c[0] >= n1 || c[1] >= n1 || c[2] >= n1) return 0;
  return S->sgn[((size_t)c[2] * (size_t)n1 + (size_t)c[1]) * (size_t)n1 + (size_t)c[0]];
}

static int src_cross(void *ctx, const int32_t a[3], int axis, double *t, double nrm[3]) {
  const src_ctx *S = (const src_ctx *)ctx;
  int32_t n1 = S->n + 1;
  if (a[0] < 0 || a[1] < 0 || a[2] < 0 || a[0] >= n1 || a[1] >= n1 || a[2] >= n1) return 1;
  size_t k = ((size_t)a[2] * (size_t)n1 + (size_t)a[1]) * (size_t)n1 + (size_t)a[0];
  double tv = S->et[3 * k + (size_t)axis];
  if (!(tv >= 0.0)) return 1;
  *t = tv;
  for (int c = 0; c < 3; c++)
    nrm[c] = S->enrm[9 * k + 3 * (size_t)axis + (size_t)c];
  return 0;
}

/* --- 2. РАЗРЕЖЁННЫЙ ИСТОЧНИК (§370) ---------------------------------------- */

typedef struct {
  int32_t child0; /* -1 = лист */
  int32_t parent; /* -1 у корня; нужен ЗАЛИВКЕ для подъёма к общему предку (§373) */
  int32_t t0;     /* начало списка треугольников листа размера 1 */
  int32_t tn;     /* длина списка — по НЕрасширенной коробке, правило плотного пути */
  uint8_t hit;    /* коробка, РАСШИРЕННАЯ на ячейку, задета треугольником */
  uint8_t occ;    /* НЕрасширенная задета: ячейка ПЕРЕКРЫВАЕТ заливку */
  uint8_t out;    /* заливка дошла */
} snode;

/* ПАЛЕЦ (§373). Запомненный ПУТЬ последнего спуска. Спуск идёт в глубину,
 * поэтому следующий запрос почти всегда лежит под одним из запомненных узлов, и
 * начинать от КОРНЯ незачем: подъём по короткому пути дешевле девяти чтений
 * вразнобой. Результат от пальца НЕ ЗАВИСИТ — принадлежность проверяется всегда;
 * ускоряется только поиск.
 * ВНИМАНИЕ (А632): это изменяемое состояние в источнике, то есть источник
 * НЕ ПОТОКОБЕЗОПАСЕН. Когда спуск пойдёт в потоки (§366 п. 4), палец обязан
 * стать потоко-локальным, иначе гонка на записи пути. */
typedef struct {
  int32_t ni[HZ_DC_MAX_LOG2SIZE + 2];
  int32_t lo[HZ_DC_MAX_LOG2SIZE + 2][3];
  int32_t size[HZ_DC_MAX_LOG2SIZE + 2];
  int len;
} spath;

/* Глубочайший запомненный узел, НАКРЫВАЮЩИЙ запрос целиком. */
static int sp_find(const spath *p, const int32_t qlo[3], int32_t qsize) {
  for (int d = p->len - 1; d >= 0; d--) {
    int in = 1;
    for (int a = 0; a < 3; a++)
      if (qlo[a] < p->lo[d][a] || qlo[a] + qsize > p->lo[d][a] + p->size[d]) in = 0;
    if (in) return d;
  }
  return -1;
}

static void sp_put(spath *p, int d, int32_t ni, const int32_t lo[3], int32_t size) {
  p->ni[d] = ni;
  p->size[d] = size;
  for (int a = 0; a < 3; a++)
    p->lo[d][a] = lo[a];
}

typedef struct {
  snode *nd;
  int32_t n, cap;
  int32_t *tri;
  int32_t ntri, captri;
  frame fr;
  const hz_objmesh *m;
  double *tlo, *thi;                     /* габариты треугольников — дешёвый отсев */
  int noflood, invert, skip0, badfinger; /* ключи контролей */
  int64_t nsign, nbox, ncross;           /* цена спуска СЧИТАЕТСЯ, а не оценивается */
  spath fcell, fbox;                     /* пальцы: поиск ячейки и поиск коробки */
} sfield;

static int32_t s_alloc8(sfield *S, int32_t ni) {
  if (S->n + 8 > S->cap) {
    int32_t nc = S->cap * 2;
    if (nc < S->n + 8) nc = S->n + 8;
    snode *nn = realloc(S->nd, (size_t)nc * sizeof *nn);
    if (nn == NULL) return -1;
    S->nd = nn;
    S->cap = nc;
  }
  int32_t base = S->n;
  memset(&S->nd[base], 0, 8 * sizeof *S->nd);
  for (int i = 0; i < 8; i++) {
    S->nd[base + i].child0 = -1;
    S->nd[base + i].parent = ni;
  }
  S->n += 8;
  return base;
}

static int s_tri_push(sfield *S, int32_t t) {
  if (S->ntri >= S->captri) {
    int32_t nc = S->captri * 2;
    int32_t *nn = realloc(S->tri, (size_t)nc * sizeof *nn);
    if (nn == NULL) return -1;
    S->tri = nn;
    S->captri = nc;
  }
  S->tri[S->ntri++] = t;
  return 0;
}

/* Коробка МЕЛКОЙ ячейки — формулой ПЛОТНОГО пути дословно (`cl + h`, а не
 * `(x+1)·h`): сложение неассоциативно, и другая формула дала бы другой набор
 * принятых треугольников на самой границе, то есть сверка мерила бы разницу
 * формул, а не разницу спусков. */
static void s_cell_box(const sfield *S, const int32_t lo[3], double cl[3], double ch[3]) {
  for (int c = 0; c < 3; c++) {
    cl[c] = S->fr.org[c] + (double)lo[c] * S->fr.h;
    ch[c] = cl[c] + S->fr.h;
  }
}

/* Коробка узла, РАСШИРЕННАЯ на одну ячейку в каждую сторону (§370 п. 3). */
static void s_pad_box(const sfield *S, const int32_t lo[3], int32_t size, double cl[3],
                      double ch[3]) {
  for (int c = 0; c < 3; c++) {
    cl[c] = S->fr.org[c] + (double)(lo[c] - 1) * S->fr.h;
    ch[c] = S->fr.org[c] + (double)(lo[c] + size + 1) * S->fr.h;
  }
}

/* Дешёвый отсев по габаритам ПЕРЕД точным отсечением. Множество принятых не
 * меняется: точное отсечение всё равно требует пересечения габаритов, и отсев
 * нестрогий. */
static int aabb_hit(const sfield *S, int32_t t, const double cl[3], const double ch[3]) {
  for (int c = 0; c < 3; c++)
    if (S->thi[3 * (size_t)t + (size_t)c] < cl[c] || S->tlo[3 * (size_t)t + (size_t)c] > ch[c])
      return 0;
  return 1;
}

static void tri_verts(const hz_objmesh *m, int32_t t, const double **A, const double **B,
                      const double **C) {
  *A = m->v + 3 * (size_t)m->f[3 * (size_t)t + 0];
  *B = m->v + 3 * (size_t)m->f[3 * (size_t)t + 1];
  *C = m->v + 3 * (size_t)m->f[3 * (size_t)t + 2];
}

static int s_build(sfield *S, int32_t ni, const int32_t lo[3], int32_t size, const int32_t *cand,
                   int32_t nc) {
  S->nd[ni].child0 = -1;
  S->nd[ni].hit = 0;
  S->nd[ni].occ = 0;
  S->nd[ni].out = 0;
  S->nd[ni].t0 = 0;
  S->nd[ni].tn = 0;
  if (nc == 0) return 0; /* пусто даже с ореолом — смены знака внутри нет */
  S->nd[ni].hit = 1;

  if (size == 1) {
    double cl[3], ch[3];
    s_cell_box(S, lo, cl, ch);
    int32_t start = S->ntri, k = 0;
    for (int32_t i = 0; i < nc; i++) {
      const double *A, *B, *C;
      if (!aabb_hit(S, cand[i], cl, ch)) continue;
      tri_verts(S->m, cand[i], &A, &B, &C);
      hz_pclip_poly P;
      if (hz_pclip_tri(A, B, C, cl, ch, &P) < 3) continue;
      if (!(hz_pclip_area(&P) > 0.0)) continue;
      if (s_tri_push(S, cand[i]) != 0) return -1;
      k++;
    }
    S->nd[ni].t0 = start;
    S->nd[ni].tn = k;
    S->nd[ni].occ = (uint8_t)(k > 0);
    return 0;
  }

  int32_t c0 = s_alloc8(S, ni);
  if (c0 < 0) return -1;
  S->nd[ni].child0 = c0;
  int32_t half = size / 2;
  int32_t *sub = malloc((size_t)nc * sizeof *sub);
  if (sub == NULL) return -1;
  for (int i = 0; i < 8; i++) {
    int32_t clo[3];
    for (int a = 0; a < 3; a++)
      clo[a] = lo[a] + (((i >> a) & 1) ? half : 0);
    double cl[3], ch[3];
    s_pad_box(S, clo, half, cl, ch);
    int32_t k = 0;
    for (int32_t j = 0; j < nc; j++) {
      const double *A, *B, *C;
      if (!aabb_hit(S, cand[j], cl, ch)) continue;
      tri_verts(S->m, cand[j], &A, &B, &C);
      hz_pclip_poly P;
      /* Здесь площадь НЕ требуется: ореол обязан быть консервативным, а
       * касание по мере нуль всё равно означает геометрию рядом. */
      if (hz_pclip_tri(A, B, C, cl, ch, &P) < 3) continue;
      sub[k++] = cand[j];
    }
    if (s_build(S, c0 + i, clo, half, sub, k) != 0) {
      free(sub);
      return -1;
    }
  }
  free(sub);
  return 0;
}

/* --- заливка по листьям разрежённого дерева -------------------------------- */

typedef struct {
  int32_t ni, lo[3], size;
} sleaf;

typedef struct {
  sleaf *v;
  int64_t n, cap;
  int oom;
} sstack;

static void st_push(sstack *st, int32_t ni, const int32_t lo[3], int32_t size) {
  if (st->n >= st->cap) {
    int64_t nc = st->cap * 2;
    sleaf *nv = realloc(st->v, (size_t)nc * sizeof *nv);
    if (nv == NULL) {
      st->oom = 1;
      return;
    }
    st->v = nv;
    st->cap = nc;
  }
  st->v[st->n].ni = ni;
  st->v[st->n].size = size;
  for (int a = 0; a < 3; a++)
    st->v[st->n].lo[a] = lo[a];
  st->n++;
}

/* Пометить и положить в стек все листья, задетые коробкой [qlo, qhi). Это и есть
 * поиск соседей по грани: у листьев разные размеры, и таблицы соседства нет —
 * есть спуск по коробке. */
static void s_spread(sfield *S, sstack *st, int32_t ni, const int32_t lo[3], int32_t size,
                     const int32_t qlo[3], const int32_t qhi[3]) {
  for (int a = 0; a < 3; a++)
    if (lo[a] >= qhi[a] || lo[a] + size <= qlo[a]) return;
  if (S->nd[ni].child0 < 0) {
    if (S->nd[ni].occ || S->nd[ni].out) return;
    S->nd[ni].out = 1;
    st_push(st, ni, lo, size);
    return;
  }
  int32_t half = size / 2;
  for (int i = 0; i < 8; i++) {
    int32_t clo[3];
    for (int a = 0; a < 3; a++)
      clo[a] = lo[a] + (((i >> a) & 1) ? half : 0);
    s_spread(S, st, S->nd[ni].child0 + i, clo, half, qlo, qhi);
  }
}

static void s_seed(sfield *S, sstack *st, int32_t ni, const int32_t lo[3], int32_t size) {
  if (S->nd[ni].child0 >= 0) {
    int32_t half = size / 2;
    for (int i = 0; i < 8; i++) {
      int32_t clo[3];
      for (int a = 0; a < 3; a++)
        clo[a] = lo[a] + (((i >> a) & 1) ? half : 0);
      s_seed(S, st, S->nd[ni].child0 + i, clo, half);
    }
    return;
  }
  if (S->nd[ni].occ || S->nd[ni].out) return;
  int touch = 0;
  for (int a = 0; a < 3; a++)
    if (lo[a] == 0 || lo[a] + size == S->fr.n) touch = 1;
  if (!touch) return;
  S->nd[ni].out = 1;
  st_push(st, ni, lo, size);
}

/* Объём в МЕЛКИХ ЯЧЕЙКАХ по трём состояниям: снаружи, недостигнуто, с геометрией. */
static void s_volume(const sfield *S, int32_t ni, int32_t size, int64_t v[3]) {
  if (S->nd[ni].child0 >= 0) {
    int32_t half = size / 2;
    for (int i = 0; i < 8; i++)
      s_volume(S, S->nd[ni].child0 + i, half, v);
    return;
  }
  int64_t c = (int64_t)size * (int64_t)size * (int64_t)size;
  if (S->nd[ni].occ)
    v[2] += c;
  else if (S->nd[ni].out)
    v[0] += c;
  else
    v[1] += c;
}

static int s_flood(sfield *S) {
  sstack st;
  st.cap = 1024;
  st.n = 0;
  st.oom = 0;
  st.v = malloc((size_t)st.cap * sizeof *st.v);
  if (st.v == NULL) return -1;
  int32_t zero[3] = {0, 0, 0};
  s_seed(S, &st, 0, zero, S->fr.n);
  while (st.n > 0 && !st.oom) {
    sleaf L = st.v[--st.n];
    for (int a = 0; a < 3; a++)
      for (int d = 0; d < 2; d++) {
        int32_t qlo[3], qhi[3];
        for (int b = 0; b < 3; b++) {
          qlo[b] = L.lo[b];
          qhi[b] = L.lo[b] + L.size;
        }
        qlo[a] = d ? L.lo[a] + L.size : L.lo[a] - 1;
        qhi[a] = qlo[a] + 1;
        if (qlo[a] < 0 || qlo[a] >= S->fr.n) continue;
        /* ПОДЪЁМ ДО ОБЩЕГО ПРЕДКА, А НЕ СПУСК ОТ КОРНЯ (§373). Коробка предка
         * восстанавливается арифметикой, а не хранением: коробка выровнена по
         * своему размеру, поэтому у родителя это `lo & ~(2·size − 1)` при
         * размере `2·size`. В худшем случае (сосед по другую сторону середины
         * сцены) поднимемся до корня и не выиграем ничего — А636. */
        int32_t ani = L.ni, alo[3] = {L.lo[0], L.lo[1], L.lo[2]}, asize = L.size;
        for (;;) {
          int cov = 1;
          for (int b = 0; b < 3; b++)
            if (qlo[b] < alo[b] || qhi[b] > alo[b] + asize) cov = 0;
          if (cov || S->nd[ani].parent < 0) break;
          ani = S->nd[ani].parent;
          asize *= 2;
          for (int b = 0; b < 3; b++)
            alo[b] &= ~(asize - 1);
        }
        s_spread(S, &st, ani, alo, asize, qlo, qhi);
      }
  }
  int bad = st.oom;
  free(st.v);
  return bad ? -1 : 0;
}

/* --- обратные вызовы разрежённого источника -------------------------------- */

static const snode *s_leaf_at(sfield *S, const int32_t cell[3]) {
  spath *p = &S->fcell;
  /* НЕГАТИВНЫЙ КОНТРОЛЬ `badfinger` (§373): палец отдаёт прошлый лист, НЕ
   * проверив, что ячейка в нём лежит. Ровно та ошибка, которую палец и может
   * внести; приёмка обязана её увидеть, иначе принята вслепую. */
  if (S->badfinger && p->len > 0) return &S->nd[p->ni[p->len - 1]];
  int d = sp_find(p, cell, 1);
  int32_t ni, size, lo[3];
  if (d < 0) {
    d = 0;
    ni = 0;
    size = S->fr.n;
    lo[0] = lo[1] = lo[2] = 0;
    sp_put(p, 0, ni, lo, size);
  } else {
    ni = p->ni[d];
    size = p->size[d];
    for (int a = 0; a < 3; a++)
      lo[a] = p->lo[d][a];
  }
  while (S->nd[ni].child0 >= 0) {
    int32_t half = size / 2;
    int bit = 0;
    for (int a = 0; a < 3; a++)
      if (cell[a] >= lo[a] + half) {
        bit |= 1 << a;
        lo[a] += half;
      }
    ni = S->nd[ni].child0 + bit;
    size = half;
    sp_put(p, ++d, ni, lo, size);
  }
  p->len = d + 1;
  return &S->nd[ni];
}

/* Знак угла: СНАРУЖИ, если хотя бы одна из восьми смежных ячеек достигнута
 * заливкой либо лежит вне сетки. Правило дословно то же, что у плотного пути. */
static int s_sign(void *ctx, const int32_t c[3]) {
  sfield *S = (sfield *)ctx;
  int32_t n = S->fr.n;
  S->nsign++;
  if (c[0] < 0 || c[1] < 0 || c[2] < 0 || c[0] > n || c[1] > n || c[2] > n) return 0;
  int isout = 0;
  for (int dz = -1; dz <= 0 && !isout; dz++)
    for (int dy = -1; dy <= 0 && !isout; dy++)
      for (int dx = -1; dx <= 0 && !isout; dx++) {
        int32_t cc[3] = {c[0] + dx, c[1] + dy, c[2] + dz};
        if (cc[0] < 0 || cc[1] < 0 || cc[2] < 0 || cc[0] >= n || cc[1] >= n || cc[2] >= n) {
          isout = 1;
          break;
        }
        if (s_leaf_at(S, cc)->out) isout = 1;
      }
  int sv = S->noflood ? 0 : (isout ? 0 : 1);
  if (S->invert) sv = S->noflood ? 0 : !sv;
  return sv;
}

static int s_cross(void *ctx, const int32_t a[3], int axis, double *t, double nrm[3]) {
  sfield *S = (sfield *)ctx;
  int32_t n = S->fr.n;
  S->ncross++;
  if (a[0] < 0 || a[1] < 0 || a[2] < 0 || a[0] > n || a[1] > n || a[2] > n) return 1;
  double P0[3], P1[3];
  for (int c = 0; c < 3; c++) {
    P0[c] = S->fr.org[c] + (double)a[c] * S->fr.h;
    P1[c] = P0[c];
  }
  P1[axis] += S->fr.h;
  double best = 2.0, bn[3] = {0.0, 0.0, 0.0};
  /* Кандидаты — из ячеек, к которым ребро принадлежит (до четырёх); порядок
   * обхода тот же, что у плотного пути: при равенстве t побеждает первый. */
  for (int dp = -1; dp <= 0; dp++)
    for (int dq = -1; dq <= 0; dq++) {
      int32_t c3[3] = {a[0], a[1], a[2]};
      c3[(axis + 1) % 3] += dp;
      c3[(axis + 2) % 3] += dq;
      if (c3[0] < 0 || c3[1] < 0 || c3[2] < 0 || c3[0] >= n || c3[1] >= n || c3[2] >= n) continue;
      const snode *L = s_leaf_at(S, c3);
      for (int32_t q = L->t0; q < L->t0 + L->tn; q++) {
        const double *A, *B, *C;
        tri_verts(S->m, S->tri[q], &A, &B, &C);
        double e1[3], e2[3], pv[3], tv[3], qv[3], dir[3];
        for (int c = 0; c < 3; c++) {
          e1[c] = B[c] - A[c];
          e2[c] = C[c] - A[c];
          dir[c] = P1[c] - P0[c];
        }
        pv[0] = dir[1] * e2[2] - dir[2] * e2[1];
        pv[1] = dir[2] * e2[0] - dir[0] * e2[2];
        pv[2] = dir[0] * e2[1] - dir[1] * e2[0];
        double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
        if (det > -1e-15 && det < 1e-15) continue;
        double inv = 1.0 / det;
        for (int c = 0; c < 3; c++)
          tv[c] = P0[c] - A[c];
        double uu = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
        if (uu < 0.0 || uu > 1.0) continue;
        qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
        qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
        qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
        double vv = (dir[0] * qv[0] + dir[1] * qv[1] + dir[2] * qv[2]) * inv;
        if (vv < 0.0 || uu + vv > 1.0) continue;
        double tt = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
        if (tt < 0.0 || tt > 1.0) continue;
        if (tt < best) {
          best = tt;
          double cr[3];
          cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
          cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
          cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
          double l2 = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
          if (!(l2 > 0.0)) l2 = 1.0;
          for (int c = 0; c < 3; c++)
            bn[c] = cr[c] / l2;
        }
      }
    }
  if (best <= 1.0) {
    *t = best;
    for (int c = 0; c < 3; c++)
      nrm[c] = bn[c];
    return 0;
  }
  /* ЗАПЛАТА — то же правило, что в сведении плотного пути (§369, А611): знак
   * сменился, а треугольника на ребре нет. Молчать нельзя, потому что
   * построитель требует «пересечение ⟺ смена знака». */
  *t = 0.5;
  for (int c = 0; c < 3; c++)
    nrm[c] = 0.0;
  nrm[axis] = 1.0;
  return 0;
}

static int s_boxq(void *ctx, const int32_t lo[3], int32_t size) {
  sfield *S = (sfield *)ctx;
  S->nbox++;
  /* НЕГАТИВНЫЙ КОНТРОЛЬ `skip0` (§370): предикат ЛЖЁТ «пусто» на одном октанте
   * корня. Отказ обязан быть fail-open — геометрия исчезает молча. */
  if (S->skip0 && size == S->fr.n / 2 && lo[0] == 0 && lo[1] == 0 && lo[2] == 0) return 0;
  spath *p = &S->fbox;
  int d = sp_find(p, lo, size);
  int32_t ni, sz, l[3];
  if (d < 0) {
    d = 0;
    ni = 0;
    sz = S->fr.n;
    l[0] = l[1] = l[2] = 0;
    sp_put(p, 0, ni, l, sz);
  } else {
    ni = p->ni[d];
    sz = p->size[d];
    for (int a = 0; a < 3; a++)
      l[a] = p->lo[d][a];
  }
  while (sz > size) {
    if (!S->nd[ni].hit) {
      p->len = d + 1;
      return 0;
    }
    if (S->nd[ni].child0 < 0) { /* лист крупнее запроса: отвечаем «может быть» */
      p->len = d + 1;
      return 1;
    }
    int32_t half = sz / 2;
    int bit = 0;
    for (int a = 0; a < 3; a++)
      if (lo[a] >= l[a] + half) {
        bit |= 1 << a;
        l[a] += half;
      }
    ni = S->nd[ni].child0 + bit;
    sz = half;
    sp_put(p, ++d, ni, l, sz);
  }
  p->len = d + 1;
  return S->nd[ni].hit ? 1 : 0;
}

/* --- 3. ОШИБКА ПОВЕРХНОСТИ И СВЕРКА ДЕРЕВЬЕВ ------------------------------- */

/* Список треугольников ячейки — у двух источников он берётся по-разному, но
 * набор кандидатов обязан быть один и тот же. */
typedef int32_t (*trilist_fn)(void *ctx, const int32_t cell[3], const int32_t **out);

typedef struct {
  const int32_t *cnt, *lst;
  int32_t n;
} dense_lists;

static int32_t dense_tris(void *ctx, const int32_t cell[3], const int32_t **out) {
  const dense_lists *D = (const dense_lists *)ctx;
  size_t ci = ((size_t)cell[2] * (size_t)D->n + (size_t)cell[1]) * (size_t)D->n + (size_t)cell[0];
  *out = D->lst + D->cnt[ci];
  return D->cnt[ci + 1] - D->cnt[ci];
}

static int32_t sparse_tris(void *ctx, const int32_t cell[3], const int32_t **out) {
  sfield *S = (sfield *)ctx;
  const snode *L = s_leaf_at(S, cell);
  *out = S->tri + L->t0;
  return L->tn;
}

typedef struct {
  int32_t lo[3], size;
  double vx[3];
} vrec;

static void collect_verts(const hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size,
                          vrec *out, int64_t *n) {
  if (t->nd[ni].flags & HZ_DC_HASVERT) {
    if (out != NULL) {
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
    collect_verts(t, t->nd[ni].child0 + i, clo, half, out, n);
  }
}

static int vrec_cmp(const void *x, const void *y) {
  const vrec *a = (const vrec *)x, *b = (const vrec *)y;
  for (int k = 0; k < 3; k++)
    if (a->lo[k] != b->lo[k]) return a->lo[k] < b->lo[k] ? -1 : 1;
  if (a->size != b->size) return a->size < b->size ? -1 : 1;
  return 0;
}

static int same_bits(double a, double b) {
  uint64_t x, y;
  memcpy(&x, &a, sizeof x);
  memcpy(&y, &b, sizeof y);
  return x == y;
}

static void surf_err(const hz_dctree *T, const frame *fr, const hz_objmesh *m, trilist_fn tl,
                     void *tctx, const char *tag) {
  int64_t nv = 0;
  collect_verts(T, 0, (int32_t[3]){0, 0, 0}, fr->n, NULL, &nv);
  if (nv == 0) return;
  double *er = malloc((size_t)nv * sizeof *er);
  vrec *vr = malloc((size_t)nv * sizeof *vr);
  if (er == NULL || vr == NULL) {
    free(er);
    free(vr);
    return;
  }
  int64_t k = 0;
  collect_verts(T, 0, (int32_t[3]){0, 0, 0}, fr->n, vr, &k);
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

/* --- 4. главная ------------------------------------------------------------ */

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr,
            "pfield ФАЙЛ.obj МАСШТАБ [lev=N] [both] [noflood] [invert] [skip0] [badfinger]\n");
    return 2;
  }
  int lev = 6, noflood = 0, invert = 0, skip0 = 0, both = 0, badfinger = 0;
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0) lev = (int)strtol(argv[i] + 4, NULL, 10);
    /* ПЛОТНЫЙ ЭТАЛОН РЯДОМ (Г42): оба пути в одном прогоне, сверка побитовая. */
    if (strcmp(argv[i], "both") == 0) both = 1;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ (§365): знак не строится вовсе — вершин обязано выйти 0. */
    if (strcmp(argv[i], "noflood") == 0) noflood = 1;
    /* СИЛЬНЫЙ КОНТРОЛЬ (А600): знак обращён — поверхность обязана выйти ТОЙ ЖЕ. */
    if (strcmp(argv[i], "invert") == 0) invert = 1;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ (§370): предикат лжёт «пусто» на октанте корня. */
    if (strcmp(argv[i], "skip0") == 0) skip0 = 1;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ (§373): палец без проверки принадлежности. */
    if (strcmp(argv[i], "badfinger") == 0) badfinger = 1;
  }
  if (lev < 1 || lev > HZ_DC_MAX_LOG2SIZE) {
    fprintf(stderr, "lev вне разрядного предела (1..%d)\n", HZ_DC_MAX_LOG2SIZE);
    return 2;
  }
  if (both && lev > HZ_DC_SAMPLE_MAX_LOG2SIZE) {
    fprintf(stderr, "both требует lev <= %d: плотного эталона выше не существует\n",
            HZ_DC_SAMPLE_MAX_LOG2SIZE);
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
  printf("== ПОЛЕ: %s, треугольников %d, сетка %d^3, ячейка %.4f м%s%s%s%s\n", argv[1], m.nt, fr.n,
         fr.h, noflood ? "  [НЕТ ЗАЛИВКИ]" : "", invert ? "  [ЗНАК ОБРАЩЁН]" : "",
         skip0 ? "  [ПРЕДИКАТ ЛЖЁТ НА ОКТАНТЕ 0]" : "", badfinger ? "  [ПАЛЕЦ БЕЗ ПРОВЕРКИ]" : "");

  int64_t dense_nv = 0;
  vrec *dense_v = NULL;

  /* ================= ПЛОТНЫЙ ПУТЬ (эталон) ================= */
  if (both) {
    double t0 = now_s();
    int32_t n = fr.n, n1 = n + 1;
    size_t ncell = (size_t)n * (size_t)n * (size_t)n;
    /* На нехватке памяти стенд УМИРАЕТ, а не возвращает код: восстанавливать
     * нечего, а «освободить и вернуть» плодит ветви, которые никогда не
     * исполняются и потому не проверены ни одним прогоном. */
    int32_t *cnt = calloc(ncell + 1, sizeof *cnt);
    unsigned char *occ = calloc(ncell, 1);
    if (cnt == NULL || occ == NULL) exit(1);
    for (int32_t t = 0; t < m.nt; t++) {
      int64_t i0[3] = {0, 0, 0}, i1[3] = {0, 0, 0};
      for (int c = 0; c < 3; c++) {
        double a = 1e300, b = -1e300;
        for (int v = 0; v < 3; v++) {
          double x = m.v[3 * (size_t)m.f[3 * (size_t)t + (size_t)v] + (size_t)c];
          if (x < a) a = x;
          if (x > b) b = x;
        }
        i0[c] = (int64_t)floor((a - fr.org[c]) / fr.h);
        i1[c] = (int64_t)floor((b - fr.org[c]) / fr.h);
        if (i0[c] < 0) i0[c] = 0;
        if (i1[c] >= n) i1[c] = n - 1;
      }
      const double *A, *B, *C;
      tri_verts(&m, t, &A, &B, &C);
      for (int64_t z = i0[2]; z <= i1[2]; z++)
        for (int64_t y = i0[1]; y <= i1[1]; y++)
          for (int64_t x = i0[0]; x <= i1[0]; x++) {
            double cl[3] = {fr.org[0] + (double)x * fr.h, fr.org[1] + (double)y * fr.h,
                            fr.org[2] + (double)z * fr.h};
            double ch[3] = {cl[0] + fr.h, cl[1] + fr.h, cl[2] + fr.h};
            hz_pclip_poly P;
            if (hz_pclip_tri(A, B, C, cl, ch, &P) < 3) continue;
            if (!(hz_pclip_area(&P) > 0.0)) continue;
            size_t ci = ((size_t)z * (size_t)n + (size_t)y) * (size_t)n + (size_t)x;
            cnt[ci + 1]++;
            occ[ci] = 1;
          }
    }
    for (size_t i = 0; i < ncell; i++)
      cnt[i + 1] += cnt[i];
    int32_t *lst = malloc((size_t)cnt[ncell] * sizeof *lst);
    int32_t *fill = calloc(ncell, sizeof *fill);
    if (lst == NULL || fill == NULL) exit(1);
    for (int32_t t = 0; t < m.nt; t++) {
      int64_t i0[3] = {0, 0, 0}, i1[3] = {0, 0, 0};
      for (int c = 0; c < 3; c++) {
        double a = 1e300, b = -1e300;
        for (int v = 0; v < 3; v++) {
          double x = m.v[3 * (size_t)m.f[3 * (size_t)t + (size_t)v] + (size_t)c];
          if (x < a) a = x;
          if (x > b) b = x;
        }
        i0[c] = (int64_t)floor((a - fr.org[c]) / fr.h);
        i1[c] = (int64_t)floor((b - fr.org[c]) / fr.h);
        if (i0[c] < 0) i0[c] = 0;
        if (i1[c] >= n) i1[c] = n - 1;
      }
      const double *A, *B, *C;
      tri_verts(&m, t, &A, &B, &C);
      for (int64_t z = i0[2]; z <= i1[2]; z++)
        for (int64_t y = i0[1]; y <= i1[1]; y++)
          for (int64_t x = i0[0]; x <= i1[0]; x++) {
            double cl[3] = {fr.org[0] + (double)x * fr.h, fr.org[1] + (double)y * fr.h,
                            fr.org[2] + (double)z * fr.h};
            double ch[3] = {cl[0] + fr.h, cl[1] + fr.h, cl[2] + fr.h};
            hz_pclip_poly P;
            if (hz_pclip_tri(A, B, C, cl, ch, &P) < 3) continue;
            if (!(hz_pclip_area(&P) > 0.0)) continue;
            size_t ci = ((size_t)z * (size_t)n + (size_t)y) * (size_t)n + (size_t)x;
            lst[cnt[ci] + fill[ci]++] = t;
          }
    }
    int64_t nocc = 0;
    for (size_t i = 0; i < ncell; i++)
      if (occ[i]) nocc++;
    printf("   [плотный] списки по ячейкам за %.2f с: занятых %lld (%.2f %%), вхождений %d\n",
           now_s() - t0, (long long)nocc, 100.0 * (double)nocc / (double)ncell, cnt[ncell]);
    free(fill);

    t0 = now_s();
    unsigned char *out = calloc(ncell, 1);
    int32_t *stk = malloc(ncell * sizeof *stk);
    if (out == NULL || stk == NULL) exit(1);
    int64_t sp = 0;
    for (int32_t z = 0; z < n; z++)
      for (int32_t y = 0; y < n; y++)
        for (int32_t x = 0; x < n; x++) {
          if (x != 0 && y != 0 && z != 0 && x != n - 1 && y != n - 1 && z != n - 1) continue;
          size_t ci = ((size_t)z * (size_t)n + (size_t)y) * (size_t)n + (size_t)x;
          if (occ[ci] || out[ci]) continue;
          out[ci] = 1;
          stk[sp++] = (int32_t)ci;
        }
    while (sp > 0) {
      int32_t ci = stk[--sp];
      int32_t x = ci % n, y = (ci / n) % n, z = ci / (n * n);
      for (int d = 0; d < 6; d++) {
        int32_t nx = x + ((d == 0) - (d == 1)), ny = y + ((d == 2) - (d == 3)),
                nz = z + ((d == 4) - (d == 5));
        if (nx < 0 || ny < 0 || nz < 0 || nx >= n || ny >= n || nz >= n) continue;
        size_t k2 = ((size_t)nz * (size_t)n + (size_t)ny) * (size_t)n + (size_t)nx;
        if (occ[k2] || out[k2]) continue;
        out[k2] = 1;
        stk[sp++] = (int32_t)k2;
      }
    }
    int64_t nout = 0, nin = 0;
    for (size_t i = 0; i < ncell; i++) {
      if (out[i])
        nout++;
      else if (!occ[i])
        nin++;
    }
    printf("   [плотный] заливка за %.2f с: снаружи %lld (%.1f %%), ВНУТРИ %lld (%.1f %%)\n",
           now_s() - t0, (long long)nout, 100.0 * (double)nout / (double)ncell, (long long)nin,
           100.0 * (double)nin / (double)ncell);
    free(stk);

    src_ctx S;
    S.n = n;
    S.sgn = calloc((size_t)n1 * (size_t)n1 * (size_t)n1, 1);
    if (S.sgn == NULL) exit(1);
    int64_t nsin = 0;
    for (int32_t z = 0; z <= n; z++)
      for (int32_t y = 0; y <= n; y++)
        for (int32_t x = 0; x <= n; x++) {
          int isout = 0;
          for (int dz = -1; dz <= 0 && !isout; dz++)
            for (int dy = -1; dy <= 0 && !isout; dy++)
              for (int dx = -1; dx <= 0 && !isout; dx++) {
                int32_t cx = x + dx, cy = y + dy, cz = z + dz;
                if (cx < 0 || cy < 0 || cz < 0 || cx >= n || cy >= n || cz >= n) {
                  isout = 1;
                  break;
                }
                if (out[((size_t)cz * (size_t)n + (size_t)cy) * (size_t)n + (size_t)cx]) isout = 1;
              }
          int sv = noflood ? 0 : (isout ? 0 : 1);
          if (invert) sv = noflood ? 0 : !sv;
          S.sgn[((size_t)z * (size_t)n1 + (size_t)y) * (size_t)n1 + (size_t)x] = (unsigned char)sv;
          if (sv) nsin++;
        }
    printf("   [плотный] углов ВНУТРИ %lld из %lld (%.1f %%)\n", (long long)nsin,
           (long long)((size_t)n1 * (size_t)n1 * (size_t)n1),
           100.0 * (double)nsin / (double)((size_t)n1 * (size_t)n1 * (size_t)n1));

    t0 = now_s();
    size_t nc1 = (size_t)n1 * (size_t)n1 * (size_t)n1;
    S.et = malloc(3 * nc1 * sizeof *S.et);
    S.enrm = malloc(9 * nc1 * sizeof *S.enrm);
    if (S.et == NULL || S.enrm == NULL) exit(1);
    for (size_t i = 0; i < 3 * nc1; i++)
      S.et[i] = -1.0;
    int64_t ncross = 0;
    for (int32_t z = 0; z <= n; z++)
      for (int32_t y = 0; y <= n; y++)
        for (int32_t x = 0; x <= n; x++) {
          size_t k = ((size_t)z * (size_t)n1 + (size_t)y) * (size_t)n1 + (size_t)x;
          for (int ax = 0; ax < 3; ax++) {
            int32_t e[3] = {x, y, z};
            if (e[ax] >= n) continue;
            double P0[3] = {fr.org[0] + (double)x * fr.h, fr.org[1] + (double)y * fr.h,
                            fr.org[2] + (double)z * fr.h};
            double P1[3] = {P0[0], P0[1], P0[2]};
            P1[ax] += fr.h;
            double best = 2.0, bn[3] = {0.0, 0.0, 0.0};
            for (int dp = -1; dp <= 0; dp++)
              for (int dq = -1; dq <= 0; dq++) {
                int32_t c3[3] = {x, y, z};
                c3[(ax + 1) % 3] += dp;
                c3[(ax + 2) % 3] += dq;
                if (c3[0] < 0 || c3[1] < 0 || c3[2] < 0 || c3[0] >= n || c3[1] >= n || c3[2] >= n)
                  continue;
                size_t ci = ((size_t)c3[2] * (size_t)n + (size_t)c3[1]) * (size_t)n + (size_t)c3[0];
                for (int32_t q = cnt[ci]; q < cnt[ci + 1]; q++) {
                  const double *A, *B, *C;
                  tri_verts(&m, lst[q], &A, &B, &C);
                  double e1[3], e2[3], pv[3], tv[3], qv[3], dir[3];
                  for (int c = 0; c < 3; c++) {
                    e1[c] = B[c] - A[c];
                    e2[c] = C[c] - A[c];
                    dir[c] = P1[c] - P0[c];
                  }
                  pv[0] = dir[1] * e2[2] - dir[2] * e2[1];
                  pv[1] = dir[2] * e2[0] - dir[0] * e2[2];
                  pv[2] = dir[0] * e2[1] - dir[1] * e2[0];
                  double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
                  if (det > -1e-15 && det < 1e-15) continue;
                  double inv = 1.0 / det;
                  for (int c = 0; c < 3; c++)
                    tv[c] = P0[c] - A[c];
                  double uu = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
                  if (uu < 0.0 || uu > 1.0) continue;
                  qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
                  qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
                  qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
                  double vv = (dir[0] * qv[0] + dir[1] * qv[1] + dir[2] * qv[2]) * inv;
                  if (vv < 0.0 || uu + vv > 1.0) continue;
                  double tt = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
                  if (tt < 0.0 || tt > 1.0) continue;
                  if (tt < best) {
                    best = tt;
                    double cr[3];
                    cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
                    cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
                    cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
                    double l2 = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
                    if (!(l2 > 0.0)) l2 = 1.0;
                    for (int c = 0; c < 3; c++)
                      bn[c] = cr[c] / l2;
                  }
                }
              }
            if (best <= 1.0) {
              S.et[3 * k + (size_t)ax] = best;
              for (int c = 0; c < 3; c++)
                S.enrm[9 * k + 3 * (size_t)ax + (size_t)c] = bn[c];
              ncross++;
            }
          }
        }
    printf("   [плотный] эрмитовы данные за %.2f с: рёбер с пересечением %lld\n", now_s() - t0,
           (long long)ncross);

    /* Сведение знака и пересечений (§369, А611): построитель требует
     * «пересечение ⟺ смена знака», а источника два. */
    int64_t nfix = 0, ndrop = 0, nkeep = 0;
    for (int32_t z = 0; z <= n; z++)
      for (int32_t y = 0; y <= n; y++)
        for (int32_t x = 0; x <= n; x++) {
          size_t k = ((size_t)z * (size_t)n1 + (size_t)y) * (size_t)n1 + (size_t)x;
          for (int ax = 0; ax < 3; ax++) {
            int32_t e[3] = {x, y, z};
            if (e[ax] >= n) continue;
            int32_t e2i[3] = {x, y, z};
            e2i[ax]++;
            size_t k2 =
                ((size_t)e2i[2] * (size_t)n1 + (size_t)e2i[1]) * (size_t)n1 + (size_t)e2i[0];
            int s0 = S.sgn[k], s1 = S.sgn[k2];
            int has = (S.et[3 * k + (size_t)ax] >= 0.0);
            if (s0 == s1) {
              if (has) {
                S.et[3 * k + (size_t)ax] = -1.0;
                ndrop++;
              }
            } else if (!has) {
              S.et[3 * k + (size_t)ax] = 0.5;
              for (int c = 0; c < 3; c++)
                S.enrm[9 * k + 3 * (size_t)ax + (size_t)c] = 0.0;
              S.enrm[9 * k + 3 * (size_t)ax + (size_t)ax] = 1.0;
              nfix++;
            } else
              nkeep++;
          }
        }
    printf("   [плотный] СВЕДЕНИЕ: пересечений при смене знака %lld; заплат %lld; отброшено %lld\n",
           (long long)nkeep, (long long)nfix, (long long)ndrop);

    t0 = now_s();
    hz_signgrid g;
    hz_htab ht;
    memset(&g, 0, sizeof g);
    if (hz_htab_init(&ht) != 0) exit(1);
    int rc = hz_dc_sample(&g, &ht, lev, src_sign, src_cross, &S);
    if (rc != HZ_DC_OK) {
      fprintf(stderr, "плотный опрос источника: код %d\n", rc);
      return 1;
    }
    hz_dctree T;
    if (hz_dc_init(&T, lev) != 0) exit(1);
    rc = hz_dc_build(&T, &g, &ht);
    printf("   [плотный] дерево DC за %.2f с: код %d, узлов %d, рёбер %d; ЗАГНАНО %d, "
           "НЕМАНИФОЛДНЫХ %d\n",
           now_s() - t0, rc, T.n, ht.n, T.nclamped, T.nmulti);
    collect_verts(&T, 0, (int32_t[3]){0, 0, 0}, fr.n, NULL, &dense_nv);
    printf("   [плотный] вершин выдано %lld\n", (long long)dense_nv);
    if (dense_nv > 0) {
      dense_v = malloc((size_t)dense_nv * sizeof *dense_v);
      if (dense_v == NULL) exit(1);
      int64_t kk = 0;
      collect_verts(&T, 0, (int32_t[3]){0, 0, 0}, fr.n, dense_v, &kk);
      qsort(dense_v, (size_t)dense_nv, sizeof *dense_v, vrec_cmp);
      /* V0 — вершины, чья коробка лежит в октанте [0, n/2)³. Это ТО САМОЕ
       * предсказанное число негативного контроля `skip0` (§370), и берётся оно
       * из НЕЗАВИСИМОГО источника — плотного дерева. */
      int64_t v0 = 0;
      for (int64_t i = 0; i < dense_nv; i++) {
        int in = 1;
        for (int a = 0; a < 3; a++)
          if (dense_v[i].lo[a] + dense_v[i].size > fr.n / 2) in = 0;
        v0 += in;
      }
      printf("   [плотный] V0 (вершин в октанте [0,n/2)³) %lld  → skip0 обязан дать %lld\n",
             (long long)v0, (long long)(dense_nv - v0));
    }
    dense_lists DL = {cnt, lst, n};
    surf_err(&T, &fr, &m, dense_tris, &DL, "[плотный]");

    hz_dc_free(&T);
    hz_htab_free(&ht);
    hz_signgrid_free(&g);
    free(S.sgn);
    free(S.et);
    free(S.enrm);
    free(out);
    free(occ);
    free(cnt);
    free(lst);
  }

  /* ================= РАЗРЕЖЁННЫЙ ПУТЬ ================= */
  sfield S;
  memset(&S, 0, sizeof S);
  S.fr = fr;
  S.m = &m;
  S.noflood = noflood;
  S.invert = invert;
  S.skip0 = skip0;
  S.badfinger = badfinger;
  S.cap = 64;
  S.nd = calloc((size_t)S.cap, sizeof *S.nd);
  S.captri = 1024;
  S.tri = malloc((size_t)S.captri * sizeof *S.tri);
  S.tlo = malloc(3 * (size_t)m.nt * sizeof *S.tlo);
  S.thi = malloc(3 * (size_t)m.nt * sizeof *S.thi);
  int32_t *root = malloc((size_t)m.nt * sizeof *root);
  if (S.nd == NULL || S.tri == NULL || S.tlo == NULL || S.thi == NULL || root == NULL) exit(1);
  S.n = 1;
  S.nd[0].child0 = -1;
  for (int32_t t = 0; t < m.nt; t++) {
    const double *A, *B, *C;
    tri_verts(&m, t, &A, &B, &C);
    for (int c = 0; c < 3; c++) {
      double a = A[c] < B[c] ? A[c] : B[c];
      double b = A[c] > B[c] ? A[c] : B[c];
      if (C[c] < a) a = C[c];
      if (C[c] > b) b = C[c];
      S.tlo[3 * (size_t)t + (size_t)c] = a;
      S.thi[3 * (size_t)t + (size_t)c] = b;
    }
    root[t] = t;
  }

  double t0 = now_s();
  int32_t zero[3] = {0, 0, 0};
  if (s_build(&S, 0, zero, fr.n, root, m.nt) != 0) {
    fprintf(stderr, "нет памяти на дерево занятости\n");
    exit(1);
  }
  free(root);
  double t_tree = now_s() - t0;
  int64_t nleaf = 0, nocc = 0;
  for (int32_t i = 0; i < S.n; i++) {
    if (S.nd[i].child0 >= 0) continue;
    nleaf++;
    if (S.nd[i].occ) nocc++;
  }
  printf("   дерево занятости за %.2f с: узлов %d, листьев %lld, ЗАНЯТЫХ мелких %lld, "
         "вхождений %d\n",
         t_tree, S.n, (long long)nleaf, (long long)nocc, S.ntri);

  t0 = now_s();
  if (s_flood(&S) != 0) {
    fprintf(stderr, "нет памяти на заливку\n");
    exit(1);
  }
  double t_flood = now_s() - t0;
  /* ОБЪЁМ, А НЕ ЧИСЛО ЛИСТЬЕВ: листья разного размера, и «листьев снаружи 40 %»
   * ничего не значит. Доля сравнима с плотной строкой «ВНУТРИ 10.4 %».
   * УГЛЫ не считаются нарочно: их (n+1)³, и перечисление вернуло бы ровно ту
   * цену O(8^L), ради снятия которой шаг и делается (П4 меряется по ячейкам). */
  int64_t vol[3] = {0, 0, 0}; /* снаружи, ВНУТРИ, с геометрией */
  s_volume(&S, 0, fr.n, vol);
  double vtot = (double)vol[0] + (double)vol[1] + (double)vol[2];
  printf("   заливка за %.2f с: снаружи %.1f %%, ВНУТРИ %.1f %%, с геометрией %.2f %% (ячеек)\n",
         t_flood, 100.0 * (double)vol[0] / vtot, 100.0 * (double)vol[1] / vtot,
         100.0 * (double)vol[2] / vtot);
  {
    double eye[3] = HZ_CFG_HALL_EYE;
    int32_t ic[3] = {0, 0, 0};
    int okc = 1;
    for (int c = 0; c < 3; c++) {
      double f = floor((eye[c] - fr.org[c]) / fr.h);
      if (!(f >= 0.0) || !(f < (double)fr.n)) okc = 0;
      ic[c] = okc ? (int32_t)f : 0;
    }
    if (okc) {
      const snode *L = s_leaf_at(&S, ic);
      printf("   ЯЧЕЙКА КАМЕРЫ §2: %s\n",
             L->occ ? "С ГЕОМЕТРИЕЙ"
                    : (L->out ? "снаружи (заливка дошла)" : "ВНУТРИ — ЗАЛИВКА НЕ ДОШЛА"));
    } else
      printf("   ЯЧЕЙКА КАМЕРЫ §2: вне сетки\n");
  }

  hz_dctree T;
  if (hz_dc_init(&T, lev) != 0) exit(1);
  t0 = now_s();
  int rc = hz_dc_shape_lazy(&T, lev, s_sign, s_boxq, &S);
  double t_shape = now_s() - t0;
  if (rc != HZ_DC_OK) {
    fprintf(stderr, "спуск за структурой: код %d\n", rc);
    return 1;
  }
  printf("   спуск за структурой за %.2f с: узлов %d\n", t_shape, T.n);

  hz_htab ht;
  if (hz_htab_init(&ht) != 0) exit(1);
  t0 = now_s();
  int32_t ndup = 0;
  rc = hz_dc_edges_lazy(&ht, &T, s_cross, &S, &ndup);
  double t_edges = now_s() - t0;
  if (rc != HZ_DC_OK) {
    fprintf(stderr, "проход за рёбрами: код %d\n", rc);
    return 1;
  }
  printf("   проход за рёбрами за %.2f с: рёбер %d, повторов схлопнуто %d\n", t_edges, ht.n, ndup);

  t0 = now_s();
  rc = hz_dc_forms_lazy(&T, &ht);
  double t_dc = now_s() - t0;
  printf("   формы и вершины за %.2f с: код %d, узлов %d; ЗАГНАНО %d, НЕМАНИФОЛДНЫХ %d\n", t_dc, rc,
         T.n, T.nclamped, T.nmulti);
  int64_t nv = 0;
  collect_verts(&T, 0, zero, fr.n, NULL, &nv);
  printf("   вершин выдано %lld\n", (long long)nv);
  printf("   ВЫЗОВОВ ИСТОЧНИКА: sign %lld, box %lld, cross %lld\n", (long long)S.nsign,
         (long long)S.nbox, (long long)S.ncross);
  printf("   ПАМЯТЬ: дерево DC %.1f МБ (%zu Б/узел), дерево занятости %.1f МБ, списки %.1f МБ, "
         "рёбра %.1f МБ\n",
         (double)T.n * (double)sizeof(hz_dcnode) / 1048576.0, sizeof(hz_dcnode),
         (double)S.n * (double)sizeof(snode) / 1048576.0, (double)S.ntri * 4.0 / 1048576.0,
         (double)ht.n * (double)sizeof(hz_hedge) / 1048576.0);
  surf_err(&T, &fr, &m, sparse_tris, &S, "");

  /* ---- КАРТИНКА ПОЛЯ: обход -> полигоны -> z-буфер -> img/ ---- */
  {
    tr3_camera cam;
    double eye[3] = HZ_CFG_HALL_EYE, at[3] = HZ_CFG_HALL_AT, up[3] = HZ_CFG_UP;
    /* 512² — разрешение, названное в бюджете кадра (§367) */
    if (tr3_camera_look(&cam, eye, at, up, HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0, 512,
                        512) == 0) {
      size_t np = 512u * 512u;
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
      snprintf(path, sizeof path, "img/pfield_hall_L%d.ppm", lev);
      int prc = hz_ppm_write(path, vb, 512, 512);
      printf("   КАРТИНКА за %.2f с: обход код %d, треугольников %lld (за камерой %lld), "
             "заполнено %lld из %zu пикселей (%.1f %%) -> %s (код %d)\n",
             t_img, wrc, (long long)R.ntri, (long long)R.nclip, (long long)nfill, np,
             100.0 * (double)nfill / (double)np, path, prc);
      free(zb);
      free(vb);
    }
  }

  /* ---- сверка двух путей: по КОРОБКЕ узла, побитово (А613) ---- */
  if (both && dense_v != NULL) {
    vrec *sv = malloc((size_t)(nv > 0 ? nv : 1) * sizeof *sv);
    if (sv == NULL) exit(1);
    int64_t k = 0;
    collect_verts(&T, 0, zero, fr.n, sv, &k);
    qsort(sv, (size_t)k, sizeof *sv, vrec_cmp);
    int64_t i = 0, j = 0, same = 0, diff = 0, onlyd = 0, onlys = 0;
    while (i < dense_nv && j < k) {
      int c = vrec_cmp(&dense_v[i], &sv[j]);
      if (c < 0) {
        onlyd++;
        i++;
      } else if (c > 0) {
        onlys++;
        j++;
      } else {
        int eq = 1;
        for (int a = 0; a < 3; a++)
          if (!same_bits(dense_v[i].vx[a], sv[j].vx[a])) eq = 0;
        same += eq;
        diff += !eq;
        i++;
        j++;
      }
    }
    onlyd += dense_nv - i;
    onlys += k - j;
    printf("   СВЕРКА С ПЛОТНЫМ: совпало ПОБИТОВО %lld из %lld; в той же коробке, но иначе %lld; "
           "только у плотного %lld; только у ленивого %lld\n",
           (long long)same, (long long)dense_nv, (long long)diff, (long long)onlyd,
           (long long)onlys);
    free(sv);
  }

  hz_dc_free(&T);
  hz_htab_free(&ht);
  free(dense_v);
  free(S.nd);
  free(S.tri);
  free(S.tlo);
  free(S.thi);
  hz_obj_free(&m);
  return 0;
}
