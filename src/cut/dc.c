/* PLAN_CUT.md Р-5б и Р-5в. Эрмитовы данные, QEF-иерархия, обход и мост к
 * разрезу. Побитовые обещания здесь те же, что в ядре (poly3.c), и держатся на
 * том же: одна формула, посчитанная один раз, и -ffp-contract=off в CFLAGS. */

#include "cut/dc.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* --- 1. эрмитовы рёбра ----------------------------------------------------- */

int hz_htab_init(hz_htab *h) {
  h->cap = 256;
  h->n = 0;
  h->e = calloc((size_t)h->cap, sizeof(hz_hedge));
  return h->e == NULL ? HZ_DC_ENOMEM : HZ_DC_OK;
}

void hz_htab_free(hz_htab *h) {
  free(h->e);
  h->e = NULL;
  h->n = h->cap = 0;
}

uint64_t hz_hedge_key(int axis, const int32_t p[3]) {
  return ((uint64_t)(uint32_t)axis << 60) | ((uint64_t)(uint32_t)p[0] << 40) |
         ((uint64_t)(uint32_t)p[1] << 20) | (uint64_t)(uint32_t)p[2];
}

int hz_htab_add(hz_htab *h, int axis, const int32_t p[3], double t, const double nrm[3],
                int in_lo) {
  if (h->n >= h->cap) {
    int32_t nc = h->cap * 2;
    hz_hedge *ne = realloc(h->e, (size_t)nc * sizeof(hz_hedge));
    if (ne == NULL) return HZ_DC_ENOMEM;
    h->e = ne;
    h->cap = nc;
  }
  hz_hedge *e = &h->e[h->n];
  e->key = hz_hedge_key(axis, p);
  for (int k = 0; k < 3; k++) {
    e->p[k] = p[k];
    e->nrm[k] = nrm[k];
  }
  e->axis = axis;
  e->t = t;
  e->in_lo = in_lo ? 1u : 0u;
  h->n++;
  return HZ_DC_OK;
}

static int hedge_cmp(const void *a, const void *b) {
  const hz_hedge *x = a, *y = b;
  if (x->key < y->key) return -1;
  return x->key > y->key ? 1 : 0;
}

int hz_htab_sort(hz_htab *h) {
  if (h->n > 1) qsort(h->e, (size_t)h->n, sizeof(hz_hedge), hedge_cmp);
  for (int32_t i = 1; i < h->n; i++)
    if (h->e[i].key == h->e[i - 1].key) return HZ_DC_EDUP;
  return HZ_DC_OK;
}

/* Побитовое равенство двух чисел. Именно битов, а не значений: сверяется
 * ОДИНАКОВОСТЬ ВЫЧИСЛЕНИЯ у четырёх ячеек, и «почти равно» тут не годится
 * (заодно -Wfloat-equal запрещает писать ==). */
static int same_bits(double a, double b) {
  uint64_t x, y;
  memcpy(&x, &a, sizeof x);
  memcpy(&y, &b, sizeof y);
  return x == y;
}

static int hedge_same(const hz_hedge *a, const hz_hedge *b) {
  if (a->axis != b->axis || a->in_lo != b->in_lo) return 0;
  for (int k = 0; k < 3; k++)
    if (a->p[k] != b->p[k] || !same_bits(a->nrm[k], b->nrm[k])) return 0;
  return same_bits(a->t, b->t);
}

int hz_htab_uniq(hz_htab *h, int32_t *ndup) {
  int32_t dup = 0;
  if (h->n > 1) qsort(h->e, (size_t)h->n, sizeof(hz_hedge), hedge_cmp);
  int32_t w = 0;
  for (int32_t i = 0; i < h->n; i++) {
    if (w > 0 && h->e[w - 1].key == h->e[i].key) {
      if (!hedge_same(&h->e[w - 1], &h->e[i])) return HZ_DC_EDUP;
      dup++;
      continue;
    }
    h->e[w++] = h->e[i];
  }
  h->n = w;
  if (ndup != NULL) *ndup = dup;
  return HZ_DC_OK;
}

/* СВЕДЕНИЕ ПО БЛИЖАЙШЕМУ ПЕРЕСЕЧЕНИЮ (Р7, §393). Беззнаковый источник заносит
 * ребро ОДИН РАЗ НА КАЖДЫЙ пересекающий его треугольник, и записи законно
 * разные — `hz_htab_uniq` здесь не годится, он требует побитового совпадения.
 *
 * ПРАВИЛО ВЫБОРА ОДНО И НАЗВАНО: остаётся пересечение с НАИМЕНЬШИМ `t`, то есть
 * ПЕРВОЕ ОТ НИЖНЕГО КОНЦА ребра. Ровно это делал и знаковый источник
 * (`best = min t`), так что правило не новое. ЧЕГО ОНО СТОИТ: если ребро
 * пересекают два листа поверхности, второй теряется — но одна вершина на ячейку
 * его и не представила бы (принятое ограничение DC), так что потеря не
 * добавляется этим правилом, а лишь не лечится им.
 *
 * ПОРЯДОК СТРОГИЙ И НЕ ЗАВИСИТ ОТ ПОРЯДКА ЗАНЕСЕНИЯ: при равных `t` сравниваются
 * координаты нормали. Иначе итог зависел бы от того, в каком порядке шли
 * треугольники, и первое же распараллеливание источника его бы изменило. */
static int hedge_cmp_mint(const void *x, const void *y) {
  const hz_hedge *a = (const hz_hedge *)x, *b = (const hz_hedge *)y;
  /* Сравнения строго через `<`/`>`, а не через `!=`: -Wfloat-equal запрещает
   * второе, и запрещает по делу — здесь нужен ПОРЯДОК, а не равенство. */
  if (a->key != b->key) return a->key < b->key ? -1 : 1;
  if (a->t < b->t) return -1;
  if (a->t > b->t) return 1;
  for (int c = 0; c < 3; c++) {
    if (a->nrm[c] < b->nrm[c]) return -1;
    if (a->nrm[c] > b->nrm[c]) return 1;
  }
  return 0;
}

int hz_htab_reduce_min(hz_htab *h, int32_t *ndrop) {
  int32_t drop = 0;
  if (h->n > 1) qsort(h->e, (size_t)h->n, sizeof(hz_hedge), hedge_cmp_mint);
  int32_t w = 0;
  for (int32_t i = 0; i < h->n; i++) {
    if (w > 0 && h->e[w - 1].key == h->e[i].key) {
      drop++;
      continue;
    }
    h->e[w++] = h->e[i];
  }
  h->n = w;
  if (ndrop != NULL) *ndrop = drop;
  return HZ_DC_OK;
}

const hz_hedge *hz_htab_find(const hz_htab *h, int axis, const int32_t p[3]) {
  uint64_t k = hz_hedge_key(axis, p);
  int32_t lo = 0, hi = h->n - 1;
  while (lo <= hi) {
    int32_t mid = lo + (hi - lo) / 2;
    if (h->e[mid].key == k) return &h->e[mid];
    if (h->e[mid].key < k)
      lo = mid + 1;
    else
      hi = mid - 1;
  }
  return NULL;
}

void hz_hedge_point(const hz_hedge *e, double p[3]) {
  for (int k = 0; k < 3; k++)
    p[k] = (double)e->p[k];
  p[e->axis] += e->t;
}

/* --- 2. сетка знаков ------------------------------------------------------- */

void hz_signgrid_free(hz_signgrid *g) {
  free(g->s);
  g->s = NULL;
  g->n1 = 0;
}

int hz_signgrid_at(const hz_signgrid *g, const int32_t c[3]) {
  return g->s[((size_t)c[0] * (size_t)g->n1 + (size_t)c[1]) * (size_t)g->n1 + (size_t)c[2]];
}

/* ОДИН опрос ребра со сменой знака. Вынесен из hz_dc_sample НЕ ради краткости:
 * ленивый спуск обязан заносить в таблицу ПОБИТОВО ту же запись, иначе сверка
 * двух путей меряет разницу двух копий формулы, а не разницу спусков. */
static int probe_edge(hz_htab *ht, hz_dc_cross cr, void *ctx, const int32_t p[3], int axis,
                      int s0) {
  double t = 0.0, nrm[3] = {0, 0, 0};
  if (cr(ctx, p, axis, &t, nrm) != 0) return HZ_DC_ETOPO;
  if (!(t >= 0.0) || !(t <= 1.0)) return HZ_DC_ETOPO;
  double m = sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
  if (!(m > 0.0)) return HZ_DC_ETOPO;
  /* Нормировка ЗДЕСЬ и один раз: вес образца в QEF есть |n|², поэтому
   * неединичная нормаль молча перевесила бы соседей. Источник вправе вернуть
   * градиент. */
  for (int k = 0; k < 3; k++)
    nrm[k] /= m;
  return hz_htab_add(ht, axis, p, t, nrm, s0);
}

int hz_dc_sample(hz_signgrid *g, hz_htab *ht, int log2size, hz_dc_sign sg, hz_dc_cross cr,
                 void *ctx) {
  if (log2size < 0 || log2size > HZ_DC_SAMPLE_MAX_LOG2SIZE) return HZ_DC_ERANGE;
  int32_t n = (int32_t)1 << log2size, n1 = n + 1;
  size_t cells = (size_t)n1 * (size_t)n1 * (size_t)n1;
  /* аллокация ПО МЕСТУ, без обёрток (CLAUDE.md: gcc-analyzer и обёртки) */
  g->s = calloc(cells, 1);
  if (g->s == NULL) return HZ_DC_ENOMEM;
  g->n1 = n1;
  g->log2size = log2size;

  for (int32_t x = 0; x < n1; x++)
    for (int32_t y = 0; y < n1; y++)
      for (int32_t z = 0; z < n1; z++) {
        int32_t c[3] = {x, y, z};
        g->s[((size_t)x * (size_t)n1 + (size_t)y) * (size_t)n1 + (size_t)z] = sg(ctx, c) ? 1u : 0u;
      }

  /* Ребро заносится РОВНО ОДИН раз: обход по ключу, а не по ячейкам. Именно это
   * и делает точку общей для всех четырёх ячеек (Г3 на уровне данных). */
  for (int axis = 0; axis < 3; axis++) {
    int32_t hi[3] = {n1, n1, n1};
    hi[axis] = n;
    for (int32_t x = 0; x < hi[0]; x++)
      for (int32_t y = 0; y < hi[1]; y++)
        for (int32_t z = 0; z < hi[2]; z++) {
          int32_t p[3] = {x, y, z}, q[3] = {x, y, z};
          q[axis]++;
          int s0 = hz_signgrid_at(g, p), s1 = hz_signgrid_at(g, q);
          if (s0 == s1) continue;
          int rc = probe_edge(ht, cr, ctx, p, axis, s0);
          if (rc != HZ_DC_OK) return rc;
        }
  }
  return hz_htab_sort(ht);
}

/* --- 3. дерево ------------------------------------------------------------- */

int hz_dc_manifold(uint8_t corner) {
  if (corner == 0u || corner == 0xFFu) return 1;
  for (int side = 0; side < 2; side++) {
    uint8_t m = side ? (uint8_t)~corner : corner;
    int seed = -1;
    for (int c = 0; c < 8; c++)
      if (m & (1u << c)) {
        seed = c;
        break;
      }
    if (seed < 0) return 0; /* пустая сторона при corner != 0, 0xFF невозможна */
    uint8_t seen = (uint8_t)(1u << seed);
    for (int iter = 0; iter < 8; iter++) /* 8 узлов: 8 проходов заведомо хватает */
      for (int c = 0; c < 8; c++) {
        if (!(seen & (1u << c))) continue;
        for (int b = 0; b < 3; b++) {
          int d = c ^ (1 << b); /* соседи куба — по одному биту */
          if (m & (1u << d)) seen = (uint8_t)(seen | (1u << d));
        }
      }
    if (seen != m) return 0;
  }
  return 1;
}

int hz_dc_init(hz_dctree *t, int log2size) {
  if (log2size < 0 || log2size > HZ_DC_MAX_LOG2SIZE) return HZ_DC_ERANGE;
  t->log2size = log2size;
  t->cap = 64;
  t->n = 0;
  t->nclamped = 0;
  t->nmulti = 0;
  t->nd = calloc((size_t)t->cap, sizeof(hz_dcnode));
  if (t->nd == NULL) return HZ_DC_ENOMEM;
  t->nd[0].child0 = -1;
  hz_qef_zero(&t->nd[0].q);
  t->n = 1;
  return HZ_DC_OK;
}

void hz_dc_free(hz_dctree *t) {
  free(t->nd);
  t->nd = NULL;
  t->n = t->cap = 0;
}

/* Блок из 8 подряд, строго вперёд — дословно node_alloc8 октодерева: на «индекс
 * ребёнка больше индекса родителя» опирается и завершаемость спуска. */
static int32_t dc_alloc8(hz_dctree *t) {
  if (t->n + 8 > t->cap) {
    int32_t nc = t->cap * 2;
    if (nc < t->n + 8) nc = t->n + 8;
    hz_dcnode *nn = realloc(t->nd, (size_t)nc * sizeof(hz_dcnode));
    if (nn == NULL) return -1;
    t->nd = nn;
    t->cap = nc;
  }
  int32_t base = t->n;
  for (int i = 0; i < 8; i++) {
    memset(&t->nd[base + i], 0, sizeof(hz_dcnode));
    t->nd[base + i].child0 = -1;
    hz_qef_zero(&t->nd[base + i].q);
  }
  t->n += 8;
  return base;
}

static uint8_t corner_mask(const hz_signgrid *g, const int32_t lo[3], int32_t size) {
  uint8_t m = 0;
  for (int c = 0; c < 8; c++) {
    int32_t p[3];
    for (int a = 0; a < 3; a++)
      p[a] = lo[a] + ((c >> a) & 1 ? size : 0);
    if (hz_signgrid_at(g, p)) m = (uint8_t)(m | (1u << c));
  }
  return m;
}

/* Есть ли в коробке смена знака. Точно, а не эвристикой: сканируется вся
 * подсетка. Это и есть та цена ПО ВХОДУ, о которой предупреждает заголовок. */
static int box_uniform(const hz_signgrid *g, const int32_t lo[3], int32_t size) {
  int32_t p0[3] = {lo[0], lo[1], lo[2]};
  int first = hz_signgrid_at(g, p0);
  for (int32_t x = 0; x <= size; x++)
    for (int32_t y = 0; y <= size; y++)
      for (int32_t z = 0; z <= size; z++) {
        int32_t p[3] = {lo[0] + x, lo[1] + y, lo[2] + z};
        if (hz_signgrid_at(g, p) != first) return 0;
      }
  return 1;
}

/* 12 рёбер единичной ячейки: ось и смещение нижнего конца по двум другим осям. */
static void unit_edge(int i, int *axis, int off[3]) {
  int a = i / 4, k = i % 4;
  int u = (a + 1) % 3, v = (a + 2) % 3;
  off[0] = off[1] = off[2] = 0;
  off[u] = k & 1;
  off[v] = (k >> 1) & 1;
  *axis = a;
}

/* Форма листа размера 1 из таблицы рёбер. ПОРЯДОК ОБХОДА РЁБЕР ФИКСИРОВАН — на
 * нём стоит побитовое совпадение починки с перестройкой (Г49). */
static void leaf_qef(hz_dctree *t, int32_t ni, const int32_t lo[3], const hz_htab *ht) {
  hz_qef_zero(&t->nd[ni].q);
  for (int i = 0; i < 12; i++) {
    int axis, off[3];
    unit_edge(i, &axis, off);
    int32_t p[3] = {lo[0] + off[0], lo[1] + off[1], lo[2] + off[2]};
    const hz_hedge *e = hz_htab_find(ht, axis, p);
    if (e == NULL) continue;
    /* Локальная координата считается ЦЕЛОЧИСЛЕННОЙ разностью плюс t, а не
     * вычитанием двух больших чисел: сокращения нет вовсе (Г41). */
    double loc[3];
    for (int k = 0; k < 3; k++)
      loc[k] = (double)(e->p[k] - lo[k]) + (k == axis ? e->t : 0.0);
    hz_qef_add_sample(&t->nd[ni].q, loc, e->nrm);
  }
}

static void solve_node(hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size) {
  hz_dcnode *nd = &t->nd[ni];
  if (nd->flags & HZ_DC_CLAMPED) t->nclamped--;
  if (nd->flags & HZ_DC_MULTI) t->nmulti--;
  nd->flags = 0;
  nd->err = 0.0;
  for (int k = 0; k < 3; k++)
    nd->vx[k] = 0.0;
  if (nd->q.n <= 0) return;
  if (!hz_dc_manifold(nd->corner)) {
    /* Г47: вершины НЕ выдаём. Одна вершина на два листа поверхности — это
     * материал там, где его нет, а «коробка ∩ полуплоскости» двух листов не
     * несёт вовсе. Отказ явный и считанный. */
    nd->flags |= HZ_DC_MULTI;
    t->nmulti++;
    return;
  }
  double blo[3] = {0.0, 0.0, 0.0}, bhi[3] = {(double)size, (double)size, (double)size};
  double x[3], r;
  int st = hz_qef_solve(&nd->q, blo, bhi, x, &r);
  if (st == HZ_QEF_EMPTY) return;
  if (st == HZ_QEF_CLAMPED) {
    nd->flags |= HZ_DC_CLAMPED;
    t->nclamped++;
  }
  for (int k = 0; k < 3; k++)
    nd->vx[k] = (double)lo[k] + x[k];
  nd->err = r;
  nd->flags |= HZ_DC_HASVERT;
}

/* Форма родителя = сумма СДВИНУТЫХ форм детей. Порядок детей 0..7 фиксирован —
 * см. Г49: побитовой является ПОЧИНКА против ПЕРЕСТРОЙКИ, и держится она
 * именно на одинаковом порядке. */
static void sum_children(hz_dctree *t, int32_t ni, int32_t size) {
  int32_t c0 = t->nd[ni].child0;
  int32_t half = size / 2;
  hz_qef acc;
  hz_qef_zero(&acc);
  for (int i = 0; i < 8; i++) {
    double sh[3];
    for (int a = 0; a < 3; a++)
      sh[a] = ((i >> a) & 1) ? (double)half : 0.0;
    hz_qef_add_shifted(&acc, &t->nd[c0 + i].q, sh);
  }
  t->nd[ni].q = acc;
}

static int build_rec(hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size,
                     const hz_signgrid *g, const hz_htab *ht) {
  t->nd[ni].corner = corner_mask(g, lo, size);
  t->nd[ni].child0 = -1;
  hz_qef_zero(&t->nd[ni].q);
  if (box_uniform(g, lo, size)) return HZ_DC_OK; /* поверхности внутри нет */

  if (size == 1) {
    leaf_qef(t, ni, lo, ht);
    solve_node(t, ni, lo, size);
    return HZ_DC_OK;
  }

  int32_t c0 = dc_alloc8(t);
  if (c0 < 0) return HZ_DC_ENOMEM;
  t->nd[ni].child0 = c0;
  int32_t half = size / 2;
  for (int i = 0; i < 8; i++) {
    int32_t clo[3];
    for (int a = 0; a < 3; a++)
      clo[a] = lo[a] + (((i >> a) & 1) ? half : 0);
    int rc = build_rec(t, c0 + i, clo, half, g, ht);
    if (rc != HZ_DC_OK) return rc;
  }
  sum_children(t, ni, size);
  solve_node(t, ni, lo, size);
  return HZ_DC_OK;
}

int hz_dc_build(hz_dctree *t, const hz_signgrid *g, const hz_htab *ht) {
  if (g->log2size != t->log2size) return HZ_DC_ERANGE;
  int32_t zero[3] = {0, 0, 0};
  int rc = build_rec(t, 0, zero, (int32_t)1 << t->log2size, g, ht);
  if (rc != HZ_DC_OK) return rc;
  return t->nmulti > 0 ? HZ_DC_EMULTI : HZ_DC_OK;
}

/* Проход ТРЕТИЙ: формы и вершины, снизу вверх, БЕЗ ЕДИНОГО обращения к
 * источнику. Порядок действий и порядок детей те же, что у плотного пути
 * (build_rec), и это не экономия строк: разойдись они, побитовое совпадение
 * вершин перестало бы что-либо значить. Остановленный узел формы не получает —
 * ровно как в плотном пути, где box_uniform возвращал управление до solve_node. */
static void lazy_forms(hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size,
                       const hz_htab *ht) {
  if (t->nd[ni].child0 < 0) {
    if (size != 1) return;
    leaf_qef(t, ni, lo, ht);
    solve_node(t, ni, lo, size);
    return;
  }
  int32_t half = size / 2;
  for (int i = 0; i < 8; i++) {
    int32_t clo[3];
    for (int a = 0; a < 3; a++)
      clo[a] = lo[a] + (((i >> a) & 1) ? half : 0);
    lazy_forms(t, t->nd[ni].child0 + i, clo, half, ht);
  }
  sum_children(t, ni, size);
  solve_node(t, ni, lo, size);
}

/* --- 3б. БЕЗЗНАКОВЫЙ ПУТЬ (Р7, §393) --------------------------------------- */

/* ЧЕМ ЭТО ОТЛИЧАЕТСЯ ОТ СНЕСЁННОГО АДАПТЕРА, И ПОЧЕМУ РАЗЛИЧИЕ НЕ СЛОВЕСНОЕ.
 * Прежний ленивый спуск задавал источнику ТРИ вопроса — знак угла, пересечение
 * ребра, «есть ли смена знака в коробке», — и первые два суть вопросы ПРО ПОЛЕ,
 * то есть про то, что источник обязан не отвечать, а СТРОИТЬ (Р1, А650).
 * Здесь остаётся один вопрос, и он про КОРОБКУ: есть ли в ней геометрия. Рёбра
 * источник ПЕРЕДАЁТ таблицей, знака не существует вовсе.
 *
 * ЗАМЕР, ИЗ КОТОРОГО ЭТО ВЫРОСЛО (§391), А НЕ ВКУС: у сцены из односторонних
 * треугольников «внутри» не определено, и заливка давала либо потерю предметов
 * (не протекла), либо потерю стен (протекла) — уровня без потери не было ни
 * одного. */
static int shape_occ_rec(hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size, hz_dc_occ oc,
                         void *ctx) {
  t->nd[ni].child0 = -1;
  t->nd[ni].corner = 0; /* знака нет; hz_dc_manifold(0) = 1, отказа Г47 не будет */
  t->nd[ni].ecross = 0;
  t->nd[ni].edir = 0;
  hz_qef_zero(&t->nd[ni].q);
  if (!oc(ctx, lo, size) || size == 1) return HZ_DC_OK;
  int32_t c0 = dc_alloc8(t);
  if (c0 < 0) return HZ_DC_ENOMEM;
  t->nd[ni].child0 = c0;
  int32_t half = size / 2;
  for (int i = 0; i < 8; i++) {
    int32_t clo[3];
    for (int a = 0; a < 3; a++)
      clo[a] = lo[a] + (((i >> a) & 1) ? half : 0);
    int rc = shape_occ_rec(t, c0 + i, clo, half, oc, ctx);
    if (rc != HZ_DC_OK) return rc;
  }
  return HZ_DC_OK;
}

int hz_dc_shape_occ(hz_dctree *t, int log2size, hz_dc_occ oc, void *ctx) {
  if (log2size < 0 || log2size > HZ_DC_MAX_LOG2SIZE) return HZ_DC_ERANGE;
  if (t->log2size != log2size) return HZ_DC_ERANGE;
  t->unsgn = 1;
  t->nbigmask = 0;
  int32_t zero[3] = {0, 0, 0};
  return shape_occ_rec(t, 0, zero, (int32_t)1 << log2size, oc, ctx);
}

/* Маски рёбер. У листа размера 1 — прямым поиском в таблице; порядок рёбер тот
 * же `unit_edge`, что у формы листа, и это не удобство: на одинаковом порядке
 * стоит побитовость починки (Г49).
 *
 * У ВНУТРЕННЕГО УЗЛА ребро `i` накрыто РОВНО ДВУМЯ детьми — теми, у кого
 * совпали биты по двум поперечным осям; их собственные рёбра с тем же номером
 * `i` и есть его половины. Поэтому `ecross` родителя есть ИЛИ по этим двум, и
 * это тождество, а не приближение. `edir` берётся у ПЕРВОГО из двух, у кого
 * пересечение есть, и при двух пересечениях сразу он НЕ ОПРЕДЕЛЁН — долг Ш3
 * (А708). При полной глубине обход читает маску только у самой мелкой из
 * четырёх ячеек (Г43), то есть всегда у ребра сетки, и долг не наступает. */
static void masks_occ_rec(hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size,
                          const hz_htab *ht) {
  hz_dcnode *nd = &t->nd[ni];
  if (nd->child0 < 0) {
    nd->ecross = 0;
    nd->edir = 0;
    if (size != 1) return;
    for (int i = 0; i < 12; i++) {
      int axis, off[3];
      unit_edge(i, &axis, off);
      int32_t p[3] = {lo[0] + off[0], lo[1] + off[1], lo[2] + off[2]};
      const hz_hedge *e = hz_htab_find(ht, axis, p);
      if (e == NULL) continue;
      nd->ecross = (uint16_t)(nd->ecross | (1u << i));
      if (e->nrm[axis] > 0.0) nd->edir = (uint16_t)(nd->edir | (1u << i));
    }
    return;
  }
  int32_t c0 = nd->child0, half = size / 2;
  for (int k = 0; k < 8; k++) {
    int32_t clo[3];
    for (int a = 0; a < 3; a++)
      clo[a] = lo[a] + (((k >> a) & 1) ? half : 0);
    masks_occ_rec(t, c0 + k, clo, half, ht);
  }
  uint16_t ec = 0, ed = 0;
  for (int i = 0; i < 12; i++) {
    int axis, off[3];
    unit_edge(i, &axis, off);
    int u = (axis + 1) % 3, v = (axis + 2) % 3;
    int base = (off[u] << u) | (off[v] << v);
    for (int h = 0; h < 2; h++) {
      const hz_dcnode *ch = &t->nd[c0 + (base | (h << axis))];
      if (!((ch->ecross >> i) & 1u)) continue;
      if (!((ec >> i) & 1u)) {
        ec = (uint16_t)(ec | (1u << i));
        if ((ch->edir >> i) & 1u) ed = (uint16_t)(ed | (1u << i));
      }
    }
  }
  t->nd[ni].ecross = ec;
  t->nd[ni].edir = ed;
}

/* А709: у недробившегося узла пересечённых рёбер быть не может — спуск идёт по
 * занятости, а ребро, задетое треугольником, лежит в задетой коробке. Это не
 * соглашение, а следствие предиката, и оно СЧИТАЕТСЯ, а не предполагается. */
static void masks_check_rec(hz_dctree *t, int32_t ni, int32_t size) {
  if (t->nd[ni].child0 < 0) {
    if (size != 1 && t->nd[ni].ecross != 0) t->nbigmask++;
    return;
  }
  for (int k = 0; k < 8; k++)
    masks_check_rec(t, t->nd[ni].child0 + k, size / 2);
}

int hz_dc_masks_occ(hz_dctree *t, const hz_htab *ht) {
  int32_t zero[3] = {0, 0, 0};
  masks_occ_rec(t, 0, zero, (int32_t)1 << t->log2size, ht);
  t->nbigmask = 0;
  masks_check_rec(t, 0, (int32_t)1 << t->log2size);
  return HZ_DC_OK;
}

int hz_dc_forms_lazy(hz_dctree *t, const hz_htab *ht) {
  int32_t zero[3] = {0, 0, 0};
  lazy_forms(t, 0, zero, (int32_t)1 << t->log2size, ht);
  return t->nmulti > 0 ? HZ_DC_EMULTI : HZ_DC_OK;
}

int hz_dc_repair(hz_dctree *t, const hz_htab *ht, const int32_t cell[3]) {
  int32_t path[HZ_DC_MAX_LOG2SIZE + 1], psize[HZ_DC_MAX_LOG2SIZE + 1];
  int32_t plo[HZ_DC_MAX_LOG2SIZE + 1][3];
  int depth = 0;
  int32_t ni = 0, size = (int32_t)1 << t->log2size, lo[3] = {0, 0, 0};
  for (int a = 0; a < 3; a++)
    if (cell[a] < 0 || cell[a] >= size) return HZ_DC_ERANGE;

  while (size > 1 && t->nd[ni].child0 >= 0) {
    path[depth] = ni;
    psize[depth] = size;
    for (int a = 0; a < 3; a++)
      plo[depth][a] = lo[a];
    depth++;
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
  /* Лист крупнее ячейки — значит правка меняет СЕТКУ ЗНАКОВ, а с ней форму
   * дерева. Это перестройка поддерева, а не починка, и молчать нельзя. */
  if (size != 1) return HZ_DC_ETOPO;

  leaf_qef(t, ni, lo, ht);
  solve_node(t, ni, lo, size);
  for (int d = depth - 1; d >= 0; d--) {
    sum_children(t, path[d], psize[d]);
    solve_node(t, path[d], plo[d], psize[d]);
  }
  return HZ_DC_OK;
}

/* --- 4. обход (Р-5в) ------------------------------------------------------- */

/* ГЕОМЕТРИЯ ВМЕСТО ТАБЛИЦ. Классический адаптивный DC выбирает детей у грани и у
 * ребра по печатным таблицам индексов. Здесь каждый ref несёт свою коробку,
 * поэтому нужный ребёнок находится ПРОБНОЙ ТОЧКОЙ в удвоенных целых координатах:
 * ошибиться в таблице нельзя, потому что таблицы нет. Цена — по одному сравнению
 * на ось, и она же снимает вопрос «а верна ли таблица», на который иначе
 * отвечать нечем (Г42: для этого и держится второй, независимый обход). */

typedef struct {
  const hz_dctree *t;
  hz_dc_stop stop;
  void *sctx;
  hz_dc_poly emit;
  void *pctx;
  int rc;
  int32_t nskip; /* пропущено полигонов у неманифолдных ячеек — СЧИТАЕТСЯ (§377) */
} walkctx;

static int leafish(const walkctx *w, const hz_dcref *r) {
  if (w->t->nd[r->ni].child0 < 0) return 1;
  return w->stop != NULL && w->stop(w->sctx, w->t, r);
}

static void child_at2(const hz_dctree *t, const hz_dcref *r, const int32_t pt2[3], hz_dcref *o) {
  int32_t half = r->size / 2;
  int bit = 0;
  o->size = half;
  for (int a = 0; a < 3; a++) {
    int hi = pt2[a] >= 2 * r->lo[a] + r->size;
    if (hi) bit |= 1 << a;
    o->lo[a] = r->lo[a] + (hi ? half : 0);
  }
  o->ni = t->nd[r->ni].child0 + bit;
}

/* Содержит ли коробка ref точку с удвоенными координатами (строго внутри). */
static int ref_holds2(const hz_dcref *r, const int32_t p2[3]) {
  for (int a = 0; a < 3; a++)
    if (p2[a] <= 2 * r->lo[a] || p2[a] >= 2 * (r->lo[a] + r->size)) return 0;
  return 1;
}

/* Четыре ячейки вокруг ребра — В ПОРЯДКЕ ПРОТИВ ЧАСОВОЙ вокруг направления e.
 * Порядок не назначается, а ВЫЧИСЛЯЕТСЯ: для каждого из четырёх квадрантов
 * ищется ячейка, его накрывающая. Так же работает и случай, когда одна крупная
 * ячейка накрывает ДВА квадранта — а он штатный на стыке уровней. */
static const int8_t QUAD_UV[4][2] = {{1, 1}, {-1, 1}, {-1, -1}, {1, -1}};

static int order_quad(const hz_dcref *in, int nin, int e, const int32_t q2[3], hz_dcref *out) {
  int u = (e + 1) % 3, v = (e + 2) % 3;
  for (int k = 0; k < 4; k++) {
    int32_t p2[3];
    for (int a = 0; a < 3; a++)
      p2[a] = q2[a];
    p2[e] += 1;
    p2[u] += QUAD_UV[k][0];
    p2[v] += QUAD_UV[k][1];
    int found = -1;
    for (int j = 0; j < nin; j++)
      if (ref_holds2(&in[j], p2)) {
        found = j;
        break;
      }
    if (found < 0) return HZ_DC_ETOPO;
    out[k] = in[found];
  }
  return HZ_DC_OK;
}

static void edge_proc(walkctx *w, const hz_dcref q[4], int e, const int32_t qlo[3], int32_t seg);

/* Минимальное ребро достигнуто: все четыре ячейки — листья среза. */
static void process_edge(walkctx *w, const hz_dcref q[4], int e, const int32_t qlo[3],
                         int32_t seg) {
  int u = (e + 1) % 3, v = (e + 2) % 3;
  int mi = 0;
  for (int k = 1; k < 4; k++)
    if (q[k].size < q[mi].size) mi = k;
  /* Инвариант обхода: сегмент И ЕСТЬ ребро самой мелкой из четырёх. Если он
   * нарушен, знак брать неоткуда — молчать здесь нельзя. */
  if (q[mi].size != seg || q[mi].lo[e] != qlo[e]) {
    w->rc = HZ_DC_ETOPO;
    return;
  }
  /* Г43: ЗНАК БЕРЁТСЯ С САМОЙ МЕЛКОЙ. У крупного узла маска углов унаследована
   * от детей, и про середину своего ребра он не знает ничего. */
  int c0 = 0;
  if (qlo[u] != q[mi].lo[u]) c0 |= 1 << u;
  if (qlo[v] != q[mi].lo[v]) c0 |= 1 << v;
  int s0;
  if (w->t->unsgn) {
    /* БЕЗЗНАКОВЫЙ КЛЮЧ (Р7, §393): поверхность опознаётся ПЕРЕСЕЧЁННЫМ РЕБРОМ,
     * а обход многоугольника — НОРМАЛЬЮ образца на нём. Номер ребра в нумерации
     * `unit_edge`: ось `e`, смещения по двум поперечным берутся из тех же битов
     * `c0`, что и знак у знакового пути, — то есть место читается одинаково, а
     * значение по-разному. `s0 = 1` значит «наружу смотрит +e», ровно как
     * «материал у нижнего конца» у знакового. */
    int i = e * 4 + (((c0 >> u) & 1) | (((c0 >> v) & 1) << 1));
    if (!((w->t->nd[q[mi].ni].ecross >> i) & 1u)) return;
    s0 = (w->t->nd[q[mi].ni].edir >> i) & 1u;
  } else {
    int c1 = c0 | (1 << e);
    s0 = (w->t->nd[q[mi].ni].corner >> c0) & 1;
    int s1 = (w->t->nd[q[mi].ni].corner >> c1) & 1;
    if (s0 == s1) return;
  }

  double vv[4][3];
  hz_dcref rr[4];
  int nv = 0;
  for (int k = 0; k < 4; k++) {
    const hz_dcnode *nd = &w->t->nd[q[k].ni];
    if (!(nd->flags & HZ_DC_HASVERT)) {
      /* Ячейка у ребра со сменой знака ОБЯЗАНА иметь вершину; её нет только у
       * неманифолдной (Г47). Полигон не выдаём — fail closed.
       *
       * НО ОБХОД ПРОДОЛЖАЕТСЯ, И ЭТО ПОПРАВКА, А НЕ ПОСЛАБЛЕНИЕ (§377). Прежде
       * здесь ставился код возврата, и вся рекурсия сворачивалась: ОДНА
       * неманифолдная ячейка отменяла ВСЮ поверхность. Замерено — на зале при
       * L = 8 таких ячеек 2 538, и обход выдавал 170 треугольников вместо
       * ~200 тысяч, то есть картинки не было вовсе. «Fail closed» значит «не
       * выдать НЕВЕРНЫЙ полигон», а не «не выдать ни одного»: дыра в месте
       * отказа честнее пустого экрана. Отказ по-прежнему НЕ МОЛЧАЛИВ — он
       * считается в w->nskip и возвращается кодом HZ_DC_EMULTI в конце. */
      w->nskip++;
      return;
    }
    if (nv > 0 && rr[nv - 1].ni == q[k].ni) continue; /* крупная ячейка на двух квадрантах */
    rr[nv] = q[k];
    for (int a = 0; a < 3; a++)
      vv[nv][a] = nd->vx[a];
    nv++;
  }
  if (nv > 1 && rr[0].ni == rr[nv - 1].ni) nv--; /* и по кругу */
  if (nv < 3) return;

  if (!s0) { /* материал у ВЕРХНЕГО конца: наружу смотрит -e, разворачиваем */
    for (int a = 0, b = nv - 1; a < b; a++, b--) {
      hz_dcref tr = rr[a];
      rr[a] = rr[b];
      rr[b] = tr;
      for (int c = 0; c < 3; c++) {
        double td = vv[a][c];
        vv[a][c] = vv[b][c];
        vv[b][c] = td;
      }
    }
  }
  if (w->emit(w->pctx, rr, vv, nv) != 0 && w->rc == HZ_DC_OK) w->rc = HZ_DC_ETOPO;
}

static void edge_proc(walkctx *w, const hz_dcref q[4], int e, const int32_t qlo[3], int32_t seg) {
  if (w->rc != HZ_DC_OK) return;
  int u = (e + 1) % 3, v = (e + 2) % 3;
  int all = 1;
  for (int k = 0; k < 4; k++)
    if (!leafish(w, &q[k])) all = 0;
  if (all) {
    process_edge(w, q, e, qlo, seg);
    return;
  }
  int32_t half = seg > 1 ? seg / 2 : seg;
  int nh = seg > 1 ? 2 : 1;
  for (int h = 0; h < nh; h++) {
    int32_t qlo2[3] = {qlo[0], qlo[1], qlo[2]};
    qlo2[e] += (int32_t)h * half;
    hz_dcref nq[4];
    for (int k = 0; k < 4; k++) {
      if (leafish(w, &q[k])) {
        nq[k] = q[k];
        continue;
      }
      /* пробная точка внутри той ячейки, со стороны которой она смотрит на ребро */
      int32_t p2[3];
      p2[e] = 2 * qlo2[e] + 1;
      p2[u] = 2 * qlo2[u] + (2 * q[k].lo[u] + q[k].size > 2 * qlo2[u] ? 1 : -1);
      p2[v] = 2 * qlo2[v] + (2 * q[k].lo[v] + q[k].size > 2 * qlo2[v] ? 1 : -1);
      child_at2(w->t, &q[k], p2, &nq[k]);
    }
    edge_proc(w, nq, e, qlo2, half);
  }
}

static void face_proc(walkctx *w, const hz_dcref *r0, const hz_dcref *r1, int d,
                      const int32_t flo[3], int32_t rect) {
  if (w->rc != HZ_DC_OK) return;
  int l0 = leafish(w, r0), l1 = leafish(w, r1);
  if (l0 && l1) return;
  int u = (d + 1) % 3, v = (d + 2) % 3;
  int32_t half = rect > 1 ? rect / 2 : rect;
  int nq = rect > 1 ? 2 : 1;

  for (int su = 0; su < nq; su++)
    for (int sv = 0; sv < nq; sv++) {
      int32_t slo[3];
      slo[d] = flo[d];
      slo[u] = flo[u] + (int32_t)su * half;
      slo[v] = flo[v] + (int32_t)sv * half;
      int32_t p2[3];
      p2[u] = 2 * slo[u] + half;
      p2[v] = 2 * slo[v] + half;
      hz_dcref a = *r0, b = *r1;
      p2[d] = 2 * slo[d] - 1;
      if (!l0) child_at2(w->t, r0, p2, &a);
      p2[d] = 2 * slo[d] + 1;
      if (!l1) child_at2(w->t, r1, p2, &b);
      face_proc(w, &a, &b, d, slo, half);
    }

  if (rect <= 1) return;
  /* Рёбра ВНУТРИ общей грани: по одному кресту, два вдоль u и два вдоль v. */
  for (int t2 = 0; t2 < 2; t2++) {
    int ee = t2 ? v : u, other = t2 ? u : v;
    for (int hh = 0; hh < 2; hh++) {
      int32_t qlo[3];
      qlo[d] = flo[d];
      qlo[ee] = flo[ee] + (int32_t)hh * half;
      qlo[other] = flo[other] + half;
      hz_dcref in[4];
      int nin = 0;
      for (int side = 0; side < 2; side++)
        for (int so = 0; so < 2; so++) {
          const hz_dcref *base = side ? r1 : r0;
          if ((side ? l1 : l0)) {
            in[nin++] = *base;
            continue;
          }
          int32_t p2[3];
          p2[d] = 2 * qlo[d] + (side ? 1 : -1);
          p2[ee] = 2 * qlo[ee] + 1;
          p2[other] = 2 * qlo[other] + (so ? 1 : -1);
          child_at2(w->t, base, p2, &in[nin++]);
        }
      int32_t q2[3] = {2 * qlo[0], 2 * qlo[1], 2 * qlo[2]};
      hz_dcref quad[4];
      int rc = order_quad(in, nin, ee, q2, quad);
      if (rc != HZ_DC_OK) {
        w->rc = rc;
        return;
      }
      edge_proc(w, quad, ee, qlo, half);
    }
  }
}

static void cell_proc(walkctx *w, const hz_dcref *r) {
  if (w->rc != HZ_DC_OK || leafish(w, r)) return;
  int32_t half = r->size / 2;
  hz_dcref c[8];
  for (int i = 0; i < 8; i++) {
    c[i].ni = w->t->nd[r->ni].child0 + i;
    c[i].size = half;
    for (int a = 0; a < 3; a++)
      c[i].lo[a] = r->lo[a] + (((i >> a) & 1) ? half : 0);
  }
  for (int i = 0; i < 8; i++)
    cell_proc(w, &c[i]);

  for (int d = 0; d < 3; d++) {
    int u = (d + 1) % 3, v = (d + 2) % 3;
    for (int bu = 0; bu < 2; bu++)
      for (int bv = 0; bv < 2; bv++) {
        int ia = (bu << u) | (bv << v), ib = ia | (1 << d);
        int32_t flo[3];
        flo[d] = c[ia].lo[d] + half;
        flo[u] = c[ia].lo[u];
        flo[v] = c[ia].lo[v];
        face_proc(w, &c[ia], &c[ib], d, flo, half);
      }
  }

  for (int e = 0; e < 3; e++) {
    int u = (e + 1) % 3, v = (e + 2) % 3;
    for (int he = 0; he < 2; he++) {
      int32_t qlo[3];
      qlo[e] = r->lo[e] + (int32_t)he * half;
      qlo[u] = r->lo[u] + half;
      qlo[v] = r->lo[v] + half;
      hz_dcref in[4];
      int nin = 0;
      for (int bu = 0; bu < 2; bu++)
        for (int bv = 0; bv < 2; bv++)
          in[nin++] = c[(he << e) | (bu << u) | (bv << v)];
      int32_t q2[3] = {2 * qlo[0], 2 * qlo[1], 2 * qlo[2]};
      hz_dcref quad[4];
      int rc = order_quad(in, nin, e, q2, quad);
      if (rc != HZ_DC_OK) {
        w->rc = rc;
        return;
      }
      edge_proc(w, quad, e, qlo, half);
    }
  }
}

int hz_dc_walk(const hz_dctree *t, hz_dc_stop stop, void *sctx, hz_dc_poly emit, void *pctx) {
  walkctx w = {t, stop, sctx, emit, pctx, HZ_DC_OK, 0};
  hz_dcref root = {0, {0, 0, 0}, (int32_t)1 << t->log2size};
  cell_proc(&w, &root);
  if (w.rc == HZ_DC_OK && w.nskip > 0) return HZ_DC_EMULTI;
  return w.rc;
}

/* Спуск к ячейке cell среза (leafish). */
static void locate(const walkctx *w, const int32_t cell[3], hz_dcref *out) {
  hz_dcref r = {0, {0, 0, 0}, (int32_t)1 << w->t->log2size};
  while (!leafish(w, &r)) {
    int32_t p2[3] = {2 * cell[0] + 1, 2 * cell[1] + 1, 2 * cell[2] + 1};
    hz_dcref ch;
    child_at2(w->t, &r, p2, &ch);
    r = ch;
  }
  *out = r;
}

int hz_dc_walk_ref(const hz_dctree *t, const hz_htab *ht, hz_dc_stop stop, void *sctx,
                   hz_dc_poly emit, void *pctx) {
  walkctx w = {t, stop, sctx, emit, pctx, HZ_DC_OK, 0};
  int32_t n = (int32_t)1 << t->log2size;
  for (int32_t i = 0; i < ht->n && w.rc == HZ_DC_OK; i++) {
    const hz_hedge *e = &ht->e[i];
    int u = (e->axis + 1) % 3, v = (e->axis + 2) % 3;
    /* Ребро на границе куба делят меньше четырёх ячеек: поверхность там
     * открыта, и рекурсивный обход таких рёбер не видит по построению. */
    if (e->p[u] == 0 || e->p[u] == n || e->p[v] == 0 || e->p[v] == n) continue;
    hz_dcref in[4];
    for (int k = 0; k < 4; k++) {
      int32_t cell[3] = {e->p[0], e->p[1], e->p[2]};
      cell[u] -= (QUAD_UV[k][0] > 0 ? 0 : 1);
      cell[v] -= (QUAD_UV[k][1] > 0 ? 0 : 1);
      locate(&w, cell, &in[k]);
    }
    int32_t q2[3] = {2 * e->p[0], 2 * e->p[1], 2 * e->p[2]};
    hz_dcref quad[4];
    int rc = order_quad(in, 4, e->axis, q2, quad);
    if (rc != HZ_DC_OK) {
      w.rc = rc;
      break;
    }
    int32_t qlo[3] = {e->p[0], e->p[1], e->p[2]};
    /* Сегмент здесь всегда единичный: таблица хранит МИНИМАЛЬНЫЕ рёбра. Знак
     * берётся из самой записи, а не из маски — второй, независимый путь. */
    int mi = 0;
    for (int k = 1; k < 4; k++)
      if (quad[k].size < quad[mi].size) mi = k;
    if (quad[mi].size != 1) {
      /* срез грубее самого мелкого уровня: у эталона отображения «ребро ->
       * полигон» больше нет, и сверять его с рабочим обходом нельзя (Г42) */
      w.rc = HZ_DC_ETOPO;
      break;
    }
    process_edge(&w, quad, e->axis, qlo, 1);
  }
  return w.rc;
}

/* --- 5. мост к разрезу ----------------------------------------------------- */

typedef struct {
  int32_t cell, fref;
} cellfacet;

typedef struct {
  hz_facettab *ft;
  cellfacet *pair;
  int32_t np, cap;
  int32_t nbad; /* плоскость, отсекающая ячейку целиком: несогласованность */
  int rc;
} facetctx;

static int pair_cmp(const void *a, const void *b) {
  const cellfacet *x = a, *y = b;
  if (x->cell != y->cell) return x->cell < y->cell ? -1 : 1;
  if (x->fref != y->fref) return x->fref < y->fref ? -1 : 1;
  return 0;
}

/* Режет ли полуплоскость {n·x <= off} коробку. Тот же ТОЧНЫЙ отбор по опорной
 * вершине, что и hz_facets_for_box: нестрогий выкинул бы режущую плоскость
 * молча, а сумма объёмов всё равно сошлась бы (Г37). */
static int plane_cuts_box(const double nn[3], double off, const int32_t lo[3], int32_t size,
                          int *outside) {
  double vmax = 0.0, vmin = 0.0;
  for (int a = 0; a < 3; a++) {
    double l = (double)lo[a], h = (double)(lo[a] + size);
    vmax += nn[a] > 0.0 ? nn[a] * h : nn[a] * l;
    vmin += nn[a] > 0.0 ? nn[a] * l : nn[a] * h;
  }
  *outside = vmin > off;
  if (*outside) return 0;
  return vmax > off;
}

static int facet_emit(void *ctx, const hz_dcref *ref, const double (*v)[3], int nv) {
  facetctx *fc = ctx;
  int ntri = nv - 2;
  for (int tri = 0; tri < ntri; tri++) {
    const double *a = v[0], *b = v[tri + 1], *c = v[tri + 2];
    double e1[3], e2[3], nn[3];
    for (int k = 0; k < 3; k++) {
      e1[k] = b[k] - a[k];
      e2[k] = c[k] - a[k];
    }
    nn[0] = e1[1] * e2[2] - e1[2] * e2[1];
    nn[1] = e1[2] * e2[0] - e1[0] * e2[2];
    nn[2] = e1[0] * e2[1] - e1[1] * e2[0];
    double m = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
    if (!(m > 0.0)) continue; /* вырожденный треугольник — плоскости нет */
    for (int k = 0; k < 3; k++)
      nn[k] /= m;
    double off = nn[0] * a[0] + nn[1] * a[1] + nn[2] * a[2];
    /* Г44: dmax НЕ ВЫЧИСЛЕН и ноль сюда ставить нельзя — невязка QEF есть
     * среднеквадратичное на образцах, а не максимум смещения по фасету. */
    int32_t fi = hz_facettab_add_units(fc->ft, nn, off, -1, HZ_FACET_DMAX_UNKNOWN);
    if (fi < 0) {
      fc->rc = HZ_DC_ENOMEM;
      return 1;
    }
    /* Плоскость раздаётся ВСЕМ ячейкам полигона, а не только вершинам этого
     * треугольника: веер вокруг вершины ячейки обязан быть полным (Г10). */
    for (int k = 0; k < nv; k++) {
      int outside = 0;
      if (!plane_cuts_box(nn, off, ref[k].lo, ref[k].size, &outside)) {
        if (outside) fc->nbad++;
        continue;
      }
      if (fc->np >= fc->cap) {
        int32_t nc = fc->cap * 2;
        cellfacet *np2 = realloc(fc->pair, (size_t)nc * sizeof(cellfacet));
        if (np2 == NULL) {
          fc->rc = HZ_DC_ENOMEM;
          return 1;
        }
        fc->pair = np2;
        fc->cap = nc;
      }
      fc->pair[fc->np].cell = ref[k].ni;
      fc->pair[fc->np].fref = fi;
      fc->np++;
    }
  }
  return 0;
}

int hz_dc_facets(const hz_dctree *t, hz_dc_stop stop, void *sctx, hz_facettab *ft, hz_cutmap *cm) {
  facetctx fc = {ft, NULL, 0, 256, 0, HZ_DC_OK};
  fc.pair = calloc((size_t)fc.cap, sizeof(cellfacet));
  if (fc.pair == NULL) return HZ_DC_ENOMEM;
  int rc = hz_dc_walk(t, stop, sctx, facet_emit, &fc);
  if (rc == HZ_DC_OK) rc = fc.rc;
  if (rc != HZ_DC_OK) {
    free(fc.pair);
    return rc;
  }
  /* Г45: ОБХОД ВЫДАЁТ ЯЧЕЙКИ В ПОРЯДКЕ ДЕРЕВА, а hz_cutmap требует строго
   * возрастающего ключа и вернул бы 2. Пересортировка обязательна, и её код
   * возврата проверяется — молча потерянные фасеты дали бы ячейку без границы,
   * а сумма объёмов при этом всё равно сошлась бы. */
  if (fc.np > 1) qsort(fc.pair, (size_t)fc.np, sizeof(cellfacet), pair_cmp);
  int32_t i = 0;
  int32_t buf[HZ_P3_MAXH];
  while (i < fc.np) {
    int32_t cell = fc.pair[i].cell, nf = 0;
    int32_t j = i;
    while (j < fc.np && fc.pair[j].cell == cell) {
      if (j == i || fc.pair[j].fref != fc.pair[j - 1].fref) {
        if (nf >= HZ_P3_MAXH) {
          free(fc.pair);
          return HZ_DC_ETOPO; /* веер не влез в ядро — предел HZ_P3_MAXH, Г34 */
        }
        buf[nf++] = fc.pair[j].fref;
      }
      j++;
    }
    if (hz_cutmap_add(cm, cell, buf, nf) != 0) {
      free(fc.pair);
      return HZ_DC_ETOPO;
    }
    i = j;
  }
  free(fc.pair);
  return HZ_DC_OK;
}
