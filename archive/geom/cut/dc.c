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

/* СТЁРТОЕ РЕБРО НЕВИДИМО ДЛЯ ВСЕХ ПОТРЕБИТЕЛЕЙ (Р8, Ш4). Правка помечает запись
 * отдельным значением `in_lo`, а не порчей `t`. Первая редакция метила именно
 * `t < 0`, и это ЗАМЕРЕНО НЕГОДНЫМ: `test_dcwalk` упал `22/28`. Причина в том,
 * что `t` есть ДАННЫЕ, их читает знаковый путь, на котором стоят Г49 и Г50, —
 * то есть метка в поле данных ломает то, что этим полем пользуется. Поле
 * `in_lo` же не читается никем (А659) и потому годится под метку целиком.
 * Запись не вырезается из массива: таблица отсортирована по ключу, и вырезание
 * сдвигало бы всё за ней, то есть стоило бы O(таблицы) вместо O(затронутых).
 * Проверка стоит ЗДЕСЬ и только здесь — иначе каждый потребитель повторял бы её,
 * и первый же забывший вернул бы стёртую поверхность. */
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
  /* ЗНАКОВЫЙ РЕЖИМ ПО УМОЛЧАНИЮ, И ЭТО НЕ МЕЛОЧЬ (найдено 08-10, §412). Поле
   * `unsgn` заведено вместе с беззнаковым путём (Р7) и ставится только в
   * `hz_dc_shape_occ`; здесь оно НЕ ИНИЦИАЛИЗИРОВАЛОСЬ, а дерево у всех
   * потребителей лежит на стеке. Тогда `process_edge` читал мусор и в половине
   * прогонов уходил в беззнаковую ветвь на ЗНАКОВОМ дереве, где маски рёбер
   * пусты, — обход выдавал ноль многоугольников. `test_dcwalk` от этого мигал
   * 50/50, и мигание было принято мною за регресс шага Ш4.
   * `nbigmask` — по той же причине: его печатают как «обязан быть 0». */
  t->unsgn = 0;
  t->nbigmask = 0;
  t->vbits = HZ_DC_VBITS;
  t->nsum = NULL;
  t->nsumcap = 0;
  t->nzeronrm = 0;
  t->nd = calloc((size_t)t->cap, sizeof(hz_dcnode));
  t->qf = calloc((size_t)t->cap, sizeof(hz_qef));
  if (t->nd == NULL || t->qf == NULL) {
    free(t->nd);
    free(t->qf);
    t->nd = NULL;
    t->qf = NULL;
    return HZ_DC_ENOMEM;
  }
  t->nd[0].child0 = -1;
  hz_qef_zero(&t->qf[0]);
  t->n = 1;
  return HZ_DC_OK;
}

void hz_dc_free(hz_dctree *t) {
  free(t->nd);
  t->nd = NULL;
  free(t->qf);
  t->qf = NULL;
  free(t->nsum);
  t->nsum = NULL;
  t->nsumcap = 0;
  t->n = t->cap = 0;
}

/* Индекс побочного массива по узлу. ПРОВЕРЯЕТСЯ, а не выводится из веры в
 * распределитель: `child0` обязан быть >= 1 и кратен 8 плюс 1 (блоки идут с
 * единицы, по восемь подряд). Не так — величины нет, и это не молчание, а
 * возврат «нет». */
static int32_t nsum_slot(const hz_dctree *t, int32_t ni) {
  if (ni < 0 || ni >= t->n) return -1;
  int32_t c0 = t->nd[ni].child0;
  if (c0 < 1 || ((c0 - 1) % 8) != 0) return -1;
  int32_t k = (c0 - 1) / 8;
  return k < t->nsumcap ? k : -1;
}

int hz_dc_nsum(const hz_dctree *t, int32_t ni, double s[3]) {
  s[0] = s[1] = s[2] = 0.0;
  int32_t k = nsum_slot(t, ni);
  if (k < 0 || t->nsum == NULL) return 0;
  for (int a = 0; a < 3; a++)
    s[a] = t->nsum[k][a];
  return 1;
}

/* Побочный массив держится под ЧИСЛО УЗЛОВ ДЕРЕВА, а не под точное число
 * внутренних: считать вторые пришлось бы отдельным проходом, а память та же с
 * точностью до восьмушки. Заводится ПОСЛЕ спуска, когда дерево уже не растёт. */
static int nsum_alloc(hz_dctree *t) {
  int32_t need = t->n / 8 + 1;
  if (t->nsum != NULL && t->nsumcap >= need) {
    memset(t->nsum, 0, (size_t)t->nsumcap * sizeof *t->nsum);
    return HZ_DC_OK;
  }
  free(t->nsum);
  t->nsum = calloc((size_t)need, sizeof *t->nsum);
  if (t->nsum == NULL) {
    t->nsumcap = 0;
    return HZ_DC_ENOMEM;
  }
  t->nsumcap = need;
  return HZ_DC_OK;
}

/* Блок из 8 подряд, строго вперёд — дословно node_alloc8 октодерева: на «индекс
 * ребёнка больше индекса родителя» опирается и завершаемость спуска. */
void hz_dc_drop_forms(hz_dctree *t) {
  free(t->qf);
  t->qf = NULL;
}

static int32_t dc_alloc8(hz_dctree *t) {
  if (t->n + 8 > t->cap) {
    int32_t nc = t->cap * 2;
    if (nc < t->n + 8) nc = t->n + 8;
    hz_dcnode *nn = realloc(t->nd, (size_t)nc * sizeof(hz_dcnode));
    if (nn == NULL) return -1;
    t->nd = nn;
    /* Формы растут ВМЕСТЕ с узлами, пока не сняты (§632). Сняты — не растут, и
     * это не молчание: дробить дерево без форм нечем, и такой путь отсекается
     * в `hz_dc_repair`. */
    if (t->qf != NULL) {
      hz_qef *nq2 = realloc(t->qf, (size_t)nc * sizeof(hz_qef));
      if (nq2 == NULL) return -1;
      t->qf = nq2;
    }
    t->cap = nc;
  }
  int32_t base = t->n;
  for (int i = 0; i < 8; i++) {
    memset(&t->nd[base + i], 0, sizeof(hz_dcnode));
    t->nd[base + i].child0 = -1;
    hz_qef_zero(&t->qf[base + i]);
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
  hz_qef_zero(&t->qf[ni]);
  for (int i = 0; i < 12; i++) {
    int axis, off[3];
    unit_edge(i, &axis, off);
    int32_t p[3] = {lo[0] + off[0], lo[1] + off[1], lo[2] + off[2]};
    const hz_hedge *e = hz_htab_find(ht, axis, p);
    if (e == NULL || e->in_lo == HZ_HEDGE_ERASED) continue;
    /* Локальная координата считается ЦЕЛОЧИСЛЕННОЙ разностью плюс t, а не
     * вычитанием двух больших чисел: сокращения нет вовсе (Г41). */
    double loc[3];
    for (int k = 0; k < 3; k++)
      loc[k] = (double)(e->p[k] - lo[k]) + (k == axis ? e->t : 0.0);
    hz_qef_add_sample(&t->qf[ni], loc, e->nrm);
  }
}

/* §635: `lo` БОЛЬШЕ НЕ НУЖЕН. Вершина хранится ОТНОСИТЕЛЬНО коробки, поэтому
 * абсолютная координата нижнего угла в решение не входит вовсе — её прибавляет
 * читатель, у которого коробка есть по построению обхода. */
static void solve_node(hz_dctree *t, int32_t ni, int32_t size) {
  hz_dcnode *nd = &t->nd[ni];
  if (nd->flags & HZ_DC_CLAMPED) t->nclamped--;
  if (nd->flags & HZ_DC_MULTI) t->nmulti--;
  nd->flags = 0;
  nd->err = 0.0f;
  for (int k = 0; k < 3; k++)
    nd->vq[k] = 0;
  /* §632: число образцов ПЕРЕПИСЫВАЕТСЯ В УЗЕЛ здесь и только здесь — это
   * единственное, что нужно от формы после постройки (`hz_dc_rms`). */
  nd->nq = t->qf[ni].n;
  if (t->qf[ni].n <= 0) return;
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
  int st = hz_qef_solve(&t->qf[ni], blo, bhi, x, &r);
  if (st == HZ_QEF_EMPTY) return;
  if (st == HZ_QEF_CLAMPED) {
    nd->flags |= HZ_DC_CLAMPED;
    t->nclamped++;
  }
  /* КВАНТОВАНИЕ ВЕРШИНЫ (§635; разрядность исправлена §636 с `8` на `16`).
   * Хранится ДОЛЯ коробки. Округление к БЛИЖАЙШЕМУ, а не отсечение: отсечение
   * дало бы систематический сдвиг к нижнему углу узла, то есть СМЕЩЕНИЕ
   * поверхности, а не шум. Зажим по краям — решатель уже держит `x` в коробке,
   * но краевое значение `size` обязано лечь ровно в `65535`.
   * ШКАЛА ХРАНЕНИЯ ВСЕГДА `65535`, СКОЛЬКО БЫ БИТ НИ БЫЛО ЗНАЧИМО (§638): при
   * `t->vbits < 16` доля сначала округляется к сетке из `qmax + 1` уровней, а
   * затем растягивается обратно в шкалу `65535`. Так `hz_dc_vertex` остаётся
   * без единой лишней операции, а негативный контроль `nodevq1` получается
   * ТОЧНЫМ: `qmax = 1` даёт `vq` ровно `0` или `65535`, то есть долю `0` или
   * `1`, а не `0` или `0.5`. */
  /* Разрядность вне `[1, HZ_DC_VBITS]` есть ошибка вызывающего, а не режим.
   * Зажим стоит здесь потому, что `qmax = 0` дал бы деление на ноль ниже, а
   * молчаливое падение хуже зажатого кванта; сам зажим ничего не «чинит»: при
   * `vbits = HZ_DC_VBITS` он тождествен. */
  int vb = t->vbits < 1 ? 1 : (t->vbits > HZ_DC_VBITS ? HZ_DC_VBITS : t->vbits);
  uint32_t qmax = (1u << (unsigned)vb) - 1u;
  for (int k = 0; k < 3; k++) {
    double f = x[k] / (double)size;
    if (!(f > 0.0)) f = 0.0;
    if (f > 1.0) f = 1.0;
    uint32_t u = (uint32_t)(f * (double)qmax + 0.5);
    if (u > qmax) u = qmax;
    nd->vq[k] = (uint16_t)((uint64_t)u * 65535u / qmax);
  }
  nd->err = (float)r;
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
    hz_qef_add_shifted(&acc, &t->qf[c0 + i], sh);
  }
  t->qf[ni] = acc;
}

static int build_rec(hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size,
                     const hz_signgrid *g, const hz_htab *ht) {
  t->nd[ni].corner = corner_mask(g, lo, size);
  t->nd[ni].child0 = -1;
  hz_qef_zero(&t->qf[ni]);
  if (box_uniform(g, lo, size)) return HZ_DC_OK; /* поверхности внутри нет */

  if (size == 1) {
    leaf_qef(t, ni, lo, ht);
    solve_node(t, ni, size);
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
  solve_node(t, ni, size);
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
    solve_node(t, ni, size);
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
  solve_node(t, ni, size);
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
  hz_qef_zero(&t->qf[ni]);
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
static void leaf_masks(hz_dctree *t, int32_t ni, const int32_t lo[3], const hz_htab *ht,
                       double nsum[3]);

static void masks_occ_rec(hz_dctree *t, int32_t ni, const int32_t lo[3], int32_t size,
                          const hz_htab *ht, double nsum[3]) {
  hz_dcnode *nd = &t->nd[ni];
  nsum[0] = nsum[1] = nsum[2] = 0.0;
  if (nd->child0 < 0) {
    /* ОДИН КОД С ПОЧИНКОЙ, А НЕ ВТОРАЯ КОПИЯ ФОРМУЛЫ (найдено 08-11, А841:
     * комментарий над `leaf_masks` обещал единый код, а копий было две — здесь
     * своя. Побитовость починки против сборки на этом и стоит). */
    nd->ecross = 0;
    nd->edir = 0;
    if (size != 1) return;
    leaf_masks(t, ni, lo, ht, nsum);
    return;
  }
  int32_t c0 = nd->child0, half = size / 2;
  for (int k = 0; k < 8; k++) {
    int32_t clo[3];
    for (int a = 0; a < 3; a++)
      clo[a] = lo[a] + (((k >> a) & 1) ? half : 0);
    double cs[3];
    masks_occ_rec(t, c0 + k, clo, half, ht, cs);
    for (int a = 0; a < 3; a++)
      nsum[a] += cs[a];
  }
  {
    int32_t sl = nsum_slot(t, ni);
    if (sl >= 0)
      for (int a = 0; a < 3; a++)
        t->nsum[sl][a] = nsum[a];
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
  int rc = nsum_alloc(t);
  if (rc != HZ_DC_OK) return rc;
  double root[3];
  masks_occ_rec(t, 0, zero, (int32_t)1 << t->log2size, ht, root);
  t->nbigmask = 0;
  masks_check_rec(t, 0, (int32_t)1 << t->log2size);
  /* А834: сумма может сократиться ТОЧНО — двусторонний лист внутри узла. Тогда
   * направления нет, и это не «нормаль (0,0,0)», а отсутствие величины. Здесь
   * оно СЧИТАЕТСЯ; молчать про остаток запрещено. */
  t->nzeronrm = 0;
  for (int32_t i = 0; i < t->n; i++) {
    int32_t sl = nsum_slot(t, i);
    if (sl < 0) continue;
    const double *s = t->nsum[sl];
    if (!(s[0] < 0.0) && !(s[0] > 0.0) && !(s[1] < 0.0) && !(s[1] > 0.0) && !(s[2] < 0.0) &&
        !(s[2] > 0.0))
      t->nzeronrm++;
  }
  return HZ_DC_OK;
}

/* --- 3в. ЛОКАЛЬНАЯ ПОЧИНКА ПОСЛЕ ПРАВКИ (Ш4, Р8) --------------------------- */

/* Маска ОДНОГО листа из таблицы — та же формула, что в , вынесена
 * затем, чтобы починка и полная сборка считали её ОДНИМ кодом: разойдись они,
 * побитовая сверка Г49 мерила бы разницу двух копий формулы. */
static void leaf_masks(hz_dctree *t, int32_t ni, const int32_t lo[3], const hz_htab *ht,
                       double nsum[3]) {
  hz_dcnode *nd = &t->nd[ni];
  nd->ecross = 0;
  nd->edir = 0;
  if (nsum != NULL) nsum[0] = nsum[1] = nsum[2] = 0.0;
  for (int i = 0; i < 12; i++) {
    int axis, off[3];
    unit_edge(i, &axis, off);
    int32_t p[3] = {lo[0] + off[0], lo[1] + off[1], lo[2] + off[2]};
    const hz_hedge *e = hz_htab_find(ht, axis, p);
    if (e == NULL || e->in_lo == HZ_HEDGE_ERASED) continue;
    nd->ecross = (uint16_t)(nd->ecross | (1u << i));
    if (e->nrm[axis] > 0.0) nd->edir = (uint16_t)(nd->edir | (1u << i));
    /* Ш9: сумма нормалей набирается из ТЕХ ЖЕ найденных записей. */
    if (nsum != NULL)
      for (int a = 0; a < 3; a++)
        nsum[a] += e->nrm[a];
  }
}

static void node_masks_from_children(hz_dctree *t, int32_t ni) {
  int32_t c0 = t->nd[ni].child0;
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

/* Сумма нормалей узла ИЗ ДЕТЕЙ — для починки за O(глубины) (Ш9). Ребёнок отдаёт
 * величину тремя разными путями, и все три законны: внутренний — из побочного
 * массива; лист размера 1 — пересчётом из таблицы (12 поисков, как при сборке);
 * крупный лист — ноль, потому что рёбер у него нет вовсе (А709). Иначе крупные
 * нормали остались бы после удара несвежими МОЛЧА: сверка Г49 сличает ячейки
 * среза, и до Ш9 бита ориентации в них не было вовсе. */
static void node_nsum_from_children(hz_dctree *t, const hz_htab *ht, int32_t ni,
                                    const int32_t lo[3], int32_t size) {
  int32_t sl = nsum_slot(t, ni);
  if (sl < 0) return;
  int32_t c0 = t->nd[ni].child0, half = size / 2;
  double acc[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 8; k++) {
    int32_t ci = c0 + k, clo[3];
    for (int a = 0; a < 3; a++)
      clo[a] = lo[a] + (((k >> a) & 1) ? half : 0);
    double cs[3] = {0.0, 0.0, 0.0};
    if (t->nd[ci].child0 >= 0) {
      hz_dc_nsum(t, ci, cs);
    } else if (half == 1) {
      hz_dcnode save = t->nd[ci];
      leaf_masks(t, ci, clo, ht, cs);
      t->nd[ci].ecross = save.ecross; /* маску трогать не наше дело: её уже */
      t->nd[ci].edir = save.edir;     /* пересчитал node_masks_from_children */
    }
    for (int a = 0; a < 3; a++)
      acc[a] += cs[a];
  }
  for (int a = 0; a < 3; a++)
    t->nsum[sl][a] = acc[a];
}

int hz_dc_fix_cell(hz_dctree *t, const hz_htab *ht, const int32_t cell[3]) {
  if (t->qf == NULL) return HZ_DC_ENOMEM;
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
  if (size != 1) return HZ_DC_OK; /* ячейка под неразделённым узлом — трогать нечего */
  leaf_masks(t, ni, lo, ht, NULL);
  leaf_qef(t, ni, lo, ht);
  solve_node(t, ni, size);
  /* Подъём: форма и маска родителя пересчитываются ИЗ ДЕТЕЙ, то есть O(глубины),
   * а не по поддереву. Это и есть обещание §383 п. 4. */
  for (int d = depth - 1; d >= 0; d--) {
    sum_children(t, path[d], psize[d]);
    node_masks_from_children(t, path[d]);
    node_nsum_from_children(t, ht, path[d], plo[d], psize[d]);
    solve_node(t, path[d], psize[d]);
  }
  return HZ_DC_OK;
}

int hz_dc_forms_lazy(hz_dctree *t, const hz_htab *ht) {
  int32_t zero[3] = {0, 0, 0};
  lazy_forms(t, 0, zero, (int32_t)1 << t->log2size, ht);
  return t->nmulti > 0 ? HZ_DC_EMULTI : HZ_DC_OK;
}

int hz_dc_repair(hz_dctree *t, const hz_htab *ht, const int32_t cell[3]) {
  /* §632: без форм починка НЕВОЗМОЖНА, и это отказ, а не тишина (Г25). */
  if (t->qf == NULL) return HZ_DC_ENOMEM;
  /* §635: путь коробок предков больше не нужен — `solve_node` работает в
   * ОТНОСИТЕЛЬНЫХ координатах, а `lo` предка ему не требуется. */
  int32_t path[HZ_DC_MAX_LOG2SIZE + 1], psize[HZ_DC_MAX_LOG2SIZE + 1];
  int depth = 0;
  int32_t ni = 0, size = (int32_t)1 << t->log2size, lo[3] = {0, 0, 0};
  for (int a = 0; a < 3; a++)
    if (cell[a] < 0 || cell[a] >= size) return HZ_DC_ERANGE;

  while (size > 1 && t->nd[ni].child0 >= 0) {
    path[depth] = ni;
    psize[depth] = size;
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
  solve_node(t, ni, size);
  for (int d = depth - 1; d >= 0; d--) {
    sum_children(t, path[d], psize[d]);
    solve_node(t, path[d], psize[d]);
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
  /* ПАМЯТЬ ОТВЕТА КРИТЕРИЯ СРЕЗА, по байту на узел (§585). NULL = считать
   * каждый раз, как было. Живёт ровно один вызов обхода. */
  unsigned char *memo;
} walkctx;

/* --- память ответа критерия среза (§585) ----------------------------------
 *
 * ЗАЧЕМ. `leafish` есть функция ОДНОГО аргумента — индекса узла: коробка
 * `(lo, size)` у узла единственная (арена, один родитель на блок из восьми,
 * индексы строго вперёд — `dc_alloc8`), а контекст среза на время обхода
 * неподвижен. Схема же `cellProc/faceProc/edgeProc` спрашивает про один и тот
 * же узел много раз: как про ячейку, как про обе стороны каждой его грани и как
 * про каждый из четырёх углов каждого его ребра, и заново на каждом уровне
 * спуска. ЗАМЕРЕНО (§585, Bistro, poly=8): `16 037 694` вызова `lod_stop` при
 * не более чем `276 904` внутренних узлах — кратность `57.9`.
 *
 * ПОЧЕМУ ЭТО НЕ МЕНЯЕТ ВЫХОД. Значение то же самое; меняется только число его
 * вычислений. Побитовость картинки — не надежда, а следствие, и проверяется она
 * прогоном с `HZ_DC_MEMO_OFF`.
 *
 * ПЕРЕКЛЮЧАТЕЛЬ — ПРИБОР, А НЕ ЧАСТЬ ДОГОВОРА. Он нужен ровно затем, чтобы
 * негативные контроли §585 существовали: без `OFF` нечем показать, что выигрыш
 * пришёл отсюда, а без `SCRAMBLE` нечем показать, что побитовое сличение вообще
 * что-нибудь заметило бы. Обход однопоточный, поэтому файловая статика здесь не
 * создаёт вопроса о потоках. */
static int g_memo_mode = HZ_DC_MEMO_ON;
/* ДВА СЧЁТЧИКА, А НЕ ОДИН, И РАЗНИЦА НЕ КОСМЕТИЧЕСКАЯ. Первая редакция считала
 * одним, и он вышел `351 113` при `276 904` внутренних узлах — то есть считал
 * ЛЮБЫЕ узлы, включая листья, а сверялся с числом ВНУТРЕННИХ. Ровно подмена
 * величины, против которой стоит §4. `evals` — сколько РАЗЛИЧНЫХ узлов вообще
 * тронул обход; `stops` — у скольких из них критерий среза вычислен полностью,
 * и только это число сравнимо с прежними `16 037 694` вызовами. */
static long long g_memo_evals, g_memo_stops, g_memo_flips;

/* Шаг порчи для `HZ_DC_MEMO_SCRAMBLE`. Простое число, и это не украшение:
 * дети выделяются блоками ПО ВОСЕМЬ, поэтому любой шаг, кратный восьми, попадал
 * бы всегда в один и тот же угол блока. Первая редакция плана портила ответ у
 * соседнего ребёнка (`ni ^ 1`), и аудит А1014 её отверг: у детей одного родителя
 * дальность и невязка близки, ответ у них чаще всего ОДИН И ТОТ ЖЕ, и порча
 * оказалась бы незаметной — то есть контроль был бы слеп ровно к тому, ради чего
 * поставлен (узор К13/К40/К94).
 *
 * КОРЕНЬ ИСКЛЮЧЁН, И ЭТО ВТОРАЯ ПРАВКА КОНТРОЛЯ, СДЕЛАННАЯ ПО ЕГО ЖЕ ПРОГОНУ
 * (А1022). `0 % 997 == 0`, поэтому первая редакция переворачивала ответ У КОРНЯ:
 * обход объявлял корень листом и кончался немедленно — `1` тронутый узел,
 * `0` многоугольников, растеризация `4.6` мс. Картинка, конечно, разошлась, но
 * доказано этим было лишь «если убить корень, кадра не будет», а проверить надо
 * другое: заметит ли побитовое сличение ЛОКАЛЬНУЮ порчу в глубине дерева. */
#define HZ_DC_MEMO_SCRAMBLE_STRIDE 997

void hz_dc_walk_memo(int mode) {
  g_memo_mode = mode;
}
long long hz_dc_walk_memo_evals(void) {
  return g_memo_evals;
}
long long hz_dc_walk_memo_stops(void) {
  return g_memo_stops;
}
long long hz_dc_walk_memo_flips(void) {
  return g_memo_flips;
}

/* СЧЁТЧИКИ ОБХОДА — ТОЛЬКО ПОД `-DHZ_DC_COUNT`, В РАБОЧЕЙ СБОРКЕ ИХ НЕТ.
 * Заведены, чтобы не гадать, где именно стоит обход: §573 уже поймал проект на
 * том, что «растеризация» оказалась спуском по дереву, а §580 — на том, что
 * лишними были многоугольники, а не обход. Считать инкремент в горячем цикле
 * рабочего пути нельзя, поэтому счёт живёт за флагом сборки. */
#ifdef HZ_DC_COUNT
long long hz_dc_n_cell, hz_dc_n_face, hz_dc_n_edge, hz_dc_n_leafish, hz_dc_n_stop, hz_dc_n_proc;
#define HZ_CNT(x) ((x)++)
#else
#define HZ_CNT(x) ((void)0)
#endif

static int leafish(const walkctx *w, const hz_dcref *r) {
  HZ_CNT(hz_dc_n_leafish);
  if (w->memo != NULL) {
    unsigned char m = w->memo[r->ni];
    if (m != 0) return m == 1;
  }
  HZ_CNT(hz_dc_n_stop);
  int haskids = w->t->nd[r->ni].child0 >= 0;
  /* СЧЁТ ПОЛНЫХ ВЫЧИСЛЕНИЙ КРИТЕРИЯ ИДЁТ В ОБОИХ РЕЖИМАХ, И ЭТО НЕ МЕЛОЧЬ.
   * Первая редакция считала его только при включённой памяти ответа — и тогда
   * «было» и «стало» оказывались разными величинами, снятыми разными приборами
   * (А1023). Сравнивать можно только это число с ним же. */
  if (haskids) g_memo_stops++;
  int lf = !haskids || (w->stop != NULL && w->stop(w->sctx, w->t, r));
  if (w->memo != NULL) {
    g_memo_evals++;
    /* ПОРЧА ТОЛЬКО У УЗЛА С ДЕТЬМИ. У листа «не лист» означало бы спуск по
     * `child0 = −1`, то есть чтение мимо массива: негативный контроль обязан
     * ломать ОТВЕТ, а не память. У узла с детьми оба значения законны. */
    if (g_memo_mode == HZ_DC_MEMO_SCRAMBLE && haskids && r->ni > 0 &&
        (r->ni % HZ_DC_MEMO_SCRAMBLE_STRIDE) == 0) {
      lf = !lf;
      g_memo_flips++;
    }
    w->memo[r->ni] = (unsigned char)(lf ? 1 : 2);
  }
  return lf;
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
  HZ_CNT(hz_dc_n_proc);
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
      vv[nv][a] = (double)q[k].lo[a] + (double)nd->vq[a] * (double)q[k].size / 65535.0;
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
  HZ_CNT(hz_dc_n_edge);
  if (w->rc != HZ_DC_OK) return;
  int u = (e + 1) % 3, v = (e + 2) % 3;
  /* Р2 (§585): ОТВЕТ ПО ЧЕТЫРЁМ ЯЧЕЙКАМ СЧИТАЕТСЯ ОДИН РАЗ. Прежде он считался
   * здесь, а потом ЗАНОВО в цикле по половинам — до восьми лишних вызовов на
   * каждый рекурсирующий `edge_proc`. Порядок обхода правка не трогает: те же
   * циклы, та же рекурсия, то же условие. */
  int lf[4], all = 1;
  for (int k = 0; k < 4; k++) {
    lf[k] = leafish(w, &q[k]);
    if (!lf[k]) all = 0;
  }
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
      if (lf[k]) {
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
  HZ_CNT(hz_dc_n_face);
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
  HZ_CNT(hz_dc_n_cell);
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

/* Выделение памяти ответа. `calloc` СТОИТ ЗДЕСЬ, А НЕ В ОБЁРТКЕ, сознательно:
 * `CLAUDE.md` называет известный класс ложных срабатываний gcc-analyzer, при
 * котором связь «ёмкость ↔ счёт» теряется через функцию-распределитель. Ноль
 * значит «не считано», и это ЕДИНСТВЕННОЕ значение по умолчанию, которое здесь
 * законно: «нет пометки» = «неизвестно», а не «лист» (А917).
 * Нехватка памяти — не ошибка: NULL возвращает прежний путь слово в слово. */
static unsigned char *memo_alloc(const hz_dctree *t) {
  g_memo_evals = 0;
  g_memo_stops = 0;
  g_memo_flips = 0;
  if (g_memo_mode == HZ_DC_MEMO_OFF || t->n <= 0) return NULL;
  return calloc((size_t)t->n, 1);
}

int hz_dc_walk_stats(const hz_dctree *t, hz_dc_stop stop, void *sctx, hz_dc_poly emit, void *pctx,
                     int32_t *nskip) {
  walkctx w = {t, stop, sctx, emit, pctx, HZ_DC_OK, 0, memo_alloc(t)};
  hz_dcref root = {0, {0, 0, 0}, (int32_t)1 << t->log2size};
  cell_proc(&w, &root);
  free(w.memo);
  if (nskip != NULL) *nskip = w.nskip;
  if (w.rc == HZ_DC_OK && w.nskip > 0) return HZ_DC_EMULTI;
  return w.rc;
}

int hz_dc_walk(const hz_dctree *t, hz_dc_stop stop, void *sctx, hz_dc_poly emit, void *pctx) {
  return hz_dc_walk_stats(t, stop, sctx, emit, pctx, NULL);
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
  /* Эталонный обход (Г42) получает ту же память ответа: он зовёт `leafish`
   * через `locate`, и считать её двумя разными приборами значило бы сличать не
   * то. Независимость эталона от этого не страдает — она в том, что он идёт ПО
   * ТАБЛИЦЕ РЁБЕР, а не рекурсией (А1015). */
  walkctx w = {t, stop, sctx, emit, pctx, HZ_DC_OK, 0, memo_alloc(t)};
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
  free(w.memo);
  return w.rc;
}

/* --- 5. мост к разрезу ----------------------------------------------------- */

typedef struct {
  int32_t cell, fref;
} cellfacet;

typedef struct {
  hz_facettab *ft;
  const hz_dctree *t; /* Р-8: раздача по куску спускается по дереву */
  cellfacet *pair;
  int32_t np, cap;
  int32_t nbad;   /* плоскость, отсекающая ячейку целиком: несогласованность */
  int bounded;    /* Р-8: раздача по КУСКУ (иначе — по плоскости, как прежде) */
  int32_t *dropc; /* Р-8/Г60: ячейки пар, отброшенных кусочным отбором */
  int32_t nd, dcap;
  int rc;
} facetctx;

static int push_pair(facetctx *fc, int32_t cell, int32_t fref) {
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
  fc->pair[fc->np].cell = cell;
  fc->pair[fc->np].fref = fref;
  fc->np++;
  return 0;
}

/* Р-8: счётчики раздачи — читаются потребителем после hz_dc_facets (стиль
 * hz_dc_walk_memo_*): молча сузившаяся раздача неотличима от бага (Г58/Г60). */
static long long g_fc_pairs = 0, g_fc_dropped = 0, g_fc_empty = 0;
long long hz_dc_facets_pairs(void) {
  return g_fc_pairs;
}
long long hz_dc_facets_dropped(void) {
  return g_fc_dropped;
}
long long hz_dc_facets_empty(void) {
  return g_fc_empty;
}

/* Р-8: пересекает ли ТРЕУГОЛЬНИК коробку [lo, lo+size) — SAT Акенина-Мёллера.
 * НЕСТРОГИЕ сравнения сознательно (Г59): касание = пересечение; ложное «да»
 * безвредно (рез даст пустую грань), ложное «нет» — потерянный кусок = дыра. */
static int tri_box_overlap(const double v0[3], const double v1[3], const double v2[3],
                           const int32_t lo[3], int32_t size) {
  double h = 0.5 * (double)size;
  double a[3][3];
  for (int k = 0; k < 3; k++) {
    double c = (double)lo[k] + h;
    a[0][k] = v0[k] - c;
    a[1][k] = v1[k] - c;
    a[2][k] = v2[k] - c;
  }
  for (int k = 0; k < 3; k++) { /* оси коробки */
    double mn = a[0][k], mx = a[0][k];
    for (int i = 1; i < 3; i++) {
      if (a[i][k] < mn) mn = a[i][k];
      if (a[i][k] > mx) mx = a[i][k];
    }
    if (mn > h || mx < -h) return 0;
  }
  double e[3][3], nn[3];
  for (int k = 0; k < 3; k++) {
    e[0][k] = a[1][k] - a[0][k];
    e[1][k] = a[2][k] - a[1][k];
    e[2][k] = a[0][k] - a[2][k];
  }
  nn[0] = e[0][1] * e[1][2] - e[0][2] * e[1][1]; /* нормаль треугольника */
  nn[1] = e[0][2] * e[1][0] - e[0][0] * e[1][2];
  nn[2] = e[0][0] * e[1][1] - e[0][1] * e[1][0];
  {
    double d = nn[0] * a[0][0] + nn[1] * a[0][1] + nn[2] * a[0][2];
    double r = h * (fabs(nn[0]) + fabs(nn[1]) + fabs(nn[2]));
    if (d > r || d < -r) return 0;
  }
  for (int i = 0; i < 3; i++) /* девять осей ребро × орт */
    for (int k = 0; k < 3; k++) {
      double ax[3] = {0.0, 0.0, 0.0};
      ax[(k + 1) % 3] = e[i][(k + 2) % 3];
      ax[(k + 2) % 3] = -e[i][(k + 1) % 3];
      double p0 = ax[0] * a[0][0] + ax[1] * a[0][1] + ax[2] * a[0][2];
      double p1 = ax[0] * a[1][0] + ax[1] * a[1][1] + ax[2] * a[1][2];
      double p2 = ax[0] * a[2][0] + ax[1] * a[2][1] + ax[2] * a[2][2];
      double mn = p0 < p1 ? (p0 < p2 ? p0 : p2) : (p1 < p2 ? p1 : p2);
      double mx = p0 > p1 ? (p0 > p2 ? p0 : p2) : (p1 > p2 ? p1 : p2);
      double r = h * (fabs(ax[0]) + fabs(ax[1]) + fabs(ax[2]));
      if (mn > r || mx < -r) return 0;
    }
  return 1;
}

static int pair_cmp(const void *a, const void *b) {
  const cellfacet *x = a, *y = b;
  if (x->cell != y->cell) return x->cell < y->cell ? -1 : 1;
  if (x->fref != y->fref) return x->fref < y->fref ? -1 : 1;
  return 0;
}

static int pair_cmp_i32(const void *a, const void *b) {
  int32_t x = *(const int32_t *)a, y = *(const int32_t *)b;
  return x < y ? -1 : (x > y ? 1 : 0);
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

/* Р-8: спуск по дереву — пара (лист, фасет) каждому листу, чей бокс пересекает
 * кусок. SAT консервативен (Г59), спуск отсекает поддеревья без пересечения.
 * Возврат 1 — ошибка памяти (rc уже выставлен push_pair). */
static int piece_pairs_rec(facetctx *fc, int32_t ni, const int32_t lo[3], int32_t size,
                           const double *a, const double *b, const double *c, int32_t fi) {
  if (!tri_box_overlap(a, b, c, lo, size)) return 0;
  if (fc->t->nd[ni].child0 < 0) {
    /* Пары — только ЛИСТЬЯМ ДНА: поверхность живёт на самом мелком уровне
     * (условие 1:1 cut3), а крупный лист занятости в сетке переноса нумеруется
     * иначе — запись на нём потерялась бы (замерено: room, se в «крупных»
     * 18 751 при раздаче любым листьям). Крупный лист сюда попадает только
     * КАСАНИЕМ границы (консервативный SAT), и пара ему не нужна. */
    if (size == 1) return push_pair(fc, ni, fi);
    return 0;
  }
  int32_t half = size / 2;
  for (int i = 0; i < 8; i++) {
    int32_t clo[3] = {lo[0] + ((i & 1) ? half : 0), lo[1] + ((i & 2) ? half : 0),
                      lo[2] + ((i & 4) ? half : 0)};
    if (piece_pairs_rec(fc, fc->t->nd[ni].child0 + i, clo, half, a, b, c, fi)) return 1;
  }
  return 0;
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
    int32_t fi;
    if (fc->bounded) { /* Р-8: фасет несёт свой кусок */
      double tv8[3][3];
      for (int k = 0; k < 3; k++) {
        tv8[0][k] = a[k];
        tv8[1][k] = b[k];
        tv8[2][k] = c[k];
      }
      fi = hz_facettab_add_units_piece(fc->ft, nn, off, -1, HZ_FACET_DMAX_UNKNOWN,
                                       (const double (*)[3])tv8);
    } else
      fi = hz_facettab_add_units(fc->ft, nn, off, -1, HZ_FACET_DMAX_UNKNOWN);
    if (fi < 0) {
      fc->rc = HZ_DC_ENOMEM;
      return 1;
    }
    if (fc->bounded) {
      /* Р-8 (вторая редакция раздачи): кусок раздаётся ВСЕМ ЛИСТЬЯМ, чей бокс
       * он пересекает, — спуском по дереву с SAT-отсечением. Раздача только по
       * ячейкам ВЕРШИН полигона теряла клетки, сквозь которые кусок проходит:
       * замерено §778 на room — 99.7 из 102.7 дыр в клетках БЕЗ записи. */
      int32_t rlo[3] = {0, 0, 0};
      if (piece_pairs_rec(fc, 0, rlo, (int32_t)1 << fc->t->log2size, a, b, c, fi)) return 1;
      continue;
    }
    /* Плоскость раздаётся ВСЕМ ячейкам полигона, а не только вершинам этого
     * треугольника: веер вокруг вершины ячейки обязан быть полным (Г10). */
    for (int k = 0; k < nv; k++) {
      int outside = 0;
      if (!plane_cuts_box(nn, off, ref[k].lo, ref[k].size, &outside)) {
        if (outside) fc->nbad++;
        continue;
      }
      if (push_pair(fc, ref[k].ni, fi)) return 1;
    }
  }
  return 0;
}

/* Общий финал раздачи (§794): сорт пар, дедуп, cmap, счётчик Г60. Буфер
 * записи ДИНАМИЧЕСКИЙ (А1280): OBJ-куски дают вееры плотных клеток (листва)
 * много больше HZ_P3_MAXH — предел здесь был бы ETOPO на первой же кроне;
 * предел РЕЗА живёт в cut3 и разводится там (§794 Р2). */
static int facets_finish(facetctx *fc, hz_cutmap *cm) {
  g_fc_pairs = fc->np;
  /* Г45: обход выдаёт ячейки в порядке дерева, hz_cutmap требует возрастания */
  if (fc->np > 1) qsort(fc->pair, (size_t)fc->np, sizeof(cellfacet), pair_cmp);
  int32_t bcap = 256;
  int32_t *buf = malloc((size_t)bcap * sizeof *buf);
  if (buf == NULL) return HZ_DC_ENOMEM;
  int32_t i = 0;
  while (i < fc->np) {
    int32_t cell = fc->pair[i].cell, nf = 0;
    int32_t j = i;
    while (j < fc->np && fc->pair[j].cell == cell) {
      if (j == i || fc->pair[j].fref != fc->pair[j - 1].fref) {
        if (nf >= bcap) {
          bcap *= 2;
          int32_t *nb = realloc(buf, (size_t)bcap * sizeof *nb);
          if (nb == NULL) {
            free(buf);
            return HZ_DC_ENOMEM;
          }
          buf = nb;
        }
        buf[nf++] = fc->pair[j].fref;
      }
      j++;
    }
    if (hz_cutmap_add(cm, cell, buf, nf) != 0) {
      free(buf);
      return HZ_DC_ETOPO;
    }
    i = j;
  }
  free(buf);
  /* Г60: опустевшие после отбора записи — счётчиком, не молчанием */
  if (fc->nd > 0) {
    qsort(fc->dropc, (size_t)fc->nd, sizeof(int32_t), pair_cmp_i32);
    for (int32_t d = 0; d < fc->nd; d++) {
      if (d > 0 && fc->dropc[d] == fc->dropc[d - 1]) continue;
      if (hz_cutmap_find(cm, fc->dropc[d]) == NULL) g_fc_empty++;
    }
  }
  return HZ_DC_OK;
}

int hz_dc_facets2(const hz_dctree *t, hz_dc_stop stop, void *sctx, hz_facettab *ft, hz_cutmap *cm,
                  int bounded) {
  facetctx fc;
  memset(&fc, 0, sizeof fc);
  fc.ft = ft;
  fc.t = t;
  fc.cap = 256;
  fc.bounded = bounded;
  fc.dcap = 256;
  fc.rc = HZ_DC_OK;
  g_fc_pairs = g_fc_dropped = g_fc_empty = 0;
  fc.pair = calloc((size_t)fc.cap, sizeof(cellfacet));
  fc.dropc = calloc((size_t)fc.dcap, sizeof(int32_t));
  if (fc.pair == NULL || fc.dropc == NULL) {
    free(fc.pair);
    free(fc.dropc);
    return HZ_DC_ENOMEM;
  }
  int rc = hz_dc_walk(t, stop, sctx, facet_emit, &fc);
  if (rc == HZ_DC_OK) rc = fc.rc;
  if (rc != HZ_DC_OK) {
    free(fc.pair);
    free(fc.dropc);
    return rc;
  }
  int frc = facets_finish(&fc, cm);
  free(fc.pair);
  free(fc.dropc);
  return frc;
}

int hz_dc_facets(const hz_dctree *t, hz_dc_stop stop, void *sctx, hz_facettab *ft, hz_cutmap *cm) {
  /* Р-8: раздача по КУСКУ — рабочее умолчание; прежнее поведение (бесконечные
   * плоскости) остаётся негативным контролем через hz_dc_facets2(..., 0). */
  return hz_dc_facets2(t, stop, sctx, ft, cm, 1);
}

/* §794: КУСКИ ИЗ ТРЕУГОЛЬНИКОВ СЦЕНЫ. Провайдер отдаёт вершины В ЕДИНИЦАХ
 * кадра (0 — треугольник пропустить); плоскость строится из вершин, кусок —
 * сам треугольник, раздача — тем же SAT-спуском по листьям дна. Ориентация:
 * авторская нормаль наружу материала → материал {n·x <= off} (fref >= 0);
 * двусторонний лист даёт слэб нулевой меры → материал пуст (корректно). */
int hz_dc_facets_tris(const hz_dctree *t, hz_dc_tri_get get, void *gctx, int32_t ntri,
                      hz_facettab *ft, hz_cutmap *cm) {
  facetctx fc;
  memset(&fc, 0, sizeof fc);
  fc.ft = ft;
  fc.t = t;
  fc.cap = 256;
  fc.bounded = 1;
  fc.dcap = 256;
  fc.rc = HZ_DC_OK;
  g_fc_pairs = g_fc_dropped = g_fc_empty = 0;
  fc.pair = calloc((size_t)fc.cap, sizeof(cellfacet));
  fc.dropc = calloc((size_t)fc.dcap, sizeof(int32_t));
  if (fc.pair == NULL || fc.dropc == NULL) {
    free(fc.pair);
    free(fc.dropc);
    return HZ_DC_ENOMEM;
  }
  int rc = HZ_DC_OK;
  for (int32_t i = 0; i < ntri && rc == HZ_DC_OK; i++) {
    double tv[3][3];
    if (!get(gctx, i, tv)) continue;
    double e1[3], e2[3], nn[3];
    for (int k = 0; k < 3; k++) {
      e1[k] = tv[1][k] - tv[0][k];
      e2[k] = tv[2][k] - tv[0][k];
    }
    nn[0] = e1[1] * e2[2] - e1[2] * e2[1];
    nn[1] = e1[2] * e2[0] - e1[0] * e2[2];
    nn[2] = e1[0] * e2[1] - e1[1] * e2[0];
    double mlen = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
    if (!(mlen > 0.0)) continue; /* вырожденный: плоскости нет */
    for (int k = 0; k < 3; k++)
      nn[k] /= mlen;
    double off = nn[0] * tv[0][0] + nn[1] * tv[0][1] + nn[2] * tv[0][2];
    /* dmax = 0 ТОЧНО, а не UNKNOWN: фасет И ЕСТЬ поверхность (§794) —
     * Г44-класс покидает перенос. */
    int32_t fi = hz_facettab_add_units_piece(ft, nn, off, -1, 0.0, (const double (*)[3])tv);
    if (fi < 0) {
      rc = HZ_DC_ENOMEM;
      break;
    }
    int32_t rlo[3] = {0, 0, 0};
    if (piece_pairs_rec(&fc, 0, rlo, (int32_t)1 << t->log2size, tv[0], tv[1], tv[2], fi)) {
      rc = fc.rc;
      break;
    }
  }
  if (rc == HZ_DC_OK) rc = fc.rc;
  if (rc == HZ_DC_OK) rc = facets_finish(&fc, cm);
  free(fc.pair);
  free(fc.dropc);
  return rc;
}

/* §800: КУСКИ-ПОЛИГОНЫ (кластеризация А1285). Провайдер отдаёт выпуклый
 * полигон в ЕДИНИЦАХ кадра (nv вершин, обход согласован с авторской нормалью;
 * 0 — пропустить) и dmax В МЕТРАХ — фактическое отклонение вершин от
 * заявляемой плоскости (А1315: ноль означает «точен»). Плоскость — по Ньюэллу
 * из всех вершин (устойчивее пары рёбер у почти-коллинеарных углов); раздача —
 * тем же SAT-спуском ВЕЕРОМ треугольников полигона: дубликаты пар (клетка,
 * фасет) снимает дедуп facets_finish. */
int hz_dc_facets_polys(const hz_dctree *t, hz_dc_poly_get get, void *gctx, int32_t npoly,
                       hz_facettab *ft, hz_cutmap *cm) {
  facetctx fc;
  memset(&fc, 0, sizeof fc);
  fc.ft = ft;
  fc.t = t;
  fc.cap = 256;
  fc.bounded = 1;
  fc.dcap = 256;
  fc.rc = HZ_DC_OK;
  g_fc_pairs = g_fc_dropped = g_fc_empty = 0;
  fc.pair = calloc((size_t)fc.cap, sizeof(cellfacet));
  fc.dropc = calloc((size_t)fc.dcap, sizeof(int32_t));
  if (fc.pair == NULL || fc.dropc == NULL) {
    free(fc.pair);
    free(fc.dropc);
    return HZ_DC_ENOMEM;
  }
  int rc = HZ_DC_OK;
  for (int32_t i = 0; i < npoly && rc == HZ_DC_OK; i++) {
    double pv[HZ_FACET_TVMAX][3];
    double dmax = 0.0;
    int nv = get(gctx, i, pv, &dmax);
    if (nv < 3 || nv > HZ_FACET_TVMAX) continue;
    /* Ньюэлл: n_x = Σ (y_j − y_k)(z_j + z_k) и циклически */
    double nn[3] = {0, 0, 0};
    for (int j = 0; j < nv; j++) {
      const double *a = pv[j], *b = pv[(j + 1) % nv];
      nn[0] += (a[1] - b[1]) * (a[2] + b[2]);
      nn[1] += (a[2] - b[2]) * (a[0] + b[0]);
      nn[2] += (a[0] - b[0]) * (a[1] + b[1]);
    }
    double mlen = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
    if (!(mlen > 0.0)) continue; /* вырожденный: плоскости нет */
    for (int k = 0; k < 3; k++)
      nn[k] /= mlen;
    double off = nn[0] * pv[0][0] + nn[1] * pv[0][1] + nn[2] * pv[0][2];
    int32_t fi = hz_facettab_add_units_poly(ft, nn, off, -1, dmax, (const double (*)[3])pv, nv);
    if (fi < 0) {
      rc = HZ_DC_ENOMEM;
      break;
    }
    int32_t rlo[3] = {0, 0, 0};
    for (int e = 1; e + 1 < nv && rc == HZ_DC_OK; e++)
      if (piece_pairs_rec(&fc, 0, rlo, (int32_t)1 << t->log2size, pv[0], pv[e], pv[e + 1], fi))
        rc = fc.rc;
    if (rc != HZ_DC_OK) break;
  }
  if (rc == HZ_DC_OK) rc = fc.rc;
  if (rc == HZ_DC_OK) rc = facets_finish(&fc, cm);
  free(fc.pair);
  free(fc.dropc);
  return rc;
}
