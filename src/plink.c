/* plink.c — сборка и решение оператора связями. Разбор — в `plink.h`. */
#include "plink.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define LK_PI 3.14159265358979323846

/* --- списки детей по лестнице ---------------------------------------------- */

/* Лестница хранит РОДИТЕЛЯ у каждого узла; для дробления нужны ДЕТИ. Строятся
 * счётной сортировкой: один проход посчитать, второй разложить. */
typedef struct {
  int32_t *start; /* nnd+1 */
  int32_t *ch;    /* дети подряд */
} lk_kids;

static int lk_kids_build(lk_kids *K, const hz_lod *L) {
  K->start = calloc((size_t)L->nnd + 1, sizeof *K->start);
  if (K->start == NULL) return 1;
  int32_t nch = 0;
  for (int32_t k = 0; k < L->nnd; k++) {
    int32_t p = L->nd[k].parent;
    if (p >= 0 && p < L->nnd) {
      K->start[p + 1]++;
      nch++;
    }
  }
  for (int32_t k = 0; k < L->nnd; k++)
    K->start[k + 1] += K->start[k];
  K->ch = malloc((size_t)(nch > 0 ? nch : 1) * sizeof *K->ch);
  if (K->ch == NULL) {
    free(K->start);
    return 1;
  }
  int32_t *pos = malloc((size_t)L->nnd * sizeof *pos);
  if (pos == NULL) {
    free(K->start);
    free(K->ch);
    return 1;
  }
  memcpy(pos, K->start, (size_t)L->nnd * sizeof *pos);
  for (int32_t k = 0; k < L->nnd; k++) {
    int32_t p = L->nd[k].parent;
    if (p >= 0 && p < L->nnd) K->ch[pos[p]++] = k;
  }
  free(pos);
  return 0;
}

static void lk_kids_free(lk_kids *K) {
  free(K->start);
  free(K->ch);
  memset(K, 0, sizeof *K);
}

/* --- геометрия связи -------------------------------------------------------- */

/* КОЭФФИЦИЕНТ ТОЧКА-ДИСК. Внутренний узел лестницы полигона не несёт — у него есть
 * опорная точка, нормаль и площадь, — поэтому излучатель приближается диском той же
 * площади. Слагаемое `A_j` в знаменателе не подгонка, а регуляризация: без него
 * выражение расходится при сближении, а с ним при `r² ≫ A_j` оно переходит в точную
 * формулу, при `r → 0` остаётся конечным. Ошибка приближения ограничена ровно тем
 * критерием, по которому идёт дробление: чем меньше угловой размер, тем точнее диск.
 *
 *     F = cosθ_i · cosθ_j · A_j / (π r² + A_j)
 */
static double lk_coef(const hz_lodnode *A, const hz_lodnode *B) {
  double d[3] = {B->cx - A->cx, B->cy - A->cy, B->cz - A->cz};
  double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  if (!(r2 > 0.0)) return 0.0;
  double r = sqrt(r2);
  double ci = (d[0] * A->n[0] + d[1] * A->n[1] + d[2] * A->n[2]) / r;
  double cj = -(d[0] * B->n[0] + d[1] * B->n[1] + d[2] * B->n[2]) / r;
  if (!(ci > 0.0) || !(cj > 0.0)) return 0.0; /* друг друга не видят по ориентации */
  return ci * cj * B->area_surf / (LK_PI * r2 + B->area_surf);
}

/* Угловой размер `B` из `A` — та же величина, что в срезе камеры, только вместо
 * глаза стоит элемент. Один допуск на всю схему. */
static double lk_subtend(const hz_lodnode *A, const hz_lodnode *B) {
  double d[3] = {B->cx - A->cx, B->cy - A->cy, B->cz - A->cz};
  double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  if (!(r2 > 0.0)) return 1e300;
  return sqrt(B->area_surf / r2);
}

/* ДОЛЯ НЕЗАСЛОНЁННОГО. Пробы ставятся на дисках узлов: точка сдвигается в
 * касательной плоскости на долю радиуса объемлющей сферы. Это грубее, чем пробы по
 * краю полигона, и намеренно: у внутреннего узла края нет вовсе, а у листа он есть,
 * но связь и так есть осреднение по паре площадей. */
static double lk_vis(const hz_scene *sc, const hz_lodnode *A, const hz_lodnode *B, int nvis,
                     int64_t *nray) {
  double eu[3], ev[3];
  int ax = 0;
  for (int c = 1; c < 3; c++)
    if (fabs(A->n[c]) < fabs(A->n[ax])) ax = c;
  double t0[3] = {0.0, 0.0, 0.0};
  t0[ax] = 1.0;
  eu[0] = A->n[1] * t0[2] - A->n[2] * t0[1];
  eu[1] = A->n[2] * t0[0] - A->n[0] * t0[2];
  eu[2] = A->n[0] * t0[1] - A->n[1] * t0[0];
  double en = sqrt(eu[0] * eu[0] + eu[1] * eu[1] + eu[2] * eu[2]);
  if (!(en > 0.0)) return 0.0;
  for (int c = 0; c < 3; c++)
    eu[c] /= en;
  ev[0] = A->n[1] * eu[2] - A->n[2] * eu[1];
  ev[1] = A->n[2] * eu[0] - A->n[0] * eu[2];
  ev[2] = A->n[0] * eu[1] - A->n[1] * eu[0];
  int hit = 0, tot = 0;
  for (int u = 0; u < nvis; u++)
    for (int v = 0; v < nvis; v++) {
      double su = (((double)u + 0.5) / nvis - 0.5) * 1.2 * A->rad;
      double sv = (((double)v + 0.5) / nvis - 0.5) * 1.2 * A->rad;
      double tu = (((double)v + 0.5) / nvis - 0.5) * 1.2 * B->rad;
      double tv = (((double)u + 0.5) / nvis - 0.5) * 1.2 * B->rad;
      double xa[3] = {A->cx + su * eu[0] + sv * ev[0], A->cy + su * eu[1] + sv * ev[1],
                      A->cz + su * eu[2] + sv * ev[2]};
      double xb[3] = {B->cx + tu * eu[0] + tv * ev[0], B->cy + tu * eu[1] + tv * ev[1],
                      B->cz + tu * eu[2] + tv * ev[2]};
      double d[3], len = 0.0;
      for (int c = 0; c < 3; c++) {
        d[c] = xb[c] - xa[c];
        len += d[c] * d[c];
      }
      len = sqrt(len);
      if (!(len > 0.0)) continue;
      double eps = 1e-6 * (1.0 + len), xo[3];
      for (int c = 0; c < 3; c++) {
        d[c] /= len;
        xo[c] = xa[c] + eps * A->n[c];
      }
      tot++;
      (*nray)++;
      if (!hz_pray_occluded(sc->g, xo, d, 0.0, len * (1.0 - 1e-5))) hit++;
    }
  return (tot > 0) ? (double)hit / (double)tot : 0.0;
}

/* --- сборка ----------------------------------------------------------------- */

typedef struct {
  int32_t a, b;
} lk_pair;

typedef struct {
  lk_pair *p;
  int64_t n, cap;
} lk_queue;

static int lk_qpush(lk_queue *q, int32_t a, int32_t b) {
  if (q->n >= q->cap) {
    int64_t nc = (q->cap > 0) ? q->cap * 2 : 4096;
    lk_pair *np = realloc(q->p, (size_t)nc * sizeof *np);
    if (np == NULL) return 1;
    q->p = np;
    q->cap = nc;
  }
  q->p[q->n].a = a;
  q->p[q->n].b = b;
  q->n++;
  return 0;
}

static int lk_push(hz_linkset *S, int32_t i, int32_t j, double f) {
  if (S->n >= S->cap) {
    int64_t nc = (S->cap > 0) ? S->cap * 2 : 4096;
    hz_link *nl = realloc(S->l, (size_t)nc * sizeof *nl);
    if (nl == NULL) return 1;
    S->l = nl;
    S->cap = nc;
  }
  S->l[S->n].i = i;
  S->l[S->n].j = j;
  S->l[S->n].f = f;
  S->n++;
  return 0;
}

/* --- подпись замкнутости `Σf` ------------------------------------------------ */

typedef struct {
  double sf, ar;
} lk_sfa;

static int lk_cmp_sf(const void *a, const void *b) {
  double x = ((const lk_sfa *)a)->sf, y = ((const lk_sfa *)b)->sf;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

/* `Σ_j f_ij` по ЛИСТУ с учётом предков: грубая связь приносит энергию всем
 * потомкам, значит и в сумму листа входит. Проценты — по ПЛОЩАДИ (А194). */
static void lk_closure(hz_linkset *S, const hz_lod *L) {
  const int32_t nn = L->nnd;
  double *sf = calloc((size_t)(nn > 0 ? nn : 1), sizeof *sf);
  lk_sfa *v = malloc((size_t)(nn > 0 ? nn : 1) * sizeof *v);
  if (sf == NULL || v == NULL) {
    free(sf);
    free(v);
    return;
  }
  for (int64_t l = 0; l < S->n; l++)
    sf[S->l[l].i] += S->l[l].f;
  for (int32_t lev = L->nlev - 1; lev >= 1; lev--)
    for (int32_t k = 0; k < nn; k++) {
      if (L->nd[k].level != lev - 1) continue;
      int32_t p = L->nd[k].parent;
      if (p >= 0 && p < nn) sf[k] += sf[p];
    }
  int32_t m = 0;
  double atot = 0.0, sum = 0.0, alow = 0.0;
  S->sf_min = 1e300;
  S->sf_max = 0.0;
  for (int32_t k = 0; k < nn; k++) {
    if (L->nd[k].level != 0) continue;
    double a = L->nd[k].area_surf;
    v[m].sf = sf[k];
    v[m].ar = a;
    m++;
    atot += a;
    sum += sf[k] * a;
    if (sf[k] < 0.9) alow += a;
    if (sf[k] < S->sf_min) S->sf_min = sf[k];
    if (sf[k] > S->sf_max) S->sf_max = sf[k];
  }
  if (m == 0 || !(atot > 0.0)) {
    S->sf_min = 0.0;
    free(sf);
    free(v);
    return;
  }
  S->sf_mean = sum / atot;
  S->sf_lowfrac = alow / atot;
  qsort(v, (size_t)m, sizeof *v, lk_cmp_sf);
  const double qs[3] = {0.10, 0.50, 0.90};
  double *out[3] = {&S->sf_p10, &S->sf_p50, &S->sf_p90};
  int qi = 0;
  double cum = 0.0;
  for (int32_t k = 0; k < m && qi < 3; k++) {
    cum += v[k].ar;
    while (qi < 3 && cum >= qs[qi] * atot) {
      *out[qi] = v[k].sf;
      qi++;
    }
  }
  while (qi < 3) {
    *out[qi] = v[m - 1].sf;
    qi++;
  }
  free(sf);
  free(v);
}

int hz_links_build(hz_linkset *S, const hz_scene *sc, double eps, int nvis, int noself) {
  memset(S, 0, sizeof *S);
  if (nvis < 1) nvis = 2;
  const hz_lod *L = sc->L;
  lk_kids K;
  if (lk_kids_build(&K, L) != 0) return 2;

  /* КОРНИ — узлы САМОГО ГРУБОГО уровня, реально занятые. Иерархия начинается
   * сверху и дробится вниз: пара грубых узлов либо порождает связь, либо
   * распадается на пары их детей. Именно так число связей и выходит `O(n log n)`
   * вместо `O(n²)` — дальние пары остаются наверху. */
  int32_t top = L->nlev - 1;
  int32_t *isroot = calloc((size_t)L->nnd, sizeof *isroot);
  if (isroot == NULL) {
    lk_kids_free(&K);
    return 2;
  }
  for (int32_t k = 0; k < L->np; k++)
    isroot[L->lab[(size_t)top * (size_t)L->np + (size_t)k]] = 1;

  lk_queue q;
  memset(&q, 0, sizeof q);
  for (int32_t a = 0; a < L->nnd; a++) {
    if (!isroot[a]) continue;
    for (int32_t b = 0; b < L->nnd; b++)
      if (isroot[b] && (a != b || !noself)) {
        if (lk_qpush(&q, a, b) != 0) {
          free(isroot);
          free(q.p);
          lk_kids_free(&K);
          return 2;
        }
      }
  }
  free(isroot);

  int64_t *cnt = calloc((size_t)L->nnd, sizeof *cnt);
  if (cnt == NULL) {
    free(q.p);
    lk_kids_free(&K);
    return 2;
  }

  for (int64_t k = 0; k < q.n; k++) {
    int32_t ia = q.p[k].a, ib = q.p[k].b;
    S->nvisit++;
    if (ia == ib) {
      /* САМОПАРА — обмен ВНУТРИ узла, то есть между его потомками. Раскрывается
       * во ВСЕ пары детей (включая их собственные самопары), а не по правилу
       * «дробить большего»: см. разбор в `plink.h`. Лист сам себя не видит. */
      int32_t s0 = K.start[ia], s1 = K.start[ia + 1];
      if (s1 <= s0) {
        S->nzero++;
        continue;
      }
      S->nrefine++;
      for (int32_t u = s0; u < s1; u++)
        for (int32_t v = s0; v < s1; v++)
          if (lk_qpush(&q, K.ch[u], K.ch[v]) != 0) {
            free(q.p);
            free(cnt);
            lk_kids_free(&K);
            return 2;
          }
      continue;
    }
    const hz_lodnode *A = &L->nd[ia], *B = &L->nd[ib];
    double f = lk_coef(A, B);
    if (!(f > 0.0)) { /* отвёрнуты друг от друга — связи нет по ориентации */
      S->nzero++;
      if (L->nd[ia].level > 0 || L->nd[ib].level > 0) S->nzero_coarse++;
      continue;
    }
    /* ДРОБИТЬ ИЛИ ЗАВОДИТЬ. Дробится тот, чей угловой размер больше: у пары
     * «стена и стул» дробить надо стену, а не стул. Если у него нет детей —
     * дробить нечем, связь заводится как есть. */
    double sa = lk_subtend(B, A), sb = lk_subtend(A, B);
    int32_t big = (sa > sb) ? ia : ib;
    double smax = (sa > sb) ? sa : sb;
    int32_t c0 = K.start[big], c1 = K.start[big + 1];
    if (smax > eps && c1 > c0) {
      S->nrefine++;
      for (int32_t c = c0; c < c1; c++) {
        int32_t ch = K.ch[c];
        int rc = (big == ia) ? lk_qpush(&q, ch, ib) : lk_qpush(&q, ia, ch);
        if (rc != 0) {
          free(q.p);
          free(cnt);
          lk_kids_free(&K);
          return 2;
        }
      }
      continue;
    }
    double v = lk_vis(sc, A, B, nvis, &S->nray);
    if (!(v > 0.0)) {
      S->nzero++;
      continue;
    }
    if (lk_push(S, ia, ib, f * v) != 0) {
      free(q.p);
      free(cnt);
      lk_kids_free(&K);
      return 2;
    }
    cnt[ia]++;
  }
  for (int32_t k = 0; k < L->nnd; k++)
    if (cnt[k] > S->nmax_node) S->nmax_node = cnt[k];
  free(cnt);
  free(q.p);
  lk_kids_free(&K);
  lk_closure(S, L);
  return 0;
}

void hz_links_free(hz_linkset *S) {
  free(S->l);
  memset(S, 0, sizeof *S);
}

/* --- решение ---------------------------------------------------------------- */

int hz_links_solve(const hz_linkset *S, const hz_lod *L, const double *Le, const double *rho,
                   double *B, double tol, int maxit, double *reshist, int nopull) {
  const int32_t nn = L->nnd;
  double *acc = malloc((size_t)nn * sizeof *acc);
  double *ar = malloc((size_t)nn * sizeof *ar);
  double *Bp = malloc((size_t)nn * sizeof *Bp); /* снимок начала итерации */
  double *sw = calloc((size_t)nn, sizeof *sw);
  double *sa = calloc((size_t)nn, sizeof *sa);
  if (acc == NULL || ar == NULL || Bp == NULL || sw == NULL || sa == NULL) {
    free(acc);
    free(ar);
    free(Bp);
    free(sw);
    free(sa);
    return -1;
  }
  /* Площадь узла для взвешивания при подъёме. У узла она уже посчитана. */
  for (int32_t k = 0; k < nn; k++) {
    ar[k] = L->nd[k].area_surf;
    B[k] = Le[k];
  }
  int it = 0;
  for (; it < maxit; it++) {
    memcpy(Bp, B, (size_t)nn * sizeof *Bp);
    for (int32_t k = 0; k < nn; k++)
      acc[k] = 0.0;
    /* 1. СБОР: одна итерация — один порядок отражения (ряд Неймана). */
    for (int64_t l = 0; l < S->n; l++) {
      const hz_link *K = &S->l[l];
      acc[K->i] += K->f * B[K->j];
    }
    /* 2. ВНИЗ: облучённость предка наследуется потомками. Обход от грубых уровней
     * к мелким, поэтому вклад проходит всю цепочку за один проход. */
    for (int lev = L->nlev - 1; lev >= 1; lev--)
      for (int32_t k = 0; k < nn; k++) {
        if (L->nd[k].level != lev - 1) continue;
        int32_t p = L->nd[k].parent;
        if (p >= 0 && p < nn) acc[k] += acc[p];
      }
    /* 3. ЗАТЕНЕНИЕ. `acc` есть Σ форм-факторов на радианс, то есть облучённость,
     * делённая на `π`; радианс отражённого потому `ρ·acc`, БЕЗ второго деления на
     * `π`. Прежняя редакция делила ещё раз и понижала альбедо с `0.5` до `0.159`
     * (А195). */
    for (int32_t k = 0; k < nn; k++)
      B[k] = Le[k] + rho[k] * acc[k];
    /* 4. ВВЕРХ: радианс предка — средневзвешенное по площади радиансов детей.
     * Обход от мелких к грубым, чтобы дети были готовы раньше родителя. */
    if (!nopull)
      for (int lev = 1; lev < L->nlev; lev++) {
        memset(sw, 0, (size_t)nn * sizeof *sw);
        memset(sa, 0, (size_t)nn * sizeof *sa);
        for (int32_t k = 0; k < nn; k++) {
          if (L->nd[k].level != lev - 1) continue;
          int32_t p = L->nd[k].parent;
          if (p >= 0 && p < nn) {
            sw[p] += B[k] * ar[k];
            sa[p] += ar[k];
          }
        }
        for (int32_t k = 0; k < nn; k++)
          if (L->nd[k].level == lev && sa[k] > 0.0) B[k] = sw[k] / sa[k];
      }
    /* НЕВЯЗКА — максимум относительного изменения за ПОЛНУЮ итерацию. Ноль на
     * неподвижной точке по построению; разбор и оговорка — в `plink.h`. */
    double dmx = 0.0, bmx = 0.0;
    for (int32_t k = 0; k < nn; k++) {
      double d = fabs(B[k] - Bp[k]), b = fabs(B[k]);
      if (d > dmx) dmx = d;
      if (b > bmx) bmx = b;
    }
    double res = (bmx > 0.0) ? dmx / bmx : 0.0;
    if (reshist != NULL) reshist[it] = res;
    if (res < tol) {
      it++;
      break;
    }
  }
  free(acc);
  free(ar);
  free(Bp);
  free(sw);
  free(sa);
  return it;
}
