/* Сборка полигонов из разметки Ш1. Разбор и оговорки — в `polygon.h`. */

#include "polygon.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* --- точный ключ пары вершин ----------------------------------------------- */

/* Ключ РОВНЫЙ, а не хеш: a·nw + b однозначен при nw < 2^31, и потому у
 * счётчика полурёбер не может быть ложного совпадения. Хешируется только НОМЕР
 * СЛОТА; сам ключ хранится и сравнивается целиком. */
static uint64_t mix64(uint64_t k) {
  k ^= k >> 33;
  k *= 0xFF51AFD7ED558CCDULL;
  k ^= k >> 33;
  k *= 0xC4CEB9FE1A85EC53ULL;
  k ^= k >> 33;
  return k;
}

typedef struct {
  uint64_t *key; /* 0 — пусто, поэтому все ключи сдвинуты на +1 */
  int32_t *val;
  int64_t nslot;
  uint64_t mask;
} hz_ht;

static int ht_init(hz_ht *h, int64_t want) {
  int64_t n = 4;
  while (n < 2 * want + 8)
    n *= 2;
  h->key = calloc((size_t)n, sizeof *h->key);
  h->val = calloc((size_t)n, sizeof *h->val);
  h->nslot = n;
  h->mask = (uint64_t)n - 1;
  return (h->key != NULL && h->val != NULL) ? 0 : 2;
}

static void ht_free(hz_ht *h) {
  free(h->key);
  free(h->val);
  memset(h, 0, sizeof *h);
}

/* Возврат слота; при insert = 0 и отсутствии ключа возвращает −1. */
static int64_t ht_slot(hz_ht *h, uint64_t k, int insert) {
  uint64_t kk = k + 1;
  uint64_t i = mix64(kk) & h->mask;
  for (;;) {
    if (h->key[i] == kk) return (int64_t)i;
    if (h->key[i] == 0) {
      if (!insert) return -1;
      h->key[i] = kk;
      return (int64_t)i;
    }
    i = (i + 1) & h->mask;
  }
}

/* --- сварка вершин по ПОЛОЖЕНИЮ --------------------------------------------- */

static uint64_t pos_key(const double p[3]) {
  uint64_t h = 0xcbf29ce484222325ULL;
  for (int a = 0; a < 3; a++) {
    double x = p[a];
    /* −0.0 и +0.0 равны как числа, но не как биты: одна и та же вершина,
     * записанная с разным знаком нуля, иначе не сварилась бы. */
    if (!(x < 0.0) && !(x > 0.0)) x = 0.0;
    uint64_t b;
    memcpy(&b, &x, sizeof b);
    h = mix64(h ^ b);
  }
  return h;
}

/* --- малая линейная алгебра -------------------------------------------------- */

/* Гаусс с частичным выбором для 3×3 с nrhs правыми частями. 0 — успех. */
static int solve3(double A[3][3], double B[3][3], int nrhs) {
  for (int c = 0; c < 3; c++) {
    int piv = c;
    for (int r = c + 1; r < 3; r++)
      if (fabs(A[r][c]) > fabs(A[piv][c])) piv = r;
    if (!(fabs(A[piv][c]) > 0.0)) return 1;
    if (piv != c) {
      for (int k = 0; k < 3; k++) {
        double t = A[c][k];
        A[c][k] = A[piv][k];
        A[piv][k] = t;
      }
      for (int k = 0; k < nrhs; k++) {
        double t = B[c][k];
        B[c][k] = B[piv][k];
        B[piv][k] = t;
      }
    }
    for (int r = 0; r < 3; r++) {
      if (r == c) continue;
      double f = A[r][c] / A[c][c];
      for (int k = 0; k < 3; k++)
        A[r][k] -= f * A[c][k];
      for (int k = 0; k < nrhs; k++)
        B[r][k] -= f * B[c][k];
    }
  }
  for (int c = 0; c < 3; c++)
    for (int k = 0; k < nrhs; k++)
      B[c][k] /= A[c][c];
  return 0;
}

/* --- край: рёбра участка ------------------------------------------------------ */

typedef struct {
  uint64_t key; /* min·nw + max */
  int32_t a, b; /* НАПРАВЛЕННОЕ ребро */
} pe_edge;

static int cmp_edge(const void *pa, const void *pb) {
  const pe_edge *x = pa, *y = pb;
  if (x->key < y->key) return -1;
  if (x->key > y->key) return 1;
  if (x->a < y->a) return -1;
  if (x->a > y->a) return 1;
  return 0;
}

static int cmp_start(const void *pa, const void *pb) {
  const pe_edge *x = pa, *y = pb;
  if (x->a < y->a) return -1;
  if (x->a > y->a) return 1;
  if (x->b < y->b) return -1;
  if (x->b > y->b) return 1;
  return 0;
}

void hz_poly_world(const hz_poly *p, double u, double v, double x[3]) {
  for (int a = 0; a < 3; a++)
    x[a] = p->org[a] + u * p->eu[a] + v * p->ev[a];
}

/* ПРАВИЛО НЕНУЛЕВОГО ОБОРОТА, А НЕ ЧЁТНО-НЕЧЁТНОЕ. Разница не стилистическая, и
 * она измерена: с чётно-нечётным при δ = 45 мм 5.41% лучей, попавших в
 * треугольник, не попадали ни в один полигон.
 *
 * Причина. Внутренность полигона есть ОБЪЕДИНЕНИЕ проекций его треугольников, а
 * они все ориентированы ПОЛОЖИТЕЛЬНО (сегментатор требует n_тр·n > 0). Значит
 * оборот равен единице там, где покрытие однократно, и ДВУМ там, где двукратно
 * — а двукратное покрытие законно возникает всякий раз, когда две параллельные
 * поверхности, смотрящие в ОДНУ сторону, оказались ближе δ друг к другу (две
 * полки, столешница и подставка). Чётно-нечётное правило объявляет такое место
 * дырой, то есть выбрасывает поверхность, которая там есть. Оборот — не
 * объявляет. Настоящая дыра обходится в ОБРАТНУЮ сторону и даёт ноль при обоих
 * правилах, поэтому ничего не теряется. */
int hz_poly_inside(const hz_polyset *ps, const hz_poly *p, double u, double v) {
  int wind = 0;
  for (int32_t l = p->l0; l < p->l0 + p->nloop; l++) {
    int32_t b = ps->loop[l], e = ps->loop[l + 1];
    int32_t n = e - b;
    for (int32_t i = 0; i < n; i++) {
      const double *A = ps->bv + (size_t)(b + i) * 2;
      const double *B = ps->bv + (size_t)(b + (i + 1) % n) * 2;
      /* Знак площади треугольника (A, B, точка): слева от ребра или справа. */
      double side = (B[0] - A[0]) * (v - A[1]) - (u - A[0]) * (B[1] - A[1]);
      if (A[1] <= v) {
        if (B[1] > v && side > 0.0) wind++;
      } else {
        if (B[1] <= v && side < 0.0) wind--;
      }
    }
  }
  return wind != 0;
}

int hz_poly_weld(const hz_objmesh *m, int32_t *wid, int32_t *nw, double *wpos) {
  hz_ht wt;
  if (ht_init(&wt, m->nv) != 0) {
    ht_free(&wt);
    return 2;
  }
  int32_t n = 0;
  for (int32_t i = 0; i < m->nv; i++) {
    const double *p = m->v + (size_t)i * 3;
    int64_t sl = ht_slot(&wt, pos_key(p), 1);
    if (wt.val[sl] == 0) {
      wt.val[sl] = n + 1;
      if (wpos != NULL)
        for (int a = 0; a < 3; a++)
          wpos[(size_t)n * 3 + (size_t)a] = p[a];
      n++;
    }
    wid[i] = wt.val[sl] - 1;
  }
  ht_free(&wt);
  *nw = n;
  return 0;
}

/* --- опорное множество (О7) ---------------------------------------------------
 *
 * Разбор и обе ветви — в `polygon.h` у поля `sup`. Здесь только устройство.
 * Точка несёт и проекцию на раму полигона (для оболочки), и ИСХОДНУЮ мировую
 * вершину: хранится в опорном множестве именно она, а не восстановленная из
 * `(u,v)`, — иначе `dmax` считался бы по точке, снесённой на плоскость, то есть
 * ровно по той ПРОЕКЦИИ, из-за которой измерение по краю и оказалось незаконным
 * (А7). */
typedef struct {
  double u, v;
  double x[3];
} sup_pt;

/* Порядок ПОЛНЫЙ и лексикографический: сперва проекция (её и требует монотонная
 * цепь), потом координаты — тогда точные дубликаты стоят рядом при любом
 * совпадении проекций. */
static int cmp_sup(const void *a, const void *b) {
  const sup_pt *p = a, *q = b;
  if (p->u < q->u) return -1;
  if (p->u > q->u) return 1;
  if (p->v < q->v) return -1;
  if (p->v > q->v) return 1;
  for (int c = 0; c < 3; c++) {
    if (p->x[c] < q->x[c]) return -1;
    if (p->x[c] > q->x[c]) return 1;
  }
  return 0;
}

static int sup_same(const sup_pt *p, const sup_pt *q) {
  for (int c = 0; c < 3; c++)
    if (p->x[c] < q->x[c] || p->x[c] > q->x[c]) return 0;
  return 1;
}

static double sup_cross(const sup_pt *o, const sup_pt *a, const sup_pt *b) {
  return (a->u - o->u) * (b->v - o->v) - (a->v - o->v) * (b->u - o->u);
}

/* Точка внутри выпуклого многоугольника, заданного против часовой стрелки?
 * Двоичный поиск по клину от `h[0]`: `O(log nh)` на точку, поэтому сторож (А48)
 * стоит `O(n log nh)` на полигон и берётся один раз при импорте.
 * Нестрогие сравнения намеренно: точка НА границе внутри. */
static int sup_inside(const sup_pt *buf, const int32_t *h, int32_t nh, const sup_pt *p) {
  if (nh < 3) return 1; /* вырожденная оболочка: цепь оставила все точки */
  if (sup_cross(&buf[h[0]], &buf[h[1]], p) < 0.0) return 0;
  if (sup_cross(&buf[h[0]], &buf[h[nh - 1]], p) > 0.0) return 0;
  int32_t lo = 1, hi = nh - 1;
  while (hi - lo > 1) {
    int32_t mid = lo + (hi - lo) / 2;
    if (sup_cross(&buf[h[0]], &buf[h[mid]], p) >= 0.0)
      lo = mid;
    else
      hi = mid;
  }
  return sup_cross(&buf[h[lo]], &buf[h[lo + 1]], p) >= 0.0;
}

/* Опорное множество полигона: собрать вершины треугольников, отсортировать,
 * различить, при плоскости ТОЧНО — свернуть до оболочки. Возвращает число
 * оставленных точек, их номера — в `idx` (номера в отсортированном `buf`).
 * `*planar` — по какой ветви пошло; `*nout` растёт на число точек, оказавшихся
 * вне оболочки (сторож А48, обязан остаться нулём). */
static int32_t sup_reduce(const hz_objmesh *m, const hz_polyset *ps, const hz_poly *P, sup_pt *buf,
                          int32_t *idx, int use_hull, int *planar, int *hulled, int64_t *nout) {
  int32_t n = 0;
  int flat = 1;
  *hulled = 0;
  for (int32_t i = 0; i < P->ntri; i++) {
    double p[3][3];
    hz_obj_tri(m, ps->tri[P->t0 + i], p);
    for (int k = 0; k < 3; k++) {
      double q[3];
      for (int a = 0; a < 3; a++)
        q[a] = p[k][a] - P->org[a];
      /* ПЛОСКОСТЬ ПРОВЕРЯЕТСЯ ЗДЕСЬ, А НЕ БЕРЁТСЯ ИЗ `P->dmax` (А44): `dmax`
       * приходит от сегментатора и означает максимум по ЕГО множеству вершин и
       * ЕГО плоскости. Ветвь (А) законна только при равенстве РОВНО нулю у
       * каждой вершины, и это дешевле проверить, чем вывести. */
      double d = p[k][0] * P->n[0] + p[k][1] * P->n[1] + p[k][2] * P->n[2] - P->off;
      if (d < 0.0 || d > 0.0) flat = 0;
      buf[n].u = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
      buf[n].v = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
      for (int a = 0; a < 3; a++)
        buf[n].x[a] = p[k][a];
      n++;
    }
  }
  *planar = flat;
  if (n <= 0) return 0;
  qsort(buf, (size_t)n, sizeof *buf, cmp_sup);
  int32_t nd = 0;
  for (int32_t i = 0; i < n; i++) {
    if (nd > 0 && sup_same(&buf[nd - 1], &buf[i])) continue;
    buf[nd++] = buf[i];
  }
  for (int32_t i = 0; i < nd; i++)
    idx[i] = i;
  if (!use_hull || !flat || nd < 4) return nd;

  /* МОНОТОННАЯ ЦЕПЬ. Точка выбрасывается только при СТРОГОМ знаке: тогда
   * коллинеарные остаются, и результат есть НАДмножество оболочки. Для
   * максимума это безопасно в ту сторону, в какую ошибаться можно, — лишняя
   * точка максимум не портит, потерянная ЗАНИЖАЕТ его (А45, класс А8). */
  int32_t nh = 0;
  for (int32_t i = 0; i < nd; i++) {
    while (nh >= 2 && sup_cross(&buf[idx[nh - 2]], &buf[idx[nh - 1]], &buf[i]) < 0.0)
      nh--;
    idx[nh++] = i;
  }
  int32_t lower = nh + 1;
  for (int32_t i = nd - 2; i >= 0; i--) {
    while (nh >= lower && sup_cross(&buf[idx[nh - 2]], &buf[idx[nh - 1]], &buf[i]) < 0.0)
      nh--;
    idx[nh++] = i;
  }
  nh--; /* последняя точка совпадает с первой */
  if (nh >= nd) {
    for (int32_t i = 0; i < nd; i++)
      idx[i] = i;
    return nd;
  }
  /* СТОРОЖ (А48): каждая точка обязана лежать внутри построенной оболочки.
   * Проверка прямая и от слепка независимая — тот ловит лишь перевернувшиеся
   * решения. Оболочка читается через `idx`, копии не заводится. */
  for (int32_t i = 0; i < nd; i++)
    if (!sup_inside(buf, idx, nh, &buf[i])) (*nout)++;
  *hulled = 1;
  return nh;
}

int hz_poly_init_empty(hz_polyset *ps) {
  memset(ps, 0, sizeof *ps);
  if (hz_facettab_init(&ps->ft) != 0) return 2;
  ps->loop = malloc(sizeof *ps->loop);
  if (ps->loop == NULL) return 2;
  ps->loop[0] = 0;
  return 0;
}

int hz_poly_add_quad(hz_polyset *ps, const double c[3], const double n[3], const double eu[3],
                     double hu, double hv, int32_t mtl) {
  if (!(hu > 0.0) || !(hv > 0.0)) return 1;
  hz_poly *np = realloc(ps->p, (size_t)(ps->np + 1) * sizeof *np);
  if (np == NULL) return 2;
  ps->p = np;
  int32_t *nl = realloc(ps->loop, (size_t)(ps->nloopall + 2) * sizeof *nl);
  if (nl == NULL) return 2;
  ps->loop = nl;
  double *nb = realloc(ps->bv, (size_t)(ps->nbv + 4) * 2 * sizeof *nb);
  if (nb == NULL) return 2;
  ps->bv = nb;
  int32_t *nbw = realloc(ps->bw, (size_t)(ps->nbv + 4) * sizeof *nbw);
  if (nbw == NULL) return 2;
  ps->bw = nbw;

  hz_poly *P = &ps->p[ps->np];
  memset(P, 0, sizeof *P);
  double nn = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  if (!(nn > 0.0)) return 1;
  for (int a = 0; a < 3; a++) {
    P->n[a] = n[a] / nn;
    P->org[a] = c[a];
  }
  double d = eu[0] * P->n[0] + eu[1] * P->n[1] + eu[2] * P->n[2];
  for (int a = 0; a < 3; a++)
    P->eu[a] = eu[a] - d * P->n[a];
  double un = sqrt(P->eu[0] * P->eu[0] + P->eu[1] * P->eu[1] + P->eu[2] * P->eu[2]);
  if (!(un > 0.0)) return 1;
  for (int a = 0; a < 3; a++)
    P->eu[a] /= un;
  P->ev[0] = P->n[1] * P->eu[2] - P->n[2] * P->eu[1];
  P->ev[1] = P->n[2] * P->eu[0] - P->n[0] * P->eu[2];
  P->ev[2] = P->n[0] * P->eu[1] - P->n[1] * P->eu[0];
  P->off = P->n[0] * c[0] + P->n[1] * c[1] + P->n[2] * c[2];
  P->area = 4.0 * hu * hv;
  P->dmax = 0.0; /* здесь ноль ЧЕСТЕН: прямоугольник плоский точно */
  P->mtl = mtl;
  P->seg = -1;
  P->t0 = 0;
  P->ntri = 0;
  P->mom[0] = 4.0 * hu * hv;
  P->mom[1] = 0.0;
  P->mom[2] = 0.0;
  P->mom[3] = P->mom[0] * hu * hu / 3.0;
  P->mom[4] = 0.0;
  P->mom[5] = P->mom[0] * hv * hv / 3.0;
  for (int a = 0; a < 3; a++) {
    P->n0[a] = P->n[a];
    P->nu[a] = 0.0;
    P->nv[a] = 0.0;
  }
  P->uvlo[0] = -hu;
  P->uvhi[0] = hu;
  P->uvlo[1] = -hv;
  P->uvhi[1] = hv;

  const hz_frame idf = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
  P->facet = hz_facettab_add_plane(&ps->ft, &idf, P->n, P->off, -1, 0.0);
  if (P->facet < 0) return 2;

  const double q[4][2] = {{-hu, -hv}, {hu, -hv}, {hu, hv}, {-hu, hv}};
  for (int i = 0; i < 4; i++) {
    ps->bv[(size_t)(ps->nbv + i) * 2 + 0] = q[i][0];
    ps->bv[(size_t)(ps->nbv + i) * 2 + 1] = q[i][1];
    /* Сварного номера у собранного полигона нет: он ни с чем не делит вершин.
     * `−1` запирает его край от упрощения, и это верно по существу — четыре
     * угла прямоугольника не избыточны. */
    ps->bw[ps->nbv + i] = -1;
  }
  P->l0 = ps->nloopall;
  P->nloop = 1;
  ps->nbv += 4;
  ps->nloopall++;
  ps->loop[ps->nloopall] = ps->nbv;

  /* ОПОРНОЕ МНОЖЕСТВО СОБРАННОГО ПОЛИГОНА — ЕГО ЧЕТЫРЕ УГЛА, и ноль здесь был
   * бы ложным. Треугольников у такого полигона нет (`ntri = 0`), поэтому
   * «максимум по треугольникам» дал бы ноль при любой плоскости; четыре угла
   * дают точный максимум, потому что прямоугольник плоский и его оболочка —
   * ровно они. В огрубление такие полигоны сейчас не попадают (`hz_merge`
   * берёт только полигоны из участков), но подразумеваемый ноль — это то, на
   * чём проект уже сидел (Ш7, А4). */
  {
    double *nsp = realloc(ps->sup, (size_t)(ps->nsupall + 4) * 3 * sizeof *nsp);
    if (nsp == NULL) return 2;
    ps->sup = nsp;
    if (ps->nsupall + 4 > INT32_MAX) return 2;
    P->s0 = (int32_t)ps->nsupall;
    P->nsup = 4;
    for (int i = 0; i < 4; i++) {
      double x[3];
      hz_poly_world(P, q[i][0], q[i][1], x);
      for (int a = 0; a < 3; a++)
        ps->sup[(size_t)(ps->nsupall + i) * 3 + (size_t)a] = x[a];
    }
    ps->nsupall += 4;
    ps->nsup_hull++;
  }
  ps->np++;
  return 0;
}

void hz_poly_free(hz_polyset *ps) {
  free(ps->p);
  free(ps->loop);
  free(ps->bv);
  free(ps->bw);
  free(ps->tri);
  free(ps->sup);
  hz_facettab_free(&ps->ft);
  memset(ps, 0, sizeof *ps);
}

/* --- сборка ------------------------------------------------------------------ */

typedef struct {
  int32_t *wid, *poff;
  double *wpos;
  hz_ht he;
  pe_edge *ed, *bd;
  int32_t *used;
  double *marea; /* площадь по материалам: полигон может накрыть шов */
} pb_work;

static void pb_free(pb_work *w) {
  free(w->wid);
  free(w->poff);
  free(w->wpos);
  ht_free(&w->he);
  free(w->ed);
  free(w->bd);
  free(w->used);
  free(w->marea);
}

int hz_poly_build(hz_polyset *ps, const hz_objmesh *m, const hz_pseglist *sg) {
  return hz_poly_build_ex(ps, m, sg, 1);
}

int hz_poly_build_ex(hz_polyset *ps, const hz_objmesh *m, const hz_pseglist *sg, int use_hull) {
  memset(ps, 0, sizeof *ps);
  if (hz_facettab_init(&ps->ft) != 0) return 2;
  const int32_t nt = m->nt, nseg = sg->nseg;
  if (nt <= 0 || nseg <= 0) return 0;

  pb_work w;
  memset(&w, 0, sizeof w);
  w.wid = malloc((size_t)m->nv * sizeof *w.wid);
  w.wpos = malloc((size_t)m->nv * 3 * sizeof *w.wpos);
  w.poff = calloc((size_t)nseg + 1, sizeof *w.poff);
  ps->tri = malloc((size_t)nt * sizeof *ps->tri);
  ps->p = calloc((size_t)nseg, sizeof *ps->p);
  if (w.wid == NULL || w.wpos == NULL || w.poff == NULL || ps->tri == NULL || ps->p == NULL) {
    pb_free(&w);
    return 2;
  }

  /* --- 1. сварка вершин по положению --- */
  {
    int32_t nwv = 0;
    if (hz_poly_weld(m, w.wid, &nwv, w.wpos) != 0) {
      pb_free(&w);
      return 2;
    }
    ps->nweld = nwv;
  }
  const int64_t nw = ps->nweld;

  /* --- 2. полурёбра ВСЕЙ сцены: манифолдность нужна силуэтам Ш6 --- */
  if (ht_init(&w.he, 3 * (int64_t)nt) != 0) {
    pb_free(&w);
    return 2;
  }
  for (int32_t t = 0; t < nt; t++)
    for (int e = 0; e < 3; e++) {
      int32_t a = w.wid[m->f[(size_t)t * 3 + (size_t)e]];
      int32_t b = w.wid[m->f[(size_t)t * 3 + (size_t)((e + 1) % 3)]];
      w.he.val[ht_slot(&w.he, (uint64_t)a * (uint64_t)nw + (uint64_t)b, 1)]++;
    }
  ps->nhe = 3 * (int64_t)nt;
  for (int32_t t = 0; t < nt; t++)
    for (int e = 0; e < 3; e++) {
      int32_t a = w.wid[m->f[(size_t)t * 3 + (size_t)e]];
      int32_t b = w.wid[m->f[(size_t)t * 3 + (size_t)((e + 1) % 3)]];
      int64_t rv = ht_slot(&w.he, (uint64_t)b * (uint64_t)nw + (uint64_t)a, 0);
      if (rv < 0 || w.he.val[rv] == 0) ps->nhe_open++;
      int64_t fw = ht_slot(&w.he, (uint64_t)a * (uint64_t)nw + (uint64_t)b, 0);
      if (fw >= 0 && w.he.val[fw] > 1) ps->nhe_dup++;
    }
  ht_free(&w.he);

  /* --- 3. треугольники по участкам --- */
  for (int32_t t = 0; t < nt; t++)
    w.poff[sg->label[t] + 1]++;
  int32_t maxntri = 0;
  for (int32_t r = 0; r < nseg; r++) {
    if (w.poff[r + 1] > maxntri) maxntri = w.poff[r + 1];
    w.poff[r + 1] += w.poff[r];
  }
  {
    int32_t *cur = malloc((size_t)nseg * sizeof *cur);
    if (cur == NULL) {
      pb_free(&w);
      return 2;
    }
    memcpy(cur, w.poff, (size_t)nseg * sizeof *cur);
    for (int32_t t = 0; t < nt; t++)
      ps->tri[cur[sg->label[t]]++] = t;
    free(cur);
  }

  /* --- 4. буферы под самый большой участок --- */
  int32_t maxe = 3 * maxntri;
  if (maxe < 1) maxe = 1;
  w.ed = malloc((size_t)maxe * sizeof *w.ed);
  w.bd = malloc((size_t)maxe * sizeof *w.bd);
  w.used = malloc((size_t)maxe * sizeof *w.used);
  w.marea = calloc((size_t)(m->nmtl > 0 ? m->nmtl : 1), sizeof *w.marea);
  int32_t *lbuf = malloc(((size_t)(maxe > 0 ? maxe : 1) + 8) * sizeof *lbuf);
  /* Буферы опорного множества — на самый большой участок, как и всё в этом
   * блоке: `3·ntri` точек и столько же номеров. Заводятся один раз, а не на
   * полигон, потому что у города максимум — 695 364 треугольника, и повторное
   * выделение на каждый из 301 430 полигонов стоило бы дороже самой работы. */
  sup_pt *sbuf = malloc((size_t)maxe * sizeof *sbuf);
  int32_t *sidx = malloc((size_t)maxe * sizeof *sidx);
  if (w.ed == NULL || w.bd == NULL || w.used == NULL || w.marea == NULL || lbuf == NULL ||
      sbuf == NULL || sidx == NULL) {
    free(lbuf);
    free(sbuf);
    free(sidx);
    pb_free(&w);
    return 2;
  }
  /* Ёмкость опорного массива растёт удвоением от четверти числа треугольников.
   * Верхняя граница ветви (Б) есть `3·nt` точек, то есть `482` МБ на городе, и
   * выделять их сразу нельзя: при работающей ветви (А) выйдут единицы
   * процентов от этого. Транзит удвоения здесь дёшев ровно потому, что итог
   * мал; на списке пар (О10) он стоил `868` МБ, и там ход обратный. */
  int64_t supcap = (int64_t)nt / 4 + 16;
  ps->sup = malloc((size_t)supcap * 3 * sizeof *ps->sup);
  if (ps->sup == NULL) {
    free(lbuf);
    free(sbuf);
    free(sidx);
    pb_free(&w);
    return 2;
  }
  /* Итоговые массивы края: длина не превосходит числа полурёбер. */
  ps->bv = malloc((size_t)maxe * 2 * sizeof *ps->bv);
  ps->bw = malloc((size_t)maxe * sizeof *ps->bw);
  ps->loop = malloc(((size_t)nt + (size_t)nseg + 2) * sizeof *ps->loop);
  int32_t bvcap = maxe, loopcap = nt + nseg + 1;
  if (ps->bv == NULL || ps->bw == NULL || ps->loop == NULL) {
    free(lbuf);
    free(sbuf);
    free(sidx);
    pb_free(&w);
    return 2;
  }
  ps->loop[0] = 0;
  ps->nloopall = 0;
  ps->nbv = 0;

  /* --- 5. полигон за полигоном --- */
  const hz_frame idf = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
  for (int32_t r = 0; r < nseg; r++) {
    hz_poly *P = &ps->p[r];
    P->seg = r;
    P->t0 = w.poff[r];
    P->ntri = w.poff[r + 1] - w.poff[r];
    /* ПУСТОЙ УЧАСТОК ЗАКОНЕН, и его надо пропустить ЯВНО. Сегментация Ш1
     * пустых не выдаёт, поэтому до огрубления (Ш7) этой ветки не бывало; там
     * же группа может остаться без треугольников, и тогда `tri[t0]` — чтение
     * за концом массива (у последнего участка `t0 == nt`). */
    if (P->ntri <= 0) {
      P->nloop = 0;
      P->l0 = ps->nloopall;
      P->facet = -1;
      /* Опорных точек ноль, и это НЕ ложный ноль: точек у полигона нет вовсе,
       * поэтому в максимум он не вносит ничего — ровно как сейчас, где у него
       * нет треугольников. Счётчик есть, чтобы молчания не было. */
      P->s0 = 0;
      P->nsup = 0;
      ps->nsup_empty++;
      continue;
    }
    P->dmax = sg->seg[r].dmax;
    P->area = sg->seg[r].area;
    for (int a = 0; a < 3; a++)
      P->n[a] = sg->seg[r].n[a];
    P->off = sg->seg[r].off;
    /* Материал — по НАИБОЛЬШЕЙ площади, а не у первого попавшегося
     * треугольника: полигон, накрывший шов двух материалов, иначе получил бы
     * альбедо по случайному признаку. Сколько таких — считается (нужно Ш8). */
    {
      int32_t bm = m->fm[ps->tri[P->t0]];
      double tota = 0.0;
      for (int32_t i = 0; i < P->ntri; i++) {
        int32_t t = ps->tri[P->t0 + i];
        double ai = hz_obj_tri_area(m, t);
        w.marea[m->fm[t]] += ai;
        tota += ai;
      }
      for (int32_t i = 0; i < P->ntri; i++) {
        int32_t t = ps->tri[P->t0 + i];
        if (w.marea[m->fm[t]] > w.marea[bm]) bm = m->fm[t];
      }
      P->mtl = bm;
      if (w.marea[bm] < tota) ps->nmixed++;
      for (int32_t i = 0; i < P->ntri; i++)
        w.marea[m->fm[ps->tri[P->t0 + i]]] = 0.0;
    }

    /* рама: ось, наименее совпадающая с нормалью — выбор детерминированный */
    int ax = 0;
    for (int a = 1; a < 3; a++)
      if (fabs(P->n[a]) < fabs(P->n[ax])) ax = a;
    double e0[3] = {0.0, 0.0, 0.0};
    e0[ax] = 1.0;
    double d = e0[0] * P->n[0] + e0[1] * P->n[1] + e0[2] * P->n[2];
    for (int a = 0; a < 3; a++)
      P->eu[a] = e0[a] - d * P->n[a];
    double en = sqrt(P->eu[0] * P->eu[0] + P->eu[1] * P->eu[1] + P->eu[2] * P->eu[2]);
    for (int a = 0; a < 3; a++)
      P->eu[a] /= en;
    P->ev[0] = P->n[1] * P->eu[2] - P->n[2] * P->eu[1];
    P->ev[1] = P->n[2] * P->eu[0] - P->n[0] * P->eu[2];
    P->ev[2] = P->n[0] * P->eu[1] - P->n[1] * P->eu[0];

    /* начало рамы — центр площади, снесённый на плоскость */
    double cw[3] = {0.0, 0.0, 0.0}, wsum = 0.0;
    for (int32_t i = 0; i < P->ntri; i++) {
      int32_t t = ps->tri[P->t0 + i];
      double p[3][3];
      hz_obj_tri(m, t, p);
      double ai = hz_obj_tri_area(m, t);
      for (int a = 0; a < 3; a++)
        cw[a] += ai * (p[0][a] + p[1][a] + p[2][a]) / 3.0;
      wsum += ai;
    }
    if (wsum > 0.0)
      for (int a = 0; a < 3; a++)
        cw[a] /= wsum;
    double cd = cw[0] * P->n[0] + cw[1] * P->n[1] + cw[2] * P->n[2] - P->off;
    for (int a = 0; a < 3; a++)
      P->org[a] = cw[a] - cd * P->n[a];

    P->facet = hz_facettab_add_plane(&ps->ft, &idf, P->n, P->off, -1, P->dmax);
    if (P->facet < 0) {
      free(lbuf);
      free(sbuf);
      free(sidx);
      pb_free(&w);
      return 2;
    }

    /* --- опорное множество (О7) --- */
    {
      int planar = 0, hulled = 0;
      int32_t ns = sup_reduce(m, ps, P, sbuf, sidx, use_hull, &planar, &hulled, &ps->nsup_out);
      if (planar) ps->nsup_flat++;
      if (hulled)
        ps->nsup_hull++;
      else
        ps->nsup_full++;
      if (ps->nsupall + ns > supcap) {
        int64_t nc = supcap;
        while (nc < ps->nsupall + ns)
          nc *= 2;
        double *nsp = realloc(ps->sup, (size_t)nc * 3 * sizeof *nsp);
        if (nsp == NULL) {
          free(lbuf);
          free(sbuf);
          free(sidx);
          pb_free(&w);
          return 2;
        }
        ps->sup = nsp;
        supcap = nc;
      }
      /* ПЕРЕПОЛНЕНИЕ ИНДЕКСА — ОТКАЗ, А НЕ ЗАВОРОТ (А50). `s0` тридцатидвух-
       * битный по образцу `t0` и `l0`, а точек в ветви (Б) до `3·nt`; на сцене
       * в сотни миллионов треугольников это выйдет за `INT32_MAX` МОЛЧА.
       * Порога здесь нет: граница ТИПА, а не выбранное число. */
      if (ps->nsupall + ns > INT32_MAX) {
        free(lbuf);
        free(sbuf);
        free(sidx);
        pb_free(&w);
        return 2;
      }
      P->s0 = (int32_t)ps->nsupall;
      P->nsup = ns;
      for (int32_t i = 0; i < ns; i++)
        for (int a = 0; a < 3; a++)
          ps->sup[(size_t)(ps->nsupall + i) * 3 + (size_t)a] = sbuf[sidx[i]].x[a];
      ps->nsupall += ns;
    }

    /* моменты и поле нормалей — один проход по треугольникам */
    double G[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}},
           RH[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    memset(P->mom, 0, sizeof P->mom);
    for (int32_t i = 0; i < P->ntri; i++) {
      int32_t t = ps->tri[P->t0 + i];
      double p[3][3], uu[3], vv[3];
      hz_obj_tri(m, t, p);
      for (int k = 0; k < 3; k++) {
        double q[3];
        for (int a = 0; a < 3; a++)
          q[a] = p[k][a] - P->org[a];
        uu[k] = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
        vv[k] = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
      }
      /* ПЛОЩАДЬ БЕРЁТСЯ В ПРОЕКЦИИ, а не в мире: моменты живут в (u,v), и
       * трёхмерная площадь наклонённого на δ треугольника им не соответствует. */
      double A = 0.5 * ((uu[1] - uu[0]) * (vv[2] - vv[0]) - (vv[1] - vv[0]) * (uu[2] - uu[0]));
      double su = uu[0] + uu[1] + uu[2], sv = vv[0] + vv[1] + vv[2];
      P->mom[0] += A;
      P->mom[1] += A * su / 3.0;
      P->mom[2] += A * sv / 3.0;
      P->mom[3] += A * (su * su + uu[0] * uu[0] + uu[1] * uu[1] + uu[2] * uu[2]) / 12.0;
      P->mom[4] += A * (su * sv + uu[0] * vv[0] + uu[1] * vv[1] + uu[2] * vv[2]) / 12.0;
      P->mom[5] += A * (sv * sv + vv[0] * vv[0] + vv[1] * vv[1] + vv[2] * vv[2]) / 12.0;

      double wq = fabs(A) / 3.0;
      for (int k = 0; k < 3; k++) {
        int32_t ni = (m->vn != NULL) ? m->fn[(size_t)t * 3 + (size_t)k] : -1;
        double nk[3];
        if (ni >= 0) {
          for (int a = 0; a < 3; a++)
            nk[a] = m->vn[(size_t)ni * 3 + (size_t)a];
        } else {
          for (int a = 0; a < 3; a++)
            nk[a] = P->n[a];
        }
        double b[3] = {1.0, uu[k], vv[k]};
        for (int i2 = 0; i2 < 3; i2++) {
          for (int j2 = 0; j2 < 3; j2++)
            G[i2][j2] += wq * b[i2] * b[j2];
          for (int a = 0; a < 3; a++)
            RH[i2][a] += wq * b[i2] * nk[a];
        }
      }
    }
    {
      double Gc[3][3], Rc[3][3];
      memcpy(Gc, G, sizeof Gc);
      memcpy(Rc, RH, sizeof Rc);
      if (solve3(Gc, Rc, 3) == 0) {
        for (int a = 0; a < 3; a++) {
          P->n0[a] = Rc[0][a];
          P->nu[a] = Rc[1][a];
          P->nv[a] = Rc[2][a];
        }
      } else {
        for (int a = 0; a < 3; a++) {
          P->n0[a] = P->n[a];
          P->nu[a] = 0.0;
          P->nv[a] = 0.0;
        }
      }
    }

    /* --- край --- */
    int32_t ne = 0;
    for (int32_t i = 0; i < P->ntri; i++) {
      int32_t t = ps->tri[P->t0 + i];
      for (int e = 0; e < 3; e++) {
        int32_t a = w.wid[m->f[(size_t)t * 3 + (size_t)e]];
        int32_t b = w.wid[m->f[(size_t)t * 3 + (size_t)((e + 1) % 3)]];
        w.ed[ne].a = a;
        w.ed[ne].b = b;
        w.ed[ne].key = (a < b) ? (uint64_t)a * (uint64_t)nw + (uint64_t)b
                               : (uint64_t)b * (uint64_t)nw + (uint64_t)a;
        ne++;
      }
    }
    qsort(w.ed, (size_t)ne, sizeof *w.ed, cmp_edge);
    int32_t nb = 0;
    for (int32_t i = 0; i < ne;) {
      int32_t j = i;
      while (j < ne && w.ed[j].key == w.ed[i].key)
        j++;
      /* Внутри группы одинакового РЕБРА: первые — направление a<b, потом a>b
       * (сортировка по `a`). Спарить можно min(f, r) штук; остаток — край. */
      int32_t f = 0;
      for (int32_t k = i; k < j; k++)
        if (w.ed[k].a < w.ed[k].b) f++;
      int32_t rr = (j - i) - f;
      int32_t keepf = (f > rr) ? f - rr : 0;
      int32_t keepr = (rr > f) ? rr - f : 0;
      for (int32_t k = i; k < j; k++) {
        int fwd = (w.ed[k].a < w.ed[k].b);
        if (fwd && keepf > 0) {
          w.bd[nb++] = w.ed[k];
          keepf--;
        } else if (!fwd && keepr > 0) {
          w.bd[nb++] = w.ed[k];
          keepr--;
        }
      }
      i = j;
    }
    qsort(w.bd, (size_t)nb, sizeof *w.bd, cmp_start);
    memset(w.used, 0, (size_t)nb * sizeof *w.used);

    P->l0 = ps->nloopall;
    P->nloop = 0;
    int64_t drop0 = ps->nloop_drop;
    for (int32_t st = 0; st < nb; st++) {
      if (w.used[st]) continue;
      int32_t start = w.bd[st].a, cur = w.bd[st].b, prev = start;
      int32_t nvloop = 0;
      w.used[st] = 1;
      lbuf[nvloop++] = start;
      int ltrunc = 0;
      for (;;) {
        if (cur == start) break;
        if (nvloop > maxe) { /* петля длиннее числа рёбер невозможна */
          ltrunc = 1;
          break;
        }
        lbuf[nvloop++] = cur;
        /* следующее неиспользованное ребро, начинающееся в cur */
        int32_t lo = 0, hi = nb;
        while (lo < hi) {
          int32_t mid = lo + (hi - lo) / 2;
          if (w.bd[mid].a < cur)
            lo = mid + 1;
          else
            hi = mid;
        }
        int32_t best = -1, cand = 0;
        double bestang = 0.0;
        double du = 0.0, dv = 0.0;
        {
          const double *A = w.wpos + (size_t)prev * 3, *B = w.wpos + (size_t)cur * 3;
          double q[3];
          for (int a = 0; a < 3; a++)
            q[a] = B[a] - A[a];
          du = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
          dv = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
        }
        for (int32_t k = lo; k < nb && w.bd[k].a == cur; k++) {
          if (w.used[k]) continue;
          cand++;
          const double *B = w.wpos + (size_t)cur * 3, *C = w.wpos + (size_t)w.bd[k].b * 3;
          double q[3];
          for (int a = 0; a < 3; a++)
            q[a] = C[a] - B[a];
          double ou = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
          double ov = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
          /* САМЫЙ ПРАВЫЙ ПОВОРОТ: внутренность слева, значит на развилке надо
           * прижиматься к ней, а не уходить в соседнюю петлю. */
          double ang = atan2(du * ov - dv * ou, du * ou + dv * ov);
          if (best < 0 || ang < bestang) {
            best = k;
            bestang = ang;
          }
        }
        if (cand > 1) ps->nfork++;
        if (best < 0) {
          ps->nopen++;
          break; /* петля замыкается насильно: край всё равно замкнутая ломаная */
        }
        w.used[best] = 1;
        prev = cur;
        cur = w.bd[best].b;
      }
      if (ltrunc) ps->nloop_trunc++;
      /* ОТБРОШЕННАЯ ПЕТЛЯ — СЧИТАЕТСЯ. Полигон при этом может СОХРАНИТЬ другие
       * петли, то есть остаться с НЕПОЛНЫМ краем, не попав в счёт «без края»
       * вовсе. Разница существенна для всего, что меряет поверхность по краю. */
      if (nvloop < 3) {
        ps->nloop_drop++;
        continue;
      }
      if (ps->nbv + nvloop > bvcap) {
        int32_t nc = (bvcap * 2 > ps->nbv + nvloop) ? bvcap * 2 : ps->nbv + nvloop;
        double *nb2 = realloc(ps->bv, (size_t)nc * 2 * sizeof *ps->bv);
        if (nb2 == NULL) {
          free(lbuf);
          free(sbuf);
          free(sidx);
          pb_free(&w);
          return 2;
        }
        ps->bv = nb2;
        int32_t *nw2 = realloc(ps->bw, (size_t)nc * sizeof *ps->bw);
        if (nw2 == NULL) {
          free(lbuf);
          free(sbuf);
          free(sidx);
          pb_free(&w);
          return 2;
        }
        ps->bw = nw2;
        bvcap = nc;
      }
      if (ps->nloopall + 1 >= loopcap) {
        int32_t nc = loopcap * 2;
        int32_t *nl = realloc(ps->loop, ((size_t)nc + 2) * sizeof *ps->loop);
        if (nl == NULL) {
          free(lbuf);
          free(sbuf);
          free(sidx);
          pb_free(&w);
          return 2;
        }
        ps->loop = nl;
        loopcap = nc;
      }
      for (int32_t i2 = 0; i2 < nvloop; i2++) {
        const double *X = w.wpos + (size_t)lbuf[i2] * 3;
        /* СВАРНОЙ НОМЕР СОХРАНЯЕТСЯ, а не выбрасывается вместе с `lbuf`. По
         * (u,v) соседние полигоны узнать общую вершину не могут — у каждого
         * своя рама, — а согласованное упрощение края без этого невозможно. */
        ps->bw[ps->nbv + i2] = lbuf[i2];
        double q[3];
        for (int a = 0; a < 3; a++)
          q[a] = X[a] - P->org[a];
        ps->bv[(size_t)(ps->nbv + i2) * 2 + 0] =
            q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
        ps->bv[(size_t)(ps->nbv + i2) * 2 + 1] =
            q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
      }
      ps->nbv += nvloop;
      ps->nloopall++;
      ps->loop[ps->nloopall] = ps->nbv;
      P->nloop++;
    }
    /* КРАЙ НЕПОЛОН, А НЕ ОТСУТСТВУЕТ: часть петель отброшена, часть осталась.
     * Именно этот класс не покрывался счётом «полигонов без края». */
    if (ps->nloop_drop > drop0 && P->nloop > 0) ps->nloop_part++;
  }
  /* габарит края в (u,v) — дешёвый отсев для луча и растеризации */
  for (int32_t r = 0; r < nseg; r++) {
    hz_poly *P = &ps->p[r];
    P->uvlo[0] = P->uvlo[1] = 1e300;
    P->uvhi[0] = P->uvhi[1] = -1e300;
    int32_t b = (P->nloop > 0) ? ps->loop[P->l0] : 0;
    int32_t e = (P->nloop > 0) ? ps->loop[P->l0 + P->nloop] : 0;
    for (int32_t i = b; i < e; i++)
      for (int a = 0; a < 2; a++) {
        double x = ps->bv[(size_t)i * 2 + (size_t)a];
        if (x < P->uvlo[a]) P->uvlo[a] = x;
        if (x > P->uvhi[a]) P->uvhi[a] = x;
      }
    if (!(P->uvlo[0] <= P->uvhi[0])) {
      P->uvlo[0] = P->uvlo[1] = 0.0;
      P->uvhi[0] = P->uvhi[1] = 0.0;
    }
  }
  ps->np = nseg;
  free(lbuf);
  free(sbuf);
  free(sidx);
  pb_free(&w);
  return 0;
}

/* Решение системы 3×3 методом Гаусса с выбором главного элемента. Заведено для
 * подгонки линейного поля по элементу (§98): та же система, что у моментов. */
int hz_solve3x3(double G[3][3], const double r[3], double out[3]) {
  double A[3][4];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++)
      A[i][j] = G[i][j];
    A[i][3] = r[i];
  }
  for (int c = 0; c < 3; c++) {
    int p = c;
    for (int i = c + 1; i < 3; i++)
      if (fabs(A[i][c]) > fabs(A[p][c])) p = i;
    if (!(fabs(A[p][c]) > 1e-300)) return 1;
    if (p != c)
      for (int j = 0; j < 4; j++) {
        double t = A[c][j];
        A[c][j] = A[p][j];
        A[p][j] = t;
      }
    for (int i = 0; i < 3; i++) {
      if (i == c) continue;
      double f = A[i][c] / A[c][c];
      for (int j = c; j < 4; j++)
        A[i][j] -= f * A[c][j];
    }
  }
  for (int i = 0; i < 3; i++)
    out[i] = A[i][3] / A[i][i];
  return 0;
}
