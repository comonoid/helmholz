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
  }
  P->l0 = ps->nloopall;
  P->nloop = 1;
  ps->nbv += 4;
  ps->nloopall++;
  ps->loop[ps->nloopall] = ps->nbv;
  ps->np++;
  return 0;
}

void hz_poly_free(hz_polyset *ps) {
  free(ps->p);
  free(ps->loop);
  free(ps->bv);
  free(ps->tri);
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
    hz_ht wt;
    if (ht_init(&wt, m->nv) != 0) {
      ht_free(&wt);
      pb_free(&w);
      return 2;
    }
    int32_t nw = 0;
    for (int32_t i = 0; i < m->nv; i++) {
      const double *p = m->v + (size_t)i * 3;
      int64_t sl = ht_slot(&wt, pos_key(p), 1);
      if (wt.val[sl] == 0) {
        wt.val[sl] = nw + 1;
        for (int a = 0; a < 3; a++)
          w.wpos[(size_t)nw * 3 + (size_t)a] = p[a];
        nw++;
      }
      w.wid[i] = wt.val[sl] - 1;
    }
    ht_free(&wt);
    ps->nweld = nw;
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
  if (w.ed == NULL || w.bd == NULL || w.used == NULL || w.marea == NULL || lbuf == NULL) {
    free(lbuf);
    pb_free(&w);
    return 2;
  }
  /* Итоговые массивы края: длина не превосходит числа полурёбер. */
  ps->bv = malloc((size_t)maxe * 2 * sizeof *ps->bv);
  ps->loop = malloc(((size_t)nt + (size_t)nseg + 2) * sizeof *ps->loop);
  int32_t bvcap = maxe, loopcap = nt + nseg + 1;
  if (ps->bv == NULL || ps->loop == NULL) {
    free(lbuf);
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
      pb_free(&w);
      return 2;
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
    for (int32_t st = 0; st < nb; st++) {
      if (w.used[st]) continue;
      int32_t start = w.bd[st].a, cur = w.bd[st].b, prev = start;
      int32_t nvloop = 0;
      w.used[st] = 1;
      lbuf[nvloop++] = start;
      for (;;) {
        if (cur == start) break;
        if (nvloop > maxe) break; /* петля длиннее числа рёбер невозможна */
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
      if (nvloop < 3) continue;
      if (ps->nbv + nvloop > bvcap) {
        int32_t nc = (bvcap * 2 > ps->nbv + nvloop) ? bvcap * 2 : ps->nbv + nvloop;
        double *nb2 = realloc(ps->bv, (size_t)nc * 2 * sizeof *ps->bv);
        if (nb2 == NULL) {
          free(lbuf);
          pb_free(&w);
          return 2;
        }
        ps->bv = nb2;
        bvcap = nc;
      }
      if (ps->nloopall + 1 >= loopcap) {
        int32_t nc = loopcap * 2;
        int32_t *nl = realloc(ps->loop, ((size_t)nc + 2) * sizeof *ps->loop);
        if (nl == NULL) {
          free(lbuf);
          pb_free(&w);
          return 2;
        }
        ps->loop = nl;
        loopcap = nc;
      }
      for (int32_t i2 = 0; i2 < nvloop; i2++) {
        const double *X = w.wpos + (size_t)lbuf[i2] * 3;
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
  pb_free(&w);
  return 0;
}
