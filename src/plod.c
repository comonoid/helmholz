/* Иерархическое огрубление и срез LOD. Разбор и оговорки — в `plod.h`. */

#include "plod.h"
#include "pmerge.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* --- малая линейная алгебра и минимаксная плоскость --------------------------
 * Дублирует по устройству решатель из `pmerge.c` и делает это СОЗНАТЕЛЬНО: там
 * он сплетён со счётчиками ворот и с конусом, и вытаскивать его наружу — правка
 * горячего кода огрубителя, которую этот шаг не заказывал. Долг записан в §58. */

static int pl_solve(double A[4][4], double b[4], double x[4], int n) {
  if (n < 1 || n > 4) return 1;
  for (int c = 0; c < n; c++) {
    int piv = c;
    for (int r = c + 1; r < n; r++)
      if (fabs(A[r][c]) > fabs(A[piv][c])) piv = r;
    if (!(fabs(A[piv][c]) > 0.0)) return 1;
    if (piv != c) {
      for (int k = 0; k < n; k++) {
        double t = A[c][k];
        A[c][k] = A[piv][k];
        A[piv][k] = t;
      }
      double t = b[c];
      b[c] = b[piv];
      b[piv] = t;
    }
    for (int r = 0; r < n; r++) {
      if (r == c) continue;
      double f = A[r][c] / A[c][c];
      for (int k = 0; k < n; k++)
        A[r][k] -= f * A[c][k];
      b[r] -= f * b[c];
    }
  }
  for (int c = 0; c < n; c++)
    x[c] = b[c] / A[c][c];
  return 0;
}

static double pl_det3(double a1, double a2, double a3, double b1, double b2, double b3, double c1,
                      double c2, double c3) {
  return a1 * (b2 * c3 - b3 * c2) - a2 * (b1 * c3 - b3 * c1) + a3 * (b1 * c2 - b2 * c1);
}

typedef struct {
  double u, v, w;
} pl_pt;

/* Знаки опорной четвёрки — из её АФФИННОЙ ЗАВИСИМОСТИ, а не перебором (А94):
 * `λ` через миноры `3×3`, величина `t = |Σ λ w| / Σ |λ|` ограничена по
 * построению, тогда как «наибольшее `t` по знакам» выбирало вырожденную систему
 * и давало наклоны `1e31`. */
static int pl_ref4(const pl_pt r[4], double *pa, double *pb, double *pc, double *pt) {
  double lam[4];
  lam[0] = +pl_det3(1, 1, 1, r[1].u, r[2].u, r[3].u, r[1].v, r[2].v, r[3].v);
  lam[1] = -pl_det3(1, 1, 1, r[0].u, r[2].u, r[3].u, r[0].v, r[2].v, r[3].v);
  lam[2] = +pl_det3(1, 1, 1, r[0].u, r[1].u, r[3].u, r[0].v, r[1].v, r[3].v);
  lam[3] = -pl_det3(1, 1, 1, r[0].u, r[1].u, r[2].u, r[0].v, r[1].v, r[2].v);
  double sa = 0.0, sw = 0.0;
  for (int k = 0; k < 4; k++) {
    sa += fabs(lam[k]);
    sw += lam[k] * r[k].w;
  }
  if (!(sa > 0.0)) return 1;
  double t = fabs(sw) / sa, sgn = (sw >= 0.0) ? 1.0 : -1.0;
  int id[4] = {0, 1, 2, 3};
  for (int i = 0; i < 4; i++)
    for (int j = i + 1; j < 4; j++)
      if (fabs(lam[id[j]]) > fabs(lam[id[i]])) {
        int tmp = id[i];
        id[i] = id[j];
        id[j] = tmp;
      }
  double A[4][4], b[4], x[4];
  memset(A, 0, sizeof A);
  for (int q = 0; q < 3; q++) {
    int k = id[q];
    double sk = (lam[k] >= 0.0) ? sgn : -sgn;
    A[q][0] = r[k].u;
    A[q][1] = r[k].v;
    A[q][2] = 1.0;
    b[q] = r[k].w - sk * t;
  }
  if (pl_solve(A, b, x, 3) != 0) return 1;
  *pa = x[0];
  *pb = x[1];
  *pc = x[2];
  *pt = t;
  return 0;
}

#define PL_EXCH 64

/* Минимаксная плоскость по набору точек в раме `(eu, ev, n)`: обмен по четвёркам.
 * Возвращает наклоны `(a, b)` и сдвиг `c` графовой формы `w = a·u + b·v + c`. */
static int pl_minimax(const pl_pt *p, int32_t n, double *pa, double *pb, double *pc) {
  if (n < 4) {
    *pa = 0.0;
    *pb = 0.0;
    if (n > 0) {
      double lo = p[0].w, hi = p[0].w;
      for (int32_t i = 1; i < n; i++) {
        if (p[i].w < lo) lo = p[i].w;
        if (p[i].w > hi) hi = p[i].w;
      }
      *pc = 0.5 * (lo + hi);
    } else
      *pc = 0.0;
    return 0;
  }
  pl_pt ref[4];
  {
    int32_t iw0 = 0, iw1 = 0, iu0 = 0, iu1 = 0;
    for (int32_t i = 1; i < n; i++) {
      if (p[i].w < p[iw0].w) iw0 = i;
      if (p[i].w > p[iw1].w) iw1 = i;
      if (p[i].u < p[iu0].u) iu0 = i;
      if (p[i].u > p[iu1].u) iu1 = i;
    }
    ref[0] = p[iw0];
    ref[1] = p[iw1];
    ref[2] = p[iu0];
    ref[3] = p[iu1];
  }
  double a = 0.0, b = 0.0, c = 0.0, t = 0.0;
  if (pl_ref4(ref, &a, &b, &c, &t) != 0) {
    a = 0.0;
    b = 0.0;
    c = 0.0;
    t = 0.0;
  }
  for (int it = 0; it < PL_EXCH; it++) {
    double worst = -1.0;
    pl_pt wp = p[0];
    for (int32_t i = 0; i < n; i++) {
      double d = fabs(p[i].w - a * p[i].u - b * p[i].v - c);
      if (d > worst) {
        worst = d;
        wp = p[i];
      }
    }
    if (!(worst > t)) break;
    double ba = a, bb = b, bc = c, bt = t;
    int best = -1;
    for (int k = 0; k < 4; k++) {
      pl_pt cand[4];
      memcpy(cand, ref, sizeof cand);
      cand[k] = wp;
      double a2, b2, c2, t2;
      if (pl_ref4(cand, &a2, &b2, &c2, &t2) != 0) continue;
      if (t2 > bt) {
        bt = t2;
        ba = a2;
        bb = b2;
        bc = c2;
        best = k;
      }
    }
    if (best < 0) break;
    ref[best] = wp;
    a = ba;
    b = bb;
    c = bc;
    t = bt;
  }
  *pa = a;
  *pb = b;
  *pc = c;
  return 0;
}

/* ПОРЯДОК ПАР — ПОЛНЫЙ, И ЭТО ТОТ ЖЕ УРОК, ЧТО `cmp_pair` В `pmerge.c`: при
 * сравнении по одному приоритету равные пары остаются на произвол `qsort`, и
 * выход перестаёт быть воспроизводимым. Поэтому доопределение по номеру.
 * Указатель файловый: `qsort` не передаёт контекста, а построение лестницы
 * однопоточное. */
static const double *pl_pri = NULL;

static int pl_cmp(const void *x, const void *y) {
  int32_t a = *(const int32_t *)x, b = *(const int32_t *)y;
  double va = pl_pri[a], vb = pl_pri[b];
  if (va > vb) return -1; /* больше — раньше */
  if (va < vb) return 1;
  return (a < b) ? -1 : ((a > b) ? 1 : 0);
}

/* --- хеш пар рёбер ----------------------------------------------------------- */

static uint64_t pl_mix(uint64_t k) {
  k ^= k >> 33;
  k *= 0xFF51AFD7ED558CCDULL;
  k ^= k >> 33;
  k *= 0xC4CEB9FE1A85EC53ULL;
  k ^= k >> 33;
  return k;
}

/* --- построение узла --------------------------------------------------------- */

/* Радиус объемлющей сферы узла вокруг точки `c` — по ОПОРНЫМ ТОЧКАМ членов, то есть
 * по тем же данным, что дают `dmax`. Оболочка уже посчитана огрубителем (О7), и
 * лишней работы здесь нет. */
static double pl_radius(const hz_polyset *ps, const int32_t *memb, int32_t nm, const double c[3]) {
  double r2 = 0.0;
  for (int32_t i = 0; i < nm; i++) {
    const hz_poly *P = &ps->p[memb[i]];
    const double *S = ps->sup + (size_t)P->s0 * 3;
    for (int32_t q = 0; q < P->nsup; q++) {
      double dx = S[(size_t)q * 3] - c[0], dy = S[(size_t)q * 3 + 1] - c[1],
             dz = S[(size_t)q * 3 + 2] - c[2];
      double d2 = dx * dx + dy * dy + dz * dz;
      if (d2 > r2) r2 = d2;
    }
  }
  return sqrt(r2);
}

/* Плоскость и `dmax` узла по ИСХОДНЫМ полигонам (А134). Ориентация — из `Σ A·n`
 * (она остаётся ОТДЕЛЬНОЙ величиной, А22), затем минимакс в её раме, затем
 * точный `dmax` по опорным точкам. */
static int pl_fit(const hz_polyset *ps, const int32_t *memb, int32_t nm, pl_pt *buf, int32_t cap,
                  hz_lodnode *nd) {
  double sA = 0.0, sn[3] = {0, 0, 0}, so[3] = {0, 0, 0};
  for (int32_t i = 0; i < nm; i++) {
    const hz_poly *P = &ps->p[memb[i]];
    sA += P->area;
    for (int c = 0; c < 3; c++) {
      sn[c] += P->area * P->n[c];
      so[c] += P->area * P->org[c];
    }
  }
  if (!(sA > 0.0)) return 1;
  double nl = sqrt(sn[0] * sn[0] + sn[1] * sn[1] + sn[2] * sn[2]);
  if (!(nl > 0.0)) return 1;
  double nn[3], org[3];
  for (int c = 0; c < 3; c++) {
    nn[c] = sn[c] / nl;
    org[c] = so[c] / sA;
  }
  /* Встречные грани — две поверхности, а не одна (А17/А22). */
  for (int32_t i = 0; i < nm; i++) {
    const hz_poly *P = &ps->p[memb[i]];
    if (!(P->n[0] * nn[0] + P->n[1] * nn[1] + P->n[2] * nn[2] > 0.0)) return 1;
  }
  double eu[3], ev[3];
  {
    int ax = 0;
    for (int a = 1; a < 3; a++)
      if (fabs(nn[a]) < fabs(nn[ax])) ax = a;
    double e0[3] = {0, 0, 0};
    e0[ax] = 1.0;
    double d = e0[0] * nn[0] + e0[1] * nn[1] + e0[2] * nn[2];
    for (int a = 0; a < 3; a++)
      eu[a] = e0[a] - d * nn[a];
    double en = sqrt(eu[0] * eu[0] + eu[1] * eu[1] + eu[2] * eu[2]);
    if (!(en > 0.0)) return 1;
    for (int a = 0; a < 3; a++)
      eu[a] /= en;
    ev[0] = nn[1] * eu[2] - nn[2] * eu[1];
    ev[1] = nn[2] * eu[0] - nn[0] * eu[2];
    ev[2] = nn[0] * eu[1] - nn[1] * eu[0];
  }
  int32_t np = 0;
  for (int32_t i = 0; i < nm; i++) {
    const hz_poly *P = &ps->p[memb[i]];
    const double *S = ps->sup + (size_t)P->s0 * 3;
    for (int32_t q = 0; q < P->nsup && np < cap; q++) {
      double dx = S[(size_t)q * 3] - org[0], dy = S[(size_t)q * 3 + 1] - org[1],
             dz = S[(size_t)q * 3 + 2] - org[2];
      buf[np].u = dx * eu[0] + dy * eu[1] + dz * eu[2];
      buf[np].v = dx * ev[0] + dy * ev[1] + dz * ev[2];
      buf[np].w = dx * nn[0] + dy * nn[1] + dz * nn[2];
      np++;
    }
  }
  if (np < 1) return 1;
  double a = 0.0, b = 0.0, c = 0.0;
  pl_minimax(buf, np, &a, &b, &c);
  double N[3];
  for (int k = 0; k < 3; k++)
    N[k] = nn[k] - a * eu[k] - b * ev[k];
  double Nl = sqrt(N[0] * N[0] + N[1] * N[1] + N[2] * N[2]);
  if (!(Nl > 0.0)) return 1;
  /* Наклон за пределы полусферы нормалей — откат к средневзвешенной (А95). */
  int okc = 1;
  for (int32_t i = 0; i < nm && okc; i++) {
    const hz_poly *P = &ps->p[memb[i]];
    if (!(P->n[0] * N[0] + P->n[1] * N[1] + P->n[2] * N[2] > 0.0)) okc = 0;
  }
  if (okc) {
    for (int k = 0; k < 3; k++)
      nd->n[k] = N[k] / Nl;
    nd->off = (N[0] * org[0] + N[1] * org[1] + N[2] * org[2] + c) / Nl;
  } else {
    for (int k = 0; k < 3; k++)
      nd->n[k] = nn[k];
    nd->off = nn[0] * org[0] + nn[1] * org[1] + nn[2] * org[2];
  }
  /* `dmax` — ТОЧНО, по тем же опорным точкам, но от ПРИНЯТОЙ плоскости. */
  double dm = 0.0;
  for (int32_t i = 0; i < nm; i++) {
    const hz_poly *P = &ps->p[memb[i]];
    const double *S = ps->sup + (size_t)P->s0 * 3;
    for (int32_t q = 0; q < P->nsup; q++) {
      double d = fabs(S[(size_t)q * 3] * nd->n[0] + S[(size_t)q * 3 + 1] * nd->n[1] +
                      S[(size_t)q * 3 + 2] * nd->n[2] - nd->off);
      if (d > dm) dm = d;
    }
  }
  nd->dmax = dm;
  nd->area_surf = sA;
  nd->cx = org[0];
  nd->cy = org[1];
  nd->cz = org[2];
  nd->rad = pl_radius(ps, memb, nm, org);
  nd->nmemb = nm;
  return 0;
}

int hz_lod_build(hz_lod *L, const hz_objmesh *m, const hz_pseglist *sg, const hz_polyset *ps,
                 double delta0, int maxlev, double eps) {
  memset(L, 0, sizeof *L);
  (void)m;
  const int32_t np = (ps->np < sg->nseg) ? ps->np : sg->nseg;
  if (np <= 0 || maxlev < 1) return 1;
  L->np = np;
  L->eps = eps;
  L->delta0 = delta0;
  L->lab = malloc((size_t)maxlev * (size_t)np * sizeof *L->lab);
  L->ndcap = 4 * np + 16;
  L->nd = malloc((size_t)L->ndcap * sizeof *L->nd);
  int32_t *memb = malloc((size_t)np * sizeof *memb);
  /* Буфер опорных точек: у самого крупного узла их может быть много, поэтому
   * берётся по всей сцене (иначе подгонка молча усечёт множество). */
  int32_t cap = (int32_t)(ps->nsupall > 0 ? ps->nsupall : 1);
  pl_pt *buf = malloc((size_t)cap * sizeof *buf);
  if (L->lab == NULL || L->nd == NULL || memb == NULL || buf == NULL) {
    free(memb);
    free(buf);
    hz_lod_free(L);
    return 2;
  }

  /* --- уровень 0: узел на полигон --- */
  for (int32_t k = 0; k < np; k++) {
    hz_lodnode *nd = &L->nd[k];
    memset(nd, 0, sizeof *nd);
    nd->level = 0;
    nd->parent = -1;
    int32_t one = k;
    if (pl_fit(ps, &one, 1, buf, cap, nd) != 0) {
      /* Полигон без опорных точек: плоскость своя, `dmax` — свой. */
      for (int c = 0; c < 3; c++)
        nd->n[c] = ps->p[k].n[c];
      nd->off = ps->p[k].off;
      nd->dmax = ps->p[k].dmax;
      nd->area_surf = ps->p[k].area;
      nd->rad = pl_radius(ps, &k, 1, ps->p[k].org);
      nd->cx = ps->p[k].org[0];
      nd->cy = ps->p[k].org[1];
      nd->cz = ps->p[k].org[2];
      nd->nmemb = 1;
    }
    nd->area_elem = ps->p[k].mom[0];
    L->lab[k] = k;
  }
  L->nnd = np;
  L->nlev = 1;

  /* --- рёбра уровня 0: пары сварных номеров и их владельцы --- */
  int32_t nbe = 0;
  for (int32_t k = 0; k < np; k++)
    for (int32_t l = ps->p[k].l0; l < ps->p[k].l0 + ps->p[k].nloop; l++)
      nbe += ps->loop[l + 1] - ps->loop[l];
  int64_t *ekey = malloc((size_t)(nbe > 0 ? nbe : 1) * sizeof *ekey);
  int32_t *eown = malloc((size_t)(nbe > 0 ? nbe : 1) * sizeof *eown);
  double *elen = malloc((size_t)(nbe > 0 ? nbe : 1) * sizeof *elen);
  if (ekey == NULL || eown == NULL || elen == NULL) {
    free(ekey);
    free(eown);
    free(elen);
    free(memb);
    free(buf);
    hz_lod_free(L);
    return 2;
  }
  int32_t ne = 0;
  for (int32_t k = 0; k < np; k++)
    for (int32_t l = ps->p[k].l0; l < ps->p[k].l0 + ps->p[k].nloop; l++) {
      int32_t b = ps->loop[l], e = ps->loop[l + 1], n = e - b;
      if (n < 2) continue;
      for (int32_t i = 0; i < n; i++) {
        int32_t w0 = ps->bw[b + i], w1 = ps->bw[b + (i + 1) % n];
        if (w0 < 0 || w1 < 0) continue;
        int32_t lo = (w0 < w1) ? w0 : w1, hi = (w0 < w1) ? w1 : w0;
        double x0[3], x1[3];
        hz_poly_world(&ps->p[k], ps->bv[(size_t)(b + i) * 2], ps->bv[(size_t)(b + i) * 2 + 1], x0);
        int32_t j = b + (i + 1) % n;
        hz_poly_world(&ps->p[k], ps->bv[(size_t)j * 2], ps->bv[(size_t)j * 2 + 1], x1);
        double dl = sqrt((x1[0] - x0[0]) * (x1[0] - x0[0]) + (x1[1] - x0[1]) * (x1[1] - x0[1]) +
                         (x1[2] - x0[2]) * (x1[2] - x0[2]));
        ekey[ne] = (int64_t)lo * 2147483647LL + hi;
        eown[ne] = k;
        elen[ne] = dl;
        ne++;
      }
    }

  /* --- уровни 1… ---
   * ВНУТРИ УРОВНЯ ИДУТ ПРОХОДЫ, и это исправление по замеру: с одним проходом
   * сокращение выходило `1.14…1.23×` вместо четырёх, потому что отвергнутая
   * воротами четвёрка больше не пробовалась, а её члены уходили на уровень
   * тождественно. Проходы дают им новых партнёров при том же допуске, и каждая
   * группа по-прежнему подгоняется ОДИН раз (§54.1). */
  int32_t *tmp = malloc((size_t)np * sizeof *tmp);
  if (tmp == NULL) {
    free(memb);
    free(buf);
    free(ekey);
    free(eown);
    free(elen);
    hz_lod_free(L);
    return 2;
  }
  for (int lev = 1; lev < maxlev; lev++) {
    int32_t *now = L->lab + (size_t)lev * (size_t)np;
    memcpy(now, L->lab + (size_t)(lev - 1) * (size_t)np, (size_t)np * sizeof *now);
    double dlev = delta0 * pow(2.0, (double)lev);
    int32_t total_made = 0;
    for (int pass = 0; pass < 8; pass++) {
      const int32_t *prev = now;
      int32_t *out = tmp;
      /* Смежность узлов предыдущего уровня по ОБЩИМ РЁБРАМ (А130): у ребра с одним
       * и тем же ключом два разных владельца-узла. Точечное касание сюда не
       * попадает вовсе, потому что ключ есть ПАРА вершин. */
      int32_t nprev = 0;
      for (int32_t k = 0; k < np; k++)
        if (prev[k] + 1 > nprev) nprev = prev[k] + 1;
      /* `grp` и `cur` индексируются ГЛОБАЛЬНЫМИ номерами узлов, а они со второго
       * уровня больше числа полигонов — поэтому выделяются здесь, по `nprev`, а не
       * снаружи по `np`. Первая редакция падала именно на этом. */
      int32_t *grp = malloc((size_t)(nprev > 0 ? nprev : 1) * sizeof *grp);
      int32_t *cur = malloc((size_t)(nprev > 0 ? nprev : 1) * sizeof *cur);
      if (grp == NULL || cur == NULL) {
        free(grp);
        free(cur);
        grp = NULL;
        cur = NULL;
        break;
      }
      /* Сортировкой по ключу найти пары: собираем массив (ключ, узел). */
      int64_t *kk = malloc((size_t)(ne > 0 ? ne : 1) * sizeof *kk);
      int32_t *kn = malloc((size_t)(ne > 0 ? ne : 1) * sizeof *kn);
      double *kl = malloc((size_t)(ne > 0 ? ne : 1) * sizeof *kl);
      if (kk == NULL || kn == NULL || kl == NULL) {
        free(kk);
        free(kn);
        free(kl);
        break;
      }
      for (int32_t i = 0; i < ne; i++) {
        kk[i] = ekey[i];
        kn[i] = prev[eown[i]];
        kl[i] = elen[i];
      }
      /* Хеш-таблица «ключ ребра → первый встреченный узел». */
      int64_t hsz = 4;
      while (hsz < 2 * (int64_t)ne + 8)
        hsz *= 2;
      int64_t *hk = calloc((size_t)hsz, sizeof *hk);
      int32_t *hv = malloc((size_t)hsz * sizeof *hv);
      /* Периметр узла: ребро, у которого владельцы РАЗНЫЕ узлы (или пары нет), —
       * граничное. Считается тут же, потому что нужен для формы `P/√A`. */
      double *per = calloc((size_t)(nprev > 0 ? nprev : 1), sizeof *per);
      if (hk == NULL || hv == NULL || per == NULL) {
        free(hk);
        free(hv);
        free(per);
        free(kk);
        free(kn);
        free(kl);
        break;
      }
      /* Пары смежности собираем в список. */
      int32_t pcap = ne + 8, pn = 0;
      int32_t *pa = malloc((size_t)pcap * sizeof *pa);
      int32_t *pb = malloc((size_t)pcap * sizeof *pb);
      if (pa == NULL || pb == NULL) {
        free(pa);
        free(pb);
        free(hk);
        free(hv);
        free(per);
        free(kk);
        free(kn);
        free(kl);
        break;
      }
      for (int32_t i = 0; i < ne; i++) {
        uint64_t h = pl_mix((uint64_t)kk[i] + 1) & (uint64_t)(hsz - 1);
        for (;;) {
          if (hk[h] == 0) {
            hk[h] = kk[i] + 1;
            hv[h] = kn[i];
            per[kn[i]] += kl[i];
            break;
          }
          if (hk[h] == kk[i] + 1) {
            int32_t o = hv[h];
            if (o != kn[i]) {
              if (pn < pcap) {
                pa[pn] = (o < kn[i]) ? o : kn[i];
                pb[pn] = (o < kn[i]) ? kn[i] : o;
                pn++;
              }
              /* Ребро внутреннее для пары, но граничное для каждого из двух узлов:
               * периметр уже учтён у первого, добавляем второму. */
              per[kn[i]] += kl[i];
            } else {
              /* Оба полурёбра одного узла — ребро ВНУТРЕННЕЕ, снимаем с периметра. */
              per[o] -= kl[i];
            }
            break;
          }
          h = (h + 1) & (uint64_t)(hsz - 1);
        }
      }
      L->npair_edge += pn;
      /* Периметр — в узлы текущего разбиения: он посчитан здесь и нужен и для
       * приоритета формы, и для доклада `P/√A`. Без этого поле оставалось нулём,
       * и печать формы показывала нули (найдено прогоном на зале). */
      for (int32_t k = 0; k < np; k++)
        L->nd[prev[k]].perim = per[prev[k]];
      /* Приоритет: сперва наибольший косинус угла, при равенстве — компактность
       * объединения по `P/√A` (§55). Сортировка простая: ключ = косинус, вниз. */
      double *pri = malloc((size_t)(pn > 0 ? pn : 1) * sizeof *pri);
      int32_t *ord = malloc((size_t)(pn > 0 ? pn : 1) * sizeof *ord);
      if (pri == NULL || ord == NULL) {
        free(pri);
        free(ord);
        free(pa);
        free(pb);
        free(hk);
        free(hv);
        free(per);
        free(kk);
        free(kn);
        free(kl);
        break;
      }
      for (int32_t i = 0; i < pn; i++) {
        const hz_lodnode *A = &L->nd[pa[i]], *B = &L->nd[pb[i]];
        double cs = A->n[0] * B->n[0] + A->n[1] * B->n[1] + A->n[2] * B->n[2];
        double as = A->area_surf + B->area_surf;
        double pp = per[pa[i]] + per[pb[i]];
        double sh = (as > 0.0 && pp > 0.0) ? pp / sqrt(as) : 1e9;
        /* Больше — раньше: косинус доминирует, форма доопределяет. */
        pri[i] = cs - 1e-3 * sh;
        ord[i] = i;
      }
      pl_pri = pri;
      qsort(ord, (size_t)pn, sizeof *ord, pl_cmp);

      /* Группировка: примерно по четыре, жадно по приоритету. */
      for (int32_t i = 0; i < nprev; i++)
        grp[i] = -1;
      int32_t *gsz = calloc((size_t)(nprev > 0 ? nprev : 1), sizeof *gsz);
      if (gsz == NULL) {
        free(pri);
        free(ord);
        free(pa);
        free(pb);
        free(hk);
        free(hv);
        free(per);
        free(kk);
        free(kn);
        free(kl);
        break;
      }
      int32_t ngrp = 0;
      for (int32_t q = 0; q < pn; q++) {
        int32_t i = ord[q], A = pa[i], B = pb[i];
        int32_t ga = grp[A], gb = grp[B];
        if (ga >= 0 && gb >= 0) continue;
        if (ga < 0 && gb < 0) {
          grp[A] = grp[B] = ngrp;
          gsz[ngrp] = 2;
          ngrp++;
          continue;
        }
        int32_t g = (ga >= 0) ? ga : gb, add = (ga >= 0) ? B : A;
        if (gsz[g] >= 4) continue; /* четвёрка — СТРУКТУРА (§54.2) */
        grp[add] = g;
        gsz[g]++;
      }
      for (int32_t i = 0; i < nprev; i++)
        if (grp[i] < 0) {
          grp[i] = ngrp;
          gsz[ngrp] = 1;
          ngrp++;
        }
      /* СПИСКИ ЧЛЕНОВ — ОДНИМ ПРОХОДОМ, А НЕ СКАНОМ НА ГРУППУ. Первая редакция
       * искала членов группы обходом всех полигонов, то есть стоила
       * `O(np · групп)`: на городе это `3e5 × 7.5e4 = 2.3e10` на проход, и
       * прогон бы не кончился. Здесь — счётная сортировка за `O(np)`. */
      for (int32_t i = 0; i < nprev; i++)
        cur[i] = -1;
      int32_t *goff = calloc((size_t)ngrp + 1, sizeof *goff);
      int32_t *gmem = malloc((size_t)np * sizeof *gmem);
      int32_t *gcur = malloc((size_t)(ngrp > 0 ? ngrp : 1) * sizeof *gcur);
      if (goff == NULL || gmem == NULL || gcur == NULL) {
        free(goff);
        free(gmem);
        free(gcur);
        free(grp);
        free(cur);
        free(gsz);
        free(pri);
        free(ord);
        free(pa);
        free(pb);
        free(hk);
        free(hv);
        free(per);
        free(kk);
        free(kn);
        free(kl);
        break;
      }
      for (int32_t k = 0; k < np; k++)
        goff[grp[prev[k]] + 1]++;
      for (int32_t g = 0; g < ngrp; g++)
        goff[g + 1] += goff[g];
      memcpy(gcur, goff, (size_t)ngrp * sizeof *gcur);
      for (int32_t k = 0; k < np; k++)
        gmem[gcur[grp[prev[k]]]++] = k;
      free(gcur);
      int32_t made = 0;
      for (int32_t g = 0; g < ngrp; g++) {
        int32_t nm = goff[g + 1] - goff[g];
        if (nm <= 0) continue;
        for (int32_t i = 0; i < nm; i++)
          memb[i] = gmem[goff[g] + i];
        /* НИЖНЯЯ ОЦЕНКА, БЕСПЛАТНАЯ И ТОЧНАЯ: `dmax` объединения не меньше `dmax`
         * любого его члена, потому что собственная плоскость члена для него
         * оптимальна. Значит член с `dmax ≥ δ` уровня делает группу безнадёжной,
         * и полную подгонку звать незачем. Отсев ОТВЕРГАЮЩИЙ, в отличие от
         * мажоранты `dev_box`, которая умеет только принимать (§54.1). */
        int hopeless = 0;
        for (int32_t i = 0; i < nm && !hopeless; i++)
          if (!(L->nd[prev[memb[i]]].dmax < dlev)) hopeless = 1;
        if (hopeless && nm > 1) {
          for (int32_t i = 0; i < nm; i++)
            out[memb[i]] = prev[memb[i]];
          L->nlow_rej++;
          continue;
        }
        if (L->nnd + 8 > L->ndcap) {
          int32_t nc = L->ndcap * 2;
          hz_lodnode *nn2 = realloc(L->nd, (size_t)nc * sizeof *nn2);
          if (nn2 == NULL) break;
          L->nd = nn2;
          L->ndcap = nc;
        }
        hz_lodnode cand;
        memset(&cand, 0, sizeof cand);
        cand.level = lev;
        cand.parent = -1;
        int ok = (pl_fit(ps, memb, nm, buf, cap, &cand) == 0);
        if (ok && !(cand.dmax < dlev)) {
          ok = 0;
          L->ngate_rej++;
        }
        if (ok) {
          int32_t id = L->nnd++;
          L->nd[id] = cand;
          for (int32_t k = 0; k < np; k++)
            if (grp[prev[k]] == g) out[k] = id;
          for (int32_t k = 0; k < np; k++)
            if (grp[prev[k]] == g) L->nd[prev[k]].parent = id;
          made++;
        } else {
          /* ТОЖДЕСТВЕННОЕ ПРОДВИЖЕНИЕ БЕЗ СОЗДАНИЯ УЗЛА: элемент остаётся собой.
           * Первая редакция дублировала узел на каждом уровне и проходе, раздувая
           * дерево и ничего не добавляя: вложенность «каждый узел предыдущего
           * уровня целиком в одном узле следующего» тождественным отображением
           * выполняется сама. */
          for (int32_t k = 0; k < np; k++)
            if (grp[prev[k]] == g) out[k] = prev[k];
          L->nident++;
        }
      }
      free(goff);
      free(gmem);
      free(grp);
      free(cur);
      free(gsz);
      free(pri);
      free(ord);
      free(pa);
      free(pb);
      free(hk);
      free(hv);
      free(per);
      free(kk);
      free(kn);
      free(kl);
      if (made > 0) memcpy(now, tmp, (size_t)np * sizeof *now);
      total_made += made;
      if (made == 0) break;
    }
    L->nlev = lev + 1;
    if (total_made == 0) break; /* слить больше нечего */
  }
  for (int32_t i = 0; i < L->nnd; i++) {
    const hz_lodnode *nd = &L->nd[i];
    if (nd->area_surf > 0.0 && nd->area_elem > 0.0) {
      double r = nd->area_elem / nd->area_surf;
      if (r > L->area_grow) L->area_grow = r;
    }
  }
  free(tmp);
  free(memb);
  free(buf);
  free(ekey);
  free(eown);
  free(elen);
  return 0;
}

void hz_lod_free(hz_lod *L) {
  free(L->nd);
  free(L->lab);
  memset(L, 0, sizeof *L);
}

int32_t hz_lod_cut(const hz_lod *L, const double eye[3], double eps, int nosin, int32_t *out) {
  /* Для каждого исходного полигона — САМЫЙ ГРУБЫЙ уровень, чей узел ещё
   * удовлетворяет `dmax < ε·R`. Расстояние берётся ЗДЕСЬ, в узле его нет (А133). */
  for (int32_t k = 0; k < L->np; k++) {
    int32_t sel = L->lab[k];
    for (int lev = 1; lev < L->nlev; lev++) {
      int32_t id = L->lab[(size_t)lev * (size_t)L->np + (size_t)k];
      const hz_lodnode *nd = &L->nd[id];
      double dx = nd->cx - eye[0], dy = nd->cy - eye[1], dz = nd->cz - eye[2];
      double Rc = sqrt(dx * dx + dy * dy + dz * dz);
      /* РАССТОЯНИЕ ДО БЛИЖАЙШЕЙ ТОЧКИ УЗЛА, А НЕ ДО ЦЕНТРА (§72). Критерий обязан
       * держаться ВЕЗДЕ на узле, а худшее место — ближайшее к глазу. */
      double R = Rc - nd->rad;
      if (!(R > 0.0)) break;
      /* Проецируемая ошибка: смещение вдоль луча экран почти не двигает (§27).
       * Косинус делится на `Rc` — длину ТОГО ЖЕ вектора `(dx,dy,dz)`, а не на
       * укороченное `R`. На `R` он зашкаливал за единицу, обрезался, `sin`
       * обращался в ноль, и критерий начинал пропускать что угодно: срез на зале
       * дал 861 элемент вместо 966 при СТРОГОМ правиле, чего быть не может. */
      double cs = (nd->n[0] * dx + nd->n[1] * dy + nd->n[2] * dz) / Rc;
      if (cs > 1.0) cs = 1.0;
      if (cs < -1.0) cs = -1.0;
      double sn = nosin ? 1.0 : sqrt(1.0 - cs * cs);
      if (nd->dmax * sn < eps * R)
        sel = id;
      else
        break;
    }
    out[k] = sel;
  }
  int32_t cnt = 0;
  {
    int32_t *seen = calloc((size_t)L->nnd, sizeof *seen);
    if (seen == NULL) return -1;
    for (int32_t k = 0; k < L->np; k++)
      if (!seen[out[k]]) {
        seen[out[k]] = 1;
        cnt++;
      }
    free(seen);
  }
  return cnt;
}

int hz_lod_seglist(const hz_lod *L, const hz_objmesh *m, const hz_pseglist *sg, const int32_t *cut,
                   hz_pseglist *so) {
  memset(so, 0, sizeof *so);
  int32_t *rank = malloc((size_t)L->nnd * sizeof *rank);
  so->label = malloc((size_t)m->nt * sizeof *so->label);
  so->seg = calloc((size_t)L->nnd, sizeof *so->seg);
  if (rank == NULL || so->label == NULL || so->seg == NULL) {
    free(rank);
    hz_seg_free(so);
    return 2;
  }
  for (int32_t i = 0; i < L->nnd; i++)
    rank[i] = -1;
  int32_t nn = 0;
  for (int32_t k = 0; k < L->np; k++) {
    int32_t id = cut[k];
    if (rank[id] < 0) {
      rank[id] = nn;
      so->seg[nn] = sg->seg[k];
      for (int c = 0; c < 3; c++)
        so->seg[nn].n[c] = L->nd[id].n[c];
      so->seg[nn].off = L->nd[id].off;
      so->seg[nn].area = 0.0;
      so->seg[nn].ntri = 0;
      so->seg[nn].dmax = 0.0;
      nn++;
    }
  }
  for (int32_t t = 0; t < m->nt; t++) {
    int32_t k = sg->label[t];
    int32_t g = (k < L->np) ? rank[cut[k]] : 0;
    so->label[t] = g;
    so->seg[g].ntri++;
    so->seg[g].area += hz_obj_tri_area(m, t);
    double p[3][3];
    hz_obj_tri(m, t, p);
    for (int i = 0; i < 3; i++) {
      double d = fabs(p[i][0] * so->seg[g].n[0] + p[i][1] * so->seg[g].n[1] +
                      p[i][2] * so->seg[g].n[2] - so->seg[g].off);
      if (d > so->seg[g].dmax) so->seg[g].dmax = d;
    }
  }
  so->nseg = nn;
  so->delta = sg->delta;
  free(rank);
  return 0;
}

/* --- ВТОРОЙ ПОСТРОИТЕЛЬ: УРОВЕНЬ = `hz_merge` НАД ПРЕДЫДУЩИМ -----------------
 * Разбор — в `plod.h`. Здесь только устройство: на каждом уровне слияние
 * работает над РАЗБИЕНИЕМ предыдущего уровня, поэтому вложенность выходит по
 * построению, а качество равно измеренному у плоского слияния — это тот же код.
 */
int hz_lod_build_merge(hz_lod *L, const hz_objmesh *m, const hz_pseglist *sg, const hz_polyset *ps0,
                       const hz_lodcfg *cf) {
  const double delta0 = cf->delta0, eps = cf->eps;
  const int maxlev = cf->maxlev;
  const double base = (cf->base > 0.0) ? cf->base : 2.0;
  memset(L, 0, sizeof *L);
  L->bands = cf->bands;
  L->bycount = cf->bycount;
  L->byangle = cf->byangle;
  L->radmul = (cf->radmul > 0.0) ? cf->radmul : 1.0;
  /* ПЕРВАЯ СТУПЕНЬ УГЛОВОЙ ЛЕСТНИЦЫ — ПАРАМЕТР, А НЕ КОНСТАНТА (А174). Стояло
   * `0.5°`, то есть НИЖЕ всего, что в сцене бывает: у города минимальный излом
   * между смежными участками ровно `90°`, пара под прямым углом даёт полураствор
   * конуса `45°`, и лестница `0.5…64°` доходила до него лишь на последнем уровне.
   * Отсюда «насыщение угловой тактики», записанное в А171 как её свойство, — а это
   * было свойство ЛЕСТНИЦЫ. Та же ошибка, что с `δ0`: начинать там, где у сцены
   * ничего нет. Ноль читается как прежние `0.5°`, чтобы старые прогоны совпали. */
  L->angle0 = (cf->angle0 > 0.0) ? cf->angle0 : 0.5;
  const int32_t np = (ps0->np < sg->nseg) ? ps0->np : sg->nseg;
  if (np <= 0 || maxlev < 1) return 1;
  L->np = np;
  L->eps = eps;
  L->delta0 = delta0;
  L->lab = malloc((size_t)maxlev * (size_t)np * sizeof *L->lab);
  L->ndcap = np + 16;
  L->nd = malloc((size_t)L->ndcap * sizeof *L->nd);
  if (L->lab == NULL || L->nd == NULL) {
    hz_lod_free(L);
    return 2;
  }
  /* Уровень 0 — узел на полигон; плоскость и `dmax` берутся у сегментации. */
  for (int32_t k = 0; k < np; k++) {
    hz_lodnode *nd = &L->nd[k];
    memset(nd, 0, sizeof *nd);
    nd->level = 0;
    nd->parent = -1;
    for (int c = 0; c < 3; c++)
      nd->n[c] = ps0->p[k].n[c];
    nd->off = ps0->p[k].off;
    nd->dmax = ps0->p[k].dmax;
    nd->area_surf = ps0->p[k].area;
    nd->area_elem = ps0->p[k].mom[0];
    nd->rad = pl_radius(ps0, &k, 1, ps0->p[k].org);
    nd->cx = ps0->p[k].org[0];
    nd->cy = ps0->p[k].org[1];
    nd->cz = ps0->p[k].org[2];
    nd->nmemb = 1;
    double per = 0.0;
    for (int32_t l = ps0->p[k].l0; l < ps0->p[k].l0 + ps0->p[k].nloop; l++) {
      int32_t b = ps0->loop[l], e = ps0->loop[l + 1], n = e - b;
      for (int32_t i = 0; i < n; i++) {
        const double *A = ps0->bv + (size_t)(b + i) * 2;
        const double *B = ps0->bv + (size_t)(b + (i + 1) % n) * 2;
        per += sqrt((B[0] - A[0]) * (B[0] - A[0]) + (B[1] - A[1]) * (B[1] - A[1]));
      }
    }
    nd->perim = per;
    L->lab[k] = k;
  }
  L->nnd = np;
  L->nlev = 1;

  /* Текущее разбиение: сегментация (уровень 0) и её полигоны. */
  /* ГАБАРИТ СЦЕНЫ — для тактики по числу: допуск, снятый с диагонали, не может
   * быть превышен ни одним отклонением внутри сцены, то есть ворота выключены
   * ДОКАЗУЕМО, а не «большим числом». */
  double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
  for (int32_t k = 0; k < np; k++)
    for (int c = 0; c < 3; c++) {
      double v = ps0->p[k].org[c];
      if (v < lo[c]) lo[c] = v;
      if (v > hi[c]) hi[c] = v;
    }
  double dscene = 0.0;
  for (int c = 0; c < 3; c++)
    dscene += (hi[c] - lo[c]) * (hi[c] - lo[c]);
  dscene = sqrt(dscene) + delta0;

  hz_pseglist cs = *sg;
  hz_polyset cp = *ps0;
  int owns = 0; /* владеем ли текущими `cs`/`cp` (уровни выше нулевого) */
  for (int lev = 1; lev < maxlev; lev++) {
    hz_mergecfg mc;
    memset(&mc, 0, sizeof mc);
    /* РАДИУС ПОИСКА ОТДЕЛЁН ОТ ДОПУСКА (§68): при `radmul > 1` кандидатами
     * становятся полигоны, отстоящие дальше допуска. Ошибка от этого не растёт —
     * её держат ВОРОТА, — а слияний становится больше: замер на зале дал
     * `202 -> 104` элемента при неизменном `dmax = 0.18 м`. */
    const double dlev = delta0 * pow(base, (double)lev);
    const double drad = dlev * L->radmul;
    if (L->byangle) {
      /* ТАКТИКА ПО УГЛУ. Ворота — полураствор конуса нормалей; ворота по `dmax`
       * сняты через `dgate` (А144), цель по числу снята. */
      mc.delta = drad;
      mc.dgate = dscene;
      /* Полураствор конуса ограничен `90°` по смыслу (полусфера), поэтому лестница
       * упирается в потолок, а не уходит за него: `89°` — последняя осмысленная
       * ступень, дальше условие вырождается в «конус существует». */
      mc.conemax = L->angle0 * pow(base, (double)(lev - 1));
      if (mc.conemax > 89.0) mc.conemax = 89.0;
      mc.target = 0;
    } else if (L->bycount) {
      /* ТАКТИКА ПО ЧИСЛУ. Цель — четверть предыдущего уровня: четвёрка здесь уже
       * не структура, а ТРЕБОВАНИЕ (§60). Допуск снят с ГАБАРИТА сцены, то есть
       * не ограничивает ничего: ни одно отклонение в сцене его не превысит. Это
       * не порог, подобранный под результат, а выключатель ворот — и достигнутый
       * `dmax` печатается уровнем как ЗАМЕР, а не проверяется. */
      /* РАДИУС кандидатов остаётся ДВОИЧНЫМ ДОПУСКОМ УРОВНЯ — иначе кандидатом
       * становится каждая пара, и город падает по памяти (А144). Выключаются
       * только ВОРОТА, через `dgate`. */
      mc.delta = drad;
      /* Ворота: либо сняты совсем (прежнее поведение), либо держатся на
       * `cgate · допуск уровня` — число ведёт, ошибка ограничена (§77). */
      mc.dgate = (cf->cgate > 0.0) ? cf->cgate * dlev : dscene;
      mc.target = (int32_t)((cs.nseg + 3) / 4);
      if (mc.target < 1) mc.target = 1;
    } else {
      mc.delta = drad;
      mc.dgate = (L->radmul > 1.0) ? dlev : 0.0;
      mc.target = 0; /* цель по числу — БЮДЖЕТ, а не критерий (§46) */
    }
    mc.use_geom = 1;
    mc.use_overlap = 1;
    mc.bands = L->bands;
    hz_pseglist so;
    hz_mergestat st;
    if (hz_merge(&so, m, &cs, &cp, &mc, &st) != 0) break;
    for (int q = 0; q < 12; q++)
      L->dih_hist[q] += st.dih_hist[q];
    L->ndih_conv += st.ndih_conv;
    L->ndih_conc += st.ndih_conc;
    L->ndih_flat += st.ndih_flat;
    if (so.nseg >= cs.nseg) {
      hz_seg_free(&so);
      break; /* слить больше нечего */
    }
    /* Разметка по ИСХОДНЫМ полигонам: у полигона все треугольники в одном узле. */
    int32_t *now = L->lab + (size_t)lev * (size_t)np;
    if (L->nnd + so.nseg + 1 > L->ndcap) {
      int32_t nc = L->nnd + so.nseg + 16;
      hz_lodnode *nn = realloc(L->nd, (size_t)nc * sizeof *nn);
      if (nn == NULL) {
        hz_seg_free(&so);
        break;
      }
      L->nd = nn;
      L->ndcap = nc;
    }
    int32_t nd0 = L->nnd; /* имя не `base`: оно занято основанием лестницы (§78) */
    for (int32_t g = 0; g < so.nseg; g++) {
      hz_lodnode *nd = &L->nd[nd0 + g];
      memset(nd, 0, sizeof *nd);
      nd->level = lev;
      nd->parent = -1;
      for (int c = 0; c < 3; c++)
        nd->n[c] = so.seg[g].n[c];
      nd->off = so.seg[g].off;
      nd->dmax = so.seg[g].dmax; /* по ВСЕМ треугольникам от плоскости (А134) */
      nd->area_surf = so.seg[g].area;
      nd->nmemb = 0;
    }
    L->nnd = nd0 + so.nseg;
    for (int32_t k = 0; k < np; k++) {
      const hz_poly *P = &ps0->p[k];
      int32_t t = (P->ntri > 0) ? ps0->tri[P->t0] : -1;
      int32_t g = (t >= 0 && t < m->nt) ? so.label[t] : 0;
      now[k] = nd0 + g;
      L->nd[nd0 + g].nmemb++;
      L->nd[L->lab[(size_t)(lev - 1) * (size_t)np + (size_t)k]].parent = nd0 + g;
    }
    /* Полигоны нового уровня: нужны и как вход следующего слияния, и ради
     * периметра с моментами (форма `P/√A`, §55). */
    hz_polyset npset;
    if (hz_poly_build(&npset, m, &so) != 0) {
      hz_seg_free(&so);
      break;
    }
    for (int32_t g = 0; g < so.nseg && g < npset.np; g++) {
      hz_lodnode *nd = &L->nd[nd0 + g];
      nd->area_elem = npset.p[g].mom[0];
      nd->rad = pl_radius(&npset, &g, 1, npset.p[g].org);
      nd->cx = npset.p[g].org[0];
      nd->cy = npset.p[g].org[1];
      nd->cz = npset.p[g].org[2];
      double per = 0.0;
      for (int32_t l = npset.p[g].l0; l < npset.p[g].l0 + npset.p[g].nloop; l++) {
        int32_t b = npset.loop[l], e = npset.loop[l + 1], n = e - b;
        for (int32_t i = 0; i < n; i++) {
          const double *A = npset.bv + (size_t)(b + i) * 2;
          const double *B = npset.bv + (size_t)(b + (i + 1) % n) * 2;
          per += sqrt((B[0] - A[0]) * (B[0] - A[0]) + (B[1] - A[1]) * (B[1] - A[1]));
        }
      }
      nd->perim = per;
    }
    if (owns) {
      hz_poly_free(&cp);
      hz_seg_free(&cs);
    }
    cs = so;
    cp = npset;
    owns = 1;
    L->nlev = lev + 1;
  }
  if (owns) {
    hz_poly_free(&cp);
    hz_seg_free(&cs);
  }
  for (int32_t i = 0; i < L->nnd; i++) {
    const hz_lodnode *nd = &L->nd[i];
    if (nd->area_surf > 0.0 && nd->area_elem > 0.0) {
      double r = nd->area_elem / nd->area_surf;
      if (r > L->area_grow) L->area_grow = r;
    }
  }
  return 0;
}
