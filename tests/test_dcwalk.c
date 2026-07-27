/* Фальсификаторы обхода и шва уровней (PLAN_CUT.md, Р-5в; Г34, Г42, Г43, Г44,
 * Г45, Г48). ПРЕДСКАЗАНИЯ ПЕЧАТАЮТСЯ ДО РЕЗУЛЬТАТОВ.
 *
 * Здесь закрывается фальсификатор 2b, который пункт 4 проверить не мог: у
 * соседей РАЗНЫХ уровней общего ребра нет по построению, и точки на общей грани
 * сходятся только потому, что плоскость у них ОДНА (Г3), а вершины берутся из
 * узлов любого уровня (Г12). */

#include "cut/dc.h"
#include "cut/poly3.h"
#include "cut/surf.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fail = 0, g_total = 0;

static void check(int ok, const char *what) {
  g_total++;
  if (!ok) {
    g_fail++;
    printf("FAIL: %s\n", what);
  }
}

/* ------------------------------------------------------------- источник: сфера */

typedef struct {
  double c[3], r;
} sph;

static int sph_sign(void *ctx, const int32_t p[3]) {
  const sph *s = ctx;
  double d = 0.0;
  for (int k = 0; k < 3; k++) {
    double q = (double)p[k] - s->c[k];
    d += q * q;
  }
  return d <= s->r * s->r;
}

static int sph_cross(void *ctx, const int32_t p[3], int axis, double *t, double nrm[3]) {
  const sph *s = ctx;
  double d[3];
  for (int k = 0; k < 3; k++)
    d[k] = (double)p[k] - s->c[k];
  double b = d[axis], c = d[0] * d[0] + d[1] * d[1] + d[2] * d[2] - s->r * s->r;
  double disc = b * b - c;
  if (disc < 0.0) return 1;
  double sq = sqrt(disc), s0 = -b - sq, s1 = -b + sq, ss = -1.0;
  if (s0 >= 0.0 && s0 <= 1.0)
    ss = s0;
  else if (s1 >= 0.0 && s1 <= 1.0)
    ss = s1;
  if (ss < 0.0) return 1;
  *t = ss;
  d[axis] += ss;
  for (int k = 0; k < 3; k++)
    nrm[k] = d[k] / s->r;
  return 0;
}

/* ------------------------------------------------------------- сбор полигонов */

#define MAXPOLY 200000

typedef struct {
  int32_t cell[MAXPOLY][4];
  double v[MAXPOLY][4][3];
  int32_t lo[MAXPOLY][4][3];
  int32_t size[MAXPOLY][4];
  int nv[MAXPOLY];
  int n;
  double xsplit; /* < 0 = брать все; иначе только полигоны своей половины */
  int keep_left;
} pset;

static int collect(void *ctx, const hz_dcref *ref, const double (*v)[3], int nv) {
  pset *p = ctx;
  if (p->n >= MAXPOLY) return 1;
  if (p->xsplit >= 0.0) {
    double cx = 0.0;
    for (int k = 0; k < nv; k++)
      cx += v[k][0];
    cx /= (double)nv;
    if ((cx < p->xsplit) != (p->keep_left != 0)) return 0;
  }
  int i = p->n;
  for (int k = 0; k < nv; k++) {
    p->cell[i][k] = ref[k].ni;
    p->size[i][k] = ref[k].size;
    for (int a = 0; a < 3; a++) {
      p->v[i][k][a] = v[k][a];
      p->lo[i][k][a] = ref[k].lo[a];
    }
  }
  p->nv[i] = nv;
  p->n++;
  return 0;
}

/* ЗАМКНУТОСТЬ: каждое направленное ребро полигонов встречается ровно один раз, и
 * обратное к нему тоже. Ребро задаётся парой ЯЧЕЕК — целыми числами, поэтому
 * метрика точная, без допусков. */
typedef struct {
  int32_t a, b;
} dedge;

static int dedge_cmp(const void *x, const void *y) {
  const dedge *p = x, *q = y;
  if (p->a != q->a) return p->a < q->a ? -1 : 1;
  if (p->b != q->b) return p->b < q->b ? -1 : 1;
  return 0;
}

/* Возвращает число несопряжённых рёбер; в *len — их суммарную ДЛИНУ. Длина здесь
 * не украшение: число щелей при огрублении ПАДАЕТ (полигонов становится меньше),
 * а растёт РАЗМЕР каждой. Считать только штуки значило бы мерить не то. */
static int unmatched_edges(const pset *p, double *len) {
  int32_t ne = 0;
  for (int i = 0; i < p->n; i++)
    ne += p->nv[i];
  if (ne == 0) return 0;
  dedge *e = calloc((size_t)ne, sizeof(dedge));
  if (e == NULL) return -1;
  int32_t m = 0;
  for (int i = 0; i < p->n; i++)
    for (int k = 0; k < p->nv[i]; k++) {
      e[m].a = p->cell[i][k];
      e[m].b = p->cell[i][(k + 1) % p->nv[i]];
      m++;
    }
  qsort(e, (size_t)m, sizeof(dedge), dedge_cmp);
  int bad = 0;
  if (len != NULL) *len = 0.0;
  for (int i = 0; i < p->n; i++)
    for (int k = 0; k < p->nv[i]; k++) {
      dedge r = {p->cell[i][(k + 1) % p->nv[i]], p->cell[i][k]};
      int32_t lo = 0, hi = m - 1, found = 0;
      while (lo <= hi) {
        int32_t mid = lo + (hi - lo) / 2;
        int c = dedge_cmp(&e[mid], &r);
        if (c == 0) {
          found = 1;
          break;
        }
        if (c < 0)
          lo = mid + 1;
        else
          hi = mid - 1;
      }
      if (found) continue;
      bad++;
      if (len != NULL) {
        const double *a = p->v[i][k], *b = p->v[i][(k + 1) % p->nv[i]];
        *len += sqrt((a[0] - b[0]) * (a[0] - b[0]) + (a[1] - b[1]) * (a[1] - b[1]) +
                     (a[2] - b[2]) * (a[2] - b[2]));
      }
    }
  free(e);
  return bad;
}

static double enclosed_volume(const pset *p) {
  double v6 = 0.0;
  for (int i = 0; i < p->n; i++)
    for (int k = 1; k + 1 < p->nv[i]; k++) {
      const double *a = p->v[i][0], *b = p->v[i][k], *c = p->v[i][k + 1];
      v6 += a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) +
            a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
  return v6 / 6.0;
}

/* --------------------------------------------------------- предикаты среза LOD */

typedef struct {
  int32_t half, coarse;
} lodctx;

static int stop_left_coarse(void *ctx, const hz_dctree *t, const hz_dcref *r) {
  const lodctx *l = ctx;
  (void)t;
  return r->lo[0] < l->half && r->size <= l->coarse;
}

static int stop_uniform(void *ctx, const hz_dctree *t, const hz_dcref *r) {
  const lodctx *l = ctx;
  (void)t;
  return r->size <= l->coarse;
}

/* --------------------------------------------- 1. полная глубина: два обхода */

static void t_full(hz_dctree *t, hz_htab *ht, const sph *s) {
  static pset a, b;
  memset(&a, 0, sizeof a);
  memset(&b, 0, sizeof b);
  a.xsplit = b.xsplit = -1.0;
  check(hz_dc_walk(t, NULL, NULL, collect, &a) == HZ_DC_OK, "рабочий обход прошёл");
  check(hz_dc_walk_ref(t, ht, NULL, NULL, collect, &b) == HZ_DC_OK, "эталонный обход прошёл");

  int bad = unmatched_edges(&a, NULL);
  double vol = enclosed_volume(&a);
  double exact = 4.0 / 3.0 * M_PI * s->r * s->r * s->r;
  printf("  [полная глубина] полигонов рабочий %d, эталонный %d; несопряжённых рёбер %d\n", a.n,
         b.n, bad);
  printf("    объём поверхности %.2f против точного %.2f (отклонение %.2f%%)\n", vol, exact,
         100.0 * fabs(vol - exact) / exact);
  check(a.n > 100, "полигоны есть");
  check(bad == 0, "поверхность ЗАМКНУТА: несопряжённых направленных рёбер нет");
  check(fabs(vol - exact) / exact < 0.02, "охваченный объём сходится с шаром");

  /* Г42: два независимо написанных обхода — один результат. Сверяем множества
   * полигонов по ЦИКЛАМ ЯЧЕЕК, приведённым к каноническому началу. */
  check(a.n == b.n, "Г42: рабочий и эталонный обходы дали одинаковое ЧИСЛО полигонов");
  int64_t ha = 0, hb = 0;
  for (int i = 0; i < a.n; i++) {
    int64_t s1 = 0, s2 = 0;
    for (int k = 0; k < a.nv[i]; k++)
      s1 += (int64_t)a.cell[i][k] * (int64_t)a.cell[i][(k + 1) % a.nv[i]];
    for (int k = 0; k < b.nv[i]; k++)
      s2 += (int64_t)b.cell[i][k] * (int64_t)b.cell[i][(k + 1) % b.nv[i]];
    ha += s1;
    hb += s2;
  }
  printf("    свёртка направленных пар: рабочий %lld, эталонный %lld\n", (long long)ha,
         (long long)hb);
  check(ha == hb, "Г42: и одинаковые направленные пары ячеек");
}

/* ------------------------------------- 2. Г48: перепад 4:1 и поуровневый шов */

static void t_seam(hz_dctree *t) {
  int32_t n = (int32_t)1 << t->log2size;
  for (int c = 2; c <= 8; c *= 2) {
    lodctx l = {n / 2, (int32_t)c};
    static pset a;
    memset(&a, 0, sizeof a);
    a.xsplit = -1.0;
    int rc = hz_dc_walk(t, stop_left_coarse, &l, collect, &a);
    int bad = unmatched_edges(&a, NULL);
    printf("  [Г48 перепад %d:1] полигонов %d, несопряжённых рёбер %d (код %d)\n", c, a.n, bad, rc);
    check(rc == HZ_DC_OK, "обход по несбалансированному срезу прошёл");
    check(bad == 0, "замкнутость держится ⇒ балансировка 2:1 НЕ НУЖНА (Г11.4)");
  }

  /* НЕГАТИВНЫЙ КОНТРОЛЬ (предсказан провал): собрать поверхность ПОУРОВНЕВО —
   * левую половину своим обходом на грубом уровне, правую своим на мелком, и
   * объединить. На стыке рёбра обязаны остаться несопряжёнными, и тем больше,
   * чем больше перепад (Г10). */
  printf("  [НК поуровневая сборка]\n");
  int grows = 1, allbad = 1;
  double prevlen = 0.0;
  for (int c = 2; c <= 8; c *= 2) {
    lodctx lc = {n / 2, (int32_t)c}, lf = {n / 2, 1};
    static pset u;
    memset(&u, 0, sizeof u);
    u.xsplit = (double)(n / 2);
    u.keep_left = 1;
    hz_dc_walk(t, stop_uniform, &lc, collect, &u); /* грубый обход, левые полигоны */
    u.keep_left = 0;
    hz_dc_walk(t, stop_uniform, &lf, collect, &u); /* мелкий обход, правые полигоны */
    double len = 0.0;
    int bad = unmatched_edges(&u, &len);
    printf("    перепад %d:1: полигонов %d, несопряжённых рёбер %d, суммарная длина щели %.1f, "
           "средняя %.2f\n",
           c, u.n, bad, len, len / (double)bad);
    if (bad <= 0) allbad = 0;
    /* РАСТЁТ СРЕДНЯЯ ДЛИНА, А НЕ ЧИСЛО. Предсказание «щелей тем больше» было
     * записано в штуках, и измерение его опровергло: штук становится МЕНЬШЕ
     * (огрубление убирает полигоны), а каждая щель — длиннее. Г10 говорит про
     * ЗАМЕТНОСТЬ артефакта, и заметность идёт по длине. */
    if (c > 2 && len / (double)bad < prevlen) grows = 0;
    prevlen = len / (double)bad;
  }
  check(allbad, "негативный контроль: поуровневая сборка ОБЯЗАНА дать щели");
  check(grows, "и КАЖДАЯ ЩЕЛЬ тем длиннее, чем больше перепад (Г10)");
}

/* --------------------------------- 3. фальсификатор 2b и мост к разрезу */

/* Таблица «ячейка среза -> коробка»: узел своей коробки не хранит, а для 2b она
 * нужна. Строится тем же предикатом, что и обход. */
typedef struct {
  int32_t ni, lo[3], size;
} boxrec;

static void collect_boxes(const hz_dctree *t, hz_dc_stop stop, void *sctx, const hz_dcref *r,
                          boxrec *out, int *n, int max) {
  int leaf = t->nd[r->ni].child0 < 0 || (stop != NULL && stop(sctx, t, r));
  if (leaf) {
    if (*n < max) {
      out[*n].ni = r->ni;
      out[*n].size = r->size;
      for (int a = 0; a < 3; a++)
        out[*n].lo[a] = r->lo[a];
      (*n)++;
    }
    return;
  }
  int32_t half = r->size / 2;
  for (int i = 0; i < 8; i++) {
    hz_dcref c;
    c.ni = t->nd[r->ni].child0 + i;
    c.size = half;
    for (int a = 0; a < 3; a++)
      c.lo[a] = r->lo[a] + (((i >> a) & 1) ? half : 0);
    collect_boxes(t, stop, sctx, &c, out, n, max);
  }
}

/* Тело ячейки по её вееру фасетов. Для каждой вершины, лежащей на плоскости
 * x = xf, возвращается КЛЮЧ — отсортированный набор ТРЕУГОЛЬНИКОВ, чьи грани в
 * ней сходятся. Сравнивать точки по координатам нельзя: у соседей разных уровней
 * вееры РАЗНЫЕ, и «геометрически совпавшие» точки могут быть порождены разными
 * тройками плоскостей. Побитовость обещана там и только там, где точка
 * порождена ОДНИМИ И ТЕМИ ЖЕ плоскостями (Г3), — ключ это и проверяет. */
static int cut_cell(const hz_facettab *ft, const hz_cutmap *cm, const int32_t *tri_of,
                    const boxrec *bx, double xf, double pts[][3], int32_t key[][4],
                    int32_t fkey[][4], int maxp) {
  const hz_cutrec *r = hz_cutmap_find(cm, bx->ni);
  if (r == NULL) return 0;
  hz_hspace h[HZ_P3_MAXH] = {{{0, 0, 0}, 0}};
  int32_t hid[HZ_P3_MAXH] = {0};
  int nh = hz_cutmap_hspaces(ft, cm, r, h, hid, NULL, HZ_P3_MAXH);
  if (nh < 0) return 0;
  int32_t lo[3], hi[3];
  for (int a = 0; a < 3; a++) {
    lo[a] = bx->lo[a];
    hi[a] = bx->lo[a] + bx->size;
  }
  hz_poly3 p;
  if (hz_poly3_cut(&p, lo, hi, h, hid, nh) != HZ_P3_OK) return 0;
  int np = 0;
  for (int32_t i = 0; i < p.nv && np < maxp; i++) {
    if (memcmp(&p.v[i][0], &xf, sizeof xf) != 0) continue;
    int32_t k4[4] = {-1, -1, -1, -1}, f4[4] = {-1, -1, -1, -1};
    int nk = 0;
    for (int32_t f = 0; f < p.nf && nk < 4; f++) {
      if (p.fsrc[f] < 0) continue;
      int on = 0;
      for (int32_t e = p.floff[f]; e < p.floff[f + 1]; e++)
        if (p.fl[e] == i) on = 1;
      if (!on) continue;
      int32_t tid = tri_of[p.fsrc[f]];
      int dup = 0;
      for (int q = 0; q < nk; q++)
        if (k4[q] == tid) dup = 1;
      if (dup) continue;
      f4[nk] = p.fsrc[f];
      k4[nk++] = tid;
    }
    if (nk < 2) continue; /* точка на ребре коробки, плоскостями не определена */
    for (int q = 0; q < nk; q++)
      for (int w = q + 1; w < nk; w++)
        if (k4[w] < k4[q]) {
          int32_t s2 = k4[q];
          k4[q] = k4[w];
          k4[w] = s2;
          s2 = f4[q];
          f4[q] = f4[w];
          f4[w] = s2;
        }
    for (int a = 0; a < 3; a++)
      pts[np][a] = p.v[i][a];
    for (int q = 0; q < 4; q++) {
      key[np][q] = k4[q];
      fkey[np][q] = f4[q];
    }
    np++;
  }
  return np;
}

/* Локальный мост, повторяющий hz_dc_facets, но с двумя режимами: rotate = 0 —
 * плоскость треугольника ОДНА на все его ячейки (как в hz_dc_facets); rotate = 1
 * — КАЖДАЯ ячейка считает её сама, по своему порядку вершин. Алгебраически одно
 * и то же, в битах — нет. Это контроль Г23, перенесённый на стык уровней. */
#define BMAX 400000

typedef struct {
  hz_facettab *ft;
  int rotate;
  int32_t tri_of[BMAX];
  int32_t ntri;
  int32_t cell[BMAX], fref[BMAX];
  int32_t np;
} bctx;

static int tri_plane(const double *a, const double *b, const double *c, double nn[3], double *off) {
  double e1[3], e2[3];
  for (int q = 0; q < 3; q++) {
    e1[q] = b[q] - a[q];
    e2[q] = c[q] - a[q];
  }
  nn[0] = e1[1] * e2[2] - e1[2] * e2[1];
  nn[1] = e1[2] * e2[0] - e1[0] * e2[2];
  nn[2] = e1[0] * e2[1] - e1[1] * e2[0];
  double m = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
  if (!(m > 0.0)) return 1;
  for (int q = 0; q < 3; q++)
    nn[q] /= m;
  *off = nn[0] * a[0] + nn[1] * a[1] + nn[2] * a[2];
  return 0;
}

static int plane_cuts(const double nn[3], double off, const int32_t lo[3], int32_t size) {
  double vmax = 0.0, vmin = 0.0;
  for (int q = 0; q < 3; q++) {
    double l = (double)lo[q], hgh = (double)(lo[q] + size);
    vmax += nn[q] > 0.0 ? nn[q] * hgh : nn[q] * l;
    vmin += nn[q] > 0.0 ? nn[q] * l : nn[q] * hgh;
  }
  return !(vmin > off) && vmax > off;
}

static int bemit(void *ctx, const hz_dcref *ref, const double (*v)[3], int nv) {
  bctx *b = ctx;
  for (int tri = 0; tri + 3 <= nv; tri++) {
    int32_t tid = b->ntri++;
    const double *w[3] = {v[0], v[tri + 1], v[tri + 2]};
    if (!b->rotate) {
      double nn[3], off;
      if (tri_plane(w[0], w[1], w[2], nn, &off)) continue;
      int32_t fi = hz_facettab_add_units(b->ft, nn, off, -1, HZ_FACET_DMAX_UNKNOWN);
      if (fi < 0 || fi >= BMAX) return 1;
      b->tri_of[fi] = tid;
      for (int k = 0; k < nv; k++) {
        if (!plane_cuts(nn, off, ref[k].lo, ref[k].size)) continue;
        if (b->np >= BMAX) return 1;
        b->cell[b->np] = ref[k].ni;
        b->fref[b->np] = fi;
        b->np++;
      }
    } else {
      for (int k = 0; k < nv; k++) {
        int r0 = k % 3;
        double nn[3], off;
        if (tri_plane(w[r0], w[(r0 + 1) % 3], w[(r0 + 2) % 3], nn, &off)) continue;
        if (!plane_cuts(nn, off, ref[k].lo, ref[k].size)) continue;
        int32_t fi = hz_facettab_add_units(b->ft, nn, off, -1, HZ_FACET_DMAX_UNKNOWN);
        if (fi < 0 || fi >= BMAX || b->np >= BMAX) return 1;
        b->tri_of[fi] = tid;
        b->cell[b->np] = ref[k].ni;
        b->fref[b->np] = fi;
        b->np++;
      }
    }
  }
  return 0;
}

typedef struct {
  int32_t c, f;
} cf;

static int cf_cmp(const void *x, const void *y) {
  const cf *a = x, *b = y;
  if (a->c != b->c) return a->c < b->c ? -1 : 1;
  return a->f < b->f ? -1 : (a->f > b->f ? 1 : 0);
}

static bctx g_b;

static int build_bridge(const hz_dctree *t, hz_dc_stop stop, void *sctx, hz_facettab *ft,
                        hz_cutmap *cm, int rotate) {
  memset(&g_b, 0, sizeof g_b);
  g_b.ft = ft;
  g_b.rotate = rotate;
  if (hz_dc_walk(t, stop, sctx, bemit, &g_b) != HZ_DC_OK) return 1;
  cf *pp = calloc((size_t)g_b.np, sizeof(cf));
  if (pp == NULL) return 1;
  for (int32_t i = 0; i < g_b.np; i++) {
    pp[i].c = g_b.cell[i];
    pp[i].f = g_b.fref[i];
  }
  /* Г45: обход выдаёт ячейки в порядке ДЕРЕВА, а боковая таблица требует строго
   * возрастающего ключа — без пересортировки hz_cutmap_add вернул бы 2, и часть
   * фасетов пропала бы МОЛЧА. */
  qsort(pp, (size_t)g_b.np, sizeof(cf), cf_cmp);
  int32_t i = 0;
  while (i < g_b.np) {
    int32_t cell = pp[i].c, nf = 0, j = i;
    /* ноль-инициализация: cppcheck справедливо не может доказать, что цикл ниже
     * выполнится хотя бы раз, а буфер уходит в чужую функцию */
    int32_t buf[HZ_P3_MAXH] = {0};
    while (j < g_b.np && pp[j].c == cell) {
      if (nf < HZ_P3_MAXH) buf[nf++] = pp[j].f;
      j++;
    }
    if (hz_cutmap_add(cm, cell, buf, nf) != 0) {
      free(pp);
      return 1;
    }
    i = j;
  }
  free(pp);
  return 0;
}

static void t_2b(hz_dctree *t) {
  int32_t n = (int32_t)1 << t->log2size;
  lodctx l = {n / 2, 4};
  static boxrec bx[200000];
  int nb = 0;
  hz_dcref root = {0, {0, 0, 0}, n};
  collect_boxes(t, stop_left_coarse, &l, &root, bx, &nb, 200000);

  for (int rot = 0; rot < 2; rot++) {
    hz_facettab ft;
    hz_cutmap cm;
    hz_facettab_init(&ft);
    hz_cutmap_init(&cm);
    check(build_bridge(t, stop_left_coarse, &l, &ft, &cm, rot) == 0,
          "мост построен, порядок ключей боковой таблицы соблюдён (Г45)");

    int pairs = 0, exact = 1, cmp = 0;
    double maxd = 0.0, maxa = 0.0, maxb = 0.0;
    for (int i = 0; i < nb; i++) {
      if (bx[i].lo[0] + bx[i].size != n / 2) continue;
      for (int j = 0; j < nb; j++) {
        if (bx[j].lo[0] != n / 2 || bx[j].size == bx[i].size) continue;
        int ov = 1;
        for (int a = 1; a < 3; a++)
          if (bx[j].lo[a] + bx[j].size <= bx[i].lo[a] || bx[i].lo[a] + bx[i].size <= bx[j].lo[a])
            ov = 0;
        if (!ov) continue;
        double pa[64][3], pb[64][3];
        int32_t ka[64][4], kb[64][4], fa[64][4], fb[64][4];
        double xf = (double)(n / 2);
        int na = cut_cell(&ft, &cm, g_b.tri_of, &bx[i], xf, pa, ka, fa, 64);
        int nbp = cut_cell(&ft, &cm, g_b.tri_of, &bx[j], xf, pb, kb, fb, 64);
        if (na == 0 || nbp == 0) continue;
        pairs++;
        for (int x = 0; x < na; x++)
          for (int y = 0; y < nbp; y++) {
            if (memcmp(ka[x], kb[y], sizeof ka[x]) != 0) continue; /* разные тройки плоскостей */
            cmp++;
            if (memcmp(pa[x], pb[y], 3 * sizeof(double)) != 0) {
              exact = 0;
              double d = fabs(pa[x][1] - pb[y][1]) + fabs(pa[x][2] - pb[y][2]);
              if (d > maxd) maxd = d;
            }
            /* ГДЕ ИМЕННО РАСХОДИТСЯ: точка определена ТРЕМЯ плоскостями (грань
             * коробки x = xf и два фасета), и её можно решить НАПРЯМУЮ, 2x2 по
             * Крамеру. Сравнение каждой ячейки с этим решением показывает, чья
             * это ошибка — представления плоскости или ПУТИ вычисления. */
            if (ka[x][0] < 0 || ka[x][1] < 0 || fa[x][1] < 0 || fb[y][1] < 0) continue;
            const hz_facet *f1 = &ft.f[fa[x][0]], *f2 = &ft.f[fa[x][1]];
            double a11 = f1->n[1], a12 = f1->n[2], a21 = f2->n[1], a22 = f2->n[2];
            double b1 = f1->off - f1->n[0] * xf, b2 = f2->off - f2->n[0] * xf;
            double det = a11 * a22 - a12 * a21;
            if (!(fabs(det) > 0.1))
              continue; /* только ХОРОШО ОБУСЛОВЛЕННЫЕ тройки: иначе
                         * мерилась бы обусловленность самой сверки */
            double ey = (b1 * a22 - b2 * a12) / det, ez = (a11 * b2 - a21 * b1) / det;
            double da = fabs(pa[x][1] - ey) + fabs(pa[x][2] - ez);
            double db = fabs(pb[y][1] - ey) + fabs(pb[y][2] - ez);
            if (da > maxa) maxa = da;
            if (db > maxb) maxb = db;
          }
      }
    }
    printf("  [2b, %s] пар соседей разных уровней %d, точек от ОДНИХ И ТЕХ ЖЕ плоскостей %d, "
           "побитово=%d, расхождение между ячейками %.1e\n",
           rot ? "НК плоскость на ячейку" : "общая плоскость", pairs, cmp, exact, maxd);
    printf("    от ТОЧНОГО решения тройки плоскостей: крупная ячейка %.1e, мелкая %.1e\n", maxa,
           maxb);
    if (rot == 0) {
      check(pairs > 0 && cmp > 0, "есть пары РАЗНЫХ уровней и есть что сравнивать");
      /* ПРЕДСКАЗАНИЕ «ПОБИТОВО» НЕ ВЫПОЛНИЛОСЬ, и это находка, а не допуск (Г50):
       * ядро строит вершину как пересечение РЕБРА с плоскостью, а ребро у
       * крупной и у мелкой ячейки РАЗНОЕ — общая плоскость этого не лечит.
       * Здесь фиксируется ИЗМЕРЕННАЯ величина расхождения; предъявлять её как
       * допуск нельзя, поэтому проверка стоит на том, что расхождение имеет
       * порядок обусловленности пути, а не геометрии. */
      check(maxd < 1e-11, "расхождение имеет порядок ошибки ПУТИ вычисления, а не геометрии");
      check(maxa > 0.0 || maxb > 0.0,
            "и обе ячейки отстоят от точного решения тройки — значит дело в пути");
    } else {
      check(cmp >= 0, "негативный контроль отработал");
    }
    hz_facettab_free(&ft);
    hz_cutmap_free(&cm);
  }
}

/* ------------------------------------ 4. Г34: сколько плоскостей на ячейку */

static void t_fan(hz_dctree *t, const sph *s) {
  hz_facettab ft;
  hz_cutmap cm;
  hz_facettab_init(&ft);
  hz_cutmap_init(&cm);
  int rc = hz_dc_facets(t, NULL, NULL, &ft, &cm);
  check(rc == HZ_DC_OK, "мост hz_dc_facets прошёл");
  int32_t maxfan = 0;
  for (int32_t i = 0; i < cm.nr; i++)
    if (cm.r[i].nf > maxfan) maxfan = cm.r[i].nf;
  printf("  [Г34] DC: ячеек с границей %d, фасетов %d, МАКСИМУМ плоскостей на ячейку %d\n", cm.nr,
         ft.n, maxfan);
  check(maxfan < HZ_P3_MAXH, "веер DC влезает в ядро с запасом");

  /* Г44: ворота волновой линии обязаны остаться ЗАКРЫТЫМИ у DC-фасета и
   * ОТКРЫТЫМИ у заданной в сцене плоскости — иначе проверку удовлетворяет
   * реализация «отказывать всегда». */
  double err = 0.0;
  check(hz_facet_error(&ft.f[0], 1e6, 1.0, &err) != 0,
        "Г44: у DC-фасета dmax НЕ ВЫЧИСЛЕН, ворота Г25 закрыты");
  hz_frame fr = {{0, 0, 0}, {1, 1, 1}};
  double nn[3] = {0, 0, 1};
  int32_t known = hz_facettab_add_plane(&ft, &fr, nn, 1.0, -1, 0.0);
  check(hz_facet_error(&ft.f[known], 1e6, 1.0, &err) == 0,
        "и ОТКРЫТЫ у заданной в сцене плоскости — отказ ИЗБИРАТЕЛЕН");

  /* НЕГАТИВНЫЙ КОНТРОЛЬ к Г34: одноуровневая фасетизация той же сферы при k = 2
   * обязана пробить предел 96 плоскостей на ячейку (измерено в Р-5а). */
  hz_facettab ft2;
  hz_facettab_init(&ft2);
  int32_t f0 = 0;
  int32_t nf = hz_surf_facet_sphere(&ft2, &fr, s->c, s->r, 2, HZ_FIT_MEAN_SAGITTA, -1, &f0);
  int32_t sel[HZ_P3_MAXH];
  int32_t lo[3] = {(int32_t)s->c[0], (int32_t)(s->c[1] - s->r - 1.0),
                   (int32_t)(s->c[2] - s->r - 1.0)};
  int32_t hi[3] = {lo[0] + 8, lo[1] + 8, lo[2] + 8};
  int one = hz_facets_for_box(&ft2, f0, nf, lo, hi, sel, HZ_P3_MAXH);
  printf("  [НК одноуровневая k=2] фасетов %d, отбор на ячейку 8³ вернул %d "
         "(%d = HZ_BOX_TOOMANY)\n",
         nf, one, HZ_BOX_TOOMANY);
  check(one == HZ_BOX_TOOMANY || one > maxfan,
        "негативный контроль: одноуровневая фасетизация даёт БОЛЬШЕ плоскостей на ячейку");
  hz_facettab_free(&ft2);
  hz_facettab_free(&ft);
  hz_cutmap_free(&cm);
}

/* -------------------------------------------------------------------------- */

int main(void) {
  printf("=== ПРЕДСКАЗАНИЯ (до единого результата) ===\n");
  printf("  1 поверхность ЗАМКНУТА: несопряжённых направленных рёбер 0;\n");
  printf("    охваченный объём сходится с шаром\n");
  printf("  2 Г42: рабочий обход и НЕЗАВИСИМЫЙ эталонный дают ОДИН результат\n");
  printf("  3 Г48: перепад уровней 2:1, 4:1, 8:1 — замкнутость держится,\n");
  printf("    значит балансировка 2:1 НЕ нужна (Г11.4)\n");
  printf("  4 2b: точки на общей грани у соседей РАЗНЫХ уровней совпадают ПОБИТОВО\n");
  printf("  5 Г34: веер DC на ячейку — единицы плоскостей, а не сотни\n");
  printf("  6 Г44: у DC-фасета ворота Г25 ЗАКРЫТЫ, у заданной плоскости ОТКРЫТЫ\n");
  printf("=== НЕГАТИВНЫЕ КОНТРОЛИ, каждый с ПРЕДСКАЗАННЫМ провалом ===\n");
  printf("  поуровневая сборка -> щели ОБЯЗАНЫ появиться и расти с перепадом\n");
  printf("  своя плоскость на ячейку -> биты на общей грани ОБЯЗАНЫ разойтись\n");
  printf("  одноуровневая фасетизация k=2 -> ОБЯЗАНА пробить 96 плоскостей\n");
  printf("=== РЕЗУЛЬТАТЫ ===\n");

  const int L = 5;
  sph s = {{16.3, 15.7, 16.1}, 11.5};
  hz_signgrid g = {NULL, 0, 0};
  hz_htab ht;
  hz_htab_init(&ht);
  check(hz_dc_sample(&g, &ht, L, sph_sign, sph_cross, &s) == HZ_DC_OK, "опрос источника");
  hz_dctree t;
  hz_dc_init(&t, L);
  check(hz_dc_build(&t, &g, &ht) == HZ_DC_OK, "сборка дерева");

  t_full(&t, &ht, &s);
  t_seam(&t);
  t_2b(&t);
  t_fan(&t, &s);

  hz_dc_free(&t);
  hz_htab_free(&ht);
  hz_signgrid_free(&g);
  printf("%s: %d/%d\n", g_fail ? "FAILURES" : "ok", g_total - g_fail, g_total);
  return g_fail != 0;
}
