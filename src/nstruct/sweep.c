/* sweep.c — свип на новых массивах: §835 (скалярный луч), §836-§838
 * (колонки/линии, материализованный обход), §841 (гигиена цикла).
 * См. sweep.h.
 */
#include "sweep.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HZ_SW_W (4.0 * 3.14159265358979323846 / 26.0) /* вес направления, изотропия 4π/N */
#define HZ_SW_SLABBITS 21 /* слэб и позиция в ключе; оси разумных сцен ≤ 2^21 */
#define HZ_SW_ND 26
#define HZ_SW_FOURPI (4.0 * 3.14159265358979323846)
#define HZ_SW_GL_EPS 1e-15 /* сходимость Ньютона для узлов Гаусса: машинный предел double */
#define HZ_SW_GL_ITMAX 100 /* предел итераций Ньютона: узел Гаусса сходится за ~5 */
#define HZ_SW_NPHI_MAX 96  /* Чебышёв по φ: шаг ≥ 3.75°, дальше вычислительно бессмысленно */
#define HZ_SW_NMU_MAX 64   /* Гаусс по μ: nd ≤ 96·64 = 6144 — потолок лестницы §843 */

/* знаки 26 направлений: 6 осей, 12 рёбер (√2), 8 углов (√3); ЦЕЛЫЕ —
 * сравнение с нулём без float-equal, нормировка отдельной таблицей */
static const int sw_sv[HZ_SW_ND][3] = {
    {1, 0, 0},  {-1, 0, 0},  {0, 1, 0},   {0, -1, 0},  {0, 0, 1},   {0, 0, -1}, {1, 1, 0},
    {1, -1, 0}, {-1, 1, 0},  {-1, -1, 0}, {1, 0, 1},   {1, 0, -1},  {-1, 0, 1}, {-1, 0, -1},
    {0, 1, 1},  {0, 1, -1},  {0, -1, 1},  {0, -1, -1}, {1, 1, 1},   {1, 1, -1}, {1, -1, 1},
    {-1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}};
static const double sw_norm[HZ_SW_ND] = {1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.4142135623730951,
                                         1.7320508075688772,
                                         1.7320508075688772,
                                         1.7320508075688772,
                                         1.7320508075688772,
                                         1.7320508075688772,
                                         1.7320508075688772,
                                         1.7320508075688772,
                                         1.7320508075688772};

typedef struct {
  int64_t key;  /* линия (А1481) */
  int64_t key2; /* ход луча внутри линии (§843/А1508) */
  int32_t li;
} hz_sw_pair;

static int sw_cmp_pair(const void *a, const void *b) {
  const hz_sw_pair *x = (const hz_sw_pair *)a, *y = (const hz_sw_pair *)b;
  if (x->key != y->key) return (x->key > y->key) - (x->key < y->key);
  return (x->key2 > y->key2) - (x->key2 < y->key2);
}

/* §843: узлы и веса Гаусса—Лагранжа на [−1,1] (Ньютон по P_n, классика
 * gauleg). Чётные n; узлы упорядочены по возрастанию. */
static void sw_gauss_legendre(int n, double *x, double *w) {
  for (int i = 0; i < n / 2; i++) {
    double z = cos(M_PI * (i + 0.75) / (n + 0.5)); /* начальное приближение */
    double pp = 0.0;
    for (int it = 0; it < HZ_SW_GL_ITMAX; it++) {
      double p0 = 1.0, p1 = z;
      for (int k = 2; k <= n; k++) {
        double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / (double)k;
        p0 = p1;
        p1 = p2;
      }
      /* P'_n через связующее тождество Лагранжа */
      pp = (double)n * (z * p1 - p0) / (z * z - 1.0);
      double dz = p1 / pp;
      z -= dz;
      if (fabs(dz) <= HZ_SW_GL_EPS) break;
    }
    x[i] = -z;
    w[i] = 2.0 / ((1.0 - z * z) * pp * pp);
    x[n - 1 - i] = z;
    w[n - 1 - i] = w[i];
  }
}

/* §843: таблица направлений. 6/26 — легаси (те же om посимвольно: om = sv/norm,
 * вес 4π/N, как HZ_SW_W); N_φ*100+N_μ — Чебышёв по φ (равные узлы, центр
 * ячейки, вес 2π/N_φ) × Гаусс—Легандр по μ (чётные N_μ). Σw = 4π точно. */
int hz_sw_dir_table(int ndirs, hz_sw_dir **out, int *nd_out) {
  hz_sw_dir *tab = NULL;
  int nd, i;
  if (!out || !nd_out || ndirs <= 0) return 1;
  if (ndirs == 6 || ndirs == 26) {
    nd = (ndirs == 26) ? HZ_SW_ND : 6;
    tab = (hz_sw_dir *)malloc((size_t)nd * sizeof *tab);
    if (!tab) return 2;
    for (i = 0; i < nd; i++) {
      tab[i].om[0] = (double)sw_sv[i][0] / sw_norm[i];
      tab[i].om[1] = (double)sw_sv[i][1] / sw_norm[i];
      tab[i].om[2] = (double)sw_sv[i][2] / sw_norm[i];
      tab[i].w = HZ_SW_FOURPI / (double)nd;
    }
  } else {
    int nphi = ndirs / 100, nmu = ndirs % 100;
    double *gx, *gw;
    if (nphi < 4 || nphi > HZ_SW_NPHI_MAX || nmu < 2 || nmu > HZ_SW_NMU_MAX || (nmu & 1) ||
        nphi * 100 + nmu != ndirs)
      return 1;
    nd = nphi * nmu;
    gx = (double *)malloc((size_t)nmu * sizeof *gx);
    gw = (double *)malloc((size_t)nmu * sizeof *gw);
    tab = (hz_sw_dir *)malloc((size_t)nd * sizeof *tab);
    if (!gx || !gw || !tab) {
      free(gx);
      free(gw);
      free(tab);
      return 2;
    }
    sw_gauss_legendre(nmu, gx, gw);
    for (i = 0; i < nphi; i++) {
      double phi = 2.0 * M_PI * ((double)i + 0.5) / (double)nphi;
      double sp = sin(phi), cp = cos(phi);
      for (int m = 0; m < nmu; m++) {
        double s = sqrt(1.0 - gx[m] * gx[m]);
        hz_sw_dir *e = &tab[i * nmu + m];
        e->om[0] = sp * s;
        e->om[1] = cp * s;
        e->om[2] = gx[m];
        e->w = (2.0 * M_PI / (double)nphi) * gw[m];
      }
    }
    free(gx);
    free(gw);
  }
  *out = tab;
  *nd_out = nd;
  return 0;
}

/* линия и слэб клетки id на направление со знаками sg (А1481; §843: sg из
 * знака om — у продуктовой квадратуры все компоненты ненулевые). */
static void sw_line_slab(const hz_pyr *py, const int sg[3], int64_t id, int64_t *line,
                         int64_t *slab) {
  const int *sv = sg;
  int64_t ii[3], len[3];
  int na[3] = {0, 0, 0}, nnz = 0, ax;
  len[0] = py->nx;
  len[1] = py->ny;
  len[2] = py->nz;
  ii[0] = id % py->nx;
  ii[1] = (id / py->nx) % py->ny;
  ii[2] = id / ((int64_t)py->nx * py->ny);
  for (ax = 0; ax < 3; ax++)
    if (sv[ax] != 0) na[nnz++] = ax;
  if (nnz == 1) {
    int c = (na[0] + 1) % 3, dd = (na[0] + 2) % 3;
    *line = ii[c] + len[c] * ii[dd];
    *slab = (int64_t)sv[na[0]] * ii[na[0]];
  } else if (nnz == 2) {
    int a = na[0], b = na[1], c = 3 - a - b;
    if (c < 0) c = 0;
    if (c > 2) c = 2;
    int64_t o1 = (int64_t)sv[b] * ii[b] - (int64_t)sv[a] * ii[a] + len[a];
    *line = o1 + (len[a] + len[b]) * ii[c];
    *slab = (int64_t)sv[a] * ii[a] + (int64_t)sv[b] * ii[b];
  } else {
    int a = na[0], b = na[1], c = na[2];
    int64_t o1 = (int64_t)sv[b] * ii[b] - (int64_t)sv[a] * ii[a] + len[a];
    int64_t o2 = (int64_t)sv[c] * ii[c] - (int64_t)sv[b] * ii[b] + len[b];
    *line = o1 + (len[a] + len[b]) * o2;
    *slab = (int64_t)sv[a] * ii[a] + (int64_t)sv[b] * ii[b] + (int64_t)sv[c] * ii[c];
  }
}

/* ПОХОДНЫЙ ЛИСТ направления (§841): cntot и T предвычислены — итерация
 * однопроходна. Только кусковые листья. */
typedef struct {
  int32_t npcs;
  int64_t first; /* офсет в wcn/wn/kd/wpid */
  int64_t line;  /* линия; смена — тёмный вход луча */
  int64_t cid;   /* §851: линейный id клетки листа (метки компонент) */
  double cntot;  /* Σ area·|ω·n| по кускам листа */
} hz_sw_wleaf;

typedef struct {
  hz_sw_wleaf *leaf;
  double *wcn, *wn, *kd;
  int32_t *wpid;
  int32_t n; /* походных листьев (с кусками) */
} hz_sw_walk;

static void sw_free_walk(hz_sw_walk *w) {
  free(w->leaf);
  free(w->wcn);
  free(w->wn);
  free(w->kd);
  free(w->wpid);
  memset(w, 0, sizeof *w);
}

/* §843/А1508: координата ХОДА луча внутри линии. Возрастает вдоль om:
 * для sg>0 — растущий индекс, для sg<0 — убывающий (len−1−ii); нулевая
 * компонента вдоль линии постоянна. Сортировка по (линия, ход) делает
 * порядок обхода полным и направленным по лучу для ЛЮБОГО знака om. */
static int64_t sw_travel_rank(const hz_pyr *py, const int sg[3], int64_t id) {
  int64_t ii[3], len[3], u = 0;
  int ax;
  len[0] = py->nx;
  len[1] = py->ny;
  len[2] = py->nz;
  ii[0] = id % py->nx;
  ii[1] = (id / py->nx) % py->ny;
  ii[2] = id / ((int64_t)py->nx * py->ny);
  for (ax = 0; ax < 3; ax++)
    u += (sg[ax] >= 0) ? ii[ax] : len[ax] - 1 - ii[ax];
  return u;
}

/* §845: поход ОБЪЁМНОГО ФРОНТА — порядок листьев из МАРША §834 (ход луча),
 * без линии и без тёмного входа на смене ключа: тёмный вход только на
 * границе домена. vindex — рабочий буфер [nleaf] (ранг марша, 1-based),
 * bits — буфер посещений [nleaf+63]/64. Линия в записях не заполняется (0)
 * — итерация mode 2 её не читает; сортировка по ОДНОМУ ключу (ранг марша). */
static int64_t sw_inv_acc; /* §845-бис: сумма инверсий марша за построение походов */

/* §845-в: СБОРКА похода из УПОРЯДОЧЕННОГО списка кусковых листьев (lis[n])
 * — общая для обоих строителей: предвычисление wcn/wn/kd/wpid и cntot. */
static int sw_walk_assemble(const hz_pyr *py, const double *om, const double *area,
                            const double *nrm, const double *kd, const int32_t *lis, int32_t n,
                            hz_sw_walk *w) {
  int32_t k;
  int64_t pos = 0;
  w->leaf = (hz_sw_wleaf *)malloc((size_t)py->nleaf * sizeof *w->leaf);
  w->wcn = (double *)malloc((size_t)py->nt * sizeof *w->wcn);
  w->wn = (double *)malloc((size_t)py->nt * sizeof *w->wn);
  w->kd = (double *)malloc((size_t)py->nt * sizeof *w->kd);
  w->wpid = (int32_t *)malloc((size_t)py->nt * sizeof *w->wpid);
  if (!w->leaf || !w->wcn || !w->wn || !w->kd || !w->wpid) {
    sw_free_walk(w);
    return 2;
  }
  w->n = n;
  for (k = 0; k < n; k++) {
    const hz_pyr_leaf *lf = &py->leaf[lis[k]];
    int32_t u;
    double cntot = 0.0;
    w->leaf[k].line = 0; /* линия в mode 2 не существует */
    w->leaf[k].cid = py->leaf_id[lis[k]];
    w->leaf[k].npcs = lf->npcs;
    w->leaf[k].first = pos;
    for (u = 0; u < lf->npcs; u++) {
      int32_t p = py->csr[lf->pcs_first + u];
      double d0 = om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                  om[2] * nrm[3 * (int64_t)p + 2];
      double an = fabs(d0);
      w->wcn[pos] = area[p] * an;
      w->wn[pos] = an;
      w->kd[pos] = kd[p];
      w->wpid[pos] = p;
      cntot += w->wcn[pos];
      pos++;
    }
    w->leaf[k].cntot = cntot;
  }
  return 0;
}

static int sw_build_walk_march(const hz_pyr *py, const double *om, const double *area,
                               const double *nrm, const double *kd, int32_t *vindex, uint64_t *bits,
                               hz_sw_walk *w) {
  hz_pyr_march_stat ms;
  hz_sw_pair *pairs;
  int32_t li, k, rc;
  hz_pyr_march(py, om, 0, &ms, bits, vindex);
  sw_inv_acc += ms.inversions;
  pairs = (hz_sw_pair *)malloc((size_t)py->nleaf * sizeof *pairs);
  if (!pairs) return 2;
  {
    int32_t np = 0;
    for (li = 0; li < py->nleaf; li++) {
      if (py->leaf[li].npcs == 0) continue; /* пустой лист — не участник переноса */
      pairs[np].key = (int64_t)vindex[li];  /* ранг марша — ход луча */
      pairs[np].key2 = 0;
      pairs[np].li = li;
      np++;
    }
    w->n = np;
  }
  qsort(pairs, (size_t)w->n, sizeof *pairs, sw_cmp_pair);
  {
    int32_t *lis = (int32_t *)malloc((size_t)w->n * sizeof *lis);
    if (!lis) {
      free(pairs);
      return 2;
    }
    for (k = 0; k < w->n; k++)
      lis[k] = pairs[k].li;
    free(pairs);
    rc = sw_walk_assemble(py, om, area, nrm, kd, lis, w->n, w);
    free(lis);
  }
  return rc;
}

/* §845-в: ПРЯМОЙ строитель — та же сортировка кусковых листьев по проекции
 * центра на ω, что даёт марш (pyr_slab), но БЕЗ обхода пирамиды и двоичных
 * поисков: O(nleaf) вычислений + одна сортировка. Инверсий нет по
 * построению (сортируем по самой величине). */
typedef struct {
  double s;
  int32_t li;
} hz_sw_sitem;

static int sw_cmp_sitem(const void *a, const void *b) {
  const hz_sw_sitem *x = (const hz_sw_sitem *)a, *y = (const hz_sw_sitem *)b;
  if (x->s > y->s) return 1;
  if (x->s < y->s) return -1;
  return (x->li > y->li) - (x->li < y->li); /* полный порядок: детерминизм */
}

static int sw_build_walk_sort(const hz_pyr *py, const double *om, const double *area,
                              const double *nrm, const double *kd, hz_sw_walk *w) {
  hz_sw_sitem *items = (hz_sw_sitem *)malloc((size_t)py->nleaf * sizeof *items);
  int32_t li, k, np = 0, rc;
  if (!items) return 2;
  for (li = 0; li < py->nleaf; li++) {
    int64_t id;
    double c[3];
    if (py->leaf[li].npcs == 0) continue;
    id = py->leaf_id[li];
    c[0] = py->lo[0] + ((double)(id % py->nx) + 0.5) * py->cell;
    c[1] = py->lo[1] + ((double)((id / py->nx) % py->ny) + 0.5) * py->cell;
    c[2] = py->lo[2] + ((double)(id / ((int64_t)py->nx * py->ny)) + 0.5) * py->cell;
    items[np].s = om[0] * c[0] + om[1] * c[1] + om[2] * c[2];
    items[np].li = li;
    np++;
  }
  qsort(items, (size_t)np, sizeof *items, sw_cmp_sitem);
  {
    int32_t *lis = (int32_t *)malloc((size_t)np * sizeof *lis);
    if (!lis) {
      free(items);
      return 2;
    }
    for (k = 0; k < np; k++)
      lis[k] = items[k].li;
    free(items);
    rc = sw_walk_assemble(py, om, area, nrm, kd, lis, np, w);
    free(lis);
  }
  return rc;
}

/* §841: поход строится ПРЯМО — фильтр кусковых листьев, сортировка ТОЛЬКО
 * их (полная сортировка 14.9M пустых не нужна вовсе), cntot/T сразу.
 * §843/А1508: сортировка по (линия, ход луча) — раньше порядок внутри
 * линии доставался Morton-порядком листьев и был ПРОТИВ луча для
 * направлений с отрицательными компонентами. */
static int sw_build_walk(const hz_pyr *py, const int sg[3], const double *om, const double *area,
                         const double *nrm, const double *kd, hz_sw_walk *w) {
  int32_t li, k;
  int64_t pos = 0;
  hz_sw_pair *pairs = (hz_sw_pair *)malloc((size_t)py->nleaf * sizeof *pairs);
  w->leaf = (hz_sw_wleaf *)malloc((size_t)py->nleaf * sizeof *w->leaf);
  w->wcn = (double *)malloc((size_t)py->nt * sizeof *w->wcn);
  w->wn = (double *)malloc((size_t)py->nt * sizeof *w->wn);
  w->kd = (double *)malloc((size_t)py->nt * sizeof *w->kd);
  w->wpid = (int32_t *)malloc((size_t)py->nt * sizeof *w->wpid);
  if (!pairs || !w->leaf || !w->wcn || !w->wn || !w->kd || !w->wpid) {
    free(pairs);
    sw_free_walk(w);
    return 2;
  }
  {
    int32_t np = 0;
    for (li = 0; li < py->nleaf; li++) {
      if (py->leaf[li].npcs == 0) continue; /* пустой лист — не участник переноса */
      sw_line_slab(py, sg, py->leaf_id[li], &pairs[np].key, &(int64_t){0});
      pairs[np].key2 = sw_travel_rank(py, sg, py->leaf_id[li]);
      pairs[np].li = li;
      np++;
    }
    w->n = np;
  }
  qsort(pairs, (size_t)w->n, sizeof *pairs, sw_cmp_pair);
  for (k = 0; k < w->n; k++) {
    const hz_pyr_leaf *lf = &py->leaf[pairs[k].li];
    int32_t u;
    double cntot = 0.0;
    {
      int64_t slab_dummy;
      sw_line_slab(py, sg, py->leaf_id[pairs[k].li], &w->leaf[k].line, &slab_dummy);
    }
    w->leaf[k].npcs = lf->npcs;
    w->leaf[k].first = pos;
    for (u = 0; u < lf->npcs; u++) {
      int32_t p = py->csr[lf->pcs_first + u];
      double d0 = om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                  om[2] * nrm[3 * (int64_t)p + 2];
      double an = fabs(d0);
      w->wcn[pos] = area[p] * an;
      w->wn[pos] = an;
      w->kd[pos] = kd[p];
      w->wpid[pos] = p;
      cntot += w->wcn[pos];
      pos++;
    }
    w->leaf[k].cntot = cntot;
  }
  free(pairs);
  return 0;
}
/* §852: ЛУЧЕВОЙ СВИП ОДНОГО НАПРАВЛЕНИЯ С ТОЧНЫМ ПЕРЕСЕЧЕНИЕМ. Трубки =
 * линии сетки вдоль om; стартуют тёмными на входных гранях; пустота
 * транслирует L точно; в клетке с кусками — БЛИЖАЙШЕЕ пересечение линии
 * с треугольником на отрезке [s_enter, s_exit]: депонирует свой инфлюкс
 * и переизлучает СВОЙ радианс. Непрозрачность тоньше/толще клетки
 * выражается геометрией — клеточных коэффициентов нет (А1557). */
typedef struct {
  double absorbed, emitted, recycled, lost;
  int64_t nvisit, ndep, ncell;
} sw_line_acc;

/* Мёллер—Трумбор: t пересечения луча (o,d) с треугольником tv[9]; RAW —
 * без отсечения позади луча (§852: стенка на границе домена даёт t<0 от
 * центра входной клетки), фильтр — отрезком марша */
static double sw_ray_tri_raw(const double o[3], const double d[3], const double tv[9]) {
  double e1[3], e2[3], p[3], t[3], q[3], det, uu, vv, tt;
  int ax;
  for (ax = 0; ax < 3; ax++) {
    e1[ax] = tv[3 + ax] - tv[ax];
    e2[ax] = tv[6 + ax] - tv[ax];
    t[ax] = o[ax] - tv[ax];
  }
  p[0] = d[1] * e2[2] - d[2] * e2[1];
  p[1] = d[2] * e2[0] - d[0] * e2[2];
  p[2] = d[0] * e2[1] - d[1] * e2[0];
  det = e1[0] * p[0] + e1[1] * p[1] + e1[2] * p[2];
  if (fabs(det) < 1e-30) return -1.0;
  uu = (t[0] * p[0] + t[1] * p[1] + t[2] * p[2]) / det;
  if (uu < 0.0 || uu > 1.0) return -1.0;
  q[0] = t[1] * e1[2] - t[2] * e1[1];
  q[1] = t[2] * e1[0] - t[0] * e1[2];
  q[2] = t[0] * e1[1] - t[1] * e1[0];
  vv = (d[0] * q[0] + d[1] * q[1] + d[2] * q[2]) / det;
  if (vv < 0.0 || uu + vv > 1.0) return -1.0;
  tt = (e2[0] * q[0] + e2[1] * q[1] + e2[2] * q[2]) / det;
  return tt;
}

/* классический: только впереди луча (§851) */
static double sw_ray_tri(const double o[3], const double d[3], const double tv[9]) {
  double tt = sw_ray_tri_raw(o, d, tv);
  return tt > 0.0 ? tt : -1.0;
}

static void sw_line_sweep_dir(const hz_pyr *py, const double om[3], double w_d, double le,
                              const double *kd, const double *area, const double *nrm,
                              const double *Eprev, double *Ed, double rho_ovr,
                              const double *trivert, int noprop, int32_t *stamp, int32_t mark,
                              sw_line_acc *acc) {
  (void)area; /* площадь нужна была клеточной модели; точное пересечение её не читает */
  int32_t ax;
  int64_t f[3];
  for (ax = 0; ax < 3; ax++) {
    int64_t n = ax == 0 ? py->nx : (ax == 1 ? py->ny : py->nz);
    if (fabs(om[ax]) < 1e-30) {
      f[ax] = -1; /* ось не входная */
      continue;
    }
    f[ax] = om[ax] > 0.0 ? 0 : n - 1;
  }
  for (ax = 0; ax < 3; ax++) {
    int64_t b1, b2;
    int32_t bxA, bxB;
    int64_t nb1, nb2;
    if (f[ax] < 0) continue;
    bxA = (ax + 1) % 3;
    bxB = (ax + 2) % 3;
    nb1 = bxA == 0 ? py->nx : (bxA == 1 ? py->ny : py->nz);
    nb2 = bxB == 0 ? py->nx : (bxB == 1 ? py->ny : py->nz);
    for (b1 = 0; b1 < nb1; b1++)
      for (b2 = 0; b2 < nb2; b2++) {
        int64_t cc[3], cell[3], stepv[3];
        double tmax[3], tdelta[3], L = 0.0, org[3], s_enter = 0.0;
        int32_t a2;
        cc[ax] = f[ax];
        cc[bxA] = b1;
        cc[bxB] = b2;
        cell[0] = cc[0];
        cell[1] = cc[1];
        cell[2] = cc[2];
        org[0] = py->lo[0] + ((double)cc[0] + 0.5) * py->cell;
        org[1] = py->lo[1] + ((double)cc[1] + 0.5) * py->cell;
        org[2] = py->lo[2] + ((double)cc[2] + 0.5) * py->cell;
        if (stamp[cell[0] + py->nx * (cell[1] + py->ny * cell[2])] == mark) continue;
        for (a2 = 0; a2 < 3; a2++) {
          if (fabs(om[a2]) < 1e-30) {
            tmax[a2] = 1e30;
            tdelta[a2] = 0.0;
            stepv[a2] = 0;
          } else {
            double boundary = om[a2] > 0.0 ? py->lo[a2] + ((double)cell[a2] + 1.0) * py->cell
                                           : py->lo[a2] + (double)cell[a2] * py->cell;
            double pos = py->lo[a2] + ((double)cell[a2] + 0.5) * py->cell;
            stepv[a2] = om[a2] > 0.0 ? 1 : -1;
            tdelta[a2] = py->cell / fabs(om[a2]);
            tmax[a2] = (boundary - pos) / om[a2];
          }
        }
        for (;;) {
          int64_t id = cell[0] + py->nx * (cell[1] + py->ny * cell[2]);
          int32_t pos = hz_pyr_leaf_pos(py, id);
          double tn, s_exit;
          acc->ncell++;
          stamp[id] = mark;
          tn = tmax[0];
          if (tmax[1] < tn) tn = tmax[1];
          if (tmax[2] < tn) tn = tmax[2];
          s_exit = tn > 1e29 ? 1e30 : tn;
          if (pos >= 0) {
            /* §852: ТОЧНОЕ пересечение линии с треугольниками клетки на
             * отрезке [s_enter, s_exit]: ближайший поглощает трубку,
             * депонирует инфлюкс и переизлучает СВОЙ радианс (А1557). */
            const hz_pyr_leaf *lf = &py->leaf[pos];
            double tbest = -1.0;
            int32_t pbest = -1, u;
            double an_best = 0.0;
            acc->nvisit++;
            for (u = 0; u < lf->npcs; u++) {
              int32_t p = py->csr[lf->pcs_first + u];
              double tt = sw_ray_tri(org, om, trivert + 9 * (int64_t)py->pcs[p].tri);
              if (tt > s_enter && tt <= s_exit && (tbest < 0.0 || tt < tbest)) {
                tbest = tt;
                pbest = p;
                an_best = fabs(om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                               om[2] * nrm[3 * (int64_t)p + 2]);
              }
            }
            if (pbest >= 0) {
              double rho = rho_ovr < 0 ? kd[pbest] : rho_ovr;
              if (L > 0.0) {
                acc->ndep++;
                Ed[pbest] += w_d * L * an_best;
                acc->absorbed += w_d * L * an_best;
              }
              acc->emitted += w_d * le * an_best;
              acc->recycled += w_d * rho * Eprev[pbest] / (2.0 * M_PI) * an_best;
              L = le + rho * Eprev[pbest] / (2.0 * M_PI); /* радианс СВОГО куска */
            }
            if (noprop) L = 0.0; /* НК: фронт не переносится */
          }
          if (s_exit > 1e29) break;
          s_enter = s_exit;
          if (tmax[0] <= tmax[1] && tmax[0] <= tmax[2]) {
            tmax[0] += tdelta[0];
            cell[0] += stepv[0];
          } else if (tmax[1] <= tmax[2]) {
            tmax[1] += tdelta[1];
            cell[1] += stepv[1];
          } else {
            tmax[2] += tdelta[2];
            cell[2] += stepv[2];
          }
          if (cell[0] < 0 || cell[1] < 0 || cell[2] < 0 || cell[0] >= py->nx || cell[1] >= py->ny ||
              cell[2] >= py->nz)
            break;
        }
      }
  }
}

/* ---- §852: МАТЕРИАЛЬНЫЙ ФРОНТ НА ПИРАМИДЕ (mode 3) ------------------------ */

/* Постановка (правлена по А1563…А1572): трубки = линии базовой сетки вдоль
 * om (как §851); марш — по узлам пирамиды. Узел ПУСТ — прыжок одним шагом
 * по tmax узла (А1569); МАТЕРИАЛЕН (max ℓ_p поддерева ≤ сторона узла) —
 * локальное взаимодействие на СВОЁМ размере клетки; ПРОМЕЖУТОЧНЫЙ — спуск,
 * продолжение DDA на детском уровне. L — РАДИАНС: сечение трубки при
 * спуске НЕ ДЕЛИТСЯ, Φ = w·L·csec с сечением БАЗОВОЙ клетки постоянен,
 * доля ребёнка 0/1 по геометрии входа, слияния нет (А1563). Клеточный
 * перехват T = 1−cntot/csec — ТОЛЬКО под предикатом «bbox куска ⊆ клетка
 * узла», вне предиката — точное пересечение, fail closed (А1566). Крат-
 * ность 1 на (кусок, направление) — штампом (А1564; чинит и самозасветку
 * тонкой оболочки vc=1: двусторонний ray-triangle больше не депонирует
 * радианс оболочки на неё же при выходе трубки). G6: на каждом спуске
 * parent_flux = Σchild_flux — по РОДИТЕЛЬСКОЙ ДОЛЕ a носителя (эмиссия —
 * новая энергия, в тождество не входит); нарушение — счётчик, на чистом
 * прогоне 0. Носитель — пара (a, b): a — родительская доля радианса,
 * b — порождённая (эмиссия кусков); депозиты делятся в той же пропорции. */
/* §863/шаг 2: статичные группы агрегации (владелец — hz_sw_run). Группа =
 * (материал × узел): представитель ray-tri, константы A/le/ρ предвычислены,
 * энергия депозита копится на группу (O(1) на попадание) и раздаётся членам
 * после итерации. Eprev-среднее пересчитывается раз в итерацию (SE). */
typedef struct {
  int64_t ng, ngcap;
  int32_t *grep;            /* представитель группы */
  double *gA, *gleA, *grA;  /* площадь, средние le и ρ (статичные) */
  double *gSE, *gacc, *gwt; /* ΣEprev·a, Σedep, Σedep·(1−cosθ) — на итерацию */
  int64_t *ggoff;
  int32_t *ggcnt, *gmem;
  int64_t gmem_n, gmem_cap;
} sw_agg;

typedef struct {
  const hz_pyr *py;
  const hz_sw_opts *o;
  const double *area, *nrm, *kd, *Eprev;
  double *Ed;
  const double *om;
  double org[3]; /* центр входной клетки трубки */
  double w_d, le, axcos;
  const double *lep; /* А1576: per-piece эмиссия (Ke); NULL — глобальный le */
  int noprop, tau0, ax;
  int64_t *pstamp; /* штамп (трубка × кусок × направление × итерация), А1564: гасит
                    * повторный депозит ОДНОЙ трубки на ту же оболочку (вход и
                    * выход тонкой оболочки); разные трубки направления
                    * депонируют каждую по-своему — пространственный сбор */
  int64_t pkey;    /* значение штампа текущей трубки: mark<<32 | номер трубки */
  int64_t ntube;   /* счётчик трубок направления (в ключе штампа) */

  int32_t mark;
  /* bbox-индекс кусков по листьям (А1566: кусок владеет ОДНОЙ клеткой по
   * центроиду, а пересекать трубку его треугольник может в СОСЕДНЕЙ —
   * взаимодействие обязано идти по bbox-перекрытию, не по владению) */
  const int32_t *bstart; /* [nleaf+1] */
  const int32_t *bpids;  /* слоты кусков */
  int32_t *pbuf;         /* куски материального узла (растёт по нужде) */
  int64_t pbufcap;
  /* §863: кэш списков кусков УЗЛОВ — gather обходил поддерево на каждую
   * трубку; теперь список узла собирается ОДИН РАЗ и живёт в nstore.
   * nstart/nlen/nstore/nstore_n принадлежат hz_sw_run (переживают fc). */
  int64_t *nstart;   /* [ncluster]: позиция в nstore или -1 — не собран */
  int32_t *nlen;     /* [ncluster]: длина списка */
  int32_t *nstore;   /* магазин списков */
  int64_t *nstore_n; /* указатель на счётчик заполнения магазина */
  int64_t nstore_cap;
  int32_t ncluster;
  int64_t noff[256];   /* смещение уровня в нумерации кластеров узлов */
  int32_t ncluster0;   /* nleaf: кластеры листов идут первыми */
  int agg;             /* §863: ключ o->agg */
  int agg_list;        /* текущий список — узловой (отсортирован, годен для групп) */
  sw_agg *ag;          /* §863/шаг 2: таблица групп (NULL — путь не активен) */
  double *abuf, *wbuf; /* скретч взаимодействия: без аллокаций на событие */
  int64_t abufcap;
  double hi[3]; /* верх сцены (клетки уровней шире домена на нечётных сетках) */
  /* §864/Б1: кэш размеров уровней (статичны на прогон; иначе level_dims —
   * 22–37 млн вызовов). Владелец — hz_sw_run, fc только ссылается. */
  int64_t (*fld)[3];
  uint8_t *fldok;
  double *pbuf_t;
  int32_t *pbuf_p;
  int64_t pbuf_n, pbuf_cap;
  double lost;       /* поток за границей домена, w·csec-единицы */
  int64_t ncellbase; /* базовых клеток полным DDA («до», прибор А1569) */
  int cbase, cfront; /* прибор считается на первой итерации (геометрия статична) */
  /* аккумуляторы */
  double depA; /* Σ родительской доли, депонированной кускам (L-единицы) */
  double absorbed, emitted, recycled;
  int64_t g6viol, njump, nmat, ndesc, ncellfront, nstamp, nlostseg, ndep;
} front_ctx;

static int32_t sw_group_make(front_ctx *fc, int64_t g0, int64_t g1); /* §863 */

/* отрезок [h0,h1] луча (org, om) в коробке, пересечённый с [tin,tout];
 * 0 — пусто */
static int front_box_seg(const double blo[3], const double bhi[3], const double org[3],
                         const double om[3], double tin, double tout, double *h0, double *h1) {
  double t0 = tin, t1 = tout;
  int ax;
  for (ax = 0; ax < 3; ax++) {
    if (fabs(om[ax]) < 1e-30) {
      if (org[ax] < blo[ax] || org[ax] > bhi[ax]) return 0;
      continue;
    }
    {
      double ta = (blo[ax] - org[ax]) / om[ax];
      double tb = (bhi[ax] - org[ax]) / om[ax];
      if (ta > tb) {
        double tt = ta;
        ta = tb;
        tb = tt;
      }
      if (ta > t0) t0 = ta;
      if (tb < t1) t1 = tb;
    }
  }
  /* t1 == t0 допустимо: вырожденная (плоская) клетка фантомной колонки —
   * стенка на mesh-границе сетки обязана быть достижима для марша (§852) */
  if (t1 < t0) return 0;
  *h0 = t0;
  *h1 = t1;
  return 1;
}

static double front_rho(const front_ctx *fc, int32_t p) {
  return fc->o->rho < 0 ? fc->kd[p] : fc->o->rho;
}

/* А1576: эмиссия куска — per-piece (Ke) плюс глобальный le; при lep=NULL
 * (умолчание) поведение побитово прежнее */
static double front_le(const front_ctx *fc, int32_t p) {
  return fc->le + (fc->lep ? (double)fc->lep[p] : 0.0);
}

static double front_cos(const front_ctx *fc, int32_t p) {
  const double *om = fc->om;
  return fabs(om[0] * fc->nrm[3 * (int64_t)p] + om[1] * fc->nrm[3 * (int64_t)p + 1] +
              om[2] * fc->nrm[3 * (int64_t)p + 2]);
}

/* §862: приращение дробного аккумулятора детальности Δ(материал,угол):
 * темнее материал (ρ) и скользящее падение (1−cosθ, cosθ — к НОРМАЛИ куска)
 * — быстрее набор этажа. Нет lpacc — нет операции (битово прежний мир). */
static void sw_accum(front_ctx *fc, int32_t p, double edep) {
  if (!fc->o->lpacc) return;
  {
    double rho = front_rho(fc, p);
    double ct = front_cos(fc, p) / (2.0 * fc->area[p]); /* |cos| к нормали куска */
    double d;
    switch (fc->o->accum_mode) {
    case 1:
      d = edep * rho * (1.0 - ct);
      break;
    case 2:
      d = edep * rho;
      break;
    case 3:
      d = edep * (1.0 - ct);
      break;
    case 4:
      d = edep;
      break;
    default:
      d = rho * (1.0 - ct);
      break; /* §862: форма за событие */
    }
    fc->o->lpacc[p] += fc->o->cdelta * d;
  }
  if (fc->o->lphits) fc->o->lphits[p] += 1.0; /* §862-диаг: ранжир (а) */
}

/* §864/Б1: размеры уровня из кэша (статичны на прогон) */
static void front_ldims(front_ctx *fc, int32_t l, int64_t d[3]) {
  if (!fc->fldok[l]) {
    hz_pyr_level_dims(fc->py, l, fc->fld[l]);
    fc->fldok[l] = 1;
  }
  d[0] = fc->fld[l][0];
  d[1] = fc->fld[l][1];
  d[2] = fc->fld[l][2];
}

/* коробка узла (уровень l, позиция pos; l<0 — лист) с ЗАЖИМОМ в сцену:
 * на нечётных сетках клетки уровней шире домена (А1567 — за доменом
 * сегментов не остаётся, остаток носителя уходит в lost) */
static void front_node_box(front_ctx *fc, int32_t l, int32_t pos, double blo[3],
                           double bhi[3]) { /* fc не const: кэш размеров §864/Б1 */
  const hz_pyr *py = fc->py;
  int64_t d[3], id;
  double side;
  int ax;
  if (l < 0) {
    d[0] = py->nx;
    d[1] = py->ny;
    d[2] = py->nz;
    side = py->cell;
    id = py->leaf_id[pos];
  } else {
    front_ldims(fc, l, d); /* §864/Б1: кэш размеров уровней */
    side = py->cell * (double)((int64_t)1 << (l + 1));
    id = py->lev[l][pos].id;
  }
  blo[0] = py->lo[0] + (double)(id % d[0]) * side;
  blo[1] = py->lo[1] + (double)((id / d[0]) % d[1]) * side;
  blo[2] = py->lo[2] + (double)(id / (d[0] * d[1])) * side;
  for (ax = 0; ax < 3; ax++) {
    bhi[ax] = blo[ax] + side;
    if (bhi[ax] > fc->hi[ax]) bhi[ax] = fc->hi[ax];
  }
}

/* собрать куски поддерева (уровень l, позиция pos; l<0 — лист) в fc->pbuf */
static int front_gather(front_ctx *fc, int32_t l, int32_t pos, int32_t *n) {
  const hz_pyr *py = fc->py;
  typedef struct {
    int32_t l, pos;
  } fp;
  fp stk[2048]; /* глубина·8: nlev ≤ 254 — стек на кадре C, без аллокации */
  int64_t sp = 0;
  if ((py->nlev + 2) * 8 > 2048) return 2;
  int32_t cnt = 0;
  int rc = 0;
  *n = 0;
  stk[sp].l = l;
  stk[sp].pos = pos;
  sp++;
  while (sp > 0) {
    int32_t cl, cpos, u;
    sp--;
    cl = stk[sp].l;
    cpos = stk[sp].pos;
    if (cl < 0) {
      int32_t nb = fc->bstart[cpos + 1] - fc->bstart[cpos];
      const int32_t *bp = fc->bpids + fc->bstart[cpos];
      if (nb > 0) {
        if (cnt + nb > fc->pbufcap) {
          int64_t nc = fc->pbufcap ? fc->pbufcap * 2 : 64;
          int32_t *nb2;
          while (nc < cnt + nb)
            nc *= 2;
          nb2 = (int32_t *)realloc(fc->pbuf, (size_t)nc * sizeof *nb2);
          if (!nb2) {
            rc = 2;
            break;
          }
          fc->pbuf = nb2;
          fc->pbufcap = nc;
        }
        for (u = 0; u < nb; u++)
          fc->pbuf[cnt++] = bp[u];
      }
      continue;
    }
    {
      int64_t pd[3], cd[3], kid = py->lev[cl][cpos].id, ci[3];
      int cx, cy, cz;
      int32_t cl2 = cl - 1;
      hz_pyr_level_dims(py, cl, pd);
      if (cl2 < 0) {
        cd[0] = py->nx;
        cd[1] = py->ny;
        cd[2] = py->nz;
      } else
        hz_pyr_level_dims(py, cl2, cd);
      ci[0] = kid % pd[0];
      ci[1] = (kid / pd[0]) % pd[1];
      ci[2] = kid / (pd[0] * pd[1]);
      for (cz = 0; cz < 2; cz++)
        for (cy = 0; cy < 2; cy++)
          for (cx = 0; cx < 2; cx++) {
            int64_t qx = 2 * ci[0] + cx, qy = 2 * ci[1] + cy, qz = 2 * ci[2] + cz;
            int64_t id = qx + cd[0] * (qy + cd[1] * qz);
            int32_t found;
            if (qx >= cd[0] || qy >= cd[1] || qz >= cd[2]) continue;
            if (cl2 < 0)
              found = hz_pyr_leaf_pos(py, id);
            else
              found = hz_pyr_node_pos(py, cl2, id);
            if (found < 0) continue; /* ПУСТ-ребёнок: кусков нет */
            stk[sp].l = cl2;
            stk[sp].pos = found;
            sp++;
          }
    }
  }
  if (rc) return rc;
  *n = cnt;
  return 0;
}

/* предикат А1566: bbox куска целиком внутри клетки (без допусков) */
static int front_contained(const double tribox[6], const double blo[3], const double bhi[3]) {
  int ax;
  for (ax = 0; ax < 3; ax++)
    if (!(tribox[ax] >= blo[ax]) || !(tribox[3 + ax] <= bhi[ax])) return 0;
  return 1;
}

/* взаимодействие трубки с материальной клеткой на отрезке [tin,tout].
 * ps/n — куски; blo/bhi — клетка. isleaf — отрезок листа (базовая клетка):
 * только там законен клеточный T-перехват (А1566 в окне «bbox ⊆ базовая
 * клетка»); на узловом отрезке (isleaf=0) — только точное пересечение,
 * А1573. Точные события (кроссирующие куски, при xi=1 — все) старше
 * клеточного перехвата: геометрия перекрывает модель. */
static void front_interact(front_ctx *fc, const int32_t *ps, int32_t n, const double blo[3],
                           const double bhi[3], double tin, double tout, double *a, double *b,
                           int isleaf) {
  const hz_pyr *py = fc->py;
  const int xi = fc->o->xint || !isleaf; /* узловой отрезок — всегда точно */
  double csec = py->cell * py->cell;     /* сечение БАЗОВОЙ клетки — инвариант А1563 */
  double csecnode = (bhi[1] - blo[1]) * (bhi[2] - blo[2]); /* сечение клетки узла: тень в T */
  double *an = NULL, *wt = NULL;
  double cntot = 0.0, swt = 0.0, Lsurf = 0.0, tbest = -1.0, anb = 0.0;
  int32_t u, pbest = -1;
  double ai = *a, bi = *b, Lin = ai + bi;
  int has_cont = 0;

  if (n <= 0) return;
  if (fc->tau0) { /* НК: слой не поглощает и не излучает — трубка не тронута */
    return;
  }
  fc->nmat++;
  if (n > fc->abufcap) {
    int64_t nc = fc->abufcap ? fc->abufcap * 2 : 64;
    double *na, *nw;
    while (nc < n)
      nc *= 2;
    na = (double *)realloc(fc->abuf, (size_t)nc * sizeof *na);
    if (!na) goto done;
    fc->abuf = na;
    nw = (double *)realloc(fc->wbuf, (size_t)nc * sizeof *nw);
    if (!nw) goto done;
    fc->wbuf = nw;
    fc->abufcap = nc;
  }
  an = fc->abuf;
  wt = fc->wbuf; /* нет памяти: сегмент проходится насквозь */
  for (u = 0; u < n; u++) {
    double cs = front_cos(fc, ps[u]);
    an[u] = fc->area[ps[u]] * cs; /* тень куска на сечении трубки */
    wt[u] = cs * an[u];           /* вес депозита: тень с косинусом */
  }
  /* 1. ТОЧНОЕ пересечение (А1566 fail closed; G1(a) при xint=1): ближайшее
   * unstamped-попадание кроссирующего куска — событие сегмента, §851-семантика. */
  for (u = 0; u < n; u++) {
    int32_t p = ps[u];
    int32_t tri = py->pcs[p].tri;
    double tt, bh0, bh1;
    if (!xi && front_contained(fc->o->tribox + 6 * (int64_t)tri, blo, bhi)) continue;
    if (fc->pstamp[p] == fc->pkey) {
      fc->nstamp++; /* кратность 1 на (кусок, направление) — А1564 */
      continue;
    }
    /* пре-фильтр bbox×сегмент (цена А1573 на узловом пути): треугольник ⊆
     * своего bbox, отсечение консервативно, результата не меняет */
    if (!front_box_seg(fc->o->tribox + 6 * (int64_t)tri, fc->o->tribox + 6 * (int64_t)tri + 3,
                       fc->org, fc->om, tin, tout, &bh0, &bh1))
      continue;
    tt = sw_ray_tri_raw(fc->org, fc->om, fc->o->trivert + 9 * (int64_t)tri);
    /* tt ≈ [tin,tout] с ДОПУСКОМ sl (единицы ulp арифметики ray-tri против
     * box-сегмента): стенка, лежащая РОВНО на конце сегмента (вход домена!),
     * иначе теряется — трубка остаётся тёмной и теряет ВЕСЬ депозит
     * (А1577: фальсификатор lev=5, E_half вдвое занижена). Двойной счёт
     * стыка гасится штампом (кусок × трубка) */
    {
      double sl = 1e-9 * (fabs(tin) + fabs(tout) + 1.0);
      if (tt >= tin - sl && tt <= tout + sl && (tbest < 0.0 || tt < tbest)) {
        tbest = tt;
        pbest = p;
        anb = front_cos(fc, p);
      }
    }
  }
  (void)anb;
  if (pbest >= 0) {
    double Lh;
    if (fc->noprop) {
      *a = 0.0;
      *b = 0.0;
      goto done;
    } /* НК: фронт не переносится */
    if (Lin > 0.0) fc->pstamp[pbest] = fc->pkey;
    if (Lin > 0.0) {
      /* Φ = w·L·csec — инвариант трубки (А1563); |cos| входит ЧЕРЕЗ ЧИСЛО
       * трубок, пересекающих кусок (A·cos/csec), поэтому G2 (ρ=0 → E=πLe)
       * держится при любой клетке */
      /* вклад трубки в ОБЛУЧЁННОСТЬ куска: Φ/A (мощность на его площадь) */
      fc->Ed[pbest] += fc->w_d * Lin * csec * fc->axcos / fc->area[pbest];
      fc->depA += ai; /* родительская доля трубки погашена целиком (G6) */
      fc->absorbed += fc->w_d * Lin * csec;
      fc->ndep++;
      Lh = front_le(fc, pbest) + front_rho(fc, pbest) * fc->Eprev[pbest] / (2.0 * M_PI);
      fc->recycled += fc->w_d * (Lh - fc->le) * csec;
    }
    fc->emitted += fc->w_d * front_le(fc, pbest) * csec;
    *a = 0.0;
    *b = front_le(fc, pbest) + front_rho(fc, pbest) * fc->Eprev[pbest] / (2.0 * M_PI);
    sw_accum(fc, pbest, fc->w_d * Lin * csec * fc->axcos / fc->area[pbest]); /* §862 */
    goto done;
  }
  /* 2. Клеточный перехват §846 — ТОЛЬКО контейнированные unstamped и ТОЛЬКО
   * на листовом отрезке (А1566+А1573: окно валидности «bbox ⊆ БАЗОВАЯ
   * клетка»; на узле T-перехват раздаёт сечение узла кускам, мимо которых
   * трубка прошла — расходимость лестницы lp>=2, замерено до правки:
   * cavity05 lp=2 E 3.43→6.05 за 6 итераций). Нормировка по wt:
   * Σ депозитов = f·вход (тождество G6), а на одной пластине f=1 депозит =
   * w·L·cos — совпадает с точной веткой и §851. */
  for (u = 0; u < n; u++) {
    int32_t p = ps[u];
    if (xi) break; /* контейнированных нет — весь сегмент точный */
    if (fc->pstamp[p] == fc->pkey) continue;
    if (!front_contained(fc->o->tribox + 6 * (int64_t)py->pcs[p].tri, blo, bhi)) continue;
    has_cont = 1;
    cntot += an[u];
    swt += wt[u];
    Lsurf += (front_le(fc, p) + front_rho(fc, p) * fc->Eprev[p] / (2.0 * M_PI)) * wt[u];
  }
  if (has_cont && swt > 0.0) {
    double T = 1.0 - cntot / csecnode, f;
    if (T < 0.0) T = 0.0;
    if (T > 1.0) T = 1.0;
    if (fc->tau0) T = 1.0; /* НК: слой не взаимодействует */
    f = 1.0 - T;
    if (f > 0.0) {
      Lsurf /= swt;
      for (u = 0; u < n; u++) {
        int32_t p = ps[u];
        if (xi) break;
        if (fc->pstamp[p] == fc->pkey) continue;
        if (!front_contained(fc->o->tribox + 6 * (int64_t)py->pcs[p].tri, blo, bhi)) continue;
        fc->Ed[p] += fc->w_d * Lin * f * csec * fc->axcos * wt[u] / swt / fc->area[p];
        sw_accum(fc, p,
                 fc->w_d * Lin * f * csec * fc->axcos * wt[u] / swt / fc->area[p]); /* §862 */
      }
      fc->depA += f * ai;
      fc->absorbed += fc->w_d * Lin * f * csec;
      fc->recycled += fc->w_d * f * (Lsurf - fc->le) * csec;
      fc->emitted += fc->w_d * fc->le * f * csec;
      *a = T * ai;
      *b = T * bi + f * Lsurf;
      if (fc->noprop) { /* НК: фронт не переносится */
        *a = 0.0;
        *b = 0.0;
      }
      goto done;
    }
  }
  /* прозрачная клетка: носитель не тронут */
  if (fc->noprop) {
    *a = 0.0;
    *b = 0.0;
  }
done:
}

/* А1576: точный МНОГОПОПАДНЫЙ проход сегмента (лист или узел): все
 * пересечения unstamped кусков с отрезком [tin,tout] по ходу трубки
 * (сортировка по t). У walk лестница ℓ_p плоская по построению — одна и та
 * же последовательность событий при любой сегментации (штамп А1564 гасит
 * дубли куска и повторы на стыках). Сключительно под ключом walk=1: модель
 * «только реальные пересечения» не совпадает с перехватной на стенках,
 * не выровненных по сетке (box: 0.88 против 4.17), — выбор за
 * фальсификатором А1576. */
static void front_seg_walk(front_ctx *fc, const int32_t *ps, int32_t n, double tin, double tout,
                           double *a, double *b) {
  const hz_pyr *py = fc->py;
  double csec = py->cell * py->cell; /* сечение БАЗОВОЙ клетки — инвариант А1563 */
  double ai = *a, bi = *b, Lin = ai + bi;
  double *ht;
  int32_t *hp;
  int32_t u;
  int64_t nh = 0, i;
  if (n <= 0) return;
  if (fc->tau0) return; /* НК: слой не взаимодействует */
  fc->nmat++;
  if (fc->noprop) { /* НК: фронт не переносится */
    *a = 0.0;
    *b = 0.0;
    return;
  }
  if (n > fc->abufcap) {
    int64_t nc = fc->abufcap ? fc->abufcap * 2 : 64;
    double *na, *nw;
    while (nc < n)
      nc *= 2;
    na = (double *)realloc(fc->abuf, (size_t)nc * sizeof *na);
    if (!na) return;
    fc->abuf = na;
    nw = (double *)realloc(fc->wbuf, (size_t)nc * sizeof *nw);
    if (!nw) return;
    fc->wbuf = nw;
    fc->abufcap = nc;
  }
  ht = fc->abuf;
  hp = (int32_t *)fc->wbuf; /* n int32 ≤ n double — места хватает */
  if (fc->agg && fc->agg_list) {
    /* §863/шаг 2: АГРЕГАЦИЯ — список узла отсортирован (mtl,id), группы
     * односоставных кусков идут подряд. ray-tri по ПРЕДСТАВИТЕЛЮ группы;
     * депозит и Lh — площадь-взвешенно по членам (диффузная доля отклика:
     * ответ группы определяется Σплощадью и средним материалом). hp хранит
     * НАЧАЛО группы; штамп — на представителе. */
    double sl = 1e-9 * (fabs(tin) + fabs(tout) + 1.0);
    double bh0, bh1;
    for (u = 0; u < n; u++) {
      int32_t gid = ps[u];
      int32_t rep = fc->ag->grep[gid];
      if (fc->pstamp[rep] != fc->pkey) {
        /* §864: префильтр bbox×сегмент — тот же консервативный отсекатель, что
         * в front_interact: треугольник ⊆ своего bbox, нет пересечения отрезка
         * с bbox — попадания в [tin,tout] не существует, результата не меняет */
        if (!front_box_seg(fc->o->tribox + 6 * (int64_t)py->pcs[rep].tri,
                           fc->o->tribox + 6 * (int64_t)py->pcs[rep].tri + 3, fc->org, fc->om, tin,
                           tout, &bh0, &bh1))
          continue;
        double tt = sw_ray_tri_raw(fc->org, fc->om, fc->o->trivert + 9 * (int64_t)py->pcs[rep].tri);
        if (tt >= tin - sl && tt <= tout + sl) {
          ht[nh] = tt;
          hp[nh] = (int32_t)u;
          nh++;
        }
      } else {
        fc->nstamp++;
      }
    }
  } else
    for (u = 0; u < n; u++) {
      int32_t p = ps[u];
      int32_t tri = py->pcs[p].tri;
      double tt, bh0, bh1;
      if (fc->pstamp[p] == fc->pkey) {
        fc->nstamp++; /* кратность 1 на (кусок, направление) — А1564 */
        continue;
      }
      /* §864: префильтр bbox×сегмент — консервативный отсекатель (как во
       * front_interact): треугольник ⊆ bbox, отрезок мимо bbox — мимо и
       * треугольника */
      if (!front_box_seg(fc->o->tribox + 6 * (int64_t)tri, fc->o->tribox + 6 * (int64_t)tri + 3,
                         fc->org, fc->om, tin, tout, &bh0, &bh1))
        continue;
      tt = sw_ray_tri_raw(fc->org, fc->om, fc->o->trivert + 9 * (int64_t)tri);
      /* Оба конца включены; от двойного события на стыке сегментов спасает
       * ШТАМП: кусок штампуется при первом же событии (проход депозита
       * ниже), повтор в следующем сегменте отбрасывается — событие
       * обрабатывается ровно один раз, поэтому лестница ℓ_p плоская. */
      {
        double sl = 1e-9 * (fabs(tin) + fabs(tout) + 1.0);
        if (tt >= tin - sl && tt <= tout + sl) {
          ht[nh] = tt;
          hp[nh] = p;
          nh++;
        }
      }
    }
  /* сортировка по ходу трубки (попаданий мало, вставками) */
  for (i = 1; i < nh; i++) {
    double t = ht[i];
    int32_t p = hp[i];
    int64_t v;
    for (v = i; v > 0 && ht[v - 1] > t; v--) {
      ht[v] = ht[v - 1];
      hp[v] = hp[v - 1];
    }
    ht[v] = t;
    hp[v] = p;
  }
  {
    int first = 1;
    for (i = 0; i < nh; i++) {
      int32_t p = hp[i];
      double Lh;
      if (fc->agg && fc->agg_list && fc->ag) { /* §863/шаг 2: O(1) на попадание */
        int32_t gid = ps[p];
        int32_t rep = fc->ag->grep[gid];
        double A = fc->ag->gA[gid], edep;
        if (fc->pstamp[rep] == fc->pkey) continue;
        fc->pstamp[rep] = fc->pkey;
        if (first) {
          fc->depA += ai;
          first = 0;
        }
        edep = fc->w_d * Lin * csec * fc->axcos / A;
        fc->ag->gacc[gid] += edep; /* раздача членам — после итерации */
        fc->ag->gwt[gid] += edep * (1.0 - front_cos(fc, rep) / (2.0 * fc->area[rep]));
        fc->absorbed += fc->w_d * Lin * csec;
        fc->emitted += fc->w_d * fc->ag->gleA[gid] * csec;
        Lh = fc->ag->gleA[gid] + fc->ag->grA[gid] * fc->ag->gSE[gid] / (A * 2.0 * M_PI);
        fc->recycled += fc->w_d * (Lh - fc->le) * csec;
        fc->ndep++;
        Lin = Lh;
        continue;
      }
      if (fc->pstamp[p] == fc->pkey) continue; /* дубль куска (список/стык) */
      /* ШТАМП при любом событии: пересечение обрабатывается один раз */
      fc->pstamp[p] = fc->pkey;
      if (first) { /* родительская доля трубки погашена первой поверхностью (G6) */
        fc->depA += ai;
        first = 0;
      }
      fc->Ed[p] += fc->w_d * Lin * csec * fc->axcos / fc->area[p];
      fc->absorbed += fc->w_d * Lin * csec;
      fc->emitted += fc->w_d * front_le(fc, p) * csec;
      Lh = front_le(fc, p) + front_rho(fc, p) * fc->Eprev[p] / (2.0 * M_PI);
      fc->recycled += fc->w_d * (Lh - fc->le) * csec;
      fc->ndep++;
      sw_accum(fc, p, fc->w_d * Lin * csec * fc->axcos / fc->area[p]); /* §862 */
      Lin = Lh;
    }
    if (nh > 0) {
      *a = 0.0;
      *b = Lin;
    }
  }
}

/* марш трубки через узел (уровень l, позиция pos; l<0 — лист) на отрезке
 * [tin,tout]; carry a/b сквозной — сечение трубки не делится (А1563) */
static void front_visit(front_ctx *fc, int32_t l, int32_t pos, double tin, double tout, double *a,
                        double *b) {
  const hz_pyr *py = fc->py;
  double blo[3], bhi[3];
  if (l < 0) {        /* лист: базовая клетка */
    fc->agg_list = 0; /* листовой список не отсортирован по материалу */
    if (fc->cfront) fc->ncellfront++;
    front_node_box(fc, -1, pos, blo, bhi);
    if (fc->o->walk) /* А1576: точный многопопадный проход */
      front_seg_walk(fc, fc->bpids + fc->bstart[pos], fc->bstart[pos + 1] - fc->bstart[pos], tin,
                     tout, a, b);
    else
      front_interact(fc, fc->bpids + fc->bstart[pos], fc->bstart[pos + 1] - fc->bstart[pos], blo,
                     bhi, tin, tout, a, b, 1);
    return;
  }
  /* МАТЕРИАЛЕН ⟺ max ℓ_p поддерева ≥ уровень узла (лист — уровень 0,
   * уровень l — parents of leaves = l+1): марш дробит фронт ДО уровня
   * ℓ_p куска; узлы глубже нужного — промежуточные (спуск), крупнее —
   * фронт уже поглощён на уровне ℓ_p и сюда не доходит. */
  uint8_t m = py->lev_lp[l][pos];
  if (m != 255 && (int32_t)m >= l + 1) {
    int32_t n = 0;
    const int32_t *ps = fc->pbuf;
    int64_t cid = -1;
    if (fc->o->walk && fc->nstart) cid = (int64_t)fc->ncluster0 + fc->noff[l] + pos;
    if (cid >= 0 && fc->nstart[cid] >= 0) {
      /* §863: список узла уже собран — gather не повторяем */
      ps = fc->nstore + fc->nstart[cid];
      n = fc->nlen[cid];
    } else {
      if (front_gather(fc, l, pos, &n) != 0) return;
      if (cid >= 0 && *fc->nstore_n + n <= fc->nstore_cap) {
        /* §863/шаг 2 (агрегация, этап 0): ДЕДУПЛИКАЦИЯ — кусок, чей bbox
         * накрывает несколько листов узла, попадал в список многократно;
         * штамп А1564 гасил дубли на депозите, но ray-tri по ним считался.
         * Порядок списка не важен: попадания сортируются по tt. */
        /* сортировка по (материал, id): группы агрегации идут подряд */
        for (int32_t q = 1; q < n; q++) {
          int32_t v = fc->pbuf[q], w2 = q;
          while (w2 > 0 &&
                 (py->pcs[fc->pbuf[w2 - 1]].mtl > py->pcs[v].mtl ||
                  (py->pcs[fc->pbuf[w2 - 1]].mtl == py->pcs[v].mtl && fc->pbuf[w2 - 1] > v))) {
            fc->pbuf[w2] = fc->pbuf[w2 - 1];
            w2--;
          }
          fc->pbuf[w2] = v;
        }
        int32_t un = 0;
        for (int32_t q = 0; q < n; q++)
          if (un == 0 || fc->pbuf[un - 1] != fc->pbuf[q]) fc->pbuf[un++] = fc->pbuf[q];
        n = un;
        int64_t st = *fc->nstore_n;
        if (fc->agg) { /* §863/шаг 2: список = ГИДЫ групп (по материалу) */
          int64_t g0 = 0;
          n = 0;
          int one = fc->o->agg >= 2; /* agg=2: одна группа на список */
          while (g0 < un) {
            int32_t mtl0 = py->pcs[fc->pbuf[g0]].mtl;
            int64_t g1 = g0;
            int32_t gid;
            while (g1 < un && (one || py->pcs[fc->pbuf[g1]].mtl == mtl0))
              g1++;
            gid = sw_group_make(fc, g0, g1);
            if (gid < 0) break; /* нет памяти: узел без кэша группы */
            fc->nstore[st + n] = gid;
            n++;
            g0 = g1;
          }
          *fc->nstore_n = st + n;
          fc->nstart[cid] = st;
          fc->nlen[cid] = n;
          ps = fc->nstore + st;
        } else {
          for (int32_t q = 0; q < n; q++)
            fc->nstore[st + q] = fc->pbuf[q];
        }
        *fc->nstore_n = st + n;
        fc->nstart[cid] = st;
        fc->nlen[cid] = n;
        ps = fc->nstore + st;
      } /* магазин полон: работаем по pbuf без кэша (прежний путь) */
      fc->agg_list = fc->agg; /* узловой список отсортирован (mtl,id) */
    }
    front_node_box(fc, l, pos, blo, bhi);
    /* А1573: на узловом отрезке клеточный T-перехват НЕзаконен — bbox куска
     * ⊆ узел не означает, что трубка задевает его тень; перехват раздаёт
     * сечение узла кускам, мимо которых трубка прошла (расходимость
     * ~×1.5/итерацию на cavity05 lp=2). Узел — только точное пересечение.
     * Предикат А1566 «bbox ⊆ клетка» валиден ТОЛЬКО на листовом отрезке
     * (bbox ⊆ БАЗОВАЯ клетка). */
    if (fc->o->walk) /* А1576: точный многопопадный проход */
      front_seg_walk(fc, ps, n, tin, tout, a, b);
    else
      front_interact(fc, fc->pbuf, n, blo, bhi, tin, tout, a, b, 0);
    return;
  }
  /* ПРОМЕЖУТОЧНЫЙ: спуск — продолжение DDA на детском уровне */
  fc->ndesc++;
  if (l == 0) {
    /* спуск на ЛИСТЬЯ: узел уровня-0 покрывает 2³ БАЗОВЫХ клетки —
     * перечисляются сами базовые клетки (id уровня-0 ≠ id базы;
     * без этого на нечётных сетках терялся весь хвост трубки) */
    int64_t pd0[3], kid0 = py->lev[0][pos].id, ci0[3];
    double c0[8], c1[8];
    int32_t cpos[8];
    int n = 0, u, v, i;
    int cx, cy, cz;
    double a_in = *a, dep_snap = fc->depA;
    hz_pyr_node *pnode = &py->lev[0][pos]; /* §864/Б1: Existence-маска детей */
    front_ldims(fc, 0, pd0);
    ci0[0] = kid0 % pd0[0];
    ci0[1] = (kid0 / pd0[0]) % pd0[1];
    ci0[2] = kid0 / (pd0[0] * pd0[1]);
    for (cz = 0; cz < 2; cz++)
      for (cy = 0; cy < 2; cy++)
        for (cx = 0; cx < 2; cx++) {
          int64_t bx = 2 * ci0[0] + cx, by = 2 * ci0[1] + cy, bz = 2 * ci0[2] + cz;
          double cblo[3], cbhi[3], h0, h1;
          int32_t found;
          int64_t id;
          if (bx >= py->nx || by >= py->ny || bz >= py->nz) continue;
          id = bx + py->nx * (by + py->ny * bz);
          /* §864/Б1: существование листа — бит-тест ДО box_seg: пустая клетка
           * обходится без слэб-теста и двоичного поиска (А1569). Маска узла
           * построена по тем же листьям — ответ совпадает с hz_pyr_leaf_pos. */
          if (!(pnode->chmask & (1u << (cx | (cy << 1) | (cz << 2))))) {
            fc->njump++; /* ПУСТ: прыжок одним шагом (А1569) */
            continue;
          }
          cblo[0] = py->lo[0] + (double)bx * py->cell;
          cblo[1] = py->lo[1] + (double)by * py->cell;
          cblo[2] = py->lo[2] + (double)bz * py->cell;
          cbhi[0] = cblo[0] + py->cell;
          cbhi[1] = cblo[1] + py->cell;
          cbhi[2] = cblo[2] + py->cell;
          for (u = 0; u < 3; u++)
            if (cbhi[u] > fc->hi[u]) cbhi[u] = fc->hi[u];
          if (!front_box_seg(cblo, cbhi, fc->org, fc->om, tin, tout, &h0, &h1)) continue;
          found = hz_pyr_leaf_pos(py, id);
          if (found < 0) {
            fc->njump++; /* ПУСТ: прыжок одним шагом по tmax клетки (А1569) */
            continue;
          }
          c0[n] = h0;
          c1[n] = h1;
          cpos[n] = found;
          n++;
        }
    for (u = 1; u < n; u++) {
      double s = c0[u], s1 = c1[u];
      int32_t pu = cpos[u];
      /* А1578: при равных c0 вперёд идёт вырожденный (h0==h1) сегмент —
       * клетка за плоской стенкой; иначе её подхват выполнялся после
       * всего марша и депонировал Lh на саму входную стенку */
      for (v = u; v > 0 && (c0[v - 1] > s || (fc->o->walk && !(c0[v - 1] < s) && !(s < c0[v - 1]) &&
                                              c1[v - 1] > s1));
           v--) {
        c0[v] = c0[v - 1];
        c1[v] = c1[v - 1];
        cpos[v] = cpos[v - 1];
      }
      c0[v] = s;
      c1[v] = s1;
      cpos[v] = pu;
    }
    for (i = 0; i < n; i++)
      front_visit(fc, -1, cpos[i], c0[i], c1[i], a, b);
    if (!fc->noprop && fabs(a_in - (fc->depA - dep_snap) - *a) > 1e-9 * (1.0 + fabs(a_in)))
      fc->g6viol++;
    return;
  }
  {
    int64_t pd[3], cd[3], kid = py->lev[l][pos].id, ci[3];
    int32_t cl = l - 1;
    double side2 = py->cell * (double)((int64_t)1 << l); /* сторона ребёнка */
    double c0[8], c1[8];
    int32_t cpos[8];
    int n = 0, u, v, i;
    int cx, cy, cz;
    double a_in = *a, dep_snap = fc->depA;
    hz_pyr_node *pnode = &py->lev[l][pos]; /* §864/Б1: Existence-маска детей */
    front_ldims(fc, l, pd);
    if (cl < 0) {
      cd[0] = py->nx;
      cd[1] = py->ny;
      cd[2] = py->nz;
    } else
      hz_pyr_level_dims(py, cl, cd);
    ci[0] = kid % pd[0];
    ci[1] = (kid / pd[0]) % pd[1];
    ci[2] = kid / (pd[0] * pd[1]);
    for (cz = 0; cz < 2; cz++)
      for (cy = 0; cy < 2; cy++)
        for (cx = 0; cx < 2; cx++) {
          int64_t qx = 2 * ci[0] + cx, qy = 2 * ci[1] + cy, qz = 2 * ci[2] + cz;
          double cblo[3], cbhi[3], h0, h1;
          int32_t found;
          int64_t id;
          if (qx >= cd[0] || qy >= cd[1] || qz >= cd[2]) continue;
          id = qx + cd[0] * (qy + cd[1] * qz);
          /* §864/Б1: существование узла-ребёнка — бит-тест ДО box_seg (А1569):
           * маска узла построена подъёмом при построении пирамиды и совпадает
           * с ответом hz_pyr_node_pos. Пустой ребёнок — один бит-тест. */
          if (!(pnode->chmask & (1u << (cx | (cy << 1) | (cz << 2))))) {
            fc->njump++; /* ПУСТ: прыжок одним шагом по tmax ребёнка (А1569) */
            continue;
          }
          cblo[0] = py->lo[0] + (double)qx * side2;
          cblo[1] = py->lo[1] + (double)qy * side2;
          cblo[2] = py->lo[2] + (double)qz * side2;
          cbhi[0] = cblo[0] + side2;
          cbhi[1] = cblo[1] + side2;
          cbhi[2] = cblo[2] + side2;
          for (u = 0; u < 3; u++)
            if (cbhi[u] > fc->hi[u]) cbhi[u] = fc->hi[u];
          if (!front_box_seg(cblo, cbhi, fc->org, fc->om, tin, tout, &h0, &h1)) continue;
          if (cl < 0)
            found = hz_pyr_leaf_pos(py, id);
          else
            found = hz_pyr_node_pos(py, cl, id);
          if (found < 0) {
            fc->njump++; /* ПУСТ: прыжок одним шагом по tmax ребёнка (А1569) */
            continue;
          }
          c0[n] = h0;
          c1[n] = h1;
          cpos[n] = found;
          n++;
        }
    /* порядок по ходу луча (n ≤ 8, вставками) */
    for (u = 1; u < n; u++) {
      double s = c0[u], s1 = c1[u];
      int32_t pu = cpos[u];
      /* А1578: при равных c0 вперёд идёт вырожденный (h0==h1) сегмент —
       * клетка за плоской стенкой; иначе её подхват выполнялся после
       * всего марша и депонировал Lh на саму входную стенку */
      for (v = u; v > 0 && (c0[v - 1] > s || (fc->o->walk && !(c0[v - 1] < s) && !(s < c0[v - 1]) &&
                                              c1[v - 1] > s1));
           v--) {
        c0[v] = c0[v - 1];
        c1[v] = c1[v - 1];
        cpos[v] = cpos[v - 1];
      }
      c0[v] = s;
      c1[v] = s1;
      cpos[v] = pu;
    }
    for (i = 0; i < n; i++) {
      front_visit(fc, cl, cpos[i], c0[i], c1[i], a, b);
    }
    /* G6: родительская доля сохраняется: a_in = Σdep_a(дети) + a_out */
    if (!fc->noprop && fabs(a_in - (fc->depA - dep_snap) - *a) > 1e-9 * (1.0 + fabs(a_in)))
      fc->g6viol++;
  }
}

/* §865/А1580 (path=1): попадания клетки в список полного пути трубки */
static void path_collect(front_ctx *fc, int32_t l, int32_t pos, double tin, double tout) {
  const hz_pyr *py = fc->py;
  const int32_t *ps;
  int32_t n = 0, u;
  double blo[3], bhi[3];
  if (l < 0) {
    ps = fc->bpids + fc->bstart[pos];
    n = fc->bstart[pos + 1] - fc->bstart[pos];
  } else {
    int64_t cid;
    front_node_box(fc, l, pos, blo, bhi);
    if (front_gather(fc, l, pos, &n) != 0) return;
    cid = (int64_t)fc->ncluster0 + fc->noff[l] + pos;
    if (cid >= 0 && fc->nstart[cid] >= 0) {
      ps = fc->nstore + fc->nstart[cid];
      n = fc->nlen[cid];
    } else
      ps = fc->pbuf;
  }
  for (u = 0; u < n; u++) {
    int32_t p = ps[u];
    double tt = sw_ray_tri_raw(fc->org, fc->om, fc->o->trivert + 9 * (int64_t)py->pcs[p].tri);
    if (!(tt >= tin) || !(tt <= tout)) continue;
    if (fc->pstamp[p] == fc->pkey) continue;
    if (fc->pbuf_n == fc->pbuf_cap) {
      int64_t nc = fc->pbuf_cap ? fc->pbuf_cap * 2 : 256;
      double *nt = (double *)realloc(fc->pbuf_t, (size_t)nc * sizeof *nt);
      int32_t *np = (int32_t *)realloc(fc->pbuf_p, (size_t)nc * sizeof *np);
      if (!nt || !np) return;
      fc->pbuf_t = nt;
      fc->pbuf_p = np;
      fc->pbuf_cap = nc;
    }
    fc->pstamp[p] = fc->pkey;
    fc->pbuf_t[fc->pbuf_n] = tt;
    fc->pbuf_p[fc->pbuf_n] = p;
    fc->pbuf_n++;
  }
}

/* §865/А1580 (path=1): DDA уровня jt (jt<0 — листья) внутри [tin,tout] */
static void path_dda(front_ctx *fc, int32_t l, double tin, double tout, const double blo[3],
                     const double bhi[3]) {
  const hz_pyr *py = fc->py;
  const double *om = fc->om;
  double side, tmax[3], tdelta[3], te = tin;
  int64_t g[3], dim[3], stepmax, s, pitch1, pitch2;
  int ax;
  if (l < 0) {
    side = py->cell;
    dim[0] = py->nx;
    dim[1] = py->ny;
    dim[2] = py->nz;
  } else {
    side = py->cell * (double)((int64_t)1 << (l + 1));
    front_ldims(fc, l, dim);
  }
  for (ax = 0; ax < 3; ax++) {
    double p = fc->org[ax] + om[ax] * tin;
    int64_t gg;
    if (p < blo[ax]) p = blo[ax];
    if (p > bhi[ax]) p = bhi[ax];
    gg = (int64_t)floor((p - py->lo[ax]) / side);
    if (gg < 0) gg = 0;
    if (gg >= dim[ax]) gg = dim[ax] - 1;
    g[ax] = gg;
    if (fabs(om[ax]) < 1e-30) {
      tmax[ax] = 1e30;
      tdelta[ax] = 0.0;
    } else {
      double plane = py->lo[ax] + (double)(g[ax] + (om[ax] > 0.0 ? 1 : 0)) * side;
      tmax[ax] = (plane - fc->org[ax]) / om[ax];
      if (tmax[ax] < tin) tmax[ax] = tin;
      tdelta[ax] = side / fabs(om[ax]);
    }
  }
  pitch1 = l < 0 ? py->nx : dim[0];
  pitch2 = l < 0 ? py->ny : dim[1]; /* id = g0 + pitch1*(g1 + pitch2*g2) */
  stepmax = dim[0] + dim[1] + dim[2] + 8;
  for (s = 0; s < stepmax; s++) {
    double tn = tmax[0];
    int axm;
    int64_t id;
    int32_t pos;
    if (tmax[1] < tn) tn = tmax[1];
    if (tmax[2] < tn) tn = tmax[2];
    axm = (tmax[0] <= tmax[1] && tmax[0] <= tmax[2]) ? 0 : (tmax[1] <= tmax[2] ? 1 : 2);
    if (tn > tout) tn = tout;
    if (g[0] >= 0 && g[1] >= 0 && g[2] >= 0 && g[0] < dim[0] && g[1] < dim[1] && g[2] < dim[2]) {
      id = g[0] + pitch1 * (g[1] + pitch2 * g[2]);
      pos = l < 0 ? hz_pyr_leaf_pos(py, id) : hz_pyr_node_pos(py, l, id);
      if (getenv("HZ_PATH") && fc->ntube == 1500 && (s < 4 || s > 60))
        fprintf(stderr, "STEP s=%lld g=(%lld,%lld,%lld) te=%.5g tn=%.5g leaf=%d\n", (long long)s,
                (long long)g[0], (long long)g[1], (long long)g[2], te, tn, pos);
      if (pos >= 0)
        path_collect(fc, l, pos, te, tn);
      else
        fc->njump++;
    }
    if (tn >= tout) break;
    te = tn;
    if (axm == 0) {
      tmax[0] += tdelta[0];
      g[0] += om[0] > 0.0 ? 1 : -1;
    } else if (axm == 1) {
      tmax[1] += tdelta[1];
      g[1] += om[1] > 0.0 ? 1 : -1;
    } else {
      tmax[2] += tdelta[2];
      g[2] += om[2] > 0.0 ? 1 : -1;
    }
  }
}

/* одна трубка: прибор «до» (полный базовый DDA), затем иерархический марш */
static void front_tube(front_ctx *fc, const int64_t cc[3], double *lostA, double *lostB) {
  const hz_pyr *py = fc->py;
  double h0, h1, a = 0.0, b = 0.0;
  int q;
  for (q = 0; q < 3; q++)
    fc->org[q] = py->lo[q] + ((double)cc[q] + 0.5) * py->cell;
  if (!front_box_seg(py->lo, fc->hi, fc->org, fc->om, -1e30, 1e30, &h0, &h1)) return;
  /* прибор «до»: число базовых клеток полным DDA (стоимость §851 vc=1) */
  {
    double tmax[3], tdelta[3], sent = h0;
    int64_t cell[3], stepv[3];
    for (q = 0; q < 3; q++) {
      double pos = fc->org[q];
      if (fabs(fc->om[q]) < 1e-30) {
        tmax[q] = 1e30;
        tdelta[q] = 0.0;
        stepv[q] = 0;
      } else {
        double boundary = fc->om[q] > 0.0 ? py->lo[q] + ((double)cc[q] + 1.0) * py->cell
                                          : py->lo[q] + (double)cc[q] * py->cell;
        stepv[q] = fc->om[q] > 0.0 ? 1 : -1;
        tdelta[q] = py->cell / fabs(fc->om[q]);
        tmax[q] = (boundary - pos) / fc->om[q];
        if (tmax[q] < 0.0) tmax[q] = 0.0;
      }
      cell[q] = cc[q];
    }
    for (;;) {
      double tn = tmax[0];
      if (tmax[1] < tn) tn = tmax[1];
      if (tmax[2] < tn) tn = tmax[2];
      if (fc->cbase) fc->ncellbase++;
      if (tn > 1e29 || sent >= h1) break;
      sent = tn;
      if (tmax[0] <= tmax[1] && tmax[0] <= tmax[2]) {
        tmax[0] += tdelta[0];
        cell[0] += stepv[0];
      } else if (tmax[1] <= tmax[2]) {
        tmax[1] += tdelta[1];
        cell[1] += stepv[1];
      } else {
        tmax[2] += tdelta[2];
        cell[2] += stepv[2];
      }
      if (cell[0] < 0 || cell[1] < 0 || cell[2] < 0 || cell[0] >= py->nx || cell[1] >= py->ny ||
          cell[2] >= py->nz)
        break;
    }
  }
  if (py->nlev == 0) { /* пирамиды нет: базовый DDA, взаимодействие на листах */
    double tmax[3], tdelta[3], sent = h0;
    int64_t cell[3], stepv[3];
    for (q = 0; q < 3; q++) {
      double pos = fc->org[q];
      if (fabs(fc->om[q]) < 1e-30) {
        tmax[q] = 1e30;
        tdelta[q] = 0.0;
        stepv[q] = 0;
      } else {
        double boundary = fc->om[q] > 0.0 ? py->lo[q] + ((double)cc[q] + 1.0) * py->cell
                                          : py->lo[q] + (double)cc[q] * py->cell;
        stepv[q] = fc->om[q] > 0.0 ? 1 : -1;
        tdelta[q] = py->cell / fabs(fc->om[q]);
        tmax[q] = (boundary - pos) / fc->om[q];
        if (tmax[q] < 0.0) tmax[q] = 0.0;
      }
      cell[q] = cc[q];
    }
    for (;;) {
      int32_t pos;
      double tn = tmax[0], te = sent;
      if (tmax[1] < tn) tn = tmax[1];
      if (tmax[2] < tn) tn = tmax[2];
      if (tn > 1e29) tn = h1;
      pos = hz_pyr_leaf_pos(py, cell[0] + py->nx * (cell[1] + py->ny * cell[2]));
      if (pos >= 0 && tn > te) front_visit(fc, -1, pos, te, tn, &a, &b);
      if (tn >= h1) break;
      sent = tn;
      if (tmax[0] <= tmax[1] && tmax[0] <= tmax[2]) {
        tmax[0] += tdelta[0];
        cell[0] += stepv[0];
      } else if (tmax[1] <= tmax[2]) {
        tmax[1] += tdelta[1];
        cell[1] += stepv[1];
      } else {
        tmax[2] += tdelta[2];
        cell[2] += stepv[2];
      }
      if (cell[0] < 0 || cell[1] < 0 || cell[2] < 0 || cell[0] >= py->nx || cell[1] >= py->ny ||
          cell[2] >= py->nz)
        break;
    }
  } else if (fc->o->path) { /* §865/А1580: марш полного пути */
    uint8_t m = py->lev_lp[py->nlev - 1][0];
    int32_t jt = (m == 255 || m == 0) ? -1 : (int32_t)m - 1;
    double csec = py->cell * py->cell;
    int64_t i;
    double Lin = 0.0;
    fc->nmat++;
    fc->pbuf_n = 0;
    path_dda(fc, jt, h0, h1, py->lo, fc->hi);
    for (i = 1; i < fc->pbuf_n; i++) {
      double tv = fc->pbuf_t[i];
      int32_t pv = fc->pbuf_p[i];
      int64_t v = i;
      while (v > 0 && fc->pbuf_t[v - 1] > tv) {
        fc->pbuf_t[v] = fc->pbuf_t[v - 1];
        fc->pbuf_p[v] = fc->pbuf_p[v - 1];
        v--;
      }
      fc->pbuf_t[v] = tv;
      fc->pbuf_p[v] = pv;
    }
    for (i = 0; i < fc->pbuf_n; i++) {
      int32_t p = fc->pbuf_p[i];
      double Lh;
      fc->Ed[p] += fc->w_d * Lin * csec * fc->axcos / fc->area[p];
      fc->absorbed += fc->w_d * Lin * csec;
      fc->emitted += fc->w_d * front_le(fc, p) * csec;
      Lh = front_le(fc, p) + front_rho(fc, p) * fc->Eprev[p] / (2.0 * M_PI);
      fc->recycled += fc->w_d * (Lh - fc->le) * csec;
      fc->ndep++;
      sw_accum(fc, p, fc->w_d * Lin * csec * fc->axcos / fc->area[p]);
      Lin = Lh;
    }
    a = 0.0;
    b = Lin;
  } else {
    front_visit(fc, py->nlev - 1, 0, h0, h1, &a, &b);
    /* G6 на корне: вошедшая (нулевая — тёмный вход) доля = Σdep + вышедшая */
    if (!fc->noprop && fabs(0.0 - fc->depA - a) > 1e-9 * (1.0 + fabs(a))) fc->g6viol++;
  }
  if (a > 0.0 || b > 0.0)
    fc->nlostseg++; /* остаток ушёл за границу домена — в lost, не исчез (А1567) */
  *lostA += a;
  *lostB += b;
}

/* направление целиком: трубки = линии базовой сетки (дедупликация штампом
 * клетки, как §851) */
static void front_dir(front_ctx *fc, int32_t *stampv) {
  const hz_pyr *py = fc->py;
  int32_t ax;
  int32_t f[3];
  for (ax = 0; ax < 3; ax++) {
    int64_t n = ax == 0 ? py->nx : (ax == 1 ? py->ny : py->nz);
    if (fabs(fc->om[ax]) < 1e-30) {
      f[ax] = -1;
      continue;
    }
    f[ax] = fc->om[ax] > 0.0 ? 0 : (int32_t)(n - 1);
  }
  for (ax = 0; ax < 3; ax++) {
    int64_t b1, b2;
    int32_t bxA = (ax + 1) % 3, bxB = (ax + 2) % 3;
    int64_t nb1, nb2;
    if (f[ax] < 0) continue;
    fc->ax = ax; /* семейство трубок: вход через грань оси ax */
    nb1 = bxA == 0 ? py->nx : (bxA == 1 ? py->ny : py->nz);
    nb2 = bxB == 0 ? py->nx : (bxB == 1 ? py->ny : py->nz);
    for (b1 = 0; b1 < nb1; b1++)
      for (b2 = 0; b2 < nb2; b2++) {
        int64_t cc[3], id;
        double lostA = 0.0, lostB = 0.0;
        cc[ax] = f[ax];
        cc[bxA] = b1;
        cc[bxB] = b2;
        id = cc[0] + py->nx * (cc[1] + py->ny * cc[2]);
        if (stampv[id] == fc->mark) continue; /* трубка уже шла (другая семья) */
        fc->axcos = fabs(fc->om[fc->ax]);
        stampv[id] = fc->mark;
        fc->pkey = ((int64_t)fc->mark << 32) | (uint32_t)fc->ntube++;
        fc->depA = 0.0; /* G6 — по РОДИТЕЛЬСКОЙ доле ЭТОЙ трубки */
        front_tube(fc, cc, &lostA, &lostB);
        if (lostA > 0.0 || lostB > 0.0) fc->lost += (lostA + lostB) * py->cell * py->cell;
      }
  }
}

/* А1578: индекс клетки по оси — ТА ЖЕ конвенция, что pyr_axis_index в pyr.c
 * (floor с клэмпом к сетке): клетки, помеченные pyr, гарантированно листья, и
 * марш их достигает (фантомную колонку у max-грани домена — вырожденным
 * сегментом [h,h], §852). Прежний span «floor…ceil−1» с клэмпом к домену клал
 * кусок плоской max-грани в слой, где листа чаще нет, — кусок выпадал из
 * bbox-индекса, и стенка теряла депозиты (А1578: 0.22 события/трубку). */
/* §863/шаг 2: создать группу из pbuf[g0..g1) (отсортирован по (mtl,id));
 * возвращает gid или -1 (нет памяти). */
static int32_t sw_group_make(front_ctx *fc, int64_t g0, int64_t g1) {
  sw_agg *ag = fc->ag;
  if (ag->ng == ag->ngcap) {
    int64_t nc = ag->ngcap ? ag->ngcap * 2 : 1024;
    int32_t *np, *nc2;
    double *na, *nl, *nr, *ns, *nc3, *nw;
    int64_t *no;
    np = (int32_t *)realloc(ag->grep, (size_t)nc * sizeof *np);
    na = (double *)realloc(ag->gA, (size_t)nc * sizeof *na);
    nl = (double *)realloc(ag->gleA, (size_t)nc * sizeof *nl);
    nr = (double *)realloc(ag->grA, (size_t)nc * sizeof *nr);
    ns = (double *)realloc(ag->gSE, (size_t)nc * sizeof *ns);
    nc3 = (double *)realloc(ag->gacc, (size_t)nc * sizeof *nc3);
    nw = (double *)realloc(ag->gwt, (size_t)nc * sizeof *nw);
    no = (int64_t *)realloc(ag->ggoff, (size_t)nc * sizeof *no);
    nc2 = (int32_t *)realloc(ag->ggcnt, (size_t)nc * sizeof *nc2);
    if (!np || !na || !nl || !nr || !ns || !nc3 || !nw || !no || !nc2) return -1;
    ag->grep = np;
    ag->gA = na;
    ag->gleA = nl;
    ag->grA = nr;
    ag->gSE = ns;
    ag->gacc = nc3;
    ag->gwt = nw;
    ag->ggoff = no;
    ag->ggcnt = nc2;
    ag->ngcap = nc;
  }
  {
    int32_t gid = (int32_t)ag->ng;
    int32_t rep = fc->pbuf[g0];
    double A = 0, leA = 0, rA = 0;
    int64_t q;
    if (ag->gmem_n + (g1 - g0) > ag->gmem_cap) {
      int64_t nc = ag->gmem_cap ? ag->gmem_cap * 2 : 4096;
      int32_t *nm;
      while (nc < ag->gmem_n + (g1 - g0))
        nc *= 2;
      nm = (int32_t *)realloc(ag->gmem, (size_t)nc * sizeof *nm);
      if (!nm) return -1;
      ag->gmem = nm;
      ag->gmem_cap = nc;
    }
    ag->ggoff[gid] = ag->gmem_n;
    for (q = g0; q < g1; q++) {
      int32_t m = fc->pbuf[q];
      double am = fc->area[m];
      ag->gmem[ag->gmem_n++] = m;
      A += am;
      leA += front_le(fc, m) * am;
      rA += front_rho(fc, m) * am;
    }
    ag->ggcnt[gid] = (int32_t)(g1 - g0);
    ag->grep[gid] = rep;
    ag->gA[gid] = A;
    ag->gleA[gid] = leA / A;
    ag->grA[gid] = rA / A;
    ag->gSE[gid] = 0.0; /* пересчитывается на итерацию */
    ag->gacc[gid] = 0.0;
    ag->gwt[gid] = 0.0;
    ag->ng++;
    return gid;
  }
}

static void sw_index_emit(int32_t pos, int32_t s, int32_t *bstart, int32_t *bpids, int32_t *fillb) {
  if (bpids)
    bpids[fillb[pos]++] = s;
  else
    bstart[pos + 1]++;
}

static void sw_index_piece(const hz_pyr *py, const hz_sw_opts *o, int32_t s, int fb,
                           int32_t *bstart, int32_t *bpids, int32_t *fillb) {
  int32_t tri = py->pcs[s].tri;
  int64_t lim[3], lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
  lim[0] = py->nx;
  lim[1] = py->ny;
  lim[2] = py->nz;
  for (int ax = 0; ax < 3; ax++) {
    double cmin = o->tribox[6 * (int64_t)tri + ax];
    double cmax = o->tribox[6 * (int64_t)tri + 3 + ax];
    double dom = o->domhi[ax] - py->lo[ax];
    lo[ax] = (int64_t)((cmin - py->lo[ax]) / py->cell);
    hi[ax] = (int64_t)ceil((cmax - py->lo[ax]) / py->cell) - 1;
    if (lo[ax] < 0) lo[ax] = 0;
    if (hi[ax] < lo[ax]) hi[ax] = lo[ax];
    if (hi[ax] > (int64_t)ceil(dom / py->cell) - 1) hi[ax] = (int64_t)ceil(dom / py->cell) - 1;
    if (lo[ax] > (int64_t)ceil(dom / py->cell) - 1) lo[ax] = (int64_t)ceil(dom / py->cell) - 1;
  }
  for (int64_t iz = lo[2]; iz <= hi[2]; iz++)
    for (int64_t iy = lo[1]; iy <= hi[1]; iy++)
      for (int64_t ix = lo[0]; ix <= hi[0]; ix++) {
        int32_t pos = hz_pyr_leaf_pos(py, ix + py->nx * (iy + py->ny * iz));
        if (pos >= 0) {
          sw_index_emit(pos, s, bstart, bpids, fillb);
          continue;
        }
        if (!fb) continue;
        for (int ax = 0; ax < 3; ax++) {
          int64_t nb[3];
          nb[0] = ix;
          nb[1] = iy;
          nb[2] = iz;
          for (int sgn = -1; sgn <= 1; sgn += 2) {
            int64_t mx = nb[ax];
            nb[ax] = mx + sgn;
            if (nb[ax] >= 0 && nb[ax] < lim[ax]) {
              int32_t npos = hz_pyr_leaf_pos(py, nb[0] + py->nx * (nb[1] + py->ny * nb[2]));
              if (npos >= 0) sw_index_emit(npos, s, bstart, bpids, fillb);
            }
            nb[ax] = mx;
          }
        }
      }
}

int hz_sw_run(hz_pyr *py, int32_t nt, const double *area, const double *nrm, const double *kd,
              const hz_sw_opts *o, hz_sw_stat *st, double *e_hist) {
  double *Ed = NULL, *Eprev = NULL;
  int32_t *order = NULL, *vindex = NULL;
  uint64_t *bits = NULL;
  hz_sw_walk *walks = NULL;
  hz_sw_dir *tab = NULL;
  int32_t *stampv = NULL; /* §851: штампы визитов линий [ncells] */
  int64_t *pstamp = NULL; /* §852/А1564: штамп (трубка × кусок) [nt] */
  int32_t *bstart = NULL, *bpids = NULL, *fillb = NULL; /* §852/А1566: bbox-индекс */
  /* §863: кэш списков кусков узлов (владелец — прогон, переживает fc) */
  int64_t *nstart = NULL;
  int32_t *nlen = NULL, *nstore = NULL;
  int64_t nstore_n = 0, nstore_cap = 0, ncluster = 0;
  int64_t noff[257] = {0};
  int64_t (*fld)[3] = NULL; /* §864/Б1: кэш размеров уровней (257 уровней max) */
  uint8_t *fldok = NULL;
  sw_agg ag; /* §863/шаг 2: таблица групп агрегации */
  memset(&ag, 0, sizeof ag);
  int nd = 0, d, it, rc = 0;
  double csec;
  uint8_t *lpflo = NULL; /* §862: этажи (когда задан lpacc) */

  memset(st, 0, sizeof *st);
  if (!py || !area || !nrm || !kd || !o) return 1;
  if (o->iters <= 0) return 1;
  if (!py->pcs || nt != py->nt) return 1;
  if (o->mode == 2 && o->vc && !o->trivert) return 1; /* §852: нужно точное пересечение */
  if (o->mode == 3 && (!o->trivert || !o->tribox || !o->lp || !o->domhi))
    return 1;                                 /* §852: фронт */
  if (o->mode == 3 && !py->leaf_lp) return 1; /* §852/А1567: per-node max ℓ_p не задан */
  /* §843: режим 0 (скалярный марш §835) определён только на легаси-наборах */
  if (o->mode == 0 && o->ndirs != 6 && o->ndirs != 26) return 1;
  rc = hz_sw_dir_table(o->ndirs, &tab, &nd);
  if (rc != 0) goto done;
  csec = py->cell * py->cell;

  Ed = (double *)calloc((size_t)nt, sizeof *Ed);
  Eprev = (double *)calloc((size_t)nt, sizeof *Eprev);
  order = (int32_t *)calloc(
      (size_t)py->nleaf,
      sizeof *order); /* calloc: анализатор видит инициализацию (FP-класс diam 07-24) */
  vindex = (int32_t *)malloc((size_t)py->nleaf * sizeof *vindex);
  bits = (uint64_t *)malloc((size_t)((py->nleaf + 63) >> 6) * sizeof *bits);
  walks = (hz_sw_walk *)calloc((size_t)nd, sizeof *walks);
  if (o->mode == 3) { /* §864/Б1: кэш размеров уровней для марша фронта */
    fld = (int64_t (*)[3])calloc(257, sizeof *fld);
    fldok = (uint8_t *)calloc(257, sizeof *fldok);
  }
  if (!Ed || !Eprev || !order || !vindex || !bits || !walks || (o->mode == 3 && (!fld || !fldok))) {
    rc = 2;
    goto done;
  }

  if (o->mode == 1 || (o->mode == 2 && !o->vc)) {
    sw_inv_acc = 0;
    for (d = 0; d < nd; d++) {
      int sg[3];
      for (int ax = 0; ax < 3; ax++)
        sg[ax] = (tab[d].om[ax] > 0.0) - (tab[d].om[ax] < 0.0);
      if (o->mode == 2) {
        if (o->build == 1)
          rc = sw_build_walk_sort(py, tab[d].om, area, nrm, kd, &walks[d]);
        else
          rc = sw_build_walk_march(py, tab[d].om, area, nrm, kd, vindex, bits, &walks[d]);
      } else
        rc = sw_build_walk(py, sg, tab[d].om, area, nrm, kd, &walks[d]);
      if (rc != 0) goto done;
      /* §845/А1525: хеш последовательности листьев — прибор против no-op:
       * порядок mode 2 обязан отличаться от mode 1, иначе G1 ничего не меряет.
       * Личность листа — его первый кусок (кусок сидит ровно в одном листе). */
      for (int32_t k = 0; k < walks[d].n; k++)
        st->order_hash = st->order_hash * 1099511628211ULL ^
                         (uint64_t)(uint32_t)walks[d].wpid[walks[d].leaf[k].first];
    }
    st->ninv = sw_inv_acc;
  } else {
    /* scalar §835: марш по всем листьям, 6 осей */
    {
      int32_t li2;
      for (li2 = 0; li2 < py->nleaf; li2++)
        order[li2] = li2; /* анализатор: полная инициализация */
    }
    hz_pyr_march_stat ms;
    double om[3] = {1.0, 0.0, 0.0};
    int32_t li;
    hz_pyr_march(py, om, 0, &ms, bits, vindex);
    for (li = 0; li < py->nleaf; li++)
      order[vindex[li] - 1] = li;
  }

  /* §851: ЛУЧЕВОЙ СВИП (vc=1): трубки = линии сетки, см. sw_line_sweep_dir.
   * Латеральная изоляция и «стороны» возникают сами — разметка компонент
   * не нужна (упрощение против первой редакции плана). Штампы исключают
   * двойные визиты клетки. */
  if ((o->mode == 2 && o->vc) || o->mode == 3) {
    int64_t ncells = (int64_t)py->nx * py->ny * py->nz;
    stampv = (int32_t *)calloc((size_t)ncells, sizeof *stampv);
    if (!stampv) {
      rc = 2;
      goto done;
    }
  }
  if (o->mode == 3) { /* §852/А1564: штамп кратности 1 на (кусок, направление) */
    pstamp = (int64_t *)calloc((size_t)nt, sizeof *pstamp);
    bstart = (int32_t *)calloc((size_t)py->nleaf + 1, sizeof *bstart);
    if (o->walk) { /* §863: кэш списков кусков узлов */
      ncluster = (int64_t)py->nleaf;
      for (int32_t l2 = 0; l2 < py->nlev; l2++) {
        noff[l2] = ncluster - (int64_t)py->nleaf;
        ncluster += py->nlev_nodes[l2];
      }
      nstart = (int64_t *)malloc((size_t)ncluster * sizeof *nstart);
      nlen = (int32_t *)malloc((size_t)ncluster * sizeof *nlen);
      nstore_cap = (int64_t)1 << 22; /* 16 млн слотов = 64 МБ; сверх — прежний путь */
      nstore = (int32_t *)malloc((size_t)nstore_cap * sizeof *nstore);
      if (!nstart || !nlen || !nstore) {
        rc = 2;
        goto done;
      }
      for (int64_t ci = 0; ci < ncluster; ci++)
        nstart[ci] = -1;
    }
    if (!pstamp || !bstart) {
      rc = 2;
      goto done;
    }
    /* bbox-индекс А1566: слот s пересекает клетку листа, если bbox его
     * ИСХОДНОГО треугольника задевает клетку (два прохода counting-sort) */
    {
      int32_t li3;
      int64_t total = 0;
      int fb = (o->walk != 0); /* А1578: фолбэк соседних листьев — только walk */
      for (int32_t s = 0; s < nt; s++)
        sw_index_piece(py, o, s, fb, bstart, NULL, NULL);
      for (li3 = 0; li3 < py->nleaf; li3++) {
        bstart[li3 + 1] += bstart[li3];
      }
      total = bstart[py->nleaf];
      bpids = (int32_t *)malloc((size_t)(total > 0 ? total : 1) * sizeof *bpids);
      fillb = (int32_t *)malloc((size_t)(py->nleaf + 1) * sizeof *fillb);
      if (!bpids || !fillb) {
        rc = 2;
        goto done;
      }
      for (li3 = 0; li3 < py->nleaf; li3++)
        fillb[li3] = bstart[li3];
      for (int32_t s = 0; s < nt; s++)
        sw_index_piece(py, o, s, fb, NULL, bpids, fillb);
      free(fillb);
      fillb = NULL;
    }
  }

  for (it = 0; it < o->iters; it++) {
    double emitted = 0, absorbed = 0, lost = 0, recycled = 0, e_sum = 0, area_sum = 0;
    int32_t p;
    memset(Ed, 0, (size_t)nt * sizeof *Ed);
    for (p = 0; p < nt; p++)
      Eprev[p] = py->pcs[p].e;
    if (o->mode == 3 && o->agg && ag.ng > 0) { /* §863: ΣEprev·a групп */
      for (int64_t g = 0; g < ag.ng; g++) {
        double s = 0;
        for (int32_t q = 0; q < ag.ggcnt[g]; q++) {
          int32_t m = ag.gmem[ag.ggoff[g] + q];
          s += Eprev[m] * area[m];
        }
        ag.gSE[g] = s;
      }
    }
    for (d = 0; d < nd; d++) {
      const double *om = tab[d].om;
      double w_d = tab[d].w;
      if (o->mode == 0 && o->ndirs == 26)
        w_d = HZ_SW_W; /* §835-паритет: скалярный путь нёс вес 4π/26 на шести осях */
      if (o->mode == 2 && o->vc) {
        /* §851: ЛУЧЕВОЙ СВИП — трубки-линии сетки, старт тёмный на
         * входных гранях; латеральная изоляция и стороны возникают сами */
        sw_line_acc la = {0, 0, 0, 0, 0, 0, 0};
        int32_t mark = (int32_t)(it * (nd + 1) + d + 1);
        sw_line_sweep_dir(py, om, w_d, o->le, kd, area, nrm, Eprev, Ed, o->rho, o->trivert,
                          o->noprop, stampv, mark, &la);
        absorbed += la.absorbed;
        emitted += la.emitted;
        recycled += la.recycled;
        lost += la.lost;
        st->nvisit += la.nvisit;
        st->ndep += la.ndep;
        st->traffic = la.ncell * 40; /* §851: грубая модель трафика лучевого свипа */
        continue;
      }
      if (o->mode == 3) {
        /* §852: МАТЕРИАЛЬНЫЙ ФРОНТ НА ПИРАМИДЕ — иерархический марш:
         * ПУСТ-прыжок по tmax узла / материальное взаимодействие на
         * LOD-уровне куска / спуск-продолжение DDA (А1563/А1566/А1569). */
        front_ctx fc;
        int32_t mark = (int32_t)(it * (nd + 1) + d + 1);
        memset(&fc, 0, sizeof fc);
        fc.py = py;
        fc.o = o;
        fc.area = area;
        fc.nrm = nrm;
        fc.kd = kd;
        fc.Eprev = Eprev;
        fc.Ed = Ed;
        fc.om = om;
        fc.w_d = w_d;
        fc.le = o->le;
        fc.lep = o->lep; /* А1576: per-piece эмиссия (NULL — прежний мир) */
        fc.noprop = o->noprop;
        fc.tau0 = o->tau0;
        fc.pstamp = pstamp;
        fc.bstart = bstart;
        fc.bpids = bpids;
        fc.mark = mark;
        fc.nstart = nstart;
        fc.nlen = nlen;
        fc.nstore = nstore;
        fc.nstore_n = &nstore_n;
        fc.nstore_cap = nstore_cap;
        fc.ncluster = (int32_t)ncluster;
        fc.ncluster0 = py->nleaf;
        fc.fld = fld; /* §864/Б1: кэш размеров уровней */
        fc.fldok = fldok;
        fc.agg = o->agg;
        fc.ag = &ag;
        for (int32_t l2 = 0; l2 < py->nlev; l2++)
          fc.noff[l2] = noff[l2];
        fc.cbase = (it == 0); /* прибор А1569 — геометрия статична, хватит раза */
        fc.cfront = (it == 0);
        fc.hi[0] = o->domhi[0]; /* меш-граница, не сеточная (§852) */
        fc.hi[1] = o->domhi[1];
        fc.hi[2] = o->domhi[2];
        front_dir(&fc, stampv);
        absorbed += fc.absorbed;
        emitted += fc.emitted;
        recycled += fc.recycled;
        lost += fc.lost;
        st->g6viol += fc.g6viol;
        st->njump += fc.njump;
        st->nmat += fc.nmat;
        st->ndesc += fc.ndesc;
        st->ncellbase += fc.ncellbase;
        st->ncellfront += fc.ncellfront;
        st->nstamp += fc.nstamp;
        st->nlostseg += fc.nlostseg;
        st->traffic = (fc.nmat + fc.ndesc + fc.njump) * 32 + fc.ncellfront * 40;
        free(fc.pbuf);
        free(fc.abuf); /* скретч растёт только на узловом пути (lp>=2, А1573) */
        free(fc.wbuf);
        free(fc.pbuf_t);
        free(fc.pbuf_p);
        continue;
      }
      double L = 0.0;
      int64_t prevline = -1;
      if (o->mode == 1) {
        const hz_sw_walk *w = &walks[d];
        int32_t k;
        int64_t run = 0; /* §843: визитов в текущей линии */
        for (k = 0; k < w->n; k++) {
          const hz_sw_wleaf *lf = &w->leaf[k];
          double Lsurf = 0.0;
          int64_t e, eend = lf->first + lf->npcs;
          int64_t pf = lf->first + 8 < nt ? lf->first + 8 : nt - 1;
          if (lf->line != prevline) { /* новая линия — тёмный вход */
            if (run == 1)
              st->nline1++;
            else if (run == 2)
              st->nline2++;
            if (run > 0) st->nline++;
            run = 0;
            L = 0.0;
            prevline = lf->line;
          }
          if (o->noprop) L = 0.0;
          if (lf->cntot <= 0.0) continue; /* куски встык лучу — слоя нет */
          st->nvisit++;                   /* §843: визит с слоем */
          __builtin_prefetch(&w->leaf[k + 8 < w->n ? k + 8 : w->n - 1]);
          __builtin_prefetch(&Eprev[w->wpid[pf]]);
          for (e = lf->first; e < eend; e++) {
            p = w->wpid[e];
            double rho = o->rho < 0 ? w->kd[e] : o->rho;
            Lsurf += (o->le + rho * Eprev[p] / (2.0 * M_PI)) * w->wcn[e];
            if (L > 0.0) st->ndep++;     /* §843: доставка света */
            Ed[p] += w_d * L * w->wn[e]; /* непрозрачный слой: инфлюкс целиком */
          }
          Lsurf /= lf->cntot;
          absorbed += w_d * L * csec;
          emitted += w_d * o->le * csec;
          recycled += w_d * (Lsurf - o->le) * csec;
          L = Lsurf; /* непрозрачный слой: луч гасится, остаётся переизлучение */
          if (o->noprop) L = 0.0;
          run++;
        }
        if (run == 1)
          st->nline1++;
        else if (run == 2)
          st->nline2++;
        if (run > 0) st->nline++;
        lost += L * csec;
      } else if (o->mode == 2) {
        /* §845/§846 (vc=0): объёмный фронт по кусковым листьям, L несётся
         * скаляром (однокомпонентная логика — только для сличения; А1551:
         * пересвечивает тонкие оболочки, латерально перемешивает). */
        const hz_sw_walk *w = &walks[d];
        int32_t k;
        L = 0.0;
        for (k = 0; k < w->n; k++) {
          const hz_sw_wleaf *lf = &w->leaf[k];
          double Lsurf = 0.0, Lin;
          int64_t e, eend = lf->first + lf->npcs;
          int64_t pf = lf->first + 8 < nt ? lf->first + 8 : nt - 1;
          if (o->noprop) L = 0.0;
          Lin = L;
          if (lf->cntot <= 0.0) continue; /* куски встык лучу — слоя нет */
          st->nvisit++;
          __builtin_prefetch(&w->leaf[k + 8 < w->n ? k + 8 : w->n - 1]);
          __builtin_prefetch(&Eprev[w->wpid[pf]]);
          for (e = lf->first; e < eend; e++) {
            p = w->wpid[e];
            double rho = o->rho < 0 ? w->kd[e] : o->rho;
            Lsurf += (o->le + rho * Eprev[p] / (2.0 * M_PI)) * w->wcn[e];
            if (Lin > 0.0) {
              st->ndep++;
              Ed[p] += w_d * Lin * w->wn[e]; /* инфлюкс клетки, доля по площади */
            }
          }
          Lsurf /= lf->cntot;
          /* §846: ПЕРВЫЙ ПОРЯДОК — перехват трубки пластинами. Тень кусков
           * T_tube = cntot/csec (для одной пластины точно; перекрытия не
           * вычитаются — консервативно, А1530); прошедшая доля луча идёт
           * сквозь клетку неперемешанной, переизлучённая — средним
           * радиансом кусков. Стена (cntot ≥ csec): T=0 → L_out=Lsurf,
           * как в §845; вакуум: T=1, луч не тронут. Депозиты выше уже
           * согласованы: их сумма = (1−T)·входной поток трубки. */
          {
            double T_cell = 1.0 - lf->cntot / csec;
            if (T_cell < 0.0) T_cell = 0.0;
            absorbed += w_d * (1.0 - T_cell) * Lin * csec;
            emitted += w_d * o->le * (1.0 - T_cell) * csec;
            recycled += w_d * (Lsurf - o->le) * (1.0 - T_cell) * csec;
            L = T_cell * Lin + (1.0 - T_cell) * Lsurf;
          }
        }
        lost += L * csec;
      } else {
        int32_t k;
        for (k = 0; k < py->nleaf; k++) {
          const hz_pyr_leaf *lf = &py->leaf[order[k]];
          double cntot = 0.0, T, Lsurf = 0.0;
          int32_t u;
          if (o->noprop) L = 0.0;
          for (u = 0; u < lf->npcs; u++) {
            p = py->csr[lf->pcs_first + u];
            double d0 = om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                        om[2] * nrm[3 * (int64_t)p + 2];
            cntot += area[p] * fabs(d0);
          }
          T = o->tau0 ? 1.0 : exp(-cntot / csec);
          if (cntot > 0.0) {
            for (u = 0; u < lf->npcs; u++) {
              p = py->csr[lf->pcs_first + u];
              double rho = o->rho < 0 ? kd[p] : o->rho;
              double cnm = fabs(om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                                om[2] * nrm[3 * (int64_t)p + 2]);
              Lsurf += (o->le + rho * Eprev[p] / (2.0 * M_PI)) * area[p] * cnm;
              Ed[p] += w_d * L * (1.0 - T) * cnm;
            }
            Lsurf /= cntot;
            absorbed += w_d * L * (1.0 - T) * csec;
            emitted += w_d * o->le * (1.0 - T) * csec;
            recycled += w_d * (Lsurf - o->le) * (1.0 - T) * csec;
            L = L * T + Lsurf * (1.0 - T);
          }
          if (o->noprop) L = 0.0;
        }
        lost += L * csec;
      }
    }
    if (o->mode == 3 && o->agg && ag.ng > 0) { /* §863: раздача вкладов групп */
      for (int64_t g = 0; g < ag.ng; g++) {
        for (int32_t q = 0; q < ag.ggcnt[g]; q++) {
          int32_t m = ag.gmem[ag.ggoff[g] + q];
          Ed[m] += ag.gacc[g]; /* edep одинаков всем членам: Σ = вклад группы */
          if (o->lpacc) o->lpacc[m] += o->cdelta * (o->rho >= 0.0 ? o->rho : kd[m]) * ag.gwt[g];
        }
        ag.gacc[g] = 0.0;
        ag.gwt[g] = 0.0;
      }
    }
    for (p = 0; p < nt; p++) {
      py->pcs[p].e = (float)Ed[p];
      e_sum += Ed[p] * area[p];
      area_sum += area[p];
    }
    st->e_avg = e_sum / area_sum;
    st->emitted = emitted;
    st->recycled = recycled;
    st->absorbed = absorbed;
    st->lost = lost;
    if (o->mode != 3) st->traffic = (int64_t)nd * ((int64_t)nt * 36 + (int64_t)walks[0].n * 40);
    if (e_hist) e_hist[it] = st->e_avg;
    if (o->mode == 3 && o->lpacc && o->lpapply) { /* §862: этаж = целая часть */
      int32_t nup = 0;
      if (!lpflo) {
        lpflo = (uint8_t *)malloc((size_t)nt);
        if (!lpflo) {
          rc = 2;
          goto done;
        }
      }
      for (p = 0; p < nt; p++) {
        double fl = o->lpacc[p];
        lpflo[p] = fl >= 255.0 ? 255 : (uint8_t)fl;
      }
      if (hz_pyr_set_lp(py, lpflo, &nup) != 0) {
        rc = 2;
        goto done;
      }
    }
  }

done:
  free(lpflo); /* §862 */
  free(fld);   /* §864/Б1 */
  free(fldok);
  free(Ed);
  free(Eprev);
  free(order);
  free(vindex);
  free(bits);
  free(tab);
  free(pstamp); /* §852/А1564 */
  free(stampv);
  free(nstart); /* §863 */
  free(nlen);
  free(nstore);
  free(bstart);
  free(bpids);
  if (walks) {
    for (d = 0; d < nd; d++)
      sw_free_walk(&walks[d]);
    free(walks);
  }
  return rc;
}

/* §843: общий кэш таблицы для сумм по направлениям. НЕ потокобезопасно —
 * вызывается из однопоточного инструмента swee3. */
static int sw_cached_table(int ndirs, const hz_sw_dir **tab, int *nd) {
  static hz_sw_dir *cache = NULL;
  static int cached_ndirs = -1, cached_nd = 0;
  if (ndirs != cached_ndirs) {
    int rc, n;
    hz_sw_dir *t = NULL;
    rc = hz_sw_dir_table(ndirs, &t, &n);
    if (rc != 0) return rc;
    free(cache);
    cache = t;
    cached_ndirs = ndirs;
    cached_nd = n;
  }
  *tab = cache;
  *nd = cached_nd;
  return 0;
}

double hz_sw_dirsum(const double n[3], int ndirs) {
  const hz_sw_dir *tab;
  int nd, rc, d;
  double s = 0.0;
  rc = sw_cached_table(ndirs, &tab, &nd);
  if (rc != 0) return 0.0;
  for (d = 0; d < nd; d++)
    s += tab[d].w * fabs(tab[d].om[0] * n[0] + tab[d].om[1] * n[1] + tab[d].om[2] * n[2]);
  return s;
}

/* §842: Σ по ВЫХОДНЫМ направлениям (n·ω > 0) w·|ω·n| — облучённость куска
 * в модели слоя при тёмном окружении. Для точного решения фикстуры. */
double hz_sw_exitwsum(const double n[3], int ndirs) {
  const hz_sw_dir *tab;
  int nd, rc, d;
  double s = 0.0;
  rc = sw_cached_table(ndirs, &tab, &nd);
  if (rc != 0) return 0.0;
  for (d = 0; d < nd; d++) {
    double c = tab[d].om[0] * n[0] + tab[d].om[1] * n[1] + tab[d].om[2] * n[2];
    if (c > 0) s += tab[d].w * c;
  }
  return s;
}
