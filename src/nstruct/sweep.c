/* sweep.c — свип на новых массивах: §835 (скалярный луч), §836-§838
 * (колонки/линии, материализованный обход), §841 (гигиена цикла).
 * См. sweep.h.
 */
#include "sweep.h"

#include <math.h>
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
/* §851: ЛУЧЕВОЙ СВИП ОДНОГО НАПРАВЛЕНИЯ. Трубки = линии сетки вдоль om:
 * стартуют тёмными (L=0) на входных гранях сетки, шаг DDA по клеткам;
 * пустота транслирует L точно (свободный пробег — точное решение на
 * вакууме); в клетке с кусками: депозит w·L·|cos| каждому куску,
 * перехват T = max(0,1−cntot/csec), переизлучение Lsurf (А1530/А1531).
 * Латеральная изоляция трубок возникает сама: линия стартует тёмной и
 * умирает на выходе — «компоненты пустоты» не нужны как разметка
 * (А1551). stamp исключает двойные визиты клетки линиями с общей
 * входной гранью. */
typedef struct {
  double absorbed, emitted, recycled, lost;
  int64_t nvisit, ndep, ncell;
} sw_line_acc;

static void sw_line_sweep_dir(const hz_pyr *py, const double om[3], double w_d, double le,
                              const double *kd, const double *area, const double *nrm,
                              const double *Eprev, double *Ed, double csec, double rho_ovr,
                              int noprop, int32_t *stamp, int32_t mark, sw_line_acc *acc) {
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
        double tmax[3], tdelta[3], L = 0.0;
        int32_t a2;
        cc[ax] = f[ax];
        cc[bxA] = b1;
        cc[bxB] = b2;
        cell[0] = cc[0];
        cell[1] = cc[1];
        cell[2] = cc[2];
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
          double tn;
          acc->ncell++;
          stamp[id] = mark;
          if (pos >= 0) {
            const hz_pyr_leaf *lf = &py->leaf[pos];
            double Lsurf = 0.0, cntot = 0.0, T_cell;
            int32_t u;
            for (u = 0; u < lf->npcs; u++) {
              int32_t p = py->csr[lf->pcs_first + u];
              double d0 = om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                          om[2] * nrm[3 * (int64_t)p + 2];
              cntot += area[p] * fabs(d0);
            }
            if (cntot > 0.0) {
              acc->nvisit++;
              for (u = 0; u < lf->npcs; u++) {
                int32_t p = py->csr[lf->pcs_first + u];
                double d0 = om[0] * nrm[3 * (int64_t)p] + om[1] * nrm[3 * (int64_t)p + 1] +
                            om[2] * nrm[3 * (int64_t)p + 2];
                double an = fabs(d0);
                double rho = rho_ovr < 0 ? kd[p] : rho_ovr;
                Lsurf += (le + rho * Eprev[p] / (2.0 * M_PI)) * area[p] * an;
                if (L > 0.0) {
                  acc->ndep++;
                  Ed[p] += w_d * L * an;
                }
              }
              Lsurf /= cntot;
              T_cell = 1.0 - cntot / csec;
              if (T_cell < 0.0) T_cell = 0.0;
              acc->absorbed += w_d * (1.0 - T_cell) * L * csec;
              acc->emitted += w_d * le * (1.0 - T_cell) * csec;
              acc->recycled += w_d * (Lsurf - le) * (1.0 - T_cell) * csec;
              L = T_cell * L + (1.0 - T_cell) * Lsurf;
              if (noprop) L = 0.0; /* НК: фронт не переносится */
            }
          }
          tn = tmax[0];
          if (tmax[1] < tn) tn = tmax[1];
          if (tmax[2] < tn) tn = tmax[2];
          if (tn > 1e29) break;
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
        if (noprop) L = 0.0;
        acc->lost += L * csec;
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
  int nd = 0, d, it, rc = 0;
  double csec;

  memset(st, 0, sizeof *st);
  if (!py || !area || !nrm || !kd || !o) return 1;
  if (o->iters <= 0) return 1;
  if (!py->pcs || nt != py->nt) return 1;
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
  if (!Ed || !Eprev || !order || !vindex || !bits || !walks) {
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
  if (o->mode == 2 && o->vc) {
    int64_t ncells = (int64_t)py->nx * py->ny * py->nz;
    stampv = (int32_t *)calloc((size_t)ncells, sizeof *stampv);
    if (!stampv) {
      rc = 2;
      goto done;
    }
  }

  for (it = 0; it < o->iters; it++) {
    double emitted = 0, absorbed = 0, lost = 0, recycled = 0, e_sum = 0, area_sum = 0;
    int32_t p;
    memset(Ed, 0, (size_t)nt * sizeof *Ed);
    for (p = 0; p < nt; p++)
      Eprev[p] = py->pcs[p].e;
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
        sw_line_sweep_dir(py, om, w_d, o->le, kd, area, nrm, Eprev, Ed, csec, o->rho, o->noprop,
                          stampv, mark, &la);
        absorbed += la.absorbed;
        emitted += la.emitted;
        recycled += la.recycled;
        lost += la.lost;
        st->nvisit += la.nvisit;
        st->ndep += la.ndep;
        st->traffic = la.ncell * 40; /* §851: грубая модель трафика лучевого свипа */
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
    st->traffic = (int64_t)nd * ((int64_t)nt * 36 + (int64_t)walks[0].n * 40);
    if (e_hist) e_hist[it] = st->e_avg;
  }

done:
  free(Ed);
  free(Eprev);
  free(order);
  free(vindex);
  free(bits);
  free(tab);
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
