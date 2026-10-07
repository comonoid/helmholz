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
/* §896-5: бюджет ног (зеркальных хопов) на направление. Число 2000000 —
 * из §896-5. Сейчас 0: механизм ног ВЫКЛЮЧЕН, и при нулевом бюджете вся нога
 * сразу уходит в lost («бюджет исчерпан — остаток в lost»). Проверено 07-10:
 * значение 2000000 присваивалось ЛОКАЛЬНОЙ переменной в hz_sw_run и в fc не
 * попадало, а fc обнуляется memset — то есть ноги были выключены НЕЯВНО, а
 * присваивание было мёртвым (отсюда -Wunused-but-set-variable). Здесь
 * состояние названо явно. Включение ног — отдельное решение владельца: оно
 * меняет физику зеркальных хопов (ksf=1 — продакшн-дефолт §875).
 * Переопределяется ТОЛЬКО для замера последствий, без правки кода:
 *   gcc ... -DHZ_SW_LEG_BUDGET=2000000 ...   (make не переопределяет) */
#ifndef HZ_SW_LEG_BUDGET
#define HZ_SW_LEG_BUDGET 0
#endif

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
/* §865/раунд 15: запись per-трубочного кэша «кусок проверен» (path=1) */
typedef struct {
  double bh0, bh1; /* bbox t-интервал куска на трубке (пустой: 1,0) */
  double tt;       /* кэш sw_ray_tri_raw */
  int64_t tstamp;  /* pkey: bh0/bh1 валидны */
  int64_t ttstamp; /* pkey: tt посчитан (ray-tri — один раз на кусок×трубку) */
  int64_t pstamp;  /* pkey: кусок уже в списке попаданий этой трубки */
} pc_rec;
typedef struct {
  int64_t ng, ngcap;
  int32_t *grep;            /* представитель группы */
  double *gA, *gleA, *grA;  /* площадь, средние le и ρ (статичные) */
  double *gSE, *gacc, *gwt; /* ΣEprev·a, Σedep, Σedep·(1−cosθ) — на итерацию */
  int64_t *ggoff;
  int32_t *ggcnt, *gmem;
  int64_t gmem_n, gmem_cap;
} sw_agg;

/* §923: списки кусков ОДНОГО уровня владения. Кусок этажа f лежит ровно в
 * одном месте: f == 0 — в списке своей листовой клетки, f ≥ 1 — в списке
 * узла-предка уровня f-1. Узел видит ТОЛЬКО свои куски — чужих списков
 * нет, итерации поддерева на визит нет (урок §922-ДОКЛАД-2). */
typedef struct {
  int32_t *off;   /* [n+1] CSR-смещения по узлам уровня (в lvl_pool) */
  int32_t *pids;  /* база пула кусков этого уровня (сдвиг внутри lvl_pool) */
  uint8_t *finer; /* [n] под узлом есть куски этажа ≤ l — нужен спуск */
  int32_t n;      /* число узлов уровня */
} sw_lvl;

typedef struct {
  const hz_pyr *py;
  const hz_sw_opts *o;
  const double *area, *nrm, *kd, *Eprev;
  double *Ed;
  const double *om;
  double omcur[3]; /* §873/T4: мутируемое направление ТЕКУЩЕЙ ноги (хоп
                    * переписывает); om указывает сюда после setup — общая
                    * таблица направлений только для чтения */
  double org[3];   /* центр входной клетки трубки */
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
  int64_t noff[256]; /* смещение уровня в нумерации кластеров узлов */
  int32_t ncluster0; /* nleaf: кластеры листов идут первыми */
  int agg;           /* §863: ключ o->agg */
  int agg_list;      /* текущий список — узловой (отсортирован, годен для групп) */
  /* §923: раскладка по уровням владения (lvls == NULL — прежний мир:
   * без lpacc/lpapply и в статике lp=; А1638 fail closed) */
  const sw_lvl *lvls;       /* [py->nlev] */
  const int32_t *leaf_loff; /* [nleaf+1] листовые списки floor-0 кусков */
  const int32_t *leaf_lpids;
  /* §924: представители декимации (reps_on — только живая лестница+walk).
   * Списки этажей держат id ≥ nt0 — расширенные tri id представителя. */
  int reps_on;
  int32_t nt0; /* исходные nt (граница слот/tri-пространств) */
  int32_t nrep;
  const int32_t *rep_koff; /* [nrep+1] */
  const int32_t *rep_kmem; /* [Σдетей] СЛОТы исходных кусков */
  const double *rep_area;  /* [nrep] */
  const double *rep_rho;   /* [nrep] */
  const double *rep_le;    /* [nrep] */
  const double *rep_atri;  /* [nrep] площадь носителя (§927) */
  const double *rep_nrm;   /* [3nrep] нормаль носителя (§927) */
  const double *rep_ep;    /* [nrep] агрегат Eprev детей (раз в итерацию) */
  int ev_src_far;          /* §931/А1666: текущее событие — дальнее (rep);
                            * фильтр этажей ближних кусков (ключ rep_near_lp) */
  /* §928: калибровка k(r,ω) — ω-бин направления, накопители кусочной ветки
   * (только последняя итерация, только мир сбора HZ_KCALDUMP), этажи для
   * отсева floor-0 событий (они в reps-мире приходят листами, не rep-хитом). */
  int dir_bin;                       /* ω-бин = уровень μ квадратуры (HZ_KCAL_NW×1) */
  double *knum;                      /* [nrep*HZ_KCAL_NW] числитель: Σ dep детей */
  double *kden;                      /* [nrep*HZ_KCAL_NW] знаменатель: Σ w·csec·axcos·Lin_entry */
  const uint8_t *kfloor;             /* [nt] этажи (слоты) или NULL */
  int64_t dbg_node, dbg_leaf, dbg_n; /* §923-дых: визиты/куски нового пути */
  sw_agg *ag;                        /* §863/шаг 2: таблица групп (NULL — путь не активен) */
  double *abuf, *wbuf;               /* скретч взаимодействия: без аллокаций на событие */
  int64_t abufcap;
  double hi[3]; /* верх сцены (клетки уровней шире домена на нечётных сетках) */
  /* §864/Б1: кэш размеров уровней (статичны на прогон; иначе level_dims —
   * 22–37 млн вызовов). Владелец — hz_sw_run, fc только ссылается. */
  int64_t (*fld)[3];
  uint8_t *fldok;
  double *pbuf_t;
  int32_t *pbuf_p;
  int64_t pbuf_n, pbuf_cap;
  /* §865/раунды 13+15: per-трубочный кэш «кусок проверен» — ОДНА запись на
   * кусок (40 Б, одна-две кэш-линии вместо 3–4 случайных по массивам):
   * bbox t-интервал по ВСЕЙ трубке + кэш ray-tri + штампы pkey.
   * Владелец hz_sw_run; NULL → прежний цикл (fail closed). */
  pc_rec *pc;        /* [nt] */
  double pth0, pth1; /* границы трубки [h0,h1] для кэша (§865/раунд 13) */
  double lost;       /* поток за границей домена, w·csec-единицы */
  int negseen;       /* §871: события Lin<0 за прогон (сентинел, лимит печати 8) */
  /* §873/T4 шаг 2: очередь зеркальных хопов текущей трубки. Глубина
   * ограничена HZ_MIRROR_BOUNCE_MAX: при ks ≤ 0.95 вклад (n+1)-го
   * отражения < 0.95ⁿ < 0.82 уже при n=4; невлезающие доли идут в lost
   * (баланс не нарушается). */
  int hop_n;                /* элементов в очереди */
  int hop_depth;            /* выполнено хопов (диагностика; усечения больше нет) */
  double hop_lost;          /* доля, вытесненная порогом/ёмкостью (в lost) */
  double hop_lost_sum;      /* §881: Σ по трубкам направления */
  int64_t hop_hops_sum;     /* §881: Σ хопов по трубкам направления */
  double Lh_last;           /* §892: диффузное продолжение клетки */
  double lh_cap;            /* §893: кап радианса итерации */
  double lh_seen;           /* §893: максимум радианса за направление */
  double hop_spawn_e;       /* §896-3: Σ ks·Lin при спавне хопов */
  double *Linmax_prev;      /* §893-b: пер-кусковый max пришедшего радианса */
  double *Linmax_cur;       /* §893-b: текущей итерации */
  int in_leg;               /* §894-c: трубка — нога (пер-хит депозит) */
  double dep_leg, dep_main; /* §894: Σ депозитов ног / основной трубки */
  int64_t *leg_pstamp;      /* §894-c-5: штампы ног — отдельное пространство */
  uint64_t leg_pkey;        /* §894-c-5: ключ штампа ноги (уникален на ногу) */
  double leg_root;          /* §894-c-5: корень цепи ноги (для порога) */
  int hop_dep[32];          /* §896-4: глубина цепи каждого хопа */
  int64_t leg_budget;       /* §896-5: бюджет ног на направление */
  int cur_leg_depth;        /* §896-4: глубина обрабатываемой ноги */
  int64_t dep_cnt;          /* §894: число депозитов */
  double t_last;            /* §932-А: t последнего материального события трубки
                             * (−1 — не было; 0 на ноге хопа); монотонный max */
  double t_cur;             /* §932-А: t текущего хита — проводка в sw_accum */
  double t_lv;              /* §932-А: ПРОБЕГ текущего события = t_cur − t_last
                             * НА МОМЕНТ хита (до обновления t_last; <0 — нет
                             * события-предшественника/инверсия — вклада нет) */
  double ldom;              /* §932-А: Ldom, диагональ меша (масштаб f=Lv/Ldom) */
  /* §932-Б: статистическая среда (агрегаты — владелец hz_sw_run,
   * пересчитываются со sw_levels_build, А1711; fc только ссылается) */
  int med_on;              /* среда активна (агрегаты построены) */
  const double *med_box;   /* [6nrep] bbox объединения замещённых детей */
  const double *med_an;    /* [3nrep] Σ A_i·n_i замещённых (σ = |an·ω|/V) */
  const double *med_v;     /* [nrep] V_r — объём bbox */
  const double *med_rho;   /* [nrep] площадь-среднее ρ_r */
  const double *med_asum;  /* [nrep] Σ площадей замещённых (нормировка) */
  const int32_t *med_nsub; /* [nrep] замещённых детей (0 — обычный rep) */
  const int32_t *med_koff; /* [nrep+1] CSR замещённых детей */
  const int32_t *med_kmem; /* [Σ замещённых] слоты детей */
  int64_t *med_stamp;      /* [nrep] pkey активации (повторы гасит) */
  double *med_pend;        /* [nrep] Σ отложенных депозитов направления */
  double *ced;             /* [nleaf] §932-Б/Ш9-ЗОНД + Ш12: клеточный
                            * приём. cellc=1 (RAW-зонд): Φ_входа звена;
                            * cellc=2: ИЗВЛЕЧЁННЫЙ поток Φ_вх·(1−e^{−Λ})
                            * (Λ — сумм. доля перехвата, §932-Б/Ш12-ПЛАН);
                            * раздача детям по площадям — в конце итерации.
                            * NULL — клеточный приём выключен */
  int cellc;               /* §932-Б/Ш12: режим клеточного приёма 0/1/2 */
  double *cabs;            /* [nleaf] §932-Б/Ш12 (только cellc=2): радианс
                            * НА ВХОДЕ клетки для этого направления (Lin в
                            * начале звена, для доли переизлучения клетки) */
  struct {
    int32_t r;
    double tin, tout, sig, rho, tprev, pend;
  } med_act[8]; /* активные среды трубки (порядок = порядок входа = t_in);
                 * pend — отложенный депозит (радианс), флеш per-region */
  int med_nact;
  double hop_lost_thr;       /* §881/П7: из hop_lost — порогом */
  double hop_lost_cap;       /* §881/П7: из hop_lost — ёмкостью */
  double row_dep;            /* §887-b: счётчик исполнений row-ветки */
  int64_t lin_pos, lin_zero; /* §887-d: депозиты с Lin>0.5 / <=0.5 */
  /* §911-6 прибор покрытия (HZ_COVDBG): Σ долей следа (Σ dep_ratio) и число
   * хитов на (кусок, направление) — сверка с леммой покрытия Σ = μ_пов */
  double *cov;         /* [nt] Σ dep_ratio куска за направление; NULL — выкл */
  uint32_t *covn;      /* [nt] число зарегистрированных хитов за направление */
  double *cov_lit;     /* §911-8: Σ dep_ratio·[Lin>0] — освещённое покрытие */
  double *cov_e;       /* §911-10: Σ dep_ratio·Lin — вклад с амплитудой света */
  int64_t lh_cap_hits; /* §911-6: срабатывания капа радианса (A3) */
  /* §881: ёмкость 32 (была 4): на зеркально-плотных сценах цепи длиннее 4 —
   * основной поток (дефект §880); усечение — порогом от корня цепи (А1617) */
  double hop_pt[32][3];
  double hop_dir[32][3];
  double hop_lin[32];
  double hop_root[32]; /* Lin цепи в точке первого зеркального удара */
  int64_t ncellbase;   /* базовых клеток полным DDA («до», прибор А1569) */
  int cbase, cfront;   /* прибор считается на первой итерации (геометрия статична) */
  /* аккумуляторы */
  double depA; /* Σ родительской доли, депонированной кускам (L-единицы) */
  double absorbed, emitted, recycled;
  int64_t g6viol, njump, nmat, ndesc, ncellfront, nstamp, nlostseg, ndep;
  /* Ш12: инкремент поглощённой СРЕДОЙ доли (радианс-единицы) — банк
   * sw_med_bank добавляет; разность снимков = доля среды на ЭТОМ звене.
   * Клеточный приём обязан её учесть: иначе среда и клетка вычтут из трубки
   * до 2× её потока (консервативность ломается). */
  double med_rem;
} front_ctx;

/* отрезок [h0,h1] луча (org, om) в коробке, пересечённый с [tin,tout];
 * 0 — пусто */
/* §930-В0: always_inline по замеру — нативная стена свипа −12%
 * (зонд always_inline против обычного O2, храм dirs=4x4 it=2, 2+2
 * прогона); вызовов 600M+, компилятор сам не инлайнит (раздувание). */
static inline __attribute__((always_inline)) int
front_box_seg(const double blo[3], const double bhi[3], const double org[3], const double om[3],
              double tin, double tout, double *h0, double *h1) {
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
    /* СТАТИКА (§925-с): t0 по осям только растёт, t1 — только сжимается
     * (max/min коммуникативны) — выход сразу после оси с t0 > t1 БИТОВО
     * равен проходу всех трёх: оставшиеся оси знак не развернут.
     * Экономит до двух делений на отвергнутом тесте. */
    if (t0 > t1) return 0;
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

/* §887-c: знаменатель депозита — «минимальная доля ряда» (форма §735):
 * вырожденный кусок (area -> 0) получал E ~ 1/area и взрывался
 * самоподкачкой Lh (room: E_avg 74.9, max куска 1e7 при it=20 против
 * сошедшегося MC 11.07). DEP_MIN_SHARE = 1/16: кусок не может принять более
 * 16 долей среднего по клетке; для обычных кусков D == area — арифметика
 * битово неизменна. */
#define FRONT_DEP_MIN_SHARE (1.0 / 16.0)
static double front_depden(const front_ctx *fc, int32_t p) {
  double cell2 = fc->py->cell * fc->py->cell;
  /* §911-13-3: у клип-полосок знаменатель — площадь родительского tri
   * (parea), не полоски: иначе все трубки tri агрегируются на одну полоску
   * со ставкой 1/D → усиление ×a_tri/D (клип-взрыв §911-3). */
  double d = fc->o->parea ? fc->o->parea[p] : fc->area[p];
  return d > FRONT_DEP_MIN_SHARE * cell2 ? d : FRONT_DEP_MIN_SHARE * cell2;
}

/* §887-d: трассировщик одного куска (HZ_TRACE_TRI) УДАЛЁН 07-10: он был
 * объявлен и никогда не вызывался (gcc -Wunused-function), а его состояние
 * tri_trace_id/tri_trace_init не читал никто, кроме него самого. Если
 * трассировка понадобится — вернуть из git-истории (коммит с ловлей
 * -Wunused-function) и вызвать в front_visit, а не держать мёртвым. */

/* §911-6: ограничитель радианса — Lh не превышает кап итерации */
#define LH_GROWTH (1.0 + 0.25) /* §893: рост радианса не быстрее +25 %/итерацию */
static int sw_covdbg_init = 0;
static int sw_covdbg_on = 0;
static int sw_covdbg(void) { /* HZ_COVDBG=1 — прибор покрытия §911-6 */
  if (!sw_covdbg_init) {
    const char *e = getenv("HZ_COVDBG");
    sw_covdbg_on = e && e[0] != '\0' && e[0] != '0';
    sw_covdbg_init = 1;
  }
  return sw_covdbg_on;
}
static double front_lh_cap(front_ctx *fc, double lh) {
  if (lh > fc->lh_cap) {
    lh = fc->lh_cap;
    fc->lh_cap_hits++; /* §911-6: счёт срабатываний (A3-дефицит цепи) */
  }
  if (lh > fc->lh_seen) fc->lh_seen = lh;
  return lh;
}

static double front_cos(const front_ctx *fc, int32_t p) {
  const double *om = fc->om;
  return fabs(om[0] * fc->nrm[3 * (int64_t)p] + om[1] * fc->nrm[3 * (int64_t)p + 1] +
              om[2] * fc->nrm[3 * (int64_t)p + 2]);
}

/* §862: приращение дробного аккумулятора детальности Δ(материал,угол):
 * темнее материал (ρ) и скользящее падение (1−cosθ, cosθ — к НОРМАЛИ куска)
 * — быстрее набор этажа. Нет lpacc — нет операции (битово прежний мир). */
static int64_t g932_n;                /* §932-А: события с пробегом (t_last≥0 ∧ Lv>0) */
static double g932_fsum;              /* §932-А: Σf (f=Lv/Ldom) — ⟨f⟩ для калибровки tvc */
static double g932_lvsum;             /* §932-А: ΣLv — ⟨Lv⟩ (стартовый tvc ~ Ldom/⟨Lv⟩) */
static int64_t g932_hist[16];         /* §932-А: гистограмма f, шаг 1/16 (А1724) */
static int g932_nomax = -1;           /* §932-А1715: env HZ_TVNOMAX — t_last=ht[i] без
                                       * max (НК с предсказанным провалом) */
static int64_t g932b_ev, g932b_big;   /* §932-Б: активаций среды; σℓ>3 (А1700-4) */
static double g932b_abs, g932b_dep;   /* §932-Б: Σ поглощённого; Σ депозитов детям */
static int64_t g932b_sact;            /* §932-Б-Ш3: активаций с σℓ>0 (статистика σℓ) */
static double g932b_ssum, g932b_smax; /* §932-Б-Ш3: Σσℓ и max σℓ по активациям */
static double g932b_edf;              /* §932-Б-Ш3: ФАКТ-Ed флеша (Σ_kid edp·A);
                                       * сверка с g932b_dep — разбор Вопроса-1:
                                       * флеш домножает Σpend на axcos последней
                                       * семьи, события несут axcos своей (Г-flux) */
static int sw_med_now = -1;           /* §932-Б-Ш3: HZ_MEDNOW — немедленная раздача
                                       * (ветка Ш1) для разбора Вопроса-1 */
static int sw_med_t1 = -1;            /* §932-Б-Ш3: HZ_MEDT1 — зонд времени:
                                       * T=1−x вместо exp (физика НЕВЕРНА, замер
                                       * доли exp в банке) */
static int sw_med_scr = -1;           /* §932-Б-Ш4: HZ_MEDSCRAMBLE — НК А1698:
                                       * детерминированный множитель σ по региону;
                                       * предсказание: баланс сохраняется
                                       * (конструктивно), поле дрейфует */
static int64_t g932b_nest;            /* §932-Б-Ш4: входов при живых средах
                                       * (вложенные регионы, Б0-в/А1720) */
static double g932c_dep;              /* §932-Б/Ш9-ЗОНД: Σ клеточного приёма */
static int64_t g932c_n;               /* §932-Б/Ш9-ЗОНД: куско-раздач за прогон */
static int64_t g932c_vis;             /* §932-Б/Ш12: клеточных изъятий за прогон */
static double g932c_fmax;             /* §932-Б/Ш12: max f изъятия (≤1 по построению) */
static int64_t g932c_ge1;             /* §932-Б/Ш12-прибор: изъятий с f ≥ 1−1e-12 */

/* §932-Б: ИНКРЕМЕНТНЫЙ БАНК СРЕДЫ (А1710) — ослабить Lin до момента t и
 * раздать поглощённое: (1−ρ_r) — замещённым детям по площадям (+depA,
 * А1733/G6), ρ_r — переизлучение в луч. Банк на каждом событии и
 * границе сегмента: Σ банков = интегралу Беера ТОЧНО (двойного счёта
 * нет: Δabs всегда от ТЕКУЩЕГО Lin). Активные среды — в порядке входа
 * (вложенные: внешняя раньше, ослабления перемножаются, Б0-в). */
static double sw_med_bank(front_ctx *fc, double Lin, double t, double csec) {
  int k, m = 0;
  for (k = 0; k < fc->med_nact; k++) {
    int32_t r = fc->med_act[k].r;
    double tout = fc->med_act[k].tout, tprev = fc->med_act[k].tprev;
    double te = tout < t ? tout : t;
    double dt = te - tprev;
    if (dt <= 0.0) { /* ещё жива (t < tprev не бывает: события сортированы) */
      fc->med_act[m++] = fc->med_act[k];
      continue;
    }
    {
      double x = fc->med_act[k].sig * dt;
      double T = sw_med_t1 ? 1.0 - x : exp(-x); /* Ш3: битово exp(-sig·dt) */
      double dabs = Lin * (1.0 - T);            /* радианс-единицы поглощённого СРЕДОЙ */
      fc->med_rem += dabs;                      /* Ш12: счётчик для клеточного приёма */
      double flux = fc->w_d * csec * fc->axcos * dabs;
      double fdep = flux * (1.0 - fc->med_act[k].rho); /* депозит детям: поток */
      g932b_abs += flux;
      g932b_dep += fdep;
      fc->absorbed += flux;
      fc->depA += dabs * (1.0 - fc->med_act[k].rho); /* А1733: радианс-единицы, как ai */
      if (sw_med_now) { /* Ш3/HZ_MEDNOW: немедленная раздача (ветка Ш1) —
                         * edp с axcos ЭТОЙ семьи, счётчику соответствует */
        double edp = fdep / fc->med_asum[r];
        int32_t q;
        for (q = fc->med_koff[r]; q < fc->med_koff[r + 1]; q++)
          fc->Ed[fc->med_kmem[q]] += edp;
      } else {
        /* Ш3/Г-flux ПРАВКА: pend копится в ЕДИНИЦАХ ПОТОКА (с axcos
         * события СВОЕЙ семьи осей) — флешу не нужно домножение на
         * axcos последней семьи (Вопрос-1: 0.5% E_avg и факт-Ed≠деп);
         * Σ_kid edp·A == g932b_dep по построению */
        fc->med_act[k].pend += fdep;
      }
      Lin = Lin * T + fc->med_act[k].rho * dabs; /* ρ_r — переизлучение в луч */
    }
    fc->med_act[k].tprev = te;
    if (t < tout)
      fc->med_act[m++] = fc->med_act[k]; /* ещё активна */
    else
      fc->med_pend[fc->med_act[k].r] += fc->med_act[k].pend; /* регион закрыт */
  }
  fc->med_nact = m;
  return Lin;
}

/* §932-Б: ВИРТУАЛЬНОЕ СОБЫТИЕ ВХОДА (А1695/А1709) — bbox региона на всю
 * трубку (front_box_seg, М2-детерминизм), σ = |Σ A n·ω|/V (А1737),
 * повторные входы гасит штамп (pkey). t — t события входа. */
static void sw_med_enter(front_ctx *fc, int32_t r, double t) {
  double h0, h1, sig, v;
  const double *an;
  if (fc->med_stamp[r] == fc->pkey) return; /* уже активирована этой трубкой */
  fc->med_stamp[r] = fc->pkey;
  v = fc->med_v[r];
  if (!(v > 1e-30)) return; /* вырожденный регион — среда не взаимодействует */
  /* Смещения региональных массивов — в 64 битах: 6·r и 3·r при nrep в сотни
   * миллионов переполнили бы int (bugprone-implicit-widening, ловля 07-10) */
  if (!front_box_seg(fc->med_box + 6 * (int64_t)r, fc->med_box + 6 * (int64_t)r + 3, fc->org,
                     fc->om, -1e30, 1e30, &h0, &h1))
    return;
  an = fc->med_an + 3 * (int64_t)r;
  sig = fabs(an[0] * fc->om[0] + an[1] * fc->om[1] + an[2] * fc->om[2]) / v;
  if (sw_med_scr) /* Ш4/А1698 НК: σ × (0.25…3.25) по региону, детерминированно */
    sig *= 0.25 + 3.0 * (double)((unsigned)r % 4u) / 4.0;
  if (!(sig > 0.0)) return;               /* просвет: T=1, поглощения нет — честно */
  if (fc->med_nact >= 8) return;          /* ёмкость: сверх — не активируем (счётчик) */
  if (fc->med_nact > 0) g932b_nest++;     /* Ш4: вход при живых средах (Б0-в/А1720) */
  if (sig * (h1 - h0) > 3.0) g932b_big++; /* А1700-4: вырождение в стену */
  { /* Ш3: статистика σℓ по активациям (битивно-нейтральный прибор) */
    double sl = sig * (h1 - h0);
    g932b_sact++;
    g932b_ssum += sl;
    if (sl > g932b_smax) g932b_smax = sl;
  }
  fc->med_act[fc->med_nact].r = r;
  fc->med_act[fc->med_nact].tin = h0;
  fc->med_act[fc->med_nact].tout = h1;
  fc->med_act[fc->med_nact].sig = sig;
  fc->med_act[fc->med_nact].rho = fc->med_rho[r];
  fc->med_act[fc->med_nact].tprev = h0 > t ? h0 : t;
  fc->med_act[fc->med_nact].pend = 0.0;
  fc->med_nact++;
  g932b_ev++;
}
/* §931: предикат «далеко» (А1663): центроид куска дальше rep_zone от
 * rep_eye. rep_zone<=0 / нет rep_cent — весь мир «далеко» = прежний §924,
 * битово. Квадрат расстояния — без sqrt (детерминизм, тот же порог, что
 * у farshare: далеко ⇔ d² > zone²). */
static int sw_far(const hz_sw_opts *o, int32_t p) {
  const double *c;
  if (o->rep_zone <= 0.0 || o->rep_cent == NULL) return 1;
  c = o->rep_cent + 3 * (int64_t)p;
  return (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]) > o->rep_zone * o->rep_zone ? 1 : 0;
}

static void sw_accum(front_ctx *fc, int32_t p, double edep) {
  if (!fc->o->lpacc) return;
  if (fc->o->rep_near_lp && fc->ev_src_far && !sw_far(fc->o, p))
    return; /* А1666: этажи БЛИЖНИХ кусков — только из ближних депозитов
             * (после А1662 rep-доли ближним не приходят вовсе — фильтр
             * страховочный, умолчание ВЫКЛ) */
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
    case 5:
      /* §918 (указание пользователя): ДИФФУЗНОСТЬ ОТСКАКА — Δ = 1−ks_eff:
       * один диффузный депозит уже повод снизить детальность, второй —
       * кардинально (floor аккумулятора растёт на ~1 за отскок);
       * зеркальный (ks→1) Δ→0 — детальность ОРИГИНАЛЬНАЯ. Нормируется
       * на единичный вклад (без edep): правило про счёт отскоков, не про
       * энергию; cdelta задаёт масштаб шага.
       * §924: ЭМИТТЕРЫ не грубеют (как зеркальные): агрегат repre-
       * зентателя разводит le по площади — источники гаснут (ловля
       * §924-дым: emitted=0, E 4.07→0.075). */
      if (fc->o->lep && fc->o->lep[p] > 0.0) {
        d = 0.0;
        break;
      }
      d = (fc->o->ks && fc->o->ks[p] >= 0.0) ? 1.0 - fc->o->ks[p]
                                             : 1.0; /* §918: ks нет/не задан — считаем диффузным */
      break;
    default:
      d = rho * (1.0 - ct);
      break; /* §862: форма за событие */
    }
    fc->o->lpacc[p] += fc->o->cdelta * d;
    if (fc->o->travel && fc->t_lv > 0.0) { /* §932-А: этаж по пустотному пробегу;
                                            * t_lv≤0 — нет предшественника
                                            * (прямой свет, А1691), инверсия
                                            * (А1702) или клеточная/path-ветка —
                                            * вклада нет автоматически */
      double lv = fc->t_lv;
      if (lv > 0.0) {
        double f = lv / fc->ldom;
        double gk;
        if (f > 1.0) f = 1.0;
        /* А1706+МЕЛОЧИ-4: гашение по ks/lep РЕБЁНКА. В rep-раздаче — просто
         * 1−ks (ks ≤ 1 по построению); в кусочной ветке — по КЛЭМПНУТОМУ
         * ks ветки депозита (ks ≤ 1−kdf). Без привязки к accum_mode. */
        if (fc->o->lep && fc->o->lep[p] > 0.0)
          gk = 0.0; /* §924-канон: эмиттеры не грубеют */
        else if (fc->o->ks && fc->o->ks[p] >= 0.0) {
          if (fc->ev_src_far) {
            gk = 1.0 - fc->o->ks[p];
          } else {
            double ks_c = fc->o->ks[p];
            double kdf = front_rho(fc, p);
            if (ks_c > 1.0 - kdf) ks_c = 1.0 - kdf > 0.0 ? 1.0 - kdf : 0.0;
            gk = 1.0 - ks_c;
          }
        } else
          gk = 1.0;
        g932_n++; /* А1724: статистика — при любом travel=1, ВНЕ guard'а вклада */
        g932_fsum += f;
        g932_lvsum += lv;
        g932_hist[f >= 1.0 ? 15 : (int)(f * 16.0)]++;
        if (fc->o->tvc > 0.0)
          fc->o->lpacc[p] += fc->o->tvc * f * gk; /* tvc — масштаб шага; tvc=0 —
                                                   * guard, битово без ключа (НК) */
      }
    }
  }
  if (fc->o->lphits) fc->o->lphits[p] += 1.0;    /* §862-диаг: ранжир (а) */
  if (fc->o->lpacc_e) fc->o->lpacc_e[p] += edep; /* §932-Б/А1737: Σдепозитов */
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
      Lh = front_lh_cap(fc, front_le(fc, pbest) +
                                front_rho(fc, pbest) * fc->Eprev[pbest] / (2.0 * M_PI));
      fc->recycled += fc->w_d * (Lh - fc->le) * csec;
    }
    fc->emitted += fc->w_d * front_le(fc, pbest) * csec;
    *a = 0.0;
    *b = front_lh_cap(fc,
                      front_le(fc, pbest) + front_rho(fc, pbest) * fc->Eprev[pbest] / (2.0 * M_PI));
    if (Lin < 0.0 && fc->negseen < 8) { /* §871: сентинел Lin<0 */
      fc->negseen++;
      fprintf(stderr, "NEG-A pbest=%d Lin=%.6g\n", pbest, Lin);
    }
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
        fc->Ed[p] += fc->w_d * Lin * f * csec * fc->axcos * wt[u] / swt / front_depden(fc, p);
        if (Lin * f < 0.0 && fc->negseen < 8) { /* §871: сентинел Lin<0 */
          fc->negseen++;
          fprintf(stderr, "NEG-T p=%d Lin=%.6g f=%.6g\n", p, Lin, f);
        }
        sw_accum(fc, p,
                 fc->w_d * Lin * f * csec * fc->axcos * wt[u] / swt /
                     front_depden(fc, p)); /* §862 */
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
static int64_t g924_hits, g924_lin;   /* §924-дых: rep-депозиты и ΣLin (прибор) */
static int64_t g924_vis, g924_geohit; /* §927-дых: визиты списков rep, гео-попадания */
static int64_t g924_kclamp;           /* §927: кулы k(r,ω) (диагностика) */
static int64_t g931_seam, g931_rh;    /* §931/А1668: шов — rep-хиты с точкой
                                       * пересечения ближе rep_zone; и все rep-хиты */
static int g924_minf = 1;             /* §924-матрика: HZ_REPMINF */

/* §932-Б/Ш12: макрос-«продолжение» цикла пер-хит событий для
 * КОНСЕРВАТИВНОГО КЛЕТОЧНОГО ПРИЁМА (cellc=2): листовая клетка набирает
 * Λ += A_p·|n_p·ω|/(csec·axcos) по кускам своего списка — это ТА ЖЕ Λ, что
 * знаменатель пер-хит депозита §894-e (там csec·axcos/A_p), только суммой по
 * клетке; извлекается Δ = Φ_вх·(1−e^{−Λ}) ≤ Φ_вх и раздаётся детям по
 * площадям (блок «КЛЕТОЧНЫЙ ПРИЁМ» в конце функции). Пер-хит ветка при
 * cellc=2 кусок НЕ обрабатывает (односчётность). `continue` внешнего цикла —
 * внутри макроса: `goto` вниз перескочил бы объявления ниже
 * (-Wjump-misses-init, ловля 07-10). */
#define SW_CELL_RECEIPT_CONTINUE()                                                                 \
  do {                                                                                             \
    if (fc->cellc == 2 && cpos >= 0) {                                                             \
      double cy = (fc->nrm != NULL) ? front_depden(fc, p) : fc->area[p];                           \
      double an2 = (fc->nrm != NULL) ? fabs(fc->om[0] * fc->nrm[3 * (int64_t)p] +                  \
                                            fc->om[1] * fc->nrm[3 * (int64_t)p + 1] +              \
                                            fc->om[2] * fc->nrm[3 * (int64_t)p + 2])               \
                                     : 1.0;                                                        \
      double c2 = (cy > FRONT_DEP_MIN_SHARE * csec && an2 > 1e-12)                                 \
                      ? (cy / an2) / (csec * fc->axcos)                                            \
                      : 1.0; /* вырожденная грань */                                               \
      fc->ced[cpos] += c2;                                                                         \
      g932c_vis++;                                                                                 \
      continue;                                                                                    \
    }                                                                                              \
  } while (0)

static void front_seg_walk(front_ctx *fc, const int32_t *ps, int32_t n, double tin, double tout,
                           double *a, double *b, int32_t cpos) {
  const hz_pyr *py = fc->py;
  double csec = py->cell * py->cell; /* сечение БАЗОВОЙ клетки — инвариант А1563 */
  double ai = *a, bi = *b, Lin = ai + bi;
  double *ht;
  int32_t *hp;
  int32_t u;
  int64_t nh = 0, i;
  double med_rem0 = fc->med_rem; /* Ш12: снимок — доля среды на ЭТОМ звене */
  if (n <= 0) return;
  if (fc->tau0) return; /* НК: слой не взаимодействует */
  fc->nmat++;
  g924_vis++; /* §927-дых: визит списка (любого) */
  /* §932-Б/Ш9-ЗОНД: КЛЕТОЧНЫЙ ПРИЁМ — накопитель на лист (Lin НА ВХОДЕ
   * сегмента; раздача детям по площадям в конце итерации). Зонд, НЕ
   * канон: сосуществует с пер-хит депозитами (двойной счёт на хитовых
   * кусках — так и было задумано в Ш9); умолчание cellc=0 — битово
   * прежний мир. cellc=2 — консервативная схема Ш12 (см. хвост функции). */
  if (fc->ced != NULL && cpos >= 0 && Lin > 0.0) /* cpos=-1: узел — вне зонда */
    fc->ced[cpos] += fc->w_d * Lin * csec * fc->axcos;
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
  for (u = 0; u < n; u++) {
    int32_t p = ps[u];
    /* §924: id ≥ nt0 — представитель (id уже РАСШИРЕННЫЙ tri id);
     * иначе — слот куска */
    int32_t tri = (fc->reps_on && p >= fc->nt0) ? p : py->pcs[p].tri;
    double tt, bh0, bh1;
    if (fc->med_on && !fc->in_leg && p >= fc->nt0 && fc->med_nsub[p - fc->nt0] > 0) {
      /* §932-Б: СРЕДА — вход в bbox региона (виртуальное событие), не
       * ray-tri носителя; повторные входы гасит med_stamp при активации */
      if (front_box_seg(fc->med_box + 6 * (int64_t)(p - fc->nt0),
                        fc->med_box + 6 * (int64_t)(p - fc->nt0) + 3, fc->org, fc->om, tin, tout,
                        &bh0, &bh1)) {
        ht[nh] = bh0;
        hp[nh] = p;
        nh++;
      }
      continue;
    }
    if (fc->agg &&
        fc->agg_list) { /* §872: серия (материал) — culling-единица:
                         * union bbox серии против [tin,tout] ОДНИМ slab-тестом; промах — вся
                         * серия мимо (строго консервативно: bbox члена ⊆ union), попадание —
                         * точная по-кусочная обработка ниже. Депозит/Lh — по-кусочно (§870:
                         * rep-ray-tri групп терял энергию — семантика снята). */
      int32_t mtl0 = py->pcs[p].mtl;
      int64_t g1 = u;
      double u0[3], u1[3];
      int ax;
      for (ax = 0; ax < 3; ax++) {
        u0[ax] = fc->o->tribox[6 * (int64_t)py->pcs[p].tri + ax];
        u1[ax] = fc->o->tribox[6 * (int64_t)py->pcs[p].tri + 3 + ax];
      }
      while (g1 < n && py->pcs[ps[g1]].mtl == mtl0) {
        const double *tb = fc->o->tribox + 6 * (int64_t)py->pcs[ps[g1]].tri;
        for (ax = 0; ax < 3; ax++) {
          if (tb[ax] < u0[ax]) u0[ax] = tb[ax];
          if (tb[3 + ax] > u1[ax]) u1[ax] = tb[3 + ax];
        }
        g1++;
      }
      if (!front_box_seg(u0, u1, fc->org, fc->om, tin, tout, &bh0, &bh1)) {
        u = (int32_t)g1 - 1; /* for-u++ перепрыгнет всю серию */
        continue;
      }
    }

    if (fc->pstamp[tri] == fc->pkey) {
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
    { /* §894-e: пер-хит депозит для ЛЮБОЙ трубки — схема согласована с
       * MC-эталоном (коридор 1.0023, одно зеркало 1.0365); линейная
       * нагрузка §892 снята: релей на ногах множил энергию (1.46) */
      int first = 1;
      for (i = 0; i < nh; i++) {
        int32_t p = hp[i];
        int32_t tri = (fc->reps_on && p >= fc->nt0) ? p : py->pcs[p].tri; /* §924 */
        double Lh;
        if (fc->pstamp[tri] == fc->pkey) continue;
        fc->pstamp[tri] = fc->pkey;
        SW_CELL_RECEIPT_CONTINUE(); /* Ш12: клеточный приём — ДО пер-хит ветки */
        if (fc->o->travel) {        /* §932-А: проводка t хита, пробега Lv и
                                     * монотонного t_last — СРАЗУ после штампа, ДО
                                     * обоих continue ниже (прозрачный rep — тоже
                                     * материальное событие, А1703); max гасит
                                     * инверсии уровней А1702. Lv — по СТАРОМУ t_last
                                     * (пробег ДО этого хита) */
          fc->t_cur = ht[i];
          fc->t_lv = fc->t_last >= 0.0 ? ht[i] - fc->t_last : -1.0;
          fc->t_last = g932_nomax ? ht[i] : (fc->t_last > ht[i] ? fc->t_last : ht[i]);
        }
        if (fc->med_on && !fc->in_leg) { /* §932-Б: ослабление средой до
                                          * события (А1710), затем — сам
                                          * случай: вход или кусковый хит */
          int32_t rr = (fc->reps_on && p >= fc->nt0) ? p - fc->nt0 : -1;
          Lin = sw_med_bank(fc, Lin, ht[i], csec);
          if (rr >= 0 && fc->med_nsub[rr] > 0) { /* ВИРТУАЛЬНЫЙ ВХОД (А1695) */
            sw_med_enter(fc, rr, ht[i]);
            if (first) {
              fc->depA += ai;
              first = 0;
            }
            continue; /* носитель-поверхность не взаимодействует */
          }
        }
        if (first) {
          fc->depA += ai;
          first = 0;
        }
        if (fc->reps_on && p >= fc->nt0) { /* §924: ПРЕДСТАВИТЕЛЬ — депозит
                                            * с раздачей детям по площадям;
                                            * штамп гасит и вложенных rep'ов
                                            * (rep_f ⊇ rep_{f+1}), и листья */
          int32_t r = p - fc->nt0;
          fc->ev_src_far = 1; /* §931/А1666: событие-источник — дальнее */
          {                   /* §931/А1668: ШОВ — rep-хит с ТОЧКОЙ пересечения ближе зоны
                               * (грубый носитель заходит в ближний шар); порог тревоги 5%
                               * от rep-хитов — тогда усиливать предикат bbox-тестом. */
            double px = fc->org[0] + fc->om[0] * ht[i] - fc->o->rep_eye[0],
                   pyy = fc->org[1] + fc->om[1] * ht[i] - fc->o->rep_eye[1],
                   pz = fc->org[2] + fc->om[2] * ht[i] - fc->o->rep_eye[2];
            g931_rh++;
            if (fc->o->rep_zone > 0.0 &&
                px * px + pyy * pyy + pz * pz <= fc->o->rep_zone * fc->o->rep_zone)
              g931_seam++;
          }
          {
            static long kcpos = 0, kcneg = 0;
            if (fc->o->kcal) {
              double kd9 = fc->o->kcal[(size_t)r * HZ_KCAL_NW + (size_t)fc->dir_bin];
              if (kd9 > 0)
                kcpos++;
              else
                kcneg++;
              if ((kcpos + kcneg) % 50000 == 0)
                fprintf(stderr, "KCDbg pos=%ld neg=%ld\n", kcpos, kcneg);
            }
          }
          if (fc->o->kcal != NULL) { /* §928: КАЛИБРОВАННАЯ таблица k(r,ωbin);
                                      * ячейка ≤ 0 — fallback на формулу §927 ниже.
                                      * extraction = k·w·csec·axcos·(Lin+Lh) — база
                                      * симметрична сбору (тёмная трубка зажигается
                                      * эмиссией самой цепочки); раздача детям по
                                      * долям СРЕДИ СВЕЖИХ: Σ депозитов = flux. */
            double kc = fc->o->kcal[(size_t)r * HZ_KCAL_NW + (size_t)fc->dir_bin];
            if (kc > 0.0) {
              double Lh9c =
                  front_lh_cap(fc, fc->rep_le[r] + fc->rep_rho[r] * fc->rep_ep[r] / (2.0 * M_PI));
              double freshA = 0.0; /* §928: свежие дети (вложенные репы делят детей —
                                    * ренормировка вместо потери доли flux) */
              int32_t nfr = 0;
              for (int32_t kq = fc->rep_koff[r]; kq < fc->rep_koff[r + 1]; kq++) {
                int32_t kid = fc->rep_kmem[kq];
                if (!sw_far(fc->o, kid))
                  continue; /* §931/А1662: ближний кусок
                             * НЕ получает rep-долю — у него
                             * свои кусочные события */
                if (fc->pstamp[py->pcs[kid].tri] == fc->pkey) continue;
                if (fc->o->lep && fc->o->lep[kid] > 0.0) continue;
                freshA += fc->area[kid];
                nfr++;
              }
              if (nfr > 0 && freshA > 0.0) {
                double flux = fc->w_d * csec * fc->axcos * (Lin + Lh9c) * kc;
                for (int32_t kq = fc->rep_koff[r]; kq < fc->rep_koff[r + 1]; kq++) {
                  int32_t kid = fc->rep_kmem[kq];
                  int32_t ktri = py->pcs[kid].tri;
                  if (!sw_far(fc->o, kid)) continue;          /* §931/А1662 */
                  if (fc->pstamp[ktri] == fc->pkey) continue; /* уже получил долю */
                  if (fc->o->lep && fc->o->lep[kid] > 0.0)
                    continue; /* §924: эмиттер — своё событие в листе */
                  fc->pstamp[ktri] = fc->pkey;
                  double ed = flux * (fc->area[kid] / freshA); /* Σ долей = 1 */
                  fc->Ed[kid] += ed;
                  sw_accum(fc, kid, ed); /* этажная логика §918 — на ребёнке */
                }
                fc->absorbed += flux; /* §928: = Σ депозитов (ренормировка) */
                fc->emitted += fc->w_d * fc->rep_le[r] * csec;
                fc->recycled += fc->w_d * (Lh9c - fc->le) * csec;
                fc->ndep++;
                g924_hits++; /* §928: как формульная ветка — честный rep_hits */
                g924_lin += (int64_t)((Lin + Lh9c) * 1000.0);
              }
              /* нет свежих детей (их забрал внешний реп) — ПРОЗРАЧНЫЙ: без
               * извлечения (иначе двойной счёт одной цепочки), физическое
               * продолжение Lin ← Lh как у любой поверхности */
              Lin = Lh9c;
              continue;
            }
          }
          /* §927: НАПРАВЛЕННЫЙ ПЕРЕХВАТ k(r,ω) = ΣДетей A_i|n_i·ω| /
           * (A_носителя·|n_r·ω|) — ожидаемые пересечения грозди на
           * пересечение носителя. Считается тем же циклом, что раздача.
           * gain петли = 1 по построению: пол (копланарные дети) → k=1,
           * колонна → k≈2 поперёк, →0 вдоль оси. Кул [0,64] — численная
           * страховка вырожденных граней (счётчик ниже). */
          double cosr = fabs(fc->om[0] * fc->rep_nrm[3 * (int64_t)r] +
                             fc->om[1] * fc->rep_nrm[3 * (int64_t)r + 1] +
                             fc->om[2] * fc->rep_nrm[3 * (int64_t)r + 2]); /* смещение 64-битное */
          double knum = 0.0;
          int32_t k;
          for (k = fc->rep_koff[r]; k < fc->rep_koff[r + 1]; k++) {
            int32_t kid = fc->rep_kmem[k];
            if (!sw_far(fc->o, kid))
              continue; /* §931/А1676: kk — только по
                         * дальним детям, ДО ffull/frest */
            const double *nk = fc->nrm + 3 * (int64_t)kid;
            knum += fc->area[kid] * fabs(fc->om[0] * nk[0] + fc->om[1] * nk[1] + fc->om[2] * nk[2]);
          }
          double kk =
              (cosr > 1e-12 && fc->rep_atri[r] > 1e-30) ? knum / (fc->rep_atri[r] * cosr) : 1.0;
          if (kk > 64.0) {
            kk = 64.0;
            g924_kclamp++;
          }
          if (kk < 1e-6) kk = 1e-6;
          /* §927-2: ЦЕПОЧКА ПЕРЕСЕЧЕНИЙ — энергосохранение строго:
           * гроздь обменивалась бы k раз: 1-е пересечение снимает
           * Φ(Lin), остальные k−1 — по Φ(Lh) (переизлучённое), дробный
           * хвост — долей. Σ депозитов = Σ извлечённого из трубки
           * (ловля §927: депозит kΦ при одном изъятии Φ = источник
           * энергии → разгон ×1.55/ит). Пол (k=1) вырождается в
           * однократное пересечение — прежняя арифметика. */
          double Lh9 =
              front_lh_cap(fc, fc->rep_le[r] + fc->rep_rho[r] * fc->rep_ep[r] / (2.0 * M_PI));
          double ffull = kk < 1.0 ? kk : 1.0;
          double frest = kk > 1.0 ? kk - 1.0 : 0.0;
          double fluxsum = fc->w_d * csec * fc->axcos * (ffull * Lin + frest * Lh9);
          double edep = fluxsum / fc->rep_area[r];
          for (int32_t kq = fc->rep_koff[r]; kq < fc->rep_koff[r + 1]; kq++) {
            int32_t kid = fc->rep_kmem[kq];
            int32_t ktri = py->pcs[kid].tri;
            if (!sw_far(fc->o, kid)) continue;          /* §931/А1662/А1676 */
            if (fc->pstamp[ktri] == fc->pkey) continue; /* уже получил долю */
            if (fc->o->lep && fc->o->lep[kid] > 0.0)
              continue; /* §924: эмиттер — своё событие в листе, штамп не
                         * вешаем: иначе его листовое событие умрёт */
            fc->pstamp[ktri] = fc->pkey;
            fc->Ed[kid] += edep;
            sw_accum(fc, kid, edep); /* этажная логика §918 — на ребёнке */
          }
          Lh = Lh9;
          fc->absorbed += fc->w_d * csec * (ffull * Lin + frest * Lh9); /* §927-2:
                                                                         * = Σ депозитов */
          fc->emitted += fc->w_d * fc->rep_le[r] * csec;
          fc->recycled += fc->w_d * (Lh - fc->le) * csec;
          fc->ndep++;
          g924_hits++;
          g924_lin += (int64_t)(Lin * 1000.0);
          {
            static int repdbg_n = 0; /* §924-дых */
            static int repdbg_on = -1;
            if (repdbg_on < 0) repdbg_on = getenv("HZ_DBG924") != NULL;
            if (repdbg_on && repdbg_n < 6) {
              int32_t nk = 0;
              for (int32_t kq = fc->rep_koff[r]; kq < fc->rep_koff[r + 1]; kq++)
                if (fc->pstamp[py->pcs[fc->rep_kmem[kq]].tri] == fc->pkey) nk++;
              fprintf(stderr,
                      "DBG924hit r=%d Lin=%.4g area=%.4g rho=%.3g le=%.3g ep=%.4g kids=%d fresh=%d "
                      "edep=%.4g Lh=%.4g\n",
                      r, Lin, fc->rep_area[r], fc->rep_rho[r], fc->rep_le[r], fc->rep_ep[r],
                      fc->rep_koff[r + 1] - fc->rep_koff[r], nk, edep, Lh);
              repdbg_n++;
            }
          }
          Lin = Lh;
          continue;
        }
        double ks = fc->o->ks ? fc->o->ks[p] : 0.0;
        fc->ev_src_far = 0; /* §931/А1666: событие-источник — кусковое */
        double kdf = front_rho(fc, p);
        if (ks > 1.0 - kdf) ks = 1.0 - kdf > 0.0 ? 1.0 - kdf : 0.0;
        double lin_hop = ks * Lin;
        int relay_ok = lin_hop > 0.0 && lin_hop >= 1e-3 * fc->leg_root; /* §894-c-5 */
        if (relay_ok && fc->hop_n < 32) {
          const double *nv2 = fc->nrm + 3 * (int64_t)p;
          double dot2 = fc->om[0] * nv2[0] + fc->om[1] * nv2[1] + fc->om[2] * nv2[2];
          double tth = ht[i];
          int q2;
          for (q2 = 0; q2 < 3; q2++) {
            fc->hop_pt[fc->hop_n][q2] = fc->org[q2] + fc->om[q2] * tth;
            fc->hop_dir[fc->hop_n][q2] = fc->om[q2] - 2.0 * dot2 * nv2[q2];
          }
          fc->hop_lin[fc->hop_n] = lin_hop;
          fc->hop_root[fc->hop_n] = fc->leg_root;
          fc->hop_n++;
        } else if (!relay_ok) {
          fc->hop_lost += lin_hop; /* порог — доля в lost (баланс цел) */
        }
        /* §911-3: перехват куском не превышает поток трубки — клэмп
         * csec·μ/depden до 1 (клип-кусок area≪csec прежде давал
         * dep/Lin до ~87 → самоподкачка → взрыв клипа). Вариант
         * «покрытие» μ·min(1,csec/area) откатан: он меняет депозиты
         * мелких кусков не-клипа (шары) без доказательства. Клип-
         * нестабильность (56→3.8e9→1.2e18, и без лампы 823 на it=2)
         * — усиление обратной связи полосок, владелец §911-4. */
        double dep_ratio = csec * fc->axcos / front_depden(fc, p);
        if (dep_ratio > 1.0) dep_ratio = 1.0;
        if (fc->cov) { /* §911-6: Σ долей следа против леммы покрытия */
          fc->cov[p] += dep_ratio;
          if (Lin > 0.0) fc->cov_lit[p] += dep_ratio; /* §911-8: с светом цепи */
          fc->cov_e[p] += dep_ratio * Lin;            /* §911-10: амплитуда */
          fc->covn[p]++;
        }
        double dep = fc->w_d * Lin * (1.0 - ks) * dep_ratio;
        {
          /* §911-2: прибор депозит/входящий поток — печать первых 20
           * фактов (piece, area, csec, axcos, dep/Lin) при HZ_DEPDBG=1;
           * здоровый прогон бесплатен */
          static int depdbg_init = 0;
          static int depdbg_on = 0;
          static int depdbg_n = 0;
          if (!depdbg_init) {
            const char *e = getenv("HZ_DEPDBG");
            depdbg_on = e && e[0] != '\0' && e[0] != '0';
            depdbg_init = 1;
          }
          if (depdbg_on && depdbg_n < 20 && Lin > 0.0 && dep > 2.0 * Lin) {
            double area_p = fc->area[p];
#pragma omp critical(depdbg)
            {
              fprintf(stderr, "DEPDBG[%d] piece=%d area=%.3g csec=%.3g axcos=%.3g dep/Lin=%.3g\n",
                      depdbg_n, p, area_p, csec, fc->axcos,
                      (csec * fc->axcos) / front_depden(fc, p));
              depdbg_n++;
            }
          }
        }
        if (fc->in_leg) {
          fc->dep_leg += dep;
        } else {
          fc->dep_main += dep;
          fc->dep_cnt++;
        }
        if (fc->o->strip_start) {
          /* §911-13-3: депозит tri распределяется на ВСЕ полоски tri —
           * E полосок однородна (гранулярность поля), Σ a_i·E_i =
           * μ·L·a_tri (консервативность) */
          int32_t t3 = py->pcs[p].tri;
          int32_t s0 = fc->o->strip_start[t3], s1 = fc->o->strip_start[t3 + 1], si;
          for (si = s0; si < s1; si++) {
            int32_t q = fc->o->strip_list[si];
            fc->Ed[q] += dep; /* §911-13-3: ровно dep (с axcos) каждой полоске */
          }
        } else {
          fc->Ed[p] += dep;
        }
        fc->absorbed += fc->w_d * Lin * (1.0 - ks) * csec;
        fc->emitted += fc->w_d * front_le(fc, p) * csec;
        Lh = front_lh_cap(fc, front_le(fc, p) + kdf * fc->Eprev[p] / (2.0 * M_PI));
        fc->recycled += fc->w_d * (Lh - fc->le) * csec;
        fc->ndep++;
        g924_hits++; /* §924-дых: событие куска (мир §923) */
        g924_lin += (int64_t)(Lin * 1000.0);
        if (fc->knum && !fc->in_leg &&
            !fc->o->strip_start) { /* §928: калибровка —
                                    * депозиты серийной цепочки региона, приписанные ОДНОМУ уровню:
                                    * этажу куска (want-семантика §924) — доставка в reps-мире
                                    * происходит ТОЛЬКО на уровне хита, считая все уровни, мы бы
                                    * завысили знаменатель и занизили k */
          int32_t fl = fc->kfloor != NULL ? (int32_t)fc->kfloor[p] : 1;
          if (fl < 1) fl = 1;
          if (fl > fc->o->rep_nlev) fl = fc->o->rep_nlev;
          int32_t r = fc->o->repof[((size_t)fl - 1u) * (size_t)fc->nt0 + (size_t)p];
          if (r >= 0) {
            if (fc->pstamp[fc->nt0 + r] != fc->pkey) { /* вход трубки в регион:
                                                        * зеркалирует штамп rep-хита.
                                                        * БАЗА (Lin+Lh): тёмная трубка
                                                        * зажигается эмиссией САМОЙ
                                                        * цепочки (le+ρEprev/2π на
                                                        * каждой поверхности) — извлечение
                                                        * на входе обязано её видеть */
              fc->pstamp[fc->nt0 + r] = fc->pkey;
              fc->kden[(size_t)r * HZ_KCAL_NW + (size_t)fc->dir_bin] +=
                  fc->w_d * csec * fc->axcos * (Lin + Lh);
            }
            fc->knum[(size_t)r * HZ_KCAL_NW + (size_t)fc->dir_bin] += dep;
          }
        }
        sw_accum(fc, p, fc->w_d * Lin * (1.0 - ks) * dep_ratio);
        Lin = Lh;
      }
      g924_geohit += nh;             /* §927-дых: гео-попадания (все списки) */
      if (fc->med_on && !fc->in_leg) /* §932-Б: банк до ГРАНИЦЫ сегмента —
                                      * состояние переживает сегмент (А1719),
                                      * инкрементность точна, двойного счёта нет */
        Lin = sw_med_bank(fc, Lin, tout, csec);
      /* Ш12: КЛЕТОЧНЫЙ ПРИЁМ (cellc=2) — извлечение из трубки и выход из
       * клетки. Λ собрана в ced[cpos] (см. макрос в цикле пер-хит событий):
       * Δ = Φ_вх·(1−e^{−Λ}) ≤ Φ_вх (консервативность по построению),
       * Δ/ΣA_cell раздаётся детям в конце итерации (Σ Ed·A = Δ).
       * Трубка через клетку: Lin ← Lin·(1−f) + L̄_cell·f, L̄_cell — средний
       * по площадям Lh кусков списка (клетка как поверхность; для клетки из
       * одного куска это пер-хит с точностью O(f²)). Носитель-поверхность
       * (виртуальные входы среды) в списке листа не встречается. */
      if (fc->cellc == 2 && cpos >= 0 && (fc->ced[cpos] > 0.0 || fc->med_rem > med_rem0)) {
        double Lin_entry = Lin; /* радианс на входе клетки — для доли
                                 * переизлучения клетки (cabs) */
        /* Λ = перехват кусками (ced) + поглощение СРЕДОЙ на звене (разность
         * снимков банка): оба слагаемых вычитаются из ОДНОЙ трубки, поэтому
         * складываются ДО экспоненты */
        double cyc = fc->ced[cpos] + (fc->med_rem - med_rem0);
        double f = -expm1(-cyc);
        double Lin_cell = Lin;
        if (f > 1.0 || !(f == f)) { /* страховка: f ≤ 1 по построению (expm1 > −1);
                                     * срабатывание = сигнал NaN/Inf в Λ */
          static int cfw = 0;
          if (cfw < 8) {
            cfw++;
            fprintf(stderr, "Ш12-CLAMP: f=%.17g cyc=%.17g ced=%.17g med=%.17g\n", f, cyc,
                    fc->ced[cpos], fc->med_rem - med_rem0);
          }
          f = 1.0;
        }
        if (f > g932c_fmax) g932c_fmax = f;
        if (f >= 1.0 - 1e-12) g932c_ge1++; /* прибор: исчерпание клетки (f→1) */
        g932c_vis++;
        if (Lin > 0.0) { /* L̄_cell — средний по площадям Lh кусков списка */
          double lsum = 0.0, asum = 0.0;
          for (int32_t q3 = fc->bstart[cpos]; q3 < fc->bstart[cpos + 1]; q3++) {
            int32_t q4 = fc->bpids[q3];
            double lh4 = front_lh_cap(fc, front_le(fc, q4) +
                                              front_rho(fc, q4) * fc->Eprev[q4] / (2.0 * M_PI));
            lsum += lh4 * fc->area[q4];
            asum += fc->area[q4];
          }
          if (asum > 0.0) Lin_cell = Lin * (1.0 - f) + (lsum / asum) * f;
        }
        fc->ced[cpos] = fc->w_d * Lin_entry * csec * fc->axcos * f;
        fc->cabs[cpos] = Lin_entry;
        fc->absorbed += fc->w_d * Lin_entry * csec * fc->axcos * f;
        fc->depA += ai;
        if (Lin_cell > 0.0)
          *b = Lin_cell;
        else
          *b = 0.0;
        *a = 0.0;
        Lin = Lin_cell;
      }
      if (nh > 0) {
        *a = 0.0;
        *b = Lin;
      }
      return;
    }
    /* §892: ЛИНЕЙНАЯ НАГРУЗКА КЛЕТКИ (T4) — поток трубки делится между
     * ВСЕМИ кусками клетки пропорционально площади: Ed_p +=
     * w_d·Lin·csec·axcos·(1−ks_p)/Σarea. Нет лотереи треугольников
     * (ray effect на уровне кусков) и 1/area-сингулярности (знаменатель
     * Σarea ограничен). Зеркальный релей — хоп от ближайшего точного
     * удара (pbest по списку ray-tri: точка, R, ks·Lin). */
    int32_t pbest = -1;
    double tb = 1e30;
    for (i = 0; i < nh; i++)
      if (ht[i] < tb) {
        tb = ht[i];
        pbest = hp[i];
      }
    /* §894-b ОТКАТ: хоп на каждый зеркальный кусок клетки ДАВАЛ пере-
     * светление ×1.5-2 (коридор 1.47, одно зеркало 1.97): хоп-нога при
     * линейной нагрузке депонирует свой поток в КАЖДОЙ клетке пути без
     * декремента — N клеток = N× энергия. Возврат к pbest-хопу до
     * проектирования корректной семантики ног (§894-c). */
    if (pbest >= 0) {
      double ks = fc->o->ks ? fc->o->ks[pbest] : 0.0;
      double kdf = front_rho(fc, pbest);
      if (ks > 1.0 - kdf) ks = 1.0 - kdf > 0.0 ? 1.0 - kdf : 0.0;
      if (ks > 0.0 && Lin > 0.0 && fc->hop_n < 32) {
        const double *nv = fc->nrm + 3 * (int64_t)pbest;
        double dot = fc->om[0] * nv[0] + fc->om[1] * nv[1] + fc->om[2] * nv[2];
        double tth = tb;
        int q2;
        for (q2 = 0; q2 < 3; q2++) {
          fc->hop_pt[fc->hop_n][q2] = fc->org[q2] + fc->om[q2] * tth;
          fc->hop_dir[fc->hop_n][q2] = fc->om[q2] - 2.0 * dot * nv[q2];
        }
        fc->hop_lin[fc->hop_n] = ks * Lin;
        fc->hop_root[fc->hop_n] = Lin;
        fc->hop_n++;
      } else if (ks > 0.0 && Lin > 0.0) {
        fc->hop_lost += ks * Lin;
      }
      {
        double lhcap_p = front_le(fc, pbest);
        if (fc->Linmax_prev[pbest] > lhcap_p) lhcap_p = fc->Linmax_prev[pbest];
        double lhraw = front_le(fc, pbest) + kdf * fc->Eprev[pbest] / (2.0 * M_PI);
        fc->Lh_last = front_lh_cap(fc, lhraw > lhcap_p ? lhcap_p : lhraw);
      }
      fc->recycled += fc->w_d * (fc->Lh_last - fc->le) * csec;
      fc->ndep++;
    }
    if (nh > 0) {
      /* §892-b: зажигание БЕЗ условия Lin>0 — тёмная трубка подхватывает
       * эмиссию поверхности (Lh = le+lep > 0); иначе она тёмная навсегда
       * (баг нулевых депозитов §892) */
      *a = 0.0;
      *b = pbest >= 0 ? fc->Lh_last : Lin;
    }
  }
}
#undef SW_CELL_RECEIPT_CONTINUE

/* §923: пересборка раскладки под текущие этажи lpflo (слот-пространство).
 * Каждый кусок — в ЕДИНСТВЕННЫЙ список уровня своего этажа: f==0 — листовая
 * клетка, f≥1 — узел-предок уровня f-1 (поиск hz_pyr_node_pos; узел-предок
 * существует гарантированно — под ним есть его кусок; отказ → лист, fail
 * closed). finer[l][pos] = «в поддереве есть куски этажа ≤ l» —
 * распространяется снизу вверх, один поиск родителя на узел. O(nt) счёт +
 * O(поиски). Вызывается раз в итерацию, после применения этажей. */
static uint8_t g_dbg923_vis[1 << 20]; /* §923-дых: отметки посещённых листьев (только HZ_DBG923) */

/* §923: клетка bbox → позиция листа; листа нет (фантомная колонка
 * max-грани, А1578) — фолбэк соседей ±1 (как fb в sw_index_piece,
 * только walk). *dst[] — до 7 позиций (центр + 6 соседей). */
static int32_t sw923_leaf_or_nb(const hz_pyr *py, int64_t ix, int64_t iy, int64_t iz, int fb,
                                int32_t dst[7]) {
  int32_t n = 0, k;
  int32_t pos = hz_pyr_leaf_pos(py, ix + py->nx * (iy + py->ny * iz));
  int64_t nb[3];
  if (pos >= 0) dst[n++] = pos;
  if (!fb) return n;
  nb[0] = ix;
  nb[1] = iy;
  nb[2] = iz;
  for (int ax = 0; ax < 3; ax++) {
    int64_t mx = nb[ax];
    for (int sgn = -1; sgn <= 1; sgn += 2) {
      int64_t lim = ax == 0 ? py->nx : (ax == 1 ? py->ny : py->nz);
      nb[ax] = mx + sgn;
      if (nb[ax] >= 0 && nb[ax] < lim) {
        pos = hz_pyr_leaf_pos(py, nb[0] + py->nx * (nb[1] + py->ny * nb[2]));
        if (pos >= 0) {
          int dup = 0;
          for (k = 0; k < n; k++)
            if (dst[k] == pos) dup = 1;
          if (!dup) dst[n++] = pos;
        }
      }
    }
    nb[ax] = mx;
  }
  return n;
}

/* §923 (v3): пересборка раскладки под текущие этажи lpflo.
 * Кусок этажа f кладётся во ВСЕ узлы уровня f-1 (f≥1) / базовые клетки
 * (f==0), которые пересекает его bbox — как А1566 для bpids (владение
 * клеткой теряет торчащие треугольники, ловля §922-дым: E 0.0112).
 * finer[l][pos] = «в bbox-перекрытии узла есть куски этажа ≤ l».
 * Итоговые числа записей известны только после счёта — пулы растятся
 * здесь (pool_cap/leaf_cap — ёмкости владельца); отказ аллокации →
 * возврат -1, вызывающий разбирает lvls и возвращается в старый мир. */
static int sw_levels_build(const hz_pyr *py, const hz_sw_opts *o, const uint8_t *lpflo, int32_t nt,
                           sw_lvl *lvls, int32_t **pool_io, int32_t *pool_cap, int32_t **leaf_io,
                           int32_t *leaf_cap, int32_t *leaf_off, const int32_t *bstart,
                           const int32_t *bpids, int32_t *want, const uint8_t *med_sub) {
  int32_t l, p, u;
  int64_t (*dims)[3];
  static int dbg = -1;
  /* §928: reps_collect — карты построены для СБОРА калибровки, но списки
   * уровней остаются кусочными (мир §923) */
  /* want != NULL в условии — НЕ украшение: при отказе calloc'а rep_want
   * вызывающий обнуляет указатель и ЖДЁТ возврата в мир без представителей
   * («fail closed», см. hz_sw_run §924). Без этой проверки ниже (2121, 2245)
   * было разыменование NULL — clang-analyzer-core.NullDereference, ловля 07-10. */
  const int reps_on =
      (want != NULL && o->nrep > 0 && o->walk && o->repof != NULL && !o->reps_collect);
  if (reps_on && g924_minf == 1) { /* §924-матрика: HZ_REPMINF (гибрид (в)) */
    const char *e = getenv("HZ_REPMINF");
    if (e) {
      int v = atoi(e);
      if (v >= 1) g924_minf = v;
    }
  }
  if (py->nlev <= 0) return 0;
  if (dbg < 0) dbg = getenv("HZ_DBG923") != NULL;
  dims = (int64_t (*)[3])malloc((size_t)py->nlev * sizeof *dims);
  if (!dims) return -1;
  for (l = 0; l < py->nlev; l++) {
    hz_pyr_level_dims(py, l, dims[l]);
    for (u = 0; u <= lvls[l].n; u++)
      lvls[l].off[u] = 0;
    memset(lvls[l].finer, 0, (size_t)lvls[l].n);
  }
  for (u = 0; u <= py->nleaf; u++)
    leaf_off[u] = 0;
  if (reps_on && want != NULL) /* §924: метки прошлой итерации гасим */
    memset(want, 0, (size_t)o->nrep * sizeof *want);
  /* 1. счёт + finer */
  for (p = 0; p < nt; p++) {
    uint8_t f = lpflo[p];
    int32_t tri = py->pcs[p].tri;
    int64_t lo[3], hi[3], cspan[6];
    if (f > (uint8_t)py->nlev) f = 0; /* клэмп отказ-на-лист (А1568) */
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
    { /* консервативный спэн ЗААНЯТОСТИ пирамиды (floor–floor, клэмп; так
       * строились листья — pyr_bbox_span): им и ТОЛЬКО им метится finer,
       * иначе листья, существующие только по широкому спэну, становятся
       * недостижимыми (ловля §923-дым: 980 листьев без визитов, E −5.5%) */
      for (int ax = 0; ax < 3; ax++) {
        double cmin = o->tribox[6 * (int64_t)tri + ax];
        double cmax = o->tribox[6 * (int64_t)tri + 3 + ax];
        int64_t lim = ax == 0 ? py->nx : (ax == 1 ? py->ny : py->nz);
        int64_t a0 = (int64_t)((cmin - py->lo[ax]) / py->cell);
        int64_t a1 = (int64_t)((cmax - py->lo[ax]) / py->cell);
        if (a0 < 0) a0 = 0;
        if (a0 > lim - 1) a0 = lim - 1;
        if (a1 < a0) a1 = a0;
        if (a1 > lim - 1) a1 = lim - 1;
        /* ax ∈ [0,3) по циклу выше — переполнения нет по построению; приведение
         * к int64_t только чтобы класс widening не шумел на 3-элементном массиве */
        cspan[(int64_t)ax * 2] = a0;
        cspan[(int64_t)ax * 2 + 1] = a1;
      }
    }
    if (reps_on && (int32_t)f >= g924_minf && (int32_t)f <= o->rep_nlev &&
        o->repof[((size_t)(int32_t)f - 1u) * (size_t)nt + (size_t)p] >= 0 && sw_far(o, p) &&
        (med_sub == NULL || med_sub[p])) { /* §932-Б: среда ЗАМЕЩАЕТ reps для
                                            * энергетически малых кусков (А1737);
                                            * med_sub==NULL — прежний §924-мир */
      /* §924: кусок заменён представителем (want[r]=f — дедуп на гроздь);
       * счёт узлов и finer — проходы B/C ниже. §931: БЛИЖНИЙ кусок не
       * заменяется (кольцо детальности — кусочный мир §923, битово). */
      want[o->repof[((size_t)(int32_t)f - 1u) * (size_t)nt + (size_t)p]] = (int32_t)f;
      goto rep_kid;
    }
    for (l = 0; l < py->nlev; l++) {
      if (f >= 1 && l + 1 < (int32_t)f) continue; /* глубже своего уровня куска нет */
      if (f >= 1 && l + 1 == (int32_t)f) {        /* СВОЙ уровень: счёт узлов-предков */
        int64_t sh = (int64_t)l + 1;
        for (int64_t az = lo[2] >> sh; az <= hi[2] >> sh; az++)
          for (int64_t ay = lo[1] >> sh; ay <= hi[1] >> sh; ay++)
            for (int64_t axx = lo[0] >> sh; axx <= hi[0] >> sh; axx++) {
              int64_t nid = axx + dims[l][0] * (ay + dims[l][1] * az);
              int32_t pos = hz_pyr_node_pos(py, l, nid);
              if (pos >= 0) {
                lvls[l].off[pos + 1]++;
              } else { /* узла нет — весь bbox в лист (fail closed, fb А1578) */
                int fb = (o->walk != 0);
                for (int64_t iz = lo[2]; iz <= hi[2]; iz++)
                  for (int64_t iy = lo[1]; iy <= hi[1]; iy++)
                    for (int64_t ix = lo[0]; ix <= hi[0]; ix++) {
                      int32_t dst[7];
                      int32_t nn = sw923_leaf_or_nb(py, ix, iy, iz, fb, dst);
                      for (int32_t k = 0; k < nn; k++)
                        leaf_off[dst[k] + 1]++;
                    }
              }
            }
      } else if (f == 0 && l == 0) {
        /* floor-0: во все пересекаемые базовые клетки (как bpids, fb А1578) */
        int fb = (o->walk != 0);
        for (int64_t iz = lo[2]; iz <= hi[2]; iz++)
          for (int64_t iy = lo[1]; iy <= hi[1]; iy++)
            for (int64_t ix = lo[0]; ix <= hi[0]; ix++) {
              int32_t dst[7];
              int32_t nn = sw923_leaf_or_nb(py, ix, iy, iz, fb, dst);
              for (int32_t k = 0; k < nn; k++)
                leaf_off[dst[k] + 1]++;
            }
      }
      if (f == 0 || l >= (int32_t)f) { /* глубже уровня l есть кусок — finer
                                        * (по консервативному спэну занятости) */
        int64_t sh = (int64_t)l + 1;
        for (int64_t az = cspan[4] >> sh; az <= cspan[5] >> sh; az++)
          for (int64_t ay = cspan[2] >> sh; ay <= cspan[3] >> sh; ay++)
            for (int64_t axx = cspan[0] >> sh; axx <= cspan[1] >> sh; axx++) {
              int64_t nid = axx + dims[l][0] * (ay + dims[l][1] * az);
              int32_t pos = hz_pyr_node_pos(py, l, nid);
              if (pos >= 0) lvls[l].finer[pos] = 1;
            }
      }
    }
    continue;
  rep_kid:; /* §924: finer rep-куска по СВОЕМУ консервативному спэну —
             * спэн rep'а (проход C) может не накрыть торчащего ребёнка,
             * а предки ребёнка обязаны звать спуск */
    for (l = (int32_t)f - 1; l < py->nlev; l++) {
      int64_t sh = (int64_t)l + 1;
      for (int64_t az = cspan[4] >> sh; az <= cspan[5] >> sh; az++)
        for (int64_t ay = cspan[2] >> sh; ay <= cspan[3] >> sh; ay++)
          for (int64_t axx = cspan[0] >> sh; axx <= cspan[1] >> sh; axx++) {
            int64_t nid = axx + dims[l][0] * (ay + dims[l][1] * az);
            int32_t pos = hz_pyr_node_pos(py, l, nid);
            if (pos >= 0) lvls[l].finer[pos] = 1;
          }
    }
  }
  if (reps_on) { /* §924, проход B: счёт want-представителей в узлы их
                  * спэна (А1578-формула по tribox rep'а, id = nt + r) */
    for (int32_t r = 0; r < o->nrep; r++) {
      if (want[r] <= 0) continue;
      int32_t li = want[r] - 1;
      int64_t rlo[3], rhi[3];
      int64_t triid = (int64_t)nt + r;
      int bad = 0;
      for (int ax = 0; ax < 3; ax++) {
        double cmin = o->tribox[6 * triid + ax];
        double cmax = o->tribox[6 * triid + 3 + ax];
        double dom = o->domhi[ax] - py->lo[ax];
        rlo[ax] = (int64_t)((cmin - py->lo[ax]) / py->cell);
        rhi[ax] = (int64_t)ceil((cmax - py->lo[ax]) / py->cell) - 1;
        if (rlo[ax] < 0) rlo[ax] = 0;
        if (rhi[ax] < rlo[ax]) rhi[ax] = rlo[ax];
        if (rhi[ax] > (int64_t)ceil(dom / py->cell) - 1)
          rhi[ax] = (int64_t)ceil(dom / py->cell) - 1;
        if (rlo[ax] > (int64_t)ceil(dom / py->cell) - 1)
          rlo[ax] = (int64_t)ceil(dom / py->cell) - 1;
        if (rhi[ax] < rlo[ax]) bad = 1;
      }
      if (bad) continue;
      int64_t sh = (int64_t)li + 1;
      for (int64_t az = rlo[2] >> sh; az <= rhi[2] >> sh; az++)
        for (int64_t ay = rlo[1] >> sh; ay <= rhi[1] >> sh; ay++)
          for (int64_t axx = rlo[0] >> sh; axx <= rhi[0] >> sh; axx++) {
            int64_t nid = axx + dims[li][0] * (ay + dims[li][1] * az);
            int32_t pos = hz_pyr_node_pos(py, li, nid);
            if (pos >= 0) lvls[li].off[pos + 1]++;
          }
    }
  }
  /* 2. префикс-суммы + рост пулов по ИТОГАМ (записей ≥ nt: bbox-множественность) */
  for (l = 0; l < py->nlev; l++)
    for (u = 0; u < lvls[l].n; u++)
      lvls[l].off[u + 1] += lvls[l].off[u];
  for (u = 0; u < py->nleaf; u++)
    leaf_off[u + 1] += leaf_off[u];
  {
    int64_t need = 0;
    int32_t base = 0;
    int32_t bases[256]; /* nlev ≤ 254 (А1568) + запас */
    for (l = 0; l < py->nlev; l++)
      need += lvls[l].off[lvls[l].n];
    if (need > *pool_cap) {
      int32_t *np = (int32_t *)realloc(*pool_io, (size_t)need * sizeof *np);
      if (!np) {
        free(dims);
        return -1;
      }
      *pool_io = np;
      *pool_cap = (int32_t)need;
    }
    if (leaf_off[py->nleaf] > *leaf_cap) {
      int32_t *np = (int32_t *)realloc(*leaf_io, (size_t)leaf_off[py->nleaf] * sizeof *np);
      if (!np) {
        free(dims);
        return -1;
      }
      *leaf_io = np;
      *leaf_cap = leaf_off[py->nleaf];
    }
    for (l = 0; l < py->nlev; l++) {
      bases[l] = base;
      base += lvls[l].off[lvls[l].n];
      lvls[l].pids = *pool_io + bases[l]; /* ДО заполнения (урок v1) */
    }
    /* 3. заполнение — тот же перебор, off/leaf_off как курсоры */
    for (p = 0; p < nt; p++) {
      uint8_t f = lpflo[p];
      int32_t tri = py->pcs[p].tri;
      int64_t lo[3], hi[3];
      if (f > (uint8_t)py->nlev) f = 0;
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
      if (reps_on && (int32_t)f >= g924_minf && (int32_t)f <= o->rep_nlev &&
          o->repof[((size_t)(int32_t)f - 1u) * (size_t)nt + (size_t)p] >= 0 && sw_far(o, p) &&
          (med_sub == NULL || med_sub[p]))
        continue; /* §924: кусок представлен — пишет проход D; §931: только
                   * дальний (ближний пишет себя, зеркально проходу 1);
                   * §932-Б: med_sub — замещённый средой кусок тоже пишет D */
      if (f == 0) {
        int fb = (o->walk != 0);
        for (int64_t iz = lo[2]; iz <= hi[2]; iz++)
          for (int64_t iy = lo[1]; iy <= hi[1]; iy++)
            for (int64_t ix = lo[0]; ix <= hi[0]; ix++) {
              int32_t dst[7];
              int32_t nn = sw923_leaf_or_nb(py, ix, iy, iz, fb, dst);
              for (int32_t k = 0; k < nn; k++)
                (*leaf_io)[leaf_off[dst[k]]++] = p;
            }
      } else {
        int32_t li = f - 1;
        int64_t sh = (int64_t)li + 1;
        for (int64_t az = lo[2] >> sh; az <= hi[2] >> sh; az++)
          for (int64_t ay = lo[1] >> sh; ay <= hi[1] >> sh; ay++)
            for (int64_t axx = lo[0] >> sh; axx <= hi[0] >> sh; axx++) {
              int64_t nid = axx + dims[li][0] * (ay + dims[li][1] * az);
              int32_t pos = hz_pyr_node_pos(py, li, nid);
              if (pos >= 0) {
                lvls[li].pids[lvls[li].off[pos]++] = p;
              } else { /* зеркало счёта: отказ узла — bbox в лист (fb А1578) */
                int fb = (o->walk != 0);
                for (int64_t iz = lo[2]; iz <= hi[2]; iz++)
                  for (int64_t iy = lo[1]; iy <= hi[1]; iy++)
                    for (int64_t ix = lo[0]; ix <= hi[0]; ix++) {
                      int32_t dst[7];
                      int32_t nn = sw923_leaf_or_nb(py, ix, iy, iz, fb, dst);
                      for (int32_t k = 0; k < nn; k++)
                        (*leaf_io)[leaf_off[dst[k]]++] = p;
                    }
              }
            }
      }
    }
    if (reps_on) { /* §924, проход D: представители в списки уровней */
      for (int32_t r = 0; r < o->nrep; r++) {
        if (want[r] <= 0) continue;
        int32_t li = want[r] - 1;
        int64_t rlo[3], rhi[3];
        int64_t triid = (int64_t)nt + r;
        int bad = 0;
        for (int ax = 0; ax < 3; ax++) {
          double cmin = o->tribox[6 * triid + ax];
          double cmax = o->tribox[6 * triid + 3 + ax];
          double dom = o->domhi[ax] - py->lo[ax];
          rlo[ax] = (int64_t)((cmin - py->lo[ax]) / py->cell);
          rhi[ax] = (int64_t)ceil((cmax - py->lo[ax]) / py->cell) - 1;
          if (rlo[ax] < 0) rlo[ax] = 0;
          if (rhi[ax] < rlo[ax]) rhi[ax] = rlo[ax];
          if (rhi[ax] > (int64_t)ceil(dom / py->cell) - 1)
            rhi[ax] = (int64_t)ceil(dom / py->cell) - 1;
          if (rlo[ax] > (int64_t)ceil(dom / py->cell) - 1)
            rlo[ax] = (int64_t)ceil(dom / py->cell) - 1;
          if (rhi[ax] < rlo[ax]) bad = 1;
        }
        if (bad) continue;
        int64_t sh = (int64_t)li + 1;
        for (int64_t az = rlo[2] >> sh; az <= rhi[2] >> sh; az++)
          for (int64_t ay = rlo[1] >> sh; ay <= rhi[1] >> sh; ay++)
            for (int64_t axx = rlo[0] >> sh; axx <= rhi[0] >> sh; axx++) {
              int64_t nid = axx + dims[li][0] * (ay + dims[li][1] * az);
              int32_t pos = hz_pyr_node_pos(py, li, nid);
              if (pos >= 0) lvls[li].pids[lvls[li].off[pos]++] = nt + r;
            }
      }
    }
    for (l = 0; l < py->nlev; l++) {
      for (u = lvls[l].n; u > 0; u--)
        lvls[l].off[u] = lvls[l].off[u - 1];
      lvls[l].off[0] = 0;
    }
    for (u = py->nleaf; u > 0; u--)
      leaf_off[u] = leaf_off[u - 1];
    leaf_off[0] = 0;
  }
  { /* §930-М2 (HZ_DBG930): занятость списков — медиана длины и доля
     * объёма мини-группы(4 подряд) от объёма списка: решает Б2 */
    static int m2 = -1;
    if (m2 < 0) m2 = getenv("HZ_DBG930") != NULL;
    if (m2) {
      for (l = 0; l < py->nlev; l++) {
        int64_t lens_sum = 0, lists = 0, grp_v = 0, lst_v = 0;
        for (u = 0; u < lvls[l].n; u++) {
          int32_t c = lvls[l].off[u + 1] - lvls[l].off[u];
          if (c <= 0) continue;
          lists++;
          lens_sum += c;
          if (c < 8) continue; /* группировать нечего */
          const int32_t *pid = lvls[l].pids + lvls[l].off[u];
          double lmin[3] = {0, 0, 0}, lmax[3] = {0, 0, 0};
          for (int32_t q = 0; q < c; q++) {
            int64_t tri = pid[q] >= nt ? pid[q] : py->pcs[pid[q]].tri;
            const double *tb = o->tribox + 6 * tri;
            for (int a = 0; a < 3; a++) {
              double mn = tb[a], mx = tb[3 + a];
              if (q == 0 || mn < lmin[a]) lmin[a] = mn;
              if (q == 0 || mx > lmax[a]) lmax[a] = mx;
            }
          }
          double lv = 1;
          for (int a = 0; a < 3; a++)
            lv *= lmax[a] - lmin[a];
          lst_v += (int64_t)lv;
          double gv_sum = 0;
          for (int32_t g0 = 0; g0 < c; g0 += 4) {
            int32_t g1 = g0 + 4 < c ? g0 + 4 : c;
            double gmin[3], gmax[3];
            for (int32_t q = g0; q < g1; q++) {
              int64_t tri = pid[q] >= nt ? pid[q] : py->pcs[pid[q]].tri;
              const double *tb = o->tribox + 6 * tri;
              for (int a = 0; a < 3; a++) {
                double mn = tb[a], mx = tb[3 + a];
                if (q == g0 || mn < gmin[a]) gmin[a] = mn;
                if (q == g0 || mx > gmax[a]) gmax[a] = mx;
              }
            }
            double gv = 1;
            for (int a = 0; a < 3; a++)
              gv *= gmax[a] - gmin[a];
            gv_sum += gv;
          }
          grp_v += (int64_t)gv_sum;
        }
        fprintf(stderr, "M2 L%d: списков=%ld ср.длина=%.1f грпп/спск=%.3f\n", l, (long)lists,
                lists ? (double)lens_sum / (double)lists : 0.0,
                lst_v > 0 ? (double)grp_v / (double)lst_v : -1.0);
      }
    }
  }
  if (dbg) {
    int64_t fmax = 0, own, fin, nl = 0, entries = 0;
    fprintf(stderr, "DBG923:");
    for (p = 0; p < nt; p++)
      if (lpflo[p] > fmax) fmax = lpflo[p];
    for (l = 0; l < py->nlev; l++) {
      own = fin = 0;
      for (u = 0; u < lvls[l].n; u++) {
        int32_t c = lvls[l].off[u + 1] - lvls[l].off[u];
        if (c > 0) own++;
        entries += c;
        if (lvls[l].finer[u]) fin++;
      }
      fprintf(stderr, " L%d:%d/%d", l, (int)own, (int)fin);
    }
    for (u = 0; u < py->nleaf; u++) {
      if (leaf_off[u + 1] - leaf_off[u] > 0) nl++;
      entries += leaf_off[u + 1] - leaf_off[u];
    }
    fprintf(stderr, " leaves=%d fmax=%d nlev=%d entries=%d\n", (int)nl, (int)fmax, py->nlev,
            (int)entries);
  }
  if (dbg && bstart && bpids) { /* сверка листовых списков с bpids (только все-0) */
    int64_t dmis = 0, dextra = 0, leafmis = 0, fmax2 = 0;
    for (p = 0; p < nt; p++)
      if (lpflo[p] > fmax2) fmax2 = lpflo[p];
    for (int32_t li2 = 0; li2 < py->nleaf && fmax2 == 0; li2++) {
      for (int32_t k = leaf_off[li2]; k < leaf_off[li2 + 1]; k++) {
        int32_t pv = (*leaf_io)[k];
        int found = 0;
        for (int32_t q = bstart[li2]; q < bstart[li2 + 1]; q++)
          if (bpids[q] == pv) found = 1;
        if (!found) dextra++;
      }
      for (int32_t q = bstart[li2]; q < bstart[li2 + 1]; q++) {
        int found = 0;
        for (int32_t k = leaf_off[li2]; k < leaf_off[li2 + 1]; k++)
          if ((*leaf_io)[k] == bpids[q]) found = 1;
        if (!found) {
          dmis++;
          if (leafmis < 8) {
            int32_t tri2 = py->pcs[bpids[q]].tri;
            fprintf(stderr, "DBG923mis leaf=%d p=%d tri=%d\n", li2, bpids[q], tri2);
            leafmis++;
          }
        }
      }
    }
    fprintf(stderr, "DBG923verify mis=%lld extra=%lld", (long long)dmis, (long long)dextra);
    {
      int64_t nv = 0;
      for (int32_t li2 = 0; li2 < py->nleaf; li2++)
        /* границу g_dbg923_vis[1<<20] проверяем ПЕРВОЙ: cppcheck
         * arrayIndexThenCheck читал порядок условий как «индекс раньше
         * проверки» (07-10). Семантика та же — && без побочных эффектов */
        if (li2 < (1 << 20) && leaf_off[li2 + 1] - leaf_off[li2] > 0 && !g_dbg923_vis[li2]) nv++;
      fprintf(stderr, " unvisited=%lld (nleaf=%d)", (long long)nv, py->nleaf);
      memset(g_dbg923_vis, 0, sizeof g_dbg923_vis);
    }
    fprintf(stderr, "\n");
  }
  free(dims);
  return 0;
}

/* марш трубки через узел (уровень l, позиция pos; l<0 — лист) на отрезке
 * [tin,tout]; carry a/b сквозной — сечение трубки не делится (А1563).
 * §923: спуск промежуточного узла вынесен в front_descend — новый путь
 * (раскладка по уровням) вызывает его же для спуска к мелким кускам. */
static void front_descend(front_ctx *fc, int32_t l, int32_t pos, double tin, double tout, double *a,
                          double *b);

static void front_visit(front_ctx *fc, int32_t l, int32_t pos, double tin, double tout, double *a,
                        double *b) {
  const hz_pyr *py = fc->py;
  double blo[3], bhi[3];
  if (fc->lvls) { /* §923: раскладка по уровням владения */
    if (l < 0) {  /* лист: только floor-0 куски этой клетки */
      int32_t n = fc->leaf_loff[pos + 1] - fc->leaf_loff[pos];
      if (fc->cfront) fc->ncellfront++;
      if (n > 0) {
        fc->dbg_leaf++;
        fc->dbg_n += n;
      }
      if (pos < (1 << 20)) g_dbg923_vis[pos] = 1;
      if (n > 0) {
        const int32_t *ps = fc->leaf_lpids + fc->leaf_loff[pos];
        if (fc->o->walk)                                   /* А1576: точный многопопадный проход */
          front_seg_walk(fc, ps, n, tin, tout, a, b, pos); /* Ш9: лист — в зонд */
        else {
          front_node_box(fc, -1, pos, blo, bhi);
          front_interact(fc, ps, n, blo, bhi, tin, tout, a, b, 1);
        }
      }
      return;
    }
    {
      const sw_lvl *L = &fc->lvls[l];
      int32_t n = L->off[pos + 1] - L->off[pos];
      if (n > 0) {
        fc->dbg_node++;
        fc->dbg_n += n;
      }
      if (n > 0) { /* СВОИ куски (этаж l+1): взаимодействие на сегменте узла.
                    * Пустота узла не обходится — она покрыта одним
                    * слэб-тестом сегмента (ход по пустоте уровнем, А1569). */
        const int32_t *ps = L->pids + L->off[pos];
        front_node_box(fc, l, pos, blo, bhi);
        if (fc->o->walk)
          front_seg_walk(fc, ps, n, tin, tout, a, b, -1); /* Ш9: узел — вне клеточного зонда */
        else
          front_interact(fc, ps, n, blo, bhi, tin, tout, a, b, 0);
      }
      if (L->finer[pos]) /* спуск — ТОЛЬКО если под узлом есть более
                          * мелкие куски; иначе марш стоит на уровне */
        front_descend(fc, l, pos, tin, tout, a, b);
      return;
    }
  }
  if (l < 0) {        /* лист: базовая клетка */
    fc->agg_list = 0; /* листовой список не отсортирован по материалу */
    if (fc->cfront) fc->ncellfront++;
    front_node_box(fc, -1, pos, blo, bhi);
    if (fc->o->walk) /* А1576: точный многопопадный проход */
      front_seg_walk(fc, fc->bpids + fc->bstart[pos], fc->bstart[pos + 1] - fc->bstart[pos], tin,
                     tout, a, b, pos); /* Ш9: лист — в зонд */
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
      /* §863/шаг 2 (этап 0): сортировка по (материал, id) + дедупликация —
       * кусок, чей bbox накрывает несколько листов узла, попадал в список
       * многократно. §872: список остаётся КУСКАМИ (не гидами) при любом
       * agg — серия по материалу это единица culling, не депозит; сортировка
       * делается ДО решения «кэшировать ли», чтобы pbuf был годен и без кэша */
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
      {
        int32_t un = 0;
        for (int32_t q = 0; q < n; q++)
          if (un == 0 || fc->pbuf[un - 1] != fc->pbuf[q]) fc->pbuf[un++] = fc->pbuf[q];
        n = un;
      }
      if (cid >= 0 && *fc->nstore_n + n <= fc->nstore_cap) {
        int64_t st = *fc->nstore_n;
        for (int32_t q = 0; q < n; q++)
          fc->nstore[st + q] = fc->pbuf[q];
        *fc->nstore_n = st + n;
        fc->nstart[cid] = st;
        fc->nlen[cid] = n;
        ps = fc->nstore + st;
      } /* магазин полон: работаем по pbuf без кэша (он тоже отсортирован) */
      fc->agg_list = fc->agg; /* список отсортирован (mtl,id): серии подряд */
    }
    front_node_box(fc, l, pos, blo, bhi);
    /* А1573: на узловом отрезке клеточный T-перехват НЕзаконен — bbox куска
     * ⊆ узел не означает, что трубка задевает его тень; перехват раздаёт
     * сечение узла кускам, мимо которых трубка прошла (расходимость
     * ~×1.5/итерацию на cavity05 lp=2). Узел — только точное пересечение.
     * Предикат А1566 «bbox ⊆ клетка» валиден ТОЛЬКО на листовом отрезке
     * (bbox ⊆ БАЗОВАЯ клетка). */
    if (fc->o->walk)                                  /* А1576: точный многопопадный проход */
      front_seg_walk(fc, ps, n, tin, tout, a, b, -1); /* Ш9: узел — вне зонда */
    else
      front_interact(fc, fc->pbuf, n, blo, bhi, tin, tout, a, b, 0);
    return;
  }
  /* ПРОМЕЖУТОЧНЫЙ: спуск — продолжение DDA на детском уровне (§923:
   * вынесен в front_descend, общий с новым путём) */
  front_descend(fc, l, pos, tin, tout, a, b);
}

static void front_descend(front_ctx *fc, int32_t l, int32_t pos, double tin, double tout, double *a,
                          double *b) {
  const hz_pyr *py = fc->py;
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

/* §865/А1580 (path=1): попадания списка кусков в список полного пути */
static void path_collect_list(front_ctx *fc, const int32_t *ps, int32_t n, double tin,
                              double tout) {
  const hz_pyr *py = fc->py;
  int32_t u;
  /* §865/раунд 13: кэш «кусок × трубка» выключен (нет памяти при постройке —
   * fail closed) → прежний путь: bbox×сегмент + ray-tri на каждый визит */
  if (!fc->pc) {
    for (u = 0; u < n; u++) {
      int32_t p = ps[u];
      double tt, bh0, bh1;
      if (fc->pstamp[p] == fc->pkey) continue; /* штамп ДО ray-tri (А1564) */
      /* §865: префильтр bbox×сегмент — консервативный, результата не меняет */
      if (!front_box_seg(fc->o->tribox + 6 * (int64_t)py->pcs[p].tri,
                         fc->o->tribox + 6 * (int64_t)py->pcs[p].tri + 3, fc->org, fc->om, tin,
                         tout, &bh0, &bh1))
        continue;
      tt = sw_ray_tri_raw(fc->org, fc->om, fc->o->trivert + 9 * (int64_t)py->pcs[p].tri);
      if (!(tt >= tin) || !(tt <= tout)) continue;
      if (fc->pbuf_n == fc->pbuf_cap) {
        int64_t nc = fc->pbuf_cap ? fc->pbuf_cap * 2 : 256;
        /* Последовательно и с записью владения ДО следующего realloc: при
         * отказе второго старый блок не теряется (clang-analyzer-unix.Malloc,
         * «potential leak of np/nt», 07-10). pbuf_cap растёт только когда
         * удались оба — меньший cap безопасен (индексы < pbuf_n ≤ cap). */
        double *nt = (double *)realloc(fc->pbuf_t, (size_t)nc * sizeof *nt);
        if (!nt) return;
        fc->pbuf_t = nt;
        int32_t *np = (int32_t *)realloc(fc->pbuf_p, (size_t)nc * sizeof *np);
        if (!np) return;
        fc->pbuf_p = np;
        fc->pbuf_cap = nc;
      }
      fc->pstamp[p] = fc->pkey;
      fc->pbuf_t[fc->pbuf_n] = tt;
      fc->pbuf_p[fc->pbuf_n] = p;
      fc->pbuf_n++;
    }
    return;
  }
  for (u = 0; u < n; u++) {
    int32_t p = ps[u];
    pc_rec *r = fc->pc + p;
    double tt, bh0, bh1;
    if (r->pstamp == fc->pkey) continue; /* уже в списке попаданий трубки */
    if (r->tstamp != fc->pkey) {         /* §865/р.13: bbox-интервал — константа трубки */
      if (!front_box_seg(fc->o->tribox + 6 * (int64_t)py->pcs[p].tri,
                         fc->o->tribox + 6 * (int64_t)py->pcs[p].tri + 3, fc->org, fc->om, fc->pth0,
                         fc->pth1, &bh0, &bh1)) {
        r->bh0 = 1.0; /* пустой интервал: [1,0] */
        r->bh1 = 0.0;
        r->tstamp = fc->pkey;
        continue;
      }
      r->bh0 = bh0;
      r->bh1 = bh1;
      r->tstamp = fc->pkey;
    } else {
      bh0 = r->bh0;
      bh1 = r->bh1;
    }
    /* клип интервала к [tin,tout] пуст ⟺ box_seg по сегменту дал бы 0 —
     * сегмент ⊂ трубки, ровно та же проверка без повторного box_seg
     * (маркер пустого интервала [1,0] отсекается первым тестом) */
    if (!(bh0 <= bh1) || bh1 < tin || bh0 > tout) continue;
    if (r->ttstamp != fc->pkey) { /* ray-tri — ОДИН РАЗ на (кусок, трубку) */
      r->tt = sw_ray_tri_raw(fc->org, fc->om, fc->o->trivert + 9 * (int64_t)py->pcs[p].tri);
      r->ttstamp = fc->pkey;
    }
    tt = r->tt;
    if (!(tt >= tin) || !(tt <= tout)) continue;
    if (fc->pbuf_n == fc->pbuf_cap) {
      int64_t nc = fc->pbuf_cap ? fc->pbuf_cap * 2 : 256;
      double *nt = (double *)realloc(fc->pbuf_t, (size_t)nc * sizeof *nt);
      if (!nt) return;
      fc->pbuf_t = nt;
      int32_t *np = (int32_t *)realloc(fc->pbuf_p, (size_t)nc * sizeof *np);
      if (!np) return;
      fc->pbuf_p = np;
      fc->pbuf_cap = nc;
    }
    r->pstamp = fc->pkey;
    fc->pbuf_t[fc->pbuf_n] = tt;
    fc->pbuf_p[fc->pbuf_n] = p;
    fc->pbuf_n++;
  }
}

/* §865/раунд 14: статическое перекрытие tribox [minx,miny,minz,maxx,maxy,
 * maxz] и бокса [blo,bhi] — без луча; касание засчитывается (разделение
 * только при строгом >), так что фильтр списка узла попаданий не теряет */
static int tri_box_overlap(const double *tb, const double blo[3], const double bhi[3]) {
  int ax;
  for (ax = 0; ax < 3; ax++)
    if (tb[ax] > bhi[ax] || tb[ax + 3] < blo[ax]) return 0;
  return 1;
}

/* §865/А1580 (path=1): выбор списка кусков клетки/узла и сбор попаданий */
static void path_collect(front_ctx *fc, int32_t l, int32_t pos, double tin, double tout) {
  const hz_pyr *py = fc->py;
  const int32_t *ps;
  int32_t n = 0;
  double blo[3], bhi[3];
  if (l < 0) {
    ps = fc->bpids + fc->bstart[pos];
    n = fc->bstart[pos + 1] - fc->bstart[pos];
  } else {
    int64_t cid;
    front_node_box(fc, l, pos, blo, bhi);
    cid = (int64_t)fc->ncluster0 + fc->noff[l] + pos;
    if (cid >= 0 && fc->nstart[cid] >= 0) { /* §863: кэш списков узлов */
      ps = fc->nstore + fc->nstart[cid];
      n = fc->nlen[cid];
    } else {
      if (front_gather(fc, l, pos, &n) != 0) return;
      /* §865/раунд 14: обрезка списка по боксу узла ДО записи в кэш — кусок
       * вне бокса узла попасть в марше этого узла не может (точка попадания
       * лежала бы в обоих боксax); фильтр статичен по геометрии — один раз */
      {
        int32_t m = 0, q;
        for (q = 0; q < n; q++) {
          int32_t p = fc->pbuf[q];
          if (tri_box_overlap(fc->o->tribox + 6 * (int64_t)py->pcs[p].tri, blo, bhi))
            fc->pbuf[m++] = p;
        }
        n = m;
      }
      ps = fc->pbuf;
      if (cid >= 0 && *fc->nstore_n + n <= fc->nstore_cap) {
        /* кэшируем собранный список (как ветка material во front_visit) */
        int64_t st = *fc->nstore_n, q;
        for (q = 0; q < n; q++)
          fc->nstore[st + q] = fc->pbuf[q];
        *fc->nstore_n = st + n;
        fc->nstart[cid] = st;
        fc->nlen[cid] = n;
        ps = fc->nstore + st;
      }
    }
  }
  path_collect_list(fc, ps, n, tin, tout);
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
  fc->hop_n = 0; /* §873/T4: состояние хопов — на трубку */
  fc->hop_depth = 0;
  fc->hop_lost = 0.0;
  fc->med_nact = 0;  /* §932-Б: активные среды — на трубку (штампы гасят повторы) */
  fc->t_last = -1.0; /* §932-А: старт трубки — события не было; здесь, а не
                      * memset'ом fc (он даёт 0.0). Прямой свет (пустотный
                      * вход трубки) детальность НЕ трогает (А1691) */
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
    fc->pth0 = h0; /* §865/раунд 13: границы трубки для кэша кусков */
    fc->pth1 = h1;
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
      fc->Ed[p] += fc->w_d * Lin * csec * fc->axcos / front_depden(fc, p);
      fc->absorbed += fc->w_d * Lin * csec;
      fc->emitted += fc->w_d * front_le(fc, p) * csec;
      {
        double lhcap_p = front_le(fc, p);
        Lh = front_lh_cap(fc,
                          front_le(fc, p) + front_rho(fc, p) * fc->Eprev[p] / (2.0 * M_PI) > lhcap_p
                              ? lhcap_p
                              : front_le(fc, p) + front_rho(fc, p) * fc->Eprev[p] / (2.0 * M_PI));
      }
      fc->recycled += fc->w_d * (Lh - fc->le) * csec;
      fc->ndep++;
      sw_accum(fc, p, fc->w_d * Lin * csec * fc->axcos / front_depden(fc, p));
      Lin = Lh;
    }
    fc->row_dep += 1.0; /* §887-b: маркер исполнения ветки (счётчик визитов) */
    a = 0.0;
    b = Lin;
  } else {
    front_visit(fc, py->nlev - 1, 0, h0, h1, &a, &b);
    /* G6 на корне: вошедшая (нулевая — тёмный вход) доля = Σdep + вышедшая.
     * §873/T4: при зеркальных хопах вход лег NOT нулевой (hop_lin), проверка
     * корня не применима — честный баланс хопов считает lost/absorbed. */
    if (!fc->noprop && fc->hop_n == 0 && fc->hop_depth == 0 &&
        fabs(0.0 - fc->depA - a) > 1e-9 * (1.0 + fabs(a)))
      fc->g6viol++;
    /* §873/T4 шаг 2: зеркальные хопы — LIFO-обработка очереди.
     * §894-c-5: ноги живут в СВОЁМ штамп-пространстве, поэтому перед циклом
     * сохраняем «основные» pstamp/pkey. Восстановление обязано стоять после
     * цикла (ниже): без него штампы основной трубки навсегда оставались бы в
     * пространстве ног — ровно то, что §894-c-5 и разделял. Сейчас не
     * срабатывает (ноги выключены, HZ_SW_LEG_BUDGET=0), но сохранение без
     * восстановления — дефект, а не экономия. */
    int64_t *main_pstamp = fc->pstamp;
    int64_t main_pkey = fc->pkey;                  /* тип как у fc->pkey: иначе -Wsign-conversion */
    fc->in_leg = 1;                                /* §894-c: дальше обрабатываются НОГИ */
    fc->leg_pkey = ((uint64_t)fc->mark << 48) | 1; /* уникальный ключ на ногу */
    while (fc->hop_n > 0) {
      /* §894-c-5: нога живёт в СВОЁМ штамп-пространстве — штампы ног не
       * блокируют депозиты основной трубки (гипотеза аудита §896-3) */
      if (fc->leg_budget <= 0) { /* §896-5: бюджет исчерпан — остаток в lost */
        fc->hop_lost += fc->hop_lin[fc->hop_n - 1];
        fc->hop_n--;
        continue;
      }
      fc->leg_budget--;
      fc->pstamp = fc->leg_pstamp;
      fc->pkey = (int64_t)fc->leg_pkey; /* ключ — битовый узор, знак не несёт:
                                         * явное приведение вместо -Wsign-conversion */
      double ai = fc->hop_lin[fc->hop_n - 1];
      const double *pt = fc->hop_pt[fc->hop_n - 1];
      const double *dir = fc->hop_dir[fc->hop_n - 1];
      fc->hop_n--;
      fc->hop_depth++;
      fc->leg_pkey++;                         /* уникальный ключ каждой ноге */
      fc->leg_root = fc->hop_root[fc->hop_n]; /* корень цепи — порог от него */
      if (fc->hop_dep[fc->hop_n] >= 32) {     /* §896-4: цепь длиннее 32 — в lost */
        fc->hop_lost += ai;
        fc->hop_n--;
        continue;
      }
      fc->cur_leg_depth = fc->hop_dep[fc->hop_n];
      /* §881/А1618: усечение по глубине снято — нога ниже порога от корня не
       * рождается при спавне; цикл зеркального коридора обрывается порогом */
      for (q = 0; q < 3; q++) {
        fc->org[q] = pt[q];
        fc->omcur[q] = dir[q]; /* om указывает на omcur — направление ноги */
      }
      fc->t_last = 0.0; /* §932-А1704: нога стартует В материальном событии —
                         * первый хит ноги = полный пробег от зеркала */
      if (!front_box_seg(py->lo, fc->hi, fc->org, fc->om, 0.0, 1e30, &h0, &h1)) continue;
      a = ai;
      b = 0.0;
      fc->depA = 0.0; /* G6 хоп-ноги: вошедшая доля = ai */
      front_visit(fc, py->nlev - 1, 0, h0, h1, &a, &b);
      if (!fc->noprop && fc->hop_n == 0 && fabs(ai - fc->depA - a - b) > 1e-9 * (1.0 + fabs(ai)))
        fc->g6viol++;
      if (b > 0.0) {
        *lostB += b; /* непогашенный хвост хоп-ноги вышел за домен */
        fc->nlostseg++;
      }
    }
    if (fc->hop_lost > 0.0) {
      *lostB += fc->hop_lost; /* вытесненные из очереди доли — в lost */
      fc->nlostseg++;
    }
    fc->in_leg = 0;
    /* §894-c-5: возвращаем основное штамп-пространство (см. сохранение выше) */
    fc->pstamp = main_pstamp;
    fc->pkey = main_pkey;
  }
  if (a > 0.0 || b > 0.0)
    fc->nlostseg++; /* остаток ушёл за границу домена — в lost, не исчез (А1567) */
  *lostA += a;
  *lostB += b;
  fc->hop_lost_sum += fc->hop_lost;
  fc->hop_hops_sum += fc->hop_depth;
}

/* направление целиком: трубки = линии базовой сетки (дедупликация штампом
 * клетки, как §851) */
static void front_dir(front_ctx *fc, int32_t *stampv, const double *odir) {
  const hz_pyr *py = fc->py;
  int32_t ax;
  int32_t f[3];
  for (ax = 0; ax < 3; ax++) {
    int64_t n = ax == 0 ? py->nx : (ax == 1 ? py->ny : py->nz);
    /* §873/T4: omcur — направление ТЕКУЩЕЙ ноги; на старте направления
     * восстанавливается из таблицы (хопы прошлой трубки его портили) */
    fc->omcur[ax] = odir[ax];
    if (fabs(odir[ax]) < 1e-30) {
      f[ax] = -1;
      continue;
    }
    f[ax] = odir[ax] > 0.0 ? 0 : (int32_t)(n - 1);
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
        /* §923-дых: печать в конце front_dir (после семейств осей) */
      }
  }
  if (getenv("HZ_DBG923"))
    fprintf(stderr,
            "DBG923dir tubes=%lld node=%lld leaf=%lld pieces=%lld absorbed=%.5g emitted=%.5g\n",
            (long long)fc->ntube, (long long)fc->dbg_node, (long long)fc->dbg_leaf,
            (long long)fc->dbg_n, fc->absorbed, fc->emitted);
  if (fc->med_on) { /* §932-Б: ФЛЕШ отложенных депозитов — один на направление:
                     * Σ_kid Ed·A = pend (Ш3: pend в единицах потока) */
    for (int mk = 0; mk < fc->med_nact; mk++) /* хвост ПОСЛЕДНЕЙ трубки направления */
      fc->med_pend[fc->med_act[mk].r] += fc->med_act[mk].pend;
    fc->med_nact = 0;
    for (int32_t r = 0; r < fc->nrep; r++) {
      if (fc->med_pend[r] != 0.0) {
        /* Ш3/Г-flux: pend уже в единицах потока (axcos события своей
         * семьи внутри банка) — флеш БЕЗ домножения; Σ_kid edp·A = pend */
        double edp = fc->med_pend[r] / fc->med_asum[r];
        g932b_edf += fc->med_pend[r]; /* Ш3: факт-Ed флеша — сверка с деп */
        for (int32_t q = fc->med_koff[r]; q < fc->med_koff[r + 1]; q++)
          fc->Ed[fc->med_kmem[q]] += edp;
        fc->med_pend[r] = 0.0;
      }
    }
  }
}

/* А1578: индекс клетки по оси — ТА ЖЕ конвенция, что pyr_axis_index в pyr.c
 * (floor с клэмпом к сетке): клетки, помеченные pyr, гарантированно листья, и
 * марш их достигает (фантомную колонку у max-грани домена — вырожденным
 * сегментом [h,h], §852). Прежний span «floor…ceil−1» с клэмпом к домену клал
 * кусок плоской max-грани в слой, где листа чаще нет, — кусок выпадал из
 * bbox-индекса, и стенка теряла депозиты (А1578: 0.22 события/трубку). */

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
  g931_rh = 0; /* §931/А1668: счётчик шва — на прогон (вызовы многократны) */
  g931_seam = 0;
  g932_n = 0; /* §932-А: приборы пробега — на прогон (А1724: печать при
               * любом travel=1, включая tvc=0-НК) */
  g932_fsum = 0.0;
  g932_lvsum = 0.0;
  memset(g932_hist, 0, sizeof g932_hist);
  g932b_ev = 0; /* §932-Б: приборы среды — на прогон */
  g932b_big = 0;
  g932b_abs = 0.0;
  g932b_dep = 0.0;
  if (g932_nomax < 0) g932_nomax = getenv("HZ_TVNOMAX") != NULL; /* §932-А1715 НК */
  double *cov = NULL;     /* §911-6: Σ долей следа на кусок (HZ_COVDBG) */
  uint32_t *covn = NULL;  /* §911-6: число хитов на кусок */
  double *cov_lit = NULL; /* §911-8: освещённое покрытие */
  double *cov_e = NULL;   /* §911-10: вклад с амплитудой света */
  int32_t *order = NULL, *vindex = NULL;
  uint64_t *bits = NULL;
  hz_sw_walk *walks = NULL;
  hz_sw_dir *tab = NULL;
  int32_t *stampv = NULL;     /* §851: штампы визитов линий [ncells] */
  int64_t *pstamp = NULL;     /* §852/А1564: штамп (трубка × кусок) [nt] */
  int64_t *leg_pstamp = NULL; /* §894-c-5: штампы ног (отдельное пространство) */
  /* §865/раунд 13: кэш «кусок × трубка» (только path=1) */
  pc_rec *pc = NULL; /* §865/раунды 13+15: кэш «кусок × трубка» (только path=1) */
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
  uint8_t *lpflo = NULL;     /* §862: этажи (когда задан lpacc) — ПРОСТРАНСТВО СЛОТОВ */
  uint8_t *lpflo_tri = NULL; /* §922/А1632: то же в пространстве ИСХОДНЫХ
                              * треугольников — hz_pyr_set_lp читает
                              * lp[pcs[u].tri] (А1491: слот ≠ треугольник);
                              * прежняя передача lpflo напрямую давала этажи
                              * не на своих кусках */
  /* §923: раскладка кусков по уровням владения (только lpacc&&lpapply;
   * agg-режим — прежний путь: списки уровней не сортированы по mtl) */
  sw_lvl *lvls = NULL;
  int32_t *lvls_loff = NULL;                    /* [nleaf+1] */
  int32_t *lvl_off_all = NULL;                  /* плоские массивы уровней: анализатор не
                                                 * держит nlev указателей (FP CWE-401) */
  uint8_t *lvl_finer_all = NULL;                /* §923 */
  int32_t *lvl_pool = NULL, *lvls_lpids = NULL; /* растятся в sw_levels_build */
  int32_t lvl_pool_cap = 0, lvls_lpids_cap = 0;
  /* §924: представители — метки want (дедуп гроздей) и агрегат Eprev */
  int32_t *rep_want = NULL;
  double *rep_ep = NULL;
  /* §928: сбор калибровки k(r,ω) — числитель/знаменатель (HZ_KCALDUMP) */
  double *knum = NULL, *kden = NULL;
  const char *kcalpath = NULL;
  int lvls_ready = 0; /* §923: списки собраны хоть раз; до того fc.lvls = NULL
                       * (ловля §923-дым: it=0 шёл по ПУСТЫМ спискам — вся
                       * эмиссия первой итерации терялась, E 0.0000 на it=1) */

  /* ВНИМАНИЕ: эти объявления перенесены сюда 07-10 и должны оставаться ВЫШЕ
   * первого `goto done` (сейчас это `rc = hz_sw_dir_table(...)` ниже).
   * Раньше они стояли в теле функции, и прыжок на метку `done:` перескакивал
   * их инициализаторы: инициализатор блочной переменной выполняется в точке
   * объявления, поэтому на пути ошибки `done:` делал free() по НЕОПРЕДЕЛЁННЫМ
   * указателям. Ловилось gcc -Wjump-misses-init и -Wmaybe-uninitialized
   * (15 указателей), то есть это была настоящая ошибка, а не шум. */
  double *Linmax_prev = NULL; /* §893-b: пер-кусковый максимум пришедшего
                               * радианса (принцип максимума); NULL вне mode 3,
                               * free(NULL) легален */
  double *Linmax_cur = NULL;
  /* §932-Б: массивы среды — владелец hz_sw_run, агрегаты пересчитываются
   * со sw_levels_build (А1711); med_sub — решение замещения куска */
  uint8_t *med_sub = NULL;
  double *med_box = NULL, *med_an = NULL, *med_v = NULL, *med_rho = NULL, *med_asum = NULL;
  int32_t *med_nsub = NULL, *med_koff = NULL, *med_kmem = NULL;
  int64_t *med_stamp = NULL;
  double *med_pend = NULL;
  int med_ready = 0;
  double *ced = NULL, *casum = NULL; /* §932-Б/Ш9-ЗОНД */
  double *cabs = NULL;               /* §932-Б/Ш12: радианс входа клетки (cellc=2) */

  memset(st, 0, sizeof *st);
  if (!py || !area || !nrm || !kd || !o) return 1;
  if (o->iters <= 0) return 1;
  if (!py->pcs || nt != py->nt) {
    fprintf(stderr, "DBG nt=%d py->nt=%d\n", nt, py->nt);
    return 1;
  }
  if (o->mode == 2 && o->vc && !o->trivert) return 1; /* §852: нужно точное пересечение */
  if (o->mode == 3 && (!o->trivert || !o->tribox || !o->lp || !o->domhi))
    return 1;                                 /* §852: фронт */
  if (o->mode == 3 && !py->leaf_lp) return 1; /* §852/А1567: per-node max ℓ_p не задан */
  /* §843: режим 0 (скалярный марш §835) определён только на легаси-наборах */
  if (o->mode == 0 && o->ndirs != 6 && o->ndirs != 26) return 1;
  if (o->travel && !o->lpacc) { /* §932-А/МЕЛОЧИ-3: fail-closed и для иных
                                 * вызывающих, не только pgather (А1705) */
    fprintf(stderr, "hz_sw_run: travel=1 требует lpacc (adapt=1)\n");
    return 2;
  }
  if (o->travel && o->tvc < 0.0) { /* §932-А1705: отрицательный масштаб — отказ */
    fprintf(stderr, "hz_sw_run: travel=1 требует tvc>=0\n");
    return 2;
  }
  if (o->medium && (o->med_tau <= 0.0 || o->lpacc_e == NULL || o->nrep <= 0 || o->repof == NULL ||
                    o->reps_collect || !o->walk)) {
    /* §932-Б/А1738: fail-closed — среда требует reps-инфраструктуру,
     * адаптив (lpacc_e) и выключена в мире kcal-сбора */
    fprintf(stderr, "hz_sw_run: medium=1 требует mtau>0, lpacc_e, nrep>0, walk, !reps_collect\n");
    return 2;
  }
  rc = hz_sw_dir_table(o->ndirs, &tab, &nd);
  if (rc != 0) goto done;
  csec = py->cell * py->cell;

  /* §928: сбор калибровки — только мир сбора (rep-карты есть, представители
   * НЕ активны: o->reps_collect). Накопление — последняя итерация. */
  kcalpath = getenv("HZ_KCALDUMP");
  if (sw_med_now < 0) sw_med_now = getenv("HZ_MEDNOW") != NULL;      /* Ш3: Вопрос-1 */
  if (sw_med_t1 < 0) sw_med_t1 = getenv("HZ_MEDT1") != NULL;         /* Ш3: зонд exp */
  if (sw_med_scr < 0) sw_med_scr = getenv("HZ_MEDSCRAMBLE") != NULL; /* Ш4: НК А1698 */
  if (kcalpath != NULL && kcalpath[0] != '\0' && o->reps_collect && o->nrep > 0 &&
      o->repof != NULL) {
    knum = (double *)calloc((size_t)o->nrep * HZ_KCAL_NW, sizeof *knum);
    kden = (double *)calloc((size_t)o->nrep * HZ_KCAL_NW, sizeof *kden);
    if (!knum || !kden) { /* fail closed: мир без калибровки */
      free(knum);
      free(kden);
      knum = NULL;
      kden = NULL;
    }
  }

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
  if (o->mode == 3) { /* §893-b: пер-кусковый max пришедшего радианса */
    Linmax_prev = (double *)calloc((size_t)nt, sizeof *Linmax_prev);
    Linmax_cur = (double *)calloc((size_t)nt, sizeof *Linmax_cur);
  }
  if (o->mode == 3 && o->lpacc && o->lpapply && py->nlev > 0 && !o->agg) {
    /* §924: буферы представителей — только живая лестница+walk;
     * §928: reps_collect — представители НЕ активны (сбор в кусочном мире) */
    if (o->nrep > 0 && o->walk && o->repof != NULL && !o->reps_collect) {
      rep_want = (int32_t *)calloc((size_t)o->nrep, sizeof *rep_want);
      rep_ep = (double *)calloc((size_t)o->nrep, sizeof *rep_ep);
      if (!rep_want || !rep_ep) {
        free(rep_want);
        free(rep_ep);
        rep_want = NULL;
        rep_ep = NULL; /* fail closed: мир без представителей */
      }
    }
    /* §923: раскладка по уровням владения. Динамическая лестница —
     * единственный гейт (А1638: статика lp= и без-ключевой мир — прежние
     * пути; agg — прежний путь, списки уровней не сортированы по mtl). */
    int32_t nlv = py->nlev;
    int64_t offcnt = 0, fincnt = 0;
    for (int32_t l2 = 0; l2 < nlv; l2++) {
      offcnt += (int64_t)py->nlev_nodes[l2] + 1;
      fincnt += py->nlev_nodes[l2];
    }
    lvls = (sw_lvl *)calloc((size_t)nlv, sizeof *lvls);
    lvl_off_all = (int32_t *)calloc((size_t)(offcnt > 0 ? offcnt : 1), sizeof *lvl_off_all);
    lvl_finer_all = (uint8_t *)calloc((size_t)(fincnt > 0 ? fincnt : 1), 1);
    lvls_loff = (int32_t *)calloc((size_t)py->nleaf + 1, sizeof *lvls_loff);
    if (lvls && lvl_off_all && lvl_finer_all && lvls_loff) {
      int64_t o2 = 0, f2 = 0;
      for (int32_t l2 = 0; l2 < nlv; l2++) {
        lvls[l2].n = py->nlev_nodes[l2];
        lvls[l2].off = lvl_off_all + o2;
        lvls[l2].finer = lvl_finer_all + f2;
        o2 += (int64_t)py->nlev_nodes[l2] + 1;
        f2 += py->nlev_nodes[l2];
      }
    } else { /* fail closed: мир без раскладки (прежние пути, §863-принцип) */
      free(lvls);
      lvls = NULL;
      free(lvl_off_all);
      lvl_off_all = NULL;
      free(lvl_finer_all);
      lvl_finer_all = NULL;
      free(lvls_loff);
      lvls_loff = NULL;
    }
  }
  /* §932-Б: массивы среды объявлены ВЫШЕ (до первого goto done) */
  if (o->medium && o->nrep > 0 && o->repof != NULL) {
    med_sub = (uint8_t *)calloc((size_t)nt, 1);
    med_box = (double *)calloc((size_t)o->nrep * 6, sizeof *med_box);
    med_an = (double *)calloc((size_t)o->nrep * 3, sizeof *med_an);
    med_v = (double *)calloc((size_t)o->nrep, sizeof *med_v);
    med_rho = (double *)calloc((size_t)o->nrep, sizeof *med_rho);
    med_asum = (double *)calloc((size_t)o->nrep, sizeof *med_asum);
    med_nsub = (int32_t *)calloc((size_t)o->nrep, sizeof *med_nsub);
    med_koff = (int32_t *)calloc((size_t)o->nrep + 1, sizeof *med_koff);
    med_kmem = (int32_t *)malloc((size_t)nt * sizeof *med_kmem);
    med_stamp = (int64_t *)calloc((size_t)o->nrep, sizeof *med_stamp);
    med_pend = (double *)calloc((size_t)o->nrep, sizeof *med_pend);
    if (!med_sub || !med_box || !med_an || !med_v || !med_rho || !med_asum || !med_nsub ||
        !med_koff || !med_kmem || !med_stamp) {
      rc = 2; /* fail closed: мир без среды (память) */
      goto done;
    }
  }
  /* §911-6: прибор покрытия (HZ_COVDBG=1) — буферы владельца hz_sw_run */
  if (o->mode == 3 && sw_covdbg()) {
    cov = (double *)calloc((size_t)nt, sizeof *cov);
    covn = (uint32_t *)calloc((size_t)nt, sizeof *covn);
    cov_lit = (double *)calloc((size_t)nt, sizeof *cov_lit);
    cov_e = (double *)calloc((size_t)nt, sizeof *cov_e);
  }
  if (!Ed || !Eprev || !order || !vindex || !bits || !walks || (o->mode == 3 && (!fld || !fldok)) ||
      (o->mode == 3 && (!Linmax_prev || !Linmax_cur)) || (cov && !covn)) {
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
    pstamp = (int64_t *)calloc((size_t)nt + (size_t)(o->nrep > 0 ? o->nrep : 0),
                               sizeof *pstamp); /* §924: +штампы представителей */
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
    leg_pstamp = (int64_t *)calloc((size_t)nt + (size_t)(o->nrep > 0 ? o->nrep : 0),
                                   sizeof *leg_pstamp); /* §924 */
    if (!pstamp || !bstart || !leg_pstamp) {
      rc = 2;
      goto done;
    }
    if (o->path) /* §865/раунды 13+15: кэш «кусок × трубка» — только path=1 */
      pc = (pc_rec *)calloc((size_t)nt, sizeof *pc);
    /* нет памяти → кэш остаётся NULL, path_collect_list идёт прежним путём
     * (fail closed: без новых массивов корректность не хуже раунда 12) */
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
  /* §932-Б/Ш9-ЗОНД + Ш12: клеточный приём — аккумулятор на лист + суммы
   * площадей детей клетки (раздача ced/ΣA в конце итерации);
   * объявлены ВЕРХом (goto done до этого блока — анализатор).
   * Ш12: cellc=2 добавляет cabs (радианс входа клетки для доли переизлучения).
   * Аллокация НА МЕСТЕ (не обёрткой): у gcc-analyzer класс FP на ёмкость↔
   * счётчик через обёртки выделения (CLAUDE.md, п. 2/S4). */
  if (o->cellc && bstart != NULL) {
    ced = (double *)calloc((size_t)py->nleaf, sizeof *ced);
    casum = (double *)calloc((size_t)py->nleaf, sizeof *casum);
    if (o->cellc == 2) cabs = (double *)calloc((size_t)py->nleaf, sizeof *cabs);
    if (!ced || !casum || (o->cellc == 2 && !cabs)) { /* fail closed: приём выкл */
      free(ced);
      free(casum);
      free(cabs);
      ced = NULL;
      casum = NULL;
      cabs = NULL;
    } else {
      for (int32_t c2 = 0; c2 < py->nleaf; c2++) {
        double s2 = 0.0;
        for (int32_t q2 = bstart[c2]; q2 < bstart[c2 + 1]; q2++)
          s2 += area[bpids[q2]];
        casum[c2] = s2;
      }
    }
  }

  for (it = 0; it < o->iters; it++) {
    double emitted = 0, absorbed = 0, lost = 0, recycled = 0, e_sum = 0, area_sum = 0;
    double lh_cap = 1e300, lh_seen = 0.0;                     /* §893 */
    double dep_main_it = 0, dep_leg_it = 0, hop_spawn_it = 0; /* §896-3 аудит */
    int32_t p;
    memset(Ed, 0, (size_t)nt * sizeof *Ed);
    for (p = 0; p < nt; p++)
      Eprev[p] = py->pcs[p].e;
    if (rep_ep != NULL) { /* §924: агрегат Eprev детей по представителям */
      for (int32_t r = 0; r < o->nrep; r++) {
        double s = 0.0;
        for (int32_t k = o->rep_koff[r]; k < o->rep_koff[r + 1]; k++) {
          int32_t kid = o->rep_kmem[k];
          s += area[kid] * Eprev[kid];
        }
        rep_ep[r] = o->rep_area[r] > 0.0 ? s / o->rep_area[r] : 0.0;
      }
    }
    if (o->mode == 3) { /* §893-b: swap пер-кускового максимума пришедшего радианса.
                         * Условие обязательно: Linmax_* выделяются ТОЛЬКО в mode 3
                         * (см. аллокацию выше), а без него memset(NULL, 0, nt·8) —
                         * настоящий null-deref для любого mode≠3
                         * (clang-analyzer-core.NonNullParamChecker, 07-10). */
      double *swp = Linmax_prev;
      Linmax_prev = Linmax_cur;
      Linmax_cur = swp;
      memset(Linmax_cur, 0, (size_t)nt * sizeof *Linmax_cur);
    }
    /* §893: Lh-ограничитель (принцип максимума радианса, §735): кап итерации
     * = LH_GROWTH·(макс радианс прошлой итерации). Рост ≤ +25 %/итерацию;
     * при альбедо<1 физическая сходимость даёт запас, расходимость — нет */
    {
      double lmax = 0.0;
      for (p = 0; p < nt; p++) {
        double rho_p = o->rho < 0 ? kd[p] : o->rho;
        double le_p = o->le + (o->lep ? o->lep[p] : 0.0);
        double lh = le_p + rho_p * Eprev[p] / (2.0 * M_PI);
        if (lh > lmax) lmax = lh;
      }
      lh_cap = LH_GROWTH * lmax;
      lh_seen = 0.0;
      /* §896-5: бюджет ног задаётся в fc (HZ_SW_LEG_BUDGET), а не локальной
       * переменной: локальное присваивание сюда не доходило (ловля 07-10) */
    }
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
        fc.lh_cap = lh_cap;
        fc.cov = cov; /* §911-6: NULL без HZ_COVDBG — прибора нет */
        fc.covn = covn;
        fc.cov_lit = cov_lit;
        fc.cov_e = cov_e;
        fc.Linmax_prev = Linmax_prev;
        fc.Linmax_cur = Linmax_cur;
        fc.py = py;
        fc.o = o;
        fc.area = area;
        fc.nrm = nrm;
        fc.kd = kd;
        fc.Eprev = Eprev;
        fc.Ed = Ed;
        fc.omcur[0] = om[0]; /* §873: хопы пишут сюда, таблица — read-only */
        fc.omcur[1] = om[1];
        fc.omcur[2] = om[2];
        fc.om = fc.omcur;
        { /* §928: ω-бин = уровень μ квадратуры, схема HZ_KCAL_NW×1
           * (таб. порядок i*nmu+m → imu = d % nmu); легаси 6/26 — бин 0 */
          int nmu2 = o->ndirs > 100 ? o->ndirs % 100 : 0;
          int imu2 = nmu2 > 0 ? d % nmu2 : 0;
          fc.dir_bin = nmu2 > 0 ? imu2 * HZ_KCAL_NW / nmu2 : 0;
          if (fc.dir_bin > HZ_KCAL_NW - 1) fc.dir_bin = HZ_KCAL_NW - 1;
        }
        fc.knum = (knum != NULL && it == o->iters - 1) ? knum : NULL; /* §928 */
        fc.kden = (fc.knum != NULL) ? kden : NULL;
        fc.kfloor = (const uint8_t *)lpflo; /* §928: NULL до первой конверсии */
        fc.w_d = w_d;
        fc.le = o->le;
        fc.lep = o->lep; /* А1576: per-piece эмиссия (NULL — прежний мир) */
        fc.noprop = o->noprop;
        fc.tau0 = o->tau0;
        fc.pstamp = pstamp;
        fc.leg_pstamp = leg_pstamp;
        fc.leg_budget = HZ_SW_LEG_BUDGET; /* §896-5: явное состояние ног (0 — выключены) */
        fc.pc = pc;
        fc.bstart = bstart;
        fc.bpids = bpids;
        fc.cellc = o->cellc; /* Ш12: режим клеточного приёма (0/1/2) */
        fc.ced = ced;        /* Ш9-зонд + Ш12: клеточный аккумулятор (NULL — выкл) */
        fc.cabs = cabs;      /* Ш12: радианс входа клетки (только cellc=2) */
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
        fc.lvls = lvls_ready ? lvls : NULL;             /* §923: до первой сборки — прежний мир */
        fc.leaf_loff = lvls_ready ? lvls_loff : NULL;   /* §923 */
        fc.leaf_lpids = lvls_ready ? lvls_lpids : NULL; /* §923 */
        /* §924: представители активны только вместе с раскладкой */
        fc.reps_on = (lvls_ready && rep_want != NULL);
        fc.nt0 = nt;
        fc.nrep = o->nrep;
        fc.rep_koff = o->rep_koff;
        fc.rep_kmem = o->rep_kmem;
        fc.rep_area = o->rep_area;
        fc.rep_rho = o->rep_rho;
        fc.rep_le = o->rep_le;
        fc.rep_atri = o->rep_atri;
        fc.rep_nrm = o->rep_nrm;
        fc.rep_ep = rep_ep;
        fc.ag = &ag;
        for (int32_t l2 = 0; l2 < py->nlev; l2++)
          fc.noff[l2] = noff[l2];
        fc.cbase = (it == 0); /* прибор А1569 — геометрия статична, хватит раза */
        fc.cfront = (it == 0);
        fc.hi[0] = o->domhi[0]; /* меш-граница, не сеточная (§852) */
        fc.hi[1] = o->domhi[1];
        fc.hi[2] = o->domhi[2];
        if (o->travel) { /* §932-А: Ldom — диагональ меша из py->lo × domhi,
                          * обе в одних координатах, >0 всегда (аудит-4) */
          double dx = o->domhi[0] - py->lo[0], dy = o->domhi[1] - py->lo[1],
                 dz = o->domhi[2] - py->lo[2];
          fc.ldom = sqrt(dx * dx + dy * dy + dz * dz);
          if (!(fc.ldom > 0.0)) fc.ldom = 1.0; /* вырожденный меш — страховка */
        }
        if (o->medium) { /* §932-Б: среда — указатели на агрегаты владельца */
          fc.med_on = med_ready;
          fc.med_box = med_box;
          fc.med_an = med_an;
          fc.med_v = med_v;
          fc.med_rho = med_rho;
          fc.med_asum = med_asum;
          fc.med_nsub = med_nsub;
          fc.med_koff = med_koff;
          fc.med_kmem = med_kmem;
          fc.med_stamp = med_stamp;
          fc.med_pend = med_pend;
        }
        front_dir(&fc, stampv, om);
        st->hop_lost += fc.hop_lost_sum;
        st->hops += fc.hop_hops_sum;
        st->hop_thr += fc.hop_lost_thr;
        st->hop_cap += fc.hop_lost_cap;
        st->row_dep += fc.row_dep;
        if (fc.lh_seen > lh_seen) lh_seen = fc.lh_seen;
        hop_spawn_it += fc.hop_spawn_e;
        st->dep_leg += fc.dep_leg;
        st->dep_main += fc.dep_main;
        dep_main_it += fc.dep_main;
        dep_leg_it += fc.dep_leg;
        st->dep_cnt += fc.dep_cnt;
        st->lin_pos += fc.lin_pos;
        st->lin_zero += fc.lin_zero;
        st->hop_cap += fc.hop_lost_cap;
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
        if (cov) { /* §911-6: отчёт покрытия направления */
          double num = 0, den = 0, numl = 0, nume = 0;
          double fuMu = 0, fuE = 0; /* §911-13-3: полные суммы по ВСЕМ кускам */
          double wr[5];
          int32_t wp[5];
          int32_t pi; /* индекс куска в отчёте: имя НЕ p — внешний p занят
                       * циклом кусков (была тень, -Wshadow) */
          int ntouch = 0, nunder = 0, nover = 0, nw = 0, q;
          for (q = 0; q < 5; q++) {
            wr[q] = 1.0;
            wp[q] = -1;
          }
          for (pi = 0; pi < nt; pi++) {
            double mu = fabs(om[0] * nrm[3 * (int64_t)pi] + om[1] * nrm[3 * (int64_t)pi + 1] +
                             om[2] * nrm[3 * (int64_t)pi + 2]);
            double rr;
            fuMu += area[pi] * mu; /* §911-13-3: полный Σa·μ без фильтра touched */
            fuE += area[pi] * cov_e[pi];
            if (covn[pi] == 0 || mu < 1e-12) continue;
            rr = cov[pi] / mu;
            num += area[pi] * cov[pi];
            numl += area[pi] * cov_lit[pi];
            nume += area[pi] * cov_e[pi];
            den += area[pi] * mu;
            ntouch++;
            if (rr < 0.95) {
              nunder++;
            } else if (rr > 1.05) {
              nover++;
            }
            if (fabs(rr - 1.0) > 0.02 && nw < 5) {
              wr[nw] = rr;
              wp[nw] = pi;
              nw++;
            }
          }
          fprintf(stderr,
                  "COVDBG d=%d nd=%d om=(%.4f,%.4f,%.4f) areacov=%.5f litcov=%.5f "
                  "raw(mu=%.4f lit=%.4f E=%.4f) full(mu=%.3f E=%.3f) ndep=%lld touched=%d "
                  "under=%d over=%d lhcaps=%lld worst:",
                  d, nd, om[0], om[1], om[2], den > 0.0 ? num / den : 0.0,
                  den > 0.0 ? numl / den : 0.0, den, numl, nume, fuMu, fuE, (long long)fc.ndep,
                  ntouch, nunder, nover, (long long)fc.lh_cap_hits);
          for (q = 0; q < nw; q++)
            fprintf(stderr, " p%d=%.3f(a=%.3g,n=%u)", wp[q], wr[q], area[wp[q]], covn[wp[q]]);
          fprintf(stderr, "\n");
          if (d == 0 && it == 0) { /* доля площади кусков ниже depden-floor */
            double sa = 0.0, sf = 0.0, floor_a = csec * FRONT_DEP_MIN_SHARE;
            for (pi = 0; pi < nt; pi++) {
              sa += area[pi];
              if (area[pi] < floor_a) sf += area[pi];
            }
            fprintf(stderr, "COVDBG small-area frac (a<%.3g): %.4f (nt=%d)\n", floor_a,
                    sa > 0.0 ? sf / sa : 0.0, nt);
          }
          memset(cov, 0, (size_t)nt * sizeof *cov);
          memset(covn, 0, (size_t)nt * sizeof *covn);
          memset(cov_lit, 0, (size_t)nt * sizeof *cov_lit);
          memset(cov_e, 0, (size_t)nt * sizeof *cov_e);
        }
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
    if (ced != NULL) { /* §932-Б/Ш9-ЗОНД + Ш12: раздача клеточного приёма
                        * детям по площадям: Σ Ed·A = ced (тождество).
                        * cellc=1: ced — Φ_входа (RAW-зонд). cellc=2: ced —
                        * ИЗВЛЕЧЁННЫЙ из трубки поток Φ_вх·(1−e^{−Λ}) */
      double cedtot = 0.0;
      int64_t cedn = 0;
      for (int32_t c2 = 0; c2 < py->nleaf; c2++) {
        double q2 = ced[c2];
        if (q2 > 0.0 && casum[c2] > 0.0) {
          double e2 = q2 / casum[c2];
          for (int32_t q3 = bstart[c2]; q3 < bstart[c2 + 1]; q3++) {
            Ed[bpids[q3]] += e2;
            cedn++;
          }
          cedtot += q2;
          ced[c2] = 0.0;
        }
      }
      g932c_dep += cedtot;
      g932c_n += cedn;
    }
    { /* §932-Б/Ш9: ИНВАРИАНТ ПОЗИТИВНОСТИ — у незамкнутого куска E>0
       * после ≥1 итерации (диффузный фон строго положителен; точные
       * нули = дыра в приёмнике — доказано Ш8-4: 96.2% нулей на UE).
       * Прибор: доля кусков с E==0; канон-гейт = 0 для открытых сцен. */
      int64_t nz = 0;
      for (p = 0; p < nt; p++) {
        py->pcs[p].e = (float)Ed[p];
        e_sum += Ed[p] * area[p];
        area_sum += area[p];
        if (!(Ed[p] > 0.0)) nz++;
      }
      st->zero_e = nz;
    }
    st->e_avg = e_sum / area_sum;
    if (getenv("HZ_AUDIT"))
      fprintf(stderr,
              "AUDIT it=%d E_avg=%.4g emitted=%.1f absorbed=%.1f lost=%.1f dep_main=%.1f "
              "dep_leg=%.1f hop_spawn=%.1f hop_lost=%.1f\n",
              it, st->e_avg, emitted, absorbed, lost, dep_main_it, dep_leg_it, hop_spawn_it,
              st->hop_lost);
    st->emitted = emitted;
    st->recycled = recycled;
    st->absorbed = absorbed;
    st->lost = lost;
    if (o->mode != 3) st->traffic = (int64_t)nd * ((int64_t)nt * 36 + (int64_t)walks[0].n * 40);
    if (e_hist) e_hist[it] = st->e_avg;
    { /* §924-дых: E по итерациям */
      static int itdbg = -1;
      if (itdbg < 0) itdbg = getenv("HZ_DBG924") != NULL;
      if (itdbg) {
        fprintf(stderr,
                "DBG924it it=%d e_avg=%.4f absorbed=%.5g emitted=%.5g rep_hits=%lld SumLin=%lld\n",
                it, st->e_avg, st->absorbed, st->emitted, (long long)g924_hits,
                (long long)g924_lin);
        if (g924_kclamp)
          fprintf(stderr, "DBG924k clamp=%lld vis=%lld geohit=%lld\n", (long long)g924_kclamp,
                  (long long)g924_vis, (long long)g924_geohit);
        g924_kclamp = 0;
        g924_vis = 0;
        g924_geohit = 0;
      }
      g924_hits = 0;
      g924_lin = 0;
    }
    if (o->mode == 3 && o->lpacc && o->lpapply) { /* §862: этаж = целая часть */
      int32_t nup = 0;
      if (!lpflo) {
        lpflo = (uint8_t *)calloc(
            (size_t)nt,
            sizeof *lpflo); /* 0: первая
                             * конверсия честно посчитает смены относительно fine-этажа */
        if (!lpflo) {
          rc = 2;
          goto done;
        }
        lpflo_tri = (uint8_t *)calloc((size_t)nt, sizeof *lpflo_tri); /* §922/А1632 */
        if (!lpflo_tri) {
          rc = 2;
          goto done;
        }
      }
      { /* §866: этаж = min(floor(lpacc), lpceil), отрицательный вклад — 0
         * (раньше (uint8_t) от отрицательного double — UB); приборы смены.
         * §868: lpnorm — относительный энерговклад fl = cdelta·lpacc/mean */
        int64_t nchg = 0;
        int32_t h[16] = {0}, q;
        double mean = 0.0;
        if (o->lpnorm)
          for (p = 0; p < nt; p++)
            mean += o->lpacc[p];
        if (o->lpnorm && nt > 0) mean /= (double)nt;
        for (p = 0; p < nt; p++) {
          double fl = o->lpacc[p];
          int32_t fi;
          if (o->lpnorm) { /* §868: база + относительный вклад, монотонно */
            fl = mean > 0.0 ? (double)o->lpbase + o->cdelta * fl / mean : (double)o->lpbase;
            fi = (int32_t)fl;
            if (fi < o->lpbase) fi = o->lpbase;
          } else {
            fi = fl >= 255.0 ? 255 : (fl < 0.0 ? 0 : (int32_t)fl);
          }
          /* §932-Б/Ш9: ХОЛОДНЫЙ СТАРТ — микрокусок (area < colda) стартует
           * на этаже coldf: реп-носитель ловит депозиты линий и раздаёт
           * детям по площадям (§924-канон); иначе яйцо-курица: без
           * депозитов этаж не набирается никогда (Ш8-4, 96.2% нулей).
           * Требует rep'а на этом этаже — иначе кусочно (fail-closed). */
          if (o->cold_f > 0 && o->cold_a > 0.0 && area[p] < o->cold_a && fi < o->cold_f) {
            int32_t fc2 = o->cold_f;
            if (fc2 > (int32_t)py->nlev) fc2 = (int32_t)py->nlev;
            if (fc2 > o->rep_nlev) fc2 = o->rep_nlev;
            if (o->repof != NULL && fc2 >= 1 &&
                o->repof[((size_t)fc2 - 1u) * (size_t)nt + (size_t)p] >= 0)
              fi = fc2;
          }
          if (o->lpceil > 0 && fi > o->lpceil) fi = o->lpceil;
          if (lpflo[p] != (uint8_t)fi) nchg++;
          lpflo[p] = (uint8_t)fi;
          lpflo_tri[py->pcs[p].tri] = (uint8_t)fi; /* §922/А1632: конверсия в tri-пространство */
          h[fi > 15 ? 15 : fi]++;
        }
        if (nchg > st->flochg_max) st->flochg_max = nchg;
        st->flochg_sum += nchg;
        for (q = 0; q < 16; q++)
          st->flhist[q] = h[q];
      }
      if (hz_pyr_set_lp(py, lpflo_tri, &nup) != 0) { /* §922/А1632: tri-пространство */
        rc = 2;
        goto done;
      }
      if (med_sub != NULL) { /* §932-Б/А1737: замещение — энергетическая
                              * малость: lpacc_e[p] < τ·mean(lpacc_e);
                              * зеркала/эмиттеры НЕ замещаются (А1732);
                              * агрегаты регионов пересчитываются вместе
                              * со списками (А1711) */
        double mean_e = 0.0;
        int32_t p2, nsub_tot = 0;
        for (p2 = 0; p2 < nt; p2++)
          mean_e += o->lpacc_e[p2];
        if (nt > 0) mean_e /= (double)nt;
        memset(med_nsub, 0, (size_t)o->nrep * sizeof *med_nsub);
        memset(med_an, 0, (size_t)o->nrep * 3 * sizeof *med_an);
        memset(med_asum, 0, (size_t)o->nrep * sizeof *med_asum);
        memset(med_rho, 0, (size_t)o->nrep * sizeof *med_rho);
        for (int32_t ar = 0; ar < o->nrep * 6; ar++) {
          med_box[ar] = ar % 6 < 3 ? 1e30 : -1e30; /* пустой union-bbox */
        }
        memset(med_koff, 0, (size_t)(o->nrep + 1) * sizeof *med_koff);
        for (p2 = 0; p2 < nt; p2++) {
          int32_t f = lpflo[p2], r = -1;
          if (f < 1 || f > o->rep_nlev) {
            med_sub[p2] = 0;
            continue;
          }
          r = o->repof[((size_t)f - 1u) * (size_t)nt + (size_t)p2];
          med_sub[p2] = (int32_t)(r >= 0 && o->lpacc_e[p2] < o->med_tau * mean_e &&
                                  !(o->ks && o->ks[p2] > 0.0) && !(o->lep && o->lep[p2] > 0.0))
                            ? 1
                            : 0;
          if (!med_sub[p2]) continue;
          { /* агрегаты: ΣA, ΣA·n, ρ-среднее, bbox-union (tribox детей) */
            int32_t tr = py->pcs[p2].tri;
            double a2 = area[p2];
            med_nsub[r]++;
            med_asum[r] += a2;
            med_rho[r] += a2 * (o->rho >= 0.0 ? o->rho : kd[p2]);
            for (int ax = 0; ax < 3; ax++) {
              med_an[3 * r + ax] += a2 * nrm[3 * (int64_t)p2 + ax];
              double lo2 = o->tribox[6 * (int64_t)tr + ax],
                     hi2 = o->tribox[6 * (int64_t)tr + 3 + ax];
              if (lo2 < med_box[6 * r + ax]) med_box[6 * r + ax] = lo2;
              if (hi2 > med_box[6 * r + 3 + ax]) med_box[6 * r + 3 + ax] = hi2;
            }
          }
        }
        for (int32_t r = 0; r < o->nrep; r++) { /* CSR: counting-sort детей */
          med_koff[r + 1] = med_koff[r] + med_nsub[r];
          if (med_asum[r] > 0.0)
            med_rho[r] /= med_asum[r];
          else
            med_rho[r] = 0.0;
          if (med_nsub[r] > 0) {
            double vx = med_box[6 * (int64_t)r + 3] - med_box[6 * (int64_t)r],
                   vy = med_box[6 * (int64_t)r + 4] - med_box[6 * (int64_t)r + 1],
                   vz = med_box[6 * (int64_t)r + 5] - med_box[6 * (int64_t)r + 2];
            med_v[r] = (vx > 0.0 ? vx : 0.0) * (vy > 0.0 ? vy : 0.0) * (vz > 0.0 ? vz : 0.0);
          }
        }
        { /* заполнение kmem по koff-курсорам (koff портится → копия) */
          int32_t *cur = (int32_t *)malloc((size_t)(o->nrep + 1) * sizeof *cur);
          if (cur != NULL) {
            memcpy(cur, med_koff, (size_t)(o->nrep + 1) * sizeof *cur);
            for (p2 = 0; p2 < nt; p2++) {
              if (!med_sub[p2]) continue;
              int32_t f = lpflo[p2];
              int32_t r = o->repof[((size_t)f - 1u) * (size_t)nt + (size_t)p2];
              med_kmem[cur[r]++] = p2;
            }
            for (int32_t r = 0; r < o->nrep; r++)
              nsub_tot += med_nsub[r];
            free(cur);
          }
        }
        med_ready = 1;
        { /* Ш3: какие регионы заместились (диагностика стенда) */
          int32_t r2, nmed = 0;
          for (r2 = 0; r2 < o->nrep; r2++)
            if (med_nsub[r2] > 0) nmed++;
          fprintf(stderr,
                  "§932-Б: замещено %d/%d кусков (tau=%.3f, mean_e=%.4g), активаций среды будет "
                  "с итерации далее; регионов с заменой=%d",
                  nsub_tot, (int)nt, o->med_tau, mean_e, nmed);
          for (r2 = 0; r2 < o->nrep && nmed > 0; r2++)
            if (med_nsub[r2] > 0)
              fprintf(stderr, "%s r%d: n=%d V=%.3g an=(%.3g,%.3g,%.3g)", r2 > 0 ? ";" : ":",
                      (int)r2, (int)med_nsub[r2], med_v[r2], med_an[3 * (int64_t)r2],
                      med_an[3 * (int64_t)r2 + 1], med_an[3 * (int64_t)r2 + 2]);
          fprintf(stderr, "\n");
        }
      }
      if (lvls) { /* §923: перекладка кусков по уровням нового этажа */
        int brc = sw_levels_build(py, o, lpflo, nt, lvls, &lvl_pool, &lvl_pool_cap, &lvls_lpids,
                                  &lvls_lpids_cap, lvls_loff, bstart, bpids, rep_want, med_sub);
        if (brc == 0) lvls_ready = 1;
        if (brc != 0) {
          /* отказ памяти: fail closed — разбор раскладки, старый мир */
          free(lvls);
          lvls = NULL;
          free(lvl_off_all);
          lvl_off_all = NULL;
          free(lvl_finer_all);
          lvl_finer_all = NULL;
          free(lvls_loff);
          lvls_loff = NULL;
        }
      }
    }
  }

  if (knum != NULL) { /* §928: дамп числителя/знаменателя калибровки */
    hz_kcal_hdr kh;
    char p2[512];
    memset(&kh, 0, sizeof kh);
    memcpy(kh.magic, "KCAL928", sizeof kh.magic);
    kh.nrep = o->nrep;
    kh.nw = HZ_KCAL_NW;
    if (snprintf(p2, sizeof p2, "%s.num", kcalpath) < (int)sizeof p2) {
      FILE *f = fopen(p2, "wb");
      if (f != NULL) {
        fwrite(&kh, sizeof kh, 1, f);
        fwrite(knum, sizeof *knum, (size_t)o->nrep * HZ_KCAL_NW, f);
        fclose(f);
      } else
        fprintf(stderr, "§928: не открылся %s\n", p2);
    }
    if (snprintf(p2, sizeof p2, "%s.den", kcalpath) < (int)sizeof p2) {
      FILE *f = fopen(p2, "wb");
      if (f != NULL) {
        fwrite(&kh, sizeof kh, 1, f);
        fwrite(kden, sizeof *kden, (size_t)o->nrep * HZ_KCAL_NW, f);
        fclose(f);
      } else
        fprintf(stderr, "§928: не открылся %s\n", p2);
    }
    fprintf(stderr, "§928: дамп калибровки %s.{num,den} (nrep=%d nw=%d)\n", kcalpath, (int)o->nrep,
            HZ_KCAL_NW);
  }
  if (o->rep_zone > 0.0 && g931_rh > 0) /* §931/А1668: ШОВ — доля rep-хитов,
                                         * чья точка пересечения ближе зоны;
                                         * >5% — усиливать предикат bbox-тестом.
                                         * (0 при it=1 — этажные списки reps
                                         * строятся после итерации 0) */
    fprintf(stderr, "§931 шов: rep-хитов=%lld, точка ближе зоны=%lld (%.2f%%)\n",
            (long long)g931_rh, (long long)g931_seam, 100.0 * (double)g931_seam / (double)g931_rh);
  g931_rh = 0;
  g931_seam = 0;
  if (o->travel) { /* §932-А: калибровочная печать (А1716/А1724) — ⟨f⟩, ⟨Lv⟩,
                    * гистограмма; стартовый tvc ~ Ldom/⟨Lv⟩ (ожидание 10–50) */
    double dx = o->domhi[0] - py->lo[0], dy = o->domhi[1] - py->lo[1], dz = o->domhi[2] - py->lo[2];
    double ldom = sqrt(dx * dx + dy * dy + dz * dz);
    if (g932_n > 0)
      fprintf(stderr,
              "§932-travel: событий с пробегом=%lld <f>=%.4f <Lv>=%.3f Ldom=%.3f "
              "стартовый tvc~Ldom/<Lv>=%.1f%s\n",
              (long long)g932_n, g932_fsum / (double)g932_n, g932_lvsum / (double)g932_n, ldom,
              ldom / (g932_lvsum / (double)g932_n), g932_nomax ? " [НК HZ_TVNOMAX]" : "");
    else
      fprintf(stderr, "§932-travel: событий с пробегом нет (t_last<0 всюду)\n");
    fprintf(stderr, "§932-travel: гистограмма f=Lv/Ldom (шаг 1/16):");
    for (int q2 = 0; q2 < 16; q2++)
      fprintf(stderr, " %d:%lld", q2, (long long)g932_hist[q2]);
    fprintf(stderr, "\n");
  }
  if (o->medium && g932b_ev >= 0) /* §932-Б: приборы среды (баланс: Σдеп =
                                   * (1−ρ̄)·Σпогл по построению банка) */
    fprintf(stderr,
            "§932-Б: активаций=%lld, σℓ>3=%lld, поглощено=%.6g, депонировано детям=%.6g "
            "(доля=%.15f)\n§932-Б-Ш3: σℓ mean=%.4f max=%.4f (n=%lld), факт-Ed=%.6g "
            "(отл./счёт=%.3g), вложенных входов=%lld%s%s%s\n",
            (long long)g932b_ev, (long long)g932b_big, g932b_abs, g932b_dep,
            g932b_abs > 0.0 ? g932b_dep / g932b_abs : 0.0,
            g932b_sact > 0 ? g932b_ssum / (double)g932b_sact : 0.0, g932b_smax,
            (long long)g932b_sact, g932b_edf, g932b_dep > 0.0 ? g932b_edf / g932b_dep : 0.0,
            (long long)g932b_nest, sw_med_now ? " [MEDNOW]" : "", sw_med_t1 ? " [MEDT1-ЗОНД]" : "",
            sw_med_scr ? " [MEDSCRAMBLE-НК]" : "");
  if (o->cellc == 1) /* §932-Б/Ш9: RAW-зонд (НЕконсервативен, НК Ш12) */
    fprintf(stderr, "§932-Б-Ш9: клеточный приём RAW-зонд: Σ=%.6g, раздач=%lld\n", g932c_dep,
            (long long)g932c_n);
  if (o->cellc == 2) /* §932-Б/Ш12: консервативная схема — Σ извлечённого
                      * из трубок потока (= Σ депозитов клеток тождество) */
    fprintf(stderr,
            "§932-Б-Ш12: клеточный приём консервативный: Σизвлечено=%.6g, раздач=%lld, "
            "изъятий=%lld, maxf=%.17g, f≥1−1e−12: %lld\n",
            g932c_dep, (long long)g932c_n, (long long)g932c_vis, g932c_fmax, (long long)g932c_ge1);

done:
  free(med_sub);       /* §932-Б */
  free(med_box);       /* §932-Б */
  free(med_an);        /* §932-Б */
  free(med_v);         /* §932-Б */
  free(med_rho);       /* §932-Б */
  free(med_asum);      /* §932-Б */
  free(med_nsub);      /* §932-Б */
  free(med_koff);      /* §932-Б */
  free(med_kmem);      /* §932-Б */
  free(med_stamp);     /* §932-Б */
  free(med_pend);      /* §932-Б (ASAN-ловля Ш3: не освобождался с Ш2) */
  free(ced);           /* §932-Б/Ш9-зонд */
  free(casum);         /* §932-Б/Ш9-зонд */
  free(cabs);          /* §932-Б/Ш12 */
  free(lpflo);         /* §862 */
  free(lpflo_tri);     /* §922/А1632 */
  free(lvls);          /* §923 */
  free(lvl_off_all);   /* §923: плоские массивы уровней */
  free(lvl_finer_all); /* §923 */
  free(lvl_pool);      /* §923 */
  free(lvls_loff);     /* §923 */
  free(lvls_lpids);    /* §923 */
  free(rep_want);      /* §924 */
  free(rep_ep);        /* §924 */
  free(knum);          /* §928 */
  free(kden);          /* §928 */
  free(fld);           /* §864/Б1 */
  free(fldok);
  free(Ed);
  free(Linmax_prev);
  free(Linmax_cur);
  free(Eprev);
  free(cov); /* §911-6 */
  free(covn);
  free(cov_lit);
  free(cov_e);
  free(order);
  free(vindex);
  free(bits);
  free(tab);
  free(pstamp);     /* §852/А1564 */
  free(leg_pstamp); /* §894-c-5 */
  free(pc);         /* §865/раунды 13+15 */
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

/* §842: Σ по ВЫХОДНЫМ направлениям (n·ω > 0) w·|ω·n| — облучённость куска
 * в модели слоя при тёмном окружении. Для точного решения фикстуры.
 *
 * Здесь был ещё hz_sw_dirsum (Σ по ВСЕМ направлениям) — УДАЛЁН 07-10: у него не
 * было ни одного вызова ни в живом слое, ни в archive/ (двойник ниже зовётся
 * tools/swee3.c). Вернуть из git-истории, если понадобится нормализация весов. */
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
