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

int hz_sw_run(hz_pyr *py, int32_t nt, const double *area, const double *nrm, const double *kd,
              const hz_sw_opts *o, hz_sw_stat *st, double *e_hist) {
  double *Ed = NULL, *Eprev = NULL;
  int32_t *order = NULL, *vindex = NULL;
  uint64_t *bits = NULL;
  hz_sw_walk *walks = NULL;
  hz_sw_dir *tab = NULL;
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

  if (o->mode == 1) {
    for (d = 0; d < nd; d++) {
      int sg[3];
      for (int ax = 0; ax < 3; ax++)
        sg[ax] = (tab[d].om[ax] > 0.0) - (tab[d].om[ax] < 0.0);
      rc = sw_build_walk(py, sg, tab[d].om, area, nrm, kd, &walks[d]);
      if (rc != 0) goto done;
    }
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
      double L = 0.0;
      int64_t prevline = -1;
      if (o->mode == 1) {
        const hz_sw_walk *w = &walks[d];
        int32_t k;
        for (k = 0; k < w->n; k++) {
          const hz_sw_wleaf *lf = &w->leaf[k];
          double Lsurf = 0.0;
          int64_t e, eend = lf->first + lf->npcs;
          int64_t pf = lf->first + 8 < nt ? lf->first + 8 : nt - 1;
          if (lf->line != prevline) { /* новая линия — тёмный вход */
            L = 0.0;
            prevline = lf->line;
          }
          if (o->noprop) L = 0.0;
          if (lf->cntot <= 0.0) continue; /* куски встык лучу — слоя нет */
          __builtin_prefetch(&w->leaf[k + 8 < w->n ? k + 8 : w->n - 1]);
          __builtin_prefetch(&Eprev[w->wpid[pf]]);
          for (e = lf->first; e < eend; e++) {
            p = w->wpid[e];
            double rho = o->rho < 0 ? w->kd[e] : o->rho;
            Lsurf += (o->le + rho * Eprev[p] / (2.0 * M_PI)) * w->wcn[e];
            Ed[p] += w_d * L * w->wn[e]; /* непрозрачный слой: инфлюкс целиком */
          }
          Lsurf /= lf->cntot;
          absorbed += w_d * L * csec;
          emitted += w_d * o->le * csec;
          recycled += w_d * (Lsurf - o->le) * csec;
          L = Lsurf; /* непрозрачный слой: луч гасится, остаётся переизлучение */
          if (o->noprop) L = 0.0;
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
