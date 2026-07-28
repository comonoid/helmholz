/* Планарная сегментация. Разбор и оговорки — в `poly_seg.h`. */

#include "poly_seg.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Уровней сетки хватает на 2^48 отношений размеров — при базовом шаге 2δ это
 * покрывает любую сцену, какую можно загрузить. */
#define SEG_MAXLEV 48
/* Проверка отпускает нарушителей; каждый проход отпускает хотя бы одного,
 * поэтому цикл конечен и без потолка. Потолок стоит затем, чтобы у худшего
 * случая была ЛИНЕЙНАЯ, а не квадратичная цена: после него плоскость
 * фиксируется плоскостью ЗАТРАВКИ, на которой участок валиден заведомо. */
#define SEG_MAXFIX 8
/* Циклический Якоби для 3×3: обмен порогом на фиксированное число проходов.
 * Порога нет сознательно (CLAUDE.md: никаких магических порогов), а 12 проходов
 * для 3×3 — заведомый запас: сходимость квадратичная. */
#define SEG_JACOBI 12
/* Раундов пересева. Каждый раунд заводит хотя бы один участок на отпущенных,
 * поэтому потолок — страховка от бесконечности, а не рабочее ограничение. */
#define SEG_MAXROUNDS 64

/* --- 3×3 симметричная: собственный вектор наименьшего значения ------------- */

static void eig3_min(const double c[6], double n[3]) {
  /* c = xx, yy, zz, xy, xz, yz */
  double a[3][3] = {{c[0], c[3], c[4]}, {c[3], c[1], c[5]}, {c[4], c[5], c[2]}};
  double v[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  for (int sweep = 0; sweep < SEG_JACOBI; sweep++)
    for (int p = 0; p < 2; p++)
      for (int q = p + 1; q < 3; q++) {
        double apq = a[p][q];
        if (!(apq < 0.0) && !(apq > 0.0)) continue; /* точный ноль — вращать нечего */
        double theta = (a[q][q] - a[p][p]) / (2.0 * apq);
        double t = (theta >= 0.0) ? 1.0 / (theta + sqrt(1.0 + theta * theta))
                                  : -1.0 / (-theta + sqrt(1.0 + theta * theta));
        double cs = 1.0 / sqrt(1.0 + t * t), sn = t * cs;
        for (int k = 0; k < 3; k++) {
          double akp = a[k][p], akq = a[k][q];
          a[k][p] = cs * akp - sn * akq;
          a[k][q] = sn * akp + cs * akq;
        }
        for (int k = 0; k < 3; k++) {
          double apk = a[p][k], aqk = a[q][k];
          a[p][k] = cs * apk - sn * aqk;
          a[q][k] = sn * apk + cs * aqk;
        }
        for (int k = 0; k < 3; k++) {
          double vkp = v[k][p], vkq = v[k][q];
          v[k][p] = cs * vkp - sn * vkq;
          v[k][q] = sn * vkp + cs * vkq;
        }
      }
  int im = 0;
  for (int i = 1; i < 3; i++)
    if (a[i][i] < a[im][im]) im = i;
  double nn = sqrt(v[0][im] * v[0][im] + v[1][im] * v[1][im] + v[2][im] * v[2][im]);
  if (!(nn > 0.0)) {
    n[0] = 0.0;
    n[1] = 0.0;
    n[2] = 1.0;
    return;
  }
  for (int k = 0; k < 3; k++)
    n[k] = v[k][im] / nn;
}

/* --- сетка поиска: иерархическая хеш-сетка --------------------------------- */

/* Уровень задаётся РАЗМЕРОМ треугольника: на своём уровне шаг ячейки не меньше
 * `габарит + 2δ`, поэтому и запись, и запрос трогают не больше восьми ячеек, а
 * число просмотренных ячеек от δ не зависит вовсе.
 *
 * СТОЛКНОВЕНИЕ КЛЮЧЕЙ БЕЗВРЕДНО, и это записано, чтобы не «чинить» его потом:
 * две разные ячейки, попавшие в один слот, дают ЛИШНИХ кандидатов, а каждый
 * кандидат всё равно проверяется точным тестом коробок. Ложных пропусков
 * столкновение не создаёт. */
typedef struct {
  uint64_t *key; /* 0 — слот пуст */
  int32_t *start;
  int32_t *idx;
  int64_t nslot, nent;
  uint64_t mask;
  double o[3], h0;
  int32_t maxlev;
  unsigned char *lev;
  unsigned char levused[SEG_MAXLEV];
} seg_grid;

static uint64_t ckey(int32_t lev, int64_t ix, int64_t iy, int64_t iz) {
  uint64_t k = (uint64_t)(uint32_t)(int32_t)ix * 0x9E3779B97F4A7C15ULL;
  k ^= (uint64_t)(uint32_t)(int32_t)iy * 0xC2B2AE3D27D4EB4FULL;
  k ^= (uint64_t)(uint32_t)(int32_t)iz * 0x165667B19E3779F9ULL;
  k ^= (uint64_t)(uint32_t)lev * 0x27D4EB2F165667C5ULL;
  /* ПЕРЕМЕШИВАНИЕ ОБЯЗАТЕЛЬНО (§5, урок 4): у произведений малых координат
   * старшие биты почти одинаковы, и хеш без него схлопывает сцену в 2–4 слота. */
  k ^= k >> 33;
  k *= 0xFF51AFD7ED558CCDULL;
  k ^= k >> 33;
  k *= 0xC4CEB9FE1A85EC53ULL;
  k ^= k >> 33;
  return k | 1ULL; /* 0 занят под «пусто» */
}

static int64_t slot_of(seg_grid *g, uint64_t k, int insert) {
  uint64_t i = k & g->mask;
  for (;;) {
    if (g->key[i] == k) return (int64_t)i;
    if (g->key[i] == 0) {
      if (!insert) return -1;
      g->key[i] = k;
      return (int64_t)i;
    }
    i = (i + 1) & g->mask;
  }
}

static void cell_range(const seg_grid *g, const double *bb, double pad, int32_t lev, int64_t c0[3],
                       int64_t c1[3]) {
  double h = ldexp(g->h0, lev);
  for (int a = 0; a < 3; a++) {
    c0[a] = (int64_t)floor((bb[a] - pad - g->o[a]) / h);
    c1[a] = (int64_t)floor((bb[3 + a] + pad - g->o[a]) / h);
  }
}

/* Конец ячейки — начало СЛЕДУЮЩЕГО слота: начала монотонны по индексу слота,
 * потому что пустой слот имеет нулевую длину. */
static void slot_span(const seg_grid *g, int64_t sl, int32_t *b, int32_t *e) {
  *b = g->start[sl];
  *e = (sl + 1 < g->nslot) ? g->start[sl + 1] : (int32_t)g->nent;
}

/* --- расстояние между коробками ≤ d ---------------------------------------- */

static int bb_near(const double *A, const double *B, double d) {
  for (int a = 0; a < 3; a++) {
    if (A[a] - B[3 + a] > d) return 0;
    if (B[a] - A[3 + a] > d) return 0;
  }
  return 1;
}

/* --- порядок затравок: ГЕОМЕТРИЧЕСКИЙ ключ, номер треугольника не участвует -- */

typedef struct {
  double area;
  double c[3];
  int32_t t;
} seg_key;

static int cmp_key(const void *pa, const void *pb) {
  const seg_key *a = pa, *b = pb;
  if (a->area > b->area) return -1;
  if (a->area < b->area) return 1;
  for (int i = 0; i < 3; i++) {
    if (a->c[i] < b->c[i]) return -1;
    if (a->c[i] > b->c[i]) return 1;
  }
  return 0;
}

/* --- подгонка плоскости ----------------------------------------------------- */

typedef struct {
  double org[3]; /* сдвиг координат = центр затравки: без него в S2 сокращаются знаки */
  double W;
  double S1[3];
  double S2[6];
  double G[3]; /* Σ площадь·нормаль — задаёт ЗНАК плоскости */
  double n[3], off;
} seg_fit;

static void fit_reset(seg_fit *f, const double org[3]) {
  memset(f, 0, sizeof *f);
  for (int a = 0; a < 3; a++)
    f->org[a] = org[a];
}

static void fit_add(seg_fit *f, const double p[3][3], const double gnr[3], double area) {
  double w = area / 3.0;
  for (int i = 0; i < 3; i++) {
    double q[3];
    for (int a = 0; a < 3; a++)
      q[a] = p[i][a] - f->org[a];
    f->W += w;
    for (int a = 0; a < 3; a++)
      f->S1[a] += w * q[a];
    f->S2[0] += w * q[0] * q[0];
    f->S2[1] += w * q[1] * q[1];
    f->S2[2] += w * q[2] * q[2];
    f->S2[3] += w * q[0] * q[1];
    f->S2[4] += w * q[0] * q[2];
    f->S2[5] += w * q[1] * q[2];
  }
  for (int a = 0; a < 3; a++)
    f->G[a] += area * gnr[a];
}

static void fit_solve(seg_fit *f) {
  if (!(f->W > 0.0)) return;
  double c[3];
  for (int a = 0; a < 3; a++)
    c[a] = f->S1[a] / f->W;
  double cov[6];
  cov[0] = f->S2[0] / f->W - c[0] * c[0];
  cov[1] = f->S2[1] / f->W - c[1] * c[1];
  cov[2] = f->S2[2] / f->W - c[2] * c[2];
  cov[3] = f->S2[3] / f->W - c[0] * c[1];
  cov[4] = f->S2[4] / f->W - c[0] * c[2];
  cov[5] = f->S2[5] / f->W - c[1] * c[2];
  double n[3];
  eig3_min(cov, n);
  if (n[0] * f->G[0] + n[1] * f->G[1] + n[2] * f->G[2] < 0.0)
    for (int a = 0; a < 3; a++)
      n[a] = -n[a];
  for (int a = 0; a < 3; a++)
    f->n[a] = n[a];
  f->off = n[0] * (c[0] + f->org[0]) + n[1] * (c[1] + f->org[1]) + n[2] * (c[2] + f->org[2]);
}

static double tri_dev(const double p[3][3], const double n[3], double off) {
  double d = 0.0;
  for (int i = 0; i < 3; i++) {
    double e = fabs(p[i][0] * n[0] + p[i][1] * n[1] + p[i][2] * n[2] - off);
    if (e > d) d = e;
  }
  return d;
}

void hz_seg_free(hz_pseglist *s) {
  free(s->label);
  free(s->seg);
  memset(s, 0, sizeof *s);
}

/* --- рабочее состояние: одна структура, чтобы освобождение было одним местом -- */

typedef struct {
  double *bb;  /* 6*nt */
  double *gn;  /* 3*nt */
  double *ar;  /* nt   */
  seg_key *ord;
  int32_t *queue;
  int32_t *rank;
  int32_t *xs;  /* nt+1: начала списков «мелкие соседи крупного» */
  int32_t *xl;  /* ncross */
  int32_t *tmp; /* nslot или nt — курсоры заполнения */
  hz_pseg *seg;
  seg_grid g;
} seg_work;

static void work_free(seg_work *w) {
  free(w->bb);
  free(w->gn);
  free(w->ar);
  free(w->ord);
  free(w->queue);
  free(w->rank);
  free(w->xs);
  free(w->xl);
  free(w->tmp);
  free(w->seg);
  free(w->g.key);
  free(w->g.start);
  free(w->g.idx);
  free(w->g.lev);
}

int hz_seg_planar(hz_pseglist *s, const hz_objmesh *m, double delta) {
  memset(s, 0, sizeof *s);
  s->delta = delta;
  const int32_t nt = m->nt;
  if (nt <= 0) return 0;
  if (!(delta > 0.0)) return 3;

  seg_work w;
  memset(&w, 0, sizeof w);
  int32_t *label = malloc((size_t)nt * sizeof *label);
  w.bb = malloc((size_t)nt * 6 * sizeof *w.bb);
  w.gn = malloc((size_t)nt * 3 * sizeof *w.gn);
  w.ar = malloc((size_t)nt * sizeof *w.ar);
  w.ord = malloc((size_t)nt * sizeof *w.ord);
  w.queue = malloc((size_t)nt * sizeof *w.queue);
  w.rank = malloc((size_t)nt * sizeof *w.rank);
  w.xs = malloc(((size_t)nt + 1) * sizeof *w.xs);
  w.g.lev = malloc((size_t)nt);
  if (label == NULL || w.bb == NULL || w.gn == NULL || w.ar == NULL || w.ord == NULL ||
      w.queue == NULL || w.rank == NULL || w.xs == NULL || w.g.lev == NULL) {
    free(label);
    work_free(&w);
    return 2;
  }

  /* --- предвычисленное по треугольникам --- */
  double lo[3] = {1e300, 1e300, 1e300};
  for (int32_t t = 0; t < nt; t++) {
    double p[3][3];
    hz_obj_tri(m, t, p);
    for (int a = 0; a < 3; a++) {
      double l = p[0][a], h = p[0][a];
      for (int i = 1; i < 3; i++) {
        if (p[i][a] < l) l = p[i][a];
        if (p[i][a] > h) h = p[i][a];
      }
      w.bb[(size_t)t * 6 + (size_t)a] = l;
      w.bb[(size_t)t * 6 + 3 + (size_t)a] = h;
      if (l < lo[a]) lo[a] = l;
    }
    double e1[3], e2[3], cr[3];
    for (int a = 0; a < 3; a++) {
      e1[a] = p[1][a] - p[0][a];
      e2[a] = p[2][a] - p[0][a];
    }
    cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
    cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
    cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
    double cn = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
    w.ar[t] = 0.5 * cn;
    for (int a = 0; a < 3; a++)
      w.gn[(size_t)t * 3 + (size_t)a] = cr[a] / cn; /* вырожденных нет: отсеяны разбором */
    label[t] = -1;
    w.ord[t].area = w.ar[t];
    w.ord[t].t = t;
    for (int a = 0; a < 3; a++)
      w.ord[t].c[a] = 0.5 * (w.bb[(size_t)t * 6 + (size_t)a] + w.bb[(size_t)t * 6 + 3 + (size_t)a]);
  }
  qsort(w.ord, (size_t)nt, sizeof *w.ord, cmp_key);
  for (int32_t i = 0; i < nt; i++)
    w.rank[w.ord[i].t] = i;

  /* --- сетка: уровни, счёт записей --- */
  seg_grid *g = &w.g;
  g->h0 = 2.0 * delta;
  for (int a = 0; a < 3; a++)
    g->o[a] = lo[a] - delta;
  int64_t nent = 0;
  int32_t maxlev = 0;
  for (int32_t t = 0; t < nt; t++) {
    double ext = 0.0;
    for (int a = 0; a < 3; a++) {
      double e = w.bb[(size_t)t * 6 + 3 + (size_t)a] - w.bb[(size_t)t * 6 + (size_t)a];
      if (e > ext) ext = e;
    }
    int32_t lv = 0;
    while (lv < SEG_MAXLEV - 1 && ldexp(g->h0, lv) < ext + 2.0 * delta)
      lv++;
    g->lev[t] = (unsigned char)lv;
    if (lv > maxlev) maxlev = lv;
    g->levused[lv] = 1;
    int64_t c0[3], c1[3];
    cell_range(g, w.bb + (size_t)t * 6, 0.0, lv, c0, c1);
    nent += (c1[0] - c0[0] + 1) * (c1[1] - c0[1] + 1) * (c1[2] - c0[2] + 1);
  }
  g->maxlev = maxlev;
  for (int32_t l = 0; l <= maxlev; l++)
    if (g->levused[l]) s->nlev++;
  if (nent > 1500000000LL) {
    free(label);
    work_free(&w);
    return 2;
  }

  int64_t ns = 4;
  while (ns < 2 * nent + 8)
    ns *= 2;
  g->nslot = ns;
  g->mask = (uint64_t)ns - 1;
  g->nent = nent;
  g->key = calloc((size_t)ns, sizeof *g->key);
  g->start = calloc((size_t)ns, sizeof *g->start);
  g->idx = malloc((size_t)nent * sizeof *g->idx);
  w.tmp = malloc((size_t)ns * sizeof *w.tmp);
  if (g->key == NULL || g->start == NULL || g->idx == NULL || w.tmp == NULL) {
    free(label);
    work_free(&w);
    return 2;
  }
  for (int32_t t = 0; t < nt; t++) {
    int32_t lv = g->lev[t];
    int64_t c0[3], c1[3];
    cell_range(g, w.bb + (size_t)t * 6, 0.0, lv, c0, c1);
    for (int64_t iz = c0[2]; iz <= c1[2]; iz++)
      for (int64_t iy = c0[1]; iy <= c1[1]; iy++)
        for (int64_t ix = c0[0]; ix <= c1[0]; ix++)
          g->start[slot_of(g, ckey(lv, ix, iy, iz), 1)]++;
  }
  int64_t acc = 0;
  for (int64_t i = 0; i < ns; i++) {
    int32_t c = g->start[i];
    g->start[i] = (int32_t)acc;
    w.tmp[i] = (int32_t)acc;
    acc += c;
  }
  for (int32_t t = 0; t < nt; t++) {
    int32_t lv = g->lev[t];
    int64_t c0[3], c1[3];
    cell_range(g, w.bb + (size_t)t * 6, 0.0, lv, c0, c1);
    for (int64_t iz = c0[2]; iz <= c1[2]; iz++)
      for (int64_t iy = c0[1]; iy <= c1[1]; iy++)
        for (int64_t ix = c0[0]; ix <= c1[0]; ix++)
          g->idx[w.tmp[slot_of(g, ckey(lv, ix, iy, iz), 0)]++] = t;
  }
  /* Порядок ВНУТРИ ячейки — геометрический: обход соседей не должен зависеть от
   * перестановки треугольников входа (заголовок). Ячейки короткие ⇒ вставками. */
  for (int64_t i = 0; i < ns; i++) {
    int32_t b, e;
    slot_span(g, i, &b, &e);
    for (int32_t x = b + 1; x < e; x++) {
      int32_t v = g->idx[x], r = w.rank[v], y = x - 1;
      while (y >= b && w.rank[g->idx[y]] > r) {
        g->idx[y + 1] = g->idx[y];
        y--;
      }
      g->idx[y + 1] = v;
    }
  }

  /* --- межуровневые пары: КРУПНЫЙ -> список МЕЛКИХ ---
   * Запрос от треугольника покрывает только уровни НЕ МЕЛЬЧЕ своего (иначе
   * крупный элемент перебирал бы тысячи мелких ячеек). Поэтому пары
   * «мелкий–крупный» собираются один раз с ДЕШЁВОЙ стороны и обращаются. */
  memset(w.xs, 0, ((size_t)nt + 1) * sizeof *w.xs);
  int64_t ncross = 0;
  for (int fill = 0; fill < 2; fill++) {
    if (fill == 1) {
      if (ncross > 1500000000LL) {
        free(label);
        work_free(&w);
        return 2;
      }
      w.xl = malloc((size_t)(ncross > 0 ? ncross : 1) * sizeof *w.xl);
      if (w.xl == NULL) {
        free(label);
        work_free(&w);
        return 2;
      }
      int64_t a2 = 0;
      for (int32_t t = 0; t < nt; t++) {
        int32_t c = w.xs[t];
        w.xs[t] = (int32_t)a2;
        w.tmp[t] = (int32_t)a2;
        a2 += c;
      }
      w.xs[nt] = (int32_t)a2;
    }
    for (int32_t b = 0; b < nt; b++)
      for (int32_t l = (int32_t)g->lev[b] + 1; l <= maxlev; l++) {
        if (!g->levused[l]) continue;
        int64_t c0[3], c1[3];
        cell_range(g, w.bb + (size_t)b * 6, delta, l, c0, c1);
        for (int64_t iz = c0[2]; iz <= c1[2]; iz++)
          for (int64_t iy = c0[1]; iy <= c1[1]; iy++)
            for (int64_t ix = c0[0]; ix <= c1[0]; ix++) {
              int64_t sl = slot_of(g, ckey(l, ix, iy, iz), 0);
              if (sl < 0) continue;
              int32_t e0, e1;
              slot_span(g, sl, &e0, &e1);
              for (int32_t k = e0; k < e1; k++) {
                int32_t a = g->idx[k];
                if (g->lev[a] != l) continue; /* запись из другой ячейки (столкновение) */
                if (!bb_near(w.bb + (size_t)a * 6, w.bb + (size_t)b * 6, delta)) continue;
                if (fill == 0) {
                  ncross++;
                  w.xs[a]++;
                } else {
                  w.xl[w.tmp[a]++] = b;
                }
              }
            }
      }
  }
  s->ncross = ncross;
  /* тот же довод, что и для ячеек: список должен быть геометрически упорядочен */
  for (int32_t a = 0; a < nt; a++) {
    int32_t b = w.xs[a], e = w.xs[a + 1];
    for (int32_t x = b + 1; x < e; x++) {
      int32_t v = w.xl[x], r = w.rank[v], y = x - 1;
      while (y >= b && w.rank[w.xl[y]] > r) {
        w.xl[y + 1] = w.xl[y];
        y--;
      }
      w.xl[y + 1] = v;
    }
  }

  /* --- рост участков --------------------------------------------------------- */
  w.seg = malloc((size_t)nt * sizeof *w.seg);
  if (w.seg == NULL) {
    free(label);
    work_free(&w);
    return 2;
  }
  int32_t nseg = 0;
  int64_t nrounds = 0, nreleased = 0, ncand = 0;
  int32_t left = nt;

  while (left > 0 && nrounds < SEG_MAXROUNDS) {
    nrounds++;
    for (int32_t oi = 0; oi < nt; oi++) {
      int32_t sd = w.ord[oi].t;
      if (label[sd] >= 0) continue;
      int32_t r = nseg++;
      double org[3];
      for (int a = 0; a < 3; a++)
        org[a] = w.ord[oi].c[a];
      seg_fit fit;
      double p[3][3];
      hz_obj_tri(m, sd, p);
      fit_reset(&fit, org);
      fit_add(&fit, p, w.gn + (size_t)sd * 3, w.ar[sd]);
      fit_solve(&fit);
      label[sd] = r;
      int32_t qs = 0, qe = 0, nextrefit = 2;
      w.queue[qe++] = sd;

      while (qs < qe) {
        int32_t x = w.queue[qs++];
        const double *bx = w.bb + (size_t)x * 6;
        for (int32_t l = (int32_t)g->lev[x]; l <= maxlev; l++) {
          if (!g->levused[l]) continue;
          int64_t c0[3], c1[3];
          cell_range(g, bx, delta, l, c0, c1);
          for (int64_t iz = c0[2]; iz <= c1[2]; iz++)
            for (int64_t iy = c0[1]; iy <= c1[1]; iy++)
              for (int64_t ix = c0[0]; ix <= c1[0]; ix++) {
                int64_t sl = slot_of(g, ckey(l, ix, iy, iz), 0);
                if (sl < 0) continue;
                int32_t e0, e1;
                slot_span(g, sl, &e0, &e1);
                for (int32_t k = e0; k < e1; k++) {
                  int32_t c = g->idx[k];
                  ncand++;
                  if (label[c] >= 0 || g->lev[c] != l) continue;
                  if (!bb_near(bx, w.bb + (size_t)c * 6, delta)) continue;
                  if (w.gn[(size_t)c * 3 + 0] * fit.n[0] + w.gn[(size_t)c * 3 + 1] * fit.n[1] +
                          w.gn[(size_t)c * 3 + 2] * fit.n[2] <=
                      0.0)
                    continue;
                  double q[3][3];
                  hz_obj_tri(m, c, q);
                  if (tri_dev(q, fit.n, fit.off) > delta) continue;
                  label[c] = r;
                  fit_add(&fit, q, w.gn + (size_t)c * 3, w.ar[c]);
                  w.queue[qe++] = c;
                  if (qe >= nextrefit) {
                    fit_solve(&fit);
                    nextrefit *= 2;
                  }
                }
              }
        }
        for (int32_t k = w.xs[x]; k < w.xs[x + 1]; k++) {
          int32_t c = w.xl[k];
          ncand++;
          if (label[c] >= 0) continue;
          if (!bb_near(bx, w.bb + (size_t)c * 6, delta)) continue;
          if (w.gn[(size_t)c * 3 + 0] * fit.n[0] + w.gn[(size_t)c * 3 + 1] * fit.n[1] +
                  w.gn[(size_t)c * 3 + 2] * fit.n[2] <=
              0.0)
            continue;
          double q[3][3];
          hz_obj_tri(m, c, q);
          if (tri_dev(q, fit.n, fit.off) > delta) continue;
          label[c] = r;
          fit_add(&fit, q, w.gn + (size_t)c * 3, w.ar[c]);
          w.queue[qe++] = c;
          if (qe >= nextrefit) {
            fit_solve(&fit);
            nextrefit *= 2;
          }
        }
      }
      fit_solve(&fit);

      /* --- ПРОВЕРКА: подгонка по среднему, критерий на МАКСИМУМЕ (Г40/Г44) --- */
      int32_t nmem = qe;
      double dmax = 0.0;
      for (int fixpass = 0;; fixpass++) {
        dmax = 0.0;
        for (int32_t i = 0; i < nmem; i++) {
          double q[3][3];
          hz_obj_tri(m, w.queue[i], q);
          double d = tri_dev(q, fit.n, fit.off);
          if (d > dmax) dmax = d;
        }
        if (!(dmax > delta)) break;

        /* Последний проход: плоскость ЗАТРАВКИ и никакого пересчёта после — иначе
         * подгонка снова уводит плоскость и цикл не сходится. */
        int last = (fixpass >= SEG_MAXFIX);
        double pn[3], poff;
        if (last) {
          for (int a = 0; a < 3; a++)
            pn[a] = w.gn[(size_t)sd * 3 + (size_t)a];
          hz_obj_tri(m, sd, p);
          poff = pn[0] * p[0][0] + pn[1] * p[0][1] + pn[2] * p[0][2];
        } else {
          for (int a = 0; a < 3; a++)
            pn[a] = fit.n[a];
          poff = fit.off;
        }
        int32_t keep = 0;
        seg_fit nf;
        fit_reset(&nf, org);
        for (int32_t i = 0; i < nmem; i++) {
          int32_t t = w.queue[i];
          double q[3][3];
          hz_obj_tri(m, t, q);
          if (t != sd && tri_dev(q, pn, poff) > delta) {
            label[t] = -1;
            nreleased++;
            continue;
          }
          w.queue[keep++] = t;
          fit_add(&nf, q, w.gn + (size_t)t * 3, w.ar[t]);
        }
        /* Нарушала только затравка ⇒ отпускать некого, и подгонка ничего не
         * изменит. Уходим на плоскость затравки немедленно. */
        int stuck = (keep == nmem);
        nmem = keep;
        if (last || stuck) {
          for (int a = 0; a < 3; a++)
            fit.n[a] = pn[a];
          fit.off = poff;
          if (stuck && !last) continue; /* пересчитать dmax на плоскости затравки */
          dmax = 0.0;
          for (int32_t i = 0; i < nmem; i++) {
            double q[3][3];
            hz_obj_tri(m, w.queue[i], q);
            double d = tri_dev(q, fit.n, fit.off);
            if (d > dmax) dmax = d;
          }
          break;
        }
        fit = nf;
        fit_solve(&fit);
      }

      double area = 0.0;
      for (int32_t i = 0; i < nmem; i++)
        area += w.ar[w.queue[i]];
      for (int a = 0; a < 3; a++)
        w.seg[r].n[a] = fit.n[a];
      w.seg[r].off = fit.off;
      w.seg[r].area = area;
      w.seg[r].ntri = nmem;
      w.seg[r].dmax = dmax;
    }
    left = 0;
    for (int32_t t = 0; t < nt; t++)
      if (label[t] < 0) left++;
  }

  hz_pseg *seg = w.seg;
  w.seg = NULL;
  work_free(&w);
  hz_pseg *shrunk = realloc(seg, (size_t)(nseg > 0 ? nseg : 1) * sizeof *seg);
  s->seg = (shrunk != NULL) ? shrunk : seg;
  s->nseg = nseg;
  s->label = label;
  s->nrounds = nrounds;
  s->nreleased = nreleased;
  s->ncand = ncand;
  return 0;
}
