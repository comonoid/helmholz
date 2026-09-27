/* kitdec.c — §916-бис: декиматор v3 — QEM + CSR-СМЕЖНОСТЬ (город).
 *
 * v2 делал O(nt)-сканов НА КОЛЛАПС (link/flip/re-add) — на San_Miguel
 * 4.96М граней это квадратично. v3: смежность вершина→треугольники
 * linked-списками (head/next + доп-слоты в хвосте при переносе) —
 * все операции коллапса O(валентность).
 *
 * Политика (указание пользователя: «слишком много этажей качества»):
 * укрупнение гроздей ×16 за уровень, target=0.25 — жирные прыжки,
 * 3–4 уровня до ×40+, а не 9 одинаковых.
 *
 * QEM: Q[v] = Σ area·pp^T; стоимость ребра = (Q_a+Q_b)(середина) —
 * АРГМИН УБРАН из горячего пути (позиция коллапса = НЕПОДВИЖНАЯ
 * вершина-приёмник: плоские кровли с вырожденными квадриками блуждали
 * цепочками оптимумов, треугольники до 529 м², §916). Запреты:
 * граничные вершины (≥2 укрупнённых гроздей) не двигаются и не
 * сливаются друг с другом; LINK CONDITION (ровно 2 общих соседа);
 * переворот нормали; рост площади треугольника ×4. Приёмка: площадь
 * уровня ±10% — жёсткий гейт записи.
 *
 * Синтаксис: kitdec IN.kit OUT.kit [target=F]   (F=0.25)
 */
#include "geom/kit.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HZ_DEC_TARGET_DEFAULT 0.25
#define HZ_DEC_AREA_TOL 0.1
#define HZ_DEC_CLUSTER_MERGE 16 /* §916-бис: жирные прыжки лестницы */
#define HZ_DEC_AREA_GROW 4.0    /* запрет роста площади треугольника */

typedef struct {
  double a[10];
} quad;

static void quad_add(quad *r, const quad *q) {
  for (int i = 0; i < 10; i++)
    r->a[i] += q->a[i];
}

static void quad_plane(quad *q, double nx, double ny, double nz, double d, double w) {
  double v[4] = {nx, ny, nz, d};
  int k = 0;
  for (int i = 0; i < 4; i++)
    for (int j = i; j < 4; j++)
      q->a[k++] = w * v[i] * v[j];
}

static double quad_eval(const quad *q, const double x[3]) {
  double v[4] = {x[0], x[1], x[2], 1.0};
  static const int idx[4][4] = {{0, 1, 2, 3}, {1, 4, 5, 6}, {2, 5, 7, 8}, {3, 6, 8, 9}};
  double s = 0;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      s += v[i] * q->a[idx[i][j]] * v[j];
  return s;
}

typedef struct {
  int32_t a, b;
  double cost;
} edge_item;

typedef struct {
  edge_item *h;
  int32_t n, cap;
} heap;

static void heap_push(heap *hp, int32_t a, int32_t b, double cost) {
  if (hp->n >= hp->cap) {
    int32_t nc = hp->cap ? hp->cap * 2 : 4096;
    edge_item *nh = realloc(hp->h, (size_t)nc * sizeof *nh);
    if (nh == NULL) return;
    hp->h = nh;
    hp->cap = nc;
  }
  int32_t i = hp->n++;
  hp->h[i].a = a;
  hp->h[i].b = b;
  hp->h[i].cost = cost;
  while (i > 0) {
    int32_t p = (i - 1) / 2;
    if (hp->h[p].cost <= hp->h[i].cost) break;
    edge_item t = hp->h[p];
    hp->h[p] = hp->h[i];
    hp->h[i] = t;
    i = p;
  }
}

static int heap_pop(heap *hp, edge_item *out) {
  if (hp->n == 0) return 0;
  *out = hp->h[0];
  hp->h[0] = hp->h[--hp->n];
  int32_t i = 0;
  for (;;) {
    int32_t l = 2 * i + 1, r = l + 1, m = i;
    if (l < hp->n && hp->h[l].cost < hp->h[m].cost) m = l;
    if (r < hp->n && hp->h[r].cost < hp->h[m].cost) m = r;
    if (m == i) break;
    edge_item t = hp->h[m];
    hp->h[m] = hp->h[i];
    hp->h[i] = t;
    i = m;
  }
  return 1;
}

static int tri_normal(const double p[3][3], double n[3]) {
  double e1[3] = {p[1][0] - p[0][0], p[1][1] - p[0][1], p[1][2] - p[0][2]};
  double e2[3] = {p[2][0] - p[0][0], p[2][1] - p[0][1], p[2][2] - p[0][2]};
  n[0] = e1[1] * e2[2] - e1[2] * e2[1];
  n[1] = e1[2] * e2[0] - e1[0] * e2[2];
  n[2] = e1[0] * e2[1] - e1[1] * e2[0];
  double nn = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  return nn > 0;
}

/* --- смежность: слоты [0,3nt) кодируют треугольник slot/3;
 * доп-слоты [3nt,...) — перенесённые вхождения, их треугольник в etri */
typedef struct {
  int32_t *head;  /* [nv] */
  int32_t *next;  /* [ёмкость] */
  uint32_t *etri; /* [доп-слоты] треугольник доп-слота */
  int32_t nslot, cap;
  int32_t nextra;
  uint32_t nt;
} adj;

static int adj_init(adj *A, uint32_t nv, uint32_t nt, const uint32_t *ta, const uint32_t *tb,
                    const uint32_t *tc) {
  A->nt = nt;
  A->nslot = (int32_t)(3.0 * (double)nt); /* nt ≤ 1e9×3 < 2^31-1: кит-кэп гарантирует; через double
                                             без sign-conversion */
  A->cap = A->nslot + 1024;
  A->nextra = 0;
  A->head = malloc((size_t)nv * sizeof(int32_t));
  A->next = malloc((size_t)A->cap * sizeof(int32_t));
  A->etri = malloc(1024 * sizeof(uint32_t));
  if (!A->head || !A->next || !A->etri) return 1;
  for (uint32_t v = 0; v < nv; v++)
    A->head[v] = -1;
  for (uint32_t t = 0; t < nt; t++) {
    uint32_t vi[3] = {ta[t], tb[t], tc[t]};
    for (int i = 0; i < 3; i++) {
      int32_t slot = (int32_t)(3) * (int32_t)t + i; /* t < 2^29: 3t < 2^31, кит-кэп */
      A->next[slot] = A->head[vi[i]];
      A->head[vi[i]] = slot;
    }
  }
  return 0;
}

/* добавить треугольник t в список вершины v (доп-слот) */
static int adj_add(adj *A, uint32_t v, uint32_t t) {
  if (A->nslot + A->nextra >= A->cap) {
    int32_t nc = A->cap * 2;
    int32_t *nn = realloc(A->next, (size_t)nc * sizeof(int32_t));
    if (nn == NULL) return 1;
    A->next = nn;
    uint32_t *ne = realloc(A->etri, (size_t)(nc - A->nslot) * sizeof(uint32_t));
    if (ne == NULL) return 1;
    A->etri = ne;
    A->cap = nc;
  }
  int32_t slot = A->nslot + A->nextra++;
  A->etri[slot - A->nslot] = t;
  A->next[slot] = A->head[v];
  A->head[v] = slot;
  return 0;
}

static uint32_t adj_tri(const adj *A, int32_t slot) {
  return slot < A->nslot ? (uint32_t)(slot / 3) : A->etri[slot - A->nslot];
}

#define ADJ_FOR(A, V, T)                                                                           \
  for (int32_t slot_ = (A)->head[(V)]; slot_ >= 0; slot_ = (A)->next[slot_])                       \
    if ((T) = adj_tri((A), slot_), !(tdel[T]))

/* общие соседи a и b по живым треугольникам (link condition) */
static uint32_t link_shared(const adj *A, uint32_t a, uint32_t b, const uint32_t *ta,
                            const uint32_t *tb, const uint32_t *tc, const unsigned char *tdel,
                            unsigned char *mark, uint32_t *stk) {
  uint32_t ns = 0;
  uint32_t t;
  ADJ_FOR(A, a, t) {
    uint32_t vi[3] = {ta[t], tb[t], tc[t]};
    for (int i = 0; i < 3; i++)
      if (vi[i] != a && !mark[vi[i]]) {
        mark[vi[i]] = 1;
        stk[ns++] = vi[i];
      }
  }
  uint32_t shared = 0;
  ADJ_FOR(A, b, t) {
    uint32_t vi[3] = {ta[t], tb[t], tc[t]};
    for (int i = 0; i < 3; i++)
      if (vi[i] != b && mark[vi[i]]) {
        mark[vi[i]] = 0; /* УНИКАЛЬНОСТЬ: сосед в двух треугольниках
          считался дважды — ссылка «234к блокировок» была этим багом */
        shared++;
      }
  }
  for (uint32_t i = 0; i < ns; i++)
    mark[stk[i]] = 0;
  return shared;
}

int main(int argc, char **argv) {
  double target = HZ_DEC_TARGET_DEFAULT;
  if (argc < 3) {
    fprintf(stderr, "kitdec IN.kit OUT.kit [target=F]\n");
    return 2;
  }
  for (int i = 3; i < argc; i++)
    if (strncmp(argv[i], "target=", 7) == 0) target = atof(argv[i] + 7);
  if (target <= 0.0 || target >= 1.0) {
    fprintf(stderr, "kitdec: 0 < target < 1\n");
    return 2;
  }

  hz_kit k;
  hz_kit_init(&k);
  FILE *f = fopen(argv[1], "rb");
  if (f == NULL) {
    fprintf(stderr, "kitdec: не читается %s\n", argv[1]);
    return 2;
  }
  int rc = hz_kit_load(&k, f);
  fclose(f);
  if (rc != HZ_KIT_OK) {
    fprintf(stderr, "kitdec: кит не читается (rc=%d)\n", rc);
    return 2;
  }
  const hz_kit_level *S = &k.lev[k.nlev - 1];
  printf("вход: nlev=%d источник lev[%d]: nv=%u nt=%u nc=%u\n", k.nlev, k.nlev - 1, S->nverts,
         S->ntris, S->nclust);
  uint32_t nt0 = S->ntris;

  uint32_t nv = S->nverts, nt = S->ntris;
  double *vx = malloc((size_t)nv * sizeof(double));
  double *vy = malloc((size_t)nv * sizeof(double));
  double *vz = malloc((size_t)nv * sizeof(double));
  uint32_t *ta = malloc(4ull * nt), *tb = malloc(4ull * nt), *tc = malloc(4ull * nt);
  uint32_t *tcl = malloc(4ull * nt), *tmtl = malloc(4ull * nt);
  quad *Q = calloc((size_t)nv, sizeof(quad));
  unsigned char *dead = calloc((size_t)nv, 1);
  unsigned char *isb = calloc((size_t)nv, 1);
  unsigned char *tdel = calloc((size_t)nt, 1);
  unsigned char *mark = calloc((size_t)nv, 1);
  uint32_t *stk = malloc((size_t)nv * sizeof(uint32_t));
  heap hp = {0};
  adj A = {0};
  if (!vx || !vy || !vz || !ta || !tb || !tc || !tcl || !tmtl || !Q || !dead || !isb || !tdel ||
      !mark || !stk)
    goto fail;
  memcpy(vx, S->vx, (size_t)nv * sizeof(double));
  memcpy(vy, S->vy, (size_t)nv * sizeof(double));
  memcpy(vz, S->vz, (size_t)nv * sizeof(double));
  memcpy(ta, S->ti0, 4ull * nt);
  memcpy(tb, S->ti1, 4ull * nt);
  memcpy(tc, S->ti2, 4ull * nt);
  memcpy(tcl, S->tcl, 4ull * nt);
  memcpy(tmtl, S->tmtl, 4ull * nt);

  /* площадь источника ДО realloc(k.lev) — S висит внутри lev */
  double area0 = 0;
  for (uint32_t t = 0; t < nt; t++) {
    double p[3][3];
    uint32_t vi[3] = {ta[t], tb[t], tc[t]};
    for (int a = 0; a < 3; a++) {
      p[a][0] = vx[vi[a]];
      p[a][1] = vy[vi[a]];
      p[a][2] = vz[vi[a]];
    }
    double n[3];
    if (tri_normal(p, n)) area0 += 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  }

  if (adj_init(&A, nv, nt, ta, tb, tc)) goto fail;

  { /* границы: вершина в треугольниках ≥2 укрупнённых гроздей */
    uint32_t *last = calloc((size_t)nv, sizeof(uint32_t));
    if (last == NULL) goto fail;
    for (uint32_t t = 0; t < nt; t++) {
      uint32_t c = tcl[t] / HZ_DEC_CLUSTER_MERGE + 1;
      uint32_t vi[3] = {ta[t], tb[t], tc[t]};
      for (int a = 0; a < 3; a++) {
        uint32_t pv = last[vi[a]];
        if (pv != 0 && pv != c) isb[vi[a]] = 1;
        last[vi[a]] = c;
      }
    }
    free(last);
  }

  for (uint32_t t = 0; t < nt; t++) {
    double p[3][3];
    uint32_t vi[3] = {ta[t], tb[t], tc[t]};
    for (int a = 0; a < 3; a++) {
      p[a][0] = vx[vi[a]];
      p[a][1] = vy[vi[a]];
      p[a][2] = vz[vi[a]];
    }
    double n[3];
    if (!tri_normal(p, n)) continue;
    double nn = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    double ar = 0.5 * nn;
    double nx = n[0] / nn, ny = n[1] / nn, nz = n[2] / nn;
    double d = -(nx * p[0][0] + ny * p[0][1] + nz * p[0][2]);
    quad q;
    memset(&q, 0, sizeof q);
    quad_plane(&q, nx, ny, nz, d, ar);
    quad_add(&Q[vi[0]], &q);
    quad_add(&Q[vi[1]], &q);
    quad_add(&Q[vi[2]], &q);
  }

  { /* начальная куча: стоимость = квадрика суммы в СЕРЕДИНЕ (аргмин
     * убран из горячего пути; позиция коллапса всё равно = приёмник) */
    for (uint32_t t = 0; t < nt; t++) {
      uint32_t e[3][2] = {{ta[t], tb[t]}, {tb[t], tc[t]}, {tc[t], ta[t]}};
      for (int i = 0; i < 3; i++) {
        if (e[i][0] >= e[i][1]) continue;
        quad qs = Q[e[i][0]];
        quad_add(&qs, &Q[e[i][1]]);
        double mid[3] = {0.5 * (vx[e[i][0]] + vx[e[i][1]]), 0.5 * (vy[e[i][0]] + vy[e[i][1]]),
                         0.5 * (vz[e[i][0]] + vz[e[i][1]])};
        double cost = quad_eval(&qs, mid);
        if (!isfinite(cost)) continue;
        heap_push(&hp, (int32_t)e[i][0], (int32_t)e[i][1], cost);
      }
    }
  }

  uint32_t want = (uint32_t)(target * (double)nt0);
  uint32_t alive = nt;
  int64_t ncol = 0, nblock = 0, nbb = 0, nlink = 0, nflip = 0;
  while (alive > want && hp.n > 0) {
    edge_item it;
    if (!heap_pop(&hp, &it)) break;
    uint32_t a = (uint32_t)it.a, b = (uint32_t)it.b;
    if (dead[a] || dead[b]) continue;
    if (isb[a] && isb[b]) {
      nblock++;
      nbb++;
      continue;
    }
    uint32_t keepv = isb[b] ? b : a;
    uint32_t killv = isb[b] ? a : b;
    /* LINK CONDITION: внутренних рёбер — ровно 2 общих соседа;
     * КРАЕВЫХ (открытая поверхность: кровля/колонны без крышек) —
     * ровно 1: коллапс краевого ребра легален, складки ловят
     * flip/area-запреты ниже */
    {
      uint32_t sh = link_shared(&A, a, b, ta, tb, tc, tdel, mark, stk);
      if (sh != 2 && sh != 1) {
        nblock++;
        nlink++;
        continue;
      }
    }
    /* переворот нормали и рост площади — по смежности killv */
    int bad = 0;
    {
      uint32_t t;
      ADJ_FOR(&A, killv, t) {
        uint32_t ia = ta[t], ib = tb[t], ic = tc[t];
        int hasp = (ia == keepv || ib == keepv || ic == keepv);
        if (hasp) continue; /* выродится — не проверяем */
        double p[3][3] = {
            {vx[ia], vy[ia], vz[ia]}, {vx[ib], vy[ib], vz[ib]}, {vx[ic], vy[ic], vz[ic]}};
        double o[3][3];
        memcpy(o, p, sizeof o);
        uint32_t ki[3] = {ia, ib, ic};
        for (int r = 0; r < 3; r++)
          if (ki[r] == killv) {
            o[r][0] = vx[keepv];
            o[r][1] = vy[keepv];
            o[r][2] = vz[keepv];
          }
        double n1[3], n2[3];
        if (!tri_normal(p, n1) || !tri_normal(o, n2)) continue;
        if (n1[0] * n2[0] + n1[1] * n2[1] + n1[2] * n2[2] <= 0) {
          bad = 1;
          break;
        }
        double s1 = sqrt(n1[0] * n1[0] + n1[1] * n1[1] + n1[2] * n1[2]);
        double s2 = sqrt(n2[0] * n2[0] + n2[1] * n2[1] + n2[2] * n2[2]);
        if (s2 > HZ_DEC_AREA_GROW * s1 + 1e-12) {
          bad = 1;
          break;
        }
      }
    }
    if (bad) {
      nblock++;
      nflip++;
      continue;
    }
    /* исполнить: killv → keepv (приёмник НЕПОДВИЖЕН) */
    dead[killv] = 1;
    quad_add(&Q[keepv], &Q[killv]);
    {
      uint32_t t;
      ADJ_FOR(&A, killv, t) {
        uint32_t *vi[3] = {&ta[t], &tb[t], &tc[t]};
        for (int i2 = 0; i2 < 3; i2++)
          if (*vi[i2] == killv) *vi[i2] = keepv;
        if (ta[t] == tb[t] || tb[t] == tc[t] || ta[t] == tc[t]) {
          tdel[t] = 1;
          alive--;
        } else if (adj_add(&A, keepv, t)) {
          goto fail;
        }
      }
    }
    ncol++;
    /* пере-добавление рёбер вокруг keepv — O(валентность) */
    {
      uint32_t t;
      ADJ_FOR(&A, keepv, t) {
        uint32_t vi[3] = {ta[t], tb[t], tc[t]};
        for (int i2 = 0; i2 < 3; i2++) {
          if (vi[i2] == keepv) continue;
          uint32_t o = vi[i2];
          if (dead[o] || (isb[o] && isb[keepv])) continue;
          uint32_t lo = keepv < o ? keepv : o, hi = keepv < o ? o : keepv;
          quad qs = Q[lo];
          quad_add(&qs, &Q[hi]);
          double mid[3] = {0.5 * (vx[lo] + vx[hi]), 0.5 * (vy[lo] + vy[hi]),
                           0.5 * (vz[lo] + vz[hi])};
          double c2 = quad_eval(&qs, mid);
          if (isfinite(c2)) heap_push(&hp, (int32_t)lo, (int32_t)hi, c2);
        }
      }
    }
  }
  printf("коллапсов=%lld запрещено=%lld (границы-оба=%lld link=%lld flip/area=%lld): живых %u из "
         "%u (цель ≤ %u)\n",
         (long long)ncol, (long long)nblock, (long long)nbb, (long long)nlink, (long long)nflip,
         alive, nt0, want);

  /* --- сборка уровня (как в v2) --- */
  uint32_t *remap = malloc((size_t)nv * sizeof(uint32_t));
  if (remap == NULL) goto fail;
  uint32_t nnv = 0;
  for (uint32_t v = 0; v < nv; v++)
    remap[v] = dead[v] ? 0xFFFFFFFFu : nnv++;
  uint32_t nnt = alive;
  if (nnt == 0 || nnt >= nt0) {
    fprintf(stderr, "kitdec: лестница нарушена (nt=%u против %u) — отказ\n", nnt, nt0);
    free(remap);
    goto fail;
  }
  hz_kit_level *L1 = realloc(k.lev, (size_t)(k.nlev + 1) * sizeof(hz_kit_level));
  if (L1 == NULL) {
    free(remap);
    goto fail;
  }
  k.lev = L1;
  memset(&k.lev[k.nlev], 0, sizeof(hz_kit_level));
  hz_kit_level *N = &k.lev[k.nlev];
  k.nlev++;
  free(k.present);
  k.present = NULL;
  N->nverts = nnv;
  N->ntris = nnt;
  N->vx = malloc((size_t)nnv * sizeof(double));
  N->vy = malloc((size_t)nnv * sizeof(double));
  N->vz = malloc((size_t)nnv * sizeof(double));
  N->ti0 = malloc(4ull * nnt);
  N->ti1 = malloc(4ull * nnt);
  N->ti2 = malloc(4ull * nnt);
  N->tcl = malloc(4ull * nnt);
  N->tmtl = malloc(4ull * nnt);
  if (!N->vx || !N->vy || !N->vz || !N->ti0 || !N->ti1 || !N->ti2 || !N->tcl || !N->tmtl) {
    free(remap);
    goto fail;
  }
  for (uint32_t v = 0; v < nv; v++)
    if (!dead[v]) {
      N->vx[remap[v]] = vx[v];
      N->vy[remap[v]] = vy[v];
      N->vz[remap[v]] = vz[v];
    }
  uint32_t w = 0, ncl = 0;
  int prev = -1;
  for (uint32_t t = 0; t < nt; t++) {
    if (tdel[t]) continue;
    uint32_t c = tcl[t] / HZ_DEC_CLUSTER_MERGE;
    if ((int)c != prev) {
      ncl++;
      prev = (int)c;
    }
    N->ti0[w] = remap[ta[t]];
    N->ti1[w] = remap[tb[t]];
    N->ti2[w] = remap[tc[t]];
    N->tmtl[w] = tmtl[t];
    N->tcl[w] = ncl - 1;
    w++;
  }
  N->nclust = ncl;
  N->cl = calloc(ncl ? ncl : 1, sizeof(hz_cluster));
  if (N->cl == NULL) {
    free(remap);
    goto fail;
  }
  {
    uint32_t ci = 0;
    for (uint32_t t = 0; t < nnt; t++) {
      if (t == 0 || N->tcl[t] != N->tcl[t - 1]) {
        if (t > 0) ci++;
        N->cl[ci].first_tri = t;
      }
      N->cl[ci].ntris++;
    }
  }
  for (uint32_t c = 0; c < ncl; c++) {
    hz_cluster *g = &N->cl[c];
    double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    int first = 1;
    for (uint32_t t = g->first_tri; t < g->first_tri + g->ntris; t++) {
      uint32_t vi[3] = {N->ti0[t], N->ti1[t], N->ti2[t]};
      for (int a = 0; a < 3; a++) {
        double p[3] = {N->vx[vi[a]], N->vy[vi[a]], N->vz[vi[a]]};
        for (int ax = 0; ax < 3; ax++) {
          if (first || p[ax] < lo[ax]) lo[ax] = p[ax];
          if (first || p[ax] > hi[ax]) hi[ax] = p[ax];
        }
        first = 0;
      }
    }
    for (int ax = 0; ax < 3; ax++) {
      g->bmin[ax] = (float)lo[ax];
      g->bmax[ax] = (float)hi[ax];
    }
    g->err = 0.0f;
  }

  double area1 = 0;
  for (uint32_t t = 0; t < nnt; t++) {
    double p[3][3];
    uint32_t vi[3] = {N->ti0[t], N->ti1[t], N->ti2[t]};
    for (int a = 0; a < 3; a++) {
      p[a][0] = N->vx[vi[a]];
      p[a][1] = N->vy[vi[a]];
      p[a][2] = N->vz[vi[a]];
    }
    double n[3];
    if (tri_normal(p, n)) area1 += 0.5 * sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  }
  printf("площадь: источник=%.4f новый=%.4f отношение=%.4f\n", area0, area1, area1 / area0);
  /* ПОЛИТИКА §916-бис: РОСТ > 1+tol — разрывы/складки, отказ; УБЫВАНИЕ
   * — легитимная грубость (колонны → рёбра), но >×0.33 потери — уже
   * не геометрия. Асимметричный допуск вместо ±10% (v2 блокировал
   * честное огрубление: 0.556 при цели 25% кусков). */
  if (area0 > 0 && (area1 / area0 > 1.0 + HZ_DEC_AREA_TOL || area1 / area0 < 0.33)) {
    fprintf(stderr, "kitdec: площадь вне допуска [0.33, %.2f]: %.3f — НЕ записан\n",
            1.0 + HZ_DEC_AREA_TOL, area1 / area0);
    free(remap);
    goto fail;
  }

  rc = hz_kit_validate(&k);
  if (rc != HZ_KIT_OK) {
    fprintf(stderr, "kitdec: итоговый кит невалиден (rc=%d)\n", rc);
    free(remap);
    goto fail;
  }
  f = fopen(argv[2], "wb");
  if (f == NULL) {
    fprintf(stderr, "kitdec: не открывается %s\n", argv[2]);
    free(remap);
    goto fail;
  }
  rc = hz_kit_save(&k, f);
  if (fclose(f) != 0) rc = HZ_KIT_E_IO;
  printf("записан %s: nlev=%d (новый: nv=%u nt=%u nc=%u) rc=%d\n", argv[2], k.nlev, N->nverts,
         N->ntris, N->nclust, rc);
  free(remap);
  free(vx);
  free(vy);
  free(vz);
  free(ta);
  free(tb);
  free(tc);
  free(tcl);
  free(tmtl);
  free(Q);
  free(dead);
  free(isb);
  free(tdel);
  free(mark);
  free(stk);
  free(hp.h);
  free(A.head);
  free(A.next);
  free(A.etri);
  hz_kit_free(&k);
  return rc == HZ_KIT_OK ? 0 : 2;
fail:
  fprintf(stderr, "kitdec: отказ\n");
  free(vx);
  free(vy);
  free(vz);
  free(ta);
  free(tb);
  free(tc);
  free(tcl);
  free(tmtl);
  free(Q);
  free(dead);
  free(isb);
  free(tdel);
  free(mark);
  free(stk);
  free(hp.h);
  free(A.head);
  free(A.next);
  free(A.etri);
  hz_kit_free(&k);
  return 2;
}
