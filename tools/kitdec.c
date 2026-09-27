/* kitdec.c — §914-Ш2'': декиматор v1, VERTEX CLUSTERING.
 *
 * Читает кит (L0), строит ГРУБЫЙ уровень L1: внутренние вершины
 * снапятся в узлы решётки (шаг h = габарит/div, именованная константа),
 * совпавшие склеиваются хешем квантованных ключей, вырожденные
 * треугольники выбрасываются. ГРАНИЧНЫЕ вершины гроздей (вершина
 * принадлежит треугольникам ≥2 гроздей) НЕ двигаются — стык уровней
 * без щелей по построению. Новый уровень дописывается САМЫМ ГРУБЫМ
 * (lev[nlev]) и сохраняется.
 *
 * Синтаксис: kitdec IN.kit OUT.kit [div=N]   (div по умолчанию 8)
 */
#include "geom/kit.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* деление габарита для шага решётки: не магия — 8 даёт умеренную
 * редукцию при запертых границах (§914-Ш2''-П1) */
#define HZ_DEC_DIV_DEFAULT 8
#define HZ_DEC_AREA_TOL 0.1
/* ГРОЗДИ УКРУПНЯЮТСЯ на уровень: 4 грозди источника → 1 нового
 * (квадродерево гроздей, STRUCTURE §2.2 v2 «смена уровня — гроздью»).
 * Без этого доля запертых граничных вершин не падает и лестница не
 * грубеет (найдено на room: 66% вершин граничные, редукция 0.3%). */
#define HZ_DEC_CLUSTER_MERGE 4

/* открытая хеш-таблица квантованных ключей: линейные зонды, степень 2 */
typedef struct {
  uint64_t *keys; /* упакованный квантованный ключ (3×21 бит) */
  int32_t *val;
  int32_t mask, used;
} ihash;

static int ih_init(ihash *h, int32_t cap2) {
  h->mask = cap2 - 1;
  h->used = 0;
  h->keys = malloc((size_t)cap2 * sizeof(int64_t));
  h->val = malloc((size_t)cap2 * sizeof(int32_t));
  if (!h->keys || !h->val) {
    free(h->keys);
    free(h->val);
    h->keys = NULL;
    h->val = NULL;
    return 1;
  }
  for (int32_t i = 0; i < cap2; i++)
    h->val[i] = -1;
  return 0;
}

static void ih_free(ihash *h) {
  free(h->keys);
  free(h->val);
}

static uint64_t mix64(uint64_t x) {
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdull;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ull;
  x ^= x >> 33;
  return x;
}

/* ключ: (ix,iy,iz) по 21 бит, беззнаковая упаковка */
static uint64_t pack_key3(int64_t ix, int64_t iy, int64_t iz) {
  uint64_t ux = (uint64_t)ix, uy = (uint64_t)iy, uz = (uint64_t)iz;
  return ((ux & 0x1fffffull) << 42) | ((uy & 0x1fffffull) << 21) | (uz & 0x1fffffull);
}

static int32_t ih_get(ihash *h, uint64_t key, int32_t *next_id) {
  uint64_t hh = mix64((uint64_t)key);
  int32_t i = (int32_t)(hh & (uint64_t)h->mask);
  for (;;) {
    if (h->val[i] < 0) {
      h->keys[i] = key;
      h->val[i] = (*next_id)++;
      h->used++;
      return h->val[i];
    }
    if (h->keys[i] == key) return h->val[i];
    i = (i + 1) & h->mask;
  }
}

typedef struct {
  double *x, *y, *z;
  uint32_t n;
} vbuf;

static int vb_push(vbuf *v, double px, double py, double pz, uint32_t cap[]) {
  if (v->n >= cap[0]) {
    cap[0] = cap[0] ? cap[0] * 2 : 1024;
    double *nx = realloc(v->x, cap[0] * sizeof(double));
    double *ny = realloc(v->y, cap[0] * sizeof(double));
    double *nz = realloc(v->z, cap[0] * sizeof(double));
    if (!nx || !ny || !nz) {
      free(nx ? nx : v->x);
      free(ny ? ny : v->y);
      free(nz ? nz : v->z);
      v->x = nx;
      v->y = ny;
      v->z = nz;
      return 1;
    }
    v->x = nx;
    v->y = ny;
    v->z = nz;
  }
  v->x[v->n] = px;
  v->y[v->n] = py;
  v->z[v->n] = pz;
  v->n++;
  return 0;
}

int main(int argc, char **argv) {
  int div = HZ_DEC_DIV_DEFAULT;
  uint32_t *vseen = NULL;
  unsigned char *isb = NULL;
  int32_t *map = NULL;
  vbuf vb = {0};
  uint32_t *t1a = NULL, *t1b = NULL, *t1c = NULL, *t1cl = NULL;
  ihash ht = {0};
  if (argc < 3) {
    fprintf(stderr, "kitdec IN.kit OUT.kit [div=N]\n");
    return 2;
  }
  for (int i = 3; i < argc; i++)
    if (strncmp(argv[i], "div=", 4) == 0) div = atoi(argv[i] + 4);
  if (div < 1) {
    fprintf(stderr, "kitdec: div>=1\n");
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
  /* источник — САМЫЙ ГРУБЫЙ существующий уровень: лестница растёт
   * вниз по грубости, а не параллельными копиями L0 */
  const hz_kit_level *L0 = &k.lev[k.nlev - 1];
  printf("вход: nlev=%d источник lev[%d]: nv=%u nt=%u nc=%u\n", k.nlev, k.nlev - 1, L0->nverts,
         L0->ntris, L0->nclust);

  /* шаг решётки от габарита L0 (границы НЕ снапятся — только внутренние) */
  double ext[3] = {0, 0, 0};
  for (uint32_t v = 0; v < L0->nverts; v++) {
    ext[0] += L0->vx[v];
    ext[1] += L0->vy[v];
    ext[2] += L0->vz[v];
  }
  double cen[3] = {ext[0] / L0->nverts, ext[1] / L0->nverts, ext[2] / L0->nverts};
  double rad = 0;
  for (uint32_t v = 0; v < L0->nverts; v++) {
    double d[3] = {L0->vx[v] - cen[0], L0->vy[v] - cen[1], L0->vz[v] - cen[2]};
    double r = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (r > rad) rad = r;
  }
  /* решётка ГРУБЕЕТ на уровень: h ×2 за шаг лестницы (геометрическая
   * прогрессия из бюджета §914); div_effective = div >> (номер уровня) */
  int diveff = div >> (k.nlev - 1);
  if (diveff < 1) diveff = 1;
  double h = 2.0 * rad / (double)diveff;
  if (h <= 0) {
    fprintf(stderr, "kitdec: нулевой габарит\n");
    hz_kit_free(&k);
    return 2;
  }

  /* граничные вершины: принадлежат треугольникам ≥2 гроздей */
  vseen = calloc(L0->nverts, sizeof(uint32_t)); /* битовая маска кластеров ≤ 31 */
  isb = calloc(L0->nverts, 1);
  if (!vseen || !isb) {
    fprintf(stderr, "kitdec: нет памяти\n");
    free(vseen);
    free(isb);
    hz_kit_free(&k);
    return 2;
  }
  for (uint32_t t = 0; t < L0->ntris; t++) {
    uint32_t c = L0->tcl[t] / HZ_DEC_CLUSTER_MERGE + 1; /* 0 = «не встречался» */
    const uint32_t vi[3] = {L0->ti0[t], L0->ti1[t], L0->ti2[t]};
    for (int a = 0; a < 3; a++) {
      uint32_t last = vseen[vi[a]];
      if (last != 0 && last != c) isb[vi[a]] = 1; /* второй кластер */
      vseen[vi[a]] = c;
    }
  }
  uint32_t nbound = 0;
  for (uint32_t v = 0; v < L0->nverts; v++)
    if (isb[v]) nbound++;
  printf("решётка h=%.6f (div=%d, эффективный %d), граничных вершин: %u из %u\n", h, div, diveff,
         nbound, L0->nverts);

  /* отображение вершин: граница — своя; внутренняя — квантованный ключ */
  if (ih_init(&ht, 1 << 16)) {
    fprintf(stderr, "kitdec: нет памяти (hash)\n");
    free(vseen);
    free(isb);
    hz_kit_free(&k);
    return 2;
  }
  map = calloc((size_t)L0->nverts, sizeof(int32_t)); /* calloc: анализатор
     теряет индукцию заполнения через ih_get (FP-класс), нули её чинят */
  uint32_t cap[1] = {0};
  int32_t next_id = 0;
  if (map == NULL) goto oom;
  for (uint32_t v = 0; v < L0->nverts; v++) {
    if (isb[v]) {
      map[v] = next_id++;
      if (vb_push(&vb, L0->vx[v], L0->vy[v], L0->vz[v], cap)) goto oom;
    } else {
      int64_t kx = (int64_t)floor(L0->vx[v] / h);
      int64_t ky = (int64_t)floor(L0->vy[v] / h);
      int64_t kz = (int64_t)floor(L0->vz[v] / h);
      int32_t id = ih_get(&ht, pack_key3(kx, ky, kz), &next_id);
      if ((uint32_t)id == vb.n) {
        /* новая ячейка: представитель — ПЕРВАЯ РЕАЛЬНАЯ вершина ячейки,
         * не узел решётки (узел растягивал слайверы: площадь 1.25 при
         * П3-пороге 1.1 — фальсификатор сработал, правка по его finding) */
        if (vb_push(&vb, L0->vx[v], L0->vy[v], L0->vz[v], cap)) goto oom;
      }
      map[v] = id;
    }
  }

  /* треугольники L1: вырожденные выбрасываем; кластер — от исходного */
  uint32_t ntr1 = 0;
  t1a = malloc(4ull * L0->ntris);
  t1b = malloc(4ull * L0->ntris);
  t1c = malloc(4ull * L0->ntris);
  t1cl = malloc(4ull * L0->ntris);
  if (!t1a || !t1b || !t1c || !t1cl) goto oom;
  for (uint32_t t = 0; t < L0->ntris; t++) {
    uint32_t a = (uint32_t)map[L0->ti0[t]], b = (uint32_t)map[L0->ti1[t]],
             c = (uint32_t)map[L0->ti2[t]];
    if (a == b || b == c || a == c) continue; /* вырожденный — в ноль */
    t1a[ntr1] = a;
    t1b[ntr1] = b;
    t1c[ntr1] = c;
    t1cl[ntr1] = L0->tcl[t] / HZ_DEC_CLUSTER_MERGE;
    ntr1++;
  }
  printf("L1: nv=%u nt=%u (сокращение %.1f%%)\n", vb.n, ntr1,
         100.0 * (1.0 - (double)ntr1 / (double)L0->ntris));
  if (ntr1 == 0 || ntr1 >= L0->ntris) {
    fprintf(stderr, "kitdec: лестница нарушена (nt1=%u ≥ nt0=%u) — отказ\n", ntr1, L0->ntris);
    ih_free(&ht);
    free(vseen);
    free(isb);
    free(map);
    free(vb.x);
    free(vb.y);
    free(vb.z);
    free(t1a);
    free(t1b);
    free(t1c);
    free(t1cl);
    hz_kit_free(&k);
    return 2;
  }

  /* суммы площадей — прибор П3 */
  double area0 = 0, area1 = 0;
  for (uint32_t t = 0; t < L0->ntris; t++) {
    uint32_t vi[3] = {L0->ti0[t], L0->ti1[t], L0->ti2[t]};
    double q[3][3];
    for (int a = 0; a < 3; a++) {
      q[a][0] = L0->vx[vi[a]];
      q[a][1] = L0->vy[vi[a]];
      q[a][2] = L0->vz[vi[a]];
    }
    double e1[3] = {q[1][0] - q[0][0], q[1][1] - q[0][1], q[1][2] - q[0][2]};
    double e2[3] = {q[2][0] - q[0][0], q[2][1] - q[0][1], q[2][2] - q[0][2]};
    double cr[3] = {e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                    e1[0] * e2[1] - e1[1] * e2[0]};
    area0 += 0.5 * sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
  }
  for (uint32_t t = 0; t < ntr1; t++) {
    double q[3][3];
    uint32_t vi[3] = {t1a[t], t1b[t], t1c[t]};
    for (int a = 0; a < 3; a++) {
      q[a][0] = vb.x[vi[a]];
      q[a][1] = vb.y[vi[a]];
      q[a][2] = vb.z[vi[a]];
    }
    double e1[3] = {q[1][0] - q[0][0], q[1][1] - q[0][1], q[1][2] - q[0][2]};
    double e2[3] = {q[2][0] - q[0][0], q[2][1] - q[0][1], q[2][2] - q[0][2]};
    double cr[3] = {e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
                    e1[0] * e2[1] - e1[1] * e2[0]};
    area1 += 0.5 * sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
  }
  printf("площадь: L0=%.4f L1=%.4f отношение=%.4f\n", area0, area1, area1 / area0);
  if (area0 > 0 && fabs(area1 / area0 - 1.0) > HZ_DEC_AREA_TOL) {
    fprintf(stderr, "kitdec: площадь вне допуска %.3f (tol %d%%) — уровень НЕ записан\n",
            area1 / area0, (int)(HZ_DEC_AREA_TOL * 100));
    ih_free(&ht);
    free(vseen);
    free(isb);
    free(map);
    free(vb.x);
    free(vb.y);
    free(vb.z);
    free(t1a);
    free(t1b);
    free(t1c);
    free(t1cl);
    hz_kit_free(&k);
    return 2;
  }

  /* проверка П2: граничные вершины битово в L1 (они вошли как есть) */
  {
    uint32_t bad = 0, checked = 0;
    int32_t *firstof = malloc((size_t)vb.n * sizeof(int32_t));
    if (firstof == NULL) goto oom;
    for (uint32_t i = 0; i < vb.n; i++)
      firstof[i] = -1;
    for (uint32_t v = 0; v < L0->nverts; v++) {
      if (!isb[v]) continue;
      int32_t id = map[v];
      if (firstof[id] < 0) firstof[id] = (int32_t)v;
      checked++;
      /* координаты отображения обязаны совпасть битово */
      if (memcmp(&vb.x[id], &L0->vx[v], sizeof(double)) != 0 ||
          memcmp(&vb.y[id], &L0->vy[v], sizeof(double)) != 0 ||
          memcmp(&vb.z[id], &L0->vz[v], sizeof(double)) != 0)
        bad++; /* битовое сравнение — НАМЕРЕННО memcmp, не == */
    }
    free(firstof);
    printf("граница: проверено %u, расхождений %u %s\n", checked, bad,
           bad == 0 ? "OK (битово)" : "FAIL");
  }

  /* собрать L1-уровень и дописать самым грубым */
  hz_kit_level *L1 = realloc(k.lev, (size_t)(k.nlev + 1) * sizeof(hz_kit_level));
  if (L1 == NULL) goto oom;
  k.lev = L1;
  memset(&k.lev[k.nlev], 0, sizeof(hz_kit_level));
  /* present-таблица рассчитана на старое nlev: после достройки кит
   * ЦЕЛИКОМ в памяти — семантика present=NULL («все уровни»).
   * ASAN поймал чтение present[nlev_old] в validate (0 bytes after). */
  free(k.present);
  k.present = NULL;
  hz_kit_level *N = &k.lev[k.nlev];
  k.nlev++;
  N->nverts = vb.n;
  N->ntris = ntr1;
  N->vx = vb.x;
  N->vy = vb.y;
  N->vz = vb.z;
  N->ti0 = t1a;
  N->ti1 = t1b;
  N->ti2 = t1c;
  N->tcl = malloc(4ull * ntr1);
  N->tmtl = malloc(4ull * ntr1);
  /* грозди L1: непрерывные прогоны одного t1cl (порядок исходный) */
  if (!N->tcl || !N->tmtl) goto oom;
  uint32_t nc1 = 0;
  for (uint32_t t = 0; t < ntr1; t++) {
    if (t == 0 || t1cl[t] != t1cl[t - 1]) nc1++;
    N->tcl[t] = nc1 - 1;
  }
  N->nclust = nc1;
  N->cl = calloc(nc1 ? nc1 : 1, sizeof(hz_cluster));
  if (N->cl == NULL) goto oom;
  uint32_t ci = 0;
  for (uint32_t t = 0; t < ntr1; t++) {
    if (t == 0 || t1cl[t] != t1cl[t - 1]) {
      if (t > 0) ci++;
      N->cl[ci].first_tri = t;
    }
    N->cl[ci].ntris++;
    N->tmtl[t] = 0; /* v1: материал уровня не отслеживается — L1 диагностический */
  }
  for (uint32_t c = 0; c < nc1; c++) {
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
    g->err = 0.0f; /* v1: ε не оценивается — прибор площади сверху */
  }
  free(t1cl); /* временный: владение уровню НЕ передаётся (ASAN нашёл) */
  t1cl = NULL;
  rc = hz_kit_validate(&k);
  if (rc != HZ_KIT_OK) {
    fprintf(stderr, "kitdec: итоговый кит невалиден (rc=%d)\n", rc);
    hz_kit_free(&k);
    ih_free(&ht);
    free(vseen);
    free(isb);
    free(map);
    return 2;
  }
  f = fopen(argv[2], "wb");
  if (f == NULL) {
    fprintf(stderr, "kitdec: не открывается %s\n", argv[2]);
    hz_kit_free(&k);
    ih_free(&ht);
    free(vseen);
    free(isb);
    free(map);
    return 2;
  }
  rc = hz_kit_save(&k, f);
  if (fclose(f) != 0) rc = HZ_KIT_E_IO;
  printf("записан %s: nlev=%d (L1: nv=%u nt=%u nc=%u) rc=%d\n", argv[2], k.nlev, N->nverts,
         N->ntris, N->nclust, rc);
  ih_free(&ht);
  free(vseen);
  free(isb);
  free(map);
  hz_kit_free(&k);
  return rc == HZ_KIT_OK ? 0 : 2;
oom:
  fprintf(stderr, "kitdec: нет памяти\n");
  ih_free(&ht);
  free(vseen);
  free(isb);
  free(map);
  free(vb.x);
  free(vb.y);
  free(vb.z);
  free(t1a);
  free(t1b);
  free(t1c);
  free(t1cl);
  hz_kit_free(&k);
  return 2;
}
