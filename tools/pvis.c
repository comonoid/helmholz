/* pvis.c — ШАГ О43 (PLAN_ELEMENTS.md, §110): НАСЛЕДУЕТСЯ ЛИ ВИДИМОСТЬ ВНИЗ ПО
 * ЛЕСТНИЦЕ LOD, замеренное НА ГОРОДЕ С КАМЕРОЙ ВНУТРИ УЛИЦЫ.
 *
 * ЗАЧЕМ ОТДЕЛЬНЫЙ СТЕНД, А НЕ ПРАВКА `prad`. Приёмники сборки сидят на ЛИСТЬЯХ
 * (А241), город даёт 446 013 участков, и полная сборка полукубами невыполнима.
 * Но проверяемое утверждение ЛОКАЛЬНО по паре «родитель — ребёнок»: полукуб
 * грубого узла запоминает узлы пространственного дерева, давшие видимый пиксель,
 * а его дети рисуют только их. Значит достаточно ВЫБОРКИ пар, а не сборки.
 *
 * ДВА ЧИСЛА, И ПОРОЗНЬ НИ ОДНО НИЧЕГО НЕ ЗНАЧИТ:
 *   ЭКОНОМИЯ — отношение работы ограниченного прогона к неограниченному, и
 *     считается ДВОЯКО (А267): по установкам треугольника и по посещённым узлам
 *     дерева. Спуск дешевеет непосещением ПОДДЕРЕВЬЕВ, а установки этого не
 *     видят;
 *   ПОТЕРЯ — доля `Σf` ребёнка, съеденная ограничением: сумма дельта-форм-
 *     факторов по пикселям, где номер полигона стал другим.
 * Экономия, равная единице, значит, что заслонять нечем (случай зала).
 *
 * ЧТО ЗДЕСЬ СОЗНАТЕЛЬНО СДЕЛАНО ТАК, А НЕ ИНАЧЕ:
 *   - ограничение ставится ПО УЗЛАМ ДЕРЕВА двумя флагами, а не одним (А276):
 *     72.6 % установок зала приходится на треугольники ВНУТРЕННИХ узлов,
 *     поэтому узел без собственных видимых треугольников обязан остаться
 *     проходимым вниз;
 *   - пары с почти пустым полукубом ребёнка исключаются (А268): у крыши `Σf`
 *     близко к нулю, отнять нечего, и потеря выходит нулевой ПО ПОСТРОЕНИЮ;
 *   - улица ищется РАНЖИРОВАНИЕМ по замеренной глубине каньона, без порогов
 *     (А272);
 *   - в каждой строке таблицы печатается радиус элемента: сравнивать зал с
 *     городом по НОМЕРУ уровня — подмена (А269).
 */

#include "phcube.h"
#include "plink.h"
#include "plod.h"
#include "poly_seg.h"
#include "polygon.h"
#include "ptree.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ИМЕНОВАННЫЕ ЧИСЛА — все с происхождением, магических порогов нет (§4). */
#define HZ_PVIS_SEED 20260802u /* зерно ЛКГ: дата шага; фиксировано ради повторяемости */
#define HZ_PVIS_EYEH 1.7       /* высота глаза над землёй, м — рост человека */
#define HZ_PVIS_CELL 4.0       /* сторона ячейки поиска улицы, м: вчетверо крупнее грани города */
#define HZ_PVIS_RING_LO 2      /* кольцо замера стен, в ячейках: 8 м */
#define HZ_PVIS_RING_HI 4      /* ... до 16 м. Это МАСШТАБ ВЫБОРКИ, а не решающий порог */
/* Граница осмысленности ОТНОШЕНИЯ, а не порог качества (А268): ниже неё элемент
 * не получает света ни при каком наследовании, и делить на неё нечего. */
#define HZ_PVIS_SFMIN 0.01
#define HZ_PVIS_CHAINS 32 /* цепочек в выборке */
#define HZ_PVIS_PSAMP_MAX 5
/* Потолок сетки поиска улицы, ячеек по стороне. Не физическая величина, а
 * граница выделения: 2048 ячеек по 4 м — 8 км, вдвенадцатеро больше Rungholt. */
#define HZ_PVIS_GRIDMAX 2048

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static uint32_t lcg(uint32_t *s) {
  *s = *s * 1664525u + 1013904223u;
  return *s;
}

static int cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

/* Квантиль по УЖЕ ОТСОРТИРОВАННОМУ массиву. */
static double quant(const double *v, int n, double q) {
  if (n <= 0) return 0.0;
  int i = (int)(q * (double)(n - 1) + 0.5);
  if (i < 0) i = 0;
  if (i >= n) i = n - 1;
  return v[i];
}

/* ---------------------------------------------------------------- ПОИСК УЛИЦЫ */

typedef struct {
  double eye[3];
  double ground;
  double canyon; /* медиана максимума высоты в кольце минус земля, м */
} pv_spot;

/* Сетка максимумов высоты по (x, z). Для клетки УЛИЦЫ максимум и есть мостовая,
 * для клетки дома — конёк; на этой разнице и стоит ранжирование. */
static int pv_heightmap(const hz_objmesh *m, double **out, int *onx, int *onz) {
  /* РАЗМЕР ВЫДЕЛЕНИЯ НЕ ДОЛЖЕН ЗАВИСЕТЬ ОТ ПЛАВАЮЩЕЙ АРИФМЕТИКИ. Габарит даёт
   * лишь ЧИСЛО ячеек, и оно тут же зажимается в целочисленные границы: иначе
   * `gcc -fanalyzer` справедливо указывает, что размер malloc выведен из
   * `double` и точной оценке не поддаётся. */
  long lx = (long)((m->hi[0] - m->lo[0]) / HZ_PVIS_CELL) + 1;
  long lz = (long)((m->hi[2] - m->lo[2]) / HZ_PVIS_CELL) + 1;
  if (lx < 1) lx = 1;
  if (lz < 1) lz = 1;
  if (lx > HZ_PVIS_GRIDMAX) lx = HZ_PVIS_GRIDMAX;
  if (lz > HZ_PVIS_GRIDMAX) lz = HZ_PVIS_GRIDMAX;
  int nx = (int)lx, nz = (int)lz;
  double *g = malloc((size_t)nx * (size_t)nz * sizeof *g);
  if (g == NULL) return 2;
  const size_t ncell = (size_t)nx * (size_t)nz;
  for (size_t i = 0; i < ncell; i++)
    g[i] = -1e300;
  for (int32_t v = 0; v < m->nv; v++) {
    const double *p = m->v + 3 * (size_t)v;
    int i = (int)((p[0] - m->lo[0]) / HZ_PVIS_CELL);
    int k = (int)((p[2] - m->lo[2]) / HZ_PVIS_CELL);
    if (i < 0) i = 0;
    if (k < 0) k = 0;
    if (i >= nx) i = nx - 1;
    if (k >= nz) k = nz - 1;
    size_t o = (size_t)k * (size_t)nx + (size_t)i;
    if (p[1] > g[o]) g[o] = p[1];
  }
  *out = g;
  *onx = nx;
  *onz = nz;
  return 0;
}

/* Лучшая по глубине каньона клетка. Порогов нет: берётся МАКСИМУМ ранга. */
static int pv_find_street(const hz_objmesh *m, pv_spot *s, pv_spot *roof) {
  double *g = NULL;
  int nx = 0, nz = 0;
  if (pv_heightmap(m, &g, &nx, &nz) != 0) return 1;
  double best = -1e300, bestroof = -1e300;
  memset(s, 0, sizeof *s);
  memset(roof, 0, sizeof *roof);
  double ring[(2 * HZ_PVIS_RING_HI + 1) * (2 * HZ_PVIS_RING_HI + 1)];
  for (int k = 0; k < nz; k++)
    for (int i = 0; i < nx; i++) {
      double h0 = g[k * nx + i];
      if (h0 < -1e299) continue; /* пустая клетка — геометрии нет */
      if (h0 > bestroof) {
        bestroof = h0;
        roof->eye[0] = m->lo[0] + ((double)i + 0.5) * HZ_PVIS_CELL;
        roof->eye[1] = h0 + HZ_PVIS_EYEH;
        roof->eye[2] = m->lo[2] + ((double)k + 0.5) * HZ_PVIS_CELL;
        roof->ground = h0;
        roof->canyon = 0.0;
      }
      int nr = 0;
      for (int dk = -HZ_PVIS_RING_HI; dk <= HZ_PVIS_RING_HI; dk++)
        for (int di = -HZ_PVIS_RING_HI; di <= HZ_PVIS_RING_HI; di++) {
          int ad = (di < 0) ? -di : di, ak = (dk < 0) ? -dk : dk;
          int ch = (ad > ak) ? ad : ak;
          if (ch < HZ_PVIS_RING_LO || ch > HZ_PVIS_RING_HI) continue;
          int i2 = i + di, k2 = k + dk;
          if (i2 < 0 || k2 < 0 || i2 >= nx || k2 >= nz) continue;
          double hv = g[k2 * nx + i2];
          if (hv < -1e299) continue;
          ring[nr++] = hv;
        }
      if (nr < 8) continue;
      qsort(ring, (size_t)nr, sizeof ring[0], cmp_d);
      double dep = quant(ring, nr, 0.5) - h0;
      if (dep > best) {
        best = dep;
        s->eye[0] = m->lo[0] + ((double)i + 0.5) * HZ_PVIS_CELL;
        s->eye[1] = h0 + HZ_PVIS_EYEH;
        s->eye[2] = m->lo[2] + ((double)k + 0.5) * HZ_PVIS_CELL;
        s->ground = h0;
        s->canyon = dep;
      }
    }
  free(g);
  return (best > -1e299) ? 0 : 2;
}

/* Свидетельство «внутри улицы» ЗАМЕРОМ: полукуб нормалью ВВЕРХ. Доля неба есть
 * `Σ dff` по пустым пикселям; дальности по горизонту берутся с НИЖНЕЙ строки
 * боковых граней — она лежит на самой касательной плоскости, и никакого порога
 * «что считать горизонтом» не требуется. Возвращает направление вдоль улицы:
 * азимут с наибольшей медианной дальностью. */
static void pv_probe(hz_hcube *h, const hz_objmesh *m, const int32_t *t2p, const hz_ptree *T,
                     int32_t *stamp, int32_t mark, const double eye[3], double *sky, double *d10,
                     double *d50, double dirout[3]) {
  const double up[3] = {0.0, 1.0, 0.0};
  hz_hcube_stat st;
  memset(&st, 0, sizeof st);
  h->draw = NULL;
  h->desc = NULL;
  hz_hcube_draw_tree(h, m, t2p, T, eye, up, -1, stamp, mark, 0, &st);
  double s = 0.0, tot = 0.0;
  for (int i = 0; i < h->npix; i++) {
    double w = h->dff[i];
    if (!(w > 0.0)) continue;
    tot += w;
    if (h->id[i] < 0) s += w;
  }
  *sky = (tot > 0.0) ? s / tot : 0.0;
  /* нижние строки четырёх боковых граней: смещение верхней грани — R*R */
  int R = h->R;
  int nside = 4 * R;
  double *dd = malloc((size_t)nside * sizeof *dd);
  double bd = -1.0;
  dirout[0] = 1.0;
  dirout[1] = 0.0;
  dirout[2] = 0.0;
  int nd = 0;
  if (dd != NULL) {
    for (int f = 0; f < 4; f++) {
      int base = R * R + f * R * (R / 2);
      double side[4096];
      int ns = 0;
      for (int x = 0; x < R && ns < 4096; x++) {
        int o = base + x; /* строка v = 0 */
        if (h->id[o] >= 0 && h->depth[o] < 1e299) {
          dd[nd++] = h->depth[o];
          side[ns++] = h->depth[o];
        }
      }
      if (ns > 0) {
        qsort(side, (size_t)ns, sizeof side[0], cmp_d);
        double md = quant(side, ns, 0.5);
        if (md > bd) {
          bd = md;
          /* середина строки — направление вдоль оси грани */
          int o = base + R / 2;
          const double *dir = h->dir + 3 * (size_t)o;
          /* рама полукуба при нормали (0,1,0): ez = y, ex/ey — горизонтальные */
          dirout[0] = dir[0];
          dirout[1] = 0.0;
          dirout[2] = dir[1];
        }
      }
    }
    if (nd > 0) {
      qsort(dd, (size_t)nd, sizeof dd[0], cmp_d);
      *d10 = quant(dd, nd, 0.10);
      *d50 = quant(dd, nd, 0.50);
    } else {
      *d10 = *d50 = -1.0;
    }
    free(dd);
  } else {
    *d10 = *d50 = -1.0;
  }
  double nn = sqrt(dirout[0] * dirout[0] + dirout[2] * dirout[2]);
  if (nn > 0.0) {
    dirout[0] /= nn;
    dirout[2] /= nn;
  }
}

/* Пять граней полукуба одной картинкой — свидетельство улицы глазами. */
static int pv_write_hcube(const hz_hcube *h, const char *path) {
  int R = h->R, W = 5 * R, H = R;
  unsigned char *px = malloc(3 * (size_t)W * (size_t)H);
  if (px == NULL) return 1;
  memset(px, 0, 3 * (size_t)W * (size_t)H);
  double dmx = 0.0;
  for (int i = 0; i < h->npix; i++)
    if (h->id[i] >= 0 && h->depth[i] < 1e299 && h->depth[i] > dmx) dmx = h->depth[i];
  if (!(dmx > 0.0)) dmx = 1.0;
  for (int f = 0; f < 5; f++) {
    int fw = R, fh = (f == 0) ? R : R / 2;
    int base = (f == 0) ? 0 : R * R + (f - 1) * R * (R / 2);
    for (int y = 0; y < fh; y++)
      for (int x = 0; x < fw; x++) {
        int o = base + y * fw + x;
        int X = f * R + x, Y = H - 1 - y;
        size_t q = 3 * ((size_t)Y * (size_t)W + (size_t)X);
        if (h->id[o] < 0) { /* небо — синим */
          px[q] = 60;
          px[q + 1] = 110;
          px[q + 2] = 220;
        } else {
          double t = 1.0 - h->depth[o] / dmx;
          if (t < 0.0) t = 0.0;
          unsigned char g = (unsigned char)(30.0 + 225.0 * t);
          px[q] = g;
          px[q + 1] = g;
          px[q + 2] = g;
        }
      }
  }
  FILE *fp = fopen(path, "wb");
  if (fp == NULL) {
    free(px);
    return 2;
  }
  fprintf(fp, "P6\n%d %d\n255\n", W, H);
  fwrite(px, 1, 3 * (size_t)W * (size_t)H, fp);
  fclose(fp);
  free(px);
  return 0;
}

/* ------------------------------------------------------------ ЗАМЕР ПАРАМИ */

typedef struct {
  int level;      /* уровень РОДИТЕЛЯ */
  double rad;     /* радиус элемента родителя, м */
  double sf_full; /* Σf ребёнка без ограничения */
  double loss;    /* доля Σf, съеденная своим родителем */
  double loss_abs;
  double rsetup, rnode; /* экономия: отношение работы */
  double keepfrac;      /* доля узлов дерева в `desc` родителя */
  double wloss;         /* потеря при ЧУЖОМ родителе; <0 — контроль не ставился */
  double wdist;         /* расстояние между своим и чужим родителем, м */
  int nchild; /* детей у родителя: `1` — ТОЖДЕСТВЕННОЕ ПРОДВИЖЕНИЕ, пара ничего не доказывает */
  int used;
} pv_rec;

int main(int argc, char **argv) {
  int city = 0, scanonly = 0, lev = 9, ptleaf = 256, R = 64, psamp = 1, w = 1920;
  double trimax = 0.0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "city") == 0) city = 1;
    if (strcmp(argv[i], "scan") == 0) scanonly = 1;
    if (strncmp(argv[i], "lev=", 4) == 0) lev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "leaf=", 5) == 0) ptleaf = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "hemi=", 5) == 0) R = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "psamp=", 6) == 0) psamp = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "tri=", 4) == 0) trimax = strtod(argv[i] + 4, NULL);
    if (strncmp(argv[i], "w=", 2) == 0) w = (int)strtol(argv[i] + 2, NULL, 10);
  }
  if (psamp < 1) psamp = 1;
  if (psamp > HZ_PVIS_PSAMP_MAX) psamp = HZ_PVIS_PSAMP_MAX;
  if (!city && !(trimax > 0.0)) trimax = 0.25; /* зал: та же геометрия, что у эталона prad */

  double t0 = now_s();
  hz_objmesh m;
  if (hz_obj_load(&m, city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ,
                  city ? HZ_CFG_CITY_SCALE : HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  double t_load = now_s() - t0;
  if (trimax > 0.0 && hz_obj_subdivide(&m, trimax) != 0) return 1;
  printf("== СЦЕНА: %s, треугольников %d, габарит %.1f × %.1f × %.1f м, разбор %.2f с\n",
         city ? "ГОРОД" : "зал", m.nt, m.hi[0] - m.lo[0], m.hi[1] - m.lo[1], m.hi[2] - m.lo[2],
         t_load);

  /* Пространственное дерево нужно и поиску улицы, и замеру. */
  t0 = now_s();
  hz_ptree T;
  if (hz_ptree_build(&T, &m, ptleaf, 0) != 0) return 1;
  double t_tree = now_s() - t0;
  printf("== ДЕРЕВО: узлов %d, листьев %lld, ссылок %lld (избыточность %.2f), за %.1f с\n", T.nnd,
         (long long)T.nleaf, (long long)T.nref, T.redundancy, t_tree);

  /* Для полукуба нужен номер полигона на треугольник. До сегментации его нет,
   * а поиску улицы он и не нужен: ставим номер самого треугольника. */
  int32_t *t2p_tri = malloc((size_t)m.nt * sizeof *t2p_tri);
  int32_t *stamp = calloc((size_t)m.nt, sizeof *stamp);
  if (t2p_tri == NULL || stamp == NULL) return 1;
  for (int32_t t = 0; t < m.nt; t++)
    t2p_tri[t] = t;

  hz_hcube hp;
  if (hz_hcube_init(&hp, R) != 0) return 1;
  printf("== ПОЛУКУБ: R = %d, пикселей %d, Σ dff = %.15f (обязана быть 1)\n", hp.R, hp.npix,
         hp.sum);

  pv_spot street, roof;
  /* УЛИЦА МОЖЕТ НЕ НАЙТИСЬ, И ЭТО НЕ ОТКАЗ. Кольцо замера — 8…16 м, а зал
   * целиком 8 × 5 м: ни у одной ячейки нет кольца, и каньона там нет ПО
   * РАЗМЕРУ СЦЕНЫ. Тогда берётся камера §2 — и это ровно тот случай, ради
   * которого зал стоит контролем: заслонять нечем. */
  int street_ok = (pv_find_street(&m, &street, &roof) == 0);
  if (!street_ok) {
    double o[3] = HZ_CFG_HALL_EYE;
    for (int c = 0; c < 3; c++)
      street.eye[c] = o[c];
    street.ground = o[1] - HZ_PVIS_EYEH;
    street.canyon = 0.0;
    printf("== УЛИЦА НЕ НАЙДЕНА: сцена меньше кольца замера (%.0f…%.0f м). Берётся камера §2\n",
           HZ_PVIS_RING_LO * HZ_PVIS_CELL, HZ_PVIS_RING_HI * HZ_PVIS_CELL);
  }
  double dirs[3];
  double sky, d10, d50;
  int32_t mk = 1;
  pv_probe(&hp, &m, t2p_tri, &T, stamp, mk++, street.eye, &sky, &d10, &d50, dirs);
  printf("== УЛИЦА (ранг по глубине каньона, кольцо %.0f…%.0f м, порогов нет):\n"
         "   глаз (%.2f, %.2f, %.2f) м, земля %.2f м, глубина каньона %.2f м\n"
         "   доля неба %.4f; дальность по горизонту p10 %.2f м, p50 %.2f м\n",
         HZ_PVIS_RING_LO * HZ_PVIS_CELL, HZ_PVIS_RING_HI * HZ_PVIS_CELL, street.eye[0],
         street.eye[1], street.eye[2], street.ground, street.canyon, sky, d10, d50);
  {
    /* ИМЯ ПО СЦЕНЕ: прогон зала затирал городскую картинку. */
    char ipath[64];
    snprintf(ipath, sizeof ipath, "img/o43_street_%s.ppm", city ? "city" : "hall");
    if (pv_write_hcube(&hp, ipath) == 0)
      printf("   картинка: %s (пять граней, синее — небо)\n", ipath);
  }
  double at[3];
  for (int c = 0; c < 3; c++)
    at[c] = street.eye[c] + 50.0 * dirs[c];
  printf("   взгляд вдоль улицы: at (%.2f, %.2f, %.2f) м\n", at[0], at[1], at[2]);

  {
    double s2, a10, a50, dd[3];
    double old[3] = HZ_CFG_CITY_EYE, hall[3] = HZ_CFG_HALL_EYE;
    const double *o = city ? old : hall;
    double oe[3] = {o[0], o[1], o[2]};
    pv_probe(&hp, &m, t2p_tri, &T, stamp, mk++, oe, &s2, &a10, &a50, dd);
    printf("== КАМЕРА §2 (%.2f, %.2f, %.2f): доля неба %.4f; горизонт p10 %.2f м, p50 %.2f м\n",
           oe[0], oe[1], oe[2], s2, a10, a50);
    pv_probe(&hp, &m, t2p_tri, &T, stamp, mk++, roof.eye, &s2, &a10, &a50, dd);
    printf("== КРЫША (%.2f, %.2f, %.2f): доля неба %.4f; горизонт p10 %.2f м, p50 %.2f м\n",
           roof.eye[0], roof.eye[1], roof.eye[2], s2, a10, a50);
  }
  if (scanonly) {
    hz_hcube_free(&hp);
    hz_ptree_free(&T);
    free(t2p_tri);
    free(stamp);
    hz_obj_free(&m);
    return 0;
  }

  /* ---------------------------------------------------- конвейер до среза */
  double dseg = city ? 0.05 : 0.045;
  t0 = now_s();
  hz_pseglist sg;
  /* ПРЕДЕЛ ГАБАРИТА УЧАСТКА — ЧИСЛО ЗАЛА, И НА ГОРОД ЕГО ПЕРЕНОСИТЬ НЕЛЬЗЯ
   * (А285, поймано прогоном). `0.5` м взято из `prad`, где зал перед этим дробят
   * до стороны `0.25` м. У Rungholt грань — единичный квадрат, и предел `0.5` м
   * мельче ВХОДНОГО треугольника: ни один участок не вмещает даже одного, и
   * сегментация выдаёт по участку на треугольник — `6 704 264` вместо `446 013`.
   * Отсюда и брались два часа лестницы. У города предел снят: входная
   * тесселяция и без того однородна (площадь грани `0.5015 ± 0.1` м², §105.2). */
  const double segcap = city ? 0.0 : 0.5;
  if (hz_seg_planar_cap2(&sg, &m, dseg, segcap, NULL) != 0) return 1;
  double t_seg = now_s() - t0;
  t0 = now_s();
  hz_polyset ps;
  if (hz_poly_build(&ps, &m, &sg) != 0) return 1;
  double t_poly = now_s() - t0;
  printf("== СЕГМЕНТАЦИЯ: участков %d (предел габарита %.2f м) за %.1f с; полигоны за %.1f с\n",
         sg.nseg, segcap, t_seg, t_poly);
  /* СТОРОЖ ВЫРОЖДЕНИЯ (А285). Участок на треугольник значит, что сегментация не
   * сработала вовсе, и всё, что ниже, будет мерить не сцену, а вход. Отказывать
   * надо здесь, а не выяснять это через два часа по профилю. */
  if (sg.nseg >= m.nt) {
    fprintf(stderr,
            "участков %d при %d треугольниках — сегментация выродилась; "
            "предел габарита %.2f м мельче входного треугольника\n",
            sg.nseg, m.nt, segcap);
    return 1;
  }

  t0 = now_s();
  hz_lod L;
  hz_lodcfg lc;
  memset(&lc, 0, sizeof lc);
  lc.delta0 = dseg;
  const double eps = (HZ_CFG_FOV_DEG * M_PI / 180.0) / (double)w;
  lc.eps = eps;
  lc.maxlev = lev;
  lc.radmul = 4.0;
  lc.base = 1.4142;
  if (hz_lod_build_merge(&L, &m, &sg, &ps, &lc) != 0) {
    fprintf(stderr, "отказ лестницы\n");
    return 1;
  }
  double t_lod = now_s() - t0;
  printf("== ЛЕСТНИЦА: уровней %d, узлов %d, за %.1f с (ε = %.3e рад)\n", L.nlev, L.nnd, t_lod,
         eps);

  int32_t *cut = malloc((size_t)L.np * sizeof *cut);
  double *pt = malloc(3 * (size_t)L.nnd * sizeof *pt);
  int32_t *t2p = malloc((size_t)m.nt * sizeof *t2p);
  if (cut == NULL || pt == NULL || t2p == NULL) return 1;
  t0 = now_s();
  int32_t ncut = hz_lod_cut(&L, street.eye, eps, 1, cut);
  if (hz_links_points(pt, &L, &ps, &m, 0) != 0) return 1;
  for (int32_t k = 0; k < ps.np; k++)
    for (int32_t t = ps.p[k].t0; t < ps.p[k].t0 + ps.p[k].ntri; t++)
      t2p[ps.tri[t]] = k;
  printf("== СРЕЗ (глаз на улице): %d узлов из %d, за %.1f с\n", ncut, L.nnd, now_s() - t0);

  /* дети узлов лестницы — CSR, чтобы `psamp > 1` не сканировал всю лестницу */
  int32_t *ccnt = calloc((size_t)L.nnd + 1, sizeof *ccnt);
  int32_t *cidx = malloc((size_t)L.nnd * sizeof *cidx);
  if (ccnt == NULL || cidx == NULL) return 1;
  for (int32_t i = 0; i < L.nnd; i++)
    if (L.nd[i].parent >= 0) ccnt[L.nd[i].parent + 1]++;
  for (int32_t i = 0; i < L.nnd; i++)
    ccnt[i + 1] += ccnt[i];
  {
    int32_t *fill = malloc((size_t)L.nnd * sizeof *fill);
    if (fill == NULL) return 1;
    memcpy(fill, ccnt, (size_t)L.nnd * sizeof *fill);
    for (int32_t i = 0; i < L.nnd; i++)
      if (L.nd[i].parent >= 0) cidx[fill[L.nd[i].parent]++] = i;
    free(fill);
  }

  /* ---------------------------------------------------------- выборка цепочек */
  int32_t *uniq = malloc((size_t)ncut * sizeof *uniq);
  unsigned char *seen = calloc((size_t)L.nnd, 1);
  if (uniq == NULL || seen == NULL) return 1;
  int32_t nu = 0;
  for (int32_t k = 0; k < L.np; k++) {
    int32_t nd = cut[k];
    if (nd >= 0 && nd < L.nnd && !seen[nd]) {
      seen[nd] = 1;
      uniq[nu++] = nd;
    }
  }
  if (nu <= 0) {
    fprintf(stderr, "срез пуст — выбирать не из чего\n");
    return 1;
  }
  const int K = HZ_PVIS_CHAINS;
  int32_t chain[HZ_PVIS_CHAINS][32];
  int clen[HZ_PVIS_CHAINS];
  uint32_t rs = HZ_PVIS_SEED;
  for (int i = 0; i < K; i++) {
    int32_t nd = uniq[(int32_t)(lcg(&rs) % (uint32_t)nu)];
    int n = 0;
    while (nd >= 0 && n < 32) {
      chain[i][n++] = nd;
      nd = L.nd[nd].parent;
    }
    clen[i] = n;
  }

  /* ПРОХОД 1: маски родителей, по цепочке и уровню. */
  int nlev = L.nlev;
  size_t msz = (size_t)T.nnd;
  unsigned char *mdraw = calloc((size_t)K * (size_t)nlev * msz, 1);
  unsigned char *mdesc = calloc((size_t)K * (size_t)nlev * msz, 1);
  unsigned char *mhave = calloc((size_t)K * (size_t)nlev, 1);
  int32_t *mnode = malloc((size_t)K * (size_t)nlev * sizeof *mnode);
  if (mdraw == NULL || mdesc == NULL || mhave == NULL || mnode == NULL) {
    fprintf(stderr, "не хватило памяти под маски\n");
    return 1;
  }
  t0 = now_s();
#pragma omp parallel
  {
    hz_hcube h;
    int32_t *st2 = calloc((size_t)m.nt, sizeof *st2);
    int ok = (st2 != NULL && hz_hcube_init(&h, R) == 0);
    hz_hcube_stat st;
#pragma omp for schedule(dynamic, 1)
    for (int i = 0; i < K; i++) {
      if (!ok) continue;
      int32_t mark = (int32_t)(i + 1) * 1000;
      for (int j = 1; j < clen[i]; j++) { /* j = 0 — сам узел среза, он ребёнок */
        int32_t p = chain[i][j];
        int lv = L.nd[p].level;
        if (lv < 0 || lv >= nlev) continue;
        size_t off = ((size_t)i * (size_t)nlev + (size_t)lv) * msz;
        for (int s = 0; s < psamp; s++) {
          const double *x = pt + 3 * (size_t)p;
          if (s > 0) {
            int32_t c0 = ccnt[p], c1 = ccnt[p + 1];
            int32_t nc = c1 - c0;
            if (nc <= 0) break;
            x = pt + 3 * (size_t)cidx[c0 + (int32_t)((int64_t)s * nc / psamp) % nc];
          }
          memset(&st, 0, sizeof st);
          h.draw = NULL;
          h.desc = NULL;
          hz_hcube_draw_tree(&h, &m, t2p, &T, x, L.nd[p].n, -1, st2, mark++, 0, &st);
          hz_hcube_mark(&h, mdraw + off);
        }
        hz_hcube_close(&T, mdraw + off, mdesc + off);
        mhave[(size_t)i * (size_t)nlev + (size_t)lv] = 1;
        mnode[(size_t)i * (size_t)nlev + (size_t)lv] = p;
      }
    }
    if (ok) hz_hcube_free(&h);
    free(st2);
  }
  double t_p1 = now_s() - t0;

  /* ПРОХОД 2: дети — полный, ограниченный своим, ограниченный чужим. */
  int npair = K * nlev;
  pv_rec *rec = calloc((size_t)npair, sizeof *rec);
  int32_t *idfull = NULL;
  int svc_diff = -1; /* служебный контроль: `keep` из одних единиц */
  if (rec == NULL) return 1;
  t0 = now_s();
#pragma omp parallel
  {
    hz_hcube h;
    int32_t *st2 = calloc((size_t)m.nt, sizeof *st2);
    int32_t *idf = malloc((size_t)hp.npix * sizeof *idf);
    int ok = (st2 != NULL && idf != NULL && hz_hcube_init(&h, R) == 0);
    hz_hcube_stat sa, sb;
#pragma omp for schedule(dynamic, 1)
    for (int i = 0; i < K; i++) {
      if (!ok) continue;
      int32_t mark = (int32_t)(i + 1) * 100000;
      for (int j = 0; j + 1 < clen[i]; j++) {
        int32_t c = chain[i][j], p = chain[i][j + 1];
        int lv = L.nd[p].level;
        if (lv < 0 || lv >= nlev || !mhave[(size_t)i * (size_t)nlev + (size_t)lv]) continue;
        size_t off = ((size_t)i * (size_t)nlev + (size_t)lv) * msz;
        pv_rec *r = &rec[i * nlev + lv];
        r->level = lv;
        r->rad = L.nd[p].rad;
        r->wloss = -1.0;
        /* полный */
        memset(&sa, 0, sizeof sa);
        h.draw = NULL;
        h.desc = NULL;
        hz_hcube_draw_tree(&h, &m, t2p, &T, pt + 3 * (size_t)c, L.nd[c].n, -1, st2, mark++, 0, &sa);
        double sf = hz_hcube_sf(&h);
        memcpy(idf, h.id, (size_t)h.npix * sizeof *idf);
        r->sf_full = sf;
        if (!(sf >= HZ_PVIS_SFMIN)) continue; /* А268: делить не на что */
        /* ограниченный СВОИМ родителем */
        memset(&sb, 0, sizeof sb);
        h.draw = mdraw + off;
        h.desc = mdesc + off;
        hz_hcube_draw_tree(&h, &m, t2p, &T, pt + 3 * (size_t)c, L.nd[c].n, -1, st2, mark++, 0, &sb);
        double ls = 0.0;
        for (int q = 0; q < h.npix; q++)
          if (idf[q] >= 0 && h.id[q] != idf[q]) ls += (double)h.dff[q];
        r->loss_abs = ls;
        r->loss = ls / sf;
        r->rsetup = (sa.nsetup > 0) ? (double)sb.nsetup / (double)sa.nsetup : 1.0;
        r->rnode = (sa.nnode > 0) ? (double)sb.nnode / (double)sa.nnode : 1.0;
        int64_t nk = 0;
        for (int32_t q = 0; q < T.nnd; q++)
          nk += (mdesc + off)[q];
        r->keepfrac = (double)nk / (double)T.nnd;
        r->nchild = (int)(ccnt[p + 1] - ccnt[p]);
        r->used = 1;
        /* НЕГАТИВНЫЙ КОНТРОЛЬ: чужой родитель того же уровня */
        int iw = (i + K / 2) % K;
        if (mhave[(size_t)iw * (size_t)nlev + (size_t)lv]) {
          size_t ow = ((size_t)iw * (size_t)nlev + (size_t)lv) * msz;
          memset(&sb, 0, sizeof sb);
          h.draw = mdraw + ow;
          h.desc = mdesc + ow;
          hz_hcube_draw_tree(&h, &m, t2p, &T, pt + 3 * (size_t)c, L.nd[c].n, -1, st2, mark++, 0,
                             &sb);
          double lw = 0.0;
          for (int q = 0; q < h.npix; q++)
            if (idf[q] >= 0 && h.id[q] != idf[q]) lw += (double)h.dff[q];
          r->wloss = lw / sf;
          int32_t pw = mnode[(size_t)iw * (size_t)nlev + (size_t)lv];
          double dsum = 0.0;
          for (size_t cc = 0; cc < 3; cc++) {
            double dv = pt[3 * (size_t)p + cc] - pt[3 * (size_t)pw + cc];
            dsum += dv * dv;
          }
          r->wdist = sqrt(dsum);
        }
      }
    }
    if (ok) hz_hcube_free(&h);
    free(st2);
    free(idf);
  }
  double t_p2 = now_s() - t0;

  /* СЛУЖЕБНЫЙ КОНТРОЛЬ (предсказан УСПЕШНЫМ, потому и не негативный): маска из
   * одних единиц обязана дать буфер номеров ПОБИТОВО тот же, что без маски. */
  {
    unsigned char *ones = malloc(msz), *od = malloc(msz);
    idfull = malloc((size_t)hp.npix * sizeof *idfull);
    int32_t *st2 = calloc((size_t)m.nt, sizeof *st2);
    if (ones != NULL && od != NULL && idfull != NULL && st2 != NULL) {
      memset(ones, 1, msz);
      memset(od, 1, msz);
      hz_hcube_stat sc;
      int32_t nd = chain[0][0];
      memset(&sc, 0, sizeof sc);
      hp.draw = NULL;
      hp.desc = NULL;
      hz_hcube_draw_tree(&hp, &m, t2p, &T, pt + 3 * (size_t)nd, L.nd[nd].n, -1, st2, 1, 0, &sc);
      memcpy(idfull, hp.id, (size_t)hp.npix * sizeof *idfull);
      memset(&sc, 0, sizeof sc);
      hp.draw = ones;
      hp.desc = od;
      hz_hcube_draw_tree(&hp, &m, t2p, &T, pt + 3 * (size_t)nd, L.nd[nd].n, -1, st2, 2, 0, &sc);
      svc_diff = 0;
      for (int q = 0; q < hp.npix; q++)
        if (hp.id[q] != idfull[q]) svc_diff++;
      hp.draw = NULL;
      hp.desc = NULL;
    }
    free(ones);
    free(od);
    free(st2);
  }

  /* --------------------------------------------------------------- доклад */
  printf("== ВРЕМЯ ЗАМЕРА: маски родителей %.1f с, дети %.1f с; цепочек %d, psamp %d\n", t_p1, t_p2,
         K, psamp);
  printf("== СЛУЖЕБНЫЙ КОНТРОЛЬ (маска из единиц против прогона без маски): "
         "различных пикселей %d (обязан быть 0)\n",
         svc_diff);

  printf("== НАСЛЕДОВАНИЕ ВИДИМОСТИ, %s, psamp = %d\n", city ? "ГОРОД (улица)" : "зал", psamp);
  printf("  ур  пар  rad p50, м  keep p50  ЭКОНОМИЯ setup p10/p50/p90   ЭКОНОМИЯ node p50   "
         "ПОТЕРЯ p50/p90/max   потеря абс p90\n");
  int nskip = 0;
  double *b1 = malloc((size_t)npair * sizeof *b1), *b2 = malloc((size_t)npair * sizeof *b2);
  double *b3 = malloc((size_t)npair * sizeof *b3), *b4 = malloc((size_t)npair * sizeof *b4);
  double *b5 = malloc((size_t)npair * sizeof *b5), *b6 = malloc((size_t)npair * sizeof *b6);
  if (b1 == NULL || b2 == NULL || b3 == NULL || b4 == NULL || b5 == NULL || b6 == NULL) return 1;
  for (int lv = 0; lv < nlev; lv++) {
    int n = 0;
    for (int i = 0; i < K; i++) {
      pv_rec *r = &rec[i * nlev + lv];
      if (!r->used) continue;
      /* ТОЖДЕСТВЕННОЕ ПРОДВИЖЕНИЕ ИСКЛЮЧАЕТСЯ (А277, найдено прогоном зала).
       * У узла с ОДНИМ ребёнком родитель и ребёнок — одна и та же поверхность с
       * одной опорной точкой, поэтому наследование там точно по построению и
       * ничего не доказывает: ложный ноль. */
      if (r->nchild <= 1) continue;
      b1[n] = r->rsetup;
      b2[n] = r->loss;
      b3[n] = r->rad;
      b4[n] = r->keepfrac;
      b5[n] = r->rnode;
      b6[n] = r->loss_abs;
      n++;
    }
    if (n == 0) continue;
    qsort(b1, (size_t)n, sizeof *b1, cmp_d);
    qsort(b2, (size_t)n, sizeof *b2, cmp_d);
    qsort(b3, (size_t)n, sizeof *b3, cmp_d);
    qsort(b4, (size_t)n, sizeof *b4, cmp_d);
    qsort(b5, (size_t)n, sizeof *b5, cmp_d);
    qsort(b6, (size_t)n, sizeof *b6, cmp_d);
    printf("  %2d  %3d   %8.2f  %8.4f      %.3f/%.3f/%.3f            %.3f          "
           "%.4f/%.4f/%.4f      %.4f\n",
           lv, n, quant(b3, n, 0.5), quant(b4, n, 0.5), quant(b1, n, 0.10), quant(b1, n, 0.50),
           quant(b1, n, 0.90), quant(b5, n, 0.50), quant(b2, n, 0.50), quant(b2, n, 0.90),
           b2[n - 1], quant(b6, n, 0.90));
  }
  int nident = 0;
  for (int i = 0; i < npair; i++) {
    if (rec[i].sf_full > 0.0 && !rec[i].used) nskip++;
    if (rec[i].used && rec[i].nchild <= 1) nident++;
  }
  printf("   пар исключено: по Σf < %.2f (А268) — %d; ТОЖДЕСТВЕННЫХ ПРОДВИЖЕНИЙ (А277) — %d\n",
         HZ_PVIS_SFMIN, nskip, nident);

  printf("== НЕГАТИВНЫЙ КОНТРОЛЬ (ЧУЖОЙ РОДИТЕЛЬ), %s\n", city ? "ГОРОД (улица)" : "зал");
  printf("  ур  пар  ПОТЕРЯ p10/p50/p90   расстояние до чужого p50, м\n");
  for (int lv = 0; lv < nlev; lv++) {
    int n = 0;
    for (int i = 0; i < K; i++) {
      pv_rec *r = &rec[i * nlev + lv];
      if (!r->used || r->wloss < 0.0 || r->nchild <= 1) continue;
      b1[n] = r->wloss;
      b2[n] = r->wdist;
      n++;
    }
    if (n == 0) continue;
    qsort(b1, (size_t)n, sizeof *b1, cmp_d);
    qsort(b2, (size_t)n, sizeof *b2, cmp_d);
    printf("  %2d  %3d      %.4f/%.4f/%.4f          %8.2f\n", lv, n, quant(b1, n, 0.10),
           quant(b1, n, 0.50), quant(b1, n, 0.90), quant(b2, n, 0.5));
  }

  free(b1);
  free(b2);
  free(b3);
  free(b4);
  free(b5);
  free(b6);
  free(rec);
  free(idfull);
  free(mdraw);
  free(mdesc);
  free(mhave);
  free(mnode);
  free(uniq);
  free(seen);
  free(ccnt);
  free(cidx);
  free(cut);
  free(pt);
  free(t2p);
  free(t2p_tri);
  free(stamp);
  hz_hcube_free(&hp);
  hz_ptree_free(&T);
  hz_lod_free(&L);
  hz_poly_free(&ps);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
