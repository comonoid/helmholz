/* pgather — СБОР ПО ПИКСЕЛЮ НА НОВОМ НОСИТЕЛЕ (§847): первая картинка
 * линии пирамиды + объёмного фронта.
 *
 * Камера — ОПЕРАТОР над решённым полем: поле считается свипом mode=2
 * (§845/§846, build=1), затем на каждый пиксель пускается луч из глаза,
 * ищется БЛИЖАЙШЕЕ пересечение с куском, яркость пикселя = L_out куска:
 *
 *     L_out = Le + rho·E_куска/(2π)     (ламберт, А1471 — двусторонняя)
 *
 * фон — 0 (тёмное окружение, §842). Пересечение — перебором ВСЕХ кусков
 * с точным ближайшим t: на синтетике (2.2–17.7 тыс. кусков) это честно;
 * марш луча по пирамиде — следующий шаг (А1537-класс, §847).
 *
 * Приёмки §847: К27 на двух касающихся шарах (лучи в диск касания — в
 * БЛИЖНИЙ); сличение сборщика со свипом (средняя яркость cavity =
 * Le + rho·E_avg/2π, ±5 %); НК le=0 / noprop → кадр тождественно 0.
 *
 * Ключи:
 *   lev=N dirs=6|26|NxM it=N rho=F le=F tau0 noprop — как в swee3;
 *   mode=2 build=1 — только объёмный фронт (умолчания);
 *   eye=X,Y,Z look=X,Y,Z fov=F W=N H=N — камера (fov в градусах, по
 *   горизонтали); out=ФАЙЛ — вывод (умолчание img/pgather.ppm, PPM P6);
 *   k27=N — К27-приёмка: N×N контрольных лучей через диск касания шаров.
 *   useke — §874/А1576: эмиссия материалов Ke попадает в свип и кадр;
 *   ksf=F — §874: дихотомия T4 (зеркало): F·ср.ks3 в свип (хоп walk) и
 *           зеркальные вторичные лучи камеры (глубина ≤ 4); 0 — прежний мир.
 *
 * КАРТИНКИ ПИШУТСЯ В img/, А НЕ В build/ (правило проекта).
 */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/resource.h>
#include <fcntl.h>

#include "nstruct/pyr.h"
#include "nstruct/sweep.h"
#include "scene_obj.h"

#define PG_FOV_MIN 5.0   /* ниже — телеобъектив вне смысла синтетических тестов */
#define PG_FOV_MAX 170.0 /* выше — кадр шире полусферы, проекция вырождается */
#define PG_WH_MAX 4096   /* потолок кадра: перебор nt на луч, больше не нужно */
#define PG_K27_RAYS 64   /* контрольных лучей на сторону диска касания */
/* §874: глубина зеркальной рекурсии камеры = HZ_MIRROR_BOUNCE_MAX из walk
 * (sweep.c, §873 шаг 2): 0.95^4 < 0.82 — глубже камера не видит того,
 * чего нет в поле. */
#define PG_MIRROR_BOUNCE_MAX 4
/* §874/А1582: сдвиг начала вторичного луча вдоль R против самопересечения,
 * в долях габарита сцены (меньше толщины стен синтетики) */
#define PG_HOP_EPS_REL 1e-6

static double now_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int parse3(const char *s, double out[3]) {
  return sscanf(s, "%lf,%lf,%lf", &out[0], &out[1], &out[2]) == 3 ? 0 : 1;
}

/* пересечение луча (o, d — единичный) с треугольником (Мёллер—Трумбор);
 * возврат t ≥ tmin или -1; p — вершины [3][3] */
static double pg_ray_tri(const double o[3], const double d[3], const double p[3][3]) {
  double e1[3], e2[3], pv[3], tv[3], qv[3];
  double det, u, v, t;
  int ax;
  for (ax = 0; ax < 3; ax++) {
    e1[ax] = p[1][ax] - p[0][ax];
    e2[ax] = p[2][ax] - p[0][ax];
  }
  pv[0] = d[1] * e2[2] - d[2] * e2[1];
  pv[1] = d[2] * e2[0] - d[0] * e2[2];
  pv[2] = d[0] * e2[1] - d[1] * e2[0];
  det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
  if (fabs(det) < 1e-30) return -1.0; /* луч в плоскости треугольника */
  u = 1.0 / det;
  tv[0] = o[0] - p[0][0];
  tv[1] = o[1] - p[0][1];
  tv[2] = o[2] - p[0][2];
  u = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * u;
  if (u < 0.0 || u > 1.0) return -1.0;
  qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
  qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
  qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
  v = (d[0] * qv[0] + d[1] * qv[1] + d[2] * qv[2]) / det;
  if (v < 0.0 || u + v > 1.0) return -1.0;
  t = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) / det;
  if (t < 1e-9) return -1.0;
  return t;
}

/* ---- §849: МАРШ ЛУЧА ПО ПИРАМИДЕ (bbox-CSR + 3D-DDA) --------------------- */

/* bbox-CSR: для каждого листа — куски, чей bbox пересекает клетку.
 * Владельческий CSR пирамиды для сбора ДЫРЯВ (А1543): треугольник краем
 * в соседней клетке. Консервативен только bbox-признак. */
typedef struct {
  int64_t *start; /* [nleaf+1] */
  int32_t *pids;  /* [ntotal] */
  int64_t ntotal;
} pg_bbox_csr;

static void pg_bbox_span(const hz_pyr *py, const double cmin[3], const double cmax[3],
                         int64_t i0[3], int64_t i1[3]) {
  int32_t ax;
  i0[0] = i0[1] = i0[2] = 0; /* анализатор теряет индукцию цикла по ax (FP-класс diam) */
  i1[0] = i1[1] = i1[2] = 0;
  for (ax = 0; ax < 3; ax++) {
    int64_t n = ax == 0 ? py->nx : (ax == 1 ? py->ny : py->nz);
    i0[ax] = (int64_t)((cmin[ax] - py->lo[ax]) / py->cell);
    i1[ax] = (int64_t)((cmax[ax] - py->lo[ax]) / py->cell);
    if (i0[ax] < 0) i0[ax] = 0;
    if (i1[ax] >= n) i1[ax] = n - 1;
  }
}

static int64_t pg_cell_id(const hz_pyr *py, int64_t i, int64_t j, int64_t k) {
  return i + py->nx * (j + py->ny * k);
}

static void pg_bbox_csr_free(pg_bbox_csr *csr) {
  free(csr->start);
  free(csr->pids);
  memset(csr, 0, sizeof *csr);
}

/* построение: два прохода counting-sort по всем клеткам bbox кусков */
static int pg_bbox_csr_build(const hz_pyr *py, const double *cmin, const double *cmax,
                             pg_bbox_csr *csr) {
  int64_t *cnt;
  int32_t p, li;
  memset(csr, 0, sizeof *csr);
  cnt = (int64_t *)calloc((size_t)py->nleaf, sizeof *cnt);
  csr->start = (int64_t *)calloc((size_t)py->nleaf + 1, sizeof *csr->start);
  if (!cnt || !csr->start) {
    free(cnt);
    pg_bbox_csr_free(csr);
    return 2;
  }
  for (p = 0; p < py->nt; p++) {
    int64_t i0[3], i1[3], i, j, k;
    int32_t tri = py->pcs[p].tri; /* cmin/cmax не переставлялись — индекс по tri */
    pg_bbox_span(py, cmin + 3 * (int64_t)tri, cmax + 3 * (int64_t)tri, i0, i1);
    for (k = i0[2]; k <= i1[2]; k++)
      for (j = i0[1]; j <= i1[1]; j++)
        for (i = i0[0]; i <= i1[0]; i++) {
          int32_t pos = hz_pyr_leaf_pos(py, pg_cell_id(py, i, j, k));
          if (pos >= 0) cnt[pos]++;
        }
  }
  for (li = 0; li < py->nleaf; li++)
    csr->start[li + 1] = csr->start[li] + cnt[li];
  csr->ntotal = csr->start[py->nleaf];
  csr->pids = (int32_t *)malloc((size_t)(csr->ntotal > 0 ? csr->ntotal : 1) * sizeof *csr->pids);
  if (!csr->pids) {
    free(cnt);
    pg_bbox_csr_free(csr);
    return 2;
  }
  for (li = 0; li < py->nleaf; li++)
    cnt[li] = csr->start[li]; /* переиспользован как курсор заполнения */
  for (p = 0; p < py->nt; p++) {
    int64_t i0[3], i1[3], i, j, k;
    int32_t tri = py->pcs[p].tri;
    pg_bbox_span(py, cmin + 3 * (int64_t)tri, cmax + 3 * (int64_t)tri, i0, i1);
    for (k = i0[2]; k <= i1[2]; k++)
      for (j = i0[1]; j <= i1[1]; j++)
        for (i = i0[0]; i <= i1[0]; i++) {
          int32_t pos = hz_pyr_leaf_pos(py, pg_cell_id(py, i, j, k));
          if (pos >= 0) csr->pids[cnt[pos]++] = p;
        }
  }
  free(cnt);
  return 0;
}

/* DDA-сбор: ближайший кусок вдоль луча (Аманатидес—Ву; срезка А1544,
 * многократная проверка куска безвредна — строгое «меньше», А1545,
 * тай-брейк по порядку осей А1546). Приборы — steps/tested.
 * Возврат куска или -1. */
static int32_t pg_dda(const hz_pyr *py, const pg_bbox_csr *csr, const hz_objmesh *m,
                      const double eye[3], const double rd[3], int64_t *steps, int64_t *tested,
                      double *thit) {
  double tlo = 0.0, thi = 1e30, tcur, tbest = -1.0;
  int64_t cell[3], stepv[3];
  double tnext[3], tdelta[3];
  int32_t ax, best = -1;
  for (ax = 0; ax < 3; ax++) {
    double n = ax == 0 ? (double)py->nx : (ax == 1 ? (double)py->ny : (double)py->nz);
    if (fabs(rd[ax]) < 1e-30) {
      if (eye[ax] < py->lo[ax] || eye[ax] > py->lo[ax] + n * py->cell) return -1;
      continue;
    }
    {
      double ta = (py->lo[ax] - eye[ax]) / rd[ax];
      double tb = (py->lo[ax] + n * py->cell - eye[ax]) / rd[ax];
      if (ta > tb) {
        double tt = ta;
        ta = tb;
        tb = tt;
      }
      if (ta > tlo) tlo = ta;
      if (tb < thi) thi = tb;
    }
  }
  if (tlo > thi) return -1; /* мимо сцены (А1544) */
  tcur = tlo;
  for (ax = 0; ax < 3; ax++) {
    int64_t nax = ax == 0 ? py->nx : (ax == 1 ? py->ny : py->nz);
    double pos = eye[ax] + tcur * rd[ax];
    int64_t idx = (int64_t)((pos - py->lo[ax]) / py->cell);
    if (idx < 0) idx = 0;
    if (idx >= nax) idx = nax - 1;
    cell[ax] = idx;
    if (fabs(rd[ax]) < 1e-30) {
      tnext[ax] = 1e30;
      tdelta[ax] = 0.0;
      stepv[ax] = 0;
    } else {
      double boundary = rd[ax] > 0.0 ? py->lo[ax] + ((double)idx + 1.0) * py->cell
                                     : py->lo[ax] + (double)idx * py->cell;
      stepv[ax] = rd[ax] > 0.0 ? 1 : -1;
      tdelta[ax] = py->cell / fabs(rd[ax]);
      tnext[ax] = tcur + (boundary - pos) / rd[ax];
    }
  }
  for (;;) {
    int32_t pos = hz_pyr_leaf_pos(py, pg_cell_id(py, cell[0], cell[1], cell[2]));
    double tn;
    (*steps)++;
    if (pos >= 0) {
      int64_t s;
      for (s = csr->start[pos]; s < csr->start[pos + 1]; s++) {
        int32_t p = csr->pids[s];
        double p3[3][3], tt;
        (*tested)++;
        hz_obj_tri(m, py->pcs[p].tri, p3);
        tt = pg_ray_tri(eye, rd, p3);
        if (tt >= 0.0 && (tbest < 0.0 || tt < tbest)) {
          tbest = tt;
          best = p;
        }
      }
    }
    /* ранний выход: вход в следующую клетку дальше ближайшего попадания */
    tn = tnext[0];
    if (tnext[1] < tn) tn = tnext[1];
    if (tnext[2] < tn) tn = tnext[2];
    if (tn > thi || (tbest >= 0.0 && tn > tbest)) break;
    if (tnext[0] <= tnext[1] && tnext[0] <= tnext[2]) {
      tcur = tnext[0];
      tnext[0] += tdelta[0];
      cell[0] += stepv[0];
    } else if (tnext[1] <= tnext[2]) {
      tcur = tnext[1];
      tnext[1] += tdelta[1];
      cell[1] += stepv[1];
    } else {
      tcur = tnext[2];
      tnext[2] += tdelta[2];
      cell[2] += stepv[2];
    }
    if (cell[0] < 0 || cell[1] < 0 || cell[2] < 0 || cell[0] >= py->nx || cell[1] >= py->ny ||
        cell[2] >= py->nz)
      break; /* страховка выхода из сетки */
  }
  if (thit) *thit = tbest; /* §874: параметр удара — вторичный луч (А1585) */
  return best;
}

/* §874/А1585: общий сборщик ближайшего попадания — базовый и вторичный лучи
 * ходят ОДНИМ кодом (gather=0 — DDA по пирамиде, gather=1 — перебор).
 * Возврат куска или -1; thit — параметр удара. */
static int32_t pg_nearest(const hz_pyr *py, const pg_bbox_csr *csr, const hz_objmesh *m,
                          const double org[3], const double rd[3], int gather_mode, double *thit,
                          int64_t *steps, int64_t *tested) {
  if (gather_mode == 0) return pg_dda(py, csr, m, org, rd, steps, tested, thit);
  {
    int32_t i, best = -1;
    double tbest = -1.0;
    for (i = 0; i < m->nt; i++) {
      double p3[3][3], tt;
      hz_obj_tri(m, py->pcs[i].tri, p3);
      tt = pg_ray_tri(org, rd, p3);
      if (tt >= 0.0 && (tbest < 0.0 || tt < tbest)) {
        tbest = tt;
        best = i;
      }
    }
    *thit = tbest;
    return best;
  }
}

/* §874: яркость по лучу С отражениями (дихотомия T4).
 * L_cam(p,d) = le + kdvis·E_p/2π + ks_eff·L_cam(отражённый, d+1),
 * ks_eff = min(ks[p], 1−kdvis) — тот же клэмп, что в депозите walk
 * (§873; А1580: kdvis по ovr-семантике, как front_rho). */
typedef struct {
  const hz_pyr *py;
  const pg_bbox_csr *csr;
  const hz_objmesh *m;
  const double *kd;  /* переставлен, индекс по куску */
  const double *lep; /* NULL — нет per-piece эмиссии */
  const double *ks;  /* NULL — дихотомия выключена (ksf=0) */
  const double *nrm;
  double le, rho;
  int gather;
  double hop_eps;
  int64_t *steps, *tested;
  int64_t *nsec; /* зеркальных вторичных лучей (прибор) */
} pg_cam;

static double pg_lcam_hit(const pg_cam *c, const double org[3], const double rd[3], int32_t p,
                          double thit, int depth) {
  double kdvis = c->rho < 0 ? c->kd[p] : c->rho;
  double L = c->le + (c->lep ? c->lep[p] : 0.0) + kdvis * (double)c->py->pcs[p].e / (2.0 * M_PI);
  if (c->ks && depth < PG_MIRROR_BOUNCE_MAX) {
    double kse = c->ks[p];
    if (kse > 1.0 - kdvis) kse = 1.0 - kdvis > 0.0 ? 1.0 - kdvis : 0.0;
    if (kse > 0.0) {
      const double *nv = c->nrm + 3 * (int64_t)p;
      double dot = rd[0] * nv[0] + rd[1] * nv[1] + rd[2] * nv[2];
      double rorg[3], rr[3], sth = -1.0;
      int32_t q;
      int ax;
      /* отражение не зависит от ориентации нормали: R = rd − 2(rd·n̂)n̂ */
      for (ax = 0; ax < 3; ax++)
        rr[ax] = rd[ax] - 2.0 * dot * nv[ax];
      for (ax = 0; ax < 3; ax++)
        rorg[ax] = org[ax] + rd[ax] * thit + rr[ax] * c->hop_eps;
      q = pg_nearest(c->py, c->csr, c->m, rorg, rr, c->gather, &sth, c->steps, c->tested);
      (*c->nsec)++;
      if (q >= 0) L += kse * pg_lcam_hit(c, rorg, rr, q, sth, depth + 1);
    }
  }
  return L;
}

/* --- §882: HBLK v1 — mmap-сбор (блок = листовая клетка, Morton) --- */
#define PG_HBLK_MAGIC 0x314B4C4248ULL
typedef struct {
  uint8_t *base; /* начало mmap */
  size_t len;
  const uint64_t *hdr;
  const uint8_t *tab; /* pb_cellrec {cell,start,cnt+pad} 24 Б */
  const int32_t *ids;
  const double *tris; /* 9 double на tri, исходный порядок */
  int64_t nocc, nx, ny, nz, nt;
  double cell;
  double lo[3];
  int nlev; /* занятых уровней бинпоиска */
} pg_hblk;

static int64_t pg_hblk_morton3(uint32_t x, uint32_t y, uint32_t z) {
  int64_t r = 0;
  int b;
  for (b = 0; b < 21; b++)
    r |= ((int64_t)(x >> b & 1) << (3 * b)) | ((int64_t)(y >> b & 1) << (3 * b + 1)) |
         ((int64_t)(z >> b & 1) << (3 * b + 2));
  return r;
}

static int pg_hblk_open(pg_hblk *h, const char *path) {
  int fd = open(path, O_RDONLY);
  if (fd < 0) return 2;
  off_t len = lseek(fd, 0, SEEK_END);
  if (len < 24 * 8) {
    close(fd);
    return 2;
  }
  h->len = (size_t)len;
  h->base = (uint8_t *)mmap(NULL, h->len, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  if (h->base == MAP_FAILED) return 2;
  madvise(h->base, h->len, MADV_RANDOM);
  h->hdr = (const uint64_t *)h->base;
  if (h->hdr[0] != PG_HBLK_MAGIC) return 2;
  h->nt = (int64_t)h->hdr[1];
  h->nx = (int64_t)h->hdr[2];
  h->ny = (int64_t)h->hdr[3];
  h->nz = (int64_t)h->hdr[4];
  memcpy(h->lo, &h->hdr[5], 3 * sizeof(double));
  memcpy(&h->cell, &h->hdr[8], sizeof(double));
  h->nocc = (int64_t)h->hdr[9];
  h->tab = h->base + h->hdr[10];
  h->ids = (const int32_t *)(h->base + h->hdr[11]);
  h->tris = (const double *)(h->base + h->hdr[12]);
  h->nlev = 0;
  while ((1LL << h->nlev) < h->nocc)
    h->nlev++;
  return 0;
}

static void pg_hblk_close(pg_hblk *h) {
  munmap(h->base, h->len);
}

/* ближайшее пересечение луча с блочным файлом; куски/вершины — из mmap */
static int32_t pg_hblk_nearest(const pg_hblk *h, const double org[3], const double rd[3],
                               int64_t *tested) {
  double tlo = 0.0, thi = 1e30, tcur, tbest = -1.0;
  int64_t cellv[3], stepv[3];
  double tnext[3], tdelta[3];
  int32_t ax, best = -1;
  for (ax = 0; ax < 3; ax++) {
    double n = ax == 0 ? (double)h->nx : (ax == 1 ? (double)h->ny : (double)h->nz);
    if (fabs(rd[ax]) < 1e-30) {
      if (org[ax] < h->lo[ax] || org[ax] > h->lo[ax] + n * h->cell) return -1;
      continue;
    }
    {
      double ta = (h->lo[ax] - org[ax]) / rd[ax];
      double tb = (h->lo[ax] + n * h->cell - org[ax]) / rd[ax];
      if (ta > tb) {
        double tt = ta;
        ta = tb;
        tb = tt;
      }
      if (ta > tlo) tlo = ta;
      if (tb < thi) thi = tb;
    }
  }
  if (tlo > thi) return -1;
  tcur = tlo;
  for (ax = 0; ax < 3; ax++) {
    int64_t nax = ax == 0 ? h->nx : (ax == 1 ? h->ny : h->nz);
    double pos = org[ax] + tcur * rd[ax];
    int64_t idx = (int64_t)((pos - h->lo[ax]) / h->cell);
    if (idx < 0) idx = 0;
    if (idx >= nax) idx = nax - 1;
    cellv[ax] = idx;
    if (fabs(rd[ax]) < 1e-30) {
      tnext[ax] = 1e30;
      tdelta[ax] = 0.0;
      stepv[ax] = 0;
    } else {
      double boundary = rd[ax] > 0.0 ? h->lo[ax] + ((double)idx + 1.0) * h->cell
                                     : h->lo[ax] + (double)idx * h->cell;
      stepv[ax] = rd[ax] > 0.0 ? 1 : -1;
      tdelta[ax] = h->cell / fabs(rd[ax]);
      tnext[ax] = tcur + (boundary - pos) / rd[ax];
    }
  }
  for (;;) {
    /* занятость клетки — бинпоиск Morton-ключа по таблице (А1623) */
    int64_t key = pg_hblk_morton3((uint32_t)cellv[0], (uint32_t)cellv[1], (uint32_t)cellv[2]);
    int64_t lo = 0, hi = h->nocc;
    while (lo < hi) {
      int64_t mid = (lo + hi) / 2;
      int64_t c;
      memcpy(&c, h->tab + (size_t)mid * 24, 8);
      if (c < key)
        lo = mid + 1;
      else
        hi = mid;
    }
    if (lo < h->nocc) {
      int64_t c;
      memcpy(&c, h->tab + (size_t)lo * 24, 8);
      if (c == key) {
        int64_t start = 0;
        int32_t cnt32 = 0; /* §882: частичное memcpy в int64 оставляет старшие
                            * байты отравленными (ASAN 0xBE) — SEGV; читать
                            * ровно int32 и расширять */
        memcpy(&start, h->tab + (size_t)lo * 24 + 8, 8);
        memcpy(&cnt32, h->tab + (size_t)lo * 24 + 16, 4);
        int64_t cnt = cnt32;
        for (int64_t s = 0; s < cnt; s++) {
          int32_t t = h->ids[start + s];
          const double *tv = h->tris + 9 * (int64_t)t;
          double p3[3][3], tt;
          (*tested)++;
          for (int q = 0; q < 9; q++)
            ((double *)p3)[q] = tv[q];
          tt = pg_ray_tri(org, rd, p3);
          if (tt >= 0.0 && (tbest < 0.0 || tt < tbest)) {
            tbest = tt;
            best = t; /* исходный tri — kd/ks/pg_ray_tri согласованы */
          }
        }
      }
    }
    double tn = tnext[0];
    if (tnext[1] < tn) tn = tnext[1];
    if (tnext[2] < tn) tn = tnext[2];
    if (tn > thi || (tbest >= 0.0 && tn > tbest)) break;
    if (tnext[0] <= tnext[1] && tnext[0] <= tnext[2]) {
      tcur = tnext[0];
      tnext[0] += tdelta[0];
      cellv[0] += stepv[0];
    } else if (tnext[1] <= tnext[2]) {
      tcur = tnext[1];
      tnext[1] += tdelta[1];
      cellv[1] += stepv[1];
    } else {
      tcur = tnext[2];
      tnext[2] += tdelta[2];
      cellv[2] += stepv[2];
    }
    if (cellv[0] < 0 || cellv[1] < 0 || cellv[2] < 0 || cellv[0] >= h->nx || cellv[1] >= h->ny ||
        cellv[2] >= h->nz)
      break;
  }
  return best;
}

int main(int argc, char **argv) {
  const char *path = NULL, *outfile = "img/pgather.ppm";
  /* §875: дефолты потребителя — ПРОДАКШН-МОДЕЛЬ (дихотомия ks §873/874,
   * per-piece эмиссия Ke А1576); useke=0 / ksf=0 — документированный выход
   * к прежнему миру. */
  double scale = 1.0, le = 1.0, rho = -1.0, fov = 60.0, ksf = 1.0;
  int useke = 1;
  double eye[3] = {0, 0, 0}, look[3] = {0, 0, 0};
  int iters = 30, lev = 6, tau0 = 0, noprop = 0, ndirs = 26, mort = 1, i, ax;
  double *lep = NULL; /* §874: per-piece эмиссия (ср. Ke); NULL — прежний мир */
  int W = 320, H = 240, k27 = 0;
  int frames = 1, have_eye2 = 0; /* §877: ходьба */
  double eye2[3] = {0, 0, 0}, look2[3] = {0, 0, 0};
  int have_delbox = 0;        /* §879: разрушаемость-прототип */
  int rgb = 0;                /* §889: RGB-рендер */
  double expmul = 1.0;        /* §890: множитель экспозиции */
  const char *blkfile = NULL; /* §882: HBLK v1, mmap-сбор */
  double delbox[6];
  hz_objmesh m;
  hz_pyr py;
  hz_sw_opts so;
  hz_sw_stat st;
  double *area = NULL, *nrm = NULL, *kd = NULL, *cent = NULL, *cmin = NULL, *cmax = NULL;
  double *ks = NULL; /* §874: per-piece зеркальная доля (ksf · ср. ks3); NULL-семантика в so */
  double maxdim = 0.0;
  int32_t *mtl = NULL;
  double cell, t0, t1, sw_time;
  double *lum = NULL;
  /* §867: данные mode=3, заполняются до memset(&so) — см. блок trivert */
  double *g_tv9 = NULL, *g_tb6 = NULL;
  uint8_t *g_lparr = NULL;
  int gather = 0; /* §849: 0 — DDA (умолчание), 1 — brute (путь верификации) */
  pg_bbox_csr csr;
  int64_t dda_steps = 0, dda_tested = 0, sec_rays = 0;

  for (i = 1; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0)
      lev = atoi(argv[i] + 4);
    else if (strncmp(argv[i], "rho=", 4) == 0)
      rho = atof(argv[i] + 4);
    else if (strncmp(argv[i], "le=", 3) == 0)
      le = atof(argv[i] + 3);
    else if (strncmp(argv[i], "it=", 3) == 0)
      iters = atoi(argv[i] + 3);
    else if (strncmp(argv[i], "dirs=", 5) == 0) {
      int a, b;
      if (sscanf(argv[i] + 5, "%dx%d", &a, &b) == 2 && a > 0 && b > 0)
        ndirs = a * 100 + b;
      else
        ndirs = atoi(argv[i] + 5);
    } else if (strncmp(argv[i], "eye=", 4) == 0)
      parse3(argv[i] + 4, eye);
    else if (strncmp(argv[i], "look=", 5) == 0)
      parse3(argv[i] + 5, look);
    else if (strncmp(argv[i], "fov=", 4) == 0)
      fov = atof(argv[i] + 4);
    else if (strncmp(argv[i], "W=", 2) == 0)
      W = atoi(argv[i] + 2);
    else if (strncmp(argv[i], "H=", 2) == 0)
      H = atoi(argv[i] + 2);
    else if (strcmp(argv[i], "useke") == 0)
      useke = 1; /* А1576; голый флаг — обратная совместимость §874 (А1587) */
    else if (strncmp(argv[i], "useke=", 6) == 0)
      useke = atoi(argv[i] + 6); /* §875: useke=0 — прежний мир */
    else if (strncmp(argv[i], "ksf=", 4) == 0)
      ksf = atof(argv[i] + 4); /* §874: дихотомия T4 в pgather */
    else if (strncmp(argv[i], "gather=", 7) == 0)
      gather = atoi(argv[i] + 7);
    else if (strncmp(argv[i], "k27=", 4) == 0)
      k27 = atoi(argv[i] + 4);
    else if (strncmp(argv[i], "expm=", 5) == 0)
      expmul = atof(argv[i] + 5); /* §890 */
    else if (strncmp(argv[i], "rgb=", 4) == 0)
      rgb = atoi(argv[i] + 4); /* §889 */
    else if (strncmp(argv[i], "blk=", 4) == 0)
      blkfile = argv[i] + 4; /* §882 */
    else if (strncmp(argv[i], "delbox=", 7) == 0) {
      if (sscanf(argv[i] + 7, "%lf,%lf,%lf:%lf,%lf,%lf", &delbox[0], &delbox[1], &delbox[2],
                 &delbox[3], &delbox[4], &delbox[5]) == 6)
        have_delbox = 1; /* §879 */
    } else if (strncmp(argv[i], "frames=", 7) == 0)
      frames = atoi(argv[i] + 7); /* §877 */
    else if (strncmp(argv[i], "eye2=", 5) == 0) {
      parse3(argv[i] + 5, eye2);
      have_eye2 = 1;
    } else if (strncmp(argv[i], "look2=", 6) == 0)
      parse3(argv[i] + 6, look2);
    else if (strncmp(argv[i], "out=", 4) == 0)
      outfile = argv[i] + 4;
    else if (strcmp(argv[i], "tau0") == 0)
      tau0 = 1;
    else if (strcmp(argv[i], "noprop") == 0)
      noprop = 1;
    else if (strncmp(argv[i], "mort=", 5) == 0)
      mort = atoi(argv[i] + 5);
    else if (strncmp(argv[i], "scale=", 6) == 0)
      scale = atof(argv[i] + 6);
    else
      path = argv[i];
  }
  if (!path || W < 1 || H < 1 || W > PG_WH_MAX || H > PG_WH_MAX || fov < PG_FOV_MIN ||
      fov > PG_FOV_MAX) {
    fprintf(stderr, "use: pgather <scene.obj> [lev=N dirs=.. it=N rho=F le=F] "
                    "[eye=X,Y,Z look=X,Y,Z fov=F W=N H=N] [k27=N] [out=ФАЙЛ]\n");
    return 2;
  }
  t0 = now_sec();
  if (!have_eye2)
    for (ax = 0; ax < 3; ax++)
      eye2[ax] = eye[ax]; /* §877: поворот на месте */
  if (frames < 1 || (have_eye2 && frames < 2)) {
    fprintf(stderr, "pgather: frames=N>=2 требует eye2=/look2= (ходьба §877)\n");
    return 2;
  }
  t0 = now_sec();
  if (hz_obj_load(&m, path, scale) != 0) {
    fprintf(stderr, "pgather: не читается %s\n", path);
    return 2;
  }
  t1 = now_sec();
  printf("СТАТЬЯ obj-загрузка: %.2f с (nt=%d, mtl=%d)\n", t1 - t0, m.nt, m.nmtl);

  if (blkfile) {
    /* --- §882 v1: mmap-сбор, самодостаточный (без пирамиды/свипа).
     * E≡0 (слот E в файле — v2): кадр = силуэт (le на попаданиях). */
    pg_hblk hb;
    struct rusage ra, rb;
    if (pg_hblk_open(&hb, blkfile) != 0) {
      fprintf(stderr, "pgather: HBLK не читается: %s\n", blkfile);
      return 2;
    }
    printf("HBLK: nt=%lld, сетка %lldx%lldx%lld, cell=%.4g, занятых клеток %lld\n",
           (long long)hb.nt, (long long)hb.nx, (long long)hb.ny, (long long)hb.nz, hb.cell,
           (long long)hb.nocc);
    getrusage(RUSAGE_SELF, &ra);
    {
      double fwd[3], right[3], up[3], tmpv[3] = {0, 0, 1};
      double tanf = tan(fov * M_PI / 360.0);
      int64_t nhit2 = 0, tested2 = 0;
      double tA = now_sec(), tB;
      for (ax = 0; ax < 3; ax++)
        fwd[ax] = look[ax] - eye[ax];
      {
        double nn = sqrt(fwd[0] * fwd[0] + fwd[1] * fwd[1] + fwd[2] * fwd[2]);
        if (nn < 1e-12) {
          fprintf(stderr, "pgather: глаз совпадает с точкой взгляда\n");
          return 2;
        }
        for (ax = 0; ax < 3; ax++)
          fwd[ax] /= nn;
      }
      for (ax = 0; ax < 3; ax++)
        right[ax] = fwd[(ax + 1) % 3] * tmpv[(ax + 2) % 3] - fwd[(ax + 2) % 3] * tmpv[(ax + 1) % 3];
      {
        double nn = sqrt(right[0] * right[0] + right[1] * right[1] + right[2] * right[2]);
        if (nn < 1e-9) {
          tmpv[0] = 1;
          tmpv[1] = tmpv[2] = 0;
          for (ax = 0; ax < 3; ax++)
            right[ax] =
                fwd[(ax + 1) % 3] * tmpv[(ax + 2) % 3] - fwd[(ax + 2) % 3] * tmpv[(ax + 1) % 3];
          nn = sqrt(right[0] * right[0] + right[1] * right[1] + right[2] * right[2]);
        }
        for (ax = 0; ax < 3; ax++)
          right[ax] /= nn;
      }
      for (ax = 0; ax < 3; ax++)
        up[ax] = fwd[(ax + 1) % 3] * right[(ax + 2) % 3] - fwd[(ax + 2) % 3] * right[(ax + 1) % 3];
      lum = (double *)malloc((size_t)W * (size_t)H * sizeof *lum);
      if (!lum) return 2;
      /* §886: параллелизм §878/§884 ВКЛЮЧЁН (указание: делать все
       * оптимизации) — строки кадра независимы, файл read-only */
#pragma omp parallel for schedule(dynamic, 16) private(ax) reduction(+ : tested2, nhit2)
      for (int iy = 0; iy < H; iy++)
        for (int ix = 0; ix < W; ix++) {
          double sx = (2.0 * (ix + 0.5) / W - 1.0) * tanf;
          double sy = (1.0 - 2.0 * (iy + 0.5) / H) * tanf * (double)H / (double)W;
          double rd[3] = {0, 0, 0}, nn;
          int32_t hit;
          int64_t tested_loc = 0;
          for (ax = 0; ax < 3; ax++)
            rd[ax] = fwd[ax] + sx * right[ax] + sy * up[ax];
          nn = sqrt(rd[0] * rd[0] + rd[1] * rd[1] + rd[2] * rd[2]);
          for (ax = 0; ax < 3; ax++)
            rd[ax] /= nn;
          hit = pg_hblk_nearest(&hb, eye, rd, &tested_loc);
          tested2 += tested_loc;
          if (hit >= 0) nhit2++;
          lum[(size_t)iy * (size_t)W + (size_t)ix] = hit >= 0 ? le : 0.0;
        }
      tB = now_sec();
      getrusage(RUSAGE_SELF, &rb);
      printf("HBLK СБОР: лучей %d, попало %" PRId64 " (%.2f %%), кусков/луч %.1f, %.2f с\n", W * H,
             nhit2, 100.0 * (double)nhit2 / ((double)W * (double)H),
             (double)tested2 / ((double)W * (double)H), tB - tA);
      printf("HBLK RSS: max %.1f МБ; фолты минорные %lld, мажорные %lld (за сбор)\n",
             (double)rb.ru_maxrss / 1024.0, (long long)(rb.ru_minflt - ra.ru_minflt),
             (long long)(rb.ru_majflt - ra.ru_majflt));
      {
        char fname[4096];
        FILE *f;
        snprintf(fname, sizeof fname, "%s", outfile);
        f = fopen(fname, "wb");
        if (!f) return 2;
        fprintf(f, "P6\n%d %d\n255\n", W, H);
        for (i = 0; i < W * H; i++) {
          unsigned char b[3];
          b[0] = b[1] = b[2] = (unsigned char)(lum[i] > 0 ? 255 : 0);
          fwrite(b, 1, 3, f);
        }
        fclose(f);
        printf("КАДР: %s записан (HBLK v1, силуэт)\n", fname);
      }
      free(lum);
    }
    pg_hblk_close(&hb);
    hz_obj_free(&m);
    return 0;
  }

  area = (double *)malloc((size_t)m.nt * sizeof *area);
  nrm = (double *)malloc((size_t)m.nt * 3 * sizeof *nrm);
  kd = (double *)malloc((size_t)m.nt * sizeof *kd);
  ks = (double *)malloc((size_t)m.nt * sizeof *ks); /* §874 */
  if (useke) lep = (double *)malloc((size_t)m.nt * sizeof *lep);
  cent = (double *)malloc((size_t)m.nt * 3 * sizeof *cent);
  cmin = (double *)malloc((size_t)m.nt * 3 * sizeof *cmin);
  cmax = (double *)malloc((size_t)m.nt * 3 * sizeof *cmax);
  mtl = (int32_t *)malloc((size_t)m.nt * sizeof *mtl);
  if (!area || !nrm || !kd || !ks || (useke && !lep) || !cent || !cmin || !cmax || !mtl) {
    fprintf(stderr, "pgather: нет памяти\n");
    return 2;
  }
  for (i = 0; i < m.nt; i++) {
    double p[3][3], nn;
    /* явная инициализация — анализатор теряет индукцию цикла по ax
     * (FP-класс diam, прецедент pg_bbox_span); все элементы далее
     * перезаписываются — арифметика не меняется */
    double e1[3] = {0, 0, 0}, e2[3] = {0, 0, 0};
    int v;
    hz_obj_tri(&m, (int32_t)i, p);
    for (ax = 0; ax < 3; ax++) {
      double lo = p[0][ax], hi = p[0][ax];
      for (v = 1; v < 3; v++) {
        if (p[v][ax] < lo) lo = p[v][ax];
        if (p[v][ax] > hi) hi = p[v][ax];
      }
      cmin[3 * (int64_t)i + ax] = lo;
      cmax[3 * (int64_t)i + ax] = hi;
      cent[3 * (int64_t)i + ax] = (p[0][ax] + p[1][ax] + p[2][ax]) / 3.0;
    }
    for (ax = 0; ax < 3; ax++) {
      e1[ax] = p[1][ax] - p[0][ax];
      e2[ax] = p[2][ax] - p[0][ax];
    }
    nrm[3 * (int64_t)i] = e1[1] * e2[2] - e1[2] * e2[1];
    nrm[3 * (int64_t)i + 1] = e1[2] * e2[0] - e1[0] * e2[2];
    nrm[3 * (int64_t)i + 2] = e1[0] * e2[1] - e1[1] * e2[0];
    nn = sqrt(nrm[3 * (int64_t)i] * nrm[3 * (int64_t)i] +
              nrm[3 * (int64_t)i + 1] * nrm[3 * (int64_t)i + 1] +
              nrm[3 * (int64_t)i + 2] * nrm[3 * (int64_t)i + 2]);
    area[i] = 0.5 * nn;
    if (nn > 0)
      for (ax = 0; ax < 3; ax++)
        nrm[3 * (int64_t)i + ax] /= nn;
    kd[i] = m.mtl[m.fm[i]].kd;
    /* §874: ks — ksf · среднее ks3 (прецедент усреднения kd, §873 шаг 1);
     * заполнение в порядке ИСХОДНЫХ треугольников ДО hz_pyr_build —
     * перестановка выровняет ks со слотами кусков (урок §873-Ф-а) */
    ks[i] = ksf * ((m.mtl[m.fm[i]].ks3[0] + m.mtl[m.fm[i]].ks3[1] + m.mtl[m.fm[i]].ks3[2]) / 3.0);
    if (lep) /* §874: front_le = le + lep[p] — в lep только ср. Ke (без le) */
      lep[i] = (m.mtl[m.fm[i]].ke3[0] + m.mtl[m.fm[i]].ke3[1] + m.mtl[m.fm[i]].ke3[2]) / 3.0;
    mtl[i] = m.fm[i];
  }

  t1 = now_sec();
  printf("СТАТЬЯ подготовка (area/nrm/kd/ks/lep/trivert): %.2f с\n", t1 - t0);
  t0 = now_sec();
  {
    for (ax = 0; ax < 3; ax++) {
      double s = m.hi[ax] - m.lo[ax];
      if (s > maxdim) maxdim = s;
    }
    cell = maxdim / (double)(1 << lev);
  }
  { /* §867: trivert/tribox/lp для mode=3 (точное пересечение + walk) —
     * в порядке ИСХОДНЫХ треугольников (pcs[].tri); так как so ниже
     * memset-ится, указатели назначаются ПОСЛЕ memset (класс бага pref) */
    double *tv9 = (double *)malloc((size_t)m.nt * 9 * sizeof *tv9);
    double *tb6 = (double *)malloc((size_t)m.nt * 6 * sizeof *tb6);
    uint8_t *lparr = (uint8_t *)malloc((size_t)m.nt);
    int32_t ti;
    if (!tv9 || !tb6 || !lparr) {
      fprintf(stderr, "нет памяти на trivert/tribox/lp\n");
      return 2;
    }
    for (ti = 0; ti < m.nt; ti++) {
      double pp[3][3];
      int32_t a2, aa;
      hz_obj_tri(&m, ti, pp);
      for (aa = 0; aa < 9; aa++)
        tv9[9 * (int64_t)ti + aa] = ((const double *)pp)[aa];
      for (a2 = 0; a2 < 3; a2++) {
        double lo = pp[0][a2], hi2 = pp[0][a2];
        int v;
        for (v = 1; v < 3; v++) {
          if (pp[v][a2] < lo) lo = pp[v][a2];
          if (pp[v][a2] > hi2) hi2 = pp[v][a2];
        }
        tb6[6 * (int64_t)ti + a2] = lo;
        tb6[6 * (int64_t)ti + 3 + a2] = hi2;
      }
      lparr[ti] = 0; /* ℓ_p = 0: рабочий уровень листьев (как в pref) */
    }
    g_tv9 = tv9;
    g_tb6 = tb6;
    g_lparr = lparr;
  }
  if (hz_pyr_build(&py, m.nt, cmin, cmax, cent, mtl, m.lo, m.hi, cell) != 0) {
    fprintf(stderr, "pgather: пирамида не построилась\n");
    return 2;
  }
  if (mort &&
      (hz_pyr_morton(&py) != 0 || hz_pyr_permute(&py, area, sizeof *area) != 0 ||
       hz_pyr_permute(&py, nrm, 3 * sizeof *nrm) != 0 || hz_pyr_permute(&py, kd, sizeof *kd) != 0 ||
       hz_pyr_permute(&py, ks, sizeof *ks) != 0 ||
       (lep && hz_pyr_permute(&py, lep, sizeof *lep) != 0))) { /* §874/А1581 */
    fprintf(stderr, "pgather: Morton не прошёл\n");
    return 2;
  }

  { /* §867: per-node max ℓ_p обязателен для mode=3 (§852) */
    int32_t nup = 0;
    if (hz_pyr_set_lp(&py, g_lparr, &nup) != 0) {
      fprintf(stderr, "pgather: per-node max lp не построился\n");
      return 2;
    }
  }
  t1 = now_sec();
  printf("СТАТЬЯ пирамида+Morton+lp: %.2f с (листьев %d)\n", t1 - t0, py.nleaf);
  t0 = now_sec();
  /* --- СВИП: поле E на кусках (фронт mode=3 + walk, §867) --- */
  if (gather == 0 && pg_bbox_csr_build(&py, cmin, cmax, &csr) != 0) {
    fprintf(stderr, "pgather: bbox-CSR не построился\n");
    return 2;
  }
  t1 = now_sec();
  if (gather == 0) printf("СТАТЬЯ bbox-CSR: %.2f с\n", t1 - t0);
  memset(&so, 0, sizeof so);
  so.le = le;
  so.rho = rho;
  so.iters = iters;
  so.tau0 = tau0;
  so.noprop = noprop;
  so.ndirs = ndirs;
  so.mode = 3; /* §867: фронт с точным пересечением (был mode=2) */
  so.build = 1;
  so.vc = 1;
  so.walk = 1;                   /* §867: продакшн-модель А1576/§861 — дефолт потребителя */
  so.ks = ksf > 0.0 ? ks : NULL; /* §874: ksf=0 — прежний мир (битово, П1) */
  so.lep = lep;                  /* §874: NULL при useke=0 — побитово прежний мир */
  so.trivert = g_tv9;
  so.tribox = g_tb6;
  so.lp = g_lparr;
  so.domhi = m.hi;
  double *Ec[3] = {NULL, NULL, NULL}; /* §889: E по каналам */
  double *kdc[3] = {NULL, NULL, NULL}, *lepc[3] = {NULL, NULL, NULL};
  if (rgb) {
    /* §889: три скалярных решения с альбедо/эмиссией канала (классика
     * радиосити). kd_c = kd3[c], lep_c = ke3[c]; le остаётся общим. */
    if (ksf > 0.0) {
      fprintf(stderr, "pgather: rgb=1 с ksf>0 не совмещается в v1 - отказ\n");
      return 2;
    }
    for (int ch = 0; ch < 3; ch++) {
      kdc[ch] = (double *)malloc((size_t)m.nt * sizeof *kdc[ch]);
      lepc[ch] = (double *)malloc((size_t)m.nt * sizeof *lepc[ch]);
      Ec[ch] = (double *)malloc((size_t)m.nt * sizeof *Ec[ch]);
      if (!kdc[ch] || !lepc[ch] || !Ec[ch]) return 2;
      for (i = 0; i < m.nt; i++) {
        kdc[ch][i] = m.mtl[m.fm[i]].kd3[ch];
        lepc[ch][i] = m.mtl[m.fm[i]].ke3[ch];
      }
    }
    t0 = now_sec();
    for (int ch = 0; ch < 3; ch++) {
      so.lep = lepc[ch]; /* kd канала идёт 3-м аргументом hz_sw_run */
      for (i = 0; i < m.nt; i++)
        py.pcs[i].e = 0.0f;
      if (hz_sw_run(&py, m.nt, area, nrm, kdc[ch], &so, &st, NULL) != 0) return 2;
      for (i = 0; i < m.nt; i++)
        Ec[ch][i] = py.pcs[i].e;
    }
    t1 = now_sec();
    sw_time = t1 - t0;
    printf("СВИП RGB: 3 канала (%.3f с), E_avg=%.4f\n", sw_time, st.e_avg);
  } else {
    for (i = 0; i < m.nt; i++)
      py.pcs[i].e = 0.0f;
    t0 = now_sec();
    if (hz_sw_run(&py, m.nt, area, nrm, kd, &so, &st, NULL) != 0) {
      fprintf(stderr, "pgather: свип не прошёл\n");
      return 2;
    }
    t1 = now_sec();
    sw_time = t1 - t0;
    printf("СВИП: E_avg=%.4f (%.3f с)\n", st.e_avg, sw_time);
  }

  if (have_delbox) {
    /* --- §879: РАЗРУШАЕМОСТЬ-1. A: решение полной сцены (выше), E по tri;
     * B: пересборка без delbox-кусков + тёплый пересвип; C: холодный (E=0),
     * те же итерации. pcs.tri остаётся исходным tri (А1604). */
    double *oldE = (double *)calloc((size_t)m.nt, sizeof *oldE);
    int *keep = (int *)calloc((size_t)m.nt, sizeof *keep);
    int32_t nt2 = 0, p2;
    double *cmin2, *cmax2, *cent2, *area2, *nrm2, *kd2, *ks2, *lep2 = NULL, *warmE;
    int32_t *mtl2, *Ltri;
    hz_sw_opts soB;
    hz_sw_stat stB, stC;
    double tB, tC, dsum = 0, dmax = 0, esum = 0;
    if (!oldE || !keep) {
      fprintf(stderr, "pgather: нет памяти на разрушение\n");
      return 2;
    }
    for (i = 0; i < m.nt; i++)
      oldE[py.pcs[i].tri] = py.pcs[i].e; /* А1608 */
    for (i = 0; i < m.nt; i++) {
      int out = cent[3 * (int64_t)i] >= delbox[0] && cent[3 * (int64_t)i] <= delbox[3] &&
                cent[3 * (int64_t)i + 1] >= delbox[1] && cent[3 * (int64_t)i + 1] <= delbox[4] &&
                cent[3 * (int64_t)i + 2] >= delbox[2] && cent[3 * (int64_t)i + 2] <= delbox[5];
      keep[i] = !out;
      if (!out) nt2++;
    }
    printf("РАЗРУШЕНИЕ: удалено %d из %d кусков\n", m.nt - nt2, m.nt);
    if (nt2 == 0) { /* А1607: keep-all допустим (НК); недопустимо удаление
                     * ВСЕЙ сцены — решать нечего */
      fprintf(stderr, "pgather: delbox удаляет всю сцену\n");
      return 2;
    }
    cmin2 = (double *)malloc((size_t)nt2 * 3 * sizeof *cmin2);
    cmax2 = (double *)malloc((size_t)nt2 * 3 * sizeof *cmax2);
    cent2 = (double *)malloc((size_t)nt2 * 3 * sizeof *cent2);
    area2 = (double *)malloc((size_t)nt2 * sizeof *area2);
    nrm2 = (double *)malloc((size_t)nt2 * 3 * sizeof *nrm2);
    kd2 = (double *)malloc((size_t)nt2 * sizeof *kd2);
    ks2 = (double *)malloc((size_t)nt2 * sizeof *ks2);
    mtl2 = (int32_t *)malloc((size_t)nt2 * sizeof *mtl2);
    Ltri = (int32_t *)malloc((size_t)nt2 * sizeof *Ltri);
    if (useke) lep2 = (double *)malloc((size_t)nt2 * sizeof *lep2);
    warmE = (double *)malloc((size_t)nt2 * sizeof *warmE);
    if (!cmin2 || !cmax2 || !cent2 || !area2 || !nrm2 || !kd2 || !ks2 || !mtl2 || !Ltri ||
        (useke && !lep2) || !warmE) {
      fprintf(stderr, "pgather: нет памяти на фазы B/C\n");
      return 2;
    }
    p2 = 0;
    for (i = 0; i < m.nt; i++) { /* компактация в порядке кусков фазы A */
      int32_t tri = py.pcs[i].tri;
      if (!keep[tri]) continue;
      cmin2[3 * p2] = cmin[3 * (int64_t)tri];
      cmin2[3 * p2 + 1] = cmin[3 * (int64_t)tri + 1];
      cmin2[3 * p2 + 2] = cmin[3 * (int64_t)tri + 2];
      cmax2[3 * p2] = cmax[3 * (int64_t)tri];
      cmax2[3 * p2 + 1] = cmax[3 * (int64_t)tri + 1];
      cmax2[3 * p2 + 2] = cmax[3 * (int64_t)tri + 2];
      cent2[3 * p2] = cent[3 * (int64_t)tri];
      cent2[3 * p2 + 1] = cent[3 * (int64_t)tri + 1];
      cent2[3 * p2 + 2] = cent[3 * (int64_t)tri + 2];
      area2[p2] = area[i];
      nrm2[3 * p2] = nrm[3 * (int64_t)i];
      nrm2[3 * p2 + 1] = nrm[3 * (int64_t)i + 1];
      nrm2[3 * p2 + 2] = nrm[3 * (int64_t)i + 2];
      kd2[p2] = kd[i];
      ks2[p2] = ks[i];
      if (lep2) lep2[p2] = lep[i];
      mtl2[p2] = mtl[i];
      Ltri[p2] = tri; /* А1604: исходный tri для pcs.tri */
      p2++;
    }
    hz_pyr_free(&py);
    if (hz_pyr_build(&py, nt2, cmin2, cmax2, cent2, mtl2, m.lo, m.hi, cell) != 0) {
      fprintf(stderr, "pgather: пирамида фазы B не построилась\n");
      return 2;
    }
    for (p2 = 0; p2 < nt2; p2++)
      py.pcs[p2].tri = Ltri[p2];
    free(Ltri);
    if (hz_pyr_morton(&py) != 0 || hz_pyr_permute(&py, area2, sizeof *area2) != 0 ||
        hz_pyr_permute(&py, nrm2, 3 * sizeof *nrm2) != 0 ||
        hz_pyr_permute(&py, kd2, sizeof *kd2) != 0 || hz_pyr_permute(&py, ks2, sizeof *ks2) != 0 ||
        (lep2 && hz_pyr_permute(&py, lep2, sizeof *lep2) != 0)) {
      fprintf(stderr, "pgather: Morton фазы B не прошёл\n");
      return 2;
    }
    {
      int32_t nup = 0;
      if (hz_pyr_set_lp(&py, g_lparr, &nup) != 0) {
        fprintf(stderr, "pgather: lp фазы B не построился\n");
        return 2;
      }
    }
    for (p2 = 0; p2 < nt2; p2++)
      py.pcs[p2].e = (float)oldE[py.pcs[p2].tri]; /* тёплый старт */
    soB = so;
    soB.ks = ksf > 0.0 ? ks2 : NULL;
    soB.lep = lep2;
    t0 = now_sec();
    if (hz_sw_run(&py, nt2, area2, nrm2, kd2, &soB, &stB, NULL) != 0) return 2;
    t1 = now_sec();
    tB = t1 - t0;
    for (p2 = 0; p2 < nt2; p2++)
      warmE[p2] = py.pcs[p2].e; /* А1606 */
    for (p2 = 0; p2 < nt2; p2++)
      py.pcs[p2].e = 0.0f;
    if (hz_sw_run(&py, nt2, area2, nrm2, kd2, &soB, &stC, NULL) != 0) return 2;
    tC = now_sec() - t1;
    for (p2 = 0; p2 < nt2; p2++) {
      double d = fabs(warmE[p2] - (double)py.pcs[p2].e);
      dsum += d;
      if (d > dmax) dmax = d;
      esum += (double)py.pcs[p2].e;
    }
    printf("РАЗРУШЕНИЕ: свип тёплый %.3f с, холодный %.3f с; |ΔE| среднее %.4g, max %.4g "
           "(E_cold среднее %.4f); E_avg warm %.4f / cold %.4f\n",
           tB, tC, dsum / nt2, dmax, esum / nt2, stB.e_avg, stC.e_avg);
    for (p2 = 0; p2 < nt2; p2++)
      py.pcs[p2].e = (float)warmE[p2]; /* кадр = тёплый */
    /* камера/сбор живут на фазе B: swap piece-indexed массивов (А1604) */
    free(kd);
    kd = kd2;
    free(nrm);
    nrm = nrm2;
    free(ks);
    ks = ks2;
    if (lep) free(lep);
    lep = lep2;
    free(oldE);
    free(keep);
    free(warmE);
    free(mtl2);
    free(cent2);
    free(cmin2);
    free(cmax2);
    free(area2);
    if (gather == 0) { /* А1605: CSR по новой пирамиде */
      pg_bbox_csr_free(&csr);
      if (pg_bbox_csr_build(&py, cmin, cmax, &csr) != 0) {
        fprintf(stderr, "pgather: CSR фазы B не построился\n");
        return 2;
      }
    }
  }

  /* --- К27: два касающихся шара — лучи через диск касания идут в БЛИЖНИЙ.
   * Шары сцены spheres.obj: центры и радиус восстанавливаются из bbox
   * кусков по знаку координаты x (левый/правый); касание — точка середины
   * между центрами. Контрольные лучи идут ИЗ глаза, стоящего на оси
   * ЛЕВОГО шара дальше от точки касания, в точки диска радиуса
   * 0.9·r вокруг направления на точку касания: обязательное попадание в
   * ближайший шар (t ближнего < t дальнего и сам луч вообще не доходит до
   * дальнего: t_ближ < t_даль − эпсилон геометрии). */
  if (k27 > 0) {
    double lo[3], hi[3];
    int v;
    for (ax = 0; ax < 3; ax++) {
      lo[ax] = hi[ax] = cent[3 * (int64_t)0 + ax];
      /* bbox по ЦЕНТРАМ кусков — куски сгруппированы вокруг двух шаров */
    }
    {
      /* кластеризация по x: порог — середина bbox по x */
      double xmid, csum[2][3] = {{0}};
      int64_t cnt[2] = {0, 0};
      double rad[2] = {0, 0};
      double c1[3], dir[3], tN, tF;
      int nk = 0, through = 0;
      xmid = 0.5 * (cmin[0] + cmax[0]);
      for (i = 0; i < m.nt; i++) {
        int s = cent[3 * (int64_t)i] < xmid ? 0 : 1;
        for (ax = 0; ax < 3; ax++)
          csum[s][ax] += cent[3 * (int64_t)i + ax];
        cnt[s]++;
      }
      for (ax = 0; ax < 3; ax++) {
        c1[ax] = csum[0][ax] / (double)cnt[0];
        lo[ax] = csum[1][ax] / (double)cnt[1];
      }
      for (i = 0; i < m.nt; i++) {
        int s = cent[3 * (int64_t)i] < xmid ? 0 : 1;
        double dd = 0;
        for (ax = 0; ax < 3; ax++) {
          double dxx = cent[3 * (int64_t)i + ax] - (s ? lo[ax] : c1[ax]);
          dd += dxx * dxx;
        }
        if (dd > rad[s]) rad[s] = dd;
      }
      rad[0] = sqrt(rad[0]);
      rad[1] = sqrt(rad[1]);
      printf("К27: центры (%.3f %.3f %.3f) r=%.3f и (%.3f %.3f %.3f) r=%.3f\n", c1[0], c1[1], c1[2],
             rad[0], lo[0], lo[1], lo[2], rad[1]);
      /* глаз — на продолжении линии центров за ЛЕВЫМ шаром */
      for (ax = 0; ax < 3; ax++)
        dir[ax] = lo[ax] - c1[ax];
      {
        double nn = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
        for (ax = 0; ax < 3; ax++)
          dir[ax] /= nn;
      }
      for (ax = 0; ax < 3; ax++)
        eye[ax] = c1[ax] - dir[ax] * (rad[0] * 4.0);
      /* точка касания — середина между центрами; лучи в диск 0.9r вокруг неё */
      for (ax = 0; ax < 3; ax++)
        look[ax] = 0.5 * (c1[ax] + lo[ax]);
      {
        /* базис диска: dir + два перпендикуляра */
        double u[3], w[3], tmp[3] = {0, 0, 1};
        double a0, ud;
        u[0] = dir[1] * tmp[2] - dir[2] * tmp[1];
        u[1] = dir[2] * tmp[0] - dir[0] * tmp[2];
        u[2] = dir[0] * tmp[1] - dir[1] * tmp[0];
        ud = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        for (ax = 0; ax < 3; ax++)
          u[ax] /= ud;
        w[0] = dir[1] * u[2] - dir[2] * u[1];
        w[1] = dir[2] * u[0] - dir[0] * u[2];
        w[2] = dir[0] * u[1] - dir[1] * u[0];
        for (int iy = 0; iy < PG_K27_RAYS; iy++)
          for (int ix = 0; ix < PG_K27_RAYS; ix++) {
            double px = (2.0 * (ix + 0.5) / PG_K27_RAYS - 1.0) * 0.9 * fmin(rad[0], rad[1]);
            double py2 = (2.0 * (iy + 0.5) / PG_K27_RAYS - 1.0) * 0.9 * fmin(rad[0], rad[1]);
            if (px * px + py2 * py2 > 0.81 * fmin(rad[0], rad[1]) * fmin(rad[0], rad[1]))
              continue; /* круг, не квадрат */
            double rd[3];
            for (ax = 0; ax < 3; ax++)
              rd[ax] = look[ax] + px * u[ax] + py2 * w[ax] - eye[ax];
            {
              double nn2 = sqrt(rd[0] * rd[0] + rd[1] * rd[1] + rd[2] * rd[2]);
              for (ax = 0; ax < 3; ax++)
                rd[ax] /= nn2;
            }
            /* ближайшие два пересечения */
            tN = -1.0;
            tF = -1.0;
            for (i = 0; i < m.nt; i++) {
              double p3[3][3], tt;
              hz_obj_tri(&m, py.pcs[i].tri, p3);
              tt = pg_ray_tri(eye, rd, p3);
              if (tt < 0.0) continue;
              if (tN < 0.0 || tt < tN) {
                tF = tN;
                tN = tt;
              } else if (tF < 0.0 || tt < tF) {
                tF = tt;
              }
            }
            nk++;
            if (tN > 0.0 && tF > 0.0 && tN < tF * (1.0 - 1e-12)) continue;
            through++;
          }
        a0 = 100.0 * (double)through / (double)(nk ? nk : 1);
        printf("К27: лучей %d, мимо ближнего %d (%.2f %%) — ОБЯЗАНА быть 0 %%\n", nk, through, a0);
      }
      (void)v;
    }
  }

  /* --- СБОР: перебор кусков, ближайшее t, L_out = le + rho·E/(2π) --- */
  lum = (double *)malloc((size_t)W * (size_t)H * sizeof *lum);
  if (!lum) {
    fprintf(stderr, "pgather: нет памяти на кадр\n");
    return 2;
  }
  for (int fr = 0; fr < frames; fr++) {
    /* §877: ходьба — линейная интерполяция (eye,look) вдоль траектории;
     * поле E решено ОДИН раз выше — кадр стоит только сбор */
    double ef[3], lf[3];
    double tk = frames > 1 ? (double)fr / (double)(frames - 1) : 0.0;
    int64_t nhit = 0;
    double lsum = 0, lmax = 0;
    for (ax = 0; ax < 3; ax++) {
      ef[ax] = eye[ax] + tk * (eye2[ax] - eye[ax]);
      lf[ax] = look[ax] + tk * (look2[ax] - look[ax]);
    }
    if (frames > 1)
      printf("ХОДЬБА: кадр %d/%d eye=(%.2f %.2f %.2f)\n", fr, frames, ef[0], ef[1], ef[2]);
    {
      /* базис камеры */
      double fwd[3], right[3], up[3], tmp[3] = {0, 0, 1};
      double tanf = tan(fov * M_PI / 360.0);
      for (ax = 0; ax < 3; ax++)
        fwd[ax] = lf[ax] - ef[ax];
      {
        double nn = sqrt(fwd[0] * fwd[0] + fwd[1] * fwd[1] + fwd[2] * fwd[2]);
        if (nn < 1e-12) {
          fprintf(stderr, "pgather: глаз совпадает с точкой взгляда\n");
          return 2;
        }
        for (ax = 0; ax < 3; ax++)
          fwd[ax] /= nn;
      }
      right[0] = fwd[1] * tmp[2] - fwd[2] * tmp[1];
      right[1] = fwd[2] * tmp[0] - fwd[0] * tmp[2];
      right[2] = fwd[0] * tmp[1] - fwd[1] * tmp[0];
      {
        double nn = sqrt(right[0] * right[0] + right[1] * right[1] + right[2] * right[2]);
        if (nn < 1e-9) { /* взгляд вдоль z — базис от x */
          tmp[0] = 1;
          tmp[1] = tmp[2] = 0;
          right[0] = fwd[1] * tmp[2] - fwd[2] * tmp[1];
          right[1] = fwd[2] * tmp[0] - fwd[0] * tmp[2];
          right[2] = fwd[0] * tmp[1] - fwd[1] * tmp[0];
          nn = sqrt(right[0] * right[0] + right[1] * right[1] + right[2] * right[2]);
        }
        for (ax = 0; ax < 3; ax++)
          right[ax] /= nn;
      }
      up[0] = fwd[1] * right[2] - fwd[2] * right[1];
      up[1] = fwd[2] * right[0] - fwd[0] * right[2];
      up[2] = fwd[0] * right[1] - fwd[1] * right[0];

      t0 = now_sec();
      pg_cam cam, cam0;
      double *Lumc[3] = {NULL, NULL, NULL};
      if (rgb)
        for (int ch = 0; ch < 3; ch++) {
          Lumc[ch] = (double *)malloc((size_t)W * (size_t)H * sizeof *Lumc[ch]);
          if (!Lumc[ch]) return 2;
        }
      memset(&cam0, 0, sizeof cam0);
      cam0.py = &py;
      cam0.csr = &csr;
      cam0.m = &m;
      cam0.kd = kd;
      cam0.lep = lep;
      cam0.ks = ksf > 0.0 ? ks : NULL; /* §874: ksf=0 — рекурсии не рождаются (А1584) */
      cam0.nrm = nrm;
      cam0.le = le;
      cam0.rho = rho;
      cam0.gather = gather;
      cam0.hop_eps = PG_HOP_EPS_REL * maxdim; /* §874/А1582 */
/* §878: лучи кадра независимы, поле read-only — строка кадра = единица
 * работы; счётчики thread-local с редукцией (А1599/А1600). */
#pragma omp parallel for schedule(dynamic, 16) private(ax, cam)                                    \
    reduction(+ : dda_steps, dda_tested, sec_rays, nhit, lsum) reduction(max : lmax)
      for (int iy = 0; iy < H; iy++) {
        int64_t st_loc = 0, te_loc = 0, sec_loc = 0; /* А1599: thread-local */
        cam = cam0;
        cam.steps = &st_loc;
        cam.tested = &te_loc;
        cam.nsec = &sec_loc;
        for (int ix = 0; ix < W; ix++) {
          double sx = (2.0 * (ix + 0.5) / W - 1.0) * tanf;
          double sy = (1.0 - 2.0 * (iy + 0.5) / H) * tanf * (double)H / (double)W;
          /* явная инициализация rd — анализатор теряет индукцию цикла по ax
           * в OMP-регионе (FP-класс diam, прецедент pg_bbox_span/e1);
           * все элементы перезаписываются — арифметика не меняется */
          double rd[3] = {0, 0, 0}, thit = -1.0, L = 0.0;
          int32_t pbest;
          for (ax = 0; ax < 3; ax++)
            rd[ax] = fwd[ax] + sx * right[ax] + sy * up[ax];
          {
            double nn = sqrt(rd[0] * rd[0] + rd[1] * rd[1] + rd[2] * rd[2]);
            for (ax = 0; ax < 3; ax++)
              rd[ax] /= nn;
          }
          /* §874/А1585: базовый и вторичный лучи — одним сборщиком */
          pbest = pg_nearest(&py, &csr, &m, ef, rd, gather, &thit, &st_loc, &te_loc);
          if (pbest >= 0) {
            if (rgb) { /* §889: поканальная яркость (зеркальный член — v2) */
              for (int ch = 0; ch < 3; ch++) {
                double kdvis_c = rho < 0 ? kdc[ch][pbest] : rho;
                double Lc = le + lepc[ch][pbest] + kdvis_c * Ec[ch][pbest] / (2.0 * M_PI);
                Lumc[ch][(size_t)iy * (size_t)W + (size_t)ix] = Lc;
                if (ch == 1) L = Lc; /* зелёный — яркостная метрика */
              }
            } else {
              L = pg_lcam_hit(&cam, ef, rd, pbest, thit, 0);
            }
            nhit++;
            lsum += L;
            if (L > lmax) lmax = L;
          }
          lum[(size_t)iy * (size_t)W + (size_t)ix] = L;
        }
        dda_steps += st_loc;
        dda_tested += te_loc;
        sec_rays += sec_loc;
      }
      t1 = now_sec();
      if (gather == 0)
        printf("DDA: клеток/луч %.1f, кусков/луч %.1f (%.2f %% от nt)\n",
               (double)dda_steps / ((double)W * (double)H),
               (double)dda_tested / ((double)W * (double)H),
               100.0 * (double)dda_tested / ((double)W * (double)H) / (double)m.nt);
      printf("СБОР: лучей %d, попало %" PRId64
             " (%.2f %%), средняя яркость %.4f, max %.4f, %.2f с\n",
             W * H, nhit, 100.0 * (double)nhit / ((double)W * (double)H),
             lsum / (nhit ? (double)nhit : 1.0), lmax, t1 - t0);
      printf("§874: вторичных зеркальных лучей %" PRId64 " (ksf=%.3g)\n", sec_rays, ksf);
      if (ksf > 0.0) { /* §881: потери hop-политики видны числом (П3) */
        printf("§881: хопов %" PRId64 ", hop_lost=%.4g (absorbed=%.4g)\n", st.hops, st.hop_lost,
               st.absorbed);
        printf("§887-b: row-ветка исполнена %.0f раз\n", st.row_dep);
        printf("§887-d: депозитов Lin>0.5: %lld, Lin<=0.5: %lld\n", (long long)st.lin_pos,
               (long long)st.lin_zero);
        printf("§894: деп_main=%.1f деп_ног=%.1f (absorbed=%.1f)\n", st.dep_main, st.dep_leg,
               st.absorbed);
      }
      {
        double rr = rho < 0 ? 0.5 : rho;
        double Lpred = le + rr * st.e_avg / (2.0 * M_PI);
        double Lmeas = lsum / (nhit ? (double)nhit : 1.0);
        printf("СЛИЧЕНИЕ: L_свипа = le + rho·E_avg/2π = %.4f; L_сбора = %.4f; отношение %.4f\n",
               Lpred, Lmeas, Lmeas / Lpred);
      }
      {
        char fname[4096];
        FILE *f;
        if (fr == 0)
          snprintf(fname, sizeof fname, "%s", outfile);
        else
          snprintf(fname, sizeof fname, "%s_f%02d.ppm", outfile, fr);
        f = fopen(fname, "wb");
        if (!f) {
          fprintf(stderr, "pgather: не открыть %s\n", fname);
          return 2;
        }
        fprintf(f, "P6\n%d %d\n255\n", W, H);
        if (rgb) { /* §889: экспозиция по средней ключевой яркости + гамма;
                    * перцентиль душил тени (динамический диапазон лампы) */
          double *lutmp = (double *)malloc((size_t)W * (size_t)H * sizeof *lutmp);
          double expk = 0.0;
          int64_t npos = 0;
          if (!lutmp) return 2;
          for (i = 0; i < W * H; i++)
            lutmp[i] = Lumc[1][i];
          for (i = 0; i < W * H; i++)
            if (lutmp[i] > 0) {
              expk += lutmp[i];
              npos++;
            }
          expk = expk * expmul / (npos > 0 ? (double)npos : 1.0) + 1e-9;
          free(lutmp);
          for (i = 0; i < W * H; i++)
            for (int ch = 0; ch < 3; ch++) {
              double v = 255.0 * pow(Lumc[ch][i] / expk < 0 ? 0 : Lumc[ch][i] / expk, 1.0 / 2.2);
              if (v > 255.0) v = 255.0;
              unsigned char bb = (unsigned char)v;
              if (fwrite(&bb, 1, 1, f) != 1) return 2;
            }
        } else
          for (i = 0; i < W * H; i++) {
            double v = lum[i] / (lmax > 0 ? lmax : 1.0) * 255.0;
            unsigned char b[3];
            if (v > 255.0) v = 255.0;
            if (v < 0.0) v = 0.0;
            b[0] = b[1] = b[2] = (unsigned char)v;
            fwrite(b, 1, 3, f);
          }
        fclose(f);
        printf("КАДР: %s записан (P6, нормировка на max кадра)\n", fname);
      }
    } /* базис камеры */
  } /* §877: ходьба */

  free(lum);
  free(g_tv9); /* §874: утечка trivert/tribox/lp (предсуществующая с §867,
                * поймана ASAN при прогоне §874) */
  free(g_tb6);
  free(g_lparr);
  if (gather == 0) pg_bbox_csr_free(&csr);
  free(area);
  free(nrm);
  free(kd);
  free(ks);
  free(lep);
  free(cent);
  free(cmin);
  free(cmax);
  free(mtl);
  hz_pyr_free(&py);
  hz_obj_free(&m);
  return 0;
}
