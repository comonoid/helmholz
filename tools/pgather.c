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
#include "geom/kit.h"
#include "scene_obj.h"
#ifdef _OPENMP /* §921: сериализация прибора-приёмника */
#include <omp.h>
#endif

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
  /* §930-Г1: мини-группы (4 подряд Мортона) с предвычисленным bbox —
   * кул одним слэб-тестом до ray-tri кусков клетки; консервативно
   * (bbox группы ⊇ bbox кусков), сбор побитово тот же */
  int64_t *gstart; /* [nleaf+1] индексы групп */
  double *gbox;    /* [6*ngroups] */
  int64_t ngroups;
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
  free(csr->gstart); /* §930-Г1 */
  free(csr->gbox);
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
  { /* §930-Г1: группы по 4 подряд (слоты мортона-упорядочены) */
    int64_t ng = 0;
    for (li = 0; li < py->nleaf; li++)
      ng += ((csr->start[li + 1] - csr->start[li]) + 3) >> 2;
    csr->gstart = (int64_t *)calloc((size_t)py->nleaf + 1, sizeof *csr->gstart);
    csr->gbox = (double *)malloc((size_t)(ng > 0 ? ng : 1) * 6 * sizeof *csr->gbox);
    if (!csr->gstart || !csr->gbox) {
      free(csr->gstart);
      free(csr->gbox);
      csr->gstart = NULL;
      csr->gbox = NULL;
      return 2; /* fail closed: без групп прежний путь не собирается */
    }
    csr->ngroups = ng;
    {
      int64_t g = 0;
      for (li = 0; li < py->nleaf; li++) {
        csr->gstart[li] = g;
        for (int64_t s = csr->start[li]; s < csr->start[li + 1]; s += 4) {
          int64_t e = s + 4 < csr->start[li + 1] ? s + 4 : csr->start[li + 1];
          double mn[3] = {0, 0, 0}, mx[3] = {0, 0, 0};
          for (int64_t q = s; q < e; q++) {
            int32_t tri = py->pcs[csr->pids[q]].tri;
            for (int a = 0; a < 3; a++) {
              double lo2 = cmin[3 * (int64_t)tri + a], hi2 = cmax[3 * (int64_t)tri + a];
              if (q == s || lo2 < mn[a]) mn[a] = lo2;
              if (q == s || hi2 > mx[a]) mx[a] = hi2;
            }
          }
          for (int a = 0; a < 3; a++) {
            csr->gbox[6 * g + a] = mn[a];
            csr->gbox[6 * g + 3 + a] = mx[a];
          }
          g++;
        }
      }
      csr->gstart[py->nleaf] = g;
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
      { /* §930-М3б: кусков в посещённой клетке (развязка Г1) */
        extern int64_t g_m3b_inlist, g_m3b_visits;
        g_m3b_inlist += csr->start[pos + 1] - csr->start[pos];
        g_m3b_visits++;
      }
      if (csr->gbox != NULL) { /* §930-Г1: кул групп (4 подряд) слэбом */
        extern int64_t g_m3c_gtests, g_m3c_ghit;
        for (int64_t g = csr->gstart[pos]; g < csr->gstart[pos + 1]; g++) {
          const double *gb = csr->gbox + 6 * g;
          double t0g = -1e30, t1g = 1e30;
          int miss = 0;
          for (int a = 0; a < 3; a++) {
            double inv = 1.0 / rd[a];
            double ta = (gb[a] - eye[a]) * inv, tb = (gb[3 + a] - eye[a]) * inv;
            double tlo2 = ta < tb ? ta : tb, thi2 = ta > tb ? ta : tb;
            if (tlo2 > t0g) t0g = tlo2;
            if (thi2 < t1g) t1g = thi2;
            if (t0g > t1g) {
              miss = 1;
              break;
            }
          }
          g_m3c_gtests++;
          if (miss) continue; /* группа мимо — её 4 куска не тестируются */
          g_m3c_ghit++;
          int64_t s0 = csr->start[pos] + ((g - csr->gstart[pos]) << 2);
          int64_t s1e = s0 + 4 < csr->start[pos + 1] ? s0 + 4 : csr->start[pos + 1];
          for (s = s0; s < s1e; s++) {
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
      } else
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
  const double *ev; /* §932-Б-Ш6: пер-вершинное E (площадь-взвешенное); NULL —
                       плоско-кусочный P0 прежнего мира, битово */
  double le, rho;
  int gather;
  double hop_eps;
  int64_t *steps, *tested;
  int64_t *nsec; /* зеркальных вторичных лучей (прибор) */
} pg_cam;

/* §932-Б-Ш6: E в точке хита — барицентрическая интерполяция пер-вершинных
 * значений (egour=1); NULL — плоско-кусочный P0, побитово прежний мир.
 * λ считаются через двойные скалярные произведения (без делений на
 * ребро); отрицательные хвосты числителя клэмпятся к 0 и сумма λ
 * нормируется — на ребре/вершине интерполяция вырождается корректно. */
static double pg_e_hit(const pg_cam *c, const double org[3], const double rd[3], int32_t p,
                       double thit) {
  int32_t t = c->py->pcs[p].tri;
  const double *A = c->m->v + 3 * (int64_t)c->m->f[3 * (int64_t)t];
  const double *B = c->m->v + 3 * (int64_t)c->m->f[3 * (int64_t)t + 1];
  const double *C = c->m->v + 3 * (int64_t)c->m->f[3 * (int64_t)t + 2];
  double e1[3], e2[3], ep[3], P[3];
  double d11, d12, d22, d1p, d2p, den, u, v, w;
  int ax;
  for (ax = 0; ax < 3; ax++)
    P[ax] = org[ax] + rd[ax] * thit;
  for (ax = 0; ax < 3; ax++) {
    e1[ax] = B[ax] - A[ax];
    e2[ax] = C[ax] - A[ax];
    ep[ax] = P[ax] - A[ax];
  }
  d11 = e1[0] * e1[0] + e1[1] * e1[1] + e1[2] * e1[2];
  d12 = e1[0] * e2[0] + e1[1] * e2[1] + e1[2] * e2[2];
  d22 = e2[0] * e2[0] + e2[1] * e2[1] + e2[2] * e2[2];
  d1p = e1[0] * ep[0] + e1[1] * ep[1] + e1[2] * ep[2];
  d2p = e2[0] * ep[0] + e2[1] * ep[1] + e2[2] * ep[2];
  den = d11 * d22 - d12 * d12;
  if (!(den > 1e-30)) return (double)c->py->pcs[p].e; /* вырожденный tri — P0 */
  u = (d22 * d1p - d12 * d2p) / den;
  v = (d11 * d2p - d12 * d1p) / den;
  if (u < 0.0) u = 0.0;
  if (v < 0.0) v = 0.0;
  w = 1.0 - u - v;
  if (w < 0.0) {
    w = 0.0;
    if (u + v > 0.0) {
      u /= u + v;
      v = 1.0 - u;
    }
  }
  return u * c->ev[c->m->f[3 * (int64_t)t + 1]] + v * c->ev[c->m->f[3 * (int64_t)t + 2]] +
         w * c->ev[c->m->f[3 * (int64_t)t]];
}

static double pg_lcam_hit(const pg_cam *c, const double org[3], const double rd[3], int32_t p,
                          double thit, int depth) {
  double kdvis = c->rho < 0 ? c->kd[p] : c->rho;
  double E = c->ev ? pg_e_hit(c, org, rd, p, thit) : (double)c->py->pcs[p].e;
  double L = c->le + (c->lep ? c->lep[p] : 0.0) + kdvis * E / (2.0 * M_PI);
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

/* --- §921: ИНВАРИАНТНЫЙ ПРИЁМНИК (метрология носителей §920-И) -----------
 * Приёмник = НЕизменный L0-меш прогона: лучи и точки попаданий одни и те
 * же для всех носителей; поле на точке попадания берётся с треугольника
 * НОСИТЕЛЯ (уровень Lk кита-лестницы, поле = sidecar §898 per-tri).
 * Знаменатель метрик — число попаданий: от носителя НЕ зависит (диагноз
 * §920-И: у E_avg знаменатель = Σ площадей кусков носителя).
 * Lookup: сетка центроидов носителя (шаг 2·√⟨S⟩), point-in-tri в плоскости
 * с запасом (носитель в ε от приёмника), fallback — ближайший центроид
 * (считается; на носителе L0 тождество с обычным сбором БИТОВОЕ). */
#define PG_RECV_CELL_CAP 512 /* кэп клеток на треугольник; крупнее — гиганты */
#define PG_RECV_RMAX 8       /* колец поиска: 8·h покрывает ε уровня и стыки */
/* строгая принадлежность: тай-полоса рёбер разбирала бы носитель L0 по
 * соседям и ломала битовое тождество; зазоры носителя закрывает
 * fallback-центроид */
#define PG_RECV_BTOL 1e-9

typedef struct {
  hz_kit kit;            /* лестница носителей; живёт до конца сбора */
  const hz_kit_level *S; /* уровень-носитель kit.lev[lv] */
  int32_t lv;
  double *E;      /* [S->ntris] поле носителя (sidecar §898) */
  int64_t *hkey;  /* хеш клеток: код клетки или -1 */
  int32_t *hhead; /* голова цепи треугольников клетки или -1 */
  int32_t *tnext; /* [nreg] следующая регистрация в цепи клетки */
  int32_t *rtri;  /* [nreg] треугольник регистрации */
  int32_t *giant; /* [ngiant] треугольники крупнее кэпа клеток */
  uint32_t ngiant;
  uint32_t hmask;
  double gorg[3], ginv; /* сетка центроидов: начало и 1/h */
  double h, ptol;
  int64_t nlook, nfall, nempty; /* приборы lookup (НК покрытия) */
} pg_recv;

static uint32_t pg_recv_hash(int64_t x, int64_t y, int64_t z) {
  uint64_t hsh = (uint64_t)x * 0x9E3779B97F4A7C15ull ^ (uint64_t)y * 0xC2B2AE3D27D4EB4Full ^
                 (uint64_t)z * 0x165667B19E3779F9ull;
  hsh ^= hsh >> 32;
  return (uint32_t)hsh;
}

/* Открыть носитель: кит + уровень lv + поле из sidecar (ровно ntris double).
 * scramble — НК: Фишер–Йетс с фикс. LCG (seed 1, минимальный стандартный),
 * применяемый ТОЛЬКО к не-эталонным носителям вызывающим. 0 — успех. */
static int pg_recv_open(pg_recv *r, const char *kitpath, int32_t lv, const char *efile,
                        int scramble) {
  memset(r, 0, sizeof *r);
  hz_kit_init(&r->kit);
  FILE *f = fopen(kitpath, "rb");
  int rc = f != NULL ? hz_kit_load(&r->kit, f) : HZ_KIT_E_IO;
  if (f != NULL) fclose(f);
  if (rc != HZ_KIT_OK) return 2;
  if (lv < 0 || lv >= r->kit.nlev) return 2;
  r->S = &r->kit.lev[lv];
  r->lv = lv;
  uint32_t nt = r->S->ntris;
  r->E = (double *)malloc((size_t)nt * sizeof *r->E);
  if (r->E == NULL) return 2;
  FILE *fe = fopen(efile, "rb");
  if (fe == NULL || fread(r->E, sizeof(double), (size_t)nt, fe) != (size_t)nt) {
    if (fe != NULL) fclose(fe);
    return 2;
  }
  fclose(fe);
  if (scramble) { /* НК: поле перемешано по треугольникам носителя */
    uint32_t st = 1;
    for (uint32_t j = nt - 1; j > 0; j--) {
      st = (uint32_t)((uint64_t)st * 48271ull % 2147483647ull);
      uint32_t kk = st % (j + 1);
      double tmp = r->E[j];
      r->E[j] = r->E[kk];
      r->E[kk] = tmp;
    }
  }
  /* bbox уровня + средняя площадь → шаг сетки центроидов */
  double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0}, meanA = 0.0;
  for (uint32_t v = 0; v < r->S->nverts; v++)
    for (int a = 0; a < 3; a++) {
      double c = a == 0 ? r->S->vx[v] : (a == 1 ? r->S->vy[v] : r->S->vz[v]);
      if (v == 0 || c < lo[a]) lo[a] = c;
      if (v == 0 || c > hi[a]) hi[a] = c;
    }
  for (uint32_t t = 0; t < nt; t++) {
    uint32_t i0 = r->S->ti0[t], i1 = r->S->ti1[t], i2 = r->S->ti2[t];
    double ax = r->S->vx[i1] - r->S->vx[i0], ay = r->S->vy[i1] - r->S->vy[i0],
           az = r->S->vz[i1] - r->S->vz[i0];
    double bx = r->S->vx[i2] - r->S->vx[i0], by = r->S->vy[i2] - r->S->vy[i0],
           bz = r->S->vz[i2] - r->S->vz[i0];
    meanA += 0.5 * sqrt((ay * bz - az * by) * (ay * bz - az * by) +
                        (az * bx - ax * bz) * (az * bx - ax * bz) +
                        (ax * by - ay * bx) * (ax * by - ay * bx));
  }
  meanA = meanA / (double)nt;
  r->h = 2.0 * sqrt(meanA > 1e-18 ? meanA : 1e-18);
  r->ginv = 1.0 / r->h;
  r->ptol = r->h; /* запас плоскости: диагональ клетки с покрытием ε */
  for (int a = 0; a < 3; a++)
    r->gorg[a] = lo[a] - r->h;
  uint32_t hs = 16;
  while (hs < 2u * nt)
    hs <<= 1;
  r->hmask = hs - 1u;
  r->hkey = (int64_t *)malloc((size_t)hs * sizeof *r->hkey);
  r->hhead = (int32_t *)malloc((size_t)hs * sizeof *r->hhead);
  if (!r->hkey || !r->hhead) return 2;
  for (uint32_t u = 0; u < hs; u++) {
    r->hkey[u] = -1;
    r->hhead[u] = -1;
  }
  r->giant = (int32_t *)malloc((size_t)nt * sizeof *r->giant); /* верхняя оценка */
  { /* пре-пас: точный бюджет регистраций (кэп на треугольник) */
    uint64_t tot = 0;
    for (uint32_t t = 0; t < nt; t++) {
      uint32_t i0 = r->S->ti0[t], i1 = r->S->ti1[t], i2 = r->S->ti2[t];
      double mn[3], mx[3], a0 = r->S->vx[i0], a1 = r->S->vx[i1], a2 = r->S->vx[i2];
      double b0 = r->S->vy[i0], b1 = r->S->vy[i1], b2 = r->S->vy[i2];
      double c0d = r->S->vz[i0], c1d = r->S->vz[i1], c2d = r->S->vz[i2];
      mn[0] = a0 < a1 ? (a0 < a2 ? a0 : a2) : (a1 < a2 ? a1 : a2);
      mx[0] = a0 > a1 ? (a0 > a2 ? a0 : a2) : (a1 > a2 ? a1 : a2);
      mn[1] = b0 < b1 ? (b0 < b2 ? b0 : b2) : (b1 < b2 ? b1 : b2);
      mx[1] = b0 > b1 ? (b0 > b2 ? b0 : b2) : (b1 > b2 ? b1 : b2);
      mn[2] = c0d < c1d ? (c0d < c2d ? c0d : c2d) : (c1d < c2d ? c1d : c2d);
      mx[2] = c0d > c1d ? (c0d > c2d ? c0d : c2d) : (c1d > c2d ? c1d : c2d);
      uint64_t span = 1;
      for (int a = 0; a < 3; a++) {
        int64_t lo_c = (int64_t)floor((mn[a] - r->gorg[a]) * r->ginv);
        int64_t hi_c = (int64_t)floor((mx[a] - r->gorg[a]) * r->ginv);
        if (hi_c < lo_c) hi_c = lo_c;
        span *= (uint64_t)(hi_c - lo_c + 1);
        if (span > (uint64_t)PG_RECV_CELL_CAP) break;
      }
      if (span <= (uint64_t)PG_RECV_CELL_CAP) tot += span;
    }
    r->tnext = (int32_t *)malloc((size_t)tot * sizeof *r->tnext);
    r->rtri = (int32_t *)malloc((size_t)tot * sizeof *r->rtri);
  }
  if (!r->giant || !r->tnext || !r->rtri) return 2;
  uint32_t nreg = 0; /* суммарные регистрации (кэп на треугольник) */
  for (uint32_t t = 0; t < nt; t++) {
    uint32_t i0 = r->S->ti0[t], i1 = r->S->ti1[t], i2 = r->S->ti2[t];
    double tx0 = r->S->vx[i0], tx1 = r->S->vx[i1], tx2 = r->S->vx[i2];
    double ty0 = r->S->vy[i0], ty1 = r->S->vy[i1], ty2 = r->S->vy[i2];
    double tz0 = r->S->vz[i0], tz1 = r->S->vz[i1], tz2 = r->S->vz[i2];
    double tmn[3], tmx[3];
    tmn[0] = tx0 < tx1 ? (tx0 < tx2 ? tx0 : tx2) : (tx1 < tx2 ? tx1 : tx2);
    tmx[0] = tx0 > tx1 ? (tx0 > tx2 ? tx0 : tx2) : (tx1 > tx2 ? tx1 : tx2);
    tmn[1] = ty0 < ty1 ? (ty0 < ty2 ? ty0 : ty2) : (ty1 < ty2 ? ty1 : ty2);
    tmx[1] = ty0 > ty1 ? (ty0 > ty2 ? ty0 : ty2) : (ty1 > ty2 ? ty1 : ty2);
    tmn[2] = tz0 < tz1 ? (tz0 < tz2 ? tz0 : tz2) : (tz1 < tz2 ? tz1 : tz2);
    tmx[2] = tz0 > tz1 ? (tz0 > tz2 ? tz0 : tz2) : (tz1 > tz2 ? tz1 : tz2);
    int64_t c0[3], c1[3];
    uint32_t span = 1;
    for (int a = 0; a < 3; a++) {
      c0[a] = (int64_t)floor((tmn[a] - r->gorg[a]) * r->ginv);
      c1[a] = (int64_t)floor((tmx[a] - r->gorg[a]) * r->ginv);
      span *= (uint32_t)(c1[a] - c0[a] + 1);
      if (c1[a] < c0[a]) c1[a] = c0[a]; /* robustness: пустой диапазон */
    }
    if (span > PG_RECV_CELL_CAP) { /* ПОЛ/кровля: линейный список гигантов */
      r->giant[r->ngiant++] = (int32_t)t;
      continue;
    }
    for (int64_t cz2 = c0[2]; cz2 <= c1[2]; cz2++)
      for (int64_t cy2 = c0[1]; cy2 <= c1[1]; cy2++)
        for (int64_t cx2 = c0[0]; cx2 <= c1[0]; cx2++) {
          int64_t key = cx2 * 73856093LL + cy2 * 19349663LL + cz2 * 83492791LL;
          uint32_t su = pg_recv_hash(cx2, cy2, cz2) & r->hmask;
          while (r->hkey[su] != -1 && r->hkey[su] != key)
            su = (su + 1u) & r->hmask;
          if (r->hkey[su] == -1) {
            r->hkey[su] = key;
            r->hhead[su] = -1;
          }
          r->tnext[nreg] = r->hhead[su];
          r->rtri[nreg] = (int32_t)t;
          r->hhead[su] = (int32_t)nreg;
          nreg++;
        }
  }
  (void)nreg; /* регистрации все разложены; цепи замкнуты через tnext */
  return 0;
}

static void pg_recv_close(pg_recv *r) {
  free(r->E);
  free(r->hkey);
  free(r->hhead);
  free(r->tnext);
  free(r->rtri);
  free(r->giant);
  hz_kit_free(&r->kit);
}

/* Поле носителя в точке q: PIP в кольцах сетки, fallback ближайший
 * центроид. Детерминировано (порядок цепей фиксирован построением). */
static double pg_recv_E(pg_recv *r, const double q[3]) {
  r->nlook++;
  int64_t c0[3];
  for (int a = 0; a < 3; a++)
    c0[a] = (int64_t)floor((q[a] - r->gorg[a]) * r->ginv);
  const hz_kit_level *S = r->S;
  double bestE = 0.0, bestD = HUGE_VAL, bestM = -HUGE_VAL, bd2 = HUGE_VAL;
  int32_t bestC = -1;
  int found = 0;
  for (int rad = 0; rad <= PG_RECV_RMAX; rad++) {
    for (int dz = -rad; dz <= rad; dz++)
      for (int dy = -rad; dy <= rad; dy++)
        for (int dx = -rad; dx <= rad; dx++) {
          int mx = dx < 0 ? -dx : dx, my = dy < 0 ? -dy : dy, mz = dz < 0 ? -dz : dz;
          int mr = mx > my ? (mx > mz ? mx : mz) : (my > mz ? my : mz);
          if (mr != rad) continue; /* только оболочка кольца rad */
          int64_t ix = c0[0] + dx, iy = c0[1] + dy, iz = c0[2] + dz;
          int64_t key = ix * 73856093LL + iy * 19349663LL + iz * 83492791LL;
          uint32_t s = pg_recv_hash(ix, iy, iz) & r->hmask;
          while (r->hkey[s] != -1) {
            if (r->hkey[s] == key) {
              for (int32_t sl = r->hhead[s]; sl >= 0; sl = r->tnext[sl]) {
                int32_t t = r->rtri[sl];
                uint32_t i0 = S->ti0[t], i1 = S->ti1[t], i2 = S->ti2[t];
                double x0 = S->vx[i0], y0 = S->vy[i0], z0 = S->vz[i0];
                double e1x = S->vx[i1] - x0, e1y = S->vy[i1] - y0, e1z = S->vz[i1] - z0;
                double e2x = S->vx[i2] - x0, e2y = S->vy[i2] - y0, e2z = S->vz[i2] - z0;
                double wx = q[0] - x0, wy = q[1] - y0, wz = q[2] - z0;
                double nx = e1y * e2z - e1z * e2y, ny = e1z * e2x - e1x * e2z,
                       nz = e1x * e2y - e1y * e2x;
                double n2 = nx * nx + ny * ny + nz * nz;
                if (n2 > 1e-30) { /* вырожденные — только в fallback-центроид */
                  /* барицентры проекции: u при v1 (e2×n), v при v2 (n×e1) */
                  double u = (wx * (e2y * nz - e2z * ny) + wy * (e2z * nx - e2x * nz) +
                              wz * (e2x * ny - e2y * nx)) /
                             n2;
                  double v = (wx * (ny * e1z - nz * e1y) + wy * (nz * e1x - nx * e1z) +
                              wz * (nx * e1y - ny * e1x)) /
                             n2;
                  double dn = wx * nx + wy * ny + wz * nz;
                  double d2 = dn * dn / n2;
                  if (u >= -PG_RECV_BTOL && v >= -PG_RECV_BTOL && u + v <= 1.0 + PG_RECV_BTOL &&
                      d2 <= r->ptol * r->ptol) {
                    /* тай-брейк копланарных (общее ребро): побеждает самый
                     * «внутренний» — на носителе L0 это треугольник попадания,
                     * тождество с обычным сбором почти битово */
                    double marg = u < v ? u : v;
                    if (1.0 - u - v < marg) marg = 1.0 - u - v;
                    double dtol = r->ptol * r->ptol * 1e-9;
                    if (d2 < bestD - dtol || (d2 <= bestD + dtol && marg > bestM)) {
                      bestD = d2;
                      bestM = marg;
                      bestE = r->E[t];
                      found = 1;
                    }
                  }
                }
                double gx = (S->vx[i0] + S->vx[i1] + S->vx[i2]) / 3.0;
                double gy = (S->vy[i0] + S->vy[i1] + S->vy[i2]) / 3.0;
                double gz = (S->vz[i0] + S->vz[i1] + S->vz[i2]) / 3.0;
                double ddx = gx - q[0], ddy = gy - q[1], ddz = gz - q[2];
                double dd2 = ddx * ddx + ddy * ddy + ddz * ddz;
                if (dd2 < bd2) {
                  bd2 = dd2;
                  bestC = t;
                }
              }
            }
            s = (s + 1u) & r->hmask;
          }
        }
    if (found && rad >= 3) break; /* содержащий дальше 3 колец хуже найденного */
  }
  { /* гиганты (ПОЛ/кровля): линейный проход, те же критерии */
    for (uint32_t gi = 0; gi < r->ngiant; gi++) {
      int32_t t = r->giant[gi];
      uint32_t i0 = S->ti0[t], i1 = S->ti1[t], i2 = S->ti2[t];
      double x0 = S->vx[i0], y0 = S->vy[i0], z0 = S->vz[i0];
      double e1x = S->vx[i1] - x0, e1y = S->vy[i1] - y0, e1z = S->vz[i1] - z0;
      double e2x = S->vx[i2] - x0, e2y = S->vy[i2] - y0, e2z = S->vz[i2] - z0;
      double wx = q[0] - x0, wy = q[1] - y0, wz = q[2] - z0;
      double nx = e1y * e2z - e1z * e2y, ny = e1z * e2x - e1x * e2z, nz = e1x * e2y - e1y * e2x;
      double n2 = nx * nx + ny * ny + nz * nz;
      if (n2 <= 1e-30) continue;
      double u =
          (wx * (e2y * nz - e2z * ny) + wy * (e2z * nx - e2x * nz) + wz * (e2x * ny - e2y * nx)) /
          n2;
      double v =
          (wx * (ny * e1z - nz * e1y) + wy * (nz * e1x - nx * e1z) + wz * (nx * e1y - ny * e1x)) /
          n2;
      double dn = wx * nx + wy * ny + wz * nz;
      double d2 = dn * dn / n2;
      if (u >= -PG_RECV_BTOL && v >= -PG_RECV_BTOL && u + v <= 1.0 + PG_RECV_BTOL &&
          d2 <= r->ptol * r->ptol) {
        double marg = u < v ? u : v;
        if (1.0 - u - v < marg) marg = 1.0 - u - v;
        double dtol = r->ptol * r->ptol * 1e-9;
        if (d2 < bestD - dtol || (d2 <= bestD + dtol && marg > bestM)) {
          bestD = d2;
          bestM = marg;
          bestE = r->E[t];
          found = 1;
        }
      }
      double gx = (S->vx[i0] + S->vx[i1] + S->vx[i2]) / 3.0;
      double gy = (S->vy[i0] + S->vy[i1] + S->vy[i2]) / 3.0;
      double gz = (S->vz[i0] + S->vz[i1] + S->vz[i2]) / 3.0;
      double ddx = gx - q[0], ddy = gy - q[1], ddz = gz - q[2];
      double dd2 = ddx * ddx + ddy * ddy + ddz * ddz;
      if (dd2 < bd2) {
        bd2 = dd2;
        bestC = t;
      }
    }
  }
  if (found) return bestE;
  if (bestC >= 0) {
    r->nfall++;
    return r->E[bestC];
  }
  r->nempty++; /* носитель не покрыл точку вовсе (контроль §914) */
  return 0.0;
}

static int pg_recv_dcmp(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return x < y ? -1 : (x > y ? 1 : 0);
}

/* --- §882: HBLK v1 — mmap-сбор (блок = листовая клетка, Morton);
 * §902: v2 — kd/lep/E в файле, освещённый кадр без OBJ;
 * §904: v3 — kd3/lep3 в файле, RGB-кадр без OBJ */
#define PG_HBLK_MAGIC 0x314B4C4248ULL  /* "HBLK1" */
#define PG_HBLK_MAGIC2 0x324B4C4248ULL /* "HBLK2" */
#define PG_HBLK_MAGIC3 0x334B4C4248ULL /* "HBLK3" */
typedef struct {
  uint8_t *base; /* начало mmap */
  size_t len;
  const uint64_t *hdr;
  const uint8_t *tab; /* pb_cellrec {cell,start,cnt+pad} 24 Б */
  const int32_t *ids;
  const double *tris; /* 9 double на tri, исходный порядок */
  const double *kd;   /* §902: v2/v3 — отклик/поле per tri (NULL в v1) */
  const double *lep;
  const double *E;
  const double *kd3; /* §904: v3 — каналы для RGB-кадра из файла */
  const double *lep3;
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
  if (len < (off_t)24 * 8) { /* приведение к типу сравнения: 24·8 — литерал заголовка HBLK */
    close(fd);
    return 2;
  }
  h->len = (size_t)len;
  h->base = (uint8_t *)mmap(NULL, h->len, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  if (h->base == MAP_FAILED) return 2;
  madvise(h->base, h->len, MADV_RANDOM);
  h->hdr = (const uint64_t *)h->base;
  if (h->hdr[0] != PG_HBLK_MAGIC && h->hdr[0] != PG_HBLK_MAGIC2 && h->hdr[0] != PG_HBLK_MAGIC3)
    return 2;
  h->kd = h->lep = h->E = h->kd3 = h->lep3 = NULL;
  if (h->hdr[0] == PG_HBLK_MAGIC2 || h->hdr[0] == PG_HBLK_MAGIC3) { /* §902: отклик/поле */
    h->kd = (const double *)(h->base + h->hdr[13]);
    h->lep = (const double *)(h->base + h->hdr[14]);
    h->E = (const double *)(h->base + h->hdr[15]);
  }
  if (h->hdr[0] == PG_HBLK_MAGIC3) { /* §904: каналы */
    h->kd3 = (const double *)(h->base + h->hdr[16]);
    h->lep3 = (const double *)(h->base + h->hdr[17]);
  }
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

/* --- §896: клип треугольника клеткой (Sutherland–Hodgman по 6 плоскостям).
 * Возврат: площадь клипа; poly/np — вершины. */
static double sh_clip_tri(const double tri[3][3], const double blo[3], const double bhi[3],
                          double out[3][3]) {
  double poly[16][3], tmp[16][3];
  int np = 3, nt2, pl; /* nt2 обнуляется в первом же витке (pl=0) до чтения:
                        * инициализатор здесь был мёртвым присваиванием (07-10) */
  for (int i = 0; i < 3; i++)
    for (int a = 0; a < 3; a++)
      poly[i][a] = tri[i][a];
  const int axv[6] = {0, 0, 1, 1, 2, 2};
  const int sgv[6] = {-1, 1, -1, 1, -1, 1};
  const double *bnd[6] = {blo, bhi, blo, bhi, blo, bhi};
  for (pl = 0; pl < 6 && np > 0; pl++) {
    int a = axv[pl], sg = sgv[pl];
    double bndv = bnd[pl][a];
    nt2 = 0;
    for (int i = 0; i < np; i++) {
      const double *cur = poly[i], *nxt = poly[(i + 1) % np];
      double c = cur[a], n = nxt[a];
      int cin = (sg < 0) ? (c >= bndv - 1e-15) : (c <= bndv + 1e-15);
      int nin = (sg < 0) ? (n >= bndv - 1e-15) : (n <= bndv + 1e-15);
      if (cin) {
        for (int q = 0; q < 3; q++)
          tmp[nt2][q] = cur[q];
        nt2++;
      }
      if ((cin && !nin) || (!cin && nin)) {
        double t = (bndv - c) / (n - c);
        for (int q = 0; q < 3; q++)
          tmp[nt2][q] = cur[q] + t * (nxt[q] - cur[q]);
        nt2++;
      }
    }
    for (int i = 0; i < nt2; i++)
      for (int q = 0; q < 3; q++)
        poly[i][q] = tmp[i][q];
    np = nt2;
  }
  if (np < 3) return 0.0;
  { /* площадь + центроид (через триангуляцию веером) */
    double area = 0.0, cx = 0, cy = 0, cz = 0;
    for (int i = 1; i < np - 1; i++) {
      double e1[3], e2[3], cr[3];
      for (int a = 0; a < 3; a++) {
        e1[a] = poly[i][a] - poly[0][a];
        e2[a] = poly[i + 1][a] - poly[0][a];
      }
      cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
      cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
      cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
      double aa = 0.5 * sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
      area += aa;
      cx += aa * (poly[0][0] + poly[i][0] + poly[i + 1][0]) / 3.0;
      cy += aa * (poly[0][1] + poly[i][1] + poly[i + 1][1]) / 3.0;
      cz += aa * (poly[0][2] + poly[i][2] + poly[i + 1][2]) / 3.0;
    }
    if (area > 1e-14) {
      out[0][0] = cx / area;
      out[0][1] = cy / area;
      out[0][2] = cz / area;
    }
    return area;
  }
}

/* ---- §924: представители декимации на динамической лестнице ----
 * Дети(R) = исходные треугольники, чей центроид лежит в R: PIP в плоскости
 * с запасом + fallback ближайший центроид (как lookup приёмника §921).
 * Всё в ПРОСТРАНСТВЕ ИСХОДНЫХ ТРЕУГОЛЬНИКОВ; слот-перестановка repof/kmem
 * делается ПОСЛЕ Мортона отдельно. Приборы О1: nfall, maxkids. */
typedef struct {
  int32_t nrep;               /* Σ треугольников уровней L1.. */
  int32_t nlev_r;             /* kk.nlev-1 */
  int32_t *repof;             /* [nlev_r*nt]: offset_f + локальный id; -1 нет */
  int32_t *koff;              /* [nrep+1] CSR детей */
  int32_t *kmem;              /* [Σдетей] tri-индексы (до Мортона) */
  double *rarea, *rrho, *rle; /* [nrep] агрегаты по детям */
  double *rle_area_skip;      /* [nrep] Σ площадей эмиттеров (§924) */
  double *ratri;              /* [nrep] площадь треугольника-носителя (§927) */
  double *rnrm;               /* [3nrep] единичная нормаль носителя (§927) */
  int64_t nfall;              /* О1: fallback-назначения */
  int32_t maxkids;            /* О1 */
} pg_reps;

#define PG924_RMAX 12 /* колец поиска; крупнее §921: центроид может уйти дальше */

typedef struct {
  int32_t *dhead, *tnext; /* §929-Х1b: ПРЯМАЯ сетка head[ncell] + цепочки;
                           * порядок цепочек = вставка по возрастанию t
                           * (тот же, что у хэш-версии), порядок обхода
                           * ячеек — таблицей оболочек: побитово то же */
  int64_t gnx, gny, gnz;  /* размеры сетки (span) */
  int64_t ncell;
  int32_t nreg;
  double gorg[3], ginv;
  double *rmax; /* §929-Х3: max |вершина − центроид| на треугольник —
                 * предфильтр ПИП (консервативный, побитово безопасный) */
  double gmaxr; /* §930-А: max rmax уровня — граница «ПИП мёртв» */
} pg924_grid;

/* §929-Х1: препcomputed-обход shells rad=0..PG924_RMAX в ТОМ же порядке
 * (dx,dy,dz по возрастанию Chebyshev-кольца), что и тройной цикл —
 * убирает ~12× мусорных итераций фильтра mr!=rad. */
int64_t g_m3b_inlist, g_m3b_visits;     /* §930-М3б */
int64_t g_m3c_gtests, g_m3c_ghit;       /* §930-Г1-дых */
static double *g_egour_ev, *g_egour_ew; /* §932-Б-Ш6: буферы Gouraud (владелец — кадр) */

/* pg924_hash удалён 07-10: не вызывался нигде (gcc -Wunused-function) —
 * вернуть из git-истории, если понадобится хеш клетки */

/* сетка по центроидам треугольников уровня; 0/память */
static int pg924_grid_build(pg924_grid *g, const hz_kit_level *S) {
  memset(g, 0, sizeof *g);
  double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
  double mean = 0.0;
  for (uint32_t t = 0; t < S->ntris; t++) {
    double cx = (S->vx[S->ti0[t]] + S->vx[S->ti1[t]] + S->vx[S->ti2[t]]) / 3.0;
    double cy = (S->vy[S->ti0[t]] + S->vy[S->ti1[t]] + S->vy[S->ti1[t]] + S->vy[S->ti2[t]]) /
                4.0; /* не используется */
    (void)cy;
    double cz = (S->vz[S->ti0[t]] + S->vz[S->ti1[t]] + S->vz[S->ti2[t]]) / 3.0;
    double e0x = S->vx[S->ti1[t]] - S->vx[S->ti0[t]], e0y = S->vy[S->ti1[t]] - S->vy[S->ti0[t]],
           e0z = S->vz[S->ti1[t]] - S->vz[S->ti0[t]];
    double e1x = S->vx[S->ti2[t]] - S->vx[S->ti0[t]], e1y = S->vy[S->ti2[t]] - S->vy[S->ti0[t]],
           e1z = S->vz[S->ti2[t]] - S->vz[S->ti0[t]];
    double nx = e0y * e1z - e0z * e1y, ny = e0z * e1x - e0x * e1z, nz = e0x * e1y - e0y * e1x;
    mean += 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
    if (t == 0 || cx < lo[0]) lo[0] = cx;
    if (t == 0 || cy < lo[1]) lo[1] = cy;
    if (t == 0 || cz < lo[2]) lo[2] = cz;
    if (t == 0 || cx > hi[0]) hi[0] = cx;
    if (t == 0 || cy > hi[1]) hi[1] = cy;
    if (t == 0 || cz > hi[2]) hi[2] = cz;
  }
  if (S->ntris == 0) return 1;
  mean /= (double)S->ntris;
  double h = 2.0 * sqrt(mean > 0 ? mean : 1e-12);
  for (int a = 0; a < 3; a++) {
    g->gorg[a] = lo[a] - h;
    g->ginv = 1.0 / h; /* одинаково по осям: h общий */
  }
  int64_t span[3];
  for (int a = 0; a < 3; a++)
    span[a] = (int64_t)((hi[a] + h - g->gorg[a]) * g->ginv) + 2;
  g->gnx = span[0] > 1 ? span[0] : 1;
  g->gny = span[1] > 1 ? span[1] : 1;
  g->gnz = span[2] > 1 ? span[2] : 1;
  g->ncell = g->gnx * g->gny * g->gnz; /* §929-Х1b: dense, без хэша */
  g->dhead = (int32_t *)malloc((size_t)g->ncell * sizeof *g->dhead);
  g->tnext = (int32_t *)malloc((size_t)S->ntris * sizeof *g->tnext);
  if (!g->dhead || !g->tnext) {
    free(g->dhead);
    free(g->tnext);
    g->dhead = NULL;
    g->tnext = NULL;
    return 2;
  }
  for (int64_t s = 0; s < g->ncell; s++)
    g->dhead[s] = -1;
  g->nreg = 0;
  g->rmax = (double *)malloc((size_t)S->ntris * sizeof *g->rmax); /* §929-Х3 */
  if (g->rmax == NULL) {
    free(g->dhead); /* §929-Х1b */
    free(g->tnext);
    g->dhead = NULL;
    g->tnext = NULL;
    return 2;
  }
  for (uint32_t t = 0; t < S->ntris; t++) {
    double cx = (S->vx[S->ti0[t]] + S->vx[S->ti1[t]] + S->vx[S->ti2[t]]) / 3.0;
    double cy = (S->vy[S->ti0[t]] + S->vy[S->ti1[t]] + S->vy[S->ti2[t]]) / 3.0;
    double cz = (S->vz[S->ti0[t]] + S->vz[S->ti1[t]] + S->vz[S->ti2[t]]) / 3.0;
    { /* §929-Х3: радиус треугольника от центроида */
      double ax = S->vx[S->ti0[t]] - cx, ay = S->vy[S->ti0[t]] - cy, az = S->vz[S->ti0[t]] - cz;
      double bx = S->vx[S->ti1[t]] - cx, by = S->vy[S->ti1[t]] - cy, bz = S->vz[S->ti1[t]] - cz;
      double ex = S->vx[S->ti2[t]] - cx, ey = S->vy[S->ti2[t]] - cy, ez = S->vz[S->ti2[t]] - cz;
      double ra = sqrt(ax * ax + ay * ay + az * az), rb = sqrt(bx * bx + by * by + bz * bz),
             rc = sqrt(ex * ex + ey * ey + ez * ez);
      g->rmax[t] = ra > rb ? (ra > rc ? ra : rc) : (rb > rc ? rb : rc);
    }
    int64_t ix = (int64_t)((cx - g->gorg[0]) * g->ginv),
            iy = (int64_t)((cy - g->gorg[1]) * g->ginv),
            iz = (int64_t)((cz - g->gorg[2]) * g->ginv);
    if (ix < 0) ix = 0; /* центроиды внутри [lo,hi] с запасом — на всякий */
    if (iy < 0) iy = 0;
    if (iz < 0) iz = 0;
    if (ix >= g->gnx) ix = g->gnx - 1;
    if (iy >= g->gny) iy = g->gny - 1;
    if (iz >= g->gnz) iz = g->gnz - 1;
    int64_t s = ix + g->gnx * (iy + g->gny * iz);
    g->tnext[t] = g->dhead[s];
    g->dhead[s] = (int32_t)t;
    g->nreg++;
  }
  g->gmaxr = 0.0; /* §930-А */
  for (uint32_t t = 0; t < S->ntris; t++)
    if (g->rmax[t] > g->gmaxr) g->gmaxr = g->rmax[t];
  return 0;
}

static void pg924_grid_free(pg924_grid *g) {
  free(g->dhead); /* §929-Х1b */
  free(g->tnext);
  free(g->rmax); /* §929-Х3 */
}

/* §929-Х1: статические таблицы обхода: смещения (dx,dy,dz) всех
 * оболочек 0..RMAX в порядке тройного цикла (mr==rad) + границы колец —
 * тот же порядок посещения ячеек, что до оптимизации (побитово). */
static const int16_t (*pg924_shells(const uint8_t **shring))[3] {
  static int16_t tab[(2 * PG924_RMAX + 1) * (2 * PG924_RMAX + 1) * (2 * PG924_RMAX + 1)][3];
  static uint8_t ring[(2 * PG924_RMAX + 1) * (2 * PG924_RMAX + 1) * (2 * PG924_RMAX + 1)];
  static int n = -1;
  if (n < 0) {
    n = 0;
    for (int rad = 0; rad <= PG924_RMAX; rad++) {
      for (int dz = -rad; dz <= rad; dz++)
        for (int dy = -rad; dy <= rad; dy++)
          for (int dx = -rad; dx <= rad; dx++) {
            int mx = dx < 0 ? -dx : dx, my = dy < 0 ? -dy : dy, mz = dz < 0 ? -dz : dz;
            int mr = mx > my ? (mx > mz ? mx : mz) : (my > mz ? my : mz);
            if (mr != rad) continue;
            tab[n][0] = (int16_t)dx;
            tab[n][1] = (int16_t)dy;
            tab[n][2] = (int16_t)dz;
            ring[n] = (uint8_t)rad;
            n++;
          }
    }
  }
  if (shring != NULL) *shring = ring;
  return tab;
}

/* треугольник уровня, содержащий q: PIP (точный) либо ближайший
 * центроид (fallback, считается). Возврат: id треугольника или -1.
 * §929: обход по таблице (тот же порядок); ПИП — только если радиус
 * треугольника достаёт до q (консервативный предфильтр, Х3). */
static int32_t pg924_find(const pg924_grid *g, const hz_kit_level *S, const double q[3],
                          double ptol, int *fallback) {
  int64_t c0[3];
  double bestD = HUGE_VAL; /* точный (PIP): min расстояния до плоскости */
  double bd2 = HUGE_VAL;   /* fallback: min расстояния до центроида */
  int32_t bestPIP = -1, bestNear = -1;
  *fallback = 0;
  for (int a = 0; a < 3; a++)
    c0[a] = (int64_t)((q[a] - g->gorg[a]) * g->ginv);
  const uint8_t *ring;
  const int16_t (*sh)[3] = pg924_shells(&ring);
  const int nsh = (2 * PG924_RMAX + 1) * (2 * PG924_RMAX + 1) * (2 * PG924_RMAX + 1);
  for (int si = 0; si < nsh; si++) {
    int64_t ix = c0[0] + sh[si][0], iy = c0[1] + sh[si][1], iz = c0[2] + sh[si][2];
    if (ix < 0 || iy < 0 || iz < 0 || ix >= g->gnx || iy >= g->gny || iz >= g->gnz) continue;
    int64_t s = ix + g->gnx * (iy + g->gny * iz); /* §929-Х1b: dense */
    {
      {
        for (int32_t t = g->dhead[s]; t >= 0; t = g->tnext[t]) {
          uint32_t i0 = S->ti0[t], i1 = S->ti1[t], i2 = S->ti2[t];
          double x0 = S->vx[i0], y0 = S->vy[i0], z0 = S->vz[i0];
          double gx = (x0 + S->vx[i1] + S->vx[i2]) / 3.0;
          double gy = (y0 + S->vy[i1] + S->vy[i2]) / 3.0;
          double gz = (z0 + S->vz[i1] + S->vz[i2]) / 3.0;
          double ddx = gx - q[0], ddy = gy - q[1], ddz = gz - q[2];
          double dd2 = ddx * ddx + ddy * ddy + ddz * ddz;
          if (dd2 < bd2) {
            bd2 = dd2;
            bestNear = t;
          }
          { /* §929-Х3: дальний по центроиду треугольник своего радиуса
             * не достанет — ПИП невозможен, пропускаем арифметику */
            double rr = g->rmax[t] + ptol;
            if (dd2 > rr * rr) continue;
          }
          double e1x = S->vx[i1] - x0, e1y = S->vy[i1] - y0, e1z = S->vz[i1] - z0;
          double e2x = S->vx[i2] - x0, e2y = S->vy[i2] - y0, e2z = S->vz[i2] - z0;
          double wx = q[0] - x0, wy = q[1] - y0, wz = q[2] - z0;
          double nx = e1y * e2z - e1z * e2y, ny = e1z * e2x - e1x * e2z, nz = e1x * e2y - e1y * e2x;
          double n2 = nx * nx + ny * ny + nz * nz;
          if (n2 > 1e-30) {
            double u = (wx * (e2y * nz - e2z * ny) + wy * (e2z * nx - e2x * nz) +
                        wz * (e2x * ny - e2y * nx)) /
                       n2;
            double v = (wx * (ny * e1z - nz * e1y) + wy * (nz * e1x - nx * e1z) +
                        wz * (nx * e1y - ny * e1x)) /
                       n2;
            double dn = wx * nx + wy * ny + wz * nz;
            double d2 = dn * dn / n2;
            if (u >= -1e-9 && v >= -1e-9 && u + v <= 1.0 + 1e-9 && d2 <= ptol * ptol &&
                d2 < bestD) {
              bestD = d2; /* первый попавшийся достаточно: тай-брейк
                           * копланарных не влияет на агрегаты */
              bestPIP = t;
            }
          }
        }
      }
    }
    /* §929: выход после ПОЛНОГО кольца — как исходный break по rad */
    if (bestPIP >= 0 && (si + 1 == nsh || ring[si + 1] != ring[si])) break;
    /* §930-А: кольцо rad ≥ pipdead — ПИП невозможен (радиус уровня не
     * достаёт); тогда если bd2 ≤ (rad−1)²h², ближе никто не найдётся —
     * выход. Оба условия доказанные отсекатели — результат побитово. */
    {
      int rad = ring[si];
      int pipdead = (int)((g->gmaxr + ptol) * g->ginv) + 2;
      /* §930-А1656: c0 считается УСЕЧЕНИЕМ к нулю (не floor) — q за
       * пределами сетки лежит до h НИЖЕ вычисленной ячейки, истинная
       * дистанция кольца rad деградирует до (rad−2)h; граница с
       * запасом (rad−2), доказано для всех q */
      if (rad >= pipdead && bd2 < HUGE_VAL &&
          bd2 <= (double)(rad - 2) * (double)(rad - 2) / (g->ginv * g->ginv))
        break;
    }
  }
  if (bestPIP >= 0) return bestPIP;
  *fallback = 1;
  return bestNear;
}

/* Сборка представителей: repof по уровням (PIP/fallback), дети-CSR,
 * агрегаты, расширение trivert/tribox (id rep = nt + j). 0 — ок. */
static int pg_reps_build(pg_reps *R, const hz_kit *kk, int32_t nt, const double *area,
                         const double *kd, const double *lep, double le, double **tv9_io,
                         double **tb6_io) {
  memset(R, 0, sizeof *R);
  R->nlev_r = kk->nlev - 1;
  if (R->nlev_r <= 0 || nt <= 0) return 1;
  int32_t *off = (int32_t *)calloc((size_t)kk->nlev, sizeof *off); /* [li] → глобальная база */
  R->repof = (int32_t *)malloc((size_t)R->nlev_r * (size_t)nt * sizeof *R->repof);
  if (!off || !R->repof) {
    free(off);
    free(R->repof);
    R->repof = NULL;
    return 2;
  }
  for (int32_t j = 0; j < R->nlev_r * nt; j++)
    R->repof[j] = -1;
  /* глобальные id: L1 с нуля (L0 — НЕ представитель) */
  for (int32_t li = 2; li < kk->nlev; li++)
    off[li] = off[li - 1] + (int32_t)(int64_t)kk->lev[li - 1].ntris;
  R->nrep = (kk->nlev >= 2) ? off[kk->nlev - 1] + (int32_t)(int64_t)kk->lev[kk->nlev - 1].ntris : 0;
  const double *tv9 = *tv9_io; /* точные центроиды исходных треугольников */
  for (int32_t li = 1; li < kk->nlev; li++) {
    const hz_kit_level *S = &kk->lev[li];
    pg924_grid g;
    if (pg924_grid_build(&g, S) != 0) {
      free(off);
      free(R->repof);
      R->repof = NULL;
      return 2;
    }
    double ptol = 0.25 / g.ginv; /* четверть клетки: PIP в плоскости;
                                  * допуск в метрах принимал ВЫСОТЫ (колонны
                                  * проецировались на пол, ловля §924-дым) */
    for (int32_t t = 0; t < nt; t++) {
      int fb = 0;
      double qc[3] = {
          (tv9[9 * (int64_t)t] + tv9[9 * (int64_t)t + 3] + tv9[9 * (int64_t)t + 6]) / 3.0,
          (tv9[9 * (int64_t)t + 1] + tv9[9 * (int64_t)t + 4] + tv9[9 * (int64_t)t + 7]) / 3.0,
          (tv9[9 * (int64_t)t + 2] + tv9[9 * (int64_t)t + 5] + tv9[9 * (int64_t)t + 8]) / 3.0};
      int32_t r = pg924_find(&g, S, qc, ptol, &fb);
      if (r >= 0) {
        R->repof[(li - 1) * nt + t] = off[li] + r;
        if (fb) R->nfall++;
      }
    }
    { /* §924-дых: раскладка уровня */
      int32_t used = 0;
      int32_t *seen = (int32_t *)calloc((size_t)S->ntris, sizeof *seen);
      if (seen == NULL) { /* без seen раскладка уровня не считается — fail closed */
        free(off);
        return 2;
      }
      for (int32_t t = 0; t < nt; t++) {
        int32_t r = R->repof[(li - 1) * nt + t];
        if (r >= 0) {
          used++;
          if (r - off[li] >= 0 && (uint32_t)(r - off[li]) < S->ntris) seen[r - off[li]]++;
        }
      }
      int32_t mx = 0, arg = -1;
      for (uint32_t t = 0; t < S->ntris; t++)
        if (seen[t] > mx) {
          mx = seen[t];
          arg = (int32_t)t;
        }
      fprintf(stderr, "DBG924 L%d: назначено %d/%d, maxkids=%d (tri %d из %u)\n", li, used, nt, mx,
              arg, S->ntris);
      if (arg >= 0 && mx > nt / 2) { /* аномалия: геометрия победителя */
        uint32_t i0 = S->ti0[arg], i1 = S->ti1[arg], i2 = S->ti2[arg];
        double blo[3] = {1e30, 1e30, 1e30}, bhi[3] = {-1e30, -1e30, -1e30};
        uint32_t vi[3] = {i0, i1, i2};
        for (int v = 0; v < 3; v++) {
          double p[3] = {S->vx[vi[v]], S->vy[vi[v]], S->vz[vi[v]]};
          for (int a = 0; a < 3; a++) {
            if (p[a] < blo[a]) blo[a] = p[a];
            if (p[a] > bhi[a]) bhi[a] = p[a];
          }
        }
        fprintf(stderr,
                "DBG924tri bbox=[%.2f..%.2f %.2f..%.2f %.2f..%.2f] "
                "v=(%.2f,%.2f,%.2f)(%.2f,%.2f,%.2f)(%.2f,%.2f,%.2f)\n",
                blo[0], bhi[0], blo[1], bhi[1], blo[2], bhi[2], S->vx[i0], S->vy[i0], S->vz[i0],
                S->vx[i1], S->vy[i1], S->vz[i1], S->vx[i2], S->vy[i2], S->vz[i2]);
      }
      free(seen);
    }
    if (getenv("HZ_DBG924")) { /* одна точка под микроскопом */
      int32_t t = nt / 2;
      double qc[3] = {
          (tv9[9 * (int64_t)t] + tv9[9 * (int64_t)t + 3] + tv9[9 * (int64_t)t + 6]) / 3.0,
          (tv9[9 * (int64_t)t + 1] + tv9[9 * (int64_t)t + 4] + tv9[9 * (int64_t)t + 7]) / 3.0,
          (tv9[9 * (int64_t)t + 2] + tv9[9 * (int64_t)t + 5] + tv9[9 * (int64_t)t + 8]) / 3.0};
      int fb2 = 0;
      int32_t r = pg924_find(&g, S, qc, ptol, &fb2);
      int64_t cc[3];
      for (int a = 0; a < 3; a++)
        cc[a] = (int64_t)((qc[a] - g.gorg[a]) * g.ginv);
      fprintf(stderr,
              "DBG924one L%d t=%d q=(%.3f,%.3f,%.3f) cell=(%lld,%lld,%lld) h=%.3f ptol=%.4f -> %d "
              "fb=%d\n",
              li, t, qc[0], qc[1], qc[2], (long long)cc[0], (long long)cc[1], (long long)cc[2],
              1.0 / g.ginv, ptol, r, fb2);
      for (int probe = 0; probe < 3; probe++) {
        int32_t t2 = probe == 0 ? 1000 : (probe == 1 ? 40000 : 90000);
        if (t2 >= nt) continue;
        double q2[3] = {
            (tv9[9 * (int64_t)t2] + tv9[9 * (int64_t)t2 + 3] + tv9[9 * (int64_t)t2 + 6]) / 3.0,
            (tv9[9 * (int64_t)t2 + 1] + tv9[9 * (int64_t)t2 + 4] + tv9[9 * (int64_t)t2 + 7]) / 3.0,
            (tv9[9 * (int64_t)t2 + 2] + tv9[9 * (int64_t)t2 + 5] + tv9[9 * (int64_t)t2 + 8]) / 3.0};
        int fb3 = 0;
        int32_t r3 = pg924_find(&g, S, q2, ptol, &fb3);
        fprintf(stderr, "DBG924one L%d t=%d q=(%.3f,%.3f,%.3f) -> %d fb=%d repof=%d\n", li, t2,
                q2[0], q2[1], q2[2], r3, fb3, R->repof[(li - 1) * nt + t2]);
      }
    }
    pg924_grid_free(&g);
  }
  /* дети-CSR + агрегаты */
  int32_t *cnt = (int32_t *)calloc((size_t)R->nrep + 1, sizeof *cnt);
  R->rarea = (double *)calloc((size_t)R->nrep, sizeof *R->rarea);
  R->rrho = (double *)calloc((size_t)R->nrep, sizeof *R->rrho);
  R->rle = (double *)calloc((size_t)R->nrep, sizeof *R->rle);
  R->rle_area_skip = (double *)calloc((size_t)R->nrep, sizeof *R->rle_area_skip);
  R->ratri = (double *)calloc((size_t)R->nrep, sizeof *R->ratri);
  R->rnrm = (double *)calloc((size_t)3 * (size_t)(R->nrep > 0 ? R->nrep : 1), sizeof *R->rnrm);
  R->kmem = (int32_t *)malloc((size_t)nt * (size_t)R->nlev_r * sizeof *R->kmem);
  if (!cnt || !R->rarea || !R->rrho || !R->rle || !R->kmem || !R->rle_area_skip || !R->ratri ||
      !R->rnrm) { /* один отказ-путь: всё снято, течей нет */
    free(off);
    free(cnt);
    free(R->repof);
    free(R->rarea);
    free(R->rrho);
    free(R->rle);
    free(R->rle_area_skip);
    free(R->ratri);
    free(R->rnrm);
    free(R->kmem);
    R->repof = NULL;
    return 2;
  }
  for (int32_t j = 0; j < R->nlev_r * nt; j++)
    if (R->repof[j] >= 0) cnt[R->repof[j] + 1]++;
  for (int32_t r = 0; r < R->nrep; r++)
    cnt[r + 1] += cnt[r];
  R->koff = cnt;
  {
    int32_t *fill = (int32_t *)malloc((size_t)R->nrep * sizeof *fill);
    if (!fill) {
      free(off);
      return 2;
    }
    for (int32_t r = 0; r < R->nrep; r++)
      fill[r] = R->koff[r];
    for (int32_t li = 1; li < kk->nlev; li++)
      for (int32_t t = 0; t < nt; t++) {
        int32_t r = R->repof[(li - 1) * nt + t];
        if (r < 0) continue;
        R->kmem[fill[r]++] = t;
        R->rarea[r] += area[t];
        R->rrho[r] += kd[t] * area[t];
        if (lep && lep[t] > 0.0)
          R->rle_area_skip[r] += area[t]; /* §924: эмиттеры —
                                           * свой le, не rep'а */
        else if (lep)
          R->rle[r] += lep[t] * area[t];
      }
    free(fill);
  }
  for (int32_t r = 0; r < R->nrep; r++) {
    if (R->rarea[r] > 0) {
      double anon = R->rarea[r] - R->rle_area_skip[r]; /* площадь не-эмиттеров */
      R->rrho[r] /= R->rarea[r];
      /* §924: как front_le = le + lep[p]: ГЛОБАЛЬНЫЙ le плюс средняя Ke
       * (ловля: rle без le гасил источники, emitted=0) */
      R->rle[r] = le + (lep ? (anon > 0 ? R->rle[r] / anon : 0.0) : 0.0);
    } else {
      R->rrho[r] = kd[0];
      R->rle[r] = le;
    }
    int32_t k = R->koff[r + 1] - R->koff[r];
    if (k > R->maxkids) R->maxkids = k;
  }
  /* расширение trivert/tribox: id rep = nt + j */
  double *xtv9 = (double *)realloc(*tv9_io, (size_t)(nt + R->nrep) * 9 * sizeof *xtv9);
  double *xtb6 = (double *)realloc(*tb6_io, (size_t)(nt + R->nrep) * 6 * sizeof *xtb6);
  if (!xtv9 || !xtb6) {
    free(off);
    return 2;
  }
  *tv9_io = xtv9;
  *tb6_io = xtb6;
  for (int32_t li = 1; li < kk->nlev; li++) {
    const hz_kit_level *S = &kk->lev[li];
    for (uint32_t t = 0; t < S->ntris; t++) {
      int64_t id = (int64_t)nt + off[li] + (int64_t)t;
      double pp[3][3];
      for (int v = 0; v < 3; v++) {
        uint32_t vi = v == 0 ? S->ti0[t] : (v == 1 ? S->ti1[t] : S->ti2[t]);
        pp[v][0] = S->vx[vi];
        pp[v][1] = S->vy[vi];
        pp[v][2] = S->vz[vi];
      }
      for (int aa = 0; aa < 9; aa++)
        xtv9[9 * id + aa] = ((const double *)pp)[aa];
      { /* §924: отношение площадей поверхность/носитель — компенсация
         * перехвата (плоский rep ловит трубку 1 раз, гроздь — N) */
        double e1x = pp[1][0] - pp[0][0], e1y = pp[1][1] - pp[0][1], e1z = pp[1][2] - pp[0][2];
        double e2x = pp[2][0] - pp[0][0], e2y = pp[2][1] - pp[0][1], e2z = pp[2][2] - pp[0][2];
        double nx = e1y * e2z - e1z * e2y, ny = e1z * e2x - e1x * e2z, nz = e1x * e2y - e1y * e2x;
        {
          int64_t rr = off[li] + (int64_t)t;
          R->ratri[rr] = 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
          if (R->ratri[rr] > 1e-30) { /* §927: нормаль носителя */
            double inv = 1.0 / (2.0 * R->ratri[rr]);
            R->rnrm[3 * rr] = nx * inv;
            R->rnrm[3 * rr + 1] = ny * inv;
            R->rnrm[3 * rr + 2] = nz * inv;
          } else {
            R->rnrm[3 * rr] = 0.0;
            R->rnrm[3 * rr + 1] = 1.0;
            R->rnrm[3 * rr + 2] = 0.0;
          }
        }
      }
      for (int a2 = 0; a2 < 3; a2++) {
        double lo = pp[0][a2], hi2 = pp[0][a2];
        for (int v = 1; v < 3; v++) {
          if (pp[v][a2] < lo) lo = pp[v][a2];
          if (pp[v][a2] > hi2) hi2 = pp[v][a2];
        }
        xtb6[6 * id + a2] = lo;
        xtb6[6 * id + 3 + a2] = hi2;
      }
    }
  }
  free(off);
  return 0;
}

int main(int argc, char **argv) {
  const char *path = NULL, *outfile = "img/pgather.ppm";
  /* §875: дефолты потребителя — ПРОДАКШН-МОДЕЛЬ (дихотомия ks §873/874,
   * per-piece эмиссия Ke А1576); useke=0 / ksf=0 — документированный выход
   * к прежнему миру. */
  double scale = 1.0, le = 1.0, rho = -1.0, fov = 60.0, ksf = 1.0;
  int useke = 1;
  double eye[3] = {0, 0, 0}, look[3] = {0, 0, 0};
  double uphint[3] = {0, 0, 0}; /* §932-Б-Ш8: up= — вертикаль сцены; {0} — Z-up (умолчание) */
  int have_up = 0;              /* флаг: up= задан (без FP-сравнений) */
  int iters = 30, lev = 6, tau0 = 0, noprop = 0, ndirs = 26, mort = 1, i, ax;
  double *lep = NULL; /* §874: per-piece эмиссия (ср. Ke); NULL — прежний мир */
  int W = 320, H = 240, k27 = 0;
  int frames = 1, have_eye2 = 0; /* §877: ходьба */
  double eye2[3] = {0, 0, 0}, look2[3] = {0, 0, 0};
  int have_delbox = 0;                            /* §879: разрушаемость-прототип */
  int rgb = 0;                                    /* §889: RGB-рендер */
  int clip = 0;                                   /* §896: кусок = (tri ∩ клетка) */
  const char *efile_out = NULL, *efile_in = NULL; /* §898: E sidecar */
  const char *recvspec = NULL; /* §921: приёмник, список kit:lev:efile через запятую */
  int recvscr = 0;             /* §921: НК — перемешать поле не-эталонных носителей */
  double expmul = 1.0;         /* §890: множитель экспозиции */
  const char *blkfile = NULL;  /* §882: HBLK v1, mmap-сбор */
  const char *kitpath = NULL;  /* §914-Ш4: геометрия из КИТА */
  double zone = -1.0;          /* §915-R3: радиус кольца детальности */
  double loderr = 0.0;         /* §933-И: экранная ошибка ε (px); >0 — кольца
                                * ГЕОМЕТРИЧЕСКИЕ (d0, 2d0, 4d0…), d0 выводится
                                * из ε/H/fov/ребра L0 (§933-ДОП) */
  const char *emipfile = NULL; /* §933-И: mip-наследование поля — L0-sidecar,
                                * E грубого tri = Σ E_i·A_i/ΣA_i по детям */
  int egour = 0;               /* §932-Б-Ш6: Gouraud-интерполяция E при чтении */
  double repzone = -1.0;       /* §931: кольцо reps — дальнее поле от eye */
  int repnearlp = 0;           /* А1666: этажи ближних только из ближних депозитов */
  double *g_rcent = NULL;      /* §931: центроиды−eye [3·nb], слоты (А1686) */
  int adapt = 0;               /* §915-R4: lpacc-адаптив (0 — битово прежний мир) */
  int travel = 0;              /* §932-А: этажи по пустотному пробегу света */
  double tvc = 0.0;            /* §932-А: масштаб travel-шага (tvc=0 — НК, битово) */
  int medium = 0;              /* §932-Б: дальнее поле как статистическая среда */
  double mtau = 0.0;           /* §932-Б/А1737: порог замещения — доля среднего */
  int coldf = 0;               /* §932-Б/Ш9: стартовый этаж микрокусков */
  int cellc = 0;               /* §932-Б/Ш9/Ш12: клеточный приём (0/1/2) */
  int tmap = 0;                /* Ш9: тон-маппинг серого кадра (эксп. по среднему) */
  double colda = 0.0;          /* §932-Б/Ш9: порог площади для старта, м² */
  int ksdiff = 0;              /* §918: Δ=1−ks (диффузный отскок грубит) */
  double cdelta = 1.0;         /* §862: вес приращения аккумулятора */
  int lpceil = 4;              /* §866: потолок этажа (из swee3-канона) */
  uint8_t *kitlvl = NULL;      /* §915-R3: уровень кита на выбранный треугольник [nt] */
  int32_t *kittri = NULL;      /* §933-И: ЛОКАЛЬНЫЙ tri кита выбранного куска [nt]
                                * (для mip-наследования E: (уровень, tri) → дети) */
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
  /* §896: куски (tri ∩ клетка) — объявления в scope функции */
  int32_t NP = 0;
  int32_t *ktri = NULL, *qmtl = NULL, *qtri = NULL;
  double *qcmin = NULL, *qcmax = NULL, *qcent = NULL;
  double *qarea = NULL, *qnrm = NULL, *qkd = NULL, *qks = NULL, *qlep = NULL;
  double *qparea = NULL; /* §911-13-3: площадь родительского tri для depden */
  /* §867: данные mode=3, заполняются до memset(&so) — см. блок trivert */
  double *g_tv9 = NULL, *g_tb6 = NULL;
  double *g_kcal_tab = NULL; /* §928: таблица kcal= (владелец, живёт до конца) */
  uint8_t *g_lparr = NULL;
  int32_t *g_strip_start = NULL, *g_strip_list = NULL; /* §911-13-3 */
  pg_reps reps;                                        /* §924 */
  hz_kit kk; /* §924: живёт до сборки представителей (use_reps) */
  hz_kit_init(&kk);
  int use_reps = 0;            /* §924: выставляется в kit-ветке (adapt, без zone) */
  int use_kit_l0 = 0;          /* §924: кит с уровнями как L0-меш (adapt, без zone) */
  int reps_collect = 0;        /* §928: HZ_KCALDUMP — карты репов для СБОРА калибровки,
                                * представители в свипе НЕ активны (мир §923) */
  const char *kcalpath = NULL; /* §928: kcal= — таблица k(r,ωbin) */
  int gather = 0;              /* §849: 0 — DDA (умолчание), 1 — brute (путь верификации) */
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
    else if (strncmp(argv[i], "up=", 3) == 0) {
      parse3(argv[i] + 3, uphint); /* §932-Б-Ш8: вертикаль сцены (Y-up: up=0,1,0) */
      have_up = 1;
    } else if (strncmp(argv[i], "W=", 2) == 0)
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
    else if (strncmp(argv[i], "Eout=", 5) == 0)
      efile_out = argv[i] + 5; /* §898 */
    else if (strncmp(argv[i], "Ein=", 4) == 0)
      efile_in = argv[i] + 4; /* §898 */
    else if (strncmp(argv[i], "recv=", 5) == 0)
      recvspec = argv[i] + 5; /* §921 */
    else if (strcmp(argv[i], "recvscramble") == 0)
      recvscr = 1; /* §921: НК */
    else if (strncmp(argv[i], "clip=", 5) == 0)
      clip = atoi(argv[i] + 5); /* §896 */
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
    else if (strncmp(argv[i], "kit=", 4) == 0)
      kitpath = argv[i] + 4; /* §914-Ш4: кит вместо OBJ (паритет — Ш4-П1) */
    else if (strncmp(argv[i], "zone=", 5) == 0)
      zone = atof(argv[i] + 5); /* §915-R3: кольца детальности вокруг eye */
    else if (strncmp(argv[i], "loderr=", 7) == 0)
      loderr = atof(argv[i] + 7); /* §933-И: экранная ошибка колец, пиксели */
    else if (strncmp(argv[i], "Emip=", 5) == 0)
      emipfile = argv[i] + 5; /* §933-И: L0-sidecar с mip-наследованием на кольца */
    else if (strncmp(argv[i], "egour=", 6) == 0)
      egour = atoi(argv[i] + 6); /* §932-Б-Ш6: барицентрическая интерполяция E
                                  * при чтении (умолчание 0 — P0 битово) */
    else if (strncmp(argv[i], "repszone=", 9) == 0)
      repzone = atof(argv[i] + 9); /* §931: reps — только дальнее поле (А1663);
                                    * читается независимо от kcal= (А1677-в) */
    else if (strncmp(argv[i], "repnearlp=", 10) == 0)
      repnearlp = atoi(argv[i] + 10); /* А1666: lpacc-фильтр ближних (умолчание 0) */
    else if (strncmp(argv[i], "ksdiff=", 7) == 0)
      ksdiff = atoi(argv[i] + 7); /* §918: диффузность отскока грубит этаж */
    else if (strncmp(argv[i], "adapt=", 6) == 0)
      adapt = atoi(argv[i] + 6); /* §915-R4: 1 — lpacc-адаптив этажей (§862) */
    else if (strncmp(argv[i], "cdelta=", 7) == 0)
      cdelta = atof(argv[i] + 7);
    else if (strncmp(argv[i], "lpceil=", 7) == 0)
      lpceil = atoi(argv[i] + 7);
    else if (strncmp(argv[i], "travel=", 7) == 0)
      travel = atoi(argv[i] + 7); /* §932-А: этажи по пустотному пробегу */
    else if (strncmp(argv[i], "tvc=", 4) == 0)
      tvc = atof(argv[i] + 4); /* §932-А: масштаб travel-шага (cdt члена) */
    else if (strncmp(argv[i], "medium=", 7) == 0)
      medium = atoi(argv[i] + 7); /* §932-Б: статистическая среда дальнего поля */
    else if (strncmp(argv[i], "mtau=", 5) == 0)
      mtau = atof(argv[i] + 5); /* §932-Б/А1737: порог замещения, доля среднего */
    else if (strncmp(argv[i], "coldf=", 6) == 0)
      coldf = atoi(argv[i] + 6); /* §932-Б/Ш9: стартовый этаж микрокусков */
    else if (strncmp(argv[i], "colda=", 6) == 0)
      colda = atof(argv[i] + 6); /* §932-Б/Ш9: порог площади, м² (0 — выкл) */
    else if (strncmp(argv[i], "cellc=", 6) == 0)
      cellc = atoi(argv[i] + 6); /* §932-Б/Ш9/Ш12: 1=RAW-зонд (НК), 2=консервативный */
    else if (strncmp(argv[i], "tmap=", 5) == 0)
      tmap = atoi(argv[i] + 5); /* Ш9: экспозиция+гамма в сером кадре */
    else if (strncmp(argv[i], "kcal=", 5) == 0)
      kcalpath = argv[i] + 5; /* §928: таблица калибровки k(r,ωbin) */
    else
      path = argv[i];
  }
  if (travel && !adapt) { /* §932/А1705: fail-closed дверь — не молчаливый no-op
                           * (прецедент «kit= с nlev>1 требует zone=») */
    fprintf(stderr, "pgather: travel=1 требует adapt=1 (fail-closed, §932)\n");
    return 2;
  }
  if (tvc < 0.0) { /* §932/А1705 */
    fprintf(stderr, "pgather: tvc= отрицательный не допускается (fail-closed, §932)\n");
    return 2;
  }
  if (medium && !adapt) { /* §932-Б/А1738: среда требует адаптив (lpacc_e) */
    fprintf(stderr, "pgather: medium=1 требует adapt=1 (fail-closed, §932-Б)\n");
    return 2;
  }
  if (medium && mtau <= 0.0) {
    fprintf(stderr, "pgather: medium=1 требует mtau>0 (порог замещения, §932-Б)\n");
    return 2;
  }
  if (!path && !blkfile) { /* §902: с blk= OBJ не обязателен — сцена в файле */
    fprintf(stderr, "use: pgather <scene.obj | blk=ФАЙЛ> [lev=N dirs=.. it=N rho=F le=F] "
                    "[eye=X,Y,Z look=X,Y,Z fov=F W=N H=N] [k27=N] [out=ФАЙЛ]\n");
    return 2;
  }
  if (W < 1 || H < 1 || W > PG_WH_MAX || H > PG_WH_MAX || fov < PG_FOV_MIN || fov > PG_FOV_MAX) {
    fprintf(stderr, "use: pgather <scene.obj | blk=ФАЙЛ> [lev=N dirs=.. it=N rho=F le=F] "
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
  memset(&m, 0, sizeof m); /* §902: без OBJ меш пуст (blk-путь его не читает) */
  if (kitpath) {
    /* §914-Ш4: кит → hz_objmesh ПОБАЙТОВО (v/f/fm/mtl); дальше весь путь
     * pgather НЕ меняется — паритет по построению. Требуется F64:
     * округление вершин сломало бы битовые якоря. */
    /* kk объявлена в scope main (§924) */
    FILE *kf = fopen(kitpath, "rb");
    t0 = now_sec();
    int krc = kf != NULL ? hz_kit_load(&kk, kf) : HZ_KIT_E_IO;
    if (kf != NULL) fclose(kf);
    if (krc != HZ_KIT_OK) {
      fprintf(stderr, "pgather: кит не читается: %s (rc=%d)\n", kitpath, krc);
      return 2;
    }
    if (!(kk.flags & HZ_KIT_FLAG_F64)) {
      fprintf(stderr, "pgather: kit= требует F64 (паритет); переиздайте kitmk\n");
      hz_kit_free(&kk);
      return 2;
    }
    if (kk.nlev > 64) {
      fprintf(stderr, "pgather: kit= nlev>64\n");
      hz_kit_free(&kk);
      return 2;
    }
    /* §924: adapt+без zone — меш из L0; HZ_NOREPS — то же, но представители
     * не передаются свипу (диагностический мир §923) */
    use_kit_l0 = (kk.nlev > 1 && zone <= 0 && adapt && !clip && loderr <= 0.0);
    /* §924-М2: до §925 (кольца) reps — ТОЛЬКО явное включение HZ_REPS=1:
     * энерго-канон плоского носителя не закрыт (см. таблицу §924-М2).
     * §928: HZ_KCALDUMP — собрать калибровку в кусочном мире (без HZ_REPS) */
    reps_collect = (use_kit_l0 && getenv("HZ_KCALDUMP") != NULL);
    use_reps = use_kit_l0 && getenv("HZ_REPS") != NULL && !reps_collect;
    if (kk.nlev > 1 && zone <= 0 && !use_kit_l0 && loderr <= 0.0) {
      fprintf(stderr,
              "pgather: kit= с nlev>1 требует zone=R или loderr=E (fail-closed, §915-R3/§933)\n");
      hz_kit_free(&kk);
      return 2;
    }
    if (emipfile != NULL && (loderr <= 0.0 || zone > 0.0)) {
      fprintf(stderr, "pgather: Emip= требует loderr= без zone= (кольцевой меш, §933)\n");
      hz_kit_free(&kk);
      return 2;
    }
    /* §933-И: d0 ЭКРАННОЙ ОШИБКИ — выводится, не подбирается: ребро s0
     * уровня ℓ занимает p = s0·H/(2·d·tan(fov/2)) пикселей; p ≤ loderr
     * даёт уровень floor(log2(d/d0)), d0 = s0·H/(2·loderr·tan(fov/2)).
     * s0 = sqrt(средней площади L0) — характерное ребро базового меша
     * (для равноупомянутых треугольников ребро ~ sqrt(A) с точностью
     * ~1.2; константаloderr=ε поглощает фактор). */
    double lod_d0 = 0.0;
    if (loderr > 0.0) {
      const hz_kit_level *L0 = &kk.lev[0];
      double asum = 0.0;
      for (uint32_t t = 0; t < L0->ntris; t++) {
        double ax0 = L0->vx[L0->ti1[t]] - L0->vx[L0->ti0[t]],
               ay0 = L0->vy[L0->ti1[t]] - L0->vy[L0->ti0[t]],
               az0 = L0->vz[L0->ti1[t]] - L0->vz[L0->ti0[t]];
        double ax1 = L0->vx[L0->ti2[t]] - L0->vx[L0->ti0[t]],
               ay1 = L0->vy[L0->ti2[t]] - L0->vy[L0->ti0[t]],
               az1 = L0->vz[L0->ti2[t]] - L0->vz[L0->ti0[t]];
        double cxp = ay0 * az1 - az0 * ay1, cyp = az0 * ax1 - ax0 * az1,
               czp = ax0 * ay1 - ay0 * ax1;
        asum += 0.5 * sqrt(cxp * cxp + cyp * cyp + czp * czp);
      }
      double s0 = L0->ntris > 0 ? sqrt(asum / (double)L0->ntris) : 1.0;
      if (!(s0 > 0.0)) s0 = 1.0;
      lod_d0 = s0 * (double)H / (2.0 * loderr * tan(fov * M_PI / 360.0));
      if (!(lod_d0 > 0.0)) lod_d0 = 1e30;
    }
    /* §915-R3: многоуровневый кит — кольца вокруг eye: треугольник
     * уровня i берётся, если ring = min(floor(d/zone), nlev-1) == i.
     * Вырожденный nlev=1 идёт прежним путём (битово). */
    uint32_t vbase[64];
    uint32_t nvtot = 0;
    for (int32_t li = 0; li < kk.nlev; li++) {
      vbase[li] = nvtot;
      nvtot += kk.lev[li].nverts;
    }
    uint32_t ntsel = 0;
    int64_t lod_cnt[64]; /* §933-И: счётчики выбранных треугольников уровней */
    memset(lod_cnt, 0, sizeof lod_cnt);
    if (use_kit_l0) { /* §924: меш = L0 целиком; уровни —
                       * представителям (после пирамиды) */
      ntsel = kk.lev[0].ntris;
    } else if (kk.nlev > 1) {
      for (int32_t li = 0; li < kk.nlev; li++) {
        const hz_kit_level *S = &kk.lev[li];
        for (uint32_t t = 0; t < S->ntris; t++) {
          double cx = (S->vx[S->ti0[t]] + S->vx[S->ti1[t]] + S->vx[S->ti2[t]]) / 3.0;
          double cy = (S->vy[S->ti0[t]] + S->vy[S->ti1[t]] + S->vy[S->ti2[t]]) / 3.0;
          double cz2 = (S->vz[S->ti0[t]] + S->vz[S->ti1[t]] + S->vz[S->ti2[t]]) / 3.0;
          double d = sqrt((cx - eye[0]) * (cx - eye[0]) + (cy - eye[1]) * (cy - eye[1]) +
                          (cz2 - eye[2]) * (cz2 - eye[2]));
          int32_t ring;
          if (loderr > 0.0) { /* §933-И: ГЕОМЕТРИЧЕСКИЕ кольца — экранная
                               * ошибка: уровень +1 на УДВОЕНИЕ дистанции
                               * (границы d0, 2d0, 4d0…), не на +zone */
            if (d <= lod_d0)
              ring = 0;
            else {
              ring = (int32_t)floor(log2(d / lod_d0));
              if (ring > kk.nlev - 1) ring = kk.nlev - 1;
            }
          } else {
            ring = (int32_t)(d / zone);
            if (ring > kk.nlev - 1) ring = kk.nlev - 1;
          }
          if (ring == li) {
            ntsel++;
            lod_cnt[li]++;
          }
        }
      }
      if (ntsel == 0) {
        fprintf(stderr, "pgather: кольца зоны пусты (zone=%.3f)\n", zone);
        hz_kit_free(&kk);
        return 2;
      }
      if (loderr > 0.0) { /* §933-И: диагностика экрана — мера П2 */
        printf("КОЛЬЦА-ЭКР: eps=%.2f px d0=%.3f м, полигонов=%u из %u [", loderr, lod_d0,
               (unsigned)ntsel, (unsigned)kk.lev[0].ntris);
        for (int32_t li = 0; li < kk.nlev; li++)
          printf("%sL%d:%lld", li ? " " : "", li, (long long)lod_cnt[li]);
        printf("]\n");
      }
    }
    const hz_kit_level *KL = &kk.lev[0];
    m.nv = (int32_t)(kk.nlev > 1 && !use_kit_l0 ? nvtot : KL->nverts);
    m.nt = (int32_t)(kk.nlev > 1 ? (int32_t)ntsel : (int32_t)KL->ntris);
    m.nmtl = (int32_t)kk.nmtl;
    m.v = malloc(3 * (size_t)m.nv * sizeof *m.v);
    m.f = malloc(3 * (size_t)m.nt * sizeof *m.f);
    m.fm = malloc((size_t)m.nt * sizeof *m.fm);
    m.mtl = calloc((size_t)(m.nmtl > 0 ? m.nmtl : 1), sizeof *m.mtl);
    m.vn = NULL;
    m.vt = NULL;
    m.ft = NULL;
    m.fn = NULL;
    if (!m.v || !m.f || !m.fm || !m.mtl) {
      fprintf(stderr, "pgather: нет памяти (kit→mesh)\n");
      hz_obj_free(&m);
      hz_kit_free(&kk);
      return 2;
    }
    for (int32_t li = 0; li < kk.nlev && (li == 0 || !use_kit_l0); li++) { /* §924: reps — L0 */
      const hz_kit_level *S = &kk.lev[li];
      for (uint32_t v = 0; v < S->nverts; v++) {
        int64_t dv = vbase[li] + v; /* без (int64_t)-обёртки: сложение уже в int64
                                     * (vbase — int64_t), обёртка была пустой
                                     * (bugprone-misplaced-widening-cast) */
        m.v[3 * dv] = S->vx[v];
        m.v[3 * dv + 1] = S->vy[v];
        m.v[3 * dv + 2] = S->vz[v];
      }
    }
    if (kk.nlev == 1 || use_kit_l0) {
      for (int32_t t = 0; t < m.nt; t++) {
        m.f[3 * (int64_t)t] = (int32_t)KL->ti0[t];
        m.f[3 * (int64_t)t + 1] = (int32_t)KL->ti1[t];
        m.f[3 * (int64_t)t + 2] = (int32_t)KL->ti2[t];
        m.fm[t] = (int32_t)KL->tmtl[t];
      }
    } else {
      /* выбранные треугольники по уровням; kitlvl[ti] = уровень;
       * §933-И: kittri[w] = ЛОКАЛЬНЫЙ tri (для mip-наследования E) */
      kitlvl = (uint8_t *)malloc((size_t)m.nt);
      kittri = (int32_t *)malloc((size_t)m.nt * sizeof *kittri);
      if (kitlvl == NULL || kittri == NULL) {
        fprintf(stderr, "pgather: нет памяти (kitlvl/kittri)\n");
        free(kittri);
        kittri = NULL;
        hz_kit_free(&kk);
        hz_obj_free(&m);
        return 2;
      }
      int32_t w = 0;
      for (int32_t li = 0; li < kk.nlev; li++) {
        const hz_kit_level *S = &kk.lev[li];
        for (uint32_t t = 0; t < S->ntris; t++) {
          double cx = (S->vx[S->ti0[t]] + S->vx[S->ti1[t]] + S->vx[S->ti2[t]]) / 3.0;
          double cy = (S->vy[S->ti0[t]] + S->vy[S->ti1[t]] + S->vy[S->ti2[t]]) / 3.0;
          double cz2 = (S->vz[S->ti0[t]] + S->vz[S->ti1[t]] + S->vz[S->ti2[t]]) / 3.0;
          double d = sqrt((cx - eye[0]) * (cx - eye[0]) + (cy - eye[1]) * (cy - eye[1]) +
                          (cz2 - eye[2]) * (cz2 - eye[2]));
          int32_t ring;
          if (loderr > 0.0) { /* §933-И: как в счётном цикле — битово та же формула */
            if (d <= lod_d0)
              ring = 0;
            else {
              ring = (int32_t)floor(log2(d / lod_d0));
              if (ring > kk.nlev - 1) ring = kk.nlev - 1;
            }
          } else {
            ring = (int32_t)(d / zone);
            if (ring > kk.nlev - 1) ring = kk.nlev - 1;
          }
          if (ring != li) continue;
          m.f[3 * (int64_t)w] = (int32_t)(S->ti0[t] + vbase[li]);
          m.f[3 * (int64_t)w + 1] = (int32_t)(S->ti1[t] + vbase[li]);
          m.f[3 * (int64_t)w + 2] = (int32_t)(S->ti2[t] + vbase[li]);
          m.fm[w] = (int32_t)S->tmtl[t];
          kitlvl[w] = (uint8_t)li;
          kittri[w] = (int32_t)t;
          w++;
        }
      }
    }
    for (int32_t mi = 0; mi < m.nmtl; mi++) {
      m.mtl[mi].kd = kk.mtl[mi].kd;
      memcpy(m.mtl[mi].kd3, kk.mtl[mi].kd3, 24);
      memcpy(m.mtl[mi].ks3, kk.mtl[mi].ks3, 24);
      memcpy(m.mtl[mi].ke3, kk.mtl[mi].ke3, 24);
    }
    /* габарит — по вершинам, как у OBJ-пути */
    for (int32_t v = 0; v < m.nv; v++)
      for (ax = 0; ax < 3; ax++) {
        double c = m.v[3 * (int64_t)v + ax];
        if (v == 0 || c < m.lo[ax]) m.lo[ax] = c;
        if (v == 0 || c > m.hi[ax]) m.hi[ax] = c;
      }
    if (!use_kit_l0 && emipfile == NULL)
      hz_kit_free(&kk); /* §924: кит живёт до сборки представителей;
                         * §933-И: Emip — до mip-наследования поля */
    t1 = now_sec();
    printf("СТАТЬЯ kit-загрузка: %.2f с (nt=%d, mtl=%d) [§914-Ш4]\n", t1 - t0, m.nt, m.nmtl);
  } else if (path) {
    t0 = now_sec();
    if (hz_obj_load(&m, path, scale) != 0) {
      fprintf(stderr, "pgather: не читается %s\n", path);
      return 2;
    }
    t1 = now_sec();
    printf("СТАТЬЯ obj-загрузка: %.2f с (nt=%d, mtl=%d)\n", t1 - t0, m.nt, m.nmtl);
  }

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
      /* явная инициализация — анализатор теряет индукцию цикла по ax
       * (FP-класс diam, прецедент pg_bbox_span/§904) */
      double fwd[3] = {0, 0, 0}, right[3] = {0, 0, 0}, up[3] = {0, 0, 0}, tmpv[3] = {0, 0, 1};
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
      double *LumcB[3] = {NULL, NULL, NULL}; /* §904: RGB из файла */
      if (rgb) {
        if (!hb.kd3) {
          fprintf(stderr, "pgather: rgb=1 требует HBLK3 (v3) - в файле нет kd3/lep3\n");
          return 2;
        }
        for (int ch = 0; ch < 3; ch++) {
          LumcB[ch] = (double *)malloc((size_t)W * (size_t)H * sizeof *LumcB[ch]);
          if (!LumcB[ch]) return 2;
        }
      }
      /* §886: параллелизм §878/§884 ВКЛЮЧЁН (указание: делать все
       * оптимизации) — строки кадра независимы, файл read-only */
#pragma omp parallel for schedule(dynamic, 16) private(ax) reduction(+ : tested2, nhit2)
      for (int iy = 0; iy < H; iy++)
        for (int ix = 0; ix < W; ix++) {
          double sx = (2.0 * (ix + 0.5) / W - 1.0) * tanf;
          double sy = (1.0 - 2.0 * (iy + 0.5) / H) * tanf * (double)H / (double)W;
          double rd[3] = {0, 0, 0}, nn, Lv = 0.0;
          double Lvc[3] = {0, 0, 0}; /* §904 */
          int32_t hit;
          int64_t tested_loc = 0;
          for (ax = 0; ax < 3; ax++)
            rd[ax] = fwd[ax] + sx * right[ax] + sy * up[ax];
          nn = sqrt(rd[0] * rd[0] + rd[1] * rd[1] + rd[2] * rd[2]);
          for (ax = 0; ax < 3; ax++)
            rd[ax] /= nn;
          hit = pg_hblk_nearest(&hb, eye, rd, &tested_loc);
          tested2 += tested_loc;
          if (hit >= 0) {
            nhit2++;
            if (rgb) { /* §904: L_c = le + lep_c + kd_c·E/2π (§889 из файла) */
              for (int ch = 0; ch < 3; ch++) {
                Lvc[ch] = le + hb.lep3[3 * (int64_t)hit + ch] +
                          hb.kd3[3 * (int64_t)hit + ch] * hb.E[hit] / (2.0 * M_PI);
                LumcB[ch][(size_t)iy * (size_t)W + (size_t)ix] = Lvc[ch];
              }
              Lv = Lvc[1]; /* зелёный — яркостная метрика (§889) */
            } else {
              /* §902: v2 — освещённый кадр из файла (та же формула, что
               * pg_lcam_hit: le + lep + kd·E/2π); v1 — силуэт */
              Lv = hb.E ? le + hb.lep[hit] + hb.kd[hit] * hb.E[hit] / (2.0 * M_PI) : le;
            }
          }
          lum[(size_t)iy * (size_t)W + (size_t)ix] = Lv;
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
        if (rgb) { /* §904: экспозиция по ключевой яркости + гамма (§889) */
          double expk = 0.0;
          int64_t npos = 0;
          for (i = 0; i < W * H; i++)
            if (LumcB[1][i] > 0) {
              expk += LumcB[1][i];
              npos++;
            }
          expk = expk * expmul / (npos > 0 ? (double)npos : 1.0) + 1e-9;
          for (i = 0; i < W * H; i++)
            for (int ch = 0; ch < 3; ch++) {
              double v = 255.0 * pow(LumcB[ch][i] / expk < 0 ? 0 : LumcB[ch][i] / expk, 1.0 / 2.2);
              if (v > 255.0) v = 255.0;
              unsigned char bb = (unsigned char)v;
              if (fwrite(&bb, 1, 1, f) != 1) return 2;
            }
        } else if (hb.E) { /* §902: v2 — серый кадр, нормировка на max (как сбор) */
          double lmaxv = 0.0;
          for (i = 0; i < W * H; i++)
            if (lum[i] > lmaxv) lmaxv = lum[i];
          for (i = 0; i < W * H; i++) {
            double v = lum[i] / (lmaxv > 0 ? lmaxv : 1.0) * 255.0;
            unsigned char b[3];
            if (v > 255.0) v = 255.0;
            if (v < 0.0) v = 0.0;
            b[0] = b[1] = b[2] = (unsigned char)v;
            fwrite(b, 1, 3, f);
          }
        } else
          for (i = 0; i < W * H; i++) {
            unsigned char b[3];
            b[0] = b[1] = b[2] = (unsigned char)(lum[i] > 0 ? 255 : 0);
            fwrite(b, 1, 3, f);
          }
        fclose(f);
        printf("КАДР: %s записан (HBLK %s)\n", fname,
               rgb    ? "v3, RGB из файла"
               : hb.E ? "v2, поле из файла"
                      : "v1, силуэт");
      }
      for (int ch = 0; ch < 3; ch++)
        free(LumcB[ch]); /* §904 */
      free(lum);
    }
    pg_hblk_close(&hb);
    if (path) hz_obj_free(&m); /* §902: blk-путь может идти без OBJ */
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
  if (clip) {
    /* §896: клип куска — кусок = (tri ∩ клетка). S-H клип по 6 плоскостям */
    int64_t cap2 = 4096, n2 = 0;
    qcmin = (double *)malloc((size_t)4096 * 6 * sizeof *qcmin);
    qcmax = (double *)malloc((size_t)4096 * 6 * sizeof *qcmax);
    qcent = (double *)malloc((size_t)4096 * 3 * sizeof *qcent);
    qmtl = (int32_t *)malloc((size_t)4096 * sizeof *qmtl);
    ktri = (int32_t *)malloc((size_t)4096 * sizeof *ktri);
    qtri = (int32_t *)malloc((size_t)4096 * sizeof *qtri);
    qarea = (double *)malloc((size_t)4096 * sizeof *qarea);
    qparea = (double *)malloc((size_t)4096 * sizeof *qparea);
    qnrm = (double *)malloc((size_t)4096 * 3 * sizeof *qnrm);
    qkd = (double *)malloc((size_t)4096 * sizeof *qkd);
    qks = (double *)malloc((size_t)4096 * sizeof *qks);
    if (useke) qlep = (double *)malloc((size_t)4096 * sizeof *qlep);
    if (!qcmin || !qcmax || !qcent || !qmtl || !qtri || !qarea || !qnrm || !qkd || !qks ||
        !qparea || (useke && !qlep))
      return 2;
    for (i = 0; i < m.nt; i++) {
      double tri[3][3];
      hz_obj_tri(&m, (int32_t)i, tri);
      int64_t i0[3], i1[3];
      for (ax = 0; ax < 3; ax++) {
        double a = 1e30, b = -1e30;
        for (int v = 0; v < 3; v++) {
          if (tri[v][ax] < a) a = tri[v][ax];
          if (tri[v][ax] > b) b = tri[v][ax];
        }
        i0[ax] = (int64_t)((a - m.lo[ax]) / cell);
        i1[ax] = (int64_t)((b - m.lo[ax]) / cell);
        if (i0[ax] < 0) i0[ax] = 0;
        if (i1[ax] > (1 << lev) - 1) i1[ax] = (1 << lev) - 1;
      }
      for (int64_t cz = i0[2]; cz <= i1[2]; cz++)
        for (int64_t cy = i0[1]; cy <= i1[1]; cy++)
          for (int64_t cx = i0[0]; cx <= i1[0]; cx++) {
            double blo[3] = {m.lo[0] + (double)cx * cell, m.lo[1] + (double)cy * cell,
                             m.lo[2] + (double)cz * cell};
            double bhi[3] = {blo[0] + cell, blo[1] + cell, blo[2] + cell};
            double outp[3][3];
            double aa = sh_clip_tri(tri, blo, bhi, outp);
            if (aa <= 1e-14) continue;
            if (n2 >= cap2) {
              cap2 *= 2;
              /* realloc НЕ в сам указатель: при отказе старый блок терялся бы
               * (cppcheck memleakOnRealloc ×11). Проверка ниже осталась */
              void *tmp;
              tmp = (double *)realloc(qcmin, (size_t)cap2 * 6 * sizeof *qcmin);
              if (tmp) qcmin = tmp;
              tmp = (double *)realloc(qcmax, (size_t)cap2 * 6 * sizeof *qcmax);
              if (tmp) qcmax = tmp;
              tmp = (double *)realloc(qcent, (size_t)cap2 * 3 * sizeof *qcent);
              if (tmp) qcent = tmp;
              tmp = (int32_t *)realloc(qmtl, (size_t)cap2 * sizeof *qmtl);
              if (tmp) qmtl = tmp;
              tmp = (int32_t *)realloc(qtri, (size_t)cap2 * sizeof *qtri);
              if (tmp) qtri = tmp;
              tmp = (double *)realloc(qarea, (size_t)cap2 * sizeof *qarea);
              if (tmp) qarea = tmp;
              tmp = (double *)realloc(qparea, (size_t)cap2 * sizeof *qparea);
              if (tmp) qparea = tmp;
              tmp = (double *)realloc(qnrm, (size_t)cap2 * 3 * sizeof *qnrm);
              if (tmp) qnrm = tmp;
              tmp = (double *)realloc(qkd, (size_t)cap2 * sizeof *qkd);
              if (tmp) qkd = tmp;
              tmp = (double *)realloc(qks, (size_t)cap2 * sizeof *qks);
              if (tmp) qks = tmp;
              if (useke) {
                tmp = (double *)realloc(qlep, (size_t)cap2 * sizeof *qlep);
                if (tmp) qlep = tmp;
              }
              if (!qcmin || !qcmax || !qcent || !qmtl || !qtri || !qarea || !qnrm || !qkd || !qks ||
                  !qparea || (useke && !qlep))
                return 2;
            }
            for (int q = 0; q < 3; q++) {
              qcmin[3 * n2 + q] = blo[q];
              qcmax[3 * n2 + q] = bhi[q];
              qcent[3 * n2 + q] = outp[0][q];
            }
            qmtl[n2] = m.fm[i];
            qtri[n2] = (int32_t)i;
            qarea[n2] = aa;
            qparea[n2] = area[i]; /* §911-13-3: depden полоски = площадь tri */
            for (int q = 0; q < 3; q++)
              qnrm[3 * n2 + q] = nrm[3 * (int64_t)i + q];
            qkd[n2] = kd[i];
            qks[n2] = ks[i];
            if (lep) qlep[n2] = lep[i];
            n2++;
          }
    }
    /* куски: cmin/cmax/cent/mtl → кусочные (для build/CSR), остальные
     * piece-массивы соберутся ниже (kd/ks/lep/nrm/area через ktri) */
    NP = (int32_t)n2;
    /* §911-13-3: прибор потери площади клипа — Σ полосок против Σ tri */
    {
      double sq = 0.0, stri = 0.0; /* stri, а не st: внешний st — hz_sw_stat */
      for (i = 0; i < n2; i++)
        sq += qarea[i];
      for (i = 0; i < m.nt; i++)
        stri += area[i];
      fprintf(stderr, "CLIPAREALOG: strips=%d Sq=%.4f Stri=%.4f ratio=%.4f\n", (int)n2, sq, stri,
              stri > 0 ? sq / stri : 0.0);
    }
    /* свап: дальше вся программа работает в терминах КУСКОВ */
    free(area);
    area = qarea; /* qparea НЕ освобождается: so.parea ссылается (§911-13-3) */
    free(nrm);
    nrm = qnrm;
    free(kd);
    kd = qkd;
    free(ks);
    ks = qks;
    free(lep);
    lep = qlep;
    ktri = qtri;
    /* §911-13-3: CSR «tri → полоски» строится ПОСЛЕ Мортона (пермутированные
     * слоты) — см. блок ниже */
    printf("§896: клип — кусков %d из %d треугольников\n", NP, m.nt);
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
      /* §915-R3: ℓ_p = уровень кита треугольника (кольца зоны);
       * вырожденный кит/OBJ — нули, битово прежний мир */
      lparr[ti] = (kitlvl != NULL && ti < m.nt) ? kitlvl[ti] : 0;
    }
    g_tv9 = tv9;
    g_tb6 = tb6;
    g_lparr = lparr;
  }
  if (use_kit_l0) { /* §924: представители (тре-пространство, до Мортона);
                     * заодно расширяет g_tv9/g_tb6 (id rep = nt + j) */
    if (pg_reps_build(&reps, &kk, m.nt, area, kd, lep, le, &g_tv9, &g_tb6) != 0) {
      fprintf(stderr, "pgather: представители не построились — мир без них (fail closed)\n");
      memset(&reps, 0, sizeof reps);
    } else {
      printf("§924 ПРЕДСТАВИТЕЛИ: %d (уровней %d), fallback=%.2f%%, maxkids=%d\n", (int)reps.nrep,
             (int)reps.nlev_r, 100.0 * (double)reps.nfall / ((double)m.nt * (double)reps.nlev_r),
             (int)reps.maxkids);
    }
    hz_kit_free(&kk);
  }
  int32_t nb = clip ? NP : m.nt;
  const double *bcmin = clip ? qcmin : cmin;
  const double *bcmax = clip ? qcmax : cmax;
  const double *bcent = clip ? qcent : cent;
  const int32_t *bmtl = clip ? qmtl : mtl;
  if (hz_pyr_build(&py, nb, bcmin, bcmax, bcent, bmtl, m.lo, m.hi, cell) != 0) {
    fprintf(stderr, "pgather: пирамида не построилась\n");
    return 2;
  }
  if (clip) { /* §896: pcs.tri = исходный tri (ray-tri/tribox по исходнику) */
    for (int32_t pi = 0; pi < nb; pi++)
      py.pcs[pi].tri = ktri[pi];
  }
  if (mort &&
      (hz_pyr_morton(&py) != 0 || hz_pyr_permute(&py, area, sizeof *area) != 0 ||
       hz_pyr_permute(&py, nrm, 3 * sizeof *nrm) != 0 || hz_pyr_permute(&py, kd, sizeof *kd) != 0 ||
       hz_pyr_permute(&py, ks, sizeof *ks) != 0 ||
       (lep && hz_pyr_permute(&py, lep, sizeof *lep) != 0) ||
       (qparea && hz_pyr_permute(&py, qparea, sizeof *qparea) != 0))) { /* §874/А1581, §911-13-3 */
    fprintf(stderr, "pgather: Morton не прошёл\n");
    return 2;
  }
  if ((use_reps || reps_collect) &&
      reps.repof !=
          NULL) { /* §924: слот-перестановка
                   * repof/kmem; §928: нужна и для сбора (repof в слот-пространстве читает свип) */
    for (int32_t li = 0; li < reps.nlev_r; li++)
      if (hz_pyr_permute(&py, reps.repof + (size_t)li * (size_t)m.nt, sizeof(int32_t)) != 0) {
        fprintf(stderr, "pgather: permute repof — не прошёл\n");
        return 2;
      }
    if (py.perm != NULL) { /* kmem: tri → слот */
      int32_t *invp = (int32_t *)malloc((size_t)m.nt * sizeof *invp);
      if (invp != NULL) {
        for (int32_t s2 = 0; s2 < m.nt; s2++)
          invp[py.perm[s2]] = s2;
        int64_t tot = reps.koff[reps.nrep];
        for (int64_t q2 = 0; q2 < tot; q2++)
          if (reps.kmem[q2] >= 0 && reps.kmem[q2] < m.nt) reps.kmem[q2] = invp[reps.kmem[q2]];
        free(invp);
      }
    }
  }
  if (clip) { /* §911-13-3: CSR «tri → полоски» ПОСЛЕ Мортона, в
               * пермутированных слотах (pcs[p].tri = исходный tri) */
    int32_t *sst = (int32_t *)calloc((size_t)m.nt + 1, sizeof *sst);
    int32_t *slst = (int32_t *)malloc((size_t)(nb ? nb : 1) * sizeof *slst);
    int32_t *fill = (int32_t *)malloc((size_t)m.nt * sizeof *fill);
    int32_t pi;
    if (!sst || !slst || !fill) return 2;
    for (pi = 0; pi < nb; pi++)
      sst[py.pcs[pi].tri + 1]++;
    for (i = 0; i < m.nt; i++)
      sst[i + 1] += sst[i];
    for (i = 0; i < m.nt; i++)
      fill[i] = sst[i];
    for (pi = 0; pi < nb; pi++)
      slst[fill[py.pcs[pi].tri]++] = pi;
    free(fill);
    g_strip_start = sst;
    g_strip_list = slst;
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
  so.cellc = cellc; /* §932-Б/Ш9/Ш12: клеточный приём (0 — битово) */
  so.rho = rho;
  so.iters = iters;
  so.tau0 = tau0;
  so.noprop = noprop;
  so.ndirs = ndirs;
  so.mode = 3; /* §867: фронт с точным пересечением (был mode=2) */
  so.build = 1;
  so.vc = 1;
  so.walk = 1;                     /* §867: продакшн-модель А1576/§861 — дефолт потребителя */
  so.ks = ksf > 0.0 ? ks : NULL;   /* §874: ksf=0 — прежний мир (битово, П1) */
  so.lep = lep;                    /* §874: NULL при useke=0 — побитово прежний мир */
  so.parea = clip ? qparea : NULL; /* §911-13-3: NULL без clip — побитово прежний мир */
  so.trivert = g_tv9;
  so.strip_start = g_strip_start; /* §911-13-3: NULL без clip */
  so.strip_list = g_strip_list;
  so.tribox = g_tb6;
  so.lp = g_lparr;
  /* §915-R4: адаптивный дробный аккумулятор §862 на РЕАЛЬНОМ свипе;
   * adapt=0 (умолчание) — so.lpacc NULL, БИТОВО прежний мир */
  if (adapt) {
    so.lpacc = (double *)calloc((size_t)m.nt, sizeof(double));
    if (so.lpacc == NULL) {
      fprintf(stderr, "pgather: нет памяти (lpacc)\n");
      return 2;
    }
    so.cdelta = cdelta;
    so.cold_f = coldf; /* §932-Б/Ш9: холодный старт этажа (0 — выкл) */
    so.cold_a = colda;
    so.lpapply = 1;
    so.lpceil = lpceil;
    if (ksdiff) {
      so.accum_mode = 5; /* §918: Δ = 1−ks_eff — счёт ДИФФУЗНЫХ отскоков */
      so.ks = ks;        /* ks[p] есть при ksf>0; при NULL — все диффузные */
    }
    so.travel = travel; /* §932-А: сюда — только при adapt (гейт выше) */
    so.tvc = tvc;
    if (medium) {      /* §932-Б: аккумулятор Σдепозитов + режим среды */
      if (!use_reps) { /* А1738: среда требует reps-инфраструктуру (kit+HZ_REPS) */
        fprintf(stderr,
                "pgather: medium=1 требует kit= с уровнями и HZ_REPS=1 (fail-closed, §932-Б)\n");
        return 2;
      }
      so.lpacc_e = (double *)calloc((size_t)m.nt, sizeof(double));
      if (so.lpacc_e == NULL) {
        fprintf(stderr, "pgather: нет памяти (lpacc_e)\n");
        return 2;
      }
      so.medium = 1;
      so.med_tau = mtau;
    }
  }
  if ((use_reps || reps_collect) && reps.repof != NULL) { /* §924: представители в свип */
    so.nrep = reps.nrep;
    so.rep_koff = reps.koff;
    so.rep_kmem = reps.kmem;
    so.rep_area = reps.rarea;
    so.rep_rho = reps.rrho;
    so.rep_le = reps.rle;
    so.rep_atri = reps.ratri;
    so.rep_nrm = reps.rnrm;
    so.repof = reps.repof;
    so.rep_nlev = reps.nlev_r;
    so.reps_collect = reps_collect && !use_reps; /* §928 */
    /* §931: КОЛЬЦА ДЕТАЛЬНОСТИ. Предикат «далеко» — центроид куска от
     * глаза КАМЕРЫ (rep_eye без нового ключа, А1663); rep_cent [3·nb] в
     * слот-порядке, координаты МИНУС глаз (А1686: предвычислено здесь,
     * в свипе только сравнение). Независимо от kcal=. */
    if (repzone > 0.0) {
      g_rcent = (double *)malloc(3u * (size_t)nb * sizeof *g_rcent);
      if (g_rcent == NULL) {
        fprintf(stderr, "pgather: нет памяти (rep_cent)\n");
        return 2;
      }
      for (int32_t pr = 0; pr < nb; pr++) {
        const double *vr = g_tv9 + 9 * (size_t)py.pcs[pr].tri;
        for (int qr = 0; qr < 3; qr++)
          g_rcent[3 * pr + qr] = (vr[qr] + vr[3 + qr] + vr[6 + qr]) / 3.0 - eye[qr];
      }
      so.rep_zone = repzone;
      so.rep_cent = g_rcent;
      for (int qr = 0; qr < 3; qr++)
        so.rep_eye[qr] = eye[qr];
      printf("§931: кольцо repszone=%g от глазa (%g,%g,%g), rep_cent=%d слотов\n", repzone, eye[0],
             eye[1], eye[2], (int)nb);
    }
    so.rep_near_lp = repnearlp;         /* А1666 */
    if (use_reps && kcalpath != NULL) { /* §928: таблица калибровки k(r,ωbin) */
      hz_kcal_hdr kh;
      FILE *kf = fopen(kcalpath, "rb");
      double *kt = NULL;
      if (kf == NULL) {
        fprintf(stderr, "pgather: kcal=%s не открылся\n", kcalpath);
        return 2;
      }
      if (fread(&kh, sizeof kh, 1, kf) != 1 || memcmp(kh.magic, "KCAL928", 8) != 0 ||
          kh.nrep != reps.nrep || kh.nw != HZ_KCAL_NW) {
        fprintf(stderr, "pgather: kcal=%s: битый заголовок (nrep %d≠%d или nw %d≠%d)\n", kcalpath,
                (int)kh.nrep, (int)reps.nrep, (int)kh.nw, HZ_KCAL_NW);
        fclose(kf);
        return 2;
      }
      kt = (double *)malloc((size_t)reps.nrep * HZ_KCAL_NW * sizeof *kt);
      if (kt == NULL || fread(kt, sizeof *kt, (size_t)reps.nrep * HZ_KCAL_NW, kf) !=
                            (size_t)reps.nrep * HZ_KCAL_NW) {
        fprintf(stderr, "pgather: kcal=%s: битый массив\n", kcalpath);
        free(kt);
        fclose(kf);
        return 2;
      }
      fclose(kf);
      g_kcal_tab = kt; /* владелец — free в конце (урок CWE-401 §874) */
      so.kcal = kt;
      printf("§928 kcal: %s загружена (nrep=%d nw=%d)\n", kcalpath, (int)reps.nrep,
             (int)HZ_KCAL_NW);
    }
  }
  so.domhi = m.hi;
  double *Ec[3] = {NULL, NULL, NULL}; /* §889: E по каналам */
  double *kdc[3] = {NULL, NULL, NULL}, *lepc[3] = {NULL, NULL, NULL};
  if (rgb) {
    /* §889: три скалярных решения с альбедо/эмиссией канала (классика
     * радиосити). kd_c = kd3[c], lep_c = ke3[c]; le остаётся общим.
     * §901: всё ПО КУСКАМ (nb): массивы каналов длины nb, kd3/ke3 берутся
     * по tri куска (py.pcs[i].tri, валиден i<nb в обоих режимах) — иначе
     * после Morton слот ≠ tri (чужой цвет), а при clip массивы короче nb. */
    if (ksf > 0.0) {
      fprintf(stderr, "pgather: rgb=1 с ksf>0 не совмещается в v1 - отказ\n");
      return 2;
    }
    for (int ch = 0; ch < 3; ch++) {
      kdc[ch] = (double *)malloc((size_t)nb * sizeof *kdc[ch]);
      lepc[ch] = (double *)malloc((size_t)nb * sizeof *lepc[ch]);
      Ec[ch] = (double *)malloc((size_t)nb * sizeof *Ec[ch]);
      if (!kdc[ch] || !lepc[ch] || !Ec[ch]) return 2;
      for (i = 0; i < nb; i++) {
        int32_t tri = py.pcs[i].tri;
        kdc[ch][i] = m.mtl[m.fm[tri]].kd3[ch];
        lepc[ch][i] = m.mtl[m.fm[tri]].ke3[ch];
      }
    }
    t0 = now_sec();
    for (int ch = 0; ch < 3; ch++) {
      so.lep = lepc[ch]; /* kd канала идёт 3-м аргументом hz_sw_run */
      for (i = 0; i < nb; i++)
        py.pcs[i].e = 0.0f;
      if (hz_sw_run(&py, nb, area, nrm, kdc[ch], &so, &st, NULL) != 0) return 2;
      for (i = 0; i < nb; i++)
        Ec[ch][i] = py.pcs[i].e;
    }
    t1 = now_sec();
    sw_time = t1 - t0;
    printf("СВИП RGB: 3 канала (%.3f с), E_avg=%.4f\n", sw_time, st.e_avg);
  } else if (emipfile) {
    /* §933-И: mip-наследование поля (§933-ДОП): E выбранного ГРУБОГО tri
     * кольцевого меша = площадь-взвешенное среднее E его L0-детей
     * (Σ E_i·A_i / Σ A_i — та же арифметика, что rep_ep в свипе).
     * Дети — PIP центроида L0-tri в грубом tri (pg924_grid, как
     * pg_reps_build); L0-куски берут своё значение напрямую. Поле
     * живёт на L0 и наследуется ВНИЗ — грубый уровень НЕ считает. */
    const hz_kit_level *L0k = &kk.lev[0];
    if (clip) {
      fprintf(stderr, "pgather: Emip= несовместим с clip= (fail-closed, §933)\n");
      return 2;
    }
    double *E0 = (double *)malloc((size_t)L0k->ntris * sizeof *E0);
    double *A0 = (double *)malloc((size_t)L0k->ntris * sizeof *A0);
    double *sumE = (double *)calloc((size_t)m.nt, sizeof *sumE);
    double *sumA = (double *)calloc((size_t)m.nt, sizeof *sumA);
    /* selr[li][t] — уровень-локальный tri → выбранный кусок w (-1 нет);
     * плоский массив по vbase-подобным смещениям уровней */
    int64_t ltri_off[HZ_KIT_MAX_NLEV];
    int32_t *selr = NULL;
    {
      int64_t tot = 0;
      for (int32_t li = 0; li < kk.nlev; li++) {
        ltri_off[li] = tot;
        tot += (int64_t)kk.lev[li].ntris;
      }
      selr = (int32_t *)malloc((size_t)tot * sizeof *selr);
    }
    if (!E0 || !A0 || !sumE || !sumA || !selr) {
      fprintf(stderr, "pgather: нет памяти (Emip)\n");
      return 2;
    }
    for (int64_t q = 0; q < ltri_off[kk.nlev - 1] + (int64_t)kk.lev[kk.nlev - 1].ntris; q++)
      selr[q] = -1;
    for (i = 0; i < m.nt; i++)
      selr[ltri_off[kitlvl[i]] + kittri[i]] = (int32_t)i;
    {
      FILE *fm = fopen(emipfile, "rb");
      if (!fm || fread(E0, sizeof(double), (size_t)L0k->ntris, fm) != (size_t)L0k->ntris) {
        fprintf(stderr, "pgather: Emip-файл не читается/короток: %s\n", emipfile);
        if (fm) fclose(fm);
        return 2;
      }
      fclose(fm);
    }
    for (uint32_t t = 0; t < L0k->ntris; t++) { /* площади L0 (кросс-произведение) */
      double ax0 = L0k->vx[L0k->ti1[t]] - L0k->vx[L0k->ti0[t]],
             ay0 = L0k->vy[L0k->ti1[t]] - L0k->vy[L0k->ti0[t]],
             az0 = L0k->vz[L0k->ti1[t]] - L0k->vz[L0k->ti0[t]];
      double ax1 = L0k->vx[L0k->ti2[t]] - L0k->vx[L0k->ti0[t]],
             ay1 = L0k->vy[L0k->ti2[t]] - L0k->vy[L0k->ti0[t]],
             az1 = L0k->vz[L0k->ti2[t]] - L0k->vz[L0k->ti0[t]];
      double cxp = ay0 * az1 - az0 * ay1, cyp = az0 * ax1 - ax0 * az1, czp = ax0 * ay1 - ay0 * ax1;
      A0[t] = 0.5 * sqrt(cxp * cxp + cyp * cyp + czp * czp);
    }
    int64_t nmiss = 0; /* L0-центроид не нашёл выбранного носителя */
    {                  /* сетки PIP — ОДНА на уровень, не на L0-tri (ловля производительности) */
      pg924_grid gm[HZ_KIT_MAX_NLEV];
      int gmok[HZ_KIT_MAX_NLEV];
      for (int32_t li = 1; li < kk.nlev; li++)
        gmok[li] = pg924_grid_build(&gm[li], &kk.lev[li]) == 0;
      for (uint32_t t = 0; t < L0k->ntris; t++) {
        double qc[3] = {(L0k->vx[L0k->ti0[t]] + L0k->vx[L0k->ti1[t]] + L0k->vx[L0k->ti2[t]]) / 3.0,
                        (L0k->vy[L0k->ti0[t]] + L0k->vy[L0k->ti1[t]] + L0k->vy[L0k->ti2[t]]) / 3.0,
                        (L0k->vz[L0k->ti0[t]] + L0k->vz[L0k->ti1[t]] + L0k->vz[L0k->ti2[t]]) / 3.0};
        int32_t own = selr[ltri_off[0] + (int32_t)t]; /* свой L0-кусок выбран? */
        if (own >= 0) { /* битово = прямой путь Ein для L0-части меша */
          sumA[own] = 1.0;
          sumE[own] = E0[t];
          continue;
        }
        /* иначе — носитель-предок на уровне выше: PIP по уровням,
         * берём БЛИЖАЙШИЙ выбранный (минимальный li ≥ 1, где tri выбран) */
        int32_t host = -1;
        for (int32_t li = 1; li < kk.nlev && host < 0; li++) {
          if (!gmok[li]) continue;
          double ptol = 0.25 / gm[li].ginv; /* как pg_reps_build (§924-ловля высот) */
          int fb;
          int32_t r = pg924_find(&gm[li], &kk.lev[li], qc, ptol, &fb);
          if (r >= 0) {
            int32_t cand = selr[ltri_off[li] + r];
            if (cand >= 0) host = cand;
          }
        }
        if (host >= 0) {
          sumE[host] += E0[t] * A0[t];
          sumA[host] += A0[t];
        } else
          nmiss++;
      }
      for (int32_t li = 1; li < kk.nlev; li++)
        if (gmok[li]) pg924_grid_free(&gm[li]);
    }
    int64_t nempty = 0;
    for (i = 0; i < nb; i++) { /* слоты мортон-переставлены — только через
                                * pcs[i].tri (как ветвь Ein, биекция не-clip) */
      int32_t w2 = py.pcs[i].tri;
      double e = kitlvl[w2] == 0 ? sumE[w2] : (sumA[w2] > 0.0 ? sumE[w2] / sumA[w2] : 0.0);
      if (kitlvl[w2] > 0 && sumA[w2] <= 0.0) nempty++;
      py.pcs[i].e = (float)e;
    }
    printf("§933-И: E mip-наследовано из %s (пропало L0=%lld, пустых грубых=%lld)\n", emipfile,
           (long long)nmiss, (long long)nempty);
    free(E0);
    free(A0);
    free(sumE);
    free(sumA);
    free(selr);
    free(kittri);
    kittri = NULL;
    hz_kit_free(&kk); /* §933-И: кит дожил до mip-наследования */
    t0 = now_sec();
  } else if (efile_in) {
    /* §899: E из sidecar — свип пропущен (свободная ходьба).
     * §902: sidecar канонически PER-TRI в ИСХОДНОМ порядке (m.nt double) —
     * чтобы pblock/E-в-файле читали его без знания о Morton; здесь
     * разворачиваем в слоты кусков по pcs[i].tri (не-clip: биекция —
     * битово; clip: куски одного tri делят E — осознанный предел v1). */
    double *Etri = (double *)malloc((size_t)m.nt * sizeof *Etri);
    if (!Etri) return 2;
    FILE *fei = fopen(efile_in, "rb");
    if (!fei) {
      fprintf(stderr, "pgather: E не читается: %s\n", efile_in);
      return 2;
    }
    if (fread(Etri, sizeof(double), (size_t)m.nt, fei) != (size_t)m.nt) {
      fprintf(stderr, "pgather: E короче nt (%s)\n", efile_in);
      fclose(fei);
      return 2;
    }
    fclose(fei);
    for (i = 0; i < nb; i++)
      py.pcs[i].e = (float)Etri[py.pcs[i].tri];
    free(Etri);
    t0 = now_sec();
    printf("§899: E из %s (свип пропущен)\n", efile_in);
  } else {
    for (i = 0; i < m.nt; i++)
      py.pcs[i].e = 0.0f;
    t0 = now_sec();
    if (hz_sw_run(&py, nb, area, nrm, kd, &so, &st, NULL) != 0) {
      fprintf(stderr, "pgather: свип не прошёл\n");
      return 2;
    }
    t1 = now_sec();
    sw_time = t1 - t0;
    printf("СВИП: E_avg=%.4f (%.3f с), кусков с E==0: %lld (инвариант позитивности, Ш9)\n",
           st.e_avg, sw_time, (long long)st.zero_e);
    /* §915-R4: приборы марша — доказательство удешевления этажами
     * (nmat/ndesc/njump) числом, адаптив (flochg) — отдельно */
    printf("МАРШ §852: материальных=%lld спусков=%lld прыжков=%lld nlpclamp=%lld визитов=%lld\n",
           (long long)st.nmat, (long long)st.ndesc, (long long)st.njump, (long long)st.nlpclamp,
           (long long)st.nvisit);
    if (so.lpacc != NULL)
      printf("АДАПТИВ §866: смен этажей Σ=%lld (макс/ит.=%lld)\n", (long long)st.flochg_sum,
             (long long)st.flochg_max);
    if (efile_out) { /* §898: дамп E — §902: per-tri в исходном порядке
                      * (Etri[tri] = E куска; не-clip — биекция) */
      double *Etri = (double *)malloc((size_t)m.nt * sizeof *Etri);
      if (!Etri) return 2;
      for (i = 0; i < m.nt; i++)
        Etri[i] = 0.0;
      for (i = 0; i < nb; i++)
        Etri[py.pcs[i].tri] = py.pcs[i].e;
      FILE *fe = fopen(efile_out, "wb");
      if (!fe) return 2;
      fwrite(Etri, sizeof(double), (size_t)m.nt, fe);
      fclose(fe);
      free(Etri);
      printf("§898: E записан в %s (%d tri)\n", efile_out, m.nt);
    }
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
      /* смещения ЛЕВОЙ части тоже 64-битные: справа уже (int64_t)tri, а 3·p2
       * в int при p2 в сотни миллионов переполнилось бы (ловля 07-10) */
      cmin2[3 * (int64_t)p2] = cmin[3 * (int64_t)tri];
      cmin2[3 * (int64_t)p2 + 1] = cmin[3 * (int64_t)tri + 1];
      cmin2[3 * (int64_t)p2 + 2] = cmin[3 * (int64_t)tri + 2];
      cmax2[3 * (int64_t)p2] = cmax[3 * (int64_t)tri];
      cmax2[3 * (int64_t)p2 + 1] = cmax[3 * (int64_t)tri + 1];
      cmax2[3 * (int64_t)p2 + 2] = cmax[3 * (int64_t)tri + 2];
      cent2[3 * (int64_t)p2] = cent[3 * (int64_t)tri];
      cent2[3 * (int64_t)p2 + 1] = cent[3 * (int64_t)tri + 1];
      cent2[3 * (int64_t)p2 + 2] = cent[3 * (int64_t)tri + 2];
      area2[p2] = area[i];
      nrm2[3 * (int64_t)p2] = nrm[3 * (int64_t)i];
      nrm2[3 * (int64_t)p2 + 1] = nrm[3 * (int64_t)i + 1];
      nrm2[3 * (int64_t)p2 + 2] = nrm[3 * (int64_t)i + 2];
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
  int nrecv = 0; /* §921: приёмник (0 — прежний мир битово) */
  pg_recv *rcv = NULL;
  double *rayL0 = NULL, **rayD = NULL;
  double rayL0sum = 0.0, *rayS = NULL, *rayAD = NULL, *rayMX = NULL;
  double ridmax = 0.0; /* §921: тождество носителя L0 с обычным сбором */
  int64_t ridcnt = 0;
  if (recvspec != NULL) {
    if (rgb) {
      fprintf(stderr, "pgather: recv= не совмещается с rgb=\n");
      return 2;
    }
    if (ksf > 0.0) {
      fprintf(stderr, "pgather: recv= требует ksf=0 (прибор диффузного поля)\n");
      return 2;
    }
    if (frames > 1) {
      fprintf(stderr, "pgather: recv= с frames>1 не поддерживается\n");
      return 2;
    }
    if (efile_in == NULL) {
      fprintf(stderr, "pgather: recv= требует Ein= (sidecar приёмника; тождество)\n");
      return 2;
    }
#ifdef _OPENMP
    omp_set_num_threads(1); /* прибор: детерминизм, счётчики lookup */
#endif
    char spec[8192];
    snprintf(spec, sizeof spec, "%s", recvspec);
    int nr = 1;
    for (const char *p = spec; *p; p++)
      if (*p == ',') nr++;
    rcv = (pg_recv *)calloc((size_t)nr, sizeof *rcv);
    if (rcv == NULL) return 2;
    char *tok = strtok(spec, ",");
    while (tok != NULL) {
      char kp[4096], ep[4096];
      int lvv;
      if (sscanf(tok, "%4095[^:]:%d:%4095s", kp, &lvv, ep) != 3) {
        fprintf(stderr, "pgather: recv= запись '%s' не вида kit:lev:efile\n", tok);
        return 2;
      }
      if (pg_recv_open(&rcv[nrecv], kp, (int32_t)lvv, ep, recvscr && nrecv > 0) != 0) {
        fprintf(stderr, "pgather: recv= носитель %s:%d:%s не открылся\n", kp, lvv, ep);
        return 2;
      }
      nrecv++;
      tok = strtok(NULL, ",");
    }
    if (nrecv < 2) {
      fprintf(stderr, "pgather: recv= нужно ≥2 записи (эталон + носитель)\n");
      return 2;
    }
    rayL0 = (double *)calloc((size_t)W * (size_t)H, sizeof *rayL0);
    rayD = (double **)calloc((size_t)nrecv, sizeof *rayD);
    rayS = (double *)calloc((size_t)nrecv, sizeof *rayS);
    rayAD = (double *)calloc((size_t)nrecv, sizeof *rayAD);
    rayMX = (double *)calloc((size_t)nrecv, sizeof *rayMX);
    if (rayL0 == NULL || rayD == NULL || rayS == NULL || rayAD == NULL || rayMX == NULL) return 2;
    for (int rk = 1; rk < nrecv; rk++) {
      rayD[rk] = (double *)calloc((size_t)W * (size_t)H, sizeof *rayD[rk]);
      if (rayD[rk] == NULL) return 2;
    }
    printf("§921 ПРИЁМНИК: %" PRId64 " носителей; lookup PIP+fallback, сериально\n",
           (int64_t)nrecv);
  }
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
      if (have_up) {
        tmp[0] = uphint[0]; /* §932-Б-Ш8: up= — подсказка вертикали сцены
                             * (Y-up сцены: up=0,1,0; умолчание Z-up битово) */
        tmp[1] = uphint[1];
        tmp[2] = uphint[2];
      }
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
      if (have_up) { /* §932-Б-Ш8: right = tmp×fwd — правый базис с вертикалью
                      * tmp; прежняя формула fwd×tmp давала зеркальный кадр
                      * (ловля по визуальной оценке пользователя) */
        right[0] = tmp[1] * fwd[2] - tmp[2] * fwd[1];
        right[1] = tmp[2] * fwd[0] - tmp[0] * fwd[2];
        right[2] = tmp[0] * fwd[1] - tmp[1] * fwd[0];
      } else {
        right[0] = fwd[1] * tmp[2] - fwd[2] * tmp[1];
        right[1] = fwd[2] * tmp[0] - fwd[0] * tmp[2];
        right[2] = fwd[0] * tmp[1] - fwd[1] * tmp[0];
      }
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
      { /* §932-Б-Ш6: пер-вершинное E — площадь-взвешенное среднее смежных
         * кусков (egour=1); отказ памяти — fail closed к P0 (ev=NULL) */
        double *ev = NULL, *ew = NULL;
        if (egour) {
          ev = (double *)calloc((size_t)m.nv, sizeof *ev);
          ew = (double *)calloc((size_t)m.nv, sizeof *ew);
          if (ev != NULL && ew != NULL) {
            int32_t p3;
            for (p3 = 0; p3 < nb; p3++) {
              int32_t t3 = py.pcs[p3].tri, k3;
              double ae = area[p3] * (double)py.pcs[p3].e;
              for (k3 = 0; k3 < 3; k3++) {
                int32_t vi = m.f[3 * (int64_t)t3 + k3];
                if (vi >= 0 && vi < m.nv) {
                  ev[vi] += ae;
                  ew[vi] += area[p3];
                }
              }
            }
            for (int32_t vi = 0; vi < m.nv; vi++)
              ev[vi] = ew[vi] > 0.0 ? ev[vi] / ew[vi] : 0.0;
            cam0.ev = ev; /* §932-Б-Ш6: чтение кадра интерполирует */
          } else {
            free(ev);
            free(ew);
          }
        }
        g_egour_ev = ev; /* освободить после петли кадра */
        g_egour_ew = ew;
      }
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
            if (rcv != NULL) { /* §921: инвариантный приёмник — поле с носителей */
              double qh[3], kdvis, Lplain, Lc;
              int a2;
              for (a2 = 0; a2 < 3; a2++)
                qh[a2] = ef[a2] + rd[a2] * thit;
              kdvis = rho < 0 ? kd[pbest] : rho;
              Lplain =
                  le + (lep ? lep[pbest] : 0.0) + kdvis * (double)py.pcs[pbest].e / (2.0 * M_PI);
              size_t pix = (size_t)iy * (size_t)W + (size_t)ix;
              for (int rk = 0; rk < nrecv; rk++) {
                double Ecr = pg_recv_E(&rcv[rk], qh);
                Lc = le + (lep ? lep[pbest] : 0.0) + kdvis * Ecr / (2.0 * M_PI);
                if (rk == 0) {
                  rayL0[pix] = Lc;
                  rayL0sum += Lc;
                  L = Lc;
                  double dd = fabs(Lc - Lplain); /* тождество (гейт П2) */
                  if (dd > ridmax) ridmax = dd;
                  if (dd > 0.0) ridcnt++;
                } else {
                  double dd = fabs(Lc - rayL0[pix]);
                  rayD[rk][pix] = dd;
                  rayS[rk] += Lc;
                  rayAD[rk] += dd;
                  if (dd > rayMX[rk]) rayMX[rk] = dd;
                }
              }
            } else if (rgb) { /* §889: поканальная яркость (зеркальный член — v2) */
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
      free(g_egour_ev); /* §932-Б-Ш6: буферы Gouraud — после петли кадра */
      free(g_egour_ew);
      if (gather == 0)
        printf("DDA: клеток/луч %.1f, кусков/луч %.1f (%.2f %% от nt)\n",
               (double)dda_steps / ((double)W * (double)H),
               (double)dda_tested / ((double)W * (double)H),
               100.0 * (double)dda_tested / ((double)W * (double)H) / (double)m.nt);
      printf("M3b: kuskov v kletke (sredn poseshchennoi) %.1f\n",
             g_m3b_visits ? (double)g_m3b_inlist / (double)g_m3b_visits : 0.0);
      printf("G1: grupp/luch %.2f, promahov %.1f%%, testov-kuskov/luch %.1f\n",
             (double)g_m3c_gtests / ((double)W * (double)H),
             g_m3c_gtests ? 100.0 * (1.0 - (double)g_m3c_ghit / (double)g_m3c_gtests) : 0.0,
             (double)dda_tested / ((double)W * (double)H));
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
      }
      printf("§894: деп_main=%.1f деп_ног=%.1f депозитов=%lld (absorbed=%.1f)\n", st.dep_main,
             st.dep_leg, (long long)st.dep_cnt, st.absorbed);
      {
        double rr = rho < 0 ? 0.5 : rho;
        double Lpred = le + rr * st.e_avg / (2.0 * M_PI);
        double Lmeas = lsum / (nhit ? (double)nhit : 1.0);
        printf("СЛИЧЕНИЕ: L_свипа = le + rho·E_avg/2π = %.4f; L_сбора = %.4f; отношение %.4f\n",
               Lpred, Lmeas, Lmeas / Lpred);
      }
      if (rcv != NULL) { /* §921: таблица носителей на ОДИНИХ лучах */
        double mL0 = rayL0sum / (nhit ? (double)nhit : 1.0);
        printf("§921 ПРИЁМНИК: попаданий %" PRId64 ", <L0>=%.6f; тождество max|dL|=%.3g"
               " (лучей с отличием %" PRId64 ")\n",
               nhit, mL0, ridmax, ridcnt);
        for (int rk = 1; rk < nrecv; rk++) {
          double *v = (double *)malloc((size_t)W * (size_t)H * sizeof *v);
          if (v == NULL) return 2;
          memcpy(v, rayD[rk], (size_t)W * (size_t)H * sizeof *v);
          qsort(v, (size_t)W * (size_t)H, sizeof *v, pg_recv_dcmp);
          double p95 = v[(size_t)(0.95 * (double)((size_t)W * (size_t)H - 1))];
          free(v);
          printf("  носитель L%-2d (%u tri): <L>=%.5f  Dmean=%.2f%%  d_ray=%.2f%%"
                 "  p95=%.2f%%  max=%.2f%%  fallback=%.2f%%  empty=%" PRId64 "\n",
                 (int)rcv[rk].lv, rcv[rk].S->ntris, rayS[rk] / (nhit ? (double)nhit : 1.0),
                 100.0 * fabs(rayS[rk] - rayL0sum) / (nhit ? (double)nhit : 1.0) /
                     (mL0 > 0 ? mL0 : 1.0),
                 100.0 * rayAD[rk] / (nhit ? (double)nhit : 1.0) / (mL0 > 0 ? mL0 : 1.0),
                 100.0 * p95 / (mL0 > 0 ? mL0 : 1.0), 100.0 * rayMX[rk] / (mL0 > 0 ? mL0 : 1.0),
                 100.0 * (double)rcv[rk].nfall / (double)(rcv[rk].nlook ? rcv[rk].nlook : 1),
                 rcv[rk].nempty);
          { /* диф-картинка |dL| (правило картинок: img/) */
            char fn[4096];
            FILE *g;
            snprintf(fn, sizeof fn, "img/recv_dL%d.ppm", (int)rcv[rk].lv);
            g = fopen(fn, "wb");
            if (g == NULL) {
              fprintf(stderr, "pgather: не открыть %s\n", fn);
              return 2;
            }
            fprintf(g, "P6\n%d %d\n255\n", W, H);
            for (size_t ip2 = 0; ip2 < (size_t)W * (size_t)H; ip2++) {
              /* видимость: 2·|dL|/<L0> -> 255 (d=50% — насыщение) */
              int b = (int)(510.0 * rayD[rk][ip2] / (mL0 > 0 ? mL0 : 1.0));
              unsigned char px[3];
              if (b > 255) b = 255;
              if (b < 0) b = 0;
              px[0] = px[1] = px[2] = (unsigned char)b;
              fwrite(px, 1, 3, g);
            }
            fclose(g);
            printf("  диф-картинка: %s\n", fn);
          }
        }
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
        } else if (tmap) {
          /* §932-Б/Ш9: экспозиция+гамма в СЕРОМ пути (как rgb): expk =
           * среднее положительных × expmul; лечит «хвост давит кадр»
           * (лампы Ke=30 делали max ~2000 и чёрный кадр); умолчание
           * выкл — прежняя нормировка на max, битово */
          double expk = 0.0;
          int64_t npos = 0;
          for (i = 0; i < W * H; i++)
            if (lum[i] > 0) {
              expk += lum[i];
              npos++;
            }
          expk = expk * expmul / (npos > 0 ? (double)npos : 1.0) + 1e-9;
          for (i = 0; i < W * H; i++) {
            double v = 255.0 * pow(lum[i] / expk < 0 ? 0 : lum[i] / expk, 1.0 / 2.2);
            unsigned char b[3];
            if (v > 255.0) v = 255.0;
            b[0] = b[1] = b[2] = (unsigned char)v;
            fwrite(b, 1, 3, f);
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
        printf("КАДР: %s записан (P6, %s)\n", fname,
               tmap ? "экспозиция по среднему + гамма (tmap)" : "нормировка на max кадра");
      }
    } /* базис камеры */
  } /* §877: ходьба */

  free(lum);
  free(reps.repof); /* §924 */
  free(reps.koff);
  free(reps.kmem);
  free(reps.rarea);
  free(reps.rrho);
  free(reps.rle);
  free(reps.rle_area_skip);
  free(reps.ratri);
  free(reps.rnrm);
  free(g_tv9); /* §874: утечка trivert/tribox/lp (предсуществующая с §867,
                * поймана ASAN при прогоне §874) */
  free(g_tb6);
  free(g_rcent);    /* §931 */
  free(g_kcal_tab); /* §928 */
  free(g_lparr);
  free(kitlvl);      /* §915-R3 */
  free(kittri);      /* §933-И (NULL-safe: освобождён раньше в Emip-ветви) */
  if (rcv != NULL) { /* §921: приёмник */
    for (int rk = 0; rk < nrecv; rk++) {
      pg_recv_close(&rcv[rk]);
      free(rayD ? rayD[rk] : NULL);
    }
    free(rcv);
    free(rayL0);
    free(rayD);
    free(rayS);
    free(rayAD);
    free(rayMX);
  }
  free(so.lpacc);   /* §915-R4 (NULL-safe) */
  free(so.lpacc_e); /* §932-Б (NULL-safe) */
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
