/* pref — ЭТАЛОН ПЕРЕНОСА (§850): итерационное решение на кусках СО
 * СЛУЧАЙНЫМИ направлениями и ТОЧНОЙ видимостью.
 *
 * ЗАЧЕМ ОТДЕЛЬНО ОТ СВИПА. Свип §845-846 решает перенос фронтом на
 * квадратуре: его зазор к непрерывной физике может быть и схемой, и
 * квадратурой. Эталон не наследует ни того, ни другого: направления
 * случайные (А1548 — ray effect не наследуется), видимость точная
 * (ray-triangle через DDA §849, побитово верифицированный против
 * перебора). Расхождение свип/эталон = ошибка СХЕМЫ ПЕРЕНОСА (А1549).
 *
 * ОЦЕНКА. Для куска p, обе полусферы (двусторонняя конвенция А1471):
 *
 *     E_p = (2π/K)·Σ_k |cos_k|·L_hit(k)
 *
 * направление k — равномерно по сфере; трассируем ЛУЧ ИЗ КУСКА вдоль k,
 * попадание q даёт входящий радианс L_out,q. Итерации Жакоби:
 * L_out = Le + ρ·E/2π. Промах мимо сцены — тёмное окружение; в замкнутой
 * сцене промахов нет (прибор замкнутости). MC-шум ~1/√K на кусок,
 * ~1/√(K·nt) на среднее (А1547).
 *
 * Запуск: pref <scene.obj> [it=N rho=F le=F dirs=6|26|NxM lev=N K=N]
 * ключи свипа те же; K — направлений на кусок (умолчание 1024).
 */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "nstruct/pyr.h"
#include "nstruct/sweep.h"
#include "scene_obj.h"

#define PREF_K_DEFAULT 1024
#define PREF_MAXIT 64
#define PREF_FACTOR_STOP 1e-3 /* ниже — ряд сошёлся в MC-шум, дальше не зачем */

static double now_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static double ref_ray_tri(const double o[3], const double d[3], const double p[3][3]) {
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
  if (fabs(det) < 1e-30) return -1.0;
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

/* ---- DDA (§849, тот же код; вынесен в pref для независимости запуска) ---- */

typedef struct {
  int64_t *start;
  int32_t *pids;
  int64_t ntotal;
} ref_bbox_csr;

static void ref_bbox_span(const hz_pyr *py, const double cmin[3], const double cmax[3],
                          int64_t i0[3], int64_t i1[3]) {
  int32_t ax;
  i0[0] = i0[1] = i0[2] = 0; /* анализатор теряет индукцию по ax (FP-класс diam) */
  i1[0] = i1[1] = i1[2] = 0;
  for (ax = 0; ax < 3; ax++) {
    int64_t n = ax == 0 ? py->nx : (ax == 1 ? py->ny : py->nz);
    i0[ax] = (int64_t)((cmin[ax] - py->lo[ax]) / py->cell);
    i1[ax] = (int64_t)((cmax[ax] - py->lo[ax]) / py->cell);
    if (i0[ax] < 0) i0[ax] = 0;
    if (i1[ax] >= n) i1[ax] = n - 1;
  }
}

static int64_t ref_cell_id(const hz_pyr *py, int64_t i, int64_t j, int64_t k) {
  return i + py->nx * (j + py->ny * k);
}

static void ref_bbox_csr_free(ref_bbox_csr *csr) {
  free(csr->start);
  free(csr->pids);
  memset(csr, 0, sizeof *csr);
}

static int ref_bbox_csr_build(const hz_pyr *py, const double *cmin, const double *cmax,
                              ref_bbox_csr *csr) {
  int64_t *cnt;
  int32_t p, li;
  memset(csr, 0, sizeof *csr);
  cnt = (int64_t *)calloc((size_t)py->nleaf, sizeof *cnt);
  csr->start = (int64_t *)calloc((size_t)py->nleaf + 1, sizeof *csr->start);
  if (!cnt || !csr->start) {
    free(cnt);
    ref_bbox_csr_free(csr);
    return 2;
  }
  for (p = 0; p < py->nt; p++) {
    int64_t i0[3], i1[3], i, j, k;
    int32_t tri = py->pcs[p].tri;
    ref_bbox_span(py, cmin + 3 * (int64_t)tri, cmax + 3 * (int64_t)tri, i0, i1);
    for (k = i0[2]; k <= i1[2]; k++)
      for (j = i0[1]; j <= i1[1]; j++)
        for (i = i0[0]; i <= i1[0]; i++) {
          int32_t pos = hz_pyr_leaf_pos(py, ref_cell_id(py, i, j, k));
          if (pos >= 0) cnt[pos]++;
        }
  }
  for (li = 0; li < py->nleaf; li++)
    csr->start[li + 1] = csr->start[li] + cnt[li];
  csr->ntotal = csr->start[py->nleaf];
  csr->pids = (int32_t *)malloc((size_t)(csr->ntotal > 0 ? csr->ntotal : 1) * sizeof *csr->pids);
  if (!csr->pids) {
    free(cnt);
    ref_bbox_csr_free(csr);
    return 2;
  }
  for (li = 0; li < py->nleaf; li++)
    cnt[li] = csr->start[li];
  for (p = 0; p < py->nt; p++) {
    int64_t i0[3], i1[3], i, j, k;
    int32_t tri = py->pcs[p].tri;
    ref_bbox_span(py, cmin + 3 * (int64_t)tri, cmax + 3 * (int64_t)tri, i0, i1);
    for (k = i0[2]; k <= i1[2]; k++)
      for (j = i0[1]; j <= i1[1]; j++)
        for (i = i0[0]; i <= i1[0]; i++) {
          int32_t pos = hz_pyr_leaf_pos(py, ref_cell_id(py, i, j, k));
          if (pos >= 0) csr->pids[cnt[pos]++] = p;
        }
  }
  free(cnt);
  return 0;
}

/* ближайший кусок вдоль луча; возврат куска или -1 */
static int32_t ref_trace(const hz_pyr *py, const ref_bbox_csr *csr, const hz_objmesh *m,
                         const double o[3], const double d[3]) {
  double tlo = 0.0, thi = 1e30, tcur, tbest = -1.0;
  int64_t cell[3], stepv[3];
  double tnext[3], tdelta[3];
  int32_t ax, best = -1;
  for (ax = 0; ax < 3; ax++) {
    double n = ax == 0 ? (double)py->nx : (ax == 1 ? (double)py->ny : (double)py->nz);
    if (fabs(d[ax]) < 1e-30) {
      if (o[ax] < py->lo[ax] || o[ax] > py->lo[ax] + n * py->cell) return -1;
      continue;
    }
    {
      double ta = (py->lo[ax] - o[ax]) / d[ax];
      double tb = (py->lo[ax] + n * py->cell - o[ax]) / d[ax];
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
    int64_t nax = ax == 0 ? py->nx : (ax == 1 ? py->ny : py->nz);
    double pos = o[ax] + tcur * d[ax];
    int64_t idx = (int64_t)((pos - py->lo[ax]) / py->cell);
    if (idx < 0) idx = 0;
    if (idx >= nax) idx = nax - 1;
    cell[ax] = idx;
    if (fabs(d[ax]) < 1e-30) {
      tnext[ax] = 1e30;
      tdelta[ax] = 0.0;
      stepv[ax] = 0;
    } else {
      double boundary = d[ax] > 0.0 ? py->lo[ax] + ((double)idx + 1.0) * py->cell
                                    : py->lo[ax] + (double)idx * py->cell;
      stepv[ax] = d[ax] > 0.0 ? 1 : -1;
      tdelta[ax] = py->cell / fabs(d[ax]);
      tnext[ax] = tcur + (boundary - pos) / d[ax];
    }
  }
  for (;;) {
    int32_t pos = hz_pyr_leaf_pos(py, ref_cell_id(py, cell[0], cell[1], cell[2]));
    double tn;
    if (pos >= 0) {
      int64_t s;
      for (s = csr->start[pos]; s < csr->start[pos + 1]; s++) {
        int32_t p = csr->pids[s];
        double p3[3][3], tt;
        hz_obj_tri(m, py->pcs[p].tri, p3);
        tt = ref_ray_tri(o, d, p3);
        if (tt >= 0.0 && (tbest < 0.0 || tt < tbest)) {
          tbest = tt;
          best = p;
        }
      }
    }
    tn = tnext[0];
    if (tnext[1] < tn) tn = tnext[1];
    if (tnext[2] < tn) tn = tnext[2];
    if (tn > thi || (tbest >= 0.0 && tn > tbest)) break;
    if (tnext[0] <= tnext[1] && tnext[0] <= tnext[2]) {
      tnext[0] += tdelta[0];
      cell[0] += stepv[0];
    } else if (tnext[1] <= tnext[2]) {
      tnext[1] += tdelta[1];
      cell[1] += stepv[1];
    } else {
      tnext[2] += tdelta[2];
      cell[2] += stepv[2];
    }
    if (cell[0] < 0 || cell[1] < 0 || cell[2] < 0 || cell[0] >= py->nx || cell[1] >= py->ny ||
        cell[2] >= py->nz)
      break;
  }
  return best;
}

/* псевдослучайные числа: xorshift64 — детерминизм прогона */
static uint64_t ref_rng = 0x9E3779B97F4A7C15ULL;
static double ref_urand(void) {
  ref_rng ^= ref_rng << 13;
  ref_rng ^= ref_rng >> 7;
  ref_rng ^= ref_rng << 17;
  return (double)(ref_rng >> 11) * (1.0 / 9007199254740992.0);
}

int main(int argc, char **argv) {
  const char *path = NULL;
  double scale = 1.0, le = 1.0, rho = -1.0;
  int iters = 30, lev = 6, mort = 1, i, ax;
  int K = PREF_K_DEFAULT;
  hz_objmesh m;
  hz_pyr py;
  hz_sw_opts so;
  hz_sw_stat st;
  double *area = NULL, *nrm = NULL, *kd = NULL, *cent = NULL, *cmin = NULL, *cmax = NULL;
  int32_t *mtl = NULL;
  double cell, *Eref = NULL, *Lout = NULL;
  ref_bbox_csr csr;

  for (i = 1; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0)
      lev = atoi(argv[i] + 4);
    else if (strncmp(argv[i], "rho=", 4) == 0)
      rho = atof(argv[i] + 4);
    else if (strncmp(argv[i], "le=", 3) == 0)
      le = atof(argv[i] + 3);
    else if (strncmp(argv[i], "it=", 3) == 0)
      iters = atoi(argv[i] + 3);
    else if (strncmp(argv[i], "K=", 2) == 0)
      K = atoi(argv[i] + 2);
    else if (strncmp(argv[i], "mort=", 5) == 0)
      mort = atoi(argv[i] + 5);
    else if (strncmp(argv[i], "scale=", 6) == 0)
      scale = atof(argv[i] + 6);
    else
      path = argv[i];
  }
  if (!path || K < 1) {
    fprintf(stderr, "use: pref <scene.obj> [it=N rho=F le=F K=N lev=N]\n");
    return 2;
  }
  if (hz_obj_load(&m, path, scale) != 0) {
    fprintf(stderr, "pref: не читается %s\n", path);
    return 2;
  }
  area = (double *)malloc((size_t)m.nt * sizeof *area);
  nrm = (double *)malloc((size_t)m.nt * 3 * sizeof *nrm);
  kd = (double *)malloc((size_t)m.nt * sizeof *kd);
  cent = (double *)malloc((size_t)m.nt * 3 * sizeof *cent);
  cmin = (double *)malloc((size_t)m.nt * 3 * sizeof *cmin);
  cmax = (double *)malloc((size_t)m.nt * 3 * sizeof *cmax);
  mtl = (int32_t *)malloc((size_t)m.nt * sizeof *mtl);
  Eref = (double *)malloc((size_t)m.nt * sizeof *Eref);
  Lout = (double *)malloc((size_t)m.nt * sizeof *Lout);
  if (!area || !nrm || !kd || !cent || !cmin || !cmax || !mtl || !Eref || !Lout) {
    fprintf(stderr, "pref: нет памяти\n");
    return 2;
  }
  for (i = 0; i < m.nt; i++) {
    double p[3][3], e1[3], e2[3], nn;
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
    mtl[i] = m.fm[i];
  }
  {
    double maxdim = 0.0;
    for (ax = 0; ax < 3; ax++) {
      double s = m.hi[ax] - m.lo[ax];
      if (s > maxdim) maxdim = s;
    }
    cell = maxdim / (double)(1 << lev);
  }
  /* §852: вершины/боксы/ℓ_p заполняются в блоке СВЕРКИ ниже (после memset so) */
  if (hz_pyr_build(&py, m.nt, cmin, cmax, cent, mtl, m.lo, m.hi, cell) != 0) {
    fprintf(stderr, "pref: пирамида не построилась\n");
    return 2;
  }
  if (mort && (hz_pyr_morton(&py) != 0 || hz_pyr_permute(&py, area, sizeof *area) != 0 ||
               hz_pyr_permute(&py, nrm, 3 * sizeof *nrm) != 0 ||
               hz_pyr_permute(&py, kd, sizeof *kd) != 0)) {
    fprintf(stderr, "pref: Morton не прошёл\n");
    return 2;
  }
  if (ref_bbox_csr_build(&py, cmin, cmax, &csr) != 0) {
    fprintf(stderr, "pref: bbox-CSR не построился\n");
    return 2;
  }

  double eavg_ref = 0.0; /* §850: живёт до СВЕРКИ */
  /* --- ЭТАЛОН: итерации Жакоби со случайными направлениями --- */
  {
    double rr = rho < 0 ? 0.5 : rho;
    double t0 = now_sec(), t1;
    double factor = 1.0, eprev = 0.0;
    int64_t miss_in = 0, miss_out = 0;
    int it;
    for (i = 0; i < m.nt; i++) {
      Lout[i] = le;
      Eref[i] = 0.0;
    }
    for (it = 0; it < iters; it++) {
      for (i = 0; i < m.nt; i++) {
        double acc = 0.0, u2[3], w2[3], nn2, org[3];
        double tmp[3] = {0, 0, 1};
        int k;
        /* базис полусферы куска — ОДИН раз на кусок */
        if (fabs(nrm[3 * (int64_t)i + 2]) > 0.9) {
          tmp[0] = 1;
          tmp[1] = tmp[2] = 0;
        }
        u2[0] = nrm[3 * (int64_t)i + 1] * tmp[2] - nrm[3 * (int64_t)i + 2] * tmp[1];
        u2[1] = nrm[3 * (int64_t)i + 2] * tmp[0] - nrm[3 * (int64_t)i + 0] * tmp[2];
        u2[2] = nrm[3 * (int64_t)i + 0] * tmp[1] - nrm[3 * (int64_t)i + 1] * tmp[0];
        nn2 = sqrt(u2[0] * u2[0] + u2[1] * u2[1] + u2[2] * u2[2]);
        u2[0] /= nn2;
        u2[1] /= nn2;
        u2[2] /= nn2;
        w2[0] = nrm[3 * (int64_t)i + 1] * u2[2] - nrm[3 * (int64_t)i + 2] * u2[1];
        w2[1] = nrm[3 * (int64_t)i + 2] * u2[0] - nrm[3 * (int64_t)i + 0] * u2[2];
        w2[2] = nrm[3 * (int64_t)i + 0] * u2[1] - nrm[3 * (int64_t)i + 1] * u2[0];
        org[0] = cent[3 * (int64_t)i];
        org[1] = cent[3 * (int64_t)i + 1];
        org[2] = cent[3 * (int64_t)i + 2];
        for (k = 0; k < K; k++) {
          double phi = 2.0 * M_PI * ref_urand();
          double mu = ref_urand(), sq = sqrt(1.0 - mu * mu);
          double s1 = sin(phi) * sq, c1 = cos(phi) * sq;
          double side = (k & 1) ? -1.0 : 1.0; /* обе полусферы (А1471) */
          double om[3];
          int32_t hit;
          om[0] = side * (mu * nrm[3 * (int64_t)i] + s1 * u2[0] + c1 * w2[0]);
          om[1] = side * (mu * nrm[3 * (int64_t)i + 1] + s1 * u2[1] + c1 * w2[1]);
          om[2] = side * (mu * nrm[3 * (int64_t)i + 2] + s1 * u2[2] + c1 * w2[2]);
          {
            double o2[3] = {org[0] + 1e-5 * om[0], org[1] + 1e-5 * om[1], org[2] + 1e-5 * om[2]};
            hit = ref_trace(&py, &csr, &m, o2, om);
          }
          if (hit >= 0)
            acc += mu * Lout[hit]; /* |cos|·L; вес 4π/K — K на ОБЕ полусферы */
          else if (side > 0)
            miss_in++;
          else
            miss_out++;
        }
        Eref[i] = (4.0 * M_PI / (double)K) * acc;
      }
      {
        double emean = 0, asum = 0;
        for (i = 0; i < m.nt; i++) {
          emean += Eref[i] * area[i];
          asum += area[i];
        }
        emean /= asum;
        if (it >= 1 && eprev > 1e-300) factor = fabs(emean - eprev) / eprev;
        eprev = emean;
        eavg_ref = emean;
        for (i = 0; i < m.nt; i++)
          Lout[i] = le + rr * Eref[i] / (2.0 * M_PI);
      }
      t1 = now_sec();
      printf("ЭТАЛОН it=%d: E_avg=%.4f (фактор %.2e, %.1f с)\n", it, eavg_ref, factor, t1 - t0);
      if (factor < PREF_FACTOR_STOP && it >= 2) break;
    }
    t1 = now_sec();
    printf("ЭТАЛОН ИТОГ: E_avg=%.4f, k непрерывной(замкн.)=%.4f, промахов внутр/внеш %" PRId64
           "/%" PRId64 ", %.1f с\n",
           eavg_ref, eavg_ref / (2.0 * M_PI * le / (1.0 - rr)), miss_in, miss_out, t1 - t0);
  }

  /* --- СВИП на том же носителе и СВЕРКА с эталоном ---
   * mode=3 (§852): полный набор точного пересечения — trivert/tribox/lp/domhi
   * + per-node max ℓ_p; раньше здесь был mode=2+vc без trivert (упал бы и со
   * старым порядком: memset ниже стирал заполненный до него so.trivert). */
  memset(&so, 0, sizeof so);
  so.le = le;
  so.rho = rho;
  so.iters = iters;
  so.ndirs = 26;
  so.mode = 3;
  so.build = 1;
  so.vc = 1;
  {
    /* §852: trivert/tribox в порядке ИСХОДНЫХ треугольников (pcs[].tri) */
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
      int32_t aa, a2;
      hz_obj_tri(&m, ti, pp);
      for (aa = 0; aa < 9; aa++)
        tv9[9 * (int64_t)ti + aa] = ((const double *)pp)[aa];
      for (a2 = 0; a2 < 3; a2++) {
        double lo = pp[0][a2], hi = pp[0][a2];
        int v;
        for (v = 1; v < 3; v++) {
          if (pp[v][a2] < lo) lo = pp[v][a2];
          if (pp[v][a2] > hi) hi = pp[v][a2];
        }
        tb6[6 * (int64_t)ti + a2] = lo;
        tb6[6 * (int64_t)ti + 3 + a2] = hi;
      }
      lparr[ti] = 0; /* ℓ_p = 0: сверка на рабочем уровне листьев */
    }
    so.trivert = tv9;
    so.tribox = tb6;
    so.lp = lparr;
    so.domhi = m.hi;
  }
  {
    int32_t nup = 0;
    if (hz_pyr_set_lp(&py, (const uint8_t *)so.lp, &nup) != 0) {
      fprintf(stderr, "pref: per-node max lp не построился\n");
      return 2;
    }
  }
  for (i = 0; i < m.nt; i++)
    py.pcs[i].e = 0.0f;
  if (hz_sw_run(&py, m.nt, area, nrm, kd, &so, &st, NULL) != 0) {
    fprintf(stderr, "pref: свип не прошёл\n");
    return 2;
  }
  {
    double esum = 0, asum = 0, dsum = 0, dmax = 0;
    int i2;
    double *rat = (double *)malloc((size_t)m.nt * sizeof *rat);
    double eavg_swp, rmed;
    int32_t u, v;
    if (!rat) return 2;
    for (i2 = 0; i2 < m.nt; i2++) {
      double e = (double)py.pcs[i2].e;
      double d = fabs(e - Eref[i2]);
      esum += e * area[i2];
      asum += area[i2];
      dsum += d * area[i2];
      if (d > dmax) dmax = d;
      rat[i2] = Eref[i2] > 1e-9 ? e / Eref[i2] : 1.0;
    }
    eavg_swp = esum / asum;
    for (u = 1; u < m.nt; u++) {
      double ru = rat[u];
      for (v = u; v > 0 && rat[v - 1] > ru; v--)
        rat[v] = rat[v - 1];
      rat[v] = ru;
    }
    rmed = rat[m.nt / 2];
    printf("СВЕРКА: E_avg свипа=%.4f, эталона=%.4f → отношение %.4f\n", eavg_swp, eavg_ref,
           eavg_swp / eavg_ref);
    printf("СВЕРКА: |Δ| средняя=%.4f, max=%.4f (по кускам); медиана отношения куска=%.4f\n",
           dsum / asum, dmax, rmed);
    free(rat);
  }
  free(Eref);
  free(Lout);
  free(area);
  free(nrm);
  free(kd);
  free(cent);
  free(cmin);
  free(cmax);
  free(mtl);
  ref_bbox_csr_free(&csr);
  hz_pyr_free(&py);
  hz_obj_free(&m);
  return 0;
}
