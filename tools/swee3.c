/* swee3 — СВИП на новых массивах (§835-§837): итерации + аналитика полости
 * + трафик-прибор §5 STRUCTURE.md + режим сличения со старым путём (cmp=).
 *
 * Ключи: it=N (умолчание 20), rho=F (<0 — kd материалов; 0 — НК),
 * le=F (умолчание 1), tau0 (НК: слой не взаимодействует),
 * noprop (НК: луч не переносится), dirs=6|26 (§837), lev=N,
 * cmp=ФАЙЛ (сличение E с дампом старого пути по индексу куска).
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

static double now_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

int main(int argc, char **argv) {
  const char *path = NULL, *cmpfile = NULL;
  double scale = 1.0, le = 1.0, rho = -1.0;
  int iters = 20, lev = 6, tau0 = 0, noprop = 0, ndirs = 6, i, ax;
  double t0, t1;
  hz_objmesh m;
  hz_pyr py;
  hz_sw_opts so;
  hz_sw_stat st;
  double *area = NULL, *nrm = NULL, *kd = NULL, *cent = NULL, *cmin = NULL, *cmax = NULL;
  int32_t *mtl = NULL;
  double *hist = NULL;
  double cell;

  for (i = 1; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0)
      lev = atoi(argv[i] + 4);
    else if (strncmp(argv[i], "rho=", 4) == 0)
      rho = atof(argv[i] + 4);
    else if (strncmp(argv[i], "le=", 3) == 0)
      le = atof(argv[i] + 3);
    else if (strncmp(argv[i], "it=", 3) == 0)
      iters = atoi(argv[i] + 3);
    else if (strcmp(argv[i], "tau0") == 0)
      tau0 = 1;
    else if (strcmp(argv[i], "noprop") == 0)
      noprop = 1;
    else if (strncmp(argv[i], "dirs=", 5) == 0)
      ndirs = atoi(argv[i] + 5);
    else if (strncmp(argv[i], "cmp=", 4) == 0)
      cmpfile = argv[i] + 4;
    else if (strncmp(argv[i], "scale=", 6) == 0)
      scale = atof(argv[i] + 6);
    else
      path = argv[i];
  }
  if (!path) {
    fprintf(stderr, "use: swee3 <scene.obj> [it=N] [rho=F] [le=F] [tau0] [noprop] [dirs=6|26] "
                    "[cmp=ФАЙЛ] [lev=N]\n");
    return 2;
  }
  if (hz_obj_load(&m, path, scale) != 0) {
    fprintf(stderr, "swee3: не читается %s\n", path);
    return 2;
  }

  area = (double *)malloc((size_t)m.nt * sizeof *area);
  nrm = (double *)malloc((size_t)m.nt * 3 * sizeof *nrm);
  kd = (double *)malloc((size_t)m.nt * sizeof *kd);
  cent = (double *)malloc((size_t)m.nt * 3 * sizeof *cent);
  cmin = (double *)malloc((size_t)m.nt * 3 * sizeof *cmin);
  cmax = (double *)malloc((size_t)m.nt * 3 * sizeof *cmax);
  mtl = (int32_t *)malloc((size_t)m.nt * sizeof *mtl);
  hist = (double *)calloc((size_t)iters, sizeof *hist);
  if (!area || !nrm || !kd || !cent || !cmin || !cmax || !mtl || !hist) {
    fprintf(stderr, "swee3: нет памяти\n");
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
  if (hz_pyr_build(&py, m.nt, cmin, cmax, cent, mtl, m.lo, m.hi, cell) != 0) {
    fprintf(stderr, "swee3: пирамида не построилась\n");
    return 2;
  }

  printf("== swee3 %s: nt=%d клетка %.4g м, it=%d le=%.3g rho=%s dirs=%d%s%s\n", path, m.nt, cell,
         iters, le, rho < 0 ? "kd" : "ovr", ndirs, tau0 ? " tau0" : "", noprop ? " noprop" : "");

  memset(&so, 0, sizeof so);
  so.le = le;
  so.rho = rho;
  so.iters = iters;
  so.tau0 = tau0;
  so.noprop = noprop;
  so.ndirs = ndirs;
  for (int mode = 1; mode < 2; mode++) { /* scalar §835 — только по требованию */
    const char *mname = mode ? "col" : "scalar";
    so.mode = mode;
    for (i = 0; i < m.nt; i++)
      py.pcs[i].e = 0.0f; /* режимы с чистого поля */
    t0 = now_sec();
    if (hz_sw_run(&py, m.nt, area, nrm, kd, &so, &st, hist) != 0) {
      fprintf(stderr, "swee3: свип не прошёл\n");
      return 2;
    }
    t1 = now_sec();
    printf("[%s nd=%d] трафик: %.2f МБ/итерацию, %.3f с → %.2f ГБ/с эффективной\n", mname,
           mode ? so.ndirs : 6, (double)st.traffic / 1048576.0, t1 - t0,
           (t1 - t0) > 1e-9 ? (double)st.traffic * (double)iters / (t1 - t0) / 1e9 : 0.0);
    printf("[%s] итерации E_avg:", mname);
    for (i = 0; i < iters; i++)
      printf(" %.4f", hist[i]);
    printf("\n");
    if (iters >= 3) {
      double f = hist[iters - 1] - hist[iters - 2];
      double f0 = hist[iters - 2] - hist[iters - 3];
      printf("[%s] фактор ряда: %.4f (ожидалось ~ρ)\n", mname, fabs(f0) > 1e-300 ? f / f0 : 0.0);
    }
    printf("[%s] баланс: излучено %.4f, поглощено %.4f, рециркуляция %.4f, потеряно %.4f\n", mname,
           st.emitted, st.absorbed, st.recycled, st.lost);
    if (fabs(rho - 1.0) > 1e-12 && !tau0 && !noprop) {
      double rr = rho < 0 ? kd[0] * 0 + 0.5 : rho; /* kd сцены печатается выше */
      printf("[%s] аналитика(ρ=%.2f): 1ст=%.4f 2ст=%.4f; E_avg=%.4f → %.2f×1ст, %.2f×2ст\n", mname,
             rr, le / (1.0 - rr), 2.0 * le / (1.0 - rr), st.e_avg, st.e_avg / (le / (1.0 - rr)),
             st.e_avg / (2.0 * le / (1.0 - rr)));
    }
  }

  /* сличение со старым путём (§837-P): дамп E по индексу куска */
  if (cmpfile) {
    FILE *f = fopen(cmpfile, "rb");
    int32_t nold = 0;
    if (!f || fread(&nold, sizeof nold, 1, f) != 1 || nold != m.nt) {
      fprintf(stderr, "swee3: дамп %s не читается или nt не совпал (%d vs %d)\n", cmpfile, nold,
              m.nt);
      return 2;
    }
    {
      double *eold = (double *)malloc((size_t)nold * sizeof *eold);
      double rmed, rsum = 0, asum = 0;
      int32_t q;
      int in2 = 0;
      double *rats;
      if (!eold || fread(eold, sizeof *eold, (size_t)nold, f) != (size_t)nold) {
        fprintf(stderr, "swee3: дамп обрезан\n");
        return 2;
      }
      fclose(f);
      rats = (double *)malloc((size_t)nold * sizeof *rats);
      for (q = 0; q < nold; q++) {
        double en = py.pcs[q].e;
        double r = (eold[q] > 1e-9) ? en / eold[q] : 1.0;
        rats[q] = r;
        rsum += r * area[q];
        asum += area[q];
        if (r > 0.5 && r < 2.0) in2++;
      }
      /* медиана отношений */
      {
        int32_t u, v;
        for (u = 1; u < nold; u++) {
          double ru = rats[u];
          for (v = u; v > 0 && rats[v - 1] > ru; v--)
            rats[v] = rats[v - 1];
          rats[v] = ru;
        }
        rmed = rats[nold / 2];
      }
      printf("ПАРИТЕТ: медиана отношения %.4f; площадь-взвешенное %.4f; доля в ×2: %.1f %%\n", rmed,
             rsum / asum, 100.0 * (double)in2 / (double)nold);
      free(rats);
      free(eold);
    }
  }

  free(area);
  free(nrm);
  free(kd);
  free(cent);
  free(cmin);
  free(cmax);
  free(mtl);
  free(hist);
  hz_pyr_free(&py);
  hz_obj_free(&m);
  return 0;
}
