/* swee3 — СВИП №1 на новых массивах (§835): итерации + аналитика полости.
 *
 * Ключи: it=N (умолчание 20), rho=F (<0 — kd материалов; 0 — НК),
 * le=F (умолчание 1), tau0 (НК: слой не взаимодействует),
 * noprop (НК: луч не переносится), lev=N (умолчание 6).
 *
 * Прогон: cavity05 — замкнутая оболочка, равновесие ряда Неймана против
 * аналитики Le/(1−ρ) (А1471: сравнение — ИЗМЕРЕНИЕ угловой погрешности
 * 6-осевой квадратуры, не допуск; допуск шага — фактор ряда и НК).
 */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "nstruct/pyr.h"
#include "nstruct/sweep.h"
#include "scene_obj.h"

int main(int argc, char **argv) {
  const char *path = NULL;
  double scale = 1.0, le = 1.0, rho = -1.0;
  int iters = 20, lev = 6, tau0 = 0, noprop = 0, i, ax;
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
    else if (strncmp(argv[i], "scale=", 6) == 0)
      scale = atof(argv[i] + 6);
    else
      path = argv[i];
  }
  if (!path) {
    fprintf(stderr, "use: swee3 <scene.obj> [it=N] [rho=F] [le=F] [tau0] [noprop] [lev=N]\n");
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

  {
    double sa = 0, sb = 0;
    for (i = 0; i < m.nt; i++) {
      sa += area[i];
      sb += hz_obj_tri_area(&m, (int32_t)i);
    }
    {
      double dmax = 0;
      int imax = 0;
      for (i = 0; i < m.nt; i++) {
        double dd = fabs(area[i] - hz_obj_tri_area(&m, (int32_t)i));
        if (dd > dmax) {
          dmax = dd;
          imax = i;
        }
      }
      printf("Σarea моя=%.4f загрузчик=%.4f; maxΔ=%.4f piece %d моя=%.5f загр=%.5f\n", sa, sb, dmax,
             imax, area[imax], hz_obj_tri_area(&m, (int32_t)imax));
    }
  }
  printf("== swee3 %s: nt=%d клетка %.4g м, it=%d le=%.3g rho=%s%s%s\n", path, m.nt, cell, iters,
         le, rho < 0 ? "kd" : "ovr", tau0 ? " tau0" : "", noprop ? " noprop" : "");
  if (rho >= 0)
    printf("   (rho-переопределение: %.3g)\n", rho);
  else {
    /* сводка kd материалов — что реально в сцене */
    double s = 0;
    for (i = 0; i < m.nt; i++)
      s += kd[i] * area[i];
    printf("   kd площади-взвешенное: %.4f\n", s);
    if (rho < 0 && le > 0 && fabs(s) < 1e-12) printf("   ВНИМАНИЕ: kd=0 — поле не родится\n");
  }

  memset(&so, 0, sizeof so);
  so.le = le;
  so.rho = rho;
  so.iters = iters;
  so.tau0 = tau0;
  so.noprop = noprop;
  for (int mode = 0; mode < 2; mode++) {
    const char *mname = mode ? "col(§836)" : "scalar(§835)";
    so.mode = mode;
    for (i = 0; i < m.nt; i++)
      py.pcs[i].e = 0.0f; /* режимы с чистого поля */
    if (hz_sw_run(&py, m.nt, area, nrm, kd, &so, &st, hist) != 0) {
      fprintf(stderr, "swee3: свип не прошёл\n");
      return 2;
    }
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
      double rr = rho < 0 ? 0.5 : rho;
      printf("[%s] аналитика: 1ст=%.4f 2ст=%.4f; E_avg=%.4f → %.2f×1ст, %.2f×2ст\n", mname,
             le / (1.0 - rr), 2.0 * le / (1.0 - rr), st.e_avg, st.e_avg / (le / (1.0 - rr)),
             st.e_avg / (2.0 * le / (1.0 - rr)));
    }
  }

  printf("итерации E_avg:");
  for (i = 0; i < iters; i++)
    printf(" %.4f", hist[i]);
  printf("\n");
  if (iters >= 3) {
    double f = (hist[iters - 1] - hist[iters - 2]);
    double f0 = (hist[iters - 2] - hist[iters - 3]);
    printf("фактор ряда ΔE_k/ΔE_(k-1): %.4f (ожидалось ~ρ)\n", fabs(f0) > 1e-300 ? f / f0 : 0.0);
  }
  printf("баланс: излучено %.4f, поглощено %.4f, рециркуляция %.4f, потеряно %.4f (ед. потока)\n",
         st.emitted, st.absorbed, st.recycled, st.lost);
  if (fabs(rho - 1.0) > 1e-12 && !tau0 && !noprop) {
    printf("аналитика полости: Le/(1−ρ)·1ст=%.4f ·2ст=%.4f; E_avg=%.4f → %.2f×1ст, %.2f×2ст\n",
           le / (1.0 - (rho < 0 ? 0.5 : rho)), 2.0 * le / (1.0 - (rho < 0 ? 0.5 : rho)), st.e_avg,
           st.e_avg / (le / (1.0 - (rho < 0 ? 0.5 : rho))),
           st.e_avg / (2.0 * le / (1.0 - (rho < 0 ? 0.5 : rho))));
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
