/* pcuts — Ш6, часть первая: СКОЛЬКО КЛИНЬЕВ ПРОХОДИТ ПРАВИЛА A–D.
 *
 * Отделено от `prender` намеренно. Вопрос «во что обходятся теневые разрезы»
 * решается ДО картинки и без неё: если правило A отвергает почти всё, то
 * измерять качество картинки уже не нужно — нечего сравнивать. Ожидание,
 * записанное до прогона: для светильника `0.3×0.3` м угловой размер с метра
 * `α ≈ 0.3`, длина клина `r_A = ε/α ≈ 15` см; с двух метров `≈30` см. Значит
 * клинья обязаны жить лоскутами у КОНТАКТОВ (ножка стула у пола, кромка
 * столешницы) и больше нигде.
 */

#include "pcut.h"
#include "pdirect.h"
#include "poly_seg.h"
#include "polygon.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

int main(int argc, char **argv) {
  double delta = (argc > 1) ? strtod(argv[1], NULL) : 0.045;
  double eps = (argc > 2) ? strtod(argv[2], NULL) : 0.045;
  double tol = (argc > 3) ? strtod(argv[3], NULL) : 0.05;
  int nsrc = (argc > 4) ? atoi(argv[4]) : 8;

  hz_objmesh m;
  if (hz_obj_load(&m, HZ_CFG_HALL_OBJ, HZ_CFG_HALL_SCALE) != 0) return 1;
  hz_pseglist sg;
  double t0 = now_s();
  if (hz_seg_planar(&sg, &m, delta) != 0) return 1;
  hz_polyset ps;
  if (hz_poly_build(&ps, &m, &sg) != 0) return 1;
  int32_t src[8];
  const double n[3] = {0.0, -1.0, 0.0}, eu[3] = {1.0, 0.0, 0.0};
  int k = 0;
  for (int i = 0; i < HZ_CFG_LAMP_NX && k < nsrc; i++)
    for (int j = 0; j < HZ_CFG_LAMP_NZ && k < nsrc; j++) {
      double c[3] = {m.lo[0] + ((double)i + 0.5) / HZ_CFG_LAMP_NX * (m.hi[0] - m.lo[0]),
                     HZ_CFG_LAMP_Y,
                     m.lo[2] + ((double)j + 0.5) / HZ_CFG_LAMP_NZ * (m.hi[2] - m.lo[2])};
      src[k] = ps.np;
      if (hz_poly_add_quad(&ps, c, n, eu, HZ_CFG_LAMP_SIDE / 2.0, HZ_CFG_LAMP_SIDE / 2.0, 0) != 0)
        return 1;
      k++;
    }
  double *Le = calloc((size_t)ps.np, sizeof *Le);
  if (Le == NULL) return 1;
  for (int i = 0; i < nsrc; i++)
    Le[src[i]] = HZ_CFG_LAMP_LE;
  double t1 = now_s();
  printf("== Ш6/клинья: зал, δ = %g м, ε = %g м, tol = %g, источников %d\n", delta, eps, tol, nsrc);
  printf("   полигонов %d, сборка %.2f с\n", ps.np, t1 - t0);

  hz_cutcfg cfg = {eps, tol, 1 << 30};
  hz_wedges ws;
  double t2 = now_s();
  if (hz_wedges_build(&ws, &m, &ps, src, nsrc, Le, &cfg) != 0) return 1;
  double t3 = now_s();
  printf("   силуэтных рёбер (ребро × источник): %lld\n", (long long)ws.nsil);
  printf("   отсеяно правилом B (контраст ≤ %g): %lld\n", tol, (long long)ws.nreject_b);
  printf("   ПРОШЛО правила: %d  (%.4f%% силуэтных), построение %.2f с\n", ws.n,
         100.0 * ws.n / (double)(ws.nsil > 0 ? ws.nsil : 1), t3 - t2);
  if (ws.n > 0) {
    double *r = malloc((size_t)ws.n * sizeof *r);
    if (r == NULL) return 1;
    for (int32_t i = 0; i < ws.n; i++)
      r[i] = ws.w[i].r;
    qsort(r, (size_t)ws.n, sizeof *r, cmp_d);
    printf("   ДЛИНА КЛИНА (правило A и C), м: p10 %.4g  p50 %.4g  p90 %.4g  max %.4g\n",
           r[ws.n / 10], r[ws.n / 2], r[ws.n * 9 / 10], r[ws.n - 1]);
    printf("   приоритет ΔΦ/Φ: max %.4g, медиана %.4g\n", ws.w[0].prio, ws.w[ws.n / 2].prio);
    free(r);
  }

  /* Во что обходится измельчение при разных бюджетах (Т2). */
  printf("\n   %-10s %10s %10s %12s %10s\n", "бюджет", "треуг.", "полигонов", "измельч., с",
         "сборка, с");
  const int32_t budg[5] = {0, 64, 256, 1024, 4096};
  for (int b = 0; b < 5; b++) {
    if (budg[b] > ws.n && b > 0 && budg[b - 1] > ws.n) break;
    hz_objmesh mo;
    hz_pseglist so;
    double ta = now_s();
    if (hz_cut_apply(&mo, &so, &m, &sg, &ws, budg[b]) != 0) {
      printf("   %-10d ОТКАЗ (нет памяти)\n", budg[b]);
      break;
    }
    double tb = now_s();
    hz_polyset po;
    if (hz_poly_build(&po, &mo, &so) != 0) return 1;
    double tc = now_s();
    printf("   %-10d %10d %10d %12.2f %10.2f\n", budg[b], mo.nt, po.np, tb - ta, tc - tb);
    hz_poly_free(&po);
    hz_seg_free(&so);
    hz_obj_free(&mo);
  }

  free(Le);
  hz_wedges_free(&ws);
  hz_poly_free(&ps);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
