/* pkernel — Ш3: ЗАМЕР ВНУТРЕННЕГО ЦИКЛА (PLAN_ELEMENTS.md).
 *
 * Что меряется: время и ТРАФИК на (направление × полигон); доля растеризации
 * против доли редукции; обе раскладки фрагментов; масштабирование по потокам.
 *
 * Запуск: build/pkernel [delta] [h ...]
 */

#include "poly_seg.h"
#include "polygon.h"
#include "psweep.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include "transport/dirs3.h"
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

/* Восемь потолочных светильников §2: сетка 2×4 при y = 1.95, квадрат 0.3×0.3,
 * нормаль вниз, L_e = 1. Сетка натягивается на габарит зала с отступом в
 * полшага, чтобы светильник не сел на стену. */
static int add_lamps(hz_polyset *ps, const hz_objmesh *m, int32_t *first) {
  const double n[3] = {0.0, -1.0, 0.0}, eu[3] = {1.0, 0.0, 0.0};
  *first = ps->np;
  for (int i = 0; i < HZ_CFG_LAMP_NX; i++)
    for (int j = 0; j < HZ_CFG_LAMP_NZ; j++) {
      double c[3];
      c[0] = m->lo[0] + ((double)i + 0.5) / HZ_CFG_LAMP_NX * (m->hi[0] - m->lo[0]);
      c[1] = HZ_CFG_LAMP_Y;
      c[2] = m->lo[2] + ((double)j + 0.5) / HZ_CFG_LAMP_NZ * (m->hi[2] - m->lo[2]);
      int rc = hz_poly_add_quad(ps, c, n, eu, HZ_CFG_LAMP_SIDE / 2.0, HZ_CFG_LAMP_SIDE / 2.0, 0);
      if (rc != 0) return rc;
    }
  return 0;
}

int main(int argc, char **argv) {
  double delta = (argc > 1) ? strtod(argv[1], NULL) : 0.045;

  hz_objmesh m;
  if (hz_obj_load(&m, HZ_CFG_HALL_OBJ, HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет %s — scripts/fetch_scene.sh assets\n", HZ_CFG_HALL_OBJ);
    return 1;
  }
  hz_pseglist sg;
  if (hz_seg_planar(&sg, &m, delta) != 0) return 1;
  hz_polyset ps;
  if (hz_poly_build(&ps, &m, &sg) != 0) return 1;
  int32_t lamp0 = 0;
  if (add_lamps(&ps, &m, &lamp0) != 0) return 1;

  hz_ptrans t;
  if (hz_ptrans_init(&t, &ps, &m) != 0) return 1;
  for (int32_t k = lamp0; k < ps.np; k++) {
    t.Le[k] = HZ_CFG_LAMP_LE;
    t.rho[k] = 0.0; /* светильник не переотражает: он источник, а не зеркало */
  }
  double area = 0.0, arho = 0.0;
  for (int32_t k = 0; k < ps.np; k++) {
    area += ps.p[k].area;
    arho += ps.p[k].area * t.rho[k];
  }
  printf("== Ш3: зал, δ = %g м; полигонов %d (из них светильников %d), площадь %.1f м²,\n"
         "   средневзвешенное альбедо %.3f\n",
         delta, ps.np, ps.np - lamp0, area, arho / area);

  /* --- ПОДПИСЬ ПОКРЫТИЯ: Σ h²/|n·ω| по фрагментам обязана дать площадь ---
   * Проверка растеризатора, а не переноса: она ловит и дыры, и двойное
   * покрытие, и неверный край — по числу, которое известно независимо. */
  {
    double w[3] = {0.3, -0.8, 0.5}, lo[3], hi[3];
    for (int a = 0; a < 3; a++) {
      lo[a] = m.lo[a];
      hi[a] = m.hi[a];
    }
    double hh = 0.01;
    hz_pview v;
    if (hz_pview_make(&v, w, lo, hi, hh) == 0) {
      hz_span *sp = NULL;
      int64_t nsp = 0, cap = 0;
      hz_rstats rs;
      if (hz_prast_spans(&sp, &nsp, &cap, &v, &ps, &rs) == 0) {
        double cov = 0.0;
        for (int64_t s = 0; s < nsp; s++) {
          const hz_poly *P = &ps.p[sp[s].poly];
          double nw = P->n[0] * v.w[0] + P->n[1] * v.w[1] + P->n[2] * v.w[2];
          cov += (double)(sp[s].i1 - sp[s].i0 + 1) * hh * hh / fabs(nw);
        }
        printf("   подпись покрытия при h = %g: Σ h²/|n·ω| = %.2f м² против площади %.2f м² "
               "(отн. %.4f); экран %d×%d, фрагментов %lld, полос %lld\n",
               hh, cov, area, cov / area - 1.0, v.W, v.H, (long long)rs.nfrag,
               (long long)rs.nspan);
      }
      free(sp);
    }
  }

  /* --- свип по шагу растра и раскладке --- */
  tr3_dirs d;
  if (tr3_dirs_product(&d, 2, 2) != 0) return 1; /* ND = 32 */
  printf("\n   ND = %d, один отскок\n", d.n);
  printf("   %-6s %-6s %5s  %9s %9s %9s %9s %9s   %8s %8s %8s\n", "h, м", "раскл", "потк",
         "полосы", "фрагм", "растр,с", "сорт,с", "редук,с", "всего,с", "ГБ/с", "нс/фр");

  for (int i = 2; i < argc; i++) {
    double h = strtod(argv[i], NULL);
    for (int lay = 0; lay < 2; lay++)
      for (int nthr = 1; nthr <= 16; nthr *= 4) {
        memset(t.E, 0, (size_t)t.np * 3 * sizeof *t.E);
        hz_pstats st;
        double ta = now_s();
        int rc = hz_psweep_bounce(&t, &d, h, nthr, 0.0, lay, &st);
        double tb = now_s();
        if (rc != 0) {
          printf("   %-6g %-6s %5d  ОТКАЗ rc=%d\n", h, lay ? "списк" : "пробг", nthr, rc);
          continue;
        }
        double tot = st.t_raster + st.t_sort + st.t_reduce;
        printf("   %-6g %-6s %5d  %9lld %9lld %9.3f %9.3f %9.3f   %8.3f %8.2f %8.1f\n", h,
               lay ? "списк" : "пробг", nthr, (long long)st.nspan, (long long)st.nfrag,
               st.t_raster, st.t_sort, st.t_reduce, tb - ta,
               (double)st.nbytes / (tot > 0.0 ? tot : 1.0) / 1e9,
               1e9 * (tb - ta) * nthr / (double)(st.nfrag > 0 ? st.nfrag : 1));
        if (st.nover != 0)
          printf("        !! ёмкость A-буфера мала: %lld фрагментов потеряно\n",
                 (long long)st.nover);
      }
  }

  /* --- цена на (направление × полигон) --- */
  if (argc > 2) {
    double h = strtod(argv[2], NULL);
    memset(t.E, 0, (size_t)t.np * 3 * sizeof *t.E);
    hz_pstats st;
    double ta = now_s();
    hz_psweep_bounce(&t, &d, h, 1, 0.0, HZ_LAYOUT_RUNS, &st);
    double tb = now_s();
    printf("\n   при h = %g, ND = %d, один поток: %.3f с на отскок,\n"
           "   %.1f нс на (направление × полигон), %.2f фрагмента на (направление × полигон),\n"
           "   решение %.4f с (%.2f%% отскока), ограничитель тронул %lld полигонов,\n"
           "   Φ_in = %.6g, Φ_out = %.6g\n",
           h, d.n, tb - ta, 1e9 * (tb - ta) / ((double)d.n * ps.np),
           (double)st.nfrag / ((double)d.n * ps.np), st.t_solve, 100.0 * st.t_solve / (tb - ta),
           (long long)st.nclip, st.phi_in, st.phi_out);
  }

  tr3_dirs_free(&d);
  hz_ptrans_free(&t);
  hz_poly_free(&ps);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
