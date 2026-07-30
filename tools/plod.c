/* plod — ИЕРАРХИЧЕСКОЕ ОГРУБЛЕНИЕ И СРЕЗ LOD (PLAN_ELEMENTS.md, §56, О16).
 *
 * Что печатается и почему именно это:
 *   - ПО УРОВНЯМ: число элементов и отношение к предыдущему. Отклонение от
 *     четырёх — ЗАМЕР, а не отклонение от нормы (§54.2): оно показывает, где
 *     сцена перестаёт быть гладкой на этом масштабе. Плюс `dmax` (обязан быть
 *     под допуском уровня) и форма `P/√A` (§55, «близость к правильному»);
 *   - ПО СРЕЗУ при камере §2: число элементов и метрика О20 через ОБЩИЙ набор
 *     лучей против нулевого уровня. Ставка шага (§56, предсказание 3): срез даёт
 *     в 2…10 раз меньше элементов, чем ЛУЧШИЙ однородный уровень при равной
 *     ошибке — поэтому рядом печатаются и однородные уровни, и срез;
 *   - СЧЁТЧИКИ ПОСТРОЙКИ: отказы ворот, тождественные продвижения, прирост
 *     площади (проверка §54.2), пары, смежные по ребру.
 *
 * ВЛОЖЕННОСТЬ ПРОВЕРЯЕТСЯ ЦЕЛОЧИСЛЕННО (А132): у каждого узла множество
 * исходных полигонов обязано быть объединением множеств детей без пересечений.
 * Проверка по МЕТКАМ, а не по площадям: «сходится по площади» даёт ложное
 * подтверждение при перекрытии.
 */

#include "pedge.h"
#include "plod.h"
#include "poly_seg.h"
#include "polygon.h"
#include "pray.h"
#include "pvfit.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include "transport/cam3.h"
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

static int cmp_dbl(const void *x, const void *y) {
  double a = *(const double *)x, b = *(const double *)y;
  return (a < b) ? -1 : ((a > b) ? 1 : 0);
}

static double pct(double *v, int64_t n, double p) {
  if (n <= 0) return 0.0;
  int64_t k = (int64_t)(p * (double)(n - 1));
  if (k < 0) k = 0;
  if (k >= n) k = n - 1;
  return v[k];
}

/* Метрика О20 через общий набор лучей: две сцены, разница по ЛУЧУ. */
static void metric(const hz_polyset *pf, const hz_polyset *pc, const tr3_camera *cam, double eps_px,
                   const char *tag) {
  hz_pray gf, gc;
  if (hz_pray_build(&gf, pf, 4.0) != 0) return;
  if (hz_pray_build(&gc, pc, 4.0) != 0) {
    hz_pray_free(&gf);
    return;
  }
  const int W = cam->w, H = cam->h;
  int64_t nb = 0, nl = 0, ng = 0;
  double *dt = malloc((size_t)W * (size_t)H * sizeof *dt);
  double *dp = malloc((size_t)W * (size_t)H * sizeof *dp);
  if (dt == NULL || dp == NULL) {
    free(dt);
    free(dp);
    hz_pray_free(&gf);
    hz_pray_free(&gc);
    return;
  }
  for (int j = 0; j < H; j++)
    for (int i = 0; i < W; i++) {
      double o[3], d[3], tf = 0.0, tc = 0.0;
      tr3_camera_ray(cam, i, j, o, d);
      int32_t hf = hz_pray_hit(&gf, o, d, 1e-6, &tf);
      int32_t hc = hz_pray_hit(&gc, o, d, 1e-6, &tc);
      if (hf < 0 && hc < 0) continue;
      if (hf >= 0 && hc < 0) {
        nl++;
        continue;
      }
      if (hf < 0) {
        ng++;
        continue;
      }
      double e = fabs(tf - tc);
      dt[nb] = e;
      const double *nf = pf->p[hf].n;
      double cr = fabs(nf[0] * d[0] + nf[1] * d[1] + nf[2] * d[2]);
      dp[nb] = e * cr / (eps_px * (tf > 0.0 ? tf : 1e-9));
      nb++;
    }
  qsort(dt, (size_t)nb, sizeof *dt, cmp_dbl);
  qsort(dp, (size_t)nb, sizeof *dp, cmp_dbl);
  int64_t tot = nb + nl + ng;
  printf("      %s: попали в обеих %.2f%%, потеряли %.2f%%, приобрели %.2f%%; |Δt| p50 %.4f p90 "
         "%.4f м; пиксели p50 %.2f p90 %.2f\n",
         tag, 100.0 * (double)nb / (double)(tot ? tot : 1),
         100.0 * (double)nl / (double)(tot ? tot : 1), 100.0 * (double)ng / (double)(tot ? tot : 1),
         pct(dt, nb, 0.5), pct(dt, nb, 0.9), pct(dp, nb, 0.5), pct(dp, nb, 0.9));
  free(dt);
  free(dp);
  hz_pray_free(&gf);
  hz_pray_free(&gc);
}

int main(int argc, char **argv) {
  int city = (argc > 1 && strcmp(argv[1], "city") == 0);
  double dseg = city ? 0.05 : 0.045;
  int simp = 0, vfit = 0, maxlev = 9, uniform = 0, viamerge = 0;
  double epsmul = 1.0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "simp") == 0) simp = 1;
    if (strcmp(argv[i], "vfit") == 0) vfit = 1;
    if (strcmp(argv[i], "uniform") == 0) uniform = 1;
    /* Второй построитель: уровень строится `hz_merge` над предыдущим (§58). */
    if (strcmp(argv[i], "viamerge") == 0) viamerge = 1;
    if (strncmp(argv[i], "lev=", 4) == 0) maxlev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "eps=", 4) == 0) epsmul = strtod(argv[i] + 4, NULL);
  }
  hz_objmesh m;
  if (hz_obj_load(&m, city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ,
                  city ? HZ_CFG_CITY_SCALE : HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  hz_pseglist sg;
  if (hz_seg_planar(&sg, &m, dseg) != 0) return 1;
  hz_polyset psf;
  if (hz_poly_build(&psf, &m, &sg) != 0) return 1;
  if (simp) {
    hz_edgestat es;
    if (hz_edge_simplify(&psf, 0.25 * dseg, HZ_EDGE_SHARED, &es) != 0) return 1;
  }
  tr3_camera cam;
  double eyeh[3] = HZ_CFG_HALL_EYE, ath[3] = HZ_CFG_HALL_AT;
  double eyec[3] = HZ_CFG_CITY_EYE, atc[3] = HZ_CFG_CITY_AT, up[3] = HZ_CFG_UP;
  const int W = 512, H = 512;
  if (tr3_camera_look(&cam, city ? eyec : eyeh, city ? atc : ath, up, HZ_CFG_FOV_DEG * M_PI / 180.0,
                      W, H) != 0)
    return 1;
  const double eps_px = (HZ_CFG_FOV_DEG * M_PI / 180.0) / (double)H;
  const double eps = eps_px * epsmul;

  double t0 = now_s();
  hz_lod L;
  int rc = viamerge ? hz_lod_build_merge(&L, &m, &sg, &psf, dseg, maxlev, eps)
                    : hz_lod_build(&L, &m, &sg, &psf, dseg, maxlev, eps);
  if (rc != 0) {
    fprintf(stderr, "отказ построения лестницы\n");
    return 1;
  }
  double t1 = now_s();
  printf("== ЛЕСТНИЦА: %s, δ0 %g м, уровней %d, узлов %d, построена за %.1f с\n",
         city ? "ГОРОД" : "зал", dseg, L.nlev, L.nnd, t1 - t0);
  printf("   постройка: пар по ребру %lld, отказов ворот %lld, тождественных продвижений %lld, "
         "худший прирост площади %.3f\n",
         (long long)L.npair_edge, (long long)L.ngate_rej, (long long)L.nident, L.area_grow);
  fflush(stdout);

  /* --- по уровням --- */
  int32_t prevn = 0;
  for (int lev = 0; lev < L.nlev; lev++) {
    const int32_t *lab = L.lab + (size_t)lev * (size_t)L.np;
    int32_t *seen = calloc((size_t)L.nnd, sizeof *seen);
    if (seen == NULL) break;
    int32_t cnt = 0;
    double *dm = malloc((size_t)L.np * sizeof *dm);
    double *sh = malloc((size_t)L.np * sizeof *sh);
    if (dm == NULL || sh == NULL) {
      free(seen);
      free(dm);
      free(sh);
      break;
    }
    int64_t nd = 0;
    for (int32_t k = 0; k < L.np; k++) {
      int32_t id = lab[k];
      if (seen[id]) continue;
      seen[id] = 1;
      cnt++;
      const hz_lodnode *n = &L.nd[id];
      dm[nd] = n->dmax;
      sh[nd] = (n->area_surf > 0.0 && n->perim > 0.0) ? n->perim / sqrt(n->area_surf) : 0.0;
      nd++;
    }
    qsort(dm, (size_t)nd, sizeof *dm, cmp_dbl);
    qsort(sh, (size_t)nd, sizeof *sh, cmp_dbl);
    printf("   уровень %d (допуск %.3f м): элементов %d; dmax p50 %.4f p90 %.4f макс %.4f; "
           "форма P/√A p50 %.2f p90 %.2f\n",
           lev, dseg * pow(2.0, lev), cnt, pct(dm, nd, 0.5), pct(dm, nd, 0.9), pct(dm, nd, 1.0),
           pct(sh, nd, 0.5), pct(sh, nd, 0.9));
    if (lev > 0 && cnt > 0)
      printf("      сокращение к предыдущему %.2f× (четвёрка — структура, не требование: §54.2)\n",
             (double)prevn / (double)cnt);
    prevn = cnt;
    fflush(stdout);
    free(seen);
    free(dm);
    free(sh);
  }

  /* --- ВЛОЖЕННОСТЬ, целочисленно (А132) --- */
  {
    int bad = 0;
    for (int lev = 1; lev < L.nlev && !bad; lev++) {
      const int32_t *p = L.lab + (size_t)(lev - 1) * (size_t)L.np;
      const int32_t *c = L.lab + (size_t)lev * (size_t)L.np;
      /* Каждый узел предыдущего уровня обязан целиком лежать в ОДНОМ узле
       * следующего: иначе разбиение не вложено. */
      int32_t *map = malloc((size_t)L.nnd * sizeof *map);
      if (map == NULL) break;
      for (int32_t i = 0; i < L.nnd; i++)
        map[i] = -1;
      for (int32_t k = 0; k < L.np; k++) {
        if (map[p[k]] < 0)
          map[p[k]] = c[k];
        else if (map[p[k]] != c[k])
          bad = 1;
      }
      free(map);
    }
    printf("   ВЛОЖЕННОСТЬ (по меткам, А132): %s\n",
           bad ? "НАРУШЕНА — разбиения уровней не вложены" : "держится на всех уровнях");
  }

  /* --- срез и однородные уровни --- */
  int32_t *cut = malloc((size_t)L.np * sizeof *cut);
  if (cut == NULL) return 1;
  const double *eye = city ? eyec : eyeh;
  int32_t ncut = hz_lod_cut(&L, eye, eps, cut);
  printf("   СРЕЗ при камере §2 (ε = %.3e рад): элементов %d\n", eps, ncut);
  {
    hz_pseglist so;
    if (hz_lod_seglist(&L, &m, &sg, cut, &so) == 0) {
      hz_polyset pc;
      if (hz_poly_build(&pc, &m, &so) == 0) {
        if (vfit) {
          hz_vfitstat vs;
          hz_poly_vfit(&pc, dseg * pow(2.0, L.nlev - 1), &vs);
          printf("      подгонка вершин: передвинуто %lld, отказов по смещению %lld\n",
                 (long long)vs.nvert_moved, (long long)vs.nvert_far);
        }
        metric(&psf, &pc, &cam, eps_px, "срез");
        hz_poly_free(&pc);
      }
      hz_seg_free(&so);
    }
  }
  /* Однородные уровни для сравнения: тот же путь, но без среза. */
  if (uniform)
    for (int lev = 1; lev < L.nlev; lev++) {
      for (int32_t k = 0; k < L.np; k++)
        cut[k] = L.lab[(size_t)lev * (size_t)L.np + (size_t)k];
      hz_pseglist so;
      if (hz_lod_seglist(&L, &m, &sg, cut, &so) != 0) continue;
      hz_polyset pc;
      if (hz_poly_build(&pc, &m, &so) == 0) {
        char tag[64];
        snprintf(tag, sizeof tag, "однородный уровень %d", lev);
        metric(&psf, &pc, &cam, eps_px, tag);
        hz_poly_free(&pc);
      }
      hz_seg_free(&so);
    }
  printf("   время: постройка %.1f с, всего %.1f с\n", t1 - t0, now_s() - t0);
  free(cut);
  hz_lod_free(&L);
  hz_poly_free(&psf);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
