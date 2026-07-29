/* pmetric — МЕТРИКА ОГРУБЛЕНИЯ, КОТОРАЯ НЕ СЛЕПА К ГЕОМЕТРИИ
 * (PLAN_ELEMENTS.md, §34, пункт О20).
 *
 * ЗАЧЕМ. Прежняя метрика (Ш7, О5) сравнивала ПОЛЕ с точным светом В ОДНОЙ И ТОЙ
 * ЖЕ точке попадания. Если геометрия уехала, уехали обе величины разом, и
 * метрика при разрушенной геометрии даже УЛУЧШАЛАСЬ: на случайном огрублении
 * (`dmax` до `101 δ`) она давала `p50 = 1.2e−7`. Такая величина не измеряет
 * огрубление, она измеряет согласованность двух расчётов между собой.
 *
 * ЧТО ЗДЕСЬ. Через ОДИН И ТОТ ЖЕ набор лучей из камеры §2 пускаются ДВЕ сцены —
 * эталонная (сегментация без огрубления) и огрублённая, — и сравнивается то, во
 * что луч попал: дальность и нормаль. Тогда уехавшая геометрия даёт расхождение
 * сразу, а не прячется.
 *
 * ТРИ ЧИСЛА О ПРОМАХАХ ПЕЧАТАЮТСЯ ОТДЕЛЬНО (А35, А107). Луч, попавший в одной
 * сцене и промахнувшийся в другой, разности не имеет; если такие молча
 * выбрасывать, метрика обнулится ровно на худших случаях — это в точности О4,
 * где число классифицированных пикселей падало `189 034 → 99 995`, а все метрики
 * становились РОВНО нулём.
 *
 * ЧЕГО ЭТА МЕТРИКА НЕ МЕРИТ, СКАЗАНО ЗАРАНЕЕ:
 *   - эталон здесь — САМ ПРИБЛИЖЕНИЕ (сегментация при `δ_seg`), поэтому меряется
 *     шаг слияния, а не полная ошибка против треугольной сетки (А105);
 *   - `|Δt|` вдоль луча — не ошибка картинки: на скользящем луче метр смещения
 *     даёт доли пикселя, на фронтальном — сотни. Поэтому печатаются ОБЕ величины,
 *     метры и пиксели (А106);
 *   - камера одна, и она смотрит туда, куда смотрит (А108);
 *   - максимум по 262 тысячам лучей неустойчив, выводы — по p50 и p90 (А109).
 */

#include "pedge.h"
#include "pmerge.h"
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

int main(int argc, char **argv) {
  int city = (argc > 1 && strcmp(argv[1], "city") == 0);
  double dseg = city ? 0.05 : 0.045;
  double dcoarse = city ? 2.0 : 0.3;
  int32_t target = city ? 20000 : 250;
  int simp = 0, random = 0;
  double conemax = 0.0, lossmax = 0.0;
  int vfit = 0;
  double vmove = 1.0; /* предел смещения вершины в долях допуска огрубления */
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "simp") == 0) simp = 1;
    if (strcmp(argv[i], "rand") == 0) random = 1;
    if (strncmp(argv[i], "d=", 2) == 0) dcoarse = strtod(argv[i] + 2, NULL);
    if (strncmp(argv[i], "cone=", 5) == 0) conemax = strtod(argv[i] + 5, NULL);
    if (strncmp(argv[i], "loss=", 5) == 0) lossmax = strtod(argv[i] + 5, NULL);
    if (strcmp(argv[i], "vfit") == 0) vfit = 1;
    /* УРОВЕНЬ ЛЕСТНИЦЫ ЗАДАЁТСЯ ДОПУСКОМ, А НЕ ЦЕЛЬЮ. С целью число полигонов на
     * грубых уровнях упирается в неё, и «сколько полигонов на уровне» меряется
     * не то: получается ответ про цель, а не про геометрию (А126). */
    if (strcmp(argv[i], "notarget") == 0) target = 0;
    if (strncmp(argv[i], "vfit=", 5) == 0) {
      vfit = 1;
      vmove = strtod(argv[i] + 5, NULL);
    }
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
  hz_mergecfg mc;
  memset(&mc, 0, sizeof mc);
  mc.delta = dcoarse;
  mc.target = target;
  mc.use_geom = 1;
  mc.use_overlap = !random;
  mc.random = random;
  mc.conemax = conemax;
  mc.lossmax = lossmax;
  hz_pseglist so;
  hz_mergestat st;
  double t0 = now_s();
  if (hz_merge(&so, &m, &sg, &psf, &mc, &st) != 0) return 1;
  hz_polyset psc;
  if (hz_poly_build(&psc, &m, &so) != 0) return 1;
  /* О17: подгонка вершин края в 3D (§40…§43). Ставится ЗДЕСЬ, после сборки
   * огрублённых полигонов: она правит только край, оставляя элемент плоским. */
  if (vfit) {
    hz_vfitstat vs;
    if (hz_poly_vfit(&psc, vmove * dcoarse, &vs) != 0) return 1;
    printf("   О17: вершин %lld, передвинуто %lld (%.2f%%), оставлено %lld, L2 вне гарантии %lld, "
           "отказов %lld; худшее t* %.4f м\n",
           (long long)vs.nvert, (long long)vs.nvert_moved,
           100.0 * (double)vs.nvert_moved / (double)(vs.nvert ? vs.nvert : 1),
           (long long)vs.nvert_fixed, (long long)vs.nvert_l2out, (long long)vs.nvert_fail, vs.tmax);
    printf("   О17: отказов по смещению %lld, худшее принятое смещение %.4f м\n",
           (long long)vs.nvert_far, vs.dmax_move);
    printf("   О17: отказов без улучшения %lld; предел смещения %.4f м\n",
           (long long)vs.nvert_noimp, vmove * dcoarse);
    printf("   О17: t* по корзинам 10^k м (от 1e-6):");
    for (int i = 0; i < 12; i++)
      if (vs.t_hist[i] > 0)
        printf(" [1e%d] %.2f%%", i - 6,
               100.0 * (double)vs.t_hist[i] / (double)(vs.nvert_moved ? vs.nvert_moved : 1));
    printf("\n");
  }
  printf("== метрика огрубления: %s, δ сегм %g, δ огр %g%s%s; участков %d -> %d, dmax_worst %.3f\n",
         city ? "ГОРОД" : "зал", dseg, dcoarse, random ? ", СЛУЧАЙНОЕ" : "",
         (conemax > 0.0) ? ", конус ограничен" : "", sg.nseg, so.nseg, st.dmax_worst);

  hz_pray gf, gc;
  if (hz_pray_build(&gf, &psf, 4.0) != 0) return 1;
  if (hz_pray_build(&gc, &psc, 4.0) != 0) return 1;

  tr3_camera cam;
  double eyeh[3] = HZ_CFG_HALL_EYE, ath[3] = HZ_CFG_HALL_AT;
  double eyec[3] = HZ_CFG_CITY_EYE, atc[3] = HZ_CFG_CITY_AT, up[3] = HZ_CFG_UP;
  const int W = 512, H = 512;
  if (tr3_camera_look(&cam, city ? eyec : eyeh, city ? atc : ath, up, HZ_CFG_FOV_DEG * M_PI / 180.0,
                      W, H) != 0)
    return 1;
  /* Угловой размер пикселя — из поля зрения и числа пикселей; он же переводит
   * смещение поверхности в пиксели (§27). */
  const double eps_px = (HZ_CFG_FOV_DEG * M_PI / 180.0) / (double)H;

  int64_t nboth = 0, nlost = 0, ngain = 0, nnone = 0;
  double *dt = malloc((size_t)W * (size_t)H * sizeof *dt);
  double *dpx = malloc((size_t)W * (size_t)H * sizeof *dpx);
  double *dang = malloc((size_t)W * (size_t)H * sizeof *dang);
  if (dt == NULL || dpx == NULL || dang == NULL) return 1;
  double t1 = now_s();
  for (int j = 0; j < H; j++)
    for (int i = 0; i < W; i++) {
      double o[3], d[3];
      tr3_camera_ray(&cam, i, j, o, d);
      double tf = 0.0, tc = 0.0;
      int32_t hf = hz_pray_hit(&gf, o, d, 1e-6, &tf);
      int32_t hc = hz_pray_hit(&gc, o, d, 1e-6, &tc);
      if (hf < 0 && hc < 0) {
        nnone++;
        continue;
      }
      if (hf >= 0 && hc < 0) {
        nlost++;
        continue;
      }
      if (hf < 0 && hc >= 0) {
        ngain++;
        continue;
      }
      double e = fabs(tf - tc);
      dt[nboth] = e;
      /* В ПИКСЕЛЯХ: смещение вдоль луча даёт видимый сдвиг только через наклон
       * поверхности; множитель — косинус между нормалью и лучом, делённый на
       * угловой след пикселя на этой дальности (§27, А106). */
      const double *nf = psf.p[hf].n;
      double cosr = fabs(nf[0] * d[0] + nf[1] * d[1] + nf[2] * d[2]);
      dpx[nboth] = e * cosr / (eps_px * (tf > 0.0 ? tf : 1e-9));
      const double *nc = psc.p[hc].n;
      double cn = nf[0] * nc[0] + nf[1] * nc[1] + nf[2] * nc[2];
      if (cn > 1.0) cn = 1.0;
      if (cn < -1.0) cn = -1.0;
      dang[nboth] = acos(fabs(cn)) * 180.0 / M_PI;
      nboth++;
    }
  double t2 = now_s();
  qsort(dt, (size_t)nboth, sizeof *dt, cmp_dbl);
  qsort(dpx, (size_t)nboth, sizeof *dpx, cmp_dbl);
  qsort(dang, (size_t)nboth, sizeof *dang, cmp_dbl);
  int64_t ntot = nboth + nlost + ngain;
  printf("   ЛУЧИ (%d×%d): попали в обеих %lld (%.2f%% от попавших хоть куда), потеряли "
         "попадание %lld (%.2f%%), приобрели %lld (%.2f%%); мимо обеих %lld\n",
         W, H, (long long)nboth, 100.0 * (double)nboth / (double)(ntot ? ntot : 1),
         (long long)nlost, 100.0 * (double)nlost / (double)(ntot ? ntot : 1), (long long)ngain,
         100.0 * (double)ngain / (double)(ntot ? ntot : 1), (long long)nnone);
  printf("   |Δt| по общим лучам, м: p50 %.4f, p90 %.4f, максимум %.4f\n", pct(dt, nboth, 0.5),
         pct(dt, nboth, 0.9), pct(dt, nboth, 1.0));
  printf("   то же В ПИКСЕЛЯХ: p50 %.2f, p90 %.2f, максимум %.2f\n", pct(dpx, nboth, 0.5),
         pct(dpx, nboth, 0.9), pct(dpx, nboth, 1.0));
  printf("   угол между нормалями, °: p50 %.2f, p90 %.2f, максимум %.2f\n", pct(dang, nboth, 0.5),
         pct(dang, nboth, 0.9), pct(dang, nboth, 1.0));
  printf("   время: слияние %.2f с, лучи %.2f с\n", t1 - t0, t2 - t1);

  free(dt);
  free(dpx);
  free(dang);
  hz_pray_free(&gf);
  hz_pray_free(&gc);
  hz_poly_free(&psc);
  hz_seg_free(&so);
  hz_poly_free(&psf);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
