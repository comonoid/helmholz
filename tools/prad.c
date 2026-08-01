/* prad — РЕНДЕРЕР БЕЗ СВИПА (PLAN_ELEMENTS.md, §84, §85).
 *
 * Вся цепочка одним проходом: иерархия полигонов с лестницей LOD → сборка
 * уравнения связями → решение итерациями → картинка.
 *
 * ЧЕМ ОН ОТЛИЧАЕТСЯ ОТ `prender`. Тот решает ординатной развёрткой (свипом) и
 * считает прямой свет отдельным механизмом. Здесь ни того, ни другого нет:
 * направление связанное (§85), поэтому оператор есть набор связей между парами
 * узлов, прямой свет — обычная связь от светящегося узла, а многократные отражения
 * — итерации того же умножения.
 *
 * ЗАЧЕМ ОБА СРАЗУ. Свип-решатель остаётся ЭТАЛОНОМ: два независимых метода на одной
 * сцене обязаны сойтись, и это сильнейшая проверка, какая тут возможна. Пока не
 * сошлись — верить нельзя ни одному.
 *
 * ИСТОЧНИК В ЭТОЙ РЕДАКЦИИ — ПОТОЛОК, и это ОСНАСТКА, а не физика. Лампы `prender`
 * суть отдельные полигоны, добавленные к набору; в лестнице их нет, потому что
 * лестница строится по геометрии сцены. Чтобы не тащить сюда всю машинерию
 * источников до того, как проверен перенос, светящимися объявляются элементы,
 * смотрящие вниз у самого верха габарита. Замена честная для проверки переноса и
 * негодная для сравнения яркостей с `prender`.
 */

#include "image.h"
#include "lodio.h"
#include "pedge.h"
#include "plink.h"
#include "plod.h"
#include "poly_seg.h"
#include "polygon.h"
#include "pray.h"
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

static int cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

int main(int argc, char **argv) {
  int city = 0, w = 512, ss = 2, nvis = 2, maxlev = 13, noself = 0, nopull = 0;
  double epsmul = 1.0, radmul = 4.0, base = 1.4142, linkmul = 1.0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "city") == 0) city = 1;
    if (strncmp(argv[i], "w=", 2) == 0) w = (int)strtol(argv[i] + 2, NULL, 10);
    if (strncmp(argv[i], "ss=", 3) == 0) ss = (int)strtol(argv[i] + 3, NULL, 10);
    if (strncmp(argv[i], "vis=", 4) == 0) nvis = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "lev=", 4) == 0) maxlev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "eps=", 4) == 0) epsmul = strtod(argv[i] + 4, NULL);
    /* Допуск СВЯЗИ в допусках среза: `1` — тот же критерий, что у камеры. */
    if (strncmp(argv[i], "link=", 5) == 0) linkmul = strtod(argv[i] + 5, NULL);
    /* ДВА ВЫКЛЮЧАТЕЛЯ — ЗАМЕРЫ, А НЕ РЕЖИМЫ (§86, А195). `noself` снимает
     * самопары и показывает, что именно они купили; `nopull` снимает подъём —
     * это негативный контроль НК1, обязанный ПРОВАЛИТЬСЯ почернением. */
    if (strcmp(argv[i], "noself") == 0) noself = 1;
    if (strcmp(argv[i], "nopull") == 0) nopull = 1;
  }
  double dseg = city ? 0.05 : 0.045;

  hz_objmesh m;
  if (hz_obj_load(&m, city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ,
                  city ? HZ_CFG_CITY_SCALE : HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  hz_pseglist sg;
  if (hz_seg_planar(&sg, &m, dseg) != 0) return 1;
  hz_polyset ps;
  if (hz_poly_build(&ps, &m, &sg) != 0) return 1;
  {
    int32_t mxl = 0;
    for (int32_t k = 0; k < ps.np; k++) {
      int32_t nb = ps.loop[ps.p[k].l0 + ps.p[k].nloop] - ps.loop[ps.p[k].l0];
      if (nb > mxl) mxl = nb;
    }
    printf("== СЦЕНА: %s, треугольников %d, участков %d; вершин края %d (на полигон %.1f, "
           "максимум %d)\n",
           city ? "ГОРОД" : "зал", m.nt, sg.nseg, ps.nbv, (double)ps.nbv / (double)(ps.np ? ps.np : 1),
           mxl);
  }

  tr3_camera cam;
  double eyeh[3] = HZ_CFG_HALL_EYE, ath[3] = HZ_CFG_HALL_AT;
  double eyec[3] = HZ_CFG_CITY_EYE, atc[3] = HZ_CFG_CITY_AT, up[3] = HZ_CFG_UP;
  double eye[3], at[3];
  for (int c = 0; c < 3; c++) {
    eye[c] = city ? eyec[c] : eyeh[c];
    at[c] = city ? atc[c] : ath[c];
  }
  for (int i = 1; i < argc; i++) {
    if (strncmp(argv[i], "eye=", 4) == 0)
      sscanf(argv[i] + 4, "%lf,%lf,%lf", &eye[0], &eye[1], &eye[2]);
    if (strncmp(argv[i], "at=", 3) == 0) sscanf(argv[i] + 3, "%lf,%lf,%lf", &at[0], &at[1], &at[2]);
  }
  if (tr3_camera_look(&cam, eye, at, up, HZ_CFG_FOV_DEG * M_PI / 180.0, w, w) != 0) return 1;
  const double eps_px = (HZ_CFG_FOV_DEG * M_PI / 180.0) / (double)w;
  const double eps = eps_px * epsmul;

  double t0 = now_s();
  hz_lod L;
  hz_lodcfg lc;
  memset(&lc, 0, sizeof lc);
  lc.delta0 = dseg;
  lc.eps = eps;
  lc.maxlev = maxlev;
  lc.radmul = radmul;
  lc.base = base;
  if (hz_lod_build_merge(&L, &m, &sg, &ps, &lc) != 0) {
    fprintf(stderr, "отказ лестницы\n");
    return 1;
  }
  printf("== ЛЕСТНИЦА: уровней %d, узлов %d, за %.1f с\n", L.nlev, L.nnd, now_s() - t0);

  /* Сетка лучей — по ИСХОДНЫМ полигонам: видимость обязана считаться по настоящей
   * геометрии, а не по огрублённой, иначе связь пройдёт сквозь стену. */
  hz_pray g;
  if (hz_pray_build(&g, &ps, 4.0) != 0) return 1;

  hz_scene sc;
  sc.L = &L;
  sc.g = &g;
  hz_linkset S;
  t0 = now_s();
  if (hz_links_build(&S, &sc, eps * linkmul, nvis, noself) != 0) {
    fprintf(stderr, "отказ сборки связей\n");
    return 1;
  }
  S.t_build = now_s() - t0;
  printf("== СВЯЗИ: %lld штук за %.1f с; пар рассмотрено %lld, дроблений %lld, отброшено %lld, "
         "проб видимости %lld, максимум связей у узла %lld\n",
         (long long)S.n, S.t_build, (long long)S.nvisit, (long long)S.nrefine, (long long)S.nzero,
         (long long)S.nray, (long long)S.nmax_node);
  printf("   из отброшенных на ГРУБОМ уровне (обрубило поддерево): %lld\n",
         (long long)S.nzero_coarse);
  printf("   на узел в среднем %.1f; связей / (n log n) = %.2f\n",
         (double)S.n / (double)(L.nnd ? L.nnd : 1),
         (double)S.n / ((double)L.nnd * log2((double)L.nnd + 2.0)));
  /* ПОДПИСЬ ЗАМКНУТОСТИ: `Σf` по листу с предками обязана быть единицей. Хвост
   * по площади, а не среднее (А194) — среднее печатается справочно. */
  printf("== Σf ПО ЛИСТУ (с предками; в замкнутой сцене = 1): p10 %.3f, p50 %.3f, p90 %.3f; "
         "среднее %.3f, min %.3f, max %.3f; площади с Σf<0.9 — %.1f %%\n",
         S.sf_p10, S.sf_p50, S.sf_p90, S.sf_mean, S.sf_min, S.sf_max, 100.0 * S.sf_lowfrac);
  fflush(stdout);

  /* ЗАМКНУТОСТЬ СЦЕНЫ — НЕЗАВИСИМАЯ ССЫЛКА ДЛЯ `Σf`. Само по себе `Σf < 1` ещё
   * не значит потери: если сцена не замкнута, это ПРАВИЛЬНЫЙ ответ. Меряется
   * лучами и потому не разделяет с `Σf` ни формулы коэффициента, ни иерархии:
   * доля косинусно-взвешенного телесного угла, упирающегося в геометрию. Для
   * замкнутой сцены она единица, и тогда `Σf` обязано быть единицей тоже.
   * Направления — решётка Хаммерсли (детерминирована, §4 запрещает случайность
   * в замерах), косинусное распределение даёт ровно вес `cosθ dω / π`. */
  {
    const int NCLO = 64; /* лучей на узел: 993·64 ≈ 64 тыс., доли секунды */
    double amiss = 0.0, atot2 = 0.0, worst = 1.0;
    for (int32_t k = 0; k < L.nnd; k++) {
      if (L.nd[k].level != 0) continue;
      const hz_lodnode *N = &L.nd[k];
      double eu[3], ev[3], t0v[3] = {0.0, 0.0, 0.0};
      int ax = 0;
      for (int c = 1; c < 3; c++)
        if (fabs(N->n[c]) < fabs(N->n[ax])) ax = c;
      t0v[ax] = 1.0;
      eu[0] = N->n[1] * t0v[2] - N->n[2] * t0v[1];
      eu[1] = N->n[2] * t0v[0] - N->n[0] * t0v[2];
      eu[2] = N->n[0] * t0v[1] - N->n[1] * t0v[0];
      double en = sqrt(eu[0] * eu[0] + eu[1] * eu[1] + eu[2] * eu[2]);
      if (!(en > 0.0)) continue;
      for (int c = 0; c < 3; c++)
        eu[c] /= en;
      ev[0] = N->n[1] * eu[2] - N->n[2] * eu[1];
      ev[1] = N->n[2] * eu[0] - N->n[0] * eu[2];
      ev[2] = N->n[0] * eu[1] - N->n[1] * eu[0];
      int hit = 0;
      for (int s = 0; s < NCLO; s++) {
        double u1 = ((double)s + 0.5) / NCLO, u2 = 0.0, f2 = 0.5;
        for (int b = s; b > 0; b >>= 1, f2 *= 0.5)
          if (b & 1) u2 += f2;
        double sn = sqrt(u1), cs = sqrt(1.0 - u1), ph = 2.0 * M_PI * u2;
        double d[3], o[3];
        for (int c = 0; c < 3; c++) {
          d[c] = sn * cos(ph) * eu[c] + sn * sin(ph) * ev[c] + cs * N->n[c];
          o[c] = ((c == 0) ? N->cx : (c == 1) ? N->cy : N->cz) + 1e-5 * N->n[c];
        }
        if (hz_pray_occluded(&g, o, d, 0.0, 1e6)) hit++;
      }
      double frac = (double)hit / NCLO;
      atot2 += N->area_surf;
      amiss += N->area_surf * (1.0 - frac);
      if (frac < worst) worst = frac;
    }
    printf("== ЗАМКНУТОСТЬ (лучами, независимо от Σf): в геометрию упирается %.3f "
           "косинусного телесного угла по площади; худший узел %.3f\n",
           (atot2 > 0.0) ? 1.0 - amiss / atot2 : 0.0, worst);
    fflush(stdout);
  }

  /* --- ИСТОЧНИК: ПОТОЛОК (оснастка, см. заголовок) --- */
  double *Le = calloc((size_t)L.nnd, sizeof *Le);
  double *rho = malloc((size_t)L.nnd * sizeof *rho);
  double *E = malloc((size_t)L.nnd * sizeof *E);
  if (Le == NULL || rho == NULL || E == NULL) return 1;
  int64_t nsrc = 0;
  double asrc = 0.0;
  double ytop = m.hi[1] - 0.05 * (m.hi[1] - m.lo[1]);
  for (int32_t k = 0; k < L.nnd; k++) {
    rho[k] = 0.5;
    /* ИСТОЧНИК ТОЛЬКО НА НУЛЕВОМ УРОВНЕ. Прежде критерий применялся к каждому
     * уровню независимо (2 узла на нулевом против 11 по лестнице), а подъём
     * затем перезаписывал `B` у всех внутренних узлов — излучение грубых просто
     * выбрасывалось. На грубые уровни оно теперь приходит подъёмом, который и
     * так есть. */
    if (L.nd[k].level == 0 && L.nd[k].n[1] < -0.9 && L.nd[k].cy > ytop) {
      Le[k] = HZ_CFG_LAMP_LE;
      rho[k] = 0.0;
      nsrc++;
      asrc += L.nd[k].area_surf;
    }
  }
  printf("== ИСТОЧНИК (оснастка): светящихся узлов уровня 0 — %lld выше y = %.2f, "
         "площадь %.3f м²\n",
         (long long)nsrc, ytop, asrc);

  const int maxit = 200;
  double *rh = calloc((size_t)maxit, sizeof *rh);
  if (rh == NULL) return 1;
  t0 = now_s();
  int it = hz_links_solve(&S, &L, Le, rho, E, 1e-4, maxit, rh, nopull);
  double tsolve = now_s() - t0;
  printf("== РЕШЕНИЕ: %d итераций за %.3f с (%.1f мкс на итерацию)%s\n", it, tsolve,
         1e6 * tsolve / (double)(it > 0 ? it : 1), nopull ? " [НК1: подъём выключен]" : "");
  printf("   невязки:");
  for (int i = 0; i < it && i < 12; i++)
    printf(" %.2e", rh[i]);
  if (it > 12) printf(" ... %.2e", rh[it - 1]);
  printf("\n   отношение соседних (измеренное сжатие):");
  for (int i = 1; i < it && i < 12; i++)
    printf(" %.3f", (rh[i - 1] > 0.0) ? rh[i] / rh[i - 1] : 0.0);
  printf("\n");
  free(rh);
  /* ПОТОК ПО НУЛЕВОМУ УРОВНЮ — величина, по которой сравниваются выключатели. */
  double flux0 = 0.0;
  for (int32_t k = 0; k < L.nnd; k++)
    if (L.nd[k].level == 0) flux0 += E[k] * L.nd[k].area_surf;
  printf("== ПОТОК: Σ B·A по уровню 0 = %.6e Вт/ср%s%s\n", flux0, noself ? " [noself]" : "",
         nopull ? " [nopull]" : "");

  /* --- СРЕЗ И КАРТИНКА --- */
  int32_t *cut = malloc((size_t)L.np * sizeof *cut);
  if (cut == NULL) return 1;
  int32_t ncut = hz_lod_cut(&L, eye, eps, 1, cut);
  printf("== СРЕЗ: %d элементов\n", ncut);

  int n = w * w;
  double *L3 = calloc((size_t)n * 3, sizeof *L3);
  unsigned char *rgb = malloc((size_t)n * 3);
  if (L3 == NULL || rgb == NULL) return 1;
  t0 = now_s();
#pragma omp parallel for schedule(dynamic, 16)
  for (int i = 0; i < n; i++) {
    double acc[3] = {0.0, 0.0, 0.0};
    for (int sy = 0; sy < ss; sy++)
      for (int sx = 0; sx < ss; sx++) {
        double o[3], d[3], t;
        tr3_camera_ray(&cam, i % w, i / w, o, d);
        int32_t k = hz_pray_hit(&g, o, d, 0.0, &t);
        if (k < 0) continue;
        /* Поле читается у УЗЛА СРЕЗА, а альбедо и цвет — у МЕЛКОГО полигона:
         * правило Т1, деталь достаётся даром и от огрубления не зависит. */
        int32_t nd = cut[k];
        double Ek = E[nd] * M_PI; /* E хранит РАДИАНС узла */
        int32_t tr = (ps.p[k].ntri > 0) ? ps.tri[ps.p[k].t0] : -1;
        int32_t mi = (tr >= 0 && m.fm != NULL) ? m.fm[tr] : 0;
        for (int c = 0; c < 3; c++) {
          double a = (mi >= 0 && mi < m.nmtl) ? m.mtl[mi].kd3[c] : 0.5;
          /* Узел несёт РАДИАНС (уже с альбедо своего уровня); цвет берётся у
           * мелкого полигона делением на среднее альбедо узла — правило Т1. */
          acc[c] += Ek * a / (0.5 * M_PI);
        }
      }
    double wgt = 1.0 / (double)(ss * ss);
    for (int c = 0; c < 3; c++)
      L3[(size_t)i * 3 + (size_t)c] = acc[c] * wgt;
  }
  double tframe = now_s() - t0;

  double *lum = malloc((size_t)n * sizeof *lum);
  if (lum == NULL) return 1;
  for (int i = 0; i < n; i++)
    lum[i] = 0.2126 * L3[3 * i] + 0.7152 * L3[3 * i + 1] + 0.0722 * L3[3 * i + 2];
  qsort(lum, (size_t)n, sizeof *lum, cmp_d);
  double mx = lum[(size_t)((double)n * 0.995)];
  free(lum);
  if (!(mx > 0.0)) mx = 1.0;
  for (int i = 0; i < 3 * n; i++) {
    double x = L3[i] / mx;
    if (x < 0.0) x = 0.0;
    if (x > 1.0) x = 1.0;
    rgb[i] = (unsigned char)(255.0 * pow(x, 1.0 / 2.2) + 0.5);
  }
  char path[96];
  snprintf(path, sizeof path, "img/prad_%s.ppm", city ? "city" : "hall");
  if (hz_ppm_write_rgb(path, rgb, w, w) != 0) return 1;
  printf("== КАДР: %s, %d×%d, суперсэмплинг %d×%d, за %.2f с\n", path, w, w, ss, ss, tframe);

  free(L3);
  free(rgb);
  free(cut);
  free(Le);
  free(rho);
  free(E);
  hz_links_free(&S);
  hz_pray_free(&g);
  hz_lod_free(&L);
  hz_poly_free(&ps);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
