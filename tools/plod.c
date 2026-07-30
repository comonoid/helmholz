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
#include "lodio.h"
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
                   const char *tag, int32_t nelem) {
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
  /* СИЛУЭТНЫЙ КАНАЛ (§69). Пиксельная метрика мерит СМЕЩЕНИЕ ПЛОСКОСТИ вдоль
   * луча, а ломается картинка на СИЛУЭТЕ — там, где огрублённый край врезается в
   * очертание или отступает от него. До сих пор это схлопывалось в проценты
   * «потеряли/приобрели», то есть в площадь, а видна глазу ШИРИНА полосы.
   * Ширина берётся без единого порога: `(потеряли + приобрели) / длина силуэта`,
   * где длина силуэта — число пикселей эталона, у которых сосед по четырём
   * направлениям УЛЕТЕЛ В НЕБО. Это внешний силуэт, определённый точно: небо есть
   * промах луча, а не значение глубины, и порога тут нет по построению. */
  uint8_t *hit = calloc((size_t)W * (size_t)H, 1);
  if (dt == NULL || dp == NULL || hit == NULL) {
    free(dt);
    free(dp);
    free(hit);
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
      if (hf >= 0) hit[(size_t)j * (size_t)W + (size_t)i] = 1;
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
  /* ХВОСТ ПЕЧАТАЕТСЯ ОБЯЗАТЕЛЬНО, И ЭТО НЕ УКРАШЕНИЕ (А158). Первый городской
   * прогон дал у среза и у ПЕРЕМЕШАННОГО контроля одинаковые `p50` и `p90` — оба
   * ноль, — из чего следовало бы, что критерий среза пуст. Но порча от
   * перестановки затрагивает малую долю пикселей, и `p90` её не видит по
   * построению: при 99.86 % целых пикселей девяностый процентиль стоит глубоко
   * внутри целой части. Различать обязаны `p99`, `p99.9`, максимум и ДОЛЯ
   * пикселей грубее одного — они и печатаются. */
  int64_t nbad = 0;
  for (int64_t i = 0; i < nb; i++)
    if (dp[i] > 1.0) nbad++;
  /* Длина внешнего силуэта эталона: попавший пиксель, у которого сосед — небо. */
  int64_t nsil = 0;
  for (int j = 0; j < H; j++)
    for (int i = 0; i < W; i++) {
      size_t q = (size_t)j * (size_t)W + (size_t)i;
      if (!hit[q]) continue;
      int edge = (i == 0 || i == W - 1 || j == 0 || j == H - 1);
      if (!edge) {
        if (!hit[q - 1] || !hit[q + 1] || !hit[q - (size_t)W] || !hit[q + (size_t)W]) edge = 1;
      }
      if (edge) nsil++;
    }
  printf("      %-28s элем %7d | попали %6.2f%% потеряли %5.2f%% приобрели %5.2f%% | "
         "ПИКСЕЛИ p50 %6.2f p90 %6.2f p99 %7.2f p99.9 %8.2f макс %8.2f | >1пкс %6.3f%% | "
         "СИЛУЭТ: длина %6lld пкс, сдвиг %6.2f пкс\n",
         tag, nelem, 100.0 * (double)nb / (double)(tot ? tot : 1),
         100.0 * (double)nl / (double)(tot ? tot : 1), 100.0 * (double)ng / (double)(tot ? tot : 1),
         pct(dp, nb, 0.5), pct(dp, nb, 0.9), pct(dp, nb, 0.99), pct(dp, nb, 0.999),
         pct(dp, nb, 1.0), 100.0 * (double)nbad / (double)(nb ? nb : 1), (long long)nsil,
         (nsil > 0) ? (double)(nl + ng) / (double)nsil : 0.0);
  free(dt);
  free(dp);
  free(hit);
  hz_pray_free(&gf);
  hz_pray_free(&gc);
}

int main(int argc, char **argv) {
  int city = (argc > 1 && strcmp(argv[1], "city") == 0);
  double dseg = city ? 0.05 : 0.045;
  int simp = 0, vfit = 0, maxlev = 9, uniform = 0, viamerge = 0, bands = 0, bycount = 0,
      byangle = 0, curve = 0;
  double epsmul = 1.0, radmul = 1.0;
  const char *save = NULL, *load = NULL;
  /* НЕИЗВЕСТНЫЙ АРГУМЕНТ — ОШИБКА, А НЕ ПРОПУСК (§60, дефект оснастки). Флаг,
   * который не совпал, молчал, и конфигурация вышла тождественной другой; поймать
   * это удалось лишь по совпадению всех семи чисел. Обрыв дешевле. */
  for (int i = 1; i < argc; i++) {
    int ok = 0;
    if (strcmp(argv[i], "city") == 0 || strcmp(argv[i], "hall") == 0) ok = 1;
    if (strcmp(argv[i], "simp") == 0) ok = simp = 1;
    if (strcmp(argv[i], "vfit") == 0) ok = vfit = 1;
    if (strcmp(argv[i], "uniform") == 0) ok = uniform = 1;
    /* Второй построитель: уровень строится `hz_merge` над предыдущим (§58). */
    if (strcmp(argv[i], "viamerge") == 0) ok = viamerge = 1;
    if (strcmp(argv[i], "bands") == 0) ok = bands = 1;
    /* ВТОРАЯ ТАКТИКА (§60): уровень задан ЧИСЛОМ (четверть), допуск не
     * ограничивает, достигнутое отклонение печатается замером. */
    if (strcmp(argv[i], "bycount") == 0) ok = bycount = 1;
    /* ТРЕТЬЯ ТАКТИКА (§67): критерий уровня — УГОЛ (полураствор конуса нормалей). */
    if (strcmp(argv[i], "byangle") == 0) ok = byangle = 1;
    /* §69: КРИВАЯ «ЭЛЕМЕНТЫ ПРОТИВ ОШИБКИ» ПО СВИПУ ε. Одна точка тактики не
     * решает ничего: на зале срез оставляет 865 полигонов из 993 на нулевом
     * уровне, то есть сравниваются сцены, совпадающие на 87 %. Сравнивать надо
     * при РАВНОЙ ЦЕНЕ, а цена задаётся ε. */
    if (strcmp(argv[i], "curve") == 0) ok = curve = 1;
    if (strncmp(argv[i], "lev=", 4) == 0) {
      maxlev = (int)strtol(argv[i] + 4, NULL, 10);
      ok = 1;
    }
    /* §68: РАДИУС ПОИСКА КАНДИДАТОВ В ДОПУСКАХ УРОВНЯ. */
    if (strncmp(argv[i], "rad=", 4) == 0) {
      radmul = strtod(argv[i] + 4, NULL);
      ok = 1;
    }
    if (strncmp(argv[i], "eps=", 4) == 0) {
      epsmul = strtod(argv[i] + 4, NULL);
      ok = 1;
    }
    /* ЛЕСТНИЦА КАК АРТЕФАКТ СЦЕНЫ (§62, пункт 2): строится один раз, дальше
     * читается. `load` заменяет постройку целиком. */
    if (strncmp(argv[i], "save=", 5) == 0) {
      save = argv[i] + 5;
      ok = 1;
    }
    if (strncmp(argv[i], "load=", 5) == 0) {
      load = argv[i] + 5;
      ok = 1;
    }
    if (!ok) {
      fprintf(
          stderr,
          "plod: неизвестный аргумент «%s»\n"
          "  ожидается: [city|hall] [simp] [vfit] [uniform] [viamerge] [bands]\n"
          "             [bycount] [byangle] [curve] [rad=X] [lev=N] [eps=X] [save=Ф] [load=Ф]\n",
          argv[i]);
      return 2;
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
  int rc;
  if (load != NULL) {
    /* ЧТЕНИЕ ВМЕСТО ПОСТРОЙКИ. Проверки — внутри `hz_lod_read`, включая полную
     * целочисленную вложенность; отказ печатается СЛОВАМИ, потому что «не ноль»
     * не различает чужой файл, порчу и несогласованность. */
    int fsimp = 0;
    rc = hz_lod_read(&L, load, &sg, m.nt, &fsimp);
    if (rc != HZ_LODIO_OK) {
      fprintf(stderr, "лестница не прочитана: %s (код %d)\n", hz_lodio_str(rc), rc);
      return 1;
    }
    /* УПРОЩЕНИЕ КРАЯ СВЕРЯЕТСЯ ОТДЕЛЬНО: сумма разметки его не видит (см. оговорку
     * в `lodio.h`), а границы полигонов от него зависят. */
    if (fsimp != (simp ? 1 : 0)) {
      fprintf(stderr, "лестница построена %s края, а запрошено %s — отказ\n",
              fsimp ? "С УПРОЩЕНИЕМ" : "БЕЗ УПРОЩЕНИЯ", simp ? "С УПРОЩЕНИЕМ" : "БЕЗ УПРОЩЕНИЯ");
      hz_lod_free(&L);
      return 1;
    }
    /* КАМЕРА СВЕРЯЕТСЯ С ФАЙЛОМ (А131). Лестница строилась под своё `ε`, и срез
     * под другим `ε` дал бы ступени, не совпадающие с построенными. Отказ, а не
     * пересчёт. */
    if (!(fabs(L.eps - eps) <= 1e-12 * (fabs(eps) + 1.0))) {
      fprintf(stderr, "лестница построена под ε = %.9g, а запрошено %.9g — отказ (А131)\n", L.eps,
              eps);
      hz_lod_free(&L);
      return 1;
    }
    dseg = L.delta0;
  } else {
    rc = viamerge ? hz_lod_build_merge(&L, &m, &sg, &psf, dseg, maxlev, eps, bands, bycount,
                                       byangle, radmul)
                  : hz_lod_build(&L, &m, &sg, &psf, dseg, maxlev, eps);
    if (rc != 0) {
      fprintf(stderr, "отказ построения лестницы\n");
      return 1;
    }
  }
  double t1 = now_s();
  printf("== ЛЕСТНИЦА: %s, δ0 %g м, уровней %d, узлов %d, %s за %.1f с\n", city ? "ГОРОД" : "зал",
         dseg, L.nlev, L.nnd, load ? "ПРОЧИТАНА" : "построена", t1 - t0);
  if (save != NULL) {
    /* ЗАПИСЬ СРАЗУ ПОСЛЕ ПОСТРОЙКИ, ДО всякого замера: если дальше что-то упадёт,
     * двести секунд постройки уже не потеряны. */
    double ts = now_s();
    int wr = hz_lod_write(&L, save, &sg, m.nt, simp);
    if (wr != HZ_LODIO_OK) {
      fprintf(stderr, "лестница не записана: %s (код %d)\n", hz_lodio_str(wr), wr);
      hz_lod_free(&L);
      return 1;
    }
    printf("== ЗАПИСАНА: %s за %.2f с\n", save, now_s() - ts);
    fflush(stdout);
  }
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
    if (L.bycount)
      printf("   уровень %d (ЦЕЛЬ по числу): элементов %d; dmax ЗАМЕР p50 %.4f p90 %.4f макс %.4f; "
             "форма P/√A p50 %.2f p90 %.2f\n",
             lev, cnt, pct(dm, nd, 0.5), pct(dm, nd, 0.9), pct(dm, nd, 1.0), pct(sh, nd, 0.5),
             pct(sh, nd, 0.9));
    else
      printf("   уровень %d (допуск %.3f м): элементов %d; dmax p50 %.4f p90 %.4f макс %.4f; "
             "форма P/√A p50 %.2f p90 %.2f\n",
             lev, dseg * pow(2.0, lev), cnt, pct(dm, nd, 0.5), pct(dm, nd, 0.9), pct(dm, nd, 1.0),
             pct(sh, nd, 0.5), pct(sh, nd, 0.9));
    if (lev > 0 && cnt > 0)
      printf("      сокращение к предыдущему %.2f× (%s)\n", (double)prevn / (double)cnt,
             L.bycount ? "четвёрка — ТРЕБОВАНИЕ: §60"
                       : "четвёрка — структура, не требование: §54.2");
    prevn = cnt;
    fflush(stdout);
    free(seen);
    free(dm);
    free(sh);
  }

  /* --- ДВУГРАННЫЙ УГОЛ ПАР-КАНДИДАТОВ (§59) ---
   * Печатается характеристика СЦЕНЫ, а не алгоритма: где стоят изломы и с какой
   * стороны. `180°` — плоское; левее выпуклое (внешний угол), правее вогнутое. */
  if (L.ndih_conv + L.ndih_conc + L.ndih_flat > 0) {
    printf("   двугранный угол кандидатов: плоских (180°) %lld, выпуклых (<180°) %lld, "
           "вогнутых (>180°) %lld\n",
           (long long)L.ndih_flat, (long long)L.ndih_conv, (long long)L.ndih_conc);
    printf("      корзины по 30°:");
    for (int q = 0; q < 12; q++)
      printf(" %d-%d:%lld", q * 30, q * 30 + 30, (long long)L.dih_hist[q]);
    printf("\n");
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
        metric(&psf, &pc, &cam, eps_px, "СРЕЗ по камере", pc.np);
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
        metric(&psf, &pc, &cam, eps_px, tag, pc.np);
        hz_poly_free(&pc);
      }
      hz_seg_free(&so);
    }
  /* --- КРИВАЯ ПО ε: ЦЕНА ПРОТИВ ОШИБКИ ДЛЯ ЭТОЙ ТАКТИКИ (§69) --- */
  if (curve) {
    const double mul[6] = {0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
    for (int i = 0; i < 6; i++) {
      hz_lod_cut(&L, eye, eps * mul[i], cut);
      hz_pseglist so;
      if (hz_lod_seglist(&L, &m, &sg, cut, &so) != 0) continue;
      hz_polyset pc;
      if (hz_poly_build(&pc, &m, &so) == 0) {
        char tag[64];
        snprintf(tag, sizeof tag, "КРИВАЯ: срез ε × %.2f", mul[i]);
        metric(&psf, &pc, &cam, eps_px, tag, pc.np);
        hz_poly_free(&pc);
      }
      hz_seg_free(&so);
    }
  }

  /* --- НЕГАТИВНЫЙ КОНТРОЛЬ: ТЕ ЖЕ УРОВНИ, РОЗДАННЫЕ ДРУГИМ ПОЛИГОНАМ (А146) ---
   *
   * Случайный ПРЕДОК не годится: он меняет цену, и провал вышел бы по причине, к
   * критерию среза отношения не имеющей. Перестановка сохраняет мультимножество
   * выбранных уровней ТОЧНО и рвёт только связь «грубее там, где дальше».
   * Перестановка ДЕТЕРМИНИРОВАННАЯ — хешем от номера: замер обязан повторяться.
   * Печатается и гистограмма уровней (обязана совпасть с гистограммой среза), и
   * итоговое число элементов — оно разойтись МОЖЕТ, и это не порок контроля, а
   * его цена: два соседних полигона на разных уровнях дают два узла там, где срез
   * давал один (А155). */
  if (uniform) {
    int32_t *lv = malloc((size_t)L.np * sizeof *lv);
    int64_t hcut[64], hperm[64];
    for (int i = 0; i < 64; i++)
      hcut[i] = hperm[i] = 0;
    if (lv != NULL) {
      /* СРЕЗ ПЕРЕСЧИТЫВАЕТСЯ ЗАНОВО: цикл по однородным уровням выше пишет в тот
       * же `cut`, и без этого контроль брал бы уровни ПОСЛЕДНЕГО однородного
       * уровня вместо среза. Поймано печатью гистограммы (А155) в первом же
       * прогоне: она показала «всё на уровне 6» там, где срез сидит на нулевом. */
      hz_lod_cut(&L, eye, eps, cut);
      for (int32_t k = 0; k < L.np; k++) {
        lv[k] = L.nd[cut[k]].level;
        if (lv[k] >= 0 && lv[k] < 64) hcut[lv[k]]++;
      }
      for (int32_t k = 0; k < L.np; k++) {
        uint64_t x = (uint64_t)(uint32_t)k * 0x9e3779b97f4a7c15ull;
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ull;
        x ^= x >> 27;
        int32_t j = (int32_t)(x % (uint64_t)(uint32_t)L.np);
        int32_t t = lv[k];
        lv[k] = lv[j];
        lv[j] = t;
      }
      for (int32_t k = 0; k < L.np; k++) {
        int32_t l = lv[k];
        if (l < 0) l = 0;
        if (l >= L.nlev) l = L.nlev - 1;
        if (l < 64) hperm[l]++;
        cut[k] = L.lab[(size_t)l * (size_t)L.np + (size_t)k];
      }
      int same = 1;
      for (int i = 0; i < 64; i++)
        if (hcut[i] != hperm[i]) same = 0;
      printf("   КОНТРОЛЬ (перестановка уровней): гистограмма уровней %s;",
             same ? "СОВПАЛА ТОЧНО — цена сохранена" : "РАЗОШЛАСЬ — контроль недоказателен");
      for (int i = 0; i < L.nlev && i < 64; i++)
        printf(" %d:%lld", i, (long long)hcut[i]);
      printf("\n");
      hz_pseglist so;
      if (hz_lod_seglist(&L, &m, &sg, cut, &so) == 0) {
        hz_polyset pc;
        if (hz_poly_build(&pc, &m, &so) == 0) {
          metric(&psf, &pc, &cam, eps_px, "КОНТРОЛЬ: уровни перемешаны", pc.np);
          hz_poly_free(&pc);
        }
        hz_seg_free(&so);
      }
      free(lv);
    }
  }

  /* --- СВИП ПО `ε`: ПОДПИСЬ АРТЕФАКТА (а) ---
   * Число элементов среза ОБЯЗАНО меняться от `ε` — от той самой величины, которой
   * его меняют. Метрика здесь НЕ считается (А157): вопрос свипа не про качество, а
   * про то, что критерий вообще работает, и ответ на него стоит секунды.
   * Рядом печатается доля полигонов, СЕВШИХ НА ПОСЛЕДНИЙ УРОВЕНЬ (А156): пока она
   * мала, свип мерит КРИТЕРИЙ; когда велика — ПОТОЛОК лестницы, и это разные вещи.
   * Оговорка: лестница откалибрована под своё `ε`, здесь меняется только срез. */
  {
    printf("   СВИП ПО ε (подпись артефакта; лестница откалибрована под ε = %.3e):\n", L.eps);
    const double mul[4] = {0.5, 1.0, 2.0, 4.0};
    for (int i = 0; i < 4; i++) {
      int32_t nc = hz_lod_cut(&L, eye, eps * mul[i], cut);
      int64_t top = 0;
      for (int32_t k = 0; k < L.np; k++)
        if (L.nd[cut[k]].level == L.nlev - 1) top++;
      printf("      ε × %.1f = %.3e: узлов среза %7d, на последнем уровне %.2f%% полигонов\n",
             mul[i], eps * mul[i], nc, 100.0 * (double)top / (double)(L.np ? L.np : 1));
    }
  }
  printf("   время: постройка %.1f с, всего %.1f с\n", t1 - t0, now_s() - t0);
  free(cut);
  hz_lod_free(&L);
  hz_poly_free(&psf);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
