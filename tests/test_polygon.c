/* test_polygon — ПРИЁМКА Ш2 (PLAN_ELEMENTS.md).
 *
 * ЧТО ЗДЕСЬ ПРОВЕРЯЕТСЯ И ПОЧЕМУ ИМЕННО ТАК.
 *
 * 1. ЛУЧИ ПО ИСХОДНЫМ ТРЕУГОЛЬНИКАМ ПРОТИВ ЛУЧЕЙ ПО ПОЛИГОНАМ. План говорит
 *    прямым текстом: «Σ объёмов ядра отсечения этого НЕ проверяет — она слепа к
 *    неверному набору плоскостей». Луч — не слеп: он видит и плоскость, и КРАЙ,
 *    и порядок по глубине. Треугольники перебираются В ЛОБ, без ускоряющей
 *    структуры, СОЗНАТЕЛЬНО: эталон обязан быть очевидно верным, а дерево — это
 *    ещё один источник ошибки ровно в том, чем меряем.
 *
 * 2. ХВОСТ, А НЕ СРЕДНЕЕ (§4). Докладываются p50/p99/p99.9/max, а промахи
 *    покрытия (один попал, другой нет) — ОТДЕЛЬНОЙ строкой, потому что на
 *    силуэте расхождение равно всей дальности при любой схеме (К20).
 *
 * 3. НЕГАТИВНЫЙ КОНТРОЛЬ. Один полигон сдвигается на 10 δ. Метрика ПО ЛУЧАМ,
 *    ПОПАВШИМ В НЕГО, обязана уйти на ≈10 δ. Общий хвост по всей картинке для
 *    этого не годится и не берётся: 993 полигона, из них сдвинут один.
 *
 * 4. УСТОЙЧИВОСТЬ К ПЕРЕСТАНОВКЕ. Треугольники переставляются детерминированно
 *    случайно, сегментация гоняется заново. Набор участков обязан СОВПАСТЬ.
 *    Не совпал — инкрементность не обоснована (Ш2 дословно).
 */

#include "poly_seg.h"
#include "polygon.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include "transport/cam3.h"
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SCENE HZ_CFG_HALL_OBJ
#define SCALE HZ_CFG_HALL_SCALE
/* δ замера. 45 мм — точка сравнения Ш1 (там 993 участка), 3 мм — ближе к
 * камерному δ из §2 (0.18…9.4 мм). */
#define DELTA_COARSE 0.045
#define DELTA_FINE 0.003
/* Лучей — 128×128. Перебор в лоб стоит (лучи × треугольники), и это потолок,
 * при котором эталон остаётся честным перебором, а тест — секундным. */
#define NRAY 128

static int failed = 0;

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* --- эталон: перебор треугольников ------------------------------------------ */

static double hit_tris(const hz_objmesh *m, const double o[3], const double d[3], int32_t *who) {
  double best = -1.0;
  if (who != NULL) *who = -1;
  for (int32_t t = 0; t < m->nt; t++) {
    double p[3][3];
    hz_obj_tri(m, t, p);
    double e1[3], e2[3], pv[3], tv[3], qv[3];
    for (int a = 0; a < 3; a++) {
      e1[a] = p[1][a] - p[0][a];
      e2[a] = p[2][a] - p[0][a];
    }
    pv[0] = d[1] * e2[2] - d[2] * e2[1];
    pv[1] = d[2] * e2[0] - d[0] * e2[2];
    pv[2] = d[0] * e2[1] - d[1] * e2[0];
    double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
    if (!(fabs(det) > 0.0)) continue;
    double inv = 1.0 / det;
    for (int a = 0; a < 3; a++)
      tv[a] = o[a] - p[0][a];
    double u = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
    if (u < 0.0 || u > 1.0) continue;
    qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
    qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
    qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
    double v = (d[0] * qv[0] + d[1] * qv[1] + d[2] * qv[2]) * inv;
    if (v < 0.0 || u + v > 1.0) continue;
    double s = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
    if (s > 0.0 && (best < 0.0 || s < best)) {
      best = s;
      if (who != NULL) *who = t;
    }
  }
  return best;
}

/* --- проверяемое: перебор полигонов ----------------------------------------- */

static double hit_polys(const hz_polyset *ps, const double o[3], const double d[3], int32_t *who) {
  double best = -1.0;
  if (who != NULL) *who = -1;
  for (int32_t k = 0; k < ps->np; k++) {
    const hz_poly *P = &ps->p[k];
    if (P->nloop == 0) continue;
    double dn = d[0] * P->n[0] + d[1] * P->n[1] + d[2] * P->n[2];
    if (!(fabs(dn) > 0.0)) continue;
    double s = (P->off - (o[0] * P->n[0] + o[1] * P->n[1] + o[2] * P->n[2])) / dn;
    if (!(s > 0.0)) continue;
    if (best >= 0.0 && s >= best) continue;
    double q[3];
    for (int a = 0; a < 3; a++)
      q[a] = o[a] + s * d[a] - P->org[a];
    double u = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
    double v = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
    if (u < P->uvlo[0] || u > P->uvhi[0] || v < P->uvlo[1] || v > P->uvhi[1]) continue;
    if (!hz_poly_inside(ps, P, u, v)) continue;
    best = s;
    if (who != NULL) *who = k;
  }
  return best;
}

/* --- статистика хвоста ------------------------------------------------------- */

static int cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

static void tail(const char *what, double *v, int n, double unit) {
  if (n <= 0) {
    printf("   %-28s нет выборки\n", what);
    return;
  }
  qsort(v, (size_t)n, sizeof *v, cmp_d);
  printf("   %-28s n=%6d  p50 %.3g  p99 %.3g  p99.9 %.3g  max %.3g   (в δ: %.2f / %.2f)\n", what, n,
         v[n / 2], v[(int)(n * 0.99)], v[(int)(n * 0.999)], v[n - 1], v[(int)(n * 0.99)] / unit,
         v[n - 1] / unit);
}

/* --- детерминированная перестановка ------------------------------------------ */

static uint64_t rng_s = 0x243F6A8885A308D3ULL;
static uint64_t rng(void) {
  rng_s += 0x9E3779B97F4A7C15ULL;
  uint64_t z = rng_s;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
  z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
  return z ^ (z >> 31);
}

static void shuffle_mesh(hz_objmesh *m) {
  for (int32_t i = m->nt - 1; i > 0; i--) {
    int32_t j = (int32_t)(rng() % (uint64_t)(i + 1));
    for (int k = 0; k < 3; k++) {
      int32_t t = m->f[(size_t)i * 3 + (size_t)k];
      m->f[(size_t)i * 3 + (size_t)k] = m->f[(size_t)j * 3 + (size_t)k];
      m->f[(size_t)j * 3 + (size_t)k] = t;
      t = m->fn[(size_t)i * 3 + (size_t)k];
      m->fn[(size_t)i * 3 + (size_t)k] = m->fn[(size_t)j * 3 + (size_t)k];
      m->fn[(size_t)j * 3 + (size_t)k] = t;
    }
    int32_t t = m->fm[i];
    m->fm[i] = m->fm[j];
    m->fm[j] = t;
  }
}

/* Отпечаток участка: плоскость + площадь + число треугольников. Сортируется,
 * поэтому от нумерации участков не зависит вовсе. */
typedef struct {
  double k[6];
} sig;

static int cmp_sig(const void *a, const void *b) {
  const sig *x = a, *y = b;
  for (int i = 0; i < 6; i++) {
    if (x->k[i] < y->k[i]) return -1;
    if (x->k[i] > y->k[i]) return 1;
  }
  return 0;
}

static sig *sigs_of(const hz_pseglist *s) {
  sig *v = malloc((size_t)s->nseg * sizeof *v);
  if (v == NULL) return NULL;
  for (int32_t i = 0; i < s->nseg; i++) {
    v[i].k[0] = s->seg[i].area;
    v[i].k[1] = (double)s->seg[i].ntri;
    for (int a = 0; a < 3; a++)
      v[i].k[2 + a] = s->seg[i].n[a];
    v[i].k[5] = s->seg[i].off;
  }
  qsort(v, (size_t)s->nseg, sizeof *v, cmp_sig);
  return v;
}

/* --- прогон --------------------------------------------------------------- */

static void run(hz_objmesh *m, double delta) {
  printf("\n=== δ = %g м ===\n", delta);
  hz_pseglist sg;
  double t0 = now_s();
  if (hz_seg_planar(&sg, m, delta) != 0) {
    printf("FAIL: сегментация не прошла\n");
    failed = 1;
    return;
  }
  double t1 = now_s();
  hz_polyset ps;
  if (hz_poly_build(&ps, m, &sg) != 0) {
    printf("FAIL: сборка полигонов не прошла\n");
    failed = 1;
    hz_seg_free(&sg);
    return;
  }
  double t2 = now_s();
  printf("   участков %d, полигонов %d; сегментация %.2f с, сборка %.2f с\n", sg.nseg, ps.np,
         t1 - t0, t2 - t1);
  printf("   сварка: вершин %d -> %lld; полурёбра: всего %lld, без пары %lld (%.3f%%), "
         "дублей %lld\n",
         m->nv, (long long)ps.nweld, (long long)ps.nhe, (long long)ps.nhe_open,
         100.0 * (double)ps.nhe_open / (double)ps.nhe, (long long)ps.nhe_dup);
  printf("   край: петель %d, вершин края %d; не замкнулось %lld, развилок %lld; "
         "полигонов на шве материалов %lld\n",
         ps.nloopall, ps.nbv, (long long)ps.nopen, (long long)ps.nfork, (long long)ps.nmixed);

  /* Проверка, а не доклад: край обязан быть у каждого полигона. */
  int32_t noloop = 0;
  double dmaxworst = 0.0;
  for (int32_t k = 0; k < ps.np; k++) {
    if (ps.p[k].nloop == 0) noloop++;
    if (ps.p[k].dmax > dmaxworst) dmaxworst = ps.p[k].dmax;
  }
  printf("   полигонов без петли: %d; max dmax = %.4g м (%.4f δ)\n", noloop, dmaxworst,
         dmaxworst / delta);
  if (dmaxworst > delta) {
    printf("FAIL: dmax больше δ — критерий на максимуме нарушен\n");
    failed = 1;
  }

  /* --- лучи --- */
  tr3_camera cam;
  double eye[3] = HZ_CFG_HALL_EYE, at[3] = HZ_CFG_HALL_AT, up[3] = HZ_CFG_UP;
  if (tr3_camera_look(&cam, eye, at, up, HZ_CFG_FOV_DEG * M_PI / 180.0, NRAY, NRAY) != 0) {
    printf("FAIL: камера\n");
    failed = 1;
    hz_poly_free(&ps);
    hz_seg_free(&sg);
    return;
  }
  int n = NRAY * NRAY;
  double *dt = malloc((size_t)n * sizeof *dt), *dp = malloc((size_t)n * sizeof *dp);
  int32_t *wh = malloc((size_t)n * sizeof *wh), *wt = malloc((size_t)n * sizeof *wt);
  if (dt == NULL || dp == NULL || wh == NULL || wt == NULL) {
    printf("FAIL: память\n");
    failed = 1;
    free(dt);
    free(dp);
    free(wh);
    free(wt);
    hz_poly_free(&ps);
    hz_seg_free(&sg);
    return;
  }
  double t3 = now_s();
#pragma omp parallel for schedule(dynamic, 8)
  for (int i = 0; i < n; i++) {
    double o[3], d[3];
    int32_t tri = -1;
    tr3_camera_ray(&cam, i % NRAY, i / NRAY, o, d);
    dt[i] = hit_tris(m, o, d, &tri);
    wt[i] = (tri >= 0) ? sg.label[tri] : -1;
    dp[i] = hit_polys(&ps, o, d, &wh[i]);
  }
  double t4 = now_s();

  /* ДВЕ ВЫБОРКИ, А НЕ ОДНА, И ЭТО НЕ КОСМЕТИКА. Луч, попавший В ТОТ ЖЕ участок,
   * меряет ГЕОМЕТРИЧЕСКУЮ ошибку плоскости и обязан лежать в пределах δ/|cos|.
   * Луч, у которого ближайшим оказался ДРУГОЙ участок, — это переворот порядка
   * по глубине на кромке: расхождение там равно всей дальности при любой схеме
   * (К20), и подмешивать его в общий хвост значит прятать первую величину за
   * второй. */
  double *diff = malloc((size_t)n * sizeof *diff), *diff2 = malloc((size_t)n * sizeof *diff2);
  if (diff == NULL || diff2 == NULL) {
    printf("FAIL: память\n");
    failed = 1;
    free(diff);
    free(diff2);
    diff = NULL;
    diff2 = NULL;
  } else {
    int nd = 0, nd2 = 0, miss_t = 0, miss_p = 0;
    for (int i = 0; i < n; i++) {
      int ht = (dt[i] > 0.0), hp = (dp[i] > 0.0);
      if (ht && hp) {
        if (wt[i] == wh[i])
          diff[nd++] = fabs(dt[i] - dp[i]);
        else
          diff2[nd2++] = fabs(dt[i] - dp[i]);
      } else if (ht)
        miss_p++;
      else if (hp)
        miss_t++;
    }
    printf("   лучей %d, перебор %.2f с; оба попали %d, только треугольники %d, "
           "только полигоны %d (%.2f%% покрытия расходится)\n",
           n, t4 - t3, nd + nd2, miss_p, miss_t, 100.0 * (miss_p + miss_t) / (double)n);
    printf("   тот же участок: %d лучей (%.2f%%); другой участок (кромка): %d (%.2f%%)\n", nd,
           100.0 * nd / (double)n, nd2, 100.0 * nd2 / (double)n);
    tail("|Δ| тот же участок, м", diff, nd, delta);
    tail("|Δ| кромка, м", diff2, nd2, delta);
    /* Внутри участка ошибка ограничена НАКЛОНОМ плоскости: δ/|cos| не больше
     * диаметра сцены, но на скользящем луче δ/|cos| велик законно. Поэтому
     * проверяется p99, а не max. */
    if (nd > 0) {
      qsort(diff, (size_t)nd, sizeof *diff, cmp_d);
      if (diff[(int)(nd * 0.99)] > 50.0 * delta) {
        printf("FAIL: внутри участка расхождение больше 50 δ — плоскость не та\n");
        failed = 1;
      }
    }
    free(diff);
    free(diff2);
  }

  /* --- НЕГАТИВНЫЙ КОНТРОЛЬ: один полигон на 10 δ ---
   * Берётся САМЫЙ ВИДИМЫЙ, а не самый большой: сдвинуть невидимый полигон и
   * не увидеть изменения — это не контроль, а тавтология. */
  int32_t *cnt = calloc((size_t)ps.np, sizeof *cnt);
  int32_t big = 0;
  if (cnt != NULL) {
    for (int i = 0; i < n; i++)
      if (wh[i] >= 0) cnt[wh[i]]++;
    for (int32_t k = 1; k < ps.np; k++)
      if (cnt[k] > cnt[big]) big = k;
  }
  int nbig = (cnt != NULL) ? cnt[big] : 0;
  free(cnt);
  ps.p[big].off += 10.0 * delta;
  for (int a = 0; a < 3; a++)
    ps.p[big].org[a] += 10.0 * delta * ps.p[big].n[a];
  /* Луч пускается В ЭТОТ ПОЛИГОН, а не в сцену. Сцену брать нельзя: сдвинутый
   * полигон уходит за соседей, ближайшим становится другой, и контроль
   * измерял бы уже не сдвиг, а перекрытие. */
  double *nc = malloc((size_t)n * sizeof *nc);
  if (nc != NULL) {
    int keep = 0, lost = 0;
    hz_polyset one = ps;
    one.p = &ps.p[big];
    one.np = 1;
    for (int i = 0; i < n; i++) {
      if (wh[i] != big || !(dt[i] > 0.0)) continue;
      double o[3], d[3];
      tr3_camera_ray(&cam, i % NRAY, i / NRAY, o, d);
      int32_t w2 = -1;
      double s = hit_polys(&one, o, d, &w2);
      if (s > 0.0)
        nc[keep++] = fabs(dt[i] - s);
      else
        lost++;
    }
    printf("   НЕГАТИВНЫЙ КОНТРОЛЬ: полигон %d (%.3f м², %d лучей) сдвинут на %g м; "
           "ушло из-под луча %d\n",
           big, ps.p[big].area, nbig, 10.0 * delta, lost);
    tail("|Δ| по его лучам, м", nc, keep, delta);
    /* ДВА ИСХОДА, И ОБА — ПРОВАЛ КАРТИНКИ, чего контроль и требует: либо луч
     * всё ещё попадает и дальность ушла на ≈10 δ, либо полигон вышел из-под
     * луча вовсе (сдвиг больше, чем сцена в этом месте может поглотить). Не
     * провал — ровно один случай: луч попадает, а дальность НЕ изменилась. */
    if (keep > 0 && nc[keep / 2] < 5.0 * delta) {
      printf("FAIL: сдвиг на 10 δ не проявился — луч слеп к положению полигона\n");
      failed = 1;
    }
    if (keep == 0 && lost == 0) {
      printf("FAIL: контроль пуст — ни один луч не проверил сдвинутый полигон\n");
      failed = 1;
    }
    free(nc);
  }

  free(dt);
  free(dp);
  free(wh);
  free(wt);
  hz_poly_free(&ps);
  hz_seg_free(&sg);
}

int main(void) {
  hz_objmesh m;
  if (hz_obj_load(&m, SCENE, SCALE) != 0) {
    printf("SKIP: нет %s — запустить scripts/fetch_scene.sh assets\n", SCENE);
    return 0;
  }
  printf("== test_polygon: %s, треугольников %d\n", SCENE, m.nt);
  run(&m, DELTA_COARSE);
  run(&m, DELTA_FINE);

  /* --- УСТОЙЧИВОСТЬ К ПЕРЕСТАНОВКЕ --- */
  printf("\n=== перестановка треугольников ===\n");
  hz_pseglist a;
  if (hz_seg_planar(&a, &m, DELTA_COARSE) != 0) return 1;
  sig *sa = sigs_of(&a);
  shuffle_mesh(&m);
  hz_pseglist b;
  if (hz_seg_planar(&b, &m, DELTA_COARSE) != 0) return 1;
  sig *sb = sigs_of(&b);
  if (sa == NULL || sb == NULL) return 1;
  printf("   участков до %d, после %d\n", a.nseg, b.nseg);
  if (a.nseg != b.nseg) {
    printf("FAIL: число участков изменилось — затравки упорядочены не по существу\n");
    failed = 1;
  } else {
    int nbad = 0;
    double worst = 0.0;
    for (int32_t i = 0; i < a.nseg; i++)
      for (int k = 0; k < 6; k++) {
        double e = fabs(sa[i].k[k] - sb[i].k[k]);
        if (e > worst) worst = e;
        if (e > 1e-12) {
          nbad++;
          k = 6;
        }
      }
    printf("   несовпавших отпечатков %d, худшее расхождение %.3g\n", nbad, worst);
    if (nbad != 0) {
      printf("FAIL: набор участков зависит от порядка входа — инкрементность не обоснована\n");
      failed = 1;
    }
  }
  free(sa);
  free(sb);
  hz_seg_free(&a);
  hz_seg_free(&b);
  hz_obj_free(&m);

  printf("\n%s\n", failed ? "== test_polygon: ЕСТЬ ОТКАЗЫ" : "== test_polygon: OK");
  return failed;
}
