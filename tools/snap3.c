/* snap3 — СНИМОК: пирамида шагов переноса на синтетике (§833).
 *
 * Новый инструмент новой структуры (STRUCTURE.md §9.1): загружает OBJ
 * существующим разборщиком (`src/scene_obj.c`, Ш2), строит пирамиду
 * (`src/nstruct/pyr.c`), печатает бюджет байтов, распределение кусков по
 * клеткам и детекторы согласованности. СВИПА НЕТ, поле — нулевой слот.
 *
 * Ключи:
 *   lev=N     клеток по максимальной оси = 2^N (умолчание 6 → 64);
 *   screwb    НК: скрестить каждый 97-й кусок (клетка +1, А1454-заворот);
 *   scrfew    НК: скрестить ровно один кусок (калибровка детектора);
 *   scale=F   масштаб сцены (умолчание 1.0).
 *
 * Чистый прогон обязан напечатать все детекторы нулями и ΔΣ посимвольным
 * нулём; прогон с НК-ключом обязан показать предсказанный провал.
 */
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "nstruct/pyr.h"
#include "scene_obj.h"

static double now_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int cmp_i32(const void *a, const void *b) {
  int32_t x = *(const int32_t *)a, y = *(const int32_t *)b;
  return (x > y) - (x < y);
}

int main(int argc, char **argv) {
  const char *path = NULL;
  double scale = 1.0, cell, t0, t1;
  int lev = 6; /* клеток по максимальной оси = 64: мелкая синтетика */
  int screwb = 0, scrfew = 0, march = 0, marchinv = 0, aggrat = 0;
  int i, ax;
  hz_objmesh m;
  hz_pyr py;
  hz_pyr_verdict vd;
  double *cmin = NULL, *cmax = NULL, *cent = NULL;
  int32_t *mtl = NULL;
  double sum_tri = 0.0, sum_pcs = 0.0;
  int64_t total, npcs_bytes, csr_bytes, leaf_bytes, lev_bytes;
  int32_t *dist = NULL, ndist = 0;
  int64_t nlevnodes = 0;

  for (i = 1; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0) {
      lev = atoi(argv[i] + 4);
    } else if (strcmp(argv[i], "screwb") == 0) {
      screwb = 1;
    } else if (strcmp(argv[i], "scrfew") == 0) {
      scrfew = 1;
    } else if (strcmp(argv[i], "march") == 0) {
      march = 1;
    } else if (strcmp(argv[i], "marchinv") == 0) {
      march = 1;
      marchinv = 1;
    } else if (strncmp(argv[i], "aggrat=", 7) == 0) {
      aggrat = atoi(argv[i] + 7);
    } else if (strncmp(argv[i], "scale=", 6) == 0) {
      scale = atof(argv[i] + 6);
    } else {
      path = argv[i];
    }
  }
  if (!path) {
    fprintf(
        stderr,
        "use: snap3 <scene.obj> [lev=N] [screwb|scrfew] [march|marchinv] [aggrat=N] [scale=F]\n");
    return 2;
  }
  if (hz_obj_load(&m, path, scale) != 0) {
    fprintf(stderr, "snap3: не читается %s\n", path);
    return 2;
  }

  /* геометрия для биннинга: bbox и центроид каждого треугольника */
  cmin = (double *)malloc((size_t)m.nt * 3 * sizeof *cmin);
  cmax = (double *)malloc((size_t)m.nt * 3 * sizeof *cmax);
  cent = (double *)malloc((size_t)m.nt * 3 * sizeof *cent);
  mtl = (int32_t *)malloc((size_t)m.nt * sizeof *mtl);
  if (!cmin || !cmax || !cent || !mtl) {
    fprintf(stderr, "snap3: нет памяти\n");
    return 2;
  }
  for (i = 0; i < m.nt; i++) {
    double p[3][3];
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
    mtl[i] = m.fm[i];
    sum_tri += hz_obj_tri_area(&m, (int32_t)i);
  }

  {
    double maxdim = 0.0;
    for (ax = 0; ax < 3; ax++) {
      double s = m.hi[ax] - m.lo[ax];
      if (s > maxdim) maxdim = s;
    }
    cell = maxdim / (double)(1 << lev);
  }

  printf("== snap3 %s\n", path);
  printf("габарит %.6g x %.6g x %.6g м; V=%d nt=%d ndegen=%" PRId64 " nquad=%" PRId64 "\n",
         m.hi[0] - m.lo[0], m.hi[1] - m.lo[1], m.hi[2] - m.lo[2], m.nv, m.nt, m.ndegen, m.nquad);
  printf("lev=%d клетка %.6g м\n", lev, cell);

  t0 = now_sec();
  if (hz_pyr_build(&py, m.nt, cmin, cmax, cent, mtl, m.lo, m.hi, cell) != 0) {
    fprintf(stderr, "snap3: построение не удалось\n");
    return 2;
  }
  t1 = now_sec();
  printf("сетка листа %dx%dx%d = %" PRId64 " клеток\n", py.nx, py.ny, py.nz, py.ncells);

  for (i = 0; i < py.nlev; i++)
    nlevnodes += py.nlev_nodes[i];

  if (screwb || scrfew) {
    int32_t ns = hz_pyr_screw(&py, scrfew ? 1 : 97);
    printf("НК: скрещено кусков %d (шаг %s)\n", ns, scrfew ? "1 (scrfew)" : "97 (screwb)");
  }

  hz_pyr_verify(&py, m.nt, cmin, cmax, cent, &vd);

  /* ΔΣ: та же сумма в том же порядке по кусочной биекции (А1456) */
  for (i = 0; i < m.nt; i++) {
    if (py.pcs[i].tri != i) break; /* биекция нарушена — ниже печатаем счётчик */
    sum_pcs += hz_obj_tri_area(&m, py.pcs[i].tri);
  }
  {
    int bijection = 1;
    for (i = 0; i < m.nt; i++)
      if (py.pcs[i].tri != i) bijection = 0;

    /* распределение кусков на клетку (по клеткам с кусками) */
    dist = (int32_t *)malloc((size_t)(py.ncentleaf > 0 ? py.ncentleaf : 1) * sizeof *dist);
    if (dist) {
      int32_t li;
      for (li = 0; li < py.nleaf; li++)
        if (py.leaf[li].npcs > 0) dist[ndist++] = py.leaf[li].npcs;
      qsort(dist, (size_t)ndist, sizeof *dist, cmp_i32);
    }

    npcs_bytes = (int64_t)m.nt * (int64_t)sizeof(hz_piece);
    csr_bytes = (int64_t)m.nt * (int64_t)sizeof(int32_t);
    leaf_bytes = (int64_t)py.nleaf * (int64_t)sizeof(hz_pyr_leaf);
    lev_bytes = nlevnodes * (int64_t)sizeof(hz_pyr_node);
    total = npcs_bytes + csr_bytes + leaf_bytes + lev_bytes;

    printf("листья: bbox-занятых %d (пометок %" PRId64 "), с кусками %" PRId64 ", пустых %" PRId64
           "\n",
           py.nleaf, (int64_t)py.nmarks, (int64_t)py.ncentleaf, (int64_t)(py.ncells - py.nleaf));
    printf("внутренних узлов %" PRId64 " на %d уровнях; агрегатов 0 (§833)\n", nlevnodes, py.nlev);
    printf("состояния: ДРОБЛЁН %d листьев + %" PRId64 " внутренних, АГРЕГАТ 0\n", py.nleaf,
           nlevnodes);
    if (ndist > 0) {
      printf("кусков на клетку (с кусками): медиана %d p99 %d max %d\n", dist[ndist / 2],
             dist[(99 * ndist) / 100], dist[ndist - 1]);
    }
    printf("бюджет: куски %" PRId64 " + CSR %" PRId64 " + листья %" PRId64 " + уровни %" PRId64
           " = %" PRId64 " Б; %.2f Б/кусок (лист %.2f, ур. %.2f)\n",
           npcs_bytes, csr_bytes, leaf_bytes, lev_bytes, total, (double)total / (double)m.nt,
           (double)leaf_bytes / (double)m.nt, (double)lev_bytes / (double)m.nt);
    printf("площадь: Σ кусков %.17g против Σ треугольников %.17g; биекция %s\n", sum_pcs, sum_tri,
           bijection ? "полная" : "НАРУШЕНА");
    printf("детекторы: (а) клетка %+" PRId64 " csr %+" PRId64 " (б) bbox %+" PRId64
           " (в) пустая %+" PRId64 "\n",
           vd.d_cell, vd.d_csr, vd.d_bbox, vd.d_empty);
    printf("время построения %.3f с\n", t1 - t0);
    if ((screwb || scrfew) &&
        (vd.d_cell == 0 && vd.d_csr == 0 && vd.d_bbox == 0 && vd.d_empty == 0))
      printf("НК НЕ СРАБОТАЛ — детекторы слепы, это провал контроля\n");
    free(dist);
  }

  /* ---- МАРШ (§834): 14 направлений = 6 осей + 8 диагоналей ---- */
  if (march) {
    static const double dirs[14][3] = {
        {1, 0, 0},  {-1, 0, 0}, {0, 1, 0},  {0, -1, 0},  {0, 0, 1},   {0, 0, -1},  {1, 1, 1},
        {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}, {-1, -1, -1}};
    const double sq3 = 0.57735026918962573; /* 1/√3, диагонали единичные */
    static const char *names[14] = {"+x",  "-x",  "+y",  "-y",  "+z",  "-z",  "+++",
                                    "++-", "+-+", "-++", "+--", "-+-", "--+", "---"};
    uint64_t *bits = (uint64_t *)malloc((size_t)((py.nleaf + 63) >> 6) * sizeof *bits);
    int32_t *vindex = (int32_t *)malloc((size_t)py.nleaf * sizeof *vindex);
    int64_t tot_l = 0, tot_n = 0, tot_inv = 0, tot_a = 0, tot_dep = 0;
    int d;
    if (aggrat > 0) {
      hz_pyr_mark_aggr(&py, (int32_t)aggrat);
      printf("МАРШ: АГРЕГАТ помечен (aggrat=%d)\n", aggrat);
    }
    for (d = 0; d < 14; d++) {
      double om[3] = {dirs[d][0], dirs[d][1], dirs[d][2]};
      hz_pyr_march_stat st;
      int miss;
      int64_t dep = 0; /* нарушения порядка СОСЕДНИХ ПО ГРАНИ клеток вдоль ω (А1468) */
      if (d >= 6) {
        om[0] *= sq3;
        om[1] *= sq3;
        om[2] *= sq3;
      }
      hz_pyr_march(&py, om, marchinv, &st, bits, vindex);
      miss = 0;
      if (bits) {
        int32_t li;
        for (li = 0; li < py.nleaf; li++)
          if (!((bits[li >> 6] >> (li & 63)) & 1)) miss++;
      }
      /* зависимые пары: ось a с ω_a ≠ 0; A=(i) → B=(i+e_a) низовой, обязана
       * быть ПОЗЖЕ. Проверяем только осевые направления — там пары полные. */
      if (d < 6) {
        int axis = d / 2, sign = (d % 2 == 0) ? 1 : -1;
        int64_t stride = axis == 0 ? 1 : (axis == 1 ? py.nx : (int64_t)py.nx * py.ny);
        int32_t li;
        for (li = 0; li < py.nleaf; li++) {
          int64_t id = py.leaf_id[li], nb;
          int64_t coord = (axis == 0)   ? id % py.nx
                          : (axis == 1) ? (id / py.nx) % py.ny
                                        : id / ((int64_t)py.nx * py.ny);
          nb = id + (sign > 0 ? stride : -stride);
          if (sign > 0 ? coord + 1 < (axis == 0   ? (int64_t)py.nx
                                      : axis == 1 ? (int64_t)py.ny
                                                  : (int64_t)py.nz)
                       : coord > 0) {
            int32_t nli = -1;
            {
              int32_t lo = 0, hi = py.nleaf;
              while (lo < hi) {
                int32_t mid = lo + (hi - lo) / 2;
                if (py.leaf_id[mid] < nb)
                  lo = mid + 1;
                else
                  hi = mid;
              }
              if (lo < py.nleaf && py.leaf_id[lo] == nb) nli = lo;
            }
            if (nli >= 0 && vindex[nli] < vindex[li]) dep++;
          }
        }
      }
      printf("  %s: листов %" PRId64 " (непосещённых %d) узлов %" PRId64 " инверсий %" PRId64
             " событий %" PRId64 " зависимых-пар-нарушено %" PRId64 "\n",
             names[d], st.leaves, miss, st.nodes, st.inversions, st.aggr, dep);
      tot_l += st.leaves;
      tot_n += st.nodes;
      tot_inv += st.inversions;
      tot_a += st.aggr;
      tot_dep += dep;
    }
    printf("МАРШ суммарно: на направление листов %" PRId64 ", цепь %" PRId64
           " (цепь/листья %.2f), инверсий %" PRId64 ", событий %" PRId64 ", зависимых-пар %" PRId64
           "\n",
           tot_l / 14, tot_n / 14, (double)(tot_n / 14) / (double)(tot_l / 14), tot_inv, tot_a,
           tot_dep);
    if (aggrat > 0) hz_pyr_unmark_aggr(&py);
    free(vindex);
    free(bits);
  }

  hz_pyr_free(&py);
  free(cmin);
  free(cmax);
  free(cent);
  free(mtl);
  hz_obj_free(&m);
  return 0;
}
