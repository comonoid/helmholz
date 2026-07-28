/* ptex — Ш8: СКОЛЬКО РАЗРЕЗОВ ТРЕБУЕТ ПЕРЕНОС ПРИ ВКЛЮЧЁННОМ Т1.
 *
 * PLAN_ELEMENTS.md, Ш8, правила Т1–Т3 (§1.9).
 *
 * Т1 ГОВОРИТ: хранится ОБЛУЧЁННОСТЬ, альбедо применяется при чтении
 * (`L_out = ρ·E`), поэтому текстура НЕ РЕЖЕТ полигоны ради изображения — вся
 * мелкая деталь сидит в `ρ` и попадает в кадр даром. Резать надо ТОЛЬКО ради
 * ПЕРЕНОСА, и только там, где подгонка `⟨ρ⟩ + градиент` уже не описывает
 * альбедо: тогда полигон отдаёт в сцену не тот поток.
 *
 * ЧТО ИМЕННО МЕРЯЕТСЯ. Для каждого полигона альбедо его треугольников
 * подгоняется двумя способами — ПОСТОЯННЫМ `⟨ρ⟩` и ЛИНЕЙНЫМ `a + b·u + c·v`, —
 * и считается взвешенная по площади невязка. Разрез нужен там, где невязка выше
 * допуска. Ответ Ш8 — сколько таких полигонов.
 *
 * ТЕКСТУРНЫХ КАРТ ЗДЕСЬ НЕТ, И ЭТО ОГОВОРКА, А НЕ УМОЛЧАНИЕ. В конференц-зале
 * их нет в `.mtl` вовсе (только `Kd`), у Rungholt они есть, но лежат в PNG, а
 * декодера в проекте нет и зависимости заводить нельзя. Поэтому мерится
 * вариация альбедо на уровне МАТЕРИАЛОВ. Это ровно та часть текстуры, о которой
 * говорит предсказание плана — «сильный контраст альбедо на КРУПНОМ полигоне
 * есть шов или кромка»; мелкая внутриматериальная деталь Т1 не режет по
 * построению, и мерить в ней нечего.
 */

#include "poly_seg.h"
#include "polygon.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

/* Гаусс 3×3 с частичным выбором. */
static int solve3(double A[3][3], double b[3]) {
  for (int c = 0; c < 3; c++) {
    int piv = c;
    for (int r = c + 1; r < 3; r++)
      if (fabs(A[r][c]) > fabs(A[piv][c])) piv = r;
    if (!(fabs(A[piv][c]) > 0.0)) return 1;
    if (piv != c) {
      for (int k = 0; k < 3; k++) {
        double t = A[c][k];
        A[c][k] = A[piv][k];
        A[piv][k] = t;
      }
      double t = b[c];
      b[c] = b[piv];
      b[piv] = t;
    }
    for (int r = 0; r < 3; r++) {
      if (r == c) continue;
      double f = A[r][c] / A[c][c];
      for (int k = 0; k < 3; k++)
        A[r][k] -= f * A[c][k];
      b[r] -= f * b[c];
    }
  }
  for (int c = 0; c < 3; c++)
    b[c] /= A[c][c];
  return 0;
}

int main(int argc, char **argv) {
  double delta = (argc > 1) ? strtod(argv[1], NULL) : 0.045;
  hz_objmesh m;
  if (hz_obj_load(&m, HZ_CFG_HALL_OBJ, HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет %s\n", HZ_CFG_HALL_OBJ);
    return 1;
  }
  printf("== Ш8: зал, материалов %d; текстурных карт в .mtl НЕТ (только Kd)\n", m.nmtl);
  printf("   альбедо материалов: ");
  double amin = 9, amax = -9;
  for (int32_t i = 0; i < m.nmtl; i++) {
    if (m.mtl[i].kd < amin) amin = m.mtl[i].kd;
    if (m.mtl[i].kd > amax) amax = m.mtl[i].kd;
  }
  printf("от %.3f до %.3f\n", amin, amax);

  const double dl[4] = {0.003, 0.01, 0.045, 0.15};
  printf("\n   %-7s %8s %9s %9s %11s %11s %11s %11s\n", "δ, м", "полиг", "на шве", "доля шва",
         "p50 конст", "p99 конст", "p50 лин", "p99 лин");
  for (int id = 0; id < 4; id++) {
    double d = (argc > 1) ? delta : dl[id];
    hz_pseglist sg;
    if (hz_seg_planar(&sg, &m, d) != 0) return 1;
    hz_polyset ps;
    if (hz_poly_build(&ps, &m, &sg) != 0) return 1;

    double *rc = malloc((size_t)ps.np * sizeof *rc);
    double *rl = malloc((size_t)ps.np * sizeof *rl);
    double *ar = malloc((size_t)ps.np * sizeof *ar);
    if (rc == NULL || rl == NULL || ar == NULL) return 1;
    int32_t nmix = 0;
    double amix = 0.0, atot = 0.0;

    for (int32_t k = 0; k < ps.np; k++) {
      const hz_poly *P = &ps.p[k];
      /* Взвешенная по площади подгонка альбедо треугольников. */
      double G[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}}, rhs[3] = {0, 0, 0};
      double sw = 0.0, sr = 0.0;
      int32_t m0 = (P->ntri > 0) ? m.fm[ps.tri[P->t0]] : 0;
      int mixed = 0;
      for (int32_t i = 0; i < P->ntri; i++) {
        int32_t t = ps.tri[P->t0 + i];
        double p[3][3], c[3] = {0, 0, 0};
        hz_obj_tri(&m, t, p);
        for (int j = 0; j < 3; j++)
          for (int a = 0; a < 3; a++)
            c[a] += p[j][a] / 3.0;
        double q[3];
        for (int a = 0; a < 3; a++)
          q[a] = c[a] - P->org[a];
        double u = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
        double v = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
        double w = hz_obj_tri_area(&m, t), rho = m.mtl[m.fm[t]].kd;
        if (m.fm[t] != m0) mixed = 1;
        double bas[3] = {1.0, u, v};
        for (int a = 0; a < 3; a++) {
          for (int b = 0; b < 3; b++)
            G[a][b] += w * bas[a] * bas[b];
          rhs[a] += w * bas[a] * rho;
        }
        sw += w;
        sr += w * rho;
      }
      ar[k] = sw;
      atot += sw;
      if (mixed) {
        nmix++;
        amix += sw;
      }
      double mean = (sw > 0.0) ? sr / sw : 0.0;
      double cc[3] = {rhs[0], rhs[1], rhs[2]};
      double GG[3][3];
      memcpy(GG, G, sizeof GG);
      int ok = (sw > 0.0) && (solve3(GG, cc) == 0);
      /* Невязки: взвешенная по площади СКО для постоянной и линейной подгонок. */
      double s2c = 0.0, s2l = 0.0;
      for (int32_t i = 0; i < P->ntri; i++) {
        int32_t t = ps.tri[P->t0 + i];
        double p[3][3], c[3] = {0, 0, 0};
        hz_obj_tri(&m, t, p);
        for (int j = 0; j < 3; j++)
          for (int a = 0; a < 3; a++)
            c[a] += p[j][a] / 3.0;
        double q[3];
        for (int a = 0; a < 3; a++)
          q[a] = c[a] - P->org[a];
        double u = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
        double v = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
        double w = hz_obj_tri_area(&m, t), rho = m.mtl[m.fm[t]].kd;
        double ec = rho - mean;
        double el = ok ? rho - (cc[0] + cc[1] * u + cc[2] * v) : ec;
        s2c += w * ec * ec;
        s2l += w * el * el;
      }
      rc[k] = (sw > 0.0) ? sqrt(s2c / sw) : 0.0;
      rl[k] = (sw > 0.0) ? sqrt(s2l / sw) : 0.0;
    }
    double *sc = malloc((size_t)ps.np * sizeof *sc);
    memcpy(sc, rc, (size_t)ps.np * sizeof *sc);
    qsort(sc, (size_t)ps.np, sizeof *sc, cmp_d);
    double *sl = malloc((size_t)ps.np * sizeof *sl);
    memcpy(sl, rl, (size_t)ps.np * sizeof *sl);
    qsort(sl, (size_t)ps.np, sizeof *sl, cmp_d);
    printf("   %-7g %8d %9d %8.2f%% %11.3e %11.3e %11.3e %11.3e\n", d, ps.np, nmix,
           100.0 * amix / atot, sc[ps.np / 2], sc[(int)(ps.np * 0.99)], sl[ps.np / 2],
           sl[(int)(ps.np * 0.99)]);

    /* СКОЛЬКО РАЗРЕЗОВ ТРЕБУЕТ ПЕРЕНОС: полигоны, где ЛИНЕЙНАЯ подгонка не
     * укладывается в допуск. Т1 снимает вопрос изображения, остаётся поток. */
    const double tols[3] = {0.02, 0.05, 0.10};
    for (int it = 0; it < 3; it++) {
      int32_t nc = 0;
      double aa = 0.0;
      for (int32_t k = 0; k < ps.np; k++)
        if (rl[k] > tols[it]) {
          nc++;
          aa += ar[k];
        }
      printf("        допуск %.2f: разрезов требует %d полигонов (%.2f%% от числа, "
             "%.2f%% площади)\n",
             tols[it], nc, 100.0 * nc / (double)ps.np, 100.0 * aa / atot);
    }
    free(sc);
    free(sl);
    free(rc);
    free(rl);
    free(ar);
    hz_poly_free(&ps);
    hz_seg_free(&sg);
    if (argc > 1) break;
  }
  hz_obj_free(&m);
  return 0;
}
