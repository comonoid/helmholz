/* pconv — Ш5: СКОЛЬКО НАПРАВЛЕНИЙ И СКОЛЬКО ОТСКОКОВ (PLAN_ELEMENTS.md).
 *
 * Три замера, и первый из них НЕ ТРОГАЕТ ПЕРЕНОС ВОВСЕ — так дешевле и так
 * честнее: Ш4 показал, что вся невязка печи угловая, значит `ND(допуск)` в
 * первом приближении есть свойство НАБОРА НАПРАВЛЕНИЙ, а не схемы.
 *
 *   A. КВАДРАТУРА: `ε(nmu, nphi) = max_n |S(n)/π − 1|`, где
 *      `S(n) = Σ_{ω·n<0} w_ω|ω·n|`. Оси `nmu` и `nphi` разделены НАМЕРЕННО:
 *      §2 задаёт свип по суммарному `ND`, а Ш4 нашёл, что `ND = 64` хуже
 *      `ND = 32`, потому что учетверение ушло не в ту ось. Суммарный `ND`
 *      этого не покажет никогда.
 *      Нормали берутся не только осевые: осевая нормаль — САМЫЙ ПЛОХОЙ случай
 *      для этого набора, и мерить только по ней значило бы завышать.
 *   B. ЗАЛ: полная развёртка при разных `ND` против самого мелкого набора.
 *      Метрика на РАДИАНСЕ (§4), взвешенная по площади, хвост а не среднее.
 *   C. ОТСКОКИ: `n(ρ, допуск)` на печи и на зале.
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

static double quad_S(const tr3_dirs *d, const double n[3]) {
  double s = 0.0;
  for (int m = 0; m < d->n; m++) {
    double c = d->ox[m] * n[0] + d->oy[m] * n[1] + d->oz[m] * n[2];
    if (c < 0.0) s += d->w[m] * (-c);
  }
  return s;
}

/* Худшая и среднеквадратичная ошибка по МНОГИМ нормалям: осевые — крайний
 * случай, и по ним одним судить нельзя. Направления берутся детерминированной
 * спиралью Фибоначчи, а не случайно: замер обязан быть воспроизводим. */
static void quad_err(const tr3_dirs *d, int nsample, double *worst, double *rms, double *axial) {
  double w = 0.0, s2 = 0.0;
  for (int i = 0; i < nsample; i++) {
    double z = 1.0 - 2.0 * ((double)i + 0.5) / (double)nsample;
    double r = sqrt(1.0 - z * z), phi = 2.399963229728653 * (double)i;
    double n[3] = {r * cos(phi), r * sin(phi), z};
    double e = quad_S(d, n) / M_PI - 1.0;
    if (fabs(e) > w) w = fabs(e);
    s2 += e * e;
  }
  *worst = w;
  *rms = sqrt(s2 / nsample);
  const double ax[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  double a = 0.0;
  for (int i = 0; i < 3; i++) {
    double e = fabs(quad_S(d, ax[i]) / M_PI - 1.0);
    if (e > a) a = e;
  }
  *axial = a;
}

/* ДВА ИСТОЧНИКА, И РАЗНИЦА МЕЖДУ НИМИ — ГЛАВНЫЙ ЗАМЕР Ш5.
 *
 * `compact = 1` — восемь светильников §2 (0.3×0.3 м, угловой размер с пола
 * `α ≈ 0.15` рад). Пропущенные ЧЕРЕЗ ОРДИНАТЫ, они дают классический лучевой
 * эффект: поверхность видит источник только вдоль тех немногих направлений,
 * которые в него попали, и между ними — полосы (К3 дословно). §1.5 прямо
 * запрещает так делать («прямой свет: аналитически, вне ординат»), и Ш5 меряет,
 * ЧЕГО этот запрет стоит, а не предполагает.
 *
 * `compact = 0` — ОДИН протяжённый источник во весь потолок, той же полной
 * мощности. Угловой размер `α ≈ 1`, лучевого эффекта нет по построению, и
 * развёртка меряет ровно то, ради чего она есть, — КОСВЕННЫЙ перенос
 * (§1.8 D: «у вторичного освещения источником служит ПРОТЯЖЁННАЯ поверхность»). */
static int add_lamps(hz_polyset *ps, const hz_objmesh *m, int32_t *first, int compact,
                     double *le_out) {
  const double n[3] = {0.0, -1.0, 0.0}, eu[3] = {1.0, 0.0, 0.0};
  *first = ps->np;
  double power = (double)(HZ_CFG_LAMP_NX * HZ_CFG_LAMP_NZ) * HZ_CFG_LAMP_SIDE * HZ_CFG_LAMP_SIDE *
                 HZ_CFG_LAMP_LE;
  if (compact) {
    for (int i = 0; i < HZ_CFG_LAMP_NX; i++)
      for (int j = 0; j < HZ_CFG_LAMP_NZ; j++) {
        double c[3];
        c[0] = m->lo[0] + ((double)i + 0.5) / HZ_CFG_LAMP_NX * (m->hi[0] - m->lo[0]);
        c[1] = HZ_CFG_LAMP_Y;
        c[2] = m->lo[2] + ((double)j + 0.5) / HZ_CFG_LAMP_NZ * (m->hi[2] - m->lo[2]);
        if (hz_poly_add_quad(ps, c, n, eu, HZ_CFG_LAMP_SIDE / 2.0, HZ_CFG_LAMP_SIDE / 2.0, 0) != 0)
          return 2;
      }
    *le_out = HZ_CFG_LAMP_LE;
  } else {
    double hx = 0.5 * (m->hi[0] - m->lo[0]) * 0.98, hz = 0.5 * (m->hi[2] - m->lo[2]) * 0.98;
    double c[3] = {0.5 * (m->lo[0] + m->hi[0]), HZ_CFG_LAMP_Y, 0.5 * (m->lo[2] + m->hi[2])};
    if (hz_poly_add_quad(ps, c, n, eu, hx, hz, 0) != 0) return 2;
    *le_out = power / (4.0 * hx * hz); /* та же полная мощность */
  }
  return 0;
}

typedef struct {
  double v;
  double a;
} wpair;

static int cmp_w(const void *x, const void *y) {
  const wpair *a = x, *b = y;
  return (a->v < b->v) ? -1 : ((a->v > b->v) ? 1 : 0);
}

/* Взвешенный по площади квантиль. Хвост, а не среднее (§4). */
static double wquant(wpair *v, int n, double q) {
  if (n <= 0) return 0.0;
  qsort(v, (size_t)n, sizeof *v, cmp_w);
  double tot = 0.0;
  for (int i = 0; i < n; i++)
    tot += v[i].a;
  double acc = 0.0;
  for (int i = 0; i < n; i++) {
    acc += v[i].a;
    if (acc >= q * tot) return v[i].v;
  }
  return v[n - 1].v;
}

int main(int argc, char **argv) {
  double delta = (argc > 1) ? strtod(argv[1], NULL) : 0.045;
  double h = (argc > 2) ? strtod(argv[2], NULL) : 0.02;

  /* --- A. КВАДРАТУРА: оси разделены -------------------------------------- */
  printf("== A. квадратура сама по себе: ε = max|S(n)/π − 1| по 4096 нормалям\n");
  printf("   %-5s %-5s %-6s %12s %12s %12s\n", "nmu", "nphi", "ND", "худшая", "СКО", "по осям");
  const int mus[6] = {1, 2, 4, 8, 16, 32}, phis[4] = {1, 2, 4, 8};
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 4; j++) {
      tr3_dirs d;
      if (tr3_dirs_product(&d, mus[i], phis[j]) != 0) return 1;
      double w, r, a;
      quad_err(&d, 4096, &w, &r, &a);
      printf("   %-5d %-5d %-6d %11.4f%% %11.4f%% %11.4f%%\n", mus[i], phis[j], d.n, 100.0 * w,
             100.0 * r, 100.0 * a);
      tr3_dirs_free(&d);
    }
  printf("\n   При РАВНОМ ND сравнить (nmu, nphi):\n");
  const int pa[4][2] = {{4, 4}, {8, 2}, {16, 1}, {2, 8}};
  for (int i = 0; i < 4; i++) {
    tr3_dirs d;
    if (tr3_dirs_product(&d, pa[i][0], pa[i][1]) != 0) return 1;
    double w, r, a;
    quad_err(&d, 4096, &w, &r, &a);
    printf("   (%2d, %2d) ND = %3d: худшая %8.4f%%, СКО %8.4f%%, по осям %8.4f%%\n", pa[i][0],
           pa[i][1], d.n, 100.0 * w, 100.0 * r, 100.0 * a);
    tr3_dirs_free(&d);
  }

  /* --- B/C. ЗАЛ ----------------------------------------------------------- */
  hz_objmesh m;
  if (hz_obj_load(&m, HZ_CFG_HALL_OBJ, HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет %s\n", HZ_CFG_HALL_OBJ);
    return 1;
  }
  hz_pseglist sg;
  if (hz_seg_planar(&sg, &m, delta) != 0) return 1;
  hz_polyset ps;
  if (hz_poly_build(&ps, &m, &sg) != 0) return 1;
  int compact = (argc > 3) ? atoi(argv[3]) : 1;
  int32_t lamp0 = 0;
  double lampLe = HZ_CFG_LAMP_LE;
  if (add_lamps(&ps, &m, &lamp0, compact, &lampLe) != 0) return 1;

  /* ЭТАЛОН ВЗЯТ СИЛЬНО МЕЛЬЧЕ СВИПА, и это не запас, а необходимость: первый
   * прогон брал эталоном `ND = 256`, а два РАЗНЫХ набора того же `ND = 256`
   * разошлись между собой на `p99 = 48%`. Эталон, который спорит сам с собой,
   * ничего не меряет. */
  const int NS = 8;
  const int smu[8] = {1, 2, 2, 4, 8, 4, 8, 8}, sphi[8] = {1, 2, 4, 4, 2, 8, 4, 16};
  double *Lref = NULL;
  wpair *buf = malloc((size_t)ps.np * sizeof *buf);
  if (buf == NULL) return 1;

  printf("\n== B/C. зал: δ = %g м, h = %g м, полигонов %d, источник %s (L_e = %g)\n", delta, h,
         ps.np, compact ? "КОМПАКТНЫЙ (8 светильников)" : "ПРОТЯЖЁННЫЙ (потолок)", lampLe);
  printf("   %-5s %-5s %-6s %6s %9s %10s %10s %9s %9s %9s\n", "nmu", "nphi", "ND", "отск",
         "время,с", "p50 отн.", "p99 отн.", "S>1%", "S>10%", "Φin/Φout");
  for (int s = NS - 1; s >= 0; s--) { /* самый мелкий набор первым — он эталон */
    hz_ptrans t;
    if (hz_ptrans_init(&t, &ps, &m) != 0) return 1;
    for (int32_t k = lamp0; k < ps.np; k++) {
      t.Le[k] = lampLe;
      t.rho[k] = 0.0;
    }
    tr3_dirs d;
    if (tr3_dirs_product(&d, smu[s], sphi[s]) != 0) return 1;
    hz_pstats st;
    double ta = now_s();
    int nb = 0;
    for (; nb < 300; nb++) {
      if (hz_psweep_bounce(&t, &d, h, 0, 0.0, HZ_LAYOUT_RUNS, &st) != 0) return 1;
      if (st.dE < 1e-10) break;
    }
    double tb = now_s();
    /* Радианс в центре полигона: L_out = L_e + ρ·E/π (§1.2, альбедо при чтении). */
    double *L = malloc((size_t)ps.np * sizeof *L);
    if (L == NULL) return 1;
    for (int32_t k = 0; k < ps.np; k++)
      L[k] = hz_ptrans_lout(&t, k, 0.0, 0.0);
    if (s == NS - 1) {
      Lref = L;
      printf("   %-5d %-5d %-6d %6d %9.2f %10s %10s %9s %9s %9.5f  ЭТАЛОН\n", smu[s], sphi[s], d.n,
             nb + 1, tb - ta, "-", "-", "-", "-", st.phi_in / st.phi_out);
    } else {
      int n = 0;
      double scale = 0.0;
      for (int32_t k = 0; k < ps.np; k++)
        if (Lref[k] > 0.0) scale += Lref[k] * ps.p[k].area;
      double atot = 0.0;
      for (int32_t k = 0; k < ps.np; k++)
        atot += ps.p[k].area;
      scale /= atot;
      for (int32_t k = 0; k < ps.np; k++) {
        /* Нормировка на СРЕДНИЙ радианс сцены, а не на местный: деление на
         * почти нулевой радианс в тени даёт бесконечность там, где глазу
         * ничего не видно (К19/К20 тот же довод). */
        buf[n].v = fabs(L[k] - Lref[k]) / scale;
        buf[n].a = ps.p[k].area;
        n++;
      }
      /* ДОЛЯ ПЛОЩАДИ ЗА ПОРОГОМ понятнее квантиля: «p99 = 0.48» не говорит,
       * пол-сцены это или три щепки, а «0.4% площади хуже 10%» — говорит. */
      double s1 = 0.0, s10 = 0.0, atot2 = 0.0;
      for (int i = 0; i < n; i++) {
        atot2 += buf[i].a;
        if (buf[i].v > 0.01) s1 += buf[i].a;
        if (buf[i].v > 0.10) s10 += buf[i].a;
      }
      printf("   %-5d %-5d %-6d %6d %9.2f %10.3e %10.3e %8.2f%% %8.2f%% %9.5f\n", smu[s], sphi[s],
             d.n, nb + 1, tb - ta, wquant(buf, n, 0.5), wquant(buf, n, 0.99), 100.0 * s1 / atot2,
             100.0 * s10 / atot2, st.phi_in / st.phi_out);
      free(L);
    }
    tr3_dirs_free(&d);
    hz_ptrans_free(&t);
  }
  free(Lref);

  /* --- C. отскоки против допуска ------------------------------------------ */
  printf("\n== C. отскоки: зал, ND = 128, против допуска по dE\n");
  printf("   %-10s %8s %14s\n", "допуск", "отскоков", "ln(доп)/ln(ρ_эфф)");
  {
    hz_ptrans t;
    if (hz_ptrans_init(&t, &ps, &m) != 0) return 1;
    for (int32_t k = lamp0; k < ps.np; k++) {
      t.Le[k] = lampLe;
      t.rho[k] = 0.0;
    }
    double arho = 0.0, atot = 0.0;
    for (int32_t k = 0; k < ps.np; k++) {
      arho += ps.p[k].area * t.rho[k];
      atot += ps.p[k].area;
    }
    double rho = arho / atot;
    tr3_dirs d;
    if (tr3_dirs_product(&d, 4, 4) != 0) return 1;
    const double tols[5] = {1e-1, 1e-2, 1e-3, 1e-5, 1e-8};
    int idx = 0, nb = 0;
    hz_pstats st;
    for (; nb < 300 && idx < 5; nb++) {
      if (hz_psweep_bounce(&t, &d, h, 0, 0.0, HZ_LAYOUT_RUNS, &st) != 0) return 1;
      while (idx < 5 && st.dE < tols[idx]) {
        printf("   %-10.0e %8d %14.1f\n", tols[idx], nb + 1, log(tols[idx]) / log(rho));
        idx++;
      }
    }
    printf("   (средневзвешенное альбедо сцены ρ_эфф = %.4f)\n", rho);
    tr3_dirs_free(&d);
    hz_ptrans_free(&t);
  }

  free(buf);
  hz_poly_free(&ps);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
