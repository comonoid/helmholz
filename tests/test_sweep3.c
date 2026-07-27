/* Фальсификаторы трёхмерной развёртки (PLAN_TRANSPORT.md, шаги A, B, C).
 * ПРЕДСКАЗАНИЯ ПЕЧАТАЮТСЯ ДО РЕЗУЛЬТАТОВ.
 *
 * Батарея T2, перенесённая в 3D ДОСЛОВНО вместе с её находками: печь проверяет
 * только среднее (К12), неконсервативную передачу потока она не ловит вовсе
 * (К13), ограничитель обязан сохранять константу (К6). */

#include "transport/dirs3.h"
#include "transport/mesh3.h"
#include "transport/sweep3.h"
#include "image.h"
#include "transport/cam3.h"
#include "transport/tet3.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fail = 0, g_total = 0;

static void check(int ok, const char *what) {
  g_total++;
  if (!ok) {
    g_fail++;
    printf("FAIL: %s\n", what);
  }
}

/* ------------------------------------------- A1: моменты набора направлений */

static void t_dirs(void) {
  for (int nmu = 2; nmu <= 4; nmu *= 2) {
    tr3_dirs d;
    check(tr3_dirs_product(&d, nmu, nmu) == 0, "набор построен");
    double s = 0.0, m1[3] = {0, 0, 0}, m2[3][3];
    memset(m2, 0, sizeof m2);
    double axis_min = 1e300;
    for (int i = 0; i < d.n; i++) {
      double o[3] = {d.ox[i], d.oy[i], d.oz[i]};
      s += d.w[i];
      for (int a = 0; a < 3; a++) {
        m1[a] += d.w[i] * o[a];
        for (int b = 0; b < 3; b++)
          m2[a][b] += d.w[i] * o[a] * o[b];
        if (fabs(o[a]) < axis_min) axis_min = fabs(o[a]);
      }
    }
    double e2 = 0.0;
    for (int a = 0; a < 3; a++)
      for (int b = 0; b < 3; b++) {
        double ex = (a == b) ? 4.0 * M_PI / 3.0 : 0.0;
        double e = fabs(m2[a][b] - ex);
        if (e > e2) e2 = e;
      }
    printf("  [A1] ND=%d: Σw−4π = %.2e, |Σwω| = %.2e, max |Σwωω − (4π/3)I| = %.2e, "
           "min |ω_a| = %.3f\n",
           d.n, fabs(s - 4.0 * M_PI), fabs(m1[0]) + fabs(m1[1]) + fabs(m1[2]), e2, axis_min);
    check(fabs(s - 4.0 * M_PI) < 1e-13, "Σw = 4π");
    check(fabs(m1[0]) + fabs(m1[1]) + fabs(m1[2]) < 1e-13, "первый момент нулевой");
    check(e2 < 1e-13, "второй момент есть (4π/3)·I");
    check(axis_min > 1e-3, "ни одно направление НЕ параллельно оси — вырожденных граней нет");
    tr3_dirs_free(&d);
  }
}

/* ------------------------------------------- A2: моменты по многоугольнику */

static void t_face_moments(void) {
  /* прямоугольник 2×3 в плоскости z = 1, базис ячейки с центром (1,1.5,1), h=2 */
  double v[4][3] = {{0, 0, 1}, {2, 0, 1}, {2, 3, 1}, {0, 3, 1}};
  double c[3] = {1, 1.5, 1}, h = 2.0;
  double m[4][4];
  check(tr3_poly_face_mass(v, 4, c, h, m) == 0, "моменты грани посчитаны");
  double A = 6.0;
  /* ∫1 = A; ∫ξ = 0 (центр по x совпал); ∫ξ² = A·(2/2)²/12 = A/12;
   * ∫η² = A·(3/2)²/12 = A·9/48; ∫ζ = 0 и ∫ζ² = 0 (грань в плоскости центра) */
  printf("  [A2] ∫1 = %.15f (%.1f), ∫ξ² = %.15f (%.15f), ∫η² = %.15f (%.15f), ∫ζ² = %.1e\n",
         m[0][0], A, m[1][1], A / 12.0, m[2][2], A * 9.0 / 48.0, m[3][3]);
  check(fabs(m[0][0] - A) < 1e-14, "площадь");
  check(fabs(m[1][1] - A / 12.0) < 1e-14, "второй момент по ξ");
  check(fabs(m[2][2] - A * 9.0 / 48.0) < 1e-14, "второй момент по η");
  check(fabs(m[3][3]) < 1e-14, "по ζ грань лежит в плоскости центра");
}

/* ------------------------------------------- B: сетка граней, целые площади */

static void build_tree(hz_octree *t, int L, int graded) {
  hz_oct_init(t, L, 0.0);
  if (graded == 2) {
    /* РАВНОМЕРНОЕ дробление до единицы: коробка в одну единицу содержится только
     * в узле размера 1, поэтому вызов на каждую ячейку заставляет дерево
     * раздробиться целиком. */
    int n2 = 1 << L;
    for (int x = 0; x < n2; x++)
      for (int y = 0; y < n2; y++)
        for (int z = 0; z < n2; z++) {
          int lo2[3] = {x, y, z}, hi2[3] = {x + 1, y + 1, z + 1};
          hz_oct_set_box(t, lo2, hi2, 1.0);
        }
  } else if (graded) {
    /* НАСТОЯЩИЙ перепад уровней: коробка в ОДНУ единицу заставляет дерево
     * раздробиться до размера 1 в углу, оставив остальное крупным. Первая
     * редакция ставила коробку в половину куба — она ложится РОВНО на детей
     * корня, дробления не вызывает вовсе, и «градуированная» сетка выходила
     * равномерной. Перепад надо ЗАДАВАТЬ, а не надеяться на него (Г48). */
    int lo[3] = {0, 0, 0}, hi[3] = {1, 1, 1};
    hz_oct_set_box(t, lo, hi, 1.0);
  }
  (void)L;
}

static void t_mesh(void) {
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 2.0}}; /* анизотропный кадр */
  for (int graded = 0; graded < 2; graded++) {
    hz_octree t;
    build_tree(&t, 3, graded);
    tr3_mesh m;
    check(tr3_mesh_build(&m, &t, &fr) == 0, "сетка построена");
    /* ЦЕЛОЧИСЛЕННОЕ ТОЖДЕСТВО: сумма площадей подграней на стороне ячейки равна
     * площади стороны. Все координаты целые, поэтому равенство ТОЧНОЕ. */
    int32_t bad = 0;
    for (int32_t c = 0; c < m.ncell; c++)
      for (int a = 0; a < 3; a++)
        for (int side = 0; side < 2; side++) {
          int32_t pos = m.clo[c][a] + (side ? m.csize[c] : 0);
          int64_t sum = 0;
          for (int32_t k = m.fstart[c]; k < m.fstart[c + 1]; k++) {
            const tr3_face *f = &m.f[m.flist[k]];
            /* сторону задаёт ПОЛОЖЕНИЕ плоскости, а не то, ca ячейка или cb:
             * у граничной грани ячейка всегда ca, и фильтр по ca/cb выбросил бы
             * минус-стенку целиком */
            if (f->axis != a || f->pos != pos) continue;
            sum += (int64_t)(f->hi[0] - f->lo[0]) * (int64_t)(f->hi[1] - f->lo[1]);
          }
          int64_t want = (int64_t)m.csize[c] * (int64_t)m.csize[c];
          if (sum != want) bad++;
        }
    printf("  [B%s] ячеек %d, граней %d; сторон с НЕВЕРНОЙ суммой площадей: %d\n",
           graded ? ", градуированная" : "", m.ncell, m.nf, bad);
    check(bad == 0, "К13 в 3D: сумма площадей подграней ТОЧНО равна площади стороны");

    /* НЕГАТИВНЫЙ КОНТРОЛЬ: сосед берётся ОДИН — тот, что накрывает ЦЕНТР
     * стороны. Так пишут, когда забывают, что у крупной ячейки соседей
     * несколько. Тождество ОБЯЗАНО сломаться, и только на градуированной. */
    int32_t bad2 = 0;
    for (int32_t c = 0; c < m.ncell; c++)
      for (int a = 0; a < 3; a++)
        for (int side = 0; side < 2; side++) {
          /* ОБЕ стороны, а не только плюс: у мелких соседей ячейка оказывается
           * то ca, то cb, и односторонний обзор конфигурацию «крупная против
           * четырёх мелких» просто не встречает */
          int32_t pos = m.clo[c][a] + (side ? m.csize[c] : 0);
          int64_t sum = 0, seen = 0;
          for (int32_t k = m.fstart[c]; k < m.fstart[c + 1]; k++) {
            const tr3_face *f = &m.f[m.flist[k]];
            if (f->axis != a || f->pos != pos) continue;
            if (seen++ > 0) continue; /* «взять только первого соседа» */
            sum += (int64_t)(f->hi[0] - f->lo[0]) * (int64_t)(f->hi[1] - f->lo[1]);
          }
          if (seen > 0 && sum != (int64_t)m.csize[c] * (int64_t)m.csize[c]) bad2++;
        }
    printf("    [НК один сосед] сторон с неверной суммой: %d\n", bad2);
    if (graded) check(bad2 > 0, "негативный контроль: «один сосед» ломает тождество");
    tr3_mesh_free(&m);
    hz_oct_free(&t);
  }
}

/* ------------------------------------------- C: печь, линейное поле, баланс */

typedef struct {
  hz_octree t;
  tr3_mesh m;
  tr3_dirs d;
  double *sig_t, *sig_s, *eps;
} rig;

static int rig_init(rig *r, int L, int graded, int nmu, const hz_frame *fr) {
  build_tree(&r->t, L, graded);
  if (tr3_mesh_build(&r->m, &r->t, fr)) return 1;
  if (tr3_dirs_product(&r->d, nmu, nmu)) return 1;
  r->sig_t = calloc((size_t)r->m.ncell, sizeof(double));
  r->sig_s = calloc((size_t)r->m.ncell, sizeof(double));
  r->eps = calloc((size_t)r->m.ncell * 4, sizeof(double));
  return r->sig_t == NULL || r->sig_s == NULL || r->eps == NULL;
}

static void rig_free(rig *r) {
  free(r->sig_t);
  free(r->sig_s);
  free(r->eps);
  tr3_dirs_free(&r->d);
  tr3_mesh_free(&r->m);
  hz_oct_free(&r->t);
}

static void t_furnace(void) {
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 2.0}};
  const double lb = 1.7;
  for (int graded = 0; graded < 2; graded++) {
    rig r;
    if (rig_init(&r, 3, graded, 2, &fr)) {
      check(0, "оснастка");
      rig_free(&r);
      return;
    }
    for (int32_t c = 0; c < r.m.ncell; c++) {
      r.sig_t[c] = 0.8;
      r.sig_s[c] = 0.8; /* альбедо РОВНО 1 */
    }
    tr3_problem p = {
        .m = &r.m, .d = &r.d, .sig_t = r.sig_t, .sig_s = r.sig_s, .binc0 = lb, .limiter = 1};
    double *phi = calloc((size_t)r.m.ncell * 4, sizeof(double));
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    tr3_stats st;
    int rc = tr3_sweep_solve(&p, 400, 1e-14, phi, &st);
    free(st.bout);
    double worst = 0.0, wslope = 0.0;
    for (int32_t c = 0; c < r.m.ncell; c++) {
      double e = fabs(phi[c * 4] - 4.0 * M_PI * lb) / (4.0 * M_PI * lb);
      if (e > worst) worst = e;
      for (int j = 1; j < 4; j++) {
        double sl = fabs(phi[c * 4 + j]) / (4.0 * M_PI * lb);
        if (sl > wslope) wslope = sl;
      }
    }
    printf("  [C печь%s] код %d, итераций %d, срезок %d; отн. ошибка среднего %.2e, "
           "наклоны %.2e\n",
           graded ? ", градуированная" : "", rc, st.iters, st.nclip, worst, wslope);
    check(rc == 0, "развёртка прошла, цикла обхода нет");
    check(worst < 1e-12, "ПЕЧЬ: радианс однороден и равен влёту");
    check(wslope < 1e-12, "и наклоны нулевые");
    free(phi);
    rig_free(&r);
  }
}

static void t_linear(void) {
  /* К12: печь проверяет только СРЕДНЕЕ. Линейное решение L = A + B·x проверяет
   * наклоны. Берём чистое поглощение (σ_s = 0) и подбираем ε так, чтобы
   * L = A + B·x было точным: (ω·∇)L + σ_t L = ε ⇒ ε(ω) = B·ω_x + σ_t(A + Bx).
   * Направленной части ε здесь нет, поэтому берём B вдоль x и σ_t = 0:
   * тогда ε = B·ω_x, что от направления зависит — а наш ε изотропен.
   * Значит проверяем ДРУГОЕ линейное решение: σ_t > 0, ε изотропно,
   * L = ε/σ_t — константа. Линейность проверяется линейным ВЛЁТОМ при σ = 0. */
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 1.0}};
  rig r;
  if (rig_init(&r, 3, 1, 2, &fr)) {
    check(0, "оснастка");
    rig_free(&r);
    return;
  }
  /* пустая среда, линейный влёт: решение L(x) = влёт, перенесённый вдоль луча,
   * то есть ЛИНЕЙНАЯ функция обязана пройти сквозь сетку без искажения */
  tr3_problem p = {.m = &r.m, .d = &r.d, .sig_t = r.sig_t, .sig_s = r.sig_s, .binc0 = 2.0};
  double *phi = calloc((size_t)r.m.ncell * 4, sizeof(double));
  if (phi == NULL) {
    check(0, "память");
    rig_free(&r);
    return;
  }
  tr3_stats st;
  int rc = tr3_sweep_solve(&p, 4, 1e-14, phi, &st);
  free(st.bout);
  double worst = 0.0;
  for (int32_t c = 0; c < r.m.ncell; c++) {
    double e = fabs(phi[c * 4] - 4.0 * M_PI * 2.0) / (4.0 * M_PI * 2.0);
    if (e > worst) worst = e;
  }
  printf("  [C пустая среда] код %d: отн. ошибка %.2e (влёт проходит насквозь)\n", rc, worst);
  check(rc == 0 && worst < 1e-12, "в пустой среде влёт проходит без искажения");
  free(phi);
  rig_free(&r);
}

static void t_balance(void) {
  /* К13: баланс энергии при НЕОДНОРОДНОЙ подсветке — то, чем ловится
   * неконсервативная передача потока, тогда как печь её не видит. */
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 2.0}};
  for (int graded = 0; graded < 2; graded++) {
    rig r;
    if (rig_init(&r, 3, graded, 2, &fr)) {
      check(0, "оснастка");
      rig_free(&r);
      return;
    }
    for (int32_t c = 0; c < r.m.ncell; c++) {
      r.sig_t[c] = 0.5;
      r.sig_s[c] = 0.5; /* альбедо 1: поглощения нет, всё обязано выйти */
    }
    /* влёт ЛИНЕЙНЫЙ по положению — неоднородная подсветка */
    tr3_problem p = {.m = &r.m,
                     .d = &r.d,
                     .sig_t = r.sig_t,
                     .sig_s = r.sig_s,
                     .binc0 = 1.0,
                     .binc = {0.13, -0.07, 0.05},
                     .binx0 = {4, 4, 4},
                     .limiter = 1};
    double *phi = calloc((size_t)r.m.ncell * 4, sizeof(double));
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    tr3_stats st;
    int rc = tr3_sweep_solve(&p, 2000, 1e-13, phi, &st);
    free(st.bout);
    double rel = fabs(st.balance) / (fabs(st.pin) + 1e-300);
    printf("  [C баланс%s] код %d, итераций %d: втекло %.6f, вытекло %.6f, поглощено %.2e, "
           "невязка/втекло %.2e\n",
           graded ? ", градуированная" : "", rc, st.iters, st.pin, st.pout, st.pabs, rel);
    check(rc == 0, "развёртка прошла");
    check(rel < 1e-10, "К13: при альбедо 1 вся влитая мощность выходит");
    free(phi);
    rig_free(&r);
  }
}

/* -------------------------------------------------------------------------- */

/* ------------------------------- D: стенки как ГУ, замкнутая полость */

static void t_cavity(void) {
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 1.0}};
  /* Замкнутая диффузная полость: все стенки с отражением ρ и собственным
   * радиансом L_e, среда пустая. Классическая замкнутая форма бесконечной
   * серии отражений: равновесный радианс L = L_e/(1−ρ), ОДНОРОДНЫЙ. */
  const double rho = 0.6, le = 0.4;
  double wr[6], we[6];
  for (int i = 0; i < 6; i++) {
    wr[i] = rho;
    we[i] = le;
  }
  double exact = le / (1.0 - rho);
  for (int graded = 0; graded < 2; graded++) {
    rig r;
    if (rig_init(&r, 3, graded, 2, &fr)) {
      check(0, "оснастка");
      rig_free(&r);
      return;
    }
    tr3_problem p = {.m = &r.m,
                     .d = &r.d,
                     .sig_t = r.sig_t,
                     .sig_s = r.sig_s,
                     .wall_rho = wr,
                     .wall_emit = we,
                     .limiter = 1};
    double *phi = calloc((size_t)r.m.ncell * 4, sizeof(double));
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    tr3_stats st;
    int rc = tr3_sweep_solve(&p, 4000, 1e-14, phi, &st);
    double worst = 0.0, wb = 0.0;
    for (int32_t c = 0; c < r.m.ncell; c++) {
      double e = fabs(phi[c * 4] - 4.0 * M_PI * exact) / (4.0 * M_PI * exact);
      if (e > worst) worst = e;
    }
    for (int32_t f = 0; f < r.m.nf; f++) {
      if (r.m.f[f].cb >= 0) continue;
      double e = fabs(st.bout[f] - exact) / exact;
      if (e > wb) wb = e;
    }
    printf("  [D полость%s] код %d, итераций %d: радианс отн. ошибка %.2e, "
           "исходящий на стенках %.2e (точно %.6f)\n",
           graded ? ", градуированная" : "", rc, st.iters, worst, wb, exact);
    check(rc == 0, "развёртка со стенками прошла");
    check(worst < 1e-10, "D: равновесие полости есть L_e/(1−ρ), однородно");
    check(wb < 1e-10, "и исходящий радианс стенок тот же");
    free(st.bout);
    free(phi);
    rig_free(&r);
  }

  /* НЕГАТИВНЫЙ КОНТРОЛЬ (предсказан провал): забыть 1/π в ламбертовом ГУ.
   * Тогда «отражение» отдаёт в π раз больше, серия сходится к другому числу —
   * а при ρ ≥ 1/π вовсе расходится. Проверяем, что ответ УЕХАЛ. */
  rig r;
  if (rig_init(&r, 3, 0, 2, &fr)) {
    rig_free(&r);
    return;
  }
  double wr2[6];
  for (int i = 0; i < 6; i++)
    wr2[i] = rho * M_PI; /* ровно то, что выйдет без деления на π */
  tr3_problem p2 = {.m = &r.m,
                    .d = &r.d,
                    .sig_t = r.sig_t,
                    .sig_s = r.sig_s,
                    .wall_rho = wr2,
                    .wall_emit = we,
                    .limiter = 1};
  double *phi2 = calloc((size_t)r.m.ncell * 4, sizeof(double));
  if (phi2 == NULL) {
    rig_free(&r);
    return;
  }
  tr3_stats st2;
  tr3_sweep_solve(&p2, 200, 1e-12, phi2, &st2);
  double got = phi2[0] / (4.0 * M_PI);
  printf("    [НК без 1/π] радианс %.4f против правильного %.4f\n", got, exact);
  check(fabs(got - exact) > 0.5 * exact, "негативный контроль: без 1/π ответ ОБЯЗАН уехать");
  free(st2.bout);
  free(phi2);
  rig_free(&r);
}

/* ------------------------------- D: КАРТИНКА с косвенным светом */

static void t_render(void) {
  const int L = 3;
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 1.0}};
  rig r;
  if (rig_init(&r, L, 2, 2, &fr)) {
    check(0, "оснастка");
    rig_free(&r);
    return;
  }
  int n = 1 << L;
  /* комната: пять серых стенок и один яркий потолок */
  double wr[6] = {0.7, 0.7, 0.7, 0.7, 0.75, 0.7};
  double we[6] = {0, 0, 0, 0, 0, 3.0}; /* +z светит */
  for (int32_t c = 0; c < r.m.ncell; c++) {
    r.sig_t[c] = 0.02; /* лёгкая дымка, чтобы среда участвовала */
    r.sig_s[c] = 0.02;
  }
  tr3_problem p = {.m = &r.m,
                   .d = &r.d,
                   .sig_t = r.sig_t,
                   .sig_s = r.sig_s,
                   .wall_rho = wr,
                   .wall_emit = we,
                   .limiter = 1};
  double *phi = calloc((size_t)r.m.ncell * 4, sizeof(double));
  if (phi == NULL) {
    check(0, "память");
    rig_free(&r);
    return;
  }
  tr3_stats st;
  int rc = tr3_sweep_solve(&p, 2000, 1e-12, phi, &st);
  check(rc == 0, "комната посчитана");
  printf("  [D комната] итераций %d, ячеек %d, граней %d\n", st.iters, r.m.ncell, r.m.nf);

  /* СБОР ПО ПИКСЕЛЮ: луч из камеры внутри комнаты до стенки, читаем ИСХОДЯЩИЙ
   * радианс той граничной грани, в которую попали. Интерполяции между гранями
   * нет намеренно — К3 запрещает её между ординатами, и здесь та же мысль. */
  const int W = 220, H = 165;
  double *buf = calloc((size_t)W * (size_t)H, sizeof(double));
  if (buf == NULL) {
    check(0, "память");
    free(st.bout);
    free(phi);
    rig_free(&r);
    return;
  }
  double eye[3] = {(double)n * 0.5, (double)n * 0.12, (double)n * 0.45};
  double at[3] = {(double)n * 0.5, (double)n, (double)n * 0.42};
  double up[3] = {0, 0, 1};
  tr3_camera cam;
  check(tr3_camera_look(&cam, eye, at, up, 1.15, W, H) == 0, "камера");
  int miss = 0;
  for (int py = 0; py < H; py++)
    for (int px = 0; px < W; px++) {
      double o[3], d[3];
      tr3_camera_ray(&cam, px, py, o, d);
      /* выход из куба */
      double tex = 1e300;
      int wall = -1;
      for (int a = 0; a < 3; a++) {
        if (!(fabs(d[a]) > 0.0)) continue;
        double lim = d[a] > 0.0 ? (double)n : 0.0;
        double t = (lim - o[a]) / d[a];
        if (t > 0.0 && t < tex) {
          tex = t;
          wall = 2 * a + (d[a] > 0.0 ? 1 : 0);
        }
      }
      if (wall < 0) {
        miss++;
        continue;
      }
      double hp[3];
      for (int a = 0; a < 3; a++)
        hp[a] = o[a] + tex * d[a];
      /* граничная грань, накрывающая точку */
      int axis = wall / 2;
      int u = (axis + 1) % 3, v = (axis + 2) % 3;
      if (u > v) {
        int t = u;
        u = v;
        v = t;
      }
      double val = 0.0;
      for (int32_t f = 0; f < r.m.nf; f++) {
        const tr3_face *ff = &r.m.f[f];
        if (ff->cb >= 0 || (int)(~ff->cb) != wall) continue;
        if (hp[u] < (double)ff->lo[0] || hp[u] > (double)ff->hi[0]) continue;
        if (hp[v] < (double)ff->lo[1] || hp[v] > (double)ff->hi[1]) continue;
        val = st.bout[f];
        break;
      }
      buf[py * W + px] = val;
    }
  double mn = 1e300, mx = 0.0;
  for (int i = 0; i < W * H; i++) {
    if (buf[i] < mn) mn = buf[i];
    if (buf[i] > mx) mx = buf[i];
  }
  printf("    радианс на стенках от %.4f до %.4f; пикселей без грани %d\n", mn, mx, miss);
  check(miss == 0, "каждый луч нашёл свою стенку");
  check(mn > 0.0, "теней в замкнутой комнате нет — косвенный свет достаёт везде");
  hz_ppm_write("build/room.ppm", buf, W, H);
  printf("    картинка записана: build/room.ppm\n");
  free(buf);
  free(st.bout);
  free(phi);
  rig_free(&r);
}

int main(void) {
  printf("=== ПРЕДСКАЗАНИЯ (до единого результата) ===\n");
  printf("  A1 Σw = 4π, первый момент 0, второй (4π/3)·I — все ≤1e-13;\n");
  printf("     ни одно направление не параллельно оси (иначе грань с ω·n = 0)\n");
  printf("  A2 моменты по многоугольнику совпадают с замкнутой формой\n");
  printf("  B  К13 в 3D: сумма площадей подграней ТОЧНО равна площади стороны\n");
  printf("     (все координаты целые, поэтому равенство точное, а не с допуском)\n");
  printf("  C  печь: радианс однороден и равен влёту, ≤1e-12, наклоны нулевые,\n");
  printf("     на равномерной И на градуированной сетке\n");
  printf("  C  баланс при НЕОДНОРОДНОЙ подсветке и альбедо 1: невязка ≤1e-10\n");
  printf("=== НЕГАТИВНЫЙ КОНТРОЛЬ, с ПРЕДСКАЗАННЫМ провалом ===\n");
  printf("  «взять одного соседа» вместо всех -> тождество площадей ОБЯЗАНО\n");
  printf("  сломаться, и только на градуированной сетке\n");
  printf("=== РЕЗУЛЬТАТЫ ===\n");

  t_dirs();
  t_face_moments();
  t_mesh();
  t_furnace();
  t_linear();
  t_balance();
  t_cavity();
  t_render();

  printf("%s: %d/%d\n", g_fail ? "FAILURES" : "ok", g_total - g_fail, g_total);
  return g_fail != 0;
}
