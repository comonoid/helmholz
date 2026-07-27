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
    free(st.sout);
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
    check(st.nclip == 0,
          "К6: на СОШЕДШЕМСЯ решении ограничитель не активен — константу он не трогает");
    check(wslope < 1e-12, "и наклоны нулевые");
    free(phi);
    rig_free(&r);
  }
}

static void t_linear_field(void) {
  /* К12 В ТРЁХМЕРИИ: печь проверяет только СРЕДНЕЕ, а наклоны — две трети
   * неизвестных DG1. Точное решение L = A + B·x достигается при σ_s = 0 и
   * ε(x,ω) = B·ω_x + σ_t·(A + B·x), причём первый член НАПРАВЛЕННЫЙ. */
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 1.0}};
  const double A = 1.3, B = 0.21;
  for (int graded = 0; graded < 2; graded++) {
    rig r;
    if (rig_init(&r, 3, graded, 2, &fr)) {
      check(0, "оснастка");
      rig_free(&r);
      return;
    }
    for (int32_t c = 0; c < r.m.ncell; c++) {
      r.sig_t[c] = 0.7;
      r.sig_s[c] = 0.0;
      double s = (double)r.m.csize[c];
      double xc = (double)r.m.clo[c][0] + 0.5 * s;
      r.eps[c * 4 + 0] = 0.7 * (A + B * xc); /* среднее по ячейке */
      r.eps[c * 4 + 1] = 0.7 * B * s;        /* наклон: b_1 = (x − x_c)/s */
    }
    tr3_problem p = {.m = &r.m,
                     .d = &r.d,
                     .sig_t = r.sig_t,
                     .sig_s = r.sig_s,
                     .eps = r.eps,
                     .eps_dir = {B, 0.0, 0.0},
                     .binc0 = A,
                     .binc = {B, 0, 0},
                     .binx0 = {0, 0, 0},
                     .limiter = 0};
    double *phi = calloc((size_t)r.m.ncell * 4, sizeof(double));
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    tr3_stats st;
    int rc = tr3_sweep_solve(&p, 4, 1e-14, phi, &st);
    free(st.bout);
    free(st.sout);
    double wm = 0.0, ws = 0.0;
    for (int32_t c = 0; c < r.m.ncell; c++) {
      double s = (double)r.m.csize[c];
      double xc = (double)r.m.clo[c][0] + 0.5 * s;
      double em = fabs(phi[c * 4] - 4.0 * M_PI * (A + B * xc));
      double es = fabs(phi[c * 4 + 1] - 4.0 * M_PI * B * s);
      if (em > wm) wm = em;
      if (es > ws) ws = es;
    }
    printf("  [К12 линейное поле%s] код %d: ошибка среднего %.2e, ошибка НАКЛОНА %.2e\n",
           graded ? ", градуированная" : "", rc, wm, ws);
    check(rc == 0 && wm < 1e-12, "К12: линейное решение воспроизведено по среднему");
    check(ws < 1e-12, "и ПО НАКЛОНУ — то, чего печь не проверяет вовсе");
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
  free(st.sout);
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
    free(st.sout);
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
    /* ХРАНИМОЕ ЕСТЬ КОЭФФИЦИЕНТЫ DG1, А НЕ ЗНАЧЕНИЕ (К39). Читать `bout[f·4]`
     * как радианс НЕЛЬЗЯ: базис ячейки центрирован на ЯЧЕЙКЕ, а грань лежит на
     * её краю, поэтому у постоянного поля c представитель минимальной нормы
     * равен (0.8c, 0.4c, 0, 0) и значение собирается только вместе с базисом.
     * Прежде это сходилось случайно — откат к константе обнулял наклоны, то
     * есть тест проходил ровно потому, что проекция была сломана. */
    for (int32_t f = 0; f < r.m.nf; f++) {
      if (r.m.f[f].cb >= 0) continue;
      double v[4][3];
      tr3_face_corners(&r.m, f, v);
      int32_t c = r.m.f[f].ca;
      double s = (double)r.m.csize[c];
      for (int q = 0; q < 4; q++) {
        double val = st.bout[f * 4];
        for (int a = 0; a < 3; a++) {
          double xu = (v[q][a] - r.m.fr.o[a]) / r.m.fr.u[a];
          val += st.bout[f * 4 + a + 1] * ((xu - ((double)r.m.clo[c][a] + 0.5 * s)) / s);
        }
        double e = fabs(val - exact) / exact;
        if (e > wb) wb = e;
      }
    }
    printf("  [D полость%s] код %d, итераций %d: радианс отн. ошибка %.2e, "
           "исходящий на стенках %.2e (точно %.6f)\n",
           graded ? ", градуированная" : "", rc, st.iters, worst, wb, exact);
    check(rc == 0, "развёртка со стенками прошла");
    check(worst < 1e-10, "D: равновесие полости есть L_e/(1−ρ), однородно");
    check(wb < 1e-10, "и исходящий радианс стенок тот же");
    free(st.bout);
    free(st.sout);
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
  free(st2.sout);
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
    free(st.sout);
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
        /* DG1 ПО ПОЛОЖЕНИЮ: радианс на грани линейный, и читать его надо В ТОЧКЕ
         * попадания, а не константой на грань. Константа и давала квадраты. */
        {
          int32_t cc = ff->ca;
          double sz = (double)r.m.csize[cc];
          double bb[4] = {1.0, 0, 0, 0};
          for (int a = 0; a < 3; a++) {
            double xu = (hp[a] - r.m.fr.o[a]) / r.m.fr.u[a];
            bb[a + 1] = (xu - ((double)r.m.clo[cc][a] + 0.5 * sz)) / sz;
          }
          val = 0.0;
          for (int i = 0; i < 4; i++)
            val += st.bout[f * 4 + i] * bb[i];
        }
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
  hz_ppm_write("img/room.ppm", buf, W, H);
  printf("    картинка записана: img/room.ppm\n");
  free(buf);
  free(st.bout);
  free(st.sout);
  free(phi);
  rig_free(&r);
}

/* ------------------------- D2: разрезанные ячейки, флюидная геометрия */

static void t_cut(void) {
  const hz_frame fr = {{0, 0, 0}, {1, 1, 1}};
  const int L = 4, N = 1 << L;
  hz_octree t;
  hz_oct_init(&t, L, 0.0);
  for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
      for (int z = 0; z < N; z++) {
        int lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        hz_oct_set_box(&t, lo, hi, 1.0);
      }
  hz_surftab stab;
  hz_facettab ftab;
  hz_cutmap cmap;
  hz_surftab_init(&stab);
  hz_facettab_init(&ftab);
  hz_cutmap_init(&cmap);
  double sc[3] = {8.3, 8.7, 7.6}, sr = 3.4;
  hz_surf sp = {HZ_SURF_SPHERE, {sc[0], sc[1], sc[2], sr, 0, 0, 0}, 1, 0};
  int32_t si = hz_surftab_add(&stab, &sp);
  int32_t f0 = 0;
  int32_t nfac = hz_surf_facet_sphere(&ftab, &fr, sc, sr, 1, HZ_FIT_MEAN_SAGITTA, si, &f0);
  check(nfac > 0, "сфера фасетизирована");
  /* Г45 В ДЕЙСТВИИ: ключ боковой таблицы обязан СТРОГО ВОЗРАСТАТЬ, а индексы
   * узлов идут порядком ВЫДЕЛЕНИЯ, а не порядком (x,y,z). Первая редакция этого
   * теста добавляла записи как попало и не смотрела на код возврата — половина
   * молча отвергалась, тело выходило дырявым, и объём материала был 105 вместо
   * 165. Ровно тот сценарий, про который Г45 и написана. */
  typedef struct {
    int32_t cell, f[HZ_P3_MAXH], nf;
  } rec_t;
  static rec_t rec[4096];
  int nrec = 0;
  for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
      for (int z = 0; z < N; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        int ns = hz_facets_for_box(&ftab, f0, nfac, lo, hi, sel, HZ_P3_MAXH);
        if (ns <= 0) continue;
        rec[nrec].cell = hz_oct_leaf(&t, x, y, z);
        rec[nrec].nf = ns;
        for (int j = 0; j < ns; j++)
          rec[nrec].f[j] = sel[j];
        nrec++;
      }
  for (int i = 1; i < nrec; i++) {
    rec_t tmp = rec[i];
    int j = i - 1;
    while (j >= 0 && rec[j].cell > tmp.cell) {
      rec[j + 1] = rec[j];
      j--;
    }
    rec[j + 1] = tmp;
  }
  int nadd = 0;
  for (int i = 0; i < nrec; i++)
    if (hz_cutmap_add(&cmap, rec[i].cell, rec[i].f, rec[i].nf) == 0) nadd++;
  check(nadd == nrec, "Г45: ВСЕ записи легли в боковую таблицу (ключи по возрастанию)");
  tr3_mesh m;
  check(tr3_mesh_build(&m, &t, &fr) == 0, "сетка");
  uint8_t *solid = calloc((size_t)m.ncell, 1);
  for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
      for (int z = 0; z < N; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        int ns = hz_facets_for_box(&ftab, f0, nfac, lo, hi, sel, HZ_P3_MAXH);
        if (ns == 0) solid[m.cellof[hz_oct_leaf(&t, x, y, z)]] = 1; /* Г38: 0 = ПОЛНАЯ */
      }
  tr3_cut cu;
  check(tr3_cut_build(&cu, &m, &ftab, &cmap, solid) == 0, "разрез");
  check(cu.nbad == 0, "условие 1:1 у разрезанных ячеек соблюдено");

  /* 1. ОБЪЁМ МАТЕРИАЛА: коробки минус флюид = объём шара, аналитически */
  double vmat = 0.0;
  for (int32_t c = 0; c < m.ncell; c++) {
    double s = (double)m.csize[c];
    vmat += s * s * s - cu.mvol[c][0][0];
  }
  free(solid);
  double vex = 4.0 / 3.0 * M_PI * sr * sr * sr;
  /* РЕШАЮЩАЯ СВЕРКА: материал считается НЕЗАВИСИМО, проверенным путём
   * hz_poly3_cut (пункт 4), и сумма «флюид + материал = коробка» обязана
   * держаться поячеечно. Если нет — виновата флюидная бухгалтерия, а не
   * фасетизация. */
  double vmat2 = 0.0, wcell = 0.0;
  for (int32_t c = 0; c < m.ncell; c++) {
    const hz_cutrec *rr = hz_cutmap_find(&cmap, m.node[c]);
    double s = (double)m.csize[c];
    if (rr == NULL) {
      vmat2 += cu.solid[c] ? s * s * s : 0.0;
      continue;
    }
    hz_hspace hh[HZ_P3_MAXH];
    int32_t hd[HZ_P3_MAXH];
    int nh = hz_cutmap_hspaces(&ftab, &cmap, rr, hh, hd, NULL, HZ_P3_MAXH);
    int32_t lo2[3] = {m.clo[c][0], m.clo[c][1], m.clo[c][2]};
    int32_t hi2[3] = {lo2[0] + m.csize[c], lo2[1] + m.csize[c], lo2[2] + m.csize[c]};
    hz_poly3 pp;
    double vm = 0.0;
    if (hz_poly3_cut(&pp, lo2, hi2, hh, hd, nh) == HZ_P3_OK) vm = hz_poly3_volume(&pp, &fr);
    vmat2 += vm;
    double e = fabs(vm + cu.mvol[c][0][0] - s * s * s);
    if (e > wcell) wcell = e;
  }
  printf("      материал независимым путём %.5f; max |флюид+материал−коробка| = %.3e\n", vmat2,
         wcell);
  check(wcell < 1e-12, "поячеечно: флюид и материал дополняют коробку ТОЧНО");
  printf("  [D2] объём материала %.5f против шара %.5f (отн. %.3e), элементов %d\n", vmat, vex,
         fabs(vmat - vex) / vex, cu.nse);
  check(fabs(vmat - vex) / vex < 0.02, "флюид дополняет материал до шара (фасетизация k=1)");

  /* 2. ЗАМКНУТОСТЬ ФЛЮИДНОЙ ОБЛАСТИ: Σ(внешняя нормаль × площадь) = 0.
   * Это и есть проверка того, что ни одна грань не потеряна и ни одна не
   * посчитана дважды — теорема о дивергенции для постоянного поля. */
  double worst = 0.0;
  int nchecked = 0;
  for (int32_t c = 0; c < m.ncell; c++) {
    if (cu.solid[c]) continue;
    if (cu.sestart[c + 1] == cu.sestart[c]) continue; /* только разрезанные */
    double acc[3] = {0, 0, 0};
    for (int32_t k = m.fstart[c]; k < m.fstart[c + 1]; k++) {
      int32_t f = m.flist[k];
      int mine_is_a = m.f[f].ca == c;
      double sgn;
      if (m.f[f].cb < 0) {
        int wall = (int)(~m.f[f].cb);
        sgn = (wall & 1) ? 1.0 : -1.0;
      } else {
        sgn = mine_is_a ? 1.0 : -1.0;
      }
      double a = mine_is_a ? cu.ffm[f][0][0] : cu.ffmb[f][0][0];
      acc[m.f[f].axis] += sgn * a;
    }
    for (int32_t k = cu.sestart[c]; k < cu.sestart[c + 1]; k++) {
      const tr3_selem *se = &cu.se[cu.selist[k]];
      for (int a = 0; a < 3; a++)
        acc[a] += -se->n[a] * se->area; /* наружу ФЛЮИДА = −n */
    }
    double e = fabs(acc[0]) + fabs(acc[1]) + fabs(acc[2]);
    if (e > worst) worst = e;
    nchecked++;
  }
  printf("      замкнутость флюида: %d разрезанных ячеек, max |Σn·A| = %.3e\n", nchecked, worst);
  check(nchecked > 50, "разрезанных ячеек достаточно");
  check(worst < 1e-12, "D2: флюидная область ЗАМКНУТА — ни одна грань не потеряна");

  /* ---------------------------------------------------------------- К39 ---
   * ВЫРОЖДЕНИЕ МАТРИЦЫ МАСС ПЛОСКОГО ЭЛЕМЕНТА. Печь и полость этого не ловят
   * ВООБЩЕ: у однородного равновесия все наклоны нули, и откат к константе там
   * безвреден. Это ровно урок К12, повторённый на поверхностях, поэтому
   * проверка ставится ПРЯМАЯ. */
  {
    double wnul = 0.0, wmom = 0.0, wgrow = 0.0, wgrow_old = 0.0;
    int nneg = 0, nneg_old = 0, ncheck2 = 0;
    /* пробная линейная функция: положительна на всей ячейке (при |b| ≤ ½
     * минимум равен 1 − ½(0.1+0.2+0.15) = 0.775 > 0), поэтому у ИСПРАВНОЙ
     * проекции повода отступать к константе нет ни у одного элемента */
    const double ex[4] = {1.0, 0.1, -0.2, 0.15};
    for (int32_t e = 0; e < cu.nse + m.nf; e++) {
      const double (*M)[4];
      double nul[4];
      if (e < cu.nse) {
        M = cu.se[e].m;
        memcpy(nul, cu.se[e].nul, sizeof nul);
      } else {
        int32_t f = e - cu.nse;
        if (m.f[f].cb >= 0) continue;
        M = cu.ffm[f];
        tr3_face_null(&m, f, m.f[f].ca, nul);
      }
      if (!(M[0][0] > 0.0)) continue;
      ncheck2++;
      /* 1. НУЛЕВОЙ ВЕКТОР ВЕРЕН И ВЫРОЖДЕНИЕ ТОЧНОЕ */
      double nm = 0.0, nv = 0.0, r = 0.0;
      for (int i = 0; i < 4; i++) {
        nv += nul[i] * nul[i];
        double acc = 0.0;
        for (int j = 0; j < 4; j++) {
          acc += M[i][j] * nul[j];
          nm += M[i][j] * M[i][j];
        }
        r += acc * acc;
      }
      double rel = sqrt(r) / (sqrt(nm) * sqrt(nv));
      if (rel > wnul) wnul = rel;
      /* 2. ПРОЕКЦИЯ ВОСПРОИЗВОДИТ ЛИНЕЙНОЕ ПОЛЕ: моменты обязаны совпасть */
      double rr2[4] = {0, 0, 0, 0}, ee[4], nr = 0.0;
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
          rr2[i] += M[i][j] * ex[j];
        nr += rr2[i] * rr2[i];
      }
      nr = sqrt(nr);
      if (tr3_project_plane(M, rr2, nul, ee) != 0) continue;
      double dm = 0.0;
      for (int i = 0; i < 4; i++) {
        double acc = -rr2[i];
        for (int j = 0; j < 4; j++)
          acc += M[i][j] * ee[j];
        dm += acc * acc;
      }
      if (nr > 0.0 && sqrt(dm) / nr > wmom) wmom = sqrt(dm) / nr;
      /* 3. ПРЕДСТАВИТЕЛЬ ОГРАНИЧЕН, И ИМЕННО ЭТО ЛЕЧИТ ОТКАТ. Добавка вдоль
       * нулевого вектора на плоскости НЕ ВИДНА, но вне её растёт без предела, и
       * проверка положительности по углам ячейки съезжает в минус. */
      double ne = 0.0, nex = 0.0;
      for (int i = 0; i < 4; i++) {
        ne += ee[i] * ee[i];
        nex += ex[i] * ex[i];
      }
      double gr = sqrt(ne / nex);
      if (gr > wgrow) wgrow = gr;
      double mn = 1e300;
      for (int k = 0; k < 8; k++) {
        double v = ee[0];
        for (int a = 0; a < 3; a++)
          v += ee[a + 1] * (((k >> a) & 1) ? 0.5 : -0.5);
        if (v < mn) mn = v;
      }
      if (mn < 0.0) nneg++;
      /* НЕГАТИВНЫЙ КОНТРОЛЬ: ПРЕЖНИЙ КОД, воспроизведённый здесь дословно.
       * `solve4` вырождения не замечает — четвёртый ведущий элемент выходит не
       * нулём, а округлением, — и возвращает произвольную добавку вдоль nul. */
      double a4[4][4], b4[4], x4[4];
      memcpy(a4, M, sizeof a4);
      memcpy(b4, rr2, sizeof b4);
      int sing = 0;
      for (int k = 0; k < 4 && !sing; k++) {
        int best = k;
        for (int i = k + 1; i < 4; i++)
          if (fabs(a4[i][k]) > fabs(a4[best][k])) best = i;
        if (!(fabs(a4[best][k]) > 0.0)) {
          sing = 1;
          break;
        }
        if (best != k) {
          for (int j = 0; j < 4; j++) {
            double tt = a4[k][j];
            a4[k][j] = a4[best][j];
            a4[best][j] = tt;
          }
          double tt = b4[k];
          b4[k] = b4[best];
          b4[best] = tt;
        }
        for (int i = k + 1; i < 4; i++) {
          double ff = a4[i][k] / a4[k][k];
          for (int j = k; j < 4; j++)
            a4[i][j] -= ff * a4[k][j];
          b4[i] -= ff * b4[k];
        }
      }
      if (sing) continue;
      for (int i = 3; i >= 0; i--) {
        double s2 = b4[i];
        for (int j = i + 1; j < 4; j++)
          s2 -= a4[i][j] * x4[j];
        x4[i] = s2 / a4[i][i];
      }
      double ne2 = 0.0, mn2 = 1e300;
      for (int i = 0; i < 4; i++)
        ne2 += x4[i] * x4[i];
      if (sqrt(ne2 / nex) > wgrow_old) wgrow_old = sqrt(ne2 / nex);
      for (int k = 0; k < 8; k++) {
        double v = x4[0];
        for (int a = 0; a < 3; a++)
          v += x4[a + 1] * (((k >> a) & 1) ? 0.5 : -0.5);
        if (v < mn2) mn2 = v;
      }
      if (mn2 < 0.0) nneg_old++;
    }
    printf("  [К39] элементов %d: max |M·nul|/(|M||nul|) = %.3e\n", ncheck2, wnul);
    printf("        проекция линейного поля: max невязка моментов %.3e, ‖ee‖/‖ex‖ ≤ %.3f, "
           "минус по углам у %d\n",
           wmom, wgrow, nneg);
    printf("    [НК прежний solve4] ‖ee‖/‖ex‖ до %.3e, минус по углам у %d из %d\n", wgrow_old,
           nneg_old, ncheck2);
    check(wnul < 1e-13,
          "К39: матрица масс плоского элемента вырождена ТОЧНО, нулевой вектор верен");
    check(wmom < 1e-9, "К39: проекция воспроизводит линейное поле — моменты совпали");
    check(wgrow < 1.0 + 1e-12, "К39: представитель МИНИМАЛЬНОЙ нормы, добавки вдоль nul нет");
    check(nneg == 0, "К39: у положительного поля откат к константе не нужен НИ РАЗУ");
    /* ПРЕДСКАЗАНИЕ «больше половины элементов уйдут в минус» НЕ СБЫЛОСЬ:
     * измерено 423 из 2212, то есть 19%. Причина в том, что здесь подаётся
     * ТОЧНО ЛИНЕЙНОЕ поле, у которого моменты представимы без остатка, и
     * величина паразитной добавки определяется только обусловленностью
     * конкретного элемента. В настоящей сцене облучённость линейной не бывает, и
     * там доля вышла 1536 из 1536 у стенок и 1372 из 1912 у поверхностей.
     * Существо контроля от этого не меняется и даже сильнее: норма представителя
     * раздувается до 1.4e17, то есть решение вырожденной системы бессмысленно
     * само по себе. Проверка переписана на ИЗМЕРЕННОЕ существо. */
    check(nneg_old > 100, "НК: прежний solve4 ОБЯЗАН уводить угловой минимум в минус");
    check(wgrow_old > 1e6, "НК: и раздувать норму представителя на порядки");
  }

  /* ------------------------------------------------------------- К47 ---
   * ФЛЮИДНАЯ ПЛОЩАДЬ НА ГРАНИ, ВТОРАЯ СТОРОНА КОТОРОЙ ПОЛНОСТЬЮ ТВЁРДАЯ.
   * У полностью твёрдой ячейки развёртка ставит L = 0 и идёт дальше, а её сосед
   * по флюиду всё равно считает выток через общую грань по СВОЕЙ флюидной
   * площади. Если та не ноль — энергия уходит в никуда, и это утечка, которую
   * баланс К40 обязан видеть. Геометрически площадь ОБЯЗАНА быть нулевой: если
   * соседняя ячейка целиком в материале, то общая грань лежит внутри материала,
   * и флюида на ней нет. Расхождение здесь есть расхождение ДВУХ ВЕЕРОВ
   * фасетов, отобранных для двух ячеек независимо. */
  {
    double leak = 0.0, total = 0.0;
    int nleak = 0;
    for (int32_t f = 0; f < m.nf; f++) {
      if (m.f[f].cb < 0) continue;
      total += cu.ffm[f][0][0] + cu.ffmb[f][0][0];
      int sa = cu.solid[m.f[f].ca], sb = cu.solid[m.f[f].cb];
      if (sa == sb) continue;
      double a = sa ? cu.ffmb[f][0][0] : cu.ffm[f][0][0]; /* площадь со стороны ФЛЮИДА */
      if (a > 0.0) {
        leak += a;
        nleak++;
      }
    }
    printf("  [К47] флюидная площадь на границе с ПОЛНОСТЬЮ ТВЁРДОЙ ячейкой: %.6e на %d гранях "
           "(вся площадь граней %.4f, доля %.3e)\n",
           leak, nleak, total, total > 0.0 ? leak / total : 0.0);
    check(!(leak > 0.0), "К47: флюид не граничит с полностью твёрдой ячейкой");
  }

  /* ------------------------------------------------------------- К48 ---
   * ДВЕ ПОЛУСФЕРНЫЕ СУММЫ У НАКЛОННОЙ НОРМАЛИ НЕ РАВНЫ. К29 велит нормировать
   * отражение на СОБСТВЕННУЮ сумму набора, а не на π из континуума. Но сумм
   * ДВЕ — по входящей полусфере и по исходящей, — и у произвольной нормали они
   * РАЗНЫЕ, потому что в полусферы попадают разные направления. Отданная
   * поверхностью мощность равна `Σ_{ω·n>0} w(ω·n) · ∫L_out`, значит делить надо
   * на ИСХОДЯЩУЮ сумму. У стенки, перпендикулярной оси, обе совпадают по
   * симметрии набора, и на равномерной сетке ошибка не проявляется ВООБЩЕ. */
  {
    tr3_dirs d2;
    if (tr3_dirs_product(&d2, 2, 2) == 0) {
      double worstr = 0.0;
      for (int32_t e = 0; e < cu.nse; e++) {
        double hin = 0.0, hout = 0.0;
        for (int mm = 0; mm < d2.n; mm++) {
          double on =
              d2.ox[mm] * cu.se[e].n[0] + d2.oy[mm] * cu.se[e].n[1] + d2.oz[mm] * cu.se[e].n[2];
          if (on < 0.0)
            hin += d2.w[mm] * (-on);
          else
            hout += d2.w[mm] * on;
        }
        if (hin > 0.0 && fabs(hout / hin - 1.0) > worstr) worstr = fabs(hout / hin - 1.0);
      }
      printf("  [К48] max |Σ_исх / Σ_вх − 1| по %d элементам: %.3e\n", cu.nse, worstr);
      tr3_dirs_free(&d2);
    }
  }

  /* ------------------------------------------------------------- К38 ---
   * СХОДИМОСТЬ ПРИ ОТРАЖАЮЩЕЙ ПОВЕРХНОСТИ ВНУТРИ ОБЛАСТИ. Прежде здесь было
   * автоколебание отката, и развёртка упиралась в maxit при невязке 1.25.
   * Ловится это ТОЛЬКО так: у полости точный ответ однороден, наклоны нулевые,
   * и переключателю не на чем колебаться. */
  {
    tr3_dirs d;
    check(tr3_dirs_product(&d, 2, 2) == 0, "ординаты");
    double *sig_t = calloc((size_t)m.ncell, sizeof(double));
    double *sig_s = calloc((size_t)m.ncell, sizeof(double));
    double *frho = calloc((size_t)ftab.n, sizeof(double));
    double *femit = calloc((size_t)ftab.n, sizeof(double));
    double *phi = calloc((size_t)m.ncell * 4, sizeof(double));
    if (sig_t == NULL || sig_s == NULL || frho == NULL || femit == NULL || phi == NULL) {
      check(0, "память");
    } else {
      for (int32_t c = 0; c < m.ncell; c++) {
        sig_t[c] = 0.015;
        sig_s[c] = 0.012;
      }
      for (int32_t i = 0; i < ftab.n; i++)
        frho[i] = 0.78;
      double wr[6] = {0.72, 0.35, 0.72, 0.72, 0.65, 0.05};
      double we[6] = {0, 0, 0, 0, 0, 6.0};
      tr3_problem p = {.m = &m,
                       .d = &d,
                       .cut = &cu,
                       .facet_rho = frho,
                       .facet_emit = femit,
                       .nfacet = ftab.n,
                       .sig_t = sig_t,
                       .sig_s = sig_s,
                       .wall_rho = wr,
                       .wall_emit = we,
                       .limiter = 1};
      /* СВИП ПО ДОПУСКУ: отставание итерации или НАСТОЯЩАЯ утечка?
       * Тождество К13 верно на любой итерации ДЛЯ ОБЪЁМА, но `sout` считается
       * из облучённости ПРЕДЫДУЩЕЙ итерации, поэтому поверхностный член отстаёт
       * на шаг. Если невязка падает вместе с допуском — это отставание и оно
       * уходит в пределе; если упирается в полку — это утечка, и её надо искать.
       * Различить их иначе нечем, и именно поэтому свип печатается (тот же приём,
       * что К22 у толщины слоя и Ф5 у числа отскоков). */
      for (int k = 0; k < 3; k++) {
        double tl = k == 0 ? 1e-7 : (k == 1 ? 1e-9 : 1e-11);
        tr3_stats s2;
        memset(phi, 0, (size_t)m.ncell * 4 * sizeof(double));
        if (tr3_sweep_solve(&p, 500, tl, phi, &s2) == 0)
          printf("        [К40 свип] допуск %.0e: итераций %3d, невязка баланса на втекшее %.3e\n",
                 tl, s2.iters, fabs(s2.balance) / s2.pin);
        free(s2.bout);
        free(s2.sout);
      }
      /* ЕДИНСТВЕННАЯ НЕЛИНЕЙНОСТЬ В ЦИКЛЕ — ОГРАНИЧИТЕЛЬ. Его замкнутая форма
       * (T2) обязана сохранять уравнение при `v = 1` ТОЧНО, то есть быть
       * консервативной. Проверяется прямо: с ним и без него. */
      for (int k = 0; k < 2; k++) {
        tr3_problem p2 = p;
        p2.limiter = k;
        tr3_stats s2;
        memset(phi, 0, (size_t)m.ncell * 4 * sizeof(double));
        if (tr3_sweep_solve(&p2, 500, 1e-10, phi, &s2) == 0)
          printf("        [К40 ограничитель %s] срезок %d, невязка баланса на втекшее %.3e\n",
                 k ? "ВКЛ" : "ВЫКЛ", s2.nclip, fabs(s2.balance) / s2.pin);
        free(s2.bout);
        free(s2.sout);
      }
      memset(phi, 0, (size_t)m.ncell * 4 * sizeof(double));
      tr3_stats st;
      int rc = tr3_sweep_solve(&p, 500, 1e-9, phi, &st);
      printf("  [К38] комната с ОТРАЖАЮЩЕЙ сферой: код %d, итераций %d, невязка %.2e, "
             "откатов проекции на последней итерации %d из %d\n",
             rc, st.iters, st.resid, st.nfallback, cu.nse + m.nf);
      /* К40: БАЛАНС С УЧЁТОМ ПОВЕРХНОСТЕЙ. Здесь ОДНО тело, то есть случая
       * «два тела в одной ячейке» нет по построению, и утечка, если она есть,
       * принадлежит самой схеме, а не построителю сцены. */
      printf("        [К40] баланс: втекло %.4f, вытекло %.4f, поглощено средой %.4f, "
             "ушло в поверхности %.4f, отдано %.4f; невязка на втекшее %.3e\n",
             st.pin, st.pout, st.pabs, st.psin, st.psout, fabs(st.balance) / st.pin);
      check(fabs(st.balance) / st.pin < 1e-9,
            "К40: баланс с поверхностями — точное дискретное тождество");
      check(rc == 0, "развёртка прошла");
      check(st.resid < 1e-9, "К38: невязка ДОСТИГЛА допуска, а не упёрлась в maxit");
      check(st.iters < 200, "К38: сходимость геометрическая — десятки итераций, а не тысячи");
      check(st.nfallback * 20 < cu.nse, "К38: откат проекции — редкое исключение, а не правило");
      free(st.bout);
      free(st.sout);
    }
    free(sig_t);
    free(sig_s);
    free(frho);
    free(femit);
    free(phi);
    tr3_dirs_free(&d);
  }

  tr3_cut_free(&cu);
  tr3_mesh_free(&m);
  hz_cutmap_free(&cmap);
  hz_facettab_free(&ftab);
  hz_surftab_free(&stab);
  hz_oct_free(&t);
}

/* --------------- ЗАМКНУТ ЛИ НАБОР ОРДИНАТ ОТНОСИТЕЛЬНО ЗЕРКАЛЬНОГО ОТРАЖЕНИЯ */

/* ОТ ЭТОГО ЗАВИСИТ, МОЖЕТ ЛИ ЗЕРКАЛО ИДТИ ЧЕРЕЗ РАЗВЁРТКУ ВОВСЕ.
 *
 * К5 записана как «зеркала через ординаты НЕ идут», и для КРИВОГО зеркала это
 * верно безоговорочно: нормаль меняется непрерывно вдоль поверхности, значит
 * отражений континуум, и ни один конечный набор под них не замкнут.
 *
 * Но для ПЛОСКОГО зеркала отражение есть ИЗОМЕТРИЯ пространства направлений, и
 * если набор переходит в себя, то зеркальная связь — ПЕРЕСТАНОВКА ординат.
 * Точная, разреженная, без интерполяции (то есть К3 не нарушается) и без
 * всякого прослеживания. Тогда бесконечная серия отражений между зеркалами
 * получается из РЕШЕНИЯ системы, а не из усечённого ряда с пределом отскоков.
 *
 * Проверяется здесь то, что нужно на деле: замкнутость относительно отражений
 * от стенок куба, то есть от плоскостей, перпендикулярных осям. Для них
 * отражение есть смена знака ОДНОЙ компоненты.
 *
 * НЕГАТИВНЫЙ КОНТРОЛЬ ЗАЛОЖЕН В ТОТ ЖЕ ТЕСТ: наклонная плоскость (нормаль по
 * диагонали) обязана замкнутость СЛОМАТЬ — иначе проверка ничего не значит,
 * потому что «замкнут под что угодно» невозможно. */
static void t_mirror_closure(void) {
  tr3_dirs d;
  if (tr3_dirs_product(&d, 2, 2)) {
    check(0, "набор направлений");
    return;
  }
  for (int ax = 0; ax < 3; ax++) {
    double worst = 0.0, wworst = 0.0;
    int missing = 0;
    for (int i = 0; i < d.n; i++) {
      double r[3] = {d.ox[i], d.oy[i], d.oz[i]};
      r[ax] = -r[ax]; /* отражение от плоскости, перпендикулярной оси ax */
      double best = 1e300;
      int bj = -1;
      for (int j = 0; j < d.n; j++) {
        double e = fabs(d.ox[j] - r[0]) + fabs(d.oy[j] - r[1]) + fabs(d.oz[j] - r[2]);
        if (e < best) {
          best = e;
          bj = j;
        }
      }
      if (best > worst) worst = best;
      if (bj >= 0 && fabs(d.w[bj] - d.w[i]) > wworst) wworst = fabs(d.w[bj] - d.w[i]);
      if (best > 1e-14) missing++;
    }
    printf(
        "  [ЗЕРКАЛО-ОСЬ %d] max промах отражённой ординаты %.3e, вес %.3e, вне набора %d из %d\n",
        ax, worst, wworst, missing, d.n);
    check(worst < 1e-14, "набор ЗАМКНУТ относительно отражения от стенки: связь есть перестановка");
    check(wworst < 1e-14, "и веса отражённых ординат совпадают");
  }
  /* ГРУППА СИММЕТРИИ НАБОРА ОКАЗАЛАСЬ ШИРЕ ОЖИДАЕМОГО, И ЭТО ЗАМЕР, А НЕ ДОГАДКА.
   * Первая редакция брала негативным контролем плоскость (1,1,0) — и он НЕ
   * ВЫСТРЕЛИЛ: набор под неё тоже замкнут, потому что азимутальные узлы
   * симметричны относительно 45°, и вся конструкция инвариантна относительно
   * ОКТАЭДРИЧЕСКОЙ группы, то есть трёх осевых плоскостей И ШЕСТИ ДИАГОНАЛЬНЫХ.
   * Девять точно представимых зеркальных плоскостей вместо трёх — это подарок,
   * но проверять надо на нормали, которая симметрией НЕ является. */
  {
    double n0[3] = {0.6, 0.48, 0.64}, nn = 0.0;
    for (int a = 0; a < 3; a++)
      nn += n0[a] * n0[a];
    nn = sqrt(nn);
    double n[3] = {n0[0] / nn, n0[1] / nn, n0[2] / nn}, worst = 0.0;
    int missing = 0;
    for (int i = 0; i < d.n; i++) {
      double o[3] = {d.ox[i], d.oy[i], d.oz[i]};
      double dn = o[0] * n[0] + o[1] * n[1] + o[2] * n[2], r[3];
      for (int a = 0; a < 3; a++)
        r[a] = o[a] - 2.0 * dn * n[a];
      double best = 1e300;
      for (int j = 0; j < d.n; j++) {
        double e = fabs(d.ox[j] - r[0]) + fabs(d.oy[j] - r[1]) + fabs(d.oz[j] - r[2]);
        if (e < best) best = e;
      }
      if (best > worst) worst = best;
      if (best > 1e-14) missing++;
    }
    printf("    [НК несимметричная нормаль (0.6,0.48,0.64)] max промах %.3e, вне набора %d из %d\n",
           worst, missing, d.n);
    check(missing > 0, "НК: под НАКЛОННОЕ зеркало набор замкнут быть НЕ должен");
  }
  tr3_dirs_free(&d);
}

/* ------------------------------------- ФОРМАТ ВЫХОДА: PFM, круговой прогон */

/* Писалки картинок врут МОЛЧА и ровно в трёх местах: порядок строк, порядок
 * байт и чередование каналов. Ни одно из трёх не видно на картинке, если
 * смотреть её тем же кодом, который писал. Поэтому файл читается ОБРАТНО
 * сырыми байтами и сверяется с тем, что клали. */
static void t_pfm(void) {
  const int W = 7, H = 5;
  double buf[35];
  for (int y = 0; y < H; y++)
    for (int x = 0; x < W; x++)
      buf[y * W + x] = 0.125 * (double)(y * W + x) + 0.0625; /* точно представимо во float */
  const char *path = "img/_pfm_roundtrip.pfm";
  check(hz_pfm_write(path, buf, buf, buf, W, H) == 0, "PFM записан");
  FILE *f = fopen(path, "rb");
  if (f == NULL) {
    check(0, "PFM открылся на чтение");
    return;
  }
  int w2 = 0, h2 = 0;
  double sc = 0.0;
  char magic[3] = {0, 0, 0};
  int ok_hdr = fscanf(f, "%2s %d %d %lf", magic, &w2, &h2, &sc) == 4;
  fgetc(f); /* один перевод строки после масштаба — часть формата */
  printf("  [PFM] заголовок: «%s» %dx%d масштаб %.1f\n", magic, w2, h2, sc);
  check(ok_hdr && magic[0] == 'P' && magic[1] == 'F', "PFM: три канала (PF), а не Pf");
  check(w2 == W && h2 == H, "PFM: размеры");
  check(sc < 0.0, "PFM: масштаб отрицателен — little-endian, как на этой машине");
  double worst = 0.0;
  int bad_order = 0;
  for (int r = 0; r < H; r++) {
    /* СТРОКИ В ФАЙЛЕ ИДУТ СНИЗУ ВВЕРХ: строка r файла есть строка H−1−r буфера */
    int y = H - 1 - r;
    for (int x = 0; x < W; x++) {
      float px[3];
      if (fread(px, sizeof(float), 3, f) != 3) {
        bad_order = 1;
        break;
      }
      float want = (float)buf[y * W + x];
      for (int c = 0; c < 3; c++) {
        double e = fabs((double)px[c] - (double)want);
        if (e > worst) worst = e;
      }
    }
  }
  long extra = 0;
  while (fgetc(f) != EOF)
    extra++;
  fclose(f);
  printf("  [PFM] круговой прогон: max |прочитано − записано| = %.3e, лишних байт %ld\n", worst,
         extra);
  check(!bad_order, "PFM: файл не оборван");
  check(!(worst > 0.0), "PFM: значения совпали ПОБИТОВО (порядок строк и байт верны)");
  check(extra == 0, "PFM: ни одного лишнего байта");
  remove(path); /* временный файл теста, а не картинка на посмотреть */
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
  t_linear_field();
  t_balance();
  t_cavity();
  t_cut();
  t_mirror_closure();
  t_pfm();
  t_render();

  printf("%s: %d/%d\n", g_fail ? "FAILURES" : "ok", g_total - g_fail, g_total);
  return g_fail != 0;
}
