/* test_rte2d.c — приёмка вехи T2 (PLAN_TRANSPORT.md): рассеяние, печная
 * задача, DG1, ограничитель. ПРЕДСКАЗАНИЯ ЗАПИСАНЫ ДО ПРОГОНА и печатаются
 * перед результатами.
 *
 * ПЕЧЬ ОДНА НЕ ЗАКРЫВАЕТ T2, и это выяснилось при постановке тестов:
 *  - однородный радианс имеет НУЛЕВОЙ наклон, поэтому ошибка в уравнениях для
 *    наклонных степеней свободы печь проходит незамеченной ⇒ нужен тест с
 *    ЛИНЕЙНЫМ точным решением;
 *  - передачу потока через разбитую грань (К8) печь тоже не ловит: любой
 *    разумный оператор передачи воспроизводит КОНСТАНТУ ⇒ нужен баланс энергии
 *    при НЕОДНОРОДНОЙ подсветке.
 * Оба добавлены рядом с печью.
 */
#include "transport/rte2d.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fail = 0;
static int g_total = 0;

static void check(int ok, const char *what) {
  g_total++;
  if (!ok) {
    g_fail++;
    printf("  FAIL: %s\n", what);
  }
}

/* Порог «точно»: величины получаются из линейных решений 3x3 и суммирования по
 * ячейкам и ординатам, то есть пол — округление, ~sqrt(число операций)·eps.
 * При тысячах ячеек и 32 ординатах это несколько 1e-14. Цель плана 1e-14. */
#define TOL_EXACT 5e-14
/* Сходимость итерации по рассеянию: относительное изменение φ. 1e-15 — пол
 * округления; ниже итерация не остановится. Ошибка установившегося решения
 * тогда ~tol/(1−ρ), при альбедо 1 и утечке ρ~0.7 это ~3e-15. */
#define TOL_ITER 1e-15
#define TOL_CONVERGED 1e-12

enum { NPER = 8 }; /* 8 узлов Гаусса на квадрант ⇒ 32 направления */

/* --------------------------------------------------------------- сетки */

typedef struct {
  const char *name;
  tr2_mesh m;
} meshcase;

/* Три сетки по требованию аудита К8: на равномерной ошибка передачи потока
 * между уровнями не проявится ВООБЩЕ. nunit=16: единица 1/16 стороны. */
static int build_meshes(meshcase *mc) {
  const double H = 1.0;
  mc[0].name = "равномерная 16x16";
  if (tr2_mesh_uniform(&mc[0].m, 16, 16, H / 16, H / 16) != 0) return -1;
  mc[1].name = "градуированная (квадродерево к центру)";
  if (tr2_mesh_graded(&mc[1].m, 16, H / 16, H / 16) != 0) return -1;
  mc[2].name = "анизотропная (4x1 слева, 1x4 справа, ux=2uy)";
  if (tr2_mesh_aniso(&mc[2].m, 16, 2.0 * H / 16, H / 16, 4) != 0) return -1;
  return 0;
}

/* ------------------------------------------------- вспомогательное */

/* Заполняет ε(x,ω) = (σ_t(x,ω) − σ_s(x))·L_b — согласованное излучение, при
 * котором ОДНОРОДНЫЙ радианс L_b есть точное решение при любом альбедо. При
 * альбедо 1 это ε = 0, то есть буквально печь из плана. */
static void fill_emit_consistent(double *emit, const tr2_mesh *m, const tr2_dirs *d,
                                 const double *sig_t, const double *sig_s, const double *gdir,
                                 double lb) {
  for (int c = 0; c < m->ncell; c++)
    for (int mi = 0; mi < d->n; mi++) {
      double gm = gdir ? gdir[mi] : 1.0;
      double *e = &emit[((size_t)c * (size_t)d->n + (size_t)mi) * 3];
      e[0] = (gm * sig_t[c] - sig_s[c]) * lb;
      e[1] = 0.0;
      e[2] = 0.0;
    }
}

int main(void) {
  printf("=== T2: рассеяние, печь, DG1, ограничитель ===\n\n");
  printf("ПРЕДСКАЗАНИЯ (записаны до прогона):\n");
  printf("  1. Sw = 2pi; p0*Sw = 1 РОВНО (p0 = 1/Sw, не 1/2pi)        <= %.0e\n", TOL_EXACT);
  printf("  2. сетки: замощение целочисленно точное; разбитых граней\n");
  printf("     БОЛЬШЕ НУЛЯ на градуированной и анизотропной, РОВНО НОЛЬ\n");
  printf("     на равномерной — иначе К8 нечем проверять\n");
  printf("  3. ПЕЧЬ, альбедо 1, eps=0, ограничитель ВКЛЮЧЁН (К6/К11):\n");
  printf("     невязка точной константы (один проход)               <= %.0e\n", TOL_EXACT);
  printf("     сходимость от нуля                                   <= %.0e\n", TOL_CONVERGED);
  printf("     на ВСЕХ ТРЁХ сетках, и число итераций от сетки почти не зависит\n");
  printf("  4. ПЕЧЬ общая: альбедо 0.5, НЕОДНОРОДНОЕ sigma_t,\n");
  printf("     eps = sigma_a*L_b ⇒ L = L_b точно                     <= %.0e\n", TOL_EXACT);
  printf("  5. ЛИНЕЙНОЕ точное решение L = A + Bx + Cy (проверяет НАКЛОНЫ,\n");
  printf("     чего печь не делает вовсе), 3 сетки                   <= %.0e\n", TOL_EXACT);
  printf("     ограничитель при этом НЕ ДОЛЖЕН срабатывать ни разу\n");
  printf("  6. ЛОВУШКА К1: направленное sigma_t + ИЗОТРОПНОЕ eps ОБЯЗАНО\n");
  printf("     поехать (>= 1e-3) у ИСПРАВНОГО кода; с согласованным\n");
  printf("     eps(omega) — точно                                    <= %.0e\n", TOL_EXACT);
  printf("  7. БАЛАНС ЭНЕРГИИ (точное дискретное тождество), альбедо 1,\n");
  printf("     подсветка ТОЛЬКО слева, 3 сетки                       <= %.0e\n", TOL_EXACT);
  printf("  8. ОГРАНИЧИТЕЛЬ: без него в задаче с тенью радианс уходит\n");
  printf("     в минус (nneg > 0); CONSERVATIVE убирает минус, срабатывает\n");
  printf("     (nlimited > 0) и СОХРАНЯЕТ баланс                     <= %.0e\n", TOL_EXACT);
  printf("     негативный контроль: режим SLOPE (без пересчёта среднего)\n");
  printf("     ОБЯЗАН сломать баланс, >= 1e-4\n");
  printf("  9. НЕГАТИВНЫЙ КОНТРОЛЬ ПЛАНА: нормировка фазовой функции x1.01.\n");
  printf("     Замкнутая форма: L/L_b = sigma_a/(sigma_a - d*sigma_s), при\n");
  printf("     альбедо 0.5 и d=0.01 это 1/0.99 ⇒ ОТКЛОНЕНИЕ 1.0101e-2,\n");
  printf("     а НЕ 1e-2. При альбедо 1 контроль вырожден (0/0) — поэтому\n");
  printf("     гоняется при 0.5, и в ОПТИЧЕСКИ ТОЛСТОЙ области: замкнутая форма\n");
  printf("     есть равновесие бесконечной среды, а в тонкой подсветка на\n");
  printf("     границе держит решение около L_b (первый прогон дал 7.7e-3 при\n");
  printf("     tau~1..5 — предсказание было неверно, поправлено)\n");
  printf(" 10. DG1 на НЕполиномиальном решении: сходится (Ричардсон), порядок\n");
  printf("     ~3 = 2p+1 — суперсходимость ФУНКЦИОНАЛА (вытекшей мощности).\n");
  printf("     Первое предсказание было «между 1 и 2, изломы по диагоналям\n");
  printf("     испортят порядок» — НЕ СБЫЛОСЬ, измерено 2.90: изломы лежат на\n");
  printf("     множестве меры нуль и интегральный функционал их не чувствует\n\n");

  tr2_dirs d;
  if (tr2_dirs_quadrant_gauss(&d, NPER) != 0) {
    printf("FATAL: dirs\n");
    return 1;
  }
  tr2_phase ph;
  tr2_phase_isotropic(&ph, &d);

  /* ---------------------------------------------------------------- 1 */
  printf("[1] направления: %d (по %d на квадрант)\n", d.n, NPER);
  {
    double sw = 0.0, sx = 0.0, sy = 0.0, minc = 1.0;
    for (int m = 0; m < d.n; m++) {
      sw += d.w[m];
      sx += d.w[m] * d.ox[m];
      sy += d.w[m] * d.oy[m];
      double a = fabs(d.ox[m]) < fabs(d.oy[m]) ? fabs(d.ox[m]) : fabs(d.oy[m]);
      if (a < minc) minc = a;
    }
    printf("    Sw-2pi = %.2e   Sw*ox = %.2e   Sw*oy = %.2e   p0*Sw-1 = %.2e\n",
           fabs(sw - 2 * M_PI), fabs(sx), fabs(sy), fabs(ph.p0 * sw - 1.0));
    printf("    min|компонента omega| = %.4f — направлений вдоль осей нет\n", minc);
    check(fabs(sw - 2 * M_PI) < TOL_EXACT, "Sw = 2pi");
    check(fabs(sx) < TOL_EXACT && fabs(sy) < TOL_EXACT, "first angular moment vanishes");
    check(fabs(ph.p0 * sw - 1.0) < TOL_EXACT, "p0*Sw = 1 exactly");
    check(minc > 1e-3, "no axis-parallel ordinate (no degenerate zero-outflow face)");
  }

  meshcase mc[3];
  if (build_meshes(mc) != 0) {
    printf("FATAL: meshes\n");
    return 1;
  }

  /* ---------------------------------------------------------------- 2 */
  printf("\n[2] сетки (замощение проверено целочисленно внутри tr2_mesh_finish)\n");
  for (int k = 0; k < 3; k++) {
    int ns = tr2_mesh_nsplit(&mc[k].m);
    printf("    %-42s ячеек %5d  граней %5d  разбитых %4d\n", mc[k].name, mc[k].m.ncell,
           mc[k].m.nface, ns);
    if (k == 0)
      check(ns == 0, "uniform mesh has no split faces (control)");
    else
      check(ns > 0, "graded/anisotropic mesh really does have split faces");
  }

  /* ---------------------------------------------------------------- 3 */
  printf("\n[3] ПЕЧЬ: альбедо 1, eps = 0, подсветка L_b со всех сторон, лимитер ON\n");
  const double LB = 2.5;
  for (int k = 0; k < 3; k++) {
    const tr2_mesh *m = &mc[k].m;
    double *st_ = calloc((size_t)m->ncell, sizeof *st_);
    double *ss_ = calloc((size_t)m->ncell, sizeof *ss_);
    double *phi0 = calloc((size_t)m->ncell * 3, sizeof *phi0);
    if (!st_ || !ss_ || !phi0) return 1;
    for (int c = 0; c < m->ncell; c++) {
      st_[c] = 3.0; /* tau по области ~3 */
      ss_[c] = 3.0; /* альбедо ровно 1 */
      phi0[3 * c] = 2.0 * M_PI * LB;
    }
    tr2_problem pr;
    memset(&pr, 0, sizeof pr);
    pr.mesh = m;
    pr.dirs = &d;
    pr.phase = &ph;
    pr.sig_t = st_;
    pr.sig_s = ss_;
    pr.lim = TR2_LIM_CONSERVATIVE;
    pr.tol = TOL_ITER;
    for (int s = 0; s < 4; s++)
      pr.lb[s] = LB;

    /* (а) невязка точного решения: один проход из точной константы */
    pr.one_pass = 1;
    pr.phi0 = phi0;
    tr2_stats a;
    if (tr2_solve(&pr, &a) != 0) return 1;
    double e_res = fmax(fabs(a.lmax / LB - 1.0), fabs(a.lmin / LB - 1.0));

    /* (б) сходимость от нуля */
    pr.one_pass = 0;
    pr.phi0 = NULL;
    tr2_stats b;
    if (tr2_solve(&pr, &b) != 0) return 1;
    double e_cnv = fmax(fabs(b.lmax / LB - 1.0), fabs(b.lmin / LB - 1.0));

    printf("    %-42s невязка %.2e   от нуля %.2e  (%d итер)\n", mc[k].name, e_res, e_cnv, b.iters);
    check(e_res < TOL_EXACT, "furnace: exact constant is a fixed point");
    check(e_cnv < TOL_CONVERGED, "furnace: converges to L_b from zero");
    check(a.nlimited == 0, "limiter does not touch the furnace (constant has no slope)");
    tr2_stats_free(&a);
    tr2_stats_free(&b);
    free(st_);
    free(ss_);
    free(phi0);
  }

  /* ---------------------------------------------------------------- 4, 9 */
  printf("\n[4] ПЕЧЬ общая: альбедо 0.5, НЕОДНОРОДНОЕ sigma_t, eps = sigma_a*L_b\n");
  double ctrl_dev = 0.0, ctrl_pred = 0.0;
  for (int k = 0; k < 3; k++) {
    const tr2_mesh *m = &mc[k].m;
    double *st_ = calloc((size_t)m->ncell, sizeof *st_);
    double *ss_ = calloc((size_t)m->ncell, sizeof *ss_);
    double *em = calloc((size_t)m->ncell * (size_t)d.n * 3, sizeof *em);
    if (!st_ || !ss_ || !em) return 1;
    for (int c = 0; c < m->ncell; c++) {
      st_[c] = 1.0 + 4.0 * (double)(c % 5) / 4.0; /* от 1 до 5, неоднородно */
      ss_[c] = 0.5 * st_[c];                      /* альбедо ровно 0.5 */
    }
    fill_emit_consistent(em, m, &d, st_, ss_, NULL, LB);
    tr2_problem pr;
    memset(&pr, 0, sizeof pr);
    pr.mesh = m;
    pr.dirs = &d;
    pr.phase = &ph;
    pr.sig_t = st_;
    pr.sig_s = ss_;
    pr.emit = em;
    pr.lim = TR2_LIM_CONSERVATIVE;
    pr.tol = TOL_ITER;
    for (int s = 0; s < 4; s++)
      pr.lb[s] = LB;
    tr2_stats a;
    if (tr2_solve(&pr, &a) != 0) return 1;
    double e = fmax(fabs(a.lmax / LB - 1.0), fabs(a.lmin / LB - 1.0));
    printf("    %-42s ошибка %.2e  (%d итер)\n", mc[k].name, e, a.iters);
    check(e < TOL_EXACT, "general furnace: constant solution exact at albedo 0.5");
    tr2_stats_free(&a);

    free(st_);
    free(ss_);
    free(em);
  }

  /* --- негативный контроль плана: нормировка фазовой функции --- */
  double ctrl_edge = 0.0;
  {
    /* ОПТИЧЕСКИ ТОЛСТАЯ среда, и мерить надо В ГЛУБИНЕ. Замкнутая форма
     * L/L_b = σ_a/(σ_a − δσ_s) есть равновесие БЕСКОНЕЧНОЙ среды, поэтому:
     *  - в тонкой области подсветка на границе держит решение около L_b и
     *    отклонение выходит МЕНЬШЕ асимптотики (измерено 7.7e-3 при τ~1..5);
     *  - а максимум по ВСЕМ ячейкам берёт ПЕРЕСВЕТ пограничного слоя: при
     *    τ на ячейку 6.25 DG1 в слое немонотонен и даёт 1.041 вместо 1.0101
     *    (ограничитель это не режет — он следит только за положительностью).
     * Поэтому: сетка 64x64 (τ на ячейку 0.47, слой разрешён), τ по области 30
     * (влияние границы e^{−15} ~ 3e-7), и отсчёт в ЦЕНТРАЛЬНОЙ ячейке. */
    tr2_mesh cm;
    if (tr2_mesh_uniform(&cm, 64, 64, 1.0 / 64, 1.0 / 64) != 0) return 1;
    const tr2_mesh *m = &cm;
    const double SIG = 30.0, DELTA = 0.01;
    double *st_ = calloc((size_t)m->ncell, sizeof *st_);
    double *ss_ = calloc((size_t)m->ncell, sizeof *ss_);
    double *em = calloc((size_t)m->ncell * (size_t)d.n * 3, sizeof *em);
    if (!st_ || !ss_ || !em) return 1;
    for (int c = 0; c < m->ncell; c++) {
      st_[c] = SIG;
      ss_[c] = 0.5 * SIG;
    }
    fill_emit_consistent(em, m, &d, st_, ss_, NULL, LB);
    tr2_problem pr;
    memset(&pr, 0, sizeof pr);
    pr.mesh = m;
    pr.dirs = &d;
    pr.phase = &ph;
    pr.sig_t = st_;
    pr.sig_s = ss_;
    pr.emit = em;
    pr.lim = TR2_LIM_CONSERVATIVE;
    pr.tol = TOL_ITER;
    for (int s = 0; s < 4; s++)
      pr.lb[s] = LB;
    /* центральная ячейка: глубина по обеим осям ~ τ/2 */
    double sw = 0.0;
    for (int mi = 0; mi < d.n; mi++)
      sw += d.w[mi];
    int cc = -1;
    for (int c = 0; c < m->ncell; c++)
      if (m->cell[c].ix0 == m->nx / 2 && m->cell[c].iy0 == m->ny / 2) cc = c;
    if (cc < 0) return 1;

    tr2_stats good;
    if (tr2_solve(&pr, &good) != 0) return 1;
    double e_good = fmax(fabs(good.lmax / LB - 1.0), fabs(good.lmin / LB - 1.0));
    tr2_stats_free(&good);
    tr2_phase bad = ph;
    bad.p0 *= 1.0 + DELTA;
    pr.phase = &bad;
    tr2_stats c;
    if (tr2_solve(&pr, &c) != 0) return 1;
    ctrl_pred = 1.0 / (1.0 - DELTA) - 1.0; /* альбедо 0.5 ⇒ σ_a/(σ_a−δσ_s) */
    ctrl_dev = c.phi[3 * cc] / sw / LB - 1.0;
    ctrl_edge = c.lmax / LB - 1.0;
    tr2_stats_free(&c);
    printf("    толстая (tau=%.0f, 64x64), НЕвозмущённая: ошибка %.2e\n", SIG, e_good);
    check(e_good < TOL_EXACT, "thick general furnace exact before the control");
    free(st_);
    free(ss_);
    free(em);
    tr2_mesh_free(&cm);
  }
  printf("\n[9] НЕГАТИВНЫЙ КОНТРОЛЬ ПЛАНА: нормировка фазовой функции x1.01\n");
  printf("    в ГЛУБИНЕ: отклонение %.8e, замкнутая форма %.8e, разница %.2e\n", ctrl_dev,
         ctrl_pred, fabs(ctrl_dev - ctrl_pred));
  printf("    максимум по всей области %.4e — здесь совпал с глубинным, потому что\n", ctrl_edge);
  printf("    tau на ячейку 0.47 и слой разрешён; при 6.25 он давал 1.041e-2, то есть\n");
  printf("    ПЕРЕСВЕТ слоя, а его ограничитель не режет — он следит только за минусом\n");
  check(fabs(ctrl_dev - ctrl_pred) < 1e-5, "phase x1.01 => furnace moves by exactly 1/0.99-1");

  /* ---------------------------------------------------------------- 5 */
  printf("\n[5] ЛИНЕЙНОЕ точное решение L = A + Bx + Cy — проверка НАКЛОНОВ\n");
  {
    const double A = 10.0, B = 1.5, C = -0.8, SIG = 2.0;
    for (int k = 0; k < 3; k++) {
      const tr2_mesh *m = &mc[k].m;
      double *st_ = calloc((size_t)m->ncell, sizeof *st_);
      double *ss_ = calloc((size_t)m->ncell, sizeof *ss_);
      double *em = calloc((size_t)m->ncell * (size_t)d.n * 3, sizeof *em);
      if (!st_ || !ss_ || !em) return 1;
      for (int c = 0; c < m->ncell; c++) {
        st_[c] = SIG;
        ss_[c] = 0.0; /* без рассеяния: проверяем чистую пространственную схему */
      }
      /* eps = omega.grad L + sigma_t*L, DG1: среднее и наклоны в локальных
       * координатах (наклон по xi равен dL/dx * hx). */
      for (int c = 0; c < m->ncell; c++) {
        const tr2_cell *kk = &m->cell[c];
        double hx = (double)(kk->ix1 - kk->ix0) * m->ux, hy = (double)(kk->iy1 - kk->iy0) * m->uy;
        double xc = 0.5 * (double)(kk->ix0 + kk->ix1) * m->ux;
        double yc = 0.5 * (double)(kk->iy0 + kk->iy1) * m->uy;
        double lmean = A + B * xc + C * yc;
        for (int mi = 0; mi < d.n; mi++) {
          double *e = &em[((size_t)c * (size_t)d.n + (size_t)mi) * 3];
          e[0] = B * d.ox[mi] + C * d.oy[mi] + SIG * lmean;
          e[1] = SIG * B * hx;
          e[2] = SIG * C * hy;
        }
      }
      tr2_problem pr;
      memset(&pr, 0, sizeof pr);
      pr.mesh = m;
      pr.dirs = &d;
      pr.phase = &ph;
      pr.sig_t = st_;
      pr.sig_s = ss_;
      pr.emit = em;
      pr.lim = TR2_LIM_CONSERVATIVE;
      pr.tol = TOL_ITER;
      /* влёт = то же линейное поле на каждой стороне */
      pr.lb[0] = A + B * 0.0;
      pr.lbt[0] = C; /* xlo: L = A + C*y */
      pr.lb[1] = A + B * (double)m->nx * m->ux;
      pr.lbt[1] = C;
      pr.lb[2] = A + C * 0.0;
      pr.lbt[2] = B; /* ylo: L = A + B*x */
      pr.lb[3] = A + C * (double)m->ny * m->uy;
      pr.lbt[3] = B;
      tr2_stats a;
      if (tr2_solve(&pr, &a) != 0) return 1;
      /* сравниваем phi: точное phi = 2pi*L(x,y), среднее по ячейке = 2pi*lmean */
      double worst = 0.0;
      for (int c = 0; c < m->ncell; c++) {
        const tr2_cell *kk = &m->cell[c];
        double xc = 0.5 * (double)(kk->ix0 + kk->ix1) * m->ux;
        double yc = 0.5 * (double)(kk->iy0 + kk->iy1) * m->uy;
        double hx = (double)(kk->ix1 - kk->ix0) * m->ux, hy = (double)(kk->iy1 - kk->iy0) * m->uy;
        double sw = 0.0;
        for (int mi = 0; mi < d.n; mi++)
          sw += d.w[mi];
        double r0 = sw * (A + B * xc + C * yc), r1 = sw * B * hx, r2 = sw * C * hy;
        double e0 = fabs(a.phi[3 * c] / r0 - 1.0);
        double e1 = fabs((a.phi[3 * c + 1] - r1) / r0);
        double e2 = fabs((a.phi[3 * c + 2] - r2) / r0);
        if (e0 > worst) worst = e0;
        if (e1 > worst) worst = e1;
        if (e2 > worst) worst = e2;
      }
      printf("    %-42s ошибка (среднее И наклоны) %.2e  срезов %ld\n", mc[k].name, worst,
             a.nlimited);
      check(worst < TOL_EXACT, "linear exact solution reproduced exactly, slopes included");
      check(a.nlimited == 0, "limiter does not fire on a smooth positive field");
      tr2_stats_free(&a);
      free(st_);
      free(ss_);
      free(em);
    }
  }

  /* ---------------------------------------------------------------- 6 */
  printf("\n[6] ЛОВУШКА К1: направленное sigma_t = G(omega)*u_L\n");
  {
    const tr2_mesh *m = &mc[0].m;
    double *gd = calloc((size_t)d.n, sizeof *gd);
    double *st_ = calloc((size_t)m->ncell, sizeof *st_);
    double *ss_ = calloc((size_t)m->ncell, sizeof *ss_);
    double *em = calloc((size_t)m->ncell * (size_t)d.n * 3, sizeof *em);
    if (!gd || !st_ || !ss_ || !em) return 1;
    double gbar = 0.0, wsum = 0.0;
    for (int mi = 0; mi < d.n; mi++) {
      gd[mi] = 1.0 + 2.0 * fabs(d.oy[mi]); /* «стебли»: поперёк видно больше */
      gbar += d.w[mi] * gd[mi];
      wsum += d.w[mi];
    }
    gbar /= wsum;
    for (int c = 0; c < m->ncell; c++) {
      st_[c] = 2.0;
      ss_[c] = 1.0; /* альбедо ~0.5 при gbar~1.6 */
    }
    tr2_problem pr;
    memset(&pr, 0, sizeof pr);
    pr.mesh = m;
    pr.dirs = &d;
    pr.phase = &ph;
    pr.sig_t = st_;
    pr.sig_s = ss_;
    pr.gdir = gd;
    pr.emit = em;
    pr.lim = TR2_LIM_CONSERVATIVE;
    pr.tol = TOL_ITER;
    for (int s = 0; s < 4; s++)
      pr.lb[s] = LB;

    /* (а) НАИВНАЯ печь: eps изотропное, по среднему G — ОБЯЗАНА поехать */
    for (int c = 0; c < m->ncell; c++)
      for (int mi = 0; mi < d.n; mi++) {
        double *e = &em[((size_t)c * (size_t)d.n + (size_t)mi) * 3];
        e[0] = (gbar * st_[c] - ss_[c]) * LB;
        e[1] = e[2] = 0.0;
      }
    tr2_stats a;
    if (tr2_solve(&pr, &a) != 0) return 1;
    double e_naive = fmax(fabs(a.lmax / LB - 1.0), fabs(a.lmin / LB - 1.0));
    tr2_stats_free(&a);

    /* (б) СОГЛАСОВАННОЕ eps(omega) — обязана быть точной */
    fill_emit_consistent(em, m, &d, st_, ss_, gd, LB);
    tr2_stats b;
    if (tr2_solve(&pr, &b) != 0) return 1;
    double e_cons = fmax(fabs(b.lmax / LB - 1.0), fabs(b.lmin / LB - 1.0));
    tr2_stats_free(&b);

    printf("    изотропное eps (наивно): %.3e   согласованное eps(omega): %.3e\n", e_naive, e_cons);
    printf("    G в [%.2f, %.2f], среднее %.3f\n", 1.0, 3.0, gbar);
    check(e_naive > 1e-3, "K1 trap fires: directional sigma_t + isotropic eps MUST deviate");
    check(e_cons < TOL_EXACT, "directional sigma_t with consistent eps(omega) is exact");
    free(gd);
    free(st_);
    free(ss_);
    free(em);
  }

  /* ---------------------------------------------------------------- 7, 8 */
  printf("\n[7] БАЛАНС ЭНЕРГИИ: альбедо 1, подсветка ТОЛЬКО слева\n");
  for (int k = 0; k < 3; k++) {
    const tr2_mesh *m = &mc[k].m;
    double *st_ = calloc((size_t)m->ncell, sizeof *st_);
    double *ss_ = calloc((size_t)m->ncell, sizeof *ss_);
    if (!st_ || !ss_) return 1;
    for (int c = 0; c < m->ncell; c++) {
      st_[c] = 2.0;
      ss_[c] = 2.0; /* альбедо 1: вся влитая энергия обязана выйти */
    }
    tr2_problem pr;
    memset(&pr, 0, sizeof pr);
    pr.mesh = m;
    pr.dirs = &d;
    pr.phase = &ph;
    pr.sig_t = st_;
    pr.sig_s = ss_;
    pr.lim = TR2_LIM_CONSERVATIVE;
    pr.tol = TOL_ITER;
    pr.lb[0] = LB; /* только xlo */
    tr2_stats a;
    if (tr2_solve(&pr, &a) != 0) return 1;
    double rel = fabs(a.balance) / a.pin;
    printf("    %-42s in %.6f  out %.6f  |баланс|/in %.2e (%d итер)\n", mc[k].name, a.pin, a.pout,
           rel, a.iters);
    check(rel < TOL_EXACT, "energy balance exact, incl. flux transfer across split faces");
    tr2_stats_free(&a);
    free(st_);
    free(ss_);
  }

  printf("\n[8] ОГРАНИЧИТЕЛЬ: задача с тенью (непрозрачный блок в центре)\n");
  {
    const tr2_mesh *m = &mc[0].m;
    double *st_ = calloc((size_t)m->ncell, sizeof *st_);
    double *ss_ = calloc((size_t)m->ncell, sizeof *ss_);
    if (!st_ || !ss_) return 1;
    for (int c = 0; c < m->ncell; c++) {
      const tr2_cell *kk = &m->cell[c];
      int in = kk->ix0 >= 6 && kk->ix1 <= 10 && kk->iy0 >= 6 && kk->iy1 <= 10;
      st_[c] = in ? 400.0 : 0.02; /* блок оптически толстый, вокруг почти пусто */
      ss_[c] = 0.0;
    }
    tr2_problem pr;
    memset(&pr, 0, sizeof pr);
    pr.mesh = m;
    pr.dirs = &d;
    pr.phase = &ph;
    pr.sig_t = st_;
    pr.sig_s = ss_;
    pr.tol = TOL_ITER;
    pr.lb[0] = LB;
    pr.lb[2] = LB;

    tr2_stats off, cons, slope;
    pr.lim = TR2_LIM_OFF;
    if (tr2_solve(&pr, &off) != 0) return 1;
    pr.lim = TR2_LIM_CONSERVATIVE;
    if (tr2_solve(&pr, &cons) != 0) return 1;
    pr.lim = TR2_LIM_SLOPE;
    if (tr2_solve(&pr, &slope) != 0) return 1;

    double b_off = fabs(off.balance) / off.pin;
    double b_cons = fabs(cons.balance) / cons.pin;
    double b_slope = fabs(slope.balance) / slope.pin;
    printf("    OFF          : min угол/L_b %+.3e  отрицательных %5ld  |баланс|/in %.2e\n",
           off.cmin / LB, off.nneg, b_off);
    printf("    CONSERVATIVE : min угол/L_b %+.3e  отрицательных %5ld  |баланс|/in %.2e"
           "  срезов %ld\n",
           cons.cmin / LB, cons.nneg, b_cons, cons.nlimited);
    printf("    SLOPE (контр): min угол/L_b %+.3e  отрицательных %5ld  |баланс|/in %.2e\n",
           slope.cmin / LB, slope.nneg, b_slope);
    check(off.nneg > 0, "unlimited DG1 really does go negative at a shadow edge (the T5 failure)");
    check(b_off < TOL_EXACT, "unlimited DG1 is conservative (baseline for the control)");
    check(cons.nlimited > 0, "limiter actually fires — otherwise the test proves nothing");
    /* Срезка ставит минимум по углам РОВНО в нуль, поэтому оставшийся минус —
     * это округление около нуля, а не ошибка схемы. Считается он относительно
     * L_b, натурального масштаба задачи, а не относительно самого нуля. */
    check(cons.cmin / LB > -TOL_EXACT, "conservative limiter removes negative radiance");
    check(b_cons < TOL_EXACT, "conservative limiter keeps the balance exact");
    check(b_slope > 1e-4, "NEGATIVE CONTROL: slope-only limiting MUST break the balance");
    tr2_stats_free(&off);
    tr2_stats_free(&cons);
    tr2_stats_free(&slope);
    free(st_);
    free(ss_);
  }

  /* ---------------------------------------------------------------- 10 */
  printf("\n[10] DG1 на НЕполиномиальном решении: Ричардсон по вытекшей мощности\n");
  {
    double p[3];
    const int nn[3] = {16, 32, 64};
    for (int i = 0; i < 3; i++) {
      tr2_mesh um;
      if (tr2_mesh_uniform(&um, nn[i], nn[i], 1.0 / nn[i], 1.0 / nn[i]) != 0) return 1;
      double *st_ = calloc((size_t)um.ncell, sizeof *st_);
      double *ss_ = calloc((size_t)um.ncell, sizeof *ss_);
      if (!st_ || !ss_) return 1;
      for (int c = 0; c < um.ncell; c++) {
        st_[c] = 3.0;
        ss_[c] = 0.0; /* чистое поглощение: решение exp(-tau), не полином */
      }
      tr2_problem pr;
      memset(&pr, 0, sizeof pr);
      pr.mesh = &um;
      pr.dirs = &d;
      pr.phase = &ph;
      pr.sig_t = st_;
      pr.sig_s = ss_;
      pr.lim = TR2_LIM_OFF;
      pr.tol = TOL_ITER;
      pr.lb[0] = LB;
      tr2_stats a;
      if (tr2_solve(&pr, &a) != 0) return 1;
      p[i] = a.pout;
      printf("    %3dx%-3d  вытекшая мощность %.10f\n", nn[i], nn[i], p[i]);
      tr2_stats_free(&a);
      free(st_);
      free(ss_);
      tr2_mesh_free(&um);
    }
    double d1 = fabs(p[1] - p[0]), d2 = fabs(p[2] - p[1]);
    double order = d2 > 0.0 ? log2(d1 / d2) : 0.0;
    printf("    разности %.3e, %.3e ⇒ наблюдаемый порядок %.2f\n", d1, d2, order);
    printf("    (на ЭТОМ решении DG1 НЕ точен — в отличие от характеристического\n");
    printf("     ядра T1, где Бугер выходил на 1e-15; см. доклад)\n");
    check(d2 < d1, "DG1 converges under refinement");
    check(order > 2.5 && order < 3.5, "outflow functional superconverges at 2p+1 = 3");
  }

  for (int k = 0; k < 3; k++)
    tr2_mesh_free(&mc[k].m);
  tr2_dirs_free(&d);
  printf("\n%d/%d checks passed%s\n", g_total - g_fail, g_total, g_fail ? " — FAILED" : " — OK");
  return g_fail ? 1 : 0;
}
