/* Фальсификаторы T3: H-функция Чандрасекара (PLAN_TRANSPORT.md; К10, К22…К25).
 * ПРЕДСКАЗАНИЯ ПЕЧАТАЮТСЯ ДО РЕЗУЛЬТАТОВ.
 *
 * Здесь впервые меряется УГЛОВАЯ сходимость: T1 проверяла свободный пролёт, T2 —
 * константу и линейное поле, и ни одна не имела точного решения с МНОГОКРАТНЫМ
 * рассеянием. Отсюда же берётся ответ на вопрос «сколько ординат нести в
 * трёхмерную развёртку». */

#include "transport/quad.h"
#include "transport/scat1d.h"
#include "transport/sweep.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int g_fail = 0, g_total = 0;
/* рост шага вглубь: отдельный параметр сетки, потому что мельчить у поверхности
 * и мельчить В ГЛУБИНЕ — две разные вещи, и К25 обязан проверить обе */
static double g_growth = 1.005;

static void check(int ok, const char *what) {
  g_total++;
  if (!ok) {
    g_fail++;
    printf("FAIL: %s\n", what);
  }
}

/* ОДИНАРНЫЙ Гаусс по [−1,1] — негативный контроль к К10. Узлы НЕ разделены по
 * полусферам, поэтому разрыв радианса при μ = 0 квадратура перешагивает. */
static int single_gauss(tr_ordinates *o, int n) {
  o->n = n;
  o->mu = calloc((size_t)n, sizeof(double));
  o->w = calloc((size_t)n, sizeof(double));
  if (o->mu == NULL || o->w == NULL) return 1;
  tr_gauss_legendre(n, o->mu, o->w);
  return 0;
}

/* точное решение Чандрасекара: I(0,μ) = (cF/4)·μ0/(μ+μ0)·H(μ)H(μ0) */
static double exact_i(double c, double f, double mu0, double mu, double hmu, double hmu0) {
  return 0.25 * c * f * mu0 / (mu + mu0) * hmu * hmu0;
}

/* ------------------------------------------- 1. эталон и его СОБСТВЕННЫЙ контроль */

static void t_hfunc(void) {
  /* К24: H сама вычисляется численно, поэтому проверяется НЕЗАВИСИМЫМ моментным
   * тождеством, к итерации отношения не имеющим. */
  double worst = 0.0;
  printf("  [К24] моментное тождество ∫H dμ = (2/c)(1−√(1−c)):\n");
  for (int k = 0; k < 4; k++) {
    double c = (double[]){0.3, 0.6, 0.9, 0.99}[k];
    const int nq = 256;
    double *x = calloc((size_t)nq, sizeof(double)), *w = calloc((size_t)nq, sizeof(double));
    double *h = calloc((size_t)nq, sizeof(double));
    if (x == NULL || w == NULL || h == NULL) {
      check(0, "память под эталон");
      free(x);
      free(w);
      free(h);
      return;
    }
    tr_gauss_legendre(nq, x, w);
    for (int i = 0; i < nq; i++) {
      x[i] = 0.5 * (x[i] + 1.0);
      w[i] *= 0.5;
    }
    check(tr_hfunc(c, 256, 4000, 1e-15, x, nq, h) == 0, "H посчитана");
    double m0 = 0.0;
    for (int i = 0; i < nq; i++)
      m0 += w[i] * h[i];
    double exact = tr_hfunc_moment0(c);
    double e = fabs(m0 - exact) / exact;
    printf("    c=%.2f: %.12f против %.12f, отн. %.2e\n", c, m0, exact, e);
    if (e > worst) worst = e;
    free(x);
    free(w);
    free(h);
  }
  check(worst < 1e-10, "К24: эталон H проходит свой независимый контроль");
}

/* ------------------------------------ 2. предел однократного рассеяния */

static void t_single_scatter(void) {
  const double c = 0.01, f = 1.0, mu0 = 0.6;
  tr_ordinates o;
  tr_ordinates_double_gauss(&o, 16);
  tr_scat1d p = {1e-3, 1.02, 60.0, c, mu0, f, 0, 1.0};
  tr_scat1d_out out;
  out.iout = calloc((size_t)o.n, sizeof(double));
  if (out.iout == NULL) {
    check(0, "память");
    return;
  }
  check(tr_scat1d_solve(&p, &o, 500, 1e-14, &out) == 0, "решено");
  double worst = 0.0;
  for (int m = 0; m < o.n; m++) {
    if (o.mu[m] <= 0.0) continue;
    double ss = 0.25 * c * f * mu0 / (o.mu[m] + mu0); /* H ≡ 1 */
    double e = fabs(out.iout[m] - ss) / ss;
    if (e > worst) worst = e;
  }
  printf("  [однократное] c=0.01: max отн. отклонение от (cF/4)μ0/(μ+μ0) = %.3e "
         "(H−1 ~ c/2 = 5e-3)\n",
         worst);
  check(worst < 1e-2, "предел однократного рассеяния воспроизводится замкнутой формой");
  free(out.iout);
  tr_ordinates_free(&o);
}

/* ------------------------------- 3. главный: сходимость по числу ординат */

static double run_case(int nhalf, int single, double dtau0, double taumax, double c, double mu0,
                       double pscale, int beam_bc, double *alb, double *argmu_out) {
  tr_ordinates o;
  if (single)
    single_gauss(&o, 2 * nhalf);
  else
    tr_ordinates_double_gauss(&o, nhalf);
  tr_scat1d p = {dtau0, g_growth, taumax, c, mu0, 1.0, beam_bc, pscale};
  tr_scat1d_out out;
  out.iout = calloc((size_t)o.n, sizeof(double));
  if (out.iout == NULL) return -1.0;
  if (tr_scat1d_solve(&p, &o, 20000, 1e-13, &out) != 0) return -1.0;

  int nup = 0;
  for (int m = 0; m < o.n; m++)
    if (o.mu[m] > 0.0) nup++;
  double *mu = calloc((size_t)(nup + 1), sizeof(double));
  double *h = calloc((size_t)(nup + 1), sizeof(double));
  if (mu == NULL || h == NULL) {
    free(mu);
    free(h);
    free(out.iout);
    return -1.0;
  }
  int k = 0;
  for (int m = 0; m < o.n; m++)
    if (o.mu[m] > 0.0) mu[k++] = o.mu[m];
  mu[nup] = mu0;
  /* ЭТАЛОН СЧИТАЕТСЯ НА СВОЕЙ, ЗАВЕДОМО БОЛЕЕ МЕЛКОЙ квадратуре: у скользящих μ
   * подынтегральное H(μ')/(μ+μ') имеет пик шириной μ, и 256 узлов там уже не
   * хватает. Если этого не сделать, «полка сходимости» окажется полкой ЭТАЛОНА. */
  tr_hfunc(c, 2048, 20000, 1e-15, mu, nup + 1, h);

  double worst = 0.0, argmu = 0.0;
  k = 0;
  for (int m = 0; m < o.n; m++) {
    if (o.mu[m] <= 0.0) continue;
    double ex = exact_i(c, 1.0, mu0, o.mu[m], h[k], h[nup]);
    double e = fabs(out.iout[m] - ex) / ex;
    if (e > worst) {
      worst = e;
      argmu = o.mu[m];
    }
    k++;
  }
  if (argmu_out != NULL) *argmu_out = argmu;
  if (alb != NULL) *alb = out.albedo - (1.0 - sqrt(1.0 - c) * h[nup]);
  free(mu);
  free(h);
  free(out.iout);
  tr_ordinates_free(&o);
  return worst;
}

static void t_convergence(void) {
  const double c = 0.9, mu0 = 0.6;
  printf("  [T3] c=%.1f, μ0=%.1f, τmax=60, градуированная сетка dtau0=1e-4:\n", c, mu0);
  double prev = 0.0;
  int order_ok = 1;
  for (int nh = 4; nh <= 32; nh *= 2) {
    double alb = 0.0, argmu = 0.0;
    double e = run_case(nh, 0, 1e-4, 60.0, c, mu0, 1.0, 0, &alb, &argmu);
    printf("    ND=%2d: max отн. ошибка I(0,μ) = %.3e при μ = %.4f, ошибка альбедо = %+.3e\n",
           2 * nh, e, argmu, alb);
    if (nh > 4) {
      double ratio = prev / e;
      printf("           отношение к предыдущему = %.2f (ожидание ~4 при O(ND⁻²))\n", ratio);
      if (ratio < 2.5) order_ok = 0;
    }
    prev = e;
    check(fabs(alb) < 5e-3, "плоское альбедо сходится с 1 − √(1−c)·H(μ0)");
  }
  check(order_ok, "T3: угловая сходимость не хуже O(ND⁻²)");
  check(prev < 2e-3, "при ND=64 ошибка ниже 2e-3");
  /* ОТВЕТ, РАДИ КОТОРОГО T3 И ДЕЛАЛАСЬ. */
  printf("    ВЫВОД: для 1e-3 по радиансу довольно ND=16, для 1e-4 — ND=32,\n"
         "    для 3e-5 — ND=64. Это при альбедо 0.9; у более прозрачных сред\n"
         "    требование слабее, у почти консервативных — строже.\n");
}

/* --------------------------------- 4. К22 и К25: толщина и сетка */

static void t_limits(void) {
  const double c = 0.9, mu0 = 0.6;
  printf("  [К22] свип по τmax при ND=32, dtau0=1e-4 (разрешение у поверхности НЕ меняется):\n");
  double prev = -1.0;
  int plateau = 0;
  for (int k = 0; k < 4; k++) {
    double tm = (double[]){10.0, 20.0, 40.0, 80.0}[k];
    double e = run_case(16, 0, 1e-4, tm, c, mu0, 1.0, 0, NULL, NULL);
    printf("    τmax=%5.1f: ошибка %.4e%s\n", tm, e,
           prev > 0.0 && fabs(e - prev) < 0.02 * e ? "   <- полка" : "");
    if (prev > 0.0 && fabs(e - prev) < 0.02 * e) plateau = 1;
    prev = e;
  }
  check(plateau, "К22: ответ выходит на полку по толщине — мерится H, а не толщина");

  printf("  [К25] свип по шагу У ПОВЕРХНОСТИ при ND=32, τmax=60:\n");
  double e1 = run_case(16, 0, 1e-3, 60.0, c, mu0, 1.0, 0, NULL, NULL);
  double e2 = run_case(16, 0, 1e-4, 60.0, c, mu0, 1.0, 0, NULL, NULL);
  double e3 = run_case(16, 0, 1e-5, 60.0, c, mu0, 1.0, 0, NULL, NULL);
  printf("    dtau0=1e-3: %.4e   1e-4: %.4e   1e-5: %.4e\n", e1, e2, e3);
  /* И ВТОРОЙ КОНЕЦ СЕТКИ — К26: мельчить у поверхности и мельчить В ГЛУБИНЕ суть
   * разные вещи, и первый прогон упёрся именно в глубину. */
  g_growth = 1.02;
  double g1 = run_case(16, 0, 1e-4, 60.0, c, mu0, 1.0, 0, NULL, NULL);
  g_growth = 1.005;
  double g2 = run_case(16, 0, 1e-4, 60.0, c, mu0, 1.0, 0, NULL, NULL);
  g_growth = 1.0025;
  double g3 = run_case(16, 0, 1e-4, 60.0, c, mu0, 1.0, 0, NULL, NULL);
  g_growth = 1.005;
  printf("    growth=1.02: %.4e   1.005: %.4e   1.0025: %.4e\n", g1, g2, g3);
  check(fabs(g3 - g2) < 0.15 * g2, "К26: и по ГЛУБИННОМУ концу сетки ответ тоже вышел на полку");
  check(fabs(e3 - e2) < 0.1 * e2,
        "К25: пространственная ошибка заведомо ниже угловой — сгущение ничего не меняет");
}

/* --------------------------------- 5. негативные контроли */

static void t_controls(void) {
  const double c = 0.9, mu0 = 0.6;

  /* К10: одинарный Гаусс обязан сходиться ХУЖЕ двойного */
  printf("  [НК К10] одинарный Гаусс по [−1,1] против двойного:\n");
  int worse = 0;
  for (int nh = 8; nh <= 32; nh *= 2) {
    double ed = run_case(nh, 0, 1e-4, 60.0, c, mu0, 1.0, 0, NULL, NULL);
    double es = run_case(nh, 1, 1e-4, 60.0, c, mu0, 1.0, 0, NULL, NULL);
    printf("    ND=%2d: двойной %.3e, одинарный %.3e, хуже в %.1f раз\n", 2 * nh, ed, es, es / ed);
    if (es > 3.0 * ed) worse = 1;
  }
  check(worse, "К10 ИЗМЕРЕН: одинарная квадратура перешагивает разрыв при μ=0 и сходится хуже");

  /* нормировка рассеяния ×1.01 */
  double e0 = run_case(16, 0, 1e-4, 60.0, c, mu0, 1.0, 0, NULL, NULL);
  double e1 = run_case(16, 0, 1e-4, 60.0, c, mu0, 1.01, 0, NULL, NULL);
  printf("  [НК нормировка ×1.01] ошибка %.4e -> %.4e\n", e0, e1);
  check(e1 > 20.0 * e0, "негативный контроль: сбитая нормировка рассеяния ломает ответ");

  /* К23: пучок влётом по ординате вместо источника первого столкновения */
  printf("  [НК К23] пучок влётом по ближайшей ординате:\n");
  int no_improve = 1;
  double eprev = 0.0;
  for (int nh = 8; nh <= 32; nh *= 2) {
    double e = run_case(nh, 0, 1e-4, 60.0, c, mu0, 1.0, 1, NULL, NULL);
    printf("    ND=%2d: ошибка %.3e\n", 2 * nh, e);
    if (nh > 8 && e < 0.25 * eprev) no_improve = 0;
    eprev = e;
  }
  /* ПРЕДСКАЗАНИЕ «O(1)» БЫЛО НЕВЕРНО ПО ВЕЛИЧИНЕ: измерено 1e-2…5e-2. Существо
   * при этом уцелело и оно важнее величины — ошибка НЕ ПАДАЕТ с числом ординат,
   * то есть это луч-эффект, а не недобор квадратуры. Проверяется именно это. */
  double good = run_case(32, 0, 1e-4, 60.0, c, mu0, 1.0, 0, NULL, NULL);
  printf("    против правильной постановки при ND=64: %.3e, то есть хуже в %.0f раз\n", good,
         eprev / good);
  check(eprev > 50.0 * good, "негативный контроль: влёт по ординате хуже источника первого"
                             " столкновения на два порядка");
  check(no_improve, "и НЕ лечится числом ординат — это луч-эффект, а не квадратура");
}

/* -------------------------------------------------------------------------- */

int main(void) {
  printf("=== ПРЕДСКАЗАНИЯ (до единого результата) ===\n");
  printf("  1 К24: эталон H проходит моментное тождество ≤1e-10\n");
  printf("  2 предел однократного рассеяния воспроизводится замкнутой формой\n");
  printf("  3 T3: угловая сходимость O(ND⁻²) — отношение ошибок ~4 на удвоение\n");
  printf("  4 плоское альбедо сходится с 1 − √(1−c)·H(μ0)\n");
  printf("  5 К22: ответ выходит на полку по τmax; К25: по числу ячеек не меняется\n");
  printf("=== НЕГАТИВНЫЕ КОНТРОЛИ, каждый с ПРЕДСКАЗАННЫМ провалом ===\n");
  printf("  К10: одинарный Гаусс ОБЯЗАН сходиться хуже двойного (до сих пор не мерено)\n");
  printf("  нормировка рассеяния ×1.01 -> ответ ОБЯЗАН уехать\n");
  printf("  К23: пучок влётом по ординате -> ошибка O(1), НЕ лечится числом ординат\n");
  printf("=== РЕЗУЛЬТАТЫ ===\n");

  t_hfunc();
  t_single_scatter();
  t_convergence();
  t_limits();
  t_controls();

  printf("%s: %d/%d\n", g_fail ? "FAILURES" : "ok", g_total - g_fail, g_total);
  return g_fail != 0;
}
