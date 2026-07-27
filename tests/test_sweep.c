/* test_sweep.c — приёмка вехи T1 (PLAN_TRANSPORT.md): ядро развёртки,
 * чистое поглощение, эталон — закон Бугера.
 *
 * ПРЕДСКАЗАНИЯ ЗАПИСАНЫ ДО ПРОГОНА и печатаются перед результатами. Условия
 * вехи: слой, σ_s = 0, ε = 0 (аудит К2 — с объёмным излучением exp(−τ)
 * перестаёт быть эталоном), σ_t линейная по ячейке, τ трапецией ТОЧНО.
 *
 * Две метрики, и у КАЖДОЙ свой негативный контроль — веса квадратуры и
 * замкнутая форма τ входят в РАЗНЫЕ числа, поэтому один контроль на оба не
 * годится:
 *   [Бугер]  относительная ошибка радианса   ← контроль Б (левый конец вместо
 *                                               трапеции ⇒ первый порядок)
 *   [Поток]  ошибка углового момента         ← контроль А (веса ×1.01 ⇒ 1%)
 * То, что контроль А НЕ шевелит метрику Бугера, записано заранее: она весов не
 * видит вовсе.
 */
#include "transport/sweep.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int g_fail = 0;
static int g_total = 0;

static void check(int ok, const char *what) {
  g_total++;
  if (!ok) {
    g_fail++;
    printf("  FAIL: %s\n", what);
  }
}

/* e^{−25} ≈ 1e-11: глубже относительная ошибка радианса определяется не схемой,
 * а округлением в самой оптической толщине (δL/L = δτ ≈ τ·n·ε), и метрика
 * перестаёт мерить то, ради чего заведена. Скользящие ординаты с τ > 25
 * проверяются отдельной метрикой — относительной ошибкой τ. */
#define TAU_FAIR 25.0

/* Пол округления характеристического ядра. L = Π exp(−τ_i) по ячейкам: на
 * ячейку ~1 ulp от exp и ~0.5 ulp от умножения, итого ~1.5·n·ε, что при n = 64
 * даёт 2e-14. Поэтому порог метрики Бугера — 5e-14, а не 1e-14 из плана;
 * РЕШАЮЩИЙ признак не порог, а ТРЕНД: округление растёт ∝ n, квадратура вместо
 * замкнутой формы падала бы степенью Δz. */
#define TOL_BOUGUER 5e-14
#define TOL_EXACT 1e-14 /* тождества квадратуры и моментов: чистая арифметика */

/* Аналитическая оптическая толщина для u_L(z) = a + b·z на пути от z0 до z1
 * в направлении с косинусом μ. Считается ЗАМКНУТОЙ ФОРМОЙ ∫(a+bz)dz, а не тем
 * же трапецеидальным циклом, что в ядре, — иначе сравнение было бы тавтологией. */
static double tau_linear(double a, double b, double z0, double z1, double amu) {
  return (a * (z1 - z0) + 0.5 * b * (z1 * z1 - z0 * z0)) / amu;
}

/* --- метрика Бугера на заданном разбиении -------------------------------- */
/* Возвращает max по ординатам с τ ≤ TAU_FAIR от |L_num/L_ref − 1|;
 * в *worst_tau — max по ВСЕМ ординатам от |τ_num/τ_ref − 1|;
 * в *min_dz — САМАЯ МЕЛКАЯ ячейка разбиения. Последнее не украшение: без него
 * нельзя отличить «ошибка не падает, потому что схема точна» от «ошибка не
 * падает, потому что сетка перестала измельчаться». */
static double bouguer_error(const tr_ordinates *o, int ncell, double stretch, double zmax, double a,
                            double b, double *worst_tau, double *min_dz) {
  tr_slab s;
  if (tr_slab_alloc(&s, ncell) != 0) {
    printf("  FAIL: slab alloc\n");
    g_fail++;
    return 1.0;
  }
  tr_slab_grade(&s, zmax, stretch);
  tr_slab_set_linear(&s, a, b);

  double mdz = zmax;
  for (int i = 0; i < ncell; i++) {
    double dz = s.z[i + 1] - s.z[i];
    if (dz < mdz) mdz = dz;
  }
  *min_dz = mdz;

  const double l0 = 1.0;
  double worst = 0.0, worst_t = 0.0;
  for (int m = 0; m < o->n; m++) {
    double amu = fabs(o->mu[m]);
    double tau_ref = tau_linear(a, b, 0.0, zmax, amu);
    double l_ref = l0 * exp(-tau_ref);
    double l_num = tr_sweep_absorb(&s, o->mu[m], 1.0, l0, NULL);
    double tau_num = -log(l_num / l0);
    double et = fabs(tau_num / tau_ref - 1.0);
    if (et > worst_t) worst_t = et;
    if (tau_ref <= TAU_FAIR) {
      double e = fabs(l_num / l_ref - 1.0);
      if (e > worst) worst = e;
    }
  }
  tr_slab_free(&s);
  *worst_tau = worst_t;
  return worst;
}

int main(void) {
  printf("=== T1: ядро развёртки, чистое поглощение ===\n\n");
  printf("ПРЕДСКАЗАНИЯ (записаны до прогона):\n");
  printf("  1. тождества квадратуры Σw=2, Σ_{μ>0}w=1, Σ_{μ>0}wμ=1/2      <= %.0e\n", TOL_EXACT);
  printf("  2. Бугер, линейная σ_t, градуированная сетка 64 ячейки        <= %.0e\n", TOL_BOUGUER);
  printf("     (цель плана 1e-14; пол округления ~1.5·n·eps = 2e-14)\n");
  printf("  3. независимость от разбиения: ncell 4..1024 — ошибка НЕ падает\n");
  printf("     степенью Dz (округление, ~n*eps). Если ~1e-8 и падает — где-то\n");
  printf("     квадратура вместо замкнутой формы. СНАЧАЛА проверяется, что\n");
  printf("     min_dz падает ~1/ncell: иначе «не изменилось» ничего не значит.\n");
  printf("  4. профиль на ВСЕХ гранях, оба знака mu                       <= %.0e\n", TOL_BOUGUER);
  printf("  5. РАЗРЫВНАЯ sigma_t (скачок на грани)                        <= %.0e\n", TOL_BOUGUER);
  printf("  6. глубокий слой tau=20, ошибка нормирована на ЭТАЛОН e^-20   <= %.0e\n", TOL_BOUGUER);
  printf("  7. направленное sigma_t = |mu|*u_L: tau не зависит от mu,\n");
  printf("     поток F+ = pi*L0*exp(-tau) ТОЧНО                           <= %.0e\n", TOL_EXACT);
  printf("  8. подпись артефакта 1: пропускание ОБЯЗАНО меняться при\n");
  printf("     изменении b, d и mu — строгая монотонность\n");
  printf("  9. негативный контроль А (веса x1.01): ошибка потока ОБЯЗАНА\n");
  printf("     стать 1.000e-2; метрика Бугера при этом НЕ меняется\n");
  printf(" 10. негативный контроль Б (левый конец sigma вместо трапеции):\n");
  printf("     ошибка >= 1e-3 и ПЕРВЫЙ порядок — при удвоении ncell /2.0\n\n");

  /* ---------------------------------------------------------------- 1 */
  enum { NHALF = 8 };
  tr_ordinates o;
  if (tr_ordinates_double_gauss(&o, NHALF) != 0) {
    printf("FATAL: ordinates alloc\n");
    return 1;
  }
  printf("[1] двойная гауссова квадратура, %d ординат\n", o.n);
  {
    double sw = 0.0, swp = 0.0, swmu = 0.0, swmu2 = 0.0;
    for (int m = 0; m < o.n; m++) {
      sw += o.w[m];
      if (o.mu[m] > 0.0) {
        swp += o.w[m];
        swmu += o.w[m] * o.mu[m];
        swmu2 += o.w[m] * o.mu[m] * o.mu[m];
      }
    }
    printf("    Sw-2 = %.2e   Sw+ -1 = %.2e   Swmu-1/2 = %.2e   Swmu2-1/3 = %.2e\n", fabs(sw - 2.0),
           fabs(swp - 1.0), fabs(swmu - 0.5), fabs(swmu2 - 1.0 / 3.0));
    check(fabs(sw - 2.0) < TOL_EXACT, "Sw = 2");
    check(fabs(swp - 1.0) < TOL_EXACT, "Sw over upper half = 1");
    check(fabs(swmu - 0.5) < TOL_EXACT, "Swmu = 1/2");
    check(fabs(swmu2 - 1.0 / 3.0) < TOL_EXACT, "Swmu^2 = 1/3");
    /* узлы строго внутри полуинтервалов ⇒ mu = 0 в наборе нет */
    int interior = 1, symmetric = 1;
    for (int m = 0; m < o.n; m++) {
      if (fabs(o.mu[m]) < 1e-12 || fabs(o.mu[m]) > 1.0) interior = 0;
      if (fabs(o.mu[m] + o.mu[o.n - 1 - m]) > TOL_EXACT) symmetric = 0;
    }
    check(interior, "0 < |mu| <= 1 for every ordinate");
    check(symmetric, "ordinate set symmetric about mu = 0");
    printf("    |mu| in [%.4f, %.4f] — mu=0 в набор не попадает\n", fabs(o.mu[NHALF]),
           fabs(o.mu[o.n - 1]));
  }

  /* ---------------------------------------------------------------- 2, 3 */
  const double A = 0.3, B = 1.4, ZMAX = 2.0, STRETCH = 8.0;
  printf("\n[2] Бугер: u_L = %.1f + %.1f z, слой [0,%.1f], градуировка (крайние ячейки) %.0f\n", A,
         B, ZMAX, STRETCH);
  {
    double wt = 0.0, mdz = 0.0;
    double e = bouguer_error(&o, 64, STRETCH, ZMAX, A, B, &wt, &mdz);
    printf("    ncell=64: max|L/L_ref-1| (tau<=%.0f) = %.3e   max|tau/tau_ref-1| = %.3e\n",
           TAU_FAIR, e, wt);
    check(e < TOL_BOUGUER, "Bouguer error at 64 cells");
    check(wt < TOL_BOUGUER, "optical-depth error at 64 cells (all ordinates)");
  }

  printf("\n[3] независимость от разбиения (градуировка %.0f, тот же слой)\n", STRETCH);
  {
    const int nn[] = {4, 16, 64, 256, 1024};
    double err[5], mdz[5];
    for (int i = 0; i < 5; i++) {
      double wt = 0.0;
      err[i] = bouguer_error(&o, nn[i], STRETCH, ZMAX, A, B, &wt, &mdz[i]);
      printf("    ncell=%5d  min_dz=%.3e  err=%.3e  err_tau=%.3e\n", nn[i], mdz[i], err[i], wt);
      check(err[i] < TOL_BOUGUER, "Bouguer error stays at rounding floor");
    }
    /* СНАЧАЛА убедиться, что параметр действует: мельчайшая ячейка обязана
     * падать ∝ 1/ncell. Без этой проверки «ошибка не изменилась» ничего не
     * значит — сетка могла просто перестать измельчаться (первая подпись). */
    double refine = mdz[0] / mdz[4];
    printf("    min_dz(4)/min_dz(1024) = %.1f  — параметр действует, ожидалось ~256\n", refine);
    check(refine > 200.0 && refine < 320.0, "refinement actually refines (min cell ~ 1/ncell)");
    /* И только теперь: 256-кратное сгущение НЕ уменьшает ошибку степенью Δz.
     * Квадратура первого порядка дала бы падение примерно в 256 раз. */
    double drop = err[0] > 0.0 ? err[4] / err[0] : 1.0;
    printf("    err(1024)/err(4) = %.2f  — при квадратуре 1-го порядка было бы ~1/256\n", drop);
    check(drop > 1.0 / 16.0, "error does NOT fall as a power of Dz (closed form, not quadrature)");
  }

  /* ---------------------------------------------------------------- 4 */
  printf("\n[4] профиль на гранях, оба направления обхода\n");
  {
    enum { NC = 48 };
    tr_slab s;
    if (tr_slab_alloc(&s, NC) != 0) return 1;
    tr_slab_grade(&s, ZMAX, STRETCH);
    tr_slab_set_linear(&s, A, B);
    double *lf = calloc((size_t)NC + 1, sizeof *lf);
    if (!lf) {
      tr_slab_free(&s);
      return 1;
    }
    double worst_up = 0.0, worst_dn = 0.0;
    for (int m = 0; m < o.n; m++) {
      double amu = fabs(o.mu[m]);
      if (tau_linear(A, B, 0.0, ZMAX, amu) > TAU_FAIR) continue;
      /* вверх: вток на z=0 */
      tr_sweep_absorb(&s, amu, 1.0, 1.0, lf);
      for (int i = 0; i <= NC; i++) {
        double ref = exp(-tau_linear(A, B, 0.0, s.z[i], amu));
        double e = fabs(lf[i] / ref - 1.0);
        if (e > worst_up) worst_up = e;
      }
      /* вниз: вток на z=zmax */
      tr_sweep_absorb(&s, -amu, 1.0, 1.0, lf);
      for (int i = 0; i <= NC; i++) {
        double ref = exp(-tau_linear(A, B, s.z[i], ZMAX, amu));
        double e = fabs(lf[i] / ref - 1.0);
        if (e > worst_dn) worst_dn = e;
      }
    }
    printf("    max по граням: вверх %.3e   вниз %.3e\n", worst_up, worst_dn);
    check(worst_up < TOL_BOUGUER, "face profile, sweep upward");
    check(worst_dn < TOL_BOUGUER, "face profile, sweep downward");
    free(lf);
    tr_slab_free(&s);
  }

  /* ---------------------------------------------------------------- 5 */
  printf("\n[5] РАЗРЫВНАЯ sigma_t: две линейные ветви со скачком на z=1\n");
  {
    enum { NC = 64 };
    const double A1 = 0.2, B1 = 0.9, A2 = 2.5, B2 = -0.4, ZCUT = 1.0;
    tr_slab s;
    if (tr_slab_alloc(&s, NC) != 0) return 1;
    tr_slab_grade(&s, ZMAX, 1.0); /* равномерная ⇒ z=1 ровно на грани NC/2 */
    for (int i = 0; i < NC; i++) {
      double zl = s.z[i], zh = s.z[i + 1];
      s.sig_lo[i] = zl < ZCUT ? A1 + B1 * zl : A2 + B2 * zl;
      s.sig_hi[i] = zh <= ZCUT ? A1 + B1 * zh : A2 + B2 * zh;
    }
    printf("    скачок на z=1: %.2f -> %.2f\n", A1 + B1 * ZCUT, A2 + B2 * ZCUT);
    double worst = 0.0;
    for (int m = 0; m < o.n; m++) {
      double amu = fabs(o.mu[m]);
      double tref = tau_linear(A1, B1, 0.0, ZCUT, amu) + tau_linear(A2, B2, ZCUT, ZMAX, amu);
      if (tref > TAU_FAIR) continue;
      double lu = tr_sweep_absorb(&s, amu, 1.0, 1.0, NULL);
      double ld = tr_sweep_absorb(&s, -amu, 1.0, 1.0, NULL);
      double ref = exp(-tref);
      double eu = fabs(lu / ref - 1.0), ed = fabs(ld / ref - 1.0);
      if (eu > worst) worst = eu;
      if (ed > worst) worst = ed;
    }
    printf("    max|L/L_ref-1| (оба направления) = %.3e\n", worst);
    check(worst < TOL_BOUGUER, "discontinuous sigma_t, both sweep directions");
    tr_slab_free(&s);
  }

  /* ---------------------------------------------------------------- 6 */
  printf("\n[6] глубокий слой: tau = 20 (подпись 3 — нормируем на эталон)\n");
  {
    enum { NC = 100 };
    const double SIG = 10.0;
    tr_slab s;
    if (tr_slab_alloc(&s, NC) != 0) return 1;
    tr_slab_grade(&s, ZMAX, 1.0);
    tr_slab_set_linear(&s, SIG, 0.0);
    double amu = o.mu[o.n - 1];
    double tref = SIG * ZMAX / amu;
    double ref = exp(-tref);
    double l = tr_sweep_absorb(&s, amu, 1.0, 1.0, NULL);
    printf("    tau=%.4f  L=%.6e  L_ref=%.6e  |L/L_ref-1| = %.3e\n", tref, l, ref,
           fabs(l / ref - 1.0));
    printf("    (нормировка на L0=1 дала бы %.1e — вот это и было бы подписью 3)\n", fabs(l - ref));
    check(fabs(l / ref - 1.0) < TOL_BOUGUER, "deep slab, error relative to the reference itself");
    tr_slab_free(&s);
  }

  /* ---------------------------------------------------------------- 7 */
  printf("\n[7] направленное sigma_t = G(mu)*u_L, G = |mu|\n");
  {
    enum { NC = 64 };
    tr_slab s;
    if (tr_slab_alloc(&s, NC) != 0) return 1;
    tr_slab_grade(&s, ZMAX, STRETCH);
    tr_slab_set_linear(&s, A, B);
    /* G = |mu| сокращает 1/|mu| в длине хорды ⇒ tau = ∫u_L dz, одинаковая
     * для всех направлений. Это НАМЕРЕННАЯ независимость от параметра, поэтому
     * ниже она подпирается проверкой, что от A значение всё-таки зависит. */
    double tref = A * ZMAX + 0.5 * B * ZMAX * ZMAX;
    double ref = exp(-tref);
    double *l = calloc((size_t)o.n, sizeof *l);
    if (!l) {
      tr_slab_free(&s);
      return 1;
    }
    double spread = 0.0;
    for (int m = 0; m < o.n; m++) {
      double v = tr_sweep_absorb(&s, o.mu[m], fabs(o.mu[m]), 1.0, NULL);
      double e = fabs(v / ref - 1.0);
      if (e > spread) spread = e;
      l[m] = o.mu[m] > 0.0 ? v : 0.0; /* поток вверх: только уходящая половина */
    }
    printf("    tau=%.4f одинаково для всех ординат, max|L/L_ref-1| = %.3e\n", tref, spread);
    check(spread < TOL_BOUGUER, "directional sigma_t: tau direction-independent, exact");

    /* Момент: F+ = 2pi * S w mu L = pi*L0*exp(-tau) ТОЧНО — двойной Гаусс
     * интегрирует mu точно, и вот здесь нормировка весов входит в число. */
    double f = tr_moment1(&o, l);
    double f_ref = M_PI * ref;
    double phi = tr_moment0(&o, l);
    double phi_ref = 2.0 * M_PI * ref;
    printf("    F+ = %.12e  ref = %.12e  err = %.3e\n", f, f_ref, fabs(f / f_ref - 1.0));
    printf("    phi+= %.12e  ref = %.12e  err = %.3e\n", phi, phi_ref, fabs(phi / phi_ref - 1.0));
    check(fabs(f / f_ref - 1.0) < TOL_EXACT, "normal flux F+ = pi L0 exp(-tau)");
    check(fabs(phi / phi_ref - 1.0) < TOL_EXACT, "half-range scalar flux = 2pi L0 exp(-tau)");

    /* --- негативный контроль А: веса x1.01 --- */
    const double BUMP = 1.01;
    for (int m = 0; m < o.n; m++)
      o.w[m] *= BUMP;
    double f_bad = tr_moment1(&o, l);
    double err_bad = fabs(f_bad / f_ref - 1.0);
    double wt = 0.0, mdz = 0.0;
    double bou_bad = bouguer_error(&o, 64, STRETCH, ZMAX, A, B, &wt, &mdz);
    for (int m = 0; m < o.n; m++)
      o.w[m] /= BUMP;
    printf("\n[7b] НЕГАТИВНЫЙ КОНТРОЛЬ А: веса x%.2f\n", BUMP);
    printf("    ошибка потока: %.3e -> %.6e   (обязана быть 1.000000e-02)\n", fabs(f / f_ref - 1.0),
           err_bad);
    printf("    метрика Бугера при этом: %.3e  (весов не видит — так и предсказано)\n", bou_bad);
    check(fabs(err_bad - (BUMP - 1.0)) < 1e-12, "weights x1.01 => flux error exactly 1%");
    check(bou_bad < TOL_BOUGUER, "weights x1.01 leaves the Bouguer metric untouched (predicted)");
    /* веса восстановлены — тождество должно вернуться */
    check(fabs(tr_moment1(&o, l) / f_ref - 1.0) < TOL_EXACT, "weights restored");
    free(l);
    tr_slab_free(&s);
  }

  /* ---------------------------------------------------------------- 8 */
  printf("\n[8] подпись артефакта 1: число обязано меняться от своих параметров\n");
  {
    enum { NC = 32 };
    tr_slab s;
    if (tr_slab_alloc(&s, NC) != 0) return 1;
    double amu = o.mu[o.n - 1];
    int mono_b = 1, mono_mu = 1, mono_d = 1;
    double prev = 2.0;
    printf("    по градиенту b:");
    for (int i = 0; i < 4; i++) {
      double b = 0.5 * (double)i;
      tr_slab_grade(&s, ZMAX, STRETCH);
      tr_slab_set_linear(&s, A, b);
      double v = tr_sweep_absorb(&s, amu, 1.0, 1.0, NULL);
      printf("  b=%.1f L=%.6f", b, v);
      if (v >= prev) mono_b = 0;
      prev = v;
    }
    printf("\n    по толщине d:");
    prev = 2.0;
    for (int i = 1; i <= 4; i++) {
      double d = 0.5 * (double)i;
      tr_slab_grade(&s, d, STRETCH);
      tr_slab_set_linear(&s, A, B);
      double v = tr_sweep_absorb(&s, amu, 1.0, 1.0, NULL);
      printf("  d=%.1f L=%.6f", d, v);
      if (v >= prev) mono_d = 0;
      prev = v;
    }
    printf("\n    по направлению mu:");
    tr_slab_grade(&s, ZMAX, STRETCH);
    tr_slab_set_linear(&s, A, B);
    prev = -1.0;
    for (int m = NHALF; m < o.n; m++) { /* mu по возрастанию */
      double v = tr_sweep_absorb(&s, o.mu[m], 1.0, 1.0, NULL);
      if (m >= o.n - 3) printf("  mu=%.3f L=%.3e", o.mu[m], v);
      if (v <= prev) mono_mu = 0;
      prev = v;
    }
    printf("\n");
    check(mono_b, "transmittance strictly falls with the sigma_t gradient b");
    check(mono_d, "transmittance strictly falls with slab thickness");
    check(mono_mu, "transmittance strictly rises with mu (shorter chord)");
    /* и направленный случай: mu-независимость намеренная, но от A зависеть обязан */
    tr_slab_set_linear(&s, A, B);
    double v1 = tr_sweep_absorb(&s, amu, fabs(amu), 1.0, NULL);
    tr_slab_set_linear(&s, 2.0 * A, B);
    double v2 = tr_sweep_absorb(&s, amu, fabs(amu), 1.0, NULL);
    printf("    направленный случай: A -> 2A даёт L %.6f -> %.6f\n", v1, v2);
    check(v2 < v1 * (1.0 - 1e-3), "directional case still responds to sigma_t magnitude");
    tr_slab_free(&s);
  }

  /* ---------------------------------------------------------------- 9 */
  printf("\n[9] НЕГАТИВНЫЙ КОНТРОЛЬ Б: левый конец sigma вместо трапеции\n");
  {
    const int nn[] = {64, 128, 256};
    double err[3];
    double amu = o.mu[o.n - 1];
    for (int i = 0; i < 3; i++) {
      tr_slab s;
      if (tr_slab_alloc(&s, nn[i]) != 0) return 1;
      tr_slab_grade(&s, ZMAX, 1.0); /* равномерная — чтобы порядок читался чисто */
      tr_slab_set_linear(&s, A, B);
      /* тот же примитив, но σ_выход подменена на σ_вход ⇒ квадратура 1-го
       * порядка. Никакого флага в рабочем коде для этого не нужно. */
      double l = 1.0;
      for (int c = 0; c < nn[i]; c++)
        l = tr_cell_absorb(s.sig_lo[c], s.sig_lo[c], (s.z[c + 1] - s.z[c]) / amu, l);
      double ref = exp(-tau_linear(A, B, 0.0, ZMAX, amu));
      err[i] = fabs(l / ref - 1.0);
      printf("    ncell=%4d  err=%.4e\n", nn[i], err[i]);
      tr_slab_free(&s);
    }
    double r1 = err[0] / err[1], r2 = err[1] / err[2];
    printf("    отношения: %.3f, %.3f  (первый порядок ⇒ 2.0)\n", r1, r2);
    check(err[0] > 1e-3, "left-endpoint control lands far above the rounding floor");
    check(fabs(r1 - 2.0) < 0.1 && fabs(r2 - 2.0) < 0.1, "left-endpoint control is first order");
  }

  tr_ordinates_free(&o);
  printf("\n%d/%d checks passed%s\n", g_total - g_fail, g_total, g_fail ? " — FAILED" : " — OK");
  return g_fail ? 1 : 0;
}
