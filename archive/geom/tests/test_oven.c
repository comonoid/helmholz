/* test_oven — Ш4: ОДИН ОТСКОК ПРОТИВ ЭТАЛОНА (PLAN_ELEMENTS.md).
 *
 * ПЕЧЬ. Замкнутая коробка, все грани с одинаковым `L_e` и одинаковым альбедо
 * `ρ`. В НЕПРЕРЫВНОМ переносе поле однородно и равно `L = L_e/(1−ρ)`, то есть
 * `E = π·L_e/(1−ρ)`. Вывод в строку: `E = ∫L cos dω = πL` по полусфере,
 * `L_out = L_e + ρE/π = L_e + ρL`, неподвижная точка `L = L_e/(1−ρ)`.
 *
 * ПЕЧЬ ЗДЕСЬ СВОЯ, А НЕ ИЗ `tests/test_sweep3.c`: та стоит на `tr3_mesh` —
 * ОБЪЁМНОМ носителе поля, отвергнутом §6. Эталон и оба инварианта (К13, К44)
 * перенесены дословно; перенести можно было только их, общего кода нет.
 *
 * ПРОВЕРКА РАЗДЕЛЕНА НА ТРИ, И ЭТО ГЛАВНОЕ РЕШЕНИЕ ЭТОГО ФАЙЛА. Первый прогон
 * дал печь мимо эталона на 1.4% при `ND = 32` и на 22…76% при `ND = 8`, а при
 * `ρ = 0.9` и `ND = 8` она РАСХОДИТСЯ (E растёт до 1e10). Смешанная проверка
 * сказала бы «перенос сломан». Он не сломан: вся невязка — УГЛОВАЯ, и её видно
 * порознь.
 *
 *   A. КВАДРАТУРА, БЕЗ ПЕРЕНОСА ВОВСЕ. `S(n) = Σ_{ω·n<0} w_ω|ω·n|` обязана
 *      равняться `π` для ЛЮБОЙ нормали. Это чистое свойство набора направлений,
 *      считается в три строки и ни одной строки переноса не трогает.
 *   B. ПЕРЕНОС ПРОТИВ КВАДРАТУРЫ. При `ρ = 0` и `L_e = 1` на всех гранях один
 *      сбор обязан дать `E_k = S(n_k)` ТОЧНО — с точностью растра и ничего
 *      больше. Это и есть проверка растеризации, редукции и решения, очищенная
 *      от угловой ошибки.
 *   C. ПЕЧЬ ЦЕЛИКОМ. Разрыв с `π·L_e/(1−ρ)` обязан совпасть с угловой ошибкой
 *      из A, усиленной рядом Неймана в `1/(1−ρ)` раз. Совпал — значит других
 *      источников ошибки нет.
 *
 * ОГРАНИЧИТЕЛЬ ВКЛЮЧЁН (К6: «ограничитель обязан сохранять константы»).
 * Проверяется не тем, что печь сошлась, а тем, что `nclip = 0`.
 *
 * ЧТО ПЕЧЬ НЕ ЛОВИТ (К44 дословно: «в однородном поле любая ошибка переноса
 * неотличима от правильного ответа») — три негативных контроля:
 *   D. РАЗНЫЕ радиансы стенок и их ПЕРЕСТАНОВКА: ответ обязан измениться на
 *      величину порядка САМОГО РАЗБРОСА по граням, а не средней яркости —
 *      среднее перестановка не трогает по построению, и метрика по нему слепа;
 *   E. сбить квадратурные веса на 1% ⇒ печь обязана поехать;
 *   F. разорвать коробку ⇒ `Φ_in < Φ_out`. Без этого «баланс сошёлся» не значит
 *      ничего: он сойдётся и на неполном наборе граней (узор К13/К40/К94).
 */

#include "polygon.h"
#include "psweep.h"
#include "transport/dirs3.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failed = 0;

/* Нормали граней коробки — ВНУТРЬ. */
static const double FN[6][3] = {{-1, 0, 0}, {1, 0, 0},  {0, -1, 0},
                                {0, 1, 0},  {0, 0, -1}, {0, 0, 1}};

static int build_oven(hz_polyset *ps, double a, int nfaces) {
  if (hz_poly_init_empty(ps) != 0) return 2;
  const double c[6][3] = {{a, 0, 0}, {-a, 0, 0}, {0, a, 0}, {0, -a, 0}, {0, 0, a}, {0, 0, -a}};
  const double e[6][3] = {{0, 1, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 1}, {1, 0, 0}, {1, 0, 0}};
  for (int i = 0; i < nfaces; i++)
    if (hz_poly_add_quad(ps, c[i], FN[i], e[i], a, a, 0) != 0) return 2;
  return 0;
}

/* S(n) = Σ_{ω·n<0} w_ω |ω·n|. Непрерывный ответ — π при любой n. */
static double quad_S(const tr3_dirs *d, const double n[3]) {
  double s = 0.0;
  for (int m = 0; m < d->n; m++) {
    double c = d->ox[m] * n[0] + d->oy[m] * n[1] + d->oz[m] * n[2];
    if (c < 0.0) s += d->w[m] * (-c);
  }
  return s;
}

static int run(hz_ptrans *t, const tr3_dirs *d, double h, double tol, hz_pstats *last, int maxb) {
  int nb = 0;
  for (; nb < maxb; nb++) {
    if (hz_psweep_bounce(t, d, h, 0, 0.0, HZ_LAYOUT_RUNS, last) != 0) return -1;
    if (last->dE < tol) break;
  }
  return nb + 1;
}

int main(void) {
  const double a = 0.5, Le = 1.0;
  const int NQ = 4;
  const int qmu[4] = {1, 2, 2, 4}, qphi[4] = {1, 2, 4, 4};
  printf("== test_oven: печь 1×1×1 м, нормали внутрь\n");

  /* --- A. КВАДРАТУРА САМА ПО СЕБЕ --------------------------------------- */
  printf("\n=== A. квадратура: S(n) = Σ w|ω·n| обязана быть π = %.6f при любой n ===\n", M_PI);
  printf("   %-6s %12s %12s %12s %12s\n", "ND", "S(ось z)", "S(ось x)", "S(диагон.)", "макс. ошиб");
  for (int q = 0; q < NQ; q++) {
    tr3_dirs d;
    if (tr3_dirs_product(&d, qmu[q], qphi[q]) != 0) return 1;
    double nz[3] = {0, 0, 1}, nx[3] = {1, 0, 0};
    double nd[3] = {0.57735026918962584, 0.57735026918962584, 0.57735026918962584};
    double sz = quad_S(&d, nz), sx = quad_S(&d, nx), sd = quad_S(&d, nd);
    double e = 0.0;
    for (int i = 0; i < 6; i++) {
      double ei = fabs(quad_S(&d, FN[i]) - M_PI) / M_PI;
      if (ei > e) e = ei;
    }
    printf("   %-6d %12.6f %12.6f %12.6f %12.4f%%\n", d.n, sz, sx, sd, 100.0 * e);
    tr3_dirs_free(&d);
  }
  printf("   Разбор: по μ квадратура ТОЧНА (∫μ dμ — многочлен первой степени),\n"
         "   по наклонной нормали — НЕТ: там под интегралом sqrt(1−μ²), у которого\n"
         "   производная на полюсе бесконечна, и Гаусс сходится алгебраически.\n"
         "   Отсюда и вся невязка печи ниже; переноса в ней нет ни грамма.\n");

  /* --- B. ПЕРЕНОС ПРОТИВ КВАДРАТУРЫ (ρ = 0, один сбор) ------------------- */
  printf("\n=== B. перенос против квадратуры: ρ = 0, E_k обязано дать S(n_k) ===\n");
  printf("   %-6s %-8s %14s %14s %12s\n", "ND", "h, м", "макс |E−S|/S", "макс |градиент|", "nclip");
  for (int q = 1; q < NQ; q++)
    for (int ih = 0; ih < 3; ih++) {
      double h = (double[]){0.05, 0.025, 0.0125}[ih];
      hz_polyset ps;
      if (build_oven(&ps, a, 6) != 0) return 1;
      hz_ptrans t;
      if (hz_ptrans_init(&t, &ps, NULL) != 0) return 1;
      for (int32_t k = 0; k < ps.np; k++) {
        t.rho[k] = 0.0;
        t.Le[k] = Le;
      }
      tr3_dirs d;
      if (tr3_dirs_product(&d, qmu[q], qphi[q]) != 0) return 1;
      hz_pstats st;
      if (hz_psweep_bounce(&t, &d, h, 0, 0.0, HZ_LAYOUT_RUNS, &st) != 0) return 1;
      double worst = 0.0, gmax = 0.0;
      for (int32_t k = 0; k < 6; k++) {
        double s = quad_S(&d, FN[k]) * Le;
        double e = fabs(t.E[(size_t)k * 3] - s) / s;
        if (e > worst) worst = e;
        double g = (fabs(t.E[(size_t)k * 3 + 1]) + fabs(t.E[(size_t)k * 3 + 2])) / s;
        if (g > gmax) gmax = g;
      }
      printf("   %-6d %-8g %14.3e %14.3e %12lld\n", d.n, h, worst, gmax, (long long)st.nclip);
      if (worst > 0.02) {
        printf("   FAIL: перенос разошёлся с собственной квадратурой больше чем на 2%%\n");
        failed = 1;
      }
      if (st.nclip != 0) {
        printf("   FAIL: ограничитель тронул однородное поле — К6 нарушен\n");
        failed = 1;
      }
      tr3_dirs_free(&d);
      hz_ptrans_free(&t);
      hz_poly_free(&ps);
    }

  /* --- C. ПЕЧЬ ЦЕЛИКОМ --------------------------------------------------- */
  printf("\n=== C. печь: E против π·L_e/(1−ρ), и сверка с усиленной угловой ошибкой ===\n");
  printf("   %-5s %-6s %6s %12s %12s %11s %13s\n", "ρ", "ND", "отск", "E ср", "эталон", "отн. ошиб",
         "ожид. ε/(1−ρ)");
  for (int ir = 0; ir < 3; ir++) {
    double rho = (double[]){0.3, 0.7, 0.9}[ir];
    for (int q = 0; q < NQ; q++) {
      hz_polyset ps;
      if (build_oven(&ps, a, 6) != 0) return 1;
      hz_ptrans t;
      if (hz_ptrans_init(&t, &ps, NULL) != 0) return 1;
      for (int32_t k = 0; k < ps.np; k++) {
        t.rho[k] = rho;
        t.Le[k] = Le;
      }
      tr3_dirs d;
      if (tr3_dirs_product(&d, qmu[q], qphi[q]) != 0) return 1;
      double eps = 0.0;
      for (int i = 0; i < 6; i++) {
        double ei = quad_S(&d, FN[i]) / M_PI - 1.0;
        if (fabs(ei) > fabs(eps)) eps = ei;
      }
      hz_pstats st;
      int nb = run(&t, &d, 0.0125, 1e-13, &st, 400);
      double ref = M_PI * Le / (1.0 - rho), sum = 0.0;
      for (int32_t k = 0; k < ps.np; k++)
        sum += t.E[(size_t)k * 3];
      double avg = sum / ps.np;
      int diverged = !(avg < 1e6);
      printf("   %-5.1f %-6d %6d %12.5f %12.5f %11.2e %13.2e%s\n", rho, d.n, nb, avg, ref,
             fabs(avg - ref) / ref, fabs(eps) / (1.0 - rho), diverged ? "  РАСХОДИТСЯ" : "");
      /* Расходимость при ρ·(1+ε) > 1 — не отказ кода, а СВОЙСТВО набора
       * направлений, и оно обязано быть видно, а не спрятано. */
      if (!diverged && rho * (1.0 + fabs(eps)) < 1.0 &&
          fabs(avg - ref) / ref > 3.0 * fabs(eps) / (1.0 - rho) + 0.01) {
        printf("   FAIL: невязка печи не объясняется угловой ошибкой — есть второй источник\n");
        failed = 1;
      }
      if (diverged && rho * (1.0 + fabs(eps)) < 1.0) {
        printf("   FAIL: печь разошлась там, где спектральный радиус меньше единицы\n");
        failed = 1;
      }
      tr3_dirs_free(&d);
      hz_ptrans_free(&t);
      hz_poly_free(&ps);
    }
  }

  /* --- F. БАЛАНС (К13) и разорванная коробка ------------------------------ */
  printf("\n=== F. баланс энергии и негативный контроль: разорванная коробка ===\n");
  for (int nf = 6; nf >= 5; nf--) {
    hz_polyset ps;
    if (build_oven(&ps, a, nf) != 0) return 1;
    hz_ptrans t;
    if (hz_ptrans_init(&t, &ps, NULL) != 0) return 1;
    for (int32_t k = 0; k < ps.np; k++) {
      t.rho[k] = 0.7;
      t.Le[k] = Le;
    }
    tr3_dirs d;
    if (tr3_dirs_product(&d, 4, 4) != 0) return 1;
    hz_pstats st;
    int nb = run(&t, &d, 0.0125, 1e-13, &st, 400);
    printf("   граней %d, отскоков %d: Φ_in = %.6f, Φ_out = %.6f, Φ_in/Φ_out = %.6f\n", nf, nb,
           st.phi_in, st.phi_out, st.phi_in / st.phi_out);
    if (nf == 6 && fabs(st.phi_in / st.phi_out - 1.0) > 0.02) {
      printf("   FAIL: замкнутая печь не сохраняет энергию\n");
      failed = 1;
    }
    if (nf == 5 && st.phi_in / st.phi_out > 0.95) {
      printf("   FAIL: разорванная коробка не потеряла энергию — баланс слеп\n");
      failed = 1;
    }
    tr3_dirs_free(&d);
    hz_ptrans_free(&t);
    hz_poly_free(&ps);
  }

  /* --- D. НЕГАТИВНЫЙ КОНТРОЛЬ (К44): РАЗНЫЕ радиансы стенок --------------- */
  printf("\n=== D. негативный контроль (К44): разные радиансы стенок и перестановка ===\n");
  {
    double Ea[6], Eb[6], gmax[2] = {0.0, 0.0};
    const double set1[6] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
    const double set2[6] = {6.0, 5.0, 4.0, 3.0, 2.0, 1.0};
    for (int pass = 0; pass < 2; pass++) {
      hz_polyset ps;
      if (build_oven(&ps, a, 6) != 0) return 1;
      hz_ptrans t;
      if (hz_ptrans_init(&t, &ps, NULL) != 0) return 1;
      for (int32_t k = 0; k < ps.np; k++) {
        t.rho[k] = 0.5;
        t.Le[k] = (pass == 0) ? set1[k] : set2[k];
      }
      tr3_dirs d;
      if (tr3_dirs_product(&d, 4, 4) != 0) return 1;
      hz_pstats st;
      run(&t, &d, 0.0125, 1e-13, &st, 400);
      for (int32_t k = 0; k < 6; k++) {
        double v = t.E[(size_t)k * 3];
        if (pass == 0)
          Ea[k] = v;
        else
          Eb[k] = v;
        double g = fabs(t.E[(size_t)k * 3 + 1]) + fabs(t.E[(size_t)k * 3 + 2]);
        if (g > gmax[pass]) gmax[pass] = g;
      }
      tr3_dirs_free(&d);
      hz_ptrans_free(&t);
      hz_poly_free(&ps);
    }
    printf("   набор 1: %.4f %.4f %.4f %.4f %.4f %.4f; max |градиент| %.4f\n", Ea[0], Ea[1], Ea[2],
           Ea[3], Ea[4], Ea[5], gmax[0]);
    printf("   набор 2: %.4f %.4f %.4f %.4f %.4f %.4f; max |градиент| %.4f\n", Eb[0], Eb[1], Eb[2],
           Eb[3], Eb[4], Eb[5], gmax[1]);
    /* МЕТРИКА ПО РАЗБРОСУ, А НЕ ПО СРЕДНЕМУ. Перестановка радиансов среднее не
     * трогает по построению (сумма та же), поэтому `Σ|ΔE|/Σ|E|` разбавлено
     * средним и в первом прогоне дало 0.086 — «не изменилось», хотя картина
     * ПЕРЕВЕРНУЛАСЬ. Сравнивать надо с тем, что перестановка и меняет. */
    double mn = Ea[0], mx = Ea[0], dmax = 0.0;
    for (int k = 0; k < 6; k++) {
      if (Ea[k] < mn) mn = Ea[k];
      if (Ea[k] > mx) mx = Ea[k];
      if (fabs(Ea[k] - Eb[k]) > dmax) dmax = fabs(Ea[k] - Eb[k]);
    }
    printf("   разброс по граням %.4f, макс |E1−E2| %.4f, отношение %.3f (обязано быть O(1))\n",
           mx - mn, dmax, dmax / (mx - mn));
    if (gmax[0] < 1e-6) {
      printf("   FAIL: при разных стенках поле осталось без градиента — линейный базис мёртв\n");
      failed = 1;
    }
    if (dmax / (mx - mn) < 0.5) {
      printf("   FAIL: перестановка радиансов не изменила ответ — перенос слеп\n");
      failed = 1;
    }
  }

  /* --- E. НЕГАТИВНЫЙ КОНТРОЛЬ: сбить квадратурные веса на 1% -------------- */
  printf("\n=== E. негативный контроль: веса квадратуры сбиты на 1%% ===\n");
  {
    double got[2];
    for (int pass = 0; pass < 2; pass++) {
      hz_polyset ps;
      if (build_oven(&ps, a, 6) != 0) return 1;
      hz_ptrans t;
      if (hz_ptrans_init(&t, &ps, NULL) != 0) return 1;
      for (int32_t k = 0; k < ps.np; k++) {
        t.rho[k] = 0.7;
        t.Le[k] = Le;
      }
      tr3_dirs d;
      if (tr3_dirs_product(&d, 4, 4) != 0) return 1;
      if (pass == 1)
        for (int i = 0; i < d.n; i++)
          d.w[i] *= 1.01;
      hz_pstats st;
      run(&t, &d, 0.0125, 1e-13, &st, 400);
      got[pass] = t.E[0];
      tr3_dirs_free(&d);
      hz_ptrans_free(&t);
      hz_poly_free(&ps);
    }
    printf(
        "   E(веса ×1) = %.6f, E(веса ×1.01) = %.6f, сдвиг %.4f%% (ожидание ≈1%%/(1−ρ) = 3.3%%)\n",
        got[0], got[1], 100.0 * (got[1] - got[0]) / got[0]);
    if (fabs(got[1] - got[0]) / got[0] < 0.005) {
      printf("   FAIL: печь не заметила сбоя весов — она измеряет не то, чем считает\n");
      failed = 1;
    }
  }

  printf("\n%s\n", failed ? "== test_oven: ЕСТЬ ОТКАЗЫ" : "== test_oven: OK");
  return failed;
}
