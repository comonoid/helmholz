/*
 * fcheck — независимый Монте-Карло эталон §911-10 / §912-R3.
 *
 * Геометрия: две параллельные пластины 1×1, зазор h = 0.2.
 * Метод: площадьная выборка точки на нижней пластине (равномерно по
 * квадрату), косинус-веса (направление полусферы с плотностью
 * cosθ/π, выборка через равномерный диск), ТОЧНЫЙ тест попадания
 * луча в квадрат верхней пластины (без трассировки).
 *
 * Оценка: F_areaavg = (1/K) Σ hit  — несмещённая оценка среднего по
 * площади коэффициента видимости (косинусная плотность даёт вес
 * cosθ именно формы_viewfactor-интеграла).
 * B_true = π · F_areaavg (§911-10: «B_true = π·F = 2.167»).
 *
 * Ожидание §912-R3: F ≈ 0.690, B_true ≈ 2.167 ± 0.005 при K ≥ 1e6.
 * НК замкнутой формой: для X=Y=a/h=5 аналитический F двух
 * соосных квадратов печатается рядом (должен совпасть с MC в 4σ).
 *
 * Ключи: K=<число> (умолчание 1000000), h=<зазор> (умолчание 0.2).
 * RNG: xorshift64* с фиксированным зерном — прогон детерминирован.
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static uint64_t rng_state = 88172645463325252ULL; /* фиксированное зерно */

static uint64_t rng_u64(void) {
  uint64_t x = rng_state;
  x ^= x >> 12;
  x ^= x << 25;
  x ^= x >> 27;
  rng_state = x;
  return x * 2685821657736338717ULL;
}

static double rng_u01(void) {
  return (double)(rng_u64() >> 11) * (1.0 / 9007199254740992.0);
}

/* Аналитический F двух соосных параллельных квадратов a×a, зазор h
 * (X=Y=a/h; замкнутая форма Howell B-1). Для НК MC. */
static double F_closed_form(double a, double h) {
  double X = a / h, Y = a / h; /* Howell B-1, соосные прямоугольники */
  double X2 = X * X, Y2 = Y * Y;
  double ln_part = 0.5 * log((1.0 + X2) * (1.0 + Y2) / (1.0 + X2 + Y2)); /* ln√ = ½ln */
  double atan_part = X * sqrt(1.0 + Y2) * atan(X / sqrt(1.0 + Y2)) +
                     Y * sqrt(1.0 + X2) * atan(Y / sqrt(1.0 + X2)) - X * atan(X) - Y * atan(Y);
  return (2.0 / (M_PI * X * Y)) * (ln_part + atan_part);
}

int main(int argc, char **argv) {
  long K = 1000000;
  double h = 0.2;
  for (int i = 1; i < argc; i++) {
    if (strncmp(argv[i], "K=", 2) == 0)
      K = atol(argv[i] + 2);
    else if (strncmp(argv[i], "h=", 2) == 0)
      h = atof(argv[i] + 2);
    else {
      fprintf(stderr, "fcheck: не читается %s\n", argv[i]);
      return 1;
    }
  }
  if (K < 1 || h <= 0) {
    fprintf(stderr, "fcheck: плохие K/h\n");
    return 1;
  }

  long hits = 0;
  for (long k = 0; k < K; k++) {
    /* точка на нижней пластине, равномерно по площади */
    double u = rng_u01(), v = rng_u01();
    /* косинус-направление: равномерный диск (dx,dy), dz=sqrt(1-r²) */
    double r = sqrt(rng_u01()), phi = 2.0 * M_PI * rng_u01();
    double dx = r * cos(phi), dy = r * sin(phi), dz = sqrt(1.0 - dx * dx - dy * dy);
    /* пересечение с плоскостью z = h и точный тест квадрата [0,1]² */
    double t = h / dz;
    double x = u + t * dx, y = v + t * dy;
    if (x >= 0.0 && x <= 1.0 && y >= 0.0 && y <= 1.0) hits++;
  }
  double F = (double)hits / (double)K;
  double sigma_F = sqrt(F * (1.0 - F) / (double)K);
  double B = M_PI * F;
  double Fcf = F_closed_form(1.0, h);
  printf("fcheck: K=%ld h=%.4f hit=%ld\n", K, h, hits);
  printf("fcheck: F_areaavg=%.6f ±%.6f (1σ)\n", F, sigma_F);
  printf("fcheck: B_true=π·F=%.6f (2.167 ± 0.005 — ожидание §912-R3)\n", B);
  printf("fcheck: НК замкнутая форма (X=Y=%.2f): F=%.6f, MC-откл=%.2fσ\n", 1.0 / h, Fcf,
         fabs(F - Fcf) / sigma_F);
  return 0;
}
