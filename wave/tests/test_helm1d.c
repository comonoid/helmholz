/* M2 experiments E1/E2: cascade vs independent FD reference (PLAN.md). */
#include "helm1d.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int g_fail = 0;
static int g_total = 0;

static void check_lt(double got, double bound, const char *what) {
  g_total++;
  if (!(got < bound)) {
    g_fail++;
    printf("FAIL: %s: got %.6g, bound %.6g\n", what, got, bound);
  }
}

enum { NCELL = 256, PER_CELL = 16 };
static const double SENS0 = 192.0, SENS1 = 208.0;
static const double XSRC = 48.0;

static double sensor_err(const hz_sol1d *sol, const double complex *ufd) {
  double num = 0.0, den = 0.0;
  for (int i = 0; i <= NCELL * PER_CELL; i++) {
    double x = (double)i / (double)PER_CELL;
    if (x < SENS0 || x > SENS1) continue;
    double complex d = hz_sol1d_eval(sol, x) - ufd[i];
    num += creal(d * conj(d));
    den += creal(ufd[i] * conj(ufd[i]));
  }
  return sqrt(num / den);
}

static void print_stats(const char *tag, const hz_lvlstat1d *st, const hz_level1d *lv, int n) {
  for (int l = 0; l < n; l++)
    printf("  %s lvl h=%-3d dim=%-4d cond_eff=%-8.3g null=%-2d res_drop=%.3g\n", tag,
           1 << lv[l].lvl, st[l].dim, st[l].cond_eff, st[l].nnull, st[l].res_drop);
}

static double run_case(const char *tag, const hz_med1d *med, const hz_level1d *lv, int nlev) {
  double complex *ufd = calloc((size_t)NCELL * PER_CELL + 1, sizeof(double complex));
  hz_lvlstat1d *st = calloc((size_t)nlev, sizeof(hz_lvlstat1d));
  if (ufd == NULL || st == NULL) {
    free(ufd);
    free(st);
    g_total++;
    g_fail++;
    printf("FAIL: %s: alloc\n", tag);
    return 1e9;
  }
  g_total++;
  if (hz_fd_reference(med, PER_CELL, XSRC, ufd) != 0) {
    g_fail++;
    printf("FAIL: %s: fd reference\n", tag);
  }
  hz_sol1d sol;
  g_total++;
  if (hz_cascade1d_solve(med, lv, nlev, XSRC, 1, &sol, st) != 0) {
    g_fail++;
    printf("FAIL: %s: cascade solve\n", tag);
    free(ufd);
    free(st);
    return 1e9;
  }
  double err = sensor_err(&sol, ufd);
  printf("%s: sensor rel L2 err = %.4f\n", tag, err);
  /* phase probe: ratio u_repr/u_fd at three sensor points (|r|, arg r) */
  for (int p = 0; p < 3; p++) {
    double x = SENS0 + (SENS1 - SENS0) * 0.25 * (double)(p + 1);
    int i = (int)(x * (double)PER_CELL);
    double complex r = hz_sol1d_eval(&sol, x) / ufd[i];
    printf("  %s probe x=%.0f ratio |%.4f| arg %+.4f rad\n", tag, x, cabs(r), carg(r));
  }
  print_stats(tag, st, lv, nlev);
  hz_sol1d_free(&sol);
  free(ufd);
  free(st);
  return err;
}

/* E3 runner: no stats, one line per variant */
static double run_e3(const char *tag, const hz_med1d *med, const hz_level1d *lv, int nlev,
                     int nsweeps) {
  double complex *ufd = calloc((size_t)NCELL * PER_CELL + 1, sizeof(double complex));
  if (ufd == NULL || hz_fd_reference(med, PER_CELL, XSRC, ufd) != 0) {
    free(ufd);
    g_total++;
    g_fail++;
    printf("FAIL: %s: setup\n", tag);
    return 1e9;
  }
  hz_sol1d sol;
  g_total++;
  if (hz_cascade1d_solve(med, lv, nlev, XSRC, nsweeps, &sol, NULL) != 0) {
    g_fail++;
    printf("FAIL: %s: solve\n", tag);
    free(ufd);
    return 1e9;
  }
  double err = sensor_err(&sol, ufd);
  printf("%s: sensor err = %.4f\n", tag, err);
  hz_sol1d_free(&sol);
  free(ufd);
  return err;
}

int main(void) {
  double k0 = 2.0 * M_PI / 16.0; /* lambda = 16 cells */
  static const hz_level1d LV_FLOOR[5] = {
      {5, -1, 9}, {4, -1, 17}, {3, -1, 33}, {2, -1, 65}, {1, -1, 129}};
  static const hz_level1d LV_FINE[6] = {{5, -1, 9},  {4, -1, 17},  {3, -1, 33},
                                        {2, -1, 65}, {1, -1, 129}, {0, -1, 257}};

  static const hz_level1d LV_ONLY1[1] = {{0, -1, 257}};
  static const hz_level1d LV_ONLY2[1] = {{1, -1, 129}};

  /* E1: empty medium, transport source -> sensor across ~9 lambda */
  hz_med1d med;
  hz_med1d_init(&med, NCELL, k0, 0.02, 32);
  /* E0: pure single-level Galerkin, no cascade — separates cascade error from
   * discretization error */
  run_case("E0 empty only h=2", &med, LV_ONLY2, 1);
  run_case("E0 empty only h=1", &med, LV_ONLY1, 1);
  double e1 = run_case("E1 empty floor(h=2)", &med, LV_FLOOR, 5);
  double e1f = run_case("E1 empty fine(h=1)", &med, LV_FINE, 6);

  /* E2: slab in the middle, n = 1.5 */
  hz_med1d_slab(&med, 96, 128, k0, 1.5, 0.02);
  double e2 = run_case("E2 slab   floor(h=2)", &med, LV_FLOOR, 5);
  double e2f = run_case("E2 slab   fine(h=1)", &med, LV_FINE, 6);
  hz_med1d_free(&med);

  /* E3 far detail: thin barrier (cells 80..82) ~7 lambda from the sensor.
   * A = floor only (Born via exact integrals), B = + h=1 at the SENSOR only
   * (old intuition — must NOT help), C = + h=1 around the barrier only
   * (local refinement — must match D), D = + h=1 everywhere (reference set). */
  static const hz_level1d LV_B[6] = {{5, -1, 9},  {4, -1, 17},  {3, -1, 33},
                                     {2, -1, 65}, {1, -1, 129}, {0, 186, 214}};
  static const hz_level1d LV_C[6] = {{5, -1, 9},  {4, -1, 17},  {3, -1, 33},
                                     {2, -1, 65}, {1, -1, 129}, {0, 72, 90}};
  /* wide patch: +lambda/2 (8 cells) each side to hold the evanescent skirt */
  static const hz_level1d LV_CW[6] = {{5, -1, 9},  {4, -1, 17},  {3, -1, 33},
                                      {2, -1, 65}, {1, -1, 129}, {0, 64, 98}};
  static const double CONTRAST[2] = {1.2, 3.0};
  for (int ci = 0; ci < 2; ci++) {
    hz_med1d mb;
    hz_med1d_init(&mb, NCELL, k0, 0.02, 32);
    hz_med1d_slab(&mb, 80, 82, k0, CONTRAST[ci], 0.02);
    printf("E3 barrier contrast %.1f:\n", CONTRAST[ci]);
    double ea = run_e3("  A floor only          ", &mb, LV_FLOOR, 5, 1);
    run_e3("  B fine at sensor      ", &mb, LV_B, 6, 1);
    double ec = run_e3("  C fine at barrier     ", &mb, LV_C, 6, 1);
    double ec2 = run_e3("  C2 barrier, 2 sweeps  ", &mb, LV_C, 6, 2);
    run_e3("  C3 wide patch, 2 sweeps", &mb, LV_CW, 6, 2);
    run_e3("  C4 wide patch, 3 sweeps", &mb, LV_CW, 6, 3);
    double ed = run_e3("  D fine everywhere     ", &mb, LV_FINE, 6, 1);
    double ed2 = run_e3("  D2 everywhere, 2 sweeps", &mb, LV_FINE, 6, 2);
    hz_med1d_free(&mb);
    /* measured finding 07-25: one-way cascade cannot feed a local refinement
     * back into the transport level — C equals A at the sensor exactly */
    check_lt(fabs(ec - ea), 0.01, "E3: one-sweep local refinement invisible (one-way cascade)");
    check_lt(ec2, ed2 * 4.0 + 0.02, "E3: 2-sweep local refinement must approach full-fine");
    check_lt(ed2, ed + 1e-12, "E3: second sweep must not hurt full-fine");
  }

  /* bounds = measured 07-25 (e1=0.0207, e1f=0.0168, e2=0.0474, e2f=0.0151)
   * plus ~50% headroom; growth here = transport/assembly regression */
  check_lt(e1, 0.032, "E1 floor err");
  check_lt(e1f, e1 + 1e-12, "E1: adding h=1 must not hurt");
  check_lt(e2, 0.072, "E2 floor err");
  check_lt(e2f, 0.024, "E2 fine err (sub-lambda refinement must pay)");
  check_lt(e2f, e2 + 1e-12, "E2: adding h=1 must not hurt");

  printf("test_helm1d: %d/%d passed\n", g_total - g_fail, g_total);
  return g_fail > 0 ? 1 : 0;
}
