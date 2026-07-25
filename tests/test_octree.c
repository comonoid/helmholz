/* M3 acceptance: queries and visits vs brute-force voxel array, roundtrip. */
#include "octree.h"
#include <complex.h>
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
    printf("FAIL: %s\n", what);
  }
}

enum { L2 = 5, S = 1 << L2 }; /* 32^3 */

static double complex ref[S][S][S];

static void ref_box(const int lo[3], const int hi[3], double complex k2) {
  for (int x = lo[0] < 0 ? 0 : lo[0]; x < (hi[0] > S ? S : hi[0]); x++)
    for (int y = lo[1] < 0 ? 0 : lo[1]; y < (hi[1] > S ? S : hi[1]); y++)
      for (int z = lo[2] < 0 ? 0 : lo[2]; z < (hi[2] > S ? S : hi[2]); z++)
        ref[x][y][z] = k2;
}

static void ref_ball(const double c[3], double r, double complex k2) {
  for (int x = 0; x < S; x++)
    for (int y = 0; y < S; y++)
      for (int z = 0; z < S; z++) {
        double dx = (double)x + 0.5 - c[0], dy = (double)y + 0.5 - c[1],
               dz = (double)z + 0.5 - c[2];
        if (dx * dx + dy * dy + dz * dz < r * r) ref[x][y][z] = k2;
      }
}

typedef struct {
  double complex sum;
  long cells;
} visit_acc;

static int acc_cb(void *ctx, const int blo[3], const int bhi[3], double complex k2) {
  visit_acc *a = ctx;
  long vol = (long)(bhi[0] - blo[0]) * (bhi[1] - blo[1]) * (bhi[2] - blo[2]);
  a->sum += k2 * (double)vol;
  a->cells += vol;
  return 0;
}

int main(void) {
  double complex bg = CMPLX(1.0, 0.02);
  hz_octree t;
  check(hz_oct_init(&t, L2, bg) == 0, "init");
  for (int x = 0; x < S; x++)
    for (int y = 0; y < S; y++)
      for (int z = 0; z < S; z++)
        ref[x][y][z] = bg;

  /* a wall, a box, a ball, an overwrite */
  int w0[3] = {0, 14, 0}, w1[3] = {S, 18, S};
  check(hz_oct_set_box(&t, w0, w1, CMPLX(4.0, 0.1)) == 0, "set wall");
  ref_box(w0, w1, CMPLX(4.0, 0.1));
  int b0[3] = {3, 3, 3}, b1[3] = {9, 7, 5};
  check(hz_oct_set_box(&t, b0, b1, CMPLX(2.25, 0.0)) == 0, "set box");
  ref_box(b0, b1, CMPLX(2.25, 0.0));
  double c[3] = {22.0, 24.0, 12.0};
  check(hz_oct_set_ball(&t, c, 6.5, CMPLX(9.0, 0.5)) == 0, "set ball");
  ref_ball(c, 6.5, CMPLX(9.0, 0.5));
  int o0[3] = {20, 22, 10}, o1[3] = {24, 26, 14};
  check(hz_oct_set_box(&t, o0, o1, bg) == 0, "overwrite");
  ref_box(o0, o1, bg);

  /* point queries vs brute force */
  int bad = 0;
  for (int x = 0; x < S; x++)
    for (int y = 0; y < S; y++)
      for (int z = 0; z < S; z++)
        if (cabs(hz_oct_at(&t, x, y, z) - ref[x][y][z]) > 1e-15) bad++;
  check(bad == 0, "all point queries match brute force");

  /* visits over random-ish boxes: cell count and k2 volume-sum must match */
  static const int probes[5][6] = {{0, 0, 0, S, S, S},
                                   {2, 13, 1, 11, 19, 6},
                                   {15, 15, 5, 32, 32, 20},
                                   {21, 23, 11, 23, 25, 13},
                                   {7, 0, 0, 8, 32, 32}};
  for (int p = 0; p < 5; p++) {
    const int *pr = probes[p];
    int lo[3] = {pr[0], pr[1], pr[2]}, hi[3] = {pr[3], pr[4], pr[5]};
    visit_acc a = {0.0, 0};
    check(hz_oct_visit(&t, lo, hi, acc_cb, &a) == 0, "visit runs");
    long cells = 0;
    double complex sum = 0.0;
    for (int x = lo[0]; x < hi[0]; x++)
      for (int y = lo[1]; y < hi[1]; y++)
        for (int z = lo[2]; z < hi[2]; z++) {
          cells++;
          sum += ref[x][y][z];
        }
    check(a.cells == cells, "visit covers exactly the box");
    check(cabs(a.sum - sum) < 1e-9 * (cabs(sum) + 1.0), "visit k2 volume-sum matches");
  }

  /* serialization roundtrip */
  FILE *f = tmpfile();
  check(f != NULL, "tmpfile");
  if (f != NULL) {
    check(hz_oct_save(&t, f) == 0, "save");
    check(fseek(f, 0, SEEK_SET) == 0, "seek");
    hz_octree t2;
    check(hz_oct_load(&t2, f) == 0, "load");
    fclose(f);
    bad = 0;
    for (int x = 0; x < S; x += 3)
      for (int y = 0; y < S; y += 3)
        for (int z = 0; z < S; z += 3)
          if (cabs(hz_oct_at(&t2, x, y, z) - hz_oct_at(&t, x, y, z)) > 0.0) bad++;
    check(bad == 0, "roundtrip queries identical");
    check(t2.n == t.n, "roundtrip node count");
    hz_oct_free(&t2);
  }
  hz_oct_free(&t);

  printf("test_octree: %d/%d passed\n", g_total - g_fail, g_total);
  return g_fail > 0 ? 1 : 0;
}
