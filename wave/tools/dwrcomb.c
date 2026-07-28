/* DUAL-WEIGHTED RESIDUAL: combine a forward run and an adjoint run of tdg2d.
 *
 * WHY THIS EXISTS. On 07-27 the prediction "the residual points at the shadow
 * and beam boundaries" was measured and FAILED: 12.1% of |r|^2 fell within one
 * cell of the six geometric lines while those bands cover 4.9% of the area, a
 * concentration of only 2.5x. But a bare residual cannot point anywhere, because
 * it does not know what reaches the observer. The estimator that does is the
 * residual weighted by the ADJOINT solution (Becker & Rannacher, dual-weighted
 * residual): contribution to the picture ~ residual x importance.
 *
 * Both runs use the same mesh and the same matrix — only the right-hand side
 * differs (physical light vs a source on the camera aperture) — so row i means
 * the same condition in both files and the two may be multiplied entry by entry.
 *
 * THE NULL MUST BE WEIGHTED TOO. Comparing a weighted concentration against an
 * unweighted area fraction would be self-fulfilling: the importance field is
 * itself a beam, so anything weighted by it concentrates. The control here is
 * the importance-weighted area fraction.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {
  double ox, oy, vx, vy;
} ray;

static int on_band(const ray *L, int nL, double x, double y, double band) {
  for (int i = 0; i < nL; i++) {
    double t = (x - L[i].ox) * L[i].vx + (y - L[i].oy) * L[i].vy;
    if (t < 0.0) continue; /* the ray starts at the corner */
    double q = fabs(-(x - L[i].ox) * L[i].vy + (y - L[i].oy) * L[i].vx);
    if (q < band) return 1;
  }
  return 0;
}

int main(int argc, char **argv) {
  if (argc != 6) {
    printf("usage: dwrcomb <theta> <orad> <band> <Lbox> <ncell_side>\n"
           "  reads img/rows_fwd.bin and img/rows_adj.bin\n");
    return 1;
  }
  double th = atof(argv[1]), a = atof(argv[2]), band = atof(argv[3]);
  double L = atof(argv[4]);
  double dx = cos(th), dy = sin(th);

  /* the six lines: two edges of the shadow slab, two per reflected beam */
  ray Ln[6];
  double pbest = -1e300, pworst = 1e300, sxb = 0, syb = 0, sxw = 0, syw = 0;
  for (int i = 0; i < 4; i++) {
    double cx = (i & 1) ? a : -a, cy = (i & 2) ? a : -a;
    double pp = -cx * dy + cy * dx;
    if (pp > pbest) {
      pbest = pp;
      sxb = cx;
      syb = cy;
    }
    if (pp < pworst) {
      pworst = pp;
      sxw = cx;
      syw = cy;
    }
  }
  Ln[0] = (ray){sxb, syb, dx, dy};
  Ln[1] = (ray){sxw, syw, dx, dy};
  Ln[2] = (ray){-a, -a, -dx, dy};
  Ln[3] = (ray){-a, a, -dx, dy};
  Ln[4] = (ray){-a, -a, dx, -dy};
  Ln[5] = (ray){a, -a, dx, -dy};

  FILE *ff = fopen("img/rows_fwd.bin", "rb");
  FILE *fa = fopen("img/rows_adj.bin", "rb");
  if (!ff || !fa) {
    printf("dwrcomb: need both img/rows_fwd.bin and img/rows_adj.bin\n");
    if (ff) fclose(ff);
    if (fa) fclose(fa);
    return 1;
  }
  double hit_r = 0.0, tot_r = 0.0; /* bare residual */
  double hit_d = 0.0, tot_d = 0.0; /* residual x importance */
  double hit_w = 0.0, tot_w = 0.0; /* importance alone, for the weighted null */
  long n = 0, mismatch = 0;
  for (;;) {
    float rf[3], ra[3];
    if (fread(rf, sizeof(float), 3, ff) != 3) break;
    if (fread(ra, sizeof(float), 3, fa) != 3) break;
    if (fabs((double)rf[0] - (double)ra[0]) + fabs((double)rf[1] - (double)ra[1]) > 1e-6 * L)
      mismatch++;
    double x = rf[0], y = rf[1];
    double r2 = (double)rf[2] * (double)rf[2];
    double w2 = (double)ra[2] * (double)ra[2];
    int on = on_band(Ln, 6, x, y, band);
    tot_r += r2;
    tot_d += r2 * w2;
    tot_w += w2;
    if (on) {
      hit_r += r2;
      hit_d += r2 * w2;
      hit_w += w2;
    }
    n++;
  }
  fclose(ff);
  fclose(fa);
  if (mismatch) {
    printf("dwrcomb: ABORT — %ld rows disagree in position; the two runs used\n"
           "         different meshes and must not be combined\n",
           mismatch);
    return 1;
  }
  printf("  rows combined: %ld\n", n);
  printf("  bare residual   on the 6 lines: %5.1f%%\n", 100.0 * hit_r / tot_r);
  printf("  DUAL-WEIGHTED   on the 6 lines: %5.1f%%\n", 100.0 * hit_d / tot_d);
  printf("  weighted NULL (importance alone): %5.1f%%   <-- the control\n", 100.0 * hit_w / tot_w);
  printf("  concentration over the weighted null: %.2fx\n",
         (hit_d / tot_d) / (hit_w / tot_w > 0.0 ? hit_w / tot_w : 1.0));
  return 0;
}
