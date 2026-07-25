/* renderer driver: build a named scene, solve, dump slices + camera images.
 * usage: render <slit|box|lens|cornell> [outdir]  (default outdir result)
 * lambda = 8 cells, floor h=1 (lambda/8); per scene two slice images:
 *   *_int.ppm  — intensity |u|^2 (inferno-ish)
 *   *_re.ppm   — Re u (blue-white-red wavefronts, medium blanked to white) */
#include "camera.h"
#include "image.h"
#include "solver3d.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double ALPHA = 0.04;

static void shell_xy(hz_octree *t, const int dom[3], double k0, int th) {
  double b = k0 * k0;
  for (int s = 0; s < th; s++) {
    double tt = (double)(th - s) / (double)th;
    double a = ALPHA + (1.5 - ALPHA) * tt * tt * tt * tt;
    int lo[3] = {s, s, 0}, hi[3] = {dom[0] - s, dom[1] - s, dom[2]};
    hz_oct_set_box(t, lo, hi, CMPLX(b, b * a));
  }
  int lo[3] = {th, th, 0}, hi[3] = {dom[0] - th, dom[1] - th, dom[2]};
  hz_oct_set_box(t, lo, hi, CMPLX(b, b * ALPHA));
}

static void shell_z(hz_octree *t, const int dom[3], double k0, int th) {
  double b = k0 * k0;
  for (int s = 0; s < th; s++) {
    double tt = (double)(th - s) / (double)th;
    double a = ALPHA + (1.5 - ALPHA) * tt * tt * tt * tt;
    int lo0[3] = {0, 0, s}, hi0[3] = {dom[0], dom[1], s + 1};
    int lo1[3] = {0, 0, dom[2] - s - 1}, hi1[3] = {dom[0], dom[1], dom[2] - s};
    hz_oct_set_box(t, lo0, hi0, CMPLX(b, b * a));
    hz_oct_set_box(t, lo1, hi1, CMPLX(b, b * a));
  }
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s <slit|box|lens|cornell> [outdir]\n", argv[0]);
    return 2;
  }
  const char *scene = argv[1];
  const char *outdir = argc > 2 ? argv[2] : "result";
  double k0 = 2.0 * M_PI / 8.0; /* lambda = 8 cells */
  double b = k0 * k0;
  double complex bg = CMPLX(b, b * ALPHA);

  hz_octree t;
  if (hz_oct_init(&t, 7, bg) != 0) return 1;
  int dom[3] = {96, 96, 16};
  hz_src3 src[1];
  int top = 3;
  double slice_at = 8.0;
  hz_camera cam;
  int use_cam = 0;

  if (strcmp(scene, "slit") == 0) {
    shell_xy(&t, dom, k0, 8);
    shell_z(&t, dom, k0, 3);
    double kw = 2.5 * k0;
    int w0[3] = {8, 46, 0}, w1[3] = {88, 50, 16};
    hz_oct_set_box(&t, w0, w1, CMPLX(kw * kw, kw * kw * 0.8)); /* absorbing wall */
    int s10[3] = {36, 46, 3}, s11[3] = {44, 50, 13};
    int s20[3] = {52, 46, 3}, s21[3] = {60, 50, 13};
    hz_oct_set_box(&t, s10, s11, bg); /* two slits, width lambda, pitch 2 lambda */
    hz_oct_set_box(&t, s20, s21, bg);
    src[0].p = (hz_pot3){0, {48, 24, 8}};
    src[0].amp = CMPLX(1.0, 0.0);
  } else if (strcmp(scene, "box") == 0) {
    shell_xy(&t, dom, k0, 8);
    shell_z(&t, dom, k0, 3);
    double kw = 2.5 * k0;
    int b0[3] = {40, 44, 0}, b1[3] = {56, 54, 16};
    hz_oct_set_box(&t, b0, b1, CMPLX(kw * kw, kw * kw * 0.8));
    src[0].p = (hz_pot3){0, {48, 20, 8}};
    src[0].amp = CMPLX(1.0, 0.0);
  } else if (strcmp(scene, "lens") == 0) {
    shell_xy(&t, dom, k0, 8);
    shell_z(&t, dom, k0, 3);
    double kl = 1.5 * k0;
    double c[3] = {48.0, 48.0, 8.0};
    hz_oct_set_ball(&t, c, 14.0, CMPLX(kl * kl, kl * kl * ALPHA)); /* r = 3.5 lambda */
    src[0].p = (hz_pot3){0, {48, 16, 8}};
    src[0].amp = CMPLX(1.0, 0.0);
  } else if (strcmp(scene, "cornell") == 0) {
    dom[0] = dom[1] = dom[2] = 64;
    slice_at = 32.0;
    double kw = 2.0 * k0;
    double complex wall = CMPLX(kw * kw, kw * kw * 0.5);
    int lo[5][3] = {{0, 0, 0}, {0, 60, 0}, {0, 0, 0}, {60, 0, 0}, {0, 0, 60}};
    int hi[5][3] = {{64, 4, 64}, {64, 64, 64}, {4, 64, 64}, {64, 64, 64}, {64, 64, 64}};
    for (int wl = 0; wl < 5; wl++)
      hz_oct_set_box(&t, lo[wl], hi[wl], wall);
    double kb = 1.6 * k0;
    double c[3] = {32.0, 40.0, 36.0};
    hz_oct_set_ball(&t, c, 10.0, CMPLX(kb * kb, kb * kb * 0.05)); /* r = 2.5 lambda */
    src[0].p = (hz_pot3){0, {32, 10, 32}};
    src[0].amp = CMPLX(1.0, 0.0);
    use_cam = 1;
    cam.cx = 32.0;
    cam.cy = 32.0;
    cam.zap = 8.0;
    cam.D = 48.0;
    cam.n = 128;
    double dobj = 36.0 - cam.zap;
    cam.di = 24.0;
    cam.f = 1.0 / (1.0 / dobj + 1.0 / cam.di);
    cam.k0 = k0;
  } else {
    fprintf(stderr, "unknown scene %s\n", scene);
    return 2;
  }

  hz_scene3 sc = {&t, {dom[0], dom[1], dom[2]}, bg, src, 1};
  hz_sol3 sol;
  hz_solver3d_set_tol(1e-3); /* picture-grade tolerance: discretization error dominates */
  printf("scene %s: dom %dx%dx%d, lambda=8, floor h=1\n", scene, dom[0], dom[1], dom[2]);
  if (hz_solve3d(&sc, top, 0, 2, 1, &sol) != 0) {
    fprintf(stderr, "solve failed\n");
    return 1;
  }

  char path[512];
  int W = dom[0] * 4, H = dom[1] * 4;
  double *ii = malloc((size_t)W * (size_t)H * sizeof(double));
  double *re = malloc((size_t)W * (size_t)H * sizeof(double));
  if (ii == NULL || re == NULL) return 1;
  for (int j = 0; j < H; j++)
    for (int i = 0; i < W; i++) {
      double p[3] = {((double)i + 0.5) / 4.0, ((double)j + 0.5) / 4.0, slice_at};
      double complex u = hz_sol3_eval(&sol, p);
      size_t idx = (size_t)j * (size_t)W + (size_t)i;
      ii[idx] = creal(u * conj(u));
      /* blank the medium (walls/objects differ in Re k2; the shell does not) */
      double complex k2 = hz_oct_at(&t, (int)(p[0]), (int)(p[1]), (int)slice_at);
      re[idx] = (fabs(creal(k2) - b) > 0.2 * b) ? 0.0 : creal(u);
    }
  snprintf(path, sizeof(path), "%s/%s_int.ppm", outdir, scene);
  if (hz_ppm_write(path, ii, W, H) == 0) printf("wrote %s\n", path);
  snprintf(path, sizeof(path), "%s/%s_re.ppm", outdir, scene);
  if (hz_ppm_write_signed(path, re, W, H) == 0) printf("wrote %s\n", path);
  free(ii);
  free(re);

  if (use_cam) {
    double *ci = malloc((size_t)cam.n * (size_t)cam.n * sizeof(double));
    if (ci != NULL && hz_camera_shoot(&sol, &cam, ci) == 0) {
      snprintf(path, sizeof(path), "%s/%s_camera.ppm", outdir, scene);
      if (hz_ppm_write(path, ci, cam.n, cam.n) == 0) printf("wrote %s\n", path);
    }
    free(ci);
  }

  hz_sol3_free(&sol);
  hz_oct_free(&t);
  return 0;
}
