/* Time-domain probe: same double-slit scene as tools/render.c, but solved by
 * explicit second-order FDTD on a uniform grid, real arithmetic throughout.
 *   u_tt + sigma*u_t = c^2 lap u + f(t) delta_src
 * The complex amplitude is recovered ONLY at the end, by accumulating one
 * Fourier coefficient over the last period — phase is never tracked or stored.
 * Purpose: measure wall-clock against the frequency-domain solver (~130 s). */
#include "image.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* scene in cells (as in render.c): 96 x 96 x 16, lambda = 8 cells */
enum { CX = 96, CY = 96, CZ = 16, PPC = 2 }; /* grid points per cell */
enum { NX = CX * PPC, NY = CY * PPC, NZ = CZ * PPC };
static const double LAM_CELLS = 8.0;

static inline size_t IDX(int i, int j, int k) {
  return ((size_t)k * NY + (size_t)j) * NX + (size_t)i;
}

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

int main(int argc, char **argv) {
  int scene = 0; /* 0 = slit, 1 = smooth ball, 2 = rough ball */
  if (argc > 1 && argv[1][0] == 'b') scene = 1;
  if (argc > 1 && argv[1][0] == 'r') scene = 2;
  size_t n = (size_t)NX * NY * NZ;
  double h = 1.0 / (double)PPC;          /* grid step in cells */
  double c0 = 1.0;                       /* background speed, cells per time unit */
  double omega = 2.0 * M_PI / LAM_CELLS; /* k0*c0 */
  double dt = 0.99 * h / (c0 * sqrt(3.0));
  double period = 2.0 * M_PI / omega;
  int steps_per_period = (int)(period / dt + 0.5);
  int nper = 40; /* enough for the transient to leave through the damping shell */
  int nsteps = nper * steps_per_period;

  float *c2 = malloc(n * sizeof(float));
  float *sig = malloc(n * sizeof(float));
  float *u0 = calloc(n, sizeof(float));
  float *u1 = calloc(n, sizeof(float));
  float *u2 = calloc(n, sizeof(float));
  if (!c2 || !sig || !u0 || !u1 || !u2) return 1;

  /* --- scene: background, absorbing shell, wall with two slits ------------ */
  const int SHELL = 8 * PPC; /* cells -> points */
  for (int k = 0; k < NZ; k++)
    for (int j = 0; j < NY; j++)
      for (int i = 0; i < NX; i++) {
        double x = (double)i * h, y = (double)j * h, z = (double)k * h; /* cells */
        double cc = c0, ss = 0.0;
        if (scene == 0) {
          /* wall y in [46,50], slits x in [36,44] and [52,60] */
          int in_wall = (y >= 46.0 && y < 50.0);
          int in_slit = (x >= 36.0 && x < 44.0) || (x >= 52.0 && x < 60.0);
          if (in_wall && !in_slit) {
            cc = c0 / 2.5;
            ss = 2.0 * omega;
          }
        } else {
          /* ball at (48,56), R = 2 lambda; scene 2 adds lambda-scale roughness
           * (deterministic harmonic sum, feature size ~ lambda) */
          double dx = x - 48.0, dy = y - 56.0;
          double rr = sqrt(dx * dx + dy * dy);
          double R = 16.0;
          if (scene == 2 && rr > 1e-9) {
            double th = atan2(dy, dx);
            R += 2.6 * sin(9.0 * th + 0.7) + 2.2 * sin(12.0 * th + 2.1) +
                 2.4 * sin(15.0 * th + 4.3) + 1.8 * sin(19.0 * th + 5.5);
          }
          if (rr < R) cc = c0 / 2.5; /* dense, index 2.5 */
        }
        /* damping shell on x,y borders and z faces (quartic ramp) */
        double d = 1e9;
        d = fmin(d, (double)i);
        d = fmin(d, (double)(NX - 1 - i));
        d = fmin(d, (double)j);
        d = fmin(d, (double)(NY - 1 - j));
        d = fmin(d, (double)k);
        d = fmin(d, (double)(NZ - 1 - k));
        if (d < SHELL) {
          double t = (double)(SHELL - d) / (double)SHELL;
          ss += 1.5 * omega * t * t * t * t;
        }
        (void)z;
        size_t id = IDX(i, j, k);
        c2[id] = (float)(cc * cc);
        sig[id] = (float)ss;
      }

  /* source: smoothly ramped sinusoid at (48,24,8) cells */
  int si = (int)(48.0 / h), sj = (int)((scene == 0 ? 24.0 : 18.0) / h), sk = (int)(8.0 / h);
  size_t sidx = IDX(si, sj, sk);

  /* --- time stepping ------------------------------------------------------ */
  double inv_h2 = 1.0 / (h * h);
  double t0 = now_s();
  int fourier_from = nsteps - steps_per_period;
  double *re = calloc((size_t)NX * NY, sizeof(double));
  double *im = calloc((size_t)NX * NY, sizeof(double));
  if (!re || !im) return 1;
  int kslice = (int)(8.0 / h);

  for (int s = 0; s < nsteps; s++) {
    double t = (double)s * dt;
    double ramp = t < 4.0 * period ? 0.5 * (1.0 - cos(M_PI * t / (4.0 * period))) : 1.0;
    double src = ramp * sin(omega * t);
#pragma omp parallel for schedule(static)
    for (int k = 1; k < NZ - 1; k++)
      for (int j = 1; j < NY - 1; j++)
        for (int i = 1; i < NX - 1; i++) {
          size_t id = IDX(i, j, k);
          /* explicit widening: the field is stored float, the stencil runs in
           * double (gate: -Wdouble-promotion forbids the implicit mix) */
          double lap = ((double)u1[id + 1] + (double)u1[id - 1] + (double)u1[id + NX] +
                        (double)u1[id - NX] + (double)u1[id + (size_t)NX * NY] +
                        (double)u1[id - (size_t)NX * NY] - 6.0 * (double)u1[id]) *
                       inv_h2;
          double a = 0.5 * (double)sig[id] * dt;
          u2[id] = (float)((2.0 * (double)u1[id] - (1.0 - a) * (double)u0[id] +
                            dt * dt * (double)c2[id] * lap) /
                           (1.0 + a));
        }
    u2[sidx] += (float)(dt * dt * src * 200.0);
    float *tmp = u0;
    u0 = u1;
    u1 = u2;
    u2 = tmp;
    if (s >= fourier_from) {
      double ph = omega * t;
      double cw = cos(ph), sw = sin(ph);
#pragma omp parallel for schedule(static)
      for (int j = 0; j < NY; j++)
        for (int i = 0; i < NX; i++) {
          double v = u1[IDX(i, j, kslice)];
          re[(size_t)j * NX + (size_t)i] += v * cw;
          im[(size_t)j * NX + (size_t)i] -= v * sw;
        }
    }
  }
  double elapsed = now_s() - t0;
  double gupd = (double)n * (double)nsteps / 1e9;
  printf("FDTD: grid %dx%dx%d = %.2fM points, %d steps (%d per period, %d periods)\n", NX, NY, NZ,
         (double)n / 1e6, nsteps, steps_per_period, nper);
  printf("FDTD: %.1f s wall clock, %.2f Gupdates (%.0f Mupd/s)\n", elapsed, gupd,
         gupd * 1e3 / elapsed);

  /* --- images: same slice as the frequency-domain renderer ---------------- */
  double norm = 2.0 / (double)steps_per_period;
  double *ii = malloc((size_t)NX * NY * sizeof(double));
  double *rr = malloc((size_t)NX * NY * sizeof(double));
  if (!ii || !rr) return 1;
  for (int j = 0; j < NY; j++)
    for (int i = 0; i < NX; i++) {
      size_t id = (size_t)j * NX + (size_t)i;
      double a = re[id] * norm, b = im[id] * norm;
      ii[id] = a * a + b * b;
      double x = (double)i * h, y = (double)j * h;
      int inside;
      if (scene == 0) {
        int in_wall = (y >= 46.0 && y < 50.0);
        int in_slit = (x >= 36.0 && x < 44.0) || (x >= 52.0 && x < 60.0);
        inside = in_wall && !in_slit;
      } else {
        double dx = x - 48.0, dy = y - 56.0;
        double rq = sqrt(dx * dx + dy * dy), R = 16.0;
        if (scene == 2 && rq > 1e-9) {
          double th = atan2(dy, dx);
          R += 2.6 * sin(9.0 * th + 0.7) + 2.2 * sin(12.0 * th + 2.1) + 2.4 * sin(15.0 * th + 4.3) +
               1.8 * sin(19.0 * th + 5.5);
        }
        inside = rq < R;
      }
      rr[id] = inside ? 0.0 : a;
    }
  const char *nm = scene == 0 ? "slit" : (scene == 1 ? "ball" : "rough");
  char p1[128], p2[128];
  snprintf(p1, sizeof(p1), "result/fdtd_%s_int.ppm", nm);
  snprintf(p2, sizeof(p2), "result/fdtd_%s_re.ppm", nm);
  if (hz_ppm_write(p1, ii, NX, NY) == 0) printf("wrote %s\n", p1);
  if (hz_ppm_write_signed(p2, rr, NX, NY) == 0) printf("wrote %s\n", p2);

  free(c2);
  free(sig);
  free(u0);
  free(u1);
  free(u2);
  free(re);
  free(im);
  free(ii);
  free(rr);
  return 0;
}
