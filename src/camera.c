#include "camera.h"
#include "fft.h"
#include "phi.h"
#include <math.h>
#include <stdlib.h>

static int next_pow2(int v) {
  int p = 1;
  while (p < v)
    p *= 2;
  return p;
}

int hz_camera_shoot(const hz_sol3 *sol, const hz_camera *cam, double *intensity) {
  int n = cam->n;
  double step = cam->D / (double)n;
  int P = next_pow2(2 * n); /* padding against wrap-around */
  double complex *fld = calloc((size_t)P * (size_t)P, sizeof(double complex));
  if (fld == NULL) return 1;

  /* sample aperture field and apply the lens mask */
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      double x = cam->cx - 0.5 * cam->D + ((double)i + 0.5) * step;
      double y = cam->cy - 0.5 * cam->D + ((double)j + 0.5) * step;
      double p[3] = {x, y, cam->zap};
      double complex u = hz_sol3_eval(sol, p);
      if (cam->f > 0.0) {
        double dx = x - cam->cx, dy = y - cam->cy;
        double ph = -cam->k0 * (dx * dx + dy * dy) / (2.0 * cam->f);
        u *= CMPLX(cos(ph), sin(ph));
      }
      fld[(size_t)j * (size_t)P + (size_t)i] = u;
    }

  /* angular spectrum: U(k) *= exp(i kz di), kz = sqrt(k0^2 - kx^2 - ky^2) */
  double complex *line = calloc((size_t)P, sizeof(double complex));
  if (line == NULL) {
    free(fld);
    return 1;
  }
  for (int j = 0; j < P; j++)
    hz_fft(fld + (size_t)j * (size_t)P, P, -1);
  for (int i = 0; i < P; i++) {
    for (int j = 0; j < P; j++)
      line[j] = fld[(size_t)j * (size_t)P + (size_t)i];
    hz_fft(line, P, -1);
    for (int j = 0; j < P; j++)
      fld[(size_t)j * (size_t)P + (size_t)i] = line[j];
  }
  double dk = 2.0 * M_PI / ((double)P * step);
  for (int j = 0; j < P; j++)
    for (int i = 0; i < P; i++) {
      int qi = i <= P / 2 ? i : i - P;
      int qj = j <= P / 2 ? j : j - P;
      double kx = (double)qi * dk, ky = (double)qj * dk;
      double complex kz = csqrt(CMPLX(cam->k0 * cam->k0 - kx * kx - ky * ky, 0.0));
      if (cimag(kz) < 0.0) kz = -kz; /* decaying evanescent branch */
      fld[(size_t)j * (size_t)P + (size_t)i] *= cexp(CMPLX(0.0, 1.0) * kz * cam->di);
    }
  for (int j = 0; j < P; j++)
    hz_fft(fld + (size_t)j * (size_t)P, P, 1);
  for (int i = 0; i < P; i++) {
    for (int j = 0; j < P; j++)
      line[j] = fld[(size_t)j * (size_t)P + (size_t)i];
    hz_fft(line, P, 1);
    for (int j = 0; j < P; j++)
      fld[(size_t)j * (size_t)P + (size_t)i] = line[j];
  }
  for (int j = 0; j < n; j++)
    for (int i = 0; i < n; i++) {
      double complex u = fld[(size_t)j * (size_t)P + (size_t)i];
      intensity[(size_t)j * (size_t)n + (size_t)i] = creal(u * conj(u));
    }
  free(line);
  free(fld);
  return 0;
}

int hz_slice_z(const hz_sol3 *sol, double z, int w, int h, double *intensity) {
  for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++) {
      double p[3] = {((double)i + 0.5) * (double)sol->dom[0] / (double)w,
                     ((double)j + 0.5) * (double)sol->dom[1] / (double)h, z};
      double complex u = hz_sol3_eval(sol, p);
      intensity[(size_t)j * (size_t)w + (size_t)i] = creal(u * conj(u));
    }
  return 0;
}
