#include "image.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* 7-stop inferno-like ramp, linear interpolation */
static const double RAMP[7][3] = {{0.0, 0.0, 0.02},   {0.11, 0.05, 0.30}, {0.42, 0.09, 0.43},
                                  {0.72, 0.21, 0.33}, {0.93, 0.44, 0.17}, {0.99, 0.75, 0.13},
                                  {0.99, 0.99, 0.75}};

static void colormap(double v, unsigned char px[3]) {
  if (v < 0.0) v = 0.0;
  if (v > 1.0) v = 1.0;
  double s = v * 6.0;
  int i = (int)s;
  if (i > 5) i = 5;
  double f = s - (double)i;
  for (int c = 0; c < 3; c++) {
    double x = RAMP[i][c] * (1.0 - f) + RAMP[i + 1][c] * f;
    px[c] = (unsigned char)(255.0 * x + 0.5);
  }
}

/* PFM: радианс как есть, без тон-маппинга. Разбор — в image.h.
 * `g` и `b` могут совпадать с `r` (один спектральный канал). */
int hz_pfm_write(const char *path, const double *r, const double *g, const double *b, int w,
                 int h) {
  if (w <= 0 || h <= 0) return 1;
  FILE *f = fopen(path, "wb");
  if (f == NULL) return 1;
  /* ПОРЯДОК БАЙТ ОПРЕДЕЛЯЕТСЯ, А НЕ ПРЕДПОЛАГАЕТСЯ: в PFM он объявлен знаком
   * масштаба, и соврать здесь значит отдать файл, который прочтётся наизнанку. */
  unsigned int probe = 1u;
  unsigned char pb;
  memcpy(&pb, &probe, 1);
  double scale = pb ? -1.0 : 1.0;
  fprintf(f, "PF\n%d %d\n%.1f\n", w, h, scale);
  const double *ch[3] = {r, g, b};
  int bad = 0;
  /* СТРОКИ СНИЗУ ВВЕРХ — это часть формата PFM, а не наш порядок хранения. */
  for (int y = h - 1; y >= 0 && !bad; y--)
    for (int x = 0; x < w; x++) {
      float px[3];
      for (int c = 0; c < 3; c++)
        px[c] = (float)ch[c][(size_t)y * (size_t)w + (size_t)x];
      if (fwrite(px, sizeof(float), 3, f) != 3) bad = 1;
    }
  if (fclose(f) != 0) bad = 1;
  return bad;
}

static int cmp_dbl(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

int hz_ppm_write(const char *path, const double *intensity, int w, int h) {
  /* white point at the 99.5th percentile: a single bright source must not
   * drown the rest of the frame */
  int n = w * h;
  double mx = 0.0;
  double *tmp = malloc((size_t)n * sizeof(double));
  if (tmp != NULL) {
    for (int i = 0; i < n; i++)
      tmp[i] = intensity[i];
    qsort(tmp, (size_t)n, sizeof(double), cmp_dbl);
    mx = tmp[(size_t)((double)n * 0.995)];
    free(tmp);
  }
  if (!(mx > 0.0)) mx = 1.0;
  FILE *f = fopen(path, "wb");
  if (f == NULL) return 1;
  fprintf(f, "P6\n%d %d\n255\n", w, h);
  unsigned char *row = malloc((size_t)w * 3);
  if (row == NULL) {
    fclose(f);
    return 1;
  }
  int rc = 0;
  for (int y = 0; y < h; y++) {
    for (int x = 0; x < w; x++) {
      double v = sqrt(intensity[(size_t)y * (size_t)w + (size_t)x] / mx); /* gamma 0.5 */
      colormap(v, &row[3 * x]);
    }
    if (fwrite(row, 3, (size_t)w, f) != (size_t)w) rc = 1;
  }
  free(row);
  fclose(f);
  return rc;
}

int hz_ppm_write_rgb(const char *path, const unsigned char *rgb, int w, int h) {
  FILE *f = fopen(path, "wb");
  if (f == NULL) return 1;
  fprintf(f, "P6\n%d %d\n255\n", w, h);
  size_t n = (size_t)w * (size_t)h * 3;
  int rc = fwrite(rgb, 1, n, f) == n ? 0 : 1;
  fclose(f);
  return rc;
}

int hz_ppm_write_signed(const char *path, const double *field, int w, int h) {
  int n = w * h;
  double mx = 0.0;
  double *tmp = malloc((size_t)n * sizeof(double));
  if (tmp != NULL) {
    for (int i = 0; i < n; i++)
      tmp[i] = fabs(field[i]);
    qsort(tmp, (size_t)n, sizeof(double), cmp_dbl);
    mx = tmp[(size_t)((double)n * 0.995)];
    free(tmp);
  }
  if (!(mx > 0.0)) mx = 1.0;
  unsigned char *rgb = malloc((size_t)n * 3);
  if (rgb == NULL) return 1;
  for (int i = 0; i < n; i++) {
    double v = field[i] / mx;
    if (v > 1.0) v = 1.0;
    if (v < -1.0) v = -1.0;
    double r, g, b;
    if (v >= 0.0) { /* white -> red */
      r = 1.0;
      g = 1.0 - v;
      b = 1.0 - v;
    } else { /* white -> blue */
      r = 1.0 + v;
      g = 1.0 + v;
      b = 1.0;
    }
    rgb[3 * i] = (unsigned char)(255.0 * r + 0.5);
    rgb[3 * i + 1] = (unsigned char)(255.0 * g + 0.5);
    rgb[3 * i + 2] = (unsigned char)(255.0 * b + 0.5);
  }
  int rc = hz_ppm_write_rgb(path, rgb, w, h);
  free(rgb);
  return rc;
}

/* Ш8 (§575): читатель PPM `P6` — см. заголовок. Комментарии по одному пробелу
 * между полями допускаются (`#` до конца строки), потому что их пишет
 * ImageMagick. */
static int ppm_tok(FILE *f, long *out) {
  int c;
  for (;;) {
    do {
      c = fgetc(f);
    } while (c == 32 || c == 9 || c == 10 || c == 13);
    if (c == 35) { /* комментарий до конца строки */
      do {
        c = fgetc(f);
      } while (c != 10 && c != 13 && c != EOF);
      continue;
    }
    break;
  }
  if (c < 48 || c > 57) return -1;
  long v = 0;
  while (c >= 48 && c <= 57) {
    if (v > 100000000L) return -1; /* битый вход, а не большая картинка */
    v = v * 10 + (c - 48);
    c = fgetc(f);
  }
  *out = v;
  return 0;
}

int hz_ppm_read(const char *path, unsigned char **rgb, int *w, int *h) {
  *rgb = NULL;
  *w = *h = 0;
  FILE *f = fopen(path, "rb");
  if (f == NULL) return 1;
  int c1 = fgetc(f), c2 = fgetc(f);
  if (c1 != 80 || c2 != 54) { /* "P6" */
    fclose(f);
    return 2;
  }
  long ww = 0, hh = 0, mx = 0;
  if (ppm_tok(f, &ww) != 0 || ppm_tok(f, &hh) != 0 || ppm_tok(f, &mx) != 0) {
    fclose(f);
    return 3;
  }
  if (ww <= 0 || hh <= 0 || mx != 255 || ww > 16384 || hh > 16384) {
    fclose(f);
    return 4;
  }
  size_t n = (size_t)ww * (size_t)hh * 3;
  unsigned char *b = malloc(n);
  if (b == NULL) {
    fclose(f);
    return 5;
  }
  if (fread(b, 1, n, f) != n) {
    free(b);
    fclose(f);
    return 6;
  }
  fclose(f);
  *rgb = b;
  *w = (int)ww;
  *h = (int)hh;
  return 0;
}
