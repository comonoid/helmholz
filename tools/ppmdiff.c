/* ppmdiff — РАЗНОСТЬ ДВУХ КАДРОВ PPM P6 (§848): ND-лестница на кадре.
 *
 * Кадры pgather записаны нормированными КАЖДЫЙ НА СВОЙ max (8 бит).
 * Сравнивать в этой шкале честно только после приведения к ОБЩЕМУ
 * максимуму: ppmdiff читает два файла, печатает отношение их максимумов
 * (масштабный множитель) и мерит разницу в единицах общей шкалы
 * [0,1] от наибольшего max. Нижний пол прибора — 8-битное квантование,
 * 1/255 ≈ 0.4 % (А1540).
 *
 * Печатается: средняя |Δ|/шкала, max |Δ|, доля пикселей с |Δ| > 10 %
 * шкалы, и max |Δ| по клеткам сетки 16×16 — где разница живёт.
 *
 * НК: сравнение кадра с самим собой даёт ровно 0; разные размеры —
 * отказ без сравнения.
 */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PPD_GRID 16 /* сетка локализации разницы */
#define PPD_MAXDIM 4096

static unsigned char *read_p6(const char *path, int *w, int *h, double *maxlum) {
  FILE *f = fopen(path, "rb");
  char magic[3] = {0};
  int W, H, maxv;
  unsigned char *buf;
  if (!f || fscanf(f, "%2s", magic) != 1 || strcmp(magic, "P6") != 0) {
    fprintf(stderr, "ppmdiff: %s не P6\n", path);
    if (f) fclose(f);
    return NULL;
  }
  if (fscanf(f, "%d %d %d", &W, &H, &maxv) != 3 || maxv != 255 || W < 1 || H < 1 ||
      W > PPD_MAXDIM || H > PPD_MAXDIM) {
    fprintf(stderr, "ppmdiff: %s — плохой заголовок\n", path);
    fclose(f);
    return NULL;
  }
  fgetc(f); /* единичный пробел после maxv */
  buf = (unsigned char *)malloc((size_t)W * (size_t)H * 3);
  if (!buf || fread(buf, 3, (size_t)W * (size_t)H, f) != (size_t)W * (size_t)H) {
    fprintf(stderr, "ppmdiff: %s обрезан\n", path);
    free(buf);
    fclose(f);
    return NULL;
  }
  fclose(f);
  *w = W;
  *h = H;
  /* кадры pgather серые (b == g == r) и нормированы на max кадра */
  {
    double m = 0;
    int64_t n = (int64_t)W * H;
    for (int64_t i = 0; i < n; i++) {
      if (buf[3 * i + 0] != buf[3 * i + 1] || buf[3 * i + 1] != buf[3 * i + 2]) {
        fprintf(stderr, "ppmdiff: %s не серый — не кадр pgather\n", path);
        free(buf);
        return NULL;
      }
      if (buf[3 * i] > m) m = buf[3 * i];
    }
    *maxlum = m;
  }
  return buf;
}

int main(int argc, char **argv) {
  int w1, h1, w2, h2;
  double m1, m2, scale;
  unsigned char *a, *b;
  if (argc != 3) {
    fprintf(stderr, "use: ppmdiff <кадр1.ppm> <кадр2.ppm>\n");
    return 2;
  }
  a = read_p6(argv[1], &w1, &h1, &m1);
  b = a ? read_p6(argv[2], &w2, &h2, &m2) : NULL;
  if (!b) {
    free(a);
    return 2;
  }
  if (w1 != w2 || h1 != h2) {
    fprintf(stderr, "ppmdiff: размеры разные (%dx%d vs %dx%d) — сравнения нет\n", w1, h1, w2, h2);
    free(a);
    free(b);
    return 2;
  }
  scale = (m1 > m2 ? m1 : m2);
  if (scale <= 0.0) {
    printf("PPMDIFF: оба кадра чёрные — разница 0\n");
    free(a);
    free(b);
    return 0;
  }
  {
    int64_t n = (int64_t)w1 * h1, nbig = 0;
    double ssum = 0, smax = 0;
    double gmax[PPD_GRID * PPD_GRID] = {0};
    for (int64_t i = 0; i < n; i++) {
      double d = fabs((double)a[3 * i] - (double)b[3 * i]) / scale;
      int gx = (int)((i % w1) * PPD_GRID / w1);
      int gy = (int)((i / w1) * PPD_GRID / h1);
      ssum += d;
      if (d > smax) smax = d;
      if (d > 0.10) nbig++;
      if (d > gmax[gy * PPD_GRID + gx]) gmax[gy * PPD_GRID + gx] = d;
    }
    printf("PPMDIFF: max кадров %.0f / %.0f (шкала = %.0f)\n", m1, m2, scale);
    printf("PPMDIFF: средняя |Δ| = %.4f шкалы (%.2f %%), max |Δ| = %.4f, пикселей >10%%: %.3f %%\n",
           ssum / (double)n, 100.0 * ssum / (double)n, smax, 100.0 * (double)nbig / (double)n);
    printf("PPMDIFF: max |Δ| по клеткам %dx%d (строки сверху вниз):\n", PPD_GRID, PPD_GRID);
    for (int gy = 0; gy < PPD_GRID; gy++) {
      for (int gx = 0; gx < PPD_GRID; gx++)
        printf("%s%4.0f", gx ? " " : "", 100.0 * gmax[gy * PPD_GRID + gx]);
      printf("\n");
    }
  }
  free(a);
  free(b);
  return 0;
}
