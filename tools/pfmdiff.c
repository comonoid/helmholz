/* pfmdiff.c — попиксельное сравнение двух буферов РАДИАНСА в формате PFM.
 *
 * ЛИНИЯ ПЕРЕНОСА, оснастка замеров К65 и К68.
 *
 * ЗАЧЕМ ОТДЕЛЬНЫЙ ИНСТРУМЕНТ, А НЕ ГЛАЗА. К19 измерила, что тон-маппинг съедает
 * метрику: та же ошибка после нормировки на максимум и гаммы 0.5 доходит до 255
 * уровней из 255. Поэтому сравнивать надо ЧИСЛА, и именно буфер радианса, а
 * PPM в метрику не входит вовсе. Отсюда же и формат: `hz_pfm_write` пишет
 * 32-битный float без нормировки, гаммы и раскраски.
 *
 * ДВЕ МЕТРИКИ, И ОБЕ ОБЯЗАТЕЛЬНЫ (К20). Пиксель — не точка, а на силуэте и на
 * кромке тени соседние пиксели дают «попал» и «не попал», то есть расхождение
 * там равно ВСЕЙ яркости независимо от того, насколько хороша схема. Мерить
 * максимум по всему кадру значит мерить кромки, а не схему. Поэтому кромки
 * выделяются и докладываются ОТДЕЛЬНО, а не подмешиваются в среднее.
 *
 * КАК ВЫДЕЛЯЕТСЯ КРОМКА, БЕЗ МАГИЧЕСКОГО ПОРОГА. Пиксель считается кромочным,
 * если локальный размах ЭТАЛОНА в его 3×3 окрестности превышает размах в этой
 * же окрестности у сглаженной картины более чем вдвое — то есть признак берётся
 * из САМОГО эталона, а не из сравниваемой пары, и от величины расхождения не
 * зависит. Порога в единицах радианса здесь нет: сравниваются две величины
 * одной размерности.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Чтение PFM. Порядок строк снизу вверх, отрицательный масштаб = little-endian —
 * часть формата, а не наша прихоть; порядок байт ОПРЕДЕЛЯЕТСЯ на месте. */
static float *pfm_read(const char *path, int *w, int *h, int *nch) {
  FILE *f = fopen(path, "rb");
  if (f == NULL) return NULL;
  char magic[3] = {0, 0, 0};
  if (fscanf(f, "%2s", magic) != 1) {
    fclose(f);
    return NULL;
  }
  int c = 3;
  if (strcmp(magic, "PF") == 0)
    c = 3;
  else if (strcmp(magic, "Pf") == 0)
    c = 1;
  else {
    fclose(f);
    return NULL;
  }
  int ww = 0, hh = 0;
  double scale = 0.0;
  if (fscanf(f, "%d %d %lf", &ww, &hh, &scale) != 3 || ww <= 0 || hh <= 0) {
    fclose(f);
    return NULL;
  }
  if (fgetc(f) == EOF) {
    fclose(f);
    return NULL;
  }
  size_t n = (size_t)ww * (size_t)hh * (size_t)c;
  float *d = malloc(n * sizeof(float));
  if (d == NULL) {
    fclose(f);
    return NULL;
  }
  if (fread(d, sizeof(float), n, f) != n) {
    free(d);
    fclose(f);
    return NULL;
  }
  /* свидетель порядка байт: scale < 0 — little-endian. На big-endian машине
   * пришлось бы переставлять; здесь проверяется, а не предполагается. */
  if (scale > 0.0) {
    unsigned char *b = (unsigned char *)d;
    for (size_t i = 0; i < n; i++) {
      unsigned char t;
      t = b[4 * i + 0];
      b[4 * i + 0] = b[4 * i + 3];
      b[4 * i + 3] = t;
      t = b[4 * i + 1];
      b[4 * i + 1] = b[4 * i + 2];
      b[4 * i + 2] = t;
    }
  }
  fclose(f);
  *w = ww;
  *h = hh;
  *nch = c;
  return d;
}

/* локальный размах в окне 3×3 */
static double span3(const float *a, int w, int h, int nch, int x, int y) {
  double lo = 1e300, hi = -1e300;
  for (int dy = -1; dy <= 1; dy++)
    for (int dx = -1; dx <= 1; dx++) {
      int xx = x + dx, yy = y + dy;
      if (xx < 0 || yy < 0 || xx >= w || yy >= h) continue;
      double v = a[((size_t)yy * (size_t)w + (size_t)xx) * (size_t)nch];
      if (v < lo) lo = v;
      if (v > hi) hi = v;
    }
  return hi - lo;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "pfmdiff ЭТАЛОН СРАВНИВАЕМЫЙ\n");
    return 2;
  }
  int w0 = 0, h0 = 0, c0 = 0, w1 = 0, h1 = 0, c1 = 0;
  float *a = pfm_read(argv[1], &w0, &h0, &c0);
  float *b = pfm_read(argv[2], &w1, &h1, &c1);
  if (a == NULL || b == NULL) {
    fprintf(stderr, "не прочитался: %s\n", a == NULL ? argv[1] : argv[2]);
    free(a);
    free(b);
    return 2;
  }
  if (w0 != w1 || h0 != h1) {
    fprintf(stderr, "размеры не совпали: %dx%d против %dx%d\n", w0, h0, w1, h1);
    free(a);
    free(b);
    return 2;
  }

  /* СРЕДНИЙ РАЗМАХ ЭТАЛОНА — масштаб, относительно которого решается, кромка
   * это или гладкое место. Берётся из эталона и от сравнения не зависит. */
  double sp_sum = 0.0;
  size_t np = (size_t)w0 * (size_t)h0;
  for (int y = 0; y < h0; y++)
    for (int x = 0; x < w0; x++)
      sp_sum += span3(a, w0, h0, c0, x, y);
  double sp_mean = sp_sum / (double)np;

  double emax_in = 0.0, esum_in = 0.0, emax_ed = 0.0;
  double rmax_in = 0.0, rsum_in = 0.0;
  size_t nin = 0, ned = 0;
  double amin = 1e300, amax = -1e300;
  for (int y = 0; y < h0; y++)
    for (int x = 0; x < w0; x++) {
      size_t i = ((size_t)y * (size_t)w0 + (size_t)x);
      double va = a[i * (size_t)c0], vb = b[i * (size_t)c1];
      if (va < amin) amin = va;
      if (va > amax) amax = va;
      double e = fabs(va - vb);
      int edge = span3(a, w0, h0, c0, x, y) > 2.0 * sp_mean;
      if (edge) {
        ned++;
        if (e > emax_ed) emax_ed = e;
      } else {
        nin++;
        esum_in += e;
        if (e > emax_in) emax_in = e;
        double r = fabs(va) > 0.0 ? e / fabs(va) : 0.0;
        rsum_in += r;
        if (r > rmax_in) rmax_in = r;
      }
    }

  printf("%dx%d, эталон %.6f…%.6f, средний размах 3x3 %.4e\n", w0, h0, amin, amax, sp_mean);
  printf("ВНЕ КРОМОК (%zu пикселей, %.1f%%): |Δ| сред %.6e, макс %.6e; отн сред %.6e, макс %.6e\n",
         nin, 100.0 * (double)nin / (double)np, esum_in / (double)nin, emax_in,
         rsum_in / (double)nin, rmax_in);
  printf("НА КРОМКАХ (%zu пикселей, %.1f%%): |Δ| макс %.6e\n", ned,
         100.0 * (double)ned / (double)np, emax_ed);
  free(a);
  free(b);
  return 0;
}
