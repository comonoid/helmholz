/* kcalmk — §928: сборка таблицы перехвата k(r,ωbin) из дампов
 * эталонного КУСОЧНОГО прогона (HZ_KCALDUMP=path → path.num/path.den).
 *
 *   k(r,ω) = num[r][ω] / den[r][ω]
 *     num — Σ депозитов кусочного мира на детях rep r по ω-бину
 *           (последняя итерация эталона);
 *     den — Σ w·csec·axcos·Lin_entry по входам трубок в регион r
 *           по ω-бину (штамп «первое событие региона»).
 *
 * Разреженность — fallback-цепочка (§928-план, порог О3):
 *   ячейка (r,ω): den > DENTOT·1e-9 → k = num/den;
 *   иначе регион-среднее (r, всё ω): Σnum[r]/Σden[r] при Σden[r] > порога;
 *   иначе 0 (= формула §927 в свипе, fail-closed).
 *
 * use: kcalmk num=F den=F out=F [scramble=1]
 *   scramble — НК §928-П5: перестановка регионов, таблица ОБЯЗАНА
 *   разрушить баланс (kcal_mismatch/E прочь из коридора).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nstruct/sweep.h"

/* Порог знаменателя: доля от полного Σden (калибровка в энергоджоулях
 * прогона — абсолютный порог непереносим). 1e-9 — «ячейка не пустая»
 * против численного шума; ужесточление = больше fallback-ступени. */
#define KCAL_DENFRAC 1e-9

static double *load_tab(const char *path, int32_t *nrep, int32_t *nw, const char *tag) {
  hz_kcal_hdr h;
  FILE *f = fopen(path, "rb");
  double *a = NULL;
  if (f == NULL) {
    fprintf(stderr, "kcalmk: %s=%s не открылся\n", tag, path);
    return NULL;
  }
  if (fread(&h, sizeof h, 1, f) != 1 || memcmp(h.magic, "KCAL928", 8) != 0 || h.nrep <= 0 ||
      h.nw != HZ_KCAL_NW) {
    fprintf(stderr, "kcalmk: %s=%s: битый заголовок (nrep=%d nw=%d)\n", tag, path, (int)h.nrep,
            (int)h.nw);
    fclose(f);
    return NULL;
  }
  a = (double *)malloc((size_t)h.nrep * (size_t)h.nw * sizeof *a);
  if (a == NULL ||
      fread(a, sizeof *a, (size_t)h.nrep * (size_t)h.nw, f) != (size_t)h.nrep * (size_t)h.nw) {
    fprintf(stderr, "kcalmk: %s=%s: битый массив\n", tag, path);
    free(a);
    fclose(f);
    return NULL;
  }
  fclose(f);
  *nrep = h.nrep;
  *nw = h.nw;
  return a;
}

int main(int argc, char **argv) {
  const char *numpath = NULL, *denpath = NULL, *outpath = NULL;
  int scramble = 0;
  int32_t nrep = 0, nw = 0, nr2 = 0, nw2 = 0;
  double *num = NULL, *den = NULL, *out = NULL, dentot = 0.0;
  long ncell = 0, nrm = 0, nfall = 0;
  int32_t r;
  int w;
  for (int i = 1; i < argc; i++) {
    if (strncmp(argv[i], "num=", 4) == 0)
      numpath = argv[i] + 4;
    else if (strncmp(argv[i], "den=", 4) == 0)
      denpath = argv[i] + 4;
    else if (strncmp(argv[i], "out=", 4) == 0)
      outpath = argv[i] + 4;
    else if (strncmp(argv[i], "scramble=", 9) == 0)
      scramble = atoi(argv[i] + 9);
  }
  if (!numpath || !denpath || !outpath) {
    fprintf(stderr, "use: kcalmk num=F den=F out=F [scramble=1]\n");
    return 2;
  }
  num = load_tab(numpath, &nrep, &nw, "num");
  den = load_tab(denpath, &nr2, &nw2, "den");
  if (!num || !den) {
    free(num);
    free(den);
    return 2;
  }
  if (nrep != nr2) {
    fprintf(stderr, "kcalmk: nrep не совпал (%d против %d) — разные прогоны?\n", (int)nrep,
            (int)nr2);
    free(num);
    free(den);
    return 2;
  }
  for (r = 0; r < nrep * nw; r++)
    dentot += den[r];
  double dthr = dentot * KCAL_DENFRAC;
  out = (double *)calloc((size_t)nrep * (size_t)nw, sizeof *out);
  if (!out) {
    free(num);
    free(den);
    return 2;
  }
  for (r = 0; r < nrep; r++) {
    double rden = 0.0, rnum = 0.0;
    for (w = 0; w < nw; w++) {
      rden += den[(size_t)r * (size_t)nw + (size_t)w];
      rnum += num[(size_t)r * (size_t)nw + (size_t)w];
    }
    for (w = 0; w < nw; w++) {
      size_t c = (size_t)r * (size_t)nw + (size_t)w;
      if (den[c] > dthr) {
        out[c] = num[c] / den[c]; /* ячейка жива */
        ncell++;
      } else if (rden > dthr) {
        out[c] = rnum / rden; /* регион-среднее по ω */
        nrm++;
      } else {
        out[c] = 0.0; /* формула §927 в свипе (fail-closed) */
        nfall++;
      }
    }
  }
  if (scramble) { /* §928-П5: НК — перестановка регионов таблицей */
    for (r = 0; r < nrep; r++) {
      int32_t r2 = (int32_t)(((long long)r * 2654435761LL) % nrep);
      for (w = 0; w < nw; w++)
        out[(size_t)r * (size_t)nw + (size_t)w] = out[(size_t)r2 * (size_t)nw + (size_t)w];
    }
  }
  {
    hz_kcal_hdr h;
    FILE *f = fopen(outpath, "wb");
    memset(&h, 0, sizeof h);
    memcpy(h.magic, "KCAL928", sizeof h.magic);
    h.nrep = nrep;
    h.nw = nw;
    if (f == NULL || fwrite(&h, sizeof h, 1, f) != 1 ||
        fwrite(out, sizeof *out, (size_t)nrep * (size_t)nw, f) != (size_t)nrep * (size_t)nw) {
      fprintf(stderr, "kcalmk: out=%s не записался\n", outpath);
      if (f != NULL) fclose(f);
      free(num);
      free(den);
      free(out);
      return 2;
    }
    fclose(f);
  }
  printf("§928 kcalmk: nrep=%d nw=%d | ячеек живых %ld (%.1f%%), регион-средних %ld "
         "(%.1f%%), формула %ld (%.1f%%) | Σden=%.6g%s\n",
         (int)nrep, (int)nw, ncell, 100.0 * (double)ncell / ((double)nrep * nw), nrm,
         100.0 * (double)nrm / ((double)nrep * nw), nfall,
         100.0 * (double)nfall / ((double)nrep * nw), dentot,
         scramble ? " | SCRAMBLE (НК §928-П5)" : "");
  free(num);
  free(den);
  free(out);
  return 0;
}
