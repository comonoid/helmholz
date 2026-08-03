/* padj.c — реализация; довод, зачем отдельный файл, — в `padj.h`.
 *
 * ПЕРЕНОС ИЗ `tools/psil.c` СДЕЛАН БЕЗ ЕДИНОЙ ПРАВКИ АРИФМЕТИКИ, и это
 * проверяется числом: `./build/psil city noocc n=32` обязан дать те же
 * `8585954` рёбер и ту же медиану `20518`.
 */
#include "padj.h"
#include <stdlib.h>

void hz_padj_free(hz_padj *A) {
  free(A->ea);
  free(A->eb);
  free(A->ef);
  free(A->enf);
  A->ea = NULL;
  A->eb = NULL;
  A->ef = NULL;
  A->enf = NULL;
}

int hz_padj_build(hz_padj *A, const hz_objmesh *m) {
  A->ea = NULL;
  A->eb = NULL;
  A->ef = NULL;
  A->enf = NULL;
  A->ne = 0;
  A->nb1 = 0;
  A->nb2 = 0;
  A->nbm = 0;
  A->nover = 0;

  int64_t nslot = (int64_t)m->nt * 3;
  int64_t nvp1 = (int64_t)m->nv + 1;
  int64_t *off = calloc((size_t)nvp1 + 1, sizeof *off);
  if (off == NULL) return 2;
  for (int32_t t = 0; t < m->nt; t++)
    for (int i = 0; i < 3; i++) {
      int32_t a = m->f[(size_t)t * 3 + (size_t)i];
      int32_t b = m->f[(size_t)t * 3 + (size_t)((i + 1) % 3)];
      int32_t lo = a < b ? a : b;
      off[(int64_t)lo + 1]++;
    }
  for (int64_t i = 0; i < nvp1; i++)
    off[i + 1] += off[i];
  int32_t *shi = malloc((size_t)nslot * sizeof *shi);
  int32_t *sfa = malloc((size_t)nslot * sizeof *sfa);
  int64_t *cur = malloc((size_t)nvp1 * sizeof *cur);
  unsigned char *mark = calloc((size_t)nslot, 1);
  if (shi == NULL || sfa == NULL || cur == NULL || mark == NULL) {
    free(off);
    free(shi);
    free(sfa);
    free(cur);
    free(mark);
    return 2;
  }
  for (int64_t i = 0; i < nvp1; i++)
    cur[i] = off[i];
  for (int32_t t = 0; t < m->nt; t++)
    for (int i = 0; i < 3; i++) {
      int32_t a = m->f[(size_t)t * 3 + (size_t)i];
      int32_t b = m->f[(size_t)t * 3 + (size_t)((i + 1) % 3)];
      int32_t lo = a < b ? a : b, hi = a < b ? b : a;
      int64_t p = cur[lo]++;
      shi[p] = hi;
      sfa[p] = t;
    }

  /* Уникальные рёбра: в каждом ведре линейная группировка по старшей вершине.
   * Вёдер `nv`, записей в ведре в среднем `3·nt/nv` (у города около шести), так
   * что `O(k²)` внутри ведра дёшево и предсказуемо. */
  A->ea = malloc((size_t)nslot * sizeof *A->ea);
  A->eb = malloc((size_t)nslot * sizeof *A->eb);
  A->ef = malloc((size_t)nslot * HZ_PADJ_MAXF * sizeof *A->ef);
  A->enf = malloc((size_t)nslot);
  if (A->ea == NULL || A->eb == NULL || A->ef == NULL || A->enf == NULL) {
    free(off);
    free(shi);
    free(sfa);
    free(cur);
    free(mark);
    return 2;
  }
  int64_t ne = 0, nb1 = 0, nb2 = 0, nbm = 0, nover = 0;
  for (int32_t a = 0; a < m->nv; a++) {
    int64_t s = off[a], e = off[a + 1];
    for (int64_t i = s; i < e; i++) {
      if (mark[i]) continue;
      int32_t h = shi[i];
      int64_t nf = 1;
      A->ea[ne] = a;
      A->eb[ne] = h;
      A->ef[ne * HZ_PADJ_MAXF] = sfa[i];
      for (int64_t j = i + 1; j < e; j++) {
        if (mark[j] || shi[j] != h) continue;
        mark[j] = 1;
        if (nf < HZ_PADJ_MAXF) A->ef[ne * HZ_PADJ_MAXF + nf] = sfa[j];
        nf++;
      }
      if (nf > HZ_PADJ_MAXF) nover++;
      A->enf[ne] = (unsigned char)(nf > 255 ? 255 : nf);
      if (nf == 1)
        nb1++;
      else if (nf == 2)
        nb2++;
      else
        nbm++;
      ne++;
    }
  }
  free(mark);
  free(cur);
  free(off);
  free(shi);
  free(sfa);
  A->ne = ne;
  A->nb1 = nb1;
  A->nb2 = nb2;
  A->nbm = nbm;
  A->nover = nover;
  return 0;
}
