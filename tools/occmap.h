/* occmap.h — ФОРМАТ КАРТЫ ЗАНЯТОСТИ (PLAN_ELEMENTS.md §385, шаг Ш0).
 *
 * ЗАЧЕМ ОТДЕЛЬНЫМ ФАЙЛОМ, А НЕ КОПИЕЙ В КАЖДОМ СТЕНДЕ. Сверяются ДВА
 * инструмента: `tools/poccref.c` (эталон, переживает Ш1б) и `tools/pfield.c`
 * (ключ `occdump`, одноразовый — уходит вместе с адаптером, А673). Разошедшийся
 * формат дал бы расхождение множеств, неотличимое от расхождения ПУТЕЙ, то есть
 * ровно ту подмену величины, ради снятия которой шаг и делается.
 *
 * РАСКЛАДКА ЯЧЕЙКИ В ЛИНЕЙНЫЙ ИНДЕКС — `(z·n + y)·n + x`, ОДНА НА ОБА
 * ИНСТРУМЕНТА. Разная раскладка даёт ТОЧНОЕ совпадение чисел занятых при
 * полностью разошедшихся множествах (А684); ловится это трёхсторонней сверкой
 * плюс СЛЕПКОМ — габаритом занятого множества в индексах ячеек, который лежит в
 * заголовке рядом с числом.
 */
#ifndef HZ_OCCMAP_H
#define HZ_OCCMAP_H

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#define HZ_OCC_MAGIC 0x31434F48u /* "HOC1" */

typedef struct {
  uint32_t magic;
  int32_t lev;
  int32_t n;
  int32_t pad; /* выравнивание под int64 ниже — явно, а не молчанием компилятора */
  int64_t count;
  int32_t blo[3], bhi[3]; /* габарит занятого множества В ИНДЕКСАХ ЯЧЕЕК */
} hz_occ_hdr;

static inline size_t hz_occ_bytes(size_t ncell) {
  return (ncell + 7u) / 8u;
}

static inline void hz_occ_set(unsigned char *b, size_t i) {
  b[i >> 3] |= (unsigned char)(1u << (i & 7u));
}

static inline int hz_occ_get(const unsigned char *b, size_t i) {
  return (b[i >> 3] >> (i & 7u)) & 1u;
}

static inline size_t hz_occ_index(int32_t n, int64_t x, int64_t y, int64_t z) {
  return ((size_t)z * (size_t)n + (size_t)y) * (size_t)n + (size_t)x;
}

/* СЛЕПОК: число и габарит. Считается ОДНОЙ функцией для обеих сторон сверки —
 * иначе разъедется он, а не множества. */
static inline void hz_occ_stamp(int lev, int32_t n, const unsigned char *bit, hz_occ_hdr *h) {
  h->magic = HZ_OCC_MAGIC;
  h->lev = lev;
  h->n = n;
  h->pad = 0;
  h->count = 0;
  for (int a = 0; a < 3; a++) {
    h->blo[a] = n;
    h->bhi[a] = -1;
  }
  for (int64_t z = 0; z < n; z++)
    for (int64_t y = 0; y < n; y++)
      for (int64_t x = 0; x < n; x++) {
        if (!hz_occ_get(bit, hz_occ_index(n, x, y, z))) continue;
        h->count++;
        int32_t c[3] = {(int32_t)x, (int32_t)y, (int32_t)z};
        for (int a = 0; a < 3; a++) {
          if (c[a] < h->blo[a]) h->blo[a] = c[a];
          if (c[a] > h->bhi[a]) h->bhi[a] = c[a];
        }
      }
  if (h->count == 0)
    for (int a = 0; a < 3; a++) {
      h->blo[a] = 0;
      h->bhi[a] = -1;
    }
}

static inline int hz_occ_write(const char *path, int lev, int32_t n, const unsigned char *bit) {
  hz_occ_hdr h;
  hz_occ_stamp(lev, n, bit, &h);
  FILE *f = fopen(path, "wb");
  if (f == NULL) return 1;
  size_t nb = hz_occ_bytes((size_t)n * (size_t)n * (size_t)n);
  int bad = (fwrite(&h, sizeof h, 1, f) != 1) || (fwrite(bit, 1, nb, f) != nb);
  if (fclose(f) != 0) bad = 1;
  return bad;
}

/* РАЗБОРЩИК НЕДОВЕРЕННОГО ВХОДА (А686). Выделен чистой функцией без
 * ввода-вывода именно затем, чтобы его брал CBMC: ничего не выделяет и никуда не
 * пишет, кроме `*out`. Возврат `0` — годен, иначе код причины. */
static inline int hz_occ_parse(const unsigned char *buf, size_t nbuf, int lev, hz_occ_hdr *out) {
  if (buf == NULL || out == NULL) return 1;
  if (nbuf < sizeof(hz_occ_hdr)) return 2;
  hz_occ_hdr h;
  memcpy(&h, buf, sizeof h);
  if (h.magic != HZ_OCC_MAGIC) return 3;
  if (h.lev < 1 || h.lev > 20) return 4; /* ДО сдвига: иначе сдвиг неопределён */
  if (h.n != (int32_t)1 << h.lev) return 5;
  if (h.lev != lev) return 6;
  size_t ncell = (size_t)h.n * (size_t)h.n * (size_t)h.n;
  if (nbuf - sizeof(hz_occ_hdr) != hz_occ_bytes(ncell)) return 7;
  if (h.count < 0 || (uint64_t)h.count > (uint64_t)ncell) return 8;
  for (int a = 0; a < 3; a++)
    if (h.count > 0 && (h.blo[a] < 0 || h.bhi[a] >= h.n || h.blo[a] > h.bhi[a])) return 9;
  *out = h;
  return 0;
}

#endif
