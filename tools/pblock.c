/* pblock — СЕРИАЛИЗАТОР СЦЕНЫ В ФОРМАТ HBLK v1 (§882): блочная раскладка
 * «блок = листовая клетка», порядок клеток Morton, CSR кусков на диске.
 *
 * Файл: заголовок + таблица занятых клеток (cell_id Morton, start, count)
 * + tri-id кусков по клеткам + треугольники (9 double, исходный порядок).
 * Читатель (pgather blk=) ходит DDA по сетке из заголовка и подтягивает
 * страницы mmap'ом — геометрия подгружается ходьбой (§2.4 STRUCTURE).
 *
 * Запуск: pblock <scene.obj> [scale=F] [lev=N] [out=ФАЙЛ]
 * lev — уровень листовой сетки (тот же смысл, что в pgather).
 */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "scene_obj.h"

#define PB_MAGIC 0x314B4C4248ULL /* "HBLK1" little-endian */

static int64_t pb_morton3(uint32_t x, uint32_t y, uint32_t z) {
  int64_t r = 0;
  int b;
  for (b = 0; b < 21; b++)
    r |= ((int64_t)(x >> b & 1) << (3 * b)) | ((int64_t)(y >> b & 1) << (3 * b + 1)) |
         ((int64_t)(z >> b & 1) << (3 * b + 2));
  return r;
}

typedef struct {
  int64_t cell; /* Morton-ключ */
  int64_t start;
  int32_t cnt;
} pb_cellrec;

static int pb_cellcmp(const void *a, const void *b) {
  const pb_cellrec *x = a, *y = b;
  return x->cell < y->cell ? -1 : (x->cell > y->cell ? 1 : 0);
}

int main(int argc, char **argv) {
  const char *path = NULL, *outfile = "scene.hblk";
  double scale = 1.0;
  int lev = 6, i, ax;
  hz_objmesh m;
  for (i = 1; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0)
      lev = atoi(argv[i] + 4);
    else if (strncmp(argv[i], "scale=", 6) == 0)
      scale = atof(argv[i] + 6);
    else if (strncmp(argv[i], "out=", 4) == 0)
      outfile = argv[i] + 4;
    else
      path = argv[i];
  }
  if (!path || lev < 1 || lev > 20) {
    fprintf(stderr, "use: pblock <scene.obj> [scale=F] [lev=N] [out=ФАЙЛ]\n");
    return 2;
  }
  if (hz_obj_load(&m, path, scale) != 0) {
    fprintf(stderr, "pblock: не читается %s\n", path);
    return 2;
  }
  double cell, lo[3];
  int64_t nx, ny, nz;
  {
    double maxdim = 0.0;
    for (ax = 0; ax < 3; ax++) {
      double s = m.hi[ax] - m.lo[ax];
      lo[ax] = m.lo[ax];
      if (s > maxdim) maxdim = s;
    }
    cell = maxdim / (double)(1 << lev);
    nx = (int64_t)(1 << lev);
    ny = nx;
    nz = nx;
  }
  /* подсчёт: клетка занята, если bbox куска её пересекает (та же математика,
   * что pg_bbox_span — А1622) */
  int64_t ncells = nx * ny * nz;
  uint8_t *occ = (uint8_t *)calloc((size_t)ncells, 1);
  int64_t *cnt = (int64_t *)calloc((size_t)ncells, sizeof *cnt);
  int32_t *cntp = (int32_t *)malloc((size_t)m.nt * sizeof *cntp);
  if (!occ || !cnt || !cntp) {
    fprintf(stderr, "pblock: нет памяти\n");
    return 2;
  }
  for (i = 0; i < m.nt; i++) {
    double p[3][3];
    hz_obj_tri(&m, (int32_t)i, p);
    int64_t i0[3], i1[3];
    for (ax = 0; ax < 3; ax++) {
      double a = p[0][ax], b = p[0][ax];
      for (int v = 1; v < 3; v++) {
        if (p[v][ax] < a) a = p[v][ax];
        if (p[v][ax] > b) b = p[v][ax];
      }
      int64_t n = ax == 0 ? nx : (ax == 1 ? ny : nz);
      i0[ax] = (int64_t)((a - lo[ax]) / cell);
      i1[ax] = (int64_t)((b - lo[ax]) / cell);
      if (i0[ax] < 0) i0[ax] = 0;
      if (i1[ax] >= n) i1[ax] = n - 1;
    }
    int32_t c = 0;
    for (int64_t k = i0[2]; k <= i1[2]; k++)
      for (int64_t j = i0[1]; j <= i1[1]; j++)
        for (int64_t g = i0[0]; g <= i1[0]; g++) {
          int64_t id = g + nx * (j + ny * k);
          if (!occ[id]) occ[id] = 1;
          cnt[id]++;
          c++;
        }
    cntp[i] = c;
  }
  /* таблица занятых клеток в Morton-порядке (А1625) */
  int64_t nocc = 0;
  for (int64_t id = 0; id < ncells; id++)
    if (occ[id]) nocc++;
  pb_cellrec *tab = (pb_cellrec *)malloc((size_t)nocc * sizeof *tab);
  if (!tab) return 2;
  {
    int64_t w = 0;
    for (int64_t id = 0; id < ncells; id++)
      if (occ[id]) {
        int64_t g = id % nx, j = (id / nx) % ny, k = id / (nx * ny);
        tab[w].cell = pb_morton3((uint32_t)g, (uint32_t)j, (uint32_t)k);
        tab[w].cnt = (int32_t)cnt[id];
        tab[w].start = 0;
        w++;
      }
  }
  qsort(tab, (size_t)nocc, sizeof *tab, pb_cellcmp);
  { /* start — префиксные суммы; внутри клетки куски в ИСХОДНОМ порядке p
     * (counting-sort, А1622) */
    int64_t acc = 0;
    for (int64_t c = 0; c < nocc; c++) {
      tab[c].start = acc;
      acc += tab[c].cnt;
    }
  }
  int64_t nids = 0;
  for (int64_t c = 0; c < nocc; c++)
    nids += tab[c].cnt;
  /* nids = Σ ссылок (кусок лежит в нескольких клетках bbox) — НЕ nt;
   * аллокация под nt = heap-buffer-overflow (ASAN, §882) */
  int32_t *ids = (int32_t *)malloc((size_t)nids * sizeof *ids);
  int64_t *fill = (int64_t *)malloc((size_t)nocc * sizeof *fill);
  if (!ids || !fill) return 2;
  for (int64_t c = 0; c < nocc; c++)
    fill[c] = tab[c].start;
  for (i = 0; i < m.nt; i++) {
    double p[3][3];
    hz_obj_tri(&m, (int32_t)i, p);
    int64_t i0[3], i1[3];
    for (ax = 0; ax < 3; ax++) {
      double a = p[0][ax], b = p[0][ax];
      for (int v = 1; v < 3; v++) {
        if (p[v][ax] < a) a = p[v][ax];
        if (p[v][ax] > b) b = p[v][ax];
      }
      int64_t n = ax == 0 ? nx : (ax == 1 ? ny : nz);
      i0[ax] = (int64_t)((a - lo[ax]) / cell);
      i1[ax] = (int64_t)((b - lo[ax]) / cell);
      if (i0[ax] < 0) i0[ax] = 0;
      if (i1[ax] >= n) i1[ax] = n - 1;
    }
    for (int64_t k = i0[2]; k <= i1[2]; k++)
      for (int64_t j = i0[1]; j <= i1[1]; j++)
        for (int64_t g = i0[0]; g <= i1[0]; g++) {
          pb_cellrec key;
          key.cell = pb_morton3((uint32_t)g, (uint32_t)j, (uint32_t)k);
          pb_cellrec *hit = (pb_cellrec *)bsearch(&key, tab, (size_t)nocc, sizeof *tab, pb_cellcmp);
          if (!hit) {
            fprintf(stderr, "pblock: клетка не найдена: key=%lld\n", (long long)key.cell);
            return 3;
          }
          ids[fill[hit - tab]++] = (int32_t)i;
        }
  }
  { /* запись */
    FILE *f = fopen(outfile, "wb");
    if (!f) {
      fprintf(stderr, "pblock: не открыть %s\n", outfile);
      return 2;
    }
    uint64_t hdr[24] = {0};
    hdr[0] = PB_MAGIC;
    hdr[1] = (uint64_t)m.nt;
    hdr[2] = (uint64_t)nx;
    hdr[3] = (uint64_t)ny;
    hdr[4] = (uint64_t)nz;
    memcpy(&hdr[5], lo, 3 * sizeof(double));
    memcpy(&hdr[8], &cell, sizeof(double));
    hdr[9] = (uint64_t)nocc;
    hdr[10] = 24 * 8;                       /* оффсет таблицы */
    hdr[11] = 24 * 8 + (uint64_t)nocc * 24; /* оффсет tri-id (int32) */
    hdr[12] = hdr[11] + (uint64_t)nids * 4; /* оффсет треугольников */
    fwrite(hdr, 8, 24, f);
    fwrite(tab, sizeof *tab, (size_t)nocc, f);
    fwrite(ids, sizeof *ids, (size_t)nids, f);
    for (i = 0; i < m.nt; i++) {
      double p[3][3];
      hz_obj_tri(&m, (int32_t)i, p);
      fwrite(p, sizeof(double), 9, f);
    }
    fclose(f);
    printf("HBLK: %s — nt=%d, занятых клеток %" PRId64 " из %" PRId64 ", ссылок %" PRId64
           " (кратность %.2f)\n",
           outfile, m.nt, nocc, ncells, nids, (double)nids / (double)m.nt);
  }
  free(occ);
  free(cnt);
  free(cntp);
  free(tab);
  free(ids);
  free(fill);
  hz_obj_free(&m);
  return 0;
}
