/* pblock — СЕРИАЛИЗАТОР СЦЕНЫ В ФОРМАТ HBLK (§882 v1, §902 v2): блочная
 * раскладка «блок = листовая клетка», порядок клеток Morton, CSR на диске.
 *
 * Файл: заголовок + таблица занятых клеток (cell_id Morton, start, count)
 * + tri-id кусков по клеткам + треугольники (9 double, исходный порядок).
 * v2 (магия «HBLK2»): после треугольников — kd[nt], lep[nt], E[nt]
 * (double): отклик и ПОЛЕ в файле — освещённый кадр из mmap без OBJ.
 * E берётся из sidecar §899 (ключ E=Ф); без E= — нули.
 * Читатель (pgather blk=) ходит DDA по сетке из заголовка и подтягивает
 * страницы mmap'ом — геометрия подгружается ходьбой (§2.4 STRUCTURE).
 *
 * Запуск: pblock <scene.obj> [scale=F] [lev=N] [E=ФАЙЛ] [out=ФАЙЛ]
 * lev — уровень листовой сетки (тот же смысл, что в pgather).
 */
#include <inttypes.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "scene_obj.h"

#define PB_MAGIC 0x334B4C4248ULL /* "HBLK3" little-endian (v3 §904) */

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
  const char *path = NULL, *outfile = "scene.hblk", *efile = NULL;
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
    else if (strncmp(argv[i], "E=", 2) == 0)
      efile = argv[i] + 2; /* §902: sidecar §899, nt doubles */
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
    /* §902: kd/lep/E per tri — отклик и поле в файле (кадр без OBJ) */
    double *kd = (double *)malloc((size_t)m.nt * sizeof *kd);
    double *lep = (double *)malloc((size_t)m.nt * sizeof *lep);
    double *E = (double *)calloc((size_t)m.nt, sizeof *E);
    if (!kd || !lep || !E) {
      fprintf(stderr, "pblock: нет памяти на kd/lep/E\n");
      return 2;
    }
    for (i = 0; i < m.nt; i++) {
      kd[i] = m.mtl[m.fm[i]].kd;
      lep[i] = (m.mtl[m.fm[i]].ke3[0] + m.mtl[m.fm[i]].ke3[1] + m.mtl[m.fm[i]].ke3[2]) / 3.0;
    }
    if (efile) { /* §899 sidecar: nt doubles, порядок ИСХОДНЫХ tri */
      FILE *fe = fopen(efile, "rb");
      if (!fe) {
        fprintf(stderr, "pblock: E не читается: %s\n", efile);
        return 2;
      }
      if (fread(E, sizeof(double), (size_t)m.nt, fe) != (size_t)m.nt) {
        fprintf(stderr, "pblock: E короче nt (%s)\n", efile);
        fclose(fe);
        return 2;
      }
      fclose(fe);
    }
    /* §904: kd3/lep3 per tri×канал — цветной кадр из файла (§889) */
    double *kd3 = (double *)malloc((size_t)m.nt * 3 * sizeof *kd3);
    double *lep3 = (double *)malloc((size_t)m.nt * 3 * sizeof *lep3);
    if (!kd3 || !lep3) {
      fprintf(stderr, "pblock: нет памяти на kd3/lep3\n");
      return 2;
    }
    for (i = 0; i < m.nt; i++)
      for (int c = 0; c < 3; c++) {
        kd3[3 * (int64_t)i + c] = m.mtl[m.fm[i]].kd3[c];
        lep3[3 * (int64_t)i + c] = m.mtl[m.fm[i]].ke3[c];
      }
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
    hdr[10] = 24 * 8;                        /* оффсет таблицы */
    hdr[11] = 24 * 8 + (uint64_t)nocc * 24;  /* оффсет tri-id (int32) */
    hdr[12] = hdr[11] + (uint64_t)nids * 4;  /* оффсет треугольников */
    hdr[13] = hdr[12] + (uint64_t)m.nt * 72; /* §902: оффсет kd */
    hdr[14] = hdr[13] + (uint64_t)m.nt * 8;  /* §902: оффсет lep */
    hdr[15] = hdr[14] + (uint64_t)m.nt * 8;  /* §902: оффсет E */
    hdr[16] = hdr[15] + (uint64_t)m.nt * 8;  /* §904: оффсет kd3 (3·nt) */
    hdr[17] = hdr[16] + (uint64_t)m.nt * 24; /* §904: оффсет lep3 (3·nt) */
    fwrite(hdr, 8, 24, f);
    fwrite(tab, sizeof *tab, (size_t)nocc, f);
    fwrite(ids, sizeof *ids, (size_t)nids, f);
    for (i = 0; i < m.nt; i++) {
      double p[3][3];
      hz_obj_tri(&m, (int32_t)i, p);
      fwrite(p, sizeof(double), 9, f);
    }
    fwrite(kd, sizeof *kd, (size_t)m.nt, f);
    fwrite(lep, sizeof *lep, (size_t)m.nt, f);
    fwrite(E, sizeof *E, (size_t)m.nt, f);
    fwrite(kd3, sizeof *kd3, (size_t)m.nt * 3, f);
    fwrite(lep3, sizeof *lep3, (size_t)m.nt * 3, f);
    fclose(f);
    printf("HBLK3: %s — nt=%d, занятых клеток %" PRId64 " из %" PRId64 ", ссылок %" PRId64
           " (кратность %.2f), E %s\n",
           outfile, m.nt, nocc, ncells, nids, (double)nids / (double)m.nt,
           efile ? "из sidecar" : "нули");
    free(kd);
    free(lep);
    free(E);
    free(kd3);
    free(lep3);
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
