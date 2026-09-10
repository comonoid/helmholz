/* pyr.c — построение пирамиды шагов переноса, шаг 1 (§833). См. pyr.h. */
#include "pyr.h"

#include <float.h>
#include <stdlib.h>
#include <string.h>

_Static_assert(sizeof(hz_piece) == 16, "бюджет STRUCTURE.md §2.3: кусок 16 Б");
_Static_assert(sizeof(hz_pyr_leaf) == 12, "бюджет STRUCTURE.md §2.1: лист 12 Б");
_Static_assert(sizeof(hz_pyr_node) == 16, "бюджет STRUCTURE.md §2.1: узел уровня 16 Б");

/* Одна ось: индекс клетки с зажимом, тотальный по входу (CBMC). NaN даёт
 * !(v >= 0) → 0; каст (int64_t) стоит только под двойным неравенством, так
 * что значение заведомо представимо. */
static int pyr_finite(double v) {
  /* Конечность одними сравнениями: CBMC не знает __builtin_isfinite, а
   * приём v−v == 0 сам порождает NaN на бесконечности; ноль проходит
   * через диапазон и без отдельного сравнения (Werror=float-equal). */
  return v > -DBL_MAX && v < DBL_MAX;
}

static int64_t pyr_axis_index(double v, double lo, double cell, int32_t n) {
  double d, t;
  if (!(cell > 0.0) || !(n > 0)) return 0;
  if (!pyr_finite(v) || !pyr_finite(lo)) return 0;
  d = v - lo;
  if (!pyr_finite(d)) return 0; /* переполнение вычитания — в клетку 0 */
  t = d / cell;                 /* конечное/положительное: NaN недостижим, ±inf ловит зажим */
  if (!(t >= 0.0)) return 0;
  if (t > (double)(n - 1)) return (int64_t)n - 1;
  return (int64_t)t;
}

int64_t hz_pyr_leaf_index(const double lo[3], double cell, int32_t nx, int32_t ny, int32_t nz,
                          const double p[3]) {
  int64_t i, j, k;
  /* NULL-аргументы — контракт вызова (hz_pyr_build проверяет), но функция
   * тотальна: возвращаем -1, вызов через pyr_find просто не найдёт клетку. */
  if (!lo || !p) return -1;
  i = pyr_axis_index(p[0], lo[0], cell, nx);
  j = pyr_axis_index(p[1], lo[1], cell, ny);
  k = pyr_axis_index(p[2], lo[2], cell, nz);
  return i + (int64_t)nx * (j + (int64_t)ny * k);
}

int64_t hz_pyr_wrap(int64_t c, int64_t n) {
  /* Контракт: 1 <= n <= INT64_MAX/2 (построение держит ncells ≤ 2^60).
   * Без guard'а CBMC находит деление на ноль и INT64_MIN % (-1). */
  int64_t r;
  if (n <= 0) return 0;
  r = c % n;
  if (r < 0) r += n; /* r ∈ (−n, n), n ≤ INT64_MAX/2 ⇒ r+n не переполняется */
  return r;
}

static int pyr_cmp_i64(const void *a, const void *b) {
  int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
  return (x > y) - (x < y);
}

/* двоичный поиск id в отсортированном массиве; -1 если нет */
static int32_t pyr_find(const int64_t *a, int32_t n, int64_t x) {
  int32_t lo = 0, hi = n;
  while (lo < hi) {
    int32_t mid = lo + (hi - lo) / 2;
    if (a[mid] < x)
      lo = mid + 1;
    else
      hi = mid;
  }
  if (lo < n && a[lo] == x) return lo;
  return -1;
}

/* Диапазон клеток сетки листа, пересекаемых bbox-ом треугольника t. */
static void pyr_bbox_span(const hz_pyr *py, int32_t t, const double *tri_min, const double *tri_max,
                          int64_t x[2], int64_t y[2], int64_t z[2]) {
  int64_t i0 =
      hz_pyr_leaf_index(py->lo, py->cell, py->nx, py->ny, py->nz, tri_min + 3 * (int64_t)t);
  int64_t i1 =
      hz_pyr_leaf_index(py->lo, py->cell, py->nx, py->ny, py->nz, tri_max + 3 * (int64_t)t);
  x[0] = i0 % py->nx;
  x[1] = i1 % py->nx;
  y[0] = (i0 / py->nx) % py->ny;
  y[1] = (i1 / py->nx) % py->ny;
  z[0] = i0 / ((int64_t)py->nx * py->ny);
  z[1] = i1 / ((int64_t)py->nx * py->ny);
}

int hz_pyr_build(hz_pyr *py, int32_t nt, const double *tri_min, const double *tri_max,
                 const double *centroid, const int32_t *mtl, const double lo[3], const double hi[3],
                 double cell) {
  int32_t t, l;
  double size[3];
  int64_t *marks = NULL; /* битовая карта bbox-занятости, ВРЕМЕННАЯ (в бюджет не входит) */
  int64_t *cnt = NULL;
  int32_t *fill = NULL;
  int64_t *cur = NULL, *par = NULL;
  int32_t ncur;
  int64_t pnx, pny, pnz;

  memset(py, 0, sizeof *py);
  if (nt <= 0 || !(cell > 0.0) || !tri_min || !tri_max || !centroid || !mtl) return 1;

  for (l = 0; l < 3; l++) {
    py->lo[l] = lo[l];
    size[l] = hi[l] - lo[l];
    if (!(size[l] > 0.0)) return 1;
  }
  py->cell = cell;
  {
    int64_t d[3];
    int ax;
    for (ax = 0; ax < 3; ax++) {
      d[ax] = (int64_t)(size[ax] / cell) + 1;
      if (d[ax] > 0x7FFFFF) return 1; /* 2^23 клеток на ось — потолок сетки листа */
    }
    py->nx = (int32_t)d[0];
    py->ny = (int32_t)d[1];
    py->nz = (int32_t)d[2];
  }
  py->ncells = (int64_t)py->nx * py->ny * py->nz;
  if (py->ncells > (int64_t)1 << 60) return 1; /* контракт hz_pyr_wrap: 2n переполнения нет */
  py->nt = nt;

  /* 1. куски В ПОРЯДКЕ ТРЕУГОЛЬНИКОВ (А1456); индекс владения появится,
   * когда будут отсортированные занятые листья */
  py->pcs = (hz_piece *)malloc((size_t)nt * sizeof *py->pcs);
  if (!py->pcs) return 2;
  for (t = 0; t < nt; t++) {
    py->pcs[t].tri = t;
    py->pcs[t].mtl = mtl[t];
    py->pcs[t].e = 0.0f;
    py->pcs[t].cell = -1;
  }

  /* 2. консервативная занятость по bbox — битовая карта */
  marks = (int64_t *)calloc((size_t)((py->ncells + 63) >> 6), sizeof *marks);
  if (!marks) {
    hz_pyr_free(py);
    return 2;
  }
  for (t = 0; t < nt; t++) {
    int64_t x[2], y[2], z[2], ix, iy, iz;
    pyr_bbox_span(py, t, tri_min, tri_max, x, y, z);
    for (iz = z[0]; iz <= z[1]; iz++)
      for (iy = y[0]; iy <= y[1]; iy++)
        for (ix = x[0]; ix <= x[1]; ix++) {
          int64_t id = ix + (int64_t)py->nx * (iy + (int64_t)py->ny * iz);
          if ((marks[id >> 6] & ((int64_t)1 << (id & 63))) == 0) {
            marks[id >> 6] |= (int64_t)1 << (id & 63);
            py->nmarkleaf++;
          }
          py->nmarks++;
        }
  }

  /* 3. занятые листья — из карты по возрастанию id */
  py->leaf_id = (int64_t *)malloc((size_t)py->nmarkleaf * sizeof *py->leaf_id);
  py->leaf = (hz_pyr_leaf *)malloc((size_t)py->nmarkleaf * sizeof *py->leaf);
  if (!py->leaf_id || !py->leaf) {
    hz_pyr_free(py);
    return 2;
  }
  {
    int32_t n = 0;
    int64_t id;
    for (id = 0; id < py->ncells; id++)
      if (marks[id >> 6] & ((int64_t)1 << (id & 63))) {
        py->leaf_id[n] = id;
        py->leaf[n].pcs_first = 0;
        py->leaf[n].npcs = 0;
        py->leaf[n].state = HZ_PYR_FINE;
        n++;
      }
    py->nleaf = n;
  }
  free(marks);
  marks = NULL;

  /* 4. владение: центроид → позиция в отсортированном leaf_id */
  cnt = (int64_t *)calloc((size_t)py->nleaf + 1, sizeof *cnt);
  if (!cnt) {
    hz_pyr_free(py);
    return 2;
  }
  for (t = 0; t < nt; t++) {
    int64_t id = hz_pyr_leaf_index(py->lo, cell, py->nx, py->ny, py->nz, centroid + 3 * (int64_t)t);
    int32_t li = pyr_find(py->leaf_id, py->nleaf, id);
    if (li < 0) {
      /* центр вне bbox-занятости — дефект входа (должно быть невозможно:
       * bbox содержит центроид), не молчим */
      free(cnt);
      hz_pyr_free(py);
      return 1;
    }
    py->pcs[t].cell = li;
    cnt[li + 1]++;
  }

  /* 5. CSR: counting-sort по клеткам; внутри клетки — порядок треугольников.
   * cnt после шага 4 — ПО-КЛЕТОЧНЫЕ счётчики (cnt[l+1] = число кусков клетки
   * l); префикс-сумма обязана быть накоплена явно, иначе pcs_first получает
   * счётчик соседней клетки — именно это и ловил детектор (а) на чистом
   * прогоне §833. */
  for (l = 0; l < py->nleaf; l++)
    cnt[l + 1] += cnt[l];
  for (l = 0; l < py->nleaf; l++) {
    py->leaf[l].pcs_first = (int32_t)cnt[l];
    py->leaf[l].npcs = (int32_t)(cnt[l + 1] - cnt[l]);
    if (py->leaf[l].npcs > 0) py->ncentleaf++;
  }
  py->csr = (int32_t *)malloc((size_t)nt * sizeof *py->csr);
  fill = (int32_t *)malloc((size_t)(py->nleaf + 1) * sizeof *fill);
  if (!py->csr || !fill) {
    free(cnt);
    hz_pyr_free(py);
    return 2;
  }
  for (l = 0; l < py->nleaf; l++)
    fill[l] = py->leaf[l].pcs_first;
  for (t = 0; t < nt; t++)
    py->csr[fill[py->pcs[t].cell]++] = t;
  free(cnt);
  cnt = NULL;
  free(fill);
  fill = NULL;

  /* 6. предки ×8 до корня: на каждом уровне отсортированные уникальные id
   * непустых клеток; уровень 1×1×1 — корень, выше не идём */
  cur = (int64_t *)malloc((size_t)py->nleaf * sizeof *cur);
  if (!cur) {
    hz_pyr_free(py);
    return 2;
  }
  memcpy(cur, py->leaf_id, (size_t)py->nleaf * sizeof *cur);
  ncur = py->nleaf;
  pnx = py->nx;
  pny = py->ny;
  pnz = py->nz;
  while (!(pnx == 1 && pny == 1 && pnz == 1)) {
    int64_t cnx = (pnx + 1) >> 1, cny = (pny + 1) >> 1, cnz = (pnz + 1) >> 1;
    int32_t u, w;
    hz_pyr_node *lvl;
    int64_t *nc;
    par = (int64_t *)malloc((size_t)ncur * sizeof *par);
    if (!par) {
      free(cur);
      hz_pyr_free(py);
      return 2;
    }
    for (u = 0; u < ncur; u++) {
      int64_t id = cur[u];
      int64_t ix = id % pnx, iy = (id / pnx) % pny, iz = id / (pnx * pny);
      par[u] = (ix >> 1) + cnx * ((iy >> 1) + cny * (iz >> 1));
    }
    qsort(par, (size_t)ncur, sizeof *par, pyr_cmp_i64);
    w = 0;
    for (u = 0; u < ncur; u++)
      if (w == 0 || par[w - 1] != par[u]) par[w++] = par[u];
    lvl = (hz_pyr_node *)malloc((size_t)w * sizeof *lvl);
    if (!lvl) {
      free(par);
      free(cur);
      hz_pyr_free(py);
      return 2;
    }
    for (u = 0; u < w; u++) {
      lvl[u].id = par[u];
      lvl[u].state = HZ_PYR_FINE; /* агрегатов в §833 нет */
    }
    py->lev = (hz_pyr_node **)realloc(py->lev, (size_t)(py->nlev + 1) * sizeof *py->lev);
    py->nlev_nodes =
        (int32_t *)realloc(py->nlev_nodes, (size_t)(py->nlev + 1) * sizeof *py->nlev_nodes);
    if (!py->lev || !py->nlev_nodes) {
      free(lvl);
      free(par);
      free(cur);
      hz_pyr_free(py);
      return 2;
    }
    py->lev[py->nlev] = lvl;
    py->nlev_nodes[py->nlev] = w;
    py->nlev++;
    /* подъём: cur ← узлы этого уровня */
    nc = (int64_t *)realloc(cur, (size_t)w * sizeof *nc);
    if (!nc) {
      free(par);
      free(cur);
      hz_pyr_free(py);
      return 2;
    }
    cur = nc;
    for (u = 0; u < w; u++)
      cur[u] = par[u];
    ncur = w;
    free(par);
    pnx = cnx;
    pny = cny;
    pnz = cnz;
  }
  free(cur);

  return 0;
}

int32_t hz_pyr_screw(hz_pyr *py, int32_t stride) {
  int32_t t, n = 0;
  if (!py || !py->pcs || py->nt <= 0 || py->nleaf <= 0) return 0;
  for (t = 0; t < py->nt; t++) {
    int32_t hit;
    if (stride == 1)
      hit = (t == 0); /* scrfew: калибровка детектора одним куском */
    else
      hit = (t % stride) == 0;
    if (hit) {
      /* клетка +1 по ЛИНЕЙНОМУ индексу с заворотом по всей сетке (А1454).
       * CSR не трогаем — разладка и есть предмет детекторов; -1 законен:
       * «кусок в объявленно пустой клетке». */
      int64_t lin = py->leaf_id[py->pcs[t].cell];
      int64_t nxt = hz_pyr_wrap(lin + 1, py->ncells);
      py->pcs[t].cell = pyr_find(py->leaf_id, py->nleaf, nxt);
      n++;
    }
  }
  return n;
}

void hz_pyr_verify(const hz_pyr *py, int32_t nt, const double *tri_min, const double *tri_max,
                   const double *centroid, hz_pyr_verdict *vd) {
  int32_t t, li;
  memset(vd, 0, sizeof *vd);
  if (!py || !py->pcs || !py->leaf || !py->leaf_id) return;

  /* (а) владение: пересчёт по центроиду … */
  for (t = 0; t < nt && t < py->nt; t++) {
    int64_t id =
        hz_pyr_leaf_index(py->lo, py->cell, py->nx, py->ny, py->nz, centroid + 3 * (int64_t)t);
    int32_t want = pyr_find(py->leaf_id, py->nleaf, id);
    if (want != py->pcs[t].cell) vd->d_cell++;
  }
  /* (а) … и встречно: CSR клетки ссылается только на её куски */
  for (li = 0; li < py->nleaf; li++) {
    int32_t u;
    for (u = 0; u < py->leaf[li].npcs; u++) {
      int32_t p = py->csr[py->leaf[li].pcs_first + u];
      if (p < 0 || p >= py->nt || py->pcs[p].cell != li) vd->d_csr++;
    }
  }

  /* (б) независимый второй проход bbox-ов */
  for (t = 0; t < nt && t < py->nt; t++) {
    int64_t x[2], y[2], z[2], ix, iy, iz;
    pyr_bbox_span(py, t, tri_min, tri_max, x, y, z);
    for (iz = z[0]; iz <= z[1]; iz++)
      for (iy = y[0]; iy <= y[1]; iy++)
        for (ix = x[0]; ix <= x[1]; ix++) {
          int64_t id = ix + (int64_t)py->nx * (iy + (int64_t)py->ny * iz);
          if (pyr_find(py->leaf_id, py->nleaf, id) < 0) vd->d_bbox++;
        }
  }

  /* (в) куски в объявленно пустых клетках (cell == -1 после скрещивания) */
  for (t = 0; t < nt && t < py->nt; t++)
    if (py->pcs[t].cell < 0 || py->pcs[t].cell >= py->nleaf) vd->d_empty++;
}

void hz_pyr_free(hz_pyr *py) {
  int32_t l;
  if (!py) return;
  free(py->pcs);
  free(py->csr);
  free(py->leaf_id);
  free(py->leaf);
  for (l = 0; l < py->nlev; l++)
    free(py->lev[l]);
  free(py->lev);
  free(py->nlev_nodes);
  memset(py, 0, sizeof *py);
}
