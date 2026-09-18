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
static int32_t pyr_find(const int64_t *a, int32_t n, int64_t x);

/* вперёд: размеры сетки уровня и поиск узла — нужны hz_pyr_set_lp (§852) */
static void pyr_level_dims(const hz_pyr *py, int32_t l, int64_t d[3]);
static int32_t pyr_find_node(const hz_pyr_node *a, int32_t n, int64_t x);

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
  int64_t *marks = NULL; /* пометки bbox-занятости, ВРЕМЕННАЯ (в бюджет не входит) */
  int64_t *cnt = NULL;
  int32_t *fill = NULL;
  int64_t *cur = NULL, *par = NULL;
  int64_t nm = 0;
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

  /* 2. консервативная занятость по bbox — МАССИВ ПОМЕТОК + сортировка
   * (§837): плотная битовая карта требует ncells/8 байт и не переживает
   * grid-сцены (габарит сетки комнат при мелкой клетке ~1e17 клеток);
   * сортировка пометок — секунды, память известна заранее. nm — число
   * пометок, живёт до конца секции 3. */
  {
    int64_t cap = 1 << 20;
    marks = (int64_t *)malloc((size_t)cap * sizeof *marks);
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
            if (nm == cap) {
              int64_t *nm2 = (int64_t *)realloc(marks, (size_t)(cap *= 2) * sizeof *marks);
              if (!nm2) {
                free(marks);
                hz_pyr_free(py);
                return 2;
              }
              marks = nm2;
            }
            marks[nm++] = id;
            py->nmarks++;
          }
    }
    qsort(marks, (size_t)nm, sizeof *marks, pyr_cmp_i64);
    py->nmarkleaf = 0;
    {
      int64_t u;
      for (u = 0; u < nm; u++)
        if (u == 0 || marks[u] != marks[u - 1]) py->nmarkleaf++;
    }
  }

  /* 3. занятые листья — из отсортированных уникальных пометок */
  py->leaf_id = (int64_t *)malloc((size_t)py->nmarkleaf * sizeof *py->leaf_id);
  py->leaf = (hz_pyr_leaf *)malloc((size_t)py->nmarkleaf * sizeof *py->leaf);
  if (!py->leaf_id || !py->leaf) {
    hz_pyr_free(py);
    return 2;
  }
  {
    int32_t n = 0;
    int64_t u;
    for (u = 0; u < nm; u++)
      if (u == 0 || marks[u] != marks[u - 1]) {
        py->leaf_id[n] = marks[u];
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

  /* (а) владение: пересчёт по центроиду ИСХОДНОГО треугольника куска
   * (pcs[t].tri — после Morton слот ≠ треугольник, А1491) … */
  for (t = 0; t < nt && t < py->nt; t++) {
    int32_t tri = py->pcs[t].tri;
    int64_t id =
        hz_pyr_leaf_index(py->lo, py->cell, py->nx, py->ny, py->nz, centroid + 3 * (int64_t)tri);
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

  /* (б) независимый второй проход bbox-ов (по исходному треугольнику) */
  for (t = 0; t < nt && t < py->nt; t++) {
    int64_t x[2], y[2], z[2], ix, iy, iz;
    int32_t tri = py->pcs[t].tri;
    pyr_bbox_span(py, tri, tri_min, tri_max, x, y, z);
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

/* §852 (А1567): per-node max ℓ_p одним подъёмом. lp — [nt] в порядке
 * ИСХОДНЫХ треугольников (кусок p смотрит треугольник py->pcs[p].tri).
 * Куски с ℓ_p выше самого крупного уровня ОБРАБАТЫВАЮТСЯ НА ЛИСТЕ
 * (ℓ_p := 0) — агрегатов нет (А1568), но счётчик обязателен: молчаливый
 * клэмп вверх = скрытая агрегация. Вызывать ПОСЛЕ hz_pyr_morton. */
void hz_pyr_clear_lp(hz_pyr *py) {
  int32_t l;
  if (!py) return;
  free(py->leaf_lp); /* NULL безопасен; повторный set_lp обязан сбросить */
  py->leaf_lp = NULL;
  if (py->lev_lp) {
    for (l = 0; l < py->nlev; l++)
      free(py->lev_lp[l]);
    free(py->lev_lp);
    py->lev_lp = NULL;
  }
}

int hz_pyr_set_lp(hz_pyr *py, const uint8_t *lp, int32_t *nup) {
  int32_t l, li, u, up = 0;
  hz_pyr_clear_lp(py); /* §862: повторный вызов (адаптивные этажи) валиден */
  if (!py || !lp || py->nt <= 0 || py->nleaf <= 0) return 1;
  if (py->nlev > 254) return 1; /* ℓ_p живёт в байте; уровней столько не бывает */

  /* лист: max ℓ_p кусков клетки; недостижимые уровни — на лист и в счётчик */
  py->leaf_lp = (uint8_t *)malloc((size_t)py->nleaf);
  if (!py->leaf_lp) return 2;
  for (li = 0; li < py->nleaf; li++)
    py->leaf_lp[li] = 255; /* нет кусков — не материален ни на каком уровне */
  for (u = 0; u < py->nt; u++) {
    uint8_t v = lp[py->pcs[u].tri]; /* ℓ_p исходного треугольника этого слота */
    if (v > (uint8_t)py->nlev) {    /* клэмп ВВЕРХ только как явный отказ-на-лист (А1568) */
      v = 0;
      up++;
    }
    li = py->pcs[u].cell;
    if (li < 0 || li >= py->nleaf) return 1;
    if (py->leaf_lp[li] == 255 || v > py->leaf_lp[li]) py->leaf_lp[li] = v;
  }
  /* подъём: узел уровня ℓ = max по существующим детям; пропуски (ПУСТ)
   * вклада не дают — состояние узла производное от поддерева */
  py->lev_lp = (uint8_t **)calloc((size_t)py->nlev, sizeof *py->lev_lp);
  if (!py->lev_lp) return 2;
  for (l = 0; l < py->nlev; l++) {
    int64_t pd[3], cd[3];
    uint8_t *arr = (uint8_t *)malloc((size_t)py->nlev_nodes[l]);
    if (!arr) return 2;
    py->lev_lp[l] = arr;
    if (l == 0) {
      cd[0] = py->nx;
      cd[1] = py->ny;
      cd[2] = py->nz;
    } else
      pyr_level_dims(py, l - 1, cd);
    pyr_level_dims(py, l, pd);
    for (u = 0; u < py->nlev_nodes[l]; u++) {
      int64_t kid = py->lev[l][u].id;
      int64_t ci[3];
      uint8_t m = 255;
      int cx, cy, cz;
      ci[0] = kid % pd[0];
      ci[1] = (kid / pd[0]) % pd[1];
      ci[2] = kid / (pd[0] * pd[1]);
      for (cz = 0; cz < 2; cz++)
        for (cy = 0; cy < 2; cy++)
          for (cx = 0; cx < 2; cx++) {
            int64_t id = 2 * ci[0] + cx + cd[0] * (2 * ci[1] + cy + cd[1] * (2 * ci[2] + cz));
            int32_t found;
            uint8_t v;
            if (2 * ci[0] + cx >= cd[0] || 2 * ci[1] + cy >= cd[1] || 2 * ci[2] + cz >= cd[2])
              continue; /* нечётная сетка: ребёнка за границей нет (А1567) */
            if (l == 0) {
              found = pyr_find(py->leaf_id, py->nleaf, id);
              if (found < 0) continue; /* ПУСТ-ребёнок: вклада нет */
              v = py->leaf_lp[found];
            } else {
              found = pyr_find_node(py->lev[l - 1], py->nlev_nodes[l - 1], id);
              if (found < 0) continue;
              v = py->lev_lp[l - 1][found];
            }
            if (v != 255 && (m == 255 || v > m)) m = v;
          }
      arr[u] = m;
    }
  }
  if (nup) *nup = up;
  return 0;
}

void hz_pyr_free(hz_pyr *py) {
  int32_t l;
  if (!py) return;
  free(py->pcs);
  free(py->csr);
  free(py->leaf_id);
  free(py->leaf);
  free(py->leaf_lp);
  if (py->lev_lp) {
    for (l = 0; l < py->nlev; l++)
      free(py->lev_lp[l]);
    free(py->lev_lp);
  }
  free(py->perm);
  for (l = 0; l < py->nlev; l++)
    free(py->lev[l]);
  free(py->lev);
  free(py->nlev_nodes);
  memset(py, 0, sizeof *py);
}

/* ---- МАРШ (§834) --------------------------------------------------------- */

/* размеры сетки уровня ℓ (0 — родители листов); сторона клетки — cell·2^(ℓ+1) */
static void pyr_level_dims(const hz_pyr *py, int32_t l, int64_t d[3]) {
  int32_t sh = (int32_t)(l + 1);
  d[0] = ((int64_t)py->nx + ((int64_t)1 << sh) - 1) >> sh;
  d[1] = ((int64_t)py->ny + ((int64_t)1 << sh) - 1) >> sh;
  d[2] = ((int64_t)py->nz + ((int64_t)1 << sh) - 1) >> sh;
}

/* проекция центра клетки на omega по СДВИГУ уровня: лист — sh=0 (базовая
 * сетка), внутренний уровень ℓ — sh=ℓ+1 (А1466-ряд: сторона cell·2^sh) */
static double pyr_slab(const hz_pyr *py, int32_t sh, int64_t id, const double omega[3],
                       int invert_x) {
  int64_t d[3];
  double side, c[3];
  int ax;
  d[0] = ((int64_t)py->nx + ((int64_t)1 << sh) - 1) >> sh;
  d[1] = ((int64_t)py->ny + ((int64_t)1 << sh) - 1) >> sh;
  d[2] = ((int64_t)py->nz + ((int64_t)1 << sh) - 1) >> sh;
  side = py->cell * (double)((int64_t)1 << sh);
  c[0] = (double)(id % d[0]) + 0.5;
  c[1] = (double)((id / d[0]) % d[1]) + 0.5;
  c[2] = (double)(id / (d[0] * d[1])) + 0.5;
  for (ax = 0; ax < 3; ax++)
    c[ax] = py->lo[ax] + c[ax] * side;
  return omega[0] * c[0] * (invert_x ? -1.0 : 1.0) + omega[1] * c[1] + omega[2] * c[2];
}

/* двоичный поиск id в отсортированном массиве узлов уровня; -1 если нет */
static int32_t pyr_find_node(const hz_pyr_node *a, int32_t n, int64_t x) {
  int32_t lo = 0, hi = n;
  while (lo < hi) {
    int32_t mid = lo + (hi - lo) / 2;
    if (a[mid].id < x)
      lo = mid + 1;
    else
      hi = mid;
  }
  if (lo < n && a[lo].id == x) return lo;
  return -1;
}

typedef struct {
  const hz_pyr *py;
  const double *omega;
  int invert_x;
  hz_pyr_march_stat *st;
  uint64_t *bits;
  double last_s;
  int have_last;
  int32_t *vindex;
} pyr_march_ctx;

static void pyr_visit(pyr_march_ctx *mc, int32_t l, int32_t pos) {
  const hz_pyr *py = mc->py;
  if (l < 0) { /* лист (уровень −1): sh = 0, базовая сетка */
    double s = pyr_slab(py, 0, py->leaf_id[pos], mc->omega, mc->invert_x);
    if (mc->have_last && s < mc->last_s) mc->st->inversions++;
    mc->last_s = s;
    mc->have_last = 1;
    if (mc->bits) mc->bits[pos >> 6] |= (uint64_t)1 << (pos & 63);
    mc->st->leaves++;
    mc->st->nodes++;
    if (mc->vindex) mc->vindex[pos] = (int32_t)mc->st->leaves; /* номер посещения, с 1 */
    return;
  }
  mc->st->nodes++;
  if (py->lev[l][pos].state == HZ_PYR_AGGR) { /* стоп: одно событие */
    mc->st->aggr++;
    return;
  }
  {
    int64_t pd[3], cd[3], ci[3], kid;
    int32_t cl = l - 1, cx, cy, cz; /* вниз к листьям: уровень 0 — предки листов */
    double s[8];
    int32_t pos_ch[8];
    int32_t n = 0, u, v;
    pyr_level_dims(py, l, pd);
    if (cl < 0) { /* листья: базовая сетка */
      cd[0] = py->nx;
      cd[1] = py->ny;
      cd[2] = py->nz;
    } else
      pyr_level_dims(py, cl, cd);
    kid = py->lev[l][pos].id;
    ci[0] = kid % pd[0];
    ci[1] = (kid / pd[0]) % pd[1];
    ci[2] = kid / (pd[0] * pd[1]);
    for (cz = 0; cz < 2; cz++)
      for (cy = 0; cy < 2; cy++)
        for (cx = 0; cx < 2; cx++) {
          int64_t qx = 2 * ci[0] + cx, qy = 2 * ci[1] + cy, qz = 2 * ci[2] + cz;
          int32_t found;
          if (qx >= cd[0] || qy >= cd[1] || qz >= cd[2]) continue;
          kid = qx + cd[0] * (qy + cd[1] * qz);
          if (cl < 0)
            found = pyr_find(py->leaf_id, py->nleaf, kid);
          else
            found = pyr_find_node(py->lev[cl], py->nlev_nodes[cl], kid);
          if (found < 0) continue; /* ПУСТ: прыжок, клетка не хранится */
          s[n] = pyr_slab(py, cl < 0 ? 0 : cl + 1, kid, mc->omega, mc->invert_x);
          pos_ch[n] = found; /* позиция в массиве своего уровня */
          n++;
        }
    /* сортировка вставками по s (n ≤ 8) */
    for (u = 1; u < n; u++) {
      double su = s[u];
      int32_t pu = pos_ch[u];
      for (v = u; v > 0 && s[v - 1] > su; v--) {
        s[v] = s[v - 1];
        pos_ch[v] = pos_ch[v - 1];
      }
      s[v] = su;
      pos_ch[v] = pu;
    }
    for (u = 0; u < n; u++)
      pyr_visit(mc, cl, pos_ch[u]);
  }
}

void hz_pyr_march(const hz_pyr *py, const double omega[3], int invert_x, hz_pyr_march_stat *st,
                  uint64_t *bits, int32_t *vindex) {
  pyr_march_ctx mc;
  if (!py || !st) {
    if (st) memset(st, 0, sizeof *st);
    return;
  }
  memset(st, 0, sizeof *st);
  if (bits) memset(bits, 0, (size_t)((py->nleaf + 63) >> 6) * sizeof *bits);
  if (vindex) memset(vindex, -1, (size_t)py->nleaf * sizeof *vindex);
  mc.py = py;
  mc.omega = omega;
  mc.invert_x = invert_x;
  mc.st = st;
  mc.vindex = vindex;
  mc.bits = bits;
  mc.last_s = 0.0;
  mc.have_last = 0;
  if (py->nleaf <= 0) return;
  if (py->nlev == 0) { /* нет внутренних уровней: занятые листья — сами по себе */
    pyr_visit(&mc, -1, 0);
    return;
  }
  pyr_visit(&mc, py->nlev - 1, 0);
}

void hz_pyr_mark_aggr(hz_pyr *py, int32_t every) {
  int32_t l, u;
  int64_t g = 0;
  if (!py || every <= 0) return;
  if (every == 1) { /* только корень (А1467) */
    py->lev[py->nlev - 1][0].state = HZ_PYR_AGGR;
    return;
  }
  for (l = 0; l < py->nlev; l++)
    for (u = 0; u < py->nlev_nodes[l]; u++)
      if (g++ % every == 0) py->lev[l][u].state = HZ_PYR_AGGR;
}

void hz_pyr_unmark_aggr(hz_pyr *py) {
  int32_t l, u;
  if (!py) return;
  for (l = 0; l < py->nlev; l++)
    for (u = 0; u < py->nlev_nodes[l]; u++)
      if (py->lev[l][u].state == HZ_PYR_AGGR) py->lev[l][u].state = HZ_PYR_FINE;
}

/* ---- §839: MORTON-РАСКЛАДКА КУСКОВ --------------------------------------- */

/* разнести 21 бит координаты через по 2 бита (63-битный код) */
static uint64_t pyr_spread(uint64_t x) {
  uint64_t r = x & 0x1FFFFF;
  r = (r | (r << 32)) & 0x7F00000000FFFF;
  r = (r | (r << 16)) & 0x7FFF0000FF0000FF;
  r = (r | (r << 8)) & 0x7F00FF00FF00FF00;
  r = (r | (r << 4)) & 0x70E38E38E38E38E3;
  r = (r | (r << 2)) & 0x1249249249249249;
  return r;
}

static uint64_t pyr_morton3(int64_t i, int64_t j, int64_t k) {
  return pyr_spread((uint64_t)i) | (pyr_spread((uint64_t)j) << 1) | (pyr_spread((uint64_t)k) << 2);
}

typedef struct {
  uint64_t key;
  int32_t p;
} pyr_mpair;

static int pyr_cmp_mpair(const void *a, const void *b) {
  uint64_t x = ((const pyr_mpair *)a)->key, y = ((const pyr_mpair *)b)->key;
  return (x > y) - (x < y);
}

int hz_pyr_morton(hz_pyr *py) {
  pyr_mpair *mp = NULL;
  hz_piece *tmp = NULL;
  int64_t *cnt = NULL;
  int32_t *fill = NULL;
  int32_t p, l;
  if (!py || !py->pcs || !py->leaf_id || py->nt <= 0) return 1;
  mp = (pyr_mpair *)malloc((size_t)py->nt * sizeof *mp);
  tmp = (hz_piece *)malloc((size_t)py->nt * sizeof *tmp);
  if (!mp || !tmp) {
    free(mp);
    free(tmp);
    return 2;
  }
  for (p = 0; p < py->nt; p++) {
    int64_t id = py->leaf_id[py->pcs[p].cell];
    mp[p].key = pyr_morton3(id % py->nx, (id / py->nx) % py->ny, id / ((int64_t)py->nx * py->ny));
    mp[p].p = p;
  }
  qsort(mp, (size_t)py->nt, sizeof *mp, pyr_cmp_mpair);
  for (p = 0; p < py->nt; p++)
    tmp[p] = py->pcs[mp[p].p];
  memcpy(py->pcs, tmp, (size_t)py->nt * sizeof *tmp);
  py->perm = (int32_t *)malloc((size_t)py->nt * sizeof *py->perm);
  if (!py->perm) {
    free(tmp);
    free(mp);
    hz_pyr_free(py);
    return 2;
  }
  for (p = 0; p < py->nt; p++)
    py->perm[p] = mp[p].p; /* слот → исходный кусок */
  free(tmp);
  tmp = NULL;
  free(mp);
  mp = NULL;

  /* CSR заново: куски сгруппированы по клеткам в новом порядке Morton */
  cnt = (int64_t *)calloc((size_t)py->nleaf + 1, sizeof *cnt);
  if (!cnt) {
    hz_pyr_free(py);
    return 2;
  }
  for (p = 0; p < py->nt; p++)
    cnt[py->pcs[p].cell + 1]++;
  for (l = 0; l < py->nleaf; l++)
    cnt[l + 1] += cnt[l];
  for (l = 0; l < py->nleaf; l++) {
    py->leaf[l].pcs_first = (int32_t)cnt[l];
    py->leaf[l].npcs = (int32_t)(cnt[l + 1] - cnt[l]);
  }
  free(py->csr); /* Morton перестраивает CSR: старый массив освобождается (LeakSanitizer §839) */
  py->csr = (int32_t *)malloc((size_t)py->nt * sizeof *py->csr);
  fill = (int32_t *)malloc((size_t)(py->nleaf + 1) * sizeof *fill);
  if (!py->csr || !fill) {
    free(cnt);
    free(fill); /* fill мог выделиться при отказе соседнего malloc */
    hz_pyr_free(py);
    return 2;
  }
  for (l = 0; l < py->nleaf; l++)
    fill[l] = py->leaf[l].pcs_first;
  for (p = 0; p < py->nt; p++)
    py->csr[fill[py->pcs[p].cell]++] = p;
  free(cnt);
  free(fill);
  return 0;
}

/* §839: применить перестановку Morton к массиву потребителя — циклами на
 * месте, без временного массива (gcc-analyzer теряет связь
 * «проверка↔выделение» на временных массивах, diam 07-24). */
int hz_pyr_permute(hz_pyr *py, void *base, size_t elem) {
  uint8_t *seen = NULL;
  int32_t i;
  char *buf, *b = (char *)base;
  if (!py || !py->perm || !base || elem == 0 || elem > 256) return 1;
  seen = (uint8_t *)calloc((size_t)py->nt, 1);
  buf = (char *)malloc(elem);
  if (!seen || !buf) {
    free(seen);
    free(buf);
    return 2;
  }
  for (i = 0; i < py->nt; i++) {
    int32_t j, k;
    if (seen[i]) continue;
    memcpy(buf, b + (size_t)i * elem, elem); /* старое значение начала цикла */
    j = i;
    while (1) {
      k = py->perm[j]; /* слот j получает старое значение из позиции perm[j] */
      seen[j] = 1;     /* ПОСЛЕДНИЙ слот цикла тоже: иначе фантомный второй проход
                        * применит перестановку дважды (ловлено §839 на полости) */
      if (k == i) {
        memcpy(b + (size_t)j * elem, buf, elem);
        break;
      }
      memcpy(b + (size_t)j * elem, b + (size_t)k * elem, elem);
      j = k;
    }
  }
  free(seen);
  free(buf);
  return 0;
}

/* §849: публичный поиск листа по линейному id клетки (для DDA-сбора
 * камеры); -1 если клетка пуста. Собственно тот же двоичный поиск, что
 * внутри марша. */
int32_t hz_pyr_leaf_pos(const hz_pyr *py, int64_t id) {
  if (!py || !py->leaf_id || py->nleaf <= 0) return -1;
  return pyr_find(py->leaf_id, py->nleaf, id);
}

int32_t hz_pyr_node_pos(const hz_pyr *py, int32_t l, int64_t id) {
  if (!py || l < 0 || l >= py->nlev) return -1;
  return pyr_find_node(py->lev[l], py->nlev_nodes[l], id);
}

void hz_pyr_level_dims(const hz_pyr *py, int32_t l, int64_t d[3]) {
  pyr_level_dims(py, l, d);
}
