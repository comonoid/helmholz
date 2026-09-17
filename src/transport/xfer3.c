/* xfer3.c — ЭТАП A, ШАГ 3: перенос хранимого DG1 между LOD-сетками. */

/* ЛОКАЛЬНЫЕ КООРДИНАТЫ И МАСШТАБ НАКЛОНА. Хранимое ячейки есть
 *   φ(ξ) = c + Σ_a g_a ξ_a,  ξ_a ∈ [−1/2, 1/2] в ЕДИНИЦАХ РАЗМЕРА ЯЧЕЙКИ,
 * то есть φ в углу = c ± g/2 — тот же договор, что читает сбор и печатает К94.
 * При смене размера ячейки h → h' безразмерная координата пересчитывается,
 * ξ' = (w − cen')/h', поэтому на общем поле g' = g·(h'/h) — наклон МЕЛЬЧЕЙ
 * ячейки в её координатах БОЛЬШЕ. (Первая редакция имела множитель перевёрнутым
 * — поймано предсказанием П2 на линейном поле до всякого прогона сцены.) Рама у обоих деревьев
 * одна, мировые множители в отношениях размеров сокращаются, и в формулах их нет. Объём — в кубе
 * единиц сетки: мировые множители стоят в СРЕДНЕЙ и в СИСТЕМЕ МОМЕНТОВ общими множителями и тоже
 * сокращаются. */

#include "transport/xfer3.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* 3×3 решение по Крамеру. Система моментов симметрична и положительно
 * определена (Σ V_i s_ia²/12 в диагонали), det строго положителен; ветка
 * det == 0 гасит наклон в ноль, оставляя ВЕРНУЮ среднюю, — отказ в закрытую
 * сторону, а не мусор. */
static void solve3(const double G[3][3], const double b[3], double g[3]) {
  double det = G[0][0] * (G[1][1] * G[2][2] - G[1][2] * G[2][1]) -
               G[0][1] * (G[1][0] * G[2][2] - G[1][2] * G[2][0]) +
               G[0][2] * (G[1][0] * G[2][1] - G[1][1] * G[2][0]);
  if (!(det > 0.0 || det < 0.0)) { /* без ==/!=: ноль здесь точный отказ, не допуск */
    g[0] = g[1] = g[2] = 0.0;
    return;
  }
  for (int k = 0; k < 3; k++) {
    double m[3][3];
    for (int r = 0; r < 3; r++)
      for (int c = 0; c < 3; c++)
        m[r][c] = c == k ? b[r] : G[r][c];
    double d = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
               m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
               m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
    g[k] = d / det;
  }
}

/* ПРОЛОНГАЦИЯ: dst-ячейка (центр cend, размер hs) лежит внутри src-листа
 * (центр csen, размер hp, хранимое cs/gs). Тождество полинома:
 * свободный член переносится вычислением φ в центре потомка; наклон
 * пересчитывается отношением размеров g' = g·(hs/hp) — из ξ = ξ'·(hs/hp)
 * (ξ безразмерен в единицах СВОЕЙ ячейки). Все размеры — в единицах сетки. */
static void prolong(const double csen[3], double hp, const double cs, const double gs[3],
                    const double cend[3], double hs, double *cd, double gd[3]) {
  double c2 = cs;
  for (int a = 0; a < 3; a++) {
    double xi0 = (cend[a] - csen[a]) / hp;
    c2 += gs[a] * xi0;
    gd[a] = gs[a] * (hs / hp);
  }
  *cd = c2;
}

/* Аккумулятор ограничения: накапливает листья src под одной dst-ячейкой.
 * G/b — точные моменты: ∫_i ξ_a φ dV и ∫_i ξ_a ξ_b dV по листу i в локальных
 * координатах dst; ∫ξ_i = 0, ∫ξ_i²/V = 1/12 — те же V/12, что в матрицах масс. */
typedef struct {
  int n;
  double vol, m0;
  double G[3][3], b[3];
} agg;

static void agg_add(agg *A, double V, double c, const double g[3], const double off[3],
                    const double s[3]) {
  /* ∫_i ξ_a φ dV = V·(c·o_a + g_a·s_a/12): перекрестных членов НЕТ, потому что
   * ∫ξ_i внутри листа i — ровно ноль. Первая редакция их добавляла — поймано
   * предсказанием П3b (линейное поле не восстановалось). */
  A->n++;
  A->vol += V;
  A->m0 += V * c;
  for (int a = 0; a < 3; a++) {
    A->b[a] += V * (c * off[a] + g[a] * s[a] / 12.0);
    for (int bb = 0; bb < 3; bb++) {
      double q = off[a] * off[bb] + (a == bb ? s[a] * s[a] / 12.0 : 0.0);
      A->G[a][bb] += V * q;
    }
  }
}

/* Обход src-поддерева, ЦЕЛИКОМ лежащего в dst-коробке: каждый лист — в
 * аккумулятор. Порядок обхода детерминирован (deterministic DFS). */
static int agg_walk(const tr3_mesh *src, const double *ps, int32_t ni, const int32_t lo[3],
                    int32_t size, const double cdst[3], double hsd, agg *A) {
  if (src->tree->nodes[ni].child0 < 0) {
    if (src->cellof[ni] < 0)
      return -1; /* А1/А2: контракт нарушен — отказ, не чтение мимо массива */
    int32_t ci = src->cellof[ni];
    double V = (double)size * (double)size * (double)size;
    double off[3], s[3];
    for (int a = 0; a < 3; a++) {
      double cen = (double)lo[a] + 0.5 * (double)size;
      off[a] = (cen - cdst[a]) / hsd;
      s[a] = (double)size / hsd;
    }
    agg_add(A, V, ps[4 * ci], &ps[4 * ci + 1], off, s);
    return 0;
  }
  int32_t h = size / 2, k0 = src->tree->nodes[ni].child0;
  for (int k = 0; k < 8; k++) {
    int32_t clo[3] = {lo[0] + ((k & 1) ? h : 0), lo[1] + ((k & 2) ? h : 0),
                      lo[2] + ((k & 4) ? h : 0)};
    if (agg_walk(src, ps, k0 + k, clo, h, cdst, hsd, A)) return -1;
  }
  return 0;
}

/* Спуск по src из корня к dst-коробке (деревья выровнены, путь единственный).
 * Возвращает 0 при успехе, −1 при структурной ошибке. */
static int descend(const tr3_mesh *dst, double *pd, const tr3_mesh *src, const double *ps,
                   int32_t dst_ci, const int32_t blo[3], int32_t bsize, int32_t ni, int32_t nsize,
                   const int32_t nlo[3]) {
  if (nsize == bsize) {
    if (src->tree->nodes[ni].child0 < 0) {
      /* EQUAL: та же коробка — та же ячейка; перенос есть тождество */
      int32_t sci = src->cellof[ni];
      if (sci < 0) return -1; /* А2: контракт mesh↔tree нарушен — отказ (fail closed) */
      for (int i = 0; i < 4; i++)
        pd[4 * dst_ci + i] = ps[4 * sci + i]; /* ПОБИТОВАЯ копия */
      return 0;
    }
    /* AGGREGATE: src мельче — средняя по объёму с точными моментами */
    agg A;
    memset(&A, 0, sizeof A);
    double cdst[3];
    for (int a = 0; a < 3; a++)
      cdst[a] = (double)blo[a] + 0.5 * (double)bsize;
    int32_t k0 = src->tree->nodes[ni].child0;
    int32_t h = nsize / 2;
    for (int k = 0; k < 8; k++) {
      int32_t clo[3] = {nlo[0] + ((k & 1) ? h : 0), nlo[1] + ((k & 2) ? h : 0),
                        nlo[2] + ((k & 4) ? h : 0)};
      if (agg_walk(src, ps, k0 + k, clo, h, cdst, (double)bsize, &A)) return -1;
    }
    if (A.n <= 0 || !(A.vol > 0.0)) return -1;
    double g[3];
    solve3(A.G, A.b, g);
    pd[4 * dst_ci] = A.m0 / A.vol;
    for (int a = 0; a < 3; a++)
      pd[4 * dst_ci + 1 + a] = g[a];
    return 0;
  }
  /* nsize > bsize: dst-коробка строго внутри src-узла */
  if (src->tree->nodes[ni].child0 < 0) {
    /* PROLONG: src-лист содержит dst-ячейку целиком */
    int32_t sci = src->cellof[ni];
    if (sci < 0) return -1; /* А2: отказ, не угадывание */
    double cend[3], csen[3];
    for (int a = 0; a < 3; a++) {
      cend[a] = (double)blo[a] + 0.5 * (double)bsize;
      csen[a] = (double)nlo[a] + 0.5 * (double)nsize;
    }
    double cd = 0.0, gd[3] = {0, 0, 0};
    prolong(csen, (double)nsize, ps[4 * sci], &ps[4 * sci + 1], cend, (double)bsize, &cd, gd);
    pd[4 * dst_ci] = cd;
    for (int a = 0; a < 3; a++)
      pd[4 * dst_ci + 1 + a] = gd[a];
    return 0;
  }
  int32_t h = nsize / 2;
  int k = (blo[0] >= nlo[0] + h ? 1 : 0) | (blo[1] >= nlo[1] + h ? 2 : 0) |
          (blo[2] >= nlo[2] + h ? 4 : 0);
  int32_t clo[3] = {nlo[0] + ((k & 1) ? h : 0), nlo[1] + ((k & 2) ? h : 0),
                    nlo[2] + ((k & 4) ? h : 0)};
  return descend(dst, pd, src, ps, dst_ci, blo, bsize, src->tree->nodes[ni].child0 + k, h, clo);
}

int tr3_xfer(const tr3_mesh *dst, double *phi_dst, const tr3_mesh *src, const double *phi_src) {
  if (dst->tree == NULL || src->tree == NULL) return -1;
  if (dst->tree->log2size != src->tree->log2size) return -1;
  /* А4: формулы используют ТОЛЬКО размеры в единицах — это честно при ОДНОЙ
   * раме (мировые множители сокращаются). Разные рамы молча дали бы мусор.
   * Сравнение ПОБИТОВОЕ и сознательное: рама обязана быть той же структурой. */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
  for (int a = 0; a < 3; a++)
    if (dst->fr.o[a] != src->fr.o[a] || dst->fr.u[a] != src->fr.u[a]) return -1;
#pragma GCC diagnostic pop
  int32_t size = 1 << src->tree->log2size;
  int32_t lo[3] = {0, 0, 0};
  for (int32_t c = 0; c < dst->ncell; c++) {
    int32_t blo[3] = {dst->clo[c][0], dst->clo[c][1], dst->clo[c][2]};
    if (descend(dst, phi_dst, src, phi_src, c, blo, dst->csize[c], 0, size, lo) != 0) return -1;
  }
  return 0;
}
