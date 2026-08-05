/* pmark.c — пометки от камеры. Разбор — в `pmark.h`. */
#include "pmark.h"
#include <math.h>

/* Целиком ли коробка снаружи хотя бы одной плоскости пирамиды. Проверяется по
 * дальнему углу: если ДАЛЬНИЙ угол снаружи, снаружи и вся коробка. */
static int box_outside(const hz_pmark_cfg *C, const double *lo, const double *hi) {
  for (int k = 0; k < 6; k++) {
    const double *P = C->fr[k];
    double s = P[3];
    for (int c = 0; c < 3; c++)
      s += (P[c] > 0.0 ? hi[c] : lo[c]) * P[c];
    if (s < 0.0) return 1;
  }
  return 0;
}

/* Предельный уровень по камерному LOD: спускаться, пока угловой размер ячейки
 * больше пола. Ячейка на уровне `L` вдвое меньше, чем на `L−1`, поэтому число
 * шагов есть логарифм отношения «нынешний угловой размер / пол». */
static int lod_level(const hz_pmark_cfg *C, const double *lo, const double *hi, int lev) {
  double r2 = 0.0, d2 = 0.0;
  for (int c = 0; c < 3; c++) {
    double h = 0.5 * (hi[c] - lo[c]);
    double dd = 0.5 * (lo[c] + hi[c]) - C->eye[c];
    r2 += h * h;
    d2 += dd * dd;
  }
  double rad = sqrt(r2), R = sqrt(d2);
  double floor_size = C->px * C->eps * R; /* допустимый размер ячейки здесь */
  if (!(floor_size > 0.0)) return HZ_PMARK_NONE;
  double k = 2.0 * rad / floor_size; /* во сколько раз ячейка крупнее нужного */
  if (!(k > 1.0)) return lev;        /* уже достаточно грубо */
  int steps = (int)ceil(log2(k));
  int L = lev + steps;
  return (L >= HZ_PMARK_NONE) ? HZ_PMARK_NONE - 1 : L;
}

static void mark_walk(const hz_ptree *T, hz_pmark_cfg *C, unsigned char *mark, int32_t nid,
                      int lev) {
  const hz_ptnode *N = &T->nd[nid];
  C->nnode++;
  int L;
  if (box_outside(C, N->lo, N->hi)) {
    /* ВНЕ КОНУСА: пол отпускается на пару уровней, а не отменяется. */
    L = lev + C->outlevels;
    C->nout++;
  } else {
    L = lod_level(C, N->lo, N->hi, lev);
    C->nin++;
  }
  mark[nid] = (unsigned char)(L > HZ_PMARK_NONE - 1 ? HZ_PMARK_NONE - 1 : L);
  if (N->child < 0 || lev >= C->level) return;
  for (int k = 0; k < 8; k++)
    mark_walk(T, C, mark, N->child + k, lev + 1);
}

void hz_pmark_build(const hz_ptree *T, hz_pmark_cfg *C, unsigned char *mark) {
  for (int32_t i = 0; i < T->nnd; i++)
    mark[i] = HZ_PMARK_NONE;
  C->nnode = 0;
  C->nout = 0;
  C->nin = 0;
  mark_walk(T, C, mark, 0, 0);
}
