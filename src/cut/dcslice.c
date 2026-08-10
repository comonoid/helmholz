/* dcslice.c — СРЕЗ. Разбор — в `dcslice.h`. */
#include "cut/dcslice.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

int hz_slice_init(hz_dcslice *s, int lev) {
  s->cap = 1024;
  s->n = 0;
  s->lev = lev;
  s->vbits = HZ_SLICE_VBITS;
  s->c = malloc((size_t)s->cap * sizeof *s->c);
  return s->c == NULL ? HZ_DC_ENOMEM : HZ_DC_OK;
}

void hz_slice_free(hz_dcslice *s) {
  free(s->c);
  s->c = NULL;
  s->n = s->cap = 0;
}

/* --- октаэдрическая нормаль ------------------------------------------------ */

/* Ошибка этого кодирования НЕ ПРЕДПОЛАГАЕТСЯ, а меряется прогоном: 8 бит на ось
 * дают угол порядка градуса, и годится это или нет — вопрос замера, а не вкуса.
 * Реализация обычная: проекция на октаэдр, нижняя полусфера отражается наружу. */
static double oct_sign(double x) {
  return x >= 0.0 ? 1.0 : -1.0;
}

uint16_t hz_oct_encode(const double n[3]) {
  double s = fabs(n[0]) + fabs(n[1]) + fabs(n[2]);
  if (!(s > 0.0)) return 0;
  double px = n[0] / s, py = n[1] / s;
  if (n[2] < 0.0) {
    double tx = (1.0 - fabs(py)) * oct_sign(px);
    double ty = (1.0 - fabs(px)) * oct_sign(py);
    px = tx;
    py = ty;
  }
  /* [-1, 1] -> [0, 255]. Округление к ближайшему, без смещения. */
  double ux = floor((px * 0.5 + 0.5) * 255.0 + 0.5);
  double uy = floor((py * 0.5 + 0.5) * 255.0 + 0.5);
  if (ux < 0.0) ux = 0.0;
  if (ux > 255.0) ux = 255.0;
  if (uy < 0.0) uy = 0.0;
  if (uy > 255.0) uy = 255.0;
  return (uint16_t)(((uint32_t)ux << 8) | (uint32_t)uy);
}

void hz_oct_decode(uint16_t o, double n[3]) {
  double px = (double)((o >> 8) & 0xFFu) / 255.0 * 2.0 - 1.0;
  double py = (double)(o & 0xFFu) / 255.0 * 2.0 - 1.0;
  double pz = 1.0 - fabs(px) - fabs(py);
  if (pz < 0.0) {
    double tx = (1.0 - fabs(py)) * oct_sign(px);
    double ty = (1.0 - fabs(px)) * oct_sign(py);
    px = tx;
    py = ty;
  }
  double m = sqrt(px * px + py * py + pz * pz);
  if (!(m > 0.0)) m = 1.0;
  n[0] = px / m;
  n[1] = py / m;
  n[2] = pz / m;
}

/* --- сборка ---------------------------------------------------------------- */

typedef struct {
  hz_dcslice *s;
  const hz_dctree *t;
  const hz_htab *ht;
  hz_dc_stop stop;
  void *sctx;
  int rc;
} sbctx;

static int sb_push(sbctx *B, int32_t ni, const int32_t lo[3], int32_t size, int lvl) {
  hz_dcslice *s = B->s;
  if (s->n >= s->cap) {
    int32_t nc = s->cap * 2;
    hz_dccell *nn = realloc(s->c, (size_t)nc * sizeof *nn);
    if (nn == NULL) return HZ_DC_ENOMEM;
    s->c = nn;
    s->cap = nc;
  }
  hz_dccell *cl = &s->c[s->n];
  memset(cl, 0, sizeof *cl);
  for (int a = 0; a < 3; a++)
    cl->lo[a] = (uint16_t)lo[a];
  cl->lvl = (uint8_t)lvl;

  /* ВЕРШИНА в долях ячейки. Хранится ПРОИЗВОДНОЕ, и потому квантуется; сколько
   * бит значимо — параметр среза, а не константа в формуле (негативный
   * контроль ставит 1). */
  double q = (double)((1u << B->s->vbits) - 1u);
  for (int a = 0; a < 3; a++) {
    double f = (hz_dc_vx(B->t, ni)[a] - (double)lo[a]) / (double)size;
    if (f < 0.0) f = 0.0;
    if (f > 1.0) f = 1.0;
    double u = floor(f * q + 0.5);
    /* Разряды выравниваются влево, чтобы декодирование было одно на все vbits. */
    unsigned v = (unsigned)u << (8 - B->s->vbits);
    cl->vx[a] = (uint8_t)(v > 255u ? 255u : v);
  }

  /* НОРМАЛЬ ЯЧЕЙКИ — сумма нормалей её пересечённых рёбер. Своего поля у узла
   * нет и заводить его незачем: величина производная, и в срез она попадает
   * ровно затем, чтобы фронт не ходил за ней в таблицу рёбер. */
  double nsum[3] = {0.0, 0.0, 0.0};
  if (size == 1) {
    for (int i = 0; i < 12; i++) {
      int a = i / 4, k = i % 4;
      int u = (a + 1) % 3, v = (a + 2) % 3;
      int32_t p[3] = {lo[0], lo[1], lo[2]};
      p[u] += k & 1;
      p[v] += (k >> 1) & 1;
      const hz_hedge *e = hz_htab_find(B->ht, a, p);
      if (e == NULL) continue;
      for (int c = 0; c < 3; c++)
        nsum[c] += e->nrm[c];
    }
  }
  double m = sqrt(nsum[0] * nsum[0] + nsum[1] * nsum[1] + nsum[2] * nsum[2]);
  if (m > 0.0) {
    for (int c = 0; c < 3; c++)
      nsum[c] /= m;
    cl->noct = hz_oct_encode(nsum);
  }
  s->n++;
  return HZ_DC_OK;
}

static void sb_rec(sbctx *B, int32_t ni, const int32_t lo[3], int32_t size, int lvl) {
  if (B->rc != HZ_DC_OK) return;
  hz_dcref r = {ni, {lo[0], lo[1], lo[2]}, size};
  int leaf = (B->t->nd[ni].child0 < 0) || (B->stop != NULL && B->stop(B->sctx, B->t, &r));
  if (leaf) {
    if (hz_dc_hasvert(B->t, ni)) B->rc = sb_push(B, ni, lo, size, lvl);
    return;
  }
  int32_t half = size / 2;
  for (int i = 0; i < 8; i++) {
    int32_t clo[3];
    for (int a = 0; a < 3; a++)
      clo[a] = lo[a] + (((i >> a) & 1) ? half : 0);
    sb_rec(B, B->t->nd[ni].child0 + i, clo, half, lvl + 1);
    if (B->rc != HZ_DC_OK) return;
  }
}

int hz_slice_build(hz_dcslice *s, const hz_dctree *t, const hz_htab *ht, hz_dc_stop stop,
                   void *sctx) {
  s->n = 0;
  sbctx B = {s, t, ht, stop, sctx, HZ_DC_OK};
  int32_t zero[3] = {0, 0, 0};
  sb_rec(&B, 0, zero, (int32_t)1 << t->log2size, 0);
  return B.rc;
}

/* --- чтение ---------------------------------------------------------------- */

/* Позиция ячейки лежит прямо, распаковывать нечего — см. разбор в заголовке:
 * мортонов код стоил 36…41 нс на ячейку против 2.4 нс потока. Порядок ячеек в
 * массиве при этом ОСТАЛСЯ мортоновым: он задаётся порядком обхода 0..7, а не
 * тем, хранится код или нет. */

void hz_slice_vertex(const hz_dcslice *s, int32_t i, double v[3]) {
  double size = (double)((int32_t)1 << (s->lev - s->c[i].lvl));
  double q = (double)((1u << s->vbits) - 1u);
  for (int a = 0; a < 3; a++) {
    double u = (double)(s->c[i].vx[a] >> (8 - s->vbits));
    v[a] = (double)s->c[i].lo[a] + size * (u / q);
  }
}

void hz_slice_normal(const hz_dcslice *s, int32_t i, double n[3]) {
  hz_oct_decode(s->c[i].noct, n);
}
