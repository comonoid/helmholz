/* phcube.c — полукуб. Разбор — в `phcube.h`. */
#include "phcube.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define HC_PI 3.14159265358979323846

/* ГРАНИ. Верхняя смотрит вдоль нормали; четыре боковых — вдоль касательных осей.
 * Каждая грань задана тройкой (вправо, вверх, вперёд) в раме полукуба, где
 * `z` — нормаль. Боковые грани занимают половину высоты: от касательной
 * плоскости до верха куба. */
typedef struct {
  double r[3], u[3], f[3];
  int w, h; /* пикселей по горизонтали и вертикали */
  double u0, u1, v0, v1;
} hc_face;

static void hc_faces(int R, hc_face fc[5]) {
  const double X[3] = {1.0, 0.0, 0.0}, Y[3] = {0.0, 1.0, 0.0}, Z[3] = {0.0, 0.0, 1.0};
  const double mX[3] = {-1.0, 0.0, 0.0}, mY[3] = {0.0, -1.0, 0.0};
  /* верхняя: вперёд = +z, полный квадрат */
  memcpy(fc[0].r, X, sizeof X);
  memcpy(fc[0].u, Y, sizeof Y);
  memcpy(fc[0].f, Z, sizeof Z);
  fc[0].w = R;
  fc[0].h = R;
  fc[0].u0 = -1.0;
  fc[0].u1 = 1.0;
  fc[0].v0 = -1.0;
  fc[0].v1 = 1.0;
  /* боковые: вперёд = ±x, ±y; «вверх» = +z, и берётся только верхняя половина */
  const double *fw[4] = {X, Y, mX, mY};
  const double *rt[4] = {mY, X, Y, mX};
  for (int k = 0; k < 4; k++) {
    memcpy(fc[k + 1].r, rt[k], 3 * sizeof(double));
    memcpy(fc[k + 1].u, Z, sizeof Z);
    memcpy(fc[k + 1].f, fw[k], 3 * sizeof(double));
    fc[k + 1].w = R;
    fc[k + 1].h = R / 2;
    fc[k + 1].u0 = -1.0;
    fc[k + 1].u1 = 1.0;
    fc[k + 1].v0 = 0.0;
    fc[k + 1].v1 = 1.0;
  }
}

int hz_hcube_init(hz_hcube *h, int R) {
  memset(h, 0, sizeof *h);
  if (R < 2 || (R & 1) != 0) return 1; /* чётное: боковые грани — ровно половина */
  hc_face fc[5];
  hc_faces(R, fc);
  int np = 0;
  for (int k = 0; k < 5; k++)
    np += fc[k].w * fc[k].h;
  h->R = R;
  h->npix = np;
  h->dff = malloc((size_t)np * sizeof *h->dff);
  h->dir = malloc(3 * (size_t)np * sizeof *h->dir);
  h->depth = malloc((size_t)np * sizeof *h->depth);
  h->id = malloc((size_t)np * sizeof *h->id);
  if (h->dff == NULL || h->dir == NULL || h->depth == NULL || h->id == NULL) {
    hz_hcube_free(h);
    return 2;
  }
  /* ДЕЛЬТА-ФОРМ-ФАКТОР ПИКСЕЛЯ — замкнутой формой, без квадратуры.
   * Пиксель грани лежит на единичном расстоянии по «вперёд»; его точка в раме
   * полукуба есть `p = f + a·r + b·u`, площадь пикселя в плоскости грани —
   * `da·db`. Телесный угол с двумя косинусами даёт
   *     ΔF = (1/π) · cosθ_приёмника · cosθ_пикселя · da·db / r²,
   * где `cosθ_пикселя = 1/|p|` (грань перпендикулярна своему «вперёд»),
   * `cosθ_приёмника = p_z/|p|`, `r² = |p|²`. Итого `ΔF = p_z·da·db/(π·|p|⁴)`. */
  int o = 0;
  double s = 0.0;
  for (int k = 0; k < 5; k++) {
    const hc_face *F = &fc[k];
    double da = (F->u1 - F->u0) / F->w, db = (F->v1 - F->v0) / F->h;
    for (int j = 0; j < F->h; j++)
      for (int i = 0; i < F->w; i++) {
        double a = F->u0 + ((double)i + 0.5) * da;
        double b = F->v0 + ((double)j + 0.5) * db;
        double p[3];
        for (int c = 0; c < 3; c++)
          p[c] = F->f[c] + a * F->r[c] + b * F->u[c];
        double q2 = p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
        double q = sqrt(q2);
        double dfv = p[2] * da * db / (HC_PI * q2 * q2);
        h->dff[o] = (float)dfv;
        for (int c = 0; c < 3; c++)
          h->dir[3 * o + c] = p[c] / q;
        s += dfv;
        o++;
      }
  }
  h->sum = s;
  return 0;
}

void hz_hcube_free(hz_hcube *h) {
  free(h->dff);
  free(h->dir);
  free(h->depth);
  free(h->id);
  memset(h, 0, sizeof *h);
}

/* Растеризация треугольника в одну грань. Вершины уже в раме полукуба
 * (приёмник в начале, нормаль по `z`). Отсечение делается ПО ПЛОСКОСТИ ГРАНИ
 * `w = p·f > 0`: без него точки позади камеры грани дают зеркальные призраки. */
static void hc_tri_face(hz_hcube *h, const hc_face *F, int base, const double v0[3],
                        const double v1[3], const double v2[3], int32_t pid, int *tb) {
  /* Инициализация не косметика: цикл отсечения пишет от 0 до 4 вершин, и
   * анализатор обязан видеть, что чтения ниже покрыты записями. */
  double P[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  int np = 0;
  const double *in[3] = {v0, v1, v2};
  double wv[3];
  for (int i = 0; i < 3; i++)
    wv[i] = in[i][0] * F->f[0] + in[i][1] * F->f[1] + in[i][2] * F->f[2];
  const double EPSW = 1e-9;
  for (int i = 0; i < 3; i++) {
    int j = (i + 1) % 3;
    if (wv[i] > EPSW) {
      memcpy(P[np++], in[i], 3 * sizeof(double));
    }
    if ((wv[i] > EPSW) != (wv[j] > EPSW)) {
      double t = (wv[i] - EPSW) / (wv[i] - wv[j]);
      for (int c = 0; c < 3; c++)
        P[np][c] = in[i][c] + t * (in[j][c] - in[i][c]);
      np++;
    }
    if (np >= 4) break;
  }
  if (np < 3) return;
  /* экранные координаты и дальность */
  double sx[4] = {0.0, 0.0, 0.0, 0.0}, sy[4] = {0.0, 0.0, 0.0, 0.0};
  double sd[4] = {0.0, 0.0, 0.0, 0.0};
  double da = (F->u1 - F->u0) / F->w, db = (F->v1 - F->v0) / F->h;
  for (int i = 0; i < np; i++) {
    double w = P[i][0] * F->f[0] + P[i][1] * F->f[1] + P[i][2] * F->f[2];
    double a = (P[i][0] * F->r[0] + P[i][1] * F->r[1] + P[i][2] * F->r[2]) / w;
    double b = (P[i][0] * F->u[0] + P[i][1] * F->u[1] + P[i][2] * F->u[2]) / w;
    sx[i] = (a - F->u0) / da - 0.5;
    sy[i] = (b - F->v0) / db - 0.5;
    sd[i] = sqrt(P[i][0] * P[i][0] + P[i][1] * P[i][1] + P[i][2] * P[i][2]);
  }
  /* веер треугольников из выпуклого многоугольника отсечения */
  for (int t = 1; t + 1 < np; t++) {
    int i0 = 0, i1 = t, i2 = t + 1;
    double xmin = sx[i0], xmax = sx[i0], ymin = sy[i0], ymax = sy[i0];
    const int ii[3] = {i0, i1, i2};
    for (int q = 1; q < 3; q++) {
      if (sx[ii[q]] < xmin) xmin = sx[ii[q]];
      if (sx[ii[q]] > xmax) xmax = sx[ii[q]];
      if (sy[ii[q]] < ymin) ymin = sy[ii[q]];
      if (sy[ii[q]] > ymax) ymax = sy[ii[q]];
    }
    int x0 = (int)floor(xmin), x1 = (int)ceil(xmax);
    int y0 = (int)floor(ymin), y1 = (int)ceil(ymax);
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 > F->w - 1) x1 = F->w - 1;
    if (y1 > F->h - 1) y1 = F->h - 1;
    if (x0 > x1 || y0 > y1) continue;
    double ax = sx[i0], ay = sy[i0], bx = sx[i1], by = sy[i1], cx = sx[i2], cy = sy[i2];
    double area = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax);
    if (!(fabs(area) > 0.0)) continue;
    double inv = 1.0 / area;
    if (tb != NULL) {
      if (x0 < tb[0]) tb[0] = x0;
      if (x1 > tb[1]) tb[1] = x1;
      if (y0 < tb[2]) tb[2] = y0;
      if (y1 > tb[3]) tb[3] = y1;
    }
    for (int y = y0; y <= y1; y++)
      for (int x = x0; x <= x1; x++) {
        double px = (double)x, py = (double)y;
        double w0 = ((bx - px) * (cy - py) - (by - py) * (cx - px)) * inv;
        double w1 = ((cx - px) * (ay - py) - (cy - py) * (ax - px)) * inv;
        double w2 = 1.0 - w0 - w1;
        if (w0 < 0.0 || w1 < 0.0 || w2 < 0.0) continue;
        double d = w0 * sd[i0] + w1 * sd[i1] + w2 * sd[i2];
        int o = base + y * F->w + x;
        if (d < h->depth[o]) {
          h->depth[o] = d;
          h->id[o] = pid;
        }
      }
  }
}

void hz_hcube_draw(hz_hcube *h, const hz_objmesh *m, const int32_t *tri2poly, const double x[3],
                   const double n[3], int32_t skip) {
  for (int i = 0; i < h->npix; i++) {
    h->depth[i] = 1e300;
    h->id[i] = -1;
  }
  /* Рама полукуба: `z` — нормаль, `x`,`y` — любая ортонормированная пара. */
  double ez[3], ex[3], ey[3], nn = 0.0;
  for (int c = 0; c < 3; c++)
    nn += n[c] * n[c];
  nn = sqrt(nn);
  if (!(nn > 0.0)) return;
  for (int c = 0; c < 3; c++)
    ez[c] = n[c] / nn;
  int ax = 0;
  for (int c = 1; c < 3; c++)
    if (fabs(ez[c]) < fabs(ez[ax])) ax = c;
  double t0[3] = {0.0, 0.0, 0.0};
  t0[ax] = 1.0;
  ex[0] = ez[1] * t0[2] - ez[2] * t0[1];
  ex[1] = ez[2] * t0[0] - ez[0] * t0[2];
  ex[2] = ez[0] * t0[1] - ez[1] * t0[0];
  double el = sqrt(ex[0] * ex[0] + ex[1] * ex[1] + ex[2] * ex[2]);
  if (!(el > 0.0)) return;
  for (int c = 0; c < 3; c++)
    ex[c] /= el;
  ey[0] = ez[1] * ex[2] - ez[2] * ex[1];
  ey[1] = ez[2] * ex[0] - ez[0] * ex[2];
  ey[2] = ez[0] * ex[1] - ez[1] * ex[0];

  hc_face fc[5];
  hc_faces(h->R, fc);
  int base[5];
  int o = 0;
  for (int k = 0; k < 5; k++) {
    base[k] = o;
    o += fc[k].w * fc[k].h;
  }
  for (int32_t t = 0; t < m->nt; t++) {
    int32_t pid = (tri2poly != NULL) ? tri2poly[t] : t;
    if (pid == skip) continue;
    double q[3][3];
    int above = 0;
    for (int i = 0; i < 3; i++) {
      const double *p = m->v + 3 * (size_t)m->f[(size_t)t * 3 + (size_t)i];
      double d[3];
      for (int c = 0; c < 3; c++)
        d[c] = p[c] - x[c];
      q[i][0] = d[0] * ex[0] + d[1] * ex[1] + d[2] * ex[2];
      q[i][1] = d[0] * ey[0] + d[1] * ey[1] + d[2] * ey[2];
      q[i][2] = d[0] * ez[0] + d[1] * ez[1] + d[2] * ez[2];
      if (q[i][2] > 0.0) above = 1;
    }
    if (!above) continue; /* целиком под касательной плоскостью — не виден */
    for (int k = 0; k < 5; k++)
      hc_tri_face(h, &fc[k], base[k], q[0], q[1], q[2], pid, NULL);
  }
}

/* --- ОБХОД ПО ПРОСТРАНСТВЕННОМУ ДЕРЕВУ, СПЕРЕДИ НАЗАД (§95) ------------------
 *
 * Три отсева, и каждый со своим счётчиком — иначе нельзя сказать, какое правило
 * работает, а какое относится к пустому множеству:
 *   1. коробка целиком ПОД касательной плоскостью приёмника;
 *   2. коробка не попала ни на одну грань (все восемь углов позади всех граней);
 *   3. ИЕРАРХИЧЕСКИЙ БУФЕР ГЛУБИНЫ: на каждой грани держится сетка плиток с
 *      МАКСИМАЛЬНОЙ глубиной; если ближняя точка коробки дальше максимума во
 *      всех накрытых плитках, поддерево заслонено ЦЕЛИКОМ.
 * Отсев обязан быть КОНСЕРВАТИВНЫМ: `Σf` и число связей не имеют права
 * измениться, и это главная приёмка §95, а не скорость. */
#define HC_TILE 8

typedef struct {
  const hz_objmesh *m;
  const int32_t *t2p;
  const hz_ptree *t;
  hz_hcube *h;
  const hc_face *fc;
  const int *base;
  double x[3], ez[3], ex[3], ey[3];
  int32_t skip;
  int32_t *stamp;
  int32_t mark;
  int nozb;
  hz_hcube_stat *st;
  /* плитки: максимум глубины и признак «пересчитать» */
  double *tmax;
  unsigned char *dirty;
  int toff[5], tw[5], th[5];
} hc_walk;

/* Максимум глубины в плитке; пересчитывается лениво, когда помечена грязной. */
static double hc_tile_max(hc_walk *W, int k, int tx, int ty) {
  int idx = W->toff[k] + ty * W->tw[k] + tx;
  if (W->dirty[idx]) {
    const hc_face *F = &W->fc[k];
    double mx = 0.0;
    int x0 = tx * HC_TILE, y0 = ty * HC_TILE;
    for (int y = y0; y < y0 + HC_TILE && y < F->h; y++)
      for (int x = x0; x < x0 + HC_TILE && x < F->w; x++) {
        double d = W->h->depth[W->base[k] + y * F->w + x];
        if (d > mx) mx = d;
      }
    W->tmax[idx] = mx;
    W->dirty[idx] = 0;
  }
  return W->tmax[idx];
}

static void hc_mark_dirty(hc_walk *W, int k, int x0, int x1, int y0, int y1) {
  for (int ty = y0 / HC_TILE; ty <= y1 / HC_TILE && ty < W->th[k]; ty++)
    for (int tx = x0 / HC_TILE; tx <= x1 / HC_TILE && tx < W->tw[k]; tx++)
      W->dirty[W->toff[k] + ty * W->tw[k] + tx] = 1;
}

/* Заслонена ли коробка узла целиком. Консервативно: «да» только если на КАЖДОЙ
 * грани, куда она попадает, все накрытые плитки ближе её ближней точки. */
static int hc_box_hidden(hc_walk *W, const double blo[3], const double bhi[3], double dnear) {
  int touched = 0;
  for (int k = 0; k < 5; k++) {
    const hc_face *F = &W->fc[k];
    double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
    int behind = 0, ahead = 0;
    double da = (F->u1 - F->u0) / F->w, db = (F->v1 - F->v0) / F->h;
    for (int c = 0; c < 8; c++) {
      double p[3];
      for (int a = 0; a < 3; a++)
        p[a] = (c & (1 << a)) ? bhi[a] : blo[a];
      double w = p[0] * F->f[0] + p[1] * F->f[1] + p[2] * F->f[2];
      if (!(w > 1e-9)) {
        behind = 1;
        continue;
      }
      ahead = 1;
      double a2 = (p[0] * F->r[0] + p[1] * F->r[1] + p[2] * F->r[2]) / w;
      double b2 = (p[0] * F->u[0] + p[1] * F->u[1] + p[2] * F->u[2]) / w;
      double sx = (a2 - F->u0) / da - 0.5, sy = (b2 - F->v0) / db - 0.5;
      if (sx < xmin) xmin = sx;
      if (sx > xmax) xmax = sx;
      if (sy < ymin) ymin = sy;
      if (sy > ymax) ymax = sy;
    }
    if (!ahead) continue; /* грань этой коробки не видит */
    if (behind) {         /* пересекает плоскость грани — считаем всю грань */
      xmin = 0.0;
      ymin = 0.0;
      xmax = F->w - 1;
      ymax = F->h - 1;
    }
    int x0 = (int)floor(xmin), x1 = (int)ceil(xmax);
    int y0 = (int)floor(ymin), y1 = (int)ceil(ymax);
    if (x0 < 0) x0 = 0;
    if (y0 < 0) y0 = 0;
    if (x1 > F->w - 1) x1 = F->w - 1;
    if (y1 > F->h - 1) y1 = F->h - 1;
    if (x0 > x1 || y0 > y1) continue;
    touched = 1;
    for (int ty = y0 / HC_TILE; ty <= y1 / HC_TILE && ty < W->th[k]; ty++)
      for (int tx = x0 / HC_TILE; tx <= x1 / HC_TILE && tx < W->tw[k]; tx++)
        if (hc_tile_max(W, k, tx, ty) > dnear) return 0; /* есть незакрытое */
  }
  return touched; /* попала хоть на одну грань и всюду закрыта */
}

static void hc_draw_tri(hc_walk *W, int32_t tr) {
  int32_t pid = (W->t2p != NULL) ? W->t2p[tr] : tr;
  if (pid == W->skip) return;
  if (W->stamp[tr] == W->mark) {
    W->st->ndup++;
    return;
  }
  W->stamp[tr] = W->mark;
  double q[3][3];
  int above = 0;
  for (int i = 0; i < 3; i++) {
    const double *p = W->m->v + 3 * (size_t)W->m->f[(size_t)tr * 3 + (size_t)i];
    double d[3];
    for (int c = 0; c < 3; c++)
      d[c] = p[c] - W->x[c];
    q[i][0] = d[0] * W->ex[0] + d[1] * W->ex[1] + d[2] * W->ex[2];
    q[i][1] = d[0] * W->ey[0] + d[1] * W->ey[1] + d[2] * W->ey[2];
    q[i][2] = d[0] * W->ez[0] + d[1] * W->ez[1] + d[2] * W->ez[2];
    if (q[i][2] > 0.0) above = 1;
  }
  if (!above) return;
  W->st->nsetup++;
  for (int k = 0; k < 5; k++) {
    const hc_face *F = &W->fc[k];
    /* ГРЯЗНЫМИ ПОМЕЧАЮТСЯ ТОЛЬКО ЗАДЕТЫЕ ПЛИТКИ. Первая редакция метила всю
     * грань на КАЖДЫЙ треугольник, отчего буфер глубины пересчитывался целиком
     * и дерево вышло медленнее прямого перебора (12.5 с против 8.3). */
    int tb[4] = {F->w, -1, F->h, -1};
    hc_tri_face(W->h, F, W->base[k], q[0], q[1], q[2], pid, W->nozb ? NULL : tb);
    if (!W->nozb && tb[1] >= tb[0]) hc_mark_dirty(W, k, tb[0], tb[1], tb[2], tb[3]);
  }
}

static void hc_walk_node(hc_walk *W, int32_t nid) {
  const hz_ptnode *N = &W->t->nd[nid];
  W->st->nnode++;
  /* коробка в раме полукуба: габарит восьми углов */
  double blo[3] = {1e300, 1e300, 1e300}, bhi[3] = {-1e300, -1e300, -1e300};
  double dnear = 1e300;
  for (int c = 0; c < 8; c++) {
    double p[3], d[3];
    for (int a = 0; a < 3; a++)
      p[a] = ((c & (1 << a)) ? N->hi[a] : N->lo[a]) - W->x[a];
    d[0] = p[0] * W->ex[0] + p[1] * W->ex[1] + p[2] * W->ex[2];
    d[1] = p[0] * W->ey[0] + p[1] * W->ey[1] + p[2] * W->ey[2];
    d[2] = p[0] * W->ez[0] + p[1] * W->ez[1] + p[2] * W->ez[2];
    for (int a = 0; a < 3; a++) {
      if (d[a] < blo[a]) blo[a] = d[a];
      if (d[a] > bhi[a]) bhi[a] = d[a];
    }
  }
  if (!(bhi[2] > 0.0)) { /* отсев 1: вся коробка под касательной плоскостью */
    W->st->ncull_h++;
    return;
  }
  /* ближняя точка коробки к началу */
  double dn2 = 0.0;
  for (int a = 0; a < 3; a++) {
    double v = (blo[a] > 0.0) ? blo[a] : ((bhi[a] < 0.0) ? -bhi[a] : 0.0);
    dn2 += v * v;
  }
  dnear = sqrt(dn2);
  if (!W->nozb && hc_box_hidden(W, blo, bhi, dnear)) { /* отсев 3 */
    W->st->ncull_z++;
    return;
  }
  /* Треугольники есть у ЛЮБОГО узла, а не только у листа: те, что не влезли
   * целиком ни в одного ребёнка, остаются здесь (см. `ptree.c`). */
  for (int32_t i = 0; i < N->ntri; i++)
    hc_draw_tri(W, W->t->ref[N->t0 + i]);
  if (N->child < 0) return;
  /* СПЕРЕДИ НАЗАД: дети по возрастанию расстояния от приёмной точки. */
  int ord[8];
  double dist[8];
  for (int k = 0; k < 8; k++) {
    const hz_ptnode *C = &W->t->nd[N->child + k];
    double s = 0.0;
    for (int a = 0; a < 3; a++) {
      double v = 0.0;
      if (W->x[a] < C->lo[a])
        v = C->lo[a] - W->x[a];
      else if (W->x[a] > C->hi[a])
        v = W->x[a] - C->hi[a];
      s += v * v;
    }
    dist[k] = s;
    ord[k] = k;
  }
  for (int a = 1; a < 8; a++) {
    int v = ord[a];
    double dv = dist[v];
    int b = a - 1;
    while (b >= 0 && dist[ord[b]] > dv) {
      ord[b + 1] = ord[b];
      b--;
    }
    ord[b + 1] = v;
    dist[v] = dv;
  }
  for (int k = 0; k < 8; k++)
    hc_walk_node(W, N->child + ord[k]);
}

void hz_hcube_draw_tree(hz_hcube *h, const hz_objmesh *m, const int32_t *tri2poly,
                        const hz_ptree *t, const double x[3], const double n[3], int32_t skip,
                        int32_t *stamp, int32_t mark, int nozb, hz_hcube_stat *st) {
  for (int i = 0; i < h->npix; i++) {
    h->depth[i] = 1e300;
    h->id[i] = -1;
  }
  hc_face fc[5];
  hc_faces(h->R, fc);
  int base[5], o = 0;
  for (int k = 0; k < 5; k++) {
    base[k] = o;
    o += fc[k].w * fc[k].h;
  }
  hc_walk W;
  memset(&W, 0, sizeof W);
  W.m = m;
  W.t2p = tri2poly;
  W.t = t;
  W.h = h;
  W.fc = fc;
  W.base = base;
  W.skip = skip;
  W.stamp = stamp;
  W.mark = mark;
  W.nozb = nozb;
  W.st = st;
  for (int c = 0; c < 3; c++)
    W.x[c] = x[c];
  double nn = 0.0;
  for (int c = 0; c < 3; c++)
    nn += n[c] * n[c];
  nn = sqrt(nn);
  if (!(nn > 0.0)) return;
  for (int c = 0; c < 3; c++)
    W.ez[c] = n[c] / nn;
  int ax = 0;
  for (int c = 1; c < 3; c++)
    if (fabs(W.ez[c]) < fabs(W.ez[ax])) ax = c;
  double t0v[3] = {0.0, 0.0, 0.0};
  t0v[ax] = 1.0;
  W.ex[0] = W.ez[1] * t0v[2] - W.ez[2] * t0v[1];
  W.ex[1] = W.ez[2] * t0v[0] - W.ez[0] * t0v[2];
  W.ex[2] = W.ez[0] * t0v[1] - W.ez[1] * t0v[0];
  double el = sqrt(W.ex[0] * W.ex[0] + W.ex[1] * W.ex[1] + W.ex[2] * W.ex[2]);
  if (!(el > 0.0)) return;
  for (int c = 0; c < 3; c++)
    W.ex[c] /= el;
  W.ey[0] = W.ez[1] * W.ex[2] - W.ez[2] * W.ex[1];
  W.ey[1] = W.ez[2] * W.ex[0] - W.ez[0] * W.ex[2];
  W.ey[2] = W.ez[0] * W.ex[1] - W.ez[1] * W.ex[0];

  int nt = 0;
  for (int k = 0; k < 5; k++) {
    W.tw[k] = (fc[k].w + HC_TILE - 1) / HC_TILE;
    W.th[k] = (fc[k].h + HC_TILE - 1) / HC_TILE;
    W.toff[k] = nt;
    nt += W.tw[k] * W.th[k];
  }
  double tmax[4096];
  unsigned char dirty[4096];
  if (nt > 4096) return;
  for (int i = 0; i < nt; i++) {
    tmax[i] = 1e300;
    dirty[i] = 0;
  }
  W.tmax = tmax;
  W.dirty = dirty;
  hc_walk_node(&W, 0);
}
