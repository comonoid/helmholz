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
                        const double v1[3], const double v2[3], int32_t pid) {
  double P[4][3];
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
  double sx[4], sy[4], sd[4];
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
      hc_tri_face(h, &fc[k], base[k], q[0], q[1], q[2], pid);
  }
}
