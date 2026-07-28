/* PLAN_TRANSPORT.md, этап A, шаг 2. Разбор и оговорки — в `raster3.h`. */

#include "transport/raster3.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Значение DG1 в мировой точке по базису ЯЧЕЙКИ. Та же формула, что в сборе:
 * нулевой коэффициент НЕ есть радианс (К39), значение собирается с базисом. */
static double dg1_at(const tr3_mesh *m, int32_t c, const double coef[4], const double x[3]) {
  double s = (double)m->csize[c], v = coef[0];
  for (int a = 0; a < 3; a++) {
    double xu = (x[a] - m->fr.o[a]) / m->fr.u[a];
    v += coef[a + 1] * ((xu - ((double)m->clo[c][a] + 0.5 * s)) / s);
  }
  return v;
}

/* Экранная координата точки. Возврат 0 — точка ЗА камерой либо на её плоскости,
 * и тогда элемент пропускается ЦЕЛИКОМ и считается: обрезать многоугольник по
 * пирамиде — отдельная работа, а молча рисовать половину нельзя. */
static int project(const tr3_camera *c, const double p[3], double *sx, double *sy, double *depth) {
  double d[3];
  for (int a = 0; a < 3; a++)
    d[a] = p[a] - c->eye[a];
  double zf = d[0] * c->fwd[0] + d[1] * c->fwd[1] + d[2] * c->fwd[2];
  if (!(zf > 0.0)) return 0;
  double xr = d[0] * c->right[0] + d[1] * c->right[1] + d[2] * c->right[2];
  double yu = d[0] * c->up[0] + d[1] * c->up[1] + d[2] * c->up[2];
  *sx = 0.5 * ((double)c->w) * (1.0 + (xr / zf) / c->tanx);
  *sy = 0.5 * ((double)c->h) * (1.0 - (yu / zf) / c->tany);
  *depth = zf;
  return 1;
}

/* Один выпуклый многоугольник в буфер. `coef` — [nch][4] коэффициенты DG1
 * этого элемента в базисе ячейки `cell`. */
static void draw_poly(const tr3_raster *r, const tr3_camera *cam, const double (*v)[3], int nv,
                      int32_t cell, const double *coef, int coefstride, double *buf, double *zbuf,
                      tr3_rstats *st) {
  const int nch = r->nch > 0 ? r->nch : 1;
  const int W = cam->w, H = cam->h;
  double sx[TR3_SE_MAXV], sy[TR3_SE_MAXV], dep[TR3_SE_MAXV];
  if (nv < 3 || nv > TR3_SE_MAXV) {
    st->nskip++;
    return;
  }
  double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
  for (int i = 0; i < nv; i++) {
    if (!project(cam, v[i], &sx[i], &sy[i], &dep[i])) {
      st->nskip++;
      return;
    }
    if (sx[i] < xmin) xmin = sx[i];
    if (sx[i] > xmax) xmax = sx[i];
    if (sy[i] < ymin) ymin = sy[i];
    if (sy[i] > ymax) ymax = sy[i];
  }
  /* ЗАЖИМ ДЕЛАЕТСЯ В `double` ДО ПРИВЕДЕНИЯ К `int`, а не после. Приведение
   * величины вне диапазона `int` есть НЕОПРЕДЕЛЁННОЕ ПОВЕДЕНИЕ, и хотя сюда
   * нельзя попасть с непроецированным многоугольником (мы вышли бы раньше),
   * полагаться на это значит держать UB за недостижимой ветвью. */
  if (xmin < 0.0) xmin = 0.0;
  if (ymin < 0.0) ymin = 0.0;
  if (xmax > (double)(W - 1)) xmax = (double)(W - 1);
  if (ymax > (double)(H - 1)) ymax = (double)(H - 1);
  if (!(xmax >= xmin) || !(ymax >= ymin)) return;
  int x0 = (int)floor(xmin), x1 = (int)ceil(xmax);
  int y0 = (int)floor(ymin), y1 = (int)ceil(ymax);
  if (x0 < 0) x0 = 0;
  if (y0 < 0) y0 = 0;
  if (x1 > W - 1) x1 = W - 1;
  if (y1 > H - 1) y1 = H - 1;
  if (x1 < x0 || y1 < y0) return;
  st->nelem++;

  /* Плоскость элемента: нормаль из первых трёх вершин, смещение по первой. */
  double e1[3], e2[3], nrm[3];
  for (int a = 0; a < 3; a++) {
    e1[a] = v[1][a] - v[0][a];
    e2[a] = v[2][a] - v[0][a];
  }
  nrm[0] = e1[1] * e2[2] - e1[2] * e2[1];
  nrm[1] = e1[2] * e2[0] - e1[0] * e2[2];
  nrm[2] = e1[0] * e2[1] - e1[1] * e2[0];
  double nn = sqrt(nrm[0] * nrm[0] + nrm[1] * nrm[1] + nrm[2] * nrm[2]);
  if (!(nn > 0.0)) return;
  for (int a = 0; a < 3; a++)
    nrm[a] /= nn;
  double off = nrm[0] * v[0][0] + nrm[1] * v[0][1] + nrm[2] * v[0][2];

  for (int py = y0; py <= y1; py++)
    for (int px = x0; px <= x1; px++) {
      double o[3], d[3];
      tr3_camera_ray(cam, px, py, o, d);
      double dn = d[0] * nrm[0] + d[1] * nrm[1] + d[2] * nrm[2];
      if (!(fabs(dn) > 0.0)) continue;
      double t = (off - (o[0] * nrm[0] + o[1] * nrm[1] + o[2] * nrm[2])) / dn;
      if (!(t > 0.0)) continue;
      double p[3];
      for (int a = 0; a < 3; a++)
        p[a] = o[a] + t * d[a];
      /* ПРИНАДЛЕЖНОСТЬ — ПО ЗНАКУ, БЕЗ ПОРОГА. Точка лежит в выпуклом
       * многоугольнике, если для всех рёбер знак `(B−A)×(P−A)·n` один и тот же.
       * Нулевой знак (ровно на ребре) принимается, и это делает границу двух
       * соседних элементов детерминированной, а не дырявой. */
      int inside = 1;
      for (int i = 0; i < nv && inside; i++) {
        const double *A = v[i], *B = v[(i + 1) % nv];
        double ab[3], ap[3], cr[3];
        for (int a = 0; a < 3; a++) {
          ab[a] = B[a] - A[a];
          ap[a] = p[a] - A[a];
        }
        cr[0] = ab[1] * ap[2] - ab[2] * ap[1];
        cr[1] = ab[2] * ap[0] - ab[0] * ap[2];
        cr[2] = ab[0] * ap[1] - ab[1] * ap[0];
        if (cr[0] * nrm[0] + cr[1] * nrm[1] + cr[2] * nrm[2] < 0.0) inside = 0;
      }
      if (!inside) continue;
      double zf = (p[0] - cam->eye[0]) * cam->fwd[0] + (p[1] - cam->eye[1]) * cam->fwd[1] +
                  (p[2] - cam->eye[2]) * cam->fwd[2];
      size_t idx = (size_t)py * (size_t)W + (size_t)px;
      if (zbuf[idx] >= 0.0 && zf >= zbuf[idx]) continue;
      zbuf[idx] = zf;
      for (int ch = 0; ch < nch; ch++) {
        double val = dg1_at(r->m, cell, coef + (size_t)ch * (size_t)coefstride, p);
        buf[(size_t)ch * (size_t)W * (size_t)H + idx] = val > 0.0 ? val : 0.0;
      }
      st->npix++;
    }
}

int tr3_raster_render(const tr3_raster *r, const tr3_camera *cam, double *buf, tr3_rstats *st) {
  memset(st, 0, sizeof *st);
  const int nch = r->nch > 0 ? r->nch : 1;
  const int W = cam->w, H = cam->h;
  double *zbuf = malloc((size_t)W * (size_t)H * sizeof(double));
  if (zbuf == NULL) return 1;
  for (size_t i = 0; i < (size_t)W * (size_t)H; i++)
    zbuf[i] = -1.0;
  memset(buf, 0, (size_t)nch * (size_t)W * (size_t)H * sizeof(double));

  /* ГРАНИЧНЫЕ ГРАНИ (стенки). Рисуются ТОЛЬКО те, чья флюидная площадь не ноль:
   * грань полностью твёрдой ячейки света не отдаёт, и рисовать её значило бы
   * закрасить стену там, где перед ней стоит материал. */
  for (int32_t f = 0; f < r->m->nf; f++) {
    if (r->m->f[f].cb >= 0) continue;
    if (r->cut != NULL && !(r->cut->farea[f] > 0.0)) continue;
    double v[4][3];
    tr3_face_corners(r->m, f, v);
    draw_poly(r, cam, v, 4, r->m->f[f].ca, r->bout + (size_t)f * 4, (int)r->m->nf * 4, buf, zbuf,
              st);
  }
  /* ПОВЕРХНОСТНЫЕ ЭЛЕМЕНТЫ (тела). */
  if (r->cut != NULL)
    for (int32_t e = 0; e < r->cut->nse; e++) {
      const tr3_selem *se = &r->cut->se[e];
      if (se->nv < 3) {
        st->nskip++;
        continue;
      }
      draw_poly(r, cam, se->v, (int)se->nv, se->cell, r->sout + (size_t)e * 4,
                (int)(r->cut->nse > 0 ? r->cut->nse : 1) * 4, buf, zbuf, st);
    }

  for (size_t i = 0; i < (size_t)W * (size_t)H; i++)
    if (zbuf[i] < 0.0) st->nempty++;
  free(zbuf);
  return 0;
}
