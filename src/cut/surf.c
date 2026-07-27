/* PLAN_CUT.md Р-3. Таблицы поверхностей, фасетов и боковая таблица. */

#include "cut/surf.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* --- 1. поверхности -------------------------------------------------------- */

int hz_surftab_init(hz_surftab *t) {
  t->cap = 16;
  t->n = 0;
  t->s = calloc((size_t)t->cap, sizeof(hz_surf));
  return t->s == NULL ? 1 : 0;
}

void hz_surftab_free(hz_surftab *t) {
  free(t->s);
  t->s = NULL;
  t->n = t->cap = 0;
}

int32_t hz_surftab_add(hz_surftab *t, const hz_surf *s) {
  if (t->n >= t->cap) {
    int32_t nc = t->cap * 2;
    hz_surf *ns = realloc(t->s, (size_t)nc * sizeof(hz_surf));
    if (ns == NULL) return -1;
    t->s = ns;
    t->cap = nc;
  }
  t->s[t->n] = *s;
  return t->n++;
}

/* --- 2. фасеты ------------------------------------------------------------- */

int hz_facettab_init(hz_facettab *t) {
  t->cap = 64;
  t->n = 0;
  t->f = calloc((size_t)t->cap, sizeof(hz_facet));
  return t->f == NULL ? 1 : 0;
}

void hz_facettab_free(hz_facettab *t) {
  free(t->f);
  t->f = NULL;
  t->n = t->cap = 0;
}

int32_t hz_facettab_add_plane(hz_facettab *t, const hz_frame *fr, const double n[3], double off,
                              int32_t surf, double dmax) {
  if (t->n >= t->cap) {
    int32_t nc = t->cap * 2;
    hz_facet *nf = realloc(t->f, (size_t)nc * sizeof(hz_facet));
    if (nf == NULL) return -1;
    t->f = nf;
    t->cap = nc;
  }
  hz_hspace h;
  hz_frame_plane(fr, n, off, &h); /* перевод в единицы — ЗДЕСЬ и ОДИН РАЗ */
  hz_facet *f = &t->f[t->n];
  memcpy(f->n, h.n, sizeof f->n);
  f->off = h.off;
  f->surf = surf;
  f->dmax = dmax;
  return t->n++;
}

int hz_facet_error(const hz_facet *f, double k, double W, double *err) {
  if (f->dmax < 0.0) return 1; /* не вычислено — потребитель ОБЯЗАН отказаться */
  if (!(W > 0.0)) return 2;
  double kd = k * f->dmax;
  *err = 0.125 * kd * kd * sqrt(f->dmax / W);
  return 0;
}

/* --- фасетизация примитива (Р-5а) ------------------------------------------ */

/* Икосаэдр: 12 вершин, 20 граней. Ориентация граней в таблице НЕ ПРОВЕРЯЕТСЯ на
 * веру — нормаль разворачивается наружу по знаку n·a во время построения, так
 * что опечатка в списке не превращается во внутрь смотрящую полуплоскость. */
static const double ICO_V[12][3] = {
    {-1, 1.6180339887498949, 0},  {1, 1.6180339887498949, 0},   {-1, -1.6180339887498949, 0},
    {1, -1.6180339887498949, 0},  {0, -1, 1.6180339887498949},  {0, 1, 1.6180339887498949},
    {0, -1, -1.6180339887498949}, {0, 1, -1.6180339887498949},  {1.6180339887498949, 0, -1},
    {1.6180339887498949, 0, 1},   {-1.6180339887498949, 0, -1}, {-1.6180339887498949, 0, 1}};
static const int8_t ICO_F[20][3] = {{0, 11, 5}, {0, 5, 1},  {0, 1, 7},   {0, 7, 10}, {0, 10, 11},
                                    {1, 5, 9},  {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                                    {3, 9, 4},  {3, 4, 2},  {3, 2, 6},   {3, 6, 8},  {3, 8, 9},
                                    {4, 9, 5},  {2, 4, 11}, {6, 2, 10},  {8, 6, 7},  {9, 8, 1}};

int32_t hz_surf_facet_sphere(hz_facettab *ft, const hz_frame *fr, const double c[3], double r,
                             int k, hz_fit fit, int32_t surf, int32_t *f0) {
  if (k < 0 || k > HZ_ICO_MAXK || !(r > 0.0)) return -1;
  /* новых вершин за подразделение = число РЁБЕР = 3F/2, а не 4F/2 */
  int32_t nf = 20, nv = 12;
  for (int i = 0; i < k; i++) {
    nv += 3 * nf / 2;
    nf *= 4;
  }
  /* аллокация по месту, без обёрток (CLAUDE.md: gcc-analyzer и обёртки) */
  double *vx = calloc((size_t)nv * 3, sizeof(double));
  int32_t *fa = calloc((size_t)nf * 3, sizeof(int32_t));
  int32_t *fb = calloc((size_t)nf * 3, sizeof(int32_t));
  int32_t *ea = calloc((size_t)nf * 3, sizeof(int32_t)); /* кэш серединных вершин */
  int32_t *eb = calloc((size_t)nf * 3, sizeof(int32_t));
  int32_t *em = calloc((size_t)nf * 3, sizeof(int32_t));
  if (vx == NULL || fa == NULL || fb == NULL || ea == NULL || eb == NULL || em == NULL) {
    free(vx);
    free(fa);
    free(fb);
    free(ea);
    free(eb);
    free(em);
    return -1;
  }

  int32_t cv = 12, cf = 20;
  for (int32_t i = 0; i < 12; i++) {
    double m =
        sqrt(ICO_V[i][0] * ICO_V[i][0] + ICO_V[i][1] * ICO_V[i][1] + ICO_V[i][2] * ICO_V[i][2]);
    for (int a = 0; a < 3; a++)
      vx[i * 3 + a] = ICO_V[i][a] / m;
  }
  for (int32_t i = 0; i < 20; i++)
    for (int a = 0; a < 3; a++)
      fa[i * 3 + a] = ICO_F[i][a];

  for (int lev = 0; lev < k; lev++) {
    int32_t ne = 0, nnf = 0;
    for (int32_t t = 0; t < cf; t++) {
      int32_t v[3] = {fa[t * 3], fa[t * 3 + 1], fa[t * 3 + 2]}, mid[3];
      for (int e = 0; e < 3; e++) {
        int32_t p = v[e], q = v[(e + 1) % 3];
        int32_t lo2 = p < q ? p : q, hi2 = p < q ? q : p; /* КАНОНИЧЕСКОЕ ребро */
        int32_t found = -1;
        for (int32_t s2 = 0; s2 < ne; s2++)
          if (ea[s2] == lo2 && eb[s2] == hi2) {
            found = em[s2];
            break;
          }
        if (found < 0) {
          double m0[3];
          double mm = 0.0;
          for (int a = 0; a < 3; a++) {
            m0[a] = 0.5 * (vx[lo2 * 3 + a] + vx[hi2 * 3 + a]);
            mm += m0[a] * m0[a];
          }
          mm = sqrt(mm);
          for (int a = 0; a < 3; a++)
            vx[cv * 3 + a] = m0[a] / mm; /* вынос на сферу */
          found = cv++;
          ea[ne] = lo2;
          eb[ne] = hi2;
          em[ne] = found;
          ne++;
        }
        mid[e] = found;
      }
      const int32_t tri[4][3] = {{v[0], mid[0], mid[2]},
                                 {v[1], mid[1], mid[0]},
                                 {v[2], mid[2], mid[1]},
                                 {mid[0], mid[1], mid[2]}};
      for (int q = 0; q < 4; q++)
        for (int a = 0; a < 3; a++)
          fb[(nnf + q) * 3 + a] = tri[q][a];
      nnf += 4;
    }
    cf = nnf;
    int32_t *sw = fa;
    fa = fb;
    fb = sw;
  }

  *f0 = ft->n;
  int32_t made = 0;
  for (int32_t t = 0; t < cf; t++) {
    const double *A = &vx[fa[t * 3] * 3], *B = &vx[fa[t * 3 + 1] * 3], *C = &vx[fa[t * 3 + 2] * 3];
    double u1[3], u2[3], n[3];
    for (int a = 0; a < 3; a++) {
      u1[a] = B[a] - A[a];
      u2[a] = C[a] - A[a];
    }
    n[0] = u1[1] * u2[2] - u1[2] * u2[1];
    n[1] = u1[2] * u2[0] - u1[0] * u2[2];
    n[2] = u1[0] * u2[1] - u1[1] * u2[0];
    double nm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
    if (!(nm > 0.0)) continue;
    double sgn = (n[0] * A[0] + n[1] * A[1] + n[2] * A[2]) < 0.0 ? -1.0 : 1.0;
    for (int a = 0; a < 3; a++)
      n[a] = sgn * n[a] / nm; /* наружу — по построению, а не по вере в таблицу */

    /* h — расстояние от центра до плоскости вписанного фасета (единичная сфера) */
    double h = n[0] * A[0] + n[1] * A[1] + n[2] * A[2];
    double hr = h * r;             /* в мире */
    double rho2 = r * r - hr * hr; /* квадрат расстояния вершины от оси фасета */
    if (rho2 < 0.0) rho2 = 0.0;
    double w = 0.0, dmax;
    if (fit == HZ_FIT_MEAN_SAGITTA) {
      /* среднее ρ² по площади правильного треугольника = ρ_v²/4, отсюда вынос
       * (3/8)ρ_v²/r обнуляет среднюю сагитту (Р-5а) */
      w = 0.375 * rho2 / r;
      double s_center = r - hr - w; /* сагитта в центре после выноса */
      double s_vertex = w;          /* и в вершинах, со знаком минус */
      dmax = fabs(s_center) > s_vertex ? fabs(s_center) : s_vertex;
    } else {
      dmax = r - hr; /* однозначное смещение вписанного */
    }
    double off = hr + w;
    off += n[0] * c[0] + n[1] * c[1] + n[2] * c[2];
    if (hz_facettab_add_plane(ft, fr, n, off, surf, dmax) < 0) {
      made = -1;
      break;
    }
    made++;
  }
  free(vx);
  free(fa);
  free(fb);
  free(ea);
  free(eb);
  free(em);
  return made;
}

int hz_facets_for_box(const hz_facettab *ft, int32_t f0, int32_t nf, const int32_t lo[3],
                      const int32_t hi[3], int32_t *out, int max) {
  int nsel = 0;
  for (int32_t j = 0; j < nf; j++) {
    const hz_facet *f = &ft->f[f0 + j];
    double vmax = 0.0, vmin = 0.0;
    for (int a = 0; a < 3; a++) {
      double l = (double)lo[a], hgh = (double)hi[a];
      vmax += f->n[a] > 0.0 ? f->n[a] * hgh : f->n[a] * l;
      vmin += f->n[a] > 0.0 ? f->n[a] * l : f->n[a] * hgh;
    }
    if (vmin > f->off) return HZ_BOX_OUTSIDE; /* коробка целиком снаружи */
    if (vmax <= f->off) continue;             /* коробка целиком внутри: лишняя */
    if (nsel >= max) return HZ_BOX_TOOMANY;
    out[nsel++] = f0 + j;
  }
  return nsel; /* 0 == коробка целиком ВНУТРИ тела (материала полна) */
}

/* --- 3. боковая таблица ---------------------------------------------------- */

int hz_cutmap_init(hz_cutmap *m) {
  m->rcap = 64;
  m->nr = 0;
  m->fcap = 256;
  m->nf = 0;
  m->r = calloc((size_t)m->rcap, sizeof(hz_cutrec));
  if (m->r == NULL) return 1;
  m->fref = calloc((size_t)m->fcap, sizeof(int32_t));
  if (m->fref == NULL) {
    free(m->r);
    m->r = NULL;
    return 1;
  }
  return 0;
}

void hz_cutmap_free(hz_cutmap *m) {
  free(m->r);
  free(m->fref);
  m->r = NULL;
  m->fref = NULL;
  m->nr = m->rcap = m->nf = m->fcap = 0;
}

int hz_cutmap_add(hz_cutmap *m, int32_t cell, const int32_t *fref, int32_t nf) {
  if (nf <= 0) return 1;
  if (m->nr > 0 && cell <= m->r[m->nr - 1].cell) return 2; /* порядок — инвариант */
  if (m->nr >= m->rcap) {
    int32_t nc = m->rcap * 2;
    hz_cutrec *nr = realloc(m->r, (size_t)nc * sizeof(hz_cutrec));
    if (nr == NULL) return 3;
    m->r = nr;
    m->rcap = nc;
  }
  if (m->nf + nf > m->fcap) {
    int32_t nc = m->fcap * 2;
    while (nc < m->nf + nf)
      nc *= 2;
    int32_t *nfr = realloc(m->fref, (size_t)nc * sizeof(int32_t));
    if (nfr == NULL) return 3;
    m->fref = nfr;
    m->fcap = nc;
  }
  for (int32_t j = 0; j < nf; j++)
    m->fref[m->nf + j] = fref[j];
  m->r[m->nr].cell = cell;
  m->r[m->nr].f0 = m->nf;
  m->r[m->nr].nf = nf;
  m->nr++;
  m->nf += nf;
  return 0;
}

const hz_cutrec *hz_cutmap_find(const hz_cutmap *m, int32_t cell) {
  int32_t lo = 0, hi = m->nr - 1;
  while (lo <= hi) {
    int32_t mid = lo + (hi - lo) / 2; /* не (lo+hi)/2: переполнение при 2e8 записях */
    if (m->r[mid].cell == cell) return &m->r[mid];
    if (m->r[mid].cell < cell)
      lo = mid + 1;
    else
      hi = mid - 1;
  }
  return NULL;
}

void hz_cutcur_begin(hz_cutcur *c, const hz_cutmap *m) {
  c->m = m;
  c->i = 0;
}

const hz_cutrec *hz_cutcur_next(hz_cutcur *c, int32_t cell) {
  const hz_cutmap *m = c->m;
  while (c->i < m->nr && m->r[c->i].cell < cell)
    c->i++;
  if (c->i < m->nr && m->r[c->i].cell == cell) return &m->r[c->i];
  return NULL;
}

/* --- мост к ядру ----------------------------------------------------------- */

int hz_cutmap_hspaces(const hz_facettab *ft, const hz_cutmap *m, const hz_cutrec *r, hz_hspace *h,
                      int32_t *hid, int32_t *hid_flipped, int max) {
  if (r->nf > max) return -1;
  for (int32_t j = 0; j < r->nf; j++) {
    int32_t ref = m->fref[r->f0 + j];
    int32_t fi = ref >= 0 ? ref : ~ref;
    if (fi < 0 || fi >= ft->n) return -1;
    const hz_facet *f = &ft->f[fi];
    if (h != NULL) {
      /* Отрицание точно в IEEE: s меняет ровно знак, а параметр пересечения
       * t = s_a/(s_a - s_b) не меняется ВОВСЕ — сосед, держащий другую сторону
       * той же плоскости, строит те же точки побитово (Р-4). */
      for (int a = 0; a < 3; a++)
        h[j].n[a] = ref >= 0 ? f->n[a] : -f->n[a];
      h[j].off = ref >= 0 ? f->off : -f->off;
    }
    if (hid != NULL) hid[j] = ref;
    if (hid_flipped != NULL) hid_flipped[j] = ~ref;
  }
  return (int)r->nf;
}
