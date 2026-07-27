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
