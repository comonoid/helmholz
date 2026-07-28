/* Ядро развёртки. Разбор и оговорки — в `psweep.h`. */

#include "psweep.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* Время фаз меряется ВНУТРИ, а не снаружи: Ш3 требует доли растеризации против
 * доли редукции, а снаружи видна только сумма. */
static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

int hz_ptrans_init(hz_ptrans *t, const hz_polyset *ps, const hz_objmesh *m) {
  memset(t, 0, sizeof *t);
  t->ps = ps;
  t->np = ps->np;
  t->rho = malloc((size_t)ps->np * sizeof *t->rho);
  t->Le = calloc((size_t)ps->np, sizeof *t->Le);
  t->E = calloc((size_t)ps->np * 3, sizeof *t->E);
  t->acc = calloc((size_t)ps->np * 3, sizeof *t->acc);
  if (t->rho == NULL || t->Le == NULL || t->E == NULL || t->acc == NULL) {
    hz_ptrans_free(t);
    return 2;
  }
  for (int32_t k = 0; k < ps->np; k++) {
    int32_t mi = ps->p[k].mtl;
    t->rho[k] = (m != NULL && mi >= 0 && mi < m->nmtl) ? m->mtl[mi].kd : 0.5;
  }
  return 0;
}

void hz_ptrans_free(hz_ptrans *t) {
  free(t->rho);
  free(t->Le);
  free(t->E);
  free(t->acc);
  memset(t, 0, sizeof *t);
}

double hz_ptrans_lout(const hz_ptrans *t, int32_t k, double u, double v) {
  const double *c = t->E + (size_t)k * 3;
  double e = c[0] + c[1] * u + c[2] * v;
  if (e < 0.0) e = 0.0; /* ограничитель уже прошёл; это защита чтения на кромке */
  return t->Le[k] + t->rho[k] * e / M_PI;
}

void hz_psweep_zero(hz_ptrans *t) {
  memset(t->acc, 0, (size_t)t->np * 3 * sizeof *t->acc);
}

/* --- 3×3: Гаусс с частичным выбором ---------------------------------------- */

static int solve3s(double A[3][3], double b[3]) {
  for (int c = 0; c < 3; c++) {
    int piv = c;
    for (int r = c + 1; r < 3; r++)
      if (fabs(A[r][c]) > fabs(A[piv][c])) piv = r;
    if (!(fabs(A[piv][c]) > 0.0)) return 1;
    if (piv != c) {
      for (int k = 0; k < 3; k++) {
        double s = A[c][k];
        A[c][k] = A[piv][k];
        A[piv][k] = s;
      }
      double s = b[c];
      b[c] = b[piv];
      b[piv] = s;
    }
    for (int r = 0; r < 3; r++) {
      if (r == c) continue;
      double f = A[r][c] / A[c][c];
      for (int k = 0; k < 3; k++)
        A[r][k] -= f * A[c][k];
      b[r] -= f * b[c];
    }
  }
  for (int c = 0; c < 3; c++)
    b[c] /= A[c][c];
  return 0;
}

/* --- ОДИН ФРАГМЕНТ: физика пары «излучатель -> приёмник» ---------------------
 *
 * Вынесено отдельно ровно потому, что редукция существует в ДВУХ раскладках
 * (Ш3 меряет, какая из них где упирается), и физика в них обязана быть одна и
 * та же. Всё, что различается у раскладок, — способ добраться до следующего
 * фрагмента; всё, что здесь, — общее.
 *
 * `*pk` / `*pu` / `*pv` — предыдущий фрагмент вдоль луча; на входе первого
 * фрагмента `*pk < 0`. */
static inline void frag_step(const hz_ptrans *t, const hz_pview *v, double wq, double h2,
                             double Lsky, double *acc, const double r[3], int32_t k, double depth,
                             int32_t *pk, double *pu, double *pv, int64_t *npair) {
  const hz_poly *P = &t->ps->p[k];
  double nw = P->n[0] * v->w[0] + P->n[1] * v->w[1] + P->n[2] * v->w[2];
  double q[3];
  for (int a = 0; a < 3; a++)
    q[a] = r[a] + depth * v->w[a] - P->org[a];
  double u = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
  double vv = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];

  if (nw < 0.0) { /* лицом навстречу ω ⇒ принимает */
    double L;
    if (*pk < 0) {
      L = Lsky; /* выше по лучу пусто: фон */
    } else {
      const hz_poly *Q = &t->ps->p[*pk];
      double nwp = Q->n[0] * v->w[0] + Q->n[1] * v->w[1] + Q->n[2] * v->w[2];
      /* Излучает вдоль ω только лицевая сторона; изнанка ЗАСЛОНЯЕТ, и это и
       * есть тень — вычисленная, а не взятая из таблицы видимости. */
      L = (nwp > 0.0) ? hz_ptrans_lout(t, *pk, *pu, *pv) : 0.0;
    }
    if (L > 0.0) {
      double c = wq * h2 * L;
      acc[(size_t)k * 3 + 0] += c;
      acc[(size_t)k * 3 + 1] += c * u;
      acc[(size_t)k * 3 + 2] += c * vv;
      (*npair)++;
    }
  }
  *pk = k;
  *pu = u;
  *pv = vv;
}

/* --- редукция, раскладка СПЛОШНЫХ ПРОБЕГОВ ---------------------------------- */

static void reduce_fbuf(const hz_ptrans *t, const hz_pview *v, const hz_fbuf *fb, double wq,
                        double Lsky, double *acc, hz_pstats *st) {
  const double h2 = v->h * v->h;
  const int64_t npix = (int64_t)fb->W * fb->H;
  for (int64_t px = 0; px < npix; px++) {
    int32_t b = fb->start[px], e = fb->start[px + 1];
    if (b == e) continue;
    double r[3];
    hz_pview_origin(v, (int32_t)(px % fb->W), (int32_t)(px / fb->W), r);
    int32_t pk = -1;
    double pu = 0.0, pv = 0.0;
    for (int32_t f = b; f < e; f++)
      frag_step(t, v, wq, h2, Lsky, acc, r, fb->poly[f], fb->depth[f], &pk, &pu, &pv, &st->npair);
  }
}

/* --- редукция, раскладка ОДНОСВЯЗНЫХ СПИСКОВ -------------------------------- */

static void reduce_abuf(const hz_ptrans *t, const hz_pview *v, const hz_abuf *ab, double wq,
                        double Lsky, double *acc, hz_pstats *st) {
  const double h2 = v->h * v->h;
  const int64_t npix = (int64_t)ab->W * ab->H;
  for (int64_t px = 0; px < npix; px++) {
    int32_t f = ab->head[px];
    if (f < 0) continue;
    double r[3];
    hz_pview_origin(v, (int32_t)(px % ab->W), (int32_t)(px / ab->W), r);
    int32_t pk = -1;
    double pu = 0.0, pv = 0.0;
    while (f >= 0) {
      frag_step(t, v, wq, h2, Lsky, acc, r, ab->poly[f], ab->depth[f], &pk, &pu, &pv, &st->npair);
      f = ab->next[f];
    }
  }
}

/* --- поток --------------------------------------------------------------- */

typedef struct {
  hz_span *sp;
  int64_t nsp, cap;
  hz_fbuf fb;
  hz_abuf ab;
  double *acc;
  hz_pstats st;
} pw_thread;

static void scene_box(const hz_polyset *ps, double lo[3], double hi[3]) {
  for (int a = 0; a < 3; a++) {
    lo[a] = 1e300;
    hi[a] = -1e300;
  }
  for (int32_t k = 0; k < ps->np; k++) {
    const hz_poly *P = &ps->p[k];
    for (int32_t l = P->l0; l < P->l0 + P->nloop; l++)
      for (int32_t b = ps->loop[l]; b < ps->loop[l + 1]; b++) {
        double x[3];
        hz_poly_world(P, ps->bv[(size_t)b * 2], ps->bv[(size_t)b * 2 + 1], x);
        for (int a = 0; a < 3; a++) {
          if (x[a] < lo[a]) lo[a] = x[a];
          if (x[a] > hi[a]) hi[a] = x[a];
        }
      }
  }
}

/* Одно направление целиком: полосы -> раскладка -> сортировка -> редукция. */
static int gather_dir(const hz_ptrans *t, const hz_pview *v, pw_thread *w, double wq, double Lsky,
                      int layout) {
  hz_rstats rs;
  double t0 = now_s();
  if (hz_prast_spans(&w->sp, &w->nsp, &w->cap, v, t->ps, &rs) != 0) return 2;
  double t1 = now_s();
  int rc;
  if (layout == HZ_LAYOUT_LIST)
    rc = hz_abuf_build(&w->ab, w->sp, w->nsp, v->W, v->H);
  else
    rc = hz_fbuf_build(&w->fb, w->sp, w->nsp, v->W, v->H);
  if (rc != 0) {
    w->st.nover += rs.nfrag;
    return 0; /* ёмкость: докладывается через nover, а не глушится */
  }
  double t2 = now_s();
  if (layout == HZ_LAYOUT_LIST)
    hz_abuf_sort(&w->ab);
  else
    hz_fbuf_sort(&w->fb);
  double t3 = now_s();
  if (layout == HZ_LAYOUT_LIST)
    reduce_abuf(t, v, &w->ab, wq, Lsky, w->acc, &w->st);
  else
    reduce_fbuf(t, v, &w->fb, wq, Lsky, w->acc, &w->st);
  double t4 = now_s();

  w->st.t_raster += (t1 - t0) + (t2 - t1);
  w->st.t_sort += t3 - t2;
  w->st.t_reduce += t4 - t3;
  w->st.nfrag += rs.nfrag;
  w->st.nspan += rs.nspan;
  /* Трафик: запись фрагмента (poly + depth) и чтение его же редукцией, плюс
   * список у списочной раскладки. Полигонные чтения сюда НЕ входят — они
   * когерентны и считать их байтами значило бы приписать схеме трафик, которого
   * кэш не производит. */
  int64_t per = (int64_t)(sizeof(int32_t) + sizeof(double)) * 2;
  if (layout == HZ_LAYOUT_LIST) per += (int64_t)sizeof(int32_t) * 2;
  w->st.nbytes += rs.nfrag * per;
  return 0;
}

int hz_psweep_gather(hz_ptrans *t, const tr3_dirs *d, double h, int nthr, double Lsky, int layout,
                     hz_pstats *st) {
  const hz_polyset *ps = t->ps;
  double lo[3], hi[3];
  scene_box(ps, lo, hi);
  if (!(lo[0] <= hi[0])) return 1;

  int nt = nthr;
#ifdef _OPENMP
  if (nt <= 0) nt = omp_get_max_threads();
#else
  nt = 1;
#endif
  if (nt < 1) nt = 1;

  /* ЭКРАН ЗАВОДИТСЯ ПО МАКСИМУМУ ПО ВСЕМ НАПРАВЛЕНИЯМ, а не по первому. Разные
   * `ω` дают разные габариты проекции, и буфер «как у направления ноль» тихо
   * выбросил бы часть развёртки — потеря, которую не поймает ни один инвариант,
   * потому что баланс сойдётся и на неполном наборе. */
  int64_t maxpix = 0;
  for (int m = 0; m < d->n; m++) {
    double w[3] = {d->ox[m], d->oy[m], d->oz[m]};
    hz_pview v;
    if (hz_pview_make(&v, w, lo, hi, h) != 0) return 1;
    if ((int64_t)v.W * v.H > maxpix) maxpix = (int64_t)v.W * v.H;
  }
  if (maxpix <= 0) return 1;

  /* Ёмкость: площадь сцены на площадь пикселя, вчетверо с запасом. Переполнение
   * НЕ глушится — `nover` докладывается, и ненулевой означает неполный ответ. */
  double area = 0.0;
  for (int32_t k = 0; k < ps->np; k++)
    area += ps->p[k].area;
  int64_t cap = (int64_t)(4.0 * area / (h * h)) + 1024;
  if (cap > 2000000000LL) cap = 2000000000LL;

  int rc = 0;
  pw_thread *w = calloc((size_t)nt, sizeof *w);
  if (w == NULL) return 2;
  for (int i = 0; i < nt; i++) {
    w[i].acc = calloc((size_t)t->np * 3, sizeof *w[i].acc);
    if (w[i].acc == NULL) rc = 2;
    if (layout == HZ_LAYOUT_LIST) {
      if (hz_abuf_init(&w[i].ab, maxpix, cap) != 0) rc = 2;
    } else {
      if (hz_fbuf_init(&w[i].fb, maxpix, cap) != 0) rc = 2;
    }
  }

  if (rc == 0) {
#pragma omp parallel for schedule(dynamic, 1) num_threads(nt)
    for (int m = 0; m < d->n; m++) {
      int id = 0;
#ifdef _OPENMP
      id = omp_get_thread_num();
#endif
      double wd[3] = {d->ox[m], d->oy[m], d->oz[m]};
      hz_pview v;
      if (hz_pview_make(&v, wd, lo, hi, h) != 0) continue;
      gather_dir(t, &v, &w[id], d->w[m], Lsky, layout);
    }
    for (int i = 0; i < nt; i++) {
      st->nfrag += w[i].st.nfrag;
      st->nspan += w[i].st.nspan;
      st->npair += w[i].st.npair;
      st->nover += w[i].st.nover;
      st->nbytes += w[i].st.nbytes;
      st->t_raster += w[i].st.t_raster;
      st->t_sort += w[i].st.t_sort;
      st->t_reduce += w[i].st.t_reduce;
      for (int64_t k = 0; k < (int64_t)t->np * 3; k++)
        t->acc[k] += w[i].acc[k];
    }
  }
  for (int i = 0; i < nt; i++) {
    free(w[i].acc);
    free(w[i].sp);
    hz_fbuf_free(&w[i].fb);
    hz_abuf_free(&w[i].ab);
  }
  free(w);
  return rc;
}

int hz_psweep_direct(hz_ptrans *t, const double w[3], double Eperp, double h, hz_pstats *st) {
  const hz_polyset *ps = t->ps;
  double lo[3], hi[3];
  scene_box(ps, lo, hi);
  hz_pview v;
  if (hz_pview_make(&v, w, lo, hi, h) != 0) return 1;
  double area = 0.0;
  for (int32_t k = 0; k < ps->np; k++)
    area += ps->p[k].area;
  pw_thread th;
  memset(&th, 0, sizeof th);
  int64_t cap = (int64_t)(4.0 * area / (h * h)) + 1024;
  if (hz_fbuf_init(&th.fb, (int64_t)v.W * v.H, cap) != 0) return 2;
  hz_rstats rs;
  int rc = hz_prast_spans(&th.sp, &th.nsp, &th.cap, &v, ps, &rs);
  if (rc == 0) rc = hz_fbuf_build(&th.fb, th.sp, th.nsp, v.W, v.H);
  if (rc != 0) {
    free(th.sp);
    hz_fbuf_free(&th.fb);
    return rc;
  }
  hz_fbuf_sort(&th.fb);
  st->nfrag += rs.nfrag;
  st->nspan += rs.nspan;

  const double h2 = v.h * v.h;
  const int64_t npix = (int64_t)v.W * v.H;
  for (int64_t px = 0; px < npix; px++) {
    int32_t b = th.fb.start[px], e = th.fb.start[px + 1];
    if (b == e) continue;
    /* ТОЛЬКО ПЕРВЫЙ фрагмент: коллимированный пучок гасится первой же
     * непрозрачной поверхностью, и это и есть тень. */
    int32_t k = th.fb.poly[b];
    const hz_poly *P = &ps->p[k];
    double nw = P->n[0] * v.w[0] + P->n[1] * v.w[1] + P->n[2] * v.w[2];
    if (!(nw < 0.0)) continue;
    double r[3], q[3];
    hz_pview_origin(&v, (int32_t)(px % v.W), (int32_t)(px / v.W), r);
    for (int a = 0; a < 3; a++)
      q[a] = r[a] + th.fb.depth[b] * v.w[a] - P->org[a];
    double u = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
    double vv = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
    double c = Eperp * h2;
    t->acc[(size_t)k * 3 + 0] += c;
    t->acc[(size_t)k * 3 + 1] += c * u;
    t->acc[(size_t)k * 3 + 2] += c * vv;
    st->npair++;
  }
  free(th.sp);
  hz_fbuf_free(&th.fb);
  return 0;
}

void hz_psweep_solve(hz_ptrans *t, hz_pstats *st) {
  const hz_polyset *ps = t->ps;
  double t0 = now_s();
  double sd = 0.0, ss = 0.0, phi_in = 0.0, phi_out = 0.0;
  int64_t nclip = 0;
  for (int32_t k = 0; k < t->np; k++) {
    const hz_poly *P = &ps->p[k];
    const double *M = P->mom;
    double A[3][3] = {{M[0], M[1], M[2]}, {M[1], M[3], M[4]}, {M[2], M[4], M[5]}};
    double b[3] = {t->acc[(size_t)k * 3], t->acc[(size_t)k * 3 + 1], t->acc[(size_t)k * 3 + 2]};
    double c[3];
    if (M[0] > 0.0 && solve3s(A, b) == 0) {
      c[0] = b[0];
      c[1] = b[1];
      c[2] = b[2];
    } else {
      /* Вырожденные моменты (щепка) — только постоянная часть. */
      c[0] = (M[0] > 0.0) ? t->acc[(size_t)k * 3] / M[0] : 0.0;
      c[1] = 0.0;
      c[2] = 0.0;
    }

    /* ОГРАНИЧИТЕЛЬ. Минимум линейного поля по габариту края; если он ниже нуля,
     * градиент сжимается, а c0 правится так, что ∫E dA не меняется.
     * У ПОСТОЯННОГО ПОЛЯ ГРАДИЕНТ НУЛЕВОЙ ⇒ ветка не срабатывает вовсе, и
     * требование К6 «ограничитель обязан сохранять константы» выполнено ПО
     * ПОСТРОЕНИЮ, а не по замеру. */
    double du = (fabs(P->uvlo[0]) > fabs(P->uvhi[0])) ? fabs(P->uvlo[0]) : fabs(P->uvhi[0]);
    double dv = (fabs(P->uvlo[1]) > fabs(P->uvhi[1])) ? fabs(P->uvlo[1]) : fabs(P->uvhi[1]);
    double swing = fabs(c[1]) * du + fabs(c[2]) * dv;
    if (swing > 0.0 && c[0] < swing) {
      double alpha = (c[0] > 0.0) ? c[0] / swing : 0.0;
      double keep = c[1] * M[1] + c[2] * M[2];
      c[1] *= alpha;
      c[2] *= alpha;
      if (M[0] > 0.0) c[0] += (1.0 - alpha) * keep / M[0];
      nclip++;
    }

    double *e = t->E + (size_t)k * 3;
    for (int a = 0; a < 3; a++) {
      sd += fabs(c[a] - e[a]);
      ss += fabs(c[a]);
      e[a] = c[a];
    }
    double ie = e[0] * M[0] + e[1] * M[1] + e[2] * M[2];
    phi_in += ie;
    phi_out += M_PI * t->Le[k] * M[0] + t->rho[k] * ie;
  }
  st->nclip = nclip;
  st->phi_in = phi_in;
  st->phi_out = phi_out;
  st->dE = (ss > 0.0) ? sd / ss : 0.0;
  st->t_solve += now_s() - t0;
}

int hz_psweep_bounce(hz_ptrans *t, const tr3_dirs *d, double h, int nthr, double Lsky, int layout,
                     hz_pstats *st) {
  memset(st, 0, sizeof *st);
  hz_psweep_zero(t);
  int rc = hz_psweep_gather(t, d, h, nthr, Lsky, layout, st);
  if (rc != 0) return rc;
  hz_psweep_solve(t, st);
  return 0;
}
