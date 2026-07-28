/* Теневые разрезы. Разбор и оговорки — в `pcut.h`. */

#include "pcut.h"
#include "pdirect.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Клиньев на один участок — потолок. Каждый добавляет бит подписи, а число
 * кусков растёт как произведение; без потолка один пол у кресла с десятком
 * ножек развалился бы на тысячи щепок. Взятые — с наибольшим контрастом (Т2). */
#define PCUT_MAXW 8
/* Кусков на треугольник — страховка от вырожденных наборов плоскостей. */
#define PCUT_MAXFRAG 64

void hz_wedges_free(hz_wedges *ws) {
  free(ws->w);
  memset(ws, 0, sizeof *ws);
}

static int cmp_wedge(const void *a, const void *b) {
  const hz_wedge *x = a, *y = b;
  return (x->prio > y->prio) ? -1 : ((x->prio < y->prio) ? 1 : 0);
}

/* --- хеш направленных полурёбер -------------------------------------------- */

typedef struct {
  uint64_t *key;
  int32_t *tri;
  int64_t n;
  uint64_t mask;
} he_tab;

static uint64_t hmix(uint64_t k) {
  k ^= k >> 33;
  k *= 0xFF51AFD7ED558CCDULL;
  k ^= k >> 33;
  k *= 0xC4CEB9FE1A85EC53ULL;
  k ^= k >> 33;
  return k;
}

static int64_t he_slot(he_tab *h, uint64_t k, int insert) {
  uint64_t kk = k + 1, i = hmix(kk) & h->mask;
  for (;;) {
    if (h->key[i] == kk) return (int64_t)i;
    if (h->key[i] == 0) {
      if (!insert) return -1;
      h->key[i] = kk;
      return (int64_t)i;
    }
    i = (i + 1) & h->mask;
  }
}

int hz_wedges_build(hz_wedges *ws, const hz_objmesh *m, const hz_polyset *ps, const int32_t *src,
                    int nsrc, const double *Le, const hz_cutcfg *cfg) {
  memset(ws, 0, sizeof *ws);
  int32_t *wid = malloc((size_t)m->nv * sizeof *wid);
  double *wpos = malloc((size_t)m->nv * 3 * sizeof *wpos);
  if (wid == NULL || wpos == NULL) {
    free(wid);
    free(wpos);
    return 2;
  }
  int32_t nw = 0;
  if (hz_poly_weld(m, wid, &nw, wpos) != 0) {
    free(wid);
    free(wpos);
    return 2;
  }

  he_tab h;
  memset(&h, 0, sizeof h);
  h.n = 4;
  while (h.n < 6 * (int64_t)m->nt)
    h.n *= 2;
  h.mask = (uint64_t)h.n - 1;
  h.key = calloc((size_t)h.n, sizeof *h.key);
  h.tri = malloc((size_t)h.n * sizeof *h.tri);
  if (h.key == NULL || h.tri == NULL) {
    free(h.key);
    free(h.tri);
    free(wid);
    free(wpos);
    return 2;
  }
  for (int32_t t = 0; t < m->nt; t++)
    for (int e = 0; e < 3; e++) {
      int32_t a = wid[m->f[(size_t)t * 3 + (size_t)e]];
      int32_t b = wid[m->f[(size_t)t * 3 + (size_t)((e + 1) % 3)]];
      h.tri[he_slot(&h, (uint64_t)a * (uint64_t)nw + (uint64_t)b, 1)] = t;
    }

  /* Средняя облучённость сцены — знаменатель правила B. Считается ДО
   * построения, из мощности источников и площади: ровно то, что план и
   * обещает («величина известна ДО построения»). */
  double Atot = 0.0, Pow = 0.0;
  for (int32_t k = 0; k < ps->np; k++)
    Atot += ps->p[k].area;
  for (int s = 0; s < nsrc; s++)
    Pow += M_PI * Le[src[s]] * ps->p[src[s]].mom[0];
  double Emean = (Atot > 0.0) ? Pow / Atot : 1.0;
  if (!(Emean > 0.0)) Emean = 1.0;

  /* геометрическая нормаль треугольника */
  double gn[3];
  for (int32_t t = 0; t < m->nt; t++) {
    for (int e = 0; e < 3; e++) {
      int32_t a = wid[m->f[(size_t)t * 3 + (size_t)e]];
      int32_t b = wid[m->f[(size_t)t * 3 + (size_t)((e + 1) % 3)]];
      if (a > b) continue; /* каждое НЕнаправленное ребро разбираем один раз */
      int64_t rv = he_slot(&h, (uint64_t)b * (uint64_t)nw + (uint64_t)a, 0);
      int32_t t2 = (rv >= 0) ? h.tri[rv] : -1;

      double p[3][3];
      hz_obj_tri(m, t, p);
      double e1[3], e2[3];
      for (int c = 0; c < 3; c++) {
        e1[c] = p[1][c] - p[0][c];
        e2[c] = p[2][c] - p[0][c];
      }
      gn[0] = e1[1] * e2[2] - e1[2] * e2[1];
      gn[1] = e1[2] * e2[0] - e1[0] * e2[2];
      gn[2] = e1[0] * e2[1] - e1[1] * e2[0];
      const double *A = wpos + (size_t)a * 3, *B = wpos + (size_t)b * 3;
      double mid[3];
      for (int c = 0; c < 3; c++)
        mid[c] = 0.5 * (A[c] + B[c]);

      for (int s = 0; s < nsrc; s++) {
        const hz_poly *S = &ps->p[src[s]];
        /* СИЛУЭТ: у двух граней ребра РАЗНЫЙ знак обращённости к источнику.
         * Ребро без пары — силуэт всегда (край дыры в меше). */
        double v1[3];
        for (int c = 0; c < 3; c++)
          v1[c] = S->org[c] - mid[c];
        double f1 = gn[0] * v1[0] + gn[1] * v1[1] + gn[2] * v1[2];
        int sil;
        if (t2 < 0) {
          sil = 1;
        } else {
          double q[3][3], f1b[3], f2b[3], g2[3];
          hz_obj_tri(m, t2, q);
          for (int c = 0; c < 3; c++) {
            f1b[c] = q[1][c] - q[0][c];
            f2b[c] = q[2][c] - q[0][c];
          }
          g2[0] = f1b[1] * f2b[2] - f1b[2] * f2b[1];
          g2[1] = f1b[2] * f2b[0] - f1b[0] * f2b[2];
          g2[2] = f1b[0] * f2b[1] - f1b[1] * f2b[0];
          double f2 = g2[0] * v1[0] + g2[1] * v1[1] + g2[2] * v1[2];
          sil = ((f1 > 0.0) != (f2 > 0.0));
        }
        if (!sil) continue;
        ws->nsil++;

        /* ПРАВИЛО A и C: угловой размер источника с ребра задаёт длину клина.
         * `α ≈ sqrt(площадь)/расстояние`; клин живёт, пока `α·d < ε`. */
        double dist = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
        if (!(dist > 0.0)) continue;
        double alpha = sqrt(S->mom[0]) / dist;
        if (!(alpha > 0.0)) continue;
        double rA = cfg->eps / alpha;

        /* ПРАВИЛО B: контраст. Прямой вклад источника у ребра против средней
         * облучённости сцены. Нормаль берём НА источник — верхняя оценка. */
        double nn[3];
        for (int c = 0; c < 3; c++)
          nn[c] = v1[c] / dist;
        double Ed = hz_direct_unocc(ps, src[s], Le[src[s]], mid, nn);
        double prio = Ed / Emean;
        if (!(prio > cfg->tol)) {
          ws->nreject_b++;
          continue;
        }
        /* ПРАВИЛО C: радиус, за которым вклад падает ниже tol по спаду 1/r².
         * Берём МЕНЬШИЙ из двух — A и B. */
        double rB = dist * sqrt(prio / cfg->tol);
        double r = (rA < rB) ? rA : rB;
        if (!(r > 0.0)) {
          ws->nreject_a++;
          continue;
        }

        /* Плоскость клина: через ребро и центр источника. */
        double ab[3], as[3], pn[3];
        for (int c = 0; c < 3; c++) {
          ab[c] = B[c] - A[c];
          as[c] = S->org[c] - A[c];
        }
        pn[0] = ab[1] * as[2] - ab[2] * as[1];
        pn[1] = ab[2] * as[0] - ab[0] * as[2];
        pn[2] = ab[0] * as[1] - ab[1] * as[0];
        double pl = sqrt(pn[0] * pn[0] + pn[1] * pn[1] + pn[2] * pn[2]);
        if (!(pl > 0.0)) continue;
        if (ws->n >= ws->cap) {
          int32_t nc = (ws->cap > 0) ? ws->cap * 2 : 4096;
          hz_wedge *nwd = realloc(ws->w, (size_t)nc * sizeof *nwd);
          if (nwd == NULL) {
            free(h.key);
            free(h.tri);
            free(wid);
            free(wpos);
            return 2;
          }
          ws->w = nwd;
          ws->cap = nc;
        }
        hz_wedge *W = &ws->w[ws->n++];
        for (int c = 0; c < 3; c++) {
          W->n[c] = pn[c] / pl;
          W->c[c] = mid[c];
        }
        W->off = W->n[0] * A[0] + W->n[1] * A[1] + W->n[2] * A[2];
        W->r = r;
        W->prio = prio;
        W->src = src[s];
      }
    }
  }
  qsort(ws->w, (size_t)ws->n, sizeof *ws->w, cmp_wedge);
  free(h.key);
  free(h.tri);
  free(wid);
  free(wpos);
  return 0;
}

/* --- применение: измельчение супа ------------------------------------------- */

typedef struct {
  double p[3][3], nrm[3][3];
  int has_n;
} frag;

/* Точка пересечения ребра с плоскостью, посчитанная ДЕТЕРМИНИРОВАННО от
 * лексикографически меньшего конца. Без этого два треугольника, делящие
 * ребро, дадут РАЗНЫЕ биты одной и той же точки, сварка их не сольёт, и в
 * крае появится щель — которую теорема Грина тут же и поймает. */
static void edge_split(const double *A, const double *B, const double *na, const double *nb,
                       const double n[3], double off, double out[3], double outn[3]) {
  const double *P = A, *Q = B, *PN = na, *QN = nb;
  int swap = 0;
  for (int a = 0; a < 3; a++) {
    if (A[a] < B[a]) break;
    if (A[a] > B[a]) {
      swap = 1;
      break;
    }
  }
  if (swap) {
    P = B;
    Q = A;
    PN = nb;
    QN = na;
  }
  double dp = n[0] * P[0] + n[1] * P[1] + n[2] * P[2] - off;
  double dq = n[0] * Q[0] + n[1] * Q[1] + n[2] * Q[2] - off;
  double s = dp / (dp - dq);
  for (int a = 0; a < 3; a++) {
    out[a] = P[a] + s * (Q[a] - P[a]);
    outn[a] = PN[a] + s * (QN[a] - PN[a]);
  }
}

/* Делит фрагмент плоскостью; кладёт куски в out, возвращает их число.
 * `sideout` получает знак стороны каждого куска. */
static int split_frag(const frag *f, const double n[3], double off, frag *out, int *sideout) {
  double d[3];
  int pos = 0, neg = 0;
  for (int i = 0; i < 3; i++) {
    d[i] = n[0] * f->p[i][0] + n[1] * f->p[i][1] + n[2] * f->p[i][2] - off;
    if (d[i] > 0.0)
      pos++;
    else
      neg++;
  }
  if (pos == 0 || neg == 0) {
    out[0] = *f;
    sideout[0] = (pos > 0) ? 1 : 0;
    return 1;
  }
  /* Одна вершина по одну сторону, две по другую. Находим одиночку. */
  int lone = 0;
  for (int i = 0; i < 3; i++)
    if ((d[i] > 0.0) != (d[(i + 1) % 3] > 0.0) && (d[i] > 0.0) != (d[(i + 2) % 3] > 0.0)) lone = i;
  int i1 = (lone + 1) % 3, i2 = (lone + 2) % 3;
  double x1[3], x2[3], n1[3], n2[3];
  edge_split(f->p[lone], f->p[i1], f->nrm[lone], f->nrm[i1], n, off, x1, n1);
  edge_split(f->p[lone], f->p[i2], f->nrm[lone], f->nrm[i2], n, off, x2, n2);
  int side_lone = (d[lone] > 0.0) ? 1 : 0;
  /* кусок с одиночкой */
  for (int a = 0; a < 3; a++) {
    out[0].p[0][a] = f->p[lone][a];
    out[0].nrm[0][a] = f->nrm[lone][a];
    out[0].p[1][a] = x1[a];
    out[0].nrm[1][a] = n1[a];
    out[0].p[2][a] = x2[a];
    out[0].nrm[2][a] = n2[a];
  }
  out[0].has_n = f->has_n;
  sideout[0] = side_lone;
  /* два куска по другую сторону — порядок вершин сохраняет ориентацию */
  for (int a = 0; a < 3; a++) {
    out[1].p[0][a] = x1[a];
    out[1].nrm[0][a] = n1[a];
    out[1].p[1][a] = f->p[i1][a];
    out[1].nrm[1][a] = f->nrm[i1][a];
    out[1].p[2][a] = f->p[i2][a];
    out[1].nrm[2][a] = f->nrm[i2][a];
    out[2].p[0][a] = x1[a];
    out[2].nrm[0][a] = n1[a];
    out[2].p[1][a] = f->p[i2][a];
    out[2].nrm[1][a] = f->nrm[i2][a];
    out[2].p[2][a] = x2[a];
    out[2].nrm[2][a] = n2[a];
  }
  out[1].has_n = out[2].has_n = f->has_n;
  sideout[1] = sideout[2] = 1 - side_lone;
  return 3;
}

int hz_cut_apply(hz_objmesh *mo, hz_pseglist *so, const hz_objmesh *mi, const hz_pseglist *si,
                 const hz_wedges *ws, int32_t ncut) {
  memset(mo, 0, sizeof *mo);
  memset(so, 0, sizeof *so);
  if (ncut > ws->n) ncut = ws->n;

  /* --- какие клинья к какому участку. Габарит участка против шара клина. --- */
  const int32_t nseg = si->nseg;
  double *sb = malloc((size_t)nseg * 6 * sizeof *sb);
  int32_t *wl = malloc((size_t)nseg * PCUT_MAXW * sizeof *wl);
  int32_t *wn = calloc((size_t)nseg, sizeof *wn);
  if (sb == NULL || wl == NULL || wn == NULL) {
    free(sb);
    free(wl);
    free(wn);
    return 2;
  }
  for (int32_t r = 0; r < nseg; r++) {
    sb[(size_t)r * 6 + 0] = sb[(size_t)r * 6 + 1] = sb[(size_t)r * 6 + 2] = 1e300;
    sb[(size_t)r * 6 + 3] = sb[(size_t)r * 6 + 4] = sb[(size_t)r * 6 + 5] = -1e300;
  }
  for (int32_t t = 0; t < mi->nt; t++) {
    int32_t r = si->label[t];
    double p[3][3];
    hz_obj_tri(mi, t, p);
    for (int i = 0; i < 3; i++)
      for (int a = 0; a < 3; a++) {
        if (p[i][a] < sb[(size_t)r * 6 + (size_t)a]) sb[(size_t)r * 6 + (size_t)a] = p[i][a];
        if (p[i][a] > sb[(size_t)r * 6 + 3 + (size_t)a])
          sb[(size_t)r * 6 + 3 + (size_t)a] = p[i][a];
      }
  }
  /* Клинья уже отсортированы по контрасту, поэтому первые попавшие участку —
   * и есть самые важные (Т2: бюджет с приоритетом по контрасту). */
  for (int32_t iw = 0; iw < ncut; iw++) {
    const hz_wedge *W = &ws->w[iw];
    for (int32_t r = 0; r < nseg; r++) {
      if (wn[r] >= PCUT_MAXW) continue;
      double dd = 0.0;
      for (int a = 0; a < 3; a++) {
        double c = W->c[a];
        double lo = sb[(size_t)r * 6 + (size_t)a], hi = sb[(size_t)r * 6 + 3 + (size_t)a];
        double q = (c < lo) ? lo - c : ((c > hi) ? c - hi : 0.0);
        dd += q * q;
      }
      if (dd > W->r * W->r) continue;
      /* Плоскость обязана участок ПЕРЕСЕКАТЬ, иначе разрез ничего не делит. */
      int pos = 0, neg = 0;
      for (int k = 0; k < 8; k++) {
        double c[3] = {(k & 1) ? sb[(size_t)r * 6 + 3] : sb[(size_t)r * 6 + 0],
                       (k & 2) ? sb[(size_t)r * 6 + 4] : sb[(size_t)r * 6 + 1],
                       (k & 4) ? sb[(size_t)r * 6 + 5] : sb[(size_t)r * 6 + 2]};
        if (W->n[0] * c[0] + W->n[1] * c[1] + W->n[2] * c[2] - W->off > 0.0)
          pos++;
        else
          neg++;
      }
      if (pos == 0 || neg == 0) continue;
      wl[(size_t)r * PCUT_MAXW + (size_t)wn[r]] = iw;
      wn[r]++;
    }
  }

  /* --- измельчение --- */
  int64_t cap = mi->nt * 2 + 1024;
  mo->f = malloc((size_t)cap * 3 * sizeof *mo->f);
  mo->fn = malloc((size_t)cap * 3 * sizeof *mo->fn);
  mo->fm = malloc((size_t)cap * sizeof *mo->fm);
  mo->v = malloc((size_t)cap * 9 * sizeof *mo->v);
  mo->vn = malloc((size_t)cap * 9 * sizeof *mo->vn);
  int32_t *lab = malloc((size_t)cap * sizeof *lab);
  int32_t *sig = malloc((size_t)cap * sizeof *sig);
  if (mo->f == NULL || mo->fn == NULL || mo->fm == NULL || mo->v == NULL || mo->vn == NULL ||
      lab == NULL || sig == NULL) {
    free(lab);
    free(sig);
    free(sb);
    free(wl);
    free(wn);
    hz_obj_free(mo);
    return 2;
  }
  int64_t nt = 0;
  frag cur[PCUT_MAXFRAG], nxt[PCUT_MAXFRAG];
  int curside[PCUT_MAXFRAG], nxtside[PCUT_MAXFRAG];

  for (int32_t t = 0; t < mi->nt; t++) {
    int32_t r = si->label[t];
    int nc = 1;
    hz_obj_tri(mi, t, cur[0].p);
    cur[0].has_n = 0;
    for (int i = 0; i < 3; i++) {
      int32_t ni = (mi->vn != NULL) ? mi->fn[(size_t)t * 3 + (size_t)i] : -1;
      for (int a = 0; a < 3; a++)
        cur[0].nrm[i][a] = (ni >= 0) ? mi->vn[(size_t)ni * 3 + (size_t)a] : 0.0;
      if (ni >= 0) cur[0].has_n = 1;
    }
    curside[0] = 0;
    for (int j = 0; j < wn[r]; j++) {
      const hz_wedge *W = &ws->w[wl[(size_t)r * PCUT_MAXW + (size_t)j]];
      int nn2 = 0;
      for (int c = 0; c < nc; c++) {
        frag out[3];
        int so2[3];
        int k = split_frag(&cur[c], W->n, W->off, out, so2);
        for (int q = 0; q < k && nn2 < PCUT_MAXFRAG; q++) {
          nxt[nn2] = out[q];
          nxtside[nn2] = curside[c] | (so2[q] << j);
          nn2++;
        }
      }
      memcpy(cur, nxt, (size_t)nn2 * sizeof *cur);
      memcpy(curside, nxtside, (size_t)nn2 * sizeof *curside);
      nc = nn2;
    }
    for (int c = 0; c < nc; c++) {
      /* РОСТ ПО ОДНОМУ УКАЗАТЕЛЮ ЗА РАЗ, с немедленной фиксацией. Пакетный
       * `realloc` семи массивов с общей проверкой в конце оставляет за собой
       * освобождённые указатели, если упал третий из семи, — gcc-analyzer это
       * и назвал use-after-free. Здесь каждый успешный `realloc` сразу
       * записывается на место, поэтому структура в любой момент согласована. */
      if (nt >= cap) {
        int64_t nc2 = cap * 2;
        void *q;
        int bad = 0;
        if ((q = realloc(mo->v, (size_t)nc2 * 9 * sizeof *mo->v)) != NULL)
          mo->v = q;
        else
          bad = 1;
        if ((q = realloc(mo->vn, (size_t)nc2 * 9 * sizeof *mo->vn)) != NULL)
          mo->vn = q;
        else
          bad = 1;
        if ((q = realloc(mo->f, (size_t)nc2 * 3 * sizeof *mo->f)) != NULL)
          mo->f = q;
        else
          bad = 1;
        if ((q = realloc(mo->fn, (size_t)nc2 * 3 * sizeof *mo->fn)) != NULL)
          mo->fn = q;
        else
          bad = 1;
        if ((q = realloc(mo->fm, (size_t)nc2 * sizeof *mo->fm)) != NULL)
          mo->fm = q;
        else
          bad = 1;
        if ((q = realloc(lab, (size_t)nc2 * sizeof *lab)) != NULL)
          lab = q;
        else
          bad = 1;
        if ((q = realloc(sig, (size_t)nc2 * sizeof *sig)) != NULL)
          sig = q;
        else
          bad = 1;
        if (bad) {
          free(lab);
          free(sig);
          free(sb);
          free(wl);
          free(wn);
          hz_obj_free(mo);
          return 2;
        }
        cap = nc2;
      }
      /* Вершины дублируются НАМЕРЕННО: сварка по положению в `hz_poly_build`
       * сольёт их обратно, а собственной таблицы вершин здесь не заводится. */
      for (int i = 0; i < 3; i++) {
        for (int a = 0; a < 3; a++) {
          mo->v[(size_t)(nt * 3 + i) * 3 + (size_t)a] = cur[c].p[i][a];
          mo->vn[(size_t)(nt * 3 + i) * 3 + (size_t)a] = cur[c].nrm[i][a];
        }
        mo->f[(size_t)nt * 3 + (size_t)i] = (int32_t)(nt * 3 + i);
        mo->fn[(size_t)nt * 3 + (size_t)i] = cur[c].has_n ? (int32_t)(nt * 3 + i) : -1;
      }
      mo->fm[nt] = mi->fm[t];
      lab[nt] = r;
      sig[nt] = curside[c];
      nt++;
    }
  }
  mo->nt = (int32_t)nt;
  mo->nv = (int32_t)(nt * 3);
  mo->nvn = (int32_t)(nt * 3);
  mo->nmtl = mi->nmtl;
  mo->mtl = malloc((size_t)mi->nmtl * sizeof *mo->mtl);
  if (mo->mtl == NULL) {
    free(lab);
    free(sig);
    free(sb);
    free(wl);
    free(wn);
    hz_obj_free(mo);
    return 2;
  }
  memcpy(mo->mtl, mi->mtl, (size_t)mi->nmtl * sizeof *mo->mtl);
  for (int a = 0; a < 3; a++) {
    mo->lo[a] = mi->lo[a];
    mo->hi[a] = mi->hi[a];
  }

  /* --- новая разметка: (участок, подпись) -> номер --- */
  so->label = malloc((size_t)nt * sizeof *so->label);
  so->seg = malloc((size_t)nt * sizeof *so->seg);
  int32_t *base = malloc((size_t)nseg * (1 << PCUT_MAXW) * sizeof *base);
  if (so->label == NULL || so->seg == NULL || base == NULL) {
    free(base);
    free(lab);
    free(sig);
    free(sb);
    free(wl);
    free(wn);
    hz_obj_free(mo);
    hz_seg_free(so);
    return 2;
  }
  for (int64_t i = 0; i < (int64_t)nseg * (1 << PCUT_MAXW); i++)
    base[i] = -1;
  int32_t nnew = 0;
  for (int64_t t = 0; t < nt; t++) {
    int64_t key = (int64_t)lab[t] * (1 << PCUT_MAXW) + sig[t];
    if (base[key] < 0) {
      base[key] = nnew;
      /* ПЛОСКОСТЬ НЕ МЕНЯЕТСЯ: разрез делит участок, а не гнёт его. */
      so->seg[nnew] = si->seg[lab[t]];
      so->seg[nnew].area = 0.0;
      so->seg[nnew].ntri = 0;
      so->seg[nnew].dmax = 0.0;
      nnew++;
    }
    int32_t g = base[key];
    so->label[t] = g;
    so->seg[g].ntri++;
    so->seg[g].area += hz_obj_tri_area(mo, (int32_t)t);
    double p[3][3];
    hz_obj_tri(mo, (int32_t)t, p);
    for (int i = 0; i < 3; i++) {
      double dv = fabs(p[i][0] * so->seg[g].n[0] + p[i][1] * so->seg[g].n[1] +
                       p[i][2] * so->seg[g].n[2] - so->seg[g].off);
      if (dv > so->seg[g].dmax) so->seg[g].dmax = dv;
    }
  }
  so->nseg = nnew;
  so->delta = si->delta;
  free(base);
  free(lab);
  free(sig);
  free(sb);
  free(wl);
  free(wn);
  return 0;
}

/* --- дробление регулярной сеткой ПЛОСКОСТЯМИ (разбор — в `pcut.h`) --------- */

#define PCUT_GRIDFRAG 8192

int hz_cut_grid(hz_objmesh *mo, hz_pseglist *so, const hz_objmesh *mi, const hz_pseglist *si,
                double L, int64_t *nover) {
  memset(mo, 0, sizeof *mo);
  memset(so, 0, sizeof *so);
  if (!(L > 0.0)) return 1;
  if (nover != NULL) *nover = 0;

  int64_t cap = (int64_t)mi->nt * 2 + 4096;
  mo->f = malloc((size_t)cap * 3 * sizeof *mo->f);
  mo->fn = malloc((size_t)cap * 3 * sizeof *mo->fn);
  mo->fm = malloc((size_t)cap * sizeof *mo->fm);
  mo->v = malloc((size_t)cap * 9 * sizeof *mo->v);
  mo->vn = malloc((size_t)cap * 9 * sizeof *mo->vn);
  int32_t *lab = malloc((size_t)cap * sizeof *lab);
  int64_t *cellid = malloc((size_t)cap * sizeof *cellid);
  frag *cur = malloc(PCUT_GRIDFRAG * sizeof *cur), *nxt = malloc(PCUT_GRIDFRAG * sizeof *nxt);
  if (mo->f == NULL || mo->fn == NULL || mo->fm == NULL || mo->v == NULL || mo->vn == NULL ||
      lab == NULL || cellid == NULL || cur == NULL || nxt == NULL) {
    free(lab);
    free(cellid);
    free(cur);
    free(nxt);
    hz_obj_free(mo);
    return 2;
  }

  const int32_t nx = (int32_t)((mi->hi[0] - mi->lo[0]) / L) + 2;
  const int32_t ny = (int32_t)((mi->hi[1] - mi->lo[1]) / L) + 2;
  int64_t nt = 0;

  for (int32_t t = 0; t < mi->nt; t++) {
    int nc = 1;
    hz_obj_tri(mi, t, cur[0].p);
    cur[0].has_n = 0;
    for (int i = 0; i < 3; i++) {
      int32_t ni = (mi->vn != NULL) ? mi->fn[(size_t)t * 3 + (size_t)i] : -1;
      for (int a = 0; a < 3; a++)
        cur[0].nrm[i][a] = (ni >= 0) ? mi->vn[(size_t)ni * 3 + (size_t)a] : 0.0;
      if (ni >= 0) cur[0].has_n = 1;
    }
    /* Режем по каждой оси всеми плоскостями сетки, попавшими в габарит. */
    for (int a = 0; a < 3 && nc > 0; a++) {
      double lo = 1e300, hi = -1e300;
      for (int c = 0; c < nc; c++)
        for (int i = 0; i < 3; i++) {
          if (cur[c].p[i][a] < lo) lo = cur[c].p[i][a];
          if (cur[c].p[i][a] > hi) hi = cur[c].p[i][a];
        }
      int64_t k0 = (int64_t)floor((lo - mi->lo[a]) / L) + 1;
      int64_t k1 = (int64_t)floor((hi - mi->lo[a]) / L);
      double n[3] = {0.0, 0.0, 0.0};
      n[a] = 1.0;
      for (int64_t k = k0; k <= k1; k++) {
        double off = mi->lo[a] + (double)k * L;
        int nn2 = 0;
        for (int c = 0; c < nc; c++) {
          frag out[3];
          int so2[3];
          int kk = split_frag(&cur[c], n, off, out, so2);
          for (int q = 0; q < kk; q++) {
            if (nn2 >= PCUT_GRIDFRAG) {
              if (nover != NULL) (*nover)++;
              continue;
            }
            nxt[nn2++] = out[q];
          }
        }
        memcpy(cur, nxt, (size_t)nn2 * sizeof *cur);
        nc = nn2;
      }
    }
    for (int c = 0; c < nc; c++) {
      if (nt >= cap) {
        int64_t nc2 = cap * 2;
        void *q;
        int bad = 0;
        if ((q = realloc(mo->v, (size_t)nc2 * 9 * sizeof *mo->v)) != NULL)
          mo->v = q;
        else
          bad = 1;
        if ((q = realloc(mo->vn, (size_t)nc2 * 9 * sizeof *mo->vn)) != NULL)
          mo->vn = q;
        else
          bad = 1;
        if ((q = realloc(mo->f, (size_t)nc2 * 3 * sizeof *mo->f)) != NULL)
          mo->f = q;
        else
          bad = 1;
        if ((q = realloc(mo->fn, (size_t)nc2 * 3 * sizeof *mo->fn)) != NULL)
          mo->fn = q;
        else
          bad = 1;
        if ((q = realloc(mo->fm, (size_t)nc2 * sizeof *mo->fm)) != NULL)
          mo->fm = q;
        else
          bad = 1;
        if ((q = realloc(lab, (size_t)nc2 * sizeof *lab)) != NULL)
          lab = q;
        else
          bad = 1;
        if ((q = realloc(cellid, (size_t)nc2 * sizeof *cellid)) != NULL)
          cellid = q;
        else
          bad = 1;
        if (bad) {
          free(lab);
          free(cellid);
          free(cur);
          free(nxt);
          hz_obj_free(mo);
          return 2;
        }
        cap = nc2;
      }
      double cc[3] = {0.0, 0.0, 0.0};
      for (int i = 0; i < 3; i++)
        for (int a = 0; a < 3; a++)
          cc[a] += cur[c].p[i][a] / 3.0;
      int64_t ix = (int64_t)((cc[0] - mi->lo[0]) / L);
      int64_t iy = (int64_t)((cc[1] - mi->lo[1]) / L);
      int64_t iz = (int64_t)((cc[2] - mi->lo[2]) / L);
      for (int i = 0; i < 3; i++) {
        for (int a = 0; a < 3; a++) {
          mo->v[(size_t)(nt * 3 + i) * 3 + (size_t)a] = cur[c].p[i][a];
          mo->vn[(size_t)(nt * 3 + i) * 3 + (size_t)a] = cur[c].nrm[i][a];
        }
        mo->f[(size_t)nt * 3 + (size_t)i] = (int32_t)(nt * 3 + i);
        mo->fn[(size_t)nt * 3 + (size_t)i] = cur[c].has_n ? (int32_t)(nt * 3 + i) : -1;
      }
      mo->fm[nt] = mi->fm[t];
      lab[nt] = si->label[t];
      cellid[nt] = (iz * ny + iy) * nx + ix;
      nt++;
    }
  }
  free(cur);
  free(nxt);
  mo->nt = (int32_t)nt;
  mo->nv = (int32_t)(nt * 3);
  mo->nvn = (int32_t)(nt * 3);
  mo->nmtl = mi->nmtl;
  mo->mtl = malloc((size_t)mi->nmtl * sizeof *mo->mtl);
  if (mo->mtl == NULL) {
    free(lab);
    free(cellid);
    hz_obj_free(mo);
    return 2;
  }
  memcpy(mo->mtl, mi->mtl, (size_t)mi->nmtl * sizeof *mo->mtl);
  for (int a = 0; a < 3; a++) {
    mo->lo[a] = mi->lo[a];
    mo->hi[a] = mi->hi[a];
  }

  /* Метка = (участок, ячейка). Ключей много, поэтому ХЕШ, а не плотный массив:
   * плотный был бы nseg × ncell и на мелкой сетке не поместился бы. */
  so->label = malloc((size_t)nt * sizeof *so->label);
  so->seg = malloc((size_t)nt * sizeof *so->seg);
  int64_t hn = 4;
  while (hn < 4 * nt)
    hn *= 2;
  uint64_t *hk = calloc((size_t)hn, sizeof *hk);
  int32_t *hv = malloc((size_t)hn * sizeof *hv);
  if (so->label == NULL || so->seg == NULL || hk == NULL || hv == NULL) {
    free(hk);
    free(hv);
    free(lab);
    free(cellid);
    hz_obj_free(mo);
    hz_seg_free(so);
    return 2;
  }
  int32_t nn = 0;
  for (int64_t i = 0; i < nt; i++) {
    uint64_t key = ((uint64_t)lab[i] << 40) ^ (uint64_t)cellid[i];
    uint64_t kk = key + 1, s = hmix(kk) & (uint64_t)(hn - 1);
    while (hk[s] != 0 && hk[s] != kk)
      s = (s + 1) & (uint64_t)(hn - 1);
    if (hk[s] == 0) {
      hk[s] = kk;
      hv[s] = nn;
      so->seg[nn] = si->seg[lab[i]];
      so->seg[nn].area = 0.0;
      so->seg[nn].ntri = 0;
      so->seg[nn].dmax = 0.0;
      nn++;
    }
    int32_t g = hv[s];
    so->label[i] = g;
    so->seg[g].ntri++;
    so->seg[g].area += hz_obj_tri_area(mo, (int32_t)i);
    double p[3][3];
    hz_obj_tri(mo, (int32_t)i, p);
    for (int j = 0; j < 3; j++) {
      double dv = fabs(p[j][0] * so->seg[g].n[0] + p[j][1] * so->seg[g].n[1] +
                       p[j][2] * so->seg[g].n[2] - so->seg[g].off);
      if (dv > so->seg[g].dmax) so->seg[g].dmax = dv;
    }
  }
  so->nseg = nn;
  so->delta = si->delta;
  free(hk);
  free(hv);
  free(lab);
  free(cellid);
  return 0;
}
