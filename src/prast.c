/* Параллельная проекция полигонов. Разбор и оговорки — в `prast.h`. */

#include "prast.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

int hz_pview_make(hz_pview *v, const double w[3], const double lo[3], const double hi[3],
                  double h) {
  memset(v, 0, sizeof *v);
  double wn = sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
  if (!(wn > 0.0) || !(h > 0.0)) return 1;
  for (int a = 0; a < 3; a++)
    v->w[a] = w[a] / wn;
  v->h = h;
  /* Ось, наименее совпадающая с ω, — выбор детерминированный и без порога. */
  int ax = 0;
  for (int a = 1; a < 3; a++)
    if (fabs(v->w[a]) < fabs(v->w[ax])) ax = a;
  double e0[3] = {0.0, 0.0, 0.0};
  e0[ax] = 1.0;
  double d = e0[0] * v->w[0] + e0[1] * v->w[1] + e0[2] * v->w[2];
  for (int a = 0; a < 3; a++)
    v->ea[a] = e0[a] - d * v->w[a];
  double an = sqrt(v->ea[0] * v->ea[0] + v->ea[1] * v->ea[1] + v->ea[2] * v->ea[2]);
  if (!(an > 0.0)) return 1;
  for (int a = 0; a < 3; a++)
    v->ea[a] /= an;
  v->eb[0] = v->w[1] * v->ea[2] - v->w[2] * v->ea[1];
  v->eb[1] = v->w[2] * v->ea[0] - v->w[0] * v->ea[2];
  v->eb[2] = v->w[0] * v->ea[1] - v->w[1] * v->ea[0];

  /* Экран накрывает проекцию габарита сцены: по восьми углам, а не по центру. */
  double amin = 1e300, amax = -1e300, bmin = 1e300, bmax = -1e300;
  for (int k = 0; k < 8; k++) {
    double c[3] = {(k & 1) ? hi[0] : lo[0], (k & 2) ? hi[1] : lo[1], (k & 4) ? hi[2] : lo[2]};
    double ca = c[0] * v->ea[0] + c[1] * v->ea[1] + c[2] * v->ea[2];
    double cb = c[0] * v->eb[0] + c[1] * v->eb[1] + c[2] * v->eb[2];
    if (ca < amin) amin = ca;
    if (ca > amax) amax = ca;
    if (cb < bmin) bmin = cb;
    if (cb > bmax) bmax = cb;
  }
  double na = ceil((amax - amin) / h) + 2.0, nb = ceil((bmax - bmin) / h) + 2.0;
  if (!(na < 2.0e9) || !(nb < 2.0e9) || na * nb > 2.0e9) return 2;
  v->W = (int32_t)na;
  v->H = (int32_t)nb;
  for (int a = 0; a < 3; a++)
    v->o[a] = amin * v->ea[a] + bmin * v->eb[a];
  return 0;
}

void hz_pview_origin(const hz_pview *v, int32_t i, int32_t j, double r[3]) {
  double sa = ((double)i + 0.5) * v->h, sb = ((double)j + 0.5) * v->h;
  for (int a = 0; a < 3; a++)
    r[a] = v->o[a] + sa * v->ea[a] + sb * v->eb[a];
}

/* --- заливка: полигоны -> полосы -------------------------------------------- */

typedef struct {
  double a0, a1, b0, b1; /* ребро в ЭКРАННЫХ пикселях, направленное */
} pe_scr;

typedef struct {
  double x;
  int dir;
} pe_cross;

/* ВСТАВКАМИ, А НЕ `qsort`, И ЭТО ЗАМЕР, А НЕ ВКУС. Пересечений в строке
 * единицы (у выпуклого края — ровно два), а `qsort` платит вызовом через
 * указатель на каждое сравнение. В первом прогоне Ш3 на этом стояло около
 * двух третей времени растеризации: 960 тысяч вызовов `qsort` на отскок при
 * средней длине списка 4. Вставки на такой длине оптимальны и без вызовов. */
static void sort_cross(pe_cross *c, int32_t n) {
  for (int32_t i = 1; i < n; i++) {
    pe_cross v = c[i];
    int32_t j = i - 1;
    while (j >= 0 && c[j].x > v.x) {
      c[j + 1] = c[j];
      j--;
    }
    c[j + 1] = v;
  }
}

int hz_prast_spans(hz_span **spp, int64_t *nsp, int64_t *cap, const hz_pview *v,
                   const hz_polyset *ps, const unsigned char *live, hz_rstats *st) {
  memset(st, 0, sizeof *st);
  *nsp = 0;
  const int32_t W = v->W, H = v->H;

  int32_t maxe = 4;
  for (int32_t k = 0; k < ps->np; k++) {
    int32_t e =
        (ps->p[k].nloop > 0) ? ps->loop[ps->p[k].l0 + ps->p[k].nloop] - ps->loop[ps->p[k].l0] : 0;
    if (e > maxe) maxe = e;
  }
  /* ТАБЛИЦА РЁБЕР ПО СТРОКАМ, а не перебор всех рёбер на каждой строке.
   * Наивная заливка стоит `O(рёбра × строки)` и штрафует ИМЕННО крупный
   * полигон — то есть ровно то, ради чего вся модель и затевалась. Замерено на
   * зале при h = 2 см: 333 полигона (по 6.8 м периметра) растеризуются ДОЛЬШЕ,
   * чем 7 244 (по 0.75 м), — 0.205 с против 0.126 при том, что фрагментов у
   * первых МЕНЬШЕ. Таблица убирает этот член: `O(рёбра + Σ активных)`. */
  pe_scr *ed = malloc((size_t)maxe * sizeof *ed);
  pe_cross *cr = malloc((size_t)maxe * sizeof *cr);
  int32_t *elo = malloc((size_t)maxe * sizeof *elo);
  int32_t *ehi = malloc((size_t)maxe * sizeof *ehi);
  int32_t *ebkt = malloc((size_t)maxe * sizeof *ebkt);
  int32_t *act = malloc((size_t)maxe * sizeof *act);
  int32_t *bcnt = calloc((size_t)H + 2, sizeof *bcnt);
  if (ed == NULL || cr == NULL || elo == NULL || ehi == NULL || ebkt == NULL || act == NULL ||
      bcnt == NULL) {
    free(ed);
    free(cr);
    free(elo);
    free(ehi);
    free(ebkt);
    free(act);
    free(bcnt);
    return 2;
  }

  for (int32_t k = 0; k < ps->np; k++) {
    const hz_poly *P = &ps->p[k];
    if (P->nloop == 0) continue;
    /* ПОЛИГОНА БОЛЬШЕ НЕТ (§124): разрушенная геометрия не рисуется вовсе —
     * значит и не заслоняет, и не светит. `live == NULL` — все живы. */
    if (live != NULL && !live[k]) continue;
    double nw = P->n[0] * v->w[0] + P->n[1] * v->w[1] + P->n[2] * v->w[2];
    /* Точный ноль: полигон РЕБРОМ к направлению — ни излучает вдоль него, ни
     * принимает, а глубина была бы делением на ноль. Порога нет и не нужно:
     * при малом, но ненулевом `nw` проекция сама вырождается в щепку. */
    if (!(nw < 0.0) && !(nw > 0.0)) {
      st->nskip++;
      continue;
    }

    double ou[3];
    for (int a = 0; a < 3; a++)
      ou[a] = P->org[a] - v->o[a];
    double A0 = (ou[0] * v->ea[0] + ou[1] * v->ea[1] + ou[2] * v->ea[2]) / v->h;
    double B0 = (ou[0] * v->eb[0] + ou[1] * v->eb[1] + ou[2] * v->eb[2]) / v->h;
    double Au = (P->eu[0] * v->ea[0] + P->eu[1] * v->ea[1] + P->eu[2] * v->ea[2]) / v->h;
    double Bu = (P->eu[0] * v->eb[0] + P->eu[1] * v->eb[1] + P->eu[2] * v->eb[2]) / v->h;
    double Av = (P->ev[0] * v->ea[0] + P->ev[1] * v->ea[1] + P->ev[2] * v->ea[2]) / v->h;
    double Bv = (P->ev[0] * v->eb[0] + P->ev[1] * v->eb[1] + P->ev[2] * v->eb[2]) / v->h;

    double no = P->n[0] * v->o[0] + P->n[1] * v->o[1] + P->n[2] * v->o[2];
    double nea = P->n[0] * v->ea[0] + P->n[1] * v->ea[1] + P->n[2] * v->ea[2];
    double neb = P->n[0] * v->eb[0] + P->n[1] * v->eb[1] + P->n[2] * v->eb[2];
    double T0 = (P->off - no - 0.5 * v->h * nea - 0.5 * v->h * neb) / nw;
    double Ti = -v->h * nea / nw, Tj = -v->h * neb / nw;

    int32_t ne = 0;
    double tmin = 1e300, tmax = -1e300;
    for (int32_t l = P->l0; l < P->l0 + P->nloop; l++) {
      int32_t b = ps->loop[l], e = ps->loop[l + 1], n = e - b;
      for (int32_t i = 0; i < n; i++) {
        const double *X = ps->bv + (size_t)(b + i) * 2;
        const double *Y = ps->bv + (size_t)(b + (i + 1) % n) * 2;
        ed[ne].a0 = A0 + Au * X[0] + Av * X[1];
        ed[ne].b0 = B0 + Bu * X[0] + Bv * X[1];
        ed[ne].a1 = A0 + Au * Y[0] + Av * Y[1];
        ed[ne].b1 = B0 + Bu * Y[0] + Bv * Y[1];
        if (ed[ne].b0 < tmin) tmin = ed[ne].b0;
        if (ed[ne].b0 > tmax) tmax = ed[ne].b0;
        ne++;
      }
    }
    if (ne < 3) continue;
    int32_t j0 = (int32_t)floor(tmin - 0.5), j1 = (int32_t)ceil(tmax + 0.5);
    if (j0 < 0) j0 = 0;
    if (j1 > H - 1) j1 = H - 1;

    /* Строчный диапазон ребра. Пересечение с центром строки `yc = j + 0.5`
     * бывает ровно при `yc ∈ [min, max)` — то же полуоткрытое условие, что и
     * ниже, только решённое относительно `j`. */
    for (int32_t i = 0; i < ne; i++) {
      double y0 = ed[i].b0, y1 = ed[i].b1;
      double ymin = (y0 < y1) ? y0 : y1, ymax = (y0 < y1) ? y1 : y0;
      int32_t a = (int32_t)ceil(ymin - 0.5), b = (int32_t)ceil(ymax - 0.5) - 1;
      if (a < j0) a = j0;
      if (b > j1) b = j1;
      elo[i] = a;
      ehi[i] = b;
    }
    const int32_t nr = j1 - j0 + 1;
    for (int32_t r = 0; r <= nr; r++)
      bcnt[r] = 0;
    for (int32_t i = 0; i < ne; i++)
      if (elo[i] <= ehi[i]) bcnt[elo[i] - j0 + 1]++;
    for (int32_t r = 0; r < nr; r++)
      bcnt[r + 1] += bcnt[r];
    for (int32_t i = 0; i < ne; i++)
      if (elo[i] <= ehi[i]) ebkt[bcnt[elo[i] - j0]++] = i;
    for (int32_t r = nr; r > 0; r--)
      bcnt[r] = bcnt[r - 1];
    bcnt[0] = 0;

    int32_t na = 0;
    for (int32_t j = j0; j <= j1; j++) {
      double yc = (double)j + 0.5;
      int32_t nc = 0, r = j - j0;
      st->nrow++;
      for (int32_t q = bcnt[r]; q < bcnt[r + 1]; q++)
        act[na++] = ebkt[q];
      for (int32_t q = 0; q < na;) {
        if (ehi[act[q]] < j)
          act[q] = act[--na];
        else
          q++;
      }
      for (int32_t q = 0; q < na; q++) {
        int32_t i = act[q];
        double y0 = ed[i].b0, y1 = ed[i].b1;
        /* Полуоткрытое условие: вершина ровно на строке считается один раз. */
        if ((y0 <= yc) == (y1 <= yc)) continue;
        st->nedge++;
        double tt = (yc - y0) / (y1 - y0);
        cr[nc].x = ed[i].a0 + tt * (ed[i].a1 - ed[i].a0);
        cr[nc].dir = (y1 > y0) ? 1 : -1;
        nc++;
      }
      if (nc < 2) continue;
      sort_cross(cr, nc);
      int wind = 0;
      for (int32_t c = 0; c + 1 < nc; c++) {
        wind += cr[c].dir;
        if (wind == 0) continue; /* правило НЕНУЛЕВОГО ОБОРОТА */
        int32_t i0 = (int32_t)ceil(cr[c].x - 0.5), i1 = (int32_t)ceil(cr[c + 1].x - 0.5) - 1;
        if (i0 < 0) i0 = 0;
        if (i1 > W - 1) i1 = W - 1;
        if (i1 < i0) continue;
        if (*nsp >= *cap) {
          int64_t nc2 = (*cap > 0) ? *cap * 2 : 4096;
          hz_span *ns = realloc(*spp, (size_t)nc2 * sizeof **spp);
          if (ns == NULL) {
            free(ed);
            free(cr);
            free(elo);
            free(ehi);
            free(ebkt);
            free(act);
            free(bcnt);
            return 2;
          }
          *spp = ns;
          *cap = nc2;
        }
        hz_span *s = &(*spp)[*nsp];
        s->j = j;
        s->i0 = i0;
        s->i1 = i1;
        s->poly = k;
        s->t0 = T0 + Tj * (double)j;
        s->ti = Ti;
        (*nsp)++;
        st->nspan++;
        st->nfrag += i1 - i0 + 1;
      }
    }
  }
  free(ed);
  free(cr);
  free(elo);
  free(ehi);
  free(ebkt);
  free(act);
  free(bcnt);
  return 0;
}

/* --- раскладка 1: сплошные пробеги ------------------------------------------ */

int hz_fbuf_init(hz_fbuf *fb, int64_t npix, int64_t cap) {
  memset(fb, 0, sizeof *fb);
  if (npix <= 0 || cap <= 0) return 1;
  fb->start = malloc((size_t)(npix + 1) * sizeof *fb->start);
  fb->poly = malloc((size_t)cap * sizeof *fb->poly);
  fb->depth = malloc((size_t)cap * sizeof *fb->depth);
  if (fb->start == NULL || fb->poly == NULL || fb->depth == NULL) {
    hz_fbuf_free(fb);
    return 2;
  }
  fb->cap = cap;
  return 0;
}

void hz_fbuf_free(hz_fbuf *fb) {
  free(fb->ord);
  free(fb->rcnt);
  free(fb->start);
  free(fb->poly);
  free(fb->depth);
  memset(fb, 0, sizeof *fb);
}

int hz_fbuf_build(hz_fbuf *fb, const hz_span *sp, int64_t nsp, int32_t W, int32_t H) {
  const int64_t npix = (int64_t)W * H;
  fb->W = W;
  fb->H = H;
  memset(fb->start, 0, (size_t)(npix + 1) * sizeof *fb->start);
  int64_t total = 0;
  for (int64_t s = 0; s < nsp; s++) {
    int64_t base = (int64_t)sp[s].j * W;
    for (int32_t i = sp[s].i0; i <= sp[s].i1; i++)
      fb->start[base + i + 1]++;
    total += sp[s].i1 - sp[s].i0 + 1;
  }
  if (total > fb->cap) return 3;
  for (int64_t p = 0; p < npix; p++)
    fb->start[p + 1] += fb->start[p];
  fb->n = total;
  /* ПОЛОСЫ ОБХОДЯТСЯ ПО ВОЗРАСТАНИЮ СТРОКИ, А НЕ В ПОРЯДКЕ ПОЛИГОНОВ (§129).
   * Замерено профилем: раскладка занимала `42.7 %` — `7` млн фрагментов на
   * направление ложились по `1.44` млн пикселей вразброс, и каждая запись
   * промахивалась мимо кэша. Полосы приходят сгруппированными по ПОЛИГОНУ, то
   * есть строки идут вперемешку; счётная сортировка по `j` стоит один проход по
   * полосам плюс `H` корзин (их тысяча, а не миллион) и делает запись почти
   * последовательной. Величина при этом не меняется вовсе: порядок фрагментов
   * внутри пикселя всё равно задаёт сортировка по глубине следом. */
  if (fb->ord == NULL || fb->ordcap < nsp) {
    int32_t *no = realloc(fb->ord, (size_t)(nsp > 0 ? nsp : 1) * sizeof *no);
    if (no == NULL) return 2;
    fb->ord = no;
    fb->ordcap = nsp;
  }
  if (fb->rcnt == NULL || fb->rcap < H + 1) {
    int32_t *nr = realloc(fb->rcnt, (size_t)(H + 1) * sizeof *nr);
    if (nr == NULL) return 2;
    fb->rcnt = nr;
    fb->rcap = H + 1;
  }
  memset(fb->rcnt, 0, (size_t)(H + 1) * sizeof *fb->rcnt);
  for (int64_t s = 0; s < nsp; s++)
    fb->rcnt[sp[s].j + 1]++;
  for (int32_t j = 0; j < H; j++)
    fb->rcnt[j + 1] += fb->rcnt[j];
  for (int64_t s = 0; s < nsp; s++)
    fb->ord[fb->rcnt[sp[s].j]++] = (int32_t)s;
  /* Курсор поверх `start`: восстанавливается сдвигом, отдельный массив на
   * миллионы пикселей не нужен. Заполняем в start[p], потом сдвигаем обратно. */
  for (int64_t q = 0; q < nsp; q++) {
    const hz_span *S = &sp[fb->ord[q]];
    int64_t base = (int64_t)S->j * W;
    for (int32_t i = S->i0; i <= S->i1; i++) {
      int32_t f = fb->start[base + i]++;
      fb->poly[f] = S->poly;
      fb->depth[f] = S->t0 + S->ti * (double)i;
    }
  }
  for (int64_t p = npix; p > 0; p--)
    fb->start[p] = fb->start[p - 1];
  fb->start[0] = 0;
  return 0;
}

void hz_fbuf_sort(hz_fbuf *fb) {
  const int64_t npix = (int64_t)fb->W * fb->H;
  for (int64_t p = 0; p < npix; p++) {
    int32_t b = fb->start[p], e = fb->start[p + 1];
    for (int32_t x = b + 1; x < e; x++) {
      double dv = fb->depth[x];
      int32_t pv = fb->poly[x], y = x - 1;
      while (y >= b && fb->depth[y] > dv) {
        fb->depth[y + 1] = fb->depth[y];
        fb->poly[y + 1] = fb->poly[y];
        y--;
      }
      fb->depth[y + 1] = dv;
      fb->poly[y + 1] = pv;
    }
  }
}

/* --- раскладка 2: односвязные списки ---------------------------------------- */

int hz_abuf_init(hz_abuf *ab, int64_t npix, int64_t cap) {
  memset(ab, 0, sizeof *ab);
  if (npix <= 0 || cap <= 0) return 1;
  ab->head = malloc((size_t)npix * sizeof *ab->head);
  ab->next = malloc((size_t)cap * sizeof *ab->next);
  ab->poly = malloc((size_t)cap * sizeof *ab->poly);
  ab->depth = malloc((size_t)cap * sizeof *ab->depth);
  if (ab->head == NULL || ab->next == NULL || ab->poly == NULL || ab->depth == NULL) {
    hz_abuf_free(ab);
    return 2;
  }
  ab->cap = cap;
  return 0;
}

void hz_abuf_free(hz_abuf *ab) {
  free(ab->head);
  free(ab->next);
  free(ab->poly);
  free(ab->depth);
  memset(ab, 0, sizeof *ab);
}

int hz_abuf_build(hz_abuf *ab, const hz_span *sp, int64_t nsp, int32_t W, int32_t H) {
  const int64_t npix = (int64_t)W * H;
  ab->W = W;
  ab->H = H;
  ab->n = 0;
  for (int64_t p = 0; p < npix; p++)
    ab->head[p] = -1;
  for (int64_t s = 0; s < nsp; s++) {
    int64_t base = (int64_t)sp[s].j * W;
    for (int32_t i = sp[s].i0; i <= sp[s].i1; i++) {
      if (ab->n >= ab->cap) return 3;
      int64_t f = ab->n++;
      ab->poly[f] = sp[s].poly;
      ab->depth[f] = sp[s].t0 + sp[s].ti * (double)i;
      ab->next[f] = ab->head[base + i];
      ab->head[base + i] = (int32_t)f;
    }
  }
  return 0;
}

void hz_abuf_sort(hz_abuf *ab) {
  const int64_t npix = (int64_t)ab->W * ab->H;
  for (int64_t px = 0; px < npix; px++) {
    int32_t h = ab->head[px];
    if (h < 0 || ab->next[h] < 0) continue; /* пусто или один фрагмент */
    int32_t sorted = -1;
    while (h >= 0) {
      int32_t nx = ab->next[h];
      if (sorted < 0 || ab->depth[h] <= ab->depth[sorted]) {
        ab->next[h] = sorted;
        sorted = h;
      } else {
        int32_t c = sorted;
        while (ab->next[c] >= 0 && ab->depth[ab->next[c]] < ab->depth[h])
          c = ab->next[c];
        ab->next[h] = ab->next[c];
        ab->next[c] = h;
      }
      h = nx;
    }
    ab->head[px] = sorted;
  }
}
