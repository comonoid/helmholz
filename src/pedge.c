/* Согласованное упрощение края. Разбор и оговорки — в `pedge.h`. */

#include "pedge.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* --- расстояние точки до отрезка В МИРЕ --------------------------------------
 * Именно в мире, а не в (u,v): у двух владельцев цепи рамы разные, и одна и та
 * же дуга в плоскости A может спроецироваться в прямую в плоскости B. */
static double seg_dist3(const double *p, const double *a, const double *b) {
  double d[3], q[3], L2 = 0.0, t = 0.0;
  for (int c = 0; c < 3; c++) {
    d[c] = b[c] - a[c];
    L2 += d[c] * d[c];
  }
  if (L2 > 0.0) {
    for (int c = 0; c < 3; c++)
      t += (p[c] - a[c]) * d[c];
    t /= L2;
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;
  }
  double s = 0.0;
  for (int c = 0; c < 3; c++) {
    q[c] = a[c] + t * d[c] - p[c];
    s += q[c] * q[c];
  }
  return sqrt(s);
}

/* --- Дуглас–Пекер, ЯВНЫЙ СТЕК ------------------------------------------------
 * Не рекурсия: глубина у неё равна числу оставленных вершин, а петли города
 * доходят до десятков тысяч. Стек растёт, а не ограничен потолком — молчаливое
 * усечение уже стоило нам отказа дробления в Ш6. */
typedef struct {
  int32_t i0, i1;
} pe_span;

typedef struct {
  pe_span *s;
  int64_t n, cap;
} pe_stack;

static int st_push(pe_stack *k, int32_t i0, int32_t i1) {
  if (k->n >= k->cap) {
    int64_t nc = (k->cap > 0) ? k->cap * 2 : 64;
    pe_span *q = realloc(k->s, (size_t)nc * sizeof *q);
    if (q == NULL) return 2;
    k->s = q;
    k->cap = nc;
  }
  k->s[k->n].i0 = i0;
  k->s[k->n].i1 = i1;
  k->n++;
  return 0;
}

/* Открытая цепь: `pid[i]` — номер точки в `pts` (по 3 double), `keep` — по
 * ЛОКАЛЬНОМУ номеру `i`. Концы помечаются всегда: они якоря, и решать за них
 * упрощение не вправе. */
static int dp_open(const double *pts, const int32_t *pid, int32_t n, double tol,
                   unsigned char *keep, pe_stack *k) {
  if (n < 2) {
    if (n == 1) keep[0] = 1;
    return 0;
  }
  keep[0] = 1;
  keep[n - 1] = 1;
  k->n = 0;
  if (st_push(k, 0, n - 1) != 0) return 2;
  while (k->n > 0) {
    k->n--;
    int32_t i0 = k->s[k->n].i0, i1 = k->s[k->n].i1;
    if (i1 <= i0 + 1) continue;
    const double *A = pts + (size_t)pid[i0] * 3, *B = pts + (size_t)pid[i1] * 3;
    double dm = -1.0;
    int32_t im = i0;
    for (int32_t i = i0 + 1; i < i1; i++) {
      double d = seg_dist3(pts + (size_t)pid[i] * 3, A, B);
      if (d > dm) {
        dm = d;
        im = i;
      }
    }
    if (dm <= tol) continue;
    keep[im] = 1;
    if (st_push(k, i0, im) != 0) return 2;
    if (st_push(k, im, i1) != 0) return 2;
  }
  return 0;
}

/* Замкнутая цепь: `pid[0..n−1]` плюс `pid[n] == pid[0]` (замыкание уже
 * дописано вызывающим). Якорей два — первая вершина и САМАЯ ДАЛЬНЯЯ от неё;
 * второй якорь берётся геометрически, а не «серединой», чтобы результат не
 * зависел от того, где петля начинается. */
static int dp_cycle(const double *pts, const int32_t *pid, int32_t n, double tol,
                    unsigned char *keep, pe_stack *k) {
  if (n < 3) {
    for (int32_t i = 0; i <= n; i++)
      keep[i] = 1;
    return 0;
  }
  const double *A = pts + (size_t)pid[0] * 3;
  int32_t far = 1;
  double dm = -1.0;
  for (int32_t i = 1; i < n; i++) {
    const double *B = pts + (size_t)pid[i] * 3;
    double d = 0.0;
    for (int c = 0; c < 3; c++)
      d += (B[c] - A[c]) * (B[c] - A[c]);
    if (d > dm) {
      dm = d;
      far = i;
    }
  }
  if (dp_open(pts, pid, far + 1, tol, keep, k) != 0) return 2;
  if (dp_open(pts, pid + far, n - far + 1, tol, keep + far, k) != 0) return 2;
  return 0;
}

/* ПОЛ В ТРИ ВЕРШИНЫ НА ПЕТЛЮ. У замкнутой цепи якорей два, и на грубом допуске
 * петля законно сжимается до них — то есть до ОТРЕЗКА, у которого площадь ноль.
 * Восстанавливать при этом петлю ЦЕЛИКОМ (первое, что приходит в голову) —
 * ошибка, и она ловится подписью из §10: число вершин перестаёт зависеть от
 * допуска и даже РАСТЁТ с ним, потому что чем грубее допуск, тем больше петель
 * восстановлено целиком. Правильный пол — ТРИ вершины: два якоря плюс самая
 * дальняя от их хорды, то есть треугольник наибольшей площади при этих якорях.
 * Выбор канонический и от полигона не зависит, поэтому в согласованном режиме
 * добавка идёт в общий `keepw` и её видят ОБА владельца. */
static void loop_floor3(const hz_polyset *ps, const double *pos, int32_t b, int32_t n,
                        int32_t out[3]) {
  const double *A = pos + (size_t)ps->bw[b] * 3;
  int32_t far = 1;
  double dm = -1.0;
  for (int32_t i = 1; i < n; i++) {
    const double *B = pos + (size_t)ps->bw[b + i] * 3;
    double d = 0.0;
    for (int c = 0; c < 3; c++)
      d += (B[c] - A[c]) * (B[c] - A[c]);
    if (d > dm) {
      dm = d;
      far = i;
    }
  }
  int32_t th = -1;
  double dm2 = -1.0;
  const double *F = pos + (size_t)ps->bw[b + far] * 3;
  for (int32_t i = 1; i < n; i++) {
    if (i == far) continue;
    double d = seg_dist3(pos + (size_t)ps->bw[b + i] * 3, A, F);
    if (d > dm2) {
      dm2 = d;
      th = i;
    }
  }
  out[0] = 0;
  out[1] = far;
  out[2] = (th >= 0) ? th : (far + 1) % n;
}

/* --- несопряжённые рёбра ------------------------------------------------------ */

static int cmp_u64(const void *a, const void *b) {
  uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

void hz_edge_unpaired(const hz_polyset *ps, int64_t *nedge, int64_t *nsingle) {
  *nedge = 0;
  *nsingle = 0;
  if (ps->nbv <= 0 || ps->bw == NULL || ps->nweld <= 0) return;
  uint64_t nw = (uint64_t)ps->nweld;
  uint64_t *key = malloc((size_t)ps->nbv * sizeof *key);
  if (key == NULL) return;
  int64_t ne = 0;
  for (int32_t k = 0; k < ps->np; k++) {
    const hz_poly *P = &ps->p[k];
    for (int32_t l = P->l0; l < P->l0 + P->nloop; l++) {
      int32_t b = ps->loop[l], e = ps->loop[l + 1], n = e - b;
      if (n < 2) continue;
      for (int32_t i = 0; i < n; i++) {
        int32_t a = ps->bw[b + i], c = ps->bw[b + (i + 1) % n];
        if (a < 0 || c < 0 || a == c) continue;
        uint64_t lo = (uint64_t)((a < c) ? a : c), hi = (uint64_t)((a < c) ? c : a);
        key[ne++] = lo * nw + hi;
      }
    }
  }
  qsort(key, (size_t)ne, sizeof *key, cmp_u64);
  int64_t single = 0, tot = 0;
  for (int64_t i = 0; i < ne;) {
    int64_t j = i;
    while (j < ne && key[j] == key[i])
      j++;
    tot++;
    if (j - i == 1) single++;
    i = j;
  }
  free(key);
  *nedge = tot;
  *nsingle = single;
}

/* --- само упрощение ----------------------------------------------------------- */

typedef struct {
  int32_t *nb0, *nb1, *chain, *occ;
  unsigned char *lock, *have, *vis, *keepw, *keepbv, *kloc;
  double *pos;
  pe_stack stk;
} pe_work;

static void pe_free(pe_work *w) {
  free(w->nb0);
  free(w->nb1);
  free(w->chain);
  free(w->occ);
  free(w->lock);
  free(w->have);
  free(w->vis);
  free(w->keepw);
  free(w->keepbv);
  free(w->kloc);
  free(w->pos);
  free(w->stk.s);
  memset(w, 0, sizeof *w);
}

int hz_edge_simplify(hz_polyset *ps, double tol, hz_edgemode mode, hz_edgestat *st) {
  memset(st, 0, sizeof *st);
  st->nbv_in = ps->nbv;
  st->nbv_out = ps->nbv;
  if (ps->nbv <= 0 || ps->bw == NULL) return 0;
  if (!(tol > 0.0)) return 0;
  int64_t nw64 = ps->nweld;
  if (nw64 <= 0 || nw64 > INT32_MAX) return 0;
  const int32_t nw = (int32_t)nw64;

  /* Рабочая длина цепи. В согласованном режиме цепь не длиннее числа сварных
   * вершин, в контрольном — не длиннее самой длинной ПЕТЛИ, а петля может
   * проходить одну вершину дважды. Берётся большее из двух: молчаливое усечение
   * буфера — ровно тот класс отказа, что стоил нам дробления в Ш6. */
  int32_t maxloop = 0;
  for (int32_t l = 0; l < ps->nloopall; l++) {
    int32_t n = ps->loop[l + 1] - ps->loop[l];
    if (n > maxloop) maxloop = n;
  }
  const size_t wcap = (size_t)((nw > maxloop) ? nw : maxloop) + 2;

  pe_work w;
  memset(&w, 0, sizeof w);
  w.nb0 = malloc((size_t)nw * sizeof *w.nb0);
  w.nb1 = malloc((size_t)nw * sizeof *w.nb1);
  w.chain = malloc(wcap * sizeof *w.chain);
  w.occ = calloc((size_t)nw, sizeof *w.occ);
  w.lock = calloc((size_t)nw, 1);
  w.have = calloc((size_t)nw, 1);
  w.vis = calloc((size_t)nw, 1);
  w.keepw = calloc((size_t)nw, 1);
  w.keepbv = calloc((size_t)ps->nbv, 1);
  w.kloc = calloc(wcap, 1);
  w.pos = malloc((size_t)nw * 3 * sizeof *w.pos);
  if (w.nb0 == NULL || w.nb1 == NULL || w.chain == NULL || w.occ == NULL || w.lock == NULL ||
      w.have == NULL || w.vis == NULL || w.keepw == NULL || w.keepbv == NULL || w.kloc == NULL ||
      w.pos == NULL) {
    pe_free(&w);
    return 2;
  }
  for (int32_t i = 0; i < nw; i++) {
    w.nb0[i] = -1;
    w.nb1[i] = -1;
  }

  /* --- 1. соседи по краю и КАНОНИЧЕСКАЯ мировая точка ---
   * Точка берётся у ПЕРВОГО владельца в порядке обхода. Владельцы дают разные
   * точки (каждый сносит вершину на свою плоскость), расходясь не больше чем на
   * `dmax ≤ δ`; выбрать надо ОДНУ, иначе цепь у двух соседей была бы разной
   * ломаной и мера отклонения — разной величиной. */
  for (int32_t k = 0; k < ps->np; k++) {
    const hz_poly *P = &ps->p[k];
    for (int32_t l = P->l0; l < P->l0 + P->nloop; l++) {
      int32_t b = ps->loop[l], e = ps->loop[l + 1], n = e - b;
      if (n < 3) continue;
      for (int32_t i = 0; i < n; i++) {
        int32_t v = ps->bw[b + i];
        if (v < 0) continue;
        w.occ[v]++;
        if (!w.have[v]) {
          w.have[v] = 1;
          hz_poly_world(P, ps->bv[(size_t)(b + i) * 2], ps->bv[(size_t)(b + i) * 2 + 1],
                        w.pos + (size_t)v * 3);
        }
        int32_t pv = ps->bw[b + (i + n - 1) % n], nx = ps->bw[b + (i + 1) % n];
        /* СТЫК ЗАПИРАЕТСЯ. Несварной сосед, совпадение соседей (защемление) и
         * петля через саму вершину — всё это случаи, где «выбросить вершину»
         * означало бы решать сразу за несколько цепей. */
        if (pv < 0 || nx < 0 || pv == nx || pv == v || nx == v) {
          w.lock[v] = 1;
          continue;
        }
        int32_t two[2] = {pv, nx};
        for (int s = 0; s < 2; s++) {
          int32_t x = two[s];
          if (w.nb0[v] < 0 || w.nb0[v] == x)
            w.nb0[v] = x;
          else if (w.nb1[v] < 0 || w.nb1[v] == x)
            w.nb1[v] = x;
          else
            w.lock[v] = 1; /* соседей больше двух — сходятся три полигона и более */
        }
      }
    }
  }
  for (int32_t v = 0; v < nw; v++) {
    if (!w.have[v]) continue;
    if (w.nb0[v] < 0 || w.nb1[v] < 0) w.lock[v] = 1;
    /* ОТКРЫТЫЙ КРАЙ: вершина входит в край ровно один раз, соседа за ней нет.
     * Согласовывать тут не с кем, и сдвиг такой вершины — не «шов уехал», а
     * дыра в силуэте. */
    if (mode == HZ_EDGE_SHARED_INNER && w.occ[v] < 2) w.lock[v] = 1;
    if (w.lock[v])
      st->nlock++;
    else
      st->nfree++;
  }

  if (mode != HZ_EDGE_INDEP) {
    /* --- 2. ЦЕПИ: упрощаются ОДИН раз, результат читают оба владельца --- */
    for (int32_t v0 = 0; v0 < nw; v0++) {
      if (!w.have[v0] || w.lock[v0] || w.vis[v0]) continue;
      /* назад до запертой вершины либо до замыкания */
      int32_t a = v0, ap = w.nb1[v0];
      int cyc = 0;
      int32_t term = -1;
      for (int32_t guard = 0; guard <= nw; guard++) {
        int32_t nx = (w.nb0[a] == ap) ? w.nb1[a] : w.nb0[a];
        if (nx == v0) {
          cyc = 1;
          break;
        }
        if (nx < 0 || w.lock[nx]) {
          term = nx;
          break;
        }
        ap = a;
        a = nx;
      }
      int32_t nc = 0;
      if (cyc) {
        int32_t prev = w.nb1[v0], cur = w.nb0[v0];
        w.chain[nc++] = v0;
        w.vis[v0] = 1;
        while (cur != v0 && nc <= nw) {
          w.chain[nc++] = cur;
          w.vis[cur] = 1;
          int32_t nx = (w.nb0[cur] == prev) ? w.nb1[cur] : w.nb0[cur];
          prev = cur;
          cur = nx;
        }
        w.chain[nc] = w.chain[0]; /* замыкание — для dp_cycle */
        memset(w.kloc, 0, (size_t)nc + 1);
        if (dp_cycle(w.pos, w.chain, nc, tol, w.kloc, &w.stk) != 0) {
          pe_free(&w);
          return 2;
        }
        for (int32_t i = 0; i <= nc; i++)
          if (w.kloc[i]) w.keepw[w.chain[i]] = 1;
        st->ncycle++;
      } else {
        if (term < 0) { /* оборванный край: запереть цепь целиком */
          w.vis[v0] = 1;
          w.keepw[v0] = 1;
          continue;
        }
        int32_t prev = term, cur = a;
        w.chain[nc++] = term;
        while (!w.lock[cur] && nc <= nw) {
          w.chain[nc++] = cur;
          w.vis[cur] = 1;
          int32_t nx = (w.nb0[cur] == prev) ? w.nb1[cur] : w.nb0[cur];
          prev = cur;
          cur = nx;
          if (cur < 0) break;
        }
        if (cur >= 0) w.chain[nc++] = cur;
        memset(w.kloc, 0, (size_t)nc);
        if (dp_open(w.pos, w.chain, nc, tol, w.kloc, &w.stk) != 0) {
          pe_free(&w);
          return 2;
        }
        for (int32_t i = 0; i < nc; i++)
          if (w.kloc[i]) w.keepw[w.chain[i]] = 1;
      }
      st->nchain++;
    }
    for (int32_t v = 0; v < nw; v++)
      if (w.have[v] && w.lock[v]) w.keepw[v] = 1;
    /* пол в три вершины — ДО раскрытия в вершины края, чтобы добавку увидели
     * ОБА владельца петли, а не тот полигон, у которого она недосчиталась */
    for (int32_t k = 0; k < ps->np; k++) {
      const hz_poly *P = &ps->p[k];
      for (int32_t l = P->l0; l < P->l0 + P->nloop; l++) {
        int32_t b = ps->loop[l], e = ps->loop[l + 1], n = e - b, cnt = 0, bad = (n < 3);
        for (int32_t i = 0; i < n && !bad; i++)
          if (ps->bw[b + i] < 0) bad = 1;
        if (bad) continue;
        for (int32_t i = 0; i < n; i++)
          if (w.keepw[ps->bw[b + i]]) cnt++;
        if (cnt >= 3) continue;
        int32_t add[3];
        loop_floor3(ps, w.pos, b, n, add);
        for (int j = 0; j < 3; j++)
          w.keepw[ps->bw[b + add[j]]] = 1;
        st->nloop_short++;
      }
    }
    for (int32_t i = 0; i < ps->nbv; i++) {
      int32_t v = ps->bw[i];
      w.keepbv[i] = (v < 0) ? 1 : w.keepw[v];
    }
  } else {
    /* --- 2'. НЕГАТИВНЫЙ КОНТРОЛЬ: каждый полигон решает за себя ---
     * Отличие от согласованного режима РОВНО одно — единица упрощения. Метод,
     * допуск, мера расстояния и канонические мировые точки те же, поэтому
     * разница в замере есть разница СХЕМ, а не реализаций. */
    for (int32_t k = 0; k < ps->np; k++) {
      const hz_poly *P = &ps->p[k];
      for (int32_t l = P->l0; l < P->l0 + P->nloop; l++) {
        int32_t b = ps->loop[l], e = ps->loop[l + 1], n = e - b;
        int bad = (n < 3);
        for (int32_t i = 0; i < n && !bad; i++)
          if (ps->bw[b + i] < 0) bad = 1;
        if (bad) {
          for (int32_t i = 0; i < n; i++)
            w.keepbv[b + i] = 1;
          continue;
        }
        for (int32_t i = 0; i < n; i++)
          w.chain[i] = ps->bw[b + i];
        w.chain[n] = w.chain[0];
        memset(w.kloc, 0, (size_t)n + 1);
        if (dp_cycle(w.pos, w.chain, n, tol, w.kloc, &w.stk) != 0) {
          pe_free(&w);
          return 2;
        }
        for (int32_t i = 0; i < n; i++)
          if (w.kloc[i]) w.keepbv[b + i] = 1;
        if (w.kloc[n]) w.keepbv[b] = 1;
        int32_t cnt = 0;
        for (int32_t i = 0; i < n; i++)
          if (w.keepbv[b + i]) cnt++;
        if (cnt < 3) {
          int32_t add[3];
          loop_floor3(ps, w.pos, b, n, add);
          for (int j = 0; j < 3; j++)
            w.keepbv[b + add[j]] = 1;
          st->nloop_short++;
        }
        st->nchain++;
      }
    }
  }

  /* --- 3. пересборка петель --- */
  /* calloc, а не malloc: пересборка заполняет только оставленные вершины, и
   * анализатор об этом не знает; нулевое дно заодно делает любую будущую
   * ошибку индексации воспроизводимой, а не мусорной. */
  double *bv2 = calloc((size_t)ps->nbv * 2, sizeof *bv2);
  int32_t *bw2 = calloc((size_t)ps->nbv, sizeof *bw2);
  int32_t *lo2 = malloc(((size_t)ps->nloopall + 2) * sizeof *lo2);
  if (bv2 == NULL || bw2 == NULL || lo2 == NULL) {
    free(bv2);
    free(bw2);
    free(lo2);
    pe_free(&w);
    return 2;
  }
  int32_t out = 0, lout = 0;
  lo2[0] = 0;
  for (int32_t k = 0; k < ps->np; k++) {
    hz_poly *P = &ps->p[k];
    int32_t l0new = lout;
    for (int32_t l = P->l0; l < P->l0 + P->nloop; l++) {
      int32_t b = ps->loop[l], e = ps->loop[l + 1], n = e - b, cnt = 0;
      for (int32_t i = 0; i < n; i++)
        if (w.keepbv[b + i]) cnt++;
      /* Петля ниже трёх вершин сюда уже не доходит — пол поставлен выше
       * (`loop_floor3`). Проверка оставлена как СТРАЖ: восстановить петлю
       * целиком безопасно, а молча выдать отрезок вместо петли — нет. */
      int whole = (cnt < 3);
      if (whole) st->nloop_short++;
      for (int32_t i = 0; i < n; i++) {
        if (!whole && !w.keepbv[b + i]) continue;
        bv2[(size_t)out * 2 + 0] = ps->bv[(size_t)(b + i) * 2 + 0];
        bv2[(size_t)out * 2 + 1] = ps->bv[(size_t)(b + i) * 2 + 1];
        bw2[out] = ps->bw[b + i];
        out++;
      }
      lout++;
      lo2[lout] = out;
    }
    P->l0 = l0new;
  }
  free(ps->bv);
  free(ps->bw);
  free(ps->loop);
  ps->bv = bv2;
  ps->bw = bw2;
  ps->loop = lo2;
  ps->nbv = out;
  ps->nloopall = lout;
  st->nbv_out = out;

  /* габарит края пересчитывается: по нему идёт дешёвый отсев луча */
  for (int32_t k = 0; k < ps->np; k++) {
    hz_poly *P = &ps->p[k];
    P->uvlo[0] = P->uvlo[1] = 1e300;
    P->uvhi[0] = P->uvhi[1] = -1e300;
    int32_t b = (P->nloop > 0) ? ps->loop[P->l0] : 0;
    int32_t e = (P->nloop > 0) ? ps->loop[P->l0 + P->nloop] : 0;
    for (int32_t i = b; i < e; i++)
      for (int a = 0; a < 2; a++) {
        double x = ps->bv[(size_t)i * 2 + (size_t)a];
        if (x < P->uvlo[a]) P->uvlo[a] = x;
        if (x > P->uvhi[a]) P->uvhi[a] = x;
      }
    if (!(P->uvlo[0] <= P->uvhi[0])) {
      P->uvlo[0] = P->uvlo[1] = 0.0;
      P->uvhi[0] = P->uvhi[1] = 0.0;
    }
  }
  pe_free(&w);
  return 0;
}
