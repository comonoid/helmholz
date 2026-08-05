/* ptree.c — пространственное дерево над треугольниками. Разбор — в `ptree.h`. */
#include "ptree.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int pt_grow_nd(hz_ptree *t, int32_t need) {
  if (t->nnd + need <= t->ndcap) return 0;
  int32_t nc = (t->ndcap > 0) ? t->ndcap * 2 : 1024;
  while (nc < t->nnd + need)
    nc *= 2;
  hz_ptnode *n = realloc(t->nd, (size_t)nc * sizeof *n);
  if (n == NULL) return 1;
  t->nd = n;
  t->ndcap = nc;
  return 0;
}

static int pt_grow_ref(hz_ptree *t, int64_t need) {
  if (t->nref + need <= t->refcap) return 0;
  int64_t nc = (t->refcap > 0) ? t->refcap * 2 : 65536;
  while (nc < t->nref + need)
    nc *= 2;
  int32_t *r = realloc(t->ref, (size_t)nc * sizeof *r);
  if (r == NULL) return 1;
  t->ref = r;
  t->refcap = nc;
  return 0;
}

/* Габарит треугольника. */
static void pt_tri_bb(const hz_objmesh *m, int32_t tr, double lo[3], double hi[3]) {
  for (int c = 0; c < 3; c++) {
    lo[c] = 1e300;
    hi[c] = -1e300;
  }
  for (int i = 0; i < 3; i++) {
    const double *p = m->v + 3 * (size_t)m->f[(size_t)tr * 3 + (size_t)i];
    for (int c = 0; c < 3; c++) {
      if (p[c] < lo[c]) lo[c] = p[c];
      if (p[c] > hi[c]) hi[c] = p[c];
    }
  }
}

/* Рекурсивное дробление. `list` — рабочий список треугольников узла. */
static int pt_split(hz_ptree *t, const hz_objmesh *m, int32_t nid, int32_t *list, int32_t n,
                    int lev, const double *tlo, const double *thi) {
  hz_ptnode *N = &t->nd[nid];
  if (n <= t->leafmax || lev >= t->maxlev) {
    if (pt_grow_ref(t, n) != 0) return 1;
    N->child = -1;
    N->t0 = (int32_t)t->nref;
    N->ntri = n;
    memcpy(t->ref + t->nref, list, (size_t)n * sizeof *list);
    t->nref += n;
    t->nleaf++;
    t->n_leaf_tri += n;
    return 0;
  }
  double mid[3];
  for (int c = 0; c < 3; c++)
    mid[c] = 0.5 * (N->lo[c] + N->hi[c]);
  if (pt_grow_nd(t, 8) != 0) return 1;
  N = &t->nd[nid];
  int32_t c0 = t->nnd;
  N->child = c0;
  N->t0 = 0;
  N->ntri = 0;
  t->nnd += 8;
  for (int k = 0; k < 8; k++) {
    hz_ptnode *C = &t->nd[c0 + k];
    for (int c = 0; c < 3; c++) {
      C->lo[c] = (k & (1 << c)) ? mid[c] : t->nd[nid].lo[c];
      C->hi[c] = (k & (1 << c)) ? t->nd[nid].hi[c] : mid[c];
    }
    C->child = -1;
    C->t0 = 0;
    C->ntri = 0;
  }
  /* ТРЕУГОЛЬНИК ЖИВЁТ В САМОМ ГЛУБОКОМ УЗЛЕ, КОТОРЫЙ СОДЕРЖИТ ЕГО ЦЕЛИКОМ.
   * Первая редакция ссылалась на него из КАЖДОГО пересекаемого листа, и замер
   * §95.1 показал катастрофу: избыточность 1193.36 ссылки на треугольник,
   * сборка 361 с вместо 8.3. Причина та же, что убила ограничение размера
   * элемента (А235) — огромные входные треугольники: пол в 18.85 м² попадал
   * во ВСЕ листья. При этом правиле избыточность равна ЕДИНИЦЕ ТОЧНО, а обход
   * обязан рисовать треугольники КАЖДОГО посещённого узла, а не только листа. */
  int32_t *sub = malloc((size_t)n * sizeof *sub);
  int32_t *keep = malloc((size_t)n * sizeof *keep);
  if (sub == NULL || keep == NULL) {
    free(sub);
    free(keep);
    return 1;
  }
  /* РЫХЛАЯ УКЛАДКА (§196, ключ `t->loose`). Строгое правило «влезает целиком»
   * оставляет у внутренних узлов 71…79 % треугольников (замерено по трём
   * сценам), и `leafmax` их размер не ограничивает. Рыхлое правило отправляет
   * треугольник вниз ВСЕГДА — по октанту его ЦЕНТРА, — а коробка ребёнка
   * расширяется, чтобы его содержать. Цена названа в §196: коробки детей
   * начинают ПЕРЕКРЫВАТЬСЯ, и обход спереди назад по порядку октантов
   * перестаёт быть точным. Поэтому правило только по ключу и НИКОГДА не
   * умолчание: у `pvis` точный порядок обхода — несущее свойство. */
  int32_t *owner = malloc((size_t)n * sizeof *owner);
  if (owner == NULL) {
    free(sub);
    free(keep);
    return 1;
  }
  int32_t nk = 0;
  for (int32_t i = 0; i < n; i++) {
    const double *bl = tlo + 3 * (size_t)list[i], *bh = thi + 3 * (size_t)list[i];
    int own = -1;
    if (t->loose) {
      own = 0;
      for (int c = 0; c < 3; c++)
        if (0.5 * (bl[c] + bh[c]) >= mid[c]) own |= (1 << c);
      hz_ptnode *C = &t->nd[c0 + own];
      for (int c = 0; c < 3; c++) {
        if (bl[c] < C->lo[c]) C->lo[c] = bl[c];
        if (bh[c] > C->hi[c]) C->hi[c] = bh[c];
      }
    } else {
      for (int k = 0; k < 8 && own < 0; k++) {
        const double *clo = t->nd[c0 + k].lo, *chi = t->nd[c0 + k].hi;
        int in = 1;
        for (int c = 0; c < 3 && in; c++)
          if (bl[c] < clo[c] || bh[c] > chi[c]) in = 0;
        if (in) own = k;
      }
    }
    owner[i] = own;
    if (own < 0) keep[nk++] = list[i];
  }
  if (pt_grow_ref(t, nk) != 0) {
    free(sub);
    free(keep);
    return 1;
  }
  t->nd[nid].t0 = (int32_t)t->nref;
  t->nd[nid].ntri = nk;
  t->n_inner += nk;
  memcpy(t->ref + t->nref, keep, (size_t)nk * sizeof *keep);
  t->nref += nk;
  free(keep);
  /* Раскладка по детям идёт ПО ТОМУ ЖЕ `owner`, что и подсчёт остатка выше.
   * Прежняя редакция перепроверяла условие вторым циклом — при рыхлом правиле
   * это дало бы РАСХОЖДЕНИЕ двух проверок (коробки уже расширены), то есть
   * ровно тот класс, где два места считают одно и то же по-разному. */
  for (int k = 0; k < 8; k++) {
    int32_t ns = 0;
    for (int32_t i = 0; i < n; i++)
      if (owner[i] == k) sub[ns++] = list[i];
    if (pt_split(t, m, c0 + k, sub, ns, lev + 1, tlo, thi) != 0) {
      free(sub);
      free(owner);
      return 1;
    }
  }
  free(sub);
  free(owner);
  return 0;
}

int hz_ptree_build(hz_ptree *t, const hz_objmesh *m, int leafmax, int maxlev) {
  return hz_ptree_build_ex(t, m, leafmax, maxlev, 0);
}

int hz_ptree_build_ex(hz_ptree *t, const hz_objmesh *m, int leafmax, int maxlev, int loose) {
  memset(t, 0, sizeof *t);
  t->loose = loose;
  t->leafmax = (leafmax > 0) ? leafmax : 16;
  t->maxlev = (maxlev > 0) ? maxlev : 12;
  if (m->nt <= 0) return 1;
  double *tlo = malloc(3 * (size_t)m->nt * sizeof *tlo);
  double *thi = malloc(3 * (size_t)m->nt * sizeof *thi);
  int32_t *list = malloc((size_t)m->nt * sizeof *list);
  if (tlo == NULL || thi == NULL || list == NULL) {
    free(tlo);
    free(thi);
    free(list);
    return 2;
  }
  for (int32_t k = 0; k < m->nt; k++) {
    pt_tri_bb(m, k, tlo + 3 * (size_t)k, thi + 3 * (size_t)k);
    list[k] = k;
  }
  if (pt_grow_nd(t, 1) != 0) {
    free(tlo);
    free(thi);
    free(list);
    return 2;
  }
  t->nnd = 1;
  for (int c = 0; c < 3; c++) {
    t->nd[0].lo[c] = m->lo[c] - 1e-6;
    t->nd[0].hi[c] = m->hi[c] + 1e-6;
  }
  t->nd[0].child = -1;
  t->nd[0].t0 = 0;
  t->nd[0].ntri = 0;
  int rc = pt_split(t, m, 0, list, m->nt, 0, tlo, thi);
  free(tlo);
  free(thi);
  free(list);
  if (rc != 0) return 2;
  t->redundancy = (double)t->nref / (double)m->nt;
  return 0;
}

/* ---- УКЛАДКА РЕЗКОЙ (§241.2, шаг О70) ------------------------------------ */

/* Спуск и градуировка меряются ПОРОЗНЬ (А476): иначе плата за структуру и
 * плата за градуировку окажутся в одном числе. */
static double pt_now(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* Рост узлов ВМЕСТЕ с массивом уровней: два массива с одним счётчиком должны
 * расти одним вызовом, иначе рано или поздно они разойдутся на одном возврате. */
static int pt_grow_cut(hz_ptree *t, int32_t need) {
  if (t->nnd + need <= t->ndcap) return 0;
  int32_t nc = (t->ndcap > 0) ? t->ndcap * 2 : 1024;
  while (nc < t->nnd + need)
    nc *= 2;
  hz_ptnode *n = realloc(t->nd, (size_t)nc * sizeof *n);
  if (n == NULL) return 1;
  t->nd = n;
  unsigned char *L = realloc(t->lev, (size_t)nc * sizeof *L);
  if (L == NULL) return 1;
  t->lev = L;
  t->ndcap = nc;
  return 0;
}

/* Завести восемь детей узла `nid`. Координата середины считается ОДИН РАЗ и
 * копируется в двух детей — на этом стоит побитовое совпадение общей грани у
 * соседей (и, через него, совпадение вершины куска, Г50). */
static int pt_children(hz_ptree *t, int32_t nid) {
  double mid[3];
  for (int c = 0; c < 3; c++)
    mid[c] = 0.5 * (t->nd[nid].lo[c] + t->nd[nid].hi[c]);
  if (pt_grow_cut(t, 8) != 0) return 1;
  int32_t c0 = t->nnd;
  t->nd[nid].child = c0;
  t->nd[nid].t0 = 0;
  t->nd[nid].ntri = 0;
  t->nnd += 8;
  for (int k = 0; k < 8; k++) {
    hz_ptnode *C = &t->nd[c0 + k];
    for (int c = 0; c < 3; c++) {
      C->lo[c] = (k & (1 << c)) ? mid[c] : t->nd[nid].lo[c];
      C->hi[c] = (k & (1 << c)) ? t->nd[nid].hi[c] : mid[c];
    }
    C->child = -1;
    C->t0 = 0;
    C->ntri = 0;
    t->lev[c0 + k] = (unsigned char)(t->lev[nid] + 1);
  }
  return 0;
}

/* Задевает ли коробка треугольника коробку узла. Соглашение `[lo, hi)`
 * распространено и сюда (А475): касание по плоскости `hi` уходит соседу. */
static int pt_hits(const double *bl, const double *bh, const double *clo, const double *chi) {
  for (int c = 0; c < 3; c++)
    if (!(bl[c] < chi[c]) || !(bh[c] >= clo[c])) return 0;
  return 1;
}

static int pt_inside(const double *bl, const double *bh, const double *clo, const double *chi) {
  for (int c = 0; c < 3; c++)
    if (bl[c] < clo[c] || bh[c] > chi[c]) return 0;
  return 1;
}

static int pt_cut_split(hz_ptree *t, int32_t nid, int32_t *list, int32_t n, int lev,
                        const double *tlo, const double *thi) {
  t->nvisited++;
  if ((int32_t)t->lev[nid] > t->depth) t->depth = (int32_t)t->lev[nid];
  if (n <= t->leafmax || lev >= t->maxlev) {
    t->nd[nid].child = -1;
    t->nd[nid].t0 = 0;
    t->nd[nid].ntri = 0;
    t->nleaf++;
    /* Сколько ССЫЛОК было бы у листьев, если бы их хранили: это и есть
     * избыточность дублирования (§95.1), измеренная без её оплаты. */
    t->n_leaf_tri += n;
    return 0;
  }
  if (pt_children(t, nid) != 0) return 1;
  int32_t c0 = t->nd[nid].child;
  double clo[8][3], chi[8][3];
  for (int k = 0; k < 8; k++)
    for (int c = 0; c < 3; c++) {
      clo[k][c] = t->nd[c0 + k].lo[c];
      chi[k][c] = t->nd[c0 + k].hi[c];
    }
  /* СЧЁТ «НЕ ВЛЕЗ НИ В ОДНОГО РЕБЁНКА» — по САМОМУ ГЛУБОКОМУ СОДЕРЖАЩЕМУ узлу,
   * то есть ровно один раз на треугольник. Это те же треугольники, что при
   * старой укладке оставались у внутренних узлов (§193: 71…79 %). */
  for (int32_t i = 0; i < n; i++) {
    const double *bl = tlo + 3 * (size_t)list[i], *bh = thi + 3 * (size_t)list[i];
    if (!pt_inside(bl, bh, t->nd[nid].lo, t->nd[nid].hi)) continue;
    int fits = 0;
    for (int k = 0; k < 8 && !fits; k++)
      if (pt_inside(bl, bh, clo[k], chi[k])) fits = 1;
    if (!fits) t->n_nofit++;
  }
  int32_t *sub = malloc((size_t)n * sizeof *sub);
  if (sub == NULL) return 1;
  for (int k = 0; k < 8; k++) {
    int32_t ns = 0;
    for (int32_t i = 0; i < n; i++) {
      const double *bl = tlo + 3 * (size_t)list[i], *bh = thi + 3 * (size_t)list[i];
      if (pt_hits(bl, bh, clo[k], chi[k])) sub[ns++] = list[i];
    }
    if (pt_cut_split(t, c0 + k, sub, ns, lev + 1, tlo, thi) != 0) {
      free(sub);
      return 1;
    }
  }
  free(sub);
  return 0;
}

/* Найти узел, содержащий точку `p`, не спускаясь глубже уровня `lmax`.
 * Возврат `-1`, если точка вне корня. */
static int32_t pt_locate(const hz_ptree *t, const double *p, int lmax) {
  for (int c = 0; c < 3; c++)
    if (!(p[c] >= t->nd[0].lo[c]) || !(p[c] < t->nd[0].hi[c])) return -1;
  int32_t nid = 0;
  while (t->nd[nid].child >= 0 && (int)t->lev[nid] < lmax) {
    int32_t c0 = t->nd[nid].child, nx = -1;
    for (int k = 0; k < 8 && nx < 0; k++) {
      int in = 1;
      for (int c = 0; c < 3 && in; c++)
        if (!(p[c] >= t->nd[c0 + k].lo[c]) || !(p[c] < t->nd[c0 + k].hi[c])) in = 0;
      if (in) nx = c0 + k;
    }
    if (nx < 0) break;
    nid = nx;
  }
  return nid;
}

/* ГРАДУИРОВКА 2:1 (§241.3, А476). Довод не общий, а от главного такта:
 * перечисление соседей через грань есть внутренний цикл фронта, и без
 * градуировки оно неограниченно (`4^Δ` вместо `4`).
 *
 * Проба — ЦЕНТР ГРАНИ, снесённый на другую сторону плоскости. Этого довольно:
 * нарушение бывает только когда сосед КРУПНЕЕ, а крупный сосед накрывает всю
 * грань целиком, значит и её центр. Сдвиг делается `nextafter`, а не малым
 * числом: порога здесь быть не должно. */
static int pt_grade(hz_ptree *t) {
  int32_t *wl = malloc((size_t)t->nnd * sizeof *wl);
  if (wl == NULL) return 1;
  int64_t nw = 0, wcap = t->nnd;
  for (int32_t i = 0; i < t->nnd; i++)
    if (t->nd[i].child < 0) wl[nw++] = i;
  while (nw > 0) {
    int32_t L = wl[--nw];
    if (t->nd[L].child >= 0) continue;
    int l = (int)t->lev[L];
    if (l < 2) continue;
    for (int a = 0; a < 3 && t->nd[L].child < 0; a++) {
      for (int s = 0; s < 2; s++) {
        double p[3];
        for (int c = 0; c < 3; c++)
          p[c] = 0.5 * (t->nd[L].lo[c] + t->nd[L].hi[c]);
        p[a] = s ? t->nd[L].hi[a] : nextafter(t->nd[L].lo[a], -HUGE_VAL);
        int32_t N = pt_locate(t, p, l - 1);
        if (N < 0 || t->nd[N].child >= 0) continue;
        if ((int)t->lev[N] > l - 2) continue;
        if (l - (int)t->lev[N] > t->ngpass) t->ngpass = (int32_t)(l - (int)t->lev[N]);
        if (pt_children(t, N) != 0) {
          free(wl);
          return 1;
        }
        t->ngrade += 8;
        t->nleaf += 7;
        int32_t c0 = t->nd[N].child;
        if (nw + 9 > wcap) {
          int64_t nc = wcap * 2 + 16;
          int32_t *w2 = realloc(wl, (size_t)nc * sizeof *w2);
          if (w2 == NULL) {
            free(wl);
            return 1;
          }
          wl = w2;
          wcap = nc;
        }
        for (int k = 0; k < 8; k++)
          wl[nw++] = c0 + k;
        wl[nw++] = L;
        break;
      }
    }
  }
  free(wl);
  return 0;
}

int hz_ptree_build_cut(hz_ptree *t, const hz_objmesh *m, int leafmax, int maxlev, int grade) {
  memset(t, 0, sizeof *t);
  t->cut = 1;
  t->leafmax = (leafmax > 0) ? leafmax : 16;
  t->maxlev = (maxlev > 0) ? maxlev : 12;
  if (m->nt <= 0) return 1;
  double *tlo = malloc(3 * (size_t)m->nt * sizeof *tlo);
  double *thi = malloc(3 * (size_t)m->nt * sizeof *thi);
  int32_t *list = malloc((size_t)m->nt * sizeof *list);
  if (tlo == NULL || thi == NULL || list == NULL) {
    free(tlo);
    free(thi);
    free(list);
    return 2;
  }
  for (int32_t k = 0; k < m->nt; k++) {
    pt_tri_bb(m, k, tlo + 3 * (size_t)k, thi + 3 * (size_t)k);
    list[k] = k;
  }
  if (pt_grow_cut(t, 1) != 0) {
    free(tlo);
    free(thi);
    free(list);
    return 2;
  }
  t->nnd = 1;
  for (int c = 0; c < 3; c++) {
    t->nd[0].lo[c] = m->lo[c] - 1e-6;
    t->nd[0].hi[c] = m->hi[c] + 1e-6;
  }
  t->nd[0].child = -1;
  t->nd[0].t0 = 0;
  t->nd[0].ntri = 0;
  t->lev[0] = 0;
  double t0 = pt_now();
  int rc = pt_cut_split(t, 0, list, m->nt, 0, tlo, thi);
  t->t_build = pt_now() - t0;
  free(tlo);
  free(thi);
  free(list);
  if (rc != 0) return 2;
  t0 = pt_now();
  if (grade && pt_grade(t) != 0) return 2;
  t->t_grade = pt_now() - t0;
  t->redundancy = (double)t->n_leaf_tri / (double)m->nt;
  return 0;
}

void hz_ptree_free(hz_ptree *t) {
  free(t->nd);
  free(t->ref);
  free(t->lev);
  memset(t, 0, sizeof *t);
}
