/* pmarch.c — марш луча по ячейкам среза. Разбор — в `pmarch.h`. */
#include "pmarch.h"
#include <math.h>

/* Сравнения с нулём и на равенство БЕЗ `==`: гейт держит `-Wfloat-equal`, а
 * равенство здесь именно точное, а не с допуском. */
static int pm_iszero(double x) {
  return !(x < 0.0) && !(x > 0.0);
}
static int pm_eq(double a, double b) {
  return !(a < b) && !(a > b);
}

/* Взят ли узел в срез. Правило пола — то же, что у О70 (§194): угловой размер
 * `(2·rad/R)/ε`, спуск пока он больше пола; глаз ВНУТРИ узла даёт спуск до
 * предела, и это не особый случай, а тот же критерий. */
static int pm_at_cut(const hz_pmarch *M, int32_t nid) {
  const hz_ptnode *N = &M->T->nd[nid];
  if (N->child < 0) return 1;
  if (M->cut.fixlev >= 0) return (int)M->T->lev[nid] >= M->cut.fixlev;
  if (!(M->cut.pxcut > 0.0)) return 0;
  double r2 = 0.0, d2 = 0.0;
  for (int c = 0; c < 3; c++) {
    double h = 0.5 * (N->hi[c] - N->lo[c]);
    double dd = 0.5 * (N->lo[c] + N->hi[c]) - M->cut.eye[c];
    r2 += h * h;
    d2 += dd * dd;
  }
  double rad = sqrt(r2), R = sqrt(d2);
  if (!(R > rad)) return 0;
  return (2.0 * rad / R) / M->cut.eps <= M->cut.pxcut;
}

/* Спуск от корня к ячейке среза, содержащей точку `p`. Равенство с серединой
 * разрешается ПО ЗНАКУ НАПРАВЛЕНИЯ (см. `pmarch.h`): это и делает покрытие луча
 * точным при `d < 0`. Возврат — узел; `M->visits` — сколько узлов пройдено. */
static int32_t pm_locate(hz_pmarch *M, const double *p) {
  int32_t nid = 0;
  M->visits = 1;
  M->byleaf = 0;
  while (!pm_at_cut(M, nid)) {
    const hz_ptnode *N = &M->T->nd[nid];
    int k = 0;
    for (int c = 0; c < 3; c++) {
      double mid = 0.5 * (N->lo[c] + N->hi[c]);
      int up = (p[c] > mid) || (pm_eq(p[c], mid) && !(M->d[c] < 0.0));
      if (up) k |= 1 << c;
    }
    nid = N->child + k;
    M->visits++;
  }
  M->byleaf = (M->T->nd[nid].child < 0) ? 1 : 0;
  return nid;
}

/* Отрезок луча внутри коробки. Возврат `0`, если не задевает. */
static int pm_slab(const double *o, const double *d, const double *lo, const double *hi, double *t0,
                   double *t1) {
  double a = *t0, b = *t1;
  for (int c = 0; c < 3; c++) {
    if (pm_iszero(d[c])) {
      if (o[c] < lo[c] || !(o[c] < hi[c])) return 0;
      continue;
    }
    double x = (lo[c] - o[c]) / d[c], y = (hi[c] - o[c]) / d[c];
    if (x > y) {
      double s = x;
      x = y;
      y = s;
    }
    if (x > a) a = x;
    if (y < b) b = y;
    if (a > b) return 0;
  }
  *t0 = a;
  *t1 = b;
  return 1;
}

/* Выход из коробки текущей ячейки: наименьшее `t` по осям. Ось выхода
 * возвращается маской — их может быть несколько сразу (луч в ребро ячейки), и
 * тогда координату надо ставить точной по КАЖДОЙ из них. */
static double pm_exit(const hz_pmarch *M, const hz_ptnode *N, int *amask) {
  double t = M->tout;
  *amask = 0;
  for (int c = 0; c < 3; c++) {
    if (pm_iszero(M->d[c])) continue;
    double b = (M->d[c] > 0.0) ? N->hi[c] : N->lo[c];
    double tb = (b - M->o[c]) / M->d[c];
    if (tb < t) t = tb;
  }
  for (int c = 0; c < 3; c++) {
    if (pm_iszero(M->d[c])) continue;
    double b = (M->d[c] > 0.0) ? N->hi[c] : N->lo[c];
    double tb = (b - M->o[c]) / M->d[c];
    if (pm_eq(tb, t)) *amask |= 1 << c;
  }
  return t;
}

/* Точка на луче при `t`; по осям выхода координата ставится РАВНОЙ границе
 * точно, а не получается счётом (иначе округление вернуло бы точку в покинутую
 * ячейку). */
static void pm_point(const hz_pmarch *M, const hz_ptnode *N, double t, int amask, double *p) {
  for (int c = 0; c < 3; c++)
    p[c] = M->o[c] + t * M->d[c];
  for (int c = 0; c < 3; c++)
    if (amask & (1 << c)) p[c] = (M->d[c] > 0.0) ? N->hi[c] : N->lo[c];
}

int hz_pmarch_begin(hz_pmarch *M, const hz_ptree *T, const hz_pmarch_cut *cut, const double *o,
                    const double *d) {
  M->T = T;
  M->cut = *cut;
  for (int c = 0; c < 3; c++) {
    M->o[c] = o[c];
    M->d[c] = d[c];
  }
  M->nstep = 0;
  M->nvisit = 0;
  M->ncap = 0;
  M->nbreak = 0;
  M->visits = 0;
  M->nid = -1;
  M->lev = 0;
  M->t0 = 0.0;
  M->t1 = 0.0;
  double a = 0.0, b = HUGE_VAL;
  if (!pm_slab(M->o, M->d, T->nd[0].lo, T->nd[0].hi, &a, &b)) return 0;
  M->tin = a;
  M->tout = b;
  double p[3];
  for (int c = 0; c < 3; c++)
    p[c] = M->o[c] + a * M->d[c];
  M->nid = pm_locate(M, p);
  M->lev = (int)T->lev[M->nid];
  M->nvisit += M->visits;
  int amask = 0;
  M->t0 = a;
  M->t1 = pm_exit(M, &T->nd[M->nid], &amask);
  M->nstep = 1;
  return 1;
}

int hz_pmarch_next(hz_pmarch *M) {
  if (M->nid < 0) return 0;
  if (!(M->t1 < M->tout)) return 0;
  if (M->nstep >= HZ_PMARCH_MAXSTEP) {
    M->ncap++;
    return 0;
  }
  const hz_ptnode *N = &M->T->nd[M->nid];
  int amask = 0;
  double t = pm_exit(M, N, &amask);
  double p[3];
  pm_point(M, N, t, amask, p);
  M->nid = pm_locate(M, p);
  M->lev = (int)M->T->lev[M->nid];
  M->nvisit += M->visits;
  /* ПОКРЫТИЕ ЛУЧА ПРОВЕРЯЕТСЯ ЗДЕСЬ, И ПРОВЕРКА ВЫБРАНА НЕ САМАЯ УДОБНАЯ, А
   * ТА, ЧТО ЛОВИТ ПРОПУСК (А498). Слабая форма — «сумма длин шагов равна длине
   * луча» — прошла бы и при пропущенной ячейке нулевой протяжённости. Здесь
   * требуется большее: у НАЙДЕННОЙ ячейки вход по её собственной коробке обязан
   * совпасть с `t` ТОЧНО. Совпасть он обязан побитово: общая грань соседей —
   * одна и та же дважды делённая пополам координата, а деление пополам двоичных
   * значений точно. Расхождение означало бы либо щель, либо перекрытие. */
  const hz_ptnode *X = &M->T->nd[M->nid];
  double a = -HUGE_VAL, b = HUGE_VAL;
  if (!pm_slab(M->o, M->d, X->lo, X->hi, &a, &b) || !pm_eq(a, t)) M->nbreak++;
  M->t0 = t;
  M->t1 = pm_exit(M, X, &amask);
  M->nstep++;
  return 1;
}
