/* plink.c — сборка и решение оператора связями. Разбор — в `plink.h`. */
#include "plink.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define LK_PI 3.14159265358979323846

/* --- списки детей по лестнице ---------------------------------------------- */

/* Лестница хранит РОДИТЕЛЯ у каждого узла; для дробления нужны ДЕТИ. Строятся
 * счётной сортировкой: один проход посчитать, второй разложить. */
typedef struct {
  int32_t *start; /* nnd+1 */
  int32_t *ch;    /* дети подряд */
} lk_kids;

static int lk_kids_build(lk_kids *K, const hz_lod *L) {
  K->start = calloc((size_t)L->nnd + 1, sizeof *K->start);
  if (K->start == NULL) return 1;
  int32_t nch = 0;
  for (int32_t k = 0; k < L->nnd; k++) {
    int32_t p = L->nd[k].parent;
    if (p >= 0 && p < L->nnd) {
      K->start[p + 1]++;
      nch++;
    }
  }
  for (int32_t k = 0; k < L->nnd; k++)
    K->start[k + 1] += K->start[k];
  K->ch = malloc((size_t)(nch > 0 ? nch : 1) * sizeof *K->ch);
  if (K->ch == NULL) {
    free(K->start);
    return 1;
  }
  int32_t *pos = malloc((size_t)L->nnd * sizeof *pos);
  if (pos == NULL) {
    free(K->start);
    free(K->ch);
    return 1;
  }
  memcpy(pos, K->start, (size_t)L->nnd * sizeof *pos);
  for (int32_t k = 0; k < L->nnd; k++) {
    int32_t p = L->nd[k].parent;
    if (p >= 0 && p < L->nnd) K->ch[pos[p]++] = k;
  }
  free(pos);
  return 0;
}

static void lk_kids_free(lk_kids *K) {
  free(K->start);
  free(K->ch);
  memset(K, 0, sizeof *K);
}

/* --- опорная точка ПЕРЕНОСА ------------------------------------------------- */

/* ЗАЧЕМ ОТДЕЛЬНАЯ ТОЧКА (§88, А210). Узел несёт ЦЕНТР ПЛОЩАДИ, и срезу камеры
 * его довольно — там нужна оценка расстояния. Переносу он не годится: участок
 * сегментации бывает невыпуклым и несвязным, и центр площади тогда лежит ВНЕ
 * поверхности. Замерено на зале: 37.4 % ПЛОЩАДИ, и видимость из таких точек
 * давала `Σf = 0.34` там, где ссылка лучами требует `0.98`.
 *
 * ЦЕНТР ПЛОЩАДИ В УЗЛЕ ПРИ ЭТОМ НЕ ТРОГАЕТСЯ: на нём сняты все замеры среза, и
 * правка обесценила бы точку отсчёта (класс А1).
 *
 * ВЫБОР — ЦЕНТР САМОГО КРУПНОГО ТРЕУГОЛЬНИКА УЧАСТКА, и он выбран ЗАМЕРОМ, а не
 * рассуждением. Первым правилом была «ближайшая к центру площади точка
 * поверхности» — она сохраняла первый момент, и А215 назвала её риск заранее:
 * у невыпуклого участка ближайшая точка ложится на КРОМКУ. Прогон это
 * подтвердил жёстче предсказанного: доля площади с точкой вне полигона упала
 * лишь 37.4 -> 34.4 %, а первая процентиль расстояния до попадания РУХНУЛА
 * `2.31e-3 -> 1.07e-4` м, то есть лучи стали упираться в геометрию мгновенно —
 * точка села на стык. Центр треугольника строго внутренний по построению и от
 * кромки отстоит соразмерно самому треугольнику. Первый момент при этом
 * приносится в жертву сознательно: точка на стыке мерит не тот свет вовсе, а
 * смещение момента ограничено размером элемента.
 *
 * Проекция на плоскость полигона с двумерным тестом принадлежности была бы
 * дешевле и НЕВЕРНА: треугольники отстоят от общей плоскости на `dmax`, и точка
 * внутри края лежала бы на ПЛОСКОСТИ, а не на поверхности (А218, тот же класс,
 * что Г40/Г44). */
/* Точки переноса для ВСЕХ узлов лестницы. Уровень 0 — ближайшая к центру площади
 * точка объединения треугольников участка; внутренний узел — точка того из его
 * полигонов, чья ближе к центру площади узла. `pt` — на `3·nnd`. */
int hz_links_points(double *pt, const hz_lod *L, const hz_polyset *ps, const hz_objmesh *m,
                    int oldpt) {
  const int32_t nn = L->nnd;
  for (int32_t k = 0; k < nn; k++) {
    pt[3 * k] = L->nd[k].cx;
    pt[3 * k + 1] = L->nd[k].cy;
    pt[3 * k + 2] = L->nd[k].cz;
  }
  if (ps == NULL || m == NULL || oldpt) return 0;
  double *pt0 = malloc(3 * (size_t)(ps->np > 0 ? ps->np : 1) * sizeof *pt0);
  if (pt0 == NULL) return 1;
  for (int32_t k = 0; k < ps->np; k++) {
    const hz_poly *p = &ps->p[k];
    pt0[3 * k] = p->org[0];
    pt0[3 * k + 1] = p->org[1];
    pt0[3 * k + 2] = p->org[2];
    double best = -1.0;
    for (int32_t t = p->t0; t < p->t0 + p->ntri; t++) {
      int32_t tr = ps->tri[t];
      const double *a = m->v + 3 * m->f[3 * tr];
      const double *b = m->v + 3 * m->f[3 * tr + 1];
      const double *c = m->v + 3 * m->f[3 * tr + 2];
      double u[3], v[3], x[3];
      for (int i = 0; i < 3; i++) {
        u[i] = b[i] - a[i];
        v[i] = c[i] - a[i];
      }
      x[0] = u[1] * v[2] - u[2] * v[1];
      x[1] = u[2] * v[0] - u[0] * v[2];
      x[2] = u[0] * v[1] - u[1] * v[0];
      double ar2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
      if (ar2 > best) {
        best = ar2;
        for (int i = 0; i < 3; i++)
          pt0[3 * k + i] = (a[i] + b[i] + c[i]) / 3.0;
      }
    }
  }
  double *bd = malloc((size_t)nn * sizeof *bd);
  if (bd == NULL) {
    free(pt0);
    return 1;
  }
  for (int32_t k = 0; k < nn; k++)
    bd[k] = 1e300;
  for (int32_t lev = 0; lev < L->nlev; lev++)
    for (int32_t k = 0; k < L->np && k < ps->np; k++) {
      int32_t nd = L->lab[(size_t)lev * (size_t)L->np + (size_t)k];
      if (nd < 0 || nd >= nn) continue;
      double d2 = 0.0;
      double q[3] = {L->nd[nd].cx, L->nd[nd].cy, L->nd[nd].cz};
      for (int i = 0; i < 3; i++)
        d2 += (pt0[3 * k + i] - q[i]) * (pt0[3 * k + i] - q[i]);
      if (d2 < bd[nd]) {
        bd[nd] = d2;
        for (int i = 0; i < 3; i++)
          pt[3 * nd + i] = pt0[3 * k + i];
      }
    }
  free(bd);
  free(pt0);
  return 0;
}

/* --- геометрия связи -------------------------------------------------------- */

/* КОЭФФИЦИЕНТ ТОЧКА-ДИСК. Внутренний узел лестницы полигона не несёт — у него есть
 * опорная точка, нормаль и площадь, — поэтому излучатель приближается диском той же
 * площади. Слагаемое `A_j` в знаменателе не подгонка, а регуляризация: без него
 * выражение расходится при сближении, а с ним при `r² ≫ A_j` оно переходит в точную
 * формулу, при `r → 0` остаётся конечным. Ошибка приближения ограничена ровно тем
 * критерием, по которому идёт дробление: чем меньше угловой размер, тем точнее диск.
 *
 *     F = cosθ_i · cosθ_j · A_j / (π r² + A_j)
 */
static double lk_coef(const hz_lodnode *A, const double xa[3], const hz_lodnode *B,
                      const double xb[3]) {
  double d[3] = {xb[0] - xa[0], xb[1] - xa[1], xb[2] - xa[2]};
  double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  if (!(r2 > 0.0)) return 0.0;
  double r = sqrt(r2);
  double ci = (d[0] * A->n[0] + d[1] * A->n[1] + d[2] * A->n[2]) / r;
  double cj = -(d[0] * B->n[0] + d[1] * B->n[1] + d[2] * B->n[2]) / r;
  if (!(ci > 0.0) || !(cj > 0.0)) return 0.0; /* друг друга не видят по ориентации */
  return ci * cj * B->area_surf / (LK_PI * r2 + B->area_surf);
}

/* ТОЧНЫЙ ФОРМ-ФАКТОР ТОЧКА–МНОГОУГОЛЬНИК (§87), контурный интеграл Ламберта:
 *
 *     F = (1/2π) · Σ_рёбра  β_k · (n_i · e_k),
 *     β_k = ∠(R_k, R_{k+1}),   e_k = (R_k × R_{k+1}) / |R_k × R_{k+1}|,
 *
 * где `R_k` — вершины края излучателя ОТНОСИТЕЛЬНО приёмной точки. Замкнутая
 * форма, точная при любой близости; квадратуры нет. Заменяет точечно-дисковую
 * оценку там, где та применялась хуже всего, — на листьях, где `√A_j / r ≈ 1`.
 *
 * ОТСЕЧЕНИЕ ПОЛУПЛОСКОСТЬЮ ПРИЁМНИКА ОБЯЗАТЕЛЬНО: формула точна лишь для части
 * излучателя НАД касательной плоскостью, а пересекающий её полигон даёт вклад
 * частично отрицательный. Рёбра, легшие в саму плоскость, вклад дают и он
 * законен — это дуга горизонта, входящая в проектированный телесный угол (А209).
 *
 * ПРИЁМНИК ОСТАЁТСЯ ТОЧКОЙ. Ошибка приёмной стороны этим не лечится, и её
 * остаток обязан быть виден как перебор `Σf` у близких пар. */
static double lk_ff_loop(const double *R, int n, const double nv[3]) {
  double s = 0.0;
  for (int k = 0; k < n; k++) {
    const double *a = R + 3 * k, *b = R + 3 * ((k + 1) % n);
    double c[3] = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
    double cl = sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
    if (!(cl > 0.0)) continue;
    double la = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
    double lb = sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2]);
    if (!(la > 0.0) || !(lb > 0.0)) continue;
    double ct = (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]) / (la * lb);
    if (ct > 1.0) ct = 1.0;
    if (ct < -1.0) ct = -1.0;
    double beta = acos(ct);
    s += beta * (nv[0] * c[0] + nv[1] * c[1] + nv[2] * c[2]) / cl;
  }
  /* ЗНАК. Внешняя петля полигона обходится против часовой стрелки В ЕГО
   * СОБСТВЕННОЙ раме (`polygon.h`, `eu × ev = n`), а приёмник смотрит на лицевую
   * сторону, то есть видит ту же петлю ПО часовой. Классическая запись формулы
   * ждёт обход против часовой ИЗ ПРИЁМНИКА, отсюда минус. Установлено
   * самопроверкой: без него оба аналитических случая выходили с точностью до
   * 15-й цифры, но отрицательными. */
  return -s / (2.0 * LK_PI);
}

/* Отсечение петли полуплоскостью `nv·R > 0` по Сазерленду–Ходжмену. Возвращает
 * число вершин на выходе; `out` обязан вмещать `n + 1`. */
static int lk_clip_half(const double *R, int n, const double nv[3], double *out) {
  int m = 0;
  for (int k = 0; k < n; k++) {
    const double *a = R + 3 * k, *b = R + 3 * ((k + 1) % n);
    double da = nv[0] * a[0] + nv[1] * a[1] + nv[2] * a[2];
    double db = nv[0] * b[0] + nv[1] * b[1] + nv[2] * b[2];
    if (da > 0.0) {
      for (int c = 0; c < 3; c++)
        out[3 * m + c] = a[c];
      m++;
    }
    if ((da > 0.0) != (db > 0.0)) {
      double t = da / (da - db);
      for (int c = 0; c < 3; c++)
        out[3 * m + c] = a[c] + t * (b[c] - a[c]);
      m++;
    }
  }
  return m;
}

/* Вклад ОДНОЙ петли края, заданной в местных `(u,v)` полигона. Мир —
 * `org + u·eu + v·ev` (`polygon.h`), приёмник вычитается, дальше отсечение и
 * контурная сумма. `w1`, `w2` — рабочие буферы на `3·(nv+1)` каждый. */
static double lk_ff_loop_uv(const double x[3], const double nrec[3], const double org[3],
                            const double eu[3], const double ev[3], const double *uv, int nv,
                            double *w1, double *w2) {
  if (nv < 3) return 0.0;
  for (int k = 0; k < nv; k++)
    for (int c = 0; c < 3; c++)
      w1[3 * k + c] = org[c] + uv[2 * k] * eu[c] + uv[2 * k + 1] * ev[c] - x[c];
  int m = lk_clip_half(w1, nv, nrec, w2);
  if (m < 3) return 0.0;
  return lk_ff_loop(w2, m, nrec);
}

/* САМОПРОВЕРКА ФОРМУЛЫ НА АНАЛИТИКЕ (А205). Проверять контурный интеграл СЦЕНОЙ
 * нельзя: ошибка в `2π`, потерянный косинус или перевёрнутый обход дали бы такое
 * же правдоподобное `Σf`, и отличить их было бы нечем. Здесь два случая с
 * известным ответом. Квадрат полуширины `a` на оси приёмника, на высоте `h`:
 *
 *     F = (4/π) · X/√(1+X²) · atan(X/√(1+X²)),   X = a/h
 *
 * (четыре угловых прямоугольника Хауэлла B-1). При `X → ∞` это ровно `1` —
 * излучатель закрывает всю полусферу. */
double hz_links_ff_selftest(double *f_big, double *f_unit, double *ref_unit) {
  double x[3] = {0.0, 0.0, 0.0}, nrec[3] = {0.0, 0.0, 1.0};
  /* Излучатель смотрит ВНИЗ, на приёмника: `n_j = (0,0,−1)`, и рама берётся с
   * `eu × ev = n_j`, как у настоящего полигона. */
  double eu[3] = {1.0, 0.0, 0.0}, ev[3] = {0.0, -1.0, 0.0};
  double w1[15], w2[18], worst = 0.0;
  const double h = 1.0;
  for (int cs = 0; cs < 2; cs++) {
    double a = (cs == 0) ? 1e6 : 1.0;
    double org[3] = {0.0, 0.0, h};
    double uv[8] = {-a, -a, a, -a, a, a, -a, a}; /* против часовой в раме (eu,ev) */
    double f = lk_ff_loop_uv(x, nrec, org, eu, ev, uv, 4, w1, w2);
    double X = a / h, t = X / sqrt(1.0 + X * X);
    double ref = (4.0 / LK_PI) * t * atan(t);
    if (cs == 0) {
      *f_big = f;
    } else {
      *f_unit = f;
      *ref_unit = ref;
    }
    double e = fabs(f - ref);
    if (e > worst) worst = e;
  }
  return worst;
}

/* Угловой размер `B` из `A` — та же величина, что в срезе камеры, только вместо
 * глаза стоит элемент. Один допуск на всю схему. */
static double lk_subtend(const double xa[3], const hz_lodnode *B, const double xb[3]) {
  double d[3] = {xb[0] - xa[0], xb[1] - xa[1], xb[2] - xa[2]};
  double r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  if (!(r2 > 0.0)) return 1e300;
  return sqrt(B->area_surf / r2);
}

/* ДОЛЯ НЕЗАСЛОНЁННОГО. Пробы ставятся на дисках узлов: точка сдвигается в
 * касательной плоскости на долю радиуса объемлющей сферы. Это грубее, чем пробы по
 * краю полигона, и намеренно: у внутреннего узла края нет вовсе, а у листа он есть,
 * но связь и так есть осреднение по паре площадей. */
static double lk_vis(const hz_scene *sc, int32_t ia, const double xa[3], int32_t ib,
                     const double xb[3], int nvis, int64_t *nray) {
  const hz_lodnode *A = &sc->L->nd[ia], *B = &sc->L->nd[ib];
  double eu[3], ev[3];
  int ax = 0;
  for (int c = 1; c < 3; c++)
    if (fabs(A->n[c]) < fabs(A->n[ax])) ax = c;
  double t0[3] = {0.0, 0.0, 0.0};
  t0[ax] = 1.0;
  eu[0] = A->n[1] * t0[2] - A->n[2] * t0[1];
  eu[1] = A->n[2] * t0[0] - A->n[0] * t0[2];
  eu[2] = A->n[0] * t0[1] - A->n[1] * t0[0];
  double en = sqrt(eu[0] * eu[0] + eu[1] * eu[1] + eu[2] * eu[2]);
  if (!(en > 0.0)) return 0.0;
  for (int c = 0; c < 3; c++)
    eu[c] /= en;
  ev[0] = A->n[1] * eu[2] - A->n[2] * eu[1];
  ev[1] = A->n[2] * eu[0] - A->n[0] * eu[2];
  ev[2] = A->n[0] * eu[1] - A->n[1] * eu[0];
  int hit = 0, tot = 0;
  for (int u = 0; u < nvis; u++)
    for (int v = 0; v < nvis; v++) {
      double su = (((double)u + 0.5) / nvis - 0.5) * 1.2 * A->rad;
      double sv = (((double)v + 0.5) / nvis - 0.5) * 1.2 * A->rad;
      double tu = (((double)v + 0.5) / nvis - 0.5) * 1.2 * B->rad;
      double tv = (((double)u + 0.5) / nvis - 0.5) * 1.2 * B->rad;
      double pa[3] = {xa[0] + su * eu[0] + sv * ev[0], xa[1] + su * eu[1] + sv * ev[1],
                      xa[2] + su * eu[2] + sv * ev[2]};
      double pb[3] = {xb[0] + tu * eu[0] + tv * ev[0], xb[1] + tu * eu[1] + tv * ev[1],
                      xb[2] + tu * eu[2] + tv * ev[2]};
      double d[3], len = 0.0;
      for (int c = 0; c < 3; c++) {
        d[c] = pb[c] - pa[c];
        len += d[c] * d[c];
      }
      len = sqrt(len);
      if (!(len > 0.0)) continue;
      double eps = 1e-6 * (1.0 + len), xo[3];
      for (int c = 0; c < 3; c++) {
        d[c] /= len;
        xo[c] = pa[c] + eps * A->n[c];
      }
      tot++;
      (*nray)++;
      /* ЭЛЕМЕНТ НЕ ЗАТЕНЯЕТ САМ СЕБЯ. Попадание в полигон, ПРИНАДЛЕЖАЩИЙ узлу
       * приёмника или узлу излучателя, заслонением не считается: внутри одного
       * элемента геометрия представлена плоской, и её подэлементный рельеф
       * (`dmax` участка, на зале до 0.045 м) относится к ОТКЛИКУ поверхности, а
       * не к взаимной видимости элементов. Без этого правила `Σf` перебором
       * давало `0.327` при ссылке лучами `0.980`; с ним — `1.058` (§88.2).
       * Проверяется через `lab`: полигон `ph` принадлежит узлу `nd` уровня
       * `lev`, если `lab[lev·np + ph] == nd`. */
      double th = 0.0;
      int32_t ph = hz_pray_hit(sc->g, xo, d, 0.0, &th);
      int blocked = (ph >= 0 && th < len * (1.0 - 1e-5));
      if (blocked && sc->L != NULL && ph < sc->L->np) {
        const hz_lod *LL = sc->L;
        int32_t la = LL->lab[(size_t)LL->nd[ia].level * (size_t)LL->np + (size_t)ph];
        int32_t lb = LL->lab[(size_t)LL->nd[ib].level * (size_t)LL->np + (size_t)ph];
        if (la == ia || lb == ib) blocked = 0;
      }
      if (!blocked) hit++;
    }
  return (tot > 0) ? (double)hit / (double)tot : 0.0;
}

/* --- сборка ----------------------------------------------------------------- */

typedef struct {
  int32_t a, b;
} lk_pair;

typedef struct {
  lk_pair *p;
  int64_t n, cap;
} lk_queue;

static int lk_qpush(lk_queue *q, int32_t a, int32_t b) {
  if (q->n >= q->cap) {
    int64_t nc = (q->cap > 0) ? q->cap * 2 : 4096;
    lk_pair *np = realloc(q->p, (size_t)nc * sizeof *np);
    if (np == NULL) return 1;
    q->p = np;
    q->cap = nc;
  }
  q->p[q->n].a = a;
  q->p[q->n].b = b;
  q->n++;
  return 0;
}

static int lk_push(hz_linkset *S, int32_t i, int32_t j, double f) {
  if (S->n >= S->cap) {
    int64_t nc = (S->cap > 0) ? S->cap * 2 : 4096;
    hz_link *nl = realloc(S->l, (size_t)nc * sizeof *nl);
    if (nl == NULL) return 1;
    S->l = nl;
    S->cap = nc;
  }
  S->l[S->n].i = i;
  S->l[S->n].j = j;
  S->l[S->n].f = f;
  S->n++;
  return 0;
}

/* --- подпись замкнутости `Σf` ------------------------------------------------ */

typedef struct {
  double sf, ar;
} lk_sfa;

static int lk_cmp_sf(const void *a, const void *b) {
  double x = ((const lk_sfa *)a)->sf, y = ((const lk_sfa *)b)->sf;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

/* `Σ_j f_ij` по ЛИСТУ с учётом предков: грубая связь приносит энергию всем
 * потомкам, значит и в сумму листа входит. Проценты — по ПЛОЩАДИ (А194). */
void hz_links_sf(const hz_linkset *S, const hz_lod *L, double *sf) {
  const int32_t nn = L->nnd;
  for (int32_t k = 0; k < nn; k++)
    sf[k] = 0.0;
  for (int64_t l = 0; l < S->n; l++)
    sf[S->l[l].i] += S->l[l].f;
  for (int32_t lev = L->nlev - 1; lev >= 1; lev--)
    for (int32_t k = 0; k < nn; k++) {
      if (L->nd[k].level != lev - 1) continue;
      int32_t p = L->nd[k].parent;
      if (p >= 0 && p < nn) sf[k] += sf[p];
    }
}

static void lk_closure(hz_linkset *S, const hz_lod *L) {
  const int32_t nn = L->nnd;
  double *sf = calloc((size_t)(nn > 0 ? nn : 1), sizeof *sf);
  lk_sfa *v = malloc((size_t)(nn > 0 ? nn : 1) * sizeof *v);
  if (sf == NULL || v == NULL) {
    free(sf);
    free(v);
    return;
  }
  hz_links_sf(S, L, sf);
  int32_t m = 0;
  double atot = 0.0, sum = 0.0, alow = 0.0;
  S->sf_min = 1e300;
  S->sf_max = 0.0;
  for (int32_t k = 0; k < nn; k++) {
    if (L->nd[k].level != 0) continue;
    double a = L->nd[k].area_surf;
    v[m].sf = sf[k];
    v[m].ar = a;
    m++;
    atot += a;
    sum += sf[k] * a;
    if (sf[k] < 0.9) alow += a;
    if (sf[k] < S->sf_min) S->sf_min = sf[k];
    if (sf[k] > S->sf_max) S->sf_max = sf[k];
  }
  if (m == 0 || !(atot > 0.0)) {
    S->sf_min = 0.0;
    free(sf);
    free(v);
    return;
  }
  S->sf_mean = sum / atot;
  S->sf_lowfrac = alow / atot;
  qsort(v, (size_t)m, sizeof *v, lk_cmp_sf);
  const double qs[3] = {0.10, 0.50, 0.90};
  double *out[3] = {&S->sf_p10, &S->sf_p50, &S->sf_p90};
  int qi = 0;
  double cum = 0.0;
  for (int32_t k = 0; k < m && qi < 3; k++) {
    cum += v[k].ar;
    while (qi < 3 && cum >= qs[qi] * atot) {
      *out[qi] = v[k].sf;
      qi++;
    }
  }
  while (qi < 3) {
    *out[qi] = v[m - 1].sf;
    qi++;
  }
  free(sf);
  free(v);
}

/* Точный форм-фактор от опорной точки узла `A` к полигону `p`: сумма по ВСЕМ
 * петлям края. Внешняя петля даёт плюс, дыра (обход по часовой) — минус, и
 * особого случая для дыр не нужно. Отрицательный итог значит «отвёрнут». */
static double lk_coef_poly(const hz_lodnode *A, const double xa[3], const hz_polyset *ps,
                           int32_t ip, int noclip, double *w1, double *w2) {
  const hz_poly *p = &ps->p[ip];
  double x[3] = {xa[0], xa[1], xa[2]};
  double s = 0.0;
  for (int32_t li = 0; li < p->nloop; li++) {
    int32_t b0 = ps->loop[p->l0 + li], b1 = ps->loop[p->l0 + li + 1];
    int nv = (int)(b1 - b0);
    if (nv < 3) continue;
    if (noclip) {
      /* НК4: без отсечения. Часть излучателя под касательной плоскостью даёт
       * вклад со своим знаком, и величина перестаёт быть форм-фактором. */
      for (int k = 0; k < nv; k++)
        for (int c = 0; c < 3; c++)
          w1[3 * k + c] = p->org[c] + ps->bv[2 * (b0 + k)] * p->eu[c] +
                          ps->bv[2 * (b0 + k) + 1] * p->ev[c] - x[c];
      s += lk_ff_loop(w1, nv, A->n);
    } else {
      s += lk_ff_loop_uv(x, A->n, p->org, p->eu, p->ev, ps->bv + 2 * b0, nv, w1, w2);
    }
  }
  return (s > 0.0) ? s : 0.0;
}

/* ЭТАЛОН ПОЛНЫМ ПЕРЕБОРОМ (§87.4). Отвечает на вопрос, которого не различают ни
 * `Σf`, ни ссылка лучами по отдельности: ТЕРЯЕТ ЛИ ЭНЕРГИЮ ИЕРАРХИЯ или сам
 * коэффициент. Для выборки листьев `Σf` считается по ВСЕМ листьям сцены без
 * всякого дробления — тем же точным форм-фактором и той же видимостью. Если
 * перебор даёт около единицы, а иерархия `0.65`, потеря в дроблении; если и
 * перебор даёт `0.65`, дело в коэффициенте либо в точечной коллокации. */
int hz_links_sf_brute(const hz_scene *sc, int stride, int nvis, int oldpt, double *sf_novis,
                      double *sf_vis, double *sf_selfok, int32_t *node, int cap) {
  const hz_lod *L = sc->L;
  if (sc->ps == NULL || stride < 1) return 0;
  if (nvis < 1) nvis = 2;
  int32_t mxv = 3;
  for (int32_t k = 0; k < sc->ps->np; k++) {
    const hz_poly *p = &sc->ps->p[k];
    for (int32_t li = 0; li < p->nloop; li++) {
      int32_t nv = sc->ps->loop[p->l0 + li + 1] - sc->ps->loop[p->l0 + li];
      if (nv > mxv) mxv = nv;
    }
  }
  double *w1 = malloc(3 * (size_t)(mxv + 2) * sizeof *w1);
  double *w2 = malloc(3 * (size_t)(mxv + 2) * sizeof *w2);
  int32_t *p2n = malloc((size_t)sc->ps->np * sizeof *p2n);
  double *pt = malloc(3 * (size_t)L->nnd * sizeof *pt);
  if (w1 == NULL || w2 == NULL || p2n == NULL || pt == NULL ||
      hz_links_points(pt, L, sc->ps, sc->m, oldpt) != 0) {
    free(w1);
    free(w2);
    free(p2n);
    free(pt);
    return 0;
  }
  for (int32_t k = 0; k < L->np && k < sc->ps->np; k++)
    p2n[k] = L->lab[k];
  int n = 0, seen = 0;
  int64_t dummy = 0;
  for (int32_t k = 0; k < L->nnd && n < cap; k++) {
    if (L->nd[k].level != 0) continue;
    if (seen++ % stride != 0) continue;
    const hz_lodnode *A = &L->nd[k];
    double s0 = 0.0, s1 = 0.0, s2 = 0.0;
    for (int32_t j = 0; j < sc->ps->np; j++) {
      if (p2n[j] == k) continue;
      double f = lk_coef_poly(A, pt + 3 * k, sc->ps, j, 0, w1, w2);
      if (!(f > 0.0)) continue;
      s0 += f;
      s1 += f * lk_vis(sc, k, pt + 3 * k, p2n[j], pt + 3 * p2n[j], nvis, &dummy);
      /* ТРЕТЬЯ СУММА: то же, но ПЕРВОЕ попадание в СВОЙ участок приёмника или в
       * сам излучатель считается НЕ заслонением. Внутри одного элемента
       * геометрия представлена плоской, и её подэлементный рельеф (до `dmax`,
       * на зале 0.045 м) относится к отклику поверхности, а не к взаимной
       * видимости элементов. Разность второй и третьей сумм и есть цена того,
       * что элемент затеняет сам себя. */
      double xa2[3] = {pt[3 * k], pt[3 * k + 1], pt[3 * k + 2]};
      double xb2[3] = {pt[3 * p2n[j]], pt[3 * p2n[j] + 1], pt[3 * p2n[j] + 2]};
      double dd[3], ln = 0.0;
      for (int c = 0; c < 3; c++) {
        dd[c] = xb2[c] - xa2[c];
        ln += dd[c] * dd[c];
      }
      ln = sqrt(ln);
      if (!(ln > 0.0)) continue;
      double ori[3];
      for (int c = 0; c < 3; c++) {
        dd[c] /= ln;
        ori[c] = xa2[c] + 1e-6 * (1.0 + ln) * A->n[c];
      }
      double th = 0.0;
      int32_t ph = hz_pray_hit(sc->g, ori, dd, 0.0, &th);
      if (ph < 0 || th >= ln * (1.0 - 1e-5) || ph == k || ph == j) s2 += f;
    }
    sf_novis[n] = s0;
    sf_vis[n] = s1;
    sf_selfok[n] = s2;
    node[n] = k;
    n++;
  }
  free(w1);
  free(w2);
  free(p2n);
  free(pt);
  return n;
}

int hz_links_build(hz_linkset *S, const hz_scene *sc, const hz_linkcfg *cfg) {
  memset(S, 0, sizeof *S);
  const double eps = cfg->eps;
  const int noself = cfg->noself;
  int nvis = cfg->nvis;
  if (nvis < 1) nvis = 2;
  const hz_lod *L = sc->L;
  lk_kids K;
  if (lk_kids_build(&K, L) != 0) return 2;

  /* КАРТА «УЗЕЛ УРОВНЯ 0 → ПОЛИГОН». Нулевой уровень лестницы и есть участки
   * сегментации, соответствие взаимно однозначно и читается из нулевой строки
   * `lab`. Рабочие буферы отсечения — в куче и один раз на сборку (А208): край
   * доходит до 1684 вершин, а `-Wvla` запрещает массив переменной длины. */
  int32_t *np2poly = NULL;
  double *w1 = NULL, *w2 = NULL;
  lk_queue q;
  memset(&q, 0, sizeof q);
  int64_t *cnt = NULL;
  /* ОПОРНЫЕ ТОЧКИ ПЕРЕНОСА (§88). Строятся до всего: на них держатся и
   * отбраковка, и критерий дробления, и видимость, и коэффициент. */
  double *pt = malloc(3 * (size_t)L->nnd * sizeof *pt);
  if (pt == NULL || hz_links_points(pt, L, sc->ps, sc->m, cfg->oldpt) != 0) {
    free(pt);
    lk_kids_free(&K);
    return 2;
  }
  if (sc->ps != NULL && !cfg->disk) {
    np2poly = malloc((size_t)L->nnd * sizeof *np2poly);
    if (np2poly == NULL) {
      lk_kids_free(&K);
      return 2;
    }
    for (int32_t k = 0; k < L->nnd; k++)
      np2poly[k] = -1;
    int32_t mxv = 3;
    for (int32_t k = 0; k < L->np && k < sc->ps->np; k++) {
      np2poly[L->lab[k]] = k;
      const hz_poly *p = &sc->ps->p[k];
      for (int32_t li = 0; li < p->nloop; li++) {
        int32_t nv = sc->ps->loop[p->l0 + li + 1] - sc->ps->loop[p->l0 + li];
        if (nv > mxv) mxv = nv;
      }
    }
    w1 = malloc(3 * (size_t)(mxv + 2) * sizeof *w1);
    w2 = malloc(3 * (size_t)(mxv + 2) * sizeof *w2);
    if (w1 == NULL || w2 == NULL) goto fail;
  }

  /* КОРНИ — узлы САМОГО ГРУБОГО уровня, реально занятые. Иерархия начинается
   * сверху и дробится вниз: пара грубых узлов либо порождает связь, либо
   * распадается на пары их детей. Именно так число связей и выходит `O(n log n)`
   * вместо `O(n²)` — дальние пары остаются наверху. */
  int32_t top = L->nlev - 1;
  int32_t *isroot = calloc((size_t)L->nnd, sizeof *isroot);
  if (isroot == NULL) goto fail;
  for (int32_t k = 0; k < L->np; k++)
    isroot[L->lab[(size_t)top * (size_t)L->np + (size_t)k]] = 1;

  for (int32_t a = 0; a < L->nnd; a++) {
    if (!isroot[a]) continue;
    for (int32_t b = 0; b < L->nnd; b++)
      if (isroot[b] && (a != b || !noself)) {
        if (lk_qpush(&q, a, b) != 0) {
          free(isroot);
          goto fail;
        }
      }
  }
  free(isroot);

  cnt = calloc((size_t)L->nnd, sizeof *cnt);
  if (cnt == NULL) goto fail;

  for (int64_t k = 0; k < q.n; k++) {
    int32_t ia = q.p[k].a, ib = q.p[k].b;
    S->nvisit++;
    if (ia == ib) {
      /* САМОПАРА — обмен ВНУТРИ узла, то есть между его потомками. Раскрывается
       * во ВСЕ пары детей (включая их собственные самопары), а не по правилу
       * «дробить большего»: см. разбор в `plink.h`. Лист сам себя не видит. */
      int32_t s0 = K.start[ia], s1 = K.start[ia + 1];
      if (s1 <= s0) {
        S->nzero++;
        continue;
      }
      S->nrefine++;
      for (int32_t u = s0; u < s1; u++)
        for (int32_t v = s0; v < s1; v++)
          if (lk_qpush(&q, K.ch[u], K.ch[v]) != 0) goto fail;
      continue;
    }
    const hz_lodnode *A = &L->nd[ia], *B = &L->nd[ib];
    const double *xa = pt + 3 * ia, *xb = pt + 3 * ib;
    double f = lk_coef(A, xa, B, xb);
    if (!(f > 0.0)) { /* отвёрнуты друг от друга — связи нет по ориентации */
      S->nzero++;
      if (L->nd[ia].level > 0 || L->nd[ib].level > 0) S->nzero_coarse++;
      continue;
    }
    /* ДРОБИТЬ ИЛИ ЗАВОДИТЬ. Дробится тот, чей угловой размер больше: у пары
     * «стена и стул» дробить надо стену, а не стул. Если у него нет детей —
     * дробить нечем, связь заводится как есть. */
    double sa = lk_subtend(xb, A, xa), sb = lk_subtend(xa, B, xb);
    int32_t big = (sa > sb) ? ia : ib;
    double smax = (sa > sb) ? sa : sb;
    int32_t c0 = K.start[big], c1 = K.start[big + 1];
    if (smax > eps && c1 > c0) {
      S->nrefine++;
      for (int32_t c = c0; c < c1; c++) {
        int32_t ch = K.ch[c];
        int rc = (big == ia) ? lk_qpush(&q, ch, ib) : lk_qpush(&q, ia, ch);
        if (rc != 0) goto fail;
      }
      continue;
    }
    double v = lk_vis(sc, ia, xa, ib, xb, nvis, &S->nray);
    if (!(v > 0.0)) {
      S->nzero++;
      continue;
    }
    /* КОЭФФИЦИЕНТ СВЯЗИ. Дробление и отбраковка выше СОЗНАТЕЛЬНО остались на
     * диске (§87): набор связей обязан выйти тем же, и тогда сдвиг `Σf`
     * приписывается коэффициенту, и больше нечему. Точная форма ставится там,
     * где у излучателя ЕСТЬ край, то есть на нулевом уровне. */
    double fc = f;
    if (np2poly != NULL && np2poly[ib] >= 0) {
      fc = lk_coef_poly(A, xa, sc->ps, np2poly[ib], cfg->noclip, w1, w2);
      S->nlink_leaf++;
      S->wleaf += fc * v;
    }
    if (lk_push(S, ia, ib, fc * v) != 0) goto fail;
    cnt[ia]++;
  }
  for (int32_t k = 0; k < L->nnd; k++)
    if (cnt[k] > S->nmax_node) S->nmax_node = cnt[k];
  {
    double wtot = 0.0;
    for (int64_t k = 0; k < S->n; k++)
      wtot += S->l[k].f;
    S->wleaf = (wtot > 0.0) ? S->wleaf / wtot : 0.0;
  }
  free(cnt);
  free(q.p);
  free(np2poly);
  free(w1);
  free(w2);
  free(pt);
  lk_kids_free(&K);
  lk_closure(S, L);
  return 0;

fail:
  free(cnt);
  free(q.p);
  free(np2poly);
  free(w1);
  free(w2);
  free(pt);
  lk_kids_free(&K);
  return 2;
}

void hz_links_free(hz_linkset *S) {
  free(S->l);
  memset(S, 0, sizeof *S);
}

/* --- решение ---------------------------------------------------------------- */

int hz_links_solve(const hz_linkset *S, const hz_lod *L, const double *Le, const double *rho,
                   double *B, double tol, int maxit, double *reshist, int nopull) {
  const int32_t nn = L->nnd;
  double *acc = malloc((size_t)nn * sizeof *acc);
  double *ar = malloc((size_t)nn * sizeof *ar);
  double *Bp = malloc((size_t)nn * sizeof *Bp); /* снимок начала итерации */
  double *sw = calloc((size_t)nn, sizeof *sw);
  double *sa = calloc((size_t)nn, sizeof *sa);
  if (acc == NULL || ar == NULL || Bp == NULL || sw == NULL || sa == NULL) {
    free(acc);
    free(ar);
    free(Bp);
    free(sw);
    free(sa);
    return -1;
  }
  /* Площадь узла для взвешивания при подъёме. У узла она уже посчитана. */
  for (int32_t k = 0; k < nn; k++) {
    ar[k] = L->nd[k].area_surf;
    B[k] = Le[k];
  }
  int it = 0;
  for (; it < maxit; it++) {
    memcpy(Bp, B, (size_t)nn * sizeof *Bp);
    for (int32_t k = 0; k < nn; k++)
      acc[k] = 0.0;
    /* 1. СБОР: одна итерация — один порядок отражения (ряд Неймана). */
    for (int64_t l = 0; l < S->n; l++) {
      const hz_link *K = &S->l[l];
      acc[K->i] += K->f * B[K->j];
    }
    /* 2. ВНИЗ: облучённость предка наследуется потомками. Обход от грубых уровней
     * к мелким, поэтому вклад проходит всю цепочку за один проход. */
    for (int lev = L->nlev - 1; lev >= 1; lev--)
      for (int32_t k = 0; k < nn; k++) {
        if (L->nd[k].level != lev - 1) continue;
        int32_t p = L->nd[k].parent;
        if (p >= 0 && p < nn) acc[k] += acc[p];
      }
    /* 3. ЗАТЕНЕНИЕ. `acc` есть Σ форм-факторов на радианс, то есть облучённость,
     * делённая на `π`; радианс отражённого потому `ρ·acc`, БЕЗ второго деления на
     * `π`. Прежняя редакция делила ещё раз и понижала альбедо с `0.5` до `0.159`
     * (А195). */
    for (int32_t k = 0; k < nn; k++)
      B[k] = Le[k] + rho[k] * acc[k];
    /* 4. ВВЕРХ: радианс предка — средневзвешенное по площади радиансов детей.
     * Обход от мелких к грубым, чтобы дети были готовы раньше родителя. */
    if (!nopull)
      for (int lev = 1; lev < L->nlev; lev++) {
        memset(sw, 0, (size_t)nn * sizeof *sw);
        memset(sa, 0, (size_t)nn * sizeof *sa);
        for (int32_t k = 0; k < nn; k++) {
          if (L->nd[k].level != lev - 1) continue;
          int32_t p = L->nd[k].parent;
          if (p >= 0 && p < nn) {
            sw[p] += B[k] * ar[k];
            sa[p] += ar[k];
          }
        }
        for (int32_t k = 0; k < nn; k++)
          if (L->nd[k].level == lev && sa[k] > 0.0) B[k] = sw[k] / sa[k];
      }
    /* НЕВЯЗКА — максимум относительного изменения за ПОЛНУЮ итерацию. Ноль на
     * неподвижной точке по построению; разбор и оговорка — в `plink.h`. */
    double dmx = 0.0, bmx = 0.0;
    for (int32_t k = 0; k < nn; k++) {
      double d = fabs(B[k] - Bp[k]), b = fabs(B[k]);
      if (d > dmx) dmx = d;
      if (b > bmx) bmx = b;
    }
    double res = (bmx > 0.0) ? dmx / bmx : 0.0;
    if (reshist != NULL) reshist[it] = res;
    if (res < tol) {
      it++;
      break;
    }
  }
  free(acc);
  free(ar);
  free(Bp);
  free(sw);
  free(sa);
  return it;
}
