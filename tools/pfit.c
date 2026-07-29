/* pfit — ЗАПАС ЧЕБЫШЕВСКОЙ (МИНИМАКСНОЙ) ПЛОСКОСТИ (PLAN_ELEMENTS.md, §23).
 *
 * ВОПРОС РОВНО ОДИН: насколько нынешняя плоскость группы — средневзвешенная по
 * площади — далека от минимаксной, и сколько слияний это стоит. О11 обещает
 * `5…20%` прироста, и обещание не подпёрто НИЧЕМ: ни одного числа о расстоянии
 * до оптимума в проекте нет. Писать линейную программу в цикл слияния под
 * непроверенное обещание — ровно то, что §5 называет «экстраполяцией без эталона
 * стоимости».
 *
 * ОГРУБИТЕЛЬ ЭТОТ СТЕНД НЕ ТРОГАЕТ. Он берёт РЕАЛЬНЫЕ группы, выданные
 * `hz_merge`, и считает по ним две величины; критерий остаётся прежним, слепок
 * разметки обязан не измениться (проверяется `pcoarse`).
 *
 * МЕРЯЮТСЯ ДВА МНОЖЕСТВА, И ВТОРОЕ ВАЖНЕЕ (А79). На ПРИНЯТЫХ группах виден
 * запас — насколько плоскость могла бы быть лучше там, где слияние уже
 * состоялось. Но прирост слияний берётся из объединений, которые критерий
 * ОТВЕРГ: у них `d_avg ≥ δ`, и вопрос в том, у скольких `d_cheb < δ`. Поэтому
 * вторым проходом перебираются пары СОСЕДНИХ групп, и считается доля пар,
 * которые минимакс пропустил бы, а средневзвешенное не пропускает.
 *
 * ВЕЛИЧИНА — ВЕРХНЯЯ ОЦЕНКА ОПТИМУМА, и это сказано заранее: любая допустимая
 * плоскость даёт оценку сверху, поэтому посчитанный запас ЗАНИЖЕН. Ошибка идёт
 * в безопасную сторону — «нашли не меньше, чем есть».
 */

#include "pedge.h"
#include "pmerge.h"
#include "poly_seg.h"
#include "polygon.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int cmp_dbl(const void *x, const void *y) {
  double a = *(const double *)x, b = *(const double *)y;
  return (a < b) ? -1 : ((a > b) ? 1 : 0);
}

/* Точки в раме группы: `(u, v)` — в плоскости, `w` — отклонение от неё. */
typedef struct {
  double u, v, w;
} fp_pt;

/* ПОЛУШИРИНА ПОЛОСЫ ПРИ ЗАДАННОМ НАКЛОНЕ. Функция `(a, b) ↦ (max − min)/2` от
 * `w − a·u − b·v` выпуклая и КУСОЧНО-ЛИНЕЙНАЯ: максимум и минимум аффинных
 * функций. Её минимум по `(a, b)` и есть чебышевское отклонение в графовой
 * форме. */
static double fp_half(const fp_pt *p, int32_t n, double a, double b, double *mid) {
  double lo = 1e300, hi = -1e300;
  for (int32_t i = 0; i < n; i++) {
    double r = p[i].w - a * p[i].u - b * p[i].v;
    if (r < lo) lo = r;
    if (r > hi) hi = r;
  }
  if (mid != NULL) *mid = 0.5 * (hi + lo);
  return 0.5 * (hi - lo);
}

/* ОДНОМЕРНЫЙ МИНИМУМ ПО НАПРАВЛЕНИЮ `(da, db)` ТЕРНАРНЫМ ПОИСКОМ.
 * ПЛАТО ОБЯЗАНО СЖИМАТЬ ОТРЕЗОК К СЕРЕДИНЕ, А НЕ ОТБРАСЫВАТЬ СТОРОНУ (А81): у
 * кусочно-линейной выпуклой функции есть плоские участки, и при `f(m1) = f(m2)`
 * минимум лежит В `[m1, m2]` — наивное «равно, значит идём влево» его теряет. */
static void fp_line(const fp_pt *p, int32_t n, double *a, double *b, double da, double db,
                    double span) {
  double lo = -span, hi = span;
  for (int it = 0; it < 60; it++) {
    double m1 = lo + (hi - lo) / 3.0, m2 = hi - (hi - lo) / 3.0;
    double f1 = fp_half(p, n, *a + m1 * da, *b + m1 * db, NULL);
    double f2 = fp_half(p, n, *a + m2 * da, *b + m2 * db, NULL);
    if (f1 < f2)
      hi = m2;
    else if (f2 < f1)
      lo = m1;
    else {
      lo = m1;
      hi = m2;
    }
  }
  double t = 0.5 * (lo + hi);
  if (fp_half(p, n, *a + t * da, *b + t * db, NULL) < fp_half(p, n, *a, *b, NULL)) {
    *a += t * da;
    *b += t * db;
  }
}

/* ВЕРХНЯЯ ОЦЕНКА ЧЕБЫШЕВСКОГО ОТКЛОНЕНИЯ. Спуск покоординатный ПЛЮС по
 * диагоналям, и диагонали здесь не украшение (А82): у негладкой выпуклой функции
 * в изломе улучшения по осям может не быть при наличии улучшения по диагонали, и
 * тогда оценка выходит рыхлой, а запас — заниженным. Насколько этого хватает,
 * проверяется ПЕРЕБОРОМ на малых группах, а не рассуждением.
 * `*orth` — то же отклонение, приведённое к ОРТОГОНАЛЬНОМУ: критерий меряет
 * расстояние до плоскости с единичной нормалью, а графовая форма — вдоль оси. */
static double fp_cheb(const fp_pt *p, int32_t n, int sweeps, double span, double *orth, double *pa,
                      double *pb, double *pc) {
  double a = 0.0, b = 0.0;
  /* НАПРАВЛЕНИЙ ВОСЕМЬ, А НЕ ЧЕТЫРЕ, И ЭТО ПО ЗАМЕРУ (А87). С четырьмя (оси плюс
   * диагонали) спуск расходился с перебором на малых группах до `12.7%`: у
   * кусочно-линейной выпуклой функции излом может не пускать ни по одной из
   * четырёх, пропуская восьмую. Каждый проход идёт по всем восьми, шаг поиска
   * сужается вдвое от прохода к проходу. */
  for (int s = 0; s < sweeps; s++) {
    double sp = span / (double)(1 << s);
    if (!(sp > 0.0)) sp = span;
    for (int k = 0; k < 8; k++) {
      double th = 3.14159265358979323846 * (double)k / 8.0;
      fp_line(p, n, &a, &b, cos(th), sin(th), sp);
    }
  }
  double mid = 0.0;
  double h = fp_half(p, n, a, b, &mid);
  /* ОРТОГОНАЛЬНОЕ ОТКЛОНЕНИЕ, А НЕ ГРАФОВОЕ: критерий меряет расстояние до
   * плоскости с ЕДИНИЧНОЙ нормалью, графовая же форма меряет вдоль оси рамы.
   * Отношение — `1/sqrt(1+a²+b²)`. Сравнивать эти две величины между собой
   * нельзя: первый прогон стенда сверял ортогональную оценку с ГРАФОВЫМ
   * перебором и получил «расхождение 13%», которого нет (А84). */
  if (orth != NULL) *orth = h / sqrt(1.0 + a * a + b * b);
  if (pa != NULL) *pa = a;
  if (pb != NULL) *pb = b;
  if (pc != NULL) *pc = mid;
  return h;
}

/* ЭТАЛОН — ПЕРЕБОР, как `brute` у сетки кандидатов. Оптимум линейной программы
 * достигается в вершине, а вершина задаётся четырьмя активными ограничениями,
 * поэтому оптимальная плоскость проходит «по четырём точкам» — перебираем все
 * четвёрки и все знаковые комбинации, решаем 4×4 и берём лучший максимум по ВСЕМ
 * точкам. Годится только на малых группах, потому и стоит эталоном, а не
 * методом. Возвращает −1, если ни одна система не решилась. */
static int fp_solve4(double A[4][4], double rhs[4], double x[4]) {
  for (int c = 0; c < 4; c++) {
    int piv = c;
    for (int r = c + 1; r < 4; r++)
      if (fabs(A[r][c]) > fabs(A[piv][c])) piv = r;
    if (!(fabs(A[piv][c]) > 0.0)) return 1;
    if (piv != c) {
      for (int k = 0; k < 4; k++) {
        double t = A[c][k];
        A[c][k] = A[piv][k];
        A[piv][k] = t;
      }
      double t = rhs[c];
      rhs[c] = rhs[piv];
      rhs[piv] = t;
    }
    for (int r = 0; r < 4; r++) {
      if (r == c) continue;
      double f = A[r][c] / A[c][c];
      for (int k = 0; k < 4; k++)
        A[r][k] -= f * A[c][k];
      rhs[r] -= f * rhs[c];
    }
  }
  for (int c = 0; c < 4; c++)
    x[c] = rhs[c] / A[c][c];
  return 0;
}

static double fp_brute(const fp_pt *p, int32_t n) {
  double best = 1e300;
  for (int32_t i = 0; i < n; i++)
    for (int32_t j = i + 1; j < n; j++)
      for (int32_t k = j + 1; k < n; k++)
        for (int32_t l = k + 1; l < n; l++) {
          int32_t id[4] = {i, j, k, l};
          for (int sg = 0; sg < 8; sg++) {
            double A[4][4], rhs[4], x[4];
            for (int r = 0; r < 4; r++) {
              A[r][0] = p[id[r]].u;
              A[r][1] = p[id[r]].v;
              A[r][2] = 1.0;
              A[r][3] = (r == 0) ? 1.0 : (((sg >> (r - 1)) & 1) ? 1.0 : -1.0);
              rhs[r] = p[id[r]].w;
            }
            if (fp_solve4(A, rhs, x) != 0) continue;
            double m = 0.0;
            for (int32_t q = 0; q < n; q++) {
              double d = fabs(p[q].w - x[0] * p[q].u - x[1] * p[q].v - x[2]);
              if (d > m) m = d;
            }
            if (m < best) best = m;
          }
        }
  return (best < 1e299) ? best : -1.0;
}

/* Набор точек группы (или объединения групп) в раме её плоскости. */
static int32_t fp_gather(fp_pt *buf, int32_t cap, const hz_polyset *ps, const int32_t *mem,
                         int32_t nm, const double n[3], double off, const double eu[3],
                         const double ev[3], const double org[3]) {
  int32_t np = 0;
  for (int32_t i = 0; i < nm; i++) {
    const hz_poly *P = &ps->p[mem[i]];
    const double *S = ps->sup + (size_t)P->s0 * 3;
    for (int32_t q = 0; q < P->nsup; q++) {
      if (np >= cap) return np;
      double d[3] = {S[(size_t)q * 3] - org[0], S[(size_t)q * 3 + 1] - org[1],
                     S[(size_t)q * 3 + 2] - org[2]};
      buf[np].u = d[0] * eu[0] + d[1] * eu[1] + d[2] * eu[2];
      buf[np].v = d[0] * ev[0] + d[1] * ev[1] + d[2] * ev[2];
      buf[np].w =
          S[(size_t)q * 3] * n[0] + S[(size_t)q * 3 + 1] * n[1] + S[(size_t)q * 3 + 2] * n[2] - off;
      np++;
    }
  }
  return np;
}

static void fp_frame(const double n[3], double eu[3], double ev[3]) {
  int ax = 0;
  for (int a = 1; a < 3; a++)
    if (fabs(n[a]) < fabs(n[ax])) ax = a;
  double e0[3] = {0.0, 0.0, 0.0};
  e0[ax] = 1.0;
  double d = e0[0] * n[0] + e0[1] * n[1] + e0[2] * n[2];
  for (int a = 0; a < 3; a++)
    eu[a] = e0[a] - d * n[a];
  double en = sqrt(eu[0] * eu[0] + eu[1] * eu[1] + eu[2] * eu[2]);
  for (int a = 0; a < 3; a++)
    eu[a] /= en;
  ev[0] = n[1] * eu[2] - n[2] * eu[1];
  ev[1] = n[2] * eu[0] - n[0] * eu[2];
  ev[2] = n[0] * eu[1] - n[1] * eu[0];
}

static uint64_t fp_rnd(uint64_t *s) {
  *s ^= *s << 13;
  *s ^= *s >> 7;
  *s ^= *s << 17;
  return *s;
}

int main(int argc, char **argv) {
  int city = (argc > 1 && strcmp(argv[1], "city") == 0);
  double dseg = (argc > 2) ? strtod(argv[2], NULL) : (city ? 0.05 : 0.045);
  double dcoarse = city ? 2.0 : 0.3;
  int32_t target = city ? 20000 : 250;
  int simp = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "simp") == 0) simp = 1;
    if (strncmp(argv[i], "d=", 2) == 0) dcoarse = strtod(argv[i] + 2, NULL);
  }
  hz_objmesh m;
  if (hz_obj_load(&m, city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ,
                  city ? HZ_CFG_CITY_SCALE : HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  hz_pseglist sg;
  if (hz_seg_planar(&sg, &m, dseg) != 0) return 1;
  hz_polyset ps;
  if (hz_poly_build(&ps, &m, &sg) != 0) return 1;
  if (simp) {
    hz_edgestat es;
    if (hz_edge_simplify(&ps, 0.25 * dseg, HZ_EDGE_SHARED, &es) != 0) return 1;
  }
  hz_mergecfg mc;
  memset(&mc, 0, sizeof mc);
  mc.delta = dcoarse;
  mc.target = target;
  mc.use_geom = 1;
  mc.use_overlap = 1;
  hz_pseglist so;
  hz_mergestat st;
  double t0 = now_s();
  if (hz_merge(&so, &m, &sg, &ps, &mc, &st) != 0) return 1;
  printf("== запас минимаксной плоскости: %s, δ сегм %g, δ огр %g, групп %d (слепок %016llx)\n",
         city ? "ГОРОД" : "зал", dseg, dcoarse, so.nseg, (unsigned long long)st.digest);

  /* полигон -> группа и списки членов по группам */
  int32_t np = (ps.np < sg.nseg) ? ps.np : sg.nseg;
  int32_t *gof = malloc((size_t)np * sizeof *gof);
  int32_t *goff = calloc((size_t)so.nseg + 1, sizeof *goff);
  if (gof == NULL || goff == NULL) return 1;
  for (int32_t k = 0; k < np; k++)
    gof[k] = -1;
  for (int32_t t = 0; t < m.nt; t++) {
    int32_t k = sg.label[t];
    if (k < np && gof[k] < 0) gof[k] = so.label[t];
  }
  for (int32_t k = 0; k < np; k++)
    if (gof[k] >= 0) goff[gof[k] + 1]++;
  for (int32_t g = 0; g < so.nseg; g++)
    goff[g + 1] += goff[g];
  int32_t *gmem = malloc((size_t)np * sizeof *gmem);
  int32_t *cur = malloc((size_t)so.nseg * sizeof *cur);
  if (gmem == NULL || cur == NULL) return 1;
  memcpy(cur, goff, (size_t)so.nseg * sizeof *cur);
  for (int32_t k = 0; k < np; k++)
    if (gof[k] >= 0) gmem[cur[gof[k]]++] = k;
  free(cur);
  /* --- А122: СКОЛЬКО ГРУПП СХОДИТСЯ В ВЕРШИНЕ ---
   * Замер ставится ПЕРВЫМ, до всякой подгонки вершин, потому что от него зависит,
   * стоит ли писать О17 вообще: двигать имеет смысл только вершины, где сходятся
   * ТРИ плоскости и больше (А121 — при одной-двух задача недоопределена в двух
   * направлениях, и вершина поедет вдоль поверхности). Если таких вершин единицы
   * процентов, вся ставка шага беспочвенна.
   * Считается по СВАРНЫМ номерам (`ps->bw`, заведены для О1): один и тот же
   * объект у всех полигонов, которые эту вершину делят. Группы обходятся ПОДРЯД
   * (по `gmem`), поэтому «последняя виденная группа» различает их точно. */
  {
    int32_t nw = (int32_t)ps.nweld;
    int32_t *lastg = malloc((size_t)(nw > 0 ? nw : 1) * sizeof *lastg);
    int32_t *ngrp = calloc((size_t)(nw > 0 ? nw : 1), sizeof *ngrp);
    if (lastg != NULL && ngrp != NULL) {
      for (int32_t i = 0; i < nw; i++)
        lastg[i] = -1;
      for (int32_t g = 0; g < so.nseg; g++)
        for (int32_t i = goff[g]; i < goff[g + 1]; i++) {
          const hz_poly *P = &ps.p[gmem[i]];
          for (int32_t l = P->l0; l < P->l0 + P->nloop; l++)
            for (int32_t b = ps.loop[l]; b < ps.loop[l + 1]; b++) {
              int32_t w = ps.bw[b];
              if (w < 0 || w >= nw) continue;
              if (lastg[w] != g) {
                lastg[w] = g;
                ngrp[w]++;
              }
            }
        }
      int64_t hist[9];
      memset(hist, 0, sizeof hist);
      int64_t tot = 0, ge3 = 0;
      for (int32_t i = 0; i < nw; i++) {
        if (ngrp[i] <= 0) continue;
        tot++;
        int k = ngrp[i] > 8 ? 8 : ngrp[i];
        hist[k]++;
        if (ngrp[i] >= 3) ge3++;
      }
      printf("   А122: групп в вершине (по сварным номерам, %lld вершин на краю):", (long long)tot);
      for (int k = 1; k <= 8; k++)
        if (hist[k] > 0)
          printf(" [%d%s] %.2f%%", k, (k == 8) ? "+" : "",
                 100.0 * (double)hist[k] / (double)(tot ? tot : 1));
      printf("; ТРИ И БОЛЬШЕ — %lld (%.2f%%)\n", (long long)ge3,
             100.0 * (double)ge3 / (double)(tot ? tot : 1));
    }
    free(lastg);
    free(ngrp);
  }

  int32_t cap = 1 << 20;
  fp_pt *buf = malloc((size_t)cap * sizeof *buf);
  double *rat = malloc((size_t)so.nseg * sizeof *rat);
  double *res = malloc((size_t)so.nseg * sizeof *res);
  if (buf == NULL || rat == NULL || res == NULL) return 1;

  /* --- ПРОХОД 1: запас на ПРИНЯТЫХ группах --- */
  double t1 = now_s();
  int32_t nrat = 0, nbord = 0, nbrute = 0;
  double wdiff_avg = 0.0, wbrute = 0.0, wbrute_lo = 0.0, qsum = 0.0;
  int32_t nq = 0;
  uint64_t rs = 0x243F6A8885A308D3ULL;
  for (int32_t g = 0; g < so.nseg; g++) {
    int32_t nm = goff[g + 1] - goff[g];
    if (nm <= 0) continue;
    double eu[3], ev[3], org[3];
    fp_frame(so.seg[g].n, eu, ev);
    for (int a = 0; a < 3; a++)
      org[a] = so.seg[g].off * so.seg[g].n[a];
    int32_t n =
        fp_gather(buf, cap, &ps, gmem + goff[g], nm, so.seg[g].n, so.seg[g].off, eu, ev, org);
    if (n < 4) continue;
    double davg = 0.0, ext = 0.0;
    for (int32_t i = 0; i < n; i++) {
      double d = fabs(buf[i].w);
      if (d > davg) davg = d;
      double e = fabs(buf[i].u) > fabs(buf[i].v) ? fabs(buf[i].u) : fabs(buf[i].v);
      if (e > ext) ext = e;
    }
    /* А80: `d_avg` по опорным точкам обязан совпасть с докладываемым `dmax`,
     * который считается по ТРЕУГОЛЬНИКАМ. После О7 это одна величина. */
    double dd = fabs(davg - so.seg[g].dmax);
    if (dd > wdiff_avg) wdiff_avg = dd;
    double span = (ext > 0.0) ? 8.0 * (davg + 1e-12) / ext : 1e-6;
    double orth = 0.0, fa = 0.0, fb = 0.0, fc = 0.0;
    double hgraph = fp_cheb(buf, n, 12, span, &orth, &fa, &fb, &fc);
    if (davg > 0.0)
      rat[nrat] = orth / davg;
    else
      rat[nrat] = 1.0;
    res[nrat] = (davg - orth) / dcoarse;
    nrat++;
    if (so.seg[g].dmax > 0.9 * dcoarse) nbord++;
    /* ЭТАЛОН на малых группах */
    if (n <= 12) {
      double bb = fp_brute(buf, n);
      if (bb >= 0.0) {
        /* СРАВНИВАТЬ ОДНОРОДНОЕ С ОДНОРОДНЫМ (А84). Перебор даёт ГРАФОВЫЙ
         * минимакс — расстояние вдоль оси рамы, — поэтому и спуск сверяется
         * своим графовым значением, а не приведённым к ортогональному. Первая
         * редакция сверяла ортогональное с графовым и напечатала «расхождение
         * 13%», которого не существует: это множитель `1/√(1+a²+b²)`. */
        double rel = (hgraph - bb) / (bb + 1e-30);
        if (rel > wbrute) wbrute = rel;
        if (rel < wbrute_lo) wbrute_lo = rel;
        nbrute++;
      }
    }
    /* НЕГАТИВНЫЙ КОНТРОЛЬ (А83): подгонка по СЛУЧАЙНОЙ ЧЕТВЕРТИ точек, максимум
     * считается по ВСЕМ. Обязана быть заметно хуже полного минимакса; отношение
     * около единицы означало бы, что крайние точки в подгонке не участвуют, и
     * вся посылка неверна. */
    if (n >= 16 && nq < 20000) {
      int32_t nqp = 0;
      for (int32_t i = 0; i < n; i++)
        if ((fp_rnd(&rs) & 3u) == 0u) {
          buf[cap - 1 - nqp] = buf[i];
          nqp++;
          if (nqp > n / 4 + 4) break;
        }
      if (nqp >= 4) {
        /* ПЛОСКОСТЬ, ПОДОГНАННАЯ ПО ЧЕТВЕРТИ, ПРИМЕНЯЕТСЯ КО ВСЕМ ТОЧКАМ.
         * Первая редакция считала здесь максимум по всем точкам от ИСХОДНОЙ
         * плоскости, то есть печатала `d_avg/d_cheb` под именем контроля
         * (А85). Контроль обязан ставить СВОЮ плоскость. */
        double qa = 0.0, qb = 0.0, qc = 0.0, o2 = 0.0;
        fp_cheb(buf + cap - nqp, nqp, 6, span, &o2, &qa, &qb, &qc);
        double full = 0.0;
        for (int32_t i = 0; i < n; i++) {
          double d = fabs(buf[i].w - qa * buf[i].u - qb * buf[i].v - qc);
          if (d > full) full = d;
        }
        full /= sqrt(1.0 + qa * qa + qb * qb);
        if (orth > 0.0) {
          qsum += full / orth;
          nq++;
        }
      }
    }
  }
  double t2 = now_s();
  qsort(rat, (size_t)nrat, sizeof *rat, cmp_dbl);
  qsort(res, (size_t)nrat, sizeof *res, cmp_dbl);
  printf("   ПРИНЯТЫЕ ГРУППЫ (%d): d_cheb/d_avg p50 %.4f, p90 %.4f, p99 %.4f, минимум %.4f\n", nrat,
         rat[nrat / 2], rat[nrat - nrat / 10 - 1], rat[nrat - nrat / 100 - 1], rat[0]);
  printf("   запас (d_avg − d_cheb)/δ: p50 %.4f, p90 %.4f, максимум %.4f; групп у границы "
         "(dmax > 0.9δ) %d (%.1f%%)\n",
         res[nrat / 2], res[nrat - nrat / 10 - 1], res[nrat - 1], nbord,
         100.0 * (double)nbord / (double)(nrat > 0 ? nrat : 1));
  printf("   СВЕРКА: |d_avg по опорным − dmax по треугольникам| максимум %.3e; эталон перебором "
         "на %d малых группах, расхождение спуска от %.3e до %.3e (спуск обязан быть НЕ НИЖЕ "
         "перебора: он даёт оценку сверху)\n",
         wdiff_avg, nbrute, wbrute_lo, wbrute);
  printf("   негативный контроль (подгонка по четверти): отношение к полному минимаксу в среднем "
         "%.3f по %d группам\n",
         (nq > 0) ? qsum / (double)nq : 0.0, nq);
  printf("   время: слияние %.2f с, проход 1 %.2f с\n", t1 - t0, t2 - t1);

  /* --- ПРОХОД 2: ОТВЕРГНУТЫЕ объединения (А79) --- */
  /* Пары СОСЕДНИХ групп: те, что делят ячейку сетки по коробкам групп. Именно на
   * них и стоит вопрос О11 — сколько объединений минимакс пропустил бы. */
  double *gb = malloc((size_t)so.nseg * 6 * sizeof *gb);
  if (gb == NULL) return 1;
  for (int32_t g = 0; g < so.nseg; g++) {
    for (int c = 0; c < 3; c++) {
      gb[(size_t)g * 6 + (size_t)c] = 1e300;
      gb[(size_t)g * 6 + 3 + (size_t)c] = -1e300;
    }
    for (int32_t i = goff[g]; i < goff[g + 1]; i++) {
      const hz_poly *P = &ps.p[gmem[i]];
      const double *S = ps.sup + (size_t)P->s0 * 3;
      for (int32_t q = 0; q < P->nsup; q++)
        for (int c = 0; c < 3; c++) {
          double x = S[(size_t)q * 3 + (size_t)c];
          if (x < gb[(size_t)g * 6 + (size_t)c]) gb[(size_t)g * 6 + (size_t)c] = x;
          if (x > gb[(size_t)g * 6 + 3 + (size_t)c]) gb[(size_t)g * 6 + 3 + (size_t)c] = x;
        }
    }
  }
  int64_t npair = 0, nwould = 0, nboth = 0;
  uint64_t rs2 = 0x13198A2E03707344ULL;
  int32_t *mem2 = malloc((size_t)np * sizeof *mem2);
  if (mem2 == NULL) return 1;
  double t3 = now_s();
  for (int32_t g = 0; g < so.nseg && npair < 20000; g++) {
    for (int32_t h = g + 1; h < so.nseg && npair < 20000; h++) {
      int near = 1;
      for (int c = 0; c < 3 && near; c++) {
        if (gb[(size_t)g * 6 + (size_t)c] - gb[(size_t)h * 6 + 3 + (size_t)c] > dcoarse) near = 0;
        if (gb[(size_t)h * 6 + (size_t)c] - gb[(size_t)g * 6 + 3 + (size_t)c] > dcoarse) near = 0;
      }
      if (!near) continue;
      /* ПРОРЕЖИВАНИЕ ТОЛЬКО ТАМ, ГДЕ ПАР МНОГО. На зале групп 250, соседних пар
       * десятки, и прореживание вчетверо-шестьдесятчетверо оставляло выборку в
       * 17 пар — числа из такой выборки не значат ничего (А86). */
      if (so.nseg > 5000 && (fp_rnd(&rs2) & 63u) != 0u) continue;
      int32_t nm = 0;
      for (int32_t i = goff[g]; i < goff[g + 1]; i++)
        mem2[nm++] = gmem[i];
      for (int32_t i = goff[h]; i < goff[h + 1]; i++)
        mem2[nm++] = gmem[i];
      /* плоскость объединения — той же формулой, что в критерии */
      double sA = 0.0, sn[3] = {0, 0, 0}, sorg[3] = {0, 0, 0};
      for (int32_t i = 0; i < nm; i++) {
        const hz_poly *P = &ps.p[mem2[i]];
        sA += P->area;
        for (int c = 0; c < 3; c++) {
          sn[c] += P->area * P->n[c];
          sorg[c] += P->area * P->org[c];
        }
      }
      double nl = sqrt(sn[0] * sn[0] + sn[1] * sn[1] + sn[2] * sn[2]);
      if (!(sA > 0.0) || !(nl > 0.0)) continue;
      double nu[3], org[3], off = 0.0;
      int ok = 1;
      for (int c = 0; c < 3; c++) {
        nu[c] = sn[c] / nl;
        org[c] = sorg[c] / sA;
      }
      for (int c = 0; c < 3; c++)
        off += nu[c] * org[c];
      for (int32_t i = 0; i < nm && ok; i++) {
        const hz_poly *P = &ps.p[mem2[i]];
        if (!(P->n[0] * nu[0] + P->n[1] * nu[1] + P->n[2] * nu[2] > 0.0)) ok = 0;
      }
      if (!ok) continue; /* встречные грани — это не про плоскость (А22) */
      double eu[3], ev[3];
      fp_frame(nu, eu, ev);
      int32_t n = fp_gather(buf, cap, &ps, mem2, nm, nu, off, eu, ev, org);
      if (n < 4) continue;
      double davg = 0.0, ext = 0.0;
      for (int32_t i = 0; i < n; i++) {
        double d = fabs(buf[i].w);
        if (d > davg) davg = d;
        double e = fabs(buf[i].u) > fabs(buf[i].v) ? fabs(buf[i].u) : fabs(buf[i].v);
        if (e > ext) ext = e;
      }
      double span = (ext > 0.0) ? 8.0 * (davg + 1e-12) / ext : 1e-6;
      double orth = 0.0;
      fp_cheb(buf, n, 6, span, &orth, NULL, NULL, NULL);
      npair++;
      if (davg >= dcoarse && orth < dcoarse) nwould++;
      if (davg < dcoarse) nboth++;
    }
  }
  printf("   ОТВЕРГНУТЫЕ ОБЪЕДИНЕНИЯ (выборка %lld пар соседних групп): минимакс пропустил бы "
         "%lld (%.2f%%); из выборки и так проходили %lld (%.2f%%)\n",
         (long long)npair, (long long)nwould, 100.0 * (double)nwould / (double)(npair ? npair : 1),
         (long long)nboth, 100.0 * (double)nboth / (double)(npair ? npair : 1));
  printf("   время: проход 2 %.2f с\n", now_s() - t3);

  free(mem2);
  free(gb);
  free(buf);
  free(rat);
  free(res);
  free(gof);
  free(goff);
  free(gmem);
  hz_seg_free(&so);
  hz_poly_free(&ps);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
