/* pfrontx.c — ШАГ О72 (Ф2, переделка): ОБХОД ФРОНТА. План §251, аудит §252.
 *
 * ЗАЧЕМ ПЕРЕДЕЛКА. О71 мерил ПУЧОК ЛУЧЕЙ, а фронт — не пучок (поправка
 * пользователя 08-05). §241.7: «нет ни активного списка, ни очереди… есть ОБХОД
 * ДЕРЕВА В ПОРЯДКЕ ОКТАНТОВ… ОДИН ВИЗИТ НА УЗЕЛ». Пучок навещает узел столько
 * раз, сколько лучей его задели (замерено: 6.7…8.6 посещения на шаг со спуском
 * от корня), и взвешивает дальние оболочки как 1/R², тогда как фронт метёт
 * оболочку целиком.
 *
 * ГЛАВНОЕ ПРАВИЛО СПУСКА — §241.4, И ОНО НЕ ТО ЖЕ, ЧТО В О70/О71: «пол
 * ограничивает снизу размер ЭЛЕМЕНТА, а НЕ размер пустого узла». Пустой узел
 * берётся ЦЕЛИКОМ, каков бы ни был его размер. Применение пола к пустоте есть
 * НЕГАТИВНЫЙ КОНТРОЛЬ (`voidfloor=1`), и он обязан взорваться: двор 70³ м при
 * h = 0.1 м даёт 3.4e11 ячеек.
 *
 * ЗАПУСК:
 *   pfrontx ФАЙЛ.obj МАСШТАБ [leaf=N] [lev=N] [grade=0|1] [src=x,y,z] [px=4]
 *           [voidfloor=1] [bboxocc=1]
 */
#include "pclip.h"
#include "ptree.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ОКТАВЫ РАССТОЯНИЯ ОТ ИСТОЧНИКА и ГИСТОГРАММА РАЗМЕРА ЯЧЕЙКИ — те же отсчёты,
 * что в О71, чтобы числа двух шагов лежали в одной шкале. */
#define OCT0 (-8)
#define NOCT 32
#define SZ0 (-24)
#define NSZ 48
/* ПРЕДЕЛ ЧИСЛА ЯЧЕЕК. Нужен ТОЛЬКО негативному контролю: при спуске по полу в
 * пустоте счёт уходит в 3.4e11 (§241.4), и без предела прогон не кончится
 * никогда. Полтораста миллионов — вчетверо больше самого большого ЗАКОННОГО
 * среза, измеренного в О70 (37.3 млн ячеек у Сан-Мигеля до листьев).
 * Срабатывание печатается: молчаливый обрыв читался бы как «столько и было». */
#define HZ_PFRONT_MAXCELL 150000000
/* ПРЕДЕЛ УРОВНЯ ДЛЯ ВИРТУАЛЬНОГО ДРОБЛЕНИЯ. Нужен потому, что у ячейки,
 * СОДЕРЖАЩЕЙ ИСТОЧНИК, угловой критерий не выполняется НИКОГДА (`R < rad` при
 * любом дроблении), и рекурсия не имеет естественного дна. В дереве дном служил
 * лист; у виртуального дробления его нет — первый прогон контроля это и показал
 * переполнением стека. Тридцать уровней: у зала это `8.12 м / 2³⁰ = 7.6` нм,
 * то есть заведомо ниже любого физического смысла. Срабатывание считается. */
#define HZ_PFRONT_MAXLEV 30

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int box_hits(const double *bl, const double *bh, const double *clo, const double *chi) {
  for (int c = 0; c < 3; c++)
    if (!(bl[c] < chi[c]) || !(bh[c] >= clo[c])) return 0;
  return 1;
}

/* ---- ЗАНЯТОСТЬ (определение О70: кусок с ненулевой площадью) --------------- */
typedef struct {
  const hz_objmesh *m;
  const hz_ptree *T;
  const double *tlo, *thi, *tpl;
  unsigned char *occ, *bocc;
  int fail;
} occctx;

static void occ_walk(occctx *X, int32_t nid, const int32_t *list, int32_t n) {
  if (X->fail) return;
  const hz_ptnode *N = &X->T->nd[nid];
  if (n > 0) X->bocc[nid] = 1; /* габаритная занятость — верхняя оценка (А511) */
  if (N->child < 0) {
    const double *lo = N->lo, *hi = N->hi;
    for (int32_t i = 0; i < n; i++) {
      int32_t t = list[i];
      const double *tp = X->tpl + 4 * (size_t)t;
      double ctr[3], hlf[3], sd = 0.0, rr = 0.0;
      for (int c = 0; c < 3; c++) {
        ctr[c] = 0.5 * (lo[c] + hi[c]);
        hlf[c] = 0.5 * (hi[c] - lo[c]);
        sd += tp[c] * ctr[c];
        rr += fabs(tp[c]) * hlf[c];
      }
      sd -= tp[3];
      if (fabs(sd) > rr) continue;
      const double *A = X->m->v + 3 * (size_t)X->m->f[3 * (size_t)t + 0];
      const double *B = X->m->v + 3 * (size_t)X->m->f[3 * (size_t)t + 1];
      const double *C = X->m->v + 3 * (size_t)X->m->f[3 * (size_t)t + 2];
      hz_pclip_poly P;
      hz_pclip_tri(A, B, C, lo, hi, &P);
      if (P.nv >= 3 && hz_pclip_area(&P) > 0.0) {
        X->occ[nid] = 1;
        return;
      }
    }
    return;
  }
  int32_t c0 = N->child;
  int32_t *sub = malloc((size_t)(n > 0 ? n : 1) * sizeof *sub);
  if (sub == NULL) {
    X->fail = 1;
    return;
  }
  for (int k = 0; k < 8; k++) {
    const double *clo = X->T->nd[c0 + k].lo, *chi = X->T->nd[c0 + k].hi;
    int32_t ns = 0;
    for (int32_t i = 0; i < n; i++) {
      const double *bl = X->tlo + 3 * (size_t)list[i], *bh = X->thi + 3 * (size_t)list[i];
      if (box_hits(bl, bh, clo, chi)) sub[ns++] = list[i];
    }
    occ_walk(X, c0 + k, sub, ns);
  }
  free(sub);
}

/* ---- ОБХОД ФРОНТА --------------------------------------------------------- */
typedef struct {
  const hz_ptree *T;
  const unsigned char *occ;
  unsigned char *seen; /* А509: «один визит на узел» ПРОВЕРЯЕТСЯ, а не заявляется */
  double src[3], eps, px;
  /* СОЛНЦЕ: направление и глубина отсчёта. Точечных источников не бывает
   * (пользователь 08-05), а у нашей же модели цены точечность есть сингулярность:
   * при alpha -> 0 величина Omega/alpha² уходит в бесконечность. У параллельного
   * источника пол задаётся УГЛОВЫМ РАЗМЕРОМ источника и пройденным путём, а не
   * расстоянием до точки: полутень растёт как alpha·(d - d_рожд), §241.11. На
   * уровне Ф2 за d_рожд берётся вход в коробку сцены — это верхняя оценка
   * размытия, потому что настоящее d_рожд лежит позже. */
  int sun;        /* 1 — параллельный фронт */
  double dir[3];  /* направление распространения, нормировано */
  double t_entry; /* проекция ближней грани корня на dir */
  double pxeps2;  /* (px·ε)² — пол в квадрате, чтобы в горячем пути не было корней */
  int lean;       /* 1 — обход БЕЗ диагностики: чистая стоимость индекса */
  int voidfloor;  /* негативный контроль: пол применяется и к пустоте */
  int64_t ncell, nempty, ngeo, nrevisit, ncap;
  int64_t nnode; /* узлов ПРОЙДЕНО, включая внутренние: цена в обращениях к памяти */
  int64_t stop_void, stop_floor, stop_leaf;
  int64_t oct_cell[NOCT], oct_empty[NOCT];
  int64_t oct_lev[NOCT]; /* сумма уровней ячеек в октаве: падает ли детальность с расстоянием */
  int64_t sib_empty[9]; /* у родителя с геометрией — сколько из 8 детей пусты (потенциал слияния) */
  int64_t sz[NSZ];
  int64_t ndeep;     /* виртуальных ячеек, упёршихся в предел уровня */
  int64_t nstraddle; /* ячеек, чьи ближняя и дальняя грани в разных октавах (А507) */
  int32_t levmin, levmax;
} frontstat;

static int oct_of(double d) {
  if (!(d > 0.0)) return 0;
  int o = (int)floor(log2(d)) - OCT0;
  if (o < 0) o = 0;
  if (o >= NOCT) o = NOCT - 1;
  return o;
}

static int sz_of(double h) {
  if (!(h > 0.0)) return 0;
  int i = (int)floor(log2(h)) - SZ0;
  if (i < 0) i = 0;
  if (i >= NSZ) i = NSZ - 1;
  return i;
}

/* Взять узел как ячейку среза. Разложение по ПРИЧИНЕ остановки — сразу (А510). */
static void take(frontstat *S, int32_t nid, int why) {
  const hz_ptnode *N = &S->T->nd[nid];
  S->ncell++;
  /* ЧИСТАЯ СТОИМОСТЬ ИНДЕКСА МЕРЯЕТСЯ ОТДЕЛЬНО ОТ СТОИМОСТИ ЕЁ ИЗМЕРЕНИЯ.
   * Всё, что ниже, — диагностика: три корня на расстояния и два `log2` на
   * гистограммы. В первом докладе О72 они попали в число «нс на ячейку», и
   * получилось, что обход медленнее, чем он есть (поправка пользователя 08-05).
   * `lean=1` выключает диагностику целиком; разность двух прогонов и есть цена
   * замера. */
  if (S->lean) {
    if (why == 0) S->stop_void++;
    if (why == 1) S->stop_floor++;
    if (why == 2) S->stop_leaf++;
    return;
  }
  if (S->seen[nid]) S->nrevisit++;
  S->seen[nid] = 1;
  int lev = (int)S->T->lev[nid];
  if (lev < S->levmin) S->levmin = lev;
  if (lev > S->levmax) S->levmax = lev;
  double dc = 0.0, dnear = 0.0, dfar = 0.0;
  for (int c = 0; c < 3; c++) {
    double ctr = 0.5 * (N->lo[c] + N->hi[c]) - S->src[c];
    dc += ctr * ctr;
    double a = N->lo[c] - S->src[c], b = N->hi[c] - S->src[c];
    double lo = (a > 0.0) ? a : ((b < 0.0) ? -b : 0.0);
    double hi = (fabs(a) > fabs(b)) ? fabs(a) : fabs(b);
    dnear += lo * lo;
    dfar += hi * hi;
  }
  int o = oct_of(sqrt(dc));
  S->oct_cell[o]++;
  S->oct_lev[o] += lev;
  if (oct_of(sqrt(dnear)) != oct_of(sqrt(dfar))) S->nstraddle++;
  S->sz[sz_of(N->hi[0] - N->lo[0])]++;
  if (S->occ[nid]) {
    S->ngeo++;
  } else {
    S->nempty++;
    S->oct_empty[o]++;
  }
  if (why == 0) S->stop_void++;
  if (why == 1) S->stop_floor++;
  if (why == 2) S->stop_leaf++;
}

/* ---- ОПЫТ: ТОТ ЖЕ ОБХОД ПО КОМПАКТНОМУ ИНДЕКСУ ---------------------------- *
 * ЗАЧЕМ. Замер показал 3.6 ГБ/с при 64 Б на узел — это латентность памяти, а не
 * счёт. Из 64 Б сорок восемь занимает КОРОБКА, а она ВЫВОДИТСЯ: коробка ребёнка
 * есть половина родительской, и деление пополам двоичных значений точно. Опыт
 * проверяет, действительно ли дело в размере узла: обход идёт по массиву одних
 * ИНДЕКСОВ ДЕТЕЙ (4 Б) плюс байт занятости, а коробка передаётся вниз по
 * рекурсии. Результат обязан СОВПАСТЬ по числу ячеек — иначе меряется не то. */
static void compact_walk(frontstat *S, const int32_t *ch, int32_t nid, const double *lo,
                         const double *hi) {
  S->nnode++;
  if (!S->occ[nid]) {
    S->ncell++;
    S->stop_void++;
    return;
  }
  if (ch[nid] < 0) {
    S->ncell++;
    S->stop_leaf++;
    return;
  }
  double r2 = 0.0, d2 = 0.0, mid[3];
  for (int c = 0; c < 3; c++) {
    double h = 0.5 * (hi[c] - lo[c]);
    mid[c] = 0.5 * (lo[c] + hi[c]);
    double dd = mid[c] - S->src[c];
    r2 += h * h;
    d2 += dd * dd;
  }
  if (d2 > r2 && 4.0 * r2 <= S->pxeps2 * d2) {
    S->ncell++;
    S->stop_floor++;
    return;
  }
  int near = 0;
  for (int c = 0; c < 3; c++)
    if (S->src[c] >= mid[c]) near |= 1 << c;
  for (int k = 0; k < 8; k++) {
    int q = k ^ near;
    double clo[3], chi[3];
    for (int c = 0; c < 3; c++) {
      clo[c] = (q & (1 << c)) ? mid[c] : lo[c];
      chi[c] = (q & (1 << c)) ? hi[c] : mid[c];
    }
    compact_walk(S, ch, ch[nid] + q, clo, chi);
  }
}

/* НЕГАТИВНЫЙ КОНТРОЛЬ, ЧЕСТНАЯ ФОРМА. Нарушение §241.4 — это дробление ПУСТОТЫ
 * по полу, а узлов там нет: дерево пустоту не дробило. Поэтому контроль ведёт
 * ВИРТУАЛЬНОЕ дробление коробки, без узлов. Он обязан упереться в предел ячеек:
 * §241.4 считает для двора 70³ м при h = 0.1 м — 3.4e11 ячеек. Если предел НЕ
 * сработает, значит правило ничего не стоит и его можно не соблюдать. */
static void void_walk(frontstat *S, const double *lo, const double *hi, int lev) {
  if (S->ncell >= HZ_PFRONT_MAXCELL) {
    S->ncap = 1;
    return;
  }
  double r2 = 0.0, d2 = 0.0;
  for (int c = 0; c < 3; c++) {
    double h = 0.5 * (hi[c] - lo[c]);
    double dd = 0.5 * (lo[c] + hi[c]) - S->src[c];
    r2 += h * h;
    d2 += dd * dd;
  }
  /* Та же запись без корней, что и в рабочем обходе: контроль обязан мерить ту
   * же цену, иначе он сравнивает не то. */
  if (lev >= HZ_PFRONT_MAXLEV || (d2 > r2 && 4.0 * r2 <= S->pxeps2 * d2)) {
    if (lev >= HZ_PFRONT_MAXLEV) S->ndeep++;
    S->ncell++;
    S->nempty++;
    S->stop_floor++;
    /* Корень остался ТОЛЬКО здесь — в гистограмме, то есть в диагностике, а не
     * в критерии спуска. */
    int o = oct_of(sqrt(d2));
    S->oct_cell[o]++;
    S->oct_empty[o]++;
    S->sz[sz_of(hi[0] - lo[0])]++;
    if (lev > S->levmax) S->levmax = lev;
    if (lev < S->levmin) S->levmin = lev;
    return;
  }
  double mid[3];
  for (int c = 0; c < 3; c++)
    mid[c] = 0.5 * (lo[c] + hi[c]);
  for (int k = 0; k < 8; k++) {
    double clo[3], chi[3];
    for (int c = 0; c < 3; c++) {
      clo[c] = (k & (1 << c)) ? mid[c] : lo[c];
      chi[c] = (k & (1 << c)) ? hi[c] : mid[c];
    }
    void_walk(S, clo, chi, lev + 1);
    if (S->ncap) return;
  }
}

static void front_walk(frontstat *S, int32_t nid) {
  if (S->ncell >= HZ_PFRONT_MAXCELL) {
    S->ncap = 1;
    return;
  }
  const hz_ptnode *N = &S->T->nd[nid];
  S->nnode++;
  /* §241.4: ПУСТОЙ УЗЕЛ БЕРЁТСЯ ЦЕЛИКОМ. */
  if (!S->occ[nid]) {
    if (S->voidfloor) {
      void_walk(S, N->lo, N->hi, (int)S->T->lev[nid]);
      return;
    }
    take(S, nid, 0);
    return;
  }
  if (N->child < 0) {
    take(S, nid, 2);
    return;
  }
  /* Пол `ε·d` ОТ ИСТОЧНИКА (А508: у Ф2 он от источника, а не от камеры).
   *
   * БЕЗ КОРНЕЙ И БЕЗ ДЕЛЕНИЯ. Условие `(2·rad/R)/ε ≤ px` равносильно
   * `4·rad² ≤ (px·ε)²·R²`, а `R > rad` — `d2 > r2`. Прежняя запись стоила два
   * `sqrt` и деление НА КАЖДЫЙ УЗЕЛ, и это была не стоимость индекса, а
   * стоимость её записи (поправка пользователя 08-05). */
  double r2 = 0.0, d2 = 0.0;
  if (S->sun) {
    /* ПАРАЛЛЕЛЬНЫЙ ФРОНТ. Расстояние — ПРОЙДЕННЫЙ ПУТЬ от входа в сцену вдоль
     * направления, а не удаление от точки. Радиус ячейки тот же. */
    double t = 0.0;
    for (int c = 0; c < 3; c++) {
      double h = 0.5 * (N->hi[c] - N->lo[c]);
      r2 += h * h;
      t += 0.5 * (N->lo[c] + N->hi[c]) * S->dir[c];
    }
    double dep = t - S->t_entry;
    d2 = dep * dep;
  } else {
    for (int c = 0; c < 3; c++) {
      double h = 0.5 * (N->hi[c] - N->lo[c]);
      double dd = 0.5 * (N->lo[c] + N->hi[c]) - S->src[c];
      r2 += h * h;
      d2 += dd * dd;
    }
  }
  if (d2 > r2 && 4.0 * r2 <= S->pxeps2 * d2) {
    take(S, nid, 1);
    return;
  }
  /* ПОРЯДОК ОКТАНТОВ ОТ ИСТОЧНИКА (§241.7): ближние дети раньше дальних. На
   * СЧЁТ он не влияет, но обход у Ф3 будет этот же, и заводить его надо здесь,
   * а не переписывать потом. */
  int32_t c0 = N->child;
  if (!S->lean) {
    int ne = 0;
    for (int k = 0; k < 8; k++)
      if (!S->occ[c0 + k]) ne++;
    S->sib_empty[ne]++;
  }
  int near = 0;
  for (int c = 0; c < 3; c++)
    if (S->src[c] >= 0.5 * (N->lo[c] + N->hi[c])) near |= 1 << c;
  for (int k = 0; k < 8; k++)
    front_walk(S, c0 + (k ^ near));
}

static void report(const frontstat *S, double secs, const char *name) {
  printf("   %s: ЯЧЕЕК %lld (пустых %lld, с геометрией %lld), уровни %d…%d\n", name,
         (long long)S->ncell, (long long)S->nempty, (long long)S->ngeo, S->levmin, S->levmax);
  printf("      УЗЛОВ ПРОЙДЕНО %lld (внутренних %lld); байт узлов %.2f ГБ\n", (long long)S->nnode,
         (long long)(S->nnode - S->ncell),
         (double)S->nnode * (double)sizeof(hz_ptnode) / 1073741824.0);
  printf("      ПОВТОРНЫХ ВИЗИТОВ %lld (обязано 0); предел ячеек сработал %lld раз\n",
         (long long)S->nrevisit, (long long)S->ncap);
  printf("      спуск остановлен: ПУСТОТОЙ %lld (%.2f %%), ПОЛОМ %lld (%.2f %%), ЛИСТОМ %lld "
         "(%.2f %%)\n",
         (long long)S->stop_void, 100.0 * (double)S->stop_void / (double)S->ncell,
         (long long)S->stop_floor, 100.0 * (double)S->stop_floor / (double)S->ncell,
         (long long)S->stop_leaf, 100.0 * (double)S->stop_leaf / (double)S->ncell);
  printf("      ЯЧЕЕК ПО ОКТАВАМ РАССТОЯНИЯ (2^k м):");
  double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, sw = 0.0;
  for (int o = 0; o < NOCT; o++) {
    if (S->oct_cell[o] <= 0) continue;
    printf(" %d:%lld", o + OCT0, (long long)S->oct_cell[o]);
    double x = (double)(o + OCT0), y = log2((double)S->oct_cell[o]);
    sx += x;
    sy += y;
    sxx += x * x;
    sxy += x * y;
    sw += 1.0;
  }
  if (sw >= 2.0) {
    double den = sw * sxx - sx * sx;
    if (den > 0.0 || den < 0.0) printf("  -> наклон в log2 = %+.3f", (sw * sxy - sx * sy) / den);
  }
  printf("\n");
  printf("      СРЕДНИЙ УРОВЕНЬ ЯЧЕЙКИ ПО ОКТАВАМ (падает ли детальность с расстоянием):");
  for (int o = 0; o < NOCT; o++)
    if (S->oct_cell[o] > 0)
      printf(" %d:%.2f", o + OCT0, (double)S->oct_lev[o] / (double)S->oct_cell[o]);
  printf("\n");
  printf("      ПУСТЫХ ДЕТЕЙ У РОДИТЕЛЯ С ГЕОМЕТРИЕЙ (0..8):");
  {
    int64_t tot = 0, sum = 0;
    for (int k = 0; k <= 8; k++) {
      tot += S->sib_empty[k];
      sum += (int64_t)k * S->sib_empty[k];
    }
    for (int k = 0; k <= 8; k++)
      if (S->sib_empty[k] > 0)
        printf(" %d:%.3f", k, (double)S->sib_empty[k] / (double)(tot > 0 ? tot : 1));
    printf("  -> в среднем %.2f из 8\n", tot > 0 ? (double)sum / (double)tot : 0.0);
  }
  printf("      ДОЛЯ ПУСТЫХ ПО ОКТАВАМ:");
  for (int o = 0; o < NOCT; o++)
    if (S->oct_cell[o] > 0)
      printf(" %d:%.2f", o + OCT0, (double)S->oct_empty[o] / (double)S->oct_cell[o]);
  printf("\n");
  printf("      РАЗМЕР ЯЧЕЙКИ, log2 м:");
  for (int i = 0; i < NSZ; i++)
    if (S->sz[i] > 0) printf(" %d:%.3f", i + SZ0, (double)S->sz[i] / (double)S->ncell);
  printf("\n");
  if (S->ndeep > 0)
    printf("      ВИРТУАЛЬНЫХ ЯЧЕЕК, УПЁРШИХСЯ В ПРЕДЕЛ УРОВНЯ %d: %lld — критерий пола у ячейки "
           "с источником не выполняется НИКОГДА\n",
           HZ_PFRONT_MAXLEV, (long long)S->ndeep);
  printf("      ячеек, задевающих две октавы, %lld (%.2f %%) — оценка А507\n",
         (long long)S->nstraddle, 100.0 * (double)S->nstraddle / (double)S->ncell);
  printf("      ВРЕМЯ %.2f с, на ячейку %.1f нс\n", secs,
         S->ncell > 0 ? 1e9 * secs / (double)S->ncell : 0.0);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: pfrontx ФАЙЛ.obj МАСШТАБ [leaf=N] [lev=N] [grade=0|1] [src=x,y,z] "
                    "[px=4] [voidfloor=1] [bboxocc=1]\n");
    return 1;
  }
  int leafmax = 0, maxlev = 0, grade = 1, voidfloor = 0, bboxocc = 0, hassrc = 0, lean = 0,
      compact = 0, sun = 0;
  double px = 4.0, src[3] = {0.0, 0.0, 0.0}, sdir[3] = {0.0, -1.0, 0.0};
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "leaf=", 5) == 0) leafmax = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "lev=", 4) == 0) maxlev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "grade=", 6) == 0) grade = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "px=", 3) == 0) px = strtod(argv[i] + 3, NULL);
    if (strncmp(argv[i], "voidfloor=", 10) == 0) voidfloor = (int)strtol(argv[i] + 10, NULL, 10);
    if (strncmp(argv[i], "bboxocc=", 8) == 0) bboxocc = (int)strtol(argv[i] + 8, NULL, 10);
    if (strncmp(argv[i], "lean=", 5) == 0) lean = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "compact=", 8) == 0) compact = (int)strtol(argv[i] + 8, NULL, 10);
    if (strncmp(argv[i], "sun=", 4) == 0) {
      char *e = NULL;
      sdir[0] = strtod(argv[i] + 4, &e);
      if (e != NULL && *e == 44) sdir[1] = strtod(e + 1, &e);
      if (e != NULL && *e == 44) sdir[2] = strtod(e + 1, NULL);
      sun = 1;
    }
    if (strncmp(argv[i], "src=", 4) == 0) {
      char *e = NULL;
      src[0] = strtod(argv[i] + 4, &e);
      if (e != NULL && *e == ',') src[1] = strtod(e + 1, &e);
      if (e != NULL && *e == ',') src[2] = strtod(e + 1, NULL);
      hassrc = 1;
    }
  }
  hz_objmesh m;
  double t0 = now_s();
  if (hz_obj_load(&m, argv[1], strtod(argv[2], NULL)) != 0) {
    fprintf(stderr, "нет сцены %s\n", argv[1]);
    return 1;
  }
  printf("== СЦЕНА %s: треугольников %d, габарит %.2f x %.2f x %.2f м, чтение %.1f с\n", argv[1],
         m.nt, m.hi[0] - m.lo[0], m.hi[1] - m.lo[1], m.hi[2] - m.lo[2], now_s() - t0);
  hz_ptree T;
  t0 = now_s();
  if (hz_ptree_build_cut(&T, &m, leafmax, maxlev, grade) != 0) {
    fprintf(stderr, "отказ дерева\n");
    hz_obj_free(&m);
    return 2;
  }
  printf("== ДЕРЕВО за %.2f с: узлов %d, листьев %lld, глубина %d\n", now_s() - t0, T.nnd,
         (long long)T.nleaf, T.depth);
  if (!hassrc)
    for (int c = 0; c < 3; c++)
      src[c] = 0.5 * (T.nd[0].lo[c] + T.nd[0].hi[c]);
  printf("== УМОЛЧАНИЯ: leafmax %d, maxlev %d (0 — 16 и 12), градуировка %d, ИСТОЧНИК %.2f,%.2f,"
         "%.2f (%s), ε %.4e, пол %g px ОТ ИСТОЧНИКА (А508), пол в пустоте %d, занятость %s\n",
         leafmax, maxlev, grade, src[0], src[1], src[2], hassrc ? "ключ src=" : "центр сцены",
         HZ_CFG_EPS, px, voidfloor, bboxocc ? "ГАБАРИТНАЯ (верхняя оценка)" : "точная (О70)");

  double *tlo = malloc(3 * (size_t)m.nt * sizeof *tlo);
  double *thi = malloc(3 * (size_t)m.nt * sizeof *thi);
  double *tpl = malloc(4 * (size_t)m.nt * sizeof *tpl);
  int32_t *list = malloc((size_t)m.nt * sizeof *list);
  unsigned char *occ = calloc((size_t)T.nnd, 1);
  unsigned char *bocc = calloc((size_t)T.nnd, 1);
  unsigned char *seen = calloc((size_t)T.nnd, 1);
  if (tlo == NULL || thi == NULL || tpl == NULL || list == NULL || occ == NULL || bocc == NULL ||
      seen == NULL) {
    free(tlo);
    free(thi);
    free(tpl);
    free(list);
    free(occ);
    free(bocc);
    free(seen);
    hz_ptree_free(&T);
    hz_obj_free(&m);
    return 2;
  }
  for (int32_t t = 0; t < m.nt; t++) {
    const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
    const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
    const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
    double u[3], w[3], nr[3];
    for (int c = 0; c < 3; c++) {
      u[c] = B[c] - A[c];
      w[c] = C[c] - A[c];
    }
    nr[0] = u[1] * w[2] - u[2] * w[1];
    nr[1] = u[2] * w[0] - u[0] * w[2];
    nr[2] = u[0] * w[1] - u[1] * w[0];
    double L = sqrt(nr[0] * nr[0] + nr[1] * nr[1] + nr[2] * nr[2]);
    if (!(L > 0.0)) L = 1.0;
    double *pl = tpl + 4 * (size_t)t;
    for (int c = 0; c < 3; c++)
      pl[c] = nr[c] / L;
    pl[3] = pl[0] * A[0] + pl[1] * A[1] + pl[2] * A[2];
    double *bl = tlo + 3 * (size_t)t, *bh = thi + 3 * (size_t)t;
    for (int c = 0; c < 3; c++) {
      bl[c] = 1e300;
      bh[c] = -1e300;
    }
    for (int i = 0; i < 3; i++) {
      const double *p = m.v + 3 * (size_t)m.f[3 * (size_t)t + (size_t)i];
      for (int c = 0; c < 3; c++) {
        if (p[c] < bl[c]) bl[c] = p[c];
        if (p[c] > bh[c]) bh[c] = p[c];
      }
    }
    list[t] = t;
  }
  occctx OX;
  memset(&OX, 0, sizeof OX);
  OX.m = &m;
  OX.T = &T;
  OX.tlo = tlo;
  OX.thi = thi;
  OX.tpl = tpl;
  OX.occ = occ;
  OX.bocc = bocc;
  t0 = now_s();
  occ_walk(&OX, 0, list, m.nt);
  for (int32_t i = T.nnd - 1; i >= 0; i--) {
    int32_t c0 = T.nd[i].child;
    if (c0 < 0) continue;
    for (int k = 0; k < 8; k++) {
      occ[i] = (unsigned char)(occ[i] | occ[c0 + k]);
      bocc[i] = (unsigned char)(bocc[i] | bocc[c0 + k]);
    }
  }
  int64_t no = 0, nb = 0;
  for (int32_t i = 0; i < T.nnd; i++) {
    if (occ[i]) no++;
    if (bocc[i]) nb++;
  }
  printf("== ЗАНЯТОСТЬ за %.2f с: точная %lld узлов (%.2f %%), габаритная %lld (%.2f %%) — "
         "верхняя оценка, А511\n",
         now_s() - t0, (long long)no, 100.0 * (double)no / (double)T.nnd, (long long)nb,
         100.0 * (double)nb / (double)T.nnd);
  free(tlo);
  free(thi);
  free(tpl);
  free(list);

  frontstat S;
  memset(&S, 0, sizeof S);
  S.T = &T;
  S.occ = bboxocc ? bocc : occ;
  S.seen = seen;
  memcpy(S.src, src, sizeof src);
  S.eps = HZ_CFG_EPS;
  S.px = px;
  S.pxeps2 = (px * HZ_CFG_EPS) * (px * HZ_CFG_EPS);
  S.lean = lean;
  S.sun = sun;
  if (sun) {
    double L = sqrt(sdir[0] * sdir[0] + sdir[1] * sdir[1] + sdir[2] * sdir[2]);
    if (!(L > 0.0)) L = 1.0;
    for (int c = 0; c < 3; c++)
      S.dir[c] = sdir[c] / L;
    /* Вход в сцену — наименьшая проекция угла коробки на направление. */
    S.t_entry = 1e300;
    for (int k = 0; k < 8; k++) {
      double t = 0.0;
      for (int c = 0; c < 3; c++)
        t += ((k & (1 << c)) ? T.nd[0].hi[c] : T.nd[0].lo[c]) * S.dir[c];
      if (t < S.t_entry) S.t_entry = t;
    }
    printf("== ПАРАЛЛЕЛЬНЫЙ ФРОНТ (солнце): направление %.3f,%.3f,%.3f; пол по ПРОЙДЕННОМУ ПУТИ от "
           "входа в сцену\n",
           S.dir[0], S.dir[1], S.dir[2]);
  }
  S.voidfloor = voidfloor;
  S.levmin = 99;
  S.levmax = -1;
  /* КОМПАКТНЫЙ ОБХОД ВЕТКИ СОЛНЦА НЕ ЗНАЕТ, И МОЛЧА ДАВАТЬ ТОЧЕЧНЫЙ ОТВЕТ ОН НЕ
   * БУДЕТ. Поймано тем, что `sun=…` при `compact=1` вернул числа точечного
   * источника ДО ПОСЛЕДНЕЙ ЦИФРЫ — та же подпись «параметр в сбор не вошёл»,
   * что в О67 и §253. Отказ вместо тихого неверного числа. */
  if (compact && sun) {
    fprintf(stderr, "compact=1 и sun=… вместе не считаются: компактный обход "
                    "ведёт пол от ТОЧКИ. Уберите compact=1.\n");
    free(occ);
    free(bocc);
    free(seen);
    hz_ptree_free(&T);
    hz_obj_free(&m);
    return 2;
  }
  t0 = now_s();
  if (compact) {
    int32_t *ch = malloc((size_t)T.nnd * sizeof *ch);
    if (ch == NULL) {
      free(occ);
      free(bocc);
      free(seen);
      hz_ptree_free(&T);
      hz_obj_free(&m);
      return 2;
    }
    for (int32_t i = 0; i < T.nnd; i++)
      ch[i] = T.nd[i].child;
    printf("== КОМПАКТНЫЙ ИНДЕКС: %d узлов по 4 Б = %.3f ГБ против %.3f ГБ полного\n", T.nnd,
           (double)T.nnd * 4.0 / 1073741824.0,
           (double)T.nnd * (double)sizeof(hz_ptnode) / 1073741824.0);
    t0 = now_s();
    compact_walk(&S, ch, 0, T.nd[0].lo, T.nd[0].hi);
    free(ch);
  } else
    front_walk(&S, 0);
  report(&S, now_s() - t0, voidfloor ? "ОБХОД С ПОЛОМ В ПУСТОТЕ (негативный контроль)" : "ОБХОД");
  free(occ);
  free(bocc);
  free(seen);
  hz_ptree_free(&T);
  hz_obj_free(&m);
  return 0;
}
