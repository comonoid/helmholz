/* pmarchx.c — ШАГ О71 (Ф2): ОБХОД И ПУСТОТА. План — §247, аудит — §248.
 *
 * ЧТО МЕРИТСЯ. Два свойства структуры О70, и оба проверяются ЗНАКОМ, а не
 * величиной (А468): касания как функция расстояния (на октаву обязано быть
 * постоянно, то есть итог логарифмичен) и распределение длины шага (обязано
 * быть двугорбым — длинные прыжки по пустоте, короткие у геометрии).
 *
 * ЧЕГО ЗДЕСЬ НЕТ СОЗНАТЕЛЬНО. Ни маски пустых детей, ни ленивого деления
 * верхних узлов, ни иной правки представления: указание пользователя 08-05 —
 * структура самоценности не имеет, точить её до замера марша нельзя (§247).
 *
 * ЗАПУСК:
 *   pmarchx ФАЙЛ.obj МАСШТАБ [leaf=N] [lev=N] [grade=0|1] [eye=x,y,z]
 *           [px=4] [fix=N] [rays=N] [uni=1]
 *   `fix=N` — срез по ФИКСИРОВАННОМУ уровню (А496: пол εR камерозависим, и при
 *   нём длинные шаги вдали выходят сами собой). `uni=1` — НЕГАТИВНЫЙ КОНТРОЛЬ:
 *   марш по равномерной сетке самого мелкого уровня, «запрет подъёма».
 */
#include "pclip.h"
#include "pmarch.h"
#include "ptree.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ОКТАВЫ РАССТОЯНИЯ. Отсчёт от 2^-8 м (≈4 мм) — ниже этого не опускается ни
 * одна ячейка ни на одной из трёх сцен (самая мелкая: зал, 8.12 м / 2^12 = 2 мм,
 * и она попадает в первую октаву). Тридцать две октавы покрывают 4 мм…1.7e7 м,
 * то есть с запасом всё, что может дать город. */
#define OCT0 (-8)
#define NOCT 32
/* ГИСТОГРАММА ДЛИНЫ ШАГА в log2 метров, тот же отсчёт и тот же довод. */
#define LEN0 (-24)
#define NLEN 48
/* ЛУЧЕЙ НА ИСТОЧНИК ПО УМОЛЧАНИЮ. Направления детерминированы (спираль по
 * золотому сечению), поэтому число задаёт не шум, а покрытие сферы: 4096 лучей
 * дают угловой шаг около 1.6°, что мельче любой особенности, которую шаг ищет. */
#define NRAY 4096
/* ДОЛЯ, НИЖЕ КОТОРОЙ МОДА НЕ СЧИТАЕТСЯ МОДОЙ (А499, назван ДО прогона). */
#define MODE_MIN 0.10

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int iszero_d(double x) {
  return !(x < 0.0) && !(x > 0.0);
}

/* ---- ЗАНЯТОСТЬ УЗЛА ------------------------------------------------------- */
/* Ячейка ПУСТА, если её коробку не задевает ни один кусок. Определение то же,
 * что у О70 («ячеек с геометрией»), поэтому числа двух шагов сравнимы. Считается
 * ОДИН раз обходом до листьев, потом поднимается вверх: если кусок есть у
 * ребёнка, он есть и у родителя, а обратное верно потому, что дети покрывают
 * родителя целиком. */
typedef struct {
  const hz_objmesh *m;
  const hz_ptree *T;
  const double *tlo, *thi, *tpl;
  unsigned char *occ;
  int fail;
} occctx;

static int box_hits(const double *bl, const double *bh, const double *clo, const double *chi) {
  for (int c = 0; c < 3; c++)
    if (!(bl[c] < chi[c]) || !(bh[c] >= clo[c])) return 0;
  return 1;
}

static void occ_walk(occctx *X, int32_t nid, const int32_t *list, int32_t n) {
  if (X->fail) return;
  const hz_ptnode *N = &X->T->nd[nid];
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

/* ---- СЧЁТЧИКИ ------------------------------------------------------------- */
typedef struct {
  int64_t nray, nstep, nvisit, nempty, ncap, nbreak, nmiss;
  int64_t nleafcut; /* шагов, где спуск остановлен ЛИСТОМ, а не полом */
  int64_t oct_step[NOCT], oct_visit[NOCT], oct_empty[NOCT], oct_ray[NOCT];
  int64_t oct_leafcut[NOCT];
  int64_t levhist[32]; /* уровни ячеек среза: где на самом деле идёт марш */
  int64_t len[NLEN];
  double covmax; /* наибольшая относительная невязка покрытия луча */
  double secs;
} mstat;

static int oct_of(double t) {
  if (!(t > 0.0)) return 0;
  int o = (int)floor(log2(t)) - OCT0;
  if (o < 0) o = 0;
  if (o >= NOCT) o = NOCT - 1;
  return o;
}

static int len_of(double L) {
  if (!(L > 0.0)) return 0;
  int i = (int)floor(log2(L)) - LEN0;
  if (i < 0) i = 0;
  if (i >= NLEN) i = NLEN - 1;
  return i;
}

/* Направления — спираль по золотому сечению. ДЕТЕРМИНИРОВАНЫ: случайности в
 * замере нет вовсе, повтор прогона обязан давать те же числа до бита. */
static void ray_dir(int i, int n, double *d) {
  double z = 1.0 - 2.0 * ((double)i + 0.5) / (double)n;
  double r = sqrt(z * z < 1.0 ? 1.0 - z * z : 0.0);
  double phi = 2.0 * 3.14159265358979323846 * (double)i * 0.6180339887498949;
  d[0] = r * cos(phi);
  d[1] = r * sin(phi);
  d[2] = z;
}

static void march_all(mstat *S, const hz_ptree *T, const unsigned char *occ,
                      const hz_pmarch_cut *cut, const double *org, int nray) {
  double t0 = now_s();
  for (int i = 0; i < nray; i++) {
    double d[3];
    ray_dir(i, nray, d);
    hz_pmarch M;
    if (!hz_pmarch_begin(&M, T, cut, org, d)) {
      S->nmiss++;
      continue;
    }
    S->nray++;
    double cov = 0.0;
    int seen[NOCT];
    memset(seen, 0, sizeof seen);
    do {
      double L = M.t1 - M.t0;
      cov += L;
      int o = oct_of(0.5 * (M.t0 + M.t1));
      S->nstep++;
      S->nvisit += M.visits;
      S->oct_step[o]++;
      S->oct_visit[o] += M.visits;
      if (!seen[o]) {
        seen[o] = 1;
        S->oct_ray[o]++;
      }
      if (occ[M.nid] == 0) {
        S->nempty++;
        S->oct_empty[o]++;
      }
      if (M.byleaf) {
        S->nleafcut++;
        S->oct_leafcut[o]++;
      }
      S->levhist[M.lev & 31]++;
      S->len[len_of(L)]++;
    } while (hz_pmarch_next(&M));
    S->ncap += M.ncap;
    S->nbreak += M.nbreak;
    double want = M.tout - M.tin;
    if (want > 0.0) {
      double e = fabs(cov - want) / want;
      if (e > S->covmax) S->covmax = e;
    }
  }
  S->secs += now_s() - t0;
}

/* НЕГАТИВНЫЙ КОНТРОЛЬ (§247): «запрет подъёма» — марш по РАВНОМЕРНОЙ сетке
 * самого мелкого уровня. Обязан дать касания, растущие ЛИНЕЙНО с расстоянием, и
 * ОДНОГОРБОЕ распределение шага. Это контроль над СЧЁТЧИКАМИ, а не над деревом:
 * он показывает, что гистограмма сообщает одну моду, когда мода одна, и что
 * счёт по октавам сообщает линейный рост, когда рост линеен. Без него
 * «двугорбо» и «логарифм» неотличимы от свойства печати. */
static void march_uniform(mstat *S, const hz_ptree *T, const double *org, int nray, int lev) {
  double t0 = now_s();
  double h[3], lo[3], hi[3];
  double s = (double)((int64_t)1 << lev);
  for (int c = 0; c < 3; c++) {
    lo[c] = T->nd[0].lo[c];
    hi[c] = T->nd[0].hi[c];
    h[c] = (hi[c] - lo[c]) / s;
  }
  for (int i = 0; i < nray; i++) {
    double d[3];
    ray_dir(i, nray, d);
    double a = 0.0, b = HUGE_VAL;
    int ok = 1;
    for (int c = 0; c < 3 && ok; c++) {
      if (iszero_d(d[c])) {
        if (org[c] < lo[c] || !(org[c] < hi[c])) ok = 0;
        continue;
      }
      double x = (lo[c] - org[c]) / d[c], y = (hi[c] - org[c]) / d[c];
      if (x > y) {
        double q = x;
        x = y;
        y = q;
      }
      if (x > a) a = x;
      if (y < b) b = y;
      if (a > b) ok = 0;
    }
    if (!ok) {
      S->nmiss++;
      continue;
    }
    S->nray++;
    int seen[NOCT];
    memset(seen, 0, sizeof seen);
    double t = a;
    int64_t guard = 0;
    while (t < b && guard < HZ_PMARCH_MAXSTEP) {
      double te = b;
      for (int c = 0; c < 3; c++) {
        if (iszero_d(d[c])) continue;
        double x = org[c] + t * d[c];
        double g = floor((x - lo[c]) / h[c]);
        double bnd = lo[c] + (d[c] > 0.0 ? g + 1.0 : g) * h[c];
        double tb = (bnd - org[c]) / d[c];
        if (!(tb > t)) tb = t + h[c] / fabs(d[c]);
        if (tb < te) te = tb;
      }
      int o = oct_of(0.5 * (t + te));
      S->nstep++;
      S->oct_step[o]++;
      if (!seen[o]) {
        seen[o] = 1;
        S->oct_ray[o]++;
      }
      S->len[len_of(te - t)]++;
      t = te;
      guard++;
    }
    if (guard >= HZ_PMARCH_MAXSTEP) S->ncap++;
  }
  S->secs += now_s() - t0;
}

/* ---- ПЕЧАТЬ --------------------------------------------------------------- */
static void report(const mstat *S, const char *name, int uni) {
  printf("   %s: лучей %lld (мимо корня %lld), ШАГОВ %lld", name, (long long)S->nray,
         (long long)S->nmiss, (long long)S->nstep);
  if (!uni) printf(", ПОСЕЩЕНИЙ УЗЛОВ %lld", (long long)S->nvisit);
  printf("\n");
  if (S->nray <= 0) return;
  printf("      на луч: шагов %.1f", (double)S->nstep / (double)S->nray);
  if (!uni)
    printf(", посещений %.1f, посещений на шаг %.2f", (double)S->nvisit / (double)S->nray,
           S->nstep > 0 ? (double)S->nvisit / (double)S->nstep : 0.0);
  printf("\n");
  if (!uni) {
    printf("      шагов по ПУСТЫМ ячейкам %lld (%.2f %%)\n", (long long)S->nempty,
           S->nstep > 0 ? 100.0 * (double)S->nempty / (double)S->nstep : 0.0);
    /* ДИАГНОСТИКА ПУСТОГО СКАНА (урок О67). Если спуск почти всегда упирается в
     * ЛИСТ, то пол среза до дела не доходит, и «скан по полу ничего не меняет»
     * есть свойство ДЕРЕВА (глубина исчерпана), а не свойство пола. */
    printf("      спуск остановлен ЛИСТОМ %lld (%.2f %%), ПОЛОМ %lld (%.2f %%)\n",
           (long long)S->nleafcut, 100.0 * (double)S->nleafcut / (double)S->nstep,
           (long long)(S->nstep - S->nleafcut),
           100.0 * (double)(S->nstep - S->nleafcut) / (double)S->nstep);
    printf("      УРОВНИ ЯЧЕЕК СРЕЗА:");
    for (int l = 0; l < 32; l++)
      if (S->levhist[l] > 0) printf(" %d:%.3f", l, (double)S->levhist[l] / (double)S->nstep);
    printf("\n");
    printf("      ДОЛЯ ШАГОВ, ОСТАНОВЛЕННЫХ ЛИСТОМ, ПО ОКТАВАМ:");
    for (int o = 0; o < NOCT; o++)
      if (S->oct_step[o] > 0)
        printf(" %d:%.2f", o + OCT0, (double)S->oct_leafcut[o] / (double)S->oct_step[o]);
    printf("\n");
  }
  printf("      предел шагов сработал %lld раз; РАЗРЫВОВ ПОКРЫТИЯ %lld (обязано 0); наибольшая "
         "невязка покрытия %.3e\n",
         (long long)S->ncap, (long long)S->nbreak, S->covmax);
  /* КАСАНИЯ ПО ОКТАВАМ РАССТОЯНИЯ — главное число шага (П1). Печатается ШАГОВ
   * НА ЛУЧ В ОКТАВЕ: делить надо на число лучей, ДОШЕДШИХ до октавы, иначе
   * дальние октавы окажутся занижены просто потому, что туда дошли не все. */
  printf("      ШАГОВ НА ЛУЧ ПО ОКТАВАМ (2^k м):");
  double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, sw = 0.0;
  for (int o = 0; o < NOCT; o++) {
    if (S->oct_ray[o] <= 0) continue;
    double f = (double)S->oct_step[o] / (double)S->oct_ray[o];
    printf(" %d:%.1f", o + OCT0, f);
    if (f > 0.0) {
      double x = (double)(o + OCT0), y = log2(f);
      sx += x;
      sy += y;
      sxx += x * x;
      sxy += x * y;
      sw += 1.0;
    }
  }
  if (sw >= 2.0) {
    double den = sw * sxx - sx * sx;
    if (!iszero_d(den)) printf("  -> наклон в log2 = %+.3f", (sw * sxy - sx * sy) / den);
  }
  printf("\n");
  if (!uni) {
    printf("      ДОЛЯ ПУСТЫХ ШАГОВ ПО ОКТАВАМ:");
    for (int o = 0; o < NOCT; o++) {
      if (S->oct_step[o] <= 0) continue;
      printf(" %d:%.2f", o + OCT0, (double)S->oct_empty[o] / (double)S->oct_step[o]);
    }
    printf("\n");
  }
  /* ДЛИНА ШАГА: гистограмма, число мод и разнос (А499 — критерий назван ДО
   * прогона: моды считаются локальными максимумами с долей не ниже MODE_MIN, и
   * приёмка требует двух таких мод, разнесённых не менее чем на порядок). */
  printf("      ДЛИНА ШАГА, log2 м:");
  for (int i = 0; i < NLEN; i++)
    if (S->len[i] > 0)
      printf(" %d:%.3f", i + LEN0, (double)S->len[i] / (double)(S->nstep > 0 ? S->nstep : 1));
  printf("\n");
  int nmode = 0, m1 = 0, m2 = 0;
  for (int i = 0; i < NLEN; i++) {
    double f = (double)S->len[i] / (double)(S->nstep > 0 ? S->nstep : 1);
    if (!(f >= MODE_MIN)) continue;
    double a = (i > 0) ? (double)S->len[i - 1] : 0.0;
    double b = (i + 1 < NLEN) ? (double)S->len[i + 1] : 0.0;
    if ((double)S->len[i] >= a && (double)S->len[i] >= b) {
      if (nmode == 0)
        m1 = i;
      else
        m2 = i;
      nmode++;
    }
  }
  printf("      МОД (доля ≥ %.0f %%) %d", 100.0 * MODE_MIN, nmode);
  if (nmode >= 2) printf(", разнос %d октав (обязано ≥ 3, то есть порядок)", m2 - m1);
  printf("\n");
  printf("      ВРЕМЯ %.2f с, на шаг %.1f нс\n", S->secs,
         S->nstep > 0 ? 1e9 * S->secs / (double)S->nstep : 0.0);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: pmarchx ФАЙЛ.obj МАСШТАБ [leaf=N] [lev=N] [grade=0|1] [eye=x,y,z] "
                    "[px=4] [fix=N] [rays=N] [uni=1]\n");
    return 1;
  }
  int leafmax = 0, maxlev = 0, grade = 1, fixlev = -1, nray = NRAY, uni = 0, haseye = 0;
  double pxcut = 4.0, eye[3] = {0.0, 0.0, 0.0};
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "leaf=", 5) == 0) leafmax = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "lev=", 4) == 0) maxlev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "grade=", 6) == 0) grade = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "px=", 3) == 0) pxcut = strtod(argv[i] + 3, NULL);
    if (strncmp(argv[i], "fix=", 4) == 0) fixlev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "rays=", 5) == 0) nray = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "uni=", 4) == 0) uni = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "eye=", 4) == 0) {
      char *e = NULL;
      eye[0] = strtod(argv[i] + 4, &e);
      if (e != NULL && *e == ',') eye[1] = strtod(e + 1, &e);
      if (e != NULL && *e == ',') eye[2] = strtod(e + 1, NULL);
      haseye = 1;
    }
  }
  const char *eyesrc = "ключ eye=";
  if (!haseye) {
    /* Зал удалён 08-11; его камера теперь у комнаты-инструмента (scene_cfg.h). */
    if (strstr(argv[1], "room") != NULL) {
      double e0[3] = HZ_CFG_HALL_EYE;
      memcpy(eye, e0, sizeof eye);
      eyesrc = "scene_cfg.h HALL_EYE";
    } else if (strstr(argv[1], "rungholt") != NULL) {
      double e0[3] = HZ_CFG_CITY_EYE;
      memcpy(eye, e0, sizeof eye);
      eyesrc = "scene_cfg.h CITY_EYE";
    } else {
      double e0[3] = HZ_CFG_MIGUEL_EYE;
      memcpy(eye, e0, sizeof eye);
      eyesrc = "scene_cfg.h MIGUEL_EYE";
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
  printf("== УМОЛЧАНИЯ: leafmax %d, maxlev %d (0 — 16 и 12), градуировка %d, глаз %.2f,%.2f,%.2f "
         "(%s), ε %.4e рад/пиксель, лучей %d, пол %g px, фикс. уровень %d\n",
         leafmax, maxlev, grade, eye[0], eye[1], eye[2], eyesrc, HZ_CFG_EPS, nray, pxcut, fixlev);

  hz_ptree T;
  t0 = now_s();
  if (hz_ptree_build_cut(&T, &m, leafmax, maxlev, grade) != 0) {
    fprintf(stderr, "отказ дерева\n");
    hz_obj_free(&m);
    return 2;
  }
  printf(
      "== ДЕРЕВО за %.2f с (спуск %.2f + градуировка %.2f): узлов %d, листьев %lld, глубина %d\n",
      now_s() - t0, T.t_build, T.t_grade, T.nnd, (long long)T.nleaf, T.depth);

  /* ЗАНЯТОСТЬ. Считается один раз; определение то же, что «ячеек с геометрией»
   * в О70, поэтому доли сравнимы между шагами. */
  double *tlo = malloc(3 * (size_t)m.nt * sizeof *tlo);
  double *thi = malloc(3 * (size_t)m.nt * sizeof *thi);
  double *tpl = malloc(4 * (size_t)m.nt * sizeof *tpl);
  int32_t *list = malloc((size_t)m.nt * sizeof *list);
  unsigned char *occ = calloc((size_t)T.nnd, 1);
  if (tlo == NULL || thi == NULL || tpl == NULL || list == NULL || occ == NULL) {
    free(tlo);
    free(thi);
    free(tpl);
    free(list);
    free(occ);
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
  t0 = now_s();
  occ_walk(&OX, 0, list, m.nt);
  for (int32_t i = T.nnd - 1; i >= 0; i--) {
    int32_t c0 = T.nd[i].child;
    if (c0 < 0) continue;
    for (int k = 0; k < 8; k++)
      occ[i] = (unsigned char)(occ[i] | occ[c0 + k]);
  }
  int64_t nocc = 0, nleafocc = 0;
  for (int32_t i = 0; i < T.nnd; i++) {
    if (occ[i]) nocc++;
    if (T.nd[i].child < 0 && occ[i]) nleafocc++;
  }
  printf("== ЗАНЯТОСТЬ за %.2f с: узлов с геометрией %lld (%.2f %%), листьев с геометрией %lld "
         "(%.2f %%)\n",
         now_s() - t0, (long long)nocc, 100.0 * (double)nocc / (double)T.nnd, (long long)nleafocc,
         T.nleaf > 0 ? 100.0 * (double)nleafocc / (double)T.nleaf : 0.0);
  free(tlo);
  free(thi);
  free(tpl);
  free(list);
  if (OX.fail) {
    fprintf(stderr, "отказ обхода занятости\n");
    free(occ);
    hz_ptree_free(&T);
    hz_obj_free(&m);
    return 2;
  }

  /* ДВА ИСТОЧНИКА ЛУЧЕЙ (А497). Луч от КАМЕРЫ — не фронт: фронт идёт от
   * источника и туда, куда камера не смотрит. Меряются оба, и расхождение само
   * по себе есть результат. */
  double ctr[3];
  for (int c = 0; c < 3; c++)
    ctr[c] = 0.5 * (T.nd[0].lo[c] + T.nd[0].hi[c]);
  const double *org[2] = {eye, ctr};
  const char *orgname[2] = {"от КАМЕРЫ", "от ЦЕНТРА сцены"};

  if (uni) {
    printf("== НЕГАТИВНЫЙ КОНТРОЛЬ: РАВНОМЕРНАЯ СЕТКА уровня %d («запрет подъёма»)\n", T.depth);
    printf("   ОБЯЗАН дать линейный рост касаний с расстоянием и ОДНУ моду длины шага.\n");
    for (int s = 0; s < 2; s++) {
      mstat S;
      memset(&S, 0, sizeof S);
      march_uniform(&S, &T, org[s], nray, T.depth);
      report(&S, orgname[s], 1);
    }
  } else {
    for (int r = 0; r < 2; r++) {
      hz_pmarch_cut C;
      memset(&C, 0, sizeof C);
      memcpy(C.eye, eye, sizeof C.eye);
      C.eps = HZ_CFG_EPS;
      if (r == 0) {
        C.pxcut = pxcut;
        C.fixlev = -1;
        printf("== СРЕЗ по полу %g пикселя (камерозависим)\n", pxcut);
      } else {
        if (fixlev < 0) continue;
        C.pxcut = 0.0;
        C.fixlev = fixlev;
        printf("== СРЕЗ по ФИКСИРОВАННОМУ уровню %d (А496: двугорбость обязана пережить и его)\n",
               fixlev);
      }
      for (int s = 0; s < 2; s++) {
        mstat S;
        memset(&S, 0, sizeof S);
        march_all(&S, &T, occ, &C, org[s], nray);
        report(&S, orgname[s], 0);
      }
    }
  }
  free(occ);
  hz_ptree_free(&T);
  hz_obj_free(&m);
  return 0;
}
