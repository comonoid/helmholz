/* psil — СЧЁТ ВЫЖИВАЮЩИХ СИЛУЭТНЫХ РЁБЕР ПРИ ЗАДАННОМ МИНИМАЛЬНОМ РАЗМЕРЕ
 * ИСТОЧНИКА (PLAN_ELEMENTS.md, §149; аудит плана до кода — §150, А333…А338;
 * зачем шаг вообще — §151).
 *
 * ВОПРОС ОДНОЙ СТРОКОЙ. Направления в схеме тратятся не на поле (у ламбертовой
 * поверхности углового содержания нет), а на ВИДИМОСТЬ. Видимость кроится
 * силуэтами. Гасит ли минимальный угловой размер источника `α` комбинаторику
 * этого раскроя — то есть ограничен ли счёт силуэтов ВЫХОДОМ (пикселями), а не
 * ВХОДОМ (содержимым сцены)?
 *
 * ЧТО СЧИТАЕТСЯ. Для набора приёмников-точек (камера города §110 плюс `N`
 * случайных треугольников по зерну) каждое ребро сетки классифицируется:
 *   - ГРАНИЧНОЕ — одна смежная грань. От приёмника НЕ зависит вовсе;
 *   - СИЛУЭТНОЕ — две смежные грани, одна лицевая, другая изнаночная
 *     относительно приёмника. От приёмника зависит.
 * Это РАЗНЫЕ множества (А336), и печатаются они порознь: до этого замера число
 * «2.47 млн граничных рёбер» цитировалось как силуэтное, а основания к тому не
 * было.
 *
 * ПОЛУТЕНЬ ЗАПИСАНА ЯВНО (А335): `w = α·|x_ребро − x_приёмник|`, то есть `d`
 * меряется ОТ ЗАСЛОНА ДО ПРИЁМНИКА, а не от источника и не от начала фронта.
 * Сравнивается `w` с ВИДИМОЙ длиной ребра — длиной, взятой ПОПЕРЁК направления
 * на приёмник: `|e × û|`, где `û` — единичное направление ребро→приёмник.
 * Ребро, глядящее на приёмник с торца, видимой длины не имеет, и это верно:
 * тени оно не кроит.
 *
 * НЕМАНИФОЛДНЫЕ СТЫКИ НЕ ЗАМАЛЧИВАЮТСЯ. У Rungholt `9.2 %` рёбер принадлежат
 * более чем двум граням (§105.2, воксельные стыки). «Лицевая и изнаночная» для
 * них определено, но неоднозначно, поэтому они идут ОТДЕЛЬНОЙ строкой, а не
 * подмешиваются в главную таблицу.
 *
 * ЧЕГО ЗДЕСЬ НЕТ И ЧТО ЭТО ЗНАЧИТ ДЛЯ ЧТЕНИЯ ЧИСЕЛ. Нет ЗАСЛОНЕНИЯ: силуэт
 * классифицируется по ориентации граней, а не по тому, виден ли он из
 * приёмника. Значит все напечатанные числа — ВЕРХНЯЯ оценка: настоящий контур
 * короче на всё, что закрыто ближней геометрией. Заслонение стоит ровно того
 * марша по сцене, ради избавления от которого шаг и делается, поэтому оно здесь
 * сознательно не считается, а оговаривается.
 *
 * ЦЕНА ЗАМЕРА ПРЕДСКАЗАНА ДО ПРОГОНА (А337): смежность `1…3` мин, счёт
 * `2…10` мин. Смежность строится ВЁДРАМИ ПО МЛАДШЕЙ ВЕРШИНЕ (CSR), а не хешем и
 * не сортировкой: у `20.1` млн рёбер-слотов это один линейный проход плюс
 * разбор вёдер по шесть записей, то есть тот класс, который в А285 давал
 * двухчасовые зависания, здесь не заводится вовсе.
 *
 * Запуск: `psil [city|hall] [n=ЧИСЛО_ПРИЁМНИКОВ]`.
 */
#include "scene_cfg.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ЗЕРНО НАБОРА ПРИЁМНИКОВ. Любое фиксированное; названо ради повторяемости
 * (§2: числа разных прогонов обязаны быть сравнимы). */
#define PSIL_SEED 0x5EED1234ULL
/* ПРИЁМНИКОВ. Первый — камера города §110, остальные случайные. Число взято
 * тем же, на котором А337 считала время замера. */
#define PSIL_NRECV 32
/* ОКТАВЫ. Младшая `2^-4 = 6.25` см — вдвое мельче грани воксельной сцены
 * (`0.5 м²` на треугольник, §105.2), старшая `2^15 = 32.8` км — заведомо
 * больше габарита обеих имеющихся сцен (655 м у города). */
#define PSIL_KMIN (-4)
#define PSIL_NBIN 20
/* СКОЛЬКО ГРАНЕЙ РЕБРА ХРАНИТСЯ. У воксельного стыка их четыре; рёбра с
 * бо́льшим числом считаются отдельно, чтобы «не поместилось» не выглядело
 * как «не бывает». */
#define PSIL_MAXF 4
/* ПРИЁМКА (А333, вывод, а не потолок). Цель — 60 кадров при 4K, это 16 мс; на
 * 16 потоках при 1e9 простых операций в секунду на поток — 2.6e8 операций на
 * кадр; приёмников (полигонов города) 446 013, то есть ~580 операций на
 * приёмник; один член контурного интеграла — десятки операций, откуда бюджет
 * 10…100 рёбер. Порог ставится 1e3 — с запасом на иерархию приёмников. */
#define PSIL_ACC_SUM 1.0e3
/* Разброс по октавам: раскрой ограничен выходом, если в октавном слое число
 * выживших держится (предсказание §149 — множитель 3, порог приёмки 5). */
#define PSIL_ACC_SPREAD 5.0
/* УБИВАЕТ (§149): дальние октавы выше ближних впятеро при α = 1° ЛИБО сумма
 * выше 1e6. */
#define PSIL_KILL_SUM 1.0e6
/* ДИАПАЗОН, НАЗВАННЫЙ ПРЕДСКАЗАНИЕМ §149 дословно: «число выживших по октавам
 * ОТ 8 М ДО 256 М держится в пределах множителя 3». 8 м = 2³, 256 м = 2⁸,
 * то есть октавы k = 3…7. */
#define PSIL_KPRED_LO 3
#define PSIL_KPRED_HI 7
/* Ниже этого расстояния приёмник считается лежащим НА ребре: направления нет,
 * ребро не классифицируется и попадает в отдельный счётчик. Число —
 * не допуск задачи, а защита от деления на ноль. */
#define PSIL_DMIN 1.0e-12

/* α = 1°, 0.5°, 0.1° по §149; α = 0 — НЕГАТИВНЫЙ КОНТРОЛЬ (точечный источник,
 * отбора нет вовсе). Контроль гоняется тем же кодом и в том же проходе: разница
 * ровно в параметре, которым величина меняется. */
static const double PSIL_ALPHA_DEG[] = {1.0, 0.5, 0.1, 0.0};
#define PSIL_NALPHA ((int)(sizeof PSIL_ALPHA_DEG / sizeof PSIL_ALPHA_DEG[0]))
/* Порядок в списке ЗНАЧИМ, и опознаётся он ИНДЕКСОМ, а не сравнением double с
 * константой: `α` — параметр перебора, а не измеренная величина, но сравнение
 * плавающих на равенство запрещено гейтом всюду, и исключений не заводится. */
#define PSIL_Q_MAIN 0                /* α = 1°: по нему сформулировано «УБИВАЕТ» */
#define PSIL_Q_CTL (PSIL_NALPHA - 1) /* α = 0: негативный контроль */

/* ЭПСИЛОН ЛУЧА — ОТНОСИТЕЛЬНЫЙ, и он не новый: ровно то соглашение, что стоит
 * в `prender.c:599` (`len * (1.0 - 1e-6)`). Нужен для СКОЛЬЗЯЩИХ лучей вдоль
 * пола, где исключение по индексу не помогает (А345). */
#define PSIL_EPSREL 1.0e-6
/* Определитель Мёллера — Трумбора ниже этого считается вырожденным: у
 * метровых треугольников он порядка единицы, `1e-12` есть луч, лежащий в
 * плоскости с точностью до `1e-12` радиана. */
#define PSIL_DET_TINY 1.0e-12
#define PSIL_DTINY 1.0e-300
/* РАЗМЕР ЯЧЕЙКИ СЕТКИ выводится из СРЕДНЕЙ площади треугольника: ячейка не
 * должна быть мельче треугольника (иначе один треугольник растекается по
 * многим ячейкам) и не крупнее нескольких его размеров (иначе в ячейке
 * очередь). `cs = 2·sqrt(2·A_ср)`; у города `A = 0.5015` м² даёт `2.0` м. */
#define PSIL_CELL_MUL 2.0
/* Потолок числа ячеек: `start[]` есть `8` байт на ячейку, `64` млн ячеек —
 * `512` МБ. При превышении ячейка укрупняется. */
#define PSIL_CELL_MAX 67108864LL
/* КОЛЛИНЕАРНОСТЬ ПРИ СШИВКЕ (А346): синус угла. Координаты до `655` м дают
 * машинный ноль порядка `1e-13` относительно, `1e-9` — три порядка над ним и
 * одновременно нанометр отклонения на метровом ребре, то есть заведомо ниже
 * любой моделируемой детали. */
#define PSIL_COLL_SIN 1.0e-9
/* Сколько рёбер в вершине разбирается при сшивке. У воксельного стыка их
 * четыре; переполнение считается и печатается, а не молчит. */
#define PSIL_VDEG_MAX 32
/* Направлений в контроле луча П1. `4096` даёт долю неба с шумом `1/√N ≈ 1.6 %`
 * при сверяемой величине `15 %` — на порядок меньше допуска сверки `±5 %`. */
#define PSIL_NPROBE 4096

/* ---- ГРУППЫ ПРИЁМНИКОВ (§156, шаг О52) ----------------------------------- */
/* Центров на радиус и приёмников в центре. `8 × 16` при пяти радиусах даёт
 * `640` приёмников — по замеру §155 (`0.44` с на приёмник) это около пяти
 * минут, то есть внутри вилки, названной планом. */
#define PSIL_CL_NC 8
#define PSIL_CL_K 16
/* Последний радиус ОТРИЦАТЕЛЕН — это метка «вся сцена», негативный контроль. */
static const double PSIL_CL_R[] = {0.5, 2.0, 8.0, 32.0, -1.0};
#define PSIL_CL_NR ((int)(sizeof PSIL_CL_R / sizeof PSIL_CL_R[0]))
/* Потолок списка кандидатов в шаре: при `R = 32` м и плотности Rungholt
 * (`6.7` млн треугольников на `2.5e7` м³) в шаре порядка `3.7e4`; миллион —
 * тридцатикратный запас, а переполнение считается и печатается. */
#define PSIL_CL_CANDMAX (1 << 20)

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int cmp_dbl(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

static int cmp_i64(const void *a, const void *b) {
  int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

/* ---- РАВНОМЕРНАЯ СЕТКА ПО ТРЕУГОЛЬНИКАМ И ЛУЧ ПО НЕЙ --------------------
 *
 * Своя, а не `hz_pray` — довод в §153, пункт 2: `hz_pray` ищет попадание по
 * ПОЛИГОНАМ сегментации, отстоящим от треугольников на `δ`, и тогда «видно/не
 * видно» решала бы поверхность, отличная от той, что породила силуэт. */
typedef struct {
  int32_t nc[3];
  double lo[3], hi[3], cs;
  int64_t *start; /* ncell+1 */
  int32_t *idx;
  int64_t ncell, nref;
} psil_grid;

static void psil_tri_box(const hz_objmesh *m, int32_t t, double blo[3], double bhi[3]) {
  for (int a = 0; a < 3; a++) {
    blo[a] = 1e300;
    bhi[a] = -1e300;
  }
  for (int i = 0; i < 3; i++) {
    const double *p = m->v + 3 * (size_t)m->f[(size_t)t * 3 + (size_t)i];
    for (int a = 0; a < 3; a++) {
      if (p[a] < blo[a]) blo[a] = p[a];
      if (p[a] > bhi[a]) bhi[a] = p[a];
    }
  }
}

static int psil_grid_build(psil_grid *g, const hz_objmesh *m) {
  double asum = 0.0;
  for (int32_t t = 0; t < m->nt; t++)
    asum += hz_obj_tri_area(m, t);
  double cs = PSIL_CELL_MUL * sqrt(2.0 * asum / (double)m->nt);
  for (int a = 0; a < 3; a++) {
    /* Коробка расширяется на ячейку: точка ровно на грани сцены обязана иметь
     * ячейку, а не выпадать из сетки. */
    g->lo[a] = m->lo[a] - cs;
    g->hi[a] = m->hi[a] + cs;
  }
  for (;;) {
    int64_t n = 1;
    for (int a = 0; a < 3; a++) {
      double w = (g->hi[a] - g->lo[a]) / cs;
      int32_t k = (int32_t)floor(w) + 1;
      g->nc[a] = k < 1 ? 1 : k;
      n *= g->nc[a];
    }
    if (n <= PSIL_CELL_MAX) {
      g->ncell = n;
      break;
    }
    cs *= 2.0;
  }
  g->cs = cs;
  g->start = calloc((size_t)g->ncell + 2, sizeof *g->start);
  if (g->start == NULL) return 2;
  for (int pass = 0; pass < 2; pass++) {
    if (pass == 1) {
      for (int64_t i = 0; i < g->ncell; i++)
        g->start[i + 1] += g->start[i];
      g->nref = g->start[g->ncell];
      g->idx = malloc((size_t)(g->nref > 0 ? g->nref : 1) * sizeof *g->idx);
      if (g->idx == NULL) return 2;
      for (int64_t i = g->ncell; i > 0; i--)
        g->start[i] = g->start[i - 1];
      g->start[0] = 0;
    }
    for (int32_t t = 0; t < m->nt; t++) {
      double blo[3], bhi[3];
      psil_tri_box(m, t, blo, bhi);
      /* Инициализация нулём не «на всякий случай»: без неё gcc-analyzer не
       * связывает цикл по трём осям с последующим использованием и даёт
       * ложное «use of uninitialized value». */
      int32_t c0[3] = {0, 0, 0}, c1[3] = {0, 0, 0};
      for (int a = 0; a < 3; a++) {
        int32_t i0 = (int32_t)floor((blo[a] - g->lo[a]) / cs);
        int32_t i1 = (int32_t)floor((bhi[a] - g->lo[a]) / cs);
        c0[a] = i0 < 0 ? 0 : (i0 >= g->nc[a] ? g->nc[a] - 1 : i0);
        c1[a] = i1 < 0 ? 0 : (i1 >= g->nc[a] ? g->nc[a] - 1 : i1);
      }
      for (int32_t z = c0[2]; z <= c1[2]; z++)
        for (int32_t y = c0[1]; y <= c1[1]; y++)
          for (int32_t x = c0[0]; x <= c1[0]; x++) {
            int64_t c = (int64_t)x + (int64_t)g->nc[0] * ((int64_t)y + (int64_t)g->nc[1] * z);
            if (pass == 0)
              g->start[c + 1]++;
            else
              g->idx[g->start[c + 1]++] = t;
          }
    }
  }
  return 0;
}

static void psil_grid_free(psil_grid *g) {
  free(g->start);
  free(g->idx);
  g->start = NULL;
  g->idx = NULL;
}

static int psil_tri_hit(const hz_objmesh *m, int32_t t, const double o[3], const double d[3],
                        double *tout) {
  const double *p0 = m->v + 3 * (size_t)m->f[(size_t)t * 3 + 0];
  const double *p1 = m->v + 3 * (size_t)m->f[(size_t)t * 3 + 1];
  const double *p2 = m->v + 3 * (size_t)m->f[(size_t)t * 3 + 2];
  double e1[3], e2[3], pv[3], tv[3], qv[3];
  for (int a = 0; a < 3; a++) {
    e1[a] = p1[a] - p0[a];
    e2[a] = p2[a] - p0[a];
  }
  pv[0] = d[1] * e2[2] - d[2] * e2[1];
  pv[1] = d[2] * e2[0] - d[0] * e2[2];
  pv[2] = d[0] * e2[1] - d[1] * e2[0];
  double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
  if (fabs(det) < PSIL_DET_TINY) return 0; /* луч лежит в плоскости */
  double inv = 1.0 / det;
  for (int a = 0; a < 3; a++)
    tv[a] = o[a] - p0[a];
  double u = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
  if (u < 0.0 || u > 1.0) return 0;
  qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
  qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
  qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
  double v = (d[0] * qv[0] + d[1] * qv[1] + d[2] * qv[2]) * inv;
  if (v < 0.0 || u + v > 1.0) return 0;
  *tout = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
  return 1;
}

/* Обход DDA. `anyhit` — вопрос «есть ли заслон» (для тени), иначе ищется
 * БЛИЖАЙШЕЕ попадание (нужно контролю по горизонту). Двусторонний: стена
 * заслоняет независимо от ориентации. */
static int psil_trace(const psil_grid *g, const hz_objmesh *m, const double o[3],
                      const double dir[3], double tmin, double tmax, const int32_t *skip, int nskip,
                      int anyhit, double *thit) {
  double t0 = tmin, t1 = tmax;
  for (int a = 0; a < 3; a++) {
    if (fabs(dir[a]) < PSIL_DTINY) {
      if (o[a] < g->lo[a] || o[a] > g->hi[a]) return 0;
    } else {
      double ta = (g->lo[a] - o[a]) / dir[a], tb = (g->hi[a] - o[a]) / dir[a];
      if (ta > tb) {
        double s = ta;
        ta = tb;
        tb = s;
      }
      if (ta > t0) t0 = ta;
      if (tb < t1) t1 = tb;
    }
  }
  if (t0 > t1) return 0;
  int32_t ix[3], st[3];
  double tnext[3], tdel[3];
  for (int a = 0; a < 3; a++) {
    double p = o[a] + t0 * dir[a];
    int32_t k = (int32_t)floor((p - g->lo[a]) / g->cs);
    ix[a] = k < 0 ? 0 : (k >= g->nc[a] ? g->nc[a] - 1 : k);
    if (dir[a] > PSIL_DTINY) {
      st[a] = 1;
      tdel[a] = g->cs / dir[a];
      tnext[a] = (g->lo[a] + (double)(ix[a] + 1) * g->cs - o[a]) / dir[a];
    } else if (dir[a] < -PSIL_DTINY) {
      st[a] = -1;
      tdel[a] = -g->cs / dir[a];
      tnext[a] = (g->lo[a] + (double)ix[a] * g->cs - o[a]) / dir[a];
    } else {
      st[a] = 0;
      tdel[a] = 1e300;
      tnext[a] = 1e300;
    }
  }
  double best = tmax;
  int found = 0;
  for (;;) {
    int64_t c = (int64_t)ix[0] + (int64_t)g->nc[0] * ((int64_t)ix[1] + (int64_t)g->nc[1] * ix[2]);
    for (int64_t i = g->start[c]; i < g->start[c + 1]; i++) {
      int32_t t = g->idx[i];
      int sk = 0;
      for (int q = 0; q < nskip; q++)
        if (skip[q] == t) sk = 1;
      if (sk) continue;
      double th;
      if (!psil_tri_hit(m, t, o, dir, &th)) continue;
      if (th <= tmin || th >= best) continue;
      if (anyhit) {
        *thit = th;
        return 1;
      }
      best = th;
      found = 1;
    }
    int a =
        (tnext[0] < tnext[1]) ? ((tnext[0] < tnext[2]) ? 0 : 2) : ((tnext[1] < tnext[2]) ? 1 : 2);
    if (tnext[a] > t1 || (found && tnext[a] > best)) break;
    ix[a] += st[a];
    if (ix[a] < 0 || ix[a] >= g->nc[a]) break;
    tnext[a] += tdel[a];
  }
  if (found) *thit = best;
  return found;
}

/* ---- СШИВКА ВЫЖИВШИХ РЁБЕР В ПРЯМЫЕ УЧАСТКИ КОНТУРА (§153, пункт 3) ------
 *
 * Сто коллинеарных метровых рёбер вдоль конька — один член контурного
 * интеграла, а не сто. Объединяются рёбра, делящие вершину и коллинеарные в
 * ней; предел объединения (А346) — только коллинеарность, поэтому печатается
 * длина самой длинной цепочки. */
typedef struct {
  int32_t *uf, *hkey, *hhead, *hnxt, *roct, *cnt;
  double *rdist;
  int64_t hcap;
} psil_scratch;

/* ---- СТАТИСТИКА ОБЩНОСТИ ГРУППЫ (§156, §159) -----------------------------
 *
 * Один и тот же счёт нужен дважды — по всей группе и по каждой октаве, — и
 * дублировать его нельзя: разойдутся определения, а не только код. */
typedef struct {
  double sumS;  /* Σ|S_i| */
  double uni;   /* |объединение| */
  double J, Jw; /* мера Жаккара: средняя и худшая по парам */
  double zero;  /* доля пар с пустым пересечением СРЕДИ ЗРЯЧИХ */
  double own;   /* доля элементов, принадлежащих ровно одному приёмнику */
  int64_t nblind;
  int64_t npair;
} psil_stat;

static void psil_group_stat(const int32_t *const *S, const int64_t *n, int keff, int64_t *key,
                            int64_t keycap, psil_stat *o) {
  memset(o, 0, sizeof *o);
  int64_t nk = 0;
  for (int j = 0; j < keff; j++) {
    o->sumS += (double)n[j];
    if (n[j] == 0) o->nblind++;
    for (int64_t i = 0; i < n[j] && nk < keycap; i++)
      key[nk++] = (int64_t)S[j][i] * PSIL_CL_K + j;
  }
  if (nk == 0) return;
  qsort(key, (size_t)nk, sizeof *key, cmp_i64);
  int64_t uni = 0, own = 0;
  for (int64_t i = 0; i < nk;) {
    int64_t e = key[i] / PSIL_CL_K, j2 = i;
    while (j2 < nk && key[j2] / PSIL_CL_K == e)
      j2++;
    uni++;
    if (j2 - i == 1) own++;
    i = j2;
  }
  o->uni = (double)uni;
  o->own = (double)own / (double)nk;
  double js = 0.0, jw = 1.0;
  int64_t np = 0, nz = 0;
  for (int a = 0; a < keff; a++)
    for (int b = a + 1; b < keff; b++) {
      /* ПАРЫ СО СЛЕПЫМ ПРИЁМНИКОМ ИСКЛЮЧАЮТСЯ (А359): у замурованного в толще
       * блоков пересечение пусто ВСЕГДА, и, попав в счёт, он выдал бы «две
       * стороны одной стены» там, где стены нет вовсе. */
      if (n[a] == 0 || n[b] == 0) continue;
      int64_t ia = 0, ib = 0, inter = 0;
      while (ia < n[a] && ib < n[b]) {
        if (S[a][ia] == S[b][ib]) {
          inter++;
          ia++;
          ib++;
        } else if (S[a][ia] < S[b][ib])
          ia++;
        else
          ib++;
      }
      double jj = (double)inter / (double)(n[a] + n[b] - inter);
      js += jj;
      if (jj < jw) jw = jj;
      if (inter == 0) nz++;
      np++;
    }
  o->npair = np;
  o->J = (np > 0) ? js / (double)np : 0.0;
  o->Jw = (np > 0) ? jw : 0.0;
  o->zero = (np > 0) ? (double)nz / (double)np : 0.0;
}

/* ---- УГЛОВОЕ ОГРУБЛЕНИЕ КОНТУРА (§162, шаг О54) --------------------------
 *
 * Из множества рёбер строятся СВЯЗНЫЕ ЦЕПИ (общая вершина, разрез в вершинах
 * степени ≠ 2), и каждая упрощается по Дугласу — Пекеру с ПЕРЕМЕННЫМ допуском
 * `w(v) = α·(|v − центр| − R_блока)`: вершина выбрасывается, если отстоит от
 * хорды меньше собственной полутени. Допуск берётся по БЛИЖАЙШЕМУ возможному
 * приёмнику блока, а не по центру (А370) — иначе на ближних октавах он вдвое
 * велик. Нового параметра здесь нет: `α` тот же, что отбирает рёбра. */
typedef struct {
  int32_t *ord;   /* рёбра объединения в порядке обхода цепей */
  int32_t *coff;  /* nchord+1 границ хорд в `ord` */
  int32_t *voct;  /* октава хорды (от центра блока) */
  int64_t nchord; /* хорд после упрощения */
  int64_t nchain; /* связных цепей — потолок при допуске ∞ */
  int64_t ncollapse;
  double nwalk; /* рёбер пройдено обходом — обязано совпасть с |U| (разбиение) */
} psil_contour;

static void psil_contour_build(psil_contour *C, const int32_t *ue, int64_t nu, const int32_t *ea,
                               const int32_t *eb, const double *v, const double cen[3], double rblk,
                               double alpha, int32_t *hkey, int32_t *hhead, int32_t *hnxt,
                               int64_t hcap, int32_t *vseq, unsigned char *keep, int32_t *stk) {
  C->nchord = C->nchain = C->ncollapse = 0;
  C->nwalk = 0.0;
  if (nu <= 0) return;
  for (int64_t i = 0; i < hcap; i++)
    hkey[i] = -1;
  /* вершина -> список концов рёбер (2·i + j) */
  for (int64_t i = 0; i < nu; i++)
    for (int j = 0; j < 2; j++) {
      int32_t vv = (j == 0) ? ea[ue[i]] : eb[ue[i]];
      uint64_t h = (uint64_t)(uint32_t)vv * 0x9E3779B97F4A7C15ULL;
      int64_t s = (int64_t)(h >> 40) & (hcap - 1);
      while (hkey[s] >= 0 && hkey[s] != vv)
        s = (s + 1) & (hcap - 1);
      if (hkey[s] < 0) {
        hkey[s] = vv;
        hhead[s] = -1;
      }
      hnxt[2 * i + j] = hhead[s];
      hhead[s] = (int32_t)(2 * i + j);
    }
  unsigned char *used = keep + nu + 1; /* хвост общего буфера: nu байт */
  for (int64_t i = 0; i < nu; i++)
    used[i] = 0;
  int64_t nord = 0;
  /* Обход: сначала от вершин степени ≠ 2 (концы и ветвления, А372), затем
   * оставшиеся замкнутые циклы. */
  for (int pass = 0; pass < 2; pass++) {
    for (int64_t i0 = 0; i0 < nu; i0++) {
      for (int j0 = 0; j0 < 2; j0++) {
        if (used[i0]) continue;
        int32_t vstart = (j0 == 0) ? ea[ue[i0]] : eb[ue[i0]];
        if (pass == 0) {
          /* степень стартовой вершины */
          uint64_t h = (uint64_t)(uint32_t)vstart * 0x9E3779B97F4A7C15ULL;
          int64_t s = (int64_t)(h >> 40) & (hcap - 1);
          while (hkey[s] >= 0 && hkey[s] != vstart)
            s = (s + 1) & (hcap - 1);
          int deg = 0;
          for (int32_t p = hhead[s]; p >= 0; p = hnxt[p])
            deg++;
          if (deg == 2) continue;
        }
        /* пройти цепь */
        int64_t vn = 0, e0 = nord;
        vseq[vn++] = vstart;
        int64_t ce = i0;
        int32_t cv = vstart;
        for (;;) {
          used[ce] = 1;
          C->ord[nord++] = ue[ce];
          int32_t nv = (ea[ue[ce]] == cv) ? eb[ue[ce]] : ea[ue[ce]];
          vseq[vn++] = nv;
          uint64_t h = (uint64_t)(uint32_t)nv * 0x9E3779B97F4A7C15ULL;
          int64_t s = (int64_t)(h >> 40) & (hcap - 1);
          while (hkey[s] >= 0 && hkey[s] != nv)
            s = (s + 1) & (hcap - 1);
          int deg = 0;
          int64_t nxt = -1;
          for (int32_t p = hhead[s]; p >= 0; p = hnxt[p]) {
            deg++;
            int64_t ei = p >> 1;
            if (ei != ce && !used[ei]) nxt = ei;
          }
          if (deg != 2 || nxt < 0) break;
          ce = nxt;
          cv = nv;
        }
        C->nchain++;
        /* Дуглас — Пекер с переменным допуском, итеративно. */
        for (int64_t i = 0; i < vn; i++)
          keep[i] = 0;
        keep[0] = 1;
        keep[vn - 1] = 1;
        int64_t sp = 0;
        stk[sp++] = 0;
        stk[sp++] = (int32_t)(vn - 1);
        while (sp > 0) {
          int32_t hi = stk[--sp], lo = stk[--sp];
          if (hi <= lo + 1) continue;
          const double *pa = v + 3 * (size_t)vseq[lo], *pb = v + 3 * (size_t)vseq[hi];
          double ab[3], la = 0.0;
          for (int a = 0; a < 3; a++) {
            ab[a] = pb[a] - pa[a];
            la += ab[a] * ab[a];
          }
          la = sqrt(la);
          double worst = 0.0;
          int32_t wi = -1;
          for (int32_t t = lo + 1; t < hi; t++) {
            const double *pt = v + 3 * (size_t)vseq[t];
            double ap[3], cr[3], dist;
            for (int a = 0; a < 3; a++)
              ap[a] = pt[a] - pa[a];
            if (la > 0.0) {
              cr[0] = ap[1] * ab[2] - ap[2] * ab[1];
              cr[1] = ap[2] * ab[0] - ap[0] * ab[2];
              cr[2] = ap[0] * ab[1] - ap[1] * ab[0];
              dist = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]) / la;
            } else {
              dist = sqrt(ap[0] * ap[0] + ap[1] * ap[1] + ap[2] * ap[2]);
            }
            double dc = 0.0;
            for (int a = 0; a < 3; a++) {
              double w2 = pt[a] - cen[a];
              dc += w2 * w2;
            }
            double tol = alpha * (sqrt(dc) - rblk);
            if (tol < 0.0) tol = 0.0;
            double ex = dist - tol;
            if (ex > worst) {
              worst = ex;
              wi = t;
            }
          }
          if (wi < 0) continue; /* вся ломаная внутри полутени — хорда одна */
          keep[wi] = 1;
          stk[sp++] = lo;
          stk[sp++] = wi;
          stk[sp++] = wi;
          stk[sp++] = hi;
        }
        /* хорды: между соседними оставленными вершинами */
        int64_t prev = 0;
        for (int64_t t = 1; t < vn; t++) {
          if (!keep[t]) continue;
          C->coff[C->nchord] = (int32_t)(e0 + prev);
          const double *pa = v + 3 * (size_t)vseq[prev], *pb = v + 3 * (size_t)vseq[t];
          double dc = 0.0, ln = 0.0;
          for (int a = 0; a < 3; a++) {
            double mid = 0.5 * (pa[a] + pb[a]) - cen[a];
            dc += mid * mid;
            ln += (pb[a] - pa[a]) * (pb[a] - pa[a]);
          }
          int kk = (dc > 0.0) ? ilogb(sqrt(dc)) - PSIL_KMIN : 0;
          if (kk < 0) kk = 0;
          if (kk >= PSIL_NBIN) kk = PSIL_NBIN - 1;
          C->voct[C->nchord] = kk;
          /* А371: цепь, схлопнувшаяся в точку, — предмет тоньше полутени.
           * Физически законно, но обязано быть посчитано, а не случиться
           * молча. */
          if (sqrt(ln) <= 0.0) C->ncollapse++;
          C->nchord++;
          prev = t;
        }
        C->nwalk += (double)(vn - 1);
      }
    }
  }
  C->coff[C->nchord] = (int32_t)nord;
}

static int32_t psil_find(int32_t *uf, int32_t x) {
  while (uf[x] != x) {
    uf[x] = uf[uf[x]];
    x = uf[x];
  }
  return x;
}

static int64_t psil_chains(int64_t n, const int32_t *eidx, const int32_t *oct, const double *dist,
                           const int32_t *ea, const int32_t *eb, const double *v, psil_scratch *sc,
                           int64_t *octcnt, int64_t *chmax, int64_t *novf) {
  if (n <= 0) return 0;
  for (int64_t i = 0; i < n; i++)
    sc->uf[i] = (int32_t)i;
  int64_t hcap = sc->hcap;
  for (int64_t i = 0; i < hcap; i++)
    sc->hkey[i] = -1;
  for (int64_t i = 0; i < n; i++)
    for (int j = 0; j < 2; j++) {
      int32_t vv = (j == 0) ? ea[eidx[i]] : eb[eidx[i]];
      uint64_t h = (uint64_t)(uint32_t)vv * 0x9E3779B97F4A7C15ULL;
      int64_t s = (int64_t)(h >> 40) & (hcap - 1);
      while (sc->hkey[s] >= 0 && sc->hkey[s] != vv)
        s = (s + 1) & (hcap - 1);
      if (sc->hkey[s] < 0) {
        sc->hkey[s] = vv;
        sc->hhead[s] = -1;
      }
      sc->hnxt[2 * i + j] = sc->hhead[s];
      sc->hhead[s] = (int32_t)(2 * i + j);
    }
  for (int64_t s = 0; s < hcap; s++) {
    if (sc->hkey[s] < 0) continue;
    int32_t buf[PSIL_VDEG_MAX];
    int nb = 0;
    for (int32_t p = sc->hhead[s]; p >= 0; p = sc->hnxt[p]) {
      if (nb < PSIL_VDEG_MAX)
        buf[nb++] = p;
      else {
        (*novf)++;
        break;
      }
    }
    if (nb < 2) continue;
    const double *pv = v + 3 * (size_t)sc->hkey[s];
    for (int i = 0; i < nb; i++)
      for (int j = i + 1; j < nb; j++) {
        int64_t i1 = buf[i] >> 1, i2 = buf[j] >> 1;
        int32_t r1 = psil_find(sc->uf, (int32_t)i1), r2 = psil_find(sc->uf, (int32_t)i2);
        if (r1 == r2) continue;
        int32_t o1 = (buf[i] & 1) ? ea[eidx[i1]] : eb[eidx[i1]];
        int32_t o2 = (buf[j] & 1) ? ea[eidx[i2]] : eb[eidx[i2]];
        double u1[3], u2[3], cr[3];
        for (int a = 0; a < 3; a++) {
          u1[a] = v[3 * (size_t)o1 + (size_t)a] - pv[a];
          u2[a] = v[3 * (size_t)o2 + (size_t)a] - pv[a];
        }
        cr[0] = u1[1] * u2[2] - u1[2] * u2[1];
        cr[1] = u1[2] * u2[0] - u1[0] * u2[2];
        cr[2] = u1[0] * u2[1] - u1[1] * u2[0];
        double n1 = sqrt(u1[0] * u1[0] + u1[1] * u1[1] + u1[2] * u1[2]);
        double n2 = sqrt(u2[0] * u2[0] + u2[1] * u2[1] + u2[2] * u2[2]);
        double cn = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
        if (cn > PSIL_COLL_SIN * n1 * n2) continue;
        sc->uf[r1] = r2;
      }
  }
  for (int64_t i = 0; i < n; i++) {
    sc->cnt[i] = 0;
    sc->rdist[i] = 1e300;
  }
  for (int64_t i = 0; i < n; i++) {
    int32_t r = psil_find(sc->uf, (int32_t)i);
    sc->cnt[r]++;
    if (dist[i] < sc->rdist[r]) {
      sc->rdist[r] = dist[i];
      sc->roct[r] = oct[i];
    }
  }
  int64_t nseg = 0;
  for (int64_t i = 0; i < n; i++) {
    if (sc->cnt[i] == 0) continue;
    nseg++;
    if (sc->cnt[i] > *chmax) *chmax = sc->cnt[i];
    if (octcnt != NULL) octcnt[sc->roct[i]]++;
  }
  return nseg;
}

static uint64_t rnd64(uint64_t *s) {
  uint64_t x = *s;
  x ^= x << 13;
  x ^= x >> 7;
  x ^= x << 17;
  *s = x;
  return x;
}

static int64_t pct_i64(const int64_t *v, int n, double p) {
  if (n <= 0) return 0;
  int k = (int)(p * (double)(n - 1));
  if (k < 0) k = 0;
  if (k >= n) k = n - 1;
  return v[k];
}

int main(int argc, char **argv) {
  int city = 1, nrecv = PSIL_NRECV, noocc = 0, clust = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "hall") == 0) city = 0;
    if (strcmp(argv[i], "city") == 0) city = 1;
    /* НК1 (§153): заслонение ВЫКЛЮЧЕНО — числа обязаны совпасть с §152 точно. */
    if (strcmp(argv[i], "noocc") == 0) noocc = 1;
    /* §156, шаг О52: приёмники группами, замер общности видимого контура. */
    if (strcmp(argv[i], "clust") == 0) clust = 1;
    if (strncmp(argv[i], "n=", 2) == 0) nrecv = (int)strtol(argv[i] + 2, NULL, 10);
  }
  if (nrecv < 1) nrecv = 1;
  if (clust) {
    noocc = 0; /* без заслонения замерять общность видимого нечего */
    nrecv = PSIL_CL_NR * PSIL_CL_NC * PSIL_CL_K;
  }

  double t0 = now_s();
  hz_objmesh m;
  if (hz_obj_load(&m, city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ,
                  city ? HZ_CFG_CITY_SCALE : HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  double t_load = now_s() - t0;

  /* ГАБАРИТ — ПЕРВЫМ ЧИСЛОМ (А334). Октавы за габаритом дадут ноль оттого, что
   * там ничего нет, а не оттого, что критерий отобрал; без этой строки приёмка
   * «разброс по октавам» выполнилась бы сама собой. */
  double ext[3], diag = 0.0;
  for (int a = 0; a < 3; a++) {
    ext[a] = m.hi[a] - m.lo[a];
    diag += ext[a] * ext[a];
  }
  diag = sqrt(diag);
  int koct_max = (int)floor(log2(diag));
  printf("== СЦЕНА: %s, треугольников %d, вершин %d\n", city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ,
         m.nt, m.nv);
  printf("== ГАБАРИТ: %.2f x %.2f x %.2f м, диагональ %.2f м -> октавы внутри сцены: "
         "k = %d…%d (%d октав, НЕ 13)\n",
         ext[0], ext[1], ext[2], diag, PSIL_KMIN, koct_max, koct_max - PSIL_KMIN + 1);
  printf("== загрузка %.1f с\n", t_load);
  fflush(stdout);

  /* ---- СМЕЖНОСТЬ РЁБЕР (пункт 1 порядка работ, А336) --------------------- */
  double t1 = now_s();
  int64_t nslot = (int64_t)m.nt * 3;
  int64_t nvp1 = (int64_t)m.nv + 1;
  int64_t *off = calloc((size_t)nvp1 + 1, sizeof *off);
  if (off == NULL) return 2;
  for (int32_t t = 0; t < m.nt; t++)
    for (int i = 0; i < 3; i++) {
      int32_t a = m.f[(size_t)t * 3 + (size_t)i];
      int32_t b = m.f[(size_t)t * 3 + (size_t)((i + 1) % 3)];
      int32_t lo = a < b ? a : b;
      off[(int64_t)lo + 1]++;
    }
  for (int64_t i = 0; i < nvp1; i++)
    off[i + 1] += off[i];
  int32_t *shi = malloc((size_t)nslot * sizeof *shi);
  int32_t *sfa = malloc((size_t)nslot * sizeof *sfa);
  int64_t *cur = malloc((size_t)nvp1 * sizeof *cur);
  unsigned char *mark = calloc((size_t)nslot, 1);
  if (shi == NULL || sfa == NULL || cur == NULL || mark == NULL) return 2;
  for (int64_t i = 0; i < nvp1; i++)
    cur[i] = off[i];
  for (int32_t t = 0; t < m.nt; t++)
    for (int i = 0; i < 3; i++) {
      int32_t a = m.f[(size_t)t * 3 + (size_t)i];
      int32_t b = m.f[(size_t)t * 3 + (size_t)((i + 1) % 3)];
      int32_t lo = a < b ? a : b, hi = a < b ? b : a;
      int64_t p = cur[lo]++;
      shi[p] = hi;
      sfa[p] = t;
    }

  /* Уникальные рёбра: в каждом ведре линейная группировка по старшей вершине.
   * Вёдер `nv`, записей в ведре в среднем `3·nt/nv` (у города около шести), так
   * что `O(k²)` внутри ведра дёшево и предсказуемо. */
  int32_t *ea = malloc((size_t)nslot * sizeof *ea);
  int32_t *eb = malloc((size_t)nslot * sizeof *eb);
  int32_t *ef = malloc((size_t)nslot * PSIL_MAXF * sizeof *ef);
  unsigned char *enf = malloc((size_t)nslot);
  if (ea == NULL || eb == NULL || ef == NULL || enf == NULL) return 2;
  int64_t ne = 0, nb1 = 0, nb2 = 0, nbm = 0, nover = 0;
  for (int32_t a = 0; a < m.nv; a++) {
    int64_t s = off[a], e = off[a + 1];
    for (int64_t i = s; i < e; i++) {
      if (mark[i]) continue;
      int32_t h = shi[i];
      int64_t nf = 1;
      ea[ne] = a;
      eb[ne] = h;
      ef[ne * PSIL_MAXF] = sfa[i];
      for (int64_t j = i + 1; j < e; j++) {
        if (mark[j] || shi[j] != h) continue;
        mark[j] = 1;
        if (nf < PSIL_MAXF) ef[ne * PSIL_MAXF + nf] = sfa[j];
        nf++;
      }
      if (nf > PSIL_MAXF) nover++;
      enf[ne] = (unsigned char)(nf > 255 ? 255 : nf);
      if (nf == 1)
        nb1++;
      else if (nf == 2)
        nb2++;
      else
        nbm++;
      ne++;
    }
  }
  free(mark);
  free(cur);
  free(off);
  free(shi);
  free(sfa);
  double t_adj = now_s() - t1;
  printf("== СМЕЖНОСТЬ (%.1f с): рёбер всего %lld\n", t_adj, (long long)ne);
  printf("   ГРАНИЧНЫХ (одна грань, от точки НЕ зависят) %lld (%.2f %%)\n", (long long)nb1,
         100.0 * (double)nb1 / (double)ne);
  printf("   манифолдных (две грани) %lld (%.2f %%), неманифолдных (больше двух) %lld (%.2f %%), "
         "из них с более чем %d гранями %lld\n",
         (long long)nb2, 100.0 * (double)nb2 / (double)ne, (long long)nbm,
         100.0 * (double)nbm / (double)ne, PSIL_MAXF, (long long)nover);
  fflush(stdout);

  /* Плоскость треугольника: нормаль (НЕ нормированная — нужен только знак) и
   * смещение. Знак `n·p − off` и есть «лицевая относительно приёмника». */
  double *pn = malloc((size_t)m.nt * 4 * sizeof *pn);
  if (pn == NULL) return 2;
  for (int32_t t = 0; t < m.nt; t++) {
    const double *p0 = m.v + 3 * (size_t)m.f[(size_t)t * 3 + 0];
    const double *p1 = m.v + 3 * (size_t)m.f[(size_t)t * 3 + 1];
    const double *p2 = m.v + 3 * (size_t)m.f[(size_t)t * 3 + 2];
    double u[3], w[3], c[3];
    for (int a = 0; a < 3; a++) {
      u[a] = p1[a] - p0[a];
      w[a] = p2[a] - p0[a];
    }
    c[0] = u[1] * w[2] - u[2] * w[1];
    c[1] = u[2] * w[0] - u[0] * w[2];
    c[2] = u[0] * w[1] - u[1] * w[0];
    pn[(size_t)t * 4 + 0] = c[0];
    pn[(size_t)t * 4 + 1] = c[1];
    pn[(size_t)t * 4 + 2] = c[2];
    pn[(size_t)t * 4 + 3] = c[0] * p0[0] + c[1] * p0[1] + c[2] * p0[2];
  }

  /* ---- СЕТКА И КОНТРОЛЬ ЛУЧА ЧУЖИМ ЧИСЛОМ (А344) ------------------------- */
  psil_grid grid;
  memset(&grid, 0, sizeof grid);
  if (!noocc) {
    double tg = now_s();
    if (psil_grid_build(&grid, &m) != 0) return 2;
    printf("== СЕТКА (%.1f с): ячейка %.3f м, %d x %d x %d = %lld ячеек, ссылок %lld (%.2f на "
           "треугольник)\n",
           now_s() - tg, grid.cs, grid.nc[0], grid.nc[1], grid.nc[2], (long long)grid.ncell,
           (long long)grid.nref, (double)grid.nref / (double)m.nt);
    /* П1. Доля неба — по КОСИНУСНОМУ весу в верхней полусфере (полукуб §110
     * считает `Σ dff`, поэтому направления берутся косинусно-взвешенными и вес
     * получается равномерным). Горизонт — по ГОРИЗОНТАЛЬНЫМ направлениям, как
     * нижняя строка боковых граней полукуба. Оба числа сняты §110 ДРУГОЙ
     * машинерией, и в этом весь смысл сверки. */
    double eyec[3] = HZ_CFG_CITY_EYE, eyeh[3] = HZ_CFG_HALL_EYE;
    const double *ey = city ? eyec : eyeh;
    double eye[3] = {ey[0], ey[1], ey[2]};
    int64_t nsky = 0;
    double *hd = malloc((size_t)PSIL_NPROBE * sizeof *hd);
    if (hd == NULL) return 2;
    double far = 0.0;
    for (int a = 0; a < 3; a++)
      far += (grid.hi[a] - grid.lo[a]) * (grid.hi[a] - grid.lo[a]);
    far = sqrt(far);
    for (int i = 0; i < PSIL_NPROBE; i++) {
      double u1 = ((double)i + 0.5) / PSIL_NPROBE;
      double u2 = fmod((double)i * 0.6180339887498949, 1.0);
      double rr = sqrt(u1), th = 2.0 * M_PI * u2;
      double d[3] = {rr * cos(th), sqrt(1.0 - u1), rr * sin(th)};
      double t;
      if (!psil_trace(&grid, &m, eye, d, 0.0, far, NULL, 0, 1, &t)) nsky++;
    }
    int nh = 0;
    for (int i = 0; i < PSIL_NPROBE; i++) {
      double th = 2.0 * M_PI * ((double)i + 0.5) / PSIL_NPROBE;
      double d[3] = {cos(th), 0.0, sin(th)};
      double t;
      if (psil_trace(&grid, &m, eye, d, 0.0, far, NULL, 0, 0, &t)) hd[nh++] = t;
    }
    qsort(hd, (size_t)nh, sizeof *hd, cmp_dbl);
    printf("== КОНТРОЛЬ ЛУЧА (П1, сверка с §110, снятым ПОЛУКУБОМ): доля неба %.4f (§110: "
           "0.1501); горизонт p10 %.2f м, p50 %.2f м (§110: 5.45 / 8.10); в небо ушло %d из "
           "%d горизонтальных\n",
           (double)nsky / PSIL_NPROBE, nh > 0 ? hd[(int)(0.10 * (nh - 1))] : -1.0,
           nh > 0 ? hd[(int)(0.50 * (nh - 1))] : -1.0, PSIL_NPROBE - nh, PSIL_NPROBE);
    free(hd);
    fflush(stdout);
  }

  /* ---- ПРИЁМНИКИ -------------------------------------------------------- */
  double *rc = malloc((size_t)nrecv * 3 * sizeof *rc);
  double *rcen = malloc((size_t)nrecv * 3 * sizeof *rcen); /* центр блока (§159) */
  int32_t *rtri = malloc((size_t)nrecv * sizeof *rtri);
  if (rc == NULL || rtri == NULL || rcen == NULL) return 2;
  uint64_t sd = PSIL_SEED;
  int32_t *clk = malloc((size_t)PSIL_CL_NR * PSIL_CL_NC * sizeof *clk);
  if (clk == NULL) return 2;
  if (!clust) {
    double eyec[3] = HZ_CFG_CITY_EYE, eyeh[3] = HZ_CFG_HALL_EYE;
    const double *ey = city ? eyec : eyeh;
    for (int a = 0; a < 3; a++)
      rc[a] = ey[a];
    rtri[0] = -1; /* камера §110 — точка в воздухе, своего треугольника нет */
    for (int r = 1; r < nrecv; r++) {
      int32_t t = (int32_t)(rnd64(&sd) % (uint64_t)m.nt);
      double p[3][3];
      hz_obj_tri(&m, t, p);
      for (int a = 0; a < 3; a++) {
        rc[(size_t)r * 3 + (size_t)a] = (p[0][a] + p[1][a] + p[2][a]) / 3.0;
        rcen[(size_t)r * 3 + (size_t)a] = rc[(size_t)r * 3 + (size_t)a];
      }
      rtri[r] = t;
    }
    for (int a = 0; a < 3; a++)
      rcen[(size_t)a] = rc[(size_t)a];
  } else {
    /* ГРУППЫ (§156): `NC` центров, в каждом `K` треугольников из шара радиуса
     * `R`. Шар набирается ПО СЕТКЕ, а не отбором наугад: при `R = 0.5` м доля
     * годных треугольников есть `1e-8`, и отбор наугад не сошёлся бы никогда. */
    int32_t *stamp = calloc((size_t)m.nt, sizeof *stamp);
    int32_t *cand = malloc((size_t)PSIL_CL_CANDMAX * sizeof *cand);
    if (stamp == NULL || cand == NULL) return 2;
    for (int i = 0; i < PSIL_CL_NR * PSIL_CL_NC; i++)
      clk[i] = PSIL_CL_K;
    int32_t cmark = 0;
    printf("\n== ГРУППЫ ПРИЁМНИКОВ: %d центров × %d приёмников на радиус\n", PSIL_CL_NC, PSIL_CL_K);
    for (int ri = 0; ri < PSIL_CL_NR; ri++) {
      double R = PSIL_CL_R[ri];
      int64_t got = 0, deficit = 0, trunc = 0;
      for (int ci = 0; ci < PSIL_CL_NC; ci++) {
        int32_t tc = (int32_t)(rnd64(&sd) % (uint64_t)m.nt);
        double p[3][3], cen[3];
        hz_obj_tri(&m, tc, p);
        for (int a = 0; a < 3; a++)
          cen[a] = (p[0][a] + p[1][a] + p[2][a]) / 3.0;
        int64_t nc2 = 0;
        if (R < 0.0) {
          /* НЕГАТИВНЫЙ КОНТРОЛЬ: вся сцена. Шар тут не нужен — берутся
           * случайные треугольники, что и есть «радиус = габарит». */
          while (nc2 < PSIL_CL_K)
            cand[nc2++] = (int32_t)(rnd64(&sd) % (uint64_t)m.nt);
        } else {
          cmark++;
          int32_t c0[3] = {0, 0, 0}, c1[3] = {0, 0, 0};
          for (int a = 0; a < 3; a++) {
            int32_t i0 = (int32_t)floor((cen[a] - R - grid.lo[a]) / grid.cs);
            int32_t i1 = (int32_t)floor((cen[a] + R - grid.lo[a]) / grid.cs);
            c0[a] = i0 < 0 ? 0 : (i0 >= grid.nc[a] ? grid.nc[a] - 1 : i0);
            c1[a] = i1 < 0 ? 0 : (i1 >= grid.nc[a] ? grid.nc[a] - 1 : i1);
          }
          for (int32_t z = c0[2]; z <= c1[2] && nc2 < PSIL_CL_CANDMAX; z++)
            for (int32_t y = c0[1]; y <= c1[1] && nc2 < PSIL_CL_CANDMAX; y++)
              for (int32_t x = c0[0]; x <= c1[0] && nc2 < PSIL_CL_CANDMAX; x++) {
                int64_t c =
                    (int64_t)x + (int64_t)grid.nc[0] * ((int64_t)y + (int64_t)grid.nc[1] * z);
                for (int64_t i = grid.start[c]; i < grid.start[c + 1]; i++) {
                  int32_t t = grid.idx[i];
                  if (stamp[t] == cmark) continue;
                  stamp[t] = cmark;
                  double q[3][3], g2[3], dd = 0.0;
                  hz_obj_tri(&m, t, q);
                  for (int a = 0; a < 3; a++) {
                    g2[a] = (q[0][a] + q[1][a] + q[2][a]) / 3.0 - cen[a];
                    dd += g2[a] * g2[a];
                  }
                  if (dd > R * R) continue;
                  if (nc2 >= PSIL_CL_CANDMAX) {
                    trunc++;
                    break;
                  }
                  cand[nc2++] = t;
                }
              }
          if (nc2 < PSIL_CL_K) deficit++;
          got += nc2;
          /* Перемешать первые `K` — выбор без возврата. */
          for (int64_t i = 0; i < PSIL_CL_K && i < nc2; i++) {
            int64_t j = i + (int64_t)(rnd64(&sd) % (uint64_t)(nc2 - i));
            int32_t s = cand[i];
            cand[i] = cand[j];
            cand[j] = s;
          }
        }
        /* РАЗМЕР ГРУППЫ — ТОЛЬКО РАЗЛИЧНЫЕ ТРЕУГОЛЬНИКИ (А353). Дублировать
         * приёмник, когда в шаре их меньше `K`, значит завысить `A`: копия
         * добавляет `|S|` в числитель и ноль в объединение. Лишние места
         * считаются, но в статистику группы НЕ входят. */
        int keff = (int)(nc2 < PSIL_CL_K ? nc2 : PSIL_CL_K);
        if (keff < 1) keff = 1;
        clk[(ri * PSIL_CL_NC) + ci] = keff;
        for (int j = 0; j < PSIL_CL_K; j++) {
          int r = ((ri * PSIL_CL_NC) + ci) * PSIL_CL_K + j;
          int32_t t = (nc2 > 0) ? cand[j < keff ? j : 0] : tc;
          double q[3][3];
          hz_obj_tri(&m, t, q);
          for (int a = 0; a < 3; a++) {
            rc[(size_t)r * 3 + (size_t)a] = (q[0][a] + q[1][a] + q[2][a]) / 3.0;
            rcen[(size_t)r * 3 + (size_t)a] = cen[a];
          }
          rtri[r] = t;
        }
      }
      /* А353: сколько треугольников НАШЛОСЬ в шаре — печатается, а не
       * подразумевается. Если их меньше `K`, совпадение множеств тривиально. */
      if (R >= 0.0)
        printf("   R = %6.2f м: в шаре в среднем %.1f треугольника; групп с недобором до %d: "
               "%lld; переполнений списка: %lld\n",
               R, (double)got / PSIL_CL_NC, PSIL_CL_K, (long long)deficit, (long long)trunc);
      else
        printf("   R = ВСЯ СЦЕНА (негативный контроль): %d случайных треугольника на группу\n",
               PSIL_CL_K);
    }
    free(cand);
    free(stamp);
    fflush(stdout);
  }

  /* ---- СЧЁТ ------------------------------------------------------------- */
  double t2 = now_s();
  int64_t *tot = calloc((size_t)nrecv * PSIL_NBIN, sizeof *tot);
  int64_t *srv = calloc((size_t)nrecv * PSIL_NALPHA * PSIL_NBIN, sizeof *srv);
  int64_t *totnm = calloc((size_t)nrecv, sizeof *totnm);
  int64_t *srvnm = calloc((size_t)nrecv * PSIL_NALPHA, sizeof *srvnm);
  int64_t *nskip = calloc((size_t)nrecv, sizeof *nskip);
  int64_t *vis = calloc((size_t)nrecv * PSIL_NALPHA * PSIL_NBIN, sizeof *vis);
  int64_t *seg = calloc((size_t)nrecv * PSIL_NALPHA * PSIL_NBIN, sizeof *seg);
  int64_t *vseg = calloc((size_t)nrecv * PSIL_NALPHA * PSIL_NBIN, sizeof *vseg);
  int64_t *vmhist = calloc((size_t)nrecv * 4, sizeof *vmhist);
  int64_t *chmax = calloc((size_t)nrecv, sizeof *chmax);
  int64_t *novf = calloc((size_t)nrecv, sizeof *novf);
  int32_t **vset = calloc((size_t)nrecv, sizeof *vset);
  unsigned char **voct = calloc((size_t)nrecv, sizeof *voct);
  int64_t *vsn = calloc((size_t)nrecv, sizeof *vsn);
  int32_t **aset = calloc((size_t)nrecv, sizeof *aset);
  int64_t *asn = calloc((size_t)nrecv, sizeof *asn);
  if (vset == NULL || vsn == NULL || voct == NULL || aset == NULL || asn == NULL) return 2;
  if (tot == NULL || srv == NULL || totnm == NULL || srvnm == NULL || nskip == NULL ||
      vis == NULL || seg == NULL || vseg == NULL || vmhist == NULL || chmax == NULL || novf == NULL)
    return 2;
  double alpha[PSIL_NALPHA];
  for (int q = 0; q < PSIL_NALPHA; q++)
    alpha[q] = PSIL_ALPHA_DEG[q] * M_PI / 180.0;

#pragma omp parallel for schedule(dynamic, 1)
  for (int r = 0; r < nrecv; r++) {
    signed char *sg = malloc((size_t)m.nt);
    if (sg == NULL) continue;
    const double px = rc[(size_t)r * 3 + 0], py = rc[(size_t)r * 3 + 1], pz = rc[(size_t)r * 3 + 2];
    for (int32_t t = 0; t < m.nt; t++) {
      double s = pn[(size_t)t * 4 + 0] * px + pn[(size_t)t * 4 + 1] * py +
                 pn[(size_t)t * 4 + 2] * pz - pn[(size_t)t * 4 + 3];
      sg[t] = (signed char)(s > 0.0);
    }
    /* ПЕРВЫЙ ПРОХОД — только счёт манифолдных силуэтов: массивы записей
     * выделяются под ТОЧНЫЙ размер, поэтому ни `realloc`, ни запаса «на всякий
     * случай» нет (тот же довод, что в `scene_obj.c`). */
    int64_t nsil = 0;
    for (int64_t ie = 0; ie < ne; ie++) {
      if (enf[ie] != 2) continue;
      if (sg[ef[ie * PSIL_MAXF]] != sg[ef[ie * PSIL_MAXF + 1]]) nsil++;
    }
    int32_t *reI = malloc((size_t)(nsil > 0 ? nsil : 1) * sizeof *reI);
    int32_t *reO = malloc((size_t)(nsil > 0 ? nsil : 1) * sizeof *reO);
    double *reD = malloc((size_t)(nsil > 0 ? nsil : 1) * sizeof *reD);
    double *reV = malloc((size_t)(nsil > 0 ? nsil : 1) * sizeof *reV);
    unsigned char *reB = calloc((size_t)(nsil > 0 ? nsil : 1), 1);
    if (reI == NULL || reO == NULL || reD == NULL || reV == NULL || reB == NULL) {
      free(reB);
      free(reV);
      free(reD);
      free(reO);
      free(reI);
      free(sg);
      continue;
    }
    int64_t nrec = 0;
    for (int64_t ie = 0; ie < ne; ie++) {
      int nf = enf[ie];
      if (nf < 2) continue;
      int nseen = nf < PSIL_MAXF ? nf : PSIL_MAXF;
      int front = 0, back = 0;
      for (int q = 0; q < nseen; q++) {
        if (sg[ef[ie * PSIL_MAXF + q]])
          front = 1;
        else
          back = 1;
      }
      if (!front || !back) continue; /* не силуэт: все грани смотрят одинаково */
      const double *va = m.v + 3 * (size_t)ea[ie];
      const double *vb = m.v + 3 * (size_t)eb[ie];
      double ev[3], dv[3];
      for (int a = 0; a < 3; a++) {
        ev[a] = vb[a] - va[a];
        dv[a] = rc[(size_t)r * 3 + (size_t)a] - 0.5 * (va[a] + vb[a]);
      }
      double d = sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]);
      if (d < PSIL_DMIN) {
        nskip[r]++;
        continue;
      }
      /* ВИДИМАЯ ДЛИНА — поперёк направления на приёмник (А335). */
      double cx = ev[1] * dv[2] - ev[2] * dv[1];
      double cy = ev[2] * dv[0] - ev[0] * dv[2];
      double cz = ev[0] * dv[1] - ev[1] * dv[0];
      double vlen = sqrt(cx * cx + cy * cy + cz * cz) / d;
      if (nf > 2) {
        totnm[r]++;
        for (int q = 0; q < PSIL_NALPHA; q++)
          if (vlen >= alpha[q] * d) srvnm[(size_t)r * PSIL_NALPHA + (size_t)q]++;
        continue;
      }
      int k = ilogb(d) - PSIL_KMIN; /* ilogb = floor(log2) точно, без ошибки округления */
      if (k < 0) k = 0;
      if (k >= PSIL_NBIN) k = PSIL_NBIN - 1;
      tot[(size_t)r * PSIL_NBIN + (size_t)k]++;
      for (int q = 0; q < PSIL_NALPHA; q++)
        if (vlen >= alpha[q] * d)
          srv[((size_t)r * PSIL_NALPHA + (size_t)q) * PSIL_NBIN + (size_t)k]++;
      if (nrec < nsil) {
        reI[nrec] = (int32_t)ie;
        reO[nrec] = k;
        reD[nrec] = d;
        reV[nrec] = vlen;
        nrec++;
      }
    }
    /* ЗАСЛОНЕНИЕ. Три точки ребра (А343): критерий — по СЕРЕДИНЕ, распределение
     * `0/3…3/3` печатается, чтобы двоичность была названа числом. */
    if (!noocc) {
      static const double sfrac[3] = {0.25, 0.5, 0.75};
      for (int64_t i = 0; i < nrec; i++) {
        int32_t ie = reI[i];
        const double *va = m.v + 3 * (size_t)ea[ie];
        const double *vb = m.v + 3 * (size_t)eb[ie];
        int32_t skip[3] = {ef[(int64_t)ie * PSIL_MAXF], ef[(int64_t)ie * PSIL_MAXF + 1], rtri[r]};
        unsigned char bits = 0;
        for (int s = 0; s < 3; s++) {
          double o[3], dd[3], len = 0.0;
          for (int a = 0; a < 3; a++) {
            o[a] = va[a] + sfrac[s] * (vb[a] - va[a]);
            dd[a] = rc[(size_t)r * 3 + (size_t)a] - o[a];
            len += dd[a] * dd[a];
          }
          len = sqrt(len);
          if (len < PSIL_DMIN) continue;
          for (int a = 0; a < 3; a++)
            dd[a] /= len;
          double th;
          if (!psil_trace(&grid, &m, o, dd, PSIL_EPSREL * len, len * (1.0 - PSIL_EPSREL), skip, 3,
                          1, &th))
            bits |= (unsigned char)(1u << s);
        }
        reB[i] = bits;
        vmhist[(size_t)r * 4 + (size_t)__builtin_popcount(bits)]++;
      }
    } else {
      for (int64_t i = 0; i < nrec; i++)
        reB[i] = 7;
      vmhist[(size_t)r * 4 + 3] += nrec;
    }
    /* СШИВКА В УЧАСТКИ КОНТУРА — по каждому `α`, отдельно для всех выживших и
     * отдельно для ВИДИМЫХ. */
    psil_scratch sc;
    int64_t hcap = 16;
    while (hcap < 4 * (nrec + 1))
      hcap <<= 1;
    sc.hcap = hcap;
    size_t nn = (size_t)(nrec > 0 ? nrec : 1);
    sc.uf = malloc(nn * sizeof *sc.uf);
    sc.hkey = malloc((size_t)hcap * sizeof *sc.hkey);
    sc.hhead = malloc((size_t)hcap * sizeof *sc.hhead);
    sc.hnxt = malloc(2 * nn * sizeof *sc.hnxt);
    sc.roct = malloc(nn * sizeof *sc.roct);
    sc.cnt = malloc(nn * sizeof *sc.cnt);
    sc.rdist = malloc(nn * sizeof *sc.rdist);
    int32_t *lstI = malloc(nn * sizeof *lstI);
    int32_t *lstO = malloc(nn * sizeof *lstO);
    double *lstD = malloc(nn * sizeof *lstD);
    if (sc.uf != NULL && sc.hkey != NULL && sc.hhead != NULL && sc.hnxt != NULL &&
        sc.roct != NULL && sc.cnt != NULL && sc.rdist != NULL && lstI != NULL && lstO != NULL &&
        lstD != NULL) {
      for (int q = 0; q < PSIL_NALPHA; q++)
        for (int onlyvis = 0; onlyvis < 2; onlyvis++) {
          int64_t n = 0;
          for (int64_t i = 0; i < nrec; i++) {
            if (reV[i] < alpha[q] * reD[i]) continue;
            if (onlyvis && !(reB[i] & 2u)) continue; /* бит 1 — середина ребра */
            lstI[n] = reI[i];
            lstO[n] = reO[i];
            lstD[n] = reD[i];
            n++;
          }
          /* §156: множество ВИДИМЫХ выживших РЁБЕР при основном `α` — то, чью
           * общность между соседями и меряет О52. Рёбра, а не сегменты: индекс
           * ребра глобален, а сегмент у разных приёмников кроится по-разному
           * (А355). */
          /* §162: ВТОРОЙ набор — выжившие силуэтные рёбра ДО фильтра видимости.
           * Заведён потому, что видимое объединение оказалось РВАНЫМ: фильтр
           * режет силуэтную кривую, и огрублять в нём нечего. Силуэт же связен
           * по построению. */
          if (clust && !onlyvis && q == PSIL_Q_MAIN && n > 0) {
            aset[r] = malloc((size_t)n * sizeof **aset);
            if (aset[r] != NULL) {
              memcpy(aset[r], lstI, (size_t)n * sizeof **aset);
              asn[r] = n;
            }
          }
          if (clust && onlyvis && q == PSIL_Q_MAIN && n > 0) {
            vset[r] = malloc((size_t)n * sizeof **vset);
            voct[r] = malloc((size_t)n);
            if (vset[r] != NULL && voct[r] != NULL) {
              memcpy(vset[r], lstI, (size_t)n * sizeof **vset);
              /* ОКТАВА ОТ ЦЕНТРА БЛОКА, а не от приёмника (А362): иначе одно и
               * то же ребро у соседей попадает в разные октавы и занижает
               * `A_k` разметкой, а не геометрией. */
              for (int64_t i = 0; i < n; i++) {
                int32_t ie = lstI[i];
                const double *va = m.v + 3 * (size_t)ea[ie];
                const double *vb = m.v + 3 * (size_t)eb[ie];
                double dc = 0.0;
                for (int a = 0; a < 3; a++) {
                  double w2 = 0.5 * (va[a] + vb[a]) - rcen[(size_t)r * 3 + (size_t)a];
                  dc += w2 * w2;
                }
                int kk = (dc > 0.0) ? ilogb(sqrt(dc)) - PSIL_KMIN : 0;
                if (kk < 0) kk = 0;
                if (kk >= PSIL_NBIN) kk = PSIL_NBIN - 1;
                voct[r][i] = (unsigned char)kk;
              }
              vsn[r] = n;
            }
          }
          int64_t *oc =
              (onlyvis ? vseg : seg) + ((size_t)r * PSIL_NALPHA + (size_t)q) * (size_t)PSIL_NBIN;
          psil_chains(n, lstI, lstO, lstD, ea, eb, m.v, &sc, oc, &chmax[r], &novf[r]);
          if (onlyvis)
            for (int64_t i = 0; i < n; i++)
              vis[((size_t)r * PSIL_NALPHA + (size_t)q) * PSIL_NBIN + (size_t)lstO[i]]++;
        }
    }
    free(lstD);
    free(lstO);
    free(lstI);
    free(sc.rdist);
    free(sc.cnt);
    free(sc.roct);
    free(sc.hnxt);
    free(sc.hhead);
    free(sc.hkey);
    free(sc.uf);
    free(reB);
    free(reV);
    free(reD);
    free(reO);
    free(reI);
    free(sg);
  }
  double t_cnt = now_s() - t2;
  printf("== СЧЁТ (%.1f с; предсказано А337: смежность 1…3 мин, счёт 2…10 мин)\n", t_cnt);

  /* ---- ДОКЛАД ----------------------------------------------------------- */
  int64_t sil0 = 0;
  for (int k = 0; k < PSIL_NBIN; k++)
    sil0 += tot[k];
  int64_t *silsum = malloc((size_t)nrecv * sizeof *silsum);
  if (silsum == NULL) return 2;
  for (int r = 0; r < nrecv; r++) {
    int64_t s = 0;
    for (int k = 0; k < PSIL_NBIN; k++)
      s += tot[(size_t)r * PSIL_NBIN + (size_t)k];
    silsum[r] = s;
  }
  int64_t *tmp = malloc((size_t)nrecv * sizeof *tmp);
  if (tmp == NULL) return 2;
  memcpy(tmp, silsum, (size_t)nrecv * sizeof *tmp);
  qsort(tmp, (size_t)nrecv, sizeof *tmp, cmp_i64);
  printf("\n== ДВА ЧИСЛА ПОРОЗНЬ (А336):\n");
  printf("   ГРАНИЧНЫХ рёбер в сцене: %lld (одно на всю сцену, от точки не зависит)\n",
         (long long)nb1);
  printf("   СИЛУЭТНЫХ рёбер (манифолдных) у камеры §110: %lld; по %d приёмникам "
         "медиана %lld, p90 %lld, max %lld\n",
         (long long)sil0, nrecv, (long long)pct_i64(tmp, nrecv, 0.5),
         (long long)pct_i64(tmp, nrecv, 0.9), (long long)tmp[nrecv - 1]);
  printf("   отношение силуэтных (медиана) к граничным: %.3f\n",
         (double)pct_i64(tmp, nrecv, 0.5) / (double)nb1);
  {
    int64_t sk = 0, nm = 0;
    for (int r = 0; r < nrecv; r++) {
      sk += nskip[r];
      nm += totnm[r];
    }
    printf("   НЕМАНИФОЛДНЫХ смешанных (>2 граней, в таблицу НЕ входят): всего по приёмникам "
           "%lld, в среднем %.0f на приёмник\n",
           (long long)nm, (double)nm / (double)nrecv);
    for (int q = 0; q < PSIL_NALPHA; q++) {
      int64_t s = 0;
      for (int r = 0; r < nrecv; r++)
        s += srvnm[(size_t)r * PSIL_NALPHA + (size_t)q];
      printf("      из них выживает при α = %.2f°: %.0f на приёмник\n", PSIL_ALPHA_DEG[q],
             (double)s / (double)nrecv);
    }
    printf("   пропущено рёбер (приёмник на ребре, d < %.0e): %lld\n", PSIL_DMIN, (long long)sk);
  }
  if (!noocc) {
    /* А343: двоичность видимости названа ЧИСЛОМ. Доля рёбер, у которых видны
     * не все три точки, и есть мера того, чего стоит выбор «по середине». */
    int64_t h[4] = {0, 0, 0, 0};
    for (int r = 0; r < nrecv; r++)
      for (int j = 0; j < 4; j++)
        h[j] += vmhist[(size_t)r * 4 + (size_t)j];
    printf("   ЧАСТИЧНАЯ ВИДИМОСТЬ ребра (А343), силуэтов по числу видимых точек из трёх: "
           "0/3 %lld, 1/3 %lld, 2/3 %lld, 3/3 %lld\n",
           (long long)h[0], (long long)h[1], (long long)h[2], (long long)h[3]);
    printf("      среди рёбер, видимых ХОТЬ ЧЕМ-ТО, частично видимы %.2f %% — вот цена выбора "
           "«по середине»\n",
           (h[1] + h[2] + h[3]) > 0 ? 100.0 * (double)(h[1] + h[2]) / (double)(h[1] + h[2] + h[3])
                                    : 0.0);
  }
  {
    int64_t cm = 0, ov = 0;
    for (int r = 0; r < nrecv; r++) {
      if (chmax[r] > cm) cm = chmax[r];
      ov += novf[r];
    }
    printf("   СШИВКА (А346): самая длинная цепочка %lld рёбер; вершин с более чем %d "
           "инцидентными силуэтами (разобрано частично) %lld\n",
           (long long)cm, PSIL_VDEG_MAX, (long long)ov);
  }

  for (int q = 0; q < PSIL_NALPHA; q++) {
    printf("\n== α = %.2f°%s\n", PSIL_ALPHA_DEG[q],
           q == PSIL_Q_CTL ? "  (НЕГАТИВНЫЙ КОНТРОЛЬ: отбора нет вовсе)" : "");
    printf("   окт |  расст., м    | силуэтных | выживших | видимых | сегм. | ВИДИМЫХ СЕГМ. |\n");
    double smin = 1e300, smax = 0.0;
    double pmin = 1e300, pmax = 0.0;
    double gmin = 1e300, gmax = 0.0;
    int nused = 0, nzero_in = 0, npred = 0, ngused = 0, ngzero = 0;
    for (int k = 0; k < PSIL_NBIN; k++) {
      int64_t st = 0, ss = 0, sv = 0, sg2 = 0, sw = 0;
      for (int r = 0; r < nrecv; r++) {
        size_t o = ((size_t)r * PSIL_NALPHA + (size_t)q) * PSIL_NBIN + (size_t)k;
        st += tot[(size_t)r * PSIL_NBIN + (size_t)k];
        ss += srv[o];
        sv += vis[o];
        sg2 += seg[o];
        sw += vseg[o];
      }
      if (st == 0) continue;
      double mt = (double)st / (double)nrecv, ms = (double)ss / (double)nrecv;
      double mv = (double)sv / (double)nrecv, mg = (double)sg2 / (double)nrecv;
      double mw = (double)sw / (double)nrecv;
      int kk = k + PSIL_KMIN;
      int inside = (kk <= koct_max);
      printf("   %3d | %6.1f…%-7.1f | %9.1f | %8.1f | %7.1f | %5.1f | %13.1f |%s\n", kk,
             ldexp(1.0, kk), ldexp(1.0, kk + 1), mt, ms, mv, mg, mw, inside ? "" : " ВНЕ ГАБАРИТА");
      if (!inside) continue;
      ngused++;
      if (mw <= 0.0) ngzero++;
      if (mw < gmin) gmin = mw;
      if (mw > gmax) gmax = mw;
      nused++;
      if (ms <= 0.0) nzero_in++;
      if (ms < smin) smin = ms;
      if (ms > smax) smax = ms;
      /* Тот же разброс, но по ОКТАВАМ, КОТОРЫЕ НАЗВАЛО ПРЕДСКАЗАНИЕ §149
       * («от 8 м до 256 м»). Печатается ВТОРЫМ числом, а не вместо первого:
       * приёмка сформулирована по всем октавам внутри габарита, и подменять её
       * задним числом нельзя. Ближние октавы у любой сцены обеднены не
       * критерием, а тем, что вплотную к приёмнику поверхности мало, — это тот
       * же класс ложного нуля, что А334, только с ближнего конца. */
      if (kk >= PSIL_KPRED_LO && kk <= PSIL_KPRED_HI) {
        npred++;
        if (ms < pmin) pmin = ms;
        if (ms > pmax) pmax = ms;
      }
    }
    int64_t *ss = malloc((size_t)nrecv * sizeof *ss);
    if (ss == NULL) return 2;
    for (int r = 0; r < nrecv; r++) {
      int64_t s = 0;
      for (int k = 0; k < PSIL_NBIN; k++)
        s += srv[((size_t)r * PSIL_NALPHA + (size_t)q) * PSIL_NBIN + (size_t)k];
      ss[r] = s;
    }
    qsort(ss, (size_t)nrecv, sizeof *ss, cmp_i64);
    double spread = (smin > 0.0) ? smax / smin : HUGE_VAL;
    double gspread = (gmin > 0.0) ? gmax / gmin : HUGE_VAL;
    printf("   СУММА выживших на приёмник: min %lld, медиана %lld, p90 %lld, max %lld\n",
           (long long)ss[0], (long long)pct_i64(ss, nrecv, 0.5), (long long)pct_i64(ss, nrecv, 0.9),
           (long long)ss[nrecv - 1]);
    /* Три ЧИСЛА ЦЕПОЧКОЙ: рёбра → сегменты → видимые сегменты. Критерий назван
     * планом §153 ДО прогона — ВИДИМЫЕ СЕГМЕНТЫ, и подменять его нельзя (А347). */
    {
      int64_t *sg3 = malloc((size_t)nrecv * sizeof *sg3);
      int64_t *sv3 = malloc((size_t)nrecv * sizeof *sv3);
      int64_t *sw3 = malloc((size_t)nrecv * sizeof *sw3);
      if (sg3 == NULL || sv3 == NULL || sw3 == NULL) return 2;
      for (int r = 0; r < nrecv; r++) {
        int64_t a1 = 0, a2 = 0, a3 = 0;
        for (int k = 0; k < PSIL_NBIN; k++) {
          size_t o = ((size_t)r * PSIL_NALPHA + (size_t)q) * PSIL_NBIN + (size_t)k;
          a1 += seg[o];
          a2 += vis[o];
          a3 += vseg[o];
        }
        sg3[r] = a1;
        sv3[r] = a2;
        sw3[r] = a3;
      }
      qsort(sg3, (size_t)nrecv, sizeof *sg3, cmp_i64);
      qsort(sv3, (size_t)nrecv, sizeof *sv3, cmp_i64);
      qsort(sw3, (size_t)nrecv, sizeof *sw3, cmp_i64);
      double medsrv = (double)pct_i64(ss, nrecv, 0.5);
      printf("   СУММА сегментов (сшивка коллинеарных): медиана %lld -> множитель сшивки %.2f\n",
             (long long)pct_i64(sg3, nrecv, 0.5),
             pct_i64(sg3, nrecv, 0.5) > 0 ? medsrv / (double)pct_i64(sg3, nrecv, 0.5) : 0.0);
      printf("   СУММА видимых рёбер: медиана %lld -> доля видимых %.1f %%\n",
             (long long)pct_i64(sv3, nrecv, 0.5),
             medsrv > 0.0 ? 100.0 * (double)pct_i64(sv3, nrecv, 0.5) / medsrv : 0.0);
      int nzr = 0;
      for (int r = 0; r < nrecv; r++)
        if (sw3[r] == 0) nzr++;
      printf("   *** СУММА ВИДИМЫХ СЕГМЕНТОВ (КРИТЕРИЙ §153): min %lld, медиана %lld, p90 %lld, "
             "max %lld; разброс по октавам %.2f (%d октав, пустых %d)\n",
             (long long)sw3[0], (long long)pct_i64(sw3, nrecv, 0.5),
             (long long)pct_i64(sw3, nrecv, 0.9), (long long)sw3[nrecv - 1], gspread, ngused,
             ngzero);
      /* Приёмник, замурованный в толще блоков, видит НОЛЬ, и такие тянут медиану
       * вниз без всякого отношения к механизму. Их число печатается.
       * СРЕДНЕЕ — не украшение: бюджет кадра есть СУММА по приёмникам, поэтому
       * А333 сравнивается со средним, а медиана говорит про типичный приёмник. */
      double mean = 0.0;
      for (int r = 0; r < nrecv; r++)
        mean += (double)sw3[r];
      mean /= (double)nrecv;
      printf("       приёмников с НУЛЁМ видимых сегментов: %d из %d; СРЕДНЕЕ %.1f; медиана среди "
             "видящих %lld\n",
             nzr, nrecv, mean, (nrecv - nzr) > 0 ? (long long)sw3[nzr + (nrecv - nzr) / 2] : 0LL);
      if (q == PSIL_Q_CTL)
        printf("   НК2: `α = 0` с заслонением обязан остаться большим (> 3e4). Медиана видимых "
               "сегментов %lld — %s\n",
               (long long)pct_i64(sw3, nrecv, 0.5),
               (double)pct_i64(sw3, nrecv, 0.5) > 3.0e4 ? "ОСТАЛСЯ" : "СЕЛ, отбирает НЕ α");
      else
        printf("   ПРИЁМКА §153: медиана ≤ %.0e — %s; разброс ≤ %.0f — %s%s\n", PSIL_ACC_SUM,
               (double)pct_i64(sw3, nrecv, 0.5) <= PSIL_ACC_SUM ? "ДА" : "НЕТ", PSIL_ACC_SPREAD,
               gspread <= PSIL_ACC_SPREAD ? "ДА" : "НЕТ",
               (q == PSIL_Q_MAIN && (double)pct_i64(sw3, nrecv, 0.5) > 1.0e4)
                   ? "; УБИВАЕТ §153: медиана выше 1e4"
                   : "");
      free(sw3);
      free(sv3);
      free(sg3);
    }
    printf("   разброс по октавам ВНУТРИ габарита (%d октав, пустых %d): %.2f\n", nused, nzero_in,
           spread);
    printf("   разброс по октавам ПРЕДСКАЗАНИЯ §149 (%d…%d м, %d октав): %.2f\n",
           1 << PSIL_KPRED_LO, 1 << (PSIL_KPRED_HI + 1), npred,
           (npred > 0 && pmin > 0.0) ? pmax / pmin : HUGE_VAL);
    /* Отношение дальней октавы к ближней — то, чем проверяется рост как d²
     * в негативном контроле (§149). Берутся крайние НЕПУСТЫЕ октавы внутри
     * габарита. */
    {
      int kf = -1, kl = -1;
      for (int k = 0; k < PSIL_NBIN; k++) {
        if (k + PSIL_KMIN > koct_max) break;
        int64_t s = 0;
        for (int r = 0; r < nrecv; r++)
          s += srv[((size_t)r * PSIL_NALPHA + (size_t)q) * PSIL_NBIN + (size_t)k];
        if (s > 0) {
          if (kf < 0) kf = k;
          kl = k;
        }
      }
      if (kf >= 0 && kl > kf) {
        int64_t sf = 0, sl = 0;
        for (int r = 0; r < nrecv; r++) {
          sf += srv[((size_t)r * PSIL_NALPHA + (size_t)q) * PSIL_NBIN + (size_t)kf];
          sl += srv[((size_t)r * PSIL_NALPHA + (size_t)q) * PSIL_NBIN + (size_t)kl];
        }
        printf("   дальняя октава (k=%d) / ближняя (k=%d): %.2f\n", kl + PSIL_KMIN, kf + PSIL_KMIN,
               (double)sl / (double)sf);
      }
    }
    double med = (double)pct_i64(ss, nrecv, 0.5);
    if (q == PSIL_Q_CTL) {
      printf("   КОНТРОЛЬ: обязан ВЗОРВАТЬСЯ (сумма > %.0e, отношение дальней к ближней ≥ 100). "
             "Сумма (медиана) %.0f — %s\n",
             PSIL_KILL_SUM, med, med > PSIL_KILL_SUM ? "ВЗОРВАЛСЯ по сумме" : "НЕ взорвался");
    } else {
      printf("   ПРИЁМКА (А333): сумма ≤ %.0e — %s; разброс ≤ %.0f — %s\n", PSIL_ACC_SUM,
             med <= PSIL_ACC_SUM ? "ДА" : "НЕТ", PSIL_ACC_SPREAD,
             spread <= PSIL_ACC_SPREAD ? "ДА" : "НЕТ");
      if (q == PSIL_Q_MAIN)
        printf("   УБИВАЕТ (§149): сумма > %.0e — %s; дальние выше ближних впятеро — по "
               "разбросу %.2f\n",
               PSIL_KILL_SUM, med > PSIL_KILL_SUM ? "ДА" : "нет", spread);
    }
    free(ss);
  }

  if (clust) {
    /* ---- ДЕЛИТСЯ ЛИ ВИДИМЫЙ КОНТУР МЕЖДУ СОСЕДЯМИ (§156) ------------------ */
    printf("\n== ОБЩНОСТЬ ВИДИМОГО КОНТУРА ПО ГРУППАМ (α = %.2f°, множества ВИДИМЫХ РЁБЕР)\n",
           PSIL_ALPHA_DEG[PSIL_Q_MAIN]);
    printf("   радиус | K факт | слепых |  |S| средн. | объединение |   A   | Жаккар ср./худш. | "
           "пустых пар СРЕДИ ЗРЯЧИХ | доля личного |\n");
    int64_t keycap = (int64_t)PSIL_CL_K * 100000;
    int64_t *key = malloc((size_t)keycap * sizeof *key);
    const int32_t **SS = malloc((size_t)PSIL_CL_K * sizeof *SS);
    int64_t *SN = malloc((size_t)PSIL_CL_K * sizeof *SN);
    int32_t *sub = malloc((size_t)keycap * sizeof *sub);
    if (key == NULL || SS == NULL || SN == NULL || sub == NULL) return 2;
    for (int ri = 0; ri < PSIL_CL_NR; ri++) {
      double sA = 0.0, sU = 0.0, sS = 0.0, sJ = 0.0, sJw = 0.0, sZ = 0.0, sOwn = 0.0;
      double sK = 0.0, sB = 0.0;
      int ncl = 0;
      for (int ci = 0; ci < PSIL_CL_NC; ci++) {
        int base = ((ri * PSIL_CL_NC) + ci) * PSIL_CL_K;
        int keff = clk[(ri * PSIL_CL_NC) + ci];
        for (int j = 0; j < keff; j++) {
          SS[j] = vset[base + j];
          SN[j] = vsn[base + j];
        }
        psil_stat st;
        psil_group_stat(SS, SN, keff, key, keycap, &st);
        if (st.sumS <= 0.0) continue;
        sA += st.sumS / st.uni;
        sU += st.uni;
        sS += st.sumS / (double)keff;
        sJ += st.J;
        sJw += st.Jw;
        sZ += st.zero;
        sOwn += st.own;
        sK += keff;
        sB += (double)st.nblind;
        ncl++;
      }
      if (ncl == 0) continue;
      double R = PSIL_CL_R[ri];
      char rn[32];
      if (R >= 0.0)
        snprintf(rn, sizeof rn, "%6.2f м", R);
      else
        snprintf(rn, sizeof rn, "  СЦЕНА");
      printf("   %s | %6.1f | %6.1f | %11.1f | %11.1f | %5.2f | %6.3f / %.3f | %19.1f %% | "
             "%11.1f %% |%s\n",
             rn, sK / ncl, sB / ncl, sS / ncl, sU / ncl, sA / ncl, sJ / ncl, sJw / ncl,
             100.0 * sZ / ncl, 100.0 * sOwn / ncl, (R < 0.0) ? " НЕГ. КОНТРОЛЬ" : "");
    }

    /* ---- РАЗЛОЖЕНИЕ ПО ОКТАВАМ ОТ ЦЕНТРА БЛОКА (§159) --------------------- */
    for (int ri = 0; ri < PSIL_CL_NR; ri++) {
      double R = PSIL_CL_R[ri];
      /* Печатаются рабочий радиус §158 и негативный контроль — остальные
       * радиусы в таблице выше уже сказали всё, что могли. */
      if (!(R < 0.0) && !(R > 1.9 && R < 2.1)) continue;
      printf("\n== ПО ОКТАВАМ ОТ ЦЕНТРА БЛОКА, %s\n",
             (R < 0.0) ? "группа ПО ВСЕЙ СЦЕНЕ (НЕГАТИВНЫЙ КОНТРОЛЬ)" : "радиус 2 м");
      printf("   окт |  расст., м    | |S_k| на приёмник | объединение_k |  A_k  | Жаккар_k |\n");
      double nearsum = 0.0, farsum = 0.0;
      for (int k = 0; k < PSIL_NBIN; k++) {
        double sA = 0.0, sU = 0.0, sS = 0.0, sJ = 0.0;
        int ncl = 0;
        for (int ci = 0; ci < PSIL_CL_NC; ci++) {
          int base = ((ri * PSIL_CL_NC) + ci) * PSIL_CL_K;
          int keff = clk[(ri * PSIL_CL_NC) + ci];
          int64_t used = 0;
          for (int j = 0; j < keff; j++) {
            int r = base + j;
            SS[j] = sub + used;
            int64_t c = 0;
            for (int64_t i = 0; i < vsn[r] && used + c < keycap; i++)
              if (voct[r][i] == k) sub[used + c++] = vset[r][i];
            SN[j] = c;
            used += c;
          }
          psil_stat st;
          psil_group_stat(SS, SN, keff, key, keycap, &st);
          if (st.sumS <= 0.0) continue;
          sA += st.sumS / st.uni;
          sU += st.uni;
          sS += st.sumS / (double)keff;
          sJ += st.J;
          ncl++;
        }
        if (ncl == 0) continue;
        int kk = k + PSIL_KMIN;
        printf("   %3d | %6.1f…%-7.1f | %17.1f | %13.1f | %5.2f | %8.3f |\n", kk, ldexp(1.0, kk),
               ldexp(1.0, kk + 1), sS / ncl, sU / ncl, sA / ncl, sJ / ncl);
        /* Граница «ближе восьми метров» — не выбранная, а предъявленная
         * таблицей (А364): перелом кривой Жаккара виден по строкам. */
        if (kk < 3)
          nearsum += sS / ncl;
        else
          farsum += sS / ncl;
      }
      printf("   ЛИЧНАЯ ЧАСТЬ (ближе 8 м): %.1f рёбер на приёмник; ОБЩАЯ (дальше): %.1f\n", nearsum,
             farsum);
    }

    /* ---- УГЛОВОЕ ОГРУБЛЕНИЕ КОНТУРА (§162, шаг О54) ----------------------- */
    for (int ri = 0; ri < PSIL_CL_NR; ri++) {
      double R = PSIL_CL_R[ri];
      if (!(R > 1.9 && R < 2.1)) continue; /* рабочий радиус §158 */
      for (int which = 0; which < 2; which++) {
        /* which = 0 — ВИДИМОЕ объединение (то, что фильтрует приёмник по §158);
         * which = 1 — ВЫЖИВШИЕ силуэты ДО фильтра видимости: кривая там связна по
         * построению, и огрублять есть что. Чистота хорды в ОБОИХ случаях
         * считается против ВИДИМОГО множества приёмника — вопрос один и тот же:
         * решается ли хорда одной проверкой. */
        int32_t **SET = which ? aset : vset;
        int64_t *SETN = which ? asn : vsn;
        printf("\n== УГЛОВОЕ ОГРУБЛЕНИЕ %s (радиус %.1f м, α = %.2f°)\n",
               which == 0 ? "ВИДИМОГО объединения блока" : "СИЛУЭТА блока (ДО видимости)", R,
               PSIL_ALPHA_DEG[PSIL_Q_MAIN]);
        double sU = 0.0, sCh[3] = {0.0, 0.0, 0.0}, sChain = 0.0, sColl = 0.0, sPart = 0.0;
        double sPure = 0.0, sMix = 0.0, sVis = 0.0;
        double octE[PSIL_NBIN], octC[PSIL_NBIN];
        for (int k = 0; k < PSIL_NBIN; k++) {
          octE[k] = 0.0;
          octC[k] = 0.0;
        }
        int ncl = 0;
        for (int ci = 0; ci < PSIL_CL_NC; ci++) {
          int base = ((ri * PSIL_CL_NC) + ci) * PSIL_CL_K;
          int keff = clk[(ri * PSIL_CL_NC) + ci];
          /* объединение блока, отсортированное и без повторов */
          int64_t nk = 0;
          for (int j = 0; j < keff; j++)
            for (int64_t i = 0; i < SETN[base + j] && nk < keycap; i++)
              key[nk++] = SET[base + j][i];
          if (nk == 0) continue;
          qsort(key, (size_t)nk, sizeof *key, cmp_i64);
          int64_t nu = 0;
          for (int64_t i = 0; i < nk; i++)
            if (i == 0 || key[i] != key[i - 1]) sub[nu++] = (int32_t)key[i];
          const double *cen = rcen + (size_t)base * 3;
          int64_t hcap2 = 16;
          while (hcap2 < 4 * (nu + 1))
            hcap2 <<= 1;
          psil_contour C;
          C.ord = malloc((size_t)nu * sizeof *C.ord);
          C.coff = malloc((size_t)(nu + 2) * sizeof *C.coff);
          C.voct = malloc((size_t)(nu + 2) * sizeof *C.voct);
          int32_t *hk = malloc((size_t)hcap2 * sizeof *hk);
          int32_t *hh = malloc((size_t)hcap2 * sizeof *hh);
          int32_t *hn = malloc((size_t)(2 * nu) * sizeof *hn);
          int32_t *vq = malloc((size_t)(nu + 2) * sizeof *vq);
          unsigned char *kp = malloc((size_t)(2 * nu + 4));
          int32_t *stk = malloc((size_t)(4 * nu + 16) * sizeof *stk);
          if (C.ord == NULL || C.coff == NULL || C.voct == NULL || hk == NULL || hh == NULL ||
              hn == NULL || vq == NULL || kp == NULL || stk == NULL) {
            free(stk);
            free(kp);
            free(vq);
            free(hn);
            free(hh);
            free(hk);
            free(C.voct);
            free(C.coff);
            free(C.ord);
            return 2;
          }
          /* Три допуска: 0 (регрессия), α (рабочий), ∞ (потолок связности, А372). */
          double tolv[3] = {0.0, PSIL_ALPHA_DEG[PSIL_Q_MAIN] * M_PI / 180.0, 1.0e9};
          for (int tv = 0; tv < 3; tv++) {
            psil_contour_build(&C, sub, nu, ea, eb, m.v, cen, R, tolv[tv], hk, hh, hn, hcap2, vq,
                               kp, stk);
            sCh[tv] += (double)C.nchord;
            if (tv == 1) {
              sChain += (double)C.nchain;
              sColl += (double)C.ncollapse;
              sPart += C.nwalk;
              for (int64_t c = 0; c < C.nchord; c++)
                octC[C.voct[c]] += 1.0;
              /* ЧИСТОТА ХОРДЫ: у приёмника все её рёбра видимы, все невидимы или
               * смешанно. Без этого «одна проверка на хорду» есть допущение. */
              for (int j = 0; j < keff; j++) {
                if (vsn[base + j] == 0) continue;
                for (int64_t i = 0; i < nu; i++)
                  kp[i] = 0;
                int64_t ia = 0, ib = 0;
                while (ia < nu && ib < vsn[base + j]) {
                  if (sub[ia] == vset[base + j][ib]) {
                    kp[ia] = 1;
                    ia++;
                    ib++;
                  } else if (sub[ia] < vset[base + j][ib])
                    ia++;
                  else
                    ib++;
                }
                for (int64_t c = 0; c < C.nchord; c++) {
                  int64_t f = 0, tt = 0;
                  for (int32_t p = C.coff[c]; p < C.coff[c + 1]; p++) {
                    /* позиция ребра `C.ord[p]` в `sub` — бинарным поиском */
                    int64_t lo2 = 0, hi2 = nu - 1, pos = -1;
                    while (lo2 <= hi2) {
                      int64_t mid = (lo2 + hi2) / 2;
                      if (sub[mid] == C.ord[p]) {
                        pos = mid;
                        break;
                      }
                      if (sub[mid] < C.ord[p])
                        lo2 = mid + 1;
                      else
                        hi2 = mid - 1;
                    }
                    tt++;
                    if (pos >= 0 && kp[pos]) f++;
                  }
                  if (f == 0) continue;
                  if (f == tt)
                    sPure += 1.0;
                  else
                    sMix += 1.0;
                  sVis += 1.0;
                }
              }
            }
          }
          for (int64_t i = 0; i < nu; i++) {
            const double *va = m.v + 3 * (size_t)ea[sub[i]];
            const double *vb = m.v + 3 * (size_t)eb[sub[i]];
            double dc = 0.0;
            for (int a = 0; a < 3; a++) {
              double w2 = 0.5 * (va[a] + vb[a]) - cen[a];
              dc += w2 * w2;
            }
            int kk = (dc > 0.0) ? ilogb(sqrt(dc)) - PSIL_KMIN : 0;
            if (kk < 0) kk = 0;
            if (kk >= PSIL_NBIN) kk = PSIL_NBIN - 1;
            octE[kk] += 1.0;
          }
          sU += (double)nu;
          ncl++;
          free(stk);
          free(kp);
          free(vq);
          free(hn);
          free(hh);
          free(hk);
          free(C.voct);
          free(C.coff);
          free(C.ord);
        }
        if (ncl == 0) continue;
        printf("   объединение блока: %.1f рёбер -> %.1f хорд при допуске α·d "
               "(множитель %.2f)\n",
               sU / ncl, sCh[1] / ncl, (sCh[1] > 0.0) ? sU / sCh[1] : 0.0);
        printf("   допуск 0 (регрессия): %.1f хорд; допуск ∞ (ПОТОЛОК связности, А372): %.1f "
               "цепей; схлопнулось в точку: %.1f\n",
               sCh[0] / ncl, sCh[2] / ncl, sColl / ncl);
        printf("   разбиение (А-проверка): пройдено %.1f рёбер против |U| = %.1f — %s\n",
               sPart / ncl, sU / ncl, (fabs(sPart - sU) < 1e-9) ? "СОВПАЛО" : "РАЗОШЛОСЬ");
        printf(
            "   ЧИСТОТА (А343 для хорд): видимых приёмнику хорд %.0f, из них СМЕШАННЫХ %.1f %%\n",
            sVis / ncl, (sVis > 0.0) ? 100.0 * sMix / sVis : 0.0);
        printf("   окт | рёбер_k | хорд_k | множитель |\n");
        for (int k = 0; k < PSIL_NBIN; k++) {
          if (octE[k] <= 0.0) continue;
          printf("   %3d | %7.1f | %6.1f | %9.2f |\n", k + PSIL_KMIN, octE[k] / ncl, octC[k] / ncl,
                 (octC[k] > 0.0) ? octE[k] / octC[k] : 0.0);
        }
        double chords = sCh[1] / ncl;
        printf("   АРИФМЕТИКА КАДРА: 446013 приёмников × %.1f хорд = %.2e операций против бюджета "
               "2.6e8 — приёмка §162 (≤ 580): %s\n",
               chords, 446013.0 * chords, (chords <= 580.0) ? "ДА" : "НЕТ");
      }
    }
    /* Регрессия К=1 (§156): объединение группы из ОДНОГО приёмника обязано
     * совпасть с его множеством, а `A` — быть равным единице ТОЧНО. */
    {
      int r = 0;
      while (r < nrecv && vsn[r] == 0)
        r++;
      if (r < nrecv) {
        int64_t uni = 0;
        for (int64_t i = 0; i < vsn[r]; i++)
          if (i == 0 || vset[r][i] != vset[r][i - 1]) uni++;
        printf("   регрессия K=1: |S| = %lld, |объединение| = %lld, A = %.6f — %s\n",
               (long long)vsn[r], (long long)uni, (double)vsn[r] / (double)uni,
               (uni == vsn[r]) ? "СОВПАЛО" : "РАЗОШЛОСЬ");
      }
    }
    free(sub);
    free(SN);
    free(SS);
    free(key);
  }
  for (int r = 0; r < nrecv; r++) {
    free(vset[r]);
    free(voct[r]);
    free(aset[r]);
  }
  free(aset);
  free(asn);
  free(voct);
  free(vset);
  free(vsn);
  free(rcen);
  free(clk);
  free(tmp);
  free(silsum);
  psil_grid_free(&grid);
  free(rtri);
  free(novf);
  free(chmax);
  free(vmhist);
  free(vseg);
  free(seg);
  free(vis);
  free(nskip);
  free(srvnm);
  free(totnm);
  free(srv);
  free(tot);
  free(rc);
  free(pn);
  free(enf);
  free(ef);
  free(eb);
  free(ea);
  hz_obj_free(&m);
  printf("\n== ВСЕГО %.1f с\n", now_s() - t0);
  return 0;
}
