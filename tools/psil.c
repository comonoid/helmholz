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

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static uint64_t rnd64(uint64_t *s) {
  uint64_t x = *s;
  x ^= x << 13;
  x ^= x >> 7;
  x ^= x << 17;
  *s = x;
  return x;
}

static int cmp_i64(const void *a, const void *b) {
  int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

static int64_t pct_i64(const int64_t *v, int n, double p) {
  if (n <= 0) return 0;
  int k = (int)(p * (double)(n - 1));
  if (k < 0) k = 0;
  if (k >= n) k = n - 1;
  return v[k];
}

int main(int argc, char **argv) {
  int city = 1, nrecv = PSIL_NRECV;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "hall") == 0) city = 0;
    if (strcmp(argv[i], "city") == 0) city = 1;
    if (strncmp(argv[i], "n=", 2) == 0) nrecv = (int)strtol(argv[i] + 2, NULL, 10);
  }
  if (nrecv < 1) nrecv = 1;

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

  /* ---- ПРИЁМНИКИ -------------------------------------------------------- */
  double *rc = malloc((size_t)nrecv * 3 * sizeof *rc);
  if (rc == NULL) return 2;
  {
    double eyec[3] = HZ_CFG_CITY_EYE, eyeh[3] = HZ_CFG_HALL_EYE;
    const double *ey = city ? eyec : eyeh;
    for (int a = 0; a < 3; a++)
      rc[a] = ey[a];
  }
  uint64_t sd = PSIL_SEED;
  for (int r = 1; r < nrecv; r++) {
    int32_t t = (int32_t)(rnd64(&sd) % (uint64_t)m.nt);
    double p[3][3];
    hz_obj_tri(&m, t, p);
    for (int a = 0; a < 3; a++)
      rc[(size_t)r * 3 + (size_t)a] = (p[0][a] + p[1][a] + p[2][a]) / 3.0;
  }

  /* ---- СЧЁТ ------------------------------------------------------------- */
  double t2 = now_s();
  int64_t *tot = calloc((size_t)nrecv * PSIL_NBIN, sizeof *tot);
  int64_t *srv = calloc((size_t)nrecv * PSIL_NALPHA * PSIL_NBIN, sizeof *srv);
  int64_t *totnm = calloc((size_t)nrecv, sizeof *totnm);
  int64_t *srvnm = calloc((size_t)nrecv * PSIL_NALPHA, sizeof *srvnm);
  int64_t *nskip = calloc((size_t)nrecv, sizeof *nskip);
  if (tot == NULL || srv == NULL || totnm == NULL || srvnm == NULL || nskip == NULL) return 2;
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
    }
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

  for (int q = 0; q < PSIL_NALPHA; q++) {
    printf("\n== α = %.2f°%s\n", PSIL_ALPHA_DEG[q],
           q == PSIL_Q_CTL ? "  (НЕГАТИВНЫЙ КОНТРОЛЬ: отбора нет вовсе)" : "");
    printf("   октава k |  расстояние, м |  силуэтных/приёмник |  выживших/приёмник |  доля |\n");
    double smin = 1e300, smax = 0.0;
    double pmin = 1e300, pmax = 0.0;
    int nused = 0, nzero_in = 0, npred = 0;
    for (int k = 0; k < PSIL_NBIN; k++) {
      int64_t st = 0, ss = 0;
      for (int r = 0; r < nrecv; r++) {
        st += tot[(size_t)r * PSIL_NBIN + (size_t)k];
        ss += srv[((size_t)r * PSIL_NALPHA + (size_t)q) * PSIL_NBIN + (size_t)k];
      }
      if (st == 0) continue;
      double mt = (double)st / (double)nrecv, ms = (double)ss / (double)nrecv;
      int kk = k + PSIL_KMIN;
      int inside = (kk <= koct_max);
      printf("   %8d | %6.2f…%-7.2f | %19.1f | %18.1f | %5.1f %% |%s\n", kk, ldexp(1.0, kk),
             ldexp(1.0, kk + 1), mt, ms, 100.0 * ms / mt, inside ? "" : " ВНЕ ГАБАРИТА");
      if (!inside) continue;
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
    printf("   СУММА выживших на приёмник: min %lld, медиана %lld, p90 %lld, max %lld\n",
           (long long)ss[0], (long long)pct_i64(ss, nrecv, 0.5), (long long)pct_i64(ss, nrecv, 0.9),
           (long long)ss[nrecv - 1]);
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

  free(tmp);
  free(silsum);
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
