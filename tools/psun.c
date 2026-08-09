/* psun.c — ШАГ О73 (Ф3): ФРОНТ ОТ СОЛНЦА С СОСТОЯНИЕМ НА ГРАНЯХ.
 * План §266, аудит §267, поправки 08-05 (огрублять, а не отсекать).
 *
 * ЗАПУСК:
 *   psun ФАЙЛ.obj МАСШТАБ [leaf=N] [lev=N] [grade=0|1] [px=9.1]
 *        [sun=dx,dy,dz] [shadow=0] [cam=1] [mark=1|2] [buf=N] [loose=1|2]
 *        [ref=N] [eta=0]
 *   `shadow=0` — НЕГАТИВНЫЙ КОНТРОЛЬ: ни заслонения, ни огрубления по ровности.
 *   Обязан воспроизвести числа О72 ДО ЕДИНИЦЫ; расхождение означает, что
 *   изменилось что-то ещё.
 *   `mark=1` — пометки невидимости буфером дальности (§278, односторонние);
 *   `mark=2` — ПРЕЖНИЙ проход фронтом, приёмку НЕ проходящий (оставлен ради
 *              воспроизводимости §271/§272 и как предмет сравнения);
 *   `buf=N`  — сторона буфера дальности, по умолчанию сторона кадра;
 *   `loose=` — негативные контроли буфера, разбор в `pcull.h`;
 *   `eta=0`  — без эталона лучами по диску: он стоит `8 × 8` лучей на выборочную
 *              ячейку перебором ВСЕХ треугольников и на Сан-Мигеле идёт
 *              десятки минут, тогда как приёмка пометок — один луч на ячейку.
 */
#include "image.h"
#include "pclip.h"
#include "pcull.h"
#include "pfront.h"
#include "pmark.h"
#include "pocc.h"
#include "ptree.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define HZ_REF_RAYS 16 /* лучей на точку выборки: диск солнца мал, шестнадцать дают шаг ~2 % */

static int cmp_dbl_psun(const void *x, const void *y) {
  double a = *(const double *)x, b = *(const double *)y;
  return (a < b) ? -1 : ((a > b) ? 1 : 0);
}

/* РИСОВАНИЕ ПО АЛГОРИТМУ ХУДОЖНИКА, БЕЗ БУФЕРА ГЛУБИНЫ (§353, довод пользователя
 * 08-09). Ячейки фронта ДИЗЪЮНКТНЫ — они разбиение пространства, — поэтому
 * порядок октантов по знакам направления на камеру есть ТОЧНЫЙ порядок
 * видимости, и достаточно рисовать сзади наперёд с простой перезаписью.
 *
 * ЧЕМ ЭТО ЛУЧШЕ ПРЕЖНЕГО. §351 сравнивал ячейки по ближнему углу коробки, а это
 * НЕ порядок видимости: коробка с более близким углом может целиком лежать за
 * другой. Здесь порядок верен по построению, буфер глубины не заводится вовсе
 * (на 1024² это 4 МБ чтения и записи за кадр — чистая полоса памяти), и
 * сравнения на пиксель тоже нет. */
static void paint_cell(const hz_ptree *T, int32_t nid, const float *cellf, const hz_pcull *C,
                       const double *eye, int side, double *out, int64_t *ndrawn) {
  const hz_ptnode *N = &T->nd[nid];
  if (cellf[nid] > -1.5f) {
    /* Экранный след коробки — по восьми углам, той же рамой, что `hz_pcull_ray`. */
    double u0 = 1e300, u1 = -1e300, v0 = 1e300, v1 = -1e300;
    for (int k = 0; k < 8; k++) {
      double W[3];
      for (int c = 0; c < 3; c++)
        W[c] = ((k & (1 << c)) ? N->hi[c] : N->lo[c]) - eye[c];
      double st = W[0] * C->fw[0] + W[1] * C->fw[1] + W[2] * C->fw[2];
      if (!(st > 1e-9)) return; /* ячейка задевает плоскость глаза — не рисуем */
      double ar = (W[0] * C->rt[0] + W[1] * C->rt[1] + W[2] * C->rt[2]) / st;
      double br = (W[0] * C->up[0] + W[1] * C->up[1] + W[2] * C->up[2]) / st;
      double hh = 0.5 * (double)side;
      double su = (ar / C->tanh_ + 1.0) * hh - 0.5;
      double sv = (br / C->tanh_ + 1.0) * hh - 0.5;
      if (su < u0) u0 = su;
      if (su > u1) u1 = su;
      if (sv < v0) v0 = sv;
      if (sv > v1) v1 = sv;
    }
    int i0 = (int)floor(u0), i1 = (int)ceil(u1);
    int j0 = (int)floor(v0), j1 = (int)ceil(v1);
    if (i0 < 0) i0 = 0;
    if (j0 < 0) j0 = 0;
    if (i1 >= side) i1 = side - 1;
    if (j1 >= side) j1 = side - 1;
    if (i1 < i0 || j1 < j0) return;
    (*ndrawn)++;
    double f = (cellf[nid] >= 0.0f) ? (double)cellf[nid] : 0.0;
    for (int jj = j0; jj <= j1; jj++)
      for (int ii = i0; ii <= i1; ii++)
        out[(size_t)jj * (size_t)side + (size_t)ii] = f;
    return;
  }
  if (N->child < 0) return;
  /* Порядок СЗАДИ НАПЕРЁД: дальний октант первым. Дальняя сторона по оси `c` —
   * та, что напротив глаза относительно середины. */
  double mid[3];
  int e[3];
  for (int c = 0; c < 3; c++) {
    mid[c] = 0.5 * (N->lo[c] + N->hi[c]);
    e[c] = (eye[c] < mid[c]) ? 1 : 0;
  }
  for (int far = 3; far >= 0; far--)
    for (int k = 0; k < 8; k++) {
      int nf = 0;
      for (int c = 0; c < 3; c++)
        if (((k >> c) & 1) == e[c]) nf++;
      if (nf != far) continue;
      paint_cell(T, N->child + k, cellf, C, eye, side, out, ndrawn);
    }
}

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "usage: psun ФАЙЛ.obj МАСШТАБ [leaf=N] [lev=N] [grade=0|1] [px=9.1] "
                    "[sun=dx,dy,dz] [shadow=0]\n");
    return 1;
  }
  int leafmax = 0, maxlev = 0, grade = 1, shadow = 1, cam = 0, nref = 0, cliplev = -1, camfull = 0,
      usemark = 0, relax = HZ_PMARK_OUT_LEVELS, refinemax = 0, nk = 1;
  double camo[3] = {0.0, 0.0, 0.0}, camf[3] = {0.0, 0.0, 1.0}, camat[3] = {0.0, 0.0, 1.0};
  int loose = 0, bufside = HZ_CFG_W, eta = 1, mlev = HZ_PMARK_LEVEL, img = 0, nrec = 0, ndbg = 0;
  /* УМОЛЧАНИЯ ПОСТАВЛЕНЫ ПО ЗАМЕРАМ 08-07, А НЕ ПО ИСТОРИИ (§302). Каждое из
   * трёх измерено и каждое улучшает ответ; правило §13 требует, чтобы умолчание
   * ошибалось в сторону ВИДИМОГО отказа, то есть было точным, а не быстрым.
   *   перекрытие ОБЪЕДИНЕНИЕМ (§292): среднее 0.0470 -> 0.0376 на большом
   *     источнике, и это единственная из трёх оценок, правильная по определению;
   *   ГЛУБИНА В ЯЧЕЙКЕ (§296) с допуском ребро/128: 6.91 % -> 0.19 % неверных
   *     пикселей при том же числе ячеек;
   *   ПОЛ ОЧЕРКА (§300): цена огрубления 7.1e-4 -> 1.6e-6.
   * Прежние числа §271…§295 воспроизводятся ключами `cover=0 indep=0 silh=0`. */
  int indep = 1, intol = 128, silh = 1;
  double ooff = 1e-5, nofrec = 0.0;
  int coverset = 0, shiftset = 0, zbuf = 1;
  int g_virt = 0;      /* §338 */
  int cellimg = 0;     /* §349: картинка из ячеек фронта, без второй видимости */
  int g_camfloor = -1; /* §340: −1 — «как virt» */
  /* УГЛОВОЙ РАДИУС ИСТОЧНИКА — ПАРАМЕТР, А НЕ КОНСТАНТА (замечание пользователя
   * 08-07). Он был зашит числом солнца в ДВУХ местах — у фронта и у эталона, — и
   * потому «протяжённый источник» из §275 означал лишь `K` выборок ТОГО ЖЕ
   * крошечного диска: полутень оставалась ýже ячейки, и §289 мерил не то.
   * Настоящая проверка — источник, у которого полутень ШИРЕ ячейки. */
  double asun = 0.5 * 9.3e-3;
  double px = 9.1, sdir[3] = {0.3, -0.9, 0.3};
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "leaf=", 5) == 0) leafmax = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "lev=", 4) == 0) maxlev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "grade=", 6) == 0) grade = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "px=", 3) == 0) px = strtod(argv[i] + 3, NULL);
    if (strncmp(argv[i], "shadow=", 7) == 0) shadow = (int)strtol(argv[i] + 7, NULL, 10);
    if (strncmp(argv[i], "cam=", 4) == 0) cam = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "ref=", 4) == 0) nref = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "clip=", 5) == 0) cliplev = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "camfull=", 8) == 0) camfull = (int)strtol(argv[i] + 8, NULL, 10);
    if (strncmp(argv[i], "mark=", 5) == 0) usemark = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "relax=", 6) == 0) relax = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "loose=", 6) == 0) loose = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "buf=", 4) == 0) bufside = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "eta=", 4) == 0) eta = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "mlev=", 5) == 0) mlev = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "img=", 4) == 0) img = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "rec=", 4) == 0) nrec = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "dbg=", 4) == 0) ndbg = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "indep=", 6) == 0) indep = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "intol=", 6) == 0) intol = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "silh=", 5) == 0) silh = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "zbuf=", 5) == 0) zbuf = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "geo=", 4) == 0) hz_ptree_geo_leaf = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "shift=", 6) == 0) hz_pfront_shift = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "ooff=", 5) == 0) ooff = strtod(argv[i] + 5, NULL);
    if (strncmp(argv[i], "nofrec=", 7) == 0) nofrec = strtod(argv[i] + 7, NULL);
    if (strncmp(argv[i], "flat=", 5) == 0) hz_ptree_flat_split = (int)strtol(argv[i] + 5, NULL, 10);
    /* ВЫБОРОК ПО ДИСКУ ИСТОЧНИКА. Механизм Ф4 (§275) был написан, но ключа не
     * имел, и `K` оставалось единицей — то есть источник точечным, а полутени в
     * постановке не было вовсе (§274). Довод пользователя 08-07: на протяжённом
     * источнике отлаживать проще, потому что поле СГЛАЖИВАЕТСЯ, а линейное
     * состояние на грани (Р1) для гладкого поля и заведено. */
    if (strncmp(argv[i], "nk=", 3) == 0) nk = (int)strtol(argv[i] + 3, NULL, 10);
    if (strncmp(argv[i], "asun=", 5) == 0) asun = strtod(argv[i] + 5, NULL);
    /* §338: ленивый спуск ниже листа; умолчание — прежнее поведение. */
    if (strncmp(argv[i], "cellimg=", 8) == 0) {
      cellimg = (int)strtol(argv[i] + 8, NULL, 10);
      continue;
    }
    if (strncmp(argv[i], "camfloor=", 9) == 0) {
      g_camfloor = (int)strtol(argv[i] + 9, NULL, 10);
      continue;
    }
    if (strncmp(argv[i], "virt=", 5) == 0) {
      g_virt = (int)strtol(argv[i] + 5, NULL, 10);
      continue;
    }
    if (strncmp(argv[i], "cover=", 6) == 0) {
      hz_pfront_cover_sum = (int)strtol(argv[i] + 6, NULL, 10);
      coverset = 1;
    }
    if (strncmp(argv[i], "sun=", 4) == 0) {
      char *e = NULL;
      sdir[0] = strtod(argv[i] + 4, &e);
      if (e != NULL && *e == ',') sdir[1] = strtod(e + 1, &e);
      if (e != NULL && *e == ',') sdir[2] = strtod(e + 1, NULL);
    }
  }
  if (!coverset) hz_pfront_cover_sum = 2; /* объединение — умолчание, см. выше */
  /* ПОЛНЫЙ ПЕРЕНОС — УМОЛЧАНИЕ (§310): средняя ошибка 0.0256 -> 0.0101 при +14 %%
   * ячеек и +47 %% времени. Правило §13: умолчание точное, а не быстрое. */
  if (!shiftset) hz_pfront_shift = 2;
  hz_objmesh m;
  double t0 = now_s();
  if (hz_obj_load(&m, argv[1], strtod(argv[2], NULL)) != 0) {
    fprintf(stderr, "нет сцены %s\n", argv[1]);
    return 1;
  }
  printf("== СЦЕНА %s: треугольников %d, габарит %.2f x %.2f x %.2f м, чтение %.1f с\n", argv[1],
         m.nt, m.hi[0] - m.lo[0], m.hi[1] - m.lo[1], m.hi[2] - m.lo[2], now_s() - t0);
  /* РАЗМЕР ТРЕУГОЛЬНИКА ПРОТИВ РАЗМЕРА ЛИСТА (§318). Избыточность разбиения
   * снизу ограничена геометрией: компактный треугольник размера `t` попадает
   * примерно в `(1 + t/L)³` ячеек ребром `L`. Значит требование по избыточности
   * ЗАДАЁТ размер листа, а тот — число треугольников в нём. Обе величины надо
   * знать в метрах, а не в уровнях. */
  {
    double sm = 0.0, smx = 0.0;
    for (int32_t t2 = 0; t2 < m.nt; t2++) {
      double e = 0.0;
      for (int i = 0; i < 3; i++) {
        const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t2 + (size_t)i];
        const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t2 + (size_t)((i + 1) % 3)];
        double d2 = 0.0;
        for (int c = 0; c < 3; c++)
          d2 += (B[c] - A[c]) * (B[c] - A[c]);
        if (d2 > e) e = d2;
      }
      e = sqrt(e);
      sm += e;
      if (e > smx) smx = e;
    }
    double tm = sm / (double)m.nt;
    double L = m.hi[0] - m.lo[0];
    printf("== ТРЕУГОЛЬНИК: наибольшее ребро в среднем %.4f м, наибольшее в сцене %.3f м; "
           "габарит %.2f м\n",
           tm, smx, L);
    printf("   лист ребром 10×треугольника (%.3f м) — это уровень %.1f, "
           "треугольников в листе ~%.0f\n",
           10.0 * tm, log2(L / (10.0 * tm)), 100.0);
  }

  hz_ptree T;
  t0 = now_s();
  if (hz_ptree_build_cut(&T, &m, leafmax, maxlev, grade) != 0) {
    fprintf(stderr, "отказ дерева\n");
    hz_obj_free(&m);
    return 2;
  }
  printf("== ДЕРЕВО за %.2f с: узлов %d, листьев %lld, глубина %d\n", now_s() - t0, T.nnd,
         (long long)T.nleaf, T.depth);
  /* ИЗБЫТОЧНОСТЬ ДУБЛИРОВАНИЯ — величина, которую построитель считает, но никто
   * не печатал (§317, подозрение пользователя 08-08, что дерево строится не так).
   * `n_leaf_tri` есть сумма кандидатов по листьям, то есть сколько ссылок было бы,
   * если бы их хранили. Отношение к числу треугольников и говорит, во сколько раз
   * геометрия размножена по ячейкам. */
  printf("   КАНДИДАТОВ ПО ЛИСТЬЯМ %lld при %d треугольниках — избыточность %.1f× ; "
         "узлов на треугольник %.1f; добавлено градуировкой %lld\n",
         (long long)T.n_leaf_tri, m.nt, (double)T.n_leaf_tri / (double)m.nt,
         (double)T.nnd / (double)m.nt, (long long)T.ngrade);

  /* ЧЕМ ОСТАНАВЛИВАЕТСЯ ДРОБЛЕНИЕ (§317): числом кандидатов или ПРЕДЕЛОМ ГЛУБИНЫ.
   * Если вторым — дерево дробится не потому, что надо, а потому, что ему не дали
   * дробиться дальше, и избыточность есть следствие, а не причина. */
  {
    int64_t byl[16] = {0};
    int64_t leafs = 0;
    for (int32_t q = 0; q < T.nnd; q++)
      if (T.nd[q].child < 0) {
        int lv = (int)T.lev[q];
        if (lv > 15) lv = 15;
        byl[lv]++;
        leafs++;
      }
    printf("   ЛИСТЬЯ ПО УРОВНЯМ:");
    for (int lv = 0; lv < 16; lv++)
      if (byl[lv] > 0) printf(" %d:%lld", lv, (long long)byl[lv]);
    printf("\n   на ПРЕДЕЛЕ ГЛУБИНЫ (%d) — %lld из %lld (%.1f %%)\n", T.maxlev,
           (long long)byl[T.maxlev < 16 ? T.maxlev : 15], (long long)leafs,
           leafs > 0 ? 100.0 * (double)byl[T.maxlev < 16 ? T.maxlev : 15] / (double)leafs : 0.0);
  }

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
  t0 = now_s();
  if (hz_pocc_build(&m, &T, tlo, thi, tpl, occ, NULL) != 0) {
    fprintf(stderr, "отказ занятости\n");
    free(tlo);
    free(thi);
    free(tpl);
    free(list);
    free(occ);
    hz_ptree_free(&T);
    hz_obj_free(&m);
    return 2;
  }
  int64_t nocc = 0;
  for (int32_t i = 0; i < T.nnd; i++)
    if (occ[i]) nocc++;
  printf("== ЗАНЯТОСТЬ за %.2f с: узлов с геометрией %lld (%.2f %%)\n", now_s() - t0,
         (long long)nocc, 100.0 * (double)nocc / (double)T.nnd);

  /* СТРУКТУРНАЯ ПРОВЕРКА ДЕРЕВА (§306). §305 показал: обход накрывает треть
   * листьев и `17 499` раз останавливается в узле повторно. Из устройства обхода
   * это не следует, значит подозрение на САМО ДЕРЕВО: у каждого узла обязан быть
   * ровно один родитель, у корня — ни одного. У `octree.c` такая проверка есть
   * (`twoparent`, `unreachable`), у `ptree` её нет ни в каком виде. */
  {
    int32_t *par = calloc((size_t)T.nnd, sizeof *par);
    if (par != NULL) {
      int64_t bad_idx = 0;
      for (int32_t q = 0; q < T.nnd; q++) {
        if (T.nd[q].child < 0) continue;
        for (int k = 0; k < 8; k++) {
          int32_t ch = T.nd[q].child + k;
          if (ch < 0 || ch >= T.nnd) {
            bad_idx++;
            continue;
          }
          par[ch]++;
        }
      }
      int64_t p0 = 0, p1 = 0, p2 = 0, pmax = 0;
      for (int32_t q = 0; q < T.nnd; q++) {
        if (par[q] == 0)
          p0++;
        else if (par[q] == 1)
          p1++;
        else
          p2++;
        if (par[q] > pmax) pmax = par[q];
      }
      printf("== СТРУКТУРА ДЕРЕВА: узлов %d; БЕЗ РОДИТЕЛЯ %lld (обязан 1 — корень), с ОДНИМ "
             "%lld, с ДВУМЯ И БОЛЕЕ %lld (обязано 0, наибольшее %lld); ссылок вне границ %lld\n",
             T.nnd, (long long)p0, (long long)p1, (long long)p2, (long long)pmax,
             (long long)bad_idx);
    }
    free(par);
  }

  hz_pfront_ctx X;
  memset(&X, 0, sizeof X);
  X.T = &T;
  X.m = &m;
  X.tlo = tlo;
  X.thi = thi;
  X.occ = occ;
  X.noshadow = !shadow;
  X.cliplev = cliplev;
  X.frustfull = camfull;
  double L = sqrt(sdir[0] * sdir[0] + sdir[1] * sdir[1] + sdir[2] * sdir[2]);
  if (!(L > 0.0)) L = 1.0;
  for (int c = 0; c < 3; c++)
    X.dir[c] = sdir[c] / L;
  /* ВЫБОРКИ ПО ДИСКУ ИСТОЧНИКА (Ф4, §275). `K = 1` — в точности О73 и он же
   * негативный контроль. Выборка — спираль Ферма, детерминированная; у эталона
   * она ДРУГАЯ и вчетверо гуще (А533), иначе сверка сойдётся по построению. */
  X.nk = (nk < 1) ? 1 : (nk > HZ_PFRONT_MAXK ? HZ_PFRONT_MAXK : nk);
  {
    double half = asun, t1[3] = {0.0, 0.0, 1.0}, q1[3], q2[3];
    if (fabs(X.dir[2]) > 0.9) {
      t1[0] = 1.0;
      t1[2] = 0.0;
    }
    q1[0] = X.dir[1] * t1[2] - X.dir[2] * t1[1];
    q1[1] = X.dir[2] * t1[0] - X.dir[0] * t1[2];
    q1[2] = X.dir[0] * t1[1] - X.dir[1] * t1[0];
    double lq = sqrt(q1[0] * q1[0] + q1[1] * q1[1] + q1[2] * q1[2]);
    if (!(lq > 0.0)) lq = 1.0;
    for (int c = 0; c < 3; c++)
      q1[c] /= lq;
    q2[0] = X.dir[1] * q1[2] - X.dir[2] * q1[1];
    q2[1] = X.dir[2] * q1[0] - X.dir[0] * q1[2];
    q2[2] = X.dir[0] * q1[1] - X.dir[1] * q1[0];
    for (int s = 0; s < X.nk; s++) {
      double rr = (X.nk == 1) ? 0.0 : half * sqrt(((double)s + 0.5) / (double)X.nk);
      double ph = 2.39996322972865332 * (double)s;
      double nn = 0.0;
      for (int c = 0; c < 3; c++) {
        X.sdir[s][c] = X.dir[c] + rr * (cos(ph) * q1[c] + sin(ph) * q2[c]);
        nn += X.sdir[s][c] * X.sdir[s][c];
      }
      nn = sqrt(nn);
      if (!(nn > 0.0)) nn = 1.0;
      for (int c = 0; c < 3; c++)
        X.sdir[s][c] /= nn;
    }
  }
  /* Умолчание камерного пола: он включается вместе с ленивым спуском, потому
   * что без спуска его просто нечему ограничивать (§340). */
  if (g_camfloor < 0) g_camfloor = g_virt;
  X.pxeps2 = (px * HZ_CFG_EPS) * (px * HZ_CFG_EPS);
  /* §340: камерный пол. Глаз тот же, что у картинки и у пирамиды. */
  for (int c = 0; c < 3; c++)
    X.eye[c] = camo[c];
  X.camfloor = g_camfloor;
  X.t_entry = 1e300;
  for (int k = 0; k < 8; k++) {
    double t = 0.0;
    for (int c = 0; c < 3; c++)
      t += ((k & (1 << c)) ? T.nd[0].hi[c] : T.nd[0].lo[c]) * X.dir[c];
    if (t < X.t_entry) X.t_entry = t;
  }
  /* ПИРАМИДА КАМЕРЫ (ключ `cam=1`). Глаз и цель — из `scene_cfg.h`, где камеры
   * НАЙДЕНЫ ЗАМЕРОМ, а не назначены (§187); поле зрения и разрешение оттуда же.
   * Плоскости строятся нормалями ВНУТРЬ. */
  /* Камера нужна и картинке (`img`), а не только пометкам: без этого условия
   * прогон с `mark=0` снимал кадр из точки (0,0,0) — поймано сравнением
   * картинок, где различие вышло во весь сигнал. */
  if (cam || usemark || img) {
    double e[3], at[3];
    if (strstr(argv[1], "conference") != NULL) {
      double a1[3] = HZ_CFG_HALL_EYE, a2[3] = HZ_CFG_HALL_AT;
      memcpy(e, a1, sizeof e);
      memcpy(at, a2, sizeof at);
    } else if (strstr(argv[1], "rungholt") != NULL) {
      double a1[3] = HZ_CFG_CITY_EYE, a2[3] = HZ_CFG_CITY_AT;
      memcpy(e, a1, sizeof e);
      memcpy(at, a2, sizeof at);
    } else {
      double a1[3] = HZ_CFG_MIGUEL_EYE, a2[3] = HZ_CFG_MIGUEL_AT;
      memcpy(e, a1, sizeof e);
      memcpy(at, a2, sizeof at);
    }
    double f[3], up[3] = {0.0, 1.0, 0.0}, r[3], u2[3];
    double ln = 0.0;
    for (int c = 0; c < 3; c++) {
      f[c] = at[c] - e[c];
      ln += f[c] * f[c];
    }
    ln = sqrt(ln);
    if (!(ln > 0.0)) ln = 1.0;
    for (int c = 0; c < 3; c++)
      f[c] /= ln;
    r[0] = f[1] * up[2] - f[2] * up[1];
    r[1] = f[2] * up[0] - f[0] * up[2];
    r[2] = f[0] * up[1] - f[1] * up[0];
    ln = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
    if (!(ln > 0.0)) ln = 1.0;
    for (int c = 0; c < 3; c++)
      r[c] /= ln;
    u2[0] = r[1] * f[2] - r[2] * f[1];
    u2[1] = r[2] * f[0] - r[0] * f[2];
    u2[2] = r[0] * f[1] - r[1] * f[0];
    double th = 0.5 * HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0;
    double cs = cos(th), sn = sin(th);
    /* Четыре боковые: нормаль внутрь есть cs·(ось) ± sn·(вперёд) с точностью до
     * знака; ближняя — плоскость через глаз вперёд. */
    double nn[6][3];
    for (int c = 0; c < 3; c++) {
      nn[0][c] = cs * r[c] + sn * f[c];
      nn[1][c] = -cs * r[c] + sn * f[c];
      nn[2][c] = cs * u2[c] + sn * f[c];
      nn[3][c] = -cs * u2[c] + sn * f[c];
      nn[4][c] = f[c];
      nn[5][c] = -f[c];
    }
    for (int k = 0; k < 6; k++) {
      double d = 0.0;
      for (int c = 0; c < 3; c++) {
        X.fr[k][c] = nn[k][c];
        d += nn[k][c] * e[c];
      }
      X.fr[k][3] = -d;
    }
    /* Дальняя плоскость: сцена целиком, то есть отодвинута за диагональ. */
    double diag = 0.0;
    for (int c = 0; c < 3; c++) {
      double s = T.nd[0].hi[c] - T.nd[0].lo[c];
      diag += s * s;
    }
    X.fr[5][3] += sqrt(diag);
    memcpy(camo, e, sizeof camo);
    memcpy(camf, f, sizeof camf);
    memcpy(camat, at, sizeof camat);
    /* Пирамида нужна и проходу ПОМЕТОК, и прямой проверке в обходе. Считается
     * она в обоих случаях, а ВКЛЮЧАЕТСЯ в обходе только по `cam=1`: иначе
     * пирамида и пометка сделали бы одну работу дважды. */
    X.usefrustum = cam;
    printf("== ПИРАМИДА КАМЕРЫ: глаз %.2f,%.2f,%.2f, цель %.2f,%.2f,%.2f, поле %.0f°;\n"
           "   ячейка ЦЕЛИКОМ снаружи берётся ОДНОЙ (огрубление), а не выбрасывается\n",
           e[0], e[1], e[2], at[0], at[1], at[2], HZ_CFG_FOV_DEG);
  }
  printf("== ФРОНТ: направление %.3f,%.3f,%.3f, пол %g px ОТ ИСТОЧНИКА, потолок огрубления %d, "
         "заслонение %s\n",
         X.dir[0], X.dir[1], X.dir[2], px, HZ_PFRONT_COARSEN_MAX, shadow ? "ВКЛ" : "ВЫКЛ");
  printf("== ОГОВОРКА ОПЕРАТОРА: %s\n", hz_pfront_note());

  /* ПРОХОД ПОМЕТОК ОТ КАМЕРЫ (ключ `mark=1`). Дёшев: идёт только по крупным
   * узлам. Кладёт каждому предельный уровень, который фронт читает третьим
   * условием остановки. */
  unsigned char *dep = calloc((size_t)T.nnd, 1);
  if (dep == NULL) return 2;
  hz_pocc_depth(&T, dep);
  X.depth = dep;
  unsigned char *mk = NULL;
  if (usemark) {
    mk = calloc((size_t)T.nnd, 1);
    if (mk == NULL) return 2;
    for (int32_t i = 0; i < T.nnd; i++)
      mk[i] = 255; /* 255 — пометки нет */
  }
  if (usemark == 1) {
    /* ПОМЕТКИ БУФЕРОМ ДАЛЬНОСТИ ОТ ГЛАЗА (§278, `pcull.h`). Здесь камера —
     * ТОЧКА, а обе границы (заслонитель по дальней точке накрытого ЦЕЛИКОМ
     * пикселя, ячейка по ближней своей точке против максимума по всей проекции)
     * повёрнуты в сторону «пометок меньше». Приёмка обязана дать НОЛЬ. */
    double diag = 0.0;
    for (int c = 0; c < 3; c++) {
      double s = T.nd[0].hi[c] - T.nd[0].lo[c];
      diag += s * s;
    }
    diag = sqrt(diag);
    double upv[3] = HZ_CFG_UP;
    hz_pcull C;
    double tmk = now_s();
    if (hz_pcull_init(&C, camo, camat, upv, HZ_CFG_FOV_DEG, bufside, diag) != 0) {
      fprintf(stderr, "отказ буфера дальности\n");
      return 2;
    }
    C.loose = loose;
    /* Пирамида кадра отдаётся `pcull` ВСЕГДА, а не по ключу `cam`: она отвечает
     * на вопрос «виден ли этот кусок пространства», а не служит послаблением
     * пола. Вклады «вне кадра» и «заслонено» печатаются порознь. */
    memcpy(C.fr, X.fr, sizeof C.fr);
    C.usefr = 1;
    hz_pcull_draw(&C, &m);
    hz_pcull_pyramid(&C);
    int64_t nset = 0;
    hz_pcull_marks(&C, &T, mlev, mk);
    /* ДОЛЯ СЦЕНЫ, ЗАСЛОНЁННАЯ ОТ КАМЕРЫ, — величина про ПРЕДСТАВЛЕНИЕ, а не про
     * фронт (замечание пользователя 08-07). §279 мерил пометки против пола
     * СОЛНЕЧНОГО прохода и получил полтора процента; но это пересечение двух
     * разных множеств, а не размер заслонённого. Здесь считается прямо: сколько
     * ЗАНЯТЫХ ЛИСТЬЕВ — то есть самой мелкой части представления — лежит под
     * заслонёнными узлами. Лист взят мерой сознательно: он и есть то, что
     * предлагается огрублять.
     *
     * ПИРАМИДА КАДРА В ЭТО ЧИСЛО НЕ ВХОДИТ: `pcull` судит невидимость ТОЛЬКО по
     * заслонению одних тел другими, а всё, что за краем кадра, он не помечает
     * вовсе. Значит доля НЕ ЗАВЫШЕНА за счёт того, что просто не попало в кадр. */
    {
      int64_t leaf_all = 0, leaf_hid = 0, occ_all = 0, occ_hid = 0, leaf_occl = 0, leaf_frust = 0;
      /* ПОЛНОТА МЕХАНИЗМА — ВТОРАЯ ПОЛОВИНА ВОПРОСА (пользователь 08-07: «а
       * наоборот мерили? сколько невидимых, но непомеченных?»). Приёмка §278
       * односторонняя и проверяет только, что среди ПОМЕЧЕННЫХ нет видимых;
       * сколько НЕВИДИМОГО механизм пропускает, не мерилось ни разу.
       *
       * ЭТАЛОН ЗДЕСЬ ДРУГОЙ И НАРОЧНО ГРУБЫЙ: ячейка считается видимой, если
       * ХОТЬ ОДНА из девяти её точек (центр и восемь углов) лежит в кадре и
       * достижима лучом от глаза. Ошибается он в сторону «невидимых больше»
       * (тонкий просвет между девятью точками он пропустит), а значит
       * ЗАНИЖАЕТ измеряемую полноту — то есть работает против нужного ответа,
       * как эталону и положено. */
      int32_t recn = 0;
      int32_t recstride = (nrec > 0) ? 1 : 1;
      int32_t *recid = NULL;
      unsigned char *recmk = NULL;
      if (nrec > 0) {
        recid = malloc((size_t)nrec * sizeof *recid);
        recmk = malloc((size_t)nrec);
        int64_t nleafocc = 0;
        for (int32_t i2 = 0; i2 < T.nnd; i2++)
          if (occ[i2] && T.nd[i2].child < 0) nleafocc++;
        recstride = (int32_t)(nleafocc / (int64_t)nrec);
        if (recstride < 1) recstride = 1;
      }
      int32_t *st = malloc((size_t)T.nnd * sizeof *st);
      unsigned char *hid = malloc((size_t)T.nnd);
      if (st != NULL && hid != NULL) {
        int32_t sp = 0;
        st[sp] = 0;
        hid[0] = (unsigned char)(mk[0] != 255 ? mk[0] : 0);
        sp = 1;
        while (sp > 0) {
          int32_t nid = st[--sp];
          int h = hid[nid];
          if (occ[nid]) {
            occ_all++;
            if (h) occ_hid++;
            if (T.nd[nid].child < 0) {
              leaf_all++;
              if (h) leaf_hid++;
              if (h == 1) leaf_occl++;
              if (h == 2) leaf_frust++;
              /* Выборка для замера ПОЛНОТЫ (см. ниже): каждый `recstride`-й
               * занятый лист вместе с его пометкой. */
              if (nrec > 0 && recn < nrec && (leaf_all % recstride) == 0) {
                recid[recn] = nid;
                recmk[recn] = (unsigned char)h;
                recn++;
              }
            }
          }
          if (T.nd[nid].child >= 0)
            for (int k = 0; k < 8; k++) {
              int32_t c = T.nd[nid].child + k;
              /* Причина наследуется вниз: заслонённый кусок пространства не
               * становится видимым от того, что его поделили. `1` — заслонено
               * телами, `2` — вне кадра; при споре побеждает та, что выше. */
              hid[c] = (unsigned char)(h != 0 ? h : (mk[c] != 255 ? mk[c] : 0));
              st[sp++] = c;
            }
        }
        printf(
            "== НЕ ВИДНО ИЗ КАМЕРЫ — ДОЛЯ ПРЕДСТАВЛЕНИЯ (занятый лист есть самое мелкое,\n"
            "   что в дереве есть, и ровно то, что предлагается огрублять):\n"
            "   ВСЕГО %lld из %lld (%.1f %%); из них ЗАСЛОНЕНО ТЕЛАМИ %lld (%.1f %%), вне "
            "кадра %lld (%.1f %%)\n"
            "   занятых узлов %lld из %lld (%.1f %%)\n",
            (long long)leaf_hid, (long long)leaf_all,
            leaf_all > 0 ? 100.0 * (double)leaf_hid / (double)leaf_all : 0.0, (long long)leaf_occl,
            leaf_all > 0 ? 100.0 * (double)leaf_occl / (double)leaf_all : 0.0,
            (long long)leaf_frust,
            leaf_all > 0 ? 100.0 * (double)leaf_frust / (double)leaf_all : 0.0, (long long)occ_hid,
            (long long)occ_all, occ_all > 0 ? 100.0 * (double)occ_hid / (double)occ_all : 0.0);
      }
      /* ЗАМЕР ПОЛНОТЫ: по выборке занятых листьев сверяем пометку с эталоном. */
      if (recid != NULL && recmk != NULL && recn > 0) {
        double trec = now_s();
        int64_t ninv = 0, ninv_mk = 0, nvis = 0, nvis_mk = 0;
#pragma omp parallel for schedule(dynamic, 4) reduction(+ : ninv, ninv_mk, nvis, nvis_mk)
        for (int32_t q = 0; q < recn; q++) {
          const hz_ptnode *N2 = &T.nd[recid[q]];
          int seen = 0;
          for (int pt = 0; pt < 9 && !seen; pt++) {
            double P0[3];
            for (int c = 0; c < 3; c++)
              P0[c] = (pt == 8) ? 0.5 * (N2->lo[c] + N2->hi[c])
                                : ((pt & (1 << c)) ? N2->hi[c] : N2->lo[c]);
            int infr = 1;
            for (int kf = 0; kf < 6 && infr; kf++) {
              const double *PL = X.fr[kf];
              double sf = PL[3];
              for (int c = 0; c < 3; c++)
                sf += P0[c] * PL[c];
              if (sf < 0.0) infr = 0;
            }
            if (!infr) continue;
            double dv[3], len = 0.0;
            for (int c = 0; c < 3; c++) {
              dv[c] = P0[c] - camo[c];
              len += dv[c] * dv[c];
            }
            len = sqrt(len);
            if (!(len > 0.0)) continue;
            for (int c = 0; c < 3; c++)
              dv[c] /= len;
            int blk = 0;
            for (int32_t t = 0; t < m.nt && !blk; t++) {
              const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
              const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
              const double *Cc = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
              double q1[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
              double q2[3] = {Cc[0] - A[0], Cc[1] - A[1], Cc[2] - A[2]};
              double pv[3], tv[3], qv[3];
              pv[0] = dv[1] * q2[2] - dv[2] * q2[1];
              pv[1] = dv[2] * q2[0] - dv[0] * q2[2];
              pv[2] = dv[0] * q2[1] - dv[1] * q2[0];
              double det = q1[0] * pv[0] + q1[1] * pv[1] + q1[2] * pv[2];
              if (det > -1e-12 && det < 1e-12) continue;
              double inv = 1.0 / det;
              tv[0] = camo[0] - A[0];
              tv[1] = camo[1] - A[1];
              tv[2] = camo[2] - A[2];
              double uu = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
              if (uu < 0.0 || uu > 1.0) continue;
              qv[0] = tv[1] * q1[2] - tv[2] * q1[1];
              qv[1] = tv[2] * q1[0] - tv[0] * q1[2];
              qv[2] = tv[0] * q1[1] - tv[1] * q1[0];
              double vv = (dv[0] * qv[0] + dv[1] * qv[1] + dv[2] * qv[2]) * inv;
              if (vv < 0.0 || uu + vv > 1.0) continue;
              double tt = (q2[0] * qv[0] + q2[1] * qv[1] + q2[2] * qv[2]) * inv;
              if (tt > 1e-6 && tt < len - 1e-6) blk = 1;
            }
            if (!blk) seen = 1;
          }
          if (seen) {
            nvis++;
            if (recmk[q]) nvis_mk++;
          } else {
            ninv++;
            if (recmk[q]) ninv_mk++;
          }
        }
        printf("== ПОЛНОТА ПОМЕТОК (эталон лучами по 9 точкам ячейки, %d занятых листьев за "
               "%.1f с):\n"
               "   НЕВИДИМЫХ %lld, из них помечено %lld (%.1f %%) — пропущено %lld\n"
               "   ВИДИМЫХ   %lld, из них помечено %lld (обязано 0)\n",
               recn, now_s() - trec, (long long)ninv, (long long)ninv_mk,
               ninv > 0 ? 100.0 * (double)ninv_mk / (double)ninv : 0.0, (long long)(ninv - ninv_mk),
               (long long)nvis, (long long)nvis_mk);
      }
      free(recid);
      free(recmk);
      free(st);
      free(hid);
    }
    for (int32_t i = 0; i < T.nnd; i++)
      if (mk[i] != 255) nset++;
    printf("== БУФЕР ДАЛЬНОСТИ ОТ ГЛАЗА за %.2f с (%d²%s): пикселей занято %lld, треугольников "
           "мимо кадра %lld, у ближней плоскости %lld\n",
           now_s() - tmk, bufside,
           loose == 1   ? ", НЕГАТИВНЫЙ КОНТРОЛЬ: по СЕРЕДИНЕ пикселя"
           : loose == 2 ? ", НЕГАТИВНЫЙ КОНТРОЛЬ: БЛИЖНЯЯ точка пикселя"
                        : "",
           (long long)C.npix, (long long)C.ntri_off, (long long)C.ntri_near);
    printf("== ПОМЕТКИ: узлов проверено %lld, НЕ ВИДНО %lld (из них ВНЕ КАДРА %lld, ЗАСЛОНЕНО "
           "%lld), пометок %lld; отказов: за краем буфера %lld, у ближней плоскости %lld; "
           "огрубление на %d уровня\n",
           (long long)C.ntest, (long long)(C.nfrust + C.nhidden), (long long)C.nfrust,
           (long long)C.nhidden, (long long)nset, (long long)C.noff, (long long)C.nnear, relax);
    hz_pcull_free(&C);
    X.mark = mk;
    X.markcur = -1;
    X.markrelax = relax;
    /* ПОЛ ОЧЕРКА (§297): угловой радиус источника и верхняя оценка расстояния до
     * приёмника — поперечник сцены. Ключ `silh=0` выключает правило. */
    if (silh) {
      X.silh_a = asun;
      double dg = 0.0;
      for (int c = 0; c < 3; c++) {
        double sdd = T.nd[0].hi[c] - T.nd[0].lo[c];
        dg += sdd * sdd;
      }
      X.silh_d = sqrt(dg);
      X.silh_local = (silh >= 2);
      X.t_exit = -1e300;
      for (int k = 0; k < 8; k++) {
        double tt = 0.0;
        for (int c = 0; c < 3; c++)
          tt += ((k & (1 << c)) ? T.nd[0].hi[c] : T.nd[0].lo[c]) * X.dir[c];
        if (tt > X.t_exit) X.t_exit = tt;
      }
    }
  } else if (usemark) {
    /* ПРЕЖНИЙ ПРОХОД ФРОНТОМ (`mark=2`) — оставлен ради воспроизводимости чисел
     * §271/§272 и как то, ЧТО ИМЕННО не проходит приёмку. Это ТОТ ЖЕ фронт,
     * пущенный от камеры: где он приходит тёмным, там ставится пометка.
     * Направление берётся взглядом камеры — приближение параллельным пучком. */
    hz_pfront_ctx CX = X;
    CX.nk = 1; /* пометки ставит ОДНА выборка: невидимость двоична */
    for (int c = 0; c < 3; c++)
      CX.sdir[0][c] = camf[c];
    CX.mark = NULL;
    CX.markout = mk;
    CX.markoutlev = mlev;
    CX.markrelax = relax;
    CX.depth = dep;
    CX.cliplev = mlev;
    CX.usefrustum = 1;
    CX.nmarkset = 0;
    for (int c = 0; c < 3; c++)
      CX.dir[c] = camf[c];
    CX.t_entry = 1e300;
    for (int k = 0; k < 8; k++) {
      double tt = 0.0;
      for (int c = 0; c < 3; c++)
        tt += ((k & (1 << c)) ? T.nd[0].hi[c] : T.nd[0].lo[c]) * CX.dir[c];
      if (tt < CX.t_entry) CX.t_entry = tt;
    }
    /* Явное обнуление: заполняется цикл ниже лишь до , и доказательства
     * этого у cppcheck нет — он прав. */
    hz_pfront_face ci[3 * HZ_PFRONT_MAXK] = {{{0}}}, co[3 * HZ_PFRONT_MAXK] = {{{0}}};
    for (int a = 0; a < 3 * CX.nk; a++) {
      hz_pfront_set(&ci[a], 1.0);
    }
    double tmk = now_s();
    CX.virt = g_virt;
    hz_pfront_walk(&CX, 0, T.nd[0].lo, T.nd[0].hi, ci, co, list, m.nt, 0, 0, 0);
    printf("== ПРОХОД ОТ КАМЕРЫ за %.2f с: ячеек %lld, ПОМЕТОК ПОСТАВЛЕНО %lld (невидимое, "
           "огрубление на %d уровня)\n",
           now_s() - tmk, (long long)CX.ncell, (long long)CX.nmarkset, relax);
    X.mark = mk;
    X.markcur = -1;
    X.markrelax = relax;
  }
  /* СРЕЗ КАДРА — ЧТО ИМЕННО ОГРУБЛЕНИЕ НЕВИДИМОГО ЭКОНОМИТ (ключ `cut=1`).
   *
   * ЗАЧЕМ ОТДЕЛЬНО ОТ ФРОНТА. §279 мерил пометки против фронта и получил
   * `1.016×`, а §283 показал, что невидимо `56…88 %` представления. Оба числа
   * верны, и противоречия нет: **фронт останавливается КРУПНЕЕ той зернистости,
   * на которой невидимость вообще доказуема.** Пометки живут на уровнях `8…12`,
   * а солнечный фронт при своём поле встаёт на `5…7` — и внутрь помеченного
   * узла просто не спускается («остановок ПО ПОМЕТКЕ `358` из `195 364`»).
   * Значит мерить надо потребителя, работающего НА ЗЕРНИСТОСТИ ЛИСТА, а таков
   * срез кадра: набор ячеек, который кадру нужен от представления.
   *
   * ЧТО СЧИТАЕТСЯ. Спуск с полом камеры `4r² ≤ (px·ε)²·d²`, где `d` — расстояние
   * ОТ ГЛАЗА (а не вдоль солнца, как у фронта). Пустой узел берётся целиком
   * (§241.4). Три варианта: без пометок, с огрублением невидимого на `relax`
   * ступеней и с полным его схлопыванием — последнее есть ПОТОЛОК, и держать
   * его рядом обязательно (§279: потолок мерить до доводки).
   *
   * ОТСЕЧЕНИЯ ЗДЕСЬ НЕТ НИ В ОДНОМ ВАРИАНТЕ: невидимая ячейка БЕРЁТСЯ ЦЕЛИКОМ,
   * а не выбрасывается — свет приходит из-за кадра и от невидимых поверхностей
   * (указание пользователя 08-05, §259 граница 2). */
  if (usemark && mk != NULL) {
    static const char *NM[3] = {"без пометок", "невидимое грубее на relax", "невидимое ЦЕЛИКОМ"};
    double eps2 = (px * HZ_CFG_EPS) * (px * HZ_CFG_EPS);
    int64_t cnt[3] = {0, 0, 0};
    double tcut = now_s();
    for (int mode = 0; mode < 3; mode++) {
      int32_t *st = malloc((size_t)T.nnd * sizeof *st);
      unsigned char *hd = malloc((size_t)T.nnd);
      if (st == NULL || hd == NULL) {
        free(st);
        free(hd);
        break;
      }
      int32_t sp = 0;
      st[sp] = 0;
      hd[0] = (unsigned char)(mk[0] != 255);
      sp = 1;
      while (sp > 0) {
        int32_t nid = st[--sp];
        int h = hd[nid];
        const hz_ptnode *N = &T.nd[nid];
        int stop = 0;
        if (!occ[nid] || N->child < 0) {
          stop = 1;
        } else if (h && mode == 2) {
          stop = 1;
        } else {
          double r2 = 0.0, d2 = 0.0;
          for (int c = 0; c < 3; c++) {
            double hh = 0.5 * (N->hi[c] - N->lo[c]);
            double dc = 0.5 * (N->lo[c] + N->hi[c]) - camo[c];
            r2 += hh * hh;
            d2 += dc * dc;
          }
          double e2 = eps2;
          if (h && mode == 1)
            for (int k = 0; k < relax && k < 16; k++)
              e2 *= 4.0;
          if (d2 > r2 && 4.0 * r2 <= e2 * d2) stop = 1;
        }
        if (stop) {
          cnt[mode]++;
          continue;
        }
        for (int k = 0; k < 8; k++) {
          int32_t c = N->child + k;
          hd[c] = (unsigned char)(h || mk[c] != 255);
          st[sp++] = c;
        }
      }
      free(st);
      free(hd);
    }
    printf("== СРЕЗ КАДРА (пол камеры %g px ОТ ГЛАЗА, ячеек в представлении, "
           "нужном кадру), за %.2f с:\n",
           px, now_s() - tcut);
    for (int mode = 0; mode < 3; mode++)
      printf("   %-28s %10lld\n", NM[mode], (long long)cnt[mode]);
    if (cnt[1] > 0 && cnt[2] > 0)
      printf("   выигрыш: огрубление на %d ступени — %.2fx, ПОТОЛОК — %.2fx\n", relax,
             (double)cnt[0] / (double)cnt[1], (double)cnt[0] / (double)cnt[2]);
  }

  /* Обнуление явное — по той же причине, что у `ci`/`co` выше. */
  hz_pfront_face in[3 * HZ_PFRONT_MAXK] = {{{0}}}, out[3 * HZ_PFRONT_MAXK] = {{{0}}};
  for (int a = 0; a < 3 * X.nk; a++) {
    hz_pfront_set(&in[a], 1.0); /* на входе в сцену диск источника открыт целиком */
  }
  double *refpt = NULL, *refvis = NULL, *refh = NULL;
  if (nref > 0) {
    refpt = malloc(3 * (size_t)nref * sizeof *refpt);
    refvis = malloc((size_t)nref * sizeof *refvis);
    refh = malloc((size_t)nref * sizeof *refh);
    X.refh = refh;
    X.refpt = refpt;
    X.refvis = refvis;
    X.refcap = nref;
    X.refstride = 997; /* простое число: выборка не попадает в такт обхода */
  }
  t0 = now_s();
  X.refinemax = refinemax;
  int32_t *vis = NULL;
  if (img) {
    vis = calloc((size_t)T.nnd, sizeof *vis);
    X.visit = vis;
  }
  /* Место, куда фронт положит свой ответ, — только если он кому-то нужен. */
  if (img) {
    X.cellf = malloc((size_t)T.nnd * sizeof *X.cellf);
    if (X.cellf != NULL)
      for (int32_t i = 0; i < T.nnd; i++)
        X.cellf[i] = -2.0f; /* сторож «не записано»: отличаем от ОТРИЦАТЕЛЬНОГО ответа */
    /* Таблица глубин заслонителя внутри ячейки (§296). Запись — 64 числа на
     * ОСТАНОВИВШУЮСЯ ячейку; их доли процента от узлов, поэтому индекс на узел
     * плюс плотный массив записей, а не 64 числа на каждый узел. */
    if (indep) {
      X.cellslot = malloc((size_t)T.nnd * sizeof *X.cellslot);
      X.capslot = 4000000;
      X.celldep = malloc(HZ_PFRONT_NCELL * (size_t)X.capslot * sizeof *X.celldep);
      X.cellval = malloc(HZ_PFRONT_NCELL * (size_t)X.capslot * sizeof *X.cellval);
      if (X.cellslot != NULL && X.celldep != NULL && X.cellval != NULL)
        for (int32_t i = 0; i < T.nnd; i++)
          X.cellslot[i] = -1;
      else {
        free(vis);
        free(X.cellslot);
        free(X.celldep);
        free(X.cellval);
        X.cellslot = NULL;
        X.celldep = NULL;
      }
    }
  }
  X.virt = g_virt;
  hz_pfront_walk(&X, 0, T.nd[0].lo, T.nd[0].hi, in, out, list, m.nt, 0, 0, 0);
  double secs = now_s() - t0;
  if (X.fail) {
    fprintf(stderr, "отказ обхода\n");
    free(tlo);
    free(thi);
    free(tpl);
    free(list);
    free(occ);
    hz_ptree_free(&T);
    hz_obj_free(&m);
    return 2;
  }
  printf("   ЯЧЕЕК %lld (с геометрией %lld, пустых %lld)\n", (long long)X.ncell,
         (long long)X.ncell_geo, (long long)X.ncell_void);
  printf("   из них ВИРТУАЛЬНЫХ (без узла, §338): %lld (%.2f %%)\n", (long long)X.nvirt,
         X.ncell > 0 ? 100.0 * (double)X.nvirt / (double)X.ncell : 0.0);
  printf("   спуск остановлен: ПУСТОТОЙ %lld (%.2f %%), ЛИСТОМ %lld (%.2f %%), ПОЛОМ %lld "
         "(%.2f %%), ТЕМНЫМ %lld (%.2f %%)\n",
         (long long)X.stop_void, 100.0 * (double)X.stop_void / (double)X.ncell,
         (long long)X.stop_leaf, 100.0 * (double)X.stop_leaf / (double)X.ncell,
         (long long)X.stop_floor, 100.0 * (double)X.stop_floor / (double)X.ncell,
         (long long)X.stop_flat, 100.0 * (double)X.stop_flat / (double)X.ncell);
  printf("   ячеек ОСВЕЩЁННЫХ %lld, В ТЕНИ %lld; потолок огрубления сработал %lld раз; предел "
         "кусков %lld раз\n",
         (long long)X.nlit, (long long)X.nshadow, (long long)X.ncap, (long long)X.npclip);
  printf("   вне пирамиды с ОТПУЩЕННЫМ полом %lld ячеек\n", (long long)X.noutside);
  printf("   остановок ПО ПОМЕТКЕ %lld; ПОЛ ОЧЕРКА запретил огрубление %lld раз\n",
         (long long)X.nmarkstop, (long long)X.nsilh);
  printf("   ЯЧЕЕК С ГЕОМЕТРИЕЙ ОДНОЙ ПЛОСКОСТИ %lld из %lld (%.2f %%) — их дробить незачем ни при "
         "каком поле\n",
         (long long)X.nflat1, (long long)X.nflatgeo,
         X.nflatgeo > 0 ? 100.0 * (double)X.nflat1 / (double)X.nflatgeo : 0.0);
  printf("   ВРЕМЯ %.2f с, на ячейку %.1f нс\n", secs,
         X.ncell > 0 ? 1e9 * secs / (double)X.ncell : 0.0);

  {
    printf("== КАНДИДАТОВ У ФРОНТА ПО УРОВНЯМ (узлов / кандидатов / в среднем):\n");
    for (int lv = 0; lv < 16; lv++)
      if (X.nodelev[lv] > 0)
        printf("   ур.%2d  узлов %8lld  кандидатов %12lld  в среднем %8.1f\n", lv,
               (long long)X.nodelev[lv], (long long)X.candlev[lv],
               (double)X.candlev[lv] / (double)X.nodelev[lv]);
  }
  if (vis != NULL) {
    int64_t v0 = 0, v1 = 0, v2 = 0, vmax = 0, v0leaf = 0, v0occ = 0;
    for (int32_t q = 0; q < T.nnd; q++) {
      if (vis[q] == 0) {
        v0++;
        if (T.nd[q].child < 0) v0leaf++;
        if (occ[q]) v0occ++;
      } else if (vis[q] == 1)
        v1++;
      else
        v2++;
      if (vis[q] > vmax) vmax = vis[q];
    }
    /* ДО КАКОЙ ГЛУБИНЫ ДЕРЕВО ВООБЩЕ ЧИТАЕТСЯ (§318). Главное требование к
     * иерархии — быстрый проход фронта и дешёвая правка при разрушении; оба
     * зависят не от того, как глубоко дерево ПОСТРОЕНО, а от того, как глубоко
     * его ЧИТАЮТ. Разница между этими двумя числами и есть чистые потери. */
    {
      int64_t vl[16] = {0};
      for (int32_t q = 0; q < T.nnd; q++)
        if (vis[q] > 0) {
          int lv = (int)T.lev[q];
          if (lv > 15) lv = 15;
          vl[lv]++;
        }
      printf("   ПОСЕЩЕНО ПО УРОВНЯМ:");
      for (int lv = 0; lv < 16; lv++)
        if (vl[lv] > 0) printf(" %d:%lld", lv, (long long)vl[lv]);
      printf("\n");
    }
    printf("== ПОСЕЩЕНИЯ ОБХОДОМ: не посещено %lld узлов (из них листьев %lld, занятых %lld), "
           "один раз %lld, БОЛЕЕ ОДНОГО %lld (наибольшее %lld)\n",
           (long long)v0, (long long)v0leaf, (long long)v0occ, (long long)v1, (long long)v2,
           (long long)vmax);
  }

  /* ПОКРЫВАЕТ ЛИ ФРОНТ ВСЁ ПРОСТРАНСТВО (§303, шаг 1). Устройство обхода
   * утверждает: остановки фронта разбивают дерево без дыр — на листе остановка
   * безусловна. Утверждение НИ РАЗУ не проверялось, а §303 нашёл `86 845`
   * пикселей, читающих пустоту. Счёт прямой: у каждого листа обязан быть предок
   * (или он сам) с записью. */
  if (img && X.cellf != NULL) {
    int64_t leaf_tot = 0, leaf_cov = 0;
    int32_t *st2 = malloc((size_t)T.nnd * sizeof *st2);
    unsigned char *cv = malloc((size_t)T.nnd);
    if (st2 != NULL && cv != NULL) {
      int32_t sp2 = 0;
      st2[sp2] = 0;
      cv[0] = (unsigned char)(X.cellf[0] >= 0.0f);
      sp2 = 1;
      while (sp2 > 0) {
        int32_t q = st2[--sp2];
        int c2 = cv[q];
        if (T.nd[q].child < 0) {
          leaf_tot++;
          if (c2) leaf_cov++;
          continue;
        }
        for (int k = 0; k < 8; k++) {
          int32_t ch = T.nd[q].child + k;
          cv[ch] = (unsigned char)(c2 || X.cellf[ch] >= 0.0f);
          st2[sp2++] = ch;
        }
      }
      /* ВТОРОЙ СЧЁТ, НАПИСАННЫЙ НЕЗАВИСИМО (§305). Первый распространяет признак
       * «накрыт» СВЕРХУ ВНИЗ; этот считает листья ПОД каждой остановкой СНИЗУ
       * ВВЕРХ и складывает. Вложенных остановок быть не может — ниже остановки
       * обход не идёт, — поэтому у разбиения сумма обязана дать ровно все листья.
       * Два счёта, устроенные по-разному, отвечают на один вопрос: если сойдутся,
       * дефект в обходе; если разойдутся, дефект в счёте. */
      int64_t *lv2 = calloc((size_t)T.nnd, sizeof *lv2);
      if (lv2 != NULL) {
        for (int32_t q = T.nnd - 1; q >= 0; q--) {
          if (T.nd[q].child < 0)
            lv2[q] = 1;
          else
            for (int k = 0; k < 8; k++)
              lv2[q] += lv2[T.nd[q].child + k];
        }
        int64_t nstop = 0, sum2 = 0, nneg = 0, sumneg = 0;
        for (int32_t q = 0; q < T.nnd; q++) {
          if (X.cellf[q] > -1.5f) { /* запись ЕСТЬ (в т.ч. отрицательная) */
            nstop++;
            sum2 += lv2[q];
            if (X.cellf[q] < 0.0f) {
              nneg++;
              sumneg += lv2[q];
            }
          }
        }
        printf("== ОТРИЦАТЕЛЬНОЕ СОСТОЯНИЕ: записей с f < 0 — %lld из %lld, листьев под ними "
               "%lld\n",
               (long long)nneg, (long long)nstop, (long long)sumneg);
        printf("== ВТОРОЙ СЧЁТ: остановок %lld (обход насчитал %lld), листьев под ними %lld из "
               "%lld (%.2f %%)\n",
               (long long)nstop, (long long)X.ncell, (long long)sum2, (long long)leaf_tot,
               leaf_tot > 0 ? 100.0 * (double)sum2 / (double)leaf_tot : 0.0);
      }
      free(lv2);
      printf("== ПОКРЫТИЕ ФРОНТОМ: листьев %lld, накрыто остановками %lld (%.2f %%), ДЫРА %lld\n",
             (long long)leaf_tot, (long long)leaf_cov,
             leaf_tot > 0 ? 100.0 * (double)leaf_cov / (double)leaf_tot : 0.0,
             (long long)(leaf_tot - leaf_cov));
    }
    free(st2);
    free(cv);
  }

  /* ПРОХОД 3 — СБОР ПО ПИКСЕЛЮ, ТО ЕСТЬ КАРТИНКА (ключ `img=1`).
   *
   * ЗАЧЕМ ОН ЗДЕСЬ И ПОЧЕМУ ТОЛЬКО ТЕПЕРЬ. У линии фронта не было ни одной
   * картинки: фронт считал долю открытого диска и ВЫБРАСЫВАЛ её, поэтому всякий
   * выигрыш оставался числом в счётчике. Замечание пользователя §277 («эти
   * оптимизации ничего не позволили нового») упирается ровно в это. Теперь фронт
   * оставляет ответ в `cellf`, и его можно посмотреть глазами.
   *
   * УСТРОЙСТВО, ПРОСТЕЙШЕЕ ИЗ ЗАКОННЫХ (§7: «камера на первое время как угодно»):
   * обычный z-буфер с номером треугольника даёт видимую поверхность в пикселе,
   * точка попадания — из дальности, ячейка — спуском по дереву до той, где фронт
   * остановился, яркость — `ρ·f·|n·ω|/π`.
   *
   * ЛУЧЕЙ ЗДЕСЬ НЕТ И БЫТЬ НЕ МОЖЕТ: дерево строится укладкой резкой и ссылок на
   * треугольники не хранит вовсе (А472), так что `ptrace` по нему не пойдёт.
   * Проекция отвечает на тот же вопрос разом для всех пикселей. */
  if (img) {
    double timg = now_s();
    size_t np = (size_t)bufside * (size_t)bufside;
    float *zb = malloc(np * sizeof *zb);
    int32_t *ib = malloc(np * sizeof *ib);
    double *pix = malloc(np * sizeof *pix);
    /* Доля открытого диска ОТДЕЛЬНО от яркости: эталон спрашивает «видно ли
     * солнце», а не «насколько ярко», и мешать в одну величину косинус с
     * альбедо значило бы сверять три вещи разом. */
    double *fpix = malloc(np * sizeof *fpix);
    double *fromcell = NULL; /* §351: доля, собранная ИЗ ЯЧЕЕК */
    if (zb != NULL && ib != NULL && pix != NULL && fpix != NULL) {
      double diag = 0.0;
      for (int c = 0; c < 3; c++) {
        double s = T.nd[0].hi[c] - T.nd[0].lo[c];
        diag += s * s;
      }
      double upv[3] = HZ_CFG_UP;
      hz_pcull C2;
      if (hz_pcull_init(&C2, camo, camat, upv, HZ_CFG_FOV_DEG, bufside, sqrt(diag)) == 0) {
        /* РИСОВАНИЕ СПЕРЕДИ НАЗАД ПО ДЕРЕВУ (§316) вместо z-буфера: пиксель,
         * однажды закрытый, больше не трогается, и когда закрыты все — обход
         * прекращается. Ключ `zbuf=1` возвращает прежний перебор ВСЕХ
         * треугольников — он же негативный контроль: картинки обязаны совпасть
         * ПОБИТОВО, иначе порядок обхода неверен. */
        int64_t sh_tri = 0, sh_cell = 0;
        /* КАРТИНКА ИЗ ЯЧЕЕК ФРОНТА (§349, вопрос пользователя «зачем нам
         * z-буфер»). Видимость считается ОДИН раз и в той же структуре, в
         * которой живёт поле: проецируется коробка ОСТАНОВИВШЕЙСЯ ячейки, в
         * пиксель пишется её состояние, глубина сравнивается ПО ЯЧЕЙКАМ.
         * Треугольники здесь не участвуют вовсе — а значит исчезает и стык двух
         * видимостей, породивший все 117 расхождений (§347).
         *
         * ЦЕНА ЭТОГО ХОДА НАЗВАНА ЗАРАНЕЕ: картинка выходит с разрешением ПОЛЯ,
         * а не пикселя, потому что мельче ячейки в ней ничего нет. */
        /* §351: ЭТАЛОНУ ГЕОМЕТРИЯ ЗАКОННА — он и есть истина, и отнимать у него
         * точную видимость было ошибкой §349 (эталон и проверяемое выродились
         * одновременно). Треугольный z-буфер остаётся ЕМУ; путь из ячеек
         * заполняет только `fpix`, со своей глубиной. Тогда сравниваются две
         * ОПРЕДЕЛЁННЫЕ величины: поле, собранное из ячеек, против точной доли
         * диска в видимой точке, и разность есть цена квантования ячейкой. */
        if (cellimg) {
          double tci = now_s();
          fromcell = malloc(np * sizeof *fromcell);
          if (fromcell == NULL) return 2;
          for (size_t p = 0; p < np; p++)
            fromcell[p] = 0.0;
          int64_t ndrawn = 0;
          paint_cell(&T, 0, X.cellf, &C2, camo, bufside, fromcell, &ndrawn);
          printf("== КАРТИНКА ИЗ ЯЧЕЕК (§349) за %.2f с: поставлено ячеек %lld из %d узлов\n",
                 now_s() - tci, (long long)ndrawn, T.nnd);
        }
        {
          double tsh = now_s();
          memcpy(C2.fr, X.fr, sizeof C2.fr);
          C2.usefr = 1;
          if (zbuf)
            hz_pcull_shot(&C2, &m, zb, ib);
          else
            hz_pcull_shot_tree(&C2, &m, &T, tlo, thi, zb, ib, &sh_tri, &sh_cell);
          printf("== РИСОВАНИЕ %s за %.2f с: треугольников поставлено %lld, ячеек обойдено %lld "
                 "(всего в сцене %d)\n",
                 zbuf ? "z-БУФЕРОМ (перебор всех)" : "СПЕРЕДИ НАЗАД по дереву", now_s() - tsh,
                 (long long)sh_tri, (long long)sh_cell, m.nt);
        }
        int64_t nsky = 0, nlit2 = 0, nsh2 = 0, nback = 0, nmiss = 0;
        int pathshown = 0;
        double sum = 0.0;
        for (size_t p = 0; p < np; p++) {
          pix[p] = 0.0;
          /* При сборке из ячеек `fpix` уже заполнено и перечитывать его по
           * точке нельзя — точки в этом пути нет вовсе (§349). */
          if (!cellimg) fpix[p] = 0.0;
          if (ib[p] < 0) {
            nsky++;
            continue;
          }
          double d[3];
          hz_pcull_ray(&C2, (int)(p % (size_t)bufside), (int)(p / (size_t)bufside), d);
          double P[3];
          for (int c = 0; c < 3; c++)
            P[c] = camo[c] + (double)zb[p] * d[c];
          /* ЯЧЕЙКА, В КОТОРОЙ ФРОНТ ОСТАНОВИЛСЯ. Точка берётся СДВИНУТОЙ К
           * ИСТОЧНИКУ, и это не мелочь: поверхность часто лежит ровно на границе
           * ячеек, гасит свет для нижней — а точка попадания геометрически
           * принадлежит именно нижней, и читается ноль там, где поверхность
           * освещена. Досье обратной популяции показало этот случай прямо: ячейки
           * уровня 9…11 (ребро 4…16 мм), кусков в них НЕТ (`глубина в ячейке =
           * нет`), нормаль смотрит на солнце (`|n·ω| = 0.93`), а `f = 0`.
           * Сдвиг — тот же приём и та же величина, что у теневого луча эталона. */
          double PL[3];
          {
            double eo = ooff * (fabs(P[0]) + fabs(P[1]) + fabs(P[2]) + 1.0);
            for (int c = 0; c < 3; c++)
              PL[c] = P[c] - eo * X.dir[c];
          }
          /* СПУСК ДО ЛИСТА С ПАМЯТЬЮ О ПОСЛЕДНЕЙ ЗАПИСИ. Раньше спуск
           * останавливался на первом узле с записью, а если её не было нигде —
           * картинка УГАДЫВАЛА темноту. Замер показал цену угадывания: `86 845`
           * пикселей (15 % поверхностных) записи не имеют, и угадывание «темно»
           * даёт `117` ложных теней, «светло» — `5 290` ложных светов (§303).
           * Правильно — брать ответ ближайшего предка, у которого запись ЕСТЬ:
           * он покрывает эту точку по построению, потому что фронт остановился
           * именно на нём. */
          int32_t nid = 0, best = (X.cellf[0] > -1.5f) ? 0 : -1;
          while (T.nd[nid].child >= 0) {
            double mid[3];
            int k = 0;
            for (int c = 0; c < 3; c++) {
              mid[c] = 0.5 * (T.nd[nid].lo[c] + T.nd[nid].hi[c]);
              if (PL[c] >= mid[c]) k |= (1 << c);
            }
            nid = T.nd[nid].child + k;
            if (X.cellf[nid] > -1.5f) best = nid;
          }
          if (best >= 0) nid = best;
          /* Одноразовая печать ПУТИ для первого пикселя без записи: где именно
           * фронт не оставил ответа (§303). */
          if (best < 0 && ndbg != 0 && !pathshown) {
            pathshown = 1;
            printf("   ПУТЬ ДЛЯ ПИКСЕЛЯ БЕЗ ЗАПИСИ (%d):\n", (int)p);
            int32_t q = 0;
            for (int lv = 0; lv < 16; lv++) {
              printf("      ур.%2d узел %9d  ребро %.5f  занят=%d  запись=%s  дети=%s\n", lv, q,
                     T.nd[q].hi[0] - T.nd[q].lo[0], (int)occ[q],
                     X.cellf[q] >= 0.0f ? "ЕСТЬ" : "нет", T.nd[q].child >= 0 ? "есть" : "ЛИСТ");
              if (T.nd[q].child < 0) break;
              double mid[3];
              int k = 0;
              for (int c = 0; c < 3; c++) {
                mid[c] = 0.5 * (T.nd[q].lo[c] + T.nd[q].hi[c]);
                if (PL[c] >= mid[c]) k |= (1 << c);
              }
              q = T.nd[q].child + k;
            }
          }
          double f = (X.cellf[nid] >= 0.0f) ? (double)X.cellf[nid] : nofrec;
          /* КЛЕТКА ПОД ТОЧКОЙ, А НЕ СРЕДНЕЕ ПО ГРАНИ (О76). Состояние — сетка,
           * и поверхность обязана читать своё место на ней: иначе освещённая
           * клетка получает усреднённую тень соседей, и вся выгода сетки до
           * картинки не доходит. */
          if (X.cellval != NULL && X.cellslot[nid] >= 0) {
            int a3 = X.celldepax;
            int p3 = (a3 + 1) % 3, q3v = (a3 + 2) % 3;
            const hz_ptnode *N3 = &T.nd[nid];
            double h3p = N3->hi[p3] - N3->lo[p3], h3q = N3->hi[q3v] - N3->lo[q3v];
            double fa3 = (X.dir[a3] > 0.0) ? N3->lo[a3] : N3->hi[a3];
            double t3 = (P[a3] - fa3) / X.dir[a3];
            double up3 = P[p3] - t3 * X.dir[p3], vp3 = P[q3v] - t3 * X.dir[q3v];
            int gi3 = (int)((double)HZ_PFRONT_COVN * (up3 - N3->lo[p3]) / h3p);
            int gj3 = (int)((double)HZ_PFRONT_COVN * (vp3 - N3->lo[q3v]) / h3q);
            if (gi3 < 0) gi3 = 0;
            if (gj3 < 0) gj3 = 0;
            if (gi3 > HZ_PFRONT_COVN - 1) gi3 = HZ_PFRONT_COVN - 1;
            if (gj3 > HZ_PFRONT_COVN - 1) gj3 = HZ_PFRONT_COVN - 1;
            f = (double)X.cellval[HZ_PFRONT_NCELL * (size_t)X.cellslot[nid] +
                                  (size_t)(gj3 * HZ_PFRONT_COVN + gi3)];
          }
          if (X.cellf[nid] < 0.0f) nmiss++;
          /* ЗАСЛОНЕНИЕ ВНУТРИ ЯЧЕЙКИ (§296): спрашиваем свой столбик сетки, есть
           * ли геометрия БЛИЖЕ к источнику, чем мы сами. Дробить ячейку для
           * этого не нужно — нужна лишь глубина, посчитанная фронтом на месте. */
          if (X.cellslot != NULL && X.cellslot[nid] >= 0) {
            int a2 = X.celldepax;
            int p2 = (a2 + 1) % 3, q2 = (a2 + 2) % 3;
            const hz_ptnode *NN = &T.nd[nid];
            double fa = (X.dir[a2] > 0.0) ? NN->lo[a2] : NN->hi[a2];
            double du2 = NN->hi[p2] - NN->lo[p2], dv2 = NN->hi[q2] - NN->lo[q2];
            /* Своя глубина: путь от входной грани до точки вдоль света. */
            double tp = (P[a2] - fa) / X.dir[a2];
            double up = P[p2] - tp * X.dir[p2], vp = P[q2] - tp * X.dir[q2];
            int gi = (int)((double)HZ_PFRONT_COVN * (up - NN->lo[p2]) / du2);
            int gj = (int)((double)HZ_PFRONT_COVN * (vp - NN->lo[q2]) / dv2);
            if (gi < 0) gi = 0;
            if (gj < 0) gj = 0;
            if (gi > HZ_PFRONT_COVN - 1) gi = HZ_PFRONT_COVN - 1;
            if (gj > HZ_PFRONT_COVN - 1) gj = HZ_PFRONT_COVN - 1;
            double dq = (double)X.celldep[HZ_PFRONT_NCELL * (size_t)X.cellslot[nid] +
                                          (size_t)(gj * HZ_PFRONT_COVN + gi)];
            /* Допуск — ширина столбика: собственная поверхность точки внутри
             * столбика имеет разброс глубины того же порядка, и без допуска она
             * заслоняла бы сама себя. */
            double tolz = (du2 > dv2 ? du2 : dv2) / (double)intol;
            if (dq < tp - tolz) f = 0.0;
          }
          /* Косинус на поверхности — из нормали видимого треугольника. */
          int32_t t = ib[p];
          const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
          const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
          const double *Cc = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
          double u1[3], u2[3], nn[3];
          for (int c = 0; c < 3; c++) {
            u1[c] = B[c] - A[c];
            u2[c] = Cc[c] - A[c];
          }
          nn[0] = u1[1] * u2[2] - u1[2] * u2[1];
          nn[1] = u1[2] * u2[0] - u1[0] * u2[2];
          nn[2] = u1[0] * u2[1] - u1[1] * u2[0];
          double ln2 = sqrt(nn[0] * nn[0] + nn[1] * nn[1] + nn[2] * nn[2]);
          if (!(ln2 > 0.0)) ln2 = 1.0;
          /* СТОРОНА ПОВЕРХНОСТИ РЕШАЕТ, А НЕ МОДУЛЬ КОСИНУСА. Сперва здесь стояло
           * `|n·ω|` с оговоркой «нормаль .obj может смотреть внутрь», и это был
           * настоящий дефект: тыльная к солнцу сторона стены освещалась как
           * лицевая. Поймано КАРТОЙ ошибки — расхождение с эталоном шло не по
           * краям теней, а сплошными поверхностями (§287).
           *
           * Правильно так: `f` есть свойство ЯЧЕЙКИ («сколько диска ей видно»), а
           * не стороны поверхности; сторону задаёт геометрия. Нормаль
           * ориентируется НА КАМЕРУ — мы видим именно эту сторону, — и если она
           * отвёрнута от солнца, света нет, сколько бы ни было `f`. */
          double vn = 0.0;
          for (int c = 0; c < 3; c++)
            vn += (nn[c] / ln2) * d[c];
          /* СКОЛЬКО ВИДИМЫХ ПИКСЕЛЕЙ ПОКАЗЫВАЮТ ТЫЛЬНУЮ ПО НОРМАЛИ ГРАНЬ.
           * Предложение пользователя 08-07 — считать невидимыми все обратные по
           * нормали поверхности — верно ТОЛЬКО для замкнутой геометрии с
           * согласованными нормалями: у открытого листа тыльная сторона видна, и
           * пометив её, мы пометим видимое, то есть сломаем односторонность.
           * Проверка прямая: если тыльных среди ВИДИМЫХ заметная доля, правило
           * на этих сценах незаконно. */
          if (vn > 0.0) nback++;
          double sgn = (vn > 0.0) ? -1.0 : 1.0; /* нормаль к глазу */
          double cs = 0.0;
          for (int c = 0; c < 3; c++)
            cs -= sgn * (nn[c] / ln2) * X.dir[c]; /* `dir` — направление ЛУЧЕЙ источника */
          if (cs < 0.0) cs = 0.0;                 /* отвёрнута от солнца — темно */
          double kd = m.mtl[m.fm[t]].kd;
          /* §351: при сборке из ячеек сравниваемая доля берётся ИЗ ЯЧЕЙКИ, а не
           * из поточечного обхода — иначе меряется прежний путь, что и вышло
           * трижды подряд. Сторона грани остаётся частью определения величины и
           * применяется к обоим путям одинаково. */
          double fv = cellimg ? fromcell[p] : f;
          pix[p] = kd * fv * cs;
          /* Эталон спрашивает «видно ли отсюда солнце», и ответ у отвёрнутой
           * стороны — НЕТ, независимо от `f`. Сверять надо ту же величину. */
          fpix[p] = (cs > 0.0) ? fv : 0.0;
          sum += pix[p];
          if (fv > 0.5)
            nlit2++;
          else
            nsh2++;
        }
        /* ЭТАЛОН ПО ПИКСЕЛЯМ — НАСКОЛЬКО КАРТИНКА НЕВЕРНА. До сих пор у линии
         * фронта сверялись ЯЧЕЙКИ (§272, выборка по обходу), а картинки не было
         * вовсе; теперь спрашивается прямо: в этом пикселе солнце видно или нет?
         * Эталон — теневой луч из точки попадания перебором ВСЕХ треугольников,
         * то есть заведомо проще проверяемого (правило эталона).
         *
         * ЛУЧ ОДИН, А НЕ ШЕСТНАДЦАТЬ, И ЭТО НЕ ЭКОНОМИЯ. При `nk = 1` фронт
         * несёт видимость ЦЕНТРА диска, то есть точечное солнце; сверять его с
         * усреднением по диску значило бы мерить полутень, которой в этой
         * постановке нет (§274). Диск вернётся вместе с `nk > 1`. */
        if (nref > 0) {
          /* ЭТАЛОН СЧИТАЕТСЯ ПО СЕТКЕ, А НЕ ПО ВЫБОРКЕ ВРАЗБРОС, И ЭТО НЕ ВКУС.
           * Две гипотезы о причине расхождения (разрешение — §286; сдвиг в
           * `face_shift`) обе ОПРОВЕРГНУТЫ замером, а третью гадать нельзя:
           * правило §271 велит смотреть, а не рассуждать. Смотреть можно только
           * КАРТОЙ ошибки, а карта требует сетки. Шаг `estep` пикселей;
           * `ref = N` читается теперь как «сторона карты N x N». */
          double tref2 = now_s();
          int32_t es = bufside / (nref > 0 ? nref : 1);
          if (es < 1) es = 1;
          int32_t ew = bufside / es;
          double *emap = malloc((size_t)ew * (size_t)ew * sizeof *emap);
          if (emap != NULL) {
            int64_t nch = 0, nlitbad = 0, nshbad = 0;
            double sumd = 0.0, maxd = 0.0;
            /* §344: разложение ошибки — без знаковых пикселей и контроль наугад. */
            /* §346: зазор до второго попадания вдоль луча камеры. */
            double *gapr = malloc((size_t)ew * (size_t)ew * sizeof *gapr);
            double *gaps = malloc((size_t)ew * (size_t)ew * sizeof *gaps);
            int64_t ngapr = 0, ngaps = 0, ncop = 0, nopp = 0, nmid = 0;
            int64_t nch_ns = 0, nch_rn = 0;
            double sumd_ns = 0.0, sumd_rn = 0.0;
/* §348: НОВЫЕ СЧЁТЧИКИ ОБЯЗАНЫ ВОЙТИ В РЕДУКЦИЮ. Первая редакция §344/§346
 * писала их из шестнадцати потоков без синхронизации — гонка, и она выдала
 * себя расхождением счётчиков между прогонами (26 554 против 26 615). Гейт
 * гонок не ловит: это дело санитайзера, а не -fanalyzer. */
#pragma omp parallel for schedule(dynamic, 8)                                                      \
    reduction(+ : nch, nlitbad, nshbad, sumd, nch_ns, sumd_ns, nch_rn, sumd_rn, ncop, nopp, nmid)  \
    reduction(max : maxd)
            for (int32_t jj = 0; jj < ew; jj++)
              for (int32_t ii = 0; ii < ew; ii++) {
                size_t p = (size_t)(jj * es) * (size_t)bufside + (size_t)(ii * es);
                emap[(size_t)jj * (size_t)ew + (size_t)ii] = 0.0;
                if (ib[p] < 0) continue;
                double d[3];
                hz_pcull_ray(&C2, (int)(ii * es), (int)(jj * es), d);
                double P[3];
                for (int c = 0; c < 3; c++)
                  P[c] = camo[c] + (double)zb[p] * d[c];
                /* СТОРОНА — ЧАСТЬ ОПРЕДЕЛЯЕМОЙ ВЕЛИЧИНЫ, И ОБЕ СТОРОНЫ СВЕРКИ
                 * ОБЯЗАНЫ СЧИТАТЬ ОДНО И ТО ЖЕ. Сперва проверяемое брало
                 * `|n·ω|`, а эталон — только луч; ошибки шли навстречу и
                 * частично гасили друг друга. Величина определяется так:
                 * «доля диска, видимая из точки С ТОЙ СТОРОНЫ, которую мы
                 * видим». У отвёрнутой от солнца стороны она НОЛЬ, и луч тут ни
                 * при чём — у тонкого треугольника он просто уходит в небо. */
                int32_t te = ib[p];
                const double *Ae = m.v + 3 * (size_t)m.f[3 * (size_t)te + 0];
                const double *Be = m.v + 3 * (size_t)m.f[3 * (size_t)te + 1];
                const double *Ce = m.v + 3 * (size_t)m.f[3 * (size_t)te + 2];
                /* Рёбра выписаны поимённо, а не циклом: `gcc -fanalyzer` не
                 * связывает заполнение массива циклом с его чтением сразу за
                 * циклом и объявляет `w1[1]` неинициализированным. Известный
                 * класс (CLAUDE.md), и лечится он не обходом вокруг, а тем,
                 * чтобы не давать анализатору повода. */
                double w1[3] = {Be[0] - Ae[0], Be[1] - Ae[1], Be[2] - Ae[2]};
                double w2[3] = {Ce[0] - Ae[0], Ce[1] - Ae[1], Ce[2] - Ae[2]};
                double ne[3];
                ne[0] = w1[1] * w2[2] - w1[2] * w2[1];
                ne[1] = w1[2] * w2[0] - w1[0] * w2[2];
                ne[2] = w1[0] * w2[1] - w1[1] * w2[0];
                double vne = 0.0, sne = 0.0;
                for (int c = 0; c < 3; c++)
                  vne += ne[c] * d[c];
                double sg2 = (vne > 0.0) ? -1.0 : 1.0;
                for (int c = 0; c < 3; c++)
                  sne -= sg2 * ne[c] * X.dir[c];
                if (!(sne > 0.0)) {
                  double dd0 = fabs(fpix[p] - 0.0);
                  emap[(size_t)jj * (size_t)ew + (size_t)ii] = dd0;
                  nch++;
                  sumd += dd0;
                  if (dd0 > maxd) maxd = dd0;
                  if (fpix[p] > 0.5) nlitbad++;
                  continue;
                }
                /* Отступ вдоль луча НА СОЛНЦЕ, а не по нормали: нормаль у .obj
                 * может смотреть внутрь, и отступ по ней уводил бы точку под
                 * поверхность ровно в половине случаев. */
                double ofs = 1e-5 * (fabs(P[0]) + fabs(P[1]) + fabs(P[2]) + 1.0);
                double O[3], S[3];
                for (int c = 0; c < 3; c++) {
                  S[c] = -X.dir[c];
                  O[c] = P[c] + ofs * S[c];
                }
                int blk = 0;
                for (int32_t t = 0; t < m.nt && !blk; t++) {
                  const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
                  const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
                  const double *Cc = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
                  /* Поимённо, а не циклом, — тот же класс ложных находок
                   * `-fanalyzer`, что и выше. */
                  double q1[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
                  double q2[3] = {Cc[0] - A[0], Cc[1] - A[1], Cc[2] - A[2]};
                  double pv[3], tv[3], qv[3];
                  pv[0] = S[1] * q2[2] - S[2] * q2[1];
                  pv[1] = S[2] * q2[0] - S[0] * q2[2];
                  pv[2] = S[0] * q2[1] - S[1] * q2[0];
                  double det = q1[0] * pv[0] + q1[1] * pv[1] + q1[2] * pv[2];
                  if (det > -1e-12 && det < 1e-12) continue;
                  double inv = 1.0 / det;
                  tv[0] = O[0] - A[0];
                  tv[1] = O[1] - A[1];
                  tv[2] = O[2] - A[2];
                  double uu = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
                  if (uu < 0.0 || uu > 1.0) continue;
                  qv[0] = tv[1] * q1[2] - tv[2] * q1[1];
                  qv[1] = tv[2] * q1[0] - tv[0] * q1[2];
                  qv[2] = tv[0] * q1[1] - tv[1] * q1[0];
                  double vv = (S[0] * qv[0] + S[1] * qv[1] + S[2] * qv[2]) * inv;
                  if (vv < 0.0 || uu + vv > 1.0) continue;
                  double tt = (q2[0] * qv[0] + q2[1] * qv[1] + q2[2] * qv[2]) * inv;
                  if (tt > 0.0) blk = 1;
                }
                double ftrue = blk ? 0.0 : 1.0;
                double four = fpix[p];
                double dd = fabs(four - ftrue);
                emap[(size_t)jj * (size_t)ew + (size_t)ii] = dd;
                nch++;
                sumd += dd;
                if (dd > maxd) maxd = dd;
                /* ВТОРОЕ ПОПАДАНИЕ ВДОЛЬ ЛУЧА КАМЕРЫ (§346, указание пользователя
                 * 08-09). Считается ЗДЕСЬ, а не в досье, и это исправление
                 * первой редакции: там классификация шла по `fpix > 0.5`, то
                 * есть по «мы даём тень», а не по расхождению с эталоном —
                 * контроль вышел пустым (0 согласных пикселей), что и выдало
                 * поломку. Здесь известны ОБЕ величины, `ftrue` и `four`.
                 *
                 * Что ищем: если сразу за видимой поверхностью стоит вторая,
                 * почти совпадающая, то точка `P` из z-буфера лежит на границе
                 * двух поверхностей, и ячейка фронта под ней перекрыта по
                 * праву — тогда дефект не в операторе, а в постановке сверки. */
                if (gapr != NULL && gaps != NULL) {
                  double t1 = 1e300, t2 = 1e300;
                  int32_t k1 = -1, k2 = -1;
                  for (int32_t tq = 0; tq < m.nt; tq++) {
                    const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)tq + 0];
                    const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)tq + 1];
                    const double *Cc = m.v + 3 * (size_t)m.f[3 * (size_t)tq + 2];
                    double q1[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
                    double q2[3] = {Cc[0] - A[0], Cc[1] - A[1], Cc[2] - A[2]};
                    double pv[3], tv[3], qv[3];
                    pv[0] = d[1] * q2[2] - d[2] * q2[1];
                    pv[1] = d[2] * q2[0] - d[0] * q2[2];
                    pv[2] = d[0] * q2[1] - d[1] * q2[0];
                    double det = q1[0] * pv[0] + q1[1] * pv[1] + q1[2] * pv[2];
                    if (det > -1e-12 && det < 1e-12) continue;
                    double inv = 1.0 / det;
                    tv[0] = camo[0] - A[0];
                    tv[1] = camo[1] - A[1];
                    tv[2] = camo[2] - A[2];
                    double uu2 = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
                    if (uu2 < 0.0 || uu2 > 1.0) continue;
                    qv[0] = tv[1] * q1[2] - tv[2] * q1[1];
                    qv[1] = tv[2] * q1[0] - tv[0] * q1[2];
                    qv[2] = tv[0] * q1[1] - tv[1] * q1[0];
                    double vv2 = (d[0] * qv[0] + d[1] * qv[1] + d[2] * qv[2]) * inv;
                    if (vv2 < 0.0 || uu2 + vv2 > 1.0) continue;
                    double tt2 = (q2[0] * qv[0] + q2[1] * qv[1] + q2[2] * qv[2]) * inv;
                    if (!(tt2 > 0.0)) continue;
                    if (tt2 < t1) {
                      t2 = t1;
                      k2 = k1;
                      t1 = tt2;
                      k1 = tq;
                    } else if (tt2 < t2) {
                      t2 = tt2;
                      k2 = tq;
                    }
                  }
                  if (k2 >= 0) {
                    double gap = t2 - t1;
                    if (ftrue > 0.5 && four < 0.5) {
/* Массивы — под критической секцией: популяция мала (сотня), цена никакая. */
#pragma omp critical
                      gapr[ngapr++] = gap;
                      /* Явные нули: `cppcheck` не прослеживает заполнение через
                       * условный указатель и ГАТИТ по нему — тот же класс, что
                       * у `gcc -fanalyzer` с обёртками (CLAUDE.md). Правится
                       * формой записи, а не подавлением. */
                      double nrm2[2][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
                      for (int qq = 0; qq < 2; qq++) {
                        int32_t tq = qq ? k2 : k1;
                        const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)tq + 0];
                        const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)tq + 1];
                        const double *Cc = m.v + 3 * (size_t)m.f[3 * (size_t)tq + 2];
                        double e1[3], e2[3], cr[3];
                        for (int c = 0; c < 3; c++) {
                          e1[c] = B[c] - A[c];
                          e2[c] = Cc[c] - A[c];
                        }
                        cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
                        cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
                        cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
                        double l2 = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
                        if (!(l2 > 0.0)) l2 = 1.0;
                        for (int c = 0; c < 3; c++)
                          nrm2[qq][c] = cr[c] / l2;
                      }
                      double cs = nrm2[0][0] * nrm2[1][0] + nrm2[0][1] * nrm2[1][1] +
                                  nrm2[0][2] * nrm2[1][2];
                      if (cs > 0.9)
                        ncop++;
                      else if (cs < -0.9)
                        nopp++;
                      else
                        nmid++;
                    } else if (dd < 1e-6) {
#pragma omp critical
                      gaps[ngaps++] = gap;
                    }
                  }
                }
                if (ftrue < 0.5 && four > 0.5) nlitbad++;
                if (ftrue > 0.5 && four < 0.5) nshbad++;
                /* §344: та же ошибка БЕЗ пикселей, разошедшихся знаком, и
                 * контроль — БЕЗ такого же числа пикселей, взятых наугад
                 * (детерминированно, хешем от номера: замер обязан
                 * повторяться). Оба накопителя идут тем же проходом. */
                int issign = (ftrue < 0.5 && four > 0.5) || (ftrue > 0.5 && four < 0.5);
                if (!issign) {
                  sumd_ns += dd;
                  nch_ns++;
                }
                uint32_t hh = (uint32_t)(jj * 73856093 ^ ii * 19349663);
                hh ^= hh >> 13;
                hh *= 2654435761u;
                if ((hh % 1000u) >=
                    4u) { /* выкинуть примерно 0.4 % — столько же, сколько знаковых */
                  sumd_rn += dd;
                  nch_rn++;
                }
              }
            if (gapr != NULL && gaps != NULL) {
              qsort(gapr, (size_t)ngapr, sizeof *gapr, cmp_dbl_psun);
              qsort(gaps, (size_t)ngaps, sizeof *gaps, cmp_dbl_psun);
              printf("     ВТОРОЕ ПОПАДАНИЕ (§346): ОБРАТНЫХ %lld, зазор p50 %.6f p90 %.6f мин "
                     "%.6f м; "
                     "нормали пары: сонаправленных %lld, встречных %lld, прочих %lld\n",
                     (long long)ngapr, ngapr ? gapr[ngapr / 2] : 0.0,
                     ngapr ? gapr[(ngapr * 9) / 10] : 0.0, ngapr ? gapr[0] : 0.0, (long long)ncop,
                     (long long)nopp, (long long)nmid);
              printf("     КОНТРОЛЬ — СОГЛАСНЫЕ: %lld, зазор p50 %.6f p90 %.6f м\n",
                     (long long)ngaps, ngaps ? gaps[ngaps / 2] : 0.0,
                     ngaps ? gaps[(ngaps * 9) / 10] : 0.0);
            }
            free(gapr);
            free(gaps);
            printf("     БЕЗ ЗНАКОВЫХ ПИКСЕЛЕЙ (§344): |Δ| среднее %.6f по %lld пикселям; "
                   "КОНТРОЛЬ (столько же наугад): %.6f по %lld\n",
                   nch_ns > 0 ? sumd_ns / (double)nch_ns : 0.0, (long long)nch_ns,
                   nch_rn > 0 ? sumd_rn / (double)nch_rn : 0.0, (long long)nch_rn);
            printf("   ЭТАЛОН ПО СЕТКЕ %dx%d (теневой луч перебором, %lld пикселей за %.1f с):\n"
                   "     |Δдоли| среднее %.4f, наибольшее %.4f; РАЗОШЛИСЬ ЗНАКОМ: у нас светло "
                   "а в тени %lld (%.2f %%), у нас тень а на свету %lld (%.2f %%)\n",
                   ew, ew, (long long)nch, now_s() - tref2, nch > 0 ? sumd / (double)nch : 0.0,
                   maxd, (long long)nlitbad, nch > 0 ? 100.0 * (double)nlitbad / (double)nch : 0.0,
                   (long long)nshbad, nch > 0 ? 100.0 * (double)nshbad / (double)nch : 0.0);
            {
              char ep[256];
              const char *b2 = strrchr(argv[1], '/');
              snprintf(ep, sizeof ep, "img/psun_%s_err.ppm", b2 ? b2 + 1 : argv[1]);
              hz_ppm_write(ep, emap, ew, ew);
              printf("     КАРТА ОШИБКИ: %s (ярко — там, где мы расходимся с эталоном)\n", ep);
            }
          }
          free(emap);
        }
        /* ДОСЬЕ НА РАСХОДЯЩИЕСЯ ПИКСЕЛИ (ключ `dbg=N`, §293). Три названные
         * мною причины подряд оказались не теми, и каждый раз это выяснялось
         * замером. Значит гадать нельзя: надо взять конкретные пиксели и
         * посмотреть, ЧТО в них происходит. Главный вопрос к каждому — ГДЕ
         * ЗАСЛОНИТЕЛЬ: если он внутри той же ячейки, промахнулось перекрытие;
         * если выше по пути — свет протёк сквозь другую ячейку, и виновата
         * передача состояния, а не оператор. */
        if (ndbg != 0) {
          printf("== ДОСЬЕ НА %d ПИКСЕЛЕЙ (%s)\n", ndbg < 0 ? -ndbg : ndbg,
                 ndbg < 0 ? "ОБРАТНАЯ популяция: мы даём ТЕНЬ, луч даёт СВЕТ"
                          : "мы даём свет, луч даёт тень");
          int shown = 0;
          /* `dbg > 0` — прямая популяция («мы светло, луч тень»); `dbg < 0` —
           * ОБРАТНАЯ («мы тень, луч светло»), та самая, что не сдвинулась ни от
           * чего за весь день. */
          int rev = (ndbg < 0);
          int nwant = rev ? -ndbg : ndbg;
          for (size_t p = 0; p < np && shown < nwant; p += 97) {
            if (ib[p] < 0) continue;
            if (!rev && !(fpix[p] > 0.5)) continue;
            if (rev && (fpix[p] > 0.5)) continue;
            double d[3];
            hz_pcull_ray(&C2, (int)(p % (size_t)bufside), (int)(p / (size_t)bufside), d);
            double P[3];
            for (int c = 0; c < 3; c++)
              P[c] = camo[c] + (double)zb[p] * d[c];
            double ofs = 1e-5 * (fabs(P[0]) + fabs(P[1]) + fabs(P[2]) + 1.0);
            double O[3], S[3];
            for (int c = 0; c < 3; c++) {
              S[c] = -X.dir[c];
              O[c] = P[c] + ofs * S[c];
            }
            int32_t blk = -1;
            double blkt = 0.0;
            for (int32_t t = 0; t < m.nt && blk < 0; t++) {
              const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
              const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
              const double *Cc = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
              double q1[3] = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
              double q2[3] = {Cc[0] - A[0], Cc[1] - A[1], Cc[2] - A[2]};
              double pv[3], tv[3], qv[3];
              pv[0] = S[1] * q2[2] - S[2] * q2[1];
              pv[1] = S[2] * q2[0] - S[0] * q2[2];
              pv[2] = S[0] * q2[1] - S[1] * q2[0];
              double det = q1[0] * pv[0] + q1[1] * pv[1] + q1[2] * pv[2];
              if (det > -1e-12 && det < 1e-12) continue;
              double inv = 1.0 / det;
              tv[0] = O[0] - A[0];
              tv[1] = O[1] - A[1];
              tv[2] = O[2] - A[2];
              double uu = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
              if (uu < 0.0 || uu > 1.0) continue;
              qv[0] = tv[1] * q1[2] - tv[2] * q1[1];
              qv[1] = tv[2] * q1[0] - tv[0] * q1[2];
              qv[2] = tv[0] * q1[1] - tv[1] * q1[0];
              double vv = (S[0] * qv[0] + S[1] * qv[1] + S[2] * qv[2]) * inv;
              if (vv < 0.0 || uu + vv > 1.0) continue;
              double tt = (q2[0] * qv[0] + q2[1] * qv[1] + q2[2] * qv[2]) * inv;
              if (tt > 0.0) {
                blk = t;
                blkt = tt;
              }
            }
            if (!rev && blk < 0) continue; /* согласие, не наш случай */
            if (rev && blk >= 0) continue; /* тень и у нас, и у эталона */
            /* Ячейка фронта, накрывшая точку. */
            int32_t nid = 0, lev = 0;
            while (X.cellf[nid] < 0.0f && T.nd[nid].child >= 0) {
              double mid[3];
              int k = 0;
              for (int c = 0; c < 3; c++) {
                mid[c] = 0.5 * (T.nd[nid].lo[c] + T.nd[nid].hi[c]);
                if (P[c] >= mid[c]) k |= (1 << c);
              }
              nid = T.nd[nid].child + k;
              lev++;
            }
            double hcell = T.nd[nid].hi[0] - T.nd[nid].lo[0];
            /* Где заслонитель: внутри этой ячейки или выше по лучу? */
            double Q[3];
            for (int c = 0; c < 3; c++)
              Q[c] = O[c] + blkt * S[c];
            int inside = 1;
            for (int c = 0; c < 3 && inside; c++)
              if (Q[c] < T.nd[nid].lo[c] || Q[c] > T.nd[nid].hi[c]) inside = 0;
            /* Размер заслонителя против размера ячейки. В ОБРАТНОЙ популяции
             * заслонителя нет вовсе (`blk = -1`), и читать его нельзя — это и
             * был SIGSEGV при первом прогоне. */
            double area = 0.0;
            if (blk >= 0) {
              const double *A2 = m.v + 3 * (size_t)m.f[3 * (size_t)blk + 0];
              const double *B2 = m.v + 3 * (size_t)m.f[3 * (size_t)blk + 1];
              const double *C2v = m.v + 3 * (size_t)m.f[3 * (size_t)blk + 2];
              double u1[3] = {B2[0] - A2[0], B2[1] - A2[1], B2[2] - A2[2]};
              double u2[3] = {C2v[0] - A2[0], C2v[1] - A2[1], C2v[2] - A2[2]};
              double cr[3];
              cr[0] = u1[1] * u2[2] - u1[2] * u2[1];
              cr[1] = u1[2] * u2[0] - u1[0] * u2[2];
              cr[2] = u1[0] * u2[1] - u1[1] * u2[0];
              area = 0.5 * sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
            }
            if (rev) {
              /* Заслонителя нет вовсе — печатать про него нечего; печатается то,
               * что могло погасить свет у НАС: состояние ячейки и её вид. */
              /* СКОЛЬЗЯЩАЯ ЛИ ПОВЕРХНОСТЬ. Подозрение: свет гасится МНОГОКРАТНО
               * одной и той же стеной — почти параллельная лучу поверхность
               * попадает в десятки ячеек подряд, и каждая применяет своё
               * `(1 − c)`. Тогда `|n·ω|` у этих пикселей обязан быть мал. */
              int32_t tv2 = ib[p];
              const double *Av = m.v + 3 * (size_t)m.f[3 * (size_t)tv2 + 0];
              const double *Bv = m.v + 3 * (size_t)m.f[3 * (size_t)tv2 + 1];
              const double *Cv2 = m.v + 3 * (size_t)m.f[3 * (size_t)tv2 + 2];
              double g1[3] = {Bv[0] - Av[0], Bv[1] - Av[1], Bv[2] - Av[2]};
              double g2[3] = {Cv2[0] - Av[0], Cv2[1] - Av[1], Cv2[2] - Av[2]};
              double gn[3];
              gn[0] = g1[1] * g2[2] - g1[2] * g2[1];
              gn[1] = g1[2] * g2[0] - g1[0] * g2[2];
              gn[2] = g1[0] * g2[1] - g1[1] * g2[0];
              double gl = sqrt(gn[0] * gn[0] + gn[1] * gn[1] + gn[2] * gn[2]);
              if (!(gl > 0.0)) gl = 1.0;
              double cw = 0.0;
              for (int c = 0; c < 3; c++)
                cw += (gn[c] / gl) * X.dir[c];
              if (cw < 0.0) cw = -cw;
              printf("   пиксель %6d: f=%.3f  ЯЧЕЙКА ур.%2d ребро %.4f м  занята=%d  "
                     "глубина в ячейке=%s  |n·ω| = %.4f\n",
                     (int)p, fpix[p], lev, hcell, (int)occ[nid],
                     (X.cellslot != NULL && X.cellslot[nid] >= 0) ? "есть" : "нет", cw);
              /* ПРОСЛЕЖИВАНИЕ ЛУЧА ПО ЯЧЕЙКАМ (§303). Ноль приходит СВЕРХУ, и
               * две гипотезы о том, почему, опровергнуты (§302). Значит надо не
               * гадать, а пройти путь от точки НАЗАД К ИСТОЧНИКУ и найти ячейку,
               * в которой свет погас: печатается место перехода `f > 0 → f = 0`
               * вместе с уровнем и занятостью этой ячейки. */
              {
                double dg2 = 0.0;
                for (int c = 0; c < 3; c++) {
                  double sd2 = T.nd[0].hi[c] - T.nd[0].lo[c];
                  dg2 += sd2 * sd2;
                }
                dg2 = sqrt(dg2);
                int32_t prev = -1;
                int nprof = 0;
                for (int st2 = 1; st2 <= 400; st2++) {
                  double Q2[3];
                  double back = dg2 * (double)st2 / 400.0;
                  for (int c = 0; c < 3; c++)
                    Q2[c] = P[c] - back * X.dir[c];
                  int out2 = 0;
                  for (int c = 0; c < 3 && !out2; c++)
                    if (Q2[c] < T.nd[0].lo[c] || Q2[c] > T.nd[0].hi[c]) out2 = 1;
                  if (out2) break;
                  int32_t q3 = 0, lv2 = 0;
                  while (X.cellf[q3] < 0.0f && T.nd[q3].child >= 0) {
                    double md[3];
                    int kk = 0;
                    for (int c = 0; c < 3; c++) {
                      md[c] = 0.5 * (T.nd[q3].lo[c] + T.nd[q3].hi[c]);
                      if (Q2[c] >= md[c]) kk |= (1 << c);
                    }
                    q3 = T.nd[q3].child + kk;
                    lv2++;
                  }
                  if (q3 == prev) continue;
                  double fq = (X.cellf[q3] >= 0.0f) ? (double)X.cellf[q3] : -1.0;
                  /* ПРОФИЛЬ ПО ПУТИ, а не только точка гашения: §307 показал,
                   * что свет входит в ячейку уже ослабленным, значит гасит не
                   * одна ячейка, а затухание НАКАПЛИВАЕТСЯ. Печатаются первые
                   * ячейки от точки назад к источнику вместе с их долей. */
                  if (nprof < 10) {
                    /* КЛЕТКА СЕТКИ, через которую идёт луч, — та же величина,
                     * что читает картинка. Среднее по грани здесь обманывает:
                     * грань может быть тёмной в среднем и светлой в нужной
                     * клетке (О76). */
                    double fcell = -1.0;
                    if (X.cellval != NULL && X.cellslot[q3] >= 0) {
                      int a4 = X.celldepax;
                      int p4 = (a4 + 1) % 3, q4 = (a4 + 2) % 3;
                      const hz_ptnode *N4 = &T.nd[q3];
                      double h4p = N4->hi[p4] - N4->lo[p4], h4q = N4->hi[q4] - N4->lo[q4];
                      double fa4 = (X.dir[a4] > 0.0) ? N4->lo[a4] : N4->hi[a4];
                      double t4 = (Q2[a4] - fa4) / X.dir[a4];
                      double u4 = Q2[p4] - t4 * X.dir[p4], v4 = Q2[q4] - t4 * X.dir[q4];
                      int gi4 = (int)((double)HZ_PFRONT_COVN * (u4 - N4->lo[p4]) / h4p);
                      int gj4 = (int)((double)HZ_PFRONT_COVN * (v4 - N4->lo[q4]) / h4q);
                      if (gi4 < 0) gi4 = 0;
                      if (gj4 < 0) gj4 = 0;
                      if (gi4 > HZ_PFRONT_COVN - 1) gi4 = HZ_PFRONT_COVN - 1;
                      if (gj4 > HZ_PFRONT_COVN - 1) gj4 = HZ_PFRONT_COVN - 1;
                      fcell = (double)X.cellval[HZ_PFRONT_NCELL * (size_t)X.cellslot[q3] +
                                                (size_t)(gj4 * HZ_PFRONT_COVN + gi4)];
                    }
                    printf("        %5.3f м назад: ур.%2d ребро %.4f занята=%d  среднее=%.4f  "
                           "КЛЕТКА=%.4f\n",
                           back, lv2, T.nd[q3].hi[0] - T.nd[q3].lo[0], (int)occ[q3], fq, fcell);
                    nprof++;
                  } else
                    break;
                  prev = q3;
                }
              }
            } else
              printf("   пиксель %6d: f=%.3f  ЯЧЕЙКА ур.%2d ребро %.4f м  ЗАСЛОНИТЕЛЬ на %.4f м "
                     "%s, площадь %.3e м² (ребро ~%.4f м)\n",
                     (int)p, fpix[p], lev, hcell, blkt, inside ? "ВНУТРИ ячейки" : "ВЫШЕ по лучу",
                     area, sqrt(2.0 * area));
            shown++;
          }
        }
        char path[256];
        const char *base = strrchr(argv[1], '/');
        snprintf(path, sizeof path, "img/psun_%s_mark%d.ppm", base ? base + 1 : argv[1], usemark);
        hz_ppm_write(path, pix, bufside, bufside);
        /* PFM рядом с PPM — чтобы сравнивать две картинки ЧИСЛАМИ (`pfmdiff`), а
         * не глазами: приёмка огрубления невидимого есть «видимое не изменилось»,
         * и глаз на такой вопрос не отвечает. */
        snprintf(path, sizeof path, "img/psun_%s_mark%d.pfm", base ? base + 1 : argv[1], usemark);
        hz_pfm_write(path, pix, pix, pix, bufside, bufside);
        printf("== КАРТИНКА (проход 3, сбор по пикселю) за %.2f с: %s\n"
               "   ПИКСЕЛЕЙ БЕЗ ЗАПИСИ ФРОНТА: %lld\n"
               "   ТЫЛЬНЫХ ПО НОРМАЛИ СРЕДИ ВИДИМЫХ: %lld из %lld (%.2f %%) — столько показало "
               "бы невидимым правило «обратные по нормали»\n"
               "   пикселей неба %lld, освещённых %lld, в тени %lld; средняя яркость %.6f\n",
               now_s() - timg, path, (long long)nmiss, (long long)nback,
               (long long)(np - (size_t)nsky),
               (np - (size_t)nsky) > 0 ? 100.0 * (double)nback / (double)(np - (size_t)nsky) : 0.0,
               (long long)nsky, (long long)nlit2, (long long)nsh2, sum / (double)np);
        hz_pcull_free(&C2);
      }
    }
    free(zb);
    free(ib);
    free(pix);
    free(fromcell);
    free(fpix);
  }

  /* ПРИЁМКА ЭТОЙ ЛИНИИ ОПТИМИЗАЦИЙ — ОДНОСТОРОННЯЯ (указание пользователя 08-06):
   * «нужно лишь, чтобы мы не огрубили ВИДИМЫЕ части; вполне допустимо, что часть
   * невидимых останется неогрублённой». Значит проверять надо НЕ расхождение с
   * эталоном по яркости — то физика тени, а у нас тень в смысле НЕВИДИМОСТИ, —
   * а одно: нет ли среди помеченных ячейки, до которой камера ДОСТАЁТ.
   * Проверка независима от фронта: прямой луч от глаза к центру ячейки перебором
   * треугольников. Нарушений обязано быть НОЛЬ. */
  if (usemark && mk != NULL && nref > 0) {
    int32_t nchk = 0, nbadv = 0, nbadbox = 0;
    double tchk = now_s();
    /* РАСПРЕДЕЛЕНИЕ НАРУШЕНИЙ ПО УГЛУ ОТ ОСИ ВЗГЛЯДА — РАЗЛИЧИТЕЛЬ ПРИЧИН.
     * Параллельное приближение камеры ошибается тем сильнее, чем дальше от оси:
     * у края кадра расхождение равно половине поля. Если нарушения сидят у края
     * — виновато оно; если разбросаны по всему полю (и есть у самой оси) —
     * причина другая, и параллельностью её не объяснить. Четыре корзины по
     * долям полуполя. */
    int32_t bin_chk[4] = {0, 0, 0, 0}, bin_bad[4] = {0, 0, 0, 0};
    double halffov = 0.5 * HZ_CFG_FOV_DEG * 3.14159265358979323846 / 180.0;
    /* Перебирать надо ПОМЕЧЕННЫЕ, а не сетку по всем узлам: пометок 47 тыс. из
     * 34 млн, и шаг по всем узлам не попал ни в одну (проверено 0). Считаем их
     * сперва, потом берём каждую `stride`-ю ИЗ НИХ. */
    int32_t nmk = 0;
    for (int32_t i = 0; i < T.nnd; i++)
      if (mk[i] != 255) nmk++;
    int32_t stride = (nmk > nref) ? nmk / nref : 1;
    int32_t seen = 0;
    for (int32_t i = 0; i < T.nnd; i++) {
      if (mk[i] == 255) continue;
      if ((seen++ % stride) != 0) continue;
      double P0[3], d2[3], len = 0.0;
      for (int c = 0; c < 3; c++) {
        P0[c] = 0.5 * (T.nd[i].lo[c] + T.nd[i].hi[c]);
        d2[c] = P0[c] - camo[c];
        len += d2[c] * d2[c];
      }
      len = sqrt(len);
      if (!(len > 0.0)) continue;
      for (int c = 0; c < 3; c++)
        d2[c] /= len;
      /* ВНУТРИ ли пирамиды — И ЗДЕСЬ ВАЖНО, ЧТО ИМЕННО ПРОВЕРЯЕТСЯ. Луч пускается
       * в ЦЕНТР ячейки, значит и на видимость проверять надо ЭТУ ТОЧКУ, а не
       * коробку: точка вне кадра лучом достижима, но не видима, и нарушением не
       * является. Прежде проверялась КОРОБКА (перекрывает ли она пирамиду), и
       * пара «коробка × луч в центр» была рассогласована.
       *
       * ЧЕСТНОСТЬ ТРЕБУЕТ ПЕЧАТАТЬ ОБА ЧИСЛА, потому что правило пометки в этом
       * же шаге стало отвечать «не видно» и за краем кадра: если бы приёмка
       * смягчалась вместе с правилом, она перестала бы что-либо проверять.
       * Поэтому считаются нарушения ПО ТОЧКЕ (строгая пара) и ПО КОРОБКЕ
       * (прежний, более широкий охват), и печатаются рядом. */
      int inbox = 1, inpt = 1;
      for (int kf = 0; kf < 6; kf++) {
        const double *PL = X.fr[kf];
        double sf = PL[3], sp = PL[3];
        for (int c = 0; c < 3; c++) {
          sf += (PL[c] > 0.0 ? T.nd[i].hi[c] : T.nd[i].lo[c]) * PL[c];
          sp += P0[c] * PL[c];
        }
        if (sf < 0.0) inbox = 0;
        if (sp < 0.0) inpt = 0;
      }
      if (!inbox) continue;
      nchk++;
      int blocked = 0;
      for (int32_t t = 0; t < m.nt && !blocked; t++) {
        const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
        const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
        const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
        double e1[3], e2[3], pv[3], tv[3], qv[3];
        for (int c = 0; c < 3; c++) {
          e1[c] = B[c] - A[c];
          e2[c] = C[c] - A[c];
        }
        pv[0] = d2[1] * e2[2] - d2[2] * e2[1];
        pv[1] = d2[2] * e2[0] - d2[0] * e2[2];
        pv[2] = d2[0] * e2[1] - d2[1] * e2[0];
        double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
        if (det > -1e-12 && det < 1e-12) continue;
        double inv = 1.0 / det;
        for (int c = 0; c < 3; c++)
          tv[c] = camo[c] - A[c];
        double uu = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
        if (uu < 0.0 || uu > 1.0) continue;
        qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
        qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
        qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
        double vv = (d2[0] * qv[0] + d2[1] * qv[1] + d2[2] * qv[2]) * inv;
        if (vv < 0.0 || uu + vv > 1.0) continue;
        double tt = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
        if (tt > 1e-6 && tt < len - 1e-6) blocked = 1;
      }
      double ca = 0.0;
      for (int c = 0; c < 3; c++)
        ca += d2[c] * camf[c];
      if (ca > 1.0) ca = 1.0;
      if (ca < -1.0) ca = -1.0;
      int bin = (int)(4.0 * acos(ca) / halffov);
      if (bin < 0) bin = 0;
      if (bin > 3) bin = 3;
      bin_chk[bin]++;
      if (!blocked) {
        nbadbox++;
        if (inpt) {
          nbadv++;
          bin_bad[bin]++;
        }
      }
    }
    printf("   ПРИЁМКА: помеченных проверено %d за %.1f с; ВИДИМЫХ СРЕДИ НИХ %d (обязано 0); "
           "по прежнему, более широкому критерию коробки — %d\n",
           nchk, now_s() - tchk, nbadv, nbadbox);
    printf("   по углу от оси взгляда (доля полуполя): ");
    for (int b = 0; b < 4; b++)
      printf("%d/4 %d из %d | ", b + 1, bin_bad[b], bin_chk[b]);
    printf("\n");
  }

  /* ЭТАЛОН ЛУЧАМИ (А519). Считает ТО ЖЕ — долю открытого диска источника, — но
   * другим способом: прямым перебором треугольников. Перебор выбран сознательно:
   * эталон обязан быть ПРОЩЕ проверяемого, иначе сверяются две одинаковые
   * ошибки. Отсюда и ограничение по сцене — на зале это 124 тыс. треугольников
   * на луч, и дальше растёт линейно.
   * Диск источника выбирается спиралью Ферма: ДЕТЕРМИНИРОВАННО, без случайности,
   * чтобы повтор давал те же числа. */
  /* ЭТАЛОН ЛУЧАМИ И ПРИЁМКА ПОМЕТОК — РАЗНЫЕ ПРОВЕРКИ, И ЦЕНА У НИХ РАЗНАЯ.
   * Эталон стреляет `8 × 8` лучей на выборочную ячейку перебором ВСЕХ
   * треугольников: на зале это две секунды, на Сан-Мигеле — `600 × 64 × 6.7e6`,
   * то есть десятки минут (замерено 08-06 прогоном, который пришлось прервать).
   * Приёмка невидимости стоит один луч на ячейку и ни от чего этого не зависит.
   * Поэтому эталон включается отдельным ключом `eta=1`. */
  if (eta && nref > 0 && X.refn > 0) {
    double half = asun; /* угловой радиус источника; у солнца 0.5·9.3e-3 (§241.4) */
    double t1[3] = {0.0, 0.0, 1.0}, r1[3], r2[3];
    if (fabs(X.dir[2]) > 0.9) {
      t1[0] = 1.0;
      t1[2] = 0.0;
    }
    r1[0] = X.dir[1] * t1[2] - X.dir[2] * t1[1];
    r1[1] = X.dir[2] * t1[0] - X.dir[0] * t1[2];
    r1[2] = X.dir[0] * t1[1] - X.dir[1] * t1[0];
    double ln2 = sqrt(r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
    if (!(ln2 > 0.0)) ln2 = 1.0;
    for (int c = 0; c < 3; c++)
      r1[c] /= ln2;
    r2[0] = X.dir[1] * r1[2] - X.dir[2] * r1[1];
    r2[1] = X.dir[2] * r1[0] - X.dir[0] * r1[2];
    r2[2] = X.dir[0] * r1[1] - X.dir[1] * r1[0];
    /* РАЗДЕЛЕНИЕ ПО СОСТОЯНИЮ (разбор 08-06). Две причины расхождения известны
     * заранее и лечатся по-разному: перекрытие по ОДНОМУ куску (А517) недооценит
     * тень ВЕЗДЕ, а линейность состояния (Р1) промахнётся только на КРАЮ тени.
     * Поэтому ошибка считается порознь в трёх корзинах: глубокая тень, открытый
     * свет и край. Если промахи только на краю — виновата линейность и это цена
     * представления; если и в тени — виновато перекрытие. */
    double sum3[3] = {0.0, 0.0, 0.0}, wor3[3] = {0.0, 0.0, 0.0};
    int64_t cnt3[3] = {0, 0, 0}, bad3[3] = {0, 0, 0};
    double sum = 0.0, worst = 0.0, tref = now_s();
    int64_t nbad = 0;
    for (int32_t s = 0; s < X.refn; s++) {
      /* СВЕРЯТЬ НАДО СРЕДНЕЕ СО СРЕДНИМ (разбор 08-06). Наше состояние есть
       * ЛИНЕЙНАЯ функция по грани, и `c0` — её среднее; эталон же стрелял из
       * ЦЕНТРА, то есть сравнивал точку с проекцией. На краю тени это давало
       * расхождение 0.248 при 100 %% ячеек вне приёмки — артефакт сверки, а не
       * ошибка фронта. Теперь эталон берёт восемь точек по ячейке (углы
       * полуразмера) и усредняет: величина становится той же, что у нас. */
      const double *Pc = refpt + 3 * (size_t)s;
      double hh = (refh != NULL) ? 0.5 * refh[s] : 0.0;
      int open = 0, ntot = 0;
      for (int pt = 0; pt < 8; pt++) {
        double P0[3];
        for (int c = 0; c < 3; c++)
          P0[c] = Pc[c] + ((pt & (1 << c)) ? hh : -hh);
        for (int k = 0; k < HZ_REF_RAYS / 2; k++) {
          ntot++;
          double rr = half * sqrt(((double)k + 0.5) / (double)(HZ_REF_RAYS / 2));
          double ph = 2.39996322972865332 * (double)k;
          double d2[3];
          for (int c = 0; c < 3; c++)
            d2[c] = -X.dir[c] + rr * (cos(ph) * r1[c] + sin(ph) * r2[c]);
          int hit = 0;
          for (int32_t t = 0; t < m.nt && !hit; t++) {
            const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
            const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
            const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
            double e1[3], e2[3], pv[3], tv[3], qv[3];
            for (int c = 0; c < 3; c++) {
              e1[c] = B[c] - A[c];
              e2[c] = C[c] - A[c];
            }
            pv[0] = d2[1] * e2[2] - d2[2] * e2[1];
            pv[1] = d2[2] * e2[0] - d2[0] * e2[2];
            pv[2] = d2[0] * e2[1] - d2[1] * e2[0];
            double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
            if (det > -1e-12 && det < 1e-12) continue;
            double inv = 1.0 / det;
            for (int c = 0; c < 3; c++)
              tv[c] = P0[c] - A[c];
            double uu = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
            if (uu < 0.0 || uu > 1.0) continue;
            qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
            qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
            qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
            double vv = (d2[0] * qv[0] + d2[1] * qv[1] + d2[2] * qv[2]) * inv;
            if (vv < 0.0 || uu + vv > 1.0) continue;
            double tt = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
            if (tt > 1e-6) hit = 1;
          }
          if (!hit) open++;
        }
      }
      double vref = (ntot > 0) ? (double)open / (double)ntot : 0.0;
      double d = refvis[s] - vref;
      if (d < 0.0) d = -d;
      sum += d;
      if (d > worst) worst = d;
      if (d > 0.02) nbad++;
      int b = (refvis[s] < 0.05) ? 0 : ((refvis[s] > 0.95) ? 1 : 2);
      cnt3[b]++;
      sum3[b] += d;
      if (d > wor3[b]) wor3[b] = d;
      if (d > 0.02) bad3[b]++;
    }
    printf("   ЭТАЛОН ЛУЧАМИ (%d лучей по диску, перебор треугольников): выборка %d ячеек за "
           "%.1f с\n",
           HZ_REF_RAYS, X.refn, now_s() - tref);
    printf("   расхождение среднее %.4f, наибольшее %.4f; выше приёмки 2 %% — %lld ячеек "
           "(%.2f %%)\n",
           sum / (double)X.refn, worst, (long long)nbad, 100.0 * (double)nbad / (double)X.refn);
    {
      static const char *nm[3] = {"глубокая тень", "открытый свет", "КРАЙ тени"};
      for (int b = 0; b < 3; b++)
        if (cnt3[b] > 0)
          printf("     %-14s ячеек %4lld, среднее %.4f, наибольшее %.4f, выше 2 %%%% — %lld "
                 "(%.1f %%%%)\n",
                 nm[b], (long long)cnt3[b], sum3[b] / (double)cnt3[b], wor3[b], (long long)bad3[b],
                 100.0 * (double)bad3[b] / (double)cnt3[b]);
    }
  }
  free(refpt);
  free(refh);
  free(refvis);
  /* Пометки и высоты поддеревьев жили до конца обхода; санитайзер 08-06
   * показал, что их не освобождали ВООБЩЕ — при выходе безвредно, но течь на
   * гейте видна, и починить дешевле, чем оговаривать. */
  free(dep);
  free(mk);
  free(X.cellf);
  free(vis);
  free(X.cellslot);
  free(X.celldep);
  free(X.cellval);
  free(tlo);
  free(thi);
  free(tpl);
  free(list);
  free(occ);
  hz_ptree_free(&T);
  hz_obj_free(&m);
  return 0;
}
