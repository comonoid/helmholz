/* pflow — ЗАМЕР «ТОЛЬКО ИТЕРАЦИИ»: сколько стоит перенос без оператора вовсе.
 *
 * ВОПРОС ПОЛЬЗОВАТЕЛЯ 08-02, дословно: «какой результат из того, что обходимся
 * ТОЛЬКО итерациями матрицы». То есть: ни таблицы связей, ни форм-факторов, ни
 * приёмников — свет протаскивается через всю сцену вдоль направлений, итерация
 * повторяется, угловое усреднение делает сама итерация.
 *
 * ЧТО СРАВНИВАЕТСЯ. Эталон уже посчитан хранимым оператором на том же городе
 * (§120): поток `Σ B·A` по уровню 0 = `2.398505e+05` Вт/ср при альбедо `0.5`,
 * небе `L = 1.0`, сборке `3102.6` с. Здесь считается то же самое той же физикой
 * и на тех же полигонах, но без сборки. Расхождение потока и есть ответ.
 *
 * ПОЧЕМУ НИЧЕГО НЕ ПИШЕТСЯ ЗАНОВО. Машинерия для этого написана давно и брошена
 * решением §85 (отозвано §121): `prast.c` — параллельная проекция со списком
 * фрагментов на пиксель, `psweep.c` — редукция «излучатель → приёмник» по парам
 * соседей вдоль луча и подгонка линейного `E` по моментам полигона. Здесь только
 * стенд: собрать сцену, погонять отскоки, напечатать числа.
 *
 * ОДНО ОТЛИЧИЕ ОТ ЭТАЛОНА, И ЕГО НАДО НАЗВАТЬ. Эталон решал на узлах СРЕЗА
 * (`152 618`), здесь решается на исходных полигонах (`446 013`): `psweep` живёт
 * на полигонах и уровней не знает. Поэтому сравнивается ПОТОК по всей
 * поверхности — величина, от разбиения не зависящая, — а не поэлементные поля.
 */

#include "plod.h"
#include "poly_seg.h"
#include "polygon.h"
#include "psweep.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include "transport/dirs3.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#endif

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* ---------------------------------------------------------------- КАМЕРА
 *
 * РАСТЕРИЗАЦИЯ, А НЕ ЛУЧ, и это не вкус: §6 отвергла ray casting для переноса
 * доводом о трафике, а замер 08-02 подтвердил его на кадре — лучевая камера
 * стоила `693` с на городе против `6.4` с решения (А301). Здесь кадр есть ещё
 * один проход по треугольникам: столько же установок, сколько у ОДНОГО
 * направления развёртки.
 *
 * Перспективно-корректная интерполяция мирового положения через `1/w`; отсечение
 * по ближней плоскости честное (иначе камера ВНУТРИ улицы теряет стены, у
 * которых часть вершин позади неё). Значение пикселя — `L_out(u, v)` того
 * полигона, что выиграл глубину, то есть ровно решённое поле, без досчёта. */
typedef struct {
  double eye[3], fwd[3], rt[3], up[3];
  double fx, fy, cx, cy;
  int W, H;
} pf_cam;

static void pf_cam_make(pf_cam *c, const double eye[3], const double at[3], const double up0[3],
                        double fov_deg, int W, int H) {
  for (int a = 0; a < 3; a++)
    c->eye[a] = eye[a];
  double d[3], n = 0.0;
  for (int a = 0; a < 3; a++) {
    d[a] = at[a] - eye[a];
    n += d[a] * d[a];
  }
  n = sqrt(n);
  for (int a = 0; a < 3; a++)
    c->fwd[a] = d[a] / n;
  c->rt[0] = c->fwd[1] * up0[2] - c->fwd[2] * up0[1];
  c->rt[1] = c->fwd[2] * up0[0] - c->fwd[0] * up0[2];
  c->rt[2] = c->fwd[0] * up0[1] - c->fwd[1] * up0[0];
  n = sqrt(c->rt[0] * c->rt[0] + c->rt[1] * c->rt[1] + c->rt[2] * c->rt[2]);
  for (int a = 0; a < 3; a++)
    c->rt[a] /= n;
  c->up[0] = c->rt[1] * c->fwd[2] - c->rt[2] * c->fwd[1];
  c->up[1] = c->rt[2] * c->fwd[0] - c->rt[0] * c->fwd[2];
  c->up[2] = c->rt[0] * c->fwd[1] - c->rt[1] * c->fwd[0];
  c->W = W;
  c->H = H;
  c->fy = 0.5 * (double)H / tan(0.5 * fov_deg * 3.14159265358979323846 / 180.0);
  c->fx = c->fy;
  c->cx = 0.5 * (double)W;
  c->cy = 0.5 * (double)H;
}

/* Одна вершина в раму камеры: `(правая, верхняя, вперёд)`. */
static void pf_view(const pf_cam *c, const double p[3], double v[3]) {
  double d[3] = {0, 0, 0};
  for (int a = 0; a < 3; a++)
    d[a] = p[a] - c->eye[a];
  v[0] = d[0] * c->rt[0] + d[1] * c->rt[1] + d[2] * c->rt[2];
  v[1] = d[0] * c->up[0] + d[1] * c->up[1] + d[2] * c->up[2];
  v[2] = d[0] * c->fwd[0] + d[1] * c->fwd[1] + d[2] * c->fwd[2];
}

/* ПЕРЕМЕСТИТЬ ТЕЛО В САМОМ МЕШЕ (§124). Не маска и не пропуск при рисовании:
 * маска оставила бы всю переднюю часть конвейера нетронутой, и «быстро» вышло бы
 * по построению (возражение пользователя 08-02, и оно верное). Здесь меняется
 * ВХОД: треугольники, чей центр попал в шар, получают СОБСТВЕННЫЕ вершины и
 * сдвигаются, после чего сегментация и полигоны строятся заново с нуля.
 *
 * Вершины дублируются, а не двигаются на месте, потому что они общие: сдвиг
 * вершины утащил бы за собой соседние треугольники и порвал бы поверхность по
 * границе выделенного тела. Существующие треугольники сохраняют свои номера —
 * новые вершины ДОПИСЫВАЮТСЯ, — и потому неизменённые участки сегментируются
 * так же, что и даёт право переносить на них прежнее поле. */
static int pf_move_body(hz_objmesh *m, const double c[3], double r, const double dv[3],
                        int32_t *nmoved) {
  int32_t nsel = 0;
  for (int32_t t = 0; t < m->nt; t++) {
    double p[3][3], g[3] = {0, 0, 0};
    hz_obj_tri(m, t, p);
    for (int i = 0; i < 3; i++)
      for (int a = 0; a < 3; a++)
        g[a] += p[i][a] / 3.0;
    double d2 = 0.0;
    for (int a = 0; a < 3; a++)
      d2 += (g[a] - c[a]) * (g[a] - c[a]);
    if (d2 <= r * r) nsel++;
  }
  *nmoved = nsel;
  if (nsel == 0) return 0;
  double *nv = realloc(m->v, 3 * (size_t)(m->nv + 3 * nsel) * sizeof *nv);
  if (nv == NULL) return 2;
  m->v = nv;
  int32_t base = m->nv;
  for (int32_t t = 0; t < m->nt; t++) {
    double p[3][3], g[3] = {0, 0, 0};
    hz_obj_tri(m, t, p);
    for (int i = 0; i < 3; i++)
      for (int a = 0; a < 3; a++)
        g[a] += p[i][a] / 3.0;
    double d2 = 0.0;
    for (int a = 0; a < 3; a++)
      d2 += (g[a] - c[a]) * (g[a] - c[a]);
    if (d2 > r * r) continue;
    for (int i = 0; i < 3; i++) {
      for (int a = 0; a < 3; a++)
        m->v[3 * (size_t)base + (size_t)a] = p[i][a] + dv[a];
      m->f[3 * (size_t)t + (size_t)i] = base;
      base++;
    }
  }
  m->nv = base;
  for (int32_t k = 0; k < m->nv; k++)
    for (int a = 0; a < 3; a++) {
      double x = m->v[3 * (size_t)k + (size_t)a];
      if (x < m->lo[a]) m->lo[a] = x;
      if (x > m->hi[a]) m->hi[a] = x;
    }
  return 0;
}

/* ПЕРЕНОС ПОЛЯ НА НОВЫЙ НАБОР ПОЛИГОНОВ — ПО ПОЛОЖЕНИЮ, А НЕ ПО НОМЕРУ.
 * Номера после пересегментации не совпадают, а опорная точка неизменённого
 * участка совпадает ПОБИТОВО: те же треугольники в том же порядке дают тот же
 * центр площади. Поэтому ключ — координаты `org`, а несовпавшие полигоны
 * стартуют с нуля. Доля совпавших ПЕЧАТАЕТСЯ: без неё «тёплый старт» неотличим
 * от холодного. */
static uint64_t pf_key(const double p[3]) {
  uint64_t h = 0xcbf29ce484222325ULL;
  for (int a = 0; a < 3; a++) {
    double x = p[a];
    if (!(x < 0.0) && !(x > 0.0)) x = 0.0;
    uint64_t b;
    memcpy(&b, &x, sizeof b);
    h ^= b;
    h *= 0x100000001b3ULL;
  }
  return h;
}

/* ПЕЧЬ — ЕДИНСТВЕННЫЙ ЭТАЛОН, НЕ ЗАВИСЯЩИЙ НИ ОТ ОДНОЙ ИЗ ДВУХ СХЕМ (§126).
 *
 * Замкнутая коробка, у всех граней одно альбедо `ρ` и одно собственное излучение
 * `Le`. Точное решение известно замкнутой формой и не требует ни геометрии, ни
 * видимости: в замкнутой полости `Σf = 1` у каждой площадки, поэтому
 *
 *     B = Le + ρ·B   ⇒   B = Le/(1 − ρ)
 *
 * во всех точках. При `Le = 1`, `ρ = 0.5` это ровно `2`. Расхождение схемы с этим
 * числом есть ЕЁ ошибка, а не спор двух подозреваемых: сегодня хранимый оператор
 * и итерации разошлись на городе в `1.93` раза, и без внешнего эталона нельзя
 * сказать, кто из них неправ.
 *
 * Грани смотрят ВНУТРЬ (нормаль к центру), иначе полость не замкнута для
 * переноса. `nsub` — на сколько делить сторону: одна грань из одного полигона
 * проверяет физику, из многих — ещё и разбиение. */
static int pf_oven(hz_objmesh *m, double side, int nsub) {
  memset(m, 0, sizeof *m);
  for (int a = 0; a < 3; a++) {
    m->lo[a] = 1e300;
    m->hi[a] = -1e300;
  }
  const double q = side / (double)nsub, hq = 0.5 * q, s2 = 0.5 * side;
  for (int ax = 0; ax < 3; ax++)
    for (int sgn = -1; sgn <= 1; sgn += 2) {
      double n[3] = {0, 0, 0}, eu[3] = {0, 0, 0};
      n[ax] = -(double)sgn; /* внутрь коробки */
      eu[(ax + 1) % 3] = 1.0;
      for (int i = 0; i < nsub; i++)
        for (int j = 0; j < nsub; j++) {
          double c[3];
          c[ax] = (double)sgn * s2;
          c[(ax + 1) % 3] = -s2 + ((double)i + 0.5) * q;
          c[(ax + 2) % 3] = -s2 + ((double)j + 0.5) * q;
          if (hz_obj_add_quad(m, "oven", c, n, eu, hq, hq) < 0) return 1;
        }
    }
  return 0;
}

int main(int argc, char **argv) {
  int city = 0, nmu = 2, nphi = 4, nb = 21, layout = HZ_LAYOUT_RUNS, imgw = 960, imgh = 540;
  int fold = 0; /* §125: угловая квадратура едет на итерации, а не вложена в неё */
  int oven = 0; /* §126: печь — аналитический эталон B = Le/(1−ρ) */
  int osub = 4; /* делений стороны печи */
  /* НЕГАТИВНЫЙ КОНТРОЛЬ §127: поправку квадратуры выключить. Ошибка печи обязана
   * вернуться к прежнему проценту; не вернулась — поправка ни при чём. */
  int noqn = 0;
  int skyf = 0; /* §130: замерить долю неба и сверить с хранимым оператором */
  int lod = 0;  /* §131: перенос по срезу лестницы вместо всех участков */
  const char *fsave = NULL, *fcmp = NULL; /* §132: поле на диск и сверка */
  int shells = 0;                         /* §133: замер цены каскада по оболочкам */
  double szmul = 0.0;                     /* §136: второй предел среза, в допусках ε */
  int link1 = 0;                          /* §137 З1: связность направлений по парам */
  /* §128: направлений за шаг свёртки. `0` читается как число потоков OpenMP —
   * иначе параллелизм по направлениям простаивает. */
  int nfold = 0;
  double oside = 4.0;
  double h = 0.5, rho = 0.5, sky = HZ_CFG_SKY_LE;
  /* Изменение сцены: шар радиуса `mvr` вокруг `mvc` сдвигается на `mvd`. */
  double mvr = 0.0, mvc[3] = {0.0, 0.0, 0.0}, mvd[3] = {0.0, 0.0, 0.0};
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "city") == 0) city = 1;
    if (strncmp(argv[i], "nmu=", 4) == 0) nmu = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "nphi=", 5) == 0) nphi = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "nb=", 3) == 0) nb = (int)strtol(argv[i] + 3, NULL, 10);
    if (strncmp(argv[i], "h=", 2) == 0) h = strtod(argv[i] + 2, NULL);
    if (strncmp(argv[i], "rho=", 4) == 0) rho = strtod(argv[i] + 4, NULL);
    if (strncmp(argv[i], "sky=", 4) == 0) sky = strtod(argv[i] + 4, NULL);
    if (strcmp(argv[i], "list") == 0) layout = HZ_LAYOUT_LIST;
    if (strcmp(argv[i], "fold") == 0) fold = 1;
    if (strcmp(argv[i], "oven") == 0) oven = 1;
    if (strcmp(argv[i], "noqn") == 0) noqn = 1;
    if (strcmp(argv[i], "skyf") == 0) skyf = 1;
    if (strcmp(argv[i], "lod") == 0) lod = 1;
    if (strncmp(argv[i], "save=", 5) == 0) fsave = argv[i] + 5;
    if (strncmp(argv[i], "cmp=", 4) == 0) fcmp = argv[i] + 4;
    if (strcmp(argv[i], "shells") == 0) shells = 1;
    if (strncmp(argv[i], "sz=", 3) == 0) szmul = strtod(argv[i] + 3, NULL);
    if (strcmp(argv[i], "link1") == 0) link1 = 1;
    if (strncmp(argv[i], "nfold=", 6) == 0) nfold = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "osub=", 5) == 0) osub = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "w=", 2) == 0) imgw = (int)strtol(argv[i] + 2, NULL, 10);
    if (strncmp(argv[i], "ih=", 3) == 0) imgh = (int)strtol(argv[i] + 3, NULL, 10);
    if (strncmp(argv[i], "mvr=", 4) == 0) mvr = strtod(argv[i] + 4, NULL);
    if (strncmp(argv[i], "mvc=", 4) == 0)
      sscanf(argv[i] + 4, "%lf,%lf,%lf", &mvc[0], &mvc[1], &mvc[2]);
    if (strncmp(argv[i], "mvd=", 4) == 0)
      sscanf(argv[i] + 4, "%lf,%lf,%lf", &mvd[0], &mvd[1], &mvd[2]);
  }

  double t0 = now_s();
  hz_objmesh m;
  if (oven) {
    if (pf_oven(&m, oside, osub) != 0) return 1;
    sky = 0.0; /* печь замкнута: неба в ней нет */
  } else if (hz_obj_load(&m, city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ,
                         city ? HZ_CFG_CITY_SCALE : HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  double t_load = now_s() - t0;

  /* Предел габарита участка снят у города по А285: `0.5` м мельче его
   * единичного треугольника, и сегментация выродилась бы в участок на
   * треугольник. У зала он остаётся прежним. */
  const double dseg = city ? 0.05 : 0.045;
  t0 = now_s();
  hz_pseglist sg;
  /* ПРЕДЕЛ ГАБАРИТА — ПО СЦЕНЕ, А НЕ ПО УМОЛЧАНИЮ (А285). У печи грань крупная,
   * и залские `0.5` м раздробили бы её в участок на треугольник — сторож ниже
   * это и поймал при первом же прогоне. */
  const double scap = (city || oven) ? 0.0 : 0.5;
  if (hz_seg_planar_cap(&sg, &m, dseg, scap) != 0) return 1;
  if (sg.nseg >= m.nt) {
    fprintf(stderr, "сегментация выродилась: участков %d при %d треугольниках\n", sg.nseg, m.nt);
    return 1;
  }
  hz_polyset ps;
  if (hz_poly_build(&ps, &m, &sg) != 0) return 1;
  double t_geo = now_s() - t0;
  printf("== СЦЕНА: %s, треугольников %d, полигонов %d; разбор %.1f с, геометрия %.1f с\n",
         city ? "ГОРОД" : "зал", m.nt, ps.np, t_load, t_geo);

  /* ---------------------------------------------- LOD В ПЕРЕНОСЕ (§131) ---
   *
   * Перенос до сих пор шёл по ВСЕМ участкам нулевого уровня — `446 013` у
   * города, — хотя лестница построена и срез камеры даёт `152 618`. Это `2.9×`
   * по числу примитивов, лежавшие нетронутыми.
   *
   * Машинерия готова целиком: `hz_lod_seglist` превращает срез в обычную
   * разметку треугольников, а `hz_poly_build` строит по ней полигоны. Ниже по
   * течению ничто не меняется — развёртка не знает, с какого уровня пришли
   * элементы.
   *
   * ЦЕНА НАЗВАНА ЗАРАНЕЕ: лестница строится `≈100` с, и это ПРЕДРАСЧЁТ, который
   * при изменении геометрии придётся повторять. Здесь мерится только выигрыш на
   * переносе; окупаемость — отдельный счёт. */
  hz_lod L;
  hz_pseglist sgc;
  hz_polyset psc;
  int lod_on = 0;
  if (lod) {
    double tl = now_s();
    hz_lodcfg lc;
    memset(&lc, 0, sizeof lc);
    lc.delta0 = dseg;
    const double eps = (HZ_CFG_FOV_DEG * M_PI / 180.0) / 1920.0;
    lc.eps = eps;
    lc.maxlev = 9;
    lc.radmul = 4.0;
    lc.base = 1.4142;
    lc.szmul = szmul;
    if (hz_lod_build_merge(&L, &m, &sg, &ps, &lc) != 0) {
      fprintf(stderr, "отказ лестницы\n");
      return 1;
    }
    double t_lad = now_s() - tl;
    double eyeh[3] = HZ_CFG_HALL_EYE, eyec[3] = HZ_CFG_CITY_EYE;
    const double *ey = city ? eyec : eyeh;
    double eye[3] = {ey[0], ey[1], ey[2]};
    int32_t *cut = malloc((size_t)L.np * sizeof *cut);
    if (cut == NULL) return 1;
    tl = now_s();
    int32_t ncut = hz_lod_cut(&L, eye, eps, 1, cut);
    if (hz_lod_seglist(&L, &m, &sg, cut, &sgc) != 0) return 1;
    if (hz_poly_build(&psc, &m, &sgc) != 0) return 1;
    free(cut);
    printf("== LOD В ПЕРЕНОСЕ: лестница %d уровней, %d узлов за %.1f с; срез %d узлов, "
           "полигонов среза %d (было %d, то есть %.2f×) за %.1f с\n",
           L.nlev, L.nnd, t_lad, ncut, psc.np, ps.np, (double)ps.np / (double)psc.np, now_s() - tl);
    /* ПОЛИГОНЫ СРЕЗА ЗАМЕЩАЮТ ИСХОДНЫЕ ЦЕЛИКОМ. Ниже по течению ничего не
     * меняется: развёртка не знает, с какого уровня пришёл элемент. */
    hz_poly_free(&ps);
    ps = psc;
    lod_on = 1;
  }
  if (lod_on && mvr > 0.0) {
    fprintf(stderr, "lod и mvr вместе не поддержаны: после правки меша лестницу надо строить "
                    "заново, и это отдельный замер\n");
    return 1;
  }

  hz_ptrans tr;
  if (hz_ptrans_init(&tr, &ps, &m) != 0) return 1;
  /* Альбедо ставится ОДИНАКОВЫМ с эталоном (§120: `ρ = 0.5` у всех), иначе
   * сравнивать потоки нельзя. Собственного излучения нет: источник — небо. */
  for (int32_t k = 0; k < ps.np; k++) {
    tr.rho[k] = rho;
    /* В печи светятся ВСЕ грани: тогда точное решение есть `Le/(1−ρ)` без
     * всякой геометрии, и сверять можно прямо с ним. */
    tr.Le[k] = oven ? 1.0 : 0.0;
  }

  tr3_dirs d;
  if (tr3_dirs_product(&d, nmu, nphi) != 0) return 1;
  printf("== НАПРАВЛЕНИЙ: %d (nmu %d × nphi %d); шаг растра %.2f м; отскоков %d; ρ %.2f; "
         "небо %.2f\n",
         d.n, nmu, nphi, h, nb, rho, sky);

  /* СВЁРНУТАЯ СХЕМА (§125): угловая квадратура НЕ вложена в итерацию, а едет на
   * ней. На каждой итерации берётся ОДНО направление, следующее по списку;
   * накопленная сумма делится на накопленный вес, то есть оценка есть
   * средневзвешенное по УЖЕ ПОСЕЩЁННЫМ направлениям. После полного цикла по
   * списку она совпадает с обычной квадратурой ТОЧНО, а до того даёт грубое, но
   * несмещённое приближение — и первые итерации, где поле всё равно сырое, не
   * платят полную угловую цену.
   *
   * Порядок обхода — с шагом, взаимно простым с длиной списка: соседние итерации
   * берут далёкие друг от друга направления, иначе первые оценки все смотрели бы
   * в одну сторону. Случайности здесь нет и не нужно: прогон обязан повторяться. */
  double *sacc = NULL, *fox = NULL, *foy = NULL, *foz = NULL, *fw = NULL;
  double wvis = 0.0, wtot = 0.0;
  tr3_dirs d1;
  int step1 = 1;
  if (fold) {
    sacc = calloc(3 * (size_t)ps.np, sizeof *sacc);
    if (sacc == NULL) return 1;
    memset(&d1, 0, sizeof d1);
    if (nfold < 1) {
      nfold = 1;
#ifdef _OPENMP
      nfold = omp_get_max_threads();
#endif
    }
    if (nfold > d.n) nfold = d.n;
    d1.n = nfold;
    fox = malloc((size_t)nfold * sizeof *fox);
    foy = malloc((size_t)nfold * sizeof *foy);
    foz = malloc((size_t)nfold * sizeof *foz);
    fw = malloc((size_t)nfold * sizeof *fw);
    if (fox == NULL || foy == NULL || foz == NULL || fw == NULL) return 1;
    for (int k = 0; k < d.n; k++)
      wtot += d.w[k];
    /* Шаг: наибольшее целое ниже d.n/φ, взаимно простое с d.n (аддитивная
     * рекурсия Вейля — стандартный способ обойти список «вразброс»). */
    step1 = (int)((double)d.n * 0.6180339887498949);
    if (step1 < 1) step1 = 1;
    while (step1 > 1) {
      int a = step1, b2 = d.n;
      while (b2 != 0) {
        int t2 = a % b2;
        a = b2;
        b2 = t2;
      }
      if (a == 1) break;
      step1--;
    }
    printf("== СВЁРНУТО: %d направлений на итерацию, шаг обхода %d из %d; Σ весов %.4f\n", nfold,
           step1, d.n, wtot);
  }

  /* ПОПРАВКА КВАДРАТУРЫ ПО НОРМАЛИ (§127) — один раз на полигон, до переноса.
   * `S(n) = Σ_{ω·n<0} w_ω|ω·n|` обязано равняться `π`; отношение и есть дефект
   * набора направлений для этой ориентации. Стоит один проход по полигонам, то
   * есть нисколько против растеризации. */
  double *qn = NULL;
  if (!noqn) {
    qn = malloc((size_t)ps.np * sizeof *qn);
    if (qn == NULL) return 1;
    double qlo = 1e300, qhi = -1e300;
    for (int32_t k = 0; k < ps.np; k++) {
      const double *nn = ps.p[k].n;
      double s = 0.0;
      for (int j = 0; j < d.n; j++) {
        double c = d.ox[j] * nn[0] + d.oy[j] * nn[1] + d.oz[j] * nn[2];
        if (c < 0.0) s += d.w[j] * (-c);
      }
      qn[k] = s / 3.14159265358979323846;
      if (qn[k] < qlo) qlo = qn[k];
      if (qn[k] > qhi) qhi = qn[k];
    }
    tr.qn = qn;
    printf("== ПОПРАВКА КВАДРАТУРЫ (§127): S(n)/π по полигонам от %.4f до %.4f "
           "(обязано быть 1; отклонение и есть угловая ошибка набора)\n",
           qlo, qhi);
  }

  /* З1 (§137, правлено §138): СВЯЗНОСТЬ НАПРАВЛЕНИЙ ПО ПАРАМ «приёмник ←
   * излучатель», а НЕ по пикселям — рама растра строится от направления, и «тот
   * же пиксель» у разных `ω` не определён (А307).
   *
   * Печатается КРИВАЯ по угловому расстоянию, а не одно число: порога `0.5` в
   * приёмке нет (А308), решение принимается по кривой. Соседство берётся по
   * ОБЕИМ осям набора порознь и выражается углом в градусах (А309).
   * Негативный контроль — противоположное направление: совпадение обязано
   * рухнуть, иначе хеш считает не то. */
  if (link1) {
    const int64_t PCAP = 40000000; /* пар на направление; при переполнении — доклад */
    tr.pair_k = malloc((size_t)PCAP * sizeof *tr.pair_k);
    tr.pair_j = malloc((size_t)PCAP * sizeof *tr.pair_j);
    tr.pair_w = malloc((size_t)PCAP * sizeof *tr.pair_w);
    if (tr.pair_k == NULL || tr.pair_j == NULL || tr.pair_w == NULL) {
      free(tr.pair_k);
      free(tr.pair_j);
      free(tr.pair_w);
      return 1;
    }
    /* ПОЛЕ СПЕРВА ЗАЖИГАЕТСЯ, иначе пар нет вовсе: на первой итерации `B ≡ 0`,
     * весь свет идёт от неба, а небесные фрагменты в пары не входят (у них нет
     * излучателя-поверхности). Три отскока полным набором дают поле, на котором
     * связность уже осмысленна. */
    for (int b = 0; b < 3; b++) {
      hz_pstats s0;
      memset(&s0, 0, sizeof s0);
      if (hz_psweep_bounce(&tr, &d, h, 0, sky, layout, &s0) != 0) return 1;
    }
    /* Опорное направление — первое; сравниваемые: соседнее по `μ`, соседнее по
     * `φ`, через одно по `φ`, и противоположное (контроль). */
    int base_i = 0;
    int cand[4] = {1, nphi * 2, nphi * 4, d.n - 1};
    int64_t nsl = 4;
    while (nsl < 4 * PCAP / 8)
      nsl *= 2;
    uint64_t *key = calloc((size_t)nsl, sizeof *key);
    double *val = calloc((size_t)nsl, sizeof *val);
    if (key == NULL || val == NULL) {
      free(key);
      free(val);
      free(tr.pair_k);
      free(tr.pair_j);
      free(tr.pair_w);
      return 1;
    }
    uint64_t msk = (uint64_t)nsl - 1;
    double wbase = 0.0;
    tr3_dirs dz;
    memset(&dz, 0, sizeof dz);
    dz.n = 1;
    for (int pass = 0; pass < 5; pass++) {
      int k = (pass == 0) ? base_i : cand[pass - 1];
      /* Сторож только для СРАВНИВАЕМЫХ: опорное имеет номер 0 и отбрасываться
       * не должно — первая редакция на этом и дала ноль пар. */
      if (k < 0 || k >= d.n || (pass > 0 && k == base_i)) continue;
      dz.ox = d.ox + k;
      dz.oy = d.oy + k;
      dz.oz = d.oz + k;
      dz.w = d.w + k;
      tr.pair_n = 0;
      tr.pair_cap = PCAP;
      hz_psweep_zero(&tr);
      hz_pstats s1;
      memset(&s1, 0, sizeof s1);
      if (hz_psweep_gather(&tr, &dz, h, 1, sky, layout, &s1) != 0) return 1;
      if (pass == 0) {
        for (int64_t q = 0; q < tr.pair_n; q++) {
          uint64_t kk = ((uint64_t)(uint32_t)tr.pair_k[q] << 32) | (uint32_t)tr.pair_j[q];
          uint64_t i = (kk * 0x9E3779B97F4A7C15ULL) & msk;
          while (key[i] != 0 && key[i] != kk + 1)
            i = (i + 1) & msk;
          key[i] = kk + 1;
          val[i] += tr.pair_w[q];
          wbase += tr.pair_w[q];
        }
        printf("== З1: опорное направление, пар %lld, суммарный вес %.6e%s\n", (long long)tr.pair_n,
               wbase, (tr.pair_n >= PCAP) ? "  [ЁМКОСТЬ ИСЧЕРПАНА]" : "");
        printf("   к чему сравнивается   угол, °   пар   совпало по счёту   совпало ПО ВЕСУ\n");
        continue;
      }
      double dot = d.ox[base_i] * d.ox[k] + d.oy[base_i] * d.oy[k] + d.oz[base_i] * d.oz[k];
      if (dot > 1.0) dot = 1.0;
      if (dot < -1.0) dot = -1.0;
      double ang = acos(dot) * 180.0 / M_PI;
      int64_t nhit = 0;
      double whit = 0.0, wall = 0.0;
      for (int64_t q = 0; q < tr.pair_n; q++) {
        uint64_t kk = ((uint64_t)(uint32_t)tr.pair_k[q] << 32) | (uint32_t)tr.pair_j[q];
        uint64_t i = (kk * 0x9E3779B97F4A7C15ULL) & msk;
        while (key[i] != 0 && key[i] != kk + 1)
          i = (i + 1) & msk;
        wall += tr.pair_w[q];
        if (key[i] == kk + 1) {
          nhit++;
          whit += tr.pair_w[q];
        }
      }
      const char *nm = (pass == 1)   ? "соседнее по mu"
                       : (pass == 2) ? "соседнее по phi"
                       : (pass == 3) ? "через одно по phi"
                                     : "ПРОТИВОПОЛОЖНОЕ [НК]";
      printf("   %-21s %7.2f  %9lld   %14.4f   %14.4f\n", nm, ang, (long long)tr.pair_n,
             (double)nhit / (double)(tr.pair_n ? tr.pair_n : 1), whit / (wall > 0.0 ? wall : 1.0));
      fflush(stdout);
    }
    free(key);
    free(val);
    free(tr.pair_k);
    free(tr.pair_j);
    free(tr.pair_w);
    tr.pair_k = NULL;
    tr.pair_j = NULL;
    tr.pair_w = NULL;
    tr.pair_cap = 0;
  }

  /* ЦЕНА КАСКАДА ПО ОБОЛОЧКАМ — АРИФМЕТИКА ПО ГОТОВЫМ ДАННЫМ (§133, пункт 1).
   *
   * Ни одной строки в библиотеке: полигоны уже есть, срез уже есть, направления
   * уже есть. Считается ровно то, что решает, стоит ли писать каскад — сумма
   * пикселей по оболочкам против нынешних `W×H` на всю сцену.
   *
   * Экран оболочки накрывает её приёмников, но РИСУЮТСЯ в него все полигоны,
   * чья проекция туда попадает, — иначе теряется дальний заслонитель. Поэтому
   * считаются два числа: пиксели (цена заливки) и число проецирующихся
   * полигонов (цена установки, и она же избыточность). */
  if (shells) {
    const int KMAX = 12;
    double eyeh[3] = HZ_CFG_HALL_EYE, eyec[3] = HZ_CFG_CITY_EYE;
    const double *ey = city ? eyec : eyeh;
    /* Оболочка полигона — по расстоянию от глаза, границы удваиваются от `R₀`.
     * `R₀` не выбирается: это шаг растра, поделённый на угловой допуск пикселя,
     * то есть дальность, на которой элемент размера `h` виден под тем же углом,
     * что и пиксель кадра. Порога здесь нет — есть уже принятые `h` и `ε`. */
    /* КРИТЕРИЙ ГЕОМЕТРИЧЕСКИЙ, А НЕ УГЛОВОЙ, и первый прогон это доказал (А305):
     * `R₀ = h/ε` дал `916.7` м, то есть весь город в одной оболочке. Смысл
     * обратный ожидаемому — при `ε = 5.45e-4` рад пиксель кадра на `10` м
     * покрывает `5.5` мм, и `h = 0.5` м уже грубее пикселя ВБЛИЗИ. Значит шаг
     * растра ведётся не дальностью, а РАЗМЕРОМ ЭЛЕМЕНТА: §130 требует, чтобы
     * `h` разрешал элемент, иначе вес `h²` переучитывает. */
    const double R0 = 0.0;
    int32_t *sh = malloc((size_t)ps.np * sizeof *sh);
    if (sh == NULL) return 1;
    int nsh = 0;
    for (int32_t k = 0; k < ps.np; k++) {
      /* Размер элемента — корень из площади: ровно та величина, которую обязан
       * разрешать шаг растра. Оболочка `s` берёт `h_s = h·2^s`, поэтому элемент
       * попадает в ту, где `h_s` впервые не крупнее его самого. */
      double dsz = sqrt(ps.p[k].area);
      int s = (dsz > h) ? (int)floor(log(dsz / h) / log(2.0)) : 0;
      if (s < 0) s = 0;
      if (s >= KMAX) s = KMAX - 1;
      sh[k] = s;
      if (s + 1 > nsh) nsh = s + 1;
    }
    (void)ey;
    (void)R0;
    printf("== ЦЕНА КАСКАДА (§133): оболочка по РАЗМЕРУ элемента, h₀ = %.2f м, оболочек %d\n", h,
           nsh);
    printf(
        "  об  граница, м   h_k, м   полигонов   габарит, м        пикселей   рисуется полиг.\n");
    double sumpix = 0.0, sumdraw = 0.0;
    for (int s = 0; s < nsh; s++) {
      double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
      int32_t ncnt = 0;
      for (int32_t k = 0; k < ps.np; k++) {
        if (sh[k] != s) continue;
        ncnt++;
        for (int a = 0; a < 3; a++) {
          if (ps.p[k].org[a] < lo[a]) lo[a] = ps.p[k].org[a];
          if (ps.p[k].org[a] > hi[a]) hi[a] = ps.p[k].org[a];
        }
      }
      if (ncnt == 0) continue;
      double hk = h * pow(2.0, (double)s);
      /* Пиксели — по НАИБОЛЬШЕЙ проекции габарита оболочки: экран каждого
       * направления не крупнее её. Оценка сверху, и она честная. */
      double e0 = hi[0] - lo[0], e1 = hi[1] - lo[1], e2 = hi[2] - lo[2];
      double q1 = (e0 > e1) ? e0 : e1,
             q2 = (e2 > ((e0 < e1) ? e0 : e1)) ? e2 : ((e0 < e1) ? e0 : e1);
      double px = (q1 / hk + 2.0) * (q2 / hk + 2.0);
      /* Сколько полигонов рисуется в этот экран: те, чья опорная точка попала в
       * габарит оболочки, растянутый вдоль всех осей (заслонители приходят
       * отовсюду вдоль `ω`, поэтому берётся полный габарит сцены по одной оси). */
      int32_t ndraw = 0;
      for (int32_t k = 0; k < ps.np; k++) {
        const double *o = ps.p[k].org;
        int in = 1;
        for (int a = 0; a < 3 && in; a++)
          if (o[a] < lo[a] - hk || o[a] > hi[a] + hk) in = 0;
        if (in) ndraw++;
      }
      sumpix += px;
      sumdraw += ndraw;
      printf("  %2d  %10.2f  %7.2f  %10d   %5.0f×%5.0f×%5.0f  %12.0f   %10d\n", s, hk, hk, ncnt, e0,
             e1, e2, px, ndraw);
    }
    double now_px = 0.0;
    {
      double e0 = m.hi[0] - m.lo[0], e1 = m.hi[1] - m.lo[1], e2 = m.hi[2] - m.lo[2];
      double q1 = (e0 > e1) ? e0 : e1,
             q2 = (e2 > ((e0 < e1) ? e0 : e1)) ? e2 : ((e0 < e1) ? e0 : e1);
      now_px = (q1 / h + 2.0) * (q2 / h + 2.0);
    }
    printf("== ИТОГ: пикселей по оболочкам %.0f против %.0f сейчас — **%.2f×**; "
           "полигонов рисуется суммарно %.0f против %d, избыточность %.2f\n",
           sumpix, now_px, now_px / (sumpix > 0.0 ? sumpix : 1.0), sumdraw, ps.np,
           sumdraw / (double)ps.np);
    free(sh);
  }

  /* ДОЛЯ ПОЛУСФЕРЫ, УШЕДШАЯ В НЕБО (§130) — БЕЗ ЕДИНОЙ СТРОКИ В БИБЛИОТЕКЕ.
   * При `ρ = 0` и `Le = 0` поле тождественно ноль, поэтому единственное, что
   * попадает в накопитель за один отскок, — небо. Тогда `acc₀[k]` есть в
   * точности `L_неба · ∫(по небу) cosθ dω` по площади полигона, а полная
   * полусфера даёт `π·A`. Отношение и есть искомая доля, сравнимая с `f_sky`
   * хранимого оператора один в один (там это `Σ dff` по пустым пикселям, а
   * `Σ dff` по всей полусфере равна единице).
   *
   * Величина эта — главный подозреваемый в расхождении `1.93` раза: у оператора
   * она снята полукубом в опорной точке ГРУБОГО узла среза, здесь — по
   * фрагментам всех `446 013` полигонов. */
  if (skyf) {
    for (int32_t k = 0; k < ps.np; k++) {
      tr.rho[k] = 0.0;
      tr.Le[k] = 0.0;
    }
    hz_pstats s0;
    memset(&s0, 0, sizeof s0);
    if (hz_psweep_bounce(&tr, &d, h, 0, sky, layout, &s0) != 0) return 1;
    double sa = 0.0, ss = 0.0, wmin = 1e300, wmax = -1e300;
    int64_t nzero = 0;
    for (int32_t k = 0; k < ps.np; k++) {
      double A = ps.p[k].area;
      if (!(A > 0.0)) continue;
      double f = tr.acc[3 * (size_t)k] / (M_PI * A * (sky > 0.0 ? sky : 1.0));
      sa += A;
      ss += f * A;
      if (f < wmin) wmin = f;
      if (f > wmax) wmax = f;
      if (!(f > 1e-6)) nzero++;
    }
    printf("== ДОЛЯ НЕБА (§130), взвешенная по площади: %.4f; по полигонам от %.4f до %.4f; "
           "не видят неба вовсе %lld из %d (%.1f %%)\n",
           (sa > 0.0) ? ss / sa : 0.0, wmin, wmax, (long long)nzero, ps.np,
           100.0 * (double)nzero / (double)(ps.np ? ps.np : 1));
    printf("   для сверки — хранимый оператор (§120): элементов с видимым небом 94 181 из "
           "152 618, средняя доля 0.1967 ПО НИМ, то есть по всему срезу 0.1214\n");
    /* Альбедо и свечение возвращаются: дальше идёт обычный прогон. */
    for (int32_t k = 0; k < ps.np; k++) {
      tr.rho[k] = rho;
      tr.Le[k] = oven ? 1.0 : 0.0;
    }
    hz_psweep_zero(&tr);
    memset(tr.E, 0, 3 * (size_t)ps.np * sizeof *tr.E);
  }

  printf("  отскок   Σ B·A, Вт/ср      dE      фрагментов    растр,с   редукция,с   всего,с\n");
  double ttot = 0.0;
  for (int b = 0; b < nb; b++) {
    hz_pstats st;
    memset(&st, 0, sizeof st);
    double tb = now_s();
    if (!fold) {
      if (hz_psweep_bounce(&tr, &d, h, 0, sky, layout, &st) != 0) {
        fprintf(stderr, "отказ итерации %d\n", b);
        return 1;
      }
    } else {
      /* НАПРАВЛЕНИЙ ЗА ШАГ — ПО ЧИСЛУ ПОТОКОВ, А НЕ ОДНО (§128, замерено).
       * Параллелизм в `psweep` устроен ПО НАПРАВЛЕНИЯМ (`psweep.c:323`), поэтому
       * одно направление на итерацию оставляет пятнадцать потоков из шестнадцати
       * стоять на барьере: профиль дал `40.4 %` в `gomp_barrier_wait_end`.
       * Свёртка при этом не теряется — за шаг берётся `nfold` направлений из
       * той же последовательности Вейля, а не весь набор. */
      for (int q = 0; q < nfold; q++) {
        int k = (int)((((int64_t)b * nfold + q) * step1) % d.n);
        fox[q] = d.ox[k];
        foy[q] = d.oy[k];
        foz[q] = d.oz[k];
        fw[q] = d.w[k];
        wvis += d.w[k];
      }
      d1.ox = fox;
      d1.oy = foy;
      d1.oz = foz;
      d1.w = fw;
      hz_psweep_zero(&tr);
      if (hz_psweep_gather(&tr, &d1, h, 0, sky, layout, &st) != 0) {
        fprintf(stderr, "отказ итерации %d\n", b);
        return 1;
      }
      double sc = wtot / wvis;
      for (int32_t q = 0; q < 3 * ps.np; q++) {
        sacc[q] += tr.acc[q];
        tr.acc[q] = sacc[q] * sc;
      }
      hz_psweep_solve(&tr, &st);
    }
    double dt = now_s() - tb;
    ttot += dt;
    /* ПОТОК ПО ВСЕЙ ПОВЕРХНОСТИ — величина, не зависящая от разбиения, и потому
     * единственная законная для сверки с эталоном, решённым на другом наборе
     * элементов. `L_out = Le + ρ·E/π`, поток на элемент — `π·L·A`. */
    double flux = 0.0;
    for (int32_t k = 0; k < ps.np; k++)
      flux += hz_ptrans_lout(&tr, k, 0.0, 0.0) * ps.p[k].area;
    printf("  %4d   %14.6e  %8.2e   %11lld   %8.2f   %8.2f   %8.2f\n", b + 1, flux, st.dE,
           (long long)st.nfrag, st.t_raster, st.t_reduce, dt);
    fflush(stdout);
  }
  if (oven) {
    /* СВЕРКА С ЗАМКНУТОЙ ФОРМОЙ: `B = Le/(1−ρ)` в каждой точке. Хвост, а не
     * среднее (§4): среднее скрыло бы, что часть граней недобирает. */
    double ex = 1.0 / (1.0 - rho), wm = 0.0, sw2 = 0.0, sa = 0.0;
    for (int32_t k = 0; k < ps.np; k++) {
      double b = hz_ptrans_lout(&tr, k, 0.0, 0.0);
      double e = fabs(b - ex) / ex;
      if (e > wm) wm = e;
      sw2 += e * ps.p[k].area;
      sa += ps.p[k].area;
    }
    printf("== ПЕЧЬ: точное B = Le/(1−ρ) = %.6f; ошибка по площади %.3e, ХУДШАЯ %.3e; "
           "полигонов %d\n",
           ex, (sa > 0.0) ? sw2 / sa : 0.0, wm, ps.np);
  }
  printf("== ВСЕГО НА ПЕРЕНОС: %.1f с на %d отскоков (%.2f с на отскок, %.3f с на направление)\n",
         ttot, nb, ttot / (double)(nb > 0 ? nb : 1),
         ttot / (double)((nb > 0 ? nb : 1) * (d.n > 0 ? d.n : 1)));
  printf("== ЭТАЛОН ХРАНИМЫМ ОПЕРАТОРОМ (§120, тот же город, ρ = 0.5, небо 1.0): "
         "поток 2.398505e+05 Вт/ср, сборка 3102.6 с\n");

  /* ------------------------------- ИЗМЕНЕНИЕ СЦЕНЫ И ТЁПЛЫЙ СТАРТ (§124) --
   *
   * Вопрос: «можно ли относительно немного считать, когда сцена не сильно
   * меняется». Проверяется ЧЕСТНО — правкой самого меша и полным повтором
   * переднего края, а не пропуском полигонов при рисовании. */
  if (mvr > 0.0) {
    int32_t np0 = ps.np;
    double *org0 = malloc(3 * (size_t)np0 * sizeof *org0);
    double *E0 = malloc(3 * (size_t)np0 * sizeof *E0);
    if (org0 == NULL || E0 == NULL) return 1;
    for (int32_t k = 0; k < np0; k++)
      for (int a = 0; a < 3; a++) {
        org0[3 * (size_t)k + (size_t)a] = ps.p[k].org[a];
        E0[3 * (size_t)k + (size_t)a] = tr.E[3 * (size_t)k + (size_t)a];
      }

    int32_t nmoved = 0;
    double t_mv = now_s();
    if (pf_move_body(&m, mvc, mvr, mvd, &nmoved) != 0) return 1;
    t_mv = now_s() - t_mv;

    double t_re = now_s();
    hz_ptrans_free(&tr);
    hz_poly_free(&ps);
    hz_seg_free(&sg);
    if (hz_seg_planar_cap(&sg, &m, dseg, scap) != 0) return 1;
    if (hz_poly_build(&ps, &m, &sg) != 0) return 1;
    if (hz_ptrans_init(&tr, &ps, &m) != 0) return 1;
    for (int32_t k = 0; k < ps.np; k++) {
      tr.rho[k] = rho;
      tr.Le[k] = 0.0;
    }
    t_re = now_s() - t_re;
    printf("== СЦЕНА ИЗМЕНЕНА: тело радиуса %.1f м сдвинуто на (%.1f, %.1f, %.1f), "
           "треугольников тронуто %d (%.3f %%); правка меша %.2f с, передний край заново %.1f с; "
           "полигонов %d -> %d\n",
           mvr, mvd[0], mvd[1], mvd[2], nmoved, 100.0 * (double)nmoved / (double)m.nt, t_mv, t_re,
           np0, ps.np);

    /* Перенос поля по ПОЛОЖЕНИЮ опорной точки. */
    int64_t nsl = 4;
    while (nsl < 2 * (int64_t)np0 + 8)
      nsl *= 2;
    uint64_t *key = calloc((size_t)nsl, sizeof *key);
    int32_t *val = malloc((size_t)nsl * sizeof *val);
    if (key == NULL || val == NULL) return 1;
    for (int64_t i = 0; i < nsl; i++)
      val[i] = -1;
    uint64_t msk = (uint64_t)nsl - 1;
    for (int32_t k = 0; k < np0; k++) {
      uint64_t kk = pf_key(org0 + 3 * (size_t)k) + 1;
      uint64_t i = kk & msk;
      while (key[i] != 0 && key[i] != kk)
        i = (i + 1) & msk;
      key[i] = kk;
      val[i] = k;
    }
    int32_t nmatch = 0;
    for (int32_t k = 0; k < ps.np; k++) {
      uint64_t kk = pf_key(ps.p[k].org) + 1;
      uint64_t i = kk & msk;
      while (key[i] != 0 && key[i] != kk)
        i = (i + 1) & msk;
      if (key[i] == kk && val[i] >= 0) {
        for (int a = 0; a < 3; a++)
          tr.E[3 * (size_t)k + (size_t)a] = E0[3 * (size_t)val[i] + (size_t)a];
        nmatch++;
      }
    }
    printf("== ПОЛЕ ПЕРЕНЕСЕНО: совпало опорных точек %d из %d (%.2f %%)\n", nmatch, ps.np,
           100.0 * (double)nmatch / (double)(ps.np > 0 ? ps.np : 1));

    /* ДОПУСК ОСТАНОВА — ИМЕНОВАННЫЙ И ОДИН НА ОБА ПРОГОНА, иначе «тёплый
     * быстрее» получилось бы выбором допуска, а не свойством схемы. */
    const double TOLW = 1e-6;
    for (int pass = 0; pass < 2; pass++) {
      if (pass == 1)
        for (int32_t k = 0; k < 3 * ps.np; k++)
          tr.E[k] = 0.0; /* холодный: с нуля */
      int it = 0;
      double tw = now_s();
      for (; it < 200; it++) {
        hz_pstats st;
        memset(&st, 0, sizeof st);
        if (hz_psweep_bounce(&tr, &d, h, 0, sky, layout, &st) != 0) return 1;
        if (st.dE < TOLW) {
          it++;
          break;
        }
      }
      double dtw = now_s() - tw;
      double flux = 0.0;
      for (int32_t k = 0; k < ps.np; k++)
        flux += hz_ptrans_lout(&tr, k, 0.0, 0.0) * ps.p[k].area;
      printf("== %s СТАРТ: %d итераций за %.1f с до dE < %.0e; поток %.6e\n",
             (pass == 0) ? "ТЁПЛЫЙ" : "ХОЛОДНЫЙ", it, dtw, TOLW, flux);
      fflush(stdout);
    }
    free(key);
    free(val);
    free(org0);
    free(E0);
  }

  /* ------------------------- ПОЛЕ НА ДИСК И СВЕРКА С ЭТАЛОНОМ (§132) ------
   *
   * Поток — ИНТЕГРАЛ, и он к местным ошибкам слеп: при `h` от `0.5` до `2` м он
   * менялся на `0.01 %`, тогда как кромка тени размывается на два метра. §4
   * требует мерить радианс и ХВОСТ распределения, а не среднее и не картинку.
   *
   * Разбиение от `h` не зависит (сегментация идёт по мешу), поэтому номера
   * полигонов у эталона и у сравниваемого совпадают, и сверка идёт один в один.
   * Проценты берутся по ПЛОЩАДИ: вклад в изображение идёт площадью. */
  if (fsave != NULL) {
    FILE *fp = fopen(fsave, "wb");
    if (fp == NULL) return 1;
    int32_t n0 = ps.np;
    fwrite(&n0, sizeof n0, 1, fp);
    for (int32_t k = 0; k < ps.np; k++) {
      double b = hz_ptrans_lout(&tr, k, 0.0, 0.0);
      fwrite(&b, sizeof b, 1, fp);
    }
    fclose(fp);
    printf("== ПОЛЕ СОХРАНЕНО: %s, полигонов %d\n", fsave, ps.np);
  }
  if (fcmp != NULL) {
    FILE *fp = fopen(fcmp, "rb");
    int32_t n0 = 0;
    if (fp == NULL || fread(&n0, sizeof n0, 1, fp) != 1 || n0 != ps.np) {
      fprintf(stderr, "эталон %s не годится (полигонов %d, ожидалось %d)\n", fcmp, n0, ps.np);
      if (fp != NULL) fclose(fp);
      return 1;
    }
    double *ref = malloc((size_t)n0 * sizeof *ref);
    double *er = malloc((size_t)n0 * sizeof *er);
    double *aw = malloc((size_t)n0 * sizeof *aw);
    if (ref == NULL || er == NULL || aw == NULL ||
        fread(ref, sizeof *ref, (size_t)n0, fp) != (size_t)n0)
      return 1;
    fclose(fp);
    /* Нормировка на СРЕДНИЙ по площади радианс эталона, а не на местное
     * значение: у тёмного элемента относительная ошибка взрывается и хвост
     * начинает мерить деление на малое, а не расхождение полей. */
    double sa = 0.0, sb = 0.0;
    for (int32_t k = 0; k < n0; k++) {
      sa += ps.p[k].area;
      sb += ref[k] * ps.p[k].area;
    }
    double bref = (sa > 0.0) ? sb / sa : 1.0;
    int nn = 0;
    for (int32_t k = 0; k < n0; k++) {
      if (!(ps.p[k].area > 0.0)) continue;
      er[nn] = fabs(hz_ptrans_lout(&tr, k, 0.0, 0.0) - ref[k]) / (bref > 0.0 ? bref : 1.0);
      aw[nn] = ps.p[k].area;
      nn++;
    }
    /* Процентили ПО ПЛОЩАДИ: сортируем по ошибке, идём по накопленной площади. */
    for (int a = 1; a < nn; a++) {
      double e = er[a], w = aw[a];
      int b = a - 1;
      while (b >= 0 && er[b] > e) {
        er[b + 1] = er[b];
        aw[b + 1] = aw[b];
        b--;
      }
      er[b + 1] = e;
      aw[b + 1] = w;
    }
    double tot = 0.0;
    for (int a = 0; a < nn; a++)
      tot += aw[a];
    double q[3] = {0.5, 0.9, 0.99}, out[3] = {0, 0, 0};
    double acc2 = 0.0;
    int qi = 0;
    for (int a = 0; a < nn && qi < 3; a++) {
      acc2 += aw[a];
      while (qi < 3 && acc2 >= q[qi] * tot)
        out[qi++] = er[a];
    }
    printf("== СВЕРКА ПОЛЯ с %s (ошибка в долях среднего радианса %.4e, по ПЛОЩАДИ): "
           "p50 %.4f, p90 %.4f, p99 %.4f, максимум %.4f\n",
           fcmp, bref, out[0], out[1], out[2], er[nn - 1]);
    free(ref);
    free(er);
    free(aw);
  }

  /* ------------------------------------------------------------------ КАДР */
  {
    double eyeh[3] = HZ_CFG_HALL_EYE, ath[3] = HZ_CFG_HALL_AT;
    double eyec[3] = HZ_CFG_CITY_EYE, atc[3] = HZ_CFG_CITY_AT, up0[3] = HZ_CFG_UP;
    pf_cam cam;
    pf_cam_make(&cam, city ? eyec : eyeh, city ? atc : ath, up0, HZ_CFG_FOV_DEG, imgw, imgh);
    size_t npx = (size_t)imgw * (size_t)imgh;
    double *zb = malloc(npx * sizeof *zb);
    float *val = malloc(npx * sizeof *val);
    if (zb == NULL || val == NULL) return 1;
    for (size_t i = 0; i < npx; i++) {
      zb[i] = 1e300;
      val[i] = -1.0f;
    }
    double tcam = now_s();
    /* Треугольник -> полигон: значение берётся у ПОЛИГОНА, поле живёт на нём. */
    int32_t *t2p = malloc((size_t)m.nt * sizeof *t2p);
    if (t2p == NULL) return 1;
    for (int32_t t = 0; t < m.nt; t++)
      t2p[t] = -1;
    for (int32_t k = 0; k < ps.np; k++)
      for (int32_t t = ps.p[k].t0; t < ps.p[k].t0 + ps.p[k].ntri; t++)
        t2p[ps.tri[t]] = k;
    const double ZNEAR = 1e-3; /* м: ближе этого вершина уходит за камеру */
    for (int32_t t = 0; t < m.nt; t++) {
      int32_t k = t2p[t];
      if (k < 0) continue;
      /* Инициализация не косметика: анализатор не прослеживает запись через вызов
       * `pf_view`, а цикл отсечения ниже читает `vv` покоординатно. */
      double w[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
      double vv[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
      hz_obj_tri(&m, t, w);
      for (int i = 0; i < 3; i++)
        pf_view(&cam, w[i], vv[i]);
      /* Отсечение по ближней плоскости: до четырёх вершин, потом веер. */
      double cw[4][3], cwd[4][3];
      int nc = 0;
      for (int i = 0; i < 3 && nc < 4; i++) {
        int j = (i + 1) % 3;
        int in0 = vv[i][2] > ZNEAR, in1 = vv[j][2] > ZNEAR;
        if (in0) {
          memcpy(cw[nc], vv[i], sizeof cw[nc]);
          memcpy(cwd[nc], w[i], sizeof cwd[nc]);
          nc++;
        }
        if (in0 != in1 && nc < 4) {
          double a = (ZNEAR - vv[i][2]) / (vv[j][2] - vv[i][2]);
          for (int c = 0; c < 3; c++) {
            cw[nc][c] = vv[i][c] + a * (vv[j][c] - vv[i][c]);
            cwd[nc][c] = w[i][c] + a * (w[j][c] - w[i][c]);
          }
          nc++;
        }
      }
      if (nc < 3) continue;
      for (int f = 1; f + 1 < nc; f++) {
        const int id[3] = {0, f, f + 1};
        double sx[3] = {0, 0, 0}, sy[3] = {0, 0, 0}, iw[3] = {0, 0, 0};
        double wx[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
        for (int i = 0; i < 3; i++) {
          double z = cw[id[i]][2];
          iw[i] = 1.0 / z;
          sx[i] = cam.cx + cw[id[i]][0] * cam.fx * iw[i];
          sy[i] = cam.cy - cw[id[i]][1] * cam.fy * iw[i];
          for (int c = 0; c < 3; c++)
            wx[i][c] = cwd[id[i]][c] * iw[i];
        }
        double area = (sx[1] - sx[0]) * (sy[2] - sy[0]) - (sy[1] - sy[0]) * (sx[2] - sx[0]);
        if (!(fabs(area) > 0.0)) continue;
        double inv = 1.0 / area;
        double xlo = sx[0], xhi = sx[0], ylo = sy[0], yhi = sy[0];
        for (int i = 1; i < 3; i++) {
          if (sx[i] < xlo) xlo = sx[i];
          if (sx[i] > xhi) xhi = sx[i];
          if (sy[i] < ylo) ylo = sy[i];
          if (sy[i] > yhi) yhi = sy[i];
        }
        int x0 = (int)floor(xlo), x1 = (int)ceil(xhi), y0 = (int)floor(ylo), y1 = (int)ceil(yhi);
        if (x0 < 0) x0 = 0;
        if (y0 < 0) y0 = 0;
        if (x1 > imgw - 1) x1 = imgw - 1;
        if (y1 > imgh - 1) y1 = imgh - 1;
        for (int y = y0; y <= y1; y++)
          for (int x = x0; x <= x1; x++) {
            double px = (double)x + 0.5, py = (double)y + 0.5;
            double l0 = ((sx[1] - px) * (sy[2] - py) - (sy[1] - py) * (sx[2] - px)) * inv;
            double l1 = ((sx[2] - px) * (sy[0] - py) - (sy[2] - py) * (sx[0] - px)) * inv;
            double l2 = 1.0 - l0 - l1;
            if (l0 < 0.0 || l1 < 0.0 || l2 < 0.0) continue;
            double q = l0 * iw[0] + l1 * iw[1] + l2 * iw[2];
            if (!(q > 0.0)) continue;
            double z = 1.0 / q;
            size_t o = (size_t)y * (size_t)imgw + (size_t)x;
            if (!(z < zb[o])) continue;
            double wp[3];
            for (int c = 0; c < 3; c++)
              wp[c] = (l0 * wx[0][c] + l1 * wx[1][c] + l2 * wx[2][c]) / q;
            double du[3];
            for (int c = 0; c < 3; c++)
              du[c] = wp[c] - ps.p[k].org[c];
            double uu = du[0] * ps.p[k].eu[0] + du[1] * ps.p[k].eu[1] + du[2] * ps.p[k].eu[2];
            double vvl = du[0] * ps.p[k].ev[0] + du[1] * ps.p[k].ev[1] + du[2] * ps.p[k].ev[2];
            zb[o] = z;
            val[o] = (float)hz_ptrans_lout(&tr, k, uu, vvl);
          }
      }
    }
    double t_cam = now_s() - tcam;
    /* ТОНОВАЯ ШКАЛА НАЗЫВАЕТСЯ, А НЕ ПОДБИРАЕТСЯ МОЛЧА: делим на `p99` непустых
     * пикселей и берём гамму `2.2`. Небо (пустой пиксель) — синим. */
    double *srt = malloc(npx * sizeof *srt);
    if (srt == NULL) return 1;
    size_t ns = 0;
    for (size_t i = 0; i < npx; i++)
      if (val[i] >= 0.0f) srt[ns++] = (double)val[i];
    double p99 = 1.0;
    if (ns > 0) {
      for (size_t a = 1; a < ns; a++) { /* частичная сортировка не нужна: один раз */
        double kv = srt[a];
        size_t b = a;
        while (b > 0 && srt[b - 1] > kv) {
          srt[b] = srt[b - 1];
          b--;
        }
        srt[b] = kv;
      }
      p99 = srt[(size_t)(0.99 * (double)(ns - 1))];
    }
    if (!(p99 > 0.0)) p99 = 1.0;
    unsigned char *px = malloc(3 * npx);
    if (px == NULL) return 1;
    for (size_t i = 0; i < npx; i++) {
      if (val[i] < 0.0f) {
        px[3 * i] = 120;
        px[3 * i + 1] = 160;
        px[3 * i + 2] = 235;
        continue;
      }
      double u = (double)val[i] / p99;
      if (u > 1.0) u = 1.0;
      unsigned char g = (unsigned char)(255.0 * pow(u, 1.0 / 2.2));
      px[3 * i] = g;
      px[3 * i + 1] = g;
      px[3 * i + 2] = g;
    }
    char path[64];
    snprintf(path, sizeof path, "img/o46_flow_%s.ppm", city ? "city" : "hall");
    FILE *fp = fopen(path, "wb");
    if (fp != NULL) {
      fprintf(fp, "P6\n%d %d\n255\n", imgw, imgh);
      fwrite(px, 1, 3 * npx, fp);
      fclose(fp);
    }
    printf("== КАДР: %s, %d×%d, за %.2f с (растеризация полигонов, не луч); "
           "шкала: делено на p99 = %.4e, гамма 2.2; пустых пикселей %.1f %%\n",
           path, imgw, imgh, t_cam, p99, 100.0 * (double)(npx - ns) / (double)npx);
    free(px);
    free(srt);
    free(zb);
    free(val);
    free(t2p);
  }

  free(sacc);
  free(qn);
  free(fox);
  free(foy);
  free(foz);
  free(fw);
  tr3_dirs_free(&d);
  hz_ptrans_free(&tr);
  hz_poly_free(&ps);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
