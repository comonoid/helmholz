/* prender — Ш6: ТЕНЕВЫЕ РАЗРЕЗЫ, КАРТИНКА И ЕЁ ЦЕНА (PLAN_ELEMENTS.md).
 *
 * ЧТО СРАВНИВАЕТСЯ И ПОЧЕМУ ИМЕННО ТАК.
 *
 * Приёмка Ш6 дословно: «при РАВНОМ числе полигонов картинка лучше, чем без
 * разрезов». Значит сравнивать надо не «с разрезами против без», а ДВА
 * представления ОДНОЙ цены: (a) грубое `δ` плюс теневые разрезы, (b) мелкое
 * `δ` без разрезов, подобранное так, чтобы полигонов вышло столько же. Иначе
 * «лучше» означало бы просто «дороже».
 *
 * ЭТАЛОН НЕ ЗАВИСИТ ОТ ПОЛИГОНОВ, и это главное решение стенда. Поле полигона
 * есть `E = E_прям + E_косв`, и обе части считаются ПОРОЗНЬ (прямая — один раз,
 * она не меняется от отскока к отскоку). Тогда:
 *
 *     проверяемое:  L = L_e + ρ·(E_прям(u,v) + E_косв(u,v))/π
 *     эталон:       L = L_e + ρ·(E_прям_ТОЧНО(x) + E_косв(u,v))/π
 *
 * где `E_прям_ТОЧНО` берётся ЗАМКНУТОЙ ФОРМОЙ в самой точке попадания луча с
 * частой пробой видимости. Разница — ровно то, что разрезы и чинят: способность
 * ПОЛИГОНАЛЬНОГО поля представить резкую границу тени. Косвенная часть в обеих
 * ветвях одна и та же и из сравнения выпадает; эталон «тот же код при мелком δ»
 * этого бы не дал — он спорил бы сам с собой (ошибка, пойманная в Ш5 на
 * эталоне `ND = 256`).
 *
 * МЕТРИКА НА РАДИАНСЕ (К19), КРОМКИ ОТДЕЛЬНО (К20): кромкой считается пиксель,
 * у которого сосед попал в ДРУГОЙ полигон. ВРЕМЯ КАДРА И ТРАФИК докладываются
 * всегда — без них выигрыш остаётся выигрышем в счётчике, а не в секундах.
 */

#include "image.h"
#include "lodio.h"
#include "pcut.h"
#include "pmerge.h"
#include "pdirect.h"
#include "poly_seg.h"
#include "polygon.h"
#include "pray.h"
#include "psweep.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include "transport/cam3.h"
#include "transport/dirs3.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* Проб видимости в ЭТАЛОНЕ. Переменная, а не константа: у ОГРУБЛЁННОЙ сцены
 * полигоны крупные, сетка лучей отсекает плохо, и эталон дорожает в разы. */
static int NVIS_REF = 6;
/* ЭТАЛОН СЧИТАЕТСЯ НЕ ВСЕГДА (§78). Он нужен МЕТРИКЕ — сравнить поле полигонов с
 * точным прямым светом в точке попадания, — и стоит `пиксели × источники × пробы`,
 * то есть на городе больше самого решения. Для КАРТИНКИ он не нужен вовсе.
 * Величина `0` не «отключает проверку», а говорит: этот прогон даёт изображение, а
 * не замер, и метрики в нём не будет. */
static int WANT_REF = 1;
/* Проб видимости в ПРОВЕРЯЕМОМ прогоне. Параметр, а не константа, и это
 * вынужденно: при 2×2 = 4 пробах полутень квантуется на пять ступеней, то есть
 * площадной источник §2 ведёт себя как четыре точечных. Тогда выигрыш разрезов
 * мог бы оказаться выигрышем над СОБСТВЕННОЙ ступенькой выборки, а не над
 * физикой. Свип по этому числу — обязательная часть приёмки Ш6. */
static int NVIS_RUN = 2;
/* Потолок выборки по приёмнику, точек на ось (§80). Значение выбирается ЗАМЕРОМ:
 * свип по нему обязан показать, с какого места оно перестаёт влиять. */
static int SAMP_CAP = 64;
/* РАДИАНС НЕБА В СВИПЕ (§80, вопрос пользователя «зачем там пробы видимости»).
 *
 * Пробы видимости — цена ОБХОДНОГО МАНЁВРА, а не свойство механизма. В свипе по
 * ординатам заслонение возникает само: свет идёт вдоль направления и гасится
 * геометрией, попарной видимости там нет как объекта. `pdirect` вынесен ИЗ
 * ординат по замеру Ш5 — компактный источник углового размера `0.15` рад виден
 * лишь вдоль немногих ординат, и при `ND = 256` это давало 11 % площади сцены с
 * ошибкой выше 10 %. Вынесли — пришлось доставать заслонение отдельно.
 *
 * К НЕБУ ЭТОТ ДОВОД НЕ ОТНОСИТСЯ: его угловой размер порядка `π`, то есть полная
 * противоположность компактному, и ординаты представляют его отлично. Поэтому у
 * города источник задаётся ФОНОМ свипа, а не площадкой через `pdirect`: ни одной
 * пробы видимости, ни одного теневого луча — тени получаются из самого свипа. */
static double SKY_L = 0.0;
/* Кадр МЕТРИКИ — 512², а не 1024² из §2. Эталон стоит (пиксели × источники ×
 * 36 лучей), и на 1024² один прогон занимает минуты, а их в свипе двадцать.
 * Для КАРТИНКИ НА ГЛАЗ это не годится, и картинка пишется отдельно в полном
 * разрешении §2; для ХВОСТА распределения 262 тысячи точек достаточно. */
static int IMGW = 512;
static int IMGH = 512;

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

typedef struct {
  hz_objmesh m; /* сетка ПОСЛЕ разрезов (или копия исходной) */
  hz_pseglist sg;
  hz_polyset ps;
  hz_ptrans t;
  hz_pray g;
  double *srcLe; /* радианс источников; в развёртке они НЕ излучают (К9) */
  double *Edir;  /* 3·np: поле ТОЛЬКО прямого света */
  int32_t src[16];
  int nsrc, ncut;
  /* ЦВЕТ И ДЕТАЛЬ БЕРУТСЯ С МЕЛКОГО ПРЕДСТАВЛЕНИЯ (§78, правило Т1). Поле `E`
   * считается на ОГРУБЛЁННЫХ элементах — это дорого и потому огрублено; альбедо
   * применяется при ЧТЕНИИ, поэтому его можно брать с какой угодно мелкой сетки
   * даром. Раньше картинка брала `rho` ЭЛЕМЕНТА, то есть среднее по слитым
   * полигонам, и вся материальная деталь стиралась вместе с геометрией — хотя
   * стираться она не обязана вовсе. */
  const hz_pray *gfine;
  const double *albf; /* 3·np мелкого набора: RGB-альбедо */
  double t_solve, t_direct, t_frame;
  int64_t nfrag, nray, nover;
  int nbounce;
} scene;

static void scene_free(scene *S) {
  free(S->srcLe);
  free(S->Edir);
  hz_pray_free(&S->g);
  hz_ptrans_free(&S->t);
  hz_poly_free(&S->ps);
  hz_seg_free(&S->sg);
  hz_obj_free(&S->m);
  memset(S, 0, sizeof *S);
}

/* Источники §2. `shift` сдвигает сетку — это НЕГАТИВНЫЙ КОНТРОЛЬ: разрезы,
 * построенные под другое положение источника. */
static int add_lamps(hz_polyset *ps, const double lo[3], const double hi[3], int32_t *out, int nsrc,
                     double shift) {
  const double n[3] = {0.0, -1.0, 0.0}, eu[3] = {1.0, 0.0, 0.0};
  int k = 0;
  for (int i = 0; i < HZ_CFG_LAMP_NX && k < nsrc; i++)
    for (int j = 0; j < HZ_CFG_LAMP_NZ && k < nsrc; j++) {
      double c[3];
      c[0] = lo[0] + ((double)i + 0.5) / HZ_CFG_LAMP_NX * (hi[0] - lo[0]) + shift;
      c[1] = HZ_CFG_LAMP_Y;
      c[2] = lo[2] + ((double)j + 0.5) / HZ_CFG_LAMP_NZ * (hi[2] - lo[2]) + shift;
      out[k] = ps->np;
      if (hz_poly_add_quad(ps, c, n, eu, HZ_CFG_LAMP_SIDE / 2.0, HZ_CFG_LAMP_SIDE / 2.0, 0) != 0)
        return 2;
      k++;
    }
  return 0;
}

/* ИСТОЧНИК ДЛЯ ГОРОДА: ОДНА БОЛЬШАЯ ИЗЛУЧАЮЩАЯ ПЛОЩАДКА НАД СЦЕНОЙ (§78).
 *
 * Лампы под потолком для города бессмысленны, а честного неба у нас нет. Взята
 * простейшая вещь, которую УЖЕ умеет замкнутая форма прямого света: площадной
 * источник — квадрат размером со сцену, поднятый над ней и смотрящий вниз. Это
 * пасмурное небо в первом приближении: свет приходит со всей верхней полусферы,
 * теней с резким краем нет, и спорить о положении солнца не приходится.
 *
 * ВЫСОТА — ПОЛОВИНА ГОРИЗОНТАЛЬНОГО РАЗМЕРА СЦЕНЫ НАД ЕЁ ВЕРХОМ. Число не
 * подобрано под картинку: при такой высоте площадка видна из середины сцены под
 * углом около `90°`, то есть закрывает верхнюю полусферу примерно так, как её
 * закрывает небо. Ниже — площадка начнёт светить как потолок, выше — как точечное
 * солнце, и оба случая пришлось бы обосновывать отдельно.
 *
 * ЯРКОСТЬ та же, что у ламп зала: сравнивать картинки между сценами всё равно
 * нельзя, а масштаб радианса в метрике сокращается. */
static int add_sky(hz_polyset *ps, const double lo[3], const double hi[3], int32_t *out) {
  const double n[3] = {0.0, -1.0, 0.0}, eu[3] = {1.0, 0.0, 0.0};
  double c[3];
  double sx = 0.5 * (hi[0] - lo[0]), sz = 0.5 * (hi[2] - lo[2]);
  double span = (sx > sz) ? sx : sz;
  c[0] = 0.5 * (lo[0] + hi[0]);
  c[1] = hi[1] + span;
  c[2] = 0.5 * (lo[2] + hi[2]);
  out[0] = ps->np;
  return hz_poly_add_quad(ps, c, n, eu, sx, sz, 0);
}

/* КОНТРОЛЬ ДИАГНОЗА: дробление полигона ПРОСТРАНСТВЕННОЙ СЕТКОЙ шага `L`.
 *
 * Нужен затем, чтобы отличить две причины ошибки, которые силуэтные разрезы
 * различить не дают. Теневой разрез кладёт линию РАЗРЫВА поля; регулярная
 * сетка просто делает полигоны МЕЛЬЧЕ, не зная ни про свет, ни про силуэты.
 * Если ошибку убирает сетка, а разрезы — нет, значит дело не в резкости поля
 * на границе тени, а в том, что ЛИНЕЙНОЕ `E` не тянет ПЛАВНУЮ вариацию `1/r²`
 * по крупному плоскому полигону. Это разные болезни и разные лекарства.
 *
 * Треугольники здесь НЕ РЕЖУТСЯ: метка получает добавку по ячейке, в которую
 * попал центр треугольника. Край выходит зубчатым, но для вопроса «помогает ли
 * уменьшение полигона» это безразлично, а кода — двадцать строк вместо двухсот. */
/* `eye != NULL` — шаг сетки РАСТЁТ С ДАЛЬНОСТЬЮ: `L(R) = L·R`, то есть ровно
 * `L = εR` из CLAUDE.md. Это и есть LOD, которого в остальном замере нет вовсе,
 * а он тут ключевой: без него пол площадью 40 м² остаётся ОДНИМ полигоном,
 * находясь в 0.4…7 м от глаза, — то есть замер наказывает модель за случай,
 * который LOD обязан не допускать. Крупный полигон законен НА УДАЛЕНИИ. */
static int subdivide_grid_eye(hz_pseglist *so, const hz_objmesh *m, const hz_pseglist *si, double L,
                              const double *eye) {
  memset(so, 0, sizeof *so);
  /* При LOD ячейка мельчает у глаза, поэтому решётка заводится по САМОЙ мелкой
   * ячейке; ключ считается делением на местный шаг, а не индексом решётки. */
  int32_t nx = (int32_t)((m->hi[0] - m->lo[0]) / L) + 1;
  int32_t ny = (int32_t)((m->hi[1] - m->lo[1]) / L) + 1;
  int32_t nz = (int32_t)((m->hi[2] - m->lo[2]) / L) + 1;
  int64_t ncell = (int64_t)nx * ny * nz;
  if (eye != NULL) {
    /* угловых бинов: (π/ε) × (2π/ε) */
    ncell = ((int64_t)(M_PI / L) + 2) * ((int64_t)(2.0 * M_PI / L) + 2);
  }
  so->label = malloc((size_t)m->nt * sizeof *so->label);
  so->seg = malloc((size_t)m->nt * sizeof *so->seg);
  int32_t *map = malloc((size_t)si->nseg * (size_t)ncell * sizeof *map);
  if (so->label == NULL || so->seg == NULL || map == NULL) {
    free(map);
    hz_seg_free(so);
    return 2;
  }
  for (int64_t i = 0; i < (int64_t)si->nseg * ncell; i++)
    map[i] = -1;
  int32_t nn = 0;
  for (int32_t t = 0; t < m->nt; t++) {
    double p[3][3], c[3] = {0, 0, 0};
    hz_obj_tri(m, t, p);
    for (int i = 0; i < 3; i++)
      for (int a = 0; a < 3; a++)
        c[a] += p[i][a] / 3.0;
    int64_t cell;
    if (eye != NULL) {
      /* УГЛОВОЙ БИН, А НЕ ПРОСТРАНСТВЕННЫЙ. `L = εR` означает ПОСТОЯННЫЙ
       * УГЛОВОЙ РАЗМЕР ячейки, и биннинг по направлению из глаза даёт его
       * буквально. Пространственная решётка с шагом `L(R)` для этого не годится:
       * индекс, посчитанный делением на РАЗНЫЙ шаг, склеил бы ячейки на разной
       * дальности, попавшие в один номер, — то есть слил бы geometрически
       * далёкие куски в один полигон. */
      double q[3], R = 0.0;
      for (int a = 0; a < 3; a++) {
        q[a] = c[a] - eye[a];
        R += q[a] * q[a];
      }
      R = sqrt(R);
      if (!(R > 0.0)) R = 1e-9;
      double th = acos(q[2] / R), ph = atan2(q[1], q[0]) + M_PI;
      int64_t it = (int64_t)(th / L), ip = (int64_t)(ph / L);
      int64_t nph = (int64_t)(2.0 * M_PI / L) + 1;
      cell = it * nph + ip;
      if (cell < 0) cell = 0;
      if (cell >= ncell) cell = cell % ncell;
    } else {
      int32_t ix = (int32_t)((c[0] - m->lo[0]) / L), iy = (int32_t)((c[1] - m->lo[1]) / L),
              iz = (int32_t)((c[2] - m->lo[2]) / L);
      if (ix < 0) ix = 0;
      if (iy < 0) iy = 0;
      if (iz < 0) iz = 0;
      if (ix >= nx) ix = nx - 1;
      if (iy >= ny) iy = ny - 1;
      if (iz >= nz) iz = nz - 1;
      cell = ((int64_t)iz * ny + iy) * nx + ix;
    }
    int64_t key = (int64_t)si->label[t] * ncell + cell;
    if (map[key] < 0) {
      map[key] = nn;
      so->seg[nn] = si->seg[si->label[t]];
      so->seg[nn].area = 0.0;
      so->seg[nn].ntri = 0;
      so->seg[nn].dmax = 0.0;
      nn++;
    }
    int32_t g = map[key];
    so->label[t] = g;
    so->seg[g].ntri++;
    so->seg[g].area += hz_obj_tri_area(m, t);
    for (int i = 0; i < 3; i++) {
      double dv = fabs(p[i][0] * so->seg[g].n[0] + p[i][1] * so->seg[g].n[1] +
                       p[i][2] * so->seg[g].n[2] - so->seg[g].off);
      if (dv > so->seg[g].dmax) so->seg[g].dmax = dv;
    }
  }
  so->nseg = nn;
  so->delta = si->delta;
  free(map);
  return 0;
}

/* Сборка сцены: сегментация, при docut — клинья и измельчение, полигоны,
 * источники, сетка лучей. `lampshift_cut` задаёт положение источника, ПОД
 * КОТОРОЕ строятся разрезы; `lampshift_lit` — где источник светит на самом
 * деле. Они различаются только в негативном контроле. */
static int scene_build(scene *S, const hz_objmesh *base, double delta, int nsrc, int docut,
                       const hz_cutcfg *cfg, double lampshift_cut, double lampshift_lit,
                       hz_wedges *wout, double gridL) {
  memset(S, 0, sizeof *S);
  hz_pseglist sg0;
  if (hz_seg_planar(&sg0, base, delta) != 0) return 1;

  const hz_objmesh *msrc = base;
  hz_pseglist *ssrc = &sg0;
  hz_objmesh mcut;
  hz_pseglist scut;
  int cut_made = 0;

  if (docut == 2) {
    /* КОНТРОЛЬ: дробление пространственной сеткой, без силуэтов и без света. */
    if (subdivide_grid_eye(&scut, base, &sg0, gridL, NULL) != 0) return 1;
    memset(&mcut, 0, sizeof mcut);
    msrc = base;
    ssrc = &scut;
    cut_made = 2;
    S->ncut = 0;
  } else if (docut == 4) {
    /* ЧИСТОЕ дробление ПЛОСКОСТЯМИ: край прямой, осколков нет. */
    int64_t ov = 0;
    if (hz_cut_grid(&mcut, &scut, base, &sg0, gridL, &ov) != 0) return 1;
    msrc = &mcut;
    ssrc = &scut;
    cut_made = 1;
    S->ncut = 0;
  } else if (docut == 3) {
    /* КОНТРОЛЬ LOD: шаг растёт с дальностью от глаза, L = εR. */
    double eye[3] = HZ_CFG_HALL_EYE;
    if (subdivide_grid_eye(&scut, base, &sg0, gridL, eye) != 0) return 1;
    memset(&mcut, 0, sizeof mcut);
    msrc = base;
    ssrc = &scut;
    cut_made = 2;
    S->ncut = 0;
  } else if (docut) {
    /* Клинья строятся на ИСХОДНОЙ геометрии с источниками на месте. */
    hz_polyset p0;
    if (hz_poly_build(&p0, base, &sg0) != 0) return 1;
    int32_t s0[16];
    if (add_lamps(&p0, base->lo, base->hi, s0, nsrc, lampshift_cut) != 0) return 1;
    double *Le0 = calloc((size_t)p0.np, sizeof *Le0);
    if (Le0 == NULL) return 1;
    for (int i = 0; i < nsrc; i++)
      Le0[s0[i]] = HZ_CFG_LAMP_LE;
    int wrc = hz_wedges_build(wout, base, &p0, s0, nsrc, Le0, cfg);
    free(Le0);
    if (wrc != 0) return 1;
    hz_poly_free(&p0);
    if (hz_cut_apply(&mcut, &scut, base, &sg0, wout, cfg->budget) != 0) return 1;
    msrc = &mcut;
    ssrc = &scut;
    cut_made = 1;
    S->ncut = (cfg->budget < wout->n) ? cfg->budget : wout->n;
  }

  if (hz_poly_build(&S->ps, msrc, ssrc) != 0) return 1;
  if (add_lamps(&S->ps, base->lo, base->hi, S->src, nsrc, lampshift_lit) != 0) return 1;
  S->nsrc = nsrc;
  if (hz_ptrans_init(&S->t, &S->ps, msrc) != 0) return 1;
  S->srcLe = calloc((size_t)S->ps.np, sizeof *S->srcLe);
  S->Edir = calloc((size_t)S->ps.np * 3, sizeof *S->Edir);
  if (S->srcLe == NULL || S->Edir == NULL) return 1;
  for (int i = 0; i < nsrc; i++) {
    S->srcLe[S->src[i]] = HZ_CFG_LAMP_LE;
    S->t.rho[S->src[i]] = 0.0;
    S->t.Le[S->src[i]] = 0.0; /* в развёртке источник НЕ излучает: К9 */
  }
  if (hz_pray_build(&S->g, &S->ps, 4.0) != 0) return 1;

  /* Сетка и разметка нужны только для сборки полигонов; дальше живут копии
   * внутри `S->ps`. Держим ту, из которой собирали, чтобы материалы были живы. */
  if (cut_made == 2) {
    S->sg = scut;
    memset(&S->m, 0, sizeof S->m);
    hz_seg_free(&sg0);
  } else if (cut_made) {
    S->m = mcut;
    S->sg = scut;
    hz_seg_free(&sg0);
  } else {
    S->sg = sg0;
    memset(&S->m, 0, sizeof S->m);
  }
  return 0;
}

/* Сцена из ГОТОВОЙ разметки: нужна огрублению, которое сперва решает задачу на
 * мелком представлении, а потом строит из него грубое. */
static int scene_from(scene *S, const hz_objmesh *m, const hz_pseglist *sg, const double lo[3],
                      const double hi[3], int nsrc) {
  memset(S, 0, sizeof *S);
  if (hz_poly_build(&S->ps, m, sg) != 0) return 1;
  /* `nsrc < 0` — ГОРОД: одна площадка-небо вместо решётки ламп (§78). */
  if (nsrc < 0) {
    if (add_sky(&S->ps, lo, hi, S->src) != 0) return 1;
    nsrc = 1;
  } else if (add_lamps(&S->ps, lo, hi, S->src, nsrc, 0.0) != 0)
    return 1;
  S->nsrc = nsrc;
  if (hz_ptrans_init(&S->t, &S->ps, m) != 0) return 1;
  S->srcLe = calloc((size_t)S->ps.np, sizeof *S->srcLe);
  S->Edir = calloc((size_t)S->ps.np * 3, sizeof *S->Edir);
  if (S->srcLe == NULL || S->Edir == NULL) return 1;
  for (int i = 0; i < nsrc; i++) {
    S->srcLe[S->src[i]] = HZ_CFG_LAMP_LE;
    S->t.rho[S->src[i]] = 0.0;
    S->t.Le[S->src[i]] = 0.0;
  }
  if (hz_pray_build(&S->g, &S->ps, 4.0) != 0) return 1;
  return 0;
}

/* Решение: прямой свет один раз (он не меняется), развёртка — до застоя. */
static int solve(scene *S, const tr3_dirs *d, double h, double tol) {
  hz_pstats st;
  memset(&st, 0, sizeof st);
  double t0 = now_s();
  hz_psweep_zero(&S->t);
  /* При небе-фоне прямой свет через `pdirect` не считается вовсе: источника-полигона
   * нет, и считать нечего. */
  if (S->nsrc > 0)
    hz_direct_add(&S->t, &S->g, S->src, S->nsrc, S->srcLe, NVIS_RUN, h, SAMP_CAP, &st);
  double *accdir = malloc((size_t)S->t.np * 3 * sizeof *accdir);
  if (accdir == NULL) return 1;
  memcpy(accdir, S->t.acc, (size_t)S->t.np * 3 * sizeof *accdir);
  S->nray += st.nfrag;
  /* Поле ТОЛЬКО прямого света — нужно эталону, чтобы вычесть его и заменить. */
  hz_psweep_solve(&S->t, &st);
  memcpy(S->Edir, S->t.E, (size_t)S->t.np * 3 * sizeof *S->Edir);
  S->t_direct = now_s() - t0;

  double t1 = now_s();
  int nb = 0;
  for (; nb < 200; nb++) {
    memset(&st, 0, sizeof st);
    memcpy(S->t.acc, accdir, (size_t)S->t.np * 3 * sizeof *accdir);
    if (hz_psweep_gather(&S->t, d, h, 0, SKY_L, HZ_LAYOUT_RUNS, &st) != 0) {
      free(accdir);
      return 1;
    }
    hz_psweep_solve(&S->t, &st);
    S->nfrag += st.nfrag;
    S->nover += st.nover;
    if (st.dE < tol) break;
  }
  free(accdir);
  S->t_solve = now_s() - t1;
  S->nbounce = nb + 1;
  return 0;
}

/* Точный прямой свет в точке — эталон. `*cls` получает класс точки:
 *   0 — освещена целиком (каждый источник виден весь),
 *   1 — ПОЛУТЕНЬ (хоть один источник виден частично),
 *   2 — тень (ни один источник не виден).
 *
 * КЛАСС НУЖЕН НЕ ДЛЯ КРАСОТЫ, А ЧТОБЫ РАЗВЕСТИ ДВЕ ПРИЧИНЫ ХВОСТА. Ш5 намерил
 * плато `p99 ≈ 3…4%`, не убираемое ни `ND`, ни шагом растра, и я приписал его
 * резкости поля на границе тени — ГИПОТЕЗА, причина не была изолирована.
 * Конкурирующее объяснение измерено в Ш4: у ЗАВЕДОМО ОДНОРОДНОГО поля печи
 * появляется паразитный градиент `0.7…10%` от `E` из-за несимметричности
 * пиксельной выборки по полигону. Это не разрыв, и разрезами не лечится.
 * Если хвост сидит в теле картинки так же, как в полутени, — причина растровая,
 * и вывод Ш6 обязан быть отрицательным. */
static double direct_exact(const scene *S, int32_t k, const double x[3], int *cls) {
  const hz_poly *P = &S->ps.p[k];
  double eps = 1e-7 * (1.0 + fabs(x[0]) + fabs(x[1]) + fabs(x[2]));
  double xo[3] = {0.0, 0.0, 0.0};
  for (int a = 0; a < 3; a++)
    xo[a] = x[a] + eps * P->n[a];
  double Ed = 0.0;
  int nfull = 0, npart = 0, nzero = 0, nabove = 0;
  for (int s = 0; s < S->nsrc; s++) {
    double Eu = hz_direct_unocc(&S->ps, S->src[s], S->srcLe[S->src[s]], x, P->n);
    if (!(Eu > 0.0)) continue;
    nabove++;
    const hz_poly *Q = &S->ps.p[S->src[s]];
    int vis = 0;
    for (int a = 0; a < NVIS_REF; a++)
      for (int b = 0; b < NVIS_REF; b++) {
        double su = Q->uvlo[0] + ((double)a + 0.5) / NVIS_REF * (Q->uvhi[0] - Q->uvlo[0]);
        double sv = Q->uvlo[1] + ((double)b + 0.5) / NVIS_REF * (Q->uvhi[1] - Q->uvlo[1]);
        double sp[3], dd[3], len = 0.0;
        hz_poly_world(Q, su, sv, sp);
        for (int c = 0; c < 3; c++) {
          dd[c] = sp[c] - xo[c];
          len += dd[c] * dd[c];
        }
        len = sqrt(len);
        if (!(len > 0.0)) continue;
        for (int c = 0; c < 3; c++)
          dd[c] /= len;
        if (!hz_pray_occluded(&S->g, xo, dd, 0.0, len * (1.0 - 1e-6))) vis++;
      }
    if (vis == 0)
      nzero++;
    else if (vis == NVIS_REF * NVIS_REF)
      nfull++;
    else
      npart++;
    Ed += Eu * (double)vis / (double)(NVIS_REF * NVIS_REF);
  }
  if (cls != NULL) *cls = (npart > 0) ? 1 : ((nfull > 0 || nabove == 0) ? 0 : 2);
  return Ed;
}

typedef struct {
  double p50, p99, pmax, p50e, p99e;
  double frac1, frac10;
  /* ХВОСТ ПО КЛАССАМ ТОЧКИ — разводит две причины (см. direct_exact) */
  double p99_lit, p99_pen, p99_shd;
  int n_lit, n_pen, n_shd;
  int npix, nedge;
  double t_frame;
} imgstat;

static int cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

/* ЗЕРКАЛЬНЫЙ И СТЕКЛЯННЫЙ ШАРЫ (§79, просьба пользователя: показать, что по
 * возможностям это сравнимо с трассировкой лучей).
 *
 * ЧТО ЭТО ЕСТЬ И ЧЕГО НЕ ЕСТЬ — СКАЗАНО ДО КАРТИНКИ, А НЕ ПОСЛЕ. Решатель хранит
 * ОБЛУЧЁННОСТЬ `E` на полигоне и читает `L = ρ·E/π`; направленности в этой
 * величине нет вовсе, поэтому зеркало в неё не помещается в принципе. Здесь
 * сделан ГИБРИД: диффузный свет берётся из решённого поля, а зеркальность и
 * преломление — трассировкой НА ЭТАПЕ КАМЕРЫ. Луч отражается или преломляется,
 * летит в сцену и читает там уже посчитанное `L`.
 *
 * ЧЕГО В НЁМ НЕТ, И ЭТО ВИДНО НА КАРТИНКЕ: шары не отдают свет обратно в сцену и
 * НЕ ОТБРАСЫВАЮТ ТЕНЕЙ — решение о них не знает. Честный путь к тому и другому
 * один: направленный радианс на элементе, то есть ординаты и BRDF с Френелем, —
 * это веха переноса, а не вечер работы.
 *
 * Преломление сферы считается двумя пересечениями (вход и выход) с законом
 * Снеллиуса на каждой границе; доля отражения — Френель по неполяризованному
 * свету (приближение Шлика здесь не нужно, точная формула столь же дёшева). */
typedef struct {
  double c[3], r;
  int kind;   /* 0 — зеркало, 1 — стекло */
  double ior; /* показатель преломления стекла */
} gball;

static int g_nball = 0;
static gball g_ball[2];

/* Ближайшее пересечение луча со сферой при `t > tmin`. */
static double ball_hit(const gball *b, const double o[3], const double d[3], double tmin) {
  double oc[3], B = 0.0, C = 0.0, dd = 0.0;
  for (int i = 0; i < 3; i++) {
    oc[i] = o[i] - b->c[i];
    B += oc[i] * d[i];
    C += oc[i] * oc[i];
    dd += d[i] * d[i];
  }
  C -= b->r * b->r;
  double disc = B * B - dd * C;
  if (disc <= 0.0) return -1.0;
  double sq = sqrt(disc);
  double t1 = (-B - sq) / dd, t2 = (-B + sq) / dd;
  if (t1 > tmin) return t1;
  if (t2 > tmin) return t2;
  return -1.0;
}

/* Диффузный радианс сцены в точке попадания луча — то самое ЧТЕНИЕ поля. */
static void shade_diffuse(const scene *S, int32_t k, const double o[3], const double d[3], double t,
                          double out[3]) {
  const hz_poly *P = &S->ps.p[k];
  double q[3];
  for (int a = 0; a < 3; a++)
    q[a] = o[a] + t * d[a] - P->org[a];
  double u = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
  double v = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
  const double *ce = S->t.E + (size_t)k * 3;
  double E = ce[0] + ce[1] * u + ce[2] * v;
  if (E < 0.0) E = 0.0;
  double alb[3] = {S->t.rho[k], S->t.rho[k], S->t.rho[k]};
  if (S->gfine != NULL && S->albf != NULL) {
    double tf = 0.0;
    int32_t kf = hz_pray_hit(S->gfine, o, d, 0.0, &tf);
    if (kf >= 0)
      for (int c = 0; c < 3; c++)
        alb[c] = S->albf[(size_t)kf * 3 + (size_t)c];
  }
  for (int c = 0; c < 3; c++)
    out[c] = S->srcLe[k] + alb[c] * E / M_PI;
}

/* Трассировка с шарами. `depth` — оставшиеся отскоки; `0` обрывает рекурсию. */
static void trace_rgb(const scene *S, const double o[3], const double d[3], int depth,
                      double out[3]) {
  out[0] = out[1] = out[2] = 0.0;
  double tp = 0.0;
  int32_t k = hz_pray_hit(&S->g, o, d, 0.0, &tp);
  int hb = -1;
  double tb = 1e300;
  for (int i = 0; i < g_nball; i++) {
    double t = ball_hit(&g_ball[i], o, d, 1e-6);
    if (t > 0.0 && t < tb) {
      tb = t;
      hb = i;
    }
  }
  if (hb >= 0 && (k < 0 || tb < tp)) {
    if (depth <= 0) return;
    const gball *b = &g_ball[hb];
    double x[3], nn[3];
    double nl = 0.0;
    for (int c = 0; c < 3; c++) {
      x[c] = o[c] + tb * d[c];
      nn[c] = (x[c] - b->c[c]) / b->r;
      nl += nn[c] * nn[c];
    }
    (void)nl;
    double dl = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    double dir[3];
    for (int c = 0; c < 3; c++)
      dir[c] = d[c] / dl;
    double cosi = -(dir[0] * nn[0] + dir[1] * nn[1] + dir[2] * nn[2]);
    double refl[3], xo[3];
    for (int c = 0; c < 3; c++) {
      refl[c] = dir[c] + 2.0 * cosi * nn[c];
      xo[c] = x[c] + 1e-5 * nn[c];
    }
    if (b->kind == 0) { /* зеркало: одно отражение, поглощения 5 % */
      double r3[3];
      trace_rgb(S, xo, refl, depth - 1, r3);
      for (int c = 0; c < 3; c++)
        out[c] = 0.95 * r3[c];
      return;
    }
    /* СТЕКЛО. Френель по неполяризованному свету, вход и выход считаются порознь. */
    double n1 = 1.0, n2 = b->ior;
    double eta = n1 / n2;
    double k2 = 1.0 - eta * eta * (1.0 - cosi * cosi);
    double R = 1.0;
    if (k2 > 0.0) {
      double cost = sqrt(k2);
      double rs = (n1 * cosi - n2 * cost) / (n1 * cosi + n2 * cost);
      double rp = (n1 * cost - n2 * cosi) / (n1 * cost + n2 * cosi);
      R = 0.5 * (rs * rs + rp * rp);
      /* Преломлённый луч внутрь, второе пересечение — выход наружу. */
      double ti[3], xin[3];
      for (int c = 0; c < 3; c++) {
        ti[c] = eta * dir[c] + (eta * cosi - cost) * nn[c];
        xin[c] = x[c] - 1e-5 * nn[c];
      }
      double t2 = ball_hit(b, xin, ti, 1e-6);
      double th[3] = {0.0, 0.0, 0.0};
      if (t2 > 0.0) {
        double y[3], n2v[3];
        for (int c = 0; c < 3; c++) {
          y[c] = xin[c] + t2 * ti[c];
          n2v[c] = (b->c[c] - y[c]) / b->r; /* нормаль ИЗНУТРИ наружу */
        }
        double cosi2 = -(ti[0] * n2v[0] + ti[1] * n2v[1] + ti[2] * n2v[2]);
        double eta2 = n2 / n1;
        double kk = 1.0 - eta2 * eta2 * (1.0 - cosi2 * cosi2);
        if (kk > 0.0) {
          double cost2 = sqrt(kk), to[3], yo[3];
          for (int c = 0; c < 3; c++) {
            to[c] = eta2 * ti[c] + (eta2 * cosi2 - cost2) * n2v[c];
            yo[c] = y[c] - 1e-5 * n2v[c];
          }
          trace_rgb(S, yo, to, depth - 1, th);
        }
      }
      double r3[3];
      trace_rgb(S, xo, refl, depth - 1, r3);
      for (int c = 0; c < 3; c++)
        out[c] = R * r3[c] + (1.0 - R) * th[c];
      return;
    }
    double r3[3]; /* полное внутреннее отражение */
    trace_rgb(S, xo, refl, depth - 1, r3);
    for (int c = 0; c < 3; c++)
      out[c] = r3[c];
    return;
  }
  if (k < 0) return;
  shade_diffuse(S, k, o, d, tp, out);
}

/* ЦВЕТНОЙ КАДР ДЛЯ ГЛАЗА (§78). Отдельный путь от `render()` СОЗНАТЕЛЬНО: тот
 * считает метрику и эталон, и мешать в него тон-маппинг с суперсэмплингом значит
 * менять замер ради красоты. Здесь наоборот — ни метрики, ни эталона, только
 * изображение.
 *
 * ТРИ ОТЛИЧИЯ ОТ ПРЕЖНЕЙ КАРТИНКИ, И ВСЕ ТРИ БЫЛИ ПРОСТО НЕ СДЕЛАНЫ:
 *   1. ЦВЕТ. `hz_ppm_write` гнал яркость через ЛОЖНУЮ РАСКРАСКУ (`colormap`), то
 *      есть фиолетово-оранжевую шкалу для чтения полей. Для глаза нужен честный
 *      серый с гаммой, а лучше — цвет материала.
 *   2. АЛЬБЕДО С МЕЛКОГО ПРЕДСТАВЛЕНИЯ, а не с элемента (правило Т1).
 *   3. СУПЕРСЭМПЛИНГ: край элемента — ступенька в один пиксель, и без сглаживания
 *      она видна на любой сцене.
 * Тон-маппинг: деление на 99.5-й процентиль яркости и гамма `2.2`. Процентиль, а
 * не максимум, — иначе одна яркая лампа гасит весь кадр. */
static int render_rgb(scene *S, const tr3_camera *cam, int ss, const char *ppm) {
  const int W = cam->w, H = cam->h;
  int n = W * H;
  double *L3 = calloc((size_t)n * 3, sizeof *L3);
  unsigned char *rgb = malloc((size_t)n * 3);
  if (L3 == NULL || rgb == NULL) {
    free(L3);
    free(rgb);
    return 1;
  }
  if (ss < 1) ss = 1;
#pragma omp parallel for schedule(dynamic, 16)
  for (int i = 0; i < n; i++) {
    int px = i % W, py = i / W;
    double acc[3] = {0.0, 0.0, 0.0};
    for (int sy = 0; sy < ss; sy++)
      for (int sx = 0; sx < ss; sx++) {
        double o[3], d[3], t;
        /* Подпиксель: камера умеет только целые пиксели, поэтому смещение
         * вносится долей пикселя в НАПРАВЛЕНИИ — этого довольно для сглаживания
         * края и не требует правки оператора камеры. */
        tr3_camera_ray(cam, px, py, o, d);
        if (ss > 1) {
          /* Шаг пикселя по обеим осям берётся у самой камеры — соседними лучами,
           * а не выводом её внутренней формулы: так подпиксельное смещение
           * остаётся верным при любом операторе камеры. */
          double ox[3], dx[3], oy[3], dy[3];
          tr3_camera_ray(cam, (px + 1 < W) ? px + 1 : px - 1, py, ox, dx);
          tr3_camera_ray(cam, px, (py + 1 < H) ? py + 1 : py - 1, oy, dy);
          double fx = ((double)sx + 0.5) / (double)ss - 0.5;
          double fy = ((double)sy + 0.5) / (double)ss - 0.5;
          if (px + 1 >= W) fx = -fx;
          if (py + 1 >= H) fy = -fy;
          for (int c = 0; c < 3; c++)
            d[c] += fx * (dx[c] - d[c]) + fy * (dy[c] - d[c]);
          double dn = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
          if (dn > 0.0)
            for (int c = 0; c < 3; c++)
              d[c] /= dn;
        }
        if (g_nball > 0) {
          /* С шарами кадр идёт через трассировку: она сама решает, что ближе —
           * шар или геометрия, — и рекурсивно читает поле за отражением. */
          double c3[3];
          trace_rgb(S, o, d, 4, c3);
          for (int c = 0; c < 3; c++)
            acc[c] += c3[c];
          continue;
        }
        int32_t k = hz_pray_hit(&S->g, o, d, 0.0, &t);
        if (k < 0) continue;
        const hz_poly *P = &S->ps.p[k];
        double q[3] = {0.0, 0.0, 0.0};
        for (int a = 0; a < 3; a++)
          q[a] = o[a] + t * d[a] - P->org[a];
        double u = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
        double v = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
        const double *ce = S->t.E + (size_t)k * 3;
        double Etot = ce[0] + ce[1] * u + ce[2] * v;
        if (Etot < 0.0) Etot = 0.0;
        double alb[3] = {S->t.rho[k], S->t.rho[k], S->t.rho[k]};
        if (S->gfine != NULL && S->albf != NULL) {
          double tf = 0.0;
          int32_t kf = hz_pray_hit(S->gfine, o, d, 0.0, &tf);
          if (kf >= 0)
            for (int c = 0; c < 3; c++)
              alb[c] = S->albf[(size_t)kf * 3 + (size_t)c];
        }
        for (int c = 0; c < 3; c++)
          acc[c] += S->srcLe[k] + alb[c] * Etot / M_PI;
      }
    double w = 1.0 / (double)(ss * ss);
    for (int c = 0; c < 3; c++)
      L3[(size_t)i * 3 + (size_t)c] = acc[c] * w;
  }
  /* Нормировка по 99.5-му процентилю ЯРКОСТИ, общая для трёх каналов: канальная
   * нормировка увела бы цвет. */
  double *lum = malloc((size_t)n * sizeof *lum);
  if (lum == NULL) {
    free(L3);
    free(rgb);
    return 1;
  }
  for (int i = 0; i < n; i++)
    lum[i] = 0.2126 * L3[3 * i] + 0.7152 * L3[3 * i + 1] + 0.0722 * L3[3 * i + 2];
  qsort(lum, (size_t)n, sizeof *lum, cmp_d);
  double mx = lum[(size_t)((double)n * 0.995)];
  free(lum);
  if (!(mx > 0.0)) mx = 1.0;
  for (int i = 0; i < 3 * n; i++) {
    double x = L3[i] / mx;
    if (x < 0.0) x = 0.0;
    if (x > 1.0) x = 1.0;
    x = pow(x, 1.0 / 2.2);
    rgb[i] = (unsigned char)(255.0 * x + 0.5);
  }
  int rc = hz_ppm_write_rgb(ppm, rgb, W, H);
  free(L3);
  free(rgb);
  return rc;
}

/* Кадр: луч в пиксель, радианс из поля; заодно эталон и метрика. */
static int render(scene *S, const tr3_camera *cam, imgstat *out, const char *ppm, const char *pfm) {
  const int W = cam->w, H = cam->h;
  int n = W * H;
  double *Lt = malloc((size_t)n * sizeof *Lt), *Lr = malloc((size_t)n * sizeof *Lr);
  int32_t *who = malloc((size_t)n * sizeof *who);
  int8_t *klass = malloc((size_t)n);
  if (Lt == NULL || Lr == NULL || who == NULL || klass == NULL) {
    free(Lt);
    free(Lr);
    free(who);
    free(klass);
    return 1;
  }

  double t0 = now_s();
#pragma omp parallel for schedule(dynamic, 16)
  for (int i = 0; i < n; i++) {
    double o[3], d[3], t;
    tr3_camera_ray(cam, i % W, i / W, o, d);
    int32_t k = hz_pray_hit(&S->g, o, d, 0.0, &t);
    who[i] = k;
    if (k < 0) {
      Lt[i] = 0.0;
      Lr[i] = 0.0;
      continue;
    }
    const hz_poly *P = &S->ps.p[k];
    double q[3] = {0.0, 0.0, 0.0};
    for (int a = 0; a < 3; a++)
      q[a] = o[a] + t * d[a] - P->org[a];
    double u = q[0] * P->eu[0] + q[1] * P->eu[1] + q[2] * P->eu[2];
    double v = q[0] * P->ev[0] + q[1] * P->ev[1] + q[2] * P->ev[2];
    /* Инициализация ПРИ ОБЪЯВЛЕНИИ, как у `q` двумя строками выше: без неё
     * `gcc -fanalyzer` теряет заполнение цикла и сообщает о чтении
     * неинициализированного в `direct_exact` (CWE-457). Находка старше правки
     * О22, но файл тронут — значит гейт на нём мой. */
    double x[3] = {0.0, 0.0, 0.0};
    for (int a = 0; a < 3; a++)
      x[a] = o[a] + t * d[a];

    const double *ce = S->t.E + (size_t)k * 3;
    const double *cd = S->Edir + (size_t)k * 3;
    double Etot = ce[0] + ce[1] * u + ce[2] * v;
    double Edf = cd[0] + cd[1] * u + cd[2] * v;
    double Eind = Etot - Edf;
    if (Eind < 0.0) Eind = 0.0;
    int cls = 0;
    double Ede = WANT_REF ? direct_exact(S, k, x, &cls) : 0.0;
    klass[i] = (int8_t)cls;
    if (Etot < 0.0) Etot = 0.0;
    Lt[i] = S->srcLe[k] + S->t.rho[k] * Etot / M_PI;
    Lr[i] = S->srcLe[k] + S->t.rho[k] * (Ede + Eind) / M_PI;
  }
  double t1 = now_s();

  /* --- метрика: тело картинки и КРОМКИ порознь (К20) --- */
  double mean = 0.0;
  int nm = 0;
  for (int i = 0; i < n; i++)
    if (who[i] >= 0) {
      mean += Lr[i];
      nm++;
    }
  mean = (nm > 0) ? mean / nm : 1.0;
  if (!(mean > 0.0)) mean = 1.0;

  double *body = malloc((size_t)n * sizeof *body), *edge = malloc((size_t)n * sizeof *edge);
  double *cl[3];
  int ncl[3] = {0, 0, 0};
  for (int c = 0; c < 3; c++)
    cl[c] = malloc((size_t)n * sizeof *cl[c]);
  if (body == NULL || edge == NULL || cl[0] == NULL || cl[1] == NULL || cl[2] == NULL) {
    free(body);
    free(edge);
    for (int c = 0; c < 3; c++)
      free(cl[c]);
    free(Lt);
    free(Lr);
    free(who);
    free(klass);
    return 1;
  }
  int nb = 0, ne = 0, n1 = 0, n10 = 0;
  for (int j = 0; j < H; j++)
    for (int i = 0; i < W; i++) {
      int p = j * W + i;
      if (who[p] < 0) continue;
      int is_edge = 0;
      if (i > 0 && who[p - 1] != who[p]) is_edge = 1;
      if (i < W - 1 && who[p + 1] != who[p]) is_edge = 1;
      if (j > 0 && who[p - W] != who[p]) is_edge = 1;
      if (j < H - 1 && who[p + W] != who[p]) is_edge = 1;
      double e = fabs(Lt[p] - Lr[p]) / mean;
      if (is_edge)
        edge[ne++] = e;
      else {
        body[nb++] = e;
        if (e > 0.01) n1++;
        if (e > 0.10) n10++;
        int c = klass[p];
        if (c >= 0 && c < 3) cl[c][ncl[c]++] = e;
      }
    }
  for (int c = 0; c < 3; c++)
    qsort(cl[c], (size_t)ncl[c], sizeof *cl[c], cmp_d);
  qsort(body, (size_t)nb, sizeof *body, cmp_d);
  qsort(edge, (size_t)ne, sizeof *edge, cmp_d);
  memset(out, 0, sizeof *out);
  out->npix = nb;
  out->nedge = ne;
  if (nb > 0) {
    out->p50 = body[nb / 2];
    out->p99 = body[(int)(nb * 0.99)];
    out->pmax = body[nb - 1];
    out->frac1 = 100.0 * n1 / nb;
    out->frac10 = 100.0 * n10 / nb;
  }
  if (ne > 0) {
    out->p50e = edge[ne / 2];
    out->p99e = edge[(int)(ne * 0.99)];
  }
  out->n_lit = ncl[0];
  out->n_pen = ncl[1];
  out->n_shd = ncl[2];
  if (ncl[0] > 0) out->p99_lit = cl[0][(int)(ncl[0] * 0.99)];
  if (ncl[1] > 0) out->p99_pen = cl[1][(int)(ncl[1] * 0.99)];
  if (ncl[2] > 0) out->p99_shd = cl[2][(int)(ncl[2] * 0.99)];
  out->t_frame = t1 - t0;
  free(body);
  free(edge);
  for (int c = 0; c < 3; c++)
    free(cl[c]);
  free(klass);

  if (ppm != NULL) hz_ppm_write(ppm, Lt, W, H);
  if (pfm != NULL) hz_pfm_write(pfm, Lt, Lt, Lt, W, H);
  free(Lt);
  free(Lr);
  free(who);
  return 0;
}

/* Один замер: сборка сцены, решение, кадр, метрика. */
static int one(const hz_objmesh *base, const tr3_camera *cam, const tr3_dirs *d, double delta,
               double h, int nsrc, int docut, const hz_cutcfg *cfg, double shift_cut,
               const char *label, const char *img, imgstat *out, int32_t *np_out, int *ncut_out,
               double gridL) {
  hz_wedges ws;
  memset(&ws, 0, sizeof ws);
  scene S;
  if (scene_build(&S, base, delta, nsrc, docut, cfg, shift_cut, 0.0, &ws, gridL) != 0) return 1;
  /* Допуск по К88: ДОЛЯ от оценённой ошибки дискретизации (она здесь ~1e−2),
   * а не абсолютная константа. 1e−6 стоило бы вчетверо больше отскоков и не
   * изменило бы ни одной цифры метрики. */
  if (solve(&S, d, h, 1e-4) != 0) return 1;
  /* Нулями ПРИ ОБЪЯВЛЕНИИ: заполняются только при , и связи между
   * этим условием и  ниже cppcheck не видит (legacyUninitvar).
   * Находка старше правки О22, но файл тронут. */
  char ppm[96] = {0}, pfm[96] = {0};
  if (img != NULL) {
    snprintf(ppm, sizeof ppm, "img/%s.ppm", img);
    snprintf(pfm, sizeof pfm, "img/%s.pfm", img);
  }
  if (render(&S, cam, out, img ? ppm : NULL, img ? pfm : NULL) != 0) return 1;
  printf("   %-30s %7d %7d %5d %9.3e %9.3e %8.2f%% %8.2f%% %9.3e %9.3e %7.2f %7.2f\n", label,
         S.ps.np, docut ? S.ncut : 0, S.nbounce, out->p50, out->p99, out->frac1, out->frac10,
         out->p50e, out->p99e, S.t_direct + S.t_solve, out->t_frame);
  /* ПОЛИГОНЫ БЕЗ КРАЯ считаются и печатаются ВСЕГДА. `hz_poly_build` отбрасывает
   * петли короче трёх вершин, и при мелком дроблении осколок вполне может не
   * дать ни одной годной петли: он остаётся в наборе своими треугольниками, но
   * для ЛУЧА его нет — `hz_poly_inside` по пустому краю всегда ложь. Это
   * ровно тот подозреваемый, на котором Ш6 списал отказ шага 0.05 м на потерю
   * геометрии; потеря вылечена (`ПОТЕРЯНО 0`), а отказ остался, значит диагноз
   * был неполон, и величину надо предъявлять, а не подразумевать. */
  int32_t noloop = 0;
  for (int32_t k = 0; k < S.ps.np; k++)
    if (S.ps.p[k].nloop == 0) noloop++;
  printf("        по классам точки  p99: свет %.3e (%d), ПОЛУТЕНЬ %.3e (%d), тень %.3e (%d)"
         "   клиньев-кандидатов %d, ПОТЕРЯНО фрагментов %lld, полигонов БЕЗ КРАЯ %d (%.2f%%)\n",
         out->p99_lit, out->n_lit, out->p99_pen, out->n_pen, out->p99_shd, out->n_shd, ws.n,
         (long long)S.nover, noloop, 100.0 * noloop / (double)(S.ps.np > 0 ? S.ps.np : 1));
  if (np_out != NULL) *np_out = S.ps.np;
  if (ncut_out != NULL) *ncut_out = S.ncut;
  fflush(stdout);
  hz_wedges_free(&ws);
  scene_free(&S);
  return 0;
}

/* δ, при котором БЕЗ разрезов выходит примерно `target` полигонов. Ищется одной
 * сегментацией на шаг: копировать собранную сцену нельзя (в `hz_pray` живёт
 * указатель на `hz_polyset`), да и незачем. */
static double match_delta(const hz_objmesh *base, double dhi, int32_t target, int nsrc) {
  double dlo = dhi * 0.01, best = dhi;
  int32_t bestnp = -1;
  for (int it = 0; it < 9; it++) {
    double dm = sqrt(dlo * dhi);
    hz_pseglist probe;
    if (hz_seg_planar(&probe, base, dm) != 0) return best;
    int32_t np = probe.nseg + nsrc;
    hz_seg_free(&probe);
    if (bestnp < 0 || labs((long)np - (long)target) < labs((long)bestnp - (long)target)) {
      bestnp = np;
      best = dm;
    }
    if (np < target)
      dhi = dm;
    else
      dlo = dm;
  }
  return best;
}

int main(int argc, char **argv) {
  double delta = (argc > 1) ? strtod(argv[1], NULL) : 0.045;
  double h = (argc > 2) ? strtod(argv[2], NULL) : 0.02;
  double eps = (argc > 3) ? strtod(argv[3], NULL) : 0.045;
  double tol = (argc > 4) ? strtod(argv[4], NULL) : 0.05;
  if (argc > 5) NVIS_RUN = atoi(argv[5]);

  /* ВЫБОР СЦЕНЫ (§78). До сих пор `prender` держал зал жёстко, и городская
   * картинка была недостижима не по существу, а по обвязке. */
  int city = 0;
  for (int i = 1; i < argc; i++)
    if (strcmp(argv[i], "city") == 0) city = 1;
  hz_objmesh base;
  if (hz_obj_load(&base, city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ,
                  city ? HZ_CFG_CITY_SCALE : HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  printf("== СЦЕНА: габарит (%.1f %.1f %.1f) … (%.1f %.1f %.1f), треугольников %d\n", base.lo[0],
         base.lo[1], base.lo[2], base.hi[0], base.hi[1], base.hi[2], base.nt);
  tr3_camera cam;
  double eyeh[3] = HZ_CFG_HALL_EYE, ath[3] = HZ_CFG_HALL_AT;
  double eyec[3] = HZ_CFG_CITY_EYE, atc[3] = HZ_CFG_CITY_AT;
  double eye[3], at[3], up[3] = HZ_CFG_UP;
  for (int c = 0; c < 3; c++) {
    eye[c] = city ? eyec[c] : eyeh[c];
    at[c] = city ? atc[c] : ath[c];
  }
  /* КАМЕРА ПЕРЕОПРЕДЕЛЯЕТСЯ С КОМАНДНОЙ СТРОКИ (§80). Городская камера из
   * конфигурации стоит на земле и смотрит горизонтально — она заводилась для
   * МЕТРИКИ, где важен разброс дальностей, а не для показа города. Положение
   * камеры на срез не влияет иначе как через дальности: `ε` лестницы сверяется по
   * разрешению кадра, а не по точке съёмки. */
  for (int i = 1; i < argc; i++) {
    if (strncmp(argv[i], "eye=", 4) == 0)
      sscanf(argv[i] + 4, "%lf,%lf,%lf", &eye[0], &eye[1], &eye[2]);
    if (strncmp(argv[i], "at=", 3) == 0) sscanf(argv[i] + 3, "%lf,%lf,%lf", &at[0], &at[1], &at[2]);
  }
  if (tr3_camera_look(&cam, eye, at, up, HZ_CFG_FOV_DEG * M_PI / 180.0, IMGW, IMGH) != 0) return 1;
  tr3_dirs d;
  if (tr3_dirs_product(&d, 4, 4) != 0) return 1;

  printf("== Ш6: зал, δ = %g м, h = %g м, ε = %g м, tol = %g, ND = %d, метрика на %d×%d\n", delta,
         h, eps, tol, d.n, IMGW, IMGH);
  printf("   %-30s %7s %7s %5s %9s %9s %9s %9s %9s %9s %7s %7s\n", "конфигурация", "полиг", "клинь",
         "отск", "p50", "p99", "S>1%", "S>10%", "кром p50", "кром p99", "реш,с", "кадр,с");

  /* РЕЖИМ «ТОЛЬКО LOD»: РАЗРЕЗЫ ПРОТИВ ДЕТАЛИЗАЦИИ, лоб в лоб.
   *
   * Вопрос узкий и он последний: единственное, чем разрезы били любую
   * светонезависимую сетку, — КРАЙНИЙ хвост (`p99 = 1.19` против `1.43`).
   * Если детализация по `L = εR` закрывает и его, теневые разрезы не нужны
   * вовсе, и правило «нужен разрез — возьми более детальный уровень»
   * заменяет весь §1.8 целиком. */
  if (argc > 6 && argv[6][0] == 'l') {
    hz_cutcfg cfg0 = {eps, tol, 0}, cfg1 = {eps, tol, 1024};
    imgstat s;
    if (one(&base, &cam, &d, delta, h, 8, 0, &cfg0, 0.0, "без разрезов, без LOD", NULL, &s, NULL,
            NULL, 0.0) != 0)
      return 1;
    if (one(&base, &cam, &d, delta, h, 8, 1, &cfg1, 0.0, "РАЗРЕЗЫ, бюджет 1024", NULL, &s, NULL,
            NULL, 0.0) != 0)
      return 1;
    for (int i = 0; i < 4; i++) {
      double L = (double[]){0.20, 0.10, 0.05, 0.025}[i];
      char lab[64], im[64];
      snprintf(lab, sizeof lab, "ЗУБЧАТОЕ L = %.3f·R", L);
      snprintf(im, sizeof im, "sh6_lod_%d", i);
      if (one(&base, &cam, &d, delta, h, 8, 3, &cfg0, 0.0, lab, im, &s, NULL, NULL, L) != 0)
        return 1;
    }
    /* СВИП ДРОБЛЕНИЯ ПЛОСКОСТЯМИ ДОВЕДЁН НИЖЕ 0.10 м (§10, открытый пункт 4).
     * Геометрия там уже проверена — площадь сохраняется ТОЧНО до шага 0.025 м
     * (обход по ЯЧЕЙКАМ вместо резки плоскостями по очереди), — а качество на
     * этих шагах снято не было, и сходимость `p99` до 0.10 м могла оказаться
     * началом плато, как это уже один раз вышло с зубчатым краем. */
    for (int i = 0; i < 6; i++) {
      double L = (double[]){1.20, 0.40, 0.20, 0.10, 0.05, 0.025}[i];
      char lab[64], im[64];
      snprintf(lab, sizeof lab, "ПЛОСКОСТЯМИ, шаг %.3f м", L);
      snprintf(im, sizeof im, "sh6_grid_%d", i);
      if (one(&base, &cam, &d, delta, h, 8, 4, &cfg0, 0.0, lab, im, &s, NULL, NULL, L) != 0)
        return 1;
    }
    tr3_dirs_free(&d);
    hz_obj_free(&base);
    return 0;
  }

  /* --- О22: РЕНДЕР СРЕЗА LOD (§62, пункт 3). Аргумент 7 — файл лестницы. ---
   *
   * ЧТО СРАВНИВАЕТСЯ. Три ветви одной сцены: НУЛЕВОЙ УРОВЕНЬ (эталон цены),
   * СРЕЗ по камере и ОДНОРОДНЫЕ УРОВНИ лестницы. Однородные берутся ГОТОВЫМИ, а
   * не подбором `δ` под равное число элементов (А148): число элементов от `δ`
   * немонотонно, подбор не гарантирован и стоит минуты на прогон, а два соседних
   * готовых уровня отвечают на тот же вопрос без единого лишнего слияния.
   *
   * ОШИБКА РАДИАНСА СЧИТАЕТСЯ ПРОТИВ ЗАМКНУТОЙ ФОРМЫ, А НЕ ПРОТИВ НУЛЕВОГО
   * УРОВНЯ, и это решение Ш6, а не удобство: эталон, зависящий от полигонов,
   * спорил бы сам с собой. Поэтому строка нулевого уровня — тоже ЗАМЕР со своей
   * ошибкой, а не ноль по определению.
   *
   * НЕГАТИВНЫЙ КОНТРОЛЬ — ПЕРЕСТАНОВКА УРОВНЕЙ (А146), а не случайный предок:
   * случайный предок меняет ЦЕНУ, и провал вышел бы по причине, к критерию среза
   * отношения не имеющей. Перестановка сохраняет набор выбранных уровней и рвёт
   * только связь «грубее там, где дальше». */
  if (argc > 6 && argv[6][0] == 'o') {
    const char *lodf = NULL;
    for (int i = 7; i < argc; i++)
      if (strncmp(argv[i], "lod=", 4) == 0) lodf = argv[i] + 4;
    if (lodf == NULL) lodf = city ? "build/lod/city_r.lod" : "build/lod/hall_r.lod";
    /* КАРТИНКА, А НЕ ЗАМЕР: эталон стоит `пиксели × источники × пробы` и на городе
     * дороже самого решения, а изображению не нужен вовсе (§78). Включается
     * словом `ref`, и тогда прогон становится замером. */
    WANT_REF = 0;
    for (int i = 7; i < argc; i++)
      if (strcmp(argv[i], "ref") == 0) WANT_REF = 1;
    NVIS_REF = 4;
    /* РАЗРЕШЕНИЕ И ПРОБЫ ВИДИМОСТИ — параметрами (§78). Умолчание 256² и четыре
     * пробы годились для ЗАМЕРА (их выбирала цена эталона), но для картинки дают
     * ступенчатый край и квантованную полутень: площадной источник при четырёх
     * пробах ведёт себя как четыре точечных. */
    IMGW = 512;
    for (int i = 7; i < argc; i++) {
      if (strncmp(argv[i], "w=", 2) == 0) IMGW = (int)strtol(argv[i] + 2, NULL, 10);
      if (strncmp(argv[i], "vis=", 4) == 0) NVIS_RUN = (int)strtol(argv[i] + 4, NULL, 10);
      if (strncmp(argv[i], "cap=", 4) == 0) SAMP_CAP = (int)strtol(argv[i] + 4, NULL, 10);
      if (strncmp(argv[i], "sky=", 4) == 0) SKY_L = strtod(argv[i] + 4, NULL);
    }
    if (IMGW < 64) IMGW = 64;
    IMGH = IMGW;
    if (tr3_camera_look(&cam, eye, at, up, HZ_CFG_FOV_DEG * M_PI / 180.0, IMGW, IMGH) != 0)
      return 1;
    const double eps_px = (HZ_CFG_FOV_DEG * M_PI / 180.0) / (double)IMGH;
    hz_pseglist sgf;
    if (hz_seg_planar(&sgf, &base, delta) != 0) return 1;
    hz_lod L;
    int fsimp = 0;
    int lrc = hz_lod_read(&L, lodf, &sgf, base.nt, &fsimp);
    if (lrc != HZ_LODIO_OK) {
      fprintf(stderr, "лестница %s не прочитана: %s (код %d)\n", lodf, hz_lodio_str(lrc), lrc);
      return 1;
    }
    /* УПРОЩЕНИЕ КРАЯ ЗДЕСЬ НЕ ПРИМЕНЯЕТСЯ, значит лестница обязана быть построена
     * БЕЗ него. Отказ, а не «примерно то же»: границы полигонов от него зависят, а
     * сумма разметки его не видит (оговорка в `lodio.h`). */
    if (fsimp) {
      fprintf(stderr, "лестница построена С УПРОЩЕНИЕМ края, а рендер его не делает — отказ\n");
      return 1;
    }
    if (!(fabs(L.eps - eps_px) <= 1e-12 * (fabs(eps_px) + 1.0))) {
      fprintf(stderr, "лестница под ε = %.9g, а кадр даёт %.9g — отказ (А131)\n", L.eps, eps_px);
      return 1;
    }
    printf("\n== О22: РЕНДЕР СРЕЗА LOD, лестница %s (уровней %d, узлов %d, ε = %.6g)\n", lodf,
           L.nlev, L.nnd, L.eps);
    printf("   %-30s %7s %9s %9s %8s %8s %7s\n", "ветвь", "элем", "p50", "p99", "S>1%", "S>10%",
           "кадр,с");
    fflush(stdout);

    /* ШАРЫ (§79). Ставятся у точки прицела камеры, радиус — от размера сцены,
     * чтобы не подбирать числа под кадр: `1/12` наибольшего горизонтального
     * размера даёт шар, занимающий заметную, но не всю картинку. */
    for (int i = 7; i < argc; i++) {
      if (strncmp(argv[i], "ball=", 5) != 0) continue;
      const char *w = argv[i] + 5;
      double sx = base.hi[0] - base.lo[0], sz = base.hi[2] - base.lo[2];
      double rr = ((sx > sz) ? sx : sz) / 16.0;
      int mir = (strcmp(w, "mirror") == 0 || strcmp(w, "both") == 0);
      int gls = (strcmp(w, "glass") == 0 || strcmp(w, "both") == 0);
      /* Центр — НЕ у точки прицела: там стена, и шары в неё упирались. Берётся
       * точка на 45 % пути от глаза к прицелу и поднимается — это середина
       * открытого пространства кадра при любой камере, а не подобранные числа. */
      /* НАД СТОЛОМ (замечание пользователя: «лучше всего, когда они над столом,
       * видно хорошо»). Точка прицела камеры и есть стол; `0.8` пути от глаза
       * ставит шары над ним и при этом не вплотную к дальней стене, в которую они
       * упирались при `1.0`. Подъём — полтора радиуса над точкой прицела. */
      double ctr[3];
      for (int c = 0; c < 3; c++)
        ctr[c] = eye[c] + 0.8 * (at[c] - eye[c]);
      double up0 = rr * 1.5;
      if (mir) {
        g_ball[g_nball].c[0] = ctr[0];
        g_ball[g_nball].c[1] = ctr[1] + up0;
        g_ball[g_nball].c[2] = ctr[2] - (gls ? rr * 1.25 : 0.0);
        g_ball[g_nball].r = rr;
        g_ball[g_nball].kind = 0;
        g_ball[g_nball].ior = 1.0;
        g_nball++;
      }
      if (gls) {
        g_ball[g_nball].c[0] = ctr[0];
        g_ball[g_nball].c[1] = ctr[1] + up0;
        g_ball[g_nball].c[2] = ctr[2] + (mir ? rr * 1.25 : 0.0);
        g_ball[g_nball].r = rr;
        g_ball[g_nball].kind = 1;
        g_ball[g_nball].ior = 1.5;
        g_nball++;
      }
      printf("== ШАРЫ: %d, радиус %.3f м, центр у точки прицела\n", g_nball, rr);
    }
    /* ТОЛЬКО СРЕЗ (§79, оснастка). Ветвь гоняла ВСЕ уровни лестницы: нулевой,
     * срез, каждый однородный и контроль — то есть семнадцать полных решений
     * свипа при пятнадцати уровнях. Для ЗАМЕРА это и нужно, для КАРТИНКИ — нет, и
     * на городе разница между «сорок минут» и «двенадцать часов». */
    int only_cut = 0;
    for (int i = 7; i < argc; i++)
      if (strcmp(argv[i], "onlycut") == 0) only_cut = 1;
    int ssaa = 2;
    for (int i = 7; i < argc; i++)
      if (strncmp(argv[i], "ss=", 3) == 0) ssaa = (int)strtol(argv[i] + 3, NULL, 10);
    int32_t *sel = malloc((size_t)L.np * sizeof *sel);
    int32_t *lev_of = malloc((size_t)L.np * sizeof *lev_of);
    if (sel == NULL || lev_of == NULL) return 1;

    /* МЕЛКОЕ ПРЕДСТАВЛЕНИЕ ДЛЯ ЦВЕТА (§78, правило Т1): строится ОДИН раз и живёт
     * всю ветвь. Поле берётся с огрублённых элементов, альбедо — отсюда, и потому
     * материальная деталь не зависит от того, насколько огрублена геометрия. */
    hz_polyset psfine;
    hz_pray gfine;
    double *albf = NULL;
    int havefine = 0;
    if (hz_poly_build(&psfine, &base, &sgf) == 0) {
      if (hz_pray_build(&gfine, &psfine, 4.0) == 0) {
        albf = malloc((size_t)psfine.np * 3 * sizeof *albf);
        if (albf != NULL) {
          for (int32_t k = 0; k < psfine.np; k++) {
            /* Материал берётся у ПЕРВОГО треугольника полигона: сегментация
             * плоская, и внутри участка материал, как правило, один. */
            int32_t t = (psfine.p[k].ntri > 0) ? psfine.tri[psfine.p[k].t0] : -1;
            int32_t mi = (t >= 0 && base.fm != NULL) ? base.fm[t] : 0;
            for (int c = 0; c < 3; c++)
              albf[(size_t)k * 3 + (size_t)c] =
                  (mi >= 0 && mi < base.nmtl) ? base.mtl[mi].kd3[c] : 0.5;
          }
          havefine = 1;
        }
      }
    }

    /* Ветвь: по выбору узлов `sel` построить сцену, решить, снять кадр. */
    for (int pass = only_cut ? 1 : 0; pass < (only_cut ? 2 : L.nlev + 2); pass++) {
      char tag[64], img[96];
      if (pass == 0) {
        /* НУЛЕВОЙ УРОВЕНЬ. */
        for (int32_t k = 0; k < L.np; k++)
          sel[k] = L.lab[k];
        snprintf(tag, sizeof tag, "уровень 0 (мелкий)");
        snprintf(img, sizeof img, "img/o22_%s_lev0.ppm", city ? "city" : "hall");
      } else if (pass == 1) {
        /* СРЕЗ. */
        /* Критерий по ПОЛНОМУ смещению (§76.4): множитель sin отпускал допуск для
         * поверхностей, видимых в лоб, и стоил порядка по качеству. */
        if (hz_lod_cut(&L, eye, L.eps, 1, sel) <= 0) return 1;
        for (int32_t k = 0; k < L.np; k++)
          lev_of[k] = L.nd[sel[k]].level;
        snprintf(tag, sizeof tag, "СРЕЗ по камере");
        snprintf(img, sizeof img, "img/o22_%s_cut.ppm", city ? "city" : "hall");
      } else {
        int lev = pass - 1;
        if (lev >= L.nlev) break;
        for (int32_t k = 0; k < L.np; k++)
          sel[k] = L.lab[(size_t)lev * (size_t)L.np + (size_t)k];
        snprintf(tag, sizeof tag, "однородный уровень %d", lev);
        snprintf(img, sizeof img, "img/o22_%s_uni%d.ppm", city ? "city" : "hall", lev);
      }
      hz_pseglist so;
      if (hz_lod_seglist(&L, &base, &sgf, sel, &so) != 0) return 1;
      scene C;
      /* `nsrc = 0` — источников-полигонов нет вовсе: город освещён ФОНОМ свипа. */
      if (scene_from(&C, &base, &so, base.lo, base.hi, city ? ((SKY_L > 0.0) ? 0 : -1) : 8) != 0)
        return 1;
      if (solve(&C, &d, h, 1e-4) != 0) return 1;
      if (havefine) {
        C.gfine = &gfine;
        C.albf = albf;
      }
      /* РАЗБИВКА ПО ФАЗАМ (§79, оснастка). «Рендер идёт пятнадцать минут» — жалоба,
       * а не адрес: тот же урок, что в `pcoarse`. Без этих трёх чисел я гадал,
       * где время — в геометрии, в прямом свете или в свипе, — и один раз уже
       * угадал неверно (винил шаг растра). */
      printf("      фазы, с: прямой свет %.1f, свип %.1f (отскоков %d), элементов %d\n", C.t_direct,
             C.t_solve, C.nbounce, C.ps.np);
      fflush(stdout);
      imgstat sc;
      if (render(&C, &cam, &sc, NULL, NULL) != 0) return 1;
      /* Цветной кадр — отдельным путём, суперсэмплинг 2×2 (§78). */
      if (render_rgb(&C, &cam, ssaa, img) != 0) return 1;
      printf("   %-30s %7d %9.3e %9.3e %8.2f%% %8.2f%% %7.2f\n", tag, C.ps.np, sc.p50, sc.p99,
             sc.frac1, sc.frac10, sc.t_frame);
      fflush(stdout);
      scene_free(&C);
      hz_seg_free(&so);
    }

    if (only_cut) {
      free(sel);
      free(lev_of);
      if (havefine) {
        hz_pray_free(&gfine);
        hz_poly_free(&psfine);
        free(albf);
      }
      hz_lod_free(&L);
      hz_seg_free(&sgf);
      tr3_dirs_free(&d);
      hz_obj_free(&base);
      return 0;
    }

    /* НЕГАТИВНЫЙ КОНТРОЛЬ: те же уровни, розданные ДРУГИМ полигонам (А146).
     * Перестановка детерминированная — хешем от номера, а не генератором: замер
     * обязан повторяться, и «случайно» здесь означало бы «неповторимо». */
    {
      for (int32_t k = 0; k < L.np; k++) {
        uint64_t x = (uint64_t)(uint32_t)k * 0x9e3779b97f4a7c15ull;
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ull;
        x ^= x >> 27;
        int32_t j = (int32_t)(x % (uint64_t)(uint32_t)L.np);
        int32_t t = lev_of[k];
        lev_of[k] = lev_of[j];
        lev_of[j] = t;
      }
      for (int32_t k = 0; k < L.np; k++) {
        int32_t lv = lev_of[k];
        if (lv < 0) lv = 0;
        if (lv >= L.nlev) lv = L.nlev - 1;
        sel[k] = L.lab[(size_t)lv * (size_t)L.np + (size_t)k];
      }
      hz_pseglist so;
      if (hz_lod_seglist(&L, &base, &sgf, sel, &so) != 0) return 1;
      scene C;
      /* `nsrc = 0` — источников-полигонов нет вовсе: город освещён ФОНОМ свипа. */
      if (scene_from(&C, &base, &so, base.lo, base.hi, city ? ((SKY_L > 0.0) ? 0 : -1) : 8) != 0)
        return 1;
      if (solve(&C, &d, h, 1e-4) != 0) return 1;
      imgstat sc;
      if (render(&C, &cam, &sc, "img/o22_perm.ppm", NULL) != 0) return 1;
      printf("   %-30s %7d %9.3e %9.3e %8.2f%% %8.2f%% %7.2f\n", "КОНТРОЛЬ: уровни перемешаны",
             C.ps.np, sc.p50, sc.p99, sc.frac1, sc.frac10, sc.t_frame);
      scene_free(&C);
      hz_seg_free(&so);
    }
    free(sel);
    free(lev_of);
    if (havefine) {
      hz_pray_free(&gfine);
      hz_poly_free(&psfine);
      free(albf);
    }
    hz_lod_free(&L);
    hz_seg_free(&sgf);
    tr3_dirs_free(&d);
    hz_obj_free(&base);
    return 0;
  }

  /* --- Ш7: ОГРУБЛЕНИЕ. Критерий ГЕОМЕТРИЧЕСКИЙ и от света НЕ зависит. --- */
  if (argc > 6 && argv[6][0] == 0x63) {
    /* Эталон удешевлён ВЧЕТВЕРО по пикселям и ВПЯТЕРО по пробам: у огрублённой
     * сцены он иначе дороже самого решения, а хвост распределения 65 тысяч
     * точек держат. */
    NVIS_REF = 4;
    IMGW = 256;
    IMGH = 256;
    if (tr3_camera_look(&cam, eye, at, up, HZ_CFG_FOV_DEG * M_PI / 180.0, IMGW, IMGH) != 0)
      return 1;
    hz_pseglist sgf;
    if (hz_seg_planar(&sgf, &base, delta) != 0) return 1;
    scene F;
    if (scene_from(&F, &base, &sgf, base.lo, base.hi, 8) != 0) return 1;
    if (solve(&F, &d, h, 1e-4) != 0) return 1;
    imgstat s0;
    if (render(&F, &cam, &s0, NULL, NULL) != 0) return 1;
    printf("   %-34s %7d %9.3e %9.3e %8.2f%% %8.2f%%\n", "исходное (мелкое)", F.ps.np, s0.p50,
           s0.p99, s0.frac1, s0.frac10);
    fflush(stdout);

    /* СЕТКА КАНДИДАТОВ РАВНА ПЕРЕБОРУ — проверяется, а не утверждается.
     * Сравниваются РАЗМЕТКИ треугольников целиком: совпадение числа участков
     * ничего не значило бы, потому что то же число достижимо разными
     * слияниями. Порядок пар при равной ошибке доопределён (`cmp_pair`),
     * поэтому совпадение обязано быть ТОЧНЫМ, а не «в пределах». */
    {
      hz_mergecfg mg;
      memset(&mg, 0, sizeof mg);
      mg.delta = eps;
      mg.target = 400;
      mg.use_geom = 1;
      mg.use_overlap = 1;
      hz_pseglist s1, s2;
      hz_mergestat t1, t2;
      double ta = now_s();
      mg.brute = 0;
      if (hz_merge(&s1, &base, &sgf, &F.ps, &mg, &t1) != 0) return 1;
      double tb = now_s();
      mg.brute = 1;
      if (hz_merge(&s2, &base, &sgf, &F.ps, &mg, &t2) != 0) return 1;
      double tc = now_s();
      int64_t bad = 0;
      for (int32_t t = 0; t < base.nt; t++)
        if (s1.label[t] != s2.label[t]) bad++;
      printf("\n   сетка кандидатов против перебора: участков %d против %d, "
             "разметка расходится на %lld треугольниках; пар %lld против %lld; "
             "%.3f с против %.3f с (×%.1f); ячеек %lld, вхождений %lld, "
             "проверено пар %lld против %lld\n",
             s1.nseg, s2.nseg, (long long)bad, (long long)t1.npair, (long long)t2.npair, tb - ta,
             tc - tb, (tc - tb) / ((tb - ta) > 0.0 ? (tb - ta) : 1e-9), (long long)t1.ngridcell,
             (long long)t1.ngrident, (long long)t1.ncand, (long long)t2.ncand);
      if (bad != 0 || s1.nseg != s2.nseg) printf("   ОТКАЗ: сетка не равна перебору\n");
      hz_seg_free(&s1);
      hz_seg_free(&s2);
      fflush(stdout);
    }

    /* СВИП ПО ЦЕЛИ. Ш7 снял две точки (708 и 258) и увидел между ними разницу
     * вчетверо по медиане; где именно ошибка становится видна — вопрос свипа, а
     * не двух точек. Критерий чисто ГЕОМЕТРИЧЕСКИЙ: члены по полю и материалу
     * убраны из кода вслед за замером Ш7 (вторая поправка). */
    printf("\n   свип по цели, критерий ГЕОМЕТРИЧЕСКИЙ (δ огрубления = %g м):\n", eps);
    const int32_t tg[9] = {850, 708, 600, 500, 400, 300, 258, 200, 150};
    for (int i = 0; i < 9; i++) {
      hz_mergecfg mc;
      memset(&mc, 0, sizeof mc);
      mc.delta = eps;
      mc.target = tg[i];
      mc.use_geom = 1;
      mc.use_overlap = 1;
      hz_pseglist sgc;
      hz_mergestat ms;
      if (hz_merge(&sgc, &base, &sgf, &F.ps, &mc, &ms) != 0) return 1;
      scene C;
      if (scene_from(&C, &base, &sgc, base.lo, base.hi, 8) != 0) return 1;
      if (solve(&C, &d, h, 1e-4) != 0) return 1;
      imgstat sc;
      char im[64];
      snprintf(im, sizeof im, "img/sh7_geom_t%d.ppm", tg[i]);
      if (render(&C, &cam, &sc, im, NULL) != 0) return 1;
      printf("   цель %4d: %7d %9.3e %9.3e %8.2f%% %8.2f%%  dmax %9.3e (%.2f δ)  "
             "пар %lld, слито %lld, отсев: геом %lld, перекр %lld\n",
             tg[i], C.ps.np, sc.p50, sc.p99, sc.frac1, sc.frac10, ms.dmax_worst,
             ms.dmax_worst / delta, (long long)ms.npair, (long long)ms.nmerged,
             (long long)ms.nrej_geom, (long long)ms.nrej_overlap);
      fflush(stdout);
      scene_free(&C);
      hz_seg_free(&sgc);
    }

    /* НЕГАТИВНЫЙ КОНТРОЛЬ: слияние по СЛУЧАЙНОЙ метрике. Ошибка обязана стать
     * O(1), и ловят её `p99` с `dmax`, а не медиана: при разрушенной геометрии
     * медиана даже УЛУЧШАЕТСЯ, потому что эталон считает точный свет в той же
     * уехавшей точке (узор К13/К40/К94, Ш7). */
    printf("\n   НЕГАТИВНЫЙ КОНТРОЛЬ (случайная метрика):\n");
    for (int i = 0; i < 2; i++) {
      int32_t t = (int32_t[]){450, 250}[i];
      hz_mergecfg mc;
      memset(&mc, 0, sizeof mc);
      mc.delta = eps;
      mc.target = t;
      mc.use_geom = 1;
      mc.random = 1;
      hz_pseglist sgc;
      hz_mergestat ms;
      if (hz_merge(&sgc, &base, &sgf, &F.ps, &mc, &ms) != 0) return 1;
      scene C;
      if (scene_from(&C, &base, &sgc, base.lo, base.hi, 8) != 0) return 1;
      if (solve(&C, &d, h, 1e-4) != 0) return 1;
      imgstat sc;
      if (render(&C, &cam, &sc, NULL, NULL) != 0) return 1;
      printf("   СЛУЧАЙНО цель %4d: %7d %9.3e %9.3e %8.2f%% %8.2f%%  dmax %9.3e (%.2f δ)\n", t,
             C.ps.np, sc.p50, sc.p99, sc.frac1, sc.frac10, ms.dmax_worst, ms.dmax_worst / delta);
      fflush(stdout);
      scene_free(&C);
      hz_seg_free(&sgc);
    }
    scene_free(&F);
    hz_seg_free(&sgf);
    tr3_dirs_free(&d);
    hz_obj_free(&base);
    return 0;
  }

  /* --- ГЛАВНАЯ ТАБЛИЦА: бюджет разрезов против равного числа полигонов --- */
  const int32_t budg[4] = {0, 256, 1024, 4096};
  for (int b = 0; b < 4; b++) {
    hz_cutcfg cfg = {eps, tol, budg[b]};
    imgstat sa, sb;
    int32_t npa = 0;
    int nc = 0;
    char lab[64], im[64];
    snprintf(lab, sizeof lab, "δ=%g + РАЗРЕЗЫ, бюджет %d", delta, budg[b]);
    snprintf(im, sizeof im, "sh6_cut_b%d", budg[b]);
    if (one(&base, &cam, &d, delta, h, 8, budg[b] > 0, &cfg, 0.0, lab, im, &sa, &npa, &nc, 0.0) !=
        0)
      return 1;
    if (budg[b] == 0) continue;
    double dm = match_delta(&base, delta, npa, 8);
    snprintf(lab, sizeof lab, "  то же числом полигонов: δ=%.4f", dm);
    snprintf(im, sizeof im, "sh6_nocut_b%d", budg[b]);
    if (one(&base, &cam, &d, dm, h, 8, 0, &cfg, 0.0, lab, im, &sb, NULL, NULL, 0.0) != 0) return 1;
    printf("        ИТОГ бюджета %d: p99 с разрезами %.3e против %.3e без — %s\n", budg[b], sa.p99,
           sb.p99, (sa.p99 < sb.p99) ? "РАЗРЕЗЫ ЛУЧШЕ" : "разрезы НЕ лучше");
  }

  /* --- КОНТРОЛЬ ДИАГНОЗА: те же деньги, но на ПРОСТРАНСТВЕННОЕ дробление ---
   * Силуэтные разрезы кладут линию РАЗРЫВА; сетка просто делает полигон МЕЛЬЧЕ.
   * Что из этого лечит ошибку — то и было её причиной. */
  printf("\n   контроль диагноза: дробление ПРОСТРАНСТВЕННОЙ сеткой (свет и силуэты не знает)\n");
  for (int i = 0; i < 4; i++) {
    double L = (double[]){2.0, 1.0, 0.5, 0.25}[i];
    hz_cutcfg cfg = {eps, tol, 0};
    imgstat sg;
    char lab[64], im[64];
    snprintf(lab, sizeof lab, "сетка %.2f м", L);
    snprintf(im, sizeof im, "sh6_grid%d", i);
    if (one(&base, &cam, &d, delta, h, 8, 2, &cfg, 0.0, lab, im, &sg, NULL, NULL, L) != 0) return 1;
  }

  /* --- КОНТРОЛЬ LOD: дробление ПО ДАЛЬНОСТИ ОТ ГЛАЗА, L = εR --- */
  printf("\n   контроль LOD: шаг дробления растёт с дальностью (L = eps·R)\n");
  for (int i = 0; i < 3; i++) {
    double L = (double[]){0.20, 0.10, 0.05}[i];
    hz_cutcfg cfg = {eps, tol, 0};
    imgstat sg;
    char lab[64], im[64];
    snprintf(lab, sizeof lab, "LOD eps = %.2f", L);
    snprintf(im, sizeof im, "sh6_lod%d", i);
    if (one(&base, &cam, &d, delta, h, 8, 3, &cfg, 0.0, lab, im, &sg, NULL, NULL, L) != 0) return 1;
  }

  /* --- НЕГАТИВНЫЙ КОНТРОЛЬ: разрезы от ДРУГОГО положения источника --- */
  printf("\n   негативный контроль: разрезы построены под источник, сдвинутый на 1.7 м\n");
  {
    hz_cutcfg cfg = {eps, tol, 1024};
    imgstat sc;
    if (one(&base, &cam, &d, delta, h, 8, 1, &cfg, 1.7, "НЕГ.КОНТРОЛЬ: чужой свет", "sh6_negctl",
            &sc, NULL, NULL, 0.0) != 0)
      return 1;
  }

  /* --- СВИП ПО ЧИСЛУ ИСТОЧНИКОВ (обязательная часть Ш6) --- */
  printf("\n   свип по числу источников при бюджете 1024\n");
  for (int is = 0; is < 4; is++) {
    int nsrc = (int[]){1, 2, 4, 8}[is];
    hz_cutcfg cfg = {eps, tol, 1024};
    imgstat sa, sb;
    int32_t npa = 0;
    char lab[64], im[64];
    snprintf(lab, sizeof lab, "источников %d, РАЗРЕЗЫ", nsrc);
    snprintf(im, sizeof im, "sh6_src%d_cut", nsrc);
    if (one(&base, &cam, &d, delta, h, nsrc, 1, &cfg, 0.0, lab, im, &sa, &npa, NULL, 0.0) != 0)
      return 1;
    double dm = match_delta(&base, delta, npa, nsrc);
    snprintf(lab, sizeof lab, "источников %d, без разрезов", nsrc);
    snprintf(im, sizeof im, "sh6_src%d_nocut", nsrc);
    if (one(&base, &cam, &d, dm, h, nsrc, 0, &cfg, 0.0, lab, im, &sb, NULL, NULL, 0.0) != 0)
      return 1;
  }

  tr3_dirs_free(&d);
  hz_obj_free(&base);
  return 0;
}
