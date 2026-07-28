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
#include "pcut.h"
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

#define NVIS_REF 6 /* проб видимости на источник в ЭТАЛОНЕ (6×6 = 36) */
/* Проб видимости в ПРОВЕРЯЕМОМ прогоне. Параметр, а не константа, и это
 * вынужденно: при 2×2 = 4 пробах полутень квантуется на пять ступеней, то есть
 * площадной источник §2 ведёт себя как четыре точечных. Тогда выигрыш разрезов
 * мог бы оказаться выигрышем над СОБСТВЕННОЙ ступенькой выборки, а не над
 * физикой. Свип по этому числу — обязательная часть приёмки Ш6. */
static int NVIS_RUN = 2;
/* Кадр МЕТРИКИ — 512², а не 1024² из §2. Эталон стоит (пиксели × источники ×
 * 36 лучей), и на 1024² один прогон занимает минуты, а их в свипе двадцать.
 * Для КАРТИНКИ НА ГЛАЗ это не годится, и картинка пишется отдельно в полном
 * разрешении §2; для ХВОСТА распределения 262 тысячи точек достаточно. */
#define IMGW 512
#define IMGH 512

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

/* Решение: прямой свет один раз (он не меняется), развёртка — до застоя. */
static int solve(scene *S, const tr3_dirs *d, double h, double tol) {
  hz_pstats st;
  memset(&st, 0, sizeof st);
  double t0 = now_s();
  hz_psweep_zero(&S->t);
  hz_direct_add(&S->t, &S->g, S->src, S->nsrc, S->srcLe, NVIS_RUN, h, &st);
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
    if (hz_psweep_gather(&S->t, d, h, 0, 0.0, HZ_LAYOUT_RUNS, &st) != 0) {
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
    double x[3];
    for (int a = 0; a < 3; a++)
      x[a] = o[a] + t * d[a];

    const double *ce = S->t.E + (size_t)k * 3;
    const double *cd = S->Edir + (size_t)k * 3;
    double Etot = ce[0] + ce[1] * u + ce[2] * v;
    double Edf = cd[0] + cd[1] * u + cd[2] * v;
    double Eind = Etot - Edf;
    if (Eind < 0.0) Eind = 0.0;
    int cls = 0;
    double Ede = direct_exact(S, k, x, &cls);
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
  char ppm[96], pfm[96];
  if (img != NULL) {
    snprintf(ppm, sizeof ppm, "img/%s.ppm", img);
    snprintf(pfm, sizeof pfm, "img/%s.pfm", img);
  }
  if (render(&S, cam, out, img ? ppm : NULL, img ? pfm : NULL) != 0) return 1;
  printf("   %-30s %7d %7d %5d %9.3e %9.3e %8.2f%% %8.2f%% %9.3e %9.3e %7.2f %7.2f\n", label,
         S.ps.np, docut ? S.ncut : 0, S.nbounce, out->p50, out->p99, out->frac1, out->frac10,
         out->p50e, out->p99e, S.t_direct + S.t_solve, out->t_frame);
  printf("        по классам точки  p99: свет %.3e (%d), ПОЛУТЕНЬ %.3e (%d), тень %.3e (%d)"
         "   клиньев-кандидатов %d, ПОТЕРЯНО фрагментов %lld\n",
         out->p99_lit, out->n_lit, out->p99_pen, out->n_pen, out->p99_shd, out->n_shd, ws.n,
         (long long)S.nover);
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

  hz_objmesh base;
  if (hz_obj_load(&base, HZ_CFG_HALL_OBJ, HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет %s\n", HZ_CFG_HALL_OBJ);
    return 1;
  }
  tr3_camera cam;
  double eye[3] = HZ_CFG_HALL_EYE, at[3] = HZ_CFG_HALL_AT, up[3] = HZ_CFG_UP;
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
    for (int i = 0; i < 4; i++) {
      double L = (double[]){1.60, 1.20, 0.80, 0.60}[i];
      char lab[64], im[64];
      snprintf(lab, sizeof lab, "ПЛОСКОСТЯМИ, шаг %.2f м", L);
      snprintf(im, sizeof im, "sh6_grid_%d", i);
      if (one(&base, &cam, &d, delta, h, 8, 4, &cfg0, 0.0, lab, im, &s, NULL, NULL, L) != 0)
        return 1;
    }
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
