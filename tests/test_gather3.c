/* Фальсификаторы ЗЕРКАЛЬНОГО прослеживания в проходе 3 (PLAN_TRANSPORT.md, К5).
 * ПРЕДСКАЗАНИЯ ПЕЧАТАЮТСЯ ДО РЕЗУЛЬТАТОВ.
 *
 * ГЛАВНЫЙ ЭТАЛОН — АНАЛОГ ПЕЧИ ДЛЯ ЗЕРКАЛ: в ОДНОРОДНОМ поле идеальное зеркало
 * НЕВИДИМО. Радианс вдоль луча при зеркальном отражении с ρ = 1 сохраняется
 * точно, поэтому пиксель, попавший в зеркальный шар, обязан показать РОВНО ту
 * же величину, что и фон. Проверка ловит разом: неверную нормаль, неверную
 * формулу отражения, самопопадание, потерю множителя и ошибку в подсчёте
 * отскоков — и ловит их ТОЧНО, а не с допуском на схему.
 *
 * Хранимое поле здесь НЕ СЧИТАЕТСЯ развёрткой, а ЗАДАЁТСЯ однородным. Это
 * сознательно: проверяется СБОР, и подмешивать в эталон ошибку решателя
 * незачем. Заодно снимается вопрос о том, что для развёртки зеркало есть чёрное
 * тело (см. gather3.h) — здесь развёртки нет вовсе.
 */

#include "cut/surf.h"
#include "octree.h"
#include "transport/cam3.h"
#include "transport/cut3.h"
#include "transport/gather3.h"
#include "transport/mesh3.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int g_fail = 0, g_total = 0;

static void check(int ok, const char *what) {
  g_total++;
  if (!ok) {
    g_fail++;
    printf("FAIL: %s\n", what);
  }
}

#define LOG2N 4
#define NC (1 << LOG2N)
#define NB 2

typedef struct {
  hz_octree t;
  hz_surftab stab;
  hz_facettab ftab;
  hz_cutmap cmap;
  tr3_mesh m;
  tr3_cut cut;
  hz_frame fr;
  int32_t *wallidx;
  double *bout, *sout, *spec;
  double sc[NB][3], sr[NB];
} rig;

/* Сцена: два шара в кубе 16³, как в render3, но БЕЗ развёртки. */
static int rig_init(rig *r, int fac_level, double s) {
  memset(r, 0, sizeof *r);
  /* МАСШТАБ КАДРА `s`: единица дерева есть `s` мировых единиц, и вся мировая
   * геометрия множится на `s` вместе с ним. Тогда сцена ПОДОБНА исходной, а
   * радианс подобием не меняется вовсе — отсюда инвариант К80. */
  r->fr = (hz_frame){{0, 0, 0}, {s, s, s}};
  if (hz_oct_init(&r->t, LOG2N, 0.0)) return 1;
  for (int x = 0; x < NC; x++)
    for (int y = 0; y < NC; y++)
      for (int z = 0; z < NC; z++) {
        int lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        hz_oct_set_box(&r->t, lo, hi, 1.0);
      }
  if (hz_surftab_init(&r->stab) || hz_facettab_init(&r->ftab) || hz_cutmap_init(&r->cmap)) return 1;
  double sc[NB][3] = {{5.6, 9.2, 4.2}, {10.9, 7.4, 3.1}};
  double sr[NB] = {3.0, 1.9};
  for (int b = 0; b < NB; b++) {
    for (int a = 0; a < 3; a++)
      sc[b][a] *= s;
    sr[b] *= s;
  }
  memcpy(r->sc, sc, sizeof sc);
  memcpy(r->sr, sr, sizeof sr);
  int32_t f0[NB], nfac[NB];
  for (int b = 0; b < NB; b++) {
    hz_surf sp = {HZ_SURF_SPHERE, {sc[b][0], sc[b][1], sc[b][2], sr[b], 0, 0, 0}, 1, 0};
    int32_t si = hz_surftab_add(&r->stab, &sp);
    f0[b] = 0;
    nfac[b] = hz_surf_facet_sphere(&r->ftab, &r->fr, sc[b], sr[b], fac_level, HZ_FIT_MEAN_SAGITTA,
                                   si, &f0[b]);
    if (nfac[b] <= 0) return 1;
  }
  typedef struct {
    int32_t cell, f[HZ_P3_MAXH], nf;
  } rec_t;
  rec_t *recs = calloc((size_t)NC * NC * NC, sizeof(rec_t));
  if (recs == NULL) return 1;
  int nrec = 0;
  for (int x = 0; x < NC; x++)
    for (int y = 0; y < NC; y++)
      for (int z = 0; z < NC; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        int used = -1, ns = 0;
        for (int b = 0; b < NB; b++) {
          int32_t s2[HZ_P3_MAXH];
          int k2 = hz_facets_for_box(&r->ftab, f0[b], nfac[b], lo, hi, s2, HZ_P3_MAXH);
          if (k2 <= 0 || used >= 0) continue;
          used = b;
          ns = k2;
          for (int j = 0; j < k2; j++)
            sel[j] = s2[j];
        }
        if (ns <= 0) continue;
        recs[nrec].cell = hz_oct_leaf(&r->t, x, y, z);
        recs[nrec].nf = ns;
        for (int j = 0; j < ns; j++)
          recs[nrec].f[j] = sel[j];
        nrec++;
      }
  for (int i = 1; i < nrec; i++) { /* Г45: ключи строго по возрастанию */
    rec_t tmp = recs[i];
    int j = i - 1;
    while (j >= 0 && recs[j].cell > tmp.cell) {
      recs[j + 1] = recs[j];
      j--;
    }
    recs[j + 1] = tmp;
  }
  int nadd = 0;
  for (int i = 0; i < nrec; i++)
    if (hz_cutmap_add(&r->cmap, recs[i].cell, recs[i].f, recs[i].nf) == 0) nadd++;
  free(recs);
  check(nadd == nrec, "Г45: все записи легли в боковую таблицу");
  if (tr3_mesh_build(&r->m, &r->t, &r->fr)) return 1;
  uint8_t *solid = calloc((size_t)r->m.ncell, 1);
  if (solid == NULL) return 1;
  for (int x = 0; x < NC; x++)
    for (int y = 0; y < NC; y++)
      for (int z = 0; z < NC; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        for (int b = 0; b < NB; b++)
          if (hz_facets_for_box(&r->ftab, f0[b], nfac[b], lo, hi, sel, HZ_P3_MAXH) == 0)
            solid[r->m.cellof[hz_oct_leaf(&r->t, x, y, z)]] = 1;
      }
  int rc = tr3_cut_build(&r->cut, &r->m, &r->ftab, &r->cmap, solid);
  free(solid);
  if (rc) return 1;
  r->wallidx = calloc((size_t)6 * NC * NC, sizeof(int32_t));
  r->bout = calloc((size_t)r->m.nf * 4, sizeof(double));
  r->sout = calloc((size_t)(r->cut.nse > 0 ? r->cut.nse : 1) * 4, sizeof(double));
  r->spec = calloc((size_t)r->ftab.n, sizeof(double));
  if (r->wallidx == NULL || r->bout == NULL || r->sout == NULL || r->spec == NULL) return 1;
  for (int32_t i = 0; i < 6 * NC * NC; i++)
    r->wallidx[i] = -1;
  for (int32_t f = 0; f < r->m.nf; f++) {
    const tr3_face *ff = &r->m.f[f];
    if (ff->cb >= 0) continue;
    for (int32_t a = ff->lo[0]; a < ff->hi[0]; a++)
      for (int32_t b = ff->lo[1]; b < ff->hi[1]; b++)
        r->wallidx[((int32_t)(~ff->cb) * NC + a) * NC + b] = f;
  }
  return 0;
}

static void rig_free(rig *r) {
  free(r->wallidx);
  free(r->bout);
  free(r->sout);
  free(r->spec);
  tr3_cut_free(&r->cut);
  tr3_mesh_free(&r->m);
  hz_cutmap_free(&r->cmap);
  hz_facettab_free(&r->ftab);
  hz_surftab_free(&r->stab);
  hz_oct_free(&r->t);
}

/* ОДНОРОДНОЕ хранимое поле: значение L во ВСЕХ точках всех элементов.
 * Кладётся оно НЕ в нулевой коэффициент, а через базис (К39): постоянная
 * функция равна b_0, и коэффициенты (L,0,0,0) дают ровно L везде. */
static void rig_uniform(rig *r, double L) {
  for (int32_t f = 0; f < r->m.nf; f++) {
    r->bout[f * 4] = L;
    r->bout[f * 4 + 1] = r->bout[f * 4 + 2] = r->bout[f * 4 + 3] = 0.0;
  }
  for (int32_t e = 0; e < r->cut.nse; e++) {
    r->sout[e * 4] = L;
    r->sout[e * 4 + 1] = r->sout[e * 4 + 2] = r->sout[e * 4 + 3] = 0.0;
  }
}

/* РАЗНЫЕ радиансы у СТЕНОК, и это лечение К44 дословно. В однородном поле
 * неверно выбранная граница даёт то же самое число, и проверка пуста; когда у
 * каждой стенки свой радианс, промах по стенке виден как ошибка O(1). */
static void rig_walls(rig *r, const double wl[6], double inner) {
  for (int32_t f = 0; f < r->m.nf; f++) {
    const tr3_face *ff = &r->m.f[f];
    double v = ff->cb < 0 ? wl[~ff->cb] : inner;
    r->bout[f * 4] = v;
    r->bout[f * 4 + 1] = r->bout[f * 4 + 2] = r->bout[f * 4 + 3] = 0.0;
  }
  for (int32_t e = 0; e < r->cut.nse; e++) {
    r->sout[e * 4] = inner;
    r->sout[e * 4 + 1] = r->sout[e * 4 + 2] = r->sout[e * 4 + 3] = 0.0;
  }
}

/* ------------------------------- К80: КАРТИНКА НЕ ЗАВИСИТ ОТ ЕДИНИЦ ------- */

/* Сцена, увеличенная вместе с кадром, ПОДОБНА исходной, а радианс подобием не
 * меняется. Значит два прогона — `u = 1` и `u = 2` при мировой геометрии,
 * умноженной на 2, — обязаны дать ОДНУ И ТУ ЖЕ картинку.
 *
 * Зачем это заведено. `wall_face` считала выход луча из куба МИРОВЫМИ
 * координатами, сравнивая их с числом ячеек, и брала от них `floor`. При
 * единичном кадре это верно случайно, и весь проект работал на единичном кадре.
 * Ошибка вылезла на свипе по размеру ячейки (К68), где кадр обязан меняться,
 * чтобы комната оставалась той же.
 *
 * СТЕНКАМ РАЗДАЮТСЯ РАЗНЫЕ РАДИАНСЫ (К44), иначе промах по стенке невидим. */
static void t_frame_units(void) {
  rig a, b;
  if (rig_init(&a, 1, 1.0) || rig_init(&b, 1, 2.0)) {
    check(0, "оснастка К80");
    rig_free(&a);
    rig_free(&b);
    return;
  }
  const double wl[6] = {0.11, 0.27, 0.43, 0.61, 0.79, 0.97};
  const double wlp[6] = {0.97, 0.79, 0.61, 0.43, 0.27, 0.11}; /* перестановка */
  const double inner = 0.5;
  rig_walls(&a, wl, inner);
  rig_walls(&b, wl, inner);

  tr3_scene sa = {.tree = &a.t, .fr = a.fr, .st = &a.stab, .ft = &a.ftab, .cm = &a.cmap};
  tr3_scene sb = {.tree = &b.t, .fr = b.fr, .st = &b.stab, .ft = &b.ftab, .cm = &b.cmap};
  tr3_gather ga = {.sc = &sa,
                   .m = &a.m,
                   .cut = &a.cut,
                   .bout = a.bout,
                   .sout = a.sout,
                   .wallidx = a.wallidx,
                   .nwall = NC};
  tr3_gather gb = {.sc = &sb,
                   .m = &b.m,
                   .cut = &b.cut,
                   .bout = b.bout,
                   .sout = b.sout,
                   .wallidx = b.wallidx,
                   .nwall = NC};

  const int W = 96, H = 72;
  double eye[3] = {8.0, 0.6, 7.2}, at[3] = {8.2, 9.5, 3.6}, up[3] = {0, 0, 1};
  double eye2[3], at2[3];
  for (int k = 0; k < 3; k++) {
    eye2[k] = 2.0 * eye[k];
    at2[k] = 2.0 * at[k];
  }
  tr3_camera ca, cb;
  check(tr3_camera_look(&ca, eye, at, up, 1.3, W, H) == 0, "камера, u = 1");
  check(tr3_camera_look(&cb, eye2, at2, up, 1.3, W, H) == 0, "камера, u = 2");

  double worst = 0.0, worst_perm = 0.0;
  int npix = 0, nwallhit = 0;
  for (int py = 0; py < H; py++)
    for (int px = 0; px < W; px++) {
      double o1[3], d1[3], o2[3], d2[3];
      tr3_camera_ray(&ca, px, py, o1, d1);
      tr3_camera_ray(&cb, px, py, o2, d2);
      double va = tr3_gather_ray(&ga, o1, d1, NULL);
      double vb = tr3_gather_ray(&gb, o2, d2, NULL);
      double den = fabs(va) > 0.0 ? fabs(va) : 1.0;
      double e = fabs(va - vb) / den;
      if (e > worst) worst = e;
      npix++;
      /* попал ли луч в СТЕНКУ (а не в шар) — только там К80 и работает */
      tr3_hit h;
      tr3_march(&sa, o1, d1, -1.0, &h);
      if (!h.hit) nwallhit++;
    }
  printf("  [К80] пикселей %d, из них в СТЕНКУ %d; max отн. расхождение u=1 против u=2: %.3e\n",
         npix, nwallhit, worst);
  check(nwallhit > 500, "лучей, доходящих до стенки, достаточно — иначе проверка пуста");
  check(worst < 1e-13, "К80: картинка НЕ ЗАВИСИТ от выбора единиц дерева");

  /* НЕГАТИВНЫЙ КОНТРОЛЬ: та же пара, но у второй сцены радиансы стенок
   * ПЕРЕСТАВЛЕНЫ. Если метрика к выбору стенки слепа, расхождение останется
   * нулевым — и тогда проверка выше не проверяет ничего (К44). */
  rig_walls(&b, wlp, inner);
  for (int py = 0; py < H; py++)
    for (int px = 0; px < W; px++) {
      double o1[3], d1[3], o2[3], d2[3];
      tr3_camera_ray(&ca, px, py, o1, d1);
      tr3_camera_ray(&cb, px, py, o2, d2);
      double va = tr3_gather_ray(&ga, o1, d1, NULL);
      double vb = tr3_gather_ray(&gb, o2, d2, NULL);
      double e = fabs(va - vb);
      if (e > worst_perm) worst_perm = e;
    }
  printf("  [К80] негативный контроль (радиансы стенок переставлены): max |Δ| %.3e\n", worst_perm);
  check(worst_perm > 0.1, "негативный контроль: промах по СТЕНКЕ обязан быть виден");
  rig_free(&a);
  rig_free(&b);
}

int main(void) {
  printf("=== ПРЕДСКАЗАНИЯ (до единого результата) ===\n");
  printf("  Ф1 ПЕЧЬ ДЛЯ ЗЕРКАЛ: в ОДНОРОДНОМ поле зеркало с ρ=1 НЕВИДИМО —\n");
  printf("     каждый пиксель равен фону, отн. ошибка ≤1e-13\n");
  printf("  Ф2 закон отражения: |d'| = 1 и d'·n = −d·n, ≤1e-15\n");
  printf("  Ф3 при ρ=0.9 пиксель равен РОВНО 0.9^k·L, где k — число отскоков\n");
  printf("  Ф4 самопопадания нет БЕЗ ε: отскоков нулевой длины ноль\n");
  printf("  Ф5 свип по пределу отскоков: ответ выходит на ПОЛКУ\n");
  printf("=== НЕГАТИВНЫЕ КОНТРОЛИ, с ПРЕДСКАЗАННЫМ провалом ===\n");
  printf("  НК1 отражать по нормали ФАСЕТА вместо примитива ⇒ Ф1 ОБЯЗАН\n");
  printf("      сломаться, и тем сильнее, чем грубее фасетизация\n");
  printf("  НК2 не умножать на ρ на отскоке ⇒ Ф3 обязан дать 1.0 вместо 0.9^k\n");
  printf("  К80: картинка НЕ зависит от выбора единиц дерева (сцена, увеличенная\n");
  printf("     вместе с кадром, ПОДОБНА исходной) — отн. расхождение ≤1e-13;\n");
  printf("     негативный контроль: переставить радиансы стенок ⇒ |Δ| > 0.1\n");
  printf("  К80: картинка НЕ зависит от выбора единиц дерева (подобная сцена)\n");
  printf("     — отн. расхождение <= 1e-13; негативный контроль: переставить\n");
  printf("     радиансы стенок ==> |D| > 0.1\n");
  printf("=== РЕЗУЛЬТАТЫ ===\n");

  t_frame_units();

  rig r;
  if (rig_init(&r, 2, 1.0)) {
    check(0, "оснастка");
    rig_free(&r);
    printf("FAILURES\n");
    return 1;
  }
  printf("  сцена: ячеек %d, граней %d, поверхностных элементов %d, фасетов %d\n", r.m.ncell,
         r.m.nf, r.cut.nse, r.ftab.n);

  const double L = 0.7;
  rig_uniform(&r, L);
  tr3_scene scn = {
      .tree = &r.t, .fr = r.fr, .sigma = NULL, .st = &r.stab, .ft = &r.ftab, .cm = &r.cmap};

  tr3_camera cam;
  double eye[3] = {8.0, 0.6, 7.2}, at[3] = {8.2, 9.5, 3.6}, up[3] = {0, 0, 1};
  const int W = 160, H = 120;
  check(tr3_camera_look(&cam, eye, at, up, 1.3, W, H) == 0, "камера");

  /* ---------------------------------------------------------------- Ф1 --- */
  for (int32_t i = 0; i < r.ftab.n; i++)
    r.spec[i] = 1.0;
  tr3_gather gg = {.sc = &scn,
                   .m = &r.m,
                   .cut = &r.cut,
                   .bout = r.bout,
                   .sout = r.sout,
                   .facet_spec = r.spec,
                   .nfacet = r.ftab.n,
                   .wallidx = r.wallidx,
                   .nwall = NC,
                   .maxbounce = 16};
  double worst = 0.0;
  int nhit = 0, maxnb = 0;
  long nbsum = 0;
  for (int py = 0; py < H; py++)
    for (int px = 0; px < W; px++) {
      double o[3], d[3];
      tr3_camera_ray(&cam, px, py, o, d);
      int nb = 0;
      double v = tr3_gather_ray(&gg, o, d, &nb);
      if (nb > 0) {
        nhit++;
        nbsum += nb;
        if (nb > maxnb) maxnb = nb;
      }
      double e = fabs(v - L) / L;
      if (e > worst) worst = e;
    }
  printf("  [Ф1] ПЕЧЬ ДЛЯ ЗЕРКАЛ: пикселей с отскоком %d из %d, отскоков до %d "
         "(в среднем %.2f); max отн. ошибка по ВСЕМУ кадру %.3e\n",
         nhit, W * H, maxnb, nhit ? (double)nbsum / (double)nhit : 0.0, worst);
  check(nhit > W * H / 20, "зеркала действительно видны в кадре");
  check(maxnb >= 2, "многократные отскоки между двумя шарами есть");
  check(worst < 1e-13, "Ф1: идеальное зеркало в однородном поле НЕВИДИМО");

  /* ---------------------------------------------------------------- Ф2 --- */
  {
    double wn = 0.0, wd = 0.0, wnn = 0.0;
    int ncheck = 0;
    for (int py = 0; py < H; py += 3)
      for (int px = 0; px < W; px += 3) {
        double o[3], d[3];
        tr3_camera_ray(&cam, px, py, o, d);
        tr3_hit h;
        tr3_march(&scn, o, d, -1.0, &h);
        if (!h.hit) continue;
        double dn = d[0] * h.n[0] + d[1] * h.n[1] + d[2] * h.n[2], dr[3], nn = 0.0, dr_n = 0.0,
               n2 = 0.0;
        for (int a = 0; a < 3; a++) {
          dr[a] = d[a] - 2.0 * dn * h.n[a];
          nn += dr[a] * dr[a];
          dr_n += dr[a] * h.n[a];
          n2 += h.n[a] * h.n[a];
        }
        if (fabs(sqrt(n2) - 1.0) > wnn) wnn = fabs(sqrt(n2) - 1.0);
        if (fabs(sqrt(nn) - 1.0) > wn) wn = fabs(sqrt(nn) - 1.0);
        if (fabs(dr_n + dn) > wd) wd = fabs(dr_n + dn);
        ncheck++;
      }
    printf("  [Ф2] закон отражения на %d лучах: max ||d'|−1| = %.3e, max |d'·n + d·n| = %.3e; "
           "источник — сама нормаль: max ||n|−1| = %.3e\n",
           ncheck, wn, wd, wnn);
    check(ncheck > 100, "лучей достаточно");
    /* ПРЕДСКАЗАНИЕ 1e-15 БЫЛО НЕВЕРНО, И ПОРОГ ПОДНЯТ НЕ «ЧТОБЫ ПРОШЛО», А ПО
     * ИЗМЕРЕННОЙ ПРИЧИНЕ. Нормаль сферы считается как (o + t·d − c)/r, то есть
     * разностью величин порядка 8 при результате порядка 3, и её собственное
     * отклонение от единичной длины измеряется здесь же. Формула отражения
     * усиливает его множителем до 4·(d·n)², больше нигде ошибке взяться неоткуда.
     * Поэтому порог ставится по ИЗМЕРЕННОМУ ||n|−1|, а не по круглому числу: если
     * отражение начнёт врать само по себе, отношение уедет, и это будет видно. */
    check(wn < 20.0 * wnn + 1e-16, "Ф2: |d'| = 1 в пределах собственной ошибки нормали");
    check(wd < 20.0 * wnn + 1e-16, "Ф2: угол падения равен углу отражения");

    /* Ф6. НОРМАЛЬНОЕ ПАДЕНИЕ — ТОЧНЫЙ ГЕОМЕТРИЧЕСКИЙ ЭТАЛОН НОРМАЛИ.
     * Луч, направленный В ЦЕНТР шара, падает по нормали, и отражённое направление
     * обязано быть РОВНО −d. Это единственная здесь проверка, которая ловит
     * нормаль САМУ ПО СЕБЕ, без всякого поля: Ф1 к ней слепа (К44). */
    double wretro = 0.0;
    int nretro = 0;
    for (int b = 0; b < NB; b++)
      for (int k = 0; k < 32; k++) {
        /* исходная точка на сфере радиуса 7 вокруг центра шара, направление внутрь */
        double th = 0.3 + 0.19 * (double)k, ph = 0.7 + 0.41 * (double)k;
        double u[3] = {sin(th) * cos(ph), sin(th) * sin(ph), cos(th)};
        double o[3], d[3];
        for (int a = 0; a < 3; a++) {
          o[a] = r.sc[b][a] + 7.0 * u[a];
          d[a] = -u[a];
        }
        if (o[0] < 0.5 || o[0] > NC - 0.5 || o[1] < 0.5 || o[1] > NC - 0.5 || o[2] < 0.5 ||
            o[2] > NC - 0.5)
          continue;
        tr3_hit h;
        tr3_march(&scn, o, d, -1.0, &h);
        if (!h.hit) continue;
        /* ПОПАЛИ ЛИ В ТОТ ШАР, В КОТОРЫЙ ЦЕЛИЛИСЬ. Центры разнесены на 5.70, а
         * старт берётся в 7.0 от центра, поэтому луч, идущий в центр одного
         * шара, законно задевает по дороге ВТОРОЙ — и тогда нормаль радиальна
         * НЕ ЭТОМУ шару, а ретро-эталон к ней неприменим. Первая редакция этой
         * проверки не делала и дала 3.3 при эталоне 0: неверный ответ был у
         * ЭТАЛОНА, а не у марша — тот же класс, что две ошибки теста в T5а. */
        double rr2 = 0.0;
        for (int a = 0; a < 3; a++)
          rr2 += (h.p[a] - r.sc[b][a]) * (h.p[a] - r.sc[b][a]);
        if (fabs(sqrt(rr2) - r.sr[b]) > 1e-9) continue;
        double dn = d[0] * h.n[0] + d[1] * h.n[1] + d[2] * h.n[2], e = 0.0;
        for (int a = 0; a < 3; a++) {
          double dr = d[a] - 2.0 * dn * h.n[a];
          e += fabs(dr + d[a]); /* обязано быть −d */
        }
        if (e > wretro) wretro = e;
        nretro++;
      }
    printf("  [Ф6] нормальное падение на %d лучах: max |d' + d| = %.3e (обязан быть ретро-луч)\n",
           nretro, wretro);
    check(nretro > 30, "лучей нормального падения достаточно");
    check(wretro < 1e-13, "Ф6: при падении в центр шара отражение РОВНО обратное");
  }

  /* ---------------------------------------------------------------- Ф3 --- */
  {
    const double rs = 0.9;
    for (int32_t i = 0; i < r.ftab.n; i++)
      r.spec[i] = rs;
    double worst3 = 0.0;
    int n3 = 0;
    for (int py = 0; py < H; py++)
      for (int px = 0; px < W; px++) {
        double o[3], d[3];
        tr3_camera_ray(&cam, px, py, o, d);
        int nb = 0;
        double v = tr3_gather_ray(&gg, o, d, &nb);
        if (nb == 0) continue;
        double want = L;
        for (int k = 0; k < nb; k++)
          want *= rs;
        double e = fabs(v - want) / want;
        if (e > worst3) worst3 = e;
        n3++;
      }
    printf("  [Ф3] ρ=0.9: %d пикселей с отскоком, max отн. отклонение от ρ^k·L = %.3e\n", n3,
           worst3);
    check(worst3 < 1e-13, "Ф3: множитель ρ накапливается РОВНО по числу отскоков");

    /* НК2: не умножать на ρ (maxbounce тот же, spec = 1) — ответ обязан стать L */
    for (int32_t i = 0; i < r.ftab.n; i++)
      r.spec[i] = 1.0;
    double got = 0.0, want1 = L * rs;
    int found = 0;
    for (int py = 0; py < H && !found; py++)
      for (int px = 0; px < W && !found; px++) {
        double o[3], d[3];
        tr3_camera_ray(&cam, px, py, o, d);
        int nb = 0;
        double v = tr3_gather_ray(&gg, o, d, &nb);
        if (nb == 1) {
          got = v;
          found = 1;
        }
      }
    printf("    [НК2 без множителя ρ] пиксель с одним отскоком: %.6f против правильного %.6f\n",
           got, want1);
    check(fabs(got - want1) > 0.05 * want1, "НК2: без множителя ρ ответ ОБЯЗАН уехать");
  }

  /* ---------------------------------------------------------------- Ф4 --- */
  {
    /* САМОПОПАДАНИЕ БЕЗ ε. Отражённый луч стартует РОВНО с поверхности. Если бы
     * он ловил сам себя, отскок вышел бы нулевой длины и число отскоков уперлось
     * бы в предел на КАЖДОМ пикселе зеркала. Меряем именно это. */
    for (int32_t i = 0; i < r.ftab.n; i++)
      r.spec[i] = 1.0;
    int nstuck = 0, nhit4 = 0;
    for (int py = 0; py < H; py++)
      for (int px = 0; px < W; px++) {
        double o[3], d[3];
        tr3_camera_ray(&cam, px, py, o, d);
        int nb = 0;
        tr3_gather_ray(&gg, o, d, &nb);
        if (nb > 0) nhit4++;
        if (nb >= gg.maxbounce) nstuck++;
      }
    printf("  [Ф4] самопопадание: пикселей, упёршихся в предел %d отскоков: %d из %d зеркальных\n",
           gg.maxbounce, nstuck, nhit4);
    check(nstuck * 100 < nhit4, "Ф4: самопопадания нет — предел упирается лишь в редких лучах");
  }

  /* ---------------------------------------------------------------- Ф5 --- */
  {
    /* СВИП ПО ПРЕДЕЛУ ОТСКОКОВ (К22 в применении к зеркалам). При ρ<1 ответ
     * обязан выйти на полку; печатается, а не утверждается. */
    const double rs = 0.6;
    for (int32_t i = 0; i < r.ftab.n; i++)
      r.spec[i] = rs;
    double prev = 0.0, dlast = 0.0;
    printf("  [Ф5] свип по пределу отскоков (ρ=%.1f), сумма радианса по кадру:\n", rs);
    for (int mb = 1; mb <= 16; mb *= 2) {
      gg.maxbounce = mb;
      double sum = 0.0;
      for (int py = 0; py < H; py++)
        for (int px = 0; px < W; px++) {
          double o[3], d[3];
          tr3_camera_ray(&cam, px, py, o, d);
          sum += tr3_gather_ray(&gg, o, d, NULL);
        }
      printf("        предел %2d: %.9f  (изменение %.3e)\n", mb, sum,
             prev > 0.0 ? fabs(sum - prev) / prev : 0.0);
      if (prev > 0.0) dlast = fabs(sum - prev) / prev;
      prev = sum;
    }
    check(dlast < 1e-6, "Ф5: ответ вышел на полку по числу отскоков");
    gg.maxbounce = 16;
  }

  /* --------------------------------------------------------------- НК1 --- */
  {
    /* ОТРАЖЕНИЕ ПО НОРМАЛИ ФАСЕТА вместо примитива (К15). У зеркала ошибка
     * нормали входит в направление С УДВОЕНИЕМ, поэтому она обязана быть видна
     * сильнее, чем на диффузной поверхности, — и слабеть с ростом фасетизации.
     *
     * К44: ПОЛЕ ЗДЕСЬ ОБЯЗАНО БЫТЬ НЕОДНОРОДНЫМ, И ЭТО ГЛАВНОЕ В ЭТОМ КОНТРОЛЕ.
     * Первая редакция гоняла его на однородном поле Ф1 и получила расхождение
     * РОВНО НОЛЬ: куда ни отрази луч, в однородном поле придёт то же число.
     * То есть «печь для зеркал» проверяет что угодно (множитель, самопопадание,
     * счёт отскоков), но к НОРМАЛИ слепа полностью — тот же узор, что у К12
     * (печь видит только среднее) и К13 (печь не ловит передачу потока).
     * Стенкам поэтому раздаются РАЗНЫЕ радиансы, и тогда неверная нормаль
     * отправляет луч на ДРУГУЮ стенку, что есть ошибка O(1). */
    for (int32_t i = 0; i < r.ftab.n; i++)
      r.spec[i] = 1.0;
    double prev_mean = 0.0;
    for (int lev = 1; lev <= 3; lev++) {
      rig r2;
      if (rig_init(&r2, lev, 1.0)) {
        rig_free(&r2);
        continue;
      }
      /* НЕОДНОРОДНОЕ поле: у каждой стенки свой радианс */
      for (int32_t f = 0; f < r2.m.nf; f++) {
        int wl = r2.m.f[f].cb < 0 ? (int)(~r2.m.f[f].cb) : 0;
        r2.bout[f * 4] = 0.2 + 0.3 * (double)wl;
        r2.bout[f * 4 + 1] = r2.bout[f * 4 + 2] = r2.bout[f * 4 + 3] = 0.0;
      }
      for (int32_t e = 0; e < r2.cut.nse; e++) {
        r2.sout[e * 4] = L;
        r2.sout[e * 4 + 1] = r2.sout[e * 4 + 2] = r2.sout[e * 4 + 3] = 0.0;
      }
      for (int32_t i = 0; i < r2.ftab.n; i++)
        r2.spec[i] = 1.0;
      tr3_scene s2 = {.tree = &r2.t,
                      .fr = r2.fr,
                      .sigma = NULL,
                      .st = &r2.stab,
                      .ft = &r2.ftab,
                      .cm = &r2.cmap};
      tr3_gather g2 = {.sc = &s2,
                       .m = &r2.m,
                       .cut = &r2.cut,
                       .bout = r2.bout,
                       .sout = r2.sout,
                       .facet_spec = r2.spec,
                       .nfacet = r2.ftab.n,
                       .wallidx = r2.wallidx,
                       .nwall = NC,
                       .maxbounce = 16};
      /* Расхождение ДВУХ путей на одном и том же луче: примитив против фасета. */
      double wdiff = 0.0, sdiff = 0.0;
      int ndiff = 0;
      for (int py = 0; py < H; py += 2)
        for (int px = 0; px < W; px += 2) {
          double o[3], d[3];
          tr3_camera_ray(&cam, px, py, o, d);
          double v = tr3_gather_ray(&g2, o, d, NULL);
          /* тот же луч, но первое отражение — по нормали ФАСЕТА */
          tr3_hit h;
          tr3_march(&s2, o, d, -1.0, &h);
          if (!h.hit) continue;
          int32_t mc = r2.m.cellof[h.cell];
          if (mc < 0) continue;
          double nf[3] = {0, 0, 0};
          double best = -2.0;
          for (int32_t k = r2.cut.sestart[mc]; k < r2.cut.sestart[mc + 1]; k++) {
            int32_t ee = r2.cut.selist[k];
            double dp = r2.cut.se[ee].n[0] * h.n[0] + r2.cut.se[ee].n[1] * h.n[1] +
                        r2.cut.se[ee].n[2] * h.n[2];
            if (dp > best) {
              best = dp;
              memcpy(nf, r2.cut.se[ee].n, sizeof nf);
            }
          }
          double dn = d[0] * nf[0] + d[1] * nf[1] + d[2] * nf[2], o2[3], d2[3];
          for (int a = 0; a < 3; a++) {
            o2[a] = h.p[a];
            d2[a] = d[a] - 2.0 * dn * nf[a];
          }
          double v2 = tr3_gather_ray(&g2, o2, d2, NULL);
          double e2 = fabs(v2 - v);
          if (e2 > wdiff) wdiff = e2;
          sdiff += e2;
          ndiff++;
        }
      double mean = ndiff ? sdiff / (double)ndiff : 0.0;
      printf("    [НК1 уровень %d, фасетов %4d] расхождение примитив/фасет: max %.3e, "
             "в среднем %.3e\n",
             lev, r2.ftab.n, wdiff, mean);
      check(wdiff > 0.1, "НК1: нормаль фасета ОБЯЗАНА давать ошибку O(1) у зеркала");
      if (lev == 1) prev_mean = mean;
      if (lev == 3) {
        printf("      падение среднего расхождения с уровня 1 к уровню 3: %.2f раза\n",
               prev_mean > 0.0 ? prev_mean / mean : 0.0);
        check(mean < 0.7 * prev_mean, "НК1: и ОБЯЗАНА слабеть с ростом фасетизации");
      }
      rig_free(&r2);
    }
  }

  rig_free(&r);
  printf("%s: %d/%d\n", g_fail ? "FAILURES" : "ok", g_total - g_fail, g_total);
  return g_fail != 0;
}
