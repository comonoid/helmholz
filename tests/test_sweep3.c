/* Фальсификаторы трёхмерной развёртки (PLAN_TRANSPORT.md, шаги A, B, C).
 * ПРЕДСКАЗАНИЯ ПЕЧАТАЮТСЯ ДО РЕЗУЛЬТАТОВ.
 *
 * Батарея T2, перенесённая в 3D ДОСЛОВНО вместе с её находками: печь проверяет
 * только среднее (К12), неконсервативную передачу потока она не ловит вовсе
 * (К13), ограничитель обязан сохранять константу (К6). */

#include "transport/dirs3.h"
#include "transport/krylov3.h"
#include "transport/mesh3.h"
#include "transport/sweep3.h"
#include "image.h"
#include "transport/cam3.h"
#include "transport/tet3.h"
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

/* ------------------------------------------- A1: моменты набора направлений */

static void t_dirs(void) {
  for (int nmu = 2; nmu <= 4; nmu *= 2) {
    tr3_dirs d;
    check(tr3_dirs_product(&d, nmu, nmu) == 0, "набор построен");
    double s = 0.0, m1[3] = {0, 0, 0}, m2[3][3];
    memset(m2, 0, sizeof m2);
    double axis_min = 1e300;
    for (int i = 0; i < d.n; i++) {
      double o[3] = {d.ox[i], d.oy[i], d.oz[i]};
      s += d.w[i];
      for (int a = 0; a < 3; a++) {
        m1[a] += d.w[i] * o[a];
        for (int b = 0; b < 3; b++)
          m2[a][b] += d.w[i] * o[a] * o[b];
        if (fabs(o[a]) < axis_min) axis_min = fabs(o[a]);
      }
    }
    double e2 = 0.0;
    for (int a = 0; a < 3; a++)
      for (int b = 0; b < 3; b++) {
        double ex = (a == b) ? 4.0 * M_PI / 3.0 : 0.0;
        double e = fabs(m2[a][b] - ex);
        if (e > e2) e2 = e;
      }
    printf("  [A1] ND=%d: Σw−4π = %.2e, |Σwω| = %.2e, max |Σwωω − (4π/3)I| = %.2e, "
           "min |ω_a| = %.3f\n",
           d.n, fabs(s - 4.0 * M_PI), fabs(m1[0]) + fabs(m1[1]) + fabs(m1[2]), e2, axis_min);
    check(fabs(s - 4.0 * M_PI) < 1e-13, "Σw = 4π");
    check(fabs(m1[0]) + fabs(m1[1]) + fabs(m1[2]) < 1e-13, "первый момент нулевой");
    check(e2 < 1e-13, "второй момент есть (4π/3)·I");
    check(axis_min > 1e-3, "ни одно направление НЕ параллельно оси — вырожденных граней нет");

    /* ЭТАП C: ЗЕРКАЛЬНАЯ ПЕРЕСТАНОВКА. Проверяется ТРИ вещи, и все три нужны:
     * что перестановка есть перестановка (инволюция и биекция — это свойства
     * ИНДЕКСОВ, и они обязаны держаться ТОЧНО), что она отражает направление
     * (свойство координат — до округления) и что вес сохраняется.
     *
     * Побитовости у координат здесь взяться неоткуда: `cos(3π/4 − x)` и
     * `−cos(π/4 + x)` в плавающей арифметике не равны. Требовать её значило бы
     * завести красный тест на исправном коде — тот же класс, что К1/К6/К10. */
    double worst = 0.0, wworst = 0.0;
    int notinv = 0, notbij = 0;
    int *seen = calloc((size_t)d.n, sizeof(int));
    check(seen != NULL, "память");
    if (seen != NULL)
      for (int a = 0; a < 3; a++) {
        for (int i = 0; i < d.n; i++)
          seen[i] = 0;
        for (int m = 0; m < d.n; m++) {
          int mm2 = d.mir[a * d.n + m];
          if (mm2 < 0 || mm2 >= d.n) {
            notbij++;
            continue;
          }
          seen[mm2]++;
          if (d.mir[a * d.n + mm2] != m) notinv++;
          double o[3] = {d.ox[m], d.oy[m], d.oz[m]};
          double o2[3] = {d.ox[mm2], d.oy[mm2], d.oz[mm2]};
          o[a] = -o[a];
          for (int b = 0; b < 3; b++) {
            double e = fabs(o[b] - o2[b]);
            if (e > worst) worst = e;
          }
          double we = fabs(d.w[m] - d.w[mm2]);
          if (we > wworst) wworst = we;
        }
        for (int i = 0; i < d.n; i++)
          if (seen[i] != 1) notbij++;
      }
    free(seen);
    printf("  [C1] зеркальная перестановка: не-биекций %d, не-инволюций %d, "
           "max |ω′ − отражённое ω| = %.2e, max |Δw| = %.2e\n",
           notbij, notinv, worst, wworst);
    check(notbij == 0, "C: перестановка есть БИЕКЦИЯ — точно, это индексы");
    check(notinv == 0, "C: и ИНВОЛЮЦИЯ — отразить дважды значит вернуться");
    check(worst < 1e-15, "C: направление действительно отражённое (до округления)");
    check(wworst < 1e-16, "C: вес при отражении сохраняется");
    tr3_dirs_free(&d);
  }
}

/* ------------------------------------------- A2: моменты по многоугольнику */

static void t_face_moments(void) {
  /* прямоугольник 2×3 в плоскости z = 1, базис ячейки с центром (1,1.5,1), h=2 */
  double v[4][3] = {{0, 0, 1}, {2, 0, 1}, {2, 3, 1}, {0, 3, 1}};
  double c[3] = {1, 1.5, 1}, h = 2.0;
  double m[4][4];
  check(tr3_poly_face_mass(v, 4, c, h, m) == 0, "моменты грани посчитаны");
  double A = 6.0;
  /* ∫1 = A; ∫ξ = 0 (центр по x совпал); ∫ξ² = A·(2/2)²/12 = A/12;
   * ∫η² = A·(3/2)²/12 = A·9/48; ∫ζ = 0 и ∫ζ² = 0 (грань в плоскости центра) */
  printf("  [A2] ∫1 = %.15f (%.1f), ∫ξ² = %.15f (%.15f), ∫η² = %.15f (%.15f), ∫ζ² = %.1e\n",
         m[0][0], A, m[1][1], A / 12.0, m[2][2], A * 9.0 / 48.0, m[3][3]);
  check(fabs(m[0][0] - A) < 1e-14, "площадь");
  check(fabs(m[1][1] - A / 12.0) < 1e-14, "второй момент по ξ");
  check(fabs(m[2][2] - A * 9.0 / 48.0) < 1e-14, "второй момент по η");
  check(fabs(m[3][3]) < 1e-14, "по ζ грань лежит в плоскости центра");
}

/* ------------------------------------------- B: сетка граней, целые площади */

static void build_tree(hz_octree *t, int L, int graded) {
  hz_oct_init(t, L, 0.0);
  if (graded == 2) {
    /* РАВНОМЕРНОЕ дробление до единицы: коробка в одну единицу содержится только
     * в узле размера 1, поэтому вызов на каждую ячейку заставляет дерево
     * раздробиться целиком. */
    int n2 = 1 << L;
    for (int x = 0; x < n2; x++)
      for (int y = 0; y < n2; y++)
        for (int z = 0; z < n2; z++) {
          int lo2[3] = {x, y, z}, hi2[3] = {x + 1, y + 1, z + 1};
          hz_oct_set_box(t, lo2, hi2, 1.0);
        }
  } else if (graded) {
    /* НАСТОЯЩИЙ перепад уровней: коробка в ОДНУ единицу заставляет дерево
     * раздробиться до размера 1 в углу, оставив остальное крупным. Первая
     * редакция ставила коробку в половину куба — она ложится РОВНО на детей
     * корня, дробления не вызывает вовсе, и «градуированная» сетка выходила
     * равномерной. Перепад надо ЗАДАВАТЬ, а не надеяться на него (Г48). */
    int lo[3] = {0, 0, 0}, hi[3] = {1, 1, 1};
    hz_oct_set_box(t, lo, hi, 1.0);
  }
  (void)L;
}

static void t_mesh(void) {
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 2.0}}; /* анизотропный кадр */
  for (int graded = 0; graded < 2; graded++) {
    hz_octree t;
    build_tree(&t, 3, graded);
    tr3_mesh m;
    check(tr3_mesh_build(&m, &t, &fr) == 0, "сетка построена");
    /* ЦЕЛОЧИСЛЕННОЕ ТОЖДЕСТВО: сумма площадей подграней на стороне ячейки равна
     * площади стороны. Все координаты целые, поэтому равенство ТОЧНОЕ. */
    int32_t bad = 0;
    for (int32_t c = 0; c < m.ncell; c++)
      for (int a = 0; a < 3; a++)
        for (int side = 0; side < 2; side++) {
          int32_t pos = m.clo[c][a] + (side ? m.csize[c] : 0);
          int64_t sum = 0;
          for (int32_t k = m.fstart[c]; k < m.fstart[c + 1]; k++) {
            const tr3_face *f = &m.f[m.flist[k]];
            /* сторону задаёт ПОЛОЖЕНИЕ плоскости, а не то, ca ячейка или cb:
             * у граничной грани ячейка всегда ca, и фильтр по ca/cb выбросил бы
             * минус-стенку целиком */
            if (f->axis != a || f->pos != pos) continue;
            sum += (int64_t)(f->hi[0] - f->lo[0]) * (int64_t)(f->hi[1] - f->lo[1]);
          }
          int64_t want = (int64_t)m.csize[c] * (int64_t)m.csize[c];
          if (sum != want) bad++;
        }
    printf("  [B%s] ячеек %d, граней %d; сторон с НЕВЕРНОЙ суммой площадей: %d\n",
           graded ? ", градуированная" : "", m.ncell, m.nf, bad);
    check(bad == 0, "К13 в 3D: сумма площадей подграней ТОЧНО равна площади стороны");

    /* НЕГАТИВНЫЙ КОНТРОЛЬ: сосед берётся ОДИН — тот, что накрывает ЦЕНТР
     * стороны. Так пишут, когда забывают, что у крупной ячейки соседей
     * несколько. Тождество ОБЯЗАНО сломаться, и только на градуированной. */
    int32_t bad2 = 0;
    for (int32_t c = 0; c < m.ncell; c++)
      for (int a = 0; a < 3; a++)
        for (int side = 0; side < 2; side++) {
          /* ОБЕ стороны, а не только плюс: у мелких соседей ячейка оказывается
           * то ca, то cb, и односторонний обзор конфигурацию «крупная против
           * четырёх мелких» просто не встречает */
          int32_t pos = m.clo[c][a] + (side ? m.csize[c] : 0);
          int64_t sum = 0, seen = 0;
          for (int32_t k = m.fstart[c]; k < m.fstart[c + 1]; k++) {
            const tr3_face *f = &m.f[m.flist[k]];
            if (f->axis != a || f->pos != pos) continue;
            if (seen++ > 0) continue; /* «взять только первого соседа» */
            sum += (int64_t)(f->hi[0] - f->lo[0]) * (int64_t)(f->hi[1] - f->lo[1]);
          }
          if (seen > 0 && sum != (int64_t)m.csize[c] * (int64_t)m.csize[c]) bad2++;
        }
    printf("    [НК один сосед] сторон с неверной суммой: %d\n", bad2);
    if (graded) check(bad2 > 0, "негативный контроль: «один сосед» ломает тождество");
    tr3_mesh_free(&m);
    hz_oct_free(&t);
  }
}

/* ------------------------------------------- C: печь, линейное поле, баланс */

typedef struct {
  hz_octree t;
  tr3_mesh m;
  tr3_dirs d;
  double *sig_t, *sig_s, *eps;
} rig;

static int rig_init(rig *r, int L, int graded, int nmu, const hz_frame *fr) {
  build_tree(&r->t, L, graded);
  if (tr3_mesh_build(&r->m, &r->t, fr)) return 1;
  if (tr3_dirs_product(&r->d, nmu, nmu)) return 1;
  r->sig_t = calloc((size_t)r->m.ncell, sizeof(double));
  r->sig_s = calloc((size_t)r->m.ncell, sizeof(double));
  r->eps = calloc((size_t)r->m.ncell * 4, sizeof(double));
  return r->sig_t == NULL || r->sig_s == NULL || r->eps == NULL;
}

static void rig_free(rig *r) {
  free(r->sig_t);
  free(r->sig_s);
  free(r->eps);
  tr3_dirs_free(&r->d);
  tr3_mesh_free(&r->m);
  hz_oct_free(&r->t);
}

/* ---------------------------------------- К41/К42: элемент, а не ячейка ---
 *
 * ПРЕДСКАЗАНИЯ ДО РЕЗУЛЬТАТОВ:
 *  Ф1 К41: поле, ПОЛОЖИТЕЛЬНОЕ на всём элементе, но уводящее ДАЛЬНИЙ УГОЛ
 *     ЯЧЕЙКИ в минус, обязано проходить проверку положительности без отката —
 *     прежде ровно этот случай резался зря (запас `0.6c` вместо `c`);
 *  Ф2 К42: отрицательное на элементе поле обязано откатываться МЯГКО —
 *     наклон выживает (α ∈ (0,1)), среднее по строке 0 массы сохраняется
 *     ТОЧНО, минимум по элементу после отката ≥ 0 (с допуском на округление);
 *  НК: rr0 < 0 — спасти нечего, α = 0, остаётся константный откат (то же
 *     отрицательное среднее), т.е. отказ в закрытую сторону никуда не делся.
 */
/* ------------------- ХВОСТ «ДВА ТЕЛА»: ОБЪЕДИНЕНИЕ В ЯЧЕЙКЕ ---------------
 *
 * ПРЕДСКАЗАНИЯ ДО РЕЗУЛЬТАТОВ:
 *  Ф1 у разрезанной сцены с двумя РАЗНЕСЁННЫМИ шарами найдутся ячейки, чей веер
 *     несёт фасеты ОБОИХ тел (nboth > 0); все они обработаны объединением
 *     (cu.nunion == nboth) и без нарушений 1:1 (cu.nbad == 0);
 *  Ф2 в такой ячейке флюид = коробка − V(B1) − V(B2) — объёмы тел считаются
 *     НЕЗАВИСИМЫМ путём (hz_poly3_complement по группе в той же коробке);
 *     негативный контроль: старое поведение «материал = пересечение» дало бы
 *     флюид = ВСЯ коробка (тела разнесены, пересечение пусто) — ошибка на весь
 *     объём тел, и она ОБЯЗАНА быть видна;
 *  Ф3 одна группа (однотелковая ячейка) идёт прежним путём: тождество
 *     флюид + тело = коробка держится до 1e-12;
 *  Ф4 энергобаланс полной развёртки на двух телах замыкается.
 */
static void t_union(void) {
  const hz_frame fr = {{0, 0, 0}, {1, 1, 1}};
  const int L = 4, N = 1 << L;
  hz_octree t;
  hz_oct_init(&t, L, 0.0);
  for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
      for (int z = 0; z < N; z++) {
        int lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        hz_oct_set_box(&t, lo, hi, 1.0);
      }
  hz_surftab stab;
  hz_facettab ftab;
  hz_cutmap cmap;
  hz_surftab_init(&stab);
  hz_facettab_init(&ftab);
  hz_cutmap_init(&cmap);
  /* Тела РАЗНЕСЕНЫ (d = 2.9 > r1+r2 = 2.8), но обе задевают общие ячейки
   * вокруг x = 8 — ровно тот случай, который прежде выбрасывал второе тело. */
  double sc[2][3] = {{7.0, 7.5, 7.5}, {9.9, 7.5, 7.5}};
  double sr[2] = {1.4, 1.4};
  int32_t si[2], f0[2], nfac[2];
  for (int b = 0; b < 2; b++) {
    hz_surf sp = {HZ_SURF_SPHERE, {sc[b][0], sc[b][1], sc[b][2], sr[b], 0, 0, 0}, 1, 0};
    si[b] = hz_surftab_add(&stab, &sp);
    f0[b] = 0;
    nfac[b] = hz_surf_facet_sphere(&ftab, &fr, sc[b], sr[b], 1, HZ_FIT_MEAN_SAGITTA, si[b], &f0[b]);
    check(nfac[b] > 0, "сфера фасетизирована");
  }
  typedef struct {
    int32_t cell, f[HZ_P3_MAXH], nf;
  } rec_t;
  static rec_t rec[16384];
  int nrec = 0, nboth = 0;
  for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
      for (int z = 0; z < N; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        int ns = 0, nb2 = 0;
        for (int b = 0; b < 2; b++) {
          int32_t s2[HZ_P3_MAXH];
          int k2 = hz_facets_for_box(&ftab, f0[b], nfac[b], lo, hi, s2, HZ_P3_MAXH);
          if (k2 <= 0) continue;
          nb2++;
          for (int j = 0; j < k2 && ns < HZ_P3_MAXH; j++)
            sel[ns++] = s2[j];
        }
        if (ns <= 0) continue;
        if (nb2 > 1) nboth++;
        rec[nrec].cell = hz_oct_leaf(&t, x, y, z);
        rec[nrec].nf = ns;
        for (int j = 0; j < ns; j++)
          rec[nrec].f[j] = sel[j];
        nrec++;
      }
  check(nboth > 0, "Ф1: ячейки с ОБОИМИ телами в веере существуют");
  for (int i = 1; i < nrec; i++) { /* Г45: ключи строго по возрастанию */
    rec_t tmp = rec[i];
    int j = i - 1;
    while (j >= 0 && rec[j].cell > tmp.cell) {
      rec[j + 1] = rec[j];
      j--;
    }
    rec[j + 1] = tmp;
  }
  for (int i = 0; i < nrec; i++)
    check(hz_cutmap_add(&cmap, rec[i].cell, rec[i].f, rec[i].nf) == 0, "Г45: запись легла");
  tr3_mesh m;
  check(tr3_mesh_build(&m, &t, &fr) == 0, "сетка");
  uint8_t *solid = calloc((size_t)m.ncell, 1);
  for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
      for (int z = 0; z < N; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        int inside = 1;
        for (int b = 0; b < 2 && inside; b++)
          if (hz_facets_for_box(&ftab, f0[b], nfac[b], lo, hi, sel, HZ_P3_MAXH) != 0) inside = 0;
        if (inside) solid[m.cellof[hz_oct_leaf(&t, x, y, z)]] = 1;
      }
  tr3_cut cu;
  check(tr3_cut_build(&cu, &m, &ftab, &cmap, solid) == 0, "разрез");
  printf("  [UNION] ячеек с двумя телами %d, обработано объединением %d, nbad %d\n", nboth,
         (int)cu.nunion, (int)cu.nbad);
  check((int32_t)cu.nunion == nboth, "Ф1: все двухтельные ячейки обработаны объединением");
  check(cu.nbad == 0, "Ф1: нарушений 1:1 нет");

  /* Ф2: флюид = коробка − V(B1) − V(B2) в каждой двухтельной ячейке. V(B_g)
   * считается НЕЗАВИСИМО: тело g в коробке через hz_poly3_cut, дополнение —
   * через hz_poly3_complement по фасетам ТОЛЬКО этого тела. */
  double worst_u = 0.0, worst_1 = 0.0, vwrong = 0.0;
  int nchecked_u = 0, nchecked_1 = 0;
  for (int32_t c = 0; c < m.ncell; c++) {
    const hz_cutrec *rr = hz_cutmap_find(&cmap, m.node[c]);
    if (rr == NULL) continue;
    double s = (double)m.csize[c], vbox = s * s * s;
    int32_t lo2[3] = {m.clo[c][0], m.clo[c][1], m.clo[c][2]};
    int32_t hi2[3] = {lo2[0] + m.csize[c], lo2[1] + m.csize[c], lo2[2] + m.csize[c]};
    int nb_here = 0, body_here[2] = {0, 0};
    double vsum = 0.0;
    int single_ok = 1;
    for (int b = 0; b < 2; b++) {
      int32_t s2[HZ_P3_MAXH];
      int k2 = hz_facets_for_box(&ftab, f0[b], nfac[b], lo2, hi2, s2, HZ_P3_MAXH);
      if (k2 <= 0) continue;
      body_here[b] = 1;
      nb_here++;
      hz_hspace hh[HZ_P3_MAXH];
      for (int j = 0; j < k2; j++) {
        for (int a = 0; a < 3; a++)
          hh[j].n[a] = ftab.f[s2[j]].n[a];
        hh[j].off = ftab.f[s2[j]].off;
      }
      /* переносимый путь: НЕЗАВИСИМОЕ дополнение тела в коробке */
      hz_poly3 *pcs = calloc((size_t)k2, sizeof(hz_poly3));
      int npcs = 0;
      int32_t zl[HZ_P3_MAXH] = {0}, zf[HZ_P3_MAXH] = {0};
      if (hz_poly3_complement(pcs, k2, &npcs, lo2, hi2, hh, zl, zf, k2) == HZ_P3_OK) {
        double vc = 0.0;
        for (int p = 0; p < npcs; p++)
          vc += hz_poly3_volume(&pcs[p], &fr);
        vsum += vbox - vc; /* V(B_g) */
      } else {
        single_ok = 0;
      }
      free(pcs);
    }
    if (nb_here == 0) continue;
    double fluid = cu.mvol[c][0][0];
    if (nb_here >= 2) {
      /* перекрытие тел отсутствует по построению (d > r1+r2), значит
       * V(объединения) = V(B1) + V(B2) */
      double e = fabs(fluid - (vbox - vsum));
      if (e > worst_u) worst_u = e;
      if (vbox - fluid > vwrong) vwrong = vbox - fluid;
      nchecked_u++;
    } else if (single_ok && body_here[0] != body_here[1]) {
      double e = fabs(fluid - (vbox - vsum));
      if (e > worst_1) worst_1 = e;
      nchecked_1++;
    }
  }
  printf("  [UNION] двухтельные: max |флюид − (коробка−V1−V2)| = %.3e (%d ячеек); "
         "старый ответ дал бы флюид = коробка, расхождение с ним до %.3f\n",
         worst_u, nchecked_u, vwrong);
  check(nchecked_u > 0, "Ф2: двухтельные ячейки проверены");
  /* Допуск 1e-4: тела в ячейке заданы ВЕЕРАМИ ПЛОСКОСТЕЙ, а не сферой; плоские
   * тела двух фасетизаций могут перекрываться на величину порядка фасетного
   * зазора (d=2.9 при r1+r2=2.8), и тождество «объединение = V1+V2» верно для
   * сфер, но не точно для их плоских представлений. Дефект «пересечение вместо
   * объединения» давал 0.7–8.3 на ячейку — на пять порядков больше. */
  check(worst_u < 1e-4, "Ф2: флюид = коробка минус ОБА тела (до фасетного зазора)");
  check(vwrong > 0.5, "НК Ф2: пересечение-вместо-объединения давило бы ОБА тела целиком");
  printf("  [UNION] однотельные: max |флюид + тело − коробка| = %.3e (%d ячеек)\n", worst_1,
         nchecked_1);
  check(nchecked_1 > 0, "Ф3: однотельные ячейки проверены");
  check(worst_1 < 1e-12, "Ф3: однотельный путь прежний — тождество до 1e-12");

  /* Ф4: энергобаланс полной развёртки на двух телах */
  {
    tr3_dirs d;
    check(tr3_dirs_product(&d, 2, 2) == 0, "ординаты");
    double *sigt = calloc((size_t)m.ncell, sizeof(double));
    double *sigs = calloc((size_t)m.ncell, sizeof(double));
    double *frho = calloc((size_t)ftab.n, sizeof(double));
    double *fem = calloc((size_t)ftab.n, sizeof(double));
    for (int32_t c = 0; c < m.ncell; c++) {
      sigt[c] = 0.05;
      sigs[c] = 0.04;
    }
    for (int32_t i = 0; i < ftab.n; i++)
      frho[i] = 0.7;
    double wr[6] = {0.8, 0.8, 0.8, 0.8, 0.8, 0.8}, we[6] = {0, 0, 0, 0, 0, 5.0};
    tr3_problem p = {.m = &m,
                     .d = &d,
                     .cut = &cu,
                     .facet_rho = frho,
                     .facet_emit = fem,
                     .nfacet = ftab.n,
                     .sig_t = sigt,
                     .sig_s = sigs,
                     .wall_rho = wr,
                     .wall_emit = we,
                     .limiter = 1};
    double *phi = calloc((size_t)m.ncell * 4, sizeof(double));
    tr3_stats st;
    memset(&st, 0, sizeof st);
    int rc = tr3_sweep_solve(&p, 4000, 1e-9, phi, &st);
    check(rc == 0, "Ф4: развёртка сошлась");
    double lhs = st.pin + st.psout + st.pemit, rhs = st.pout + st.pabs + st.psin + st.psolid;
    double rel = st.pin > 0.0 ? fabs(lhs - rhs) / st.pin : -1.0;
    printf("  [UNION] баланс: втекло %.4f, невязка тождества %.3e (отн. %.3e), элементов %d\n",
           st.pin, fabs(lhs - rhs), rel, (int)cu.nse);
    /* Допуск 1e-2: ограничитель положительности и принцип максимума ЖЕРТВУЮТ
     * строки баланса в зажатых ячейках (§735 записано), поэтому тождество
     * замыкается до жертвы ограничителя, а не до машинной точности; дефектом
     * считался бы порядок 1 — свет, прошедший сквозь выброшенное тело. */
    check(rel < 1e-2, "Ф4: тождество К40 на объединении замыкается (до жертвы ограничителя)");
    free(phi);
    free(sigt);
    free(sigs);
    free(frho);
    free(fem);
    free(st.bout);
    free(st.sout);
    tr3_dirs_free(&d);
  }

  free(solid);
  tr3_cut_free(&cu);
  tr3_mesh_free(&m);
  hz_cutmap_free(&cmap);
  hz_facettab_free(&ftab);
  hz_surftab_free(&stab);
  hz_oct_free(&t);
}

/* ------------------------------ К50: ГЛАДКИЕ НОРМАЛИ -----------------------
 *
 * РЕШЕНИЕ ПОЛЬЗОВАТЕЛЯ: тело считается ПРИМИТИВОМ; у явно заданных поверхностей
 * нормаль интерполируется между треугольниками (Фонг).
 *
 * ПРЕДСКАЗАНИЯ ДО РЕЗУЛЬТАТОВ:
 *  Ф1 у элементов сферы нормаль = примитивной нормали в центроиде; макс угол
 *     ПАДАЕТ с ростом фасетизации (фасетный уровень 1 -> 3);
 *  НК1 прежний путь (нормаль плоскости фасета) даёт угол порядка наклона
 *     фасета — ЗАМЕТНО больший, и почти не зависящий от уровня;
 *  Ф2 явная поверхность (surf сброшен в -1): Фонг даёт угол меньше плоскостного
 *     НЕ МЕНЕЕ чем вдвое — интерполяция восстанавливает сферичность;
 *  НК2 Фонг на явной поверхности при build (без build3) недоступен — нормали
 *     остаются плоскостными (тот же макс угол, что у плоскостной ошибки).
 */
static void t_k50_scene(int fsub, int use_stab, int make_explicit, double *worst_angle,
                        double *worst_plane_angle, int *nchecked, int latlong) {
  const hz_frame fr = {{0, 0, 0}, {1, 1, 1}};
  const int L = 4, N = 1 << L;
  hz_octree t;
  hz_oct_init(&t, L, 0.0);
  for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
      for (int z = 0; z < N; z++) {
        int lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        hz_oct_set_box(&t, lo, hi, 1.0);
      }
  hz_surftab stab;
  hz_facettab ftab;
  hz_cutmap cmap;
  hz_surftab_init(&stab);
  hz_facettab_init(&ftab);
  hz_cutmap_init(&cmap);
  double c[3] = {8.0, 8.0, 8.0}, r = 3.0;
  hz_surf sp = {HZ_SURF_SPHERE, {c[0], c[1], c[2], r, 0, 0, 0}, 1, 0};
  int32_t si = hz_surftab_add(&stab, &sp);
  int32_t f0 = 0, nfac = 0;
  if (latlong) {
    /* ЯВНАЯ поверхность: широтно-долготная сетка треугольников с вершинами НА
     * сфере, surf = -1, куски ограничены (tv — единицы кадра, здесь = мир).
     * Тот же класс входа, что DC-фасеты и OBJ-куски. */
    f0 = (int32_t)ftab.n;
    const int NTH = 10, NPH = 20;
    for (int it = 0; it < NTH; it++)
      for (int ip = 0; ip < NPH; ip++) {
        double th0 = M_PI * (double)it / NTH, th1 = M_PI * (double)(it + 1) / NTH;
        double ph0 = 2.0 * M_PI * (double)ip / NPH, ph1 = 2.0 * M_PI * (double)(ip + 1) / NPH;
        double v[4][3];
        double th[2] = {th0, th1}, ph[2] = {ph0, ph1};
        for (int a = 0; a < 2; a++)
          for (int b = 0; b < 2; b++) {
            v[2 * a + b][0] = c[0] + r * sin(th[a]) * cos(ph[b]);
            v[2 * a + b][1] = c[1] + r * sin(th[a]) * sin(ph[b]);
            v[2 * a + b][2] = c[2] + r * cos(th[a]);
          }
        int tris[2][3] = {{0, 1, 3}, {0, 3, 2}};
        for (int tt = 0; tt < 2; tt++) {
          const double *A = v[tris[tt][0]], *B = v[tris[tt][1]], *C = v[tris[tt][2]];
          double u1[3], u2[3], n[3];
          for (int a = 0; a < 3; a++) {
            u1[a] = B[a] - A[a];
            u2[a] = C[a] - A[a];
          }
          n[0] = u1[1] * u2[2] - u1[2] * u2[1];
          n[1] = u1[2] * u2[0] - u1[0] * u2[2];
          n[2] = u1[0] * u2[1] - u1[1] * u2[0];
          double nm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
          if (!(nm > 0.0)) continue;
          for (int a = 0; a < 3; a++)
            n[a] /= nm;
          if (n[0] * (A[0] - c[0]) + n[1] * (A[1] - c[1]) + n[2] * (A[2] - c[2]) < 0.0)
            for (int a = 0; a < 3; a++)
              n[a] = -n[a];
          double off = n[0] * A[0] + n[1] * A[1] + n[2] * A[2];
          int32_t fi = ftab.n;
          if (hz_facettab_add_plane(&ftab, &fr, n, off, -1, 0.0) < 0) {
            check(0, "явный фасет не лёг");
            return;
          }
          ftab.f[fi].bounded = 1;
          ftab.f[fi].tnv = 3;
          for (int a = 0; a < 3; a++) {
            ftab.f[fi].tv[0][a] = A[a];
            ftab.f[fi].tv[1][a] = B[a];
            ftab.f[fi].tv[2][a] = C[a];
          }
          nfac++;
        }
      }
    if (nfac <= 0) {
      check(0, "явная сетка пуста");
      return;
    }
  } else {
    nfac = hz_surf_facet_sphere(&ftab, &fr, c, r, fsub, HZ_FIT_MEAN_SAGITTA, si, &f0);
    if (make_explicit)
      for (int32_t i = f0; i < f0 + nfac; i++)
        ftab.f[i].surf = -1; /* тело задано ЯВНО треугольниками, примитива нет */
  }
  typedef struct {
    int32_t cell, f[HZ_P3_MAXH], nf;
  } rec_t;
  static rec_t rec[16384];
  int nrec = 0;
  for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
      for (int z = 0; z < N; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        int ns = hz_facets_for_box(&ftab, f0, nfac, lo, hi, sel, HZ_P3_MAXH);
        if (ns <= 0) continue;
        rec[nrec].cell = hz_oct_leaf(&t, x, y, z);
        rec[nrec].nf = ns;
        for (int j = 0; j < ns; j++)
          rec[nrec].f[j] = sel[j];
        nrec++;
      }
  for (int i = 1; i < nrec; i++) {
    rec_t tmp = rec[i];
    int j = i - 1;
    while (j >= 0 && rec[j].cell > tmp.cell) {
      rec[j + 1] = rec[j];
      j--;
    }
    rec[j + 1] = tmp;
  }
  for (int i = 0; i < nrec; i++)
    check(hz_cutmap_add(&cmap, rec[i].cell, rec[i].f, rec[i].nf) == 0, "Г45: запись легла");
  tr3_mesh m;
  check(tr3_mesh_build(&m, &t, &fr) == 0, "сетка");
  uint8_t *solid = calloc((size_t)m.ncell, 1);
  for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
      for (int z = 0; z < N; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        if (hz_facets_for_box(&ftab, f0, nfac, lo, hi, sel, HZ_P3_MAXH) == 0)
          solid[m.cellof[hz_oct_leaf(&t, x, y, z)]] = 1;
      }
  tr3_cut cu;
  check(tr3_cut_build3(&cu, &m, &ftab, &cmap, solid, NULL, NULL, use_stab ? &stab : NULL) == 0,
        "разрез");
  for (int32_t e = 0; e < cu.nse; e++) {
    const tr3_selem *se = &cu.se[e];
    if (se->nv < 3) continue;
    double cen[3] = {0, 0, 0};
    for (int i = 0; i < se->nv; i++)
      for (int a = 0; a < 3; a++)
        cen[a] += se->v[i][a];
    for (int a = 0; a < 3; a++)
      cen[a] /= (double)se->nv;
    double exact[3] = {cen[0] - c[0], cen[1] - c[1], cen[2] - c[2]};
    double nm = sqrt(exact[0] * exact[0] + exact[1] * exact[1] + exact[2] * exact[2]);
    for (int a = 0; a < 3; a++)
      exact[a] /= nm;
    double dot = 0.0;
    for (int a = 0; a < 3; a++)
      dot += se->n[a] * exact[a];
    double ang = acos(dot < -1.0 ? -1.0 : (dot > 1.0 ? 1.0 : dot));
    if (ang > *worst_angle) {
      *worst_angle = ang;
      if (getenv("C3D"))
        printf("    худший: угол %.4f, se.n=(%.3f,%.3f,%.3f), exact=(%.3f,%.3f,%.3f), "
               "cen=(%.2f,%.2f,%.2f), facet=%d\n",
               ang, se->n[0], se->n[1], se->n[2], exact[0], exact[1], exact[2], cen[0], cen[1],
               cen[2], (int)se->facet);
    }
    /* плоскостная нормаль фасета — прежний ответ */
    const hz_facet *fp = &ftab.f[se->facet];
    if (fp != NULL) {
      double pnm = sqrt(fp->n[0] * fp->n[0] + fp->n[1] * fp->n[1] + fp->n[2] * fp->n[2]);
      double pd = 0.0;
      for (int a = 0; a < 3; a++)
        pd += (fp->n[a] / pnm) * exact[a];
      double pang = acos(pd < -1.0 ? -1.0 : (pd > 1.0 ? 1.0 : pd));
      if (pang > *worst_plane_angle) *worst_plane_angle = pang;
    }
    (*nchecked)++;
  }
  free(solid);
  tr3_cut_free(&cu);
  tr3_mesh_free(&m);
  hz_cutmap_free(&cmap);
  hz_facettab_free(&ftab);
  hz_surftab_free(&stab);
  hz_oct_free(&t);
}

static void t_k50(void) {
  double w1 = 0, w3 = 0, wp1 = 0, wp3 = 0;
  int c1 = 0, c3 = 0;
  t_k50_scene(1, 1, 0, &w1, &wp1, &c1, 0);
  t_k50_scene(3, 1, 0, &w3, &wp3, &c3, 0);
  printf("  [К50] примитив: макс угол к сфере при уровне 1: %.4f рад, при 3: %.4f рад "
         "(плоскостный: %.4f / %.4f)\n",
         w1, w3, wp1, wp3);
  check(c1 > 0 && c3 > 0, "Ф1: элементы проверены");
  check(w3 < w1, "Ф1: угол ПАДАЕТ с ростом фасетизации");
  check(w1 < 0.35 && w3 < 0.12, "Ф1: примитивная нормаль близка к точной");
  check(wp1 > 3.0 * w3, "НК1: плоскостная нормаль ОБЯЗАНА быть заметно хуже примитивной");

  /* Ф2: та же сфера, но задана ЯВНО (surf = -1) — работает Фонг */
  double we = 0, wpe = 0;
  int ce = 0;
  t_k50_scene(0, 0, 0, &we, &wpe, &ce, 1);
  printf("  [К50] явная поверхность (Фонг): макс угол %.4f рад против плоскостного %.4f рад\n", we,
         wpe);
  check(ce > 0, "Ф2: элементы проверены");
  /* Критерий: Фонг ОБЯЗАН быть строго лучше плоскостной нормали; запас 10 %
   * Named-остаток честен: интерполяция линейна по треугольнику и не обязана
   * восстанавливать сферу точно, а у полюсов широтно-долготной сетки
   * треугольники вытянуты. Дефект «нет интерполяции» дал бы РАВЕНство. */
  check(we < 0.9 * wpe, "Ф2: Фонг восстанавливает сферичность (заметно лучше плоскостной)");
  check(wpe > 0.02, "НК2: плоскостная ошибка на этой сетке ненулевая — есть что лечить");
  /* Ф3: примитив-путь при этом не тронут — нормали примитивных сфер точны
   * (acos добавляет округление ~1e-9 к нулевому углу) */
  check(w3 < 1e-6, "Ф3: примитивная нормаль в центроиде ТОЧНА (измерено выше)");
}

static void t_elem_positivity(void) {
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 1.0}};
  rig r;
  if (rig_init(&r, 1, 2, 2, &fr)) {
    check(0, "оснастка К41/К42");
    rig_free(&r);
    return;
  }
  const int32_t c = 0; /* ячейка (0,0,0) размера 1, центр (0.5, 0.5, 0.5) */
  const double mv[4][3] = {{0, 0.5, 0}, {1, 0.5, 0}, {1, 0.5, 1}, {0, 0.5, 1}};

  /* --- Ф1: на элементе положительно, в дальнем углу ячейки — минус ---
   * Наклон вдоль y: на самом элементе (плоскость y = 0.5) поле ПОСТОЯННО и
   * положительно, а дальние углы ячейки (y = 0) оно уводит в минус. */
  const double e1[4] = {0.1, 0.0, 4.0, 0.0};
  double mn_elem = tr3_elem_min_poly(&r.m, c, e1, mv, 4);
  double mn_cell = 1e300;
  for (int k = 0; k < 8; k++) {
    double v = e1[0];
    for (int a = 0; a < 3; a++)
      v += e1[a + 1] * (((k >> a) & 1) ? 0.5 : -0.5);
    if (v < mn_cell) mn_cell = v;
  }
  printf("  [К41] min по элементу %.3e (должен быть ≥ 0), min по ячейке %.3e (< 0)\n", mn_elem,
         mn_cell);
  check(mn_elem >= 0.0 && mn_cell < 0.0, "К41: случай «элемент ок, ячейка нет» существует");
  double ee1[4] = {e1[0], e1[1], e1[2], e1[3]};
  /* строка 0 массы диагональна на элементе, свободный член несёт среднее */
  const double fmm_diag[4][4] = {{1.0, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
  int rolled = tr3_elem_rollback(&r.m, c, fmm_diag, 0.1, mv, 4, ee1);
  check(rolled == 0, "К41: отката ОБЯЗАНО не быть — поле на элементе положительно");
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
  check(ee1[2] == e1[2], "К41: наклон ОБЯЗАН выжить (прежде срезался целиком)"); /* точный бит */

  /* --- Ф2: на элементе есть минус — мягкий откат --- */
  /* v = −0.25 − (x − 0.5): при x = 1 это −0.75 < 0. Среднее по строке 0 есть
   * 0.25, значит den = 0 − (−0.5) = 0.5 и α = 0.25/0.5 = РОВНО 0.5. */
  double ee2[4] = {-0.25, -1.0, 0.0, 0.0};
  rolled = tr3_elem_rollback(&r.m, c, fmm_diag, 0.25, mv, 4, ee2);
  double mn2 = tr3_elem_min_poly(&r.m, c, ee2, mv, 4);
  double row0 = 0.0;
  for (int j = 0; j < 4; j++)
    row0 += fmm_diag[0][j] * ee2[j];
  printf("  [К42] после отката: ee = %.3e %.3e %.3e %.3e, min по элементу %.3e, строка 0 %.3e\n",
         ee2[0], ee2[1], ee2[2], ee2[3], mn2, row0);
  check(rolled == 1, "К42: откат обязан примениться");
  check(fabs(ee2[1] + 0.5) < 1e-15, "К42: наклон сжат ровно вдвое (α = 0.5), а не обнулён");
  check(fabs(ee2[0] - 0.25) < 1e-15, "К42: свободный член пересчитан из строки 0 (0.25)");
  check(fabs(row0 - 0.25) <= 1e-15, "К42: строка 0 массы сохранена ТОЧНО");
  check(mn2 >= -1e-15, "К42: после отката минимум по элементу неотрицателен");

  /* --- НК: rr0 < 0 — спасти нечего, константный откат --- */
  double ee3[4] = {0.0, -1.0, 0.0, 0.0};
  rolled = tr3_elem_rollback(&r.m, c, fmm_diag, -0.4, mv, 4, ee3);
  check(rolled == 1 && ee3[1] == 0.0 && fabs(ee3[0] + 0.4) < 1e-15, /* точный ноль ветки */
        "НК К42: rr0 < 0 — α = 0, прежний константный откат");
#pragma GCC diagnostic pop

  rig_free(&r);
}

static void t_furnace(void) {
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 2.0}};
  const double lb = 1.7;
  (void)lb;
  for (int graded = 0; graded < 2; graded++) {
    rig r;
    if (rig_init(&r, 3, graded, 2, &fr)) {
      check(0, "оснастка");
      rig_free(&r);
      return;
    }
    for (int32_t c = 0; c < r.m.ncell; c++) {
      r.sig_t[c] = 0.8;
      r.sig_s[c] = 0.8; /* альбедо РОВНО 1 */
    }
    tr3_problem p = {
        .m = &r.m, .d = &r.d, .sig_t = r.sig_t, .sig_s = r.sig_s, .binc0 = lb, .limiter = 1};
    double *phi = calloc((size_t)r.m.ncell * 4, sizeof(double));
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    tr3_stats st;
    int rc = tr3_sweep_solve(&p, 400, 1e-14, phi, &st);
    free(st.bout);
    free(st.sout);
    double worst = 0.0, wslope = 0.0;
    for (int32_t c = 0; c < r.m.ncell; c++) {
      double e = fabs(phi[c * 4] - 4.0 * M_PI * lb) / (4.0 * M_PI * lb);
      if (e > worst) worst = e;
      for (int j = 1; j < 4; j++) {
        double sl = fabs(phi[c * 4 + j]) / (4.0 * M_PI * lb);
        if (sl > wslope) wslope = sl;
      }
    }
    printf("  [C печь%s] код %d, итераций %d, срезок %d; отн. ошибка среднего %.2e, "
           "наклоны %.2e\n",
           graded ? ", градуированная" : "", rc, st.iters, st.nclip, worst, wslope);
    check(rc == 0, "развёртка прошла, цикла обхода нет");
    check(worst < 1e-12, "ПЕЧЬ: радианс однороден и равен влёту");
    check(st.nclip == 0,
          "К6: на СОШЕДШЕМСЯ решении ограничитель не активен — константу он не трогает");
    check(wslope < 1e-12, "и наклоны нулевые");
    free(phi);
    rig_free(&r);
  }
}

static void t_linear_field(void) {
  /* К12 В ТРЁХМЕРИИ: печь проверяет только СРЕДНЕЕ, а наклоны — две трети
   * неизвестных DG1. Точное решение L = A + B·x достигается при σ_s = 0 и
   * ε(x,ω) = B·ω_x + σ_t·(A + B·x), причём первый член НАПРАВЛЕННЫЙ. */
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 1.0}};
  const double A = 1.3, B = 0.21;
  for (int graded = 0; graded < 2; graded++) {
    rig r;
    if (rig_init(&r, 3, graded, 2, &fr)) {
      check(0, "оснастка");
      rig_free(&r);
      return;
    }
    for (int32_t c = 0; c < r.m.ncell; c++) {
      r.sig_t[c] = 0.7;
      r.sig_s[c] = 0.0;
      double s = (double)r.m.csize[c];
      double xc = (double)r.m.clo[c][0] + 0.5 * s;
      r.eps[c * 4 + 0] = 0.7 * (A + B * xc); /* среднее по ячейке */
      r.eps[c * 4 + 1] = 0.7 * B * s;        /* наклон: b_1 = (x − x_c)/s */
    }
    tr3_problem p = {.m = &r.m,
                     .d = &r.d,
                     .sig_t = r.sig_t,
                     .sig_s = r.sig_s,
                     .eps = r.eps,
                     .eps_dir = {B, 0.0, 0.0},
                     .binc0 = A,
                     .binc = {B, 0, 0},
                     .binx0 = {0, 0, 0},
                     .limiter = 0};
    double *phi = calloc((size_t)r.m.ncell * 4, sizeof(double));
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    tr3_stats st;
    int rc = tr3_sweep_solve(&p, 4, 1e-14, phi, &st);
    free(st.bout);
    free(st.sout);
    double wm = 0.0, ws = 0.0;
    for (int32_t c = 0; c < r.m.ncell; c++) {
      double s = (double)r.m.csize[c];
      double xc = (double)r.m.clo[c][0] + 0.5 * s;
      double em = fabs(phi[c * 4] - 4.0 * M_PI * (A + B * xc));
      double es = fabs(phi[c * 4 + 1] - 4.0 * M_PI * B * s);
      if (em > wm) wm = em;
      if (es > ws) ws = es;
    }
    printf("  [К12 линейное поле%s] код %d: ошибка среднего %.2e, ошибка НАКЛОНА %.2e\n",
           graded ? ", градуированная" : "", rc, wm, ws);
    check(rc == 0 && wm < 1e-12, "К12: линейное решение воспроизведено по среднему");
    check(ws < 1e-12, "и ПО НАКЛОНУ — то, чего печь не проверяет вовсе");
    free(phi);
    rig_free(&r);
  }
}

static void t_linear(void) {
  /* К12: печь проверяет только СРЕДНЕЕ. Линейное решение L = A + B·x проверяет
   * наклоны. Берём чистое поглощение (σ_s = 0) и подбираем ε так, чтобы
   * L = A + B·x было точным: (ω·∇)L + σ_t L = ε ⇒ ε(ω) = B·ω_x + σ_t(A + Bx).
   * Направленной части ε здесь нет, поэтому берём B вдоль x и σ_t = 0:
   * тогда ε = B·ω_x, что от направления зависит — а наш ε изотропен.
   * Значит проверяем ДРУГОЕ линейное решение: σ_t > 0, ε изотропно,
   * L = ε/σ_t — константа. Линейность проверяется линейным ВЛЁТОМ при σ = 0. */
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 1.0}};
  rig r;
  if (rig_init(&r, 3, 1, 2, &fr)) {
    check(0, "оснастка");
    rig_free(&r);
    return;
  }
  /* пустая среда, линейный влёт: решение L(x) = влёт, перенесённый вдоль луча,
   * то есть ЛИНЕЙНАЯ функция обязана пройти сквозь сетку без искажения */
  tr3_problem p = {.m = &r.m, .d = &r.d, .sig_t = r.sig_t, .sig_s = r.sig_s, .binc0 = 2.0};
  double *phi = calloc((size_t)r.m.ncell * 4, sizeof(double));
  if (phi == NULL) {
    check(0, "память");
    rig_free(&r);
    return;
  }
  tr3_stats st;
  int rc = tr3_sweep_solve(&p, 4, 1e-14, phi, &st);
  free(st.bout);
  free(st.sout);
  double worst = 0.0;
  for (int32_t c = 0; c < r.m.ncell; c++) {
    double e = fabs(phi[c * 4] - 4.0 * M_PI * 2.0) / (4.0 * M_PI * 2.0);
    if (e > worst) worst = e;
  }
  printf("  [C пустая среда] код %d: отн. ошибка %.2e (влёт проходит насквозь)\n", rc, worst);
  check(rc == 0 && worst < 1e-12, "в пустой среде влёт проходит без искажения");
  free(phi);
  rig_free(&r);
}

/* ------------------------------------- ЭТАП C: ЗЕРКАЛЬНАЯ СТЕНКА ---------- */

/* ФАЛЬСИФИКАТОР ЗДЕСЬ НЕ ПЕЧЬ, И ЭТО К55.
 *
 * План предлагал проверять зеркала замкнутой полостью с равновесием
 * `L = L_e/(1−ρ)`. Такой эталон СЛЕП: точное решение однородно, а в однородном
 * поле ЛЮБАЯ перестановка ординат — включая тождественную и заведомо неверную —
 * даёт тот же ответ. Это ровно К44, только применённая ДО кода, а не после.
 *
 * Проверяется СИММЕТРИЯ, и она неоднородна по построению. Источник несёт
 * направленную часть `ε += eps_dir·ω`, а отражение в плоскости `z` меняет знак
 * `ω_z`. Значит задача с зеркалом на стенке `−z` при `eps_dir = (0,0,+c)`
 * обязана дать РОВНО то же, что задача с зеркалом на `+z` при `(0,0,−c)`,
 * отражённая по `z`. У скалярного потока при этом меняет знак наклон по `z`, а
 * прочие коэффициенты сохраняются.
 *
 * Неверная перестановка отправляет свет не туда и ломает это на `O(1)`. */
static void t_mirror_wall(void) {
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 1.0}};
  const int L = 3, N = 1 << L;
  double res[2] = {1.0, 0.0};
  for (int broken = 0; broken < 2; broken++) {
    double *keep = NULL;
    for (int side = 0; side < 2; side++) {
      rig r;
      if (rig_init(&r, L, 0, 2, &fr)) {
        check(0, "оснастка C");
        free(keep);
        return;
      }
      /* НЕГАТИВНЫЙ КОНТРОЛЬ — И ОН ЛОМАЕТСЯ ТОЛЬКО В ОДНОЙ ИЗ ДВУХ ЗАДАЧ.
       *
       * Первая редакция подменяла перестановку тождественной В ОБЕИХ, и
       * контроль НЕ ВЫСТРЕЛИЛ: `4.6e-16` вместо `O(1)`. Причина содержательная
       * и стоит того, чтобы её записать: тождественная перестановка САМА
       * симметрична относительно `z`, поэтому сломанные одинаково задачи
       * по-прежнему переходят друг в друга, и метрика к поломке слепа. Это К55
       * в четвёртый раз — и на этот раз на СВОЁМ тесте.
       *
       * Ломается поэтому ОДНА сторона: тогда вопрос звучит так, как и должен, —
       * «зависит ли ответ от перестановки вообще». */
      if (broken && side == 0)
        for (int a = 0; a < 3; a++)
          for (int mm = 0; mm < r.d.n; mm++)
            r.d.mir[a * r.d.n + mm] = mm;
      for (int32_t c = 0; c < r.m.ncell; c++) {
        r.sig_t[c] = 0.30;
        r.sig_s[c] = 0.15;
        r.eps[c * 4] = 1.0;
      }
      double ws[6] = {0, 0, 0, 0, 0, 0};
      ws[side ? 5 : 4] = 1.0; /* зеркало на −z (side=0) либо на +z (side=1) */
      double wrho[6] = {0, 0, 0, 0, 0, 0};
      tr3_problem p = {.m = &r.m,
                       .d = &r.d,
                       .sig_t = r.sig_t,
                       .sig_s = r.sig_s,
                       .eps = r.eps,
                       .eps_dir = {0.0, 0.0, side ? -0.4 : 0.4},
                       .wall_rho = wrho,
                       .wall_spec = ws,
                       .limiter = 0};
      double *ph = calloc((size_t)r.m.ncell * 4, sizeof(double));
      if (ph == NULL) {
        check(0, "память");
        rig_free(&r);
        free(keep);
        return;
      }
      tr3_stats st;
      int rc = tr3_sweep_solve(&p, 2000, 1e-12, ph, &st);
      check(rc == 0, "C: развёртка с зеркальной стенкой прошла");
      /* К57: ЗЕРКАЛО ПЕРЕСТАЁТ БЫТЬ ЧЁРНЫМ, И ПРИБОР ОБЯЗАН ЭТО ПЕРЕЖИТЬ.
       * Прежде тождество `pin + psout = pout + pabs + psin` держалось отчасти
       * потому, что перехваченное зеркалом просто ПРОПАДАЛО. Теперь оно
       * возвращается, и если путь возврата не учтён, невязка станет ЗАКОННОЙ —
       * то есть прибор снова ослепнет, как до К40. */
      if (!broken) {
        double rel = fabs(st.balance) / (fabs(st.pin) + 1e-300);
        printf("    [C3] баланс при зеркальной стенке: втекло %.4f, вытекло %.4f, "
               "поглощено %.4f, невязка/втекло %.2e\n",
               st.pin, st.pout, st.pabs, rel);
        check(rel < 1e-10, "К57: тождество баланса держится и с зеркальной стенкой");
      }
      free(st.bout);
      free(st.sout);
      if (side == 0) {
        keep = ph;
      } else {
        double worst = 0.0, scale = 0.0;
        for (int32_t c = 0; c < r.m.ncell; c++) {
          int32_t zz = N - 1 - r.m.clo[c][2];
          int32_t c2 = r.m.cellof[hz_oct_leaf(&r.t, r.m.clo[c][0], r.m.clo[c][1], zz)];
          const double sgn[4] = {1.0, 1.0, 1.0, -1.0};
          for (int j = 0; j < 4; j++) {
            double e = fabs(keep[c * 4 + j] - sgn[j] * ph[c2 * 4 + j]);
            if (e > worst) worst = e;
            if (fabs(keep[c * 4 + j]) > scale) scale = fabs(keep[c * 4 + j]);
          }
        }
        res[broken] = worst / (scale > 0.0 ? scale : 1.0);
        printf("  [C2%s] зеркало −z против +z, отражённое: max отн. расхождение %.3e\n",
               broken ? ", НЕГАТИВНЫЙ КОНТРОЛЬ (перестановка тождественна)" : "", res[broken]);
        free(keep);
        keep = NULL;
        free(ph);
      }
      rig_free(&r);
    }
  }
  check(res[0] < 1e-12, "C: зеркальная стенка симметрична — перестановка верна");
  check(res[1] > 1e-2, "негативный контроль: тождественная перестановка ЛОМАЕТ симметрию");
}

static void t_balance(void) {
  /* К13: баланс энергии при НЕОДНОРОДНОЙ подсветке — то, чем ловится
   * неконсервативная передача потока, тогда как печь её не видит. */
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 2.0}};
  for (int graded = 0; graded < 2; graded++) {
    rig r;
    if (rig_init(&r, 3, graded, 2, &fr)) {
      check(0, "оснастка");
      rig_free(&r);
      return;
    }
    for (int32_t c = 0; c < r.m.ncell; c++) {
      r.sig_t[c] = 0.5;
      r.sig_s[c] = 0.5; /* альбедо 1: поглощения нет, всё обязано выйти */
    }
    /* влёт ЛИНЕЙНЫЙ по положению — неоднородная подсветка */
    tr3_problem p = {.m = &r.m,
                     .d = &r.d,
                     .sig_t = r.sig_t,
                     .sig_s = r.sig_s,
                     .binc0 = 1.0,
                     .binc = {0.13, -0.07, 0.05},
                     .binx0 = {4, 4, 4},
                     .limiter = 1};
    double *phi = calloc((size_t)r.m.ncell * 4, sizeof(double));
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    tr3_stats st;
    int rc = tr3_sweep_solve(&p, 2000, 1e-13, phi, &st);
    free(st.bout);
    free(st.sout);
    double rel = fabs(st.balance) / (fabs(st.pin) + 1e-300);
    printf("  [C баланс%s] код %d, итераций %d: втекло %.6f, вытекло %.6f, поглощено %.2e, "
           "невязка/втекло %.2e\n",
           graded ? ", градуированная" : "", rc, st.iters, st.pin, st.pout, st.pabs, rel);
    check(rc == 0, "развёртка прошла");
    check(rel < 1e-10, "К13: при альбедо 1 вся влитая мощность выходит");
    free(phi);
    rig_free(&r);
  }
}

/* -------------------------------------------------------------------------- */

/* ------------------------------- D: стенки как ГУ, замкнутая полость */

static void t_cavity(void) {
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 1.0}};
  /* Замкнутая диффузная полость: все стенки с отражением ρ и собственным
   * радиансом L_e, среда пустая. Классическая замкнутая форма бесконечной
   * серии отражений: равновесный радианс L = L_e/(1−ρ), ОДНОРОДНЫЙ. */
  const double rho = 0.6, le = 0.4;
  double wr[6], we[6];
  for (int i = 0; i < 6; i++) {
    wr[i] = rho;
    we[i] = le;
  }
  double exact = le / (1.0 - rho);
  for (int graded = 0; graded < 2; graded++) {
    rig r;
    if (rig_init(&r, 3, graded, 2, &fr)) {
      check(0, "оснастка");
      rig_free(&r);
      return;
    }
    tr3_problem p = {.m = &r.m,
                     .d = &r.d,
                     .sig_t = r.sig_t,
                     .sig_s = r.sig_s,
                     .wall_rho = wr,
                     .wall_emit = we,
                     .limiter = 1};
    double *phi = calloc((size_t)r.m.ncell * 4, sizeof(double));
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    if (phi == NULL) {
      check(0, "память");
      rig_free(&r);
      return;
    }
    tr3_stats st;
    int rc = tr3_sweep_solve(&p, 4000, 1e-14, phi, &st);
    double worst = 0.0, wb = 0.0;
    for (int32_t c = 0; c < r.m.ncell; c++) {
      double e = fabs(phi[c * 4] - 4.0 * M_PI * exact) / (4.0 * M_PI * exact);
      if (e > worst) worst = e;
    }
    /* ХРАНИМОЕ ЕСТЬ КОЭФФИЦИЕНТЫ DG1, А НЕ ЗНАЧЕНИЕ (К39). Читать `bout[f·4]`
     * как радианс НЕЛЬЗЯ: базис ячейки центрирован на ЯЧЕЙКЕ, а грань лежит на
     * её краю, поэтому у постоянного поля c представитель минимальной нормы
     * равен (0.8c, 0.4c, 0, 0) и значение собирается только вместе с базисом.
     * Прежде это сходилось случайно — откат к константе обнулял наклоны, то
     * есть тест проходил ровно потому, что проекция была сломана. */
    for (int32_t f = 0; f < r.m.nf; f++) {
      if (r.m.f[f].cb >= 0) continue;
      double v[4][3];
      tr3_face_corners(&r.m, f, v);
      int32_t c = r.m.f[f].ca;
      double s = (double)r.m.csize[c];
      for (int q = 0; q < 4; q++) {
        double val = st.bout[f * 4];
        for (int a = 0; a < 3; a++) {
          double xu = (v[q][a] - r.m.fr.o[a]) / r.m.fr.u[a];
          val += st.bout[f * 4 + a + 1] * ((xu - ((double)r.m.clo[c][a] + 0.5 * s)) / s);
        }
        double e = fabs(val - exact) / exact;
        if (e > wb) wb = e;
      }
    }
    printf("  [D полость%s] код %d, итераций %d: радианс отн. ошибка %.2e, "
           "исходящий на стенках %.2e (точно %.6f)\n",
           graded ? ", градуированная" : "", rc, st.iters, worst, wb, exact);
    check(rc == 0, "развёртка со стенками прошла");
    check(worst < 1e-10, "D: равновесие полости есть L_e/(1−ρ), однородно");
    check(wb < 1e-10, "и исходящий радианс стенок тот же");
    free(st.bout);
    free(st.sout);
    free(phi);
    rig_free(&r);
  }

  /* НЕГАТИВНЫЙ КОНТРОЛЬ (предсказан провал): забыть 1/π в ламбертовом ГУ.
   * Тогда «отражение» отдаёт в π раз больше, серия сходится к другому числу —
   * а при ρ ≥ 1/π вовсе расходится. Проверяем, что ответ УЕХАЛ. */
  rig r;
  if (rig_init(&r, 3, 0, 2, &fr)) {
    rig_free(&r);
    return;
  }
  double wr2[6];
  for (int i = 0; i < 6; i++)
    wr2[i] = rho * M_PI; /* ровно то, что выйдет без деления на π */
  tr3_problem p2 = {.m = &r.m,
                    .d = &r.d,
                    .sig_t = r.sig_t,
                    .sig_s = r.sig_s,
                    .wall_rho = wr2,
                    .wall_emit = we,
                    .limiter = 1};
  double *phi2 = calloc((size_t)r.m.ncell * 4, sizeof(double));
  if (phi2 == NULL) {
    rig_free(&r);
    return;
  }
  tr3_stats st2;
  tr3_sweep_solve(&p2, 200, 1e-12, phi2, &st2);
  free(st2.sout);
  double got = phi2[0] / (4.0 * M_PI);
  printf("    [НК без 1/π] радианс %.4f против правильного %.4f\n", got, exact);
  check(fabs(got - exact) > 0.5 * exact, "негативный контроль: без 1/π ответ ОБЯЗАН уехать");
  free(st2.bout);
  free(phi2);
  rig_free(&r);
}

/* ------------------------------- D: КАРТИНКА с косвенным светом */

static void t_render(void) {
  const int L = 3;
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 1.0}};
  rig r;
  if (rig_init(&r, L, 2, 2, &fr)) {
    check(0, "оснастка");
    rig_free(&r);
    return;
  }
  int n = 1 << L;
  /* комната: пять серых стенок и один яркий потолок */
  double wr[6] = {0.7, 0.7, 0.7, 0.7, 0.75, 0.7};
  double we[6] = {0, 0, 0, 0, 0, 3.0}; /* +z светит */
  for (int32_t c = 0; c < r.m.ncell; c++) {
    r.sig_t[c] = 0.02; /* лёгкая дымка, чтобы среда участвовала */
    r.sig_s[c] = 0.02;
  }
  tr3_problem p = {.m = &r.m,
                   .d = &r.d,
                   .sig_t = r.sig_t,
                   .sig_s = r.sig_s,
                   .wall_rho = wr,
                   .wall_emit = we,
                   .limiter = 1};
  double *phi = calloc((size_t)r.m.ncell * 4, sizeof(double));
  if (phi == NULL) {
    check(0, "память");
    rig_free(&r);
    return;
  }
  tr3_stats st;
  int rc = tr3_sweep_solve(&p, 2000, 1e-12, phi, &st);
  check(rc == 0, "комната посчитана");
  printf("  [D комната] итераций %d, ячеек %d, граней %d\n", st.iters, r.m.ncell, r.m.nf);

  /* ТЁПЛЫЙ СТАРТ — ИНКРЕМЕНТНОСТЬ МЕЖДУ КАДРАМИ, ИЗМЕРЕННАЯ, А НЕ ЗАЯВЛЕННАЯ.
   *
   * Архитектура держится на утверждении «поле не зависит от камеры, кадр есть
   * переизлучение уже известного», а рассеяние связывает сцену целиком: сдвиньте
   * что-нибудь — и поле, строго говоря, меняется везде. Локального приращения тут
   * нет и быть не может. Но оно и не нужно: `φ` прошлого кадра есть начальное
   * приближение, и `krylov3.h` его ПРИНИМАЕТ («на входе начальное приближение»).
   * Значит инкрементность выражается ЧИСЛОМ ПРОХОДОВ, а не числом тронутых
   * ячеек, и это ровно то, что здесь меряется.
   *
   * ВОЗМУЩЕНИЕ — ИЗМЕНЕНИЕ СЦЕНЫ, А НЕ КАМЕРЫ, и это сознательно: от камеры поле
   * не зависит вовсе, поэтому «сдвинули камеру» дало бы ноль проходов по
   * построению и ничего не проверило бы. Двигается яркость потолка на 5 %.
   *
   * ПРЕДСКАЗАНИЕ ДО ПРОГОНА: холодный старт здесь и так дёшев (BiCGStab берёт
   * единицы проходов), поэтому тёплый может не дать НИЧЕГО — и это будет
   * законный отрицательный ответ, а не неудача замера. Ценность его в том, что
   * он отделяет «инкрементность полезна» от «инкрементность не нужна, потому что
   * холодный старт дешёв». */
  {
    double we2[6] = {0, 0, 0, 0, 0, 3.15}; /* потолок ярче на 5 % — «следующий кадр» */
    tr3_problem p2 = p;
    p2.wall_emit = we2;
    int32_t nphi = r.m.ncell * 4, nb = r.m.nf * 4;
    double *pc = calloc((size_t)nphi, sizeof(double));
    double *pw = calloc((size_t)nphi, sizeof(double));
    double *b0 = calloc((size_t)(nb > 0 ? nb : 1), sizeof(double));
    double *bc = calloc((size_t)(nb > 0 ? nb : 1), sizeof(double));
    double *bw = calloc((size_t)(nb > 0 ? nb : 1), sizeof(double));
    double *p0 = calloc((size_t)nphi, sizeof(double));
    if (pc != NULL && pw != NULL && b0 != NULL && bc != NULL && bw != NULL && p0 != NULL) {
      tr3_kstats k0, kc, kw;
      /* кадр 1 */
      tr3_krylov_solve(&p, 200, 1e-10, p0, b0, NULL, &k0);
      /* кадр 2 ХОЛОДНЫЙ: с нуля */
      tr3_krylov_solve(&p2, 200, 1e-10, pc, bc, NULL, &kc);
      /* кадр 2 ТЁПЛЫЙ: от решения кадра 1 */
      memcpy(pw, p0, (size_t)nphi * sizeof(double));
      memcpy(bw, b0, (size_t)(nb > 0 ? nb : 1) * sizeof(double));
      tr3_krylov_solve(&p2, 200, 1e-10, pw, bw, NULL, &kw);
      double dmax = 0.0, sc = 0.0;
      for (int32_t i = 0; i < nphi; i++) {
        double a = fabs(pc[i]);
        if (a > sc) sc = a;
      }
      if (!(sc > 0.0)) sc = 1.0;
      for (int32_t i = 0; i < nphi; i++) {
        double e = fabs(pw[i] - pc[i]) / sc;
        if (e > dmax) dmax = e;
      }
      printf("  [ТЁПЛЫЙ СТАРТ] кадр 1: %ld проходов; кадр 2 С НУЛЯ: %ld; кадр 2 ОТ ПРОШЛОГО: %ld; "
             "расхождение ответов %.2e\n",
             (long)k0.npass, (long)kc.npass, (long)kw.npass, dmax);
      check(dmax < 1e-6, "ТЁПЛЫЙ СТАРТ: ответ ТОТ ЖЕ, иначе экономия куплена другим ответом");
    }
    free(pc);
    free(pw);
    free(b0);
    free(bc);
    free(bw);
    free(p0);
  }

  /* СБОР ПО ПИКСЕЛЮ: луч из камеры внутри комнаты до стенки, читаем ИСХОДЯЩИЙ
   * радианс той граничной грани, в которую попали. Интерполяции между гранями
   * нет намеренно — К3 запрещает её между ординатами, и здесь та же мысль. */
  const int W = 220, H = 165;
  double *buf = calloc((size_t)W * (size_t)H, sizeof(double));
  if (buf == NULL) {
    check(0, "память");
    free(st.bout);
    free(st.sout);
    free(phi);
    rig_free(&r);
    return;
  }
  double eye[3] = {(double)n * 0.5, (double)n * 0.12, (double)n * 0.45};
  double at[3] = {(double)n * 0.5, (double)n, (double)n * 0.42};
  double up[3] = {0, 0, 1};
  tr3_camera cam;
  check(tr3_camera_look(&cam, eye, at, up, 1.15, W, H) == 0, "камера");
  int miss = 0;
  for (int py = 0; py < H; py++)
    for (int px = 0; px < W; px++) {
      double o[3], d[3];
      tr3_camera_ray(&cam, px, py, o, d);
      /* выход из куба */
      double tex = 1e300;
      int wall = -1;
      for (int a = 0; a < 3; a++) {
        if (!(fabs(d[a]) > 0.0)) continue;
        double lim = d[a] > 0.0 ? (double)n : 0.0;
        double t = (lim - o[a]) / d[a];
        if (t > 0.0 && t < tex) {
          tex = t;
          wall = 2 * a + (d[a] > 0.0 ? 1 : 0);
        }
      }
      if (wall < 0) {
        miss++;
        continue;
      }
      double hp[3];
      for (int a = 0; a < 3; a++)
        hp[a] = o[a] + tex * d[a];
      /* граничная грань, накрывающая точку */
      int axis = wall / 2;
      int u = (axis + 1) % 3, v = (axis + 2) % 3;
      if (u > v) {
        int t = u;
        u = v;
        v = t;
      }
      double val = 0.0;
      for (int32_t f = 0; f < r.m.nf; f++) {
        const tr3_face *ff = &r.m.f[f];
        if (ff->cb >= 0 || (int)(~ff->cb) != wall) continue;
        if (hp[u] < (double)ff->lo[0] || hp[u] > (double)ff->hi[0]) continue;
        if (hp[v] < (double)ff->lo[1] || hp[v] > (double)ff->hi[1]) continue;
        /* DG1 ПО ПОЛОЖЕНИЮ: радианс на грани линейный, и читать его надо В ТОЧКЕ
         * попадания, а не константой на грань. Константа и давала квадраты. */
        {
          int32_t cc = ff->ca;
          double sz = (double)r.m.csize[cc];
          double bb[4] = {1.0, 0, 0, 0};
          for (int a = 0; a < 3; a++) {
            double xu = (hp[a] - r.m.fr.o[a]) / r.m.fr.u[a];
            bb[a + 1] = (xu - ((double)r.m.clo[cc][a] + 0.5 * sz)) / sz;
          }
          val = 0.0;
          for (int i = 0; i < 4; i++)
            val += st.bout[f * 4 + i] * bb[i];
        }
        break;
      }
      buf[py * W + px] = val;
    }
  double mn = 1e300, mx = 0.0;
  for (int i = 0; i < W * H; i++) {
    if (buf[i] < mn) mn = buf[i];
    if (buf[i] > mx) mx = buf[i];
  }
  printf("    радианс на стенках от %.4f до %.4f; пикселей без грани %d\n", mn, mx, miss);
  check(miss == 0, "каждый луч нашёл свою стенку");
  check(mn > 0.0, "теней в замкнутой комнате нет — косвенный свет достаёт везде");
  hz_ppm_write("img/room.ppm", buf, W, H);
  printf("    картинка записана: img/room.ppm\n");
  free(buf);
  free(st.bout);
  free(st.sout);
  free(phi);
  rig_free(&r);
}

/* ------------------------- D2: разрезанные ячейки, флюидная геометрия */

static void t_cut(void) {
  const hz_frame fr = {{0, 0, 0}, {1, 1, 1}};
  const int L = 4, N = 1 << L;
  hz_octree t;
  hz_oct_init(&t, L, 0.0);
  for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
      for (int z = 0; z < N; z++) {
        int lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        hz_oct_set_box(&t, lo, hi, 1.0);
      }
  hz_surftab stab;
  hz_facettab ftab;
  hz_cutmap cmap;
  hz_surftab_init(&stab);
  hz_facettab_init(&ftab);
  hz_cutmap_init(&cmap);
  double sc[3] = {8.3, 8.7, 7.6}, sr = 3.4;
  hz_surf sp = {HZ_SURF_SPHERE, {sc[0], sc[1], sc[2], sr, 0, 0, 0}, 1, 0};
  int32_t si = hz_surftab_add(&stab, &sp);
  int32_t f0 = 0;
  int32_t nfac = hz_surf_facet_sphere(&ftab, &fr, sc, sr, 1, HZ_FIT_MEAN_SAGITTA, si, &f0);
  check(nfac > 0, "сфера фасетизирована");
  /* Г45 В ДЕЙСТВИИ: ключ боковой таблицы обязан СТРОГО ВОЗРАСТАТЬ, а индексы
   * узлов идут порядком ВЫДЕЛЕНИЯ, а не порядком (x,y,z). Первая редакция этого
   * теста добавляла записи как попало и не смотрела на код возврата — половина
   * молча отвергалась, тело выходило дырявым, и объём материала был 105 вместо
   * 165. Ровно тот сценарий, про который Г45 и написана. */
  typedef struct {
    int32_t cell, f[HZ_P3_MAXH], nf;
  } rec_t;
  static rec_t rec[4096];
  int nrec = 0;
  for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
      for (int z = 0; z < N; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        int ns = hz_facets_for_box(&ftab, f0, nfac, lo, hi, sel, HZ_P3_MAXH);
        if (ns <= 0) continue;
        rec[nrec].cell = hz_oct_leaf(&t, x, y, z);
        rec[nrec].nf = ns;
        for (int j = 0; j < ns; j++)
          rec[nrec].f[j] = sel[j];
        nrec++;
      }
  for (int i = 1; i < nrec; i++) {
    rec_t tmp = rec[i];
    int j = i - 1;
    while (j >= 0 && rec[j].cell > tmp.cell) {
      rec[j + 1] = rec[j];
      j--;
    }
    rec[j + 1] = tmp;
  }
  int nadd = 0;
  for (int i = 0; i < nrec; i++)
    if (hz_cutmap_add(&cmap, rec[i].cell, rec[i].f, rec[i].nf) == 0) nadd++;
  check(nadd == nrec, "Г45: ВСЕ записи легли в боковую таблицу (ключи по возрастанию)");
  tr3_mesh m;
  check(tr3_mesh_build(&m, &t, &fr) == 0, "сетка");
  uint8_t *solid = calloc((size_t)m.ncell, 1);
  for (int x = 0; x < N; x++)
    for (int y = 0; y < N; y++)
      for (int z = 0; z < N; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        int ns = hz_facets_for_box(&ftab, f0, nfac, lo, hi, sel, HZ_P3_MAXH);
        if (ns == 0) solid[m.cellof[hz_oct_leaf(&t, x, y, z)]] = 1; /* Г38: 0 = ПОЛНАЯ */
      }
  tr3_cut cu;
  check(tr3_cut_build(&cu, &m, &ftab, &cmap, solid) == 0, "разрез");
  check(cu.nbad == 0, "условие 1:1 у разрезанных ячеек соблюдено");

  /* 1. ОБЪЁМ МАТЕРИАЛА: коробки минус флюид = объём шара, аналитически */
  double vmat = 0.0;
  for (int32_t c = 0; c < m.ncell; c++) {
    double s = (double)m.csize[c];
    vmat += s * s * s - cu.mvol[c][0][0];
  }
  free(solid);
  double vex = 4.0 / 3.0 * M_PI * sr * sr * sr;
  /* РЕШАЮЩАЯ СВЕРКА: материал считается НЕЗАВИСИМО, проверенным путём
   * hz_poly3_cut (пункт 4), и сумма «флюид + материал = коробка» обязана
   * держаться поячеечно. Если нет — виновата флюидная бухгалтерия, а не
   * фасетизация. */
  double vmat2 = 0.0, wcell = 0.0;
  for (int32_t c = 0; c < m.ncell; c++) {
    const hz_cutrec *rr = hz_cutmap_find(&cmap, m.node[c]);
    double s = (double)m.csize[c];
    if (rr == NULL) {
      vmat2 += cu.solid[c] ? s * s * s : 0.0;
      continue;
    }
    hz_hspace hh[HZ_P3_MAXH];
    int32_t hd[HZ_P3_MAXH];
    int nh = hz_cutmap_hspaces(&ftab, &cmap, rr, hh, hd, NULL, HZ_P3_MAXH);
    int32_t lo2[3] = {m.clo[c][0], m.clo[c][1], m.clo[c][2]};
    int32_t hi2[3] = {lo2[0] + m.csize[c], lo2[1] + m.csize[c], lo2[2] + m.csize[c]};
    hz_poly3 pp;
    double vm = 0.0;
    if (hz_poly3_cut(&pp, lo2, hi2, hh, hd, nh) == HZ_P3_OK) vm = hz_poly3_volume(&pp, &fr);
    vmat2 += vm;
    double e = fabs(vm + cu.mvol[c][0][0] - s * s * s);
    if (e > wcell) wcell = e;
  }
  printf("      материал независимым путём %.5f; max |флюид+материал−коробка| = %.3e\n", vmat2,
         wcell);
  check(wcell < 1e-12, "поячеечно: флюид и материал дополняют коробку ТОЧНО");
  printf("  [D2] объём материала %.5f против шара %.5f (отн. %.3e), элементов %d\n", vmat, vex,
         fabs(vmat - vex) / vex, cu.nse);
  check(fabs(vmat - vex) / vex < 0.02, "флюид дополняет материал до шара (фасетизация k=1)");

  /* 2. ЗАМКНУТОСТЬ ФЛЮИДНОЙ ОБЛАСТИ: Σ(внешняя нормаль × площадь) = 0.
   * Это и есть проверка того, что ни одна грань не потеряна и ни одна не
   * посчитана дважды — теорема о дивергенции для постоянного поля. */
  double worst = 0.0;
  int nchecked = 0;
  for (int32_t c = 0; c < m.ncell; c++) {
    if (cu.solid[c]) continue;
    if (cu.sestart[c + 1] == cu.sestart[c]) continue; /* только разрезанные */
    double acc[3] = {0, 0, 0};
    for (int32_t k = m.fstart[c]; k < m.fstart[c + 1]; k++) {
      int32_t f = m.flist[k];
      int mine_is_a = m.f[f].ca == c;
      double sgn;
      if (m.f[f].cb < 0) {
        int wall = (int)(~m.f[f].cb);
        sgn = (wall & 1) ? 1.0 : -1.0;
      } else {
        sgn = mine_is_a ? 1.0 : -1.0;
      }
      double a = mine_is_a ? cu.ffm[f][0][0] : cu.ffmb[f][0][0];
      acc[m.f[f].axis] += sgn * a;
    }
    for (int32_t k = cu.sestart[c]; k < cu.sestart[c + 1]; k++) {
      const tr3_selem *se = &cu.se[cu.selist[k]];
      for (int a = 0; a < 3; a++)
        acc[a] += -se->n[a] * se->area; /* наружу ФЛЮИДА = −n */
    }
    double e = fabs(acc[0]) + fabs(acc[1]) + fabs(acc[2]);
    if (e > worst) worst = e;
    nchecked++;
  }
  printf("      замкнутость флюида: %d разрезанных ячеек, max |Σn·A| = %.3e\n", nchecked, worst);
  check(nchecked > 50, "разрезанных ячеек достаточно");
  check(worst < 1e-12, "D2: флюидная область ЗАМКНУТА — ни одна грань не потеряна");

  /* ---------------------------------------------------------------- К39 ---
   * ВЫРОЖДЕНИЕ МАТРИЦЫ МАСС ПЛОСКОГО ЭЛЕМЕНТА. Печь и полость этого не ловят
   * ВООБЩЕ: у однородного равновесия все наклоны нули, и откат к константе там
   * безвреден. Это ровно урок К12, повторённый на поверхностях, поэтому
   * проверка ставится ПРЯМАЯ. */
  {
    double wnul = 0.0, wmom = 0.0, wgrow = 0.0, wgrow_old = 0.0;
    int nneg = 0, nneg_old = 0, ncheck2 = 0;
    /* пробная линейная функция: положительна на всей ячейке (при |b| ≤ ½
     * минимум равен 1 − ½(0.1+0.2+0.15) = 0.775 > 0), поэтому у ИСПРАВНОЙ
     * проекции повода отступать к константе нет ни у одного элемента */
    const double ex[4] = {1.0, 0.1, -0.2, 0.15};
    for (int32_t e = 0; e < cu.nse + m.nf; e++) {
      const double (*M)[4];
      double nul[4];
      if (e < cu.nse) {
        M = cu.se[e].m;
        memcpy(nul, cu.se[e].nul, sizeof nul);
      } else {
        int32_t f = e - cu.nse;
        if (m.f[f].cb >= 0) continue;
        M = cu.ffm[f];
        tr3_face_null(&m, f, m.f[f].ca, nul);
      }
      if (!(M[0][0] > 0.0)) continue;
      ncheck2++;
      /* 1. НУЛЕВОЙ ВЕКТОР ВЕРЕН И ВЫРОЖДЕНИЕ ТОЧНОЕ */
      double nm = 0.0, nv = 0.0, r = 0.0;
      for (int i = 0; i < 4; i++) {
        nv += nul[i] * nul[i];
        double acc = 0.0;
        for (int j = 0; j < 4; j++) {
          acc += M[i][j] * nul[j];
          nm += M[i][j] * M[i][j];
        }
        r += acc * acc;
      }
      double rel = sqrt(r) / (sqrt(nm) * sqrt(nv));
      if (rel > wnul) wnul = rel;
      /* 2. ПРОЕКЦИЯ ВОСПРОИЗВОДИТ ЛИНЕЙНОЕ ПОЛЕ: моменты обязаны совпасть */
      double rr2[4] = {0, 0, 0, 0}, ee[4], nr = 0.0;
      for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++)
          rr2[i] += M[i][j] * ex[j];
        nr += rr2[i] * rr2[i];
      }
      nr = sqrt(nr);
      if (tr3_project_plane(M, rr2, nul, ee) != 0) continue;
      double dm = 0.0;
      for (int i = 0; i < 4; i++) {
        double acc = -rr2[i];
        for (int j = 0; j < 4; j++)
          acc += M[i][j] * ee[j];
        dm += acc * acc;
      }
      if (nr > 0.0 && sqrt(dm) / nr > wmom) wmom = sqrt(dm) / nr;
      /* 3. ПРЕДСТАВИТЕЛЬ ОГРАНИЧЕН, И ИМЕННО ЭТО ЛЕЧИТ ОТКАТ. Добавка вдоль
       * нулевого вектора на плоскости НЕ ВИДНА, но вне её растёт без предела, и
       * проверка положительности по углам ячейки съезжает в минус. */
      double ne = 0.0, nex = 0.0;
      for (int i = 0; i < 4; i++) {
        ne += ee[i] * ee[i];
        nex += ex[i] * ex[i];
      }
      double gr = sqrt(ne / nex);
      if (gr > wgrow) wgrow = gr;
      double mn = 1e300;
      for (int k = 0; k < 8; k++) {
        double v = ee[0];
        for (int a = 0; a < 3; a++)
          v += ee[a + 1] * (((k >> a) & 1) ? 0.5 : -0.5);
        if (v < mn) mn = v;
      }
      if (mn < 0.0) nneg++;
      /* НЕГАТИВНЫЙ КОНТРОЛЬ: ПРЕЖНИЙ КОД, воспроизведённый здесь дословно.
       * `solve4` вырождения не замечает — четвёртый ведущий элемент выходит не
       * нулём, а округлением, — и возвращает произвольную добавку вдоль nul. */
      double a4[4][4], b4[4], x4[4];
      memcpy(a4, M, sizeof a4);
      memcpy(b4, rr2, sizeof b4);
      int sing = 0;
      for (int k = 0; k < 4 && !sing; k++) {
        int best = k;
        for (int i = k + 1; i < 4; i++)
          if (fabs(a4[i][k]) > fabs(a4[best][k])) best = i;
        if (!(fabs(a4[best][k]) > 0.0)) {
          sing = 1;
          break;
        }
        if (best != k) {
          for (int j = 0; j < 4; j++) {
            double tt = a4[k][j];
            a4[k][j] = a4[best][j];
            a4[best][j] = tt;
          }
          double tt = b4[k];
          b4[k] = b4[best];
          b4[best] = tt;
        }
        for (int i = k + 1; i < 4; i++) {
          double ff = a4[i][k] / a4[k][k];
          for (int j = k; j < 4; j++)
            a4[i][j] -= ff * a4[k][j];
          b4[i] -= ff * b4[k];
        }
      }
      if (sing) continue;
      for (int i = 3; i >= 0; i--) {
        double s2 = b4[i];
        for (int j = i + 1; j < 4; j++)
          s2 -= a4[i][j] * x4[j];
        x4[i] = s2 / a4[i][i];
      }
      double ne2 = 0.0, mn2 = 1e300;
      for (int i = 0; i < 4; i++)
        ne2 += x4[i] * x4[i];
      if (sqrt(ne2 / nex) > wgrow_old) wgrow_old = sqrt(ne2 / nex);
      for (int k = 0; k < 8; k++) {
        double v = x4[0];
        for (int a = 0; a < 3; a++)
          v += x4[a + 1] * (((k >> a) & 1) ? 0.5 : -0.5);
        if (v < mn2) mn2 = v;
      }
      if (mn2 < 0.0) nneg_old++;
    }
    printf("  [К39] элементов %d: max |M·nul|/(|M||nul|) = %.3e\n", ncheck2, wnul);
    printf("        проекция линейного поля: max невязка моментов %.3e, ‖ee‖/‖ex‖ ≤ %.3f, "
           "минус по углам у %d\n",
           wmom, wgrow, nneg);
    printf("    [НК прежний solve4] ‖ee‖/‖ex‖ до %.3e, минус по углам у %d из %d\n", wgrow_old,
           nneg_old, ncheck2);
    check(wnul < 1e-13,
          "К39: матрица масс плоского элемента вырождена ТОЧНО, нулевой вектор верен");
    check(wmom < 1e-9, "К39: проекция воспроизводит линейное поле — моменты совпали");
    check(wgrow < 1.0 + 1e-12, "К39: представитель МИНИМАЛЬНОЙ нормы, добавки вдоль nul нет");
    check(nneg == 0, "К39: у положительного поля откат к константе не нужен НИ РАЗУ");
    /* ПРЕДСКАЗАНИЕ «больше половины элементов уйдут в минус» НЕ СБЫЛОСЬ:
     * измерено 423 из 2212, то есть 19%. Причина в том, что здесь подаётся
     * ТОЧНО ЛИНЕЙНОЕ поле, у которого моменты представимы без остатка, и
     * величина паразитной добавки определяется только обусловленностью
     * конкретного элемента. В настоящей сцене облучённость линейной не бывает, и
     * там доля вышла 1536 из 1536 у стенок и 1372 из 1912 у поверхностей.
     * Существо контроля от этого не меняется и даже сильнее: норма представителя
     * раздувается до 1.4e17, то есть решение вырожденной системы бессмысленно
     * само по себе. Проверка переписана на ИЗМЕРЕННОЕ существо. */
    check(nneg_old > 100, "НК: прежний solve4 ОБЯЗАН уводить угловой минимум в минус");
    check(wgrow_old > 1e6, "НК: и раздувать норму представителя на порядки");
  }

  /* ------------------------------------------------------------- К47 ---
   * ФЛЮИДНАЯ ПЛОЩАДЬ НА ГРАНИ, ВТОРАЯ СТОРОНА КОТОРОЙ ПОЛНОСТЬЮ ТВЁРДАЯ.
   * У полностью твёрдой ячейки развёртка ставит L = 0 и идёт дальше, а её сосед
   * по флюиду всё равно считает выток через общую грань по СВОЕЙ флюидной
   * площади. Если та не ноль — энергия уходит в никуда, и это утечка, которую
   * баланс К40 обязан видеть. Геометрически площадь ОБЯЗАНА быть нулевой: если
   * соседняя ячейка целиком в материале, то общая грань лежит внутри материала,
   * и флюида на ней нет. Расхождение здесь есть расхождение ДВУХ ВЕЕРОВ
   * фасетов, отобранных для двух ячеек независимо. */
  {
    double leak = 0.0, total = 0.0;
    int nleak = 0;
    for (int32_t f = 0; f < m.nf; f++) {
      if (m.f[f].cb < 0) continue;
      total += cu.ffm[f][0][0] + cu.ffmb[f][0][0];
      int sa = cu.solid[m.f[f].ca], sb = cu.solid[m.f[f].cb];
      if (sa == sb) continue;
      double a = sa ? cu.ffmb[f][0][0] : cu.ffm[f][0][0]; /* площадь со стороны ФЛЮИДА */
      if (a > 0.0) {
        leak += a;
        nleak++;
      }
    }
    printf("  [К47] флюидная площадь на границе с ПОЛНОСТЬЮ ТВЁРДОЙ ячейкой: %.6e на %d гранях "
           "(вся площадь граней %.4f, доля %.3e)\n",
           leak, nleak, total, total > 0.0 ? leak / total : 0.0);
    check(!(leak > 0.0), "К47: флюид не граничит с полностью твёрдой ячейкой");
  }

  /* ------------------------------------------------------------- К48 ---
   * ДВЕ ПОЛУСФЕРНЫЕ СУММЫ У НАКЛОННОЙ НОРМАЛИ НЕ РАВНЫ. К29 велит нормировать
   * отражение на СОБСТВЕННУЮ сумму набора, а не на π из континуума. Но сумм
   * ДВЕ — по входящей полусфере и по исходящей, — и у произвольной нормали они
   * РАЗНЫЕ, потому что в полусферы попадают разные направления. Отданная
   * поверхностью мощность равна `Σ_{ω·n>0} w(ω·n) · ∫L_out`, значит делить надо
   * на ИСХОДЯЩУЮ сумму. У стенки, перпендикулярной оси, обе совпадают по
   * симметрии набора, и на равномерной сетке ошибка не проявляется ВООБЩЕ. */
  {
    tr3_dirs d2;
    if (tr3_dirs_product(&d2, 2, 2) == 0) {
      double worstr = 0.0;
      for (int32_t e = 0; e < cu.nse; e++) {
        double hin = 0.0, hout = 0.0;
        for (int mm = 0; mm < d2.n; mm++) {
          double on =
              d2.ox[mm] * cu.se[e].n[0] + d2.oy[mm] * cu.se[e].n[1] + d2.oz[mm] * cu.se[e].n[2];
          if (on < 0.0)
            hin += d2.w[mm] * (-on);
          else
            hout += d2.w[mm] * on;
        }
        if (hin > 0.0 && fabs(hout / hin - 1.0) > worstr) worstr = fabs(hout / hin - 1.0);
      }
      printf("  [К48] max |Σ_исх / Σ_вх − 1| по %d элементам: %.3e\n", cu.nse, worstr);
      tr3_dirs_free(&d2);
    }
  }

  /* ------------------------------------------------------------- К38 ---
   * СХОДИМОСТЬ ПРИ ОТРАЖАЮЩЕЙ ПОВЕРХНОСТИ ВНУТРИ ОБЛАСТИ. Прежде здесь было
   * автоколебание отката, и развёртка упиралась в maxit при невязке 1.25.
   * Ловится это ТОЛЬКО так: у полости точный ответ однороден, наклоны нулевые,
   * и переключателю не на чем колебаться. */
  {
    tr3_dirs d;
    check(tr3_dirs_product(&d, 2, 2) == 0, "ординаты");
    double *sig_t = calloc((size_t)m.ncell, sizeof(double));
    double *sig_s = calloc((size_t)m.ncell, sizeof(double));
    double *frho = calloc((size_t)ftab.n, sizeof(double));
    double *femit = calloc((size_t)ftab.n, sizeof(double));
    double *phi = calloc((size_t)m.ncell * 4, sizeof(double));
    if (sig_t == NULL || sig_s == NULL || frho == NULL || femit == NULL || phi == NULL) {
      check(0, "память");
    } else {
      for (int32_t c = 0; c < m.ncell; c++) {
        sig_t[c] = 0.015;
        sig_s[c] = 0.012;
      }
      for (int32_t i = 0; i < ftab.n; i++)
        frho[i] = 0.78;
      double wr[6] = {0.72, 0.35, 0.72, 0.72, 0.65, 0.05};
      double we[6] = {0, 0, 0, 0, 0, 6.0};
      tr3_problem p = {.m = &m,
                       .d = &d,
                       .cut = &cu,
                       .facet_rho = frho,
                       .facet_emit = femit,
                       .nfacet = ftab.n,
                       .sig_t = sig_t,
                       .sig_s = sig_s,
                       .wall_rho = wr,
                       .wall_emit = we,
                       .limiter = 1};
      /* СВИП ПО ДОПУСКУ: отставание итерации или НАСТОЯЩАЯ утечка?
       * Тождество К13 верно на любой итерации ДЛЯ ОБЪЁМА, но `sout` считается
       * из облучённости ПРЕДЫДУЩЕЙ итерации, поэтому поверхностный член отстаёт
       * на шаг. Если невязка падает вместе с допуском — это отставание и оно
       * уходит в пределе; если упирается в полку — это утечка, и её надо искать.
       * Различить их иначе нечем, и именно поэтому свип печатается (тот же приём,
       * что К22 у толщины слоя и Ф5 у числа отскоков). */
      for (int k = 0; k < 3; k++) {
        double tl = k == 0 ? 1e-7 : (k == 1 ? 1e-9 : 1e-11);
        tr3_stats s2;
        memset(phi, 0, (size_t)m.ncell * 4 * sizeof(double));
        if (tr3_sweep_solve(&p, 500, tl, phi, &s2) == 0)
          printf("        [К40 свип] допуск %.0e: итераций %3d, невязка баланса на втекшее %.3e\n",
                 tl, s2.iters, fabs(s2.balance) / s2.pin);
        free(s2.bout);
        free(s2.sout);
      }
      /* ЕДИНСТВЕННАЯ НЕЛИНЕЙНОСТЬ В ЦИКЛЕ — ОГРАНИЧИТЕЛЬ. Его замкнутая форма
       * (T2) обязана сохранять уравнение при `v = 1` ТОЧНО, то есть быть
       * консервативной. Проверяется прямо: с ним и без него. */
      for (int k = 0; k < 2; k++) {
        tr3_problem p2 = p;
        p2.limiter = k;
        tr3_stats s2;
        memset(phi, 0, (size_t)m.ncell * 4 * sizeof(double));
        if (tr3_sweep_solve(&p2, 500, 1e-10, phi, &s2) == 0)
          printf("        [К40 ограничитель %s] срезок %d, невязка баланса на втекшее %.3e\n",
                 k ? "ВКЛ" : "ВЫКЛ", s2.nclip, fabs(s2.balance) / s2.pin);
        free(s2.bout);
        free(s2.sout);
      }
      memset(phi, 0, (size_t)m.ncell * 4 * sizeof(double));
      tr3_stats st;
      int rc = tr3_sweep_solve(&p, 500, 1e-9, phi, &st);
      printf("  [К38] комната с ОТРАЖАЮЩЕЙ сферой: код %d, итераций %d, невязка %.2e, "
             "откатов проекции на последней итерации %d из %d\n",
             rc, st.iters, st.resid, st.nfallback, cu.nse + m.nf);
      /* К40: БАЛАНС С УЧЁТОМ ПОВЕРХНОСТЕЙ. Здесь ОДНО тело, то есть случая
       * «два тела в одной ячейке» нет по построению, и утечка, если она есть,
       * принадлежит самой схеме, а не построителю сцены. */
      printf("        [К40] баланс: втекло %.4f, вытекло %.4f, поглощено средой %.4f, "
             "ушло в поверхности %.4f, отдано %.4f; невязка на втекшее %.3e\n",
             st.pin, st.pout, st.pabs, st.psin, st.psout, fabs(st.balance) / st.pin);
      check(fabs(st.balance) / st.pin < 1e-9,
            "К40: баланс с поверхностями — точное дискретное тождество");
      check(rc == 0, "развёртка прошла");
      /* ИСХОДОВ ТРИ, А НЕ ДВА — К81/К84. Прежде здесь стояло `st.resid < 1e-9`,
       * то есть «сошлось по допуску», и это было верно ровно до тех пор, пока
       * невязка мерила ОДНО `φ`. Теперь она мерит всё состояние `(φ, bout, sout)`
       * — картинка делается из последних двух, и не смотреть на них значило не
       * измерять то, ради чего решение считается, — а у поверхностной части свой
       * ПОЛ. Поэтому законных исходов два: допуск и ЗАСТОЙ. Незаконен один —
       * упереться в `maxit`, и ровно это ловила К38.
       *
       * Сливать их нельзя (Г38): «сошлось» и «упёрлось в пол» — разные ответы. */
      check(st.iters < 499, "К38: остановка ПО СУЩЕСТВУ, а не упор в maxit");
      check(st.resid < 1e-9 || st.stalled, "К38: исход назван — либо допуск, либо застой");
      check(st.iters < 200, "К38: сходимость геометрическая — десятки итераций, а не тысячи");
      check(st.nfallback * 20 < cu.nse, "К38: откат проекции — редкое исключение, а не правило");
      free(st.bout);
      free(st.sout);
    }
    free(sig_t);
    free(sig_s);
    free(frho);
    free(femit);
    free(phi);
    tr3_dirs_free(&d);
  }

  tr3_cut_free(&cu);
  tr3_mesh_free(&m);
  hz_cutmap_free(&cmap);
  hz_facettab_free(&ftab);
  hz_surftab_free(&stab);
  hz_oct_free(&t);
}

/* ------------------------- КРЫЛОВСКОЕ УСКОРЕНИЕ ИТЕРАЦИИ ПО РАССЕЯНИЮ ------ */

/* ИТЕРАЦИЯ ПО РАССЕЯНИЮ ЕСТЬ РЯД НЕЙМАНА — САМЫЙ МЕДЛЕННЫЙ СПОСОБ РЕШИТЬ
 * `(I − T)φ = b`. Он просто складывает порядки рассеяния, и его скорость равна
 * альбедо: `‖φⁿ − φ‖ ~ ρⁿ`. Крыловский метод на том же операторе идёт по `√κ`
 * вместо `κ`, где `κ = 1/(1−ρ)`, и выигрыш тем больше, чем ближе альбедо к
 * единице — то есть ровно там, где С1 и предупреждал про сотни итераций.
 *
 * ОПЕРАТОР СТРОИТСЯ ИЗ САМОЙ РАЗВЁРТКИ, а не выписывается заново:
 *
 *     S(φ) — один проход развёртки (при `warm_start`, maxit = 1)
 *     b = S(0)                       — вклад предписанного влёта
 *     T·v = S(v) − b                 — линейная часть
 *     A·v = v − T·v = v − S(v) + b
 *
 * Одно умножение на `A` есть ОДИН проход развёртки, поэтому честное сравнение
 * идёт по ЧИСЛУ ПРОХОДОВ, а не по числу итераций метода: BiCGStab делает два
 * прохода на итерацию.
 *
 * ЛОВУШКА, НАЗВАННАЯ ДО КОДА: Крылову нужен ЛИНЕЙНЫЙ оператор, а в проходе
 * сидит ограничитель положительности. Крыловские векторы физическим радиансом
 * не являются, у них есть отрицательные компоненты, и ограничитель на них
 * сработает — оператор перестанет быть линейным. Поэтому опыт ставится с
 * ВЫКЛЮЧЕННЫМ ограничителем, и законно это ровно на печи: там на сошедшемся
 * решении `nclip = 0` (измерено, К31), значит выключение ответа не меняет.
 * Отражающих границ и поверхностных элементов в печи нет, других нелинейностей
 * не остаётся. Для сцен с отражением так делать НЕЛЬЗЯ — там ограничитель
 * активен и на сходимости (458 срезок в рендере), и нужен другой ход. */

typedef struct {
  tr3_problem p;
  double *work;
  long npass;
  int32_t n;
} krylov_op;

/* v_out = A·v = v − S(v) + b. Разрушает `work`. */
static void kr_apply(krylov_op *k, const double *v, const double *b, double *out) {
  tr3_stats st;
  memcpy(k->work, v, (size_t)k->n * sizeof(double));
  if (tr3_sweep_solve(&k->p, 1, 0.0, k->work, &st) != 0) return;
  k->npass++;
  free(st.bout);
  free(st.sout);
  for (int32_t i = 0; i < k->n; i++)
    out[i] = v[i] - k->work[i] + b[i];
}

static double kr_dot(const double *a, const double *b, int32_t n) {
  double s = 0.0;
  for (int32_t i = 0; i < n; i++)
    s += a[i] * b[i];
  return s;
}

static double kr_nrm(const double *a, int32_t n) {
  double m = 0.0;
  for (int32_t i = 0; i < n; i++)
    if (fabs(a[i]) > m) m = fabs(a[i]);
  return m;
}

static void t_krylov(void) {
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 2.0}};
  const double lb = 1.7;
  (void)lb;
  /* альбедо РОВНО 1 — тот самый случай, ради которого С1 заводил разгон */
  /* СВИП ПО ОПТИЧЕСКОЙ ТОЛЩИНЕ ЯЧЕЙКИ, потому что от неё всё и зависит.
   * При `σ_t·h ~ 1` свободный пробег равен ячейке, рассеянный свет почти не
   * уходит из неё, оператор рассеяния ПОЧТИ ДИАГОНАЛЕН, спектр сжат — и Крылов
   * сносит его за один шаг. Это измерено: на случайном векторе `cos(v, A·v)`
   * вышел `0.997`. Такой режим ускорение ЗАВЫШАЕТ, и мерить только на нём
   * нельзя. В нашем рендере `σ_t·h = 0.015` — другой конец шкалы. */
  const double sigh[7] = {0.05, 0.8, 4.0, 0.8, 0.05, 0.8, 4.0};
  for (int cs = 0; cs < 7; cs++) {
    rig r;
    if (rig_init(&r, 3, 0, 2, &fr)) {
      check(0, "оснастка");
      rig_free(&r);
      return;
    }
    /* cs = 0 — НЕОДНОРОДНАЯ среда при альбедо 1. Однородную печь для замера
     * УСКОРЕНИЯ брать нельзя: у неё точное решение однородно, оно почти
     * параллельно правой части, и BiCGStab сходится за один шаг по причине,
     * не имеющей отношения к скорости метода. Это тот же урок, что К12 и К44 —
     * эталон с однородным решением проверяет только однородные степени свободы.
     * cs = 1 — НЕГАТИВНЫЙ КОНТРОЛЬ без рассеяния. */
    double alb = cs == 3 ? 0.0 : 1.0; /* cs = 3 — НК; cs = 4 — СЛУЧАЙНАЯ правая часть */
    for (int32_t c = 0; c < r.m.ncell; c++) {
      /* σ меняется по ячейкам в 6 раз — решение заведомо неоднородно */
      double f = 0.3 + 1.7 * (double)((c * 37) % 11) / 10.0;
      r.sig_t[c] = sigh[cs] * f;
      r.sig_s[c] = sigh[cs] * f * alb;
    }
    /* ИСТОЧНИК ОБЯЗАН ВОЗБУЖДАТЬ МНОГО МОД, иначе замер вырожден.
     * Первая редакция брала ОДНОРОДНЫЙ влёт на границе, и незатенённое поле от
     * него оказалось почти точно ОСНОВНОЙ МОДОЙ оператора рассеяния: измерено
     * `cos(b, A·b) = 1.000000` и `T·b = 0.524·b`. Крылов такую задачу решает за
     * один шаг тождественно, и ускорение на ней измерить нельзя.
     * Здесь источник ЛОКАЛЬНЫЙ — светится одна ячейка из восьми по каждой оси,
     * и в правой части оказывается широкий спектр мод. */
    for (int32_t c = 0; c < r.m.ncell; c++)
      r.eps[c * 4] = ((r.m.clo[c][0] + r.m.clo[c][1] * 3 + r.m.clo[c][2] * 7) % 8 == 0) ? 5.0 : 0.0;
    int32_t n = r.m.ncell * 4;
    tr3_problem base = {
        .m = &r.m, .d = &r.d, .sig_t = r.sig_t, .sig_s = r.sig_s, .eps = r.eps, .limiter = 0};

    /* --- эталон: обычная итерация по рассеянию --- */
    double *ref = calloc((size_t)n, sizeof(double));
    tr3_stats st0;
    tr3_sweep_solve(&base, 2000, 1e-12, ref, &st0);
    free(st0.bout);
    free(st0.sout);
    int passes_ref = st0.iters + 1;

    /* --- крыловский путь --- */
    krylov_op k = {.p = base, .work = calloc((size_t)n, sizeof(double)), .npass = 0, .n = n};
    k.p.warm_start = 1;
    double *b = calloc((size_t)n, sizeof(double));
    double *x = calloc((size_t)n, sizeof(double));
    double *rr = calloc((size_t)n, sizeof(double));
    double *rh = calloc((size_t)n, sizeof(double));
    double *pv = calloc((size_t)n, sizeof(double));
    double *vv = calloc((size_t)n, sizeof(double));
    double *ss = calloc((size_t)n, sizeof(double));
    double *tt = calloc((size_t)n, sizeof(double));
    if (k.work == NULL || b == NULL || x == NULL || rr == NULL || rh == NULL || pv == NULL ||
        vv == NULL || ss == NULL || tt == NULL) {
      check(0, "память");
      free(ref);
      free(k.work);
      free(b);
      free(x);
      free(rr);
      free(rh);
      free(pv);
      free(vv);
      free(ss);
      free(tt);
      rig_free(&r);
      return;
    }
    /* b = S(0): один проход от нулевого поля */
    memset(k.work, 0, (size_t)n * sizeof(double));
    {
      tr3_stats st1;
      tr3_sweep_solve(&k.p, 1, 0.0, k.work, &st1);
      k.npass++;
      free(st1.bout);
      free(st1.sout);
      memcpy(b, k.work, (size_t)n * sizeof(double));
    }
    /* РЕШАЮЩИЙ ОПЫТ: правая часть СЛУЧАЙНАЯ, оператор ТОТ ЖЕ.
     * `boff` — аффинный сдвиг внутри оператора, его подменять нельзя; подменяем
     * только правую часть системы. Если и на случайной метод сходится за шаг —
     * дело в операторе; если берёт много — значит быстрая сходимость была
     * свойством ЭТОЙ правой части, и цитировать ускорение по ней нельзя. */
    double *rhs = calloc((size_t)n, sizeof(double));
    if (rhs == NULL) {
      check(0, "память");
      free(ref);
      free(k.work);
      free(b);
      free(x);
      free(rr);
      free(rh);
      free(pv);
      free(vv);
      free(ss);
      free(tt);
      free(rhs);
      rig_free(&r);
      return;
    }
    if (cs >= 4) {
      uint32_t sd = 777u;
      for (int32_t i = 0; i < n; i++) {
        sd = sd * 1664525u + 1013904223u;
        rhs[i] = (double)(sd >> 8) / 8388608.0 - 1.0;
      }
    } else {
      memcpy(rhs, b, (size_t)n * sizeof(double));
    }
    /* BiCGStab по A·x = rhs, старт с нуля: r = rhs */
    memcpy(rr, rhs, (size_t)n * sizeof(double));
    memcpy(rh, rr, (size_t)n * sizeof(double));
    double rho = 1.0, alpha = 1.0, omega = 1.0, bn = kr_nrm(rhs, n);
    const char *why = "предел итераций";
    int it = 0, maxit = 400;
    for (; it < maxit; it++) {
      if (kr_nrm(rr, n) <= 1e-12 * bn) {
        why = "сошлось";
        break;
      }
      double rho1 = kr_dot(rh, rr, n);
      if (!(fabs(rho1) > 0.0)) {
        why = "обрыв: rho = 0";
        break;
      }
      if (it == 0) {
        memcpy(pv, rr, (size_t)n * sizeof(double));
      } else {
        double beta = (rho1 / rho) * (alpha / omega);
        for (int32_t i = 0; i < n; i++)
          pv[i] = rr[i] + beta * (pv[i] - omega * vv[i]);
      }
      kr_apply(&k, pv, b, vv);
      if (it == 0) {
        /* ДИАГНОЗ ОПЕРАТОРА, А НЕ ЗАДАЧИ: если `A·v ∥ v` и для СЛУЧАЙНОГО `v`,
         * значит оператор построен неверно и вырожден в кратный единичному. */
        double *rv = calloc((size_t)n, sizeof(double)), *av = calloc((size_t)n, sizeof(double));
        if (rv != NULL && av != NULL) {
          uint32_t sd = 12345u;
          for (int32_t i = 0; i < n; i++) {
            sd = sd * 1664525u + 1013904223u;
            rv[i] = (double)(sd >> 8) / 8388608.0 - 1.0;
          }
          kr_apply(&k, rv, b, av);
          k.npass--; /* диагностический проход в счёт не идёт */
          printf(
              "        [диагноз] на b: ‖A·b‖/‖b‖ = %.4f, cos = %.6f;  на СЛУЧАЙНОМ: cos = %.6f\n",
              kr_nrm(vv, n) / kr_nrm(b, n),
              kr_dot(b, vv, n) / (sqrt(kr_dot(b, b, n)) * sqrt(kr_dot(vv, vv, n))),
              kr_dot(rv, av, n) / (sqrt(kr_dot(rv, rv, n)) * sqrt(kr_dot(av, av, n))));
        }
        free(rv);
        free(av);
      }
      double den = kr_dot(rh, vv, n);
      if (!(fabs(den) > 0.0)) {
        why = "обрыв: (r0,v) = 0";
        break;
      }
      alpha = rho1 / den;
      for (int32_t i = 0; i < n; i++)
        ss[i] = rr[i] - alpha * vv[i];
      if (it == 0)
        printf("        [диагноз2] после 1 шага: ‖s‖/‖b‖ = %.3e, α = %.6f\n", kr_nrm(ss, n) / bn,
               alpha);
      if (kr_nrm(ss, n) <= 1e-12 * bn) {
        for (int32_t i = 0; i < n; i++)
          x[i] += alpha * pv[i];
        /* НАСТОЯЩАЯ невязка на этой ветке есть `s`, а не старое `r`: шаг сделан
         * половинный, и `r` ещё не обновлена. Печатать старое `r` — врать. */
        memcpy(rr, ss, (size_t)n * sizeof(double));
        why = "сошлось на половинном шаге";
        break;
      }
      kr_apply(&k, ss, b, tt);
      double tsq = kr_dot(tt, tt, n);
      if (!(tsq > 0.0)) {
        why = "обрыв: (t,t) = 0";
        break;
      }
      omega = kr_dot(tt, ss, n) / tsq;
      for (int32_t i = 0; i < n; i++) {
        x[i] += alpha * pv[i] + omega * ss[i];
        rr[i] = ss[i] - omega * tt[i];
      }
      rho = rho1;
      if (!(fabs(omega) > 0.0)) {
        why = "обрыв: omega = 0";
        break;
      }
    }
    /* сверка с эталоном */
    double worst = 0.0, scale = kr_nrm(ref, n);
    if (!(scale > 0.0)) scale = 1.0;
    for (int32_t i = 0; i < n; i++) {
      double e = fabs(x[i] - ref[i]) / scale;
      if (e > worst) worst = e;
    }
    if (cs >= 4) {
      printf("  [КРЫЛОВ, СЛУЧАЙНАЯ правая часть] σ_t·h = %.2f: ряд Неймана %4d, BiCGStab %3ld "
             "проходов; невязка %.2e; останов: %s\n",
             sigh[cs], passes_ref, k.npass, kr_nrm(rr, n) / bn, why);
      check(k.npass > 4, "РЕШАЮЩИЙ ОПЫТ: на общей правой части метод обязан взять БОЛЬШЕ шагов");
      check(k.npass < 40, "и всё же ЦЕНА ОСТАЁТСЯ ОГРАНИЧЕННОЙ при любой толщине");
    } else
      printf("  [КРЫЛОВ] σ_t·h = %.2f, альбедо %.0f: ряд Неймана %4d проходов, BiCGStab %3ld "
             "(×%.1f); расхождение %.2e\n",
             sigh[cs], alb, passes_ref, k.npass, (double)passes_ref / (double)k.npass, worst);
    if (cs < 3) {
      check(worst < 1e-9, "КРЫЛОВ: ответ тот же, что у итерации по рассеянию");
      check(k.npass < passes_ref, "КРЫЛОВ: проходов МЕНЬШЕ, чем у ряда Неймана");
    } else if (cs == 3) {
      /* НЕГАТИВНЫЙ КОНТРОЛЬ: без рассеяния T = 0, и обе схемы обязаны решить
       * задачу за ОДИН проход. Если Крылов возьмёт больше — ошибка в постановке
       * оператора, а не в скорости. */
      check(k.npass <= 3, "НК: при σ_s = 0 оператор тождественный, проходов должно быть ~1");
      check(worst < 1e-12, "НК: и ответ совпадает точно");
    }
    free(k.work);
    free(b);
    free(rhs);
    free(x);
    free(rr);
    free(rh);
    free(pv);
    free(vv);
    free(ss);
    free(tt);
    free(ref);
    rig_free(&r);
  }
}

/* --------------- ЗАМКНУТ ЛИ НАБОР ОРДИНАТ ОТНОСИТЕЛЬНО ЗЕРКАЛЬНОГО ОТРАЖЕНИЯ */

/* ОТ ЭТОГО ЗАВИСИТ, МОЖЕТ ЛИ ЗЕРКАЛО ИДТИ ЧЕРЕЗ РАЗВЁРТКУ ВОВСЕ.
 *
 * К5 записана как «зеркала через ординаты НЕ идут», и для КРИВОГО зеркала это
 * верно безоговорочно: нормаль меняется непрерывно вдоль поверхности, значит
 * отражений континуум, и ни один конечный набор под них не замкнут.
 *
 * Но для ПЛОСКОГО зеркала отражение есть ИЗОМЕТРИЯ пространства направлений, и
 * если набор переходит в себя, то зеркальная связь — ПЕРЕСТАНОВКА ординат.
 * Точная, разреженная, без интерполяции (то есть К3 не нарушается) и без
 * всякого прослеживания. Тогда бесконечная серия отражений между зеркалами
 * получается из РЕШЕНИЯ системы, а не из усечённого ряда с пределом отскоков.
 *
 * Проверяется здесь то, что нужно на деле: замкнутость относительно отражений
 * от стенок куба, то есть от плоскостей, перпендикулярных осям. Для них
 * отражение есть смена знака ОДНОЙ компоненты.
 *
 * НЕГАТИВНЫЙ КОНТРОЛЬ ЗАЛОЖЕН В ТОТ ЖЕ ТЕСТ: наклонная плоскость (нормаль по
 * диагонали) обязана замкнутость СЛОМАТЬ — иначе проверка ничего не значит,
 * потому что «замкнут под что угодно» невозможно. */
static void t_mirror_closure(void) {
  tr3_dirs d;
  if (tr3_dirs_product(&d, 2, 2)) {
    check(0, "набор направлений");
    return;
  }
  for (int ax = 0; ax < 3; ax++) {
    double worst = 0.0, wworst = 0.0;
    int missing = 0;
    for (int i = 0; i < d.n; i++) {
      double r[3] = {d.ox[i], d.oy[i], d.oz[i]};
      r[ax] = -r[ax]; /* отражение от плоскости, перпендикулярной оси ax */
      double best = 1e300;
      int bj = -1;
      for (int j = 0; j < d.n; j++) {
        double e = fabs(d.ox[j] - r[0]) + fabs(d.oy[j] - r[1]) + fabs(d.oz[j] - r[2]);
        if (e < best) {
          best = e;
          bj = j;
        }
      }
      if (best > worst) worst = best;
      if (bj >= 0 && fabs(d.w[bj] - d.w[i]) > wworst) wworst = fabs(d.w[bj] - d.w[i]);
      if (best > 1e-14) missing++;
    }
    printf(
        "  [ЗЕРКАЛО-ОСЬ %d] max промах отражённой ординаты %.3e, вес %.3e, вне набора %d из %d\n",
        ax, worst, wworst, missing, d.n);
    check(worst < 1e-14, "набор ЗАМКНУТ относительно отражения от стенки: связь есть перестановка");
    check(wworst < 1e-14, "и веса отражённых ординат совпадают");
  }
  /* ГРУППА СИММЕТРИИ НАБОРА ОКАЗАЛАСЬ ШИРЕ ОЖИДАЕМОГО, И ЭТО ЗАМЕР, А НЕ ДОГАДКА.
   * Первая редакция брала негативным контролем плоскость (1,1,0) — и он НЕ
   * ВЫСТРЕЛИЛ: набор под неё тоже замкнут, потому что азимутальные узлы
   * симметричны относительно 45°, и вся конструкция инвариантна относительно
   * ОКТАЭДРИЧЕСКОЙ группы, то есть трёх осевых плоскостей И ШЕСТИ ДИАГОНАЛЬНЫХ.
   * Девять точно представимых зеркальных плоскостей вместо трёх — это подарок,
   * но проверять надо на нормали, которая симметрией НЕ является. */
  {
    double n0[3] = {0.6, 0.48, 0.64}, nn = 0.0;
    for (int a = 0; a < 3; a++)
      nn += n0[a] * n0[a];
    nn = sqrt(nn);
    double n[3] = {n0[0] / nn, n0[1] / nn, n0[2] / nn}, worst = 0.0;
    int missing = 0;
    for (int i = 0; i < d.n; i++) {
      double o[3] = {d.ox[i], d.oy[i], d.oz[i]};
      double dn = o[0] * n[0] + o[1] * n[1] + o[2] * n[2], r[3];
      for (int a = 0; a < 3; a++)
        r[a] = o[a] - 2.0 * dn * n[a];
      double best = 1e300;
      for (int j = 0; j < d.n; j++) {
        double e = fabs(d.ox[j] - r[0]) + fabs(d.oy[j] - r[1]) + fabs(d.oz[j] - r[2]);
        if (e < best) best = e;
      }
      if (best > worst) worst = best;
      if (best > 1e-14) missing++;
    }
    printf("    [НК несимметричная нормаль (0.6,0.48,0.64)] max промах %.3e, вне набора %d из %d\n",
           worst, missing, d.n);
    check(missing > 0, "НК: под НАКЛОННОЕ зеркало набор замкнут быть НЕ должен");
  }
  tr3_dirs_free(&d);
}

/* ------------------------------------- ФОРМАТ ВЫХОДА: PFM, круговой прогон */

/* Писалки картинок врут МОЛЧА и ровно в трёх местах: порядок строк, порядок
 * байт и чередование каналов. Ни одно из трёх не видно на картинке, если
 * смотреть её тем же кодом, который писал. Поэтому файл читается ОБРАТНО
 * сырыми байтами и сверяется с тем, что клали. */
static void t_pfm(void) {
  const int W = 7, H = 5;
  double buf[35];
  for (int y = 0; y < H; y++)
    for (int x = 0; x < W; x++)
      buf[y * W + x] = 0.125 * (double)(y * W + x) + 0.0625; /* точно представимо во float */
  const char *path = "img/_pfm_roundtrip.pfm";
  check(hz_pfm_write(path, buf, buf, buf, W, H) == 0, "PFM записан");
  FILE *f = fopen(path, "rb");
  if (f == NULL) {
    check(0, "PFM открылся на чтение");
    return;
  }
  int w2 = 0, h2 = 0;
  double sc = 0.0;
  char magic[3] = {0, 0, 0};
  int ok_hdr = fscanf(f, "%2s %d %d %lf", magic, &w2, &h2, &sc) == 4;
  fgetc(f); /* один перевод строки после масштаба — часть формата */
  printf("  [PFM] заголовок: «%s» %dx%d масштаб %.1f\n", magic, w2, h2, sc);
  check(ok_hdr && magic[0] == 'P' && magic[1] == 'F', "PFM: три канала (PF), а не Pf");
  check(w2 == W && h2 == H, "PFM: размеры");
  check(sc < 0.0, "PFM: масштаб отрицателен — little-endian, как на этой машине");
  double worst = 0.0;
  int bad_order = 0;
  for (int r = 0; r < H; r++) {
    /* СТРОКИ В ФАЙЛЕ ИДУТ СНИЗУ ВВЕРХ: строка r файла есть строка H−1−r буфера */
    int y = H - 1 - r;
    for (int x = 0; x < W; x++) {
      float px[3];
      if (fread(px, sizeof(float), 3, f) != 3) {
        bad_order = 1;
        break;
      }
      float want = (float)buf[y * W + x];
      for (int c = 0; c < 3; c++) {
        double e = fabs((double)px[c] - (double)want);
        if (e > worst) worst = e;
      }
    }
  }
  long extra = 0;
  while (fgetc(f) != EOF)
    extra++;
  fclose(f);
  printf("  [PFM] круговой прогон: max |прочитано − записано| = %.3e, лишних байт %ld\n", worst,
         extra);
  check(!bad_order, "PFM: файл не оборван");
  check(!(worst > 0.0), "PFM: значения совпали ПОБИТОВО (порядок строк и байт верны)");
  check(extra == 0, "PFM: ни одного лишнего байта");
  remove(path); /* временный файл теста, а не картинка на посмотреть */
}

/* ------------------------------------ К102: линейность сравнения с Крыловым ---
 *
 * РЕГРЕССИЯ НАХОДКИ К102: стенд ЭТАП B сравнивал развёртку (с НЕЛИНЕЙНЫМИ
 * §735-проекциями maxp) с линейным оператором Крылова — ответы разошлись на
 * сотни, прибор стенда показывал невязку 14–114. Замок: ЛИНЕЙНАЯ развёртка
 * (limiter=0, maxp_off=1) обязана удовлетворять оператору Крылова, Крылов
 * обязан дать тот же ответ, и возмущённое поле обязано невязку ПОДНЯТЬ
 * (иначе прибор ничего не мерит). */
static void t_k102(void) {
  const hz_frame fr = {{0, 0, 0}, {1.0, 1.0, 2.0}};
  for (int graded = 0; graded < 2; graded++) {
    rig r;
    if (rig_init(&r, 3, graded, 2, &fr)) {
      check(0, "оснастка");
      rig_free(&r);
      return;
    }
    for (int32_t c = 0; c < r.m.ncell; c++) {
      r.sig_t[c] = 0.5;
      r.sig_s[c] = 0.5;
    }
    tr3_problem p = {.m = &r.m,
                     .d = &r.d,
                     .sig_t = r.sig_t,
                     .sig_s = r.sig_s,
                     .binc0 = 1.0,
                     .binc = {0.13, -0.07, 0.05},
                     .binx0 = {4, 4, 4},
                     .limiter = 0,
                     .maxp_off = 1 /* ТОТ ЖЕ линейный оператор, что у Крылова */};
    double *phi_lin = calloc((size_t)r.m.ncell * 4, sizeof(double));
    double *phi_kry = calloc((size_t)r.m.ncell * 4, sizeof(double));
    double *bk = calloc((size_t)r.m.nf * 4, sizeof(double));
    if (!phi_lin || !phi_kry || !bk) {
      check(0, "память");
      free(phi_lin);
      free(phi_kry);
      free(bk);
      rig_free(&r);
      return;
    }
    tr3_stats st;
    int rcA = tr3_sweep_solve(&p, 4000, 1e-12, phi_lin, &st);
    double relA = -1.0;
    int rrA = tr3_krylov_residual(&p, phi_lin, st.bout, st.sout, &relA);
    printf("  [К102 %s] развёртка: код %d, итераций %d; прибор оператора: %.3e\n",
           graded ? "градуир." : "равном.", rcA, st.iters, relA);
    check(rcA == 0 && rrA == 0, "К102: оба решателя прошли");
    check(relA <= 1e-8, "К102: линейная развёртка удовлетворяет оператору Крылова");
    /* НЕГАТИВНЫЙ КОНТРОЛЬ (предсказан провал): возмущение одной ячейки
     * обязано поднять невязку — иначе прибор не чувствителен */
    phi_lin[0] += 0.5;
    double relP = -1.0;
    tr3_krylov_residual(&p, phi_lin, st.bout, st.sout, &relP);
    phi_lin[0] -= 0.5;
    check(relP > relA * 10.0 && relP > 1e-6,
          "К102 НК: возмущённое поле ОБЯЗАНО поднять невязку оператора");
    /* Крылов на ТОЙ ЖЕ линейной задаче: ответы обязаны совпасть */
    tr3_kstats kst;
    int rcB = tr3_krylov_solve(&p, 200, 1e-11, phi_kry, bk, NULL, &kst);
    check(rcB == 0, "К102: Крылов прошёл");
    double dm = 0.0, pm = 0.0;
    for (int32_t i = 0; i < r.m.ncell * 4; i++) {
      double e = fabs(phi_lin[i] - phi_kry[i]);
      if (e > dm) dm = e;
      if (fabs(phi_kry[i]) > pm) pm = fabs(phi_kry[i]);
    }
    printf("    Крылов: проходов %ld, итераций %d; max|Δφ| развёртка/Крылов %.3e при "
           "|φ|max %.3e\n",
           kst.npass, kst.iters, dm, pm);
    check(dm <= 1e-6, "К102: Крылов даёт ТОТ ЖЕ ответ, что линейная развёртка");
    free(st.bout);
    free(st.sout);
    free(st.eirr);
    free(phi_lin);
    free(phi_kry);
    free(bk);
    rig_free(&r);
  }
}

/* -------------------------------------------------------------------------- */

int main(void) {
  printf("=== ПРЕДСКАЗАНИЯ (до единого результата) ===\n");
  printf("  A1 Σw = 4π, первый момент 0, второй (4π/3)·I — все ≤1e-13;\n");
  printf("     ни одно направление не параллельно оси (иначе грань с ω·n = 0)\n");
  printf("  A2 моменты по многоугольнику совпадают с замкнутой формой\n");
  printf("  B  К13 в 3D: сумма площадей подграней ТОЧНО равна площади стороны\n");
  printf("     (все координаты целые, поэтому равенство точное, а не с допуском)\n");
  printf("  C  печь: радианс однороден и равен влёту, ≤1e-12, наклоны нулевые,\n");
  printf("     на равномерной И на градуированной сетке\n");
  printf("  C  баланс при НЕОДНОРОДНОЙ подсветке и альбедо 1: невязка ≤1e-10\n");
  printf("=== НЕГАТИВНЫЙ КОНТРОЛЬ, с ПРЕДСКАЗАННЫМ провалом ===\n");
  printf("  «взять одного соседа» вместо всех -> тождество площадей ОБЯЗАНО\n");
  printf("  сломаться, и только на градуированной сетке\n");
  printf("=== РЕЗУЛЬТАТЫ ===\n");

  t_dirs();
  t_face_moments();
  t_mesh();
  t_elem_positivity();
  t_union();
  t_k102();
  t_k50();
  t_furnace();
  t_linear();
  t_linear_field();
  t_mirror_wall();
  t_balance();
  t_cavity();
  t_cut();
  t_krylov();
  t_mirror_closure();
  t_pfm();
  t_render();

  printf("%s: %d/%d\n", g_fail ? "FAILURES" : "ok", g_total - g_fail, g_total);
  return g_fail != 0;
}
