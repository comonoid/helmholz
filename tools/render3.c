/* render3.c — ПОЛНОЦЕННЫЙ РЕНДЕР линии переноса: сфера в освещённой комнате,
 * глобальное освещение развёрткой по ординатам.
 *
 * PLAN_TRANSPORT.md, все три прохода вместе:
 *   развёртка   — полное решение уравнения переноса с многократными отражениями;
 *   поверхности — граничное условие, исходящий радианс хранится DG1 по положению;
 *   сбор        — луч на пиксель до первого попадания, чтение хранимого В ТОЧКЕ.
 *
 * СИЛУЭТ СФЕРЫ ТОЧЕН: попадание считается по АНАЛИТИЧЕСКОМУ ПРИМИТИВУ (К15),
 * фасеты работают только пространственным индексом и режут ячейку. Фасетной
 * огранки на картинке поэтому нет вовсе — она была бы, если бы луч спрашивал
 * фасеты.
 *
 * ДВОЙНОГО УЧЁТА НЕТ ПО ПОСТРОЕНИЮ (К9/К17): источник здесь — ИЗЛУЧАЮЩАЯ
 * ПОВЕРХНОСТЬ, развёртка считает всё, а сбор только ЧИТАЕТ. Прямой свет
 * отдельным проходом не добавляется, потому что он уже внутри.
 */

#include "cut/surf.h"
#include "image.h"
#include "octree.h"
#include "transport/cam3.h"
#include "transport/cut3.h"
#include "transport/dirs3.h"
#include "transport/mesh3.h"
#include "transport/ray3.h"
#include "transport/sweep3.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

#define LOG2N 4
#define NC (1 << LOG2N)

int main(int argc, char **argv) {
  int W = argc > 1 ? atoi(argv[1]) : 640;
  int H = argc > 2 ? atoi(argv[2]) : 480;
  int nmu = argc > 3 ? atoi(argv[3]) : 2;
  const char *out = argc > 4 ? argv[4] : "img/room_sphere.ppm";

  hz_frame fr = {{0, 0, 0}, {1, 1, 1}};
  hz_octree t;
  if (hz_oct_init(&t, LOG2N, 0.0)) return 1;
  /* РАВНОМЕРНОЕ дробление: условие 1:1 у разрезанных ячеек (cut3) требует, чтобы
   * грань сетки совпадала с гранью коробки. Градуированную сетку у поверхности
   * пришлось бы ещё и обрезать прямоугольником — это отдельная работа. */
  for (int x = 0; x < NC; x++)
    for (int y = 0; y < NC; y++)
      for (int z = 0; z < NC; z++) {
        int lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        hz_oct_set_box(&t, lo, hi, 1.0);
      }

  /* сфера как АНАЛИТИЧЕСКИЙ примитив плюс её фасеты для разреза */
  hz_surftab stab;
  hz_facettab ftab;
  hz_cutmap cmap;
  if (hz_surftab_init(&stab) || hz_facettab_init(&ftab) || hz_cutmap_init(&cmap)) return 1;
  /* ДВА ТЕЛА: тени одного на другом и переотражения между ними — то, ради чего
   * развёртка и нужна. Материал в модели разреза есть ПЕРЕСЕЧЕНИЕ полуплоскостей,
   * то есть ОДНО выпуклое тело на ячейку; два тела в одной ячейке дали бы
   * пересечение вместо объединения. Тела разнесены, и это ПРОВЕРЯЕТСЯ. */
  const int NB = 2;
  double sc[2][3] = {{5.6, 9.2, 4.2}, {10.9, 7.4, 3.1}};
  double sr[2] = {3.0, 1.9};
  int32_t f0[2], nfac[2];
  for (int b = 0; b < NB; b++) {
    hz_surf sp = {HZ_SURF_SPHERE, {sc[b][0], sc[b][1], sc[b][2], sr[b], 0, 0, 0}, 1, 0};
    int32_t si = hz_surftab_add(&stab, &sp);
    f0[b] = 0;
    nfac[b] = hz_surf_facet_sphere(&ftab, &fr, sc[b], sr[b], 2, HZ_FIT_MEAN_SAGITTA, si, &f0[b]);
    if (nfac[b] <= 0) return 1;
  }

  /* боковая таблица: ключи ОБЯЗАНЫ идти по возрастанию (Г45) */
  typedef struct {
    int32_t cell, f[HZ_P3_MAXH], nf;
  } rec_t;
  rec_t *recs = calloc((size_t)NC * NC * NC, sizeof(rec_t));
  if (recs == NULL) return 1;
  int nrec = 0, nboth = 0;
  for (int x = 0; x < NC; x++)
    for (int y = 0; y < NC; y++)
      for (int z = 0; z < NC; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        int used = -1, ns = 0;
        for (int b = 0; b < NB; b++) {
          int32_t s2[HZ_P3_MAXH];
          int k2 = hz_facets_for_box(&ftab, f0[b], nfac[b], lo, hi, s2, HZ_P3_MAXH);
          if (k2 <= 0) continue;
          if (used >= 0) {
            nboth++;
            continue;
          }
          used = b;
          ns = k2;
          for (int j = 0; j < k2; j++)
            sel[j] = s2[j];
        }
        if (ns <= 0) continue;
        int32_t ni = hz_oct_leaf(&t, x, y, z);
        recs[nrec].cell = ni;
        recs[nrec].nf = ns;
        for (int j = 0; j < ns; j++)
          recs[nrec].f[j] = sel[j];
        nrec++;
      }
  for (int i = 1; i < nrec; i++) { /* сортировка вставками по ключу */
    rec_t tmp = recs[i];
    int j = i - 1;
    while (j >= 0 && recs[j].cell > tmp.cell) {
      recs[j + 1] = recs[j];
      j--;
    }
    recs[j + 1] = tmp;
  }
  for (int i = 0; i < nrec; i++)
    if (hz_cutmap_add(&cmap, recs[i].cell, recs[i].f, recs[i].nf) != 0) return 1;
  free(recs);

  tr3_mesh mesh;
  if (tr3_mesh_build(&mesh, &t, &fr)) return 1;
  /* Г38: маска полных ячеек строится ОТДЕЛЬНО — в боковой таблице их нет */
  uint8_t *solid = calloc((size_t)mesh.ncell, 1);
  if (solid == NULL) return 1;
  for (int x = 0; x < NC; x++)
    for (int y = 0; y < NC; y++)
      for (int z = 0; z < NC; z++) {
        int32_t lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        int32_t sel[HZ_P3_MAXH];
        for (int b = 0; b < NB; b++)
          if (hz_facets_for_box(&ftab, f0[b], nfac[b], lo, hi, sel, HZ_P3_MAXH) == 0)
            solid[mesh.cellof[hz_oct_leaf(&t, x, y, z)]] = 1;
      }
  tr3_cut cut;
  if (tr3_cut_build(&cut, &mesh, &ftab, &cmap, solid)) return 1;
  tr3_dirs dirs;
  if (tr3_dirs_product(&dirs, nmu, nmu)) return 1;

  double *sig_t = calloc((size_t)mesh.ncell, sizeof(double));
  double *sig_s = calloc((size_t)mesh.ncell, sizeof(double));
  if (sig_t == NULL || sig_s == NULL) return 1;
  for (int32_t c = 0; c < mesh.ncell; c++) {
    sig_t[c] = 0.015; /* лёгкая дымка: среда участвует, но не глушит */
    sig_s[c] = 0.012; /* альбедо 0.8, а не 1: при 1 ряд по рассеянию сходится как
                       * ρⁿ и упирается ровно в то, ради чего С1 и заводил DSA */
    ;
  }
  /* комната: потолок светит, пол и стены серые, одна стена красноватая по
   * яркости (цвета нет — считается один спектральный канал) */
  double wr[6] = {0.72, 0.35, 0.72, 0.72, 0.65, 0.05};
  double we[6] = {0, 0, 0, 0, 0, 6.0};
  double *frho = calloc((size_t)ftab.n, sizeof(double));
  double *femit = calloc((size_t)ftab.n, sizeof(double));
  if (frho == NULL || femit == NULL) return 1;
  for (int32_t i = 0; i < ftab.n; i++)
    frho[i] = 0.78; /* светлая диффузная сфера */

  tr3_problem prob = {.m = &mesh,
                      .d = &dirs,
                      .cut = &cut,
                      .facet_rho = frho,
                      .facet_emit = femit,
                      .nfacet = ftab.n,
                      .sig_t = sig_t,
                      .sig_s = sig_s,
                      .wall_rho = wr,
                      .wall_emit = we,
                      .limiter = 1};
  double *phi = calloc((size_t)mesh.ncell * 4, sizeof(double));
  if (phi == NULL) return 1;
  tr3_stats st;
  printf("ячеек %d, граней %d, поверхностных элементов %d, направлений %d, 1:1 нарушений %d, "
         "ячеек с ДВУМЯ телами %d (обязано быть 0)\n",
         mesh.ncell, mesh.nf, cut.nse, dirs.n, cut.nbad, nboth);
  double t0 = now();
  int rc = tr3_sweep_solve(&prob, 4000, 1e-9, phi, &st);
  double t_sweep = now() - t0;
  printf("развёртка (ХОЛОДНЫЙ СТАРТ): код %d, итераций %d, невязка %.2e, срезок %d, %.2f с "
         "(%.2f мс на итерацию)\n",
         rc, st.iters, st.resid, st.nclip, t_sweep, 1e3 * t_sweep / (double)(st.iters + 1));
  printf("баланс: втекло %.4f, вытекло %.4f, поглощено %.4f, невязка %.2e\n", st.pin, st.pout,
         st.pabs, st.balance);
  if (rc != 0) return 1;

  /* --- сбор по пикселю --- */
  tr3_scene scn = {.tree = &t, .fr = fr, .sigma = NULL, .st = &stab, .ft = &ftab, .cm = &cmap};
  tr3_camera cam;
  double eye[3] = {8.0, 0.6, 7.2}, at[3] = {8.2, 9.5, 3.6}, up[3] = {0, 0, 1};
  if (tr3_camera_look(&cam, eye, at, up, 1.3, W, H)) return 1;
  double *buf = calloc((size_t)W * (size_t)H, sizeof(double));
  if (buf == NULL) return 1;

  /* ИНДЕКС ГРАНИЧНЫХ ГРАНЕЙ ПО СЕТКЕ. Первая редакция искала грань ЛИНЕЙНЫМ
   * перебором всех 13056 граней на КАЖДЫЙ пиксель — это 2.7e10 сравнений на
   * кадр, и сбор стоил 3930 нс на луч при разумных 100-200. Грани лежат на
   * ЦЕЛОЧИСЛЕННОЙ сетке, поэтому индекс прямой: (стенка, u, v) -> грань. */
  int32_t *wallidx = calloc((size_t)6 * NC * NC, sizeof(int32_t));
  if (wallidx == NULL) return 1;
  for (int32_t i = 0; i < 6 * NC * NC; i++)
    wallidx[i] = -1;
  for (int32_t f = 0; f < mesh.nf; f++) {
    const tr3_face *ff = &mesh.f[f];
    if (ff->cb >= 0) continue;
    int wl = (int)(~ff->cb);
    for (int32_t a = ff->lo[0]; a < ff->hi[0]; a++)
      for (int32_t b = ff->lo[1]; b < ff->hi[1]; b++)
        wallidx[((int32_t)wl * NC + a) * NC + b] = f;
  }

  double t1 = now();
  for (int py = 0; py < H; py++)
    for (int px = 0; px < W; px++) {
      double o[3], d[3];
      tr3_camera_ray(&cam, px, py, o, d);
      tr3_hit h;
      tr3_march(&scn, o, d, -1.0, &h);
      double val = 0.0;
      int skip = 0;
      if (h.hit) {
        /* ИНДЕКС ЯЧЕЙКИ, А НЕ УЗЛА: марш возвращает индекс узла октодерева, а
         * таблицы сетки и разреза живут по КОМПАКТНОМУ индексу ячейки. Перевод
         * делает cellof, и без него поиск не находил ничего — сфера выходила
         * чёрной при идеально круглом силуэте. */
        int32_t mc = mesh.cellof[h.cell];
        if (mc < 0) skip = 1;
        /* ЭЛЕМЕНТ ИЩЕТСЯ ПО НОРМАЛИ, А НЕ ПО НОМЕРУ ФАСЕТА. В примитивном пути
         * (К15) номер фасета не проставляется вовсе — попадание считает ПРИМИТИВ,
         * и фасет к ответу отношения не имеет. Ближайший по нормали элемент и
         * есть тот кусок поверхности, в который попал луч. */
        double best = -2.0;
        int32_t bi = -1;
        if (skip)
          bi = -1;
        else
          for (int32_t k = cut.sestart[mc]; k < cut.sestart[mc + 1]; k++) {
            int32_t e = cut.selist[k];
            double dp = cut.se[e].n[0] * h.n[0] + cut.se[e].n[1] * h.n[1] + cut.se[e].n[2] * h.n[2];
            if (dp > best) {
              best = dp;
              bi = e;
            }
          }
        if (bi >= 0) {
          int32_t e = bi;
          double sz = (double)mesh.csize[mc], b[4] = {1, 0, 0, 0};
          for (int a = 0; a < 3; a++) {
            double xu = (h.p[a] - fr.o[a]) / fr.u[a];
            b[a + 1] = (xu - ((double)mesh.clo[mc][a] + 0.5 * sz)) / sz;
          }
          for (int i = 0; i < 4; i++)
            val += st.sout[e * 4 + i] * b[i];
        }
      } else {
        /* СТЕНКА: выход из куба, затем грань, накрывающая точку */
        double tex = 1e300;
        int wall = -1;
        for (int a = 0; a < 3; a++) {
          if (!(fabs(d[a]) > 0.0)) continue;
          double lim = d[a] > 0.0 ? (double)NC : 0.0;
          double tt = (lim - o[a]) / d[a];
          if (tt > 0.0 && tt < tex) {
            tex = tt;
            wall = 2 * a + (d[a] > 0.0 ? 1 : 0);
          }
        }
        if (wall >= 0) {
          double hp[3];
          for (int a = 0; a < 3; a++)
            hp[a] = o[a] + tex * d[a];
          int axis = wall / 2, u = (axis + 1) % 3, v = (axis + 2) % 3;
          if (u > v) {
            int s2 = u;
            u = v;
            v = s2;
          }
          int32_t iu = (int32_t)floor(hp[u]), iv = (int32_t)floor(hp[v]);
          if (iu < 0) iu = 0;
          if (iu >= NC) iu = NC - 1;
          if (iv < 0) iv = 0;
          if (iv >= NC) iv = NC - 1;
          int32_t f = wallidx[((int32_t)wall * NC + iu) * NC + iv];
          if (f >= 0) {
            const tr3_face *ff = &mesh.f[f];
            int32_t cc = ff->ca;
            double sz = (double)mesh.csize[cc], b[4] = {1, 0, 0, 0};
            for (int a = 0; a < 3; a++) {
              double xu = (hp[a] - fr.o[a]) / fr.u[a];
              b[a + 1] = (xu - ((double)mesh.clo[cc][a] + 0.5 * sz)) / sz;
            }
            for (int i = 0; i < 4; i++)
              val += st.bout[f * 4 + i] * b[i];
          }
        }
      }
      buf[(size_t)py * (size_t)W + (size_t)px] = val > 0.0 ? val : 0.0;
    }

  double t_gather = now() - t1;
  /* РАЗДЕЛЬНЫЙ ДОКЛАД ХОЛОДНОГО СТАРТА И КАДРА — требование CLAUDE.md. Поле от
   * камеры не зависит, поэтому поворот камеры стоит ТОЛЬКО сбора. */
  printf("СБОР ПО ПИКСЕЛЮ (кадр): %.3f с на %dx%d = %.2f млн лучей, %.0f нс на луч\n", t_gather, W,
         H, 1e-6 * (double)W * (double)H, 1e9 * t_gather / ((double)W * (double)H));
  if (hz_ppm_write(out, buf, W, H) != 0) {
    fprintf(stderr, "не записалось: %s\n", out);
    return 1;
  }
  printf("картинка: %s (%dx%d)\n", out, W, H);
  free(wallidx);
  free(buf);
  free(phi);
  free(st.bout);
  free(st.sout);
  free(sig_t);
  free(sig_s);
  free(frho);
  free(femit);
  tr3_dirs_free(&dirs);
  tr3_cut_free(&cut);
  tr3_mesh_free(&mesh);
  hz_cutmap_free(&cmap);
  hz_facettab_free(&ftab);
  hz_surftab_free(&stab);
  hz_oct_free(&t);
  return 0;
}
