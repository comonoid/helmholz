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
  double d[3];
  for (int a = 0; a < 3; a++)
    d[a] = p[a] - c->eye[a];
  v[0] = d[0] * c->rt[0] + d[1] * c->rt[1] + d[2] * c->rt[2];
  v[1] = d[0] * c->up[0] + d[1] * c->up[1] + d[2] * c->up[2];
  v[2] = d[0] * c->fwd[0] + d[1] * c->fwd[1] + d[2] * c->fwd[2];
}

int main(int argc, char **argv) {
  int city = 0, nmu = 2, nphi = 4, nb = 21, layout = HZ_LAYOUT_RUNS, imgw = 960, imgh = 540;
  double h = 0.5, rho = 0.5, sky = HZ_CFG_SKY_LE;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "city") == 0) city = 1;
    if (strncmp(argv[i], "nmu=", 4) == 0) nmu = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "nphi=", 5) == 0) nphi = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "nb=", 3) == 0) nb = (int)strtol(argv[i] + 3, NULL, 10);
    if (strncmp(argv[i], "h=", 2) == 0) h = strtod(argv[i] + 2, NULL);
    if (strncmp(argv[i], "rho=", 4) == 0) rho = strtod(argv[i] + 4, NULL);
    if (strncmp(argv[i], "sky=", 4) == 0) sky = strtod(argv[i] + 4, NULL);
    if (strcmp(argv[i], "list") == 0) layout = HZ_LAYOUT_LIST;
    if (strncmp(argv[i], "w=", 2) == 0) imgw = (int)strtol(argv[i] + 2, NULL, 10);
    if (strncmp(argv[i], "ih=", 3) == 0) imgh = (int)strtol(argv[i] + 3, NULL, 10);
  }

  double t0 = now_s();
  hz_objmesh m;
  if (hz_obj_load(&m, city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ,
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
  if (hz_seg_planar_cap(&sg, &m, dseg, city ? 0.0 : 0.5) != 0) return 1;
  if (sg.nseg >= m.nt) {
    fprintf(stderr, "сегментация выродилась: участков %d при %d треугольниках\n", sg.nseg, m.nt);
    return 1;
  }
  hz_polyset ps;
  if (hz_poly_build(&ps, &m, &sg) != 0) return 1;
  double t_geo = now_s() - t0;
  printf("== СЦЕНА: %s, треугольников %d, полигонов %d; разбор %.1f с, геометрия %.1f с\n",
         city ? "ГОРОД" : "зал", m.nt, ps.np, t_load, t_geo);

  hz_ptrans tr;
  if (hz_ptrans_init(&tr, &ps, &m) != 0) return 1;
  /* Альбедо ставится ОДИНАКОВЫМ с эталоном (§120: `ρ = 0.5` у всех), иначе
   * сравнивать потоки нельзя. Собственного излучения нет: источник — небо. */
  for (int32_t k = 0; k < ps.np; k++) {
    tr.rho[k] = rho;
    tr.Le[k] = 0.0;
  }

  tr3_dirs d;
  if (tr3_dirs_product(&d, nmu, nphi) != 0) return 1;
  printf("== НАПРАВЛЕНИЙ: %d (nmu %d × nphi %d); шаг растра %.2f м; отскоков %d; ρ %.2f; "
         "небо %.2f\n",
         d.n, nmu, nphi, h, nb, rho, sky);

  printf("  отскок   Σ B·A, Вт/ср      dE      фрагментов    растр,с   редукция,с   всего,с\n");
  double ttot = 0.0;
  for (int b = 0; b < nb; b++) {
    hz_pstats st;
    memset(&st, 0, sizeof st);
    double tb = now_s();
    if (hz_psweep_bounce(&tr, &d, h, 0, sky, layout, &st) != 0) {
      fprintf(stderr, "отказ отскока %d\n", b);
      return 1;
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
  printf("== ВСЕГО НА ПЕРЕНОС: %.1f с на %d отскоков (%.2f с на отскок, %.3f с на направление)\n",
         ttot, nb, ttot / (double)(nb > 0 ? nb : 1),
         ttot / (double)((nb > 0 ? nb : 1) * (d.n > 0 ? d.n : 1)));
  printf("== ЭТАЛОН ХРАНИМЫМ ОПЕРАТОРОМ (§120, тот же город, ρ = 0.5, небо 1.0): "
         "поток 2.398505e+05 Вт/ср, сборка 3102.6 с\n");

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
        double sx[3], sy[3], iw[3];
        double wx[3][3];
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

  tr3_dirs_free(&d);
  hz_ptrans_free(&tr);
  hz_poly_free(&ps);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
