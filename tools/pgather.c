/* pgather — СБОР ПО ПИКСЕЛЮ НА НОВОМ НОСИТЕЛЕ (§847): первая картинка
 * линии пирамиды + объёмного фронта.
 *
 * Камера — ОПЕРАТОР над решённым полем: поле считается свипом mode=2
 * (§845/§846, build=1), затем на каждый пиксель пускается луч из глаза,
 * ищется БЛИЖАЙШЕЕ пересечение с куском, яркость пикселя = L_out куска:
 *
 *     L_out = Le + rho·E_куска/(2π)     (ламберт, А1471 — двусторонняя)
 *
 * фон — 0 (тёмное окружение, §842). Пересечение — перебором ВСЕХ кусков
 * с точным ближайшим t: на синтетике (2.2–17.7 тыс. кусков) это честно;
 * марш луча по пирамиде — следующий шаг (А1537-класс, §847).
 *
 * Приёмки §847: К27 на двух касающихся шарах (лучи в диск касания — в
 * БЛИЖНИЙ); сличение сборщика со свипом (средняя яркость cavity =
 * Le + rho·E_avg/2π, ±5 %); НК le=0 / noprop → кадр тождественно 0.
 *
 * Ключи:
 *   lev=N dirs=6|26|NxM it=N rho=F le=F tau0 noprop — как в swee3;
 *   mode=2 build=1 — только объёмный фронт (умолчания);
 *   eye=X,Y,Z look=X,Y,Z fov=F W=N H=N — камера (fov в градусах, по
 *   горизонтали); out=ФАЙЛ — вывод (умолчание img/pgather.ppm, PPM P6);
 *   k27=N — К27-приёмка: N×N контрольных лучей через диск касания шаров.
 *
 * КАРТИНКИ ПИШУТСЯ В img/, А НЕ В build/ (правило проекта).
 */
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "nstruct/pyr.h"
#include "nstruct/sweep.h"
#include "scene_obj.h"

#define PG_FOV_MIN 5.0   /* ниже — телеобъектив вне смысла синтетических тестов */
#define PG_FOV_MAX 170.0 /* выше — кадр шире полусферы, проекция вырождается */
#define PG_WH_MAX 4096   /* потолок кадра: перебор nt на луч, больше не нужно */
#define PG_K27_RAYS 64   /* контрольных лучей на сторону диска касания */

static double now_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static int parse3(const char *s, double out[3]) {
  return sscanf(s, "%lf,%lf,%lf", &out[0], &out[1], &out[2]) == 3 ? 0 : 1;
}

/* пересечение луча (o, d — единичный) с треугольником (Мёллер—Трумбор);
 * возврат t ≥ tmin или -1; p — вершины [3][3] */
static double pg_ray_tri(const double o[3], const double d[3], const double p[3][3]) {
  double e1[3], e2[3], pv[3], tv[3], qv[3];
  double det, u, v, t;
  int ax;
  for (ax = 0; ax < 3; ax++) {
    e1[ax] = p[1][ax] - p[0][ax];
    e2[ax] = p[2][ax] - p[0][ax];
  }
  pv[0] = d[1] * e2[2] - d[2] * e2[1];
  pv[1] = d[2] * e2[0] - d[0] * e2[2];
  pv[2] = d[0] * e2[1] - d[1] * e2[0];
  det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
  if (fabs(det) < 1e-30) return -1.0; /* луч в плоскости треугольника */
  u = 1.0 / det;
  tv[0] = o[0] - p[0][0];
  tv[1] = o[1] - p[0][1];
  tv[2] = o[2] - p[0][2];
  u = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * u;
  if (u < 0.0 || u > 1.0) return -1.0;
  qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
  qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
  qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
  v = (d[0] * qv[0] + d[1] * qv[1] + d[2] * qv[2]) / det;
  if (v < 0.0 || u + v > 1.0) return -1.0;
  t = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) / det;
  if (t < 1e-9) return -1.0;
  return t;
}

int main(int argc, char **argv) {
  const char *path = NULL, *outfile = "img/pgather.ppm";
  double scale = 1.0, le = 1.0, rho = -1.0, fov = 60.0;
  double eye[3] = {0, 0, 0}, look[3] = {0, 0, 0};
  int iters = 30, lev = 6, tau0 = 0, noprop = 0, ndirs = 26, mort = 1, i, ax;
  int W = 320, H = 240, k27 = 0;
  hz_objmesh m;
  hz_pyr py;
  hz_sw_opts so;
  hz_sw_stat st;
  double *area = NULL, *nrm = NULL, *kd = NULL, *cent = NULL, *cmin = NULL, *cmax = NULL;
  int32_t *mtl = NULL;
  double cell, t0, t1, sw_time;
  double *lum = NULL;

  for (i = 1; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0)
      lev = atoi(argv[i] + 4);
    else if (strncmp(argv[i], "rho=", 4) == 0)
      rho = atof(argv[i] + 4);
    else if (strncmp(argv[i], "le=", 3) == 0)
      le = atof(argv[i] + 3);
    else if (strncmp(argv[i], "it=", 3) == 0)
      iters = atoi(argv[i] + 3);
    else if (strncmp(argv[i], "dirs=", 5) == 0) {
      int a, b;
      if (sscanf(argv[i] + 5, "%dx%d", &a, &b) == 2 && a > 0 && b > 0)
        ndirs = a * 100 + b;
      else
        ndirs = atoi(argv[i] + 5);
    } else if (strncmp(argv[i], "eye=", 4) == 0)
      parse3(argv[i] + 4, eye);
    else if (strncmp(argv[i], "look=", 5) == 0)
      parse3(argv[i] + 5, look);
    else if (strncmp(argv[i], "fov=", 4) == 0)
      fov = atof(argv[i] + 4);
    else if (strncmp(argv[i], "W=", 2) == 0)
      W = atoi(argv[i] + 2);
    else if (strncmp(argv[i], "H=", 2) == 0)
      H = atoi(argv[i] + 2);
    else if (strncmp(argv[i], "k27=", 4) == 0)
      k27 = atoi(argv[i] + 4);
    else if (strncmp(argv[i], "out=", 4) == 0)
      outfile = argv[i] + 4;
    else if (strcmp(argv[i], "tau0") == 0)
      tau0 = 1;
    else if (strcmp(argv[i], "noprop") == 0)
      noprop = 1;
    else if (strncmp(argv[i], "mort=", 5) == 0)
      mort = atoi(argv[i] + 5);
    else if (strncmp(argv[i], "scale=", 6) == 0)
      scale = atof(argv[i] + 6);
    else
      path = argv[i];
  }
  if (!path || W < 1 || H < 1 || W > PG_WH_MAX || H > PG_WH_MAX || fov < PG_FOV_MIN ||
      fov > PG_FOV_MAX) {
    fprintf(stderr, "use: pgather <scene.obj> [lev=N dirs=.. it=N rho=F le=F] "
                    "[eye=X,Y,Z look=X,Y,Z fov=F W=N H=N] [k27=N] [out=ФАЙЛ]\n");
    return 2;
  }
  if (hz_obj_load(&m, path, scale) != 0) {
    fprintf(stderr, "pgather: не читается %s\n", path);
    return 2;
  }

  area = (double *)malloc((size_t)m.nt * sizeof *area);
  nrm = (double *)malloc((size_t)m.nt * 3 * sizeof *nrm);
  kd = (double *)malloc((size_t)m.nt * sizeof *kd);
  cent = (double *)malloc((size_t)m.nt * 3 * sizeof *cent);
  cmin = (double *)malloc((size_t)m.nt * 3 * sizeof *cmin);
  cmax = (double *)malloc((size_t)m.nt * 3 * sizeof *cmax);
  mtl = (int32_t *)malloc((size_t)m.nt * sizeof *mtl);
  if (!area || !nrm || !kd || !cent || !cmin || !cmax || !mtl) {
    fprintf(stderr, "pgather: нет памяти\n");
    return 2;
  }
  for (i = 0; i < m.nt; i++) {
    double p[3][3], e1[3], e2[3], nn;
    int v;
    hz_obj_tri(&m, (int32_t)i, p);
    for (ax = 0; ax < 3; ax++) {
      double lo = p[0][ax], hi = p[0][ax];
      for (v = 1; v < 3; v++) {
        if (p[v][ax] < lo) lo = p[v][ax];
        if (p[v][ax] > hi) hi = p[v][ax];
      }
      cmin[3 * (int64_t)i + ax] = lo;
      cmax[3 * (int64_t)i + ax] = hi;
      cent[3 * (int64_t)i + ax] = (p[0][ax] + p[1][ax] + p[2][ax]) / 3.0;
    }
    for (ax = 0; ax < 3; ax++) {
      e1[ax] = p[1][ax] - p[0][ax];
      e2[ax] = p[2][ax] - p[0][ax];
    }
    nrm[3 * (int64_t)i] = e1[1] * e2[2] - e1[2] * e2[1];
    nrm[3 * (int64_t)i + 1] = e1[2] * e2[0] - e1[0] * e2[2];
    nrm[3 * (int64_t)i + 2] = e1[0] * e2[1] - e1[1] * e2[0];
    nn = sqrt(nrm[3 * (int64_t)i] * nrm[3 * (int64_t)i] +
              nrm[3 * (int64_t)i + 1] * nrm[3 * (int64_t)i + 1] +
              nrm[3 * (int64_t)i + 2] * nrm[3 * (int64_t)i + 2]);
    area[i] = 0.5 * nn;
    if (nn > 0)
      for (ax = 0; ax < 3; ax++)
        nrm[3 * (int64_t)i + ax] /= nn;
    kd[i] = m.mtl[m.fm[i]].kd;
    mtl[i] = m.fm[i];
  }

  {
    double maxdim = 0.0;
    for (ax = 0; ax < 3; ax++) {
      double s = m.hi[ax] - m.lo[ax];
      if (s > maxdim) maxdim = s;
    }
    cell = maxdim / (double)(1 << lev);
  }
  if (hz_pyr_build(&py, m.nt, cmin, cmax, cent, mtl, m.lo, m.hi, cell) != 0) {
    fprintf(stderr, "pgather: пирамида не построилась\n");
    return 2;
  }
  if (mort && (hz_pyr_morton(&py) != 0 || hz_pyr_permute(&py, area, sizeof *area) != 0 ||
               hz_pyr_permute(&py, nrm, 3 * sizeof *nrm) != 0 ||
               hz_pyr_permute(&py, kd, sizeof *kd) != 0)) {
    fprintf(stderr, "pgather: Morton не прошёл\n");
    return 2;
  }

  /* --- СВИП: поле E на кусках (только объёмный фронт §845/§846) --- */
  memset(&so, 0, sizeof so);
  so.le = le;
  so.rho = rho;
  so.iters = iters;
  so.tau0 = tau0;
  so.noprop = noprop;
  so.ndirs = ndirs;
  so.mode = 2;
  so.build = 1;
  for (i = 0; i < m.nt; i++)
    py.pcs[i].e = 0.0f;
  t0 = now_sec();
  if (hz_sw_run(&py, m.nt, area, nrm, kd, &so, &st, NULL) != 0) {
    fprintf(stderr, "pgather: свип не прошёл\n");
    return 2;
  }
  t1 = now_sec();
  sw_time = t1 - t0;
  printf("СВИП: E_avg=%.4f (%.3f с)\n", st.e_avg, sw_time);

  /* --- К27: два касающихся шара — лучи через диск касания идут в БЛИЖНИЙ.
   * Шары сцены spheres.obj: центры и радиус восстанавливаются из bbox
   * кусков по знаку координаты x (левый/правый); касание — точка середины
   * между центрами. Контрольные лучи идут ИЗ глаза, стоящего на оси
   * ЛЕВОГО шара дальше от точки касания, в точки диска радиуса
   * 0.9·r вокруг направления на точку касания: обязательное попадание в
   * ближайший шар (t ближнего < t дальнего и сам луч вообще не доходит до
   * дальнего: t_ближ < t_даль − эпсилон геометрии). */
  if (k27 > 0) {
    double lo[3], hi[3];
    int v;
    for (ax = 0; ax < 3; ax++) {
      lo[ax] = hi[ax] = cent[3 * (int64_t)0 + ax];
      /* bbox по ЦЕНТРАМ кусков — куски сгруппированы вокруг двух шаров */
    }
    {
      /* кластеризация по x: порог — середина bbox по x */
      double xmid, csum[2][3] = {{0}};
      int64_t cnt[2] = {0, 0};
      double rad[2] = {0, 0};
      double c1[3], dir[3], tN, tF;
      int nk = 0, through = 0;
      xmid = 0.5 * (cmin[0] + cmax[0]);
      for (i = 0; i < m.nt; i++) {
        int s = cent[3 * (int64_t)i] < xmid ? 0 : 1;
        for (ax = 0; ax < 3; ax++)
          csum[s][ax] += cent[3 * (int64_t)i + ax];
        cnt[s]++;
      }
      for (ax = 0; ax < 3; ax++) {
        c1[ax] = csum[0][ax] / (double)cnt[0];
        lo[ax] = csum[1][ax] / (double)cnt[1];
      }
      for (i = 0; i < m.nt; i++) {
        int s = cent[3 * (int64_t)i] < xmid ? 0 : 1;
        double dd = 0;
        for (ax = 0; ax < 3; ax++) {
          double dxx = cent[3 * (int64_t)i + ax] - (s ? lo[ax] : c1[ax]);
          dd += dxx * dxx;
        }
        if (dd > rad[s]) rad[s] = dd;
      }
      rad[0] = sqrt(rad[0]);
      rad[1] = sqrt(rad[1]);
      printf("К27: центры (%.3f %.3f %.3f) r=%.3f и (%.3f %.3f %.3f) r=%.3f\n", c1[0], c1[1], c1[2],
             rad[0], lo[0], lo[1], lo[2], rad[1]);
      /* глаз — на продолжении линии центров за ЛЕВЫМ шаром */
      for (ax = 0; ax < 3; ax++)
        dir[ax] = lo[ax] - c1[ax];
      {
        double nn = sqrt(dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]);
        for (ax = 0; ax < 3; ax++)
          dir[ax] /= nn;
      }
      for (ax = 0; ax < 3; ax++)
        eye[ax] = c1[ax] - dir[ax] * (rad[0] * 4.0);
      /* точка касания — середина между центрами; лучи в диск 0.9r вокруг неё */
      for (ax = 0; ax < 3; ax++)
        look[ax] = 0.5 * (c1[ax] + lo[ax]);
      {
        /* базис диска: dir + два перпендикуляра */
        double u[3], w[3], tmp[3] = {0, 0, 1};
        double a0, ud;
        u[0] = dir[1] * tmp[2] - dir[2] * tmp[1];
        u[1] = dir[2] * tmp[0] - dir[0] * tmp[2];
        u[2] = dir[0] * tmp[1] - dir[1] * tmp[0];
        ud = sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        for (ax = 0; ax < 3; ax++)
          u[ax] /= ud;
        w[0] = dir[1] * u[2] - dir[2] * u[1];
        w[1] = dir[2] * u[0] - dir[0] * u[2];
        w[2] = dir[0] * u[1] - dir[1] * u[0];
        for (int iy = 0; iy < PG_K27_RAYS; iy++)
          for (int ix = 0; ix < PG_K27_RAYS; ix++) {
            double px = (2.0 * (ix + 0.5) / PG_K27_RAYS - 1.0) * 0.9 * fmin(rad[0], rad[1]);
            double py2 = (2.0 * (iy + 0.5) / PG_K27_RAYS - 1.0) * 0.9 * fmin(rad[0], rad[1]);
            if (px * px + py2 * py2 > 0.81 * fmin(rad[0], rad[1]) * fmin(rad[0], rad[1]))
              continue; /* круг, не квадрат */
            double rd[3];
            for (ax = 0; ax < 3; ax++)
              rd[ax] = look[ax] + px * u[ax] + py2 * w[ax] - eye[ax];
            {
              double nn2 = sqrt(rd[0] * rd[0] + rd[1] * rd[1] + rd[2] * rd[2]);
              for (ax = 0; ax < 3; ax++)
                rd[ax] /= nn2;
            }
            /* ближайшие два пересечения */
            tN = -1.0;
            tF = -1.0;
            for (i = 0; i < m.nt; i++) {
              double p3[3][3], tt;
              hz_obj_tri(&m, py.pcs[i].tri, p3);
              tt = pg_ray_tri(eye, rd, p3);
              if (tt < 0.0) continue;
              if (tN < 0.0 || tt < tN) {
                tF = tN;
                tN = tt;
              } else if (tF < 0.0 || tt < tF) {
                tF = tt;
              }
            }
            nk++;
            if (tN > 0.0 && tF > 0.0 && tN < tF * (1.0 - 1e-12)) continue;
            through++;
          }
        a0 = 100.0 * (double)through / (double)(nk ? nk : 1);
        printf("К27: лучей %d, мимо ближнего %d (%.2f %%) — ОБЯЗАНА быть 0 %%\n", nk, through, a0);
      }
      (void)v;
    }
  }

  /* --- СБОР: перебор кусков, ближайшее t, L_out = le + rho·E/(2π) --- */
  lum = (double *)malloc((size_t)W * (size_t)H * sizeof *lum);
  if (!lum) {
    fprintf(stderr, "pgather: нет памяти на кадр\n");
    return 2;
  }
  {
    /* базис камеры */
    double fwd[3], right[3], up[3], tmp[3] = {0, 0, 1};
    double tanf = tan(fov * M_PI / 360.0);
    int64_t nhit = 0;
    double lsum = 0, lmax = 0;
    for (ax = 0; ax < 3; ax++)
      fwd[ax] = look[ax] - eye[ax];
    {
      double nn = sqrt(fwd[0] * fwd[0] + fwd[1] * fwd[1] + fwd[2] * fwd[2]);
      if (nn < 1e-12) {
        fprintf(stderr, "pgather: глаз совпадает с точкой взгляда\n");
        return 2;
      }
      for (ax = 0; ax < 3; ax++)
        fwd[ax] /= nn;
    }
    right[0] = fwd[1] * tmp[2] - fwd[2] * tmp[1];
    right[1] = fwd[2] * tmp[0] - fwd[0] * tmp[2];
    right[2] = fwd[0] * tmp[1] - fwd[1] * tmp[0];
    {
      double nn = sqrt(right[0] * right[0] + right[1] * right[1] + right[2] * right[2]);
      if (nn < 1e-9) { /* взгляд вдоль z — базис от x */
        tmp[0] = 1;
        tmp[1] = tmp[2] = 0;
        right[0] = fwd[1] * tmp[2] - fwd[2] * tmp[1];
        right[1] = fwd[2] * tmp[0] - fwd[0] * tmp[2];
        right[2] = fwd[0] * tmp[1] - fwd[1] * tmp[0];
        nn = sqrt(right[0] * right[0] + right[1] * right[1] + right[2] * right[2]);
      }
      for (ax = 0; ax < 3; ax++)
        right[ax] /= nn;
    }
    up[0] = fwd[1] * right[2] - fwd[2] * right[1];
    up[1] = fwd[2] * right[0] - fwd[0] * right[2];
    up[2] = fwd[0] * right[1] - fwd[1] * right[0];

    t0 = now_sec();
    for (int iy = 0; iy < H; iy++)
      for (int ix = 0; ix < W; ix++) {
        double sx = (2.0 * (ix + 0.5) / W - 1.0) * tanf;
        double sy = (1.0 - 2.0 * (iy + 0.5) / H) * tanf * (double)H / (double)W;
        double rd[3], tbest = -1.0;
        int pbest = -1;
        for (ax = 0; ax < 3; ax++)
          rd[ax] = fwd[ax] + sx * right[ax] + sy * up[ax];
        {
          double nn = sqrt(rd[0] * rd[0] + rd[1] * rd[1] + rd[2] * rd[2]);
          for (ax = 0; ax < 3; ax++)
            rd[ax] /= nn;
        }
        for (i = 0; i < m.nt; i++) {
          double p3[3][3], tt;
          hz_obj_tri(&m, py.pcs[i].tri, p3);
          tt = pg_ray_tri(eye, rd, p3);
          if (tt >= 0.0 && (tbest < 0.0 || tt < tbest)) {
            tbest = tt;
            pbest = i;
          }
        }
        {
          double L = 0.0;
          if (pbest >= 0) {
            int q = py.pcs[pbest].tri; /* яркость — по ИСХОДНОМУ индексу для отчёта,
                                          поле — у куска pbest */
            (void)q;
            L = le + (rho < 0 ? kd[pbest] : rho) * (double)py.pcs[pbest].e / (2.0 * M_PI);
            nhit++;
            lsum += L;
            if (L > lmax) lmax = L;
          }
          lum[(size_t)iy * (size_t)W + (size_t)ix] = L;
        }
      }
    t1 = now_sec();
    printf("СБОР: лучей %d, попало %" PRId64 " (%.2f %%), средняя яркость %.4f, max %.4f, %.2f с\n",
           W * H, nhit, 100.0 * (double)nhit / ((double)W * (double)H),
           lsum / (nhit ? (double)nhit : 1.0), lmax, t1 - t0);
    {
      double rr = rho < 0 ? 0.5 : rho;
      double Lpred = le + rr * st.e_avg / (2.0 * M_PI);
      double Lmeas = lsum / (nhit ? (double)nhit : 1.0);
      printf("СЛИЧЕНИЕ: L_свипа = le + rho·E_avg/2π = %.4f; L_сбора = %.4f; отношение %.4f\n",
             Lpred, Lmeas, Lmeas / Lpred);
    }
    {
      FILE *f = fopen(outfile, "wb");
      if (!f) {
        fprintf(stderr, "pgather: не открыть %s\n", outfile);
        return 2;
      }
      fprintf(f, "P6\n%d %d\n255\n", W, H);
      for (i = 0; i < W * H; i++) {
        double v = lum[i] / (lmax > 0 ? lmax : 1.0) * 255.0;
        unsigned char b[3];
        if (v > 255.0) v = 255.0;
        if (v < 0.0) v = 0.0;
        b[0] = b[1] = b[2] = (unsigned char)v;
        fwrite(b, 1, 3, f);
      }
      fclose(f);
      printf("КАДР: %s записан (P6, нормировка на max кадра)\n", outfile);
    }
  }

  free(lum);
  free(area);
  free(nrm);
  free(kd);
  free(cent);
  free(cmin);
  free(cmax);
  free(mtl);
  hz_pyr_free(&py);
  hz_obj_free(&m);
  return 0;
}
