/* swee3 — СВИП на новых массивах (§835-§837): итерации + аналитика полости
 * + трафик-прибор §5 STRUCTURE.md + режим сличения со старым путём (cmp=).
 *
 * Ключи: it=N (умолчание 20), rho=F (<0 — kd материалов; 0 — НК),
 * le=F (умолчание 1), tau0 (НК: слой не взаимодействует),
 * noprop (НК: луч не переносится), dirs=6|26 (§837) или dirs=NxM —
 * продуктовая квадратура §843 (Чебышёв×Гаусс, nd = N·M), lev=N,
 * cmp=ФАЙЛ (сличение E с дампом старого пути по индексу куска).
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

static double now_sec(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

int main(int argc, char **argv) {
  const char *path = NULL, *cmpfile = NULL;
  double scale = 1.0, le = 1.0, rho = -1.0;
  int iters = 20, lev = 6, tau0 = 0, noprop = 0, ndirs = 6, mort = 1, i, ax;
  int mode = 1;          /* §845: 2 — объёмный фронт (L на клетку, марш-порядок) */
  int build = 0;         /* §845-в: строитель похода (1 — прямой сортировочный) */
  int vc = 1;            /* §851: компоненты пустоты по умолчанию включены */
  int lp = 0;            /* §852: LOD-уровень кусков ℓ_p (единый; лестница G1) */
  int adapt = 0;         /* §862: 1 — дробный аккумулятор детальности */
  double *lpacc = NULL;  /* §862: аккумулятор [nt], живёт до конца main */
  double *lphits = NULL; /* §862-диаг: счётчик событий [nt] */
  int dumpE = 0;         /* §862-диаг: дамп per-piece E последней итерации */
  double cdelta = 1.0;   /* §862: вес приращения Δ(материал,угол) */
  int amode = 0;         /* §862-диаг: форма Δ (см. sweep.h accum_mode) */
  int agg = 0;           /* §863/шаг 2: 1 — группы по материалу; 2 — одна
                          * группа на список (диффузное приближение) */
  int xint = 0;          /* §852/G1(a): 1 — всегда точное пересечение */
  int screw = 0;         /* НК А1572: скрестить куски с шагом screw (детекторы>0) */
  int walk = 0;          /* А1576: 1 — точный многопопадный проход (модель «реальные
                          * пересечения»); 0 — схема §852 (перехват под предикатом) */
  int useke = 0;         /* А1576: 1 — per-piece le из Ke материалов (плюс глобальный
                          * le); 0 — прежний мир (глобальный le) */
  int nphi = 0, nmu = 0;
  int has_recv = 0; /* blocked-beam (А1566): recv=X:<x0>:<x1> — приёмник-слэб */
  double recv0 = 0.0, recv1 = 0.0;
  double t0, t1;
  hz_objmesh m;
  hz_pyr py;
  hz_sw_opts so;
  hz_sw_stat st;
  double *area = NULL, *nrm = NULL, *kd = NULL, *cent = NULL, *cmin = NULL, *cmax = NULL;
  int32_t *mtl = NULL;
  double *hist = NULL;
  double *lep9 = NULL; /* А1576: per-piece эмиссия (Ke) в порядке ИСХОДНЫХ */
  double cell;

  for (i = 1; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0)
      lev = atoi(argv[i] + 4);
    else if (strncmp(argv[i], "rho=", 4) == 0)
      rho = atof(argv[i] + 4);
    else if (strncmp(argv[i], "le=", 3) == 0)
      le = atof(argv[i] + 3);
    else if (strncmp(argv[i], "it=", 3) == 0)
      iters = atoi(argv[i] + 3);
    else if (strcmp(argv[i], "tau0") == 0)
      tau0 = 1;
    else if (strcmp(argv[i], "noprop") == 0)
      noprop = 1;
    else if (strncmp(argv[i], "dirs=", 5) == 0) {
      /* §843: dirs=6|26 — легаси; dirs=NxM — продуктовая квадратура
       * (Чебышёв по φ × Гаусс по μ, nd = N·M), кодируется N*100+M */
      if (sscanf(argv[i] + 5, "%dx%d", &nphi, &nmu) == 2 && nphi > 0 && nmu > 0)
        ndirs = nphi * 100 + nmu;
      else
        ndirs = atoi(argv[i] + 5);
    } else if (strncmp(argv[i], "build=", 6) == 0) {
      build = atoi(argv[i] + 6); /* §845-в: 1 — прямой сортировочный строитель */
    } else if (strncmp(argv[i], "vc=", 3) == 0) {
      vc = atoi(argv[i] + 3); /* §851: 0 — старая однокомпонентная логика */
    } else if (strncmp(argv[i], "lp=", 3) == 0) {
      lp = atoi(argv[i] + 3); /* §852: ℓ_p кусков */
    } else if (strncmp(argv[i], "adapt=", 6) == 0) {
      adapt = atoi(argv[i] + 6); /* §862 */
    } else if (strncmp(argv[i], "cdelta=", 7) == 0) {
      cdelta = atof(argv[i] + 7); /* §862 */
    } else if (strncmp(argv[i], "agg=", 4) == 0) {
      agg = atoi(argv[i] + 4); /* §863/шаг 2 */
    } else if (strncmp(argv[i], "amode=", 6) == 0) {
      amode = atoi(argv[i] + 6); /* §862-диаг */
    } else if (strncmp(argv[i], "dumpE=", 6) == 0) {
      dumpE = atoi(argv[i] + 6); /* §862-диаг */
    } else if (strncmp(argv[i], "xint=", 5) == 0) {
      xint = atoi(argv[i] + 5); /* §852/G1(a): точное пересечение везде */
    } else if (strncmp(argv[i], "screw=", 6) == 0) {
      screw = atoi(argv[i] + 6); /* НК А1572 */
    } else if (strncmp(argv[i], "walk=", 5) == 0) {
      walk = atoi(argv[i] + 5); /* А1576: модель «реальные пересечения» */
    } else if (strncmp(argv[i], "ke=", 3) == 0) {
      useke = atoi(argv[i] + 3); /* А1576: per-piece le из Ke */
    } else if (strncmp(argv[i], "recv=", 5) == 0) {
      /* А1566-фальсификатор: средняя E по кускам, чей ИСХОДНЫЙ центроид в
       * слэбе по оси X (приёмник за перегородкой) */
      if (sscanf(argv[i] + 5, "X:%lf:%lf", &recv0, &recv1) == 2) has_recv = 1;
    } else if (strncmp(argv[i], "mode=", 5) == 0) {
      mode = atoi(argv[i] + 5); /* §845: 2 — объёмный фронт */
    } else if (strncmp(argv[i], "cmp=", 4) == 0)
      cmpfile = argv[i] + 4;
    else if (strncmp(argv[i], "mort=", 5) == 0)
      mort = atoi(argv[i] + 5);
    else if (strncmp(argv[i], "scale=", 6) == 0)
      scale = atof(argv[i] + 6);
    else
      path = argv[i];
  }
  if (!path) {
    fprintf(stderr, "use: swee3 <scene.obj> [it=N] [rho=F] [le=F] [tau0] [noprop] "
                    "[dirs=6|26|NxM] [cmp=ФАЙЛ] [lev=N]\n");
    return 2;
  }
  if (hz_obj_load(&m, path, scale) != 0) {
    fprintf(stderr, "swee3: не читается %s\n", path);
    return 2;
  }

  { /* §843-G4: моменты квадратуры — Σw = 4π точно; Σw|n·ω̂| ≈ 2π (излом |μ|
     * в нуле даёт ошибку порядка квадрата шага по μ, см. А1510) */
    hz_sw_dir *tab = NULL;
    int nd2 = 0;
    static const double raw[3][3] = {{1, 2, 3}, {-5, 1, 2}, {3, 7, 11}};
    if (hz_sw_dir_table(ndirs, &tab, &nd2) != 0) {
      fprintf(stderr, "swee3: таблица направлений не строится (dirs=%d)\n", ndirs);
      return 2;
    }
    double sw = 0.0;
    for (i = 0; i < nd2; i++)
      sw += tab[i].w;
    printf("КВАДРАТУРА: dirs=%d nd=%d Σw=%.15g (4π=%.15g)", ndirs, nd2, sw,
           4.0 * 3.14159265358979323846);
    for (ax = 0; ax < 3; ax++) {
      double nn = 0.0, sm = 0.0;
      int q;
      for (q = 0; q < 3; q++)
        nn += raw[ax][q] * raw[ax][q];
      nn = sqrt(nn);
      for (i = 0; i < nd2; i++)
        sm += tab[i].w * fabs(tab[i].om[0] * raw[ax][0] / nn + tab[i].om[1] * raw[ax][1] / nn +
                              tab[i].om[2] * raw[ax][2] / nn);
      printf("; Σw|cos|=%.12g", sm);
    }
    printf(" (2π=%.12g)\n", 2.0 * 3.14159265358979323846);
    free(tab);
  }

  area = (double *)malloc((size_t)m.nt * sizeof *area);
  nrm = (double *)malloc((size_t)m.nt * 3 * sizeof *nrm);
  kd = (double *)malloc((size_t)m.nt * sizeof *kd);
  cent = (double *)malloc((size_t)m.nt * 3 * sizeof *cent);
  cmin = (double *)malloc((size_t)m.nt * 3 * sizeof *cmin);
  cmax = (double *)malloc((size_t)m.nt * 3 * sizeof *cmax);
  mtl = (int32_t *)malloc((size_t)m.nt * sizeof *mtl);
  hist = (double *)calloc((size_t)iters, sizeof *hist);
  if (!area || !nrm || !kd || !cent || !cmin || !cmax || !mtl || !hist) {
    fprintf(stderr, "swee3: нет памяти\n");
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
  memset(&so, 0, sizeof so); /* рано: дальше блоки §852 заполняют so сами */
  {                          /* §852: вершины и bbox-ы треугольников (в порядке ИСХОДНЫХ) — точное
                              * пересечение лучевого свипа и предикат А1566 */
    double *tv9 = (double *)malloc((size_t)m.nt * 9 * sizeof *tv9);
    double *tb6 = (double *)malloc((size_t)m.nt * 6 * sizeof *tb6);
    uint8_t *lparr = (uint8_t *)malloc((size_t)m.nt);
    if (adapt) lpacc = (double *)calloc((size_t)m.nt, sizeof *lpacc);   /* §862 */
    if (adapt) lphits = (double *)calloc((size_t)m.nt, sizeof *lphits); /* §862-диаг */
    if (adapt && !lphits) {
      fprintf(stderr, "нет памяти на lphits\n");
      return 2;
    }
    int32_t ti;
    if (!tv9 || !tb6 || !lparr || (adapt && !lpacc)) {
      fprintf(stderr, "нет памяти на trivert/tribox/lp/lpacc\n");
      return 2;
    }
    for (ti = 0; ti < m.nt; ti++) {
      double pp[3][3];
      int aa;
      hz_obj_tri(&m, ti, pp);
      for (int v = 0; v < 3; v++)
        for (aa = 0; aa < 3; aa++)
          tv9[9 * (int64_t)ti + 3 * v + aa] = pp[v][aa];
      for (ax = 0; ax < 3; ax++) {
        double lo = pp[0][ax], hi = pp[0][ax];
        for (int v = 1; v < 3; v++) {
          if (pp[v][ax] < lo) lo = pp[v][ax];
          if (pp[v][ax] > hi) hi = pp[v][ax];
        }
        tb6[6 * (int64_t)ti + ax] = lo;
        tb6[6 * (int64_t)ti + 3 + ax] = hi;
      }
      lparr[ti] = (uint8_t)lp; /* §852: единый ℓ_p прогоном (лестница G1);
                                * поэлементная политика §844 — отдельно */
    }
    so.trivert = tv9;
    so.tribox = tb6;
    so.domhi = m.hi;
    so.lp = lparr;
    if (adapt) { /* §862: дробный аккумулятор детальности */
      for (ti = 0; ti < (int)m.nt; ti++)
        lpacc[ti] = (double)lp;
      so.lpacc = lpacc;
      so.cdelta = cdelta;
      so.accum_mode = amode;
      so.agg = agg;
      so.lphits = lphits;
      so.lpapply = (adapt == 1); /* adapt=2 — пассивная диагностика */
    }
    so.xint = xint;
    so.walk = walk;
    if (useke) {
      /* А1576: per-piece эмиссия = СРЕДНЕЕ Ke материала (прецедент
       * скаляризации pfield.c) ПЛЮС глобальный le — на сценах с Ke=0
       * побитово прежний мир. Порядок — ИСХОДНЫЕ треугольники, переставится
       * вместе с area/nrm/kd под Morton ниже. */
      lep9 = (double *)malloc((size_t)m.nt * sizeof *lep9);
      if (!lep9) {
        fprintf(stderr, "нет памяти на lep\n");
        return 2;
      }
      for (int32_t tq = 0; tq < m.nt; tq++) {
        const double *kq = m.mtl[m.fm[tq]].ke3;
        lep9[tq] = (kq[0] + kq[1] + kq[2]) / 3.0 + le;
      }
      so.lep = lep9;
    }
  }
  if (hz_pyr_build(&py, m.nt, cmin, cmax, cent, mtl, m.lo, m.hi, cell) != 0) {
    fprintf(stderr, "swee3: пирамида не построилась\n");
    return 2;
  }
  if (mort) {
    if (hz_pyr_morton(&py) != 0) {
      fprintf(stderr, "swee3: Morton не прошёл\n");
      return 2;
    }
    /* параллельные массивы инструмента — той же перестановкой, циклами
     * на месте (А1491: соответствие кусок↔треугольник — pcs[].tri) */
    if (hz_pyr_permute(&py, area, sizeof *area) != 0 ||
        hz_pyr_permute(&py, nrm, 3 * sizeof *nrm) != 0 ||
        hz_pyr_permute(&py, kd, sizeof *kd) != 0 ||
        (lep9 && hz_pyr_permute(&py, lep9, sizeof(double)) != 0)) {
      fprintf(stderr, "swee3: перестановка не прошла\n");
      return 2;
    }
    {
      double dmax = 0;
      int imax = -1;
      for (i = 0; i < m.nt; i++) {
        double dd = fabs(area[i] - hz_obj_tri_area(&m, py.pcs[i].tri));
        if (dd > dmax) {
          dmax = dd;
          imax = i;
        }
      }
      printf("ИНВАРИАНТ: max|area[i]-tri_area(pcs[i].tri)| = %.3g на куске %d\n", dmax, imax);
    }
  }
  { /* §852/А1567: per-node max ℓ_p одним подъёмом (после Morton — CSR свежий) */
    int32_t nup = 0;
    if (hz_pyr_set_lp(&py, (const uint8_t *)so.lp, &nup) != 0) {
      fprintf(stderr, "swee3: per-node max lp не построился\n");
      return 2;
    }
    /* А1568: клэмп вверх — отказ с печатаемым счётчиком (обработка на листе) */
    printf("LOD: lp=%d, кусков с ℓ_p выше доступного уровня (обработаны на листе): %d\n", lp,
           (int)nup);
  }
  if (screw) { /* НК А1572: калибровка детекторов — при screw они ОБЯЗАНЫ
                * сработать (d_cell/d_empty > 0) */
    int32_t ns = hz_pyr_screw(&py, (int32_t)screw);
    printf("НК: screw=%d — скрещено %d кусков\n", screw, (int)ns);
  }
  { /* НК А1572: детекторы пирамиды; на ЧИСТОМ прогоне — все нули */
    hz_pyr_verdict vd;
    hz_pyr_verify(&py, m.nt, cmin, cmax, cent, &vd);
    printf("НК hz_pyr_verify: d_cell=%" PRId64 " d_csr=%" PRId64 " d_bbox=%" PRId64
           " d_empty=%" PRId64 " (чистый прогон: все нули; при screw — ненули)\n",
           vd.d_cell, vd.d_csr, vd.d_bbox, vd.d_empty);
  }

  printf("== swee3 %s: nt=%d клетка %.4g м, it=%d le=%.3g rho=%s dirs=%d%s%s\n", path, m.nt, cell,
         iters, le, rho < 0 ? "kd" : "ovr", ndirs, tau0 ? " tau0" : "", noprop ? " noprop" : "");
  printf("   листья: занятых %d, с кусками %d — пустых в обходе %d (по 24 Б на визит)\n", py.nleaf,
         (int)py.ncentleaf, py.nleaf - (int)py.ncentleaf);
  { /* §843/А1513: кратность куска по листьям — свип посещает кусок в КАЖДОМ
     * листе его bbox; если кратность > 1, площадь куска участвует в переносе
     * многократно, и дискретизация может не сходиться с измельчением */
    int *pcnt = (int *)calloc((size_t)m.nt, sizeof *pcnt);
    int64_t slots = 0;
    int32_t li2, u2, pmax = 0;
    double pavg;
    if (!pcnt) {
      fprintf(stderr, "swee3: нет памяти\n");
      return 2;
    }
    for (li2 = 0; li2 < py.nleaf; li2++)
      for (u2 = 0; u2 < py.leaf[li2].npcs; u2++)
        pcnt[py.csr[py.leaf[li2].pcs_first + u2]]++;
    for (i = 0; i < m.nt; i++) {
      slots += pcnt[i];
      if (pcnt[i] > pmax) pmax = pcnt[i];
    }
    pavg = (double)slots / (double)m.nt;
    printf("   КУСКИ×ЛИСТЬЯ: слотов %" PRId64
           " на %d кусков — средняя кратность %.2f, максимум %d\n",
           slots, m.nt, pavg, pmax);
    free(pcnt);
  }

  so.build = build; /* so уже обнулён и заполнен §852-блоками выше */
  so.vc = vc;
  so.le = le;
  so.rho = rho;
  so.iters = iters;
  so.tau0 = tau0;
  so.noprop = noprop;
  so.ndirs = ndirs;
  for (int md = mode; md < mode + 1; md++) { /* один режим за прогон (§845: mode=2) */
    const char *mname = (md == 3) ? "pyrfront" : ((md == 2) ? "front" : (md ? "col" : "scalar"));
    so.mode = md;
    for (i = 0; i < m.nt; i++)
      py.pcs[i].e = 0.0f; /* режимы с чистого поля */
    t0 = now_sec();
    if (hz_sw_run(&py, m.nt, area, nrm, kd, &so, &st, hist) != 0) {
      fprintf(stderr, "swee3: свип не прошёл\n");
      return 2;
    }
    t1 = now_sec();
    printf("[%s dirs=%d] трафик: %.2f МБ/итерацию, %.3f с → %.2f ГБ/с эффективной\n", mname,
           md == 0 ? 6 : so.ndirs, (double)st.traffic / 1048576.0, t1 - t0,
           (t1 - t0) > 1e-9 ? (double)st.traffic * (double)iters / (t1 - t0) / 1e9 : 0.0);
    printf("[%s] итерации E_avg:", mname);
    for (i = 0; i < iters; i++)
      printf(" %.4f", hist[i]);
    printf("\n");
    if (iters >= 3) {
      double f = hist[iters - 1] - hist[iters - 2];
      double f0 = hist[iters - 2] - hist[iters - 3];
      printf("[%s] фактор ряда: %.4f (ожидалось ~ρ)\n", mname, fabs(f0) > 1e-300 ? f / f0 : 0.0);
    }
    printf("[%s] баланс: излучено %.4f, поглощено %.4f, рециркуляция %.4f, потеряно %.4f\n", mname,
           st.emitted, st.absorbed, st.recycled, st.lost);
    printf("[%s] доставки: визитов %" PRId64 ", с светом %" PRId64 " (%.1f %%)\n", mname, st.nvisit,
           st.ndep, st.nvisit ? 100.0 * (double)st.ndep / (double)st.nvisit : 0.0);
    printf("[%s] линии: %" PRId64 ", из 1 визита %.1f %%, из 2 — %.1f %%\n", mname, st.nline,
           st.nline ? 100.0 * (double)st.nline1 / (double)st.nline : 0.0,
           st.nline ? 100.0 * (double)st.nline2 / (double)st.nline : 0.0);
    printf("[%s] хеш порядка похода: %016" PRIx64 " (§845/А1525: режимы обязаны различаться)\n",
           mname, st.order_hash);
    printf("[%s] инверсии марша: %" PRId64 " (§845-бис: 0 ожидается на целочисленных)\n", mname,
           st.ninv);
    if (md == 3) { /* §852: приборы фронта */
      double vfrac = st.ncellbase > 0 ? 1.0 - (double)st.ncellfront / (double)st.ncellbase : 0.0;
      printf("[pyrfront] G6 нарушений: %" PRId64 " (чистый прогон: 0); спусков %" PRId64
             ", материальных событий %" PRId64 ", ПУСТ-прыжков %" PRId64 "\n",
             st.g6viol, st.ndesc, st.nmat, st.njump);
      printf("[pyrfront] прибор пустоты (А1569): клеток полным DDA %" PRId64
             ", посещено фронтом %" PRId64 " — пустых %.1f %% (предсказание > 90 %%)\n",
             st.ncellbase, st.ncellfront, 100.0 * vfrac);
      printf("[pyrfront] штамп А1564: пресечено повторных депозитов %" PRId64
             "; сегментов за границей домена %" PRId64 "\n",
             st.nstamp, st.nlostseg);
    }
    if (fabs(rho - 1.0) > 1e-12 && !tau0 && !noprop) {
      double rr = rho < 0 ? 0.5 : rho;
      /* §842: точное решение МОДЕЛИ слоя: E = W·Le/(A − ρ·W/(2π)),
       * W = Σ_p area_p·(выходная сумма весов) — квадратура учтена точно */
      double W = 0, A = 0;
      for (i = 0; i < m.nt; i++) {
        W += hz_sw_exitwsum(nrm + 3 * (int64_t)i, ndirs) * area[i];
        A += area[i];
      }
      double emodel = W * le / (A - rr * W / (2.0 * M_PI));
      printf("[%s] модель слоя: E = %.4f (W=%.3f A=%.3f); свип: %.4f → %.3f×модели\n", mname,
             emodel, W, A, st.e_avg, st.e_avg / emodel);
      printf("[%s] справочно, непрерывная физика: 2πLe/(1−ρ) = %.3f\n", mname,
             2.0 * M_PI * le / (1.0 - rr));
    }
    if (has_recv) {
      /* А1566-фальсификатор: приёмник — куски с ИСХОДНЫМ центроидом в слэбе
       * recv0<=x<=recv1 (за перегородкой); средняя E площадь-взвешенная */
      double esum = 0.0, asum = 0.0;
      int nrecv = 0;
      for (i = 0; i < m.nt; i++) {
        int32_t tri = py.pcs[i].tri;
        double cx = cent[3 * (int64_t)tri];
        if (cx >= recv0 && cx <= recv1) {
          esum += (double)py.pcs[i].e * area[i];
          asum += area[i];
          nrecv++;
        }
      }
      printf("[pyrfront] ПРИЁМНИК x∈[%.4g,%.4g]: кусков %d, E_avg = %.6g (Σ=%.6g)\n", recv0, recv1,
             nrecv, asum > 0 ? esum / asum : 0.0, esum);
    }
  }

  /* сличение со старым путём (§837-P): дамп E по индексу куска */
  if (cmpfile) {
    FILE *f = fopen(cmpfile, "rb");
    int32_t nold = 0;
    if (!f || fread(&nold, sizeof nold, 1, f) != 1 || nold != m.nt) {
      fprintf(stderr, "swee3: дамп %s не читается или nt не совпал (%d vs %d)\n", cmpfile, nold,
              m.nt);
      return 2;
    }
    {
      double *eold = (double *)malloc((size_t)nold * sizeof *eold);
      double rmed, rsum = 0, asum = 0;
      int32_t q;
      int in2 = 0;
      double *rats;
      if (!eold || fread(eold, sizeof *eold, (size_t)nold, f) != (size_t)nold) {
        fprintf(stderr, "swee3: дамп обрезан\n");
        return 2;
      }
      fclose(f);
      rats = (double *)malloc((size_t)nold * sizeof *rats);
      if (!rats) {
        fprintf(stderr, "swee3: нет памяти\n");
        free(eold);
        return 2;
      }
      for (q = 0; q < nold; q++) {
        int32_t tri = py.pcs[q].tri; /* сличение по ИСХОДНОМУ треугольнику (А1491) */
        double en = py.pcs[q].e;
        double r = (eold[tri] > 1e-9) ? en / eold[tri] : 1.0;
        rats[q] = r;
        rsum += r * area[q];
        asum += area[q];
        if (r > 0.5 && r < 2.0) in2++;
      }
      /* медиана отношений */
      {
        int32_t u, v;
        for (u = 1; u < nold; u++) {
          double ru = rats[u];
          for (v = u; v > 0 && rats[v - 1] > ru; v--)
            rats[v] = rats[v - 1];
          rats[v] = ru;
        }
        rmed = rats[nold / 2];
      }
      printf("ПАРИТЕТ: медиана отношения %.4f; площадь-взвешенное %.4f; доля в ×2: %.1f %%\n", rmed,
             rsum / asum, 100.0 * (double)in2 / (double)nold);
      free(rats);
      free(eold);
    }
  }

  free(area);
  free(nrm);
  free(kd);
  free(cent);
  free(cmin);
  free(cmax);
  free(mtl);
  free(hist);
  free((void *)(uintptr_t)so.trivert);
  free((void *)(uintptr_t)so.tribox);
  free((void *)(uintptr_t)so.lp);
  if (dumpE) /* §862-диаг: нагрузка куска, последняя итерация */
    for (i = 0; i < (int)m.nt; i++)
      printf("E1 p=%d tri=%d e=%.8g hits=%.0f sdelta=%.8g\n", i, py.pcs[i].tri, (double)py.pcs[i].e,
             lphits ? lphits[i] : 0.0, lpacc ? lpacc[i] : 0.0);
  if (lpacc) { /* §862: гистограмма этажей на последней итерации */
    int fl, cnt[16] = {0}, maxfl = 0;
    for (i = 0; i < (int)m.nt; i++) {
      fl = (int)lpacc[i];
      if (fl < 0)
        fl = 0; /* §862-диаг: вклад бывает < 0 (отрицательные
                 * Lin-депозиты) — вне гистограммы, не в стёк */
      if (fl > 15) fl = 15;
      cnt[fl]++;
      if (fl > maxfl) maxfl = fl;
    }
    printf("ЛОД-АДАПТ: cdelta=%.3g этажи", cdelta);
    for (fl = 0; fl <= maxfl; fl++)
      printf(" %d:%d", fl, cnt[fl]);
    printf(")\n");
    free(lpacc);
    free(lphits);
  }
  free(lep9); /* А1576: per-piece эмиссия */
  hz_pyr_free(&py);
  hz_obj_free(&m);
  return 0;
}
