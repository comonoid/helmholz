/* pfield — ИСТОЧНИК ЭРМИТОВА ПОЛЯ ИЗ МЕША (PLAN_ELEMENTS.md, §363, §365, §366).
 *
 * ЗАЧЕМ. Решение §363: представление сцены — эрмитово поле на восьмидереве, а не
 * индекс над треугольниками. Готовая машинерия (`src/cut/dc.c`) лежит без дела
 * ровно по одной причине: ей нужен ЗНАК «внутри/снаружи», а наши сцены —
 * незамкнутый суп треугольников. Этот стенд знак строит и проверяет, ЖИВА ЛИ
 * машинерия, — до всякой перестройки конвейера (указание пользователя 08-09).
 *
 * ЗНАК — ЗАЛИВКОЙ ПУСТОТЫ ОТ ГРАНИЦЫ, И ПЕРЕКРЫВАЕТСЯ ОН ЯЧЕЙКОЙ С ГЕОМЕТРИЕЙ,
 * А НЕ ПЕРЕСЕЧЁННЫМ РЕБРОМ (А596). Поверхность, целиком лежащая внутри ячейки и
 * не задевающая её рёбер (столешница, спинка стула — медианная грань зала
 * `0.131` м при ячейке `0.153` м), рёбер не пересекает и заливку не остановит:
 * она протечёт сквозь стену и объявит комнату наружей. Это отказ в пользу
 * пустоты, что запрещено (Г25). Занятость же считается честным отсечением.
 *
 * ЧЕГО ЭТОТ СТЕНД НЕ ПОКАЗЫВАЕТ (А597). Сетка `64³` даёт на зале ячейку
 * `0.153` м — ГРУБЕЕ медианной грани. Значит доля неманифолдных ячеек здесь
 * говорит о РАЗРЕШЕНИИ, а не о пригодности конструкции, и при срабатывании
 * условия «убивает» первым делом проверяется сетка, а не замысел.
 */
#include "cut/dc.h"
#include "pclip.h"
#include "scene_cfg.h"
#include "scene_obj.h"
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

static int cmp_d(const void *x, const void *y) {
  double a = *(const double *)x, b = *(const double *)y;
  return (a < b) ? -1 : ((a > b) ? 1 : 0);
}

/* Контекст источника: всё предвычислено, обратные вызовы только читают. */
typedef struct {
  int n;              /* сторона сетки ячеек */
  unsigned char *sgn; /* (n+1)^3: 1 — угол ВНУТРИ тела */
  /* эрмитовы данные минимальных рёбер: для оси a и нижнего конца c */
  double *et;   /* 3*(n+1)^3: доля вдоль ребра, <0 — пересечения нет */
  double *enrm; /* 9*(n+1)^3 */
} src_ctx;

static int src_sign(void *ctx, const int32_t c[3]) {
  const src_ctx *S = (const src_ctx *)ctx;
  int n1 = S->n + 1;
  if (c[0] < 0 || c[1] < 0 || c[2] < 0 || c[0] >= n1 || c[1] >= n1 || c[2] >= n1) return 0;
  return S->sgn[((size_t)c[2] * (size_t)n1 + (size_t)c[1]) * (size_t)n1 + (size_t)c[0]];
}

static int src_cross(void *ctx, const int32_t a[3], int axis, double *t, double nrm[3]) {
  const src_ctx *S = (const src_ctx *)ctx;
  int n1 = S->n + 1;
  if (a[0] < 0 || a[1] < 0 || a[2] < 0 || a[0] >= n1 || a[1] >= n1 || a[2] >= n1) return 1;
  size_t k = ((size_t)a[2] * (size_t)n1 + (size_t)a[1]) * (size_t)n1 + (size_t)a[0];
  double tv = S->et[3 * k + (size_t)axis];
  if (!(tv >= 0.0)) return 1;
  *t = tv;
  for (int c = 0; c < 3; c++)
    nrm[c] = S->enrm[9 * k + 3 * (size_t)axis + (size_t)c];
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "pfield ФАЙЛ.obj МАСШТАБ [lev=N] [noflood]\n");
    return 2;
  }
  int lev = 6, noflood = 0, invert = 0;
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0) lev = (int)strtol(argv[i] + 4, NULL, 10);
    /* НЕГАТИВНЫЙ КОНТРОЛЬ (§365): знак не строится вовсе — вершин обязано выйти 0. */
    if (strcmp(argv[i], "noflood") == 0) noflood = 1;
    /* СИЛЬНЫЙ КОНТРОЛЬ (А600): знак обращён — поверхность обязана выйти ТОЙ ЖЕ. */
    if (strcmp(argv[i], "invert") == 0) invert = 1;
  }
  if (lev < 1 || lev > HZ_DC_SAMPLE_MAX_LOG2SIZE) {
    fprintf(stderr, "lev вне предела построителя (1..%d)\n", HZ_DC_SAMPLE_MAX_LOG2SIZE);
    return 2;
  }
  hz_objmesh m;
  if (hz_obj_load(&m, argv[1], strtod(argv[2], NULL)) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  const int n = 1 << lev, n1 = n + 1;
  double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
  for (int32_t k = 0; k < m.nv; k++)
    for (int c = 0; c < 3; c++) {
      double x = m.v[3 * (size_t)k + (size_t)c];
      if (x < lo[c]) lo[c] = x;
      if (x > hi[c]) hi[c] = x;
    }
  /* Кубическая рама: сетка одна на все оси, иначе ключ ребра и координаты DC
   * перестают быть однородными. Габарит расширяется до куба с запасом в ячейку. */
  double side = 0.0;
  for (int c = 0; c < 3; c++)
    if (hi[c] - lo[c] > side) side = hi[c] - lo[c];
  double h = side / (double)(n - 2);
  for (int c = 0; c < 3; c++) {
    double mid = 0.5 * (lo[c] + hi[c]);
    lo[c] = mid - 0.5 * (double)n * h;
  }
  printf("== ПОЛЕ: %s, треугольников %d, сетка %d^3, ячейка %.4f м\n", argv[1], m.nt, n, h);

  /* --- 1. списки треугольников по ячейкам (А601: в дереве ссылок нет) --- */
  double t0 = now_s();
  size_t ncell = (size_t)n * (size_t)n * (size_t)n;
  int32_t *cnt = calloc(ncell + 1, sizeof *cnt);
  unsigned char *occ = calloc(ncell, 1);
  if (cnt == NULL || occ == NULL) return 1;
  for (int32_t t = 0; t < m.nt; t++) {
    int64_t i0[3] = {0, 0, 0}, i1[3] = {0, 0, 0};
    for (int c = 0; c < 3; c++) {
      double a = 1e300, b = -1e300;
      for (int v = 0; v < 3; v++) {
        double x = m.v[3 * (size_t)m.f[3 * (size_t)t + (size_t)v] + (size_t)c];
        if (x < a) a = x;
        if (x > b) b = x;
      }
      i0[c] = (int64_t)floor((a - lo[c]) / h);
      i1[c] = (int64_t)floor((b - lo[c]) / h);
      if (i0[c] < 0) i0[c] = 0;
      if (i1[c] >= n) i1[c] = n - 1;
    }
    const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
    const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
    const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
    for (int64_t z = i0[2]; z <= i1[2]; z++)
      for (int64_t y = i0[1]; y <= i1[1]; y++)
        for (int64_t x = i0[0]; x <= i1[0]; x++) {
          double cl[3] = {lo[0] + (double)x * h, lo[1] + (double)y * h, lo[2] + (double)z * h};
          double ch[3] = {cl[0] + h, cl[1] + h, cl[2] + h};
          hz_pclip_poly P;
          if (hz_pclip_tri(A, B, C, cl, ch, &P) < 3) continue;
          if (!(hz_pclip_area(&P) > 0.0)) continue;
          size_t ci = ((size_t)z * (size_t)n + (size_t)y) * (size_t)n + (size_t)x;
          cnt[ci + 1]++;
          occ[ci] = 1;
        }
  }
  for (size_t i = 0; i < ncell; i++)
    cnt[i + 1] += cnt[i];
  int32_t *lst = malloc((size_t)cnt[ncell] * sizeof *lst);
  int32_t *fill = calloc(ncell, sizeof *fill);
  if (lst == NULL || fill == NULL) return 1;
  for (int32_t t = 0; t < m.nt; t++) {
    int64_t i0[3] = {0, 0, 0}, i1[3] = {0, 0, 0};
    for (int c = 0; c < 3; c++) {
      double a = 1e300, b = -1e300;
      for (int v = 0; v < 3; v++) {
        double x = m.v[3 * (size_t)m.f[3 * (size_t)t + (size_t)v] + (size_t)c];
        if (x < a) a = x;
        if (x > b) b = x;
      }
      i0[c] = (int64_t)floor((a - lo[c]) / h);
      i1[c] = (int64_t)floor((b - lo[c]) / h);
      if (i0[c] < 0) i0[c] = 0;
      if (i1[c] >= n) i1[c] = n - 1;
    }
    const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
    const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
    const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
    for (int64_t z = i0[2]; z <= i1[2]; z++)
      for (int64_t y = i0[1]; y <= i1[1]; y++)
        for (int64_t x = i0[0]; x <= i1[0]; x++) {
          double cl[3] = {lo[0] + (double)x * h, lo[1] + (double)y * h, lo[2] + (double)z * h};
          double ch[3] = {cl[0] + h, cl[1] + h, cl[2] + h};
          hz_pclip_poly P;
          if (hz_pclip_tri(A, B, C, cl, ch, &P) < 3) continue;
          if (!(hz_pclip_area(&P) > 0.0)) continue;
          size_t ci = ((size_t)z * (size_t)n + (size_t)y) * (size_t)n + (size_t)x;
          lst[cnt[ci] + fill[ci]++] = t;
        }
  }
  int64_t nocc = 0;
  for (size_t i = 0; i < ncell; i++)
    if (occ[i]) nocc++;
  printf("   списки по ячейкам за %.2f с: занятых ячеек %lld (%.2f %%), вхождений %d\n",
         now_s() - t0, (long long)nocc, 100.0 * (double)nocc / (double)ncell, cnt[ncell]);
  free(fill);

  /* --- 2. ЗАЛИВКА ПУСТОТЫ ОТ ГРАНИЦЫ (А596: перекрывает ЗАНЯТАЯ ячейка) --- */
  t0 = now_s();
  unsigned char *out = calloc(ncell, 1);
  int32_t *stk = malloc(ncell * sizeof *stk);
  if (out == NULL || stk == NULL) return 1;
  int64_t sp = 0;
  for (int z = 0; z < n; z++)
    for (int y = 0; y < n; y++)
      for (int x = 0; x < n; x++) {
        if (x != 0 && y != 0 && z != 0 && x != n - 1 && y != n - 1 && z != n - 1) continue;
        size_t ci = ((size_t)z * (size_t)n + (size_t)y) * (size_t)n + (size_t)x;
        if (occ[ci] || out[ci]) continue;
        out[ci] = 1;
        stk[sp++] = (int32_t)ci;
      }
  while (sp > 0) {
    int32_t ci = stk[--sp];
    int x = ci % n, y = (ci / n) % n, z = ci / (n * n);
    for (int d = 0; d < 6; d++) {
      int nx = x + ((d == 0) - (d == 1)), ny = y + ((d == 2) - (d == 3)),
          nz = z + ((d == 4) - (d == 5));
      if (nx < 0 || ny < 0 || nz < 0 || nx >= n || ny >= n || nz >= n) continue;
      size_t k2 = ((size_t)nz * (size_t)n + (size_t)ny) * (size_t)n + (size_t)nx;
      if (occ[k2] || out[k2]) continue;
      out[k2] = 1;
      stk[sp++] = (int32_t)k2;
    }
  }
  int64_t nout = 0, nin = 0;
  for (size_t i = 0; i < ncell; i++) {
    if (out[i])
      nout++;
    else if (!occ[i])
      nin++;
  }
  printf("   заливка за %.2f с: снаружи %lld (%.1f %%), ВНУТРИ %lld (%.1f %%), с геометрией %lld\n",
         now_s() - t0, (long long)nout, 100.0 * (double)nout / (double)ncell, (long long)nin,
         100.0 * (double)nin / (double)ncell, (long long)nocc);
  {
    /* КУДА ПОПАЛА КАМЕРА (возражение пользователя 08-09, и оно решает всё).
     * Комната — тоже ЗАМКНУТАЯ полость: заливка от границы габарита может
     * пометить как «внутри» ровно то, на что мы смотрим. Тогда «даровое
     * удаление невидимого» окажется удалением ВИДИМОГО, и вывод перевёрнут. */
    double eye[3] = HZ_CFG_HALL_EYE;
    int64_t ic[3] = {0, 0, 0};
    int okc = 1;
    for (int c = 0; c < 3; c++) {
      ic[c] = (int64_t)floor((eye[c] - lo[c]) / h);
      if (ic[c] < 0 || ic[c] >= n) okc = 0;
    }
    if (okc) {
      size_t ci = ((size_t)ic[2] * (size_t)n + (size_t)ic[1]) * (size_t)n + (size_t)ic[0];
      printf("   ЯЧЕЙКА КАМЕРЫ §2: %s\n",
             occ[ci] ? "С ГЕОМЕТРИЕЙ"
                     : (out[ci] ? "снаружи (заливка дошла)" : "ВНУТРИ — ЗАЛИВКА НЕ ДОШЛА"));
    } else
      printf("   ЯЧЕЙКА КАМЕРЫ §2: вне сетки\n");
  }
  free(stk);

  /* Знак угла: СНАРУЖИ, если хотя бы одна смежная ячейка снаружи. Тогда у стены
   * внешние углы снаружи, внутренние внутри, и смена знака приходится на неё. */
  src_ctx S;
  S.n = n;
  S.sgn = calloc((size_t)n1 * (size_t)n1 * (size_t)n1, 1);
  if (S.sgn == NULL) {
    free(out);
    free(occ);
    free(cnt);
    free(lst);
    return 1;
  }
  int64_t nsin = 0;
  for (int z = 0; z <= n; z++)
    for (int y = 0; y <= n; y++)
      for (int x = 0; x <= n; x++) {
        int isout = 0;
        for (int dz = -1; dz <= 0 && !isout; dz++)
          for (int dy = -1; dy <= 0 && !isout; dy++)
            for (int dx = -1; dx <= 0 && !isout; dx++) {
              int cx = x + dx, cy = y + dy, cz = z + dz;
              if (cx < 0 || cy < 0 || cz < 0 || cx >= n || cy >= n || cz >= n) {
                isout = 1;
                break;
              }
              if (out[((size_t)cz * (size_t)n + (size_t)cy) * (size_t)n + (size_t)cx]) isout = 1;
            }
        int sv = noflood ? 0 : (isout ? 0 : 1);
        if (invert) sv = noflood ? 0 : !sv;
        S.sgn[((size_t)z * (size_t)n1 + (size_t)y) * (size_t)n1 + (size_t)x] = (unsigned char)sv;
        if (sv) nsin++;
      }
  printf("   углов ВНУТРИ %lld из %lld (%.1f %%)%s\n", (long long)nsin,
         (long long)((size_t)n1 * (size_t)n1 * (size_t)n1),
         100.0 * (double)nsin / (double)((size_t)n1 * (size_t)n1 * (size_t)n1),
         noflood ? "  [НЕГАТИВНЫЙ КОНТРОЛЬ: заливки нет]" : (invert ? "  [ЗНАК ОБРАЩЁН]" : ""));

  /* --- 3. эрмитовы данные на рёбрах сетки --- */
  t0 = now_s();
  size_t nc1 = (size_t)n1 * (size_t)n1 * (size_t)n1;
  S.et = malloc(3 * nc1 * sizeof *S.et);
  S.enrm = malloc(9 * nc1 * sizeof *S.enrm);
  if (S.et == NULL || S.enrm == NULL) {
    free(S.sgn);
    free(S.et);
    free(S.enrm);
    free(out);
    free(occ);
    free(cnt);
    free(lst);
    return 1;
  }
  for (size_t i = 0; i < 3 * nc1; i++)
    S.et[i] = -1.0;
  int64_t ncross = 0;
  for (int z = 0; z <= n; z++)
    for (int y = 0; y <= n; y++)
      for (int x = 0; x <= n; x++) {
        size_t k = ((size_t)z * (size_t)n1 + (size_t)y) * (size_t)n1 + (size_t)x;
        for (int ax = 0; ax < 3; ax++) {
          int e[3] = {x, y, z};
          if (e[ax] >= n) continue;
          double P0[3] = {lo[0] + (double)x * h, lo[1] + (double)y * h, lo[2] + (double)z * h};
          double P1[3] = {P0[0], P0[1], P0[2]};
          P1[ax] += h;
          /* Кандидаты — из ячеек, к которым ребро принадлежит (до четырёх). */
          double best = 2.0;
          double bn[3] = {0.0, 0.0, 0.0};
          for (int dp = -1; dp <= 0; dp++)
            for (int dq = -1; dq <= 0; dq++) {
              int c3[3] = {x, y, z};
              c3[(ax + 1) % 3] += dp;
              c3[(ax + 2) % 3] += dq;
              if (c3[0] < 0 || c3[1] < 0 || c3[2] < 0 || c3[0] >= n || c3[1] >= n || c3[2] >= n)
                continue;
              size_t ci = ((size_t)c3[2] * (size_t)n + (size_t)c3[1]) * (size_t)n + (size_t)c3[0];
              for (int32_t q = cnt[ci]; q < cnt[ci + 1]; q++) {
                int32_t t = lst[q];
                const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
                const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
                const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
                double e1[3], e2[3], pv[3], tv[3], qv[3], dir[3];
                for (int c = 0; c < 3; c++) {
                  e1[c] = B[c] - A[c];
                  e2[c] = C[c] - A[c];
                  dir[c] = P1[c] - P0[c];
                }
                pv[0] = dir[1] * e2[2] - dir[2] * e2[1];
                pv[1] = dir[2] * e2[0] - dir[0] * e2[2];
                pv[2] = dir[0] * e2[1] - dir[1] * e2[0];
                double det = e1[0] * pv[0] + e1[1] * pv[1] + e1[2] * pv[2];
                if (det > -1e-15 && det < 1e-15) continue;
                double inv = 1.0 / det;
                for (int c = 0; c < 3; c++)
                  tv[c] = P0[c] - A[c];
                double uu = (tv[0] * pv[0] + tv[1] * pv[1] + tv[2] * pv[2]) * inv;
                if (uu < 0.0 || uu > 1.0) continue;
                qv[0] = tv[1] * e1[2] - tv[2] * e1[1];
                qv[1] = tv[2] * e1[0] - tv[0] * e1[2];
                qv[2] = tv[0] * e1[1] - tv[1] * e1[0];
                double vv = (dir[0] * qv[0] + dir[1] * qv[1] + dir[2] * qv[2]) * inv;
                if (vv < 0.0 || uu + vv > 1.0) continue;
                double tt = (e2[0] * qv[0] + e2[1] * qv[1] + e2[2] * qv[2]) * inv;
                if (tt < 0.0 || tt > 1.0) continue;
                if (tt < best) {
                  best = tt;
                  double cr[3];
                  cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
                  cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
                  cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
                  double l2 = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
                  if (!(l2 > 0.0)) l2 = 1.0;
                  for (int c = 0; c < 3; c++)
                    bn[c] = cr[c] / l2;
                }
              }
            }
          if (best <= 1.0) {
            S.et[3 * k + (size_t)ax] = best;
            for (int c = 0; c < 3; c++)
              S.enrm[9 * k + 3 * (size_t)ax + (size_t)c] = bn[c];
            ncross++;
          }
        }
      }
  printf("   эрмитовы данные за %.2f с: рёбер с пересечением %lld\n", now_s() - t0,
         (long long)ncross);

  /* --- 3b. СВЕДЕНИЕ ЗНАКА И ПЕРЕСЕЧЕНИЙ (найдено прогоном: hz_dc_sample вернул
   * ETOPO). Построитель требует инварианта: пересечение есть ТОГДА И ТОЛЬКО
   * ТОГДА, когда у концов ребра разный знак. У меня же два НЕЗАВИСИМЫХ
   * источника — заливка по ячейкам и пересечение с треугольником, — и совпадать
   * они не обязаны. Расхождение печатается ЗАМЕРОМ: оно и есть мера того,
   * насколько грубая сетка расходится с геометрией. Ребро со сменой знака без
   * найденного треугольника получает заплату в середине с нормалью соседа. */
  {
    int64_t nfix = 0, ndrop = 0, nkeep = 0;
    for (int z = 0; z <= n; z++)
      for (int y = 0; y <= n; y++)
        for (int x = 0; x <= n; x++) {
          size_t k = ((size_t)z * (size_t)n1 + (size_t)y) * (size_t)n1 + (size_t)x;
          for (int ax = 0; ax < 3; ax++) {
            int e[3] = {x, y, z};
            if (e[ax] >= n) continue;
            int e2i[3] = {x, y, z};
            e2i[ax]++;
            size_t k2 =
                ((size_t)e2i[2] * (size_t)n1 + (size_t)e2i[1]) * (size_t)n1 + (size_t)e2i[0];
            int s0 = S.sgn[k], s1 = S.sgn[k2];
            int has = (S.et[3 * k + (size_t)ax] >= 0.0);
            if (s0 == s1) {
              if (has) {
                S.et[3 * k + (size_t)ax] = -1.0;
                ndrop++;
              }
            } else if (!has) {
              S.et[3 * k + (size_t)ax] = 0.5;
              double nn[3] = {0.0, 0.0, 0.0};
              nn[ax] = 1.0;
              for (int c = 0; c < 3; c++)
                S.enrm[9 * k + 3 * (size_t)ax + (size_t)c] = nn[c];
              nfix++;
            } else
              nkeep++;
          }
        }
    printf("   СВЕДЕНИЕ: пересечений при СМЕНЕ знака %lld; заплат (знак сменился, треугольника "
           "нет) %lld; отброшено (треугольник есть, знак не сменился) %lld\n",
           (long long)nkeep, (long long)nfix, (long long)ndrop);
  }

  /* --- 4. дерево DC --- */
  t0 = now_s();
  hz_signgrid g;
  hz_htab ht;
  memset(&g, 0, sizeof g);
  if (hz_htab_init(&ht) != 0) return 1;
  int rc = hz_dc_sample(&g, &ht, lev, src_sign, src_cross, &S);
  if (rc != HZ_DC_OK) {
    fprintf(stderr, "опрос источника: код %d\n", rc);
    return 1;
  }
  hz_dctree T;
  if (hz_dc_init(&T, lev) != 0) return 1;
  rc = hz_dc_build(&T, &g, &ht);
  printf("   дерево DC за %.2f с: код %d, узлов %d, рёбер в таблице %d; ЗАГНАНО %d, "
         "НЕМАНИФОЛДНЫХ %d\n",
         now_s() - t0, rc, T.n, ht.n, T.nclamped, T.nmulti);
  int64_t nvert = 0;
  for (int32_t i = 0; i < T.n; i++)
    if (T.nd[i].flags & HZ_DC_HASVERT) nvert++;
  printf("   вершин выдано %lld\n", (long long)nvert);

  /* --- 5. ошибка поверхности (ОДНОСТОРОННЯЯ, А598) --- */
  if (nvert > 0) {
    double *er = malloc((size_t)nvert * sizeof *er);
    if (er != NULL) {
      int64_t ne = 0;
      for (int32_t i = 0; i < T.n; i++) {
        if (!(T.nd[i].flags & HZ_DC_HASVERT)) continue;
        double V[3];
        for (int c = 0; c < 3; c++)
          V[c] = lo[c] + T.nd[i].vx[c] * h;
        int64_t ic[3];
        int okc = 1;
        for (int c = 0; c < 3; c++) {
          ic[c] = (int64_t)floor((V[c] - lo[c]) / h);
          if (ic[c] < 0 || ic[c] >= n) okc = 0;
        }
        if (!okc) continue;
        double bd = 1e300;
        for (int dz = -1; dz <= 1; dz++)
          for (int dy = -1; dy <= 1; dy++)
            for (int dx = -1; dx <= 1; dx++) {
              int64_t c3[3] = {ic[0] + dx, ic[1] + dy, ic[2] + dz};
              if (c3[0] < 0 || c3[1] < 0 || c3[2] < 0 || c3[0] >= n || c3[1] >= n || c3[2] >= n)
                continue;
              size_t ci = ((size_t)c3[2] * (size_t)n + (size_t)c3[1]) * (size_t)n + (size_t)c3[0];
              for (int32_t q = cnt[ci]; q < cnt[ci + 1]; q++) {
                const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)lst[q] + 0];
                const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)lst[q] + 1];
                const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)lst[q] + 2];
                double e1[3], e2[3], cr[3];
                for (int c = 0; c < 3; c++) {
                  e1[c] = B[c] - A[c];
                  e2[c] = C[c] - A[c];
                }
                cr[0] = e1[1] * e2[2] - e1[2] * e2[1];
                cr[1] = e1[2] * e2[0] - e1[0] * e2[2];
                cr[2] = e1[0] * e2[1] - e1[1] * e2[0];
                double l2 = sqrt(cr[0] * cr[0] + cr[1] * cr[1] + cr[2] * cr[2]);
                if (!(l2 > 0.0)) continue;
                double d =
                    fabs((V[0] - A[0]) * cr[0] + (V[1] - A[1]) * cr[1] + (V[2] - A[2]) * cr[2]) /
                    l2;
                if (d < bd) bd = d;
              }
            }
        if (bd < 1e299) er[ne++] = bd;
      }
      if (ne > 0) {
        qsort(er, (size_t)ne, sizeof *er, cmp_d);
        printf("   ОШИБКА ПОВЕРХНОСТИ (вершина DC → ближайшая плоскость треугольника; "
               "ОДНОСТОРОННЯЯ, А598): p50 %.4f p90 %.4f max %.4f м, по %lld вершинам\n",
               er[ne / 2], er[(ne * 9) / 10], er[ne - 1], (long long)ne);
      }
      free(er);
    }
  }

  hz_dc_free(&T);
  hz_htab_free(&ht);
  hz_signgrid_free(&g);
  free(S.sgn);
  free(S.et);
  free(S.enrm);
  free(out);
  free(occ);
  free(cnt);
  free(lst);
  hz_obj_free(&m);
  return 0;
}
