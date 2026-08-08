/* pcull.c — односторонние пометки невидимости. Разбор и оговорки — в `pcull.h`. */
#include "pcull.h"
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* ДАЛЬНОСТЬ ВЕЗДЕ СЧИТАЕТСЯ ВДОЛЬ ОСИ ВЗГЛЯДА, и это не соглашение об удобстве:
 * направление луча пикселя нормировано условием `D · fw = 1`, поэтому дальность
 * вдоль оси есть параметр вдоль ЛУЧА, и порядок «ближе — дальше» у неё тот же.
 * Только поэтому записанное в буфер сравнимо с ближней точкой ячейки.
 *
 * Проекция точки в пиксельные координаты. Возвращает `0`, если точка не перед
 * глазом (тогда числа не значат ничего и звать их нельзя). */
static int project(const hz_pcull *C, const double *X, double *sx, double *sy, double *z) {
  double v[3], a = 0.0, b = 0.0, zz = 0.0;
  for (int c = 0; c < 3; c++) {
    v[c] = X[c] - C->eye[c];
    a += v[c] * C->rt[c];
    b += v[c] * C->up[c];
    zz += v[c] * C->fw[c];
  }
  if (!(zz > C->zmin)) return 0;
  double h = 0.5 * (double)C->side;
  *sx = h * (1.0 + (a / zz) / C->tanh_);
  *sy = h * (1.0 + (b / zz) / C->tanh_);
  *z = zz;
  return 1;
}

int hz_pcull_init(hz_pcull *C, const double eye[3], const double at[3], const double up[3],
                  double fovdeg, int side, double diag) {
  memset(C, 0, sizeof *C);
  if (side < 2) return 1;
  double ln = 0.0;
  for (int c = 0; c < 3; c++) {
    C->eye[c] = eye[c];
    C->fw[c] = at[c] - eye[c];
    ln += C->fw[c] * C->fw[c];
  }
  ln = sqrt(ln);
  if (!(ln > 0.0)) return 1;
  for (int c = 0; c < 3; c++)
    C->fw[c] /= ln;
  C->rt[0] = C->fw[1] * up[2] - C->fw[2] * up[1];
  C->rt[1] = C->fw[2] * up[0] - C->fw[0] * up[2];
  C->rt[2] = C->fw[0] * up[1] - C->fw[1] * up[0];
  ln = sqrt(C->rt[0] * C->rt[0] + C->rt[1] * C->rt[1] + C->rt[2] * C->rt[2]);
  if (!(ln > 0.0)) return 1;
  for (int c = 0; c < 3; c++)
    C->rt[c] /= ln;
  C->up[0] = C->rt[1] * C->fw[2] - C->rt[2] * C->fw[1];
  C->up[1] = C->rt[2] * C->fw[0] - C->rt[0] * C->fw[2];
  C->up[2] = C->rt[0] * C->fw[1] - C->rt[1] * C->fw[0];
  C->tanh_ = tan(0.5 * fovdeg * 3.14159265358979323846 / 180.0);
  C->side = side;
  /* БЛИЖНЯЯ ПЛОСКОСТЬ — ОТНОСИТЕЛЬНАЯ, А НЕ НАЗНАЧЕННАЯ. Ближе неё проекция
   * вырастает настолько, что рёберная функция теряет разряды, на которых стоит
   * односторонность. Берётся `2⁻²⁰` диагонали сцены: при этом пиксельные
   * координаты не превосходят `10⁶`, а произведения в рёберной функции —
   * `10¹²`, то есть до конца двойной точности остаётся три порядка. */
  C->zmin = ldexp(diag > 0.0 ? diag : 1.0, -20);
  int n = 0;
  for (int s = side; s >= 1; s = (s + 1) / 2) {
    n++;
    if (s == 1) break;
  }
  C->nlev = n;
  C->lvl = calloc((size_t)n, sizeof *C->lvl);
  C->lw = calloc((size_t)n, sizeof *C->lw);
  if (C->lvl == NULL || C->lw == NULL) {
    hz_pcull_free(C);
    return 2;
  }
  int s = side;
  for (int k = 0; k < n; k++) {
    C->lw[k] = s;
    C->lvl[k] = malloc((size_t)s * (size_t)s * sizeof *C->lvl[k]);
    if (C->lvl[k] == NULL) {
      hz_pcull_free(C);
      return 2;
    }
    for (size_t i = 0; i < (size_t)s * (size_t)s; i++)
      C->lvl[k][i] = INFINITY;
    s = (s + 1) / 2;
  }
  return 0;
}

void hz_pcull_free(hz_pcull *C) {
  if (C->lvl != NULL)
    for (int k = 0; k < C->nlev; k++)
      free(C->lvl[k]);
  free(C->lvl);
  free(C->lw);
  C->lvl = NULL;
  C->lw = NULL;
  C->nlev = 0;
}

/* Дальность до плоскости треугольника вдоль луча пиксельной точки `(px, py)`.
 * `num`, `nr`, `nu`, `nf` — числитель и разложение нормали по раме камеры.
 * Возвращает `0`, если знаменатель не положителен: тогда точка луча плоскости не
 * встречает впереди, и пиксель считать нельзя. */
static int plane_z(const hz_pcull *C, double px, double py, double num, double nr, double nu,
                   double nf, double *z) {
  double h = 0.5 * (double)C->side;
  double as = C->tanh_ * (px / h - 1.0), bs = C->tanh_ * (py / h - 1.0);
  double den = nf + as * nr + bs * nu;
  if (!(den > 0.0)) return 0;
  double zz = num / den;
  if (!(zz > 0.0) || !(zz < HUGE_VAL)) return 0;
  *z = zz;
  return 1;
}

void hz_pcull_draw(hz_pcull *C, const hz_objmesh *m) {
  float *Z = C->lvl[0];
  const int S = C->side;
  for (int32_t t = 0; t < m->nt; t++) {
    const double *V[3];
    for (int i = 0; i < 3; i++)
      V[i] = m->v + 3 * (size_t)m->f[3 * (size_t)t + (size_t)i];
    double sx[3], sy[3], sz[3];
    int ok = 1;
    for (int i = 0; i < 3 && ok; i++)
      ok = project(C, V[i], &sx[i], &sy[i], &sz[i]);
    if (!ok) {
      C->ntri_near++; /* пересекает ближнюю плоскость — пропускается ЦЕЛИКОМ */
      continue;
    }
    double xmin = sx[0], xmax = sx[0], ymin = sy[0], ymax = sy[0];
    for (int i = 1; i < 3; i++) {
      if (sx[i] < xmin) xmin = sx[i];
      if (sx[i] > xmax) xmax = sx[i];
      if (sy[i] < ymin) ymin = sy[i];
      if (sy[i] > ymax) ymax = sy[i];
    }
    if (!(xmax > 0.0) || !(ymax > 0.0) || !(xmin < (double)S) || !(ymin < (double)S)) {
      C->ntri_off++;
      continue;
    }
    int i0 = (int)floor(xmin), i1 = (int)ceil(xmax), j0 = (int)floor(ymin), j1 = (int)ceil(ymax);
    if (i0 < 0) i0 = 0;
    if (j0 < 0) j0 = 0;
    if (i1 > S) i1 = S;
    if (j1 > S) j1 = S;
    if (i1 - i0 < 1 || j1 - j0 < 1) continue; /* целого пикселя не накрыть */
    /* Ориентация и нормаль в раме камеры. */
    double a2 = (sx[1] - sx[0]) * (sy[2] - sy[0]) - (sy[1] - sy[0]) * (sx[2] - sx[0]);
    if (!(a2 > 0.0) && !(a2 < 0.0)) continue;
    double sg = (a2 > 0.0) ? 1.0 : -1.0;
    double e1[3], e2[3], nw[3];
    for (int c = 0; c < 3; c++) {
      e1[c] = V[1][c] - V[0][c];
      e2[c] = V[2][c] - V[0][c];
    }
    nw[0] = e1[1] * e2[2] - e1[2] * e2[1];
    nw[1] = e1[2] * e2[0] - e1[0] * e2[2];
    nw[2] = e1[0] * e2[1] - e1[1] * e2[0];
    double num = 0.0;
    for (int c = 0; c < 3; c++)
      num += nw[c] * (V[0][c] - C->eye[c]);
    if (num < 0.0) {
      num = -num;
      for (int c = 0; c < 3; c++)
        nw[c] = -nw[c];
    }
    if (!(num > 0.0)) continue; /* плоскость проходит через глаз */
    double nr = 0.0, nu = 0.0, nf = 0.0;
    for (int c = 0; c < 3; c++) {
      nr += nw[c] * C->rt[c];
      nu += nw[c] * C->up[c];
      nf += nw[c] * C->fw[c];
    }
    /* Запас на округление рёберной функции: масштаб координат × длина ребра. */
    double scale = fabs(xmin) + fabs(xmax) + fabs(ymin) + fabs(ymax);
    for (int j = j0; j < j1; j++) {
      for (int i = i0; i < i1; i++) {
        double zz;
        if (C->loose == 1) {
          /* НЕГАТИВНЫЙ КОНТРОЛЬ: обычный z-буфер по середине пикселя. */
          double cx = (double)i + 0.5, cy = (double)j + 0.5;
          int in = 1;
          for (int k = 0; k < 3 && in; k++) {
            int k2 = (k + 1) % 3;
            double ex = sx[k2] - sx[k], ey = sy[k2] - sy[k];
            double ev = (cx - sx[k]) * ey - (cy - sy[k]) * ex;
            if (!(sg * ev < 0.0)) in = 0;
          }
          if (!in) continue;
          if (!plane_z(C, cx, cy, num, nr, nu, nf, &zz)) continue;
        } else {
          /* НАКРЫТ ЛИ ПИКСЕЛЬ ЦЕЛИКОМ: все четыре угла внутри всех трёх рёбер, и
           * с запасом. Иначе луч, проходящий сквозь непокрытый угол, ничем не
           * заслонён, а буфер объявил бы его заслонённым. */
          int full = 1;
          for (int k = 0; k < 3 && full; k++) {
            int k2 = (k + 1) % 3;
            double ex = sx[k2] - sx[k], ey = sy[k2] - sy[k];
            double tol = HZ_PCULL_ULP * DBL_EPSILON * (fabs(ex) + fabs(ey)) * scale;
            for (int q = 0; q < 4 && full; q++) {
              double cx = (double)(i + (q & 1)), cy = (double)(j + ((q >> 1) & 1));
              double ev = (cx - sx[k]) * ey - (cy - sy[k]) * ex;
              if (!(sg * ev < -tol)) full = 0;
            }
          }
          if (!full) continue;
          /* ДАЛЬНЯЯ точка треугольника в этом пикселе. Дальность есть отношение
           * постоянной к аффинной функции экранных координат, а такая функция
           * достигает крайних значений на выпуклом множестве В ЕГО ВЕРШИНАХ —
           * поэтому четырёх углов достаточно, а не «достаточно приближённо». */
          zz = (C->loose == 2) ? HUGE_VAL : 0.0;
          int good = 1;
          for (int q = 0; q < 4; q++) {
            double cz;
            if (!plane_z(C, (double)(i + (q & 1)), (double)(j + ((q >> 1) & 1)), num, nr, nu, nf,
                         &cz)) {
              good = 0;
              break;
            }
            if (C->loose == 2 ? (cz < zz) : (cz > zz)) zz = cz;
          }
          if (!good) continue;
        }
        /* Округление ВВЕРХ при укладке во `float`: записанное обязано быть не
         * меньше вычисленного, иначе часть заслонённого объёма окажется «дальше
         * буфера» лишь из-за усечения мантиссы. */
        float fz = nextafterf((float)zz, INFINITY);
        size_t idx = (size_t)j * (size_t)S + (size_t)i;
        if (fz < Z[idx]) {
          Z[idx] = fz;
          C->npix++;
        }
      }
    }
  }
}

void hz_pcull_pyramid(hz_pcull *C) {
  for (int k = 1; k < C->nlev; k++) {
    const float *A = C->lvl[k - 1];
    float *B = C->lvl[k];
    int wa = C->lw[k - 1], wb = C->lw[k];
    for (int j = 0; j < wb; j++)
      for (int i = 0; i < wb; i++) {
        /* МАКСИМУМ, а не минимум: грубый уровень отвечает на вопрос «что здесь
         * заведомо не заслонено», и отсутствующий (за краем) отпрыск читается
         * как `+∞` — то есть незаслонённый. Обе оговорки в сторону «пометок
         * меньше». */
        float mx = -INFINITY;
        for (int q = 0; q < 4; q++) {
          int ii = 2 * i + (q & 1), jj = 2 * j + ((q >> 1) & 1);
          float v = (ii < wa && jj < wa) ? A[(size_t)jj * (size_t)wa + (size_t)ii] : INFINITY;
          if (v > mx) mx = v;
        }
        B[(size_t)j * (size_t)wb + (size_t)i] = mx;
      }
  }
}

/* Максимум буфера по прямоугольнику пикселей — с ЗАПАСОМ: берётся уровень, на
 * котором прямоугольник укладывается в две клетки по стороне, и максимум по ним.
 * Область при этом НАКРЫВАЕТСЯ с избытком, максимум выходит завышенным, и
 * ошибка снова идёт в сторону «пометок меньше». */
static float rect_max(const hz_pcull *C, int i0, int i1, int j0, int j1) {
  int L = 0;
  while (L + 1 < C->nlev && ((i1 >> L) - (i0 >> L) > 1 || (j1 >> L) - (j0 >> L) > 1))
    L++;
  const float *A = C->lvl[L];
  int w = C->lw[L];
  float mx = -INFINITY;
  for (int j = j0 >> L; j <= (j1 >> L); j++)
    for (int i = i0 >> L; i <= (i1 >> L); i++) {
      if (i < 0 || j < 0 || i >= w || j >= w) return INFINITY;
      float v = A[(size_t)j * (size_t)w + (size_t)i];
      if (v > mx) mx = v;
    }
  return mx;
}

int hz_pcull_hidden(hz_pcull *C, const double lo[3], const double hi[3]) {
  C->ntest++;
  /* ВНЕ КАДРА — ЭТО ТОЖЕ «НЕ ВИДНО», И ОТВЕЧАТЬ НА ЭТО НАДО ЗДЕСЬ. Раньше такие
   * куски пространства уходили в отказ («за краем буфера», «у ближней
   * плоскости»), и на Сан-Мигеле отказом кончались девять проверок из десяти:
   * `26.8` млн узлов у ближней плоскости и `10.7` млн за краем буфера из `41.8`
   * млн. Для пола фронта это было безразлично (§279), для КЛАССИФИКАЦИИ
   * ПРОСТРАНСТВА — решающе.
   *
   * Проверка обычная и консервативная: коробка целиком снаружи хотя бы одной
   * плоскости. Ближняя плоскость входит в тот же набор, поэтому «позади глаза»
   * ловится ею же, а не отдельным правилом. */
  if (C->usefr) {
    for (int k = 0; k < 6; k++) {
      const double *P = C->fr[k];
      double s = P[3];
      for (int c = 0; c < 3; c++)
        s += (P[c] > 0.0 ? hi[c] : lo[c]) * P[c];
      if (s < 0.0) {
        C->nfrust++;
        return 2;
      }
    }
  }
  /* Глаз внутри коробки — о заслонении речи нет. */
  int inside = 1;
  for (int c = 0; c < 3 && inside; c++)
    if (C->eye[c] < lo[c] || C->eye[c] > hi[c]) inside = 0;
  if (inside) return 0;
  double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300, znear = 1e300;
  for (int k = 0; k < 8; k++) {
    double P[3] = {(k & 1) ? hi[0] : lo[0], (k & 2) ? hi[1] : lo[1], (k & 4) ? hi[2] : lo[2]};
    double px, py, pz;
    if (!project(C, P, &px, &py, &pz)) {
      C->nnear++;
      return 0; /* хоть один угол не перед глазом — не помечаем */
    }
    if (px < xmin) xmin = px;
    if (px > xmax) xmax = px;
    if (py < ymin) ymin = py;
    if (py > ymax) ymax = py;
    if (pz < znear) znear = pz;
  }
  /* Проекция коробки есть выпуклая оболочка проекций её углов (все они перед
   * глазом), поэтому габарит по углам накрывает её целиком. Отступ в пиксель —
   * запас на округление.
   *
   * ЗА КРАЕМ БУФЕРА — НЕ «НЕИЗВЕСТНО», А «НЕ ВИДНО». Сперва тут стоял отказ, и
   * он был лишней осторожностью: буфер накрывает РОВНО поле кадра, значит
   * пиксели за его краем суть направления ВНЕ КАДРА, а вне кадра видимого нет.
   * Поэтому прямоугольник запроса ОБРЕЗАЕТСЯ по буферу, и решается вопрос о той
   * части куска пространства, которая в кадр попадает. Если не попадает ничего —
   * кусок целиком вне кадра, и это ловится пирамидой выше. */
  int i0 = (int)floor(xmin) - 1, i1 = (int)ceil(xmax), j0 = (int)floor(ymin) - 1,
      j1 = (int)ceil(ymax);
  if (i0 < 0) i0 = 0;
  if (j0 < 0) j0 = 0;
  if (i1 > C->side - 1) i1 = C->side - 1;
  if (j1 > C->side - 1) j1 = C->side - 1;
  if (i0 > i1 || j0 > j1) {
    C->noff++;
    return 0;
  }
  float q = rect_max(C, i0, i1, j0, j1);
  double zn = znear * (1.0 - HZ_PCULL_ULP * DBL_EPSILON);
  if ((double)q < zn) {
    C->nhidden++;
    return 1;
  }
  return 0;
}

static void marks_walk(hz_pcull *C, const hz_ptree *T, int32_t nid, int maxlev,
                       unsigned char *mark) {
  int h = hz_pcull_hidden(C, T->nd[nid].lo, T->nd[nid].hi);
  if (h != 0) {
    mark[nid] = (unsigned char)h;
    return;
  }
  if (T->nd[nid].child < 0 || (int)T->lev[nid] >= maxlev) return;
  for (int k = 0; k < 8; k++)
    marks_walk(C, T, T->nd[nid].child + k, maxlev, mark);
}

void hz_pcull_marks(hz_pcull *C, const hz_ptree *T, int maxlev, unsigned char *mark) {
  if (T->nnd > 0) marks_walk(C, T, 0, maxlev, mark);
}

/* --- КАМЕРА: ОБЫЧНЫЙ z-БУФЕР С НОМЕРОМ ТРЕУГОЛЬНИКА ------------------------ */
/* Разбор — в `pcull.h`. Здесь важно одно: правила осторожности, на которых стоит
 * односторонность пометок, тут НЕ ДЕЙСТВУЮТ и действовать не должны. */

void hz_pcull_ray(const hz_pcull *C, int i, int j, double d[3]) {
  double h = 0.5 * (double)C->side;
  double as = C->tanh_ * (((double)i + 0.5) / h - 1.0);
  double bs = C->tanh_ * (((double)j + 0.5) / h - 1.0);
  for (int c = 0; c < 3; c++)
    d[c] = C->fw[c] + as * C->rt[c] + bs * C->up[c];
}

void hz_pcull_shot(hz_pcull *C, const hz_objmesh *m, float *z, int32_t *id) {
  const int S = C->side;
  for (size_t p = 0; p < (size_t)S * (size_t)S; p++) {
    z[p] = INFINITY;
    id[p] = -1;
  }
  for (int32_t t = 0; t < m->nt; t++) {
    const double *V[3];
    for (int i = 0; i < 3; i++)
      V[i] = m->v + 3 * (size_t)m->f[3 * (size_t)t + (size_t)i];
    double sx[3], sy[3], sz[3];
    int ok = 1;
    for (int i = 0; i < 3 && ok; i++)
      ok = project(C, V[i], &sx[i], &sy[i], &sz[i]);
    if (!ok) continue; /* пересекает ближнюю плоскость — не рисуется */
    double xmin = sx[0], xmax = sx[0], ymin = sy[0], ymax = sy[0];
    for (int i = 1; i < 3; i++) {
      if (sx[i] < xmin) xmin = sx[i];
      if (sx[i] > xmax) xmax = sx[i];
      if (sy[i] < ymin) ymin = sy[i];
      if (sy[i] > ymax) ymax = sy[i];
    }
    if (!(xmax > 0.0) || !(ymax > 0.0) || !(xmin < (double)S) || !(ymin < (double)S)) continue;
    int i0 = (int)floor(xmin), i1 = (int)ceil(xmax), j0 = (int)floor(ymin), j1 = (int)ceil(ymax);
    if (i0 < 0) i0 = 0;
    if (j0 < 0) j0 = 0;
    if (i1 > S) i1 = S;
    if (j1 > S) j1 = S;
    double a2 = (sx[1] - sx[0]) * (sy[2] - sy[0]) - (sy[1] - sy[0]) * (sx[2] - sx[0]);
    if (!(a2 > 0.0) && !(a2 < 0.0)) continue;
    double sg = (a2 > 0.0) ? 1.0 : -1.0;
    double e1[3], e2[3], nw[3];
    for (int c = 0; c < 3; c++) {
      e1[c] = V[1][c] - V[0][c];
      e2[c] = V[2][c] - V[0][c];
    }
    nw[0] = e1[1] * e2[2] - e1[2] * e2[1];
    nw[1] = e1[2] * e2[0] - e1[0] * e2[2];
    nw[2] = e1[0] * e2[1] - e1[1] * e2[0];
    double num = 0.0;
    for (int c = 0; c < 3; c++)
      num += nw[c] * (V[0][c] - C->eye[c]);
    if (num < 0.0) {
      num = -num;
      for (int c = 0; c < 3; c++)
        nw[c] = -nw[c];
    }
    if (!(num > 0.0)) continue;
    double nr = 0.0, nu = 0.0, nf = 0.0;
    for (int c = 0; c < 3; c++) {
      nr += nw[c] * C->rt[c];
      nu += nw[c] * C->up[c];
      nf += nw[c] * C->fw[c];
    }
    for (int j = j0; j < j1; j++)
      for (int i = i0; i < i1; i++) {
        double cx = (double)i + 0.5, cy = (double)j + 0.5;
        int in = 1;
        for (int k = 0; k < 3 && in; k++) {
          int k2 = (k + 1) % 3;
          double ex = sx[k2] - sx[k], ey = sy[k2] - sy[k];
          double ev = (cx - sx[k]) * ey - (cy - sy[k]) * ex;
          if (!(sg * ev < 0.0)) in = 0;
        }
        if (!in) continue;
        double zz;
        if (!plane_z(C, cx, cy, num, nr, nu, nf, &zz)) continue;
        size_t p = (size_t)j * (size_t)S + (size_t)i;
        if ((float)zz < z[p]) {
          z[p] = (float)zz;
          id[p] = t;
        }
      }
  }
}

/* --- РИСОВАНИЕ ПО ДЕРЕВУ, СПЕРЕДИ НАЗАД (§316) ------------------------------ */

typedef struct {
  hz_pcull *C;
  const hz_objmesh *m;
  const hz_ptree *T;
  const double *tlo, *thi;
  float *z;
  int32_t *id;
  unsigned char *done; /* закрыт ли пиксель — один байт вместо глубины */
  int64_t nopen;       /* сколько пикселей ещё открыто: ранний выход */
  int64_t ntri, ncell;
} shot_ctx;

/* Нарисовать один треугольник в уже закрытый-незакрытый буфер: пишем только в
 * ОТКРЫТЫЕ пиксели, потому что всё, что впереди, уже нарисовано. */
static void shot_tri(shot_ctx *S, int32_t t) {
  hz_pcull *C = S->C;
  const hz_objmesh *m = S->m;
  const int SD = C->side;
  const double *V[3];
  for (int i = 0; i < 3; i++)
    V[i] = m->v + 3 * (size_t)m->f[3 * (size_t)t + (size_t)i];
  /* Обнуление явное: `project` может вернуть отказ, и анализатор прав, что
   * доказательства заполненности у него нет (тот же класс, что в `psun`). */
  double sx[3] = {0.0, 0.0, 0.0}, sy[3] = {0.0, 0.0, 0.0}, sz[3] = {0.0, 0.0, 0.0};
  for (int i = 0; i < 3; i++)
    if (!project(C, V[i], &sx[i], &sy[i], &sz[i])) return;
  double xmin = sx[0], xmax = sx[0], ymin = sy[0], ymax = sy[0];
  for (int i = 1; i < 3; i++) {
    if (sx[i] < xmin) xmin = sx[i];
    if (sx[i] > xmax) xmax = sx[i];
    if (sy[i] < ymin) ymin = sy[i];
    if (sy[i] > ymax) ymax = sy[i];
  }
  if (!(xmax > 0.0) || !(ymax > 0.0) || !(xmin < (double)SD) || !(ymin < (double)SD)) return;
  int i0 = (int)floor(xmin), i1 = (int)ceil(xmax), j0 = (int)floor(ymin), j1 = (int)ceil(ymax);
  if (i0 < 0) i0 = 0;
  if (j0 < 0) j0 = 0;
  if (i1 > SD) i1 = SD;
  if (j1 > SD) j1 = SD;
  double a2 = (sx[1] - sx[0]) * (sy[2] - sy[0]) - (sy[1] - sy[0]) * (sx[2] - sx[0]);
  if (!(a2 > 0.0) && !(a2 < 0.0)) return;
  double sg = (a2 > 0.0) ? 1.0 : -1.0;
  double e1[3], e2[3], nw[3];
  for (int c = 0; c < 3; c++) {
    e1[c] = V[1][c] - V[0][c];
    e2[c] = V[2][c] - V[0][c];
  }
  nw[0] = e1[1] * e2[2] - e1[2] * e2[1];
  nw[1] = e1[2] * e2[0] - e1[0] * e2[2];
  nw[2] = e1[0] * e2[1] - e1[1] * e2[0];
  double num = 0.0;
  for (int c = 0; c < 3; c++)
    num += nw[c] * (V[0][c] - C->eye[c]);
  if (num < 0.0) {
    num = -num;
    for (int c = 0; c < 3; c++)
      nw[c] = -nw[c];
  }
  if (!(num > 0.0)) return;
  double nr = 0.0, nu = 0.0, nf = 0.0;
  for (int c = 0; c < 3; c++) {
    nr += nw[c] * C->rt[c];
    nu += nw[c] * C->up[c];
    nf += nw[c] * C->fw[c];
  }
  S->ntri++;
  for (int j = j0; j < j1; j++)
    for (int i = i0; i < i1; i++) {
      size_t p = (size_t)j * (size_t)SD + (size_t)i;
      if (S->done[p]) continue; /* впереди уже нарисовано — это и есть порядок */
      double cx = (double)i + 0.5, cy = (double)j + 0.5;
      int in = 1;
      for (int k = 0; k < 3 && in; k++) {
        int k2 = (k + 1) % 3;
        double ex = sx[k2] - sx[k], ey = sy[k2] - sy[k];
        double ev = (cx - sx[k]) * ey - (cy - sy[k]) * ex;
        if (!(sg * ev < 0.0)) in = 0;
      }
      if (!in) continue;
      double zz;
      if (!plane_z(C, cx, cy, num, nr, nu, nf, &zz)) continue;
      S->z[p] = (float)zz;
      S->id[p] = t;
      S->done[p] = 1;
      S->nopen--;
    }
}

/* Обход узла: дети в порядке БЛИЖНИЙ-К-ГЛАЗУ ПЕРВЫМ. */
static void shot_walk(shot_ctx *S, int32_t nid, const int32_t *list, int32_t n) {
  if (S->nopen <= 0 || n <= 0) return;
  const hz_ptnode *N = &S->T->nd[nid];
  /* Ячейка целиком вне кадра — пропускается; это не отсечение света, а лишь
   * рисование, и вне кадра рисовать нечего. */
  if (S->C->usefr) {
    for (int k = 0; k < 6; k++) {
      const double *P = S->C->fr[k];
      double s = P[3];
      for (int c = 0; c < 3; c++)
        s += (P[c] > 0.0 ? N->hi[c] : N->lo[c]) * P[c];
      if (s < 0.0) return;
    }
  }
  S->ncell++;
  if (N->child < 0) {
    for (int32_t i = 0; i < n; i++)
      shot_tri(S, list[i]);
    return;
  }
  double mid[3];
  for (int c = 0; c < 3; c++)
    mid[c] = 0.5 * (N->lo[c] + N->hi[c]);
  /* Ближний октант задаётся положением ГЛАЗА относительно середины: это и есть
   * «порядок октантов от источника» (§241.7), только источник здесь — камера. */
  int e[3];
  for (int c = 0; c < 3; c++)
    e[c] = (S->C->eye[c] >= mid[c]) ? 1 : 0;
  int32_t *sub = malloc(((size_t)n + 1) * sizeof *sub);
  if (sub == NULL) return;
  for (int far = 0; far <= 3 && S->nopen > 0; far++)
    for (int k = 0; k < 8 && S->nopen > 0; k++) {
      int b[3], nf = 0;
      for (int c = 0; c < 3; c++) {
        b[c] = (k >> c) & 1;
        if (b[c] != e[c]) nf++;
      }
      if (nf != far) continue;
      int32_t ch = N->child + k;
      const hz_ptnode *Cn = &S->T->nd[ch];
      int32_t ns = 0;
      for (int32_t i = 0; i < n; i++) {
        const double *bl = S->tlo + 3 * (size_t)list[i], *bh = S->thi + 3 * (size_t)list[i];
        int hit = 1;
        for (int c = 0; c < 3 && hit; c++)
          if (!(bl[c] < Cn->hi[c]) || !(bh[c] >= Cn->lo[c])) hit = 0;
        if (hit) sub[ns++] = list[i];
      }
      shot_walk(S, ch, sub, ns);
    }
  free(sub);
}

void hz_pcull_shot_tree(hz_pcull *C, const hz_objmesh *m, const hz_ptree *T, const double *tlo,
                        const double *thi, float *z, int32_t *id, int64_t *ntri_done,
                        int64_t *ncell_done) {
  const int SD = C->side;
  size_t np = (size_t)SD * (size_t)SD;
  unsigned char *done = calloc(np, 1);
  int32_t *list = malloc((size_t)m->nt * sizeof *list);
  if (done == NULL || list == NULL) {
    free(done);
    free(list);
    return;
  }
  for (size_t p = 0; p < np; p++) {
    z[p] = INFINITY;
    id[p] = -1;
  }
  for (int32_t t = 0; t < m->nt; t++)
    list[t] = t;
  shot_ctx S;
  memset(&S, 0, sizeof S);
  S.C = C;
  S.m = m;
  S.T = T;
  S.tlo = tlo;
  S.thi = thi;
  S.z = z;
  S.id = id;
  S.done = done;
  S.nopen = (int64_t)np;
  shot_walk(&S, 0, list, m->nt);
  if (ntri_done != NULL) *ntri_done = S.ntri;
  if (ncell_done != NULL) *ncell_done = S.ncell;
  free(done);
  free(list);
}
