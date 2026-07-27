/* PLAN_TRANSPORT.md, шаги C и D. Трёхмерная развёртка: топологический порядок,
 * DG1 с противопотоком, итерация по рассеянию, ограничитель положительности. */

#include "transport/sweep3.h"
#include "transport/tet3.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Базис ячейки: {1, ξ, η, ζ}, ξ = (x_ед − c_ед)/s. Координаты ЕДИНИЧНЫЕ, потому
 * что коробка есть куб именно в них; в мир всё переносится множителями кадра. */
static void cell_basis(const tr3_mesh *m, int32_t c, const double xw[3], double b[4]) {
  b[0] = 1.0;
  double s = (double)m->csize[c];
  for (int k = 0; k < 3; k++) {
    double xu = (xw[k] - m->fr.o[k]) / m->fr.u[k];
    b[k + 1] = (xu - ((double)m->clo[c][k] + 0.5 * s)) / s;
  }
}

/* ∫ b^A_i b^B_j dA по многоугольнику в МИРЕ: базисы РАЗНЫХ ячеек (на стыке
 * уровней они разные), поэтому матрица не симметрична и считается смешанной. */
static void face_mass2(const tr3_mesh *m, const double (*v)[4][3], int32_t ca, int32_t cb,
                       double out[4][4]) {
  memset(out, 0, 16 * sizeof(double));
  const double (*vv)[3] = *v;
  for (int e = 1; e + 1 < 4; e++) {
    double tri[3][3];
    for (int k = 0; k < 3; k++) {
      tri[0][k] = vv[0][k];
      tri[1][k] = vv[e][k];
      tri[2][k] = vv[e + 1][k];
    }
    double area = tr3_poly_face_area(tri, 3);
    double ba[3][4], bb[3][4];
    for (int i = 0; i < 3; i++) {
      cell_basis(m, ca, tri[i], ba[i]);
      if (cb >= 0)
        cell_basis(m, cb, tri[i], bb[i]);
      else
        memcpy(bb[i], ba[i], sizeof ba[i]);
    }
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        double s1 = 0.0, sf = 0.0, sg = 0.0;
        for (int q = 0; q < 3; q++) {
          s1 += ba[q][i] * bb[q][j];
          sf += ba[q][i];
          sg += bb[q][j];
        }
        out[i][j] += area * (s1 + sf * sg) / 12.0;
      }
  }
}

/* Решение 4x4 с частичным выбором. Возврат 1 при вырождении. */
static int solve4(double a[4][4], double b[4], double x[4]) {
  for (int k = 0; k < 4; k++) {
    int best = k;
    for (int i = k + 1; i < 4; i++)
      if (fabs(a[i][k]) > fabs(a[best][k])) best = i;
    if (!(fabs(a[best][k]) > 0.0)) return 1;
    if (best != k) {
      for (int j = 0; j < 4; j++) {
        double t = a[k][j];
        a[k][j] = a[best][j];
        a[best][j] = t;
      }
      double t = b[k];
      b[k] = b[best];
      b[best] = t;
    }
    for (int i = k + 1; i < 4; i++) {
      double f = a[i][k] / a[k][k];
      for (int j = k; j < 4; j++)
        a[i][j] -= f * a[k][j];
      b[i] -= f * b[k];
    }
  }
  for (int i = 3; i >= 0; i--) {
    double s = b[i];
    for (int j = i + 1; j < 4; j++)
      s -= a[i][j] * x[j];
    x[i] = s / a[i][i];
  }
  return 0;
}

/* Минимум линейной функции по УГЛАМ ячейки: у DG1 экстремум всегда в углу. */
static double corner_min(const double c[4]) {
  double mn = 1e300;
  for (int k = 0; k < 8; k++) {
    double v = c[0];
    for (int a = 0; a < 3; a++)
      v += c[a + 1] * (((k >> a) & 1) ? 0.5 : -0.5);
    if (v < mn) mn = v;
  }
  return mn;
}

int tr3_sweep_solve(const tr3_problem *p, int maxit, double tol, double *phi, tr3_stats *st) {
  const tr3_mesh *m = p->m;
  const tr3_dirs *d = p->d;
  int32_t nc = m->ncell, nd = d->n;
  double *L = calloc((size_t)nc * 4, sizeof(double));
  double *phin = calloc((size_t)nc * 4, sizeof(double));
  int32_t *indeg = calloc((size_t)nc, sizeof(int32_t));
  int32_t *order = calloc((size_t)nc, sizeof(int32_t));
  int32_t *queue = calloc((size_t)nc, sizeof(int32_t));
  double (*fm)[4][4] = calloc((size_t)m->nf, sizeof(double[4][4]));
  double (*fmb)[4][4] = calloc((size_t)m->nf, sizeof(double[4][4]));
  double (*fmx)[4][4] = calloc((size_t)m->nf, sizeof(double[4][4]));
  double *farea = calloc((size_t)m->nf, sizeof(double));
  double *bout = calloc((size_t)m->nf, sizeof(double));
  double *binf = calloc((size_t)m->nf, sizeof(double));
  if (L == NULL || phin == NULL || indeg == NULL || order == NULL || queue == NULL || fm == NULL ||
      fmx == NULL || farea == NULL) {
    free(L);
    free(phin);
    free(indeg);
    free(order);
    free(queue);
    free(fm);
    free(fmb);
    free(fmx);
    free(farea);
    free(bout);
    free(binf);
    return 1;
  }

  /* Матрицы граней считаются ОДИН раз: они от направления не зависят. */
  for (int32_t f = 0; f < m->nf; f++) {
    double v[4][3];
    tr3_face_corners(m, f, v);
    farea[f] = tr3_face_area(m, f);
    face_mass2(m, (const double (*)[4][3]) & v, m->f[f].ca, m->f[f].ca, fm[f]);
    /* И ОТДЕЛЬНО матрица В БАЗИСЕ ЯЧЕЙКИ cb: у выточного члена стоит СВОЙ базис,
     * а на стыке уровней он у соседей разный. Одна матрица на обе стороны была
     * ошибкой, и на равномерной сетке она бы не проявилась вовсе. */
    if (m->f[f].cb >= 0)
      face_mass2(m, (const double (*)[4][3]) & v, m->f[f].cb, m->f[f].cb, fmb[f]);
    face_mass2(m, (const double (*)[4][3]) & v, m->f[f].ca, m->f[f].cb, fmx[f]);
  }

  memset(phi, 0, (size_t)nc * 4 * sizeof(double));
  /* ПОЛУСФЕРНЫЙ ПОТОК СЧИТАЕТСЯ ПО САМОМУ НАБОРУ, А НЕ БЕРЁТСЯ РАВНЫМ π.
   * В континууме ∫_{ω·N>0}(ω·N)dω = π точно, а на дискретном наборе — нет: у
   * боковых стенок под интегралом стоят sqrt(1−μ²) и cos φ, которые ни гауссова
   * квадратура по μ, ни квадрантная по азимуту точно не берут. Делить на π
   * значило бы терять энергию на каждом отражении: измерено 1.3%% на замкнутой
   * полости, где точный ответ известен. Нормировка на СОБСТВЕННУЮ сумму делает
   * диффузное отражение энергетически точным на любом наборе — дословно тот же
   * приём, что p0 = 1/Σw для фазовой функции в T2. */
  double hsum[6];
  for (int wl = 0; wl < 6; wl++) {
    int ax = wl / 2;
    double sgnw = (wl & 1) ? 1.0 : -1.0;
    hsum[wl] = 0.0;
    for (int mm = 0; mm < nd; mm++) {
      double on = (ax == 0 ? d->ox[mm] : (ax == 1 ? d->oy[mm] : d->oz[mm])) * sgnw;
      if (on > 0.0) hsum[wl] += d->w[mm] * on;
    }
  }
  if (p->wall_rho != NULL)
    for (int32_t f = 0; f < m->nf; f++)
      if (m->f[f].cb < 0) bout[f] = p->wall_emit != NULL ? p->wall_emit[(int)(~m->f[f].cb)] : 0.0;
  st->nclip = 0;
  int it = 0;
  double resid = 0.0;
  for (it = 0; it < maxit; it++) {
    memset(phin, 0, (size_t)nc * 4 * sizeof(double));
    memset(binf, 0, (size_t)m->nf * sizeof(double));
    st->pin = st->pout = st->pabs = 0.0;

    for (int mm = 0; mm < nd; mm++) {
      double om[3] = {d->ox[mm], d->oy[mm], d->oz[mm]};
      /* --- топологический порядок для этого направления --- */
      memset(indeg, 0, (size_t)nc * sizeof(int32_t));
      for (int32_t f = 0; f < m->nf; f++) {
        if (m->f[f].cb < 0) continue;
        double on = om[m->f[f].axis]; /* нормаль ca->cb есть +axis */
        if (on > 0.0)
          indeg[m->f[f].cb]++;
        else if (on < 0.0)
          indeg[m->f[f].ca]++;
      }
      int32_t qh = 0, qt = 0, no = 0;
      for (int32_t c = 0; c < nc; c++)
        if (indeg[c] == 0) queue[qt++] = c;
      while (qh < qt) {
        int32_t c = queue[qh++];
        order[no++] = c;
        for (int32_t k = m->fstart[c]; k < m->fstart[c + 1]; k++) {
          int32_t f = m->flist[k];
          if (m->f[f].cb < 0) continue;
          double on = om[m->f[f].axis];
          int32_t down = -1;
          if (on > 0.0 && m->f[f].ca == c)
            down = m->f[f].cb;
          else if (on < 0.0 && m->f[f].cb == c)
            down = m->f[f].ca;
          if (down >= 0 && --indeg[down] == 0) queue[qt++] = down;
        }
      }
      if (no != nc) { /* цикл обхода — не считать молча */
        free(L);
        free(phin);
        free(indeg);
        free(order);
        free(queue);
        free(fm);
        free(fmb);
        free(fmx);
        free(farea);
        free(bout);
        free(binf);
        return 2;
      }

      /* --- проход --- */
      for (int32_t oi = 0; oi < no; oi++) {
        int32_t c = order[oi];
        double A[4][4], rhs[4];
        memset(A, 0, sizeof A);
        memset(rhs, 0, sizeof rhs);
        double s = (double)m->csize[c];
        double vol = s * s * s * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
        double M[4] = {vol, vol / 12.0, vol / 12.0, vol / 12.0};

        /* объёмные члены: σ_t·∫b_i b_j и −∫ b_i (ω·∇b_j) */
        for (int j = 0; j < 4; j++)
          A[j][j] += p->sig_t[c] * M[j];
        for (int j = 1; j < 4; j++) {
          double g = om[j - 1] / (s * m->fr.u[j - 1]); /* (ω·∇b_j), постоянная */
          /* −∫ b_i (ω·∇b_j) dV. Градиент постоянен, поэтому интеграл есть
           * g_j·∫b_i dV, а ∫b_i dV равен НУЛЮ при i > 0 (ячейка симметрична
           * относительно своего центра). Значит член связывает тестовую функцию
           * ТОЛЬКО с постоянной модой, и лишний вклад при i == j был ошибкой. */
          A[j][0] -= g * M[0];
        }
        /* источник: σ_s/(4π)·φ + ε */
        for (int j = 0; j < 4; j++) {
          double q = p->sig_s[c] / (4.0 * M_PI) * phi[c * 4 + j];
          if (p->eps != NULL) q += p->eps[c * 4 + j];
          rhs[j] += q * M[j];
        }

        for (int32_t k = m->fstart[c]; k < m->fstart[c + 1]; k++) {
          int32_t f = m->flist[k];
          int mine_is_a = m->f[f].ca == c;
          /* ВНЕШНЯЯ НОРМАЛЬ ЯЧЕЙКИ. У внутренней грани её задаёт сторона: ca лежит
           * со стороны меньших координат, значит у неё нормаль +axis. У ГРАНИЧНОЙ
           * ячейка всегда записана в ca, и знак берётся из номера стенки —
           * иначе минус-стенка получила бы нормаль плюсовой. Ровно это и поймала
           * печь: решение уехало на 17%. */
          double sgn;
          if (m->f[f].cb < 0) {
            int wall = (int)(~m->f[f].cb);
            sgn = (wall & 1) ? 1.0 : -1.0;
          } else {
            sgn = mine_is_a ? 1.0 : -1.0;
          }
          double on = om[m->f[f].axis] * sgn;
          if (!(fabs(on) > 0.0)) continue;
          if (on > 0.0) { /* ВЫТОК: своя же функция, в СВОЁМ базисе */
            const double (*fs)[4] = mine_is_a ? fm[f] : fmb[f];
            for (int i = 0; i < 4; i++)
              for (int j = 0; j < 4; j++)
                A[j][i] += on * fs[i][j];
          } else { /* ВТОК: значение верхней по потоку ячейки */
            if (m->f[f].cb < 0) {
              /* граница: предписанный влёт, линейный по положению */
              double v[4][3];
              tr3_face_corners(m, f, v);
              double bl[4][4];
              face_mass2(m, (const double (*)[4][3]) & v, c, -1, bl);
              /* ∫ L_b b_j: L_b линейна, поэтому берём её значения в вершинах и
               * пользуемся тем же тождеством через смешанную матрицу с b_0 */
              for (int j = 0; j < 4; j++) {
                double acc = 0.0;
                for (int e = 1; e + 1 < 4; e++) {
                  double tri[3][3];
                  for (int q = 0; q < 3; q++)
                    for (int a = 0; a < 3; a++)
                      tri[q][a] = v[q == 0 ? 0 : (q == 1 ? e : e + 1)][a];
                  double area = tr3_poly_face_area(tri, 3);
                  double fv[3], bv[3][4];
                  for (int q = 0; q < 3; q++) {
                    if (p->wall_rho != NULL) {
                      fv[q] = bout[f]; /* стенка светит СВОИМ исходящим радиансом */
                    } else {
                      fv[q] = p->binc0;
                      for (int a = 0; a < 3; a++)
                        fv[q] += p->binc[a] * (tri[q][a] - p->binx0[a]);
                    }
                    cell_basis(m, c, tri[q], bv[q]);
                  }
                  double s1 = 0.0, sf = 0.0, sg = 0.0;
                  for (int q = 0; q < 3; q++) {
                    s1 += fv[q] * bv[q][j];
                    sf += fv[q];
                    sg += bv[q][j];
                  }
                  acc += area * (s1 + sf * sg) / 12.0;
                }
                rhs[j] -= on * acc;
                /* ВТЕКШАЯ МОЩНОСТЬ — из того же интеграла, что и правая часть,
                 * а не из площади на константу: при линейном влёте второе
                 * неверно, и баланс К13 померил бы не то. */
                if (j == 0) st->pin += -on * d->w[mm] * acc;
              }
            } else {
              int32_t up = mine_is_a ? m->f[f].cb : m->f[f].ca;
              /* fmx[f] есть ∫ b^ca_i b^cb_j; нужна ∫ b^up_i b^my_j */
              for (int j = 0; j < 4; j++) {
                double acc = 0.0;
                for (int i = 0; i < 4; i++)
                  acc += L[up * 4 + i] * (mine_is_a ? fmx[f][j][i] : fmx[f][i][j]);
                rhs[j] -= on * acc;
              }
            }
          }
        }

        /* СТРОКА БАЛАНСА СОХРАНЯЕТСЯ ДО РЕШЕНИЯ: ограничителю она нужна целой, а
         * solve4 матрицу разрушает. Без неё пришлось бы либо считать систему
         * дважды, либо срезать наклоны БЕЗ пересчёта среднего — а это ровно
         * режим SLOPE, у которого в T2 измерена утечка 4.3e-3 против 5e-16. */
        double a0row[4], rhs0 = rhs[0];
        for (int j = 0; j < 4; j++)
          a0row[j] = A[0][j];
        double cf[4];
        if (solve4(A, rhs, cf) != 0) {
          for (int j = 0; j < 4; j++)
            cf[j] = 0.0;
        }
        /* ОГРАНИЧИТЕЛЬ ПОЛОЖИТЕЛЬНОСТИ (К14/К6): наклоны срезаются множителем α,
         * а среднее ПЕРЕСЧИТЫВАЕТСЯ из строки баланса — иначе энергия течёт
         * (в T2 измерено 4.3e-3 против 5e-16). α берётся замкнутой формой. */
        if (p->limiter && corner_min(cf) < 0.0 && fabs(a0row[0]) > 0.0) {
          /* ЗАМКНУТАЯ ФОРМА (перенесена из T2): среднее есть АФФИННАЯ функция
           * коэффициента срезки, поэтому положительность и баланс достигаются
           * ОБА и за один шаг, а порога при этом не появляется — α берётся из
           * уравнения. */
          double kk = 0.0;
          for (int j = 1; j < 4; j++)
            kk += a0row[j] * cf[j];
          double dev = corner_min((double[4]){0.0, cf[1], cf[2], cf[3]});
          double c0h = rhs0 / a0row[0], kh = kk / a0row[0];
          double alpha = 1.0, den = kh - dev;
          if (den > 0.0) {
            alpha = c0h / den;
            if (alpha > 1.0) alpha = 1.0;
            if (alpha < 0.0) alpha = 0.0;
          }
          if (alpha < 1.0) {
            for (int j = 1; j < 4; j++)
              cf[j] *= alpha;
            cf[0] = (rhs0 - alpha * kk) / a0row[0];
            st->nclip++;
          }
        }
        for (int j = 0; j < 4; j++)
          L[c * 4 + j] = cf[j];
      }

      /* вклад в скалярный поток и в баланс */
      for (int32_t c = 0; c < nc; c++)
        for (int j = 0; j < 4; j++)
          phin[c * 4 + j] += d->w[mm] * L[c * 4 + j];
      for (int32_t f = 0; f < m->nf; f++) {
        if (m->f[f].cb >= 0) continue;
        /* внешняя нормаль граничной грани: +axis у «плюс»-стенки, −axis у «минус» */
        int wall = (int)(~m->f[f].cb);
        double on = om[m->f[f].axis] * ((wall & 1) ? 1.0 : -1.0);
        if (on > 0.0) {
          double v[4][3];
          tr3_face_corners(m, f, v);
          double bl[4][4];
          face_mass2(m, (const double (*)[4][3]) & v, m->f[f].ca, -1, bl);
          double acc = 0.0;
          for (int i = 0; i < 4; i++)
            acc += L[m->f[f].ca * 4 + i] * bl[i][0];
          st->pout += on * d->w[mm] * acc;
          /* ОБЛУЧЁННОСТЬ СТЕНКИ: то, что из ячейки вытекло, для стенки есть
           * ВХОДЯЩЕЕ. E = ∫L|ω·n|dω, и делится потом на площадь. */
          binf[f] += on * d->w[mm] * acc;
        }
      }
    }

    /* СТЕНКИ: новый исходящий радианс из накопленной облучённости. Это и есть
     * внешняя итерация «развёртка ↔ отражение», и она идёт в том же цикле, что
     * и итерация по рассеянию — они сходятся вместе. */
    if (p->wall_rho != NULL) {
      for (int32_t f = 0; f < m->nf; f++) {
        if (m->f[f].cb >= 0) continue;
        int wall = (int)(~m->f[f].cb);
        double e = binf[f] / farea[f];
        bout[f] =
            (p->wall_emit != NULL ? p->wall_emit[wall] : 0.0) + p->wall_rho[wall] * e / hsum[wall];
      }
    }
    resid = 0.0;
    for (int32_t i = 0; i < nc * 4; i++) {
      double dd = fabs(phin[i] - phi[i]);
      if (dd > resid) resid = dd;
      phi[i] = phin[i];
    }
    if (resid < tol) break;
  }
  /* поглощено = ∫(σ_t − σ_s)·φ */
  st->pabs = 0.0;
  for (int32_t c = 0; c < nc; c++) {
    double s = (double)m->csize[c];
    double vol = s * s * s * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
    st->pabs += (p->sig_t[c] - p->sig_s[c]) * phi[c * 4] * vol;
  }
  st->balance = st->pin - st->pout - st->pabs;
  st->iters = it;
  st->resid = resid;
  st->bout = bout; /* владение переходит вызывающему: сбор по пикселю читает это */

  free(L);
  free(phin);
  free(indeg);
  free(order);
  free(queue);
  free(fm);
  free(fmb);
  free(fmx);
  free(farea);
  free(binf);
  return 0;
}
