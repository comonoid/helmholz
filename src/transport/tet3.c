/* PLAN_TRANSPORT.md T4. Тетраэдризация разрезанной ячейки и точные моменты
 * линейного базиса — интегратор ПЕРЕНОСА над общей геометрией. */

#include "transport/tet3.h"
#include <math.h>
#include <string.h>

double tr3_tet_volume(const tr3_tet *t) {
  double a[3], b[3], c[3];
  for (int k = 0; k < 3; k++) {
    a[k] = t->v[1][k] - t->v[0][k];
    b[k] = t->v[2][k] - t->v[0][k];
    c[k] = t->v[3][k] - t->v[0][k];
  }
  double cr[3] = {b[1] * c[2] - b[2] * c[1], b[2] * c[0] - b[0] * c[2], b[0] * c[1] - b[1] * c[0]};
  return (a[0] * cr[0] + a[1] * cr[1] + a[2] * cr[2]) / 6.0;
}

int tr3_tet_from_poly(const hz_poly3 *p, tr3_tet *out, int max) {
  if (p->nv < 4 || p->nf < 4) return 0;
  int32_t apex = p->fl[0]; /* вершина веера: первая в петле первой грани */
  int n = 0;
  for (int32_t f = 0; f < p->nf; f++) {
    int32_t b0 = p->floff[f], k = p->floff[f + 1] - b0;
    if (k < 3) continue;
    int has = 0;
    for (int32_t e = 0; e < k; e++)
      if (p->fl[b0 + e] == apex) has = 1;
    if (has) continue; /* грань через вершину веера даёт нулевой объём */
    for (int32_t e = 1; e + 1 < k; e++) {
      if (n >= max) return -1;
      for (int c = 0; c < 3; c++) {
        out[n].v[0][c] = p->v[apex][c];
        out[n].v[1][c] = p->v[p->fl[b0]][c];
        out[n].v[2][c] = p->v[p->fl[b0 + e]][c];
        out[n].v[3][c] = p->v[p->fl[b0 + e + 1]][c];
      }
      n++;
    }
  }
  return n;
}

/* Базис {1, ξ, η, ζ} в точке x (единицы) при центре c и масштабе h. */
static void basis(const double x[3], const double c[3], double h, double b[4]) {
  b[0] = 1.0;
  for (int k = 0; k < 3; k++)
    b[k + 1] = (x[k] - c[k]) / h;
}

int tr3_mass_matrix(const hz_poly3 *p, const hz_frame *fr, const double c[3], double h,
                    double m[4][4]) {
  tr3_tet tet[TR3_MAXTET];
  int nt = tr3_tet_from_poly(p, tet, TR3_MAXTET);
  if (nt < 0) return 1;
  memset(m, 0, 16 * sizeof(double));
  double jac = fr->u[0] * fr->u[1] * fr->u[2];
  for (int t = 0; t < nt; t++) {
    double vol = tr3_tet_volume(&tet[t]);
    /* ЗНАК СОХРАНЯЕТСЯ: петли граней ориентированы наружу, поэтому у правильного
     * разбиения все тетраэдры одного знака, и сумма со знаком равна объёму. */
    double b[4][4];
    for (int i = 0; i < 4; i++)
      basis(tet[t].v[i], c, h, b[i]);
    /* ТОЧНО: ∫_T f g dV = V/20 · (Σ f_i g_i + (Σf_i)(Σg_i)) для линейных f, g.
     * Это тождество, а не квадратура: следствие ∫λ^α = 6V α!/(|α|+3)!. */
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        double s1 = 0.0, sf = 0.0, sg = 0.0;
        for (int q = 0; q < 4; q++) {
          s1 += b[q][i] * b[q][j];
          sf += b[q][i];
          sg += b[q][j];
        }
        m[i][j] += vol * (s1 + sf * sg) / 20.0;
      }
  }
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      m[i][j] *= jac;
  return 0;
}

double tr3_int_linear(const hz_poly3 *p, const hz_frame *fr, const double c[3], double h,
                      const double a[4]) {
  tr3_tet tet[TR3_MAXTET];
  int nt = tr3_tet_from_poly(p, tet, TR3_MAXTET);
  if (nt < 0) return 0.0;
  double acc = 0.0;
  for (int t = 0; t < nt; t++) {
    double vol = tr3_tet_volume(&tet[t]);
    /* ∫_T f dV = V·(f1+f2+f3+f4)/4 для линейной f — тоже тождество */
    double s = 0.0;
    for (int i = 0; i < 4; i++) {
      double b[4];
      basis(tet[t].v[i], c, h, b);
      double f = a[0];
      for (int k = 1; k < 4; k++)
        f += a[k] * b[k];
      s += f;
    }
    acc += vol * s / 4.0;
  }
  return acc * fr->u[0] * fr->u[1] * fr->u[2];
}
