/* pclip.c — отсечение треугольника коробкой ячейки. Разбор — в `pclip.h`. */
#include "pclip.h"
#include <assert.h>
#include <math.h>
#include <string.h>

/* Сравнение с нулём БЕЗ `==`: гейт держит `-Wfloat-equal`, а проверка «делитель
 * ноль» здесь не допуск, а именно точный ноль. */
static int pc_iszero(double x) {
  return !(x < 0.0) && !(x > 0.0);
}

/* Рабочая вершина. Помимо координат и происхождения несёт ЗНАК НЕСУЩЕЙ ПРЯМОЙ
 * звена, идущего ОТ неё к следующей: либо ребро треугольника, либо осевая
 * плоскость. Именно этот знак и делает вершину независимой от ячейки — по нему
 * выбирается формула, а не по тому, что лежит в соседней ячейке буфера. */
typedef struct {
  double v[3];
  hz_pclip_prov p;
  unsigned char lk; /* 0 — ребро треугольника, 1 — осевая плоскость */
  unsigned char le; /* lk = 0: номер ребра */
  unsigned char la; /* lk = 1: ось */
  double lc;        /* lk = 1: координата */
} pc_vtx;

/* Точка на отрезке PQ с координатой `c` по оси `a`. Координата по этой оси
 * ставится ТОЧНО, а не получается счётом: на ней держится и полуоткрытое
 * соглашение, и совпадение вершин соседей. */
static void pc_lerp(double *x, const double *P, const double *Q, int a, double c) {
  double t = (c - P[a]) / (Q[a] - P[a]);
  for (int k = 0; k < 3; k++)
    x[k] = P[k] + t * (Q[k] - P[k]);
  x[a] = c;
}

/* Новая вершина на пересечении звена (P → Q) с плоскостью (a, c).
 * Схема выбирается по несущей прямой звена, а не по его концам. */
static void pc_isect(pc_vtx *X, const pc_vtx *P, const pc_vtx *Q, const double tri[3][3],
                     const double *nrm, double dpl, int a, double c, int naive, int *nfall) {
  int canon = 0;
  memset(&X->p, 0, sizeof X->p);
  if (P->lk == 0) {
    int e = P->le;
    const double *Pa = tri[e], *Pb = tri[(e + 1) % 3];
    X->p.kind = 1;
    X->p.e = (unsigned char)e;
    X->p.a1 = (unsigned char)a;
    X->p.c1 = c;
    if (!pc_iszero(Pb[a] - Pa[a])) {
      canon = 1;
      if (!naive) pc_lerp(X->v, Pa, Pb, a, c);
    }
  } else {
    int aq = P->la;
    double cq = P->lc;
    /* `aq == a` означало бы звено, лежащее в ТОЙ ЖЕ плоскости, по которой идёт
     * отсечение: у него оба конца по эту сторону, и сюда управление не
     * приходит. Ключ помечается родом 3, чтобы вырожденная вершина не выдала
     * себя за законную при сверке соседей. */
    if (aq != a) {
      int a1 = (a < aq) ? a : aq, a2 = (a < aq) ? aq : a;
      double c1 = (a1 == a) ? c : cq, c2 = (a2 == a) ? c : cq;
      int a3 = 3 - a1 - a2;
      X->p.kind = 2;
      X->p.a1 = (unsigned char)a1;
      X->p.c1 = c1;
      X->p.a2 = (unsigned char)a2;
      X->p.c2 = c2;
      if (!pc_iszero(nrm[a3])) {
        canon = 1;
        if (!naive) {
          X->v[a1] = c1;
          X->v[a2] = c2;
          X->v[a3] = (dpl - nrm[a1] * c1 - nrm[a2] * c2) / nrm[a3];
        }
      }
    } else {
      X->p.kind = 3;
    }
  }
  /* ОТСТУПЛЕНИЕ. Каноническая формула невозможна только в вырожденной
   * обстановке (ребро параллельно секущей плоскости; плоскость треугольника
   * содержит направление третьей оси). Отступление СЧИТАЕТСЯ и печатается:
   * молчаливое отступление вернуло бы зависимость вершины от ячейки, то есть
   * ровно Г50, и ноль расхождений стал бы неправдой. */
  if (!canon) {
    if (!naive) (*nfall)++;
    pc_lerp(X->v, P->v, Q->v, a, c);
  } else if (naive) {
    pc_lerp(X->v, P->v, Q->v, a, c);
  }
}

/* Выдать вершину. Совпавшую ПОБИТОВО с предыдущей не удваивать: вырожденное
 * звено нулевой длины не несёт ничего, а счёт вершин ограничен постусловием. */
static int pc_emit(pc_vtx *out, int n, const pc_vtx *x) {
  if (n > 0 && memcmp(out[n - 1].v, x->v, sizeof x->v) == 0) {
    out[n - 1].lk = x->lk;
    out[n - 1].le = x->le;
    out[n - 1].la = x->la;
    out[n - 1].lc = x->lc;
    return n;
  }
  /* ПОСТУСЛОВИЕ §241.3, ОНО ЖЕ ГРАНИЦА БУФЕРА. Обрезания здесь нет
   * сознательно (§244 п. 1): обрезание спрятало бы ошибку отсечения. */
  assert(n < HZ_PCLIP_MAXV);
  out[n] = *x;
  return n + 1;
}

/* Один проход Сазерленда—Ходжмана по одной осевой полуплоскости.
 * `islo`: внутри значит `>= c`; иначе внутри значит `< c` (соглашение А475). */
static int pc_stage(const pc_vtx *in, int n, pc_vtx *out, const double tri[3][3], const double *nrm,
                    double dpl, int a, double c, int islo, int naive, int *nfall) {
  int m = 0;
  for (int i = 0; i < n; i++) {
    const pc_vtx *P = &in[i], *Q = &in[(i + 1) % n];
    int ip = islo ? (P->v[a] >= c) : (P->v[a] < c);
    int iq = islo ? (Q->v[a] >= c) : (Q->v[a] < c);
    if (ip) m = pc_emit(out, m, P);
    if (ip != iq) {
      pc_vtx X;
      pc_isect(&X, P, Q, tri, nrm, dpl, a, c, naive, nfall);
      if (ip) {
        /* уходим наружу: следующее звено идёт ПО секущей плоскости */
        X.lk = 1;
        X.le = 0;
        X.la = (unsigned char)a;
        X.lc = c;
      } else {
        /* возвращаемся внутрь: звено продолжает прежнюю несущую прямую */
        X.lk = P->lk;
        X.le = P->le;
        X.la = P->la;
        X.lc = P->lc;
      }
      m = pc_emit(out, m, &X);
    }
  }
  return m;
}

int hz_pclip_tri_ex(const double *A, const double *B, const double *C, const double *lo,
                    const double *hi, hz_pclip_poly *out, int naive) {
  double tri[3][3];
  for (int k = 0; k < 3; k++) {
    tri[0][k] = A[k];
    tri[1][k] = B[k];
    tri[2][k] = C[k];
  }
  /* Плоскость треугольника: нормаль НЕ нормируется — корень ничего не добавляет
   * к точности, а лишняя операция есть лишний источник расхождения. Формула
   * одна на все ячейки, потому и вершина рода 2 одна на все ячейки. */
  double u[3], w[3], nrm[3];
  for (int k = 0; k < 3; k++) {
    u[k] = tri[1][k] - tri[0][k];
    w[k] = tri[2][k] - tri[0][k];
  }
  nrm[0] = u[1] * w[2] - u[2] * w[1];
  nrm[1] = u[2] * w[0] - u[0] * w[2];
  nrm[2] = u[0] * w[1] - u[1] * w[0];
  double dpl = nrm[0] * tri[0][0] + nrm[1] * tri[0][1] + nrm[2] * tri[0][2];

  pc_vtx buf[2][HZ_PCLIP_MAXV];
  int nfall = 0, n = 3, cur = 0;
  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < 3; k++)
      buf[0][i].v[k] = tri[i][k];
    memset(&buf[0][i].p, 0, sizeof buf[0][i].p);
    buf[0][i].p.kind = 0;
    buf[0][i].p.e = (unsigned char)i;
    buf[0][i].lk = 0;
    buf[0][i].le = (unsigned char)i;
    buf[0][i].la = 0;
    buf[0][i].lc = 0.0;
  }
  out->nv = 0;
  out->nfall = 0;
  for (int a = 0; a < 3; a++) {
    for (int s = 0; s < 2; s++) {
      n = pc_stage(buf[cur], n, buf[1 - cur], tri, nrm, dpl, a, s ? hi[a] : lo[a], s == 0, naive,
                   &nfall);
      cur = 1 - cur;
      if (n < 3) {
        out->nfall = nfall;
        return 0;
      }
    }
  }
  /* Замыкание: последняя вершина, совпавшая с первой, — то же вырожденное
   * звено, что снимает `pc_emit` внутри прохода. */
  while (n >= 2 && memcmp(buf[cur][n - 1].v, buf[cur][0].v, sizeof buf[cur][0].v) == 0)
    n--;
  out->nfall = nfall;
  if (n < 3) return 0;
  assert(n <= HZ_PCLIP_MAXV);
  for (int i = 0; i < n; i++) {
    for (int k = 0; k < 3; k++)
      out->v[i][k] = buf[cur][i].v[k];
    out->p[i] = buf[cur][i].p;
  }
  out->nv = n;
  return n;
}

int hz_pclip_tri(const double *A, const double *B, const double *C, const double *lo,
                 const double *hi, hz_pclip_poly *out) {
  return hz_pclip_tri_ex(A, B, C, lo, hi, out, 0);
}

double hz_pclip_area(const hz_pclip_poly *P) {
  if (P->nv < 3) return 0.0;
  double s[3] = {0.0, 0.0, 0.0};
  for (int i = 1; i + 1 < P->nv; i++) {
    double u[3], w[3];
    for (int k = 0; k < 3; k++) {
      u[k] = P->v[i][k] - P->v[0][k];
      w[k] = P->v[i + 1][k] - P->v[0][k];
    }
    s[0] += u[1] * w[2] - u[2] * w[1];
    s[1] += u[2] * w[0] - u[0] * w[2];
    s[2] += u[0] * w[1] - u[1] * w[0];
  }
  return 0.5 * sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
}
