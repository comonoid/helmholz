/* PLAN_TRANSPORT.md, этап B. BiCGStab по всему состоянию `(φ, bout, sout)`.
 * Разбор и оговорки — в `krylov3.h`, здесь только механика. */

#include "transport/krylov3.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
  tr3_problem lin; /* оператор: источники ЗАНУЛЕНЫ, ограничитель выключен */
  int32_t nphi, nb, ns, n;
  double *w;  /* рабочее состояние прохода */
  long npass; /* ПРОХОДОВ развёртки — единица сравнения с рядом Неймана */
  int failed; /* проход вернул ошибку: дальше считать нельзя */
} kop;

static double kdot(const double *a, const double *b, int32_t n) {
  double s = 0.0;
  for (int32_t i = 0; i < n; i++)
    s += a[i] * b[i];
  return s;
}

static double knrm(const double *a, int32_t n) {
  double m = 0.0;
  for (int32_t i = 0; i < n; i++)
    if (fabs(a[i]) > m) m = fabs(a[i]);
  return m;
}

/* out = A·v = v − T·v.
 *
 * `T·v` есть ОДИН проход развёртки от состояния `v` при выключенных источниках.
 * Состояние передаётся целиком: `φ` через `warm_start`, поверхности через
 * `bout_in`/`sout_in` — без последних половина состояния терялась бы между
 * вызовами, и оператор перестал бы быть тем, чем притворяется (К84). */
static void kapply(kop *k, const double *v, double *out) {
  if (k->failed) return;
  tr3_stats st;
  memcpy(k->w, v, (size_t)k->nphi * sizeof(double));
  k->lin.bout_in = v + k->nphi;
  k->lin.sout_in = k->ns > 0 ? v + k->nphi + k->nb : NULL;
  if (tr3_sweep_solve(&k->lin, 1, 0.0, k->w, &st) != 0) {
    k->failed = 1;
    return;
  }
  k->npass++;
  for (int32_t i = 0; i < k->nphi; i++)
    out[i] = v[i] - k->w[i];
  for (int32_t i = 0; i < k->nb; i++)
    out[k->nphi + i] = v[k->nphi + i] - st.bout[i];
  for (int32_t i = 0; i < k->ns; i++)
    out[k->nphi + k->nb + i] = v[k->nphi + k->nb + i] - st.sout[i];
  free(st.bout);
  free(st.sout);
}

int tr3_krylov_solve(const tr3_problem *p, int maxit, double tol, double *phi, double *bout,
                     double *sout, tr3_kstats *st) {
  memset(st, 0, sizeof *st);
  st->why = "не запускалось";
  int32_t ncell = p->m->ncell, nf = p->m->nf;
  int32_t nse = p->cut != NULL ? p->cut->nse : 0;
  kop k;
  memset(&k, 0, sizeof k);
  k.nphi = ncell * 4;
  k.nb = nf * 4;
  k.ns = nse * 4;
  k.n = k.nphi + k.nb + k.ns;

  /* ОПЕРАТОР: тот же `p`, но БЕЗ источников и БЕЗ ограничителя (см. заголовок).
   * Массивы нулей заводятся здесь же, а не подставляются `NULL`: `NULL` у
   * `wall_emit` значит «нет излучения», но у `wall_rho` — «нет отражения», и
   * путать эти два случая нельзя. */
  double wz[6] = {0, 0, 0, 0, 0, 0};
  double *fz = calloc((size_t)(p->nfacet > 0 ? p->nfacet : 1), sizeof(double));
  double *b = calloc((size_t)k.n, sizeof(double));
  double *x = calloc((size_t)k.n, sizeof(double));
  double *rr = calloc((size_t)k.n, sizeof(double));
  double *rh = calloc((size_t)k.n, sizeof(double));
  double *pv = calloc((size_t)k.n, sizeof(double));
  double *vv = calloc((size_t)k.n, sizeof(double));
  double *ss = calloc((size_t)k.n, sizeof(double));
  double *tt = calloc((size_t)k.n, sizeof(double));
  k.w = calloc((size_t)(k.nphi > 0 ? k.nphi : 1), sizeof(double));
  /* НУЛЕВОЕ поверхностное состояние для правой части: `S` применяется к нулю
   * ТЕМ ЖЕ путём, что и оператор, поэтому нули подаются ЯВНО, а не через
   * холодный старт, который сначала кладёт `bout = emit`. */
  double *zb = calloc((size_t)(k.nb > k.ns ? k.nb : k.ns) + 1, sizeof(double));
  if (fz == NULL || b == NULL || x == NULL || rr == NULL || rh == NULL || pv == NULL ||
      vv == NULL || ss == NULL || tt == NULL || k.w == NULL || zb == NULL) {
    free(fz);
    free(b);
    free(x);
    free(rr);
    free(rh);
    free(pv);
    free(vv);
    free(ss);
    free(tt);
    free(k.w);
    free(zb);
    return 1;
  }
  k.lin = *p;
  k.lin.limiter = 0;
  k.lin.warm_start = 1;
  k.lin.trace = 0;
  k.lin.wall_emit = p->wall_rho != NULL ? wz : NULL;
  k.lin.facet_emit = fz;
  k.lin.eps = NULL;
  k.lin.eps_dir[0] = k.lin.eps_dir[1] = k.lin.eps_dir[2] = 0.0;
  k.lin.binc0 = 0.0;
  k.lin.binc[0] = k.lin.binc[1] = k.lin.binc[2] = 0.0;

  /* ПРАВАЯ ЧАСТЬ `b = S(0)` — И СЧИТАТЬ ЕЁ НАДО ТЕМ ЖЕ ПУТЁМ, ЧТО И ОПЕРАТОР.
   *
   * Соблазн: взять `warm_start = 0`, то есть обычный холодный старт. НЕВЕРНО:
   * холодный старт СНАЧАЛА кладёт `bout = emit` и лишь потом идёт развёрткой,
   * то есть считает `S` применённым к состоянию `(0, emit)`, а не к нулевому.
   * Это лишние полшага, и из-за них система решалась с чужой правой частью:
   * Крылов сходился к состоянию с `φ ≈ 0` и невязкой `7e-10` — то есть
   * безупречно решал НЕ ТУ задачу.
   *
   * Правильно: тот же путь, что у оператора (`warm_start = 1`, поверхности
   * поданы явно), только источники ЖИВЫЕ, а состояние на входе НУЛЕВОЕ. */
  {
    tr3_problem p0 = *p;
    p0.limiter = 0;
    p0.warm_start = 1;
    p0.trace = 0;
    p0.bout_in = zb;
    p0.sout_in = k.ns > 0 ? zb : NULL;
    tr3_stats s0;
    if (tr3_sweep_solve(&p0, 1, 0.0, b, &s0) != 0) {
      free(fz);
      free(b);
      free(x);
      free(rr);
      free(rh);
      free(pv);
      free(vv);
      free(ss);
      free(tt);
      free(k.w);
      free(zb);
      return 1;
    }
    k.npass++;
    memcpy(b + k.nphi, s0.bout, (size_t)k.nb * sizeof(double));
    if (k.ns > 0) memcpy(b + k.nphi + k.nb, s0.sout, (size_t)k.ns * sizeof(double));
    free(s0.bout);
    free(s0.sout);
  }

  /* СТАРТ С ЗАДАННОГО ПРИБЛИЖЕНИЯ, а не обязательно с нуля: тёплый старт по
   * К85 даёт немного, но он бесплатен и складывается с остальным. */
  memcpy(x, phi, (size_t)k.nphi * sizeof(double));
  memcpy(x + k.nphi, bout, (size_t)k.nb * sizeof(double));
  if (k.ns > 0) memcpy(x + k.nphi + k.nb, sout, (size_t)k.ns * sizeof(double));

  kapply(&k, x, rr);
  for (int32_t i = 0; i < k.n; i++)
    rr[i] = b[i] - rr[i];
  memcpy(rh, rr, (size_t)k.n * sizeof(double));
  double bn = knrm(b, k.n);
  if (!(bn > 0.0)) bn = 1.0;
  double rho = 1.0, alpha = 1.0, omega = 1.0;
  st->why = "предел итераций";
  int it = 0;
  for (; it < maxit; it++) {
    double rn = knrm(rr, k.n);
    if (rn <= tol * bn) {
      st->why = "сошлось";
      st->converged = 1;
      break;
    }
    double rho1 = kdot(rh, rr, k.n);
    if (!(fabs(rho1) > 0.0)) {
      st->why = "обрыв: rho = 0";
      break;
    }
    if (it == 0) {
      memcpy(pv, rr, (size_t)k.n * sizeof(double));
    } else {
      double beta = (rho1 / rho) * (alpha / omega);
      for (int32_t i = 0; i < k.n; i++)
        pv[i] = rr[i] + beta * (pv[i] - omega * vv[i]);
    }
    kapply(&k, pv, vv);
    if (k.failed) {
      st->why = "обрыв: проход вернул ошибку";
      break;
    }
    double den = kdot(rh, vv, k.n);
    if (!(fabs(den) > 0.0)) {
      st->why = "обрыв: (rh, v) = 0";
      break;
    }
    alpha = rho1 / den;
    for (int32_t i = 0; i < k.n; i++)
      ss[i] = rr[i] - alpha * vv[i];
    /* ПОЛОВИННЫЙ ШАГ: если он уже дал невязку, настоящая невязка есть `s`, а не
     * старое `r`. Печать старого `r` однажды выглядела как обрыв метода, хотя
     * метод работал (доклад 07-28). */
    if (knrm(ss, k.n) <= tol * bn) {
      for (int32_t i = 0; i < k.n; i++)
        x[i] += alpha * pv[i];
      memcpy(rr, ss, (size_t)k.n * sizeof(double));
      st->why = "сошлось на половинном шаге";
      st->converged = 1;
      it++;
      break;
    }
    kapply(&k, ss, tt);
    if (k.failed) {
      st->why = "обрыв: проход вернул ошибку";
      break;
    }
    double tsd = kdot(tt, tt, k.n);
    if (!(tsd > 0.0)) {
      st->why = "обрыв: (t, t) = 0";
      break;
    }
    omega = kdot(tt, ss, k.n) / tsd;
    for (int32_t i = 0; i < k.n; i++) {
      x[i] += alpha * pv[i] + omega * ss[i];
      rr[i] = ss[i] - omega * tt[i];
    }
    if (!(fabs(omega) > 0.0)) {
      st->why = "обрыв: omega = 0";
      break;
    }
    rho = rho1;
  }
  st->iters = it;
  st->npass = k.npass;
  st->resid = knrm(rr, k.n) / bn;

  memcpy(phi, x, (size_t)k.nphi * sizeof(double));
  memcpy(bout, x + k.nphi, (size_t)k.nb * sizeof(double));
  if (k.ns > 0) memcpy(sout, x + k.nphi + k.nb, (size_t)k.ns * sizeof(double));

  free(fz);
  free(b);
  free(x);
  free(rr);
  free(rh);
  free(pv);
  free(vv);
  free(ss);
  free(tt);
  free(k.w);
  free(zb);
  return 0;
}

/* ПРОВЕРКА ПРИБОРА (см. заголовок): невязка ЗАДАННОГО состояния под тем же
 * оператором. Собирается той же машинерией, что и решатель, — иначе она
 * проверяла бы другой оператор, и смысл терялся бы целиком. */
int tr3_krylov_residual(const tr3_problem *p, const double *phi, const double *bout,
                        const double *sout, double *rel) {
  int32_t nphi = p->m->ncell * 4, nb = p->m->nf * 4;
  int32_t ns = (p->cut != NULL ? p->cut->nse : 0) * 4;
  int32_t n = nphi + nb + ns;
  double *x = calloc((size_t)n, sizeof(double));
  double *ax = calloc((size_t)n, sizeof(double));
  double *b = calloc((size_t)n, sizeof(double));
  double *fz = calloc((size_t)(p->nfacet > 0 ? p->nfacet : 1), sizeof(double));
  double *w = calloc((size_t)(nphi > 0 ? nphi : 1), sizeof(double));
  if (x == NULL || ax == NULL || b == NULL || fz == NULL || w == NULL) {
    free(x);
    free(ax);
    free(b);
    free(fz);
    free(w);
    return 1;
  }
  memcpy(x, phi, (size_t)nphi * sizeof(double));
  memcpy(x + nphi, bout, (size_t)nb * sizeof(double));
  if (ns > 0) memcpy(x + nphi + nb, sout, (size_t)ns * sizeof(double));

  double wz[6] = {0, 0, 0, 0, 0, 0};
  kop k;
  memset(&k, 0, sizeof k);
  k.nphi = nphi;
  k.nb = nb;
  k.ns = ns;
  k.n = n;
  k.w = w;
  k.lin = *p;
  k.lin.limiter = 0;
  k.lin.warm_start = 1;
  k.lin.trace = 0;
  k.lin.wall_emit = p->wall_rho != NULL ? wz : NULL;
  k.lin.facet_emit = fz;
  k.lin.eps = NULL;
  k.lin.eps_dir[0] = k.lin.eps_dir[1] = k.lin.eps_dir[2] = 0.0;
  k.lin.binc0 = 0.0;
  k.lin.binc[0] = k.lin.binc[1] = k.lin.binc[2] = 0.0;

  /* правая часть считается ТЕМ ЖЕ путём, что в решателе (см. там разбор) */
  double *zb = calloc((size_t)(nb > ns ? nb : ns) + 1, sizeof(double));
  if (zb == NULL) {
    free(x);
    free(ax);
    free(b);
    free(fz);
    free(w);
    return 1;
  }
  tr3_problem p0 = *p;
  p0.limiter = 0;
  p0.warm_start = 1;
  p0.trace = 0;
  p0.bout_in = zb;
  p0.sout_in = ns > 0 ? zb : NULL;
  tr3_stats s0;
  if (tr3_sweep_solve(&p0, 1, 0.0, b, &s0) != 0) {
    free(x);
    free(ax);
    free(b);
    free(fz);
    free(w);
    return 1;
  }
  memcpy(b + nphi, s0.bout, (size_t)nb * sizeof(double));
  if (ns > 0) memcpy(b + nphi + nb, s0.sout, (size_t)ns * sizeof(double));
  free(s0.bout);
  free(s0.sout);

  kapply(&k, x, ax);
  double num = 0.0;
  for (int32_t i = 0; i < n; i++) {
    double d = fabs(ax[i] - b[i]);
    if (d > num) num = d;
  }
  double den = knrm(b, n);
  *rel = den > 0.0 ? num / den : num;
  free(zb);
  free(x);
  free(ax);
  free(b);
  free(fz);
  free(w);
  return k.failed ? 1 : 0;
}
