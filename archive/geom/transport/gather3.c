/* PLAN_TRANSPORT.md, проход 3 со ЗЕРКАЛАМИ (К5 и К3). Разбор — в gather3.h. */

#include "transport/gather3.h"
#include <math.h>
#include <string.h>

/* Значение DG1-поля, хранимого в базисе ячейки `c`, в мировой точке `x`.
 * ЧИТАТЬ НУЛЕВОЙ КОЭФФИЦИЕНТ КАК РАДИАНС НЕЛЬЗЯ (К39): базис центрирован на
 * ЯЧЕЙКЕ, а элемент лежит на её краю, поэтому у постоянного поля `c`
 * представитель минимальной нормы равен (0.8c, 0.4c, 0, 0). */
static double dg1_at(const tr3_mesh *m, int32_t c, const double coef[4], const double x[3]) {
  double s = (double)m->csize[c], v = coef[0];
  for (int a = 0; a < 3; a++) {
    double xu = (x[a] - m->fr.o[a]) / m->fr.u[a];
    v += coef[a + 1] * ((xu - ((double)m->clo[c][a] + 0.5 * s)) / s);
  }
  return v;
}

/* Поверхностный элемент, в который попал луч: БЛИЖАЙШИЙ ПО НОРМАЛИ.
 * Искать по номеру фасета нельзя — в примитивном пути (К15) номер не
 * проставляется вовсе, потому что попадание считает примитив, и фасет к ответу
 * отношения не имеет. */
static int32_t pick_selem(const tr3_cut *cut, int32_t mc, const double n[3]) {
  double best = -2.0;
  int32_t bi = -1;
  for (int32_t k = cut->sestart[mc]; k < cut->sestart[mc + 1]; k++) {
    int32_t e = cut->selist[k];
    double dp = cut->se[e].n[0] * n[0] + cut->se[e].n[1] * n[1] + cut->se[e].n[2] * n[2];
    if (dp > best) {
      best = dp;
      bi = e;
    }
  }
  return bi;
}

/* Граничная грань, накрывающая точку выхода луча из куба.
 *
 * К80. СЧЁТ ИДЁТ В ЕДИНИЦАХ ДЕРЕВА, А ЛУЧ ПРИХОДИТ В МИРЕ, И ПЕРЕВОД
 * ОБЯЗАТЕЛЕН. Первая редакция сравнивала МИРОВУЮ координату с `nwall` и брала
 * от неё `floor`, то есть молча предполагала, что единица дерева равна единице
 * мира. Пока кадр был `u = (1,1,1)`, это выполнялось, и ошибка не проявлялась
 * ВОВСЕ; при `u = 0.5` (комната того же размера на вдвое более мелкой сетке)
 * луч уводился к координате `nwall` в МИРЕ — то есть за пределы комнаты, — а
 * индекс грани брался от мировой координаты. Картинка при этом не ломалась
 * заметно, а МЕНЯЛАСЬ: радианс уходил с `0.45…3.38` на `0.0…20.9`.
 *
 * Это ровно тот класс, что `PLAN_CUT.md` Г21 и Г46: величина, посчитанная в
 * единичных координатах вместо мировых, не проявляется на кубическом единичном
 * кадре и вылезает только тогда, когда кадр перестаёт быть единичным. Здесь его
 * вскрыл свип по размеру ячейки (К68), который для того и заводился. */
static int32_t wall_face(const tr3_gather *g, const double o[3], const double d[3], double hp[3]) {
  if (g->wallidx == NULL) return -1;
  const hz_frame *fr = &g->m->fr;
  double ou[3], du[3];
  for (int a = 0; a < 3; a++) {
    ou[a] = (o[a] - fr->o[a]) / fr->u[a];
    du[a] = d[a] / fr->u[a];
  }
  double tex = 1e300;
  int wall = -1;
  for (int a = 0; a < 3; a++) {
    if (!(fabs(du[a]) > 0.0)) continue;
    double lim = du[a] > 0.0 ? (double)g->nwall : 0.0;
    double tt = (lim - ou[a]) / du[a];
    if (tt > 0.0 && tt < tex) {
      tex = tt;
      wall = 2 * a + (du[a] > 0.0 ? 1 : 0);
    }
  }
  if (wall < 0) return -1;
  double hu[3];
  for (int a = 0; a < 3; a++) {
    hu[a] = ou[a] + tex * du[a];
    /* наружу точка отдаётся В МИРЕ: `dg1_at` переводит её обратно сам, и
     * возвращать единицы значило бы применить кадр дважды (Г46 наизнанку) */
    hp[a] = fr->o[a] + fr->u[a] * hu[a];
  }
  int axis = wall / 2, u = (axis + 1) % 3, v = (axis + 2) % 3;
  if (u > v) {
    int t = u;
    u = v;
    v = t;
  }
  int32_t iu = (int32_t)floor(hu[u]), iv = (int32_t)floor(hu[v]);
  if (iu < 0) iu = 0;
  if (iu >= g->nwall) iu = g->nwall - 1;
  if (iv < 0) iv = 0;
  if (iv >= g->nwall) iv = g->nwall - 1;
  return g->wallidx[((int32_t)wall * g->nwall + iu) * g->nwall + iv];
}

/* К45: НАБЛЮДАТЕЛЬ СЕГМЕНТОВ СРЕДЫ. На каждый средовый отрезок [ta, tb] (те же,
 * по которым марш копит τ) считает точный интеграл
 *
 *     I = ∫_0^ℓ e^{−(τa + σ_t·u)} · (A + B·u) du,   A + B·u = φ(s),
 *
 * где φ — DG1 (по отрезку линейна), и прибавляет (σ_s/4π)·I к каналу с текущим
 * зеркальным пропусканием. Замкнутая форма, не квадратура: при σ_t = 0 она
 * вырождается в трапецию той же точности, что и в T1, — частный случай
 * проверенного. τa ведётся КОНТЕКСТОМ от глаза через все зеркальные отскоки:
 * марш между отскоками τ не переносит. */
typedef struct {
  const tr3_gather *g;
  const double *o, *d;
  double thr;  /* произведение зеркальных долей ПРОЙДЕННЫХ отскоков */
  double tau;  /* накопленная τ ОТ ГЛАЗА до начала текущего сегмента */
  double *val; /* [nch] канальный аккумулятор */
  int nch;
} seg45;

/* 4π — фазовая функция изотропного рассеяния, та же, что в шаге C развёртки */
#define SEG45_FOUR_PI 12.566370614359172

static void seg45_add(void *vctx, int32_t ni, double ta, double tb) {
  seg45 *sx = vctx;
  const tr3_gather *g = sx->g;
  if (g->phi == NULL || g->sig_t == NULL || g->m == NULL) return;
  const int32_t c = g->m->cellof[ni];
  if (c < 0) return;
  const double st = g->sig_t[c], ss = g->sig_s != NULL ? g->sig_s[c] : 0.0;
  const double ell = tb - ta;
  if (!(ell > 0.0)) return;
  double p0[3], p1[3], f0[4] = {0, 0, 0, 0}, f1[4] = {0, 0, 0, 0};
  for (int a = 0; a < 3; a++) {
    p0[a] = sx->o[a] + ta * sx->d[a];
    p1[a] = sx->o[a] + tb * sx->d[a];
  }
  f0[0] = dg1_at(g->m, c, g->phi + (size_t)c * 4, p0);
  f1[0] = dg1_at(g->m, c, g->phi + (size_t)c * 4, p1);
  const double A = f0[0], B = (f1[0] - f0[0]) / ell;
  /* имя `I` занято мнимой единицей complex.h из общего слоя — не трогать */
  double intg;
  if (st > 0.0) {
    const double q = exp(-st * ell);
    intg = A * (1.0 - q) / st + B * (1.0 - q * (1.0 + st * ell)) / (st * st);
  } else {
    intg = ell * (A + 0.5 * B * ell); /* вакуум внутри ячейки: трапеция точна */
  }
  const double add = sx->thr * (ss / SEG45_FOUR_PI) * exp(-sx->tau) * intg;
  for (int c2 = 0; c2 < sx->nch; c2++)
    sx->val[c2] += add;
  sx->tau += st * ell;
}

double tr3_gather_ray(const tr3_gather *g, const double o[3], const double d[3], int *nbounce,
                      double *out) {
  /* ОДИН МАРШ НА ВСЕ КАНАЛЫ (этап E): геометрия у них общая, различаются только
   * ЗНАЧЕНИЯ хранимого. Марш есть почти вся цена кадра, поэтому звать сбор по
   * разу на канал значило бы утроить кадр там, где честная цена — проценты. */
  const int nch = g->nch > 0 ? g->nch : 1;
  const int32_t nse4 = (g->cut != NULL && g->cut->nse > 0 ? g->cut->nse : 1) * 4;
  double p[3], dir[3], thr = 1.0;
  double val[TR3_MAXCH];
  for (int c = 0; c < nch; c++)
    val[c] = 0.0;
  memcpy(p, o, sizeof p);
  memcpy(dir, d, sizeof dir);
  int nb = 0;
  /* К45: наблюдатель сегментов ставится, только если среда задана; иначе —
   * прежний марш, мир побитово прежний (встроенный негативный контроль). */
  seg45 sx = {.g = g, .o = p, .d = dir, .thr = 1.0, .tau = 0.0, .val = val, .nch = nch};
  const int withmed = g->phi != NULL && g->sig_t != NULL;
  for (;;) {
    tr3_hit h;
    if (withmed)
      tr3_march_sink(g->sc, p, dir, -1.0, &h, seg45_add, &sx);
    else
      tr3_march(g->sc, p, dir, -1.0, &h);
    if (!h.hit) { /* ушёл в стенку куба */
      double hp[3];
      int32_t f = wall_face(g, p, dir, hp);
      if (f >= 0) {
        /* К45: стенка видна ЧЕРЕЗ среду — на T = exp(−τ) от глаза до неё */
        double T = withmed ? exp(-sx.tau) : 1.0;
        for (int c = 0; c < nch; c++)
          val[c] += thr * T *
                    dg1_at(g->m, g->m->f[f].ca,
                           g->bout + (size_t)c * (size_t)g->m->nf * 4 + (size_t)f * 4, hp);
      }
      break;
    }
    int32_t mc = g->m->cellof[h.cell];
    if (mc < 0) break;
    int32_t e = pick_selem(g->cut, mc, h.n);
    if (e < 0) break;
    double spec = 0.0;
    if (g->facet_spec != NULL && g->cut->se[e].facet < g->nfacet)
      spec = g->facet_spec[g->cut->se[e].facet];
    if (spec > 0.0 && nb < g->maxbounce) {
      /* ЗАКОН ОТРАЖЕНИЯ по нормали ПРИМИТИВА: d' = d − 2(d·n)n.
       * Точка старта берётся РОВНО на поверхности, без смещения по нормали:
       * самопопадание исключено строгостью неравенств в марше (ray3.h). */
      double dn = dir[0] * h.n[0] + dir[1] * h.n[1] + dir[2] * h.n[2];
      for (int a = 0; a < 3; a++) {
        p[a] = h.p[a];
        dir[a] -= 2.0 * dn * h.n[a];
      }
      thr *= spec;
      sx.thr = thr; /* К45: сегменты ПОСЛЕ отскока фильтруются им же */
      nb++;
      continue;
    }
    /* К45: поверхность тоже видна через среду — тот же множитель T */
    double Tsurf = withmed ? exp(-sx.tau) : 1.0;
    for (int c = 0; c < nch; c++)
      val[c] +=
          thr * Tsurf * dg1_at(g->m, mc, g->sout + (size_t)c * (size_t)nse4 + (size_t)e * 4, h.p);
    break;
  }
  if (nbounce != NULL) *nbounce = nb;
  for (int c = 0; c < nch; c++) {
    if (!(val[c] > 0.0)) val[c] = 0.0;
    if (out != NULL) out[c] = val[c];
  }
  return val[0];
}
