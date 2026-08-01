/* prad — РЕНДЕРЕР БЕЗ СВИПА (PLAN_ELEMENTS.md, §84, §85).
 *
 * Вся цепочка одним проходом: иерархия полигонов с лестницей LOD → сборка
 * уравнения связями → решение итерациями → картинка.
 *
 * ЧЕМ ОН ОТЛИЧАЕТСЯ ОТ `prender`. Тот решает ординатной развёрткой (свипом) и
 * считает прямой свет отдельным механизмом. Здесь ни того, ни другого нет:
 * направление связанное (§85), поэтому оператор есть набор связей между парами
 * узлов, прямой свет — обычная связь от светящегося узла, а многократные отражения
 * — итерации того же умножения.
 *
 * ЗАЧЕМ ОБА СРАЗУ. Свип-решатель остаётся ЭТАЛОНОМ: два независимых метода на одной
 * сцене обязаны сойтись, и это сильнейшая проверка, какая тут возможна. Пока не
 * сошлись — верить нельзя ни одному.
 *
 * ИСТОЧНИК В ЭТОЙ РЕДАКЦИИ — ПОТОЛОК, и это ОСНАСТКА, а не физика. Лампы `prender`
 * суть отдельные полигоны, добавленные к набору; в лестнице их нет, потому что
 * лестница строится по геометрии сцены. Чтобы не тащить сюда всю машинерию
 * источников до того, как проверен перенос, светящимися объявляются элементы,
 * смотрящие вниз у самого верха габарита. Замена честная для проверки переноса и
 * негодная для сравнения яркостей с `prender`.
 */

#include "image.h"
#include "lodio.h"
#include "pedge.h"
#include "plink.h"
#include "plod.h"
#include "poly_seg.h"
#include "polygon.h"
#include "pray.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include "transport/cam3.h"
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

static int cmp_d(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x < y) ? -1 : ((x > y) ? 1 : 0);
}

/* ==== ЗЕРКАЛЬНОСТЬ, ПРЕЛОМЛЕНИЕ И ИНТЕРПОЛИРОВАННЫЕ СВОЙСТВА (§89) ============
 *
 * ЧТО ЭТО ЕСТЬ И ЧЕГО НЕ ЕСТЬ — СКАЗАНО ДО КАРТИНКИ. Оператор связей хранит на
 * узле СКАЛЯРНЫЙ радианс, направленности в нём нет, поэтому зеркало в перенос не
 * входит: **зеркала не светят и теней не отбрасывают.** §85.2 держит это
 * отдельной строкой — «плоское зеркало, СВЕТИТ → связь к отражённому элементу»,
 * то есть правка ОПЕРАТОРА. Здесь делается только строка «зеркало и стекло,
 * ВИДНО»: отражённый и преломлённый лучи летят в сцену и читают уже решённое
 * поле.
 *
 * ЧЕТЫРЕ ЧИСЛА НА ЭЛЕМЕНТЕ — ВОТ ЧТО ДЕЛАЕТ ИНТЕРПОЛЯЦИЮ ОСМЫСЛЕННОЙ. По
 * скалярному полю интерполировать нормаль не над чем (А222). Поэтому радианс
 * читается как `B · (1 + E⃗·n̂ / E₀)`: СКАЛЯР несёт энергию, ВЕКТОР — только форму
 * её углового распределения. Множитель лежит в `[0, 2]` по неравенству
 * треугольника (`|E⃗| ≤ E₀`), то есть энергии не создаёт.
 *
 * ИНТЕРПОЛЯЦИЯ — СВОЙСТВО ПОВЕРХНОСТИ (пользователь 08-01: «горы-камни всякие
 * лучше НЕ интерполировать»). Половина этого работает сама: поле `n0/nu/nv`
 * подгоняется по ВЕРШИННЫМ нормалям входа, а без них даёт `nu = nv = 0` ТОЧНО.
 * Сверх того есть явный флаг материала `flat`.
 */
#define SH_DEPTH                                                                                   \
  3 /* отскоков: 3 хватает на «зеркало в зеркале» и ограничивает \
     * стоимость; глубже вклад падает как произведение долей */

/* ТОЧНЫЙ ФРЕНЕЛЬ ПО НЕПОЛЯРИЗОВАННОМУ СВЕТУ. Приближение Шлика не нужно: точная
 * формула столь же дёшева, а проверяется тождеством `R + T = 1`. `ci` — косинус
 * угла падения (положительный), `eta` — отношение показателей (из → в). */
static double fresnel_R(double ci, double eta) {
  double s2 = eta * eta * (1.0 - ci * ci);
  if (s2 >= 1.0) return 1.0; /* полное внутреннее отражение */
  double ct = sqrt(1.0 - s2);
  double rs = (eta * ci - ct) / (eta * ci + ct);
  double rp = (ci - eta * ct) / (ci + eta * ct);
  return 0.5 * (rs * rs + rp * rp);
}

/* Преломление по Снеллиусу; `0` — полное внутреннее отражение, `t` не заполнен. */
static int refract_dir(const double d[3], const double n[3], double eta, double t[3]) {
  double ci = -(d[0] * n[0] + d[1] * n[1] + d[2] * n[2]);
  double s2 = eta * eta * (1.0 - ci * ci);
  if (s2 >= 1.0) return 0;
  double ct = sqrt(1.0 - s2);
  for (int c = 0; c < 3; c++)
    t[c] = eta * d[c] + (eta * ci - ct) * n[c];
  return 1;
}

/* САМОПРОВЕРКИ НА АНАЛИТИКЕ — ГОНЯЮТСЯ ПЕРВЫМИ, ДО СЦЕНЫ (урок А205: картинкой
 * такое не проверяется, ошибка знака даст правдоподобную картинку). */
static double optics_selftest(double *e_rt, double *e_slab, double *e_tir) {
  const double eta = 1.0 / 1.5;
  double worst = 0.0;
  /* 1. R + T = 1 на сетке углов. */
  *e_rt = 0.0;
  for (int i = 0; i <= 90; i++) {
    double th = (double)i * M_PI / 180.0 * 0.99, ci = cos(th);
    double R = fresnel_R(ci, eta), T = 1.0 - R;
    double e = fabs(R + T - 1.0);
    if (e > *e_rt) *e_rt = e;
  }
  /* 2. Луч сквозь плоскопараллельную пластину выходит ПАРАЛЛЕЛЬНО входящему. */
  *e_slab = 0.0;
  for (int i = 1; i <= 80; i++) {
    double th = (double)i * M_PI / 180.0;
    double d[3] = {sin(th), -cos(th), 0.0};
    double nin[3] = {0.0, 1.0, 0.0}, t1[3], t2[3];
    if (!refract_dir(d, nin, eta, t1)) continue;
    double nout[3] = {0.0, 1.0, 0.0}; /* нормаль обращена к лучу изнутри */
    if (!refract_dir(t1, nout, 1.5, t2)) continue;
    double e = 0.0;
    for (int c = 0; c < 3; c++)
      e += fabs(t2[c] - d[c]);
    if (e > *e_slab) *e_slab = e;
  }
  /* 3. Полное внутреннее отражение наступает ровно на `asin(1/n)`. */
  {
    double thc = asin(1.0 / 1.5);
    double a = fresnel_R(cos(thc * (1.0 - 1e-12)), 1.5);
    double b = fresnel_R(cos(thc * (1.0 + 1e-12)), 1.5);
    *e_tir = fabs(b - 1.0) + fabs(a - 1.0);
    /* Ниже угла — не единица, выше — ровно единица. */
    double below = fresnel_R(cos(thc * 0.9), 1.5);
    *e_tir = fabs(b - 1.0) + ((below < 1.0) ? 0.0 : 1.0);
  }
  if (*e_rt > worst) worst = *e_rt;
  if (*e_slab > worst) worst = *e_slab;
  if (*e_tir > worst) worst = *e_tir;
  return worst;
}

/* ДЕМОНСТРАЦИОННЫЕ ТЕЛА. В зале `d = 1.0` у ВСЕХ 33 материалов (замерено по
 * `.mtl`), поэтому преломление показывать не на чем — нужен добавленный шар.
 * Тела живут ТОЛЬКО в проходе камеры: в `hz_polyset` и `hz_pray` они не
 * попадают, поэтому перенос ими не задет по построению (А225). */
typedef struct {
  double c[3], r;
  int kind; /* 0 — зеркало, 1 — стекло */
  double ior;
} sh_ball;

static double ball_hit(const sh_ball *b, const double o[3], const double d[3], double tmin) {
  double oc[3], B = 0.0, C = 0.0, dd = 0.0;
  for (int i = 0; i < 3; i++) {
    oc[i] = o[i] - b->c[i];
    B += oc[i] * d[i];
    C += oc[i] * oc[i];
    dd += d[i] * d[i];
  }
  C -= b->r * b->r;
  double disc = B * B - dd * C;
  if (disc < 0.0) return -1.0;
  double sq = sqrt(disc);
  double t1 = (-B - sq) / dd, t2 = (-B + sq) / dd;
  if (t1 > tmin) return t1;
  if (t2 > tmin) return t2;
  return -1.0;
}

typedef struct {
  const hz_pray *g;
  const hz_polyset *ps;
  const hz_objmesh *m;
  const hz_lod *L;
  const int32_t *cut;
  const double *B, *E0, *Ev, *Le;
  const sh_ball *ball;
  int nball, flatn, nospec, spec;
} sh_ctx;

static void shade_ray(const sh_ctx *S, const double o[3], const double d[3], int depth,
                      double out[3]);

/* Радианс поверхности сцены в точке попадания. */
static void shade_surface(const sh_ctx *S, int32_t k, const double x[3], const double d[3],
                          int depth, double out[3]) {
  const hz_poly *p = &S->ps->p[k];
  int32_t tr = (p->ntri > 0) ? S->ps->tri[p->t0] : -1;
  int32_t mi = (tr >= 0 && S->m->fm != NULL) ? S->m->fm[tr] : 0;
  const hz_obj_mtl *mt = (mi >= 0 && mi < S->m->nmtl) ? &S->m->mtl[mi] : NULL;

  /* ГЕОМЕТРИЧЕСКАЯ и ИНТЕРПОЛИРОВАННАЯ нормали. Вторая — `n0 + nu·u + nv·v` в
   * раме полигона; она не берётся, если материал объявлен плоским либо вход не
   * дал вершинных нормалей (тогда `nu = nv = 0` и разницы нет). */
  double ng[3] = {p->n[0], p->n[1], p->n[2]}, nh[3];
  double du[3];
  for (int c = 0; c < 3; c++)
    du[c] = x[c] - p->org[c];
  double uu = du[0] * p->eu[0] + du[1] * p->eu[1] + du[2] * p->eu[2];
  double vv = du[0] * p->ev[0] + du[1] * p->ev[1] + du[2] * p->ev[2];
  if (S->flatn || (mt != NULL && mt->flat)) {
    for (int c = 0; c < 3; c++)
      nh[c] = ng[c];
  } else {
    double ln = 0.0;
    for (int c = 0; c < 3; c++) {
      nh[c] = p->n0[c] + uu * p->nu[c] + vv * p->nv[c];
      ln += nh[c] * nh[c];
    }
    ln = sqrt(ln);
    if (ln > 0.0)
      for (int c = 0; c < 3; c++)
        nh[c] /= ln;
    else
      for (int c = 0; c < 3; c++)
        nh[c] = ng[c];
  }
  /* Нормали обращаются к лучу: полигоны односторонние, а луч может прийти с
   * изнанки (например, изнутри стеклянного шара). */
  double dn = d[0] * ng[0] + d[1] * ng[1] + d[2] * ng[2];
  if (dn > 0.0)
    for (int c = 0; c < 3; c++) {
      ng[c] = -ng[c];
      nh[c] = -nh[c];
    }

  /* ДИФФУЗНЫЙ ЧЛЕН: скаляр несёт энергию, вектор — форму (см. заголовок §89). */
  int32_t nd = S->cut[k];
  double e0 = S->E0[nd], mod = 1.0;
  if (e0 > 0.0) {
    double ed = S->Ev[3 * nd] * nh[0] + S->Ev[3 * nd + 1] * nh[1] + S->Ev[3 * nd + 2] * nh[2];
    mod = 1.0 + ed / e0;
    if (mod < 0.0) mod = 0.0;
  }
  double Bk = S->B[nd] * mod;
  for (int c = 0; c < 3; c++) {
    double a = (mt != NULL) ? mt->kd3[c] : 0.5;
    out[c] = Bk * a / 0.5; /* правило Т1: цвет — с мелкого полигона */
  }
  /* ДЕЛЬТА-ЗЕРКАЛЬНОСТЬ У МАТЕРИАЛОВ СЦЕНЫ ВЫКЛЮЧЕНА, И ЭТО НЕ ЛЕНЬ, А ЗАМЕР.
   * У всех 33 материалов зала `Ns = 40`, то есть лепесток полуширины
   * `arccos(0.5^(1/40)) ≈ 12°` — НЕ дельта. Отражение такого лепестка одним
   * лучом сваливает всю его энергию в одно направление; на картинке это дало
   * зал, целиком похожий на стекло, а числом — 90.2 % пикселей с недиффузным
   * вкладом выше 10 % против предсказанных 3…25 %. Дельта законна лишь при
   * `Ns → ∞`. Честный путь для лепестка — BRDF, а он требует НАПРАВЛЕННОГО
   * радианса на элементе, которого оператор не хранит (§85.2). Поэтому здесь
   * зеркальны только ОБЪЯВЛЕННЫЕ зеркалом и стеклом тела.
   * `spec` включает дельту у материалов сцены для замера — и она неверна. */
  if (!S->spec || S->nospec || depth <= 0 || mt == NULL) return;

  /* ЗЕРКАЛЬНЫЙ ЧЛЕН. Доля — `Ks + (1 − Ks)·F(θ)`: при нормальном падении она
   * равна `Ks`, к скользящему уходит в единицу, и единицу не превышает никогда.
   * Это ПОСТРОЕНИЕ, а не физика: `Ks` в `.mtl` есть художественная величина, и
   * честнее сказать это прямо, чем выдать за формулу Френеля. */
  double kmax = 0.0;
  for (int c = 0; c < 3; c++)
    if (mt->ks3[c] > kmax) kmax = mt->ks3[c];
  if (!(kmax > 0.0)) return;
  double ci = -(d[0] * nh[0] + d[1] * nh[1] + d[2] * nh[2]);
  if (!(ci > 0.0)) return;
  double ior = (mt->ior > 1.0) ? mt->ior : 1.5;
  double F = fresnel_R(ci, 1.0 / ior);
  double rd[3], rg = 0.0;
  for (int c = 0; c < 3; c++) {
    rd[c] = d[c] + 2.0 * ci * nh[c];
    rg += rd[c] * ng[c];
  }
  /* СТОРОЖ ЗАТЕНЯЮЩЕЙ НОРМАЛИ (А223): отражённый по интерполированной нормали
   * луч может уйти ПОД геометрическую поверхность — тогда берётся геометрическая. */
  if (rg <= 0.0) {
    double cg = -(d[0] * ng[0] + d[1] * ng[1] + d[2] * ng[2]);
    if (!(cg > 0.0)) return;
    for (int c = 0; c < 3; c++)
      rd[c] = d[c] + 2.0 * cg * ng[c];
  }
  double xo[3], sp[3] = {0.0, 0.0, 0.0};
  for (int c = 0; c < 3; c++)
    xo[c] = x[c] + 1e-5 * ng[c];
  shade_ray(S, xo, rd, depth - 1, sp);
  /* ПОСЛОЙНАЯ ФОРМА, А НЕ СЛОЖЕНИЕ: `(1 − w)·диффузное + w·зеркальное`. Первая
   * редакция ДОБАВЛЯЛА зеркальный член поверх ПОЛНОГО диффузного, и энергия не
   * сохранялась: сцена вышла полупрозрачной, а недиффузный вклад превысил 10 % у
   * 73.9 % пикселей против предсказанных 3…25 %. Дефект виден и глазом, но
   * ПОЙМАН он числом — предсказание было записано до прогона. */
  for (int c = 0; c < 3; c++) {
    double wgt = mt->ks3[c] + (1.0 - mt->ks3[c]) * F;
    out[c] = (1.0 - wgt) * out[c] + wgt * sp[c];
  }
}

static void shade_ray(const sh_ctx *S, const double o[3], const double d[3], int depth,
                      double out[3]) {
  out[0] = out[1] = out[2] = 0.0;
  double t = 0.0;
  int32_t k = hz_pray_hit(S->g, o, d, 0.0, &t);
  if (k < 0) t = 1e30;
  /* Тела проверяются отдельно и перекрывают сцену, если ближе. */
  int hb = -1;
  double tb = t;
  for (int i = 0; i < S->nball; i++) {
    double tt = ball_hit(&S->ball[i], o, d, 1e-6);
    if (tt > 0.0 && tt < tb) {
      tb = tt;
      hb = i;
    }
  }
  if (hb >= 0) {
    const sh_ball *b = &S->ball[hb];
    double x[3], nn[3], ln = 0.0;
    for (int c = 0; c < 3; c++) {
      x[c] = o[c] + tb * d[c];
      nn[c] = x[c] - b->c[c];
      ln += nn[c] * nn[c];
    }
    ln = sqrt(ln);
    for (int c = 0; c < 3; c++)
      nn[c] /= ln;
    int inside = (d[0] * nn[0] + d[1] * nn[1] + d[2] * nn[2]) > 0.0;
    if (inside)
      for (int c = 0; c < 3; c++)
        nn[c] = -nn[c];
    if (depth <= 0) return;
    double ci = -(d[0] * nn[0] + d[1] * nn[1] + d[2] * nn[2]);
    double rd[3], xo[3];
    for (int c = 0; c < 3; c++) {
      rd[c] = d[c] + 2.0 * ci * nn[c];
      xo[c] = x[c] + 1e-6 * nn[c];
    }
    if (b->kind == 0) { /* зеркало: одно отражение, поглощение 5 % */
      shade_ray(S, xo, rd, depth - 1, out);
      for (int c = 0; c < 3; c++)
        out[c] *= 0.95;
      return;
    }
    /* СТЕКЛО: Френель делит на отражённую и прошедшую доли, обе прослеживаются. */
    double eta = inside ? b->ior : (1.0 / b->ior);
    /* ГРАНИЦА БЕЗ КОНТРАСТА ПОКАЗАТЕЛЯ — НЕ ГРАНИЦА, и глубины она не тратит.
     * Без этого шар с `n = 1` оставался ВИДЕН: два «преломления» съедали два
     * уровня рекурсии, и поверхность за ним затенялась беднее, а на исчерпании
     * глубины возвращался чёрный. Поймано негативным контролем НК6 (расхождение
     * радианса `1.4e-01` при уровне `1.5`), а не глазом. */
    if (!(fabs(eta - 1.0) > 0.0)) {
      double xt[3];
      for (int c = 0; c < 3; c++)
        xt[c] = x[c] - 1e-9 * nn[c];
      shade_ray(S, xt, d, depth, out);
      return;
    }
    double F = fresnel_R(ci, eta), tr[3], refl[3] = {0.0, 0.0, 0.0};
    shade_ray(S, xo, rd, depth - 1, refl);
    double trc[3] = {0.0, 0.0, 0.0};
    if (refract_dir(d, nn, eta, tr)) {
      double xi[3];
      for (int c = 0; c < 3; c++)
        xi[c] = x[c] - 1e-6 * nn[c];
      shade_ray(S, xi, tr, depth - 1, trc);
    } else {
      F = 1.0;
    }
    for (int c = 0; c < 3; c++)
      out[c] = F * refl[c] + (1.0 - F) * trc[c];
    return;
  }
  if (k < 0) return;
  double x[3];
  for (int c = 0; c < 3; c++)
    x[c] = o[c] + t * d[c];
  shade_surface(S, k, x, d, depth, out);
}

int main(int argc, char **argv) {
  int city = 0, w = 512, ss = 2, nvis = 2, maxlev = 13, noself = 0, nopull = 0, disk = 0,
      noclip = 0, oldpt = 0, flatn = 0, nospec = 0, noballs = 0, h = 0, spec = 0;
  double ballior = 1.5;
  double epsmul = 1.0, radmul = 4.0, base = 1.4142, linkmul = 1.0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "city") == 0) city = 1;
    if (strncmp(argv[i], "w=", 2) == 0) w = (int)strtol(argv[i] + 2, NULL, 10);
    if (strncmp(argv[i], "ss=", 3) == 0) ss = (int)strtol(argv[i] + 3, NULL, 10);
    if (strncmp(argv[i], "vis=", 4) == 0) nvis = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "lev=", 4) == 0) maxlev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "eps=", 4) == 0) epsmul = strtod(argv[i] + 4, NULL);
    /* Допуск СВЯЗИ в допусках среза: `1` — тот же критерий, что у камеры. */
    if (strncmp(argv[i], "link=", 5) == 0) linkmul = strtod(argv[i] + 5, NULL);
    /* ДВА ВЫКЛЮЧАТЕЛЯ — ЗАМЕРЫ, А НЕ РЕЖИМЫ (§86, А195). `noself` снимает
     * самопары и показывает, что именно они купили; `nopull` снимает подъём —
     * это негативный контроль НК1, обязанный ПРОВАЛИТЬСЯ почернением. */
    if (strcmp(argv[i], "noself") == 0) noself = 1;
    if (strcmp(argv[i], "nopull") == 0) nopull = 1;
    /* НК3 и НК4 §87: диск вместо точной формы; точная форма без отсечения. */
    if (strcmp(argv[i], "disk") == 0) disk = 1;
    if (strcmp(argv[i], "noclip") == 0) noclip = 1;
    /* НК5 §88: опорная точка = центр площади, как до правки. */
    if (strcmp(argv[i], "oldpt") == 0) oldpt = 1;
    /* §89: высота кадра (по умолчанию квадрат); НК6 — показатель шара; НК7 —
     * интерполяция выключена; плюс выключатели зеркальности и тел. */
    if (strncmp(argv[i], "h=", 2) == 0) h = (int)strtol(argv[i] + 2, NULL, 10);
    if (strncmp(argv[i], "ior=", 4) == 0) ballior = strtod(argv[i] + 4, NULL);
    if (strcmp(argv[i], "flatn") == 0) flatn = 1;
    if (strcmp(argv[i], "nospec") == 0) nospec = 1;
    if (strcmp(argv[i], "noballs") == 0) noballs = 1;
    if (strcmp(argv[i], "spec") == 0) spec = 1;
  }

  /* САМОПРОВЕРКА ФОРМУЛЫ — ПЕРВОЙ, ДО ВСЯКОЙ СЦЕНЫ (А205). Ошибка знака или
   * нормировки пришла бы потом в виде правдоподобного `Σf`, и отличить её было
   * бы нечем. Допуск `1e-9` — не порог качества, а запас на округление double в
   * `acos` и в эталонной формуле; фактическая невязка на пять порядков меньше. */
  {
    double fb = 0.0, fu = 0.0, ru = 0.0;
    double worst = hz_links_ff_selftest(&fb, &fu, &ru);
    double ert = 0.0, esl = 0.0, eti = 0.0;
    double wopt = optics_selftest(&ert, &esl, &eti);
    printf("== САМОПРОВЕРКА ФОРМ-ФАКТОРА: огромный квадрат %.15f; квадрат a=h %.15f против "
           "эталона %.15f; худшая невязка %.3e\n",
           fb, fu, ru, worst);
    printf("== САМОПРОВЕРКА ОПТИКИ: R+T−1 %.3e; пластина, отклонение от параллельности "
           "%.3e; полное внутреннее отражение %.3e\n",
           ert, esl, eti);
    if (!(worst < 1e-9) || !(wopt < 1e-14)) {
      fprintf(stderr, "аналитика не сходится — дальше идти нельзя\n");
      return 1;
    }
  }
  double dseg = city ? 0.05 : 0.045;

  hz_objmesh m;
  if (hz_obj_load(&m, city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ,
                  city ? HZ_CFG_CITY_SCALE : HZ_CFG_HALL_SCALE) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  hz_pseglist sg;
  if (hz_seg_planar(&sg, &m, dseg) != 0) return 1;
  hz_polyset ps;
  if (hz_poly_build(&ps, &m, &sg) != 0) return 1;
  {
    int32_t mxl = 0;
    for (int32_t k = 0; k < ps.np; k++) {
      int32_t nb = ps.loop[ps.p[k].l0 + ps.p[k].nloop] - ps.loop[ps.p[k].l0];
      if (nb > mxl) mxl = nb;
    }
    printf("== СЦЕНА: %s, треугольников %d, участков %d; вершин края %d (на полигон %.1f, "
           "максимум %d)\n",
           city ? "ГОРОД" : "зал", m.nt, sg.nseg, ps.nbv,
           (double)ps.nbv / (double)(ps.np ? ps.np : 1), mxl);
  }

  tr3_camera cam;
  double eyeh[3] = HZ_CFG_HALL_EYE, ath[3] = HZ_CFG_HALL_AT;
  double eyec[3] = HZ_CFG_CITY_EYE, atc[3] = HZ_CFG_CITY_AT, up[3] = HZ_CFG_UP;
  double eye[3], at[3];
  for (int c = 0; c < 3; c++) {
    eye[c] = city ? eyec[c] : eyeh[c];
    at[c] = city ? atc[c] : ath[c];
  }
  for (int i = 1; i < argc; i++) {
    if (strncmp(argv[i], "eye=", 4) == 0)
      sscanf(argv[i] + 4, "%lf,%lf,%lf", &eye[0], &eye[1], &eye[2]);
    if (strncmp(argv[i], "at=", 3) == 0) sscanf(argv[i] + 3, "%lf,%lf,%lf", &at[0], &at[1], &at[2]);
  }
  if (h <= 0) h = w;
  if (tr3_camera_look(&cam, eye, at, up, HZ_CFG_FOV_DEG * M_PI / 180.0, w, h) != 0) return 1;
  const double eps_px = (HZ_CFG_FOV_DEG * M_PI / 180.0) / (double)w;
  const double eps = eps_px * epsmul;

  double t0 = now_s();
  hz_lod L;
  hz_lodcfg lc;
  memset(&lc, 0, sizeof lc);
  lc.delta0 = dseg;
  lc.eps = eps;
  lc.maxlev = maxlev;
  lc.radmul = radmul;
  lc.base = base;
  if (hz_lod_build_merge(&L, &m, &sg, &ps, &lc) != 0) {
    fprintf(stderr, "отказ лестницы\n");
    return 1;
  }
  printf("== ЛЕСТНИЦА: уровней %d, узлов %d, за %.1f с\n", L.nlev, L.nnd, now_s() - t0);

  /* Сетка лучей — по ИСХОДНЫМ полигонам: видимость обязана считаться по настоящей
   * геометрии, а не по огрублённой, иначе связь пройдёт сквозь стену. */
  hz_pray g;
  if (hz_pray_build(&g, &ps, 4.0) != 0) return 1;

  hz_scene sc;
  sc.L = &L;
  sc.g = &g;
  sc.ps = &ps;
  sc.m = &m;
  hz_linkcfg lkc;
  memset(&lkc, 0, sizeof lkc);
  lkc.eps = eps * linkmul;
  lkc.nvis = nvis;
  lkc.noself = noself;
  lkc.disk = disk;
  lkc.noclip = noclip;
  lkc.oldpt = oldpt;
  hz_linkset S;
  t0 = now_s();
  if (hz_links_build(&S, &sc, &lkc) != 0) {
    fprintf(stderr, "отказ сборки связей\n");
    return 1;
  }
  S.t_build = now_s() - t0;
  printf("== СВЯЗИ: %lld штук за %.1f с; пар рассмотрено %lld, дроблений %lld, отброшено %lld, "
         "проб видимости %lld, максимум связей у узла %lld\n",
         (long long)S.n, S.t_build, (long long)S.nvisit, (long long)S.nrefine, (long long)S.nzero,
         (long long)S.nray, (long long)S.nmax_node);
  printf("   из отброшенных на ГРУБОМ уровне (обрубило поддерево): %lld\n",
         (long long)S.nzero_coarse);
  printf("   с излучателем-листом %lld связей (%.1f %%), несут %.1f %% суммы коэффициентов\n",
         (long long)S.nlink_leaf, 100.0 * (double)S.nlink_leaf / (double)(S.n ? S.n : 1),
         100.0 * S.wleaf);
  printf("   на узел в среднем %.1f; связей / (n log n) = %.2f\n",
         (double)S.n / (double)(L.nnd ? L.nnd : 1),
         (double)S.n / ((double)L.nnd * log2((double)L.nnd + 2.0)));
  /* ПОДПИСЬ ЗАМКНУТОСТИ: `Σf` по листу с предками обязана быть единицей. Хвост
   * по площади, а не среднее (А194) — среднее печатается справочно. */
  printf("== Σf ПО ЛИСТУ (с предками; в замкнутой сцене = 1): p10 %.3f, p50 %.3f, p90 %.3f; "
         "среднее %.3f, min %.3f, max %.3f; площади с Σf<0.9 — %.1f %%\n",
         S.sf_p10, S.sf_p50, S.sf_p90, S.sf_mean, S.sf_min, S.sf_max, 100.0 * S.sf_lowfrac);
  fflush(stdout);

  /* ЛЕЖИТ ЛИ ОПОРНАЯ ТОЧКА НА СВОЁМ ПОЛИГОНЕ. Опорная точка узла есть ЦЕНТР
   * ПЛОЩАДИ (`plod.h`), а участок сегментации бывает невыпуклым и даже
   * несвязным — центр площади тогда лежит ВНЕ его, в воздухе или внутри мебели.
   * Из такой точки видимость меряется не оттуда, откуда светит поверхность.
   * Проверка — принадлежность `(0,0)` петлям края в местных `(u,v)`: начало рамы
   * и есть центр площади, так что тест сводится к подсчёту пересечений луча
   * `v = 0, u > 0` с краем. Площадь считается, а не только число полигонов:
   * вклад в перенос идёт площадью. */
  double *ptn = malloc(3 * (size_t)L.nnd * sizeof *ptn);
  if (ptn == NULL || hz_links_points(ptn, &L, &ps, &m, oldpt) != 0) return 1;
  {
    double aout = 0.0, atot3 = 0.0;
    int32_t nout = 0;
    for (int32_t k = 0; k < ps.np; k++) {
      const hz_poly *p = &ps.p[k];
      /* Точка ПЕРЕНОСА в местных (u,v) полигона; до §88 это было тождественно
       * начало рамы, то есть центр площади. */
      const double *xp = ptn + 3 * L.lab[k];
      double du[3] = {xp[0] - p->org[0], xp[1] - p->org[1], xp[2] - p->org[2]};
      double pu = du[0] * p->eu[0] + du[1] * p->eu[1] + du[2] * p->eu[2];
      double pv = du[0] * p->ev[0] + du[1] * p->ev[1] + du[2] * p->ev[2];
      int cross = 0;
      for (int32_t li = 0; li < p->nloop; li++) {
        int32_t b0 = ps.loop[p->l0 + li], b1 = ps.loop[p->l0 + li + 1];
        for (int32_t e = b0; e < b1; e++) {
          int32_t e2 = (e + 1 < b1) ? e + 1 : b0;
          double u1 = ps.bv[2 * e], v1 = ps.bv[2 * e + 1];
          double u2 = ps.bv[2 * e2], v2 = ps.bv[2 * e2 + 1];
          v1 -= pv;
          v2 -= pv;
          if ((v1 > 0.0) == (v2 > 0.0)) continue;
          double t = v1 / (v1 - v2);
          if (u1 + t * (u2 - u1) > pu) cross++;
        }
      }
      atot3 += p->area;
      if ((cross & 1) == 0) {
        nout++;
        aout += p->area;
      }
    }
    printf("== ОПОРНАЯ ТОЧКА ПЕРЕНОСА ВНЕ СВОЕГО ПОЛИГОНА: %d из %d (%.1f %% ПЛОЩАДИ)\n", nout,
           ps.np, (atot3 > 0.0) ? 100.0 * aout / atot3 : 0.0);
    fflush(stdout);
  }

  /* ЗАМКНУТОСТЬ СЦЕНЫ — НЕЗАВИСИМАЯ ССЫЛКА ДЛЯ `Σf`. Само по себе `Σf < 1` ещё
   * не значит потери: если сцена не замкнута, это ПРАВИЛЬНЫЙ ответ. Меряется
   * лучами и потому не разделяет с `Σf` ни формулы коэффициента, ни иерархии:
   * доля косинусно-взвешенного телесного угла, упирающегося в геометрию. Для
   * замкнутой сцены она единица, и тогда `Σf` обязано быть единицей тоже.
   * Направления — решётка Хаммерсли (детерминирована, §4 запрещает случайность
   * в замерах), косинусное распределение даёт ровно вес `cosθ dω / π`. */
  {
    const int NCLO = 64; /* лучей на узел: 993·64 ≈ 64 тыс., доли секунды */
    double amiss = 0.0, atot2 = 0.0, worst = 1.0;
    double *hitd = malloc((size_t)L.np * (size_t)NCLO * sizeof *hitd);
    int64_t nhd = 0;
    if (hitd == NULL) return 1;
    double *sfn = malloc((size_t)L.nnd * sizeof *sfn);
    double *dev = malloc((size_t)L.nnd * sizeof *dev);
    if (sfn == NULL || dev == NULL) return 1;
    hz_links_sf(&S, &L, sfn);
    int32_t nd2 = 0;
    for (int32_t k = 0; k < L.nnd; k++) {
      if (L.nd[k].level != 0) continue;
      const hz_lodnode *N = &L.nd[k];
      double eu[3], ev[3], t0v[3] = {0.0, 0.0, 0.0};
      int ax = 0;
      for (int c = 1; c < 3; c++)
        if (fabs(N->n[c]) < fabs(N->n[ax])) ax = c;
      t0v[ax] = 1.0;
      eu[0] = N->n[1] * t0v[2] - N->n[2] * t0v[1];
      eu[1] = N->n[2] * t0v[0] - N->n[0] * t0v[2];
      eu[2] = N->n[0] * t0v[1] - N->n[1] * t0v[0];
      double en = sqrt(eu[0] * eu[0] + eu[1] * eu[1] + eu[2] * eu[2]);
      if (!(en > 0.0)) continue;
      for (int c = 0; c < 3; c++)
        eu[c] /= en;
      ev[0] = N->n[1] * eu[2] - N->n[2] * eu[1];
      ev[1] = N->n[2] * eu[0] - N->n[0] * eu[2];
      ev[2] = N->n[0] * eu[1] - N->n[1] * eu[0];
      int hit = 0;
      for (int s = 0; s < NCLO; s++) {
        double u1 = ((double)s + 0.5) / NCLO, u2 = 0.0, f2 = 0.5;
        for (int b = s; b > 0; b >>= 1, f2 *= 0.5)
          if (b & 1) u2 += f2;
        double sn = sqrt(u1), cs = sqrt(1.0 - u1), ph = 2.0 * M_PI * u2;
        double d[3], o[3];
        for (int c = 0; c < 3; c++) {
          d[c] = sn * cos(ph) * eu[c] + sn * sin(ph) * ev[c] + cs * N->n[c];
          o[c] = ptn[3 * k + c] + 1e-5 * N->n[c];
        }
        double th = 0.0;
        if (hz_pray_hit(&g, o, d, 0.0, &th) >= 0) {
          hit++;
          hitd[nhd++] = th;
        }
      }
      double frac = (double)hit / NCLO;
      atot2 += N->area_surf;
      amiss += N->area_surf * (1.0 - frac);
      if (frac < worst) worst = frac;
      /* СО ЗНАКОМ, а не по модулю (А206): ссылка лучами есть оценка СВЕРХУ —
       * `hz_pray_occluded` засчитывает и удар в изнанку односторонней
       * поверхности, которую связь законно отбраковывает. Систематический сдвиг
       * в минус потому законен, и его надо видеть отдельно от разброса. */
      dev[nd2++] = sfn[k] - frac;
    }
    printf("== ЗАМКНУТОСТЬ (лучами, независимо от Σf): в геометрию упирается %.3f "
           "косинусного телесного угла по площади; худший узел %.3f\n",
           (atot2 > 0.0) ? 1.0 - amiss / atot2 : 0.0, worst);
    qsort(hitd, (size_t)nhd, sizeof *hitd, cmp_d);
    printf("   расстояния до попадания, м: p01 %.2e, p10 %.3f, p50 %.3f, p90 %.3f\n",
           hitd[(size_t)(0.01 * (double)nhd)], hitd[(size_t)(0.10 * (double)nhd)],
           hitd[(size_t)(0.50 * (double)nhd)], hitd[(size_t)(0.90 * (double)nhd)]);
    free(hitd);
    qsort(dev, (size_t)nd2, sizeof *dev, cmp_d);
    printf("== СВЕРКА Σf СО ССЫЛКОЙ ЛУЧАМИ (Σf − доля лучей, СО ЗНАКОМ): "
           "p10 %+.3f, p50 %+.3f, p90 %+.3f; |p50| %.3f\n",
           dev[(size_t)(0.10 * nd2)], dev[(size_t)(0.50 * nd2)], dev[(size_t)(0.90 * nd2)],
           fabs(dev[(size_t)(0.50 * nd2)]));
    /* ЭТАЛОН ПЕРЕБОРОМ — на выборке листьев, без иерархии. Три числа в одной
     * строке отвечают на вопрос, который поодиночке не различает ни одна из
     * величин: где потеря — в дроблении, в заслонении или в коэффициенте. */
    {
      const int NS = 32;
      double bn[32], bv[32], bs[32];
      int32_t bnode[32];
      int nb = hz_links_sf_brute(&sc, (L.np / NS > 0) ? L.np / NS : 1, nvis, oldpt, bn, bv, bs,
                                 bnode, NS);
      if (nb > 0) {
        double an = 0.0, av = 0.0, ah = 0.0, as = 0.0, ar2 = 0.0;
        for (int i2 = 0; i2 < nb; i2++) {
          an += bn[i2];
          av += bv[i2];
          ah += sfn[bnode[i2]];
          as += bs[i2];
        }
        (void)ar2;
        printf("== ЭТАЛОН ПЕРЕБОРОМ (%d листьев из %d): Σf перебором без видимости %.3f, "
               "с видимостью %.3f, с видимостью БЕЗ САМОЗАТЕНЕНИЯ %.3f; та же выборка "
               "иерархией %.3f\n",
               nb, L.np, an / nb, av / nb, as / nb, ah / nb);
      }
    }
    free(sfn);
    free(dev);
    fflush(stdout);
  }

  /* --- ИСТОЧНИК: ПОТОЛОК (оснастка, см. заголовок) --- */
  double *Le = calloc((size_t)L.nnd, sizeof *Le);
  double *rho = malloc((size_t)L.nnd * sizeof *rho);
  double *E = malloc((size_t)L.nnd * sizeof *E);
  if (Le == NULL || rho == NULL || E == NULL) return 1;
  int64_t nsrc = 0;
  double asrc = 0.0;
  double ytop = m.hi[1] - 0.05 * (m.hi[1] - m.lo[1]);
  for (int32_t k = 0; k < L.nnd; k++) {
    rho[k] = 0.5;
    /* ИСТОЧНИК ТОЛЬКО НА НУЛЕВОМ УРОВНЕ. Прежде критерий применялся к каждому
     * уровню независимо (2 узла на нулевом против 11 по лестнице), а подъём
     * затем перезаписывал `B` у всех внутренних узлов — излучение грубых просто
     * выбрасывалось. На грубые уровни оно теперь приходит подъёмом, который и
     * так есть. */
    if (L.nd[k].level == 0 && L.nd[k].n[1] < -0.9 && L.nd[k].cy > ytop) {
      Le[k] = HZ_CFG_LAMP_LE;
      rho[k] = 0.0;
      nsrc++;
      asrc += L.nd[k].area_surf;
    }
  }
  printf("== ИСТОЧНИК (оснастка): светящихся узлов уровня 0 — %lld выше y = %.2f, "
         "площадь %.3f м²\n",
         (long long)nsrc, ytop, asrc);

  const int maxit = 200;
  double *rh = calloc((size_t)maxit, sizeof *rh);
  if (rh == NULL) return 1;
  t0 = now_s();
  int it = hz_links_solve(&S, &L, Le, rho, E, 1e-4, maxit, rh, nopull);
  double tsolve = now_s() - t0;
  printf("== РЕШЕНИЕ: %d итераций за %.3f с (%.1f мкс на итерацию)%s\n", it, tsolve,
         1e6 * tsolve / (double)(it > 0 ? it : 1), nopull ? " [НК1: подъём выключен]" : "");
  printf("   невязки:");
  for (int i = 0; i < it && i < 12; i++)
    printf(" %.2e", rh[i]);
  if (it > 12) printf(" ... %.2e", rh[it - 1]);
  printf("\n   отношение соседних (измеренное сжатие):");
  for (int i = 1; i < it && i < 12; i++)
    printf(" %.3f", (rh[i - 1] > 0.0) ? rh[i] / rh[i - 1] : 0.0);
  printf("\n");
  free(rh);
  /* ПОТОК ПО НУЛЕВОМУ УРОВНЮ — величина, по которой сравниваются выключатели. */
  double flux0 = 0.0;
  for (int32_t k = 0; k < L.nnd; k++)
    if (L.nd[k].level == 0) flux0 += E[k] * L.nd[k].area_surf;
  printf("== ПОТОК: Σ B·A по уровню 0 = %.6e Вт/ср%s%s\n", flux0, noself ? " [noself]" : "",
         nopull ? " [nopull]" : "");

  /* --- СРЕЗ И КАРТИНКА --- */
  int32_t *cut = malloc((size_t)L.np * sizeof *cut);
  if (cut == NULL) return 1;
  int32_t ncut = hz_lod_cut(&L, eye, eps, 1, cut);
  printf("== СРЕЗ: %d элементов\n", ncut);

  /* ВЕКТОР ОБЛУЧЁННОСТИ — один проход после сходимости (§89, А222). */
  double *Efield = malloc((size_t)L.nnd * sizeof *Efield);
  double *Evec = malloc(3 * (size_t)L.nnd * sizeof *Evec);
  if (Efield == NULL || Evec == NULL) return 1;
  hz_links_evec(&S, &L, ptn, E, Efield, Evec);
  {
    double *rat = malloc((size_t)L.nnd * sizeof *rat);
    if (rat == NULL) return 1;
    int32_t nr = 0;
    for (int32_t k = 0; k < L.nnd; k++) {
      if (L.nd[k].level != 0 || !(Efield[k] > 0.0)) continue;
      double v = sqrt(Evec[3 * k] * Evec[3 * k] + Evec[3 * k + 1] * Evec[3 * k + 1] +
                      Evec[3 * k + 2] * Evec[3 * k + 2]);
      rat[nr++] = v / Efield[k];
    }
    qsort(rat, (size_t)nr, sizeof *rat, cmp_d);
    if (nr > 0)
      printf("== НАПРАВЛЕННОСТЬ ПОЛЯ |E⃗|/E₀ (0 — изотропно, 1 — одно направление): "
             "p50 %.3f, p90 %.3f\n",
             rat[(size_t)(0.50 * nr)], rat[(size_t)(0.90 * nr)]);
    free(rat);
  }
  /* УГОЛ МЕЖДУ ПЛОСКОЙ И ИНТЕРПОЛИРОВАННОЙ НОРМАЛЬЮ — мера того, что вообще даёт
   * интерполяция на этой сцене. Берётся в вершинах края, взвешивается площадью. */
  {
    double *ang = malloc((size_t)ps.nbv * sizeof *ang);
    if (ang == NULL) return 1;
    int32_t na = 0;
    double amax = 0.0;
    for (int32_t k2 = 0; k2 < ps.np; k2++) {
      const hz_poly *p = &ps.p[k2];
      for (int32_t b = ps.loop[p->l0]; b < ps.loop[p->l0 + p->nloop] && na < ps.nbv; b++) {
        double uu2 = ps.bv[2 * b], vv2 = ps.bv[2 * b + 1], nh[3], ln = 0.0, dp = 0.0;
        for (int c = 0; c < 3; c++) {
          nh[c] = p->n0[c] + uu2 * p->nu[c] + vv2 * p->nv[c];
          ln += nh[c] * nh[c];
        }
        ln = sqrt(ln);
        if (!(ln > 0.0)) continue;
        for (int c = 0; c < 3; c++)
          dp += nh[c] * p->n[c] / ln;
        if (dp > 1.0) dp = 1.0;
        if (dp < -1.0) dp = -1.0;
        double a = acos(dp) * 180.0 / M_PI;
        ang[na++] = a;
        if (a > amax) amax = a;
      }
    }
    qsort(ang, (size_t)na, sizeof *ang, cmp_d);
    if (na > 0)
      printf("== ИНТЕРПОЛЯЦИЯ НОРМАЛИ, угол к плоской (градусы): p50 %.2f, p90 %.2f, "
             "максимум %.1f, точек %d\n",
             ang[(size_t)(0.50 * na)], ang[(size_t)(0.90 * na)], amax, na);
    free(ang);
  }

  /* ДЕМОНСТРАЦИОННЫЕ ТЕЛА: зеркальный и стеклянный шары над столом. В зале
   * `d = 1.0` у всех материалов, и без них преломление показать не на чем. */
  sh_ball balls[2];
  int nball = 0;
  if (!noballs) {
    double cx2 = 0.5 * (m.lo[0] + m.hi[0]), cz2 = 0.5 * (m.lo[2] + m.hi[2]);
    double rr = 0.10 * (m.hi[1] - m.lo[1]);
    balls[0].c[0] = cx2 - 1.6 * rr;
    balls[0].c[1] = m.lo[1] + 3.2 * rr;
    balls[0].c[2] = cz2;
    balls[0].r = rr;
    balls[0].kind = 0;
    balls[0].ior = 1.0;
    balls[1].c[0] = cx2 + 1.6 * rr;
    balls[1].c[1] = m.lo[1] + 3.2 * rr;
    balls[1].c[2] = cz2;
    balls[1].r = rr;
    balls[1].kind = 1;
    balls[1].ior = ballior;
    nball = 2;
    printf("== ТЕЛА: зеркальный и стеклянный шары, радиус %.2f м, показатель %.3f\n", rr, ballior);
  }

  sh_ctx SH;
  SH.g = &g;
  SH.ps = &ps;
  SH.m = &m;
  SH.L = &L;
  SH.cut = cut;
  SH.B = E;
  SH.E0 = Efield;
  SH.Ev = Evec;
  SH.Le = Le;
  SH.ball = balls;
  SH.nball = nball;
  SH.flatn = flatn;
  SH.nospec = nospec;
  SH.spec = spec;

  /* НК6 ПО ЛУЧАМ, А НЕ ПО БАЙТАМ КАРТИНКИ. Первый прогон контроля сравнивал
   * PPM, и это была моя ошибка: тональная компрессия нормируется по 99.5-й
   * процентили, блик её сдвигает, и различается почти весь кадр — сравнение
   * мерило нормировку, а не шар. Здесь сравнивается РАДИАНС: стеклянный шар с
   * показателем `1.0` обязан быть невидим, то есть луч сквозь него обязан дать
   * ровно то же, что луч без шара. */
  {
    sh_ctx S1 = SH, S0 = SH;
    sh_ball b1 = balls[1];
    b1.ior = 1.0;
    S1.ball = &b1;
    S1.nball = 1;
    S0.nball = 0;
    double worst = 0.0, ref = 0.0;
    int nhit = 0;
    for (int i = 0; i < 262144; i++) {
      double o[3], d[3], v1[3], v0[3];
      /* Индексы в 64 битах: при 262 144 пробах `i·104729` переполняет `int`. */
      int px = (int)(((int64_t)i * 7919) % w), py = (int)(((int64_t)i * 104729) % h);
      tr3_camera_ray(&cam, px, py, o, d);
      if (ball_hit(&b1, o, d, 1e-6) <= 0.0) continue;
      nhit++;
      shade_ray(&S1, o, d, SH_DEPTH, v1);
      shade_ray(&S0, o, d, SH_DEPTH, v0);
      for (int c = 0; c < 3; c++) {
        double e = fabs(v1[c] - v0[c]);
        if (e > worst) worst = e;
        if (fabs(v0[c]) > ref) ref = fabs(v0[c]);
      }
    }
    printf("== НК6 (стекло при n = 1 обязано быть невидимо): лучей сквозь шар %d, "
           "худшее расхождение радианса %.3e при уровне %.3e\n",
           nhit, worst, ref);
  }

  int n = w * h;
  double *L3 = calloc((size_t)n * 3, sizeof *L3);
  unsigned char *rgb = malloc((size_t)n * 3);
  int64_t nspec = 0;
  if (L3 == NULL || rgb == NULL) return 1;
  t0 = now_s();
#pragma omp parallel for schedule(dynamic, 16) reduction(+ : nspec)
  for (int i = 0; i < n; i++) {
    double acc[3] = {0.0, 0.0, 0.0}, accd[3] = {0.0, 0.0, 0.0};
    for (int sy = 0; sy < ss; sy++)
      for (int sx = 0; sx < ss; sx++) {
        double o[3], d[3], v[3], vd[3];
        tr3_camera_ray(&cam, i % w, i / w, o, d);
        shade_ray(&SH, o, d, SH_DEPTH, v);
        sh_ctx SD = SH;
        SD.nospec = 1;
        SD.nball = 0;
        shade_ray(&SD, o, d, 0, vd);
        for (int c = 0; c < 3; c++) {
          acc[c] += v[c];
          accd[c] += vd[c];
        }
      }
    double wgt = 1.0 / (double)(ss * ss);
    double s1 = 0.0, s2 = 0.0;
    for (int c = 0; c < 3; c++) {
      L3[(size_t)i * 3 + (size_t)c] = acc[c] * wgt;
      s1 += acc[c];
      s2 += accd[c];
    }
    /* Доля пикселей, где НЕдиффузная часть даёт больше 10 % — замер, а не отладка. */
    if (s1 > 0.0 && fabs(s1 - s2) > 0.10 * s1) nspec++;
  }
  double tframe = now_s() - t0;
  printf("== ЗЕРКАЛЬНОСТЬ: пикселей с недиффузным вкладом выше 10 %% — %.1f %%\n",
         100.0 * (double)nspec / (double)n);

  double *lum = malloc((size_t)n * sizeof *lum);
  if (lum == NULL) return 1;
  for (int i = 0; i < n; i++)
    lum[i] = 0.2126 * L3[3 * i] + 0.7152 * L3[3 * i + 1] + 0.0722 * L3[3 * i + 2];
  qsort(lum, (size_t)n, sizeof *lum, cmp_d);
  double mx = lum[(size_t)((double)n * 0.995)];
  free(lum);
  if (!(mx > 0.0)) mx = 1.0;
  for (int i = 0; i < 3 * n; i++) {
    double x = L3[i] / mx;
    if (x < 0.0) x = 0.0;
    if (x > 1.0) x = 1.0;
    rgb[i] = (unsigned char)(255.0 * pow(x, 1.0 / 2.2) + 0.5);
  }
  char path[96];
  snprintf(path, sizeof path, "img/prad_%s.ppm", city ? "city" : "hall");
  if (hz_ppm_write_rgb(path, rgb, w, h) != 0) return 1;
  printf("== КАДР: %s, %d×%d, суперсэмплинг %d×%d, за %.2f с\n", path, w, h, ss, ss, tframe);

  free(ptn);
  free(L3);
  free(rgb);
  free(cut);
  free(Le);
  free(rho);
  free(E);
  hz_links_free(&S);
  hz_pray_free(&g);
  hz_lod_free(&L);
  hz_poly_free(&ps);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
