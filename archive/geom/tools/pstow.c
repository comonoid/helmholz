/* pstow.c — ШАГ О70 (Ф1): УКЛАДКА РЕЗКОЙ. План — §244, «ПОПРАВЛЕННЫЙ ПЛАН».
 *
 * ЧТО МЕРИТСЯ. Не фронт и не кадр, а УКЛАДКА, на которой фронт стоит: во что
 * обходится резать входную геометрию по коробке ячейки вместо того, чтобы
 * оставлять не влезшее у внутренних узлов (§193: 71…79 %) или ссылаться на него
 * из каждого листа (§95.1: избыточность 1193.36, сборка 361.3 с).
 *
 * ТРИ ВЕЩИ, КОТОРЫЕ ЗДЕСЬ РАЗДЕЛЕНЫ НАМЕРЕННО.
 *   ДЕРЕВО   строится КАМЕРОНЕЗАВИСИМО (А470): спуск до `maxlev`/`leafmax`,
 *            плюс градуировка 2:1. Дерево, построенное под камеру, пришлось бы
 *            перестраивать на каждое её движение, и вся инкрементность,
 *            ради которой структура выбрана, исчезла бы.
 *   СРЕЗ     выбирается ПРИ ЧТЕНИИ, по полу в пикселях (А133): узел берётся,
 *            когда его угловой размер упал до `px` пикселей.
 *   КУСКИ    выводятся отсечением НА СРЕЗЕ и НЕ ХРАНЯТСЯ (А472).
 *
 * ЗАПУСК:
 *   pstow ФАЙЛ.obj МАСШТАБ [leaf=N] [lev=N] [grade=0|1] [eye=x,y,z]
 *         [px=8,4,2,1,0] [naive=1] [strict=1] [vtri=N]
 *   `px=0` — срез ДО ЛИСТЬЕВ; `strict=1` — негативный контроль (старая укладка
 *   «влезает целиком»); `naive=1` — негативный контроль над счётчиком побитовых
 *   расхождений вершины.
 */
#include "pclip.h"
#include "ptree.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ДОПУСК СЛИЯНИЯ КОМПЛАНАРНЫХ (А473), В ДОЛЯХ РЕБРА ЯЧЕЙКИ СРЕЗА. Абсолютного
 * числа здесь быть не может: ячейка среза по построению имеет размер в
 * несколько ПИКСЕЛЕЙ (§241.4), и различимость плоскостей задаётся ею, а не
 * метром. Шестнадцатая доля ячейки — заведомо ниже пикселя (то есть ниже того,
 * что фронт вообще разрешает) и на четырнадцать порядков выше шума подгонки
 * плоскости в двойной точности. Число НЕ подпирает вывод в одиночку: `k`
 * печатается сразу при трёх допусках, вчетверо грубее и вчетверо тоньше. */
#define HZ_PSTOW_PLANE_REL 0.0625
/* ПРЕДЕЛ ЧИСЛА ГРУПП НА ЯЧЕЙКУ. Слияние идёт перебором по представителям, то
 * есть O(k·g); ячейка с тысячей РАЗНЫХ плоскостей есть скопление по любому
 * критерию, и точное `g` для неё ничего не решает. Срабатывания СЧИТАЮТСЯ и
 * печатаются: молчаливый предел читался бы как «столько и было». */
#define HZ_PSTOW_GMAX 1024
/* СКОЛЬКО ТРЕУГОЛЬНИКОВ ИДЁТ В СВЕРКУ ВЕРШИН. Полная сверка требует таблицы на
 * все порождённые вершины; при избыточности 10³ это сотни миллионов ключей.
 * Берётся ДЕТЕРМИНИРОВАННАЯ выборка `t % шаг == 0`, шаг печатается. */
#define HZ_PSTOW_VTRI 200000
#define HZ_PSTOW_VBITS 23 /* 2^23 гнёзд таблицы сверки (≈470 МБ) */
/* СКОЛЬКО ЯЧЕЕК УРОВНЯ ДОЛЖНО НАБРАТЬСЯ, ЧТОБЫ ДОЛЯ С НЕГО ШЛА В НАКЛОН. */
#define HZ_PSTOW_LEVMIN 64

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* ---- ТАБЛИЦА СВЕРКИ ВЕРШИН (А474) ---------------------------------------- */
/* Ключ — ПРОИСХОЖДЕНИЕ вершины (треугольник плюс набор секущих плоскостей).
 * Утверждение, которое проверяется: одинаковый ключ обязан давать ПОБИТОВО
 * одинаковые координаты, из какой бы ячейки и с какого бы уровня вершина ни
 * считалась. Это и есть лечение Г50 в самом лёгком виде. */
typedef struct {
  int32_t tri;
  unsigned char used, kind, e, a1, a2, lev;
  double c1, c2, v[3];
} vent;

typedef struct {
  vent *e;
  int64_t cap, n;
  int64_t nshared;     /* ключей, встреченных больше одного раза */
  int64_t nshared_lev; /* из них — с РАЗНЫХ уровней дерева */
  int64_t nmis;        /* побитовых расхождений */
  int64_t nmis_lev;    /* из них на перепаде уровней */
  int64_t nfull;       /* вершин, не поместившихся в таблицу */
  int32_t stride;
} vmap;

static uint64_t v_mix(uint64_t x) {
  x ^= x >> 33;
  x *= 0xff51afd7ed558ccdULL;
  x ^= x >> 33;
  x *= 0xc4ceb9fe1a85ec53ULL;
  x ^= x >> 33;
  return x;
}

static int vmap_init(vmap *M, int32_t stride) {
  memset(M, 0, sizeof *M);
  M->cap = (int64_t)1 << HZ_PSTOW_VBITS;
  M->e = calloc((size_t)M->cap, sizeof *M->e);
  M->stride = stride;
  return M->e == NULL ? 1 : 0;
}

static void vmap_free(vmap *M) {
  free(M->e);
  memset(M, 0, sizeof *M);
}

static void vmap_put(vmap *M, int32_t tri, const hz_pclip_prov *p, const double *v, int lev) {
  if (M->e == NULL || p->kind == 0 || p->kind == 3) return;
  if (M->n * 10 > M->cap * 7) {
    M->nfull++;
    return;
  }
  uint64_t c1b, c2b;
  memcpy(&c1b, &p->c1, sizeof c1b);
  memcpy(&c2b, &p->c2, sizeof c2b);
  uint64_t h = v_mix((uint64_t)(uint32_t)tri * 0x9e3779b97f4a7c15ULL + p->kind * 131ULL +
                     p->e * 17ULL + p->a1 * 7ULL + p->a2 * 3ULL) ^
               v_mix(c1b) ^ v_mix(c2b + 0x1234567ULL);
  int64_t i = (int64_t)(h & (uint64_t)(M->cap - 1));
  for (;;) {
    vent *E = &M->e[i];
    if (!E->used) {
      E->used = 1;
      E->tri = tri;
      E->kind = p->kind;
      E->e = p->e;
      E->a1 = p->a1;
      E->a2 = p->a2;
      E->c1 = p->c1;
      E->c2 = p->c2;
      E->lev = (unsigned char)lev;
      memcpy(E->v, v, 3 * sizeof *v);
      M->n++;
      return;
    }
    if (E->tri == tri && E->kind == p->kind && E->e == p->e && E->a1 == p->a1 && E->a2 == p->a2 &&
        memcmp(&E->c1, &p->c1, sizeof E->c1) == 0 && memcmp(&E->c2, &p->c2, sizeof E->c2) == 0) {
      M->nshared++;
      int dl = (E->lev != (unsigned char)lev);
      if (dl) M->nshared_lev++;
      if (memcmp(E->v, v, 3 * sizeof *v) != 0) {
        M->nmis++;
        if (dl) M->nmis_lev++;
      }
      return;
    }
    i = (i + 1) & (M->cap - 1);
  }
}

/* ---- СЧЁТЧИКИ СРЕЗА ------------------------------------------------------- */
#define KHIST 4096

typedef struct {
  int64_t ncell, ncell_geo, ncand, npiece, nempty, ndegen, nfall;
  int maxnv;
  int64_t nvhist[HZ_PCLIP_MAXV + 1];
  double area;
  int64_t cell_lev[32], multi_lev[32], piece_lev[32];
  int64_t khist[KHIST + 1], phist[KHIST + 1]; /* по ПЛОСКОСТЯМ и по КУСКАМ */
  int64_t kmax, pmax, ksum, psum;
  int64_t kmax_lo, kmax_hi; /* тот же max при допуске вчетверо грубее/тоньше */
  int64_t ksum_lo, ksum_hi;
  int64_t nmulti, nsimple, nsimple_ix, gcap_hit;
  int32_t levmin, levmax;
} cutstat;

/* ---- РАБОЧИЕ БУФЕРЫ ЯЧЕЙКИ ------------------------------------------------ */
typedef struct {
  int32_t *tri;      /* исходный треугольник куска */
  double *pl;        /* 4 числа на кусок: единичная нормаль и смещение */
  int32_t *grp[3];   /* группа куска при трёх допусках */
  int32_t *rep;      /* представители групп (номер куска) */
  int32_t *uf, *uf2; /* объединение-поиск по группам: по КООРДИНАТАМ и по ИНДЕКСАМ */
  uint64_t *ekey;    /* хеш рёбер: ключ */
  int32_t *eval;     /* хеш рёбер: кусок */
  int32_t cap, ecap;
} scratch;

static void sc_free(scratch *S);

/* НА ЛЮБОМ ОТКАЗЕ РОСТА БУФЕР ОСВОБОЖДАЕТСЯ ЦЕЛИКОМ. Это не аккуратность ради
 * аккуратности: `sc_reserve` есть ОБЁРТКА НАД ВЫДЕЛЕНИЕМ, а на них
 * gcc-анализатор теряет связь «указатель уложен в структуру» и сообщает о
 * фантомной утечке (класс описан в CLAUDE.md по разбору diam 07-24). Полная
 * очистка на отказе снимает находку ПО СУЩЕСТВУ, а не подавлением: отказ здесь
 * смертелен для замера, и держать половину буферов незачем. */
static int sc_reserve(scratch *S, int32_t n) {
  if (n <= S->cap) return 0;
  int32_t nc = S->cap > 0 ? S->cap : 1024;
  while (nc < n)
    nc *= 2;
  int32_t ec = nc * 8;
  int fail = 0;
  /* Указатель кладётся в структуру СРАЗУ, а признак отказа копится отдельно:
   * так у каждого выделенного блока ровно один владелец на всех путях, и
   * рассуждать о нём (и анализатору, и человеку) можно по одному месту.
   * Старый блок при неудачном `realloc` остаётся жив и лежит там же. */
  void *q;
  q = realloc(S->tri, (size_t)nc * sizeof *S->tri);
  if (q == NULL)
    fail = 1;
  else
    S->tri = q;
  q = realloc(S->pl, 4 * (size_t)nc * sizeof *S->pl);
  if (q == NULL)
    fail = 1;
  else
    S->pl = q;
  for (int j = 0; j < 3; j++) {
    q = realloc(S->grp[j], (size_t)nc * sizeof *S->grp[j]);
    if (q == NULL)
      fail = 1;
    else
      S->grp[j] = q;
  }
  /* по представителю на группу И НА КАЖДЫЙ ИЗ ТРЁХ ДОПУСКОВ */
  q = realloc(S->rep, 3 * (size_t)nc * sizeof *S->rep);
  if (q == NULL)
    fail = 1;
  else
    S->rep = q;
  q = realloc(S->uf, (size_t)nc * sizeof *S->uf);
  if (q == NULL)
    fail = 1;
  else
    S->uf = q;
  q = realloc(S->uf2, (size_t)nc * sizeof *S->uf2);
  if (q == NULL)
    fail = 1;
  else
    S->uf2 = q;
  q = realloc(S->ekey, (size_t)ec * sizeof *S->ekey);
  if (q == NULL)
    fail = 1;
  else
    S->ekey = q;
  q = realloc(S->eval, (size_t)ec * sizeof *S->eval);
  if (q == NULL)
    fail = 1;
  else
    S->eval = q;
  if (fail) {
    sc_free(S);
    return 1;
  }
  S->cap = nc;
  S->ecap = ec;
  return 0;
}

static void sc_free(scratch *S) {
  free(S->tri);
  free(S->pl);
  for (int q = 0; q < 3; q++)
    free(S->grp[q]);
  free(S->rep);
  free(S->uf);
  free(S->uf2);
  free(S->ekey);
  free(S->eval);
  memset(S, 0, sizeof *S);
}

static int32_t uf_find(int32_t *u, int32_t x) {
  while (u[x] != x)
    x = u[x] = u[u[x]];
  return x;
}

/* Наибольшее расхождение двух плоскостей ПО КОРОБКЕ ЯЧЕЙКИ. Критерий один и
 * связывает нормаль со смещением (poly_seg.h, урок 1): раздельного допуска по
 * нормали нет, потому что при фиксированном допуске по нормали кривая
 * поверхность не сливалась бы ни на какой дальности. */
static double pl_gap(const double *A, const double *B, const double *ctr, const double *hlf) {
  double dn[3], dd = A[3] - B[3], s = 0.0;
  for (int c = 0; c < 3; c++) {
    dn[c] = A[c] - B[c];
    s += fabs(dn[c]) * hlf[c];
  }
  double at = dn[0] * ctr[0] + dn[1] * ctr[1] + dn[2] * ctr[2] - dd;
  return fabs(at) + s;
}

/* Задевает ли отрезок коробку — параметрическое отсечение по трём слоям. Нужен
 * критерию связности: рёбра считаются общими ВНУТРИ ячейки (А130), иначе два
 * куска, чей общий стык лежит за её пределами, объявятся сросшимися. */
static int seg_hits_box(const double *P, const double *Q, const double *lo, const double *hi) {
  double t0 = 0.0, t1 = 1.0;
  for (int c = 0; c < 3; c++) {
    double d = Q[c] - P[c];
    if (fabs(d) < 1e-300) {
      if (P[c] < lo[c] || P[c] > hi[c]) return 0;
      continue;
    }
    double a = (lo[c] - P[c]) / d, b = (hi[c] - P[c]) / d;
    if (a > b) {
      double s = a;
      a = b;
      b = s;
    }
    if (a > t0) t0 = a;
    if (b < t1) t1 = b;
    if (t0 > t1) return 0;
  }
  return 1;
}

/* ---- ОБРАБОТКА ОДНОЙ ЯЧЕЙКИ СРЕЗА ---------------------------------------- */
static int cell_do(const hz_objmesh *m, const double *tpl, const double *lo, const double *hi,
                   const int32_t *list, int32_t n, int lev, cutstat *S, scratch *W, vmap *M,
                   int naive) {
  S->ncell++;
  if (lev < S->levmin) S->levmin = lev;
  if (lev > S->levmax) S->levmax = lev;
  S->ncand += n;
  if (n <= 0) return 0;
  if (sc_reserve(W, n) != 0) return 1;
  double ctr[3], hlf[3];
  for (int c = 0; c < 3; c++) {
    ctr[c] = 0.5 * (lo[c] + hi[c]);
    hlf[c] = 0.5 * (hi[c] - lo[c]);
  }
  double tol[3];
  tol[0] = HZ_PSTOW_PLANE_REL * 0.25 * (hi[0] - lo[0]);
  tol[1] = HZ_PSTOW_PLANE_REL * (hi[0] - lo[0]);
  tol[2] = HZ_PSTOW_PLANE_REL * 4.0 * (hi[0] - lo[0]);
  int32_t np = 0, ng[3] = {0, 0, 0};
  for (int32_t i = 0; i < n; i++) {
    int32_t t = list[i];
    const double *A = m->v + 3 * (size_t)m->f[3 * (size_t)t + 0];
    const double *B = m->v + 3 * (size_t)m->f[3 * (size_t)t + 1];
    const double *C = m->v + 3 * (size_t)m->f[3 * (size_t)t + 2];
    /* ДЕШЁВЫЙ ОТКАЗ ПО ПЛОСКОСТИ. Коробка треугольника задевает ячейку у ста
     * кандидатов из ста одного (замерено на зале: 802 млн пустых отсечений
     * против 7 млн кусков), а плоскость отсекает их одним сравнением. Это не
     * приближение: если коробка целиком по одну сторону плоскости, куска нет
     * ТОЧНО. */
    const double *tp = tpl + 4 * (size_t)t;
    double sd = tp[0] * ctr[0] + tp[1] * ctr[1] + tp[2] * ctr[2] - tp[3];
    double rr = fabs(tp[0]) * hlf[0] + fabs(tp[1]) * hlf[1] + fabs(tp[2]) * hlf[2];
    if (fabs(sd) > rr) {
      S->nempty++;
      continue;
    }
    hz_pclip_poly P;
    hz_pclip_tri_ex(A, B, C, lo, hi, &P, naive);
    S->nfall += P.nfall;
    if (P.nv < 3) {
      S->nempty++;
      continue;
    }
    double ar = hz_pclip_area(&P);
    if (!(ar > 0.0)) {
      /* Вырожденный кусок отбрасывается ЯВНО и считается (А475): молчаливый
       * отброс здесь и есть потерянный заслон, который не находится картинкой. */
      S->ndegen++;
      continue;
    }
    if (P.nv > S->maxnv) S->maxnv = P.nv;
    S->nvhist[P.nv]++;
    S->area += ar;
    S->npiece++;
    if (M != NULL && M->stride > 0 && t % M->stride == 0)
      for (int q = 0; q < P.nv; q++)
        vmap_put(M, t, &P.p[q], P.v[q], lev);
    /* плоскость куска — плоскость ИСХОДНОГО треугольника, единичная */
    double *pl = W->pl + 4 * (size_t)np;
    memcpy(pl, tp, 4 * sizeof *pl);
    W->tri[np] = t;
    /* СЛИЯНИЕ КОМПЛАНАРНЫХ. Знак нормали проверяется, величина — нет: иначе
     * передняя и задняя грани тонкой панели слились бы в одну плоскость. */
    for (int q = 0; q < 3; q++) {
      int32_t g = -1;
      for (int32_t j = 0; j < ng[q] && g < 0; j++) {
        const double *R = W->pl + 4 * (size_t)W->rep[j * 3 + q];
        if (R[0] * pl[0] + R[1] * pl[1] + R[2] * pl[2] <= 0.0) continue;
        if (pl_gap(R, pl, ctr, hlf) <= tol[q]) g = j;
      }
      if (g < 0) {
        if (ng[q] < HZ_PSTOW_GMAX) {
          W->rep[ng[q] * 3 + q] = np;
          g = ng[q]++;
        } else {
          if (q == 1) S->gcap_hit++;
          g = ng[q] - 1;
        }
      }
      W->grp[q][np] = g;
    }
    np++;
  }
  if (np == 0) return 0;
  S->ncell_geo++;
  S->cell_lev[lev & 31]++;
  S->piece_lev[lev & 31] += np;
  int32_t k = ng[1];
  S->ksum += k;
  S->ksum_lo += ng[0];
  S->ksum_hi += ng[2];
  if (k > S->kmax) S->kmax = k;
  if (ng[0] > S->kmax_lo) S->kmax_lo = ng[0];
  if (ng[2] > S->kmax_hi) S->kmax_hi = ng[2];
  S->khist[k < KHIST ? k : KHIST]++;
  S->psum += np;
  if (np > S->pmax) S->pmax = np;
  S->phist[np < KHIST ? np : KHIST]++;
  if (k < 2) return 0;
  /* МНОГОПЛОСКОСТНАЯ ЯЧЕЙКА. Считать «многокусковость» по числу ФРАГМЕНТОВ
   * нельзя (А473): плоская стена из тысячи треугольников дала бы k = 1000 и
   * прошла бы за излом. */
  S->nmulti++;
  S->multi_lev[lev & 31]++;
  for (int32_t j = 0; j < k; j++) {
    W->uf[j] = j;
    W->uf2[j] = j;
  }
  int32_t eh = 1;
  while (eh < np * 8)
    eh *= 2;
  if (eh > W->ecap) eh = W->ecap;
  for (int32_t j = 0; j < eh; j++)
    W->eval[j] = -1;
  /* СВЯЗНОСТЬ ПО ОБЩИМ РЁБРАМ ВНУТРИ ЯЧЕЙКИ (§241.6, критерий А130: по рёбрам,
   * а не по вершинам — точечное касание не есть связность).
   *
   * КЛЮЧ РЕБРА — КООРДИНАТЫ, А НЕ ИНДЕКСЫ ВЕРШИН, и это не мелочь. Урок 2
   * `poly_seg.h`: смежность по общим индексам меряет НАРЕЗКУ МЕША, а не сцену —
   * у конференц-зала вход распадается на 3 603 связных компоненты, потому что
   * шов записан удвоенными вершинами. Ключ по координатам ТОЧНЫЙ (побитовый,
   * без допуска) и через такой шов проходит. Тот же счёт ПО ИНДЕКСАМ ведётся
   * рядом: разница двух чисел и есть цена нарезки входа. */
  for (int32_t j = 0; j < np; j++) {
    int32_t t = W->tri[j];
    const int32_t *fv = m->f + 3 * (size_t)t;
    for (int q = 0; q < 3; q++) {
      int32_t va = fv[q], vb = fv[(q + 1) % 3];
      const double *PA = m->v + 3 * (size_t)va, *PB = m->v + 3 * (size_t)vb;
      /* канонический порядок концов — по битам координат, а не по номеру */
      if (memcmp(PA, PB, 3 * sizeof *PA) > 0) {
        const double *S0 = PA;
        PA = PB;
        PB = S0;
      }
      uint64_t key = 0;
      for (int c = 0; c < 3; c++) {
        uint64_t b1, b2;
        memcpy(&b1, PA + c, sizeof b1);
        memcpy(&b2, PB + c, sizeof b2);
        key = v_mix(key ^ v_mix(b1)) + v_mix(b2);
      }
      uint64_t kx = (uint64_t)(uint32_t)(va < vb ? va : vb) * 2654435761ULL +
                    (uint64_t)(uint32_t)(va < vb ? vb : va);
      int32_t h = (int32_t)(v_mix(key) & (uint64_t)(eh - 1));
      for (;;) {
        if (W->eval[h] < 0) {
          W->eval[h] = j;
          W->ekey[h] = key;
          break;
        }
        if (W->ekey[h] == key) {
          int32_t o = W->eval[h];
          if (W->grp[1][o] != W->grp[1][j] && seg_hits_box(PA, PB, lo, hi)) {
            int32_t ra = uf_find(W->uf, W->grp[1][o]), rb = uf_find(W->uf, W->grp[1][j]);
            if (ra != rb) W->uf[ra] = rb;
            /* тот же союз, но только если ребро общее И ПО ИНДЕКСАМ */
            const int32_t *fo = m->f + 3 * (size_t)W->tri[o];
            int same = 0;
            for (int z = 0; z < 3; z++) {
              int32_t oa = fo[z], ob = fo[(z + 1) % 3];
              uint64_t ko = (uint64_t)(uint32_t)(oa < ob ? oa : ob) * 2654435761ULL +
                            (uint64_t)(uint32_t)(oa < ob ? ob : oa);
              if (ko == kx) same = 1;
            }
            if (same) {
              int32_t r2a = uf_find(W->uf2, W->grp[1][o]), r2b = uf_find(W->uf2, W->grp[1][j]);
              if (r2a != r2b) W->uf2[r2a] = r2b;
            }
          }
          break;
        }
        h = (int32_t)((h + 1) & (eh - 1));
      }
    }
  }
  int32_t ncomp = 0, ncomp2 = 0;
  for (int32_t j = 0; j < k; j++) {
    if (uf_find(W->uf, j) == j) ncomp++;
    if (uf_find(W->uf2, j) == j) ncomp2++;
  }
  if (ncomp == 1) S->nsimple++;
  if (ncomp2 == 1) S->nsimple_ix++;
  return 0;
}

/* ---- ОБХОД СО СРЕЗОМ ------------------------------------------------------ */
typedef struct {
  const hz_objmesh *m;
  const hz_ptree *T;
  const double *tlo, *thi, *tpl;
  double eye[3], eps;
  double pxcut; /* 0 — срез ДО ЛИСТЬЕВ */
  cutstat *S;
  scratch *W;
  vmap *M;
  int naive;
  int fail;
} walkctx;

static int box_hits(const double *bl, const double *bh, const double *clo, const double *chi) {
  for (int c = 0; c < 3; c++)
    if (!(bl[c] < chi[c]) || !(bh[c] >= clo[c])) return 0;
  return 1;
}

/* Взят ли узел в срез. Правило §194: угловой размер `px = (2·rad/R)/ε`, спуск
 * пока `px` больше пола. Глаз ВНУТРИ узла даёт бесконечный размер, то есть
 * спуск до предела — это не особый случай, а тот же критерий. */
static int at_cut(const walkctx *X, int32_t nid) {
  const hz_ptnode *N = &X->T->nd[nid];
  if (N->child < 0) return 1;
  if (!(X->pxcut > 0.0)) return 0;
  double r2 = 0.0, d2 = 0.0;
  for (int c = 0; c < 3; c++) {
    double h = 0.5 * (N->hi[c] - N->lo[c]);
    double dd = 0.5 * (N->lo[c] + N->hi[c]) - X->eye[c];
    r2 += h * h;
    d2 += dd * dd;
  }
  double rad = sqrt(r2), R = sqrt(d2);
  if (!(R > rad)) return 0;
  return (2.0 * rad / R) / X->eps <= X->pxcut;
}

static void walk(walkctx *X, int32_t nid, int32_t *list, int32_t n) {
  if (X->fail) return;
  if (at_cut(X, nid)) {
    if (cell_do(X->m, X->tpl, X->T->nd[nid].lo, X->T->nd[nid].hi, list, n, (int)X->T->lev[nid],
                X->S, X->W, X->M, X->naive) != 0)
      X->fail = 1;
    return;
  }
  int32_t c0 = X->T->nd[nid].child;
  int32_t *sub = malloc((size_t)(n > 0 ? n : 1) * sizeof *sub);
  if (sub == NULL) {
    X->fail = 1;
    return;
  }
  for (int k = 0; k < 8; k++) {
    const double *clo = X->T->nd[c0 + k].lo, *chi = X->T->nd[c0 + k].hi;
    int32_t ns = 0;
    for (int32_t i = 0; i < n; i++) {
      const double *bl = X->tlo + 3 * (size_t)list[i], *bh = X->thi + 3 * (size_t)list[i];
      if (box_hits(bl, bh, clo, chi)) sub[ns++] = list[i];
    }
    walk(X, c0 + k, sub, ns);
  }
  free(sub);
}

/* ---- ПЕЧАТЬ --------------------------------------------------------------- */
static int64_t pct_of(const int64_t *h, int64_t tot, double q) {
  int64_t acc = 0, want = (int64_t)ceil(q * (double)tot);
  if (want < 1) want = 1;
  for (int i = 0; i <= KHIST; i++) {
    acc += h[i];
    if (acc >= want) return i;
  }
  return KHIST;
}

static void report(const cutstat *S, const hz_objmesh *m, double area_ref, double secs,
                   const vmap *M) {
  printf("   ячеек среза %lld (с геометрией %lld), уровни %d…%d, кандидатов %lld\n",
         (long long)S->ncell, (long long)S->ncell_geo, S->levmin, S->levmax, (long long)S->ncand);
  printf("   КУСКОВ %lld, ИНФЛЯЦИЯ кусков/треугольников = %.3f; кандидатов на кусок %.1f\n",
         (long long)S->npiece, (double)S->npiece / (double)m->nt,
         S->npiece > 0 ? (double)S->ncand / (double)S->npiece : 0.0);
  printf("   отсечений пустых %lld, ВЫРОЖДЕННЫХ (нулевая площадь) %lld, отступлений формулы %lld\n",
         (long long)S->nempty, (long long)S->ndegen, (long long)S->nfall);
  printf("   max вершин куска %d (обязано ≤ %d); гистограмма nv:", S->maxnv, HZ_PCLIP_MAXV);
  for (int i = 3; i <= HZ_PCLIP_MAXV; i++)
    printf(" %d:%lld", i, (long long)S->nvhist[i]);
  printf("\n");
  printf("   ПЛОЩАДЬ кусков %.6f м² против %.6f у сцены, отн. невязка %.3e\n", S->area, area_ref,
         area_ref > 0.0 ? fabs(S->area - area_ref) / area_ref : 0.0);
  if (S->ncell_geo > 0) {
    printf("   k по ПЛОСКОСТЯМ: среднее %.3f, 99-й проц. %lld, max %lld",
           (double)S->ksum / (double)S->ncell_geo, (long long)pct_of(S->khist, S->ncell_geo, 0.99),
           (long long)S->kmax);
    printf("   [допуск ×0.25: среднее %.3f max %lld; ×4: среднее %.3f max %lld]\n",
           (double)S->ksum_lo / (double)S->ncell_geo, (long long)S->kmax_lo,
           (double)S->ksum_hi / (double)S->ncell_geo, (long long)S->kmax_hi);
    printf("   k по КУСКАМ:     среднее %.3f, 99-й проц. %lld, max %lld; предел групп сработал "
           "%lld раз\n",
           (double)S->psum / (double)S->ncell_geo, (long long)pct_of(S->phist, S->ncell_geo, 0.99),
           (long long)S->pmax, (long long)S->gcap_hit);
    printf("   МНОГОПЛОСКОСТНЫХ ячеек %lld (%.2f %%), из них ОДНОСВЯЗНЫХ %lld (%.2f %%); тот же "
           "счёт по ИНДЕКСАМ вершин %lld (%.2f %%)\n",
           (long long)S->nmulti, 100.0 * (double)S->nmulti / (double)S->ncell_geo,
           (long long)S->nsimple,
           S->nmulti > 0 ? 100.0 * (double)S->nsimple / (double)S->nmulti : 0.0,
           (long long)S->nsimple_ix,
           S->nmulti > 0 ? 100.0 * (double)S->nsimple_ix / (double)S->nmulti : 0.0);
    printf("   доля многоплоскостных ПО УРОВНЯМ:");
    double sx = 0.0, sy = 0.0, sxx = 0.0, sxy = 0.0, sw = 0.0;
    for (int l = 0; l < 32; l++) {
      if (S->cell_lev[l] <= 0) continue;
      double f = (double)S->multi_lev[l] / (double)S->cell_lev[l];
      printf(" %d:%.3f", l, f);
      if (S->cell_lev[l] >= HZ_PSTOW_LEVMIN && f > 0.0) {
        double x = l, y = log2(f);
        sx += x;
        sy += y;
        sxx += x * x;
        sxy += x * y;
        sw += 1.0;
      }
    }
    if (sw >= 2.0) {
      double den = sw * sxx - sx * sx;
      printf("  -> наклон в log2 = %.3f (по %.0f уровням, где ячеек ≥ %d)",
             (den > 0.0 || den < 0.0) ? (sw * sxy - sx * sy) / den : 0.0, sw, HZ_PSTOW_LEVMIN);
    }
    printf("\n");
  }
  if (M != NULL && M->e != NULL)
    printf("   СВЕРКА ВЕРШИН: ключей %lld (каждый %d-й треугольник), общих %lld, из них с РАЗНЫХ "
           "уровней "
           "%lld;\n                  ПОБИТОВЫХ РАСХОЖДЕНИЙ %lld (на перепаде уровней %lld), не "
           "поместилось %lld\n",
           (long long)M->n, M->stride, (long long)M->nshared, (long long)M->nshared_lev,
           (long long)M->nmis, (long long)M->nmis_lev, (long long)M->nfull);
  printf("   ВРЕМЯ среза %.2f с, на кусок %.1f нс\n", secs,
         S->npiece > 0 ? 1e9 * secs / (double)S->npiece : 0.0);
}

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "pstow ФАЙЛ.obj МАСШТАБ [leaf=N] [lev=N] [grade=0|1] [eye=x,y,z]\n"
                    "      [px=8,4,2,1,0] [naive=1] [strict=1] [vtri=N]\n");
    return 1;
  }
  int leafmax = 0, maxlev = 0, grade = 1, naive = 0, strict = 0, vtri = HZ_PSTOW_VTRI;
  double eye[3] = {0.0, 0.0, 0.0};
  int haseye = 0;
  double pxl[8] = {8.0, 4.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0};
  int npx = 5;
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "leaf=", 5) == 0) leafmax = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "lev=", 4) == 0) maxlev = (int)strtol(argv[i] + 4, NULL, 10);
    if (strncmp(argv[i], "grade=", 6) == 0) grade = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "naive=", 6) == 0) naive = (int)strtol(argv[i] + 6, NULL, 10);
    if (strncmp(argv[i], "strict=", 7) == 0) strict = (int)strtol(argv[i] + 7, NULL, 10);
    if (strncmp(argv[i], "vtri=", 5) == 0) vtri = (int)strtol(argv[i] + 5, NULL, 10);
    if (strncmp(argv[i], "eye=", 4) == 0) {
      char *e = NULL;
      eye[0] = strtod(argv[i] + 4, &e);
      if (e != NULL && *e == ',') eye[1] = strtod(e + 1, &e);
      if (e != NULL && *e == ',') eye[2] = strtod(e + 1, NULL);
      haseye = 1;
    }
    if (strncmp(argv[i], "px=", 3) == 0) {
      npx = 0;
      const char *s = argv[i] + 3;
      while (*s != '\0' && npx < 8) {
        char *e = NULL;
        pxl[npx++] = strtod(s, &e);
        if (e == NULL || *e != ',') break;
        s = e + 1;
      }
    }
  }
  /* КАЖДОЕ УМОЛЧАНИЕ ПЕЧАТАЕТСЯ (§241.13): неназванное умолчание превращает
   * грубую ошибку в незаметную. Глаз берётся из `scene_cfg.h` по имени файла —
   * камеры там НАЙДЕНЫ ЗАМЕРОМ, а не назначены. */
  const char *eyesrc = "ключ eye=";
  if (!haseye) {
    /* Зал удалён 08-11; его камера теперь у комнаты-инструмента. */
    if (strstr(argv[1], "room") != NULL) {
      double e0[3] = HZ_CFG_HALL_EYE;
      memcpy(eye, e0, sizeof eye);
      eyesrc = "scene_cfg.h HALL_EYE";
    } else if (strstr(argv[1], "rungholt") != NULL) {
      double e0[3] = HZ_CFG_CITY_EYE;
      memcpy(eye, e0, sizeof eye);
      eyesrc = "scene_cfg.h CITY_EYE";
    } else {
      double e0[3] = HZ_CFG_MIGUEL_EYE;
      memcpy(eye, e0, sizeof eye);
      eyesrc = "scene_cfg.h MIGUEL_EYE";
    }
  }

  hz_objmesh m;
  double t0 = now_s();
  if (hz_obj_load(&m, argv[1], strtod(argv[2], NULL)) != 0) {
    fprintf(stderr, "нет сцены %s\n", argv[1]);
    return 1;
  }
  printf("== СЦЕНА %s: треугольников %d, вершин %d, габарит %.2f x %.2f x %.2f м, чтение %.1f с\n",
         argv[1], m.nt, m.nv, m.hi[0] - m.lo[0], m.hi[1] - m.lo[1], m.hi[2] - m.lo[2],
         now_s() - t0);
  printf("== УМОЛЧАНИЯ: leafmax %d, maxlev %d (0 — 16 и 12), градуировка %d, глаз %.2f,%.2f,%.2f "
         "(%s), ε %.4e рад/пиксель\n",
         leafmax, maxlev, grade, eye[0], eye[1], eye[2], eyesrc, HZ_CFG_EPS);
  double area_ref = 0.0;
  for (int32_t t = 0; t < m.nt; t++) {
    const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
    const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
    const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
    double u[3] = {0.0, 0.0, 0.0}, w[3] = {0.0, 0.0, 0.0}, nr[3];
    for (int c = 0; c < 3; c++) {
      u[c] = B[c] - A[c];
      w[c] = C[c] - A[c];
    }
    nr[0] = u[1] * w[2] - u[2] * w[1];
    nr[1] = u[2] * w[0] - u[0] * w[2];
    nr[2] = u[0] * w[1] - u[1] * w[0];
    area_ref += 0.5 * sqrt(nr[0] * nr[0] + nr[1] * nr[1] + nr[2] * nr[2]);
  }
  printf("== ПЛОЩАДЬ СЦЕНЫ %.6f м² (инвариант дробления, §193)\n", area_ref);

  hz_ptree T;
  t0 = now_s();
  if (strict) {
    /* НЕГАТИВНЫЙ КОНТРОЛЬ (А478): та же машинерия, СТАРАЯ укладка «влезает
     * целиком». Обязан воспроизвести опубликованные 71 / 72.5 / 79 % (§193) и
     * десятки тысяч треугольников в одном узле (39 854 на городе, §217). */
    if (hz_ptree_build_ex(&T, &m, leafmax, maxlev, 0) != 0) {
      fprintf(stderr, "отказ дерева\n");
      return 2;
    }
    double tt = now_s() - t0;
    int64_t mx = 0, mxi = 0;
    for (int32_t i = 0; i < T.nnd; i++) {
      if (T.nd[i].ntri > mx) mx = T.nd[i].ntri;
      if (T.nd[i].child >= 0 && T.nd[i].ntri > mxi) mxi = T.nd[i].ntri;
    }
    printf("== СТАРАЯ УКЛАДКА «ВЛЕЗАЕТ ЦЕЛИКОМ» за %.2f с: узлов %d, листьев %lld, ссылок %lld, "
           "избыточность %.3f\n",
           tt, T.nnd, (long long)T.nleaf, (long long)T.nref, T.redundancy);
    printf("   n_inner = %lld (%.2f %% сцены), у листьев %lld; max треугольников в узле %lld, "
           "у ВНУТРЕННЕГО %lld, в КОРНЕ %d\n",
           (long long)T.n_inner, 100.0 * (double)T.n_inner / (double)m.nt, (long long)T.n_leaf_tri,
           (long long)mx, (long long)mxi, T.nd[0].ntri);
    hz_ptree_free(&T);
    hz_obj_free(&m);
    return 0;
  }
  if (hz_ptree_build_cut(&T, &m, leafmax, maxlev, grade) != 0) {
    fprintf(stderr, "отказ дерева\n");
    return 2;
  }
  double t_tree = now_s() - t0;
  printf("== ДЕРЕВО УКЛАДКИ РЕЗКОЙ за %.2f с (спуск %.2f + градуировка %.2f, камеронезависимо): "
         "узлов %d, листьев %lld, глубина %d\n",
         t_tree, T.t_build, T.t_grade, T.nnd, (long long)T.nleaf, T.depth);
  printf("   n_inner = %lld (тождество укладки: «не влезло» невыразимо), НЕ ВЛЕЗАЕТ целиком ни в "
         "одного ребёнка %lld (%.2f %% сцены)\n",
         (long long)T.n_inner, (long long)T.n_nofit, 100.0 * (double)T.n_nofit / (double)m.nt);
  printf("   ГРАДУИРОВКА 2:1 добавила узлов %lld (%.2f %% от всех), наибольший перепад ДО неё %d\n",
         (long long)T.ngrade, T.nnd > 0 ? 100.0 * (double)T.ngrade / (double)T.nnd : 0.0, T.ngpass);
  printf("   узлов на входной треугольник %.4f; ссылок у листьев было бы %lld (избыточность "
         "%.3f)\n",
         (double)T.nnd / (double)m.nt, (long long)T.n_leaf_tri, T.redundancy);

  double *tlo = malloc(3 * (size_t)m.nt * sizeof *tlo);
  double *thi = malloc(3 * (size_t)m.nt * sizeof *thi);
  double *tpl = malloc(4 * (size_t)m.nt * sizeof *tpl);
  int32_t *list = malloc((size_t)m.nt * sizeof *list);
  if (tlo == NULL || thi == NULL || tpl == NULL || list == NULL) {
    free(tlo);
    free(thi);
    free(tpl);
    free(list);
    return 2;
  }
  for (int32_t t = 0; t < m.nt; t++) {
    const double *A = m.v + 3 * (size_t)m.f[3 * (size_t)t + 0];
    const double *B = m.v + 3 * (size_t)m.f[3 * (size_t)t + 1];
    const double *C = m.v + 3 * (size_t)m.f[3 * (size_t)t + 2];
    double u[3] = {0.0, 0.0, 0.0}, w[3] = {0.0, 0.0, 0.0}, nr[3];
    for (int c = 0; c < 3; c++) {
      u[c] = B[c] - A[c];
      w[c] = C[c] - A[c];
    }
    nr[0] = u[1] * w[2] - u[2] * w[1];
    nr[1] = u[2] * w[0] - u[0] * w[2];
    nr[2] = u[0] * w[1] - u[1] * w[0];
    double L = sqrt(nr[0] * nr[0] + nr[1] * nr[1] + nr[2] * nr[2]);
    if (!(L > 0.0)) L = 1.0;
    double *pl = tpl + 4 * (size_t)t;
    for (int c = 0; c < 3; c++)
      pl[c] = nr[c] / L;
    pl[3] = pl[0] * A[0] + pl[1] * A[1] + pl[2] * A[2];
    double *bl = tlo + 3 * (size_t)t, *bh = thi + 3 * (size_t)t;
    for (int c = 0; c < 3; c++) {
      bl[c] = 1e300;
      bh[c] = -1e300;
    }
    for (int i = 0; i < 3; i++) {
      const double *p = m.v + 3 * (size_t)m.f[3 * (size_t)t + (size_t)i];
      for (int c = 0; c < 3; c++) {
        if (p[c] < bl[c]) bl[c] = p[c];
        if (p[c] > bh[c]) bh[c] = p[c];
      }
    }
    list[t] = t;
  }

  scratch W;
  memset(&W, 0, sizeof W);
  int rc = 0;
  double infl_prev = 4.0;
  for (int q = 0; q < npx; q++) {
    cutstat S;
    memset(&S, 0, sizeof S);
    S.levmin = 99;
    S.levmax = -1;
    vmap M;
    /* ШАГ ВЫБОРКИ ДЛЯ СВЕРКИ ВЕРШИН выбирается ПОД ЁМКОСТЬ ТАБЛИЦЫ, а не
     * назначается: каждый кусок даёт до четырёх новых вершин, а кусков на
     * треугольник — измеренная инфляция предыдущего (более грубого) среза.
     * Так сверка накрывает выборку ЦЕЛИКОМ, а не обрывается на середине. */
    double est = 4.0 * 2.0 * infl_prev;
    int32_t stride =
        (int32_t)ceil((double)m.nt * est / (0.6 * (double)((int64_t)1 << HZ_PSTOW_VBITS)));
    if (stride < 1) stride = 1;
    if (vtri > 0 && m.nt / stride > vtri) stride = m.nt / vtri;
    if (vmap_init(&M, stride) != 0) {
      rc = 2;
      break;
    }
    walkctx X;
    memset(&X, 0, sizeof X);
    X.m = &m;
    X.T = &T;
    X.tlo = tlo;
    X.thi = thi;
    X.tpl = tpl;
    memcpy(X.eye, eye, sizeof eye);
    X.eps = HZ_CFG_EPS;
    X.pxcut = pxl[q];
    X.S = &S;
    X.W = &W;
    X.M = &M;
    X.naive = naive;
    t0 = now_s();
    walk(&X, 0, list, m.nt);
    double tt = now_s() - t0;
    if (X.fail) {
      fprintf(stderr, "отказ обхода\n");
      rc = 2;
      vmap_free(&M);
      break;
    }
    if (pxl[q] > 0.0)
      printf("== СРЕЗ по полу %g пикселя%s\n", pxl[q], naive ? " [НАИВНАЯ вершина]" : "");
    else
      printf("== СРЕЗ ДО ЛИСТЬЕВ%s\n", naive ? " [НАИВНАЯ вершина]" : "");
    report(&S, &m, area_ref, tt, &M);
    if (S.npiece > 0) infl_prev = (double)S.npiece / (double)m.nt;
    vmap_free(&M);
  }
  sc_free(&W);
  free(tlo);
  free(thi);
  free(tpl);
  free(list);
  hz_ptree_free(&T);
  hz_obj_free(&m);
  return rc;
}
