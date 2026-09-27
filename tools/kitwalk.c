/* kitwalk.c — §914-Ш2': ДВИЖЕНИЕ ФРОНТА ПО ГЕОМЕТРИЧЕСКОЙ СТРУКТУРЕ.
 * Синтетический пол 16×16 кводов, три уровня лестницы (L0=512 тр,
 * L1=32, L2=2), одна пирамида (занятость общая, ℓ_p тр = уровень).
 *
 * Кусок с этажом P (сторона ~cell·2^P; P=0 — листовой) идёт сверху:
 *   ПУСТ узел            → прыжок (jumps), стоимость 0;
 *   узел уровня l, l+1==P → СОБЫТИЕ: кусок дошёл до своей детальности;
 *   иначе                → спуск к детям (порядок детей — Ш3, формулой).
 * Клэмп: P больше высоты пирамиды → обслуживание в корне.
 */
#include "nstruct/pyr.h"

#include <stdio.h>
#include <string.h>

#define NQ 16
#define NLEV 3
#define NEV 8 /* этажей событий: лист + 7 уровней вверх */

typedef struct {
  int64_t visits;  /* узлов+листьев посещено (без прыжков) */
  int64_t jumps;   /* пустых детей пропущено */
  int64_t ev[NEV]; /* ev[0]=лист, ev[l+1]=узлы уровня l */
  int64_t ev_total;
} walk_stat;

static void serve(walk_stat *st, int slot) {
  if (slot >= 0 && slot < NEV) st->ev[slot]++;
  st->ev_total++;
}

static void walk(const hz_pyr *py, int32_t l, int64_t id, int P, const double om[3],
                 walk_stat *st) {
  (void)om; /* порядок детей формулой направления — Ш3; множество то же */
  st->visits++;
  if (l + 1 == P) { /* свой этаж — событие, спуск прекращён */
    serve(st, l + 1);
    return;
  }
  if (l == 0) { /* дети — листья; сюда попадает только P=0 (P=1 обслужен выше).
                    id — узел сетки уровня 0: разлагаем ЕЁ размерами, дети —
                    в листовой сетке (nx, ny, nz) */
    int64_t dl[3];
    hz_pyr_level_dims(py, 0, dl);
    int64_t ci[3];
    ci[0] = id % dl[0];
    ci[1] = (id / dl[0]) % dl[1];
    ci[2] = id / (dl[0] * dl[1]);
    for (int cz = 0; cz < 2; cz++)
      for (int cy = 0; cy < 2; cy++)
        for (int cx = 0; cx < 2; cx++) {
          int64_t ccx = ci[0] * 2 + cx, ccy = ci[1] * 2 + cy, ccz = ci[2] * 2 + cz;
          if (ccx >= py->nx || ccy >= py->ny || ccz >= py->nz) {
            st->jumps++; /* за пределами листовой сетки: без проверки id
                            АЛИАСИТСЯ в соседний ряд (урок этого прогона) */
            continue;
          }
          int64_t cid = ccx + (ccy + ccz * py->ny) * py->nx;
          if (hz_pyr_leaf_pos(py, cid) < 0) {
            st->jumps++;
            continue;
          }
          st->visits++;
          serve(st, 0);
        }
    return;
  }
  int64_t dl[3], cd[3];
  hz_pyr_level_dims(py, l, dl);
  hz_pyr_level_dims(py, l - 1, cd);
  int64_t ci[3];
  ci[0] = id % dl[0];
  ci[1] = (id / dl[0]) % dl[1];
  ci[2] = id / (dl[0] * dl[1]);
  if (ci[2] >= dl[2] || ci[1] >= dl[1] || ci[0] >= dl[0]) {
    printf("BADID l=%d id=%lld dims=%lldx%lldx%lld\n", l, (long long)id, (long long)dl[0],
           (long long)dl[1], (long long)dl[2]);
    return;
  }
  for (int cz = 0; cz < 2; cz++)
    for (int cy = 0; cy < 2; cy++)
      for (int cx = 0; cx < 2; cx++) {
        int64_t ccx = ci[0] * 2 + cx, ccy = ci[1] * 2 + cy, ccz = ci[2] * 2 + cz;
        if (ccx >= cd[0] || ccy >= cd[1] || ccz >= cd[2]) {
          st->jumps++; /* срез за пределами размерностей уровня НЕ существует —
                          без проверки его линейный id АЛИАСИТСЯ в соседний ряд
                          (найдено дублями листьев в этом же прогоне) */
          continue;
        }
        int64_t cid = ccx + (ccy + ccz * cd[1]) * cd[0];
        if (hz_pyr_node_pos(py, l - 1, cid) < 0) {
          st->jumps++;
          continue;
        }
        walk(py, l - 1, cid, P, om, st);
      }
}

/* Три-массивы синтетики — в файловой области: их видит и верификатор
 * pyr (детекторы (а)(б) требуют пересчёта центроидов). */
enum { NTMAX = 2 * (16 * 16 + 4 * 4 + 1) };
static double tmin[3 * NTMAX], tmax[3 * NTMAX], cen[3 * NTMAX];
static int32_t mtl[NTMAX];
static uint8_t lp[NTMAX];
static int nt_global = 0;

static int build_floor(hz_pyr *py) {
  int nt = 0;
  for (int li = 0; li < NLEV; li++) {
    int nq = NQ >> (2 * li); /* сторона квода: 16 → 4 → 1 (площадь ×16) */
    double q = 16.0 / nq;
    for (int j = 0; j < nq; j++)
      for (int i = 0; i < nq; i++) {
        double x0 = i * q, y0 = j * q;
        for (int h = 0; h < 2; h++) {
          double *mn = &tmin[3 * nt], *mx = &tmax[3 * nt], *cn = &cen[3 * nt];
          mn[0] = x0;
          mn[1] = y0;
          mn[2] = -1e-3;
          mx[0] = x0 + q;
          mx[1] = y0 + q;
          mx[2] = 1e-3;
          cn[0] = x0 + q / 2;
          cn[1] = y0 + q / 2;
          cn[2] = 0;
          mtl[nt] = 0;
          lp[nt] = (uint8_t)li;
          nt++;
        }
      }
  }
  double lo[3] = {0, 0, -0.5}, hi[3] = {16, 16, 0.5};
  if (hz_pyr_build(py, nt, tmin, tmax, cen, mtl, lo, hi, 1.0)) return 1;
  int32_t nup = 0;
  if (hz_pyr_set_lp(py, lp, &nup)) return 1;
  nt_global = nt;
  printf("floor: nt=%d (L0=%d L1=%d L2=%d) grid=%dx%dx%d nleaf=%d ncentleaf=%d nup=%d nlev=%d\n",
         nt, 2 * 16 * 16, 2 * 4 * 4, 2, py->nx, py->ny, py->nz, py->nleaf, (int)py->ncentleaf, nup,
         py->nlev);
  return 0;
}

static void run_walk(const hz_pyr *py, int32_t ltop, int P, const double om[3], const char *tag) {
  walk_stat st;
  memset(&st, 0, sizeof st);
  if (P > ltop + 1) { /* клэмп: кусок крупнее корня — сервис в корне */
    st.visits = 1;
    serve(&st, ltop + 1);
  } else {
    walk(py, ltop, 0, P, om, &st);
  }
  printf("%s P=%d: visits=%lld jumps=%lld events=%lld (leaf=%lld", tag, P, (long long)st.visits,
         (long long)st.jumps, (long long)st.ev_total, (long long)st.ev[0]);
  for (int l = 1; l < NEV; l++)
    if (st.ev[l]) printf(" u%d=%lld", l - 1, (long long)st.ev[l]);
  printf(")\n");
}

int main(void) {
  hz_pyr py;
  memset(&py, 0, sizeof py);
  if (build_floor(&py)) {
    printf("build FAIL\n");
    return 1;
  }
  enum { NT = 2 * (16 * 16 + 4 * 4 + 1) };
  (void)NT;
  hz_pyr_verdict vd;
  hz_pyr_verify(&py, nt_global, tmin, tmax, cen, &vd);
  printf("verify: d_cell=%lld d_csr=%lld d_bbox=%lld d_empty=%lld\n", (long long)vd.d_cell,
         (long long)vd.d_csr, (long long)vd.d_bbox, (long long)vd.d_empty);

  int32_t ltop = py.nlev - 1;
  int64_t dtop[3];
  hz_pyr_level_dims(&py, ltop, dtop);
  printf("top: level=%d dims=%lldx%lldx%lld\n", ltop, (long long)dtop[0], (long long)dtop[1],
         (long long)dtop[2]);

  const double om[3] = {0.48, 0.6, 0.64};
  run_walk(&py, ltop, 0, om, "walk");
  run_walk(&py, ltop, 2, om, "walk");
  /* П5: ω → −ω: события неизменны */
  const double om2[3] = {-0.48, -0.6, -0.64};
  walk_stat s1, s2;
  memset(&s1, 0, sizeof s1);
  memset(&s2, 0, sizeof s2);
  walk(&py, ltop, 0, 0, om, &s1);
  walk(&py, ltop, 0, 0, om2, &s2);
  printf("NK omega-flip: events %lld vs %lld %s\n", (long long)s1.ev_total, (long long)s2.ev_total,
         s1.ev_total == s2.ev_total ? "OK" : "FAIL");
  hz_pyr_free(&py);
  return 0;
}
