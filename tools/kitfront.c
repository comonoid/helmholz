/* kitfront.c — §915: многоуровневый свип-прототип.
 *
 * Кусок = (этаж P, энергия E, аккумулятор d). Ходит по ПИРАМИДЕ
 * СВОЕГО уровня (Ш3): ПУСТ → прыжок, занятый лист → СОБЫТИЕ:
 * поглощение E·(1−ρ), отражение E·ρ в повёрнутое направление,
 * d += 0.5 (полностью диффузный отскок; d ≥ 1 → этаж грубее).
 * Хопы = итерации свипа: куски перевыпускаются от клеток-приёмников.
 *
 * Закон сохранения — двойная бухгалтерия: Σиспущено = Σпоглощено +
 * E_живых на обрыве; невязка обязана быть ≤ 1e-12 относительная.
 */
#include "nstruct/pyr.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NQ 16
#define NLEV 2
#define NHOPS 12
#define RHO 0.5
#define ECUT 1e-4 /* обрыв куска: энергия ниже — поглощается целиком */
#define NZ 2      /* ожидаемая высота пирамид (сетка 17×17×2 / 5×5×1) */
#define MAXPIECES (1 << 14)

typedef struct {
  double x, y, z; /* позиция (центр клетки-приёмника последнего события) */
  double om[3];
  double E;
  int P;
  double d;
} piece;

typedef struct {
  int64_t visits;
  int64_t events_l0;
  int64_t events_l1;
} hop_stat;

static int build_pyr(hz_pyr *py, int nq, double cell, uint8_t lpv) {
  enum { NTMAX = 2 * 256 };
  static double tmin[3 * NTMAX], tmax[3 * NTMAX], cen[3 * NTMAX];
  static int32_t mtl[NTMAX];
  static uint8_t lp[NTMAX];
  int nt = 0;
  double q = 16.0 / nq;
  for (int j = 0; j < nq; j++)
    for (int i = 0; i < nq; i++)
      for (int h = 0; h < 2; h++) {
        double *mn = &tmin[3 * nt], *mx = &tmax[3 * nt], *cn = &cen[3 * nt];
        mn[0] = i * q;
        mn[1] = j * q;
        mn[2] = -1e-3;
        mx[0] = i * q + q;
        mx[1] = j * q + q;
        mx[2] = 1e-3;
        cn[0] = i * q + q / 2;
        cn[1] = j * q + q / 2;
        cn[2] = 0;
        mtl[nt] = 0;
        lp[nt] = lpv;
        nt++;
      }
  double lo[3] = {0, 0, -0.5}, hi[3] = {16, 16, 0.5};
  int32_t nup = 0;
  if (hz_pyr_build(py, nt, tmin, tmax, cen, mtl, lo, hi, cell)) return 1;
  if (hz_pyr_set_lp(py, lp, &nup)) return 1;
  return 0;
}

/* обход пирамиды сверху до ПЕРВОГО занятого листа по направлению om
 * (v1: ближайший по ходу лист без точного DDA — прототип энергии);
 * возвращает лист или -1 */
static int64_t first_hit(const hz_pyr *py, const double om[3], const double p[3], int64_t *visits) {
  /* v1: последовательный марш по листовой сетке вдоль om от p —
   * простая линейная трассировка по cells (прототип, не DDA) */
  int64_t i = (int64_t)floor(p[0] / py->cell);
  int64_t j = (int64_t)floor(p[1] / py->cell);
  double side = py->cell;
  if (fabs(om[0]) < 1e-12 && fabs(om[1]) < 1e-12) return -1; /* вертикальных нет */
  for (int step = 0; step < 3 * (py->nx + py->ny); step++) {
    if (i < 0 || j < 0 || i >= py->nx || j >= py->ny) return -1;
    (*visits)++;
    int64_t id = i + py->nx * j;
    if (hz_pyr_leaf_pos(py, id) >= 0) {
      int32_t pos = hz_pyr_leaf_pos(py, id);
      if (py->leaf[pos].npcs > 0) return id; /* материален */
    }
    /* шаг к следующей клетке: пересечение луча с границами клетки */
    double ci = (double)i, cj = (double)j;
    double tx = om[0] > 0 ? ((ci + 1.0) * side - p[0]) / om[0]
                          : (om[0] < 0 ? (ci * side - p[0]) / om[0] : 1e30);
    double ty = om[1] > 0 ? ((cj + 1.0) * side - p[1]) / om[1]
                          : (om[1] < 0 ? (cj * side - p[1]) / om[1] : 1e30);
    double t = tx < ty ? tx : ty;
    if (t >= 1e29) return -1;
    if (tx < ty)
      i += om[0] > 0 ? 1 : -1;
    else
      j += om[1] > 0 ? 1 : -1;
  }
  return -1;
}

int main(void) {
  hz_pyr plev[NLEV];
  memset(plev, 0, sizeof plev);
  if (build_pyr(&plev[0], NQ, 1.0, 0) || build_pyr(&plev[1], 4, 4.0, 0)) {
    printf("build FAIL\n");
    return 1;
  }
  printf("pyr0: %dx%d nleaf=%d; pyr1: %dx%d nleaf=%d\n", plev[0].nx, plev[0].ny, plev[0].nleaf,
         plev[1].nx, plev[1].ny, plev[1].nleaf);

  /* источник: лампа-квадрат в центре, 512 кусков по кругу направлений */
  piece pcs[MAXPIECES];
  int np = 0;
  double emitted = 0;
  for (int a = 0; a < 512; a++) {
    double ang = 2.0 * M_PI * a / 512.0;
    pcs[np].x = 8.0;
    pcs[np].y = 8.0;
    pcs[np].z = 0;
    pcs[np].om[0] = cos(ang);
    pcs[np].om[1] = sin(ang);
    pcs[np].om[2] = 0;
    pcs[np].E = 1.0 / 512.0;
    pcs[np].P = 0;
    pcs[np].d = 0;
    emitted += pcs[np].E;
    np++;
  }
  double absorbed = 0;
  double alive = 0;

  for (int hop = 1; hop <= NHOPS; hop++) {
    hop_stat st;
    memset(&st, 0, sizeof st);
    double absorbed_hop = 0;
    int born = 0;
    piece next[MAXPIECES]; /* cppcheck не видит запись через q-копию — явный нуль */
    memset(next, 0, sizeof next);
    int nn = 0;
    for (int u = 0; u < np; u++) {
      const hz_pyr *py = &plev[pcs[u].P];
      int64_t visits = 0;
      double pos[3] = {pcs[u].x, pcs[u].y, pcs[u].z};
      int64_t leaf = first_hit(py, pcs[u].om, pos, &visits);
      st.visits += visits;
      if (leaf < 0) {
        /* улетел: энергия в живых (обрыв прототипа — границы) */
        alive += pcs[u].E;
        continue;
      }
      if (pcs[u].P == 0)
        st.events_l0++;
      else
        st.events_l1++;
      double eabs = pcs[u].E * (1.0 - RHO);
      absorbed_hop += eabs;
      double enew = pcs[u].E - eabs;
      if (enew < ECUT) {
        absorbed_hop += enew;
        continue;
      }
      /* отражение: поворот на 137.5° (золотой угол, детерминированно) */
      piece q = pcs[u];
      double c = cos(2.399963), s = sin(2.399963);
      q.om[0] = pcs[u].om[0] * c - pcs[u].om[1] * s;
      q.om[1] = pcs[u].om[0] * s + pcs[u].om[1] * c;
      q.E = enew;
      q.x = (floor(pcs[u].x / py->cell) + 0.5) * py->cell;
      q.y = (floor(pcs[u].y / py->cell) + 0.5) * py->cell;
      q.d += 0.5;
      if (q.d >= 1.0 && q.P < NLEV - 1) {
        q.d -= 1.0;
        q.P++;
      }
      if (nn < MAXPIECES) next[nn++] = q;
      born++;
    }
    absorbed += absorbed_hop;
    memcpy(pcs, next, (size_t)nn * sizeof(piece));
    np = nn;
    printf("хоп%2d: кусков=%d визитов=%lld событий L0=%lld L1=%lld поглощено=%.5f\n", hop, np,
           (long long)st.visits, (long long)st.events_l0, (long long)st.events_l1, absorbed_hop);
    if (np == 0) break;
  }
  alive += np > 0 ? 0 : 0;
  for (int u = 0; u < np; u++)
    alive += pcs[u].E; /* хвост живых */
  double bal = emitted - absorbed - alive;
  printf("БАЛАНС: испущено=%.6f поглощено=%.6f живые=%.6f невязка=%.3e (отн %.3e) %s\n", emitted,
         absorbed, alive, bal, fabs(bal) / emitted, fabs(bal) / emitted <= 1e-12 ? "OK" : "FAIL");
  hz_pyr_free(&plev[0]);
  hz_pyr_free(&plev[1]);
  return 0;
}
