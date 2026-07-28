/* M9c (part 1): HOW MANY CARRIER DIRECTIONS DOES ONE ELEMENT ACTUALLY NEED?
 *
 * This is the gate the whole architecture hangs on. The count of unknowns is
 * (elements) x (directions per element); LOD bounds the first factor, and the
 * second is assumed small. "Assumed" is the operative word — it has never been
 * measured.
 *
 * WHY PROJECTION AND NOT A SOLVE. Every solver bench in this project so far was
 * spoiled by its boundary or its ladder — five separate artefacts, each of which
 * first looked like a property of the method. Best-approximation error has no
 * boundary, no cascade and no conditioning of its own, and it is a LOWER BOUND
 * on any Galerkin solve: if the directions are not there, no solver can invent
 * them. A pass here does not promise a good solve; a failure here is fatal.
 *
 * WHAT IS PREDICTED. An element of width W resolves angles no finer than
 * lambda/W (that is just diffraction from an aperture of size W). So a field
 * arriving inside an angular cone of width dtheta should need about
 *     N ~ max(1, dtheta * W / lambda)
 * carriers — the number of resolvable angular slots — and NOT the number of
 * plane waves actually present. If that holds, directional sparsity is real and
 * quantitative; if N tracks the number of incident waves instead, the element
 * is not compressing anything and the (kL)^(d-1) capacity bound is also the
 * cost.
 *
 * The target is the WINDOWED field phi(x/W)phi(y/W)*u — an element's job is to
 * carry its own share of the field under the partition of unity, not to
 * reproduce u on an unweighted square. With that target a single plane wave is
 * exactly representable, which makes the control case meaningful.
 * All fields are finite sums of plane waves, so every integral factorises into
 * two 1D closed forms and nothing is sampled. */
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum { MDIR = 96, MAXWAVE = 256 };
static const double LAM = 16.0;
static const double RCOND = 1e-12;
/* Error level at which a direction count is declared "enough". 1% is the
 * accuracy class the 1D benches reached at their best, so asking for more here
 * would measure the tail of the projection rather than the angular content. */
static const double TARGET = 0.01;

typedef struct {
  double kx, ky;
  double complex amp;
} wave;

/* Int phi(x/W)^2 e^{i w x} dx — the windowed-field weight. Both factors are the
 * same potential, so this is the ordinary two-factor oscillatory integral. */
static double complex win_osc(double W, double w) {
  hz_phi_factor f = {W, 0.0, 0};
  return hz_phi_prod_integral_osc(-2.0 * W, 2.0 * W, f, f, w);
}

int main(void) {
  double k = 2.0 * M_PI / LAM;

  printf("M9c/angle: directions needed by ONE element vs the angular width of\n");
  printf("the incident field. Prediction: N ~ max(1, dtheta*W/lambda).\n");
  printf("Target error %.0f%%, %d candidate directions spread over 2*pi.\n\n", TARGET * 100.0,
         MDIR);

  static double complex G[MDIR * MDIR], rhs[MDIR], c[MDIR];
  static double sv[MDIR];
  static wave wv[MAXWAVE];

  printf("  %10s %10s %8s %12s %8s %10s\n", "W/lambda", "dtheta", "waves", "dtheta*W/lam",
         "N needed", "err at N");
  static const double WL[4] = {1.0, 10.0, 100.0, 1000.0};
  static const double DTH[6] = {0.0, 0.001, 0.01, 0.05, 0.2, 1.0};
  for (int iw = 0; iw < 4; iw++) {
    double W = WL[iw] * LAM;
    for (int it = 0; it < 6; it++) {
      double dth = DTH[it];
      /* incident field: a fan of plane waves filling the cone. Deliberately
       * MORE waves than the predicted slot count, so that a result of
       * "N = slots" is a genuine compression and not just "N = waves". */
      int nw = 40;
      double phase0 = 0.7; /* fixed, so runs are reproducible */
      for (int p = 0; p < nw; p++) {
        double t = (nw == 1) ? 0.0 : -0.5 + (double)p / (double)(nw - 1);
        double th = t * dth;
        wv[p].kx = k * cos(th);
        wv[p].ky = k * sin(th);
        double ph = phase0 * (double)((p * 37) % 101);
        wv[p].amp = CMPLX(cos(ph), sin(ph));
      }

      /* CANDIDATE FAN. It must not be the binding constraint: a uniform sweep
       * over 2*pi puts only ONE candidate inside a narrow incident cone, and
       * then the measurement reports the poverty of the candidate set instead
       * of the angular content of the field. (That is exactly what the first
       * run of this bench did — everything past W = 10 lambda "needed" all 96.)
       * So the fan spans a few times the incident cone, and never less than a
       * few element resolutions lambda/W, sampled MDIR times inside it. */
      double res = LAM / W; /* angular resolution of this element */
      double span = 4.0 * dth;
      if (span < 8.0 * res) span = 8.0 * res;
      double dstep = span / (double)(MDIR - 1);
      static double ang[MDIR];
      for (int m = 0; m < MDIR; m++)
        ang[m] = -0.5 * span + dstep * (double)m;

      /* Gram and rhs over the candidate directions */
      for (int m = 0; m < MDIR; m++) {
        double am = ang[m];
        double mx = k * cos(am), my = k * sin(am);
        for (int l = 0; l < MDIR; l++) {
          double al = ang[l];
          double lx = k * cos(al), ly = k * sin(al);
          G[m * MDIR + l] = win_osc(W, lx - mx) * win_osc(W, ly - my);
        }
        double complex t = 0.0;
        for (int p = 0; p < nw; p++)
          t += wv[p].amp * win_osc(W, wv[p].kx - mx) * win_osc(W, wv[p].ky - my);
        rhs[m] = t;
      }
      /* ||windowed u||^2 */
      double nu = 0.0;
      for (int p = 0; p < nw; p++)
        for (int q = 0; q < nw; q++)
          nu += creal(wv[p].amp * conj(wv[q].amp) * win_osc(W, wv[p].kx - wv[q].kx) *
                      win_osc(W, wv[p].ky - wv[q].ky));

      /* full solve, then refit with the strongest n directions */
      static double complex Gw[MDIR * MDIR], rw[MDIR];
      for (int i = 0; i < MDIR * MDIR; i++)
        Gw[i] = G[i];
      for (int i = 0; i < MDIR; i++)
        rw[i] = rhs[i];
      lapack_int rank = 0;
      if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, MDIR, MDIR, 1, Gw, MDIR, rw, 1, sv, RCOND, &rank) != 0)
        continue;
      for (int i = 0; i < MDIR; i++)
        c[i] = rw[i];
      static int ord[MDIR];
      for (int i = 0; i < MDIR; i++)
        ord[i] = i;
      for (int i = 1; i < MDIR; i++) {
        int key = ord[i];
        int j = i - 1;
        while (j >= 0 && cabs(c[ord[j]]) < cabs(c[key])) {
          ord[j + 1] = ord[j];
          j--;
        }
        ord[j + 1] = key;
      }

      int nneed = MDIR;
      double errn = -1.0;
      for (int n = 1; n <= MDIR; n++) {
        static double complex Gs[MDIR * MDIR], rs[MDIR];
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++)
            Gs[i * n + j] = G[ord[i] * MDIR + ord[j]];
          rs[i] = rhs[ord[i]];
        }
        lapack_int rk = 0;
        if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, n, n, 1, Gs, n, rs, 1, sv, RCOND, &rk) != 0) continue;
        double vv = 0.0, uv = 0.0;
        for (int i = 0; i < n; i++) {
          uv += creal(conj(rs[i]) * rhs[ord[i]]);
          for (int j = 0; j < n; j++)
            vv += creal(rs[i] * conj(rs[j]) * G[ord[i] * MDIR + ord[j]]);
        }
        double e2 = nu - 2.0 * uv + vv;
        if (e2 < 0.0) e2 = 0.0;
        double e = sqrt(e2 / nu);
        if (e <= TARGET) {
          nneed = n;
          errn = e;
          break;
        }
        errn = e;
      }
      printf("  %10.0f %10.4f %8d %12.3f %8d %10.3e\n", WL[iw], dth, nw, dth * W / LAM, nneed,
             errn);
    }
    printf("\n");
  }
  printf("READ: if 'N needed' tracks the dtheta*W/lam column and NOT the 'waves'\n");
  printf("column, an element compresses its angular content and the architecture\n");
  printf("has its second factor bounded. If it tracks 'waves', it does not.\n");

  /* ------------------------------------------------------------------------
   * PART 2 — AMPLITUDE, NOT JUST WIDTH.
   * Part 1 used equal amplitudes across the cone, which is the WORST case, and
   * taken literally it breaks the architecture: a 1 cm element at 10 m lit over
   * 0.01 rad has ~180 resolvable angular slots, so ~200 carriers each, and with
   * 2e8 elements the count is gone. Real illumination is not like that: one
   * bright source plus weak diffuse light from everywhere. Since the target is
   * a RELATIVE error, directions carrying less energy than that target cannot
   * matter. This measures whether N is governed by the amplitude distribution
   * rather than by the angular width alone — if it is, the architecture's count
   * survives and "angular sparsity" means sparsity IN ENERGY, not in support. */
  printf("\n[2] one dominant direction + weak diffuse background over 1 rad\n");
  printf("  %10s %12s %12s %8s %10s\n", "W/lambda", "bg amp each", "bg energy", "N needed",
         "err at N");
  static const double BG[6] = {1.0, 0.3, 0.1, 0.03, 0.01, 0.003};
  for (int iw = 0; iw < 3; iw++) {
    double W = WL[iw + 1] * LAM; /* 10, 100, 1000 lambda */
    for (int ib = 0; ib < 6; ib++) {
      int nw = 40;
      double dth = 1.0;
      for (int p = 0; p < nw; p++) {
        double t = -0.5 + (double)p / (double)(nw - 1);
        double th = t * dth;
        wv[p].kx = k * cos(th);
        wv[p].ky = k * sin(th);
        double ph = 0.7 * (double)((p * 37) % 101);
        double a = BG[ib];
        wv[p].amp = CMPLX(a * cos(ph), a * sin(ph));
      }
      /* the dominant wave: straight ahead, unit amplitude */
      wv[nw].kx = k;
      wv[nw].ky = 0.0;
      wv[nw].amp = 1.0;
      nw++;
      double bg_energy = BG[ib] * BG[ib] * 40.0;

      /* CANDIDATES = THE INCIDENT DIRECTIONS THEMSELVES. Part 1 already measured
       * the RESOLUTION question (how finely an element can distinguish angles);
       * this phase asks the ENERGY question (how many of the directions that are
       * actually present carry enough energy to matter), and the two must not be
       * confounded. A uniform fan cannot serve here: a 1 rad cone seen by a
       * 1000-lambda element has ~1000 resolvable slots, so 96 candidates miss
       * the dominant direction entirely — which is exactly what the first run of
       * this phase reported as "N = 96 always". With candidates on the incident
       * directions the field lies exactly in their span, and N is purely a
       * statement about energy. */
      int nc = nw;
      static double angx[MDIR], angy[MDIR];
      for (int m = 0; m < nc; m++) {
        angx[m] = wv[m].kx;
        angy[m] = wv[m].ky;
      }
      for (int m = 0; m < nc; m++) {
        for (int l = 0; l < nc; l++)
          G[m * nc + l] = win_osc(W, angx[l] - angx[m]) * win_osc(W, angy[l] - angy[m]);
        double complex t = 0.0;
        for (int p = 0; p < nw; p++)
          t += wv[p].amp * win_osc(W, wv[p].kx - angx[m]) * win_osc(W, wv[p].ky - angy[m]);
        rhs[m] = t;
      }
      double nu = 0.0;
      for (int p = 0; p < nw; p++)
        for (int q = 0; q < nw; q++)
          nu += creal(wv[p].amp * conj(wv[q].amp) * win_osc(W, wv[p].kx - wv[q].kx) *
                      win_osc(W, wv[p].ky - wv[q].ky));

      static double complex Gw2[MDIR * MDIR], rw2[MDIR];
      for (int i = 0; i < nc * nc; i++)
        Gw2[i] = G[i];
      for (int i = 0; i < nc; i++)
        rw2[i] = rhs[i];
      lapack_int rk2 = 0;
      if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, nc, nc, 1, Gw2, nc, rw2, 1, sv, RCOND, &rk2) != 0)
        continue;
      for (int i = 0; i < nc; i++)
        c[i] = rw2[i];
      static int ord2[MDIR];
      for (int i = 0; i < nc; i++)
        ord2[i] = i;
      for (int i = 1; i < nc; i++) {
        int key = ord2[i];
        int j = i - 1;
        while (j >= 0 && cabs(c[ord2[j]]) < cabs(c[key])) {
          ord2[j + 1] = ord2[j];
          j--;
        }
        ord2[j + 1] = key;
      }
      int nneed = nc;
      double errn = -1.0;
      for (int n = 1; n <= nc; n++) {
        static double complex Gs[MDIR * MDIR], rs[MDIR];
        for (int i = 0; i < n; i++) {
          for (int j = 0; j < n; j++)
            Gs[i * n + j] = G[ord2[i] * nc + ord2[j]];
          rs[i] = rhs[ord2[i]];
        }
        lapack_int rk = 0;
        if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, n, n, 1, Gs, n, rs, 1, sv, RCOND, &rk) != 0) continue;
        double vv = 0.0, uv = 0.0;
        for (int i = 0; i < n; i++) {
          uv += creal(conj(rs[i]) * rhs[ord2[i]]);
          for (int j = 0; j < n; j++)
            vv += creal(rs[i] * conj(rs[j]) * G[ord2[i] * nc + ord2[j]]);
        }
        double e2 = nu - 2.0 * uv + vv;
        if (e2 < 0.0) e2 = 0.0;
        double e = sqrt(e2 / nu);
        if (e <= TARGET) {
          nneed = n;
          errn = e;
          break;
        }
        errn = e;
      }
      printf("  %10.0f %12.4f %12.4f %8d %10.3e\n", WL[iw + 1], BG[ib], bg_energy, nneed, errn);
    }
    printf("\n");
  }
  printf("READ: if N collapses as the background energy drops below the target,\n");
  printf("then angular sparsity is sparsity IN ENERGY and the count survives.\n");
  return 0;
}
