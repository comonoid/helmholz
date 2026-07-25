/* M9c: angular sparsity on a GENERATED SCENE with real multiple scattering.
 *
 * Parts 1 and 2 used superpositions of plane waves, which answered the
 * resolution question and the energy question but left the one that matters:
 * where does a REAL scene sit between "one direction" and "speckle"? PLAN
 * finding A3 says angular sparsity and speckle are one condition seen twice, so
 * a bench on a smooth scene would report sparsity and create false confidence.
 * The scene therefore has to scatter, repeatedly, and be exact.
 *
 * SCENE MODEL: N point scatterers, solved by Foldy-Lax. Each scatterer has an
 * exciting field that includes the field re-radiated by all the others:
 *     u_ex,i = u_inc(x_i) + sum_{j!=i} f G(x_i - x_j) u_ex,j
 * and the total field is u = u_inc + sum_j f G(x - x_j) u_ex,j. That is EXACT
 * for point scatterers (no Born truncation, all orders of multiple scattering),
 * costs one N x N solve, needs no grid, and is valid at any lambda — the same
 * rule every reference in this project has followed. G = (i/4) H_0^(1)(k r),
 * built from the libm Bessel functions.
 * Density is the single knob: from one scatterer (a clean directional field) up
 * to a dense cloud (a speckle field), which is exactly the axis A3 predicts the
 * compression to die along.
 *
 * MEASUREMENT: put a probe element of width W somewhere in the field, project
 * the windowed field onto M candidate carrier directions spread over 2*pi, and
 * count how many are needed for 1%. The Gram is closed-form (the carriers are
 * plane waves and the window is a tensor product of potentials); only the
 * right-hand side is sampled, because Hankel functions do not factorise. That
 * caps W at a few tens of lambda, which is where the interesting range is. */
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum { MAXSC = 400, MDIR = 256, NSPP = 8 };
static const double LAM = 16.0;
static const double RCOND = 1e-12;
static const double TARGET = 0.01;
/* Point-scatterer strength. A model parameter, not a material: it sets how much
 * each scatterer re-radiates, hence how many orders of scattering matter. */
static const double FSTRENGTH = 0.35;

static double complex hankel0(double z) {
  if (z < 1e-12) z = 1e-12;
  return CMPLX(j0(z), y0(z));
}

/* 2D free-space Green's function of (lap + k^2) */
static double complex green2d(double k, double dx, double dy) {
  double r = sqrt(dx * dx + dy * dy);
  return CMPLX(0.0, 0.25) * hankel0(k * r);
}

/* deterministic scatterer layout: a fixed low-discrepancy-ish sequence, so runs
 * are reproducible and the "density" knob is the only thing that changes */
static void place(double *xs, double *ys, int n, double box) {
  double gx = 0.7548776662466927, gy = 0.5698402909980532; /* R2 sequence */
  for (int i = 0; i < n; i++) {
    double u = fmod(0.5 + gx * (double)(i + 1), 1.0);
    double v = fmod(0.5 + gy * (double)(i + 1), 1.0);
    xs[i] = (u - 0.5) * box;
    ys[i] = (v - 0.5) * box;
  }
}

int main(void) {
  double k = 2.0 * M_PI / LAM;

  printf("M9c on a generated scene: point scatterers, exact Foldy-Lax, all\n");
  printf("orders of multiple scattering. Probe element projected onto %d\n", MDIR);
  printf("candidate directions over 2*pi; target error %.0f%%.\n\n", TARGET * 100.0);
  printf("  %8s %10s %8s %8s %8s %10s %10s\n", "scatt", "dens/lam^2", "W/lam", "slots", "N needed",
         "frac slots", "err at N");

  static double xs[MAXSC], ys[MAXSC];
  static double complex uex[MAXSC], M[MAXSC * MAXSC];
  static double complex G[MDIR * MDIR], rhs[MDIR], c[MDIR];
  static double sv[MDIR];
  static int ord[MDIR];

  static const int NSC[5] = {1, 8, 40, 120, 300};
  static const double WL[3] = {2.0, 8.0, 32.0};
  double box = 120.0 * LAM;

  for (int is = 0; is < 5; is++) {
    int ns = NSC[is];
    place(xs, ys, ns, box);
    /* Foldy-Lax: (I - f G_offdiag) u_ex = u_inc, incident wave along +x */
    for (int i = 0; i < ns; i++) {
      for (int j = 0; j < ns; j++)
        M[i * ns + j] = (i == j) ? 1.0 : -FSTRENGTH * green2d(k, xs[i] - xs[j], ys[i] - ys[j]);
      uex[i] = cexp(CMPLX(0.0, 1.0) * k * xs[i]);
    }
    lapack_int ipiv[MAXSC];
    if (LAPACKE_zgesv(LAPACK_ROW_MAJOR, ns, 1, M, ns, ipiv, uex, 1) != 0) continue;

    for (int iw = 0; iw < 3; iw++) {
      double W = WL[iw] * LAM;
      double xc = 0.0, yc = 0.0; /* probe at the centre of the cloud */
      int npt = (int)(4.0 * W / LAM * (double)NSPP);
      if (npt > 1400) npt = 1400;
      double h = 4.0 * W / (double)npt;

      /* CANDIDATE COUNT IS TIED TO WHAT THE ELEMENT CAN RESOLVE. An element of
       * width W separates directions no finer than lambda/W, so 2*pi/(lambda/W)
       * = 2*pi*W/lambda slots exist and a denser fan is pure redundancy: the
       * min-norm solve then spreads energy over near-identical directions and
       * the greedy count is inflated. The first run of this bench reported
       * N = 83 at W = 2 lambda, where only 13 slots exist at all. */
      int mdir = (int)(2.0 * M_PI * W / LAM) + 2;
      if (mdir > MDIR) mdir = MDIR;
      static double ang[MDIR];
      for (int m = 0; m < mdir; m++)
        ang[m] = 2.0 * M_PI * (double)m / (double)mdir;
      hz_phi_factor fx = {W, xc / W, 0};
      for (int m = 0; m < mdir; m++) {
        double mx = k * cos(ang[m]), my = k * sin(ang[m]);
        for (int l = 0; l < mdir; l++) {
          double lx = k * cos(ang[l]), ly = k * sin(ang[l]);
          G[m * mdir + l] = hz_phi_prod_integral_osc(xc - 2.0 * W, xc + 2.0 * W, fx, fx, lx - mx) *
                            hz_phi_prod_integral_osc(yc - 2.0 * W, yc + 2.0 * W, fx, fx, ly - my);
        }
        rhs[m] = 0.0;
      }

      /* sampled right-hand side and |windowed u|^2 */
      double nu = 0.0;
      for (int px = 0; px < npt; px++)
        for (int py = 0; py < npt; py++) {
          double x = xc - 2.0 * W + ((double)px + 0.5) * h;
          double y = yc - 2.0 * W + ((double)py + 0.5) * h;
          double wgt = hz_phi(x / W - fx.n) * hz_phi(y / W - fx.n);
          if (wgt <= 0.0) continue;
          double complex u = cexp(CMPLX(0.0, 1.0) * k * x);
          for (int j = 0; j < ns; j++)
            u += FSTRENGTH * green2d(k, x - xs[j], y - ys[j]) * uex[j];
          double complex wu = wgt * u;
          nu += creal(wu * conj(wu)) * h * h;
          for (int m = 0; m < mdir; m++) {
            double complex e = cexp(CMPLX(0.0, -1.0) * k * (cos(ang[m]) * x + sin(ang[m]) * y));
            rhs[m] += wgt * wgt * u * e * h * h;
          }
        }

      static double complex Gw[MDIR * MDIR], rw[MDIR];
      for (int i = 0; i < mdir * mdir; i++)
        Gw[i] = G[i];
      for (int i = 0; i < mdir; i++)
        rw[i] = rhs[i];
      lapack_int rank = 0;
      if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, mdir, mdir, 1, Gw, mdir, rw, 1, sv, RCOND, &rank) != 0)
        continue;
      for (int i = 0; i < mdir; i++) {
        c[i] = rw[i];
        ord[i] = i;
      }
      for (int i = 1; i < mdir; i++) {
        int key = ord[i];
        int j2 = i - 1;
        while (j2 >= 0 && cabs(c[ord[j2]]) < cabs(c[key])) {
          ord[j2 + 1] = ord[j2];
          j2--;
        }
        ord[j2 + 1] = key;
      }

      int nneed = mdir;
      double errn = -1.0;
      for (int n = 1; n <= mdir; n++) {
        static double complex Gs[MDIR * MDIR], rs[MDIR];
        for (int i = 0; i < n; i++) {
          for (int j2 = 0; j2 < n; j2++)
            Gs[i * n + j2] = G[ord[i] * mdir + ord[j2]];
          rs[i] = rhs[ord[i]];
        }
        lapack_int rk = 0;
        if (LAPACKE_zgelsd(LAPACK_ROW_MAJOR, n, n, 1, Gs, n, rs, 1, sv, RCOND, &rk) != 0) continue;
        double vv = 0.0, uv = 0.0;
        for (int i = 0; i < n; i++) {
          uv += creal(conj(rs[i]) * rhs[ord[i]]);
          for (int j2 = 0; j2 < n; j2++)
            vv += creal(rs[i] * conj(rs[j2]) * G[ord[i] * mdir + ord[j2]]);
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
      printf("  %8d %10.4f %8.0f %8d %8d %10.3f %10.3e\n", ns,
             (double)ns / ((box / LAM) * (box / LAM)), WL[iw], mdir, nneed,
             (double)nneed / (double)mdir, errn);
    }
    printf("\n");
  }
  printf("READ: N should stay small and roughly flat while the cloud is sparse,\n");
  printf("then climb towards the full candidate set as the field turns into\n");
  printf("speckle. Where it climbs is where the architecture's second factor\n");
  printf("stops being bounded — and that boundary is what A3 predicts exists.\n");
  return 0;
}
