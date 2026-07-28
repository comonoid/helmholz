/* PLAN question 3: how continuous LOD (L = eps*R) maps onto the DYADIC level
 * ladder that the two-scale relation needs — and whether that mapping can be
 * done without the degeneracy that wrecked the earlier ladder.
 *
 * THE RULE BEING TESTED. Level j has width W_j = W0 * 2^j. Shell j covers
 * distances R in [a_j, a_{j+1}) with a_j = W_j / eps; since W_{j+1} = 2 W_j,
 * a_{j+1} = 2 a_j, so the shells are DYADIC ANNULI and each doubling of
 * distance is one level coarser. Membership is by element centre, widened by
 * MARGIN elements on each side so the partition of unity (four overlapping
 * translates) holds throughout the shell rather than dipping at its edges.
 *
 * WHY THIS IS THE DECISIVE MEASUREMENT. The earlier ladder reached cond ~1e10
 * at 27 levels and its accuracy collapsed. The suspicion is that this was NOT
 * inherent to having many levels: at eps = 1/4 a shell is only 2W wide while
 * the support is 4W, so with a +/-2W margin the levels overlapped almost
 * completely, and phi is REFINABLE — a coarse function inside a region also
 * covered by a fine level is exactly reproducible by it. Redundancy at a seam
 * is unavoidable and local; redundancy everywhere is a construction error.
 * So: does conditioning grow with the NUMBER OF SHELLS, or only with their
 * OVERLAP FRACTION? The answer decides whether question 3 has a solution.
 *
 * Also measured here, because it is pure bookkeeping and needs no solver:
 * how many elements change when the camera moves — the claim the whole
 * real-time verdict rests on (PLAN A2). */
#include "carrier.h"
#include "phi.h"
#include <complex.h>
#include <lapacke.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

enum { MAXDIM = 3000, MARGIN = 2 };
static const double LAM = 16.0; /* cells per wavelength */
static const double W0 = 2.0;   /* finest element = lambda/8 */

/* One shell-quantised active set. Elements are identified by (level, node) so
 * that two sets taken at different camera positions can be compared directly. */
typedef struct {
  int lev;
  int node;
  double W;
  double kx;
} elem;

static int build_shells(elem *e, int cap, double xcam, double dom, double eps, int *nshell) {
  int d = 0;
  *nshell = 0;
  double kx = 2.0 * M_PI / LAM;
  for (int j = 0; j < 40; j++) {
    double W = W0 * pow(2.0, (double)j);
    if (W > dom) break;
    double alo = W / eps, ahi = 2.0 * alo;
    /* the finest level also owns everything nearer than its shell, the coarsest
     * everything beyond — otherwise the near field and the far field have no
     * representation at all */
    if (j == 0) alo = 0.0;
    if (W0 * pow(2.0, (double)(j + 1)) > dom) ahi = 2.0 * dom;
    alo -= (double)MARGIN * W;
    ahi += (double)MARGIN * W;
    if (alo < 0.0) alo = 0.0;
    int added = 0;
    int n0 = (int)((xcam - ahi) / W) - 2, n1 = (int)((xcam + ahi) / W) + 2;
    for (int n = n0; n <= n1; n++) {
      double xc = (double)n * W;
      if (xc < -2.0 * W || xc > dom + 2.0 * W) continue;
      double dist = fabs(xc - xcam);
      if (dist < alo || dist >= ahi) continue;
      if (d + 2 > cap) return d;
      e[d++] = (elem){j, n, W, kx};
      e[d++] = (elem){j, n, W, -kx};
      added = 1;
    }
    if (added) (*nshell)++;
  }
  return d;
}

/* Fraction of the shell that is also covered by a neighbouring shell: the
 * quantity the hypothesis says conditioning should track. Shell width is
 * a_j = W/eps; the margins add MARGIN*W on each side and the support reaches
 * 2W further, so the doubly-covered strip is ~2*(MARGIN+2)*W per seam. */
static double overlap_fraction(double eps) {
  double shell = 1.0 / eps;                    /* in units of W */
  double strip = 2.0 * ((double)MARGIN + 2.0); /* in units of W */
  return strip / shell;
}

/* L2 Gram of the directional carrier basis: <B_i,B_j> = Int phi_i phi_j
 * e^{i(kx_j - kx_i)x} dx  (CONJUGATED — this is an inner product, not the
 * bilinear Galerkin form). Returns the singular-value spread. */
static double gram_cond(const elem *e, int dim, double dom, int *rank_out) {
  double complex *G = calloc((size_t)dim * (size_t)dim, sizeof(double complex));
  double *sv = calloc((size_t)dim, sizeof(double));
  double complex *work = calloc((size_t)dim, sizeof(double complex));
  if (!G || !sv || !work) {
    free(G);
    free(sv);
    free(work);
    return -1.0;
  }
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++) {
      hz_phi_factor fi = {e[i].W, (double)e[i].node, 0};
      hz_phi_factor fj = {e[j].W, (double)e[j].node, 0};
      double lo = e[i].W * ((double)e[i].node - 2.0), hi = e[i].W * ((double)e[i].node + 2.0);
      double lo2 = e[j].W * ((double)e[j].node - 2.0), hi2 = e[j].W * ((double)e[j].node + 2.0);
      if (lo2 > lo) lo = lo2;
      if (hi2 < hi) hi = hi2;
      if (lo < 0.0) lo = 0.0;
      if (hi > dom) hi = dom;
      if (lo >= hi) continue;
      G[(size_t)i * (size_t)dim + (size_t)j] =
          hz_phi_prod_integral_osc(lo, hi, fi, fj, e[j].kx - e[i].kx);
    }
  /* singular values via zgelsd on a zero right-hand side: it factors G anyway
   * and hands back the spectrum, which is all we want here */
  lapack_int rank = 0;
  lapack_int info =
      LAPACKE_zgelsd(LAPACK_ROW_MAJOR, dim, dim, 1, G, dim, work, 1, sv, 1e-14, &rank);
  double c = (info == 0 && sv[0] > 0.0 && sv[dim - 1] > 0.0) ? sv[0] / sv[dim - 1] : -1.0;
  if (rank_out != NULL) *rank_out = (int)rank;
  free(G);
  free(sv);
  free(work);
  return c;
}

/* How many elements differ between two active sets, matched on (level, node). */
static int changed(const elem *a, int na, const elem *b, int nb) {
  int diff = 0;
  for (int i = 0; i < na; i++) {
    int found = 0;
    for (int j = 0; j < nb; j++)
      if (a[i].lev == b[j].lev && a[i].node == b[j].node) {
        found = 1;
        break;
      }
    if (!found) diff++;
  }
  for (int j = 0; j < nb; j++) {
    int found = 0;
    for (int i = 0; i < na; i++)
      if (b[j].lev == a[i].lev && b[j].node == a[i].node) {
        found = 1;
        break;
      }
    if (!found) diff++;
  }
  return diff;
}

int main(void) {
  static elem e1[MAXDIM], e2[MAXDIM];

  printf("Shell-quantised LOD: does degeneracy track SHELL COUNT or OVERLAP?\n");
  printf("shell j spans distance [W_j/eps, 2 W_j/eps), W_j = %.0f*2^j cells, margin %d elems\n\n",
         W0, MARGIN);

  printf("[1] element count vs domain — expect logarithmic growth\n");
  printf("  %12s %8s %8s\n", "domain(lam)", "shells", "dim");
  for (int i = 0; i < 6; i++) {
    double dom = 100.0 * pow(10.0, (double)i) * LAM;
    int ns = 0;
    int d = build_shells(e1, MAXDIM, 0.15 * dom, dom, 0.0625, &ns);
    printf("  %12.0f %8d %8d\n", dom / LAM, ns, d);
  }

  printf("\n[2] DECISIVE: conditioning vs number of shells, at several eps\n");
  printf("  eps gives overlap fraction = 2*(margin+2)*W / (W/eps) = %s\n", "2*(MARGIN+2)*eps");
  printf("  %8s %10s %8s %8s %10s %8s\n", "eps", "overlap", "shells", "dim", "cond", "rank/dim");
  static const double EPSV[4] = {0.25, 0.125, 0.0625, 0.03125};
  static const double DOMV[3] = {1e3, 1e5, 1e7};
  for (int a = 0; a < 4; a++) {
    for (int b = 0; b < 3; b++) {
      double dom = DOMV[b] * LAM;
      int ns = 0;
      int d = build_shells(e1, MAXDIM, 0.15 * dom, dom, EPSV[a], &ns);
      if (d >= MAXDIM) {
        printf("  %8.4f %10.3f %8d %8d   (over cap)\n", EPSV[a], overlap_fraction(EPSV[a]), ns, d);
        continue;
      }
      int rank = 0;
      double c = gram_cond(e1, d, dom, &rank);
      printf("  %8.4f %10.3f %8d %8d %10.2e %5d/%-5d\n", EPSV[a], overlap_fraction(EPSV[a]), ns, d,
             c, rank, d);
    }
  }
  printf("  READ: if cond tracks the overlap column and NOT the shell column,\n");
  printf("  question 3 has a solution — shells may be stacked as deep as needed.\n");

  printf("\n[3] camera motion: how many elements change (PLAN A2)\n");
  printf("  %10s %10s %8s %8s %10s\n", "shift/lam", "shift/R_near", "dim", "changed", "fraction");
  double dom = 1e5 * LAM, eps = 0.0625, xcam = 0.15 * dom;
  int ns = 0;
  int d1 = build_shells(e1, MAXDIM, xcam, dom, eps, &ns);
  double rnear = W0 / eps; /* inner radius of the first shell, in cells */
  for (int i = 0; i < 6; i++) {
    double shift = W0 * pow(2.0, (double)i);
    int n2 = 0;
    int d2 = build_shells(e2, MAXDIM, xcam + shift, dom, eps, &n2);
    int ch = changed(e1, d1, e2, d2);
    printf("  %10.3f %10.3f %8d %8d %10.4f\n", shift / LAM, shift / rnear, d2, ch,
           (double)ch / (double)d1);
  }
  printf("  READ: fraction should scale like shift/R_near, not like a constant.\n");
  return 0;
}
