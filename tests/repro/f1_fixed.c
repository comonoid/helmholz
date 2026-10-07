/* F1 probe: sweep.c:140 "The left operand of '*' is a garbage value".
 * Shape = hz_sw_dir_table(): guarded nmu (even, >=2) -> sw_gauss_legendre()
 * fills gx[i]/gx[n-1-i] for i<n/2 -> caller reads gx[m] for m<nmu. */
#include <math.h>
#include <stdlib.h>

static void gl(int n, double *x, double *w) {
  for (int i = 0; i < n; i++) { /* VARIANT: one write per index */
    double z = 0.5;
    x[i] = (i < n / 2) ? -z : z;
    w[i] = 2.0;
  }
}

static int table(int ndirs) {
  int nphi = ndirs / 100, nmu = ndirs % 100;
  double *gx, *gw;
  if (nphi < 4 || nmu < 2 || (nmu & 1) || nphi * 100 + nmu != ndirs) return 1;
  gx = (double *)malloc((size_t)nmu * sizeof *gx);
  gw = (double *)malloc((size_t)nmu * sizeof *gw);
  if (!gx || !gw) {
    free(gx);
    free(gw);
    return 2;
  }
  gl(nmu, gx, gw);
  double acc = 0.0;
  for (int i = 0; i < nphi; i++) {
    double sp = 0.5, cp = 0.5;
    for (int m = 0; m < nmu; m++) {
      double s = sqrt(1.0 - gx[m] * gx[m]);
      acc += (sp * s + cp * s + gx[m]) * gw[m];
    }
  }
  free(gx);
  free(gw);
  return acc > 0.0;
}

int main(int argc, char **argv) {
  static int cached_ndirs = -1;
  int ndirs = 404;
  if (argc > 1) ndirs = atoi(argv[1]);
  if (ndirs != cached_ndirs) return table(ndirs);
  return 0;
}
