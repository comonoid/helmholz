
#include <math.h>
#include <stdlib.h>
typedef struct { int nt; } mesh;
double now_sec(void);                 /* opaque, no args */
void hz_pyr_build(void *py, int nt, double *area, double *nrm, double *kd);
double probe_eref(const mesh *m) {
  int i, i2;
  double *Eref = (double *)malloc((size_t)m->nt * sizeof *Eref);
  double *area = (double *)malloc((size_t)m->nt * sizeof *area);
  if (!Eref || !area) return 2.0;
  for (i = 0; i < m->nt; i++) area[i] = 1.0;
  double t0 = now_sec();                 /* opaque call before the zero loop */
  for (i = 0; i < m->nt; i++) {
    area[i] = area[i] + t0;
    Eref[i] = 0.0;
  }
  hz_pyr_build((void *)m, m->nt, area, area, area); /* opaque call between */
  double d = 0.0;
  for (i2 = 0; i2 < m->nt; i2++) d += fabs(area[i2] - Eref[i2]);
  return d;
}
int main(void) { return 0; }
