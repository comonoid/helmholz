#include "fft.h"
#include <math.h>

void hz_fft(double complex *a, int n, int sign) {
  for (int i = 1, j = 0; i < n; i++) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1)
      j ^= bit;
    j |= bit;
    if (i < j) {
      double complex t = a[i];
      a[i] = a[j];
      a[j] = t;
    }
  }
  for (int len = 2; len <= n; len <<= 1) {
    double ang = (double)sign * 2.0 * M_PI / (double)len;
    double complex wl = CMPLX(cos(ang), sin(ang));
    for (int i = 0; i < n; i += len) {
      double complex w = 1.0;
      for (int k = 0; k < len / 2; k++) {
        double complex u = a[i + k];
        double complex v = a[i + k + len / 2] * w;
        a[i + k] = u + v;
        a[i + k + len / 2] = u - v;
        w *= wl;
      }
    }
  }
  if (sign > 0)
    for (int i = 0; i < n; i++)
      a[i] /= (double)n;
}
