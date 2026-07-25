#ifndef HZ_FFT_H
#define HZ_FFT_H

#include <complex.h>

/* In-place iterative radix-2 FFT, n = power of two.
 * sign = -1 forward, +1 inverse (inverse includes the 1/n factor). */
void hz_fft(double complex *a, int n, int sign);

#endif
