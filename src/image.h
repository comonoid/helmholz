#ifndef HZ_IMAGE_H
#define HZ_IMAGE_H

/* Minimal image output: binary PPM (P6), zero dependencies. Values are
 * intensities >= 0; tone mapping = normalize to max, gamma 0.5, inferno-ish
 * tiny colormap. */

int hz_ppm_write(const char *path, const double *intensity, int w, int h);

/* signed field (e.g. Re u): blue-white-red, symmetric percentile scaling */
int hz_ppm_write_signed(const char *path, const double *field, int w, int h);

/* raw rgb writer for composed images */
int hz_ppm_write_rgb(const char *path, const unsigned char *rgb, int w, int h);

#endif
