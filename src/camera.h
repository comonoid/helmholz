#ifndef HZ_CAMERA_H
#define HZ_CAMERA_H

/* M7: camera operator. The field is sampled on an aperture patch (a constant-z
 * plane inside the clean zone), multiplied by a thin-lens phase mask
 * exp(-i k0 rho^2 / 2f), then propagated to the sensor by the angular-spectrum
 * method (exact in homogeneous space, evanescent-safe) — two FFTs, zero extra
 * unknowns. Output: sensor intensity |u|^2. */

#include "solver3d.h"

typedef struct {
  double cx, cy; /* aperture center (cells) */
  double zap;    /* aperture plane z */
  double D;      /* aperture side */
  int n;         /* samples per side (grid step = D/n) */
  double f;      /* thin-lens focal length (cells); f<=0 disables the lens */
  double di;     /* aperture->sensor distance */
  double k0;     /* free-space wavenumber for the camera path */
} hz_camera;

/* intensity out: n*n row-major; returns 0 on ok */
int hz_camera_shoot(const hz_sol3 *sol, const hz_camera *cam, double *intensity);

/* debug/showcase: |u|^2 slice at fixed z over [0,dom]x[0,dom], w x h samples */
int hz_slice_z(const hz_sol3 *sol, double z, int w, int h, double *intensity);

#endif
