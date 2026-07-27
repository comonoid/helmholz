/* PLAN_TRANSPORT.md T5а. Камера-обскура, прямой свет, сбор по пикселю. */

#include "transport/cam3.h"
#include <math.h>
#include <string.h>

static void cross3(const double a[3], const double b[3], double o[3]) {
  o[0] = a[1] * b[2] - a[2] * b[1];
  o[1] = a[2] * b[0] - a[0] * b[2];
  o[2] = a[0] * b[1] - a[1] * b[0];
}

static double norm3(double v[3]) {
  double m = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
  if (m > 0.0)
    for (int k = 0; k < 3; k++)
      v[k] /= m;
  return m;
}

int tr3_camera_look(tr3_camera *c, const double eye[3], const double at[3], const double up[3],
                    double fov, int w, int h) {
  if (w <= 0 || h <= 0 || !(fov > 0.0) || !(fov < M_PI)) return 1;
  memset(c, 0, sizeof *c);
  for (int k = 0; k < 3; k++) {
    c->eye[k] = eye[k];
    c->fwd[k] = at[k] - eye[k];
  }
  if (!(norm3(c->fwd) > 0.0)) return 1;
  cross3(c->fwd, up, c->right);
  if (!(norm3(c->right) > 0.0)) return 1; /* up параллелен взгляду */
  cross3(c->right, c->fwd, c->up);
  norm3(c->up);
  c->tanx = tan(0.5 * fov);
  /* пиксель КВАДРАТНЫЙ: вертикаль выводится из соотношения сторон */
  c->tany = c->tanx * (double)h / (double)w;
  c->w = w;
  c->h = h;
  return 0;
}

void tr3_camera_ray(const tr3_camera *c, int px, int py, double o[3], double d[3]) {
  /* ЦЕНТР пикселя (К20): +0.5, а не угол. Ось y экрана вниз. */
  double sx = (2.0 * ((double)px + 0.5) / (double)c->w - 1.0) * c->tanx;
  double sy = (1.0 - 2.0 * ((double)py + 0.5) / (double)c->h) * c->tany;
  for (int k = 0; k < 3; k++) {
    o[k] = c->eye[k];
    d[k] = c->fwd[k] + sx * c->right[k] + sy * c->up[k];
  }
  norm3(d);
}

double tr3_shade(const tr3_scene *sc, const tr3_sun *sun, const tr3_material *mat,
                 const double o[3], const double d[3], int shadow) {
  tr3_hit h;
  tr3_march(sc, o, d, -1.0, &h);
  if (!h.hit) return 0.0; /* фона нет: T5а не рисует небо */

  double cosi = 0.0;
  for (int k = 0; k < 3; k++)
    cosi += h.n[k] * sun->dir[k];
  if (!(cosi > 0.0)) return 0.0; /* обратная сторона: сама себя и затеняет */

  double vis = 1.0;
  if (shadow) {
    /* Теневой луч выходит РОВНО с поверхности и не смещается по нормали: вход в
     * полупространство требует ds/dt < 0 строго, а при взгляде наружу
     * производная положительна (ray3.h). ε не нужен. */
    tr3_hit s;
    tr3_march(sc, h.p, sun->dir, -1.0, &s);
    /* НАКОПЛЕННАЯ толщина, а не булева видимость (К16/С4) */
    vis = s.hit ? 0.0 : exp(-s.tau);
  }
  double lout = mat->albedo / M_PI * sun->e * cosi * vis;
  return exp(-h.tau) * lout;
}

int tr3_render(const tr3_scene *sc, const tr3_camera *cam, const tr3_sun *sun,
               const tr3_material *mat, int shadow, double *buf) {
  if (buf == NULL) return 1;
  for (int py = 0; py < cam->h; py++)
    for (int px = 0; px < cam->w; px++) {
      double o[3], d[3];
      tr3_camera_ray(cam, px, py, o, d);
      buf[(size_t)py * (size_t)cam->w + (size_t)px] = tr3_shade(sc, sun, mat, o, d, shadow);
    }
  return 0;
}
