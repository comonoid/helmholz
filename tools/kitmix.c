/* kitmix.c — §914-П3 / §915-R2: смешанный кит по ЗОНЕ.
 *
 * Два режима:
 *   split=S     — плоскость z: L0 ниже, L1 выше (П3-эксперимент);
 *   eye=,zone=R — КОЛЬЦА вокруг глаза: треугольник уровня i берётся,
 *                 если ring = min(floor(d/R), nlev-1) == i, d —
 *                 расстояние центроида до глаза (камерная зона §913).
 * Выход — одноуровневый кит для обычного pgather (свип не меняется).
 *
 * Синтаксис: kitmix IN.kit OUT.kit split=S | eye=X,Y,Z zone=R
 */
#include "geom/kit.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
  double split = 0.0;
  double eye[3] = {0, 0, 0}, zone = -1.0; /* zone<0 — режим split */
  if (argc < 3) {
    fprintf(stderr, "kitmix IN.kit OUT.kit [split=S | eye=X,Y,Z zone=R]\n");
    return 2;
  }
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "split=", 6) == 0)
      split = atof(argv[i] + 6);
    else if (strncmp(argv[i], "zone=", 5) == 0)
      zone = atof(argv[i] + 5);
    else if (strncmp(argv[i], "eye=", 4) == 0) {
      if (sscanf(argv[i] + 4, "%lf,%lf,%lf", &eye[0], &eye[1], &eye[2]) != 3) {
        fprintf(stderr, "kitmix: eye=X,Y,Z\n");
        return 2;
      }
    }
  }
  int radial = zone > 0;

  hz_kit k;
  hz_kit_init(&k);
  FILE *f = fopen(argv[1], "rb");
  if (f == NULL) {
    fprintf(stderr, "kitmix: не читается %s\n", argv[1]);
    return 2;
  }
  int rc = hz_kit_load(&k, f);
  fclose(f);
  if (rc != HZ_KIT_OK || k.nlev < 1) {
    fprintf(stderr, "kitmix: кит не читается (rc=%d nlev=%d)\n", rc, k.nlev);
    hz_kit_free(&k);
    return 2;
  }
  uint32_t vbase[64];
  uint32_t nv = 0;
  for (int32_t li = 0; li < k.nlev; li++) {
    vbase[li] = nv;
    nv += k.lev[li].nverts;
  }

  /* Счёт: ring(i) == уровень i (радиальный) либо плоскость split. */
  uint32_t nt = 0;
  uint32_t took[64] = {0};
  if (k.nlev > 64) {
    fprintf(stderr, "kitmix: nlev>64\n");
    hz_kit_free(&k);
    return 2;
  }
  for (int32_t li = 0; li < k.nlev; li++) {
    const hz_kit_level *S = &k.lev[li];
    for (uint32_t t = 0; t < S->ntris; t++) {
      double cx = (S->vx[S->ti0[t]] + S->vx[S->ti1[t]] + S->vx[S->ti2[t]]) / 3.0;
      double cy = (S->vy[S->ti0[t]] + S->vy[S->ti1[t]] + S->vy[S->ti2[t]]) / 3.0;
      double cz = (S->vz[S->ti0[t]] + S->vz[S->ti1[t]] + S->vz[S->ti2[t]]) / 3.0;
      int take;
      if (radial) {
        double d = sqrt((cx - eye[0]) * (cx - eye[0]) + (cy - eye[1]) * (cy - eye[1]) +
                        (cz - eye[2]) * (cz - eye[2]));
        int32_t ring = (int32_t)(d / zone);
        if (ring > k.nlev - 1) ring = k.nlev - 1;
        take = (ring == li);
      } else {
        take = (li == 0) ? (cz < split) : (cz >= split);
      }
      if (take) {
        took[li]++;
        nt++;
      }
    }
  }
  printf("смешивание (nlev=%d, %s): ", k.nlev, radial ? "кольца" : "плоскость");
  for (int32_t li = 0; li < k.nlev; li++)
    printf("L%d=%u ", li, took[li]);
  printf("итого=%u\n", nt);
  if (nt == 0) {
    fprintf(stderr, "kitmix: пустой результат\n");
    hz_kit_free(&k);
    return 2;
  }

  /* Выходной кит: вершины ВСЕХ уровней (смещения vbase), треугольники
   * подряд по уровням, грозди-прогоны по источнику. */
  hz_kit out;
  hz_kit_init(&out);
  out.nlev = 1;
  out.flags = HZ_KIT_FLAG_F64;
  out.lev = calloc(1, sizeof(hz_kit_level));
  out.nmtl = k.nmtl;
  out.mtl = malloc((size_t)k.nmtl * sizeof(hz_kit_mtl));
  if (out.lev == NULL || out.mtl == NULL) {
    fprintf(stderr, "kitmix: нет памяти\n");
    hz_kit_free(&out);
    hz_kit_free(&k);
    return 2;
  }
  memcpy(out.mtl, k.mtl, (size_t)k.nmtl * sizeof(hz_kit_mtl));
  hz_kit_level *M = &out.lev[0];
  M->nverts = nv;
  M->ntris = nt;
  M->vx = malloc((size_t)nv * sizeof(double));
  M->vy = malloc((size_t)nv * sizeof(double));
  M->vz = malloc((size_t)nv * sizeof(double));
  M->ti0 = malloc(4ull * nt);
  M->ti1 = malloc(4ull * nt);
  M->ti2 = malloc(4ull * nt);
  M->tcl = malloc(4ull * nt);
  M->tmtl = malloc(4ull * nt);
  if (!M->vx || !M->vy || !M->vz || !M->ti0 || !M->ti1 || !M->ti2 || !M->tcl || !M->tmtl) {
    fprintf(stderr, "kitmix: нет памяти (уровень)\n");
    hz_kit_free(&out);
    hz_kit_free(&k);
    return 2;
  }
  for (int32_t li = 0; li < k.nlev; li++) {
    const hz_kit_level *S = &k.lev[li];
    memcpy(M->vx + vbase[li], S->vx, (size_t)S->nverts * sizeof(double));
    memcpy(M->vy + vbase[li], S->vy, (size_t)S->nverts * sizeof(double));
    memcpy(M->vz + vbase[li], S->vz, (size_t)S->nverts * sizeof(double));
  }

  /* Грозди заполняются В ЭТОМ ЖЕ цикле: массив размером nt (гроздь ≥1
   * треугольника ⇒ ncl ≤ nt), после — усечение realloc. Индукция по ci
   * в отдельном цикле была невидима анализатору (CWE-122 FP-класс). */
  M->cl = calloc(nt, sizeof(hz_cluster));
  if (M->cl == NULL) {
    fprintf(stderr, "kitmix: нет памяти (грозди)\n");
    hz_kit_free(&out);
    hz_kit_free(&k);
    return 2;
  }
  uint32_t w = 0;
  uint32_t ncl = 0;
  int prev_src = -1;
  for (int32_t li = 0; li < k.nlev; li++) {
    const hz_kit_level *S = &k.lev[li];
    for (uint32_t t = 0; t < S->ntris; t++) {
      double cx = (S->vx[S->ti0[t]] + S->vx[S->ti1[t]] + S->vx[S->ti2[t]]) / 3.0;
      double cy = (S->vy[S->ti0[t]] + S->vy[S->ti1[t]] + S->vy[S->ti2[t]]) / 3.0;
      double cz = (S->vz[S->ti0[t]] + S->vz[S->ti1[t]] + S->vz[S->ti2[t]]) / 3.0;
      int take;
      if (radial) {
        double d = sqrt((cx - eye[0]) * (cx - eye[0]) + (cy - eye[1]) * (cy - eye[1]) +
                        (cz - eye[2]) * (cz - eye[2]));
        int32_t ring = (int32_t)(d / zone);
        if (ring > k.nlev - 1) ring = k.nlev - 1;
        take = (ring == li);
      } else {
        take = (li == 0) ? (cz < split) : (cz >= split);
      }
      if (!take) continue;
      M->ti0[w] = S->ti0[t] + vbase[li];
      M->ti1[w] = S->ti1[t] + vbase[li];
      M->ti2[w] = S->ti2[t] + vbase[li];
      M->tmtl[w] = S->tmtl[t];
      if (prev_src != li) {
        if (ncl > 0) M->cl[ncl - 1].ntris = w - M->cl[ncl - 1].first_tri;
        M->cl[ncl].first_tri = w;
        ncl++;
        prev_src = li;
      }
      M->tcl[w] = ncl - 1;
      w++;
    }
  }
  if (ncl > 0) M->cl[ncl - 1].ntris = w - M->cl[ncl - 1].first_tri;
  M->nclust = ncl;
  {
    hz_cluster *shr = realloc(M->cl, (size_t)ncl * sizeof(hz_cluster));
    if (shr != NULL) M->cl = shr; /* усечение вниз: отказ невозможен при ncl>0 */
  }
  for (uint32_t c = 0; c < ncl; c++) {
    hz_cluster *g = &M->cl[c];
    double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    int first = 1;
    for (uint32_t t = g->first_tri; t < g->first_tri + g->ntris; t++) {
      uint32_t vi[3] = {M->ti0[t], M->ti1[t], M->ti2[t]};
      for (int a = 0; a < 3; a++) {
        double p[3] = {M->vx[vi[a]], M->vy[vi[a]], M->vz[vi[a]]};
        for (int ax = 0; ax < 3; ax++) {
          if (first || p[ax] < lo[ax]) lo[ax] = p[ax];
          if (first || p[ax] > hi[ax]) hi[ax] = p[ax];
        }
        first = 0;
      }
    }
    for (int ax = 0; ax < 3; ax++) {
      g->bmin[ax] = (float)lo[ax];
      g->bmax[ax] = (float)hi[ax];
    }
    g->err = 0.0f;
  }
  rc = hz_kit_validate(&out);
  if (rc != HZ_KIT_OK) {
    fprintf(stderr, "kitmix: итог невалиден (rc=%d)\n", rc);
    hz_kit_free(&out);
    hz_kit_free(&k);
    return 2;
  }
  f = fopen(argv[2], "wb");
  if (f == NULL) {
    fprintf(stderr, "kitmix: не открывается %s\n", argv[2]);
    hz_kit_free(&out);
    hz_kit_free(&k);
    return 2;
  }
  rc = hz_kit_save(&out, f);
  if (fclose(f) != 0) rc = HZ_KIT_E_IO;
  printf("записан %s: nv=%u nt=%u nc=%u (rc=%d)\n", argv[2], M->nverts, M->ntris, M->nclust, rc);
  hz_kit_free(&out);
  hz_kit_free(&k);
  return rc == HZ_KIT_OK ? 0 : 2;
}
