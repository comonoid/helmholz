/* kitmk.c — §914-Ш4: OBJ → КИТ (вырожденный, nlev=1, F64).
 *
 * Паритет (§914-Ш4-П1): вершины хранятся double (флаг F64) и читаются
 * из OBJ побайтово; треугольники — в ПОРЯДКЕ OBJ; материалы — все поля
 * hz_obj_mtl, которые перенос использует (kd/kd3/ks3/ke3; ns/ior/alpha/
 * flat/tex в ките v2 НЕ хранятся — pgather их не читает, записано
 * ограничением). Грозди — чанки по HZ_KITMK_CHUNK треугольников,
 * покрытие ровное (А1583), bbox по вершинам чанка, ε=0 (L0).
 *
 * Синтаксис: kitmk IN.obj OUT.kit [scale=S]
 */
#include "geom/kit.h"
#include "scene_obj.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Размер грозди: целевой размер из бюджета §914 (128 тр); не магия —
 * масштаб чанка для чтения грубых уровней будущими Ш5+. */
#define HZ_KITMK_CHUNK 128

int main(int argc, char **argv) {
  double scale = 1.0;
  if (argc < 3) {
    fprintf(stderr, "kitmk IN.obj OUT.kit [scale=S]\n");
    return 2;
  }
  for (int i = 3; i < argc; i++)
    if (strncmp(argv[i], "scale=", 6) == 0) scale = atof(argv[i] + 6);

  hz_objmesh m;
  memset(&m, 0, sizeof m);
  if (hz_obj_load(&m, argv[1], scale) != 0) {
    fprintf(stderr, "kitmk: не читается %s\n", argv[1]);
    return 2;
  }
  printf("obj: nv=%d nt=%d mtl=%d degen=%lld quad=%lld\n", m.nv, m.nt, m.nmtl, (long long)m.ndegen,
         (long long)m.nquad);
  if (m.nt < 1 || m.nv < 3) {
    fprintf(stderr, "kitmk: пустой меш\n");
    hz_obj_free(&m);
    return 2;
  }

  hz_kit k;
  hz_kit_init(&k);
  k.nlev = 1;
  k.lev = calloc(1, sizeof(hz_kit_level));
  k.nmtl = (uint32_t)m.nmtl;
  k.mtl = calloc((size_t)m.nmtl, sizeof(hz_kit_mtl));
  if (k.lev == NULL || k.mtl == NULL) {
    fprintf(stderr, "kitmk: нет памяти\n");
    hz_kit_free(&k);
    hz_obj_free(&m);
    return 2;
  }
  k.flags = HZ_KIT_FLAG_F64; /* паритет: вершины без округления */

  hz_kit_level *L = &k.lev[0];
  L->nverts = (uint32_t)m.nv;
  L->ntris = (uint32_t)m.nt;
  L->nclust = (uint32_t)((m.nt + HZ_KITMK_CHUNK - 1) / HZ_KITMK_CHUNK);
  L->vx = malloc((size_t)m.nv * sizeof(double));
  L->vy = malloc((size_t)m.nv * sizeof(double));
  L->vz = malloc((size_t)m.nv * sizeof(double));
  L->ti0 = malloc(sizeof(uint32_t) * (size_t)m.nt);
  L->ti1 = malloc(sizeof(uint32_t) * (size_t)m.nt);
  L->ti2 = malloc(sizeof(uint32_t) * (size_t)m.nt);
  L->tcl = malloc(sizeof(uint32_t) * (size_t)m.nt);
  L->tmtl = malloc(sizeof(uint32_t) * (size_t)m.nt);
  L->cl = calloc(L->nclust, sizeof(hz_cluster));
  if (!L->vx || !L->vy || !L->vz || !L->ti0 || !L->ti1 || !L->ti2 || !L->tcl || !L->tmtl ||
      !L->cl) {
    fprintf(stderr, "kitmk: нет памяти (уровень)\n");
    hz_kit_free(&k);
    hz_obj_free(&m);
    return 2;
  }
  for (int32_t v = 0; v < m.nv; v++) {
    L->vx[v] = m.v[3 * (int64_t)v];
    L->vy[v] = m.v[3 * (int64_t)v + 1];
    L->vz[v] = m.v[3 * (int64_t)v + 2];
  }
  for (int32_t t = 0; t < m.nt; t++) {
    L->ti0[t] = (uint32_t)m.f[3 * (int64_t)t];
    L->ti1[t] = (uint32_t)m.f[3 * (int64_t)t + 1];
    L->ti2[t] = (uint32_t)m.f[3 * (int64_t)t + 2];
    L->tcl[t] = (uint32_t)(t / HZ_KITMK_CHUNK);
    L->tmtl[t] = (uint32_t)m.fm[t];
  }
  for (uint32_t c = 0; c < L->nclust; c++) {
    hz_cluster *g = &L->cl[c];
    g->first_tri = c * HZ_KITMK_CHUNK;
    g->ntris = L->ntris - g->first_tri;
    if (g->ntris > HZ_KITMK_CHUNK) g->ntris = HZ_KITMK_CHUNK;
    g->err = 0.0f;
    double lo[3] = {0, 0, 0}, hi[3] = {0, 0, 0};
    int first = 1;
    for (uint32_t t = g->first_tri; t < g->first_tri + g->ntris; t++) {
      const uint32_t idx[3] = {L->ti0[t], L->ti1[t], L->ti2[t]};
      for (int vi = 0; vi < 3; vi++) {
        double p[3] = {L->vx[idx[vi]], L->vy[idx[vi]], L->vz[idx[vi]]};
        for (int a = 0; a < 3; a++) {
          if (first || p[a] < lo[a]) lo[a] = p[a];
          if (first || p[a] > hi[a]) hi[a] = p[a];
        }
        first = 0;
      }
    }
    for (int a = 0; a < 3; a++) {
      g->bmin[a] = (float)lo[a];
      g->bmax[a] = (float)hi[a];
    }
  }
  for (int32_t mi = 0; mi < m.nmtl; mi++) {
    k.mtl[mi].kd = m.mtl[mi].kd;
    memcpy(k.mtl[mi].kd3, m.mtl[mi].kd3, 24);
    memcpy(k.mtl[mi].ks3, m.mtl[mi].ks3, 24);
    memcpy(k.mtl[mi].ke3, m.mtl[mi].ke3, 24);
  }

  int rc = hz_kit_validate(&k);
  if (rc != HZ_KIT_OK) {
    fprintf(stderr, "kitmk: кит не прошёл валидацию: %d\n", rc);
    hz_kit_free(&k);
    hz_obj_free(&m);
    return 2;
  }
  FILE *f = fopen(argv[2], "wb");
  if (f == NULL) {
    fprintf(stderr, "kitmk: не открывается %s\n", argv[2]);
    hz_kit_free(&k);
    hz_obj_free(&m);
    return 2;
  }
  rc = hz_kit_save(&k, f);
  if (fclose(f) != 0) rc = HZ_KIT_E_IO;
  printf("kit: nlev=1 (L0) verts=%u tris=%u clusters=%u mtl=%u → %s (%s)\n", L->nverts, L->ntris,
         L->nclust, k.nmtl, argv[2], rc == HZ_KIT_OK ? "ok" : "FAIL");
  long fsz = 0;
  f = fopen(argv[2], "rb");
  if (f != NULL) {
    fseek(f, 0, SEEK_END);
    fsz = ftell(f);
    fclose(f);
  }
  printf("размер файла: %ld Б (%.2f Б/треугольник)\n", fsz, m.nt > 0 ? (double)fsz / m.nt : 0.0);
  hz_kit_free(&k);
  hz_obj_free(&m);
  return rc == HZ_KIT_OK ? 0 : 2;
}
