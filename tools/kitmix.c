/* kitmix.c — §914-П3: смешанный кит — L0 ближней половины + L1 дальней.
 *
 * Вход: 2-уровневый кит (после kitdec). Выход: ОДНОуровневый кит,
 * геометрия = L0-треугольники с центроидом z < split плюс
 * L1-треугольники с центроидом z >= split. Аппроксимация камерной
 * зоны для измерения ошибки грубости переносом (без многоуровневого
 * свипа — изоляция ошибки).
 *
 * Синтаксис: kitmix IN.kit OUT.kit split=S
 */
#include "geom/kit.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
  double split = 0.0;
  if (argc < 3) {
    fprintf(stderr, "kitmix IN.kit OUT.kit split=S\n");
    return 2;
  }
  for (int i = 3; i < argc; i++)
    if (strncmp(argv[i], "split=", 6) == 0) split = atof(argv[i] + 6);

  hz_kit k;
  hz_kit_init(&k);
  FILE *f = fopen(argv[1], "rb");
  if (f == NULL) {
    fprintf(stderr, "kitmix: не читается %s\n", argv[1]);
    return 2;
  }
  int rc = hz_kit_load(&k, f);
  fclose(f);
  if (rc != HZ_KIT_OK || k.nlev < 2) {
    fprintf(stderr, "kitmix: нужен кит с nlev>=2 (rc=%d nlev=%d)\n", rc, k.nlev);
    hz_kit_free(&k);
    return 2;
  }
  const hz_kit_level *A = &k.lev[0]; /* детальный */
  const hz_kit_level *B = &k.lev[1]; /* грубый */

  /* Счёт треугольников по половинам. */
  uint32_t na = 0, nb = 0;
  for (uint32_t t = 0; t < A->ntris; t++) {
    double cz = (A->vz[A->ti0[t]] + A->vz[A->ti1[t]] + A->vz[A->ti2[t]]) / 3.0;
    if (cz < split) na++;
  }
  for (uint32_t t = 0; t < B->ntris; t++) {
    double cz = (B->vz[B->ti0[t]] + B->vz[B->ti1[t]] + B->vz[B->ti2[t]]) / 3.0;
    if (cz >= split) nb++;
  }
  uint32_t nt = na + nb;
  printf("смешивание: L0-близко=%u L1-далеко=%u итого=%u (split=%.3f)\n", na, nb, nt, split);
  if (nt == 0) {
    fprintf(stderr, "kitmix: пустой результат\n");
    hz_kit_free(&k);
    return 2;
  }

  /* Выходной кит: вершины обоих уровней (без переиспользования:
   * смещение индексов B на A->nverts), треугольники подряд,
   * грозди-прогоны по источнику. */
  uint32_t nv = A->nverts + B->nverts;
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
  memcpy(M->vx, A->vx, (size_t)A->nverts * sizeof(double));
  memcpy(M->vy, A->vy, (size_t)A->nverts * sizeof(double));
  memcpy(M->vz, A->vz, (size_t)A->nverts * sizeof(double));
  memcpy(M->vx + A->nverts, B->vx, (size_t)B->nverts * sizeof(double));
  memcpy(M->vy + A->nverts, B->vy, (size_t)B->nverts * sizeof(double));
  memcpy(M->vz + A->nverts, B->vz, (size_t)B->nverts * sizeof(double));

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
  int prev_src = -1; /* 0 = L0, 1 = L1 */
  for (int pass = 0; pass < 2; pass++) {
    const hz_kit_level *S = pass == 0 ? A : B;
    uint32_t base = pass == 0 ? 0 : A->nverts;
    for (uint32_t t = 0; t < S->ntris; t++) {
      double cz = (S->vz[S->ti0[t]] + S->vz[S->ti1[t]] + S->vz[S->ti2[t]]) / 3.0;
      if (pass == 0 ? !(cz < split) : !(cz >= split)) continue;
      M->ti0[w] = S->ti0[t] + base;
      M->ti1[w] = S->ti1[t] + base;
      M->ti2[w] = S->ti2[t] + base;
      M->tmtl[w] = S->tmtl[t]; /* материал пережил декимацию (Ш2-фикс) */
      if (prev_src != pass) {
        if (ncl > 0) M->cl[ncl - 1].ntris = w - M->cl[ncl - 1].first_tri;
        M->cl[ncl].first_tri = w;
        ncl++;
        prev_src = pass;
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
