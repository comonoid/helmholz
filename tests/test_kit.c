/* test_kit.c — Ш1 §914: формат кита v1. Позитив: save→load→битовое
 * тождество; порядок блоков (грубые уровни РАНЬШЕ L0 в файле, §2.2 v2).
 * НК: каждый инвариант формата мутируется по байтовому офсету и обязан
 * дать СВОЙ код отказа (урок Г26: тест, видящий только «отказ», ничего
 * не доказывает). */
#include "geom/kit.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int nfail = 0;

#define CHECK(cond, name)                                                                          \
  do {                                                                                             \
    if (cond)                                                                                      \
      printf("  ok   %s\n", name);                                                                 \
    else {                                                                                         \
      printf("  FAIL %s\n", name);                                                                 \
      nfail++;                                                                                     \
    }                                                                                              \
  } while (0)

/* Октаэдр: 6 вершин, 8 треугольников, 2 грозди по 4. */
static void fill_octa(hz_kit_level *L) {
  L->nverts = 6;
  L->ntris = 8;
  L->nclust = 2;
  static const double X[6] = {1, -1, 0, 0, 0, 0};
  static const double Y[6] = {0, 0, 1, -1, 0, 0};
  static const double Z[6] = {0, 0, 0, 0, 1, -1};
  /* 8 граней: верхние 4 + нижние 4, наружная ориентация */
  static const uint32_t T[8][3] = {{0, 2, 4}, {2, 1, 4}, {1, 3, 4}, {3, 0, 4},
                                   {2, 0, 5}, {1, 2, 5}, {3, 1, 5}, {0, 3, 5}};
  L->vx = malloc(6 * sizeof(double));
  L->vy = malloc(6 * sizeof(double));
  L->vz = malloc(6 * sizeof(double));
  L->ti0 = malloc((size_t)8 * 4);
  L->ti1 = malloc((size_t)8 * 4);
  L->ti2 = malloc((size_t)8 * 4);
  L->tcl = malloc((size_t)8 * 4);
  L->tmtl = malloc((size_t)8 * 4);
  L->cl = malloc(2 * sizeof(hz_cluster));
  if (!L->vx || !L->vy || !L->vz || !L->ti0 || !L->ti1 || !L->ti2 || !L->tcl || !L->tmtl || !L->cl)
    exit(2);
  memcpy(L->vx, X, sizeof X);
  memcpy(L->vy, Y, sizeof Y);
  memcpy(L->vz, Z, sizeof Z);
  for (int t = 0; t < 8; t++)
    L->tmtl[t] = 0;
  for (int t = 0; t < 8; t++) {
    L->ti0[t] = T[t][0];
    L->ti1[t] = T[t][1];
    L->ti2[t] = T[t][2];
    L->tcl[t] = (uint32_t)(t / 4);
  }
  for (int c = 0; c < 2; c++) {
    hz_cluster *g = &L->cl[c];
    g->first_tri = (uint32_t)(4 * c);
    g->ntris = 4;
    g->err = 0.0f;
    g->bmin[0] = -1;
    g->bmin[1] = -1;
    g->bmin[2] = -1;
    g->bmax[0] = 1;
    g->bmax[1] = 1;
    g->bmax[2] = 1;
  }
}

/* Грубый уровень L1: октаэдр из 4 «четвертных» треугольников (замещение). */
static void fill_octa_coarse(hz_kit_level *L) {
  fill_octa(L); /* затем прореживаем: 8 → 4 треугольника, 1 гроздь */
  L->ntris = 4;
  L->nclust = 1;
  L->tcl[0] = L->tcl[1] = L->tcl[2] = L->tcl[3] = 0;
  L->cl[0].first_tri = 0;
  L->cl[0].ntris = 4;
}

static int same_level(const hz_kit_level *a, const hz_kit_level *b) {
  if (a->nverts != b->nverts || a->ntris != b->ntris || a->nclust != b->nclust) return 0;
  if (memcmp(a->vx, b->vx, 4ull * a->nverts) || memcmp(a->vy, b->vy, 4ull * a->nverts) ||
      memcmp(a->vz, b->vz, 4ull * a->nverts))
    return 0;
  if (memcmp(a->ti0, b->ti0, 4ull * a->ntris) || memcmp(a->ti1, b->ti1, 4ull * a->ntris) ||
      memcmp(a->ti2, b->ti2, 4ull * a->ntris) || memcmp(a->tcl, b->tcl, 4ull * a->ntris))
    return 0;
  for (uint32_t c = 0; c < a->nclust; c++)
    if (memcmp(&a->cl[c], &b->cl[c], sizeof(hz_cluster))) return 0;
  return 1;
}

static const char *path = "build/kit_tmp.bin";

static int save_kit(const hz_kit *k) {
  FILE *f = fopen(path, "wb");
  if (f == NULL) return -1;
  int rc = hz_kit_save(k, f);
  if (fclose(f) != 0) rc = HZ_KIT_E_IO;
  return rc;
}

/* Мутирует байтовый офсет в сохранённом файле и перезагружает. Байт
 * ВОССТАВЛИВАЕТСЯ: порчи не накапливаются, каждая проверка видит СВОЮ
 * мутацию (иначе всё падает на первой же — урок этой сессии). Каждый
 * FILE* живёт в своём блоке с безусловным закрытием — так gcc-analyzer
 * видит закрытие на всех путях (смешивание в одном указателе теряет
 * его, известный класс FPs). */
static int reload_mut(uint32_t off, uint32_t val32) {
  int rc = -1;
  uint32_t orig = 0;
  int mutated = 0;
  FILE *f = fopen(path, "r+b");
  if (f != NULL) {
    if (fseek(f, (long)off, SEEK_SET) == 0 && fread(&orig, 4, 1, f) == 1 &&
        fseek(f, (long)off, SEEK_SET) == 0 && fwrite(&val32, 4, 1, f) == 1)
      mutated = 1;
    fclose(f);
  }
  if (!mutated) return rc;
  {
    FILE *g = fopen(path, "rb");
    if (g != NULL) {
      hz_kit k;
      rc = hz_kit_load(&k, g);
      fclose(g);
      if (rc == HZ_KIT_OK) hz_kit_free(&k); /* неожиданно валиден — утечки не нужно */
    }
  }
  {
    FILE *h = fopen(path, "r+b");
    if (h != NULL) {
      if (fseek(h, (long)off, SEEK_SET) != 0 || fwrite(&orig, 4, 1, h) != 1)
        rc = -1; /* файл не восстановлен — дальше не годен */
      fclose(h);
    } else {
      rc = -1;
    }
  }
  return rc;
}

int main(void) {
  printf("kit: позитив — save/load, битовое тождество\n");
  hz_kit k;
  hz_kit_init(&k);
  k.nlev = 2;
  k.lev = calloc(2, sizeof(hz_kit_level));
  if (k.lev == NULL) return 2;
  k.nmtl = 1;
  k.mtl = calloc(1, sizeof(hz_kit_mtl));
  if (k.mtl == NULL) {
    hz_kit_free(&k);
    return 2;
  }
  k.mtl[0].kd = 0.5;
  for (int c = 0; c < 3; c++) {
    k.mtl[0].kd3[c] = 0.5;
    k.mtl[0].ks3[c] = 0.1;
    k.mtl[0].ke3[c] = 0.0;
  }
  fill_octa(&k.lev[0]);
  fill_octa_coarse(&k.lev[1]);
  CHECK(hz_kit_validate(&k) == HZ_KIT_OK, "валидность 2-уровневого кита");
  CHECK(save_kit(&k) == HZ_KIT_OK, "save");

  hz_kit l;
  FILE *f = fopen(path, "rb");
  CHECK(f != NULL, "файл открыт");
  if (f == NULL) return 2;
  int rc = hz_kit_load(&l, f);
  fclose(f);
  CHECK(rc == HZ_KIT_OK, "load");
  if (rc == HZ_KIT_OK) {
    CHECK(l.nlev == 2 && l.nmtl == 1, "nlev/nmtl");
    CHECK(l.mtl != NULL && memcmp(&l.mtl[0].kd, &k.mtl[0].kd, sizeof(hz_kit_mtl)) == 0,
          "материал битово");
    CHECK(same_level(&k.lev[0], &l.lev[0]) && same_level(&k.lev[1], &l.lev[1]),
          "битовое тождество уровней");
  }

  printf("kit: порядок блоков — грубый раньше L0\n");
  {
    /* оглавление: L0 (i=0) на офсете 64, L1 (i=1) на 64+32; офсет — u64
     * на +24 внутри записи. */
    unsigned char hdr[128];
    FILE *g = fopen(path, "rb");
    if (g == NULL || fread(hdr, 1, 128, g) != 128) {
      printf("  FAIL чтение оглавления\n");
      if (g != NULL) fclose(g);
      return 2;
    }
    fclose(g);
    uint64_t off0, off1;
    memcpy(&off0, hdr + 64 + 24, 8);
    memcpy(&off1, hdr + 64 + 32 + 24, 8);
    CHECK(off1 < off0, "данные L1 лежат в файле раньше L0 (§2.2 v2)");
  }

  printf("kit: раздельность в памяти — load_part (грубые без L0)\n");
  {
    hz_kit p;
    FILE *g = fopen(path, "rb");
    if (g == NULL) return 2;
    int prc = hz_kit_load_part(&p, g, 1); /* только уровни [1, nlev) */
    fclose(g);
    CHECK(prc == HZ_KIT_OK, "load_part(min_lev=1)");
    if (prc == HZ_KIT_OK && p.present != NULL) {
      CHECK(p.nlev == 2, "present-таблица есть, nlev цел");
      CHECK(p.present[1] == 1, "L1 присутствует");
      CHECK(p.present[0] == 0, "L0 НЕ загружен (раздельность §914-R2)");
      CHECK(p.lev[0].vx == NULL && p.lev[0].cl == NULL, "память L0 не занята");
      CHECK(same_level(&k.lev[1], &p.lev[1]), "L1 битово тот же");
    } else if (prc == HZ_KIT_OK) {
      CHECK(0, "present-таблица есть");
    }
    if (prc == HZ_KIT_OK) hz_kit_free(&p);
  }

  printf("kit: Ш5 — mmap-путь битово равен fread-пути\n");
  {
    hz_kit p2;
    int mrc = hz_kit_load_range(&p2, path, 0);
    CHECK(mrc == HZ_KIT_OK, "load_mmap(все уровни)");
    if (mrc == HZ_KIT_OK) {
      CHECK(p2.nlev == 2 && same_level(&k.lev[0], &p2.lev[0]) && same_level(&k.lev[1], &p2.lev[1]),
            "mmap: уровни битово те же");
      CHECK(p2.nmtl == k.nmtl && memcmp(p2.mtl, k.mtl, k.nmtl * sizeof(hz_kit_mtl)) == 0,
            "mmap: материалы битово те же");
      hz_kit_free(&p2);
    }
    mrc = hz_kit_load_range(&p2, path, 1); /* грубые без L0 */
    CHECK(mrc == HZ_KIT_OK && p2.present[0] == 0 && p2.present[1] == 1, "load_mmap(min_lev=1)");
    if (mrc == HZ_KIT_OK) hz_kit_free(&p2);
  }

  printf("kit: НК — каждый инвариант даёт СВОЙ код\n");
  save_kit(&k);
  CHECK(reload_mut(HZ_KIT_OFF_MAGIC, 0xdeadbeefu) == HZ_KIT_E_MAGIC, "magic → E_MAGIC");
  CHECK(reload_mut(HZ_KIT_OFF_VERSION, HZ_KIT_VERSION + 1u) == HZ_KIT_E_VERSION,
        "версия → E_VERSION");
  CHECK(reload_mut(HZ_KIT_OFF_FLAGS, 0x80000000u) == HZ_KIT_E_FLAGS, "флаги → E_FLAGS");
  CHECK(reload_mut(HZ_KIT_OFF_ENDIAN, 0x76543210u) == HZ_KIT_E_ENDIAN, "endian → E_ENDIAN");
  CHECK(reload_mut(HZ_KIT_OFF_CANARY, 0) == HZ_KIT_E_CANARY, "канарейка → E_CANARY");
  CHECK(reload_mut(HZ_KIT_OFF_RESERVED, 1u) == HZ_KIT_E_RESERVED, "reserved → E_RESERVED");
  CHECK(reload_mut(HZ_KIT_OFF_NLEV, 0u) == HZ_KIT_E_NLEV, "nlev=0 → E_NLEV");
  CHECK(reload_mut(HZ_KIT_OFF_NLEV, 65u) == HZ_KIT_E_NLEV, "nlev=65 → E_NLEV");

  /* счётчик ntris уровня L0 (первая запись оглавления, +4): 0 → E_COUNT */
  CHECK(reload_mut(64 + 4, 0u) == HZ_KIT_E_COUNT, "ntris=0 → E_COUNT");
  CHECK(reload_mut(64 + 4, HZ_KIT_MAX_COUNT + 1u) == HZ_KIT_E_COUNT, "ntris>кэп → E_COUNT");
  /* офсет данных L0 (запись 0, +24..32) за EOF → E_RANGE */
  CHECK(reload_mut(64 + 24, 0xfffffff0u) == HZ_KIT_E_RANGE, "офсет за EOF → E_RANGE");

  /* структурные НК — через повторный save мутированного кита */
  {
    hz_kit bad;
    hz_kit_init(&bad);
    bad.nlev = 1;
    bad.nmtl = 1;
    bad.mtl = calloc(1, sizeof(hz_kit_mtl));
    if (bad.mtl == NULL) {
      hz_kit_free(&bad);
      return 2;
    }
    bad.mtl[0].kd = 0.5;
    bad.lev = calloc(1, sizeof(hz_kit_level));
    if (bad.lev == NULL) {
      hz_kit_free(&bad);
      return 2;
    }
    fill_octa(&bad.lev[0]);
    bad.lev[0].ti0[3] = 6; /* ≥ nverts */
    CHECK(hz_kit_validate(&bad) == HZ_KIT_E_TRIDX, "индекс вершины → E_TRIDX");
    bad.lev[0].ti0[3] = 0;
    bad.lev[0].tmtl[3] = 1; /* ≥ nmtl=1 */
    CHECK(hz_kit_validate(&bad) == HZ_KIT_E_MTLIDX, "tmtl ≥ nmtl → E_MTLIDX");
    bad.lev[0].tmtl[3] = 0;
    bad.lev[0].cl[0].ntris = 3; /* покрытие дырявое */
    CHECK(hz_kit_validate(&bad) == HZ_KIT_E_COVER, "дыра в покрытии → E_COVER");
    bad.lev[0].cl[0].ntris = 4;
    bad.lev[0].vx[2] = NAN;
    CHECK(hz_kit_validate(&bad) == HZ_KIT_E_NAN, "NaN-вершина → E_NAN");
    bad.lev[0].vx[2] = 0;
    bad.nlev = 2;
    hz_kit_level *lev2 = realloc(bad.lev, 2 * sizeof(hz_kit_level));
    if (lev2 == NULL) return 2;
    bad.lev = lev2;
    memset(&bad.lev[1], 0, sizeof(hz_kit_level));
    fill_octa(&bad.lev[1]); /* L1 такой же ДЕТАЛЬНЫЙ+ — ntris не убывает? */
    /* равный ntris допустим (не возрастает), ломаем индексом грозди */
    bad.lev[1].tcl[0] = 2;
    CHECK(hz_kit_validate(&bad) == HZ_KIT_E_CLIDX, "tcl ≥ nclust → E_CLIDX");
    bad.lev[1].tcl[0] = 0;
    /* лестница: L0 грубее L1 нельзя — пересобираем L0 начисто */
    free(bad.lev[0].vx);
    free(bad.lev[0].vy);
    free(bad.lev[0].vz);
    free(bad.lev[0].ti0);
    free(bad.lev[0].ti1);
    free(bad.lev[0].ti2);
    free(bad.lev[0].tcl);
    free(bad.lev[0].tmtl);
    free(bad.lev[0].cl);
    fill_octa(&bad.lev[0]);
    bad.lev[0].ntris = 4;
    bad.lev[0].nclust = 1;
    bad.lev[0].cl[0].first_tri = 0;
    bad.lev[0].cl[0].ntris = 4;
    CHECK(hz_kit_validate(&bad) == HZ_KIT_E_LADDER, "L0 грубее L1 → E_LADDER");
    hz_kit_free(&bad);
  }

  hz_kit_free(&k);
  if (rc == HZ_KIT_OK) hz_kit_free(&l);
  remove(path);
  printf(nfail ? "kit: %d FAIL\n" : "kit: ALL OK\n", nfail);
  return nfail ? 1 : 0;
}
