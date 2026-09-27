/* kitpages.c — §914-Ш5, П6: touched-страницы mmap при чтении грубых
 * уровней НЕ включают страницы детальных. Синтетика: пол 16×16 кводов,
 * L0=512 тр + L1=32 (укрупнённая сетка) — как в kitwalk. Прибор:
 * mincore до/после hz_kit_load_range(min_lev).
 */
#include "geom/kit.h"

#include <fcntl.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/mman.h>
#include <unistd.h>

#define NQ 64 /* L0 крупнее страницы: блок L0 ~= 30 страниц */

static int fill_level(hz_kit_level *L, int nq) {
  L->nverts = (uint32_t)((nq + 1) * (nq + 1));
  L->ntris = (uint32_t)(2 * nq * nq);
  L->nclust = (uint32_t)((2 * nq * nq + 127) / 128);
  L->vx = malloc(L->nverts * sizeof(double));
  L->vy = malloc(L->nverts * sizeof(double));
  L->vz = malloc(L->nverts * sizeof(double));
  L->ti0 = malloc(4ull * L->ntris);
  L->ti1 = malloc(4ull * L->ntris);
  L->ti2 = malloc(4ull * L->ntris);
  L->tcl = malloc(4ull * L->ntris);
  L->tmtl = malloc(4ull * L->ntris);
  L->cl = calloc(L->nclust, sizeof(hz_cluster));
  if (!L->vx || !L->vy || !L->vz || !L->ti0 || !L->ti1 || !L->ti2 || !L->tcl || !L->tmtl || !L->cl)
    return 1;
  double q = 16.0 / nq;
  uint32_t nv = 0;
  for (int j = 0; j <= nq; j++)
    for (int i = 0; i <= nq; i++) {
      L->vx[nv] = i * q;
      L->vy[nv] = j * q;
      L->vz[nv] = 0;
      nv++;
    }
  uint32_t nt = 0;
  for (int j = 0; j < nq; j++)
    for (int i = 0; i < nq; i++) {
      uint32_t a = (uint32_t)(j * (nq + 1) + i), b = a + 1, c = a + (uint32_t)(nq + 1), d = c + 1;
      L->ti0[nt] = a;
      L->ti1[nt] = b;
      L->ti2[nt] = d;
      nt++;
      L->ti0[nt] = a;
      L->ti1[nt] = d;
      L->ti2[nt] = c;
      nt++;
    }
  for (uint32_t t = 0; t < L->ntris; t++) {
    L->tcl[t] = t / 128;
    L->tmtl[t] = 0;
  }
  for (uint32_t c = 0; c < L->nclust; c++) {
    hz_cluster *g = &L->cl[c];
    g->first_tri = c * 128;
    g->ntris = L->ntris - g->first_tri;
    if (g->ntris > 128) g->ntris = 128;
    double lo[3] = {0, 0, 0}, hi[3] = {16, 16, 0};
    for (int a = 0; a < 3; a++) {
      g->bmin[a] = (float)lo[a];
      g->bmax[a] = (float)hi[a];
    }
    g->err = 0.0f;
  }
  return 0;
}

static long count_resident(const void *addr, size_t len, long page) {
  size_t np = (len + (size_t)page - 1) / (size_t)page;
  unsigned char *vec = calloc(np, 1);
  if (vec == NULL) return -1;
  if (mincore((void *)((uintptr_t)addr & ~(uintptr_t)(page - 1)), np * (size_t)page, vec) != 0) {
    free(vec);
    return -1;
  }
  long r = 0;
  for (size_t i = 0; i < np; i++)
    if (vec[i] & 1) r++;
  free(vec);
  return r;
}

int main(void) {
  hz_kit k;
  hz_kit_init(&k);
  k.nlev = 2;
  k.flags = HZ_KIT_FLAG_F64;
  k.lev = calloc(2, sizeof(hz_kit_level));
  k.nmtl = 1;
  k.mtl = calloc(1, sizeof(hz_kit_mtl));
  if (k.lev == NULL || k.mtl == NULL) {
    hz_kit_free(&k);
    return 2;
  }
  k.mtl[0].kd = 0.5;
  if (fill_level(&k.lev[0], NQ) || fill_level(&k.lev[1], 4)) {
    hz_kit_free(&k);
    return 2;
  }
  if (hz_kit_validate(&k) != HZ_KIT_OK) {
    printf("kit invalid\n");
    hz_kit_free(&k);
    return 1;
  }
  const char *path = "build/pg2.kit";
  FILE *f = fopen(path, "wb");
  if (f == NULL || hz_kit_save(&k, f) != HZ_KIT_OK) {
    printf("save FAIL\n");
    hz_kit_free(&k);
    return 1;
  }
  fclose(f);
  FILE *g = fopen(path, "rb");
  if (g == NULL) return 2;
  fseek(g, 0, SEEK_END);
  long fsz = ftell(g);
  fclose(g);
  /* сбросить page cache после собственной записи — иначе база горячая
   * и прибор меряет кеш от save, не от загрузки */
  {
    int dfd = open(path, O_RDONLY);
    if (dfd >= 0) {
      /* грязные страницы fadvise не сбрасывает: сначала выsyncить */
      fsync(dfd);
      posix_fadvise(dfd, 0, 0, POSIX_FADV_DONTNEED);
      close(dfd);
    }
  }
  printf("kit: L0=%u тр, L1=%u тр, файл %ld Б\n", k.lev[0].ntris, k.lev[1].ntris, fsz);
  hz_kit_free(&k);

  /* ПРИНЦИП ПРИБОРА: загрузчик имеет собственное отображение, но page
   * cache ОБЩИЙ: страницы, которых коснулся загрузчик, resident и в
   * нашем контрольном отображении. База снимается до загрузки. */
  long page = sysconf(_SC_PAGESIZE);
  int fd = open(path, O_RDONLY);
  if (fd < 0) return 2;
  unsigned char *map = mmap(NULL, (size_t)fsz, PROT_READ, MAP_SHARED, fd, 0);
  close(fd);
  if (map == MAP_FAILED) return 2;
  madvise(map, (size_t)fsz, MADV_DONTNEED);
  long before = count_resident(map, (size_t)fsz, page);

  hz_kit kk;
  int rc = hz_kit_load_range(&kk, path, 1); /* только L1 */
  if (rc != HZ_KIT_OK) {
    printf("load_mmap(min_lev=1) rc=%d FAIL\n", rc);
    return 1;
  }
  long after = count_resident(map, (size_t)fsz, page);
  printf("П6: страниц resident: до=%ld после=%ld (страниц в файле %ld)\n", before, after,
         (fsz + page - 1) / page);
  printf("загружен L1: nv=%u nt=%u; L0 present=%d\n", kk.lev[1].nverts, kk.lev[1].ntris,
         kk.present[0]);
  hz_kit_free(&kk);

  /* Какие страницы тронуты: L0-блок начинается после заголовка+dir+mtl+L1.
   * Возьмём офсет L0 из dir-записи. */
  {
    unsigned char hdr[128];
    FILE *h = fopen(path, "rb");
    if (h == NULL) return 2;
    if (fread(hdr, 1, 128, h) != 128) {
      fclose(h);
      return 2;
    }
    fclose(h);
    uint64_t off0;
    memcpy(&off0, hdr + 64 + 24, 8);
    printf("блок L0 начинается на байте %llu (страница %lld)\n", (unsigned long long)off0,
           (long long)(off0 / (uint64_t)page));
    /* страницы, ЦЕЛИКОМ принадлежащие блоку L0: первая такая —
     * ceil(off0/page)+0, если L0 занимает её целиком; страница 0
     * общая (заголовок+оглавление+L1) и корректно resident. */
    size_t np = (size_t)((fsz + page - 1) / page);
    unsigned char *vec = calloc(np, 1);
    if (vec == NULL) return 2;
    if (mincore(map, np * (size_t)page, vec) != 0) return 2;
    size_t first_full = (size_t)((off0 + (uint64_t)page - 1) / (uint64_t)page);
    long res_l0 = 0, np_l0 = 0;
    for (size_t i = first_full; i < np; i++) {
      np_l0++;
      if (vec[i] & 1) res_l0++;
    }
    free(vec);
    printf("страниц ЦЕЛИКОМ L0: %ld, из них resident: %ld %s\n", np_l0, res_l0,
           res_l0 == 0 ? "OK (L0 не тронут)" : "FAIL (L0 тронут)");
  }
  munmap(map, (size_t)fsz);
  remove(path);
  return 0;
}
