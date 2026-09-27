/* kit.c — ГЕОМЕТРИЯ v2, Ш1 §914: кит с лестницей замещающих уровней.
 * Формат и инварианты — kit.h; аудит до кода — §914-А (А1580–А1585).
 *
 * Гейт: clang-format, ccheck CLEAN, санитайзеры (tests/test_kit.c).
 * Горячих путей нет (формат/сериализация); аллокации инлайном у места
 * использования — без обёрток (урок gcc-analyzer, CLAUDE.md гейт 2).
 */
#include "kit.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Свидетели формата: те же идеи, что octree Г28 — версия ловит только
 * изменения, о которых ПОМНИЛИ; byte-order и IEEE754 не трогают версию. */
#define HZ_KIT_ENDIAN_WITNESS 0x01234567u
/* Канарейка: double с "неудобным" полным представлением бит. */
#define HZ_KIT_CANARY_BITS 0x4005BF0A8B145769ull /* ≈ 2.718281828… */
#define HZ_KIT_ALIGN 64u

static const char HZ_KIT_MAGIC[8] = {'H', 'Z', 'K', 'I', 'T', '1', 0, 0};

static double kit_canary(void) {
  double d;
  uint64_t b = HZ_KIT_CANARY_BITS;
  memcpy(&d, &b, sizeof d);
  return d;
}

void hz_kit_init(hz_kit *k) {
  k->nlev = 0;
  k->lev = NULL;
  k->flags = 0;
  k->present = NULL;
}

void hz_kit_free(hz_kit *k) {
  if (k == NULL) return;
  if (k->lev != NULL) {
    for (int32_t i = 0; i < k->nlev; i++) {
      if (k->present != NULL && !k->present[i]) continue; /* не загружен */
      hz_kit_level *L = &k->lev[i];
      free(L->vx);
      free(L->vy);
      free(L->vz);
      free(L->ti0);
      free(L->ti1);
      free(L->ti2);
      free(L->tcl);
      free(L->tmtl);
      free(L->cl);
    }
    free(k->lev);
  }
  free(k->present);
  free(k->mtl);
  k->nlev = 0;
  k->lev = NULL;
  k->present = NULL;
  k->mtl = NULL;
  k->nmtl = 0;
}

/* --- проверка структур в памяти ---------------------------------------- */

static int level_validate(const hz_kit_level *L, uint32_t nmtl) {
  if (L->nverts < 3u || L->nverts > HZ_KIT_MAX_COUNT) return HZ_KIT_E_COUNT;
  if (L->ntris == 0u || L->ntris > HZ_KIT_MAX_COUNT) return HZ_KIT_E_COUNT;
  if (L->nclust == 0u || L->nclust > HZ_KIT_MAX_COUNT) return HZ_KIT_E_COUNT;
  for (uint32_t t = 0; t < L->ntris; t++) {
    if (L->ti0[t] >= L->nverts || L->ti1[t] >= L->nverts || L->ti2[t] >= L->nverts)
      return HZ_KIT_E_TRIDX;
    if (L->tcl[t] >= L->nclust) return HZ_KIT_E_CLIDX;
    if (L->tmtl[t] >= nmtl) return HZ_KIT_E_MTLIDX;
  }
  /* Грозди покрывают [0, ntris) РОВНО, упорядочены (А1583). */
  uint32_t next = 0;
  for (uint32_t c = 0; c < L->nclust; c++) {
    const hz_cluster *g = &L->cl[c];
    if (g->first_tri != next) return HZ_KIT_E_COVER;
    if (g->ntris == 0u || g->ntris > L->ntris - g->first_tri) return HZ_KIT_E_COVER;
    next = g->first_tri + g->ntris;
  }
  if (next != L->ntris) return HZ_KIT_E_COVER;
  /* tcl согласован с диапазоном своей грозди. */
  for (uint32_t c = 0; c < L->nclust; c++)
    for (uint32_t t = L->cl[c].first_tri; t < L->cl[c].first_tri + L->cl[c].ntris; t++)
      if (L->tcl[t] != c) return HZ_KIT_E_COVER;
  for (uint32_t v = 0; v < L->nverts; v++)
    if (!isfinite(L->vx[v]) || !isfinite(L->vy[v]) || !isfinite(L->vz[v])) return HZ_KIT_E_NAN;
  for (uint32_t c = 0; c < L->nclust; c++) {
    const hz_cluster *g = &L->cl[c];
    for (int a = 0; a < 3; a++)
      if (!isfinite(g->bmin[a]) || !isfinite(g->bmax[a])) return HZ_KIT_E_NAN;
    if (!isfinite(g->err)) return HZ_KIT_E_NAN;
  }
  return HZ_KIT_OK;
}

int hz_kit_validate(const hz_kit *k) {
  if (k == NULL || k->lev == NULL) return HZ_KIT_E_ARG;
  if (k->nlev < 1 || k->nlev > HZ_KIT_MAX_NLEV) return HZ_KIT_E_NLEV;
  int rc = HZ_KIT_OK;
  int have_prev = 0;
  uint32_t prev_ntris = 0;
  for (int32_t i = 0; i < k->nlev; i++) {
    if (k->present != NULL && !k->present[i]) continue; /* не загружен */
    rc = level_validate(&k->lev[i], k->nmtl);
    if (rc != HZ_KIT_OK) return rc;
    /* Замещение: вверх по лестнице НЕ грубеет (А1582); по присутствующим. */
    if (!have_prev)
      have_prev = 1;
    else if (k->lev[i].ntris > prev_ntris)
      return HZ_KIT_E_LADDER;
    prev_ntris = k->lev[i].ntris;
  }
  return rc;
}

/* --- сериализация ------------------------------------------------------- */

static int wr(FILE *f, const void *p, size_t n) {
  return fwrite(p, 1, n, f) == n ? 0 : HZ_KIT_E_IO;
}

static int wr_u32(FILE *f, uint32_t v) {
  return wr(f, &v, 4);
}
static int wr_u64(FILE *f, uint64_t v) {
  return wr(f, &v, 8);
}

static int pad_to(FILE *f, uint64_t from, uint64_t align) {
  uint64_t n = (align - (from % align)) % align;
  if (n == 0) return 0;
  unsigned char z[64] = {0};
  if (n > sizeof z) return HZ_KIT_E_ARG; /* align ≤ 64 — аргумент, а не файл */
  return wr(f, z, (size_t)n);
}

static uint64_t level_bytes(const hz_kit_level *L, int f64) {
  uint64_t nv = L->nverts, nt = L->ntris, nc = L->nclust;
  uint64_t vb = f64 ? 8ull : 4ull;
  return vb * nv * 3ull + 4ull * nt * 5ull + 40ull * nc; /* verts, idx/tcl/tmtl, грозди */
}

int hz_kit_save(const hz_kit *k, FILE *f) {
  if (k == NULL || f == NULL || k->lev == NULL) return HZ_KIT_E_ARG;
  int rc = hz_kit_validate(k);
  if (rc != HZ_KIT_OK) return rc;

  unsigned char h[HZ_KIT_HDRSIZE] = {0};
  memcpy(h + HZ_KIT_OFF_MAGIC, HZ_KIT_MAGIC, 8);
  uint32_t ver = HZ_KIT_VERSION, fl = k->flags & (HZ_KIT_FLAG_PAD64 | HZ_KIT_FLAG_F64),
           en = HZ_KIT_ENDIAN_WITNESS, nl = (uint32_t)k->nlev;
  memcpy(h + HZ_KIT_OFF_VERSION, &ver, 4);
  memcpy(h + HZ_KIT_OFF_FLAGS, &fl, 4);
  memcpy(h + HZ_KIT_OFF_ENDIAN, &en, 4);
  memcpy(h + HZ_KIT_OFF_NLEV, &nl, 4);
  double ca = kit_canary();
  memcpy(h + HZ_KIT_OFF_CANARY, &ca, 8);
  uint32_t nm = k->nmtl;
  memcpy(h + HZ_KIT_OFF_NMTL, &nm, 4); /* [44..64) остаётся нулевым */
  rc = wr(f, h, HZ_KIT_HDRSIZE);
  if (rc) return rc;

  /* Оглавление: уровни по порядку i=0(L0)..nlev-1; офсет считается от
   * начала файла, данные пишутся ОТ ГРУБОГО К L0 (kit.h, §2.2 v2). */
  uint64_t dir_end = (uint64_t)HZ_KIT_HDRSIZE + (uint64_t)HZ_KIT_DIRSIZE * (uint64_t)k->nlev;
  uint64_t mtl_off = (dir_end + HZ_KIT_ALIGN - 1) / HZ_KIT_ALIGN * HZ_KIT_ALIGN;
  uint64_t mtl_bytes = (uint64_t)sizeof(hz_kit_mtl) * (uint64_t)k->nmtl; /* 80 Б/запись */
  uint64_t data_start = (mtl_off + mtl_bytes + HZ_KIT_ALIGN - 1) / HZ_KIT_ALIGN * HZ_KIT_ALIGN;
  uint64_t off = data_start;
  uint64_t offs[HZ_KIT_MAX_NLEV];
  for (int32_t i = k->nlev - 1; i >= 0; i--) { /* грубые раньше */
    offs[i] = off;
    off += level_bytes(&k->lev[i], (k->flags & HZ_KIT_FLAG_F64) != 0);
    off = (off + HZ_KIT_ALIGN - 1) / HZ_KIT_ALIGN * HZ_KIT_ALIGN;
  }
  for (int32_t i = 0; i < k->nlev; i++) {
    const hz_kit_level *L = &k->lev[i];
    uint32_t pad0 = 0;
    /* запись оглавления — ровно HZ_KIT_DIRSIZE=32 Б: 3 счётчика +
     * 3 паддинга + офсет u64 (ABI-независимость, урок Г5). */
    rc = wr_u32(f, L->nverts) || wr_u32(f, L->ntris) || wr_u32(f, L->nclust) || wr_u32(f, pad0) ||
         wr_u32(f, pad0) || wr_u32(f, pad0) || wr_u64(f, offs[i]);
    if (rc) return rc;
  }
  rc = pad_to(f, (uint64_t)HZ_KIT_HDRSIZE + HZ_KIT_DIRSIZE * (uint64_t)k->nlev, HZ_KIT_ALIGN);
  if (rc) return rc;
  /* таблица материалов (v2): сразу после выровненного оглавления */
  for (uint32_t mi = 0; mi < k->nmtl && rc == 0; mi++) {
    const hz_kit_mtl *mm = &k->mtl[mi];
    rc = wr(f, &mm->kd, 8) || wr(f, mm->kd3, 24) || wr(f, mm->ks3, 24) || wr(f, mm->ke3, 24);
  }
  if (rc) return rc;

  for (int32_t i = k->nlev - 1; i >= 0; i--) {
    const hz_kit_level *L = &k->lev[i];
    long here = ftell(f);
    if (here < 0) return HZ_KIT_E_IO;
    if (pad_to(f, (uint64_t)here, HZ_KIT_ALIGN)) return HZ_KIT_E_IO;
    if (k->flags & HZ_KIT_FLAG_F64) {
      rc = wr(f, L->vx, 8ull * L->nverts) || wr(f, L->vy, 8ull * L->nverts) ||
           wr(f, L->vz, 8ull * L->nverts);
      if (rc) return rc;
    } else {
      /* хранение float, в памяти double — конверсия при записи */
      float *tmp = malloc((size_t)L->nverts * sizeof *tmp);
      if (tmp == NULL) return HZ_KIT_E_MEM;
      for (int c = 0; c < 3; c++) {
        const double *src = c == 0 ? L->vx : (c == 1 ? L->vy : L->vz);
        for (uint32_t v = 0; v < L->nverts; v++)
          tmp[v] = (float)src[v];
        if (wr(f, tmp, 4ull * L->nverts)) {
          free(tmp);
          return HZ_KIT_E_IO;
        }
      }
      free(tmp);
    }
    rc = wr(f, L->ti0, 4ull * L->ntris) || wr(f, L->ti1, 4ull * L->ntris) ||
         wr(f, L->ti2, 4ull * L->ntris) || wr(f, L->tcl, 4ull * L->ntris) ||
         wr(f, L->tmtl, 4ull * L->ntris);
    if (rc) return rc;
    for (uint32_t c = 0; c < L->nclust && rc == 0; c++) {
      const hz_cluster *g = &L->cl[c];
      rc = wr(f, g->bmin, 12) || wr(f, g->bmax, 12) || wr(f, &g->err, 4) ||
           wr_u32(f, g->first_tri) || wr_u32(f, g->ntris) || wr(f, h, 4);
    }
    if (rc) return rc;
  }
  return HZ_KIT_OK;
}

/* --- загрузка ------------------------------------------------------------ */

int hz_kit_hdr_decode(const unsigned char h[HZ_KIT_HDRSIZE], int32_t *nlev, uint32_t *flags) {
  if (memcmp(h + HZ_KIT_OFF_MAGIC, HZ_KIT_MAGIC, 8) != 0) return HZ_KIT_E_MAGIC;
  uint32_t ver, fl, en, nl;
  memcpy(&ver, h + HZ_KIT_OFF_VERSION, 4);
  memcpy(&fl, h + HZ_KIT_OFF_FLAGS, 4);
  memcpy(&en, h + HZ_KIT_OFF_ENDIAN, 4);
  memcpy(&nl, h + HZ_KIT_OFF_NLEV, 4);
  if (ver != HZ_KIT_VERSION) return HZ_KIT_E_VERSION;
  if (fl & ~(HZ_KIT_FLAG_PAD64 | HZ_KIT_FLAG_F64)) return HZ_KIT_E_FLAGS;
  /* v2: [44,64) зарезервированы нулями (nmtl читает загрузчик) */
  for (int b = HZ_KIT_OFF_NMTL + 4; b < HZ_KIT_HDRSIZE; b++)
    if (h[b] != 0) return HZ_KIT_E_RESERVED;
  if (en != HZ_KIT_ENDIAN_WITNESS) return HZ_KIT_E_ENDIAN;
  double ca, want = kit_canary();
  memcpy(&ca, h + HZ_KIT_OFF_CANARY, 8);
  if (memcmp(&ca, &want, 8) != 0) return HZ_KIT_E_CANARY;
  uint64_t r;
  memcpy(&r, h + HZ_KIT_OFF_RESERVED, 8);
  if (r != 0) return HZ_KIT_E_RESERVED;
  if (nl < 1u || nl > (uint32_t)HZ_KIT_MAX_NLEV) return HZ_KIT_E_NLEV;
  *nlev = (int32_t)nl;
  *flags = fl;
  return HZ_KIT_OK;
}

static int rd(FILE *f, void *p, size_t n) {
  return fread(p, 1, n, f) == n ? 0 : HZ_KIT_E_IO;
}

static void level_clear(hz_kit_level *L) {
  free(L->vx);
  free(L->vy);
  free(L->vz);
  free(L->ti0);
  free(L->ti1);
  free(L->ti2);
  free(L->tcl);
  free(L->tmtl);
  free(L->cl);
  memset(L, 0, sizeof *L);
}

static int level_load(hz_kit_level *L, FILE *f, uint64_t fsize, int f64, uint32_t nmtl) {
  uint32_t hdr[6];
  if (rd(f, hdr, 24)) return HZ_KIT_E_IO;
  uint64_t off;
  if (rd(f, &off, 8)) return HZ_KIT_E_IO;
  L->nverts = hdr[0];
  L->ntris = hdr[1];
  L->nclust = hdr[2];
  if (L->nverts < 3u || L->nverts > HZ_KIT_MAX_COUNT || L->ntris == 0u ||
      L->ntris > HZ_KIT_MAX_COUNT || L->nclust == 0u || L->nclust > HZ_KIT_MAX_COUNT)
    return HZ_KIT_E_COUNT;
  uint64_t nb = level_bytes(L, f64);
  if (off > fsize || nb > fsize - off) return HZ_KIT_E_RANGE;
  /* Аллокации инлайном, без обёрток (CLAUDE.md гейт 2); при любом отказе
   * уровень чистит сам — вызывающий не знает о частичных аллокациях.
   * В памяти вершины ВСЕГДА double; флаг файла задаёт точность хранения. */
  L->vx = malloc((size_t)L->nverts * sizeof(double));
  L->vy = malloc((size_t)L->nverts * sizeof(double));
  L->vz = malloc((size_t)L->nverts * sizeof(double));
  L->ti0 = malloc(4ull * L->ntris);
  L->ti1 = malloc(4ull * L->ntris);
  L->ti2 = malloc(4ull * L->ntris);
  L->tcl = malloc(4ull * L->ntris);
  L->tmtl = malloc(4ull * L->ntris);
  L->cl = malloc(40ull * L->nclust);
  if (!L->vx || !L->vy || !L->vz || !L->ti0 || !L->ti1 || !L->ti2 || !L->tcl || !L->tmtl ||
      !L->cl) {
    level_clear(L);
    return HZ_KIT_E_MEM;
  }
  if (fseek(f, (long)off, SEEK_SET) != 0) {
    level_clear(L);
    return HZ_KIT_E_IO;
  }
  {
    int bad = 0;
    for (int c = 0; c < 3 && !bad; c++) {
      double *dst = c == 0 ? L->vx : (c == 1 ? L->vy : L->vz);
      if (f64) {
        bad = rd(f, dst, (size_t)L->nverts * sizeof(double));
      } else {
        float *tmp = malloc((size_t)L->nverts * sizeof *tmp);
        if (tmp == NULL) {
          level_clear(L);
          return HZ_KIT_E_MEM;
        }
        bad = rd(f, tmp, (size_t)L->nverts * sizeof(float));
        for (uint32_t v = 0; v < L->nverts && !bad; v++)
          dst[v] = (double)tmp[v];
        free(tmp);
      }
    }
    if (bad || rd(f, L->ti0, 4ull * L->ntris) || rd(f, L->ti1, 4ull * L->ntris) ||
        rd(f, L->ti2, 4ull * L->ntris) || rd(f, L->tcl, 4ull * L->ntris) ||
        rd(f, L->tmtl, 4ull * L->ntris)) {
      level_clear(L);
      return HZ_KIT_E_IO;
    }
  }
  for (uint32_t c = 0; c < L->nclust; c++) {
    hz_cluster *g = &L->cl[c];
    unsigned char raw[40];
    if (rd(f, raw, 40)) {
      level_clear(L);
      return HZ_KIT_E_IO;
    }
    memcpy(g->bmin, raw, 12);
    memcpy(g->bmax, raw + 12, 12);
    memcpy(&g->err, raw + 24, 4);
    memcpy(&g->first_tri, raw + 28, 4);
    memcpy(&g->ntris, raw + 32, 4);
  }
  int rc = level_validate(L, nmtl);
  if (rc != HZ_KIT_OK) level_clear(L);
  return rc;
}

int hz_kit_load(hz_kit *k, FILE *f) {
  return hz_kit_load_part(k, f, 0);
}

int hz_kit_load_part(hz_kit *k, FILE *f, int32_t min_lev) {
  if (k == NULL || f == NULL) return HZ_KIT_E_ARG;
  hz_kit_init(k);
  if (min_lev < 0) return HZ_KIT_E_ARG;
  unsigned char h[HZ_KIT_HDRSIZE];
  if (rd(f, h, HZ_KIT_HDRSIZE)) return HZ_KIT_E_IO;
  int32_t nlev;
  uint32_t flags;
  int rc = hz_kit_hdr_decode(h, &nlev, &flags);
  if (rc) return rc;
  k->flags = flags;
  uint32_t nmtl = 0;
  memcpy(&nmtl, h + HZ_KIT_OFF_NMTL, 4);
  if (nmtl > HZ_KIT_MAX_COUNT) return HZ_KIT_E_COUNT;
  k->nmtl = nmtl;
  if (fseek(f, 0, SEEK_END) != 0) return HZ_KIT_E_IO;
  long fsz = ftell(f);
  if (fsz < 0) return HZ_KIT_E_IO;
  uint64_t fsize = (uint64_t)fsz;
  uint64_t need = (uint64_t)HZ_KIT_HDRSIZE + (uint64_t)HZ_KIT_DIRSIZE * (uint64_t)nlev;
  if (need > fsize) return HZ_KIT_E_RANGE;
  k->lev = calloc((size_t)nlev, sizeof(hz_kit_level));
  if (k->lev == NULL) return HZ_KIT_E_MEM;
  k->present = calloc((size_t)nlev, 1);
  if (k->present == NULL) {
    free(k->lev);
    k->lev = NULL;
    return HZ_KIT_E_MEM;
  }
  k->nlev = nlev;
  if (nmtl > 0) {
    k->mtl = malloc((size_t)nmtl * sizeof(hz_kit_mtl));
    if (k->mtl == NULL) {
      hz_kit_free(k);
      return HZ_KIT_E_MEM;
    }
    uint64_t moff =
        ((uint64_t)HZ_KIT_HDRSIZE + (uint64_t)HZ_KIT_DIRSIZE * (uint64_t)nlev + HZ_KIT_ALIGN - 1) /
        HZ_KIT_ALIGN * HZ_KIT_ALIGN;
    if (moff + (uint64_t)sizeof(hz_kit_mtl) * (uint64_t)nmtl > fsize) {
      hz_kit_free(k);
      return HZ_KIT_E_RANGE;
    }
    if (fseek(f, (long)moff, SEEK_SET) != 0) {
      hz_kit_free(k);
      return HZ_KIT_E_IO;
    }
    for (uint32_t mi = 0; mi < nmtl; mi++) {
      hz_kit_mtl *mm = &k->mtl[mi];
      if (rd(f, &mm->kd, 8) || rd(f, mm->kd3, 24) || rd(f, mm->ks3, 24) || rd(f, mm->ke3, 24)) {
        hz_kit_free(k);
        return HZ_KIT_E_IO;
      }
    }
  }
  if (fseek(f, HZ_KIT_HDRSIZE, SEEK_SET) != 0) return HZ_KIT_E_IO;
  for (int32_t i = 0; i < nlev; i++) {
    if (i < min_lev) continue; /* страницы детальных уровней НЕ читаются */
    /* Записи оглавления читаются по СВОИМ офсетам: после чтения данных
     * предыдущего уровня поток стоит в конце данных, а не в оглавлении. */
    if (fseek(f, (long)(HZ_KIT_HDRSIZE + HZ_KIT_DIRSIZE * (uint64_t)i), SEEK_SET) != 0) {
      hz_kit_free(k);
      return HZ_KIT_E_IO;
    }
    rc = level_load(&k->lev[i], f, fsize, (flags & HZ_KIT_FLAG_F64) != 0, nmtl);
    if (rc) {
      hz_kit_free(k); /* present[i] ещё 0 — освобождается только загруженное */
      return rc;
    }
    k->present[i] = 1;
  }
  return hz_kit_validate(k);
}
