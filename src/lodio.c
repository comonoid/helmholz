/* lodio.c — чтение и запись файла лестницы LOD. Разбор — в `lodio.h`. */
#include "lodio.h"
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define HZ_LODIO_MAGIC "HZLOD\x01\x00\x00"
#define HZ_LODIO_VERSION 3u
/* Образец порядка байт: все восемь байт различны, поэтому любая перестановка
 * видна, а не только обращение. */
#define HZ_LODIO_ENDIAN 0x0102030405060708ull

typedef struct {
  char magic[8];
  uint32_t version;
  uint32_t sz_double; /* свидетель: sizeof(double) */
  uint32_t sz_node;   /* свидетель: sizeof(hz_lodnode) — ловит выравнивание */
  uint32_t pad;       /* до кратности восьми; пишется нулём и проверяется */
  uint64_t endian;    /* образец порядка байт */
  double one;         /* образец самого double; сверяется ПО БИТАМ, а не `==` */
  int32_t nlev, np;   /* слоёв и исходных полигонов */
  int32_t nnd;        /* узлов всего */
  int32_t bands;      /* происхождение: порядок обхода */
  int32_t bycount;    /* происхождение: тактика уровня */
  int32_t simp;       /* базовая разметка прошла hz_edge_simplify */
  int32_t nt;         /* треугольников в базовой разметке */
  int32_t pad2;       /* до кратности восьми */
  uint64_t sum_seg;   /* сумма по `sg->label`: С КАКОЙ сегментацией лестница парна */
  double eps, delta0; /* под какую камеру и с каким нулевым допуском */
  uint64_t sum;       /* FNV-1a 64 по НАГРУЗКЕ (узлы, затем метки) */
  int64_t npair_edge; /* замеры постройки — чтобы доклад воспроизводился */
  int64_t ngate_rej, nident, nlow_rej;
  double area_grow;
} hz_lodhdr;

const char *hz_lodio_str(int err) {
  switch (err) {
  case HZ_LODIO_OK:
    return "ок";
  case HZ_LODIO_ERR_IO:
    return "ввод-вывод (нет файла либо короткое чтение)";
  case HZ_LODIO_ERR_MAGIC:
    return "не файл лестницы (магия)";
  case HZ_LODIO_ERR_VERSION:
    return "другая версия формата";
  case HZ_LODIO_ERR_ABI:
    return "свидетели представления не совпали (чужая платформа)";
  case HZ_LODIO_ERR_SUM:
    return "контрольная сумма нагрузки не сошлась (порча)";
  case HZ_LODIO_ERR_RANGE:
    return "метка либо родитель вне диапазона";
  case HZ_LODIO_ERR_LEVEL:
    return "уровень узла не равен номеру слоя";
  case HZ_LODIO_ERR_NEST:
    return "вложенность нарушена";
  case HZ_LODIO_ERR_MEM:
    return "нет памяти";
  case HZ_LODIO_ERR_SEG:
    return "лестница построена на ДРУГОЙ сегментации (сумма разметки)";
  default:
    return "неизвестный код";
  }
}

/* FNV-1a 64. Детектор ПОРЧИ, не подпись: одиночный битовый переворот меняет
 * сумму заведомо, а против подделки она и не заявлена. */
static uint64_t fnv(const void *p, size_t n, uint64_t h) {
  const unsigned char *b = (const unsigned char *)p;
  for (size_t i = 0; i < n; i++) {
    h ^= (uint64_t)b[i];
    h *= 1099511628211ull;
  }
  return h;
}

int hz_lod_check(const hz_lod *L) {
  if (L->nlev < 1 || L->np < 1 || L->nnd < 1) return HZ_LODIO_ERR_RANGE;
  if (L->lab == NULL || L->nd == NULL) return HZ_LODIO_ERR_RANGE;
  /* 1. МЕТКИ И РОДИТЕЛИ В ДИАПАЗОНЕ; УРОВЕНЬ УЗЛА РАВЕН НОМЕРУ СЛОЯ. */
  for (int32_t lev = 0; lev < L->nlev; lev++) {
    const int32_t *lab = L->lab + (size_t)lev * (size_t)L->np;
    for (int32_t k = 0; k < L->np; k++) {
      int32_t id = lab[k];
      if (id < 0 || id >= L->nnd) return HZ_LODIO_ERR_RANGE;
      if (L->nd[id].level != lev) return HZ_LODIO_ERR_LEVEL;
      int32_t par = L->nd[id].parent;
      if (par != -1 && (par < 0 || par >= L->nnd)) return HZ_LODIO_ERR_RANGE;
    }
  }
  /* 2. ВЛОЖЕННОСТЬ ПО МЕТКАМ, ЦЕЛОЧИСЛЕННО (А132). Все полигоны, лежащие в одном
   * узле уровня `L`, обязаны лежать в ОДНОМ узле уровня `L+1`. Проверяется через
   * таблицу «узел -> первый увиденный родитель», то есть за один проход и без
   * сортировки. */
  int32_t *seen = malloc((size_t)L->nnd * sizeof *seen);
  if (seen == NULL) return HZ_LODIO_ERR_MEM;
  int rc = HZ_LODIO_OK;
  for (int32_t lev = 1; lev < L->nlev && rc == HZ_LODIO_OK; lev++) {
    const int32_t *lo = L->lab + (size_t)(lev - 1) * (size_t)L->np;
    const int32_t *hi = L->lab + (size_t)lev * (size_t)L->np;
    for (int32_t i = 0; i < L->nnd; i++)
      seen[i] = -1;
    for (int32_t k = 0; k < L->np; k++) {
      int32_t a = lo[k], b = hi[k];
      if (seen[a] == -1)
        seen[a] = b;
      else if (seen[a] != b) {
        rc = HZ_LODIO_ERR_NEST;
        break;
      }
    }
  }
  free(seen);
  return rc;
}

int hz_lod_write(const hz_lod *L, const char *path, const hz_pseglist *sg, int32_t nt, int simp) {
  int rc = hz_lod_check(L);
  if (rc != HZ_LODIO_OK) return rc; /* негодную лестницу на диск не пускаем */
  hz_lodhdr h;
  memset(&h, 0, sizeof h);
  memcpy(h.magic, HZ_LODIO_MAGIC, 8);
  h.version = HZ_LODIO_VERSION;
  h.sz_double = (uint32_t)sizeof(double);
  h.sz_node = (uint32_t)sizeof(hz_lodnode);
  h.endian = HZ_LODIO_ENDIAN;
  h.one = 1.0;
  h.nlev = L->nlev;
  h.np = L->np;
  h.nnd = L->nnd;
  h.bands = L->bands;
  h.bycount = L->bycount;
  h.eps = L->eps;
  h.delta0 = L->delta0;
  h.npair_edge = L->npair_edge;
  h.ngate_rej = L->ngate_rej;
  h.nident = L->nident;
  h.nlow_rej = L->nlow_rej;
  h.area_grow = L->area_grow;
  h.simp = simp ? 1 : 0;
  h.nt = nt;
  h.sum_seg = fnv(sg->label, (size_t)nt * sizeof *sg->label, 1469598103934665603ull);
  size_t nnode = (size_t)L->nnd * sizeof *L->nd;
  size_t nlabb = (size_t)L->nlev * (size_t)L->np * sizeof *L->lab;
  h.sum = fnv(L->lab, nlabb, fnv(L->nd, nnode, 1469598103934665603ull));
  FILE *f = fopen(path, "wb");
  if (f == NULL) return HZ_LODIO_ERR_IO;
  int ok = (fwrite(&h, sizeof h, 1, f) == 1) && (fwrite(L->nd, nnode, 1, f) == 1) &&
           (fwrite(L->lab, nlabb, 1, f) == 1);
  if (fclose(f) != 0) ok = 0;
  return ok ? HZ_LODIO_OK : HZ_LODIO_ERR_IO;
}

int hz_lod_read(hz_lod *L, const char *path, const hz_pseglist *sg, int32_t nt, int *simp) {
  memset(L, 0, sizeof *L);
  FILE *f = fopen(path, "rb");
  if (f == NULL) return HZ_LODIO_ERR_IO;
  hz_lodhdr h;
  if (fread(&h, sizeof h, 1, f) != 1) {
    fclose(f);
    return HZ_LODIO_ERR_IO;
  }
  if (memcmp(h.magic, HZ_LODIO_MAGIC, 8) != 0) {
    fclose(f);
    return HZ_LODIO_ERR_MAGIC;
  }
  if (h.version != HZ_LODIO_VERSION) {
    fclose(f);
    return HZ_LODIO_ERR_VERSION;
  }
  /* Свидетели ДО всякого выделения памяти: чужой файл не должен успеть
   * попросить гигабайт по своему `nnd`. Образец `double` сверяется ПО БИТАМ:
   * `==` на числах с плавающей точкой здесь не только вызывает предупреждение
   * гейта, но и слабее — оно приняло бы иное представление, дающее то же
   * значение, а свидетель заведён именно про ПРЕДСТАВЛЕНИЕ. */
  const double one = 1.0;
  if (h.sz_double != (uint32_t)sizeof(double) || h.sz_node != (uint32_t)sizeof(hz_lodnode) ||
      h.endian != HZ_LODIO_ENDIAN || memcmp(&h.one, &one, sizeof one) != 0 || h.pad != 0 ||
      h.pad2 != 0) {
    fclose(f);
    return HZ_LODIO_ERR_ABI;
  }
  if (h.nlev < 1 || h.np < 1 || h.nnd < 1) {
    fclose(f);
    return HZ_LODIO_ERR_RANGE;
  }
  /* Произведение слоёв на полигоны считается в 64 битах и сверяется с пределом
   * `int32` — иначе испорченный заголовок даёт переполнение размера, а не отказ. */
  int64_t nlab = (int64_t)h.nlev * (int64_t)h.np;
  if (nlab > 2147483000ll) {
    fclose(f);
    return HZ_LODIO_ERR_RANGE;
  }
  /* ВЕРХНЯЯ ГРАНИЦА РАЗМЕРА БЕРЁТСЯ ИЗ ДЛИНЫ ФАЙЛА, А НЕ ИЗ ВЫДУМАННОГО ПРЕДЕЛА
   * (найдено гейтом: `nnd` из файла — величина, подконтрольная тому, кто файл
   * подсунул, CWE-789). Нагрузка ОБЯЗАНА совпасть с объявленной ПОБАЙТОВО: тогда
   * просить память можно ровно под то, что на диске и правда лежит, а обрезание и
   * приписка ловятся ДО всякого выделения. Порога здесь нет — есть тождество. */
  if (fseek(f, 0, SEEK_END) != 0) {
    fclose(f);
    return HZ_LODIO_ERR_IO;
  }
  long fend = ftell(f);
  if (fend < 0 || fseek(f, (long)sizeof h, SEEK_SET) != 0) {
    fclose(f);
    return HZ_LODIO_ERR_IO;
  }
  /* Произведения считаются в 64 битах, и каждый множитель уже ограничен выше:
   * `nlab` проверен, `nnd` есть `int32`, значит `nnd * 104 < 2.3e11` — переполнения
   * нет ни в одном слагаемом. */
  int64_t body = (int64_t)fend - (int64_t)sizeof h;
  if (body < 0) {
    fclose(f);
    return HZ_LODIO_ERR_IO;
  }
  /* ЯВНАЯ ВЕРХНЯЯ ГРАНИЦА НА КАЖДЫЙ ТАИНТЕД-МНОЖИТЕЛЬ, а не только тождество
   * размеров: тождества `gcc -fanalyzer` не отслеживает через арифметику (это его
   * известный класс, CLAUDE.md), и без прямого сравнения он справедливо считает
   * `nnd` неограниченным. Граница снова не выдумана — она есть длина файла,
   * делённая на размер записи. */
  if ((int64_t)h.nnd > body / (int64_t)sizeof(hz_lodnode)) {
    fclose(f);
    return HZ_LODIO_ERR_RANGE;
  }
  if (nlab > body / 4) {
    fclose(f);
    return HZ_LODIO_ERR_RANGE;
  }
  /* И тождество: нагрузка совпадает с объявленной ПОБАЙТОВО. Обрезание и приписка
   * ловятся ДО выделения памяти. */
  int64_t need = (int64_t)h.nnd * (int64_t)sizeof(hz_lodnode) + nlab * 4;
  if (need != body) {
    fclose(f);
    return HZ_LODIO_ERR_IO;
  }
  hz_lodnode *nd = malloc((size_t)h.nnd * sizeof *nd);
  int32_t *lab = malloc((size_t)nlab * sizeof *lab);
  if (nd == NULL || lab == NULL) {
    free(nd);
    free(lab);
    fclose(f);
    return HZ_LODIO_ERR_MEM;
  }
  size_t nnode = (size_t)h.nnd * sizeof *nd;
  size_t nlabb = (size_t)nlab * sizeof *lab;
  int ok = (fread(nd, nnode, 1, f) == 1) && (fread(lab, nlabb, 1, f) == 1);
  fclose(f);
  if (!ok) {
    free(nd);
    free(lab);
    return HZ_LODIO_ERR_IO;
  }
  if (fnv(lab, nlabb, fnv(nd, nnode, 1469598103934665603ull)) != h.sum) {
    free(nd);
    free(lab);
    return HZ_LODIO_ERR_SUM;
  }
  L->nd = nd;
  L->lab = lab;
  L->nnd = h.nnd;
  L->ndcap = h.nnd;
  L->nlev = h.nlev;
  L->np = h.np;
  L->eps = h.eps;
  L->delta0 = h.delta0;
  L->bands = h.bands;
  L->bycount = h.bycount;
  L->npair_edge = h.npair_edge;
  L->ngate_rej = h.ngate_rej;
  L->nident = h.nident;
  L->nlow_rej = h.nlow_rej;
  L->area_grow = h.area_grow;
  if (simp != NULL) *simp = h.simp;
  /* ПАРНОСТЬ С СЕГМЕНТАЦИЕЙ — ПРОВЕРКА, А НЕ ДОГОВОРЁННОСТЬ. Сверяется и число
   * треугольников, и сумма по меткам: одного числа мало, оно совпадает у любых
   * двух разметок одной сцены. */
  if (sg != NULL) {
    if (nt != h.nt ||
        fnv(sg->label, (size_t)nt * sizeof *sg->label, 1469598103934665603ull) != h.sum_seg) {
      hz_lod_free(L);
      return HZ_LODIO_ERR_SEG;
    }
  }
  int rc = hz_lod_check(L);
  if (rc != HZ_LODIO_OK) {
    hz_lod_free(L); /* половина лестницы потребителю не достаётся */
    return rc;
  }
  return HZ_LODIO_OK;
}
