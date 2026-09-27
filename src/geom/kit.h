#ifndef HZ_GEOM_KIT_H
#define HZ_GEOM_KIT_H

/* kit.h — ГЕОМЕТРИЯ v2, Ш1 §914 (STRUCTURE.md §2.2 ревизия v2 от 12-09).
 *
 * КИТ = лестница ЗАМЕЩАЮЩИХ уровней: lev[0] = L0 (исходные треугольники,
 * их никто не переписывает), lev[i>0] — грубее (декимация — Ш2). Уровни
 * не дублируют, а заменяют: ntris монотонно не возрастает вверх (А1582).
 *
 * ГРОЗДЬ — единица выбора уровня и резидентности: связная группа
 * треугольников уровня; грозди покрывают треугольники уровня РОВНО,
 * упорядочены по first_tri (А1583). Граничные вершины грозди декимация
 * не двигает — это уговор Ш2, формат только резервирует ε (err).
 *
 * РАСКЛАДКА ФАЙЛА (§2.2 v2, «уровни отдельно»): [заголовок 64 Б]
 * [оглавление уровней, 32 Б/уровень] [данные: от САМОГО ГРУБОГО к L0],
 * блоки выровнены по 64 Б. Страница 2 МБ — Ш5 (А1584). Чтение грубых
 * уровней не обязано касаться страниц L0 — порядок блоков это даёт.
 *
 * Формат v1, свидетели — класс Г5/Г28 (А1585): magic, версия, byte-order
 * witness, IEEE754-канарейка, unknown flags rejected, reserved = 0.
 * Сериализация SoA: отдельные массивы vx/vy/vz/ti0/ti1/ti2/tcl + грозди;
 * сырой дамп структур не пишется (переносимость ABI, урок Г5).
 */
#include <stdint.h>
#include <stdio.h>

/* Кэпы счётчиков (А1581, урок HZ_OCT_MAX_NODES / CWE-789): атакующий файл
 * не должен уметь спросить гигабайт до первой проверки. 1e9 недостижимо
 * для реальных китов и покрывает "храм" (1e7) с запасом ×100. */
#define HZ_KIT_MAX_NLEV ((int32_t)64)
#define HZ_KIT_MAX_COUNT ((uint32_t)1000000000u)

/* Формат файла v2 (Ш4 §914): + материалы, + флаг F64. v1-файлов в
 * природе нет — версия бампнута, v1 отклоняется. */
#define HZ_KIT_HDRSIZE 64
#define HZ_KIT_DIRSIZE 32
#define HZ_KIT_VERSION 2u
#define HZ_KIT_FLAG_PAD64 1u /* выравнивание блоков уровней по 64 Б */
#define HZ_KIT_FLAG_F64                                                                            \
  2u /* вершины double (8 Б): паритет/архив;                                   \
      * без флага — float (рабочая лестница) */

/* байтовые офсеты заголовка — тесты мутируют файл по ним */
#define HZ_KIT_OFF_MAGIC 0
#define HZ_KIT_OFF_VERSION 8
#define HZ_KIT_OFF_FLAGS 12
#define HZ_KIT_OFF_ENDIAN 16
#define HZ_KIT_OFF_NLEV 20
#define HZ_KIT_OFF_CANARY 24
#define HZ_KIT_OFF_RESERVED 32 /* v2: [32,40) нули; [40,44) nmtl; [44,64) нули */
#define HZ_KIT_OFF_NMTL 40

typedef struct {
  float bmin[3]; /* bbox грозди, метры */
  float bmax[3];
  float err;          /* ошибка уровня ε на грозди (Ш2); L0: 0 */
  uint32_t first_tri; /* диапазон [first_tri, first_tri+ntris) уровня */
  uint32_t ntris;
} hz_cluster; /* 36 Б как в памяти (7×float + 2×u32, без атрибутов) */

/* Сериализация грозди — фиксированная запись 40 Б (7 float + 2 u32 +
 * 4 Б паддинга), от компиляторной упаковки не зависит. */

typedef struct {
  uint32_t nverts, ntris, nclust;
  /* SoA; индексы вершин — uint32 (кит ≤ 1e9 вершин, 24-бит уже мало).
   * В ПАМЯТИ вершины всегда double; флаг F64 файла задаёт ТОЧНОСТЬ
   * хранения (паритет/архив — double, рабочая лестница — float). */
  double *vx, *vy, *vz;      /* [nverts] */
  uint32_t *ti0, *ti1, *ti2; /* [ntris] */
  uint32_t *tcl;             /* [ntris] гроздь треугольника */
  uint32_t *tmtl;            /* [ntris] индекс материала */
  hz_cluster *cl;            /* [nclust], упорядочены по first_tri */
} hz_kit_level;

/* Материал кита (§2.2 v2: «источники — те же киты, Ke в материале»).
 * Поля соответствуют hz_obj_mtl (kd, kd3, ks3, ke3) — конверсия в
 * hz_objmesh побайтовая, паритет pgather по построению. 80 Б
 * (8 + 3×24); СЕРИЛИЗУЕТСЯ ПОЛЯМИ (kd, kd3, ks3, ke3), не дампом. */
typedef struct {
  double kd, kd3[3], ks3[3], ke3[3];
} hz_kit_mtl;

typedef struct {
  int32_t nlev;      /* ≥ 1; lev[0] = L0 (самый детальный) */
  hz_kit_level *lev; /* [nlev]; lev[nlev-1] — самый грубый */
  hz_kit_mtl *mtl;   /* [nmtl]; [0] — умолчание (как hz_obj_mtl) */
  uint32_t nmtl;
  uint32_t flags; /* PAD64 | F64 */
  /* Раздельность уровней В ПАМЯТИ (§914-R2): present == NULL — все
   * уровни в памяти (программное построение); иначе present[i] = 1
   * только для загруженных. Грубые уровни резидентны БЕЗ детальных:
   * hz_kit_load_part читает оглавление + уровни [min_lev, nlev). */
  uint8_t *present; /* [nlev] или NULL (= все) */
} hz_kit;

void hz_kit_init(hz_kit *k);
void hz_kit_free(hz_kit *k);

/* Проверка структур в памяти (те же инварианты, что у загрузчика);
 * отсутствующие уровни (present[i] == 0) пропускаются. Инварианта
 * лестницы проверяется по ПРИСУТСТВУЮЩИМ уровням. */
int hz_kit_validate(const hz_kit *k);

int hz_kit_save(const hz_kit *k, FILE *f);
int hz_kit_load(hz_kit *k, FILE *f); /* в неинициализированный kit; все уровни */

/* Загрузить ТОЛЬКО грубые уровни [min_lev, nlev): pages L0 не читаются
 * вовсе (требование §2.2 v2 + §914-R2). min_lev = 0 эквивалентно
 * hz_kit_load. Возврат — код hz_kit_status. */
int hz_kit_load_part(hz_kit *k, FILE *f, int32_t min_lev);

/* Чистая половина загрузки — декодер заголовка (под CBMC: без FILE*). */
int hz_kit_hdr_decode(const unsigned char h[HZ_KIT_HDRSIZE], int32_t *nlev, uint32_t *flags);

typedef enum {
  HZ_KIT_OK = 0,
  HZ_KIT_E_IO = 1, /* короткое чтение/запись */
  HZ_KIT_E_MAGIC = 2,
  HZ_KIT_E_VERSION = 3,
  HZ_KIT_E_FLAGS = 4, /* неизвестный бит флагов */
  HZ_KIT_E_ENDIAN = 5,
  HZ_KIT_E_CANARY = 6, /* IEEE754 представление отличается */
  HZ_KIT_E_RESERVED = 7,
  HZ_KIT_E_NLEV = 8,    /* nlev вне [1, 64] */
  HZ_KIT_E_COUNT = 9,   /* счётчик уровня вне (0, 1e9] / nverts < 3 */
  HZ_KIT_E_LADDER = 10, /* ntris возрастает вверх (А1582) */
  HZ_KIT_E_TRIDX = 11,  /* индекс вершины ≥ nverts */
  HZ_KIT_E_CLIDX = 12,  /* tcl ≥ nclust */
  HZ_KIT_E_COVER = 13,  /* грозди не покрывают [0,ntris) ровно (А1583) */
  HZ_KIT_E_RANGE = 14,  /* офсет/размер блока вне файла (u64 overflow) */
  HZ_KIT_E_NAN = 15,    /* не-finite в вершинах/bbox/err */
  HZ_KIT_E_MEM = 16,    /* аллокация */
  HZ_KIT_E_ARG = 17,    /* NULL-аргументы в конструкторах */
  HZ_KIT_E_MTLIDX = 18  /* tmtl ≥ nmtl */
} hz_kit_status;

#endif
