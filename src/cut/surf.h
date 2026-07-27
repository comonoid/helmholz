#ifndef HZ_CUT_SURF_H
#define HZ_CUT_SURF_H

/* ОБЩИЙ СЛОЙ (PLAN_CUT.md, Р-3). Три таблицы, а не одна: у них разные времена
 * жизни и разные владельцы.
 *
 *   ПОВЕРХНОСТИ  аналитические примитивы — ИСТОЧНИК ИСТИНЫ, LOD не подлежат
 *   ФАСЕТЫ       плоскости, которыми режут; у кривого примитива их много
 *   БОКОВАЯ      ячейка -> набор фасетов со СТОРОНОЙ
 *
 * ПОЧЕМУ ССЫЛКА, А НЕ ЗНАЧЕНИЯ (Г3). Соседние ячейки ссылаются на ТОТ ЖЕ
 * объект, поэтому водонепроницаемость точна ПО ПОСТРОЕНИЮ — рёберных
 * пересечений, согласования и щелей не существует как темы. Любой интерполянт
 * по углам грубил бы геометрию до O(h²), что запрет на LOD по геометрии
 * (CLAUDE.md) не разрешает вовсе.
 *
 * ПОЧЕМУ КЛЮЧ — АБСТРАКТНАЯ ЯЧЕЙКА, А НЕ УЗЕЛ ДЕРЕВА (Г20). Волновая линия
 * живёт на hz_octree (ключ = индекс узла, hz_oct_leaf), перенос на tr2_mesh
 * (ключ = индекс ячейки). Таблица со словом «узел» в ключе была бы непригодна
 * переносу — это ровно класс ошибки Г2. Здесь ключ int32 и ничей.
 *
 * ЧЕГО ЗДЕСЬ НЕТ. Материалов: hz_surf несёт две НЕПРОЗРАЧНЫЕ метки сторон, и
 * общий слой в них не смотрит — волновая линия положит туда индекс своего k²,
 * перенос индекс своих σ. Иначе физика одной линии просочилась бы в общий
 * слой (Г2, Г6). */

#include "cut/poly3.h"
#include <stdint.h>

/* --- 1. поверхности: аналитические примитивы ------------------------------ */

typedef enum {
  HZ_SURF_PLANE = 0,  /* p: n[3], off */
  HZ_SURF_SPHERE = 1, /* p: c[3], r   */
  HZ_SURF_CYL = 2     /* p: a[3], d[3], r */
} hz_surf_kind;

typedef struct {
  int32_t kind;
  double p[7];
  int32_t tag_le; /* НЕПРОЗРАЧНАЯ метка стороны {n·x <= off} примитива */
  int32_t tag_ge; /* то же для {n·x >= off} */
} hz_surf;

typedef struct {
  hz_surf *s;
  int32_t n, cap;
} hz_surftab;

int hz_surftab_init(hz_surftab *t);
void hz_surftab_free(hz_surftab *t);
int32_t hz_surftab_add(hz_surftab *t, const hz_surf *s); /* индекс или -1 */

/* --- 2. фасеты: плоскости, которыми режут --------------------------------- */

/* «Смещение не вычислено». НЕ NaN: NaN в поле данных разъедается арифметикой и
 * ссорится с feenableexcept(FE_INVALID) из гейта, а отрицательное смещение не
 * имеет смысла и потому свободно как признак. */
#define HZ_FACET_DMAX_UNKNOWN (-1.0)

typedef struct {
  double n[3], off; /* полупространство n·x <= off, В ЕДИНИЦАХ кадра */
  int32_t surf;     /* примитив-родитель; -1 = плоскость задана прямо */
  double dmax;      /* макс |смещение| фасет<->примитив, МЕТРЫ; <0 = не вычислено */
} hz_facet;

typedef struct {
  hz_facet *f;
  int32_t n, cap;
} hz_facettab;

int hz_facettab_init(hz_facettab *t);
void hz_facettab_free(hz_facettab *t);

/* Плоскость задаётся в МИРОВЫХ координатах и переводится в единицы кадра ЗДЕСЬ,
 * ОДИН РАЗ. В этом весь смысл: две ячейки, сославшиеся на этот фасет, получат
 * ПОБИТОВО одну и ту же плоскость, потому что перевод был один. */
int32_t hz_facettab_add_plane(hz_facettab *t, const hz_frame *fr, const double n[3], double off,
                              int32_t surf, double dmax);

/* Ошибка волновой линии от смещения фасета: измеренный закон M14 [3],
 * 0.125 (k d)² sqrt(d / W). ВОЗВРАЩАЕТ НЕНОЛЬ, ЕСЛИ dmax НЕ ВЫЧИСЛЕН, и *err
 * не трогает (Г25): между пунктом 3, который поле заводит, и пунктом 5, который
 * его заполняет, потребитель обязан ОТКАЗАТЬСЯ, а не подставить ноль. Ноль тут
 * самое опасное значение — он означает «фасет точен». */
int hz_facet_error(const hz_facet *f, double k, double W, double *err);

/* --- 3. боковая таблица: ЯЧЕЙКА -> набор фасетов со стороной --------------- */

/* fref >= 0 — фасет fref, держим {n·x <= off};
 * fref <  0 — фасет ~fref, держим {n·x >= off}.
 * ДОПОЛНЕНИЕ, А НЕ МИНУС: -0 == 0, и нулевой фасет потерял бы сторону. */
typedef struct {
  int32_t cell, f0, nf;
} hz_cutrec;

typedef struct {
  hz_cutrec *r;
  int32_t nr, rcap;
  int32_t *fref;
  int32_t nf, fcap;
} hz_cutmap;

int hz_cutmap_init(hz_cutmap *m);
void hz_cutmap_free(hz_cutmap *m);
/* cell обязан СТРОГО ВОЗРАСТАТЬ: массив держится сортированным, и на этом стоит
 * как двоичный поиск, так и курсор. Нарушение — ошибка, а не пересортировка. */
int hz_cutmap_add(hz_cutmap *m, int32_t cell, const int32_t *fref, int32_t nf);
const hz_cutrec *hz_cutmap_find(const hz_cutmap *m, int32_t cell); /* O(log n) */

/* КУРСОР. Точечного поиска мало: сборка обходит ячейки ПО ВОЗРАСТАНИЮ, и
 * логарифм сел бы во внутренний цикл на всех 1e7…1e8 элементах. Слияние даёт
 * O(1) на ячейку. cell от вызова к вызову не убывает. */
typedef struct {
  const hz_cutmap *m;
  int32_t i;
} hz_cutcur;

void hz_cutcur_begin(hz_cutcur *c, const hz_cutmap *m);
const hz_cutrec *hz_cutcur_next(hz_cutcur *c, int32_t cell);

/* --- мост к ядру ----------------------------------------------------------- */

/* Разворачивает запись в ОРИЕНТИРОВАННЫЕ полупространства, готовые для
 * hz_poly3_cut. hid[j] = сам fref, поэтому fsrc грани в hz_poly3 говорит, каким
 * фасетом и с какой стороны она порождена; hid_flipped[j] = ~fref для
 * hz_poly3_complement. Любой из трёх выходных массивов может быть NULL. */
int hz_cutmap_hspaces(const hz_facettab *ft, const hz_cutmap *m, const hz_cutrec *r, hz_hspace *h,
                      int32_t *hid, int32_t *hid_flipped, int max);

#endif
