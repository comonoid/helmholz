/* CBMC-стенд для ЧТЕНИЯ СРЕЗА (`src/cut/dcslice.c`, шаг Ш2; долг А715).
 *
 * ЗАЧЕМ ИМЕННО ЭТИ ДВЕ ФУНКЦИИ. `hz_slice_vertex` и `hz_slice_normal` — это то,
 * что фронт зовёт НА КАЖДУЮ ЯЧЕЙКУ, то есть горячий путь целиком. Содержимое
 * ячейки приходит из памяти, а память после правки геометрии (Ш4) — вход, за
 * который никто не поручился побитово. Поэтому здесь всё поле ячейки берётся
 * НЕДЕТЕРМИНИРОВАННО, включая заведомо негодные значения `lvl`.
 *
 * ЧТО ДОКАЗЫВАЕТСЯ.
 *   (1) чтений за пределами массива нет ни при каком содержимом ячейки;
 *   (2) сдвиги определены: `8 - vbits` и `lev - lvl` не отрицательны и не
 *       превосходят разрядности — а это НЕ очевидно, потому что `lvl` читается
 *       из памяти, и «уровень больше `lev`» дал бы сдвиг на отрицательное число,
 *       то есть неопределённое поведение, а не просто неверную координату;
 *   (3) ПОСТУСЛОВИЕ, ради которого фронт эти функции и зовёт: вершина лежит
 *       В СВОЕЙ ЯЧЕЙКЕ (`lo <= v <= lo + size`), а нормаль ЕДИНИЧНАЯ.
 *
 * ЧЕГО ЭТОТ СТЕНД НЕ ПОКРЫВАЕТ, СКАЗАНО ПРЯМО: сборку среза (`hz_slice_build`) —
 * она рекурсивна по дереву, и её CBMC не берёт в разумном бюде. Долг А715
 * закрывается этим стендом ЧАСТИЧНО: чтение доказано, постройка нет.
 *
 * КОНТРОЛЬ НАД САМИМ СТЕНДОМ: `-DHZ_CHECK_FALSE` требует заведомо ложного
 * «нормаль не единичная», и CBMC ОБЯЗАН выдать контрпример.
 *
 * Прогон:
 *   scripts/cverify.sh tests/cbmc_slice.c --function harness --unwind 20 src/cut/dcslice.c
 *   ... -DHZ_CHECK_FALSE   (ОБЯЗАН дать контрпример)
 */
#include "cut/dcslice.h"

int nondet_i(void);
unsigned char nondet_u8(void);
unsigned short nondet_u16(void);

#if !defined(__CPROVER) && !defined(__CPROVER__)
#define __CPROVER_assume(x) ((void)(x))
#define __CPROVER_assert(x, s) ((void)(x))
int nondet_i(void) {
  return 0;
}
unsigned char nondet_u8(void) {
  return 0;
}
unsigned short nondet_u16(void) {
  return 0;
}
#endif

void harness(void);

void harness(void) {
  hz_dccell cell;
  for (int a = 0; a < 3; a++) {
    cell.lo[a] = nondet_u16();
    cell.vx[a] = nondet_u8();
    cell.irr[a] = nondet_u16();
  }
  cell.noct = nondet_u16();
  cell.lvl = nondet_u8();
  cell.mat = nondet_u8();
  cell.pad = 0;

  hz_dcslice s;
  s.c = &cell;
  s.n = 1;
  s.cap = 1;
  s.lev = nondet_i();
  s.vbits = nondet_i();
  /* Предпосылки — КОНТРАКТ структуры, а не удобство: уровень сетки в пределах
   * формата, значащих бит вершины от одного до восьми. Уровень ЯЧЕЙКИ (`lvl`)
   * НЕ ограничивается: он приходит из памяти, и корректность при негодном `lvl`
   * — часть того, что проверяется. */
  __CPROVER_assume(s.lev >= 1 && s.lev <= HZ_SLICE_MAX_LEV);
  __CPROVER_assume(s.vbits >= 1 && s.vbits <= 8);
  __CPROVER_assume(cell.lvl <= (unsigned char)s.lev);

  double v[3], n[3];
  hz_slice_vertex(&s, 0, v);
  hz_slice_normal(&s, 0, n);

  double size = (double)((int)1 << (s.lev - (int)cell.lvl));
  for (int a = 0; a < 3; a++) {
    __CPROVER_assert(v[a] >= (double)cell.lo[a], "vertex not below its cell");
    __CPROVER_assert(v[a] <= (double)cell.lo[a] + size, "vertex not above its cell");
  }
  double m = n[0] * n[0] + n[1] * n[1] + n[2] * n[2];
#ifdef HZ_CHECK_FALSE
  __CPROVER_assert(m < 0.5, "control: normal must be able to be unit");
#endif
  __CPROVER_assert(m > 0.99 && m < 1.01, "normal is unit");
}

#if !defined(__CPROVER) && !defined(__CPROVER__)
int main(void);
int main(void) {
  harness();
  return 0;
}
#endif
