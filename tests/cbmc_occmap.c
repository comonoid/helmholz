/* CBMC-стенд для РАЗБОРЩИКА КАРТЫ ЗАНЯТОСТИ (`tools/occmap.h`, шаг Ш0, А686).
 *
 * ЗАЧЕМ. §383 называет CBMC только в Ш2 и Ш4, но CLAUDE.md требует формальный
 * слой на РАЗБОРЩИКАХ и на всём, что читает недоверенный вход, — а `hz_occ_parse`
 * читает файл с диска и отдаёт вызывающему число ячеек, по которому тот потом
 * ИНДЕКСИРУЕТ битовую карту. Если разбор пропустит несогласованный заголовок
 * (`n`, не равное `2^lev`; длина, не отвечающая `n³`; `count` больше числа
 * ячеек), сверка полезет за пределы буфера — и выглядеть это будет как
 * расхождение МНОЖЕСТВ, то есть подменит величину, ради которой шаг делается.
 *
 * ЧТО ДОКАЗЫВАЕТСЯ.
 *   (1) сам разбор безопасен при ЛЮБОМ содержимом буфера длины `nbuf`: чтение
 *       заголовка происходит только после проверки `nbuf >= sizeof`, и никаких
 *       записей, кроме `*out`, нет;
 *   (2) ПОСТУСЛОВИЕ, на которое опирается вызывающий: при `rc == 0` длина
 *       полезной части РОВНО `(n³+7)/8`, поэтому `hz_occ_get` по любому индексу
 *       `< n³` лежит в границах. Это проверяется не рассуждением, а обращением
 *       по НЕДЕТЕРМИНИРОВАННОМУ индексу под `--bounds-check`.
 *
 * ПОЧЕМУ `lev = 1`. Буфер в стенде обязан быть КОНСТАНТНОГО размера, а полезная
 * часть растёт как `8^lev / 8`: уже при `lev = 4` это `512` байт, и путь
 * «заголовок годен» стал бы недостижим внутри малого буфера, то есть стенд
 * доказывал бы пустое. При `lev = 1` годная длина равна `48 + 1 = 49`, и обе
 * ветви — и отказ, и приём — достижимы. Ветвление по `lev` в самом разборщике
 * при этом покрыто: `lev` берётся недетерминированно, а вход сравнивается с ним.
 *
 * КОНТРОЛЬ НАД САМИМ СТЕНДОМ (образец А431 и `cbmc_pclip.c`). Ключ
 * `-DHZ_CHECK_FALSE` требует заведомо ложного «разбор всегда отказывает», и CBMC
 * ОБЯЗАН выдать контрпример. Без него «SUCCESS» неотличим от стенда, чьи
 * утверждения не дошли до решателя.
 *
 * Прогон:
 *   scripts/cverify.sh tests/cbmc_occmap.c --function harness --unwind 64
 *   ... -DHZ_CHECK_FALSE   (ОБЯЗАН дать контрпример)
 */
#include "occmap.h"

int nondet_i(void);
unsigned char nondet_u8(void);
unsigned long nondet_ul(void);

/* СТОРОЖ ПРОВЕРЯЕТ ОБА ИМЕНИ: CBMC 6.9.0 определяет `__CPROVER__`, а не
 * `__CPROVER` (проверено прямым прогоном `#error` — см. `cbmc_pclip.c`). Под
 * одним только `__CPROVER` предпосылки подменялись бы пустышками, а
 * «недетерминированный» вход оказался бы константой. */
#if !defined(__CPROVER) && !defined(__CPROVER__)
/* чтобы гейт качества (gcc/clang-tidy/cppcheck) разбирал и этот файл */
#define __CPROVER_assume(x) ((void)(x))
#define __CPROVER_assert(x, s) ((void)(x))
int nondet_i(void) {
  return 0;
}
unsigned char nondet_u8(void) {
  return 0;
}
unsigned long nondet_ul(void) {
  return 0;
}
#endif

/* 48 байт заголовка плюс 8 байт полезной части — с запасом на `lev = 1`, где
 * годной является длина 49. Больше не нужно и вредно: размер буфера прямо
 * задаёт стоимость перебора. */
#define HZ_OCCBUF 56

void harness(void);

void harness(void) {
  unsigned char buf[HZ_OCCBUF];
  for (int i = 0; i < HZ_OCCBUF; i++)
    buf[i] = nondet_u8();

  unsigned long nb = nondet_ul();
  __CPROVER_assume(nb <= (unsigned long)HZ_OCCBUF);
  int lev = nondet_i();
  __CPROVER_assume(lev == 1);

  hz_occ_hdr h;
  int rc = hz_occ_parse(buf, (size_t)nb, lev, &h);

#ifdef HZ_CHECK_FALSE
  /* КОНТРОЛЬ НАД СТЕНДОМ: заведомо ложно — приём достижим. Обязан упасть. */
  __CPROVER_assert(rc != 0, "control: parse must be able to accept");
#endif

  if (rc != 0) return;

  /* Постусловие, на которое опирается вызывающий. */
  __CPROVER_assert(h.n == 2, "lev=1 => n=2");
  __CPROVER_assert(h.count >= 0 && h.count <= 8, "count in [0, n^3]");
  __CPROVER_assert(nb == sizeof(hz_occ_hdr) + 1u, "payload length is exactly (n^3+7)/8");

  /* ГЛАВНОЕ: индексация полезной части по ЛЮБОМУ законному индексу ячейки не
   * выходит за буфер. Проверяется обращением, а не рассуждением. */
  unsigned long ci = nondet_ul();
  __CPROVER_assume(ci < 8u);
  int bit = hz_occ_get(buf + sizeof(hz_occ_hdr), (size_t)ci);
  __CPROVER_assert(bit == 0 || bit == 1, "bit is 0 or 1");
}

#if !defined(__CPROVER) && !defined(__CPROVER__)
int main(void);
int main(void) {
  harness();
  return 0;
}
#endif
