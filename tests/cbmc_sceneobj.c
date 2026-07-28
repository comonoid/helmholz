/* Оснастка CBMC для разборщика OBJ (PLAN_ELEMENTS.md, Ш2: «Гейт качества
 * целиком, включая CBMC — разбор недоверенного входа»).
 *
 * ПОЧЕМУ ОСНАСТКА, А НЕ `cverify.sh src/scene_obj.c`. У CBMC нет тела ни для
 * `fread`, ни для `FILE*`, поэтому «CBMC на hz_obj_load» недостижимо как
 * задача. Разбор поэтому разделён (`hz_obj_parse` берёт буфер), а здесь
 * буфер — НЕДЕТЕРМИНИРОВАННЫЙ.
 *
 * ФАЙЛ ВКЛЮЧАЕТСЯ ЦЕЛИКОМ, потому что опасное — статическое. Лексика
 * (`skip_sp`, `skip_line`, `skip_token`, `count_tokens`) и разбор индексов
 * (`parse_int`, `parse_fvert`) — это и есть хождение по недоверенному буферу;
 * они `static`, и снаружи их не позвать.
 *
 * ЧТО ДОКАЗЫВАЕТСЯ, И ПОЧЕМУ ИМЕННО ЭТО.
 *
 *   1. ЛЕКСИКА НЕ ПЕРЕШАГИВАЕТ NUL. Все четыре функции идут по буферу без
 *      знания его длины — конца как отдельного условия в них нет вовсе, его
 *      роль играет `'\0'`. Если хоть одна пропускает терминатор, `--bounds-check`
 *      это ловит, и ловит на ЛЮБОМ содержимом, а не на том, что я придумал.
 *
 *   2. `parse_fvert` ПРИ `ok` ОТДАЁТ ИНДЕКС В ПРЕДЕЛАХ. Это несущее
 *      утверждение всего разбора: `iv` идёт прямо в `m->f`, а оттуда в
 *      `m->v[3·iv]`. Отрицательные (относительные) индексы OBJ делают проверку
 *      неочевидной — `-1` значит «последняя виденная», и арифметика идёт в
 *      `long`, чтобы `nv_seen + a` не переполнилось до сравнения. Здесь это
 *      утверждение, а не комментарий.
 *
 *   3. `count_tokens` НЕ ОТРИЦАТЕЛЕН. От него зависит `ntri += k - 2`, то есть
 *      РАЗМЕР выделения. Отрицательное `k` дало бы отрицательное `ntri`.
 *
 * РАЗМЕР ЗАДАЧИ — ЛЕСТНИЦА, КАК В `cbmc_octree.c`: `HZ_BUFN` компилируемая
 * константа, каждая ступень — отдельный дешёвый прогон.
 *
 * Запуск (лестница, выполненная в докладе Ш2):
 *   scripts/cverify.sh tests/cbmc_sceneobj.c --function harness \
 *       --unwind 14 --unwinding-assertions -DHZ_BUFN=12
 */

/* NOLINTBEGIN */
#include "../src/scene_obj.c"
/* NOLINTEND */
#include <assert.h>

char nondet_char(void);
int32_t nondet_i32(void);

#ifndef __CPROVER
/* чтобы гейт качества (gcc/clang-tidy/cppcheck) тоже разбирал этот файл; под
 * CBMC имя встроенное и это объявление не компилируется */
void __CPROVER_assume(int);
#endif

#ifndef HZ_BUFN
#define HZ_BUFN 12
#endif

void harness(void);

void harness(void) {
  char buf[HZ_BUFN];
  for (int i = 0; i < HZ_BUFN - 1; i++)
    buf[i] = nondet_char();
  /* ЕДИНСТВЕННОЕ ОГРАНИЧЕНИЕ ОСНАСТКИ, и оно есть договор с `slurp`, а не
   * удобство: `slurp` всегда дописывает `'\0'`, поэтому буфер без терминатора
   * разборщику прийти не может. Всё остальное содержимое произвольно. */
  buf[HZ_BUFN - 1] = '\0';

  /* 1. лексика не перешагивает NUL: адрес обязан остаться в буфере */
  const char *a = skip_sp(buf);
  assert(a >= buf && a <= buf + HZ_BUFN - 1);
  const char *b = skip_token(buf);
  assert(b >= buf && b <= buf + HZ_BUFN - 1);
  const char *c = skip_line(buf);
  assert(c >= buf && c <= buf + HZ_BUFN);

  /* 3. счётчик слов неотрицателен — от него зависит РАЗМЕР выделения */
  int k = count_tokens(buf);
  assert(k >= 0);

#ifdef HZ_CBMC_NEGCTL
  /* НЕГАТИВНЫЙ КОНТРОЛЬ САМОЙ ОСНАСТКИ, с предсказанным провалом. Буфер из
   * одних `'\0'` даёт `k = 0`, значит `k >= 1` ОБЯЗАНО не доказаться. Прогон,
   * который это проходит, означает, что недетерминированность до буфера не
   * доходит и все SUCCESS выше — пустые. */
  assert(k >= 1);
#endif

  /* 2. индексы вершин в пределах виденного */
  int32_t nv_seen = nondet_i32(), nvn_seen = nondet_i32();
  __CPROVER_assume(nv_seen >= 0 && nv_seen <= 4);
  __CPROVER_assume(nvn_seen >= 0 && nvn_seen <= 4);
  int32_t iv = -7, ivn = -7;
  int ok = 0;
  const char *e = parse_fvert(buf, nv_seen, nvn_seen, &iv, &ivn, &ok);
  assert(e >= buf && e <= buf + HZ_BUFN - 1);
  if (ok) {
    assert(iv >= 0);
    assert(iv < nv_seen);
    assert(ivn == -1 || (ivn >= 0 && ivn < nvn_seen));
  }

  /* parse_int: только выход за буфер; ЗНАЧЕНИЕ не ограничивается — потолок
   * 2e8 внутри есть отказ от битого файла, а не свойство безопасности */
  long v = 0;
  int ok2 = 0;
  const char *g = parse_int(buf, &v, &ok2);
  assert(g >= buf && g <= buf + HZ_BUFN - 1);
}
