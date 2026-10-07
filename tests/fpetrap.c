/* fpetrap.c — ловушка FP-исключений для РАНТАЙМ-слоя гейта (07-10).
 *
 * Зачем. CLAUDE.md, пункт 5 гейта: «численные прогоны — раз с feenableexcept
 * на FE_INVALID|FE_DIVBYZERO, чтобы NaN-источники ловились рано». До сих пор
 * это было только прозой. NaN/деление на ноль в поле обязано падать СРАЗУ, а
 * не расползаться в тёмные куски: у проекта свежий прецедент — расхождение
 * клеточного зонда (Ш10) искали долго и дорого.
 *
 * Как. Это НЕ отдельный прибор и не правка продуктового кода: библиотека
 * подгружается через LD_PRELOAD и в своём конструкторе включает ловушки на
 * процесс. Продуктовые файлы остаются нетронутыми.
 *
 *   make check-fpe
 *
 * Выключить на один прогон (например, чтобы отличить источник): HZ_FPETRAP=0.
 *
 * Почему SIGFPE-обработчик, а не голый feenableexcept: без обработчика процесс
 * умирает молча, и в отчёте гейта не видно, ЧТО случилось. Здесь печатается
 * строка и код возврата 77 — «гейт-сбой», отличимый от обычного падения.
 * Трассировки стека нет намеренно: в конструкторе preload-библиотеки нет
 * надёжного доступа к символам, а частичный backtrace хуже честной строки.
 */
#define _GNU_SOURCE
#include <fenv.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/* `#pragma STDC FENV_ACCESS ON` здесь СОЗНАТЕЛЬНО нет: gcc его игнорирует
 * (-Werror=unknown-pragmas, ловля гейта 07-10), а требование «не переставлять
 * FP-операции» в проекте держится флагами сборки: -ffp-contract=off и запрет
 * -ffast-math (CLAUDE.md, Г31). Прагма ничего не добавляла, кроме сломанного S1. */

static void hz_fpe_handler(int sig) {
  (void)sig;
  /* write(), а не fprintf: обработчик сигнала обязан быть reentrant */
  static const char msg[] =
      ">>> ГЕЙТ-СБОЙ: FP-исключение (FE_INVALID или FE_DIVBYZERO) — см. check-fpe\n";
  if (write(2, msg, sizeof msg - 1) < 0) {
    /* ничего не поделать: сообщение — единственный выход */
  }
  _exit(77);
}

__attribute__((constructor)) static void hz_fpetrap_on(void) {
  const char *e = getenv("HZ_FPETRAP");
  if (e != NULL && e[0] == '0') return; /* выключатель на разбор одного прогона */
  struct sigaction sa;
  sa.sa_handler = hz_fpe_handler;
  sigemptyset(&sa.sa_mask);
  sa.sa_flags = SA_NODEFER; /* исключение не должно теряться в обработчике */
  (void)sigaction(SIGFPE, &sa, NULL);
  (void)feenableexcept(FE_INVALID | FE_DIVBYZERO);
}
