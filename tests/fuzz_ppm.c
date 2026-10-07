/* fuzz_ppm.c — фаззер НЕДОВЕРЕННОГО ВХОДА: читатель PPM `P6` (Ш8 §575).
 *
 * Читатель берёт ПУТЬ, а не буфер, поэтому драйвер пишет вход в файл с
 * фиксированным именем. Итерация от этого дороже, чем у fuzz_obj (файловый
 * ввод), поэтому цель отдельная: смешивать дешёвую и дорогую в одной сборке
 * значит платить за обе по цене дорогой.
 *
 * Возвращённый буфер освобождается ВСЕГДА (free(NULL) легален): иначе утечка
 * драйвера замаскировала бы утечку читателя.
 *
 * Запуск: make fuzz   (FUZZ_TIME=секунды на цель, по умолчанию 60)
 */
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "image.h"

#define HZ_FUZZ_PPM_PATH "build/fuzz/ppm_input"

/* Прототип — иначе gcc -Wmissing-prototypes считает определение сиротой */
int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size);

/* Точку входа зовёт фреймворк libFuzzer, а не наш код — потому «не используется»;
 * пояснение отдельной строкой: cppcheck не терпит текст в самой директиве */
/* cppcheck-suppress unusedFunction */
int LLVMFuzzerTestOneInput(const uint8_t *data, size_t size) {
  FILE *f = fopen(HZ_FUZZ_PPM_PATH, "wb");
  if (f == NULL) return 0; /* каталога нет — молча: это забота цели Makefile */
  if (size > 0) {
    if (fwrite(data, 1, size, f) != size) {
      fclose(f);
      return 0;
    }
  }
  fclose(f);

  unsigned char *rgb = NULL;
  int w = 0, h = 0;
  (void)hz_ppm_read(HZ_FUZZ_PPM_PATH, &rgb, &w, &h);
  free(rgb);
  return 0;
}
