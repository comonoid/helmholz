#ifndef HZ_IMAGE_H
#define HZ_IMAGE_H

/* Minimal image output: binary PPM (P6), zero dependencies. Values are
 * intensities >= 0; tone mapping = normalize to max, gamma 0.5, inferno-ish
 * tiny colormap. */

int hz_ppm_write(const char *path, const double *intensity, int w, int h);

/* КАНОНИЧЕСКИЙ ВЫХОД: РАДИАНС КАК ЕСТЬ, 32-битный float, RGB, PFM.
 *
 * ЗАЧЕМ ОТДЕЛЬНЫЙ ФОРМАТ. `hz_ppm_write` нормирует на 99.5-й процентиль, берёт
 * гамму 0.5 и раскрашивает ЛОЖНЫМ цветом. Все три операции необратимы, и К19
 * измерила, чем это оборачивается: та же ошибка после них доходит до 255
 * уровней из 255, то есть по картинке нельзя сверить ничего. Метрика поэтому
 * берётся на буфере радианса, а не на PPM, — и хранить надо ровно тот буфер.
 *
 * PFM выбран потому, что пишется двадцатью строками без единой библиотеки
 * (проект на C и без зависимостей), несёт полный float без обрезки и читается
 * ImageMagick, GIMP, Krita, OpenCV. Порядок строк в PFM — СНИЗУ ВВЕРХ, это
 * часть формата, а не наша прихоть.
 *
 * ЦВЕТА У НАС ПОКА НЕТ, И ФАЙЛ ОБ ЭТОМ НЕ ВРЁТ. Считается ОДИН спектральный
 * канал, поэтому R = G = B, и записанное есть честный серый радианс. Формат
 * взят трёхканальный (`PF`, а не `Pf`) сознательно: когда каналов станет три,
 * поменяется вызывающий код, а не формат файла и не читатели.
 *
 * scale < 0 в заголовке означает little-endian — записывается порядок ЭТОЙ
 * машины, определяемый на месте, а не предполагаемый. */
int hz_pfm_write(const char *path, const double *r, const double *g, const double *b, int w, int h);

/* signed field (e.g. Re u): blue-white-red, symmetric percentile scaling */
int hz_ppm_write_signed(const char *path, const double *field, int w, int h);

/* raw rgb writer for composed images */
int hz_ppm_write_rgb(const char *path, const unsigned char *rgb, int w, int h);

#endif
