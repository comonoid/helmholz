# scripts/cflags.mk — ЕДИНЫЙ источник флагов живого слоя (src, tools, tests).
#
# Зачем: до 07-10 список предупреждений жил в трёх копиях (Makefile, ccheck.sh,
# cverify.sh) и уже разошёлся — `-Wstrict-overflow=2` был только в ccheck.sh,
# поэтому сборка и гейт проверяли РАЗНОЕ. Теперь копия одна: Makefile делает
# `include scripts/cflags.mk`, а bash-скрипты читают её через load_cflags()
# (scripts/cgate.sh).
#
# ПРАВИЛА ФАЙЛА (иначе сломается bash-загрузчик):
#   * только плоские значения: `ИМЯ := значение` или `ИМЯ ?= значение`;
#   * без переносов строк, без ссылок $(...) и без кавычек внутри значений;
#   * комментарии — целой строкой.
# Составные строки (база + оптимизация + предупреждения) собирают потребители.

ARCH ?= -march=skylake

# --- включения и база компиляции -------------------------------------------
HZ_INC := -I src
HZ_BASE := -std=gnu11 -fopenmp -ffp-contract=off
HZ_LIBS := -llapacke -llapack -lblas -lm

# --- предупреждения ---------------------------------------------------------
# Общий набор: годится и gcc, и clang (его же читает clang-tidy через
# build/compile_commands.json).
HZ_WARN := -Wall -Wextra -Wshadow -Wconversion -Wsign-conversion -Wpointer-arith -Wnull-dereference -Wcast-qual -Wwrite-strings -Wvla -Wformat=2 -Wundef -Wstrict-prototypes -Wold-style-definition -Wmissing-prototypes -Wstrict-overflow=2 -Wdouble-promotion -Wfloat-equal
# Только gcc. -Wjump-misses-init — класс «goto перескакивает инициализатор»,
# тот самый, что даёт free() по мусору. Берётся В ОДИНОЧКУ: обычно его
# включают группой -Wc++-compat, но та добавляет 204 предупреждения о
# неявном приведении void* (в C это норма), то есть шум. Проверено 07-10:
# -Wjump-misses-init отдельно даёт ровно 98 находок и ничего лишнего.
HZ_WARN_GCC := -Wjump-misses-init

# --- санитайзеры (рантайм-слой, см. Makefile) -------------------------------
HZ_SAN := -fsanitize=address,undefined -fno-omit-frame-pointer -g
# Целочисленная топология: clang-only (implicit-conversion, unsigned overflow).
HZ_SAN_INT := -fsanitize=integer,implicit-conversion -fno-omit-frame-pointer -g
# Детерминированный мусор в автоматических переменных вместо случайного:
# ловит чтение неинициализированного как повторяемый результат.
HZ_AUTOINIT := -ftrivial-auto-var-init=pattern
