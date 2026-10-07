# helmholz — сборка живого слоя (nstruct + общий слой) и ГЕЙТ качества.
# Toolchain — nix (shell.nix, пин инструментов; см. CLAUDE.md).
#
# 12-09: старый геометрический путь (src/p*, src/cut, pfield/pwalk/render3)
# перенесён в archive/geom — исходники §914 создаются заново, старые якоря
# сверяются через git-историю.
# Живая линия: src/nstruct/ (пирамида, свип), pgather/pref (сбор и эталон),
# общий слой: octree, scene_obj, image, transport/dirs3 + quad (test_bounce).
#
# 07-10: флаги вынесены в scripts/cflags.mk (ЕДИНЫЙ источник для сборки, гейта
# и compile_commands.json), гейт — в scripts/cgate.sh (ступени S0…S4), добавлен
# рантайм-слой (санитайзеры, integer-санитайзер, FPE-ловушка, valgrind).

# Флаги: один источник. Значения — scripts/cflags.mk (его же читает bash).
include scripts/cflags.mk

PKGS = gcc lapack blas pkg-config
# Внутри пин-шелла (nix-shell shell.nix) вложенный nix-shell не нужен.
ifeq ($(IN_NIX_SHELL),)
RUN = nix-shell --run
else
RUN =
endif
CFLAGS = $(HZ_BASE) -O2 $(ARCH) $(HZ_WARN) $(HZ_INC)
LIBS   = $(HZ_LIBS)

all: build/test_octree build/pgather build/pref

build:
	mkdir -p build

#-- -- общий слой : октодерево -- --
build/test_octree: tests/test_octree.c src/octree.c src/octree.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_octree.c src/octree.c -lm'

build/test_octfmt: tests/test_octfmt.c src/octree.c src/octree.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_octfmt.c src/octree.c -lm'

#-- -- О75 : нормировка первого отскока(А537) — живой гейт §875 /§912 -- --
build/test_bounce: tests/test_bounce.c src/transport/dirs3.c src/transport/dirs3.h \
                   src/transport/quad.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_bounce.c src/transport/dirs3.c \
	  src/transport/quad.c -lm'

#-- -- НОВАЯ СТРУКТУРА -- --
build/snap3: tools/snap3.c src/nstruct/pyr.c src/nstruct/pyr.h src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/snap3.c src/nstruct/pyr.c src/scene_obj.c -lm'

build/swee3: tools/swee3.c src/nstruct/pyr.c src/nstruct/pyr.h src/nstruct/sweep.c \
	src/nstruct/sweep.h src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/swee3.c src/nstruct/pyr.c src/nstruct/sweep.c src/scene_obj.c -lm'

build/pblock: tools/pblock.c src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pblock.c src/scene_obj.c -lm'

build/pgather: tools/pgather.c src/nstruct/pyr.c src/nstruct/pyr.h src/nstruct/sweep.c \
	src/nstruct/sweep.h src/scene_obj.c src/geom/kit.c src/geom/kit.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pgather.c src/nstruct/pyr.c src/nstruct/sweep.c \
	  src/scene_obj.c src/geom/kit.c -lm'

#-- -- §928: сборка таблицы k(r,ω) из дампов эталона -- --
build/kcalmk: tools/kcalmk.c src/nstruct/sweep.h | build
	$(RUN) 'gcc $(CFLAGS) -I. -o $@ tools/kcalmk.c -lm'

#-- -- §931 Ш1: доли колец детальности D(R) по зрелому E-дампу -- --
build/farshare: tools/farshare.c src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/farshare.c src/scene_obj.c -lm'

build/pref: tools/pref.c src/nstruct/pyr.c src/nstruct/pyr.h src/nstruct/sweep.c \
	src/nstruct/sweep.h src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pref.c src/nstruct/pyr.c src/nstruct/sweep.c src/scene_obj.c -lm'

build/ppmdiff: tools/ppmdiff.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/ppmdiff.c -lm'

build/pfmdiff: tools/pfmdiff.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pfmdiff.c -lm'

build/scenechk: tools/scenechk.c src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/scenechk.c src/scene_obj.c -lm'

build/pview: tools/pview.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pview.c -lm'

#-- -- ГЕОМЕТРИЯ v2(§914, Ш1) -- --
build/test_kit: tests/test_kit.c src/geom/kit.c src/geom/kit.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_kit.c src/geom/kit.c -lm'

# §914 - Ш2': движение фронта по геометрической структуре (кусок с этажом P).
build/kitwalk: tools/kitwalk.c src/nstruct/pyr.c src/nstruct/pyr.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/kitwalk.c src/nstruct/pyr.c -lm'

# §914 - Ш4 : OBJ → вырожденный кит(F64) для проверки паритета pgather.
build/kitmk: tools/kitmk.c src/geom/kit.c src/geom/kit.h src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/kitmk.c src/geom/kit.c src/scene_obj.c -lm'

# §914 - Ш5 / П6 : touched - страницы mmap при чтении грубых уровней.
build/kitpages: tools/kitpages.c src/geom/kit.c src/geom/kit.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/kitpages.c src/geom/kit.c -lm'

#fast tests(seconds)
FAST_TESTS = build/test_octree build/test_octfmt build/test_bounce build/test_kit

test: $(FAST_TESTS)
	./build/test_octree
	./build/test_octfmt
	./build/test_bounce
	./build/test_kit

# ================================ ГЕЙТ ======================================
# Ступени статики (S0…S4) — scripts/cgate.sh; база разобранных находок —
# scripts/cgate-baseline.txt; отчёт — build/cgate/report.tsv.
#
#   make check-staged   S0–S2 по изменённым файлам (для pre-commit: секунды)
#   make check-fast     S0–S2 по всему живому слою (секунды)
#   make check          S0–S4 по всему живому слою (главный статический гейт)
#   make gate           всё: статика + lean + CBMC-стенды + рантайм-слой
#
# Порядок ступеней — по цене, с ранним выходом (см. шапку cgate.sh).

compile-db:
	$(RUN) 'scripts/gen_compile_db.sh'

check-fast:
	scripts/cgate.sh --all --stage=0,1,2

check-staged:
	scripts/cgate.sh --staged --stage=0,1,2

check:
	scripts/cgate.sh --all --stage=0,1,2,3,4

check-deep: check

# Разбор находок: дописать новые в базу (вердикт «НОВОЕ» требует разбора).
gate-baseline:
	scripts/cgate.sh --all --stage=0,1,2,3,4 --update-baseline

# Таблица всех находок — вход для триажа (в т.ч. субагентом).
gate-digest:
	scripts/cgate.sh --all --stage=0,1,2,3,4 --strict --digest

gate-versions:
	scripts/cgate.sh --help >/dev/null
	@nix-shell --run 'gcc -dumpfullversion; clang-tidy --version | sed -n 2p; \
	  cppcheck --version; clang-format --version; cbmc --version; \
	  valgrind --version'

gate: check
	$(MAKE) lean
	$(MAKE) cbmc
	$(MAKE) test-asan
	$(MAKE) test-uninit
	$(MAKE) check-int
	$(MAKE) check-fpe
	$(MAKE) valgrind
	@echo ">>> gate: статика + lean + CBMC + рантайм отработали"

# S5: межфайловая мелочь (неиспользуемые функции/поля) — целиком, на коммит.
lean:
	$(RUN) 'scripts/lean.sh'

# S6: формальный слой (CBMC) — по функциям, только на опасном.
# Обе оснастки — с точкой входа `harness` (в них нет main), параметры взяты из
# их шапок: HZ_BUFN=12 для разборщика (лестница Ш2), HZ_N=17/HZ_L2=2 для
# октодерева. ВАЖНО: у оснастки октодерева в списке файлов обязан быть
# src/octree.c — иначе CBMC не имеет тела вызываемого hz_oct_validate и его
# свойства тривиально ПРОВАЛИВАЮТСЯ («no body for callee», 8 из 35; ловля 07-10).
cbmc:
	$(RUN) 'scripts/cverify.sh tests/cbmc_sceneobj.c --function harness \
	  --unwind 14 --unwinding-assertions -DHZ_BUFN=12'
	$(RUN) 'scripts/cverify.sh tests/cbmc_octree.c --function harness \
	  --unwind 24 --unwinding-assertions -DHZ_N=17 -DHZ_L2=2 src/octree.c'

# --- S7: рантайм-слой -------------------------------------------------------
SAN = build/san
$(SAN):
	mkdir -p $(SAN)

$(SAN)/test_octree: tests/test_octree.c src/octree.c | $(SAN)
	$(RUN) 'gcc $(HZ_BASE) -O1 $(HZ_SAN) $(HZ_WARN) $(HZ_INC) -o $@ tests/test_octree.c src/octree.c -lm'
$(SAN)/test_octfmt: tests/test_octfmt.c src/octree.c | $(SAN)
	$(RUN) 'gcc $(HZ_BASE) -O1 $(HZ_SAN) $(HZ_WARN) $(HZ_INC) -o $@ tests/test_octfmt.c src/octree.c -lm'
$(SAN)/test_bounce: tests/test_bounce.c src/transport/dirs3.c src/transport/quad.c | $(SAN)
	$(RUN) 'gcc $(HZ_BASE) -O1 $(HZ_SAN) $(HZ_WARN) $(HZ_INC) -o $@ tests/test_bounce.c \
	  src/transport/dirs3.c src/transport/quad.c -lm'
$(SAN)/test_kit: tests/test_kit.c src/geom/kit.c | $(SAN)
	$(RUN) 'gcc $(HZ_BASE) -O1 $(HZ_SAN) $(HZ_WARN) $(HZ_INC) -o $@ tests/test_kit.c src/geom/kit.c -lm'

# ASan+UBSan на реальном входе (быстрые тесты). detect_leaks=1 — намеренно.
test-asan: $(SAN)/test_octree $(SAN)/test_octfmt $(SAN)/test_bounce $(SAN)/test_kit
	ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 $(SAN)/test_octree
	ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 $(SAN)/test_octfmt
	ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 $(SAN)/test_bounce
	ASAN_OPTIONS=detect_leaks=1 UBSAN_OPTIONS=print_stacktrace=1:halt_on_error=1 $(SAN)/test_kit
	@echo ">>> test-asan: ASan+UBSan чисто"

# Чтение неинициализированного: детерминированный мусор вместо случайного.
UNINIT = build/uninit
$(UNINIT):
	mkdir -p $(UNINIT)
$(UNINIT)/test_octree: tests/test_octree.c src/octree.c | $(UNINIT)
	$(RUN) 'gcc $(HZ_BASE) -O1 $(HZ_AUTOINIT) $(HZ_WARN) $(HZ_INC) -o $@ tests/test_octree.c src/octree.c -lm'
$(UNINIT)/test_octfmt: tests/test_octfmt.c src/octree.c | $(UNINIT)
	$(RUN) 'gcc $(HZ_BASE) -O1 $(HZ_AUTOINIT) $(HZ_WARN) $(HZ_INC) -o $@ tests/test_octfmt.c src/octree.c -lm'
$(UNINIT)/test_kit: tests/test_kit.c src/geom/kit.c | $(UNINIT)
	$(RUN) 'gcc $(HZ_BASE) -O1 $(HZ_AUTOINIT) $(HZ_WARN) $(HZ_INC) -o $@ tests/test_kit.c src/geom/kit.c -lm'

test-uninit: $(UNINIT)/test_octree $(UNINIT)/test_octfmt $(UNINIT)/test_kit
	$(UNINIT)/test_octree
	$(UNINIT)/test_octfmt
	$(UNINIT)/test_kit
	@echo ">>> test-uninit: чтений неинициализированного не проявилось"

# Целочисленная топология: clang -fsanitize=integer,implicit-conversion.
INT = build/int
$(INT):
	mkdir -p $(INT)
$(INT)/test_octree: tests/test_octree.c src/octree.c | $(INT)
	$(RUN) 'clang $(HZ_BASE) -O1 $(HZ_SAN_INT) $(HZ_WARN) $(HZ_INC) -o $@ tests/test_octree.c src/octree.c -lm'
$(INT)/test_octfmt: tests/test_octfmt.c src/octree.c | $(INT)
	$(RUN) 'clang $(HZ_BASE) -O1 $(HZ_SAN_INT) $(HZ_WARN) $(HZ_INC) -o $@ tests/test_octfmt.c src/octree.c -lm'
$(INT)/test_kit: tests/test_kit.c src/geom/kit.c | $(INT)
	$(RUN) 'clang $(HZ_BASE) -O1 $(HZ_SAN_INT) $(HZ_WARN) $(HZ_INC) -o $@ tests/test_kit.c src/geom/kit.c -lm'

check-int: $(INT)/test_octree $(INT)/test_octfmt $(INT)/test_kit
	UBSAN_OPTIONS=print_stacktrace=1 $(INT)/test_octree
	UBSAN_OPTIONS=print_stacktrace=1 $(INT)/test_octfmt
	UBSAN_OPTIONS=print_stacktrace=1 $(INT)/test_kit
	@echo ">>> check-int: целочисленных нарушений не проявилось"

# FP-ловушка (FE_INVALID|FE_DIVBYZERO) на быстрых тестах, без правок продукта.
build/libfpetrap.so: tests/fpetrap.c | build
	$(RUN) 'gcc -shared -fPIC -O1 -g $(HZ_WARN) -o $@ tests/fpetrap.c -lm'

check-fpe: build/libfpetrap.so $(FAST_TESTS)
	LD_PRELOAD=$(PWD)/build/libfpetrap.so ./build/test_octree
	LD_PRELOAD=$(PWD)/build/libfpetrap.so ./build/test_octfmt
	LD_PRELOAD=$(PWD)/build/libfpetrap.so ./build/test_bounce
	LD_PRELOAD=$(PWD)/build/libfpetrap.so ./build/test_kit
	@echo ">>> check-fpe: FP-исключений не было"

valgrind: $(FAST_TESTS)
	$(RUN) 'valgrind --quiet --error-exitcode=1 --leak-check=full ./build/test_octree'
	$(RUN) 'valgrind --quiet --error-exitcode=1 --leak-check=full ./build/test_octfmt'
	$(RUN) 'valgrind --quiet --error-exitcode=1 --leak-check=full ./build/test_bounce'
	$(RUN) 'valgrind --quiet --error-exitcode=1 --leak-check=full ./build/test_kit'
	@echo ">>> valgrind: утечек и ошибок памяти нет"

# --- S7-фаззинг: НЕДОВЕРЕННЫЙ ВХОД (tests/fuzz_obj.c, tests/fuzz_ppm.c) ------
# Фаззер дополняет CBMC, а не заменяет его: CBMC доказывает свойства в границах
# раскрутки, фаззер ищет то, что не придёт в голову автору теста. Находки
# однозначны (крэш/утечка под ASan+UBSan), поэтому в отчёт, а не в базу.
# Длительность — параметр: по умолчанию 60 с на цель (для гейта), длинные
# прогоны — осознанно (FUZZ_TIME=900 make fuzz).
FUZZ_TIME ?= 60
FUZZ = build/fuzz
$(FUZZ):
	mkdir -p $(FUZZ)/corpus_obj $(FUZZ)/corpus_ppm $(FUZZ)/artifacts

build/fuzz_obj: tests/fuzz_obj.c src/scene_obj.c src/scene_obj.h | $(FUZZ)
	$(RUN) 'clang $(HZ_BASE) -O1 $(HZ_SAN_FUZZ) $(HZ_WARN) $(HZ_INC) -o $@ \
	  tests/fuzz_obj.c src/scene_obj.c -lm'

build/fuzz_ppm: tests/fuzz_ppm.c src/image.c src/image.h | $(FUZZ)
	$(RUN) 'clang $(HZ_BASE) -O1 $(HZ_SAN_FUZZ) $(HZ_WARN) $(HZ_INC) -o $@ \
	  tests/fuzz_ppm.c src/image.c -lm'

# Зёрна: наши мелкие (tests/fuzz_seeds, в git) + реальные сцены, если assets/
# скачан (269 МБ, в git не лежит; см. scripts/fetch_scene.sh).
seed-corpus: | $(FUZZ)
	@cp -n tests/fuzz_seeds/obj/* $(FUZZ)/corpus_obj/ 2>/dev/null || true
	@cp -n tests/fuzz_seeds/ppm/* $(FUZZ)/corpus_ppm/ 2>/dev/null || true
	@if [ -d assets/synth ]; then \
	  find assets/synth -name '*.obj' -size -64k -exec cp -n {} $(FUZZ)/corpus_obj/ \; 2>/dev/null || true; \
	fi
	@echo "  корпус: obj=$$(ls $(FUZZ)/corpus_obj | wc -l) ppm=$$(ls $(FUZZ)/corpus_ppm | wc -l)"

fuzz: build/fuzz_obj build/fuzz_ppm seed-corpus
	@echo "=== фаззинг OBJ: $(FUZZ_TIME) с ==="
	./build/fuzz_obj -max_total_time=$(FUZZ_TIME) -print_final_stats=1 \
	  -artifact_prefix=$(FUZZ)/artifacts/ $(FUZZ)/corpus_obj
	@echo "=== фаззинг PPM: $(FUZZ_TIME) с ==="
	./build/fuzz_ppm -max_total_time=$(FUZZ_TIME) -print_final_stats=1 \
	  -artifact_prefix=$(FUZZ)/artifacts/ $(FUZZ)/corpus_ppm
	@echo ">>> fuzz: крэши — $(FUZZ)/artifacts/ (пусто = не было)"

# Хук на дешёвые ступени: «закоммичено» = «S0–S2 пройдены».
hooks:
	chmod +x .githooks/pre-commit
	git config core.hooksPath .githooks
	@echo ">>> hooks: core.hooksPath=.githooks (снять: git config --unset core.hooksPath)"

#Г31 - страж для старого ядра реза уехал вместе с ним(archive / geom);
#для живого слоя контракт FMA не был уговором — страж снят 12-09.

#СТРАЖ НАБОРА ИНСТРУКЦИЙ : ymm / xmm / fma в горячих объектниках живой линии.
check-simd: | build
	$(RUN) 'set -e; \
	  for f in src/nstruct/pyr.c src/nstruct/sweep.c src/scene_obj.c; do \
	    o=build/simd_$$(basename $$f .c).o; \
	    gcc $(CFLAGS) -c $$f -o $$o; \
	    y=$$(objdump -d $$o | grep -c "%ymm" || true); \
	    x=$$(objdump -d $$o | grep -c "%xmm" || true); \
	    m=$$(objdump -d $$o | grep -cE "vfmadd|vfmsub" || true); \
	    printf "   %-20s ymm %5s  xmm %5s  fma %3s\n" $$f $$y $$x $$m; \
	    [ "$$m" -eq 0 ] || { echo "   Г31 НАРУШЕН: fma в $$f"; exit 1; }; \
	  done; \
	  echo; echo "== потоки:"; \
	  nproc | sed "s/^/   ядер: /"; \
	  gcc $(CFLAGS) -o build/simd_omp tools/ompinfo.c && ./build/simd_omp'

.PHONY: all test check check-deep check-fast check-staged check-simd compile-db \
	gate gate-baseline gate-digest gate-versions lean cbmc \
	test-asan test-uninit check-int check-fpe valgrind hooks fuzz seed-corpus
