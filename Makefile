# helmholz — НОВАЯ СТРУКТУРА (nstruct) + общий слой. Toolchain via nix-shell
# (see CLAUDE.md).
#
# 12-09: старый геометрический путь (src/p* conveyor, src/cut, transport-солвер,
# pfield/pwalk/render3 и их стенды) ПЕРЕНЕСЁН В archive/geom по указанию
# пользователя — исходники §914 создаются заново; старые якоря сверяются через
# git-историю. Волновая линия в wave/ — без изменений (запрет действует).
# Живая линия: src/nstruct/ (пирамида, свип), pgather/pref (сбор и эталон),
# общий слой: octree, scene_obj, image, transport/dirs3+quad (test_bounce).
PKGS = gcc lapack blas pkg-config
RUN  = nix-shell -p $(PKGS) --run
WARN = -Wall -Wextra -Wshadow -Wconversion -Wsign-conversion -Wpointer-arith \
       -Wnull-dereference -Wcast-qual -Wwrite-strings -Wvla -Wformat=2 -Wundef \
       -Wstrict-prototypes -Wold-style-definition -Wmissing-prototypes \
       -Wdouble-promotion -Wfloat-equal
# -ffp-contract=off: Г31 (история в archive/geom) — побитовая воспроизводимость
# чисел; -ffast-math ЗАПРЕЩЁН по той же причине.
ARCH   ?= -march=skylake
CFLAGS = -std=gnu11 -O2 $(ARCH) -fopenmp -ffp-contract=off $(WARN) -I src
LIBS   = -llapacke -llapack -lblas -lm

all: build/test_octree build/pgather build/pref

build:
	mkdir -p build

# ---- общий слой: октодерево ----
build/test_octree: tests/test_octree.c src/octree.c src/octree.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_octree.c src/octree.c -lm'

build/test_octfmt: tests/test_octfmt.c src/octree.c src/octree.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_octfmt.c src/octree.c -lm'

# ---- О75: нормировка первого отскока (А537) — живой гейт §875/§912 ----
build/test_bounce: tests/test_bounce.c src/transport/dirs3.c src/transport/dirs3.h \
                   src/transport/quad.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_bounce.c src/transport/dirs3.c \
	  src/transport/quad.c -lm'

# ---- НОВАЯ СТРУКТУРА ----
build/snap3: tools/snap3.c src/nstruct/pyr.c src/nstruct/pyr.h src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/snap3.c src/nstruct/pyr.c src/scene_obj.c -lm'

build/swee3: tools/swee3.c src/nstruct/pyr.c src/nstruct/pyr.h src/nstruct/sweep.c \
	src/nstruct/sweep.h src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/swee3.c src/nstruct/pyr.c src/nstruct/sweep.c src/scene_obj.c -lm'

build/pblock: tools/pblock.c src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pblock.c src/scene_obj.c -lm'

build/pgather: tools/pgather.c src/nstruct/pyr.c src/nstruct/pyr.h src/nstruct/sweep.c \
	src/nstruct/sweep.h src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pgather.c src/nstruct/pyr.c src/nstruct/sweep.c src/scene_obj.c -lm'

build/pref: tools/pref.c src/nstruct/pyr.c src/nstruct/pyr.h src/nstruct/sweep.c \
	src/nstruct/sweep.h src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pref.c src/nstruct/pyr.c src/nstruct/sweep.c src/scene_obj.c -lm'

build/ppmdiff: tools/ppmdiff.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/ppmdiff.c -lm'

build/pfmdiff: tools/pfmdiff.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pfmdiff.c -lm'

build/scenechk: tools/scenechk.c src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/scenechk.c src/scene_obj.c -lm'

# ---- ГЕОМЕТРИЯ v2 (§914, Ш1) ----
build/test_kit: tests/test_kit.c src/geom/kit.c src/geom/kit.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_kit.c src/geom/kit.c -lm'

# §914-Ш2': движение фронта по геометрической структуре (кусок с этажом P).
build/kitwalk: tools/kitwalk.c src/nstruct/pyr.c src/nstruct/pyr.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/kitwalk.c src/nstruct/pyr.c -lm'

# fast tests (seconds)
test: build/test_octree build/test_octfmt build/test_bounce build/test_kit
	./build/test_octree
	./build/test_octfmt
	./build/test_bounce
	./build/test_kit

check:
	scripts/ccheck.sh src/octree.c src/image.c src/scene_obj.c \
	  tests/test_octree.c tests/test_octfmt.c tests/cbmc_octree.c tests/cbmc_sceneobj.c \
	  src/nstruct/pyr.c src/nstruct/sweep.c \
	  src/geom/kit.c tests/test_kit.c \
	  tools/snap3.c tools/swee3.c tools/pblock.c tools/pgather.c tools/pref.c \
	  tools/scenechk.c tools/ppmdiff.c tools/pfmdiff.c tools/fcheck.c

# Г31-страж для старого ядра реза уехал вместе с ним (archive/geom);
# для живого слоя контракт FMA не был уговором — страж снят 12-09.

# СТРАЖ НАБОРА ИНСТРУКЦИЙ: ymm/xmm/fma в горячих объектниках живой линии.
check-simd: | build
	nix-shell -p gcc binutils --run 'set -e; \
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

.PHONY: all test check check-simd
