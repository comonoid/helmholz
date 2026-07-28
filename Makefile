# helmholz — build/test/check. Toolchain via nix-shell (see CLAUDE.md).
#
# ЛИНИЯ ПЕРЕНОСА И ОБЩИЙ СЛОЙ. Волновая линия ОТСТАВЛЕНА 07-28 и вынесена
# в wave/ вместе со своим Makefile; отсюда она не собирается и не тестируется
# СОЗНАТЕЛЬНО. Условия отставки и возврата — в CLAUDE.md, блок в начале файла.
# Возврат сборки: собирать в wave/ (`make -C wave`), ничего сюда не возвращая,
# пока запрет не снят явно.
PKGS = gcc lapack blas pkg-config
RUN  = nix-shell -p $(PKGS) --run
WARN = -Wall -Wextra -Wshadow -Wconversion -Wsign-conversion -Wpointer-arith \
       -Wnull-dereference -Wcast-qual -Wwrite-strings -Wvla -Wformat=2 -Wundef \
       -Wstrict-prototypes -Wold-style-definition -Wmissing-prototypes \
       -Wdouble-promotion -Wfloat-equal
# -ffp-contract=off: PLAN_CUT.md Г31. Слияние умножения-сложения меняет последний
# бит скалярного произведения в зависимости от инлайнинга, а на побитовом
# совпадении держится водонепроницаемость разреза. #pragma STDC FP_CONTRACT
# gcc не реализует — механизм только флагом.
CFLAGS = -std=gnu11 -O2 -fopenmp -ffp-contract=off $(WARN) -I src
LIBS   = -llapacke -llapack -lblas -lm

all: build/test_octree build/test_sweep3 build/test_gather3 build/render3

build:
	mkdir -p build

# ---- общий слой: октодерево ----
build/test_octree: tests/test_octree.c src/octree.c src/octree.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_octree.c src/octree.c -lm'

build/test_octfmt: tests/test_octfmt.c src/octree.c src/octree.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_octfmt.c src/octree.c -lm'

# ---- общий слой: ГЕОМЕТРИЯ РАЗРЕЗА (PLAN_CUT.md) ----
build/test_poly3: tests/test_poly3.c src/cut/poly3.c src/cut/poly3.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_poly3.c src/cut/poly3.c -lm'

build/test_surf: tests/test_surf.c src/cut/surf.c src/cut/surf.h src/cut/poly3.c src/octree.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_surf.c src/cut/surf.c src/cut/poly3.c src/octree.c src/transport/cam3.c src/transport/ray3.c src/cut/surf.c src/image.c -lm'

build/test_facet: tests/test_facet.c src/cut/surf.c src/cut/surf.h src/cut/poly3.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_facet.c src/cut/surf.c src/cut/poly3.c -lm'

build/test_qef: tests/test_qef.c src/cut/qef.c src/cut/qef.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_qef.c src/cut/qef.c -lm'

build/test_dc: tests/test_dc.c src/cut/dc.c src/cut/dc.h src/cut/qef.c src/cut/surf.c src/cut/poly3.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_dc.c src/cut/dc.c src/cut/qef.c src/cut/surf.c src/cut/poly3.c -lm'

build/test_dcwalk: tests/test_dcwalk.c src/cut/dc.c src/cut/dc.h src/cut/qef.c src/cut/surf.c src/cut/poly3.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_dcwalk.c src/cut/dc.c src/cut/qef.c src/cut/surf.c src/cut/poly3.c -lm'

# ---- линия ПЕРЕНОСА (PLAN_TRANSPORT.md) ----
build/test_sweep: tests/test_sweep.c src/transport/sweep.c src/transport/sweep.h \
                  src/transport/quad.c src/transport/quad.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_sweep.c src/transport/sweep.c \
	  src/transport/quad.c -lm'

build/test_rte2d: tests/test_rte2d.c src/transport/rte2d.c src/transport/rte2d.h \
                  src/transport/quad.c src/transport/quad.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_rte2d.c src/transport/rte2d.c \
	  src/transport/quad.c -lm'

# ---- T5а: марш по лучу, камера, первая трёхмерная картинка ----
build/test_ray3: tests/test_ray3.c src/transport/ray3.c src/transport/ray3.h \
                 src/transport/cam3.c src/transport/cam3.h src/cut/surf.c src/cut/poly3.c \
                 src/octree.c src/image.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_ray3.c src/transport/ray3.c src/transport/cam3.c \
	  src/cut/surf.c src/cut/poly3.c src/octree.c src/image.c -lm'

build/test_scat1d: tests/test_scat1d.c src/transport/scat1d.c src/transport/scat1d.h \
                   src/transport/sweep.c src/transport/quad.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_scat1d.c src/transport/scat1d.c \
	  src/transport/sweep.c src/transport/quad.c -lm'

build/test_tet3: tests/test_tet3.c src/transport/tet3.c src/transport/tet3.h \
                 src/cut/poly3.c src/cut/surf.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_tet3.c src/transport/tet3.c src/cut/poly3.c \
	  src/cut/surf.c -lm'

build/test_sweep3: tests/test_sweep3.c src/transport/sweep3.c src/transport/mesh3.c \
                   src/transport/dirs3.c src/transport/tet3.c src/transport/cut3.c src/transport/quad.c \
                   src/cut/poly3.c src/octree.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_sweep3.c src/transport/sweep3.c \
	  src/transport/mesh3.c src/transport/dirs3.c src/transport/tet3.c src/transport/cut3.c src/transport/quad.c \
	  src/cut/poly3.c src/octree.c src/transport/cam3.c src/transport/ray3.c src/cut/surf.c src/image.c -lm'

build/test_gather3: tests/test_gather3.c src/transport/gather3.c src/transport/gather3.h \
                    src/transport/mesh3.c src/transport/cut3.c src/transport/tet3.c \
                    src/transport/ray3.c src/transport/cam3.c src/cut/poly3.c src/cut/surf.c \
                    src/octree.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_gather3.c src/transport/gather3.c \
	  src/transport/mesh3.c src/transport/cut3.c src/transport/tet3.c src/transport/ray3.c \
	  src/transport/cam3.c src/cut/poly3.c src/cut/surf.c src/octree.c -lm'

build/lod3: tools/lod3.c src/transport/ray3.c src/cut/surf.c src/cut/poly3.c src/octree.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/lod3.c src/transport/ray3.c src/cut/surf.c src/cut/poly3.c src/octree.c -lm'

# ---- ПОЛИГОНАЛЬНАЯ МОДЕЛЬ ПОЛЯ (PLAN_ELEMENTS.md) ----
# Ш1: свип по δ на скачанных сценах. Не входит в `make test` — нужны assets/,
# которых в git нет; запускается руками.
build/segstat: tools/segstat.c src/scene_obj.c src/poly_seg.c src/polygon.c \
               src/cut/surf.c src/cut/poly3.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/segstat.c src/scene_obj.c src/poly_seg.c \
	  src/polygon.c src/cut/surf.c src/cut/poly3.c -lm'

# Ш2: приёмка импорта. Сама СНИМАЕТСЯ, если сцены нет (печатает SKIP), поэтому
# в `make test` входить может.
build/test_polygon: tests/test_polygon.c src/scene_obj.c src/poly_seg.c src/polygon.c \
                    src/cut/surf.c src/cut/poly3.c src/transport/cam3.c src/transport/ray3.c \
                    src/octree.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_polygon.c src/scene_obj.c src/poly_seg.c \
	  src/polygon.c src/cut/surf.c src/cut/poly3.c src/transport/cam3.c \
	  src/transport/ray3.c src/octree.c -lm'

PELEM = src/scene_obj.c src/poly_seg.c src/polygon.c src/prast.c src/psweep.c \
        src/transport/dirs3.c src/transport/quad.c src/cut/surf.c src/cut/poly3.c

# Ш3: замер внутреннего цикла. Не в `make test` — нужны assets/.
build/pkernel: tools/pkernel.c $(PELEM) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pkernel.c $(PELEM) -lm'

# Ш4: ПЕЧЬ. Сцена СОБИРАЕТСЯ в тесте, assets не нужны ⇒ в `make test` входит.
build/test_oven: tests/test_oven.c $(PELEM) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_oven.c $(PELEM) -lm'

# СРАВНЕНИЕ ДВУХ БУФЕРОВ РАДИАНСА (PFM), оснастка замеров К65 и К68.
# Метрика берётся на РАДИАНСЕ, а не на картинке: К19 измерила, что тон-маппинг
# ошибку съедает. Кромки (К20) докладываются ОТДЕЛЬНО, а не подмешиваются в
# среднее — на силуэте расхождение равно всей яркости при любой схеме.
build/pfmdiff: tools/pfmdiff.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pfmdiff.c -lm'

# ПОЛНОЦЕННЫЙ РЕНДЕР линии переноса. Картинки пишутся в img/, а НЕ в build/:
# build — только артефакты сборки, и мусорить в нём нельзя.
build/render3: tools/render3.c src/transport/gather3.c src/transport/krylov3.c src/transport/raster3.c src/transport/sweep3.c src/transport/mesh3.c src/transport/cut3.c \
               src/transport/dirs3.c src/transport/tet3.c src/transport/quad.c \
               src/transport/ray3.c src/transport/cam3.c src/cut/poly3.c src/cut/surf.c \
               src/octree.c src/image.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/render3.c src/transport/gather3.c src/transport/krylov3.c src/transport/raster3.c src/transport/sweep3.c src/transport/mesh3.c \
	  src/transport/cut3.c src/transport/dirs3.c src/transport/tet3.c src/transport/quad.c \
	  src/transport/ray3.c src/transport/cam3.c src/cut/poly3.c src/cut/surf.c src/octree.c \
	  src/image.c -lm'

# fast tests (seconds..minutes)
test: check-fp build/test_octree build/test_octfmt build/test_poly3 build/test_surf \
      build/test_facet build/test_qef build/test_dc build/test_dcwalk \
      build/test_sweep build/test_rte2d build/test_ray3 build/test_scat1d build/test_tet3 \
      build/test_sweep3 build/test_gather3 build/test_polygon build/test_oven
	./build/test_octree
	./build/test_octfmt
	./build/test_poly3
	./build/test_surf
	./build/test_facet
	./build/test_qef
	./build/test_dc
	./build/test_dcwalk
	./build/test_sweep
	./build/test_rte2d
	./build/test_ray3
	./build/test_scat1d
	./build/test_tet3
	./build/test_sweep3
	./build/test_gather3
	./build/test_polygon
	./build/test_oven

check:
	scripts/ccheck.sh src/octree.c src/image.c \
	  tests/test_octree.c tests/test_octfmt.c \
	  src/cut/poly3.c tests/test_poly3.c src/cut/surf.c tests/test_surf.c tests/test_facet.c \
	  src/cut/qef.c tests/test_qef.c src/cut/dc.c tests/test_dc.c tests/test_dcwalk.c \
	  src/transport/sweep.c tests/test_sweep.c src/transport/quad.c \
	  src/transport/rte2d.c tests/test_rte2d.c \
	  src/transport/ray3.c src/transport/cam3.c tests/test_ray3.c \
	  src/transport/scat1d.c tests/test_scat1d.c \
	  src/transport/tet3.c tests/test_tet3.c \
	  src/transport/dirs3.c src/transport/mesh3.c src/transport/sweep3.c tests/test_sweep3.c \
	  src/transport/gather3.c tests/test_gather3.c \
	  src/transport/krylov3.c src/transport/raster3.c src/transport/cut3.c \
	  src/scene_obj.c src/poly_seg.c src/polygon.c src/prast.c src/psweep.c \
	  tests/test_polygon.c tests/test_oven.c tests/cbmc_sceneobj.c \
	  tools/lod3.c tools/pfmdiff.c tools/render3.c tools/segstat.c

# Г31: СТРАЖ КОНФИГУРАЦИИ СБОРКИ, А НЕ ЧИСЕЛ. Побитовое совпадение выходов ядра
# контракцию НЕ ловит — измерено: с 18 fma-инструкциями внутри test_poly3
# остаётся зелёным. Причина в устройстве теста: он сравнивает выход ОДНОГО
# скомпилированного кода с самим собой, а контракция ломает согласие лишь тогда,
# когда ОДНО скалярное произведение скомпилировано по-разному в двух местах
# (разный инлайнинг). Наблюдаемая метрика поэтому — сами инструкции.
# Базовый x86-64 без -mfma их и так не даёт, поэтому проверять надо С -mfma.
# Р-5б добавил ТРЕТИЙ побитовый уговор — «починка равна перестройке» (Г49), и он
# держится на тех же строгих IEEE, что и водонепроницаемость ядра. Поэтому страж
# накрывает qef.c и dc.c тоже, а не один poly3.c.
check-fp: | build
	nix-shell -p gcc binutils --run 'for f in poly3 qef dc; do \
	  gcc $(CFLAGS) -mfma -c src/cut/$$f.c -o build/$${f}_fma.o || exit 1; \
	  n=$$(objdump -d build/$${f}_fma.o | grep -cE "vfmadd|vfmsub" || true); \
	  echo "Г31: fma-инструкций в $$f.o при -mfma = $$n (обязано быть 0)"; \
	  [ "$$n" -eq 0 ] || exit 1; done'

.PHONY: all test check check-fp
