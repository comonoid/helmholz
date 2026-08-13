# helmholz —src/transport/cam3.c src/transport/ray3.c src/transport/mesh3.c src/transport/cut3.c src/transport/quad.c src/transport/tet3.c src/octree.c | build/test/check. Toolchain via nix-shell (see CLAUDE.md).
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
#
# НАБОР ИНСТРУКЦИЙ. До 07-29 сборка шла на базовом x86-64: SSE2 (он обязателен в
# ABI для double) есть, AVX2 нет ВОВСЕ — проверено счётом регистров `ymm` в
# дизассемблере, их было НОЛЬ во всех стендах.
#
# ПОЧЕМУ `skylake`, А НЕ `native`. `-march=native` на этой машине определяет
# процессор НЕВЕРНО и разворачивается в базовый `x86-64` — проверено прямо:
# пробный цикл `a[i] = b[i]*c[i] + a[i]` даёт `ymm 0` при `native` и `ymm 4` при
# `-march=skylake`. Процессор — i9-9880H (Coffee Lake): AVX2 и FMA есть,
# AVX-512 нет, и `skylake` — ровно его набор. Ставить сюда `native`, увидев его
# в чужом Makefile, значит молча остаться без вектора: тот самый класс, где
# величина перестаёт зависеть от параметра, которым её меняют.
#
# ТРИ ВЕЩИ, КОТОРЫЕ ПРИ ЭТОМ НАДО ДЕРЖАТЬ В ГОЛОВЕ:
#   1. IEEE-СЕМАНТИКА НЕ МЕНЯЕТСЯ. Без `-ffast-math`/`-fassociative-math` gcc не
#      переставляет слагаемые в редукциях, поэтому побитовые уговоры Г31/Г49
#      остаются в силе. `-ffast-math` в этом проекте ЗАПРЕЩЁН по той же причине.
#   2. FMA. `skylake` включает и его, но `-ffp-contract=off` запрещает
#      контракцию. Замерено на том же пробном цикле: с `-ffp-contract=off` —
#      `fma 0`, без него — `fma 2`. То есть оговорка Г31 держится флагом, а не
#      надеждой, и `check-fp`/`check-simd` это проверяют.
#   3. ПЕРЕНОСИМОСТЬ БИНАРНИКА теряется. Цена нулевая: считаем здесь. На другой
#      машине — `make ARCH=-march=<её>`.
# Проверка, что набор реально используется, — `make check-simd`, а не вера.
ARCH   ?= -march=skylake
CFLAGS = -std=gnu11 -O2 $(ARCH) -fopenmp -ffp-contract=off $(WARN) -I src
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
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_surf.c src/cut/surf.c src/cut/poly3.c src/octree.c src/transport/cam3.c src/transport/ray3.c src/image.c -lm'

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

# О75: НОРМИРОВКА ПЕРВОГО ОТСКОКА (А537). Эталоны аналитические, сцена не нужна
# вовсе — стенд обязан идти в `make test` и обязан идти БЫСТРО.
build/test_bounce: tests/test_bounce.c src/transport/dirs3.c src/transport/dirs3.h \
                   src/transport/quad.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_bounce.c src/transport/dirs3.c src/transport/quad.c -lm'

# Слой ЛУЧА И КРАЯ: нужен стендам, которые спрашивают у сцены «что видно».
PGEOM = src/scene_obj.c src/poly_seg.c src/polygon.c src/pedge.c src/pray.c \
        src/cut/surf.c src/cut/poly3.c src/transport/cam3.c src/transport/ray3.c src/octree.c

# §10, пункт 1: СОГЛАСОВАННОЕ УПРОЩЕНИЕ КРАЯ против независимого. Не в
# `make test` — нужны assets/.
build/psimp: tools/psimp.c $(PGEOM) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/psimp.c $(PGEOM) -lm'

# Полный набор полигональной модели: развёртка, луч, прямой свет, разрезы,
# огрубление. Стенды Ш5…Ш8 и §10 собираются из него.
PFULL = src/scene_obj.c src/poly_seg.c src/polygon.c src/pedge.c src/prast.c src/psweep.c \
        src/pray.c src/phcube.c src/ptree.c src/pvert.c src/pdirect.c src/pcut.c src/pmerge.c src/pvfit.c src/plod.c src/lodio.c src/plink.c src/pff.c src/image.c \
        src/pclip.c \
        src/cut/surf.c src/cut/poly3.c src/cut/qef.c \
        src/transport/dirs3.c src/transport/quad.c src/transport/cam3.c src/transport/ray3.c \
        src/octree.c

build/prender: tools/prender.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/prender.c $(PFULL) -lm'

build/pconv: tools/pconv.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pconv.c $(PFULL) -lm'

build/pcuts: tools/pcuts.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pcuts.c $(PFULL) -lm'

build/ptex: tools/ptex.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/ptex.c $(PFULL) -lm'

# §10, пункты 2 и 6: ЦЕНА ОГРУБЛЕНИЯ и ГОРОД. Картинки не строит СОЗНАТЕЛЬНО —
# на городе решать поле нечем и незачем, а вопрос здесь про число полигонов,
# `dmax` и время.
build/pcoarse: tools/pcoarse.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pcoarse.c $(PFULL) -lm'

# plod — иерархическое огрубление и срез LOD (§56, О16).
# РЕНДЕРЕР БЕЗ СВИПА (§84, §85): лестница -> связи -> решение -> картинка.
build/scenechk: tools/scenechk.c src/scene_obj.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/scenechk.c src/scene_obj.c -lm'

# psil — §149, шаг О50: сколько силуэтных рёбер переживает минимальный угловой
# размер источника. Кроме сетки не нужно НИЧЕГО (ни сегментации, ни лестницы),
# поэтому и собирается из одного `scene_obj.c`: замер про рёбра ВХОДА.
build/psil: tools/psil.c src/scene_obj.c src/pgrid.c src/pgrid.h src/padj.c src/padj.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/psil.c src/scene_obj.c src/pgrid.c src/padj.c -lm'

# pcell — сегментация по ячейкам против глобальной (замер О57, §191).
build/pcell: tools/pcell.c src/scene_obj.c src/ptree.c src/ptree.h src/poly_seg.c src/poly_seg.h src/pgrid.c src/pgrid.h src/ptrace.c src/ptrace.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pcell.c src/scene_obj.c src/ptree.c src/poly_seg.c src/pgrid.c src/ptrace.c -lm'

# pstow — шаг О70 (Ф1, план §244): УКЛАДКА РЕЗКОЙ. Дерево камеронезависимо,
# срез по полу в пикселях, куски выводятся отсечением и не хранятся.
build/pstow: tools/pstow.c src/scene_obj.c src/ptree.c src/ptree.h src/pclip.c src/pclip.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pstow.c src/scene_obj.c src/ptree.c src/pclip.c -lm'

# ОТДЕЛЬНОЕ ИМЯ ДЛЯ РАБОТЫ ПОВЕРХ ИДУЩЕГО ПРОГОНА — тот же довод, что у `plodx`.
build/pstowx: tools/pstow.c src/scene_obj.c src/ptree.c src/pclip.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pstow.c src/scene_obj.c src/ptree.c src/pclip.c -lm'

# pmarchx — шаг О71 (Ф2, план §247): ОБХОД И ПУСТОТА. Марш луча по ячейкам среза;
# касания как функция расстояния и распределение длины шага.
build/pmarchx: tools/pmarchx.c src/scene_obj.c src/ptree.c src/ptree.h src/pclip.c src/pclip.h src/pmarch.c src/pmarch.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pmarchx.c src/scene_obj.c src/ptree.c src/pclip.c src/pmarch.c -lm'

# pfrontx — шаг О72 (Ф2, переделка по §251): ОБХОД ФРОНТА в порядке октантов,
# один визит на узел. Пустой узел берётся целиком (§241.4).
build/pfrontx: tools/pfrontx.c src/scene_obj.c src/ptree.c src/ptree.h src/pclip.c src/pclip.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pfrontx.c src/scene_obj.c src/ptree.c src/pclip.c -lm'

# psun — шаг О73 (Ф3, план §266): ФРОНТ ОТ СОЛНЦА с состоянием на гранях.
# Огрубление вместо отсечения; shadow=0 — негативный контроль.
# Пометки невидимости с mark=1 ставит `pcull` (§278, односторонние по
# построению); mark=2 — прежний проход фронтом, оставленный ради сверки.
build/psun: tools/psun.c src/scene_obj.c src/ptree.c src/pclip.c src/pfront.c src/pfront.h src/pocc.c src/pocc.h src/pmark.c src/pmark.h src/pcull.c src/pcull.h src/image.c src/image.h src/scene_cfg.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/psun.c src/scene_obj.c src/ptree.c src/pclip.c src/pfront.c src/pocc.c src/pmark.c src/pcull.c src/image.c -lm'

# pcam — камера сцены находится ЗАМЕРОМ, а не назначается (§187). Обе прежние
# камеры были назначены на глаз и обе оказались негодными.
build/pcam: tools/pcam.c src/scene_obj.c src/pgrid.c src/pgrid.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pcam.c src/scene_obj.c src/pgrid.c -lm'

# pfront — §172, шаг О55: насыщается ли сложность фронта при ходе. Собирается
# из того же общего слоя, что psil (сетка, смежность, форм-фактор) плюс запись
# картинки — ни сегментации, ни лестницы шагу не нужно.
build/pfront: tools/pfront.c src/scene_obj.c src/pgrid.c src/pgrid.h src/padj.c src/padj.h \
              src/pff.c src/pff.h src/image.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pfront.c src/scene_obj.c src/pgrid.c src/padj.c \
	  src/pff.c src/image.c -lm'

build/prad: tools/prad.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/prad.c $(PFULL) -lm'

# ОТДЕЛЬНОЕ ИМЯ ДЛЯ РАБОТЫ ПОВЕРХ ИДУЩЕЙ ОЧЕРЕДИ. Пересборка build/plod, пока
# очередь запускает его по одному прогону, дала бы РАЗНЫЕ версии на разных
# конфигурациях замера — этот дефект уже кусался (затёртый эталон, §36).
build/plodx: tools/plod.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/plod.c $(PFULL) -lm'

build/plod: tools/plod.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/plod.c $(PFULL) -lm'

# pmetric — метрика огрубления: две сцены через ОБЩИЙ набор лучей (§34, О20).
# Стенд ЗАМЕРА: огрубитель не трогает.
build/pmetric: tools/pmetric.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pmetric.c $(PFULL) -lm'

# pfit — запас минимаксной (чебышевской) плоскости против средневзвешенной (§23).
# Стенд ЗАМЕРА: огрубитель не трогает, критерий не меняет.
build/pfit: tools/pfit.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pfit.c $(PFULL) -lm'

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
      build/test_sweep3 build/test_gather3 build/test_polygon build/test_oven build/test_bounce
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
	./build/test_bounce

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
	  src/pgrid.c src/ptrace.c tests/cbmc_ptrace.c \
	  src/pcull.c tests/cbmc_pcull.c \
	  tests/test_polygon.c tests/test_oven.c tests/test_bounce.c tests/cbmc_sceneobj.c \
	  tools/lod3.c tools/pfmdiff.c tools/render3.c tools/segstat.c
	CCHECK_EXTRA=-I$(CURDIR)/tools scripts/ccheck.sh tools/poccref.c tests/cbmc_occmap.c

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

# СТРАЖ НАБОРА ИНСТРУКЦИЙ И ПОТОКОВ. Отвечает на вопрос «а точно ли это
# используется» ИЗМЕРЕНИЕМ, а не флагом в командной строке: флаг можно передать
# и не получить ничего, если цикл не векторизуется.
#   ymm — 256-битные регистры AVX; ноль означает, что AVX не используется вовсе;
#   xmm — 128-битные (SSE2), они есть всегда, потому что ABI x86-64 считает
#         double именно в них, и их наличие НИЧЕГО не доказывает;
#   fma — обязано быть НОЛЬ (Г31), см. check-fp.
# Отдельно печатается, что gcc реально включил (`-Q --help=target`) и сколько
# потоков видит OpenMP.
check-simd: | build
	nix-shell -p gcc binutils --run 'set -e; \
	  echo "== что включил компилятор:"; \
	  gcc $(ARCH) -Q --help=target 2>/dev/null | grep -E "^ *-m(avx2|avx512f|fma|sse2) " || true; \
	  echo; echo "== регистры в горячих объектниках:"; \
	  for f in src/psweep.c src/prast.c src/pray.c src/pmerge.c src/pdirect.c src/cut/qef.c; do \
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

.PHONY: all test check check-fp check-simd

# pvis — §110, шаг О43: наследуется ли видимость вниз по лестнице LOD. Замер
# ставится на ГОРОДЕ с камерой ВНУТРИ УЛИЦЫ; на зале он слеп по построению.
build/pvis: tools/pvis.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pvis.c $(PFULL) -lm'

# ОТДЕЛЬНОЕ ИМЯ ДЛЯ РАБОТЫ ПОВЕРХ ИДУЩЕГО ПРОГОНА — тот же довод, что у `plodx`.
build/pvisx: tools/pvis.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pvis.c $(PFULL) -lm'

# pflow — §123: перенос БЕЗ оператора, только итерации (parast + psweep).
build/pflow: tools/pflow.c $(PFULL) | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/pflow.c $(PFULL) -lm'

# pfield — источник эрмитова поля из меша (§365). Общий слой `cut/` плюс отсечение.
build/pfield: tools/pfield.c src/scene_obj.c src/pclip.c src/cut/dcslice.c src/cut/dc.c src/cut/qef.c src/cut/poly3.c src/cut/surf.c src/image.c src/transport/cam3.c src/transport/ray3.c src/transport/mesh3.c src/transport/sweep3.c src/transport/dirs3.c src/transport/cut3.c src/transport/quad.c src/transport/tet3.c src/octree.c | build
	$(RUN) 'gcc $(CFLAGS) -I tools -o $@ tools/pfield.c src/scene_obj.c src/pclip.c src/cut/dcslice.c src/cut/dc.c src/cut/qef.c src/cut/poly3.c src/cut/surf.c src/image.c src/transport/cam3.c src/transport/ray3.c src/transport/mesh3.c src/transport/sweep3.c src/transport/dirs3.c src/transport/cut3.c src/transport/quad.c src/transport/tet3.c src/octree.c -lm'

# pwalk — ТОТ ЖЕ `pfield`, но собранный с окном (SDL) и потому понимающий ключ
# `walk` (§589): цикл кадров, камера с клавиатуры и мыши.
#
# ОТДЕЛЬНАЯ ЦЕЛЬ, А НЕ ФЛАГ ПО УМОЛЧАНИЮ, И ЭТО НЕ ВКУСОВЩИНА. `build/pfield` —
# ЗАМЕРНЫЙ инструмент: на нём сняты все числа `PLAN_ELEMENTS.md`. Появись у него
# внешняя зависимость — и условия прежних прогонов изменились бы задним числом,
# а сравнивать стало бы не с чем. Исходник при этом ОДИН: растеризатор, срез и
# прямой свет не копируются, окно живёт под `#ifdef HZ_SDL`.
build/pwalk: tools/pfield.c src/scene_obj.c src/pclip.c src/cut/dcslice.c src/cut/dc.c src/cut/qef.c src/cut/poly3.c src/cut/surf.c src/image.c src/transport/cam3.c src/transport/ray3.c src/transport/mesh3.c src/transport/sweep3.c src/transport/dirs3.c src/transport/cut3.c src/transport/quad.c src/transport/tet3.c src/octree.c | build
	nix-shell -p gcc SDL2 pkg-config --run 'gcc $(CFLAGS) -DHZ_SDL -I tools -o $@ tools/pfield.c src/scene_obj.c src/pclip.c src/cut/dcslice.c src/cut/dc.c src/cut/qef.c src/cut/poly3.c src/cut/surf.c src/image.c src/transport/cam3.c src/transport/ray3.c src/transport/mesh3.c src/transport/sweep3.c src/transport/dirs3.c src/transport/cut3.c src/transport/quad.c src/transport/tet3.c src/octree.c $(shell nix-shell -p SDL2 pkg-config --run "pkg-config --cflags --libs sdl2" 2>/dev/null) -lm'

# poccref — ЭТАЛОН ЗАНЯТОСТИ (§385, шаг Ш0). Ни деревьев, ни ускорителей:
# только меш и отсечение. Он обязан пережить удаление адаптера в Ш1б, поэтому и
# стоит отдельным инструментом, а не ключом `pfield`.
build/poccref: tools/poccref.c tools/occmap.h src/scene_obj.c src/pclip.c | build
	$(RUN) 'gcc $(CFLAGS) -I tools -o $@ tools/poccref.c src/scene_obj.c src/pclip.c -lm'
