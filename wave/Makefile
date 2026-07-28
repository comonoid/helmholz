# wave/Makefile — ВОЛНОВАЯ ЛИНИЯ, ОТСТАВЛЕНА 07-28.
#
# ЭТОТ МАКЕФАЙЛ НЕ ВЫЗЫВАЕТСЯ ИЗ КОРНЕВОГО СОЗНАТЕЛЬНО. Пока в CLAUDE.md стоит
# блок отставки, работать с этой линией нельзя — ни собирать по своей инициативе,
# ни править, ни «прибирать». Файл существует, чтобы линия ВОССТАНАВЛИВАЛАСЬ
# одной командой, когда запрет снимут явно:
#
#     make -C wave test
#
# Общий слой (src/image.*, src/octree.*) остался в корне и подключается
# через -I ../src; дублировать его сюда нельзя — он живой и принадлежит переносу.
PKGS = gcc lapack blas pkg-config
RUN  = nix-shell -p $(PKGS) --run
WARN = -Wall -Wextra -Wshadow -Wconversion -Wsign-conversion -Wpointer-arith \
       -Wnull-dereference -Wcast-qual -Wwrite-strings -Wvla -Wformat=2 -Wundef \
       -Wstrict-prototypes -Wold-style-definition -Wmissing-prototypes \
       -Wdouble-promotion -Wfloat-equal
CFLAGS = -std=gnu11 -O2 -fopenmp -ffp-contract=off $(WARN) -I src -I ../src
CORE3D = src/solver3d.c src/assemble3d.c ../src/octree.c src/phi.c src/fft.c
LIBS   = -llapacke -llapack -lblas -lm

all: build/test_phi build/test_helm1d build/test_m2forms build/test_asm3d \
     build/test_solver3d build/render

build:
	mkdir -p build

build/test_phi: tests/test_phi.c src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_phi.c src/phi.c -lm'

build/test_helm1d: tests/test_helm1d.c src/helm1d.c src/helm1d.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tests/test_helm1d.c src/helm1d.c src/phi.c $(LIBS)'

build/test_m2forms: tests/test_m2forms.c src/helm1d.c src/helm1d.h src/phi.c src/phi.h \
                    src/fft.c src/fft.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tests/test_m2forms.c src/helm1d.c src/phi.c src/fft.c $(LIBS)'

build/test_mg3d: tests/test_mg3d.c $(CORE3D) src/solver3d.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tests/test_mg3d.c $(CORE3D) $(LIBS)'

build/test_asm3d: tests/test_asm3d.c src/assemble3d.c src/assemble3d.h ../src/octree.c src/phi.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_asm3d.c src/assemble3d.c ../src/octree.c src/phi.c -lm'

build/test_solver3d: tests/test_solver3d.c $(CORE3D) src/solver3d.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tests/test_solver3d.c $(CORE3D) $(LIBS)'

build/render: tools/render.c src/camera.c src/camera.h ../src/image.c ../src/image.h $(CORE3D) | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/render.c src/camera.c ../src/image.c $(CORE3D) $(LIBS)'

# M9: несущий (Трефтц) базис и перекрёстная проверка FDTD
build/carrier1d: tools/carrier1d.c src/helm1d.c src/helm1d.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier1d.c src/helm1d.c src/phi.c $(LIBS)'

build/fdtd: tools/fdtd.c ../src/image.c ../src/image.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/fdtd.c ../src/image.c -lm'

build/carrier_scale: tools/carrier_scale.c src/carrier.c src/carrier.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier_scale.c src/carrier.c src/phi.c $(LIBS)'

build/carrier_proj: tools/carrier_proj.c src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier_proj.c src/phi.c $(LIBS)'

build/test_carrier: tests/test_carrier.c src/carrier.c src/carrier.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tests/test_carrier.c src/carrier.c src/phi.c $(LIBS)'

build/carrier_shell: tools/carrier_shell.c src/carrier.c src/carrier.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier_shell.c src/carrier.c src/phi.c $(LIBS)'

build/carrier_angle: tools/carrier_angle.c src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier_angle.c src/phi.c $(LIBS)'

build/carrier_cascade: tools/carrier_cascade.c src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier_cascade.c src/phi.c $(LIBS)'

build/test_carrier_op: tests/test_carrier_op.c src/carrier.c src/carrier.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_carrier_op.c src/carrier.c src/phi.c -lm'

build/carrier_term: tools/carrier_term.c src/carrier.c src/carrier.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier_term.c src/carrier.c src/phi.c $(LIBS)'

build/scene2d: tools/scene2d.c src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/scene2d.c src/phi.c $(LIBS)'

build/carrier_2x2: tools/carrier_2x2.c src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier_2x2.c src/phi.c $(LIBS)'

build/carrier_solve: tools/carrier_solve.c src/carrier.c src/carrier.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier_solve.c src/carrier.c src/phi.c $(LIBS)'

build/carrier_iter: tools/carrier_iter.c src/carrier.c src/carrier.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier_iter.c src/carrier.c src/phi.c $(LIBS)'

build/carrier_incr: tools/carrier_incr.c src/carrier.c src/carrier.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier_incr.c src/carrier.c src/phi.c $(LIBS)'

build/test_carrier2d: tests/test_carrier2d.c src/carrier2d.c src/carrier2d.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_carrier2d.c src/carrier2d.c src/phi.c -lm'

build/test_cut2d: tests/test_cut2d.c src/cut2d.c src/cut2d.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_cut2d.c src/cut2d.c src/phi.c -lm'

build/test_nitsche2d: tests/test_nitsche2d.c src/nitsche2d.c src/nitsche2d.h src/cut2d.c \
                      src/cut2d.h src/carrier2d.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_nitsche2d.c src/nitsche2d.c src/cut2d.c \
	  src/phi.c -lm'

build/carrier_cut2d: tools/carrier_cut2d.c src/cut2d.c src/cut2d.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier_cut2d.c src/cut2d.c src/phi.c $(LIBS)'

build/test_bessel: tests/test_bessel.c src/bessel.c src/bessel.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_bessel.c src/bessel.c -lm'

build/test_mie2d: tests/test_mie2d.c src/mie2d.c src/mie2d.h src/bessel.c src/bessel.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_mie2d.c src/mie2d.c src/bessel.c -lm'

build/test_dtn2d: tests/test_dtn2d.c src/dtn2d.c src/dtn2d.h src/bessel.c src/carrier2d.c src/phi.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_dtn2d.c src/dtn2d.c src/bessel.c src/carrier2d.c src/phi.c -lm'

build/slice2d: tools/slice2d.c src/cut2d.c src/nitsche2d.c src/dtn2d.c src/bessel.c src/mie2d.c src/carrier2d.c src/phi.c | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/slice2d.c src/cut2d.c src/nitsche2d.c src/dtn2d.c src/bessel.c src/mie2d.c src/carrier2d.c src/phi.c $(LIBS)'

build/slab2d: tools/slab2d.c src/cut2d.c src/nitsche2d.c src/carrier2d.c src/phi.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/slab2d.c src/cut2d.c src/nitsche2d.c src/carrier2d.c src/phi.c -lm'

build/tdg2d: tools/tdg2d.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/tdg2d.c -lm'

build/dwrcomb: tools/dwrcomb.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/dwrcomb.c -lm'

test: build/test_phi build/test_carrier build/test_carrier_op build/test_carrier2d \
      build/test_cut2d build/test_nitsche2d build/test_bessel build/test_mie2d \
      build/test_dtn2d build/test_helm1d build/test_m2forms build/test_asm3d build/test_mg3d
	./build/test_phi
	./build/test_carrier
	./build/test_carrier_op
	./build/test_carrier2d
	./build/test_cut2d
	./build/test_nitsche2d
	./build/test_bessel
	./build/test_mie2d
	./build/test_dtn2d
	./build/test_helm1d
	./build/test_m2forms
	./build/test_asm3d
	./build/test_mg3d

test-slow: build/test_solver3d
	./build/test_solver3d

.PHONY: all test test-slow
