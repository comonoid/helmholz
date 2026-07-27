# helmholz — build/test/check. Toolchain via nix-shell (see CLAUDE.md).
# LAPACK = reference lapack+blas, NOT openblas (stack-smash landmine, CLAUDE.md).
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
CORE3D = src/solver3d.c src/assemble3d.c src/octree.c src/phi.c src/fft.c
LIBS   = -llapacke -llapack -lblas -lm

all: build/test_phi build/test_helm1d build/test_m2forms build/test_octree \
     build/test_asm3d build/test_solver3d build/render

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

build/test_octree: tests/test_octree.c src/octree.c src/octree.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_octree.c src/octree.c -lm'

build/test_octfmt: tests/test_octfmt.c src/octree.c src/octree.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_octfmt.c src/octree.c -lm'

build/test_mg3d: tests/test_mg3d.c $(CORE3D) src/solver3d.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tests/test_mg3d.c $(CORE3D) $(LIBS)'

build/test_asm3d: tests/test_asm3d.c src/assemble3d.c src/assemble3d.h src/octree.c src/phi.c \
                  | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_asm3d.c src/assemble3d.c src/octree.c src/phi.c -lm'

build/test_solver3d: tests/test_solver3d.c $(CORE3D) src/solver3d.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tests/test_solver3d.c $(CORE3D) $(LIBS)'

build/render: tools/render.c src/camera.c src/camera.h src/image.c src/image.h $(CORE3D) | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/render.c src/camera.c src/image.c $(CORE3D) $(LIBS)'

# M9 experiments: carrier (Trefftz) basis and the FDTD cross-check
build/carrier1d: tools/carrier1d.c src/helm1d.c src/helm1d.h src/phi.c src/phi.h | build
	$(RUN) 'gcc $(CFLAGS) $$(pkg-config --cflags lapacke) -o $@ \
	  tools/carrier1d.c src/helm1d.c src/phi.c $(LIBS)'

build/fdtd: tools/fdtd.c src/image.c src/image.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tools/fdtd.c src/image.c -lm'

# fast tests (seconds..minutes); test_solver3d is the slow validation (~10 min)
test: check-fp build/test_phi build/test_carrier build/test_carrier_op build/test_carrier2d build/test_cut2d build/test_nitsche2d build/test_bessel build/test_mie2d build/test_dtn2d build/test_helm1d build/test_m2forms build/test_octree build/test_octfmt build/test_poly3 build/test_surf build/test_facet build/test_asm3d \
      build/test_mg3d build/test_sweep build/test_rte2d
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
	./build/test_octree
	./build/test_octfmt
	./build/test_poly3
	./build/test_surf
	./build/test_facet
	./build/test_asm3d
	./build/test_mg3d
	./build/test_sweep
	./build/test_rte2d

test-slow: build/test_solver3d
	./build/test_solver3d

check:
	scripts/ccheck.sh src/phi.c src/helm1d.c src/fft.c src/octree.c src/assemble3d.c \
	  src/solver3d.c src/camera.c src/image.c \
	  tests/test_phi.c tests/test_helm1d.c tests/test_m2forms.c tests/test_octree.c \
	  tests/test_octfmt.c src/cut/poly3.c tests/test_poly3.c src/cut/surf.c tests/test_surf.c tests/test_facet.c \
	  tests/test_asm3d.c tests/test_solver3d.c tests/test_mg3d.c tools/render.c \
	  tools/carrier1d.c tools/fdtd.c src/carrier.c tools/carrier_scale.c tools/carrier_proj.c \
  tests/test_carrier.c tests/test_carrier_op.c tools/carrier_term.c tools/scene2d.c tools/carrier_shell.c tools/carrier_angle.c tools/carrier_cascade.c tools/carrier_solve.c tools/carrier_iter.c tools/carrier_incr.c src/carrier2d.c tests/test_carrier2d.c src/cut2d.c tests/test_cut2d.c tools/carrier_cut2d.c src/bessel.c tests/test_bessel.c src/mie2d.c tests/test_mie2d.c src/dtn2d.c tests/test_dtn2d.c tools/slab2d.c src/nitsche2d.c tests/test_nitsche2d.c tools/slab2d.c tools/tdg2d.c tools/slice2d.c src/transport/sweep.c tests/test_sweep.c src/transport/quad.c src/transport/rte2d.c tests/test_rte2d.c

.PHONY: all test test-slow check check-fp

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

# ---- линия ПЕРЕНОСА (PLAN_TRANSPORT.md) — отдельно от волновой ----
build/test_sweep: tests/test_sweep.c src/transport/sweep.c src/transport/sweep.h \
                  src/transport/quad.c src/transport/quad.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_sweep.c src/transport/sweep.c \
	  src/transport/quad.c -lm'

build/test_rte2d: tests/test_rte2d.c src/transport/rte2d.c src/transport/rte2d.h \
                  src/transport/quad.c src/transport/quad.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_rte2d.c src/transport/rte2d.c \
	  src/transport/quad.c -lm'

build/test_poly3: tests/test_poly3.c src/cut/poly3.c src/cut/poly3.h | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_poly3.c src/cut/poly3.c -lm'

# Г31: СТРАЖ КОНФИГУРАЦИИ СБОРКИ, А НЕ ЧИСЕЛ. Побитовое совпадение выходов ядра
# контракцию НЕ ловит — измерено: с 18 fma-инструкциями внутри test_poly3
# остаётся зелёным. Причина в устройстве теста: он сравнивает выход ОДНОГО
# скомпилированного кода с самим собой, а контракция ломает согласие лишь тогда,
# когда ОДНО скалярное произведение скомпилировано по-разному в двух местах
# (разный инлайнинг). Наблюдаемая метрика поэтому — сами инструкции.
# Базовый x86-64 без -mfma их и так не даёт, поэтому проверять надо С -mfma.
check-fp:
	nix-shell -p gcc binutils --run 'gcc $(CFLAGS) -mfma -c src/cut/poly3.c -o build/poly3_fma.o && \
	  n=$$(objdump -d build/poly3_fma.o | grep -cE "vfmadd|vfmsub" || true); \
	  echo "Г31: fma-инструкций в poly3.o при -mfma = $$n (обязано быть 0)"; \
	  [ "$$n" -eq 0 ]'

build/test_surf: tests/test_surf.c src/cut/surf.c src/cut/surf.h src/cut/poly3.c src/octree.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_surf.c src/cut/surf.c src/cut/poly3.c src/octree.c -lm'

build/test_facet: tests/test_facet.c src/cut/surf.c src/cut/surf.h src/cut/poly3.c | build
	$(RUN) 'gcc $(CFLAGS) -o $@ tests/test_facet.c src/cut/surf.c src/cut/poly3.c -lm'
