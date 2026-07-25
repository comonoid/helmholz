# helmholz — build/test/check. Toolchain via nix-shell (see CLAUDE.md).
# LAPACK = reference lapack+blas, NOT openblas (stack-smash landmine, CLAUDE.md).
PKGS = gcc lapack blas pkg-config
RUN  = nix-shell -p $(PKGS) --run
WARN = -Wall -Wextra -Wshadow -Wconversion -Wsign-conversion -Wpointer-arith \
       -Wnull-dereference -Wcast-qual -Wwrite-strings -Wvla -Wformat=2 -Wundef \
       -Wstrict-prototypes -Wold-style-definition -Wmissing-prototypes \
       -Wdouble-promotion -Wfloat-equal
CFLAGS = -std=gnu11 -O2 -fopenmp $(WARN) -I src
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
test: build/test_phi build/test_carrier build/test_carrier_op build/test_helm1d build/test_m2forms build/test_octree build/test_asm3d \
      build/test_mg3d
	./build/test_phi
	./build/test_carrier
	./build/test_carrier_op
	./build/test_helm1d
	./build/test_m2forms
	./build/test_octree
	./build/test_asm3d
	./build/test_mg3d

test-slow: build/test_solver3d
	./build/test_solver3d

check:
	scripts/ccheck.sh src/phi.c src/helm1d.c src/fft.c src/octree.c src/assemble3d.c \
	  src/solver3d.c src/camera.c src/image.c \
	  tests/test_phi.c tests/test_helm1d.c tests/test_m2forms.c tests/test_octree.c \
	  tests/test_asm3d.c tests/test_solver3d.c tests/test_mg3d.c tools/render.c \
	  tools/carrier1d.c tools/fdtd.c src/carrier.c tools/carrier_scale.c tools/carrier_proj.c \
  tests/test_carrier.c tests/test_carrier_op.c tools/carrier_term.c tools/scene2d.c tools/carrier_shell.c tools/carrier_angle.c tools/carrier_cascade.c

.PHONY: all test test-slow check

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
