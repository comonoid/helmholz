#!/usr/bin/env bash
set -u
S=assets/San_Miguel/san-miguel.obj
R=result/O64_gather_elem.txt
: > $R
run() { echo "### $*" | tee -a $R; ./build/pcell $S 1.0 "$@" 2>&1 | tee -a $R; echo | tee -a $R; }
run gather=pixel
cp -f img/pcell_lit.ppm img/O64_pixel.ppm
run gather=elem
cp -f img/pcell_lit.ppm img/O64_elem.ppm
echo '##### НЕГАТИВНЫЙ КОНТРОЛЬ: косвенное берётся у ЧУЖОГО элемента' | tee -a $R
run gather=elem shift=1
cp -f img/pcell_lit.ppm img/O64_shift.ppm
