#!/usr/bin/env bash
set -u
S=assets/San_Miguel/san-miguel.obj
R=result/O67_agg_scan.txt
: > $R
run() { echo "### $*" | tee -a $R; ./build/pcell $S 1.0 "$@" 2>&1 | grep -E 'ВТОРИЧНЫХ ИСТОЧНИКОВ|СБОР НА ЭЛЕМЕНТЕ|ВИДИМОСТЬ В КОСВЕННОМ|ВТОРАЯ ПРОВЕРКА|ИНВАРИАНТ МОНОТОННОСТИ' | tee -a $R; echo | tee -a $R; }
for L in 2 3 4 5 6 7 8; do run gather=elem vis=on agg=$L; done
echo '##### НЕГАТИВНЫЙ КОНТРОЛЬ: вся сцена одним источником в корне' | tee -a $R
run gather=elem vis=on agg=0
