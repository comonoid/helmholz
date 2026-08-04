#!/usr/bin/env bash
set -u
S=assets/San_Miguel/san-miguel.obj
R=result/O68_srck_scan.txt
: > $R
run() { echo "### $*" | tee -a $R; ./build/pcell $S 1.0 "$@" 2>&1 | grep -E 'СБОР НА ЭЛЕМЕНТЕ|ВИДИМОСТЬ В КОСВЕННОМ|ВТОРАЯ ПРОВЕРКА|ИНВАРИАНТ МОНОТОННОСТИ' | tee -a $R; echo | tee -a $R; }
for K in 1 2 4 8 16; do run gather=elem vis=on srck=$K; done
echo '##### НЕГАТИВНЫЙ КОНТРОЛЬ: спуска нет, источник один — корень' | tee -a $R
run gather=elem vis=on srck=0
