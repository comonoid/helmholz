#!/usr/bin/env bash
set -u
S=assets/San_Miguel/san-miguel.obj
R=result/O66_vis.txt
: > $R
run() { echo "### $*" | tee -a $R; ./build/pcell $S 1.0 "$@" 2>&1 | grep -E 'СБОР НА ЭЛЕМЕНТЕ|ВИДИМОСТЬ В КОСВЕННОМ|ВТОРАЯ ПРОВЕРКА|ИНВАРИАНТ МОНОТОННОСТИ|КАРТИНКА img|СТАДИИ КАДРА|ПРОХОД ВИДИМОСТИ' | tee -a $R; echo | tee -a $R; }
run gather=elem vis=off
cp -f img/pcell_lit.ppm img/O66_novis.ppm
run gather=elem vis=on
cp -f img/pcell_lit.ppm img/O66_vis.ppm
echo '##### НЕГАТИВНЫЙ КОНТРОЛЬ: инверсия теста видимости' | tee -a $R
run gather=elem vis=inverse
cp -f img/pcell_lit.ppm img/O66_inv.ppm
