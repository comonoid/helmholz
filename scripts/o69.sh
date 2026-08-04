#!/usr/bin/env bash
set -u
S=assets/San_Miguel/san-miguel.obj
R=result/O69_self.txt
: > $R
run() { echo "### $*" | tee -a $R; ./build/pcell $S 1.0 "$@" 2>&1 | grep -E 'СБОР НА ЭЛЕМЕНТЕ|ВИДИМОСТЬ В КОСВЕННОМ|ГДЕ УПИРАЕТСЯ|ВТОРАЯ ПРОВЕРКА|ИНВАРИАНТ МОНОТ|КАРТИНКА img' | tee -a $R; echo | tee -a $R; }
run gather=elem vis=on self=off
cp -f img/pcell_lit.ppm img/O69_off.ppm
run gather=elem vis=on self=own
cp -f img/pcell_lit.ppm img/O69_own.ppm
echo '##### НЕГАТИВНЫЙ КОНТРОЛЬ: исключается СОСЕДНИЙ элемент, а не свой' | tee -a $R
run gather=elem vis=on self=next
cp -f img/pcell_lit.ppm img/O69_next.ppm
