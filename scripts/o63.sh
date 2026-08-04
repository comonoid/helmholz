#!/usr/bin/env bash
# Прогоны шага О63 (план §215). Порядок: сперва разброс (А428), потом обе
# конфигурации порознь, потом негативный контроль.
set -u
S=assets/San_Miguel/san-miguel.obj
R=result/O63_trace.txt
: > $R
run() { echo "### $*" | tee -a $R; ./build/pcell $S 1.0 "$@" 2>&1 | tee -a $R; echo | tee -a $R; }

echo "##### РАЗБРОС ПРОГОНОВ (А428): три раза одна конфигурация" | tee -a $R
for i in 1 2 3; do run trace=grid sun=64 sky=24; done
cp -f img/pcell_lit.ppm img/O63_grid.ppm

echo "##### ДЕРЕВО, та же конфигурация" | tee -a $R
run trace=tree sun=64 sky=24
cp -f img/pcell_lit.ppm img/O63_tree.ppm

echo "##### ДОЛИ ВЫЧИТАНИЕМ (протокол §210)" | tee -a $R
run trace=grid sun=1 sky=24
run trace=tree sun=1 sky=24
run trace=grid sun=1 sky=1
run trace=tree sun=1 sky=1

echo "##### НЕГАТИВНЫЙ КОНТРОЛЬ: дерево из ОДНОГО листа, кадр 64x64" | tee -a $R
run trace=tree img=64 sun=1 sky=1
run trace=tree img=64 sun=1 sky=1 leaf=10000000

echo "##### СЛИЧЕНИЕ КАРТИНОК" | tee -a $R
cmp img/O63_grid.ppm img/O63_tree.ppm >> $R 2>&1 && echo "картинки СОВПАЛИ побитово" | tee -a $R
