#!/usr/bin/env bash
set -u
S=assets/San_Miguel/san-miguel.obj
R=result/O65_stages.txt
: > $R
run() { echo "### $*" | tee -a $R; ./build/pcell $S 1.0 "$@" 2>&1 | grep -E 'СТАДИИ|доли:|ИНВАРИАНТ ПОКРЫТИЯ: Σ|КАРТИНКА img|ПРОХОД ВИДИМОСТИ|СБОР НА ЭЛЕМЕНТЕ' | tee -a $R; echo | tee -a $R; }
run gather=pixel
run gather=elem
echo '##### НАКЛАДНЫЕ ПРИБОРА (А444): те же прогоны без часов' | tee -a $R
run gather=pixel notimer=all
run gather=elem notimer=all
echo '##### НЕГАТИВНЫЙ КОНТРОЛЬ: таймер солнца намеренно не копится' | tee -a $R
run gather=pixel notimer=sun
echo '##### СВЕРКА С РАЗБИВКОЙ §210 В ЛОБ' | tee -a $R
run gather=pixel sun=1 sky=1
