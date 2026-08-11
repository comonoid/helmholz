#!/usr/bin/env bash
# Прогоны шага О70 (план §244, «ПОПРАВЛЕННЫЙ ПЛАН»). Порядок: сперва оба
# негативных контроля, потом замер по трём сценам, потом контроль над
# счётчиком побитовых расхождений.
set -u
R=result/O70_stow.txt
: > $R
SM=assets/San_Miguel/san-miguel.obj
HALL=assets/synth/room.obj  # зал удалён 08-11, роль у комнаты-инструмента
CITY=assets/rungholt/rungholt.obj

# ПОТОЛОК ПАМЯТИ ОБЯЗАТЕЛЕН (CLAUDE.md: «ulimit -v a sane cap in any long-running
# launcher script»). Прежняя редакция его не ставила, и 08-05 глобальный OOM-киллер
# снял не тот процесс, который виноват. 32 ГБ — впятеро выше самого тяжёлого
# законного прогона (Сан-Мигель); больше означает не «мало памяти», а ошибку.
CAP=${CAP:-32G}
run() { echo "### $*" | tee -a $R; scripts/cap.sh $CAP ./build/pstow "$@" 2>&1 | tee -a $R; echo | tee -a $R; }

echo '##### НЕГАТИВНЫЙ КОНТРОЛЬ 1 (А478): СТАРАЯ УКЛАДКА «ВЛЕЗАЕТ ЦЕЛИКОМ».' | tee -a $R
echo '##### ОБЯЗАНА воспроизвести 71 / 72.5 / 79 % (§193) и десятки тысяч в узле (§217).' | tee -a $R
run $HALL 0.003 strict=1
run $SM 1.0 strict=1
run assets/San_Miguel/san-miguel-low-poly.obj 1.0 strict=1
run $CITY 1.0 strict=1

echo '##### НЕГАТИВНЫЙ КОНТРОЛЬ 2 (НК13): leafmax > nt — кусков РОВНО nt.' | tee -a $R
run $HALL 0.003 leaf=100000000 px=0
run $CITY 1.0 leaf=100000000 px=0
run $SM 1.0 leaf=100000000 px=0

echo '##### ЗАМЕР: три сцены, срез 8/4/2/1 пикселя и до листьев.' | tee -a $R
run $HALL 0.003
run $CITY 1.0
run $SM 1.0

echo '##### ГРАДУИРОВКА ВЫКЛЮЧЕНА — плата за неё отдельным числом (А476).' | tee -a $R
run $HALL 0.003 grade=0
run $SM 1.0 grade=0 px=8

echo '##### КОНТРОЛЬ НАД СЧЁТЧИКОМ ПОБИТОВЫХ РАСХОЖДЕНИЙ: наивная вершина.' | tee -a $R
echo '##### ОБЯЗАН дать НЕНУЛЕВОЕ число расхождений, в том числе на перепаде уровней.' | tee -a $R
run $HALL 0.003 naive=1 px=2
run $SM 1.0 naive=1 px=8
