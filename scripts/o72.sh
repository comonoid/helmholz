#!/usr/bin/env bash
# Прогоны шага О72 (Ф2 переделанный, план §251, аудит §252). Обход фронта в
# порядке октантов, один визит на узел. Потолок памяти обязателен (А492).
set -u
R=${R:-result/O72_front.txt}
CAP=${CAP:-32G}
: > $R
SM=assets/San_Miguel/san-miguel.obj
HALL=assets/conference/conference.obj
CITY=assets/rungholt/rungholt.obj

run() { echo "### $*" | tee -a $R; scripts/cap.sh $CAP ./build/pfrontx "$@" 2>&1 | tee -a $R; echo | tee -a $R; }

echo '##### НЕГАТИВНЫЙ КОНТРОЛЬ 1: ПОЛ В ПУСТОТЕ (нарушение §241.4).' | tee -a $R
echo '##### ОБЯЗАН ВЗОРВАТЬСЯ: §241.4 считает 3.4e11 ячеек на двор 70³ м.' | tee -a $R
run $HALL 0.003 px=4 voidfloor=1

echo '##### НЕГАТИВНЫЙ КОНТРОЛЬ 2: дерево из одного узла — ячеек РОВНО 1.' | tee -a $R
run $HALL 0.003 px=4 lev=0

echo '##### ЗАМЕР: три сцены, пол 4 пикселя ОТ ИСТОЧНИКА (центр сцены).' | tee -a $R
run $HALL 0.003 px=4
run $CITY 1.0 px=4
run $SM 1.0 px=4

echo '##### СКАН ПО ПОЛУ: 1, 16, 64 пикселя — входит ли пол в сбор (урок А502).' | tee -a $R
run $CITY 1.0 px=1
run $CITY 1.0 px=16
run $CITY 1.0 px=64

echo '##### ЗАНЯТОСТЬ ГАБАРИТНАЯ (верхняя оценка, А511): обязана дать БОЛЬШЕ ячеек.' | tee -a $R
run $CITY 1.0 px=4 bboxocc=1
