#!/usr/bin/env bash
# night9.sh — О23, опыт «ворота стоят, радиус растёт» НА ГОРОДЕ (§67).
# На зале уже измерено: при том же допуске 0.18 м счёт падает 202 -> 104, то есть
# ВДВОЕ, просто от того, что кандидатами становятся далёкие компланарные куски.
# Предел памяти ставится явно: при δ = 8 м на городе через ворота проходило 57 млн
# пар, и пик списка выше итога в полтора раза (удвоение realloc).
set -u
cd "$(dirname "$0")/.."
LOG=build/night/log9.txt
: > "$LOG"
echo "=== город, ворота 0.4 м, радиус растёт : $(date +%H:%M:%S)" >> "$LOG"
ulimit -v 40000000 || true
OMP_NUM_THREADS=16 timeout 21600 ./build/pcoarse city 0.05 simp radexp=0.4 \
  > build/night/o23_radexp_city.txt 2>&1
echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
grep -E "сетка +[0-9]+ |ОПЫТ" build/night/o23_radexp_city.txt >> "$LOG" 2>/dev/null
echo "ОЧЕРЕДЬ 9 ЗАКОНЧЕНА $(date +%H:%M:%S)" >> "$LOG"
