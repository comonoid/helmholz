#!/usr/bin/env bash
# night13.sh — §71: ГОРОД С ДАЛЬНЕЙ КАМЕРЫ. Там, где LOD и должен работать.
# Предсказание записано до прогона: при удалении в k раз ε·R растёт в k раз, и с
# k >= 2 (ε·R > 0.7 м = половина медианной грани) срез обязан начать падать резко,
# а однородные уровни — отставать от него всё сильнее.
set -u
cd "$(dirname "$0")/.."
LOG=build/night/log13.txt
: > "$LOG"
while [ ! -f build/night/log12.txt ] || ! grep -q "ОЧЕРЕДЬ 12 ЗАКОНЧЕНА" build/night/log12.txt; do
  echo "ждём: очередь 12 (силуэт против радиуса) ещё идёт, $(date +%H:%M:%S)" >> "$LOG"
  sleep 300
done
run() {
  local name=$1; shift
  echo "=== $name : $* : $(date +%H:%M:%S)" >> "$LOG"
  ulimit -v 60000000 || true
  OMP_NUM_THREADS=16 timeout 43200 ./build/plod "$@" > "build/night/$name.txt" 2>&1
  echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
  grep -E "ОТОДВИНУТА|СРЕЗ по камере|однородный уровень" "build/night/$name.txt" >> "$LOG" 2>/dev/null
}
run o25_city_eye2  city simp lev=8 viamerge rad=4 uniform eye=2
run o25_city_eye4  city simp lev=8 viamerge rad=4 uniform eye=4
run o25_city_eye8  city simp lev=8 viamerge rad=4 uniform eye=8
echo "ОЧЕРЕДЬ 13 ЗАКОНЧЕНА $(date +%H:%M:%S)" >> "$LOG"
