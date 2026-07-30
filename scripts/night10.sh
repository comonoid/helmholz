#!/usr/bin/env bash
# night10.sh — §68: лестница города с РАЗВЯЗАННЫМ радиусом поиска и по УГЛУ.
# Замер pcoarse уже показал: при том же допуске 0.4 м счёт падает 217706 -> 62798,
# то есть в 3.5 раза, если радиус кандидатов растёт до 6.4 м. Здесь то же самое
# проверяется на ЛЕСТНИЦЕ целиком и по срезу.
set -u
cd "$(dirname "$0")/.."
LOG=build/night/log10.txt
: > "$LOG"
run() {
  local name=$1; shift
  echo "=== $name : $* : $(date +%H:%M:%S)" >> "$LOG"
  ulimit -v 60000000 || true
  OMP_NUM_THREADS=16 timeout 21600 ./build/plod "$@" > "build/night/$name.txt" 2>&1
  echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
  grep -E "уровень|СРЕЗ по камере|КОНТРОЛЬ:" "build/night/$name.txt" >> "$LOG" 2>/dev/null
}
run o23_city_rad4   city simp lev=8 viamerge rad=4
run o23_city_rad16  city simp lev=8 viamerge rad=16
run o23_city_angle  city simp lev=8 viamerge byangle rad=4
echo "ОЧЕРЕДЬ 10 ЗАКОНЧЕНА $(date +%H:%M:%S)" >> "$LOG"
