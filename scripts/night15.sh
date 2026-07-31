#!/usr/bin/env bash
# night15.sh — §72: срез с расстоянием ДО БЛИЖАЙШЕЙ ТОЧКИ узла, дальняя камера.
# Прежний замер (очередь 13) дал при ×4 срез 75984 элемента с p90 = 49 пикселей —
# критерий обещал единицу. Причина найдена: расстояние бралось до ЦЕНТРА узла, а
# грубый узел тянется на сотни метров. Теперь у узла есть радиус объемлющей сферы.
set -u
cd "$(dirname "$0")/.."
LOG=build/night/log15.txt
: > "$LOG"
run() {
  local name=$1; shift
  echo "=== $name : $* : $(date +%H:%M:%S)" >> "$LOG"
  ulimit -v 60000000 || true
  OMP_NUM_THREADS=16 timeout 43200 ./build/plod "$@" > "build/night/$name.txt" 2>&1
  echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
  grep -E "СРЕЗ по камере|однородный уровень" "build/night/$name.txt" >> "$LOG" 2>/dev/null
}
run o27_city_eye1 city simp lev=8 viamerge rad=4 uniform
run o27_city_eye2 city simp lev=8 viamerge rad=4 uniform eye=2
run o27_city_eye4 city simp lev=8 viamerge rad=4 uniform eye=4
run o27_city_eye8 city simp lev=8 viamerge rad=4 uniform eye=8
echo "ОЧЕРЕДЬ 15 ЗАКОНЧЕНА $(date +%H:%M:%S)" >> "$LOG"
