#!/usr/bin/env bash
# night11.sh — §69: КРИВЫЕ по ε на ГОРОДЕ, четыре тактики при равной цене.
# На зале кривая уже показала: угловая сортировка даёт вдевятеро меньше испорченных
# пикселей при меньшем числе элементов. Город — та сцена, ради которой всё и делается.
set -u
cd "$(dirname "$0")/.."
LOG=build/night/log11.txt
: > "$LOG"
run() {
  local name=$1; shift
  echo "=== $name : $* : $(date +%H:%M:%S)" >> "$LOG"
  ulimit -v 60000000 || true
  OMP_NUM_THREADS=16 timeout 43200 ./build/plod "$@" > "build/night/$name.txt" 2>&1
  echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
  grep -E "КРИВАЯ" "build/night/$name.txt" >> "$LOG" 2>/dev/null
}
run o23c_delta       city simp lev=8 viamerge curve
run o23c_delta_rad4  city simp lev=8 viamerge curve rad=4
run o23c_angle       city simp lev=8 viamerge byangle curve
run o23c_angle_rad4  city simp lev=8 viamerge byangle curve rad=4
echo "ОЧЕРЕДЬ 11 ЗАКОНЧЕНА $(date +%H:%M:%S)" >> "$LOG"
