#!/usr/bin/env bash
# night8.sh — О23: город, срез против однородных уровней ПО ХВОСТУ метрики.
# Первый прогон (p50/p90) не различил срез и перемешанный контроль — оба ноль;
# А158: порча сидит в малой доле пикселей, и p90 её не видит по построению.
set -u
cd "$(dirname "$0")/.."
LOG=build/night/log8.txt
: > "$LOG"
echo "=== О23, город, хвост метрики : $(date +%H:%M:%S)" >> "$LOG"
OMP_NUM_THREADS=16 timeout 21600 ./build/plod city simp lev=8 viamerge uniform \
  load=build/lod/city.lod > build/night/o23_city.txt 2>&1
echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
grep -E "элем|КОНТРОЛЬ|ε ×" build/night/o23_city.txt >> "$LOG" 2>/dev/null
echo "ОЧЕРЕДЬ 8 ЗАКОНЧЕНА $(date +%H:%M:%S)" >> "$LOG"
