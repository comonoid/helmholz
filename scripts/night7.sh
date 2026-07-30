#!/usr/bin/env bash
# night7.sh — ГОРОД, тактика по числу ПОСЛЕ развода допуска и радиуса (А144).
# Прежний прогон падал (выход 139): снятие ворот делало кандидатом каждую пару.
set -u
cd "$(dirname "$0")/.."
LOG=build/night/log7.txt
: > "$LOG"
echo "=== город, тактика по ЧИСЛУ, после А144 : $(date +%H:%M:%S)" >> "$LOG"
OMP_NUM_THREADS=16 timeout 21600 ./build/plod city simp lev=8 viamerge bycount \
  > build/night/plod_city_bycount.txt 2>&1
echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
grep -E "уровень|СРЕЗ|срез:" build/night/plod_city_bycount.txt >> "$LOG" 2>/dev/null
echo "ОЧЕРЕДЬ 7 ЗАКОНЧЕНА $(date +%H:%M:%S)" >> "$LOG"
