#!/usr/bin/env bash
# night6.sh — ГОРОД, тактика по числу (§60). Отдельной очередью, потому что первый
# прогон был потерян: `pkill -f` совпал с собственной командной строкой.
#
# ОЖИДАНИЕ ПЕЧАТАЕТ, ЧЕГО ЖДЁТ (§60, дефект оснастки), и ждёт по ФАЙЛУ-ПРИЗНАКУ, а
# не по `pgrep`: шаблон `pgrep` в командной строке сам себя и находит.
set -u
cd "$(dirname "$0")/.."
OUT=build/night
LOG="$OUT/log6.txt"
: > "$LOG"
while ! grep -q "ОЧЕРЕДЬ 5 ЗАКОНЧЕНА" "$OUT/log5.txt" 2>/dev/null; do
  echo "ждём: очередь 5 (город: полосы, затем число+полосы) ещё идёт, $(date +%H:%M:%S)" >> "$LOG"
  sleep 120
done
echo "=== город, тактика по ЧИСЛУ : $(date +%H:%M:%S)" >> "$LOG"
OMP_NUM_THREADS=16 timeout 21600 ./build/plod city simp lev=8 viamerge bycount \
  > "$OUT/plod_city_bycount.txt" 2>&1
echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
grep -E "уровень|СРЕЗ|срез:" "$OUT/plod_city_bycount.txt" >> "$LOG" 2>/dev/null
echo "ОЧЕРЕДЬ 6 ЗАКОНЧЕНА $(date +%H:%M:%S)" >> "$LOG"
