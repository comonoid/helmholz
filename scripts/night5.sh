#!/usr/bin/env bash
# night5.sh — ГОРОД: четыре конфигурации лестницы (§60), последовательно.
# ЗАЧЕМ: таблица §60 снята на ЗАЛЕ, а вердикт «А основная, Б под бюджет» обязан
# проверяться на сцене, где уровни 0…3 бесполезны (А136) — то есть на городе.
set -u
cd "$(dirname "$0")/.."
OUT=build/night
mkdir -p "$OUT"
LOG="$OUT/log5.txt"
: > "$LOG"
run() {
  local name=$1; shift
  if [ -s "$OUT/$name.txt" ]; then echo "SKIP $name" >> "$LOG"; return; fi
  echo "=== $name : $* : $(date +%H:%M:%S)" >> "$LOG"
  OMP_NUM_THREADS=16 timeout 10800 ./build/plod "$@" > "$OUT/$name.txt" 2>&1
  echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
  grep -E "уровень|СРЕЗ|срез:" "$OUT/$name.txt" >> "$LOG" 2>/dev/null
}
run plod_city_bands      city simp lev=8 viamerge bands
run plod_city_bycount_b  city simp lev=8 viamerge bycount bands
echo "ОЧЕРЕДЬ 5 ЗАКОНЧЕНА $(date +%H:%M:%S)" >> "$LOG"
