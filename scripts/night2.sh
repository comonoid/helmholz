#!/usr/bin/env bash
# night2.sh — ЛЕСТНИЦА УРОВНЕЙ БЕЗ ЦЕЛИ ПО ЧИСЛУ УЧАСТКОВ (§45).
#
# ЗАЧЕМ ОТДЕЛЬНО ОТ `night.sh`. Там прогоны идут с целью (город — 20 000), и на
# грубых уровнях число участков упирается в неё: при δ = 8 м выходит ровно 20 000
# независимо от того, что дал бы допуск. Для вопроса «сколько полигонов на уровне»
# это неверная величина — она про цель, а не про геометрию (А126). Здесь цель
# снята, и уровень задаётся ТОЛЬКО допуском: δ, 2δ, 4δ… — двоичная лестница §27.
set -u
cd "$(dirname "$0")/.."
OUT=build/night2
mkdir -p "$OUT"
LOG="$OUT/log.txt"
: > "$LOG"

# Ждём, пока освободится машина: параллельные прогоны портят и время, и память.
# Ожидание снято 02:08 — при трёх параллельных очередях машина не освобождается,
# и очередь простояла полчаса. Дефект оснастки, а не замысла.

run() {
  local name=$1
  shift
  [ -s "$OUT/$name.txt" ] && { echo "SKIP $name" >> "$LOG"; return; }
  echo "=== $name : $* : $(date +%H:%M:%S)" >> "$LOG"
  OMP_NUM_THREADS=16 timeout 7200 ./build/pmetric "$@" > "$OUT/$name.txt" 2>&1
  echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
  grep -E "участков|ПИКСЕЛЯХ|попали в обеих" "$OUT/$name.txt" >> "$LOG" 2>/dev/null
}

# Зал: лестница от δ сегментации вверх, цель снята.
for D in 0.09 0.18 0.36 0.72 1.44; do
  run "hall_L${D}_plain" hall simp notarget d=$D
  run "hall_L${D}_vfit" hall simp notarget d=$D vfit=0.25
  run "hall_L${D}_loss15" hall simp notarget d=$D loss=1.5
  run "hall_L${D}_vfit_loss" hall simp notarget d=$D vfit=0.25 loss=1.5
done

# Город: двоичная лестница 0.5 -> 8, цель снята. ГЛАВНЫЙ замер.
for D in 0.5 1.0 2.0 4.0 8.0; do
  run "city_L${D}_plain" city simp notarget d=$D
  run "city_L${D}_vfit" city simp notarget d=$D vfit=0.25
  run "city_L${D}_loss15" city simp notarget d=$D loss=1.5
  run "city_L${D}_vfit_loss" city simp notarget d=$D vfit=0.25 loss=1.5
done
echo "ЛЕСТНИЦА ЗАКОНЧЕНА $(date +%H:%M:%S)" >> "$LOG"
