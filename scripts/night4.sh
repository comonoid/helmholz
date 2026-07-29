#!/usr/bin/env bash
# night4.sh — ВЫТЯНУТОСТЬ КАК КРИТЕРИЙ (А129), а не как распределение.
#
# ЗАЧЕМ. Предложение ограничивать вытянутость (§37) было ИЗМЕРЕНО распределением,
# но ни разу не ПРИМЕНЕНО: ни одного прогона с ограничением нет, и эффект на
# качество и на число элементов неизвестен. Это недостающий замер, а не
# «объективные обстоятельства».
#
# ЧТО МЕРИТСЯ. Порог вытянутости 2, 4, 8 (отношение полуосей ковариации группы),
# сам по себе и вместе с сохранением площади. Плюс на цилиндрических кусках
# ограничение вытянутости ДОЛЖНО быть вредно (§49: у односторонней кривизны длинная
# полоска законна и дешевле круглого куска) — на городе это увидим по числу
# элементов.
set -u
cd "$(dirname "$0")/.."
OUT=build/night4
mkdir -p "$OUT"
LOG="$OUT/log.txt"
: > "$LOG"
while pgrep -x pmetric > /dev/null; do sleep 60; done

run() {
  local name=$1
  shift
  [ -s "$OUT/$name.txt" ] && { echo "SKIP $name" >> "$LOG"; return; }
  echo "=== $name : $* : $(date +%H:%M:%S)" >> "$LOG"
  OMP_NUM_THREADS=16 timeout 7200 ./build/pmetric "$@" > "$OUT/$name.txt" 2>&1
  echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
  grep -E "участков|ПИКСЕЛЯХ|LOD:" "$OUT/$name.txt" >> "$LOG" 2>/dev/null
}

for A in 2 4 8; do
  run "hall_asp${A}" hall simp notarget d=0.36 asp=$A
  run "hall_asp${A}_loss" hall simp notarget d=0.36 asp=$A loss=1.5
done
for A in 2 4 8; do
  run "city_asp${A}" city simp notarget d=2.0 asp=$A
  run "city_asp${A}_loss" city simp notarget d=2.0 asp=$A loss=1.5
done
run "city_asp4_vfit" city simp notarget d=2.0 asp=4 vfit=0.25
run "hall_asp4_vfit" hall simp notarget d=0.36 asp=4 vfit=0.25
echo "ВЫТЯНУТОСТЬ ЗАКОНЧЕНА $(date +%H:%M:%S)" >> "$LOG"
