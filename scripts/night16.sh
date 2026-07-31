#!/usr/bin/env bash
# night16.sh — §76: три долга на ГОРОДЕ. Очередь ВОЗОБНОВЛЯЕМАЯ: прогон, чей файл
# уже содержит строку среза, пропускается. Заведено после того, как перезагрузка
# машины 07-31 убила очередь на первом из четырёх прогонов и заставила начинать
# всё заново.
set -u
cd "$(dirname "$0")/.."
LOG=build/night/log16.txt
: > "$LOG"
run() {
  local name=$1; shift
  if grep -q "СРЕЗ по камере" "build/night/$name.txt" 2>/dev/null; then
    echo "ПРОПУСК $name (уже досчитан)" >> "$LOG"
    return
  fi
  echo "=== $name : $* : $(date +%H:%M:%S)" >> "$LOG"
  ulimit -v 60000000 || true
  OMP_NUM_THREADS=16 timeout 43200 ./build/plodx "$@" > "build/night/$name.txt" 2>&1
  echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
  grep -E "АВТОКАЛИБ|КАЛИБРОВКА|уровень|СРЕЗ по камере|КОНТРОЛЬ:" "build/night/$name.txt" >> "$LOG" 2>/dev/null
}
run o28_city_nosin   city simp lev=8 viamerge rad=4 nosin uniform
run o28_city_auto    city simp lev=8 viamerge rad=4 nosin auto0
run o28_city_sqrt2   city simp lev=15 viamerge rad=4 nosin base=1.4142
run o28_city_nosin8  city simp lev=8 viamerge rad=4 nosin eye=8 uniform
echo "ОЧЕРЕДЬ 16 ЗАКОНЧЕНА $(date +%H:%M:%S)" >> "$LOG"
