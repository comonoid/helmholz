#!/usr/bin/env bash
# city_pic.sh — §78: первая ГОРОДСКАЯ картинка со среза LOD.
# Лестница строится под пиксель кадра 512² (ε сверяется при срезе, А131), затем
# рендер: прямой свет замкнутой формой от площадки-неба плюс свип косвенного.
# Цена свипа на ста тысячах элементов не мерена ни разу — это и есть главный риск.
set -u
cd "$(dirname "$0")/.."
LOG=build/night/city_pic.txt
: > "$LOG"
if [ ! -s build/lod/city_512.lod ]; then
  echo "=== лестница под 512² : $(date +%H:%M:%S)" >> "$LOG"
  OMP_NUM_THREADS=16 ./build/plodx city lev=15 viamerge rad=4 base=1.4142 nosin \
    save=build/lod/city_512.lod >> "$LOG" 2>&1
  echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
fi
echo "=== рендер : $(date +%H:%M:%S)" >> "$LOG"
OMP_NUM_THREADS=16 timeout 43200 ./build/prender 0.05 0.02 0.05 0.05 2 o city \
  lod=build/lod/city_512.lod w=512 vis=4 ss=2 >> "$LOG" 2>&1
echo "    выход $? : $(date +%H:%M:%S)" >> "$LOG"
echo "ГОТОВО $(date +%H:%M:%S)" >> "$LOG"
