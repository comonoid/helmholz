#!/usr/bin/env bash
# gen_grid.sh N — сетка N экземпляров room.obj (сдвиг X/Z, зазор 10 м)
# → assets/synth/gridN.obj. room.obj использует АБСОЛЮТНЫЕ индексы вершин,
# поэтому f-строки копии сдвигаются на вершинную базу копии (k·V), иначе все
# копии ссылаются на вершины первой — сцена схлопывается (урок §839).
set -e
N=${1:?use: gen_grid.sh N}
cols=$(awk -v n="$N" 'BEGIN{printf "%d", int(sqrt(n)+0.999)}')
vc=$(grep -c '^v ' assets/synth/room.obj)
out="assets/synth/grid$N.obj"
echo "mtllib grid.mtl" > "$out"
i=0
while [ "$i" -lt "$N" ]; do
  ix=$((i % cols)); iz=$((i / cols))
  ox=$(awk -v a="$ix" 'BEGIN{printf "%.6f", a*10}')
  oz=$(awk -v a="$iz" 'BEGIN{printf "%.6f", a*10}')
  base=$((i * vc))
  awk -v ox="$ox" -v oz="$oz" -v base="$base" '
    /^v /  { printf "v %.6f %s %.6f\n", $2+ox, $3, $4+oz; next }
    /^f /  { printf "f %d %d %d\n", $2+base, $3+base, $4+base; next }
    /^mtllib/ { next }
    { print }' assets/synth/room.obj >> "$out"
  i=$((i + 1))
done
echo "$out: $N комнат, $(( $(grep -c '^f ' assets/synth/room.obj) * N )) граней, сдвиг f на $vc/копию"
