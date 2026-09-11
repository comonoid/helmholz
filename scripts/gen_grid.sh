#!/usr/bin/env bash
# gen_grid.sh N — сетка N экземпляров room.obj (сдвиг X/Z, зазор 10 м)
# → assets/synth/gridN.obj. Нормали и грани копируются как есть (нормали —
# направления, сдвигу не подлежат); сдвигаются только вершины (x += ox, z += oz).
set -e
N=${1:?use: gen_grid.sh N}
cols=$(awk -v n="$N" 'BEGIN{printf "%d", int(sqrt(n)+0.999)}')
out="assets/synth/grid$N.obj"
tmp=$(mktemp -d assets/synth/.gridXXXX)
cp assets/synth/room.mtl "$tmp/grid.mtl"
echo "mtllib grid.mtl" > "$out"
i=0
while [ "$i" -lt "$N" ]; do
  ix=$((i % cols)); iz=$((i / cols))
  ox=$(awk -v a="$ix" 'BEGIN{printf "%.6f", a*10}')
  oz=$(awk -v a="$iz" 'BEGIN{printf "%.6f", a*10}')
  awk -v ox="$ox" -v oz="$oz" '
    /^v /  { printf "v %.6f %s %.6f\n", $2+ox, $3, $4+oz; next }
    /^mtllib/ { next }
    { print }' assets/synth/room.obj >> "$out"
  i=$((i + 1))
done
rm -rf "$tmp"
echo "$out: $N комнат, $(( $(grep -c '^f ' assets/synth/room.obj) * N )) граней"
