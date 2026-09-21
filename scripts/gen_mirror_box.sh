#!/usr/bin/env bash
# MIRROR-BOX (§873/T4-фальсификатор Ф-б): замкнутая коробка [0,1]^3.
# Стена x=0 — ЭМИТТЕР (Kd 0, Ks 0, Ke 1); остальные пять стен — ЗЕРКАЛА
# (Kd 0, Ks 1). Ожидания (walk, ksf=1, ke=1):
#   - депозиты ТОЛЬКО на эмиттерной стене (зеркала отражают всё: 1-Ks=0);
#   - баланс: излучено = поглощено + потеряно (хвост хопов уходит в lost
#     после HZ-глубины 4 — поток в замкнутой зеркальной коробке не гаснет).
set -eu
out="${1:-assets/synth/mirror_box.obj}"
mtl="${out%.obj}.mtl"
cat > "$mtl" <<EOF
newmtl emit
Kd 0.000000
Ke 1.000000 1.000000 1.000000
Ks 0.000000 0.000000 0.000000
newmtl mir
Kd 0.000000
Ks 1.000000 1.000000 1.000000
EOF
awk -v mtlname="$(basename "$mtl")" 'BEGIN {
  printf "# MIRROR-BOX (§873/T4-фальсификатор): эмиттер x=0, зеркала — остальные стены.\n";
  printf "# Сгенерировано scripts/gen_mirror_box.sh — не править руками.\n";
  printf "mtllib %s\n", mtlname;
}
function vert(x, y, z) { printf "v %.9f %.9f %.9f\n", x, y, z; nv++; return nv; }
function quad(a, b, c, d) { printf "f %d %d %d\n", a, b, c; printf "f %d %d %d\n", a, c, d; }
BEGIN {
  nv = 0;
  printf "usemtl emit\n";
  # эмиттер x=0 (два треугольника)
  a1=vert(0,0,0); a2=vert(0,1,0); a3=vert(0,1,1); a4=vert(0,0,1); quad(a1,a2,a3,a4);
  printf "usemtl mir\n";
  # x=1
  b1=vert(1,0,0); b2=vert(1,1,0); b3=vert(1,1,1); b4=vert(1,0,1); quad(b1,b2,b3,b4);
  # y=0
  c1=vert(0,0,0); c2=vert(1,0,0); c3=vert(1,0,1); c4=vert(0,0,1); quad(c1,c2,c3,c4);
  # y=1
  d1=vert(0,1,0); d2=vert(1,1,0); d3=vert(1,1,1); d4=vert(0,1,1); quad(d1,d2,d3,d4);
  # z=0
  e1=vert(0,0,0); e2=vert(1,0,0); e3=vert(1,1,0); e4=vert(0,1,0); quad(e1,e2,e3,e4);
  # z=1
  f1=vert(0,0,1); f2=vert(1,0,1); f3=vert(1,1,1); f4=vert(0,1,1); quad(f1,f2,f3,f4);
}' > "$out"
echo "written: $out"
