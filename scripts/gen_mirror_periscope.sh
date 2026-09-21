#!/usr/bin/env bash
# PERISCOPE (§873/T4-фальсификатор Ф-а): свет +x поворачивается ДВУМЯ
# 45°-зеркалами (x+z=3 -> -z, x+z=2.99 -> +x) и попадает на приёмник x=3.
# Эмиттер — плавающий квадрат x=0.5 (Ke 1, Kd 0, Ks 0); зеркала и приёмник
# без собственного излучения. Ожидания (walk, le=0, ksf=1):
#   - e>0 ТОЛЬКО на приёмнике (2 куска x=3);
#   - зеркала e=0 (всё отражают), эмиттер e=0 (тёмный вход);
#   - без ksf (ksf=0) свет умирает на первом зеркале: e≡0 всюду.
set -eu
out="${1:-assets/synth/mirror_periscope.obj}"
mtl="${out%.obj}.mtl"
cat > "$mtl" <<EOF
newmtl emit
Kd 0.000000 0.000000 0.000000
Ks 0.000000 0.000000 0.000000
Ke 1.000000 1.000000 1.000000
newmtl mir
Kd 0.000000 0.000000 0.000000
Ks 1.000000 1.000000 1.000000
newmtl recv
Kd 1.000000 1.000000 1.000000
Ks 0.000000 0.000000 0.000000
EOF
awk -v mtlname="$(basename "$mtl")" 'BEGIN {
  printf "# PERISCOPE (§873/T4-Ф-а): эмиттер x=0.5, зеркала x+z=3 и x+z=2.99, приёмник x=3.\n";
  printf "# Сгенерировано scripts/gen_mirror_periscope.sh — не править руками.\n";
  printf "mtllib %s\n", mtlname;
}
function vert(x, y, z) { printf "v %.9f %.9f %.9f\n", x, y, z; nv++; return nv; }
function quad(a, b, c, d) { printf "f %d %d %d\n", a, b, c; printf "f %d %d %d\n", a, c, d; }
BEGIN {
  nv = 0;
  printf "usemtl emit\n";
  a1=vert(0.5,0,0); a2=vert(0.5,1,0); a3=vert(0.5,1,1); a4=vert(0.5,0,1); quad(a1,a2,a3,a4);
  printf "usemtl mir\n";
  # зеркало 1: плоскость x+z=3, полоса z in [0.5,1] (трубки верхней половины);
  # малое пятно — иначе занавес x+z=3 блокирует +x-луч после второго отражения
  b1=vert(2,0,0.5); b2=vert(2.5,0,1); b3=vert(2.5,1,1); b4=vert(2,1,0.5); quad(b1,b2,b3,b4);
  # зеркало 2: плоскость x+z=2.75, полоса z in [0.25,0.5]: ловит луч -z,
  # разворачивает в +x; луч проходит ПОД зеркалом 1 (z<0.5) к приёмнику
  c1=vert(2.25,0,0.5); c2=vert(2.5,0,0.25); c3=vert(2.5,1,0.25); c4=vert(2.25,1,0.5); quad(c1,c2,c3,c4);
  printf "usemtl recv\n";
  # приёмник x=3
  d1=vert(3,0,0); d2=vert(3,1,0); d3=vert(3,1,1); d4=vert(3,0,1); quad(d1,d2,d3,d4);
}' > "$out"
echo "written: $out"
