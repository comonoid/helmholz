#!/usr/bin/env bash
# BLOCKED-BEAM (§852/А1566-фальсификатор) — закрытая коробка [0,2]×[0,1]×[0,1],
# перегородка на x = XW (вне плоскостей сетки) делит её на два отсека.
# Приёмник — все панели дальнего отсека (центроид x > XW). Перегородка
# герметична: при it=1 депозиты на приёмнике обязаны быть РОВНО 0 при любом
# ℓ_p (обмен светом между отсеками идёт только сквозь перегородку).
# hole=1 — вариант с прямоугольной дырой y,z ∈ [HY0,HY1]×[HZ0,HZ1] в
# перегородке: НК прибора, свет обязан пройти (E_recv(it=1) > 0).
# Панели КОМПАРТМЕНТНЫЕ (5 граней на отсек), а не общие длинные — иначе
# центроиды приёмных панелей попадают в ближний отсек.
set -eu
xw="${1:-1.003}"
hole="${2:-0}"
out="${3:-assets/synth/blocked_beam.obj}"
mtl="${out%.obj}.mtl"
cat > "$mtl" <<EOF
newmtl beam
Kd 0.500000
EOF
awk -v xw="$xw" -v hole="$hole" -v mtlname="$(basename "$mtl")" 'BEGIN {
  printf "# BLOCKED-BEAM (А1566-фальсификатор): перегородка x=%.4f, hole=%d.\n", xw, hole;
  printf "# Сгенерировано scripts/gen_blocked_beam.sh — не править руками.\n";
  printf "mtllib %s\nusemtl beam\n", mtlname;
  nv = 0;
}
function vert(x, y, z) {
  printf "v %.9f %.9f %.9f\n", x, y, z;
  nv++;
  return nv;
}
function quad(a, b, c, d) {
  printf "f %d %d %d\n", a, b, c;
  printf "f %d %d %d\n", a, c, d;
}
BEGIN {
  # ---- ближний отсек [0,xw] ----
  # грань x=0
  v1 = vert(0, 0, 0); v2 = vert(0, 1, 0); v3 = vert(0, 1, 1); v4 = vert(0, 0, 1);
  quad(v1, v2, v3, v4);
  # грань y=0
  v1 = vert(0, 0, 0); v2 = vert(xw, 0, 0); v3 = vert(xw, 0, 1); v4 = vert(0, 0, 1);
  quad(v1, v2, v3, v4);
  # грань y=1
  v1 = vert(0, 1, 0); v2 = vert(xw, 1, 0); v3 = vert(xw, 1, 1); v4 = vert(0, 1, 1);
  quad(v1, v2, v3, v4);
  # грань z=0
  v1 = vert(0, 0, 0); v2 = vert(xw, 0, 0); v3 = vert(xw, 1, 0); v4 = vert(0, 1, 0);
  quad(v1, v2, v3, v4);
  # грань z=1
  v1 = vert(0, 0, 1); v2 = vert(xw, 0, 1); v3 = vert(xw, 1, 1); v4 = vert(0, 1, 1);
  quad(v1, v2, v3, v4);
  # ---- дальний отсек [xw,2] ----
  # грань x=2
  v1 = vert(2, 0, 0); v2 = vert(2, 1, 0); v3 = vert(2, 1, 1); v4 = vert(2, 0, 1);
  quad(v1, v2, v3, v4);
  # грань y=0
  v1 = vert(xw, 0, 0); v2 = vert(2, 0, 0); v3 = vert(2, 0, 1); v4 = vert(xw, 0, 1);
  quad(v1, v2, v3, v4);
  # грань y=1
  v1 = vert(xw, 1, 0); v2 = vert(2, 1, 0); v3 = vert(2, 1, 1); v4 = vert(xw, 1, 1);
  quad(v1, v2, v3, v4);
  # грань z=0
  v1 = vert(xw, 0, 0); v2 = vert(2, 0, 0); v3 = vert(2, 1, 0); v4 = vert(xw, 1, 0);
  quad(v1, v2, v3, v4);
  # грань z=1
  v1 = vert(xw, 0, 1); v2 = vert(2, 0, 1); v3 = vert(2, 1, 1); v4 = vert(xw, 1, 1);
  quad(v1, v2, v3, v4);
  # ---- перегородка x=xw ----
  if (hole < 0.5) {
    v1 = vert(xw, 0, 0); v2 = vert(xw, 1, 0); v3 = vert(xw, 1, 1); v4 = vert(xw, 0, 1);
    quad(v1, v2, v3, v4);
  } else {
    # дыра y,z ∈ [0.4,0.6]: рама из 4 панелей
    v1 = vert(xw, 0, 0); v2 = vert(xw, 1, 0); v3 = vert(xw, 1, 0.4); v4 = vert(xw, 0, 0.4);
    quad(v1, v2, v3, v4);                       # z ∈ [0,0.4]
    v1 = vert(xw, 0, 0.6); v2 = vert(xw, 1, 0.6); v3 = vert(xw, 1, 1); v4 = vert(xw, 0, 1);
    quad(v1, v2, v3, v4);                       # z ∈ [0.6,1]
    v1 = vert(xw, 0, 0.4); v2 = vert(xw, 0.4, 0.4); v3 = vert(xw, 0.4, 0.6); v4 = vert(xw, 0, 0.6);
    quad(v1, v2, v3, v4);                       # y ∈ [0,0.4]
    v1 = vert(xw, 0.6, 0.4); v2 = vert(xw, 1, 0.4); v3 = vert(xw, 1, 0.6); v4 = vert(xw, 0.6, 0.6);
    quad(v1, v2, v3, v4);                       # y ∈ [0.6,1]
  }
}' > "$out"
echo "OK: $out (xw=$xw hole=$hole)"
