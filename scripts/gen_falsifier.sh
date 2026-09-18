#!/usr/bin/env bash
# FALSIFIER (А1576) — две параллельные пластины в пустоте, ρ=0 (итог точен на
# it=1): эмиттер на x=XE (Ke=1, Kd=0), приёмник на x=XR во всё сечение (Ke=0).
# YE — ширина эмиттера по y (0.47 — край ВНЕ плоскостей сетки; 1.0 — полная).
# Трубки, не пересёкшие эмиттер, несут ноль: у модели «реальные пересечения»
# E_приёмника пропорционально A_E — ЗАМКНУТЫЙ ответ E(YE)/E(1.0) = YE ±
# погрешность квадратуры. Перехватная модель добавляет фантомные f-депозиты
# из клеток, смежных с краем эмиттера, — отношение уплывает вверх от YE и
# сползает к YE с ростом lev. Прибор выбора модели А1576.
set -eu
xe="${1:-0.31}"
xr="${2:-0.83}"
ye="${3:-0.47}"
out="${4:-assets/synth/fals.obj}"
mtl="${out%.obj}.mtl"
cat > "$mtl" <<EOF
newmtl emit
Kd 0.000000
Ke 1.000000 1.000000 1.000000
newmtl recv
Kd 0.000000
EOF
awk -v xe="$xe" -v xr="$xr" -v ye="$ye" -v mtlname="$(basename "$mtl")" 'BEGIN {
  printf "# FALSIFIER (А1576): эмиттер x=%.4f y<=%.4f Ke=1; приёмник x=%.4f; ρ=0.\n", xe, ye, xr;
  printf "# Сгенерировано scripts/gen_falsifier.sh — не править руками.\n";
  printf "mtllib %s\n", mtlname;
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
  printf "usemtl emit\n";
  v1 = vert(xe, 0, 0); v2 = vert(xe, ye, 0); v3 = vert(xe, ye, 1); v4 = vert(xe, 0, 1);
  quad(v1, v2, v3, v4);
  printf "usemtl recv\n";
  v1 = vert(xr, 0, 0); v2 = vert(xr, 1, 0); v3 = vert(xr, 1, 1); v4 = vert(xr, 0, 1);
  quad(v1, v2, v3, v4);
  # ИНЕРТНЫЙ ЭКРАН за приёмником: делает домен шире, чтобы приёмник был
  # ВНУТРЕННЕЙ стенкой (пересечения строго внутри сегментов, без ulp-гонки
  # на выходной грани домена). Kd=0, Ke=0 — на метрику приёмника не влияет.
  printf "usemtl recv\n";
  v1 = vert(1.5, 0, 0); v2 = vert(1.5, 1, 0); v3 = vert(1.5, 1, 1); v4 = vert(1.5, 0, 1);
  quad(v1, v2, v3, v4);
}' > "$out"
echo "OK: $out (XE=$xe XR=$xr YE=$ye)"
