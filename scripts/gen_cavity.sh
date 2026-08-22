#!/usr/bin/env bash
# ИНТЕГРИРУЮЩАЯ СФЕРА (§764) — внешний аналитический эталон отскока.
# Замкнутая полость: UV-сфера с нормалями ВНУТРЬ (обмотка перевёрнута против
# gen_spheres.sh), материал светится (Ke = 1) и отражает (Kd = ρ, параметр).
# Точное решение: косвенная доля psin = ρ/(1−ρ); ряд отскоков — ρⁿ.
# Полюса — одной вершиной (приём gen_spheres.sh: копии не сшивают рёбра).
set -eu
rho="${1:-0.5}"
out="${2:-assets/synth/cavity.obj}"
mtl="${out%.obj}.mtl"
awk -v rho="$rho" -v mtlname="$(basename "$mtl")" 'BEGIN {
  PI = 3.14159265358979323846;
  cx = 0; cy = 0; cz = 0; r = 1.6;
  nu = 48; nv = 24;
  printf "# ИНТЕГРИРУЮЩАЯ СФЕРА (§764): нормали внутрь, Ke=1, Kd=%.3f.\n", rho;
  printf "# Сгенерировано scripts/gen_cavity.sh — не править руками.\n";
  printf "mtllib %s\n", mtlname;
  printf "usemtl cavity\n";
  # вершины: полюса одной вершиной
  printf "v %.9f %.9f %.9f\n", cx, cy + r, cz;              # 1: северный
  for (j = 1; j < nv; j++) {
    th = PI * j / nv;
    for (i = 0; i < nu; i++) {
      ph = 2 * PI * i / nu;
      printf "v %.9f %.9f %.9f\n", cx + r*sin(th)*cos(ph), cy + r*cos(th), cz + r*sin(th)*sin(ph);
    }
  }
  printf "v %.9f %.9f %.9f\n", cx - 0 + cx, cy - r, cz;     # южный
  south = 1 + (nv - 1) * nu + 1;
  # грани: ОБМОТКА ВНУТРЬ (порядок против часовой снаружи -> по часовой)
  for (i = 0; i < nu; i++) {
    a = 2 + i; b = 2 + (i + 1) % nu;
    printf "f 1 %d %d\n", a, b;                              # шапка севера
  }
  for (j = 0; j < nv - 2; j++) {
    for (i = 0; i < nu; i++) {
      a = 2 + j*nu + i;         b = 2 + j*nu + (i + 1) % nu;
      c = 2 + (j+1)*nu + i;     d = 2 + (j+1)*nu + (i + 1) % nu;
      printf "f %d %d %d\n", a, c, d;
      printf "f %d %d %d\n", a, d, b;
    }
  }
  for (i = 0; i < nu; i++) {
    a = 2 + (nv-2)*nu + i; b = 2 + (nv-2)*nu + (i + 1) % nu;
    printf "f %d %d %d\n", south, b, a;                      # шапка юга
  }
}' > "$out"
cat > "$mtl" << EOF
# Материал полости §764: излучает и отражает; Kd — параметр эталона.
newmtl cavity
Kd $rho $rho $rho
Ke 1.0 1.0 1.0
EOF
echo "OK: $out (rho=$rho), $mtl"
