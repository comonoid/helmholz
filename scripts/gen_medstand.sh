#!/usr/bin/env bash
# КОМНАТА-СТЕНД СРЕДЫ (§932-Б-ИСП Ш3, Б0-а; v2 — сеточная). Изолированная
# среда с КОНТРОЛИРУЕМЫМ σℓ. v1 (крупные квады) вскрыла: регионы = кластеры
# УРОВНЯ кита (декимация), и при 2 чанках облако мешалось со стенами; v2
# дробит ВСЁ сеткой 4×4, облако = 64 квада × 32 три = 2048 три = РОВНО 16
# чankов kitmk подряд — кластер уровня кита = облако целиком:
#   V_r = 4·4·2 = 32 (футпринт 4×4, y∈[1,3]), ΣA·n = 64·q²·ŷ,
#   σℓ(ŷ) = ΣA/16 = 4q², q = √σℓ/2 → номинал σℓ(ŷ) = {5,1,0.1}.
# Обходы граней выведены (u×v = наружная нормаль), толщина стен 0.2 (§470).
# Лампа — сетка 4×4 Ke на потолке (lep исключает из замещения, А1732).
# Сгенерировано scripts/gen_medstand.sh — не править руками.
#
#   gen_medstand.sh [sigma_ном ∈ {5,1,0.1}] [выход.obj]
set -eu
sig="${1:-1}"
out="${2:-assets/synth/medstand_s1.obj}"
case "$sig" in
  5) ;; 1) ;; 0.1) ;;
  *) echo "gen_medstand: sigma_ном ∈ {5,1,0.1}" >&2; exit 2 ;;
esac
mtl="${out%.obj}.mtl"
awk -v SIG="$sig" -v OUT="$out" -v MTL="$mtl" 'BEGIN {
  G = 4;                              # сетка дробления (квад → G×G клеток)
  q = sqrt(SIG) / 2.0;               # сторона квада облака: σℓ(ŷ) = 4q² = SIG
  nv = 0;
  printf "mtllib %s\n", substr(MTL, match(MTL, /[^\/]*$/)) > OUT;
  printf "# КОМНАТА-СТЕНД СРЕДЫ §932-Б-ИСП Ш3 (Б0-а) v2, sigma=%.3f\n", SIG > OUT;
  printf "# Сгенерировано scripts/gen_medstand.sh — не править руками.\n" > OUT;

  # ---- ОБЛАКО: 64 квада (каждый G×G), футпринт [−2,2]², y∈[1,3] ----
  printf "usemtl cloud\n" > OUT;
  hq = q / 2.0; cell = q / G;
  for (k = 0; k < 64; k++) {
    i = k % 8; j = int(k / 8);
    x0 = -2.0 + i * (4.0 - q) / 7.0;   # bbox облака = ровно [−2,2]²
    z0 = -2.0 + j * (4.0 - q) / 7.0;
    y = 1.0 + 2.0 * k / 63.0;          # bbox y = ровно [1,3]
    gridquad(x0, y, z0, 0, 0, cell, cell, 0, 0, G, G); # u=+z, v=+x → нормаль +y
  }

  # ---- ОБОЛОЧКА: 6 плит, грани сеткой G×G, нормали наружу ----
  printf "usemtl wall\n" > OUT;
  vbox[0] = "-4.2 -0.2 -4.2  4.2 0.0  4.2";  # пол: нормаль −y → u=+x, v=+z
  vbox[1] = "-4.2  4.0 -4.2  4.2 4.2  4.2";  # потолок: +y → u=+z, v=+x
  vbox[2] = "-4.2 -0.2 -4.2 -4.0 4.2  4.2";  # x−: u=+z, v=+y
  vbox[3] = " 4.0 -0.2 -4.2  4.2 4.2  4.2";  # x+: u=+y, v=+z
  vbox[4] = "-4.2 -0.2 -4.2  4.2 4.2 -4.0";  # z−: u=+y, v=+x
  vbox[5] = "-4.2 -0.2  4.0  4.2 4.2  4.2";  # z+: u=+x, v=+y
  for (bi = 0; bi < 6; bi++) {
    split(vbox[bi], V, " ");
    gridbox(V[1], V[2], V[3], V[4], V[5], V[6]);
  }

  # ---- ЛАМПА: сетка 2×2 на y=3.9, нормаль −y (u=+x, v=+z) ----
  printf "usemtl lamp\n" > OUT;
  gridquad(-1.0, 3.9, -1.0, 1.0, 0, 0, 0, 0, 1.0, 2, 2);

  printf "newmtl cloud\nKd 0.5 0.5 0.5\n" > MTL;
  printf "newmtl wall\nKd 0.5 0.5 0.5\n" > MTL;
  printf "newmtl lamp\nKd 0.0 0.0 0.0\nKe 3.0 3.0 3.0\n" > MTL;
}
# сетка квадратов: угол (x,y,z), шаги (ux,uy,uz) по u и (vx,vy,vz) по v,
# nu×nv клеток; обход клетки CCW в (u,v) → нормаль u×v
function gridquad(x, y, z, ux, uy, uz, vx, vy, vz, nu, nv,   iu, jv, a, b2, c, d) {
  for (iu = 0; iu < nu; iu++)
    for (jv = 0; jv < nv; jv++) {
      a = nv0(x + iu * ux + jv * vx, y + iu * uy + jv * vy, z + iu * uz + jv * vz);
      b2 = nv0(x + (iu + 1) * ux + jv * vx, y + (iu + 1) * uy + jv * vy, z + (iu + 1) * uz + jv * vz);
      c = nv0(x + (iu + 1) * ux + (jv + 1) * vx, y + (iu + 1) * uy + (jv + 1) * vy,
              z + (iu + 1) * uz + (jv + 1) * vz);
      d = nv0(x + iu * ux + (jv + 1) * vx, y + iu * uy + (jv + 1) * vy, z + iu * uz + (jv + 1) * vz);
      printf "f %d %d %d\nf %d %d %d\n", a, b2, c, a, c, d > OUT;
    }
}
function nv0(x, y, z) { printf "v %.9f %.9f %.9f\n", x, y, z > OUT; return ++nv; }
# коробка гранями-сетками: u×v = наружная нормаль каждой грани
function gridbox(x0, y0, z0, x1, y1, z1) {
  gridquad(x0, y0, z0, x1 - x0, 0, 0, 0, 0, z1 - z0, G, G); # y0: −y
  gridquad(x0, y1, z0, 0, 0, z1 - z0, x1 - x0, 0, 0, G, G); # y1: +y
  gridquad(x0, y0, z0, 0, 0, z1 - z0, 0, y1 - y0, 0, G, G); # x0: −x
  gridquad(x1, y0, z0, 0, y1 - y0, 0, 0, 0, z1 - z0, G, G); # x1: +x
  gridquad(x0, y0, z0, 0, y1 - y0, 0, x1 - x0, 0, 0, G, G); # z0: −z
  gridquad(x0, y0, z1, x1 - x0, 0, 0, 0, y1 - y0, 0, G, G); # z1: +z
}'
echo "gen_medstand v2: sigma=$sig → $out (+ .mtl)"
