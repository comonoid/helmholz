#!/usr/bin/env bash
# lodneg.sh — НЕГАТИВНЫЙ КОНТРОЛЬ ФОРМАТА .lod (§62, приёмка; А147).
#
# ЗАЧЕМ. В приёмке стояло «порча одного байта ловится», и А147 показал, что без
# КОНТРОЛЬНОЙ СУММЫ это ложное обещание: проверка вложенности не видит обмена двух
# меток внутри одного родителя. Сумма заведена — и вот проверка, что она РАБОТАЕТ.
# Три порчи, три РАЗНЫХ кода: иначе «отказал» неотличимо от «отказал не тому».
set -u
cd "$(dirname "$0")/.."
SRC=build/lod/hall.lod
T=build/lod/neg
mkdir -p "$T"
if [ ! -s "$SRC" ]; then echo "нет $SRC — сперва save="; exit 2; fi
fail=0
chk() { # chk <имя> <ожидаемый код> <файл>
  ./build/plod hall simp lev=7 viamerge "load=$3" > /dev/null 2>"$T/$1.err"
  local got
  got=$(grep -o "код [0-9]*" "$T/$1.err" | grep -o "[0-9]*" | head -1)
  got=${got:-нет}
  if [ "$got" = "$2" ]; then
    echo "  $1: код $got — ОЖИДАЛСЯ $2, СХОДИТСЯ"
  else
    echo "  $1: код $got — ОЖИДАЛСЯ $2, НЕ СХОДИТСЯ; $(cat "$T/$1.err")"
    fail=1
  fi
}
SZ=$(stat -c%s "$SRC")
# 1. ПОРЧА НАГРУЗКИ: один бит в середине файла, то есть заведомо в узлах либо
#    метках, а не в заголовке. Ожидается код 5 (сумма).
cp "$SRC" "$T/sum.lod"
printf "\xff" | dd of="$T/sum.lod" bs=1 seek=$((SZ/2)) count=1 conv=notrunc status=none
chk sum 5 "$T/sum.lod"
# 2. ОБРЕЗАНИЕ: файл короче нагрузки. Ожидается код 1 (короткое чтение).
head -c $((SZ-64)) "$SRC" > "$T/trunc.lod"
chk trunc 1 "$T/trunc.lod"
# 3. ЧУЖОЙ ФАЙЛ: порча магии. Ожидается код 2.
cp "$SRC" "$T/magic.lod"
printf "X" | dd of="$T/magic.lod" bs=1 seek=1 count=1 conv=notrunc status=none
chk magic 2 "$T/magic.lod"
# 4. ЧУЖАЯ ПЛАТФОРМА: порча свидетеля sizeof(hz_lodnode). Ожидается код 4.
cp "$SRC" "$T/abi.lod"
printf "\x63" | dd of="$T/abi.lod" bs=1 seek=16 count=1 conv=notrunc status=none
chk abi 4 "$T/abi.lod"
if [ $fail -eq 0 ]; then echo "НЕГАТИВНЫЙ КОНТРОЛЬ ФОРМАТА: ВСЕ ЧЕТЫРЕ ПОРЧИ ПОЙМАНЫ СВОИМ КОДОМ"; else echo "НЕГАТИВНЫЙ КОНТРОЛЬ ФОРМАТА: ПРОВАЛ"; fi
exit $fail
