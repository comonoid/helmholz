#!/usr/bin/env bash
# cap.sh — ЗАПУСК ТЯЖЁЛОГО ПРОГОНА ПОД ПОТОЛКОМ ПАМЯТИ.
#
# ЗАЧЕМ. 08-05 в 05:50 ядро сняло `cbmc` глобальным OOM-киллером: 123.7 ГБ
# резидентно при 125 ГБ машины. Умер не тот процесс, который виноват, а тот,
# который ядро выбрало; прогон при этом не оставил ни следа, ни числа. CLAUDE.md
# требует «ulimit -v a sane cap in any long-running launcher script», и требование
# было не исполнено — `scripts/o70.sh` потолка не ставил вовсе.
#
# ПОЧЕМУ CGROUP, А НЕ `ulimit -v`. `ulimit -v` ограничивает АДРЕСНОЕ пространство:
# он режет и mmap, которым ничего не занято, поэтому для солверов его приходится
# ставить втрое выше нужного — и тогда он не защищает. Потолок cgroup v2
# (`MemoryMax`) считает РЕЗИДЕНТНОЕ и убивает ровно наш процесс, не трогая
# остальную машину. `MemorySwapMax=0` дописан затем, чтобы прогон не уползал в
# своп: 31 ГБ свопа превратили бы отказ по памяти в многочасовое перемалывание.
#
# ПОТОЛКА НА ПРОГОН НЕ ХВАТАЕТ, И ЭТО ВЫЯСНИЛОСЬ РАЗБОРОМ ТОГО ЖЕ СЛУЧАЯ:
# `cbmc` в тот раз шло ДВА, из двух сессий сразу. Два процесса по 60 ГБ каждый
# укладываются в свой потолок и вместе кладут машину. Поэтому прогон помещается
# ещё и в ОБЩУЮ долю `hzheavy.slice` (описание — `scripts/hzheavy.slice`), у
# которой свой потолок на ВСЕ прогоны всех сессий разом.
#
# ВЫХОД: код возврата команды; 137 значит «снят по потолку» — это ЗАКОНОМЕРНЫЙ
# исход, а не сбой стенда, и читать его надо как «задача в потолок не влезла».
#
# ПРИМЕР: scripts/cap.sh 48G ./build/pstow assets/synth/room.obj 1.0
set -u
if [ $# -lt 2 ]; then
  echo "usage: cap.sh <MAX, напр. 48G> <команда> [аргументы...]" >&2
  exit 2
fi
CAP=$1
shift

# Общая доля ставится из файла в репозитории, а не правкой ~/.config руками:
# иначе потолок есть на этой машине и его нет в истории.
SLICE_SRC=$(dirname "$0")/hzheavy.slice
SLICE_DST=${XDG_CONFIG_HOME:-$HOME/.config}/systemd/user/hzheavy.slice
if [ ! -f "$SLICE_DST" ] || ! cmp -s "$SLICE_SRC" "$SLICE_DST"; then
  mkdir -p "$(dirname "$SLICE_DST")"
  cp "$SLICE_SRC" "$SLICE_DST"
  systemctl --user daemon-reload
  echo "== общая доля hzheavy.slice установлена/обновлена из репозитория" >&2
fi

# Полный путь: systemd-run исполняет команду сам, PATH оболочки ему не наследуется.
CMD=$1
shift
CMD=$(command -v "$CMD") || {
  echo "cap.sh: не найдено: $CMD" >&2
  exit 127
}

UNIT=hzcap-$$
SLICE_MAX=$(systemctl --user show hzheavy.slice -p MemoryMax --value 2>/dev/null)
echo "== ПОТОЛОК ПАМЯТИ $CAP на прогон, общая доля hzheavy.slice ($SLICE_MAX Б) на все" >&2
echo "== своп запрещён; единица $UNIT" >&2
systemd-run --user --scope -q --unit="$UNIT" --slice=hzheavy.slice \
  -p MemoryMax="$CAP" -p MemorySwapMax=0 -p MemoryAccounting=yes \
  -- "$CMD" "$@"
RC=$?
if [ $RC -eq 137 ] || [ $RC -eq 9 ]; then
  echo "== СНЯТ ПО ПОТОЛКУ ПАМЯТИ ($CAP). Это исход замера, а не поломка." >&2
fi
exit $RC
