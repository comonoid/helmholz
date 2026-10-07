#!/usr/bin/env bash
# lean.sh — whole-project "deadweight" pass: unused functions, unused/unread
# variables, redundant assignments. Эти проверки требуют контекста ВСЕГО проекта
# (символ, использованный в другом файле, должен быть виден), поэтому это не
# пофайловая ступень S1, а отдельный проход по живому набору целиком.
#
# Usage:
#   scripts/lean.sh                 # default: весь живой слой из git
#   scripts/lean.sh FILE.c ...      # явный набор (анализируется вместе)
#
# 07-10: список файлов берётся из git (а не рукописный) — прежний рукописный
# список в Makefile пропускал tools/pgather.c и транспортный слой. Шелл — пин
# shell.nix (версии движка влияют на находки).
#
# ОСОЗНАННО НЕИСПОЛЬЗУЕМОЕ вычитается по scripts/lean-allow.txt: публичный API,
# который зовёт АРХИВНАЯ линия (archive/geom/ не собирается, но остаётся
# эталоном), и точки входа внешних фреймворков (CBMC/libFuzzer). Без вычета
# вывод lean превращался в список известного, и настоящее мёртвое в нём тонуло.
# Настоящий мёртвый код (ноль вызовов ГДЕ УГОДНО) — удалять, а не вписывать сюда.
#
# ЧТО ГЕЙТИТ, А ЧТО НЕТ. Гейтят unusedFunction / unusedVariable / unreadVariable /
# redundantAssignment — эти проверки в этом проекте работают. НЕ гейтит
# unusedStructMember: cppcheck разбирает заголовок В ОТРЫВЕ от пользователей и
# объявляет неиспользуемыми ВСЕ поля структур (проверено 07-10: в списке
# «never used» стоят hz_objmesh::nt, ::v, ::vn, ::mtl, ::lo, ::hi — а их читает
# и tools/scenechk.c, и весь разбор/солвер). Такие строки печатаются отдельной
# сводкой «к сведению»: слепота прибора не должна ни валить гейт, ни прятаться.
set -uo pipefail
HZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HZ"

if [ $# -gt 0 ]; then
  FILES=("$@")
else
  # tests/repro/ исключён: это фикстуры базы гейта, а не код (см. cgate.sh,
  # live_paths) — часть намеренно не форматирована и зовёт __CPROVER_assume
  mapfile -t FILES < <(git ls-files 'src/*.c' 'src/*.h' 'tools/*.c' 'tests/*.c' |
    grep -v '^tests/repro/')
fi
[ ${#FILES[@]} -ge 1 ] || { echo "lean: no C sources found"; exit 0; }

# Символы из allowlist: первое поле строки (само имя).
ALLOW=""
if [ -r "$HZ/scripts/lean-allow.txt" ]; then
  ALLOW="$(grep -vE '^[[:space:]]*(#|$)' "$HZ/scripts/lean-allow.txt" | awk '{print $1}' | paste -sd'|' -)"
fi

echo "=== lean: deadweight over ${#FILES[@]} files (whole-project) ==="
nix-shell "$HZ/shell.nix" --run "
INC=\"-I$HZ/src\"
out=\$(cppcheck --enable=unusedFunction,style --inline-suppr \
  --suppress=missingIncludeSystem --suppress=normalCheckLevelMaxBranches \
  \$INC ${FILES[*]} 2>&1 \
  | grep -E 'unusedFunction|unusedStructMember|unusedVariable|unreadVariable|redundantAssignment')
if [ -z \"\$out\" ]; then echo '  (no deadweight found)'; exit 0; fi

# 1. НЕ гейтит: слепота прибора на поля структур — отдельной сводкой.
struct=\$(printf '%s\n' \"\$out\" | grep -c 'unusedStructMember' || true)
if [ \"\$struct\" -gt 0 ]; then
  echo \"  К СВЕДЕНИЮ: \$struct строк unusedStructMember — cppcheck разбирает заголовок в отрыве\"
  echo '    от пользователей и зовёт «неиспользуемыми» все поля (в т.ч. hz_objmesh::nt/v/vn/mtl).'
  echo '    Гейт по ним НЕ валится; смотреть руками, если поле действительно нигде не читается.'
fi

# 2. Гейтит: мёртвые функции и переменные, минус осознанный API.
gate=\$(printf '%s\n' \"\$out\" | grep -v 'unusedStructMember')
if [ -n '$ALLOW' ]; then
  rest=\$(printf '%s\n' \"\$gate\" | grep -vE \"'($ALLOW)'\" || true)
  skipped=\$(printf '%s\n' \"\$gate\" | grep -cE \"'($ALLOW)'\" || true)
  [ \"\$skipped\" -gt 0 ] && echo \"  вычтено осознанного API (scripts/lean-allow.txt): \$skipped\"
else
  rest=\"\$gate\"
fi
if [ -n \"\$rest\" ]; then printf '%s\n' \"\$rest\"; exit 1; fi
echo '  мёртвого кода нет'
exit 0
"
