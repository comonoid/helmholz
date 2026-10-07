#!/usr/bin/env bash
# lean.sh — whole-project "deadweight" pass: unused functions, never-used struct
# members, unused/unread variables. These need whole-project context (a symbol
# used in another file must be visible) so this is NOT a per-file ccheck check —
# run it across the active source set at once.
#
# Usage:
#   scripts/lean.sh                 # default: all of src/
#   scripts/lean.sh FILE.c ...      # explicit set (still analyzed together)
#
# Per-file unused locals + dead static functions are already caught by the gate's
# S1 (gcc -Wall -Wextra); this adds the cross-file deadweight cppcheck can only
# see with all files in scope.
#
# 07-10: список файлов берётся из git (а не рукописный) — прежний рукописный
# список в Makefile пропускал tools/pgather.c и транспортный слой. Шелл — пин
# shell.nix (версии движка влияют на находки).
set -uo pipefail
HZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$HZ"

if [ $# -gt 0 ]; then
  FILES=("$@")
else
  mapfile -t FILES < <(git ls-files 'src/*.c' 'src/*.h' 'tools/*.c' 'tests/*.c')
fi
[ ${#FILES[@]} -ge 1 ] || { echo "lean: no C sources found"; exit 0; }

echo "=== lean: deadweight over ${#FILES[@]} files (whole-project) ==="
nix-shell "$HZ/shell.nix" --run "
INC=\"-I$HZ/src\"
cppcheck --enable=unusedFunction,style --inline-suppr \
  --suppress=missingIncludeSystem --suppress=normalCheckLevelMaxBranches \
  \$INC ${FILES[*]} 2>&1 \
  | grep -E 'unusedFunction|unusedStructMember|unusedVariable|unreadVariable|redundantAssignment' \
  || echo '  (no deadweight found)'
"
