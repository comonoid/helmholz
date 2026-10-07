#!/usr/bin/env bash
# gen_compile_db.sh — build/compile_commands.json для живого слоя.
#
# Зачем база: clang-tidy (S3) переразбирает ВСЕ системные заголовки в каждом
# вызове (~4 с фиксированной цены на файл), а без базы флаги приходится
# дублировать в вызове. База строится из тех же scripts/cflags.mk, что и гейт,
# поэтому третьей копии флагов не возникает.
#
# Живой слой = git-отслеживаемые src/*.c, tools/*.c, tests/*.c. archive/ — архив,
# wave/ — вне работы (см. CLAUDE.md), ни то, ни другое в базу не попадает.
#
# Usage: scripts/gen_compile_db.sh [ВЫХОД.json]     (по умолчанию build/compile_commands.json)
set -euo pipefail

HZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
# shellcheck source=scripts/cflags.sh
. "$HZ/scripts/cflags.sh"
hz_load_cflags "$HZ/scripts/cflags.mk"
cd "$HZ"

OUT="${1:-build/compile_commands.json}"
mkdir -p "$(dirname "$OUT")"

OMP=""
[ -n "${HZ_OMP_INC:-}" ] && OMP=" -isystem $HZ_OMP_INC"

{
  printf '[\n'
  first=1
  while IFS= read -r f; do
    [ -n "$f" ] || continue
    [ "$first" = 1 ] || printf ',\n'
    first=0
    # Каталог/путь в команде абсолютные: clang-tidy запускается из любого cwd.
    printf '  {"directory": "%s", "file": "%s/%s", "command": "gcc %s -O2 %s %s %s%s -c %s -o /dev/null"}' \
      "$HZ" "$HZ" "$f" "$HZ_BASE" "$ARCH" "$HZ_WARN" "$HZ_INC" "$OMP" "$f"
  done < <(git ls-files 'src/*.c' 'tools/*.c' 'tests/*.c')
  printf '\n]\n'
} >"$OUT"

n=$(grep -c '"file"' "$OUT" || true)
echo "compile_commands.json: $n записей → $OUT"
