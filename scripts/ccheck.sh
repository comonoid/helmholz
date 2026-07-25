#!/usr/bin/env bash
# ccheck.sh — layered static analysis for one C file (or several).
# Runs three independent engines so every generated file is vetted for
# null-deref / use-after-free / double-free / leak / UB before it lands.
#
#   gcc -fanalyzer   path-sensitive: null deref, double-free, UAF, leak, taint
#   clang-tidy       2nd model: clang-analyzer + bugprone + cert checks
#   cppcheck         independent parser: leaks, bounds, uninit, portability
#
# Usage:
#   scripts/ccheck.sh FILE.c [FILE2.c ...]
# Env:
#   CCHECK_EXTRA="-I foo"   extra compile flags for gcc/clang-tidy
#
# Exit: non-zero if any engine reports a problem. Designed to gate generated
# code: run this on every .c written before considering it done.
set -uo pipefail
[ $# -ge 1 ] || { echo "usage: $0 FILE.c [...]"; exit 2; }

HZ=/home/n/helmholz
FILES=("$@")
EXTRA="${CCHECK_EXTRA:-}"

# All engines + project deps (openblas/lapacke) in one nix-shell so includes
# resolve and nothing is re-fetched between tools.
nix-shell -p gcc clang-tools cppcheck lapack blas pkg-config --run "
set -uo pipefail
# INC = include/define flags only (safe for cppcheck); CF adds -std for gcc/clang.
INC=\"\$(pkg-config --cflags lapacke 2>/dev/null) -I$HZ/src $EXTRA\"
CF=\"\$INC -std=gnu11 -fopenmp\"
# cppcheck accepts only -I/-D from the flag soup (pkg-config may emit -fopenmp)
CPPINC=\$(printf '%s\n' \$INC | grep -E '^-(I|D)' | tr '\n' ' ')
WARN='-Wall -Wextra -Wshadow -Wconversion -Wsign-conversion -Wpointer-arith \
-Wnull-dereference -Wcast-qual -Wwrite-strings -Wvla -Wformat=2 -Wundef \
-Wstrict-prototypes -Wold-style-definition -Wmissing-prototypes -Wstrict-overflow=2 \
-Wdouble-promotion -Wfloat-equal'
rc=0
for f in ${FILES[*]}; do
  echo '================================================================'
  echo \"=== ccheck: \$f\"
  echo '================================================================'

  echo '--- [1/3] gcc: warnings (-O2) + -fanalyzer (-O0) ---'
  # Two passes, both -Werror (gcc exits 0 on plain warnings, so -Werror is what
  # makes this engine actually gate):
  #  -O2: optimizer-dependent warnings (-Wnull-dereference, -Wstrict-overflow);
  #  -O0 -fanalyzer: at -O1+ allocation-DCE deletes dead malloc/free pairs
  #      BEFORE the analyzer runs and leaks go unseen — analyzer needs -O0.
  # NIX_HARDENING_ENABLE='': nix cc-wrapper force-defines _FORTIFY_SOURCE,
  # which #warns (=fails under -Werror) at -O0; -U can't override the wrapper.
  # Analyzer pass first (superset of plain warnings at -O0), then only the
  # O2-unique lines (optimizer-dependent) — avoids printing shared warnings twice.
  o2log=\$(mktemp); anlog=\$(mktemp)
  NIX_HARDENING_ENABLE='' gcc \$CF \$WARN -Werror -O0 -fanalyzer -c \"\$f\" -o /dev/null >\"\$anlog\" 2>&1 || rc=1
  gcc \$CF \$WARN -Werror -O2 -c \"\$f\" -o /dev/null >\"\$o2log\" 2>&1 || rc=1
  cat \"\$anlog\"
  grep -Fxvf \"\$anlog\" \"\$o2log\" || true
  rm -f \"\$o2log\" \"\$anlog\"

  echo '--- [2/3] clang-tidy (analyzer+bugprone+cert) ---'
  # advisory: clang-tidy counts system-header warnings and exits non-zero on
  # them, so its findings are displayed (and acted on) but do not gate rc.
  # The file-level findings are what matter; rc is gated by gcc -fanalyzer
  # and cppcheck, which have clean per-file exit semantics.
  clang-tidy --quiet \"\$f\" -- \$CF 2>&1 | grep -vE '[0-9]+ warnings generated|warnings? treated as errors' || true

  echo '--- [3/3] cppcheck ---'
  cppcheck --enable=warning,performance,portability --std=c11 \
    --inline-suppr --suppress=missingIncludeSystem \
    --suppress=normalCheckLevelMaxBranches \
    --error-exitcode=1 \$CPPINC \"\$f\" 2>&1 || rc=1
done
echo
if [ \$rc -eq 0 ]; then echo '>>> ccheck: CLEAN'; else echo '>>> ccheck: FINDINGS ABOVE'; fi
exit \$rc
"
