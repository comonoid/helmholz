#!/usr/bin/env bash
# cverify.sh — heavy formal layer (CBMC bounded model checking) for ONE file.
# Exhaustively proves absence of pointer/bounds/UAF/overflow/leak defects within
# the loop-unwinding bound — no false negatives inside the bound, unlike the
# heuristic engines in ccheck.sh. Slow + bounded: use on safety-critical
# functions (octree traversal/mutation, scene-file parsers, anything on
# untrusted input), NOT on every file.
#
# Usage:
#   scripts/cverify.sh FILE.c [--function NAME] [--unwind N] [extra cbmc args]
# Examples:
#   scripts/cverify.sh src/octree.c                    # entry = main()
#   scripts/cverify.sh src/octree.c --function insert  # params = nondeterministic
#
# Exit: 0 = verified within bound; 10 = counterexample found; other = tool error.
set -uo pipefail
[ $# -ge 1 ] || { echo "usage: $0 FILE.c [--function NAME] [--unwind N] [cbmc args]"; exit 2; }
FILE="$1"; shift
HZ=/home/n/helmholz

# Default to --unwind 10 unless the caller passes their own.
case " $* " in *" --unwind "*) UNW="";; *) UNW="--unwind 10 --unwinding-assertions";; esac

nix-shell -p cbmc --run "
# NIX_HARDENING_ENABLE='': the nix cc-wrapper force-defines _FORTIFY_SOURCE at
# preprocessing, and fortified glibc wrappers (bits/stdio2.h printf) use
# __builtin_va_arg_pack, which CBMC has no body for — spurious FAILURE.
export NIX_HARDENING_ENABLE=''
INC=\"-I$HZ/src\"
echo '=== CBMC bounded verification: $FILE ==='
cbmc \"$FILE\" \$INC \
  --pointer-check --bounds-check --memory-leak-check \
  --signed-overflow-check --unsigned-overflow-check \
  --div-by-zero-check --pointer-overflow-check \
  --nan-check \
  $UNW $* 2>&1
"
