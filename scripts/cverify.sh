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

# HARD MEMORY CAP. CBMC's SAT back end can grow without bound on a harness that
# is one parameter too large, and an unbounded run does not fail — it invites
# the OOM killer, which leaves NO output at all and looks exactly like "the tool
# said nothing". A cap turns that into a fast, legible failure (CLAUDE.md: put a
# ulimit in the launcher rather than trusting the OOM killer). Override with
# CVERIFY_MEM_KB; the default is deliberately far below installed RAM, because
# a run that needs more than this is a harness to shrink, not a run to feed.
CVERIFY_MEM_KB="${CVERIFY_MEM_KB:-16000000}"

nix-shell -p cbmc --run "
ulimit -v $CVERIFY_MEM_KB || true
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
