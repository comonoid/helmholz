/* CBMC harness for the tree-file topology check (PLAN_CUT.md, Г5 / Г29 / Р-2).
 *
 * WHY A HARNESS AND NOT `cverify.sh src/octree.c`. CBMC has no body for
 * fread/FILE*, so "CBMC on hz_oct_load" is unreachable as literally stated in
 * the work order. The load was therefore split, and this drives the pure,
 * integer-only half with a nondeterministic node array.
 *
 * WHAT IS PROVED, AND WHY IT IS THE CONSUMER AND NOT THE FORMULA. After
 * hz_oct_validate accepts, the harness runs the REAL descent hz_oct_leaf on
 * nondeterministic coordinates. Its `while (child0 >= 0)` loop has no floor on
 * `size`, so a tree with a cycle or with a chain deeper than log2size makes it
 * spin; --unwinding-assertions turns that into a failed proof. So the property
 * proved is the one the invariant exists for, not a restatement of the check.
 *
 * THE SIZE OF THE PROBLEM IS A LADDER, AND IT IS PARAMETERS, NOT A REWRITE.
 * n and log2size are COMPILE-TIME here rather than nondeterministic: a symbolic
 * n makes the malloc inside hz_oct_validate a symbolic-size dynamic object, and
 * that alone dominates the cost. Each rung is a separate, cheap run, and the
 * union of the rungs covers what one symbolic run would have.
 *   HZ_N  = 1 + 8k, the node count           (1, 9, 17)
 *   HZ_L2 = log2size                          (0, 1, 2)
 *   HZ_CHILD_LO/HI = range the nondeterministic child0 is drawn from.
 * The child0 range is a HARNESS RESTRICTION AND IS DECLARED AS ONE: every test
 * hz_oct_validate makes on c is a comparison against -1, 0, i, n-8 or 8, so
 * values beyond a small window past n exercise no path that the window does
 * not. Widening it costs SAT time and proves nothing new — but it is a
 * restriction, not a theorem, and the default window is deliberately wider than
 * n on both sides.
 *
 * Run (see PLAN_CUT.md Р-2 for the ladder actually executed):
 *   scripts/cverify.sh tests/cbmc_octree.c --function harness \
 *       --unwind 24 --unwinding-assertions -DHZ_N=17 -DHZ_L2=2 src/octree.c
 */
#include "octree.h"
#include <assert.h>
#include <complex.h>
#include <stdint.h>

/* undefined on purpose: CBMC treats a body-less function as nondeterministic */
int32_t nondet_i32(void);
int nondet_i(void);

#ifndef __CPROVER
/* so the quality gate (gcc/clang-tidy/cppcheck) can parse this file too; under
 * CBMC the name is a builtin and this declaration is not compiled */
void __CPROVER_assume(int);
#endif

#ifndef HZ_N
#define HZ_N 9
#endif
#ifndef HZ_L2
#define HZ_L2 1
#endif
#ifndef HZ_CHILD_LO
#define HZ_CHILD_LO (-3)
#endif
#ifndef HZ_CHILD_HI
#define HZ_CHILD_HI (HZ_N + 8)
#endif

void harness(void);

void harness(void) {
  /* k2 is left at its static zero: the whole check is integer, and CMPLX
   * expands to __builtin_complex, for which CBMC has no body — assigning it
   * here would fail the run on the harness rather than on the code. */
  static hz_onode nodes[HZ_N];
  for (int32_t i = 0; i < HZ_N; i++) {
    int32_t c = nondet_i32();
    __CPROVER_assume(c >= HZ_CHILD_LO && c <= HZ_CHILD_HI);
    nodes[i].child0 = c;
  }

  hz_octree t;
  t.nodes = nodes;
  t.n = HZ_N;
  t.cap = HZ_N;
  t.log2size = HZ_L2;

  if (hz_oct_validate(&t) != HZ_OCT_OK) return;

  /* the invariants, re-derived rather than re-read from the checker */
  for (int32_t i = 0; i < HZ_N; i++) {
    int32_t c = nodes[i].child0;
    assert(c == -1 || (c > i && c + 8 <= HZ_N && (c - 1) % 8 == 0));
  }

  /* THE GEOMETRIC INVARIANT, ASSERTED SEPARATELY BECAUSE NOTHING ELSE SEES IT.
   * Depth <= log2size is what stops `size` in hz_oct_leaf underflowing past 1
   * to 0, where every `x >= lo + size` turns true and the octant index sticks
   * at 7 (the mechanism in Г27). Memory safety does NOT depend on it — the
   * descent stays in bounds either way — so without this loop the depth guard
   * could be deleted and every gate would stay green. Re-derived here rather
   * than re-read from hz_oct_validate. */
  unsigned char d[HZ_N];
  for (int32_t i = 0; i < HZ_N; i++)
    d[i] = 0xFFu;
  d[0] = 0;
  for (int32_t i = 0; i < HZ_N; i++) {
    if (d[i] == 0xFFu) continue;
    int32_t c = nodes[i].child0;
    if (c < 0) continue;
    assert((int)d[i] + 1 <= HZ_L2);
    for (int j = 0; j < 8; j++)
      d[c + j] = (unsigned char)(d[i] + 1);
  }

  /* the consumer: must stay in bounds AND terminate within the unwind bound.
   * hz_oct_leaf, not hz_oct_at, because it is the integer descent that the
   * invariants protect — hz_oct_at is that descent plus one array read, and its
   * double complex return aborts CBMC's symex before anything is proved. */
  int x = nondet_i(), y = nondet_i(), z = nondet_i();
  int32_t leaf = hz_oct_leaf(&t, x, y, z);
  assert(leaf >= -1 && leaf < HZ_N);
  assert(leaf < 0 || nodes[leaf].child0 == -1); /* it really is a leaf */
}
