/* Tree file format v1 (PLAN_CUT.md, Г5 / Р-2): magic, version, witnesses and
 * the topology invariants hz_oct_at depends on.
 *
 * PREDICTIONS ARE PRINTED BEFORE THE RESULTS. Every corrupted file below
 * differs from the good one in EXACTLY ONE field, and the good one is required
 * to load and to compare equal cell by cell — without that pair the whole table
 * is satisfied by `return 1;` (Г26).
 *
 * The negative control is the SIGALRM one at the end: the cyclic file must not
 * merely be rejected, the cycle must actually hang hz_oct_at when nobody
 * rejects it (Г27). If it does not hang, the invariant guards nothing. */
#include "octree.h"
#include <complex.h>
#include <setjmp.h>
#include <signal.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

static int g_fail = 0;
static int g_total = 0;

static void check(int ok, const char *what) {
  g_total++;
  if (!ok) {
    g_fail++;
    printf("FAIL: %s\n", what);
  }
}

/* --- file in a buffer ----------------------------------------------------- */

static unsigned char *save_buf(const hz_octree *t, long *len) {
  FILE *f = tmpfile();
  if (f == NULL) return NULL;
  if (hz_oct_save(t, f) != HZ_OCT_OK) {
    fclose(f);
    return NULL;
  }
  long n = ftell(f);
  if (n <= 0) {
    fclose(f);
    return NULL;
  }
  rewind(f);
  unsigned char *b = malloc((size_t)n);
  if (b == NULL) {
    fclose(f);
    return NULL;
  }
  if (fread(b, 1, (size_t)n, f) != (size_t)n) {
    free(b);
    fclose(f);
    return NULL;
  }
  fclose(f);
  *len = n;
  return b;
}

static int load_buf(const unsigned char *b, long n, hz_octree *t) {
  FILE *f = tmpfile();
  if (f == NULL) return -1;
  if (n > 0 && fwrite(b, 1, (size_t)n, f) != (size_t)n) {
    fclose(f);
    return -1;
  }
  rewind(f);
  int rc = hz_oct_load(t, f);
  fclose(f);
  return rc;
}

/* copy, patch plen bytes at off, expect exactly `want` back */
static void expect_reject(const unsigned char *base, long n, long off, const void *patch,
                          size_t plen, int want, const char *what) {
  unsigned char *b = malloc((size_t)n);
  if (b == NULL) {
    check(0, what);
    return;
  }
  memcpy(b, base, (size_t)n);
  memcpy(b + off, patch, plen);
  hz_octree t;
  int rc = load_buf(b, n, &t);
  if (rc == HZ_OCT_OK) hz_oct_free(&t);
  if (rc != want) printf("  (%s: got %d, wanted %d)\n", what, rc, want);
  check(rc == want, what);
  free(b);
}

static void patch_i32(const unsigned char *base, long n, long off, int32_t v, int want,
                      const char *what) {
  expect_reject(base, n, off, &v, sizeof v, want, what);
}

/* --- the tree under test -------------------------------------------------- */

enum { L2 = 2, S = 1 << L2 }; /* 4^3 cube, so depth 2 is the legal maximum */

static int build(hz_octree *t) {
  if (hz_oct_init(t, L2, CMPLX(1.0, 0.02)) != 0) return 1;
  int lo[3] = {0, 0, 0}, hi[3] = {1, 1, 1}; /* one unit cell: forces depth 2 */
  if (hz_oct_set_box(t, lo, hi, CMPLX(2.5, -0.1)) != 0) return 1;
  int lo2[3] = {2, 0, 0}, hi2[3] = {4, 2, 2}; /* one octant: depth 1 */
  return hz_oct_set_box(t, lo2, hi2, CMPLX(-3.0, 0.5));
}

static int same_cells(const hz_octree *a, const hz_octree *b) {
  for (int x = 0; x < S; x++)
    for (int y = 0; y < S; y++)
      for (int z = 0; z < S; z++) {
        /* bitwise, via memcmp: the file is required to be exact, and == on
         * doubles would both warn and quietly accept a changed representation */
        double complex va = hz_oct_at(a, x, y, z), vb = hz_oct_at(b, x, y, z);
        if (memcmp(&va, &vb, sizeof va) != 0) return 0;
      }
  return 1;
}

/* --- negative control: does the cycle really hang? ------------------------ */

static sigjmp_buf g_jb;
static volatile sig_atomic_t g_hit;

static void on_alarm(int s) {
  (void)s;
  g_hit = 1;
  siglongjmp(g_jb, 1);
}

/* Build IN MEMORY a cyclic tree that the OLD check (child0 + 8 <= n, and
 * nothing else) accepts, and call the consumer. Returns 1 if it did not come
 * back within the alarm.
 *
 * THE FIRST VERSION OF THIS CONTROL WAS WRONG AND SAID SO. A bare self-loop
 * child0[i] = i does NOT hang: once `size` in hz_oct_at underflows past 1 to 0,
 * every `x >= lo + size` test is true, the octant index sticks at 7 and the
 * walk moves FORWARD out of the loop. The cycle has to close after that
 * underflow, i.e. through the +7: node 8 pointing back at block 1 sends
 * ni = 1 + 7 = 8 forever. That is the shape below, and the old check passes it
 * because 1 + 8 <= 9. */
static int cycle_hangs(void) {
  hz_onode nodes[9];
  for (int i = 0; i < 9; i++) {
    nodes[i].child0 = -1;
    nodes[i].k2 = 0.0;
  }
  nodes[0].child0 = 1; /* root subdivides */
  nodes[1].child0 = 1; /* drives size 2 -> 1 -> 0 */
  nodes[8].child0 = 1; /* closes the cycle through the stuck octant index 7 */
  hz_octree t;
  t.nodes = nodes;
  t.n = 9;
  t.cap = 9;
  t.log2size = L2;

  struct sigaction sa;
  memset(&sa, 0, sizeof sa);
  sa.sa_handler = on_alarm;
  sigemptyset(&sa.sa_mask);
  if (sigaction(SIGALRM, &sa, NULL) != 0) return -1;
  g_hit = 0;
  if (sigsetjmp(g_jb, 1) == 0) {
    alarm(1);
    double complex v = hz_oct_at(&t, 0, 0, 0);
    alarm(0);
    (void)v;
    return 0; /* returned => no hang */
  }
  return 1;
}

/* -------------------------------------------------------------------------- */

int main(void) {
  printf("=== PREDICTIONS (before any result) ===\n");
  printf("  1 save->load->save is byte-identical, and the loaded tree equals the\n");
  printf("    original cell by cell (the pair that stops `return 1;` passing)\n");
  printf("  2 the same tree built twice gives byte-identical files\n");
  printf("    (weak on its own: fresh mmap pages are zero, so valgrind is the gate)\n");
  printf("  3 a v0 file (no magic) is rejected as E_MAGIC=%d\n", HZ_OCT_E_MAGIC);
  printf("  4 each single-field corruption is rejected by ITS OWN code:\n");
  printf("    version=%d flags=%d endian=%d canary=%d reserved=%d log2size=%d\n", HZ_OCT_E_VERSION,
         HZ_OCT_E_FLAGS, HZ_OCT_E_ENDIAN, HZ_OCT_E_CANARY, HZ_OCT_E_RESERVED, HZ_OCT_E_LOG2SIZE);
  printf("    ncount=%d childneg=%d cycle=%d range=%d align=%d depth=%d twoparent=%d\n",
         HZ_OCT_E_NCOUNT, HZ_OCT_E_CHILD_NEG, HZ_OCT_E_CYCLE, HZ_OCT_E_RANGE, HZ_OCT_E_ALIGN,
         HZ_OCT_E_DEPTH, HZ_OCT_E_TWOPARENT);
  printf("    unreachable=%d, truncation=%d\n", HZ_OCT_E_UNREACHABLE, HZ_OCT_E_IO);
  printf("  5 NEGATIVE CONTROL, predicted to FAIL to return: the cyclic tree hangs\n");
  printf("    hz_oct_at when it is not rejected\n");
  printf("=== RESULTS ===\n");

  hz_octree t;
  check(build(&t) == 0, "build");
  long n = 0;
  unsigned char *b = save_buf(&t, &n);
  check(b != NULL, "save");
  if (b == NULL) return 1;
  printf("  tree: log2size=%d n=%d, file=%ld bytes (header %d + %d*20)\n", t.log2size, t.n, n,
         HZ_OCT_HDRSIZE, t.n);
  check(n == HZ_OCT_HDRSIZE + 20L * t.n, "file length is header + SoA arrays");
  check(t.n == 17, "the test tree is the two-level one the mutations assume");

  /* 1 — roundtrip, and the good file really loads */
  hz_octree t2;
  int rc = load_buf(b, n, &t2);
  check(rc == HZ_OCT_OK, "good file loads");
  if (rc == HZ_OCT_OK) {
    check(same_cells(&t, &t2), "loaded tree equals the original cell by cell");
    long n2 = 0;
    unsigned char *b2 = save_buf(&t2, &n2);
    check(b2 != NULL && n2 == n && memcmp(b, b2, (size_t)n) == 0, "roundtrip is byte-identical");
    free(b2);
    hz_oct_free(&t2);
  }

  /* 2 — determinism: a second, independently allocated copy of the same tree */
  hz_octree t3;
  check(build(&t3) == 0, "build again");
  long n3 = 0;
  unsigned char *b3 = save_buf(&t3, &n3);
  check(b3 != NULL && n3 == n && memcmp(b, b3, (size_t)n) == 0,
        "two independent builds give identical bytes (no padding in the file)");
  free(b3);
  hz_oct_free(&t3);

  /* 3 — the v0 file: raw {log2size, n} header, then the node array */
  {
    long n0 = 8 + (long)sizeof(hz_onode) * t.n;
    unsigned char *b0 = calloc((size_t)n0, 1);
    check(b0 != NULL, "alloc v0");
    if (b0 != NULL) {
      int32_t hdr[2] = {(int32_t)t.log2size, t.n};
      memcpy(b0, hdr, sizeof hdr);
      memcpy(b0 + 8, t.nodes, sizeof(hz_onode) * (size_t)t.n);
      hz_octree tv0;
      int r0 = load_buf(b0, n0, &tv0);
      if (r0 == HZ_OCT_OK) hz_oct_free(&tv0);
      check(r0 == HZ_OCT_E_MAGIC, "v0 file is rejected as E_MAGIC");
      free(b0);
    }
  }

  /* 4 — one corrupted field at a time */
  {
    uint32_t u = 2;
    expect_reject(b, n, HZ_OCT_OFF_VERSION, &u, sizeof u, HZ_OCT_E_VERSION, "version 2");
    u = 1;
    expect_reject(b, n, HZ_OCT_OFF_FLAGS, &u, sizeof u, HZ_OCT_E_FLAGS, "unknown flag bit");
    u = 0x01020304u;
    expect_reject(b, n, HZ_OCT_OFF_ENDIAN, &u, sizeof u, HZ_OCT_E_ENDIAN, "byte order witness");
    double c = 2.0;
    expect_reject(b, n, HZ_OCT_OFF_CANARY, &c, sizeof c, HZ_OCT_E_CANARY, "IEEE754 canary");
  }
  patch_i32(b, n, HZ_OCT_OFF_RESERVED, 1, HZ_OCT_E_RESERVED, "reserved word");
  patch_i32(b, n, HZ_OCT_OFF_LOG2SIZE, HZ_OCT_MAX_LOG2SIZE + 1, HZ_OCT_E_LOG2SIZE, "log2size 21");
  patch_i32(b, n, HZ_OCT_OFF_N, 2, HZ_OCT_E_NCOUNT, "n = 2 is not 1 + 8k");
  patch_i32(b, n, HZ_OCT_OFF_LOG2SIZE, 1, HZ_OCT_E_DEPTH, "depth 2 under log2size 1");

  {
    long c0 = HZ_OCT_HDRSIZE; /* child0[i] lives at c0 + 4*i */
    patch_i32(b, n, c0 + 4 * 5, -2, HZ_OCT_E_CHILD_NEG, "child0 = -2");
    patch_i32(b, n, c0 + 4 * 1, 1, HZ_OCT_E_CYCLE, "child0[1] = 1 (self-loop)");
    patch_i32(b, n, c0 + 4 * 0, 17, HZ_OCT_E_RANGE, "child block past the arena");
    patch_i32(b, n, c0 + 4 * 0, 2, HZ_OCT_E_ALIGN, "child0 not at a block start");
    patch_i32(b, n, c0 + 4 * 2, 9, HZ_OCT_E_TWOPARENT, "block 9 claimed twice");
    patch_i32(b, n, c0 + 4 * 1, -1, HZ_OCT_E_UNREACHABLE, "block 9 orphaned");
  }

  {
    hz_octree tt;
    int r = load_buf(b, n - 1, &tt);
    if (r == HZ_OCT_OK) hz_oct_free(&tt);
    check(r == HZ_OCT_E_IO, "truncated file");
    r = load_buf(b, HZ_OCT_HDRSIZE - 1, &tt);
    if (r == HZ_OCT_OK) hz_oct_free(&tt);
    check(r == HZ_OCT_E_IO, "file shorter than the header");
  }

  /* 5 — negative control with a PREDICTED failure to return */
  {
    int h = cycle_hangs();
    if (h == 0)
      printf("FAIL: negative control: hz_oct_at RETURNED on the cyclic tree, so\n"
             "      E_CYCLE guards nothing\n");
    check(h == 1, "negative control: the cycle hangs hz_oct_at when not rejected");
  }

  free(b);
  hz_oct_free(&t);
  printf("%s: %d/%d\n", g_fail ? "FAILURES" : "ok", g_total - g_fail, g_total);
  return g_fail != 0;
}
