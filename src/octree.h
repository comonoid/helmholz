#ifndef HZ_OCTREE_H
#define HZ_OCTREE_H

/* Stage-1 octree: cubic voxels, piecewise-constant complex k^2 per leaf.
 * Cube [0, 2^log2size)^3 in unit cells. Nodes live in one arena, children
 * referenced by index (8 contiguous) — no per-node malloc (CLAUDE.md: alloc
 * wrappers vs gcc-analyzer, plus locality). Node: child0 < 0 = leaf. */

#include <complex.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
  int32_t child0; /* index of first of 8 children, or -1 for a leaf */
  double complex k2;
} hz_onode;

typedef struct {
  hz_onode *nodes;
  int32_t n, cap;
  int log2size;
} hz_octree;

int hz_oct_init(hz_octree *t, int log2size, double complex background);
void hz_oct_free(hz_octree *t);

/* set k2 on the cell box [lo, hi) (clipped to the cube); returns 0 on ok */
int hz_oct_set_box(hz_octree *t, const int lo[3], const int hi[3], double complex k2);
/* set k2 on cells whose centers lie inside the ball */
int hz_oct_set_ball(hz_octree *t, const double c[3], double r, double complex k2);

/* Index of the leaf containing (x,y,z), or -1 outside the cube. The descent is
 * pure integer, which is what makes it reachable by CBMC (hz_oct_at returns a
 * double complex, and CBMC's symex aborts on the real->complex conversion).
 * It is also the value the cut side table is keyed by: PLAN_CUT.md Г4/Г20 hangs
 * boundary records off an abstract int32 CELL index, and for this tree the cell
 * index is the node index. */
int32_t hz_oct_leaf(const hz_octree *t, int x, int y, int z);
double complex hz_oct_at(const hz_octree *t, int x, int y, int z);

/* visit leaves intersecting [lo, hi): cb gets the CLIPPED box [blo, bhi).
 * Nonzero from cb aborts and is returned. */
typedef int (*hz_oct_cb)(void *ctx, const int blo[3], const int bhi[3], double complex k2);
int hz_oct_visit(const hz_octree *t, const int lo[3], const int hi[3], hz_oct_cb cb, void *ctx);

/* --- serialization (PLAN_CUT.md, Г5 / Р-2) -------------------------------- */

/* Format v1. Header of HZ_OCT_HDRSIZE bytes, then THREE SEPARATE arrays
 * (child0[n], re(k2)[n], im(k2)[n]) — structure-of-arrays, NOT a raw dump of
 * the node array.
 *   why not the raw dump: hz_onode is int32 + 4 padding bytes + double complex,
 *   and nodes handed out by realloc never touch the padding, so a raw fwrite
 *   ships INDETERMINATE bytes — the file is irreproducible and valgrind flags
 *   the write. SoA also decouples the file from the compiler's layout, so the
 *   Г4 side tables can be added without a node-format version bump.
 *   why two witnesses on top of the version number: a version only catches
 *   changes somebody REMEMBERED to bump. double->float, another endianness or
 *   another ABI do not touch it, so the header also carries an integer
 *   byte-order witness and an IEEE754 canary, both compared for exact equality
 *   (Г28). Unknown flag bits are REJECTED (fail closed). */
#define HZ_OCT_HDRSIZE 40
#define HZ_OCT_VERSION 1u
/* byte offsets inside the header — the tests mutate the file by these */
#define HZ_OCT_OFF_MAGIC 0
#define HZ_OCT_OFF_VERSION 8
#define HZ_OCT_OFF_FLAGS 12
#define HZ_OCT_OFF_ENDIAN 16
#define HZ_OCT_OFF_LOG2SIZE 20
#define HZ_OCT_OFF_N 24
#define HZ_OCT_OFF_RESERVED 28
#define HZ_OCT_OFF_CANARY 32

/* 2^20 units per axis: 1e18 leaf cells, far past the arena bound below, and it
 * keeps tree depth inside the uint8 scratch hz_oct_validate uses. */
#define HZ_OCT_MAX_LOG2SIZE 20
/* 2^26 nodes = 1.5 GiB of arena — beyond any sane scene, and it bounds the
 * attacker-controlled allocation (CWE-789). */
#define HZ_OCT_MAX_NODES (1 << 26)

/* Status codes. NOT a plain 0/1: a test that only sees "rejected" cannot show
 * that the file was rejected by the invariant that was corrupted in it, and
 * then the whole falsifier table is satisfied by `return 1;` (Г26). */
typedef enum {
  HZ_OCT_OK = 0,
  HZ_OCT_E_IO = 1,          /* short read / short write */
  HZ_OCT_E_MAGIC = 2,       /* not a tree file (this is the Г5 hole) */
  HZ_OCT_E_VERSION = 3,     /* known kind, unknown version */
  HZ_OCT_E_FLAGS = 4,       /* unknown flag bit set */
  HZ_OCT_E_ENDIAN = 5,      /* integer byte order differs */
  HZ_OCT_E_CANARY = 6,      /* floating-point representation differs */
  HZ_OCT_E_RESERVED = 7,    /* reserved word not zero */
  HZ_OCT_E_LOG2SIZE = 8,    /* cube size out of range */
  HZ_OCT_E_NCOUNT = 9,      /* n out of range, or (n-1) not a multiple of 8 */
  HZ_OCT_E_CHILD_NEG = 10,  /* a leaf is marked ONLY by -1 */
  HZ_OCT_E_CYCLE = 11,      /* child0 <= own index: hz_oct_at would hang (Г27) */
  HZ_OCT_E_RANGE = 12,      /* child block does not fit the arena */
  HZ_OCT_E_ALIGN = 13,      /* child0 is not the start of an 8-block */
  HZ_OCT_E_DEPTH = 14,      /* deeper than log2size: a size-1 node has children */
  HZ_OCT_E_TWOPARENT = 15,  /* a block claimed by two parents */
  HZ_OCT_E_UNREACHABLE = 16 /* node reachable from no parent = junk in the file */
} hz_oct_status;

int hz_oct_save(const hz_octree *t, FILE *f);
int hz_oct_load(hz_octree *t, FILE *f); /* into an uninitialized tree */

/* Pure, integer-only halves of the load — this is what CBMC can reach, since it
 * has no body for fread/FILE* (Г29). */
int hz_oct_hdr_decode(const unsigned char *h, int *log2size, int32_t *n);
/* Every invariant hz_oct_at and the recursive walkers rely on, in ONE forward
 * pass without recursion: children always sit at a HIGHER index than their
 * parent (node_alloc8 only ever appends), so "c > i" alone gives acyclicity,
 * and by the time the pass reaches i every possible parent of i is behind it. */
int hz_oct_validate(const hz_octree *t);

#endif
