# CLAUDE.md — helmholz project instructions

## What this is
3D renderer via numerical solution of the inhomogeneous Helmholtz equation
(∇²u + k(x)²u = f) in an inhomogeneous medium — wave rendering, not ray tracing.
- Scene: two octrees — one for the medium, a separate one for the light
  sources f(x) (tactically simpler than one shared tree). k is piecewise-constant
  per material (later: linear gradient). A node at ANY level may, instead of
  subdividing, store a linear function (corner interpolation) whose zero set is
  the material boundary — sloped walls/floors as one coarse node, no staircase.
  Large empty regions = one empty node.
  Staging: STAGE 1 = plain hierarchical (cubic) voxels, no boundary
  interpolation — then every operator integral is separable and reduces to a
  small precomputed table of 1D dyadic integrals (no quadrature anywhere).
  STAGE 2 = planar cuts via cube clipping + tetrahedra + exact degree-6
  cubature, swapped in behind the same integrate_cell() interface. Tree node
  format carries an optional boundary-corner payload from day one.
- Solution basis: "potentials" (wavelet-like). 3D potential = tensor product of
  1D potentials; 1D: φ(x)=2−x² on [0,1], (x−2)² on (1,2], mirrored for x<0,
  zero outside [−2,2]. C¹ piecewise-quadratic, integer knots, φ'' piecewise
  constant (∓2) ⇒ Helmholtz-operator matrix elements are analytic; translates
  sum to 4 (partition of unity); φ = 2·(two adjacent quadratic B-splines) ⇒
  refinable — coarse octree levels expand exactly into finer ones.
  The solution is a weighted sum over potentials of many widths at once.
  Deliberately overcomplete (a multilevel generating system, not a basis) ⇒ the
  system is rank-deficient by construction — hence min-norm solvers only (LSQR,
  zgelsd); never form AᴴA. Medium and sources enter through the equation
  residual on the potentials' supports.
- DISCRETIZATION IS SET BY GEOMETRY, NOT BY WAVELENGTH (user directive 07-25;
  this SUPERSEDES the earlier "λ-resolving floor active everywhere"). Element
  size comes from object detail, and at the camera from PIXEL COUNT; λ leaves
  the dimension count entirely. A plain smooth bump cannot do this — it must
  resolve λ even in vacuum (plus Babuška–Sauter pollution). Only a CARRIER
  (Trefftz/PUM/UWVF) basis φ(x/W−n)·e^{ik·x} can: in (d²/dx²+k²)[φ cos] the k²
  term cancels analytically, which is exactly what makes coarse cells legal.
  The carrier is a CONDITION for the directive, not an optimization. Plane-wave
  carriers specifically (not spherical harmonics) because e^{i(kₓx+k_yy+k_zz)}
  FACTORIZES — the tensor separability the assembly depends on survives.
  Unknowns = (elements) × (modes per element). Elements: LOD L = εR plus
  empty-space collapse by the octree (~1e7 for a 10 m room at 1000×1000 px,
  against 3e24 for a λ/8 grid — but that 1e7 assumes empty-space collapse AND
  smooth surfaces; it is conditional, not established).
  Modes: (kL)^(d−1) is a CAPACITY BOUND — how many independent wave modes a
  region of size L can hold — NOT a cost. What is actually needed is the angular
  content of the real field: one carrier toward the camera in free space, two at
  a smooth surface (incident + reflected, Fresnel is analytic, independent of L),
  more only where the surface is structured or multiple scattering is strong.
  The actual requirement is UNKNOWN and is precisely what M9c measures — do not
  quote the bound as if it were the cost.
  For the element→aperture channel the communication-mode count is
  N = L²D²/(λR)² = (εD/λ)², exactly 1 when ε is the diffraction limit λ/D — one
  carrier per element, independent of distance. That is the CAMERA CHANNEL, not
  the local field.
  WHERE λ LEGITIMATELY REMAINS (it left the GRID, not the physics — it is in the
  equation): in the operator, carried analytically by e^{ik·x} at zero grid cost;
  as the evanescence cutoff (sub-λ detail does not radiate ⇒ when to stop
  refining); and as the floor λ/D under ε — ε itself is set by PIXEL COUNT,
  L = R·(field of view / pixels). λ sets no grid spacing anywhere.
  Three regulators replace "λ-floor everywhere": (1) LOD L = εR — CEILING on
  element size, set by aperture angular resolution; (2) refinement at geometry —
  FLOOR on element size, driven by k CONTRAST, not by distance; (3) the medium
  octree stays fine INDEPENDENTLY — sub-element detail enters through the exact
  k² integral over the support, not by refining the basis.
  The camera is an OPERATOR (thin-lens phase mask + Fresnel, both FFT): the
  field is solved on the APERTURE PATCH, no sub-λ layers at the sensor at all.
- CUTS, NOT AVOIDANCE. The basis must be DISCONTINUOUS AT a material boundary,
  and that is achieved by CUTTING the element, never by shrinking it. "Elements
  must not cross a boundary" would bound element size by distance-to-surface and
  destroy LOD (a long wall would force refinement along its whole length).
  A large element stays large and its function is cut by the surface:
  φ·[left]·e^{ik₁x} and φ·[right]·e^{ik₂x} — one extra degree of freedom per
  (element × material region), so the count grows with MATERIAL complexity, not
  with the number of wavelengths. Standard technique (cut-cell, XFEM, Nitsche);
  the geometry for it is already the Stage 2 plan (linear function in the node,
  cube clipping) — which therefore moves onto the critical path, not "later".
  Carrier and cut are PARTNERS, not alternatives: the cut is a local repair at a
  finite number of surfaces (~1e7 in a room), the carrier handles the
  oscillation of the interior (~1e22 wavelengths in the same room). Neither
  substitutes for the other. Mechanically: THE CARRIER IS A SHARED PHASE
  REFERENCE — one continuous function per region, so phase continues across
  element joins automatically instead of being renegotiated; THE CUT IS WHERE
  THE REFERENCE CHANGES. The join rests on the partition of unity
  (Σφ(x/W−n) = 4); break it and spurious amplitude modulation appears (measured:
  30% error in the scale bench).
- LOD ON THE BASIS IS FREE, LOD ON THE GEOMETRY IS FORBIDDEN. Geometric data
  precision and element size are independent: a plane is three doubles, and
  specifying it to 1e-9 m costs exactly what 1 cm costs. An element 1e6 λ wide
  cut by an exactly specified plane is legal — S = ∫k is computed from ANALYTIC
  geometry. The sub-λ requirement (δx ≲ λ/(2π·Δn/n), ~0.3λ at contrast 0.5)
  bites only when the geometry itself is coarsened: voxelised boundaries, or
  classic mesh-simplification LOD. The medium octree must keep exact boundaries
  EVERYWHERE, including far away. This upgrades the third regulator from an
  optimisation to a CORRECTNESS REQUIREMENT.
- Related law, not yet measured quantitatively: the carrier phase must be right
  to ~1 rad across an element, i.e. δk/k ≲ λ/(2π·W). Large elements are bought
  with accurate knowledge of the phase. For gradient media the carrier becomes
  the WKB phase e^{iS}, S' = k(x) — the k² term still cancels exactly, leaving
  one term ∝ k'. Caustics are not fatal there: only the PHASE is borrowed, the
  amplitude is a solved coefficient, so superposed carriers stay finite.
- THE FIELD DOES NOT DEPEND ON THE CAMERA. u solves the same equation whatever
  the observer does; only its REPRESENTATION is camera-dependent (through LOD).
  So a camera move is a RE-PROJECTION of an already-known field into a new basis
  (two-scale restriction/prolongation), not a re-solve; a pure rotation costs
  the camera operator alone. Path tracing has no equivalent — its sampling is
  camera-anchored, so every frame restarts. This is an architectural asset, not
  a late optimization. Consequence for honesty: report cold-start and
  steady-state frame cost SEPARATELY; quoting only the latter is cheating.
  The OPERATOR is incremental too: an entry ⟨B_i, L B_j⟩ depends only on the
  element pair and the medium — no camera in the formula — so surviving elements
  keep their entries verbatim.
  LOD stays ANCHORED TO THE CAMERA: that is what makes the cost OUTPUT-bounded
  (∝ pixels, only logarithmic in scene extent). A scene anchor would make it
  INPUT-bounded and is not an option. Per octave of distance the count is
  independent of distance (c/ε²), so total ≈ (c/ε²)·log2(R_max/R_min); the
  active set must cover the FULL SPHERE and occluded layers, since light arrives
  from outside the frustum and through invisible surfaces — frustum culling is
  not available. Room ≈ 2e8 elements, whole city ≈ 5e8: a city costs 2-3x a room
  for 1e4x the linear size. That ratio is the entire point of LOD.
  Reuse rides on the SAME level hierarchy: LEVEL = SPATIAL SCALE AND UPDATE
  RATE. Coarse levels form a persistent, scene-anchored core computed once
  (global transport); fine levels are the per-frame camera shell, updated
  incrementally. The cascade was a solver device; it doubles as the temporal
  reuse structure.
- COHERENCE IS THE SAME AXIS AS GEOMETRIC PRECISION: both ask how much phase
  information the scene actually determines. Real sources have finite coherence
  length (LED ~3 um, sunlight ~1 um, laser metres), so the strict phase
  requirements apply only WITHIN a coherence volume — which is the only reason
  rendering ordinary light is possible at all. Conversely a single monochromatic
  solve across a scene much larger than that produces interference that does not
  exist in nature: the v1 speckle is a broken assumption, not noise to be
  averaged later. The coherence regime must be a SCENE PARAMETER, not a global
  assumption; strictly coherent treatment is meaningful only where the geometry
  is analytic AND the source narrowband (lens, mirror, grating).
- STATUS: architecture accepted, NUMBERS NOT CONFIRMED. Critical path is
  M9a/M9b/M9c (result/M9_CARRIER_AUDIT.md — audited before running, falsifiers
  fixed in advance; M9a done). Two architectural pieces are NOT specified yet
  and must not be improvised: (i) how continuous LOD sizing maps onto the
  DYADIC level ladder that M1's two-scale relation needs, (ii) what solves the
  system now that there is neither a floor nor a dense cascade head. PLAN.md
  open questions 12-13.
  Everything λ-calibrated in PLAN.md (λ/8 floor, Toeplitz/FFT floor, "domain
  ≥5-6λ", shell thickness and ramp profile, transport percentages, BiCGStab
  verdict) is v1: measured, superseded — DO NOT carry those numbers over.
- Target: GPU. Linear systems (least-squares origin → CG/LSQR-family iterative
  solvers) are the core workload.
- Language: C.

## C CODE QUALITY GATE — MANDATORY, ALWAYS, WHOLE PROJECT
Every `.c`/`.h` written or edited goes through this sequence before it is
"done". Not optional, not per-task — all C in the repo, every time.

1. **Format**: `clang-format -i FILE.c` — config `.clang-format` (LLVM base,
   2-space, no tabs, col 100, K&R).
2. **Static-analysis gate**: `scripts/ccheck.sh FILE.c [...]` — three engines
   in one nix-shell:
   - `gcc -fanalyzer` (path-sensitive: null-deref / UAF / double-free / leak /
     taint) — **gates**;
   - `clang-tidy` (clang-analyzer + bugprone + cert; config `.clang-tidy`) —
     advisory;
   - `cppcheck` (bounds / uninit / leak / realloc / portability) — **gates**.
   Fix until it prints `>>> ccheck: CLEAN`. (gcc-analyzer + cppcheck non-zero =
   must fix; clang-tidy findings are advisory but read them.)
   Known gcc-analyzer false-positive class (diam audit 07-24): it loses the
   capacity↔count link through allocation WRAPPER functions (xmalloc/xcalloc
   style) and reports phantom heap overflows; guards/if-forms do NOT cure it,
   inlining the calloc at the use site does. So: no alloc wrappers in
   hot/indexed-buffer code paths — allocate inline; if a finding looks like
   this class, build the minimal repro before "fixing" real code around it.
3. **Formal layer** (safety-critical / pointer-heavy code: octree
   traversal/mutation, scene-file parsers, anything on untrusted input):
   `scripts/cverify.sh FILE.c [--function NAME] [--unwind N]` — CBMC bounded
   model checking (pointer/bounds/UAF/overflow/leak/div-zero/NaN, no false-neg
   inside the bound). Per-function with `__CPROVER_assume` preconditions. Not
   on every file; on the dangerous ones.
4. **Deadweight** (before commits): `scripts/lean.sh` — whole-project unused
   funcs/members/vars (cross-file; per-file unused already caught by ccheck's
   `-Wall -Wextra`).
5. **Runtime layer** (code that runs on real input): build with
   `-fsanitize=address,undefined -g` and run, or `valgrind --leak-check=full`.
   Numerics extra: run once with `-ffpe-trap`-style checks (feenableexcept on
   FE_INVALID|FE_DIVBYZERO) on a small case to catch NaN sources early.
6. **GPU layer**: kernels (.cu/.cl/shaders) are NOT seen by ccheck — the gate
   covers host C only. For CUDA kernels use `compute-sanitizer` (memcheck +
   racecheck + initcheck) on a small case; keep a CPU reference implementation
   of every kernel and diff results (bitwise for int paths, tolerance for
   float) before trusting GPU output.

## Build / run
- Compiler via nix: `nix-shell -p gcc --run 'gcc -O2 -o build/NAME src/NAME.c -lm'`.
- **LANDMINE (07-25): nixpkgs openblas 0.3.33 LAPACKE zgelsd/zgelss SMASH THE
  STACK on rectangular complex matrices (m != n), both row- and col-major
  (minimal repro verified; square sizes work). Link LAPACK as
  `-llapacke -llapack -lblas` (reference, from nix `lapack blas` packages),
  NOT `$(pkg-config --libs openblas)`. openblas stays only as a header source
  in ccheck. Re-test before ever switching back.**
- Analyzer tools come from nix inside the scripts (gcc, clang-tools, cppcheck,
  cbmc).
- Heavy runs (large grids/solves) — mind memory: `ulimit -v` a sane cap in any
  long-running launcher script rather than trusting the OOM killer.

## Conventions
- All geometry/tree indexing in integer arithmetic where possible; float
  belongs in field values and solver math, not in tree topology.
- No magic thresholds in numerics — tolerances and iteration caps are named
  constants with a comment stating where the number comes from.
