# CLAUDE.md — helmholz project instructions

<!-- ==================== ОТСЮДА ==================== -->
## ВОЛНОВАЯ ЛИНИЯ ОТСТАВЛЕНА С 07-28. НЕ ПОМИНАТЬ ВОВСЕ

**Действует, пока этот блок присутствует в файле. Чтобы вернуть волновую линию —
УДАЛИТЬ БЛОК ЦЕЛИКОМ, от маркера «ОТСЮДА» до маркера «ДОСЮДА». Одна правка, и
всё, что ниже, снова читается как написано.**

Пока блок на месте, **этот запрет ПЕРЕКРЫВАЕТ всё остальное в файле, включая
раздел «ДВА ПРОДУКТА» ниже**:

- **не упоминать волновую линию ни в ответах, ни в докладах, ни в планах** — ни
  Гельмгольца, ни несущий/Трефтц-базис, ни потенциалы φ, ни фазу, ни голограммы,
  ни когерентность, ни `PLAN.md`, ни вехи M-серии;
- **не предлагать** волновых ходов, не сравнивать с ними, не приводить их как
  довод и не «для контекста»;
- **не докладывать состояние обоих продуктов.** Докладывается ТОЛЬКО перенос.
  Требование «называть состояние ОБОИХ» ниже — недействительно;
- если сказанное пользователем можно прочесть и как волновое, и как переносное —
  читать как ПЕРЕНОСНОЕ и не оговариваться.

**ВСЁ ВОЛНОВОЕ ВЫНЕСЕНО В `wave/` 07-28 (по требованию пользователя).** Там
`wave/PLAN.md`, `wave/SLICE_PLAN.md`, `wave/src/`, `wave/tools/`, `wave/tests/` и
свой `wave/Makefile`. Перенос сделан через `git mv`, обратим одной командой.

**ЗАПРЕТ НА РАБОТУ: в `wave/` не входить.** Не собирать, не запускать, не править,
не рефакторить, не «прибирать», не чинить предупреждения, не обновлять под
изменения общего слоя. **Не читать оттуда даже «для контекста»** — именно чтение
и приводило к тому, что волновые доводы просачивались в ответы. Единственное
исключение: пользователь ЯВНО просит что-то из `wave/`.

`wave/` из корневого `Makefile` не собирается СОЗНАТЕЛЬНО; `make test` в корне
гоняет только перенос и общий слой (проверено 07-28: `EXIT=0`). Восстановление
линии — `make -C wave test`, но только после снятия запрета.

**Волновой код НЕ МЁРТВЫЙ, а отложенный** — удалять `wave/` нельзя.

**Оговорки волновой линии в общем слое (Г25, Г40/Г44 про `dmax = UNKNOWN`)
СОХРАНЯЮТСЯ в коде и в тексте** — они там как ворота «fail closed», и снимать их
из-за отставки нельзя: вернётся линия — вернётся и требование.
<!-- ==================== ДОСЮДА ==================== -->

## ЧТО В ПРОЕКТЕ (обновлено 11-09 после чистки; прежняя редакция — в git)

Одна действующая линия — **перенос излучения** (объёмное уравнение переноса,
дискретные ординаты, DG1) с целью «по сцене можно ходить»: интерактивно,
разрушаемо, корректное диффузное переотражение (бюджет ≥5 к/с при 512²,
город и тяжёлые сцены докладываются отдельно).

Документы — всё живое в корне, всё старое в `archive/` (указатель
`archive/README.md`):

- **STRUCTURE.md** — действующее ПРЕДСТАВЛЕНИЕ (киты × размещения, поле на
  поверхности, марш формулей порядка направления; согласовано 09-10). Старый
  путь (конвейер `pfield`/swee3) остаётся ПРОДЮСЕРОМ и ЭТАЛОНОМ до паритета.
  Читать ДО любой структурной работы.
- **PLAN_ELEMENTS.md** — живой append-only ЖУРНАЛ (§833+, эпоха новой
  структуры). Читать С КОНЦА. Дисциплина шага — в его шапке, ОБЯЗАТЕЛЬНА:
  план и предсказания до кода, аудит А-серии до кода, гейт, негативный
  контроль с предсказанным провалом, доклад, аудит исполнения.
- **PLAN_TRANSPORT.md** — вехи T1…T7 (живой хвост T5б → T6 → T7), канон чисел,
  отложенные долги с владельцами.
- **START.md** — точка входа: где мы, порядок дальше, канонические команды,
  действующие ключи, сцены.

В `archive/`: журнал §1–§832 (эпоха старого представления), PLAN_CUT.md
(аудит геометрии реза Г1…Г62 — читать перед ЛЮБОЙ работой с разрезом/DC: там
записаны отвергнутые ходы, и без этого их естественно предложить заново),
new-structure.md (правила и замеры-предшественники STRUCTURE.md: закон формы
§2, свет §6, замеры §9bis), старые артефакты result/.

**Порядок чтения новой сессии:** START.md → хвост PLAN_ELEMENTS.md →
STRUCTURE.md (при структурной работе) → PLAN_TRANSPORT.md (при работе вех).

Живой код: `src/transport/`, `src/cut/`, `src/octree.*`, `src/image.*`,
`src/nstruct/`, `src/p*.c` (пространственное дерево и старый конвейер —
эталонный путь), `tools/`, `tests/`. Удалять код старого пути нельзя до
паритета новой структуры — он эталон.

## What this is
3D renderer via numerical solution of the radiative transport equation
(discrete ordinates, discontinuous Galerkin) — light is propagated as a front
over scene elements, not ray tracing per pixel. Language: C. Solver: sweep
(wavefront traversal) on CPU today; GPU later. Octrees, LOD and the camera
operator are shared machinery. Target scene scale: a city.

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
- **КАРТИНКИ ПИШУТСЯ В `img/`, А НЕ В `build/`.** `build` — только артефакты
  сборки, и мусорить в нём нельзя: он чистится вместе с объектниками, и то, на
  что смотрят глазами, там теряется. Всякий инструмент, выдающий изображение,
  пишет в `img/` по умолчанию (`tools/render3.c` — так).

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
- Битовые слепки/якоря — только при `OMP_NUM_THREADS=1` (А1313).
- Line-based чистилки по grep — только с `git diff` перед коммитом (урок §873:
  авто-чистка съела ~270 строк swee3.c).
