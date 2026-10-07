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

Живой код (12-09, после переноса старого пути в `archive/geom/`):
`src/nstruct/` (пирамида, свип), `src/octree.*`, `src/image.*`,
`src/scene_obj.*`, `src/transport/dirs3+quad` (test_bounce), `tools/`
(pgather — продакшн, pref — эталон, swee3/snap3), `tests/`. Старый
геометрический путь (src/p*, src/cut, pfield/pwalk/render3) — в
`archive/geom/`: поглядывать можно, якоря сверять через git-историю.
Геометрия v2 (§914) пишется с чистого листа в `src/geom/`.

## What this is
3D renderer via numerical solution of the radiative transport equation
(discrete ordinates, discontinuous Galerkin) — light is propagated as a front
over scene elements, not ray tracing per pixel. Language: C. Solver: sweep
(wavefront traversal) on CPU today; GPU later. Octrees, LOD and the camera
operator are shared machinery. Target scene scale: a city.

## C CODE QUALITY GATE — MANDATORY, ALWAYS, WHOLE PROJECT
Every `.c`/`.h` written or edited goes through this sequence before it is
"done". Not optional, not per-task — all C in the repo, every time.
**Процедура целиком, с командами и контрактом триажа — в скилле `.dsh/skills/c-gate`.**

Исполняемая часть — `scripts/cgate.sh` (ступени S0…S4, порядок по цене, ранний
выход, кэш, база разобранных находок, отчёт `build/cgate/report.tsv`):

    make check-staged   S0–S1 по изменённым файлам (то же гоняет pre-commit)
    make check-fast     S0–S2 по всему живому слою (~1 мин, S2 — дорогая)
    make check          S0–S4 по всему живому слою — главный статический гейт
    make gate           всё: статика + lean + CBMC + санитайзеры + FPE + valgrind

1. **Format** (S0): `clang-format -i FILE.c` — config `.clang-format` (LLVM
   base, 2-space, no tabs, col 100, K&R). Проверка не мутирует файл:
   `clang-format --dry-run --Werror`.
2. **Static-analysis gate**:
   - **S1** `gcc -O2 -c -Werror` с общим набором `scripts/cflags.mk`
     (включая `-Wjump-misses-init` — класс «goto перескакивает инициализатор»;
     07-10 он дал настоящий `free()` по мусору в `src/nstruct/sweep.c:3338`) —
     **gates**;
   - **S2** `cppcheck` (bounds / uninit / leak / realloc / portability) —
     **gates**;
   - **S3** `clang-tidy` (clang-analyzer + bugprone + cert; config `.clang-tidy`,
     флаги из `build/compile_commands.json`) — гейтит **по политике `tidy_blocks`
     в `scripts/cgate.sh`**, а не по коду возврата (при 46 находках rc=0;
     система шапок даёт ~4 с фиксированной цены на файл). Блокируют классы с
     прямым дефектом (`clang-analyzer-core.*`, `unix.Stream`, утечка при
     `realloc`, целочисленное деление в FP-контексте, сравнение объектного
     представления, смещённое расширение; `unix.Malloc` — для `src/**`). В
     ОТЧЁТ (`build/cgate/advisory/`, 224 находки на 07-10) уходят гигиена
     (`cert-err34-c` на `atoi` в разборе argv, `misc-include-cleaner`) и
     «optin»-классы; `bugprone-implicit-widening` — полезен, но требует планового
     перевода индексации на 64 бита, поэтому до решения он в отчёте. Отчёт
     печатается сводкой по классам — «не блокирует» не значит «не существует»;
   - **S4** `gcc -O0 -fanalyzer` — **gates**, но идёт ПОСЛЕДНИМ и под таймаутом
     (44 с на `sweep.c`; на дешёвой ступени красный файл до него не доходит).
   Известный FP-класс gcc-analyzer (аудит diam 07-24): он теряет связь
   «ёмкость↔счётчик» через ОБЁРТКИ выделения (xmalloc/xcalloc) и рапортует
   фантомные переполнения кучи; защитные условия это НЕ лечат, лечит calloc на
   месте использования. Отсюда правило: в горячих/индексных буферах обёрток
   выделения нет — выделять на месте; похоже на этот класс — сначала минимальный
   репро, потом правка настоящего кода.
   **База и ратчет**: `scripts/cgate-baseline.txt` — разобранные находки с
   вердиктом (`FP (механизм …)` / `ОСОЗНАННО (…§…)`). Гейт валит только НОВОЕ.
   Вердикт «НОВОЕ — требует вердикта» — долг, а не разрешение.
3. **Formal layer** (safety-critical / pointer-heavy code: octree
   traversal/mutation, scene-file parsers, anything on untrusted input):
   `make cbmc` или `scripts/cverify.sh FILE.c [--function NAME] [--unwind N]` —
   CBMC bounded model checking (pointer/bounds/UAF/overflow/leak/div-zero/NaN, no
   false-neg inside the bound). Per-function with `__CPROVER_assume`
   preconditions. Not on every file; on the dangerous ones. **07-10: оснастка
   `tests/cbmc_sceneobj.c` два месяца не компилировалась** (сигнатура
   `parse_fvert` разошлась) — то есть слой существовал только на бумаге; стенд
   починен и теперь проверяется воротами: `make cbmc` даёт 4591 свойство
   (разборщик, VERIFICATION SUCCESSFUL) и 1645 свойств (октодерево,
   VERIFICATION SUCCESSFUL), включая несущее свойство текстурного индекса.
   Оснастке октодерева обязателен `src/octree.c` в списке файлов: без него у
   CBMC нет тела `hz_oct_validate`, и свойства проваливаются ЛОЖНО.
4. **Deadweight** (before commits): `make lean` (`scripts/lean.sh`) —
   whole-project unused funcs/members/vars (cross-file; per-file unused already
   caught by S1's `-Wall -Wextra`). Список файлов — из `git ls-files`, не
   рукописный (рукописный список пропускал `tools/pgather.c` и `src/transport/`).
5. **Runtime layer**: `make test-asan` (ASan+UBSan, `detect_leaks=1`),
   `make test-uninit` (`-ftrivial-auto-var-init=pattern` — чтение
   неинициализированного), `make check-int` (clang
   `-fsanitize=integer,implicit-conversion` по целочисленной топологии),
   `make check-fpe` (FE_INVALID|FE_DIVBYZERO через `tests/fpetrap.c` +
   `LD_PRELOAD`, без правок продуктового кода), `make valgrind`.
6. **GPU layer**: kernels (.cu/.cl/shaders) are NOT seen by the gate — the gate
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
- **ДИРЕКТИВА ПОЛЬЗОВАТЕЛЯ (01-10): ВСЕ примеры и бенчмарки рендера —
  в разрешении FullHD (1920×1080)**, если явно не сказано иное. Мелкие
  разрешения (320×240 и т.п.) допустимы только для диагностических
  приборов, не для отчётных чисел. Цель-рамка «≥5 к/с» измеряется на
  FullHD.

## ДИРЕКТИВА ПОЛЬЗОВАТЕЛЯ (11-09, изм. тем же днём): ДЕЛАТЬ ВСЕ ОПТИМИЗАЦИИ; ВОЗМОЖНОСТИ ПАРАЛЛЕЛИЗМА И GPU ОСТАЮТСЯ

**УТОЧНЕНИЕ ПОЛЬЗОВАТЕЛЯ (01-10, ПОСЛЕ ПОВТОРНЫХ ЗАМЕЧАНИЙ — ЖЁСТКО):
ВСЕ прогоны, бенчмарки и отчётные числа — ТОЛЬКО `OMP_NUM_THREADS=1`.
НЕ предлагать, НЕ запускать и НЕ публиковать многопоточные замеры
(«16 ядер» и т.п.) — даже как вариант. Параллелизм и GPU существуют
ИСКЛЮЧИТЕЛЬНО как проектное требование к новому коду (пункты 1–3
ниже): hot path stateless, данные SoA — чтобы ПОТОМ можно было
распараллелить на 16 потоков и перенести на видеокарту. Но режимом
измерения это никогда не становится.**

Цель — ВСЕ оптимизации (указание: «запаса по скорости ещё чуть ли не в 5
раз»). Параллелизм blk-сбора ВКЛЮЧЁН (§886); требования к любому новому коду:

1. **Горячий путь обязан оставаться stateless на уровне независимой единицы
   работы** (луч/строка кадра/трубка): никакой мутабельной памяти, общей
   между единицами — распараллеливание обязано сводиться к разбиению
   диапазона (прецеденты: §878 RAM-сбор ×8–14, §884 blk-сбор ×7.5 — оба
   включаются одной строкой `#pragma omp`).
2. **Данные — SoA/блочно, без указателей наружу блока** (формат HBLK §882
   этому соответствует): GPU-порт = загрузка блока в видеопамять как есть.
3. Однопоточный прогон (OMP_NUM_THREADS=1) остаётся обязательной точкой
   отсчёта и НК; переключатель параллелизма не удалять и не ломать
   (проверка битовости OMP=1 vs N — обязательная гейт-процедура).

## ДИРЕКТИВА ПОЛЬЗОВАТЕЛЯ (27-09, после §912): ПРОТОКОЛ НЕ РАЗДУВАТЬ; КАЖДЫЙ ПРИБОР ОПРАВДАТЬ

Повод: сессии §911–§912 потратили часы на якорные/аудитные прогоны и
эталонные построения при нулевом движении продукта («как сдрейфовали от
0.1 с до 40 минут» — стакан диагностических режимов: OMP=1 × 128 ординат
× it до сходимости × clip без эффекта).

1. **НИ ОДНОГО прогона длиннее ~30 с без явного разрешения пользователя.**
   Перед запуском назвать: что меряем, что ждём, сколько займёт, что будет
   сделано иначе в зависимости от результата. Нет ответа — не запускать.
2. **Каждая новая «эталонная штука» (инструмент, лестница K, переиздание
   якорей, кадр-серия) — ДО создания записывается: какую product-задачу
   закрывает, ожидаемая польза, цена. Пользователь вправе потребовать
   обоснование и отклонить. Польза «подтвердить уже записанное число»
   обоснованием НЕ считается.**
3. Продукт-цель главнее метрологии: интерактивный проход по сцене —
   мерило; если раунд не приближает кадр/с — он не запускается.
