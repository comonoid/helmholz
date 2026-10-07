#!/usr/bin/env bash
# cgate.sh — СТУПЕНЧАТЫЙ C-гейт живого слоя (src, tools, tests).
#
# Почему ступени, а не «три движка по кругу» (как было в ccheck.sh):
# замерено 07-10 на 437/660/4408 строк — clang-tidy 5.4/4.5/15.6 с, gcc
# -fanalyzer 0.9/1.7/44.2 с, при этом gcc -O2 0.3/0.5/3.2 с, cppcheck
# 0.2/2.4/3.0 с. Прежний порядок запускал САМЫЙ дорогой движок первым, без
# раннего выхода: на файл тратилось 44 с прежде, чем выяснялось, что он не
# проходит дешёвую ступень. Здесь порядок обратный цене и есть ранний выход.
#
# Ступени:
#   S0 clang-format --dry-run --Werror      ~10 мс   формат (НЕ мутирует файл)
#   S1 gcc -O2 -c -Werror (+ -Wjump-misses-init)     предупреждения; среди них
#       класс «goto перескакивает инициализатор» — именно он давал free() по
#       мусору в src/nstruct/sweep.c:3338 (ловля 07-10)
#   S2 cppcheck                             независимый парсер: bounds/uninit/leak
#   S3 clang-tidy -p build                  вторая модель (флаги из compile_commands)
#   S4 gcc -O0 -fanalyzer -c                path-sensitive, самый дорогой, последним
#
# НЕ здесь (не про один файл; вызываются целями Makefile):
#   S5 scripts/lean.sh (межфайловая мелочь)  make lean
#   S6 scripts/cverify.sh (CBMC)             make cbmc
#   S7 санитайзеры/FPE/valgrind              make test-asan | check-int | check-fpe | valgrind
#
# База и ратчет: находка, записанная в scripts/cgate-baseline.txt, гейт НЕ
# валит (она разобрана и имеет вердикт). Новых находок гейт не прощает. Так
# «>>> cgate: CLEAN» снова означает что-то: без базы гейт был красен всегда, и
# новую настоящую ошибку было не отличить от стоячего шума.
#
# Usage:
#   scripts/cgate.sh [--all | --staged | --changed=REF | ФАЙЛ...]
#                    [--stage=0,1,2] [--strict] [--digest]
#                    [--jobs=N] [--report=DIR] [--baseline=FILE]
#                    [--update-baseline] [--no-cache] [--quiet]
# По умолчанию: ступени 0,1,2 по всем файлам живого слоя.
# --strict   : база не учитывается, валит любая находка.
# --digest   : напечатать ВСЕ находки таблицей (файл, ступень, класс, текст) —
#              это вход для триажа субагентом.
# --update-baseline: дописать новые находки в базу с вердиктом «НОВОЕ».
#
# Артефакты: build/cgate/report.tsv (итог), logs/ (полные логи движков),
# findings/ (новые), all/ (все), stamps/ (кэш по хэшу файла+заголовков+пина).
#
# Exit: 0 — новых находок нет; 1 — есть новые или ступень не отработала; 2 — ошибка вызова.
set -uo pipefail

HZ="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# ---------------------------------------------------------------- аргументы --
# Умолчания ставятся через `:=` (только если переменной НЕТ в окружении):
# воркер — отдельный процесс, и он обязан получить ступени/каталог отчёта/базу
# от родителя. Жёсткое присваивание здесь УЖЕ съело это однажды: воркер
# подставлял свои «0,1,2», и ступени S3/S4 не запускались ни разу, а гейт при
# этом печатал «CLEAN» (ловля 07-10 — ровно тот молчаливый отказ, ради
# которого гейт и делался).
: "${STAGES:=0,1,2}"
: "${JOBS:=$(nproc 2>/dev/null || echo 4)}"
: "${REPORT:=$HZ/build/cgate}"
: "${BASELINE:=$HZ/scripts/cgate-baseline.txt}"
MODE=""
REF=""
CACHE=1
STRICT=0
DIGEST=0
QUIET=0
UPDATE_BASELINE=0
WORKER_FILE=""
FILES=()
# Аргументы сохраняются ДО разбора: перезапуск в пине (ниже) обязан получить их
# все — иначе --stage/--quiet/список файлов молча теряются (ловля 07-10).
ORIG_ARGS=("$@")
FANALYZER_TIMEOUT="${CGATE_FANALYZER_TIMEOUT:-300}"
FANALYZER_MEM_KB="${CGATE_FANALYZER_MEM_KB:-8000000}"
STAGE_TIMEOUT="${CGATE_STAGE_TIMEOUT:-180}"

usage() { sed -n '2,45p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 2; }

while [ $# -gt 0 ]; do
  case "$1" in
    --stage=*) STAGES="${1#*=}" ;;
    --all) MODE=all ;;
    --staged) MODE=staged ;;
    --changed=*) MODE=changed; REF="${1#*=}" ;;
    --jobs=*) JOBS="${1#*=}" ;;
    --report=*) REPORT="${1#*=}" ;;
    --baseline=*) BASELINE="${1#*=}" ;;
    --no-cache) CACHE=0 ;;
    --strict) STRICT=1 ;;
    --digest) DIGEST=1 ;;
    --quiet) QUIET=1 ;;
    --update-baseline) UPDATE_BASELINE=1 ;;
    --worker) WORKER_FILE="${2:-}"; shift ;;
    -h | --help) usage ;;
    -*) echo "cgate: неизвестный ключ $1" >&2; usage ;;
    *) FILES+=("$1") ;;
  esac
  shift
done

# ------------------------------------------------------ пин инструментов -----
# Гейт всегда работает в ОДНОМ шелле (shell.nix): версии движков версие-
# зависимы, а база и кэш ступеней без пина бессмысленны.
if [ "${CGATE_IN_SHELL:-}" != "1" ]; then
  export CGATE_IN_SHELL=1
  exec nix-shell "$HZ/shell.nix" --run "exec bash $(printf '%q ' "$HZ/scripts/cgate.sh" "${ORIG_ARGS[@]}")"
fi

# shellcheck source=scripts/cflags.sh
. "$HZ/scripts/cflags.sh"
hz_load_cflags "$HZ/scripts/cflags.mk" || exit 2
export HZ HZ_INC HZ_BASE HZ_LIBS HZ_WARN HZ_WARN_GCC ARCH CGATE_IN_SHELL

# ------------------------------------------------------------- служебное -----
mangle() { printf '%s' "$1" | sed 's#^\./##; s#[/.]#_#g'; }

tool_versions() {
  printf 'gcc\t%s\n' "$(gcc -dumpfullversion 2>/dev/null || gcc --version 2>/dev/null | head -1)"
  printf 'clang\t%s\n' "$(clang-tidy --version 2>/dev/null | sed -n 's/.*LLVM version \([0-9][0-9.]*\).*/\1/p' | head -1)"
  printf 'cppcheck\t%s\n' "$(cppcheck --version 2>/dev/null | awk '{print $2}')"
  printf 'clang-format\t%s\n' "$(clang-format --version 2>/dev/null | awk '{print $NF}')"
}

live_paths() {
  # tests/repro/ — ФИКСТУРЫ базы гейта: урезанные файлы той же формы, что и
  # места находок, часть намеренно не форматирована и зовёт __CPROVER_assume.
  # Это НЕ код проекта, и гейт их не судит. Исключение явное, потому что
  # `git ls-files 'tests/*.c'` их берёт: в pathspec git `*` пересекает `/`
  # (проверено 07-10 — без этого фильтра хук заблокировал коммит на 15 находках).
  git ls-files 'src/*.c' 'src/*.h' 'tools/*.c' 'tools/*.h' 'tests/*.c' 'tests/*.h' |
    grep -v '^tests/repro/'
}

tool_key() {
  {
    tool_versions
    cat "$HZ/scripts/cflags.mk"
    cat "$HZ/scripts/cgate.sh"
  } | sha256sum | cut -c1-12
}

header_hash() {
  local p
  for p in $(live_paths | grep '\.h$'); do
    printf '%s %s\n' "$p" "$(sha256sum "$p" | cut -c1-16)"
  done | sha256sum | cut -c1-16
}

# ------------------------------------------------------------ нормализация ---
# Ключ находки = (файл, ступень, класс, текст) БЕЗ номера строки: иначе сдвиг
# строк обнулял бы базу на каждой правке. Текст сохраняется дословно —
# потерять имена переменных значит потерять различие между находками.
#
# Класс берётся из ХВОСТА В КВАДРАТНЫХ СКОБКАХ и бывает двух видов:
#   gcc        [-Wmaybe-uninitialized], [-Werror=…], [-Wanalyzer-…]
#   clang-tidy [cert-err34-c], [bugprone-…,cert-…] — БЕЗ ведущего -W
# Прежняя версия признавала только [-W…] и потому МОЛЧА ВЫБРАСЫВАЛА все находки
# clang-tidy (S3 показывал 0 при 3571 предупреждении в логе). Ловля 07-10:
# любая нормализация обязана пропускать неизвестный класс ГРОМКО, а не тихо.
norm_diag() {
  awk '
    /: (error|warning): / {
      line = $0; cls = ""
      if (match(line, /\[[A-Za-z0-9_.,+=:-]+\]$/)) {
        cls = substr(line, RSTART + 1, RLENGTH - 2)
        sub(/ \[[A-Za-z0-9_.,+=:-]+\]$/, "", line)
      }
      sub(/^[^ ]*:[0-9]+:[0-9]+: /, "", line)
      sub(/^(error|warning): /, "", line)
      if (cls == "") { printf "БЕЗ-КЛАССА\t%s\n", line }
      else printf "%s\t%s\n", cls, line
    }'
}

norm_cppcheck() {
  awk '
    /: (error|warning|style|performance|portability): / {
      line = $0; cls = ""
      if (match(line, /\[[a-zA-Z0-9_]+\]$/)) {
        cls = substr(line, RSTART + 1, RLENGTH - 2)
        sub(/ \[[a-zA-Z0-9_]+\]$/, "", line)
      }
      sub(/^\[[^]]*\]: /, "", line)
      sub(/^[^ ]*:[0-9]+:[0-9]+: /, "", line)
      sub(/^\([a-z]+\) /, "", line)
      sub(/^(error|warning|style|performance|portability): /, "", line)
      if (cls != "") printf "%s\t%s\n", cls, line
    }'
}

declare -A BASE
load_baseline() {
  BASE=()
  [ -r "$BASELINE" ] || return 0
  local f st cls msg verdict
  while IFS=$'\t' read -r f st cls msg verdict; do
    case "$f" in '' | '#'*) continue ;; esac
    BASE["$f|$st|$cls|$msg"]="$verdict"
  done <"$BASELINE"
}

# ---------------------------------------------------------------- ступени ----
stage_cmd() {
  local st="$1" f="$2" log="$3"
  case "$st" in
    0)
      clang-format --dry-run --Werror "$f" >"$log" 2>&1
      ;;
    1)
      case "$f" in
        *.h) gcc $HZ_INC $HZ_BASE -O2 $ARCH $HZ_WARN $HZ_WARN_GCC -Werror -fsyntax-only -x c-header "$f" >"$log" 2>&1 ;;
        *) gcc $HZ_INC $HZ_BASE -O2 $ARCH $HZ_WARN $HZ_WARN_GCC -Werror -c "$f" -o /dev/null >"$log" 2>&1 ;;
      esac
      ;;
    2)
      local cppinc="-I$HZ/src" lang=""
      [ -n "${HZ_OMP_INC:-}" ] && cppinc="$cppinc -I$HZ_OMP_INC"
      case "$f" in *.h) lang="--language=c" ;; esac
      # Потолок на ступень: cppcheck на tools/pgather.c считается ~38 с
      # (замерено 07-10, все конфигурации одинаково — цена файла, не флагов), а
      # на патологическом файле может расти без предела. Потолок превращает это
      # в ГЕЙТ-СБОЙ, а не в «гейт молчит».
      timeout "$STAGE_TIMEOUT" cppcheck --enable=warning,performance,portability --std=c11 \
        --inline-suppr --suppress=missingIncludeSystem --suppress=normalCheckLevelMaxBranches \
        --error-exitcode=1 $lang $cppinc "$f" >"$log" 2>&1
      [ $? -eq 124 ] && echo "ГЕЙТ-СБОЙ: ступень S2 не уложилась в ${STAGE_TIMEOUT}s" >>"$log"
      ;;
    3)
      timeout "$STAGE_TIMEOUT" clang-tidy --quiet -p "$HZ/build" "$f" >"$log" 2>&1
      [ $? -eq 124 ] && echo "ГЕЙТ-СБОЙ: ступень S3 не уложилась в ${STAGE_TIMEOUT}s" >>"$log"
      ;;
    4)
      (
        ulimit -v "$FANALYZER_MEM_KB" 2>/dev/null || true
        NIX_HARDENING_ENABLE='' timeout "$FANALYZER_TIMEOUT" \
          gcc $HZ_INC $HZ_BASE -O0 $HZ_WARN -fanalyzer -c "$f" -o /dev/null
      ) >"$log" 2>&1
      local rc=$?
      [ "$rc" -eq 124 ] && echo "ГЕЙТ-СБОЙ: ступень S4 не уложилась в ${FANALYZER_TIMEOUT}s" >>"$log"
      ;;
  esac
  return 0
}

stages_for() {
  case "$1" in *.h) echo "0 1 2" ;; *) echo "0 1 2 3 4" ;; esac
}

# S3: какие классы clang-tidy БЛОКИРУЮТ, а какие уходят в отчёт.
# Решение 07-10 по замеру на 44 файлах (S3 был включён и разобран впервые):
#   блокируют — классы с прямым дефектом: разыменование/UB/неинициализированное
#     (clang-analyzer-core.*), FILE*-потоки (unix.Stream), утечка при realloc,
#     целочисленное деление в FP-контексте, сравнение объектного представления,
#     смещённое расширение;
#   в отчёт — гигиена (cert-err34-c на atoi в разборе argv, misc-include-cleaner)
#     и «optin»-классы анализатора (portability/performance/taint), плюс
#     bugprone-implicit-widening* — класс признан полезным, но требует планового
#     перевода индексной арифметики на 64 бита; до этого решения он в отчёте,
#     а не в гейте (иначе гейт красный без действия).
# unix.Malloc: для src/** БЛОКИРУЕТ (солвер живёт долго), для tools/tests —
# отчёт: утечка на пути отказа умирает вместе с короткоживущим процессом.
# implicit-widening: для src/** БЛОКИРУЕТ (проверено 07-10: в продакшне такие
# сайты существуют и починены — 6·r/3·r в региональных массивах, смещения
# изображений), для tools/tests — отчёт: там индексы-литералы или счётчики ≤512.
tidy_blocks() { # $1 = файл, $2 = класс (clang-tidy может отдать список через запятую)
  case "$2" in
    *clang-analyzer-core.* | *clang-analyzer-unix.Stream* | *bugprone-suspicious-realloc-usage* | \
      *bugprone-integer-division* | *bugprone-suspicious-memory-comparison* | \
      *bugprone-misplaced-widening-cast*) return 0 ;;
  esac
  case "$1" in
    src/*)
      case "$2" in
        *clang-analyzer-unix.Malloc* | *bugprone-implicit-widening-of-multiplication-result*) return 0 ;;
      esac
      ;;
  esac
  return 1
}

# ------------------------------------------------------------------ воркер ---
worker() {
  local f="$WORKER_FILE"
  [ -n "$f" ] || { echo "cgate: --worker без файла" >&2; exit 2; }
  load_baseline

  local -a want=()
  local s
  for s in $(stages_for "$f"); do
    case ",$STAGES," in *",$s,"*) want+=("$s") ;; esac
  done

  local m rows allfd findings
  m="$(mangle "$f")"
  rows="$REPORT/rows/$m.tsv"
  allfd="$REPORT/all/$m.tsv"
  findings="$REPORT/findings/$m.tsv"
  advisory="$REPORT/advisory/$m.tsv"
  : >"$rows"; : >"$allfd"; : >"$findings"; : >"$advisory"

  if [ "${#want[@]}" -eq 0 ]; then
    printf '%s\tS-ПРОПУСК\tSKIP\t0\t0\t0\t0\tни одна ступень не подходит файлу\n' "$f" >>"$rows"
    exit 0
  fi

  local fhash st log t0 t1 ms raw new base status cls msg key line
  fhash="$(sha256sum "$f" | cut -c1-16)"
  local hashkey="${HDRHASH:-x}-$fhash"

  for st in "${want[@]}"; do
    log="$REPORT/logs/$m.S$st.log"
    local stamp="$REPORT/stamps/$TOOLKEY/S$st/$hashkey"
    : >"$allfd.raw"
    if [ "$CACHE" = 1 ] && [ -s "$stamp.row" ]; then
      # Кэш хранит НАХОДКИ ступени (СЫРЫЕ), а не вердикт: вердикт зависит от базы
      # и от политики блокировки S3, а они правятся без смены ключа кэша (ловля
      # 07-10: кэшированная ступень продолжала считаться NEW после внесения
      # находки в базу). Фильтр и классификация выполняются всегда заново.
      [ -s "$stamp.raw" ] && cat "$stamp.raw" >"$allfd.raw"
      ms="$(cut -f7 "$stamp.row")"
      toolfail="$(cut -f3 "$stamp.row")"
    else
      t0=$(date +%s%N)
      stage_cmd "$st" "$f" "$log"
      t1=$(date +%s%N)
      ms=$(((t1 - t0) / 1000000))
      toolfail=""
      raw=""
      case "$st" in
        0 | 1 | 3 | 4) raw="$(norm_diag <"$log")" ;;
        2) raw="$(norm_cppcheck <"$log")" ;;
      esac
      if [ -n "$raw" ]; then
        while IFS=$'\t' read -r cls msg; do
          [ -n "$cls" ] || continue
          printf '%s\tS%s\t%s\t%s\t1\n' "$f" "$st" "$cls" "$msg" >>"$allfd.raw"
        done < <(printf '%s\n' "$raw" | sort -u)
      fi
      grep -q 'ГЕЙТ-СБОЙ' "$log" 2>/dev/null && toolfail="TOOLFAIL"
    fi

    # Политика S3: не блокирующие классы уходят в ОТЧЁТ (см. tidy_blocks) —
    # видимы, но гейт не валят; молча не выбрасывается ничего.
    : >"$allfd.stage"
    if [ -s "$allfd.raw" ]; then
      while IFS=$'\t' read -r _f _st cls msg _one; do
        [ -n "$cls" ] || continue
        if [ "$st" = 3 ] && ! tidy_blocks "$f" "$cls"; then
          printf '%s\tS3-СОВЕТ\t%s\t%s\n' "$f" "$cls" "$msg" >>"$advisory"
          continue
        fi
        printf '%s\tS%s\t%s\t%s\t1\n' "$f" "$st" "$cls" "$msg" >>"$allfd.stage"
      done <"$allfd.raw"
    fi

    # Классификация — всегда СВЕЖАЯ, по текущей базе; материал ступени пишется в
    # свои файлы (не «всё, что накопилось»), иначе кэш-попадание дублировало бы
    # в отчёте строки предыдущих ступеней.
    : >"$findings.stage"
    new=0
    base=0
    local rawcount=0
    if [ -s "$allfd.stage" ]; then
      while IFS=$'\t' read -r _f _st cls msg _one; do
        [ -n "$cls" ] || continue
        rawcount=$((rawcount + 1))
        key="$f|S$st|$cls|$msg"
        if [ -n "${BASE[$key]:-}" ] && [ "$STRICT" = 0 ]; then
          base=$((base + 1))
        else
          new=$((new + 1))
          printf '%s\tS%s\t%s\t%s\t1\n' "$f" "$st" "$cls" "$msg" >>"$findings.stage"
        fi
      done <"$allfd.stage"
    fi

    if [ "$toolfail" = "TOOLFAIL" ]; then
      status="TOOLFAIL"
    elif [ "$new" -gt 0 ]; then
      status="NEW"
    elif [ "$base" -gt 0 ]; then
      status="BASELINED"
    else
      status="OK"
    fi

    row="$(printf '%s\tS%s\t%s\t%s\t%s\t%s\t%s' "$f" "$st" "$status" "$new" "$base" "$rawcount" "$ms")"
    printf '%s\n' "$row" >>"$rows"
    cat "$allfd.stage" >>"$allfd"
    cat "$findings.stage" >>"$findings"

    if [ "$CACHE" = 1 ]; then
      mkdir -p "$(dirname "$stamp")"
      printf '%s\n' "$row" >"$stamp.row"
      cp "$allfd.raw" "$stamp.raw"
    fi
    rm -f "$allfd.stage" "$findings.stage"

    # РАННИЙ ВЫХОД: дешёвая ступень нашла НОВОЕ — дорогие не запускаем.
    if [ "$status" = "NEW" ] || [ "$status" = "TOOLFAIL" ]; then
      local rest=""
      for s in "${want[@]}"; do
        [ "$s" -gt "$st" ] && rest="$rest S$s"
      done
      [ -n "$rest" ] && printf '%s\tS-ПРОПУСК\tSKIP\t0\t0\t0\t0\tне запускались:%s (S%s красная)\n' "$f" "$rest" "$st" >>"$rows"
      break
    fi
  done
  exit 0
}

# --------------------------------------------------------------- родитель ----
if [ -n "$WORKER_FILE" ]; then
  : "${TOOLKEY:?воркер запускается только из родителя}"
  worker
fi

if [ -z "$MODE" ] && [ "${#FILES[@]}" -eq 0 ]; then MODE=all; fi

case "$MODE" in
  all) mapfile -t FILES < <(live_paths) ;;
  staged) mapfile -t FILES < <({ git diff --name-only --cached --diff-filter=ACMR; git diff --name-only --diff-filter=ACMR; } | sort -u) ;;
  changed) mapfile -t FILES < <(git diff --name-only --diff-filter=ACMR "$REF") ;;
esac

sel=()
for f in "${FILES[@]}"; do
  case "$f" in tests/repro/*) continue ;; esac # фикстуры базы, не код (см. live_paths)
  case "$f" in
    src/*.c | src/*.h | tools/*.c | tools/*.h | tests/*.c | tests/*.h) [ -f "$f" ] && sel+=("$f") ;;
  esac
done
FILES=("${sel[@]:-}")

if [ "${#FILES[@]}" -eq 0 ]; then
  echo "cgate: файлов для проверки нет (режим ${MODE:-список})"
  exit 0
fi

STAGES="$(printf '%s' "$STAGES" | tr -d ' ')"

case ",$STAGES," in
  *,3,*) bash "$HZ/scripts/gen_compile_db.sh" >/dev/null || echo "cgate: compile_commands.json не собрался" ;;
esac

TOOLKEY="$(tool_key)"
HDRHASH="$(header_hash)"
export TOOLKEY HDRHASH REPORT BASELINE STAGES CACHE STRICT FANALYZER_TIMEOUT FANALYZER_MEM_KB STAGE_TIMEOUT

mkdir -p "$REPORT"/{logs,rows,all,findings,advisory,stamps}
# БЛОКИРОВКА: два прогона гейта в одном каталоге отчёта затирают друг другу
# rows/all/advisory и печатают ЛОЖНУЮ сводку (ловля 07-10: параллельный прогон
# по двум файлам обнулил отчёт полного прогона, а тот отчитался «46 советов»
# вместо реальных). Ждём освобождения, а не портим отчёт.
exec 9>"$REPORT/.lock"
if ! flock -w 3600 9; then
  echo "cgate: не дождался блокировки $REPORT/.lock (другой прогон больше часа)" >&2
  exit 2
fi
rm -f "$REPORT"/rows/*.tsv "$REPORT"/all/*.tsv "$REPORT"/findings/*.tsv "$REPORT"/advisory/*.tsv 2>/dev/null

if [ "$QUIET" != 1 ]; then
  echo "=== cgate: ${#FILES[@]} файлов, ступени S${STAGES//,/, S}, потоков $JOBS"
  echo "=== инструменты: $(tool_versions | tr '\n' ' ')"
  echo "=== пин/кэш: $TOOLKEY (заголовки $HDRHASH)"
fi

printf '%s\n' "${FILES[@]}" | xargs -P "$JOBS" -r -I{} bash "$HZ/scripts/cgate.sh" --worker {}

cat "$REPORT"/rows/*.tsv 2>/dev/null | sort >"$REPORT/report.tsv"

if [ "$UPDATE_BASELINE" = 1 ]; then
  newfile="$BASELINE.tmp"
  {
    printf '# scripts/cgate-baseline.txt — разобранные находки гейта (база/ратчет).\n'
    printf '# Формат: файл<TAB>ступень<TAB>класс<TAB>текст<TAB>вердикт\n'
    printf '# Инструменты: %s\n' "$(tool_versions | tr '\n' ' ')"
    [ -r "$BASELINE" ] && grep -E '^[^#]' "$BASELINE"
    cat "$REPORT"/all/*.tsv 2>/dev/null | cut -f1-4 |
      while IFS=$'\t' read -r f st cls msg; do
        printf '%s\t%s\t%s\t%s\tНОВОЕ — требует вердикта\n' "$f" "$st" "$cls" "$msg"
      done
  } | awk -F'\t' '!seen[$1 FS $2 FS $3 FS $4]++' >"$newfile"
  mv "$newfile" "$BASELINE"
  echo "база обновлена: $BASELINE ($(grep -cE '^[^#]' "$BASELINE") записей; вердикт «НОВОЕ» требует разбора)"
fi

if [ "$DIGEST" = 1 ]; then
  cat "$REPORT"/all/*.tsv 2>/dev/null | awk -F'\t' '{printf "%s\t%s\t%s\t%s\n", $1, $2, $3, $4}' | sort -u
  # советы S3 — тоже вход для триажа: помечаются, чтобы не путать с гейтом
  cat "$REPORT"/advisory/*.tsv 2>/dev/null | awk -F'\t' '{printf "%s\t%s\t%s\t%s\n", $1, $2, $3, $4}' | sort -u
  exit 0
fi

awk -F'\t' '
  { st=$2; seen[st]++
    if ($3=="NEW") new[st]+=$4
    if ($3=="BASELINED") base[st]+=$5
    if ($3=="OK") ok[st]++
    if ($3=="TOOLFAIL") fail[st]++
    if (st=="S-ПРОПУСК") skip[st]++ }
  END {
    printf "%-10s %8s %8s %8s %8s %8s\n", "ступень", "OK", "БАЗА", "НОВЫХ", "СБОЙ", "ПРОПУСК"
    split("S0 S1 S2 S3 S4", order, " ")
    for (i=1; i<=5; i++) { st=order[i]; if (seen[st]=="") continue
      printf "%-10s %8d %8d %8d %8d %8d\n", st, ok[st], base[st], new[st], fail[st], skip[st] }
    if (seen["S-ПРОПУСК"] != "") printf "%-10s %8d %8d %8d %8d %8d\n", "S-ПРОПУСК", 0,0,0,0, skip["S-ПРОПУСК"]
  }' "$REPORT/report.tsv"

# Советы S3: НЕ гейтят (политика в tidy_blocks), но обязаны быть видны —
# иначе «не блокирует» превратится в «не существует».
adv_n=$(cat "$REPORT"/advisory/*.tsv 2>/dev/null | wc -l)
if [ "$adv_n" -gt 0 ]; then
  echo
  echo "--- S3-СОВЕТЫ (не гейтят; плановый разбор): $adv_n ---"
  cat "$REPORT"/advisory/*.tsv | cut -f3 | sort | uniq -c | sort -rn | head -8 |
    sed 's/^/      /'
  echo "      полный список: $REPORT/advisory/"
fi

if awk -F'\t' '$3=="NEW"||$3=="TOOLFAIL"{f=1} END{exit !f}' "$REPORT/report.tsv"; then
  echo
  echo "--- новые находки (каждая строка = находка; разбирать с минимальным репро) ---"
  cat "$REPORT"/findings/*.tsv 2>/dev/null | sort -u | head -40
  n=$(cat "$REPORT"/findings/*.tsv 2>/dev/null | sort -u | wc -l)
  [ "$n" -gt 40 ] && echo "... всего $n (полный список: $REPORT/findings/)"
  echo
  echo ">>> cgate: НОВЫЕ НАХОДКИ ($n) — гейт красный"
  exit 1
fi

echo ">>> cgate: CLEAN (новых находок нет)"
exit 0
