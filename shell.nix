# shell.nix — ПИН инструментов гейта (07-10).
#
# Зачем: до этого каждый скрипт делал свой `nix-shell -p gcc clang-tools ...`,
# то есть версии движков плавали вместе с каналом nixpkgs. Находки анализатора
# версие-зависимы, а база известных срабатываний (scripts/cgate-baseline.txt)
# и кэш ступеней — бессмысленны, если движок под тем же именем меняется.
# Здесь инструменты объявлены один раз и берутся из ОДНОГО набора nixpkgs.
#
# Использование:
#   nix-shell                 # интерактивно
#   make gate                 # гейт сам перезапустится внутри этого шелла
#
# Проверка версий, с которыми снимались числа в build/cgate/report.tsv:
#   nix-shell --run 'make gate-versions'

{ pkgs ? import <nixpkgs> { } }:

pkgs.mkShell {
  packages = with pkgs; [
    gcc # S1, S4 (gcc -fanalyzer), сборка
    clang-tools # S0 (clang-format), S3 (clang-tidy), check-int
    clang # check-int: -fsanitize=integer,implicit-conversion
    cppcheck # S2
    cbmc # S6 (cverify.sh, формальный слой)
    valgrind # S7 (make valgrind)
    binutils # objdump для make check-simd
    lapack
    blas
    pkg-config
    llvmPackages.openmp # omp.h для clang-tidy (иначе tools/*.c не разбираются)
  ];

  shellHook = ''
    # Явный путь к omp.h нужен генератору compile_commands.json: clang-tidy
    # получает флаги из базы, а не из NIX_CFLAGS_COMPILE.
    export HZ_OMP_INC="${pkgs.llvmPackages.openmp.dev}/include"
    export HZ_SHELL_PIN="nixpkgs-${pkgs.gcc.version}"
  '';
}
