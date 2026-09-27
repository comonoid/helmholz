/* poccref — ЭТАЛОН ЗАНЯТОСТИ (PLAN_ELEMENTS.md §385, шаг Ш0; аудит §386).
 *
 * ЗАЧЕМ. А665: приёмка Ш1 («потеряно полем») считает занятые ячейки ПО ДЕРЕВУ
 * ЗАНЯТОСТИ, а Ш1б это дерево сносит — метрика испарится ровно тогда, когда
 * понадобится. Здесь занятость получает определение, не зависящее НИ ОТ ОДНОЙ
 * структуры, которую мы правим: ни дерева занятости, ни дерева DC, ни ускорителей.
 *
 * ОПРЕДЕЛЕНИЕ ОДНО И НЕ НОВОЕ (`src/pocc.h`, О70): ячейка `[lo, lo+h)` занята
 * тогда и только тогда, когда найдётся треугольник, у которого кусок, отсечённый
 * коробкой ячейки, имеет ПОЛОЖИТЕЛЬНУЮ ПЛОЩАДЬ. Габаритный суррогат запрещён и
 * замерен негодным: 97.28 % занятых узлов против 12.79 % настоящих.
 *
 * ТРИ СПОСОБА СЧИТАТЬ ОДНО МНОЖЕСТВО, И У КАЖДОГО СВОЯ РОЛЬ.
 *   А РАССЕИВАНИЕ  по каждому треугольнику пройти ячейки его габарита. Габарит
 *                  здесь ТОЧНЫЙ ОТСЕВ, а не ускоритель: ячейка вне габарита не
 *                  пересекает треугольник ни при каком раскладе. Рабочая форма.
 *   Б ТУПОЙ        по каждой ячейке — ВСЕ треугольники. Буквальное определение,
 *                  `O(8^L·nt)`, только малые `L`. Существует затем, чтобы А не
 *                  приходилось брать на веру.
 *   В SAT          13 осей, ЗАМКНУТАЯ коробка. Другой алгоритм, но МАЖОРАНТА:
 *                  засчитывает касание нулевой площади. За точный не выдаётся
 *                  (класс Г40/Г44).
 *
 * ЧЕГО ЭТОТ ЭТАЛОН НЕ ПРОВЕРЯЕТ — А675, и умалчивать это нельзя. Способы А и Б
 * зовут ТОТ ЖЕ `hz_pclip_tri`, что и дерево занятости в `tools/pfield.c:239`.
 * Значит «А против Б = 0» доказывает РОВНО ОДНО: габаритный отсев не теряет
 * ячеек. О верности самого отсечения оно не говорит ничего. Независимость даёт
 * только способ В, и он мажоранта; негативный контроль самого отсечения — ключ
 * `naive` в `hz_pclip_tri_ex` (Г50), и он не здесь.
 *
 * ЭПСИЛОНОВ НЕТ НИ В ОДНОМ СПОСОБЕ (А682): `площадь > 0.0` — строгое неравенство
 * с нулём, а не порог; в SAT сравнения по осям тоже строгие. Допуск превратил бы
 * мажоранту в величину с неизвестной стороной ошибки.
 */
#include "occmap.h"

#include "pclip.h"
#include "scene_obj.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_s(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* --- рама сцены: формула ДОСЛОВНО из tools/pfield.c:923…940 ---------------- */

/* Другая формула дала бы другой набор ячеек на самой границе, и сверка мерила бы
 * разницу ФОРМУЛ, а не разницу путей. Поэтому здесь копия, а не «эквивалент». */
typedef struct {
  double org[3];
  double h;
  int32_t n;
  int lev;
} frame;

static void frame_of(frame *fr, const hz_objmesh *m, int lev) {
  fr->lev = lev;
  fr->n = (int32_t)1 << lev;
  double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
  for (int32_t k = 0; k < m->nv; k++)
    for (int c = 0; c < 3; c++) {
      double x = m->v[3 * (size_t)k + (size_t)c];
      if (x < lo[c]) lo[c] = x;
      if (x > hi[c]) hi[c] = x;
    }
  double side = 0.0;
  for (int c = 0; c < 3; c++)
    if (hi[c] - lo[c] > side) side = hi[c] - lo[c];
  fr->h = side / (double)(fr->n - 2);
  /* ПОЛКЛЕТКИ СМЕЩЕНИЯ — ДОСЛОВНО ТА ЖЕ ФОРМУЛА, ЧТО В `pfield.c` (§451). Без
   * него `h = сторона/(n−2)` ставит габарит РОВНО на границы ячеек, и осевая
   * стена попадает в плоскость сеточных вершин: её пересекают ОБА смежных ребра,
   * поверхность выдаётся дважды, площадь выходит в `1.5` раза больше истинной.
   * Формула обязана совпадать в обоих инструментах ДОСЛОВНО — иначе эталон Ш0
   * меряет другую сетку, чем рабочий путь, и сверка теряет смысл (§444). */
  for (int c = 0; c < 3; c++)
    fr->org[c] = 0.5 * (lo[c] + hi[c]) - 0.5 * (double)fr->n * fr->h + 0.5 * fr->h;
}

/* Коробка ячейки — тоже дословно (`cl + h`, а не `(x+1)·h`): сложение
 * неассоциативно, и другая запись даёт другой набор принятых треугольников. */
static void cell_box(const frame *fr, int64_t x, int64_t y, int64_t z, double cl[3], double ch[3]) {
  cl[0] = fr->org[0] + (double)x * fr->h;
  cl[1] = fr->org[1] + (double)y * fr->h;
  cl[2] = fr->org[2] + (double)z * fr->h;
  for (int c = 0; c < 3; c++)
    ch[c] = cl[c] + fr->h;
}

static void tri_verts(const hz_objmesh *m, int32_t t, const double **A, const double **B,
                      const double **C) {
  *A = m->v + 3 * (size_t)m->f[3 * (size_t)t + 0];
  *B = m->v + 3 * (size_t)m->f[3 * (size_t)t + 1];
  *C = m->v + 3 * (size_t)m->f[3 * (size_t)t + 2];
}

/* Диапазон ячеек габарита треугольника — формула из pfield.c:960…972. */
static void tri_cells(const hz_objmesh *m, const frame *fr, int32_t t, int64_t i0[3],
                      int64_t i1[3]) {
  for (int c = 0; c < 3; c++) {
    double a = 1e300, b = -1e300;
    for (int v = 0; v < 3; v++) {
      double x = m->v[3 * (size_t)m->f[3 * (size_t)t + (size_t)v] + (size_t)c];
      if (x < a) a = x;
      if (x > b) b = x;
    }
    /* РАСШИРЕНИЕ НА ОДНУ ЯЧЕЙКУ В КАЖДУЮ СТОРОНУ (§443). Прежний диапазон
     * назывался ТОЧНЫМ отсевом, и на зале это подтверждалось; на ОСЕВОЙ
     * геометрии он терял 52 % занятых ячеек: габарит треугольника вырожден по
     * оси, диапазон сжимается в один слой, а плоскость касается двух. Точность
     * обеспечивает точное отсечение, которое всё равно вызывается; габарит
     * обязан лишь НЕ ТЕРЯТЬ, и теперь он консервативен. */
    i0[c] = (int64_t)floor((a - fr->org[c]) / fr->h) - 1;
    i1[c] = (int64_t)floor((b - fr->org[c]) / fr->h) + 1;
    if (i0[c] < 0) i0[c] = 0;
    if (i1[c] >= fr->n) i1[c] = fr->n - 1;
  }
}

/* --- битовая карта занятости ----------------------------------------------- */

/* Раскладка, слепок, формат файла и разборщик — в `tools/occmap.h`, ОДНИ на оба
 * инструмента сверки (иначе разъедется формат, а выглядеть это будет как
 * расхождение путей). Здесь только сокращение под раму. */
static size_t cell_index(const frame *fr, int64_t x, int64_t y, int64_t z) {
  return hz_occ_index(fr->n, x, y, z);
}

/* --- СПОСОБ А: рассеивание -------------------------------------------------- */

static int64_t method_a(const hz_objmesh *m, const frame *fr, const unsigned char *drop,
                        unsigned char *bit, int64_t *nclip) {
  int64_t nocc = 0;
  for (int32_t t = 0; t < m->nt; t++) {
    if (drop != NULL && drop[t]) continue;
    int64_t i0[3], i1[3];
    tri_cells(m, fr, t, i0, i1);
    const double *A, *B, *C;
    tri_verts(m, t, &A, &B, &C);
    for (int64_t z = i0[2]; z <= i1[2]; z++)
      for (int64_t y = i0[1]; y <= i1[1]; y++)
        for (int64_t x = i0[0]; x <= i1[0]; x++) {
          double cl[3], ch[3];
          cell_box(fr, x, y, z, cl, ch);
          hz_pclip_poly P;
          (*nclip)++;
          if (hz_pclip_tri(A, B, C, cl, ch, &P) < 3) continue;
          if (!(hz_pclip_area(&P) > 0.0)) continue;
          size_t ci = cell_index(fr, x, y, z);
          if (!hz_occ_get(bit, ci)) {
            hz_occ_set(bit, ci);
            nocc++;
          }
        }
  }
  return nocc;
}

/* --- СПОСОБ Б: тупой полный перебор ---------------------------------------- */

/* Буквальное определение: для КАЖДОЙ ячейки спрашиваются ВСЕ треугольники.
 * Габаритного отсева нет вовсе — в этом весь смысл. */
static int64_t method_b(const hz_objmesh *m, const frame *fr, const unsigned char *drop,
                        unsigned char *bit, int64_t *nclip) {
  int64_t nocc = 0;
  for (int64_t z = 0; z < fr->n; z++)
    for (int64_t y = 0; y < fr->n; y++)
      for (int64_t x = 0; x < fr->n; x++) {
        double cl[3], ch[3];
        cell_box(fr, x, y, z, cl, ch);
        for (int32_t t = 0; t < m->nt; t++) {
          if (drop != NULL && drop[t]) continue;
          const double *A, *B, *C;
          tri_verts(m, t, &A, &B, &C);
          hz_pclip_poly P;
          (*nclip)++;
          if (hz_pclip_tri(A, B, C, cl, ch, &P) < 3) continue;
          if (!(hz_pclip_area(&P) > 0.0)) continue;
          hz_occ_set(bit, cell_index(fr, x, y, z));
          nocc++;
          break;
        }
      }
  return nocc;
}

/* --- СПОСОБ В: SAT, 13 осей, замкнутая коробка ------------------------------ */

/* НЕЗАВИСИМ ПО ПРЕДИКАТУ, а не по перечислению: ячейки берутся тем же габаритным
 * диапазоном, что и в А. Независимость заявляется только на алгоритм проверки —
 * ни `hz_pclip_tri`, ни площадь здесь не участвуют.
 * МАЖОРАНТА: замкнутая коробка засчитывает касание нулевой площади, а
 * полуоткрытое соглашение `[lo, hi)` отсечения — нет. Поэтому ожидается
 * `А ⊆ В`, и `|А \ В| = 0` есть проверка, а `|В \ А| > 0` — законно. */
static int sat_axis(const double p[3], const double a[3], const double e[3]) {
  double r = e[0] * fabs(a[0]) + e[1] * fabs(a[1]) + e[2] * fabs(a[2]);
  double mn = p[0], mx = p[0];
  for (int k = 1; k < 3; k++) {
    if (p[k] < mn) mn = p[k];
    if (p[k] > mx) mx = p[k];
  }
  return !(mn > r || mx < -r);
}

static int sat_tri_box(const double *A, const double *B, const double *C, const double ctr[3],
                       const double e[3]) {
  double v[3][3], f[3][3];
  for (int c = 0; c < 3; c++) {
    v[0][c] = A[c] - ctr[c];
    v[1][c] = B[c] - ctr[c];
    v[2][c] = C[c] - ctr[c];
  }
  for (int c = 0; c < 3; c++) {
    f[0][c] = v[1][c] - v[0][c];
    f[1][c] = v[2][c] - v[1][c];
    f[2][c] = v[0][c] - v[2][c];
  }
  /* девять осей e_i × f_j */
  for (int j = 0; j < 3; j++)
    for (int i = 0; i < 3; i++) {
      double a[3] = {0.0, 0.0, 0.0};
      int i1 = (i + 1) % 3, i2 = (i + 2) % 3;
      a[i1] = -f[j][i2];
      a[i2] = f[j][i1];
      double p[3];
      for (int k = 0; k < 3; k++)
        p[k] = a[0] * v[k][0] + a[1] * v[k][1] + a[2] * v[k][2];
      if (!sat_axis(p, a, e)) return 0;
    }
  /* три оси коробки */
  for (int c = 0; c < 3; c++) {
    double mn = v[0][c], mx = v[0][c];
    for (int k = 1; k < 3; k++) {
      if (v[k][c] < mn) mn = v[k][c];
      if (v[k][c] > mx) mx = v[k][c];
    }
    if (mn > e[c] || mx < -e[c]) return 0;
  }
  /* плоскость треугольника */
  double nrm[3];
  nrm[0] = f[0][1] * f[1][2] - f[0][2] * f[1][1];
  nrm[1] = f[0][2] * f[1][0] - f[0][0] * f[1][2];
  nrm[2] = f[0][0] * f[1][1] - f[0][1] * f[1][0];
  double d = nrm[0] * v[0][0] + nrm[1] * v[0][1] + nrm[2] * v[0][2];
  double r = e[0] * fabs(nrm[0]) + e[1] * fabs(nrm[1]) + e[2] * fabs(nrm[2]);
  return !(fabs(d) > r);
}

static int64_t method_c(const hz_objmesh *m, const frame *fr, const unsigned char *drop,
                        unsigned char *bit) {
  int64_t nocc = 0;
  double e[3] = {0.5 * fr->h, 0.5 * fr->h, 0.5 * fr->h};
  for (int32_t t = 0; t < m->nt; t++) {
    if (drop != NULL && drop[t]) continue;
    int64_t i0[3], i1[3];
    tri_cells(m, fr, t, i0, i1);
    const double *A, *B, *C;
    tri_verts(m, t, &A, &B, &C);
    for (int64_t z = i0[2]; z <= i1[2]; z++)
      for (int64_t y = i0[1]; y <= i1[1]; y++)
        for (int64_t x = i0[0]; x <= i1[0]; x++) {
          double cl[3], ch[3], ctr[3];
          cell_box(fr, x, y, z, cl, ch);
          for (int c = 0; c < 3; c++)
            ctr[c] = 0.5 * (cl[c] + ch[c]);
          if (!sat_tri_box(A, B, C, ctr, e)) continue;
          size_t ci = cell_index(fr, x, y, z);
          if (!hz_occ_get(bit, ci)) {
            hz_occ_set(bit, ci);
            nocc++;
          }
        }
  }
  return nocc;
}

/* РАЗБОР РАСХОЖДЕНИЯ А\В, а не отмахивание от него. `|А\В| = 0` предсказано
 * ДОКАЗАТЕЛЬСТВОМ (положительная площадь в полуоткрытой коробке влечёт
 * пересечение с замкнутой), поэтому всякий случай — либо дефект одного из двух
 * способов, либо ГРАНИЧНАЯ АРИФМЕТИКА. Различает их ПЛОЩАДЬ куска: доля от `h²`
 * порядка `1e-16` означает касание на последнем бите (коробка SAT строится как
 * `ctr ± h/2`, и `ctr = 0.5·(cl + ch)` округляется), а доля порядка единицы —
 * настоящий дефект. Печатается сама доля, а не вердикт. */
static void sat_diag(const hz_objmesh *m, const frame *fr, const unsigned char *drop,
                     const unsigned char *ba, const unsigned char *bc, int64_t nshow) {
  int64_t shown = 0;
  for (int64_t z = 0; z < fr->n && shown < nshow; z++)
    for (int64_t y = 0; y < fr->n && shown < nshow; y++)
      for (int64_t x = 0; x < fr->n && shown < nshow; x++) {
        size_t ci = cell_index(fr, x, y, z);
        if (!hz_occ_get(ba, ci) || hz_occ_get(bc, ci)) continue;
        double cl[3], ch[3];
        cell_box(fr, x, y, z, cl, ch);
        double amax = 0.0;
        for (int32_t t = 0; t < m->nt; t++) {
          if (drop != NULL && drop[t]) continue;
          int64_t i0[3], i1[3];
          tri_cells(m, fr, t, i0, i1);
          if (x < i0[0] || x > i1[0] || y < i0[1] || y > i1[1] || z < i0[2] || z > i1[2]) continue;
          const double *A, *B, *C;
          tri_verts(m, t, &A, &B, &C);
          hz_pclip_poly P;
          if (hz_pclip_tri(A, B, C, cl, ch, &P) < 3) continue;
          double a = hz_pclip_area(&P);
          if (a > amax) amax = a;
        }
        printf("      А\\В: ячейка [%lld %lld %lld], МАКС ПЛОЩАДЬ КУСКА %.3e м² = %.3e от h²\n",
               (long long)x, (long long)y, (long long)z, amax, amax / (fr->h * fr->h));
        shown++;
      }
}

/* --- чтение чужой карты ----------------------------------------------------- */

/* Формат, слепок и разборщик — в `tools/occmap.h`. Здесь только ввод-вывод: сам
 * разбор обязан остаться ЧИСТОЙ функцией, иначе его не возьмёт CBMC (А686). */
static unsigned char *occ_read(const char *path, int lev, hz_occ_hdr *h, size_t *nbuf) {
  FILE *f = fopen(path, "rb");
  if (f == NULL) return NULL;
  if (fseek(f, 0, SEEK_END) != 0) {
    fclose(f);
    return NULL;
  }
  long sz = ftell(f);
  if (sz <= 0 || fseek(f, 0, SEEK_SET) != 0) {
    fclose(f);
    return NULL;
  }
  unsigned char *buf = malloc((size_t)sz);
  if (buf == NULL) {
    fclose(f);
    return NULL;
  }
  size_t got = fread(buf, 1, (size_t)sz, f);
  fclose(f);
  if (got != (size_t)sz) {
    free(buf);
    return NULL;
  }
  int rc = hz_occ_parse(buf, got, lev, h);
  if (rc != 0) {
    fprintf(stderr, "дамп %s не годен: код %d\n", path, rc);
    free(buf);
    return NULL;
  }
  *nbuf = got;
  return buf;
}

/* --- трёхсторонняя сверка --------------------------------------------------- */

/* ЧИСЛОМ СВЕРЯТЬ НЕЛЬЗЯ (А6, А684): одинаковое число получается и при
 * разошедшихся множествах. Отсюда три величины, а не одна. */
typedef struct {
  int64_t both, only_l, only_r;
  int64_t only_r_out; /* «только справа» ВНЕ октанта [0, n/2)³ — для контроля drop0 */
} cmp3;

static void occ_cmp(const frame *fr, const unsigned char *L, const unsigned char *R, cmp3 *o) {
  memset(o, 0, sizeof *o);
  int32_t hn = fr->n / 2;
  for (int64_t z = 0; z < fr->n; z++)
    for (int64_t y = 0; y < fr->n; y++)
      for (int64_t x = 0; x < fr->n; x++) {
        size_t ci = cell_index(fr, x, y, z);
        int a = hz_occ_get(L, ci), b = hz_occ_get(R, ci);
        if (a && b)
          o->both++;
        else if (a)
          o->only_l++;
        else if (b) {
          o->only_r++;
          if (x >= hn || y >= hn || z >= hn) o->only_r_out++;
        }
      }
}

/* --- главная ---------------------------------------------------------------- */

int main(int argc, char **argv) {
  if (argc < 3) {
    fprintf(stderr, "poccref ФАЙЛ.obj МАСШТАБ [lev=N] [dumb] [sat] [dump=ПУТЬ] [cmp=ПУТЬ] "
                    "[drop0] [half]\n");
    return 2;
  }
  int lev = 6, dumb = 0, sat = 0, drop0 = 0, half = 0;
  const char *dump = NULL, *cmp = NULL;
  for (int i = 3; i < argc; i++) {
    if (strncmp(argv[i], "lev=", 4) == 0) lev = (int)strtol(argv[i] + 4, NULL, 10);
    /* СВЕРКА ЭКВИВАЛЕНТНОСТИ, а не негативный контроль (А10): способ Б есть
     * буквальное определение, и А обязан ему равняться. */
    if (strcmp(argv[i], "dumb") == 0) dumb = 1;
    /* НЕЗАВИСИМЫЙ ПРЕДИКАТ-МАЖОРАНТА (А675). */
    if (strcmp(argv[i], "sat") == 0) sat = 1;
    if (strncmp(argv[i], "dump=", 5) == 0) dump = argv[i] + 5;
    if (strncmp(argv[i], "cmp=", 4) == 0) cmp = argv[i] + 4;
    /* НЕГАТИВНЫЙ КОНТРОЛЬ (§385): отбросить треугольники, чей диапазон ячеек
     * целиком в октанте [0, n/2)³. Занятость обязана УПАСТЬ. */
    if (strcmp(argv[i], "drop0") == 0) drop0 = 1;
    /* ЗАПАСНОЙ КОНТРОЛЬ (А679), если drop0 окажется слаб на этих данных. */
    if (strcmp(argv[i], "half") == 0) half = 1;
  }
  if (lev < 1 || lev > 12) {
    fprintf(stderr, "lev вне 1..12 (выше карта не влезает в память по замыслу шага)\n");
    return 2;
  }
  hz_objmesh m;
  if (hz_obj_load(&m, argv[1], strtod(argv[2], NULL)) != 0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  frame fr;
  frame_of(&fr, &m, lev);
  size_t ncell = (size_t)fr.n * (size_t)fr.n * (size_t)fr.n;
  size_t nb = hz_occ_bytes(ncell);
  printf(
      "== ЭТАЛОН ЗАНЯТОСТИ: %s, треугольников %d, сетка %d^3, ячейка %.6f м, карта %.2f МБ%s%s\n",
      argv[1], m.nt, fr.n, fr.h, (double)nb / 1048576.0, drop0 ? "  [КОНТРОЛЬ drop0]" : "",
      half ? "  [КОНТРОЛЬ half]" : "");

  unsigned char *drop = NULL;
  if (drop0 || half) {
    drop = calloc((size_t)m.nt, 1);
    if (drop == NULL) exit(1);
    int64_t nd = 0;
    for (int32_t t = 0; t < m.nt; t++) {
      int d = 0;
      if (half && (t & 1)) d = 1;
      if (drop0) {
        int64_t i0[3], i1[3];
        tri_cells(&m, &fr, t, i0, i1);
        int in = 1;
        for (int c = 0; c < 3; c++)
          if (i1[c] >= fr.n / 2) in = 0;
        if (in) d = 1;
      }
      drop[t] = (unsigned char)d;
      nd += d;
    }
    printf("   ОТБРОШЕНО ТРЕУГОЛЬНИКОВ: %lld из %d (%.2f %%)\n", (long long)nd, m.nt,
           100.0 * (double)nd / (double)m.nt);
  }

  unsigned char *ba = calloc(nb, 1);
  if (ba == NULL) exit(1);
  int64_t nclip_a = 0;
  double t0 = now_s();
  int64_t na = method_a(&m, &fr, drop, ba, &nclip_a);
  double ta = now_s() - t0;
  printf("   [А рассеивание] ЗАНЯТЫХ ЯЧЕЕК %lld из %zu (%.3f %%) за %.2f с, отсечений %lld\n",
         (long long)na, ncell, 100.0 * (double)na / (double)ncell, ta, (long long)nclip_a);

  if (dumb) {
    unsigned char *bb = calloc(nb, 1);
    if (bb == NULL) exit(1);
    int64_t nclip_b = 0;
    t0 = now_s();
    int64_t nbo = method_b(&m, &fr, drop, bb, &nclip_b);
    double tb = now_s() - t0;
    cmp3 c;
    occ_cmp(&fr, ba, bb, &c);
    printf("   [Б тупой] ЗАНЯТЫХ %lld за %.2f с, отсечений %lld\n", (long long)nbo, tb,
           (long long)nclip_b);
    printf("   А ПРОТИВ Б (сверка ГАБАРИТНОГО ОТСЕВА, не отсечения — А675): общих %lld, только "
           "А %lld, только Б %lld\n",
           (long long)c.both, (long long)c.only_l, (long long)c.only_r);
    free(bb);
  }

  if (sat) {
    unsigned char *bc = calloc(nb, 1);
    if (bc == NULL) exit(1);
    t0 = now_s();
    int64_t nc = method_c(&m, &fr, drop, bc);
    double tc = now_s() - t0;
    cmp3 c;
    occ_cmp(&fr, ba, bc, &c);
    printf("   [В SAT, МАЖОРАНТА] ЗАНЯТЫХ %lld за %.2f с\n", (long long)nc, tc);
    printf("   А ПРОТИВ В: общих %lld, |А\\В| %lld (обязан быть 0), |В\\А| %lld (%.2f %% от |А|, "
           "законно: касание нулевой площади)\n",
           (long long)c.both, (long long)c.only_l, (long long)c.only_r,
           100.0 * (double)c.only_r / (double)(na ? na : 1));
    if (c.only_l > 0) sat_diag(&m, &fr, drop, ba, bc, 5);
    free(bc);
  }

  if (dump != NULL) {
    int rc = hz_occ_write(dump, fr.lev, fr.n, ba);
    printf("   ДАМП -> %s (код %d)\n", dump, rc);
  }

  if (cmp != NULL) {
    hz_occ_hdr h;
    size_t nbuf = 0;
    unsigned char *buf = occ_read(cmp, lev, &h, &nbuf);
    if (buf == NULL) {
      free(ba);
      free(drop);
      hz_obj_free(&m);
      return 1;
    }
    hz_occ_hdr mine;
    hz_occ_stamp(fr.lev, fr.n, ba, &mine);
    printf("   СЛЕПОК ЭТАЛОНА: занятых %lld, габарит [%d..%d][%d..%d][%d..%d]\n",
           (long long)mine.count, mine.blo[0], mine.bhi[0], mine.blo[1], mine.bhi[1], mine.blo[2],
           mine.bhi[2]);
    printf("   СЛЕПОК ДАМПА:   занятых %lld, габарит [%d..%d][%d..%d][%d..%d]\n",
           (long long)h.count, h.blo[0], h.bhi[0], h.blo[1], h.bhi[1], h.blo[2], h.bhi[2]);
    cmp3 c;
    occ_cmp(&fr, ba, buf + sizeof(hz_occ_hdr), &c);
    printf("   СВЕРКА С ДЕРЕВОМ: общих %lld, ТОЛЬКО ЭТАЛОН %lld, ТОЛЬКО ДЕРЕВО %lld (из них вне "
           "октанта [0,n/2)³ %lld)\n",
           (long long)c.both, (long long)c.only_l, (long long)c.only_r, (long long)c.only_r_out);
    free(buf);
  }

  free(ba);
  free(drop);
  hz_obj_free(&m);
  return 0;
}
