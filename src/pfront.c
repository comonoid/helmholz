/* pfront.c — обход фронта с переносом состояния. Разбор — в `pfront.h`. */
#include "pfront.h"
#include "pclip.h"
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* УСТРОЙСТВО РЕКУРСИИ (§241.7: один визит на узел, порядок октантов).
 *
 * У ячейки три ВХОДЯЩИЕ грани — по одной на ось, с той стороны, откуда идёт свет
 * (`lo`, если `d[c] > 0`, иначе `hi`), и три ИСХОДЯЩИЕ. При дроблении ребёнок
 * берёт входящую грань либо У РОДИТЕЛЯ (тогда это ЧЕТВЕРТЬ родительской грани, и
 * состояние получается вычислением полинома — Р3), либо У СОСЕДА-РЕБЁНКА,
 * посчитанного раньше. Порядок «раньше» задаётся числом осей, по которым ребёнок
 * стоит с ДАЛЬНЕЙ стороны: сперва ближний угол, потом три соседних, потом три,
 * потом дальний. Это и есть октантный порядок, и он же даёт зависимости без
 * циклов.
 *
 * Исходящие грани родителя собираются из детских: на каждую приходится четыре
 * детские грани, и они складываются СЛОЖЕНИЕМ МОМЕНТОВ (Р2, точная L²-проекция).
 */

/* Индексы осей грани: для оси `a` местные оси грани — две другие, по возрастанию. */
static void face_axes(int a, int *p, int *q) {
  *p = (a + 1) % 3;
  *q = (a + 2) % 3;
}

/* Р3: ЧЕТВЕРТЬ грани. `su`, `sv` — знаки смещения четверти (−1 или +1). Центр
 * четверти лежит в `±1/2` местных координат, полуразмер вдвое меньше — отсюда
 * наклоны делятся пополам, а центр смещается на половину наклона. */
static hz_pfront_face face_quarter(const hz_pfront_face *F, double su, double sv) {
  hz_pfront_face R;
  R.c0 = F->c0 + 0.5 * su * F->cu + 0.5 * sv * F->cv;
  R.cu = 0.5 * F->cu;
  R.cv = 0.5 * F->cv;
  return R;
}

/* Р2: сборка грани из четырёх четвертей — сложение моментов. Точная L²-проекция
 * линейной функции с четырёх подквадратов на один: среднее даёт `c0`, а наклон
 * складывается из наклонов четвертей и из разности их средних. */
static hz_pfront_face face_join(const hz_pfront_face *Q) {
  /* Q[0] = (−,−), Q[1] = (+,−), Q[2] = (−,+), Q[3] = (+,+) */
  hz_pfront_face R;
  R.c0 = 0.25 * (Q[0].c0 + Q[1].c0 + Q[2].c0 + Q[3].c0);
  R.cu = 0.25 * (Q[0].cu + Q[1].cu + Q[2].cu + Q[3].cu) +
         0.5 * (Q[1].c0 + Q[3].c0 - Q[0].c0 - Q[2].c0) * 0.5;
  R.cv = 0.25 * (Q[0].cv + Q[1].cv + Q[2].cv + Q[3].cv) +
         0.5 * (Q[2].c0 + Q[3].c0 - Q[0].c0 - Q[1].c0) * 0.5;
  return R;
}

/* ОПЕРАТОР ПУСТОГО УЗЛА — ТОЖДЕСТВО (§241.7). Доля открытого диска постоянна
 * ВДОЛЬ ЛУЧА, поэтому исходящая грань несёт то же, что входящая, с точностью до
 * пересчёта местных координат. Для осевых граней и общего направления пересчёт
 * этот аффинный, а значит линейность сохраняется ТОЧНО — то самое свойство,
 * ради которого состояние взято линейным (Р1). */
static hz_pfront_face face_shift(const hz_pfront_face *F) {
  return *F;
}

/* КОНСЕРВАТИВНАЯ ОЦЕНКА ПЕРЕКРЫТИЯ (А517). Складывать доли отдельных кусков
 * нельзя: их проекции пересекаются, и сумма объявит закрытым угол, который
 * открыт, — ошибка несимметрична и гасит свет там, где он есть (класс А475).
 * Поэтому берётся МАКСИМУМ по кускам, то есть заведомая НИЖНЯЯ оценка
 * перекрытия. Следствие: тени выходят слабее настоящих, огрубления МЕНЬШЕ, и
 * измеренный выигрыш есть НИЖНЯЯ граница. Ошибаться в эту сторону безопасно. */
static double cover_max(double m) {
  return m > 1.0 ? 1.0 : m;
}

/* Площадь проекции куска на грань оси `a`, делённая на площадь грани. Проекция
 * идёт ВДОЛЬ НАПРАВЛЕНИЯ, а не по нормали: заслоняет то, что стоит на пути. */
static double piece_cover(const hz_pclip_poly *P, const double *lo, const double *hi,
                          const double *dir, int a) {
  int p, q;
  face_axes(a, &p, &q);
  /* Постусловие `nv ≤ 9` доказано в `pclip`, но СЮДА кусок приходит извне, и
   * анализатор справедливо не берёт его на веру. Проверка явная — она же граница
   * буферов ниже. */
  int nv = P->nv;
  if (nv < 3 || nv > HZ_PCLIP_MAXV) return 0.0;
  double du = hi[p] - lo[p], dv = hi[q] - lo[q];
  if (!(du > 0.0) || !(dv > 0.0)) return 0.0;
  /* Снос вдоль направления на грань `a`: параметр берётся так, чтобы точка легла
   * в плоскость грани. При `dir[a] == 0` кусок вдоль грани не сносится вовсе. */
  if (dir[a] < 1e-300 && dir[a] > -1e-300) return 0.0;
  double fa = (dir[a] > 0.0) ? hi[a] : lo[a];
  /* Явная инициализация: анализатор не связывает заполнение в первом цикле с
   * чтением во втором (известный класс его слабостей, ср. разбор diam 07-24 в
   * CLAUDE.md). Девять записей на кусок против шести делений — цена никакая, а
   * значение определено по построению, а не по рассуждению. */
  double u[HZ_PCLIP_MAXV] = {0}, v[HZ_PCLIP_MAXV] = {0};
  for (int i = 0; i < nv; i++) {
    double t = (fa - P->v[i][a]) / dir[a];
    u[i] = (P->v[i][p] + t * dir[p] - lo[p]) / du;
    v[i] = (P->v[i][q] + t * dir[q] - lo[q]) / dv;
  }
  double s = 0.0;
  for (int i = 0; i < nv; i++) {
    int j = (i + 1) % nv;
    s += u[i] * v[j] - u[j] * v[i];
  }
  s = 0.5 * (s < 0.0 ? -s : s);
  return s > 1.0 ? 1.0 : s;
}

/* Ослабить грань перекрытием `c` (доля [0,1]), сохраняя линейность. */
static hz_pfront_face face_block(const hz_pfront_face *F, double c) {
  hz_pfront_face R;
  double k = 1.0 - c;
  R.c0 = F->c0 * k;
  R.cu = F->cu * k;
  R.cv = F->cv * k;
  return R;
}

int hz_pfront_step(const hz_pfront_face in[3], hz_pfront_face out[3], const hz_pclip_poly *pieces,
                   int npiece, const double *lo, const double *hi, const double *dir) {
  double cov[3];
  for (int a = 0; a < 3; a++) {
    /* Максимум копится на месте: временный массив здесь не нужен, а cppcheck на
     * нём справедливо ругался — при `npiece == 0` он оставался незаполненным. */
    double m = 0.0;
    for (int i = 0; i < npiece; i++) {
      double c = piece_cover(&pieces[i], lo, hi, dir, a);
      if (c > m) m = c;
    }
    cov[a] = cover_max(m);
  }
  for (int a = 0; a < 3; a++) {
    hz_pfront_face f = face_shift(&in[a]);
    out[a] = face_block(&f, cov[a]);
  }
  return 0;
}

const char *hz_pfront_note(void) {
  return "состояние линейное (Р1); перекрытие — МАКСИМУМ по кускам, то есть нижняя "
         "оценка (А517), значит и выигрыш от огрубления — нижняя граница";
}

/* ---- ОБХОД С ПЕРЕНОСОМ СОСТОЯНИЯ ------------------------------------------ */
/* Три причины остановить спуск, и третья — ради чего шаг затевался:
 *   ЛИСТ      дерево кончилось;
 *   ПОЛ       угловой размер ячейки упал ниже пола прохода;
 *   РОВНОЕ    все три входящие грани представимы линейно в пределах допуска —
 *             это и есть ОГРУБЛЕНИЕ ЗАСЛОНЁННОГО (и освещённого) вместо
 *             отсечения. Ограничено потолком `HZ_PFRONT_COARSEN_MAX`, и
 *             срабатывание потолка считается отдельно. */
static int box_hits3(const double *bl, const double *bh, const double *clo, const double *chi) {
  for (int c = 0; c < 3; c++)
    if (!(bl[c] < chi[c]) || !(bh[c] >= clo[c])) return 0;
  return 1;
}

static void pf_pieces(hz_pfront_ctx *X, const double *lo, const double *hi, const int32_t *list,
                      int32_t n, hz_pclip_poly *out, int *np) {
  int k = 0;
  for (int32_t i = 0; i < n && k < HZ_PFRONT_MAXPIECE; i++) {
    int32_t t = list[i];
    const double *A = X->m->v + 3 * (size_t)X->m->f[3 * (size_t)t + 0];
    const double *B = X->m->v + 3 * (size_t)X->m->f[3 * (size_t)t + 1];
    const double *C = X->m->v + 3 * (size_t)X->m->f[3 * (size_t)t + 2];
    hz_pclip_poly P;
    hz_pclip_tri(A, B, C, lo, hi, &P);
    if (P.nv >= 3 && hz_pclip_area(&P) > 0.0) out[k++] = P;
  }
  if (n > HZ_PFRONT_MAXPIECE) X->npclip++;
  *np = k;
}

void hz_pfront_walk(hz_pfront_ctx *X, int32_t nid, const double *lo, const double *hi,
                    const hz_pfront_face in[3], hz_pfront_face out[3], const int32_t *list,
                    int32_t n, int coarsened) {
  const hz_ptnode *N = &X->T->nd[nid];
  int stop = 0, why = 0;
  if (N->child < 0) {
    stop = 1;
    why = 0;
  } else {
    double r2 = 0.0, t = 0.0;
    for (int c = 0; c < 3; c++) {
      double h = 0.5 * (hi[c] - lo[c]);
      r2 += h * h;
      t += 0.5 * (lo[c] + hi[c]) * X->dir[c];
    }
    double dep = t - X->t_entry, d2 = dep * dep;
    if (d2 > r2 && 4.0 * r2 <= X->pxeps2 * d2) {
      stop = 1;
      why = 1;
    } else if (hz_pfront_flat(&in[0]) && hz_pfront_flat(&in[1]) && hz_pfront_flat(&in[2])) {
      if (coarsened < HZ_PFRONT_COARSEN_MAX) {
        stop = 1;
        why = 2;
      } else {
        X->ncap++;
      }
    }
  }
  if (stop) {
    X->ncell++;
    if (why == 0) X->stop_leaf++;
    if (why == 1) X->stop_floor++;
    if (why == 2) X->stop_flat++;
    hz_pclip_poly pc[HZ_PFRONT_MAXPIECE];
    int np = 0;
    if (n > 0) pf_pieces(X, lo, hi, list, n, pc, &np);
    if (np > 0)
      X->ncell_geo++;
    else
      X->ncell_void++;
    hz_pfront_step(in, out, pc, np, lo, hi, X->dir);
    double f = out[0].c0 + out[1].c0 + out[2].c0;
    if (f <= 3.0 * HZ_PFRONT_FRAC_TOL)
      X->nshadow++;
    else if (f >= 3.0 * (1.0 - HZ_PFRONT_FRAC_TOL))
      X->nlit++;
    return;
  }

  int32_t c0 = N->child;
  double mid[3];
  for (int c = 0; c < 3; c++)
    mid[c] = 0.5 * (lo[c] + hi[c]);
  int e[3];
  for (int c = 0; c < 3; c++)
    e[c] = (X->dir[c] > 0.0) ? 0 : 1;
  int32_t *sub = malloc((size_t)(n > 0 ? n : 1) * sizeof *sub);
  if (sub == NULL) {
    X->fail = 1;
    return;
  }
  hz_pfront_face cout[8][3];
  /* Порядок «ближние раньше дальних»: по числу осей, где ребёнок с дальней
   * стороны. Он делает зависимости ациклическими (§241.7). */
  for (int far = 0; far <= 3; far++) {
    for (int k = 0; k < 8; k++) {
      int b[3], nf = 0;
      for (int c = 0; c < 3; c++) {
        b[c] = (k >> c) & 1;
        if (b[c] != e[c]) nf++;
      }
      if (nf != far) continue;
      double clo[3], chi[3];
      for (int c = 0; c < 3; c++) {
        clo[c] = b[c] ? mid[c] : lo[c];
        chi[c] = b[c] ? hi[c] : mid[c];
      }
      hz_pfront_face cin[3];
      for (int a = 0; a < 3; a++) {
        int p, q;
        face_axes(a, &p, &q);
        if (b[a] == e[a]) {
          double su = b[p] ? 1.0 : -1.0, sv = b[q] ? 1.0 : -1.0;
          cin[a] = face_quarter(&in[a], su, sv);
        } else {
          cin[a] = cout[k ^ (1 << a)][a];
        }
      }
      int32_t ns = 0;
      for (int32_t i = 0; i < n; i++) {
        const double *bl = X->tlo + 3 * (size_t)list[i], *bh = X->thi + 3 * (size_t)list[i];
        if (box_hits3(bl, bh, clo, chi)) sub[ns++] = list[i];
      }
      hz_pfront_walk(X, c0 + k, clo, chi, cin, cout[k], sub, ns, coarsened);
      if (X->fail) {
        free(sub);
        return;
      }
    }
  }
  free(sub);
  /* Исходящие грани родителя — сложение моментов по четырём детским (Р2). */
  for (int a = 0; a < 3; a++) {
    int p, q;
    face_axes(a, &p, &q);
    hz_pfront_face Q[4];
    for (int j = 0; j < 4; j++) {
      int b[3];
      b[a] = 1 - e[a];
      b[p] = j & 1;
      b[q] = (j >> 1) & 1;
      Q[j] = cout[(b[0]) | (b[1] << 1) | (b[2] << 2)][a];
    }
    out[a] = face_join(Q);
  }
}
