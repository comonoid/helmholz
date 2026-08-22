/* PLAN_TRANSPORT.md, шаги C и D. Трёхмерная развёртка: топологический порядок,
 * DG1 с противопотоком, итерация по рассеянию, ограничитель положительности. */

#include "transport/sweep3.h"
#include "transport/tet3.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Базис ячейки: {1, ξ, η, ζ}, ξ = (x_ед − c_ед)/s. Координаты ЕДИНИЧНЫЕ, потому
 * что коробка есть куб именно в них; в мир всё переносится множителями кадра. */
static void cell_basis(const tr3_mesh *m, int32_t c, const double xw[3], double b[4]) {
  b[0] = 1.0;
  double s = (double)m->csize[c];
  for (int k = 0; k < 3; k++) {
    double xu = (xw[k] - m->fr.o[k]) / m->fr.u[k];
    b[k + 1] = (xu - ((double)m->clo[c][k] + 0.5 * s)) / s;
  }
}

/* ∫ b^A_i b^B_j dA по многоугольнику в МИРЕ: базисы РАЗНЫХ ячеек (на стыке
 * уровней они разные), поэтому матрица не симметрична и считается смешанной. */
static void face_mass2(const tr3_mesh *m, const double (*v)[4][3], int32_t ca, int32_t cb,
                       double out[4][4]) {
  memset(out, 0, 16 * sizeof(double));
  const double (*vv)[3] = *v;
  for (int e = 1; e + 1 < 4; e++) {
    double tri[3][3];
    for (int k = 0; k < 3; k++) {
      tri[0][k] = vv[0][k];
      tri[1][k] = vv[e][k];
      tri[2][k] = vv[e + 1][k];
    }
    double area = tr3_poly_face_area(tri, 3);
    double ba[3][4], bb[3][4];
    for (int i = 0; i < 3; i++) {
      cell_basis(m, ca, tri[i], ba[i]);
      if (cb >= 0)
        cell_basis(m, cb, tri[i], bb[i]);
      else
        memcpy(bb[i], ba[i], sizeof ba[i]);
    }
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++) {
        double s1 = 0.0, sf = 0.0, sg = 0.0;
        for (int q = 0; q < 3; q++) {
          s1 += ba[q][i] * bb[q][j];
          sf += ba[q][i];
          sg += bb[q][j];
        }
        out[i][j] += area * (s1 + sf * sg) / 12.0;
      }
  }
}

/* Решение 4x4 с частичным выбором. Возврат 1 при вырождении. */
static int solve4(double a[4][4], double b[4], double x[4]) {
  for (int k = 0; k < 4; k++) {
    int best = k;
    for (int i = k + 1; i < 4; i++)
      if (fabs(a[i][k]) > fabs(a[best][k])) best = i;
    if (!(fabs(a[best][k]) > 0.0)) return 1;
    if (best != k) {
      for (int j = 0; j < 4; j++) {
        double t = a[k][j];
        a[k][j] = a[best][j];
        a[best][j] = t;
      }
      double t = b[k];
      b[k] = b[best];
      b[best] = t;
    }
    for (int i = k + 1; i < 4; i++) {
      double f = a[i][k] / a[k][k];
      for (int j = k; j < 4; j++)
        a[i][j] -= f * a[k][j];
      b[i] -= f * b[k];
    }
  }
  for (int i = 3; i >= 0; i--) {
    double s = b[i];
    for (int j = i + 1; j < 4; j++)
      s -= a[i][j] * x[j];
    x[i] = s / a[i][i];
  }
  return 0;
}

/* §731: ТРАССА ЦЕПОЧКИ — константы прибора.
 * HZ_CHAIN_MAX: одно направление пересекает не больше 3·2^lev плоскостей сетки
 * (на lev = 7 это 384); 512 — с запасом, а не подбор.
 * HZ_CHAIN_FLOOR: пол обрыва прогулки — уровень ЕДИНИЧНОГО входа `xunit`;
 * ниже него ячейка не усилена, идти дальше нечего (А1114: звено ниже пола
 * записывается последним, чтобы скачок через пол не потерялся).
 * HZ_CHAIN_AMP: порог «усиливающего шага» для сводки — здоровый разовый
 * перелёт DG1 замерен §729 на box: 1.267; 1.5 лежит выше него. */
#define HZ_CHAIN_MAX 512
static const double HZ_CHAIN_FLOOR = 1.0;
static const double HZ_CHAIN_AMP = 1.5;

/* §735: ПРИНЦИП МАКСИМУМА — допуск клипа. Сама граница уже ЗАВЫШЕНА (мажоранта
 * по углам коробки ⊇ флюида, хорда — диагональю), допуск покрывает только
 * накопление округления в цепочке сборки: 2^-38 ≈ 3.6e-12 при сотнях операций
 * double. Не подбор: у больных превышение в десятки раз, у точных решений —
 * равенство, и изменение допуска на порядки в обе стороны исход не меняет. */
static const double HZ_MAXP_TOL = 1.0 + 3.6e-12;

/* §735: точный максимум |полинома DG1| по восьми углам коробки. */
static double corner_amax(const double c[4]) {
  return fabs(c[0]) + 0.5 * (fabs(c[1]) + fabs(c[2]) + fabs(c[3]));
}

/* §690: ОБУСЛОВЛЕННОСТЬ МАТРИЦЫ МАСС ЯЧЕЙКИ, `‖M‖_F · ‖M⁻¹‖_F`.
 *
 * Обратная берётся ЧЕТЫРЬМЯ решениями `M x = e_i` тем же `solve4`, каким
 * пользуется рабочий путь: тогда мерится обусловленность ТОГО решателя, а не
 * абстрактная. `solve4` портит вход, поэтому копия обязательна.
 *
 * ЧИСЛО, ПО КОТОРОМУ ЭТО ПРОВЕРЯЕТСЯ, ИЗВЕСТНО ТОЧНО: у ЦЕЛОЙ кубической ячейки
 * матрица диагональна, `diag(V, V/12, V/12, V/12)`, и `κ = ‖M‖_F·‖M⁻¹‖_F`
 * от `V` не зависит — это `12·√(1+3/144)·√(1+3·144)/12`… считать не нужно:
 * важно, что величина ОДНА И ТА ЖЕ у всех целых ячеек любого размера. Разброс у
 * целых и есть проверка прибора. */
static double mass_cond(const double M[4][4]) {
  double nf = 0.0;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      nf += M[i][j] * M[i][j];
  nf = sqrt(nf);
  if (!(nf > 0.0)) return -1.0;
  double ni = 0.0;
  for (int k = 0; k < 4; k++) {
    double A[4][4], b[4] = {0, 0, 0, 0}, x[4];
    memcpy(A, M, sizeof A);
    b[k] = 1.0;
    if (solve4(A, b, x) != 0) return -2.0; /* вырождена: отдельный исход, не число */
    for (int i = 0; i < 4; i++)
      ni += x[i] * x[i];
  }
  return nf * sqrt(ni);
}

static int cmp_dbl_dbg(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return x < y ? -1 : (x > y ? 1 : 0);
}

static int cmp_i32_dbg(const void *a, const void *b) {
  int32_t x = *(const int32_t *)a, y = *(const int32_t *)b;
  return x < y ? -1 : (x > y ? 1 : 0);
}

/* §727: медиана значений > 0 по набору индексов. Значения ≤ 0 исключаются:
 * нуль у усиления проектора значит «не считано» (А1097, ложный ноль), а
 * отрицательная κ — флаг вырождения из mass_cond, не число. Возврат −1:
 * положительных значений в наборе нет. */
static double med_idx_dbg(const double *v, const int32_t *idx, int64_t n, double *scratch) {
  int64_t k = 0;
  for (int64_t i = 0; i < n; i++)
    if (v[idx[i]] > 0.0) scratch[k++] = v[idx[i]];
  if (k == 0) return -1.0;
  qsort(scratch, (size_t)k, sizeof *scratch, cmp_dbl_dbg);
  return scratch[k / 2];
}

/* §727: число различных значений; буфер сортируется на месте */
static int64_t nuniq_i32_dbg(int32_t *a, int64_t n) {
  if (n == 0) return 0;
  qsort(a, (size_t)n, sizeof *a, cmp_i32_dbg);
  int64_t u = 1;
  for (int64_t i = 1; i < n; i++)
    if (a[i] != a[i - 1]) u++;
  return u;
}

/* Решение 3x3 с частичным выбором. Возврат 1 при вырождении. */
static int solve3(double a[3][3], double b[3], double x[3]) {
  for (int k = 0; k < 3; k++) {
    int best = k;
    for (int i = k + 1; i < 3; i++)
      if (fabs(a[i][k]) > fabs(a[best][k])) best = i;
    if (!(fabs(a[best][k]) > 0.0)) return 1;
    if (best != k) {
      for (int j = 0; j < 3; j++) {
        double t = a[k][j];
        a[k][j] = a[best][j];
        a[best][j] = t;
      }
      double t = b[k];
      b[k] = b[best];
      b[best] = t;
    }
    for (int i = k + 1; i < 3; i++) {
      double f = a[i][k] / a[k][k];
      for (int j = k; j < 3; j++)
        a[i][j] -= f * a[k][j];
      b[i] -= f * b[k];
    }
  }
  for (int i = 2; i >= 0; i--) {
    double s = b[i];
    for (int j = i + 1; j < 3; j++)
      s -= a[i][j] * x[j];
    x[i] = s / a[i][i];
  }
  return 0;
}

/* ПРОЕКЦИЯ НА ЛИНЕЙНЫЕ ФУНКЦИИ, ЖИВУЩИЕ НА ПЛОСКОСТИ (К39).
 *
 * Задача: найти DG1-разложение поля, заданного своими моментами `rr_j = ∫ b_j E`,
 * по ПЛОСКОМУ элементу — грани сетки или куску поверхности. Наивно это
 * `solve4(M, rr)`, и ровно так было написано. НО МАТРИЦА `M` ЗДЕСЬ ВЫРОЖДЕНА
 * ТОЧНО: элемент плоский, а базис ячейки {1, ξ, η, ζ} на плоскости связан
 * тождеством `nul·b ≡ 0` (вывод — в `cut3.h`), поэтому ранг равен трём, а не
 * четырём. Измерено: `|M·nul| / (|M|·|nul|) = 1.4e-16` на ВСЕХ 1912 элементах
 * сцены, то есть вырождение не «почти», а по построению.
 *
 * Чем это было плохо. `solve4` вырождения НЕ ОБНАРУЖИВАЕТ — четвёртый ведущий
 * элемент выходит не нулём, а округлением, — и возвращает решение с ПРОИЗВОЛЬНОЙ
 * добавкой вдоль `nul`. На самой плоскости добавка не видна (там `nul·b = 0`),
 * зато ВНЕ плоскости она огромна. А проверка положительности `corner_min` считает
 * как раз по углам ЯЧЕЙКИ, то есть вне плоскости, и потому срабатывала почти
 * всегда: 1536 стенок из 1536 и 1372 поверхностных элементов из 1912. Хранимый
 * радианс падал до КОНСТАНТЫ (обещанный «DG1 по положению» не работал ни разу —
 * отсюда блоки на стенах), а исход переключателя решали последние биты, и он
 * МЕНЯЛСЯ от итерации к итерации: 1372 → 1330. Это и есть автоколебание, из-за
 * которого развёртка не сходилась 4000 итераций (К38).
 *
 * Как правильно. Решать в ТРЁХ неизвестных, а не в четырёх: тождество `nul·b = 0`
 * позволяет выразить одну степень свободы через остальные. Исключается та, у
 * которой коэффициент в `nul` НАИБОЛЬШИЙ ПО МОДУЛЮ, — тогда деление устойчиво, и
 * порога здесь не появляется, потому что выбор идёт по argmax, а не по сравнению
 * с числом. Индекс 0 (константа) не исключается никогда: без него не представить
 * постоянное поле.
 *
 * Возврат 1 — вырождение НАСТОЯЩЕЕ (элемент выродился в отрезок или точку), и
 * тогда вызывающий обязан отступить к константе. Это отказ «в закрытую сторону».
 */
int tr3_project_plane(const double M[4][4], const double rr[4], const double nul[4], double ee[4]) {
  int drop = 1;
  for (int k = 2; k < 4; k++)
    if (fabs(nul[k]) > fabs(nul[drop])) drop = k;
  if (!(fabs(nul[drop]) > 0.0)) return 1; /* нулевого вектора нет — не плоскость */
  int id[3], nid = 0;
  for (int k = 0; k < 4; k++)
    if (k != drop) id[nid++] = k;
  /* P: приведённые коэффициенты -> полные. ee[drop] = −Σ nul_k·g_k / nul_drop */
  double P[4][3];
  memset(P, 0, sizeof P);
  for (int k = 0; k < 3; k++) {
    P[id[k]][k] = 1.0;
    P[drop][k] = -nul[id[k]] / nul[drop];
  }
  double G[3][3], g[3], rg[3];
  for (int i = 0; i < 3; i++) {
    rg[i] = 0.0;
    for (int a = 0; a < 4; a++)
      rg[i] += P[a][i] * rr[a];
    for (int j = 0; j < 3; j++) {
      double s = 0.0;
      for (int a = 0; a < 4; a++)
        for (int b = 0; b < 4; b++)
          s += P[a][i] * M[a][b] * P[b][j];
      G[i][j] = s;
    }
  }
  if (solve3(G, rg, g) != 0) return 1;
  memset(ee, 0, 4 * sizeof(double));
  for (int a = 0; a < 4; a++)
    for (int k = 0; k < 3; k++)
      ee[a] += P[a][k] * g[k];
  return 0;
}

/* Нулевой вектор матрицы масс ГРАНИ СЕТКИ: та же формула, что в `cut3.h`, при
 * нормали вдоль оси и смещении `pos`. Грань тоже плоская, и её матрица масс
 * вырождена ровно так же — 1536 стенок из 1536 отступали к константе. */
void tr3_face_null(const tr3_mesh *m, int32_t f, int32_t c, double nul[4]) {
  double s = (double)m->csize[c];
  int a = (int)m->f[f].axis;
  memset(nul, 0, 4 * sizeof(double));
  nul[0] = ((double)m->clo[c][a] + 0.5 * s) - (double)m->f[f].pos;
  nul[a + 1] = s;
}

/* Минимум линейной функции по УГЛАМ ячейки: у DG1 экстремум всегда в углу. */
static double corner_min(const double c[4]) {
  double mn = 1e300;
  for (int k = 0; k < 8; k++) {
    double v = c[0];
    for (int a = 0; a < 3; a++)
      v += c[a + 1] * (((k >> a) & 1) ? 0.5 : -0.5);
    if (v < mn) mn = v;
  }
  return mn;
}

int tr3_sweep_solve(const tr3_problem *p, int maxit, double tol, double *phi, tr3_stats *st) {
  const tr3_mesh *m = p->m;
  const tr3_dirs *d = p->d;
  const tr3_cut *cu = p->cut;
  int32_t nc = m->ncell, nd = d->n, nse = cu != NULL ? cu->nse : 0;
  /* ВЛАДЕНИЕ ОБНУЛЯЕТСЯ ПЕРВОЙ СТРОКОЙ (К30): на путях раннего возврата поля
   * иначе остались бы неинициализированными, и вызывающий освободил бы мусор. */
  st->bout = NULL;
  st->eirr = NULL;
  st->sout = NULL;

  /* ПАРАЛЛЕЛЬНО ПО ОРДИНАТАМ. Внутри одной итерации направления НЕЗАВИСИМЫ:
   * каждое строит свой порядок обхода и своё угловое поле, а связывает их только
   * φ предыдущей итерации. Поэтому распараллеливание здесь не оптимизация с
   * риском, а прямое следствие устройства метода. Каждому потоку — своя копия
   * транзиентных массивов; Р1 при этом не нарушается: копий столько, сколько
   * ПОТОКОВ, а не сколько НАПРАВЛЕНИЙ. */
  /* ОДНОПОТОЧНО СОЗНАТЕЛЬНО: пока чинится алгоритм, параллелизм только мешает
   * читать. Структура под него уже готова — направления внутри итерации
   * независимы, у каждого свой транзиентный массив, — и включается он одной
   * строкой #pragma omp parallel for над циклом по mm. */
  int nth = 1;
  double *Lall = calloc((size_t)nth * (size_t)nc * 4, sizeof(double));
  double *phinall = calloc((size_t)nth * (size_t)nc * 4, sizeof(double));
  int32_t *indegall = calloc((size_t)nth * (size_t)nc, sizeof(int32_t));
  int32_t *orderall = calloc((size_t)nth * (size_t)nc, sizeof(int32_t));
  int32_t *queueall = calloc((size_t)nth * (size_t)nc, sizeof(int32_t));
  double *phin = calloc((size_t)nc * 4, sizeof(double));
  double *bout = calloc((size_t)m->nf * 4, sizeof(double));
  double *binf = calloc((size_t)m->nf * 4, sizeof(double));
  double *binfall = NULL, *sinfall = NULL;
  double *sout = calloc((size_t)(nse > 0 ? nse : 1) * 4, sizeof(double));
  double *sinf = calloc((size_t)(nse > 0 ? nse : 1) * 4, sizeof(double));
  /* §690: снимок `φ` прошлого такта — только под `trace`, рабочий путь не платит. */
  double *phiprev_dbg = NULL;
  /* §721: усиление проектора по элементам — только под `trace`. */
  double *gdbg = NULL;
  /* §723: обусловленность матрицы обновления ячейки, максимум по направлениям. */
  double *kadbg = NULL;
  /* §731: трасса цепочки — лучшая (по `|L|` в цели) цепочка такта. */
  int32_t *chnc = NULL;
  double *chnl = NULL;
  /* §735: макс мажоранты L по направлениям на элемент — граница облучённости */
  /* §746: то же для ГРАНЕЙ (стенки и стык) — граница обновления bout */
  double *bfmax = calloc((size_t)(0 < m->nf ? m->nf : 1), sizeof(double));
  double *semax = calloc((size_t)(nse > 0 ? nse : 1), sizeof(double));
  double *hs_se = calloc((size_t)(nse > 0 ? nse : 1), sizeof(double));
  double *hs_out = calloc((size_t)(nse > 0 ? nse : 1), sizeof(double));
  /* ЭТАП C: ЗЕРКАЛЬНЫЕ ГРАНИ. Индекс `mfid[f]` есть номер грани среди
   * зеркальных или −1; хранимое на них НАПРАВЛЕННОЕ, `nmf · nd · 4`. Память
   * платится только за зеркальные грани, а не за все. */
  int32_t *mfid = calloc((size_t)m->nf, sizeof(int32_t));
  int32_t nmf = 0;
  if (mfid != NULL) {
    for (int32_t f = 0; f < m->nf; f++) {
      mfid[f] = -1;
      if (m->f[f].cb >= 0) continue;
      int wl = (int)(~m->f[f].cb);
      if (p->wall_spec != NULL && p->wall_spec[wl] > 0.0) mfid[f] = nmf++;
    }
  }
  double *mspec = calloc((size_t)(nmf > 0 ? nmf : 1) * (size_t)nd * 4, sizeof(double));
  double *mspin = calloc((size_t)(nmf > 0 ? nmf : 1) * (size_t)nd * 4, sizeof(double));
  double *mprev = calloc((size_t)(nmf > 0 ? nmf : 1) * (size_t)nd * 4, sizeof(double));

  /* предыдущее ПОВЕРХНОСТНОЕ состояние — для невязки по всему состоянию (К84) */
  double *bprev = calloc((size_t)m->nf * 4, sizeof(double));
  double *sprev = calloc((size_t)(nse > 0 ? nse : 1) * 4, sizeof(double));
  binfall = calloc((size_t)nth * (size_t)m->nf * 4, sizeof(double));
  sinfall = calloc((size_t)nth * (size_t)(nse > 0 ? nse : 1) * 4, sizeof(double));
  if (Lall == NULL || phinall == NULL || indegall == NULL || orderall == NULL || queueall == NULL ||
      phin == NULL || binfall == NULL || sinfall == NULL || bout == NULL || binf == NULL ||
      sout == NULL || sinf == NULL || hs_se == NULL || hs_out == NULL || semax == NULL ||
      bfmax == NULL || bprev == NULL || sprev == NULL) {
    free(Lall);
    free(phinall);
    free(indegall);
    free(orderall);
    free(queueall);
    free(phin);
    free(binfall);
    free(sinfall);
    free(bout);
    free(binf);
    free(sout);
    free(sinf);
    free(kadbg);
    free(gdbg);
    free(phiprev_dbg);
    free(chnc);
    free(chnl);
    free(hs_se);
    free(semax);
    free(bfmax);
    free(hs_out);
    free(bprev);
    free(sprev);
    free(mfid);
    free(mspec);
    free(mspin);
    free(mprev);
    return 1;
  }

  /* К29 ДЛЯ ПРОИЗВОЛЬНОЙ НОРМАЛИ. Полусферный поток на дискретном наборе не
   * равен π (у наклонной поверхности тем более), и делить надо на СОБСТВЕННУЮ
   * сумму набора — иначе каждое отражение теряет энергию. */
  double hsum[6];
  for (int wl = 0; wl < 6; wl++) {
    int ax = wl / 2;
    double sgnw = (wl & 1) ? 1.0 : -1.0;
    hsum[wl] = 0.0;
    for (int mm = 0; mm < nd; mm++) {
      double on = (ax == 0 ? d->ox[mm] : (ax == 1 ? d->oy[mm] : d->oz[mm])) * sgnw;
      if (on > 0.0) hsum[wl] += d->w[mm] * on;
    }
  }
  for (int32_t e = 0; e < nse; e++) {
    double s = 0.0, so = 0.0;
    for (int mm = 0; mm < nd; mm++) {
      double on =
          d->ox[mm] * cu->se[e].n[0] + d->oy[mm] * cu->se[e].n[1] + d->oz[mm] * cu->se[e].n[2];
      if (on < 0.0)
        s += d->w[mm] * (-on); /* приходящие НА поверхность */
      else
        so += d->w[mm] * on; /* уходящие С поверхности — нужны балансу (К40) */
    }
    hs_se[e] = s;
    hs_out[e] = so;
  }

  /* §709 (Р1): СОГЛАСОВАНЫ ЛИ МАТРИЦЫ ГРАНИ. `ffmx[i][0]` и `ffm[i][0]` суть
   * один и тот же интеграл `∫b^ca_i dA` (базисный член `b_0 ≡ 1`), и то же для
   * `ffmx[0][j]` против `ffmb[j][0]`. Расхождение означало бы, что вылет из
   * одной ячейки и влёт в другую считаются по РАЗНЫМ величинам. Проверяется
   * один раз, до итераций. */
  if (p->trace > 0 && cu != NULL) {
    double wmax = 0.0, wmax_j = 0.0;
    int64_t nbadf = 0, nfl2 = 0;
    for (int32_t f = 0; f < m->nf; f++) {
      int32_t ca = m->f[f].ca, cb = m->f[f].cb;
      if (cb < 0 || cu->solid[ca] || cu->solid[cb]) continue;
      nfl2++;
      for (int i = 0; i < 4; i++) {
        double a = cu->ffmx[f][i][0], b = cu->ffm[f][i][0];
        double sc = fabs(a) + fabs(b);
        double e = sc > 0.0 ? fabs(a - b) / sc : 0.0;
        if (e > wmax) wmax = e;
        double a2 = cu->ffmx[f][0][i], b2 = cu->ffmb[f][i][0];
        double sc2 = fabs(a2) + fabs(b2);
        double e2 = sc2 > 0.0 ? fabs(a2 - b2) / sc2 : 0.0;
        if (e2 > wmax_j) wmax_j = e2;
        if (e > 1e-12 || e2 > 1e-12) {
          nbadf++;
          break;
        }
      }
    }
    printf("    §709 МАТРИЦЫ ГРАНИ: флюид-флюид граней %lld; max |ffmx[i][0] − ffm[i][0]| отн. "
           "%.3e; max |ffmx[0][j] − ffmb[j][0]| отн. %.3e; НАРУШЕНИЙ выше 1e-12: %lld\n",
           (long long)nfl2, wmax, wmax_j, (long long)nbadf);
  }
  /* §690 (ПРИЧИНА): распределение обусловленности матрицы масс ФЛЮИДНОЙ части,
   * отдельно у ЩЕПОК и у ЦЕЛЫХ. Печатается один раз, до итераций. */
  if (p->trace > 0 && cu != NULL) {
    double *kw = malloc((size_t)nc * sizeof *kw), *kt = malloc((size_t)nc * sizeof *kt);
    if (kw != NULL && kt != NULL) {
      int64_t nw = 0, nt = 0, nsing = 0;
      for (int32_t ci = 0; ci < nc; ci++) {
        double s3 = (double)m->csize[ci];
        double vcell = s3 * s3 * s3 * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
        double vfl = cu->mvol[ci][0][0];
        if (!(vfl > 0.0) || !(vcell > 0.0)) continue;
        double k = mass_cond(cu->mvol[ci]);
        if (k < 0.0) {
          nsing++;
          continue;
        }
        if (vfl < 0.01 * vcell)
          kt[nt++] = k;
        else
          kw[nw++] = k;
      }
      for (int64_t i = 1; i < nw; i++) {
        double x = kw[i];
        int64_t j = i - 1;
        while (j >= 0 && kw[j] > x) {
          kw[j + 1] = kw[j];
          j--;
        }
        kw[j + 1] = x;
      }
      for (int64_t i = 1; i < nt; i++) {
        double x = kt[i];
        int64_t j = i - 1;
        while (j >= 0 && kt[j] > x) {
          kt[j + 1] = kt[j];
          j--;
        }
        kt[j + 1] = x;
      }
      printf("    §690 ОБУСЛОВЛЕННОСТЬ МАТРИЦЫ МАСС: ЦЕЛЫХ %lld (медиана %.4g, p99 %.4g, макс "
             "%.4g); ЩЕПОК %lld (медиана %.4g, p99 %.4g, макс %.4g); ВЫРОЖДЕННЫХ %lld\n",
             (long long)nw, nw ? kw[nw / 2] : -1.0, nw ? kw[(nw * 99) / 100] : -1.0,
             nw ? kw[nw - 1] : -1.0, (long long)nt, nt ? kt[nt / 2] : -1.0,
             nt ? kt[(nt * 99) / 100] : -1.0, nt ? kt[nt - 1] : -1.0, (long long)nsing);
    }
    free(kw);
    free(kt);
  }

  /* §682: КОЭФФИЦИЕНТ ПЕРЕДАЧИ ЭЛЕМЕНТА. Две суммы выше считаются по ОДНОМУ
   * набору ординат, но по РАЗНЫМ полусферам, и точная квадратура дала бы обе
   * равными π. Исходящий радианс делится на `hs_se`, а мощность с элемента
   * набирается по `hs_out`, поэтому отношение `rho·hs_out/hs_se` есть
   * коэффициент передачи: больше единицы — оператор усиливает при ЛЮБОМ альбедо.
   * Печатается один раз, до итераций, и только при `trace` — это прибор, а не
   * рабочий путь. */
  if (p->trace > 0 && nse > 0) {
    double *v = malloc((size_t)nse * sizeof *v);
    if (v != NULL) {
      double gmax = 0.0, rmax = 0.0, axmax = 0.0;
      int64_t nax = 0, ngt1 = 0;
      for (int32_t e = 0; e < nse; e++) {
        v[e] = hs_se[e];
        double r = hs_se[e] > 0.0 ? hs_out[e] / hs_se[e] : 0.0;
        double rho = (p->facet_rho != NULL && cu->se[e].facet < p->nfacet)
                         ? p->facet_rho[cu->se[e].facet]
                         : 0.0;
        if (r > rmax) rmax = r;
        if (rho * r > gmax) gmax = rho * r;
        if (rho * r > 1.0) ngt1++;
        int axial = 0;
        for (int a = 0; a < 3; a++)
          if (fabs(cu->se[e].n[a]) > 0.996) axial = 1;
        if (axial) {
          nax++;
          if (fabs(r - 1.0) > axmax) axmax = fabs(r - 1.0);
        }
      }
      for (int32_t i = 1; i < nse; i++) { /* сортировка вставками мала: один раз */
        double x = v[i];
        int32_t j = i - 1;
        while (j >= 0 && v[j] > x) {
          v[j + 1] = v[j];
          j--;
        }
        v[j + 1] = x;
      }
      printf("    §682 ПЕРЕДАЧА ЭЛЕМЕНТА: элементов %d; hs_se мин %.4f p1 %.4f медиана %.4f "
             "макс %.4f (пи = %.4f); max hs_out/hs_se %.4f; MAX rho*hs_out/hs_se %.4f, у %lld "
             "элементов оно > 1; ОСЕВЫХ %lld, у них max |hs_out/hs_se - 1| = %.3e\n",
             nse, v[0], v[nse / 100], v[nse / 2], v[nse - 1], 3.14159265358979323846, rmax, gmax,
             (long long)ngt1, (long long)nax, axmax);
      free(v);
    }
  }
  /* ПРИ `warm_start` ВХОДНОЕ ПОЛЕ СОХРАНЯЕТСЯ — тогда один проход есть
   * применение ОПЕРАТОРА к заданному вектору, а не итерация от нуля (см.
   * `sweep3.h`). Ограничение про `bout`/`sout` там же. */
  if (p->trace > 0 && nse > 0) gdbg = calloc((size_t)nse, sizeof *gdbg);
  if (p->trace > 0) kadbg = calloc((size_t)nc, sizeof *kadbg);
  /* §731: цель трассы хранится в задаче со сдвигом +1 (см. sweep3.h) */
  const int32_t chain_tgt = p->chain_cell1 - 1;
  if (p->trace > 0 && p->chain_cell1 > 0 && p->chain_cell1 <= nc) {
    chnc = calloc(HZ_CHAIN_MAX, sizeof *chnc);
    chnl = calloc(HZ_CHAIN_MAX, sizeof *chnl);
  }
  if (p->trace > 0) {
    phiprev_dbg = calloc((size_t)nc, sizeof *phiprev_dbg);
    if (phiprev_dbg != NULL)
      for (int32_t ci = 0; ci < nc; ci++)
        phiprev_dbg[ci] = phi[4 * (size_t)ci];
  }
  if (!p->warm_start) memset(phi, 0, (size_t)nc * 4 * sizeof(double));
  /* К76: ПОВЕРХНОСТНОЕ СОСТОЯНИЕ ТОЖЕ МОЖЕТ ПРИЙТИ ИЗВНЕ. Без него «тёплый
   * старт» на сцене с отражением тёплым не является: `bout`/`sout` есть вторая
   * половина состояния итерации, и заводить их из одного излучения значит
   * начинать отражённую часть с нуля. */
  if (p->bout_in != NULL) {
    memcpy(bout, p->bout_in, (size_t)m->nf * 4 * sizeof(double));
    /* §750: при тёплом поверхностном старте ПРОШЛОЕ состояние — прокинутое, а
     * не ноль: демпфер §748 в первом такте смешивает с bprev, и нулевой bprev
     * гасил бы (1−ω) состояния к нулю КАЖДЫЙ вызов — прокидка maxit=1 вставала
     * на ложный пол (замерено: 9.34e-2 у γ=0). То же для sprev ниже. */
    memcpy(bprev, p->bout_in, (size_t)m->nf * 4 * sizeof(double));
  } else if (p->wall_rho != NULL)
    for (int32_t f = 0; f < m->nf; f++)
      if (m->f[f].cb < 0)
        bout[f * 4] = p->wall_emit != NULL ? p->wall_emit[(int)(~m->f[f].cb)] : 0.0;
  if (p->sout_in != NULL && nse > 0) {
    memcpy(sout, p->sout_in, (size_t)nse * 4 * sizeof(double));
    memcpy(sprev, p->sout_in, (size_t)nse * 4 * sizeof(double)); /* §750 */
  } else
    for (int32_t e = 0; e < nse; e++) {
      /* §670: излучение по ЭЛЕМЕНТУ имеет приоритет; фасетное осталось для
       * воспроизведения прежних прогонов. Оба места (здесь и в рабочем цикле)
       * правятся ВМЕСТЕ — иначе первая итерация разойдётся со следующими. */
      if (p->elem_emit != NULL)
        sout[e * 4] = p->elem_emit[e];
      else if (p->facet_emit != NULL && cu->se[e].facet < p->nfacet)
        sout[e * 4] = p->facet_emit[cu->se[e].facet];
    }

  int nclip_last = 0, it = 0, nfb = 0;
  /* §735: счётчики клипа принципа максимума — ячейки и элементы */
  int nmaxp_last = 0, nmaxpe_last = 0, nmaxpf_last = 0;
  /* §677: ОТКАТ РАСЩЕПЛЁН НА ЧЕТЫРЕ. К86 говорит, что условий два и природа у
   * них разная: вырождение элемента — свойство ГЕОМЕТРИИ (множество постоянно),
   * `corner_min < 0` — свойство ПОЛЯ, то есть ПЕРЕКЛЮЧАТЕЛЬ. Пока они считались
   * одним числом, сказать, ЧТО дрожит, было нечем. `_thin` — доля на ЩЕПКАХ. */
  int nfb_fg = 0, nfb_fp = 0, nfb_eg = 0, nfb_ep = 0, nfb_thin = 0;
  /* §731: лучшая цепочка такта — длина, направление, `|L|` в цели. */
  int chn_n = 0, chn_mm = -1;
  double chn_bestlt = -1.0;
  double resid = 0.0;
  /* К81: пол невязки ловится ЗАСТОЕМ, а не порогом на её величину */
  double best = 1e300;
  int nstall = 0;
  st->stalled = 0;
  for (it = 0; it < maxit; it++) {
    memset(phin, 0, (size_t)nc * 4 * sizeof(double));
    memset(binf, 0, (size_t)m->nf * 4 * sizeof(double));
    memset(sinf, 0, (size_t)(nse > 0 ? nse : 1) * 4 * sizeof(double));
    if (nmf > 0) memset(mspin, 0, (size_t)nmf * (size_t)nd * 4 * sizeof(double));
    memset(phinall, 0, (size_t)nth * (size_t)nc * 4 * sizeof(double));
    memset(binfall, 0, (size_t)nth * (size_t)m->nf * 4 * sizeof(double));
    memset(sinfall, 0, (size_t)nth * (size_t)(nse > 0 ? nse : 1) * 4 * sizeof(double));
    memset(semax, 0, (size_t)(nse > 0 ? nse : 1) * sizeof(double));
    memset(bfmax, 0, (size_t)m->nf * sizeof(double));
    nclip_last = 0;
    nmaxp_last = 0;
    nmaxpe_last = 0;
    nmaxpf_last = 0;
    nfb = 0;
    nfb_fg = nfb_fp = nfb_eg = nfb_ep = nfb_thin = 0;
    st->pin = st->pout = st->pabs = 0.0;
    double pin_acc = 0.0, pout_acc = 0.0;
    /* §705: поток, ушедший в СПЛОШНЫЕ ячейки и прежде нигде не учтённый. */
    double psolid_acc = 0.0;
    int fail = 0;

    for (int mm = 0; mm < nd; mm++) {
      int tid = 0;
      double *L = Lall + (size_t)tid * (size_t)nc * 4;
      double *phit = phinall + (size_t)tid * (size_t)nc * 4;
      int32_t *indeg = indegall + (size_t)tid * (size_t)nc;
      int32_t *order = orderall + (size_t)tid * (size_t)nc;
      int32_t *queue = queueall + (size_t)tid * (size_t)nc;
      double *binft = binfall + (size_t)tid * (size_t)m->nf * 4;
      double *sinft = sinfall + (size_t)tid * (size_t)(nse > 0 ? nse : 1) * 4;
      double om[3] = {d->ox[mm], d->oy[mm], d->oz[mm]};
      /* --- топологический порядок для этого направления ---
       * ПЕРЕСТРАИВАЕТСЯ КАЖДУЮ ИТЕРАЦИЮ, И ЭТО СОЗНАТЕЛЬНО: кэш стоил бы
       * ND×ncell целых, а это ровно то, что запрещает Р1. Цена O(ncell+nf) —
       * тот же порядок, что у самого прохода. */
      memset(indeg, 0, (size_t)nc * sizeof(int32_t));
      for (int32_t f = 0; f < m->nf; f++) {
        if (m->f[f].cb < 0) continue;
        double on = om[m->f[f].axis];
        if (on > 0.0)
          indeg[m->f[f].cb]++;
        else if (on < 0.0)
          indeg[m->f[f].ca]++;
      }
      int32_t qh = 0, qt = 0, no = 0;
      for (int32_t c = 0; c < nc; c++)
        if (indeg[c] == 0) queue[qt++] = c;
      while (qh < qt) {
        int32_t c = queue[qh++];
        order[no++] = c;
        for (int32_t k = m->fstart[c]; k < m->fstart[c + 1]; k++) {
          int32_t f = m->flist[k];
          if (m->f[f].cb < 0) continue;
          double on = om[m->f[f].axis];
          int32_t down = -1;
          if (on > 0.0 && m->f[f].ca == c)
            down = m->f[f].cb;
          else if (on < 0.0 && m->f[f].cb == c)
            down = m->f[f].ca;
          if (down >= 0 && --indeg[down] == 0) queue[qt++] = down;
        }
      }
      if (no != nc) {
        fail = 1; /* цикл обхода: наружу выходим флагом, а не return из omp */
        continue;
      }

      for (int32_t oi = 0; oi < no; oi++) {
        int32_t c = order[oi];
        if (cu != NULL && cu->solid[c]) { /* ячейка целиком в материале */
          /* §705/§707: НАКОПЛЕНИЕ ПЕРЕЕХАЛО В ЦИКЛ ПО ГРАНЯМ. Здесь оно стояло,
           * пока стык не имел граничного условия; теперь тот же поток считается
           * там же, где считается вылет через грани куба, и оставлять его тут
           * значит считать дважды — регрессия §707 П1 это и поймала (`5.6072 ->
           * 11.214`, ровно вдвое). */
          for (int j = 0; j < 4; j++)
            L[c * 4 + j] = 0.0;
          continue;
        }
        double A[4][4], rhs[4];
        memset(A, 0, sizeof A);
        memset(rhs, 0, sizeof rhs);
        double s = (double)m->csize[c];
        /* БЕЗ КОПИРОВАНИЯ: матрица берётся ССЫЛКОЙ. memcpy по 128 байт на ячейку
         * и на каждую грань — это килобайт лишней памяти на одно обновление
         * ячейки, а обновлений здесь сотни миллионов. */
        double MMbox[4][4];
        const double (*MM)[4];
        if (cu != NULL) {
          MM = cu->mvol[c];
        } else {
          memset(MMbox, 0, sizeof MMbox);
          double vol = s * s * s * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
          MMbox[0][0] = vol;
          for (int k = 1; k < 4; k++)
            MMbox[k][k] = vol / 12.0;
          MM = MMbox;
        }
        for (int j = 0; j < 4; j++)
          for (int i = 0; i < 4; i++)
            A[j][i] += p->sig_t[c] * MM[j][i];
        /* −∫ b_i (ω·∇b_j): градиент постоянен, поэтому член есть g_j·∫b_i dV, а
         * ∫b_i dV на ФЛЮИДНОЙ области уже не ноль (симметрии нет) — берётся
         * первая строка матрицы масс. */
        for (int j = 1; j < 4; j++) {
          double g = om[j - 1] / (s * m->fr.u[j - 1]);
          for (int i = 0; i < 4; i++)
            A[j][i] -= g * MM[0][i];
        }
        for (int j = 0; j < 4; j++) {
          double q = 0.0;
          for (int i = 0; i < 4; i++)
            q += p->sig_s[c] / (4.0 * M_PI) * phi[c * 4 + i] * MM[j][i];
          if (p->eps != NULL)
            for (int i = 0; i < 4; i++)
              q += p->eps[c * 4 + i] * MM[j][i];
          if (j == 0)
            q += (p->eps_dir[0] * om[0] + p->eps_dir[1] * om[1] + p->eps_dir[2] * om[2]) * MM[0][0];
          rhs[j] += q;
        }
        /* §735: граница принципа максимума собирается вместе со сборкой:
         * bin — максимум угловых мажорант ВЛЁТА, qc — полином источника. */
        double bin = 0.0;
        double qc[4];
        for (int i = 0; i < 4; i++)
          qc[i] = p->sig_s[c] / (4.0 * M_PI) * phi[c * 4 + i] +
                  (p->eps != NULL ? p->eps[c * 4 + i] : 0.0);
        qc[0] += p->eps_dir[0] * om[0] + p->eps_dir[1] * om[1] + p->eps_dir[2] * om[2];
        double qmax = corner_amax(qc);

        for (int32_t k = m->fstart[c]; k < m->fstart[c + 1]; k++) {
          int32_t f = m->flist[k];
          int mine_is_a = m->f[f].ca == c;
          double sgn;
          if (m->f[f].cb < 0) {
            int wall = (int)(~m->f[f].cb);
            sgn = (wall & 1) ? 1.0 : -1.0;
          } else {
            sgn = mine_is_a ? 1.0 : -1.0;
          }
          double on = om[m->f[f].axis] * sgn;
          if (!(fabs(on) > 0.0)) continue;
          const double (*fmine)[4] = cu != NULL ? (mine_is_a ? cu->ffm[f] : cu->ffmb[f]) : NULL;
          double fbox[4][4];
          if (cu == NULL) {
            double v[4][3];
            tr3_face_corners(m, f, v);
            face_mass2(m, (const double (*)[4][3]) & v, c, c, fbox);
            fmine = fbox;
          }
          if (on > 0.0) {
            for (int i = 0; i < 4; i++)
              for (int j = 0; j < 4; j++)
                A[j][i] += on * fmine[i][j];
          } else {
            if (m->f[f].cb < 0) {
              double lb[4];
              /* ЭТАП C: у ЗЕРКАЛЬНОЙ грани влёт зависит от ОРДИНАТЫ — хранимое
               * там направленное, и берётся оно по индексу `mm` напрямую.
               * Перестановка уже применена при записи, поэтому здесь никакого
               * поиска нет: чтение по тому же `mm`, что и всё остальное. */
              if (mfid != NULL && mfid[f] >= 0) {
                for (int i = 0; i < 4; i++)
                  lb[i] = mspec[((size_t)mfid[f] * (size_t)nd + (size_t)mm) * 4 + (size_t)i];
              } else if (p->wall_rho != NULL) {
                for (int i = 0; i < 4; i++)
                  lb[i] = bout[f * 4 + i];
              } else {
                /* предписанный влёт, линейный по положению: раскладываем по
                 * базису ЯЧЕЙКИ через значения в углах грани */
                /* Разложение ЛИНЕЙНОГО влёта по базису ЯЧЕЙКИ. Центр берётся
                 * ЯЧЕЙКИ, а не грани: базис b_i центрирован на ячейке, и подмена
                 * центра сдвигает постоянную часть на binc·(центр грани − центр
                 * ячейки) — то есть ровно на полклетки, что К12 и поймал. */
                double cen[3];
                for (int a = 0; a < 3; a++)
                  cen[a] = m->fr.o[a] + m->fr.u[a] * ((double)m->clo[c][a] + 0.5 * s);
                lb[0] = p->binc0;
                for (int a = 0; a < 3; a++)
                  lb[0] += p->binc[a] * (cen[a] - p->binx0[a]);
                for (int a = 0; a < 3; a++)
                  lb[a + 1] = p->binc[a] * s * m->fr.u[a];
              }
              double accj[4] = {0, 0, 0, 0};
              for (int j = 0; j < 4; j++)
                for (int i = 0; i < 4; i++)
                  accj[j] += lb[i] * fmine[i][j];
              for (int j = 0; j < 4; j++)
                rhs[j] -= on * accj[j];
              pin_acc += -on * d->w[mm] * accj[0];
              double bm = corner_amax(lb); /* §735 */
              if (bm > bin) bin = bm;
            } else {
              int32_t up = mine_is_a ? m->f[f].cb : m->f[f].ca;
              if (cu != NULL && cu->solid[up]) {
                /* §707: ВЛЁТ СО СТЫКА «ФЛЮИД — СПЛОШНОЕ». Радианс сплошной
                 * ячейки нулевой по построению (материал непрозрачен), поэтому
                 * брать надо не его, а ИСХОДЯЩИЙ радианс стены, хранимый на
                 * грани, — ровно как у граней куба. В `pin_acc` это НЕ идёт:
                 * втекло — про поток ИЗВНЕ области, а здесь возвращает своя же
                 * стена. */
                const double (*fmine2)[4] = mine_is_a ? cu->ffm[f] : cu->ffmb[f];
                double accs[4] = {0, 0, 0, 0};
                for (int j = 0; j < 4; j++)
                  for (int i = 0; i < 4; i++)
                    accs[j] += bout[f * 4 + i] * fmine2[i][j];
                for (int j = 0; j < 4; j++)
                  rhs[j] -= on * accs[j];
                double bm = corner_amax(&bout[f * 4]); /* §735 */
                if (bm > bin) bin = bm;
                continue;
              }
              double fxb[4][4];
              const double (*fx)[4];
              if (cu != NULL) {
                fx = cu->ffmx[f];
              } else {
                double v[4][3];
                tr3_face_corners(m, f, v);
                face_mass2(m, (const double (*)[4][3]) & v, m->f[f].ca, m->f[f].cb, fxb);
                fx = fxb;
              }
              for (int j = 0; j < 4; j++) {
                double acc = 0.0;
                for (int i = 0; i < 4; i++)
                  acc += L[up * 4 + i] * (mine_is_a ? fx[j][i] : fx[i][j]);
                rhs[j] -= on * acc;
              }
              double bm = corner_amax(&L[up * 4]); /* §735 */
              if (bm > bin) bin = bm;
            }
          }
        }

        /* ПОВЕРХНОСТНЫЕ ЭЛЕМЕНТЫ: та же роль, что у граничной грани, но нормаль
         * произвольная. ω·n > 0 — поверхность светит В ячейку (вток), < 0 — луч
         * упирается в неё (выток, и он же облучённость поверхности). */
        if (cu != NULL)
          for (int32_t k = cu->sestart[c]; k < cu->sestart[c + 1]; k++) {
            int32_t e = cu->selist[k];
            const tr3_selem *se = &cu->se[e];
            double on = om[0] * se->n[0] + om[1] * se->n[1] + om[2] * se->n[2];
            if (!(fabs(on) > 0.0)) continue;
            if (on < 0.0) { /* выток на поверхность */
              for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                  A[j][i] += (-on) * se->m[i][j];
            } else { /* вток с поверхности */
              for (int j = 0; j < 4; j++) {
                double acc = 0.0;
                for (int i = 0; i < 4; i++)
                  acc += sout[e * 4 + i] * se->m[i][j];
                rhs[j] += on * acc;
              }
              double bm = corner_amax(&sout[e * 4]); /* §735 */
              if (bm > bin) bin = bm;
            }
          }

        /* §733: снимок системы до solve4 — он портит вход. */
        int isdump = 0;
        if (p->trace > 0 && cu != NULL && p->dump_dir1 == mm + 1)
          for (int q = 0; q < 4; q++)
            if (p->dump_cell1[q] == c + 1) isdump = 1;
        double Adump[4][4], rdump[4];
        if (isdump) {
          memcpy(Adump, A, sizeof Adump);
          memcpy(rdump, rhs, sizeof rdump);
        }
        double a0row[4], rhs0 = rhs[0];
        for (int j = 0; j < 4; j++)
          a0row[j] = A[0][j];
        double cf[4];
        if (p->trace > 0 && kadbg != NULL) {
          /* §723: обусловленность матрицы ОБНОВЛЕНИЯ ячейки. Это не `mvol`:
           * `A` собирается на лету из объёма, граней и поверхностных элементов,
           * и именно её решает `solve4`. Берётся МАКСИМУМ по направлениям. */
          double ka = mass_cond(A);
          if (ka > kadbg[c]) kadbg[c] = ka;
        }
        if (solve4(A, rhs, cf) != 0)
          for (int j = 0; j < 4; j++)
            cf[j] = 0.0;
        if (p->limiter && corner_min(cf) < 0.0 && fabs(a0row[0]) > 0.0) {
          double kk = 0.0;
          for (int j = 1; j < 4; j++)
            kk += a0row[j] * cf[j];
          double dev = corner_min((double[4]){0.0, cf[1], cf[2], cf[3]});
          double c0h = rhs0 / a0row[0], kh = kk / a0row[0];
          double alpha = 1.0, den = kh - dev;
          if (den > 0.0) {
            alpha = c0h / den;
            if (alpha > 1.0) alpha = 1.0;
            if (alpha < 0.0) alpha = 0.0;
          }
          if (alpha < 1.0) {
            for (int j = 1; j < 4; j++)
              cf[j] *= alpha;
            cf[0] = (rhs0 - alpha * kk) / a0row[0];
            nclip_last++;
          }
        }
        /* §735: ПРИНЦИП МАКСИМУМА. Вдоль луча L′ + σ_t·L = q, поэтому max L в
         * ячейке не превышает максимума ВЛЁТА и q/σ_t (при σ_t = 0 вдоль
         * хорды прибавляется не больше q_max·хорда). Решение за границей
         * заменяется КОНСТАНТНЫМ БАЛАНСНЫМ rhs0/A00: строка 0 (баланс потока)
         * выполняется ТОЧНО, жертвуются наклоны — консервативность выше
         * поточечной границы (§734, К40). Механизм и адреса — §734: полином
         * ±1e+09 при влёте 1.9e+06 (расщепление), среднее 365 при потоке 0.73
         * (щепка). */
        if (!p->maxp_off) {
          double s3d = (double)m->csize[c];
          double chord = s3d * sqrt(m->fr.u[0] * m->fr.u[0] + m->fr.u[1] * m->fr.u[1] +
                                    m->fr.u[2] * m->fr.u[2]);
          double Bmax = bin + qmax * chord;
          if (p->sig_t[c] > 0.0) Bmax += qmax / p->sig_t[c];
          double BT = Bmax * HZ_MAXP_TOL;
          if (corner_amax(cf) > BT) {
            /* НЕПРЕРЫВНАЯ ФОРМА, НЕ ПЕРЕКЛЮЧАТЕЛЬ: жёсткий сброс в константу
             * дал на печи §699 предельный цикл (невязка стояла на 4.4e-2 при
             * дрожащем maxp 2346↔2348 — ровно предсказание А1125). Наклоны
             * сжимаются коэффициентом α, а cf0(α) = (rhs0 − α·kk)/A00 держит
             * строку 0 точно — тот же приём, что у штатного ограничителя
             * положительности. α ищется бисекцией: 30 шагов, детерминированно,
             * точность 1e-9 по α — не допуск качества, а шаг сетки поиска. */
            if (fabs(a0row[0]) > 0.0) {
              double kk2 = 0.0;
              for (int j = 1; j < 4; j++)
                kk2 += a0row[j] * cf[j];
              double S2 = fabs(cf[1]) + fabs(cf[2]) + fabs(cf[3]);
              double lo = 0.0, hi = 1.0;
              for (int b2 = 0; b2 < 30; b2++) {
                double mid = 0.5 * (lo + hi);
                double c0m = (rhs0 - mid * kk2) / a0row[0];
                if (fabs(c0m) + 0.5 * mid * S2 > BT)
                  hi = mid;
                else
                  lo = mid;
              }
              /* даже α = 0 может быть вне границы (балансное среднее выше
               * неё): поток дороже поточечной границы — остаёмся на α = 0 */
              cf[0] = (rhs0 - lo * kk2) / a0row[0];
              for (int j = 1; j < 4; j++)
                cf[j] *= lo;
            } else {
              cf[0] = cf[1] = cf[2] = cf[3] = 0.0;
            }
            nmaxp_last++;
          }
        }
        for (int j = 0; j < 4; j++)
          L[c * 4 + j] = cf[j];

        /* §733: ПЕЧАТЬ ВСКРЫТИЯ. Раскладка строки 0 (баланс потоков)
         * пересчитывается повторной прогулкой по граням и элементам — все
         * данные (L верховых, bout, sout, матрицы) в этой точке ещё живы.
         * Отношения печатаются и к |rhs0|, и к Σ|влётов| — rhs0 может быть
         * мал из-за компенсации влётов, и деление на него врёт (А1120). */
        if (isdump) {
          printf("    §733 ВСКРЫТИЕ ячейки %d, направление %d (%.3f, %.3f, %.3f), такт %d\n", c, mm,
                 om[0], om[1], om[2], it);
          for (int j = 0; j < 4; j++)
            printf("      §733 A[%d] = %14.6e %14.6e %14.6e %14.6e   rhs[%d] = %14.6e\n", j,
                   Adump[j][0], Adump[j][1], Adump[j][2], Adump[j][3], j, rdump[j]);
          double cmin = 1e300, cmax = -1e300;
          for (int k2 = 0; k2 < 8; k2++) {
            double v = cf[0];
            for (int a = 0; a < 3; a++)
              v += cf[a + 1] * (((k2 >> a) & 1) ? 0.5 : -0.5);
            if (v < cmin) cmin = v;
            if (v > cmax) cmax = v;
          }
          double r0 = 0.0;
          for (int i = 0; i < 4; i++)
            r0 += Adump[0][i] * cf[i];
          printf("      §733 решение cf = %.6e %.6e %.6e %.6e; κ(A) %.4g; углы полинома мин "
                 "%.6e макс %.6e; невязка строки 0 %.3e (отн. %.3e)\n",
                 cf[0], cf[1], cf[2], cf[3], mass_cond(Adump), cmin, cmax, fabs(r0 - rdump[0]),
                 fabs(rdump[0]) > 0.0 ? fabs(r0 - rdump[0]) / fabs(rdump[0]) : -1.0);
          double sumin = 0.0, sumout = 0.0, absin = 0.0, maxout = 0.0;
          for (int32_t k2 = m->fstart[c]; k2 < m->fstart[c + 1]; k2++) {
            int32_t f = m->flist[k2];
            int mia = m->f[f].ca == c;
            double sgn2;
            if (m->f[f].cb < 0)
              sgn2 = ((int)(~m->f[f].cb) & 1) ? 1.0 : -1.0;
            else
              sgn2 = mia ? 1.0 : -1.0;
            double on2 = om[m->f[f].axis] * sgn2;
            if (!(fabs(on2) > 0.0)) continue;
            const double (*fmine2)[4] = mia ? cu->ffm[f] : cu->ffmb[f];
            if (on2 > 0.0) {
              double v0 = 0.0;
              for (int i = 0; i < 4; i++)
                v0 += cf[i] * fmine2[i][0];
              v0 *= on2;
              int32_t dn = (m->f[f].cb < 0) ? -1 : (mia ? m->f[f].cb : m->f[f].ca);
              printf("      §733 НИЗОВАЯ  грань %7d ось %d сосед %7d%s: вылет %14.6e\n", f,
                     m->f[f].axis, dn, dn < 0 ? " (КУБ)" : ((cu->solid[dn]) ? " (СТЫК)" : ""), v0);
              sumout += v0;
              if (fabs(v0) > maxout) maxout = fabs(v0);
            } else {
              int32_t up2 = (m->f[f].cb < 0) ? -1 : (mia ? m->f[f].cb : m->f[f].ca);
              double v0 = 0.0;
              const char *kind = "";
              if (up2 < 0) {
                /* куб: на этом стенде wall_rho = NULL и влёт binc; повторяем
                 * формулу сборки для строки 0 */
                double s2 = (double)m->csize[c];
                double cen[3];
                for (int a = 0; a < 3; a++)
                  cen[a] = m->fr.o[a] + m->fr.u[a] * ((double)m->clo[c][a] + 0.5 * s2);
                double lb0 = p->binc0;
                for (int a = 0; a < 3; a++)
                  lb0 += p->binc[a] * (cen[a] - p->binx0[a]);
                double lb[4] = {lb0, 0.0, 0.0, 0.0};
                for (int a = 0; a < 3; a++)
                  lb[a + 1] = p->binc[a] * s2 * m->fr.u[a];
                for (int i = 0; i < 4; i++)
                  v0 += lb[i] * fmine2[i][0];
                v0 *= -on2;
                kind = " (КУБ)";
              } else if (cu->solid[up2]) {
                for (int i = 0; i < 4; i++)
                  v0 += bout[f * 4 + i] * fmine2[i][0];
                v0 *= -on2;
                kind = " (СТЫК)";
              } else {
                const double (*fx2)[4] = cu->ffmx[f];
                for (int i = 0; i < 4; i++)
                  v0 += L[up2 * 4 + i] * (mia ? fx2[0][i] : fx2[i][0]);
                v0 *= -on2;
              }
              printf("      §733 ВЕРХОВАЯ грань %7d ось %d сосед %7d%s: влёт  %14.6e\n", f,
                     m->f[f].axis, up2, kind, v0);
              sumin += v0;
              absin += fabs(v0);
            }
          }
          for (int32_t k2 = cu->sestart[c]; k2 < cu->sestart[c + 1]; k2++) {
            int32_t e = cu->selist[k2];
            const tr3_selem *se2 = &cu->se[e];
            double on2 = om[0] * se2->n[0] + om[1] * se2->n[1] + om[2] * se2->n[2];
            if (!(fabs(on2) > 0.0)) continue;
            if (on2 < 0.0) {
              double v0 = 0.0;
              for (int i = 0; i < 4; i++)
                v0 += cf[i] * se2->m[i][0];
              v0 *= -on2;
              printf("      §733 ЭЛЕМЕНТ-ВЫТОК %7d (площадь %.3e): вылет %14.6e\n", e, se2->area,
                     v0);
              sumout += v0;
              if (fabs(v0) > maxout) maxout = fabs(v0);
            } else {
              double v0 = 0.0;
              for (int i = 0; i < 4; i++)
                v0 += sout[e * 4 + i] * se2->m[i][0];
              v0 *= on2;
              printf("      §733 ЭЛЕМЕНТ-ВТОК  %7d (площадь %.3e): влёт  %14.6e\n", e, se2->area,
                     v0);
              sumin += v0;
              absin += fabs(v0);
            }
          }
          printf("    §733 СВОДКА ячейки %d: Σвлёт %14.6e, Σвылет %14.6e, Σ|влёт| %14.6e; МАКС "
                 "|вылет| %14.6e; макс/|rhs0| %.4g; макс/Σ|влёт| %.4g\n",
                 c, sumin, sumout, absin, maxout,
                 fabs(rdump[0]) > 0.0 ? maxout / fabs(rdump[0]) : -1.0,
                 absin > 0.0 ? maxout / absin : -1.0);
        }

        /* облучённость поверхностных элементов — из того же выточного члена */
        if (cu != NULL)
          for (int32_t k = cu->sestart[c]; k < cu->sestart[c + 1]; k++) {
            int32_t e = cu->selist[k];
            const tr3_selem *se = &cu->se[e];
            double on = om[0] * se->n[0] + om[1] * se->n[1] + om[2] * se->n[2];
            if (on >= 0.0) continue;
            for (int j = 0; j < 4; j++) {
              double acc = 0.0;
              for (int i = 0; i < 4; i++)
                acc += cf[i] * se->m[i][j];
              sinft[e * 4 + j] += (-on) * d->w[mm] * acc;
            }
            /* §735: граница облучённости — максимум мажоранты ИТОГОВОГО (после
             * клипа) радианса по направлениям, дающим вклад в элемент */
            double cm5 = corner_amax(cf);
            if (cm5 > semax[e]) semax[e] = cm5;
          }
      }

      for (int32_t c = 0; c < nc; c++)
        for (int j = 0; j < 4; j++)
          phit[c * 4 + j] += d->w[mm] * L[c * 4 + j];
      /* §731: ОБРАТНАЯ ПРОГУЛКА ОТ ЦЕЛИ. Угловое поле `L` направления в этой
       * точке полно; шаг — к наибольшему по `|L|` верховому соседу (прокси
       * вклада, А1113). Однопоточно, как и весь цикл (nth = 1). */
      if (chnc != NULL && chnl != NULL && fabs(L[(size_t)chain_tgt * 4]) > chn_bestlt) {
        chn_bestlt = fabs(L[(size_t)chain_tgt * 4]);
        chn_mm = mm;
        int32_t c2 = chain_tgt;
        int nch = 0;
        while (nch < HZ_CHAIN_MAX) {
          chnc[nch] = c2;
          chnl[nch] = fabs(L[(size_t)c2 * 4]);
          nch++;
          int32_t bu = -1;
          double bl = -1.0;
          for (int32_t q = m->fstart[c2]; q < m->fstart[c2 + 1]; q++) {
            int32_t f = m->flist[q];
            double on = om[m->f[f].axis];
            int32_t u = -1;
            if (on > 0.0) {
              if (m->f[f].cb == c2) u = m->f[f].ca;
            } else if (on < 0.0) {
              if (m->f[f].ca == c2) u = m->f[f].cb;
            }
            if (u < 0) continue;
            if (cu != NULL && cu->solid[u]) continue;
            double lu = fabs(L[(size_t)u * 4]);
            if (lu > bl) {
              bl = lu;
              bu = u;
            }
          }
          if (bu < 0) break;
          if (bl < HZ_CHAIN_FLOOR) {
            /* А1114: звено ниже пола записывается последним — скачок ЧЕРЕЗ пол
             * иначе потерялся бы вместе с обрывом */
            if (nch < HZ_CHAIN_MAX) {
              chnc[nch] = bu;
              chnl[nch] = bl;
              nch++;
            }
            break;
          }
          c2 = bu;
        }
        chn_n = nch;
      }
      for (int32_t f = 0; f < m->nf; f++) {
        /* §707: ЭТОТ ЦИКЛ ТЕПЕРЬ ОБСЛУЖИВАЕТ ДВА РОДА ГРАНИЦ, А НЕ ОДИН.
         * Прежде здесь были только грани КУБА (`cb < 0`), а поток, уходивший из
         * флюидной ячейки в СПЛОШНУЮ, не накапливался нигде и терялся — это и
         * есть измеренные `76 %` (§706). Стык «флюид — сплошное» есть такая же
         * стена: осевая грань с известной внешней нормалью, и потому обходится
         * тем же кодом с той же нормировкой `hsum`. */
        int32_t c;
        double on;
        const double (*fmm)[4];
        double fmmb[4][4];
        int at_solid = 0;
        if (m->f[f].cb < 0) {
          int wall0 = (int)(~m->f[f].cb);
          on = om[m->f[f].axis] * ((wall0 & 1) ? 1.0 : -1.0);
          if (!(on > 0.0)) continue;
          c = m->f[f].ca;
          if (cu != NULL) {
            fmm = cu->ffm[f];
          } else {
            double v[4][3];
            tr3_face_corners(m, f, v);
            face_mass2(m, (const double (*)[4][3]) & v, c, c, fmmb);
            fmm = fmmb;
          }
        } else if (cu != NULL && (cu->solid[m->f[f].ca] != cu->solid[m->f[f].cb])) {
          int a_solid = cu->solid[m->f[f].ca] ? 1 : 0;
          c = a_solid ? m->f[f].cb : m->f[f].ca;
          on = a_solid ? -om[m->f[f].axis] : om[m->f[f].axis];
          if (!(on > 0.0)) continue;
          fmm = a_solid ? cu->ffmb[f] : cu->ffm[f];
          at_solid = 1;
        } else {
          continue;
        }
        double accj[4] = {0, 0, 0, 0};
        for (int j = 0; j < 4; j++)
          for (int i = 0; i < 4; i++)
            accj[j] += L[c * 4 + i] * fmm[i][j];
        if (at_solid)
          psolid_acc += on * d->w[mm] * accj[0]; /* пришло В СТЕНУ */
        else
          pout_acc += on * d->w[mm] * accj[0]; /* вышло ИЗ ОБЛАСТИ */
        for (int j = 0; j < 4; j++)
          binft[f * 4 + j] += on * d->w[mm] * accj[j];
        /* §746: граница обновления bout — макс мажоранты L вкладчиков (как
         * semax у элементов; однопоточно, nth = 1) */
        double bm7 = corner_amax(&L[c * 4]);
        if (bm7 > bfmax[f]) bfmax[f] = bm7;
        /* ЭТАП C: у зеркальной грани копится момент ПО ОРДИНАТЕ, без веса —
         * зеркало не интегрирует по полусфере, оно переставляет направление. */
        if (!at_solid && mfid != NULL && mfid[f] >= 0)
          for (int j = 0; j < 4; j++)
            mspin[((size_t)mfid[f] * (size_t)nd + (size_t)mm) * 4 + (size_t)j] += accj[j];
      }
    }

    /* сведение потоковых накоплений */
    if (fail) {
      free(Lall);
      free(phinall);
      free(indegall);
      free(orderall);
      free(queueall);
      free(phin);
      free(bout);
      free(binf);
      free(sout);
      free(sinf);
      free(kadbg);
      free(gdbg);
      free(phiprev_dbg);
      free(chnc);
      free(chnl);
      free(hs_se);
      free(semax);
      free(bfmax);
      free(hs_out);
      free(bprev);
      free(sprev);
      free(mfid);
      free(mspec);
      free(mspin);
      free(mprev);
      free(binfall);
      free(sinfall);
      return 2;
    }
    for (int th = 0; th < nth; th++) {
      const double *pt = phinall + (size_t)th * (size_t)nc * 4;
      for (int32_t i = 0; i < nc * 4; i++)
        phin[i] += pt[i];
      const double *bt = binfall + (size_t)th * (size_t)m->nf * 4;
      for (int32_t i = 0; i < m->nf * 4; i++)
        binf[i] += bt[i];
      const double *stt = sinfall + (size_t)th * (size_t)(nse > 0 ? nse : 1) * 4;
      for (int32_t i = 0; i < nse * 4; i++)
        sinf[i] += stt[i];
    }
    st->pin = pin_acc;
    st->pout = pout_acc;
    /* §707: в баланс идёт ПОГЛОЩЁННОЕ стеной, а не весь пришедший поток:
     * доля `solid_rho` возвращается во флюид и учтена влётом с грани. */
    st->psolid = (1.0 - (p->solid_rho > 0.0 ? p->solid_rho : 0.0)) * psolid_acc;

    /* --- стенки и поверхности: новый исходящий радианс, DG1 по положению --- */
    if (p->wall_rho != NULL || (cu != NULL && p->solid_rho > 0.0))
      for (int32_t f = 0; f < m->nf; f++) {
        /* §707: та же развилка, что при накоплении, — грань КУБА или СТЫК. */
        int32_t cown;
        int wall = -1;
        double rho_w;
        if (m->f[f].cb < 0) {
          if (p->wall_rho == NULL) continue;
          wall = (int)(~m->f[f].cb);
          rho_w = p->wall_rho[wall];
          cown = m->f[f].ca;
        } else if (cu != NULL && p->solid_rho > 0.0 &&
                   (cu->solid[m->f[f].ca] != cu->solid[m->f[f].cb])) {
          int a_solid = cu->solid[m->f[f].ca] ? 1 : 0;
          cown = a_solid ? m->f[f].cb : m->f[f].ca;
          /* Номер «стенки» по ВНЕШНЕЙ нормали флюидной ячейки: у грани оси `a`
           * это `2a+1`, если материал с плюс-стороны, и `2a` — если с минус.
           * Та же нумерация, что у граней куба, поэтому `hsum` берётся готовым. */
          wall = 2 * m->f[f].axis + (a_solid ? 1 : 0);
          rho_w = p->solid_rho;
        } else {
          continue;
        }
        double fmm[4][4], rr[4], ee[4];
        if (cu != NULL) {
          memcpy(fmm, (cown == m->f[f].ca) ? cu->ffm[f] : cu->ffmb[f], sizeof fmm);
        } else {
          double v[4][3];
          tr3_face_corners(m, f, v);
          face_mass2(m, (const double (*)[4][3]) & v, cown, cown, fmm);
        }
        if (!(fmm[0][0] > 0.0)) continue;
        for (int j = 0; j < 4; j++)
          rr[j] = binf[f * 4 + j];
        /* ПРОЕКЦИЯ ОБЛУЧЁННОСТИ НА DG1 — В ТРЁХ НЕИЗВЕСТНЫХ, А НЕ В ЧЕТЫРЁХ (К39).
         * Здесь стояло `solve4(fmm, rr, ee)`, а матрица масс ГРАНИ вырождена
         * точно: грань плоская, и базис ячейки на ней связан. Измерено, чем это
         * оборачивалось: откат к константе срабатывал на 1536 гранях из 1536,
         * то есть DG1 по положению не работал НИ РАЗУ, а его исход менялся от
         * итерации к итерации и не давал развёртке сойтись. */
        double nul[4];
        tr3_face_null(m, f, cown, nul);
        /* ДВА РАЗНЫХ УСЛОВИЯ, И ТОЛЬКО ВТОРОЕ ЕСТЬ НЕЛИНЕЙНОСТЬ — К86.
         * Вырождение элемента есть свойство ГЕОМЕТРИИ: множество вырожденных
         * элементов от поля не зависит вовсе, и откат на них ЛИНЕЕН. А
         * `corner_min < 0` есть свойство ПОЛЯ, то есть переключатель, и он
         * обязан подчиняться `p->limiter` наравне с объёмным. Прежде он стоял
         * безусловно, и потому выключение ограничителя оператор линейным НЕ
         * ДЕЛАЛО: Крылов решал не ту систему и приходил к другому ответу
         * (невязка сошедшегося ответа Неймана под его оператором была 1.000). */
        int fb_geom_f = tr3_project_plane(fmm, rr, nul, ee) != 0;
        if (fb_geom_f || (p->limiter && corner_min(ee) < 0.0)) {
          nfb++;
          /* §677: КАКОГО РОДА ОТКАТ. Геометрический считается отдельно от
           * полевого — только второй способен дрожать от такта к такту, и
           * только он может давать предельный цикл. ЩЕПКА определяется по
           * флюидному объёму ЯЧЕЙКИ, к которой грань принадлежит с наветренной
           * стороны: порог `1 %` целой ячейки — не магический, он ровно тот, по
           * которому §675 назвал максимум `φ` сидящим на щепке. */
          if (fb_geom_f)
            nfb_fg++;
          else {
            nfb_fp++;
            if (cu != NULL && m->f[f].ca >= 0) {
              double vfl = cu->mvol[m->f[f].ca][0][0];
              double s3 = (double)m->csize[m->f[f].ca];
              double vcell = s3 * s3 * s3 * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
              if (vcell > 0.0 && vfl < 0.01 * vcell) nfb_thin++;
            }
          }
          /* ОТКАЗ В ЗАКРЫТУЮ СТОРОНУ: либо элемент выродился геометрически, либо
           * проекция ушла в минус. Тогда остаётся ТОЧНОЕ среднее по грани —
           * константа, представимая в любом случае. Порога в критерии нет. */
          ee[0] = binf[f * 4] / fmm[0][0];
          ee[1] = ee[2] = ee[3] = 0.0;
        }
        /* §746: ПРИНЦИП МАКСИМУМА ДЛЯ ГРАНИ — та же граница, что у элементов
         * (§735): E(x) ≤ hsum·max|L| по вкладчикам (bfmax). Проекция DG1 на
         * грани не ограничена так же, как на элементах; взрыв §714 (bout до
         * 4.9e+34) жил здесь. Форма непрерывная (А1125/А1128). */
        if (!p->maxp_off) {
          double bnd = bfmax[f] * hsum[wall] * HZ_MAXP_TOL;
          double S3 = 0.5 * (fabs(ee[1]) + fabs(ee[2]) + fabs(ee[3]));
          if (fabs(ee[0]) + S3 > bnd) {
            if (fabs(ee[0]) < bnd && S3 > 0.0) {
              double t3 = (bnd - fabs(ee[0])) / S3;
              for (int j = 1; j < 4; j++)
                ee[j] *= t3;
            } else {
              ee[0] = binf[f * 4] / fmm[0][0];
              ee[1] = ee[2] = ee[3] = 0.0;
            }
            nmaxpf_last++;
          }
        }
        for (int j = 0; j < 4; j++)
          bout[f * 4 + j] = (hsum[wall] > 0.0 ? rho_w * ee[j] / hsum[wall] : 0.0);
        if (m->f[f].cb < 0 && p->wall_emit != NULL) bout[f * 4] += p->wall_emit[wall];
      }

    /* ЭТАП C: ЗЕРКАЛЬНЫЕ ГРАНИ — СВЯЗЬ `m → m′` ТОЧНОЙ ПЕРЕСТАНОВКОЙ.
     *
     * Радианс, пришедший на грань по ординате `mm`, уходит с неё по ординате
     * `mir`, и никакой свёртки по полусфере здесь нет: зеркало не усредняет, а
     * переставляет. Интерполировать между ординатами запрещает К3, и она здесь
     * не нужна — набор замкнут относительно осевых плоскостей ТОЧНО.
     *
     * Проекция моментов на DG1 берётся та же, что у ламбертовой грани, вместе с
     * разделением К86: вырождение элемента — геометрия, положительность —
     * переключатель. */
    if (nmf > 0)
      for (int32_t f = 0; f < m->nf; f++) {
        if (mfid[f] < 0) continue;
        int wall = (int)(~m->f[f].cb);
        int ax = (int)m->f[f].axis;
        double fmmb[4][4];
        const double (*fmm)[4];
        if (cu != NULL) {
          fmm = cu->ffm[f];
        } else {
          double v[4][3];
          tr3_face_corners(m, f, v);
          face_mass2(m, (const double (*)[4][3]) & v, m->f[f].ca, m->f[f].ca, fmmb);
          fmm = fmmb;
        }
        if (!(fmm[0][0] > 0.0)) continue;
        double nul[4];
        tr3_face_null(m, f, m->f[f].ca, nul);
        for (int mm2 = 0; mm2 < nd; mm2++) {
          double on2 = (d->ox[mm2] * (ax == 0) + d->oy[mm2] * (ax == 1) + d->oz[mm2] * (ax == 2)) *
                       ((wall & 1) ? 1.0 : -1.0);
          if (!(on2 > 0.0)) continue; /* только УХОДЯЩИЕ от ячейки в стенку */
          double rr2[4], ee2[4];
          for (int j = 0; j < 4; j++)
            rr2[j] = mspin[((size_t)mfid[f] * (size_t)nd + (size_t)mm2) * 4 + (size_t)j];
          if (tr3_project_plane(fmm, rr2, nul, ee2) != 0 || (p->limiter && corner_min(ee2) < 0.0)) {
            ee2[0] = rr2[0] / fmm[0][0];
            ee2[1] = ee2[2] = ee2[3] = 0.0;
          }
          int mr = d->mir[(size_t)ax * (size_t)nd + (size_t)mm2];
          for (int j = 0; j < 4; j++)
            mspec[((size_t)mfid[f] * (size_t)nd + (size_t)mr) * 4 + (size_t)j] =
                p->wall_spec[wall] * ee2[j];
        }
      }
    for (int32_t e = 0; e < nse; e++) {
      const tr3_selem *se = &cu->se[e];
      double rho = (p->facet_rho != NULL && se->facet < p->nfacet) ? p->facet_rho[se->facet] : 0.0;
      double em =
          p->elem_emit != NULL
              ? p->elem_emit[e]
              : ((p->facet_emit != NULL && se->facet < p->nfacet) ? p->facet_emit[se->facet] : 0.0);
      double fmm[4][4], rr[4], ee[4];
      memcpy(fmm, se->m, sizeof fmm);
      for (int j = 0; j < 4; j++)
        rr[j] = sinf[e * 4 + j];
      if (!(fmm[0][0] > 0.0)) continue;
      /* та же правка К39: нулевой вектор у поверхностного элемента хранится
       * рядом с матрицей масс и выведен из уравнения плоскости (cut3.h) */
      /* К86, то же самое на поверхностных элементах */
      int fb_geom_e = tr3_project_plane(fmm, rr, se->nul, ee) != 0;
      /* §721: БЕЗРАЗМЕРНОЕ УСИЛЕНИЕ ПРОЕКТОРА, считается ПОСЛЕ проекции и ДО
       * отката: `g` показывает, во сколько раз решение DG1 больше безопасного
       * константного ответа `rr[0]/M[0][0]`, которым откат и пользуется. У
       * здорового элемента `g` порядка единицы. */
      if (p->trace > 0 && gdbg != NULL) {
        double nee = 0.0;
        for (int j = 0; j < 4; j++)
          nee += ee[j] * ee[j];
        nee = sqrt(nee);
        double safe = fabs(rr[0]) > 0.0 ? fabs(rr[0]) / fmm[0][0] : 0.0;
        gdbg[e] = safe > 0.0 ? nee / safe : 0.0;
      }
      if (fb_geom_e || (p->limiter && corner_min(ee) < 0.0)) {
        ee[0] = sinf[e * 4] / fmm[0][0];
        ee[1] = ee[2] = ee[3] = 0.0;
        nfb++;
        if (fb_geom_e)
          nfb_eg++;
        else {
          nfb_ep++;
          double vfl = cu->mvol[se->cell][0][0];
          double s3 = (double)m->csize[se->cell];
          double vcell = s3 * s3 * s3 * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
          if (vcell > 0.0 && vfl < 0.01 * vcell) nfb_thin++;
        }
      }
      /* §735: ПРИНЦИП МАКСИМУМА ДЛЯ ОБЛУЧЁННОСТИ: E(x) = ∫L|ω·n|dω ≤
       * hs_se · max|L| по вкладывавшим направлениям (semax). Проекция DG1 на
       * элементе в принципе не ограничена (§722: p99 усиления 92). Форма
       * непрерывная (А1125): наклоны сжимаются до границы; если и среднее
       * выше — константный откат, как у прочих отказов. */
      if (!p->maxp_off) {
        double bnd = semax[e] * hs_se[e] * HZ_MAXP_TOL;
        double S3 = 0.5 * (fabs(ee[1]) + fabs(ee[2]) + fabs(ee[3]));
        if (fabs(ee[0]) + S3 > bnd) {
          if (fabs(ee[0]) < bnd && S3 > 0.0) {
            double t3 = (bnd - fabs(ee[0])) / S3;
            for (int j = 1; j < 4; j++)
              ee[j] *= t3;
          } else {
            ee[0] = sinf[e * 4] / fmm[0][0];
            ee[1] = ee[2] = ee[3] = 0.0;
          }
          nmaxpe_last++;
        }
      }
      for (int j = 0; j < 4; j++)
        sout[e * 4 + j] = (hs_se[e] > 0.0 ? rho * ee[j] / hs_se[e] : 0.0);
      sout[e * 4] += em;
    }

    /* НЕВЯЗКА МЕРИТ ВСЁ СОСТОЯНИЕ ИТЕРАЦИИ, А НЕ ТОЛЬКО `φ` — К84.
     *
     * Состояние здесь есть тройка `(φ, bout, sout)`, и обновляются они ПО
     * ОЧЕРЕДИ: развёртка считает `φ` по СТАРЫМ поверхностям, и лишь потом
     * поверхности пересчитываются. Значит изменение, вошедшее через поверхности,
     * доходит до `φ` только на СЛЕДУЮЩЕЙ итерации, и критерий, глядящий на одно
     * `φ`, объявляет сходимость раньше времени.
     *
     * Поймано негативным контролем К76: при тёплом старте с изменённым светом
     * развёртка останавливалась на НУЛЕВОЙ итерации и возвращала СТАРОЕ решение
     * (`max|Δφ| = 5.04` при `|φ| = 55`, строго пропорционально возмущению), то
     * есть давала ложное ускорение `×47`. Без встроенной проверки «оба ответа
     * обязаны совпасть» это выглядело бы как блестящий результат.
     *
     * И ВТОРОЙ ДОВОД, СИЛЬНЕЕ ПЕРВОГО: КАРТИНКА ДЕЛАЕТСЯ ИЗ `bout`/`sout`, а не
     * из `φ`. Критерий, который на них не смотрит, не измеряет то, ради чего
     * решение и считается. */
    /* §748: демпфирование x ← (1−ω)·x + ω·S(x). Неподвижная точка та же при
     * любом ω > 0; гасится ПРЕДЕЛЬНЫЙ ЦИКЛ переключателей (диагноз §748:
     * невязка падала до ~1e-2 и дрожала скачками при nclip, дрожащем
     * тысячами). Невязка — по НЕДЕМПФИРОВАННОМУ шагу |S(x) − x| (А1164):
     * старое состояние в момент замера ещё лежит в phi/bprev/sprev/mprev. */
    const double w748 = p->relax > 0.0 ? p->relax : 1.0;
    resid = 0.0;
    for (int32_t i = 0; i < nc * 4; i++) {
      double dd = fabs(phin[i] - phi[i]);
      if (dd > resid) resid = dd;
      phi[i] = (1.0 - w748) * phi[i] + w748 * phin[i];
    }
    for (int32_t i = 0; i < m->nf * 4; i++) {
      double dd = fabs(bout[i] - bprev[i]);
      if (dd > resid) resid = dd;
      bout[i] = (1.0 - w748) * bprev[i] + w748 * bout[i];
      bprev[i] = bout[i];
    }
    /* §731: ПЕЧАТЬ ТРАССЫ — цепочка направления с максимальным `|L|` в цели.
     * «Разрезанность» — отношение флюидного объёма к коробке меньше единицы;
     * страховка 1e-9 покрывает разные порядки умножения при равных величинах. */
    if (p->trace > 0 && chnc != NULL && chn_n > 0) {
      const double cuteps = 1e-9;
      printf("    §731 ТРАССА к ячейке %d (такт %d): направление %d (%.3f, %.3f, %.3f), |L| в "
             "цели %.4g, звеньев %d\n",
             chain_tgt, it, chn_mm, d->ox[chn_mm], d->oy[chn_mm], d->oz[chn_mm], chn_bestlt, chn_n);
      double gmax = -1.0, glogsum = 0.0;
      int32_t gmaxc = -1;
      int64_t nratio = 0, namp = 0, namp_cut = 0, ncut_all = 0;
      for (int k = 0; k < chn_n; k++) {
        int32_t c2 = chnc[k];
        double s3 = (double)m->csize[c2];
        double vcell = s3 * s3 * s3 * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
        double frac = (cu != NULL && vcell > 0.0) ? cu->mvol[c2][0][0] / vcell : 1.0;
        int iscut = frac < 1.0 - cuteps;
        if (iscut) ncut_all++;
        double ratio = (k + 1 < chn_n && chnl[k + 1] > 0.0) ? chnl[k] / chnl[k + 1] : -1.0;
        printf("      §731 звено %3d: ячейка %7d  размер %3d  доля флюида %.4g  |L| %.6g  "
               "множитель %.6g\n",
               k, c2, m->csize[c2], frac, chnl[k], ratio);
        if (ratio > 0.0) {
          nratio++;
          glogsum += log(ratio);
          if (ratio > gmax) {
            gmax = ratio;
            gmaxc = c2;
          }
          if (ratio > HZ_CHAIN_AMP) {
            namp++;
            if (iscut) namp_cut++;
          }
        }
      }
      printf("    §731 СВОДКА: звеньев %d; МАКС МНОЖИТЕЛЬ %.6g (ячейка %d); геом. среднее "
             "%.4g по %lld шагам; шагов > %.2g: %lld, из них РАЗРЕЗАННЫХ %lld; разрезанных "
             "вдоль всей цепочки %lld из %d\n",
             chn_n, gmax, gmaxc, nratio > 0 ? exp(glogsum / (double)nratio) : -1.0,
             (long long)nratio, HZ_CHAIN_AMP, (long long)namp, (long long)namp_cut,
             (long long)ncut_all, chn_n);
      chn_bestlt = -1.0;
      chn_n = 0;
      chn_mm = -1;
    }
    /* §725: НАБЛЮДАЕМЫЙ КОЭФФИЦИЕНТ УСИЛЕНИЯ ЭЛЕМЕНТА ЗА ТАКТ. Величина уже
     * течёт через развёртку: `sout` этого такта против `sprev` прошлого.
     * Считается ДО перезаписи `sprev` в блоке невязки. */
    if (p->trace > 0 && nse > 0) {
      double *gi = malloc((size_t)nse * sizeof *gi);
      if (gi != NULL) {
        int64_t ng2 = 0, ngt = 0;
        double gmx2 = 0.0;
        int32_t igm = -1;
        for (int32_t e2 = 0; e2 < nse; e2++) {
          double a = fabs(sprev[e2 * 4]), b = fabs(sout[e2 * 4]);
          if (!(a > 0.0)) continue;
          double g2 = b / a;
          gi[ng2++] = g2;
          if (g2 > 1.0) ngt++;
          if (g2 > gmx2) {
            gmx2 = g2;
            igm = e2;
          }
        }
        if (ng2 > 0) {
          qsort(gi, (size_t)ng2, sizeof *gi, cmp_dbl_dbg);
          printf("    §725 УСИЛЕНИЕ ЭЛЕМЕНТА ЗА ТАКТ %d: элементов %lld, медиана %.4g, p99 %.4g, "
                 "МАКСИМУМ %.4g (элемент %d); с g > 1: %lld (%.2f %%); у 98531 g = %.4g\n",
                 it, (long long)ng2, gi[ng2 / 2], gi[(ng2 * 99) / 100], gmx2, igm, (long long)ngt,
                 100.0 * (double)ngt / (double)ng2,
                 (98531 < nse && fabs(sprev[98531 * 4]) > 0.0)
                     ? fabs(sout[98531 * 4]) / fabs(sprev[98531 * 4])
                     : -1.0);
        }
        free(gi);
      }
    }
    /* §727: ЧТО ОБЩЕГО У ХВОСТА. Медианы шести уже стоящих величин по
     * множеству S (усиление за такт выше порога), по всей популяции §725 и по
     * контрольной выборке того же размера — без неё различия не читаются
     * (А1063), а если различия появятся и у неё, врёт сам отбор. Считается ДО
     * перезаписи `sprev`; `kadbg` в этой точке ещё держит значения ТЕКУЩЕГО
     * такта (обнуляется ниже, в блоке §723). */
    if (p->trace > 0 && nse > 0 && cu != NULL && gdbg != NULL && kadbg != NULL) {
      /* §726: медиана g_it 1.05, хвост 6e+06…1.1e+07. Два порога: 100 —
       * исходное задание (отстоит от обеих мод на два порядка), 1e+06 выделяет
       * сам «миллионный» хвост. А1101: порог 100 ловит ещё и широкое плечо
       * (10…24 % популяции), разбавляющее медианный портрет. */
      static const double gthr2[2] = {100.0, 1e+06};
      enum { NQ = 6 };
      double *q6 = malloc((size_t)NQ * (size_t)nse * sizeof *q6);
      int32_t *pe = malloc((size_t)nse * sizeof *pe);
      int32_t *te = malloc((size_t)nse * sizeof *te);
      int32_t *ke = malloc((size_t)nse * sizeof *ke);
      int32_t *cbuf = malloc((size_t)nse * sizeof *cbuf);
      double *scr = malloc((size_t)nse * sizeof *scr);
      if (q6 != NULL && pe != NULL && te != NULL && ke != NULL && cbuf != NULL && scr != NULL) {
        double h2 = m->fr.u[0] * m->fr.u[1];
        int64_t npop = 0;
        for (int32_t e2 = 0; e2 < nse; e2++) {
          const tr3_selem *se2 = &cu->se[e2];
          int32_t c2 = se2->cell;
          double *qq = q6; /* строка q — величина, столбец — элемент */
          qq[0 * (size_t)nse + (size_t)e2] = se2->area / h2;
          qq[1 * (size_t)nse + (size_t)e2] = hs_se[e2];
          qq[2 * (size_t)nse + (size_t)e2] = gdbg[e2];
          double ffr3 = -1.0, kmv = -1.0, ka = -1.0;
          if (c2 >= 0 && c2 < nc) {
            double s3 = (double)m->csize[c2];
            double vcell = s3 * s3 * s3 * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
            ffr3 = vcell > 0.0 ? cu->mvol[c2][0][0] / vcell : -1.0;
            kmv = mass_cond(cu->mvol[c2]);
            ka = kadbg[c2];
          }
          qq[3 * (size_t)nse + (size_t)e2] = ffr3;
          qq[4 * (size_t)nse + (size_t)e2] = kmv;
          qq[5 * (size_t)nse + (size_t)e2] = ka;
          if (fabs(sprev[e2 * 4]) > 0.0) pe[npop++] = e2;
        }
        for (int t2 = 0; t2 < 2 && npop > 0; t2++) {
          const double gtail = gthr2[t2];
          int64_t ntl = 0;
          for (int64_t i = 0; i < npop; i++) {
            int32_t e2 = pe[i];
            if (fabs(sout[e2 * 4]) / fabs(sprev[e2 * 4]) > gtail) te[ntl++] = e2;
          }
          if (ntl == 0) continue;
          /* А1100: контроль — хеш-детерминированный, НЕ стрижка. Элементы
           * одной ячейки лежат в массиве подряд, и выборка «каждый k-й» при
           * k ≥ 2 не может взять двух из одной ячейки — её доля различных
           * ячеек равна 1 по построению, как база кластеризации она слепа.
           * Кнутов множитель даёт воспроизводимую псевдослучайную выборку
           * ожидаемого размера ntl без rand(). */
          int64_t nk2 = 0;
          for (int64_t i = 0; i < npop; i++) {
            uint32_t hsh = (uint32_t)pe[i] * 2654435761u;
            if ((int64_t)(hsh % (uint32_t)npop) < ntl) ke[nk2++] = pe[i];
          }
          double mt[NQ], mp[NQ], mk[NQ];
          for (int q = 0; q < NQ; q++) {
            const double *vq = q6 + (size_t)q * (size_t)nse;
            mt[q] = med_idx_dbg(vq, te, ntl, scr);
            mp[q] = med_idx_dbg(vq, pe, npop, scr);
            mk[q] = nk2 > 0 ? med_idx_dbg(vq, ke, nk2, scr) : -1.0;
          }
          for (int64_t i = 0; i < ntl; i++)
            cbuf[i] = cu->se[te[i]].cell;
          int64_t uct = nuniq_i32_dbg(cbuf, ntl);
          for (int64_t i = 0; i < nk2; i++)
            cbuf[i] = cu->se[ke[i]].cell;
          int64_t uck = nuniq_i32_dbg(cbuf, nk2);
          printf("    §727 ХВОСТ g > %g (такт %d): элементов %lld из %lld (%.2f %%), различных "
                 "ячеек %lld из %lld (контроль: %lld из %lld)\n",
                 gtail, it, (long long)ntl, (long long)npop, 100.0 * (double)ntl / (double)npop,
                 (long long)uct, (long long)ntl, (long long)uck, (long long)nk2);
          printf("    §727   медианы ХВОСТ/ВСЕ/КОНТРОЛЬ: площадь h² %.4g/%.4g/%.4g; hs_se "
                 "%.4g/%.4g/%.4g; g_proj %.4g/%.4g/%.4g; доля флюида %.4g/%.4g/%.4g; κ(mvol) "
                 "%.4g/%.4g/%.4g; κ(A) %.4g/%.4g/%.4g\n",
                 mt[0], mp[0], mk[0], mt[1], mp[1], mk[1], mt[2], mp[2], mk[2], mt[3], mp[3], mk[3],
                 mt[4], mp[4], mk[4], mt[5], mp[5], mk[5]);
        }
      }
      free(q6);
      free(pe);
      free(te);
      free(ke);
      free(cbuf);
      free(scr);
    }
    for (int32_t i = 0; i < nse * 4; i++) {
      double dd = fabs(sout[i] - sprev[i]);
      if (dd > resid) resid = dd;
      sout[i] = (1.0 - w748) * sprev[i] + w748 * sout[i]; /* §748 */
      sprev[i] = sout[i];
    }
    /* ЗЕРКАЛЬНОЕ хранимое — тоже часть состояния (К84), и без него критерий
     * объявил бы сходимость, пока зеркала ещё не установились. */
    for (int32_t i = 0; i < nmf * nd * 4; i++) {
      double dd = fabs(mspec[i] - mprev[i]);
      if (dd > resid) resid = dd;
      mspec[i] = (1.0 - w748) * mprev[i] + w748 * mspec[i]; /* §748 */
      mprev[i] = mspec[i];
    }
    /* ИСТОРИЯ НЕВЯЗКИ ПЕЧАТАЕТСЯ ПО ТРЕБОВАНИЮ, И ЭТО НЕ ОТЛАДКА (К38).
     *
     * «Не сошлось за N итераций» и «сходится медленно» — РАЗНЫЕ диагнозы, и
     * различает их не число итераций, а история невязки. Один раз эта разница
     * уже стоила проекту целого пункта плана: 4000 итераций были прочитаны как
     * медленная сходимость и приписаны отсутствию DSA, а оказались
     * автоколебанием нелинейного отката поверх вырожденной проекции. Тогда
     * историю печатали разово и руками; здесь она есть штатный вывод.
     *
     * Вместе с невязкой печатаются `nclip` и `nfallback` ПОСЛЕДНЕЙ итерации:
     * если невязка стоит, а эти два числа МЕНЯЮТСЯ — это переключатель, а не
     * скорость, и лечится оно не разгоном. */
    if (p->trace > 0 && (it % p->trace == 0 || resid < tol))
      printf("    it %5d  resid %.3e  nclip %d  maxp %d+%d+%d  nfb %d = грани(геом %d, поле %d) + "
             "элементы(геом %d, поле %d); из полевых на ЩЕПКАХ %d\n",
             it, resid, nclip_last, nmaxp_last, nmaxpe_last, nmaxpf_last, nfb, nfb_fg, nfb_fp,
             nfb_eg, nfb_ep, nfb_thin);
    /* §693: РОСТ КАК ФУНКЦИЯ ДОЛИ ФЛЮИДА, БЕЗ ПОРОГОВ. Декада — способ показать
     * кривую: каждая ячейка попадает ровно в одну, ни одна не отбрасывается.
     * Печатается медиана (а не среднее — А1090), p99 и максимум. */
    if (p->trace > 0 && cu != NULL && phiprev_dbg != NULL) {
      typedef struct {
        double frac, g, k;
        int32_t ci;
      } grec;
      grec *gr = malloc((size_t)nc * sizeof *gr);
      if (gr != NULL) {
        int64_t ng = 0;
        for (int32_t ci = 0; ci < nc; ci++) {
          double a = fabs(phiprev_dbg[ci]), b = fabs(phi[4 * (size_t)ci]);
          if (!(a > 0.0)) continue;
          double s3 = (double)m->csize[ci];
          double vcell = s3 * s3 * s3 * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
          if (!(vcell > 0.0)) continue;
          gr[ng].frac = cu->mvol[ci][0][0] / vcell;
          gr[ng].g = b / a;
          gr[ng].k = mass_cond(cu->mvol[ci]);
          gr[ng].ci = ci;
          ng++;
        }
        double *tmp = malloc((size_t)(ng > 0 ? ng : 1) * sizeof *tmp);
        if (tmp != NULL && ng > 0) {
          for (int64_t i = 0; i < ng; i++)
            tmp[i] = gr[i].g;
          qsort(tmp, (size_t)ng, sizeof *tmp, cmp_dbl_dbg);
          printf("    §693 РОСТ ПО СЦЕНЕ (такт %d): ячеек %lld, МЕДИАНА %.4g, p99 %.4g, макс "
                 "%.4g\n",
                 it, (long long)ng, tmp[ng / 2], tmp[(ng * 99) / 100], tmp[ng - 1]);
          printf("    §693 ЗАВИСИМОСТЬ ОТ ДОЛИ ФЛЮИДА (декады, без порогов):\n");
          for (int dec = 9; dec >= 0; dec--) {
            double lo = (dec == 0) ? 0.0 : pow(10.0, -(double)dec);
            double hi = (dec == 0) ? 1.0e30 : pow(10.0, -(double)(dec - 1));
            int64_t cnt = 0;
            for (int64_t i = 0; i < ng; i++)
              if (gr[i].frac >= lo && gr[i].frac < hi) tmp[cnt++] = gr[i].g;
            if (cnt == 0) continue;
            qsort(tmp, (size_t)cnt, sizeof *tmp, cmp_dbl_dbg);
            double kmed = 0.0;
            {
              int64_t c2 = 0;
              double *tk = malloc((size_t)cnt * sizeof *tk);
              if (tk != NULL) {
                for (int64_t i = 0; i < ng; i++)
                  if (gr[i].frac >= lo && gr[i].frac < hi) tk[c2++] = gr[i].k;
                qsort(tk, (size_t)c2, sizeof *tk, cmp_dbl_dbg);
                kmed = tk[c2 / 2];
                free(tk);
              }
            }
            printf("      доля [%.0e, %.0e): ячеек %7lld  медиана роста %.4g  p99 %.4g  макс "
                   "%.4g  медиана κ %.4g\n",
                   lo, hi > 1.0e29 ? 1.0 : hi, (long long)cnt, tmp[cnt / 2], tmp[(cnt * 99) / 100],
                   tmp[cnt - 1], kmed);
          }
        }
        /* §695: СВОЙСТВА НАСЕЛЕНИЯ ВЫШЕ `p99` ПО РОСТУ. Популяция задана РАНГОМ,
         * а не порогом (урок А1089): берутся ячейки, чей рост больше значения
         * `p99` того же распределения. Те же величины считаются по ВСЕЙ сцене —
         * без знаменателя доля не читается (А1063), и совпадение сценной доли с
         * уже напечатанными 22.19 % служит проверкой прибора. */
        if (ng > 0) {
          double gp99 = tmp[(ng * 99) / 100];
          for (int pass2 = 0; pass2 < 2; pass2++) {
            int64_t ncel = 0, nbnd = 0, fj = 0, fa = 0;
            double *fr2 = malloc((size_t)ng * sizeof *fr2);
            double *kk2 = malloc((size_t)ng * sizeof *kk2);
            double *nf2 = malloc((size_t)ng * sizeof *nf2);
            if (fr2 == NULL || kk2 == NULL || nf2 == NULL) {
              free(fr2);
              free(kk2);
              free(nf2);
              break;
            }
            for (int64_t i = 0; i < ng; i++) {
              if (pass2 == 0 && !(gr[i].g > gp99)) continue;
              int32_t ci = gr[i].ci;
              int bnd = 0, nfc = 0;
              for (int32_t q = m->fstart[ci]; q < m->fstart[ci + 1]; q++) {
                int32_t fi = m->flist[q];
                int32_t nb = (m->f[fi].ca == ci) ? m->f[fi].cb : m->f[fi].ca;
                nfc++;
                fa++;
                if (nb < 0)
                  bnd = 1;
                else if (m->csize[nb] != m->csize[ci])
                  fj++;
              }
              fr2[ncel] = gr[i].frac;
              kk2[ncel] = gr[i].k;
              nf2[ncel] = (double)nfc;
              ncel++;
              if (bnd) nbnd++;
            }
            if (ncel > 0) {
              qsort(fr2, (size_t)ncel, sizeof *fr2, cmp_dbl_dbg);
              qsort(kk2, (size_t)ncel, sizeof *kk2, cmp_dbl_dbg);
              qsort(nf2, (size_t)ncel, sizeof *nf2, cmp_dbl_dbg);
              printf("    §695 %s: ячеек %lld; граней %lld, с ПЕРЕПАДОМ %lld (%.2f %%); медиана "
                     "граней %.0f; на границе области %lld (%.2f %%); медиана доли флюида %.4g; "
                     "медиана κ %.4g\n",
                     pass2 == 0 ? "ВЫШЕ p99 ПО РОСТУ" : "ВСЯ СЦЕНА (тот же счёт)", (long long)ncel,
                     (long long)fa, (long long)fj, 100.0 * (double)fj / (double)(fa ? fa : 1),
                     nf2[ncel / 2], (long long)nbnd, 100.0 * (double)nbnd / (double)ncel,
                     fr2[ncel / 2], kk2[ncel / 2]);
            }
            free(fr2);
            free(kk2);
            free(nf2);
          }
        }
        free(tmp);
        free(gr);
      }
      /* ПЕРЕПАД УРОВНЕЙ У ЯЧЕЙКИ: сколько её граней имеют соседа ДРУГОГО
       * размера. Печатается и по сцене — иначе сравнивать не с чем (А1063). */
      int64_t fall = 0, fjump = 0;
      for (int32_t ci = 0; ci < nc; ci++)
        for (int32_t q = m->fstart[ci]; q < m->fstart[ci + 1]; q++) {
          int32_t fi = m->flist[q];
          int32_t nb = (m->f[fi].ca == ci) ? m->f[fi].cb : m->f[fi].ca;
          fall++;
          if (nb >= 0 && m->csize[nb] != m->csize[ci]) fjump++;
        }
      printf("    §693 ПЕРЕПАД УРОВНЕЙ ПО СЦЕНЕ: граней %lld, из них с соседом другого размера "
             "%lld (%.2f %%)\n",
             (long long)fall, (long long)fjump, 100.0 * (double)fjump / (double)(fall ? fall : 1));
    }
    /* §690 (СЛЕДСТВИЕ): ГДЕ ИМЕННО РАСТЁТ ПОЛЕ. Взрыв на шесть порядков за такт
     * обязан где-то сидеть; если он размазан — виноват не разрез. Печатается
     * десятка худших вместе со свойствами ячейки, чтобы причина и следствие
     * стояли в одной строке. */
    if (p->trace > 0 && cu != NULL && phiprev_dbg != NULL) {
      double gtop[10];
      int32_t itop[10];
      for (int k = 0; k < 10; k++) {
        gtop[k] = -1.0;
        itop[k] = -1;
      }
      double gsum = 0.0;
      int64_t ngr = 0;
      for (int32_t ci = 0; ci < nc; ci++) {
        double a = fabs(phiprev_dbg[ci]), b = fabs(phi[4 * (size_t)ci]);
        if (!(a > 0.0)) continue;
        double g = b / a;
        gsum += g;
        ngr++;
        for (int k = 0; k < 10; k++)
          if (g > gtop[k]) {
            for (int q = 9; q > k; q--) {
              gtop[q] = gtop[q - 1];
              itop[q] = itop[q - 1];
            }
            gtop[k] = g;
            itop[k] = ci;
            break;
          }
      }
      printf("    §690 РОСТ ЗА ТАКТ %d: средний %.3e по %lld ячейкам; ДЕСЯТКА ХУДШИХ:\n", it,
             ngr ? gsum / (double)ngr : 0.0, (long long)ngr);
      for (int k = 0; k < 10 && itop[k] >= 0; k++) {
        int32_t ci = itop[k];
        double s3 = (double)m->csize[ci];
        double vcell = s3 * s3 * s3 * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
        double vfl = cu->mvol[ci][0][0];
        printf("      рост %.3e  ячейка %d  размер %d  доля флюида %.3e  κ(M) %.4g\n", gtop[k], ci,
               m->csize[ci], vcell > 0.0 ? vfl / vcell : -1.0, mass_cond(cu->mvol[ci]));
      }
    }
    if (phiprev_dbg != NULL)
      for (int32_t ci = 0; ci < nc; ci++)
        phiprev_dbg[ci] = phi[4 * (size_t)ci];
    /* §723: распределение обусловленности матрицы обновления ячейки. */
    if (p->trace > 0 && kadbg != NULL) {
      double *kv = malloc((size_t)nc * sizeof *kv);
      if (kv != NULL) {
        int64_t nk = 0;
        double kmx = 0.0;
        int32_t ikmx = -1;
        for (int32_t ci = 0; ci < nc; ci++) {
          if (!(kadbg[ci] > 0.0)) continue;
          kv[nk++] = kadbg[ci];
          if (kadbg[ci] > kmx) {
            kmx = kadbg[ci];
            ikmx = ci;
          }
        }
        if (nk > 0) {
          qsort(kv, (size_t)nk, sizeof *kv, cmp_dbl_dbg);
          printf("    §723 κ(A) ОБНОВЛЕНИЯ ЯЧЕЙКИ (такт %d): ячеек %lld, медиана %.4g, p99 %.4g, "
                 "МАКСИМУМ %.4g (ячейка %d); у ячейки 81060 κ(A) = %.4g\n",
                 it, (long long)nk, kv[nk / 2], kv[(nk * 99) / 100], kmx, ikmx,
                 (81060 < nc) ? kadbg[81060] : -1.0);
        }
        free(kv);
      }
      memset(kadbg, 0, (size_t)nc * sizeof *kadbg);
    }
    /* §721: распределение усиления проектора по элементам. Медиана печатается
     * рядом с максимумом — без неё «максимум велик» не читается (А1063). */
    if (p->trace > 0 && gdbg != NULL && nse > 0) {
      double *gg = malloc((size_t)nse * sizeof *gg);
      if (gg != NULL) {
        int64_t ngg = 0;
        double gmx = 0.0;
        int32_t igmx = -1;
        for (int32_t e2 = 0; e2 < nse; e2++) {
          if (!(gdbg[e2] > 0.0)) continue;
          gg[ngg++] = gdbg[e2];
          if (gdbg[e2] > gmx) {
            gmx = gdbg[e2];
            igmx = e2;
          }
        }
        if (ngg > 0) {
          qsort(gg, (size_t)ngg, sizeof *gg, cmp_dbl_dbg);
          printf("    §721 УСИЛЕНИЕ ПРОЕКТОРА (такт %d): элементов %lld, медиана %.4g, p99 %.4g, "
                 "МАКСИМУМ %.4g на элементе %d; у элемента 98531 g = %.4g\n",
                 it, (long long)ngg, gg[ngg / 2], gg[(ngg * 99) / 100], gmx, igmx,
                 (98531 < nse) ? gdbg[98531] : -1.0);
        }
        free(gg);
      }
    }
    /* §713: ГДЕ СИДИТ МАКСИМУМ. `resid` есть макс-норма по ТРЁМ носителям сразу,
     * и по ней нельзя сказать, что именно взорвалось. Разделяем и печатаем
     * врозь, с адресом худшего. */
    if (p->trace > 0) {
      double mphi = 0.0, mbout = 0.0, msout = 0.0;
      int32_t iphi = -1, ibout = -1, isout = -1;
      for (int32_t i = 0; i < nc * 4; i++)
        if (fabs(phin[i]) > mphi) {
          mphi = fabs(phin[i]);
          iphi = i / 4;
        }
      for (int32_t i = 0; i < m->nf * 4; i++)
        if (fabs(bout[i]) > mbout) {
          mbout = fabs(bout[i]);
          ibout = i / 4;
        }
      for (int32_t i = 0; i < nse * 4; i++)
        if (fabs(sout[i]) > msout) {
          msout = fabs(sout[i]);
          isout = i / 4;
        }
      double ffr = -1.0;
      if (iphi >= 0 && cu != NULL) {
        double s3 = (double)m->csize[iphi];
        double V = s3 * s3 * s3 * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
        ffr = V > 0.0 ? cu->mvol[iphi][0][0] / V : -1.0;
      }
      const char *bkind = "нет";
      if (ibout >= 0)
        bkind = (m->f[ibout].cb < 0)
                    ? "ГРАНЬ КУБА"
                    : ((cu != NULL && (cu->solid[m->f[ibout].ca] || cu->solid[m->f[ibout].cb]))
                           ? "СТЫК"
                           : "внутренняя");
      /* §719: СВОЙСТВА ХУДШИХ ОБЪЕКТОВ, А НЕ ТОЛЬКО ИХ НОМЕРА. Рядом печатается
       * ТИПИЧНЫЙ элемент (медианный по площади) — иначе «κ = 1e4» не с чем
       * сравнить (А1063). */
      if (isout >= 0 && cu != NULL) {
        double h2 = m->fr.u[0] * m->fr.u[1];
        const tr3_selem *sw = &cu->se[isout];
        double kk = mass_cond(sw->m);
        int32_t cw = sw->cell;
        double ffr2 = -1.0, kc = -1.0;
        int nfc2 = 0;
        if (cw >= 0 && cw < nc) {
          double s3 = (double)m->csize[cw];
          double V = s3 * s3 * s3 * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
          ffr2 = V > 0.0 ? cu->mvol[cw][0][0] / V : -1.0;
          kc = mass_cond(cu->mvol[cw]);
          nfc2 = (int)(m->fstart[cw + 1] - m->fstart[cw]);
        }
        /* типичный элемент: берётся первый, чья площадь ближе всего к медиане
         * из уже напечатанного распределения — здесь просто площадь около 0.124 h² */
        int32_t itypic = -1;
        double bestd = 1e300;
        for (int32_t e2 = 0; e2 < nse; e2++) {
          if (!(cu->se[e2].area > 0.0)) continue;
          double d2 = fabs(cu->se[e2].area / h2 - 0.124);
          if (d2 < bestd) {
            bestd = d2;
            itypic = e2;
          }
        }
        printf("    §719 ХУДШИЙ ЭЛЕМЕНТ %d: площадь %.4e h², hs_se %.4f, κ(M_e) %.4g, ячейка %d "
               "(доля флюида %.4g, граней %d, κ(mvol) %.4g)\n",
               isout, sw->area / h2, hs_se[isout], kk, cw, ffr2, nfc2, kc);
        if (itypic >= 0)
          printf("    §719 ТИПИЧНЫЙ ЭЛЕМЕНТ %d (для сравнения): площадь %.4e h², hs_se %.4f, "
                 "κ(M_e) %.4g\n",
                 itypic, cu->se[itypic].area / h2, hs_se[itypic], mass_cond(cu->se[itypic].m));
      }
      printf("    §713 ТАКТ %d, МАКСИМУМ ПО НОСИТЕЛЯМ: |φ| %.4e (ячейка %d, доля флюида %.4g); "
             "|bout| %.4e (грань %d, %s); |sout| %.4e (элемент %d)\n",
             it, mphi, iphi, ffr, mbout, ibout, bkind, msout, isout);
    }
    /* §709 (Р2): СРЕДНЕЕ ПРОТИВ НАКЛОНОВ. Баланс проверяет только нулевой момент
     * (К12), наклоны переносят ноль энергии и потому невидимы для него. Если
     * растут именно они — расходимость сидит в DG1, а не в переносе. */
    if (p->trace > 0) {
      double m0 = 0.0, m1 = 0.0;
      double wf = -1.0;
      int32_t wc2 = -1;
      for (int32_t ci = 0; ci < nc; ci++) {
        double a0 = fabs(phi[4 * (size_t)ci]);
        if (a0 > m0) m0 = a0;
        for (int j = 1; j < 4; j++) {
          double aj = fabs(phi[4 * (size_t)ci + (size_t)j]);
          if (aj > m1) {
            m1 = aj;
            wc2 = ci;
          }
        }
      }
      if (wc2 >= 0 && cu != NULL) {
        double s3 = (double)m->csize[wc2];
        double V = s3 * s3 * s3 * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
        wf = V > 0.0 ? cu->mvol[wc2][0][0] / V : -1.0;
      }
      printf("    §709 ТАКТ %d: max|СРЕДНЕЕ| %.4e, max|НАКЛОН| %.4e, отношение %.4e; худший по "
             "наклону — ячейка %d, доля флюида %.4g\n",
             it, m0, m1, m0 > 0.0 ? m1 / m0 : -1.0, wc2, wf);
    }
    if (resid < tol) break;

    /* ОСТАНОВКА ПО ЗАСТОЮ, А НЕ ТОЛЬКО ПО ДОПУСКУ — К81.
     *
     * Допуск здесь АБСОЛЮТНЫЙ, а `φ` имеет масштаб задачи (в комнате ~16), то
     * есть `1e-9` есть `6e-11` относительных. У невязки при этом ЕСТЬ ПОЛ:
     * измерено, что на сетке 64³ она падает геометрически до `4e-8` и дальше
     * болтается немонотонно, тогда как решение уже стоит в десятом знаке.
     * Критерий тогда не выполним НИКОГДА, и цикл упирается в `maxit`: 4000
     * итераций вместо пятидесяти, 4 ч 41 мин вместо четырёх минут.
     *
     * ПОЧЕМУ ЭТО НЕ МЕДЛЕННАЯ СХОДИМОСТЬ, А ПОЛ (различает их история, К38):
     * скорость та же 0.63, что на 16³ и 32³, а `nclip` и `nfallback` ПОСТОЯННЫ
     * с двадцатой итерации — то есть нелинейный переключатель стоит, и дело не
     * в нём. Механизм самого пола НЕ УСТАНОВЛЕН и записан открытым вопросом.
     *
     * КРИТЕРИЙ БЕЗ МАГИЧЕСКОГО ПОРОГА: следим не за ВЕЛИЧИНОЙ невязки, а за тем,
     * УЛУЧШАЕТСЯ ли она. `HZ_SWEEP_STALL` итераций подряд без улучшения лучшего
     * достигнутого значения означают пол, и вот почему это выведено, а не
     * подобрано: ряд по рассеянию сходится как `ρⁿ`, и даже при `ρ = 0.9`
     * десять итераций обязаны дать множитель `0.35`. Ноль улучшения за десять
     * итераций медленной сходимостью быть не может. */
    if (resid < best) {
      best = resid;
      nstall = 0;
    } else if (++nstall >= HZ_SWEEP_STALL) {
      st->stalled = 1;
      break;
    }
  }

  /* ПОГЛОЩЕНИЕ БЕРЁТ ВСЕ ЧЕТЫРЕ МОМЕНТА, А НЕ ОДНО СРЕДНЕЕ (К48).
   *
   * В уравнении при тестовой функции `v = 1` член поглощения есть
   * `∫(σ_t−σ_s)·L dV = (σ_t−σ_s)·Σ_i φ_i·∫b_i dV`, то есть `Σ_i φ_i·M[0][i]`.
   * У КОРОБКИ `∫b_i dV = 0` при `i ≥ 1` по симметрии, и остаётся ровно
   * `φ_0·объём` — так здесь и было написано. Но У РАЗРЕЗАННОЙ ЯЧЕЙКИ СИММЕТРИИ
   * НЕТ, и `M[0][i]` при `i ≥ 1` не ноль (это уже отмечено выше, в члене
   * переноса). Наклонные члены выбрасывались, и баланс терял 4.5e-9 от
   * втекшего.
   *
   * Признак, по которому это нашлось: без ограничителя баланс был ХУЖЕ в 800
   * раз (3.6e-6 против 4.5e-9). Ограничитель гасит наклоны, а выброшен был
   * именно наклонный вклад — то есть «утечка» была пропорциональна тому,
   * насколько разболтано поле. Ни геометрия (флюид на границе с твёрдой
   * ячейкой — ноль), ни две полусферные суммы (совпадают до 4.4e-16) виноваты
   * не были, и обе гипотезы проверены ЧИСЛОМ, а не рассуждением. */
  st->pabs = 0.0;
  for (int32_t c = 0; c < nc; c++) {
    double sa = p->sig_t[c] - p->sig_s[c];
    if (cu != NULL) {
      for (int i = 0; i < 4; i++)
        st->pabs += sa * phi[c * 4 + i] * cu->mvol[c][0][i];
    } else {
      double vol = (double)m->csize[c] * (double)m->csize[c] * (double)m->csize[c] * m->fr.u[0] *
                   m->fr.u[1] * m->fr.u[2];
      st->pabs += sa * phi[c * 4] * vol; /* у коробки ∫b_i = 0 при i ≥ 1 */
    }
  }

  /* ПОВЕРХНОСТИ ВНУТРИ ОБЛАСТИ — ТОЖЕ СТАТЬЯ БАЛАНСА (К40).
   *
   * Прежде баланс считался как `pin − pout − pabs`, и в сцене БЕЗ разреза это
   * было точное дискретное тождество (измерено 1.5e-14). Но как только внутри
   * появилась поверхность, она стала брать энергию себе, а в тождестве её не
   * было: невязка выходила `1.66e2` при пропускной способности `1.0e4`, и она
   * была ЗАКОННОЙ — сфера с ρ = 0.78 поглощает 22% того, что на неё падает.
   *
   * ЧЕМ ЭТО БЫЛО ПЛОХО: проверка, у которой законная невязка составляет
   * полтора процента, не отличит от неё ошибку схемы в полтора процента. То
   * есть главный инвариант метода (К13: тождество верно на ЛЮБОЙ итерации, а
   * не только после сходимости) в разрезанной сцене НЕ ПРОВЕРЯЛ НИЧЕГО. Это
   * третья подпись артефакта в чистом виде — метрика, нечувствительная к той
   * ошибке, которую она якобы стережёт.
   *
   * Теперь считаются обе стороны отдельно: `psin` — мощность, УШЕДШАЯ из
   * объёма в поверхности, `psout` — мощность, отданная поверхностями обратно
   * (отражение плюс собственное излучение). Тождество:
   *
   *     pin + psout = pout + pabs + psin
   *
   * `psin` берётся из того же накопителя `sinf`, которым считается облучённость,
   * поэтому это не независимая оценка, а ровно та величина, что вошла в схему.
   * `psout` — интеграл хранимого исходящего радианса по элементу, умноженный на
   * СОБСТВЕННУЮ полусферную сумму набора (К29): нормировка из континуума здесь
   * потеряла бы энергию ровно так же, как теряла её в отражении. */
  st->psin = st->psout = 0.0;
  for (int32_t e = 0; e < nse; e++) {
    st->psin += sinf[e * 4];
    double io = 0.0; /* ∫ L_out dA = Σ_i sout_i · ∫b_i dA */
    for (int i = 0; i < 4; i++)
      io += sout[e * 4 + i] * cu->se[e].m[i][0];
    st->psout += hs_out[e] * io;
  }
  /* К91: ОБЪЁМНЫЙ ИСТОЧНИК В ТОЖДЕСТВЕ ОТСУТСТВОВАЛ.
   *
   * Тождество читалось `pin + psout = pout + pabs + psin` и было верно ровно до
   * тех пор, пока `ε = 0`. Ни один тест с объёмным источником баланс не
   * проверял, поэтому дыра держалась с шага C: в задаче с `ε` невязка выходила
   * равной ВСЕЙ излучённой мощности (измерено: `1.14e+01` от втекшего), и это
   * читалось бы как ошибка схемы.
   *
   * Излучённая мощность есть `∫∫ ε dΩ dV`. Изотропная часть даёт `4π·∫ε dV`, а
   * НАПРАВЛЕННАЯ `eps_dir·ω` при интегрировании по сфере даёт ноль — поэтому в
   * тождество она не входит вовсе, и это надо назвать, а не молча опустить.
   * Объёмный интеграл берётся ВСЕМИ ЧЕТЫРЬМЯ моментами, а не одним средним:
   * у РАЗРЕЗАННОЙ ячейки `∫b_i dV ≠ 0` при `i ≥ 1`. Это дословно К48, и то, что
   * та же поправка понадобилась второй раз, — довод считать её правилом. */
  st->pemit = 0.0;
  if (p->eps != NULL)
    for (int32_t c = 0; c < nc; c++) {
      if (cu != NULL) {
        for (int i = 0; i < 4; i++)
          st->pemit += 4.0 * M_PI * p->eps[c * 4 + i] * cu->mvol[c][0][i];
      } else {
        double vol = (double)m->csize[c] * (double)m->csize[c] * (double)m->csize[c] * m->fr.u[0] *
                     m->fr.u[1] * m->fr.u[2];
        st->pemit += 4.0 * M_PI * p->eps[c * 4] * vol;
      }
    }
  st->balance = st->pin + st->psout + st->pemit - st->pout - st->pabs - st->psin - st->psolid;
  st->iters = it;
  st->resid = resid;
  st->nclip = nclip_last;
  st->nfallback = nfb;
  st->bout = bout;
  st->sout = sout;
  /* §744: экспорт средней облучённости элементов финального такта — чистая E
   * без альбедо (sinf0/M00). Финальный такт, не сошедшийся предел (А1152). */
  if (cu != NULL && nse > 0) {
    st->eirr = malloc((size_t)nse * sizeof *st->eirr);
    if (st->eirr != NULL)
      for (int32_t e = 0; e < nse; e++)
        st->eirr[e] = cu->se[e].m[0][0] > 0.0 ? sinf[e * 4] / cu->se[e].m[0][0] : 0.0;
  }

  free(Lall);
  free(phinall);
  free(indegall);
  free(orderall);
  free(queueall);
  free(phin);
  free(binf);
  free(sinf);
  free(kadbg);
  free(gdbg);
  free(phiprev_dbg);
  free(chnc);
  free(chnl);
  free(hs_se);
  free(semax);
  free(bfmax);
  free(hs_out);
  free(bprev);
  free(sprev);
  free(mfid);
  free(mspec);
  free(mspin);
  free(mprev);
  free(binfall);
  free(sinfall);
  return 0;
}
