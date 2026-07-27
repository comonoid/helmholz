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
  double *hs_se = calloc((size_t)(nse > 0 ? nse : 1), sizeof(double));
  binfall = calloc((size_t)nth * (size_t)m->nf * 4, sizeof(double));
  sinfall = calloc((size_t)nth * (size_t)(nse > 0 ? nse : 1) * 4, sizeof(double));
  if (Lall == NULL || phinall == NULL || indegall == NULL || orderall == NULL || queueall == NULL ||
      phin == NULL || binfall == NULL || sinfall == NULL || bout == NULL || binf == NULL ||
      sout == NULL || sinf == NULL || hs_se == NULL) {
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
    free(hs_se);
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
    double s = 0.0;
    for (int mm = 0; mm < nd; mm++) {
      double on =
          d->ox[mm] * cu->se[e].n[0] + d->oy[mm] * cu->se[e].n[1] + d->oz[mm] * cu->se[e].n[2];
      if (on < 0.0) s += d->w[mm] * (-on); /* приходящие НА поверхность */
    }
    hs_se[e] = s;
  }

  memset(phi, 0, (size_t)nc * 4 * sizeof(double));
  if (p->wall_rho != NULL)
    for (int32_t f = 0; f < m->nf; f++)
      if (m->f[f].cb < 0)
        bout[f * 4] = p->wall_emit != NULL ? p->wall_emit[(int)(~m->f[f].cb)] : 0.0;
  for (int32_t e = 0; e < nse; e++)
    if (p->facet_emit != NULL && cu->se[e].facet < p->nfacet)
      sout[e * 4] = p->facet_emit[cu->se[e].facet];

  int nclip_last = 0, it = 0, nfb = 0;
  double resid = 0.0;
  for (it = 0; it < maxit; it++) {
    memset(phin, 0, (size_t)nc * 4 * sizeof(double));
    memset(binf, 0, (size_t)m->nf * 4 * sizeof(double));
    memset(sinf, 0, (size_t)(nse > 0 ? nse : 1) * 4 * sizeof(double));
    memset(phinall, 0, (size_t)nth * (size_t)nc * 4 * sizeof(double));
    memset(binfall, 0, (size_t)nth * (size_t)m->nf * 4 * sizeof(double));
    memset(sinfall, 0, (size_t)nth * (size_t)(nse > 0 ? nse : 1) * 4 * sizeof(double));
    nclip_last = 0;
    nfb = 0;
    st->pin = st->pout = st->pabs = 0.0;
    double pin_acc = 0.0, pout_acc = 0.0;
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
              if (p->wall_rho != NULL) {
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
            } else {
              int32_t up = mine_is_a ? m->f[f].cb : m->f[f].ca;
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
            }
          }

        double a0row[4], rhs0 = rhs[0];
        for (int j = 0; j < 4; j++)
          a0row[j] = A[0][j];
        double cf[4];
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
        for (int j = 0; j < 4; j++)
          L[c * 4 + j] = cf[j];

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
          }
      }

      for (int32_t c = 0; c < nc; c++)
        for (int j = 0; j < 4; j++)
          phit[c * 4 + j] += d->w[mm] * L[c * 4 + j];

      for (int32_t f = 0; f < m->nf; f++) {
        if (m->f[f].cb >= 0) continue;
        int wall = (int)(~m->f[f].cb);
        double on = om[m->f[f].axis] * ((wall & 1) ? 1.0 : -1.0);
        if (!(on > 0.0)) continue;
        int32_t c = m->f[f].ca;
        double fmmb[4][4];
        const double (*fmm)[4];
        if (cu != NULL) {
          fmm = cu->ffm[f];
        } else {
          double v[4][3];
          tr3_face_corners(m, f, v);
          face_mass2(m, (const double (*)[4][3]) & v, c, c, fmmb);
          fmm = fmmb;
        }
        double accj[4] = {0, 0, 0, 0};
        for (int j = 0; j < 4; j++)
          for (int i = 0; i < 4; i++)
            accj[j] += L[c * 4 + i] * fmm[i][j];
        pout_acc += on * d->w[mm] * accj[0];
        for (int j = 0; j < 4; j++)
          binft[f * 4 + j] += on * d->w[mm] * accj[j];
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
      free(hs_se);
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

    /* --- стенки и поверхности: новый исходящий радианс, DG1 по положению --- */
    if (p->wall_rho != NULL)
      for (int32_t f = 0; f < m->nf; f++) {
        if (m->f[f].cb >= 0) continue;
        int wall = (int)(~m->f[f].cb);
        double fmm[4][4], rr[4], ee[4];
        if (cu != NULL) {
          memcpy(fmm, cu->ffm[f], sizeof fmm);
        } else {
          double v[4][3];
          tr3_face_corners(m, f, v);
          face_mass2(m, (const double (*)[4][3]) & v, m->f[f].ca, m->f[f].ca, fmm);
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
        tr3_face_null(m, f, m->f[f].ca, nul);
        if (tr3_project_plane(fmm, rr, nul, ee) != 0 || corner_min(ee) < 0.0) {
          nfb++;
          /* ОТКАЗ В ЗАКРЫТУЮ СТОРОНУ: либо элемент выродился геометрически, либо
           * проекция ушла в минус. Тогда остаётся ТОЧНОЕ среднее по грани —
           * константа, представимая в любом случае. Порога в критерии нет. */
          ee[0] = binf[f * 4] / fmm[0][0];
          ee[1] = ee[2] = ee[3] = 0.0;
        }
        for (int j = 0; j < 4; j++)
          bout[f * 4 + j] = (hsum[wall] > 0.0 ? p->wall_rho[wall] * ee[j] / hsum[wall] : 0.0);
        if (p->wall_emit != NULL) bout[f * 4] += p->wall_emit[wall];
      }
    for (int32_t e = 0; e < nse; e++) {
      const tr3_selem *se = &cu->se[e];
      double rho = (p->facet_rho != NULL && se->facet < p->nfacet) ? p->facet_rho[se->facet] : 0.0;
      double em = (p->facet_emit != NULL && se->facet < p->nfacet) ? p->facet_emit[se->facet] : 0.0;
      double fmm[4][4], rr[4], ee[4];
      memcpy(fmm, se->m, sizeof fmm);
      for (int j = 0; j < 4; j++)
        rr[j] = sinf[e * 4 + j];
      if (!(fmm[0][0] > 0.0)) continue;
      /* та же правка К39: нулевой вектор у поверхностного элемента хранится
       * рядом с матрицей масс и выведен из уравнения плоскости (cut3.h) */
      if (tr3_project_plane(fmm, rr, se->nul, ee) != 0 || corner_min(ee) < 0.0) {
        ee[0] = sinf[e * 4] / fmm[0][0];
        ee[1] = ee[2] = ee[3] = 0.0;
        nfb++;
      }
      for (int j = 0; j < 4; j++)
        sout[e * 4 + j] = (hs_se[e] > 0.0 ? rho * ee[j] / hs_se[e] : 0.0);
      sout[e * 4] += em;
    }

    resid = 0.0;
    for (int32_t i = 0; i < nc * 4; i++) {
      double dd = fabs(phin[i] - phi[i]);
      if (dd > resid) resid = dd;
      phi[i] = phin[i];
    }
    if (resid < tol) break;
  }

  st->pabs = 0.0;
  for (int32_t c = 0; c < nc; c++) {
    double vol = cu != NULL ? cu->mvol[c][0][0]
                            : (double)m->csize[c] * (double)m->csize[c] * (double)m->csize[c] *
                                  m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
    st->pabs += (p->sig_t[c] - p->sig_s[c]) * phi[c * 4] * vol;
  }
  st->balance = st->pin - st->pout - st->pabs;
  st->iters = it;
  st->resid = resid;
  st->nclip = nclip_last;
  st->nfallback = nfb;
  st->bout = bout;
  st->sout = sout;

  free(Lall);
  free(phinall);
  free(indegall);
  free(orderall);
  free(queueall);
  free(phin);
  free(binf);
  free(sinf);
  free(hs_se);
  free(binfall);
  free(sinfall);
  return 0;
}
