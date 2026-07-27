/* PLAN_CUT.md Р-4. Последовательное отсечение выпуклого многогранника
 * полупространствами (Сазерленд–Ходжмен в 3D) с сохранением топологии. */

/* Г31: ВОДОНЕПРОНИЦАЕМОСТЬ МОЖЕТ БЫТЬ СЛОМАНА ФЛАГАМИ, А НЕ КОДОМ. Всё, что
 * ниже написано про побитовое совпадение, верно для строгого IEEE и ложно при
 * слиянии умножения-сложения: n·v в двух местах программы дало бы разные
 * последние биты просто из-за разного инлайнинга, появились бы щели, и ни один
 * тест исходника их не объяснил бы.
 *
 * МЕХАНИЗМ — ФЛАГ СБОРКИ, А НЕ PRAGMA. Здесь стояла `#pragma STDC FP_CONTRACT
 * OFF`, как и записано в Р-4; замер показал, что GCC её НЕ РЕАЛИЗУЕТ и молча
 * пропускает (-Wunknown-pragmas), то есть защита была бы декоративной. Поэтому
 * `-ffp-contract=off` стоит в CFLAGS проекта, а `-ffast-math` запрещён; сверка
 * сборок -O0 и -O2 на побитовое совпадение — в `make test`. */

#include "poly3.h"
#include <math.h>
#include <stdint.h>

/* Классификация вершины относительно плоскости: ровно ноль есть третий случай
 * «НА плоскости», а не «внутри с погрешностью». -Wfloat-equal гасится здесь
 * точечно, потому что предупреждение описывает ровно то, что делается
 * сознательно: допуска в классификации нет и быть не должно (Г21/Р-4).
 * memcmp тут негоден — -0.0 обязан считаться нулём, а бит знака у него другой. */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wfloat-equal"
static int on_plane(double s) {
  return s == 0.0;
}
#pragma GCC diagnostic pop

void hz_frame_plane(const hz_frame *fr, const double n[3], double off, hz_hspace *h) {
  /* мир = o + u∘i, значит n·x <= off  <=>  Σ(n_a u_a) i_a <= off - n·o */
  double d = 0.0;
  for (int a = 0; a < 3; a++) {
    h->n[a] = n[a] * fr->u[a];
    d += n[a] * fr->o[a];
  }
  h->off = off - d;
}

/* --- коробка --------------------------------------------------------------- */

/* Петли шести граней по углам куба, занумерованным битами (1=x, 2=y, 4=z).
 * Каждая — ПРОТИВ ЧАСОВОЙ СТРЕЛКИ, ЕСЛИ СМОТРЕТЬ СНАРУЖИ; проверено тем, что
 * (v1-v0) x (v2-v1) смотрит вдоль внешней нормали грани. */
static const int8_t box_face[6][4] = {
    {0, 4, 6, 2}, /* -x */ {1, 3, 7, 5}, /* +x */
    {0, 1, 5, 4}, /* -y */ {2, 6, 7, 3}, /* +y */
    {0, 2, 3, 1}, /* -z */ {4, 5, 7, 6}, /* +z */
};

int hz_poly3_box(hz_poly3 *p, const int32_t lo[3], const int32_t hi[3]) {
  for (int a = 0; a < 3; a++)
    if (hi[a] <= lo[a]) return HZ_P3_EMPTY;
  p->nv = 8;
  for (int c = 0; c < 8; c++)
    for (int a = 0; a < 3; a++)
      p->v[c][a] = (double)((c >> a) & 1 ? hi[a] : lo[a]);
  p->nf = 6;
  for (int f = 0; f < 6; f++) {
    p->fsrc[f] = ~(int32_t)f;
    p->floff[f] = 4 * (int32_t)f;
    for (int k = 0; k < 4; k++)
      p->fl[4 * f + k] = box_face[f][k];
  }
  p->floff[6] = 24;
  return HZ_P3_OK;
}

/* --- отсечение ------------------------------------------------------------- */

/* Лексикографический порядок точек. Только < и >, без сравнений на равенство:
 * -Wfloat-equal ругался бы по делу, а порядок от этого не зависит. */
static int lex_less(const double a[3], const double b[3]) {
  for (int k = 0; k < 3; k++) {
    if (a[k] < b[k]) return 1;
    if (b[k] < a[k]) return 0;
  }
  return 0;
}

/* Пересечение ребра с плоскостью, ВСЕГДА от канонически упорядоченной пары
 * концов. Иначе две стороны посчитали бы t против 1-t и разошлись бы в
 * последнем бите — та самая щель Г10.
 *
 * ПОРЯДОК ПО КООРДИНАТАМ, А НЕ ПО ИНДЕКСАМ (Г39, найдено тестом). Сначала здесь
 * стояло (min, max) ИНДЕКСОВ вершин. Этого хватает внутри одной ячейки, но
 * индексы назначаются в порядке обхода, и у соседа они ДРУГИЕ: одна и та же
 * точка на общей стене получалась из пары (p,q) у одного и (q,p) у другого.
 * На одной плоскости совпадало случайно — нумерация углов коробки давала тот же
 * относительный порядок, — и тест пункта 4 это пропустил; на фасетах сферы
 * разошлось. Лексикографический порядок КООРДИНАТ от ячейки не зависит вовсе,
 * а координаты концов у соседей уже совпадают побитово по индукции от целых
 * углов коробки. */
static void edge_cross(const hz_poly3 *in, const double *s, int32_t a, int32_t b, double out[3]) {
  int32_t p = a, q = b;
  if (!lex_less(in->v[p], in->v[q])) {
    p = b;
    q = a;
  }
  double sp = s[p], sq = s[q];
  double t = sp / (sp - sq);
  for (int k = 0; k < 3; k++)
    out[k] = in->v[p][k] + t * (in->v[q][k] - in->v[p][k]);
}

int hz_poly3_clip(const hz_poly3 *in, const hz_hspace *h, int32_t hid, hz_poly3 *out) {
  double s[HZ_P3_MAXV];
  for (int32_t i = 0; i < in->nv; i++) {
    double d = 0.0;
    for (int a = 0; a < 3; a++)
      d += h->n[a] * in->v[i][a];
    s[i] = d - h->off; /* <0 внутри, >0 снаружи, ==0 НА плоскости */
  }

  int32_t vmap[HZ_P3_MAXV]; /* старая вершина -> новая, -1 = не перенесена */
  for (int32_t i = 0; i < in->nv; i++)
    vmap[i] = -1;
  /* склейка пересечений по КАНОНИЧЕСКОМУ ребру: вершина на ребре, общем для
   * двух граней, создаётся ОДИН раз (линейный поиск, <= HZ_P3_MAXV, детерминирован) */
  int32_t xa[HZ_P3_MAXV], xb[HZ_P3_MAXV], xi[HZ_P3_MAXV];
  int32_t nx = 0;
  uint8_t on[HZ_P3_MAXV]; /* новая вершина лежит на плоскости */

  out->nv = 0;
  out->nf = 0;
  out->floff[0] = 0;

  int32_t cap_a[HZ_P3_MAXF], cap_b[HZ_P3_MAXF];
  int32_t ncap = 0;
  int coplanar = 0; /* нашлась грань, целиком лежащая на плоскости */

  for (int32_t f = 0; f < in->nf; f++) {
    int32_t b0 = in->floff[f], b1 = in->floff[f + 1], k = b1 - b0;
    int32_t loop[HZ_P3_MAXFL];
    int32_t nl = 0;
    for (int32_t e = 0; e < k; e++) {
      int32_t a = in->fl[b0 + e], b = in->fl[b0 + (e + 1 == k ? 0 : e + 1)];
      if (s[a] <= 0.0) { /* внутри или НА плоскости — переносится */
        if (vmap[a] < 0) {
          if (out->nv >= HZ_P3_MAXV) return HZ_P3_ECAPACITY;
          vmap[a] = out->nv;
          for (int c = 0; c < 3; c++)
            out->v[out->nv][c] = in->v[a][c];
          on[out->nv] = on_plane(s[a]) ? 1u : 0u;
          out->nv++;
        }
        if (nl >= HZ_P3_MAXFL) return HZ_P3_ECAPACITY;
        loop[nl++] = vmap[a];
      }
      /* пересечение — только при СТРОГОЙ смене знака; вершина на плоскости уже
       * перенесена выше, и лишней точки не создаётся (Г8 №4: куска нулевого
       * объёма не рождается) */
      if ((s[a] < 0.0 && s[b] > 0.0) || (s[a] > 0.0 && s[b] < 0.0)) {
        int32_t p = a < b ? a : b, q = a < b ? b : a, found = -1;
        for (int32_t t = 0; t < nx; t++)
          if (xa[t] == p && xb[t] == q) {
            found = xi[t];
            break;
          }
        if (found < 0) {
          if (out->nv >= HZ_P3_MAXV || nx >= HZ_P3_MAXV) return HZ_P3_ECAPACITY;
          edge_cross(in, s, a, b, out->v[out->nv]);
          on[out->nv] = 1u;
          found = out->nv++;
          xa[nx] = p;
          xb[nx] = q;
          xi[nx] = found;
          nx++;
        }
        if (nl >= HZ_P3_MAXFL) return HZ_P3_ECAPACITY;
        loop[nl++] = found;
      }
    }

    /* убрать подряд идущие повторы (замкнуто по кругу) */
    int32_t m = 0;
    for (int32_t e = 0; e < nl; e++)
      if (loop[e] != loop[(e + nl - 1) % nl]) loop[m++] = loop[e];
    if (m < 3) continue; /* выродилась в точку или отрезок — грани нет */

    if (out->nf >= HZ_P3_MAXF) return HZ_P3_ECAPACITY;
    int32_t base = out->floff[out->nf];
    if (base + m > HZ_P3_MAXFL) return HZ_P3_ECAPACITY;
    for (int32_t e = 0; e < m; e++)
      out->fl[base + e] = loop[e];
    out->fsrc[out->nf] = in->fsrc[f];
    out->nf++;
    out->floff[out->nf] = base + m;

    /* Рёбра будущей крышки. Грань отдаёт ребро (y,x) для своей пары соседних
     * вершин, обе из которых лежат на плоскости: смежные грани замкнутого тела
     * проходят общее ребро в ПРОТИВОПОЛОЖНЫХ направлениях, поэтому крышка
     * берёт его развёрнутым. */
    int32_t allon = 1;
    for (int32_t e = 0; e < m; e++)
      if (!on[loop[e]]) {
        allon = 0;
        break;
      }
    if (allon) {
      coplanar = 1; /* эта грань И ЕСТЬ крышка, второй заводить нельзя */
      continue;
    }
    for (int32_t e = 0; e < m; e++) {
      int32_t x = loop[e], y = loop[(e + 1) % m];
      if (on[x] && on[y]) {
        if (ncap >= HZ_P3_MAXF) return HZ_P3_ECAPACITY;
        cap_a[ncap] = y;
        cap_b[ncap] = x;
        ncap++;
      }
    }
  }

  if (out->nf == 0) return HZ_P3_EMPTY;

  if (!coplanar && ncap >= 3) {
    /* сцепить рёбра в петлю: следующее начинается там, где кончилось прошлое */
    int32_t used[HZ_P3_MAXF];
    for (int32_t t = 0; t < ncap; t++)
      used[t] = 0;
    int32_t base = out->floff[out->nf];
    if (out->nf >= HZ_P3_MAXF || base + ncap > HZ_P3_MAXFL) return HZ_P3_ECAPACITY;
    int32_t start = cap_a[0], cur = cap_b[0], nl = 0;
    used[0] = 1;
    out->fl[base + nl++] = cap_a[0];
    while (cur != start) {
      int32_t nxt = -1;
      for (int32_t t = 0; t < ncap; t++)
        if (!used[t] && cap_a[t] == cur) {
          nxt = t;
          break;
        }
      if (nxt < 0) return HZ_P3_ECAP; /* не замкнулась: вырождение, молчать нельзя */
      used[nxt] = 1;
      if (base + nl >= HZ_P3_MAXFL) return HZ_P3_ECAPACITY;
      out->fl[base + nl++] = cap_a[nxt];
      cur = cap_b[nxt];
    }
    if (nl >= 3) {
      out->fsrc[out->nf] = hid;
      out->nf++;
      out->floff[out->nf] = base + nl;
    }
  }
  return out->nf >= 4 ? HZ_P3_OK : HZ_P3_EMPTY;
}

int hz_poly3_cut(hz_poly3 *out, const int32_t lo[3], const int32_t hi[3], const hz_hspace *h,
                 const int32_t *hid, int nh) {
  /* Пинг-понг между out и tmp: отсечение читает один буфер и пишет другой,
   * поэтому вершины сжимаются на каждом шаге и мёртвых не остаётся.
   * НАЧАЛЬНЫЙ БУФЕР ВЫБИРАЕТСЯ ПО ЧЁТНОСТИ nh, чтобы последний результат осел
   * ровно в out. Так было не сразу: сначала коробка всегда клалась в out, а в
   * конце при нечётном nh делалось `*out = tmp` — копия ЦЕЛОЙ структуры, то
   * есть ~8 КБ, из которых прописана только начальная часть массивов. Хвосты
   * при этом копируются неинициализированными (cppcheck это и поймал), а
   * заодно копия стоит на ровном месте. Выбор буфера убирает и то и другое. */
  hz_poly3 tmp;
  hz_poly3 *buf[2];
  int s = nh & 1;
  buf[0] = out;
  buf[1] = &tmp;
  int rc = hz_poly3_box(buf[s], lo, hi);
  if (rc != HZ_P3_OK) return rc;
  for (int j = 0; j < nh; j++) {
    rc = hz_poly3_clip(buf[s], &h[j], hid[j], buf[1 - s]);
    s = 1 - s;
    if (rc != HZ_P3_OK) {
      out->nf = 0; /* ничего не осталось (или вырождение) — не отдавать полуфабрикат */
      return rc;
    }
  }
  return HZ_P3_OK;
}

int hz_poly3_complement(hz_poly3 *out, int maxout, int *nout, const int32_t lo[3],
                        const int32_t hi[3], const hz_hspace *h, const int32_t *hid,
                        const int32_t *hid_flipped, int nh) {
  hz_hspace hs[HZ_P3_MAXH];
  int32_t ids[HZ_P3_MAXH];
  *nout = 0;
  if (nh > HZ_P3_MAXH) return HZ_P3_ECAPACITY;
  for (int j = 0; j < nh; j++) {
    if (*nout >= maxout) return HZ_P3_ECAPACITY;
    for (int t = 0; t < j; t++) {
      hs[t] = h[t];
      ids[t] = hid[t];
    }
    for (int a = 0; a < 3; a++)
      hs[j].n[a] = -h[j].n[a]; /* отрицание точно: биты пересечений не поедут */
    hs[j].off = -h[j].off;
    ids[j] = hid_flipped[j];
    int rc = hz_poly3_cut(&out[*nout], lo, hi, hs, ids, j + 1);
    if (rc == HZ_P3_OK)
      (*nout)++;
    else if (rc != HZ_P3_EMPTY)
      return rc;
  }
  return HZ_P3_OK;
}

/* --- метрика --------------------------------------------------------------- */

double hz_poly3_volume(const hz_poly3 *p, const hz_frame *fr) {
  /* Теорема о дивергенции: сумма знаковых тетраэдров с вершиной в начале
   * координат. Для замкнутой поверхности с петлями против часовой снаружи она
   * не зависит от выбора начала. */
  double v6 = 0.0;
  for (int32_t f = 0; f < p->nf; f++) {
    int32_t b0 = p->floff[f], k = p->floff[f + 1] - b0;
    const double *a = p->v[p->fl[b0]];
    for (int32_t e = 1; e + 1 < k; e++) {
      const double *b = p->v[p->fl[b0 + e]], *c = p->v[p->fl[b0 + e + 1]];
      v6 += a[0] * (b[1] * c[2] - b[2] * c[1]) + a[1] * (b[2] * c[0] - b[0] * c[2]) +
            a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
  }
  /* объём — единственная величина, которая переносится в мир умножением */
  return v6 / 6.0 * fr->u[0] * fr->u[1] * fr->u[2];
}

/* Σ (w_i x w_{i+1}) по петле в МИРОВЫХ координатах: длина = 2*площадь,
 * направление = внешняя нормаль. В единицах это не площадь и не нормаль. */
static void face_vec(const hz_poly3 *p, int f, const hz_frame *fr, double acc[3]) {
  int32_t b0 = p->floff[f], k = p->floff[f + 1] - b0;
  acc[0] = acc[1] = acc[2] = 0.0;
  for (int32_t e = 0; e < k; e++) {
    const double *vi = p->v[p->fl[b0 + e]], *vj = p->v[p->fl[b0 + (e + 1 == k ? 0 : e + 1)]];
    double a[3], b[3];
    for (int c = 0; c < 3; c++) {
      a[c] = fr->o[c] + fr->u[c] * vi[c];
      b[c] = fr->o[c] + fr->u[c] * vj[c];
    }
    acc[0] += a[1] * b[2] - a[2] * b[1];
    acc[1] += a[2] * b[0] - a[0] * b[2];
    acc[2] += a[0] * b[1] - a[1] * b[0];
  }
}

double hz_poly3_area(const hz_poly3 *p, int f, const hz_frame *fr) {
  double a[3];
  face_vec(p, f, fr, a);
  return 0.5 * sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

void hz_poly3_normal(const hz_poly3 *p, int f, const hz_frame *fr, double nout[3]) {
  double a[3];
  face_vec(p, f, fr, a);
  double m = sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
  for (int c = 0; c < 3; c++)
    nout[c] = m > 0.0 ? a[c] / m : 0.0;
}
