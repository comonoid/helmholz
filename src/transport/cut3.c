/* PLAN_TRANSPORT.md T4. Стык развёртки с геометрией разреза: флюидная область
 * ячейки и поверхностные элементы на фасетах. */

#include "transport/cut3.h"
#include "transport/tet3.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

/* метки для hz_poly3_complement: неперевёрнутые плоскости помечаются как
 * ВНУТРЕННИЕ (по ту сторону — другой кусок той же флюидной области), а
 * перевёрнутая — как поверхность. Диапазоны не пересекаются с ~(0..5) у коробки. */
#define CUT3_INNER (-1000)

/* §772: уникальные ссылки плоскостей ПОДДЕРЕВА узла. Обе ориентации одной
 * плоскости ЛЕГАЛЬНЫ (материал-слэб тонкой стены выпукл); невыпуклый материал
 * ловится не здесь, а объёмной сверкой в предикате огрубления. Возврат: число
 * ссылок; −1 — не влезло в max (кандидат не грубится). */
int tr3_cut_subtree_refs(const hz_octree *t, const hz_cutmap *cm, int32_t ni, int32_t *refs,
                         int max) {
  /* стек: 7·глубина + 1 ≤ 7·16 + 1 = 113 при пределе решётки uint16 */
  int32_t stack[128];
  int sp = 0, n = 0;
  stack[sp++] = ni;
  while (sp > 0) {
    int32_t cur = stack[--sp];
    const hz_cutrec *r = hz_cutmap_find(cm, cur);
    if (r != NULL)
      for (int32_t j = 0; j < r->nf; j++) {
        int32_t ref = cm->fref[r->f0 + j];
        int dup = 0;
        for (int q = 0; q < n; q++)
          if (refs[q] == ref) {
            dup = 1;
            break;
          }
        if (dup) continue;
        if (n >= max) return -1;
        refs[n++] = ref;
      }
    if (t->nodes[cur].child0 >= 0) {
      if (sp + 8 > (int)(sizeof stack / sizeof stack[0])) return -1;
      for (int k = 0; k < 8; k++)
        stack[sp++] = t->nodes[cur].child0 + k;
    }
  }
  return n;
}

/* §772: флюидный объём коробки узла, резанной набором ссылок — тем же
 * дополнением, что рабочий путь сборки. Пустой набор — объём коробки тем же
 * ядром (одна арифметика на обе стороны сверки). Возврат < 0 — отказ ядра. */
double tr3_cut_refs_fluid_vol(const hz_frame *fr, const hz_facettab *ft, const int32_t *refs,
                              int nrefs, const int32_t lo[3], int32_t size) {
  int32_t hi[3] = {lo[0] + size, lo[1] + size, lo[2] + size};
  if (nrefs > HZ_P3_MAXH || nrefs < 0) return -1.0;
  if (nrefs == 0) {
    hz_poly3 box;
    if (hz_poly3_cut(&box, lo, hi, NULL, NULL, 0) != HZ_P3_OK) return -1.0;
    return hz_poly3_volume(&box, fr);
  }
  hz_hspace h[HZ_P3_MAXH];
  int32_t hid[HZ_P3_MAXH], hflip[HZ_P3_MAXH];
  for (int j = 0; j < nrefs; j++) {
    int32_t ref = refs[j], fi = ref >= 0 ? ref : ~ref;
    if (fi < 0 || fi >= ft->n) return -1.0;
    const hz_facet *f = &ft->f[fi];
    for (int a = 0; a < 3; a++)
      h[j].n[a] = ref >= 0 ? f->n[a] : -f->n[a];
    h[j].off = ref >= 0 ? f->off : -f->off;
    hid[j] = CUT3_INNER;
    hflip[j] = CUT3_INNER;
  }
  hz_poly3 *pieces = calloc((size_t)nrefs, sizeof(hz_poly3));
  if (pieces == NULL) return -1.0;
  int npc = 0;
  if (hz_poly3_complement(pieces, nrefs, &npc, lo, hi, h, hid, hflip, nrefs) != HZ_P3_OK) {
    free(pieces);
    return -1.0;
  }
  double v = 0.0;
  for (int p = 0; p < npc; p++)
    v += hz_poly3_volume(&pieces[p], fr);
  free(pieces);
  return v;
}
#define CUT3_SURF0 (1000)

/* Р-8 (Д3): клип многоугольника элемента РЕБЁРНЫМИ полуплоскостями куска.
 * Кусок — треугольник фасета (единицы кадра), клип идёт В МИРЕ; сторона
 * каждой ребёрной полуплоскости берётся по третьей вершине (знак, не допуск —
 * порога здесь нет). Точки ровно на ребре держат обе смежные полуплоскости:
 * перекрытие меры нуль, а щель не возникает (Г61: биты не обещаются, площади
 * складываются в кусок с плавучей точностью). Ёмкость CUT3_CLIPV: вход ≤ 32,
 * каждая из трёх полуплоскостей добавляет ≤ 1 вершину. */
#define CUT3_CLIPV 16
/* Р-8 (Д3, вторая редакция): элемент bounded-фасета = КУСОК ∩ КОРОБКА, прямым
 * клипом в МИРЕ. Первая редакция клиповала куском грань «материального
 * многогранника» — у тонкой поверхности материал вырожден, его грани стоят не
 * там (класс «избыток» §779), и пересечение с куском исчезало (замерено:
 * дефицит 98.7 % на cavity). Чужие полуплоскости элементу с куском не нужны:
 * кусок сам несёт свои границы, а покрытие «кусок ∩ коробка» точно по
 * построению (куски поверхности не перекрываются, коробки — разбиение).
 * Ёмкость: треугольник × 6 полуплоскостей ≤ 3 + 6 вершин. */
static int cut3_piece_in_box(const tr3_mesh *m, const hz_facet *fp, const int32_t lo[3],
                             const int32_t hi[3], double (*vw)[3]) {
  double buf[2][CUT3_CLIPV][3];
  for (int i = 0; i < 3; i++)
    for (int a = 0; a < 3; a++)
      buf[0][i][a] = m->fr.o[a] + m->fr.u[a] * fp->tv[i][a];
  int cur = 0, ncp = 3;
  for (int ax = 0; ax < 3 && ncp >= 3; ax++)
    for (int side = 0; side < 2 && ncp >= 3; side++) {
      double lim = m->fr.o[ax] + m->fr.u[ax] * (double)(side ? hi[ax] : lo[ax]);
      double sgn = side ? -1.0 : 1.0; /* внутри: sgn·(x − lim) ≥ 0 */
      int oth = 1 - cur, nxt = 0;
      for (int e = 0; e < ncp; e++) {
        const double *P = buf[cur][e], *Q = buf[cur][(e + 1) % ncp];
        double dp = sgn * (P[ax] - lim), dq = sgn * (Q[ax] - lim);
        if (dp >= 0.0 && nxt < CUT3_CLIPV) {
          /* поэлементно, не memcpy: cur != oth, но анализатор связи не видит */
          for (int a = 0; a < 3; a++)
            buf[oth][nxt][a] = P[a];
          nxt++;
        }
        if ((dp > 0.0 && dq < 0.0) || (dp < 0.0 && dq > 0.0)) {
          double t = dp / (dp - dq);
          if (nxt < CUT3_CLIPV) {
            for (int a = 0; a < 3; a++)
              buf[oth][nxt][a] = P[a] + t * (Q[a] - P[a]);
            nxt++;
          }
        }
      }
      cur = oth;
      ncp = nxt;
    }
  if (ncp < 3) return 0;
  memcpy(vw, buf[cur], (size_t)ncp * 3 * sizeof(double));
  return ncp;
}

void tr3_cut_free(tr3_cut *cu) {
  free(cu->mvol);
  free(cu->ffm);
  free(cu->ffmb);
  free(cu->ffmx);
  free(cu->farea);
  free(cu->solid);
  free(cu->se);
  free(cu->sestart);
  free(cu->selist);
  memset(cu, 0, sizeof *cu);
}

/* базис ячейки в мировой точке */
static void cbasis(const tr3_mesh *m, int32_t c, const double xw[3], double b[4]) {
  b[0] = 1.0;
  double s = (double)m->csize[c];
  for (int k = 0; k < 3; k++) {
    double xu = (xw[k] - m->fr.o[k]) / m->fr.u[k];
    b[k + 1] = (xu - ((double)m->clo[c][k] + 0.5 * s)) / s;
  }
}

/* ∫ b^A_i b^B_j dA по МИРОВОМУ многоугольнику, плюс площадь */
static double poly_mass2(const tr3_mesh *m, const double (*v)[3], int nv, int32_t ca, int32_t cb,
                         double out[4][4]) {
  memset(out, 0, 16 * sizeof(double));
  double atot = 0.0;
  for (int e = 1; e + 1 < nv; e++) {
    double tri[3][3];
    for (int k = 0; k < 3; k++) {
      tri[0][k] = v[0][k];
      tri[1][k] = v[e][k];
      tri[2][k] = v[e + 1][k];
    }
    double area = tr3_poly_face_area(tri, 3);
    atot += area;
    double ba[3][4], bb[3][4];
    for (int i = 0; i < 3; i++) {
      cbasis(m, ca, tri[i], ba[i]);
      if (cb >= 0)
        cbasis(m, cb, tri[i], bb[i]);
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
  return atot;
}

static int push_se(tr3_cut *cu, const tr3_selem *e) {
  if (cu->nse >= cu->secap) {
    int32_t nc = cu->secap ? cu->secap * 2 : 64;
    tr3_selem *ns = realloc(cu->se, (size_t)nc * sizeof(tr3_selem));
    if (ns == NULL) return 1;
    cu->se = ns;
    cu->secap = nc;
  }
  cu->se[cu->nse++] = *e;
  return 0;
}

/* §782: mvol ГРУБОЙ ячейки — ТОЧНОЙ АГРЕГАЦИЕЙ по листьям поддерева: моменты
 * аддитивны, интегралы листовых кусков берутся сразу в базисе грубой ячейки
 * (центр c3 и размер hh — параметры tr3_mass_matrix, пересчёта матриц нет).
 * Лист без записи — полный бокс либо ноль по листовой сплошности. Вырожденный
 * листовой рез (union съел лист) у кусочной записи без сплошности — полный
 * бокс, ровно правило degen8 листовой сетки (§781). Внутренний узел с записью
 * нарушает контракт Р-8 (записи — листьям дна) и считается в nbad. */
static int aggr_mvol_rec(tr3_cut *cu, const tr3_mesh *m, const hz_facettab *ft, const hz_cutmap *cm,
                         tr3_leaf_solid_fn leaf_solid, void *lsctx, int32_t ni, const int32_t lo[3],
                         int32_t size, const double c3[3], double hh, int32_t c) {
  if (m->tree->nodes[ni].child0 >= 0) {
    if (hz_cutmap_find(cm, ni) != NULL) {
      cu->nbad++;
      return 0;
    }
    int32_t half = size / 2;
    for (int i = 0; i < 8; i++) {
      int32_t clo[3] = {lo[0] + ((i & 1) ? half : 0), lo[1] + ((i & 2) ? half : 0),
                        lo[2] + ((i & 4) ? half : 0)};
      if (aggr_mvol_rec(cu, m, ft, cm, leaf_solid, lsctx, m->tree->nodes[ni].child0 + i, clo, half,
                        c3, hh, c))
        return 1;
    }
    return 0;
  }
  int32_t hi2[3] = {lo[0] + size, lo[1] + size, lo[2] + size};
  const hz_cutrec *r = hz_cutmap_find(cm, ni);
  double mm[4][4];
  if (r == NULL) {
    if (leaf_solid(lsctx, lo, size)) return 0;
    hz_poly3 box;
    if (hz_poly3_cut(&box, lo, hi2, NULL, NULL, 0) != HZ_P3_OK) return 1;
    if (tr3_mass_matrix(&box, &m->fr, c3, hh, mm) == 0)
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          cu->mvol[c][i][j] += mm[i][j];
    return 0;
  }
  hz_hspace h2[HZ_P3_MAXH];
  int32_t hid2[HZ_P3_MAXH], hfl2[HZ_P3_MAXH];
  int nh2 = hz_cutmap_hspaces(ft, cm, r, h2, hid2, hfl2, HZ_P3_MAXH);
  if (nh2 <= 0) return 0; /* как в листовом пути: ячейка остаётся полной */
  for (int j = 0; j < nh2; j++) {
    hid2[j] = CUT3_INNER;
    hfl2[j] = CUT3_INNER;
  }
  hz_poly3 *pcs = calloc((size_t)nh2, sizeof(hz_poly3));
  if (pcs == NULL) return 1;
  int np2 = 0;
  if (hz_poly3_complement(pcs, nh2, &np2, lo, hi2, h2, hid2, hfl2, nh2) != HZ_P3_OK) np2 = 0;
  if (np2 == 0) {
    int allb = 1;
    for (int32_t j = 0; j < r->nf && allb; j++) {
      int32_t rf = cm->fref[r->f0 + j];
      int32_t fi2 = rf >= 0 ? rf : ~rf;
      if (!(fi2 >= 0 && fi2 < ft->n && ft->f[fi2].bounded)) allb = 0;
    }
    if (allb && !leaf_solid(lsctx, lo, size)) {
      hz_poly3 box;
      if (hz_poly3_cut(&box, lo, hi2, NULL, NULL, 0) == HZ_P3_OK &&
          tr3_mass_matrix(&box, &m->fr, c3, hh, mm) == 0)
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
            cu->mvol[c][i][j] += mm[i][j];
    }
    free(pcs);
    return 0;
  }
  for (int p = 0; p < np2; p++)
    if (tr3_mass_matrix(&pcs[p], &m->fr, c3, hh, mm) == 0)
      for (int i = 0; i < 4; i++)
        for (int j = 0; j < 4; j++)
          cu->mvol[c][i][j] += mm[i][j];
  free(pcs);
  return 0;
}

int tr3_cut_build(tr3_cut *cu, const tr3_mesh *m, const hz_facettab *ft, const hz_cutmap *cm,
                  const uint8_t *solid_in) {
  return tr3_cut_build2(cu, m, ft, cm, solid_in, NULL, NULL);
}

int tr3_cut_build2(tr3_cut *cu, const tr3_mesh *m, const hz_facettab *ft, const hz_cutmap *cm,
                   const uint8_t *solid_in, tr3_leaf_solid_fn leaf_solid, void *lsctx) {
  memset(cu, 0, sizeof *cu);
  cu->m = m;
  cu->mvol = calloc((size_t)m->ncell, sizeof(double[4][4]));
  cu->ffm = calloc((size_t)m->nf, sizeof(double[4][4]));
  cu->ffmb = calloc((size_t)m->nf, sizeof(double[4][4]));
  cu->ffmx = calloc((size_t)m->nf, sizeof(double[4][4]));
  cu->farea = calloc((size_t)m->nf, sizeof(double));
  cu->solid = calloc((size_t)m->ncell, 1);
  if (cu->mvol == NULL || cu->ffm == NULL || cu->ffmb == NULL || cu->ffmx == NULL ||
      cu->farea == NULL || cu->solid == NULL) {
    tr3_cut_free(cu);
    return 1;
  }

  /* --- сперва всё как у коробок: полная ячейка и полные грани --- */
  for (int32_t c = 0; c < m->ncell; c++) {
    double s = (double)m->csize[c];
    double vol = s * s * s * m->fr.u[0] * m->fr.u[1] * m->fr.u[2];
    cu->mvol[c][0][0] = vol;
    for (int k = 1; k < 4; k++)
      cu->mvol[c][k][k] = vol / 12.0;
  }
  for (int32_t f = 0; f < m->nf; f++) {
    double v[4][3];
    tr3_face_corners(m, f, v);
    cu->farea[f] = poly_mass2(m, v, 4, m->f[f].ca, m->f[f].ca, cu->ffm[f]);
    poly_mass2(m, v, 4, m->f[f].ca, m->f[f].cb, cu->ffmx[f]);
    if (m->f[f].cb >= 0) poly_mass2(m, v, 4, m->f[f].cb, m->f[f].cb, cu->ffmb[f]);
  }

  /* ЯЧЕЙКИ ЦЕЛИКОМ В МАТЕРИАЛЕ: у них ни граней разреза, ни флюида. */
  if (solid_in != NULL)
    for (int32_t c = 0; c < m->ncell; c++) {
      if (!solid_in[c]) continue;
      cu->solid[c] = 1;
      memset(cu->mvol[c], 0, 16 * sizeof(double));
    }

  if (ft == NULL || cm == NULL) return 0;

  /* Грани полностью твёрдых ячеек обнуляются ПОСЛЕ разбора разрезанных, иначе
   * порядок обхода решал бы, чей вклад уцелеет. */
  /* --- ячейки с границей: флюидная часть и поверхностные элементы --- */
  for (int32_t c = 0; c < m->ncell; c++) {
    const hz_cutrec *rec = hz_cutmap_find(cm, m->node[c]);
    int32_t refs772[HZ_P3_MAXH]; /* §772: ссылки грубой ячейки — элементам нужен фасет */
    hz_hspace h[HZ_P3_MAXH];
    int32_t hid[HZ_P3_MAXH], hflip[HZ_P3_MAXH];
    int nh = 0, aggr8 = 0, over8 = 0;
    if (rec != NULL) {
      nh = hz_cutmap_hspaces(ft, cm, rec, h, hid, hflip, HZ_P3_MAXH);
      if (nh < 0) {
        /* §794: веер записи больше предела ЯДРА (листва OBJ-кусков) — рез
         * невозможен; «первые 96» не режутся (неверный рез хуже полного
         * флюида), элементы строятся ниже ПО ЗАПИСИ без предела. Ячейка
         * остаётся полной коробкой; сплошную решает solid_in. Считается. */
        over8 = 1;
        nh = 0;
        cu->nrezover++;
        if (solid_in != NULL && solid_in[c]) continue; /* сплошная остаётся */
      }
    } else if (m->tree->nodes[m->node[c]].child0 >= 0) {
      /* §772/§782: ГРУБАЯ ячейка (внутренний узел без своей записи) — куски
       * собираются с ПОДДЕРЕВА; при заданном leaf_solid объём считается
       * АГРЕГАЦИЕЙ по листьям, union-рез большой коробки не гоняется. */
      int nr = tr3_cut_subtree_refs(m->tree, cm, m->node[c], refs772, HZ_P3_MAXH);
      for (int j = 0; j < (nr > 0 ? nr : 0); j++) {
        int32_t ref = refs772[j], fi2 = ref >= 0 ? ref : ~ref;
        if (fi2 < 0 || fi2 >= ft->n) {
          nr = 0;
          break;
        }
        const hz_facet *fp = &ft->f[fi2];
        for (int a = 0; a < 3; a++)
          h[j].n[a] = ref >= 0 ? fp->n[a] : -fp->n[a];
        h[j].off = ref >= 0 ? fp->off : -fp->off;
      }
      nh = nr > 0 ? nr : 0;
      aggr8 = leaf_solid != NULL && nh > 0;
    }
    if (nh <= 0 && !over8) continue;
    /* УСЛОВИЕ 1:1 (записано в заголовке): грань сетки у разрезанной ячейки
     * обязана совпадать с гранью коробки, иначе флюидную часть пришлось бы ещё
     * и обрезать прямоугольником. Проверяется, а не предполагается. */
    for (int32_t k = m->fstart[c]; k < m->fstart[c + 1]; k++) {
      const tr3_face *fa = &m->f[m->flist[k]];
      if (fa->hi[0] - fa->lo[0] != m->csize[c] || fa->hi[1] - fa->lo[1] != m->csize[c]) {
        cu->nbad++;
      }
    }
    for (int j = 0; j < nh; j++) {
      hid[j] = CUT3_INNER;
      hflip[j] = CUT3_INNER;
    }
    int32_t lo[3], hi[3];
    for (int a = 0; a < 3; a++) {
      lo[a] = m->clo[c][a];
      hi[a] = lo[a] + m->csize[c];
    }
    /* ФЛЮИД = ДОПОЛНЕНИЕ материала, точное разбиение на выпуклые куски */
    hz_poly3 *pieces = calloc(
        (size_t)(nh > 0 ? nh : 1),
        sizeof(
            hz_poly3)); /* nh=0 при over8: не нулевой буфер — анализатор строил ложный over-read */
    if (pieces == NULL) {
      tr3_cut_free(cu);
      return 1;
    }
    double c3[3];
    for (int a = 0; a < 3; a++)
      c3[a] = (double)m->clo[c][a] + 0.5 * (double)m->csize[c];
    double hh = (double)m->csize[c];
    if (aggr8) {
      /* §782: ОБЪЁМ ГРУБОЙ ЯЧЕЙКИ — ТОЧНОЙ АГРЕГАЦИЕЙ по листьям; грани
       * остаются ПОЛНЫМИ (дефолт — упрощение §782, названо в плане: листовая
       * флюидная обрезка граней не агрегируется, вклад дальней зоны < 0.01 %
       * по А1212). Элементы кусков строит piece-цикл ниже; union-рез большой
       * коробки не гоняется вовсе (предел §773). */
      memset(cu->mvol[c], 0, 16 * sizeof(double));
      if (aggr_mvol_rec(cu, m, ft, cm, leaf_solid, lsctx, m->node[c], lo, m->csize[c], c3, hh, c)) {
        free(pieces);
        tr3_cut_free(cu);
        return 1;
      }
    }
    int npc = 0;
    if (!aggr8 && !over8 &&
        hz_poly3_complement(pieces, nh, &npc, lo, hi, h, hid, hflip, nh) != HZ_P3_OK)
      npc = 0;
    if (npc > nh)
      npc = 0; /* невозможно по контракту complement (max = nh);
                * страховка от ложного over-read анализатора */
    /* Р-8: запись целиком из ОГРАНИЧЕННЫХ кусков? (нужно ниже дважды) */
    int allb8 = nh > 0 ? 1 : 0;
    for (int j8 = 0; j8 < nh && allb8; j8++) {
      int32_t r8 = rec != NULL ? cm->fref[rec->f0 + j8] : refs772[j8];
      int32_t f8 = r8 >= 0 ? r8 : ~r8;
      if (!(f8 >= 0 && f8 < ft->n && ft->f[f8].bounded)) allb8 = 0;
    }
    int degen8 = 0;
    if (!aggr8 && !over8 && npc == 0) {
      if (allb8 && (solid_in == NULL || !solid_in[c])) {
        /* Р-8: union-рез КУСОЧНОЙ записи съел коробку целиком — у кусков
         * объёма нет, «сплошная» здесь ложь того же класса, что дыры §779
         * (замерено: топ-дыры cavity с долей флюида ровно 0.000). Ячейка
         * остаётся ПОЛНОЙ коробкой (дефолты mvol/ffm не перезаписываются),
         * элементы кусков строятся ниже. Толстые тела сюда не попадают —
         * их метит solid_in. */
        degen8 = 1;
      } else { /* флюида нет: ячейка целиком в материале */
        cu->solid[c] = 1;
        memset(cu->mvol[c], 0, 16 * sizeof(double));
        for (int32_t k = m->fstart[c]; k < m->fstart[c + 1]; k++) {
          int32_t fi = m->flist[k];
          cu->farea[fi] = 0.0;
          memset(cu->ffm[fi], 0, 16 * sizeof(double));
          memset(cu->ffmb[fi], 0, 16 * sizeof(double));
          memset(cu->ffmx[fi], 0, 16 * sizeof(double));
        }
        free(pieces);
        continue;
      }
    }

    if (!degen8 && !aggr8 && !over8) {
      memset(cu->mvol[c], 0, 16 * sizeof(double));
      for (int32_t k = m->fstart[c]; k < m->fstart[c + 1]; k++) {
        int32_t fi = m->flist[k];
        cu->farea[fi] = 0.0;
        memset(cu->ffm[fi], 0, 16 * sizeof(double));
        memset(cu->ffmb[fi], 0, 16 * sizeof(double));
        memset(cu->ffmx[fi], 0, 16 * sizeof(double));
      }
    }

    /* ПОВЕРХНОСТЬ БЕРЁТСЯ ИЗ САМОГО МАТЕРИАЛА, А НЕ ИЗ ДОПОЛНЕНИЯ. Это правка по
     * измерению: грань куска j дополнения лежит на плоскости h_j, но граничит
     * она с множеством {h_0..h_j}, куда входят и ПОЗДНИЕ куски дополнения, а не
     * только материал. Поэтому такие грани — НЕ граница материала, и площадь
     * поверхности выходила завышенной: замкнутость флюида ломалась на 9.9 при
     * площади грани порядка 1. Материал же даёт свои грани прямо: у hz_poly3_cut
     * грань с fsrc >= 0 И ЕСТЬ кусок фасета. */
    if (!aggr8 && !over8) { /* §782: у агрегата все фасеты кусочные — mat-путь пуст */
      hz_poly3 mat;
      int32_t hid2[HZ_P3_MAXH] = {0};
      for (int j = 0; j < nh; j++)
        hid2[j] = j;
      if (hz_poly3_cut(&mat, lo, hi, h, hid2, nh) == HZ_P3_OK) {
        for (int32_t fj = 0; fj < mat.nf; fj++) {
          int32_t src = mat.fsrc[fj];
          if (src < 0) continue; /* грань коробки: это не поверхность */
          int32_t b0 = mat.floff[fj], nv = mat.floff[fj + 1] - b0;
          if (nv < 3 || nv > 32) continue;
          double vw[CUT3_CLIPV][3];
          for (int32_t e = 0; e < nv; e++)
            for (int a = 0; a < 3; a++)
              vw[e][a] = m->fr.o[a] + m->fr.u[a] * mat.v[mat.fl[b0 + e]][a];
          int32_t ref = rec != NULL ? cm->fref[rec->f0 + src] : refs772[src]; /* §772 */
          int32_t fidx = ref >= 0 ? ref : ~ref;
          /* Р-8 (Д3): элементы bounded-фасетов строятся ПРЯМЫМ клипом куска к
           * коробке (цикл ниже) — грань вырожденного материала пропускается */
          if (fidx >= 0 && fidx < ft->n && ft->f[fidx].bounded) continue;
          tr3_selem se;
          memset(&se, 0, sizeof se);
          se.cell = c;
          se.facet = fidx;
          se.area = poly_mass2(m, vw, (int)nv, c, c, se.m);
          if (!(se.area > 0.0)) continue;
          /* вершины — проекционному сбору; переполнение ПОМЕЧАЕТСЯ, а не режется */
          if (nv <= TR3_SE_MAXV) {
            se.nv = nv;
            for (int32_t e = 0; e < nv; e++)
              for (int a = 0; a < 3; a++)
                se.v[e][a] = vw[e][a];
          } else {
            se.nv = 0;
            cu->nsebig++;
          }
          /* НУЛЕВОЙ ВЕКТОР матрицы масс — прямо из уравнения плоскости (К39).
           * Не вычисляется численно и не угадывается: он ЗАДАН геометрией. */
          se.nul[0] = -h[src].off;
          for (int a = 0; a < 3; a++) {
            se.nul[0] += h[src].n[a] * ((double)m->clo[c][a] + 0.5 * hh);
            se.nul[a + 1] = hh * h[src].n[a];
          }
          /* Материал держит {n·x ≤ off}, значит НАРУЖУ материала (в флюид)
           * смотрит +n. В мир — делением на u и пере-нормировкой (Г21). */
          double nw[3], nm2 = 0.0;
          for (int a = 0; a < 3; a++) {
            nw[a] = h[src].n[a] / m->fr.u[a];
            nm2 += nw[a] * nw[a];
          }
          nm2 = sqrt(nm2);
          for (int a = 0; a < 3; a++)
            se.n[a] = nw[a] / nm2;
          if (push_se(cu, &se)) {
            free(pieces);
            tr3_cut_free(cu);
            return 1;
          }
        }
      }
    }

    /* Р-8 (Д3): элементы ОГРАНИЧЕННЫХ кусков — прямым клипом куска к коробке.
     * Дубль фасета в записи (обе ориентации, §772-слэб) даёт ОДИН элемент.
     * §794: цикл идёт ПО ЗАПИСИ (не по h[]) — предела HZ_P3_MAXH у элементов
     * нет; ориентированные нормаль и нулевой вектор строятся из ФАСЕТА той же
     * формулой, что hz_cutmap_hspaces (К39/Г21). */
    int32_t nref8 = rec != NULL ? rec->nf : nh;
    for (int32_t j8 = 0; j8 < nref8; j8++) {
      int32_t ref8 = rec != NULL ? cm->fref[rec->f0 + j8] : refs772[j8];
      int32_t fi8 = ref8 >= 0 ? ref8 : ~ref8;
      if (fi8 < 0 || fi8 >= ft->n || !ft->f[fi8].bounded) continue;
      int dup8 = 0;
      for (int32_t q8 = 0; q8 < j8 && !dup8; q8++) {
        int32_t r2 = rec != NULL ? cm->fref[rec->f0 + q8] : refs772[q8];
        if ((r2 >= 0 ? r2 : ~r2) == fi8) dup8 = 1;
      }
      if (dup8) continue;
      double vw8[CUT3_CLIPV][3];
      int nv8 = cut3_piece_in_box(m, &ft->f[fi8], lo, hi, vw8);
      if (nv8 < 3) continue;
      tr3_selem se;
      memset(&se, 0, sizeof se);
      se.cell = c;
      se.facet = fi8;
      se.area = poly_mass2(m, vw8, nv8, c, c, se.m);
      if (!(se.area > 0.0)) continue;
      if (nv8 <= TR3_SE_MAXV) {
        se.nv = nv8;
        for (int e = 0; e < nv8; e++)
          for (int a = 0; a < 3; a++)
            se.v[e][a] = vw8[e][a];
      } else {
        se.nv = 0;
        cu->nsebig++;
      }
      double hn8[3] = {0, 0, 0}, hoff8; /* нули — ложный класс анализатора (тернарник в цикле) */
      const hz_facet *fp8 = &ft->f[fi8];
      for (int a = 0; a < 3; a++)
        hn8[a] = ref8 >= 0 ? fp8->n[a] : -fp8->n[a];
      hoff8 = ref8 >= 0 ? fp8->off : -fp8->off;
      se.nul[0] = -hoff8;
      for (int a = 0; a < 3; a++) {
        se.nul[0] += hn8[a] * ((double)m->clo[c][a] + 0.5 * hh);
        se.nul[a + 1] = hh * hn8[a];
      }
      double nw8[3], nm8 = 0.0;
      for (int a = 0; a < 3; a++) {
        nw8[a] = hn8[a] / m->fr.u[a];
        nm8 += nw8[a] * nw8[a];
      }
      nm8 = sqrt(nm8);
      for (int a = 0; a < 3; a++)
        se.n[a] = nw8[a] / nm8;
      if (push_se(cu, &se)) {
        free(pieces);
        tr3_cut_free(cu);
        return 1;
      }
    }

    for (int pi = 0; pi < npc; pi++) {
      hz_poly3 *p = &pieces[pi];
      double mm[4][4];
      if (tr3_mass_matrix(p, &m->fr, c3, hh, mm) == 0)
        for (int i = 0; i < 4; i++)
          for (int j = 0; j < 4; j++)
            cu->mvol[c][i][j] += mm[i][j];

      for (int32_t fj = 0; fj < p->nf; fj++) {
        int32_t src = p->fsrc[fj];
        if (src == CUT3_INNER) continue; /* внутренняя грань флюида */
        int32_t b0 = p->floff[fj], nv = p->floff[fj + 1] - b0;
        if (nv < 3 || nv > 32) continue;
        double vw[32][3];
        for (int32_t e = 0; e < nv; e++)
          for (int a = 0; a < 3; a++)
            vw[e][a] = m->fr.o[a] + m->fr.u[a] * p->v[p->fl[b0 + e]][a];

        if (src >= 0) continue; /* не грань коробки — к поверхности отношения не имеет */
        {                       /* грань КОРОБКИ: ищем соответствующую грань сетки */
          int bf = (int)(~src);
          if (bf < 0 || bf > 5) continue;
          int axis = bf / 2;
          int32_t pos = m->clo[c][axis] + ((bf & 1) ? m->csize[c] : 0);
          /* УСЛОВИЕ 1:1 СНЯТО — К96: СТОРОНА КОРОБКИ МОЖЕТ БЫТЬ РАЗБИТА НА
           * НЕСКОЛЬКО ГРАНЕЙ, И ФЛЮИД ОБРЕЗАЕТСЯ ПО КАЖДОЙ.
           *
           * Прежде здесь стоял `break` на ПЕРВОЙ подходящей грани: весь флюидный
           * многоугольник стороны сваливался в неё одну, а остальные подграни
           * получали ноль. На равномерной сетке подгрань ровно одна, и дефект не
           * проявлялся ВОВСЕ; на градуированной он ломал поток, и `cut3` честно
           * считал такие ячейки в `nbad`, но не исправлял.
           *
           * Обрезка — Сазерленд—Ходжмен по четырём полуплоскостям прямоугольника
           * грани, в ДВУХ поперечных координатах (третья постоянна на грани).
           * Многоугольник флюида выпуклый, прямоугольник выпуклый, значит
           * пересечение выпукло и алгоритма достаточно. Порога здесь нет:
           * стороны прямоугольника — целые координаты сетки, и сравнение идёт с
           * ними, а не с подобранным числом. */
          int ua = (axis + 1) % 3, wa = (axis + 2) % 3;
          if (ua > wa) {
            int tt = ua;
            ua = wa;
            wa = tt;
          }
          for (int32_t k = m->fstart[c]; k < m->fstart[c + 1]; k++) {
            int32_t fi = m->flist[k];
            if (m->f[fi].axis != axis || m->f[fi].pos != pos) continue;
            /* прямоугольник грани в МИРЕ по двум поперечным осям */
            double rlo[2] = {m->fr.o[ua] + m->fr.u[ua] * (double)m->f[fi].lo[0],
                             m->fr.o[wa] + m->fr.u[wa] * (double)m->f[fi].lo[1]};
            double rhi[2] = {m->fr.o[ua] + m->fr.u[ua] * (double)m->f[fi].hi[0],
                             m->fr.o[wa] + m->fr.u[wa] * (double)m->f[fi].hi[1]};
            double cp[2][64][3];
            int cur = 0, ncp = (int)nv;
            for (int e = 0; e < ncp; e++)
              for (int a = 0; a < 3; a++)
                cp[0][e][a] = vw[e][a];
            for (int side = 0; side < 4 && ncp >= 3; side++) {
              int ax = (side < 2) ? ua : wa;
              double lim = (side & 1) ? rhi[side < 2 ? 0 : 1] : rlo[side < 2 ? 0 : 1];
              double sgn = (side & 1) ? -1.0 : 1.0; /* внутри: sgn*(x−lim) ≥ 0 */
              int nxt = 0, oth = 1 - cur;
              for (int e = 0; e < ncp; e++) {
                const double *A = cp[cur][e], *B = cp[cur][(e + 1) % ncp];
                double da = sgn * (A[ax] - lim), db = sgn * (B[ax] - lim);
                if (da >= 0.0 && nxt < 64) {
                  for (int a = 0; a < 3; a++)
                    cp[oth][nxt][a] = A[a];
                  nxt++;
                }
                if ((da > 0.0 && db < 0.0) || (da < 0.0 && db > 0.0)) {
                  double tt = da / (da - db);
                  if (nxt < 64) {
                    for (int a = 0; a < 3; a++)
                      cp[oth][nxt][a] = A[a] + tt * (B[a] - A[a]);
                    nxt++;
                  }
                }
              }
              cur = oth;
              ncp = nxt;
            }
            if (ncp < 3) continue;
            double t1[4][4];
            cu->farea[fi] += poly_mass2(m, cp[cur], ncp, m->f[fi].ca, m->f[fi].ca, t1);
            for (int i = 0; i < 4; i++)
              for (int j = 0; j < 4; j++)
                cu->ffm[fi][i][j] += t1[i][j];
            poly_mass2(m, cp[cur], ncp, m->f[fi].ca, m->f[fi].cb, t1);
            for (int i = 0; i < 4; i++)
              for (int j = 0; j < 4; j++)
                cu->ffmx[fi][i][j] += t1[i][j];
            if (m->f[fi].cb >= 0) {
              poly_mass2(m, cp[cur], ncp, m->f[fi].cb, m->f[fi].cb, t1);
              for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                  cu->ffmb[fi][i][j] += t1[i][j];
            }
          }
        }
      }
    }
    free(pieces);
  }

  if (solid_in != NULL)
    for (int32_t c = 0; c < m->ncell; c++) {
      if (!cu->solid[c]) continue;
      for (int32_t k = m->fstart[c]; k < m->fstart[c + 1]; k++) {
        int32_t fi = m->flist[k];
        int mine_is_a = m->f[fi].ca == c;
        /* обнуляется ТОЛЬКО своя сторона: у соседа-флюида его сторона грани
         * остаётся, а поток через неё запирает поверхность в его же ячейке */
        if (mine_is_a)
          memset(cu->ffm[fi], 0, 16 * sizeof(double));
        else
          memset(cu->ffmb[fi], 0, 16 * sizeof(double));
        memset(cu->ffmx[fi], 0, 16 * sizeof(double));
      }
    }

  /* CSR по ячейкам */
  cu->sestart = calloc((size_t)m->ncell + 1, sizeof(int32_t));
  if (cu->sestart == NULL) {
    tr3_cut_free(cu);
    return 1;
  }
  for (int32_t i = 0; i < cu->nse; i++)
    cu->sestart[cu->se[i].cell + 1]++;
  for (int32_t c = 0; c < m->ncell; c++)
    cu->sestart[c + 1] += cu->sestart[c];
  cu->selist = calloc((size_t)(cu->nse > 0 ? cu->nse : 1), sizeof(int32_t));
  int32_t *fill = calloc((size_t)m->ncell, sizeof(int32_t));
  if (cu->selist == NULL || fill == NULL) {
    free(fill);
    tr3_cut_free(cu);
    return 1;
  }
  for (int32_t i = 0; i < cu->nse; i++) {
    int32_t c = cu->se[i].cell;
    cu->selist[cu->sestart[c] + fill[c]++] = i;
  }
  free(fill);
  return 0;
}
