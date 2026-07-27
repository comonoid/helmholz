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
#define CUT3_SURF0 (1000)

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

int tr3_cut_build(tr3_cut *cu, const tr3_mesh *m, const hz_facettab *ft, const hz_cutmap *cm,
                  const uint8_t *solid_in) {
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
    if (rec == NULL) continue;
    hz_hspace h[HZ_P3_MAXH];
    int32_t hid[HZ_P3_MAXH], hflip[HZ_P3_MAXH];
    int nh = hz_cutmap_hspaces(ft, cm, rec, h, hid, hflip, HZ_P3_MAXH);
    if (nh <= 0) continue;
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
    hz_poly3 *pieces = calloc((size_t)nh, sizeof(hz_poly3));
    if (pieces == NULL) {
      tr3_cut_free(cu);
      return 1;
    }
    int npc = 0;
    if (hz_poly3_complement(pieces, nh, &npc, lo, hi, h, hid, hflip, nh) != HZ_P3_OK) npc = 0;
    if (npc == 0) { /* флюида нет: ячейка целиком в материале */
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

    memset(cu->mvol[c], 0, 16 * sizeof(double));
    for (int32_t k = m->fstart[c]; k < m->fstart[c + 1]; k++) {
      int32_t fi = m->flist[k];
      cu->farea[fi] = 0.0;
      memset(cu->ffm[fi], 0, 16 * sizeof(double));
      memset(cu->ffmb[fi], 0, 16 * sizeof(double));
      memset(cu->ffmx[fi], 0, 16 * sizeof(double));
    }

    double c3[3];
    for (int a = 0; a < 3; a++)
      c3[a] = (double)m->clo[c][a] + 0.5 * (double)m->csize[c];
    double hh = (double)m->csize[c];

    /* ПОВЕРХНОСТЬ БЕРЁТСЯ ИЗ САМОГО МАТЕРИАЛА, А НЕ ИЗ ДОПОЛНЕНИЯ. Это правка по
     * измерению: грань куска j дополнения лежит на плоскости h_j, но граничит
     * она с множеством {h_0..h_j}, куда входят и ПОЗДНИЕ куски дополнения, а не
     * только материал. Поэтому такие грани — НЕ граница материала, и площадь
     * поверхности выходила завышенной: замкнутость флюида ломалась на 9.9 при
     * площади грани порядка 1. Материал же даёт свои грани прямо: у hz_poly3_cut
     * грань с fsrc >= 0 И ЕСТЬ кусок фасета. */
    {
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
          double vw[32][3];
          for (int32_t e = 0; e < nv; e++)
            for (int a = 0; a < 3; a++)
              vw[e][a] = m->fr.o[a] + m->fr.u[a] * mat.v[mat.fl[b0 + e]][a];
          tr3_selem se;
          memset(&se, 0, sizeof se);
          se.cell = c;
          int32_t ref = cm->fref[rec->f0 + src];
          se.facet = ref >= 0 ? ref : ~ref;
          se.area = poly_mass2(m, vw, (int)nv, c, c, se.m);
          if (!(se.area > 0.0)) continue;
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
          for (int32_t k = m->fstart[c]; k < m->fstart[c + 1]; k++) {
            int32_t fi = m->flist[k];
            if (m->f[fi].axis != axis || m->f[fi].pos != pos) continue;
            double t1[4][4];
            cu->farea[fi] += poly_mass2(m, vw, (int)nv, m->f[fi].ca, m->f[fi].ca, t1);
            for (int i = 0; i < 4; i++)
              for (int j = 0; j < 4; j++)
                cu->ffm[fi][i][j] += t1[i][j];
            poly_mass2(m, vw, (int)nv, m->f[fi].ca, m->f[fi].cb, t1);
            for (int i = 0; i < 4; i++)
              for (int j = 0; j < 4; j++)
                cu->ffmx[fi][i][j] += t1[i][j];
            if (m->f[fi].cb >= 0) {
              poly_mass2(m, vw, (int)nv, m->f[fi].cb, m->f[fi].cb, t1);
              for (int i = 0; i < 4; i++)
                for (int j = 0; j < 4; j++)
                  cu->ffmb[fi][i][j] += t1[i][j];
            }
            break;
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
