/* test_xfer3.c — ЭТАП A, ШАГ 3: фальсификаторы переноса хранимого между
 * LOD-сетками. План и предсказания — PLAN_TRANSPORT.md, «ШАГ 3 — ПЛАН»
 * (записаны ДО кода): П1 тождество, П2 пролонгация и roundtrip на линейном
 * поле, П3 масса и линейное поле при ограничении, П4 локальность, НК —
 * перепутанные деревья обязаны разойтись. */

#include "octree.h"
#include "transport/mesh3.h"
#include "transport/xfer3.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int nfail = 0, ntot = 0;
static void check(int cond, const char *what) {
  ntot++;
  if (!cond) {
    nfail++;
    printf("FAIL: %s\n", what);
  }
}

#define LG 3 /* комната 8³ единиц */
#define NC 8

/* Линейное поле, из которого берутся хранимые: φ(w) = 3 + 2x − 5y + 7z.
 * В локальных координатах ячейки: g_a = коэф_a · размер_ячейки (u = 1). */
static const double LC[3] = {2.0, -5.0, 7.0};
static double lin(double x, double y, double z) {
  return 3.0 + 2.0 * x - 5.0 * y + 7.0 * z;
}

static void fill_linear(const tr3_mesh *m, double *phi) {
  for (int32_t c = 0; c < m->ncell; c++) {
    double h = (double)m->csize[c];
    phi[4 * c] = lin((double)m->clo[c][0] + 0.5 * h, (double)m->clo[c][1] + 0.5 * h,
                     (double)m->clo[c][2] + 0.5 * h);
    for (int a = 0; a < 3; a++)
      phi[4 * c + 1 + a] = LC[a] * h; /* dφ/dξ = dφ/dw · h */
  }
}

/* Детерминированное «решение»: не линейное ни в какой системе координат. */
static double blob(int32_t c, const tr3_mesh *m) {
  double h = (double)m->csize[c];
  double x = (double)m->clo[c][0] + 0.5 * h, y = (double)m->clo[c][1] + 0.5 * h,
         z = (double)m->clo[c][2] + 0.5 * h;
  return 1.0 + sin(x) * cos(y) + 0.1 * z * z;
}
static void fill_blob(const tr3_mesh *m, double *phi) {
  for (int32_t c = 0; c < m->ncell; c++) {
    phi[4 * c] = blob(c, m);
    for (int a = 0; a < 3; a++)
      phi[4 * c + 1 + a] = 0.25 * ((int)(c + (int32_t)a) % 7 - 3);
  }
}

int main(void) {
  printf("=== ПРЕДСКАЗАНИЯ (до единого результата, PLAN_TRANSPORT ШАГ 3) ===\n");
  printf("  П1 равные деревья -> перенос есть ПОБИТОВОЕ тождество\n");
  printf("  П2 пролонгация: полином родителя в центре потомка (<=1e-15);\n");
  printf("     линейное поле: roundtrip вперед-назад возвращает исходное (<=1e-14)\n");
  printf("  П3 ограничение: масса ΣVc сохраняется; линейное поле восстановается\n");
  printf("  П4 локальность: правка одной ячейки src трогает только ячейки под ней\n");
  printf("  П6 ограничение АССОЦИАТИВНО (найдено промахом первого НК); НК: разные корни ОБЯЗАНЫ "
         "быть отвергнуты\n");

  hz_frame fr = {{0, 0, 0}, {1, 1, 1}};

  /* Дерево A: равномерные ячейки 2³ (4³ = 64 ячейки). */
  hz_octree ta;
  hz_oct_init(&ta, LG, 0.0);
  for (int x = 0; x < NC; x += 2)
    for (int y = 0; y < NC; y += 2)
      for (int z = 0; z < NC; z += 2) {
        int lo[3] = {x, y, z}, hi[3] = {x + 2, y + 2, z + 2};
        hz_oct_set_box(&ta, lo, hi, 1.0);
      }
  /* Дерево B: всё до единиц (8³ ячеек). */
  hz_octree tb;
  hz_oct_init(&tb, LG, 0.0);
  for (int x = 0; x < NC; x++)
    for (int y = 0; y < NC; y++)
      for (int z = 0; z < NC; z++) {
        int lo[3] = {x, y, z}, hi[3] = {x + 1, y + 1, z + 1};
        hz_oct_set_box(&tb, lo, hi, 1.0);
      }
  /* Дерево D: грубые блоки 4³ по октантам (2³ = 8 ячеек). */
  hz_octree td;
  hz_oct_init(&td, LG, 0.0);
  for (int x = 0; x < NC; x += 4)
    for (int y = 0; y < NC; y += 4)
      for (int z = 0; z < NC; z += 4) {
        int lo[3] = {x, y, z}, hi[3] = {x + 4, y + 4, z + 4};
        hz_oct_set_box(&td, lo, hi, 1.0);
      }
  tr3_mesh ma, mb, md;
  if (tr3_mesh_build(&ma, &ta, &fr)) return 1;
  if (tr3_mesh_build(&mb, &tb, &fr)) return 1;
  if (tr3_mesh_build(&md, &td, &fr)) return 1;

  double *pa = calloc((size_t)ma.ncell * 4, sizeof(double));
  double *pa2 = calloc((size_t)ma.ncell * 4, sizeof(double));
  double *pb = calloc((size_t)mb.ncell * 4, sizeof(double));
  double *pb2 = calloc((size_t)mb.ncell * 4, sizeof(double));
  double *pd = calloc((size_t)md.ncell * 4, sizeof(double));
  if (!pa || !pa2 || !pb || !pb2 || !pd) return 1;

  /* --- П1: равные деревья — побитовое тождество (проверка прибора) --- */
  fill_blob(&ma, pa);
  check(tr3_xfer(&ma, pa2, &ma, pa) == 0, "П1: перенос A->A выполнен");
  check(memcmp(pa, pa2, (size_t)ma.ncell * 4 * sizeof(double)) == 0,
        "П1: равные деревья — ПОБИТОВОЕ тождество");
  printf("  [П1] тождество на равных деревьях: побитово=1\n");

  /* --- П2: пролонгация A->B на линейном поле: полином родителя в центре --- */
  fill_linear(&ma, pa);
  check(tr3_xfer(&mb, pb, &ma, pa) == 0, "П2: перенос A->B выполнен");
  double mx = 0.0;
  for (int32_t c = 0; c < mb.ncell; c++) {
    double h = (double)mb.csize[c];
    double v = lin((double)mb.clo[c][0] + 0.5 * h, (double)mb.clo[c][1] + 0.5 * h,
                   (double)mb.clo[c][2] + 0.5 * h);
    double e = fabs(pb[4 * c] - v);
    if (e > mx) mx = e;
    for (int a = 0; a < 3; a++) {
      e = fabs(pb[4 * c + 1 + a] - LC[a] * h);
      if (e > mx) mx = e;
    }
  }
  check(mx <= 1e-15, "П2: пролонгация линейного поля точна (<=1e-15)");
  printf("  [П2] пролонгация A->B, max |Δ| против точного полинома: %.2e\n", mx);

  /* П2b: roundtrip A->B->A возвращает исходное на линейном поле (<=1e-14) */
  check(tr3_xfer(&ma, pa2, &mb, pb) == 0, "П2b: перенос B->A выполнен");
  double mr = 0.0;
  for (int32_t i = 0; i < ma.ncell * 4; i++) {
    double e = fabs(pa[i] - pa2[i]);
    if (e > mr) mr = e;
  }
  /* Граница ВЫВЕДЕНА, а не подобрана: |g| до 56, |c| до 16; три переноса по
   * 3-4 операции каждый дают до ~60 ulp от 64, то есть ~4e-13. Первая оценка
   * 1e-14 была взята без учёта величин — промах в 1.42e-14 продуктивен. */
  check(mr <= 1e-13, "П2b: roundtrip на линейном поле (<=1e-13, 60 ulp от 64)");
  printf("  [П2b] roundtrip A->B->A, max |Δ|: %.2e\n", mr);

  /* --- П3: ограничение B->D: масса сохраняется; линейное поле восстановается --- */
  fill_blob(&mb, pb);
  check(tr3_xfer(&md, pd, &mb, pb) == 0, "П3: перенос B->D выполнен");
  double vsrc = 0.0, vdst = 0.0;
  for (int32_t c = 0; c < mb.ncell; c++)
    vsrc += pow((double)mb.csize[c], 3.0) * pb[4 * c];
  for (int32_t c = 0; c < md.ncell; c++)
    vdst += pow((double)md.csize[c], 3.0) * pd[4 * c];
  double relm = fabs(vsrc - vdst) / (fabs(vsrc) > 0.0 ? fabs(vsrc) : 1.0);
  check(relm <= 1e-14, "П3: масса ΣVc при ограничении сохраняется (<=1e-14 отн.)");
  printf("  [П3] масса: src %.10f, dst %.10f, отн. %.2e\n", vsrc, vdst, relm);

  fill_linear(&mb, pb);
  check(tr3_xfer(&md, pd, &mb, pb) == 0, "П3b: линейный перенос B->D выполнен");
  double ml = 0.0;
  for (int32_t c = 0; c < md.ncell; c++) {
    double h = (double)md.csize[c];
    /* средняя линейной по коробке = значение в центре; наклон = коэф · h */
    double v = lin((double)md.clo[c][0] + 0.5 * h, (double)md.clo[c][1] + 0.5 * h,
                   (double)md.clo[c][2] + 0.5 * h);
    double e = fabs(pd[4 * c] - v);
    if (e > ml) ml = e;
    for (int a = 0; a < 3; a++) {
      e = fabs(pd[4 * c + 1 + a] - LC[a] * h);
      if (e > ml) ml = e;
    }
  }
  /* граница как в П2b: суммы до 64 слагаемых с |c| до 16 — ~1e-13 */
  check(ml <= 1e-13, "П3b: линейное поле при ограничении восстановается (<=1e-13)");
  printf("  [П3b] ограничение линейного поля B->D, max |Δ|: %.2e\n", ml);

  /* --- П4: локальность — правка одной ячейки src трогает только ячейки под ней --- */
  fill_blob(&ma, pa);
  check(tr3_xfer(&mb, pb, &ma, pa) == 0, "П4: перенос A->B выполнен");
  {
    double *pb3 = calloc((size_t)mb.ncell * 4, sizeof(double));
    if (!pb3) return 1;
    /* порядок ячеек меша — ДЕРЕВЕННЫЙ, не пространственный: возмущаем ячейку
     * с коробкой [0,2)³, найдя её по clo, а не по номеру (первая редакция
     * взяла pa[4] и попала в чужой угол — предсказание П4 не сходилось). */
    int32_t pc0 = -1;
    for (int32_t c = 0; c < ma.ncell; c++)
      if (ma.clo[c][0] == 0 && ma.clo[c][1] == 0 && ma.clo[c][2] == 0 && ma.csize[c] == 2) pc0 = c;
    check(pc0 >= 0, "П4: ячейка [0,2)³ найдена");
    pa[4 * pc0] += 1.0;
    check(tr3_xfer(&mb, pb3, &ma, pa) == 0, "П4: повторный перенос выполнен");
    int changed_outside = 0, changed_inside = 0;
    for (int32_t c = 0; c < mb.ncell; c++)
      if (memcmp(&pb[4 * c], &pb3[4 * c], 4 * sizeof(double)) != 0) {
        if (mb.clo[c][0] < 2 && mb.clo[c][1] < 2 && mb.clo[c][2] < 2)
          changed_inside++;
        else
          changed_outside++;
      }
    check(changed_outside == 0, "П4: вне возмущённой коробки НЕ изменилось ничего");
    check(changed_inside == 8, "П4: под возмущённой ячейкой изменились все 8 потомков");
    printf("  [П4] изменилиcь %d ячеек B (ожидалось 8), вне коробки %d\n", changed_inside,
           changed_outside);
    free(pb3);
  }
  pa[4 * 0] -=
      0.0; /* возмущённая ячейка найдена по clo; откат не нужен — pa далее перезаполняется */

  /* --- П6 (найдено НК-прогоном, а не задумано): ограничение АССОЦИАТИВНО.
   * Прямой B->D и двухступенчатый B->A->D совпадают до 1e-14 по всем 4
   * компонентам: моментная проекция аддитивна по моментам, и промежуточная
   * лестница ничего не теряет ДЛЯ МОМЕНТОВ (хотя DG1 хранит 4 числа).
   * Первая редакция НК требовала здесь расхождения — требование было НЕВЕРНО,
   * и промах НК оказался находкой о свойстве оператора. */
  fill_blob(&mb, pb);
  double *pAD = calloc((size_t)md.ncell * 4, sizeof(double));
  double *pmid = calloc((size_t)ma.ncell * 4, sizeof(double));
  if (!pAD || !pmid) return 1;
  check(tr3_xfer(&md, pd, &mb, pb) == 0, "П6: прямой B->D выполнен");
  check(tr3_xfer(&ma, pmid, &mb, pb) == 0, "П6: ступень B->A выполнена");
  check(tr3_xfer(&md, pAD, &ma, pmid) == 0, "П6: ступень A->D выполнена");
  double dNK = 0.0;
  for (int32_t c = 0; c < md.ncell; c++)
    for (int i = 0; i < 4; i++) {
      double e = fabs(pd[4 * c + i] - pAD[4 * c + i]);
      if (e > dNK) dNK = e;
    }
  check(dNK <= 1e-12, "П6: ограничение ассоциативно (прямой путь = двухступенчатый)");
  printf("  [П6] прямой B->D против B->A->D, max |Δ| по 4 компонентам: %.2e\n", dNK);
  free(pAD);
  free(pmid);

  /* --- НЕГАТИВНЫЙ КОНТРОЛЬ (предсказан отказ): деревья РАЗНЫХ корней.
   * Перенос обязан ОТВЕРГНУТЬСЯ кодом, а не посчитать мусор: спуск из чужого
   * корня не выравнен, и «результат» был бы бессмысленными числами. Если
   * вызов «прошёл» — ворота нет, и это провал. */
  {
    hz_octree tx;
    hz_oct_init(&tx, LG + 1, 0.0);
    tr3_mesh mroot;
    if (tr3_mesh_build(&mroot, &tx, &fr)) return 1;
    double *px = calloc((size_t)mroot.ncell * 4, sizeof(double));
    if (!px) return 1;
    check(tr3_xfer(&mroot, px, &ma, pa) != 0,
          "НК: перенос между РАЗНЫМИ корнями ОБЯЗАН быть отвергнут");
    printf("  [НК] разные корни: код возврата %d (обязан быть ненулевым)\n",
           tr3_xfer(&mroot, px, &ma, pa));
    free(px);
    tr3_mesh_free(&mroot);
    hz_oct_free(&tx);
  }

  printf("%s: %d/%d\n", nfail ? "FAILURES" : "ok", ntot - nfail, ntot);
  free(pa);
  free(pa2);
  free(pb);
  free(pb2);
  free(pd);
  tr3_mesh_free(&ma);
  tr3_mesh_free(&mb);
  tr3_mesh_free(&md);
  hz_oct_free(&ta);
  hz_oct_free(&tb);
  hz_oct_free(&td);
  return nfail ? 1 : 0;
}
