/* plod — ИЕРАРХИЧЕСКОЕ ОГРУБЛЕНИЕ И СРЕЗ LOD (PLAN_ELEMENTS.md, §56, О16).
 *
 * Что печатается и почему именно это:
 *   - ПО УРОВНЯМ: число элементов и отношение к предыдущему. Отклонение от
 *     четырёх — ЗАМЕР, а не отклонение от нормы (§54.2): оно показывает, где
 *     сцена перестаёт быть гладкой на этом масштабе. Плюс `dmax` (обязан быть
 *     под допуском уровня) и форма `P/√A` (§55, «близость к правильному»);
 *   - ПО СРЕЗУ при камере §2: число элементов и метрика О20 через ОБЩИЙ набор
 *     лучей против нулевого уровня. Ставка шага (§56, предсказание 3): срез даёт
 *     в 2…10 раз меньше элементов, чем ЛУЧШИЙ однородный уровень при равной
 *     ошибке — поэтому рядом печатаются и однородные уровни, и срез;
 *   - СЧЁТЧИКИ ПОСТРОЙКИ: отказы ворот, тождественные продвижения, прирост
 *     площади (проверка §54.2), пары, смежные по ребру.
 *
 * ВЛОЖЕННОСТЬ ПРОВЕРЯЕТСЯ ЦЕЛОЧИСЛЕННО (А132): у каждого узла множество
 * исходных полигонов обязано быть объединением множеств детей без пересечений.
 * Проверка по МЕТКАМ, а не по площадям: «сходится по площади» даёт ложное
 * подтверждение при перекрытии.
 */

/* §327: сколько размеров ячейки в лестнице касаний — двоичная, 12 октав от
 * габарита сцены. Двенадцать, а не «сколько-нибудь»: у города габарит 79 м, и
 * 12 октав доводят ячейку до 19 мм, то есть ниже входного треугольника. */
#define HZ_NH 12

#include "pedge.h"
#include "lodio.h"
#include "plod.h"
#include "poly_seg.h"
#include "polygon.h"
#include "pray.h"
#include "pvfit.h"
#include "scene_cfg.h"
#include "scene_obj.h"
#include "transport/cam3.h"
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

static int cmp_dbl(const void *x, const void *y) {
  double a = *(const double *)x, b = *(const double *)y;
  return (a < b) ? -1 : ((a > b) ? 1 : 0);
}

static double pct(double *v, int64_t n, double p) {
  if (n <= 0) return 0.0;
  int64_t k = (int64_t)(p * (double)(n - 1));
  if (k < 0) k = 0;
  if (k >= n) k = n - 1;
  return v[k];
}

/* Метрика О20 через общий набор лучей: две сцены, разница по ЛУЧУ. */
static void metric(const hz_polyset *pf, const hz_polyset *pc, const tr3_camera *cam, double eps_px,
                   const char *tag, int32_t nelem) {
  hz_pray gf, gc;
  if (hz_pray_build(&gf, pf, 4.0) != 0) return;
  if (hz_pray_build(&gc, pc, 4.0) != 0) {
    hz_pray_free(&gf);
    return;
  }
  const int W = cam->w, H = cam->h;
  int64_t nb = 0, nl = 0, ng = 0;
  double *dt = malloc((size_t)W * (size_t)H * sizeof *dt);
  double *dp = malloc((size_t)W * (size_t)H * sizeof *dp);
  /* СИЛУЭТНЫЙ КАНАЛ (§69). Пиксельная метрика мерит СМЕЩЕНИЕ ПЛОСКОСТИ вдоль
   * луча, а ломается картинка на СИЛУЭТЕ — там, где огрублённый край врезается в
   * очертание или отступает от него. До сих пор это схлопывалось в проценты
   * «потеряли/приобрели», то есть в площадь, а видна глазу ШИРИНА полосы.
   * Ширина берётся без единого порога: `(потеряли + приобрели) / длина силуэта`,
   * где длина силуэта — число пикселей эталона, у которых сосед по четырём
   * направлениям УЛЕТЕЛ В НЕБО. Это внешний силуэт, определённый точно: небо есть
   * промах луча, а не значение глубины, и порога тут нет по построению. */
  uint8_t *hit = calloc((size_t)W * (size_t)H, 1);
  double *tref = calloc((size_t)W * (size_t)H, sizeof *tref);
  int32_t *pref = malloc((size_t)W * (size_t)H * sizeof *pref);
  if (dt == NULL || dp == NULL || hit == NULL || tref == NULL || pref == NULL) {
    free(dt);
    free(dp);
    free(hit);
    free(tref);
    free(pref);
    hz_pray_free(&gf);
    hz_pray_free(&gc);
    return;
  }
  for (int j = 0; j < H; j++)
    for (int i = 0; i < W; i++) {
      double o[3], d[3], tf = 0.0, tc = 0.0;
      tr3_camera_ray(cam, i, j, o, d);
      int32_t hf = hz_pray_hit(&gf, o, d, 1e-6, &tf);
      int32_t hc = hz_pray_hit(&gc, o, d, 1e-6, &tc);
      if (hf >= 0) {
        hit[(size_t)j * (size_t)W + (size_t)i] = 1;
        tref[(size_t)j * (size_t)W + (size_t)i] = tf;
        pref[(size_t)j * (size_t)W + (size_t)i] = hf;
      }
      if (hf < 0 && hc < 0) continue;
      if (hf >= 0 && hc < 0) {
        nl++;
        continue;
      }
      if (hf < 0) {
        ng++;
        continue;
      }
      double e = fabs(tf - tc);
      dt[nb] = e;
      const double *nf = pf->p[hf].n;
      double cr = fabs(nf[0] * d[0] + nf[1] * d[1] + nf[2] * d[2]);
      dp[nb] = e * cr / (eps_px * (tf > 0.0 ? tf : 1e-9));
      nb++;
    }
  qsort(dt, (size_t)nb, sizeof *dt, cmp_dbl);
  qsort(dp, (size_t)nb, sizeof *dp, cmp_dbl);
  int64_t tot = nb + nl + ng;
  /* ХВОСТ ПЕЧАТАЕТСЯ ОБЯЗАТЕЛЬНО, И ЭТО НЕ УКРАШЕНИЕ (А158). Первый городской
   * прогон дал у среза и у ПЕРЕМЕШАННОГО контроля одинаковые `p50` и `p90` — оба
   * ноль, — из чего следовало бы, что критерий среза пуст. Но порча от
   * перестановки затрагивает малую долю пикселей, и `p90` её не видит по
   * построению: при 99.86 % целых пикселей девяностый процентиль стоит глубоко
   * внутри целой части. Различать обязаны `p99`, `p99.9`, максимум и ДОЛЯ
   * пикселей грубее одного — они и печатаются. */
  int64_t nbad = 0;
  for (int64_t i = 0; i < nb; i++)
    if (dp[i] > 1.0) nbad++;
  /* ДЛИНА СИЛУЭТА (А166, вторая редакция). Первая считала соседство с промахом
   * луча и дала ровно `4·512 − 4 = 2044` на ОБЕИХ сценах — то есть периметр кадра:
   * сцены заполняют кадр, неба нет. Вторая брала скачок глубины больше `ε·t` и
   * пометила 65 % кадра: на СКОЛЬЗЯЩЕМ взгляде непрерывная поверхность меняет
   * глубину быстрее пиксельного следа, и порог верен лишь при взгляде в лоб.
   * Третья, эта: сосед считается разрывом, если его точка попадания лежит ВНЕ
   * ПЛОСКОСТИ текущего элемента дальше пиксельного следа `ε·t`. Тогда наклон
   * поверхности учтён точно — по её собственной плоскости, — и порога нет: `ε·t`
   * есть предел разрешения кадра. */
  int64_t nsil = 0;
  for (int j = 0; j < H; j++)
    for (int i = 0; i < W; i++) {
      size_t q = (size_t)j * (size_t)W + (size_t)i;
      if (!hit[q]) continue;
      const hz_poly *P = &pf->p[pref[q]];
      double o[3], d[3];
      tr3_camera_ray(cam, i, j, o, d);
      int edge = 0;
      const int di[4] = {-1, 1, 0, 0}, dj[4] = {0, 0, -1, 1};
      for (int k = 0; k < 4 && !edge; k++) {
        int ii = i + di[k], jj = j + dj[k];
        if (ii < 0 || ii >= W || jj < 0 || jj >= H) continue;
        size_t r = (size_t)jj * (size_t)W + (size_t)ii;
        if (!hit[r]) {
          edge = 1;
          break;
        }
        double o2[3], d2[3];
        tr3_camera_ray(cam, ii, jj, o2, d2);
        double dev = 0.0;
        for (int c = 0; c < 3; c++)
          dev += P->n[c] * (o2[c] + tref[r] * d2[c]);
        dev = fabs(dev - P->off);
        if (dev > eps_px * tref[q]) edge = 1;
      }
      if (edge) nsil++;
    }
  printf("      %-28s элем %7d | попали %6.2f%% потеряли %5.2f%% приобрели %5.2f%% | "
         "ПИКСЕЛИ p50 %6.2f p90 %6.2f p99 %7.2f p99.9 %8.2f макс %8.2f | >1пкс %6.3f%% | "
         "СИЛУЭТ: длина %6lld пкс, сдвиг %6.2f пкс\n",
         tag, nelem, 100.0 * (double)nb / (double)(tot ? tot : 1),
         100.0 * (double)nl / (double)(tot ? tot : 1), 100.0 * (double)ng / (double)(tot ? tot : 1),
         pct(dp, nb, 0.5), pct(dp, nb, 0.9), pct(dp, nb, 0.99), pct(dp, nb, 0.999),
         pct(dp, nb, 1.0), 100.0 * (double)nbad / (double)(nb ? nb : 1), (long long)nsil,
         (nsil > 0) ? (double)(nl + ng) / (double)nsil : 0.0);
  /* СБРОС БУФЕРА ПОСЛЕ КАЖДОЙ СТРОКИ. Без него вывод, перенаправленный в файл,
   * копится блоками, и прогон, убитый на середине (перезагрузка машины 07-31),
   * оставляет ноль метрических строк при часе работы. Строка метрики стоит
   * минуты — сброс на её фоне бесплатен. */
  fflush(stdout);
  free(dt);
  free(dp);
  free(hit);
  free(tref);
  free(pref);
  hz_pray_free(&gf);
  hz_pray_free(&gc);
}

/* ПОТОЛОК КОМПЛАНАРНОГО СЛИЯНИЯ (§70). Сколько РАЗЛИЧНЫХ плоскостей в сцене —
 * это нижняя граница числа элементов, достижимая слиянием С НУЛЕВОЙ ошибкой, если
 * бы радиус поиска не мешал. Считается сортировкой по квантованной плоскости, а не
 * кластеризацией: квантование РАЗРЕЗАЕТ группы на границах корзин, поэтому число
 * получается ЗАВЫШЕННЫМ, то есть оценка потолка КОНСЕРВАТИВНА, и это её главное
 * свойство. Печатается по нескольким допускам сразу — единственного «правильного»
 * тут нет, и подбирать его под ответ нельзя. */
typedef struct {
  int64_t k0, k1, k2, k3;
} plkey;

static int cmp_plkey(const void *x, const void *y) {
  const plkey *a = (const plkey *)x, *b = (const plkey *)y;
  if (a->k0 != b->k0) return (a->k0 < b->k0) ? -1 : 1;
  if (a->k1 != b->k1) return (a->k1 < b->k1) ? -1 : 1;
  if (a->k2 != b->k2) return (a->k2 < b->k2) ? -1 : 1;
  if (a->k3 != b->k3) return (a->k3 < b->k3) ? -1 : 1;
  return 0;
}

/* УГОЛ МЕЖДУ СОСЕДНИМИ УЧАСТКАМИ (§70.1, требование пользователя 07-30: «угол
 * считать только между соседними полигонами, объединять имеет смысл только
 * соседние»). Прежний замер считал КОМПЛАНАРНОСТЬ и дал ноль: смежные
 * компланарные грани сегментация уже слила, их не осталось. Но вопрос был не про
 * компланарность, а про УГОЛ — и распределение углов по СМЕЖНЫМ парам не мерилось
 * ни разу. Гистограмма ниже и есть ответ на «есть ли чем работать угловой
 * сортировке»: сколько смежных пар лежит ниже каждого порога.
 *
 * Двугранный угол печатается в шкале пользователя: `180°` — плоско, вниз к `90°`
 * — выпуклый излом, вверх к `270°` — вогнутый. Знак берётся признаком
 * `(n_A − n_B)·(c_B − c_A)`. */
static int32_t uf_find(int32_t *par, int32_t x) {
  while (par[x] != x) {
    par[x] = par[par[x]];
    x = par[x];
  }
  return x;
}

typedef struct {
  int64_t key;
  int32_t tri;
} edgerec;

static int cmp_edge(const void *x, const void *y) {
  const edgerec *a = (const edgerec *)x, *b = (const edgerec *)y;
  return (a->key < b->key) ? -1 : ((a->key > b->key) ? 1 : 0);
}

static double adj_angles(const hz_objmesh *m, const hz_pseglist *sg, const hz_polyset *ps) {
  edgerec *er = malloc((size_t)m->nt * 3 * sizeof *er);
  int32_t *par = malloc((size_t)ps->np * sizeof *par);
  if (er == NULL || par == NULL) {
    free(er);
    free(par);
    return 0.0;
  }
  for (int32_t i = 0; i < ps->np; i++)
    par[i] = i;
  int64_t ne = 0;
  for (int32_t t = 0; t < m->nt; t++)
    for (int e = 0; e < 3; e++) {
      int32_t v0 = m->f[(size_t)t * 3 + (size_t)e],
              v1 = m->f[(size_t)t * 3 + (size_t)((e + 1) % 3)];
      int32_t lo = (v0 < v1) ? v0 : v1, hi = (v0 < v1) ? v1 : v0;
      er[ne].key = (int64_t)lo * 2147483647ll + (int64_t)hi;
      er[ne].tri = t;
      ne++;
    }
  qsort(er, (size_t)ne, sizeof *er, cmp_edge);
  /* Пороги — та же двоично-десятичная шкала, что у полос §59, и выбраны они до
   * замера, а не под него. */
  const double thr[8] = {0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 45.0, 91.0};
  int64_t below[8] = {0, 0, 0, 0, 0, 0, 0, 0};
  /* ПОЛНАЯ ГИСТОГРАММА `0…180°` ПО 5° (А173). Прежние пороги кончались на `91°`, и
   * весь хвост выше сваливался в одну корзину: отличить `135°` (скат к стене) от
   * `180°` (лицо к изнанке) было нечем, а я при этом утверждал, что в сцене «всё
   * либо 90, либо 180». Утверждение шло сверх измеренного; корзины заведены. */
  int64_t hist[36];
  for (int q = 0; q < 36; q++)
    hist[q] = 0;
  int64_t nadjpair = 0, nconv = 0, nconc = 0;
  double amin = 999.0;
  for (int64_t i = 1; i < ne; i++) {
    if (er[i].key != er[i - 1].key) continue;
    int32_t la = sg->label[er[i].tri], lb = sg->label[er[i - 1].tri];
    if (la < 0 || lb < 0 || la >= ps->np || lb >= ps->np || la == lb) continue;
    const hz_poly *A = &ps->p[la], *B = &ps->p[lb];
    double cs = A->n[0] * B->n[0] + A->n[1] * B->n[1] + A->n[2] * B->n[2];
    if (cs > 1.0) cs = 1.0;
    if (cs < -1.0) cs = -1.0;
    double th = acos(cs) * 180.0 / 3.14159265358979323846;
    double dd = 0.0;
    for (int c = 0; c < 3; c++)
      dd += (A->n[c] - B->n[c]) * (B->org[c] - A->org[c]);
    nadjpair++;
    if (dd < 0.0)
      nconv++;
    else
      nconc++;
    if (th < amin) amin = th;
    for (int q = 0; q < 8; q++)
      if (th <= thr[q]) below[q]++;
    {
      int hb = (int)(th / 5.0);
      if (hb < 0) hb = 0;
      if (hb > 35) hb = 35;
      hist[hb]++;
    }
    /* Заодно: во сколько групп схлопнулась бы сцена при пороге 45°, если сливать
     * ТОЛЬКО соседей. Это верхняя оценка для угловой тактики при связности. */
    if (th <= 45.0) {
      int32_t ra = uf_find(par, la), rb = uf_find(par, lb);
      if (ra != rb) par[ra] = rb;
    }
  }
  int64_t ng = 0;
  for (int32_t i = 0; i < ps->np; i++)
    if (uf_find(par, i) == i) ng++;
  printf("   УГОЛ МЕЖДУ СМЕЖНЫМИ УЧАСТКАМИ (§70.1): смежных пар %lld, из них выпуклых %lld, "
         "вогнутых %lld; минимальный излом %.4f°\n",
         (long long)nadjpair, (long long)nconv, (long long)nconc, amin);
  printf("      доля пар с изломом не более:");
  for (int q = 0; q < 8; q++)
    printf(" %.1f°:%.2f%%", thr[q], 100.0 * (double)below[q] / (double)(nadjpair ? nadjpair : 1));
  printf("\n      ПОЛНАЯ гистограмма по 5°, %% (только непустые):");
  for (int q = 0; q < 36; q++)
    if (hist[q] > 0)
      printf(" [%d-%d°] %.2f", q * 5, q * 5 + 5,
             100.0 * (double)hist[q] / (double)(nadjpair ? nadjpair : 1));
  printf("\n      если слить ВСЕ смежные пары с изломом до 45°: %lld групп (сокращение %.1fx)\n",
         (long long)ng, (double)ps->np / (double)(ng ? ng : 1));
  /* КАЛИБРОВКА УГЛОВОЙ ЛЕСТНИЦЫ ПО СЦЕНЕ (§76). Возвращается нижняя граница первой
   * корзины, где лежит хотя бы `1 %` смежных пар, — то есть излом, с которого у
   * сцены НАЧИНАЕТСЯ масса. Полураствор конуса у пары с изломом `θ` равен `θ/2`,
   * поэтому первая ступень лестницы есть половина этой величины. Один процент —
   * не подобранный порог, а требование «корзина не пуста статистически»: при
   * меньшей доле ступень обслуживала бы единицы пар из миллионов. */
  double a1 = 0.0;
  {
    /* ПЕРВЫЙ ПРОЦЕНТИЛЬ ПО НАКОПЛЕНИЮ, А НЕ ПЕРВАЯ КОРЗИНА С МАССОЙ (§76, поправка
     * по замеру). Первая редакция брала корзину, где лежит хотя бы процент пар, —
     * то есть ОСНОВНУЮ МАССУ, — и на зале выбрала `70°`, проскочив все пологие
     * стыки, которых там 0.44 % ниже пяти градусов. Лестница обязана начинаться у
     * НИЖНЕГО края распределения, иначе повторяется ровно та ошибка, ради которой
     * калибровка и заводится. Один процент — требование статистической
     * непустоты: ниже него ступень обслуживала бы единицы пар из миллионов. */
    int64_t need = (nadjpair + 99) / 100, acc = 0;
    for (int q = 0; q < 36; q++) {
      acc += hist[q];
      if (acc >= need) {
        a1 = (double)q * 5.0;
        break;
      }
    }
  }
  printf("      КАЛИБРОВКА: масса углов начинается с %.1f°, первая ступень лестницы %.2f°\n", a1,
         0.5 * a1);
  free(er);
  free(par);
  return 0.5 * a1;
}

/* СКОЛЬКО ГРАНЕЙ СМОТРИТ В СТЕНУ (§70.2). Третий путь для города, где пологих
 * стыков нет вовсе: грань, упирающаяся в соседний блок, не видна ниоткуда и не
 * участвует ни в переносе, ни в картинке — её удаление стоит РОВНО НОЛЬ ошибки.
 * Меряется прямо: из центра элемента вдоль его ВНЕШНЕЙ нормали пускается луч, и
 * печатается распределение расстояния до первого попадания. Порога здесь нет —
 * печатается кривая, а «сколько удалить» решается по ней, а не до неё. */
static void buried_faces(const hz_polyset *ps) {
  hz_pray g;
  if (hz_pray_build(&g, ps, 4.0) != 0) return;
  const double thr[6] = {0.01, 0.05, 0.2, 1.0, 5.0, 20.0};
  int64_t below[6] = {0, 0, 0, 0, 0, 0};
  int64_t nfree = 0, ntot = 0;
  for (int32_t i = 0; i < ps->np; i++) {
    const hz_poly *P = &ps->p[i];
    if (!(P->area > 0.0)) continue;
    double o[3], d[3];
    for (int c = 0; c < 3; c++) {
      d[c] = P->n[c];
      o[c] = P->org[c] + 1e-6 * P->n[c];
    }
    double t = 0.0;
    int32_t h = hz_pray_hit(&g, o, d, 1e-9, &t);
    ntot++;
    if (h < 0) {
      nfree++;
      continue;
    }
    for (int q = 0; q < 6; q++)
      if (t <= thr[q]) below[q]++;
  }
  printf("   ГРАНИ, СМОТРЯЩИЕ В СТЕНУ (§70.2): элементов %lld; луч вдоль внешней нормали "
         "не встретил ничего у %lld (%.2f%%)\n",
         (long long)ntot, (long long)nfree, 100.0 * (double)nfree / (double)(ntot ? ntot : 1));
  printf("      доля с препятствием ближе:");
  for (int q = 0; q < 6; q++)
    printf(" %5.2f м:%.2f%%", thr[q], 100.0 * (double)below[q] / (double)(ntot ? ntot : 1));
  printf("\n");
  hz_pray_free(&g);
}

/* ЗАПЕРТЫЕ ГРАНИ (§70.3, предложение пользователя 07-30: «находясь на расстоянии
 * от закрытой коробки, ты в принципе не можешь видеть, что там внутри»).
 *
 * ЧТО МЕРИТСЯ И ЧТО НЕТ. Из центра элемента пускается `K` лучей по полусфере его
 * нормали; если НИ ОДИН не ушёл наружу, элемент заперт — снаружи он не виден ни
 * при какой камере. Это ЗАМЕР ПОТЕНЦИАЛА, а не критерий удаления, и разница
 * принципиальна: выборка по направлениям может ПРОПУСТИТЬ узкую щель (окно,
 * дверь, зазор), и тогда грань, на самом деле видимую, мы объявим запертой. Для
 * УДАЛЕНИЯ нужен консервативный тест, который щель пропустить не может, — заливка
 * пустого пространства от внешней границы с шагом мельче самого узкого проёма.
 * Пользователь предупредил об этом прямо («чтобы не было там полупрозрачных стен
 * или окон»), и предупреждение записано здесь, а не в докладе задним числом. */
static void enclosed_faces(const hz_polyset *ps, int K) {
  hz_pray g;
  if (hz_pray_build(&g, ps, 4.0) != 0) return;
  int64_t nlock = 0, ntot = 0, nopen = 0;
  double alock = 0.0, atot = 0.0;
  for (int32_t i = 0; i < ps->np; i++) {
    const hz_poly *P = &ps->p[i];
    if (!(P->area > 0.0)) continue;
    ntot++;
    atot += P->area;
    int esc = 0;
    for (int k = 0; k < K && !esc; k++) {
      /* Направления — детерминированной решёткой Фибоначчи по полусфере: замер
       * обязан повторяться, поэтому никакого генератора. */
      double u = ((double)k + 0.5) / (double)K;
      double z = u;                 /* косинусное распределение не нужно: вопрос */
      double r = sqrt(1.0 - z * z); /* бинарный — ушёл или нет */
      double ph = 2.39996322972865332 * (double)k;
      double lu = r * cos(ph), lv = r * sin(ph);
      double eu[3], ev[3];
      int ax = 0;
      for (int c = 1; c < 3; c++)
        if (fabs(P->n[c]) < fabs(P->n[ax])) ax = c;
      double t0[3] = {0.0, 0.0, 0.0};
      t0[ax] = 1.0;
      eu[0] = P->n[1] * t0[2] - P->n[2] * t0[1];
      eu[1] = P->n[2] * t0[0] - P->n[0] * t0[2];
      eu[2] = P->n[0] * t0[1] - P->n[1] * t0[0];
      double en = sqrt(eu[0] * eu[0] + eu[1] * eu[1] + eu[2] * eu[2]);
      if (!(en > 0.0)) continue;
      for (int c = 0; c < 3; c++)
        eu[c] /= en;
      ev[0] = P->n[1] * eu[2] - P->n[2] * eu[1];
      ev[1] = P->n[2] * eu[0] - P->n[0] * eu[2];
      ev[2] = P->n[0] * eu[1] - P->n[1] * eu[0];
      double d[3], o[3];
      for (int c = 0; c < 3; c++) {
        d[c] = z * P->n[c] + lu * eu[c] + lv * ev[c];
        o[c] = P->org[c] + 1e-5 * P->n[c];
      }
      double t = 0.0;
      if (hz_pray_hit(&g, o, d, 1e-9, &t) < 0) esc = 1;
    }
    if (esc)
      nopen++;
    else {
      nlock++;
      alock += P->area;
    }
  }
  printf("   ЗАПЕРТЫЕ ГРАНИ (§70.3, %d луча по полусфере): заперто %lld из %lld (%.2f%%), "
         "по площади %.2f%%\n",
         K, (long long)nlock, (long long)ntot, 100.0 * (double)nlock / (double)(ntot ? ntot : 1),
         100.0 * alock / (atot > 0.0 ? atot : 1.0));
  printf("      ЭТО ПОТЕНЦИАЛ, А НЕ КРИТЕРИЙ: выборка направлений может пропустить узкую щель, "
         "и тогда видимая грань объявлена запертой. Удалять — только по заливке пустоты.\n");
  hz_pray_free(&g);
}

/* РАЗМЕР ГРАНИ И ПРЕДСКАЗАНИЕ ПУСТЫХ УРОВНЕЙ (§70.4, вывод пользователя 07-30:
 * «это просто геометрия такая; предсказуемо, что первые два уровня работать не
 * будут»).
 *
 * ПОЧЕМУ ЭТО ПРЕДСКАЗАНИЕ, А НЕ НАБЛЮДЕНИЕ. Замер §70.1 дал минимальный излом
 * между СМЕЖНЫМИ участками ровно `90°` на городе. Две перпендикулярные грани
 * размера `s`, слитые в одну плоскость, дают отклонение около `s/2`: подогнанная
 * плоскость режет угол по диагонали. Значит уровень с допуском `δ` не может слить
 * НИЧЕГО, пока `δ < s/2` для типичной грани, и число пустых уровней считается
 * заранее — по распределению размеров, а не по прогону лестницы.
 * Размер берётся как `sqrt(площадь)`: у почти квадратной грани это её сторона, а
 * подгонять более хитрую меру не под что. */
static double face_sizes(const hz_polyset *ps, double delta0) {
  double *sz = malloc((size_t)ps->np * sizeof *sz);
  if (sz == NULL) return 0.0;
  int64_t n = 0;
  for (int32_t i = 0; i < ps->np; i++)
    if (ps->p[i].area > 0.0) sz[n++] = sqrt(ps->p[i].area);
  if (n == 0) {
    free(sz);
    return 0.0;
  }
  qsort(sz, (size_t)n, sizeof *sz, cmp_dbl);
  double p10 = sz[(int64_t)(0.10 * (double)(n - 1))], p50 = sz[(int64_t)(0.50 * (double)(n - 1))];
  double p90 = sz[(int64_t)(0.90 * (double)(n - 1))];
  printf("   РАЗМЕР ГРАНИ (§70.4): sqrt(площадь) p10 %.3f, p50 %.3f, p90 %.3f м\n", p10, p50, p90);
  /* Уровень L имеет допуск `delta0·2^L`; первый работающий — тот, где допуск
   * дошёл до половины МЕДИАННОЙ грани. */
  int lfirst = 0;
  while (delta0 * pow(2.0, (double)lfirst) < 0.5 * p50 && lfirst < 32)
    lfirst++;
  int lp10 = 0;
  while (delta0 * pow(2.0, (double)lp10) < 0.5 * p10 && lp10 < 32)
    lp10++;
  printf("      ПРЕДСКАЗАНИЕ: при δ0 = %.3f м первые %d уровней обязаны быть ПУСТЫМИ "
         "(допуск ниже половины медианной грани %.3f м); самые мелкие грани (p10) "
         "начнут сливаться с уровня %d\n",
         delta0, lfirst, 0.5 * p50, lp10);
  free(sz);
  return 0.5 * p50;
}

static void plane_ceiling(const hz_polyset *ps) {
  const double ang[3] = {0.1, 0.5, 2.0};    /* градусы */
  const double off[3] = {0.01, 0.05, 0.20}; /* метры */
  plkey *k = malloc((size_t)ps->np * sizeof *k);
  if (k == NULL) return;
  printf("   ПОТОЛОК КОМПЛАНАРНОГО СЛИЯНИЯ (§70): элементов %d; различных плоскостей —\n", ps->np);
  for (int a = 0; a < 3; a++)
    for (int o = 0; o < 3; o++) {
      double qa = ang[a] * 3.14159265358979323846 / 180.0;
      for (int32_t i = 0; i < ps->np; i++) {
        const hz_poly *P = &ps->p[i];
        k[i].k0 = (int64_t)floor(P->n[0] / qa);
        k[i].k1 = (int64_t)floor(P->n[1] / qa);
        k[i].k2 = (int64_t)floor(P->n[2] / qa);
        k[i].k3 = (int64_t)floor(P->off / off[o]);
      }
      qsort(k, (size_t)ps->np, sizeof *k, cmp_plkey);
      int64_t nd = (ps->np > 0) ? 1 : 0;
      for (int32_t i = 1; i < ps->np; i++)
        if (cmp_plkey(&k[i], &k[i - 1]) != 0) nd++;
      printf("      угол %4.1f°, смещение %5.2f м: %8lld  (сокращение до %5.1f×)\n", ang[a], off[o],
             (long long)nd, (double)ps->np / (double)(nd ? nd : 1));
    }
  free(k);
}

int main(int argc, char **argv) {
  int city = 0, miguel = 0, miglow = 0, levels = 0;
  for (int i = 1; i < argc; i++) {
    if (strcmp(argv[i], "city") == 0) city = 1;
    if (strcmp(argv[i], "miguel") == 0) miguel = 1;
    if (strcmp(argv[i], "miglow") == 0) miguel = miglow = 1;
  }
  /* Допуск сегментации. Сан-Мигелю берётся ГОРОДСКОЙ `0.05` м, а не зальный: он
   * в метрах, как и город, тогда как зальные `45` мм привязаны к масштабу зала
   * `0.003` (тот же довод, что в `tools/pcoarse.c`). */
  double dseg = (city || miguel) ? 0.05 : 0.045;
  int simp = 0, vfit = 0, maxlev = 9, uniform = 0, viamerge = 0, bands = 0, bycount = 0,
      byangle = 0, curve = 0, auto0 = 0, nosin = 0;
  double epsmul = 1.0, radmul = 1.0, eyemul = 1.0, ang0 = 0.0, vfitlim = 0.0;
  double cgate = 0.0, ladbase = 0.0;
  const char *save = NULL, *load = NULL;
  /* НЕИЗВЕСТНЫЙ АРГУМЕНТ — ОШИБКА, А НЕ ПРОПУСК (§60, дефект оснастки). Флаг,
   * который не совпал, молчал, и конфигурация вышла тождественной другой; поймать
   * это удалось лишь по совпадению всех семи чисел. Обрыв дешевле. */
  for (int i = 1; i < argc; i++) {
    int ok = 0;
    if (strcmp(argv[i], "city") == 0 || strcmp(argv[i], "hall") == 0) ok = 1;
    /* ТРЕТЬЯ СЦЕНА В ЛЕСТНИЧНОМ СТЕНДЕ (§321). До сих пор Сан-Мигель знали
     * только `pcoarse`/`pcell`/`pstow`, и лестница на нём не мерилась ни разу —
     * а §319 требует условие постройки грубых уровней на ВСЕХ ТРЁХ сценах. */
    if (strcmp(argv[i], "miguel") == 0 || strcmp(argv[i], "miglow") == 0) ok = 1;
    /* ТОЛЬКО ТАБЛИЦА УРОВНЕЙ (§321). Срез, лучевая метрика и зондирование сцены
     * лучами (§70.2/§70.3) к выходу этого шага отношения не имеют, а на
     * Сан-Мигеле стоят дороже самого замера: метрика гоняет лучи против
     * `9.96` млн треугольников дважды. */
    if (strcmp(argv[i], "levels") == 0) ok = levels = 1;
    if (strcmp(argv[i], "simp") == 0) ok = simp = 1;
    if (strcmp(argv[i], "vfit") == 0) ok = vfit = 1;
    /* §73: ПРЕДЕЛ СМЕЩЕНИЯ ВЕРШИНЫ, МЕТРЫ. Прежде он брался как допуск САМОГО
     * ГРУБОГО уровня (`0.045·2⁶ = 2.88 м` на зале) — то есть вершине разрешалось
     * уехать на два метра в комнате, и подгонка портила качество в восемь раз:
     * `p99` с 0.08 до 35.01 пикселя. Величина обязана быть параметром и
     * выбираться замером, а не наследоваться от глубины лестницы. */
    if (strncmp(argv[i], "vfit=", 5) == 0) {
      vfitlim = strtod(argv[i] + 5, NULL);
      ok = vfit = 1;
    }
    if (strcmp(argv[i], "uniform") == 0) ok = uniform = 1;
    /* Второй построитель: уровень строится `hz_merge` над предыдущим (§58). */
    if (strcmp(argv[i], "viamerge") == 0) ok = viamerge = 1;
    if (strcmp(argv[i], "bands") == 0) ok = bands = 1;
    /* ВТОРАЯ ТАКТИКА (§60): уровень задан ЧИСЛОМ (четверть), допуск не
     * ограничивает, достигнутое отклонение печатается замером. */
    if (strcmp(argv[i], "bycount") == 0) ok = bycount = 1;
    /* ТРЕТЬЯ ТАКТИКА (§67): критерий уровня — УГОЛ (полураствор конуса нормалей). */
    if (strcmp(argv[i], "byangle") == 0) ok = byangle = 1;
    /* §76: АВТОКАЛИБРОВКА ОБЕИХ ЛЕСТНИЦ ПО РАСПРЕДЕЛЕНИЯМ СЦЕНЫ. Обе величины —
     * первый допуск и первый угол — задавались руками, и оба раза мимо всего, что
     * в сцене есть: у города четыре нижних уровня пусты по построению, а угловая
     * лестница доходила до первого кандидата лишь на последнем уровне. Следствия я
     * оба раза записал как свойство сцены или тактики. Здесь они берутся из
     * измеренных распределений: половина медианной грани и половина угла, с
     * которого начинается масса. */
    if (strcmp(argv[i], "auto0") == 0) ok = auto0 = 1;
    /* §79: критерий среза по ПОЛНОМУ смещению, без множителя sin угла взгляда. */
    if (strcmp(argv[i], "nosin") == 0) ok = nosin = 1;
    /* §77: ворота при цели по числу, в допусках уровня. */
    if (strncmp(argv[i], "cgate=", 6) == 0) {
      cgate = strtod(argv[i] + 6, NULL);
      ok = 1;
    }
    /* §78: основание лестницы (было жёстко 2). */
    if (strncmp(argv[i], "base=", 5) == 0) {
      ladbase = strtod(argv[i] + 5, NULL);
      ok = 1;
    }
    /* §69: КРИВАЯ «ЭЛЕМЕНТЫ ПРОТИВ ОШИБКИ» ПО СВИПУ ε. Одна точка тактики не
     * решает ничего: на зале срез оставляет 865 полигонов из 993 на нулевом
     * уровне, то есть сравниваются сцены, совпадающие на 87 %. Сравнивать надо
     * при РАВНОЙ ЦЕНЕ, а цена задаётся ε. */
    if (strcmp(argv[i], "curve") == 0) ok = curve = 1;
    if (strncmp(argv[i], "lev=", 4) == 0) {
      maxlev = (int)strtol(argv[i] + 4, NULL, 10);
      ok = 1;
    }
    /* §71: ПЕРВАЯ СТУПЕНЬ УГЛОВОЙ ЛЕСТНИЦЫ, ГРАДУСЫ. */
    if (strncmp(argv[i], "ang0=", 5) == 0) {
      ang0 = strtod(argv[i] + 5, NULL);
      ok = 1;
    }
    /* §71: ОТОДВИНУТЬ КАМЕРУ. Глаз уезжает от точки прицела в `eyemul` раз, поле
     * зрения то же. Нужно потому, что при штатной камере `ε·R` по кадру равно
     * размеру грани города (p50 = 1.414 м), то есть LOD меряется там, где ему
     * нечего делать: огрубление законно лишь с `ε·R > s/2`. Свип по `ε` этого не
     * заменяет — он меняет допуск среза, но не то, ЧТО видно и с какого удаления. */
    if (strncmp(argv[i], "eye=", 4) == 0) {
      eyemul = strtod(argv[i] + 4, NULL);
      ok = 1;
    }
    /* §68: РАДИУС ПОИСКА КАНДИДАТОВ В ДОПУСКАХ УРОВНЯ. */
    if (strncmp(argv[i], "rad=", 4) == 0) {
      radmul = strtod(argv[i] + 4, NULL);
      ok = 1;
    }
    if (strncmp(argv[i], "eps=", 4) == 0) {
      epsmul = strtod(argv[i] + 4, NULL);
      ok = 1;
    }
    /* ЛЕСТНИЦА КАК АРТЕФАКТ СЦЕНЫ (§62, пункт 2): строится один раз, дальше
     * читается. `load` заменяет постройку целиком. */
    if (strncmp(argv[i], "save=", 5) == 0) {
      save = argv[i] + 5;
      ok = 1;
    }
    if (strncmp(argv[i], "load=", 5) == 0) {
      load = argv[i] + 5;
      ok = 1;
    }
    if (!ok) {
      fprintf(stderr,
              "plod: неизвестный аргумент «%s»\n"
              "  ожидается: [city|hall|miguel|miglow] [levels] [simp] [vfit] [uniform] [viamerge] "
              "[bands]\n"
              "             [bycount] [byangle] [curve] [rad=X] [eye=X] [ang0=X] [vfit=X] [auto0] "
              "[lev=N] "
              "[eps=X] "
              "[save=Ф] "
              "[load=Ф]\n",
              argv[i]);
      return 2;
    }
  }
  hz_objmesh m;
  const char *scene =
      miguel ? (miglow ? "САН-МИГЕЛЬ low-poly" : "САН-МИГЕЛЬ") : (city ? "ГОРОД" : "зал");
  if (hz_obj_load(&m,
                  miguel ? (miglow ? HZ_CFG_MIGUEL_LOW_OBJ : HZ_CFG_MIGUEL_OBJ)
                         : (city ? HZ_CFG_CITY_OBJ : HZ_CFG_HALL_OBJ),
                  miguel ? HZ_CFG_MIGUEL_SCALE : (city ? HZ_CFG_CITY_SCALE : HZ_CFG_HALL_SCALE)) !=
      0) {
    fprintf(stderr, "нет сцены\n");
    return 1;
  }
  hz_pseglist sg;
  if (hz_seg_planar(&sg, &m, dseg) != 0) return 1;
  hz_polyset psf;
  if (hz_poly_build(&psf, &m, &sg) != 0) return 1;
  if (simp) {
    hz_edgestat es;
    if (hz_edge_simplify(&psf, 0.25 * dseg, HZ_EDGE_SHARED, &es) != 0) return 1;
  }
  plane_ceiling(&psf);
  /* Зондирование сцены лучами (§70.2/§70.3) и угловая калибровка стоят дорого и
   * к таблице уровней отношения не имеют — под `levels` они снимаются. Угол
   * остаётся, если его просит `auto0`: он там вход, а не доклад. */
  double cal_ang = (!levels || auto0) ? adj_angles(&m, &sg, &psf) : 0.0;
  double cal_delta = face_sizes(&psf, dseg);
  if (!levels) {
    buried_faces(&psf);
    enclosed_faces(&psf, 32);
  }
  tr3_camera cam;
  double eyeh[3] = HZ_CFG_HALL_EYE, ath[3] = HZ_CFG_HALL_AT;
  double eyec[3] = HZ_CFG_CITY_EYE, atc[3] = HZ_CFG_CITY_AT, up[3] = HZ_CFG_UP;
  double eyem[3] = HZ_CFG_MIGUEL_EYE, atm[3] = HZ_CFG_MIGUEL_AT;
  double *eyes = miguel ? eyem : (city ? eyec : eyeh);
  const double *ats = miguel ? atm : (city ? atc : ath);
  if (!(eyemul >= 1.0 && eyemul <= 1.0)) {
    for (int c = 0; c < 3; c++)
      eyes[c] = ats[c] + (eyes[c] - ats[c]) * eyemul;
    printf("== КАМЕРА ОТОДВИНУТА в %.1f раз: глаз (%.1f %.1f %.1f)\n", eyemul, eyes[0], eyes[1],
           eyes[2]);
  }
  const int W = 512, H = 512;
  if (tr3_camera_look(&cam, eyes, ats, up, HZ_CFG_FOV_DEG * M_PI / 180.0, W, H) != 0) return 1;
  const double eps_px = (HZ_CFG_FOV_DEG * M_PI / 180.0) / (double)H;
  const double eps = eps_px * epsmul;

  /* Ладдер начинается со СВОЕГО допуска, отдельного от допуска СЕГМЕНТАЦИИ:
   * сегментация задаёт точность базового представления и мельчить её нельзя, а
   * лестнице ниже половины медианной грани делать нечего. */
  double lad0 = dseg;
  if (auto0) {
    if (cal_delta > dseg) lad0 = cal_delta;
    if (ang0 <= 0.0 && cal_ang > 0.0) ang0 = cal_ang;
    printf("== АВТОКАЛИБРОВКА: допуск лестницы %.4f м (сегментация %.4f м), первый угол %.2f°\n",
           lad0, dseg, ang0);
  }
  double t0 = now_s();
  hz_lod L;
  int rc;
  if (load != NULL) {
    /* ЧТЕНИЕ ВМЕСТО ПОСТРОЙКИ. Проверки — внутри `hz_lod_read`, включая полную
     * целочисленную вложенность; отказ печатается СЛОВАМИ, потому что «не ноль»
     * не различает чужой файл, порчу и несогласованность. */
    int fsimp = 0;
    rc = hz_lod_read(&L, load, &sg, m.nt, &fsimp);
    if (rc != HZ_LODIO_OK) {
      fprintf(stderr, "лестница не прочитана: %s (код %d)\n", hz_lodio_str(rc), rc);
      return 1;
    }
    /* УПРОЩЕНИЕ КРАЯ СВЕРЯЕТСЯ ОТДЕЛЬНО: сумма разметки его не видит (см. оговорку
     * в `lodio.h`), а границы полигонов от него зависят. */
    if (fsimp != (simp ? 1 : 0)) {
      fprintf(stderr, "лестница построена %s края, а запрошено %s — отказ\n",
              fsimp ? "С УПРОЩЕНИЕМ" : "БЕЗ УПРОЩЕНИЯ", simp ? "С УПРОЩЕНИЕМ" : "БЕЗ УПРОЩЕНИЯ");
      hz_lod_free(&L);
      return 1;
    }
    /* КАМЕРА СВЕРЯЕТСЯ С ФАЙЛОМ (А131). Лестница строилась под своё `ε`, и срез
     * под другим `ε` дал бы ступени, не совпадающие с построенными. Отказ, а не
     * пересчёт. */
    if (!(fabs(L.eps - eps) <= 1e-12 * (fabs(eps) + 1.0))) {
      fprintf(stderr, "лестница построена под ε = %.9g, а запрошено %.9g — отказ (А131)\n", L.eps,
              eps);
      hz_lod_free(&L);
      return 1;
    }
    dseg = L.delta0;
  } else {
    hz_lodcfg lc;
    memset(&lc, 0, sizeof lc);
    lc.delta0 = lad0;
    lc.eps = eps;
    lc.maxlev = maxlev;
    lc.bands = bands;
    lc.bycount = bycount;
    lc.byangle = byangle;
    lc.radmul = radmul;
    lc.angle0 = ang0;
    lc.cgate = cgate;
    lc.base = ladbase;
    rc = viamerge ? hz_lod_build_merge(&L, &m, &sg, &psf, &lc)
                  : hz_lod_build(&L, &m, &sg, &psf, dseg, maxlev, eps);
    if (rc != 0) {
      fprintf(stderr, "отказ построения лестницы\n");
      return 1;
    }
  }
  double t1 = now_s();
  printf("== ЛЕСТНИЦА: %s, δ0 %g м, уровней %d, узлов %d, %s за %.1f с\n", scene, dseg, L.nlev,
         L.nnd, load ? "ПРОЧИТАНА" : "построена", t1 - t0);
  if (save != NULL) {
    /* ЗАПИСЬ СРАЗУ ПОСЛЕ ПОСТРОЙКИ, ДО всякого замера: если дальше что-то упадёт,
     * двести секунд постройки уже не потеряны. */
    double ts = now_s();
    int wr = hz_lod_write(&L, save, &sg, m.nt, simp);
    if (wr != HZ_LODIO_OK) {
      fprintf(stderr, "лестница не записана: %s (код %d)\n", hz_lodio_str(wr), wr);
      hz_lod_free(&L);
      return 1;
    }
    printf("== ЗАПИСАНА: %s за %.2f с\n", save, now_s() - ts);
    fflush(stdout);
  }
  printf("   постройка: пар по ребру %lld, отказов ворот %lld, тождественных продвижений %lld, "
         "худший прирост площади %.3f\n",
         (long long)L.npair_edge, (long long)L.ngate_rej, (long long)L.nident, L.area_grow);
  fflush(stdout);

  /* --- по уровням --- */
  /* ДОПУСК УРОВНЯ ПЕЧАТАЕТСЯ ФАКТИЧЕСКИЙ, А НЕ `dseg·2^L` (§321, А555). Основание
   * лестницы — параметр (§78), и при `base=1` (негативный контроль этого шага)
   * прежняя строка показывала бы удвоение там, где допуск не растёт вовсе. */
  const double ladb = (ladbase > 0.0) ? ladbase : 2.0;
  int32_t *ncnt = malloc((size_t)(L.nlev > 0 ? L.nlev : 1) * sizeof *ncnt);
  double *nper = malloc((size_t)(L.nlev > 0 ? L.nlev : 1) * sizeof *nper);
  if (ncnt == NULL || nper == NULL) {
    free(ncnt);
    free(nper);
    return 1;
  }
  for (int lev = 0; lev < L.nlev; lev++) {
    ncnt[lev] = 0;
    nper[lev] = 0.0;
  }
  int32_t prevn = 0;
  for (int lev = 0; lev < L.nlev; lev++) {
    const int32_t *lab = L.lab + (size_t)lev * (size_t)L.np;
    int32_t *seen = calloc((size_t)L.nnd, sizeof *seen);
    if (seen == NULL) break;
    int32_t cnt = 0;
    double *dm = malloc((size_t)L.np * sizeof *dm);
    double *sh = malloc((size_t)L.np * sizeof *sh);
    if (dm == NULL || sh == NULL) {
      free(seen);
      free(dm);
      free(sh);
      break;
    }
    int64_t nd = 0;
    for (int32_t k = 0; k < L.np; k++) {
      int32_t id = lab[k];
      if (seen[id]) continue;
      seen[id] = 1;
      cnt++;
      const hz_lodnode *n = &L.nd[id];
      dm[nd] = n->dmax;
      sh[nd] = (n->area_surf > 0.0 && n->perim > 0.0) ? n->perim / sqrt(n->area_surf) : 0.0;
      nper[lev] += n->perim;
      nd++;
    }
    qsort(dm, (size_t)nd, sizeof *dm, cmp_dbl);
    qsort(sh, (size_t)nd, sizeof *sh, cmp_dbl);
    if (L.bycount)
      printf("   уровень %d (ЦЕЛЬ по числу): элементов %d; dmax ЗАМЕР p50 %.4f p90 %.4f макс %.4f; "
             "форма P/√A p50 %.2f p90 %.2f\n",
             lev, cnt, pct(dm, nd, 0.5), pct(dm, nd, 0.9), pct(dm, nd, 1.0), pct(sh, nd, 0.5),
             pct(sh, nd, 0.9));
    else
      printf("   уровень %d (допуск %.4f м): элементов %d; dmax p50 %.4f p90 %.4f макс %.4f; "
             "форма P/√A p50 %.2f p90 %.2f\n",
             lev, L.delta0 * pow(ladb, lev), cnt, pct(dm, nd, 0.5), pct(dm, nd, 0.9),
             pct(dm, nd, 1.0), pct(sh, nd, 0.5), pct(sh, nd, 0.9));
    if (lev > 0 && cnt > 0)
      printf("      сокращение к предыдущему %.2f× (%s)\n", (double)prevn / (double)cnt,
             L.bycount ? "четвёрка — ТРЕБОВАНИЕ: §60"
                       : "четвёрка — структура, не требование: §54.2");
    prevn = cnt;
    ncnt[lev] = cnt;
    fflush(stdout);
    free(seen);
    free(dm);
    free(sh);
  }

  /* --- ЦЕНА ПИРАМИДЫ УРОВНЕЙ (§319, поправка; шаг §321) ---
   *
   * ЧТО ЭТО ЗА ЧИСЛО. Грубые уровни несут СВОЮ геометрию, то есть хранится
   * больше, чем `N_0`; §319 считает цену пирамиды при сокращении вчетверо на
   * уровень равной `1 + 1/4 + 1/16 + … = 4/3`, и это ровно тот предел
   * «20…40 %», который назван допустимым. Здесь та же сумма берётся ЗАМЕРОМ:
   * `Σ_L N_L / N_0`. Все слагаемые считались и раньше, но не печатались ни разу
   * (правило §317).
   *
   * ПОДЛЕСТНИЦЫ — СЕМЕЙСТВО, А НЕ ПОРОГ (А554). Уровень, на котором сливать
   * нечего, копию геометрии заводить не обязан; вместо порога «полезный
   * уровень» печатается цена лестниц с основанием `2`, `4` и `8` (шаг по
   * уровням 1, 2, 3). САМЫЙ ГРУБЫЙ УРОВЕНЬ ВХОДИТ ВО ВСЕ ТРИ — иначе они
   * покрывали бы разный диапазон масштабов и были бы несравнимы. */
  if (L.nlev > 0 && ncnt[0] > 0) {
    printf("   ЦЕНА ПИРАМИДЫ (§319: при сокращении вчетверо на уровень было бы 1.33×)\n");
    printf("      уровней %d, N_0 = %d, самый грубый уровень %d: N = %d, допуск %.4f м\n", L.nlev,
           ncnt[0], L.nlev - 1, ncnt[L.nlev - 1], L.delta0 * pow(ladb, L.nlev - 1));
    printf("      доли N_L/N_0:");
    for (int lev = 0; lev < L.nlev; lev++)
      printf(" %.3f", (double)ncnt[lev] / (double)ncnt[0]);
    printf("\n");
    /* ВТОРАЯ ОСЬ ЦЕНЫ — ПЕРИМЕТР КРАЯ (А550). Счёт КУСКОВ занижает складскую
     * цену: край грубого элемента рваный (`P/√A` растёт по лестнице), а хранится
     * именно край. Периметр к тому же РАЗЛИЧАЕТ два разных слияния, которых счёт
     * не различает: при настоящем огрублении смежные куски сливаются и общий
     * край ИСЧЕЗАЕТ (Σ P падает), а при сборе разрозненной мелочи в одну
     * плоскость границы остаются все до одной (Σ P стоит на месте). */
    if (nper[0] > 0.0) {
      printf("      доли ΣP_L/ΣP_0 (периметр края):");
      for (int lev = 0; lev < L.nlev; lev++)
        printf(" %.3f", nper[lev] / nper[0]);
      printf("\n");
    }
    /* ВЕРШИНЫ КРАЯ — САМА ХРАНИМАЯ ВЕЛИЧИНА, А НЕ ПРОКСИ (А550, А558). Счёт
     * кусков и периметр — две оценки с разных сторон, и на зале они разошлись в
     * три раза (`1.59×` против `5.01×`). Элемент хранится как плоскость плюс
     * КРАЙ, значит спор решает число вершин края, и оно тут считается прямо:
     * разметка уровня → полигоны → `nbv`. Цена — пересборка полигонов на каждом
     * уровне (на Сан-Мигеле около двух минут при прогоне в 38), и она заплачена
     * сознательно: без этого числа вердикт §319 держался бы на прокси. */
    /* ЛЕСТНИЦА РАЗМЕРОВ ЯЧЕЙКИ — ДВОИЧНАЯ, ОТ ГАБАРИТА СЦЕНЫ ВНИЗ (§327). Якоря
     * «своя ячейка уровня» нет СОЗНАТЕЛЬНО: он был бы порогом, подобранным под
     * ответ. Печатается вся матрица, диагональ читает потребитель. */
    double slo[3] = {1e300, 1e300, 1e300}, shi[3] = {-1e300, -1e300, -1e300};
    for (int32_t k = 0; k < m.nv; k++)
      for (int c = 0; c < 3; c++) {
        double x = m.v[3 * (size_t)k + (size_t)c];
        if (x < slo[c]) slo[c] = x;
        if (x > shi[c]) shi[c] = x;
      }
    double sdiag =
        sqrt((shi[0] - slo[0]) * (shi[0] - slo[0]) + (shi[1] - slo[1]) * (shi[1] - slo[1]) +
             (shi[2] - slo[2]) * (shi[2] - slo[2]));
    double hcell[HZ_NH];
    for (int j = 0; j < HZ_NH; j++)
      hcell[j] = sdiag / pow(2.0, (double)j);
    double *touch = calloc((size_t)L.nlev * HZ_NH, sizeof *touch);
    int32_t *lcut = malloc((size_t)L.np * sizeof *lcut);
    int64_t *nbv = malloc((size_t)L.nlev * sizeof *nbv);
    int64_t *nbvs = malloc((size_t)L.nlev * sizeof *nbvs);
    int64_t *nloop = malloc((size_t)L.nlev * sizeof *nloop);
    int64_t *nouter = malloc((size_t)L.nlev * sizeof *nouter);
    int64_t *nhole = malloc((size_t)L.nlev * sizeof *nhole);
    if (lcut != NULL && nbv != NULL && nbvs != NULL && nloop != NULL && nouter != NULL &&
        nhole != NULL && touch != NULL) {
      for (int lev = 0; lev < L.nlev; lev++) {
        nbv[lev] = 0;
        nbvs[lev] = 0;
        nloop[lev] = 0;
        nouter[lev] = 0;
        nhole[lev] = 0;
        for (int32_t k = 0; k < L.np; k++)
          lcut[k] = L.lab[(size_t)lev * (size_t)L.np + (size_t)k];
        hz_pseglist so;
        if (hz_lod_seglist(&L, &m, &sg, lcut, &so) != 0) continue;
        hz_polyset pc;
        if (hz_poly_build(&pc, &m, &so) == 0) {
          nbv[lev] = pc.nbv;
          /* КРАЙ УРОВНЯ УПРОЩАЕТСЯ ДОПУСКОМ ЭТОГО УРОВНЯ (§324). Прежний замер
           * (§323) считал вершины края БЕЗ упрощения на всех уровнях сразу и
           * получил `Σ V/V₀ = 6.4…7.7×` — но это цена не пирамиды, а лестницы,
           * которая огрубляет ГРУППИРОВКУ и не огрубляет КРАЙ: на городе
           * уровень 8 нёс 4280 вершин на элемент. Допуск берётся ТОТ ЖЕ, что
           * гейтит поверхность уровня (`δ_L`), а не новый порог: у грубого
           * уровня обе ошибки — смещения поверхности и спрямления края —
           * ограничены одним числом. */
          /* ПЕТЛИ СЧИТАЮТСЯ ОТДЕЛЬНО, И ЭТО НЕ ЛЮБОПЫТСТВО (§324). Спрямление
           * края не может убрать петлю целиком: у всякой петли остаётся не
           * менее трёх вершин при любом допуске. Значит `V ≥ 3 · (число
           * НЕСВЯЗНЫХ кусков)`, и если лестница число кусков не сокращает, а
           * лишь переклеивает на них ярлык группы, край упрётся в этот пол. */
          /* ВНЕШНИЙ КОНТУР ОТДЕЛЬНО ОТ ДЫРКИ (§326). Петля есть граничный ЦИКЛ,
           * и это ЛИБО внешний обвод куска, ЛИБО дырка в нём; §324 назвал петли
           * «несвязными кусками» и ошибся — на нулевом уровне зала их `3` на
           * элемент, где элемент связен по построению. Различает ЗНАК площади в
           * местной раме: у внешнего обвода и у дырки обходы противоположны.
           * Число внешних обводов и есть число СВЯЗНЫХ КУСКОВ элемента. */
          for (int32_t q = 0; q < pc.np; q++) {
            const hz_poly *P = &pc.p[q];
            nloop[lev] += P->nloop;
            for (int32_t l = 0; l < P->nloop; l++) {
              int32_t a = pc.loop[P->l0 + l], b = pc.loop[P->l0 + l + 1];
              double s2 = 0.0;
              for (int32_t v = a; v < b; v++) {
                int32_t w = (v + 1 < b) ? v + 1 : a;
                s2 += pc.bv[2 * (size_t)v] * pc.bv[2 * (size_t)w + 1] -
                      pc.bv[2 * (size_t)w] * pc.bv[2 * (size_t)v + 1];
              }
              if (s2 > 0.0)
                nouter[lev]++;
              else if (s2 < 0.0)
                nhole[lev]++;
            }
          }
          /* КАСАНИЯ (ЯЧЕЙКА, КАНДИДАТ) ПО РАВНОМЕРНОЙ СЕТКЕ (§327, А546).
           * Элемент с габаритом `[lo,hi]` при ячейке `h` встречает ровно
           * `∏(⌈hi/h⌉ − ⌊lo/h⌋)` ячеек — считается точно, без дерева и без
           * фронта: правило дробления, пол прохода и камера к вопросу не
           * относятся, а отношение `T(L,h)/T(0,h)` при одной `h` от общего
           * множителя сетки не зависит. Габарит элемента берётся по его
           * ТРЕУГОЛЬНИКАМ — то же множество, по которому считан `dmax`. */
          for (int32_t q = 0; q < pc.np; q++) {
            const hz_poly *P = &pc.p[q];
            double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
            for (int32_t t = 0; t < P->ntri; t++) {
              int32_t ti = pc.tri[P->t0 + t];
              for (int cc = 0; cc < 3; cc++) {
                const double *vv = m.v + (size_t)m.f[3 * (size_t)ti + (size_t)cc] * 3;
                for (int c = 0; c < 3; c++) {
                  if (vv[c] < lo[c]) lo[c] = vv[c];
                  if (vv[c] > hi[c]) hi[c] = vv[c];
                }
              }
            }
            if (!(lo[0] <= hi[0])) continue;
            for (int j = 0; j < HZ_NH; j++) {
              double h = hcell[j];
              double cells = 1.0;
              for (int c = 0; c < 3; c++)
                cells *= ceil(hi[c] / h) - floor(lo[c] / h);
              touch[(size_t)lev * HZ_NH + (size_t)j] += cells;
            }
          }
          hz_edgestat es;
          if (hz_edge_simplify(&pc, L.delta0 * pow(ladb, lev), HZ_EDGE_SHARED, &es) == 0)
            nbvs[lev] = es.nbv_out;
          hz_poly_free(&pc);
        }
        hz_seg_free(&so);
      }
      if (nbv[0] > 0) {
        printf("      доли V_L/V_0 (ВЕРШИН КРАЯ, край НЕ огрублён); V_0 = %lld:",
               (long long)nbv[0]);
        for (int lev = 0; lev < L.nlev; lev++)
          printf(" %.3f", (double)nbv[lev] / (double)nbv[0]);
        printf("\n      вершин края НА ЭЛЕМЕНТ:");
        for (int lev = 0; lev < L.nlev; lev++)
          printf(" %.0f", (double)nbv[lev] / (double)(ncnt[lev] ? ncnt[lev] : 1));
        printf("\n");
      }
      if (nbvs[0] > 0) {
        printf("      доли Vs_L/Vs_0 (край ОГРУБЛЁН допуском уровня); Vs_0 = %lld:",
               (long long)nbvs[0]);
        for (int lev = 0; lev < L.nlev; lev++)
          printf(" %.3f", (double)nbvs[lev] / (double)nbvs[0]);
        printf("\n      вершин края НА ЭЛЕМЕНТ после огрубления:");
        for (int lev = 0; lev < L.nlev; lev++)
          printf(" %.0f", (double)nbvs[lev] / (double)(ncnt[lev] ? ncnt[lev] : 1));
        printf("\n");
      }
      if (touch[0] > 0.0) {
        printf("      КАСАНИЙ (ячейка, кандидат) ПО РАВНОМЕРНОЙ СЕТКЕ, §327:\n");
        printf("         ячейка, м:   ");
        for (int j = 0; j < HZ_NH; j++)
          printf(" %9.3g", hcell[j]);
        printf("\n         УРОВЕНЬ 0, штук:");
        for (int j = 0; j < HZ_NH; j++)
          printf(" %9.3g", touch[j]);
        printf("\n");
        for (int lev = 1; lev < L.nlev; lev++) {
          printf("         ур.%d, T_L/T_0: ", lev);
          for (int j = 0; j < HZ_NH; j++)
            printf(" %9.3f",
                   (touch[j] > 0.0) ? touch[(size_t)lev * HZ_NH + (size_t)j] / touch[j] : 0.0);
          printf("\n");
        }
      }
      if (nloop[0] > 0) {
        printf("      ПЕТЕЛЬ края (граничных циклов) всего:");
        for (int lev = 0; lev < L.nlev; lev++)
          printf(" %lld", (long long)nloop[lev]);
        printf("\n      из них ВНЕШНИХ ОБВОДОВ (= СВЯЗНЫХ КУСКОВ):");
        for (int lev = 0; lev < L.nlev; lev++)
          printf(" %lld", (long long)nouter[lev]);
        printf("\n      СВЯЗНЫХ КУСКОВ НА ЭЛЕМЕНТ:");
        for (int lev = 0; lev < L.nlev; lev++)
          printf(" %.1f", (double)nouter[lev] / (double)(ncnt[lev] ? ncnt[lev] : 1));
        printf("\n      петель НА ЭЛЕМЕНТ:");
        for (int lev = 0; lev < L.nlev; lev++)
          printf(" %.0f", (double)nloop[lev] / (double)(ncnt[lev] ? ncnt[lev] : 1));
        printf("\n      вершин на ПЕТЛЮ после огрубления (пол — три):");
        for (int lev = 0; lev < L.nlev; lev++)
          printf(" %.1f", (double)nbvs[lev] / (double)(nloop[lev] ? nloop[lev] : 1));
        printf("\n");
      }
    }
    free(nloop);
    free(touch);
    free(nouter);
    free(nhole);
    for (int step = 1; step <= 3; step++) {
      double sum = 0.0, sump = 0.0, sumv = 0.0, sumvs = 0.0;
      int nkeep = 0;
      for (int lev = 0; lev < L.nlev; lev += step) {
        sum += (double)ncnt[lev] / (double)ncnt[0];
        if (nper[0] > 0.0) sump += nper[lev] / nper[0];
        if (nbv != NULL && nbv[0] > 0) sumv += (double)nbv[lev] / (double)nbv[0];
        if (nbvs != NULL && nbvs[0] > 0) sumvs += (double)nbvs[lev] / (double)nbvs[0];
        nkeep++;
      }
      /* самый грубый уровень обязателен во всякой подлестнице */
      if ((L.nlev - 1) % step != 0) {
        sum += (double)ncnt[L.nlev - 1] / (double)ncnt[0];
        if (nper[0] > 0.0) sump += nper[L.nlev - 1] / nper[0];
        if (nbv != NULL && nbv[0] > 0) sumv += (double)nbv[L.nlev - 1] / (double)nbv[0];
        if (nbvs != NULL && nbvs[0] > 0) sumvs += (double)nbvs[L.nlev - 1] / (double)nbvs[0];
        nkeep++;
      }
      printf("      шаг по уровням %d (основание %.0f): уровней %d, Σ N_L/N_0 = %.2f×, "
             "Σ P_L/P_0 = %.2f×, Σ V_L/V_0 = %.2f× (край не огрублён), "
             "Σ Vs_L/Vs_0 = %.2f× (ОГРУБЛЁН)\n",
             step, pow(ladb, step), nkeep, sum, sump, sumv, sumvs);
    }
    free(lcut);
    free(nbv);
    free(nbvs);
    if (L.nlev > 1 && ncnt[L.nlev - 1] > 0)
      printf("      среднее геометрическое сокращение %.2f× на уровень (от N_0 к самому "
             "грубому)\n",
             pow((double)ncnt[0] / (double)ncnt[L.nlev - 1], 1.0 / (double)(L.nlev - 1)));
    fflush(stdout);
  }
  free(ncnt);
  free(nper);

  /* --- ДВУГРАННЫЙ УГОЛ ПАР-КАНДИДАТОВ (§59) ---
   * Печатается характеристика СЦЕНЫ, а не алгоритма: где стоят изломы и с какой
   * стороны. `180°` — плоское; левее выпуклое (внешний угол), правее вогнутое. */
  if (L.ndih_conv + L.ndih_conc + L.ndih_flat > 0) {
    printf("   двугранный угол кандидатов: плоских (180°) %lld, выпуклых (<180°) %lld, "
           "вогнутых (>180°) %lld\n",
           (long long)L.ndih_flat, (long long)L.ndih_conv, (long long)L.ndih_conc);
    printf("      корзины по 30°:");
    for (int q = 0; q < 12; q++)
      printf(" %d-%d:%lld", q * 30, q * 30 + 30, (long long)L.dih_hist[q]);
    printf("\n");
  }

  /* --- ВЛОЖЕННОСТЬ, целочисленно (А132) --- */
  {
    int bad = 0;
    for (int lev = 1; lev < L.nlev && !bad; lev++) {
      const int32_t *p = L.lab + (size_t)(lev - 1) * (size_t)L.np;
      const int32_t *c = L.lab + (size_t)lev * (size_t)L.np;
      /* Каждый узел предыдущего уровня обязан целиком лежать в ОДНОМ узле
       * следующего: иначе разбиение не вложено. */
      int32_t *map = malloc((size_t)L.nnd * sizeof *map);
      if (map == NULL) break;
      for (int32_t i = 0; i < L.nnd; i++)
        map[i] = -1;
      for (int32_t k = 0; k < L.np; k++) {
        if (map[p[k]] < 0)
          map[p[k]] = c[k];
        else if (map[p[k]] != c[k])
          bad = 1;
      }
      free(map);
    }
    printf("   ВЛОЖЕННОСТЬ (по меткам, А132): %s\n",
           bad ? "НАРУШЕНА — разбиения уровней не вложены" : "держится на всех уровнях");
  }

  /* ТОЛЬКО ТАБЛИЦА УРОВНЕЙ (§321): всё ниже — срез, лучевая метрика и негативный
   * контроль СРЕЗА, то есть другой вопрос и другая цена. */
  if (levels) {
    hz_lod_free(&L);
    hz_poly_free(&psf);
    hz_seg_free(&sg);
    hz_obj_free(&m);
    return 0;
  }

  /* --- срез и однородные уровни --- */
  int32_t *cut = malloc((size_t)L.np * sizeof *cut);
  if (cut == NULL) return 1;
  const double *eye = eyes;
  int32_t ncut = hz_lod_cut(&L, eye, eps, nosin, cut);
  printf("   СРЕЗ при камере §2 (ε = %.3e рад): элементов %d\n", eps, ncut);
  {
    hz_pseglist so;
    if (hz_lod_seglist(&L, &m, &sg, cut, &so) == 0) {
      hz_polyset pc;
      if (hz_poly_build(&pc, &m, &so) == 0) {
        if (vfit) {
          hz_vfitstat vs;
          hz_poly_vfit(&pc, dseg * pow(2.0, L.nlev - 1), &vs);
          printf("      подгонка вершин: передвинуто %lld, отказов по смещению %lld\n",
                 (long long)vs.nvert_moved, (long long)vs.nvert_far);
        }
        metric(&psf, &pc, &cam, eps_px, "СРЕЗ по камере", pc.np);
        hz_poly_free(&pc);
      }
      hz_seg_free(&so);
    }
  }
  /* Однородные уровни для сравнения: тот же путь, но без среза. */
  if (uniform)
    for (int lev = 1; lev < L.nlev; lev++) {
      for (int32_t k = 0; k < L.np; k++)
        cut[k] = L.lab[(size_t)lev * (size_t)L.np + (size_t)k];
      hz_pseglist so;
      if (hz_lod_seglist(&L, &m, &sg, cut, &so) != 0) continue;
      hz_polyset pc;
      if (hz_poly_build(&pc, &m, &so) == 0) {
        char tag[64];
        snprintf(tag, sizeof tag, "однородный уровень %d", lev);
        metric(&psf, &pc, &cam, eps_px, tag, pc.np);
        hz_poly_free(&pc);
      }
      hz_seg_free(&so);
    }
  /* --- КРИВАЯ ПО ε: ЦЕНА ПРОТИВ ОШИБКИ ДЛЯ ЭТОЙ ТАКТИКИ (§69) --- */
  if (curve) {
    const double mul[6] = {0.25, 0.5, 1.0, 2.0, 4.0, 8.0};
    for (int i = 0; i < 6; i++) {
      hz_lod_cut(&L, eye, eps * mul[i], nosin, cut);
      hz_pseglist so;
      if (hz_lod_seglist(&L, &m, &sg, cut, &so) != 0) continue;
      hz_polyset pc;
      if (hz_poly_build(&pc, &m, &so) == 0) {
        /* ПОДГОНКА ВЕРШИН ВНУТРИ КРИВОЙ (§73). Прежде она стояла только в ветви
         * «срез», и прогон с `vfit` давал кривую, ТОЖДЕСТВЕННУЮ прогону без него —
         * до последней цифры. Молчаливо недействительный замер того же класса, что
         * и `bycount bands` в zsh: конфигурация выглядит новой и повторяет старую. */
        if (vfit) {
          hz_vfitstat vs;
          hz_poly_vfit(&pc, (vfitlim > 0.0) ? vfitlim : dseg * pow(2.0, L.nlev - 1), &vs);
          /* ПОЛНЫЙ ЗАМЕР ПОДГОНКИ, а не два счётчика (§73). Вопрос пользователя:
           * как оптимизация может ухудшать, если у неё есть проверка улучшения?
           * Ответ обязан быть числом: `t*` — невязка ЕЁ цели (вершина на всех
           * своих плоскостях), и она обязана быть мала, тогда как лучевая метрика
           * при этом портится. Расхождение двух чисел и есть механизм: цель одна,
           * мерится другое, а между ними — проекция в двумерный край. */
          printf("      подгонка: вершин %lld, двинуто %lld, оставлено %lld, отказов "
                 "(далеко %lld, без улучшения %lld); t* макс %.4f м, смещение макс %.4f м, "
                 "расхождение площади %.3e\n",
                 (long long)vs.nvert, (long long)vs.nvert_moved, (long long)vs.nvert_fixed,
                 (long long)vs.nvert_far, (long long)vs.nvert_noimp, vs.tmax, vs.dmax_move,
                 vs.dmax_area);
        }
        char tag[64];
        snprintf(tag, sizeof tag, "КРИВАЯ: срез ε × %.2f", mul[i]);
        metric(&psf, &pc, &cam, eps_px, tag, pc.np);
        hz_poly_free(&pc);
      }
      hz_seg_free(&so);
    }
  }

  /* --- НЕГАТИВНЫЙ КОНТРОЛЬ: ТЕ ЖЕ УРОВНИ, РОЗДАННЫЕ ДРУГИМ ПОЛИГОНАМ (А146) ---
   *
   * Случайный ПРЕДОК не годится: он меняет цену, и провал вышел бы по причине, к
   * критерию среза отношения не имеющей. Перестановка сохраняет мультимножество
   * выбранных уровней ТОЧНО и рвёт только связь «грубее там, где дальше».
   * Перестановка ДЕТЕРМИНИРОВАННАЯ — хешем от номера: замер обязан повторяться.
   * Печатается и гистограмма уровней (обязана совпасть с гистограммой среза), и
   * итоговое число элементов — оно разойтись МОЖЕТ, и это не порок контроля, а
   * его цена: два соседних полигона на разных уровнях дают два узла там, где срез
   * давал один (А155). */
  if (uniform) {
    int32_t *lv = malloc((size_t)L.np * sizeof *lv);
    int64_t hcut[64], hperm[64];
    for (int i = 0; i < 64; i++)
      hcut[i] = hperm[i] = 0;
    if (lv != NULL) {
      /* СРЕЗ ПЕРЕСЧИТЫВАЕТСЯ ЗАНОВО: цикл по однородным уровням выше пишет в тот
       * же `cut`, и без этого контроль брал бы уровни ПОСЛЕДНЕГО однородного
       * уровня вместо среза. Поймано печатью гистограммы (А155) в первом же
       * прогоне: она показала «всё на уровне 6» там, где срез сидит на нулевом. */
      hz_lod_cut(&L, eye, eps, nosin, cut);
      for (int32_t k = 0; k < L.np; k++) {
        lv[k] = L.nd[cut[k]].level;
        if (lv[k] >= 0 && lv[k] < 64) hcut[lv[k]]++;
      }
      for (int32_t k = 0; k < L.np; k++) {
        uint64_t x = (uint64_t)(uint32_t)k * 0x9e3779b97f4a7c15ull;
        x ^= x >> 30;
        x *= 0xbf58476d1ce4e5b9ull;
        x ^= x >> 27;
        int32_t j = (int32_t)(x % (uint64_t)(uint32_t)L.np);
        int32_t t = lv[k];
        lv[k] = lv[j];
        lv[j] = t;
      }
      for (int32_t k = 0; k < L.np; k++) {
        int32_t l = lv[k];
        if (l < 0) l = 0;
        if (l >= L.nlev) l = L.nlev - 1;
        if (l < 64) hperm[l]++;
        cut[k] = L.lab[(size_t)l * (size_t)L.np + (size_t)k];
      }
      int same = 1;
      for (int i = 0; i < 64; i++)
        if (hcut[i] != hperm[i]) same = 0;
      printf("   КОНТРОЛЬ (перестановка уровней): гистограмма уровней %s;",
             same ? "СОВПАЛА ТОЧНО — цена сохранена" : "РАЗОШЛАСЬ — контроль недоказателен");
      for (int i = 0; i < L.nlev && i < 64; i++)
        printf(" %d:%lld", i, (long long)hcut[i]);
      printf("\n");
      hz_pseglist so;
      if (hz_lod_seglist(&L, &m, &sg, cut, &so) == 0) {
        hz_polyset pc;
        if (hz_poly_build(&pc, &m, &so) == 0) {
          metric(&psf, &pc, &cam, eps_px, "КОНТРОЛЬ: уровни перемешаны", pc.np);
          hz_poly_free(&pc);
        }
        hz_seg_free(&so);
      }
      free(lv);
    }
  }

  /* --- СВИП ПО `ε`: ПОДПИСЬ АРТЕФАКТА (а) ---
   * Число элементов среза ОБЯЗАНО меняться от `ε` — от той самой величины, которой
   * его меняют. Метрика здесь НЕ считается (А157): вопрос свипа не про качество, а
   * про то, что критерий вообще работает, и ответ на него стоит секунды.
   * Рядом печатается доля полигонов, СЕВШИХ НА ПОСЛЕДНИЙ УРОВЕНЬ (А156): пока она
   * мала, свип мерит КРИТЕРИЙ; когда велика — ПОТОЛОК лестницы, и это разные вещи.
   * Оговорка: лестница откалибрована под своё `ε`, здесь меняется только срез. */
  {
    printf("   СВИП ПО ε (подпись артефакта; лестница откалибрована под ε = %.3e):\n", L.eps);
    const double mul[4] = {0.5, 1.0, 2.0, 4.0};
    for (int i = 0; i < 4; i++) {
      int32_t nc = hz_lod_cut(&L, eye, eps * mul[i], nosin, cut);
      int64_t top = 0;
      for (int32_t k = 0; k < L.np; k++)
        if (L.nd[cut[k]].level == L.nlev - 1) top++;
      printf("      ε × %.1f = %.3e: узлов среза %7d, на последнем уровне %.2f%% полигонов\n",
             mul[i], eps * mul[i], nc, 100.0 * (double)top / (double)(L.np ? L.np : 1));
    }
  }
  printf("   время: постройка %.1f с, всего %.1f с\n", t1 - t0, now_s() - t0);
  free(cut);
  hz_lod_free(&L);
  hz_poly_free(&psf);
  hz_seg_free(&sg);
  hz_obj_free(&m);
  return 0;
}
