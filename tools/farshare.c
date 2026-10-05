/* farshare — §931 Ш1 (А1664/А1672/А1681/А1685): офлайн-таблица колец
 * детальности по ЗРЕЛОМУ E-дампу эталона (tmp/E_928A.bin, it=8).
 *
 *   ближний кусок  — центроид ближе R от глаза (кусковой мир §923, битово);
 *   дальний кусок  — центроид дальше R (мир reps, недобор терпим, К2).
 *
 * Для каждого R сетки печатает ДОЛЕ ЭНЕРГИИ (ΣEd·A):
 *   D_near(R) — ближняя доля, растёт по R; санити А1685: D_near(0)=0,
 *     D_near(∞)=1 (тождественно — ловит перепутанные индексы/единицы);
 *   D_far(R) = 1 − D_near(R) — доля дальних, гейт К2: дрейф ≤ (1−k_far)·D_far,
 *     пригодность храма А1684: ∃R с D_far ≥ 10%.
 *
 * Маска int8[nt] (А1672/А1677(а)): 1 = ближний (центроид ≤ R), 0 = дальний,
 * в TRI-пространстве OBJ — та же нумерация, что у E-файла (§902).
 *
 * use: farshare E=файл obj=файл eye=X,Y,Z R=2,3,4,5,6,8 [mask=шаблон{R}]
 *   все ключи обязательны, кроме mask=; в шаблоне маски {R} заменяется на
 *   значение R (без {R} маска пишется только для последнего R сетки).
 * Поле Ed от глаза не зависит (это не кадр); глаз нужен только для маски и
 * долей. Глаз по умолчанию НЕ угадывается: eye= обязателен (А1672).
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "scene_obj.h"

/* Порог пригодности полигонa А1684: доля дальнего поля, ниже которой
 * кольцам нечего мерить (k_far по малому множеству шумит, А1674). */
#define FS_FARMIN 0.10

static int parse3(const char *s, double o[3]) {
  return sscanf(s, "%lf,%lf,%lf", &o[0], &o[1], &o[2]) == 3 ? 0 : 1;
}

/* Шаблон пути: {R} → значение R (%g). Без {R} — копия как есть. */
static int mkpath(char *dst, size_t cap, const char *tpl, double r) {
  const char *p = strstr(tpl, "{R}");
  int n;
  if (p == NULL)
    n = snprintf(dst, cap, "%s", tpl);
  else
    n = snprintf(dst, cap, "%.*s%g%s", (int)(p - tpl), tpl, r, p + 3);
  return (n > 0 && (size_t)n < cap) ? 0 : 1;
}

int main(int argc, char **argv) {
  const char *epath = NULL, *opath = NULL, *mtpl = NULL, *rspec = NULL, *e2path = NULL;
  double *E2 = NULL; /* §931 A2: опциональный мир колец для К1-L1/k_far */
  double eye[3] = {0, 0, 0};
  double rs[64];
  int nr = 0, i;
  hz_objmesh m;
  double *Ed = NULL, *cent = NULL, *area = NULL, tot = 0, dmaxc = 0;
  FILE *f;
  for (i = 1; i < argc; i++) {
    if (strncmp(argv[i], "E=", 2) == 0)
      epath = argv[i] + 2;
    else if (strncmp(argv[i], "obj=", 4) == 0)
      opath = argv[i] + 4;
    else if (strncmp(argv[i], "eye=", 4) == 0) {
      if (parse3(argv[i] + 4, eye)) goto usage;
    } else if (strncmp(argv[i], "R=", 2) == 0)
      rspec = argv[i] + 2;
    else if (strncmp(argv[i], "mask=", 5) == 0)
      mtpl = argv[i] + 5;
    else if (strncmp(argv[i], "E2=", 3) == 0)
      e2path = argv[i] + 3; /* §931 A2: мир колец против E=REF — К1-L1/k_far */
    else
      goto usage;
  }
  if (epath == NULL || opath == NULL || rspec == NULL) goto usage;
  { /* R-сетка: запятые, без пробелов */
    char *buf = strdup(rspec), *tok;
    if (buf == NULL) return 2;
    for (tok = strtok(buf, ","); tok != NULL && nr < 64; tok = strtok(NULL, ","))
      rs[nr++] = atof(tok);
    free(buf);
    if (nr <= 0) goto usage;
  }
  if (hz_obj_load(&m, opath, 1.0) != 0) {
    fprintf(stderr, "farshare: obj не читается: %s\n", opath);
    return 2;
  }
  /* КОНТРОЛЬ ВХОДА А1681: E-файл обязан иметь РОВНО nt double — fail closed. */
  f = fopen(epath, "rb");
  if (f == NULL) {
    fprintf(stderr, "farshare: E не открылся: %s\n", epath);
    return 2;
  }
  if (fseek(f, 0, SEEK_END) != 0) return 2;
  {
    long sz = ftell(f);
    if (sz != (long)m.nt * (long)sizeof(double)) {
      fprintf(stderr, "farshare: nt E-файла (%ld) != nt OBJ (%d) — ОТКАЗ\n",
              sz / (long)sizeof(double), (int)m.nt);
      return 2;
    }
    if (fseek(f, 0, SEEK_SET) != 0) return 2;
  }
  Ed = (double *)malloc((size_t)m.nt * sizeof *Ed);
  cent = (double *)malloc(3u * (size_t)m.nt * sizeof *cent);
  area = (double *)malloc((size_t)m.nt * sizeof *area);
  if (Ed == NULL || cent == NULL || area == NULL) return 2;
  if (fread(Ed, sizeof(double), (size_t)m.nt, f) != (size_t)m.nt) return 2;
  fclose(f);
  if (e2path != NULL) { /* §931 A2: мир колец, тот же контроль nt (А1681) */
    FILE *f2 = fopen(e2path, "rb");
    long sz2;
    if (f2 == NULL) {
      fprintf(stderr, "farshare: E2 не открылся: %s\n", e2path);
      return 2;
    }
    if (fseek(f2, 0, SEEK_END) != 0) return 2;
    sz2 = ftell(f2);
    if (sz2 != (long)m.nt * (long)sizeof(double)) {
      fprintf(stderr, "farshare: nt E2-файла (%ld) != nt OBJ (%d) — ОТКАЗ\n",
              sz2 / (long)sizeof(double), (int)m.nt);
      return 2;
    }
    if (fseek(f2, 0, SEEK_SET) != 0) return 2;
    E2 = (double *)malloc((size_t)m.nt * sizeof *E2);
    if (E2 == NULL) return 2;
    if (fread(E2, sizeof(double), (size_t)m.nt, f2) != (size_t)m.nt) return 2;
    fclose(f2);
  }
  for (i = 0; i < m.nt; i++) {
    const double *a = m.v + 3 * (size_t)m.f[3 * i], *b = m.v + 3 * (size_t)m.f[3 * i + 1],
                 *c = m.v + 3 * (size_t)m.f[3 * i + 2];
    double u[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]},
           v[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]}, cx = u[1] * v[2] - u[2] * v[1],
           cy = u[2] * v[0] - u[0] * v[2], cz = u[0] * v[1] - u[1] * v[0];
    int k;
    area[i] = 0.5 * sqrt(cx * cx + cy * cy + cz * cz);
    for (k = 0; k < 3; k++)
      cent[3 * i + k] = (a[k] + b[k] + c[k]) / 3.0 - eye[k];
    {
      double d2 = cent[3 * i] * cent[3 * i] + cent[3 * i + 1] * cent[3 * i + 1] +
                  cent[3 * i + 2] * cent[3 * i + 2];
      if (sqrt(d2) > dmaxc) dmaxc = sqrt(d2);
    }
    tot += Ed[i] * area[i];
  }
  printf("farshare: nt=%d eye=(%g,%g,%g) max|c−eye|=%.3f ΣEd·A=%.6g\n", (int)m.nt, eye[0], eye[1],
         eye[2], dmaxc, tot);
  /* САНИТИ А1685 (по ближней доле, до использования маски в гейтах).
   * Тождественность — целочисленными счётчиками, без float-сравнений:
   * R=0 «ближний» ⇔ d≤0 (ждём ноль таких), R=∞ — все (ждём nt). */
  {
    int32_t cn0 = 0, cninf = 0;
    for (i = 0; i < m.nt; i++) {
      double d = sqrt(cent[3 * i] * cent[3 * i] + cent[3 * i + 1] * cent[3 * i + 1] +
                      cent[3 * i + 2] * cent[3 * i + 2]);
      if (d <= 0.0) cn0++;
      if (d <= 1e9 * (1.0 + dmaxc)) cninf++;
    }
    printf("САНИТИ А1685: near(0)=%d (ждём 0)  near(∞)=%d (ждём %d)  %s\n", (int)cn0, (int)cninf,
           (int)m.nt, (cn0 == 0 && cninf == m.nt) ? "PASS" : "FAIL — маску в гейты НЕ ДАВАТЬ");
    if (cn0 != 0 || cninf != m.nt) return 3;
  }
  printf("%6s %12s %12s %10s %6s\n", "R", "D_far", "D_near", "near_tri", "годен");
  for (i = 0; i < nr; i++) {
    double r = rs[i], dnf = 0;
    int32_t tn = 0, t;
    int8_t *msk = NULL;
    for (t = 0; t < m.nt; t++) {
      double d = sqrt(cent[3 * (int64_t)t] * cent[3 * (int64_t)t] +
                      cent[3 * t + 1] * cent[3 * t + 1] + cent[3 * t + 2] * cent[3 * t + 2]);
      if (d > r)
        dnf += Ed[t] * area[t];
      else
        tn++;
    }
    printf("%6g %12.6f %12.6f %10d %6s\n", r, dnf / tot, 1.0 - dnf / tot, (int)tn,
           (dnf / tot >= FS_FARMIN) ? "ДА" : "нет");
    if (E2 != NULL) { /* §931 A2: гейты точки против REF (Ed) — К1-L1/k_far
                       * (А1667/А1669); ближнее множество = маска d<=R */
      double l1n = 0, sref_n = 0, sworld_n = 0, sref_f = 0, sworld_f = 0, sw_a = 0, ref_a = 0;
      for (t = 0; t < m.nt; t++) {
        double d = sqrt(cent[3 * (int64_t)t] * cent[3 * (int64_t)t] +
                        cent[3 * t + 1] * cent[3 * t + 1] + cent[3 * t + 2] * cent[3 * t + 2]);
        double wa = E2[t] * area[t], ra = Ed[t] * area[t];
        sw_a += wa;
        ref_a += ra;
        if (d <= r) {
          l1n += fabs(E2[t] - Ed[t]) * area[t];
          sref_n += ra;
          sworld_n += wa;
        } else {
          sref_f += ra;
          sworld_f += wa;
        }
      }
      printf("  К1: L1_near=%.6f (гейт<=0.01)  E_avg_near=%+.4f%% (гейт ±1%%)\n"
             "  К2: k_far=%.4f  D_far=%.4f  дрейф E_avg=%+.4f%%  граница(конс. k=0.21)=%.2f%%\n",
             l1n / sref_n, 100.0 * (sworld_n / sref_n - 1.0), sworld_f / sref_f, dnf / tot,
             100.0 * (sw_a / ref_a - 1.0), 100.0 * 0.79 * dnf / tot);
    }
    if (mtpl != NULL && (strstr(mtpl, "{R}") != NULL || i == nr - 1)) {
      char path[512];
      if (mkpath(path, sizeof path, mtpl, r) != 0) return 2;
      msk = (int8_t *)malloc((size_t)m.nt);
      if (msk == NULL) return 2;
      for (t = 0; t < m.nt; t++) {
        double d = sqrt(cent[3 * (int64_t)t] * cent[3 * (int64_t)t] +
                        cent[3 * t + 1] * cent[3 * t + 1] + cent[3 * t + 2] * cent[3 * t + 2]);
        msk[t] = (d <= r) ? 1 : 0;
      }
      f = fopen(path, "wb");
      if (f == NULL) return 2;
      if (fwrite(msk, 1, (size_t)m.nt, f) != (size_t)m.nt) return 2;
      fclose(f);
      free(msk);
      printf("  маска ближних: %s (int8[%d])\n", path, (int)m.nt);
    }
  }
  free(E2); /* §931 A2 */
  free(Ed);
  free(cent);
  free(area);
  hz_obj_free(&m);
  return 0;
usage:
  fprintf(stderr, "use: farshare E=файл obj=файл eye=X,Y,Z R=2,3,4 [mask=шаблон{R}]\n");
  return 1;
}
