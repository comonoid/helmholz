/* gentemple.c — §915-R5: синтетический ХРАМ-МИНИ (~1.2e5 треугольников).
 *
 * Колоннада 24×24 (12-гранные призмы, 4 сегмента высоты), пол сеткой
 * 100×100, кровля сеткой 100×100 (поднята), ступени базы. Всё — v/f
 * строки OBJ; материалы: default (mtl не пишем — pgather берёт kd=0.5
 * по умолчанию, как у синтетик gen_room).
 *
 * Синтаксис: gentemple OUT.obj [side=N grid=M]  (side=24 grid=100)
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv) {
  int side = 24, grid = 100;
  if (argc < 2) {
    fprintf(stderr, "gentemple OUT.obj [side=N grid=M]\n");
    return 2;
  }
  for (int i = 2; i < argc; i++) {
    if (strncmp(argv[i], "side=", 5) == 0)
      side = atoi(argv[i] + 5);
    else if (strncmp(argv[i], "grid=", 5) == 0)
      grid = atoi(argv[i] + 5);
  }
  if (side < 2 || grid < 4 || side > 64 || grid > 400) {
    fprintf(stderr, "gentemple: 2<=side<=64, 4<=grid<=400\n");
    return 2;
  }
  FILE *f = fopen(argv[1], "w");
  if (f == NULL) {
    fprintf(stderr, "gentemple: не открывается %s\n", argv[1]);
    return 2;
  }
  long nv = 0, nt = 0;
  const double S = 40.0;    /* габарит храма, м */
  const double colR = 0.35; /* радиус колонны */
  const double colH = 6.0;  /* высота колонны */
  const int SIDES = 12;     /* грани колонны */
  const int SEGS = 4;       /* сегменты высоты */
  fprintf(f, "# gentemple side=%d grid=%d\n", side, grid);

  /* пол: сетка grid×grid на y=0 */
  double q = S / grid;
  for (int j = 0; j <= grid; j++)
    for (int i = 0; i <= grid; i++)
      fprintf(f, "v %.6f 0.000000 %.6f\n", i * q - S / 2, j * q - S / 2), nv++;
  for (int j = 0; j < grid; j++)
    for (int i = 0; i < grid; i++) {
      long a = 1 + (long)j * (grid + 1) + i, b = a + 1, c = a + grid + 1, d = c + 1;
      fprintf(f, "f %ld %ld %ld\nf %ld %ld %ld\n", a, d, b, a, c, d);
      nt += 2;
    }

  /* кровля: та же сетка, поднята на colH+0.5, нормали вниз (порядок) */
  long roof_base = nv;
  for (int j = 0; j <= grid; j++)
    for (int i = 0; i <= grid; i++)
      fprintf(f, "v %.6f %.6f %.6f\n", i * q - S / 2, colH + 0.5, j * q - S / 2), nv++;
  for (int j = 0; j < grid; j++)
    for (int i = 0; i < grid; i++) {
      long a = 1 + roof_base + (long)j * (grid + 1) + i, b = a + 1, c = a + grid + 1, d = c + 1;
      fprintf(f, "f %ld %ld %ld\nf %ld %ld %ld\n", a, b, d, a, d, c);
      nt += 2;
    }

  /* колоннада side×side */
  double step = (S - 2.0) / (side - 1);
  for (int cj = 0; cj < side; cj++)
    for (int ci = 0; ci < side; ci++) {
      double cx = -S / 2 + 1.0 + ci * step, cz = -S / 2 + 1.0 + cj * step;
      long base = nv;
      for (int sg = 0; sg <= SEGS; sg++) {
        double y = colH * sg / SEGS;
        for (int a = 0; a < SIDES; a++) {
          double ang = 2.0 * 3.14159265358979323846 * a / SIDES;
          fprintf(f, "v %.6f %.6f %.6f\n", cx + colR * cos(ang), y, cz + colR * sin(ang));
          nv++;
        }
      }
      for (int sg = 0; sg < SEGS; sg++)
        for (int a = 0; a < SIDES; a++) {
          long v00 = 1 + base + (long)sg * SIDES + a;
          long v01 = 1 + base + (long)sg * SIDES + (a + 1) % SIDES;
          long v10 = v00 + SIDES, v11 = v01 + SIDES;
          fprintf(f, "f %ld %ld %ld\nf %ld %ld %ld\n", v00, v01, v11, v00, v11, v10);
          nt += 2;
        }
      /* крышка колонны (венчающая плита) — quad в кровлю не считаем:
       * колонна закрыта кровлей сверху и полом снизу; боковины только */
    }

  /* ступени базы: 3 кольца вокруг (quad-стены по периметру) */
  for (int st = 0; st < 3; st++) {
    double ext = S / 2 + 1.0 + st, y0 = 0.15 * st, y1 = 0.15 * (st + 1);
    long vb = nv;
    for (int k = 0; k < 4; k++) {
      /* углы прямоугольника, дважды (низ/верх) */
      double sx = (k == 0 || k == 3) ? -ext : ext;
      double sz = (k < 2) ? -ext : ext;
      fprintf(f, "v %.6f %.6f %.6f\n", sx, y0, sz);
      fprintf(f, "v %.6f %.6f %.6f\n", sx, y1, sz);
      nv += 2;
    }
    /* 4 стены по периметру: углы 0-1, 1-2, 2-3, 3-0 */
    static const int ed[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    for (int k = 0; k < 4; k++) {
      long a0 = 1 + vb + 2 * ed[k][0], a1 = 1 + vb + 2 * ed[k][1];
      long b0 = a0 + 1, b1 = a1 + 1;
      fprintf(f, "f %ld %ld %ld\nf %ld %ld %ld\n", a0, a1, b1, a0, b1, b0);
      nt += 2;
    }
  }

  fprintf(f, "# gentemple: %ld вершин, %ld треугольников\n", nv, nt);
  fclose(f);
  printf("gentemple: nv=%ld nt=%ld → %s (side=%d grid=%d)\n", nv, nt, argv[1], side, grid);
  return 0;
}
