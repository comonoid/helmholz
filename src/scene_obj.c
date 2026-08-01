/* Разбор OBJ/MTL. Разбор и оговорки — в `scene_obj.h`. */

#include "scene_obj.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --- лексика --------------------------------------------------------------- */

static int is_sp(char c) {
  return c == ' ' || c == '\t' || c == '\r';
}
static int is_eol(char c) {
  return c == '\n' || c == '\0';
}

static const char *skip_sp(const char *p) {
  while (is_sp(*p))
    p++;
  return p;
}

static const char *skip_line(const char *p) {
  while (!is_eol(*p))
    p++;
  return (*p == '\n') ? p + 1 : p;
}

static const char *skip_token(const char *p) {
  while (!is_sp(*p) && !is_eol(*p))
    p++;
  return p;
}

/* Сколько слов до конца строки. Нужно ровно затем, чтобы посчитать, сколько
 * треугольников даст грань, ДО выделения памяти. */
static int count_tokens(const char *p) {
  int k = 0;
  for (;;) {
    p = skip_sp(p);
    if (is_eol(*p)) return k;
    k++;
    p = skip_token(p);
  }
}

/* Целое со знаком. `*ok` = 0, если цифр не было вовсе. */
static const char *parse_int(const char *p, long *out, int *ok) {
  int neg = 0;
  long v = 0;
  int any = 0;
  if (*p == '-') {
    neg = 1;
    p++;
  } else if (*p == '+') {
    p++;
  }
  while (*p >= '0' && *p <= '9') {
    /* Переполнение здесь означает битый файл, а не большую сцену: индексов
     * больше 2^31 в OBJ не бывает. Обрываем, а не заворачиваем по модулю. */
    if (v > 200000000L) {
      *ok = 0;
      return p;
    }
    v = v * 10 + (*p - '0');
    p++;
    any = 1;
  }
  *out = neg ? -v : v;
  *ok = any;
  return p;
}

/* Один угол грани: `v`, `v/vt`, `v//vn`, `v/vt/vn`. Индексы 1-базовые, а
 * отрицательные отсчитываются от ЧИСЛА ВИДЕННЫХ К ЭТОМУ МОМЕНТУ (см. заголовок).
 * `*ivn` = −1, если нормали у угла нет. */
static const char *parse_fvert(const char *p, int32_t nv_seen, int32_t nvn_seen, int32_t *iv,
                               int32_t *ivn, int *ok) {
  long a = 0;
  int got = 0;
  *ivn = -1;
  p = parse_int(p, &a, &got);
  if (!got) {
    *ok = 0;
    return p;
  }
  long idx = (a > 0) ? a - 1 : (long)nv_seen + a;
  if (idx < 0 || idx >= (long)nv_seen) {
    *ok = 0;
    return p;
  }
  *iv = (int32_t)idx;
  if (*p == '/') {
    p++;
    if (*p != '/') {
      long vt = 0;
      int g2 = 0;
      p = parse_int(p, &vt, &g2); /* координата текстуры: до Ш8 не нужна */
      (void)vt;
      (void)g2;
    }
    if (*p == '/') {
      p++;
      long b = 0;
      int g3 = 0;
      p = parse_int(p, &b, &g3);
      if (g3) {
        long jdx = (b > 0) ? b - 1 : (long)nvn_seen + b;
        if (jdx >= 0 && jdx < (long)nvn_seen) *ivn = (int32_t)jdx;
      }
    }
  }
  *ok = 1;
  return p;
}

/* --- чтение файла целиком -------------------------------------------------- */

/* Возврат NULL — не открыть, не прочесть либо не хватило памяти. Буфер всегда
 * NUL-терминирован, и на этом стоит вся лексика выше: конца буфера как
 * отдельного условия в ней нет. */
static char *slurp(const char *path, size_t *len) {
  FILE *fp = fopen(path, "rb");
  if (fp == NULL) return NULL;
  if (fseek(fp, 0, SEEK_END) != 0) {
    fclose(fp);
    return NULL;
  }
  long sz = ftell(fp);
  if (sz < 0) {
    fclose(fp);
    return NULL;
  }
  if (fseek(fp, 0, SEEK_SET) != 0) {
    fclose(fp);
    return NULL;
  }
  char *buf = malloc((size_t)sz + 1);
  if (buf == NULL) {
    fclose(fp);
    return NULL;
  }
  size_t got = fread(buf, 1, (size_t)sz, fp);
  fclose(fp);
  buf[got] = '\0';
  *len = got;
  return buf;
}

/* --- материалы ------------------------------------------------------------- */

static int mtl_find(const hz_objmesh *m, const char *name, size_t n) {
  if (n >= HZ_OBJ_MTLNAME) n = HZ_OBJ_MTLNAME - 1;
  for (int32_t i = 0; i < m->nmtl; i++)
    if (strncmp(m->mtl[i].name, name, n) == 0 && m->mtl[i].name[n] == '\0') return (int)i;
  return -1;
}

/* Добавление материала. РАСТУЩИЙ массив здесь законен: материалов десятки, это
 * не индексируемый горячий буфер, о котором предупреждает CLAUDE.md. */
static int mtl_add(hz_objmesh *m, const char *name, size_t n, double kd) {
  int have = mtl_find(m, name, n);
  if (have >= 0) return have;
  hz_obj_mtl *nm = realloc(m->mtl, (size_t)(m->nmtl + 1) * sizeof *nm);
  if (nm == NULL) return -1;
  m->mtl = nm;
  if (n >= HZ_OBJ_MTLNAME) n = HZ_OBJ_MTLNAME - 1;
  memcpy(m->mtl[m->nmtl].name, name, n);
  m->mtl[m->nmtl].name[n] = '\0';
  m->mtl[m->nmtl].kd = kd;
  for (int c = 0; c < 3; c++) {
    m->mtl[m->nmtl].kd3[c] = kd;
    m->mtl[m->nmtl].ks3[c] = 0.0;
  }
  m->mtl[m->nmtl].ns = 0.0;
  m->mtl[m->nmtl].ior = 1.0;
  m->mtl[m->nmtl].alpha = 1.0;
  m->mtl[m->nmtl].flat = 0;
  m->nmtl++;
  return (int)(m->nmtl - 1);
}

/* .mtl рядом с .obj. Отсутствие файла — НЕ ошибка: сцена без материалов
 * считается с альбедо по умолчанию, и это лучше, чем отказ грузить геометрию. */
static void mtl_load(hz_objmesh *m, const char *objpath, const char *name, size_t nlen) {
  char path[1024];
  const char *slash = strrchr(objpath, '/');
  size_t dirlen = (slash != NULL) ? (size_t)(slash - objpath) + 1 : 0;
  if (dirlen + nlen + 1 > sizeof path) return;
  memcpy(path, objpath, dirlen);
  memcpy(path + dirlen, name, nlen);
  path[dirlen + nlen] = '\0';

  size_t len = 0;
  char *buf = slurp(path, &len);
  if (buf == NULL) return;

  const char *p = buf;
  int cur = -1;
  while (*p != '\0') {
    p = skip_sp(p);
    if (strncmp(p, "newmtl", 6) == 0 && is_sp(p[6])) {
      const char *q = skip_sp(p + 6);
      const char *e = skip_token(q);
      cur = mtl_add(m, q, (size_t)(e - q), 0.5);
    } else if (p[0] == 'K' && p[1] == 'd' && is_sp(p[2]) && cur >= 0) {
      char *e = NULL;
      double c0 = strtod(p + 2, &e);
      double c1 = strtod(e, &e);
      double c2 = strtod(e, &e);
      /* ОДИН КАНАЛ — СРЕДНЕЕ, а не яркость: сохраняется доля энергии (заголовок).
       * Три канала хранятся рядом и нужны только картинке (§78). */
      m->mtl[cur].kd = (c0 + c1 + c2) / 3.0;
      m->mtl[cur].kd3[0] = c0;
      m->mtl[cur].kd3[1] = c1;
      m->mtl[cur].kd3[2] = c2;
    } else if (p[0] == 0x4B && p[1] == 0x73 && is_sp(p[2]) && cur >= 0) {
      char *e = NULL;
      m->mtl[cur].ks3[0] = strtod(p + 2, &e);
      m->mtl[cur].ks3[1] = strtod(e, &e);
      m->mtl[cur].ks3[2] = strtod(e, &e);
    } else if (p[0] == 0x4E && p[1] == 0x73 && is_sp(p[2]) && cur >= 0) {
      m->mtl[cur].ns = strtod(p + 2, NULL);
    } else if (p[0] == 0x4E && p[1] == 0x69 && is_sp(p[2]) && cur >= 0) {
      m->mtl[cur].ior = strtod(p + 2, NULL);
    } else if (p[0] == 0x64 && is_sp(p[1]) && cur >= 0) {
      m->mtl[cur].alpha = strtod(p + 1, NULL);
    } else if (strncmp(p, "hz_flat", 7) == 0 && is_sp(p[7]) && cur >= 0) {
      m->mtl[cur].flat = (strtod(p + 7, NULL) > 0.5);
    }
    p = skip_line(p);
  }
  free(buf);
}

/* --- загрузка -------------------------------------------------------------- */

void hz_obj_free(hz_objmesh *m) {
  free(m->v);
  free(m->vn);
  free(m->f);
  free(m->fn);
  free(m->fm);
  free(m->mtl);
  memset(m, 0, sizeof *m);
}

void hz_obj_tri(const hz_objmesh *m, int32_t t, double p[3][3]) {
  for (int i = 0; i < 3; i++) {
    int32_t vi = m->f[(size_t)t * 3 + (size_t)i];
    for (int a = 0; a < 3; a++)
      p[i][a] = m->v[(size_t)vi * 3 + (size_t)a];
  }
}

double hz_obj_tri_area(const hz_objmesh *m, int32_t t) {
  double p[3][3], e1[3], e2[3], c[3];
  hz_obj_tri(m, t, p);
  for (int a = 0; a < 3; a++) {
    e1[a] = p[1][a] - p[0][a];
    e2[a] = p[2][a] - p[0][a];
  }
  c[0] = e1[1] * e2[2] - e1[2] * e2[1];
  c[1] = e1[2] * e2[0] - e1[0] * e2[2];
  c[2] = e1[0] * e2[1] - e1[1] * e2[0];
  return 0.5 * sqrt(c[0] * c[0] + c[1] * c[1] + c[2] * c[2]);
}

int hz_obj_load(hz_objmesh *m, const char *path, double scale) {
  size_t len = 0;
  char *buf = slurp(path, &len);
  if (buf == NULL) {
    memset(m, 0, sizeof *m);
    return 1;
  }
  int rc = hz_obj_parse(m, buf, scale, path);
  free(buf);
  return rc;
}

int hz_obj_parse(hz_objmesh *m, const char *buf, double scale, const char *objpath) {
  memset(m, 0, sizeof *m);
  if (buf == NULL) return 3;

  /* --- проход 1: счёт. Растущих массивов нет ровно поэтому. --- */
  int64_t nv = 0, nvn = 0, ntri = 0, nquad = 0;
  for (const char *p = buf; *p != '\0'; p = skip_line(p)) {
    const char *q = skip_sp(p);
    if (q[0] == 'v' && is_sp(q[1])) {
      nv++;
    } else if (q[0] == 'v' && q[1] == 'n' && is_sp(q[2])) {
      nvn++;
    } else if (q[0] == 'f' && is_sp(q[1])) {
      int k = count_tokens(q + 1);
      if (k >= 3) {
        ntri += k - 2;
        if (k > 3) nquad++;
      }
    }
  }
  if (nv > 2000000000LL || ntri > 700000000LL) {

    return 3;
  }

  m->v = malloc((size_t)(nv > 0 ? nv : 1) * 3 * sizeof *m->v);
  m->vn = (nvn > 0) ? malloc((size_t)nvn * 3 * sizeof *m->vn) : NULL;
  m->f = malloc((size_t)(ntri > 0 ? ntri : 1) * 3 * sizeof *m->f);
  m->fn = malloc((size_t)(ntri > 0 ? ntri : 1) * 3 * sizeof *m->fn);
  m->fm = malloc((size_t)(ntri > 0 ? ntri : 1) * sizeof *m->fm);
  if (m->v == NULL || m->f == NULL || m->fn == NULL || m->fm == NULL ||
      (nvn > 0 && m->vn == NULL)) {

    hz_obj_free(m);
    return 2;
  }
  /* Материал по умолчанию — индекс 0, поэтому fm валиден всегда (заголовок). */
  if (mtl_add(m, "__default", 9, 0.5) < 0) {

    hz_obj_free(m);
    return 2;
  }

  /* --- проход 2: заполнение --- */
  for (int a = 0; a < 3; a++) {
    m->lo[a] = 1e300;
    m->hi[a] = -1e300;
  }
  int32_t cv = 0, cvn = 0, cur_mtl = 0;
  int64_t ct = 0, ndeg = 0;
  int bad = 0;

  for (const char *p = buf; *p != '\0' && !bad; p = skip_line(p)) {
    const char *q = skip_sp(p);
    if (q[0] == 'v' && is_sp(q[1])) {
      char *e = NULL;
      double x = strtod(q + 1, &e) * scale;
      double y = strtod(e, &e) * scale;
      double z = strtod(e, &e) * scale;
      m->v[(size_t)cv * 3 + 0] = x;
      m->v[(size_t)cv * 3 + 1] = y;
      m->v[(size_t)cv * 3 + 2] = z;
      if (x < m->lo[0]) m->lo[0] = x;
      if (y < m->lo[1]) m->lo[1] = y;
      if (z < m->lo[2]) m->lo[2] = z;
      if (x > m->hi[0]) m->hi[0] = x;
      if (y > m->hi[1]) m->hi[1] = y;
      if (z > m->hi[2]) m->hi[2] = z;
      cv++;
    } else if (q[0] == 'v' && q[1] == 'n' && is_sp(q[2]) && m->vn != NULL) {
      /* `m->vn != NULL` недостижимо ложно — проход 1 считал `vn` тем же
       * предикатом, — но проверка стоит явно: без неё запись в NULL прячется за
       * согласованностью двух проходов, а это не то, на что стоит опираться в
       * разборщике недоверенного входа. */
      char *e = NULL;
      double x = strtod(q + 2, &e);
      double y = strtod(e, &e);
      double z = strtod(e, &e);
      double n = sqrt(x * x + y * y + z * z);
      if (!(n > 0.0)) n = 1.0;
      m->vn[(size_t)cvn * 3 + 0] = x / n;
      m->vn[(size_t)cvn * 3 + 1] = y / n;
      m->vn[(size_t)cvn * 3 + 2] = z / n;
      cvn++;
    } else if (strncmp(q, "usemtl", 6) == 0 && is_sp(q[6])) {
      const char *s = skip_sp(q + 6);
      const char *e = skip_token(s);
      int idx = mtl_find(m, s, (size_t)(e - s));
      cur_mtl = (idx >= 0) ? idx : 0;
    } else if (strncmp(q, "mtllib", 6) == 0 && is_sp(q[6])) {
      const char *s = skip_sp(q + 6);
      const char *e = skip_token(s);
      if (objpath != NULL) mtl_load(m, objpath, s, (size_t)(e - s));
    } else if (q[0] == 'f' && is_sp(q[1])) {
      /* Веер (0, i, i+1): грани входа плоские и выпуклые (заголовок). */
      int32_t first = -1, prev = -1, fnfirst = -1, fnprev = -1;
      const char *s = q + 1;
      for (;;) {
        s = skip_sp(s);
        if (is_eol(*s)) break;
        int32_t iv = 0, ivn = -1;
        int ok = 0;
        s = parse_fvert(s, cv, cvn, &iv, &ivn, &ok);
        if (!ok) {
          bad = 1;
          break;
        }
        s = skip_token(s);
        if (first < 0) {
          first = iv;
          fnfirst = ivn;
        } else if (prev < 0) {
          prev = iv;
          fnprev = ivn;
        } else {
          m->f[(size_t)ct * 3 + 0] = first;
          m->f[(size_t)ct * 3 + 1] = prev;
          m->f[(size_t)ct * 3 + 2] = iv;
          m->fn[(size_t)ct * 3 + 0] = fnfirst;
          m->fn[(size_t)ct * 3 + 1] = fnprev;
          m->fn[(size_t)ct * 3 + 2] = ivn;
          m->fm[ct] = cur_mtl;
          if (hz_obj_tri_area(m, (int32_t)ct) > 0.0) {
            ct++;
          } else {
            ndeg++;
          }
          prev = iv;
          fnprev = ivn;
        }
      }
    }
  }

  if (bad) {
    hz_obj_free(m);
    return 3;
  }
  m->nv = cv;
  m->nvn = cvn;
  m->nt = (int32_t)ct;
  m->ndegen = ndeg;
  m->nquad = nquad;
  if (cv == 0) {
    for (int a = 0; a < 3; a++) {
      m->lo[a] = 0.0;
      m->hi[a] = 0.0;
    }
  }
  return 0;
}

int hz_obj_add_quad(hz_objmesh *m, const char *mtlname, const double c[3], const double n[3],
                    const double eu[3], double hu, double hv) {
  int mi = mtl_add(m, mtlname, strlen(mtlname), 0.05);
  if (mi < 0) return -1;
  double nn = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
  if (!(nn > 0.0) || !(hu > 0.0) || !(hv > 0.0)) return -1;
  double un[3], ev[3], eul = 0.0;
  for (int a = 0; a < 3; a++)
    un[a] = n[a] / nn;
  double dp = eu[0] * un[0] + eu[1] * un[1] + eu[2] * un[2];
  for (int a = 0; a < 3; a++) {
    ev[a] = eu[a] - dp * un[a];
    eul += ev[a] * ev[a];
  }
  eul = sqrt(eul);
  if (!(eul > 0.0)) return -1;
  double e1[3], e2[3];
  for (int a = 0; a < 3; a++)
    e1[a] = ev[a] / eul;
  /* `e2 = n × e1`, чтобы обход шёл против часовой стрелки со стороны нормали. */
  e2[0] = un[1] * e1[2] - un[2] * e1[1];
  e2[1] = un[2] * e1[0] - un[0] * e1[2];
  e2[2] = un[0] * e1[1] - un[1] * e1[0];

  double *nv = realloc(m->v, (size_t)(m->nv + 4) * 3 * sizeof *nv);
  if (nv == NULL) return -1;
  m->v = nv;
  int32_t *nf = realloc(m->f, (size_t)(m->nt + 2) * 3 * sizeof *nf);
  if (nf == NULL) return -1;
  m->f = nf;
  int32_t *nm = realloc(m->fm, (size_t)(m->nt + 2) * sizeof *nm);
  if (nm == NULL) return -1;
  m->fm = nm;
  /* Нормалей вершин у лампы НЕТ СОЗНАТЕЛЬНО: она плоская, и подгонка `n0/nu/nv`
   * даст `nu = nv = 0` точно (§89.3). Если массив нормалей в сцене есть, углам
   * лампы ставится `-1` — «нормали нет». */
  if (m->fn != NULL) {
    int32_t *nfn = realloc(m->fn, (size_t)(m->nt + 2) * 3 * sizeof *nfn);
    if (nfn == NULL) return -1;
    m->fn = nfn;
    for (int i = 0; i < 6; i++)
      m->fn[(size_t)m->nt * 3 + (size_t)i] = -1;
  }
  const double su[4] = {-1.0, 1.0, 1.0, -1.0}, sv[4] = {-1.0, -1.0, 1.0, 1.0};
  int32_t v0 = m->nv;
  for (int k = 0; k < 4; k++)
    for (int a = 0; a < 3; a++)
      m->v[(size_t)(v0 + k) * 3 + (size_t)a] = c[a] + su[k] * hu * e1[a] + sv[k] * hv * e2[a];
  m->nv += 4;
  int32_t t0 = m->nt;
  const int32_t idx[6] = {0, 1, 2, 0, 2, 3};
  for (int i = 0; i < 6; i++)
    m->f[(size_t)t0 * 3 + (size_t)i] = v0 + idx[i];
  m->fm[t0] = mi;
  m->fm[t0 + 1] = mi;
  m->nt += 2;
  /* ГАБАРИТ — ПО ВЕРШИНАМ, а не по `c ± (hu+hv)`: площадка ПЛОСКАЯ, вдоль
   * нормали протяжённости нет вовсе, и раздувание коробки по всем осям подняло
   * потолок сцены с 2.07 до 2.19 м, сорвав замер базы (§91.2). */
  for (int k = 0; k < 4; k++)
    for (int a = 0; a < 3; a++) {
      double p = m->v[(size_t)(v0 + k) * 3 + (size_t)a];
      if (m->lo[a] > p) m->lo[a] = p;
      if (m->hi[a] < p) m->hi[a] = p;
    }
  return mi;
}
