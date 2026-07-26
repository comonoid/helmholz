/* TREFFTZ-DG at the element size the project actually targets, with a material
 * interface. Derivation and numbers: PLAN.md, "ПОСТАНОВКА СМЕНЕНА НА ТРЕФТЦ-DG".
 *
 * WHAT slab2d PROVED AND WHY THIS EXISTS. With a CONTINUOUS (partition-of-unity)
 * carrier basis the volume Galerkin residual is VACUOUS at k W >> 1: the carrier
 * cancels k^2 analytically, so a plane wave — the whole solution — sits in the
 * kernel of the volume operator, while a C^1 basis has continuity built in and so
 * yields no equations either. That leaves NE^2 * ND unknowns against O(NE * ND)
 * boundary conditions. Measured: the exact solution satisfied the assembled
 * system to 3e-13 while LSQR drove the residual to 9e-8 and walked AWAY from it.
 *
 * HERE cells do not overlap and the basis is PURE plane waves, so every basis
 * function solves the Helmholtz equation exactly on a region of ANY shape and
 * there is no volume integral anywhere. Everything lives on the skeleton:
 *     interior face:  [u] = 0   and   [dn u]/(i k) = 0
 *     box wall:       dn u - i k u = h     (impedance; the data come from the
 *                                           exact solution, so termination is
 *                                           NOT what this bench tests)
 * The two interior conditions carry the same units, so there is no penalty
 * parameter anywhere in the formulation.
 *
 * THE MATERIAL INTERFACE IS NOT A SPECIAL CASE. A cell the interface crosses is
 * split into two sub-cells carrying their own wavenumbers, and the interface
 * becomes an ORDINARY FACE between them with the very same two conditions —
 * which are exactly the physical transmission conditions. Nothing is "cut": the
 * apparatus a continuous basis needed (exact phi x phi x plane-wave integration
 * over polygons, a Nitsche term on the cut locus) is not required at all, because
 * there are no volume integrals to restrict. A cut costs only the clipping of
 * SEGMENTS, since every integral is a face integral.
 *
 * TEST FUNCTIONS ARE CONJUGATED, unlike the project's bilinear volume assembly:
 * the tests on a face are traces of the same plane waves that make up the jump,
 * and a bilinear test would give a symmetric rather than Hermitian Gram, so "all
 * rows zero" would stop meaning "the jump is zero". */
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static double now_s(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + 1e-9 * (double)t.tv_nsec;
}

enum { MAXCELL = 40000, MAXSUB = 48000, MAXDIR = 160, MAXDIM = 200000, MAXROW = 2000000 };
static const double LAM = 1.0;

typedef struct {
  double W, k0, theta, n, alpha, dth;
  int nd, ne, oracle, lev, spec, abc;
  double lodk;
  double nhx, nhy, thx, thy;
  /* TWO incident waves with their Fresnel partners: 0,1,2 = incident, reflected,
   * transmitted of the first; 3,4,5 of the second. The BASIS is given only the
   * first triple, so with amp2 != 0 the field contains a component that is
   * genuinely NOT representable — which is how this bench stops measuring an
   * exactly-representable field and starts measuring approximation. */
  double complex amp[256];
  double px[256], py[256];
  double amp2, theta2;
  double ex[4]; /* extra carrier directions, angles in radians, magnitude from the medium */
  int nex;
  int nwave;
  /* MIRROR SCENES. 1 = corridor: mirrors on the y walls, characteristic ends in
   * x, exact field = the guided mode sin(kap(y+L)) exp(i bet x). 2 = dead end:
   * the -x wall is a mirror too, exact field = sin(kap(y+L)) sin(bet(x+L)).
   * The mode index sets how many times the ray bounces along the guide
   * (bounces = kap/bet) WITHOUT changing the number of directions, which is 2
   * and 4. That is the whole claim being measured. */
  int mirror, mode, img, drop;
  double beam, by0; /* beam=0 means the pure mode; beam>0 localises the drive */
  /* A BEAM HAS ANGULAR WIDTH. Two directions can only carry a field of constant
   * amplitude, because the face conditions demand continuity; a localised beam
   * needs the angular spread lam/w around each mode direction, and fan/spread
   * are how much of it the basis is given. */
  int fan, nmode, obj, ss;
  double orad, oxc, oyc;
  double spread;
  double Lbox, kap, bet;
} cfg;

static void fresnel(cfg *c) {
  c->nhx = cos(c->alpha);
  c->nhy = sin(c->alpha);
  c->thx = -sin(c->alpha);
  c->thy = cos(c->alpha);
  double k1 = c->k0, k2 = c->k0 * c->n;
  /* Written UNROLLED with literal indices. The loop form (3*w + j) is correct —
   * nw is 1 or 2 by construction — but gcc's analyzer cannot carry that bound
   * through the index expression and reports a phantom overflow. Six literal
   * assignments cost nothing and leave the gate clean, which is the trade
   * CLAUDE.md prescribes for this class of finding. */
  double ka = k1 * sin(c->theta), kb = k1 * sin(c->theta2);
  double a1 = sqrt(k1 * k1 - ka * ka), a2 = sqrt(k2 * k2 - ka * ka);
  double b1 = sqrt(k1 * k1 - kb * kb), b2 = sqrt(k2 * k2 - kb * kb);
  c->nwave = (fabs(c->amp2) > 0.0) ? 6 : 3;
  c->amp[0] = 1.0;
  c->amp[1] = (a1 - a2) / (a1 + a2);
  c->amp[2] = 2.0 * a1 / (a1 + a2);
  c->px[0] = ka * c->thx - a1 * c->nhx;
  c->py[0] = ka * c->thy - a1 * c->nhy;
  c->px[1] = ka * c->thx + a1 * c->nhx;
  c->py[1] = ka * c->thy + a1 * c->nhy;
  c->px[2] = ka * c->thx - a2 * c->nhx;
  c->py[2] = ka * c->thy - a2 * c->nhy;
  c->amp[3] = c->amp2;
  c->amp[4] = c->amp2 * (b1 - b2) / (b1 + b2);
  c->amp[5] = c->amp2 * 2.0 * b1 / (b1 + b2);
  c->px[3] = kb * c->thx - b1 * c->nhx;
  c->py[3] = kb * c->thy - b1 * c->nhy;
  c->px[4] = kb * c->thx + b1 * c->nhx;
  c->py[4] = kb * c->thy + b1 * c->nhy;
  c->px[5] = kb * c->thx - b2 * c->nhx;
  c->py[5] = kb * c->thy - b2 * c->nhy;
}

static double sdist(const cfg *c, double x, double y) {
  return x * c->nhx + y * c->nhy;
}

/* AN OBJECT, not an infinite plane. obj=1 is a square block (the cross-section of
 * a cube), obj=2 a disc (the cross-section of a sphere). Both sit inside the box
 * and are lit by a plane wave arriving at angle th, so the picture has a stated
 * light direction and casts a shadow.
 * The block is placed on cell boundaries, so NO cell is cut and the material
 * boundary is exactly a set of ordinary interior faces — which in this
 * formulation already carry the transmission conditions. The disc is a staircase
 * at cell resolution and is therefore APPROXIMATE geometry: it is here to be
 * looked at, not to be measured against Mie. */
static int inside_obj(const cfg *c, double x, double y) {
  double dx = x - c->oxc, dy = y - c->oyc;
  if (c->obj == 1) return fabs(dx) < c->orad && fabs(dy) < c->orad;
  if (c->obj == 2) return dx * dx + dy * dy < c->orad * c->orad;
  return 0;
}

/* The two or four plane waves of the mirror scene, in the same amp/px/py slots
 * the rest of the bench already reads: the wall term, the column diagnostic and
 * the oracle basis then need no special case at all. */
static void mirror_setup(cfg *c) {
  double complex i1 = CMPLX(0.0, 1.0);
  double L = c->Lbox;
  c->kap = (double)c->mode * M_PI / (2.0 * L);
  double r2 = c->k0 * c->k0 - c->kap * c->kap;
  if (!(r2 > 0.0)) {
    printf("  ABORT: mode %d is below cutoff in a guide of %.4g lam\n", c->mode, 2.0 * L / LAM);
    exit(1);
  }
  c->bet = sqrt(r2);
  double complex ep = cexp(i1 * c->kap * L), em = cexp(-i1 * c->kap * L);
  if (c->nmode > 1) {
    /* A BEAM AS A SUM OF MODES. Every mode satisfies the mirrors exactly, so the
     * sum does too, and the drive is an exact trace rather than the staircase a
     * per-cell envelope produced (measured: that staircase pinned the residual
     * at W/w, 0.16 for W=2 w=8, and no amount of extra directions moved it).
     * A packet of M modes centred on m0 is localised to about (guide height)/M
     * across, so the direction count of a LOCALISED field is 2 * height/width —
     * it is set by how tight the beam is, not by how far it travels or how many
     * times it bounces. The weight sin(kap (y0+L)) projects a point source at
     * y0; the Gaussian in m keeps the packet smooth. */
    int M = c->nmode;
    c->nwave = 0;
    double sw = (double)M / 3.0;
    for (int j = 0; j < M; j++) {
      int mm = c->mode - M / 2 + j;
      if (mm < 1) continue;
      double kp = (double)mm * M_PI / (2.0 * L);
      double r2m = c->k0 * c->k0 - kp * kp;
      if (!(r2m > 0.0)) continue;
      double bt = sqrt(r2m);
      double dm = ((double)mm - (double)c->mode) / sw;
      double wgt = exp(-dm * dm) * sin(kp * (c->by0 + L));
      if (c->nwave + 2 > 256) break;
      c->amp[c->nwave] = wgt * cexp(i1 * kp * L) / (2.0 * i1);
      c->px[c->nwave] = bt;
      c->py[c->nwave] = kp;
      c->nwave++;
      c->amp[c->nwave] = -wgt * cexp(-i1 * kp * L) / (2.0 * i1);
      c->px[c->nwave] = bt;
      c->py[c->nwave] = -kp;
      c->nwave++;
    }
    printf("  [mirror] BEAM of %d modes around m=%d at y0=%.1f: %d directions, waist ~%.1f lam, "
           "ray angle %.1f deg, bounces %.1f\n",
           M, c->mode, c->by0, c->nwave, 2.0 * L / (double)M, atan2(c->kap, c->bet) * 180.0 / M_PI,
           c->kap / c->bet);
    return;
  }
  if (c->mirror == 1) {
    c->nwave = 2;
    c->amp[0] = ep / (2.0 * i1);
    c->amp[1] = -em / (2.0 * i1);
    c->px[0] = c->bet;
    c->py[0] = c->kap;
    c->px[1] = c->bet;
    c->py[1] = -c->kap;
  } else {
    double complex fp = cexp(i1 * c->bet * L), fm = cexp(-i1 * c->bet * L);
    c->nwave = 4;
    c->amp[0] = ep * fp / (2.0 * i1 * 2.0 * i1);
    c->amp[1] = -em * fp / (2.0 * i1 * 2.0 * i1);
    c->amp[2] = -ep * fm / (2.0 * i1 * 2.0 * i1);
    c->amp[3] = em * fm / (2.0 * i1 * 2.0 * i1);
    c->px[0] = c->bet;
    c->py[0] = c->kap;
    c->px[1] = c->bet;
    c->py[1] = -c->kap;
    c->px[2] = -c->bet;
    c->py[2] = c->kap;
    c->px[3] = -c->bet;
    c->py[3] = -c->kap;
  }
  printf("  [mirror] scene %d  mode m=%d  kap=%.4f bet=%.4f  ray bounces along the guide = %.1f  "
         "directions = %d\n",
         c->mirror, c->mode, c->kap, c->bet, c->kap / c->bet, c->nwave);
}

static double complex uexact(const cfg *c, double x, double y) {
  double complex i1 = CMPLX(0.0, 1.0), u = 0.0;
  if (c->mirror) {
    for (int q = 0; q < c->nwave; q++)
      u += c->amp[q] * cexp(i1 * (c->px[q] * x + c->py[q] * y));
    return u;
  }
  int out = sdist(c, x, y) > 0.0;
  for (int q = 0; q < c->nwave; q++) {
    if (((q % 3) == 2) == (out != 0)) continue; /* transmitted lives inside only */
    u += c->amp[q] * cexp(i1 * (c->px[q] * x + c->py[q] * y));
  }
  return u;
}

/* Integral over the segment P0->P1 of exp(i g.x) ds, closed form. */
static double complex seg_exp(double gx, double gy, double x0, double y0, double x1, double y1) {
  double dx = x1 - x0, dy = y1 - y0;
  double len = sqrt(dx * dx + dy * dy);
  if (!(len > 0.0)) return 0.0;
  double gam = gx * dx + gy * dy;
  double complex e0 = cexp(CMPLX(0.0, 1.0) * (gx * x0 + gy * y0));
  if (fabs(gam) < 1e-8) { /* (e^{ig}-1)/(ig) loses everything at small g */
    double complex s = 1.0, t = 1.0;
    for (int m = 1; m < 12; m++) {
      t *= CMPLX(0.0, 1.0) * gam / (double)(m + 1);
      s += t;
    }
    return len * e0 * s;
  }
  return len * e0 * (cexp(CMPLX(0.0, 1.0) * gam) - 1.0) / (CMPLX(0.0, 1.0) * gam);
}

/* Clip a segment to one side of the interface: side 0 = medium (s<0), 1 = the
 * rest. This is the ONLY geometry a cut costs here, because every integral in the
 * formulation is a face integral. */
static int clip_side(const cfg *c, int side, double *ax, double *ay, double *bx, double *by) {
  double s0 = sdist(c, *ax, *ay), s1 = sdist(c, *bx, *by);
  int in0 = (side == 0) ? (s0 < 0.0) : (s0 > 0.0);
  int in1 = (side == 0) ? (s1 < 0.0) : (s1 > 0.0);
  if (!in0 && !in1) return 0;
  if (in0 && in1) return 1;
  double t = s0 / (s0 - s1);
  double mx = *ax + t * (*bx - *ax), my = *ay + t * (*by - *ay);
  if (in0) {
    *bx = mx;
    *by = my;
  } else {
    *ax = mx;
    *ay = my;
  }
  return 1;
}

typedef struct {
  int cell, side; /* side 0 = medium, 1 = the rest */
  double cx, cy;  /* phase reference: the cell centre */
  int base, nd;   /* first unknown, and how many directions */
  double kx[MAXDIR], ky[MAXDIR];
} sub;

static sub S[MAXSUB];
static int subof[MAXCELL][2];
static double cellx[MAXCELL], celly[MAXCELL], cellh[MAXCELL];
static int celld[MAXCELL];
static double skey[2][MAXCELL];
static int ord[2][MAXCELL];

int main(int argc, char **argv) {
  double t_start = now_s();
  static const char *const KEYS[] = {"W",   "nd",     "ne",   "th",   "it",   "n",      "alpha",
                                     "dth", "oracle", "lev",  "lodk", "amp2", "th2",    "spec",
                                     "ex1", "ex2",    "ex3",  "ex4",  "abc",  "mirror", "mode",
                                     "img", "drop",   "beam", "by0",  "fan",  "spread", "nmode",
                                     "obj", "orad",   "oxc",  "oyc",  "ss",   NULL};
  for (int i = 1; i < argc; i++) {
    const char *eq = strchr(argv[i], '=');
    int ok = 0;
    if (eq)
      for (int k = 0; KEYS[k]; k++)
        if ((size_t)(eq - argv[i]) == strlen(KEYS[k]) &&
            !strncmp(argv[i], KEYS[k], strlen(KEYS[k])))
          ok = 1;
    if (!ok) {
      printf("tdg2d: unknown argument '%s'; keys are", argv[i]);
      for (int k = 0; KEYS[k]; k++)
        printf(" %s", KEYS[k]);
      printf("\n");
      return 1;
    }
    /* THE VALUE MUST PARSE WHOLE. Validating only the key let "ne=8 lev=3"
     * through as ONE argument (zsh does not word-split an unquoted expansion),
     * atof read 8 and silently dropped lev — a sweep that looked like it varied
     * a parameter and did not. That is the project's most frequent artefact
     * signature, arriving this time through the command line. */
    char *end = NULL;
    strtod(eq + 1, &end);
    while (end && (*end == ' ' || *end == '\t'))
      end++;
    if (!end || *end != '\0') {
      printf("tdg2d: argument '%s' has trailing junk after the number\n", argv[i]);
      return 1;
    }
  }
  cfg c;
  memset(&c, 0, sizeof c);
  c.W = 1.0 * LAM;
  c.k0 = 2.0 * M_PI / LAM;
  c.theta = 0.3;
  c.n = 1.0;
  c.alpha = 0.4;
  c.nd = 8;
  c.ne = 8;
  int itmax = 4000;
  for (int i = 1; i < argc; i++) {
    double v = atof(strchr(argv[i], '=') + 1);
    if (!strncmp(argv[i], "W=", 2)) c.W = v * LAM;
    if (!strncmp(argv[i], "nd=", 3)) c.nd = (int)v;
    if (!strncmp(argv[i], "ne=", 3)) c.ne = (int)v;
    if (!strncmp(argv[i], "th=", 3)) c.theta = v;
    if (!strncmp(argv[i], "n=", 2)) c.n = v;
    if (!strncmp(argv[i], "alpha=", 6)) c.alpha = v;
    if (!strncmp(argv[i], "dth=", 4)) c.dth = v;
    if (!strncmp(argv[i], "it=", 3)) itmax = (int)v;
    if (!strncmp(argv[i], "oracle=", 7)) c.oracle = (int)v;
    if (!strncmp(argv[i], "lev=", 4)) c.lev = (int)v;
    if (!strncmp(argv[i], "lodk=", 5)) c.lodk = v;
    if (!strncmp(argv[i], "amp2=", 5)) c.amp2 = v;
    if (!strncmp(argv[i], "th2=", 4)) c.theta2 = v;
    if (!strncmp(argv[i], "spec=", 5)) c.spec = (int)v;
    if (!strncmp(argv[i], "abc=", 4)) c.abc = (int)v;
    if (!strncmp(argv[i], "mirror=", 7)) c.mirror = (int)v;
    if (!strncmp(argv[i], "mode=", 5)) c.mode = (int)v;
    if (!strncmp(argv[i], "img=", 4)) c.img = (int)v;
    if (!strncmp(argv[i], "drop=", 5)) c.drop = (int)v;
    if (!strncmp(argv[i], "beam=", 5)) c.beam = v;
    if (!strncmp(argv[i], "by0=", 4)) c.by0 = v;
    if (!strncmp(argv[i], "fan=", 4)) c.fan = (int)v;
    if (!strncmp(argv[i], "spread=", 7)) c.spread = v;
    if (!strncmp(argv[i], "nmode=", 6)) c.nmode = (int)v;
    if (!strncmp(argv[i], "obj=", 4)) c.obj = (int)v;
    if (!strncmp(argv[i], "orad=", 5)) c.orad = v;
    if (!strncmp(argv[i], "oxc=", 4)) c.oxc = v;
    if (!strncmp(argv[i], "oyc=", 4)) c.oyc = v;
    if (!strncmp(argv[i], "ss=", 3)) c.ss = (int)v;
    for (int e = 0; e < 4; e++) {
      char key[8];
      snprintf(key, sizeof key, "ex%d=", e + 1);
      if (!strncmp(argv[i], key, 4)) {
        c.ex[e] = v;
        if (e + 1 > c.nex) c.nex = e + 1;
      }
    }
  }
  if (c.nd > MAXDIR || c.ne * c.ne > MAXCELL) return 1;
  fresnel(&c);
  int ncell = c.ne * c.ne;
  double L = 0.5 * (double)c.ne * c.W;
  int contrast = fabs(c.n - 1.0) > 1e-12;
  int planecut = contrast && !c.obj; /* the object sits on cell boundaries: nothing to clip */
  if (c.mirror) {
    c.Lbox = L;
    c.abc = 1; /* the open ends are characteristic, and only they */
    mirror_setup(&c);
  }
  if (c.obj) {
    /* ONE incident plane wave — a distant light of stated direction. Everything
     * else in the picture (reflection, refraction, shadow) is solved, not
     * prescribed, so there is no exact reference and the quality measure is the
     * LSQR residual alone. */
    c.abc = 1;
    c.nwave = 1;
    c.amp[0] = 1.0;
    c.px[0] = c.k0 * cos(c.theta);
    c.py[0] = c.k0 * sin(c.theta);
    printf("  [scene] %s of n=%.2f, half-size %.4g lam at (%.4g,%.4g); light from %.1f deg\n",
           c.obj == 1 ? "block" : "disc", c.n, c.orad, c.oxc, c.oyc, c.theta * 180.0 / M_PI);
  }
  printf("tdg2d: W=%.4g lam  kW=%.4g  ND=%d  NE=%d  box=%.4g lam  n=%.3f alpha=%.2f th=%.2f\n",
         c.W / LAM, c.k0 * c.W, c.nd, c.ne, 2.0 * L / LAM, c.n, c.alpha, c.theta);

  /* --- CELLS: a graded ladder, not a grid --------------------------------
   * The architecture sizes elements by LOD, so the mesh has to carry cells of
   * different sizes side by side. In DG that costs nothing structurally: a face
   * between a coarse and a fine cell is just the segment they share, and a
   * coarse cell simply has several faces along the side a fine neighbour
   * touches. "Hanging nodes" are not a concept here.
   * Refinement rule: split while the cell is larger than lodk times its distance
   * to the interface — regulator 2 of the plan (a FLOOR on size driven by
   * geometry), capped by lev. lev = 0 gives back the uniform grid. */
  int ncellv = 0;
  {
    int head = 0;
    for (int a = 0; a < c.ne; a++)
      for (int bb = 0; bb < c.ne; bb++) {
        cellx[ncellv] = -L + ((double)a + 0.5) * c.W;
        celly[ncellv] = -L + ((double)bb + 0.5) * c.W;
        cellh[ncellv] = c.W;
        celld[ncellv] = 0;
        ncellv++;
      }
    while (head < ncellv) {
      int q = head++;
      if (celld[q] >= c.lev) continue;
      double d = fabs(sdist(&c, cellx[q], celly[q]));
      if (!(cellh[q] > c.lodk * d)) continue;
      if (ncellv + 3 >= MAXCELL) break;
      double h2 = 0.5 * cellh[q], x0 = cellx[q], y0 = celly[q];
      int dep = celld[q] + 1;
      cellh[q] = h2;
      celld[q] = dep;
      cellx[q] = x0 - 0.5 * h2;
      celly[q] = y0 - 0.5 * h2;
      const double ox[3] = {0.5, -0.5, 0.5}, oy[3] = {-0.5, 0.5, 0.5};
      for (int t = 0; t < 3; t++) {
        cellx[ncellv] = x0 + ox[t] * h2;
        celly[ncellv] = y0 + oy[t] * h2;
        cellh[ncellv] = h2;
        celld[ncellv] = dep;
        ncellv++;
      }
      head = 0; /* the split cell and its siblings must be reconsidered */
    }
  }
  ncell = ncellv;
  {
    double hmin = 1e300, hmax = 0.0;
    for (int q = 0; q < ncell; q++) {
      if (cellh[q] < hmin) hmin = cellh[q];
      if (cellh[q] > hmax) hmax = cellh[q];
    }
    printf("  cells=%d  sizes %.4g..%.4g lam  (ratio %.0f)\n", ncell, hmin / LAM, hmax / LAM,
           hmax / hmin);
  }

  /* --- sub-cells ---------------------------------------------------------- */
  int nsub = 0, dim = 0;
  {
    for (int q = 0; q < ncell; q++) {
      subof[q][0] = subof[q][1] = -1;
      double smin = 1e300, smax = -1e300;
      for (int i = 0; i < 4; i++) {
        double x = cellx[q] + ((i & 1) ? 0.5 : -0.5) * cellh[q];
        double y = celly[q] + ((i & 2) ? 0.5 : -0.5) * cellh[q];
        double s = sdist(&c, x, y);
        if (s < smin) smin = s;
        if (s > smax) smax = s;
      }
      int cut = contrast && !c.obj && smin < 0.0 && smax > 0.0;
      int own = (contrast && sdist(&c, cellx[q], celly[q]) < 0.0) ? 0 : 1;
      if (c.obj) {
        /* the object is placed on cell boundaries, so a cell is wholly in or
         * wholly out and the material boundary is a set of ordinary faces */
        cut = 0;
        own = inside_obj(&c, cellx[q], celly[q]) ? 0 : 1;
      }
      for (int sd = 0; sd < 2; sd++) {
        if (!cut && sd != own) continue;
        if (nsub >= MAXSUB) return 1;
        sub *s = &S[nsub];
        s->cell = q;
        s->side = sd;
        s->cx = cellx[q];
        s->cy = celly[q];
        double km = (sd == 0 && contrast) ? c.k0 * c.n : c.k0;
        if (c.oracle) {
          /* the field's OWN directions: transmitted inside, incident + reflected
           * outside; with no contrast the two outside waves coincide */
          if (c.mirror) {
            /* two directions for the corridor, four for the dead end — and this
             * number does NOT depend on the bounce count */
            /* drop = the negative control: take away directions the mode needs
             * and it must become unrepresentable. Direction 0 is kept because it
             * is the one the characteristic wall drives through — the wall can
             * only prescribe an incoming amplitude for a direction that IS in
             * the basis, so dropping it would silence the source instead of
             * breaking the approximation, and the control would be vacuous. */
            s->nd = c.nwave - c.drop;
            if (s->nd < 1) s->nd = 1;
            for (int d = 0; d < s->nd; d++) {
              s->kx[d] = c.px[d];
              s->ky[d] = c.py[d];
            }
            /* the fan, appended AFTER the exact directions so the wall drive
             * still finds direction 0 by exact match */
            int nb = s->nd;
            for (int b = 0; b < nb; b++) {
              double a0 = atan2(c.py[b], c.px[b]);
              for (int j = 1; j <= c.fan && s->nd < MAXDIR - 1; j++)
                for (int sg = -1; sg <= 1; sg += 2) {
                  double aa = a0 + (double)sg * c.spread * (double)j / (double)c.fan;
                  s->kx[s->nd] = c.k0 * cos(aa);
                  s->ky[s->nd] = c.k0 * sin(aa);
                  s->nd++;
                }
            }
          } else if (!contrast) {
            /* the incident wave's OWN direction, which fresnel() defines against
             * the interface normal — not the raw angle theta. dth rotates it
             * deliberately, which is how the direction tolerance is measured. */
            double a0 = atan2(c.py[0], c.px[0]) + c.dth;
            s->nd = 1;
            s->kx[0] = c.k0 * cos(a0);
            s->ky[0] = c.k0 * sin(a0);
          } else if (sd == 0) {
            s->nd = 1;
            s->kx[0] = c.px[2];
            s->ky[0] = c.py[2];
          } else {
            s->nd = 2;
            for (int d = 0; d < 2; d++) {
              s->kx[d] = c.px[d];
              s->ky[d] = c.py[d];
            }
          }
        } else {
          s->nd = c.nd;
          for (int d = 0; d < c.nd; d++) {
            /* With an object the fan is ANCHORED to the illumination direction,
             * so d=0 is exactly the incident wave. The characteristic wall can
             * only prescribe an amplitude for a direction that is in the basis,
             * and an unanchored fan therefore leaves the scene unlit. */
            double th = c.obj ? c.theta + 2.0 * M_PI * (double)d / (double)c.nd
                              : 2.0 * M_PI * ((double)d + 0.5) / (double)c.nd;
            s->kx[d] = km * cos(th);
            s->ky[d] = km * sin(th);
          }
        }
        /* EXTRA DIRECTIONS, appended to every sub-cell. This is the mechanism the
         * residual-driven search needs: propose an angle, add it everywhere,
         * see what the residual does. The magnitude follows the sub-cell's own
         * medium, so one angle serves both sides of an interface. */
        for (int e = 0; e < c.nex && s->nd < MAXDIR; e++) {
          s->kx[s->nd] = km * cos(c.ex[e]);
          s->ky[s->nd] = km * sin(c.ex[e]);
          s->nd++;
        }
        s->base = dim;
        dim += s->nd;
        subof[q][sd] = nsub++;
      }
    }
  }
  if (dim > MAXDIM) {
    printf("tdg2d: dim %d too large\n", dim);
    return 1;
  }
  printf("  cells=%d  sub-cells=%d  unknowns=%d\n", ncell, nsub, dim);

  size_t cap = 48u << 20;
  int *ja = malloc(cap * sizeof(int));
  double complex *va = malloc(cap * sizeof(double complex));
  size_t *rp = malloc((MAXROW + 1) * sizeof(size_t));
  double complex *rhs = calloc(MAXROW, sizeof(double complex));
  if (!ja || !va || !rp || !rhs) return 1;
  size_t nnz = 0;
  int nrow = 0;
  double complex i1 = CMPLX(0.0, 1.0);

  /* One face. sb < 0 marks a box wall, where the impedance condition and its
   * analytic drive replace the two jump conditions; wmask selects which exact
   * waves live on that side. */
#define ADDFACE(SA, SB, X0, Y0, X1, Y1, NX, NY, WMASK)                                             \
  do {                                                                                             \
    int scs[2] = {(SA), (SB)};                                                                     \
    int isint = ((SB) >= 0);                                                                       \
    for (int ti = 0; ti < 2; ti++) {                                                               \
      if (scs[ti] < 0) continue;                                                                   \
      const sub *ST = &S[scs[ti]];                                                                 \
      for (int tm = 0; tm < ST->nd; tm++) {                                                        \
        double tgx = -ST->kx[tm], tgy = -ST->ky[tm];                                               \
        double complex tph = cexp(-i1 * (tgx * ST->cx + tgy * ST->cy));                            \
        for (int cond = 0; cond < (isint ? 2 : 1); cond++) {                                       \
          if (nrow >= MAXROW) return 1;                                                            \
          rp[nrow] = nnz;                                                                          \
          for (int t2 = 0; t2 < 2; t2++) {                                                         \
            if (scs[t2] < 0) continue;                                                             \
            const sub *SU = &S[scs[t2]];                                                           \
            double sgn = isint ? (t2 == 0 ? -1.0 : 1.0) : 1.0;                                     \
            for (int d = 0; d < SU->nd; d++) {                                                     \
              double gx = SU->kx[d] + tgx, gy = SU->ky[d] + tgy;                                   \
              double complex ph = cexp(-i1 * (SU->kx[d] * SU->cx + SU->ky[d] * SU->cy));           \
              double complex Iv = seg_exp(gx, gy, (X0), (Y0), (X1), (Y1)) * ph * tph;              \
              double dn = SU->kx[d] * (NX) + SU->ky[d] * (NY);                                     \
              double complex w;                                                                    \
              if (!isint)                                                                          \
                w = (i1 * dn - i1 * c.k0) * Iv;                                                    \
              else if (cond == 0)                                                                  \
                w = sgn * Iv;                                                                      \
              else                                                                                 \
                w = sgn * i1 * dn * Iv / (i1 * c.k0);                                              \
              if (!(cabs(w) > 0.0) || nnz >= cap) continue;                                        \
              ja[nnz] = SU->base + d;                                                              \
              va[nnz++] = w;                                                                       \
            }                                                                                      \
          }                                                                                        \
          double complex dr = 0.0;                                                                 \
          if (!isint)                                                                              \
            for (int q2 = 0; q2 < c.nwave; q2++) {                                                 \
              if (!((WMASK) & (1 << q2))) continue;                                                \
              double gx = c.px[q2] + tgx, gy = c.py[q2] + tgy;                                     \
              double complex Iv = seg_exp(gx, gy, (X0), (Y0), (X1), (Y1)) * tph;                   \
              dr += c.amp[q2] * (i1 * (c.px[q2] * (NX) + c.py[q2] * (NY)) - i1 * c.k0) * Iv;       \
            }                                                                                      \
          rhs[nrow] = dr;                                                                          \
          nrow++;                                                                                  \
        }                                                                                          \
      }                                                                                            \
    }                                                                                              \
  } while (0)

  /* THE TERMINATION, AND WHY KNOWING THE DIRECTIONS MAKES IT EXACT.
   * The impedance condition dn u - i k u = h is only first order: it absorbs a
   * normally incident wave and reflects an oblique one, and its data h had to
   * carry the OUTGOING field as well — which a real scene does not know.
   * With a directional basis the split is exact instead of approximate: on a
   * wall with outward normal n every basis direction is unambiguously OUTGOING
   * (d.n > 0) or INCOMING (d.n < 0). Outgoing waves must leave untouched, so
   * they get no equation at all; incoming ones are exactly the illumination
   * entering the domain, so their amplitudes are prescribed. Nothing reflects,
   * and the only datum needed is what a scene actually knows: what comes IN. */
/* A PERFECT MIRROR is one more kind of wall and costs one equation per test
 * function: Int_F u conj(T) ds = 0, i.e. Dirichlet. Nothing else changes — and
 * that is the point of the scene family below: a wave bouncing between mirrors
 * an unbounded number of times is still only a FEW directions, so the cost does
 * not grow with the number of bounces at all. A path tracer pays depth and
 * variance per bounce; here the bounces are not traversed, they are solved. */
#define ADDMIRROR(SA, X0, Y0, X1, Y1)                                                              \
  do {                                                                                             \
    const sub *SM = &S[(SA)];                                                                      \
    for (int tm = 0; tm < SM->nd; tm++) {                                                          \
      double tgx = -SM->kx[tm], tgy = -SM->ky[tm];                                                 \
      double complex tph = cexp(-i1 * (tgx * SM->cx + tgy * SM->cy));                              \
      if (nrow >= MAXROW) return 1;                                                                \
      rp[nrow] = nnz;                                                                              \
      for (int d = 0; d < SM->nd; d++) {                                                           \
        double gx = SM->kx[d] + tgx, gy = SM->ky[d] + tgy;                                         \
        double complex ph = cexp(-i1 * (SM->kx[d] * SM->cx + SM->ky[d] * SM->cy));                 \
        double complex w = seg_exp(gx, gy, (X0), (Y0), (X1), (Y1)) * ph * tph;                     \
        if (!(cabs(w) > 0.0) || nnz >= cap) continue;                                              \
        ja[nnz] = SM->base + d;                                                                    \
        va[nnz++] = w;                                                                             \
      }                                                                                            \
      rhs[nrow] = 0.0;                                                                             \
      nrow++;                                                                                      \
    }                                                                                              \
  } while (0)

#define ADDWALL(SA, X0, Y0, X1, Y1, NX, NY, WMASK)                                                 \
  do {                                                                                             \
    if (c.mirror && (fabs(NY) > 0.0 || (c.mirror == 2 && (NX) < 0.0))) {                           \
      ADDMIRROR((SA), (X0), (Y0), (X1), (Y1)); /* the guide walls, and the dead end */             \
      break;                                                                                       \
    }                                                                                              \
    if (!c.abc) {                                                                                  \
      ADDFACE((SA), -1, (X0), (Y0), (X1), (Y1), (NX), (NY), (WMASK));                              \
      break;                                                                                       \
    }                                                                                              \
    const sub *SW = &S[(SA)];                                                                      \
    for (int d = 0; d < SW->nd; d++) {                                                             \
      if (SW->kx[d] * (NX) + SW->ky[d] * (NY) >= 0.0) continue; /* outgoing: free */               \
      if (nrow >= MAXROW || nnz >= cap) return 1;                                                  \
      rp[nrow] = nnz;                                                                              \
      ja[nnz] = SW->base + d;                                                                      \
      va[nnz++] = 1.0;                                                                             \
      double complex vin = 0.0;                                                                    \
      for (int q2 = 0; q2 < c.nwave; q2++) {                                                       \
        /* the mirror scene has no inside/outside split, so no mask applies */                     \
        if (!c.mirror && !c.obj && !((WMASK) & (1 << q2))) continue;                               \
        if (fabs(SW->kx[d] - c.px[q2]) + fabs(SW->ky[d] - c.py[q2]) > 1e-9 * c.k0) continue;       \
        vin = c.amp[q2] * cexp(i1 * (c.px[q2] * SW->cx + c.py[q2] * SW->cy));                      \
        break;                                                                                     \
      }                                                                                            \
      rhs[nrow] = vin;                                                                             \
      nrow++;                                                                                      \
    }                                                                                              \
  } while (0)

  /* INTER-CELL FACES BY GEOMETRIC ADJACENCY, not by grid indices. A face is the
   * OVERLAP of two cell sides that lie on the same line, so a coarse cell facing
   * several fine ones simply yields several faces along that side, each of the
   * fine cell's length. Emitted once per pair: only the +x and +y sides of the
   * cell behind the normal are scanned. */
  int nfint = 0, nfwall = 0;
  for (int ax2 = 0; ax2 < 2; ax2++) {
    for (int q = 0; q < ncell; q++) {
      skey[ax2][q] = (ax2 == 0 ? cellx[q] : celly[q]) - 0.5 * cellh[q];
      ord[ax2][q] = q;
    }
    for (int i = 1; i < ncell; i++) { /* insertion sort on a nearly sorted list */
      int v = ord[ax2][i];
      double kv = skey[ax2][v];
      int j = i - 1;
      while (j >= 0 && skey[ax2][ord[ax2][j]] > kv) {
        ord[ax2][j + 1] = ord[ax2][j];
        j--;
      }
      ord[ax2][j + 1] = v;
    }
  }
  for (int qa = 0; qa < ncell; qa++)
    for (int ax2 = 0; ax2 < 2; ax2++) {
      double ha = 0.5 * cellh[qa];
      double coord = (ax2 == 0) ? cellx[qa] + ha : celly[qa] + ha; /* the +x / +y side */
      double lo = (ax2 == 0) ? celly[qa] - ha : cellx[qa] - ha, hi = lo + cellh[qa];
      double nx = (ax2 == 0) ? 1.0 : 0.0, ny = (ax2 == 0) ? 0.0 : 1.0;
      /* the box wall: no neighbour there by construction */
      if (fabs(coord - L) < 1e-9 * c.W) {
        for (int sd = 0; sd < 2; sd++) {
          if (c.obj ? sd != 0 : !contrast && sd == 0) continue; /* obj: one sub per cell */
          double px0 = (ax2 == 0) ? coord : lo, py0 = (ax2 == 0) ? lo : coord;
          double px1 = (ax2 == 0) ? coord : hi, py1 = (ax2 == 0) ? hi : coord;
          if (planecut && !clip_side(&c, sd, &px0, &py0, &px1, &py1)) continue;
          int sa = c.obj ? (subof[qa][0] >= 0 ? subof[qa][0] : subof[qa][1]) : subof[qa][sd];
          if (sa < 0) continue;
          int mask = (sd == 0) ? 0x24 : 0x1b; /* inside: q=2,5   outside: q=0,1,3,4 */
          nfwall++;
          ADDWALL(sa, px0, py0, px1, py1, nx, ny, mask);
        }
        continue;
      }
      /* Neighbours by binary search on the sorted "-x / -y side" coordinate.
       * A quadratic scan caps the mesh at a couple of thousand cells, which is
       * two orders below what the scaling question needs. */
      int lo_i = 0, hi_i = ncell;
      while (lo_i < hi_i) {
        int mid = (lo_i + hi_i) / 2;
        if (skey[ax2][ord[ax2][mid]] < coord - 1e-9 * c.W)
          lo_i = mid + 1;
        else
          hi_i = mid;
      }
      for (int ii = lo_i; ii < ncell; ii++) {
        int qb = ord[ax2][ii];
        if (skey[ax2][qb] > coord + 1e-9 * c.W) break;
        if (qb == qa) continue;
        double hb = 0.5 * cellh[qb];
        double blo = (ax2 == 0) ? celly[qb] - hb : cellx[qb] - hb, bhi = blo + cellh[qb];
        double olo = lo > blo ? lo : blo, ohi = hi < bhi ? hi : bhi;
        if (!(ohi - olo > 1e-9 * c.W)) continue;
        for (int sd = 0; sd < 2; sd++) {
          if (c.obj ? sd != 0 : !contrast && sd == 0) continue;
          double px0 = (ax2 == 0) ? coord : olo, py0 = (ax2 == 0) ? olo : coord;
          double px1 = (ax2 == 0) ? coord : ohi, py1 = (ax2 == 0) ? ohi : coord;
          if (planecut && !clip_side(&c, sd, &px0, &py0, &px1, &py1)) continue;
          int sa = c.obj ? (subof[qa][0] >= 0 ? subof[qa][0] : subof[qa][1]) : subof[qa][sd];
          int sb = c.obj ? (subof[qb][0] >= 0 ? subof[qb][0] : subof[qb][1]) : subof[qb][sd];
          int mask = (sd == 0) ? 0x24 : 0x1b; /* inside: q=2,5   outside: q=0,1,3,4 */
          if (sa >= 0 && sb >= 0) {
            nfint++;
            ADDFACE(sa, sb, px0, py0, px1, py1, nx, ny, mask);
          }
        }
      }
    }
  /* the -x and -y sides that lie on the box wall */
  for (int qa = 0; qa < ncell; qa++)
    for (int ax2 = 0; ax2 < 2; ax2++) {
      double ha = 0.5 * cellh[qa];
      double coord = (ax2 == 0) ? cellx[qa] - ha : celly[qa] - ha;
      if (fabs(coord + L) > 1e-9 * c.W) continue;
      double lo = (ax2 == 0) ? celly[qa] - ha : cellx[qa] - ha, hi = lo + cellh[qa];
      double nx = (ax2 == 0) ? -1.0 : 0.0, ny = (ax2 == 0) ? 0.0 : -1.0;
      for (int sd = 0; sd < 2; sd++) {
        if (c.obj ? sd != 0 : !contrast && sd == 0) continue;
        double px0 = (ax2 == 0) ? coord : lo, py0 = (ax2 == 0) ? lo : coord;
        double px1 = (ax2 == 0) ? coord : hi, py1 = (ax2 == 0) ? hi : coord;
        if (planecut && !clip_side(&c, sd, &px0, &py0, &px1, &py1)) continue;
        int sa = c.obj ? (subof[qa][0] >= 0 ? subof[qa][0] : subof[qa][1]) : subof[qa][sd];
        if (sa < 0) continue;
        int mask = (sd == 0) ? 0x24 : 0x1b;
        nfwall++;
        ADDWALL(sa, px0, py0, px1, py1, nx, ny, mask);
      }
    }

  /* THE MATERIAL INTERFACE, as an ordinary face between the two sub-cells of a
   * cut cell. Same two conditions as every interior face — they ARE the physical
   * transmission conditions. */
  int nif = 0;
  if (planecut)
    for (int q = 0; q < ncell; q++) {
      if (subof[q][0] < 0 || subof[q][1] < 0) continue;
      double sc = sdist(&c, cellx[q], celly[q]);
      double ox = cellx[q] - sc * c.nhx, oy = celly[q] - sc * c.nhy;
      double t0 = -1e300, t1 = 1e300;
      const double lo[2] = {cellx[q] - 0.5 * c.W, celly[q] - 0.5 * c.W};
      const double hi[2] = {cellx[q] + 0.5 * c.W, celly[q] + 0.5 * c.W};
      const double dd[2] = {c.thx, c.thy}, oo[2] = {ox, oy};
      int empty = 0;
      for (int ax2 = 0; ax2 < 2; ax2++) {
        if (fabs(dd[ax2]) < 1e-300) {
          if (oo[ax2] < lo[ax2] || oo[ax2] > hi[ax2]) empty = 1;
          continue;
        }
        double ta = (lo[ax2] - oo[ax2]) / dd[ax2], tb = (hi[ax2] - oo[ax2]) / dd[ax2];
        if (ta > tb) {
          double sw = ta;
          ta = tb;
          tb = sw;
        }
        if (ta > t0) t0 = ta;
        if (tb < t1) t1 = tb;
      }
      if (empty || !(t0 < t1)) continue;
      double ax = ox + c.thx * t0, ay = oy + c.thy * t0;
      double bx = ox + c.thx * t1, by = oy + c.thy * t1;
      nif++;
      ADDFACE(subof[q][0], subof[q][1], ax, ay, bx, by, c.nhx, c.nhy, 0);
    }
  rp[nrow] = nnz;
  printf("  rows=%d  (%.1fx overdetermined)  nnz=%zu  interface faces=%d\n", nrow,
         (double)nrow / (double)dim, nnz, nif);

  for (int r = 0; r < nrow; r++) {
    double s = 0.0;
    for (size_t p = rp[r]; p < rp[r + 1]; p++)
      s += cabs(va[p]) * cabs(va[p]);
    s = sqrt(s);
    if (!(s > 0.0)) continue;
    for (size_t p = rp[r]; p < rp[r + 1]; p++)
      va[p] /= s;
    rhs[r] /= s;
  }

  /* With the field's own directions and no deliberate angular error the exact
   * solution is IN THE SPAN exactly, so the assembled rows must reproduce it. */
  if (c.oracle && !(fabs(c.dth) > 0.0)) {
    double complex *ce = calloc((size_t)dim, sizeof(double complex));
    if (ce) {
      /* Match each basis direction to a wave of the exact field BY ITS VECTOR,
       * not by its index. With extra directions the index correspondence is
       * gone, and this is the only way the diagnostic survives enrichment: a
       * direction that matches a wave gets that wave's amplitude, one that
       * matches nothing gets zero. If the exact field is in the span, this
       * vector satisfies the system. */
      int matched = 0;
      for (int s = 0; s < nsub; s++) {
        const sub *SU = &S[s];
        for (int d = 0; d < SU->nd; d++) {
          ce[SU->base + d] = 0.0;
          for (int q2 = 0; q2 < c.nwave; q2++) {
            int wave_inside = ((q2 % 3) == 2);
            if (contrast && wave_inside != (SU->side == 0)) continue;
            if (fabs(SU->kx[d] - c.px[q2]) + fabs(SU->ky[d] - c.py[q2]) > 1e-9 * c.k0) continue;
            ce[SU->base + d] = c.amp[q2] * cexp(i1 * (c.px[q2] * SU->cx + c.py[q2] * SU->cy));
            matched++;
            break;
          }
        }
      }
      printf("  [diag] %d of %d columns matched a wave of the exact field\n", matched, dim);
      double e2 = 0.0, nb = 0.0;
      for (int r = 0; r < nrow; r++) {
        double complex s = 0.0;
        for (size_t p = rp[r]; p < rp[r + 1]; p++)
          s += va[p] * ce[ja[p]];
        e2 += cabs(s - rhs[r]) * cabs(s - rhs[r]);
        nb += cabs(rhs[r]) * cabs(rhs[r]);
      }
      printf("  [diag] exact solution in assembled rows: |Ac-b| = %.3e of |b| = %.3e\n", sqrt(e2),
             sqrt(nb));
      free(ce);
    }
  }

  /* --- LSQR --------------------------------------------------------------- */
  double t_asm = now_s() - t_start, t_solve0 = now_s();
  double complex *xs = calloc((size_t)dim, sizeof(double complex));
  double complex *uu = calloc((size_t)nrow, sizeof(double complex));
  double complex *vv = calloc((size_t)dim, sizeof(double complex));
  double complex *ww = calloc((size_t)dim, sizeof(double complex));
  double complex *tv = calloc((size_t)dim, sizeof(double complex));
  if (!xs || !uu || !vv || !ww || !tv) return 1;
  for (int r = 0; r < nrow; r++)
    uu[r] = rhs[r];
  double beta = 0.0;
  for (int r = 0; r < nrow; r++)
    beta += cabs(uu[r]) * cabs(uu[r]);
  beta = sqrt(beta);
  if (!(beta > 0.0)) {
    printf("  ABORT: zero drive\n");
    return 1;
  }
  for (int r = 0; r < nrow; r++)
    uu[r] /= beta;
  for (int j = 0; j < dim; j++)
    vv[j] = 0.0;
  for (int r = 0; r < nrow; r++)
    for (size_t p = rp[r]; p < rp[r + 1]; p++)
      vv[ja[p]] += conj(va[p]) * uu[r];
  double alpha = 0.0;
  for (int j = 0; j < dim; j++)
    alpha += cabs(vv[j]) * cabs(vv[j]);
  alpha = sqrt(alpha);
  for (int j = 0; j < dim; j++) {
    vv[j] /= alpha;
    ww[j] = vv[j];
  }
  double phibar = beta, rhobar = alpha, b0 = beta;
  int n2 = 0, n3 = 0, n4 = 0, n6 = 0, n8 = 0, brk = 0, itdone = 0;
  for (int it = 0; it < itmax; it++) {
    itdone = it + 1;
    for (int r = 0; r < nrow; r++) {
      double complex s = 0.0;
      for (size_t p = rp[r]; p < rp[r + 1]; p++)
        s += va[p] * vv[ja[p]];
      uu[r] = s - alpha * uu[r];
    }
    beta = 0.0;
    for (int r = 0; r < nrow; r++)
      beta += cabs(uu[r]) * cabs(uu[r]);
    beta = sqrt(beta);
    if (!(beta > 0.0)) {
      brk = 1;
      break;
    }
    for (int r = 0; r < nrow; r++)
      uu[r] /= beta;
    for (int j = 0; j < dim; j++)
      tv[j] = 0.0;
    for (int r = 0; r < nrow; r++)
      for (size_t p = rp[r]; p < rp[r + 1]; p++)
        tv[ja[p]] += conj(va[p]) * uu[r];
    alpha = 0.0;
    for (int j = 0; j < dim; j++) {
      vv[j] = tv[j] - beta * vv[j];
      alpha += cabs(vv[j]) * cabs(vv[j]);
    }
    alpha = sqrt(alpha);
    if (!(alpha > 0.0)) {
      brk = 2;
      break;
    }
    for (int j = 0; j < dim; j++)
      vv[j] /= alpha;
    double rho = sqrt(rhobar * rhobar + beta * beta);
    double cs = rhobar / rho, sn = beta / rho;
    double th2 = sn * alpha;
    rhobar = -cs * alpha;
    double phi = cs * phibar;
    phibar = sn * phibar;
    for (int j = 0; j < dim; j++) {
      xs[j] += (phi / rho) * ww[j];
      ww[j] = vv[j] - (th2 / rho) * ww[j];
    }
    double rr = phibar / b0;
    if (!n2 && rr <= 1e-2) n2 = it + 1;
    if (!n3 && rr <= 1e-3) n3 = it + 1;
    if (!n4 && rr <= 1e-4) n4 = it + 1;
    if (!n6 && rr <= 1e-6) n6 = it + 1;
    if (!n8 && rr <= 1e-8) n8 = it + 1;
  }
  printf("  LSQR |r|/|b| = %.3e   iters to 1e-2/1e-3/1e-4/1e-6/1e-8: %d/%d/%d/%d/%d\n", phibar / b0,
         n2, n3, n4, n6, n8);
  double t_solve = now_s() - t_solve0;
  printf("  LSQR stopped after %d iters: %s\n", itdone,
         brk == 1   ? "BETA breakdown (Krylov space exhausted)"
         : brk == 2 ? "ALPHA breakdown"
                    : "iteration cap");

  /* --- WHICH DIRECTION IS MISSING? A MATCHED FILTER ON THE RESIDUAL --------
   * A direction absent from the basis cannot be reproduced, so it survives in
   * the residual of the boundary condition; along a wall its phase runs as
   * k (d.t) s, and referenced to absolute position it is exp(i k d.x). So
   * projecting the sampled residual onto exp(i k d.x) over all walls and
   * scanning d over the circle is a matched filter whose peak NAMES the missing
   * direction.
   * If this works, directions need not be supplied by geometry at all: they can
   * be DISCOVERED from a solve and added where the residual asks for them. That
   * is the whole point of doing it this way rather than building an a priori
   * estimator — an a priori rule cannot know about edges, caustics or multiple
   * scattering, and the residual does. */
  if (c.spec > 0) {
    /* TWO wavenumber circles, not one. A wave living inside the denser medium
     * has |k| = n k0 and is INVISIBLE to a filter that scans only |k| = k0 —
     * measured: with only the transmitted wave missing, the peak landed nowhere
     * near it. */
    enum { NPHI = 1440, NSMP = 64, NKC = 2 };
    static double complex acc[NKC][NPHI];
    const double kc[NKC] = {c.k0, c.k0 * c.n};
    for (int m = 0; m < NKC; m++)
      for (int p = 0; p < NPHI; p++)
        acc[m][p] = 0.0;
    double norm = 0.0;
    for (int q = 0; q < ncell; q++)
      for (int ax2 = 0; ax2 < 2; ax2++)
        for (int sgn = -1; sgn <= 1; sgn += 2) {
          double h = 0.5 * cellh[q];
          double coord = (ax2 == 0) ? cellx[q] + sgn * h : celly[q] + sgn * h;
          if (fabs(fabs(coord) - L) > 1e-9 * c.W) continue;
          double nx = (ax2 == 0) ? (double)sgn : 0.0, ny = (ax2 == 0) ? 0.0 : (double)sgn;
          double lo = (ax2 == 0) ? celly[q] - h : cellx[q] - h;
          for (int sd = 0; sd < 2; sd++) {
            int sa = subof[q][sd];
            if (sa < 0) continue;
            const sub *SU = &S[sa];
            for (int t = 0; t < NSMP; t++) {
              double u = lo + cellh[q] * ((double)t + 0.5) / (double)NSMP;
              double x = (ax2 == 0) ? coord : u, y = (ax2 == 0) ? u : coord;
              if (planecut && ((SU->side == 0) != (sdist(&c, x, y) < 0.0))) continue;
              double complex uu2 = 0.0, du = 0.0;
              for (int d = 0; d < SU->nd; d++) {
                double complex e = cexp(i1 * (SU->kx[d] * (x - SU->cx) + SU->ky[d] * (y - SU->cy)));
                uu2 += xs[SU->base + d] * e;
                du += xs[SU->base + d] * i1 * (SU->kx[d] * nx + SU->ky[d] * ny) * e;
              }
              double complex hh = 0.0;
              for (int qq = 0; qq < c.nwave; qq++) {
                if (((qq % 3) == 2) == (sdist(&c, x, y) > 0.0)) continue;
                hh += c.amp[qq] * (i1 * (c.px[qq] * nx + c.py[qq] * ny) - i1 * c.k0) *
                      cexp(i1 * (c.px[qq] * x + c.py[qq] * y));
              }
              double complex r = (du - i1 * c.k0 * uu2) - hh;
              double w = cellh[q] / (double)NSMP;
              norm += cabs(r) * cabs(r) * w;
              for (int m = 0; m < NKC; m++)
                for (int p = 0; p < NPHI; p++) {
                  double ph = 2.0 * M_PI * (double)p / (double)NPHI;
                  double gx = kc[m] * cos(ph), gy = kc[m] * sin(ph);
                  acc[m][p] += w * r * cexp(-i1 * (gx * x + gy * y));
                }
            }
          }
        }
    int best = 0, bestm = 0;
    double bv = 0.0, tot = 0.0;
    for (int m = 0; m < NKC; m++)
      for (int p = 0; p < NPHI; p++) {
        double v = cabs(acc[m][p]);
        tot += v * v;
        if (v > bv) {
          bv = v;
          best = p;
          bestm = m;
        }
      }
    double bphi = 2.0 * M_PI * (double)best / (double)NPHI;
    printf("  [spectrum] wall residual %.3e; peak at %.5f rad on |k|=%s, sharpness %.3f\n",
           sqrt(norm), bphi, bestm ? "n*k0" : "k0", bv / sqrt(tot / (double)(NKC * NPHI)));
    /* what the peak SHOULD be if it is finding the wave the basis lacks */
    for (int q = 3; q < c.nwave; q++)
      printf("             missing wave %d true direction %.15f rad\n", q, atan2(c.py[q], c.px[q]));
  }

  /* --- error against Fresnel, sampled inside each sub-cell's own region ---- */
  double num = 0.0, den = 0.0;
  int nsamp = 0;
  for (int s = 0; s < nsub; s++) {
    const sub *SU = &S[s];
    for (int p = 0; p < 25; p++) {
      double ox = ((double)(p % 5) - 2.0) * 0.2 * c.W, oy = ((double)(p / 5) - 2.0) * 0.2 * c.W;
      double x = SU->cx + ox, y = SU->cy + oy;
      if (planecut && ((SU->side == 0) != (sdist(&c, x, y) < 0.0))) continue;
      double complex got = 0.0;
      for (int d = 0; d < SU->nd; d++)
        got += xs[SU->base + d] * cexp(i1 * (SU->kx[d] * (x - SU->cx) + SU->ky[d] * (y - SU->cy)));
      double complex ex = uexact(&c, x, y);
      num += cabs(got - ex) * cabs(got - ex);
      den += cabs(ex) * cabs(ex);
      nsamp++;
    }
  }
  printf("  FIELD ERROR = %.4e   (%d samples)\n", sqrt(num / den), nsamp);

  /* --- RENDER --------------------------------------------------------------
   * The picture is not a second solve. The coefficients already ARE the field
   * everywhere, in closed form, so a pixel costs one evaluation of its own
   * cell's ND plane waves — nothing else, at any resolution. That is the
   * architecture's "the field does not depend on the camera" made literal, and
   * it is why solve time and render time are reported SEPARATELY here: quoting
   * one as the other is exactly the kind of cheating PLAN.md forbids.
   * The loop is over CELLS, not pixels: a pixel lookup would need a tree walk,
   * while scattering each cell into the pixels it covers is O(pixels) flat. */
  double t_img = 0.0;
  if (c.img > 0) {
    int N = c.img;
    double t_img0 = now_s();
    float *fre = malloc((size_t)N * (size_t)N * sizeof(float));
    float *fab = malloc((size_t)N * (size_t)N * sizeof(float));
    if (!fre || !fab) {
      free(fre);
      free(fab);
      printf("  ABORT: image buffers\n");
      exit(1);
    }
    for (size_t p = 0; p < (size_t)N * (size_t)N; p++) {
      fre[p] = 0.0f;
      fab[p] = 0.0f;
    }
    double px = 2.0 * L / (double)N;
    for (int s = 0; s < nsub; s++) {
      const sub *SU = &S[s];
      int q = SU->cell;
      double h = 0.5 * cellh[q];
      int i0 = (int)floor((cellx[q] - h + L) / px), i1x = (int)ceil((cellx[q] + h + L) / px);
      int j0 = (int)floor((celly[q] - h + L) / px), j1 = (int)ceil((celly[q] + h + L) / px);
      if (i0 < 0) i0 = 0;
      if (j0 < 0) j0 = 0;
      if (i1x > N) i1x = N;
      if (j1 > N) j1 = N;
      /* A PIXEL IS AN AREA, NOT A POINT. A sensor integrates |u|^2 over its own
       * cell, and a coherent field oscillates at lam/2 — sampling it at one
       * point per pixel produces the moire that made the first pictures
       * unreadable. ss x ss subsamples per pixel is that integration; the phase
       * picture keeps the centre sample, since averaging a signed oscillation
       * would erase it. */
      int ns = c.ss > 0 ? c.ss : 1;
      for (int jj = j0; jj < j1; jj++)
        for (int ii = i0; ii < i1x; ii++) {
          double xc = -L + ((double)ii + 0.5) * px, yc = -L + ((double)jj + 0.5) * px;
          if (planecut && ((SU->side == 0) != (sdist(&c, xc, yc) < 0.0))) continue;
          double acc = 0.0;
          double complex ctr = 0.0;
          for (int sy = 0; sy < ns; sy++)
            for (int sx = 0; sx < ns; sx++) {
              double x = -L + ((double)ii + ((double)sx + 0.5) / (double)ns) * px;
              double y = -L + ((double)jj + ((double)sy + 0.5) / (double)ns) * px;
              double complex got = 0.0;
              for (int d = 0; d < SU->nd; d++)
                got += xs[SU->base + d] *
                       cexp(i1 * (SU->kx[d] * (x - SU->cx) + SU->ky[d] * (y - SU->cy)));
              acc += cabs(got) * cabs(got);
              if (sx == ns / 2 && sy == ns / 2) ctr = got;
            }
          size_t p = (size_t)(N - 1 - jj) * (size_t)N + (size_t)ii;
          fre[p] = (float)creal(ctr);
          fab[p] = (float)sqrt(acc / (double)(ns * ns));
        }
    }
    t_img = now_s() - t_img0;
    double amax = 0.0;
    for (size_t p = 0; p < (size_t)N * (size_t)N; p++)
      if ((double)fab[p] > amax) amax = (double)fab[p];
    if (!(amax > 0.0)) amax = 1.0;
    char nm[64];
    for (int which = 0; which < 2; which++) {
      snprintf(nm, sizeof nm, "img/tdg_%s_%d%d.ppm", which ? "re" : "abs", c.mirror, c.obj);
      FILE *f = fopen(nm, "wb");
      if (!f) continue;
      fprintf(f, "P6\n%d %d\n255\n", N, N);
      unsigned char *row = malloc(3u * (size_t)N);
      if (!row) {
        fclose(f);
        continue;
      }
      for (int jj = 0; jj < N; jj++) {
        for (int ii = 0; ii < N; ii++) {
          size_t p = (size_t)jj * (size_t)N + (size_t)ii;
          if (which) { /* Re u: blue - white - red, so the phase is visible */
            double v = (double)fre[p] / amax;
            double a = v > 0.0 ? v : -v;
            if (a > 1.0) a = 1.0;
            double lo = 1.0 - a;
            row[3 * ii + 0] = (unsigned char)(255.0 * (v > 0.0 ? 1.0 : lo));
            row[3 * ii + 1] = (unsigned char)(255.0 * lo);
            row[3 * ii + 2] = (unsigned char)(255.0 * (v > 0.0 ? lo : 1.0));
          } else { /* |u|^2, what a sensor integrates, on a display gamma */
            double v = (double)fab[p] / amax;
            v = v * v;
            if (v > 1.0) v = 1.0;
            unsigned char g = (unsigned char)(255.0 * pow(v, 1.0 / 2.2));
            row[3 * ii + 0] = g;
            row[3 * ii + 1] = g;
            row[3 * ii + 2] = g;
            /* THE GEOMETRY DRAWN IN. Without it the viewer cannot tell object
             * from interference, which is precisely how the first pictures
             * failed: green = the material boundary, yellow = the incoming
             * light, drawn along its direction of travel. */
            if (c.obj) {
              double xw = -L + ((double)ii + 0.5) * px;
              double yw = -L + ((double)(N - 1 - jj) + 0.5) * px;
              int in0 = inside_obj(&c, xw, yw);
              if (in0 != inside_obj(&c, xw + px, yw) || in0 != inside_obj(&c, xw, yw + px)) {
                row[3 * ii + 0] = 0;
                row[3 * ii + 1] = 255;
                row[3 * ii + 2] = 0;
              }
              double dxl = cos(c.theta), dyl = sin(c.theta);
              double ax = -0.92 * L * dxl, ay = -0.92 * L * dyl;
              double tt = (xw - ax) * dxl + (yw - ay) * dyl;
              if (tt > 0.0 && tt < 0.35 * L) {
                double dperp = fabs(-(xw - ax) * dyl + (yw - ay) * dxl);
                if (dperp < 1.5 * px) {
                  row[3 * ii + 0] = 255;
                  row[3 * ii + 1] = 220;
                  row[3 * ii + 2] = 0;
                }
              }
            }
          }
        }
        fwrite(row, 3u, (size_t)N, f);
      }
      free(row);
      fclose(f);
      printf("  wrote %s\n", nm);
    }
    free(fre);
    free(fab);
  }
  printf("  [time] assemble %.3f s   solve %.3f s (%d iters, %.2f ms/iter)   render %.3f s", t_asm,
         t_solve, itdone, 1e3 * t_solve / (double)(itdone > 0 ? itdone : 1), t_img);
  if (c.img > 0) printf("  = %.0f ns/pixel", 1e9 * t_img / ((double)c.img * (double)c.img));
  printf("\n");

  free(ja);
  free(va);
  free(rp);
  free(rhs);
  free(xs);
  free(uu);
  free(vv);
  free(ww);
  free(tv);
  return 0;
}
