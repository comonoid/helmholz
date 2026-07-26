#include "cut2d.h"
#include "phi.h"
#include <complex.h>
#include <math.h>

/* Degrees: phi*phi is 4 per axis; where the oscillation is absorbed (see ABSORB)
 * the primitive axis gains ABSORB_TERMS-1 = 18 and the antiderivative one more,
 * so the worst product Qx*Fy is 23 + 4 = 27. 30 leaves head-room. */
enum { MAXDEG = 30, MAXPT = 12, MAXV = HZ_CUT_MAXHP + 8, ABSORB_TERMS = 19 };

/* |omega| * (piece half-width) below which exp(i*omega*s) is absorbed into the
 * polynomial instead of being carried.
 * WHY THIS EXISTS AND WHY THE THRESHOLD IS 1 AND NOT SMALL. The antiderivative
 * recursion divides by i*omega, so its top coefficient reaches deg!/omega^(deg+1)
 * while the piece integral is only ~q*w: the cancellation is deg!/(omega w)^deg,
 * i.e. 1e7 at omega*w = 0.05. The first run of test_cut2d showed exactly that —
 * 7.6e-10 on the om=0 row with a threshold of 1e-2. Absorbing instead costs a
 * truncation of THRESH^M/M!, which at THRESH = 1 and M = 19 is 8e-18.
 * The decision is PER PIECE, and only on the axis the primitive is taken along:
 * the other axis never divides by its omega, so it has nothing to lose. */
static const double ABSORB = 1.0;
/* Series/recurrence seam for the moments. Series cancellation is bounded by the
 * peak term e^{|a|} (3e3 at 8, i.e. under four digits); the upward recurrence is
 * stable once |a| >= n, and n <= 8 whenever the exponential is NOT absorbed. */
static const double MOM_SERIES_MAX = 8.0;

/* ------------------------------------------------------------------ moments */
/* mom[n] = integral over [-1,1] of u^n exp(i a u) du, n = 0..N. */
static void osc_moments_n(double complex a, int N, double complex *mom) {
  double aa = cabs(a);
  double lim = MOM_SERIES_MAX > (double)N ? MOM_SERIES_MAX : (double)N;
  if (aa <= lim) {
    for (int n = 0; n <= N; n++)
      mom[n] = 0.0;
    double complex coef = 1.0; /* (i a)^m / m! */
    for (int m = 0; m < 80; m++) {
      if (m > 0) coef *= CMPLX(0.0, 1.0) * a / (double)m;
      for (int n = m & 1; n <= N; n += 2)
        mom[n] += 2.0 * coef / (double)(n + m + 1);
      if ((double)m > aa && cabs(coef) < 1e-20) break;
    }
    return;
  }
  double complex ia = CMPLX(0.0, 1.0) * a;
  double complex ep = cexp(CMPLX(0.0, 1.0) * a), em = cexp(CMPLX(0.0, -1.0) * a);
  mom[0] = (ep - em) / ia;
  for (int n = 1; n <= N; n++) {
    double complex sgn = (n & 1) ? -1.0 : 1.0;
    mom[n] = (ep - sgn * em) / ia - (double)n * mom[n - 1] / ia;
  }
}

void hz_cut2d_osc_moments(double complex a, int N, double complex *mom) {
  if (N < 0 || N > MAXDEG) return;
  osc_moments_n(a, N, mom);
}

/* --------------------------------------------------------------- polynomials */
/* d(t) = c(a + b t), degree preserved. */
static void poly_compose(const double complex *c, int deg, double complex a, double complex b,
                         double complex *d) {
  for (int i = 0; i <= deg; i++)
    d[i] = 0.0;
  for (int k = deg; k >= 0; k--) {
    for (int i = deg; i >= 1; i--)
      d[i] = d[i] * a + d[i - 1] * b;
    d[0] = d[0] * a + c[k];
  }
}

static int poly_mul(const double complex *a, int da, const double complex *b, int db,
                    double complex *out) {
  for (int i = 0; i <= da + db; i++)
    out[i] = 0.0;
  for (int i = 0; i <= da; i++)
    for (int j = 0; j <= db; j++)
      out[i + j] += a[i] * b[j];
  return da + db;
}

/* Integral over [-1,1] of P(u) exp(i gam u).
 * The caller composes into this unit interval rather than passing [t0,t1]: with
 * the absorbed branch the polynomial reaches degree ~27, and composing from a
 * coordinate whose origin sits a few units away raises intermediate
 * coefficients by offset^27 and loses them again to cancellation. Measured: the
 * V3 additivity rows read 1.1e-13 that way and 1e-15 after the change. */
static double complex poly_exp_unit(const double complex *c, int deg, double complex gam) {
  double complex mom[MAXDEG + 1];
  osc_moments_n(gam, deg, mom);
  double complex acc = 0.0;
  for (int n = 0; n <= deg; n++)
    acc += c[n] * mom[n];
  return acc;
}

/* Q with Q' + i*om*Q = q (degree preserved), or the plain antiderivative with
 * Q(0) = 0 when the oscillation was absorbed (degree +1). */
static int antideriv(const double complex *q, int deg, double complex om, int absorbed,
                     double complex *Q) {
  if (absorbed) {
    Q[0] = 0.0;
    for (int k = 0; k <= deg; k++)
      Q[k + 1] = q[k] / (double)(k + 1);
    return deg + 1;
  }
  double complex io = CMPLX(0.0, 1.0) * om;
  Q[deg] = q[deg] / io;
  for (int k = deg - 1; k >= 0; k--)
    Q[k] = (q[k] - (double)(k + 1) * Q[k + 1]) / io;
  return deg;
}

/* ---------------------------------------------------------------- one axis */
typedef struct {
  double pt[MAXPT]; /* breakpoints, npt of them => npt-1 pieces */
  int npt;
  double mid[MAXPT];                   /* piece midpoints */
  double complex q[MAXPT][MAXDEG + 1]; /* factor polynomial in s = x - mid */
  int deg[MAXPT];                      /* its degree (higher where absorbed) */
  int absorb[MAXPT];                   /* exp(i om s) folded into q[] */
} axis_pieces;

static void axis_clip(hz_axis2 f, double *lo, double *hi) {
  for (int i = 0; i < f.nf; i++) {
    double s0 = f.f[i].h * (f.f[i].n - 2.0), s1 = f.f[i].h * (f.f[i].n + 2.0);
    if (s0 > *lo) *lo = s0;
    if (s1 < *hi) *hi = s1;
  }
}

/* Breakpoints and per-piece polynomials. Where allow_absorb is set and the
 * piece is short enough (see ABSORB), exp(i om (x-xref)) is folded into the
 * polynomial as its truncated series and the piece is marked. */
static void axis_build(hz_axis2 f, double lo, double hi, double complex om, int allow_absorb,
                       double xref, axis_pieces *ap) {
  static const double knot_t[4] = {-2.0, -1.0, 1.0, 2.0};
  int np = 0;
  ap->pt[np++] = lo;
  ap->pt[np++] = hi;
  for (int i = 0; i < f.nf; i++)
    for (int k = 0; k < 4; k++) {
      double x = f.f[i].h * (f.f[i].n + knot_t[k]);
      if (x > lo && x < hi && np < MAXPT) ap->pt[np++] = x;
    }
  for (int i = 1; i < np; i++) { /* insertion sort */
    double v = ap->pt[i];
    int j = i - 1;
    while (j >= 0 && ap->pt[j] > v) {
      ap->pt[j + 1] = ap->pt[j];
      j--;
    }
    ap->pt[j + 1] = v;
  }
  ap->npt = np;
  for (int p = 0; p + 1 < np; p++) {
    double m = 0.5 * (ap->pt[p] + ap->pt[p + 1]);
    double w = 0.5 * (ap->pt[p + 1] - ap->pt[p]);
    int absorbed = allow_absorb && cabs(om) * w < ABSORB;
    /* Only as many series terms as the piece actually needs: (om w)^M/M! < 1e-17.
     * Equal carriers give om == 0 exactly, and those pairs are the MAJORITY of a
     * Gram matrix — carrying 19 terms for them made the bench 5x slower for
     * nothing. */
    int nser = 1;
    if (absorbed) {
      double z = cabs(om) * w, t = 1.0;
      while (nser < ABSORB_TERMS && t > 1e-17) {
        t *= z / (double)nser;
        nser++;
      }
    }
    ap->mid[p] = m;
    ap->absorb[p] = absorbed;
    ap->deg[p] = 2 * f.nf + (absorbed ? nser - 1 : 0);
    double complex base[MAXDEG + 1];
    for (int i = 0; i <= ap->deg[p]; i++)
      base[i] = 0.0;
    base[0] = 1.0;
    int db = 0;
    for (int i = 0; i < f.nf; i++) {
      double qq[3];
      hz_phi_local_poly(f.f[i], m, qq);
      double complex fac[3] = {qq[0], qq[1], qq[2]}, tmp[MAXDEG + 1];
      db = poly_mul(base, db, fac, 2, tmp);
      for (int k = 0; k <= db; k++)
        base[k] = tmp[k];
    }
    if (absorbed) {
      /* exp(i om x) = exp(i om (x-xref)) * exp(i om xref); the second factor is
       * pulled out of the whole integral by the caller. */
      double complex ser[ABSORB_TERMS] = {0},
                     tmp[MAXDEG + 1]; /* zeroed: cppcheck cannot see nser >= 1 */
      double complex c = cexp(CMPLX(0.0, 1.0) * om * (m - xref));
      double complex t = 1.0;
      for (int k = 0; k < nser; k++) {
        ser[k] = c * t;
        t *= CMPLX(0.0, 1.0) * om / (double)(k + 1);
      }
      db = poly_mul(base, db, ser, nser - 1, tmp);
      for (int k = 0; k <= db; k++)
        base[k] = tmp[k];
    }
    for (int k = 0; k <= ap->deg[p]; k++)
      ap->q[p][k] = k <= db ? base[k] : 0.0;
  }
}

static int piece_of(const axis_pieces *ap, double x) {
  for (int p = 0; p + 2 < ap->npt; p++)
    if (x < ap->pt[p + 1]) return p;
  return ap->npt - 2;
}

/* ------------------------------------------------------------------ polygon */
static int clip_half(const double *vx, const double *vy, int nv, double ca, double sa, double c,
                     int keep_le, double *ox, double *oy) {
  int no = 0;
  for (int i = 0; i < nv; i++) {
    int j = (i + 1) % nv;
    double si = vx[i] * ca + vy[i] * sa - c, sj = vx[j] * ca + vy[j] * sa - c;
    double di = keep_le ? -si : si, dj = keep_le ? -sj : sj; /* >=0 means inside */
    if (di >= 0.0) {
      ox[no] = vx[i];
      oy[no++] = vy[i];
    }
    if ((di >= 0.0) != (dj >= 0.0)) {
      double t = di / (di - dj);
      ox[no] = vx[i] + t * (vx[j] - vx[i]);
      oy[no++] = vy[i] + t * (vy[j] - vy[i]);
    }
    if (no >= MAXV) break;
  }
  return no;
}

/* --------------------------------------------------------------------- core */
/* Primitive taken along x: F = (W(x) Fy(y) e^{i omy y}, 0), so only edges with
 * dy != 0 contribute and horizontal edges drop out. */
static double complex core(double xlo, double xhi, double ylo, double yhi, hz_axis2 fx, hz_axis2 fy,
                           double complex omx, double complex omy, const hz_half2 *hp, int nhp) {
  double xref = 0.5 * (xlo + xhi), yref = 0.5 * (ylo + yhi);
  axis_pieces px, py;
  axis_build(fx, xlo, xhi, omx, 1, xref, &px);
  axis_build(fy, ylo, yhi, omy, 0, yref, &py);

  /* antiderivatives along x, made continuous piece to piece */
  double complex Q[MAXPT][MAXDEG + 1], C[MAXPT];
  int degQ[MAXPT];
  double complex val = 0.0; /* W at the left end of the current piece */
  for (int p = 0; p + 1 < px.npt; p++) {
    int ab = px.absorb[p];
    degQ[p] = antideriv(px.q[p], px.deg[p], omx, ab, Q[p]);
    double sL = px.pt[p] - px.mid[p], sR = px.pt[p + 1] - px.mid[p];
    double complex eL = ab ? 1.0 : cexp(CMPLX(0.0, 1.0) * omx * (px.pt[p] - xref));
    double complex eR = ab ? 1.0 : cexp(CMPLX(0.0, 1.0) * omx * (px.pt[p + 1] - xref));
    double complex qL = 0.0, qR = 0.0;
    for (int k = degQ[p]; k >= 0; k--) {
      qL = qL * sL + Q[p][k];
      qR = qR * sR + Q[p][k];
    }
    C[p] = val - qL * eL;
    val = C[p] + qR * eR;
  }

  double vx[MAXV], vy[MAXV], wx[MAXV], wy[MAXV];
  vx[0] = xlo, vy[0] = ylo;
  vx[1] = xhi, vy[1] = ylo;
  vx[2] = xhi, vy[2] = yhi;
  vx[3] = xlo, vy[3] = yhi;
  int nv = 4;
  for (int hpi = 0; hpi < nhp && nv >= 3; hpi++) {
    nv = clip_half(vx, vy, nv, hp[hpi].ca, hp[hpi].sa, hp[hpi].c, hp[hpi].keep_le, wx, wy);
    for (int i = 0; i < nv; i++)
      vx[i] = wx[i], vy[i] = wy[i];
  }
  if (nv < 3) return 0.0;
  double area2 = 0.0;
  for (int i = 0; i < nv; i++) {
    int j = (i + 1) % nv;
    area2 += vx[i] * vy[j] - vx[j] * vy[i];
  }
  if (area2 < 0.0) { /* the divergence theorem wants counter-clockwise */
    for (int i = 0; i < nv / 2; i++) {
      double tx = vx[i], ty = vy[i];
      vx[i] = vx[nv - 1 - i], vy[i] = vy[nv - 1 - i];
      vx[nv - 1 - i] = tx, vy[nv - 1 - i] = ty;
    }
  }

  double complex total = 0.0;
  for (int e = 0; e < nv; e++) {
    double x0 = vx[e], y0 = vy[e];
    double dx = vx[(e + 1) % nv] - x0, dy = vy[(e + 1) % nv] - y0;
    if (!(fabs(dy) > 0.0)) continue; /* n_x dl = dy: horizontal edges drop out */
    double ts[2 * MAXPT + 2];
    int nt = 0;
    ts[nt++] = 0.0;
    ts[nt++] = 1.0;
    if (fabs(dx) > 0.0)
      for (int p = 1; p + 1 < px.npt; p++) {
        double t = (px.pt[p] - x0) / dx;
        if (t > 0.0 && t < 1.0) ts[nt++] = t;
      }
    for (int p = 1; p + 1 < py.npt; p++) {
      double t = (py.pt[p] - y0) / dy;
      if (t > 0.0 && t < 1.0) ts[nt++] = t;
    }
    for (int i = 1; i < nt; i++) {
      double v = ts[i];
      int j = i - 1;
      while (j >= 0 && ts[j] > v) {
        ts[j + 1] = ts[j];
        j--;
      }
      ts[j + 1] = v;
    }
    for (int s = 0; s + 1 < nt; s++) {
      double ta = ts[s], tb = ts[s + 1];
      if (tb - ta <= 0.0) continue;
      double tc = 0.5 * (ta + tb), hw = 0.5 * (tb - ta);
      double xm = x0 + dx * tc, ym = y0 + dy * tc;
      int ip = piece_of(&px, xm), jp = piece_of(&py, ym);
      double complex Py[MAXDEG + 1], Qx[MAXDEG + 1], prod[MAXDEG + 1];
      double complex t1[MAXDEG + 1] = {0}; /* zeroed: cppcheck cannot see deg >= 0 */
      int ab = px.absorb[ip];
      /* everything is composed into u in [-1,1] about the SUB-EDGE midpoint, so
       * every offset stays inside one piece (see poly_exp_unit) */
      poly_compose(py.q[jp], py.deg[jp], ym - py.mid[jp], dy * hw, Py);
      poly_compose(Q[ip], degQ[ip], xm - px.mid[ip], dx * hw, Qx);
      int dp = poly_mul(Qx, degQ[ip], Py, py.deg[jp], prod);
      for (int k = 0; k <= py.deg[jp]; k++)
        t1[k] = C[ip] * Py[k];

      /* the y oscillation is always carried; the x one only where the piece did
       * not absorb it */
      double complex gy = omy * dy * hw;
      double complex cy = cexp(CMPLX(0.0, 1.0) * omy * (ym - yref));
      double complex g2 = ab ? gy : (omx * dx + omy * dy) * hw;
      double complex c2 = ab ? cy : cy * cexp(CMPLX(0.0, 1.0) * omx * (xm - xref));
      total +=
          dy * hw * (cy * poly_exp_unit(t1, py.deg[jp], gy) + c2 * poly_exp_unit(prod, dp, g2));
    }
  }
  return total * cexp(CMPLX(0.0, 1.0) * (omx * xref + omy * yref));
}

double complex hz_cut2d_poly(double xlo, double xhi, double ylo, double yhi, hz_axis2 fx,
                             hz_axis2 fy, double complex omx, double complex omy,
                             const hz_half2 *hp, int nhp) {
  axis_clip(fx, &xlo, &xhi);
  axis_clip(fy, &ylo, &yhi);
  if (!(xlo < xhi) || !(ylo < yhi) || nhp > HZ_CUT_MAXHP) return 0.0;
  /* Take the primitive along the axis with the larger |omega|: the recursion in
   * antideriv() divides by i*omega, so the weaker axis is the wrong one. */
  if (cabs(omy) > cabs(omx)) {
    hz_half2 sw[HZ_CUT_MAXHP];
    for (int i = 0; i < nhp; i++) {
      sw[i].ca = hp[i].sa;
      sw[i].sa = hp[i].ca;
      sw[i].c = hp[i].c;
      sw[i].keep_le = hp[i].keep_le;
    }
    return core(ylo, yhi, xlo, xhi, fy, fx, omy, omx, sw, nhp);
  }
  return core(xlo, xhi, ylo, yhi, fx, fy, omx, omy, hp, nhp);
}

double complex hz_cut2d_integral(double xlo, double xhi, double ylo, double yhi, hz_axis2 fx,
                                 hz_axis2 fy, double complex omx, double complex omy,
                                 hz_strip2 st) {
  hz_half2 hp[2];
  int nhp = 0;
  if (!(st.slo < st.shi)) return 0.0;
  if (st.slo > -HZ_CUT_INF) {
    hp[nhp].ca = st.ca;
    hp[nhp].sa = st.sa;
    hp[nhp].c = st.slo;
    hp[nhp].keep_le = 0;
    nhp++;
  }
  if (st.shi < HZ_CUT_INF) {
    hp[nhp].ca = st.ca;
    hp[nhp].sa = st.sa;
    hp[nhp].c = st.shi;
    hp[nhp].keep_le = 1;
    nhp++;
  }
  return hz_cut2d_poly(xlo, xhi, ylo, yhi, fx, fy, omx, omy, hp, nhp);
}

int hz_cut2d_disc(double cx, double cy, double r, int m, int fit, hz_half2 *hp) {
  if (m < 3) m = 3;
  if (m > HZ_CUT_MAXHP) m = HZ_CUT_MAXHP;
  /* Facet i is the line normal to angle th_i at distance d from the centre.
   * Inscribed takes d = r cos(pi/m), the chord through two vertices ON the
   * circle, so the polygon is entirely inside and the displacement is
   * one-signed; mid-fit splits the sagitta so the polygon crosses the circle
   * and the displacement has zero mean. M14 [3] measured the error law for a
   * UNIFORM (one-signed) displacement, so inscribed is the case that law
   * covers and mid-fit is the one that should beat it. */
  double c = cos(M_PI / (double)m);
  double d = fit ? 0.5 * r * (1.0 + c) : r * c;
  for (int i = 0; i < m; i++) {
    double th = 2.0 * M_PI * ((double)i + 0.5) / (double)m;
    hp[i].ca = cos(th);
    hp[i].sa = sin(th);
    hp[i].c = d + cx * cos(th) + cy * sin(th);
    hp[i].keep_le = 1;
  }
  return m;
}

int hz_cut2d_disc_verts(double cx, double cy, double r, int m, int fit, double *vx, double *vy) {
  if (m < 3) m = 3;
  if (m > HZ_CUT_MAXHP) m = HZ_CUT_MAXHP;
  /* Facet normals sit at (i+0.5)*2pi/m, so two consecutive facets meet at
   * i*2pi/m; the distance d of the facet from the centre is the vertex radius
   * times cos(pi/m). Same d as hz_cut2d_disc, taken from the same expression. */
  double c = cos(M_PI / (double)m);
  double d = fit ? 0.5 * r * (1.0 + c) : r * c;
  double rv = d / c;
  for (int i = 0; i < m; i++) {
    double th = 2.0 * M_PI * (double)i / (double)m;
    vx[i] = cx + rv * cos(th);
    vy[i] = cy + rv * sin(th);
  }
  return m;
}
