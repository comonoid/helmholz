#ifndef HZ_NITSCHE2D_H
#define HZ_NITSCHE2D_H

#include "carrier2d.h"
#include <complex.h>

/* THE THREE INTERFACE INTEGRALS A CUT COSTS IN A GALERKIN FORM.
 *
 * WHY THEY EXIST AT ALL. A cut basis function has a JUMP across the cut locus.
 * The strong form Int B_i (Lap + k^2) B_j integrates each side SEPARATELY and
 * therefore silently drops the distributional terms that live on the locus (a
 * delta and a delta-prime): the discrete operator does not see the jump at all,
 * the jump is free, and every combination a*B_in + b*B_out with a != b sits in
 * the kernel WITH a non-zero field. Measured (vacuum, slice step 3): 0.197..0.405
 * with the cut on against 0.0073..0.0166 with it off.
 *
 * WHAT REPAIRS IT, DERIVED RATHER THAN QUOTED. Supports are compact and phi is
 * C^1, so Green's identity over the two regions is exact and gives
 *
 *   Int B_i (Lap + k^2) B_j = -a_broken(B_i,B_j)
 *                             - Int_G ( [B_i]{dn B_j} + {B_i}[dn B_j] )
 *   a_broken(u,v) = Int (grad u . grad v - k^2 u v)          (symmetric)
 *
 * i.e. THE STRONG FORM ALREADY CONTAINS THE WEAK ONE, integrated by parts
 * analytically, and only two surface terms are missing. Adding them and the
 * standard Nitsche pair turns the assembled matrix into (minus) the symmetric
 * Nitsche form for the transmission problem [u] = 0, [dn u] = 0:
 *
 *   a_N(u,v) = a_broken(u,v) + Int_G ({dn u}[v] + {dn v}[u]) + beta Int_G [u][v]
 *
 * WITH THAT SIGN, AND IT IS NOT A CONVENTION. The right-hand side is untouched
 * only if a_N vanishes on the exact solution, and it is NOT enough that every
 * added term contain a jump: for the exact solution [u] = [dn u] = 0, but the
 * TEST function still jumps, and Green's identity leaves
 * a_broken(u,v) = -Int_G [v]{dn u} behind. The interface pair must CANCEL that,
 * not double it. Written with the other sign the assembled operator missed the
 * known plane wave by |Ac| = 0.40 while still coming out symmetric to 8e-16 —
 * the symmetry check cannot see this, because flipping a symmetric pair leaves
 * a symmetric matrix.
 * Unlike point constraints there is no sampling parameter left; beta is the only
 * coefficient and it is bounded below by an inverse inequality (M15 audit §2).
 *
 * BILINEAR, NOT SESQUILINEAR - the test function is not conjugated, exactly as
 * in the volume assembly, so the two exponents simply ADD.
 *
 * The 3D case is the same object with the line integral replaced by an integral
 * over a polygonal face; not implemented here. */

typedef struct {
  double complex t0;  /* Int_G B_i B_j ds */
  double complex t1j; /* Int_G B_i (n . grad B_j) ds */
  double complex t1i; /* Int_G B_j (n . grad B_i) ds */
} hz_nit2d;

/* One straight facet, from (x0,y0) to (x1,y1), with the unit normal (nx,ny)
 * pointing the way the caller's jump convention counts as "plus". Exact to
 * rounding: the integrand is (piecewise polynomial of degree <= 8 in the arc
 * length) x (plane wave), split at every knot the four phi factors cross. */
hz_nit2d hz_nitsche2d_seg(hz_carrier2d bi, hz_carrier2d bj, double x0, double y0, double x1,
                          double y1, double nx, double ny);

/* Closed counter-clockwise polygon (the faceted material boundary). The OUTWARD
 * normal of edge e is (dy,-dx)/L, which for a counter-clockwise ring points
 * away from the interior - i.e. from the scatterer outwards. */
hz_nit2d hz_nitsche2d_poly(hz_carrier2d bi, hz_carrier2d bj, const double *vx, const double *vy,
                           int nv);

#endif
