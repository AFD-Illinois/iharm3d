
/* M1: no change */

/**

 advance particle positions using half-step primitives 

 CFG 7/2/12: interpolation may be half-zone off.  Needs testing.

**/

#include "decs.h"

void advance_particles(grid_prim_type ph, double Dt)
{
	int i, j, k, l, d;
	double f1, f2, f3, ucon[NDIM], uconll[NDIM], uconlh[NDIM],
	    uconhl[NDIM], uconhh[NDIM];
	struct of_geom *geom;

	for (l = 0; l < NPTOT; l++) {
		/* linear interpolate half-step velocity to position of particle */
		f1 = (xp[l][1] - startx[1]) / dx[1];
		f2 = (xp[l][2] - startx[2]) / dx[2];
		f3 = (xp[l][3] - startx[3]) / dx[3];
		i = (int) f1;	/* location in integral zones */
		j = (int) f2;
		k = (int) f3;
		f1 -= i;	/* location in fractional zones */
		f2 -= j;
		f3 -= k;

		/* watch for particles that leave the grid...
		   don't update them. */
		if (i >= 0 && i < N1 && j >= 0 && j < N2) {

			geom = get_geometry(i, j, CORN);
			ucon_calc(ph[i][j][k], geom, uconll);
			geom = get_geometry(i + 1, j, CORN);
			ucon_calc(ph[i + 1][j][k], geom, uconhl);
			geom = get_geometry(i, j + 1, CORN);
			ucon_calc(ph[i][j + 1][k], geom, uconlh);
			geom = get_geometry(i + 1, j + 1, CORN);
			ucon_calc(ph[i + 1][j + 1][k], geom, uconhh);

			for (d = 0; d < NDIM; d++) {
				ucon[d] =
				    f1 * f2 * uconll[d] + (1. -
							   f1) * f2 *
				    uconhl[d]
				    + f1 * (1. - f2) * uconlh[d] + (1. -
								    f1) *
				    (1. - f2) * uconhh[d];
			}

			/* step particle position forward */
			for (d = 1; d < NDIM; d++)
				xp[l][d] += Dt * (ucon[d]/ucon[0]);

		}

	}
	/* done! */
}
