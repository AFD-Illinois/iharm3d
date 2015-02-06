
/* M1: no changes required */

/*
 *
 * initialize Lagrangian tracer particles
 *
 * cfg 4 feb 09
 *
 */

#include "decs.h"

#define SAMPLE_FACTOR	4.

void init_particles()
{
	int i, j, k, Np;
	double Nt, Nexp, sample_factor, X[NDIM];

	/* global variables */
	pdump_cnt = 0;
	DTp = 100.;

	/* assign particles according to restmass density */

	/* first total up a quantity that is proportional
	   to rest mass density */
	Nt = 0.;
	ZLOOP {
		Nt += ggeom[i][j][CENT].g * p[i][j][k][RHO];
	}

	/* normalization factor so that we get the # of
	   particles we want */
	sample_factor = NPTOT / Nt;

	/* now sample, keeping running total of particles
	   created */
	Np = 0;
	Nexp = 0.;
	ZLOOP {
		Nexp += ggeom[i][j][CENT].g * p[i][j][k][RHO] * sample_factor;

		/* assign particle to random position in cell */
		while (Nexp >= 1.) {

			/* lower left corner of zone is here */
			coord(i, j, CORN, X);
			X[3] = startx[3] + (k-START3 + global_start[2])*dx[3];

			xp[Np][0] = 0.;
			xp[Np][1] = X[1] + ranc() * dx[1];
			xp[Np][2] = X[2] + ranc() * dx[2];
			xp[Np][3] = X[3] + ranc() * dx[3];

			Np++;
			Nexp -= 1.;
		}
	}

	/* report back! */
	fprintf(stderr, "made %d particles of %d\n", Np, NPTOT);

	/* done! */
}
