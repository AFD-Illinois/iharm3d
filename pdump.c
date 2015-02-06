
/*
	M1: no changes required 
*/

#include "decs.h"

void pdump(FILE * fp)
{
	int l;
	double r, th;

	for (l = 0; l < NPTOT; l++) {
		bl_coord(xp[l], &r, &th);
		fprintf(fp, "%d %g %g %g %g\n", l, xp[l][1], xp[l][2], r,
			th);
	}
}
