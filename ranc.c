
#include "decs.h"

void init_ranc(int seed)
{
        r = gsl_rng_alloc(gsl_rng_mt19937);     /* use Mersenne twister */
        gsl_rng_set(r, seed);
}

/* return pseudo-random value between 0 and 1 */
double ranc()
{
        return (gsl_rng_uniform(r));
}

