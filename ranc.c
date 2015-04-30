
#include <stdint.h>

static uint64_t xorx;
static uint64_t xors[16];

static int p_;

uint64_t xorshift_init(void)
{
	xorx ^= xorx >> 12; // a
	xorx ^= xorx << 25; // b
	xorx ^= xorx >> 27; // c
	return xorx * UINT64_C(2685821657736338717);
}

void init_ranc(int seed)
{
	xorx = seed;
	for(int i=0;i<16;i++) xors[i] = xorshift_init();
}

double ranc()
{
	uint64_t s0 = xors[ p_ ];
	uint64_t s1 = xors[ p_ = ( p_ + 1 ) & 15 ];
	s1 ^= s1 << 31;
	s1 ^= s1 >> 11;
	s0 ^= s0 >> 30;
	return (( xors[ p_ ] = s0 ^ s1 ) * UINT64_C(1181783497276652981))/((double)UINT64_MAX);
}

/*
void init_ranc(int seed)
{
        r = gsl_rng_alloc(gsl_rng_mt19937);     // use Mersenne twister
        gsl_rng_set(r, seed);
}

// return pseudo-random value between 0 and 1
double ranc()
{
        return (gsl_rng_uniform(r));
}
*/
