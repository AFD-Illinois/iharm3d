/******************************************************************************  
 *                                                                            *  
 * RANDOM.C                                                                   *  
 *                                                                            *  
 * WRAPPER FOR RANDOM NUMBER GENERATOR                                        *  
 *                                                                            *  
 ******************************************************************************/

#include "decs.h"

void init_random(int seed)
{
  // Use Mersenne twister
  r = gsl_rng_alloc(gsl_rng_mt19937);
  gsl_rng_set(r, seed);
}

// Return pseudorandom value between 0 and 1
double get_random()
{
  return (gsl_rng_uniform(r));
}

