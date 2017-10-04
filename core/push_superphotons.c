/******************************************************************************
 *                                                                            *
 * PUSH_SUPERPHOTONS.C                                                        *
 *                                                                            *
 * INTEGRATE SUPERPHOTON GEODESICS                                            *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
int push(double t, double dt, double tmpX[NDIM], tmpK[NDIM], 
  struct of_photon *ph)

// Don't do work when sources are trivial
#if METRIC == MINKOWSKI
const int killing[] = {1, 1, 1, 1};
#elif METRIC == MKS
const int killing[] = {1, 0, 0, 1};
#endif

#define MAX_SUBDIV (7)
void push_superphotons(double t, double dt)
{
  #pragma omp parallel
  {
    struct of_photon *ph = photon_lists[omp_get_thread_num()];
    struct of_photon *prev = NULL;
    struct of_photon *head = ph;
    //double tmpX[NDIM], tmpKcov[NDIM];
    int failure, nsteps;

    // Push each photon from n to n+1
    while (ph != NULL)
    {
      failure = 1;
      nsteps = 0;
      
      // Store X and K at step n
      for (int mu = 0; mu < NDIM; mu++) {
        ph->Xprev[mu] = ph->X[mu];
        ph->Kcovprev[mu] = ph->Kcov[mu];
      }

      // If integration fails, repeat with progressively smaller timesteps
      while (failure == 1 && nsteps < MAX_SUBDIV) {
        for (int mu = 0; mu < NDIM; mu++) {
          ph->X[mu] = ph->Xprev[mu];
          ph->Kcov[mu] = ph->Kcovprev[mu];
          //tmpX[mu] = ph->X[mu];
          //tmpKcov[mu] = ph->K[mu];
        }
        
        for (int n = 0; n <= nsteps; n++) {
          failure = push(dt/pow(2,nsteps), ph);
        }
        nsteps++;
      }
      
      printf("ph->X[1] = %e\n", ph->X[1]);

      if (failure == 1) {
        list_remove(&ph, &head, &prev);
        #pragma omp atomic
        step_lost++;
        continue;
      } // failure

      prev = ph;
      ph = ph->next;

      #pragma omp atomic
      step_tot++;
    } // ph != NULL
  } // omp parallel
}
#undef MAX_SUBDIV

// Second order update to X^{\mu}, K_{\mu} from t to t + dt
int push(double dt, struct of_photon *ph)
{
  /* Heun's method:
   *   x_n+1 = x_n + dt*(0.5*c1 + 0.5*c2)
   *     c1  = dxdt(t_n, x_n)
   *     c2  = dxdt(t_n + 0.5*dt, x_n + 0.5*dt*c1)
   *   y_n+1 = y_n + dt*(0.5*d1 + 0.5*d2)
   *     d1  = dydt(t_n, y_n)
   *     d2  = dydt(t_n + 0.5*dt, y_n + 0.5*dt*d1)
   *
   *   where
   *     x = X^{\mu}, \mu = [1, 2, 3] (X[0] known)
   *     y = K_{\mu}, \mu = [1, 2]    (K_0, K_3 conserved)
   *     dydt = -1/(2*g^{0\nu}*k_{\nu})*k_b*k_c*(d g^{bc} / dx^{\mu})
   *     dxdt = k^{\mu} / k^0
   */

   double c1[NDIM], c2[NDIM], d1[NDIM], d2[NDIM];
   double gcov[NDIM][NDIM], gcon[NDIM][NDIM];
   double X[NDIM], Kcov[NDIM], Kcon[NDIM];

   // First stage
   memcpy(X,    ph->X,    NDIM*sizeof(double));
   memcpy(Kcov, ph->Kcov, NDIM*sizeof(double));

   gcov_func(X, gcov);
   gcon_func(gcov, gcon);
   for (int mu = 0; mu < NDIM; mu++) {
     Kcon[mu] = 0.;
     for (int nu = 0; nu < NDIM; nu++) {
       Kcon[mu] += gcon[mu][nu]*Kcov[nu]
     }
   }
   //for (int mu = 0; mu < NDIM; mu++) {
     
   //}

   get_X_source(Kcon, c1);
   get_K_source(X, Kcov, Kcon, d1);

    
}

void get_X_source(double Kcon[NDIM], double c[NDIM])
{
  for (int i = 1; i < NDIM; i++) {
    c1[i] = Kcon[i]/Kcon[0];
  }
}

void get_K_source(double X[NDIM], double Kcov[NDIM], double Kcon[NDIM], 
  double d[NDIM])
{
  for (int mu = 0; mu < NDIM; mu++) {
    if (killing[mu] == 1) {
      d[mu] = 0.;
    } else {

    }
  }
}
#endif // RADIATION
