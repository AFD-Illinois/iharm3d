/******************************************************************************
 *                                                                            *
 * MAKE_SUPER_PHOTONS.C                                                       *
 *                                                                            *
 * EMISSION OF MONTE CARLO SAMPLES                                            *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION

static double lnu_min, lnu_max, dlnu, nusamp[NU_BINS+1], Ns;
static double dnzs[N1+2*NG][N2+2*NG][N3+2*NG];

void sample_photon(int i, int j, int k, double t, double dt, 
  double dndlnu[NU_BINS+1], struct of_photon *tmp, double Econ[NDIM][NDIM], 
  double Ecov[NDIM][NDIM], double Ne, double Thetae, double Bmag, 
  double Bcon[NDIM]);
void get_dndlnu(int i, int j, int k, double dt, double dndlnu[NU_BINS+1], 
  double Ne, double Thetae, double Bmag);

void make_superphotons(grid_prim_type Prad, double t, double dt)
{
  get_dnz(Prad);

  #pragma omp parallel
  {
    struct of_photon *tmp, *head = photon_lists[omp_get_thread_num()];
    double dndlnu[NU_BINS+1];
    double Econ[NDIM][NDIM], Ecov[NDIM][NDIM];
    double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    double Ne, Thetae, Bmag;
    double X[NDIM];
    int nz;

    #pragma omp for collapse(3)
    ZLOOP
    {
      nz = (int)dnzs[i][j][k];
      if (dnzs[i][j][k] - nz > get_random()) nz++;
      
      if (nz > 0) {
        // Set up zone
        coord(i, j, k, CENT, X);
        get_fluid_zone(i, j, k, Prad, &Ne, &Thetae, &Bmag, Ucon, Ucov, Bcon, 
          Bcov);
        make_tetrad(i, j, k, Ucon, Bcon, ggeom[i][j][k].gcov, Econ, Ecov);
        get_dndlnu(i, j, k, dt, dndlnu, Ne, Thetae, Bmag);

        // Create superphotons in pairs
        for (int n = 0; n < nz; n++) {
          tmp = my_malloc(1, sizeof(struct of_photon));
          tmp->next = my_malloc(1, sizeof(struct of_photon));

          sample_photon(i, j, k, t, dt, dndlnu, tmp, Econ, Ecov, Ne,
            Thetae, Bmag, Bcon);

          (tmp->next)->next = head;
          head = tmp;
        } // n < nz

        #pragma omp atomic
        step_made += 2*nz;
      } // nz > 0
    } // ZLOOP

    // Prepend created superphotons to each thread's global list
    photon_lists[omp_get_thread_num()] = head;
  } // omp parallel
}
       
void sample_photon(int i, int j, int k, double t, double dt, 
  double dndlnu[NU_BINS+1], struct of_photon *ph, double Econ[NDIM][NDIM], 
  double Ecov[NDIM][NDIM], double Ne, double Thetae, double Bmag, 
  double Bcon[NDIM])
{
  double nu, th, cth[2], sth[2], phi, sphi[2], cphi[2];
  double Kcov_tetrad[NDIM];
  struct of_photon *tmp[2];
  tmp[0] = ph;
  tmp[1] = ph->next;
  tmp[1]->next = NULL;

  // Sample emissivity to get frequency
  do {
    nu = exp(get_random()*(lnu_max - lnu_min) + lnu_min);
  } while (get_random() > linear_interp_log(nu, dndlnu, lnu_min, dlnu));

  // Get weight from global weight parameter
  double weight = wgtC/nu;

  // Sample emissivity in solid angle
  double jmax = jnu(nu, Ne, Thetae, Bmag, 0.5*M_PI);
  do {
    cth[0] = 2.*get_random() - 1.;
    th = acos(cth[0]);
  } while (get_random() > jnu(nu, Ne, Thetae, Bmag, th)/jmax);
  
  sth[0] = sqrt(1. - cth[0]*cth[0]);
  phi = 2.*M_PI*get_random();
  cphi[0] = cos(phi);
  sphi[0] = sin(phi);

  // Second photon antiparallel in fluid frame
  cth[1]  = -cth[0];
  sth[1]  =  sth[0];
  cphi[1] = -cphi[0];
  sphi[1] = -sphi[0];

  double E = nu*HPL/(ME*CL*CL);

  for (int n = 0; n < 2; n++) {
    // Set position
    tmp[n]->X[0] = t;
    coord(i, j, k, CENT, tmp[n]->X);

    // Get coordinate frame wavevector
    Kcov_tetrad[0] = -E;
    Kcov_tetrad[1] = E*cth[n];
    Kcov_tetrad[2] = E*cphi[n]*sth[n];
    Kcov_tetrad[3] = E*sphi[n]*sth[n];
    tetrad_to_coord(Ecov, Kcov_tetrad, tmp[n]->Kcov);

    // Set superphoton weight
    tmp[n]->w = 0.5*weight;

    // Diagnostics
    tmp[n]->origin[0] = nstep;
    tmp[n]->origin[1] = i;
    tmp[n]->origin[2] = j;
    tmp[n]->origin[3] = k;

    // Record radiation four-force
    for (int mu = 0; mu < NDIM; mu++) {
      #pragma omp atomic
      radG[i][j][k][mu] -= 1/(dt*dx[1]*dx[2]*dx[3])*kphys_to_num*tmp[n]->w*
                  tmp[n]->Kcov[mu];
    }
  }
}

#define TINY (1.e-200)
void get_dndlnu(int i, int j, int k, double dt, double dndlnu[NU_BINS+1], 
  double Ne, double Thetae, double Bmag)
{
  for (int n = 0; n < NU_BINS; n++) {
    dndlnu[n] = 0.;
  }

  double dndlnu_max = -1.e100;
  for (int n = 0; n <= NU_BINS; n++) {
    double Jsamp = Jnu(nusamp[n], Ne, Thetae, Bmag);
    Jsamp *= dx[1]*dx[2]*dx[3]*pow(L_unit,3.)*ggeom[i][j][CENT].g;

    dndlnu[n] = Jsamp/(wgtC/nusamp[n]*HPL + TINY);

    if (dndlnu[n] > dndlnu_max) {
      dndlnu_max = dndlnu[n];
    }
  }

  for (int n = 0; n <= NU_BINS; n++) {
    dndlnu[n] /= dndlnu_max;
  }
}
#undef TINY

void set_weight(grid_prim_type Prad)
{
  double Jtot;
  double zoneVol = dV*L_unit*L_unit*L_unit;

  // Set static variables
  Ns = tune_emiss;
  lnu_min = log(numin);
  lnu_max = log(numax);
  dlnu = (lnu_max - lnu_min)/NU_BINS;
  for (int n = 0; n <= NU_BINS; n++) {
    nusamp[n] = exp(n*dlnu + lnu_min);
  }
  Jtot = 0.;

  #pragma omp parallel
  {
    #pragma omp for collapse(3) reduction(+:Jtot)
    ZLOOP
    {
      double Ne, Thetae, Bmag;
      double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
      get_fluid_zone(i, j, k, Prad, &Ne, &Thetae, &Bmag, Ucon, Ucov, Bcon,
        Bcov);

      for (int n = 0; n <= NU_BINS; n++) {
        Jtot += Jnu(nusamp[n], Ne, Thetae, Bmag)*zoneVol*ggeom[i][j][CENT].g;
      }
    } // ZLOOP
  } // omp parallel

  wgtC = Jtot/(HPL*Ns)*nusamp[0];
  printf("wgtC = %e\n", wgtC);
}

// Use gsl's Gauss-Kronrod to integrate dNs/dlnu in each zone
void get_dnz(grid_prim_type Prad)
{
  #pragma omp parallel
  {
    gsl_integration_workspace *w = gsl_integration_workspace_alloc(1000);

    double f (double x, void *params) {
      double nu     = exp(x);
      double Ne     = ((double*)params)[0];
      double Thetae = ((double*)params)[1];
      double Bmag   = ((double*)params)[2];
      double Jsamp =  Jnu(nu, Ne, Thetae, Bmag)*nu;
      if (isinf(Jsamp)) {
        return 0.;
      } else {
        return Jsamp;
      }
    }

    double result, error;
    gsl_function F;
    F.function = &f;

    double zoneVol = dV*L_unit*L_unit*L_unit;
    double Ne, Thetae, Bmag;
    double Ucon[NDIM], Ucov[NDIM], Bcon[NDIM], Bcov[NDIM];
    
    #pragma omp for collapse(3)
    ZLOOP {
      // Ignore emission outside region of interest
      double X[NDIM];
      coord(i, j, k, CENT, X);
      if (X[1] < startx_rad[1] || X[2] > stopx_rad[1]) {
        dnzs[i][j][k] = 0.;
        continue;
      }

      get_fluid_zone(i, j, k, Prad, &Ne, &Thetae, &Bmag, Ucon, Ucov, Bcon, Bcov);

      // Get number of superphotons to be emitted
      double params[3] = {Ne, Thetae, Bmag};
      F.params = params;
      gsl_integration_qags(&F, lnu_min, lnu_max, 1.e100, 1.e-4, 1000, w,
        &result, &error);
      result /= HPL;
      result /= wgtC;
      result *= zoneVol;
      result *= ggeom[i][j][CENT].g;
      result *= dt*T_unit;
      if (isnan(result/(nthreads*mpi_nprocs()))) {
        dnzs[i][j][k] = 0.;
      } else {
        dnzs[i][j][k] = result/(nthreads*mpi_nprocs());
      }
    } // ZLOOP
    gsl_integration_workspace_free(w);
  } // pragma omp parallel
}

void *my_malloc(int cnt, int size)
{
  void *A = malloc(cnt*size);
  if (A == NULL) {
    fprintf(stderr, "Failed to malloc\n");
    exit(-1);
  }
  return A;
}
#endif // RADIATION

