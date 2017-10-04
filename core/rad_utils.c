/******************************************************************************
 *                                                                            *
 * RAD_UTILS.C                                                                *
 *                                                                            *
 * HELPER FUNCTIONS FOR RADIATION INFRASTRUCTURE                              *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

#if RADIATION
void init_rad(grid_prim_type Prad)
{
  set_units();
  
  photon_lists = malloc(nthreads*sizeof(struct of_photon*));
  #pragma omp parallel
  {
    photon_lists[omp_get_thread_num()] = NULL;
  }

  init_emissivity();

  set_weight(Prad);

  init_superphoton_resolution();
}

void init_superphoton_resolution()
{
}

double linear_interp_log(double x, double *table, double lx_min, double dlx)
{
  double lx = log(x);
  double dn = (lx - lx_min)/dlx;
  int n = (int)dn;
  dn = dn - n;

  return (1. - dn)*table[n] + dn*table[n+1];
}

void set_units()
{
  #if METRIC == MKS
  L_unit = GNEWT*mbh/(CL*CL);
  #endif
  T_unit = L_unit/CL;
  RHO_unit = M_unit*pow(L_unit,-3.);
  U_unit = RHO_unit*CL*CL;
  B_unit = CL*sqrt(4.*M_PI*RHO_unit);
  Ne_unit = RHO_unit/MP;
  #if ELECTRONS
  Thetae_unit = MP/ME;
  #else
  Thetae_unit = (gam-1.)*MP/ME/(1. + tp_over_te);
  #endif
  kphys_to_num = ME/M_unit;
}

// Remove superphoton from list and release memory
void list_remove(struct of_photon **ph, struct of_photon **ph_head,
  struct of_photon **ph_prev)
{
  if(*ph_prev != NULL)
  {
    (*ph_prev)->next = (*ph)->next;
    free(*ph);
    *ph = (*ph_prev)->next;
  } else {
    *ph_head = (*ph)->next;
    free(*ph);
    *ph = *ph_head;
  }
}

void get_fluid_zone(int i, int j, int k, grid_prim_type Prad, double *Ne,
  double *Thetae, double *B, double Ucon[NDIM], double Ucov[NDIM],
  double Bcon[NDIM], double Bcov[NDIM])
{
  double Bp[NDIM], Vcon[NDIM], Vfac, VdotV, UdotBp;

  *Ne = Prad[i][j][k][RHO]*Ne_unit;

  #if ELECTRONS
  *Thetae = Prad[i][j][k][KEL]*pow(Prad[i][j][k][RHO],game-1.)*Thetae_unit;
  #else
  *Thetae = Prad[i][j][k][UU]/Prad[i][j][k][RHO]*Thetae_unit;
  #endif

  *Thetae = MY_MIN(*Thetae, THETAE_MAX);

  Bp[1] = Prad[i][j][k][B1]*B_unit;
  Bp[2] = Prad[i][j][k][B2]*B_unit;
  Bp[3] = Prad[i][j][k][B3]*B_unit;

  Vcon[1] = Prad[i][j][k][U1];
  Vcon[2] = Prad[i][j][k][U2];
  Vcon[3] = Prad[i][j][k][U3];

  // Get Ucov
  VdotV = 0.;
  for(int l = 1; l < NDIM; l++) {
    for(int m = 1; m < NDIM; m++) {
      VdotV += ggeom[i][j][CENT].gcov[l][m]*Vcon[l]*Vcon[m];
    }
  }
  Vfac = sqrt(-1./ggeom[i][j][CENT].gcon[0][0]*(1. + fabs(VdotV)));
  Ucon[0] = -Vfac*ggeom[i][j][CENT].gcon[0][0];
  for(int l = 1; l < NDIM; l++) 
    Ucon[l] = Vcon[l] - Vfac*ggeom[i][j][CENT].gcon[0][l];
  lower(Ucon, &(ggeom[i][j][CENT]), Ucov);

  // Get Bcon, Bcov, and B
  UdotBp = 0.;
  for(int l = 1; l < NDIM; l++) 
    UdotBp += Ucov[l]*Bp[l];
  Bcon[0] = UdotBp;
  for(int l = 1; l < NDIM; l++) 
    Bcon[l] = (Bp[l] + Ucon[l]*UdotBp)/Ucon[0];
  lower(Bcon, &(ggeom[i][j][CENT]), Bcov);
  *B = sqrt(Bcon[0]*Bcov[0] + Bcon[1]*Bcov[1] + Bcon[2]*Bcov[2]
    + Bcon[3]*Bcov[3]);

  // Prevent highly magnetized regions from emitting
  double sigma = pow(*B/B_unit,2.)/(*Ne/Ne_unit);
  if (sigma > SIGMA_MAX) {
    *Thetae = SMALL;
  }
}
#endif // RADIATION

