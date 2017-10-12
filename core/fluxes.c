/******************************************************************************
 *                                                                            *
 * FLUXES.C                                                                   *
 *                                                                            *
 * CALCULATES FLUID FLUXES                                                    *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

static struct FluidState *Sl, *Sr;//, *Sl_shift;
static struct FluidEMF *emf;
int first_call = 1;

double get_flux(struct GridGeom *G, struct FluidState *S, struct FluidFlux *F)
{
  if (first_call) {
    Sl  = (struct FluidState*)malloc(sizeof(struct FluidState));
    //Sl_shift  = (struct FluidState*)malloc(sizeof(struct FluidState));
    Sr  = (struct FluidState*)malloc(sizeof(struct FluidState));
    emf = (struct FluidEMF*)malloc(sizeof(struct FluidEMF));
    first_call = 0;
  }

  //static GridPrim Pl, Pr;
  //static struct FluidState Sl, Sr;
  static GridDouble ctop;
  double cmax[NDIM], ndts[NDIM];
  memset(cmax, 0, NDIM*sizeof(double));
  memset(ndts, 0, NDIM*sizeof(double));

  // reconstruct X1
  reconstruct(S, Sl->P, Sr->P, 1);
  PLOOP ZLOOPALL if (isnan(Sl->P[ip][k][j][i])) {printf(" %i %i %i Sl->P\n", i, j, k); exit(-1);}
  PLOOP ZLOOPALL if (isnan(Sr->P[ip][k][j][i])) {printf(" %i %i %i Sr->P\n", i, j, k); exit(-1);}

  // lr_to_flux X1
  lr_to_flux(G, Sl, Sr, 1, FACE1, F->X1, ctop);

  PLOOP ZLOOPALL if (isnan(F->X1[ip][k][j][i])) {printf(" %i %i %i AOSIiSjjJD\n", i, j, k); exit(-1);}

  // cmax1
  timer_start(TIMER_CMAX);
  ZLOOP {
    cmax[1] = (ctop[k][j][i] > cmax[1] ? ctop[k][j][i] : cmax[1]);
  }
  timer_stop(TIMER_CMAX);

  // reconstruct X2
  reconstruct(S, Sl->P, Sr->P, 2);

  //printf("[U1] LR = %e %e\n", Sl->P[U1][3][4][6], Sr->P[U1][3][4][6]);

  // lr_to_flux X2
  lr_to_flux(G, Sl, Sr, 2, FACE2, F->X2, ctop);
  PLOOP ZLOOP if (isnan(F->X1[ip][k][j][i])) {printf("FX2\n"); exit(-1);}

  // cmax2
  timer_start(TIMER_CMAX);
  ZLOOP {
    cmax[2] = (ctop[k][j][i] > cmax[2] ? ctop[k][j][i] : cmax[2]);
  }
  timer_stop(TIMER_CMAX);

  // reconstruct X3
  reconstruct(S, Sl->P, Sr->P, 3);

  // lr_to_flux X3
  lr_to_flux(G, Sl, Sr, 3, FACE3, F->X3, ctop);
  PLOOP ZLOOP if (isnan(F->X1[ip][k][j][i])) {printf("FX3\n"); exit(-1);}

  // cmax3
  timer_start(TIMER_CMAX);
  ZLOOP {
    cmax[3] = (ctop[k][j][i] > cmax[3] ? ctop[k][j][i] : cmax[3]);
  }
  timer_stop(TIMER_CMAX);

  for (int mu = 1; mu < NDIM; mu++) {
    ndts[mu] = cour*dx[mu]/cmax[mu];
  }

  return 1./(1./ndts[1] + 1./ndts[2] + 1./ndts[3]);
}

// Note that the sense of L/R flips from zone to interface during function call
void lr_to_flux(struct GridGeom *G, struct FluidState *Sr,
  struct FluidState *Sl, int dir, int loc, GridPrim flux, GridDouble ctop)
{
  timer_start(TIMER_LR_TO_F);
  static GridPrim fluxL, fluxR;
  static GridDouble cmaxL, cmaxR, cminL, cminR, cmax, cmin;//, ctop;

  /*int max_indices[] = {0, N1-1, N2-1, N3-1};
  if (dir == 1) {
    max_indices[1] += 1;
  } else if (dir == 2) {
    max_indices[2] += 1;
  } else if (dir == 3) {
    max_indices[3] += 1;
  }*/
  int max_indices[] = {0, N1, N2, N3};

  //if (dir == 2)
  //printf("LRFLUX Pr[U1] = %e\n", Sr->P[U1][3][4][6]);
  //printf("DIR = %i: MAX_INDICES = %i %i %i %i\n", dir,
  //max_indices[0], max_indices[1], max_indices[2], max_indices[3]);

  // Properly offset left face
  PLOOP {
    if (dir == 1) {
      #pragma omp parallel for collapse(2)
      ZSLOOP_REVERSE(-1, max_indices[3], -1, max_indices[2], -1, max_indices[1])
        Sl->P[ip][k][j][i] = Sl->P[ip][k][j][i-1];
    } else if (dir == 2) {
      #pragma omp parallel for collapse(1)
      ZSLOOP_REVERSE(-1, max_indices[3], -1, max_indices[2], -1, max_indices[1])
        Sl->P[ip][k][j][i] = Sl->P[ip][k][j-1][i];
    } else if (dir == 3) {
      ZSLOOP_REVERSE(-1, max_indices[3], -1, max_indices[2], -1, max_indices[1])
        Sl->P[ip][k][j][i] = Sl->P[ip][k-1][j][i];
    }
  }
  //if (dir == 2)
  //printf("LRFLUX Pr[U1] = %e\n", Sr->P[U1][3][4][6]);

  // IF GDET < SMALL?!?!??!?! FLUX = 0!!!  MAXSPEED = 0!!!
  // MAYBE SING_FIX REMOVES THIS

  timer_start(TIMER_LR_STATE);

  get_state_vec(G, Sl, loc, -1, N3, -1, N2, -1, N1);

  get_state_vec(G, Sr, loc, -1, N3, -1, N2, -1, N1);

  //if (dir == 2)
  //printf("LRFLUX Pr[U1] = %e\n", Sr->P[U1][3][4][6]);
  /*#pragma omp parallel for collapse(2)
  KSLOOP(0, max_indices[3]) {
    JSLOOP(0, max_indices[2]) {
      #pragma omp simd
      ISLOOP(0, max_indices[1]) {
        get_state(G, Sl, i, j, k, loc);
      }
    }
  }
  #pragma omp parallel for collapse(2)
  KSLOOP(0, max_indices[3]) {
    JSLOOP(0, max_indices[2]) {
      #pragma omp simd
      ISLOOP(0, max_indices[1]) {
        get_state(G, Sr, i, j, k, loc);
      }
    }
  }*/
  timer_stop(TIMER_LR_STATE);

  timer_start(TIMER_LR_PTOF);

  prim_to_flux_vec(G, Sl, 0,   loc, -1, N3, -1, N2, -1, N1, Sl->U);
  prim_to_flux_vec(G, Sl, dir, loc, -1, N3, -1, N2, -1, N1, fluxL);

  prim_to_flux_vec(G, Sr, 0,   loc, -1, N3, -1, N2, -1, N1, Sr->U);
  prim_to_flux_vec(G, Sr, dir, loc, -1, N3, -1, N2, -1, N1, fluxR);

  timer_stop(TIMER_LR_PTOF);

  timer_start(TIMER_LR_VCHAR);
  #pragma omp parallel for simd collapse(2)
  ZSLOOP(-1, max_indices[3], -1, max_indices[2], -1, max_indices[1]) {
    mhd_vchar(G, Sl, i, j, k, loc, dir, cmaxL, cminL);
    mhd_vchar(G, Sr, i, j, k, loc, dir, cmaxR, cminR);
  }
  timer_stop(TIMER_LR_VCHAR);

  timer_start(TIMER_LR_CMAX);
  #pragma omp parallel for collapse(3)
  ZSLOOP(-1, max_indices[3], -1, max_indices[2], -1, max_indices[1]) {
    cmax[k][j][i] = fabs(MY_MAX(MY_MAX(0., cmaxL[k][j][i]), cmaxR[k][j][i]));
    cmin[k][j][i] = fabs(MY_MAX(MY_MAX(0., -cminL[k][j][i]), -cminR[k][j][i]));
    ctop[k][j][i] = MY_MAX(cmax[k][j][i], cmin[k][j][i]);
    if (isnan(1./ctop[k][j][i])) {printf("%i %i %i (%i)\n", k,j,i,dir); exit(-1);}
  }
  timer_stop(TIMER_LR_CMAX);

  timer_start(TIMER_LR_FLUX);

  PLOOP {
    //#pragma omp parallel for collapse(3)
    /*ZSLOOP(0, max_indices[3], 0, max_indices[2], 0, max_indices[1]) {
      flux[ip][k][j][i] = 0.5*(fluxL[ip][k][j][i] + fluxR[ip][k][j][i] -
        ctop[k][j][i]*(Sr->U[ip][k][j][i] - Sl->U[ip][k][j][i]));
    }*/
    #pragma omp parallel for simd collapse(2)
    ZSLOOP(-1, max_indices[3], -1, max_indices[2], -1, max_indices[1]) {
    //ZSLOOP(0, N3-1, 0, N2-1, 0, N1) {
      flux[ip][k][j][i] = 0.5*(fluxL[ip][k][j][i] + fluxR[ip][k][j][i] -
        ctop[k][j][i]*(Sr->U[ip][k][j][i] - Sl->U[ip][k][j][i]));
      if (isnan(flux[ip][k][j][i])) {
        printf("%e %e %e %e %e\n", fluxL[ip][k][j][i], fluxR[ip][k][j][i],
        ctop[k][j][i], Sr->U[ip][k][j][i], Sl->U[ip][k][j][i]);
      }
    }
  }
  timer_stop(TIMER_LR_FLUX);
  timer_stop(TIMER_LR_TO_F);
}

void flux_ct(struct FluidFlux *F)
{
  timer_start(TIMER_FLUX_CT);
  #pragma omp parallel
  {
    //printf("BEFORE F->X3[B1][3][64][64] = %e\n", F->X3[B1][3][64][64]);
    //printf("BEFORE F->X1[B3][3][64][64] = %e\n", F->X1[B3][3][64][64]);
    //printf("BEFORE F->X1[B3][2][64][64] = %e\n", F->X1[B3][2][64][64]);
    //printf("F->X1[B3][4][64][64] = %e\n", F->X1[B3][4][64][64]);
    // Calculate EMFs via average to corners (Toth approach)
    #pragma omp for collapse(3)
    ZSLOOP(0, N3, 0, N2, 0, N1) {
      emf->X3[k][j][i] =  0.25*(F->X1[B2][k][j][i] + F->X1[B2][k][j-1][i]
                              - F->X2[B1][k][j][i] - F->X2[B1][k][j][i-1]);
      emf->X2[k][j][i] = -0.25*(F->X1[B3][k][j][i] + F->X1[B3][k-1][j][i]
                              - F->X3[B1][k][j][i] - F->X3[B1][k][j][i-1]);
      if (isnan(emf->X2[k][j][i])) printf("%i %i %i %g %g %g %g\n", k,j,i,
        F->X1[B3][k][j][i], F->X1[B3][k-1][j][i], F->X3[B1][k][j][i], F->X3[B1][k][j][i-1]);
      emf->X1[k][j][i] =  0.25*(F->X2[B3][k][j][i] + F->X2[B3][k-1][j][i]
                              - F->X3[B2][k][j][i] - F->X3[B2][k][j-1][i]);
    }

    //printf("emf->X2[3][64][64] = %e\n", emf->X2[3][64][64]);
    //printf("emf->X2[4][64][64] = %e\n", emf->X2[4][65][64]);

    // Rewrite EMFs as fluxes, after Toth
    #pragma omp for collapse(3) nowait
    ZSLOOP(0, N3 - 1, 0, N2 - 1, 0, N1) {
      F->X1[B1][k][j][i] =  0.;
      F->X1[B2][k][j][i] =  0.5*(emf->X3[k][j][i] + emf->X3[k][j+1][i]);
      F->X1[B3][k][j][i] = -0.5*(emf->X2[k][j][i] + emf->X2[k+1][j][i]);
      //if (isnan(F->X1[B3][k][j][i])) printf("1\n");
    }
    #pragma omp for collapse(3) nowait
    ZSLOOP(0, N3 - 1, 0, N2, 0, N1 - 1) {
      F->X2[B1][k][j][i] = -0.5*(emf->X3[k][j][i] + emf->X3[k][j][i+1]);
      F->X2[B2][k][j][i] =  0.;
      F->X2[B3][k][j][i] =  0.5*(emf->X1[k][j][i] + emf->X1[k+1][j][i]);
    }
    #pragma omp for collapse(3)
    ZSLOOP(0, N3, 0, N2 - 1, 0, N1 - 1) {
      F->X3[B1][k][j][i] =  0.5*(emf->X2[k][j][i] + emf->X2[k][j][i+1]);
      F->X3[B2][k][j][i] = -0.5*(emf->X1[k][j][i] + emf->X1[k][j+1][i]);
      F->X3[B3][k][j][i] =  0.;
    }

    //printf("AFTER F->X1[B3][3][64][64] = %e\n", F->X1[B3][3][64][64]);
  } // omp parallel
  timer_stop(TIMER_FLUX_CT);
}

