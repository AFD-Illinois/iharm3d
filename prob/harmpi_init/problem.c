//Modified by Alexander Tchekhovskoy: MPI+3D
/***********************************************************************************
    Copyright 2006 Charles F. Gammie, Jonathan C. McKinney, Scott C. Noble, 
                   Gabor Toth, and Luca Del Zanna

                        HARM  version 1.0   (released May 1, 2006)

    This file is part of HARM.  HARM is a program that solves hyperbolic 
    partial differential equations in conservative form using high-resolution
    shock-capturing techniques.  This version of HARM has been configured to 
    solve the relativistic magnetohydrodynamic equations of motion on a 
    stationary black hole spacetime in Kerr-Schild coordinates to evolve
    an accretion disk model. 

    You are morally obligated to cite the following two papers in his/her 
    scientific literature that results from use of any part of HARM:

    [1] Gammie, C. F., McKinney, J. C., \& Toth, G.\ 2003, 
        Astrophysical Journal, 589, 444.

    [2] Noble, S. C., Gammie, C. F., McKinney, J. C., \& Del Zanna, L. \ 2006, 
        Astrophysical Journal, 641, 626.

   
    Further, we strongly encourage you to obtain the latest version of 
    HARM directly from our distribution website:
    http://rainman.astro.uiuc.edu/codelib/


    HARM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    HARM is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with HARM; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

***********************************************************************************/

/*
 *
 * generates initial conditions for a fishbone & moncrief disk 
 * with exterior at minimum values for density & internal energy.
 *
 * cfg 8-10-01
 *
 */

#include "decs.h"
#include <float.h>

struct of_geom {
	double gcon[NDIM][NDIM];
	double gcov[NDIM][NDIM];
	double g;
	double alpha;
};

void coord_transform(struct GridGeom *G, struct FluidState *S, int i, int j,int k) ;
double compute_Amax( double (*A)[N2+2*NG][N3+2*NG] );
double compute_B_from_A( struct GridGeom *G, struct FluidState *S, double (*A)[N2+2*NG][N3+2*NG]);
double normalize_B_by_maxima_ratio(struct GridGeom *G, struct FluidState *S, double beta_target, double *norm_value);
double normalize_B_by_beta(struct GridGeom *G, struct FluidState *S, double beta_target, double rmax, double *norm_value);
void blgset(int i, int j, struct of_geom *geom);
double bl_gdet_func(double r, double th);
void bl_gcov_func(double r, double th, double gcov[][NDIM]);
void bl_gcon_func(double r, double th, double gcon[][NDIM]);

/////////////////////
//magnetic field geometry and normalization
#define NORMALFIELD (0)
#define MADFIELD (1)

#define WHICHFIELD MADFIELD

#define NORMALIZE_FIELD_BY_MAX_RATIO (1)
#define NORMALIZE_FIELD_BY_BETAMIN (2)
#define WHICH_FIELD_NORMALIZATION NORMALIZE_FIELD_BY_MAX_RATIO
//end magnetic field
//////////////////////

//////////////////////
//torus density normalization
#define THINTORUS_NORMALIZE_DENSITY (1)
#define DOAUTOCOMPUTEENK0 (1)

#define NORMALIZE_BY_TORUS_MASS (1)
#define NORMALIZE_BY_DENSITY_MAX (2)

#define DENSITY_NORMALIZATION NORMALIZE_BY_DENSITY_MAX
//torus density normalization
//////////////////////

double rmax = 0.;
double rhomax = 1.;

void set_problem_params() {}

void init(struct GridGeom *G, struct FluidState *S)
{
  int i,j,k ;
  double r,th,sth,cth ;
  double ur,uh,up,u,rho ;
  double X[NDIM] ;

  /* for disk interior */
  double l,rin,lnh,expm2chi,up1 ;
  double DD,AA,SS,thin,sthin,cthin,DDin,AAin,SSin ;
  double kappa,hm1 ;

  /* for magnetic field */
  double A[N1+2*NG][N2+2*NG][N3+2*NG] ;
  double rho_av,umax,beta,bsq_max,norm,q,beta_act ;
  double lfish_calc(double rmax) ;
  
  int iglob, jglob, kglob;
  double rancval;
  
  double aphipow;
  
  /* some physics parameters */
  gam = 5./3. ;

  /* disk parameters (use fishbone.m to select new solutions) */
  a = 0.9 ;
  rin = 6. ;
  rmax = 13. ;
  l = lfish_calc(rmax) ;

  kappa =1.e-3;
  beta = 100. ;

  /* some numerical parameters */
  lim = MC ;
  failed = 0 ;	/* start slow */
  cour = .8 ;
  dt = 1.e-5 ;
  R0 = 0.0 ;
  Rin = 0.87*(1. + sqrt(1. - a*a)) ;  //.98
  Rout = 1e5;
  
  //Extra grid params VHARM uses
  Rhor = (1. + sqrt(1. - a*a));
  double z1 = 1 + pow(1 - a*a,1./3.)*(pow(1+a,1./3.) + pow(1-a,1./3.));
  double z2 = sqrt(3*a*a + z1*z1);
  Risco = 3 + z2 - sqrt((3-z1)*(3 + z1 + 2*z2));
  
//  rbr = 400.;
//  npow2=4.0; //power exponent
//  cpow2=1.0; //exponent prefactor (the larger it is, the more hyperexponentiation is)


  t = 0. ;
  hslope = 0.3 ;

//  if(N2!=1) {
//    //2D problem, use full pi-wedge in theta
//    fractheta = 1.;
//  }
//  else{
//    //1D problem (since only 1 cell in theta-direction), use a restricted theta-wedge
//    fractheta = 1.e-2;
//  }
//
//  fracphi = 0.5;

//  set_arrays() ;
  zero_arrays();
  set_grid(G) ;

  // TODO isn't this just bl_coord?
  // Reads from phys_coords which is set by bl_coord_vec
  //get_phys_coord(5,0,0,&r,&th) ;
  coord(NG+5,NG,NG,CENT,X) ;
  bl_coord(X,&r,&th);

  if(mpi_io_proc()) {
    fprintf(stderr,"r[5]: %g\n",r) ;
    fprintf(stderr,"r[5]/rhor: %g",r/(1. + sqrt(1. - a*a))) ;
    if( r > 1. + sqrt(1. - a*a) ) {
      fprintf(stderr, ": INSUFFICIENT RESOLUTION, ADD MORE CELLS INSIDE THE HORIZON\n" );
    }
    else {
      fprintf(stderr, "\n");
    }
  }

  /* output choices */
  tf = 10000.0 ;

  DTd = 10.; /* dumping frequency, in units of M */
  DTl = 10. ;	/* logfile frequency, in units of M */
  DTr = 10. ; /* restart file frequ., in units of M */

  /* start diagnostic counters */
  dump_cnt = 0 ;
  rdump_cnt = 0 ;

  rhomax = 0. ;
  umax = 0. ;
  //ZSLOOP(0,N1-1,0,N2-1,0,N3-1) {
  for(iglob=0;iglob<N1TOT;iglob++) {
    for(jglob=0;jglob<N2TOT;jglob++) {
      for(kglob=0;kglob<N3TOT;kglob++) {
        
        rancval = get_random();
        i = iglob - global_start[0];
        j = jglob - global_start[1];
        k = kglob - global_start[2];
        if(i<0 ||
           j<0 ||
           k<0 ||
           i>=N1 ||
           j>=N2 ||
           k>=N3){
          continue;
        }

        //get_phys_coord(i,j,k,&r,&th) ;
        coord(i,j,k,CENT,X) ;
        bl_coord(X,&r,&th);

        sth = sin(th) ;
        cth = cos(th) ;

        /* calculate lnh */
        DD = r*r - 2.*r + a*a ;
        AA = (r*r + a*a)*(r*r + a*a) - DD*a*a*sth*sth ;
        SS = r*r + a*a*cth*cth ;

        thin = M_PI/2. ;
        sthin = sin(thin) ;
        cthin = cos(thin) ;
        DDin = rin*rin - 2.*rin + a*a ;
        AAin = (rin*rin + a*a)*(rin*rin + a*a) 
                - DDin*a*a*sthin*sthin ;
        SSin = rin*rin + a*a*cthin*cthin ;

        if(r >= rin) {
          lnh = 0.5*log((1. + sqrt(1. + 4.*(l*l*SS*SS)*DD/
                  (AA*sth*AA*sth)))/(SS*DD/AA)) 
                  - 0.5*sqrt(1. + 4.*(l*l*SS*SS)*DD/(AA*AA*sth*sth))
                  - 2.*a*r*l/AA 
                  - (0.5*log((1. + sqrt(1. + 4.*(l*l*SSin*SSin)*DDin/
                  (AAin*AAin*sthin*sthin)))/(SSin*DDin/AAin)) 
                  - 0.5*sqrt(1. + 4.*(l*l*SSin*SSin)*DDin/
                          (AAin*AAin*sthin*sthin)) 
                  - 2.*a*rin*l/AAin ) ;
        }
        else
          lnh = 1. ;


        /* regions outside torus */
        if(lnh < 0. || r < rin) {
          //reset density and internal energy to zero outside torus
          rho = 0.; //1.e-7*RHOMIN ;
          u = 0.; //1.e-7*UUMIN ;

          ur = 0. ;
          uh = 0. ;
          up = 0. ;

          S->P[RHO][k][j][i] = rho ;
          S->P[UU][k][j][i] = u ;
          S->P[U1][k][j][i] = ur ;
          S->P[U2][k][j][i] = uh ;
          S->P[U3][k][j][i] = up ;
        }
        /* region inside magnetized torus; u^i is calculated in
         * Boyer-Lindquist coordinates, as per Fishbone & Moncrief,
         * so it needs to be transformed at the end */
        else { 
          hm1 = exp(lnh) - 1. ;
          rho = pow(hm1*(gam - 1.)/(kappa*gam),
                                  1./(gam - 1.)) ; 
          u = kappa*pow(rho,gam)/(gam - 1.) ;
          ur = 0. ;
          uh = 0. ;

          /* calculate u^phi */
          expm2chi = SS*SS*DD/(AA*AA*sth*sth) ;
          up1 = sqrt((-1. + sqrt(1. + 4.*l*l*expm2chi))/2.) ;
          up = 2.*a*r*sqrt(1. + up1*up1)/sqrt(AA*SS*DD) +
                  sqrt(SS/AA)*up1/sth ;


          S->P[RHO][k][j][i] = rho ;
          if(rho > rhomax) rhomax = rho ;
          S->P[UU][k][j][i] = u*(1. + 4.e-2*(rancval-0.5)) ;
          if(u > umax && r > rin) umax = u ;
          S->P[U1][k][j][i] = ur ;
          S->P[U2][k][j][i] = uh ;

          S->P[U3][k][j][i] = up ;

          /* convert from 4-vel in BL coords to relative 4-vel in code coords */
          coord_transform(G,S,i,j,k) ;
        }

        S->P[B1][k][j][i] = 0. ;
        S->P[B2][k][j][i] = 0. ;
        S->P[B3][k][j][i] = 0. ;
      }
    }
  }

#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&rhomax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&umax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  
  /* Normalize the densities so that max(rho) = 1 */
  if(mpi_io_proc()) fprintf(stderr,"rhomax: %g\n",rhomax) ;
  ZSLOOP(0,N3-1,0,N2-1,0,N1-1) {
          S->P[RHO][k][j][i] /= rhomax ;
          S->P[UU][k][j][i]  /= rhomax ;
  }
  umax /= rhomax ;
  kappa *= pow(rhomax,gam-1);
//  global_kappa = kappa;
  rhomax = 1. ;

  //(fixup here)?
  set_bounds(G, S);

  if (WHICHFIELD == NORMALFIELD) {
    aphipow = 0.;
  } else if (WHICHFIELD == MADFIELD) {
    aphipow = 5.;
  } else {
    fprintf(stderr, "Unknown field type: %d\n", (int)WHICHFIELD);
    exit(321);
  }

  /* first find corner-centered vector potential */
  ZSLOOP(0,N3-1+NG,0,N2-1+NG,0,N1-1+NG) A[i][j][k] = 0. ;
  ZSLOOP(0,N3-1+NG,0,N2-1+NG,0,N1-1+NG) {
          /* radial field version */
          /*
          coord(i,j,k,CORN,X) ;
          bl_coord(X,&r,&th) ;
         
          A[i][j][k] = (1-cos(th)) ;
          */

    
          /* vertical field version */
          /*
          coord(i,j,k,CORN,X) ;
          bl_coord(X,&r,&th) ;

          A[i][j][k] = r*r*sin(th)*sin(th) ;
          */
    
    

          /* field-in-disk version */
          /* flux_ct */
    
      //cannot use get_phys_coords() here because it can only provide coords at CENT
      coord(i,j,k,CORN,X) ;
      bl_coord(X,&r,&th) ;


          rho_av = 0.25*(
                  S->P[RHO][k][j][i] +
                  S->P[RHO][k][j][i-1] +
                  S->P[RHO][k][j-1][i] +
                  S->P[RHO][k][j-1][i-1]) ;

          q = pow(r,aphipow)*rho_av/rhomax ;
          if (WHICHFIELD == NORMALFIELD) {
            q -= 0.2;
          }
          if(q > 0.) A[i][j][k] = q ;

  }
  
  fixup(G, S);

  /* now differentiate to find cell-centered B,
     and begin normalization */
  
  bsq_max = compute_B_from_A(G, S, A);
  
    if(mpi_io_proc())
      fprintf(stderr,"initial bsq_max: %g\n",bsq_max) ;

    /* finally, normalize to set field strength */
    beta_act =(gam - 1.)*umax/(0.5*bsq_max) ;

    if(mpi_io_proc())
      fprintf(stderr,"initial beta: %g (should be %g)\n",beta_act,beta) ;
    
    if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_BETAMIN) {
      beta_act = normalize_B_by_beta(G, S, beta, 10*rmax, &norm);
    } else if(WHICH_FIELD_NORMALIZATION == NORMALIZE_FIELD_BY_MAX_RATIO) {
      beta_act = normalize_B_by_maxima_ratio(G, S, beta, &norm);
    } else {
      if(mpi_io_proc()) {
        fprintf(stderr, "Unknown magnetic field normalization %d\n",
                WHICH_FIELD_NORMALIZATION);
        MPI_Finalize();
        exit(2345);
      }

    if(mpi_io_proc())
      fprintf(stderr,"final beta: %g (should be %g)\n",beta_act,beta) ;
    }

    
  /* enforce boundary conditions */
  fixup(G, S);
  set_bounds(G, S); //was bound_prim


#if( DO_FONT_FIX )
  set_Katm();
#endif 


}

//note that only axisymmetric A is supported
double compute_Amax( double (*A)[N2+2*NG][N3+2*NG] )
{
  double Amax = 0.;
  
  ZLOOP {
    if(A[i][j][k] > Amax) Amax = A[i][j][k];
  }
  
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&Amax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  
  return(Amax);
}


//note that only axisymmetric A is supported
double compute_B_from_A(struct GridGeom *G, struct FluidState *S, double (*A)[N2+2*NG][N3+2*NG])
{
  double bsq_max = 0., bsq_ij ;
  
  ZLOOP {
    
    /* flux-ct */
    S->P[B1][k][j][i] = -(A[i][j][k] - A[i][j+1][k]
                       + A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[2]*G->gdet[CENT][j][i]) ;
    S->P[B2][k][j][i] = (A[i][j][k] + A[i][j+1][k]
                      - A[i+1][j][k] - A[i+1][j+1][k])/(2.*dx[1]*G->gdet[CENT][j][i]) ;
    
    S->P[B3][k][j][i] = 0. ;
    
    get_state(G, S, i, j, k, CENT);
    bsq_ij = bsq_calc(S, i,j,k) ;
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

  return(bsq_max);
}

double normalize_B_by_maxima_ratio(struct GridGeom *G, struct FluidState *S, double beta_target, double *norm_value)
{
  double beta_act, bsq_ij, u_ij, umax = 0., bsq_max = 0.;
  double norm;
  
  ZLOOP {
    get_state(G, S, i, j, k, CENT);
    bsq_ij = bsq_calc(S, i, j, k);
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
    u_ij = S->P[UU][k][j][i];
    if(u_ij > umax) umax = u_ij;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&umax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif

  /* finally, normalize to set field strength */
  beta_act =(gam - 1.)*umax/(0.5*bsq_max) ;
  
  norm = sqrt(beta_act/beta_target) ;
  bsq_max = 0. ;
  ZLOOP {
    S->P[B1][k][j][i] *= norm ;
    S->P[B2][k][j][i] *= norm ;
    
    get_state(G, S, i, j, k, CENT);
    bsq_ij = bsq_calc(S, i, j, k);
    if(bsq_ij > bsq_max) bsq_max = bsq_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&bsq_max,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#endif
  
  beta_act = (gam - 1.)*umax/(0.5*bsq_max) ;

  if(norm_value) {
    *norm_value = norm;
  }
  return(beta_act);
}

//normalize the magnetic field using the values inside r < rmax
double normalize_B_by_beta(struct GridGeom *G, struct FluidState *S, double beta_target, double rmax, double *norm_value)
{
  double beta_min = 1e100, beta_ij, beta_act, bsq_ij, u_ij;
  double norm;
  double X[NDIM], r, th;
  
  ZLOOP {
    coord(i, j, k, CENT, X);
    bl_coord(X, &r, &th);
    if (r>rmax) {
      continue;
    }
    get_state(G, S, i, j, k, CENT);
    bsq_ij = bsq_calc(S, i, j, k);

    u_ij = S->P[UU][k][j][i];
    beta_ij = (gam - 1.)*u_ij/(0.5*(bsq_ij+SMALL)) ;
    if(beta_ij < beta_min) beta_min = beta_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&beta_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
  
  /* finally, normalize to set field strength */
  beta_act = beta_min;
  
  norm = sqrt(beta_act/beta_target) ;
  beta_min = 1e100;
  ZLOOP {
    S->P[B1][k][j][i] *= norm ;
    S->P[B2][k][j][i] *= norm ;
    S->P[B3][k][j][i] *= norm ;
    get_state(G, S, i, j, k, CENT);
    bsq_ij = bsq_calc(S, i, j, k);
    u_ij = S->P[UU][k][j][i];
    beta_ij = (gam - 1.)*u_ij/(0.5*(bsq_ij+SMALL)) ;
    if(beta_ij < beta_min) beta_min = beta_ij ;
  }
#ifdef MPI
  //exchange the info between the MPI processes to get the true max
  MPI_Allreduce(MPI_IN_PLACE,&beta_min,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
#endif
  
  beta_act = beta_min;

  if(norm_value) {
    *norm_value = norm;
  }

  return(beta_act);
}

// Convert Boyer-Lindquist four-velocity to MKS 3-velocity
void coord_transform(struct GridGeom *G, struct FluidState *S, int i, int j,
  int k)
{
  double X[NDIM], r, th, ucon[NDIM], trans[NDIM][NDIM], tmp[NDIM];
  double AA, BB, CC, discr;
  double alpha, gamma, beta[NDIM];
  struct blgeom;
  struct of_geom blgeom;

  coord(i, j, k, CENT, X);
  bl_coord(X, &r, &th);
  blgset(i, j, &blgeom);

  ucon[1] = S->P[U1][k][j][i];
  ucon[2] = S->P[U2][k][j][i];
  ucon[3] = S->P[U3][k][j][i];

  AA = blgeom.gcov[0][0];
  BB = 2.*(blgeom.gcov[0][1]*ucon[1] +
           blgeom.gcov[0][2]*ucon[2] +
           blgeom.gcov[0][3]*ucon[3]);
  CC = 1. +
      blgeom.gcov[1][1]*ucon[1]*ucon[1] +
      blgeom.gcov[2][2]*ucon[2]*ucon[2] +
      blgeom.gcov[3][3]*ucon[3]*ucon[3] +
      2.*(blgeom.gcov[1][2]*ucon[1]*ucon[2] +
          blgeom.gcov[1][3]*ucon[1]*ucon[3] +
          blgeom.gcov[2][3]*ucon[2]*ucon[3]);

  discr = BB*BB - 4.*AA*CC;
  ucon[0] = (-BB - sqrt(discr))/(2.*AA);
  // This is ucon in BL coords

  // transform to Kerr-Schild
  // Make transform matrix
  memset(trans, 0, 16*sizeof(double));
  for (int mu = 0; mu < NDIM; mu++) {
    trans[mu][mu] = 1.;
  }
  trans[0][1] = 2.*r/(r*r - 2.*r + a*a);
  trans[3][1] = a/(r*r - 2.*r + a*a);

  // Transform ucon
  for (int mu = 0; mu < NDIM; mu++) {
    tmp[mu] = 0.;
  }
  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      tmp[mu] += trans[mu][nu]*ucon[nu];
    }
  }
  for (int mu = 0; mu < NDIM; mu++) {
    ucon[mu] = tmp[mu];
  }

  // This is ucon in KS coords

  // Transform to KS' coords
  //ucon[1] *= (1. / (r - R0));
  ucon[1] /= dr_dx(X[1]);
  //ucon[2] *=
  //    (1. / (M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2])));
  ucon[2] /= dth_dx(X[2]);

  // Solve for v. Use same u^t, unchanged under KS -> KS'
  //geom = get_geometry(ii, jj, CENT);
  alpha = G->lapse[CENT][j][i];
  //alpha = 1.0 / sqrt( -geom->gcon[0][0] ) ;
  gamma = ucon[0]*alpha;

  beta[1] = alpha*alpha*G->gcon[CENT][0][1][j][i];//geom->gcon[0][1];
  beta[2] = alpha*alpha*G->gcon[CENT][0][2][j][i];//geom->gcon[0][2];
  beta[3] = alpha*alpha*G->gcon[CENT][0][3][j][i];//geom->gcon[0][3];

  S->P[U1][k][j][i] = ucon[1] + beta[1]*gamma/alpha;
  S->P[U2][k][j][i] = ucon[2] + beta[2]*gamma/alpha;
  S->P[U3][k][j][i] = ucon[3] + beta[3]*gamma/alpha;
}


double lfish_calc(double r)
{
	return(
   ((pow(a,2) - 2.*a*sqrt(r) + pow(r,2))*
      ((-2.*a*r*(pow(a,2) - 2.*a*sqrt(r) + pow(r,2)))/
         sqrt(2.*a*sqrt(r) + (-3. + r)*r) +
        ((a + (-2. + r)*sqrt(r))*(pow(r,3) + pow(a,2)*(2. + r)))/
         sqrt(1 + (2.*a)/pow(r,1.5) - 3./r)))/
    (pow(r,3)*sqrt(2.*a*sqrt(r) + (-3. + r)*r)*(pow(a,2) + (-2. + r)*r))
	) ;
}

// SUPPORT FUNCTIONS FROM HARM/HARMPI


void blgset(int i, int j, struct of_geom *geom)
{
  double r, th, X[NDIM];

  coord(i, j, 0, CENT, X);
  bl_coord(X, &r, &th);

  if (th < 0)
    th *= -1.;
  if (th > M_PI)
    th = 2. * M_PI - th;

  geom->g = bl_gdet_func(r, th);
  bl_gcov_func(r, th, geom->gcov);
  bl_gcon_func(r, th, geom->gcon);
}

double bl_gdet_func(double r, double th)
{
  double a2, r2;

  a2 = a * a;
  r2 = r * r;
  return (r * r * fabs(sin(th)) *
    (1. + 0.5 * (a2 / r2) * (1. + cos(2. * th))));
}

void bl_gcov_func(double r, double th, double gcov[NDIM][NDIM])
{
  double sth, cth, s2, a2, r2, DD, mu;

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      gcov[mu][nu] = 0.;
    }
  }

  sth = fabs(sin(th));
  s2 = sth*sth;
  cth = cos(th);
  a2 = a*a;
  r2 = r*r;
  DD = 1. - 2./r + a2/r2;
  mu = 1. + a2*cth*cth/r2;

  gcov[0][0] = -(1. - 2./(r*mu));
  gcov[0][3]  = -2.*a*s2/(r*mu);
  gcov[3][0]  = gcov[0][3];
  gcov[1][1]   = mu/DD;
  gcov[2][2]   = r2*mu;
  gcov[3][3]   = r2*sth*sth*(1. + a2/r2 + 2.*a2*s2/(r2*r*mu));

}

void bl_gcon_func(double r, double th, double gcon[NDIM][NDIM])
{
  double sth, cth, a2, r2, r3, DD, mu;

  for (int mu = 0; mu < NDIM; mu++) {
    for (int nu = 0; nu < NDIM; nu++) {
      gcon[mu][nu] = 0.;
    }
  }

  sth = sin(th);
  cth = cos(th);

#if(COORDSINGFIX)
  if (fabs(sth) < SINGSMALL) {
    if (sth >= 0)
      sth = SINGSMALL;
    if (sth < 0)
      sth = -SINGSMALL;
  }
#endif

  a2 = a*a;
  r2 = r*r;
  r3 = r2*r;
  DD = 1. - 2./r + a2/r2;
  mu = 1. + a2*cth*cth/r2;

  gcon[0][0] = -1. - 2.*(1. + a2/r2)/(r*DD*mu);
  gcon[0][3] = -2.*a/(r3*DD*mu);
  gcon[3][0] = gcon[0][3];
  gcon[1][1] = DD/mu;
  gcon[2][2] = 1./(r2*mu);
  gcon[3][3] = (1. - 2./(r*mu))/(r2*sth*sth*DD);
}
