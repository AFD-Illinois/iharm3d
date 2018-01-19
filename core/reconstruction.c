/******************************************************************************
 *                                                                            *
 * RECONSTRUCTION.C                                                           *
 *                                                                            *
 * RECONSTRUCTION ALGORITHMS                                                  *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

// Performs the slope-limiting for the numerical flux calculation
double slope_lim(double y1, double y2, double y3)
{
  double Dqm, Dqp, Dqc, s;

  // Woodward, or monotonized central, slope limiter
  if (lim == MC) {
    Dqm = 2. * (y2 - y1);
    Dqp = 2. * (y3 - y2);
    Dqc = 0.5 * (y3 - y1);
    s = Dqm * Dqp;
    if (s <= 0.)
      return 0.;
    else {
      if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
        return (Dqm);
      else if (fabs(Dqp) < fabs(Dqc))
        return (Dqp);
      else
        return (Dqc);
    }
  }

  // van Leer slope limiter
  else if (lim == VANL) {
    Dqm = (y2 - y1);
    Dqp = (y3 - y2);
    s = Dqm * Dqp;
    if (s <= 0.)
      return 0.;
    else
      return (2. * s / (Dqm + Dqp));
  }

  // Minmod slope limiter (crude but robust)
  else if (lim == MINM) {
    Dqm = (y2 - y1);
    Dqp = (y3 - y2);
    s = Dqm * Dqp;
    if (s <= 0.)
      return 0.;
    else if (fabs(Dqm) < fabs(Dqp))
      return Dqm;
    else
      return Dqp;
  }

  fprintf(stderr, "unknown slope limiter\n");
  exit(10);

  return (0.);
}

void linear_mc(double x1, double x2, double x3, double *lout, double *rout)
{
  double Dqm,Dqp,Dqc,s;

  Dqm = 2. * (x2 - x1);
  Dqp = 2. * (x3 - x2);
  Dqc = 0.5 * (x3 - x1);

  s = Dqm * Dqp;

  if (s <= 0.)
    s = 0.;
  else {
    if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
      s = Dqm;
    else if (fabs(Dqp) < fabs(Dqc))
      s = Dqp;
    else
      s = Dqc;
  }

  // Reconstruct left, right
  *lout = x2 - 0.5*s;
  *rout = x2 + 0.5*s;
}

// Parabolic interpolation (see Colella & Woodward 1984; CW)
// Implemented by Xiaoyue Guan
void para(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{
  double y[5], dq[5];
  double Dqm, Dqc, Dqp, aDqm,aDqp,aDqc,s,l,r,qa, qd, qe;

  y[0] = x1;
  y[1] = x2;
  y[2] = x3;
  y[3] = x4;
  y[4] = x5;

  // CW 1.7
  for(int i = 1; i < 4; i++) {
    Dqm = 2. *(y[i  ] - y[i-1]);
    Dqp = 2. *(y[i+1] - y[i  ]);
    Dqc = 0.5*(y[i+1] - y[i-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

    // CW 1.8
    if (s <= 0.) {
      dq[i] = 0.;
    } else {
      dq[i] = MY_MIN(aDqc,MY_MIN(aDqm,aDqp))*MY_SIGN(Dqc);
    }
  }

  // CW 1.6
  l = 0.5*(y[2] + y[1]) - (dq[2] - dq[1])/6.0;
  r = 0.5*(y[3] + y[2]) - (dq[3] - dq[2])/6.0;

  qa = (r - y[2])*(y[2] - l);
  qd = (r - l);
  qe = 6.0*(y[2] - 0.5*(l + r));

  if (qa <= 0.) {
    l = y[2];
    r = y[2];
  }

  if (qd*(qd - qe) < 0.0) {
    l = 3.0*y[2] - 2.0*r;
  } else if (qd*(qd + qe) < 0.0) {
    r = 3.0*y[2] - 2.0*l;
  }

  *lout = l;
  *rout = r;
}

// WENO interpolation. See Tchekhovskoy et al. 2007 (T07), Shu 2011 (S11)
// Implemented by Monika Moscibrodzka
void weno(double x1, double x2, double x3, double x4, double x5, double *lout, double *rout)
{
  // S11 1, 2, 3
  double vr[3], vl[3];
  vr[0] =  (3./8.)*x1 - (5./4.)*x2 + (15./8.)*x3;
  vr[1] = (-1./8.)*x2 + (3./4.)*x3 + (3./8.)*x4;
  vr[2] =  (3./8.)*x3 + (3./4.)*x4 - (1./8.)*x5;

  vl[0] =  (3./8.)*x5 - (5./4.)*x4 + (15./8.)*x3;
  vl[1] = (-1./8.)*x4 + (3./4.)*x3 + (3./8.)*x2;
  vl[2] =  (3./8.)*x3 + (3./4.)*x2 - (1./8.)*x1;

  // Smoothness indicators, T07 A18 or S11 8
  double beta[3];
  beta[0] = (13./12.)*pow(x1 - 2.*x2 + x3, 2) +
            (1./4.)*pow(x1 - 4.*x2 + 3.*x3, 2);
  beta[1] = (13./12.)*pow(x2 - 2.*x3 + x4, 2) +
            (1./4.)*pow(x4 - x2, 2);
  beta[2] = (13./12.)*pow(x3 - 2.*x4 + x5, 2) +
            (1./4.)*pow(x5 - 4.*x4 + 3.*x3, 2);

  // Nonlinear weights S11 9
  double den, wtr[3], Wr, wr[3], wtl[3], Wl, wl[3], eps;
  eps=1.e-26;

  den = eps + beta[0]; den *= den; wtr[0] = (1./16.)/den;
  den = eps + beta[1]; den *= den; wtr[1] = (5./8. )/den;
  den = eps + beta[2]; den *= den; wtr[2] = (5./16.)/den;
  Wr = wtr[0] + wtr[1] + wtr[2];
  wr[0] = wtr[0]/Wr ;
  wr[1] = wtr[1]/Wr ;
  wr[2] = wtr[2]/Wr ;

  den = eps + beta[2]; den *= den; wtl[0] = (1./16.)/den;
  den = eps + beta[1]; den *= den; wtl[1] = (5./8. )/den;
  den = eps + beta[0]; den *= den; wtl[2] = (5./16.)/den;
  Wl = wtl[0] + wtl[1] + wtl[2];
  wl[0] = wtl[0]/Wl;
  wl[1] = wtl[1]/Wl;
  wl[2] = wtl[2]/Wl;

  *lout = vl[0]*wl[0] + vl[1]*wl[1] + vl[2]*wl[2];
  *rout = vr[0]*wr[0] + vr[1]*wr[1] + vr[2]*wr[2];
}

// MP5 reconstruction from PLUTO
// Imported by Mani Chandra
#define MINMOD(a, b) ((a)*(b) > 0.0 ? (fabs(a) < fabs(b) ? (a):(b)):0.0)
double Median(double a, double b, double c)
{
  return (a + MINMOD(b - a, c - a));
}
double mp5_subcalc(double Fjm2, double Fjm1, double Fj, double Fjp1, double Fjp2)
{
  double f, d2, d2p, d2m;
  double dMMm, dMMp;
  double scrh1,scrh2, Fmin, Fmax;
  double fAV, fMD, fLC, fUL, fMP;
  static double alpha = 4.0, epsm = 1.e-12;

  f  = 2.0*Fjm2 - 13.0*Fjm1 + 47.0*Fj + 27.0*Fjp1 - 3.0*Fjp2;
  f /= 60.0;

  fMP = Fj + MINMOD(Fjp1 - Fj, alpha*(Fj - Fjm1));

  if ((f - Fj)*(f - fMP) <= epsm)
    return f;

  d2m = Fjm2 + Fj   - 2.0*Fjm1;              // Eqn. 2.19
  d2  = Fjm1 + Fjp1 - 2.0*Fj;
  d2p = Fj   + Fjp2 - 2.0*Fjp1;              // Eqn. 2.19

  scrh1 = MINMOD(4.0*d2 - d2p, 4.0*d2p - d2);
  scrh2 = MINMOD(d2, d2p);
  dMMp  = MINMOD(scrh1,scrh2);               // Eqn. 2.27

  scrh1 = MINMOD(4.0*d2m - d2, 4.0*d2 - d2m);
  scrh2 = MINMOD(d2, d2m);
  dMMm  = MINMOD(scrh1,scrh2);               // Eqn. 2.27

  fUL = Fj + alpha*(Fj - Fjm1);              // Eqn. 2.8
  fAV = 0.5*(Fj + Fjp1);                     // Eqn. 2.16
  fMD = fAV - 0.5*dMMp;                      // Eqn. 2.28
  fLC = 0.5*(3.0*Fj - Fjm1) + 4.0/3.0*dMMm;  // Eqn. 2.29

  scrh1 = fmin(Fj, Fjp1); scrh1 = fmin(scrh1, fMD);
  scrh2 = fmin(Fj, fUL);    scrh2 = fmin(scrh2, fLC);
  Fmin  = fmax(scrh1, scrh2);                // Eqn. (2.24a)

  scrh1 = fmax(Fj, Fjp1); scrh1 = fmax(scrh1, fMD);
  scrh2 = fmax(Fj, fUL);    scrh2 = fmax(scrh2, fLC);
  Fmax  = fmin(scrh1, scrh2);                // Eqn. 2.24b

  f = Median(f, Fmin, Fmax);                 // Eqn. 2.26
  return f;
}

void mp5(double x1, double x2, double x3, double x4, double x5, double *lout,
  double *rout)
{
  *rout = mp5_subcalc(x1, x2, x3, x4, x5);
  *lout = mp5_subcalc(x5, x4, x3, x2, x1);
}
#undef MINMOD

//void reconstruct_lr_lin(double Ptmp[NMAX+2*NG][NVAR], int N,
//  double P_l[NMAX+2*NG][NVAR], double P_r[NMAX+2*NG][NVAR])
void reconstruct_lr_lin(struct FluidState *S, GridPrim Pl, GridPrim Pr, int dir)
{
  PLOOP {
	if (dir == 1) {
	  #pragma omp parallel for simd collapse(2)
	  KSLOOP(-1, N3) {
		JSLOOP(-1, N2) {
		  ISLOOP(-1, N1) {
			linear_mc(S->P[ip][k][j][i-1], S->P[ip][k][j][i], S->P[ip][k][j][i+1],
					  &(Pl[ip][k][j][i]), &(Pr[ip][k][j][i]));
		  }
		}
	  }
	} else if (dir == 2) {
	  #pragma omp parallel for simd collapse(2)
	  KSLOOP(-1, N3) {
		JSLOOP(-1, N2) {
		  ISLOOP(-1, N1) {
			linear_mc(S->P[ip][k][j-1][i], S->P[ip][k][j][i], S->P[ip][k][j+1][i],
					  &(Pl[ip][k][j][i]), &(Pr[ip][k][j][i]));
		  }
		}
	  }
	} else if (dir == 3) {
	  #pragma omp parallel for simd collapse(2)
	  KSLOOP(-1, N3) {
		JSLOOP(-1, N2) {
		  ISLOOP(-1, N1) {
			linear_mc(S->P[ip][k-1][j][i], S->P[ip][k][j][i], S->P[ip][k+1][j][i],
					  &(Pl[ip][k][j][i]), &(Pr[ip][k][j][i]));
		  }
		}
	  }
	}
  } // PLOOP
}


//void reconstruct_lr_par(double Ptmp[NMAX+2*NG][NVAR], int N,
//  double P_l[NMAX+2*NG][NVAR], double P_r[NMAX+2*NG][NVAR])
void reconstruct_lr_par(struct FluidState *S, GridPrim Pl, GridPrim Pr, int dir)
{
  PLOOP {
    if (dir == 1) {
      #pragma omp parallel for simd collapse(2)
      KSLOOP(-1, N3) {
        JSLOOP(-1, N2) {
          ISLOOP(-1, N1) {
            para(S->P[ip][k][j][i-2], S->P[ip][k][j][i-1], S->P[ip][k][j][i],
                 S->P[ip][k][j][i+1], S->P[ip][k][j][i+2], &(Pl[ip][k][j][i]),
                 &(Pr[ip][k][j][i]));
          }
        }
      }
    } else if (dir == 2) {
      #pragma omp parallel for simd collapse(2)
      KSLOOP(-1, N3) {
        JSLOOP(-1, N2) {
          ISLOOP(-1, N1) {
            para(S->P[ip][k][j-2][i], S->P[ip][k][j-1][i], S->P[ip][k][j][i],
                 S->P[ip][k][j+1][i], S->P[ip][k][j+2][i], &(Pl[ip][k][j][i]),
                 &(Pr[ip][k][j][i]));
          }
        }
      }
    } else if (dir == 3) {
      #pragma omp parallel for simd collapse(2)
      KSLOOP(-1, N3) {
        JSLOOP(-1, N2) {
          ISLOOP(-1, N1) {
            para(S->P[ip][k-2][j][i], S->P[ip][k-1][j][i], S->P[ip][k][j][i],
                 S->P[ip][k+1][j][i], S->P[ip][k+2][j][i], &(Pl[ip][k][j][i]),
                 &(Pr[ip][k][j][i]));
          }
        }
      }
    }
  } // PLOOP
}

//void reconstruct_lr_weno(double Ptmp[NMAX+2*NG][NVAR], int N,
//  double P_l[NMAX+2*NG][NVAR], double P_r[NMAX+2*NG][NVAR])
//{
//  ISLOOP(-1,N) {
//    PLOOP {
//      weno(Ptmp[i-2][ip],
//           Ptmp[i-1][ip],
//           Ptmp[i][ip],
//           Ptmp[i+1][ip],
//           Ptmp[i+2][ip],
//           &P_l[i][ip],
//           &P_r[i][ip]);
//    }
//  }
//}

void reconstruct_lr_weno(struct FluidState *S, GridPrim Pl, GridPrim Pr, int dir)
{
  PLOOP {
    if (dir == 1) {
      #pragma omp parallel for simd collapse(2)
      KSLOOP(-1, N3) {
        JSLOOP(-1, N2) {
          ISLOOP(-1, N1) {
            weno(S->P[ip][k][j][i-2], S->P[ip][k][j][i-1], S->P[ip][k][j][i],
                 S->P[ip][k][j][i+1], S->P[ip][k][j][i+2], &(Pl[ip][k][j][i]),
                 &(Pr[ip][k][j][i]));
          }
        }
      }
    } else if (dir == 2) {
      #pragma omp parallel for simd collapse(2)
      KSLOOP(-1, N3) {
        JSLOOP(-1, N2) {
          ISLOOP(-1, N1) {
            weno(S->P[ip][k][j-2][i], S->P[ip][k][j-1][i], S->P[ip][k][j][i],
                 S->P[ip][k][j+1][i], S->P[ip][k][j+2][i], &(Pl[ip][k][j][i]),
                 &(Pr[ip][k][j][i]));
          }
        }
      }
    } else if (dir == 3) {
      #pragma omp parallel for simd collapse(2)
      KSLOOP(-1, N3) {
        JSLOOP(-1, N2) {
          ISLOOP(-1, N1) {
            weno(S->P[ip][k-2][j][i], S->P[ip][k-1][j][i], S->P[ip][k][j][i],
                 S->P[ip][k+1][j][i], S->P[ip][k+2][j][i], &(Pl[ip][k][j][i]),
                 &(Pr[ip][k][j][i]));
          }
        }
      }
    }
  } // PLOOP
}

void reconstruct_lr_mp5(double Ptmp[NMAX+2*NG][NVAR], int N,
  double P_l[NMAX+2*NG][NVAR], double P_r[NMAX+2*NG][NVAR])
{
  ISLOOP(-1,N) {
    PLOOP {
      mp5(Ptmp[i-2][ip],
          Ptmp[i-1][ip],
          Ptmp[i][ip],
          Ptmp[i+1][ip],
          Ptmp[i+2][ip],
          &P_l[i][ip],
          &P_r[i][ip]);
    }
  }
}

//void reconstruct(double Ptmp[NMAX+2*NG][NVAR], int N,
//  double P_l[NMAX+2*NG][NVAR], double P_r[NMAX+2*NG][NVAR])
void reconstruct(struct FluidState *S, GridPrim Pl, GridPrim Pr, int dir)
{
  timer_start(TIMER_RECON);
  #if RECONSTRUCTION == LINEAR
    //reconstruct_lr_lin(Ptmp, N, P_l, P_r);
    reconstruct_lr_lin(S, Pl, Pr, dir);
  #elif RECONSTRUCTION == PPM
    //reconstruct_lr_par(Ptmp, N, P_l, P_r);
    reconstruct_lr_par(S, Pl, Pr, dir);
  #elif RECONSTRUCTION == WENO
    reconstruct_lr_weno(S, Pl, Pr, dir);
  #elif RECONSTRUCTION == MP5
    reconstruct_lr_mp5(Ptmp, N, P_l, P_r);
  #else
    fprintf(stderr, "Reconstruction algorithm not specified!\n");
    exit(-1);
  #endif
  timer_stop(TIMER_RECON);
}

