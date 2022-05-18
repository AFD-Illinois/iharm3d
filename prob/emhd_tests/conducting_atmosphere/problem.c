/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR HYDROSTATIC CONDUCTING ATMOSPHERE                   *
 *                                                                            *
 ******************************************************************************/

 #include "decs.h"
 #include "hdf5_utils.h"

 void set_problem_params()
{
  if (mpi_io_proc()) fprintf(stdout, "No problem specific params to set from param.dat!\n\n");
}

void save_problem_data(hid_t string_type){
  hdf5_write_single_val("emhd/conducting_atmosphere", "PROB", string_type);
}

// Set chi, nu, tau. Problem dependent
void set_emhd_parameters(struct GridGeom *G, struct FluidState *S, int i, int j, int k){
     
  //Initializations
  double rho = S->P[RHO][k][j][i];

  //set EMHD parameters based on closure relations
  double tau   = 10.;
  double kappa = 0.1;

  S->tau[k][j][i]      = tau;
  S->chi_emhd[k][j][i] = kappa / MY_MAX(SMALL, rho);
  S->nu_emhd[k][j][i]  = 0.;
 }

void init(struct GridGeom *G, struct FluidState *S) {
    
  double X[NDIM];

  // Grid parameters
  Rhor = (1. + sqrt(1. - a*a));

  // Set grid
  set_grid(G);
  LOG("Set grid");

  // load conducting atmosphere solution
  char fname_rCoords[STRLEN];
  sprintf(fname_rCoords, "atmosphere_soln_rCoords.txt");
  char fname_rho[STRLEN];
  sprintf(fname_rho, "atmosphere_soln_rho.txt");
  char fname_u[STRLEN];
  sprintf(fname_u, "atmosphere_soln_u.txt");
  char fname_q[STRLEN];
  sprintf(fname_q, "atmosphere_soln_phi.txt");

  FILE *fp_rCoords, *fp_rho, *fp_u, *fp_q;
  fp_rCoords = fopen(fname_rCoords, "r");
  fp_rho     = fopen(fname_rho, "r");
  fp_u       = fopen(fname_u, "r");
  fp_q       = fopen(fname_q, "r");

  // Load coordinates 'r' and compare against grid values
  double rCoords[N1 + 2*NG] = {0}, error = 0.;
  ILOOPALL {

    fscanf(fp_rCoords, "%lf", &(rCoords[i]));
    coord(i, NG, NG, CENT, X); // j and k don't matter since we need to compare only the radial coordinate
    error = fabs(X[1] - log(rCoords[i]));
    if (error > 1.e-10) { 
      fprintf(stdout, "Error at radial zone i = %d, iharm3d: %g, sage nb: %g\n", i, X[1], rCoords[i]);
    }
  }

  if (error > 1.e-10) exit(-1);

  // Loop over physical zones and initialize primitives
  // Temporaries
  double rho_temp, u_temp, q_temp;
  ILOOPALL {

    fscanf(fp_rho, "%lf", &(rho_temp));
    fscanf(fp_u,   "%lf", &(u_temp));
    fscanf(fp_q,   "%lf", &(q_temp));

    KLOOP {
      JLOOP {

        coord(i, j, k, CENT, X);
        S->P[RHO][k][j][i]           = rho_temp;
        S->P[UU][k][j][i]            = u_temp;
        S->P[U1][k][j][i]            = 0.;
        S->P[U2][k][j][i]            = 0.;
        S->P[U3][k][j][i]            = 0.;
        S->P[B1][k][j][i]            = 0.;
        S->P[B2][k][j][i]            = 0.;
        S->P[B3][k][j][i]            = 0.;
        S->P[Q_TILDE][k][j][i]       = q_temp;
        S->P[DELTA_P_TILDE][k][j][i] = 0.;

        // Note that the  velocity primitives defined up there isn't quite right.
        // For a fluid at rest wrt. the normal observer, ucon = {-1/g_tt,0,0,0}. We need to use this info
        // to obtain the correct values for U1, U2 and U3
        
        // ucon in BL
        double ucon[NDIM] = {0};
        ucon[0] = 1./sqrt(-G->gcov[CENT][0][0][j][i]);
        ucon[1] = 0.;
        ucon[2] = 0.;
        ucon[3] = 0.;

        double r = exp(X[1]);

        double trans[NDIM][NDIM], tmp[NDIM];
        double alpha, gamma, beta[NDIM];

        // Solve for v
        alpha = G->lapse[CENT][j][i];
        gamma = ucon[0]*alpha;

        beta[1] = alpha*alpha*G->gcon[CENT][0][1][j][i];
        beta[2] = alpha*alpha*G->gcon[CENT][0][2][j][i];
        beta[3] = alpha*alpha*G->gcon[CENT][0][3][j][i];

        S->P[U1][k][j][i] = ucon[1] + beta[1]*gamma/alpha;
        S->P[U2][k][j][i] = ucon[2] + beta[2]*gamma/alpha;
        S->P[U3][k][j][i] = ucon[3] + beta[3]*gamma/alpha;

       if (higher_order_terms == 1) {

         double rho   = S->P[RHO][k][j][i];
          double u     = S->P[UU][k][j][i];
          double Theta = (gam - 1.) * u / rho;

          set_emhd_parameters(G, S, i, j, k);
          double tau      = S->tau[k][j][i];
          double chi_emhd = S->chi_emhd[k][j][i];

          S->P[Q_TILDE][k][j][i] *= sqrt(tau / (chi_emhd * rho * pow(Theta, 2)));
       }
        
      }
    }

    // Save boundary values for Dirichlet boundary conditions
    if (i < NG) {
      S->P_BOUND[RHO][i]           = S->P[RHO][NG][NG][i];
      S->P_BOUND[UU][i]            = S->P[UU][NG][NG][i];
      S->P_BOUND[U1][i]            = S->P[U1][NG][NG][i]; // Since U1 is independent of j and k, we can use the value at any of the zones in those directions
      S->P_BOUND[U2][i]            = S->P[U2][NG][NG][i];
      S->P_BOUND[U3][i]            = S->P[U3][NG][NG][i];
      S->P_BOUND[B1][i]            = 0.;
      S->P_BOUND[B2][i]            = 0.;
      S->P_BOUND[B3][i]            = 0.;
      S->P_BOUND[Q_TILDE][i]       = S->P[Q_TILDE][NG][NG][i];
      S->P_BOUND[DELTA_P_TILDE][i] = 0.;
    }
    if (i > N1 + NG - 1) {
      S->P_BOUND[RHO][i-N1]           = S->P[RHO][NG][NG][i];
      S->P_BOUND[UU][i-N1]            = S->P[UU][NG][NG][i];
      S->P_BOUND[U1][i-N1]            = S->P[U1][NG][NG][i]; // Since U1 is independent of j and k, we can use the value at any of the zones in those directions
      S->P_BOUND[U2][i-N1]            = S->P[U2][NG][NG][i];
      S->P_BOUND[U3][i-N1]            = S->P[U3][NG][NG][i];
      S->P_BOUND[B1][i-N1]            = 0.;
      S->P_BOUND[B2][i-N1]            = 0.;
      S->P_BOUND[B3][i-N1]            = 0.;
     S->P_BOUND[Q_TILDE][i-N1]       = S->P[Q_TILDE][NG][NG][i];
     S->P_BOUND[DELTA_P_TILDE][i-N1] = 0.;
    }

  }

  fclose(fp_rCoords);
  fclose(fp_rho);
  fclose(fp_u);
  fclose(fp_q);

  //Enforce boundary conditions
  set_bounds(G, S);

}
