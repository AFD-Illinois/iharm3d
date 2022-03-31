/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR AN EMHD SHOCK                                       *
 *                                                                            *
 ******************************************************************************/

 #include "decs.h"
 #include "hdf5_utils.h"

void set_problem_params() {
 if (mpi_io_proc()) fprintf(stdout, "No problem specific params to set from param.dat!\n\n");
}

void save_problem_data(hid_t string_type){
  hdf5_write_single_val("emhd/shock_test", "PROB", string_type);
}

// Set chi, nu, tau. Problem dependent
void set_emhd_parameters(struct GridGeom *G, struct FluidState *S, int i, int j, int k){
     
  // Initializations
  double rho = S->P[RHO][k][j][i];
  double u   = S->P[UU][k][j][i];
  double P   = (gam - 1.) * u;

  // sound speed
  double cs2 = (gam * P) / (rho + (gam * u));

  // set EMHD parameters based on closure relations
  double tau = 0.1;
  S->tau[k][j][i]      = tau;
  S->chi_emhd[k][j][i] = conduction_alpha * cs2 * tau;
  S->nu_emhd[k][j][i]  = viscosity_alpha * cs2 * tau;
}

void init(struct GridGeom *G, struct FluidState *S) {
    
  double X[NDIM];

  // Set grid
  set_grid(G);
  LOG("Set grid");

  // Grid center
  double x1_center = (x1Min + x1Max) / 2.;

  // Left and right states
  double rhoL = 1.,     rhoR = 3.08312999;
  double uL   = 1.,     uR   = 4.94577705;
  double u1L  = 1.,     u1R  = 0.32434571;
  double u2L  = 0.,     u2R  = 0.;
  double u3L  = 0.,     u3R  = 0.;
  double B1L  = 1.e-5,  B1R  = 1.e-5;
  double B2L  = 0,      B2R  = 0.;
  double B3L  = 0.,     B3R  = 0.;

  // load BVP solution

  char fname_rho[STRLEN];
  sprintf(fname_rho, "shock_soln_rho.txt");
  char fname_u[STRLEN];
  sprintf(fname_u, "shock_soln_u.txt");
  char fname_u1[STRLEN];
  sprintf(fname_u1, "shock_soln_u1.txt");
  char fname_q[STRLEN];
  sprintf(fname_q, "shock_soln_q.txt");
  char fname_dP[STRLEN];
  sprintf(fname_dP, "shock_soln_dP.txt");

  FILE *fp_rho, *fp_u, *fp_u1, *fp_q, *fp_dP;
  fp_rho = fopen(fname_rho, "r");
  fp_u   = fopen(fname_u, "r");
  fp_u1  = fopen(fname_u1, "r");
  fp_q   = fopen(fname_q, "r");
  fp_dP  = fopen(fname_dP, "r");

  // Loop over physical zones and initialize primitives
  ZLOOP {
    coord(i, j, k, CENT, X);

    bool lhs = X[1] < x1_center;

    // Initialize primitives
    fscanf(fp_rho, "%lf", &(S->P[RHO][k][j][i]));
    fscanf(fp_u, "%lf", &(S->P[UU][k][j][i]));
    fscanf(fp_u1, "%lf", &(S->P[U1][k][j][i]));
    S->P[U2][k][j][i]  = (lhs) ? u2L  : u2R;
    S->P[U3][k][j][i]  = (lhs) ? u3L  : u3R;
    S->P[B1][k][j][i]  = (lhs) ? B1L  : B1R;
    S->P[B2][k][j][i]  = (lhs) ? B2L  : B2R;
    S->P[B3][k][j][i]  = (lhs) ? B3L  : B3R;
    fscanf(fp_q, "%lf", &(S->P[Q_TILDE][k][j][i]));
    fscanf(fp_dP, "%lf", &(S->P[DELTA_P_TILDE][k][j][i]));

    if (higher_order_terms == 1) {

      double rho = S->P[RHO][k][j][i];
      double u   = S->P[UU][k][j][i];
      double Theta = (gam - 1.) * u / rho;

      set_emhd_parameters(G, S, i, j, k);
      double tau = S->tau[k][j][i];
      double chi_emhd = S->chi_emhd[k][j][i];
      double nu_emhd  = S->nu_emhd[k][j][i];

      S->P[Q_TILDE][k][j][i]       *= sqrt(tau / (chi_emhd * rho * pow(Theta, 2)));
      S->P[DELTA_P_TILDE][k][j][i] *= sqrt(tau / (nu_emhd * rho * Theta));
    }
  }

  fclose(fp_rho);
  fclose(fp_u);
  fclose(fp_u1);
  fclose(fp_q);
  fclose(fp_dP);

  //Enforce boundary conditions
  set_bounds(G, S);
}
