/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR AN EMHD SHOCK                                       *
 *                                                                            *
 ******************************************************************************/

 #include "decs.h"
 #include <complex.h>
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

    // Loop over physical zones and initialize primitives
    ZLOOP {
        coord(i, j, k, CENT, X);

        bool lhs = X[1] < x1_center;

        // Initialize primitives
        S->P[RHO][k][j][i] = (lhs) ? rhoL : rhoR;
        S->P[UU][k][j][i]  = (lhs) ? uL   : uR;
        S->P[U1][k][j][i]  = (lhs) ? u1L  : u1R;
        S->P[U2][k][j][i]  = (lhs) ? u2L  : u2R;
        S->P[U3][k][j][i]  = (lhs) ? u3L  : u3R;
        S->P[B1][k][j][i]  = (lhs) ? B1L  : B1R;
        S->P[B2][k][j][i]  = (lhs) ? B2L  : B2R;
        S->P[B3][k][j][i]  = (lhs) ? B3L  : B3R;
        S->P[Q_TILDE][k][j][i]       = 0.;
        S->P[DELTA_P_TILDE][k][j][i] = 0.;
    }

    //Enforce boundary conditions
    set_bounds(G, S);
}
