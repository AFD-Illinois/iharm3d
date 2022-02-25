/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR 2D EMHD WAVE                                        *
 *                                                                            *
 ******************************************************************************/

 #include "decs.h"
 #include <complex.h>
 #include "hdf5_utils.h"

 void set_problem_params()
{
  set_param("rho_floor_fluid_element", &rho_floor_fluid_element);
  set_param("uu_floor_fluid_element", &uu_floor_fluid_element);
  set_param("bsq_floor_fluid_element", &bsq_floor_fluid_element);
  set_param("Theta_floor_fluid_element", &Theta_floor_fluid_element);
}

void save_problem_data(hid_t string_type){
  hdf5_write_single_val("emhd/anisotropic_conduction", "PROB", string_type);
  hdf5_write_single_val(&rho_floor_fluid_element, "rho_floor_fluid_element", H5T_IEEE_F64LE);
  hdf5_write_single_val(&uu_floor_fluid_element, "uu_floor_fluid_element", H5T_IEEE_F64LE);
  hdf5_write_single_val(&bsq_floor_fluid_element, "bsq_floor_fluid_element", H5T_IEEE_F64LE);
  hdf5_write_single_val(&Theta_floor_fluid_element, "Theta_floor_fluid_element", H5T_IEEE_F64LE);
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
    S->chi_emhd[k][j][i] = 0.01;
    S->nu_emhd[k][j][i]  = 0.;
 }

void init(struct GridGeom *G, struct FluidState *S) {
    
    double cos_phi, sin_phi, X[NDIM];

    // Set grid
    set_grid(G);
    LOG("Set grid");

    // Anisotropic conduction parameters
    double A = 0.2;
    double R = sqrt(0.005);
    double B0 = 1.e-4;
    double k = 4.;

    // Loop over physical zones and initialize primitives
    ZLOOP {
        coord(i, j, k, CENT, X);

        // Local variables to improve readability later
        double r = sqrt(pow((X[1] - 0.5), 2) + pow((X[2] - 0.5), 2));

        // Initialize primitives
        S->P[RHO][k][j][i] = 1 - (A * exp(-pow(r,2) / pow(R,2)));
        S->P[UU][k][j][i]  = 1.;
        S->P[U1][k][j][i]  = 0.;
        S->P[U2][k][j][i]  = 0.;
        S->P[U3][k][j][i]  = 0.;
        S->P[B1][k][j][i]  = B0;
        S->P[B2][k][j][i]  = B0 * sin(2*M_PI*k*X[1]);
        S->P[B3][k][j][i]  = 0.;
        S->P[Q_TILDE][k][j][i]       = 0.;
        S->P[DELTA_P_TILDE][k][j][i] = 0.;
    }

    //Enforce boundary conditions
    set_bounds(G, S);
}
