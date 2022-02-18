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

 static double amp, k1, k2, real_omega, imag_omega;

 void set_problem_params()
{
    set_param("amp", &amp);
    set_param("k1", &k1);
    set_param("k2", &k2);
    set_param("real_omega", &real_omega);
    set_param("imag_omega", &imag_omega);
}

void save_problem_data(hid_t string_type){
    hdf5_write_single_val("emhd/linear_modes", "PROB", string_type);
    hdf5_write_single_val(&amp, "amp", H5T_IEEE_F64LE);
    hdf5_write_single_val(&k1, "k1", H5T_IEEE_F64LE);
    hdf5_write_single_val(&k2, "k2", H5T_IEEE_F64LE);   
    hdf5_write_single_val(&real_omega, "real_omega", H5T_IEEE_F64LE);
    hdf5_write_single_val(&imag_omega, "imag_omega", H5T_IEEE_F64LE);
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
    double tau = 1.;
    S->tau[k][j][i]      = tau;
    S->chi_emhd[k][j][i] = cs2 * tau;
    S->nu_emhd[k][j][i]  = cs2 * tau;
 }

void init(struct GridGeom *G, struct FluidState *S) {
    
    double cos_phi, sin_phi, X[NDIM];

    // Set grid
    set_grid(G);
    LOG("Set grid");

    // Mean state
    double rho0 = 1.;
    double u0   = 2.;
    double u10  = 0.;
    double u20  = 0.;
    double u30  = 0.;
    double B10  = 0.1;
    double B20  = 0.3;
    double B30  = 0.;
    // NOTE: The following are NOT the relaxed values!! These are just background/mean state values for linear analysis
    double q0   = 0.;
    double delta_p0 = 0.;

    // Loop over physical zones and initialize primitives
    ZLOOP {
        coord(i, j, k, CENT, X);

        // Local variables to improve readability later
        cos_phi = cos(k1*X[1] + k2*X[2]);
        sin_phi = sin(k1*X[1] + k2*X[2]);

        // Perturbations
        double drho     = amp * (((-0.518522524082246)*cos_phi) + ((0.1792647678001878)*sin_phi));
        double du       = amp * ((0.5516170736393813)*cos_phi);
        double du1      = amp * (((0.008463122479547856)*cos_phi) + ((-0.011862022608466367)*sin_phi));
        double du2      = amp * (((-0.16175466371870734)*cos_phi) + ((0.034828080823603294)*sin_phi));
        double du3      = 0.;
        double dB1      = amp * (((-0.05973794979640743)*cos_phi) + ((0.03351707506150924)*sin_phi));
        double dB2      = amp * (((0.02986897489820372)*cos_phi) - ((0.016758537530754618)*sin_phi));
        double dB3      = 0.;
        double dq       = amp * (((0.5233486841539436)*cos_phi) - ((0.04767672501939603)*sin_phi));
        double ddelta_p = amp * (((0.2909106062057657)*cos_phi) - ((0.02159452055336572)*sin_phi));

        // NOTE: Will have to edit this when higher order terms are included

        // Initialize primitives
        S->P[RHO][k][j][i] = rho0 + drho;
        S->P[UU][k][j][i]  = u0 + du;
        S->P[U1][k][j][i]  = u10 + du1;
        S->P[U2][k][j][i]  = u20 + du2;
        S->P[U3][k][j][i]  = u30 + du3;
        S->P[B1][k][j][i]  = B10 + dB1;
        S->P[B2][k][j][i]  = B20 + dB2;
        S->P[B3][k][j][i]  = B30 + dB3;
        S->P[Q_TILDE][k][j][i]       = q0 + dq;
        S->P[DELTA_P_TILDE][k][j][i] = delta_p0 + ddelta_p;
    }

    //Enforce boundary conditions
    set_bounds(G, S);
}