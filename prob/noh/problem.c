/******************************************************************************
 *                                                                            *
 * PROBLEM.C                                                                  *
 *                                                                            *
 * INITIAL CONDITIONS FOR 1D (NOH) SHOCK TEST                                 *
 *                                                                            *
 ******************************************************************************/

 #include "decs.h"
 #include "hdf5_utils.h"

 static double mach;
 static double rhoL;
 static double rhoR;
 static double PL;
 static double PR;
 static int set_tlim;

 void set_problem_params() {
     set_param("mach", &mach);
     set_param("rhoL", &rhoL);
     set_param("rhoR", &rhoR);
     set_param("PL", &PL);
     set_param("PR", &PR);
     set_param("set_tlim", &set_tlim);
 }

 void save_problem_data(hid_t string_type) {
    hdf5_write_single_val(&fel_constant, "fel_constant", H5T_IEEE_F64LE);
    hdf5_write_single_val(&mach, "mach", H5T_IEEE_F64LE);
    hdf5_write_single_val(&rhoL, "rhoL", H5T_IEEE_F64LE);
    hdf5_write_single_val(&rhoR, "rhoR", H5T_IEEE_F64LE);
    hdf5_write_single_val(&PL, "PL", H5T_IEEE_F64LE);
    hdf5_write_single_val(&PR, "PR", H5T_IEEE_F64LE);
    hdf5_write_single_val(&set_tlim, "set_tlim", H5T_STD_I32LE);
 }

 void init(struct GridGeom *G, struct FluidState *S) {
    double X[NDIM];

    // Modify final time according to Ressler+ (2015)
    double cs2 = (gam * (gam - 1) * PL) / rhoL;
    double v1 = mach * sqrt(cs2);
    if (set_tlim == 1) tf = 0.6*(x1Max - x1Min)/v1;

    set_grid(G);
    LOG("Set grid");

    double center = (x1Min + x1Max) / 2.;

    // Compute gamma 
    double gamma = 1. / sqrt(1. - v1 * v1); // Since we are in flat space

    ZLOOP {
        coord(i, j, k, CENT, X);
        bool lhs = X[1] < center;
        S->P[RHO][k][j][i] = (lhs) ? rhoL : rhoR;
        S->P[UU][k][j][i] = ((lhs) ? PL : PR)/(gam - 1.);
        S->P[U1][k][j][i] = ((lhs) ? v1 : -v1) * gamma;
        S->P[U2][k][j][i] = 0.;
        S->P[U3][k][j][i] = 0.;
        S->P[B1][k][j][i] = 0.;
        S->P[B2][k][j][i] = 0.;
        S->P[B3][k][j][i] = 0.;
    }

    #if ELECTRONS
        init_electrons(G,S);
    #endif

    set_bounds(G, S);
    LOG("Finished init()");
 }