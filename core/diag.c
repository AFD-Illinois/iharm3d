/******************************************************************************
 *                                                                            *
 * DIAG.C                                                                     *
 *                                                                            *
 * DIAGNOSTIC OUTPUT                                                          *
 *                                                                            *
 ******************************************************************************/

#include "decs.h"

void reset_log_variables()
{
  #if RADIATION
  step_tot = step_made = step_abs = step_scatt = step_rec = 0;
  #endif
}

// Evaluate flux based diagnostics; put results in global variables
void diag_flux(struct FluidFlux *F)
{
  mdot = edot = ldot = 0.;
  if (global_start[1] == 0) {
    #pragma omp parallel for \
      reduction(+:mdot) reduction(-:edot) reduction(+:ldot) \
      collapse(2)
    JSLOOP(0, N2 - 1) {
      KSLOOP(0, N3 - 1) {
        mdot -= F->X1[RHO][k][j][NG]*dx[2]*dx[3];
        edot -= (F->X1[UU][k][j][NG] - F->X1[RHO][k][j][NG])*dx[2]*dx[3];
        ldot += F->X1[U3][k][j][NG]*dx[2]*dx[3];
      }
    }
  }
  mdot = mpi_io_reduce(mdot);
  edot = mpi_io_reduce(edot);
  ldot = mpi_io_reduce(ldot);
}

void diag(struct GridGeom *G, struct FluidState *S, int call_code)
{
  static FILE *ener_file;

  if (call_code == DIAG_INIT) {
    // Set things up
    if(mpi_io_proc()) {
      ener_file = fopen("dumps/log.out", "a"); //TODO make this more flexible
      if (ener_file == NULL) {
        fprintf(stderr,
          "error opening energy output file\n");
        exit(1);
      }
    }
  }

  double pp = 0.;
  double divbmax = 0.;
  double rmed = 0.;
  double e = 0.;
  double Phi = 0;

  // Calculate conserved quantities
  if ((call_code == DIAG_INIT || call_code == DIAG_LOG ||
       call_code == DIAG_FINAL) && !failed) {

    get_state_vec(G, S, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1);
    prim_to_flux_vec(G, S, 0, CENT, 0, N3 - 1, 0, N2 - 1, 0, N1 - 1, S->U);

    ZLOOP {
      rmed += S->U[RHO][k][j][i]*dV;
      pp += S->U[U3][k][j][i]*dV;
      e += S->U[UU][k][j][i]*dV;

      double divb = flux_ct_divb(G, S, i, j, k);

      if (divb > divbmax) {
        divbmax = divb;
      }

      if (global_start[1] == 0 && i == 5+NG) {
        Phi += 0.5*fabs(S->P[B1][k][j][i])*dx[2]*dx[3]*G->gdet[j][i][CENT];
      }
    }
  }

  rmed = mpi_io_reduce(rmed);
  pp = mpi_io_reduce(pp);
  e = mpi_io_reduce(e);
  divbmax = mpi_io_max(divbmax);
  if (mdot > SMALL) {
    Phi = mpi_io_reduce(Phi);
    Phi /= sqrt(mdot);
  } else {
      Phi = 0;
  }

  if ((call_code == DIAG_INIT && !is_restart) ||
    call_code == DIAG_DUMP || call_code == DIAG_FINAL) {
    dump(G, S);
    dump_cnt++;
  }

  //if (call_code == DIAG_FINAL) {
  //  dump(G, S);
  //}

  if (call_code == DIAG_INIT || call_code == DIAG_LOG ||
      call_code == DIAG_FINAL) {
    if(mpi_io_proc()) {
      fprintf(stdout, "LOG      t=%g \t divbmax: %g\n",
        t,divbmax);
      fprintf(ener_file, "%10.5g %10.5g %10.5g %10.5g ", t, rmed, pp, e);
      //fprintf(ener_file, "%10.5g %10.5g %10.5g %10.5g %15.8g %15.8g ",
      //  t, rmed, pp, e,
      //  P[N1/2][N2/2][N3/2][UU]*pow(P[N1/2][N2/2][N3/2][RHO], -gam),
      //  P[N1/2][N2/2][N3/2][UU]);
      // REEVALUATE MDOT LOCALLY!
      fprintf(ener_file, "%15.8g %15.8g %15.8g %15.8g ", mdot, edot, ldot, Phi);
      fprintf(ener_file, "%15.8g ", divbmax);
      fprintf(ener_file, "\n");
      fflush(ener_file);
    }
  }
}

// Diagnostic routines
void fail(int fail_type, int i, int j, int k)
{
  failed = 1;

  fprintf(stderr, "\n\n[%d %d %d] FAIL: %d\n", i, j, k, fail_type);

  //area_map(i, j, k, P);

  //diag(DIAG_FINAL);

  exit(0);
}

// Map out region around failure point
/*void area_map(int i, int j, int k, grid_prim_type prim)
{
  fprintf(stderr, "*** AREA MAP ***\n");

  PLOOP {
    fprintf(stderr, "variable %d \n", ip);
    fprintf(stderr, "i = \t %12d %12d %12d\n", i - 1, i,
      i + 1);
    fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j + 1,
      prim[i - 1][j + 1][k][ip], prim[i][j + 1][k][ip],
      prim[i + 1][j + 1][k][ip]);
    fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j,
      prim[i - 1][j][k][ip], prim[i][j][k][ip],
      prim[i + 1][j][k][ip]);
    fprintf(stderr, "j = %d \t %12.5g %12.5g %12.5g\n", j - 1,
      prim[i - 1][j - 1][k][ip], prim[i][j - 1][k][ip],
      prim[i + 1][j - 1][k][ip]);
  }

  fprintf(stderr, "****************\n");
}*/

// Evaluate flux based diagnostics; put results in global variables
/*void diag_flux(grid_prim_type F1, grid_prim_type F2, grid_prim_type F3)
{
  mdot = edot = ldot = 0.;
  #pragma omp parallel for \
    reduction(+:mdot) reduction(-:edot) reduction(+:ldot) \
    collapse(2)
  JSLOOP(0, N2 - 1) {
    KSLOOP(0, N3 - 1) {
      mdot += F1[0][j][k][RHO]*dx[2]*dx[3];
      edot -= (F1[0][j][k][UU] - F1[0][j][k][RHO])*dx[2]*dx[3];
      ldot += F1[0][j][k][U3]*dx[2]*dx[3];
    }
  }
}*/

double flux_ct_divb(struct GridGeom *G, struct FluidState *S, int i, int j,
  int k)
{
  #if N3 > 1
  if(i > 0 + NG && j > 0 + NG && k > 0 + NG &&
     i < N1 + NG && j < N2 + NG && k < N3 + NG) {
  #elif N2 > 1
  if(i > 0 + NG && j > 0 + NG &&
     i < N1 + NG && j < N2 + NG) {
  #elif N1 > 1
  if(i > 0 + NG &&
     i < N1 + NG) {
  #else
  if (0) {
  #endif
    return fabs(0.25*(
      S->P[B1][k][j][i]*G->gdet[CENT][j][i]
      + S->P[B1][k][j-1][i]*G->gdet[CENT][j-1][i]
      + S->P[B1][k-1][j][i]*G->gdet[CENT][j][i]
      + S->P[B1][k-1][j-1][i]*G->gdet[CENT][j-1][i]
      - S->P[B1][k][j][i-1]*G->gdet[CENT][j][i-1]
      - S->P[B1][k][j-1][i-1]*G->gdet[CENT][j-1][i-1]
      - S->P[B1][k-1][j][i-1]*G->gdet[CENT][j][i-1]
      - S->P[B1][k-1][j-1][i-1]*G->gdet[CENT][j-1][i-1]
      )/dx[1] +
      0.25*(
      S->P[B2][k][j][i]*G->gdet[CENT][j][i]
      + S->P[B2][k][j][i-1]*G->gdet[CENT][j][i-1]
      + S->P[B2][k-1][j][i]*G->gdet[CENT][j][i]
      + S->P[B2][k-1][j][i-1]*G->gdet[CENT][j][i-1]
      - S->P[B2][k][j-1][i]*G->gdet[CENT][j-1][i]
      - S->P[B2][k][j-1][i-1]*G->gdet[CENT][j-1][i-1]
      - S->P[B2][k-1][j-1][i]*G->gdet[CENT][j-1][i]
      - S->P[B2][k-1][j-1][i-1]*G->gdet[CENT][j-1][i-1]
      )/dx[2] +
      0.25*(
      S->P[B3][k][j][i]*G->gdet[CENT][j][i]
      + S->P[B3][k][j-1][i]*G->gdet[CENT][j-1][i]
      + S->P[B3][k][j][i-1]*G->gdet[CENT][j][i-1]
      + S->P[B3][k][j-1][i-1]*G->gdet[CENT][j-1][i-1]
      - S->P[B3][k-1][j][i]*G->gdet[CENT][j][i]
      - S->P[B3][k-1][j-1][i]*G->gdet[CENT][j-1][i]
      - S->P[B3][k-1][j][i-1]*G->gdet[CENT][j][i-1]
      - S->P[B3][k-1][j-1][i-1]*G->gdet[CENT][j-1][i-1]
      )/dx[3]);
  } else {
    return 0.;
  }
}
