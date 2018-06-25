/*
 * pack.c
 *
 *  Created on: Jun 25, 2018
 *      Author: bprather
 */

#include "decs.h"

// Following double/float defs differ only in signature. C is a bear.
void pack_scalar_double(double in[N3+2*NG][N2+2*NG][N1+2*NG], double *out)
{
  int ind = 0;
  ZLOOP_OUT {
    out[ind] = in[k][j][i];
    ind++;
  }
}

void pack_scalar_float(double in[N3+2*NG][N2+2*NG][N1+2*NG], float *out)
{
  int ind = 0;
  ZLOOP_OUT {
    out[ind] = in[k][j][i];
    ind++;
  }
}

void pack_scalar_int(int in[N3+2*NG][N2+2*NG][N1+2*NG], int *out)
{
  int ind = 0;
  ZLOOP_OUT {
    out[ind] = in[k][j][i];
    ind++;
  }
}

void pack_vector_double(double in[][N3+2*NG][N2+2*NG][N1+2*NG], double *out, int vector_len)
{
  int ind = 0;
  for (int mu=0; mu < vector_len; mu++) {
    ZLOOP_OUT {
      out[ind] = in[mu][k][j][i];
      ind++;
    }
  }
}

void pack_vector_float(double in[][N3+2*NG][N2+2*NG][N1+2*NG], float *out, int vector_len)
{
  int ind = 0;
  ZLOOP_OUT {
    for (int mu=0; mu < vector_len; mu++) {
      out[ind] = in[mu][k][j][i];
      ind++;
    }
  }
}

// TODO more functions if I want to output F or other tensors
// Also: this wastes quite a bit of space writing all phi
void pack_Gtensor_double(double in[NDIM][NDIM][N2+2*NG][N1+2*NG], double *out)
{
  int ind = 0;
  ZLOOP_OUT {
    DLOOP2 {
	  out[ind] = in[nu][mu][j][i];
	  ind++;
    }
  }
}

void pack_Gtensor_float(double in[NDIM][NDIM][N2+2*NG][N1+2*NG], float *out)
{
  int ind = 0;
  ZLOOP_OUT {
    DLOOP2 {
	  out[ind] = in[nu][mu][j][i];
	  ind++;
    }
  }
}
