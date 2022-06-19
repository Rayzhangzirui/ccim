#ifndef PBDATA_H
#define PBDATA_H

#include "grid.h"
#include "matrix.h"

struct PBData
{
   double epsilonp, epsilonm;
   double gamma0;
   double ***psi;

   double **x;
   int N;
   double *Q;
   double beta;
   double *c;
   double *epsilone;
   int dim;

   double rho0;
   double *epsilonlj;
   double *sigmalj;
   double ***LJ;
};


void init_PBData(PBData &pb, GridData &grid, double epsp, double epsm);

#endif