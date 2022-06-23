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


double getpsivac(double *x, PBData &pb);
void getDpsivac(double *Dpsi, double *x, PBData &pb);
double getDpsivacn(int *index, double ***S, PBData &pb, GridData &grid);
void getD2psivac(double **D2psi, double *x, PBData &pb);
double Bval(double s, PBData &pb);
double Bprime(double s, PBData &pb);
double B2prime(double s, PBData &pb);


#endif