#ifndef INPUT_H
#define INPUT_H

#include "grid.h"

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


double getf(double *x, int thesign, PBData &pb, GridData &grid);
double getf(int *index, int rstar, int sstar, double alpha, int thesign, PBData &pb, 
            GridData &grid);
double getu(double *x, int thesign, GridData &grid);
double getu(int *index, int rstar, int sstar, double alpha, int thesign, GridData &grid);
double getDu(double *x, int s, int thesign, GridData &grid);
double getDu(int *index, int s, int rstar, int sstar, double alpha, int thesign, 
             GridData &grid);
double getD2u(double *x, int r, int s, int thesign, GridData &grid);
double getD2u(int *index, int r, int s, int rstar, int sstar, double alpha, 
              int thesign, GridData &grid);
//double getf(double ***f, int index, int rstar, int sstar, double alpha, int thesign, 
//GridData &grid);
double geta(double *x, double thesign, PBData& pb, GridData& grid);
void geta(double ***a, double ***u, double ***S, PBData &pb, GridData &grid);
double gettau(double *x, PBData &pb, GridData &grid);
void gettau(double &tau, int *index, int rstar, int sstar, double alpha, 
            GridData &grid);
void getDtau(double *Dtau, double *x, PBData &pb, GridData &grid);
void getDtau(double *Dtau, int *index, int rstar, int sstar, double alpha, 
             GridData &grid);
void getD2tau(double **D2tau, double *x, PBData &pb, GridData &grid);
void getD2tau(double **D2tau, int *index, int rstar, int sstar, double alpha, 
              GridData &grid);
double getsigma(double *x, double *normal, PBData &pb, GridData &grid);
void getsigma(double &sigma, int *index, int rstar, int sstar, double alpha, 
              double *normal, PBData &pb, GridData &grid);
void getDsigma(double *Dsigma, double *x, double *normal, double **Dnormal, PBData &pb, 
               GridData &grid);
void getDsigma(double *Dsigma, int *index, int rstar, int sstar, double alpha, 
              double *normal, double **Dnormal, PBData &pb, GridData &grid);

double getexactnormalTest(double *x, int r);

#endif