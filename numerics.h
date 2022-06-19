#ifndef NUMERICS_H
#define NUMERICS_H

#include "grid.h"

int gecp0(double **U, int PLR[], int PLC[], double **A, int n, int m);
void forwardbacksub0(double *x, double *b, double **U, int *PLR, int *PLC, int n);


double newtonweno(double ***u, int *index, int rstar, int sstar, double alpha, 
                  GridData &grid);
double newtonweno(double *u, double x0, double a, double dx, double nx);
double weno6interp(double ***u, int *index, int rstar, int sstar, double alpha, 
                   GridData &grid);
void dweno6interp1d(double &val, double &dval, double x, double *u, double a,
                    double dx, int nx);
double weno6interp1d(double x, double *u, double a, double dx, int nx);
double regulafalsiweno(double ***u, int *index, int rstar, int sstar, GridData &grid);
double regulafalsiweno(double *u, int i, int j, double a, double dx, double nx);
double lagrangeinterp1d(double x, double *u, double a, double dx, int nx);
double regulafalsi(double ***u, int *index, int rstar, int sstar, GridData &grid);
double regulafalsi(double *u, int i, int j, double a, double dx, double nx);
double regulafalsiexact(int *index, int rstar, int sstar, GridData &grid);
double bisection(double ***u, int *index, int rstar, int sstar, GridData &grid);
double bisection(double *u, int i, int j, double a, double dx, double nx);


void weno(double& dfp, double& dfn, double* fi, const long i,
          const double dx, const int order);

#endif