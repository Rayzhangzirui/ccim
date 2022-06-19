#ifndef INTERFACE_H
#define INTERFACE_H
#include "pb.h"
#include "grid.h"

// vector operation
double getdotprod(double *v, double *w, int thedim);
void normalize(double *v, double *w, int thedim);
void getcrossprod(double *v, double *w, double *z);
void getunitcrossprod(double *v, double *w, double *z);
void project(double *w, double *normal, double *v, int dim);


// interface info
char getinterfaceinfo(double &alpha, double *tangent, double *normal, double ***S,
                      int *index, int rstar, int sstar, GridData &grid);
char getinterfaceinfo(double &alpha, double *tangent, double *normal, double ***S,
                      int *index1, int *index2, GridData &grid);

char getinterfaceinfo(double &alpha, double *tangent1, double *tangent2, double *normal, 
                      double &sigma, double **Dn, double *Dsigma, double &jumpfe, 
                      double ***S, int *index1, int *index2, PBData &pb, GridData &grid);

void getinterfaceinfo(double *tangent1, double *tangent2, double &sigma, double **Dn,
                      double *Dsigma, double &jumpfe, int *index, int rstar, int sstar,
                      double alpha, double ***S, PBData &pb, GridData &grid);

void getinterfaceinfo(double *tangent1, double *tangent2, double &tau, double &sigma, 
                      double **Dn, double *Dsigma, double &jumpfe, double &aehere,
                      double &aethere, int *index, int rstar, int sstar, double alpha, 
                      double ***a, double ***S, PBData &pb, GridData &grid);


#endif