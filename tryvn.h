#ifndef TRYVN_H
#define TRYVN_H

#include <math.h>
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <string>
#include <vector>

#include "grid.h"
#include "matrix.h"

#include "input.h"
#include "solver.h"
#include "hypresolver.h"
#include "amg.h"
#include "march.h"
#include "storage.h"
#include "numerics.h"
#include "ccim.h"
#include "icim.h"
#include "finitediff.h"

// using namespace std;
using std::ofstream;
static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2))









void output(ofstream& outfile, double ***phi, int *nx);
void output(ofstream& outfile, double **phi, int *nx);
void output(ofstream& outfile, double **phi, char **tube, int *nx);
void output(ofstream& outfile, double **u, double **S, int *nx);
void outputneighbors(double ***S, int *index, GridData &grid);

// void getinit(double*** &S, PBData &pb, MarchStruct &march, TempStruct &tmp, 
//              GridData &grid);


void interiorpt2(SparseElt2**** &A, double ***b, int *index, double ***S, PBData &pb, 
                 GridData &grid);
void interiordirection2(SparseElt2**** &A, double ***b, int *index, int r, double ***S, 
                        PBData &pb, GridData &grid);
void interiorptsmall(SparseElt2**** &A, double ***b, int *index, double ***S, 
                     PBData &pb, GridData &grid);
void cim1again2(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
                PBData &pb, GridData &grid);
void cim1again3(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
                PBData &pb, GridData &grid);



void perturbstatus(double ***S, double tol, int maxsteps, PBData &pb, GridData &grid);


char perturbstatuselt(double ***S, int *index, double tol, PBData &pb, GridData &grid);





void outputstatus(ofstream& outfile, double ***phi, GridData &grid);
void cim2again3(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
                PBData &pb, GridData &grid);

void cim3(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
          PBData &pb, GridData &grid);
void cim3again(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
               PBData &pb, GridData &grid);
void cim4(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
          PBData &pb, GridData &grid);
void cim5(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
          PBData &pb, GridData &grid);
void cim345(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
            PBData &pb, GridData &grid);
void cim345(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &smallsize, 
            int *index, int gamma[][2], double ***S, PBData &pb, GridData &grid);


void getD1(double ****D1, int *sk, GridData grid);
void getD2(double ***D2[][3], int sk2[][3][4], GridData grid);
void getD2(double ***D2, int m, int n, int sk2[][3][4], GridData grid);
void getD2(double ***D2, int m, int n, int sk2[][3][4], int *tindex, GridData grid);
void getD2(double ***D2, double &jumpuxxcoef, int m, int n, int rstar, int sstar, 
           double thesign, int sk2[][3][4], GridData grid);


void getjumpux(double& u0, double ***ucoef, double *uxxcoef, int *index, int rstar, 
               int *sk, double alpha, int thesign, double *normal, double *tangent, 
               double ****D1ucoef, double **D1uxxcoef, double ***S, PBData &pb, 
               GridData &grid);
void getjumpux(double& u0, double *uxcoef, double **uxxcoef, int *index, int rstar, 
               int sk, double alpha, int thesign, double *normal, double *tangent, 
               double **D1uxcoef, double ***D1uxxcoef, double ***S, PBData &pb, 
               GridData &grid);
void getjumpux(double& u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
               int *index, int rstar, int sk, double alpha, int thesign, 
               double *normal, double *tangent, int mid, double ****D1ucoef, 
               double **D1uxcoef, double **D1uxxcoef, double ***S, PBData &pb, 
               GridData &grid);
void recast(double &u0, double ***ucoef, double *uxcoef, double **uxxcoef, 
            double **jumpD2u, double *****jumpD2ucoef, double ***jumpD2uxcoef, 
            double ***jumpD2uxxcoef, double ***D2[][3], double D2jumpuxxcoef[][3],
            GridData &grid);
void recast(double ***ucoef, double *uxcoef, double **uxxcoef, double ***D2[][3], 
            double *D2uxcoef[][3], double *D2uxxcoef[][3], int mid, GridData &grid);
void getjumpuxx(double& u0, double ***ucoef, double *uxxcoef, int *index, int rstar, 
                int *sk, double alpha, int thesign, double *normal, double ****D1ucoef,
                double **D1uxxcoef, double ****D1, double*** D2[][3], double ***S, 
                PBData &pb, GridData &grid);
void getjumpuxx(double& u0, double ***ucoef, double *uxxcoef, int *index, int rstar, 
                int *sk, double alpha, int thesign, double *normal, double ****D1ucoef,
                double **D1uxxcoef, double*** D2[][3], double ***S, PBData &pb, 
                GridData &grid);
void getjumpuxx(double& u0, double ***ucoef, double *uxxcoef, int *index, int rstar, 
                int *sk, double alpha, int thesign, double *normal, double ****D1,
                double*** D2[][3], double ***S, PBData &pb, GridData &grid);
void getjumpuxx(double **u0, double *****ucoef, double ***uxcoef, double ***uxxcoef, 
                int *index, int rstar, int sk, double alpha, int thesign, 
                double *normal, double **D1uxcoef, double ***D1uxxcoef, 
                double*** D2[][3], double D2jumpuxxcoef[][3], double ***S, 
                PBData &pb, GridData &grid);
void getjumpuxx(double &u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                int *index, int rstar, int sk, double alpha, int thesign, 
                double *normal, int mid, double ****D1ucoef, double **D1uxcoef, 
                double **D1uxxcoef, double*** D2[][3], double *D2uxcoef[][3], 
                double *D2uxxcoef[][3], double ***S, PBData &pb, GridData &grid);


void getDu(double ****ucoef, double **uxxcoef, int *index, int rstar,
           int sstar, double alpha, int thesign, int *sk, double ****D1, 
           double*** D2[][3], GridData grid);
void getDu(double **uxcoef, double ***uxxcoef, int *index, int rstar, int sstar, 
           double alpha, int thesign, GridData grid);
void getDu(double *Du, int *index, int rstar, int sstar, double alpha, int thesign, 
           int *sk, double ****D1, double*** D2[][3], GridData grid);

int getstatus4(double ***S, int *index, GridData grid);
// double getinterfacegrad3(double *grad, double ***u, double ***S, int *index, int rstar, 
//                          int sstar, PBData &pb, GridData grid);


void linearsystembasic(SparseElt2**** &A, double ***b, double ***S, GridData &grid);
void linearsystem3(SparseElt2**** &A, double ***b, double ***S, PBData &pb, 
                   GridData &grid);
void linearsystem4(SparseElt2**** &A, double ***b, double ***S, PBData &pb, 
                   GridData &grid);
void linearsystem5(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, 
                   int &smallsize, double ***S, PBData &pb, GridData &grid);
void linearsystem6(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, 
                   int &smallsize, double ***a, double ***S, PBData &pb, 
                   GridData &grid);
void linearsystemZL(SparseElt2**** &A, double ***b, double ***S, PBData &pb, 
                    GridData &grid);
int checkcim3(double ***S, int *index, GridData &grid);

char checkcimstatus(double ***S, GridData &grid);



void getRHS(double ***rhs, double ***u, PBData &pb, MarchStruct &march, TempStruct &tmp, 
            GridData &grid);
void getRHS(double ***rhs, double ***u, double ***a, PBData &pb, MarchStruct &march, 
            TempStruct &tmp, GridData &grid);
void getRHShyp(double ***rhs, double ***u, PBData &pb, MarchStruct &march, 
               TempStruct &tmp, GridData &grid);
void getRHSexactvn(double ***rhs, double ***u, PBData &pb, MarchStruct &march, 
                   TempStruct &tmp, GridData &grid);
void getRHScurv(double ***rhs, double ***u, TempStruct &tmp, GridData &grid);




void checkDanswer(double ***u, double ***S, PBData &pb, GridData &grid);
void checkDanswer(double *uxerr, int *index, int gamma[][2], double ***u, double ***S, 
                  PBData &pb, GridData &grid);
void checkDanswer2(double ***u, double ***S, PBData &pb, GridData &grid);
void checkDanswer2(double *uxerr, int *index, int gamma[][2], double ***u, double ***S, 
                   PBData &pb, GridData &grid);
void flipvector(double *v, double *w, int *index, int rstar, int sstar, double alpha,
                double *normal, double ***S, PBData &pb, GridData &grid);
void flipvector2(double *v, double *w, int *index, int rstar, int sstar, double alpha,
                 double *normal, double ***S, PBData &pb, GridData &grid);
void getcim3Du(double &uint, double *Du, int *zindex, int rstar, int sstar, 
               double ***u, double ***S, PBData &pb, GridData &grid);
void checkcim3Du(double ***u, double ***S, PBData &pb, GridData &grid);
void checkallDu(double ***u, double ***S, PBData &pb, GridData &grid);
void checkcim5Du(double ***u, double ***S, PBData &pb, GridData &grid);
void checkcim345Du(double ***u, double ***S, PBData &pb, GridData &grid);

void checkcim345Dupm(double ***u, double ***S, PBData &pb, GridData &grid);
void checkcim345jumpDu(double ***u, double ***S, PBData &pb, GridData &grid);
void getcim4Du(double &uint, double *Du, int *index, int rstar, int sstar, double ***u, 
               double ***S, PBData &pb, GridData &grid);
void getcim5Du(double &uint, double *Du, int *index, int rstar, int sstar, double ***u, 
               double ***S, PBData &pb, GridData &grid);
void getinterfaceDu(double &uint, double *Du, int *index, int rstar, int sstar,
                    double ***u, double ***S, PBData &pb, GridData &grid);
void cimexact(SparseElt2**** &A, double ***b, int *index, double ***S, PBData &pb, 
              GridData &grid);



void getcim1Du(double &uint, double *Du, int *index, int rstar, int sstar, double ***u,
               int gamma[][2], double ***S, PBData &pb, GridData &grid);
void getcim2Du(double &uint, double *Du, int *index, int rstar, int sstar, double ***u,
               int gamma[][2], double ***S, PBData &pb, GridData &grid);
void getcim12Du(double &uint, double *Du, int *index, int rstar, int sstar,
                double ***u, double ***S, PBData &pb, GridData &grid);
void checkcim12Du(double ***u, double ***S, PBData &pb, GridData &grid);





double getexactresidual(double ***r, SparseElt2 ****A, double ***b, GridData &grid,
                        double ***S, PBData &pb);
double getexactnormal(int *index, int r, int rstar, int sstar, double alpha, 
                      GridData &grid);
double getexactDnormal(int *index, int r, int s, int rstar, int sstar, double alpha, 
                       GridData &grid);
void checkexact(double ***S, double radius, PBData &pb, GridData &grid);


#endif