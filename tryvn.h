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
#include "amg.h"
#include "march.h"
#include "storage.h"
#include "numerics.h"

// using namespace std;
using std::ofstream;
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ?\
(dmaxarg1) : (dmaxarg2))
static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
(iminarg1) : (iminarg2))








struct ZLHeapElt
{
   int index[3];
   ZLHeapElt* child[2];
   ZLHeapElt* parent;
   int num;
   double val;
};

struct ZLHeapStruct
{
   ZLHeapElt *head;
   ZLHeapElt *tail;
   char ***mark;
   int dim;
};




void output(ofstream& outfile, double ***phi, int *nx);
void output(ofstream& outfile, double **phi, int *nx);
void output(ofstream& outfile, double **phi, char **tube, int *nx);
void output(ofstream& outfile, double **u, double **S, int *nx);
void outputneighbors(double ***S, int *index, GridData &grid);

// void getinit(double*** &S, PBData &pb, MarchStruct &march, TempStruct &tmp, 
//              GridData &grid);


void project(double *w, double *normal, double *v, int dim);

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
void interiorpt2(SparseElt2**** &A, double ***b, int *index, double ***S, PBData &pb, 
                 GridData &grid);
void interiordirection2(SparseElt2**** &A, double ***b, int *index, int r, double ***S, 
                        PBData &pb, GridData &grid);
void interiorptsmall(SparseElt2**** &A, double ***b, int *index, double ***S, 
                     PBData &pb, GridData &grid);
void cim1again2(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
                PBData &pb, GridData &grid);
void cim1(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize,
          int *index, int gamma[][2], double ***S, PBData &pb, GridData &grid);
void cim1(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize,
          int *index, double ***a, int gamma[][2], double ***S, PBData &pb, 
          GridData &grid);
void cim1again3(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
                PBData &pb, GridData &grid);


void perturb(double ***S, double tol, GridData &grid);
void perturb(double ***S, double tol, PBData &pb, GridData &grid);
void perturbstatus(double ***S, double tol, int maxsteps, PBData &pb, GridData &grid);
void perturbelt(double ***S, int *index, double tol);
void perturbelt(double ***S, int *index, double tol, PBData &pb);
char perturbstatuselt(double ***S, int *index, double tol, PBData &pb, GridData &grid);



double getpsivn(double ***u, double ***S, int *index, int rstar, PBData &pb, 
                GridData grid);
double getpsivn2(double ***u, double ***S, int *index, int rstar, PBData &pb, 
                 GridData grid);
double getLJvn(double ***u, double ***S, int *index, int rstar, PBData &pb, 
               GridData grid);
double getvn(double ***psi, double ***LJ, double ***S, int *index, int rstar, 
             PBData &pb, GridData grid);
double getheatvn(double ***u, double ***S, int *index, int rstar, int sstar, 
                 PBData &pb, GridData &grid);
double getheatvn(double ***u, double ***S, int *index, int rstar, int sstar, 
                 StorageStruct *Dusmall, int smallsize, PBData &pb, GridData &grid);
double getrishuvn(double ***u, double ***S, int *index, int rstar, int sstar, 
                  StorageStruct *Dusmall, int smallsize, PBData &pb, GridData &grid);
double getexactvn(int *index, int r, int s, double ***S, GridData &grid);
double formvn(double *Dup, double *Dum, double *normal, GridData &grid);
double LJwater(int *index, PBData &pb, GridData &grid);


void outputstatus(ofstream& outfile, double ***phi, GridData &grid);
void cim2again3(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
                PBData &pb, GridData &grid);
void cim2(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize,
          int *index, int gamma[][2], double ***S, PBData &pb, GridData &grid);
void cim2(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize,
          int *index, double ***aa, int gamma[][2], double ***S, PBData &pb, 
          GridData &grid);
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
void cim345(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize, 
            int *index, double ***a, int gamma[][2], double ***S, PBData &pb, 
            GridData &grid);
double evalcoef(double u0, double ***ucoef, double *uxxcoef, int *index, double ***S, 
                GridData grid);
double evalcoef(double u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                double **jumpuxxcoef, int *index, int rstar, int sstar, double alpha, 
                double ***S, GridData grid);
double evalcoef(double u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                int *index, int rstar, int sstar, double alpha, int mid, 
                double thesign, GridData grid);
double evalcoef(double u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                int *index, int rstar, int sstar, double alpha, int mid, 
                double ***S, GridData grid);
double evalcoef(double u0, double ***ucoef, double *uxcoef, double **uxxcoef, 
                int *index, int rstar, int sstar, double alpha, int mid, 
                double thesign, GridData grid);
void getD1(double ****D1, int *sk, GridData grid);
void getD2(double ***D2[][3], int sk2[][3][4], GridData grid);
void getD2(double ***D2, int m, int n, int sk2[][3][4], GridData grid);
void getD2(double ***D2, int m, int n, int sk2[][3][4], int *tindex, GridData grid);
void getD2(double ***D2, double &jumpuxxcoef, int m, int n, int rstar, int sstar, 
           double thesign, int sk2[][3][4], GridData grid);
void getD2(double ***D2, int m, int n, int *sk2, int *tindex, int *N, GridData grid);
void getcim345D2u(double &u0, double ***u, double *uxcoef, double *uxxcoef, 
                  double &jumpuxxcoef, char &perm, int m, int n, int *index, int rstar, 
                  int sstar, int mid, double ***S, GridData &grid);
void getcim345D2udist(double &u0, double ***u, double *uxcoef, double *uxxcoef, 
                      double &jumpuxxcoef, char &perm, int m, int n, int *index, 
                      int rstar, int sstar, int mid, double ***S, GridData &grid);
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
void getcim345jumpux(double& u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                     double **jumpuxxcoef, int *index, int rstar, int sk, double alpha, 
                     int thesign, double *normal, double *tangent, int mid, 
                     double ****D1ucoef, double **D1uxcoef, double **D1uxxcoef, 
                     double ***D1jumpuxxcoef, double ***S, PBData &pb, GridData &grid);
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
void getcim345jumpuxx(double &u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                      int *index, int rstar, int sk, double alpha, int thesign, 
                      double *normal, int mid, double *D1u, double ****D1ucoef, 
                      double **D1uxcoef, double **D1uxxcoef, double ***D1jumpuxxcoef, 
                      double **D2u, double *****D2ucoef, double ***D2uxcoef, 
                      double ***D2uxxcoef, double **D2jumpuxxcoef, double &jumpD1u,
                      double ***jumpD1ucoef, double *jumpD1uxcoef, 
                      double *jumpD1uxxcoef, double **jumpD1jumpuxxcoef, double ***S, 
                      PBData &pb, GridData &grid);
void getcim345jumpuxx(double &u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                      int *index, int rstar, int sk, double alpha, int thesign, 
                      double *normal, int mid, double ***a, double *D1u, 
                      double ****D1ucoef, double **D1uxcoef, double **D1uxxcoef, 
                      double ***D1jumpuxxcoef, double **D2u, double *****D2ucoef, 
                      double ***D2uxcoef, double ***D2uxxcoef, double **D2jumpuxxcoef, 
                      double &jumpD1u, double ***jumpD1ucoef, double *jumpD1uxcoef, 
                      double *jumpD1uxxcoef, double **jumpD1jumpuxxcoef, double ***S, 
                      PBData &pb, GridData &grid);
void getDu(double ****ucoef, double **uxxcoef, int *index, int rstar,
           int sstar, double alpha, int thesign, int *sk, double ****D1, 
           double*** D2[][3], GridData grid);
void getDu(double **uxcoef, double ***uxxcoef, int *index, int rstar, int sstar, 
           double alpha, int thesign, GridData grid);
void getDu(double *Du, int *index, int rstar, int sstar, double alpha, int thesign, 
           int *sk, double ****D1, double*** D2[][3], GridData grid);
void getcim345Du(double *u0, double ****ucoef, double **uxcoef, double **uxxcoef,
                 double ***jumpuxxcoef, int *index, int rstar, int sstar, double alpha,
                 double thesign, double **D2u, double *****D2ucoef, double ***D2uxcoef, 
                 double ***D2uxxcoef, double **D2jumpuxxcoef, int mid, GridData &grid);
char getsk2(int *sk2, int i, int j, int *index, double ***S, GridData &grid);
char yessk2(int *sk2, int i, int j, int *index, double ***S, GridData &grid);
int getstatus3(double ***S, int *index, GridData grid);
int getstatus4(double ***S, int *index, GridData grid);
double getinterfacegrad3(double *grad, double ***u, double ***S, int *index, int rstar, 
                         int sstar, PBData &pb, GridData grid);
double getinterfacegrad4(double *grad, double ***u, double ***S, int *index, int rstar,
                         int sstar, PBData &pb, GridData grid);
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
char getcim5D2(double ***D2ucoef, double *D2uxcoef, double *D2uxxcoef, int m, int n, 
               int *index, int mid, double ***S, GridData &grid);
char yescim5D2(double ***D2ucoef, double *D2uxcoef, double *D2uxxcoef, int m, int n, 
               int *index, int mid, double ***S, GridData &grid);
int getstatus5(double ***S, int *index, GridData &grid);
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


double getpsivac(double *x, PBData &pb);
void getDpsivac(double *Dpsi, double *x, PBData &pb);
double getDpsivacn(int *index, double ***S, PBData &pb, GridData &grid);
void getD2psivac(double **D2psi, double *x, PBData &pb);
double Bval(double s, PBData &pb);
double Bprime(double s, PBData &pb);
double B2prime(double s, PBData &pb);
double getdotprod(double *v, double *w, int thedim);
void normalize(double *v, double *w, int thedim);
void getcrossprod(double *v, double *w, double *z);
void getunitcrossprod(double *v, double *w, double *z);
void checkanswer(double ***u, double ***S, PBData &pb, GridData &grid);
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
void checkcim345Du(double ***u, StorageStruct *Dusmall, int smallsize, double ***S, 
                   PBData &pb, GridData &grid);
void checkcim345Dupm(double ***u, double ***S, PBData &pb, GridData &grid);
void checkcim345jumpDu(double ***u, double ***S, PBData &pb, GridData &grid);
void getcim4Du(double &uint, double *Du, int *index, int rstar, int sstar, double ***u, 
               double ***S, PBData &pb, GridData &grid);
void getcim5Du(double &uint, double *Du, int *index, int rstar, int sstar, double ***u, 
               double ***S, PBData &pb, GridData &grid);
void getcim345Du(double &uint, double *Du, int *index, int rstar, int sstar,
                 double ***u, double ***S, PBData &pb, GridData &grid);
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
double getexactradius(double thetime, double radius0, double x0, double tol, int Nstep,
                      GridData &grid);
void checkwithexact(double ***S, double radius, GridData &grid);


double getexactresidual(double ***r, SparseElt2 ****A, double ***b, GridData &grid,
                        double ***S, PBData &pb);
double getexactnormal(int *index, int r, int rstar, int sstar, double alpha, 
                      GridData &grid);
double getexactDnormal(int *index, int r, int s, int rstar, int sstar, double alpha, 
                       GridData &grid);
void checkexact(double ***S, double radius, PBData &pb, GridData &grid);
void getinterfaceinfo(double *normal, double *tangent1, double *tangent2, double **Dn,
                      double &tau, double *Dtau, double **D2tau, double &sigma, 
                      double *Dsigma, double &jumpfe, double *x, double ***S,
                      PBData &pb, GridData &grid);
void getiimjumps(double &up0, double &upum, double *uxp0, double **uxpuxm, 
                 double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm,
                 int *index, int rstar, int sk, double alpha, int thesign, 
                 double *normal, double ***S, PBData &pb, GridData &grid);
void getiimjumps(double &up0, double &upum, double *uxp0, double **uxpuxm, 
                 double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm,
                 double *x, int thesign, double ***S, PBData &pb, GridData &grid);
void getiimstencilmx(double **A, int M, int N, double *x, int thesign, int **index, 
                     double up0, double upum, double *uxp0, double **uxpuxm, 
                     double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm, 
                     double ***S, PBData &pb, GridData &grid);
void getiimstencilmx(double **A, int M, int N, double *x, int thesign, int **index, 
                     char yesC, double up0, double upum, double *uxp0, 
                     double **uxpuxm, double **uxxp0, double ***uxxpuxm, 
                     double ****uxxpuxxm, double ***S, PBData &pb, GridData &grid);
void getiimCstencilmx(double **A, int M, int N, double *x, int thesign, int **index, 
                      double up0, double upum, double *uxp0, double **uxpuxm, 
                      double **uxxp0, double ***uxxpuxm, double ****uxxpuxxm, 
                      double ***S, PBData &pb, GridData &grid);
void getiimgridstencilmx(double **A, int M, int N, double *x, int *gindex, int thesign, 
                         int **index, char yesC, double up0, double upum, double *uxp0, 
                         double **uxpuxm, double **uxxp0, double ***uxxpuxm, 
                         double ****uxxpuxxm, double ***S, PBData &pb, GridData &grid);
void iim(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S,
         PBData &pb, GridData &grid);
void iim(SparseElt2**** &A, double ***b, int *index, int add, double ***S, 
         char ***tube, PBData &pb, GridData &grid);
void iim(SparseElt2**** &A, double ***b, int *index, char yesC, int add, double ***S, 
         char ***tube, PBData &pb, GridData &grid);
void iimghost(SparseElt2**** &A, double ***b, int *index, char yesC, int add, 
              double ***S, char ***tube, PBData &pb, GridData &grid);
void iimC(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S,
          PBData &pb, GridData &grid);
void iimC(SparseElt2**** &A, double ***b, int *index, int add, double ***S, 
          char ***tube, PBData &pb, GridData &grid);
void addtoheap(ZLHeapStruct &heap, int *index, double val);
void fixheapeltup(ZLHeapStruct &heap, ZLHeapElt *fix);
ZLHeapElt* fixheapeltempty(ZLHeapStruct &heap, ZLHeapElt *fix);
void fixheapeltdelete(ZLHeapStruct &heap, ZLHeapElt *del);
ZLHeapElt* heapgoto(ZLHeapStruct &heap, int num);
void readZLHeap(ZLHeapStruct heap);
void getnearestgrid(int **index, int &N, double *x, int maxpts, double maxdist, 
                    char ***tube, GridData &grid);
double getdist(double *z, double *x, int thedim);
double bilinearinterp(double *x, double ***u, GridData &grid);
double weno6interpdirect(double *x, double ***u, GridData &grid);
void getallgrad(double ****grad, double ***S, GridData &grid);
int newtondir(double *x, double *x0, double *grad, double tol, double ***S, 
               double ****allgrad, GridData &grid);
int regulafalsidir(double *x, double *x0, double *grad, double tol, double ***S, 
                   GridData &grid);
void getnearestinterface(double *x, int *index, double ***S, GridData &grid) ;

double pythag(double a, double b);
void svdcmp(double **a, int m, int n, double w[], double **v);
double svbksb(double **u, double w[], double **v, int m, int n, double b[], double x[],
              double thresh);

char GSarnoldismall(double**** &Vcol, double** &Hcol, int &k, int maxk, double ***x0, 
                    SparseElt2**** &A, double ***b, double ***S, PBData &pb, 
                    GridData &grid);
char GSarnoldipreleftsmall(double**** &Vcol, double** &Hcol, int &k, int maxk, 
                           double ***x0, SparseElt2**** &A, double ***b, 
                           SparseElt2**** &M, double ***S, PBData &pb, GridData &grid);
char GMRESsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
                int numsteps, double tol, double ***S, PBData &pb, TempStruct &tmp);
char GMRESpreleftsmall(double ***x, SparseElt2**** &A, double ***b, SparseElt2 &M, 
                       GridData &grid, int numsteps, double tol, double ***S, 
                       PBData &pb, TempStruct &tmp);
char GMRESpreleftsmall2(double ***x, SparseElt2**** &A, double ***b, SparseElt2**** &M, 
                        GridData &grid, int numsteps, double tol, double ***S, 
                        PBData &pb, TempStruct &tmp);
void GMRESrestartsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
                       int numsteps, double tol, double ***S, PBData &pb, 
                       TempStruct &tmp);
void GMRESpreleftrestartsmall(double ***x, SparseElt2**** &A, double ***b, 
                              SparseElt2**** &M, GridData &grid, int numsteps, 
                              double tol, double ***S, PBData &pb, TempStruct &tmp);
void testZL(char yesC, int add, double ***S, char ***tube, PBData &pb, 
            GridData &grid);
void testZLatx(double *theerr, double *x, double thesign, char yesC, int add, 
               double ***S, char ***tube, PBData &pb, GridData &grid);



void amgsolve (double ***sx, SparseElt2**** &sA, double ***sb, GridData &grid, double ***S, PBData &pb);

#endif