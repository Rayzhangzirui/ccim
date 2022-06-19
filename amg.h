#ifndef AMG_H
#define AMG_H


#include "grid.h"
#include "sparse.h"
#include <iostream>

struct SparseAMGElt
{
   int roc;
   double *val;
   char *strong;
};

struct SparseAMGMx
{
   int *nr;
   SparseAMGElt **row;
   int *nc;
   SparseAMGElt **col;
   int n;
   int m;
};

struct SparseAMGList
{
   SparseAMGMx A;
   SparseAMGMx P;
   SparseAMGList *next;
   SparseAMGList *prev;
   double **B;
   int *PLR, *PLC;
};





struct StrHeapStruct
{
   int *value;
   int *index;
   int *loc;
   int num;
};


void startamg(int *strength, SparseAMGMx &A, double theta);
void startamg(int *strength, int &strhead, int* &strnext, int* &strprev, 
              SparseAMGMx &A, double theta);
void startamg(int *strength, StrHeapStruct &heap, SparseAMGMx &A, double theta);
void amgfirstpass(int *coarse, int *strength, SparseAMGMx &A);
void amgfirstpass(int *coarse, int *strength, int &strhead, int* &strnext, 
                  int* &strprev, SparseAMGMx &A);
void amgfirstpass(int *coarse, int *strength, StrHeapStruct &heap, SparseAMGMx &A);
void amgsecondpass(int *coarse, SparseAMGMx &A);
void makecoarsemx(int *coarse, int *strength, SparseAMGMx &A);
void makecoarsemx(int *coarse, int *strength, int &strhead, int *strnext, int *strprev, 
                  SparseAMGMx &A);
void makecoarsemx(int *coarse, int *strength, StrHeapStruct &heap, SparseAMGMx &A);
char yesinfluencedby(int i, int j, SparseAMGMx &A);

void removeelt(SparseElt* &current, SparseElt* &head);
SparseAMGMx createP(int *coarse, SparseAMGMx &A);
double getdotprod(SparseAMGElt *v, int numv, SparseAMGElt *w, int numw);
double getdotprod(SparseAMGElt *v, int numv, SparseElt *w);
double getdotprod(SparseAMGElt *v, int numv, double *w);
SparseAMGMx multPtP(SparseAMGMx &A, SparseAMGMx &P);
SparseAMGMx multPtP2(SparseAMGMx &A, SparseAMGMx &P);
void clearsparse(SparseAMGMx &A);
void mxvecmult(double *a, SparseAMGMx A, double *b);
void mxtvecmult(double *a, SparseAMGMx A, double *b);
void interpolate(SparseAMGMx &B, double* &c, SparseAMGMx &P, int *coarse, 
                 SparseAMGMx &A, double *b);
void interpolate(SparseAMGMx &B, SparseAMGMx &P, int *coarse, SparseAMGMx &A);
void gaussseidel(double *x, SparseAMGMx &A, double *b, double *x0, int nu);
double Vcycle(double *x, SparseAMGMx &A, double *b, double *x0, int nu1, int nu2, 
              int small);
double Vcycle(double *x, SparseAMGList* &L, double *b, double *x0, int nu1, int nu2, 
              int small);
double getresidual(SparseAMGMx &A, double *x, double *b);
double getresidualinf(SparseAMGMx &A, double *x, double *b);
double getrelresidual(SparseAMGMx &A, double *x, double *b);
double getrelresidualinf(SparseAMGMx &A, double *x, double *b);
void Vcycleloop(double *x, SparseAMGMx &A, double *b, double *x0, int nu1, int nu2, 
                int small, int numsteps, double tol);
void Vcycleloop(double *x, SparseAMGList* &L, double *b, double *x0, int nu1, int nu2, 
                int small, int numsteps, double tol);
void AMGsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, int nu1, 
              int nu2, int small, int numsteps, double tol, double ***S, PBData &pb);
void AMGsmall2(double ***x, SparseElt2**** &A, double ***b, GridData &grid, int nu1, 
               int nu2, int small, int numsteps, double tol, double ***S, PBData &pb);
void AMGsmall3(double ***x, SparseElt2**** &A, double ***b, GridData &grid, int nu1, 
               int nu2, int small, int numsteps, double tol, double ***S, PBData &pb);
void AMGsmall3(double ***x, SparseElt2**** &A, double ***b, GridData &grid, int nu1, 
               int nu2, int small, int numsteps, double tol, double ***a, double ***S, 
               PBData &pb);
void checkmx(SparseAMGMx A);
void peekmx(SparseAMGMx A);
void getpoisson(SparseAMGMx &A, double *b, int *nx, double *dx, int thedim);



char closertohead(StrHeapStruct &heap, int r, int s);
void fixheap(StrHeapStruct &heap, int num);
void fixheapdown(StrHeapStruct &heap, int num);
void fixheapup(StrHeapStruct &heap, int num);
void addtoheap(StrHeapStruct &heap, int index);
void popfromheap(StrHeapStruct &heap, int num);
void switchheapelt(StrHeapStruct &heap, int r, int s);
void sortfromheap(StrHeapStruct &heap);
char checkheap(StrHeapStruct &heap);


void outputAMGmx(std::ofstream& outfile, SparseAMGMx A);
void OutputAmgMxByRow(std::string filename, SparseAMGMx A);
void OutputAmgMxByCol(std::string filename, SparseAMGMx A);

#endif