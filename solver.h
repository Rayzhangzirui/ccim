#ifndef SOLVER_H
#define SOLVER_H


#include "grid.h"
#include "matrix.h"
#include "sparse.h"
#include "input.h"

void gaussseidelsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
                      int numsteps, double ***S, PBData &pb);
void gaussseidelsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
                      int numsteps, double ***a, double ***S, PBData &pb);
void BICGSTAB(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
              int numsteps, double tol, TempStruct &tmp);
void BICGSTABsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
                   int numsteps, double tol, double ***S, PBData &pb, TempStruct &tmp);
void BICGSTABsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
                   int numsteps, double tol, double ***a, double ***S, PBData &pb, 
                   TempStruct &tmp);
void prerightBICGSTABsmall(double ***x, SparseElt2**** &A, double ***b, 
                           SparseElt2**** &M, GridData &grid, int numsteps, double tol, 
                           double ***S, PBData &pb, TempStruct &tmp);
void prerightBICGSTABsmall(double ***x, SparseElt2**** &A, double ***b, 
                           SparseElt2**** &M, GridData &grid, int numsteps, double tol, 
                           double ***a, double ***S, PBData &pb, TempStruct &tmp);
void preBICGSTABsmall(double ***x, SparseElt2**** &A, double ***b, SparseElt2**** M, 
                      GridData &grid, int numsteps, double tol, double ***S, 
                      PBData &pb, TempStruct &tmp);

void ILU(SparseElt2**** &M, SparseElt2**** &A, GridData &grid);
void Atransposesmall(SparseElt2**** &B, SparseElt2**** &A, GridData &grid, double ***S,
                     PBData &pb);
void copyfromsmall(SparseElt2**** &B, SparseElt2**** &A, GridData &grid, double ***S,
                   PBData &pb);
void copyfromsmall(SparseElt2**** &B, SparseElt2**** &A, GridData &grid, double ***a, 
                   double ***S, PBData &pb);
void ILUsmall(SparseElt2**** &M, SparseElt2**** &A, GridData &grid, double ***S, 
              PBData &pb);
void ILUsmall(SparseElt2**** &M, SparseElt2**** &A, GridData &grid, double ***a, 
              double ***S, PBData &pb);


#endif