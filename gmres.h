#ifndef GMRES_H
#define GMRES_H

#include "sparse.h"


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

#endif