#ifndef SPARSE_H
#define SPARSE_H

#include "grid.h"
#include "input.h"
#include "matrix.h"
#include <iostream>

struct SparseElt2
{
   double val;
   SparseElt2 *next;
   int *cindex;
};

struct TempStruct
{
   double ****fourd;
   char *fdstatus;
   int Nfd;
   SparseElt2 ****A;
   SparseElt2 ****M;
//   SparseElt2 ****B;
   double ***b;
};




SparseElt2 ****sparseelt2ptrmatrix(int row, int col, int fr);
inline SparseElt2 *evalarray(SparseElt2**** &A, int *index);
inline void setvalarray(SparseElt2**** &A, int *index, SparseElt2 *value);
inline SparseElt2 *evalarray(SparseElt2**** &A, int *index)
{
   return A[index[0]][index[1]][index[2]];
}
inline void setvalarray(SparseElt2**** &A, int *index, SparseElt2 *value)
{
   A[index[0]][index[1]][index[2]] = value;
}


void sparse2(int *rindex, int *cindex, SparseElt2**** &A, double value, GridData &grid);
void sparse2(int *cindex, SparseElt2* &A, double value, GridData &grid);
void sparse2(int *index, SparseElt2* &R, SparseElt2**** &A, GridData &grid);
void sparseorder(int *rindex, int *cindex, SparseElt2**** &A, double value, 
                 GridData &grid);
void sparseorder(int *cindex, SparseElt2* &A, double value, GridData &grid);
void multsparserow2(SparseElt2 ****A, int *index, double value);
void multsparserow2(SparseElt2 *A, double value);
void clearsparse(SparseElt2* &A);
void clearsparse(SparseElt2**** &A, GridData &grid);

void outputsparserow2(SparseElt2 *R, GridData &grid);
void outputsparsesmall(std::ofstream& outfile, SparseElt2 ****A, double ***S, PBData &pb, 
                       GridData &grid);


double ***setfourd(double ****fourd, char *fdstatus, int Nfd, GridData &grid);
void removefourd(double ***r, double ****fourd, char *fdstatus, int Nfd);
void clearfourd(double**** &fourd, char* &fdstatus, int Nfd, GridData &grid);

#endif