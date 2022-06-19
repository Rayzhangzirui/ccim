#ifndef MARCH_H
#define MARCH_H

#include "grid.h"
#include "input.h"
#include "storage.h"

struct HeapElt
{
   int index[3];
   HeapElt* child[3];
   HeapElt* parent;
   int num;
};

struct HeapStruct
{
   HeapElt *head;
   HeapElt *tail;
   HeapElt ****loc;
   int dim;
};


struct MarchStruct
{
   double ***dist;
   char ***status;
   double ****extend;
   int nval;
   double ***dorig;
   double ****eorig;
   HeapStruct heap;
   double tol;
   PBData pb;

   StorageStruct *Dusmall;
   int smallsize;
};

// initialize marching structure, S is surface
void init_march(MarchStruct &march, double*** S, PBData &pb, GridData&grid);


void fmarch(MarchStruct &march, double phitube, GridData &grid);
void fmarchstep(MarchStruct &march, int *index, GridData &grid);
char fmarchdistance(MarchStruct &march, char *yes, double *xvalue, double **xvalvalue,
                    double *thedx, int *index, GridData &grid);
void fmarchdirection(char *yes, double *value, double **valvalue, double *thedx,
                     MarchStruct &march, int *index, GridData &grid);
void addtoheap(HeapStruct &heap, int *index, double ***phi);
void fixheapeltup(HeapStruct &heap, HeapElt *fix, double ***phi);
HeapElt* fixheapeltempty(HeapStruct &heap, HeapElt *fix, double ***phi);
void fixheapeltdelete(HeapStruct &heap, HeapElt *del, double ***phi);
void fixheapeltreplace(HeapStruct &heap, int *index, double ***phi);
HeapElt* heapgoto(HeapStruct &heap, int num);
HeapElt* evalarray(HeapElt ***A, int *index);
HeapElt* evalarray(HeapElt ****A, int *index);
void setvalarray(HeapElt ***A, int *index, HeapElt *value);
void setvalarray(HeapElt ****A, int *index, HeapElt *value);
HeapElt ***heapmatrix(int row, int col);
HeapElt ****heapmatrix(int row, int col, int frb);
char checkcompat(HeapStruct heap);
void outputheap(HeapStruct heap, double ***value, double ***ovalue);


#endif