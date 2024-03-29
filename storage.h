#ifndef STORAGE_H
#define STORAGE_H

#include "grid.h"
#include "input.h"
#include "sparse.h"


// Dusmall array of StorageStruct
//

struct StorageStruct
{
   int *info;
   SparseElt *head;
   int N;
   int mid;
};


void newstorage(StorageStruct &Dusmall);
int getstoragesize(double ***S, GridData &grid);
int findinstorage(int *info, StorageStruct *Dusmall, int smallsize);
void evalfromstorage(double &uint, double *Du, int *index, int rstar, int sstar,
                     StorageStruct *Dusmall, int smallsize, double ***u, double ***S, 
                     PBData &pb, GridData &grid);

void clearstorage(StorageStruct* &Dusmall, int &smallsize);

void checkDuStorage(double ***u, StorageStruct *Dusmall, int smallsize, double ***S, PBData &pb, GridData &grid);

#endif