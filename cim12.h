#ifndef CIM12_H
#define CIM12_H

#include "sparse.h"
#include "storage.h"

void cim1(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize, int *index, 
        int gamma[][2], double ***S, PBData &pb, GridData &grid);
void cim1(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize, int *index, double ***a,
         int gamma[][2], double ***S, PBData &pb,GridData &grid);

void cim2(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize, int *index,
         int gamma[][2], double ***S, PBData &pb, GridData &grid);
void cim2(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize,int *index, double ***aa, 
        int gamma[][2], double ***S, PBData &pb, GridData &grid);

int getstatus3(double ***S, int *index, GridData grid);

void getcim1Du(double &uint, double *Du, int *index, int rstar, int sstar, double ***u,
               int gamma[][2], double ***S, PBData &pb, GridData &grid);
void getcim2Du(double &uint, double *Du, int *index, int rstar, int sstar, double ***u,
               int gamma[][2], double ***S, PBData &pb, GridData &grid);
void getcim12Du(double &uint, double *Du, int *index, int rstar, int sstar,
                double ***u, double ***S, PBData &pb, GridData &grid);
void checkcim12Du(double ***u, double ***S, PBData &pb, GridData &grid);


double getinterfacegrad4(double *grad, double ***u, double ***S, int *index, int rstar,int sstar, PBData &pb, GridData grid);
#endif