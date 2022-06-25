#ifndef CIM345COND_H
#define CIM345COND_H


#include "input.h"
#include "storage.h"
#include "sparse.h"
#include "grid.h"
#include <vector>


void checkcim345DuAll(double ***u, double ***S, PBData &pb, GridData &grid);


int getstatus5debug(double ***S, int *index, GridData &grid);

void cim345cond(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize, int *index, double ***a, int gamma[][2], double ***S, PBData &pb, GridData &grid);

char yescim5D2All(std::vector<double***> &D2ucoefvec, std::vector<double*> &D2uxcoefvec, std::vector<double*> &D2uxxcoefvec, std::vector<std::vector<int>> &offsetvec, int m, int n, int *index, int mid, double ***S, GridData &grid);




#endif