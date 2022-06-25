#ifndef FINITEDIFF_H
#define FINITEDIFF_H

#include "grid.h"


char yessk2(int *sk2, int i, int j, int *index, double ***S, GridData &grid);

void getD2(double ***D2, int m, int n, int *sk2, int *tindex, int *N, GridData grid);




#endif