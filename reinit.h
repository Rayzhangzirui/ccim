#ifndef REINIT_H
#define REINIT_H

#include "grid.h"

void reinitfull(double ***u, int numsteps, double ***rhs, double ***der, double ***tmp,
                GridData &grid);
void getReRHSfull(double ***rhs, double ***u, double ***u0, GridData &grid);

#endif