#ifndef HYPRESOLVER_H
#define HYPRESOLVER_H

#include "grid.h"
#include "pb.h"
#include "sparse.h"

void HypreSolve(double ***sx, SparseElt2**** &sA, double ***sb, GridData &grid, double ***S, PBData &pb, double ***a);

#endif