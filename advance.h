#ifndef ADVANCE_H
#define ADVANCE_H

#include "grid.h"
#include "sparse.h"
#include "input.h"
#include "march.h"

void advance(double ***u, PBData &pb, MarchStruct &march, TempStruct &tmp, 
             GridData &grid);
void advance(double ***u, double ***a, PBData &pb, MarchStruct &march, TempStruct &tmp, 
             GridData &grid);

void advanceheat(double ***u, TempStruct &tmp, GridData &grid);
void advanceheat(double ***u, double finaltime, TempStruct &tmp, GridData &grid);


#endif