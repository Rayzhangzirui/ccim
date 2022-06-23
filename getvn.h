#ifndef GETVN_H
#define GETVN_H

#include "storage.h"
#include "grid.h"


double getpsivn(double ***u, double ***S, int *index, int rstar, PBData &pb, 
                GridData grid);
double getpsivn2(double ***u, double ***S, int *index, int rstar, PBData &pb, 
                 GridData grid);
double getLJvn(double ***u, double ***S, int *index, int rstar, PBData &pb, 
               GridData grid);
double getvn(double ***psi, double ***LJ, double ***S, int *index, int rstar, 
             PBData &pb, GridData grid);
double getheatvn(double ***u, double ***S, int *index, int rstar, int sstar, 
                 PBData &pb, GridData &grid);
double getheatvn(double ***u, double ***S, int *index, int rstar, int sstar, 
                 StorageStruct *Dusmall, int smallsize, PBData &pb, GridData &grid);
double getrishuvn(double ***u, double ***S, int *index, int rstar, int sstar, 
                  StorageStruct *Dusmall, int smallsize, PBData &pb, GridData &grid);
double getexactvn(int *index, int r, int s, double ***S, GridData &grid);
double formvn(double *Dup, double *Dum, double *normal, GridData &grid);
double LJwater(int *index, PBData &pb, GridData &grid);

#endif