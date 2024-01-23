#ifndef GRID_H
#define GRID_H
#include <cmath>
#include "global.h"


#include <random>


struct GridData
{
   static const int dim = 3;
   int nx[dim];
   int N;
   double a[dim], dx[dim];
   double dt, mindx, t, radius, maxvn, tfinal, radius0;
   double tol;
};


inline int sub2ind(int *index, const int *n, int dim);
inline void ind2sub(int *index, int i, const int *n, int dim);
inline void sub2coord(double *x, int *index, GridData &grid);

//initialize grid

inline void init_grid(GridData &grid, int nx, double a){
   grid.nx[0] = nx; //number of cell
   grid.nx[1] = grid.nx[0];
   grid.nx[2] = grid.nx[0];
   grid.N = grid.nx[0] * grid.nx[1] * grid.nx[2];
   grid.a[0] = -a;
   grid.a[1] = -a;
   grid.a[2] = -a;
   grid.dx[0] = 2.0*fabs(grid.a[0])/grid.nx[0];
   grid.dx[1] = 2.0*fabs(grid.a[1])/grid.nx[1];
   grid.dx[2] = 2.0*fabs(grid.a[2])/grid.nx[2];
   grid.mindx = fmin(fmin(grid.dx[0],grid.dx[1]),grid.dx[2]);
   grid.dt = grid.mindx*grid.mindx/(2.0*grid.dim);
   grid.tol = 1.0e-14;
   grid.t = 0.0;
   grid.tfinal = globtime;

   grid.radius = RADIUS; // current radius for sphere
   grid.radius0 = RADIUS; // initial radius for sphere




  if(globgridperturb){
    std::default_random_engine generator(globgridperturb);
    std::uniform_real_distribution<double> distribution(0.0, grid.mindx);
    double rand_shift[3] = {distribution(generator), distribution(generator), distribution(generator)};
    grid.a[0] += rand_shift[0];
    grid.a[1] += rand_shift[1];
    grid.a[2] += rand_shift[2];
    printf("randomly perturb grid by (%f,%f,%f)\n",rand_shift[0],rand_shift[1],rand_shift[2]);
  }
}


inline int sub2ind(int *index, const int *n, int dim)
{
   int r, s, i;

   s = 0;
   i = 1;
   for (r = dim-1; r >= 0; r--)
   {
      s += index[r]*i;
      i *= n[r]+1;
   }

   return s;
}

inline void ind2sub(int *index, int i, const int *n, int dim)
{
   int r, j;

   j = i;
   for (r = dim-1; r >= 0; r--)
   {
      index[r] = j%(n[r]+1);
      j = j/(n[r]+1);
   }
}

inline void sub2coord(double *x, int *index, GridData &grid)
{
   int r;

   for (r = 0; r < grid.dim; r++)
      x[r] = grid.a[r]+index[r]*grid.dx[r];
}

#endif