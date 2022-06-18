#ifndef GRID_H
#define GRID_H

struct GridData
{
   static const int dim = 3;
   int nx[dim];
   double a[dim], dx[dim];
   double dt, mindx, t, radius, maxvn, tfinal, radius0;
   double tol;
};


inline int sub2ind(int *index, const int *n, int dim);
inline void ind2sub(int *index, int i, const int *n, int dim);
inline void sub2coord(double *x, int *index, GridData &grid);


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