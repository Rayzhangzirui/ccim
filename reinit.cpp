#include "reinit.h"
#include "matrix.h"
#include "numerics.h"
#include <iostream>
#include <cmath>

using namespace std;


void reinitfull(double ***u, int numsteps, double ***rhs, double ***der, double ***tmp,
                GridData &grid)
{
   int i, step;
   int tindex[grid.dim];
   double mindx, dt;

   mindx = grid.dx[0];
   for (i = 1; i < grid.dim; i++)
      if (grid.dx[i] < mindx)
         mindx = grid.dx[i];
   dt = mindx/2.0;

   for (step = 1; step <= numsteps; step++)
   {
      cout << "reinit step = " << step << endl;
      getReRHSfull(rhs,u,u,grid);
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(tmp,tindex,evalarray(u,tindex)+dt*evalarray(rhs,tindex));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      getReRHSfull(der,tmp,tmp,grid);
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(rhs,tindex,evalarray(rhs,tindex)+evalarray(der,tindex));
         setvalarray(tmp,tindex,evalarray(u,tindex)+0.25*dt*evalarray(rhs,tindex));

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      getReRHSfull(der,tmp,tmp,grid);
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(rhs,tindex,evalarray(rhs,tindex)+4.0*evalarray(der,tindex));
         setvalarray(u,tindex,evalarray(u,tindex)+dt*evalarray(rhs,tindex)/6.0);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
   }
}

void getReRHSfull(double ***rhs, double ***u, double ***u0, GridData &grid)
{
   double dfp[grid.dim], dfn[grid.dim], df[grid.dim];
   double mindx, eps, sign, grad;
   int deg = 3;
   double u1d[2*deg+1];
   double tol = 1.0e-14;
   int i, r, s;
   int tindex[grid.dim], rindex[grid.dim];

   mindx = grid.dx[0];
   for (i = 1; i < grid.dim; i++)
      if (grid.dx[i] < mindx)
         mindx = grid.dx[i];
   eps = mindx;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
//      u1d[deg] = evalarray(u,tindex);
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         for (s = 0, rindex[r] = max(tindex[r]-deg,0); s <= 2*deg;
              s++, rindex[r] = min(max(tindex[r]-deg+s,0),grid.nx[r]))
            u1d[s] = evalarray(u,rindex);
         weno(dfp[r],dfn[r],u1d,deg,grid.dx[r],2*deg-1);
//         weno(dfp[r],dfn[r],u1d,deg,grid.dx[r],3);
         rindex[r] = tindex[r];
      }
      if (fabs(evalarray(u0,tindex)) < tol)
         setvalarray(rhs,tindex,0.0);
      else
      {
         if (evalarray(u0,tindex) >= tol)
            for (r = 0; r < grid.dim; r++)
            {
               dfp[r] = (dfp[r] < 0.0) ? -dfp[r] : 0.0;
               dfn[r] = (dfn[r] > 0.0 ) ? dfn[r] : 0.0;
               df[r]  = (dfp[r] > dfn[r]) ? dfp[r] : dfn[r];
            }
         else
            for (r = 0; r < grid.dim; r++)
            {
               dfp[r] = (dfp[r] > 0.0) ? dfp[r] : 0.0;
               dfn[r] = (dfn[r] < 0.0) ? -dfn[r] : 0.0;
               df[r]  = (dfp[r] > dfn[r]) ? dfp[r] : dfn[r];
            }
         grad = 0.0;
         for (r = 0; r < grid.dim; r++)
            grad += df[r]*df[r];
         grad = sqrt(grad);
         sign = evalarray(u0,tindex)/sqrt(evalarray(u0,tindex)*
                                          evalarray(u0,tindex)+eps);
         setvalarray(rhs,tindex,sign*(1.0-grad));
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
}
