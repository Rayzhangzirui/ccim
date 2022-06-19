#include "numerics.h"
#include "matrix.h"
#include <cmath>


// PLU factorization, n,m should be dim-1
int gecp0(double **U, int PLR[], int PLC[], double **A, int n, int m)
{
   int i, j, k, I, J, kl, minn;
   double Amax, temp, Z, eps = 2.2204e-16;
   int tempi;

   for (i = 0; i <= n; i++)
      for (j = 0; j <= m; j++)
         U[i][j] = A[i][j];
   for (j = 0; j <= n; j++)
      PLR[j] = j;
   for (j = 0; j <= m; j++)
      PLC[j] = j;
   Amax = 0.0;
   for (i = 0; i <= n; i++)
      for (j = 0; j <= m; j++)
         if (fabs(A[i][j]) > Amax)
            Amax = fabs(A[i][j]);

   minn = (m < n) ? m : n;
   for (k = 0; k <= minn-1; k++)// k-th step of GE
   {
      Z = 0.0;
      I = k;  J = k;
      for (j = k; j <= m; j++)
         for (i = k; i <= n; i++)
            if (fabs(U[i][j]) > Z)
            {
               Z = fabs(U[i][j]); // max of U(k:n,k:m) at (I,J)
               I = i;
               J = j;
            }
      kl = J;
      for (i = 0; i <= n; i++)//swap column J and kl of U(k:n,k:m)
      {
         temp = U[i][kl];
         U[i][kl] = U[i][k];
         U[i][k] = temp;
      }
      tempi = PLC[kl];//record column permutation in PLC
      PLC[kl] = PLC[k];
      PLC[k] = tempi;
      kl = I;
      for (j = 0; j <= m; j++)//swap row I and row k of U(k:n,k:m)
      {
         temp = U[kl][j];
         U[kl][j] = U[k][j];
         U[k][j] = temp;
      }
      tempi = PLR[kl];//record row permutation in PLR
      PLR[kl] = PLR[k];
      PLR[k] = tempi;

      if (fabs(U[k][k]) <= (eps * Amax))
      {
//         cout << "Input matrix is numerically singular; abort computation" <<
//                 endl;
         return k-1;
      }
      for (i = k+1; i <= n; i++) //Gaussian elimination
      {
         U[i][k] = U[i][k]/U[k][k];
         for (j = k+1; j <= m; j++)
            U[i][j] = U[i][j]-U[i][k]*U[k][j];
      }
   }

   if (n > m)//tall matrix
   {
      Z = 0.0;
      J = m;
      for (i = m; i <= n; i++)
         if (fabs(U[i][m]) > Z)

         {
            Z = fabs(U[i][m]);
            J = i;
         }
      kl = J;
      for (j = 0; j <= m; j++)
      {
         temp = U[kl][j];
         U[kl][j] = U[m][j];
         U[m][j] = temp;
      }
      tempi = PLR[kl];
      PLR[kl] = PLR[m];
      PLR[m] = tempi;

      if (fabs(U[m][m]) <= (eps * Amax))
      {
//         cout << "Input matrix is numerically singular; abort computation" <<
//                 endl;
         return 0;
      }
      for (i = m+1; i <= n; i++)
         U[i][m] = U[i][m]/U[m][m];
   }
   if (m > n) //wide matrix
   {
      Z = 0.0;
      J = n;
      for (j = n; j <= m; j++)
         if (fabs(U[n][j]) > Z)
         {
            Z = fabs(U[n][j]);
            J = j;
         }
      kl = J;
      for (i = 0; i <= n; i++)
      {
         temp = U[i][kl];
         U[i][kl] = U[i][n];
         U[i][n] = temp;
      }
      tempi = PLC[kl];
      PLC[kl] = PLC[n];
      PLC[n] = tempi;
   }
   return 1;
}

void forwardbacksub0(double *x, double *b, double **U, int *PLR, int *PLC, int n)
{
   int j, k;
   double y[n+1];

// back solve L part
   for (j = 0; j <= n; j++)
      y[j] = 0.0;
   for (j = 0; j <= n; j++)
   {
      for (k = 0; k < j; k++)
         y[j] -= U[j][k]*y[k];
      y[j] += b[PLR[j]];
   }

// back solve U part
   for (j = 0; j <= n; j++)
      x[j] = 0.0;
   for (j = n; j >= 0; j--)
   {
      for (k = j+1; k <= n; k++)
         x[PLC[j]] -= U[j][k]*x[PLC[k]];
      x[PLC[j]] += y[j];
      x[PLC[j]] /= U[j][j];
   }
}




double newtonweno(double ***u, int *index, int rstar, int sstar, double alpha, 
                  GridData &grid)
{
   int r, deg = 3, rindex[grid.dim];
   double temp[2*deg];

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < 2*deg; r++)
   {
      rindex[rstar] = fmin(fmax(index[rstar]+sstar*(-deg+1+r),0),grid.nx[rstar]);
      temp[r] = evalarray(u,rindex);
   }
   
   return newtonweno(temp,alpha,1.0-deg,1.0,2*deg-1);
}

double newtonweno(double *u, double x0, double a, double dx, double nx)
{
   int step;
   double x, xold, uval, duval, tol = 1.0e-14;
   char notstop = 1;

   x = x0;
   for (step = 1; notstop; step++)
   {
      dweno6interp1d(uval,duval,x,u,a,dx,nx);
      xold = x;
      x -= uval/duval;
      if (fabs(x-xold) < tol)
         notstop = 0;
   }

   return x;
}

void dweno6interp1d(double &val, double &dval, double x, double *u, double a,
                    double dx, int nx)
{
   int i, s;
   int deg = 3;
   double xi, thesum, C, IS, alpha, p, dp;
   double epsilon = 1.0e-6;

   i = (int) ((x-a)/dx);
   xi = a+i*dx;
  
   val = 0.0;
   dval = 0.0;
   thesum = 0.0;
   s = 0;
   p = u[i-2+s]+(x-(xi-(2.0-s)*dx))*((u[i-1+s]-u[i-2+s])/dx+(x-(xi-(1.0-s)*dx))*
       ((u[i+s]-2.0*u[i-1+s]+u[i-2+s])/(2.0*dx*dx)+(x-(xi+s*dx))*
        (u[i+1+s]-3.0*u[i+s]+3.0*u[i-1+s]-u[i-2+s])/(6.0*dx*dx*dx)));
   dp = (u[i-1+s]-u[i-2+s])/dx+
        (u[i+s]-2.0*u[i-1+s]+u[i-2+s])/(2.0*dx*dx)*
        ((x-(xi-(2.0-s)*dx))+(x-(xi-(1.0-s)*dx)))+
        (u[i+1+s]-3.0*u[i+s]+3.0*u[i-1+s]-u[i-2+s])/(6.0*dx*dx*dx)*
        ((x-(xi-(2.0-s)*dx))*(x-(xi-(1.0-s)*dx))+(x-(xi-(2.0-s)*dx))*(x-(xi+s*dx))+
         (x-(xi-(1.0-s)*dx))*(x-(xi+s*dx)));
   C = (xi+2.0*dx-x)*(xi+3.0*dx-x)/(20.0*dx*dx);
   IS = (-3579.0*u[i+1]*u[i]+2634.0*u[i+1]*u[i-1]-683.0*u[i+1]*u[i-2]-
         6927.0*u[i]*u[i-1]+1854.0*u[i]*u[i-2]-1659.0*u[i-1]*u[i-2]+
         814.0*u[i+1]*u[i+1]+4326.0*u[i]*u[i]+2976.0*u[i-1]*u[i-1]+
         244.0*u[i-2]*u[i-2])/180.0;
   alpha = C/((epsilon+IS)*(epsilon+IS));
   thesum += alpha;
   val += alpha*p;
   dval += alpha*dp;
   s++;
   p = u[i-2+s]+(x-(xi-(2.0-s)*dx))*((u[i-1+s]-u[i-2+s])/dx+(x-(xi-(1.0-s)*dx))*
       ((u[i+s]-2.0*u[i-1+s]+u[i-2+s])/(2.0*dx*dx)+(x-(xi+s*dx))*
        (u[i+1+s]-3.0*u[i+s]+3.0*u[i-1+s]-u[i-2+s])/(6.0*dx*dx*dx)));
   dp = (u[i-1+s]-u[i-2+s])/dx+
        (u[i+s]-2.0*u[i-1+s]+u[i-2+s])/(2.0*dx*dx)*
        ((x-(xi-(2.0-s)*dx))+(x-(xi-(1.0-s)*dx)))+
        (u[i+1+s]-3.0*u[i+s]+3.0*u[i-1+s]-u[i-2+s])/(6.0*dx*dx*dx)*
        ((x-(xi-(2.0-s)*dx))*(x-(xi-(1.0-s)*dx))+(x-(xi-(2.0-s)*dx))*(x-(xi+s*dx))+
         (x-(xi-(1.0-s)*dx))*(x-(xi+s*dx)));
   C = (xi+2.0*dx-x)*(xi+3.0*dx-x)/(10.0*dx*dx);
   IS = (-3777.0*u[i+1]*u[i]+1074.0*u[i+1]*u[i-1]-1269.0*u[i]*u[i-1]+
         1986.0*u[i+1]*u[i+1]+1986.0*u[i]*u[i]+244.0*u[i-1]*u[i-1]+
         244.0*u[i+2]*u[i+2]-1269.0*u[i+2]*u[i+1]+1074.0*u[i+2]*u[i]-
         293.0*u[i+2]*u[i-1])/180.0;
   alpha = C/((epsilon+IS)*(epsilon+IS));
   thesum += alpha;
   val += alpha*p;
   dval += alpha*dp;
   s++;
   p = u[i-2+s]+(x-(xi-(2.0-s)*dx))*((u[i-1+s]-u[i-2+s])/dx+(x-(xi-(1.0-s)*dx))*
       ((u[i+s]-2.0*u[i-1+s]+u[i-2+s])/(2.0*dx*dx)+(x-(xi+s*dx))*
        (u[i+1+s]-3.0*u[i+s]+3.0*u[i-1+s]-u[i-2+s])/(6.0*dx*dx*dx)));
   dp = (u[i-1+s]-u[i-2+s])/dx+
        (u[i+s]-2.0*u[i-1+s]+u[i-2+s])/(2.0*dx*dx)*
        ((x-(xi-(2.0-s)*dx))+(x-(xi-(1.0-s)*dx)))+
        (u[i+1+s]-3.0*u[i+s]+3.0*u[i-1+s]-u[i-2+s])/(6.0*dx*dx*dx)*
        ((x-(xi-(2.0-s)*dx))*(x-(xi-(1.0-s)*dx))+(x-(xi-(2.0-s)*dx))*(x-(xi+s*dx))+
         (x-(xi-(1.0-s)*dx))*(x-(xi+s*dx)));
   C = (x-(xi-2.0*dx))*(x-(xi-dx))/(20.0*dx*dx);
   IS = (-3579.0*u[i+1]*u[i]+4326.0*u[i+1]*u[i+1]+814.0*u[i]*u[i]+
         2976.0*u[i+2]*u[i+2]+244.0*u[i+3]*u[i+3]-683.0*u[i+3]*u[i]-
         6927.0*u[i+2]*u[i+1]+2634.0*u[i+2]*u[i]-1659.0*u[i+3]*u[i+2]+
         1854.0*u[i+3]*u[i+1])/180.0;
   alpha = C/((epsilon+IS)*(epsilon+IS));
   thesum += alpha;
   val += alpha*p;
   dval += alpha*dp;

   val /= thesum;
   dval /= thesum;
}

double weno6interp(double ***u, int *index, int rstar, int sstar, double alpha, 
                   GridData &grid)
{
   int r, deg = 3, rindex[grid.dim];
   double temp[2*deg];

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < 2*deg; r++)
   {
      rindex[rstar] = fmin(fmax(index[rstar]+sstar*(-deg+1+r),0),grid.nx[rstar]);
      temp[r] = evalarray(u,rindex);
   }
   
   return weno6interp1d(alpha,temp,1.0-deg,1.0,2*deg-1);
}

double weno6interp1d(double x, double *u, double a, double dx, int nx)
{
   int i, s;
   int deg = 3;
   double xi, tvalue, thesum, p, C, IS, alpha;
   double epsilon = 1.0e-6;

   i = (int) ((x-a)/dx);
   xi = a+i*dx;
  
   tvalue = 0.0;
   thesum = 0.0;
   s = 0;
   p = u[i-2+s]+(x-(xi-(2.0-s)*dx))*((u[i-1+s]-u[i-2+s])/dx+(x-(xi-(1.0-s)*dx))*
       ((u[i+s]-2.0*u[i-1+s]+u[i-2+s])/(2.0*dx*dx)+(x-(xi+s*dx))*
        (u[i+1+s]-3.0*u[i+s]+3.0*u[i-1+s]-u[i-2+s])/(6.0*dx*dx*dx)));
   C = (xi+2.0*dx-x)*(xi+3.0*dx-x)/(20.0*dx*dx);
   IS = (-3579.0*u[i+1]*u[i]+2634.0*u[i+1]*u[i-1]-683.0*u[i+1]*u[i-2]-
         6927.0*u[i]*u[i-1]+1854.0*u[i]*u[i-2]-1659.0*u[i-1]*u[i-2]+
         814.0*u[i+1]*u[i+1]+4326.0*u[i]*u[i]+2976.0*u[i-1]*u[i-1]+
         244.0*u[i-2]*u[i-2])/180.0;
   alpha = C/((epsilon+IS)*(epsilon+IS));
   thesum += alpha;
   tvalue += alpha*p;
   s++;
   p = u[i-2+s]+(x-(xi-(2.0-s)*dx))*((u[i-1+s]-u[i-2+s])/dx+(x-(xi-(1.0-s)*dx))*
       ((u[i+s]-2.0*u[i-1+s]+u[i-2+s])/(2.0*dx*dx)+(x-(xi+s*dx))*
        (u[i+1+s]-3.0*u[i+s]+3.0*u[i-1+s]-u[i-2+s])/(6.0*dx*dx*dx)));
   C = (xi+2.0*dx-x)*(xi+3.0*dx-x)/(10.0*dx*dx);
   IS = (-3777.0*u[i+1]*u[i]+1074.0*u[i+1]*u[i-1]-1269.0*u[i]*u[i-1]+
         1986.0*u[i+1]*u[i+1]+1986.0*u[i]*u[i]+244.0*u[i-1]*u[i-1]+
         244.0*u[i+2]*u[i+2]-1269.0*u[i+2]*u[i+1]+1074.0*u[i+2]*u[i]-
         293.0*u[i+2]*u[i-1])/180.0;
   alpha = C/((epsilon+IS)*(epsilon+IS));
   thesum += alpha;
   tvalue += alpha*p;
   s++;
   p = u[i-2+s]+(x-(xi-(2.0-s)*dx))*((u[i-1+s]-u[i-2+s])/dx+(x-(xi-(1.0-s)*dx))*
       ((u[i+s]-2.0*u[i-1+s]+u[i-2+s])/(2.0*dx*dx)+(x-(xi+s*dx))*
        (u[i+1+s]-3.0*u[i+s]+3.0*u[i-1+s]-u[i-2+s])/(6.0*dx*dx*dx)));
   C = (x-(xi-2.0*dx))*(x-(xi-dx))/(20.0*dx*dx);
   IS = (-3579.0*u[i+1]*u[i]+4326.0*u[i+1]*u[i+1]+814.0*u[i]*u[i]+
         2976.0*u[i+2]*u[i+2]+244.0*u[i+3]*u[i+3]-683.0*u[i+3]*u[i]-
         6927.0*u[i+2]*u[i+1]+2634.0*u[i+2]*u[i]-1659.0*u[i+3]*u[i+2]+
         1854.0*u[i+3]*u[i+1])/180.0;
   alpha = C/((epsilon+IS)*(epsilon+IS));
   thesum += alpha;
   tvalue += alpha*p;

   tvalue /= thesum;

   return tvalue;
}

double regulafalsiweno(double ***u, int *index, int rstar, int sstar, GridData &grid)
{
   int r, deg = 3, rindex[grid.dim];
   double temp[2*deg];

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < 2*deg; r++)
   {
      rindex[rstar] = fmin(fmax(index[rstar]+sstar*(-deg+1+r),0),grid.nx[rstar]);//-2, -1, 0, 1, 2, 3
      temp[r] = evalarray(u,rindex);
   }
   
   return regulafalsiweno(temp,deg-1,deg,1.0-deg,1.0,2*deg-1);
}

double regulafalsiweno(double *u, int i, int j, double a, double dx, double nx)
{
   int step;
   double b, c, s, ub, uc, us, tol = 1.0e-14;
   char notstop = 1;

   b = a+i*dx;
   c = a+j*dx;
   ub = u[i];
   uc = u[j];
   for (step = 1; notstop; step++)
   {
      s = b-ub*(b-c)/(ub-uc);
      us = weno6interp1d(s,u,a,dx,nx);
      if (fmin(fabs(s-b),fabs(s-c)) < tol || us == 0.0)
         notstop = 0;
      else
         if (us*ub < 0.0)
         {
            c = s;
            uc = us;
         }
         else
         {
            b = s;
            ub = us;
         }
   }

   return s;
}

double lagrangeinterp1d(double x, double *u, double a, double dx, int nx)
{
   int i, j;
   double tvalue, theprod;

   tvalue = 0.0;
   for (i = 0; i <= nx; i++)
   {
      theprod = 1.0;
      for (j = 0; j <= nx; j++)
         if (i != j)
            theprod *= (x-(a+j*dx))/((i-j)*dx);
      tvalue += u[i]*theprod;
   }

   return tvalue;
}

double regulafalsi(double ***u, int *index, int rstar, int sstar, GridData &grid)
{
   int r, deg = 3, rindex[grid.dim];
   double temp[2*deg];

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];

//   rindex[rstar] = index[rstar]+sstar;
//   return -evalarray(u,index)/(evalarray(u,rindex)-evalarray(u,index));

   for (r = 0; r < 2*deg; r++)
   {
      rindex[rstar] = fmin(fmax(index[rstar]+sstar*(-deg+1+r),0),grid.nx[rstar]);
      temp[r] = evalarray(u,rindex);
   }
   
   return regulafalsi(temp,deg-1,deg,1.0-deg,1.0,2*deg-1);
/*
//   return regulafalsi(temp,deg-1,deg,(1.0-deg)*grid.dx[rstar],grid.dx[rstar],2*deg-1)/
//          grid.dx[rstar];
   if (index[0] == 20 && index[1] == 26 && index[2] == 34)
      globdebug2 = 1;
   double val = regulafalsi(temp,deg-1,deg,grid.a[rstar]+(index[rstar]+
                                           (1.0-deg))*grid.dx[rstar],grid.dx[rstar],
                            2*deg-1);
   double x = -1.0+2.0/70.0*20.0, y = -1.0+2.0/70.0*26, z;
   if (index[0] == 20 && index[1] == 26 && index[2] == 34)
   {
      cout << "IN regula falsi " << val << " " << sqrt(x*x+y*y+val*val)-0.5 << endl;
      cout << "   " << rstar << " " << sstar << endl;
//      for (r = 0; r < 2*deg; r++)
//      {
//         z = -1.0+(index[2]+sstar*(-deg+1+r))*2.0/70.0;
//         cout << "   " << temp[r]-(sqrt(x*x+y*y+z*z)-0.5) << endl;
//      }
      double *temp2;
      temp2 = new double[71];
      for (r = 0; r <= 70; r++)
      {
         z = -1.0+2.0/70.0*r;
         temp2[r] = sqrt(x*x+y*y+z*z)-0.5;
      }
      for (r = 0; r < 2*deg; r++)
      {
         cout << "   " << temp[r]-temp2[34-deg+1+r] << endl;
         cout << "   " << -1+2.0/70.0*(34-deg+1+r)-
                          grid.a[rstar]-(index[rstar]+1.0-deg+r)*grid.dx[rstar] << endl;
      }
      int p = 34;
      cout << "   " << "SECOND TIME" << endl;
      double val2 = regulafalsi(temp2,p,p+1,-1.0,2.0/70.0,70);
      cout << "   " << val2 << endl;
      globdebug2 = 0;

      delete [] temp2;
   }

   double xx[grid.dim];
   sub2coord(xx,index,grid);
   return (val-xx[rstar])/grid.dx[rstar];
*/
}

double regulafalsi(double *u, int i, int j, double a, double dx, double nx)
{
   int step;
   double b, c, s, ub, uc, us, tol = 1.0e-15;
   char notstop = 1;

   b = a+i*dx;
   c = a+j*dx;
   ub = u[i];
   uc = u[j];
   for (step = 1; notstop; step++)
   {
      s = b-ub*(b-c)/(ub-uc);
      us = lagrangeinterp1d(s,u,a,dx,nx);
      if (fmin(fabs(s-b),fabs(s-c)) < tol || us == 0.0)
         notstop = 0;
      else
         if (us*ub < 0.0)
         {
            c = s;
            uc = us;
         }
         else
         {
            b = s;
            ub = us;
         }
   }

   return s;
}

double regulafalsiexact(int *index, int rstar, int sstar, GridData &grid)
{
   int r, step;
   double b, c, s, ub, uc, us, tol = 1.0e-15, x[grid.dim];
   char notstop = 1;

   sub2coord(x,index,grid);
   if (sstar > 0)
   {
      b = x[rstar];
      c = x[rstar]+sstar*grid.dx[rstar];
   }
   else
   {
      b = x[rstar]+sstar*grid.dx[rstar];
      c = x[rstar];
   }
   ub = 0.0;
   uc = 0.0;
   for (r = 0; r < grid.dim; r++)
      if (r != rstar)
      {
         ub += x[r]*x[r]; 
         uc += x[r]*x[r]; 
      }
      else
      {
         ub += b*b; 
         uc += c*c; 
      }
   ub = sqrt(ub)-0.5;
   uc = sqrt(uc)-0.5;
   
   for (step = 1; notstop; step++)
   {
      s = b-ub*(b-c)/(ub-uc);
      us = 0.0;
      for (r = 0; r < grid.dim; r++)
         if (r != rstar)
            us += x[r]*x[r]; 
         else
            us += s*s; 
      us = sqrt(us)-0.5;
      if (fmin(fabs(s-b),fabs(s-c)) < tol || us == 0.0)
         notstop = 0;
      else
         if (us*ub < 0.0)
         {
            c = s;
            uc = us;
         }
         else
         {
            b = s;
            ub = us;
         }
   }

   return fabs(s-x[rstar])/grid.dx[rstar];
}

double bisection(double ***u, int *index, int rstar, int sstar, GridData &grid)
{
   int r, deg = 3, rindex[grid.dim];
   double temp[2*deg];

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < 2*deg; r++)
   {
      rindex[rstar] = fmin(fmax(index[rstar]+sstar*(-deg+1+r),0),grid.nx[rstar]);
      temp[r] = evalarray(u,rindex);
   }
   
   return bisection(temp,deg-1,deg,1.0-deg,1.0,2*deg-1);
//   return bisection(temp,deg-1,deg,(1.0-deg)*grid.dx[rstar],grid.dx[rstar],2*deg-1)/
//          grid.dx[rstar];
//   double val = regulafalsi(temp,deg-1,deg,grid.a[rstar]+(index[rstar]+
//                                           (1.0-deg))*grid.dx[rstar],grid.dx[rstar],
//                            2*deg-1);
//   double xx[grid.dim];
//   sub2coord(xx,index,grid);
//   return (val-xx[rstar])/grid.dx[rstar];
}

double bisection(double *u, int i, int j, double a, double dx, double nx)
{
   int step;
   double b, c, s, ub, uc, us, tol = 1.0e-15;
   char notstop = 1;

   b = a+i*dx;
   c = a+j*dx;
   ub = u[i];
   uc = u[j];
   for (step = 1; notstop; step++)
   {
      s = b+(c-b)/2.0;
      us = lagrangeinterp1d(s,u,a,dx,nx);
      if (fmin(fabs(s-b),fabs(s-c)) < tol || us == 0.0)
         notstop = 0;
      else
         if ((us < 0.0)+(ub < 0.0) == 1)
         {
            c = s;
            uc = us;
         }
         else
         {
            b = s;
            ub = us;
         }
   }

   return s;
}



// fi is function value, dfp is forward diff, dfn is backward diff
void weno(double& dfp, double& dfn, double* fi, const long i,
          const double dx, const int order)
{
  double aa, bb, cc, dd, ee, ff, df, ih;
  double IS0, IS1, IS2, a0, a1, a2, w0, w1, w2;
  const double eps = 1.e-6;

  ih  = 1. / dx;

   if (order == 1)
   {
      dfp = (fi[i+1]-fi[i])/dx;
      dfn = (fi[i]-fi[i-1])/dx;
   }
  else if ( order == 3 ) {                // WENO-3
    aa = fi[i-1] - fi[i-2];
    bb = fi[i]   - fi[i-1];
    cc = fi[i+1] - fi[i];
    dd = fi[i+2] - fi[i+1];

    df = bb + cc;

    a0 = bb - aa;
    a1 = cc - bb;
    a2 = dd - cc;

    IS0 = ( eps + a2*a2 ) / (eps + a1*a1);
    w0  = 1. / (1. + 2.*IS0*IS0);

    dfp = .5 * ih * ( df - w0*(dd - 2.*cc + bb) );

    IS0 = ( eps + a0*a0 ) / (eps + a1*a1);
    w0  = 1 / (1 + 2.*IS0*IS0);

    dfn = .5 * ih * ( df - w0*(aa - 2.*bb + cc) );
  }

  else {                             // WENO-5
    aa = (fi[i-2] - fi[i-3]) * ih;
    bb = (fi[i-1] - fi[i-2]) * ih;
    cc = (fi[i]   - fi[i-1]) * ih;
    dd = (fi[i+1] - fi[i]  ) * ih;
    ee = (fi[i+2] - fi[i+1]) * ih;
    ff = (fi[i+3] - fi[i+2]) * ih;

    df = (-bb + 7.* (cc + dd) - ee) / 12.;

    aa = bb - aa;
    bb = cc - bb;
    cc = dd - cc;
    dd = ee - dd;
    ee = ff - ee;

    //  dfp = df + fweno(ee, dd, cc, bb);
    //  dfn = df - fweno(aa, bb, cc, dd);

    IS0 = 13. * (ee-dd) * (ee-dd) + 3. * (ee-3.*dd) * (ee-3.*dd);
    IS1 = 13. * (dd-cc) * (dd-cc) + 3. * (dd +  cc) * (dd +  cc);
    IS2 = 13. * (cc-bb) * (cc-bb) + 3. * (3.*cc-bb) * (3.*cc-bb);

    a0 = 1. / ( (eps + IS0) * (eps + IS0) );
    a1 = 6. / ( (eps + IS1) * (eps + IS1) );
    a2 = 3. / ( (eps + IS2) * (eps + IS2) );

    w1 = 1. / ( a0 + a1 + a2 );
    w0 = a0 * w1;
    w2 = a2 * w1;

    dfp = df + w0 * ( ee - 2.*dd + cc ) / 3.  +
                    ( w2 - .5 ) * ( dd - 2.*cc + bb ) / 6.;

    IS0 = 13. * (aa-bb) * (aa-bb) + 3. * (aa-3.*bb) * (aa-3.*bb);
    IS1 = 13. * (bb-cc) * (bb-cc) + 3. * (bb +  cc) * (bb +  cc);
    IS2 = 13. * (cc-dd) * (cc-dd) + 3. * (3.*cc-dd) * (3.*cc-dd);

    a0 = 1. / ( (eps + IS0) * (eps + IS0) );
    a1 = 6. / ( (eps + IS1) * (eps + IS1) );
    a2 = 3. / ( (eps + IS2) * (eps + IS2) );

    w1 = 1. / ( a0 + a1 + a2 );
    w0 = a0 * w1;
    w2 = a2 * w1;

    dfn = df - w0 * ( aa - 2.*bb + cc ) / 3.  -
                    ( w2 - .5 ) * ( bb - 2.*cc + dd ) / 6.;
  }
}