#include "gmres.h"
#include "numerics.h"
#include <iostream>
#include <cmath>
using namespace std;

char GSarnoldismall(double**** &Vcol, double** &Hcol, int &k, int maxk, double ***x0, 
                    SparseElt2**** &A, double ***b, double ***S, PBData &pb, 
                    GridData &grid)
// if k < maxk, creates Vcol[k] and, if k > 0, Hcol[k-1][0:k].
// note: creates memory for Vcol[k], though not correct one, even if 
//    Hcol[k-1][k] <= grid.tol
// returns 0 if Hcol[k-1][k] <= grid.tol or if k >= maxk
{
// create v_{k+1}
   int i, j, m, n, tindex[grid.dim], rindex[grid.dim];
   double temp, ehere;
   SparseElt2* current;

//   cout << "in GS arnoldi" << endl;
   k++;
   Vcol[k] = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   if (k > 0 && k < maxk)
   {
// create H column based on v_1,...,v_k
//      cout << "in main part" << endl;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         for (m = 0; m < grid.dim; m++)
            rindex[m] = tindex[m];
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;

         temp = 0.0;
         if (evalarray(A,tindex) == NULL)
            for (m = 0; m < grid.dim; m++)
            {
               temp += 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(Vcol[k-1],rindex);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp += -ehere/(grid.dx[m]*grid.dx[m])*evalarray(Vcol[k-1],rindex);
               }
               rindex[m] = tindex[m];
            }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               temp += (*current).val*evalarray(Vcol[k-1],(*current).cindex);
         setvalarray(Vcol[k],tindex,temp);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
//      cout << "finished Av_k" << endl;

      Hcol[k-1] = new double[k+1];
      for (j = 0; j < k; j++)
      {
         Hcol[k-1][j] = 0.0;
         for (i = 0; i < grid.dim; i++)
            tindex[i] = 0;
         while (tindex[0] <= grid.nx[0])
         {
            Hcol[k-1][j] += evalarray(Vcol[j],tindex)*evalarray(Vcol[k],tindex);

            (tindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
            {
               tindex[i] = 0;
               (tindex[i-1])++;
            }
         }
      }
//      cout << "finished most of H" << endl;

// create v_{k+1}
      for (j = 0; j < k; j++)
      {
         for (i = 0; i < grid.dim; i++)
            tindex[i] = 0;
         while (tindex[0] <= grid.nx[0])
         {
            setvalarray(Vcol[k],tindex,evalarray(Vcol[k],tindex)-
                                       Hcol[k-1][j]*evalarray(Vcol[j],tindex));

            (tindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
            {
               tindex[i] = 0;
               (tindex[i-1])++;
            }
         }
      }

// create rest of hat H based on v_{k+1}
      Hcol[k-1][k] = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         Hcol[k-1][k] += evalarray(Vcol[k],tindex)*evalarray(Vcol[k],tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      Hcol[k-1][k] = sqrt(Hcol[k-1][k]);
      cout << "last H value = " << Hcol[k-1][k] << endl;

      if (Hcol[k-1][k] > grid.tol)
      {
         for (i = 0; i < grid.dim; i++)
            tindex[i] = 0;
         while (tindex[0] <= grid.nx[0])
         {
            setvalarray(Vcol[k],tindex,evalarray(Vcol[k],tindex)/Hcol[k-1][k]);

            (tindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
            {
               tindex[i] = 0;
               (tindex[i-1])++;
            }
         }

/*
         double theerr = 0.0;
//         cout << "should be hat H = " << endl;
         for (int ii = 0; ii < k+1; ii++)
         {
            for (int jj = 0; jj < k; jj++)
            {
               double ***a = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   
               for (i = 0; i < grid.dim; i++)
                  tindex[i] = 0;
               while (tindex[0] <= grid.nx[0])
               {
                  for (m = 0; m < grid.dim; m++)
                     rindex[m] = tindex[m];
                  if (evalarray(S,tindex) < 0.0)
                     ehere = pb.epsilonm;
                  else
                     ehere = pb.epsilonp;
   
                  double temp = 0.0;
                  if (evalarray(A,tindex) == NULL)
                     for (m = 0; m < grid.dim; m++)
                     {
                        temp += 2.0*ehere/(grid.dx[m]*grid.dx[m])*
                                evalarray(Vcol[jj],rindex);
                        for (n = -1; n <= 1; n += 2)
                        {
                           rindex[m] = tindex[m]+n;
                           if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                              temp += -ehere/(grid.dx[m]*grid.dx[m])*
                                      evalarray(Vcol[jj],rindex);
                        }
                        rindex[m] = tindex[m];
                     }
                  else
                     for (current = evalarray(A,tindex); current != NULL; 
                          current = (*current).next)
                        temp += (*current).val*evalarray(Vcol[jj],(*current).cindex);
   
                  setvalarray(a,tindex,temp);
   
                  (tindex[grid.dim-1])++;
                  for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
                  {
                     tindex[i] = 0;
                     (tindex[i-1])++;
                  }
               }
   
               double val = 0.0;
               for (i = 0; i < grid.dim; i++)
                  tindex[i] = 0;
               while (tindex[0] <= grid.nx[0])
               {
                  val += evalarray(Vcol[ii],tindex)*evalarray(a,tindex);
   
                  (tindex[grid.dim-1])++;
                  for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
                  {
                     tindex[i] = 0;
                     (tindex[i-1])++;
                  }
               }
               if (ii <= jj+1 && fabs(val-Hcol[jj][ii]) > theerr)
                  theerr = fabs(val-Hcol[jj][ii]);
               else if (ii > jj+1 && fabs(val) > theerr)
                  theerr = fabs(val);

//               cout << val << " ";
               free_matrix(a,grid.nx[0],grid.nx[1],grid.nx[2]);
            }
//            cout << endl;
         }
         cout << "   THEERR = " << theerr << endl;
*/

/*
         cout << "is hat H = " << endl;
         for (int ii = 0; ii < k+1; ii++)
         {
            for (int jj = 0; jj < k; jj++)
               if (ii <= jj+1)
                  cout << Hcol[jj][ii] << " ";
               else
                  cout << "0.0 ";
            cout << endl;
         }
*/
      }
      else
         return 0;
   }
   else if (k == 0)
   {
// create v_1 based on initial guess
      double r0norm = 0.0;

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         for (m = 0; m < grid.dim; m++)
            rindex[m] = tindex[m];
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;

         temp = evalarray(b,tindex);
         if (evalarray(A,tindex) == NULL)
            for (m = 0; m < grid.dim; m++)
            {
               temp -= 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(x0,rindex);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp -= -ehere/(grid.dx[m]*grid.dx[m])*evalarray(x0,rindex);
               }
               rindex[m] = tindex[m];
            }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               temp -= (*current).val*evalarray(x0,(*current).cindex);

         setvalarray(Vcol[k],tindex,temp);
         r0norm += temp*temp;

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      r0norm = sqrt(r0norm);

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(Vcol[k],tindex,evalarray(Vcol[k],tindex)/r0norm);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
   }
   else
      return 0;
   
/*
   double theerr = 0.0;
//   cout << "check newest = " << endl;
   for (m = 0; m < k; m++)
   {
      for (n = 0; n < k; n++)
      {
         double r0norm = 0.0;
         for (i = 0; i < grid.dim; i++)
            tindex[i] = 0;
         while (tindex[0] <= grid.nx[0])
         {
            r0norm += evalarray(Vcol[m],tindex)*evalarray(Vcol[n],tindex);

            (tindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
            {
               tindex[i] = 0;
               (tindex[i-1])++;
            }
         }
         if (m == n && fabs(r0norm-1.0) > theerr)
            theerr = fabs(r0norm-1.0);
         else if (m != n && fabs(r0norm) > theerr)
            theerr = fabs(r0norm);
//         cout << r0norm << " ";
      }   
//      cout << endl;
   }
   cout << "   THEERR2 = " << theerr << endl;
*/
   return 1;
}
 
char GSarnoldipreleftsmall(double**** &Vcol, double** &Hcol, int &k, int maxk, 
                           double ***x0, SparseElt2**** &A, double ***b, 
                           SparseElt2**** &M, double ***S, PBData &pb, GridData &grid)
// if k < maxk, creates Vcol[k] and, if k > 0, Hcol[k-1][0:k].
// note: creates memory for Vcol[k], though not correct one, even if 
//    Hcol[k-1][k] <= grid.tol
// returns 0 if Hcol[k-1][k] <= grid.tol or if k >= maxk
{
// create v_{k+1}
   int i, j, m, n, tindex[grid.dim], rindex[grid.dim];
   double temp, ehere;
   SparseElt2* current;

//   cout << "in GS arnoldi" << endl;
   k++;
   Vcol[k] = matrix(grid.nx[0],grid.nx[1],grid.nx[2]);
   if (k > 0 && k < maxk)
   {
// create H column based on v_1,...,v_k
//      cout << "in main part" << endl;
      leftmultmxsmall(Vcol[k],A,Vcol[k-1],grid,S,pb);
      leftmultILUinv(Vcol[k],M,Vcol[k],grid);
//      cout << "finished Av_k" << endl;

      Hcol[k-1] = new double[k+1];
      for (j = 0; j < k; j++)
      {
         Hcol[k-1][j] = 0.0;
         for (i = 0; i < grid.dim; i++)
            tindex[i] = 0;
         while (tindex[0] <= grid.nx[0])
         {
            Hcol[k-1][j] += evalarray(Vcol[j],tindex)*evalarray(Vcol[k],tindex);

            (tindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
            {
               tindex[i] = 0;
               (tindex[i-1])++;
            }
         }
      }
//      cout << "finished most of H" << endl;

// create v_{k+1}
      for (j = 0; j < k; j++)
      {
         for (i = 0; i < grid.dim; i++)
            tindex[i] = 0;
         while (tindex[0] <= grid.nx[0])
         {
            setvalarray(Vcol[k],tindex,evalarray(Vcol[k],tindex)-
                                       Hcol[k-1][j]*evalarray(Vcol[j],tindex));

            (tindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
            {
               tindex[i] = 0;
               (tindex[i-1])++;
            }
         }
      }

// create rest of hat H based on v_{k+1}
      Hcol[k-1][k] = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         Hcol[k-1][k] += evalarray(Vcol[k],tindex)*evalarray(Vcol[k],tindex);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      Hcol[k-1][k] = sqrt(Hcol[k-1][k]);
      cout << "last H value = " << Hcol[k-1][k] << endl;

      if (Hcol[k-1][k] > grid.tol)
      {
         for (i = 0; i < grid.dim; i++)
            tindex[i] = 0;
         while (tindex[0] <= grid.nx[0])
         {
            setvalarray(Vcol[k],tindex,evalarray(Vcol[k],tindex)/Hcol[k-1][k]);

            (tindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
            {
               tindex[i] = 0;
               (tindex[i-1])++;
            }
         }
      }
      else
         return 0;
   }
   else if (k == 0)
   {
// create v_1 based on initial guess
      double r0norm = 0.0;

      leftmultmxsmall(Vcol[k],A,x0,grid,S,pb);
      leftmultILUinv(Vcol[k],M,Vcol[k],grid);
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         temp = evalarray(b,tindex)-evalarray(Vcol[k],tindex);
         setvalarray(Vcol[k],tindex,temp);
         r0norm += temp*temp;

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      r0norm = sqrt(r0norm);

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(Vcol[k],tindex,evalarray(Vcol[k],tindex)/r0norm);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
   }
   else
      return 0;
   
   return 1;
}
 
char GMRESsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
                int numsteps, double tol, double ***S, PBData &pb, TempStruct &tmp)
{
   int i, j, k, s, m, n;
   int tindex[grid.dim], rindex[grid.dim];
   double residual, r0norm, bnorm, temp, ehere; 
   double ***z = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   SparseElt2 *current;

// if maxk == numsteps, then last k is numsteps-1, which creates Vcol[numsteps-1]
//    and Hcol[numsteps-2][0...numsteps-1]
   double **Hcol = new double*[numsteps-1];
   double ****Vcol = new double***[numsteps];

// get two norm of b
   bnorm = 0.0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      bnorm += evalarray(b,tindex)*evalarray(b,tindex);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   bnorm = sqrt(bnorm);

// get residual and its two norm
   r0norm = 0.0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (m = 0; m < grid.dim; m++)
         rindex[m] = tindex[m];
      if (evalarray(S,tindex) < 0.0)
         ehere = pb.epsilonm;
      else
         ehere = pb.epsilonp;

      temp = evalarray(b,tindex);
      if (evalarray(A,tindex) == NULL)
         for (m = 0; m < grid.dim; m++)
         {
            temp -= 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
            for (n = -1; n <= 1; n += 2)
            {
               rindex[m] = tindex[m]+n;
               if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                  temp -= -ehere/(grid.dx[m]*grid.dx[m])*evalarray(x,rindex);
            }
            rindex[m] = tindex[m];
         }
      else
         for (current = evalarray(A,tindex); current != NULL; 
              current = (*current).next)
            temp -= (*current).val*evalarray(x,(*current).cindex);

      r0norm += temp*temp;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   r0norm = sqrt(r0norm);

// get Vcol[0]
   k = -1;
   GSarnoldismall(Vcol,Hcol,k,numsteps,x,A,b,S,pb,grid);
   cout << "Set initial column of V" << endl;
   residual = r0norm/bnorm;
   cout << "   residual = " << residual << endl;
// if maxk == numsteps, then last k is numsteps-1, which creates Vcol[numsteps-1]
//    and Hcol[numsteps-2][0...numsteps-1]
   double **B = matrix(numsteps-2,numsteps-2);
   double *c = new double[numsteps-1];
   double *y = new double[numsteps-1];
   int *PLR = new int[numsteps-1];
   int *PLC = new int[numsteps-1];
   double **LU = matrix(numsteps-2,numsteps-2);

   char notstop = 1;
   for ( ; k+1 < numsteps && residual >= tol && notstop; )
   {
      notstop = GSarnoldismall(Vcol,Hcol,k,numsteps,x,A,b,S,pb,grid);
      cout << "Finished " << k << "th column of V" << endl;
//      cout << "B = " << endl;
      for (i = 0; i < k; i++) 
      {
         for (j = 0; j < k; j++) 
         {
            B[i][j] = 0.0;
            for (s = 0; s <= min(i+1,j+1); s++) 
               B[i][j] += Hcol[i][s]*Hcol[j][s];
//            cout << B[i][j] << " ";
         }
//         cout << endl;
         c[i] = Hcol[i][0]*r0norm;
      }
      gecp0(LU,PLR,PLC,B,k-1,k-1);
      forwardbacksub0(y,c,LU,PLR,PLC,k-1);
/*
      double **BB = matrix(k,k-1);
      double *cc = new double[k+1];
      double *w = new double[k];
      double **CC = matrix(k-1,k-1);
      for (i = 0; i < k+1; i++) 
      {
         for (j = 0; j < k; j++) 
            if (i <= j+1)
               BB[i][j] = Hcol[j][i];
            else
               BB[i][j] = 0.0;
         cc[i] = 0.0;
      }
      cc[0] = r0norm;
      svdcmp(BB,k+1,k,w,CC);
      svbksb(BB,w,CC,k+1,k,cc,y,1.0e-14);
      free_matrix(BB,k,k-1);
      free_matrix(CC,k-1,k-1);
      delete [] cc;
      delete [] w;
*/

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(z,tindex,evalarray(x,tindex));
         for (s = 0; s < k; s++)
            setvalarray(z,tindex,evalarray(z,tindex)+evalarray(Vcol[s],tindex)*y[s]);
       
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      residual = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         for (m = 0; m < grid.dim; m++)
            rindex[m] = tindex[m];
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;
   
         temp = evalarray(b,tindex);
         if (evalarray(A,tindex) == NULL)
            for (m = 0; m < grid.dim; m++)
            {
               temp -= 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(z,rindex);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp -= -ehere/(grid.dx[m]*grid.dx[m])*evalarray(z,rindex);
               }
               rindex[m] = tindex[m];
            }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               temp -= (*current).val*evalarray(z,(*current).cindex);
   
         residual += temp*temp;
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      residual = sqrt(residual)/bnorm;
      cout << "   residual = " << residual << endl;
//      getchar();
   }

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(x,tindex,evalarray(z,tindex));

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   free_matrix(B,numsteps-2,numsteps-2);
   free_matrix(LU,numsteps-2,numsteps-2);
   delete [] c;
   delete [] y;
   delete [] PLR;
   delete [] PLC;

   for (i = 0; i <= min(k,numsteps-1); i++)
      free_matrix(Vcol[i],grid.nx[0],grid.nx[1],grid.nx[2]);
   delete [] Vcol;
   for (i = 0; i < min(k,numsteps-1); i++)
      delete [] Hcol[i];
   delete [] Hcol;

   removefourd(z,tmp.fourd,tmp.fdstatus,tmp.Nfd);

   getchar();

   if (residual >= tol)
      return 0;
   else
      return 1;
}

char GMRESpreleftsmall(double ***x, SparseElt2**** &A, double ***b, SparseElt2**** &M, 
                       GridData &grid, int numsteps, double tol, double ***S, 
                       PBData &pb, TempStruct &tmp)
{
   int i, j, k, s, m, n;
   int tindex[grid.dim], rindex[grid.dim];
   double residual, r0norm, bnorm, temp, ehere; 
   double ***z = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   double ***p = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   double ***r = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   SparseElt2 *current;

// if maxk == numsteps, then last k is numsteps-1, which creates Vcol[numsteps-1]
//    and Hcol[numsteps-2][0...numsteps-1]
   double **Hcol = new double*[numsteps-1];
   double ****Vcol = new double***[numsteps];

   leftmultILUinv(p,M,b,grid);
// get two norm of b
   bnorm = 0.0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      bnorm += evalarray(p,tindex)*evalarray(p,tindex);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   bnorm = sqrt(bnorm);

// get residual and its two norm
   leftmultmxsmall(r,A,x,grid,S,pb);
   residual = 0.0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp = evalarray(b,tindex)-evalarray(r,tindex);

      residual += temp*temp;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   residual = sqrt(residual);

   leftmultILUinv(r,M,r,grid);
   r0norm = 0.0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp = evalarray(p,tindex)-evalarray(r,tindex);

      r0norm += temp*temp;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   r0norm = sqrt(r0norm);

// get Vcol[0]
   k = -1;
   GSarnoldipreleftsmall(Vcol,Hcol,k,numsteps,x,A,p,M,S,pb,grid);
   cout << "Set initial column of V" << endl;
   residual /= bnorm;
   cout << "   residual = " << residual << endl;
// if maxk == numsteps, then last k is numsteps-1, which creates Vcol[numsteps-1]
//    and Hcol[numsteps-2][0...numsteps-1]
   double **B = matrix(numsteps-2,numsteps-2);
   double *c = new double[numsteps-1];
   double *y = new double[numsteps-1];
   int *PLR = new int[numsteps-1];
   int *PLC = new int[numsteps-1];
   double **LU = matrix(numsteps-2,numsteps-2);

   char notstop = 1;
   for ( ; k+1 < numsteps && residual >= tol && notstop; )
   {
      notstop = GSarnoldipreleftsmall(Vcol,Hcol,k,numsteps,x,A,p,M,S,pb,grid);
      cout << "Finished " << k << "th column of V" << endl;
//      cout << "B = " << endl;
      for (i = 0; i < k; i++) 
      {
         for (j = 0; j < k; j++) 
         {
            B[i][j] = 0.0;
            for (s = 0; s <= min(i+1,j+1); s++) 
               B[i][j] += Hcol[i][s]*Hcol[j][s];
         }
         c[i] = Hcol[i][0]*r0norm;
      }
      gecp0(LU,PLR,PLC,B,k-1,k-1);
      forwardbacksub0(y,c,LU,PLR,PLC,k-1);

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(z,tindex,evalarray(x,tindex));
         for (s = 0; s < k; s++)
            setvalarray(z,tindex,evalarray(z,tindex)+evalarray(Vcol[s],tindex)*y[s]);
       
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      residual = 0.0;
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         for (m = 0; m < grid.dim; m++)
            rindex[m] = tindex[m];
         if (evalarray(S,tindex) < 0.0)
            ehere = pb.epsilonm;
         else
            ehere = pb.epsilonp;
   
         temp = evalarray(b,tindex);
         if (evalarray(A,tindex) == NULL)
            for (m = 0; m < grid.dim; m++)
            {
               temp -= 2.0*ehere/(grid.dx[m]*grid.dx[m])*evalarray(z,rindex);
               for (n = -1; n <= 1; n += 2)
               {
                  rindex[m] = tindex[m]+n;
                  if (rindex[m] >= 0 && rindex[m] <= grid.nx[m])
                     temp -= -ehere/(grid.dx[m]*grid.dx[m])*evalarray(z,rindex);
               }
               rindex[m] = tindex[m];
            }
         else
            for (current = evalarray(A,tindex); current != NULL; 
                 current = (*current).next)
               temp -= (*current).val*evalarray(z,(*current).cindex);
   
         residual += temp*temp;
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      residual = sqrt(residual)/bnorm;
      cout << "   residual = " << residual << endl;
   }

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(x,tindex,evalarray(z,tindex));

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   free_matrix(B,numsteps-2,numsteps-2);
   free_matrix(LU,numsteps-2,numsteps-2);
   delete [] c;
   delete [] y;
   delete [] PLR;
   delete [] PLC;

   for (i = 0; i <= min(k,numsteps-1); i++)
      free_matrix(Vcol[i],grid.nx[0],grid.nx[1],grid.nx[2]);
   delete [] Vcol;
   for (i = 0; i < min(k,numsteps-1); i++)
      delete [] Hcol[i];
   delete [] Hcol;

   removefourd(z,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(p,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(r,tmp.fourd,tmp.fdstatus,tmp.Nfd);

   getchar();

   if (residual >= tol)
      return 0;
   else
      return 1;
}

char GMRESpreleftsmall2(double ***x, SparseElt2**** &A, double ***b, SparseElt2**** &M, 
                        GridData &grid, int numsteps, double tol, double ***S, 
                        PBData &pb, TempStruct &tmp)
{
   int i, j, k, s, m, n;
   int tindex[grid.dim], rindex[grid.dim];
   double residual, r0norm, bnorm, temp, ehere; 
   double ***z = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   double ***p = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   double ***r = setfourd(tmp.fourd,tmp.fdstatus,tmp.Nfd,grid);
   SparseElt2 *current;

// if maxk == numsteps, then last k is numsteps-1, which creates Vcol[numsteps-1]
//    and Hcol[numsteps-2][0...numsteps-1]
   double **Hcol = new double*[numsteps-1];
   double ****Vcol = new double***[numsteps];

// get two norm of b
   leftmultILUinv(p,M,b,grid);
   
   bnorm = 0.0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      bnorm += evalarray(p,tindex)*evalarray(p,tindex);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   bnorm = sqrt(bnorm);

// get residual and its two norm
   r0norm = 0.0;
   leftmultmxsmall(r,A,x,grid,S,pb);
   leftmultILUinv(r,M,r,grid);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      temp = evalarray(p,tindex)-evalarray(r,tindex);

      r0norm += temp*temp;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   r0norm = sqrt(r0norm);

// get Vcol[0]
   k = -1;
   GSarnoldipreleftsmall(Vcol,Hcol,k,numsteps,x,A,p,M,S,pb,grid);
   cout << "Set initial column of V" << endl;
   residual = r0norm/bnorm;
   cout << "   residual = " << residual << endl;
// if maxk == numsteps, then last k is numsteps-1, which creates Vcol[numsteps-1]
//    and Hcol[numsteps-2][0...numsteps-1]
   double **B = matrix(numsteps-2,numsteps-2);
   double *c = new double[numsteps-1];
   double *y = new double[numsteps-1];
   int *PLR = new int[numsteps-1];
   int *PLC = new int[numsteps-1];
   double **LU = matrix(numsteps-2,numsteps-2);

   char notstop = 1;
   for ( ; k+1 < numsteps && residual >= tol && notstop; )
   {
      notstop = GSarnoldipreleftsmall(Vcol,Hcol,k,numsteps,x,A,p,M,S,pb,grid);
      cout << "Finished " << k << "th column of V" << endl;
      for (i = 0; i < k; i++) 
      {
         for (j = 0; j < k; j++) 
         {
            B[i][j] = 0.0;
            for (s = 0; s <= min(i+1,j+1); s++) 
               B[i][j] += Hcol[i][s]*Hcol[j][s];
         }
         c[i] = Hcol[i][0]*r0norm;
      }
      gecp0(LU,PLR,PLC,B,k-1,k-1);
      forwardbacksub0(y,c,LU,PLR,PLC,k-1);

      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         setvalarray(z,tindex,evalarray(x,tindex));
         for (s = 0; s < k; s++)
            setvalarray(z,tindex,evalarray(z,tindex)+evalarray(Vcol[s],tindex)*y[s]);
       
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

      residual = 0.0;
      leftmultmxsmall(r,A,z,grid,S,pb);
      leftmultILUinv(r,M,r,grid);
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= grid.nx[0])
      {
         temp = evalarray(p,tindex)-evalarray(r,tindex);
   
         residual += temp*temp;
   
         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }
      residual = sqrt(residual)/bnorm;
      cout << "   residual = " << residual << endl;
   }

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      setvalarray(x,tindex,evalarray(z,tindex));

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   free_matrix(B,numsteps-2,numsteps-2);
   free_matrix(LU,numsteps-2,numsteps-2);
   delete [] c;
   delete [] y;
   delete [] PLR;
   delete [] PLC;

   for (i = 0; i <= min(k,numsteps-1); i++)
      free_matrix(Vcol[i],grid.nx[0],grid.nx[1],grid.nx[2]);
   delete [] Vcol;
   for (i = 0; i < min(k,numsteps-1); i++)
      delete [] Hcol[i];
   delete [] Hcol;

   removefourd(z,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(p,tmp.fourd,tmp.fdstatus,tmp.Nfd);
   removefourd(r,tmp.fourd,tmp.fdstatus,tmp.Nfd);

//   getchar();

   if (residual >= tol)
      return 0;
   else
      return 1;
}

void GMRESrestartsmall(double ***x, SparseElt2**** &A, double ***b, GridData &grid, 
                       int numsteps, double tol, double ***S, PBData &pb, 
                       TempStruct &tmp)
{
   while (!GMRESsmall(x,A,b,grid,numsteps,tol,S,pb,tmp));
}

void GMRESpreleftrestartsmall(double ***x, SparseElt2**** &A, double ***b, 
                              SparseElt2**** &M, GridData &grid, int numsteps, 
                              double tol, double ***S, PBData &pb, TempStruct &tmp)
{
   while (!GMRESpreleftsmall2(x,A,b,M,grid,numsteps,tol,S,pb,tmp));
//   while (!GMRESpreleftsmall(x,A,b,M,grid,numsteps,tol,S,pb,tmp));
}
