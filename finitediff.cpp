#include "finitediff.h"
#include "matrix.h"
#include "input.h"

// look for 4 neighboring points that does not change sign in i,j-plane to approx cross derivative.
// if exist, return as sk2, else return true, else false
// tempsk[0,1] and tempsk[2,3] each go through (-1,0) (0,1) (-1,1)
// the four points are in i-dim, temp[0] temp[1], in j-dim, temp[2] temp[3]
// e.g. index(i,j) = (a,b), 4 points (a+temp[0], b+temp[2]), (a+temp[1], b+temp[2]),(a+temp[0], b+temp[3]),(a+temp[1], b+temp[3])
// tempsk2 goes through the following sequence
// -1 1 -1 1                                                                                                             
// -1 1 -1 0                                                                                                             
// -1 1 0 1                                                                                                              
// -1 0 -1 1                                                                                                             
// -1 0 -1 0                                                                                                             
// -1 0 0 1                                                                                                              
// 0 1 -1 1                                                                                                              
// 0 1 -1 0                                                                                                              
// 0 1 0 1    
// prefer central diff
char yessk2(int *sk2, int i, int j, int *index, double ***S, GridData &grid)
{
  int combimation[9][4] = {
{-1, 1, -1, 1},
{-1, 1, -1, 0},
{-1, 1,  0, 1},
{-1, 0, -1, 1},
{ 0, 1, -1, 1},
{-1, 0, -1, 0},
{-1, 0,  0, 1},
{ 0, 1, -1, 0},
{ 0, 1,  0, 1} };

   int r, s, t, m, n, temp, rindex[grid.dim], tempsk2[4], sk2stat[4], tempsk2stat[4];
   char bad;

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];

   bad = 1;
   for(int k = 0; k < 9; k ++)
            {
// check the four points of tempsk2
               tempsk2[0] = combimation[k][0];
               tempsk2[1] = combimation[k][1];
               tempsk2[2] = combimation[k][2];
               tempsk2[3] = combimation[k][3];

               bad = 0;
               for (r = 0; r < 2 && !bad; r++)
               {
                  rindex[i] = index[i]+tempsk2[r];
                  for (s = 2; s < 4 && !bad; s++)
                  {
                     rindex[j] = index[j]+tempsk2[s];
                     if ((evalarray(S,index) < 0.0)+
                         (evalarray(S,rindex) < 0.0) == 1)
                        bad = 1;
                  }
                  rindex[j] = index[j];
               }
               rindex[i] = index[i];
               
               if (!bad){
                  for (r = 0; r < 4; r++)
                     sk2[r] = tempsk2[r];
                   return (!bad);
               }

            }

   return (!bad);
}



// used in getcim345Du, approx D2 of tindex in m,n-plane,
// write approximation from get sk2 into NxNxN coeff of u-value
void getD2(double ***D2, int m, int n, int *sk2, int *tindex, int *N, GridData grid)
{
   int i, s, t, sindex[grid.dim];

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N[0])
   {
      setvalarray(D2,sindex,0.0);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N[i]; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }

   for (s = 0; s < grid.dim; s++)
      sindex[s] = tindex[s];
   for (s = 0; s < 2; s++)
   {
      sindex[m] = tindex[m]+sk2[s];
      for (t = 0; t < 2; t++)
      {
         sindex[n] = tindex[n]+sk2[2+t];
         setvalarray(D2,sindex,(2*s-1)*(2*t-1)/((sk2[1]-sk2[0])*(sk2[3]-sk2[2])*grid.dx[m]*grid.dx[n]));
      }
      sindex[n] = tindex[n];
   }
   sindex[m] = tindex[m];
}



// approximate tion of u value
double evalcoef(double u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                int *index, int rstar, int sstar, double alpha, int mid, 
                double ***S,
                GridData grid)
{
   int i, s, N = 2*mid, tindex[grid.dim], sindex[grid.dim];
   double value = u0, thesign;

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-mid+sindex[s];
      if (evalarray(S,tindex) < 0.0)
         thesign = -1.0;
      else
         thesign = 1.0;
      value += evalarray(ucoef,sindex)*getu(tindex,0,0,0.0,thesign,grid);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

   if (evalarray(S,index) < 0.0)
      thesign = -1.0;
   else
      thesign = 1.0;
   for (s = 0; s < grid.dim; s++)
   {
      value += uxcoef[s]*getDu(index,s,rstar,sstar,alpha,thesign,grid);
      value += uxxcoef[s]*getD2u(index,s,s,rstar,sstar,alpha,thesign,grid);
   }

   return value;
}