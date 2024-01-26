#include "grid.h"
#include "global.h"
#include "helper.h"
#include "ccim.h"
#include "interface.h"
#include "numerics.h"
#include <cmath>
#include "math.h"
#include <algorithm>
#include <iostream>
#include "finitediff.h"
using namespace std;



// get mixed derivative at index in i,j-plane, stored in sk2,
// e.g. index(i,j) = (a,b), 4 points (a+sk2[0], b+sk2[2]), (a+sk2[1], b+sk2[2]),(a+sk2[0], b+sk2[3]),(a+sk2[1], b+sk2[3])
// return char 1 if found, 0 otherwise
char getsk2(int *sk2, int i, int j, int *index, double ***S, GridData &grid)
{
   int r, s, t, m, n, temp, rindex[grid.dim], tempsk2[4], sk2stat[4], tempsk2stat[4];
   char bad;

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];

   if (!globbiasstat)
      return yessk2(sk2,i,j,index,S,grid);
   else
   {
      char foundone = 0;
      // m is increment in i-dim, m=2 is central diff, m=1 is forward or backward
      // n is increment in j-dim, 
      // tempsk2[0] is first point in i-dim, tempsk2[0] = -1 means x_(i-1)
      // tempsk2[1] is second point in i-dim, tempsk2[1] = tempsk2[0] + m
      // order (-1,1), (-1,0), (0,1)
      // tempsk2[2] and tempsk2[3] are in j-dim
      for (m = 2; m >= 1 && !foundone; m--)
         for (n = 2; n >= 1 && !foundone; n--)
            for (tempsk2[0] = -1,tempsk2[1] = tempsk2[0]+m; tempsk2[1] <= 1; 
                 (tempsk2[0])++,tempsk2[1] = tempsk2[0]+m)
               for (tempsk2[2] = -1,tempsk2[3] = tempsk2[2]+n; tempsk2[3] <= 1; 
                    (tempsk2[2])++,tempsk2[3] = tempsk2[2]+n)
               {
// check the four points of tempsk2
                  bad = 0;
                  for (r = 0; r < 2 && !bad; r++) // r = 0 or 1
                  {
                     rindex[i] = index[i]+tempsk2[r];
                     for (s = 2; s < 4 && !bad; s++) // r = 2 or 3
                     {
                        rindex[j] = index[j]+tempsk2[s];
                        // if rindex goes out of bound || rindex and index have different sign, bad = 1
                        if (rindex[i] < 0 || rindex[i] > grid.nx[i] || 
                            rindex[j] < 0 || rindex[j] > grid.nx[j] || 
                            (evalarray(S,index) < 0.0)+
                            (evalarray(S,rindex) < 0.0) == 1)
                           bad = 1;
                        else
                        // if for current combimation of tempsk2, rindex is not bad(change sign)
                        // record status of rindex, sort the 4 points from large to small
                        // that is, tempsk2stat[0] > ... tempsk2stat[3]. 
                        {
// getstatus uses yessk2
                           tempsk2stat[2*r+(s-2)] = getstatus5(S,rindex,grid);
                           //goes down the list, swap if tempsk2stat[t]>tempsk2stat[t-1]
                           for (t = 2*r+(s-2); t > 0 && tempsk2stat[t] > tempsk2stat[t-1]; t--)
                           {
                              temp = tempsk2stat[t-1];
                              tempsk2stat[t-1] = tempsk2stat[t];
                              tempsk2stat[t] = temp;
                           }
                        }
                     }
                     rindex[j] = index[j];
                  }
                  rindex[i] = index[i];
                  

                  // 
                  if (!bad)
                  {
                     if (!foundone)//if not bad and haven't found 1, record current 1
                     {
                        foundone = 1;
                        for (r = 0; r < 4; r++)
                        {
                           sk2[r] = tempsk2[r];
                           sk2stat[r] = tempsk2stat[r];
                        }
                     }
                     else //if not bad but alread found 1,
                     {
                        for (r = 0; r < 4 && tempsk2stat[r] == sk2stat[r]; r++);// find the disagreement
//                        if (r < 4 && tempsk2stat[r] < sk2stat[r])
                        if (r < 4 && tempsk2stat[r] > sk2stat[r])// prefer largers status?
                           for (r = 0; r < 4; r++)
                           {
                              sk2[r] = tempsk2[r];
                              sk2stat[r] = tempsk2stat[r];
                           }
                     }
                  }
               }

      return foundone;
   }
}






//cim12.pdf p35, approximate mixed derivtive in (m,n) plane at index
//triangle stencil make use of u and uxx 
char yescim5D2(double ***D2ucoef, double *D2uxcoef, double *D2uxxcoef, int m, int n, 
               int *index, int mid, double ***S, GridData &grid)
{
   int i, r, s, thesign, N = 2*mid, ntheta = 8, signs[ntheta], sk[2], offset[2][2];
   int tindex[grid.dim], sindex[grid.dim];
   double theta0, theta, dtheta, themax, value;

   for (i = 0; i < grid.dim; i++)
      sindex[i] = 0;
   while (sindex[0] <= N)
   {
      setvalarray(D2ucoef,sindex,0.0);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

   for (r = 0; r < grid.dim; r++)
   {
      D2uxcoef[r] = 0.0;
      D2uxxcoef[r] = 0.0;
   }

   for (r = 0; r < grid.dim; r++)
      tindex[r] = index[r];
   thesign = 2*(evalarray(S,tindex) >= 0.0)-1;//1=outside, -1 indside, sign at current index
   dtheta = 2.0*M_PI/ntheta;
   theta0 = -3.0*M_PI/4.0;// start from bottom left corner
// computes signs of points surrounding node, start from lower left corner, couter-clockwise
   for (r = 0; r < ntheta; r++)
   {
      theta = theta0+r*dtheta;
      themax = max(fabs(cos(theta)),fabs(sin(theta)));
      offset[0][0] = round(cos(theta)/themax);
      offset[0][1] = round(sin(theta)/themax);
      tindex[m] = index[m]+offset[0][0];
      tindex[n] = index[n]+offset[0][1];
      signs[r] = 2*(evalarray(S,tindex) >= 0.0)-1;
   }
   tindex[m] = index[m];
   tindex[n] = index[n];

   // looks for central differencing possibility around node, "large" triangle stencil
   // r = 0, 2 ,4, 6 corresponds to bot-left, bot-right, top-right, top-left
   for (r = 0; r < ntheta; r += 2)
      // signs[r] and signs[r+2] same side as index
      if ((thesign < 0)+(signs[r] < 0) != 1 && (thesign < 0)+(signs[(r+2)%ntheta] < 0) != 1)
      {
         theta = theta0+r*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[0][0] = round(cos(theta)/themax);
         offset[0][1] = round(sin(theta)/themax);
         theta = theta0+((r+2)%ntheta)*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[1][0] = round(cos(theta)/themax);
         offset[1][1] = round(sin(theta)/themax);
         sk[0] = offset[0][0]-offset[1][0];
         sk[1] = offset[0][1]-offset[1][1];
         if (sk[0] == 0)
         {
            sk[0] = offset[0][0];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[n] = -1.0/(sk[0]*grid.dx[m]);
         }
         else
         {
            sk[1] = offset[0][1];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[m] = -1.0/(sk[1]*grid.dx[n]);
         }
         sindex[m] = mid;
         sindex[n] = mid;

         return 1;
      }

   // looks for forward differencing possibility around node, "small" triangle stencil
   for (r = 0; r < ntheta; r++)
      if ((thesign < 0)+(signs[r] < 0) != 1 && 
          (thesign < 0)+(signs[(r+1)%ntheta] < 0) != 1)
      {
         theta = theta0+r*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[0][0] = round(cos(theta)/themax);
         offset[0][1] = round(sin(theta)/themax);
         theta = theta0+((r+1)%ntheta)*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[1][0] = round(cos(theta)/themax);
         offset[1][1] = round(sin(theta)/themax);
         sk[0] = offset[0][0]-offset[1][0];
         sk[1] = offset[0][1]-offset[1][1];
         if (sk[0] == 0)
         {
            sk[0] = offset[0][0];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[n] = -1.0/(sk[0]*grid.dx[m]);
//            D2uxxcoef[n] = -0.5*sk[1]*grid.dx[n]/(sk[0]*grid.dx[m]);
            D2uxxcoef[n] = -0.5*(offset[0][1]+offset[1][1])*grid.dx[n]/
                            (sk[0]*grid.dx[m]);
         }
         else
         {
            sk[1] = offset[0][1];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[m] = -1.0/(sk[1]*grid.dx[n]);
//            D2uxxcoef[m] = -0.5*sk[0]*grid.dx[m]/(sk[1]*grid.dx[n]);
            D2uxxcoef[m] = -0.5*(offset[0][0]+offset[1][0])*grid.dx[m]/
                           (sk[1]*grid.dx[n]);
         }
         sindex[m] = mid;
         sindex[n] = mid;

         return 1;
      }

   return 0;
}


// get status of point index: 0 = boundary, 1 = interior, 2 = cim3, 3 = cim5, 4 = cim4, 5 = cim1
int getstatus5(double ***S, int *index, GridData &grid)
{
   int r, sk, m, n, mid = 1, rindex[grid.dim], sk2[4];
   int yesinterior, thecim, cimstatus[grid.dim][grid.dim];
   double ***junk1, junk2[grid.dim], junk3[grid.dim];

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
      {
         rindex[r] = index[r]+sk;
         if (rindex[r] < 0 || rindex[r] > grid.nx[r])
            return 0;//if index is boundary
         rindex[r] = index[r];
      }

   yesinterior = 1;
   for (r = 0; r < grid.dim; r++)
   {
      for (sk = -1; sk <= 1; sk += 2)
      {
         rindex[r] = index[r]+sk;
         if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1) // if near interface
            yesinterior = 0;
      }
      rindex[r] = index[r];
   }

   if (yesinterior)
      return 1;


   junk1 = matrix(2,2,2);
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         cimstatus[m][n] = 3; 
         if (!yessk2(sk2,m,n,index,S,grid))//if do not have usual mixed derivative
         {
            //if has cim5 mixed derivative: triangle stencil with ux and uxx
            if (yescim5D2(junk1,junk2,junk3,m,n,index,mid,S,grid))
               cimstatus[m][n] = 5; 
            // do not have stencil from nearby 8 direct nbrs, look for D2u approximation at neighoring grid points
            // if exist, cimstatus = 4. otherwise cimstatus = 50
            else 
            {
               cimstatus[m][n] = 4; 
               int count = 0;
               // look at each dim, front and back
               for (r = 0; r < grid.dim; r++)
                  for (sk = -1; sk <= 1; sk += 2)
                  {
                     rindex[r] = index[r]+sk;
                     // if same side and no sk2, increment count
                     if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1 &&
                         !yessk2(sk2,m,n,rindex,S,grid)){
                          ++count;
                          // cimstatus[m][n] = 0;
                      }
                     rindex[r] = index[r];
                  }
               // all up-down-left-right-front-back 6 nbrs do not have D2u approximation
                if (count==grid.dim*2){
                  cimstatus[m][n] = 0;
                }
            }
         }
      }
   free_matrix(junk1,2,2,2);

   thecim = 3;
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
         if (cimstatus[m][n] == 0)
            return 5;
         else if (cimstatus[m][n] == 4)
            thecim = 4;
         else if (thecim == 3 && cimstatus[m][n] == 5)
            thecim = 5;


   if (thecim == 3)
   {
// second derivatives exhibit cim 3 structure but have to still check first derivatives
      int count = 0;

      for (r = 0; r < grid.dim; r++)
         rindex[r] = index[r];
      for (r = 0; r < grid.dim && count < 2; r++)
      {
         count = 0;
         for (sk = -1; sk <= 1; sk += 2)
         {
            rindex[r] = index[r]+sk;
            if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1)
               count++;
         }
         rindex[r] = index[r];
      }
      if (count < 2)
         return 2;
      else
         thecim = 5;
   }

   if (thecim == 5)//
      return 3;
   else if (thecim == 4)
      return 4;
   else
      return 5;
}

/*
approximate D2u in (m,n) plane at index in terms of constant u0, NxNxN u coeff, 3x1 uxcoeff, 3x1 uxx coef, 1x1 jumpuxx coef 
In the first pass, rstar=0, sstar=0, return perm=1 if can be approx with same side, otherwise perm=0. Record coeff.
In the second pass,
(1) if the first pass return 0, in the second pass, provide extra info on interface location, see if use uxy on the other side (cim4)
(2) if the first pass return 1, the coeffs are already recorded, in the second pass, nothing needs to be done, so the whole block is skipped
*/
void getcim345D2u(double &u0, double ***u, double *uxcoef, double *uxxcoef, 
                  double &jumpuxxcoef, char &perm, int m, int n, int *index, int rstar, 
                  int sstar, int mid, double ***S, GridData &grid)
{
   int t, N = 2*mid, sindex[grid.dim], rindex[grid.dim], nindex[grid.dim], sk2[4];
   double thesign;

   if (!perm)
   {
      if (getsk2(sk2,m,n,index,S,grid))//if D2u in (m,n) plane at index has approximation in same side
      {
         for (t = 0; t < grid.dim; t++)
         {
            sindex[t] = mid;
            nindex[t] = N;
         }
         getD2(u,m,n,sk2,sindex,nindex,grid);// approximate D2u in terms of u-value only
         u0 = 0.0;
         for (t = 0; t < grid.dim; t++)
         {
            uxcoef[t] = 0.0;
            uxxcoef[t] = 0.0;
         }
         jumpuxxcoef = 0.0;
         perm = 1;
      }
      else if (getcim5D2(u,uxcoef,uxxcoef,m,n,index,mid,S,grid)) // use triangle stencil with ux and uxx
      {
         u0 = 0.0;
         jumpuxxcoef = 0.0;
         perm = 1;
      }
      else 
      {
         int r, s;
         int **anosk2 = imatrix(1,3);//2x4 matrix, anosk2[0] is D2u stencil for s=-1, anosk2[1] for s=1
         char yes[2];

         if (globdirectD2) // approximate D2 from nbr of same side, prefer out of plane
         {
            for (t = 0; t < grid.dim; t++)
               rindex[t] = index[t];
            for (r = 0; r < grid.dim && (r == m || r == n); r++); // find r no equal to m or n
            for (s = -1; s <= 1; s += 2)
            {
               rindex[r] = index[r]+s; //rindex is nbr out of m,n-plane
               if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) != 1)//if change sign
                  yes[(s+1)/2] = getsk2(anosk2[(s+1)/2],m,n,rindex,S,grid); 
                  //yes[0] = 1 if s=-1, rindex has approx of mixed derivative
               else
                  yes[(s+1)/2] = 0;
            }
         }

         // if we could found s = +/-1, same side, out-of-plane, have D2u approx with u-value
         if (globdirectD2 && (yes[0] || yes[1])) 
         {  
            // with use s=-1 if 
            // (1) s=-1 has sk2 but s=1 no sk2, or 
            // (2) s=-1 has sk2, s=1 has sk2, s=-1 has larger stencil,i.e. prefer central differencing
            if ((yes[0] && !(yes[1])) ||
               (yes[0] && yes[1] && 
               abs(anosk2[0][0]-anosk2[0][1])+abs(anosk2[0][2]-anosk2[0][3])>
               abs(anosk2[1][0]-anosk2[1][1])+abs(anosk2[1][2]-anosk2[1][3])))              
               s = -1;// use s=-1
            else
               s = 1;
            for (t = 0; t < grid.dim; t++)
            {
               sindex[t] = mid;
               nindex[t] = N;
            }
            sindex[r] = mid+s;
            getD2(u,m,n,anosk2[(s+1)/2],sindex,nindex,grid);
            u0 = 0.0;
            for (t = 0; t < grid.dim; t++)
            {
               uxcoef[t] = 0.0;
               uxxcoef[t] = 0.0;
            }
            jumpuxxcoef = 0.0;
            perm = 1;
         }
         // cannot find any approx in the same side, but rstar and sstar is provided, i.e. know where is interface
         // approxiamte D2u at points on different side, make use of jumpuxx
         else if (rstar >= 0 && rstar < grid.dim && sstar != 0)
         {
            cout << "Using cim4" << endl;//use mix derivative from the other side
            for (t = 0; t < grid.dim; t++)
               rindex[t] = index[t];
            rindex[rstar] = index[rstar]+sstar;
            if (getsk2(sk2,m,n,rindex,S,grid))
            {
               if (evalarray(S,index) < 0.0)
                  thesign = -1;
               else
                  thesign = 1;
               for (t = 0; t < grid.dim; t++)
               {
                  sindex[t] = mid;
                  nindex[t] = N;
               }
               sindex[rstar] = mid+sstar;
               getD2(u,m,n,sk2,sindex,nindex,grid);
               jumpuxxcoef = thesign;
               u0 = 0.0;
               for (t = 0; t < grid.dim; t++)
               {
                  uxcoef[t] = 0.0;
                  uxxcoef[t] = 0.0;
               }
               perm = 0;
            }
            else
            {
               cout << "bad status in getcim345D2" << endl;
               exit(1);
            }
            rindex[rstar] = index[rstar];
         }
         //cannot find any approx in the same side, sstar==0, i.e. don't know where is interface
         else
         {
            perm = 0;
         }
         free_matrix(anosk2,1,3);
      }


   }
/*
   if (!perm)
      if (getsk2(sk2,m,n,index,S,grid))
      {
         for (t = 0; t < grid.dim; t++)
         {
            sindex[t] = mid;
            nindex[t] = N;
         }
         getD2(u,m,n,sk2,sindex,nindex,grid);
         u0 = 0.0;
         for (t = 0; t < grid.dim; t++)
         {
            uxcoef[t] = 0.0;
            uxxcoef[t] = 0.0;
         }
         jumpuxxcoef = 0.0;
         perm = 1;
      }
      else if (getcim5D2(u,uxcoef,uxxcoef,m,n,index,mid,S,grid))
      {
         u0 = 0.0;
         jumpuxxcoef = 0.0;
         perm = 1;
      }
      else if (rstar >= 0 && rstar < grid.dim && sstar != 0)
      {
         for (t = 0; t < grid.dim; t++)
            rindex[t] = index[t];
         rindex[rstar] = index[rstar]+sstar;
         if (getsk2(sk2,m,n,rindex,S,grid))
         {
            if (evalarray(S,index) < 0.0)
               thesign = -1;
            else
               thesign = 1;
            for (t = 0; t < grid.dim; t++)
            {
               sindex[t] = mid;
               nindex[t] = N;
            }
            sindex[rstar] = mid+sstar;
            getD2(u,m,n,sk2,sindex,nindex,grid);
            jumpuxxcoef = thesign;
            u0 = 0.0;
            for (t = 0; t < grid.dim; t++)
            {
               uxcoef[t] = 0.0;
               uxxcoef[t] = 0.0;
            }
            perm = 0;
         }
         else
         {
            cout << "bad status in getcim345D2" << endl;
            exit(1);
         }
         rindex[rstar] = index[rstar];
      }
      else
         perm = 0;
*/
}


// cim5D2 with bias,
char getcim5D2(double ***D2ucoef, double *D2uxcoef, double *D2uxxcoef, int m, int n, 
               int *index, int mid, double ***S, GridData &grid)
{
   int i, r, s, thesign, N = 2*mid, ntheta = 8, signs[ntheta], sk[2], offset[2][2];
   int tindex[grid.dim], sindex[grid.dim];
   double theta0, theta, dtheta, themax, value;

   if (!globbiasstat)
      yescim5D2(D2ucoef,D2uxcoef,D2uxxcoef,m,n,index,mid,S,grid);
   else
   {
      int t, temp, foundone = 0, stat[2], tmpstat[2], tmpoffset[2][2];
      int rindex[grid.dim];
      for (r = 0; r < grid.dim; r++)
         rindex[r] = index[r];

      for (i = 0; i < grid.dim; i++)
         sindex[i] = 0;
      while (sindex[0] <= N)
      {
         setvalarray(D2ucoef,sindex,0.0);
   
         (sindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
         {
            sindex[i] = 0;
            (sindex[i-1])++;
         }
      }
      for (s = 0; s < grid.dim; s++)
         sindex[s] = mid;

      for (r = 0; r < grid.dim; r++)
      {
         D2uxcoef[r] = 0.0;
         D2uxxcoef[r] = 0.0;
      }

      for (r = 0; r < grid.dim; r++)
         tindex[r] = index[r];
      thesign = 2*(evalarray(S,tindex) >= 0.0)-1;
      dtheta = 2.0*M_PI/ntheta;
      theta0 = -3.0*M_PI/4.0;
// computes signs of points surrounding node, start with lower left, move counter-clockwise,
      for (r = 0; r < ntheta; r++)
      {
         theta = theta0+r*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[0][0] = round(cos(theta)/themax);
         offset[0][1] = round(sin(theta)/themax);
         tindex[m] = index[m]+offset[0][0];
         tindex[n] = index[n]+offset[0][1];
         signs[r] = 2*(evalarray(S,tindex) >= 0.0)-1;
      }
      tindex[m] = index[m];
      tindex[n] = index[n];
   
// looks for central and forward differencing possibilities around node
      for (t = 2; t >= 1 && !foundone; t--)// t is step size, if t=2 have found 1, no need to consider t = 1
         for (r = 0; r < ntheta; r += t) //r is starting index, when t=2, the kernel is 1/4 triangle, when t = 1, the kernel is 1/8 triangle
            if ((thesign < 0)+(signs[r] < 0) != 1 && 
                (thesign < 0)+(signs[(r+t)%ntheta] < 0) != 1) // consider valid combinations
            {
               theta = theta0+r*dtheta;
               themax = max(fabs(cos(theta)),fabs(sin(theta)));
               tmpoffset[0][0] = round(cos(theta)/themax);
               tmpoffset[0][1] = round(sin(theta)/themax);
               theta = theta0+((r+t)%ntheta)*dtheta;
               themax = max(fabs(cos(theta)),fabs(sin(theta)));
               tmpoffset[1][0] = round(cos(theta)/themax);
               tmpoffset[1][1] = round(sin(theta)/themax);

               for (s = 0; s < 2; s++)
               {
                  rindex[m] = tindex[m]+tmpoffset[s][0];
                  rindex[n] = tindex[n]+tmpoffset[s][1];
                  tmpstat[s] = getstatus5(S,rindex,grid); // tempstat is the status of the two points
               }
               if (tmpstat[1] > tmpstat[0]) // order tempstat
               {
                  temp = tmpstat[0];
                  tmpstat[0] = tmpstat[1];
                  tmpstat[1] = temp;
               }
               rindex[m] = tindex[m];
               rindex[n] = tindex[n];
               if (!foundone) // if first time found available stencil
               {
                  foundone = t;
                  for (s = 0; s < 2; s++)
                  {
                     offset[s][0] = tmpoffset[s][0];
                     offset[s][1] = tmpoffset[s][1];
                     stat[s] = tmpstat[s];
                  }
               }
               else // if alread have stencil,
               {
                  for (s = 0; s < 2 && tmpstat[s] == stat[s]; s++); // tempstat and stat are ordered, s stops at fist difference
//                  if (s < 2 && tmpstat[s] < stat[s])
                  if (s < 2 && tmpstat[s] > stat[s])// if tempstat is larger, use tempstat
                     for (s = 0; s < 2; s++)
                     {
                        offset[s][0] = tmpoffset[s][0];
                        offset[s][1] = tmpoffset[s][1];
                        stat[s] = tmpstat[s];
                     }
               }
            }
      if (foundone)
      {
         sk[0] = offset[0][0]-offset[1][0];
         sk[1] = offset[0][1]-offset[1][1];
         if (sk[0] == 0)// when the two points are along this dimension
         {
            sk[0] = offset[0][0];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[n] = -1.0/(sk[0]*grid.dx[m]);
            if (foundone == 1)
               D2uxxcoef[n] = -0.5*(offset[0][1]+offset[1][1])*grid.dx[n]/
                               (sk[0]*grid.dx[m]);
         }
         else
         {
            sk[1] = offset[0][1];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[m] = -1.0/(sk[1]*grid.dx[n]);
            if (foundone == 1)
               D2uxxcoef[m] = -0.5*(offset[0][0]+offset[1][0])*grid.dx[m]/
                              (sk[1]*grid.dx[n]);
         }
         sindex[m] = mid;
         sindex[n] = mid;

         return 1;
      }
      else
         return 0;
   }

   cout << "getcim5D2 NOT SUPPOSED TO BE HERE" << endl;
   return 0;
}


/*

Input: 
Du- in terms of (const, u-val, Du, principal D2u, jumpD2u (can have mixed))
Du-_r = D1ucoef[r][:sindex] u_:sindex +  D1uxcoef[r][:s] Du_:s + D1uxxcoef[r][:s] D2u_:ss + D1jumpuxxcoef[r][:m][:n] jumpD2u_:m:n

Output:
get jumpDu_rstar at interface in dim rstar direction sstar
jumpDu_rstar = u0 + ucoef[:sindex] u_:sindex +  uxcoef[:s] Du_:s + uxxoef[:s] D2u_:ss + jumpuxxcoef[:m][:n] jumpD2u_:m:n


Method:
jumpDu depends on Du- and interface condition, cim12.pdf p27

*/
void getcim345jumpux(double& u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                     double **jumpuxxcoef, int *index, int rstar, int sk, double alpha, 
                     int thesign, double *normal, double *tangent, int mid, 
                     double ****D1ucoef, double **D1uxcoef, double **D1uxxcoef, 
                     double ***D1jumpuxxcoef, double ***S, PBData &pb, GridData &grid)
{
   int i, s, t, m, n, N = 2*mid, sindex[grid.dim];
   double ehere, ethere, sigma, Dtau[grid.dim], value;

   if (evalarray(S,index) < 0.0)
   {
      ehere = pb.epsilonm;
      ethere = pb.epsilonp;
   }
   else
   {
      ehere = pb.epsilonp;
      ethere = pb.epsilonm;
   }
   getsigma(sigma,index,rstar,sk,alpha,normal,pb,grid);
   getDtau(Dtau,index,rstar,sk,alpha,grid);
   
   // constant term, only one tangent vector is needed, as the other tangent vector (tangent2) can be chosen
   // to be in the plane with normal e[rstar], so tangent2[rstar] = 0
   u0 = sigma/ethere*normal[rstar] + getdotprod(Dtau,tangent,grid.dim)*tangent[rstar];

   // recast Du-, [epsilon]/epsilonp dot(grad(Du-),n) n[star]
   // thesign*(ethere-ehere): if thesign = 1, then epsm-epsp. if thesign = -1, then -1 (epsp-epsm), therefore awlays - (epsp-epsm) 
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      value = 0.0;
      for (t = 0; t < grid.dim; t++)
         value += evalarray(D1ucoef[t],sindex)*normal[t];
      setvalarray(ucoef,sindex,thesign*(ethere-ehere)/ethere*value*normal[rstar]);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }

   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

   for (s = 0; s < grid.dim; s++)
   {
      value = 0.0;
      for (t = 0; t < grid.dim; t++)
         value += D1uxcoef[t][s]*normal[t];
      uxcoef[s] = thesign*(ethere-ehere)/ethere*value*normal[rstar];
      value = 0.0;
      for (t = 0; t < grid.dim; t++)
         value += D1uxxcoef[t][s]*normal[t];
      uxxcoef[s] = thesign*(ethere-ehere)/ethere*value*normal[rstar];
   }

   for (m = 0; m < grid.dim; m++)
      for (n = 0; n < grid.dim; n++)
         jumpuxxcoef[m][n] = 0.0;
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         value = 0.0;
         for (t = 0; t < grid.dim; t++)
            value += D1jumpuxxcoef[t][m][n]*normal[t];
         jumpuxxcoef[m][n] = thesign*(ethere-ehere)/ethere*value*normal[rstar];
      }
}


// without a
void getcim345jumpuxx(double &u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                      int *index, int rstar, int sk, double alpha, int thesign, 
                      double *normal, int mid, double *D1u, double ****D1ucoef, 
                      double **D1uxcoef, double **D1uxxcoef, double ***D1jumpuxxcoef, 
                      double **D2u, double *****D2ucoef, double ***D2uxcoef, 
                      double ***D2uxxcoef, double **D2jumpuxxcoef, double &jumpD1u,
                      double ***jumpD1ucoef, double *jumpD1uxcoef, 
                      double *jumpD1uxxcoef, double **jumpD1jumpuxxcoef, double ***S, 
                      PBData &pb, GridData &grid)
{
   int i, r, s, n, m, N = 2*mid;
   int rindex[grid.dim], sindex[grid.dim];
   double b0[2*grid.dim], bxcoef[2*grid.dim][grid.dim], bxxcoef[2*grid.dim][grid.dim]; 
   double ***bcoef[2*grid.dim];
   
   double **Dn; 
   double Dndott[grid.dim], Dndots[grid.dim], dotDndot[grid.dim];
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double tangent1[grid.dim], tangent2[grid.dim];
   double ethere, ehere, jumpfe;
   double value, temp[2*grid.dim];
   double sigma, tau, Dsigma[grid.dim];
   double Dtau[grid.dim], **D2tau, dotD2taudot[grid.dim], Dtaudot[grid.dim];
   int two2one[grid.dim][grid.dim];
   char theorder = 2;

   double **LU = matrix(2*grid.dim-1,2*grid.dim-1);
   double **B = matrix(2*grid.dim-1,2*grid.dim-1);

   Dn = matrix(grid.dim-1,grid.dim-1);
   D2tau = matrix(grid.dim-1,grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      bcoef[r] = matrix(N,N,N);

   if (evalarray(S,index) < 0.0)
   {
      ehere = pb.epsilonm;
      ethere = pb.epsilonp;
   }
   else
   {
      ehere = pb.epsilonp;
      ethere = pb.epsilonm;
   }

   getinterfaceinfo(tangent1,tangent2,sigma,Dn,Dsigma,jumpfe,index,rstar,sk,alpha,S,
                    pb,grid);
// form Dn dot product with various vectors
   getDtau(Dtau,index,rstar,sk,alpha,grid);
   getD2tau(D2tau,index,rstar,sk,alpha,grid);
   for (n = 0; n < grid.dim; n++)
   {
      Dndott[n] = 0.0;
      Dndots[n] = 0.0;
      dotDndot[n] = 0.0;
      dotD2taudot[n] = 0.0;
      Dtaudot[n] = 0.0;
   }
   for (n = 0; n < grid.dim; n++)
   {
      for (m = 0; m < grid.dim; m++)
      {
         Dndott[n] += Dn[n][m]*tangent1[m]; //Dndott = Dn tangent1
         Dndots[n] += Dn[n][m]*tangent2[m]; //Dndots = Dn tangent2
         dotD2taudot[0] += tangent1[n]*D2tau[n][m]*tangent1[m]; //dotD2taudot[0] = tangent1 D2tau tangent1
         dotD2taudot[1] += tangent2[n]*D2tau[n][m]*tangent2[m]; //dotD2taudot[1] = tangent2 D2tau tangent2
         dotD2taudot[2] += tangent1[n]*D2tau[n][m]*tangent2[m]; //dotD2taudot[2] = tangent1 D2tau tangent2
      }
      dotDndot[0] += tangent1[n]*Dndott[n]; //dotDndot[0] = tangent1 Dn tangent1
      dotDndot[1] += tangent2[n]*Dndots[n]; //dotDndot[1] = tangent2 Dn tangent2
      dotDndot[2] += tangent1[n]*Dndots[n]; //dotDndot[2] = tangent1 Dn tangent2
      Dtaudot[0] += Dtau[n]*normal[n]; //Dtaudot[0] = Dtau n
      Dtaudot[1] += Dtau[n]*tangent1[n]; //Dtaudot[1] = Dtau tangent1
      Dtaudot[2] += Dtau[n]*tangent2[n]; //Dtaudot[2] = Dtau tangent2
   }
// form matrix
   for (n = 0; n < grid.dim; n++)
      two2one[n][n] = n;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (m = 0; m < grid.dim; m++)
      for (n = 0; n < m; n++)
         two2one[m][n] = two2one[n][m];

   for (n = 0; n < grid.dim; n++)
   {
      B[0][n] = tangent1[n]*tangent1[n];
      B[1][n] = tangent2[n]*tangent2[n];
      B[2][n] = tangent1[n]*tangent2[n];
      B[3][n] = normal[n]*tangent1[n];
      B[4][n] = normal[n]*tangent2[n];
      B[5][n] = 1.0;
   }
   for (n = grid.dim; n < 2*grid.dim; n++)
   {
      m = n-grid.dim;
      s = (m+1)%grid.dim;
      B[0][n] = 2.0*tangent1[m]*tangent1[s];
      B[1][n] = 2.0*tangent2[m]*tangent2[s];
      B[2][n] = tangent1[m]*tangent2[s]+tangent1[s]*tangent2[m];
      B[3][n] = normal[m]*tangent1[s]+normal[s]*tangent1[m];
      B[4][n] = normal[m]*tangent2[s]+normal[s]*tangent2[m];
      B[5][n] = 0.0;
   }
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
// get contributions of Du to jump data
         value = 0.0;
         for (r = 0; r < grid.dim; r++)
            value += D1jumpuxxcoef[r][m][n]*normal[r];// dot(Du-,n)
         for (r = 0; r < grid.dim; r++)
            B[r][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*value*dotDndot[r];
         value = 0.0;
         for (r = 0; r < grid.dim; r++)
            value += D1jumpuxxcoef[r][m][n]*Dndott[r]; // Du- Dn tangent1
         B[grid.dim][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*value;
         value = 0.0;
         for (r = 0; r < grid.dim; r++)
            value += D1jumpuxxcoef[r][m][n]*Dndots[r]; // Du- Dn tangent2
         B[grid.dim+1][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*value;
// get contributions of D2u to jump data
         B[grid.dim][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*
                                       (normal[m]*tangent1[n]+normal[n]*tangent1[m])*
                                       D2jumpuxxcoef[m][n]; //n D2u- tangent1
         B[grid.dim+1][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*
                                         (normal[m]*tangent2[n]+normal[n]*tangent2[m])*
                                         D2jumpuxxcoef[m][n];//n D2u- tangent2
      }
  // constant term rhs
   for (n = 0; n < grid.dim; n++)
      b0[n] = (sigma/ethere-Dtaudot[0])*dotDndot[n]+dotD2taudot[n];//?
   b0[grid.dim] = getdotprod(Dsigma,tangent1,grid.dim)/ethere-
                  Dtaudot[1]*dotDndot[0]-Dtaudot[2]*dotDndot[2];
   b0[grid.dim+1] = getdotprod(Dsigma,tangent2,grid.dim)/ethere-
                    Dtaudot[1]*dotDndot[2]-Dtaudot[2]*dotDndot[1];
   b0[grid.dim+2] = -jumpfe;

// get b coefs, coef of u-value, b[r] = NxNxN matrix for r = 0,...5
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
// initialize to zero
      for (m = 0; m < 2*grid.dim; m++)
         setvalarray(bcoef[m],sindex,0.0);

// get contributions of Du at interface
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += evalarray(D1ucoef[n],sindex)*normal[n]; //Du- dot n
      for (m = 0; m < grid.dim; m++)
          // [eps]/epsp (Du- dot n) (si Dn sj)
         setvalarray(bcoef[m],sindex,evalarray(bcoef[m],sindex)+
                                     thesign*(ethere-ehere)/ethere*value*dotDndot[m]); 
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += evalarray(D1ucoef[n],sindex)*Dndott[n]; //Du- Dn s
       //[eps]/epsp Du- Dn s
      setvalarray(bcoef[grid.dim],sindex,evalarray(bcoef[grid.dim],sindex)+
                                  thesign*(ethere-ehere)/ethere*value);
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += evalarray(D1ucoef[n],sindex)*Dndots[n];
      setvalarray(bcoef[grid.dim+1],sindex,evalarray(bcoef[grid.dim+1],sindex)+
                                    thesign*(ethere-ehere)/ethere*value);

// get contributions of D2u
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            // -[eps]/epsp n D2u- s1
            setvalarray(bcoef[grid.dim],sindex,
                        evalarray(bcoef[grid.dim],sindex)+
                        thesign*(ethere-ehere)/ethere*
                        (normal[m]*tangent1[n]+normal[n]*tangent1[m])*
                         evalarray(D2ucoef[m][n],sindex));
            // -[eps]/epsp n D2u- s2
            setvalarray(bcoef[grid.dim+1],sindex,
                        evalarray(bcoef[grid.dim+1],sindex)+
                        thesign*(ethere-ehere)/ethere*
                        (normal[m]*tangent2[n]+normal[n]*tangent2[m])*
                         evalarray(D2ucoef[m][n],sindex));
         }

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

// get bx and bxx coefs, bxcoef[r] = 1x3 matrix for Du_{1,2,3}, bxxcoef[r] = 1x3 matrix for DDu_{11,22,33}
   for (s = 0; s < grid.dim; s++)
   {
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxcoef[n][s]*normal[n]; // Du- n
      for (m = 0; m < grid.dim; m++)
         bxcoef[m][s] = thesign*(ethere-ehere)/ethere*value*dotDndot[m];
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxcoef[n][s]*Dndott[n];
      bxcoef[grid.dim][s] = thesign*(ethere-ehere)/ethere*value;
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxcoef[n][s]*Dndots[n];
      bxcoef[grid.dim+1][s] = thesign*(ethere-ehere)/ethere*value;
      bxcoef[grid.dim+2][s] = 0.0;
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            bxcoef[grid.dim][s] += thesign*(ethere-ehere)/ethere*
                                   (normal[m]*tangent1[n]+normal[n]*tangent1[m])*
                                   D2uxcoef[m][n][s];
            bxcoef[grid.dim+1][s] += thesign*(ethere-ehere)/ethere*
                                     (normal[m]*tangent2[n]+normal[n]*tangent2[m])*
                                     D2uxcoef[m][n][s];
         }

      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxxcoef[n][s]*normal[n];
      for (m = 0; m < grid.dim; m++)
         bxxcoef[m][s] = thesign*(ethere-ehere)/ethere*value*dotDndot[m];
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxxcoef[n][s]*Dndott[n];
      bxxcoef[grid.dim][s] = thesign*(ethere-ehere)/ethere*
                             (value+normal[s]*tangent1[s]);
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxxcoef[n][s]*Dndots[n];
      bxxcoef[grid.dim+1][s] = thesign*(ethere-ehere)/ethere*
                               (value+normal[s]*tangent2[s]);
      bxxcoef[grid.dim+2][s] = 0.0;
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            bxxcoef[grid.dim][s] += thesign*(ethere-ehere)/ethere*
                                    (normal[m]*tangent1[n]+normal[n]*tangent1[m])*
                                    D2uxxcoef[m][n][s];
            bxxcoef[grid.dim+1][s] += thesign*(ethere-ehere)/ethere*
                                      (normal[m]*tangent2[n]+normal[n]*tangent2[m])*
                                      D2uxxcoef[m][n][s];
         }
   }

   // check uxx matrix
   // double realjumpuxx[6], realux[3],dotrealuxxdot[2];
   //  for (n = 0; n < grid.dim; n++)
   //    realjumpuxx[n] = getD2u(index,n,n,rstar,sk,alpha,1,grid)-
   //                     getD2u(index,n,n,rstar,sk,alpha,-1,grid);
   // realjumpuxx[grid.dim] = getD2u(index,0,1,rstar,sk,alpha,1,grid)-
   //                         getD2u(index,0,1,rstar,sk,alpha,-1,grid);
   // realjumpuxx[grid.dim+1] = getD2u(index,1,2,rstar,sk,alpha,1,grid)-
   //                           getD2u(index,1,2,rstar,sk,alpha,-1,grid);
   // realjumpuxx[grid.dim+2] = getD2u(index,0,2,rstar,sk,alpha,1,grid)-
   //                           getD2u(index,0,2,rstar,sk,alpha,-1,grid);
   // for (n = 0; n < grid.dim; n++)
   //    realux[n] = getDu(index,n,rstar,sk,alpha,thesign,grid);
   // dotrealuxxdot[0] = 0.0;
   // dotrealuxxdot[1] = 0.0;
   // for (n = 0; n < grid.dim; n++)
   //    for (m = 0; m < grid.dim; m++)
   //    {
   //       dotrealuxxdot[0] += normal[n]*
   //                           getD2u(index,n,m,rstar,sk,alpha,thesign,grid)*
   //                           tangent1[m];
   //       dotrealuxxdot[1] += normal[n]*
   //                           getD2u(index,n,m,rstar,sk,alpha,thesign,grid)*
   //                           tangent2[m];
   //    }
   // cout << index[0] << " " << index[1] << " " << index[2] << endl;
   
   // cout << "   checking " << getdotprod(B[0],realjumpuxx,2*grid.dim) << " "
   //      << (sigma/ethere + thesign*(ethere-ehere)/ethere*getdotprod(realux,normal,grid.dim) - Dtaudot[0])*dotDndot[0]+dotD2taudot[0]<<" "<<dotD2taudot[0] << endl;
   
   // cout << "   checking " << getdotprod(B[1],realjumpuxx,2*grid.dim) << " "
   //      << (sigma/ethere + thesign*(ethere-ehere)/ethere*getdotprod(realux,normal,grid.dim) - Dtaudot[0])*dotDndot[1]+dotD2taudot[1]<<" "<<dotD2taudot[1] << endl;
   
   // cout << "   checking " << getdotprod(B[2],realjumpuxx,2*grid.dim) << " "
   //      << (sigma/ethere + thesign*(ethere-ehere)/ethere*getdotprod(realux,normal,grid.dim) - Dtaudot[0])*dotDndot[2]+dotD2taudot[2]<<" "<<dotD2taudot[2] << endl;
   
   // cout << "   checking " << getdotprod(B[3],realjumpuxx,2*grid.dim) << " "
   //      << (getdotprod(Dsigma,tangent1,grid.dim) + thesign*(ethere-ehere)*(dotrealuxxdot[0] + getdotprod(realux,Dndott,grid.dim)))/ethere - Dtaudot[1]*dotDndot[0] - Dtaudot[2]*dotDndot[2] << endl;
   
   // cout << "   checking " << getdotprod(B[4],realjumpuxx,2*grid.dim) << " "
   //      << (getdotprod(Dsigma,tangent2,grid.dim) + thesign*(ethere-ehere)*(dotrealuxxdot[1] + getdotprod(realux,Dndots,grid.dim)))/ethere - Dtaudot[1]*dotDndot[2] - Dtaudot[2]*dotDndot[1] << endl;
   
   // cout << "   checking " << getdotprod(B[5],realjumpuxx,2*grid.dim) << " "
   //      << -jumpfe << endl;

   // getchar();

   gecp0(LU,PLR,PLC,B,2*grid.dim-1,2*grid.dim-1);

// form jumpuxx in rstar direction and also recast Du and D2u
   forwardbacksub0(temp,b0,LU,PLR,PLC,2*grid.dim-1);
   u0 = temp[two2one[rstar][rstar]];
   //recast constant term of uxx in D1jumpuxxcoef to D1u
   for (s = 0; s < grid.dim; s++)
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
            D1u[s] += D1jumpuxxcoef[s][m][n]*temp[two2one[m][n]];
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         D2u[m][n] += D2jumpuxxcoef[m][n]*temp[two2one[m][n]];//recast D2jumpuxxcoef to D2u
         jumpD1u += jumpD1jumpuxxcoef[m][n]*temp[two2one[m][n]];//recast jumpD1jumpuxxcoef to jumpD1u
      }

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      for (m = 0; m < 2*grid.dim; m++)
         temp[m] = evalarray(bcoef[m],sindex); 
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1); // inv(B) temp[:] u[sindex]
      setvalarray(ucoef,sindex,temp[two2one[rstar][rstar]]); //two2one[rstar][rstar] = rstar
      // recast D1jumpuxxcoef to D1ucoef
      for (s = 0; s < grid.dim; s++)
         for (m = 0; m < grid.dim; m++)
            for (n = m+1; n < grid.dim; n++)
               setvalarray(D1ucoef[s],sindex,evalarray(D1ucoef[s],sindex)+
                                             D1jumpuxxcoef[s][m][n]*
                                             temp[two2one[m][n]]);
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            // recast D2jumpuxxcoef to D2ucoef
            setvalarray(D2ucoef[m][n],sindex,evalarray(D2ucoef[m][n],sindex)+
                                             D2jumpuxxcoef[m][n]*temp[two2one[m][n]]);
            // recat jumpD1jumpuxxcoef to jumpD1ucoef
            setvalarray(jumpD1ucoef,sindex,evalarray(jumpD1ucoef,sindex)+
                                           jumpD1jumpuxxcoef[m][n]*temp[two2one[m][n]]);
         }

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;
    
   for (s = 0; s < grid.dim; s++)
   {
      for (m = 0; m < 2*grid.dim; m++)
         temp[m] = bxcoef[m][s]; // temp[:] Du_{s}
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      uxcoef[s] = temp[two2one[rstar][rstar]]; //[DDu_{rstar}] = ... uxcoef[s] Du{s} ...
      // recast D1jumpuxxcoef to D1uxcoef
      for (r = 0; r < grid.dim; r++)
         for (m = 0; m < grid.dim; m++)
            for (n = m+1; n < grid.dim; n++)
               D1uxcoef[r][s] += D1jumpuxxcoef[r][m][n]*temp[two2one[m][n]];
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            D2uxcoef[m][n][s] += D2jumpuxxcoef[m][n]*temp[two2one[m][n]];
            jumpD1uxcoef[s] += jumpD1jumpuxxcoef[m][n]*temp[two2one[m][n]];
         }
      for (m = 0; m < 2*grid.dim; m++)
         temp[m] = bxxcoef[m][s]; // temp[:] DDu_{ss}
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      uxxcoef[s] = temp[two2one[rstar][rstar]];//[DDu_{rstar}] = ... uxxcoef[s] DDu_{s} ...
      for (r = 0; r < grid.dim; r++)
         for (m = 0; m < grid.dim; m++)
            for (n = m+1; n < grid.dim; n++)
               D1uxxcoef[r][s] += D1jumpuxxcoef[r][m][n]*temp[two2one[m][n]];
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            D2uxxcoef[m][n][s] += D2jumpuxxcoef[m][n]*temp[two2one[m][n]];
            jumpD1uxxcoef[s] += jumpD1jumpuxxcoef[m][n]*temp[two2one[m][n]];
         }
   }
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         for (s = 0; s < grid.dim; s++)
            D1jumpuxxcoef[s][m][n] = 0.0;
         D2jumpuxxcoef[m][n] = 0.0;
         jumpD1jumpuxxcoef[m][n] = 0.0;
      }

   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(B,2*grid.dim-1,2*grid.dim-1);
   free_matrix(Dn,grid.dim-1,grid.dim-1);
   free_matrix(D2tau,grid.dim-1,grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(bcoef[r],N,N,N);
}




/*
solve the linear system of (mixed) jumpD2u, rhs only involve const, u-val, Du, principal D2u
only jumpD2u_rstar,rstar is needed as output
recast jumpuxxcoeff to ucoef, uxcoef and uxxcoef

jumpD1u_rstar might depends on (mixed) jumpD2u, back sub to eliminate it

output:
(1) jumpD2u_rstar,rstar = u0 + ucoef[:sindex] u_sindex + uxcoef[:s] Du_:s + uxxcoef[:s] D2u_:ss

(3) back substitution to eliminate  D1jumpuxxcoef, D2jumpuxxcoef, jumpD1jumpuxxcoef in the input

after this call, every thing depends on const, u, D1u, principle D2u,                      

called with
getcim345jumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,
                             alpha[r][(sk+1)/2],thesign,normal,mid,a,D1u[r][(sk+1)/2],
                             D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                             D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,D2u,D2ucoef,D2uxcoef,
                             D2uxxcoef,D2jumpuxxcoef,jumpD1u,jumpD1ucoef,jumpD1uxcoef,
                             jumpD1uxxcoef,jumpD1jumpuxxcoef,S,pb,grid);


Input: all in terms of const, u, D1u, principle D2u, (mixed) jumpD2u

(1) approximation of Du- at interface in dim rstar direction sstar
Du-_r = D1u[r] + D1ucoef[r][:sindex] u_sindex +  D1uxcoef[r][:s] D1u_:s +  D1uxxcoef[r][:s] D2u_:ss +  D1jumpuxxcoef[r][:m][:n] jumpD2u_:m:n

(2) approximation of D2u- at grid point for all pairs of m,n
D2u_mn = D2u[m][n] + D2ucoef[m][n][:sindex] u_:sindex + D2uxcoef[m][n][:s] Du_:s + D2uxxcoef[m][n][:s] D2u_:ss + D2jumpuxxcoef[m][n] jumpD2u_mn

(3) approximation of jumpD1u at interface 
jumpD1u_rstar = jumpD1u + jumpD1ucoef[:sindex] u_:sindex + jumpD1uxcoef[:s] Du_:s + jumpD1uxxcoef[:s] D2u_:ss + jumpD1jumpuxxcoef[:m][:n] jumpD2u_:m:n
*/



void getcim345jumpuxx(double &u0, double ***ucoef, double *uxcoef, double *uxxcoef, 
                      int *index, int rstar, int sk, double alpha, int thesign, 
                      double *normal, int mid, double ***a, double *D1u, 
                      double ****D1ucoef, double **D1uxcoef, double **D1uxxcoef, 
                      double ***D1jumpuxxcoef, double **D2u, double *****D2ucoef, 
                      double ***D2uxcoef, double ***D2uxxcoef, double **D2jumpuxxcoef, 
                      double &jumpD1u, double ***jumpD1ucoef, double *jumpD1uxcoef, 
                      double *jumpD1uxxcoef, double **jumpD1jumpuxxcoef, double ***S, 
                      PBData &pb, GridData &grid)
{
   int i, r, s, n, m, N = 2*mid;
   int rindex[grid.dim], sindex[grid.dim];
   double b0[2*grid.dim], bxcoef[2*grid.dim][grid.dim], bxxcoef[2*grid.dim][grid.dim]; 
   double ***bcoef[2*grid.dim];
   double **LU, **B, **Dn; 
   double Dndott[grid.dim], Dndots[grid.dim], dotDndot[grid.dim];
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double tangent1[grid.dim], tangent2[grid.dim];
   double ethere, ehere, jumpfe;
   double value, temp[2*grid.dim];
   double sigma, tau, Dsigma[grid.dim];
   double Dtau[grid.dim], **D2tau, dotD2taudot[grid.dim], Dtaudot[grid.dim];
   int two2one[grid.dim][grid.dim];
   char theorder = 2;
   double aehere, aethere, jumpae;

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   B = matrix(2*grid.dim-1,2*grid.dim-1);
   Dn = matrix(grid.dim-1,grid.dim-1);
   D2tau = matrix(grid.dim-1,grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      bcoef[r] = matrix(N,N,N);

   if (evalarray(S,index) < 0.0)
   {
      ehere = pb.epsilonm;
      ethere = pb.epsilonp;
   }
   else
   {
      ehere = pb.epsilonp;
      ethere = pb.epsilonm;
   }

   getinterfaceinfo(tangent1,tangent2,tau,sigma,Dn,Dsigma,jumpfe,aehere,aethere,index,
                    rstar,sk,alpha,a,S,pb,grid);
   jumpae = thesign*(aehere-aethere); // jumpae defined as a+e+ - a-e-
// form Dn dot product with various vectors
   getDtau(Dtau,index,rstar,sk,alpha,grid);
   getD2tau(D2tau,index,rstar,sk,alpha,grid);


   for (n = 0; n < grid.dim; n++)
   {
      Dndott[n] = 0.0;
      Dndots[n] = 0.0;
      dotDndot[n] = 0.0;
      dotD2taudot[n] = 0.0;
      Dtaudot[n] = 0.0;
   }
   for (n = 0; n < grid.dim; n++)
   {
      for (m = 0; m < grid.dim; m++)
      {
         Dndott[n] += Dn[n][m]*tangent1[m]; //Dndott = Dn tangent1
         Dndots[n] += Dn[n][m]*tangent2[m]; //Dndots = Dn tangent2
         dotD2taudot[0] += tangent1[n]*D2tau[n][m]*tangent1[m]; //dotD2taudot[0] = tangent1 D2tau tangent1
         dotD2taudot[1] += tangent2[n]*D2tau[n][m]*tangent2[m]; //dotD2taudot[1] = tangent2 D2tau tangent2
         dotD2taudot[2] += tangent1[n]*D2tau[n][m]*tangent2[m]; //dotD2taudot[2] = tangent1 D2tau tangent2
      }
      dotDndot[0] += tangent1[n]*Dndott[n]; //dotDndot[0] = tangent1 Dn tangent1
      dotDndot[1] += tangent2[n]*Dndots[n]; //dotDndot[1] = tangent2 Dn tangent2
      dotDndot[2] += tangent1[n]*Dndots[n]; //dotDndot[2] = tangent1 Dn tangent2
      Dtaudot[0] += Dtau[n]*normal[n]; //Dtaudot[0] = Dtau n
      Dtaudot[1] += Dtau[n]*tangent1[n]; //Dtaudot[1] = Dtau tangent1
      Dtaudot[2] += Dtau[n]*tangent2[n]; //Dtaudot[2] = Dtau tangent2
   }

// form matrix
   for (n = 0; n < grid.dim; n++)
      two2one[n][n] = n;
   two2one[0][1] = grid.dim;
   two2one[1][2] = grid.dim+1;
   two2one[0][2] = grid.dim+2;
   for (m = 0; m < grid.dim; m++)
      for (n = 0; n < m; n++)
         two2one[m][n] = two2one[n][m];

   // === LHS of linear system
   // original LHS, only involve tangent and normal, not yet consdier d2u in RHS might include jumpuxx
   // first 3 row, s1[D2u]s1, s2[D2u]s2, s1[D2u]s2
   // next 2 row, s1[D2u]n, s2[D2u]n
   // last row, [lap(u)]
   for (n = 0; n < grid.dim; n++)
   {
      B[0][n] = tangent1[n]*tangent1[n];
      B[1][n] = tangent2[n]*tangent2[n];
      B[2][n] = tangent1[n]*tangent2[n];
      B[3][n] = normal[n]*tangent1[n];
      B[4][n] = normal[n]*tangent2[n];
      B[5][n] = 1.0;
   }
   for (n = grid.dim; n < 2*grid.dim; n++)
   {
      m = n-grid.dim;
      s = (m+1)%grid.dim;
      B[0][n] = 2.0*tangent1[m]*tangent1[s];
      B[1][n] = 2.0*tangent2[m]*tangent2[s];
      B[2][n] = tangent1[m]*tangent2[s]+tangent1[s]*tangent2[m];
      B[3][n] = normal[m]*tangent1[s]+normal[s]*tangent1[m];
      B[4][n] = normal[m]*tangent2[s]+normal[s]*tangent2[m];
      B[5][n] = 0.0;
   }

   // === contribution from Du and D2u (from RHS) to jump(LHS)
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         // get contributions of Du to jump data
         // first 3 rows, rhs include [eps]/epsm dot(Du-,n) sn.Dn.sm
         value = 0.0;
         for (r = 0; r < grid.dim; r++)
            value += D1jumpuxxcoef[r][m][n]*normal[r]; // dot(Du-,n)
         for (r = 0; r < grid.dim; r++)
            B[r][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*value*dotDndot[r];

         // in the 4th row, rhs include [eps]/epsm s1.Dn.Du-
         value = 0.0;
         for (r = 0; r < grid.dim; r++)
            value += D1jumpuxxcoef[r][m][n]*Dndott[r]; // Du- Dn tangent1
         B[grid.dim][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*value;
         
         // in the 5th row, rhs include [eps]/epsm s2.Dn.Du-
         value = 0.0;
         for (r = 0; r < grid.dim; r++)
            value += D1jumpuxxcoef[r][m][n]*Dndots[r]; // Du- Dn tangent2
         B[grid.dim+1][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*value;

         // in our application, we assume eps is piecewise constant, so no Deps term
         // no correction in the 6th row
         
         // get contributions of D2u to jump data
         // in the 4th row, rhs include [eps]/epsm s1.D2u-.n 
         B[grid.dim][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*
                                       (normal[m]*tangent1[n]+normal[n]*tangent1[m])*
                                       D2jumpuxxcoef[m][n];
         // in the 5th row, rhs include [eps]/epsm s1.D2u-.n 
         B[grid.dim+1][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*
                                         (normal[m]*tangent2[n]+normal[n]*tangent2[m])*
                                         D2jumpuxxcoef[m][n];
      }


   //=== constant term
   // in the first 3 row, sn.D2tau.sm + (sigma/epsm - Dtau.n )sn.Dn.sm
   for (n = 0; n < grid.dim; n++)
      b0[n] = (sigma/ethere-Dtaudot[0])*dotDndot[n]+dotD2taudot[n];

   b0[grid.dim] = getdotprod(Dsigma,tangent1,grid.dim)/ethere-
                  Dtaudot[1]*dotDndot[0]-Dtaudot[2]*dotDndot[2];
   b0[grid.dim+1] = getdotprod(Dsigma,tangent2,grid.dim)/ethere-
                    Dtaudot[1]*dotDndot[2]-Dtaudot[2]*dotDndot[1];
   // last row [f/eps] + ap/epsp tau
   b0[grid.dim+2] = -jumpfe+aethere*tau;

   // === The RHS include u-val, Du, D2u, with bcoef, bxcoef, bxxcoef
   // m = 1:6, s = 1:3
   // bcoef[m][sindex]: coef of u[sindx] in m-row
   // bxcoef[m][s]: coef of Du_s in m-row
   // bxxcoef[m][s]: coef of D2u_ss in m-row

   // === get b coefs
   // coef of u-value, b[r] = NxNxN matrix for r = 0,...5
   // all coeff of u-val in N-nbr of the r-th row of equation, 
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
// initialize to zero
      for (m = 0; m < 2*grid.dim; m++)
         setvalarray(bcoef[m],sindex,0.0);

      // get contributions from Du-
      // [eps]/epsp (Du- dot n) (si Dn sj)
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += evalarray(D1ucoef[n],sindex)*normal[n]; //Du- dot n
      for (m = 0; m < grid.dim; m++)
         setvalarray(bcoef[m],sindex,evalarray(bcoef[m],sindex)+
                                     thesign*(ethere-ehere)/ethere*value*dotDndot[m]);

       
      //[eps]/epsp Du- Dn s
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += evalarray(D1ucoef[n],sindex)*Dndott[n]; //Du- Dn s
      setvalarray(bcoef[grid.dim],sindex,evalarray(bcoef[grid.dim],sindex)+
                                  thesign*(ethere-ehere)/ethere*value);
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += evalarray(D1ucoef[n],sindex)*Dndots[n];
      setvalarray(bcoef[grid.dim+1],sindex,evalarray(bcoef[grid.dim+1],sindex)+
                                    thesign*(ethere-ehere)/ethere*value);

      // get contributions from D2u in row 4 and 5
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            setvalarray(bcoef[grid.dim],sindex,
                        evalarray(bcoef[grid.dim],sindex)+
                        thesign*(ethere-ehere)/ethere*
                        (normal[m]*tangent1[n]+normal[n]*tangent1[m])*
                         evalarray(D2ucoef[m][n],sindex));
            setvalarray(bcoef[grid.dim+1],sindex,
                        evalarray(bcoef[grid.dim+1],sindex)+
                        thesign*(ethere-ehere)/ethere*
                        (normal[m]*tangent2[n]+normal[n]*tangent2[m])*
                         evalarray(D2ucoef[m][n],sindex));
         }

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

   // get contributions from u- in the last row
   setvalarray(bcoef[grid.dim+2],sindex,jumpae);



   //=== get bx and bxx coefs
   for (s = 0; s < grid.dim; s++)
   {
      //=== get bx coefs
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxcoef[n][s]*normal[n];
      for (m = 0; m < grid.dim; m++)
         bxcoef[m][s] = thesign*(ethere-ehere)/ethere*value*dotDndot[m];

      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxcoef[n][s]*Dndott[n];
      bxcoef[grid.dim][s] = thesign*(ethere-ehere)/ethere*value;
      
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxcoef[n][s]*Dndots[n];
      bxcoef[grid.dim+1][s] = thesign*(ethere-ehere)/ethere*value;
      bxcoef[grid.dim+2][s] = 0.0;
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            bxcoef[grid.dim][s] += thesign*(ethere-ehere)/ethere*
                                   (normal[m]*tangent1[n]+normal[n]*tangent1[m])*
                                   D2uxcoef[m][n][s];
            bxcoef[grid.dim+1][s] += thesign*(ethere-ehere)/ethere*
                                     (normal[m]*tangent2[n]+normal[n]*tangent2[m])*
                                     D2uxcoef[m][n][s];
         }


      //=== get bxx coefs
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxxcoef[n][s]*normal[n];
      for (m = 0; m < grid.dim; m++)
         bxxcoef[m][s] = thesign*(ethere-ehere)/ethere*value*dotDndot[m];
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxxcoef[n][s]*Dndott[n];
      bxxcoef[grid.dim][s] = thesign*(ethere-ehere)/ethere*
                             (value+normal[s]*tangent1[s]);
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += D1uxxcoef[n][s]*Dndots[n];
      bxxcoef[grid.dim+1][s] = thesign*(ethere-ehere)/ethere*
                               (value+normal[s]*tangent2[s]);
      bxxcoef[grid.dim+2][s] = 0.0;
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            bxxcoef[grid.dim][s] += thesign*(ethere-ehere)/ethere*
                                    (normal[m]*tangent1[n]+normal[n]*tangent1[m])*
                                    D2uxxcoef[m][n][s];
            bxxcoef[grid.dim+1][s] += thesign*(ethere-ehere)/ethere*
                                      (normal[m]*tangent2[n]+normal[n]*tangent2[m])*
                                      D2uxxcoef[m][n][s];
         }
   }
// addition
   bxcoef[grid.dim+2][rstar] = jumpae*sk*alpha*grid.dx[rstar];
   bxxcoef[grid.dim+2][rstar] = jumpae*0.5*alpha*grid.dx[rstar]*alpha*grid.dx[rstar];//possible with minus sign?

   // === GE of LHS matrix
   gecp0(LU,PLR,PLC,B,2*grid.dim-1,2*grid.dim-1);


   // === for RHS cont, u-val, Du, D2u, solve the linear system,

   // === solve for RHS const
   // form jumpuxx in rstar direction and also recast Du and D2u
   // temp is constant term of jumpD2u[m][n]
   forwardbacksub0(temp,b0,LU,PLR,PLC,2*grid.dim-1);
   u0 = temp[two2one[rstar][rstar]];

   
   // D1u (Du-r) might depends on jumpD2u by D1jumpuxxcoef
   // recast constant term of jumpD2u in D1jumpuxxcoef to D1u
   for (s = 0; s < grid.dim; s++)
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
            D1u[s] += D1jumpuxxcoef[s][m][n]*temp[two2one[m][n]];
   
   // D2u might depends on jumpD2u by D2jumpuxxcoef
   // jumpD1u might depends on jumpD2u by jumpD1jumpuxxcoef
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         D2u[m][n] += D2jumpuxxcoef[m][n]*temp[two2one[m][n]];
         jumpD1u += jumpD1jumpuxxcoef[m][n]*temp[two2one[m][n]];
      }

   // === solve for RHS u-val
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      for (m = 0; m < 2*grid.dim; m++)
         temp[m] = evalarray(bcoef[m],sindex);
      // temp is coeff of u[sindex]
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      setvalarray(ucoef,sindex,temp[two2one[rstar][rstar]]);

      //D1u- depends on jumpD2u, which depends on u, so recast D1jumpuxxcoef to D1ucoef
      for (s = 0; s < grid.dim; s++)
         for (m = 0; m < grid.dim; m++)
            for (n = m+1; n < grid.dim; n++)
               setvalarray(D1ucoef[s],sindex,evalarray(D1ucoef[s],sindex)+
                                             D1jumpuxxcoef[s][m][n]*
                                             temp[two2one[m][n]]);

      //D2u- depends on jumpD2u, which depends on u, so recast D2jumpuxxcoef to D2ucoef
      //jumpD1u- depends on jumpD2u, which depends on u, so recast jumpD1jumpuxxcoef to jumpD1ucoef
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            setvalarray(D2ucoef[m][n],sindex,evalarray(D2ucoef[m][n],sindex)+
                                             D2jumpuxxcoef[m][n]*temp[two2one[m][n]]);
            setvalarray(jumpD1ucoef,sindex,evalarray(jumpD1ucoef,sindex)+
                                           jumpD1jumpuxxcoef[m][n]*temp[two2one[m][n]]);
         }

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }


   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

   // === solve for RHS Du- D2u-,  
   for (s = 0; s < grid.dim; s++)
   {
      for (m = 0; m < 2*grid.dim; m++)
         temp[m] = bxcoef[m][s];
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      uxcoef[s] = temp[two2one[rstar][rstar]];
      for (r = 0; r < grid.dim; r++)
         for (m = 0; m < grid.dim; m++)
            for (n = m+1; n < grid.dim; n++)
               D1uxcoef[r][s] += D1jumpuxxcoef[r][m][n]*temp[two2one[m][n]];
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            D2uxcoef[m][n][s] += D2jumpuxxcoef[m][n]*temp[two2one[m][n]];
            jumpD1uxcoef[s] += jumpD1jumpuxxcoef[m][n]*temp[two2one[m][n]];
         }
      for (m = 0; m < 2*grid.dim; m++)
         temp[m] = bxxcoef[m][s];
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      uxxcoef[s] = temp[two2one[rstar][rstar]];
      for (r = 0; r < grid.dim; r++)
         for (m = 0; m < grid.dim; m++)
            for (n = m+1; n < grid.dim; n++)
               D1uxxcoef[r][s] += D1jumpuxxcoef[r][m][n]*temp[two2one[m][n]];
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            D2uxxcoef[m][n][s] += D2jumpuxxcoef[m][n]*temp[two2one[m][n]];
            jumpD1uxxcoef[s] += jumpD1jumpuxxcoef[m][n]*temp[two2one[m][n]];
         }
   }
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         for (s = 0; s < grid.dim; s++)
            D1jumpuxxcoef[s][m][n] = 0.0;
         D2jumpuxxcoef[m][n] = 0.0;
         jumpD1jumpuxxcoef[m][n] = 0.0;
      }

   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(B,2*grid.dim-1,2*grid.dim-1);
   free_matrix(Dn,grid.dim-1,grid.dim-1);
   free_matrix(D2tau,grid.dim-1,grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(bcoef[r],N,N,N);
}




/*
after getting the coupling matrx
if  gamma[r][s] == 1, there is interface in dim r, s side

to reconstruct u- 
u- = u_index + (alpha[r][s] h) u_r + 1/2 (alpha[r][s] h) u_rr + O(h^3)
Du_r = ux[r] + uxcoef[:sindex] u_:sindex
D2u_rr = D2u[r][r] + D2ucoef[r][r][:sindex] u_:sindex

to reconstruct Du_t- at dim r direction s
after getcim345jumpuxx, we have
Du-_t = D1u[r][s][t] + D1ucoef[r][s][t][:sindex] u_:sindex + D1uxcoef[r][s][t][:s] Du_:s + D1uxxcoef[r][s][t][:s] D2u_:ss
backsub Du_:s and D2u_:ss by u-val

*/
void addtostorage( StorageStruct* &Dusmall, int &buildsize, 
   int mid, int *index,int gamma[][2], double**alpha, 
   double* ux, double**** uxcoef, 
   double*** D1u,double*** *** D1ucoef, double**** D1uxcoef, double**** D1uxxcoef,
   double** D2u, double***** D2ucoef, GridData &grid)
{
   int sindex[grid.dim];
   int rindex[grid.dim];
   int N = 2*mid;
   int Narray[grid.dim];

   for (int r = 0; r < grid.dim; r++)
      Narray[r] = N;

   double tol = 1.0e-14;
   for (int r = 0; r < grid.dim; r++)
      for (int sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
            newstorage(Dusmall[buildsize]);
            Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
            Dusmall[buildsize].info[1] = r;
            Dusmall[buildsize].info[2] = sk;
            Dusmall[buildsize].info[3] = -1;
            Dusmall[buildsize].mid = mid;
            double value = sk*alpha[r][(sk+1)/2]*grid.dx[r]*ux[r];
            if (globintorder == 3)
               value += 0.5*(alpha[r][(sk+1)/2]*grid.dx[r])*
                            (alpha[r][(sk+1)/2]*grid.dx[r])*D2u[r][r];
            sparseorder(-1,Dusmall[buildsize].head,value);
            
            for (int s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               int t;
               for (t = 0; t < grid.dim && sindex[t] == mid; t++);
               if (t >= grid.dim)
                  value = 1.0;
               else
                  value = 0.0;
               value += sk*alpha[r][(sk+1)/2]*grid.dx[r]*evalarray(uxcoef[r],sindex);
               if (globintorder == 3)
                  value += 0.5*(alpha[r][(sk+1)/2]*grid.dx[r])*
                               (alpha[r][(sk+1)/2]*grid.dx[r])*
                               evalarray(D2ucoef[r][r],sindex);
               if (fabs(value) > tol)
               {
                  if (globsmall == 1)
                     sparseorder(sub2ind(sindex,Narray,grid.dim),Dusmall[buildsize].head,
                                 value);
                  else if (globsmall == 2)
                  {
                     for (int s = 0; s < grid.dim; s++)
                        rindex[s] = index[s]+sindex[s]-mid;
                     sparseorder(sub2ind(rindex,grid.nx,grid.dim),
                                 Dusmall[buildsize].head,value);
                  }
               }
   
               (sindex[grid.dim-1])++;
               for (int i = grid.dim-1; i > 0 && sindex[i] > N; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (int s = 0; s < grid.dim; s++)
               sindex[s] = mid;

            buildsize++;

// storing Du info
            for (int t = 0; t < grid.dim; t++)
            {
               newstorage(Dusmall[buildsize]);
               Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
               Dusmall[buildsize].info[1] = r;
               Dusmall[buildsize].info[2] = sk;
               Dusmall[buildsize].info[3] = t; // which component of gradient
               Dusmall[buildsize].mid = mid;

               // constant term
               value = D1u[r][(sk+1)/2][t];
               for (int m = 0; m < grid.dim; m++)
                  value += D1uxcoef[r][(sk+1)/2][t][m]*ux[m]+
                           D1uxxcoef[r][(sk+1)/2][t][m]*D2u[m][m];
               if (fabs(value) > tol)
                  sparseorder(-1,Dusmall[buildsize].head,value);

               // u-val term
               for (int s = 0; s < grid.dim; s++)
                  sindex[s] = 0;
               while (sindex[0] <= N)
               {
                  value = evalarray(D1ucoef[r][(sk+1)/2][t],sindex);
                  for (int m = 0; m < grid.dim; m++)
                     value += D1uxcoef[r][(sk+1)/2][t][m]*evalarray(uxcoef[m],sindex)+
                              D1uxxcoef[r][(sk+1)/2][t][m]*
                              evalarray(D2ucoef[m][m],sindex);
                  if (fabs(value) > tol)
                  {
                     if (globsmall == 1)
                        sparseorder(sub2ind(sindex,Narray,grid.dim),
                                    Dusmall[buildsize].head,value);
                     else if (globsmall == 2)
                     {
                        for (int s = 0; s < grid.dim; s++)
                           rindex[s] = index[s]+sindex[s]-mid;
                        sparseorder(sub2ind(rindex,grid.nx,grid.dim),
                                    Dusmall[buildsize].head,value);
                     }
                  }

                  (sindex[grid.dim-1])++;
                  for (int i = grid.dim-1; i > 0 && sindex[i] > N; i--)
                  {
                     sindex[i] = 0;
                     (sindex[i-1])++;
                  }
               }
               for (int s = 0; s < grid.dim; s++)
                  sindex[s] = mid;

               buildsize++;
            }
         }
}

/*
approximate (mixed) D2u at grid point in terms of const, u-val, Du, principle D2u, (mixed) jumpD2u
D2u_mn = u0 + u[sindex] u_sindex + uxcoef[:s] Du_s + uxxcoef[:s] D2u_ss + jumpuxxcoef jumpD2u_mn

acutally u0 always 0

used when with globdist = 1, otherwise use [getcim345D2u]

in the first pass, perm=0, rstart = sstar = 0, 
return perm=2 if there is usual mix Du, then second pass is skipped
return perm=1 if have to use the other side, then in the second pass, location for interface is provided,  use rstar and sstar to find jumpuxxcoeff

globdistvar = 1: 
getsk2 with central diff (includeing rectangle) > cim5 (triangle stencil) > getsk2 without central diff > sk2 on out-of-plane nbr point(prefer the side with central diff) > sk2 on the other side
globdistvar = 0: this is the same as [getcim345D2u]
getsk2 (with or without central diff) > cim5 (triangle stencil) > sk2 on out-of-plane point(prefer the side with central diff) > sk2 on the other side
*/

void getcim345D2udist(double &u0, double ***u, double *uxcoef, double *uxxcoef, 
                      double &jumpuxxcoef, char &perm, int m, int n, int *index, 
                      int rstar, int sstar, int mid, double ***S, GridData &grid)
{
   if (perm != 2)
   {
      int r, s, t, N = 2*mid, sindex[grid.dim], rindex[grid.dim], nindex[grid.dim]; 
      int sk2reg[4], sk2other[4], **sk2near = imatrix(1,3);  //sk2near is 2 x 4 matrix
      char regD2, useuxD2, otherD2, nearD2[2];
      double thesign, tol = 1.0e-14;

      for (t = 0; t < grid.dim; t++)
         rindex[t] = index[t];
      if (perm != 1)
      {
         regD2 = getsk2(sk2reg,m,n,index,S,grid);
         useuxD2 = getcim5D2(u,uxcoef,uxxcoef,m,n,index,mid,S,grid);
         for (r = 0; r < grid.dim && (r == m || r == n); r++);// find out of plane dimension
         for (s = -1; s <= 1; s += 2)
         {
            rindex[r] = index[r]+s;
            if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) != 1)
               nearD2[(s+1)/2] = getsk2(sk2near[(s+1)/2],m,n,rindex,S,grid);//nearD2 is approx of D2 of nbr point in same side, out of plane
            else
               nearD2[(s+1)/2] = 0;
         }
         rindex[r] = index[r];
      }
      rindex[rstar] = index[rstar]+sstar; // get the other side
      otherD2 = getsk2(sk2other,m,n,rindex,S,grid); 
      rindex[rstar] = index[rstar];
/*
in either case below , use sk2reg, set perm=2
case 1: perm != 1 && !globdistvar && regD2 && ( [1] || [2] )
   where [1] abs(sk2reg[1]-sk2reg[0])+abs(sk2reg[3]-sk2reg[2]) >= 3 means one direction of cross derivative use central difference
         [2] (useuxD2 && (fabs(uxxcoef[m]) > tol || fabs(uxxcoef[n]) > tol) means cim5D2 use uxxcoef in m n plane
case 2: perm != 1 &&  globdistvar && regD2 &&  [1]
if globdistvar = 0, we choose sk2reg even when sk2reg is not central diff,
if globdistvar = 1, we choose sk2reg only when sk2reg is central diff     
*/
      if ((perm != 1 && !globdistvar && regD2 && (abs(sk2reg[1]-sk2reg[0])+abs(sk2reg[3]-sk2reg[2]) >= 3 ||(useuxD2 && (fabs(uxxcoef[m]) > tol || fabs(uxxcoef[n]) > tol))) ) ||
          (perm != 1 && globdistvar && regD2 && abs(sk2reg[1]-sk2reg[0])+abs(sk2reg[3]-sk2reg[2]) >= 3))
      {
         for (t = 0; t < grid.dim; t++)
         {
            sindex[t] = mid;
            nindex[t] = N;
         }
         getD2(u,m,n,sk2reg,sindex,nindex,grid);
         u0 = 0.0;
         for (t = 0; t < grid.dim; t++)
         {
            uxcoef[t] = 0.0;
            uxxcoef[t] = 0.0;
         }
         jumpuxxcoef = 0.0;
         perm = 2;
      }
      else if (perm != 1 && useuxD2) // if (regD2 failed) or (regD2 but no central diff) and (has getcim5D2), then use cim5
      {
         u0 = 0.0;
         jumpuxxcoef = 0.0;
         perm = 2;
      }
      else if (perm != 1 && (nearD2[0] || nearD2[1])) //if (regD2 and getcim5D2 failed), use out-of-plane points, prefer central diff
      {
         if ((nearD2[0] && !(nearD2[1])) || (nearD2[0] && nearD2[1] && 
            abs(sk2near[0][0]-sk2near[0][1])+abs(sk2near[0][2]-sk2near[0][3]) >abs(sk2near[1][0]-sk2near[1][1])+abs(sk2near[1][2]-sk2near[1][3])))
            s = -1;
         else
            s = 1;
         for (t = 0; t < grid.dim; t++)
         {
            sindex[t] = mid;
            nindex[t] = N;
         }
         sindex[r] = mid+s;
         getD2(u,m,n,sk2near[(s+1)/2],sindex,nindex,grid);
         u0 = 0.0;
         for (t = 0; t < grid.dim; t++)
         {
            uxcoef[t] = 0.0;
            uxxcoef[t] = 0.0;
         }
         jumpuxxcoef = 0.0;
         perm = 2;
      }
      // every method that use same side fails. Use the other side, set perm = 1
      // this is skiped in the first pass, sstar = 0
      // entered in the second pass with sstar !=1, know interface location, use the other side
      else if (rstar >= 0 && rstar < grid.dim && sstar != 0 && otherD2) 
      {
         if (evalarray(S,index) < 0.0)
            thesign = -1;
         else
            thesign = 1;
         for (t = 0; t < grid.dim; t++)
         {
            sindex[t] = mid;
            nindex[t] = N;
         }
         sindex[rstar] = mid+sstar;
         getD2(u,m,n,sk2other,sindex,nindex,grid);
         jumpuxxcoef = thesign;
         u0 = 0.0;
         for (t = 0; t < grid.dim; t++)
         {
            uxcoef[t] = 0.0;
            uxxcoef[t] = 0.0;
         }
         perm = 1;
      }
      // also skipped in the first pass as sstar=0
      else if (rstar >= 0 && rstar < grid.dim && sstar != 0)
      {
         cout << "bad status" << endl;
         perm = 0;
      }
      else
         perm = 1;

      free_matrix(sk2near,1,3); 
   }
}


// calculate Du- at interface between index and index[rstar]+start
// full process: construct and solve coupling equation.
void getcim345Du(double &uint, double *Du, int *index, int rstar, int sstar,
                 double ***u, double ***S, PBData &pb, GridData &grid)
{
   int r, s, t, i, m, n, sk, mid = 2, N = 2*mid;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim];
   double d0[2*grid.dim];
   double ***dcoef[2*grid.dim];
   double **LU, **G;
   int gamma[grid.dim][2];
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double alpha, beta, normal[grid.dim], tangent[grid.dim];
   double ethere, ehere;
   double temp[2*grid.dim];
   double sigma, tau, Dsigma[grid.dim], Dtau[grid.dim], thesign;

   double D1u[grid.dim], ****D1ucoef, **D1uxcoef, **D1uxxcoef, 
          ***D1jumpuxxcoef;
   D1ucoef = new double ***[grid.dim];
   for (r = 0; r < grid.dim; r++)
      D1ucoef[r] = matrix(N,N,N);
   D1uxcoef = matrix(grid.dim-1,grid.dim-1);
   D1uxxcoef = matrix(grid.dim-1,grid.dim-1);
   D1jumpuxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);

   double finalalpha, finalD1u[grid.dim], ****finalD1ucoef, **finalD1uxcoef, 
          **finalD1uxxcoef;
   finalD1ucoef = new double ***[grid.dim];
   for (r = 0; r < grid.dim; r++)
      finalD1ucoef[r] = matrix(N,N,N);
   finalD1uxcoef = matrix(grid.dim-1,grid.dim-1);
   finalD1uxxcoef = matrix(grid.dim-1,grid.dim-1);

   double **D2u, *****D2ucoef, ***D2uxcoef, ***D2uxxcoef, **D2jumpuxxcoef;
   D2u = matrix(grid.dim-1,grid.dim-1);
   D2ucoef = new double ****[grid.dim];
   for (r = 0; r < grid.dim; r++)
   {
      D2ucoef[r] = new double ***[grid.dim];
      for (s = 0; s < grid.dim; s++)
         D2ucoef[r][s] = matrix(N,N,N);
   }
   D2uxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   D2uxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   D2jumpuxxcoef = matrix(grid.dim-1,grid.dim-1);

   double jumpD1u = 0.0, ***jumpD1ucoef, jumpD1uxcoef[grid.dim], 
          jumpD1uxxcoef[grid.dim], **jumpD1jumpuxxcoef;
   jumpD1ucoef = matrix(N,N,N);
   jumpD1jumpuxxcoef = matrix(grid.dim-1,grid.dim-1);

   double jumpD2u = 0.0, ***jumpD2ucoef, jumpD2uxcoef[grid.dim], 
          jumpD2uxxcoef[grid.dim];
   jumpD2ucoef = matrix(N,N,N);

   char yesD2[grid.dim][grid.dim];


   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
   {
      for (s = -1; s <= 1; s += 2)
      {
         rindex[r] = min(max(index[r]+s,0),grid.nx[r]);
         if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1)
            gamma[r][(s+1)/2] = 1;
         else
            gamma[r][(s+1)/2] = 0;//gamma[r][s=0] = dim r, s=-1, has interface
      }
      rindex[r] = index[r];
   }

   LU = matrix(2*grid.dim-1,2*grid.dim-1);
   G = matrix(2*grid.dim-1,2*grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      dcoef[r] = matrix(N,N,N);

   if (evalarray(S,index) < 0.0)
   {
      ehere = pb.epsilonm;
      ethere = pb.epsilonp;
      thesign = -1.0;
   }
   else
   {
      ehere = pb.epsilonp;
      ethere = pb.epsilonm;
      thesign = 1.0;
   }
   // get approximation of D2u at index
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         yesD2[m][n] = 0;
         if (!globdist)
            getcim345D2u(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],
                         D2jumpuxxcoef[m][n],yesD2[m][n],m,n,index,0,0,mid,S,grid);
         else
            getcim345D2udist(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],
                             D2jumpuxxcoef[m][n],yesD2[m][n],m,n,index,0,0,mid,S,grid);
      }

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
//            cout << "r = " << r << " " << sk << endl;
//            cout << "getting D2u" << endl;
            for (m = 0; m < grid.dim; m++)
               for (n = m+1; n < grid.dim; n++)
                  if (!globdist)
                     getcim345D2u(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],
                                  D2uxxcoef[m][n],D2jumpuxxcoef[m][n],yesD2[m][n],
                                  m,n,index,r,sk,mid,S,grid);
                  else
                     getcim345D2udist(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],
                                      D2uxxcoef[m][n],D2jumpuxxcoef[m][n],yesD2[m][n],
                                      m,n,index,r,sk,mid,S,grid);
            
            rindex[r] = index[r]+sk;
// getting interface info
            getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
            beta = 1.0-alpha;
// get derivatives
//            cout << "getting Du" << endl;
            getcim345Du(D1u,D1ucoef,D1uxcoef,D1uxxcoef,D1jumpuxxcoef,index,r,sk,alpha,
                        thesign,D2u,D2ucoef,D2uxcoef,D2uxxcoef,D2jumpuxxcoef,mid,grid);
// getting jump of uxx in r
            gettau(tau,index,r,sk,alpha,grid);
//            cout << "getting jump Du" << endl;
            getcim345jumpux(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,
                            jumpD1jumpuxxcoef,index,r,sk,alpha,thesign,normal,tangent,
                            mid,D1ucoef,D1uxcoef,D1uxxcoef,D1jumpuxxcoef,S,pb,grid);
//            cout << "getting jump D2u" << endl;
            getcim345jumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,
                             alpha,thesign,normal,mid,D1u,D1ucoef,D1uxcoef,D1uxxcoef,
                             D1jumpuxxcoef,D2u,D2ucoef,D2uxcoef,D2uxxcoef,D2jumpuxxcoef,
                             jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,
                             jumpD1jumpuxxcoef,S,pb,grid);
            if (r == rstar && sk == sstar)
            {
               finalalpha = alpha;

               for (s = 0; s < grid.dim; s++)
                  finalD1u[s] = D1u[s];

               for (s = 0; s < grid.dim; s++)
                  sindex[s] = 0;
               while (sindex[0] <= N)
               {
                  for (s = 0; s < grid.dim; s++)
                     setvalarray(finalD1ucoef[s],sindex,evalarray(D1ucoef[s],sindex));

                  (sindex[grid.dim-1])++;
                  for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
                  {
                     sindex[i] = 0;
                     (sindex[i-1])++;
                  }
               }
               for (s = 0; s < grid.dim; s++)
                  sindex[s] = mid;

               for (s = 0; s < grid.dim; s++)
                  for (t = 0; t < grid.dim; t++)
                  {
                     finalD1uxcoef[s][t] = D1uxcoef[s][t];
                     finalD1uxxcoef[s][t] = D1uxxcoef[s][t];
                  }
            }
// form d0 and dcoef's rth entry 
            d0[grid.dim*(sk+1)/2+r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+
                                               sk*beta*jumpD1u/grid.dx[r]+
                                               0.5*beta*beta*jumpD2u)-
                                      sk*D1u[r]/grid.dx[r];

            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                           thesign*(sk*beta*evalarray(jumpD1ucoef,sindex)/grid.dx[r]+
                                    0.5*beta*beta*evalarray(jumpD2ucoef,sindex))-
                           sk*evalarray(D1ucoef[r],sindex)/grid.dx[r]);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = mid;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                        evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)-
                        1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = mid+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                        evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)+
                        1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = mid;
// form matrix's rth row
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m] = -thesign*
                                           (sk*beta*jumpD1uxxcoef[m]/grid.dx[r]+
                                            0.5*beta*beta*jumpD2uxxcoef[m])+
                                           sk*D1uxxcoef[r][m]/grid.dx[r];
            G[grid.dim*(sk+1)/2+r][r] += 0.5*(beta*beta-alpha*alpha);
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m+grid.dim] = -thesign*
                                                    (sk*beta*jumpD1uxcoef[m]/grid.dx[r]+
                                                     0.5*beta*beta*jumpD2uxcoef[m])+
                                                    sk*D1uxcoef[r][m]/grid.dx[r];
            rindex[r] = index[r];
         }
         else
         {
// if interior point
            d0[grid.dim*(sk+1)/2+r] = 0.0;
            for (s = 0; s < 2*grid.dim; s++)
               G[grid.dim*(sk+1)/2+r][s] = 0.0;
            G[grid.dim*(sk+1)/2+r][r] = 0.5;
            G[grid.dim*(sk+1)/2+r][grid.dim+r] = sk/grid.dx[r];
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,0.0);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = mid;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,-1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = mid+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = mid;
         }

// solve for uxx and put in matrix
   int j;
   double uxx[grid.dim], ux[grid.dim];

   gecp0(LU,PLR,PLC,G,2*grid.dim-1,2*grid.dim-1);
   forwardbacksub0(temp,d0,LU,PLR,PLC,2*grid.dim-1);
   for (j = 0; j < grid.dim; j++)
      uxx[j] = temp[j];
   for (j = 0; j < grid.dim; j++)
      ux[j] = temp[j+grid.dim];

   for (s = 0; s < grid.dim; s++)
      tindex[s] = index[s];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      for (n = 0; n < 2*grid.dim; n++)
         temp[n] = evalarray(dcoef[n],sindex);
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-mid+sindex[s];
      for (j = 0; j < grid.dim; j++)
      {
         uxx[j] += temp[j]*evalarray(u,tindex);
         ux[j] += temp[j+grid.dim]*evalarray(u,tindex);
      }

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

   if (globintorder == 3)
      uint = evalarray(u,index)+sstar*finalalpha*grid.dx[rstar]*ux[rstar]+
             0.5*(finalalpha*grid.dx[rstar])*(finalalpha*grid.dx[rstar])*uxx[rstar];
   else
      uint = evalarray(u,index)+sstar*finalalpha*grid.dx[rstar]*ux[rstar];

   for (s = 0; s < grid.dim; s++)
      Du[s] = finalD1u[s];

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-mid+sindex[s];
      for (s = 0; s < grid.dim; s++)
         Du[s] += evalarray(finalD1ucoef[s],sindex)*evalarray(u,tindex);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

   for (s = 0; s < grid.dim; s++)
      for (t = 0; t < grid.dim; t++)
         Du[s] += finalD1uxcoef[s][t]*ux[t];

   for (s = 0; s < grid.dim; s++)
      for (t = 0; t < grid.dim; t++)
         Du[s] += finalD1uxxcoef[s][t]*uxx[t];

/*
   char yescim4 = 0;
   for (m = 0; m < grid.dim && !yescim4; m++)
      for (n = m+1; n < grid.dim && !yescim4; n++)
         if (!(yesD2[m][n]))
            yescim4 = 1;
   if (yescim4)
   {
      rindex[rstar] = index[rstar]+sstar;
      getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
      cout << "CIM 4 Du Error at " << index[0] << " " << index[1] << " " << index[2] 
           << " in direction " << rstar << " " << sstar << ": " 
           << Du[0]-getDu(index,0,rstar,sstar,alpha,thesign,grid) << " "
           << Du[1]-getDu(index,1,rstar,sstar,alpha,thesign,grid) << " "
           << Du[2]-getDu(index,2,rstar,sstar,alpha,thesign,grid) << endl;
      cout << "   official status is " << getstatus5(S,index,grid) << endl;
      cout << "   ";
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
            cout << (int) yesD2[m][n] << " ";
      cout << endl;
      rindex[rstar] = index[rstar];
   }
*/
/*
   if (index[0] == 34 && index[1] == 31 && index[2] == 56)
   {
      cout << rstar << " " << sstar << endl;
      for (m = 0; m < grid.dim; m++)
         rindex[m] = index[m];
      rindex[rstar] = index[rstar]+sstar;
      getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
      cout << Du[0] << " " << Du[1] << " " << Du[2] << endl;
      cout << getDu(index,0,rstar,sstar,alpha,thesign,grid) << " "
           << getDu(index,1,rstar,sstar,alpha,thesign,grid) << " "
           << getDu(index,2,rstar,sstar,alpha,thesign,grid) << endl;
      cout << endl;
      cout << ux[0] << " " << ux[1] << " " << ux[2] << endl;
      cout << getDu(index,0,0,0,0.0,thesign,grid) << " "
           << getDu(index,1,0,0,0.0,thesign,grid) << " "
           << getDu(index,2,0,0,0.0,thesign,grid) << endl;
      cout << endl;
      cout << uxx[0] << " " << uxx[1] << " " << uxx[2] << endl;
      cout << getD2u(index,0,0,0,0,0.0,thesign,grid) << " "
           << getD2u(index,1,1,0,0,0.0,thesign,grid) << " "
           << getD2u(index,2,2,0,0,0.0,thesign,grid) << endl;
      cout << endl;
   }
*/

   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(G,2*grid.dim-1,2*grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(dcoef[r],N,N,N);

   for (r = 0; r < grid.dim; r++)
      free_matrix(D1ucoef[r],N,N,N);
   delete [] D1ucoef;
   free_matrix(D1uxcoef,grid.dim-1,grid.dim-1);
   free_matrix(D1uxxcoef,grid.dim-1,grid.dim-1);
   free_matrix(D1jumpuxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);

   for (r = 0; r < grid.dim; r++)
      free_matrix(finalD1ucoef[r],N,N,N);
   delete [] finalD1ucoef;
   free_matrix(finalD1uxcoef,grid.dim-1,grid.dim-1);
   free_matrix(finalD1uxxcoef,grid.dim-1,grid.dim-1);

   free_matrix(D2u,grid.dim-1,grid.dim-1);
   for (r = 0; r < grid.dim; r++)
   {
      for (s = 0; s < grid.dim; s++)
         free_matrix(D2ucoef[r][s],N,N,N);
      delete [] D2ucoef[r];
   }
   delete [] D2ucoef;
   free_matrix(D2uxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(D2uxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(D2jumpuxxcoef,grid.dim-1,grid.dim-1);

   free_matrix(jumpD1ucoef,N,N,N);
   free_matrix(jumpD1jumpuxxcoef,grid.dim-1,grid.dim-1);

   free_matrix(jumpD2ucoef,N,N,N);
}


/*

Notation: :s means summation over s, Du_s means partial derivative

Input:
called after getcim345D2u, in which we obtained
D2u_mn = D2u[m][n] + D2ucoef[m][n][sindex] u_sindex + D2uxcoef[m][n][:s] Du_s + D2uxxcoef[m][n][:s] D2u_ss + D2jumpuxxcoef[m][n] jumpD2u_mn

Output:
get Du-, i.e. Du at interface dim rstar direction star in terms of (const, u-val, Du, principal D2u, jumpD2u (can have mixed))
component r
Du-_r = u0[r] +  ucoef[r][:index] u_:index +  uxcoef[r][:s] Du_:s + uxxoef[r][:s] D2u_:ss + jumpuxxcoef[r][:m][:n] jumpD2u_:m:n

Method:
start with taylor, t-component of Du-, interface location in dim rstar direction star
Du-_t = Du_t + alpha h D2u_rstar,t
then substitute D2u_rstar,t

actually u0 always 0, because D2u always 0 (no constant term involved when approximating D2u)

called with
getcim345Du(D1u[r][(sk+1)/2],D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
            D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,index,r,sk,
            alpha[r][(sk+1)/2],thesign,D2u,D2ucoef,D2uxcoef,D2uxxcoef,
            D2jumpuxxcoef,mid,grid);
*/

void getcim345Du(double *u0, double ****ucoef, double **uxcoef, double **uxxcoef,
                 double ***jumpuxxcoef, int *index, int rstar, int sstar, double alpha,
                 double thesign, double **D2u, double *****D2ucoef, double ***D2uxcoef, 
                 double ***D2uxxcoef, double **D2jumpuxxcoef, int mid, GridData &grid)
{
// getDu and recast
   int i, r, s, t, m, n, N = 2*mid, sindex[grid.dim];
   double ***uxxcoeflarge = matrix(grid.dim-1,grid.dim-1,grid.dim-1);//uxxcoeflarg[r][i][j] = Du-_r in terms of DDu_ij

   // set 0 uxxcoeflarge jumpuxxcoef u0 uxcoef
   for (r = 0; r < grid.dim; r++)
   {
      u0[r] = 0.0;
      for (s = 0; s < grid.dim; s++)
      {
         uxcoef[r][s] = 0.0;
         for (t = 0; t < grid.dim; t++)
         {
            uxxcoeflarge[r][s][t] = 0.0;
            jumpuxxcoef[r][s][t] = 0.0;
         }
      }
   }

   // set 0: ucoef
   for (i = 0; i < grid.dim; i++)
      sindex[i] = 0;
   while (sindex[0] <= N)
   {
      for (r = 0; r < grid.dim; r++)
         setvalarray(ucoef[r],sindex,0.0);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }

   // recast 
   // first order taylor Du-[s] = Du[s] + sstar  alpha h D2u[s][rstar]
   // so uxxcoef[s][s] = 1, uxxcoeflarge[s][s][rstar] = sstar alpha h
   // 
   for (s = 0; s < grid.dim; s++)
   {
      uxcoef[s][s] = 1.0;
      uxxcoeflarge[s][min(rstar,s)][max(rstar,s)] = sstar*alpha*grid.dx[rstar];
   }
   // recast uxxcoeflarge, 
   for (r = 0; r < grid.dim; r++)
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
         {
            u0[r] += uxxcoeflarge[r][m][n]*D2u[m][n];

            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               setvalarray(ucoef[r],sindex,evalarray(ucoef[r],sindex)+
                                           uxxcoeflarge[r][m][n]*evalarray(D2ucoef[m][n],sindex)); 
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }

            for (s = 0; s < grid.dim; s++)
               uxcoef[r][s] += uxxcoeflarge[r][m][n]*D2uxcoef[m][n][s]; 
            for (s = 0; s < grid.dim; s++)
               uxxcoeflarge[r][s][s] += uxxcoeflarge[r][m][n]*D2uxxcoef[m][n][s]; 
   
            jumpuxxcoef[r][m][n] += uxxcoeflarge[r][m][n]*D2jumpuxxcoef[m][n];
   
            uxxcoeflarge[r][m][n] = 0.0;
            uxxcoeflarge[r][n][m] = 0.0;
         }

   for (r = 0; r < grid.dim; r++)
      for (s = 0; s < grid.dim; s++)
         uxxcoef[r][s] = uxxcoeflarge[r][s][s];

   free_matrix(uxxcoeflarge,grid.dim-1,grid.dim-1,grid.dim-1);
}


// with a, use dusmall
void linearsystemcim345(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, 
                   int &smallsize, double ***a, double ***S, PBData &pb, 
                   GridData &grid)
{
// can use cim1, cim2, cim3, cim4, cim5 and saves on speed in calculating Du
   int gamma[grid.dim][2];
   int i, r, s, j, k, thestatus;
   int tindex[grid.dim], rindex[grid.dim], numcount = 7, count[numcount];
   double Sval1, Sval2;
   double x[grid.dim];
   char except;
   double temp;
   int tempsk2[4];

   cout << "in linearsystemcim345" << endl;

   clearsparse(A,grid);

   int buildsize = 0;
   smallsize = getstoragesize(S,grid);
   Dusmall = new StorageStruct[smallsize];
   cout << "smallsize = " << smallsize << endl;
   
   for (i = 0; i < numcount; i++)
      count[i] = 0;

   // Set b vector by getf
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      if (evalarray(S,tindex) < 0.0)
         setvalarray(b,tindex,getf(tindex,0,0,0.0,-1,pb,grid));
      else
         setvalarray(b,tindex,getf(tindex,0,0,0.0,1,pb,grid));

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   // Assemble A matrix
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         for (s = -1; s <= 1; s += 2)
         {
            rindex[r] = min(max(tindex[r]+s,0),grid.nx[r]);
            if (evalarray(S,tindex)*evalarray(S,rindex) < 0.0)
               gamma[r][(s+1)/2] = 1;
            else
               gamma[r][(s+1)/2] = 0;
         }
         rindex[r] = tindex[r];
      }
		thestatus = getstatus5(S,tindex,grid);

      if (thestatus == 1)
      {
         // interiorptsmall(A,b,tindex,S,pb,grid); // temporary: interiorptsmall should do nothing
         (count[thestatus])++;
      }
      else if (thestatus >= 2 && thestatus <=4)
      {        
         cim345cond(A,b,Dusmall,buildsize,tindex,a,gamma,S,pb,grid);  
         (count[thestatus])++;
      } 
      else if (thestatus == 5)
      {
         // cim1(A,b,Dusmall,buildsize,tindex,a,gamma,S,pb,grid);
         cout<<"cim1 not yet separated, not suppose to be here"<<endl;
         exit(1);
         (count[thestatus])++;
      }
      else if (thestatus == 0)
      {
         sparse2(tindex,tindex,A,1.0,grid);
         for (r = 0; r < grid.dim; r++)
            x[r] = grid.a[r]+tindex[r]*grid.dx[r];
         setvalarray(b,tindex,DBC(x,grid.dim,0.0));
         (count[thestatus])++;
      }
      else
         (count[6])++;

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   int thecount = 0, thecount2 = 0.0;
   SparseElt2 *current;
   cout << sizeof(SparseElt2) << " " << sizeof(int *) << " " << sizeof(double *) << " "
        << sizeof(int) << " " << sizeof(double) << endl;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      for (current = evalarray(A,tindex); current != NULL; current = (*current).next)
      {
          thecount++;
      }
      
      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
   cout << "found " << count[0] << " boundary pts, " << count[1] << " interior pts, ";
   cout << count[2] << " cim3 pts, " << count[3] << " cim5 pts," << endl
           << "   and " << count[4] << " cim4 pts, " << count[5] << " cim1 pts." 
           << endl;


// Dusmall will hold more memory: smallsize; but actually use less: buildsize
   smallsize = buildsize;
   cout << "smallsize reduced to " << smallsize << endl;

   
   thecount2 = 0;
   for (i = 0; i < smallsize; i++)
      for (SparseElt *current = Dusmall[i].head; current != NULL; 
           current = (*current).next)
         thecount2++;
   cout << "Du storage structure has " << thecount2 << " elements." << endl;

}





// to use without a(x), pass ***a as nullptr, only difference is getcim345jumpuxx
void gmatrix(double **G, double** LU, int PLR[], int PLC[],
  double d0[], double**** dcoef,
   double*** D1u,double*** *** D1ucoef, double**** D1uxcoef, double**** D1uxxcoef,
   double** D2u, double***** D2ucoef, 
   int *index,
   double** alpha, int gamma[][2], double ***S, double***a, PBData &pb, GridData &grid){
   


   int r, s, t, i, j, m, n, sk, mid = 2, N = 2*mid;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim], Narray[grid.dim];
   double value;
   double beta, normal[grid.dim], tangent[grid.dim];
   double ethere, ehere;
   double temp[2*grid.dim];
   double sigma, tau, Dsigma[grid.dim], Dtau[grid.dim], thesign;

   double exactd[2*grid.dim], exactres, tempsign;
   int tempindex[grid.dim];

   for (r = 0; r < grid.dim; r++)
      Narray[r] = N;

   
   double ***D1jumpuxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   
   

   double ***D2uxcoef, ***D2uxxcoef, **D2jumpuxxcoef;
   D2uxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   D2uxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   D2jumpuxxcoef = matrix(grid.dim-1,grid.dim-1);

   double jumpD1u = 0.0, ***jumpD1ucoef, jumpD1uxcoef[grid.dim], 
          jumpD1uxxcoef[grid.dim], **jumpD1jumpuxxcoef;
   jumpD1ucoef = matrix(N,N,N);
   jumpD1jumpuxxcoef = matrix(grid.dim-1,grid.dim-1);

   double jumpD2u = 0.0, ***jumpD2ucoef, jumpD2uxcoef[grid.dim], 
          jumpD2uxxcoef[grid.dim];
   jumpD2ucoef = matrix(N,N,N);

   char yesD2[grid.dim][grid.dim];
   
   // print geometry
   if (equal(index,eindex))
   {
      cout << "S" << endl;
      print_surf(S, index, 1);
      cout << "u" << endl;
      for (t = -1; t <= 1; t++) 
      {
         for (s = 1; s >= -1; s--) 
         {
            for (r = -1; r <= 1; r++) 
            {
               tindex[0] = index[0]+r;
               tindex[1] = index[1]+s;
               tindex[2] = index[2]+t;
               thesign = (evalarray(S,index) < 0.0)?-1:1;
               printf("%4.16f ",getu(tindex,0,0,0.0,thesign,grid));
            }
            cout << endl;
         }
         cout << endl;
      }
   }

   value = 0.0;

   
   for (r = 0; r < 2*grid.dim; r++)
      dcoef[r] = matrix(N,N,N);


   // actual calculation
   if (evalarray(S,index) < 0.0)
   {
      ehere = pb.epsilonm;
      ethere = pb.epsilonp;
      thesign = -1.0;
   }
   else
   {
      ehere = pb.epsilonp;
      ethere = pb.epsilonm;
      thesign = 1.0;
   }

   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
// yesD2 only defined for m < n
         yesD2[m][n] = 0;
         if (!globdist)
            getcim345D2u(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],
                         D2jumpuxxcoef[m][n],yesD2[m][n],m,n,index,0,0,mid,S,grid);
         else
            getcim345D2udist(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],D2uxxcoef[m][n],
                             D2jumpuxxcoef[m][n],yesD2[m][n],m,n,index,0,0,mid,S,grid);
      }



   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
         // one row of couping metrix, if interface in r, s
         if (gamma[r][(sk+1)/2] == 1)
         {
            D1u[r][(sk+1)/2] = new double[grid.dim];
            D1ucoef[r][(sk+1)/2] = new double ***[grid.dim];
            for (t = 0; t < grid.dim; t++)
               D1ucoef[r][(sk+1)/2][t] = matrix(N,N,N);
            D1uxcoef[r][(sk+1)/2] = matrix(grid.dim-1,grid.dim-1);
            D1uxxcoef[r][(sk+1)/2] = matrix(grid.dim-1,grid.dim-1);

            if (equal(index,eindex))
            {
               cout<<endl<<"dim = "<<r <<" sk = "<<sk << endl;
            }

            for (m = 0; m < grid.dim; m++)
               for (n = m+1; n < grid.dim; n++)
               {
                  if (!globdist)
                     getcim345D2u(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],
                                  D2uxxcoef[m][n],D2jumpuxxcoef[m][n],yesD2[m][n],
                                  m,n,index,r,sk,mid,S,grid);
                  else
                     getcim345D2udist(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],
                                      D2uxxcoef[m][n],D2jumpuxxcoef[m][n],yesD2[m][n],
                                      m,n,index,r,sk,mid,S,grid);

                  // #ifdef FIXBANANA 
                  //surgical fix for banana
                  // if ( m==0 && n==2  && GRIDNUM==110 && SURFOPT==13 && index[0]==18 && (index[1]==54 || index[1]==56) && (index[2]==11||index[2]==99))
                  if (globfixbanana)
                  {
                        // printf("fix (%d,%d) plane at index(%d,%d,%d)\n",m,n,index[0],index[1],index[2]);
                        // cout<<"fix "<<m<<" "<<n<<" at "<<index[0]<<","<<index[1]<<","<<index[2]<<endl;
                        double D2ueItf = getD2u(index,m,n,r,sk,alpha[r][(sk+1)/2],thesign,grid);//exact D2u at interface
                        vector<double***> D2ucoefvec;
                        vector<double*> D2uxcoefvec, D2uxxcoefvec;
                        vector<vector<int> > offsetvec;
                        vector<double> err;
                        bool yes = yescim5D2All(D2ucoefvec,D2uxcoefvec,D2uxxcoefvec,offsetvec,m,n,index,mid,S,grid);
                        D2u[m][n] = 0;
                        if(yes){
                           for(int i = 0; i < D2ucoefvec.size(); i++){
                              double D2uApx = evalcoef(D2u[m][n], D2ucoefvec[i],D2uxcoefvec[i],D2uxxcoefvec[i],index,0,0,0.0,mid,thesign,grid);
                              err.push_back(abs(D2uApx - D2ueItf));
                           }

                           int minIdx = std::min_element(err.begin(),err.end()) - err.begin();
                        
                           copyMat(D2ucoef[m][n],D2ucoefvec[minIdx],N,N,N);
                           copy(D2uxcoefvec[minIdx],D2uxcoefvec[minIdx]+3, D2uxcoef[m][n]);
                           copy(D2uxxcoefvec[minIdx],D2uxxcoefvec[minIdx]+3, D2uxxcoef[m][n]);

                        }
                        
                        for(int i = 0; i < D2ucoefvec.size(); i++){
                           free_matrix(D2ucoefvec[i],N,N,N);
                           delete [] D2uxcoefvec[i];
                           delete [] D2uxxcoefvec[i];
                        }
                     }
                  // #endif
                  if (equal(index,eindex)){
                     // should be O(h) approximate, exact is at interface, approx is at grid
                     cout <<"computed D2u in ("<< m<<","<<n<<") plane" << endl;
                     cout << " apprx = "
                          << evalcoef(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],
                                      D2uxxcoef[m][n],index,0,0,0.0,mid,thesign,grid) 
                          << " exact = "
                          << getD2u(index,m,n,r,sk,alpha[r][(sk+1)/2],thesign,grid)
                          << " error = "
                          << evalcoef(D2u[m][n],D2ucoef[m][n],D2uxcoef[m][n],
                                      D2uxxcoef[m][n],index,0,0,0.0,mid,thesign,grid)-
                             getD2u(index,m,n,r,sk,alpha[r][(sk+1)/2],thesign,grid) 
                          << endl;
                  }
               }
            
            rindex[r] = index[r]+sk;
// getting interface info
            getinterfaceinfo(alpha[r][(sk+1)/2],tangent,normal,S,index,rindex,grid);
            beta = 1.0-alpha[r][(sk+1)/2];

// get derivatives
//            cout << "getting Du" << endl;
            getcim345Du(D1u[r][(sk+1)/2],D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                        D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,index,r,sk,
                        alpha[r][(sk+1)/2],thesign,D2u,D2ucoef,D2uxcoef,D2uxxcoef,
                        D2jumpuxxcoef,mid,grid);
            double rhs = 0.0;

            if (equal(index,eindex))
            {
               cout <<endl<< "computed Du" << endl;
               for (m = 0; m < grid.dim; m++)
                  cout<<"dim "<< m 
                     <<" apprx = " << evalcoef(D1u[r][(sk+1)/2][m],D1ucoef[r][(sk+1)/2][m],
                                   D1uxcoef[r][(sk+1)/2][m],D1uxxcoef[r][(sk+1)/2][m],
                                   index,0,0,0.0,mid,thesign,grid) 
                      << ", exact ="
                      << getDu(index,m,r,sk,alpha[r][(sk+1)/2],thesign,grid)
                      << ", error ="
                      << evalcoef(D1u[r][(sk+1)/2][m],D1ucoef[r][(sk+1)/2][m],
                                   D1uxcoef[r][(sk+1)/2][m],D1uxxcoef[r][(sk+1)/2][m],
                                   index,0,0,0.0,mid,thesign,grid)-
                          getDu(index,m,r,sk,alpha[r][(sk+1)/2],thesign,grid) << endl;
               
               rhs += sk*grid.dx[r]*
                      evalcoef(D1u[r][(sk+1)/2][r],D1ucoef[r][(sk+1)/2][r],
                               D1uxcoef[r][(sk+1)/2][r],D1uxxcoef[r][(sk+1)/2][r],
                               index,0,0,0.0,mid,thesign,grid);
               cout << "rhs 1 = " << rhs << ", exact = "
                    << sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid) //rhs 1 = h ux-
                    << endl;
            }
// getting jump of uxx in r
            gettau(tau,index,r,sk,alpha[r][(sk+1)/2],grid);
//            cout << "getting jump Du" << endl;
            getcim345jumpux(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,
                            jumpD1jumpuxxcoef,index,r,sk,alpha[r][(sk+1)/2],thesign,
                            normal,tangent,mid,D1ucoef[r][(sk+1)/2],
                            D1uxcoef[r][(sk+1)/2],D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,
                            S,pb,grid);
            if (equal(index,eindex))
            {
               double x[grid.dim];
               sub2coord(x,index,grid);
               x[r] += sk*alpha[r][(sk+1)/2]*grid.dx[r];

               double computed_jump = evalcoef(jumpD1u,jumpD1ucoef,jumpD1uxcoef,
                                               jumpD1uxxcoef,jumpD1jumpuxxcoef,index,0,0,0.0,S,grid);
               double exact_jump = getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                                   getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid);

               cout<<endl;
               cout << "computed jump in D1u = "
                    << computed_jump
                    << " "
                    <<", exact = "
                    << exact_jump 
                    << "another exact "
                  //   << -(pb.epsilonp-pb.epsilonm)/ethere*
                  //       (getDu(index,0,r,sk,alpha[r][(sk+1)/2],thesign,grid)*normal[0]+
                  //        getDu(index,1,r,sk,alpha[r][(sk+1)/2],thesign,grid)*normal[1]+
                  //        getDu(index,2,r,sk,alpha[r][(sk+1)/2],thesign,grid)*normal[2])*
                  //       normal[r] << " "
                  //   << 2.0*x[r]*(1.0-pb.epsilonp/pb.epsilonm) 
                    << " error =  "
                    << computed_jump-exact_jump
                    << endl;

               rhs += -thesign*sk*beta*grid.dx[r]*computed_jump;
               cout << "rhs 2 = " << rhs << " exact =  "
                    << sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       thesign*sk*beta*grid.dx[r]*(exact_jump) << endl; //rhs2 = h ux- + beta h [ux], refer to eq 6 in ccim paper
            }
//            cout << "getting jump D2u" << endl;
            
            // if(a){
              // if a is not nullptr
              getcim345jumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,
                             alpha[r][(sk+1)/2],thesign,normal,mid,a,D1u[r][(sk+1)/2],
                             D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                             D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,D2u,D2ucoef,D2uxcoef,
                             D2uxxcoef,D2jumpuxxcoef,jumpD1u,jumpD1ucoef,jumpD1uxcoef,
                             jumpD1uxxcoef,jumpD1jumpuxxcoef,S,pb,grid);
            // }else{
            //   getcim345jumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,
            //                  alpha[r][(sk+1)/2],thesign,normal,mid,D1u[r][(sk+1)/2],
            //                  D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
            //                  D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,D2u,D2ucoef,D2uxcoef,
            //                  D2uxxcoef,D2jumpuxxcoef,jumpD1u,jumpD1ucoef,jumpD1uxcoef,
            //                  jumpD1uxxcoef,jumpD1jumpuxxcoef,S,pb,grid);  
            // }
            
            // error check
            if (equal(index,eindex))
            {
               cout << "computed jump in D2u = "
                    << evalcoef(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,0,0,
                                0.0,mid,thesign,grid) 
                    << ", exact = "
                    << getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                       getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],-1,grid)
                    << ", error = "
                    << evalcoef(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,0,0,
                                0.0,mid,thesign,grid)-
                       (getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],-1,grid)) 
                    << endl;
               rhs += -thesign*0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                      evalcoef(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,0,0,
                               0.0,mid,thesign,grid); 
               cout << "rhs 3 = " << rhs << " " //rhs4 = h ux- + beta h [ux] + (1/2) (betta h)^2  [uxx]
                    << ", exact = "
                    << -thesign*sk*beta*grid.dx[r]*
                        (getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                         getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       thesign*0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                       (getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],-1,grid)) << endl;
               
               rhs += 0.5*(beta*beta-alpha[r][(sk+1)/2]*alpha[r][(sk+1)/2])*grid.dx[r]*grid.dx[r]*
                      getD2u(index,r,r,0,0,0.0,thesign,grid); // uxx in dim r
               
               cout << "rhs 4 = " << rhs << " " //rhs4 = h ux- + beta h [ux] + (1/2) (beta h)^2  [uxx] + (1/2) h^2 (beta^2-alpha^2) (uxx-)
                    << ", exact = "
                    << -thesign*sk*beta*grid.dx[r]*
                        (getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                         getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       thesign*0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                       (getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       0.5*(beta*beta-alpha[r][(sk+1)/2]*alpha[r][(sk+1)/2])*grid.dx[r]*grid.dx[r]*
                       getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],thesign,grid) << endl;
               
               cout << 0.5*(beta*beta-alpha[r][(sk+1)/2]*alpha[r][(sk+1)/2])*grid.dx[r]*grid.dx[r]*getD2u(index,r,r,0,0,0.0,thesign,grid) << " "
                    << 0.5*(beta*beta-alpha[r][(sk+1)/2]*alpha[r][(sk+1)/2])*grid.dx[r]*grid.dx[r]*getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],thesign,grid) << " "
                    << getD2u(index,r,r,0,0,0.0,thesign,grid) << " "
                    << getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],thesign,grid) << endl;
               
               int tempindex[grid.dim];
               for (m = 0; m < grid.dim; m++)
                  tempindex[m] = index[m];
               tempindex[r] += sk; //tempindex at the other side of interface
               cout << "rhs = " << rhs << endl;
                      // asume cim12 p25, (u-) - alpha (ux-) + 0.5 (alpha h)^2 uxx-
               cout << getu(index,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       sk*alpha[r][(sk+1)/2]*grid.dx[r]*
                       getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)+
                       0.5*alpha[r][(sk+1)/2]*alpha[r][(sk+1)/2]*grid.dx[r]*grid.dx[r]*
                       getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       
                       (getu(tempindex,r,-sk,beta,-thesign,grid)+
                        sk*beta*grid.dx[r]*
                        getDu(tempindex,r,r,-sk,beta,-thesign,grid)+
                        0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                        getD2u(tempindex,r,r,r,-sk,beta,-thesign,grid));
               cout << ", should ="
                    << thesign*tau
                       -thesign*sk*beta*grid.dx[r]*
                        (getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                         getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       thesign*0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                       (getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       0.5*grid.dx[r]*grid.dx[r]*(beta*beta-alpha[r][(sk+1)/2]*
                                                            alpha[r][(sk+1)/2])*
                       getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],thesign,grid) << endl;
              
               cout << "thesign = " << thesign << " and dim = " << r << " and direction = " << sk << endl;
               cout << " move to same side should be zero "//move everything to one side
                    << getu(index,0,0,0.0,thesign,grid)-
                       getu(index,r,sk,1.0,-thesign,grid)
                       + (-thesign)*tau
                       -thesign*sk*beta*grid.dx[r]*
                        (getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                         getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       thesign*0.5*beta*beta*grid.dx[r]*grid.dx[r]*
                       (getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],-1,grid))+
                       0.5*grid.dx[r]*grid.dx[r]*(beta*beta-alpha[r][(sk+1)/2]*
                                                            alpha[r][(sk+1)/2])*
                       getD2u(index,r,r,r,sk,alpha[r][(sk+1)/2],thesign,grid) << endl;
               
               double uhere = getu(index,0,0,0.0,thesign,grid);
               double uthere = getu(index,r,sk,1.0,-thesign,grid);
               cout << "(at grid) u_here - u_there "
                    << " = " << uhere 
                    << " - " << uthere
                    << " = "<<uhere - uthere << endl;
               
               cout << "u jump (u+ - u-) = "
                    << getu(index,r,sk,alpha[r][(sk+1)/2],1,grid) - getu(index,r,sk,alpha[r][(sk+1)/2],-1,grid) << endl;
               
               double x[grid.dim];
               sub2coord(x,index,grid);
               x[r] += sk*alpha[r][(sk+1)/2]*grid.dx[r];
               cout << "error in radius = " <<sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-grid.radius0 << endl;
            }
// form d0 and dcoef's rth entry, d0 + dcoeff u = G [ux uxx]
            d0[grid.dim*(sk+1)/2+r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+
                                               sk*beta*jumpD1u/grid.dx[r]+
                                               0.5*beta*beta*jumpD2u)-
                                      sk*D1u[r][(sk+1)/2][r]/grid.dx[r];
            
            exactd[grid.dim*(sk+1)/2+r] = d0[grid.dim*(sk+1)/2+r];//constant part is exact

            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            { 
                //coefficient of u[sindex], move rhs to lhs
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                           thesign*(sk*beta*evalarray(jumpD1ucoef,sindex)/grid.dx[r]+
                                    0.5*beta*beta*evalarray(jumpD2ucoef,sindex))-
                           sk*evalarray(D1ucoef[r][(sk+1)/2][r],sindex)/grid.dx[r]);
               for (m = 0; m < grid.dim; m++)
                  tempindex[m] = index[m]+sindex[m]-mid;
               if (evalarray(S,tempindex) < 0.0)
                  tempsign = -1.0;
               else
                  tempsign = 1.0;

               if (equal(index,eindex))
               {
                  exactd[grid.dim*(sk+1)/2+r] += evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)*
                                                             getu(tempindex,0,0,0.0,tempsign,grid);
               }
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = mid;
             // coeff due to ( u_{i,j+1} - u_{i,j} )/ h^2 
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                        evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)-
                        1.0/(grid.dx[r]*grid.dx[r]));
            for (m = 0; m < grid.dim; m++)
               tempindex[m] = index[m]+sindex[m]-mid;
            if (evalarray(S,tempindex) < 0.0)
               tempsign = -1.0;
            else
               tempsign = 1.0;

            if (equal(index,eindex)){
               exactd[grid.dim*(sk+1)/2+r] += -1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(tempindex,0,0,0.0,tempsign,grid);   
            }
            
            sindex[r] = mid+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,
                        evalarray(dcoef[grid.dim*(sk+1)/2+r],sindex)+
                        1.0/(grid.dx[r]*grid.dx[r]));
            for (m = 0; m < grid.dim; m++)
               tempindex[m] = index[m]+sindex[m]-mid;
            
            if (evalarray(S,tempindex) < 0.0)
               tempsign = -1.0;
            else
               tempsign = 1.0;

            if (equal(index,eindex)){
               exactd[grid.dim*(sk+1)/2+r] += 1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(tempindex,0,0,0.0,tempsign,grid);
            }

            if (equal(index,eindex))
            {
               double somezeros[grid.dim];
               for (m = 0; m < grid.dim; m++)
                  somezeros[m] = 0.0;
               cout << "Row = " <<grid.dim*(sk+1)/2+r
                    << ", exactd = "
                    << exactd[grid.dim*(sk+1)/2+r]
                    << ", approx = "
                    << evalcoef(d0[grid.dim*(sk+1)/2+r],dcoef[grid.dim*(sk+1)/2+r],
                                somezeros,somezeros,index,0,0,0.0,mid,S,grid)
                    << endl;
            }
            sindex[r] = mid;
// form matrix's rth row
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m] = -thesign*
                                           (sk*beta*jumpD1uxxcoef[m]/grid.dx[r]+
                                            0.5*beta*beta*jumpD2uxxcoef[m])+
                                           sk*D1uxxcoef[r][(sk+1)/2][r][m]/grid.dx[r]; // coeff of u_{xx} from [u_xx] [u_x]

            G[grid.dim*(sk+1)/2+r][r] += 0.5*(beta*beta-alpha[r][(sk+1)/2]*alpha[r][(sk+1)/2]); // 1/2 (beta^2-alpha^2)u_{xx}
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m+grid.dim] = -thesign*
                                                    (sk*beta*jumpD1uxcoef[m]/grid.dx[r]+
                                                     0.5*beta*beta*jumpD2uxcoef[m])+
                                                    sk*D1uxcoef[r][(sk+1)/2][r][m]/
                                                    grid.dx[r];// coeff of u_x from [u_xx] [u_x]
            
            // check r-th row of G
            if (equal(index,eindex))
            {  
               cout<<endl<<"check r-th row of G"<<endl;
               exactres = exactd[grid.dim*(sk+1)/2+r];
               for (m = 0; m < 2*grid.dim; m++)
                  if (m < grid.dim)
                     exactres -= G[grid.dim*(sk+1)/2+r][m]*
                                 getD2u(index,m,m,0,0,0.0,thesign,grid);
                  else
                     exactres -= G[grid.dim*(sk+1)/2+r][m]*
                                 getDu(index,m-grid.dim,0,0,0.0,thesign,grid);
               cout << "d_exact - G_apprx DDu_exact = " << exactres << endl;
               cout << "tau = "
                    << getu(index,r,sk,alpha[r][(sk+1)/2],1,grid)-
                       getu(index,r,sk,alpha[r][(sk+1)/2],-1,grid)
                    << ", tau/h^2  ="
                    << (getu(index,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getu(index,r,sk,alpha[r][(sk+1)/2],-1,grid))/
                       (grid.dx[r]*grid.dx[r]) << endl;
               cout << "-ehere lap(u) - f = "
                    << -ehere*(getD2u(index,0,0,0,0,0.0,thesign,grid)+
                               getD2u(index,1,1,0,0,0.0,thesign,grid)+
                               getD2u(index,2,2,0,0,0.0,thesign,grid))-
                       getf(index,0,0,0.0,thesign,pb,grid) 
                    << " =  "
                    << -ehere*(getD2u(index,0,0,0,0,0.0,thesign,grid)+
                               getD2u(index,1,1,0,0,0.0,thesign,grid)+
                               getD2u(index,2,2,0,0,0.0,thesign,grid)) 
                    << " - "
                    << getf(index,0,0,0.0,thesign,pb,grid) 
                    << endl
                    << "uxx = "<<getD2u(index,0,0,0,0,0.0,thesign,grid) << " "
                    << "uyy = "<<getD2u(index,1,1,0,0,0.0,thesign,grid) << " "
                    << "uzz = "<<getD2u(index,2,2,0,0,0.0,thesign,grid) << " "
                    << "ehere = "<<ehere 
                    << endl;

               double x[grid.dim];
               sub2coord(x,index,grid);
               x[r] += sk*alpha[r][(sk+1)/2]*grid.dx[r];
               cout <<"at interface x = ("<< x[0] << ", " << x[1] << ", " << x[2] <<")"<< endl;
               cout <<"error of radius" << sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])-grid.radius
                    << ", alpha = "<< alpha[r][(sk+1)/2] << endl;
            }
            rindex[r] = index[r];
         }
         else
         {
// if interior point
            d0[grid.dim*(sk+1)/2+r] = 0.0;
            exactd[grid.dim*(sk+1)/2+r] = d0[grid.dim*(sk+1)/2+r];
            for (s = 0; s < 2*grid.dim; s++)
               G[grid.dim*(sk+1)/2+r][s] = 0.0;
            G[grid.dim*(sk+1)/2+r][r] = 0.5;
            G[grid.dim*(sk+1)/2+r][grid.dim+r] = sk/grid.dx[r];
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
               setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,0.0);
   
               (sindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
               {
                  sindex[i] = 0;
                  (sindex[i-1])++;
               }
            }
            for (s = 0; s < grid.dim; s++)
               sindex[s] = mid;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,-1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = mid+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,1.0/(grid.dx[r]*grid.dx[r]));
            sindex[r] = mid;

            // compute exact rhs
            if (equal(index,eindex))
            {
               exactd[grid.dim*(sk+1)/2+r] += -1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(index,0,0,0.0,thesign,grid);
               exactd[grid.dim*(sk+1)/2+r] += 1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(index,r,sk,1.0,thesign,grid);
            }
            
         }

   // print coupling matrix, exact rhs
   if (equal(index,eindex))
   {
      cout<<endl<<"G matrix"<<endl;
      for (m = 0; m < 2*grid.dim; m++)
      {
         for (n = 0; n < 2*grid.dim; n++)
            cout << setw(10)<< G[m][n] << " ";
         cout << endl;
      }
      cout<<"exactd"<<endl;
      for (m = 0; m < 2*grid.dim; m++)
         cout << exactd[m] << endl;
      gecp0(LU,PLR,PLC,G,2*grid.dim-1,2*grid.dim-1);
      forwardbacksub0(exactd,exactd,LU,PLR,PLC,2*grid.dim-1);

      cout<<"solved var with exact rhs | exact var "<<endl;
      for (m = 0; m < 2*grid.dim; m++)
      {
         cout  << exactd[m] << " | ";
         if (m < grid.dim)
            cout << getD2u(index,m,m,0,0,0.0,thesign,grid) << endl;
         else
            cout << getDu(index,m-grid.dim,0,0,0.0,thesign,grid) << endl;
      }
      cout <<" solved -eps lap(u) - f =" <<-ehere*(exactd[0]+exactd[1]+exactd[2])-getf(index,0,0,0.0,thesign,pb,grid)
           << " =  "<< -ehere*(exactd[0]+exactd[1]+exactd[2]) 
           << " - " << getf(index,0,0,0.0,thesign,pb,grid) << endl;
      
   }

   free_matrix(D1jumpuxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);

   
   free_matrix(D2uxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(D2uxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);
   free_matrix(D2jumpuxxcoef,grid.dim-1,grid.dim-1);

   free_matrix(jumpD1ucoef,N,N,N);
   free_matrix(jumpD1jumpuxxcoef,grid.dim-1,grid.dim-1);

   free_matrix(jumpD2ucoef,N,N,N);

}

// estimate condition number
double condest(const vector<vector<double>> &A){
   int n = A.size();
   double n1 = norm1(A);
   vector<double> x(n, 1.0/n);
   int i1 = -1, i2;
   double c1 = 0, c2 = 0;
   

   while(true){
      x = gesolve(A,x);

      c2 = sum(abs(x));

      x = sign(x);

      x = gesolve(transpose(A),x);

      i2 = max_abs_idx(x);

      if(0 <= i1){
         if(i1==i2 || c2<=c1){
            break;
         }
      }

      i1 = i2;
      c1 = c2;

      fill(x.begin(),x.end(),0);
      x[i1] = 1.0;
   }

   return c2*n1;
}

double condest(double **G, int n){
   return condest(Pointer2dToVector(G,n,n));
}





//
void cim345cond(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize, 
            int *index, double***a, int gamma[][2], double ***S, PBData &pb, GridData &grid){
   int r, s, t, i, j, m, n;

   int mid = 2;
   int N = 2*mid;
   
   double **alpha = matrix(grid.dim-1,1);


   double ***D1u, ******D1ucoef, ****D1uxcoef, ****D1uxxcoef;
   D1u = new double **[grid.dim];
   D1ucoef = new double *****[grid.dim];
   D1uxcoef = new double ***[grid.dim];
   D1uxxcoef = new double ***[grid.dim];
   for (r = 0; r < grid.dim; r++)
   {
      D1u[r] = new double *[2];
      D1ucoef[r] = new double ****[2];
      D1uxcoef[r] = new double **[2];
      D1uxxcoef[r] = new double **[2];
      for (s = 0; s <= 1; s++)
      {
         D1u[r][s] = NULL;
         D1ucoef[r][s] = NULL;
         D1uxcoef[r][s] = NULL;
         D1uxxcoef[r][s] = NULL;
      }
   }

   double **D2u, *****D2ucoef;
   D2u = matrix(grid.dim-1,grid.dim-1);
   D2ucoef = new double ****[grid.dim];
   for (r = 0; r < grid.dim; r++)
   {
      D2ucoef[r] = new double ***[grid.dim];
      for (s = 0; s < grid.dim; s++)
         D2ucoef[r][s] = matrix(N,N,N);
   }


   // coupling matrix, G, d0, dcoef
   double** LU = matrix(2*grid.dim-1,2*grid.dim-1);
   double** G = matrix(2*grid.dim-1,2*grid.dim-1);
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double d0[2*grid.dim];
   double ***dcoef[2*grid.dim];
   for (int r = 0; r < 2*grid.dim; r++)
      dcoef[r] = matrix(N,N,N);

   // choose coupling matrix based on estimated condition number
   char tempglobdist = globdist;

  
   gmatrix(G, LU, PLR, PLC,
    d0, dcoef, 
    D1u, D1ucoef,  D1uxcoef,  D1uxxcoef,
    D2u,  D2ucoef,
    index,
    alpha,  gamma, S, a, pb, grid);

   double estcond0 = condest(G, 6);

   // if globcim == 6, i.e. with condition number, compute the other globdist,
   if (globcim == 6){

      globdist = 1 - globdist;   
      gmatrix(G, LU, PLR, PLC,
       d0, dcoef, 
       D1u, D1ucoef,  D1uxcoef,  D1uxxcoef,
       D2u,  D2ucoef,
       index,
       alpha,  gamma, S, a, pb, grid);
      double estcond1= condest(G, 6);

      globdist = (estcond0<estcond1)?0:1;
      
      gmatrix(G, LU, PLR, PLC,
       d0, dcoef, 
       D1u, D1ucoef,  D1uxcoef,  D1uxxcoef,
       D2u,  D2ucoef,
       index,
       alpha,  gamma, S, a, pb, grid);   


      // reset
      globdist = tempglobdist;
   }
   




   double temp[2*grid.dim];

   // === figure out coeff
   double ehere, ethere, thesign;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim];


   if (evalarray(S,index) < 0.0)
   {
      ehere = pb.epsilonm;
      ethere = pb.epsilonp;
      thesign = -1.0;
   }
   else
   {
      ehere = pb.epsilonp;
      ethere = pb.epsilonm;
      thesign = 1.0;
   }

   // store D2u_rr as function of u-value in Dusmall
   // ux = const term, uxcoef = coef of u
   double ux[grid.dim], ****uxcoef;
   uxcoef = new double ***[grid.dim];
   for (r = 0; r < grid.dim; r++)
      uxcoef[r] = matrix(N,N,N);

   // ==== solve coupling matrix
   gecp0(LU,PLR,PLC,G,2*grid.dim-1,2*grid.dim-1);
   forwardbacksub0(temp,d0,LU,PLR,PLC,2*grid.dim-1);
   double value = 0.0;
   // constant term in D2u[n][n] (first dim term), put to rhs
   for (int n = 0; n < grid.dim; n++)
      value += ehere*temp[n];
   setvalarray(b,index,getf(index,0,0,0.0,thesign,pb,grid)+value);

   // const term for D2u[j][j]
   for (j = 0; j < grid.dim; j++)
      D2u[j][j] = temp[j];

   // const term for Du[j]
   for (j = 0; j < grid.dim; j++)
      ux[j] = temp[j+grid.dim];

   for (s = 0; s < grid.dim; s++)
      tindex[s] = index[s];
   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      // temp is vector for u[sindex]
      for (n = 0; n < 2*grid.dim; n++)
         temp[n] = evalarray(dcoef[n],sindex);

      // A(index, tindex) = value
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      for (s = 0; s < grid.dim; s++)
         tindex[s] = index[s]-mid+sindex[s];
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value -= ehere*temp[n];
      if (value != 0.0)
         sparse2(index,tindex,A,value,grid);

      // D2u[n][n] in terms of u-val
      for (n = 0; n < grid.dim; n++)
         setvalarray(D2ucoef[n][n],sindex,temp[n]);

      // Du[n] in terms of u-val
      for (n = 0; n < grid.dim; n++)
         setvalarray(uxcoef[n],sindex,temp[n+grid.dim]);

      (sindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && sindex[i] > N; i--)
      {
         sindex[i] = 0;
         (sindex[i-1])++;
      }
   }
   for (s = 0; s < grid.dim; s++)
      sindex[s] = mid;

// add value of a to diagonals
   if(a){
      value = evalarray(a,index);
      sparse2(index,index,A,value,grid);   
   }
   

// storing uint info
   addtostorage( Dusmall, buildsize, mid, index, gamma, alpha, 
    ux,  uxcoef, 
    D1u, D1ucoef,  D1uxcoef,  D1uxxcoef,
    D2u,  D2ucoef, grid);


   free_matrix(alpha,grid.dim-1,1);
   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(G,2*grid.dim-1,2*grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(dcoef[r],N,N,N);



   for (r = 0; r < grid.dim; r++)
      for (int sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
            delete [] D1u[r][(sk+1)/2];
            for (t = 0; t < grid.dim; t++)
               free_matrix(D1ucoef[r][(sk+1)/2][t],N,N,N);
            delete [] D1ucoef[r][(sk+1)/2];
            free_matrix(D1uxcoef[r][(sk+1)/2],grid.dim-1,grid.dim-1);
            free_matrix(D1uxxcoef[r][(sk+1)/2],grid.dim-1,grid.dim-1);
         }
   for (r = 0; r < grid.dim; r++)
   {
      delete [] D1u[r];
      delete [] D1ucoef[r];
      delete [] D1uxcoef[r];
      delete [] D1uxxcoef[r];
   }
   delete [] D1u;
   delete [] D1ucoef;
   delete [] D1uxcoef;
   delete [] D1uxxcoef;


   free_matrix(D2u,grid.dim-1,grid.dim-1);
   for (r = 0; r < grid.dim; r++)
   {
      for (s = 0; s < grid.dim; s++)
         free_matrix(D2ucoef[r][s],N,N,N);
      delete [] D2ucoef[r];
   }
   delete [] D2ucoef;


  
}








//return a vector of all D2u scheme
char yescim5D2All(vector<double***> &D2ucoefvec, vector<double*> &D2uxcoefvec, vector<double*> &D2uxxcoefvec, vector<vector<int>> &offsetvec, int m, int n, 
               int *index, int mid, double ***S, GridData &grid)
{
   int i, r, s, thesign, N = 2*mid, ntheta = 8, signs[ntheta], sk[2], offset[2][2];
   int tindex[grid.dim], sindex[grid.dim];
   double theta0, theta, dtheta, themax, value;

   for (i = 0; i < grid.dim; i++)
      sindex[i] = mid;
   for (r = 0; r < grid.dim; r++)
      tindex[r] = index[r];

   thesign = 2*(evalarray(S,tindex) >= 0.0)-1;//1=outside, -1 indside
   dtheta = 2.0*M_PI/ntheta;
   theta0 = -3.0*M_PI/4.0;// start from bottom left corner
// computes signs of points surrounding node, start from (-1,-1) corner, couter-clockwise
   for (r = 0; r < ntheta; r++)
   {
      theta = theta0+r*dtheta;
      themax = max(fabs(cos(theta)),fabs(sin(theta)));
      offset[0][0] = round(cos(theta)/themax);
      offset[0][1] = round(sin(theta)/themax);
      tindex[m] = index[m]+offset[0][0];
      tindex[n] = index[n]+offset[0][1];
      signs[r] = 2*(evalarray(S,tindex) >= 0.0)-1;
   }
   tindex[m] = index[m];
   tindex[n] = index[n];

// looks for central differencing possibility around node
   for (r = 0; r < ntheta; r += 2)
      if ((thesign < 0)+(signs[r] < 0) != 1 && 
          (thesign < 0)+(signs[(r+2)%ntheta] < 0) != 1)
      {
         double ***D2ucoef = matrix(N,N,N);
         double *D2uxcoef = new double[3]();
         double *D2uxxcoef = new double[3]();

         theta = theta0+r*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[0][0] = round(cos(theta)/themax);
         offset[0][1] = round(sin(theta)/themax);
         theta = theta0+((r+2)%ntheta)*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[1][0] = round(cos(theta)/themax);
         offset[1][1] = round(sin(theta)/themax);
         sk[0] = offset[0][0]-offset[1][0];
         sk[1] = offset[0][1]-offset[1][1];
         if (sk[0] == 0)
         {
            sk[0] = offset[0][0];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[n] = -1.0/(sk[0]*grid.dx[m]);
         }
         else
         {
            sk[1] = offset[0][1];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[m] = -1.0/(sk[1]*grid.dx[n]);
         }
         sindex[m] = mid;
         sindex[n] = mid;

         offsetvec.push_back(vector<int>{offset[0][0],offset[0][1],offset[1][0],offset[1][1]});
         D2ucoefvec.push_back(D2ucoef);
         D2uxcoefvec.push_back(D2uxcoef);
         D2uxxcoefvec.push_back(D2uxxcoef);
      }

// looks for forward differencing possibility around node
   for (r = 0; r < ntheta; r++)
      if ((thesign < 0)+(signs[r] < 0) != 1 && 
          (thesign < 0)+(signs[(r+1)%ntheta] < 0) != 1)
      {
         double ***D2ucoef = matrix(N,N,N);
         double *D2uxcoef = new double[3]();
         double *D2uxxcoef = new double[3]();

         theta = theta0+r*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[0][0] = round(cos(theta)/themax);
         offset[0][1] = round(sin(theta)/themax);
         theta = theta0+((r+1)%ntheta)*dtheta;
         themax = max(fabs(cos(theta)),fabs(sin(theta)));
         offset[1][0] = round(cos(theta)/themax);
         offset[1][1] = round(sin(theta)/themax);
         sk[0] = offset[0][0]-offset[1][0];
         sk[1] = offset[0][1]-offset[1][1];
         if (sk[0] == 0)
         {
            sk[0] = offset[0][0];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[n] = -1.0/(sk[0]*grid.dx[m]);
//            D2uxxcoef[n] = -0.5*sk[1]*grid.dx[n]/(sk[0]*grid.dx[m]);
            D2uxxcoef[n] = -0.5*(offset[0][1]+offset[1][1])*grid.dx[n]/
                            (sk[0]*grid.dx[m]);
         }
         else
         {
            sk[1] = offset[0][1];
            value = 1.0/(sk[0]*sk[1]*grid.dx[m]*grid.dx[n]);
            sindex[m] = mid+offset[0][0];
            sindex[n] = mid+offset[0][1];
            setvalarray(D2ucoef,sindex,value);
            sindex[m] = mid+offset[1][0];
            sindex[n] = mid+offset[1][1];
            setvalarray(D2ucoef,sindex,-value);
            D2uxcoef[m] = -1.0/(sk[1]*grid.dx[n]);
//            D2uxxcoef[m] = -0.5*sk[0]*grid.dx[m]/(sk[1]*grid.dx[n]);
            D2uxxcoef[m] = -0.5*(offset[0][0]+offset[1][0])*grid.dx[m]/
                           (sk[1]*grid.dx[n]);
         }
         sindex[m] = mid;
         sindex[n] = mid;

         offsetvec.push_back(vector<int>{offset[0][0],offset[0][1],offset[1][0],offset[1][1]});
         D2ucoefvec.push_back(D2ucoef);
         D2uxcoefvec.push_back(D2uxcoef);
         D2uxxcoefvec.push_back(D2uxxcoef);
      }
   if (offsetvec.size() == 0){
      return 0;
   }
   return 1;
   
}


