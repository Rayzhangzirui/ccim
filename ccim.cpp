#include "grid.h"
#include "global.h"
#include "ccim.h"
#include "interface.h"
#include "numerics.h"
#include <cmath>
#include "math.h"
#include <algorithm>
#include <iostream>
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
              copy(combimation[k],combimation[k]+4,tempsk2);
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



// get jump in Du in rstar dim, cim12.pdf p27
// Du-[r] = NxNxN D1ucoeff[r] u + 1x3 D1uxcoef[r] ux + 1x3 D1uxxcoef[r] uxx + 3x3 D1jumpuxxcoef [uxx]
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
   
   u0 = sigma/ethere*normal[rstar]+getdotprod(Dtau,tangent,grid.dim)*tangent[rstar];//constant term

   // recast Du-, [epsilon]/epsilonp dot(grad(Du-),n) n[star]
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


// without a, [DDu_{rstar}] = u0 + ucoef NxNxN u-value + uxcoef Du_{1,2,3} + uxxcoef 1x3 DDu_{1,2,3}
// recast jumpuxxcoeff to ucoef, uxcoef and uxxcoef
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

   // free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   // free_matrix(B,2*grid.dim-1,2*grid.dim-1);
   free_matrix(Dn,grid.dim-1,grid.dim-1);
   free_matrix(D2tau,grid.dim-1,grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(bcoef[r],N,N,N);
}

// with a
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
         Dndott[n] += Dn[n][m]*tangent1[m];
         Dndots[n] += Dn[n][m]*tangent2[m];
         dotD2taudot[0] += tangent1[n]*D2tau[n][m]*tangent1[m];
         dotD2taudot[1] += tangent2[n]*D2tau[n][m]*tangent2[m];
         dotD2taudot[2] += tangent1[n]*D2tau[n][m]*tangent2[m];
      }
      dotDndot[0] += tangent1[n]*Dndott[n];
      dotDndot[1] += tangent2[n]*Dndots[n];
      dotDndot[2] += tangent1[n]*Dndots[n];
      Dtaudot[0] += Dtau[n]*normal[n];
      Dtaudot[1] += Dtau[n]*tangent1[n];
      Dtaudot[2] += Dtau[n]*tangent2[n];
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
            value += D1jumpuxxcoef[r][m][n]*normal[r];
         for (r = 0; r < grid.dim; r++)
            B[r][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*value*dotDndot[r];
         value = 0.0;
         for (r = 0; r < grid.dim; r++)
            value += D1jumpuxxcoef[r][m][n]*Dndott[r];
         B[grid.dim][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*value;
         value = 0.0;
         for (r = 0; r < grid.dim; r++)
            value += D1jumpuxxcoef[r][m][n]*Dndots[r];
         B[grid.dim+1][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*value;
// get contributions of D2u to jump data
         B[grid.dim][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*
                                       (normal[m]*tangent1[n]+normal[n]*tangent1[m])*
                                       D2jumpuxxcoef[m][n];
         B[grid.dim+1][two2one[m][n]] -= thesign*(ethere-ehere)/ethere*
                                         (normal[m]*tangent2[n]+normal[n]*tangent2[m])*
                                         D2jumpuxxcoef[m][n];
      }

   for (n = 0; n < grid.dim; n++)
      b0[n] = (sigma/ethere-Dtaudot[0])*dotDndot[n]+dotD2taudot[n];
   b0[grid.dim] = getdotprod(Dsigma,tangent1,grid.dim)/ethere-
                  Dtaudot[1]*dotDndot[0]-Dtaudot[2]*dotDndot[2];
   b0[grid.dim+1] = getdotprod(Dsigma,tangent2,grid.dim)/ethere-
                    Dtaudot[1]*dotDndot[2]-Dtaudot[2]*dotDndot[1];
   b0[grid.dim+2] = -jumpfe+aethere*tau;
//   b0[grid.dim+2] = getD2u(index,0,0,rstar,sk,alpha,1.0,grid)-
//                    getD2u(index,0,0,rstar,sk,alpha,-1.0,grid)+
//                    getD2u(index,1,1,rstar,sk,alpha,1.0,grid)-
//                    getD2u(index,1,1,rstar,sk,alpha,-1.0,grid)+
//                    getD2u(index,2,2,rstar,sk,alpha,1.0,grid)-
//                    getD2u(index,2,2,rstar,sk,alpha,-1.0,grid);

// get b coefs
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
         value += evalarray(D1ucoef[n],sindex)*normal[n];
      for (m = 0; m < grid.dim; m++)
         setvalarray(bcoef[m],sindex,evalarray(bcoef[m],sindex)+
                                     thesign*(ethere-ehere)/ethere*value*dotDndot[m]);
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value += evalarray(D1ucoef[n],sindex)*Dndott[n];
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
// addition
//   double x[grid.dim];
//   sub2coord(x,index,grid);
//   x[rstar] += sk*alpha*grid.dx[rstar];
//   jumpae = cos(x[grid.dim-1])-sin(x[0]);
   setvalarray(bcoef[grid.dim+2],sindex,jumpae);

// get bx and bxx coefs
   for (s = 0; s < grid.dim; s++)
   {
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
//   cout << index[0] << " " << index[1] << " " << index[2] << " "
//        << evalcoef(b0[grid.dim+2],bcoef[grid.dim+2],bxcoef[grid.dim+2],
//                    bxxcoef[grid.dim+2],index,0,0,0.0,mid,S,grid) << " ";
//   cout << -jumpfe+jumpae*getu(index,rstar,sk,alpha,thesign,grid) << endl;
//   getchar();
//   cout << index[0] << " " << index[1] << " " << index[2] << " "
//        << evalcoef(b0[grid.dim+2],bcoef[grid.dim+2],bxcoef[grid.dim+2],
//                    bxxcoef[grid.dim+2],index,0,0,0.0,mid,S,grid) << " "
//        << getD2u(index,0,0,rstar,sk,alpha,1.0,grid)-
//           getD2u(index,0,0,rstar,sk,alpha,-1.0,grid)+
//           getD2u(index,1,1,rstar,sk,alpha,1.0,grid)-
//           getD2u(index,1,1,rstar,sk,alpha,-1.0,grid)+
//           getD2u(index,2,2,rstar,sk,alpha,1.0,grid)-
//           getD2u(index,2,2,rstar,sk,alpha,-1.0,grid) << endl;
//   getchar();

   gecp0(LU,PLR,PLC,B,2*grid.dim-1,2*grid.dim-1);

// form jumpuxx in rstar direction and also recast Du and D2u
   forwardbacksub0(temp,b0,LU,PLR,PLC,2*grid.dim-1);
   u0 = temp[two2one[rstar][rstar]];
   for (s = 0; s < grid.dim; s++)
      for (m = 0; m < grid.dim; m++)
         for (n = m+1; n < grid.dim; n++)
            D1u[s] += D1jumpuxxcoef[s][m][n]*temp[two2one[m][n]];
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         D2u[m][n] += D2jumpuxxcoef[m][n]*temp[two2one[m][n]];
         jumpD1u += jumpD1jumpuxxcoef[m][n]*temp[two2one[m][n]];
      }

   for (s = 0; s < grid.dim; s++)
      sindex[s] = 0;
   while (sindex[0] <= N)
   {
      for (m = 0; m < 2*grid.dim; m++)
         temp[m] = evalarray(bcoef[m],sindex);
      forwardbacksub0(temp,temp,LU,PLR,PLC,2*grid.dim-1);
      setvalarray(ucoef,sindex,temp[two2one[rstar][rstar]]);
      for (s = 0; s < grid.dim; s++)
         for (m = 0; m < grid.dim; m++)
            for (n = m+1; n < grid.dim; n++)
               setvalarray(D1ucoef[s],sindex,evalarray(D1ucoef[s],sindex)+
                                             D1jumpuxxcoef[s][m][n]*
                                             temp[two2one[m][n]]);
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

// with ***a, 
void cim345(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize, 
            int *index, double ***a, int gamma[][2], double ***S, PBData &pb, 
            GridData &grid)
{
   int r, s, t, i, j, m, n, sk, mid = 2, N = 2*mid;
   int rindex[grid.dim], tindex[grid.dim], sindex[grid.dim], Narray[grid.dim];
   double d0[2*grid.dim];
   double ***dcoef[2*grid.dim];
   double **LU, **G, value;
   int PLR[2*grid.dim], PLC[2*grid.dim];
   double **alpha = matrix(grid.dim-1,1), beta, normal[grid.dim], tangent[grid.dim];
   double ethere, ehere;
   double temp[2*grid.dim];
   double sigma, tau, Dsigma[grid.dim], Dtau[grid.dim], thesign;

            double exactd[2*grid.dim], exactres, tempsign;
            int tempindex[grid.dim];

   for (r = 0; r < grid.dim; r++)
      Narray[r] = N;

   double ***D1u, ******D1ucoef, ****D1uxcoef, ****D1uxxcoef, ***D1jumpuxxcoef;
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
   D1jumpuxxcoef = matrix(grid.dim-1,grid.dim-1,grid.dim-1);
   
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

   double ux[grid.dim], ****uxcoef;
   uxcoef = new double ***[grid.dim];
   for (r = 0; r < grid.dim; r++)
      uxcoef[r] = matrix(N,N,N);

   char yesD2[grid.dim][grid.dim];
   

   value = 0.0;

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
         if (gamma[r][(sk+1)/2] == 1)
         {
            D1u[r][(sk+1)/2] = new double[grid.dim];
            D1ucoef[r][(sk+1)/2] = new double ***[grid.dim];
            for (t = 0; t < grid.dim; t++)
               D1ucoef[r][(sk+1)/2][t] = matrix(N,N,N);
            D1uxcoef[r][(sk+1)/2] = matrix(grid.dim-1,grid.dim-1);
            D1uxxcoef[r][(sk+1)/2] = matrix(grid.dim-1,grid.dim-1);

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
// getting jump of uxx in r
            gettau(tau,index,r,sk,alpha[r][(sk+1)/2],grid);
//            cout << "getting jump Du" << endl;
            getcim345jumpux(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,
                            jumpD1jumpuxxcoef,index,r,sk,alpha[r][(sk+1)/2],thesign,
                            normal,tangent,mid,D1ucoef[r][(sk+1)/2],
                            D1uxcoef[r][(sk+1)/2],D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,
                            S,pb,grid);
//            cout << "getting jump D2u" << endl;
            getcim345jumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,
                             alpha[r][(sk+1)/2],thesign,normal,mid,a,D1u[r][(sk+1)/2],
                             D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                             D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,D2u,D2ucoef,D2uxcoef,
                             D2uxxcoef,D2jumpuxxcoef,jumpD1u,jumpD1ucoef,jumpD1uxcoef,
                             jumpD1uxxcoef,jumpD1jumpuxxcoef,S,pb,grid);
// form d0 and dcoef's rth entry 
            d0[grid.dim*(sk+1)/2+r] = thesign*(tau/(grid.dx[r]*grid.dx[r])+
                                               sk*beta*jumpD1u/grid.dx[r]+
                                               0.5*beta*beta*jumpD2u)-
                                      sk*D1u[r][(sk+1)/2][r]/grid.dx[r];
            exactd[grid.dim*(sk+1)/2+r] = d0[grid.dim*(sk+1)/2+r];

            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
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
               exactd[grid.dim*(sk+1)/2+r] += evalarray(dcoef[grid.dim*(sk+1)/2+r],
                                                        sindex)*
                                              getu(tempindex,0,0,0.0,tempsign,grid);
   
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
            for (m = 0; m < grid.dim; m++)
               tempindex[m] = index[m]+sindex[m]-mid;
            if (evalarray(S,tempindex) < 0.0)
               tempsign = -1.0;
            else
               tempsign = 1.0;
            exactd[grid.dim*(sk+1)/2+r] += -1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(tempindex,0,0,0.0,tempsign,grid);
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
            exactd[grid.dim*(sk+1)/2+r] += 1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(tempindex,0,0,0.0,tempsign,grid);
            sindex[r] = mid;
// form matrix's rth row
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m] = -thesign*
                                           (sk*beta*jumpD1uxxcoef[m]/grid.dx[r]+
                                            0.5*beta*beta*jumpD2uxxcoef[m])+
                                           sk*D1uxxcoef[r][(sk+1)/2][r][m]/grid.dx[r];

            G[grid.dim*(sk+1)/2+r][r] += 0.5*(beta*beta-alpha[r][(sk+1)/2]*
                                                        alpha[r][(sk+1)/2]);
            for (m = 0; m < grid.dim; m++)
               G[grid.dim*(sk+1)/2+r][m+grid.dim] = -thesign*
                                                    (sk*beta*jumpD1uxcoef[m]/grid.dx[r]+
                                                     0.5*beta*beta*jumpD2uxcoef[m])+
                                                    sk*D1uxcoef[r][(sk+1)/2][r][m]/
                                                    grid.dx[r];
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
            exactd[grid.dim*(sk+1)/2+r] += -1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(index,0,0,0.0,thesign,grid);
            sindex[r] = mid+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,1.0/(grid.dx[r]*grid.dx[r]));
            exactd[grid.dim*(sk+1)/2+r] += 1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(index,r,sk,1.0,thesign,grid);
            sindex[r] = mid;
         }

// solve for uxx and put in matrix
   gecp0(LU,PLR,PLC,G,2*grid.dim-1,2*grid.dim-1);
   forwardbacksub0(temp,d0,LU,PLR,PLC,2*grid.dim-1);
   value = 0.0;
   for (n = 0; n < grid.dim; n++)
      value += ehere*temp[n];
//   setvalarray(b,index,evalarray(b,index)+value);
   setvalarray(b,index,getf(index,0,0,0.0,thesign,pb,grid)+value);
   for (j = 0; j < grid.dim; j++)
      D2u[j][j] = temp[j];
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
      value = 0.0;
      for (n = 0; n < grid.dim; n++)
         value -= ehere*temp[n];
      if (value != 0.0)
         sparse2(index,tindex,A,value,grid);

      for (n = 0; n < grid.dim; n++)
         setvalarray(D2ucoef[n][n],sindex,temp[n]);
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
   value = evalarray(a,index);
   if (value != 0.0)
      sparse2(index,index,A,value,grid);

// storing uint info
   double tol = 1.0e-14;
   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
         if (gamma[r][(sk+1)/2] == 1)
         {
            newstorage(Dusmall[buildsize]);
            Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
            Dusmall[buildsize].info[1] = r;
            Dusmall[buildsize].info[2] = sk;
            Dusmall[buildsize].info[3] = -1;
            Dusmall[buildsize].mid = mid;
            value = sk*alpha[r][(sk+1)/2]*grid.dx[r]*ux[r];
            if (globintorder == 3)
               value += 0.5*(alpha[r][(sk+1)/2]*grid.dx[r])*
                            (alpha[r][(sk+1)/2]*grid.dx[r])*D2u[r][r];
            sparseorder(-1,Dusmall[buildsize].head,value);
            for (s = 0; s < grid.dim; s++)
               sindex[s] = 0;
            while (sindex[0] <= N)
            {
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
                     for (s = 0; s < grid.dim; s++)
                        rindex[s] = index[s]+sindex[s]-mid;
                     sparseorder(sub2ind(rindex,grid.nx,grid.dim),
                                 Dusmall[buildsize].head,value);
                  }
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
            buildsize++;

// storing Du info
            for (t = 0; t < grid.dim; t++)
            {
               newstorage(Dusmall[buildsize]);
               Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
               Dusmall[buildsize].info[1] = r;
               Dusmall[buildsize].info[2] = sk;
               Dusmall[buildsize].info[3] = t;
               Dusmall[buildsize].mid = mid;

               value = D1u[r][(sk+1)/2][t];
               for (m = 0; m < grid.dim; m++)
                  value += D1uxcoef[r][(sk+1)/2][t][m]*ux[m]+
                           D1uxxcoef[r][(sk+1)/2][t][m]*D2u[m][m];
               if (fabs(value) > tol)
                  sparseorder(-1,Dusmall[buildsize].head,value);

               for (s = 0; s < grid.dim; s++)
                  sindex[s] = 0;
               while (sindex[0] <= N)
               {
                  value = evalarray(D1ucoef[r][(sk+1)/2][t],sindex);
                  for (m = 0; m < grid.dim; m++)
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
                        for (s = 0; s < grid.dim; s++)
                           rindex[s] = index[s]+sindex[s]-mid;
                        sparseorder(sub2ind(rindex,grid.nx,grid.dim),
                                    Dusmall[buildsize].head,value);
                     }
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

               buildsize++;
            }
         }

   free_matrix(alpha,grid.dim-1,1);
   free_matrix(LU,2*grid.dim-1,2*grid.dim-1);
   free_matrix(G,2*grid.dim-1,2*grid.dim-1);
   for (r = 0; r < 2*grid.dim; r++)
      free_matrix(dcoef[r],N,N,N);

   for (r = 0; r < grid.dim; r++)
      for (sk = -1; sk <= 1; sk += 2)
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
   free_matrix(D1jumpuxxcoef,grid.dim-1,grid.dim-1,grid.dim-1);

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

   for (r = 0; r < grid.dim; r++)
      free_matrix(uxcoef[r],N,N,N);
   delete [] uxcoef;
}

/*
used when with globdist = 1, otherwise use [getcim345D2u]

in the first pass, perm=0, rstart = sstar = 0, 
return perm=2 if there is usual mix Du, then second pass is skipped
return perm=1 if have to use the other side, then in the second pass, location for interface is provided,  use rstar and sstar to find jumpuxxcoeff

globdistvar = 1: 
getsk2 with central diff > cim5 (triangle stencil) > getsk2 without central diff > sk2 on out-of-plane nbr point(prefer the side with central diff) > sk2 on the other side
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
   int eindex[grid.dim];
   eindex[0] = 17;
   eindex[1] = 5;
   eindex[2] = 19;

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


// get Du- between index[][][] and index[rstart+sstart][][] 
// in dim-r, Du-[r] = u0[r] +  sum ucoef[r][index] u[index] + sum uxcoef[r][0:2] Du + sum uxxoef[r][0:2] DDu
// uxcoeff[r][i] = Du-_r in terms of Du_i, i = 1,2,3
// uxxoef[r][i] = Du-_r in terms of DDu_i, i = 1,2,3
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
                                           uxxcoeflarge[r][m][n]*
                                           evalarray(D2ucoef[m][n],sindex)); 
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
       cim345(A,b,Dusmall,buildsize,tindex,a,gamma,S,pb,grid);
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
