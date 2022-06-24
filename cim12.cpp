#include "cim12.h"
#include "numerics.h"
#include "interface.h"
#include "helper.h"
#include "finitediff.h"
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

void cim1(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize,
          int *index, int gamma[][2], double ***S, PBData &pb, GridData &grid)
{
   int r, s, j, k, t, m, mid = 1, N = 2*mid;
   int Narray[grid.dim], rindex[grid.dim], sindex[grid.dim];
   for (r = 0; r < grid.dim; r++)
      Narray[r] = N;
   double **c = matrix(grid.dim-1,1),
          ***f = matrix(grid.dim-1,1,grid.dim-1);
   double **LU, **M, value;
   int PLR[grid.dim], PLC[grid.dim];
   double alpha[grid.dim][2], beta, tangent[grid.dim], normal[grid.dim],
          temp[grid.dim];
   double ethere, ehere, ebar;
// ADDING
   double sigma, thesign, tau, Dtau[grid.dim], **d = matrix(grid.dim-1,1);

   LU = matrix(grid.dim-1,grid.dim-1);
   M = matrix(grid.dim-1,grid.dim-1);

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

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (s = -1; s <= 1; s += 2)
   {
      for (r = 0; r < grid.dim; r++)
      {
         for (j = 0; j < grid.dim; j++)
            f[r][(s+1)/2][j] = 0.0;
         c[r][(s+1)/2] = 0.0;
         d[r][(s+1)/2] = 0.0;
      }
      for (r = 0; r < grid.dim; r++)
      {
         rindex[r] = index[r]+s;
         getinterfaceinfo(alpha[r][(s+1)/2],tangent,normal,S,index,rindex,grid);
         beta = 1.0-alpha[r][(s+1)/2];
         ebar = alpha[r][(s+1)/2]*ethere+beta*ehere;
         if (globcim1order == 1)
            for (j = 0; j < grid.dim; j++)
               M[r][j] = -gamma[r][(s+1)/2]*gamma[j][(s+1)/2]*beta*(ehere-ethere)/ebar*
                          tangent[r]*tangent[j];
         else
            for (j = 0; j < grid.dim; j++)
               M[r][j] = -gamma[r][(s+1)/2]*beta*(ehere-ethere)/ebar*tangent[r]*
                          tangent[j];
         M[r][r] += 1.0;

         if (globcim1order == 1)
         {
            for (j = 0; j < grid.dim; j++)
               f[r][(s+1)/2][j] = (gamma[r][(s+1)/2]*beta*(ehere-ethere)*tangent[r]*
                                  (1-gamma[j][(s+1)/2])*s*tangent[j])/(grid.dx[r]*ebar);
            f[r][(s+1)/2][r] += (s*ethere)/(grid.dx[r]*ebar);
            c[r][(s+1)/2] = s*ethere;
            for (j = 0; j < grid.dim; j++)
               c[r][(s+1)/2] += gamma[r][(s+1)/2]*beta*(ehere-ethere)*tangent[r]*
                                (1-gamma[j][(s+1)/2])*s*tangent[j];
         }
         else
         {
            f[r][(s+1)/2][r] += (s*ethere)/(grid.dx[r]*ebar);
            c[r][(s+1)/2] = s*ethere;
         }
         c[r][(s+1)/2] /= -grid.dx[r]*ebar;
// ADDING
         getsigma(sigma,index,r,s,alpha[r][(s+1)/2],normal,pb,grid);
         gettau(tau,index,r,s,alpha[r][(s+1)/2],grid);
         getDtau(Dtau,index,r,s,alpha[r][(s+1)/2],grid);
         d[r][(s+1)/2] = thesign*gamma[r][(s+1)/2]*
                         (beta/ebar*normal[r]*sigma+s*ethere/(grid.dx[r]*ebar)*tau+
                          ethere*beta/ebar*tangent[r]*getdotprod(Dtau,tangent,grid.dim));
         rindex[r] = index[r];
      }
      gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);

      for (r = 0; r < grid.dim; r++)
      {
         for (j = 0; j < grid.dim; j++)
            temp[j] = f[j][(s+1)/2][r];
         forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
         rindex[r] = index[r]+s;
         value = 0.0;
         for (j = 0; j < grid.dim; j++)
            value += -s*ehere*temp[j]/grid.dx[r];
//            value += temp[j];
//         value *= -s*ehere/grid.dx[r];
         sparse2(index,rindex,A,value,grid);
         rindex[r] = index[r];

         for (j = 0; j < grid.dim; j++)
            f[j][(s+1)/2][r] = temp[j];
      }
      for (j = 0; j < grid.dim; j++)
         temp[j] = c[j][(s+1)/2];
      forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
      value = 0.0;
      for (j = 0; j < grid.dim; j++)
         value += temp[j]/grid.dx[j];
      value *= -s*ehere;
      sparse2(index,index,A,value,grid);

      for (j = 0; j < grid.dim; j++)
         c[j][(s+1)/2] = temp[j];
// ADDING
      for (j = 0; j < grid.dim; j++)
         temp[j] = d[j][(s+1)/2];
      forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
      value = 0.0;
      for (j = 0; j < grid.dim; j++)
         value += temp[j]/grid.dx[j];
      value *= s*ehere;
      setvalarray(b,index,evalarray(b,index)+value);
//      setvalarray(b,index,getf(index,0,0,0.0,thesign,pb,grid)+value);

      for (j = 0; j < grid.dim; j++)
         d[j][(s+1)/2] = temp[j];
   }

   if (globsmall == 1)
      for (r = 0; r < grid.dim; r++)
         sindex[r] = mid;
   else if (globsmall == 2)
      for (r = 0; r < grid.dim; r++)
         rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
      for (s = -1; s <= 1; s += 2)
         if (gamma[r][(s+1)/2] == 1)
         {
            newstorage(Dusmall[buildsize]);
            Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
            Dusmall[buildsize].info[1] = r;
            Dusmall[buildsize].info[2] = s;
            Dusmall[buildsize].info[3] = -1;
            Dusmall[buildsize].mid = mid;
            sparseorder(-1,Dusmall[buildsize].head,
                        d[r][(s+1)/2]*s*alpha[r][(s+1)/2]*grid.dx[r]);
            if (globsmall == 1)
               sparseorder(sub2ind(sindex,Narray,grid.dim),Dusmall[buildsize].head,
                           1.0+c[r][(s+1)/2]*s*alpha[r][(s+1)/2]*grid.dx[r]);
            else if (globsmall == 2)
               sparseorder(sub2ind(rindex,grid.nx,grid.dim),Dusmall[buildsize].head,
                           1.0+c[r][(s+1)/2]*s*alpha[r][(s+1)/2]*grid.dx[r]);
            for (m = 0; m < grid.dim; m++)
            {
               if (globsmall == 1)
               {
                  sindex[m] = mid+s;
                  sparseorder(sub2ind(sindex,Narray,grid.dim),Dusmall[buildsize].head,
                              f[r][(s+1)/2][m]*s*alpha[r][(s+1)/2]*grid.dx[r]);
                  sindex[m] = mid;
               }
               else if (globsmall == 2)
               {
                  rindex[m] = index[m]+s;
                  sparseorder(sub2ind(rindex,grid.nx,grid.dim),Dusmall[buildsize].head,
                              f[r][(s+1)/2][m]*s*alpha[r][(s+1)/2]*grid.dx[r]);
                  rindex[m] = index[m];
               }
            }
            buildsize++;

            for (t = 0; t < grid.dim; t++)
            {
               newstorage(Dusmall[buildsize]);
               Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
               Dusmall[buildsize].info[1] = r;
               Dusmall[buildsize].info[2] = s;
               Dusmall[buildsize].info[3] = t;
               Dusmall[buildsize].mid = mid;
               sparseorder(-1,Dusmall[buildsize].head,d[t][(s+1)/2]);
               if (globsmall == 1)
                  sparseorder(sub2ind(sindex,Narray,grid.dim),Dusmall[buildsize].head,
                              c[t][(s+1)/2]);
               else if (globsmall == 2)
                  sparseorder(sub2ind(rindex,grid.nx,grid.dim),Dusmall[buildsize].head,
                              c[t][(s+1)/2]);
               for (m = 0; m < grid.dim; m++)
               {
                  if (globsmall == 1)
                  {
                     sindex[m] = mid+s;
                     sparseorder(sub2ind(sindex,Narray,grid.dim),Dusmall[buildsize].head,
                                 f[t][(s+1)/2][m]);
                     sindex[m] = mid;
                  }
                  else if (globsmall == 2)
                  {
                     rindex[m] = index[m]+s;
                     sparseorder(sub2ind(rindex,grid.nx,grid.dim),
                                 Dusmall[buildsize].head,f[t][(s+1)/2][m]);
                     rindex[m] = index[m];
                  }
               }
               buildsize++;
            }
         }


   free_matrix(M,grid.dim-1,grid.dim-1);
   free_matrix(LU,grid.dim-1,grid.dim-1);
   free_matrix(c,grid.dim-1,1);
   free_matrix(d,grid.dim-1,1);
   free_matrix(f,grid.dim-1,1,grid.dim-1);
}

void cim1(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize,
          int *index, double ***a, int gamma[][2], double ***S, PBData &pb, 
          GridData &grid)
{
   int r, s, j, k, t, m, mid = 1, N = 2*mid;
   int Narray[grid.dim], rindex[grid.dim], sindex[grid.dim];
   for (r = 0; r < grid.dim; r++)
      Narray[r] = N;
   double **c = matrix(grid.dim-1,1),
          ***f = matrix(grid.dim-1,1,grid.dim-1);
   double **LU, **M, value;
   int PLR[grid.dim], PLC[grid.dim];
   double alpha[grid.dim][2], beta, tangent[grid.dim], normal[grid.dim],
          temp[grid.dim];
   double ethere, ehere, ebar;
// ADDING
   double sigma, thesign, tau, Dtau[grid.dim], **d = matrix(grid.dim-1,1);

   LU = matrix(grid.dim-1,grid.dim-1);
   M = matrix(grid.dim-1,grid.dim-1);

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

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (s = -1; s <= 1; s += 2)
   {
      for (r = 0; r < grid.dim; r++)
      {
         for (j = 0; j < grid.dim; j++)
            f[r][(s+1)/2][j] = 0.0;
         c[r][(s+1)/2] = 0.0;
         d[r][(s+1)/2] = 0.0;
      }
      for (r = 0; r < grid.dim; r++)
      {
         rindex[r] = index[r]+s;
         getinterfaceinfo(alpha[r][(s+1)/2],tangent,normal,S,index,rindex,grid);
         beta = 1.0-alpha[r][(s+1)/2];
         ebar = alpha[r][(s+1)/2]*ethere+beta*ehere;
         if (globcim1order == 1)
            for (j = 0; j < grid.dim; j++)
               M[r][j] = -gamma[r][(s+1)/2]*gamma[j][(s+1)/2]*beta*(ehere-ethere)/ebar*
                          tangent[r]*tangent[j];
         else
            for (j = 0; j < grid.dim; j++)
               M[r][j] = -gamma[r][(s+1)/2]*beta*(ehere-ethere)/ebar*tangent[r]*
                          tangent[j];
         M[r][r] += 1.0;

         if (globcim1order == 1)
         {
            for (j = 0; j < grid.dim; j++)
               f[r][(s+1)/2][j] = (gamma[r][(s+1)/2]*beta*(ehere-ethere)*tangent[r]*
                                  (1-gamma[j][(s+1)/2])*s*tangent[j])/(grid.dx[r]*ebar);
            f[r][(s+1)/2][r] += (s*ethere)/(grid.dx[r]*ebar);
            c[r][(s+1)/2] = s*ethere;
            for (j = 0; j < grid.dim; j++)
               c[r][(s+1)/2] += gamma[r][(s+1)/2]*beta*(ehere-ethere)*tangent[r]*
                                (1-gamma[j][(s+1)/2])*s*tangent[j];
         }
         else
         {
            f[r][(s+1)/2][r] += (s*ethere)/(grid.dx[r]*ebar);
            c[r][(s+1)/2] = s*ethere;
         }
         c[r][(s+1)/2] /= -grid.dx[r]*ebar;
// ADDING
         getsigma(sigma,index,r,s,alpha[r][(s+1)/2],normal,pb,grid);
         gettau(tau,index,r,s,alpha[r][(s+1)/2],grid);
         getDtau(Dtau,index,r,s,alpha[r][(s+1)/2],grid);
         d[r][(s+1)/2] = thesign*gamma[r][(s+1)/2]*
                         (beta/ebar*normal[r]*sigma+s*ethere/(grid.dx[r]*ebar)*tau+
                          ethere*beta/ebar*tangent[r]*getdotprod(Dtau,tangent,grid.dim));
         rindex[r] = index[r];
      }
      gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);

      for (r = 0; r < grid.dim; r++)
      {
         for (j = 0; j < grid.dim; j++)
            temp[j] = f[j][(s+1)/2][r];
         forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
         rindex[r] = index[r]+s;
         value = 0.0;
         for (j = 0; j < grid.dim; j++)
            value += -s*ehere*temp[j]/grid.dx[r];
//            value += temp[j];
//         value *= -s*ehere/grid.dx[r];
         sparse2(index,rindex,A,value,grid);
         rindex[r] = index[r];

         for (j = 0; j < grid.dim; j++)
            f[j][(s+1)/2][r] = temp[j];
      }
      for (j = 0; j < grid.dim; j++)
         temp[j] = c[j][(s+1)/2];
      forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
      value = 0.0;
      for (j = 0; j < grid.dim; j++)
         value += temp[j]/grid.dx[j];
      value *= -s*ehere;
      sparse2(index,index,A,value,grid);

      for (j = 0; j < grid.dim; j++)
         c[j][(s+1)/2] = temp[j];
// ADDING
      for (j = 0; j < grid.dim; j++)
         temp[j] = d[j][(s+1)/2];
      forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
      value = 0.0;
      for (j = 0; j < grid.dim; j++)
         value += temp[j]/grid.dx[j];
      value *= s*ehere;
      setvalarray(b,index,evalarray(b,index)+value);
//      setvalarray(b,index,getf(index,0,0,0.0,thesign,pb,grid)+value);

      for (j = 0; j < grid.dim; j++)
         d[j][(s+1)/2] = temp[j];
   }
// addition of a
   sparse2(index,index,A,evalarray(a,index),grid);

   if (globsmall == 1)
      for (r = 0; r < grid.dim; r++)
         sindex[r] = mid;
   else if (globsmall == 2)
      for (r = 0; r < grid.dim; r++)
         rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
      for (s = -1; s <= 1; s += 2)
         if (gamma[r][(s+1)/2] == 1)
         {
            newstorage(Dusmall[buildsize]);
            Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
            Dusmall[buildsize].info[1] = r;
            Dusmall[buildsize].info[2] = s;
            Dusmall[buildsize].info[3] = -1;
            Dusmall[buildsize].mid = mid;
            sparseorder(-1,Dusmall[buildsize].head,
                        d[r][(s+1)/2]*s*alpha[r][(s+1)/2]*grid.dx[r]);
            if (globsmall == 1)
               sparseorder(sub2ind(sindex,Narray,grid.dim),Dusmall[buildsize].head,
                           1.0+c[r][(s+1)/2]*s*alpha[r][(s+1)/2]*grid.dx[r]);
            else if (globsmall == 2)
               sparseorder(sub2ind(rindex,grid.nx,grid.dim),Dusmall[buildsize].head,
                           1.0+c[r][(s+1)/2]*s*alpha[r][(s+1)/2]*grid.dx[r]);
            for (m = 0; m < grid.dim; m++)
            {
               if (globsmall == 1)
               {
                  sindex[m] = mid+s;
                  sparseorder(sub2ind(sindex,Narray,grid.dim),Dusmall[buildsize].head,
                              f[r][(s+1)/2][m]*s*alpha[r][(s+1)/2]*grid.dx[r]);
                  sindex[m] = mid;
               }
               else if (globsmall == 2)
               {
                  rindex[m] = index[m]+s;
                  sparseorder(sub2ind(rindex,grid.nx,grid.dim),Dusmall[buildsize].head,
                              f[r][(s+1)/2][m]*s*alpha[r][(s+1)/2]*grid.dx[r]);
                  rindex[m] = index[m];
               }
            }
            buildsize++;

            for (t = 0; t < grid.dim; t++)
            {
               newstorage(Dusmall[buildsize]);
               Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
               Dusmall[buildsize].info[1] = r;
               Dusmall[buildsize].info[2] = s;
               Dusmall[buildsize].info[3] = t;
               Dusmall[buildsize].mid = mid;
               sparseorder(-1,Dusmall[buildsize].head,d[t][(s+1)/2]);
               if (globsmall == 1)
                  sparseorder(sub2ind(sindex,Narray,grid.dim),Dusmall[buildsize].head,
                              c[t][(s+1)/2]);
               else if (globsmall == 2)
                  sparseorder(sub2ind(rindex,grid.nx,grid.dim),Dusmall[buildsize].head,
                              c[t][(s+1)/2]);
               for (m = 0; m < grid.dim; m++)
               {
                  if (globsmall == 1)
                  {
                     sindex[m] = mid+s;
                     sparseorder(sub2ind(sindex,Narray,grid.dim),Dusmall[buildsize].head,
                                 f[t][(s+1)/2][m]);
                     sindex[m] = mid;
                  }
                  else if (globsmall == 2)
                  {
                     rindex[m] = index[m]+s;
                     sparseorder(sub2ind(rindex,grid.nx,grid.dim),
                                 Dusmall[buildsize].head,f[t][(s+1)/2][m]);
                     rindex[m] = index[m];
                  }
               }
               buildsize++;
            }
         }


   free_matrix(M,grid.dim-1,grid.dim-1);
   free_matrix(LU,grid.dim-1,grid.dim-1);
   free_matrix(c,grid.dim-1,1);
   free_matrix(d,grid.dim-1,1);
   free_matrix(f,grid.dim-1,1,grid.dim-1);
}

void cim1again3(SparseElt2**** &A, double ***b, int *index, int gamma[][2], double ***S, 
                PBData &pb, GridData &grid)
{
// uses 2nd order central differencing when possible
   int r, s, j, k;
   int rindex[grid.dim];
   double c[grid.dim], f[grid.dim][grid.dim][2];
   double **LU, **M, value; 
   int PLR[grid.dim], PLC[grid.dim];
   double alpha, beta, tangent[grid.dim], normal[grid.dim], temp[grid.dim];
   double ethere, ehere, ebar;
// ADDING
   double sigma, d[grid.dim], thesign, tau, Dtau[grid.dim];

   LU = matrix(grid.dim-1,grid.dim-1);
   M = matrix(grid.dim-1,grid.dim-1);

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

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (s = -1; s <= 1; s += 2)
   {
      for (r = 0; r < grid.dim; r++)
      {
         for (j = 0; j < grid.dim; j++)
            for (k = 0; k <= 1; k++)
               f[r][j][k] = 0.0;
         c[r] = 0.0;
         d[r] = 0.0;
      }
      for (r = 0; r < grid.dim; r++)
      {
         rindex[r] = index[r]+s;
         getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
         beta = 1.0-alpha;
         ebar = alpha*ethere+beta*ehere;
         for (j = 0; j < grid.dim; j++)
            M[r][j] = -gamma[r][(s+1)/2]*gamma[j][(s+1)/2]*beta*(ehere-ethere)/ebar*
                       tangent[r]*tangent[j];
         M[r][r] += 1.0;
            
         for (j = 0; j < grid.dim; j++)
            if (gamma[j][0] == 0 && gamma[j][1] == 0)
               for (k = 0; k <= 1; k++)
                  f[r][j][k] = (gamma[r][(s+1)/2]*beta*(ehere-ethere)*tangent[r]*
                               (1-gamma[j][k])*(2*k-1)*tangent[j])/(2.0*grid.dx[r]*ebar);
            else
               f[r][j][(s+1)/2] = (gamma[r][(s+1)/2]*beta*(ehere-ethere)*tangent[r]*
                                  (1-gamma[j][(s+1)/2])*s*tangent[j])/(grid.dx[r]*ebar);
         f[r][r][(s+1)/2] += (s*ethere)/(grid.dx[r]*ebar);
         c[r] = s*ethere;
         for (j = 0; j < grid.dim; j++)
            if (gamma[j][0] == 0 && gamma[j][1] == 0);
            else
               c[r] += gamma[r][(s+1)/2]*beta*(ehere-ethere)*tangent[r]*
                       (1-gamma[j][(s+1)/2])*s*tangent[j];
         c[r] /= -grid.dx[r]*ebar;
// ADDING
         getsigma(sigma,index,r,s,alpha,normal,pb,grid);
         gettau(tau,index,r,s,alpha,grid);
         getDtau(Dtau,index,r,s,alpha,grid);
//         d[r] = thesign*beta/ebar*normal[r]*sigma;
         d[r] = thesign*gamma[r][(s+1)/2]*
                (beta/ebar*normal[r]*sigma+s*ethere/(grid.dx[r]*ebar)*tau+
                 ethere*beta/ebar*tangent[r]*getdotprod(Dtau,tangent,grid.dim));
         rindex[r] = index[r];
      }
      gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);
      for (k = 0; k <= 1; k++)
         for (r = 0; r < grid.dim; r++)
         {
            for (j = 0; j < grid.dim; j++)
               temp[j] = f[j][r][k];
            forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
            rindex[r] = index[r]+2*k-1;
            value = 0.0;
            for (j = 0; j < grid.dim; j++)
               value += temp[j];
            value *= -s*ehere/grid.dx[r];
            sparse2(index,rindex,A,value,grid);
            rindex[r] = index[r];
         }
      forwardbacksub0(c,c,LU,PLR,PLC,grid.dim-1);
      value = 0.0;
      for (j = 0; j < grid.dim; j++)
         value += c[j]/grid.dx[j];
      value *= -s*ehere;
      sparse2(index,index,A,value,grid);
// ADDING
      forwardbacksub0(d,d,LU,PLR,PLC,grid.dim-1);
      value = 0.0;
      for (j = 0; j < grid.dim; j++)
         value += d[j]/grid.dx[j];
      value *= s*ehere;
      setvalarray(b,index,evalarray(b,index)+value);
//      setvalarray(b,index,getf(index,0,0,0.0,thesign,pb,grid)+value);
   }
   
   free_matrix(M,grid.dim-1,grid.dim-1);
   free_matrix(LU,grid.dim-1,grid.dim-1);
}



//cim2 with Dusmall
void cim2(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize,
          int *index, int gamma[][2], double ***S, PBData &pb, GridData &grid)
{
   int r, s, i, j, k, tmps, t, m, mid = 2, N = 2*mid;
   int rindex[grid.dim], tindex[grid.dim], Narray[grid.dim];
   for (r = 0; r < grid.dim; r++)
      Narray[r] = N;
   double ****f;
   double **LU, **M, value, Sval;
   int PLR[grid.dim], PLC[grid.dim];
   double alpha[grid.dim], beta, tangent[grid.dim], normal[grid.dim];
   double bk, rthere, rhere, ethere, ehere, ehat, C;
   double temp[grid.dim], a[4];
   int sk[grid.dim];
// ADDING
   double sigma, d[grid.dim], thesign, tau, Dtau[grid.dim];
   double uxxcoef[grid.dim][grid.dim]; // uxxcoeff[r][j] is coeff of uxx for Du-[j] at x_r (interface along r)
   double *****Du, ****uint;

   LU = matrix(grid.dim-1,grid.dim-1);
   M = matrix(grid.dim-1,grid.dim-1);
   f = matrix(grid.dim-1,N,N,N); // rhs of coupling linear system,  r-the row M urr =  d[r] + f[r][idx] u[idx]

   for (r = 0; r < grid.dim; r++)
      sk[r] = gamma[r][1]-gamma[r][0]; // 1 if one interface along dim r, 0  if no interface or two interface

   uint = matrix(grid.dim-1,N,N,N);
   Du = new double ****[grid.dim]; // Du[r][j] is coef of u derivative in j dim at x_r ( interface along r)
   for (r = 0; r < grid.dim; r++)
      Du[r] = matrix(grid.dim-1,N,N,N);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= N)
   {
      for (r = 0; r < grid.dim; r++)
      {
         setvalarray(uint[r],tindex,0.0);
         for (m = 0; m < grid.dim; m++)
            setvalarray(Du[r][m],tindex,0.0);
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > N; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   Sval = evalarray(S,index);
   if (Sval < 0.0)
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

   for (r = 0; r < grid.dim; r++)
   {
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= N)
      {
         setvalarray(f[r],tindex,0.0);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > N; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

// ADDING
      d[r] = 0.0;
   }

   #if 1
   if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
        print_surf(S, index, 2);    
      }
   #endif

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
      tindex[r] = mid;
   for (r = 0; r < grid.dim; r++)
   {
      rindex[r] = index[r]+sk[r];
      getinterfaceinfo(alpha[r],tangent,normal,S,index,rindex,grid);
      rindex[r] = index[r];
      if (abs(sk[r]) == 1) //has only 1 interface in dim r
      {
         beta = 1.0-alpha[r];
         ehat = (beta+beta*beta)*(0.5+alpha[r])*ehere+
                (alpha[r]+alpha[r]*alpha[r])*(0.5+beta)*ethere;
         rhere = ehere/ehat;
         rthere = ethere/ehat;
         a[0] = (beta+beta*beta)*rhere+alpha[r]*(1.0+2.0*beta)*rthere; //a_{i,-1} to a_{i,2} eqn(10) cim paper
         a[1] = -(beta+beta*beta)*rhere-(1.0+alpha[r])*(1.0+2.0*beta)*rthere;
         a[2] = (1.0+beta)*(1.0+beta)*rthere;
         a[3] = -beta*beta*rthere;
         bk = -(beta+beta*beta)*(rthere-rhere);
         C = bk*tangent[r]; // C = b_k (t_k . e_k), eqn(34), coeff for grau(u-)
         for (s = 0; s < 4; s++)
         {
            tindex[r] = mid+(s-1)*sk[r];
            setvalarray(f[r],tindex,evalarray(f[r],tindex)+a[s]); //u coeff from L operator eqn(32)
         }
         tindex[r] = mid-sk[r];
         setvalarray(f[r],tindex,evalarray(f[r],tindex)-sk[r]*C*tangent[r]*sk[r]);
         setvalarray(Du[r][r],tindex,evalarray(Du[r][r],tindex)-sk[r]); // eqn(37), eqn(30) coeff of u(x + s_k h e_k)
         tindex[r] = mid;
         setvalarray(f[r],tindex,evalarray(f[r],tindex)+sk[r]*C*tangent[r]*sk[r]);
         setvalarray(Du[r][r],tindex,evalarray(Du[r][r],tindex)+sk[r]); //eqn(37), eqn(30) coeff of u(x)

         M[r][r] = 1.0-abs(sk[r])*(0.5+alpha[r])*bk*tangent[r]*tangent[r];
         uxxcoef[r][r] = sk[r]*(0.5+alpha[r]); //eqn(37) case 1

         for (j = 0; j < grid.dim; j++)
            if (j != r) //look at dimension other than r
            {
               tmps = sk[j];
               if (tmps == 0)// both sides have interface or both sides do not have
                  for (s = -1; s <= 1 && tmps == 0; s += 2)
                  {
                     rindex[j] = index[j]+s;
                     if (evalarray(S,rindex)*Sval <= 0.0)
                        tmps = s;// find s on the other side,
                  }
               if (tmps == 0)//can not find such s, both sides do not have interface, take a step back along r
               {
                  rindex[r] = index[r]-sk[r];
                  for (s = -1; s <= 1 && tmps == 0; s += 2)
                  {
                     rindex[j] = index[j]+s;
                     if (evalarray(S,rindex)*Sval <= 0.0)
                        tmps = s;
                  }
               }
               rindex[r] = index[r];
               rindex[j] = index[j];

               M[r][j] = -0.5*tmps*sk[r]*bk*tangent[j]*tangent[r];
               uxxcoef[r][j] = 0.5*tmps; //eqn(37) case two

               if (abs(tmps) == 1) //eqn(39) operator (D^{s_j=+/-1}_j). Approx Du[r]- in dim j with forward/backward diff in dim j with index and index[r]-1. tmps = s_j
               {
                  // dim j, step back from interface, on u(x - s_j h e_j)
                  tindex[j] = mid-tmps;
                  rindex[j] = index[j]-tmps;
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)-
                                          sk[r]*C*tangent[j]*(1.0+alpha[r])*tmps); // eqn(39), eqn(30)
                  setvalarray(Du[r][j],tindex,evalarray(Du[r][j],tindex)-
                                              tmps*(1.0+alpha[r])); // eqn(37), eqn(30), coeff is  -s_j (1 + alpha_k)
                  // dim r and j, step back from interface, on u(x - s_j h e_j - s_k h e_k)
                  tindex[r] = mid-sk[r];
                  rindex[r] = index[r]-sk[r];
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)+
                                          sk[r]*C*tangent[j]*alpha[r]*tmps);
                  setvalarray(Du[r][j],tindex,evalarray(Du[r][j],tindex)+tmps*alpha[r]); // eqn(37), eqn(30), coef is s_j alpha_k

                  // dim r step back from interface, on u(x - s_k h e_k)
                  tindex[j] = mid;
                  rindex[j] = index[j];
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)-
                                          sk[r]*C*tangent[j]*alpha[r]*tmps);
                  setvalarray(Du[r][j],tindex,evalarray(Du[r][j],tindex)-tmps*alpha[r]);// eqn(37), eqn(30), coef is - s_j alpha_k

                  // back to index on u(x)
                  tindex[r] = mid;
                  rindex[r] = index[r];
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)+
                                          sk[r]*C*tangent[j]*(1.0+alpha[r])*tmps);
                  setvalarray(Du[r][j],tindex,evalarray(Du[r][j],tindex)+
                                              tmps*(1.0+alpha[r])); // eqn(37), eqn(30), coef is s_j (1 + alpha_k)
               }
               else // eqn(39) operator (D^{s_j=0}_j). Approx Du[r]- in dim j with central diff in dim j with index and index[r]-1
               {
                  for (s = -1; s <= 1; s += 2)
                  {
                     tindex[j] = mid+s;
                     rindex[j] = index[j]+s;
                     setvalarray(f[r],tindex,evalarray(f[r],tindex)+
                                             sk[r]*C*tangent[j]*(1.0+alpha[r])*0.5*s);
                     setvalarray(Du[r][j],tindex,evalarray(Du[r][j],tindex)+
                                                 s*(1.0+alpha[r])*0.5);
                  }
                  tindex[r] = mid-sk[r];
                  rindex[r] = index[r]-sk[r];
                  for (s = -1; s <= 1; s += 2)
                  {
                     tindex[j] = mid+s;
                     rindex[j] = index[j]+s;
                     setvalarray(f[r],tindex,evalarray(f[r],tindex)-
                                             sk[r]*C*tangent[j]*alpha[r]*0.5*s);
                     setvalarray(Du[r][j],tindex,evalarray(Du[r][j],tindex)-
                                                 s*alpha[r]*0.5);
                  }
                  tindex[r] = mid;
                  tindex[j] = mid;
                  rindex[r] = index[r];
                  rindex[j] = index[j];
               }
            }

// ADDING
         getsigma(sigma,index,r,sk[r],alpha[r],normal,pb,grid);
         gettau(tau,index,r,sk[r],alpha[r],grid);
         getDtau(Dtau,index,r,sk[r],alpha[r],grid);
//         d[r] = thesign*sk[r]*(beta+beta*beta)*grid.dx[r]/ehat*normal[r]*sigma;
         d[r] = thesign*(abs(sk[r])*(1.0+2.0*beta)*rthere*tau+
                         sk[r]*(beta+beta*beta)*grid.dx[r]*
                         (sigma/ehat*normal[r]+
                          rthere*getdotprod(Dtau,tangent,grid.dim)*tangent[r]));// J_k in eqn(36) cim paper
      }
      else
      {
         for (j = 0; j < grid.dim; j++)
            if (j == r)
               M[r][r] = 1.0;
            else
               M[r][j] = 0.0;
         for (s = -1; s <= 1; s += 2)
         {
            tindex[r] = mid+s;
            setvalarray(f[r],tindex,evalarray(f[r],tindex)+1.0);
         }
         tindex[r] = mid;
         setvalarray(f[r],tindex,evalarray(f[r],tindex)-2.0);
      }
   }

   gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);
   #if 1
   
  if (index[0] == eindex[0] && index[1] == eindex[1] && index[2] == eindex[2]){
    printf("\nAt (%i,%i,%i)\n",index[0],index[1],index[2]);
    double v[3] = {0};
    double tempuxcoef[3] = {0};
    double tempuxxcoef[3] = {0};
    for(int j = 0; j < 3; j++){
      printf("d[%i] = %.6f\n ", j, d[j]);
      for( int iloc = 0; iloc <= N; iloc ++){
        for( int jloc = 0; jloc <= N; jloc ++){
          for( int kloc = 0; kloc <= N; kloc ++){
            if (abs(f[j][iloc][jloc][kloc])>1e-6){
              printf("%.6f  u(%i,%i,%i)\n", f[j][iloc][jloc][kloc],iloc,jloc,kloc);
            }
          }
        }
      }
      //double evalcoef(double u0, double ***ucoef, double *uxcoef, double *uxxcoef, int *index, int rstar, int sstar, double alpha, int mid, double thesign, GridData grid)
      v[j] = 1/pow(grid.dx[0],2) * evalcoef(d[j],f[j], tempuxcoef, tempuxxcoef, index, 0, 0, 0, mid, thesign, grid );
    }

    
    cout<<"M matrix"<<endl;
    printMat(3,3,M);
    cout<<"exact rhs"<<endl;
    printMat(3,v);
    forwardbacksub0(v,v,LU,PLR,PLC,grid.dim-1);
    

    double exactD2u[3] = {0};
    for(int j = 0; j < 3; j++){
      exactD2u[j] = getD2u(index,j,j,0,0,0.0,thesign,grid);  
    }
    cout<<"solved, exact, error"<<endl;
    for(int j = 0; j < 3; j++){
      cout<<std::setw(12)<< v[j]<<", "<<std::setw(12)<<exactD2u[j]<<", "<<std::setw(12)<<v[j] - exactD2u[j]<<endl;
    }
    cout <<" -eps lap(u) - f = " <<- ehere * (v[0]+v[1]+v[2])-getf(index,0,0,0.0,thesign,pb,grid)<<endl;
  }
   // end of local truncation error
  #endif


   // Setup A and b
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= N)
   {
      for (r = 0; r < grid.dim; r++)
         temp[r] = evalarray(f[r],tindex); //temp = f[:][tindex]
      forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1); // temp[1:3]= [uxx uyy uzz] in term of u[index] = inv(M) f[:][tindex] 
      value = 0.0;
      for (r = 0; r < grid.dim; r++)
         value += -ehere*temp[r]/(grid.dx[r]*grid.dx[r]); // - ehere lap(u) = - ehere sum(inv(M)f[:][tindex])
      if (value != 0.0)
      {
         for (r = 0; r < grid.dim; r++)
            rindex[r] = index[r]-mid+tindex[r];
         sparse2(index,rindex,A,value,grid);
      }

      for (r = 0; r < grid.dim; r++)
         setvalarray(f[r],tindex,temp[r]); // f <- inv(M) f

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > N; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

// ADDING
   forwardbacksub0(d,d,LU,PLR,PLC,grid.dim-1);
   value = 0.0;
   for (r = 0; r < grid.dim; r++)
      value += ehere*d[r]/(grid.dx[r]*grid.dx[r]);
   setvalarray(b,index,evalarray(b,index)+value);

// Save Dusmall
   double tol = 1.0e-14;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = mid;
   for (r = 0; r < grid.dim; r++)
   {
      setvalarray(uint[r],tindex,1.0+alpha[r]); // eqn(7), u-(x_r) = (1 + alpha) u_i - alpha u_i-1 + 1/2 h^2 alpha (1 + alpha ) u_i''
      tindex[r] = mid-sk[r];
      setvalarray(uint[r],tindex,-alpha[r]);
      tindex[r] = mid;
   }
   for (r = 0; r < grid.dim; r++)
      if (abs(sk[r]) == 1)
      {
         s = sk[r];
         newstorage(Dusmall[buildsize]);// store coeff to reconstruct u at interface
         Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
         Dusmall[buildsize].info[1] = r;
         Dusmall[buildsize].info[2] = s;
         Dusmall[buildsize].info[3] = -1; // indicator for u at interface
         Dusmall[buildsize].mid = mid;
         value = 0.5*grid.dx[r]*grid.dx[r]*d[r]*alpha[r]*(1.0+alpha[r]); // constant part, 1/2 h^2 alpha (1 + alpha ) inv(M) f
         if (fabs(value) > tol)
            sparseorder(-1,Dusmall[buildsize].head,value);
         for (i = 0; i < grid.dim; i++)
            tindex[i] = 0;
         while (tindex[0] <= N)
         {
            value = evalarray(uint[r],tindex)+
                    0.5*grid.dx[r]*grid.dx[r]*evalarray(f[r],tindex)*alpha[r]*
                    (1.0+alpha[r]); // u part, 
            if (fabs(value) > tol)
            {
               if (globsmall == 1)
                  sparseorder(sub2ind(tindex,Narray,grid.dim),Dusmall[buildsize].head,
                              value);
               else if (globsmall == 2)
               {
                  for (i = 0; i < grid.dim; i++)
                     rindex[i] = index[i]+tindex[i]-mid;
                  sparseorder(sub2ind(rindex,grid.nx,grid.dim),Dusmall[buildsize].head,
                              value);
               }
            }

            (tindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && tindex[i] > N; i--)
            {
               tindex[i] = 0;
               (tindex[i-1])++;
            }
         }
         buildsize++;

         for (t = 0; t < grid.dim; t++)
         {
            newstorage(Dusmall[buildsize]);
            Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
            Dusmall[buildsize].info[1] = r;
            Dusmall[buildsize].info[2] = s;
            Dusmall[buildsize].info[3] = t; //indicator for Du[t]
            Dusmall[buildsize].mid = mid;
            value = uxxcoef[r][t]*d[t]/grid.dx[r]; // coeff in eqn(37) inv(M) d / h 
            if (fabs(value) > tol)
               sparseorder(-1,Dusmall[buildsize].head,value);

            for (i = 0; i < grid.dim; i++)
               tindex[i] = 0;
            while (tindex[0] <= N)
            {
               value = (evalarray(Du[r][t],tindex)+
                        uxxcoef[r][t]*evalarray(f[t],tindex))/grid.dx[r]; // 
               if (fabs(value) > tol)
               {
                  if (globsmall == 1)
                     sparseorder(sub2ind(tindex,Narray,grid.dim),Dusmall[buildsize].head,
                                 value);
                  else if (globsmall == 2)
                  {
                     for (i = 0; i < grid.dim; i++)
                        rindex[i] = index[i]+tindex[i]-mid;
                     sparseorder(sub2ind(rindex,grid.nx,grid.dim),
                                 Dusmall[buildsize].head,value);
                  }
               }

               (tindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && tindex[i] > N; i--)
               {
                  tindex[i] = 0;
                  (tindex[i-1])++;
               }
            }
            buildsize++;
         }
      }

   free_matrix(M,grid.dim-1,grid.dim-1);
   free_matrix(LU,grid.dim-1,grid.dim-1);
   free_matrix(f,grid.dim-1,4,4,4);
   for (r = 0; r < grid.dim; r++)
      free_matrix(Du[r],grid.dim-1,N,N,N);
   delete [] Du;
   free_matrix(uint,grid.dim-1,4,4,4);
}
// Only difference is  ***aa
void cim2(SparseElt2**** &A, double ***b, StorageStruct* &Dusmall, int &buildsize,
          int *index, double ***aa, int gamma[][2], double ***S, PBData &pb, 
          GridData &grid)
{
   int r, s, i, j, k, tmps, t, m, mid = 2, N = 2*mid;
   int rindex[grid.dim], tindex[grid.dim], Narray[grid.dim];
   for (r = 0; r < grid.dim; r++)
      Narray[r] = N;
   double ****f;
   double **LU, **M, value, Sval;
   int PLR[grid.dim], PLC[grid.dim];
   double alpha[grid.dim], beta, tangent[grid.dim], normal[grid.dim];
   double bk, rthere, rhere, ethere, ehere, ehat, C;
   double temp[grid.dim], a[4];
   int sk[grid.dim];
// ADDING
   double sigma, d[grid.dim], thesign, tau, Dtau[grid.dim];
   double uxxcoef[grid.dim][grid.dim];
   double *****Du, ****uint;

   LU = matrix(grid.dim-1,grid.dim-1);
   M = matrix(grid.dim-1,grid.dim-1);
   f = matrix(grid.dim-1,N,N,N);

   for (r = 0; r < grid.dim; r++)
      sk[r] = gamma[r][1]-gamma[r][0];

   uint = matrix(grid.dim-1,N,N,N);
   Du = new double ****[grid.dim];
   for (r = 0; r < grid.dim; r++)
      Du[r] = matrix(grid.dim-1,N,N,N);
   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= N)
   {
      for (r = 0; r < grid.dim; r++)
      {
         setvalarray(uint[r],tindex,0.0);
         for (m = 0; m < grid.dim; m++)
            setvalarray(Du[r][m],tindex,0.0);
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > N; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   Sval = evalarray(S,index);
   if (Sval < 0.0)
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

   for (r = 0; r < grid.dim; r++)
   {
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] <= N)
      {
         setvalarray(f[r],tindex,0.0);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] > N; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

// ADDING
      d[r] = 0.0;
   }

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
      tindex[r] = mid;
   for (r = 0; r < grid.dim; r++)
   {
      rindex[r] = index[r]+sk[r];
      getinterfaceinfo(alpha[r],tangent,normal,S,index,rindex,grid);
      rindex[r] = index[r];
      if (abs(sk[r]) == 1)
      {
         beta = 1.0-alpha[r];
         ehat = (beta+beta*beta)*(0.5+alpha[r])*ehere+
                (alpha[r]+alpha[r]*alpha[r])*(0.5+beta)*ethere;
         rhere = ehere/ehat;
         rthere = ethere/ehat;
         a[0] = (beta+beta*beta)*rhere+alpha[r]*(1.0+2.0*beta)*rthere;
         a[1] = -(beta+beta*beta)*rhere-(1.0+alpha[r])*(1.0+2.0*beta)*rthere;
         a[2] = (1.0+beta)*(1.0+beta)*rthere;
         a[3] = -beta*beta*rthere;
         bk = -(beta+beta*beta)*(rthere-rhere);
         C = bk*tangent[r];
         for (s = 0; s < 4; s++)
         {
            tindex[r] = mid+(s-1)*sk[r];
            setvalarray(f[r],tindex,evalarray(f[r],tindex)+a[s]);
         }
         tindex[r] = mid-sk[r];
         setvalarray(f[r],tindex,evalarray(f[r],tindex)-sk[r]*C*tangent[r]*sk[r]);
         setvalarray(Du[r][r],tindex,evalarray(Du[r][r],tindex)-sk[r]);
         tindex[r] = mid;
         setvalarray(f[r],tindex,evalarray(f[r],tindex)+sk[r]*C*tangent[r]*sk[r]);
         setvalarray(Du[r][r],tindex,evalarray(Du[r][r],tindex)+sk[r]);

         M[r][r] = 1.0-abs(sk[r])*(0.5+alpha[r])*bk*tangent[r]*tangent[r];
         uxxcoef[r][r] = sk[r]*(0.5+alpha[r]);

         for (j = 0; j < grid.dim; j++)
            if (j != r)
            {
               tmps = sk[j];
               if (tmps == 0)
                  for (s = -1; s <= 1 && tmps == 0; s += 2)
                  {
                     rindex[j] = index[j]+s;
                     if (evalarray(S,rindex)*Sval <= 0.0)
                        tmps = s;
                  }
               if (tmps == 0)
               {
                  rindex[r] = index[r]-sk[r];
                  for (s = -1; s <= 1 && tmps == 0; s += 2)
                  {
                     rindex[j] = index[j]+s;
                     if (evalarray(S,rindex)*Sval <= 0.0)
                        tmps = s;
                  }
               }
               rindex[r] = index[r];
               rindex[j] = index[j];

               M[r][j] = -0.5*tmps*sk[r]*bk*tangent[j]*tangent[r];
               uxxcoef[r][j] = 0.5*tmps;

               if (abs(tmps) == 1)
               {
                  tindex[j] = mid-tmps;
                  rindex[j] = index[j]-tmps;
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)-
                                          sk[r]*C*tangent[j]*(1.0+alpha[r])*tmps);
                  setvalarray(Du[r][j],tindex,evalarray(Du[r][j],tindex)-
                                              tmps*(1.0+alpha[r]));
                  tindex[r] = mid-sk[r];
                  rindex[r] = index[r]-sk[r];
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)+
                                          sk[r]*C*tangent[j]*alpha[r]*tmps);
                  setvalarray(Du[r][j],tindex,evalarray(Du[r][j],tindex)+tmps*alpha[r]);
                  tindex[j] = mid;
                  rindex[j] = index[j];
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)-
                                          sk[r]*C*tangent[j]*alpha[r]*tmps);
                  setvalarray(Du[r][j],tindex,evalarray(Du[r][j],tindex)-tmps*alpha[r]);
                  tindex[r] = mid;
                  rindex[r] = index[r];
                  setvalarray(f[r],tindex,evalarray(f[r],tindex)+
                                          sk[r]*C*tangent[j]*(1.0+alpha[r])*tmps);
                  setvalarray(Du[r][j],tindex,evalarray(Du[r][j],tindex)+
                                              tmps*(1.0+alpha[r]));
               }
               else
               {
                  for (s = -1; s <= 1; s += 2)
                  {
                     tindex[j] = mid+s;
                     rindex[j] = index[j]+s;
                     setvalarray(f[r],tindex,evalarray(f[r],tindex)+
                                             sk[r]*C*tangent[j]*(1.0+alpha[r])*0.5*s);
                     setvalarray(Du[r][j],tindex,evalarray(Du[r][j],tindex)+
                                                 s*(1.0+alpha[r])*0.5);
                  }
                  tindex[r] = mid-sk[r];
                  rindex[r] = index[r]-sk[r];
                  for (s = -1; s <= 1; s += 2)
                  {
                     tindex[j] = mid+s;
                     rindex[j] = index[j]+s;
                     setvalarray(f[r],tindex,evalarray(f[r],tindex)-
                                             sk[r]*C*tangent[j]*alpha[r]*0.5*s);
                     setvalarray(Du[r][j],tindex,evalarray(Du[r][j],tindex)-
                                                 s*alpha[r]*0.5);
                  }
                  tindex[r] = mid;
                  tindex[j] = mid;
                  rindex[r] = index[r];
                  rindex[j] = index[j];
               }
            }

// ADDING
         getsigma(sigma,index,r,sk[r],alpha[r],normal,pb,grid);
         gettau(tau,index,r,sk[r],alpha[r],grid);
         getDtau(Dtau,index,r,sk[r],alpha[r],grid);
//         d[r] = thesign*sk[r]*(beta+beta*beta)*grid.dx[r]/ehat*normal[r]*sigma;
         d[r] = thesign*(abs(sk[r])*(1.0+2.0*beta)*rthere*tau+
                         sk[r]*(beta+beta*beta)*grid.dx[r]*
                         (sigma/ehat*normal[r]+
                          rthere*getdotprod(Dtau,tangent,grid.dim)*tangent[r]));
      }
      else
      {
         for (j = 0; j < grid.dim; j++)
            if (j == r)
               M[r][r] = 1.0;
            else
               M[r][j] = 0.0;
         for (s = -1; s <= 1; s += 2)
         {
            tindex[r] = mid+s;
            setvalarray(f[r],tindex,evalarray(f[r],tindex)+1.0);
         }
         tindex[r] = mid;
         setvalarray(f[r],tindex,evalarray(f[r],tindex)-2.0);
      }
   }

   gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= N)
   {
      for (r = 0; r < grid.dim; r++)
         temp[r] = evalarray(f[r],tindex);
      forwardbacksub0(temp,temp,LU,PLR,PLC,grid.dim-1);
      value = 0.0;
      for (r = 0; r < grid.dim; r++)
         value += -ehere*temp[r]/(grid.dx[r]*grid.dx[r]);
      if (value != 0.0)
      {
         for (r = 0; r < grid.dim; r++)
            rindex[r] = index[r]-mid+tindex[r];
         sparse2(index,rindex,A,value,grid);
      }

      for (r = 0; r < grid.dim; r++)
         setvalarray(f[r],tindex,temp[r]);

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > N; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }
// addition of a
   sparse2(index,index,A,evalarray(aa,index),grid);

// ADDING
   forwardbacksub0(d,d,LU,PLR,PLC,grid.dim-1);
   value = 0.0;
   for (r = 0; r < grid.dim; r++)
      value += ehere*d[r]/(grid.dx[r]*grid.dx[r]);
   setvalarray(b,index,evalarray(b,index)+value);

   double tol = 1.0e-14;
//   double tol = 0.0;
   for (i = 0; i < grid.dim; i++)
      tindex[i] = mid;
   for (r = 0; r < grid.dim; r++)
   {
      setvalarray(uint[r],tindex,1.0+alpha[r]);
      tindex[r] = mid-sk[r];
      setvalarray(uint[r],tindex,-alpha[r]);
      tindex[r] = mid;
   }
   for (r = 0; r < grid.dim; r++)
      if (abs(sk[r]) == 1)
      {
         s = sk[r];
         newstorage(Dusmall[buildsize]);
         Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
         Dusmall[buildsize].info[1] = r;
         Dusmall[buildsize].info[2] = s;
         Dusmall[buildsize].info[3] = -1;
         Dusmall[buildsize].mid = mid;
         value = 0.5*grid.dx[r]*grid.dx[r]*d[r]*alpha[r]*(1.0+alpha[r]);
         if (fabs(value) > tol)
            sparseorder(-1,Dusmall[buildsize].head,value);
         for (i = 0; i < grid.dim; i++)
            tindex[i] = 0;
         while (tindex[0] <= N)
         {
            value = evalarray(uint[r],tindex)+
                    0.5*grid.dx[r]*grid.dx[r]*evalarray(f[r],tindex)*alpha[r]*
                    (1.0+alpha[r]);
            if (fabs(value) > tol)
            {
               if (globsmall == 1)
                  sparseorder(sub2ind(tindex,Narray,grid.dim),Dusmall[buildsize].head,
                              value);
               else if (globsmall == 2)
               {
                  for (i = 0; i < grid.dim; i++)
                     rindex[i] = index[i]+tindex[i]-mid;
                  sparseorder(sub2ind(rindex,grid.nx,grid.dim),Dusmall[buildsize].head,
                              value);
               }
            }

            (tindex[grid.dim-1])++;
            for (i = grid.dim-1; i > 0 && tindex[i] > N; i--)
            {
               tindex[i] = 0;
               (tindex[i-1])++;
            }
         }
         buildsize++;

         for (t = 0; t < grid.dim; t++)
         {
            newstorage(Dusmall[buildsize]);
            Dusmall[buildsize].info[0] = sub2ind(index,grid.nx,grid.dim);
            Dusmall[buildsize].info[1] = r;
            Dusmall[buildsize].info[2] = s;
            Dusmall[buildsize].info[3] = t;
            Dusmall[buildsize].mid = mid;
            value = uxxcoef[r][t]*d[t]/grid.dx[r];
            if (fabs(value) > tol)
               sparseorder(-1,Dusmall[buildsize].head,value);

            for (i = 0; i < grid.dim; i++)
               tindex[i] = 0;
            while (tindex[0] <= N)
            {
               value = (evalarray(Du[r][t],tindex)+
                        uxxcoef[r][t]*evalarray(f[t],tindex))/grid.dx[r];
               if (fabs(value) > tol)
               {
                  if (globsmall == 1)
                     sparseorder(sub2ind(tindex,Narray,grid.dim),Dusmall[buildsize].head,
                                 value);
                  else if (globsmall == 2)
                  {
                     for (i = 0; i < grid.dim; i++)
                        rindex[i] = index[i]+tindex[i]-mid;
                     sparseorder(sub2ind(rindex,grid.nx,grid.dim),
                                 Dusmall[buildsize].head,value);
                  }
               }

               (tindex[grid.dim-1])++;
               for (i = grid.dim-1; i > 0 && tindex[i] > N; i--)
               {
                  tindex[i] = 0;
                  (tindex[i-1])++;
               }
            }
            buildsize++;
         }
      }

   free_matrix(M,grid.dim-1,grid.dim-1);
   free_matrix(LU,grid.dim-1,grid.dim-1);
   free_matrix(f,grid.dim-1,4,4,4);
   for (r = 0; r < grid.dim; r++)
      free_matrix(Du[r],grid.dim-1,N,N,N);
   delete [] Du;
   free_matrix(uint,grid.dim-1,4,4,4);
}


/*
check max du in each dimension of cim12
*/
void checkcim12Du(double ***u, double ***S, PBData &pb, GridData &grid)
{
   int i, r, s, t, m, thestatus;
   int tindex[grid.dim], rindex[grid.dim], **zindex = imatrix(grid.dim-1,grid.dim-1);
   double uint, Du[grid.dim];
   double alpha, thesign, tangent[grid.dim], normal[grid.dim];
   double tmperr, theerr[grid.dim];
   for (t = 0; t < grid.dim; t++)
      theerr[t] = 0.0;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      if (evalarray(S,tindex) < 0.0)
         thesign = -1.0;
      else
         thesign = 1.0;
      for (r = 0; r < grid.dim; r++)
         rindex[r] = tindex[r];
      for (r = 0; r < grid.dim; r++)
      {
         for (s = -1; s <= 1; s += 2)
         {
            rindex[r] = tindex[r]+s;
            if (rindex[r] >= 0 && rindex[r] <= grid.nx[r] &&
                (evalarray(S,tindex) < 0.0)+(evalarray(S,rindex) < 0.0) == 1)
            {
               thestatus = getstatus3(S,tindex,grid);
               getcim12Du(uint,Du,tindex,r,s,u,S,pb,grid);

               getinterfaceinfo(alpha,tangent,normal,S,tindex,r,s,grid);
//               if (thestatus == 2)
               if (1)
                  for (t = 0; t < grid.dim; t++)
                  {
                     tmperr = fabs(Du[t]-getDu(tindex,t,r,s,alpha,thesign,grid));
                     if (tmperr > theerr[t])
                     {
                        theerr[t] = tmperr;
                        for (m = 0; m < grid.dim; m++)
                           zindex[t][m] = tindex[m];
                     }
                  }
            }
         }
         rindex[r] = tindex[r];
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

   cout << "[checkcim12Du]Error in cim12 derivative is " << theerr[0] << " " << theerr[1] << " "
        << theerr[2] << endl;
   cout << "   at";
   for (t = 0; t < grid.dim; t++)
   {
      for (m = 0; m < grid.dim; m++)
         cout << " " << zindex[t][m];
      cout << " (" << getstatus3(S,zindex[t],grid) << ")";
      if (t < grid.dim-1)
         cout << ", ";
      else
         cout << "." << endl;
   }

   free_matrix(zindex,grid.dim-1,grid.dim-1);
}


void getcim1Du(double &uint, double *Du, int *index, int rstar, int sstar, double ***u,
               int gamma[][2], double ***S, PBData &pb, GridData &grid)
{
   int r, s, j, k;
   int rindex[grid.dim];
   double **LU, **M, value;
   int PLR[grid.dim], PLC[grid.dim];
   double alpha, beta, tangent[grid.dim], normal[grid.dim], temp[grid.dim];
   double ethere, ehere, ebar;
// ADDING
   double sigma, thesign, tau, Dtau[grid.dim], finalalpha;

   LU = matrix(grid.dim-1,grid.dim-1);
   M = matrix(grid.dim-1,grid.dim-1);

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

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];

   s = sstar;
   for (r = 0; r < grid.dim; r++)
   {
      Du[r] = 0.0;
      rindex[r] = index[r]+s;
      getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
      rindex[r] = index[r];
      if (r == rstar)
         finalalpha = alpha;
      beta = 1.0-alpha;
      ebar = alpha*ethere+beta*ehere;
      for (j = 0; j < grid.dim; j++)
         M[r][j] = -gamma[r][(s+1)/2]*gamma[j][(s+1)/2]*beta*(ehere-ethere)/ebar*
                    tangent[r]*tangent[j];
      M[r][r] += 1.0;

      for (j = 0; j < grid.dim; j++)
      {
         rindex[j] = index[j]+s;
         Du[r] += (gamma[r][(s+1)/2]*beta*(ehere-ethere)*tangent[r]*
                   (1.0-gamma[j][(s+1)/2])*s*tangent[j])/(grid.dx[r]*ebar)*
                  evalarray(u,rindex);
         rindex[j] = index[j];
      }
      rindex[r] = index[r]+s;
      Du[r] += (s*ethere)/(grid.dx[r]*ebar)*evalarray(u,rindex);
      rindex[r] = index[r];
      Du[r] += -s*ethere/(grid.dx[r]*ebar)*evalarray(u,index);
      for (j = 0; j < grid.dim; j++)
         Du[r] += -gamma[r][(s+1)/2]*beta*(ehere-ethere)*tangent[r]*
                  (1.0-gamma[j][(s+1)/2])*s*tangent[j]/(grid.dx[r]*ebar)*
                  evalarray(u,index);
// ADDING
      getsigma(sigma,index,r,s,alpha,normal,pb,grid);
      gettau(tau,index,r,s,alpha,grid);
      getDtau(Dtau,index,r,s,alpha,grid);
      Du[r] += thesign*gamma[r][(s+1)/2]*
               (beta/ebar*normal[r]*sigma+s*ethere/(grid.dx[r]*ebar)*tau+
                ethere*beta/ebar*tangent[r]*getdotprod(Dtau,tangent,grid.dim));
   }
   gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);
   forwardbacksub0(Du,Du,LU,PLR,PLC,grid.dim-1);

   uint = evalarray(u,index)+Du[rstar]*sstar*finalalpha*grid.dx[rstar];

   free_matrix(M,grid.dim-1,grid.dim-1);
   free_matrix(LU,grid.dim-1,grid.dim-1);
}

void getcim2Du(double &uint, double *Du, int *index, int rstar, int sstar, double ***u,
               int gamma[][2], double ***S, PBData &pb, GridData &grid)
{
   int r, s, i, j, k, tmps;
   int rindex[grid.dim], tindex[grid.dim];
   double ****f;
   double **LU, **M, value, Sval;
   int PLR[grid.dim], PLC[grid.dim];
   double alpha, beta, tangent[grid.dim], normal[grid.dim];
   double bk, rthere, rhere, ethere, ehere, ehat, C;
   double temp[grid.dim], a[4];
   int sk[grid.dim];
// ADDING
   double sigma, d[grid.dim], thesign, tau, Dtau[grid.dim], finalalpha;

   LU = matrix(grid.dim-1,grid.dim-1);
   M = matrix(grid.dim-1,grid.dim-1);
   f = matrix(grid.dim-1,4,4,4);

   Sval = evalarray(S,index);
   if (Sval < 0.0)
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

   for (r = 0; r < grid.dim; r++)
      sk[r] = gamma[r][1]-gamma[r][0];
// assuming sk[rstar] == sstar

   for (r = 0; r < grid.dim; r++)
   {
      for (i = 0; i < grid.dim; i++)
         tindex[i] = 0;
      while (tindex[0] < 5)
      {
         setvalarray(f[r],tindex,0.0);

         (tindex[grid.dim-1])++;
         for (i = grid.dim-1; i > 0 && tindex[i] >= 5; i--)
         {
            tindex[i] = 0;
            (tindex[i-1])++;
         }
      }

// ADDING
      d[r] = 0.0;
   }

   double uxx[grid.dim], uxxcoef[grid.dim];
   for (r = 0; r < grid.dim; r++)
   {
      Du[r] = 0.0;
      uxx[r] = 0.0;
      uxxcoef[r] = 0.0;
   }
   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
      tindex[r] = 2;
   for (r = 0; r < grid.dim; r++)
   {
      rindex[r] = index[r]+sk[r];
      getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
      if (r == rstar)
         finalalpha = alpha;
      rindex[r] = index[r];
      if (abs(sk[r]) == 1)
      {
         beta = 1.0-alpha;
         ehat = (beta+beta*beta)*(0.5+alpha)*ehere+(alpha+alpha*alpha)*(0.5+beta)*ethere;
         rhere = ehere/ehat;
         rthere = ethere/ehat;
         a[0] = (beta+beta*beta)*rhere+alpha*(1.0+2.0*beta)*rthere;
         a[1] = -(beta+beta*beta)*rhere-(1.0+alpha)*(1.0+2.0*beta)*rthere;
         a[2] = (1.0+beta)*(1.0+beta)*rthere;
         a[3] = -beta*beta*rthere;
         bk = -(beta+beta*beta)*(rthere-rhere);
         C = bk*tangent[r];
         for (s = 0; s < 4; s++)
         {
            rindex[r] = index[r]+(s-1)*sk[r];
            uxx[r] += a[s]*evalarray(u,rindex);
         }
         rindex[r] = index[r]-sk[r];
         uxx[r] += -sk[r]*C*tangent[r]*sk[r]*evalarray(u,rindex);
         if (r == rstar)
            Du[r] += -sk[r]*evalarray(u,rindex);
         rindex[r] = index[r];
         uxx[r] += sk[r]*C*tangent[r]*sk[r]*evalarray(u,rindex);
         if (r == rstar)
            Du[r] += sk[r]*evalarray(u,rindex);

         M[r][r] = 1.0-abs(sk[r])*(0.5+alpha)*bk*tangent[r]*tangent[r];
         if (r == rstar)
            uxxcoef[r] = sk[r]*(0.5+alpha);

         for (j = 0; j < grid.dim; j++)
            if (j != r)
            {
               tmps = sk[j];
               if (tmps == 0)
                  for (s = -1; s <= 1 && tmps == 0; s += 2)
                  {
                     rindex[j] = index[j]+s;
                     if (evalarray(S,rindex)*Sval <= 0.0)
                        tmps = s;
                  }
               if (tmps == 0)
               {
                  rindex[r] = index[r]-sk[r];
                  for (s = -1; s <= 1 && tmps == 0; s += 2)
                  {
                     rindex[j] = index[j]+s;
                     if (evalarray(S,rindex)*Sval <= 0.0)
                        tmps = s;
                  }
               }
               rindex[r] = index[r];
               rindex[j] = index[j];

               M[r][j] = -0.5*tmps*sk[r]*bk*tangent[j]*tangent[r];
               if (r == rstar)
                  uxxcoef[j] = 0.5*tmps;

               if (abs(tmps) == 1)
               {
                  rindex[j] = index[j]-tmps;
                  uxx[r] += -sk[r]*C*tangent[j]*(1.0+alpha)*tmps*evalarray(u,rindex);
                  if (r == rstar)
                     Du[j] += -tmps*(1.0+alpha)*evalarray(u,rindex);
                  rindex[r] = index[r]-sk[r];
                  uxx[r] += sk[r]*C*tangent[j]*alpha*tmps*evalarray(u,rindex);
                  if (r == rstar)
                     Du[j] += tmps*alpha*evalarray(u,rindex);
                  rindex[j] = index[j];
                  uxx[r] += -sk[r]*C*tangent[j]*alpha*tmps*evalarray(u,rindex);
                  if (r == rstar)
                     Du[j] += -tmps*alpha*evalarray(u,rindex);
                  rindex[r] = index[r];
                  uxx[r] += sk[r]*C*tangent[j]*(1.0+alpha)*tmps*evalarray(u,rindex);
                  if (r == rstar)
                     Du[j] += tmps*(1.0+alpha)*evalarray(u,rindex);
               }
               else
               {
                  for (s = -1; s <= 1; s += 2)
                  {
                     rindex[j] = index[j]+s;
                     uxx[r] += sk[r]*C*tangent[j]*(1.0+alpha)*0.5*s*evalarray(u,rindex);
                     if (r == rstar)
                        Du[j] += s*(1.0+alpha)*0.5*evalarray(u,rindex);
                  }
                  rindex[r] = index[r]-sk[r];
                  for (s = -1; s <= 1; s += 2)
                  {
                     rindex[j] = index[j]+s;
                     uxx[r] += -sk[r]*C*tangent[j]*alpha*0.5*s*evalarray(u,rindex);
                     if (r == rstar)
                        Du[j] += -s*alpha*0.5*evalarray(u,rindex);
                  }
                  rindex[r] = index[r];
                  rindex[j] = index[j];
               }
            }

// ADDING
         getsigma(sigma,index,r,sk[r],alpha,normal,pb,grid);
         gettau(tau,index,r,sk[r],alpha,grid);
         getDtau(Dtau,index,r,sk[r],alpha,grid);

         uxx[r] += thesign*(abs(sk[r])*(1.0+2.0*beta)*rthere*tau+
                            sk[r]*(beta+beta*beta)*grid.dx[r]*
                            (sigma/ehat*normal[r]+
                             rthere*getdotprod(Dtau,tangent,grid.dim)*tangent[r]));
      }
      else
      {
         for (j = 0; j < grid.dim; j++)
            if (j == r)
               M[r][r] = 1.0;
            else
               M[r][j] = 0.0;
         for (s = -1; s <= 1; s += 2)
         {
            rindex[r] = index[r]+s;
            uxx[r] += evalarray(u,rindex);
         }
         rindex[r] = index[r];
         uxx[r] += -2.0*evalarray(u,rindex);
      }
   }

   gecp0(LU,PLR,PLC,M,grid.dim-1,grid.dim-1);
   forwardbacksub0(uxx,uxx,LU,PLR,PLC,grid.dim-1);
   for (r = 0; r < grid.dim; r++)
      Du[r] += uxxcoef[r]*uxx[r];
   for (r = 0; r < grid.dim; r++)
      Du[r] /= grid.dx[rstar];

   rindex[rstar] = index[rstar]-sstar;
   uint = evalarray(u,index)+
          ((evalarray(u,index)-evalarray(u,rindex))/grid.dx[rstar]+
           0.5*grid.dx[rstar]*uxx[rstar])*finalalpha*grid.dx[rstar]+
          0.5*uxx[rstar]*finalalpha*finalalpha*grid.dx[rstar]*grid.dx[rstar];

   free_matrix(M,grid.dim-1,grid.dim-1);
   free_matrix(LU,grid.dim-1,grid.dim-1);
   free_matrix(f,grid.dim-1,4,4,4);
}

void getcim12Du(double &uint, double *Du, int *index, int rstar, int sstar,
                double ***u, double ***S, PBData &pb, GridData &grid)
{
   int r, s, rindex[grid.dim], gamma[grid.dim][2];

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
            gamma[r][(s+1)/2] = 0;
      }
      rindex[r] = index[r];
   }

   if (getstatus3(S,index,grid) == 3)
      getcim1Du(uint,Du,index,rstar,sstar,u,gamma,S,pb,grid);
   else if (getstatus3(S,index,grid) == 2)
      getcim2Du(uint,Du,index,rstar,sstar,u,gamma,S,pb,grid);
}


// status 0 = boundary, 1 = interior, 2 = cim2, 3 = exceptional
int getstatus3(double ***S, int *index, GridData grid)
{
   int r, s, j, k, sk;
   int tindex[grid.dim];
   double Sval1, Sval2, Sval3, Sval4, Sval5;
   char interior = 1, except = 0, failed;

   for (r = 0; r < grid.dim && index[r] > 0 && index[r] < grid.nx[r]; r++); // r=3 if interior
   if (r < grid.dim)
      return 0;

   for (r = 0; r < grid.dim; r++)
      tindex[r] = index[r];
   for (r = 0; r < grid.dim; r++)
   {
      sk = 0;
      for (s = -1; s <= 1; s += 2)
      {
         tindex[r] = index[r]+s;
         Sval1 = evalarray(S,index);//Sval1 is current point
         Sval2 = evalarray(S,tindex);//Sval2 is nbr +/- 1
         if (Sval1*Sval2 < 0.0)//change sign in r dim
         {
            sk += 1;
            if (sk == 2) // if sk = 2, then both front and back are in different side
               return 3;
            tindex[r] = index[r]+2*s;
            Sval3 = evalarray(S,tindex); //Sval3 = dimr + 2s
            if (Sval2*Sval3 < 0.0) // less then 2 points in the other side
               return 3;
            interior = 0;
         }
         else// not yet return, 2 points in both side, look for cross derivative
            for (j = 0; j < grid.dim; j++) 
               if (j != r) // j is a dimension different thant r
               {
                  failed = 0;
                  for (k = -1; k <= 1; k += 2)
                  {
                     tindex[r] = index[r];
                     tindex[j] = index[j]+k; // tindex is diag
                     Sval3 = evalarray(S,tindex);
                     if (Sval1*Sval3 >= 0.0)// if diagonal has same sign
                     {
                        tindex[r] = index[r]+s;
                        Sval4 = evalarray(S,tindex);
                        if (Sval3*Sval4 < 0.0)
                           failed += 1; //fail to apprx cross derivative
                     }
                     else 
                        failed += 1;//change sign in dim=j
                     tindex[j] = index[j];
                  }
                  if (failed == 2)
                     except = 1; // can not approximate cross derivative
               }
      }
      tindex[r] = index[r];
   }

   if (interior)
      return 1;// interior point
   else if (except)
      return 3;//exceptional
   else
      return 2;// can use CIM2
}



double getinterfacegrad4(double *grad, double ***u, double ***S, int *index, int rstar,
                         int sstar, PBData &pb, GridData grid)
{                  
   int r, s, rindex[grid.dim], gamma[grid.dim][2];
   double uint;
   
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
            gamma[r][(s+1)/2] = 0;
      }
      rindex[r] = index[r];
   }

   if (getstatus3(S,index,grid) == 3)
      getcim1Du(uint,grad,index,rstar,sstar,u,gamma,S,pb,grid);
   else if (getstatus3(S,index,grid) == 2)
      getcim2Du(uint,grad,index,rstar,sstar,u,gamma,S,pb,grid);

   return uint;
}