#include "cim345cond.h"
#include "global.h"
#include "helper.h"
#include "interface.h"
#include "ccim.h"
#include "numerics.h"
#include "tryvn.h"
#include <iostream>
#include <cmath>
using namespace std;





// output ux at grid point and Du at interface
// modifed from getcim345Du, add *ux as argument
void getcim345DuAll(double &uint, double *ux, double *Du, int *index, int rstar, int sstar,
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
   double uxx[grid.dim];

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

// check Du at gridpoint
void checkcim345DuAll(double ***u, double ***S, PBData &pb, GridData &grid)
{   
   double uint, ux[grid.dim], Du[grid.dim]; //results from getcim345
   
   array<int,3> maxErrIdx;
   array<int,3> maxErrItfIdx;
   double maxErr; // maxErr at all grid point
   double maxErrItf; // maxErr at grid point near interface
   vector<double> sqrdCompErrAll; // vector of error component at all grid point
   vector<double> sqrdCompErrItf; // vector of error component at grid point near interface
   
   double max_Dup_itf_err = 0.0;// error at interface plus side
   array<int,3> max_Dup_itf_err_idx;

   double max_Dup_itf_rel_err = 0.0;// relative error at interface plus side
   array<int,3> max_Dup_itf_rel_err_idx;
   
   double max_Dum_itf_err = 0.0;// error at interface minus side
   array<int,3> max_Dum_itf_err_idx;

   double max_Dum_itf_rel_err = 0.0;// relative error at interface minus side
   array<int,3> max_Dum_itf_rel_err_idx;
   
   // do not go through boundary
   for(int i = 1; i < grid.nx[0]; i++){
      for(int j = 1; j < grid.nx[1]; j++){
         for(int k = 1; k < grid.nx[2]; k++){
            int tindex[3] = {i,j,k};
            int rindex[3] = {i,j,k};

            int thestatus = getstatus5(S,tindex,grid);
            double thesign = (evalarray(S,tindex) < 0.0)? -1.0 : 1.0;
            
            double alpha, tangent[grid.dim], normal[grid.dim];
            if (thestatus>=2){ 
               //cim point, find rstar and sstar
               for (int r = 0; r < grid.dim; r++){
                  for (int s = -1; s <= 1; s += 2){
                     rindex[r] = tindex[r]+s;

                     if ((evalarray(S,tindex) < 0.0)+(evalarray(S,rindex) < 0.0) == 1){
                        
                        getcim345DuAll(uint,ux,Du,tindex,r,s,u,S,pb,grid);
                        
                        getinterfaceinfo(alpha,tangent,normal,S,tindex,r,s,grid);
                        
                        //Du exact at interface
                        array<double,3> du_itf_exact {getDu(tindex,0,r,s,alpha,thesign,grid),
                                                getDu(tindex,1,r,s,alpha,thesign,grid),
                                                getDu(tindex,2,r,s,alpha,thesign,grid)};
                        double max_Du_itf_exact = max(max(fabs(du_itf_exact[0]),fabs(du_itf_exact[1])),fabs(du_itf_exact[2]));
                        //Du at interface

                        array<double,3> Du_itf_err =  { fabs(Du[0] - du_itf_exact[0]),
                                                   fabs(Du[1] - du_itf_exact[1]),
                                                   fabs(Du[2] - du_itf_exact[2])};
                        double max_Du_itf_err = max(max(Du_itf_err[0],Du_itf_err[1]),Du_itf_err[2]);
                        
                        double Du_itf_rel_err = max_Du_itf_err/max_Du_itf_exact; //relative error
                        
                        if(thesign<0){
                           max_Dum_itf_err = (max_Du_itf_err > max_Dum_itf_err)? max_Du_itf_err : max_Dum_itf_err;

                           max_Dum_itf_err_idx = array<int,3>{i,j,k};

                           max_Dum_itf_rel_err = (Du_itf_rel_err > max_Dum_itf_rel_err)? Du_itf_rel_err : max_Dum_itf_rel_err;
                           max_Dum_itf_rel_err_idx = array<int,3>{i,j,k};
                        }else{
                           max_Dup_itf_err = (max_Du_itf_err > max_Dup_itf_err)? max_Du_itf_err : max_Dup_itf_err;
                           max_Dup_itf_err_idx = array<int,3>{i,j,k};

                           max_Dup_itf_rel_err = (Du_itf_rel_err > max_Dup_itf_rel_err)? Du_itf_rel_err : max_Dup_itf_rel_err;
                           max_Dup_itf_rel_err_idx = array<int,3>{i,j,k};
                        }
                     }else{
                        continue;
                     }
                  }//end of s loop
                  rindex[r] = tindex[r];
               }// end of r loop


            }else if (thestatus == 1){ 
               //interior, central differencing
               ux[0] = (u[i+1][j][k] - u[i-1][j][k]) / (2*grid.dx[0]);
               ux[1] = (u[i][j+1][k] - u[i][j-1][k]) / (2*grid.dx[1]);
               ux[2] = (u[i][j][k+1] - u[i][j][k-1]) / (2*grid.dx[2]);

            }else{
               cerr<<"check Du unclassified"<<endl;
               exit(1);
            }

            

            double xerr = abs(ux[0] - getDu(tindex,0, 0, 0, 0.0, thesign, grid));
            double yerr = abs(ux[1] - getDu(tindex,1, 0, 0, 0.0, thesign, grid));
            double zerr = abs(ux[2] - getDu(tindex,2, 0, 0, 0.0, thesign, grid));
            double tmpMaxErr = max(max(xerr,yerr),zerr);

            if (tmpMaxErr>maxErr){
               maxErr = tmpMaxErr;
               maxErrIdx = array<int,3>{i,j,k};
            }
            sqrdCompErrAll.push_back(xerr);
            sqrdCompErrAll.push_back(yerr);
            sqrdCompErrAll.push_back(zerr);

            if(thestatus>=2){
               if (tmpMaxErr>maxErrItf){
                  maxErrItf = tmpMaxErr;
                  maxErrItfIdx = array<int,3>{i,j,k};
               }
               sqrdCompErrItf.push_back(xerr);
               sqrdCompErrItf.push_back(yerr);
               sqrdCompErrItf.push_back(zerr);
               
            }
         }//end k
      }// end j
   }// end i
    
   cout<<"[checkcim345DuAll]"<<endl;
   printf("max Du err All = %f at (%d,%d,%d)\n",maxErr,maxErrIdx[0],maxErrIdx[1],maxErrIdx[2]);
   printf("max Du err Interface grid point = %f at (%d,%d,%d)\n",maxErrItf,maxErrItfIdx[0],maxErrItfIdx[1],maxErrItfIdx[2]);
   printf("rmse Du err All = %f \n",rmse(sqrdCompErrAll));
   printf("rmse Du err Interface grid point = %f \n",rmse(sqrdCompErrItf));

   printf("max Du err negative interface = %f at (%d,%d,%d)\n", max_Dum_itf_err, max_Dum_itf_err_idx[0], max_Dum_itf_err_idx[1], max_Dum_itf_err_idx[2]);
   printf("max Du err positive interface = %f at (%d,%d,%d)\n", max_Dup_itf_err, max_Dup_itf_err_idx[0], max_Dup_itf_err_idx[1], max_Dup_itf_err_idx[2]);

   printf("max Du rel err negative interface = %f at (%d,%d,%d)\n", max_Dum_itf_rel_err, max_Dum_itf_rel_err_idx[0], max_Dum_itf_rel_err_idx[1], max_Dum_itf_rel_err_idx[2]);
   printf("max Du rel err positive interface = %f at (%d,%d,%d)\n", max_Dup_itf_rel_err, max_Dup_itf_rel_err_idx[0], max_Dup_itf_rel_err_idx[1], max_Dup_itf_rel_err_idx[2]);
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
   
   if (equal(index,eindex))
   {
      cout << "S" << endl;
      for (t = -1; t <= 1; t++) 
      {
         for (s = 1; s >= -1; s--) 
         {
            for (r = -1; r <= 1; r++) 
               printf("%4.16f ",S[index[0]+r][index[1]+s][index[2]+t]);
            cout << endl;
         }
         cout << endl;
      }
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
                  if (equal(index,eindex))
                  {
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
               cout << "computed Du" << endl;
               for (m = 0; m < grid.dim; m++)
                  cout<<"dim "<< m <<" apprx = " 
                      << evalcoef(D1u[r][(sk+1)/2][m],D1ucoef[r][(sk+1)/2][m],
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
               cout << "computed jump in D1u = "
                    << evalcoef(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef, 
                                jumpD1jumpuxxcoef,index,0,0,0.0,S,grid) << " "
                    <<", exact = "
                    << getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                       getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid) 
                    // << "another exact "
                    // << -(pb.epsilonp-pb.epsilonm)/ethere*
                    //     (getDu(index,0,r,sk,alpha[r][(sk+1)/2],thesign,grid)*normal[0]+
                    //      getDu(index,1,r,sk,alpha[r][(sk+1)/2],thesign,grid)*normal[1]+
                    //      getDu(index,2,r,sk,alpha[r][(sk+1)/2],thesign,grid)*normal[2])*
                    //     normal[r] << " "
                    // << 2.0*x[r]*(1.0-pb.epsilonp/pb.epsilonm) 
                    << " error =  "
                    << evalcoef(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,
                                jumpD1jumpuxxcoef,index,0,0,0.0,S,grid)-
                       (getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                        getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid)) << endl;

               rhs += -thesign*sk*beta*grid.dx[r]*
                      evalcoef(jumpD1u,jumpD1ucoef,jumpD1uxcoef,jumpD1uxxcoef,
                               jumpD1jumpuxxcoef,index,0,0,0.0,S,grid);
               cout << "rhs 2 = " << rhs << " error =  "
                    << sk*grid.dx[r]*getDu(index,r,r,sk,alpha[r][(sk+1)/2],thesign,grid)-
                       thesign*sk*beta*grid.dx[r]*
                        (getDu(index,r,r,sk,alpha[r][(sk+1)/2],1,grid)-
                         getDu(index,r,r,sk,alpha[r][(sk+1)/2],-1,grid)) << endl; //rhs2 = h ux- + beta h [ux]
            }
//            cout << "getting jump D2u" << endl;
            
            if(a){
              // if a is not nullptr
              getcim345jumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,
                             alpha[r][(sk+1)/2],thesign,normal,mid,a,D1u[r][(sk+1)/2],
                             D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                             D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,D2u,D2ucoef,D2uxcoef,
                             D2uxxcoef,D2jumpuxxcoef,jumpD1u,jumpD1ucoef,jumpD1uxcoef,
                             jumpD1uxxcoef,jumpD1jumpuxxcoef,S,pb,grid);
            }else{
              getcim345jumpuxx(jumpD2u,jumpD2ucoef,jumpD2uxcoef,jumpD2uxxcoef,index,r,sk,
                             alpha[r][(sk+1)/2],thesign,normal,mid,D1u[r][(sk+1)/2],
                             D1ucoef[r][(sk+1)/2],D1uxcoef[r][(sk+1)/2],
                             D1uxxcoef[r][(sk+1)/2],D1jumpuxxcoef,D2u,D2ucoef,D2uxcoef,
                             D2uxxcoef,D2jumpuxxcoef,jumpD1u,jumpD1ucoef,jumpD1uxcoef,
                             jumpD1uxxcoef,jumpD1jumpuxxcoef,S,pb,grid);  
            }
            
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
               cout << "(at grid) u_here - u_there "
                    << " = " << getu(index,0,0,0.0,thesign,grid) 
                    << " - " << getu(index,r,sk,1.0,-thesign,grid)
                    << " = "<<getu(index,0,0,0.0,thesign,grid) - getu(index,r,sk,1.0,-thesign,grid) << endl;
               cout << "u jump = "
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
            
            if (equal(index,eindex))
            {
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
            exactd[grid.dim*(sk+1)/2+r] += -1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(index,0,0,0.0,thesign,grid);
            sindex[r] = mid+sk;
            setvalarray(dcoef[grid.dim*(sk+1)/2+r],sindex,1.0/(grid.dx[r]*grid.dx[r]));
            exactd[grid.dim*(sk+1)/2+r] += 1.0/(grid.dx[r]*grid.dx[r])*
                                           getu(index,r,sk,1.0,thesign,grid);
            sindex[r] = mid;
         }

   // print G matrix, exact rhs
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

      cout<<"solved | exact"<<endl;
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

   globdist = 0;
   gmatrix(G, LU, PLR, PLC,
    d0, dcoef, 
    D1u, D1ucoef,  D1uxcoef,  D1uxxcoef,
    D2u,  D2ucoef,
    index,
    alpha,  gamma, S, a, pb, grid);
   double estcond0 = condest(G, 6);


   globdist = 1;   
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


   // reset
   globdist = tempglobdist;
  
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
         free_matrix(D2ucoef,N,N,N);
         delete [] D2uxcoef;
         delete [] D2uxxcoef;
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
         free_matrix(D2ucoef,N,N,N);
         delete [] D2uxcoef;
         delete [] D2uxxcoef;
      }
   if (offsetvec.size() == 0){
      return 0;
   }
   return 1;
   
}


// get status of point index: 0 = boundary, 1 = interior, 2 = cim3, 3 = cim5, 4 = cim4, 5 = cim1
int getstatus5debug(double ***S, int *index, GridData &grid)
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

// old incorrect version returned cim 5 if cim 5 for any m and n 
   junk1 = matrix(2,2,2);
//   thecim = 3;
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
      {
         cimstatus[m][n] = 3; 
         if (!yessk2(sk2,m,n,index,S,grid))//if do not have usual mixed derivative
         {

            if (yescim5D2(junk1,junk2,junk3,m,n,index,mid,S,grid))//if has cim5 mixed derivative
//               if (thecim == 3)
//                  thecim = 5;
               cimstatus[m][n] = 5; 
            else 
            {
//               if (thecim == 3 || thecim == 5)
//                  thecim = 4;
               cimstatus[m][n] = 4; 
               for (r = 0; r < grid.dim; r++)
                  for (sk = -1; sk <= 1; sk += 2)
                  {
                     rindex[r] = index[r]+sk;
                     if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1 &&
                         !yessk2(sk2,m,n,rindex,S,grid))
//                        thecim = 0;
                        cout<<"m = "<<m<<", n = "<<n<<", r = "<<r<<", sk = "<<sk<<", cimstatus = "<<cimstatus[m][n]<<endl;
                        cimstatus[m][n] = 0;
                     rindex[r] = index[r];
                  }
            }
         }
      }
   free_matrix(junk1,2,2,2);

   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++){
         cout<<"m = "<<m<<", n = "<<n<<", cimstatus = "<<cimstatus[m][n]<<endl;
      }



   thecim = 3;
   for (m = 0; m < grid.dim; m++)
      for (n = m+1; n < grid.dim; n++)
         if (cimstatus[m][n] == 0)
            return 5;
         else if (cimstatus[m][n] == 4)
            thecim = 4;
         else if (thecim == 3 && cimstatus[m][n] == 5)
            thecim = 5;
   
   cout<<"thecim = "<<thecim<<endl;

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
            if ((evalarray(S,index) < 0.0)+(evalarray(S,rindex) < 0.0) == 1)//if change sign, count++
               count++;
         }
         rindex[r] = index[r];
      }
      if (count < 2)
         return 2;
      else
         thecim = 5;
   }

   if (thecim == 5)
      return 3;
   else if (thecim == 4)
      return 4;
   else
      return 5; 
}

