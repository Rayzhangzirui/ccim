#include "getvn.h"
#include "cim12.h"
#include "ccim.h"
#include "interface.h"
#include "numerics.h"
#include <iostream>

using namespace std;

double getpsivn2(double ***u, double ***S, int *index, int rstar, PBData &pb, 
                 GridData grid)
{
   int r;
   int rindex[grid.dim], tindex[grid.dim];
   double grad[grid.dim], normal[grid.dim], tangent[grid.dim];
   double alpha, ehere, vn, temp;
   double uint;
   int status1, status2;

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   rindex[rstar] = index[rstar]+1;
   status1 = getstatus3(S,index,grid);
   status2 = getstatus3(S,rindex,grid);
   if (status1 < status2 || 
       (status1 == status2 && 
        ((evalarray(S,index) >= 0.0 && pb.epsilonp >= pb.epsilonm) || 
         (evalarray(S,index) <= 0.0 && pb.epsilonm >= pb.epsilonp))))
      for (r = 0; r < grid.dim; r++)
         tindex[r] = index[r];
   else
      for (r = 0; r < grid.dim; r++)
      {
         tindex[r] = rindex[r];
         rindex[r] = index[r];
      }
   getinterfaceinfo(alpha,tangent,normal,S,tindex,rindex,grid);
   if (evalarray(S,tindex) < 0.0)
      ehere = pb.epsilonm;
   else
      ehere = pb.epsilonp;

   uint = getinterfacegrad4(grad,u,S,tindex,rstar,rindex[rstar]-tindex[rstar],pb,grid);
   temp = 0.0;
   for (r = 0; r < grid.dim; r++)
      temp += grad[r]*normal[r];
   temp *= ehere;
   vn = -0.5*(1.0/pb.epsilonm-1.0/pb.epsilonp)*temp*temp;
   project(tangent,normal,grad,grid.dim);
   temp = 0.0;
   for (r = 0; r < grid.dim; r++)
      temp += tangent[r]*tangent[r];
   vn += -0.5*(pb.epsilonp-pb.epsilonm)*temp;
   vn += -Bval(uint,pb);

   return vn;
}

double getLJvn(double ***u, double ***S, int *index, int rstar, PBData &pb, 
               GridData grid)
{
   int r;
   int rindex[grid.dim];
   double grad[grid.dim], normal[grid.dim], tangent[grid.dim];
   double alpha, ehere, ethere, thesign, v1, v2, vn, temp;

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   rindex[rstar] = index[rstar]+1;
   getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);

   return pb.rho0*weno6interp(u,index,rstar,1,alpha,grid);
}

double getvn(double ***psi, double ***LJ, double ***S, int *index, int rstar, 
             PBData &pb, GridData grid)
{
   return getpsivn2(psi,S,index,rstar,pb,grid)+getLJvn(LJ,S,index,rstar,pb,grid);
}


// get interface vn at index, rstar is which dimension, sstar is +1 or -1 change sign in dim r
double getheatvn(double ***u, double ***S, int *index, int rstar, int sstar, 
                 PBData &pb, GridData &grid)
{
   int r, rindex[grid.dim];
   double thesign, uhere, uthere, vn, Duhere[grid.dim], Duthere[grid.dim];
   double alpha, tangent[grid.dim], normal[grid.dim];

   if (evalarray(S,index) < 0.0)
      thesign = -1.0;
   else
      thesign = 1.0;

   getcim345Du(uhere,Duhere,index,rstar,sstar,u,S,pb,grid);
   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   rindex[rstar] = index[rstar]+sstar;
   getcim345Du(uthere,Duthere,rindex,rstar,-sstar,u,S,pb,grid);

   getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);

   if (thesign < 0.0)
      vn = formvn(Duthere,Duhere,normal,grid);// if index in inside, 'there' is Du in outside,
   else
      vn = formvn(Duhere,Duthere,normal,grid);// if index in outside, 'here' is Du in outside,

   return vn;
}

// getvn using Dusmall StorageStruct, 
double getheatvn(double ***u, double ***S, int *index, int rstar, int sstar, 
                 StorageStruct *Dusmall, int smallsize, PBData &pb, GridData &grid)
{
   int r, rindex[grid.dim];
   double thesign, uhere, uthere, vn, Duhere[grid.dim], Duthere[grid.dim];
   double alpha, tangent[grid.dim], normal[grid.dim];

   if (evalarray(S,index) < 0.0)
      thesign = -1.0;
   else
      thesign = 1.0;


   evalfromstorage(uhere,Duhere,index,rstar,sstar,Dusmall,smallsize,u,S,pb,grid);

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   rindex[rstar] = index[rstar]+sstar;

   evalfromstorage(uthere,Duthere,rindex,rstar,-sstar,Dusmall,smallsize,u,S,pb,grid);

   getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);

   if (thesign < 0.0)
      vn = formvn(Duthere,Duhere,normal,grid);
   else
      vn = formvn(Duhere,Duthere,normal,grid);

   double exactvn, x[grid.dim];
   sub2coord(x,index,grid);
   x[rstar] += sstar*alpha*grid.dx[rstar];
   if (globtestnum == 0)
      exactvn = 4.0*grid.radius/((1.0+grid.radius*grid.radius)*
                                 (1.0+grid.radius*grid.radius));
   else if (globtestnum == 1)
      exactvn = 2.0*(1.0-pb.epsilonp/pb.epsilonm)*0.5*
                exp(2.0*(1.0-pb.epsilonp/pb.epsilonm)*grid.t);
   else if (globtestnum == 2)
      exactvn = 2.0*grid.radius*exp(grid.radius*grid.radius)*
                (1.0-pb.epsilonp/pb.epsilonm);
   else if (globtestnum == 3 || globtestnum == 4)
      exactvn = (1.0-pb.epsilonp/pb.epsilonm)*grid.radius/
                sqrt(grid.radius*grid.radius+fabs(1.0-pb.epsilonp/pb.epsilonm));
   if (fabs(exactvn-vn) > grid.maxvn)
   {
      grid.maxvn = fabs(exactvn-vn);
   }

   return vn;
}

double getrishuvn(double ***u, double ***S, int *index, int rstar, int sstar, 
                  StorageStruct *Dusmall, int smallsize, PBData &pb, GridData &grid)
{
   int r, rindex[grid.dim];
   double thesign, uhere, uthere, vn, Duhere[grid.dim], Duthere[grid.dim];
   double alpha, tangent[grid.dim], normal[grid.dim]; 
   double x[grid.dim], Dpsivac[grid.dim];

   if (evalarray(S,index) < 0.0)
      thesign = -1.0;
   else
      thesign = 1.0;
   getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
   sub2coord(x,index,grid);
   x[rstar] += sstar*alpha*grid.dx[rstar];

   evalfromstorage(uhere,Duhere,index,rstar,sstar,Dusmall,smallsize,u,S,pb,grid);
   uhere += getpsivac(x,pb);
   getDpsivac(Dpsivac,x,pb);
   for (r = 0; r < grid.dim; r++)
      Duhere[r] += Dpsivac[r];

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   rindex[rstar] = index[rstar]+sstar;
   evalfromstorage(uthere,Duthere,rindex,rstar,-sstar,Dusmall,smallsize,u,S,pb,grid);
   uthere += getpsivac(x,pb);
   getDpsivac(Dpsivac,x,pb);
   for (r = 0; r < grid.dim; r++)
      Duthere[r] += Dpsivac[r];

   if (thesign < 0.0)
      vn = formvn(Duthere,Duhere,normal,grid);
   else
      vn = formvn(Duhere,Duthere,normal,grid);

   return vn;
}

double getexactvn(int *index, int r, int s, double ***S, GridData &grid)
{
   int t, rindex[grid.dim];
   double vn, alpha, tangent[grid.dim], normal[grid.dim], Dup[grid.dim], Dum[grid.dim];

   for (t = 0; t < grid.dim; t++)
      rindex[t] = index[t];
   rindex[r] = index[r]+s;
   getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);

   for (t = 0; t < grid.dim; t++)
   {
      Dup[t] = getDu(index,t,r,s,alpha,1,grid);
      Dum[t] = getDu(index,t,r,s,alpha,-1,grid);
   }

   vn = formvn(Dup,Dum,normal,grid);

   return vn;
}
// vn = dot(Dup,normal)- dot(Dum,normal),
double formvn(double *Dup, double *Dum, double *normal, GridData &grid)
{
   int r;
   double vn = 0.0;

   for (r = 0; r < grid.dim; r++)
      vn += Dup[r]*normal[r];
   for (r = 0; r < grid.dim; r++)
      vn -= Dum[r]*normal[r];

   return vn;
//   return -vn;
}

double LJwater(int *index, PBData &pb, GridData &grid)
{
   double temp, val = 0.0;
   double epsilon, sigma, tol, x[pb.dim];
   int r, s;

   sub2coord(x,index,grid);

   for (r = 0; r < pb.N; r++)
   {
      epsilon = pb.epsilonlj[r];
      sigma = pb.sigmalj[r];
      temp = 0.0;
      for (s = 0; s < grid.dim; s++)
         temp += (x[s]-pb.x[r][s])*(x[s]-pb.x[r][s]);
      temp = sqrt(temp);
      tol = 0.0;
      for (s = 0; s < grid.dim; s++)
         tol += grid.dx[s]*grid.dx[s];
      tol = sqrt(tol);
      if (temp < tol)
         temp = tol;
      val += 4.0*epsilon*(exp(12.0*log(sigma/temp))-exp(6.0*log(sigma/temp)));
   }

   return val;
}