#include "interface.h"
#include "numerics.h"
#include "input.h"
#include "global.h"
#include "helper.h"
#include <cmath>
#include <iostream>

using namespace std;


double getdotprod(double *v, double *w, int thedim)
{
   double value;
   int r;

   value = 0.0;
   for (r = 0; r < thedim; r++)
      value += v[r]*w[r];

   return value;
}

void normalize(double *v, double *w, int thedim)
{
// can use same array for v and w.
   double length;
   int r;

   length = sqrt(getdotprod(w,w,thedim));
   for (r = 0; r < thedim; r++)
      v[r] = w[r]/length;
}

void getcrossprod(double *v, double *w, double *z)
{
// can use same array for v and w or z.
   double temp[3];

   temp[0] = w[1]*z[2]-w[2]*z[1];
   temp[1] = -w[0]*z[2]+w[2]*z[0];
   temp[2] = w[0]*z[1]-w[1]*z[0];

   v[0] = temp[0];
   v[1] = temp[1];
   v[2] = temp[2];
}

void getunitcrossprod(double *v, double *w, double *z)
{
// can use same array for v and w or z.
   getcrossprod(v,w,z);
   normalize(v,v,3);
}

// w = v - proj(v to normal)
// project v to the plane with normal
void project(double *w, double *normal, double *v, int dim)
{
   int r;
   double dotprod, length2;

   dotprod = 0.0;
   for (r = 0; r < dim; r++)
      dotprod += v[r]*normal[r]; //v dot normal
   length2 = 0.0;
   for (r = 0; r < dim; r++)
      length2 += normal[r]*normal[r]; //length of normal
   for (r = 0; r < dim; r++)
      w[r] = v[r]-dotprod*normal[r]/length2;
}

char getinterfaceinfo(double &alpha, double *tangent, double *normal, double ***S,
                      int *index, int rstar, int sstar, GridData &grid)
{
   int r;
   int rindex[grid.dim];

   for (r = 0; r < grid.dim; r++)
      rindex[r] = index[r];
   rindex[rstar] = index[rstar]+sstar;
   return getinterfaceinfo(alpha,tangent,normal,S,index,rindex,grid);
}

// index1 = index of current point, index2 = index across interface, 
// tangent is projection of (index2-index1) to normal
char getinterfaceinfo(double &alpha, double *tangent, double *normal, double ***S,
                      int *index1, int *index2, GridData &grid)
{
   int r, s, rstar, sstar, rindex1[grid.dim], rindex2[grid.dim];
   double length, tol = 1.0e-14, normal1[grid.dim], normal2[grid.dim];

   for (r = 0; r < grid.dim; r++)// rindex1 copy index1
   {
      rindex1[r] = index1[r];
      rindex2[r] = index2[r];
   }
   for (r = 0; r < grid.dim; r++)
   {
      normal1[r] = 0.0;
      normal2[r] = 0.0;
   }
   if (evalarray(S,index1)*evalarray(S,index2) < 0.0)//if one of index1/2 is inside
   {
      for (rstar = 0; rstar < grid.dim && index1[rstar] == index2[rstar]; rstar++);//rstar = which dimension to diffrence
      sstar = index2[rstar]-index1[rstar];

      if (globregfalweno)
         alpha = regulafalsiweno(S,index1,rstar,sstar,grid);
      else
         alpha = regulafalsi(S,index1,rstar,sstar,grid);
// USING EXACT
//      alpha = bisection(S,index1,rstar,sstar,grid);
//      alpha = regulafalsiexact(index1,rstar,sstar,grid);
      for (r = 0; r < grid.dim; r++)
      {
         normal1[r] = 0.0;
         normal2[r] = 0.0;
         for (s = -1; s <= 1; s += 2)
         {
            
            rindex1[r] = index1[r]+s;
            rindex2[r] = index2[r]+s;
            normal1[r] += s*evalarray(S,rindex1);
            normal2[r] += s*evalarray(S,rindex2);
         }
         normal1[r] /= 2.0*grid.dx[r];//central differencing in dimension r
         normal2[r] /= 2.0*grid.dx[r];
         normal[r] = (1.0-alpha)*normal1[r]+alpha*normal2[r];//linear interpolate the first derivative
         rindex1[r] = index1[r];
         rindex2[r] = index2[r];
      }

      //normalization of normal vector
      length = 0.0;
      for (r = 0; r < grid.dim; r++)
         length += normal[r]*normal[r];
      length = sqrt(length);
      for (r = 0; r < grid.dim; r++)
         normal[r] /= length;
// USING EXACT
//      for (r = 0; r < grid.dim; r++)
//         normal[r] = getexactnormal(index1,r,rstar,sstar,alpha,grid);

      for (r = 0; r < grid.dim; r++)
         tangent[r] = index2[r]-index1[r];
      project(tangent,normal,tangent,grid.dim); //tangent is projection of u(nit vector crossing interface) to the plane with normal
      
      length = 0.0;
      for (r = 0; r < grid.dim; r++)
         length += tangent[r]*tangent[r];
      length = sqrt(length);
      
      if (length > tol)
         for (r = 0; r < grid.dim; r++)
            tangent[r] /= length;
      else
         for (r = 0; r < grid.dim; r++){
            tangent[r] = 0.0; 
            if(abs(normal[r])<tol) {// if tangent = 0, that means normal is parallel to tindex2-index1, make no differece as long as tk ek = 0 in cim2
              tangent[r] = 1.0;
              break;
            }
         }
      return 1;
   }
   else
   {
      alpha = 1.0;
      for (r = 0; r < grid.dim; r++)
      {
         normal[r] = 0.0;
         tangent[r] = 0.0;
      }
      return 0;
   }
}

char getinterfaceinfo(double &alpha, double *tangent1, double *tangent2, double *normal, 
                      double &sigma, double **Dn, double *Dsigma, double &jumpfe, 
                      double ***S, int *index1, int *index2, PBData &pb, GridData &grid)
{
   int r, s, t, m, n, rstar, sstar, rindex1[grid.dim], rindex2[grid.dim];
   double length, tol = 1.0e-14, normal1[grid.dim], normal2[grid.dim], x[grid.dim];
   double Dn2[grid.dim][grid.dim];
   double length2;

   double realDn[grid.dim][grid.dim], realD2[grid.dim][grid.dim];
   double realnormal[grid.dim];

   for (r = 0; r < grid.dim; r++)
   {
      rindex1[r] = index1[r];
      rindex2[r] = index2[r];
   }
   for (r = 0; r < grid.dim; r++)
   {
      normal1[r] = 0.0;
      normal2[r] = 0.0;
   }
   if (evalarray(S,index1)*evalarray(S,index2) < 0.0)
   {
      for (rstar = 0; rstar < grid.dim && index1[rstar] == index2[rstar]; rstar++);
      sstar = index2[rstar]-index1[rstar];

      if (globregfalweno)
         alpha = regulafalsiweno(S,index1,rstar,sstar,grid);
      else
         alpha = regulafalsi(S,index1,rstar,sstar,grid);
// USING EXACT
//      alpha = bisection(S,index1,rstar,sstar,grid);
//      alpha = regulafalsiexact(index1,rstar,sstar,grid);
      for (r = 0; r < grid.dim; r++)
      {
         normal1[r] = 0.0;
         normal2[r] = 0.0;
         for (s = -1; s <= 1; s += 2)
         {
            rindex1[r] = index1[r]+s;
            rindex2[r] = index2[r]+s;
            normal1[r] += s*evalarray(S,rindex1);
            normal2[r] += s*evalarray(S,rindex2);
         }
         normal1[r] /= 2.0*grid.dx[r];
         normal2[r] /= 2.0*grid.dx[r];
         normal[r] = (1.0-alpha)*normal1[r]+alpha*normal2[r];
         rindex1[r] = index1[r];
         rindex2[r] = index2[r];
      }
      length = 0.0;
      for (r = 0; r < grid.dim; r++)
         length += normal[r]*normal[r];
      length = sqrt(length);
      for (r = 0; r < grid.dim; r++)
      {
         normal[r] /= length;
         normal2[r] = normal[r];
      }
      length2 = length;
      for (r = 0; r < grid.dim; r++)
         tangent1[r] = index2[r]-index1[r];
      project(tangent1,normal,tangent1,grid.dim);
      length = 0.0;
      for (r = 0; r < grid.dim; r++)
         length += tangent1[r]*tangent1[r];
      length = sqrt(length);
      if (length > tol)
         for (r = 0; r < grid.dim; r++)
            tangent1[r] /= length;
      else
         for (r = 0; r < grid.dim; r++)
            tangent1[r] = 0.0;
      getunitcrossprod(tangent2,normal,tangent1);

      for (r = 0; r < grid.dim; r++)
      {
         Dn[r][r] = 0.0;
         Dn2[r][r] = 0.0;
         for (s = -1; s <= 1; s += 2)
         {
            Dn[r][r] -= evalarray(S,index1);
            rindex1[r] = index1[r]+s;
            Dn[r][r] += evalarray(S,rindex1);
            Dn2[r][r] -= evalarray(S,index2);
            rindex2[r] = index2[r]+s;
            Dn2[r][r] += evalarray(S,rindex2);
         }
         rindex1[r] = index1[r];
         rindex2[r] = index2[r];
         Dn[r][r] /= grid.dx[r]*grid.dx[r];
         Dn2[r][r] /= grid.dx[r]*grid.dx[r];
/*
         sub2coord(x,index1,grid);
         cout << Dn[r][r] << " " << (1.0-x[r]*x[r]/(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))/
                                    sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) << endl;
         sub2coord(x,index2,grid);
         cout << Dn2[r][r] << " " << (1.0-x[r]*x[r]/(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))/
                                     sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]) << endl;
*/
         Dn2[r][r] = (1.0-alpha)*Dn[r][r]+alpha*Dn2[r][r];
         for (m = r+1; m < grid.dim; m++)
         {
            Dn[r][m] = 0.0;
            Dn2[r][m] = 0.0;
            for (s = -1; s <= 1; s += 2)
            {
               rindex1[r] = index1[r]+s;
               rindex2[r] = index2[r]+s;
               for (t = -1; t <= 1; t += 2)
               {
                  rindex1[m] = index1[m]+t;
                  rindex2[m] = index2[m]+t;
                  Dn[r][m] += s*t*evalarray(S,rindex1);
                  Dn2[r][m] += s*t*evalarray(S,rindex2);
               }
               rindex1[m] = index1[m];
               rindex2[m] = index2[m];
            }
            rindex1[r] = index1[r];
            rindex2[r] = index2[r];
            Dn[r][m] /= 4.0*grid.dx[r]*grid.dx[m];
            Dn2[r][m] /= 4.0*grid.dx[r]*grid.dx[m];
            Dn2[r][m] = (1.0-alpha)*Dn[r][m]+alpha*Dn2[r][m];
            Dn2[m][r] = Dn2[r][m];
         }
      }
      sub2coord(x,index1,grid);
      x[rstar] += sstar*alpha*grid.dx[rstar];
/*
      for (r = 0; r < grid.dim; r++)
         for (m = 0; m < grid.dim; m++)
         {
            cout << "   in Dn " << r << " " << m << " " << Dn2[r][m] << " ";
            if (r != m)
               realD2[r][m] = -x[r]*x[m]/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2])/
                              (x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
            else
               realD2[r][m] = (1.0-x[r]*x[m]/(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))/
                              sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
            cout << realD2[r][m] << endl;
         }
*/
/*
      cout << "   in Dn " << normal2[0] << " " << normal2[1] << " " << normal2[2] << endl;
      for (r = 0; r < grid.dim; r++)
         realnormal[r] = x[r]/sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
      cout << "   in Dn " << realnormal[0] << " " << realnormal[1] << " " 
           << realnormal[2] << endl;
*/
      for (r = 0; r < grid.dim; r++)
         for (m = 0; m < grid.dim; m++)
         {
            Dn[r][m] = Dn2[r][m]/length2;
            for (n = 0; n < grid.dim; n++)
               Dn[r][m] -= normal2[r]*Dn2[m][n]*normal2[n]/length2;
//            realDn[r][m] = realD2[r][m];
//            for (n = 0; n < grid.dim; n++)
//               realDn[r][m] -= realnormal[r]*realD2[m][n]*realnormal[n];
//            cout << "in Dn " << realDn[r][m] << endl;
         }
      jumpfe = getf(index1,rstar,sstar,alpha,1,pb,grid)/pb.epsilonp-
               getf(index1,rstar,sstar,alpha,-1,pb,grid)/pb.epsilonm;
      getsigma(sigma,index1,rstar,sstar,alpha,normal,pb,grid);
      getDsigma(Dsigma,index1,rstar,sstar,alpha,normal,Dn,pb,grid);
      return 1;
   }
   else
   {
      alpha = 1.0;
      for (r = 0; r < grid.dim; r++)
      {
         normal[r] = 0.0;
         tangent1[r] = 0.0;
      }
      return 0;
   }
}

void getinterfaceinfo(double *tangent1, double *tangent2, double &sigma, double **Dn,
                      double *Dsigma, double &jumpfe, int *index, int rstar, int sstar,
                      double alpha, double ***S, PBData &pb, GridData &grid)
{
   int r, s, t, m, n, index2[grid.dim], rindex1[grid.dim], rindex2[grid.dim];
   double grad, length, tol = 1.0e-14, normal[grid.dim], normal2[grid.dim];
   double Dn2[grid.dim][grid.dim];

   for (r = 0; r < grid.dim; r++)
   {
      rindex1[r] = index[r]; // rindex1 = index = current
      index2[r] = index[r];
   }
   index2[rstar] += sstar;
   for (r = 0; r < grid.dim; r++)
      rindex2[r] = index2[r];//rindex2 = index2 = index on the other side
   for (r = 0; r < grid.dim; r++)
   {
      normal[r] = 0.0;
      normal2[r] = 0.0;
      for (s = -1; s <= 1; s += 2)
      {
         rindex1[r] = index[r]+s;
         rindex2[r] = index2[r]+s;
         normal[r] += s*evalarray(S,rindex1);
         normal2[r] += s*evalarray(S,rindex2);
      }
      normal[r] /= 2.0*grid.dx[r];
      normal2[r] /= 2.0*grid.dx[r]; //calculate late normal at both side by central diff
      normal[r] = (1.0-alpha)*normal[r]+alpha*normal2[r];//linear interpolate to get interface
      rindex1[r] = index[r];
      rindex2[r] = index2[r];
   }
   grad = sqrt(getdotprod(normal,normal,grid.dim));
   for (r = 0; r < grid.dim; r++)
      normal[r] /= grad;
// USING EXACT
//   for (r = 0; r < grid.dim; r++)
//      normal[r] = getexactnormal(index,r,rstar,sstar,alpha,grid);
   length = 0.0;
   t = -1;
   for (r = 0; r < grid.dim; r++) 
   {   
      for (s = 0; s < grid.dim; s++)
         tangent1[s] = 0.0;
      tangent1[r] = 1.0; //go through basis vector and find orthonal projection as tangent
      project(tangent1,normal,tangent1,grid.dim);
      if (sqrt(getdotprod(tangent1,tangent1,grid.dim)) > length)//find the one with the largest length, ie, angle closest to 90
      {
         length = sqrt(getdotprod(tangent1,tangent1,grid.dim));
         t = r;
      }
   }
// can change this
//   for (r = 0; r < grid.dim; r++)
//      tangent1[r] = 0.0;
//   tangent1[rstar] = sstar;
   for (s = 0; s < grid.dim; s++)
      tangent1[s] = 0.0;
   tangent1[t] = 1.0;
   project(tangent1,normal,tangent1,grid.dim);
   length = sqrt(getdotprod(tangent1,tangent1,grid.dim));
   for (s = 0; s < grid.dim; s++)
      tangent1[s] /= length;
   getunitcrossprod(tangent2,normal,tangent1);

   for (r = 0; r < grid.dim; r++)
   {
      Dn[r][r] = 0.0;
      Dn2[r][r] = 0.0;
      for (s = -1; s <= 1; s += 2)
      {
         Dn[r][r] -= evalarray(S,index);
         rindex1[r] = index[r]+s;
         Dn[r][r] += evalarray(S,rindex1);
         Dn2[r][r] -= evalarray(S,index2);
         rindex2[r] = index2[r]+s;
         Dn2[r][r] += evalarray(S,rindex2);
      }
      rindex1[r] = index[r];
      rindex2[r] = index2[r];
      Dn[r][r] /= grid.dx[r]*grid.dx[r]; //Dn[r][r] is DDphi/DDr at index by central diff
      Dn2[r][r] /= grid.dx[r]*grid.dx[r];//Dn2[r][r] is DDphi/DDr at index2 the other side
      Dn2[r][r] = (1.0-alpha)*Dn[r][r]+alpha*Dn2[r][r];//linear interpolate at interface
      for (m = r+1; m < grid.dim; m++)
      {
         Dn[r][m] = 0.0;
         Dn2[r][m] = 0.0;
         for (s = -1; s <= 1; s += 2)
         {
            rindex1[r] = index[r]+s;
            rindex2[r] = index2[r]+s;
            for (t = -1; t <= 1; t += 2)
            {
               rindex1[m] = index[m]+t;
               rindex2[m] = index2[m]+t;
               Dn[r][m] += s*t*evalarray(S,rindex1);
               Dn2[r][m] += s*t*evalarray(S,rindex2);
            }
            rindex1[m] = index[m];
            rindex2[m] = index2[m];
         }
         rindex1[r] = index[r];
         rindex2[r] = index2[r];
         Dn[r][m] /= 4.0*grid.dx[r]*grid.dx[m]; //Dn[r][m] is DDphi/DrDm at index by central diff
         Dn2[r][m] /= 4.0*grid.dx[r]*grid.dx[m];
         Dn2[r][m] = (1.0-alpha)*Dn[r][m]+alpha*Dn2[r][m];
         Dn2[m][r] = Dn2[r][m];
      }
   }
   for (r = 0; r < grid.dim; r++)
      for (m = 0; m < grid.dim; m++)
      {
         Dn[r][m] = Dn2[r][m]/grad;
         for (n = 0; n < grid.dim; n++)
            Dn[r][m] -= normal[r]*Dn2[m][n]*normal[n]/grad; //grad n
// USING EXACT
//         Dn[r][m] = getexactDnormal(index,r,m,rstar,sstar,alpha,grid);
      }
   jumpfe = getf(index,rstar,sstar,alpha,1,pb,grid)/pb.epsilonp-
            getf(index,rstar,sstar,alpha,-1,pb,grid)/pb.epsilonm;
   getsigma(sigma,index,rstar,sstar,alpha,normal,pb,grid);
   getDsigma(Dsigma,index,rstar,sstar,alpha,normal,Dn,pb,grid);
}

// with jumpfe = [f/eps], aehere, 
void getinterfaceinfo(double *tangent1, double *tangent2, double &tau, double &sigma, 
                      double **Dn, double *Dsigma, double &jumpfe, double &aehere,
                      double &aethere, int *index, int rstar, int sstar, double alpha, 
                      double ***a, double ***S, PBData &pb, GridData &grid)
{
   int r, s, t, m, n, index2[grid.dim], rindex1[grid.dim], rindex2[grid.dim];
   double grad, length, tol = 1.0e-14, normal[grid.dim], normal2[grid.dim];
   double Dn2[grid.dim][grid.dim];

   for (r = 0; r < grid.dim; r++)
   {
      rindex1[r] = index[r];
      index2[r] = index[r];
   }
   index2[rstar] += sstar;
   for (r = 0; r < grid.dim; r++)
      rindex2[r] = index2[r];
   for (r = 0; r < grid.dim; r++)
   {
      normal[r] = 0.0;
      normal2[r] = 0.0;
      for (s = -1; s <= 1; s += 2)
      {
         rindex1[r] = index[r]+s;
         rindex2[r] = index2[r]+s;
         normal[r] += s*evalarray(S,rindex1);
         normal2[r] += s*evalarray(S,rindex2);
      }
      normal[r] /= 2.0*grid.dx[r];
      normal2[r] /= 2.0*grid.dx[r];
      normal[r] = (1.0-alpha)*normal[r]+alpha*normal2[r];
      rindex1[r] = index[r];
      rindex2[r] = index2[r];
   }
   grad = sqrt(getdotprod(normal,normal,grid.dim));
   for (r = 0; r < grid.dim; r++)
      normal[r] /= grad;
// USING EXACT
//   for (r = 0; r < grid.dim; r++)
//      normal[r] = getexactnormal(index,r,rstar,sstar,alpha,grid);
   length = 0.0;
   t = -1;
   for (r = 0; r < grid.dim; r++)
   {   
      for (s = 0; s < grid.dim; s++)
         tangent1[s] = 0.0;
      tangent1[r] = 1.0;
      project(tangent1,normal,tangent1,grid.dim);
      if (sqrt(getdotprod(tangent1,tangent1,grid.dim)) > length)
      {
         length = sqrt(getdotprod(tangent1,tangent1,grid.dim));
         t = r;
      }
   }
// can change this
//   for (r = 0; r < grid.dim; r++)
//      tangent1[r] = 0.0;
//   tangent1[rstar] = sstar;
   for (s = 0; s < grid.dim; s++)
      tangent1[s] = 0.0;
   tangent1[t] = 1.0;
   project(tangent1,normal,tangent1,grid.dim);
   length = sqrt(getdotprod(tangent1,tangent1,grid.dim));
   for (s = 0; s < grid.dim; s++)
      tangent1[s] /= length;
   getunitcrossprod(tangent2,normal,tangent1);

   for (r = 0; r < grid.dim; r++)
   {
      Dn[r][r] = 0.0;
      Dn2[r][r] = 0.0;
      for (s = -1; s <= 1; s += 2)
      {
         Dn[r][r] -= evalarray(S,index);
         rindex1[r] = index[r]+s;
         Dn[r][r] += evalarray(S,rindex1);
         Dn2[r][r] -= evalarray(S,index2);
         rindex2[r] = index2[r]+s;
         Dn2[r][r] += evalarray(S,rindex2);
      }
      rindex1[r] = index[r];
      rindex2[r] = index2[r];
      Dn[r][r] /= grid.dx[r]*grid.dx[r];
      Dn2[r][r] /= grid.dx[r]*grid.dx[r];
      Dn2[r][r] = (1.0-alpha)*Dn[r][r]+alpha*Dn2[r][r];
      for (m = r+1; m < grid.dim; m++)
      {
         Dn[r][m] = 0.0;
         Dn2[r][m] = 0.0;
         for (s = -1; s <= 1; s += 2)
         {
            rindex1[r] = index[r]+s;
            rindex2[r] = index2[r]+s;
            for (t = -1; t <= 1; t += 2)
            {
               rindex1[m] = index[m]+t;
               rindex2[m] = index2[m]+t;
               Dn[r][m] += s*t*evalarray(S,rindex1);
               Dn2[r][m] += s*t*evalarray(S,rindex2);
            }
            rindex1[m] = index[m];
            rindex2[m] = index2[m];
         }
         rindex1[r] = index[r];
         rindex2[r] = index2[r];
         Dn[r][m] /= 4.0*grid.dx[r]*grid.dx[m];
         Dn2[r][m] /= 4.0*grid.dx[r]*grid.dx[m];
         Dn2[r][m] = (1.0-alpha)*Dn[r][m]+alpha*Dn2[r][m];
         Dn2[m][r] = Dn2[r][m];
      }
   }
   for (r = 0; r < grid.dim; r++)
      for (m = 0; m < grid.dim; m++)
      {
         Dn[r][m] = Dn2[r][m]/grad;
         for (n = 0; n < grid.dim; n++)
            Dn[r][m] -= normal[r]*Dn2[m][n]*normal[n]/grad;
// USING EXACT
//         Dn[r][m] = getexactDnormal(index,r,m,rstar,sstar,alpha,grid);
      }
   jumpfe = getf(index,rstar,sstar,alpha,1,pb,grid)/pb.epsilonp-
            getf(index,rstar,sstar,alpha,-1,pb,grid)/pb.epsilonm;
   gettau(tau,index,rstar,sstar,alpha,grid);
   getsigma(sigma,index,rstar,sstar,alpha,normal,pb,grid);
   getDsigma(Dsigma,index,rstar,sstar,alpha,normal,Dn,pb,grid);
   if (evalarray(S,index) < 0.0)
   {
      aehere = evalarray(a,index)/pb.epsilonm;
      aethere = evalarray(a,index2)/pb.epsilonp;
   }
   else
   {
      aehere = evalarray(a,index)/pb.epsilonp;
      aethere = evalarray(a,index2)/pb.epsilonm;
   }
}



// go through all grid points, calculated radius, display maximum error
void checkwithexact(double ***S, double radius, GridData &grid)
{
   int i, r, s, t, tindex[grid.dim], rindex[grid.dim];
   double value, theerr = 0.0, alpha, x[grid.dim], tangent[grid.dim], normal[grid.dim];
   double y[grid.dim], morerr[grid.dim];
   vector<double> rvec; 
   for (t = 0; t < grid.dim; t++)
      morerr[t] = 0.0;

   for (i = 0; i < grid.dim; i++)
      tindex[i] = 0;
   while (tindex[0] <= grid.nx[0])
   {
      s = 1;
      for (r = 0; r < grid.dim; r++)
      {
         for (i = 0; i < grid.dim; i++)
            rindex[i] = tindex[i];
         rindex[r] = tindex[r]+s;
         if (rindex[r] <= grid.nx[r] && 
             (evalarray(S,tindex) < 0.0)+(evalarray(S,rindex) < 0.0) == 1)
         {
            sub2coord(x,tindex,grid);
//            cout << tindex[0] << " " << tindex[1] << " " << tindex[2] << endl;
//            cout << rindex[0] << " " << rindex[1] << " " << rindex[2] << endl;
            getinterfaceinfo(alpha,tangent,normal,S,tindex,rindex,grid);
            x[r] += s*alpha*grid.dx[r];//x is interface location
            value = 0.0;
            for (i = 0; i < grid.dim; i++)
               value += x[i]*x[i];
            value = sqrt(value);// value = norm(x) = radius
            if (tindex[0] == grid.nx[0]/2 && tindex[1] == grid.nx[1]/2 && 
                tindex[2] >= grid.nx[2]/2 && r == 2)// error of radius along positive z axis if grid.nx is even
               cout << "Especially z-axis " << value << " " << radius << " " 
                    << fabs(value-radius) << " " << x[2] << endl;
//            mu2d << value-radius << " ";
            rvec.push_back(fabs(value));
            if (fabs(value-radius) > theerr) 
            {
               theerr = fabs(value-radius);
               for (t = 0; t < grid.dim; t++)
                  y[t] = x[t];
            }
            for (t = 0; t < grid.dim; t++)
               if (fabs(normal[t]-x[t]/value) > morerr[t])//suppose surface is a sphere, then normal should be same as x/norm(x)
                  morerr[t] = fabs(normal[t]-x[t]/value);
         }
      }

      (tindex[grid.dim-1])++;
      for (i = grid.dim-1; i > 0 && tindex[i] > grid.nx[i]; i--)
      {
         tindex[i] = 0;
         (tindex[i-1])++;
      }
   }

//   mu2d << endl;
   cout << "radius mean = " << mean(rvec) << " variance = " << variance(rvec) << " rmse =" << rmse(rvec, radius) <<endl;
   cout << "Max error for radius = " << theerr << " at " << y[0] << " " << y[1] << " "
        << y[2] << endl;
   cout << "Max error for normal = " << morerr[0] << " " << morerr[1] << " " 
        << morerr[2] << endl; 
   // epsilon2d << grid.t << " " << theerr << " " << morerr[0] << " " << morerr[1] << " "
   //           << morerr[2] << endl;
}


